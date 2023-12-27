// Original Author: Rad1ance
// This version: Jeff Luo

#include "analyze_trans.h"
#include <fstream>

AnalyzeTrans::AnalyzeTrans(Equation *equation, CompMapType *compMap, std::set<int> *nodeList, float deltaT, float tFinal)
{
    // initialize the time related parameters
    initialTimePace(deltaT);
    initialFinalTime(tFinal);
    initialCurrentT();
    initialDeltaTUnchanged(deltaT);

    initialDeltaTOrigin(deltaT);

    // initialize the initial value of the transient analysis and transInductor and transCapacitor
    std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it1;

    int flagC = 0;
    int flagI = 0;
    for (it1 = equation->transformedComponents.begin(); it1 != equation->transformedComponents.end(); it1++)
    {
        // initialize std::map<const NonlinearComponent *, std::vector<BasicComponent *>>
        // focused on the Inductors and Conductors
        // in order to implement the first iteration of newton method
        if (it1->first->getType() == enumCapacitor)
        {
            flagC++;
            Capacitor *tempC = (Capacitor *)it1->first;
            TransCapacitor *c_tran = new TransCapacitor(tempC->getValue(), 0, 0,
                                                        tempC->getNodeIn(), tempC->getNodeOut(), equation->applyNewNode());
            Capacitors.push_back(c_tran);
            updateCapacitor(it1->first, it1->second, equation, Capacitors[flagC - 1], Capacitors[flagC - 1]->getNewNode());
        }
        if (it1->first->getType() == enumInductor)
        {
            flagI++;
            Inductor *tempI = (Inductor *)it1->first;
            TransInductor *l_tran = new TransInductor(tempI->getValue(), 0, 0,
                                                      tempI->getNodeIn(), tempI->getNodeOut());
            Inductors.push_back(l_tran);
            updateInductor(it1->first, it1->second, equation, Inductors[flagI - 1]);
        }
    }
    // update the MNA equation and go into newton method
    equation->equationAddTransformed();
    newtonMethod(nodeList, equation, compMap);
}

// calculate the recommended time pace of certain capacitor
void AnalyzeTrans::timeStepControlC(Equation *equation, float tPre, float &tNow, int flag, std::vector<float> Voltages)
{

    // Implementation of Dynamic Step Control
    int nin = std::distance(equation->nodeList->begin(), equation->nodeList->find(Capacitors[flag - 1]->getNodeIn())) - 1;
    int nout = std::distance(equation->nodeList->begin(), equation->nodeList->find(Capacitors[flag - 1]->getNodeOut())) - 1;
    int nnew = std::distance(equation->nodeList->begin(), equation->nodeList->find(Capacitors[flag - 1]->getNewNode())) - 1;

    float rc = parameters.back() / (2 * Capacitors[flag - 1]->getValue());
    float i_now = (Voltages[nin] - Voltages[nnew]) / rc;
    float i_pre = Capacitors[flag - 1]->getCurrent();

    float tol = (Capacitors[flag - 1]->getValue() / 1000) +
                (Capacitors[flag - 1]->getValue() / 100000) * std::max(fabs(i_now), fabs(i_pre)) * Capacitors[flag - 1]->getValue();

    // truncate error
    epsilon = -(deltaT * deltaT * (i_now - i_pre)) / (12 * Capacitors[flag - 1]->getValue()) * Capacitors[flag - 1]->getValue();

    // control factor
    float q = fabs(epsilon / tol);

    // while (q > 1)
    // {
    //     tNow = 0.9 * tPre * std::max(sqrt(1 / q), 0.3);
    //     epsilon = -(tNow * tNow * (i_now - i_pre)) / (12 * Capacitors[flag - 1]->getValue()) * Capacitors[flag - 1]->getValue();
    //     q = fabs(epsilon / tol);
    //     tPre = tNow;
    // }

    // if (q <= 1)
    //     tNow = tPre;
    if (q > q_max) q_max = q;
    if (q <= 1)
        tNow = tPre;
    else
        tNow = tPre * std::max(sqrt(1/q), 0.3);
}

// calculate the recommended time pace of certain inductor
void AnalyzeTrans::timeStepControlL(Equation *equation, float tPre, float &tNow, int flag, std::vector<float> Voltages)
{

    // Implementation of Dynamic Step Control
    int nin = std::distance(equation->nodeList->begin(), equation->nodeList->find(Inductors[flag - 1]->getNodeIn())) - 1;
    int nout = std::distance(equation->nodeList->begin(), equation->nodeList->find(Inductors[flag - 1]->getNodeOut())) - 1;

    float v_now = Voltages[nin] - Voltages[nout];
    float v_pre = Inductors[flag - 1]->getVoltage();

    // error tolerance (to avoid inf, multiply the value of itself)
    float tol = ((Inductors[flag - 1]->getValue() / 1000) +
                 (Inductors[flag - 1]->getValue() / 100000) * std::max(fabs(v_now), fabs(v_pre))) *
                Inductors[flag - 1]->getValue();

    // update epsilon
    epsilonL = -(deltaT * deltaT * (v_now - v_pre)) / (12 * Inductors[flag - 1]->getValue()) * Inductors[flag - 1]->getValue();

    float q = fabs(epsilonL / tol);


    // while (q > 1)
    // {
    //     tNow = 0.9 * tPre * std::max(sqrt(1 / q), 0.3);
    //     epsilonL = -(tNow * tNow * (v_now - v_pre)) / (12 * Inductors[flag - 1]->getValue()) * Inductors[flag - 1]->getValue();
    //     q = fabs(epsilonL / tol);
    //     tPre = tNow;
    // }

    // if (q <= 1)
    //     tNow = tPre;
    if (q > q_max) q_max = q;
    if (q <= 1)
        tNow = tPre;
    else
        tNow = tPre * std::max(sqrt(1/q), 0.3);
}

void AnalyzeTrans::newtonMethod(std::set<int> *nodeList, Equation *equation, CompMapType *compMap)
{
    int iteration = 0;
    int flag = 0;
    parameters.push_back(deltaTUnchanged);

    int solveflag = 1;
    std::vector<float> Voltages;

    // for output 
    int currentcnt = 0; 
    probeNameList = equation->getProbeNode();
    
    for (int i = 0; i < probeNameList.size(); i++)
    {
        int nodepos = std::distance(equation->nodeList->begin(), 
                        equation->nodeList->find(probeNameList[i])) - 1;
        probeList.push_back(nodepos);
    }
    

    while (currentT < tFinal)
    {
        iteration++;

        if (iteration % 10000 == 0)
            std::cout << "iteration reached: " << iteration << std::endl;

        if (solveflag)
        {
            nonlinearSolver = new NonlinearSolver(equation, compMap, nodeList, finalSolution);
            Voltages = nonlinearSolver->getFinalSolution();
            delete nonlinearSolver;
            nonlinearSolver = NULL;
        }
        

        // the first loop to implement dynamic time pace control
        std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it1;

        int flagC = 0, flagI = 0; // index of C & L. initial value is 0
        bool flagPaceChanged = 0;
        for (it1 = equation->transformedComponents.begin(); it1 != equation->transformedComponents.end(); it1++)
        {
            if (it1->first->getType() == enumCapacitor)
            {
                flagC++;
                timeStepControlC(equation, deltaTUnchanged, deltaT, flagC, Voltages);
                if (deltaT != deltaTUnchanged)
                {
                    flagPaceChanged = 1;
                    deltaTUnchanged = deltaT;
                    solveflag = 0;
                    break;
                }
            }
            if (it1->first->getType() == enumInductor)
            {
                flagI++;
                timeStepControlL(equation, deltaTUnchanged, deltaT, flagI, Voltages);
                if (deltaT != deltaTUnchanged)
                {
                    flagPaceChanged = 1;
                    deltaTUnchanged = deltaT;
                    solveflag = 0;
                    break;
                }
            }
        }

        // if the time pace is changed, give up results in this iteration
        if (flagPaceChanged)
        {
            continue;
        }

        flagC = 0;
        flagI = 0;
        solveflag = 1;
        parameters.pop_back(); // parameters.push_back(deltaTUnchanged);

        // if not, update all the results and parameters
        for (it1 = equation->transformedComponents.begin(); it1 != equation->transformedComponents.end(); it1++)
        {
            if (it1->first->getType() == enumCapacitor)
            {
                flagC++;
                updateVC_C(Capacitors[flagC - 1], Voltages, equation->nodeList, equation);
                updateCapacitor(it1->first, it1->second, equation, Capacitors[flagC - 1], Capacitors[flagC - 1]->getNewNode());
            }
            if (it1->first->getType() == enumInductor)
            {
                flagI++;
                updateVC_I(Inductors[flagI - 1], Voltages, equation->nodeList, equation);
                updateInductor(it1->first, it1->second, equation, Inductors[flagI - 1]);
            }
        }
        finalSolution = Voltages;
        currentT += deltaT;


        //output to file
        if (currentT >= currentcnt * deltaTOrigin)
        {
            // std::cout << "flag = " << flag << std::endl;
            // std::cout << "deltaT = " << deltaT << std::endl;
            currentcnt++;
            std::cout << "currentT: " << currentT << std::endl;

            // probeNameList and probeList have the same size
            for (int i = 0; i < probeList.size() ; i++)
                std::cout << "node " << probeNameList[i] <<": " << Voltages[probeList[i]] << std::endl;
        }

        // increase the time length each time control successs, 
        if (q_max < 0.2){
            deltaT = std::min(2 * deltaT, deltaTOrigin);
            deltaTUnchanged = deltaT;
        }
        q_max = 0;

        initialTimePace(deltaTUnchanged);
        equation->equationAddTransformed();
        parameters.push_back(deltaTUnchanged);
    }

    // output
    std::cout << std::endl;
    std::cout << "Transient analysis setting:" << std::endl;
    std::cout << "deltaT: " << deltaTOrigin << " s" <<std::endl;
    std::cout << "tFinal: " << tFinal << " s" <<std::endl;
    std::cout << std::endl;
    std::cout << "After " << iteration << " iterations, Newton Method was ended. " << std::endl;
    std::cout << "The voltages of the nodes:" << std::endl;
    std::set<int>::iterator temp = ++nodeList->begin();
    for (float i : finalSolution)
    {
        if (temp != nodeList->end())
        {
            std::cout << "node." << (*temp) << "\t" << i << std::endl;
            temp++;
        }
    }
}

// calculate the model according to the nonlinear solver outputs
void AnalyzeTrans::getParametersC(TransCapacitor *transCapacitor)
{

    //...
    float vc;
    float rc;
    vc = transCapacitor->getVoltage() + deltaT *
                                            transCapacitor->getCurrent() / (2 * transCapacitor->getValue());
    rc = deltaT / (2 * transCapacitor->getValue());
    resistorC = rc;
    voltageC = vc;
    parameters.push_back(vc);
    parameters.push_back(rc);
}

void AnalyzeTrans::getParametersI(TransInductor *transInductor)
{

    //...
    float currentl;
    float gl;
    currentl = transInductor->getCurrent() + deltaT *
                                                 transInductor->getVoltage() / (2 * transInductor->getValue());
    gl = deltaT / (2 * transInductor->getValue());
    resistorI = 1 / gl;
    currentI = currentl;
    parameters.push_back(currentl);
    parameters.push_back(gl);
}

// get the node voltages
float AnalyzeTrans::getVoltage(int index, const std::vector<float> &Voltages)
{
    return index == -1 ? 0 : Voltages[index];
}

// update the voltage and current of the C and I components
void AnalyzeTrans::updateVC_C(TransCapacitor *transCapacitor, const std::vector<float> Voltages, std::set<int> *nodeList,
                              Equation *equation)
{
    transCapacitor->setVoltage(
        getVoltage(equation->getNodeIndex(transCapacitor->getNodeIn()), Voltages) -
        getVoltage(equation->getNodeIndex(transCapacitor->getNodeOut()), Voltages));
    transCapacitor->setCurrent(
        (getVoltage(equation->getNodeIndex(transCapacitor->getNodeIn()), Voltages) -
         getVoltage(equation->getNodeIndex(transCapacitor->getNewNode()), Voltages)) /
        resistorC);
}

void AnalyzeTrans::updateVC_I(TransInductor *transInductor, const std::vector<float> Voltages, std::set<int> *nodeList,
                              Equation *equation)
{
    transInductor->setCurrent(
        (getVoltage(equation->getNodeIndex(transInductor->getNodeIn()), Voltages) -
         getVoltage(equation->getNodeIndex(transInductor->getNodeOut()), Voltages)) /
            resistorI +
        currentI);

    transInductor->setVoltage(
        getVoltage(equation->getNodeIndex(transInductor->getNodeIn()), Voltages) -
        getVoltage(equation->getNodeIndex(transInductor->getNodeOut()), Voltages));
}

// update the std::map<const NonlinearComponent *, std::vector<BasicComponent *>>
void AnalyzeTrans::updateInductor(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                                  Equation *equation, TransInductor *transInductor)
{

    // ...
    getParametersI(transInductor);
    Inductor *temp = (Inductor *)component;

    // correspond to Gl
    BasicComponent *rl = new Resistor(temp->getNodeIn(), temp->getNodeOut(), 1 / parameters.back());
    parameters.pop_back();

    // caution: the direction of currentsource
    BasicComponent *currentl = new CurrentSource(temp->getNodeOut(), temp->getNodeIn(), parameters.back());
    parameters.pop_back();

    basicComponents.clear();
    basicComponents.push_back(rl);
    basicComponents.push_back(currentl);
}

void AnalyzeTrans::updateCapacitor(const NonlinearComponent *component,
                                   std::vector<BasicComponent *> &basicComponents, Equation *equation,
                                   TransCapacitor *transCapacitor, int newNode)
{

    // ...
    getParametersC(transCapacitor);
    Capacitor *temp = (Capacitor *)component;

    // NOTICE: CREATED NEW NODE
    BasicComponent *rc = new Resistor(temp->getNodeIn(), newNode, parameters.back());
    parameters.pop_back();

    BasicComponent *vc = new VoltageSource(newNode, temp->getNodeOut(), parameters.back());
    parameters.pop_back();

    basicComponents.clear();
    basicComponents.push_back(rc);
    basicComponents.push_back(vc);
}

// insert the parasitic capacitance
void AnalyzeTrans::insertCapacitor(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                                   Equation *equation, CompMapType *compMap)
{
    if (component->getType() == enumMOSFET)
    {
        MOSFET *temp = (MOSFET *)component;
        std::map<std::string, float> tempMap = (*compMap)[temp->getID()];
        parameters.push_back(0.5 * tempMap["COX"] * temp->getWidth() * temp->getLength());
        parameters.push_back(tempMap["CJ0"]);

        float test = parameters[1];

        Capacitor *capacitor1 = new Capacitor(temp->getNodeGate(), temp->getNodeOut(), 1e-17);

        equation->transformedComponents.insert(
            std::make_pair((const NonlinearComponent *)capacitor1, std::vector<BasicComponent *>()));

        Capacitor *capacitor2 = new Capacitor(temp->getNodeGate(), temp->getNodeIn(), 1e-17);

        equation->transformedComponents.insert(
            std::make_pair((const NonlinearComponent *)capacitor2, std::vector<BasicComponent *>()));
        Capacitor *capacitor3 = new Capacitor(temp->getNodeOut(), 0, 1e-14);

        equation->transformedComponents.insert(
            std::make_pair((const NonlinearComponent *)capacitor3, std::vector<BasicComponent *>()));
        Capacitor *capacitor4 = new Capacitor(temp->getNodeIn(), 0, 1e-14);

        equation->transformedComponents.insert(
            std::make_pair((const NonlinearComponent *)capacitor4, std::vector<BasicComponent *>()));
        parameters.clear();
    }
}
