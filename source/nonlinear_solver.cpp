// Original Author: Rad1ance
// This version: Jeff Luo

#include "nonlinear_solver.h"

NonlinearSolver::NonlinearSolver(Equation *equation, CompMapType *compMap, std::set<int> *nodeList,
                 std::vector<float> InitialVoltages)
{
    linearSolver = new LinearSolver();
    std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it1;
    std::map<int, std::map<std::string, float>> it2 = *compMap;
    int flag = 0; // 数组mosComponent的索引
    for (it1 = equation->transformedComponents.begin(); it1 != equation->transformedComponents.end(); it1++)
    {
        //...
        // initialize std::map<const NonlinearComponent *, std::vector<BasicComponent *>>
        // focused on the MOSFETs
        // in order to implement the first iteration of newton method
        if (it1->first->getType() == enumMOSFET)
        {
            flag++;
            MOSFET *temp = (MOSFET *)it1->first;
            std::map<std::string, float> tempMap = it2[temp->getID()];
            nonlinearMOS *mos = new nonlinearMOS(0, 0.5, 0, tempMap["VT"],
                                                 temp->getWidth() / temp->getLength(), tempMap["COX"], 
                                                 tempMap["MU"], tempMap["LAMBDA"], temp->getType());
            if (!InitialVoltages.empty())
                initialValue(temp, *mos, nodeList, InitialVoltages); //?
            mosComponents.push_back(mos);
            updateMatrix(it1->first, it1->second, mosComponents[flag - 1]);
        }
    }
    equation->equationAddTransformed();
    newtonMethod(nodeList, *equation);

    // outputResults(nodeList);
}

// output
void NonlinearSolver::outputResults(std::set<int> *nodeList)
{
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

// initial the node voltage of the MOSFET
// void NonlinearSolver::initialValue(MOSFET *temp, nonlinearMOS &mos)
// {

// }

void NonlinearSolver::initialValue(MOSFET *temp, nonlinearMOS &mos, std::set<int> *nodeList, std::vector<float> Voltages)
{

    //...
    int ng = std::distance(nodeList->begin(), nodeList->find(temp->getNodeGate())) - 1;
    int nd = std::distance(nodeList->begin(), nodeList->find(temp->getNodeIn())) - 1;
    int ns = std::distance(nodeList->begin(), nodeList->find(temp->getNodeOut())) - 1;

    mos.setVg(Voltages[ng]);
    mos.setVd(Voltages[nd]);
    mos.setVs(Voltages[ns]);
}

// calculate the error parameter of the newton method
void NonlinearSolver::checkError(std::vector<float> &Voltages, float &error)
{
    // check when the Newton Method is over
    error = 0;
    for (size_t i = 0; i < Voltages.size(); i++)
    {
        error += (Voltages[i] - finalSolution[i]) * (Voltages[i] - finalSolution[i]);
    }
    error = sqrt(error);
}

void NonlinearSolver::newtonMethod(std::set<int> *nodeList, Equation &equation)
{
    float error = 1;
    float errorThreshold = 1e-5;
    int iteration = 0;
    while (error > errorThreshold)
    {
        // Calling the linear solver
        // update the MNA Matrix
        // in this .cpp, I call the function updateParameter() and updateMatrix()
        std::vector<float> Voltages;
        
        //notice that matrix has changed
        linearSolver->solve(Voltages, equation, SolverOption::DirectPrecise, 1);
        if (iteration > 0)
            checkError(Voltages, error);


        // 更新参数
        int flag = 0;
        std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it1;
        
        // cannot use 'auto' here. otherwise equation won't be changed when it1 is changed
        for (it1 = equation.transformedComponents.begin(); it1 != equation.transformedComponents.end(); it1++)
        {
            if (it1->first->getType() == enumMOSFET)
            {
                /* code */
                flag++;
                updateParameter(it1->first, mosComponents[flag - 1], Voltages, nodeList);
                updateMatrix(it1->first, it1->second, mosComponents[flag - 1]); 
            }
        }
        // 更新方程
        equation.equationAddTransformed();

        // 更新结果
        finalSolution = Voltages;
        iteration++;
        if (iteration >= 21)
            break; // 设定最大迭代次数
    }
    // std::cout << "After " << iteration << " iterations, Newton Method was ended. " << std::endl;
}

// set the MOSFET Parameters
//  Only the nmos situation is left
//  此处计算得到gm, gds, ieq，并把它们插入向量parameters中
void NonlinearSolver::getParameter(nonlinearMOS *mos)
{
    switch (mos->getType())
    {
    case (true):
    {
        if (mos->getVg() - mos->getVs() <= mos->getVt())
        {
            mos->setCutoffGds();
            mos->setCutoffGm();
            mos->setCutoffIeq();
            mos->setCutoffIds();
            parameters.push_back(mos->getGm());
            parameters.push_back(mos->getGds());
            parameters.push_back(-mos->getIeq()); //in: d, out:s
        }
        else
        {
            if (mos->getVg() - mos->getVt() >= mos->getVd())
            {
                mos->setLinearIds(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                  mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                mos->setLinearGds(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                  mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                mos->setLinearGm(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                 mos->getCox(),
                                 mos->getMu(), mos->getLamda());
                mos->setLinearIeq(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                  mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                parameters.push_back(mos->getGm());
                parameters.push_back(mos->getGds());
                parameters.push_back(-mos->getIeq());
            }
            else if (mos->getVg() - mos->getVt() < mos->getVd())
            {
                mos->setSaturationIds(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                      mos->getCox(),
                                      mos->getMu(), mos->getLamda());
                mos->setSaturationGds(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                      mos->getCox(),
                                      mos->getMu(), mos->getLamda());
                mos->setSaturationGm(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                     mos->getCox(),
                                     mos->getMu(), mos->getLamda());
                mos->setSaturationIeq(mos->getVg(), mos->getVd(), mos->getVs(), mos->getVt(), mos->getW_L(),
                                      mos->getCox(),
                                      mos->getMu(), mos->getLamda());
                parameters.push_back(mos->getGm());
                parameters.push_back(mos->getGds());
                parameters.push_back(-mos->getIeq());
            }
        }
        break;
    }

    case (false):
    {
        if (mos->getVs() - mos->getVg() <= -mos->getVt())
        {
            mos->setCutoffGds();
            mos->setCutoffGm();
            mos->setCutoffIeq();
            mos->setCutoffIds();
            parameters.push_back(-mos->getGm());
            parameters.push_back(-mos->getGds());
            parameters.push_back(mos->getIeq());
        }
        else
        {
            if (mos->getVs() - mos->getVg() - fabs(mos->getVt()) >= mos->getVs()- mos->getVd())
            {
                mos->setLinearIds(mos->getVg(), mos->getVd(), mos->getVs(),
                                  mos->getVt(), mos->getW_L(), mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                mos->setLinearGds(mos->getVg(), mos->getVd(), mos->getVs(),
                                  mos->getVt(), mos->getW_L(), mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                mos->setLinearGm(mos->getVg(), mos->getVd(), mos->getVs(),
                                 mos->getVt(), mos->getW_L(), mos->getCox(),
                                 mos->getMu(), mos->getLamda());
                mos->setLinearIeq(mos->getVg(), mos->getVd(), mos->getVs(),
                                  mos->getVt(), mos->getW_L(), mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                parameters.push_back(-mos->getGm());
                parameters.push_back(-mos->getGds());
                parameters.push_back(mos->getIeq());
            }
            else if (mos->getVs() - mos->getVg() - fabs(mos->getVt()) < mos->getVs()- mos->getVd())
            {
                mos->setSaturationIds(mos->getVg(), mos->getVd(), mos->getVs(),
                                      mos->getVt(), mos->getW_L(), mos->getCox(),
                                      mos->getMu(), mos->getLamda());
                mos->setSaturationGds(mos->getVg(), mos->getVd(), mos->getVs(),
                                      mos->getVt(), mos->getW_L(), mos->getCox(),
                                      mos->getMu(), mos->getLamda());
                mos->setSaturationGm(mos->getVg(), mos->getVd(), mos->getVs(),
                                     mos->getVt(), mos->getW_L(), mos->getCox(),
                                     mos->getMu(), mos->getLamda());
                mos->setSaturationIeq(mos->getVg(), mos->getVd(), mos->getVs(),
                                  mos->getVt(), mos->getW_L(), mos->getCox(),
                                  mos->getMu(), mos->getLamda());
                parameters.push_back(-mos->getGm());
                parameters.push_back(mos->getGds());
                parameters.push_back(mos->getIeq());
            }
        }
    }
    }
}

// update the parameters of the nonlinear components (Gm Gds Ieq) (这里的注释好像和文档的不一样？)
// 文档注释：存储前一次牛顿迭代节点电压
void NonlinearSolver::updateParameter(const NonlinearComponent *nonlinearComponent, nonlinearMOS *mos,
                                      std::vector<float> Voltages, std::set<int> *nodeList)
{
    MOSFET *temp = (MOSFET *)nonlinearComponent;
    int ng = std::distance(nodeList->begin(), nodeList->find(temp->getNodeGate())) - 1;
    int nd = std::distance(nodeList->begin(), nodeList->find(temp->getNodeIn())) - 1;
    int ns = std::distance(nodeList->begin(), nodeList->find(temp->getNodeOut())) - 1;
    if (ng >= 0)
        mos->setVg(Voltages[ng]);
    if (nd >= 0)
        mos->setVd(Voltages[nd]);
    if (ns >= 0)
        mos->setVs(Voltages[ns]);
}

// transmit the linear components to the basic components
// In other words, update the std::pair<const NonlinearComponent*, std::vector<BasicComponent*>>
// 对每一个mosConponent中的元素分别进行操作
void NonlinearSolver::updateMatrix(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                                   nonlinearMOS *mos)
{
    getParameter(mos);

    MOSFET *temp = (MOSFET*)component;
    
    BasicComponent *ieq = new CurrentSource(temp->getNodeIn(), temp->getNodeOut(),
                                            parameters.back());
    parameters.pop_back();

    //cannot push gds directly. Instead, make it resistor (i.e. push 1/gds)
    BasicComponent *rds = new Resistor(temp->getNodeIn(), temp->getNodeOut(),
                                       1/parameters.back());
    parameters.pop_back();

    //notice that gm is not a resistor between D&S. Instead, it's a voltage-control-current source
    BasicComponent *gm_src = new VCCS(temp->getNodeIn(), temp->getNodeOut(), parameters.back(),
                                    temp->getNodeGate(), temp->getNodeOut());
    parameters.pop_back();

    basicComponents.clear();
    basicComponents.push_back(gm_src);
    basicComponents.push_back(rds);
    basicComponents.push_back(ieq);
}
