// Author: Rad1ance

#include "analyze_dc.h"

AnalyzeDC::AnalyzeDC(Equation *equation, CompMapType *compMap, std::set<int> *nodeList) {
//    std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it3;
//    std::map<int, std::map<std::string, float>> it2 = *compMap;
//    for (it3 = equation->transformedComponents.begin(); it3 != equation->transformedComponents.end(); it3++) {
//        if (it3->first->getType() == enumMOSFET) {
//            MOSFET *temp = (MOSFET *) it3->first;
//            Capacitor *capacitor = new Capacitor( temp->getNodeIn(), temp->getNodeOut(), 2);
//            equation->transformedComponents.insert(std::make_pair((const NonlinearComponent *)capacitor, std::vector<BasicComponent *>()));
//        }
//    }
    std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it1;
    for (it1 = equation->transformedComponents.begin(); it1 != equation->transformedComponents.end(); it1++) {
        if (it1->first->getType() == enumInductor) {
            VoltageSource *voltageSource = new VoltageSource(it1->first->getNodeIn(), it1->first->getNodeOut(), 0);
            it1->second.push_back(voltageSource);
        }
    }

    std::vector<float> InitialVoltages;
    NonlinearSolver NonlinearSolver(equation, compMap, nodeList, InitialVoltages);
    NonlinearSolver.outputResults(nodeList);
}
