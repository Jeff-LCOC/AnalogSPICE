// Author: JefferyLi0903

#ifndef ANALOG_SOLVER_H_
#define ANALOG_SOLVER_H_

#include <map>
#include <string>
#include <vector>
#include "analyze_dc.h"
#include "component.h"
#include "equation.h"
#include "nonlinear_solver.h"
#include "parser.h"

class AnalogSolver {
public:
    AnalogSolver(std::string filename);
    void run();
private:
    std::string filename;
    std::vector<Component*> components;
    std::set<int> nodeList;
    std::vector<Variable*> outputList;
    Parser parser;
    CompMapType componentMap;
    int startTime, endTime;

    std::string outfilename;
};

#endif // ANALOG_SOLVER_H_
