// Author: Rad1ance

#ifndef ANALOGSPICE_ANALYZE_DC_H_
#define ANALOGSPICE_ANALYZE_DC_H_

#include "equation.h"
#include "linear_solver.h"
#include "nonlinear_solver.h"
#include "parser.h"

class AnalyzeDC {
public:
    AnalyzeDC(Equation *equation, CompMapType *compMap, std::set<int> *nodeList);
//    ~AnalyzeDC(){ delete nonlinearSolver; }

private:
    NonlinearSolver *nonlinearSolver;
};

#endif //ANALOGSPICE_ANALYZE_DC_H_
