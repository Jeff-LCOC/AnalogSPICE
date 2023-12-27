// Author: ShaoqunLi

#ifndef RELAX_LINEAR_SOLVER_H_
#define RELAX_LINEAR_SOLVER_H_

#include <cmath>
#include <iostream>
#include <vector>
#include "equation.h"
#include "matrix.h"
#include "simple_factory_mat.h"

class RelaxLinearSolver {
public:
    RelaxLinearSolver();
    ~RelaxLinearSolver();
    void solve(std::vector<float>& solution, const Equation& eq, RelaxOption option = RelaxOption::GJ);
private:
    void getDelta(const std::vector<float>& x1, const std::vector<float>& x2, float& delta) const;
    SimpleFactoryMat factory;
    Matrix* matACalc;
    std::vector<float> vecb;
};

#endif // RELAX_LINEAR_SOLVER_H_
