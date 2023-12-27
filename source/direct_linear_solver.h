// Author: ShaoqunLi

#ifndef DIRECT_LINEAR_SOLVER_H_
#define DIRECT_LINEAR_SOLVER_H_

#include <vector>
#include "equation.h"
#include "matrix.h"
#include "simple_factory_mat.h"

// DirectLinearSolver:
// - LU decomposition: implemented
// - Cholesky decomposition(not yet implemented)

enum class Priority {
    Basic,
    Fast,
    Precise,
    MorePrecise,    // L and U can't be reused
    None = -1
};

class DirectLinearSolver {
public:
    DirectLinearSolver();
    ~DirectLinearSolver();
    void solve(std::vector<float>& solution, const Equation& eq, Priority priority = Priority::Precise, bool matChanged = true);
    void setPriorityNone(); // when using relax solver, set perPriority back to None
private:
    SimpleFactoryMat factory;
    Matrix* matACalc;   // for solving, will be changed
    // Matrix* matACopy;   // copy of sorted matrix A
    Matrix* matLower;
    Matrix* matUpper;
    std::vector<ChangeRowIdx> changeRowRec;
    Priority prePriority;
};

#endif // DIRECT_LINEAR_SOLVER_H_
