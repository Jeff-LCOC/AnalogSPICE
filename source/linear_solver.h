// Original Author: ShaoqunLi
// this version: Jeff Luo

#ifndef LINEAR_SOLVER_H_
#define LINEAR_SOLVER_H_

#include <vector>
#include "equation.h"
#include "direct_linear_solver.h"
#include "relax_linear_solver.h"
#include "simple_factory_mat.h"

// for this version:
// - direct basic        available
// - direct fast         available
// - direct precise      available
// - direct more precise available
// - relax G-J           available, but not always useful
// - (relax G-S)         not available

enum class SolverOption {
    DirectBasic,        // Basic LU decomposition
    DirectFast,         // Markowitz
    DirectPrecise,
    DirectMorePrecise,  // No L and U reused
    RelaxGJ,            // Gauss-Jacobi
    RelaxGS
};

class LinearSolver {
public:
    LinearSolver();
    void solve(std::vector<float>& solution, const Equation& eq, SolverOption option = SolverOption::DirectPrecise, bool matChanged = true);
private:
    DirectLinearSolver directSolver;
    RelaxLinearSolver relaxSolver;
};

#endif // LINEAR_SOLVER_H_
