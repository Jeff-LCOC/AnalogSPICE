// Author: ShaoqunLi

#include "linear_solver.h"

LinearSolver::LinearSolver()
        : directSolver(), relaxSolver() {}

void LinearSolver::solve(std::vector<float>& solution, const Equation& eq, SolverOption option, bool matChanged) {
    if (option == SolverOption::DirectBasic) {
        directSolver.solve(solution, eq, Priority::Basic, matChanged);
    }
    else if (option == SolverOption::DirectFast) {
        directSolver.solve(solution, eq, Priority::Fast, matChanged);
    }
    else if (option == SolverOption::DirectPrecise) {
        directSolver.solve(solution, eq, Priority::Precise, matChanged);
    }
    else if (option == SolverOption::DirectMorePrecise) {
        directSolver.solve(solution, eq, Priority::MorePrecise);
    }
    else if (option == SolverOption::RelaxGJ) {
        directSolver.setPriorityNone();
        relaxSolver.solve(solution, eq, RelaxOption::GJ);
    }
    else {

    }
};
