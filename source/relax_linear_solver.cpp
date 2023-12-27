// Author: ShaoqunLi

#include "relax_linear_solver.h"

RelaxLinearSolver::RelaxLinearSolver() {
    std::string matType = "SparseMatrix"; // temp
    matACalc = factory.getMatrix(matType);
}

RelaxLinearSolver::~RelaxLinearSolver() {
    delete matACalc;
}

void RelaxLinearSolver::solve(std::vector<float>& solution, const Equation& eq, RelaxOption option) {
    // vecb.clear();
    solution.clear();

    matACalc->copyMat(*(eq.getA())); // sorted
    for (int i = 0; i < eq.getASize(); i++) {
        vecb.push_back(eq.getb()->getElement(i, 0));
    }

    matACalc->makeDiagNonZero(vecb); // sorted too

    float e;
    float eThreshold = 1e-10;
    int iteration = 0;
    std::vector<float> newX(eq.getASize(), 0);
    do {
        matACalc->relaxIterate(newX, vecb, option);
        getDelta(solution, newX, e);
        solution = newX;
        iteration++;
        if (iteration > 10000) {
            std::cout << "Relaxation Solver(Gauss-Jacobi) failed." << std::endl;
            break;
        }
    } while (e > eThreshold);
}

void RelaxLinearSolver::getDelta(const std::vector<float>& x1, const std::vector<float>& x2, float& delta) const {
    float norm = 0;
    for (int i = 0; i < x1.size(); i++) {
        norm += (x1[i] - x2[i]) * (x1[i] - x2[i]);
    }
    delta = sqrt(norm);
}
