// Author: ShaoqunLi

#include "direct_linear_solver.h"

DirectLinearSolver::DirectLinearSolver() : prePriority(Priority::None) {
    std::string matType = "SparseMatrix"; // temp
    matACalc = factory.getMatrix(matType);
    // matACopy = factory.getMatrix(matType);
    matLower = factory.getMatrix(matType);
    matUpper = factory.getMatrix(matType);
}

DirectLinearSolver::~DirectLinearSolver() {
    delete matACalc;
    // delete matACopy;
    delete matLower;
    delete matUpper;
}

void DirectLinearSolver::solve(std::vector<float>& solution, const Equation& eq, Priority priority, bool matChanged) {
    solution.clear();

    // temp for this version
    if (priority == Priority::Basic) {
        if (matChanged || prePriority == Priority::None || prePriority == Priority::MorePrecise) {
            matACalc->copyMat(*(eq.getA())); // will be changed later
            // matACopy->copyMat(*(eq.getA())); // not changed
            matACalc->decompLU(*matLower, *matUpper, changeRowRec, DecompConfig::DecompBasic); // all inputs will be cleared before solving
        }

        // get original vector b
        for (int i = 0; i < eq.getASize(); i++) {
            solution.push_back(eq.getb()->getElement(i, 0));
        }

        // generate transformed vector b(row transformation)
        for (int i = 0; i < changeRowRec.size(); i++) {
            int row1 = changeRowRec[i].row1;
            int row2 = changeRowRec[i].row2;
            float temp = solution[row1];
            solution[row1] = solution[row2];
            solution[row2] = temp;
        }

        // forward substitution and backward substitution
        // if matLower and matUpper are SparseMatrix objects, they will not be changed
        matLower->forwardSub(solution);
        matUpper->backwardSub(solution);

        prePriority = Priority::Basic;
    }
    else if (priority == Priority::Fast) {
        if (matChanged || prePriority == Priority::None || prePriority == Priority::MorePrecise) {
            matACalc->copyMat(*(eq.getA())); // will be changed later
            // matACopy->copyMat(*(eq.getA())); // not changed
            matACalc->decompLU(*matLower, *matUpper, changeRowRec, DecompConfig::DecompFast); // all inputs will be cleared before solving
        }

        // get original vector b
        for (int i = 0; i < eq.getASize(); i++) {
            solution.push_back(eq.getb()->getElement(i, 0));
        }

        // generate transformed vector b(row transformation)
        for (int i = 0; i < changeRowRec.size(); i++) {
            int row1 = changeRowRec[i].row1;
            int row2 = changeRowRec[i].row2;
            float temp = solution[row1];
            solution[row1] = solution[row2];
            solution[row2] = temp;
        }

        // forward substitution and backward substitution
        // if matLower and matUpper are SparseMatrix objects, they will not be changed
        matLower->forwardSub(solution);
        matUpper->backwardSub(solution);

        prePriority = Priority::Fast;
    }
    else if (priority == Priority::Precise) {
        if (matChanged || prePriority == Priority::None  || prePriority == Priority::Basic || prePriority == Priority::Fast || prePriority == Priority::MorePrecise) {
            matACalc->copyMat(*(eq.getA())); // will be changed later
            // matACopy->copyMat(*(eq.getA())); // not changed
            matACalc->decompLU(*matLower, *matUpper, changeRowRec, DecompConfig::DecompPrecise);
        }

        // get original vector b
        for (int i = 0; i < eq.getASize(); i++) {
            solution.push_back(eq.getb()->getElement(i, 0));
        }

        // generate transformed vector b(row transformation)
        for (int i = 0; i < changeRowRec.size(); i++) {
            int row1 = changeRowRec[i].row1;
            int row2 = changeRowRec[i].row2;
            float temp = solution[row1];
            solution[row1] = solution[row2];
            solution[row2] = temp;
        }

        // forward substitution and backward substitution
        // if matLower and matUpper are SparseMatrix objects, they will not be changed
        matLower->forwardSub(solution);
        matUpper->backwardSub(solution);

        prePriority = Priority::Precise;
    }
    else if(priority == Priority::MorePrecise) {
        matACalc->copyMat(*(eq.getA())); // matACalc will be changed later

        // get original vector b
        for (int i = 0; i < eq.getASize(); i++) {
            solution.push_back(eq.getb()->getElement(i, 0));
        }
        // normalization
        matACalc->normalizeByJ(solution);

        matACalc->decompLU(*matLower, *matUpper, changeRowRec, DecompConfig::DecompPrecise);

        // generate transformed vector b(row transformation)
        for (int i = 0; i < changeRowRec.size(); i++) {
            int row1 = changeRowRec[i].row1;
            int row2 = changeRowRec[i].row2;
            float temp = solution[row1];
            solution[row1] = solution[row2];
            solution[row2] = temp;
        }

        matLower->forwardSub(solution);
        matUpper->backwardSub(solution);

        prePriority = Priority::MorePrecise;
    }

}

void DirectLinearSolver::setPriorityNone() {
    prePriority = Priority::None;
}
