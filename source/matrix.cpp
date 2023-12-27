// Author: ShaoqunLi

#include "matrix.h"

bool compElement (const Element& a, const Element& b) {
    if (a.row < b.row) return true;
    if (a.row > b.row) return false;
    return a.col < b.col;
}

float SparseMatrix::getElement(int row, int col) const {
    for (int i = 0; i < elements.size(); i++) {
        if (elements[i].row == row && elements[i].col == col) {
            return elements[i].value;
        }
    }
    return 0;
}

// public though, use this method carefully!
void SparseMatrix::setElement(int row, int col, float value) {
    if (row < 0 || col < 0) return;
    if (fabs(value) <= std::numeric_limits<float>::epsilon()) return;

    elements.push_back({row, col, value});
    numRows = row < numRows ? numRows : row + 1;
    numCols = col < numCols ? numCols : col + 1;

    sorted = false;
}

void SparseMatrix::updateElement(int row, int col, float addend) {
    if (row < 0 || col < 0) return;

    for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); it++) {
        if (it->row == row && it->col == col) {
            float newValue = it->value + addend;
            if (fabs(newValue) <= std::numeric_limits<float>::epsilon()) {
                elements.erase(it);
                return;
            }
            it->value = newValue;
            return;
        }
    }
    setElement(row, col, addend);
}

void SparseMatrix::clear() {
    elements.clear();
    numRows = 0;
    numCols = 0;
    sorted = true;
}

// after copying, the SparseMatrix is sorted
void SparseMatrix::copyMat(const Matrix& src) {
    clear(); // sorted = true
    numRows = src.getNumRows();
    numCols = src.getNumCols();
    float element = 0;
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            element = src.getElement(i, j);
            if (fabs(element) > std::numeric_limits<float>::epsilon()) {
                elements.push_back({i, j, element});
            }
        }
    }
}

// not recommended now
bool SparseMatrix::isSameMat(const Matrix& other) {
    if (this->getNumRows() != other.getNumRows() || this->getNumCols() != other.getNumCols()) {
        return false;
    }
    for (int i = 0; i < this->getNumRows(); i++) {
        for (int j = 0; j < this->getNumCols(); j++) {
            if (fabs(this->getElement(i, j) - other.getElement(i, j)) >= std::numeric_limits<float>::epsilon()) {
                return false;
            }
        }
    }
    return true;
}

void SparseMatrix::changeRow(int row1, int row2) {
    for (int i = 0; i < elements.size(); i++) {
        if (elements[i].row == row1) {
            elements[i].row = row2;
        }
        else if (elements[i].row == row2) {
            elements[i].row = row1;
        }
    }

    sorted = false;
}

// this matrix will be sorted if not
void SparseMatrix::decompLU(Matrix& lower, Matrix& upper, std::vector<ChangeRowIdx>& rec, DecompConfig config) {
    if (numRows != numCols) return;

    makeMatSorted();

    lower.clear();
    upper.clear();
    rec.clear();
    
    for (int i = 0; i < numRows; i++) {
        float pivot = 0;
        std::vector<luRecord> uRec;
        std::vector<luRecord> lRec;

        // get pivot
        if (config == DecompConfig::DecompBasic) {
            if (!getElementSorted(i, i, pivot)) {
                for (int j = i + 1; j < numRows; j++) {
                    if (getElementSorted(j, i, pivot)) {
                        changeRow(i, j);
                        lower.changeRow(i, j);
                        rec.push_back({i, j});
                        sortMat();  // A is always sorted
                        break;

                        if (j == numRows - 1) {
                            /* error */
                        }
                    }
                }
            }
        }
        else if (config == DecompConfig::DecompFast) {
            // get (NZUR - 1)s
            std::vector<int> nzur_list(numRows - i, 0);
            std::vector<Element>::iterator it = elements.begin();

            // assume that each row has at least one non-zero element
            while (it != elements.end() && it->row != i) {
                it++;
            }
            for (int j = i; j < numRows; j++) {
                while (it != elements.end() && it->col < i) {
                    it++;
                }

                // get NZUR(only when mat(j, i) is non-zero element)
                if (it->col == i) {
                    do {
                        nzur_list[j - i]++;
                        it++;
                    } while (it != elements.end() && it->row == j);
                }
                else {
                    do {
                        it++;
                    } while (it != elements.end() && it->row == j);
                }
            }

            // find the minimum NZUR in list
            {
                int min_nzur_row = i;
                int min_nzur = 0;
                for (int j = 0; j < numRows - i; j++) {
                    if (nzur_list[j] != 0) {
                        if (min_nzur == 0 || min_nzur > nzur_list[j]) {
                            min_nzur_row = i + j;
                            min_nzur = nzur_list[j];
                        }
                    }
                }
                getElementSorted(min_nzur_row, i, pivot);

                if (min_nzur_row != i) {
                    changeRow(i, min_nzur_row);
                    lower.changeRow(i, min_nzur_row);
                    rec.push_back({i, min_nzur_row});
                    sortMat();
                }
            }
        }
        else if (config == DecompConfig::DecompPrecise) {
            int pivot_row = i;
            float temp_pivot = 0;
            getElementSorted(i, i, pivot);
            for (int j = i + 1; j < numRows; j++) {
                if (getElementSorted(j, i, temp_pivot) && fabs(temp_pivot) > fabs(pivot)) {
                    pivot = temp_pivot;
                    pivot_row = j;
                }
            }
            if(pivot_row != i) {
                changeRow(i, pivot_row);
                lower.changeRow(i, pivot_row);
                rec.push_back({i, pivot_row});
                sortMat();  // A is always sorted
            }
        }

        // update upper matrix
        upper.setElement(i, i, pivot);
        for (int j = i + 1; j < numRows; j++) {
            float u = 0;
            if (getElementSorted(i, j, u)) {
                upper.setElement(i, j, u);
                uRec.push_back({j, u});
            }
        }
        
        // update lower matrix
        lower.setElement(i, i, 1);
        for (int j = i + 1; j < numRows; j++) {
            float m = 0;
            if (getElementSorted(j, i, m)) {
                float l = m / pivot;
                lower.setElement(j, i, l);
                lRec.push_back({j, l});
            }
        }

        // update matrix A
        for (std::vector<Element>::iterator it = elements.begin(); it != elements.end(); ) {
            if (it->row == i || it->col == i) {
                it = elements.erase(it);
            } else {
                it++;
            }
        }
        for (int j = 0; j < uRec.size(); j++) {
            for (int k = 0; k < lRec.size(); k++) {
                updateElement(lRec[k].index, uRec[j].index, - lRec[k].value * uRec[j].value);
            }
        }
        sortMat();    
    }

    // note that lower matrix is not sorted if it's a SparseMatrix object
}

// only for lower matrix
// make sure the diagonal elements of the lower matrix are all 1
void SparseMatrix::forwardSub(std::vector<float>& b) {
    makeMatSorted();

    std::vector<Element>::iterator it = elements.begin();

    // for each row of the lower matrix
    for (int i = 0; i < numRows; i++) {
        while (it->col < i) {
            b[i] -= it->value * b[it->col];
            it++;
        }

        // when reach element(i,i), do nothing but it++
        it++;

        if (it == elements.end()) break;
    }
}

// only for upper matrix
void SparseMatrix::backwardSub(std::vector<float>& c) {
    // makeMatSorted();

    std::vector<Element>::iterator it = elements.end() - 1;

    // for each row of the lower matrix
    for (int i = numRows - 1; i > -1; i--) {
        while (it->col > i) {
            c[i] -= it->value * c[it->col];
            it--;
        }

        c[i] /= it->value;
        // it--;

        // if (it + 1 == elements.begin()) break;
        if (it == elements.begin()) break;
        
        it--;
    }
}

// for more-precise solving
// make sure the SparseMatrix is sorted before
void SparseMatrix::normalizeByJ(std::vector<float>& originJ) {
    makeMatSorted();

    std::vector<Element>::iterator it = elements.begin();
    float epsilon = 1e-6; // temp
    for (int i = 0; i < originJ.size(); i++) {
        if (fabs(originJ[i]) > epsilon) {
            while (it != elements.end() && it->row == i) {
                it->value /= originJ[i];
                it++;
            }
            originJ[i] = 1;
        }
        else {
            while (it != elements.end() && it->row == i) {
                it++;
            }
        }
    }
}

// for relaxation solver
// it's better to have the SparseMatrix sorted
// vector b will be transformed too
void SparseMatrix::makeDiagNonZero(std::vector<float>& b) {
    BipartiteGraph bipartiteGraph(numRows);
    getBipartiteGraph(bipartiteGraph);
    if (bipartiteGraph.hungarianForPerfectMatch()) {
        bipartiteGraph.transPreToResult();
        std::vector<int> result = bipartiteGraph.getResult();
        changeRowByHungarian(result); // then not sorted
        // sortMat(); // necessary for relaxIterate()

        // generate transformed vector b(row transformation)
        std::vector<float> originalB(b);
        for (int i = 0; i < numRows; i++) {
            b[result[i]] = originalB[i];
        }
    }
}

// before using this method, make sure that the SparseMatrix is sorted before
void SparseMatrix::relaxIterate(std::vector<float>& vecX, const std::vector<float>& b, RelaxOption config) {
    makeMatSorted();

    if (config == RelaxOption::GJ) {
        std::vector<Element>::iterator it = elements.begin();
        std::vector<float> tempVec(numRows, 0);
        for (int i = 0; i < numRows; i++) {
            float aii = 0; // if aii = 0 in the end, there will be an error
            tempVec[i] = b[i];
            while (it != elements.end() && it->row == i) {
                if (it->col == i) {
                    aii = it->value;
                }
                else {
                    tempVec[i] -= it->value * vecX[it->col];
                }
                it++;
            }
            tempVec[i] /= aii;
        }

        // copy
        for (int i = 0; i < numRows; i++) {
            vecX[i] = tempVec[i];
        }
    }
}

void SparseMatrix::sortMat() {
    sort(elements.begin(), elements.end(), compElement);
    sorted = true;
}

bool SparseMatrix::getElementSorted(int row, int col) const {
    for (int i = 0; i < elements.size(); i++) {
        if (elements[i].row == row) {
            if (elements[i].col == col) {
                return true;
            }
            if (elements[i].col > col) {
                return false;
            }
        }
        if (elements[i].row > row) {
            return false;
        }
    }
    return false;
}

bool SparseMatrix::getElementSorted(int row, int col, float& value) const {
    for (int i = 0; i < elements.size(); i++) {
        if (elements[i].row == row) {
            if (elements[i].col == col) {
                value = elements[i].value;
                return true;
            }
            if (elements[i].col > col) {
                value = 0;
                return false;
            }
        }
        if (elements[i].row > row) {
            value = 0;
            return false;
        }
    }
    value = 0;
    return false;
}

void SparseMatrix::getBipartiteGraph(BipartiteGraph& biGraph) const {
    biGraph.clear();
    for (auto i : elements) {
        biGraph.addEdge(i.row, i.col);
    }
}

void SparseMatrix::changeRowByHungarian(const std::vector<int>& result) {
    for (auto& i : elements) {
        i.row = result[i.row];
    }
    sorted = false;
}
