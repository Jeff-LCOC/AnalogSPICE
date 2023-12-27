// Author: ShaoqunLi

#ifndef MATRIX_H_
#define MATRIX_H_

#include <algorithm>
#include <cmath> // for fabs()
#include <vector>
#include "bipartite_graph.h"

enum class DecompConfig {
    DecompBasic,
    DecompFast,
    DecompPrecise
};

enum class RelaxOption {
    GJ
};

struct ChangeRowIdx {
    int row1;
    int row2;
};

class Matrix {
public:
    virtual ~Matrix() {};
    virtual int getNumRows() const = 0;
    virtual int getNumCols() const = 0;
    virtual float getElement(int row, int col) const = 0;
    virtual void setElement(int row, int col, float value) = 0;     // not safe, try not to use this public method
    virtual void updateElement(int row, int col, float addend) = 0; // safe
    virtual void clear() = 0;
    virtual void copyMat(const Matrix& src) = 0; // copy from source matrix
    virtual bool isSameMat(const Matrix& other) = 0; // too complex, not recommended
    virtual void changeRow(int row1, int row2) = 0;

    virtual void decompLU(Matrix& lower, Matrix& upper, std::vector<ChangeRowIdx>& rec, DecompConfig config = DecompConfig::DecompPrecise) = 0;
    virtual void forwardSub(std::vector<float>& b) = 0;  // for lower matrix(special)
    virtual void backwardSub(std::vector<float>& c) = 0; // for upper matrix
    virtual void normalizeByJ(std::vector<float>& originJ) = 0; // for more-precise decomposition

    virtual void makeDiagNonZero(std::vector<float>& b) = 0;
    virtual void relaxIterate(std::vector<float>& vecX, const std::vector<float>& b, RelaxOption config = RelaxOption::GJ) = 0;
};

struct Element {
    int row;
    int col;
    float value;
};

// for LU decomposition
struct luRecord {
    int index;
    float value;
};

bool compElement (const Element& a, const Element& b);

class SparseMatrix : public Matrix {
public:
    SparseMatrix() : numRows(0), numCols(0), sorted(true) {};
    ~SparseMatrix() {};
    int getNumRows() const override { return numRows; };
    int getNumCols() const override { return numCols; };
    float getElement(int row, int col) const override;
    void setElement(int row, int col, float value) override;
    void updateElement(int row, int col, float addend) override;
    void clear() override;
    void copyMat(const Matrix& src) override;
    bool isSameMat(const Matrix& other) override;
    void changeRow(int row1, int row2) override;

    void decompLU(Matrix& lower, Matrix& upper, std::vector<ChangeRowIdx>& rec, DecompConfig config = DecompConfig::DecompPrecise) override;
    void forwardSub(std::vector<float>& b) override; // make sure the diagonal elements of the lower matrix are all 1
    void backwardSub(std::vector<float>& c) override;
    void normalizeByJ(std::vector<float>& originJ) override;

    void makeDiagNonZero(std::vector<float>& b) override;
    void relaxIterate(std::vector<float>& vecX, const std::vector<float>& b, RelaxOption config = RelaxOption::GJ) override;
private:
    void sortMat();
    bool getElementSorted(int row, int col) const;
    bool getElementSorted(int row, int col, float& value) const;
    void getBipartiteGraph(BipartiteGraph& biGraph) const;
    void changeRowByHungarian(const std::vector<int>& result);
    void makeMatSorted() { if (!sorted) sortMat(); }
    
    std::vector<Element> elements;
    int numRows;
    int numCols;

    bool sorted;
};

#endif // MATRIX_H_
