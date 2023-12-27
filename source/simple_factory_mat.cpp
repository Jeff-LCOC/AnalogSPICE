// Author: ShaoqunLi

#include "simple_factory_mat.h" 

Matrix* SimpleFactoryMat::getMatrix(std::string inputString) {
    if (inputString == "SparseMatrix") {
        return new SparseMatrix();
    }
    return nullptr;
}
