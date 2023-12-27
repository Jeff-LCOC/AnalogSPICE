// Author: ShaoqunLi

#ifndef SIMPLE_FACTORY_MAT_H_
#define SIMPLE_FACTORY_MAT_H_

#include <string>
#include "matrix.h"

class SimpleFactoryMat {
public:
    Matrix* getMatrix(std::string inputString);
};

#endif // SIMPLE_FACTORY_MAT_H_
