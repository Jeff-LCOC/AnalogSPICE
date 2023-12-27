// Author: JefferyLi0903

#include <iostream>
#include <string>
#include "analog_solver.h"

int main(int argc, char** argv) {

// read the first argument and pass it to construct an AnalogSolver Object
    if (argc > 1) {
        std::string filename = argv[1];
        AnalogSolver solver(filename);
        solver.run();
    }
    else {
        std::cout << "Please provide a file name" << std::endl;
    }

    return 0;
}
