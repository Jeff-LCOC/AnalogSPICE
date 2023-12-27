// Author: Rad1ance

#include "analog_solver.h" 
#include "analyze_trans.h"
#include <ctime>

AnalogSolver::AnalogSolver(std::string filename) : filename(filename),
                                                   parser(filename, &componentMap, components, nodeList, outputList) {
    parser.run();
}

void AnalogSolver::run() {
    Equation equation(&components, &outputList, &nodeList);
    // test
    std::cout << "test started" << std::endl;
    int modeSelect;
    std::cout << "Please choose the mode, 1 for DC Analyze and 2 for Transient Analyze" << std::endl;
    std::cin >> modeSelect;


    //file operation
    //delete the postfix 
    size_t lastDot = filename.find_last_of(".");

    if (lastDot != std::string::npos){
        outfilename = filename.substr(0, lastDot) + ".lis";
    }
    else
        outfilename = filename + ".lis";

    std::ofstream outFile(outfilename);

    if (!outFile.is_open()) {
        std::cerr << "Failed to open the file！" << std::endl;
    }

    //considering the memory leak, the while(1) is not provided
    if (modeSelect == 1) {
        startTime = clock();

        // 将输出重定向到文件流
        std::streambuf *coutbuf = std::cout.rdbuf(); 
        std::cout.rdbuf(outFile.rdbuf()); 

        std::cout << "DC analysis started..." << std::endl;
        AnalyzeDC analyzeDC(&equation, &componentMap, &nodeList);
        endTime = clock();
        std::cout << "Running time = " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
        
        std::cout.rdbuf(coutbuf); //reset cout 
    } else if (modeSelect == 2) {
        float deltaT, tFinal;
        std::cout << "Please input the deltaT:" << std::endl;
        std::cin >> deltaT;
        std::cout << "Please input the tFinal:" << std::endl;
        std::cin >> tFinal;
        startTime = clock();

        // 将输出重定向到文件流
        std::streambuf *coutbuf = std::cout.rdbuf(); 
        std::cout.rdbuf(outFile.rdbuf()); 

        std::cout << "Transient analysis started..." << std::endl;
        
        AnalyzeTrans analyzeTrans(&equation, &componentMap, &nodeList, deltaT, tFinal);
        endTime = clock();
        std::cout << "Running time = " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;

        std::cout.rdbuf(coutbuf); //reset cout 
    } else
        std::cout << "Invalid input" << std::endl;
    
    // std::cout << "Running time = " << (double) (endTime - startTime) / CLOCKS_PER_SEC << "s" << std::endl;
    
}

