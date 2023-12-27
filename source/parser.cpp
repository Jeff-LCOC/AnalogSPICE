// Author: JefferyLi0903

#include <cassert>
#include "parser.h"

Parser::Parser(std::string filename, CompMapType *compMap, std::vector<Component *> &components,
               std::set<int> &nodeList, std::vector<Variable *> &outputList)
        : compMap(compMap), components(components), nodeList(nodeList), outputList(outputList) {
    infile.open(filename);
}

void Parser::run() {
    assert(components.empty());
    std::string inputBuffer;
    std::vector<std::string> lineString;
    std::vector<std::string> compList;
    int lineNum(0);
    bool exitLoop = false;
    while (!exitLoop && getline(infile, inputBuffer)) {
        if (inputBuffer.empty()) {
            continue;
        }
        lineString = split(inputBuffer);
        switch (lineString[0][0]) {
            case '*':
                continue;
            case '.':
                if (lineString[0] == ".end" || lineString[0] == ".END") {
                    exitLoop = true;
                }
                else if (lineString[0] == ".model" || lineString[0] == ".MODEL") {
                    std::map<std::string, float> tempMap;
                    for (int i = 2; i < lineString.size(); i += 2) {
                        tempMap[lineString[i]] = std::stof(lineString[i + 1]);
                    }
                    compMap->insert(std::pair<int, std::map<std::string, float> >(std::stoi(lineString[1]), tempMap));
                }
                else if (lineString[0] == ".PLOTBI") {
                    std::vector<std::string>::iterator it = std::find(compList.begin(), compList.end(), lineString[1]);
                    if (it == compList.end()) {
                        std::cerr << "No Such Component! " << std::endl;
                    } else {
                        outputList.push_back(new Current(0, (*components[std::distance(compList.begin(), it)])));
                    }
                } else if (lineString[0] == ".PLOTNV") {
                    std::set<int>::iterator it = nodeList.find(std::stoi(lineString[1]));
                    if (it == nodeList.end()) {
                        std::cerr << "No Such Node! " << std::endl;
                    } else {
                        outputList.push_back(new Voltage(0, std::stoi(lineString[1])));
                    }
                } else {
                    std::cout << "Error: Unknown command" << std::endl;
                    exitLoop = true;
                }
                continue;
            case 'R':
                components.push_back(
                        new Resistor(std::stoi(lineString[1]), std::stoi(lineString[2]), std::stof(lineString[3])));
                nodeList.insert(std::stoi(lineString[1]));
                nodeList.insert(std::stoi(lineString[2]));
                continue;
            case 'V':
                if (lineString[3] == "DC") {
                    components.push_back(new VoltageSource(std::stoi(lineString[1]), std::stoi(lineString[2]),
                                                           std::stof(lineString[4])));
                    nodeList.insert(std::stoi(lineString[1]));
                    nodeList.insert(std::stoi(lineString[2]));
                } else if (lineString[3] == "SIN") {
                    components.push_back(new ACVoltageSource(std::stoi(lineString[1]), std::stoi(lineString[2]),
                                                             std::stof(lineString[4]), std::stof(lineString[5]),
                                                             std::stof(lineString[6]), std::stof(lineString[7])));
                    nodeList.insert(std::stoi(lineString[1]));
                    nodeList.insert(std::stoi(lineString[2]));
                }
                continue;
            case 'I':
                if (lineString[3] == "DC") {
                    components.push_back(new CurrentSource(std::stoi(lineString[1]), std::stoi(lineString[2]),
                                                           std::stof(lineString[4])));
                    nodeList.insert(std::stoi(lineString[1]));
                    nodeList.insert(std::stoi(lineString[2]));
                } else if (lineString[3] == "SIN") {
                    // TODO
                }
                continue;
            case 'M':
                components.push_back(
                        new MOSFET(std::stoi(lineString[1]), std::stoi(lineString[3]), std::stoi(lineString[2]),
                                   lineString[4] == "n" ? true : false, std::stof(lineString[5]),
                                   std::stof(lineString[6]), std::stoi(lineString[7])));
                nodeList.insert(std::stoi(lineString[1]));
                nodeList.insert(std::stoi(lineString[3]));
                nodeList.insert(std::stoi(lineString[2]));
                continue;
            case 'C':
                components.push_back(
                        new Capacitor(std::stoi(lineString[1]), std::stoi(lineString[2]), std::stof(lineString[3])));
                nodeList.insert(std::stoi(lineString[1]));
                nodeList.insert(std::stoi(lineString[2]));
                continue;
            case 'L':
                components.push_back(
                        new Inductor(std::stoi(lineString[1]), std::stoi(lineString[2]), std::stof(lineString[3])));
                nodeList.insert(std::stoi(lineString[1]));
                nodeList.insert(std::stoi(lineString[2]));
                continue;
        }
    }
}

// if (inputBuffer == "R") {
//     int nodeIn, nodeOut;
//     float value;
//     infile >> nodeIn >> nodeOut >> value;
//     components.push_back(Resistor(nodeIn, nodeOut, value));
// }

std::vector<std::string> split(std::string str) {
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream.
    std::string tok;

    while (ss >> tok) {
        internal.push_back(tok);
    }

    return internal;
}

Parser::~Parser() {
    infile.close();
    for (auto comp: components) {
        delete comp;
    }
    for (auto item: outputList) {
        delete item;
    }
}
