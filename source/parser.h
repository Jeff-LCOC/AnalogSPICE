// Author: JefferyLi0903

#ifndef PARSER_H_
#define PARSER_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include "component.h"
#include "equation.h"

typedef std::map<int, std::map<std::string, float>> CompMapType;

class Parser {
public:
    Parser(std::string filename, CompMapType* compMap, std::vector<Component*>& components, std::set<int>& nodeList, std::vector<Variable*>& outputList);
    ~Parser();
    void run();
private:
    std::ifstream infile;
    CompMapType* compMap;
    std::set<int>& nodeList;
    std::vector<Component*>& components;
    std::vector<Variable*>& outputList;
};

std::vector<std::string> split(std::string str);

#endif // PARSER_H_
