// Author: JefferyLi0903

#ifndef EQUATION_H_
#define EQUATION_H_

#include <map>
#include <set>
#include <vector>
#include "component.h"
#include "simple_factory_mat.h"

// class EquationInterface {
// public:
//     virtual int getASize() const = 0;
//     virtual float getA(int i, int j) const = 0;
//     virtual float getb(int i) const = 0;
// };

enum VariableType {
    enumVoltage,
    enumCurrent
};

struct Variable {
    Variable(float value, VariableType type) : value(value), type(type) {};
    float value;
    VariableType type;
};

struct Voltage : public Variable {
    Voltage(float value, int node);
    int node;
};

struct Current : public Variable {
    Current(float value, const Component& comp);
    const Component& comp;
};

class Equation {
public:
    Equation(const std::vector<Component*>* components, std::vector<Variable*>* outputList, std::set<int>* nodeList);
    ~Equation();
    int getASize() const;
    Matrix* getA() const;
    Matrix* getb() const;
    void equationAddTransformed();
    std::vector<Variable*> wrap(std::vector<float> solution);
    std::map<const NonlinearComponent*, std::vector<BasicComponent*>> transformedComponents;
    std::set<int>* nodeList;

    std::vector<int> getProbeNode();

    // Add nodes
    int applyNewNode();
    int getNodeIndex(int nodeTag);
    int getCurrentIndex(VoltageSource* ptr);
    // AC Sweep
    void ACTransform(float time);
private:
    int size;
    std::vector<Variable*>* outputList;
    std::vector<const BasicComponent*> basicComponents;
    std::vector<Variable*> appendedVariables;
    Matrix* basicMat;
    Matrix* finalMat;
    Matrix* basicb;
    Matrix* finalb;
    SimpleFactoryMat factory;
    void equationAdd(Matrix* mat, Matrix* b, BasicComponent* ptr, bool isSpecial = false);
    // New nodes set
    int nonSpecialSize;
    std::set<int>* newNodeList; // All negative tags
    std::vector<BasicComponent*>* specialComponents; // All components containing negative tags

    std::vector<int> probeNodeList;
};

#endif // EQUATION_H_