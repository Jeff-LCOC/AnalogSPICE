// Original Author: Rad1ance
// This version: Jeff Luo

#ifndef ANALOGSPICE_ANALYZE_TRANS_H_
#define ANALOGSPICE_ANALYZE_TRANS_H_

#include "equation.h"
#include "linear_solver.h"
#include "nonlinear_solver.h"
#include "parser.h"
#include "analyze_dc.h"
#include <cmath>

class TransCapacitor {
public:
    TransCapacitor(float value, float current, float voltage, int nodeIn, int nodeOut, int newNode) {
        this->value = value;
        this->current = current;
        this->voltage = voltage;
        this->nodeIn = nodeIn;
        this->nodeOut = nodeOut;
        this->newNode = newNode;
    }

    float getCurrent() const { return current; }

    float getVoltage() const { return voltage; }

    float getValue() const { return value; }

    int getNewNode() const { return newNode; }

    int getNodeIn() const { return nodeIn; }

    int getNodeOut() const { return nodeOut; }

    void setCurrent(float current) { this->current = current; }

    void setVoltage(float voltage) { this->voltage = voltage; }

private:
    float value;
    float current;
    float voltage;
    int nodeIn;
    int nodeOut;
    int newNode;
};

class TransInductor {
public:
    TransInductor(float value, float current, float voltage, int nodeIn, int nodeOut) {
        this->value = value;
        this->current = current;
        this->voltage = voltage;
        this->nodeIn = nodeIn;
        this->nodeOut = nodeOut;
    }

    float getCurrent() const { return current; }

    float getVoltage() const { return voltage; }

    float getValue() const { return value; }

    int getNodeIn() const { return nodeIn; }

    int getNodeOut() const { return nodeOut; }

    void setCurrent(float current) { this->current = current; }

    void setVoltage(float voltage) { this->voltage = voltage; }

private:
    float value;
    float current;
    float voltage;
    int nodeIn;
    int nodeOut;

};

class AnalyzeTrans {
public:
    AnalyzeTrans(Equation *equation, CompMapType *compMap, std::set<int> *nodeList, float deltaT, float tFinal);
    float getVoltageC() const { return voltageC; }
    float getCurrentI() const { return currentI; }
    float getResistorC() const { return resistorC; }
    float getResistorI() const { return resistorI; }
    float getDeltaT() const { return deltaT; }
private:
    void updateInductor(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                        Equation *equation, TransInductor *transInductor);
    void updateCapacitor(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                         Equation *equation, TransCapacitor *transCapacitor, int newNode);
    void newtonMethod(std::set<int> *nodeList, Equation *equation, CompMapType *compMap);
    void insertCapacitor(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                         Equation *equation, CompMapType *compMap);
    void initialTimePace(float deltaT) { this->deltaT = deltaT; }
    void initialFinalTime(float tFinal) { this->tFinal = tFinal; }
    void initialCurrentT() { this->currentT = 0; }
//    void timePaceGen();
    void getParametersC(TransCapacitor *transCapacitor);
    void getParametersI(TransInductor *transInductor);
    void updateVC_C(TransCapacitor *transCapacitor, const std::vector<float> Voltages, std::set<int> *nodeList,
                    Equation *equation);
    void updateVC_I(TransInductor *transInductor, const std::vector<float> Voltages, std::set<int> *nodeList,
                    Equation *equation);
    float getVoltage(int index, const std::vector<float> &Voltages);
    void timeStepControlC(Equation *equation, float tPre, float &tNow, int flag, std::vector<float> Voltages);
    void timeStepControlL(Equation *equation, float tPre, float &tNow, int flag, std::vector<float> Voltages);
    void initialDeltaTUnchanged(float t) { deltaTUnchanged = t ;}

    void initialDeltaTOrigin(float t) {deltaTOrigin = t;}

    std::vector<float> parameters;
    float epsilonL{};
    float epsilon{};
    float tFinal{};
    float deltaT{};
    float currentT{};
    float voltageC{};
    float currentI{};
    float resistorC{};
    float resistorI{};
    float deltaTUnchanged{};
    
    float q_max{};
    float deltaTOrigin{};

    std::vector<int> probeList;
    std::vector<int> probeNameList;

    std::vector<float> finalSolution;
    std::vector<TransCapacitor *> Capacitors;
    std::vector<TransInductor *> Inductors;
    NonlinearSolver *nonlinearSolver{};
    AnalyzeDC *analyzeDC{};
};

#endif //ANALOGSPICE_ANALYZE_TRANS_H_
