// Original Author: Rad1ance
// This version: Jeff Luo

#ifndef NONLINEAR_SOLVER_H_
#define NONLINEAR_SOLVER_H_

#include "equation.h"
#include "linear_solver.h"
#include "parser.h"

class nonlinearMOS {
public:
    //the data of last iteration
    nonlinearMOS(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda, bool type) {
        this->Vg = Vg;
        this->Vd = Vd;
        this->Vs = Vs;
        this->Vt = Vt;
        this->W_L = W_L;
        this->Cox = Cox;
        this->Mu = Mu;
        this->Lamda = Lamda;
        this->type = type;
    }

    float getVg() const;

    float getVd() const;

    float getVs() const;

    float getVt() const;

    float getW_L() const;

    float getCox() const;

    float getMu() const;

    float getLamda() const;

    bool getType() const;

    void setVg(float Vg);

    void setVd(float Vd);

    void setVs(float Vs);

    void setVt(float Vt);

    void setW_L(float W_L);

    void setCox(float Cox);

    void setMu(float Mu);

    void setLamda(float Lamda);

    float getIds() const;

    float getGm() const;

    float getGds() const;

    float getIeq() const;

    void setLinearIds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setLinearGm(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setLinearGds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setLinearIeq(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setCutoffIds();

    void setCutoffGm();

    void setCutoffGds();

    void setCutoffIeq();

    void setSaturationIds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setSaturationGm(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setSaturationGds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

    void setSaturationIeq(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda);

private:
    float Vg;
    float Vd;
    float Vs;
    float Vt;
    float W_L;
    float Cox;
    float Mu;
    float Lamda;
    bool type;
    float Ids;
    float Gm;
    float Gds;
    float Ieq;
};

inline bool nonlinearMOS::getType() const {
    return type;
}

inline void
nonlinearMOS::setLinearIds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Ids = Mu * Cox * W_L * ((Vg - Vs - Vt) * (Vd - Vs) - (Vd - Vs) * (Vd - Vs) / 2);
}

inline void
nonlinearMOS::setLinearGm(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Gm = Mu * Cox * W_L * (Vd - Vs);
}

inline void
nonlinearMOS::setLinearGds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Gds = Mu * Cox * W_L * ((Vg - Vs - Vt) - (Vd - Vs));
}

inline void
nonlinearMOS::setLinearIeq(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Ieq = this->Ids - this->Gm * (Vg - Vs) - this->Gds * (Vd - Vs);
}

inline void
nonlinearMOS::setSaturationIds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Ids = Mu * Cox * W_L * (Vg - Vs - Vt) * (Vg - Vs - Vt) * (1 + Lamda * (Vd - Vs)) / 2;
}

inline void
nonlinearMOS::setSaturationGm(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Gm = Mu * Cox * W_L * (Vg - Vs - Vt) * (1 + Lamda * (Vd - Vs));
}

inline void
nonlinearMOS::setSaturationGds(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Gds = Mu * Cox * W_L * (Vg - Vs - Vt) * (Vg - Vs - Vt) * Lamda / 2;
}

inline void
nonlinearMOS::setSaturationIeq(float Vg, float Vd, float Vs, float Vt, float W_L, float Cox, float Mu, float Lamda) {
    this->Ieq = this->Ids - this->Gm * (Vg - Vs) - this->Gds * (Vd - Vs);
}

inline void nonlinearMOS::setCutoffGds() {
    this->Gds = 0;
}

inline void nonlinearMOS::setCutoffGm() {
    this->Gm = 0;
}

inline void nonlinearMOS::setCutoffIds() {
    this->Ids = 0;
}

inline void nonlinearMOS::setCutoffIeq() {
    this->Ieq = 0;
}

inline float nonlinearMOS::getVt() const {
    return Vt;
}

inline float nonlinearMOS::getVd() const {
    return Vd;
}

inline float nonlinearMOS::getVg() const {
    return Vg;
}

inline float nonlinearMOS::getVs() const {
    return Vs;
}

inline float nonlinearMOS::getW_L() const {
    return W_L;
}

inline float nonlinearMOS::getCox() const {
    return Cox;
}

inline float nonlinearMOS::getMu() const {
    return Mu;
}

inline float nonlinearMOS::getLamda() const {
    return Lamda;
}

inline float nonlinearMOS::getIds() const {
    return Ids;
}

inline float nonlinearMOS::getGm() const {
    return Gm;
}

inline float nonlinearMOS::getGds() const {
    return Gds;
}

inline float nonlinearMOS::getIeq() const {
    return Ieq;
}

inline void nonlinearMOS::setVt(float Vt) {
    this->Vt = Vt;
}

inline void nonlinearMOS::setVd(float Vd) {
    this->Vd = Vd;
}

inline void nonlinearMOS::setVg(float Vg) {
    this->Vg = Vg;
}

inline void nonlinearMOS::setVs(float Vs) {
    this->Vs = Vs;
}

inline void nonlinearMOS::setW_L(float W_L) {
    this->W_L = W_L;
}

inline void nonlinearMOS::setCox(float Cox) {
    this->Cox = Cox;
}

inline void nonlinearMOS::setMu(float Mu) {
    this->Mu = Mu;
}

inline void nonlinearMOS::setLamda(float Lamda) {
    this->Lamda = Lamda;
}

class NonlinearSolver {
public:
    // NonlinearSolver(Equation *equation, CompMapType *compMap, std::set<int> *nodeList);
    NonlinearSolver(Equation *equation, CompMapType *compMap, std::set<int> *nodeList, std::vector<float> InitialVoltages);
    
    ~NonlinearSolver() { delete linearSolver; }
    void outputResults(std::set<int> *nodeList);
    std::vector<float> getFinalSolution() const { return finalSolution; }

private:
    void updateMatrix(const NonlinearComponent *component, std::vector<BasicComponent *> &basicComponents,
                      nonlinearMOS *mos);
    void getParameter(nonlinearMOS *mos);
    void updateParameter(const NonlinearComponent *nonlinearComponent, nonlinearMOS *mos, std::vector<float> Voltages,
                         std::set<int> *nodeList);
    void newtonMethod(std::set<int> *nodeList, Equation &equation);
    void checkError(std::vector<float> &Voltages, float &error);

    // void initialValue(MOSFET *temp, nonlinearMOS &mos);
    void initialValue(MOSFET *temp, nonlinearMOS &mos, std::set<int> *nodeList, std::vector<float> Voltages);

    std::vector<float> parameters;
    std::vector<nonlinearMOS *> mosComponents;
    LinearSolver *linearSolver;
    std::vector<float> finalSolution;

};

#endif // NONLINEAR_SOLVER_H_
