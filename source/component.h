// Author: JefferyLi0903

#ifndef COMPONENT_H_
#define COMPONENT_H_
#include <cmath>
constexpr double pi = 3.14159265358979323846;

// Components: 
// - BasicComponent : Resistor, VoltageSource, CurrentSource, VCCS, VCVS, CCCS, CCVS
// - NonlinearComponent : Capacitor, Inductor, Diode, BJT, MOSFET

enum ComponentType {
    enumResistor,
    enumVoltageSource,
    enumCurrentSource,
    enumVCCS,
    enumVCVS,
    enumCCCS,
    enumCCVS,
    enumCapacitor,
    enumInductor,
    enumDiode,
    enumBJT,
    enumMOSFET,
    enumACVoltageSource
};

class Component {
public:
    Component(int nodeIn, int nodeOut, enum ComponentType type);
    enum ComponentType getType() const;
    int getNodeIn() const { return nodeIn; }
    int getNodeOut() const { return nodeOut; }
    bool isBasicComponent;
protected:
    const int nodeIn, nodeOut;
    enum ComponentType type;
};

class BasicComponent : public Component {
public:
    BasicComponent(int nodeIn, int nodeOut, float value, enum ComponentType type)
        : Component(nodeIn, nodeOut, type), value(value) {};
    void setValue(const float value);
    float getValue() const;
private:
    float value;
};

class Resistor : public BasicComponent {
public:
    Resistor(int nodeIn, int nodeOut, float value)
        : BasicComponent(nodeIn, nodeOut, value, enumResistor) {};
};

class VoltageSource : public BasicComponent {
public:
    VoltageSource(int nodeIn, int nodeOut, float value)
        : BasicComponent(nodeIn, nodeOut, value, enumVoltageSource) {};
};

class CurrentSource : public BasicComponent {
public:
    CurrentSource(int nodeIn, int nodeOut, float value)
        : BasicComponent(nodeIn, nodeOut, value, enumCurrentSource) {};
};

class VCCS : public BasicComponent {
public:
    VCCS(int nodeIn, int nodeOut, float value, int nodeVpos, int nodeVneg)
        : BasicComponent(nodeIn, nodeOut, value, enumVCCS), nodeVpos(nodeVpos), nodeVneg(nodeVneg) {};
    const int nodeVpos;
    const int nodeVneg;
};

class VCVS : public BasicComponent {
public:
    VCVS(int nodeIn, int nodeOut, float value, int nodeVpos, int nodeVneg)
        : BasicComponent(nodeIn, nodeOut, value, enumVCVS), nodeVpos(nodeVpos), nodeVneg(nodeVneg) {};
    const int nodeVpos; 
    const int nodeVneg;
};

class CCCS : public BasicComponent {
public:
    CCCS(int nodeIn, int nodeOut, float value, const Component& comp)
        : BasicComponent(nodeIn, nodeOut, value, enumCCCS), comp(comp) {};
    const Component& comp;
};

class CCVS : public BasicComponent {
public:
    CCVS(int nodeIn, int nodeOut, float value, const Component& comp)
        : BasicComponent(nodeIn, nodeOut, value, enumCCVS), comp(comp) {};
    const Component& comp;
};

class NonlinearComponent : public Component {
public:
    NonlinearComponent(int nodeIn, int nodeOut, enum ComponentType type)
        : Component(nodeIn, nodeOut, type) {};
};

// type : true = N-type, false = P-type
class MOSFET : public NonlinearComponent {
public:
    MOSFET(int nodeIn, int nodeOut, int nodeGate, bool type, float width, float length, int id)
        : NonlinearComponent(nodeIn, nodeOut, enumMOSFET), nodeGate(nodeGate), type(type), width(width), length(length), id(id) {};
    inline int getNodeGate() {
        return nodeGate;
    }
    inline bool getType() {
        return type;
    }
    inline float getWidth() {
        return width;
    }
    inline float getLength() {
        return length;
    }
    inline int getID() {
        return id;
    }
private:
    int nodeGate;
    bool type;
    float width;
    float length;
    int id;
};

class Capacitor : public NonlinearComponent {
public:
    Capacitor(int nodeIn, int nodeOut, float value)
        : NonlinearComponent(nodeIn, nodeOut, enumCapacitor), value(value) {};
    float getValue() const;
private:
    float value;
};

class Inductor : public NonlinearComponent {
public:
    Inductor(int nodeIn, int nodeOut, float value)
        : NonlinearComponent(nodeIn, nodeOut, enumInductor), value(value) {};
    float getValue() const;
private:
    float value;
};

// timeDelay: nanoseconds
class ACVoltageSource : public NonlinearComponent {
public:
    ACVoltageSource(int nodeIn, int nodeOut, float DCMagnitude, float ACMagnitude, float freq, float timeDelay)
        : NonlinearComponent(nodeIn, nodeOut, enumACVoltageSource), DCMagnitude(DCMagnitude), ACMagnitude(ACMagnitude), freq(freq), timeDelay(timeDelay) {};
    float getValue(float time) const {
        return DCMagnitude + ACMagnitude * sin(2 * pi * freq * (time - timeDelay * 1e-9));
    }
    float getDCMagnitude () const {
        return DCMagnitude;
    }
    float getACMagnitude () const {
        return ACMagnitude;
    }
    float getFreq () const {
        return freq;
    }
    float getTimeDelay () const {
        return timeDelay;
    }
private:
    float DCMagnitude;
    float ACMagnitude;
    float freq;
    float timeDelay;
};

#endif // COMPONENT_H_
