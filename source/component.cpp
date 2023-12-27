// Author: JefferyLi0903

#include "component.h"

Component::Component(int nodeIn, int nodeOut, enum ComponentType type) : nodeIn(nodeIn), nodeOut(nodeOut), type(type) {
    if (type == enumResistor || type == enumVoltageSource || type == enumCurrentSource || type == enumVCCS || type == enumCCCS || type == enumVCVS || type == enumCCVS)
        isBasicComponent = true;
    else 
        isBasicComponent = false;
}

float BasicComponent::getValue() const {
    return value;
}

void BasicComponent::setValue(const float value) {
    this->value = value;
}

enum ComponentType Component::getType() const {
    return type;
}

float Capacitor::getValue() const {
    return value;
}

float Inductor::getValue() const {
    return value;
}
