// Author: JefferyLi0903

#include <cassert>
#include "equation.h"
#include <iostream>

Voltage::Voltage(float value, int node) : Variable(value, enumVoltage), node(node) {}
Current::Current(float value, const Component &comp) : Variable(value, enumCurrent), comp(comp) {}

Equation::Equation(const std::vector<Component *> *components, std::vector<Variable *> *outputList, std::set<int> *nodeList)
    : outputList(outputList), size(nodeList->size() - 1), nodeList(nodeList), appendedVariables(0)
{
    // For every component in components, if it is basic, push it to BasicComponents
    // else push it to NonBasicComponents
    basicMat = factory.getMatrix("SparseMatrix");
    basicb = factory.getMatrix("SparseMatrix");
    for (int i = 0; i < components->size(); i++)
    {
        const Component *comp = (*components)[i];
        if (comp->isBasicComponent)
        {
            BasicComponent *tempPtr = (BasicComponent *)comp;
            basicComponents.push_back(tempPtr);
            equationAdd(basicMat, basicb, tempPtr);
        }
        else
        {
            transformedComponents.insert(std::pair<const NonlinearComponent *, std::vector<BasicComponent *>>(static_cast<const NonlinearComponent *>(comp), std::vector<BasicComponent *>()));
        }
    }
    finalMat = factory.getMatrix("SparseMatrix");
    finalb = factory.getMatrix("SparseMatrix");
    finalMat->copyMat(*basicMat);
    finalb->copyMat(*basicb);

    // Initialize all memeber variables
    newNodeList = new std::set<int>();
    specialComponents = new std::vector<BasicComponent *>();
}

Equation::~Equation()
{
    delete basicMat;
    delete basicb;
    delete finalMat;
    delete finalb;
    delete newNodeList;

    for (auto pair : transformedComponents)
    {
        for (auto item : pair.second)
        {
            delete item;
        }
    }
    for (auto item : (*specialComponents))
    {
        delete item;
    }
    delete specialComponents;
}

void Equation::equationAddTransformed()
{
    nonSpecialSize = -1; // Initialize nonSpecialSize
    finalMat->copyMat(*basicMat);
    finalb->copyMat(*basicb);
    size = basicMat->getNumRows();
    assert(size == basicMat->getNumCols());
    for (auto item : transformedComponents)
    {
        for (BasicComponent *comp : item.second)
            equationAdd(finalMat, finalb, comp);
    }
    nonSpecialSize = size;
    size = size + newNodeList->size();
    // Handle the Special Components
    for (auto item : *specialComponents)
    {
        equationAdd(finalMat, finalb, item, true);
    }
    //    for(int i = 0; i < size; i++) {
    //        for(int j = 0; j < size; j++) {
    //            std::cout << finalMat->getElement(i, j) << "\t";
    //        }
    //        std::cout << "=" << finalb->getElement(i, 0) << std::endl;
    //    }
    //    std::cout << "----------" <<std::endl;
    // Clear the addedNodes
    specialComponents->clear();
}

int Equation::applyNewNode()
{
    int newNode = -newNodeList->size() - 1;
    newNodeList->insert(newNode);
    return newNode;
}

int Equation::getNodeIndex(int nodeTag)
{
    if (nodeTag >= 0)
        return std::distance(nodeList->begin(), nodeList->find(nodeTag)) - 1;
    else
    {
        assert(nonSpecialSize != -1);
        // WRONG IF OTHER NON-SPECIAL ITEM HAS NOT BEEN HANDLED
        return nonSpecialSize + std::distance(newNodeList->begin(), newNodeList->find(nodeTag));
    }
}

int Equation::getCurrentIndex(VoltageSource *ptr)
{
    // Find ptr in the appendedVariables
    for (int i = 0; i < appendedVariables.size(); i++)
    {
        if (appendedVariables[i]->type == enumCurrent && &(((Current *)appendedVariables[i])->comp) == (Component *)ptr)
            return i + nodeList->size() - 1;
    }
    return -1;
}

void Equation::ACTransform(float time)
{
    std::map<const NonlinearComponent *, std::vector<BasicComponent *>>::iterator it;
    for (it = transformedComponents.begin(); it != transformedComponents.end(); it++)
    {
        if (it->first->getType() == enumACVoltageSource)
        {
            VoltageSource *voltageSource = new VoltageSource(it->first->getNodeIn(), it->first->getNodeOut(), ((ACVoltageSource *)(it->first))->getValue(time));
            it->second.push_back(voltageSource);
        }
    }
}

int Equation::getASize() const
{
    return size;
}

Matrix *Equation::getb() const
{
    return finalb;
}

Matrix *Equation::getA() const
{
    return finalMat;
}

void Equation::equationAdd(Matrix *matPtr, Matrix *b, BasicComponent *ptr, bool isSpecial)
{
    Matrix &mat = (*matPtr);
    //    int n1 = std::distance(nodeList->begin(), nodeList->find(ptr->getNodeIn())) - 1;
    //    int n2 = std::distance(nodeList->begin(), nodeList->find(ptr->getNodeOut())) - 1;

    // Handle Special Component
    if ((ptr->getNodeIn() < 0 || ptr->getNodeOut() < 0) && !isSpecial)
    {
        specialComponents->push_back(ptr);
        return;
    }

    int n1 = getNodeIndex(ptr->getNodeIn());
    int n2 = getNodeIndex(ptr->getNodeOut());
    assert(n1 != n2);

    if (ptr->getType() == enumResistor)
    {
        float G = 1 / ptr->getValue();
        mat.updateElement(n1, n1, G);
        mat.updateElement(n2, n2, G);
        mat.updateElement(n1, n2, -G);
        mat.updateElement(n2, n1, -G);
    }
    else if (ptr->getType() == enumCurrentSource)
    {
        float I = ptr->getValue();
        b->updateElement(n1, 0, I);
        b->updateElement(n2, 0, -I);
    }
    else if (ptr->getType() == enumVoltageSource)
    {
        float V = ptr->getValue();
        appendedVariables.push_back(new Current(0, *ptr));
        mat.setElement(size, n1, 1);
        mat.setElement(size, n2, -1);
        mat.setElement(n1, size, 1);
        mat.setElement(n2, size, -1);
        b->setElement(size, 0, V);
        size++;
    }
    else if (ptr->getType() == enumVCCS)
    {
        int nVpos = std::distance(nodeList->begin(), nodeList->find(((VCCS *)ptr)->nodeVpos)) - 1;
        int nVneg = std::distance(nodeList->begin(), nodeList->find(((VCCS *)ptr)->nodeVneg)) - 1;
        float G_m = ptr->getValue();
        mat.updateElement(n1, nVpos, G_m);
        mat.updateElement(n2, nVpos, -G_m);
        mat.updateElement(n1, nVneg, -G_m);
        mat.updateElement(n2, nVneg, G_m);
    }
    else if (ptr->getType() == enumVCVS)
    {
        int nVpos = std::distance(nodeList->begin(), nodeList->find(((VCCS *)ptr)->nodeVpos)) - 1;
        int nVneg = std::distance(nodeList->begin(), nodeList->find(((VCCS *)ptr)->nodeVneg)) - 1;
        float A_v = ptr->getValue();
        appendedVariables.push_back(new Current(0, *ptr));
        mat.setElement(size, n1, 1);
        mat.setElement(size, n2, -1);
        mat.setElement(n1, size, 1);
        mat.setElement(n2, size, -1);
        mat.setElement(size, nVpos, -A_v);
        mat.setElement(size, nVneg, A_v);
        size++;
    }
    else if (ptr->getType() == enumCCCS)
    { // Need revision
        float beta = ptr->getValue();
        const Component *temp = (Component *)(&(((CCCS *)ptr)->comp));
        int loc = 0;
        if (temp->getType() == enumVoltageSource || temp->getType() == enumVCVS || temp->getType() == enumCCVS)
        {
            for (Variable *item : appendedVariables)
            {
                if (item->type == enumCurrent && &(((Current *)item)->comp) == temp)
                {
                    loc = std::distance(appendedVariables.begin(), std::find(appendedVariables.begin(), appendedVariables.end(), item));
                    loc += nodeList->size();
                    break;
                }
            }
            mat.updateElement(n1, loc, beta);
            mat.updateElement(n2, loc, -beta);
        }
        else if (temp->getType() == enumCurrentSource)
        {
            CurrentSource transFormedComp = CurrentSource(ptr->getNodeIn(), ptr->getNodeOut(), ((CurrentSource *)temp)->getValue() * beta);
            BasicComponent *compPtr = (BasicComponent *)&transFormedComp;
            equationAdd(matPtr, b, compPtr);
        }
        else if (temp->getType() == enumVCCS)
        {
            VCCS transFormedComp = VCCS(ptr->getNodeIn(), ptr->getNodeOut(), ((VCCS *)temp)->nodeVpos, ((VCCS *)temp)->nodeVneg, ((VCCS *)temp)->getValue() * beta);
            BasicComponent *compPtr = (BasicComponent *)&transFormedComp;
            equationAdd(matPtr, b, compPtr);
        }
        else if (temp->getType() == enumCCCS)
        {
            CCCS transFormedComp = CCCS(ptr->getNodeIn(), ptr->getNodeOut(), ((CCCS *)temp)->getValue() * beta, ((CCCS *)temp)->comp);
            BasicComponent *compPtr = (BasicComponent *)&transFormedComp;
            equationAdd(matPtr, b, compPtr);
        }
        else if (temp->getType() == enumResistor)
        {
            BasicComponent *compr = (BasicComponent *)temp;
            int nr1 = compr->getNodeIn() == 0 ? -1 : std::distance(nodeList->begin(), nodeList->find(ptr->getNodeIn())) - 1;
            int nr2 = compr->getNodeOut() == 0 ? -1 : std::distance(nodeList->begin(), nodeList->find(ptr->getNodeOut())) - 1;
            float rG = 1 / compr->getValue();
            mat.setElement(n1, size, beta);
            mat.setElement(n2, size, -beta);
            mat.setElement(size, nr1, rG);
            mat.setElement(size, nr2, -rG);
            size++;
            appendedVariables.push_back(new Current(0, *compr));
        }
    }
    else if (ptr->getType() == enumCCVS)
    { // Need revision
        float z = ptr->getValue();
        const Component *temp = (Component *)(&(((CCCS *)ptr)->comp));
        int loc = 0;
        mat.setElement(n1, size, 1);
        mat.setElement(n2, size, -1);
        if (temp->getType() == enumVoltageSource || temp->getType() == enumVCVS || temp->getType() == enumCCVS)
        {
            mat.setElement(size, n1, 1);
            mat.setElement(size, n2, -1);
            for (Variable *item : appendedVariables)
            {
                if (item->type == enumCurrent && &(((Current *)item)->comp) == temp)
                {
                    loc = std::distance(appendedVariables.begin(), std::find(appendedVariables.begin(), appendedVariables.end(), item)) - 1;
                    loc += nodeList->size();
                    break;
                }
            }
            mat.updateElement(size, loc, -z);
        }
        else if (temp->getType() == enumCurrentSource)
        {
            mat.setElement(size, n1, 1);
            mat.setElement(size, n2, -1);
            b->setElement(size, 0, z);
        }
        else if (temp->getType() == enumCCCS)
        {
            mat.setElement(size, n1, 1);
            mat.setElement(size, n2, -1);
            CCVS transFormedComp = CCVS(ptr->getNodeIn(), ptr->getNodeOut(), ((CCCS *)temp)->getValue() * z, ((CCCS *)temp)->comp);
            BasicComponent *compPtr = (BasicComponent *)&transFormedComp;
            equationAdd(matPtr, b, compPtr);
        }
        else if (temp->getType() == enumVCCS)
        {
            int nVposvccs = std::distance(nodeList->begin(), nodeList->find(((VCCS *)ptr)->nodeVpos)) - 1;
            int nVnegvccs = std::distance(nodeList->begin(), nodeList->find(((VCCS *)ptr)->nodeVneg)) - 1;
            float G_mvccs = ptr->getValue();
            mat.updateElement(size, nVposvccs, G_mvccs);
            mat.updateElement(size, nVposvccs, G_mvccs);
            mat.updateElement(size, size, -z);
        }
    }
}

std::vector<int> Equation::getProbeNode()
{
    if (outputList->empty())
    {
        std::cout << "Please provide the node to be measured (Hint: Using .PLOTNV)" << std::endl;
        exit(0);
    }
    else
    {
        for (int i = 0; i < outputList->size(); i++)
        {
            Voltage *temp = (Voltage*)(*outputList)[i];
            probeNodeList.push_back(temp->node);            
        }
        return probeNodeList;
    }
}