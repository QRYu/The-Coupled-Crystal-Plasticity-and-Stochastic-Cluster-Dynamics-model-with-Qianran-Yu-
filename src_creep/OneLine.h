// OneLine.h -- class OneLine
#pragma once
#include<iostream>
#include<unordered_map>
#include"Object.h"

class OneLine {
private:
    // private data member
    int sinks[LEVELS+1][POINTS];
    double diffRToF;            // diffusion rate from jth element to (j-1)th element
    double diffRToB;            // diffusion rate from jth element to (j+1)th latter one
    double sinkR;                 // go to sink reaction rate
    long double dissociationR[LEVELS];
    // dissociation reaction rate
    long double totalRate;             // totalRate of the line;
    std::unordered_map<int64, long double> secondR;  // second reaction rate
    
    // private functions
    void setOneLine(const Object* const, const int&, unordered_map<int64, Object*>&); // build one line
    void computeDiffReaction(const Object* const, const int&);  // compute diffusion rate
    void computeDissReaction(const Object* const, const int&, const int&); // compute dissociation reaction rate
    double computeCombReaction(const Object* const, const Object* const, const int&);  // compute 2nd order reaction rate
    double computeDimensionTerm(const double&, const Object* const, const Object* const, const int&);
    
public:
    OneLine(const Object* const, const int&, unordered_map<int64, Object*>&);
    Reaction selectReaction(const Object* const,int64&, double&);
    void addReaction(const Object* const, const Object* const, const int&); /* add one reaction when a new object is been created */
    void removeReaction(const int64&); /* delete one reaction when an object is deleted */
    void updateReaction(const Object* const, const Object* const, const int&);/* when the number of another object has changed, one rate in this line should be changed */
    void updateLine(const Object* const, const int&, unordered_map<int64, Object*>&); /* when number of this object has changed, rates in this line should be updated */
    void computeSinkReaction(const Object* const, const int&);         // compute absorption reaction rate
    const long double computeTotalRate();
    const double getDiffRateF() const;
    const double getDiffRateB() const;
    void display(const Object* const);
};
