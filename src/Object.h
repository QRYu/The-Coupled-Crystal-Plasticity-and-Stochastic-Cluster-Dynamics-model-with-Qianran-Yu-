// Object.h -- declaration of Object class that help establish an object
#pragma once
#include<iostream>
#include<cmath>
#include"constant.h"
using namespace std;
// some constant values
const double avol = ALATT*ALATT*ALATT / 2;
const double jumped = sqrt(3.0)* ALATT / 2.0;
typedef long long int  int64;

class Object {
private:
    // private data member
    int64 oKey;
    int attributes[LEVELS];
    int number[POINTS];  // number of this object in several elements
    int totalNumber; // number of this object in the whole system
    unsigned int dimensionality;
    double diffusivity;
    long double bind[LEVELS];
    double bindSH ; /* bind energy of H-H at surface */
    double r1, r1e;
    double sinkStrength;
    double DISLOCATION;
    
    // private functions
    void setKey(); /* use attributes to get key */
    void setAttributes(const int64&); /* use key to get attributes */
    void setNumber(); /* set zeros to every element */
    int setDimensionality();
    void computeR1R1e();
    void computeDiffCoeff();
    void computeBindTerm();
    
    void setProperties(const int&, const int&);
public:
    Object(const int64&, const int&, double&,const int& n=1);  /* constructor one, establishing object by key */
    Object(const int*, const int&, double&,const int& n=1);/* constructor two, establishing object by arrtibutes */
    Object(const int64&, double&, const int*); /* designed for restart, add object by knowing numbers in every element*/
    void addNumber(const int&, const int& n=1); /* default, add 1, but can add whatever n */
    void reduceNumber(const int&);
    int signof(const int64&) const;
    double zero(const int&);
    // functions that get access to private data memeber
    int64 getKey() const;  // get access to object key;
    double getDiff() const;  // get access to diffusivity
    int getNumber(const int&) const;  // get access to the number of object in this element
    int getTotalNumber() const;
    int getAttri(const int&) const;   // get access to one attribute
    double getSink() const;  // get access to sink strength
    long double getBind(const int&) const;  // get access to one bind
    double getR1()const;  // get access to r1
    double getR1e() const;  // get access to r1e
    int getDim() const;     // get access to dimensionality
    double getBindSH() const; // get binding energy of surface hydrogen
    void getThreeNumber(const int&, int*) const;// get access to object number (number[POINTS])
    void computeSinks();
    void updateDislocation(double&);
    void display();
};

