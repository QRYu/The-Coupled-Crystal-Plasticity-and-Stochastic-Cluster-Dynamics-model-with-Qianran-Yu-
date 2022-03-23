//Damage.h -- Damage class to hold damage
#pragma once
#include<iostream>
#include<fstream>
#include<sstream>
#include "constant.h"
using namespace std;

class Damage {
private:
    // private data member:
    double DPA_RATE[POINTS];
    double NRT[POINTS];
    double damage[POINTS][CHANNELS];
    double totalRate[POINTS];
    // private functions:
    void readFile();
    void computeDamageZero(const int&);
    void computeDamageOne(const int&);
    void computeDamageTwo(const int&);
    void computeDamageOther(const int&, const int&);
public:
    Damage(); // constructor
    Reaction selectDamage(const int&, double&);
    void display(const int&) const; // display damage in one element
    const double getTotalDamage(const int&) const;
    double getDpaRate(const int&);
};

