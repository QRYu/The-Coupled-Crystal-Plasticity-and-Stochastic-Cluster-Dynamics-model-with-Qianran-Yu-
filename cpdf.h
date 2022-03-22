// cpdf.h -- cpdf class to read cpdf files
#pragma once
#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<unordered_map>
#include"constant.h"
using namespace std;

class Cpdf {
private:
    unordered_map<int, vector<double> > recoilEnergy;
    unordered_map<int, vector<double> > cumulPossibility;
    int size[POINTS];
    double maxPossibility[POINTS];
public:
    Cpdf(); /* constructor */
    double samplePkaEnergy(const double&, const int&) const;
    double getMaxPossibility(const int&);
};

