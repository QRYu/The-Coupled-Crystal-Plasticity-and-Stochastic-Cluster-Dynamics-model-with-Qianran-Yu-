//
//  slipSystem.hpp
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//

#pragma once
#include <iostream>
#include <stdio.h>
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include <ctime>
#include <unistd.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <math.h>
#include "SCDWrapper.h"
#include "constant.h"
using namespace std;
using namespace Eigen;

class SlipSystem{
private:
    double theta; // [0, pi/4] polar coordinate for standard triangle
    double phi; // [0, 54.7*pi/180.0] polar coordinate for standard triangle
    /* slip system */
    Vector3d slipDirection[12];
    /* slip plane normal */
    Vector3d slipPlaneNormal[12];
    /* slip direction -60 degree of slip plane (non-schmid) */
    //Vector3d slipDirection_60[12];
    Matrix3d sCCn[12]; /* s circle_cross n */
    //Matrix3d sCCn1[12]; /* s circle_cross n1 (non-schmid)*/
    Matrix3d nCsCCn[12]; /* n cross s circle_cross n */
    //Matrix3d n1CsCCn1[12]; /* n1 cross s circle_cross n1 (non-schmid) */
    Matrix3d L; // transformation lattice
    void computeTransfromMatrix(); // coordinate transformation
    void transferSlipSystems();
    void computeSlipSystems();
public:
    SlipSystem(double&, double&);
    Vector3d getSlipDirection(int&);
    Vector3d getSlipPlaneNormal(int&);
    //Vector3d getSlipDirection_60(int&);
    Matrix3d getSCCN(int&);
    //Matrix3d getSCCN1(int&);
    Matrix3d getNCSCCN(int&);
    void updateSlipSystems(double&, double&);
    //Matrix3d getN1CSCCN1(int&);
    ~SlipSystem(){};
};
