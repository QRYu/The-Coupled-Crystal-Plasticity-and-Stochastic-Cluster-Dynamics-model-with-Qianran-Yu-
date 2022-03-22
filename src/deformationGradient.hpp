//
//  deformationGradient.hpp
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//
#pragma once

#include "slipSystem.hpp"
#include<ctime>


typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 6, 6> Matrix6d;

class DeformationGradient{
private:
    SCDWrapper srscd;
    SlipSystem slip_systems;
    /* declare necessary constants in this class -- for tungsten*/
    int maxiter;
    double h; // [cm] kink height
    //double T0; // [K] 0.8*Tm
    double xi[12][12]; //dislocation junction strength matrix
    //double lambda[6]; // dislocation-defect interaction coefficients
    /* declare other parameters */
    //double l; // dislocation length = |121|*ALATT
    double totalStrainRate; //[s^{-1}]total strain rate, assume constant strain rate
    double TOL;
    double t; // total time
    double tempPS; // temporary plastic strain
    int yield;
    double acctime;
    double dt;
    int unstable_time;
    int disU[12]; // wether to update dislocation
    int total_step;
    int round;
    double totalWallTime, totalWalltimeCP, totalWalltimeSCD;
    double dSCDTime;
    double dCPTime;
    int searchStepCP;
    int executeStepSCD;
    double Th; // Theta
    double Ph; // Phi
    //Matrix3d F, F_L, F_P, F_theta, F_E;
    /**
     * F : total deformation gradient
     * F_L : lattice deformation gradient
     * F_P : lattice preserving plastic deformation gradient
     * F_theta : thermal expansion component
     * F_E : elastic deformation gradient
     * F = F_L * F_P = F_E * F_theta * F_P
     **/
    //Matrix3d L_P; //plastic velocity (correspond to plastic strain-rate) gradient --> plastic strain rate
    //Matrix3d dF_P, dF_theta;
    /**
     * dF_P : plastic deformation rate
     * dF_theta : thermal expansion rate
     **/
    double dGamma[12], V[12], B_drag[12]; //, deltaG_kp[12];
    double one_over_L; // 1/L [cm^-1] reverse of point-like defects with effective average spacing L
    /**
     * dGamma : slip rate
     * V : [cm/s] dislocation velocity
     * B_drag : [Pa*s] drag coefficient
     * deltaG_kp : kink-pair nucleation free energy
     **/
    double Rho[12], dRho[12], Rho_f[12], tau[12]; // Lambda[12],  , L, Xi[12];
    double tempRho[12];
    double Rho_int[12];
    /**
     * Rho : [cm^-2] dislocation density
     * dRho : dislocation density rate (small increment of dislocation density)
     * Lambda : dislocation mean free path
     * Rho_f : forest density
     * tau: [GPa] resolved shear stress 1 GPa = 1e+9 Pa
     * L: effective average spacing of point-like defects
     * Xi: ratio of resolved shear stress and the Peierls stress
     **/
    /* declaration of necessary matrices */
    Matrix3d strain_tot; // total strain matrix
    Matrix3d strain_elastic; // elastic strain matrix
    Matrix3d strain_plastic; // plastic strain matrix
    Matrix3d Cauchy_stress; // cauchy stress tensor
    Vector6d stress_vec; //stress vector
    Vector6d strain_vec; //strain vector
    Matrix6d stiffness; // stiffness matrix for Hookie's law
    //Matrix3d dStrainT; // total strain increament
    //Matrix3d dStrainE; // elastic strain increament
    fstream fss;
    fstream fsr;
    fstream fdd;
    fstream ft;
    int sign(double&); //sign function
    void computeStiffness(); // compute stiffness matrix
    Matrix3d computeDStrainT(); // compute total strain
    Matrix3d computeDPlasticStrain(double&);
    void computeStrainE(); // compute elastic strain
    bool computeOneRSS(int&); // compute one resolved applied stress
    void computeOneRho_f(int&); // compute one forest dislocation density
    //void computeOneLambda(int&); // compute one mean free path
    void computeOneDRho(int&); // compute one dislocation density rate
    //void computeOneXi(int&); // compute one Xi
    void computeOneB_drag(int&); // compute one drag coefficient
    void computeOneV(int&); // compute one dislocation velocity
    void computeOneDGamma(int&); // compute one slip rate
    //void computeLP(); // compute L_P
    //void computeL(); // compute average spacing L
    Matrix3d computeDCauchyStress(Matrix3d&);
    void updateTempDislocationDensity(int&);
    double putTempToRho(); // put temporary dislocation density to official dislocation density
    void plotStressStrainTensors();
    
public:
    DeformationGradient(double&, double&);
    /* get functions */
    Matrix6d getStiffness(); // return stiffness matrix
    //Matrix3d getDStrainT();
    int deformation(); // process whole deformation
    double getMainStrain();
    void initializeParameters(double&, double&);
    ~DeformationGradient(){};
};

