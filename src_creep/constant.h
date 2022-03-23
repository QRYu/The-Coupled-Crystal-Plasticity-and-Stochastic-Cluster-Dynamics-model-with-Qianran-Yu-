//
//  constant.h
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//
#pragma once
#define KB 8.617e-05       // [ev/K] Boltzmann's constant.
#define ALATT 2.90e-8     // [cm] Lattice parameter of Fe.
#define BURGER 2.51147e-8 // [cm] burgers vector
#define mu 78.0 // [GPa] shear modulus (Fe 78 - 84)
#define Self 0.009
#define Coplanar 0.009
#define Collinear 0.72
#define Orthorgonal 0.05
#define Glissile 0.09
#define Sessile 0.06 //dislocation interaction types
#define E 204.0 // [GPa] Young's Modulus (Fe 204 - 212 GPa)
#define poisson 0.29 // Poisson ratio (Fe 0.29 - 0.30)
#define tau_p 0.3 // [GPa] Peierls barrier
#define Delta_H0 0.57 // [eV] change of enthalpy at 0 K ??
#define p0 0.67
#define q0 1.18 // fitting parameters for calculating kink-pair nucleation free energy
#define nu0 1e11 // [s^-1]
#define TEMPERATURE 300 //[K] temperature 763 K - 878 K
#define Tm 1805.0 // [K] melting temperature
#define STRESS 0.0 // [MPa] constant stress
//#define DT 0.0000001 // time interval
//#define TotalStrain 0.03

//#define Scenario 1 // loading from [11-2] works
//#define Scenario 8 // loading from [-112] works
//#define Scenario 32 // loading from [1-12] works
//#define Scenario 33 // loading from [211] doesn't work
//#define Scenario 34 // loading from [-211] doesn't work
//#define Scenario 35 // loading from [2-11] doesn't work
//#define Scenario 36 // loading from [21-1] doesn't work
//#define Scenario 37 // loading from [121] doesn't work
//#define Scenario 38 // loading from [-121] doesn't work
//#define Scenario 39 // loading from [1-21] doesn't work
//#define Scenario 40 // loading from [12-1] works
//#define Scenario 41 // loading from [-1-12] works
//#define Scenario 42 // loading from [-1-1-2] works
//#define Scenario 43 // loading from [112] works
//#define Scenario 44 // loading from [1-1-2] works
//#define Scenario 45 // loading from [-2-11] doesn't work
//#define Scenario 46 // loading from [-21-1] doesn't work
//#define Scenario 47 // loading from [2-1-1] doesn't work
//#define Scenario 48 // loading from [-2-1-1] doesn't work
//#define Scenario 49 // loading from [-1-21] work
//#define Scenario 50 // loading from [1-2-1] doesn't work
//#define Scenario 51 // loading from [-12-1] works
//#define Scenario 52 // loading from [-1-2-1] doesn't work

//#define Scenario 2 // loading from [111] works
//#define Scenario 11 // loading from [-111] works
//#define Scenario 12 // loading from [1-11] works
//#define Scenario 13 // loading from [11-1] doesn't work
//#define Scenario 14 // loading from [-1-11] works
//#define Scenario 15 // loading from [-11-1] works
//#define Scenario 16 // loading from [1-1-1] works
//#define Scenario 17 // loading from [-1-1-1] works

//#define Scenario 3 // loading from [-110] doesn't (don't yield)
//#define Scenario 9 // loading from [110] doesn't (don't yield)
//#define Scenario 22 // loading from [1-10] doesn't work
//#define Scenario 23 // loading from [101] works
//#define Scenario 24 // loading from [-101] works
//#define Scenario 25 // loading from [10-1] works
//#define Scenario 26 // loading from [-10-1] works
//#define Scenario 27 // loading from [011] works
//#define Scenario 28 // loading from [0-11] works
//#define Scenario 29 // loading from [01-1] works
//#define Scenario 30 // loading from [0-1-1] works
//#define Scenario 31 // loading from [-1-10] doesn't work


//#define Scenario 4 // loading from [100] doesn't (don't yield)
//#define Scenario 10 //loading from [001] works
//#define Scenario 18  // loading from [010] doesn't work
//#define Scenario 19  // loading from [-100] doesn't work
//#define Scenario 20  // loading from [00-1] works
//#define Scenario 21  // loading from [0-10] doesn't work

//#define Scenario 5 // loading from [123] works
//#define Scenario 6 // loading from [331] works
//#define Scenario 7 // loading from [-149] works
/* below are factors all related to non-schmid effect

/* directions that will be tested */
//#define Scenario 43 // loading from [112] works
//#define Scenario 23 // loading from [101] works
#define Scenario 10 //loading from [001] works
//#define Scenario 2 // loading from [111] works
//#define Scenario 53 // loading from [117] works
//#define Scenario 54 // loading from [115] works
//#define Scenario 55 // loading from [113] works
//#define Scenario 56 // loading from [335] works
//#define Scenario 57 // loading from [535] unstable
//#define Scenario 58 // loading from [212] unstable
//#define Scenario 59 // loading from [313] unstable
//#define Scenario 60 // loading from [315] works
//#define Scenario 61 // loading from [213] works
//#define Scenario 62 // loading from [102] works
//#define Scenario 63 // loading from [103] works
//#define Scenario 64 // loading from (thita, phi)

 
 
 
 //#define B_0 9.8e-4
 //#define B_1 0
 //#define B_k 8.3e-5 /* [Pa*s], [Pa*s/K], [Pa*s] for calculating drag coefficient */



//#define a_0 1.50
//#define a_1 1.15
//#define a_2 2.32
//#define a_3 4.29 // fitting parameters for non-schmid effect criterion

//#define d_g 2.72e-8 // [cm] grain size
//#define d_a 2.72e-8 // [cm] critical dipole annihlation spacing

/* Parameters from SCD */
//#define ION
#define NEUTRON
#define PI 3.1415926
#define DENSITY 8.20e+22   // [atoms/cm^3] Atomic density Fe.
#define ATOMICVOLUME 1.21945e-23 //[cm^3] Atomic volume of Fe.
#define ODS_R 2.5e-07       // [cm] ODS-particle radius.
#define ODS_DENSITY 0.0
#define GRAIN_SIZE 0.0002 //[cm]
#define FOIL_THICKNESS 0.0002 //[cm] Foil thickness from UCSD (2000nm)
#define SURFACE_THICKNESS 0.544 //[nm] thickness of surface (conrresponds to two monolayers of tungsten)
#define NU0 6.1e+12           // [Hz] Attempt frequency.
#define C_DENSITY 10        // [appm] C-atom density
#define GAMMA 1.0           // Fraction of surface emission.
#define TDE 25              // [eV] Threshold displacement energy for W.
#define TOTAL_TIME 51724137.9  // [s] Total simulated time.
#define VOLUME 1.0e-14 // [cm^3] System volume.
#define RATIO_HE 0       // [appm/dpa] He-to-dpa ratio.
#define RATIO_H 0
#define CHANNELS 3         // Irradiation channels used (1:W, 2:He, 3:H,...). the number of different particle insertion(irradiation) process.
#define PSTEPS 50000 // Print data every so many.
#define TSTEPS 50000 // Run these many steps.
#define LEVELS 3
#define EXP10 3
#define POINTS 1       // number of elements: one surface(Point 0), other bulk elements(NO.1,2,3,4...)
#define dpa_rate 2.9e-7 //dpa rate [dpa/s]
// Auxiliary definitions:
enum Reaction { DIFFUSETOF, DIFFUSETOB, SINK, DISSOCIATION, COMBINATION, NONE, PARTICLE, HE, H, DISSV, DISSH, ERROR};
const double AVG_ION_EN[POINTS] = {131.304};
#define AVG_NEUTRON_EN 18.8 // (from SPECTER) Total damage energy in keV produced by a neutron in ITER.
//#define AVG_NEUTRON_EN 100
#define RESTART            // Do restart. */
