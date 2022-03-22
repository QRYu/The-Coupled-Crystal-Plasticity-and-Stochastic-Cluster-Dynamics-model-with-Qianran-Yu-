//
//  constant.h
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//
#pragma once
#define PI 3.1415926
#define KB 8.617e-05       // [ev/K] Boltzmann's constant.
#define ALATT 3.165e-8 // [cm] lattice parameter of
#define BURGER 2.80e-8 // [cm] burgers vector
#define mu 130.0 // [GPa] shear modulus (tungsten 130 - 160 GPa) -- 88 is froma paper
#define Self 0.009
#define Coplanar 0.009
#define Collinear 0.72
#define Orthorgonal 0.05
#define Glissile 0.09
#define Sessile 0.06 //dislocation interaction types
#define E 340.0 // [GPa] Young's Modulus (tungsten 340 - 405 GPa) -- calculated from mu = E/2/(1+Poisson)
#define poisson 0.27 // Poisson ratio (tungsten 0.27 - 0.29)
#define tau_p 2.03 // [GPa] Peierls barrier
//#define dt 1e-5 // [s] time increament
#define Delta_H0 1.63 // [eV] change of enthalpy at 0 K ??
#define p0 0.86
#define q0 1.69 // fitting parameters for calculating kink-pair nucleation free energy
#define nu0 1e11 // [s^-1]
#define TEMPERATURE 600 //[K] temperature
#define Tm 3695.0 // [K] melting temperature
#define TotalStrain 0.03
#define DT 1e-5 // [s]

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

//#define Scenario 44 // loading from [1-1-2] works
//#define Scenario 45 // loading from [-2-11] doesn't work
//#define Scenario 46 // loading from [-21-1] doesn't work
//#define Scenario 47 // loading from [2-1-1] doesn't work
//#define Scenario 48 // loading from [-2-1-1] doesn't work
//#define Scenario 49 // loading from [-1-21] work
//#define Scenario 50 // loading from [1-2-1] doesn't work
//#define Scenario 51 // loading from [-12-1] works!!
//#define Scenario 52 // loading from [-1-2-1] doesn't work

//#define Scenario 2 // loading from [111] works
//#define Scenario 12 // loading from [1-11] works
//#define Scenario 13 // loading from [11-1] doesn't work
//#define Scenario 14 // loading from [-1-11] works
//#define Scenario 15 // loading from [-11-1] works!!
//#define Scenario 16 // loading from [1-1-1] works
//#define Scenario 17 // loading from [-1-1-1] works

//#define Scenario 3 // loading from [-110] doesn't (don't yield)
//#define Scenario 9 // loading from [110] doesn't (don't yield)
//#define Scenario 22 // loading from [1-10] doesn't work

//#define Scenario 24 // loading from [-101] doesn't works
//#define Scenario 26 // loading from [-10-1] works
//#define Scenario 27 // loading from [011] works!!
//#define Scenario 28 // loading from [0-11] works
//#define Scenario 30 // loading from [0-1-1] works!!
//#define Scenario 31 // loading from [-1-10] doesn't work


//#define Scenario 4 // loading from [100] doesn't (don't yield)
//#define Scenario 20  // loading from [00-1] works
//#define Scenario 18  // loading from [010] doesn't work
//#define Scenario 19  // loading from [-100] doesn't work
//#define Scenario 21  // loading from [0-10] doesn't work
//#define Scenario 5 // loading from [123] doesn't works
//#define Scenario 6 // loading from [331] works!!


/* directions that will be tested */
//#define Scenario 23 // loading from [101] doesn't works
//#define Scenario 10 //loading from [001] doesn't works
//#define Scenario 2 // loading from [111] doesn't works
//#define Scenario 53 // loading from [117] doesn't works
//#define Scenario 54 // loading from [115] doesn't works
//#define Scenario 62 // loading from [102] doesn't works
//#define Scenario 63 // loading from [103] doesn't works




/* directions that will be tested_B */
//#define Scenario 64 // loading from (thita, phi)
//#define Scenario 65 // loading from [-103] doesn't work
//#define Scenario 66 // loading from [-102] doesn't work
//#define Scenario 67 //loading from [-313] doesn't work
//#define Scenario 68 // loading from [-212] doesn't work
//#define Scenario 69 // loading from [-117] doesn't work
//#define Scenario 70 // loading from [-535] doesn't work
//#define Scenario 71 // loading from [-213] doesn't work
//#define Scenario 72 // loading from [-315] doesn't work
//#define Scenario 74 // loading from [-113] doesn't work
//#define Scenario 8 // loading from [-112] works!!
//#define Scenario 75 // loading from [-335] doesn't work
//#define Scenario 24 // loading from [-101] doesn't works
//#define Scenario 55 // loading from [113] works!!
//#define Scenario 73 // loading from [-115] works!!


//#define Scenario 59 // loading from [313] works!!
//#define Scenario 60 // loading from [315] works!!
#define Scenario 61 // loading from [213] works!!
//#define Scenario 7 // loading from [-149] works
//#define Scenario 25 // loading from [10-1] works


//#define Scenario 43 // loading from [112] work!!

//#define Scenario 57 // loading from [535] works!!

//#define Scenario 56 // loading from [335] works!!
//#define Scenario 11 // loading from [-111] works!!
//#define Scenario 58 // loading from [212] works
//#define Scenario 29 // loading from [01-1] works!!



/* below are factors all related to non-schmid effect
 
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
//#define NEUTRON
#define ION
#define PI 3.1415926
#define DENSITY 6.30705e+22   // [atoms/cm^3] Atomic density for W.
#define ATOMICVOLUME 0.0158 //[nm^3] Atomic volume of W.
#define ODS_R 2.5e-07       // [cm] ODS-particle radius.
#define ODS_DENSITY 0.0
#define GRAIN_SIZE 0.0002 //[cm]
#define FOIL_THICKNESS 0.0002 //[cm] Foil thickness from UCSD (2000nm)
#define SURFACE_THICKNESS 0.544 //[nm] thickness of surface (conrresponds to two monolayers of tungsten)
#define NU0 6.1e+12           // [Hz] Attempt frequency.
#define C_DENSITY 10        // [appm] C-atom density
#define GAMMA 1.0           // Fraction of surface emission.
#define TDE 90              // [eV] Threshold displacement energy for W.
#define TOTAL_TIME 200  // [s] Total simulated time.
#define VOLUME 1.0e-13 // [cm^3] System volume.
#define RATIO_HE 1.1       // [appm/dpa] He-to-dpa ratio.
#define RATIO_H 0
#define CHANNELS 3         // Irradiation channels used (1:W, 2:He, 3:H,...). the number of different particle insertion(irradiation) process.
#define PSTEPS 50000 // Print data every so many.
#define TSTEPS 50000 // Run these many steps.
#define LEVELS 3
#define EXP10 3
#define POINTS 1       // number of elements: one surface(Point 0), other bulk elements(NO.1,2,3,4...)
// Auxiliary definitions:
enum Reaction { DIFFUSETOF, DIFFUSETOB, SINK, DISSOCIATION, COMBINATION, NONE, PARTICLE, HE, H, DISSV, DISSH, ERROR};
const double AVG_ION_EN[POINTS] = {765.591};
#define AVG_NEUTRON_EN 40.6 // (from SPECTER) Total damage energy in keV produced by a neutron in ITER.
#define RESTART            // Do restart. */
