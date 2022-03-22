//
//  main.cpp
//  CP
//
//  Created by fx on 7/7/20.
//  Copyright © 2020 fx. All rights reserved.
//  This code is developed based on Algorithm 1 created by Jaime.
//  Equations of some parameters is not accurate and used the simplified numbers. (i.e. drag coefficient)
//

#include <iostream>
#include "deformationGradient.hpp"

int main(int argc, const char * argv[]) {
    double theta = 0.0;
    double phi = 0.0;
    int unstable;
    int file = 0;
    fstream ftp;
    DeformationGradient F(theta, phi);
    double strainT = F.getMainStrain();
    double oneDegree = 1.0/180.0*PI;
    if(Scenario == 64){
        while(theta <= PI/4.0){
            phi = 0.0;
            while(phi <= PI/4.0){
                /*
                file++;
                ftp.open("theta_phi.txt", ios::app);
                ftp << file++ <<" " << cos(theta)*sin(phi)*(1.0/(cos(phi)+1.0)) << "   " << sin(theta)*sin(phi)*(1.0/(cos(phi)+1.0)) <<endl;
                ftp.close();
                */
                F.initializeParameters(theta, phi);
                strainT = F.getMainStrain();
                while(strainT <= TotalStrain){
                    unstable = F.deformation();
                    strainT = F.getMainStrain();
                    if(unstable == -1){
                        break;
                    }
                }
                phi += oneDegree;
            }
            theta += oneDegree;
        }
        return 0;
    }else{
        while(strainT <= TotalStrain){
            F.deformation();
            strainT = F.getMainStrain();
            if(unstable == -1){
                break;
            }
        }
    }
    return 0;
}

