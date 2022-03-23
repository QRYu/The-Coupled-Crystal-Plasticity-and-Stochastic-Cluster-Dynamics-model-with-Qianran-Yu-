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
    double tot_dpa = 0.2;
    int unstable;
    //int file = 0;
    //fstream ftp;
    DeformationGradient F(theta, phi);
    double strainT = F.getMainStrain();
    if(Scenario == 64){
        while(theta <= PI/4.0){
            phi = 0.0;
            while(phi <= 54.7*PI/180.0){
                //file++;
                //ftp.open("theta_phi.txt", ios::app);
                //ftp <<file << " " << theta << "   " << phi <<endl;
                //ftp.close();
                F.initializeParameters(theta, phi);
                strainT = F.getMainStrain();
                while(F.getTime() <= TOTAL_TIME){
                    unstable = F.deformation();
                    strainT = F.getMainStrain();
                    if(unstable == -1){
                        break;
                    }
                }
                phi += 0.02*PI;
            }
            theta += 0.025*PI;
        }
    }else{
        //while(F.getDPA() <= tot_dpa){
        while(F.getTime()<=TOTAL_TIME){
            F.deformation();
            //strainT = F.getMainStrain();
            if(unstable == -1){
                break;
            }
        }
    }
    
    return 0;
}

