//
//  slipSystem.cpp
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//

#include "slipSystem.hpp"
/* public functions */
SlipSystem::SlipSystem(double& Theta, double& Phi){
    
    /* slip direction -60 degree of slip plane */
    /*
     slipDirection_60[0]<< -1.0, 0.0, 1.0;
     slipDirection_60[1]<< -1.0, -1.0, 0.0;
     slipDirection_60[2]<< 1.0, -1.0, 0.0;
     slipDirection_60[3]<< 1.0, 0.0, 1.0;
     slipDirection_60[4]<< 1.0, 1.0, 0.0;
     slipDirection_60[5]<< 0.0, 1.0, 1.0;
     slipDirection_60[6]<< 0.0, -1.0, 1.0;
     slipDirection_60[7]<< -1.0, -1.0, 0.0;
     slipDirection_60[8]<< 0.0, 1.0, -1.0;
     slipDirection_60[9]<< 1.0, 0.0, -1.0;
     slipDirection_60[10]<< -1.0, 0.0, 1.0;
     slipDirection_60[11]<< 0.0, 1.0, 1.0;
     */
    /* compute various relationships */
    updateSlipSystems(Theta, Phi);
    /*
     for(int i = 0; i<12; i++){
     cout << i+1 <<" = "<<slipDirection[i].transpose()*Cauchy_stress*slipPlaneNormal[i]<<endl;
     }
     */
}

Vector3d SlipSystem::getSlipDirection(int& count){
    
    return slipDirection[count];
    
}

Vector3d SlipSystem::getSlipPlaneNormal(int& count){
    
    return slipPlaneNormal[count];
    
}

/*
 Vector3d SlipSystem::getSlipDirection_60(int& count){
 
 return slipDirection_60[count];
 
 }
 */

Matrix3d SlipSystem::getSCCN(int& count){
    
    return sCCn[count];
    
}
/*
 Matrix3d SlipSystem::getSCCN1(int& count){
 
 return sCCn1[count];
 
 }
 */
Matrix3d SlipSystem::getNCSCCN(int& count){
    
    return nCsCCn[count];
    
}
/*
 Matrix3d SlipSystem::getN1CSCCN1(int& count){
 
 return n1CsCCn1[count];
 
 }
 */
/* private fuctions */
void SlipSystem::computeSlipSystems(){
    for(int i = 0; i < 12; i++){
        sCCn[i] = slipDirection[i]*slipPlaneNormal[i].transpose();
        //sCCn1[i] = slipDirection[i]*slipDirection_60[i].transpose();
        nCsCCn[i] = slipPlaneNormal[i].cross(slipDirection[i])*slipPlaneNormal[i].transpose();
        //n1CsCCn1[i] = slipDirection_60[i].cross(slipDirection[i])*slipDirection_60[i].transpose();
    }
}

void SlipSystem::updateSlipSystems(double& Theta, double& Phi){
    double magnitude_Pla = 1.0/sqrt(2.0);
    double magnitude_Dir = 1.0/sqrt(3.0);
    slipDirection[0]<< 1.0*magnitude_Dir, -1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[1]<< -1.0*magnitude_Dir, -1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[2]<< 1.0*magnitude_Dir, 1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[3]<< -1.0*magnitude_Dir, 1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[4]<< -1.0*magnitude_Dir, 1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[5]<< -1.0*magnitude_Dir, -1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[6]<< 1.0*magnitude_Dir, 1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[7]<< 1.0*magnitude_Dir, -1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[8]<< -1.0*magnitude_Dir, 1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[9]<< -1.0*magnitude_Dir, 1.0*magnitude_Dir, -1.0*magnitude_Dir;
    slipDirection[10]<< 1.0*magnitude_Dir, 1.0*magnitude_Dir, 1.0*magnitude_Dir;
    slipDirection[11]<< 1.0*magnitude_Dir, 1.0*magnitude_Dir, -1.0*magnitude_Dir;
    /* initialize slip plane normal */
    slipPlaneNormal[0]<< 0.0, 1.0*magnitude_Pla, 1.0*magnitude_Pla;
    slipPlaneNormal[1]<< 0.0, 1.0*magnitude_Pla, 1.0*magnitude_Pla;
    slipPlaneNormal[2]<< 0.0, -1.0*magnitude_Pla, 1.0*magnitude_Pla;
    slipPlaneNormal[3]<< 0.0, -1.0*magnitude_Pla, 1.0*magnitude_Pla;
    slipPlaneNormal[4]<< 1.0*magnitude_Pla, 0.0, 1.0*magnitude_Pla;
    slipPlaneNormal[5]<< 1.0*magnitude_Pla, 0.0, 1.0*magnitude_Pla;
    slipPlaneNormal[6]<< -1.0*magnitude_Pla, 0.0, 1.0*magnitude_Pla;
    slipPlaneNormal[7]<< -1.0*magnitude_Pla, 0.0, 1.0*magnitude_Pla;
    slipPlaneNormal[8]<< 1.0*magnitude_Pla, 1.0*magnitude_Pla, 0.0;
    slipPlaneNormal[9]<< 1.0*magnitude_Pla, 1.0*magnitude_Pla, 0.0;
    slipPlaneNormal[10]<< -1.0*magnitude_Pla, 1.0*magnitude_Pla, 0.0;
    slipPlaneNormal[11]<< -1.0*magnitude_Pla, 1.0*magnitude_Pla, 0.0;
    theta = Theta;
    phi = Phi;
    computeTransfromMatrix();
    transferSlipSystems();
    computeSlipSystems();
    Matrix3d Cauchy_stress;
    Cauchy_stress<< 0,0,0,
    0,1,0,
    0,0,0;
    
}

void SlipSystem::computeTransfromMatrix(){
    // coordinate transformation
    //original reference: x = [100], y = [010], z = [001]
    // declare original reference
    Vector3d x(1.0,0.0,0.0);
    Vector3d y(0.0,1.0,0.0);
    Vector3d z(0.0,0.0,1.0);
    Vector3d x1;
    Vector3d y1;
    Vector3d z1;
    if(Scenario == 1){
        /** Scenario 1
         * current referenct: x1 = [1-10], y1 = [11-2], z = [111]
         **/
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2), 0.0;
        y1 << 1.0/sqrt(6.0), 1.0/sqrt(6.0),-2.0*1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0),1.0/sqrt(3.0);
    }else if(Scenario == 2){
        /** Scenario 2
         * current referenct: x1 = [11-2], y1 = [111], z = [1-10]
         **/
        x1 << 1.0/sqrt(6.0), 1.0/sqrt(6.0),-2.0*1.0/sqrt(6.0);
        y1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0),1.0/sqrt(3.0);
        z1 << 1.0/sqrt(2.0), -1.0/sqrt(2), 0.0;
    }else if(Scenario == 3){
        /** Scenario 3
         * current referenct: x1 = [11-2], y1 = [-110], z = [111]
         **/
        x1 << 1.0/sqrt(6.0), 1.0/sqrt(6.0),-2.0*1.0/sqrt(6.0);
        y1 << -1.0/sqrt(2.0), 1.0/sqrt(2), 0.0;
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0),1.0/sqrt(3.0);
    }else if(Scenario == 4){
        x1 << 0.0, 0.0, 1.0;
        y1 << 1.0, 0.0, 0.0;
        z1 << 0.0, 1.0, 0.0;
    }else if(Scenario == 5){
        x1 << 1.0/sqrt(21.0)*(-4.0), 1.0/sqrt(21.0)*(-1.0),1.0/sqrt(21.0)*2.0;
        y1 << 1.0/sqrt(14.0)*1.0, 1.0/sqrt(14.0)*2.0, 1.0/sqrt(14.0)*3.0;
        z1 << 1.0/sqrt(6.0)*(-1.0),1.0/sqrt(6.0)*2.0, 1.0/sqrt(6.0)*(-1.0);
    }else if(Scenario == 6){
        x1 << 1.0/sqrt(38.0)*(-1.0),1.0/sqrt(38.0)*(-1.0),1.0/sqrt(38.0)*(6.0);
        y1 << 1.0/sqrt(19.0)*3.0, 1.0/sqrt(19.0)*3.0, 1.0/sqrt(19.0)*1.0;
        z1 << 1.0/sqrt(2.0)*1.0,1.0/sqrt(2.0)*(-1.0), 1.0/sqrt(2.0)*(0.0);
    }else if(Scenario == 7){
        x1 << 1.0/sqrt(147.0)*(-11.0),1.0/sqrt(147.0)*(-5.0),1.0/sqrt(147.0)*(1.0);
        y1 << 1.0/sqrt(98.0)*(-1.0), 1.0/sqrt(98.0)*4.0, 1.0/sqrt(98.0)*9.0;
        z1 << 1.0/sqrt(6.0)*(-1.0),1.0/sqrt(6.0)*(2.0), 1.0/sqrt(6.0)*(-1.0);
    }else if(Scenario == 8){
        y1 << -1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
        x1 << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        z1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 9){
        y1 << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        x1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
        z1 << -1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 10){
        y1 << 0.0, 0.0, 1.0;
        x1 << 1.0, 0.0, 0.0;
        z1 << 0.0, -1.0, 0.0;
    }else if(Scenario == 11){
        y1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
        x1 << 1.0/sqrt(2), 1.0/sqrt(2.0), 0.0;
        z1 << 1.0/sqrt(6.0), -1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 12){
        x1 << -1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
        z1 << -1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 13){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0);
        z1 << -1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 14){
        x1 << -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
        z1 << 1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 15){
        x1 << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0);
        z1 << -1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 16){
        x1 << -1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0);
        z1 << 1.0/sqrt(6.0), -1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 17){
        x1 << -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0);
        z1 << -1.0/sqrt(6.0), -1.0/sqrt(6.0), 2.0/sqrt(6.0);
    }else if(Scenario == 18){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, 1.0, 0.0;
        z1 << 0.0, 0.0, 1.0;
    }else if(Scenario == 19){
        x1 << 0.0, 1.0, 0.0;
        y1 << -1.0, 0.0, 0.0;
        z1 << 0.0, 0.0, 1.0;
    }else if(Scenario == 20){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, 0.0, -1.0;
        z1 << 0.0, 1.0, 0.0;
    }else if(Scenario == 21){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, -1.0, 0.0;
        z1 << 0.0, 0.0, -1.0;
    }else if(Scenario == 22){
        x1 << -1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        z1 << 0.0, 0.0, 1.0;
    }else if(Scenario == 23){
        x1 << -1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        z1 << 0.0, 1.0, 0.0;
    }else if(Scenario == 24){
        x1 << 0.0, 1.0, 0.0;
        y1 << -1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        z1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
    }else if(Scenario == 25){
        x1 << 0.0, -1.0, 0.0;
        y1 << 1.0/sqrt(2.0), 0.0, -1.0/sqrt(2.0);
        z1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
    }else if(Scenario == 26){
        x1 << 0.0, 1.0, 0.0;
        y1 << -1.0/sqrt(2.0), 0.0, -1.0/sqrt(2.0);
        z1 << -1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
    }else if(Scenario == 27){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0);
        z1 << 0.0, -1.0/sqrt(2.0),1.0/sqrt(2.0);
    }else if(Scenario == 28){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, -1.0/sqrt(2.0), 1.0/sqrt(2.0);
        z1 << 0.0, -1.0/sqrt(2.0), -1.0/sqrt(2.0);
    }else if(Scenario == 29){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0);
        z1 << 0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0);
    }else if(Scenario == 30){
        x1 << 1.0, 0.0, 0.0;
        y1 << 0.0, -1.0/sqrt(2.0), -1.0/sqrt(2.0);
        z1 << 0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0);
    }else if(Scenario == 31){
        x1 << -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        z1 << 0.0, 0.0, 1.0;
    }else if(Scenario == 32){
        x1 << -1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(6.0), -1.0/sqrt(6.0), 2.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 33){
        x1 << 0.0, -1.0/sqrt(2.0), 1.0/sqrt(2.0);
        y1 << 2.0/sqrt(6.0), 1.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 34){
        x1 << 0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0);
        y1 << -2.0/sqrt(6.0), 1.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 35){
        x1 << 0.0, -1.0/sqrt(2.0), -1.0/sqrt(2.0);
        y1 << 2.0/sqrt(6.0), -1.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 36){
        x1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
        y1 << 2.0/sqrt(6.0), 1.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << 0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0);
    }else if(Scenario == 37){
        x1 << 1.0/sqrt(2.0), 0.0, -1.0/sqrt(2.0);
        y1 << 1.0/sqrt(6.0), 2.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 38){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << -1.0/sqrt(6.0), 2.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 39){
        x1 << -1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << 1.0/sqrt(6.0), -2.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 40){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << 1.0/sqrt(6.0), 2.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 41){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(6.0), -1.0/sqrt(6.0), 2.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 42){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(6.0), -1.0/sqrt(6.0), -2.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 43){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(6.0), 1.0/sqrt(6.0), 2.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 44){
        x1 << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(6.0), -1.0/sqrt(6.0), -2.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 45){
        x1 << 0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0);
        y1 << -2.0/sqrt(6.0), -1.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 46){
        x1 << 0.0, 1.0/sqrt(2.0), 1.0/sqrt(2.0);
        y1 << -2.0/sqrt(6.0), 1.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 47){
        x1 << 0.0, -1.0/sqrt(2.0), 1.0/sqrt(2.0);
        y1 << 2.0/sqrt(6.0), -1.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 48){
        x1 << 0.0, 1.0/sqrt(2.0), -1.0/sqrt(2.0);
        y1 << -2.0/sqrt(6.0), -1.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 49){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << -1.0/sqrt(6.0), -2.0/sqrt(6.0), 1.0/sqrt(6.0);
        z1 << 1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 50){
        x1 << -1.0/sqrt(2.0), 0.0, -1.0/sqrt(2.0);
        y1 << 1.0/sqrt(6.0), -2.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0/sqrt(3.0);
    }else if(Scenario == 51){
        x1 << -1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << -1.0/sqrt(6.0), 2.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 52){
        x1 << 1.0/sqrt(2.0), 0.0, -1.0/sqrt(2.0);
        y1 << -1.0/sqrt(6.0), -2.0/sqrt(6.0), -1.0/sqrt(6.0);
        z1 << -1.0/sqrt(3.0), 1.0/sqrt(3.0), -1.0/sqrt(3.0);
    }else if(Scenario == 53){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(51.0), 1.0/sqrt(51.0), 7.0/sqrt(51.0);
        z1 << -7.0/sqrt(102.0), -7.0/sqrt(102.0), 2.0/sqrt(102.0);
    }else if(Scenario == 54){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(27.0), 1.0/sqrt(27.0), 5.0/sqrt(27.0);
        z1 << -5.0/sqrt(54.0), -5.0/sqrt(54.0), 2.0/sqrt(54.0);
    }else if(Scenario == 55){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 1.0/sqrt(11.0), 1.0/sqrt(11.0), 3.0/sqrt(11.0);
        z1 << -3.0/sqrt(22.0), -3.0/sqrt(22.0), 2.0/sqrt(22.0);
    }else if(Scenario == 56){
        x1 << 1.0/sqrt(2.0), -1.0/sqrt(2.0), 0.0;
        y1 << 3.0/sqrt(43.0), 3.0/sqrt(43.0), 5.0/sqrt(43.0);
        z1 << -5.0/sqrt(86.0), -5.0/sqrt(86.0), 6.0/sqrt(86.0);
    }else if(Scenario == 57){
        x1 << 3.0/sqrt(34.0), -5.0/sqrt(34.0), 0.0;
        y1 << 5.0/sqrt(59.0), 3.0/sqrt(59.0), 5.0/sqrt(59.0);
        z1 << -5.0/sqrt(80.24), -3.0/sqrt(80.24), 6.8/sqrt(80.24);
    }else if(Scenario == 58){
        x1 << 1.0/sqrt(5.0), -2.0/sqrt(5.0), 0.0;
        y1 << 2.0/sqrt(9.0), 1.0/sqrt(9.0), 2.0/sqrt(59.0);
        z1 << -2.0/sqrt(11.25), -1.0/sqrt(11.25), 2.5/sqrt(11.25);
    }else if(Scenario == 59){
        x1 << 1.0/sqrt(10.0), -3.0/sqrt(10.0), 0.0;
        y1 << 3.0/sqrt(19.0), 1.0/sqrt(19.0), 3.0/sqrt(19.0);
        z1 << -3.0/sqrt(21.11), -1.0/sqrt(21.11), 3.3333/sqrt(21.11);
    }else if(Scenario == 60){
        x1 << 1.0/sqrt(10.0), -3.0/sqrt(10.0), 0.0;
        y1 << 3.0/sqrt(35.0), 1.0/sqrt(35.0), 5.0/sqrt(35.0);
        z1 << -3.0/sqrt(14.0), -1.0/sqrt(14.0), 2.0/sqrt(14.0);
    }else if(Scenario == 61){
        x1 << 1.0/sqrt(5.0), -2.0/sqrt(5.0), 0.0;
        y1 << 2.0/sqrt(14.0), 1.0/sqrt(14.0), 3.0/sqrt(14.0);
        z1 << -2.0/sqrt(7.78), -1.0/sqrt(7.78), 1.66667/sqrt(7.78);
    }else if(Scenario == 62){
        x1 << 0.0, -1.0, 0.0;
        y1 << 1.0/sqrt(5.0), 0.0, 2.0/sqrt(5.0);
        z1 << -2.0/sqrt(5.0), 0.0, 1.0/sqrt(5.0);
    }else if(Scenario == 63){
        x1 << 0.0, -1.0, 0.0;
        y1 << 1.0/sqrt(10.0), 0.0, 3.0/sqrt(10.0);
        z1 << -3.0/sqrt(10.0), 0.0, 1.0/sqrt(10.0);
    }else if(Scenario == 64){
        x1 << sin(theta), -cos(theta), 0.0;
        y1 << sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi);
        z1 << cos(theta), sin(theta), -tan(phi);
    }else if(Scenario == 65){
        x1 << 3.0/sqrt(10.0), 0.0, 1.0/sqrt(10.0);
        y1 << -1.0/sqrt(10.0), 0.0, 3.0/sqrt(10.0);
        z1 << 0.0, -1.0, 0.0;
    }else if(Scenario == 66){
        x1 << 2.0/sqrt(5.0), 0.0, 1.0/sqrt(5.0);
        y1 << -1.0/sqrt(5.0), 0.0, 2.0/sqrt(5.0);
        z1 << 0.0, -1.0, 0.0;
    }else if(Scenario == 67){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << -3.0/sqrt(19.0), 1.0/sqrt(19.0), 3.0/sqrt(19.0);
        z1 << -1.0/sqrt(38.0), -6.0/sqrt(38.0), 1.0/sqrt(38.0);
    }else if(Scenario == 68){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << -2.0/3.0, 1.0/3.0, 2.0/3.0;
        z1 << -1.0/sqrt(18.0), -4.0/sqrt(18.0), 1.0/sqrt(18.0);
    }else if(Scenario == 69){
        x1 << 1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0;
        y1 << -1.0/sqrt(51.0), 1.0/sqrt(51.0), 7.0/sqrt(51.0);
        z1 << 7.0/sqrt(102.0), -7.0/sqrt(102.0), -2.0/sqrt(102.0);
    }else if(Scenario == 70){
        x1 << 1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0);
        y1 << -5.0/sqrt(59.0), 3.0/sqrt(59.0), 5.0/sqrt(59.0);
        z1 << -10.0/sqrt(5206.0/9.0), -65.0/3.0/sqrt(5206.0/9.0), 3.0/sqrt(5206.0/9.0);
    }else if(Scenario == 71){
        x1 << 3.0/sqrt(13.0), 0.0, 2.0/sqrt(13.0);
        y1 << -2.0/sqrt(14.0), 1.0/sqrt(14.0), 3.0/sqrt(14.0);
        z1 << -2.0/sqrt(182.0), -13.0/sqrt(182.0), 3.0/sqrt(182.0);
    }else if(Scenario == 72){
        x1 << 5.0/sqrt(34.0), 0.0, 3.0/sqrt(34.0);
        y1 << -3.0/sqrt(35.0), 1.0/sqrt(35.0), 5.0/sqrt(35.0);
        z1 << -3.0/sqrt(1190.0), -34.0/sqrt(1190.0), 5.0/sqrt(1190.0);
    }else if(Scenario == 73){
        x1 << 5.0/sqrt(26.0), 0.0, 1.0/sqrt(26.0);
        y1 << -1.0/sqrt(27.0), 1.0/sqrt(27.0), 5.0/sqrt(27.0);
        z1 << -1.0/sqrt(702.0), -26.0/sqrt(702.0), 5.0/sqrt(702.0);
    }else if(Scenario == 74){
        x1 << 3.0/sqrt(10.0), 0.0, 1.0/sqrt(10.0);
        y1 << -1.0/sqrt(11.0), 1.0/sqrt(11.0), 3.0/sqrt(11.0);
        z1 << -3.0/sqrt(46.0), -6.0/sqrt(46.0), 1.0/sqrt(46.0);
    }else if(Scenario == 75){
        x1 << 5.0/sqrt(34.0), 0.0, 3.0/sqrt(34.0);
        y1 << -3.0/sqrt(43.0), 3.0/sqrt(43.0), 5.0/sqrt(43.0);
        z1 << -3.0/sqrt(132.22), -34.0/3.0/sqrt(132.22), 1.0/sqrt(132.22);
    }
    
    L<< x1.dot(x), x1.dot(y), x1.dot(z),
    y1.dot(x), y1.dot(y), y1.dot(z),
    z1.dot(x), z1.dot(y), z1.dot(z);
    
}

void SlipSystem::transferSlipSystems(){
    for(int i = 0; i<12; i++){
        slipDirection[i] = L*slipDirection[i];
        slipPlaneNormal[i] = L*slipPlaneNormal[i];
        //slipDirection_60[i] = L*slipDirection_60[i];
    }
}

