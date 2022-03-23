//
//  deformationGradient.cpp
//  CP_Algorithm2
//
//  Created by fx on 7/15/20.
//  Copyright © 2020 fx. All rights reserved.
//

#include "deformationGradient.hpp"

/* public functions */
DeformationGradient::DeformationGradient(double& theta, double& phi):slip_systems(theta, phi),srscd(){
    //maxiter = 50;
    
    //totalStrainRate = 10.0;
    //TOL = 1e-6*dt*totalStrainRate;
    //tempPS = 0.0;
    //unstable_time = 0;
    /*
     lambda[0] = 0.15;
     lambda[1] = 0.25;
     lambda[2] = 0.30;
     lambda[3] = 0.35;
     lambda[4] = 0.40;
     lambda[5] = 0.60;
     */
    computeStiffness();
    srand(time(0));
    xi[0][0] = xi[1][1] = xi[2][2] = xi[3][3] = xi[4][4] =xi[5][5] = xi[6][6] = xi[7][7] = xi[8][8] = xi[9][9] = xi[11][11] = Self;
    xi[1][0] = xi[3][2] = xi[5][4] = xi[7][6] = xi[9][8] = xi[11][10] = xi[0][1] = xi[2][3] = xi[4][5] = xi[6][7] = xi[8][9] = xi[10][11] = Coplanar;
    xi[10][10] = xi[2][0] = xi[2][1] = xi[3][0] = xi[3][1] = xi[6][4] = xi[6][5] = xi[7][4] = xi[7][5] = xi[10][8] = xi[10][9] = xi[11][8] = xi[11][9] = xi[0][2] = xi[1][2] = xi[0][3] = xi[1][3] = xi[4][6] = xi[5][6] = xi[4][7] = xi[5][7] = xi[8][10] = xi[9][10] = xi[8][11] = xi[9][11] = Sessile;
    xi[4][3] = xi[5][1] = xi[6][2] = xi[7][0] = xi[8][3] = xi[8][4] = xi[9][0] = xi[9][7] = xi[10][2] = xi[10][6] = xi[11][1] = xi[11][5] = xi[3][4] = xi[1][5] = xi[2][6] = xi[0][7] = xi[3][8] = xi[4][8] = xi[0][9] = xi[7][9] = xi[2][10] = xi[6][10] = xi[1][11] = xi[5][11] = Collinear;
    xi[4][0] = xi[5][2] = xi[6][1] = xi[7][3] = xi[8][1] = xi[8][6] = xi[9][2] = xi[9][5] = xi[10][0] = xi[10][4] = xi[11][3] = xi[11][7] = xi[0][4] = xi[2][5] = xi[1][6] = xi[3][7] = xi[1][8] = xi[6][8] = xi[2][9] = xi[5][9] = xi[0][10] = xi[4][10] = xi[3][11] = xi[7][11] = Glissile;
    xi[4][1] = xi[1][4] = xi[4][2] = xi[2][4] = xi[5][0] = xi[0][5] = xi[5][3] = xi[3][5] = xi[6][3] = xi[3][6] = xi[6][0] = xi[0][6] = xi[7][1] = xi[1][7] = xi[7][2] = xi[2][7] = xi[8][0] = xi[0][8] = xi[8][2] = xi[2][8] = xi[8][5] = xi[5][8] = xi[8][7] = xi[7][8] = xi[9][1] = xi[1][9] = xi[9][3] = xi[3][9] = xi[9][4] = xi[4][9] = xi[9][6] = xi[6][9] = xi[10][1] = xi[1][10] = xi[10][3] = xi[3][10] = xi[10][5] = xi[5][10] = xi[10][7] = xi[7][10] = xi[11][0] = xi[0][11] = xi[11][2] = xi[2][11] = xi[11][4] = xi[4][11] = xi[11][6] = xi[6][11] = Orthorgonal;
    /* initialize stress/strain tensors */
    initializeParameters(theta, phi);
    volume_strain = 0.0;
    effective_strain = 0.0;
    volume_strain_old = 0.0;
    effective_strain_old = 0.0;
    DT = 1.0;
    dNabs = 0.0;
    ABSRho = 0.0;
    //L_P = Matrix3d::Zero();
}

Matrix6d DeformationGradient::getStiffness(){
    return stiffness;
}

/*
 Matrix3d DeformationGradient::getDStrainT(){
 return dStrainT;
 }
*/

int DeformationGradient::sign(double& a){
    if(a >= 0.0){
        return 1;
    }else{
        return -1;
    }
}

int DeformationGradient::deformation(){
    
    fstream fy;
    fstream ff;
    fstream fc;
    fstream fv;
    char strainStressFile[30], DDFile[30];
    double incre_ec = 0.0;
    double sumRho1 = 0.0;
    Matrix3d delta_strain_P = Matrix3d::Zero(); //small plastic strain
    Matrix3d small_eg = Matrix3d::Zero();
    Matrix3d small_ec = Matrix3d::Zero();
    if(t == 0.0){
        char strainStressFile[30];
        cout <<  "round = "<< round << endl;
        sprintf(strainStressFile, "strain_stress%d.txt", round);
        fss.open(strainStressFile, ios::app);
        fss << "theta = " << Th << "     phi = "<< Ph << endl;
        acctime = 0.0;
        plotStressStrainTensors();
        // write into files
        sprintf(DDFile, "dd%d.txt", round);
        sprintf(strainStressFile, "strain_stress%d.txt", round);
        fss.open(strainStressFile, ios::app);
        fss << strain_tot(1,1) << "   " << Cauchy_stress(1,1) << endl;
        fss.close();
        fdd.open(DDFile, ios::app);
        fdd << t << "   " << volume_strain<< " "<< effective_strain<<"  "<<sumRho <<"    ";
        for(int i = 0; i<12; i++){
            fdd << Rho_e[i]<<"  "<< Rho_s[i] << "  ";
        }
        fdd<<endl;
        fdd.close();
        fc.open("creep.txt", ios::app);
        //fc << t << "    "<<strain_elastic(1,1) <<"   "<< strain_climb <<"    "<<strain_plastic(1,1)<< "  " << incre_ec/DT <<"   "<<(delta_strain_P-incre_ec)/DT<<"  "<<delta_strain_P/DT<<endl;
        fc << t << "    "<<strain_elastic(1,1) <<"   "<< volume_strain <<"    "<<effective_strain << "  "<< (volume_strain - volume_strain_old)/DT << "  " << (effective_strain-effective_strain_old)/DT <<endl;
        fc.close();
        //fv.open("velocity.txt", ios::app);
        //for(int i = 0; i<12; i++){
            //fv << V_c[i] << "   " << V_e[i] << "    ";
        //}
        //fv << endl;
        //fv.close();
    }
    srscd.doClusterDynamics(one_over_L, Nv_Ni, dNabs, avg_r, DT, ABSRho);
    /*
    while(one_over_L == 0.0){
        t += DT;
        srscd.doClusterDynamics(one_over_L, Nv_Ni, dNabs, avg_r, DT, ABSRho);
    }
    */
    /*
    fstream fN;
    fN.open("bias.txt", ios::app);
    if(Nv_Ni != 0){
        fN<<t<<"    "<<abs(Nv_Ni)/DT<<endl;
    }
    fN.close();
    */
    bool whether;
    for(int i = 0; i < 12; i++){
        //double small_eg = 0.0; // strain caused by glide
        //double small_ec = 0.0; // strain caused by climb
        computeOneRho_f(i);
        //computeOneEffectiveRSS(i); // compute stress
        computeOneVc(i); // compute climb velocity
        computeOneV(i); // compute edge/screw dislocation slip velocity
        computeOneDGamma(i); // compute slip rate
        computeOneDRho(i); // compute dislocation density increase
        updateDislocationDensity(i); // update dislocation density
        small_eg += (abs(dGamma[i])*slip_systems.getSCCN(i)*DT);
        //small_eg += (dGamma[i]*slip_systems.getSCCN(i)*DT);
        computeOneDBeta(i);
        small_ec += (abs(dBeta[i])*slip_systems.getBCCB(i)*DT);
        //small_ec += (dBeta[i]*slip_systems.getBCCB(i)*DT);
        sumRho1 += 2*(Rho_e[i] + Rho_s[i]);
    }
    delta_strain_P += small_eg + small_ec;
    sumRho = sumRho1;
    //sumRho_edge = sumRho_e1;
    Matrix3d hydro = Matrix3d::Zero();
    Matrix3d devi = Matrix3d::Zero();
    t += DT;
    acctime += DT;
    srscd.updateRates(sumRho);
    strain_plastic += delta_strain_P;
    strain_tot += delta_strain_P;
    /* calculate creep strain and swelling */
    hydro << strain_plastic(0,0), 0.0, 0.0,
    0.0, strain_plastic(1,1), 0.0,
    0.0, 0.0, strain_plastic(2,2);
    devi = strain_plastic-1.0/3.0*hydro;  //strain caused by swelling
    //cout << "H^devi = "<< devi << endl;
    //getchar();
    volume_strain = hydro.trace();
    //cout << "plastic strain = " <<strain_plastic<<endl;
    //getchar();
    double e11 = devi(0,0);
    double e12 = devi(0,1);
    double e13 = devi(0,2);
    double e21 = devi(1,0);
    double e22 = devi(1,1);
    double e23 = devi(1,2);
    double e31 = devi(2,0);
    double e32 = devi(2,1);
    double e33 = devi(2,2);
    
    //effective_strain =2.0/3.0*sqrt(1.5*(e11*e11+e22*e22+e33*e33)+0.75*(e12*e12+e23*e23+e13*e13));  // strain caused by creep
    effective_strain =sqrt(2.0/3.0*(e11*e11+e12*e12+e13*e13+e21*e21+e22*e22+e23*e23+e31*e31+e32*e32+e33*e33));
    //double v_eff = effective_strain/sumRho/BURGER/0.5/DT;
    //fstream fv;
    //fv.open("effective.txt", ios::app);
    //if(Nv_Ni != 0){
       // fv<<t<<"    "<<v_eff<<endl;
    //}
    //fv.close();
    if(acctime >= 500*DT){
        acctime = 0.0;
        //cout<<"delta_strain_Pp = "<< delta_strain_Pp << endl;
        plotStressStrainTensors();
        // write into files
        
        sprintf(DDFile, "dd%d.txt", round);
        sprintf(strainStressFile, "strain_stress%d.txt", round);
        fss.open(strainStressFile, ios::app);
        fss << strain_tot(1,1) << "   " << Cauchy_stress(1,1) << endl;
        fss.close();
        fdd.open(DDFile, ios::app);
        fdd << t << "   " << volume_strain << " "<< effective_strain <<"  "<<sumRho <<"    ";
        for(int i = 0; i<12; i++){
            fdd << Rho_e[i]<<"  "<< Rho_s[i] << "  ";
        }
        fdd<<endl;
        fdd.close();
        fc.open("creep.txt", ios::app);
        fc << t << "    "<<strain_elastic(1,1) <<"   "<< volume_strain <<"    "<<effective_strain << "  "<< (volume_strain - volume_strain_old)/DT << "  " << (effective_strain-effective_strain_old)/DT <<endl;
        fc.close();
        effective_strain_old = effective_strain;
        volume_strain_old = volume_strain;
        if(DT < 100){
            DT = DT*1.02;
        }
        cout<<"DT = "<<DT<<endl;
        //fv.open("velocity.txt", ios::app);
        //for(int i = 0; i<12; i++){
           // fv << V_c[i] << "   " << V_e[i] << "    ";
        //}
        //fv << endl;
        //fv.close();
    }
}

double DeformationGradient::getMainStrain(){
    return strain_tot(1,1);
}

/* private functions */
void DeformationGradient::computeStiffness(){
    double prefactor = E/(1+poisson)/(1-2*poisson);
    double element = (1.0-2.0*poisson)/2.0;
    stiffness << prefactor*(1-poisson), prefactor*poisson, prefactor*poisson, 0.0, 0.0, 0.0,
    prefactor*poisson, prefactor*(1-poisson), prefactor*poisson, 0.0, 0.0, 0.0,
    prefactor*poisson, prefactor*poisson, prefactor*(1-poisson), 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, prefactor*element, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, prefactor*element, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, prefactor*element;
}
/*
Matrix3d DeformationGradient::computeDStrainT(){
    Matrix3d d_strain;
    d_strain<< -0.5*totalStrainRate*dt, 0.0, 0.0,
    0.0, totalStrainRate*dt, 0.0,
    0.0, 0.0, -0.5*totalStrainRate*dt;
    return d_strain;
}
*/
Matrix3d DeformationGradient::computeDPlasticStrain(double& principal){
    Matrix3d d_strain;
    d_strain <<-0.5*principal, 0.0, 0.0,
    0.0, principal, 0.0,
    0.0, 0.0, -0.5*principal;
    return d_strain;
}

void DeformationGradient::computeStrainE(){
    
    strain_elastic = strain_tot - strain_plastic;
}

bool DeformationGradient::computeOneEffectiveRSS(int& count){
    
    double dtau_f = mu*BURGER*sqrt(Rho_int[count]);
    double ddis_d = mu*BURGER*one_over_L;
    tau[count] = tau_RSS[count] - sign(tau_RSS[count])*(dtau_f+ddis_d);
    
    if(abs(tau_RSS[count])>abs(dtau_f)+abs(ddis_d)){
        return true;
    }else{
        return false;
    }
    
}

void DeformationGradient::computeOneRho_f(int& count){
    double sum = 0.0;
    double sum_1 = 0.0;
    Vector3d n_1 = slip_systems.getSlipPlaneNormal(count); // slip normal alpha
    for(int i = 0; i < 12; i++){
        Vector3d s = slip_systems.getSlipDirection(i);
        Vector3d n = slip_systems.getSlipPlaneNormal(i);
        sum += 2*(Rho_e[i]*abs(n_1.dot(s.cross(n)))+Rho_s[i]*abs(n_1.dot(s)));
        sum_1 += 2*xi[count][i]*(Rho_e[i] + Rho_s[i]);
    }
    Rho_f[count] = sum;
    Rho_int[count] = sum_1;
    //cout<<"Rho_f = "<< Rho_f[count]<<endl;
}
/*
 void DeformationGradient::computeOneLambda(int& count){
 Lambda[count] = 1.0/(1.0/d_g + sqrt(Rho_f[count]));
 }
 */
void DeformationGradient::computeOneDRho(int& count){
    Rho_tot = 0.0;
    for(int i = 0; i < 12; i++){
        Rho_tot += Rho_e[i];
        Rho_tot += Rho_s[i];
    }
    //dRho_e[count] = (2*Rho_s[count]*abs(V_s[count])/lambda_alpha[count]-Rho_e[count]*Rho_e[count]*BURGER*(abs(V_e[count])+abs(V_c[count])))*DT;
    //dRho_s[count] = (2*Rho_e[count]*abs(V_e[count])/lambda_alpha[count]-Rho_s[count]*Rho_s[count]*BURGER*abs(V_s[count]))*DT;
    dRho_s[count] = ((ABSRho/DT/12.0+2*Rho_e[count]/lambda_alpha[count]*abs(V_e[count]))-Rho_s[count]*Rho_s[count]*BURGER*abs(V_s[count])/4.0)*DT;
    dRho_e[count] = ((ABSRho/DT/12.0+abs(V_e[count])*sqrt(Rho_e[count])/lambda_alpha[count]/lambda_alpha[count])-Rho_e[count]*Rho_e[count]*BURGER*(abs(V_e[count])+abs(V_c[count])))*DT;
    //dRho_e[count] = (ABSRho/DT/12.0-Rho_e[count]*Rho_e[count]*2*BURGER*(abs(V_e[count])+abs(V_c[count])))*DT;
    //dRho_e[count] = ((ABSRho/DT/12.0+abs(V_e[count])*sqrt(Rho_e[count])*Rho_e[count]*BURGER/lambda_alpha[count])-Rho_e[count]*Rho_e[count]*BURGER*(abs(V_e[count])))*DT;
    //dRho_e[count] = ((ABSRho/DT/12.0+abs(V_e[count])*sqrt(Rho_tot)*Rho_e[count]*BURGER/lambda_alpha[count])-Rho_e[count]*Rho_e[count]*4*BURGER*(abs(V_e[count])))*DT;
    //cout<<lambda_alpha[count]<<endl;
    //cout<<"***********************"<<endl;
    //cout<<"Delta edge dislocation = "<< dRho_e[count] << endl;
    //cout<<"Rho_s = "<<Rho_s[count]<<endl;
    //cout <<"V_s = "<<V_s[count]<<endl;
    //cout <<"lambda = "<<lambda_alpha[count]<<endl;
    //cout <<"Rho_e = "<<Rho_e[count]<<endl;
    //cout<<"***********************"<<endl;
    //getchar();
    //cout<<"screw dislocation"<<dRho_s[count]<<endl;
    
}

/*
 void DeformationGradient::computeOneXi(int& count){
 double tau_101 = (Cauchy_stress.array()*slip_systems.getSCCN(count).array()).sum();
 double tau0_11 = (Cauchy_stress.array()*slip_systems.getSCCN1(count).array()).sum();
 double ptau_101= (Cauchy_stress.array()*slip_systems.getNCSCCN(count).array()).sum();
 double ptau0_11 = (Cauchy_stress.array()*slip_systems.getN1CSCCN1(count).array()).sum();
 Xi[count] = (tau0_11 + a_1*tau0_11)/(a_0*tau_p-a_2*ptau0_11-a_3*ptau0_11);
 if(Xi[count] <= 0.0){
 Xi[count] = abs(Xi[count]);
 }else if(Xi[count] > 1.0){
 Xi[count] = Xi[count]-1.0;
 }
 }
 */
void DeformationGradient::computeOneB_drag(int& count){
    /*
     deltaG_kp[count] = Delta_H0*(pow(abs(1-pow(Xi[count],p)),q)-T/T0);
     if(deltaG_kp[count]>0){
     
     B_drag[count] = B_k*ALATT*(2*ALATT*exp(deltaG_kp[count]/2/KB/T)+l)/2/h/l;
     
     }else{
     
     B_drag[count] = B_0 + B_1*T;
     }
     */
    B_drag[count] = 1e-4; // Pa*s
}

void DeformationGradient::computeOneV(int& count){
    double dtau_f = mu*BURGER*sqrt(Rho_int[count]);
    double ddis_d = mu*BURGER*one_over_L;
    double dirr=1.0/one_over_L;
    lambda_alpha[count] = 1.0/(sqrt(Rho_f[count])+one_over_L);
    tau[count] = tau_RSS[count] - sign(tau_RSS[count])*(dtau_f+ddis_d);
    double B = 6.7e-7*TEMPERATURE; // Pa.s
    //cout << "d_irr = " << dirr << endl;
    //cout << "rbar = " << avg_r << endl;
    //cout << "vclimb = " << V_c[count]<<endl;
    //getchar();
    //double B = 8e-7*TEMPERATURE;
    if(abs(tau_RSS[count])>abs(dtau_f)+abs(ddis_d)){
        double w = 11.0*BURGER;
        double v_0 = sign(tau_RSS[count])*nu0*h/BURGER*(lambda_alpha[count]-w);
        /* compute glide velocity */
        //cout<<"v0 = "<<v_0<<endl;
        if(lambda_alpha[count] - w > 0){
            V_s[count] = v_0*exp(-Delta_H0/KB/TEMPERATURE*(pow((1-pow(abs(tau[count]/tau_p),p0)),q0)));
            
        }else{
            V_s[count] = 0.0;
            
        }
        if(tau_RSS[count] == 0.0){
            V_s[count] = 0.0;
        }
        double v1 = sign(tau_RSS[count])*BURGER*tau[count]*1e9/B;
        double v2 = 0.0;
        V_e[count] = sign(tau_RSS[count])*BURGER*tau[count]*1e9/B;
        
    }else{
        V_s[count] = 0.0;
        if(tau_RSS[count] == 0){
            V_e[count] = 0;
            
        }else{
            double numerator = BURGER*tau_RSS[count]*(lambda_alpha[count]+2*BURGER+avg_r)*V_c[count];
            double denomenator = (B*1e-9*V_c[count]*(lambda_alpha[count]+2*BURGER+avg_r)+BURGER*tau_RSS[count]*(2*BURGER+avg_r));
            if(numerator == 0.0){
                V_e[count] = 0.0;
            }else if(denomenator == 0.0){
                V_e[count] = 0.0;
            }else{
                V_e[count] = numerator/denomenator;
            }
        }
        
    }
    //cout<<"v_e = "<<V_e[count]<<endl;
    //fstream fve;
    //fve.open("ve.txt", ios::app);
    //fve << tau_RSS[count] << "  " << V_e[count] <<endl;
    //fve.close();
    //getchar();
}
/*
void DeformationGradient::computeOneVc(int& count){
    // compute climb velocity
    //N_0 = exp(-Hf/KB/TEMPERATURE)*VOLUME*DENSITY;
    //N_0 = 1;
    if(N_v == 0){
        N_v = 1;
    }
    double y = N_0/N_v; //number ratio
    //double y = N_v/N_0;
    double A = 1.9e+5-1.5e+3*y;
    double Q = 1.3; //[eV]
    cout<<"sigma_xx = "<<sigma_xx[count]<<endl;
    //V_c[count] = sign(y)*A*exp(-Q/KB/TEMPERATURE)*(KB*TEMPERATURE*1.602e-19/(ATOMICVOLUME*1e-27)*(1-N_v/N_0)-sigma_xx[count]*1e+9)*BURGER; //[cm/s]
    V_c[count] = sign(y)*A*exp(-Q/KB/TEMPERATURE)*(KB*TEMPERATURE*1.602e-19/(ATOMICVOLUME*1e-27)*(1-N_v/N_0)-sigma_xx[count]*1e+9)*BURGER; //[cm/s]
    //cout << "C0 = "<< N_0  <<endl;
    //cout << "A = " << A << endl;
    //cout<<KB*TEMPERATURE*1.602e-19/(ATOMICVOLUME*1e-27)<<endl;
    //cout<<"N_v/N_0 = "<<N_v/N_0<<endl;
    //cout<<"V_c = "<<V_c[count]<<endl;
}
*/
void DeformationGradient::computeOneVc(int& count){
    double totalRhoEdge = 0.0;
    V_c[count] = abs(Nv_Ni)*ATOMICVOLUME*Rho_e[count]/VOLUME/sumRho/DT/BURGER/sumRho;
    /*
    Rho_tot = sumRho;
    for(int i = 0; i<12; i++){
        totalRhoEdge += 2*Rho_e[i];
    }
    V_c[count] = (ATOMICVOLUME/VOLUME)*(abs(Nv_Ni)/BURGER)*(totalRhoEdge/Rho_tot)*(1.0/Rho_e[count]);
    */
    //cout << "dNabs = " << dNabs << endl;
    //V_c[count] = abs(dNabs)*ATOMICVOLUME*Rho_e[count]/BURGER/DT/VOLUME/sumRho;
    //cout<<"Rho_e = "<<Rho_e[count]<<endl;
    //cout<<"sumRho = "<<sumRho<<endl;
    //cout<<"Vc["<<count<<"] = "<<V_c[count]<<endl;
    //getchar();
}

void DeformationGradient::computeOneDGamma(int& count){
    dGamma[count] = (Rho_e[count]*V_e[count]+Rho_s[count]*V_s[count])*BURGER;
}

void DeformationGradient::computeOneDBeta(int& count){
    dBeta[count] = Rho_e[count]*V_c[count]*BURGER;
}
/*
 void DeformationGradient::computeLP(){
 L_P = Matrix3d::Zero();
 for(int i = 0; i<12; i++){
 Matrix3d tempo = dGamma[i]*slip_systems.getSCCN(i);
 if(tempo(1,1)<0){
 tempo(1,1) = abs(tempo(1,1));
 }
 L_P += tempo;
 }
 }
 */
//void DeformationGradient::computeL(){} // not until coupled with SCD

Matrix3d DeformationGradient::computeDCauchyStress(Matrix3d& strain_elastic){
    Matrix3d d_CauchyStress;
    strain_vec << strain_elastic(0,0), strain_elastic(1, 1), strain_elastic(2,2), 2.0*strain_elastic(1, 2), 2.0*strain_elastic(0, 2),2.0*strain_elastic(0,1);
    stress_vec = stiffness*strain_vec;
    d_CauchyStress << stress_vec(0), stress_vec(5), stress_vec(4),
    stress_vec(5), stress_vec(1), stress_vec(3),
    stress_vec(4), stress_vec(3), stress_vec(2);
    return d_CauchyStress;
}

void DeformationGradient::updateDislocationDensity(int& count){
    
    Rho_e[count] = Rho_e[count] + dRho_e[count];
    Rho_s[count] = Rho_s[count] + dRho_s[count];
}
/*
double DeformationGradient::putTempToRho(){
    
    double sum = 0.0;
    for(int i = 0; i< 12; i++){
        Rho[i] = tempRho[i];
        sum += tempRho[i];
    }
    return sum;
    
}
*/
void DeformationGradient::plotStressStrainTensors(){
    cout << "time = " << t << endl;
    cout << "dpa = " << 2.9e-7*t << endl;
    cout << "sigma = " << Cauchy_stress(1,1)<<endl;
    cout << "e_E = " << strain_elastic(1,1) << endl;
    //cout << "e_C = " << strain_climb << endl;
    //cout << "e_P = " << strain_plastic(1,1) << endl;
    //cout << "e_T = " << strain_tot(1,1) << endl << endl;
    cout << "e_swelling = " << volume_strain << endl;
    cout << "e_creep = " << effective_strain << endl << endl;
}

void DeformationGradient::initializeParameters(double& theta, double& phi){
    acctime = 0.0;
    h = ALATT*sqrt(3.0/2.0);
    round=0;
    t = 0.0;
    sumRho = 0.0;
    Cauchy_stress << 0.0, 0.0, 0.0,
                    0.0, STRESS*1e-3, 0.0,
                    0.0, 0.0, 0.0; //[GPa]
    //strain_elastic << 0.0, 0.0, 0.0,
    //                0.0, STRESS*1e-3/E, 0.0,
    //                0.0, 0.0, 0.0;
    stress_vec << Cauchy_stress(0,0), Cauchy_stress(1, 1), Cauchy_stress(2,2), 2.0*Cauchy_stress(1, 2), 2.0*Cauchy_stress(0, 2), 2.0*Cauchy_stress(0,1);
    strain_vec = stiffness.colPivHouseholderQr().solve(stress_vec);
    stress_vec = stiffness*strain_vec;
    Cauchy_stress << stress_vec(0), stress_vec(5), stress_vec(4),
    stress_vec(5), stress_vec(1), stress_vec(3),
    stress_vec(4), stress_vec(3), stress_vec(2);
    strain_elastic << strain_vec(0), strain_vec(5), strain_vec(4),
    strain_vec(5), strain_vec(1), strain_vec(3),
    strain_vec(4), strain_vec(3), strain_vec(2);
    cout<<"Cauchy_stress = "<<Cauchy_stress<<endl;
    cout<<"strain_elastic = "<<strain_elastic<<endl;
    getchar();
    strain_climb = 0.0;
    strain_plastic = Matrix3d::Zero();
    strain_tot = strain_elastic + strain_plastic;
    Th = theta;
    Ph = phi;
    slip_systems.updateSlipSystems(theta, phi);
    acctime = 0.0;
    //Hf = 1.7+STRESS*1e6*(1.673+1.0)*0.01179e-27*6.2415e+18;
    //cout<<"Hf = 1.7 + "<< STRESS*1e6*(1.673+1.0)*0.01179e-27*6.2415e+18 <<"= "<<Hf<<endl;
    // [eV] formation enthalpy Hf = Ef+p(Omega_rel + Omega_a)
    //double volume_dis = 1.5e11*VOLUME*PI*BURGER*BURGER;
    /* initialize dislocation densities etc. */
    //sumRho_edge = 0.0;
    for(int i = 0; i<12; i++){
        V_c[i] = 0.0;
        /* initialize dislocation densities */
        //Rho_e[i] = 1.0e10/12.0/4.0*1.1; // [cm^-2] initial dislocation density
        //Rho_s[i] = 1.0e10/12.0/4.0*1.1;
        Rho_e[i] = 5e10/12.0/4.0; // [cm^-2] initial dislocation density
        Rho_s[i] = 5e10/12.0/4.0;
        Rho_f[i] = 0.0;
        Rho_int[i] = 0.0;
        dRho_e[i] = 0.0;
        dRho_s[i] = 0.0;
        /* initialize velocity */
        V_e[i] = 0.0;
        V_s[i] = 0.0;
        /* initialize slip or climb rate */
        dGamma[i] = 0.0;
        dBeta[i] = 0.0;
        /* initialize resolved shear stress */
        Vector3d s = slip_systems.getSlipDirection(i);
        Vector3d n = slip_systems.getSlipPlaneNormal(i);
        tau_RSS[i] = s.transpose()*Cauchy_stress*n;
        sigma_xx[i] = s.transpose()*Cauchy_stress*s; //[GPa]
        tau[i] = 0.0;
        sumRho +=2*(Rho_e[i]+Rho_s[i]);
        //sumRho_edge += Rho_e[i];
    }
    //N_0 = exp(-Hf/KB/TEMPERATURE)*DENSITY*VOLUME;
    //cout<<"N_0 = "<<N_0<<endl;
    //getchar();
    Nv_Ni = 0;
    srscd.clearSRSCD();
    avg_r = 1.0;
    /* initialize one_over_L */
    one_over_L = 0.0; //0.0
    total_step = 0;
    //cout<< "Rho_0 = " << Rho[0] << endl << endl;
    //cout<<"total step = " << TotalStrain/dt/totalStrainRate<<endl;
    //cin.get();
}

double DeformationGradient::getTime(){
    return t;
}

double DeformationGradient::getDPA(){
    return srscd.getTotalDpa(t);
}
