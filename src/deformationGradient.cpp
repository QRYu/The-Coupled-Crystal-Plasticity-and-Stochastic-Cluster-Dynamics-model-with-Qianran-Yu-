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
    
    round = -1;
    maxiter = 50;
    acctime = 0.0;
    h = ALATT*sqrt(3.0/2.0);
    dt = DT;
    totalStrainRate = 0.0001;
    TOL = 1e-6*dt*totalStrainRate;
    tempPS = 0.0;
    unstable_time = 0;
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
    //L_P = Matrix3d::Zero();
    ft.open("time.txt", ios::app);
    ft<<"wall-clock(tot)    wal-clock(SCD)  wal-clock(CP)   simulation-time   step(scd)   step(cp)    scdrate(s/step)    cprate(s/step)"<<endl;
    ft.close();
    totalWallTime = 0.0;
    totalWalltimeCP = 0.0;
    totalWalltimeSCD = 0.0;
    dSCDTime = 0.0;
    dCPTime = 0.0;
    searchStepCP = 0;
    executeStepSCD = 0;
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
    if(t == 0.0){
        char strainStressFile[30];
        cout <<  "round = "<< round << endl;
        sprintf(strainStressFile, "strain_stress%d.txt", round);
        fss.open(strainStressFile, ios::app);
        fss << "theta = " << Th << "     phi = "<< Ph << endl;
    }
    dSCDTime = clock();
    executeStepSCD = srscd.doClusterDynamics(one_over_L, dt);
    dSCDTime = clock()-dSCDTime;
    Matrix3d stress_store = Cauchy_stress;
    Matrix3d delta_strain_T = computeDStrainT();
    double delta_strain_Pp = 0.0 ;
    Matrix3d delta_strain_P = delta_strain_T;
    //double delta_strain_Pp = delta_strain_T(1,1);
    //Matrix3d delta_strain_P = Matrix3d::Zero();
    Matrix3d delta_strain_E;
    Matrix3d delta_stress;
    bool whether;
    int n=1;
    dCPTime = clock();
    while(abs(delta_strain_P(1,1) - delta_strain_Pp)> TOL && n < maxiter){
        searchStepCP++;
        //cout<<"1"<<endl;
        delta_strain_P = computeDPlasticStrain(delta_strain_Pp);
        delta_strain_Pp = 0.0;
        delta_strain_E = delta_strain_T - delta_strain_P;
        delta_stress = computeDCauchyStress(delta_strain_E);
        Cauchy_stress = stress_store+delta_stress;
        for(int i = 0; i < 12; i++){
            computeOneRho_f(i);
            whether = computeOneRSS(i);
            //cout << "whether = " << whether << endl;
            if(whether == false){
                continue;
            }
            computeOneV(i);
            computeOneDGamma(i);
            computeOneDRho(i);
            updateTempDislocationDensity(i);
            double small_ep = abs((abs(dGamma[i])*slip_systems.getSCCN(i)*dt)(1,1));
            //double small_ep = (abs(dGamma[i])*slip_systems.getSCCN(i)*dt)(1,1);
            delta_strain_Pp += small_ep;
        }
        n +=1;
    }
    dCPTime = clock()-dCPTime;
    totalWalltimeSCD += dSCDTime;
    totalWalltimeCP += dCPTime;
    totalWallTime = totalWalltimeSCD + totalWalltimeCP;
    if(n == maxiter && delta_strain_Pp >(delta_strain_T(1,1)+TOL)){
        unstable_time++;
        //exit();
    }
    double sumRho = putTempToRho();
    t += dt;
    acctime +=dt;
    strain_tot += delta_strain_T;
    strain_plastic += delta_strain_P;
    strain_elastic += delta_strain_E;
    if(Cauchy_stress(1,1) > 5.0 || Cauchy_stress(1,1) < -1.0){
        ff.open("fracture_toughness.txt", ios::app);
        double x0 = (-cos(Th)+sqrt(cos(Th)*cos(Th)+1))*Ph*4.0/PI;
        ff << round << "    " << x0*cos(Th) << " "<< sin(Th)*x0 << " " << strain_tot(1,1) <<" " << stress_store(1,1) << endl;
        ff.close();
        return -1;
    }
    if((Cauchy_stress(1,1)-268.0*(strain_tot(1,1)-0.002))<0.001 && yield == 0){
        yield = 1;
        fy.open("yield_stress.txt", ios::app);
        double x0 = (-cos(Th)+sqrt(cos(Th)*cos(Th)+1))*Ph*4.0/PI;
        fy << round << "    " << x0*cos(Th) << "  " << sin(Th)*x0 << " " << strain_tot(1,1) <<" " << Cauchy_stress(1,1) << endl;
        fy.close();
    }
    srscd.updateRates(sumRho);
    if(acctime >= 0.1){
        acctime = 0.0;
        //cout<<"delta_strain_Pp = "<< delta_strain_Pp << endl;
        //plotStressStrainTensors();
        // write into files
        char strainStressFile[30], DDFile[30];
        sprintf(DDFile, "dd%d.txt", round);
        sprintf(strainStressFile, "strain_stress%d.txt", round);
        fss.open(strainStressFile, ios::app);
        fss << strain_tot(1,1) << "   " << Cauchy_stress(1,1) << endl;
        fss.close();
        fdd.open(DDFile, ios::app);
        fdd << t << "   " << strain_plastic(1,1)<< " "<< strain_tot(1,1)<<"  "<<sumRho <<"    ";
        for(int i = 0; i<12; i++){
            fdd << Rho[i] << "  ";
        }
        ft.open("time.txt", ios::app);
        ft << totalWallTime/CLOCKS_PER_SEC << "    " << totalWalltimeSCD/CLOCKS_PER_SEC << "    " << totalWalltimeCP/CLOCKS_PER_SEC << "  " << t << " " << executeStepSCD << "    " << searchStepCP<<"  " << totalWalltimeSCD/CLOCKS_PER_SEC/executeStepSCD << "  " << totalWalltimeCP/CLOCKS_PER_SEC/searchStepCP << endl;
        ft.close();
        //fdd<<endl;
        //fdd.close();
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

Matrix3d DeformationGradient::computeDStrainT(){
    Matrix3d d_strain;
    //d_strain<< -poisson*totalStrainRate*dt, 0.0, 0.0,
    //0.0, totalStrainRate*dt, 0.0,
    //0.0, 0.0, -poisson*totalStrainRate*dt;
    d_strain<< 0.0, 0.0, 0.0,
    0.0, totalStrainRate*dt, 0.0,
    0.0, 0.0, 0.0;
    return d_strain;
}

Matrix3d DeformationGradient::computeDPlasticStrain(double& principal){
    Matrix3d d_strain;
    //d_strain <<-poisson*principal, 0.0, 0.0,
    //0.0, principal, 0.0,
    //0.0, 0.0, -poisson*principal;
    d_strain <<0.0, 0.0, 0.0,
    0.0, principal, 0.0,
    0.0, 0.0, 0.0;
    return d_strain;
}

void DeformationGradient::computeStrainE(){
    
    strain_elastic = strain_tot - strain_plastic;
}

bool DeformationGradient::computeOneRSS(int& count){
    
    Vector3d s = slip_systems.getSlipDirection(count);
    Vector3d n = slip_systems.getSlipPlaneNormal(count);
    double RSS = s.transpose()*Cauchy_stress*n;
    double dtau_f = mu*BURGER*sqrt(Rho_int[count]);
    double ddis_d = mu*BURGER*one_over_L;
    tau[count] = RSS - sign(RSS)*(dtau_f+ddis_d);
    
    if(abs(RSS)>abs(dtau_f)+abs(ddis_d)){
        return true;
    }else{
        return false;
    }
    
}

void DeformationGradient::computeOneRho_f(int& count){
    double sum = 0.0;
    double sum_1 = 0.0;
    for(int i = 0; i < 12; i++){
        sum += Rho[i]*abs(slip_systems.getSlipPlaneNormal(count).dot(slip_systems.getSlipDirection(i)));
        sum_1 += xi[count][i]*Rho[i];
        //cout<< "Rho = "<<Rho[i]<<endl;
    }
    Rho_f[count] = sum;
    Rho_int[count] = sum_1;
    //cout<< "count = "<< count << "  Rho_f = "<< Rho_f[count] <<"    Rho_int = " << Rho_int[count] << endl;
}
/*
 void DeformationGradient::computeOneLambda(int& count){
 Lambda[count] = 1.0/(1.0/d_g + sqrt(Rho_f[count]));
 }
 */
void DeformationGradient::computeOneDRho(int& count){
    double lambda = sqrt(Rho_f[count]) + one_over_L;
    dRho[count] = abs(dGamma[count])*dt/BURGER*(lambda-2*BURGER*Rho[count]);
    //dRho[count] = dGamma[count]*dt/BURGER*(sqrt(Rho_f[count])-2*BURGER*Rho[count]);
    //cout<<"dRho = "<< dRho[count] << endl;
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
    
    Vector3d s = slip_systems.getSlipDirection(count);
    Vector3d n = slip_systems.getSlipPlaneNormal(count);
    double RSS = s.transpose()*Cauchy_stress*n;
    double w = 11.0*BURGER;
    double lambda_alpha = 1.0/(sqrt(Rho_f[count])+one_over_L);
    
    double v_0 = sign(tau[count])*nu0*h/BURGER*(lambda_alpha-w);
    //cout << "p0 = " << p0 << " , q0 = " << q0 <<endl;
    if(abs(tau[count])<=tau_p){
        V[count] = v_0*exp(-Delta_H0/KB/TEMPERATURE*(pow((1-pow(abs(tau[count]/tau_p),p0)),q0)));
    }else{
        computeOneB_drag(count);
        V[count] = (tau[count]-sign(RSS)*tau_p)*1e9*BURGER/B_drag[count];
    }
    
    
}

void DeformationGradient::computeOneDGamma(int& count){
    dGamma[count] = Rho[count]*BURGER*V[count];
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

void DeformationGradient::updateTempDislocationDensity(int& count){
    
    tempRho[count] = Rho[count] + dRho[count];
    
}

double DeformationGradient::putTempToRho(){
    
    double sum = 0.0;
    for(int i = 0; i< 12; i++){
        Rho[i] = tempRho[i];
        sum += tempRho[i];
    }
    return sum;
    
}

void DeformationGradient::plotStressStrainTensors(){
    cout << "sigma = " << Cauchy_stress(1,1)<<endl;
    cout << "e_E = " << strain_elastic(1,1) << endl;
    cout << "e_P = " << strain_plastic(1,1) << endl;
    cout << "e_T = " << strain_tot(1,1) << endl << endl;
}

void DeformationGradient::initializeParameters(double& theta, double& phi){
    round++;
    yield = 0;
    t = 0.0;
    strain_tot = Matrix3d::Zero();
    strain_elastic = Matrix3d::Zero();
    strain_plastic = Matrix3d::Zero();
    Cauchy_stress = Matrix3d::Zero();
    Th = theta;
    Ph = phi;
    slip_systems.updateSlipSystems(theta, phi);
    acctime = 0.0;
    /* initialize dislocation densities etc. */
    for(int i = 0; i<12; i++){
        Rho[i] = 1.0e10/12.0; // [cm^-2] initial dislocation density
        tempRho[i] = 1.0e10/12.0;
        Rho_f[i] = 0.0;
        V[i] = 0.0;
        dGamma[i] = 0.0;
        dRho[i] = 0.0;
        disU[i] = 0;
    }
    srscd.clearSRSCD();
    /* initialize one_over_L */
    one_over_L = 0.0; //0.0
    total_step = 0;
    cout<< "Rho_0 = " << Rho[0] << endl << endl;
    cout<<"total step = " << TotalStrain/dt/totalStrainRate<<endl;
    //cin.get();
}

