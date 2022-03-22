#include "Object.h"
#include<cmath>
// Object.cpp -- implementations of Object class

/* public method implementations */
Object::Object(
               const int64 & key,
               const int & count,
               double& dd,
               const int& n) :oKey(key), totalNumber(0), bindSH(0.0)
{
    DISLOCATION = dd;
    setAttributes(key);
    setProperties(count, n);
}

Object::Object(const int * attr,
               const int& count,
               double& dd,
               const int& n): totalNumber(0)
{
    DISLOCATION = dd;
    for (int i = 0; i < LEVELS; i++) {
        attributes[i] = attr[i];
    }
    oKey = 0;
    setKey();
    setProperties(count, n);
}

Object::Object(const int64 &key, double& dd, const int *number):oKey(key), totalNumber(0)
{
    DISLOCATION = dd;
    setAttributes(key);
    dimensionality = setDimensionality();
    computeDiffCoeff();
    computeBindTerm();
    computeR1R1e();
    computeSinks();
    setNumber();
    for (int i = 0; i < POINTS; i++) {
        addNumber(i, number[i]);
    }
    
}


void Object::addNumber(const int & count, const int& n)
{
    number[count] += n;
    totalNumber += n;
}

void Object::reduceNumber(const int & count)
{
    --number[count];
    --totalNumber;
}

int Object::signof(const int64 & key) const
{
    return (key < 0) ? -1 : 1;
}

double Object::zero(const int & defect)
{
    return (abs(defect) > 1) ? 1.0 : 0.0;
}

int64 Object::getKey() const
{
    return oKey;
}

double Object::getDiff() const
{
    return diffusivity;
}


int Object::getNumber(const int & count) const
{
    return number[count];
}

int Object::getTotalNumber() const
{
    return totalNumber;
}

int Object::getAttri(const int & index) const
{
    return attributes[index];
}

double Object::getSink() const
{
    return sinkStrength;
}

long double Object::getBind(const int & index) const
{
    return bind[index];
}

double Object::getR1() const
{
    return r1;
}

double Object::getR1e() const
{
    return r1e;
}

int Object::getDim() const
{
    return dimensionality;
}

double Object::getBindSH() const
{
	return bindSH;
}

void Object::getThreeNumber(const int & count, int* objectN) const
{
    /*
     * objectN[0] = object number in this element
     * objectN[1] = object number in the previous element
     * objectN[2] = object number in the next element
     */
    objectN[0] = number[count];
    if (count == 0) {
        /* when this is the surface element */
        objectN[1] = 0; /* Question, ask Jaime */
        objectN[2] = number[count + 1];
    }
    else if (count == POINTS - 1) {
        /* when this is the last element*/
        objectN[1] = number[count - 1];
        objectN[2] = 0; /* Question, ask Jaime */
    }
    else {
        objectN[1] = number[count - 1];
        objectN[2] = number[count + 1];
    }
}

void Object::display()
{
    cout << "Information of Object " << oKey << ": " << endl;
    cout << "Attributes: ";
    for (int i = 0; i < LEVELS; ++i) {
        cout << attributes[i] << "    ";
    }
    cout << endl;
    cout << "number: ";
    for (int i = 0; i < POINTS; ++i) {
        cout << number[i] << "    ";
    }
    cout << endl;
    cout << "total number: " << totalNumber << endl;
    cout << "dimensionality: " << dimensionality << endl;
    cout << "diffusivity: " << diffusivity << endl;
    cout << "bind term: ";
    for (int i = 0; i < LEVELS; ++i) {
        cout << bind[i] << "    ";
    }
    cout << endl;
    cout << "r1 = " << r1 << endl;
    cout << "r1e = " << r1e << endl;
    cout << "sink strength: " << sinkStrength << endl;
    cout << endl;
}

/* private method impementations */
void Object::setKey()
{
    for (int i = 0; i < LEVELS; ++i) {
        oKey += labs(attributes[i])*((int64)pow(10.0, (double)EXP10*(LEVELS - 1 - i)));
    }
    oKey *= signof(attributes[0]);
}

void Object::setAttributes(const int64 & key)
{
    int64 tempKey = abs(key);
    for (int i = 0; i < LEVELS; i++) {
        attributes[i] = double(tempKey) / pow(10.0, (double)EXP10*(LEVELS - 1 - i));
        tempKey -= ((int64)attributes[i])*((int64)pow(10.0, (double)EXP10*(LEVELS - 1 - i)));
    }
    attributes[0] *= signof(key);
}

void Object::setNumber()
{
    for (int i = 0; i < POINTS; i++) {
        number[i] = 0;
        totalNumber += number[i];
    }
}

int Object::setDimensionality()
{
    return attributes[0] > 4 ? 1 : 3;
}

void Object::computeR1R1e()
{
    int ndef = attributes[0];
    if (ndef <= 0) {
        r1 = zero(ndef)*pow(3.0*fabs((double)ndef)*avol / 4.0 / PI, 0.333333333333333333) + jumped;
        if (ndef != 0)
            r1e = pow(3.0*(fabs((double)ndef) - 1)*avol / 4.0 / PI, 0.333333333333333333) + jumped;
        else
            r1e = jumped;
    }
    else if (ndef > 0) {
        r1 = zero(ndef)*sqrt((double)ndef*avol / jumped / PI) + jumped;
        r1e = zero(ndef)*sqrt(((double)ndef - 1)*avol / jumped / PI) + jumped;
    }
}

void Object::computeDiffCoeff()
{
    const double fi = 0.9, fv = 0.7; // Diffusion correlationm factors.
    const double gi = 0.5, gv = 0.125; // Geometric factor for diffusion.
    double prefactor = 0, energy_m = 0;
    int check_all = 0;
    int check_He = 0;
    int check_H = 0;
    // int check_C = 0;
    int i;
    
    for (i = 1; i < LEVELS; i++) {
        check_all |= attributes[i];
        if (i >= 2) check_He |= attributes[i];
        if (i >= 3) check_H |= attributes[i];
    }
    /* check_all: 0 when this is a pure defect cluster.
     *  check_He:  0 when this is a cluster with He.
     *  check_H:   0 when this is a cluster with H.
     */
    
    // Pure defect clusters:
    if (!check_all) {
        
        if (attributes[0] > 0) { // SIAs
            if (abs(attributes[0]) == 1) { // 1I
                prefactor = 8.744e-4;
                energy_m = 0.009;
            }else if (abs(attributes[0]) == 2) { // 1I
                prefactor = 7.97e-4;
                energy_m = 0.024;
            }else if (abs(attributes[0]) == 3) { // 1I
                prefactor = 0.0;
                //prefactor = 3.92e-4;
                energy_m = 0.033;
            }
            else if(abs(attributes[0]) > 3) { // >1I
                //prefactor = gi*jumped*jumped*fi*NU0*pow(fabs(attributes[0]), -0.5);
                prefactor = 0.0;
                energy_m = 0.013;
            }
        }
        else if (attributes[0] < 0) { // Vacancies.
            if (abs(attributes[0]) == 1) { // 1V
                
                prefactor = 1.77e-2;
                energy_m = 1.29;
            }
            else if (abs(attributes[0]) == 2) { // >1V
                
                prefactor = 2.91e-5;
                energy_m = 1.66;
            }else if (abs(attributes[0]) > 2) { // >1V
                
                //prefactor = gv*jumped*jumped*fv*NU0*pow(0.001, fabs(attributes[0]) - 1.0);
                prefactor = 0.0;
                energy_m = 1.66;
            }
        }
    }
    else if (!check_He) {
        
        if (attributes[0] != 0) { // SIA-He and V-He clusters immobile.
            prefactor = 0.0;
        }
        else if (attributes[1] == 1) { // He1.
            prefactor = gi*jumped*jumped*NU0;
            energy_m = 0.01;
        }
        else  if (attributes[1] == 2) { // He2.
            prefactor = gi*jumped*jumped*NU0*0.01;
            energy_m = 0.03;
        }
        else  if (attributes[1] == 3) { // He3.
            prefactor = gi*jumped*jumped*NU0*0.01;
            energy_m = 0.05;
        }
        else  if (attributes[1] > 3) { // He>3.
            prefactor = gi*jumped*jumped*NU0*0.01;
            energy_m = 0.06;
        }
        else prefactor = 0.0;
    }
    else if (!check_H) {
        if (attributes[0] != 0) { // SIA-H and V-H clusters immobile.
            prefactor = 0.0;
        }
        else if (attributes[2] == 1 || attributes[2] == 2) { // H1, H2.
            /*
             //number below is from ppt Liu(2015) at 300K D_H=10e-9 m^2/s = 10e-5 cm^2/s
             prefactor = 3.8e-3;
             energy_m = 0.41;
             */
            //prefactor = 1.58e-3;
            prefactor = 0.0;
            energy_m = 0.25;
        }
        else
            prefactor = 0.0;
    }
    /* All data from [CS Becquart et al., J Nucl Mater 403 (2010) 75] */
    diffusivity = prefactor*exp(-energy_m / KB / TEMPERATURE);
}

void Object::computeBindTerm()
{
    long double energy_d[LEVELS] = { 0.0 };
    long double energy_b = 0.0, energy_m = 0.0;
    double attfreq = 1.0;
    double efi = 9.96, emi = 0.013; // Ab initio migration and formation energies of V and SIA in pure W.
    double efv = 3.23, emv = 1.66;
    double eb2v = -0.1, eb2i = 2.12, eb2he = 1.03;
    double efhe = 4.0, emhe = 0.01;
    double emh = 0.39;
    int check_all = 0, check_He = 0, check_H = 0;
    int i;
    
    for (i = 0; i < LEVELS; i++) {
        bind[i] = 0.0;
        if (i >= 1) check_all |= attributes[i];
        if (i >= 2) check_He |= attributes[i];
        if (i >= 3) check_H |= attributes[i];
    }
    /**
     * bind energy positive(eg. 5 eV) means easy to get together, hard to dissociate, when dissociate, absorb 5eV energy
     * bind energy negative(eg. -5 eV) means hard to get together, easy to dissociate, when dissociate, release 5eV energy
     **/
    // Pure defect clusters:
    if (!check_all) {
        if (attributes[0]>0) { // SIAs
            energy_m = emi;
            if (abs(attributes[0]) == 1) { // 1I
                attfreq = 0.0;
            }
            else if (abs(attributes[0]) == 2) { // 2I
                energy_b = 2.12;
            }
            else if (abs(attributes[0]) == 3) { // 3I
                energy_b = 3.02;
            }
            else if (abs(attributes[0]) == 4) { // 4I
                energy_b = 3.60;
            }
            else if (abs(attributes[0]) == 5) { // 5I
                energy_b = 3.98;
            }
            else if (abs(attributes[0]) == 6) { // 6I
                energy_b = 4.27;
            }
            else if (abs(attributes[0]) == 7) { // 7I
                energy_b = 5.39;
            }
            else if (abs(attributes[0])>7) // > 7I
                energy_b = efi + (eb2i - efi)*(pow(fabs((double)attributes[0]), 0.6666667) - pow((fabs((double)attributes[0]) - 1.0), 0.6666667)) / 0.5847;
        }
        else if (attributes[0]<0) { // Vacancies.
            energy_m = emv;
            if (abs(attributes[0]) == 1) { // 1V
                attfreq = 0.0;
            }
            else if (abs(attributes[0]) == 2) { // 2V
                energy_b = eb2v;
            }
            else if (abs(attributes[0]) == 3) { // 3V
                energy_b = 0.04;
            }
            else if (abs(attributes[0]) == 4) { // 4V
                energy_b = 0.64;
            }
            else if (abs(attributes[0]) == 5) { // 5V
                energy_b = 0.72;
            }
            else if (abs(attributes[0]) == 6) { // 6V
                energy_b = 0.89;
            }
            else if (abs(attributes[0]) == 7) { // 7V
                energy_b = 0.72;
            }
            else if (abs(attributes[0]) == 8) { // 8V
                energy_b = 0.88;
            }
            else if (abs(attributes[0])>8) // > 8V
                energy_b = efv + (eb2v - efv)*(pow(fabs((double)attributes[0]), 0.6666667) - pow((fabs((double)attributes[0]) - 1.0), 0.6666667)) / 0.5874;
        }
        energy_d[0] = energy_b; // + energy_m
        bind[0] = attfreq*exp(-energy_d[0] / KB / TEMPERATURE);
    }
    // He-defect clusters:
    else if (!check_He) {
        
        if (attributes[0]<0) { // He-V clusters:
            double ratio = fabs(((double)attributes[1]) / ((double)attributes[0]));
            //printf("%dV - %dHe\n", abs(attributes[0]), abs(attributes[1]));
            if (abs(attributes[0]) == 1 && abs(attributes[1]) == 1) { // 1V-1He
                energy_d[0] = 4.6;
                energy_d[1] = energy_d[0];
            }
            else {
                energy_d[0] = 2.4 + 3.5*log10(ratio) + 1.7*log10(ratio)*log10(ratio);
                // binding energy of V to cluster.
                energy_d[1] = 4.6 - 1.1*log10(ratio) - 0.3*log10(ratio)*log10(ratio);
                // binding energy of He to cluster.
            }
            bind[0] = attfreq*exp(-energy_d[0] / KB / TEMPERATURE);
            bind[1] = attfreq*exp(-energy_d[1] / KB / TEMPERATURE);
        }
        else if (attributes[0]>0) // He-SIA clusters.
            attfreq = 0.0; // No dissociation between He and SIA clusters.
        else if (attributes[0] == 0) { // pure He clusters.
            if (attributes[1] == 1) { // He1.
                attfreq = 0;
            }
            else if (attributes[1] == 2) { // He2.
                energy_b = 1.03;
            }
            else if (attributes[1] == 3) { // He3.
                energy_b = 1.36;
            }
            else if (attributes[1] == 4) { // He4.
                energy_b = 1.52;
            }
            else { // He>4.
                energy_b = efhe + (eb2he - efhe)*(pow(fabs((double)attributes[0]), 0.6666667) - pow((fabs((double)attributes[0]) - 1.0), 0.6666667)) / 0.5874;
            }
            energy_d[1] = energy_b + emhe;
            bind[1] = attfreq*exp(-energy_d[1] / KB / TEMPERATURE);
        }
    }
    
    // H-defect clusters:
    // H-defect clusters:
    else if (!check_H){
        if (attributes[0]<0) { // H-V clusters:
            double ratio = fabs( ((double) attributes[2])/((double) attributes[0]) );
            attfreq = 0.0;
            //printf("%dV - %dH\n", abs(attr[0]), abs(attr[2]));
            /*this part is for binding energy of mV-nH that try to dissociate a V from Li Xiaochun(2015)*/
            if (ratio==1 /*&& attributes[0]> -5 && attributes[2]< 15*/) { // V-H.
                energy_b = 1.60;
            } else if (ratio==2) { // V-H2
                energy_b = 2.10;
            } else if (ratio==3) { // V-H3
                energy_b = 2.40;
            } else if (ratio==4) { // V-H4
                energy_b = 3.90;
            } else if (ratio==5) { // V-H5
                energy_b = 4.35;
            } else if (ratio==6) { // V-H6
                energy_b = 5.80;
            } else if (ratio==7) { // V-H7
                energy_b = 7.0;
            } else if(ratio == 8){ // V-H8
                energy_b = 8.25;
            } else if(ratio == 9){ // V-H9
                energy_b = 9.80;
            } else if(ratio == 10){ // V-H10
                energy_b = 11.25;
            } else if(ratio>10){
                
                energy_d[0] = 1.91 + 0.0974 * ratio * ratio; //extrapolation
                // if happened form this, dissociate as soon as possible
            }
            energy_d[0] = energy_b + emv;
            /*this part is for binding energy of mV-nH+1H from Ohsawa(2015)*/
            if (ratio == 1 ) { // V-H.
                energy_b = 1.223;
            } else if (ratio==2 ) { // V-H2
                energy_b = 1.192;
            } else if (ratio==3 ) { // V-H3
                energy_b = 1.077;
            } else if (ratio==4 ) { // V-H4
                energy_b = 0.958;
            } else if (ratio==5 ) { // V-H5
                energy_b = 0.862;
            } else if (ratio==6 ) { // V-H6
                energy_b = 0.645;
            } else if (ratio==7 ) { // V-H7
                energy_b = 0.243;
            } else if(ratio == 8 ){ // V-H8
                energy_b = 0.306;
            } else if(ratio == 9 ){ // V-H9
                energy_b = 0.191;
            } else if(ratio == 10){ // V-H10
                energy_b = 0.184;
            } else if(ratio == 11){ // V-H11
                energy_b = 0.004;
            } else if(ratio == 12){ // V-H12
                energy_b = 0.367;
            } else if(ratio == 13){ // V-H13
                energy_b = -0.889;
            } else if(ratio == 14){ // V-H14
                energy_b = 0.747;
            } else if(ratio == 15){ // V-H15
                energy_b = -1.242;
            } else{
                energy_b = 1.45 - 0.125*ratio-0.00159*ratio*ratio; //extrapolation
                // energy_b = (numeric_limits<double>::min)(); // if happenedly get this large make it to disappear
            }
            energy_d[2] = energy_b + emh;
            bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
            bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
        }
        else if (attributes[0]>0){ // H-SIA clusters.
            energy_d[0] = 0.67 + emi;
            attfreq = 0.0;
            if (attributes[0]==1 && attributes[2]==1){ // SIA-H
                energy_b = 0.67;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[0]==1 && attributes[2]==2){ // SIA-H2
                energy_b = 0.40;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[0]==1 && attributes[2]==3){ // SIA-H3
                energy_b = 0.40;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[0]==1 && attributes[2]==4){ // SIA-H4
                energy_b = 0.05;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
                
            } else if (attributes[0]==1 && attributes[2]==5){ // SIA-H5
                energy_b = 0.20;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[0]==2 && attributes[2]==1){ // SIA2-H
                energy_d[0] = 2.12;
                energy_b = 0.57;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
                
            } else if (attributes[0]==2 && attributes[2]==2){ // SIA2-H2
                energy_d[0] = 2.12;
                energy_b = 0.45;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[0]==2 && attributes[2]==3){ // SIA2-H3
                energy_d[0] = 2.12;
                energy_b = 0.1;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[0]==2 && attributes[2]==4){ // SIA2-H4
                energy_d[0] = 2.12;
                energy_b = 0.3;
                energy_d[2] = energy_b + emh;
                bind[0] = attfreq*exp(-energy_d[0]/KB/TEMPERATURE);
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else{
                int H = attributes[2];
                energy_d[0] = 2.12;
                energy_d[2] = H + 0.013 * H * H * H * H - 0.44 * H * H;
                /* this part is unknown */
                
            }
            
        }else if (attributes[0]== 0) { // nH clusters, this is binding energy of nH cluster dissociating 1 H from Qin(2015)
            /*
             if (attr[2]==1) { // H
             energy_b = 0;
             } else if (attr[2]==2) { // 2H
             energy_b = 0.02;
             } else if (attr[2]==3) { // 3H
             energy_b = 0.08;
             } else if (attr[2]==4) { // 4H
             energy_b = 0.20;
             } else if (attr[2]==5) { // 5H
             energy_b = 0.27;
             }
             energy_d[2] = energy_b + emh;
             if (attr[2] > 5) // nH(n>5) or higher
             energy_d[2] = 0.0;
             */
            bind[0] = 0; //because there's no V/SIA in cluster
            attfreq = 0.0;
            if (attributes[2]==1) { // H   data from Xiaochun Li(2015)
                bind[2] = 0;
                
            } else if (attributes[2]==2) { // 2H
                bindSH = attfreq * exp(-(0.01 + emh)/KB/TEMPERATURE);
                /* on average, binding energy at surface is 0.01 eV*/
                energy_b = -0.12;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[2]==3) { // 3H
                energy_b = -0.1;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[2]==4) { // 4H
                energy_b = 0.20;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if (attributes[2]==5) { // 5H
                energy_b = -0.20;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            } else if(attributes[2]==6){
                energy_b = -0.30;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            }else if(attributes[2]==7){
                energy_b = -0.45;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            }else if(attributes[2]==8){
                energy_b = -0.15;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            }else if(attributes[2]==9){
                energy_b = 0.2;
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            }else if(attributes[2]==10){
                energy_b = -1.0; /* at this point SAV happen */
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                
            }else{
                int a = attributes[2];
                energy_b = 0.119 + 0.0407/(sin(7.94*a)) + 0.004*a*a*sin(7.93*a) - 0.104*a*sin(7.99*a)*sin(7.94*a);
                energy_d[2] = energy_b + emh;
                bind[2] = attfreq*exp(-energy_d[2]/KB/TEMPERATURE);
                /*extrapolation*/
                
            }
            
        }
    }
}

void Object::computeSinks()
{
    /* The total sink strength for all defects are stored in the array s.
     [0] for vacancies; [1] for SIAs; */
    double Zdv = 1.0, Zdi = 1.1;
    double Zodsv = 0.0, Zodsi = 0.0;
    double Zgbv = 1.0, Zgbi = 1.0;
    double Zsv = 1.0, Zsi = 1.1;
    double Sd, Sods, Sgbv, Sgbi, Sf;
    double s[2] = { 0 };
    /* 1. Dislocation sink strength: */
    Sd = DISLOCATION*exp(-1*KB*(TEMPERATURE - 300.0));
    
    /* 2. ODS-particle sink strength: */
    Sods = 4 * PI * ODS_R * ODS_DENSITY;
    
    /* 3. Grain boundary sink strength: */
    Sgbv = 6*sqrt(Zdv*Sd + Zodsv*Sods)/GRAIN_SIZE;
    Sgbi = 6*sqrt(Zdi*Sd + Zodsi*Sods)/GRAIN_SIZE;
    
    /* when internal sinks are weak */
    Sgbi = 24 / GRAIN_SIZE / GRAIN_SIZE;
    Sgbv = 24 / GRAIN_SIZE / GRAIN_SIZE;
    
    /* 4. Thin foil/interface: */
    // Sf = (2*sqrt(Zdi*5d + Zodsi*Sods)/FOIL_THICKNESS)*coth(sqrt(Zdi*5d + Zodsi*Sods)*FOIL_THINKNESS/4); */
    /* When internal sinks are weak */
    //Sf = 8 / FOIL_THICKNESS / FOIL_THICKNESS;
    //s[0] = Zdv * Sd + Zgbv * Sgbv + Zodsv * Sods + Zsv * Sf;
    //s[1] = Zdi * Sd + Zgbi * Sgbi + Zodsi * Sods + Zsi * Sf;
    
    //s[0] = Zdv * Sd + Zgbv * Sgbv + Zodsv * Sods;
    //s[1] = Zdi * Sd + Zgbi * Sgbi + Zodsi * Sods;
    s[0] = Zdv * Sd;
    s[1] = Zdi * Sd;
    sinkStrength = attributes[0] < 0 ? s[0] : s[1];
}

void Object::setProperties(const int & count, const int & n)
{
    setNumber();
    addNumber(count, n);
    dimensionality = setDimensionality();
    computeDiffCoeff();
    computeBindTerm();
    computeR1R1e();
    computeSinks();
}

void Object::updateDislocation(double& dd){
    DISLOCATION = dd;
}
