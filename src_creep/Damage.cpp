#include "Damage.h"
// Damage.cpp --  implementations of the damage class

Damage::Damage()
{
    int index;
    for (index = 0; index < POINTS; ++index) {
        totalRate[index] = 0.0;
    }
    readFile();
    //for(index=0; index<1; ++index){
    for (index = 0; index < POINTS; ++index) {
        computeDamageZero(index);
        if (CHANNELS > 1) {
            computeDamageOne(index);
        }
        if (CHANNELS > 2) {
            computeDamageTwo(index);
        }
        if (CHANNELS > 3) {
            int iindex;
            for (iindex = 3; iindex < CHANNELS; ++iindex) {
                computeDamageOther(index, iindex);
            }
        }
    }
}


Reaction Damage::selectDamage(const int & n, double & randRate)
{
    int index = 0;
    double tempRate = randRate;
    if (totalRate[n] < randRate) {
        randRate -= totalRate[n];
        return NONE;
        /* this block will never be excuted */
    }
    while (index < CHANNELS) {
        if (damage[n][index] >= tempRate) {
            if (index == 0) {
                return PARTICLE;
            }
            else if (index == 1) {
                return HE;
            }
            else if (index == 2) {
                return H;
            }
        }
        else {
            tempRate -= damage[n][index];
            ++index;
        }
    }
    return NONE;
}

void Damage::display(const int &count) const
{
    cout << "Damage information in element " << count << endl;
    for (int i = 0; i < CHANNELS; i++) {
        cout << "Damage[" << i << "]: " << damage[count][i] << endl;
    }
    cout << "Total Damage Rate: " << totalRate[count] << endl;
    //getchar();
}

const double Damage::getTotalDamage(const int & n) const
{
    return totalRate[n];
}

/* private method implementation */
void Damage::readFile()
{
    int index = 0;
    ifstream fd;
    fd.open("damage.txt");
    string oneLine;
    stringstream lineHold;
    while (getline(fd, oneLine)&& index < POINTS) {
        lineHold.str(oneLine);
        lineHold >> DPA_RATE[index] >> NRT[index];
        lineHold.clear();
        ++index;
    }
    fd.close();
}

void Damage::computeDamageZero(const int & n)
{
    //damage[n][0] = 0.0; /* after damage */
    /*
    if(n != 0){
        if( NRT[n] == 0.0 ){
            
            damage[n][0] = 0.0;
        }else{
            damage[n][0] = (DPA_RATE[n] * DENSITY*VOLUME / NRT[n]);
        }
        
    }else{
        damage[n][0] = 0.0;
    }*/
    double r_pka = 4.33e+14;
    //double nrt = 153.82; //Gilbert et al(2015)
    damage[n][0] = r_pka*VOLUME;
    //damage[n][0] = DPA_RATE[n] * DENSITY*VOLUME / NRT[n];
    totalRate[n] += damage[n][0];
}

void Damage::computeDamageOne(const int & n)
{
    damage[n][1]= RATIO_HE*1.0e-06*DPA_RATE[n]*DENSITY*VOLUME;
    totalRate[n] += damage[n][1];
}

void Damage::computeDamageTwo(const int & n)
{
    
    double concentration_H = 1.34;
    double flux_H = 4.00e+16;
    /*if (n == 0) {
        //damage[n][2] = concentration_H * flux_H * volume;
        damage[n][2] = 0.0; // no hydrogen insertion
    }
    else {
        damage[n][2] = 0.0;
    }
    */
    damage[n][2] = RATIO_H*1.0e-06*DPA_RATE[n] * DENSITY*VOLUME;
    //damage[n][2] = 0.0;
    totalRate[n] += damage[n][2];
}

void Damage::computeDamageOther(const int & n, const int & m)
{
    damage[n][m] = 0.0;
    totalRate[n] += damage[n][m];
}

double Damage::getDpaRate(const int& n){
    return DPA_RATE[n];
}
