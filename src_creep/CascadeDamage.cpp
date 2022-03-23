#include "CascadeDamage.h"
#include<iostream>
CascadeDamage::~CascadeDamage()
{
    cleanDamage();
}

void CascadeDamage::generateNeutronDamage(const double &energy, int & ndef)
{
    const double fcli = 0.55; // Fraction of interstitials in clusters
    const double fclv = 0.25; // Fraction of vacancies in clusters.
    /* L. Malerbra, ...C. S. Bequart (2010) */
    const double fmd = 1.0; // kMC escape probabilty.
    const int types = 2; // Only vacancies and SIAs in a cascade.
    
    double nrt;
    double psia, pv;
    int n = 0;
    int i, j, k = 0;
    double eta = 0.75*(exp(-0.3*energy)+0.4); //see L. Malerba (2006)
    while (n < 1){ // Repeat until at least one stable Frenkel pair produced.
    //n = rint(1.49*pow(energy, 0.82)); //for tungsten
        n = rint(eta*0.8*energy/2.0/0.04);// see L. Malerba (2006)
        if(n<1){
            n = 1;
        }
    }
    /* Dynamically allocate first dimension of pointer array: */
    damage.resize(types, nullptr); // [0]: vacancy defects; [1]: SIA defects.
    /* Dynamically allocate second dimension of pointer array: */
    for (i = 0; i<types; i++) {
        damage[i] = new int[n];
    }
    // Initialize array:
    for (j = 0; j<types; j++)
        for (i = 0; i<n; i++)
            damage[j][i] = 0;
    /* damage[2][n] is an integer array containing the counts (from 0 to the total number of defects)
     of defects and clusters ([0]: vacancy; [1]: SIA) produced by the cascade
     cascade. For example, damage[1][3]=2 means that two I3 clusters were formed. */
    ndef = n;
    if (n > 1) {
        n *= fmd;
        psia = 1.0 - pow((1 - fcli), 1.0 / ((double)types * (double)n - 1.0));
        // Parameter for the binomial sampling of SIA-clusters.
        pv = 1.0 - pow((1 - fclv), 1.0 / ((double)types * (double)n - 1.0));
        // Parameter for the binomial sampling of V-clusters.
        
        // Randomly sample SIA-clusters from the Binomial distribution:
        int nv = 0;
        while (nv <= n) {
            k = Binomial((n - 1), pv) + 1;
            damage[0][k - 1]++;
            nv += k;
            //std::cout<<"a something"<<std::endl;
        }
        int nsia = 0;
        while (nsia <= n) {
            k = Binomial((n - 1), psia) + 1;
            damage[1][k - 1]++;
            nsia += k;
            //std::cout<<"a neutron"<<std::endl;
        }
        // Complete with remaining point defects:
        if (nv < nsia) damage[0][0] += (nsia - nv);
        if (nsia < nv) damage[1][0] += (nv - nsia);
    }
    else {
        double xi = (double)rand()/RAND_MAX;
        if (xi<fmd) {
            damage[0][0] = 1;
            damage[1][0] = 1;
        }
    }
}

void CascadeDamage::generateIonDamage(const double & energy, int & ndef)
{
    const double fcli = 0.55; // Fraction of interstitials in clusters
    const double fclv = 0.25; // Fraction of vacancies in clusters, both from: [L Malerba, JNM 351 (2006) 28].
    const double fmd = 0.65; // kMC escape probabilty
    // from [Soneda and Diaz de la Rubia, Phil Mag A 78 (1998) 995].
    const int types = 2; // Only vacancies and SIAs in a cascade.
    
    double nrt, eta;
    double psia, pv;
    int n = 0;
    int k, j, i;
    
    while (n < 1) { // Repeat until at least one stable Frenkel pair produced.
        eta = exp( -3.57 * energy ) + 0.3; // from: [L Malerba, JNM 351 (2006) 28].
        nrt = eta * 0.8 * energy / TDE / 2;
        n = rint(Poisson(nrt)); // Given the average number of defects, sample  from some statistical
        // distribution to give some variability (Poisson distribution only.
    }
    
    // Dynamically allocate first dimension of pointer array:
    damage.resize(types, nullptr); // [0]: vacancy defects; [1]: SIA defects.
    
    // Dynamically allocate second dimension of pointer array:
    for (i = 0; i<types; i++) {
        damage[i] = new int[n];
    }
    
    // This allocation overdimensions damage[][].
    
    // Initialize array:
    for (j = 0; j < types; ++j)
        for (i = 0; i<n; ++i)
            damage[j][i] = 0;
    
    /* damage[2][n] is an integer array containing the counts (from 0 to the total number of defects)
     of defects and clusters ([0]: vacancy; [1]: SIA) produced by the cascade
     cascade. For example, damage[1][3]=2 means that two I3 clusters were formed. */
    
    ndef = n;
    if (n > 1) {
        n *= fmd;
        psia = 1.0 - pow((1.0 - fcli), 1.0 / ((double)types * (double)n - 1.0));
        // Parameter for the binomial sampling of SIA-clusters.
        pv = 1.0 - pow((1.0 - fclv), 1.0 / ((double)types * (double)n - 1.0));
        // Parameter for the binomial sampling of V-clusters.
        
        // Randomly sample SIA-clusters from the Binomial distribution:
        int nv = 0;
        while (nv < n) {
            k = Binomial((n - 1), pv) + 1;
            damage[0][k - 1]++;
            nv += k;
        }
        int nsia = 0;
        while (nsia < n) {
            k = Binomial((n - 1), psia) + 1;
            damage[1][k - 1]++;
            nsia += k;
        }
        
        // Complete with remaining point defects:
        if (nv < nsia) damage[0][0] += (nsia - nv);
        else if (nsia < nv) damage[1][0] += (nv - nsia);
    }
    else {
        double xi = (double)rand()/RAND_MAX;
        if (xi < fmd) {
            damage[0][0] = 1;
            damage[1][0] = 1;
        }
    }
}

const int CascadeDamage::getDamage(const int & n, const int & m) const
{
    return damage[n][m];
}

void CascadeDamage::cleanDamage()
{
    if (!damage.empty()) {
        for (int i = 0; i < damage.size(); ++i) {
            delete[] damage[i];
        }
        damage.clear();
    }
}

int CascadeDamage::size()
{
    return damage.size();
}
