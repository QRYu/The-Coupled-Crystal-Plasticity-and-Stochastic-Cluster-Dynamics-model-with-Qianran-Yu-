// CascadeDamage.h -- class CascadeDamage
#pragma once
#include<vector>
#include<cstdlib>
#include<cmath>
#include"rvgs.h"
#include"constant.h"

using namespace std;
class CascadeDamage {
private:
    vector<int*> damage;
public:
    CascadeDamage() {};/* constructor */
    ~CascadeDamage(); /* destructor */
    /* generate cascade damage functions */
    void generateNeutronDamage(const double&, int&);
    void generateIonDamage(const double&, int&);
    const int getDamage(const int&, const int&) const;
    void cleanDamage();
    int size();
};

