// rvgs.cpp -- cpp file of rvgs.c
#include<cmath>
#include<ctime>
#include<cstdlib>
#include<cstdio>
#include"rvgs.h"


int Bernoulli(double p)
/* ========================================================
 * Returns 1 with probability p or 0 with probability 1 - p.
 * NOTE: use 0.0 < p < 1.0
 * ========================================================
 */
{
    return (((double)rand() / RAND_MAX < (1 - p)) ? 0 : 1);
}

int Binomial(int n, double p)
/* ================================================================
 * Returns a binomial distributed integer between 0 and n inclusive.
 * NOTE: use n > 0 and 0.0 < p < 1.0
 * ================================================================
 */
{
    int i, x = 0;
    
    for (i = 0; i < n; i++)
        x += Bernoulli(p);
    return (x);
}

int Equilikely(int a, int b)
/*===================================================================
 * Returns an equilikely distributed integer between a and b inclusive.
 * NOTE: use a < b
 * ===================================================================
 */
{
    return (a + (int)((b - a + 1) * (double)rand() / RAND_MAX));
}

int Geometric(double p)
/* ====================================================
 * Returns a geometric distributed non-negative integer.
 * NOTE: use 0.0 < p < 1.0
 * ====================================================
 */
{
    return ((int)(log(1.0 - (double)rand() / RAND_MAX) / log(p)));
}

int Pascal(int n, double p)
/* =================================================
 * Returns a Pascal distributed non-negative integer.
 * NOTE: use n > 0 and 0.0 < p < 1.0
 * =================================================
 */
{
    int i, x = 0;
    
    for (i = 0; i < n; i++)
        x += Geometric(p);
    return (x);
}

int Poisson(double m)
/* ==================================================
 * Returns a Poisson distributed non-negative integer.
 * NOTE: use m > 0
 * ==================================================
 */
{
    double t = 0.0;
    int   x = 0;
    
    while (t < m) {
        t += Exponential(1.0);
        x++;
    }
    return (x - 1);
}

double Uniform(double a, double b)
/* ===========================================================
 * Returns a uniformly distributed real number between a and b.
 * NOTE: use a < b
 * ===========================================================
 */
{
    return (a + (b - a) * (double)rand() / RAND_MAX);
}

double Exponential(double m)
/* =========================================================
 * Returns an exponentially distributed positive real number.
 * NOTE: use m > 0.0
 * =========================================================
 */
{
    double k = -m * log(1.0 - (double)rand() / RAND_MAX);
    return k;
}

double Erlang(int n, double b)
/* ==================================================
 * Returns an Erlang distributed positive real number.
 * NOTE: use n > 0 and b > 0.0
 * ==================================================
 */
{
    int   i;
    double x = 0.0;
    
    for (i = 0; i < n; i++)
        x += Exponential(b);
    return (x);
}

double Normal(double m, double s)
/* ========================================================================
 * Returns a normal (Gaussian) distributed real number.
 * NOTE: use s > 0.0
 *
 * Uses a very accurate approximation of the normal idf due to Odeh & Evans,
 * J. Applied Statistics, 1974, vol 23, pp 96-97.
 * ========================================================================
 */
{
    const double p0 = 0.322232431088;     const double q0 = 0.099348462606;
    const double p1 = 1.0;                const double q1 = 0.588581570495;
    const double p2 = 0.342242088547;     const double q2 = 0.531103462366;
    const double p3 = 0.204231210245e-1;  const double q3 = 0.103537752850;
    const double p4 = 0.453642210148e-4;  const double q4 = 0.385607006340e-2;
    double u, t, p, q, z;
    
    u = (double)rand()/RAND_MAX;
    if (u < 0.5)
        t = sqrt(-2.0 * log(u));
    else
        t = sqrt(-2.0 * log(1.0 - u));
    p = p0 + t * (p1 + t * (p2 + t * (p3 + t * p4)));
    q = q0 + t * (q1 + t * (q2 + t * (q3 + t * q4)));
    if (u < 0.5)
        z = (p / q) - t;
    else
        z = t - (p / q);
    return (m + s * z);
}

double Lognormal(double a, double b)
/* ====================================================
 * Returns a lognormal distributed positive real number.
 * NOTE: use b > 0.0
 * ====================================================
 */
{
    return (exp(a + b * Normal(0.0, 1.0)));
}

double Chisquare(int n)
/* =====================================================
 * Returns a chi-square distributed positive real number.
 * NOTE: use n > 0
 * =====================================================
 */
{
    int   i;
    double z, x = 0.0;
    
    for (i = 0; i < n; i++) {
        z = Normal(0.0, 1.0);
        x += z * z;
    }
    return (x);
}

double Student(int n)
/* ===========================================
 * Returns a student-t distributed real number.
 * NOTE: use n > 0
 * ===========================================
 */
{
    return (Normal(0.0, 1.0) / sqrt(Chisquare(n) / n));
}
