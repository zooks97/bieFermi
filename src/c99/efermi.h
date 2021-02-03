#ifndef EFERMI_H_
#define EFERMI_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// parameters
#define XACC 1e-6      // desired accuracy on the Fermi energy
#define JMAX 10000     // maximum number of bisections
#define NMAX 100000    // maximum number of xe1, xe2 redefinitions
#define FDCUT 30.0     // Fermi-Dirac cutoff
#define HMCUT 10.0     // Hermite cutoff
#define POSHMA -0.5634 // Positive Hermite (cold I) `a`

// math constants
#define SQRT2 1.414213562373095     // sqrt(2)
#define ISQRT2 0.707106781186548    // sqrt(2) / 2,  1 / sqrt(2)
#define ISQRTPI 0.564189583547756   // 1 / sqrt(pi)
#define I2SQRTPI 0.2820947917738781 // 1 / (2 * sqrt(pi))
#define HSQRTE 0.824360635350064    // sqrt(e) / 2

// smearing function pointer type
typedef double SMEARING_FUNC(double);

// bisection
double efermi(size_t nkpts, size_t nbands,
              double bands[nkpts][nbands], double weights[nkpts],
              int nelecs, double smearing_width, int smearing_type); // Fermi energy
double smear(size_t nkpts, size_t nbands,
             double bands[nkpts][nbands], double weights[nkpts],
             double x, int nelecs, float smearing_width, int smearing_type); // calculate smeared band

// smearing
double gaussian(double x); // Gaussian
double fermid(double x);   // Fermi-Dirac
double delthm(double x);   // Hermite delta
double spline(double x);   // Gaussian spline
double poshm(double x);    // Positive Hermite (cold I)
double poshm2(double x);   // Positive Hermite (cold II)

#endif // EFERMI_H_