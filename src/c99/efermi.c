#include "efermi.h"

double efermi(size_t nkpt, size_t nbnd,
              double bands[nkpt][nbnd], double weights[nkpt],
              int nelec,
              double swidth, int stype)
// calculate Fermi energy using bisection
{
    double eigmin, eigmax, x0, x1, x2, dx, xmid, f, fmid, rtb;
    int i, j, k, n;

    // Get min, max eigenvalue and set as initial bounds
    x1 = bands[0][0];
    x2 = bands[0][0];
    for (register int i = 0; i < nkpt; i++)
    {
        for (register int j = 0; j < nbnd; ++j)
        {
            if (bands[i][j] < x1)
                x1 = bands[i][j];
            if (bands[i][j] > x2)
                x2 = bands[i][j];
        }
    }
    x0 = (x1 + x2) * 0.5;

    // Calculate initial f, fmid
    f = smear(nkpt, nbnd, bands, weights, x1, nelec, swidth, stype);
    fmid = smear(nkpt, nbnd, bands, weights, x2, nelec, swidth, stype);

    // Find bounds which bracket the Fermi energy
    for (register int n; n < NMAX; ++n)
    {
        if (f * fmid >= 0.0)
        {
            x1 = x0 - (double)n * swidth;
            x2 = x0 + ((double)n - 0.5) * swidth;
            f = smear(nkpt, nbnd, bands, weights, x1, nelec, swidth, stype);
            fmid = smear(nkpt, nbnd, bands, weights, x2, nelec, swidth, stype);
        }
        else
        {
            break;
        }
    }
    if (f * fmid >= 0.0)
    {
        printf("Could not bracket Fermi energy. Smearing too small?\n");
        return 0.0;
    }

    // Set initial Fermi energy guess
    rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);

    // Do bisection
    for (register int j = 0; j < JMAX; ++j)
    {
        if ((fabs(dx) <= XACC) || (fmid == 0.0))
        {
            return rtb;
        }
        dx = dx * 0.5;
        xmid = rtb + dx;
        fmid = smear(nkpt, nbnd, bands, weights, xmid, nelec, swidth, stype);
        if (fmid <= 0.0)
        {
            rtb = xmid;
        }
    }
    printf("Reached maximum number of bisections.\n");
    return rtb;
}

double smear(size_t nkpt, size_t nbnd,
             double bands[nkpt][nbnd], double weights[nkpt],
             double xe,
             int nelec, float swidth, int stype)
// calculate smeared value used for bisection
{
    double x, z;
    SMEARING_FUNC *smearing_funcs[6] = {gaussian, fermid, delthm, spline, poshm, poshm2};

    z = 0.0;
    for (register int i = 0; i < nkpt; ++i)
    {
        for (register int j = 0; j < nbnd; ++j)
        {
            x = (xe - bands[i][j]) / swidth;
            z = z + weights[i] * smearing_funcs[stype - 1](x);
        }
    }
    return z - (double)nelec;
}

// Gaussian
static inline double gaussian(double x)
{
    return 2.0 - erfc(x);
}

// Fermi-Dirac
static inline double fermid(double x)
{
    // AZ: flip around x-axis so all smearing func.s called the same way
    x = -x;
    if (x > FDCUT)
        return 0.0;
    else if (x < -FDCUT)
        return 2.0;
    else
        return 2.0 / (1.0 + exp(x));
}

// Hermite delta expansion
static inline double delthm(double x)
{
    if (x > HMCUT)
        return 2.0;
    else if (x < -HMCUT)
        return 0.0;
    else
        return (2.0 - erfc(x)) + x * exp(-x * x) * ISQRTPI;
}

// Gaussian spline
static inline double spline(double x)
{
    // AZ: flip around x-axis so all smearing func.s called the same way
    x = -x;
    if (x > 0.0)
        return 2.0 * (HSQRTE * exp(-(x + ISQRT2) * (x + ISQRT2)));
    else
        return 2.0 * (1.0 - HSQRTE * exp(-(x - ISQRT2) * (x - ISQRT2)));
}

// Positive Hermite (cold I)
static inline double poshm(double x)
{
    // NM: NOTE: g's are all intended to be normalized to 1!
    // NM: function = 2 * int_minf^x [g(t)] dt
    if (x > HMCUT)
        return 2.0;
    else if (x < -HMCUT)
        return 0.0;
    else
        return (2.0 - erfc(x)) + (-2.0 * POSHMA * x * x + 2.0 * x + POSHMA) * exp(-x * x) * I2SQRTPI;
}

// Positive Hermite (cold II)
static inline double poshm2(double x)
{
    // NM: NOTE: g's are all intended to be normalized to 1!
    // NM: function = 2 * int_minf^x [g(t)] dt
    if (x > HMCUT)
        return 2.0;
    else if (x < -HMCUT)
        return 0.0;
    else
        return (2.0 - erfc(x - ISQRT2)) + SQRT2 * exp(-x * x + SQRT2 * x - 0.5) * ISQRTPI;
}