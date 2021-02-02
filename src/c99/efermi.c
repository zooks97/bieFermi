#include "efermi.h"

double efermi(size_t nkpts, size_t nbands,
              double bands[nkpts][nbands], double weights[nkpts],
              int nelecs,
              double smearing_width, int smearing_type)
// calculate Fermi energy using bisection
{
    double eigmin, eigmax, xe0, xe1, xe2, dx, xmid, f, fmid, ef;
    int i, j, k, n;

    // get eigenvalue extrema
    eigmin = bands[0][0];
    eigmax = bands[0][0];
    for (i = 0; i < nkpts; i++)
    {
        for (j = 0; j < nbands; ++j)
        {
            if (bands[i][j] < eigmin)
                eigmin = bands[i][j];
            if (bands[i][j] > eigmax)
                eigmax = bands[i][j];
        }
    }

    // set initial bounds and values
    xe1 = eigmin;
    xe2 = eigmax;
    xe0 = (xe1 + xe2) * 0.5;

    // calculate initial f, fmid
    f = smear(nkpts, nbands, bands, weights, xe1, nelecs, smearing_width, smearing_type);
    fmid = smear(nkpts, nbands, bands, weights, xe2, nelecs, smearing_width, smearing_type);

    // find xe1, xe2 which bracket the Fermi energy
    n = 1;
    while (f * fmid >= 0.0)
    {
        if (n >= NMAX)
        {
            printf("Could not bracket Fermi energy. Is the electronic temperture too small?\n");
            return 0.0;
        }

        n++;
        xe1 = xe0 - (double)n * smearing_width;
        xe2 = xe0 + ((double)n - 0.5) * smearing_width;

        f = smear(nkpts, nbands, bands, weights, xe1, nelecs, smearing_width, smearing_type);
        fmid = smear(nkpts, nbands, bands, weights, xe2, nelecs, smearing_width, smearing_type);
    }

    // set initial guess
    if (f < 0.0)
    {
        ef = xe1;
        dx = xe2 - xe1;
    }
    else
    {
        ef = xe2;
        dx = xe1 - xe2;
    }

    // do bisection
    j = 0;
    while ((fabs(dx) > XACC) && (fmid != 0.0))
    {
        if (j >= JMAX)
        {
            printf("Reached maximum number of bisections.\n");
            return 0.0;
        }
        dx *= 0.5;
        xmid = ef + dx;

        fmid = smear(nkpts, nbands, bands, weights, xmid, nelecs, smearing_width, smearing_type);
        if (fmid <= 0.0)
            ef = xmid;

        j++;
    }

    return ef;
}

double smear(size_t nkpts, size_t nbands,
             double bands[nkpts][nbands], double weights[nkpts],
             double xe,
             int nelecs, float smearing_width, int smearing_type)
// calculate smeared value used for bisection
{
    register int i, j;
    double x, z;
    SMEARING_FUNC *smearing_funcs[6] = {gaussian, fermid, delthm, spline, poshm, poshm2};

    z = 0.0;
    for (i = 0; i < nkpts; ++i)
    {
        for (j = 0; j < nbands; ++j)
        {
            x = (xe - bands[i][j]) / smearing_width;
            z += weights[i] * smearing_funcs[smearing_type - 1](x);
        }
    }
    return z - (double)nelecs;
}

double *occupations(int_t nkpts, int_t nbands,
                    double bands[nkpts][nbands], double ef,
                    double smearing_width, int smearing_type)
// construct occupations given ef, smearing_width, and smearing_type
{
    register int i, j;
    double x;
    double occ[nkpts][nbands];
    SMEARING_FUNC *smearing_funcs[6] = {gaussian, fermid, delthm, spline, poshm, poshm2};

    for (i = 0; i < nkpts; ++i)
    {
        for (j = 0; j < nbands; ++j)
        {
            x = (ef - bands[i, j]) / smearing_width;
            occ[i, j] = smearing_funcs[smearing_type](x) * 0.5;
        }
    }

    return occ;
}

// TODO: implement correction functions for each smearing type,
// TODO: and implement similarly to `smear`
// double smearing_correction(int_t nkpts, int_t nbands,
//                            double bands[nkpts][nbands], double weights[nkpts],
//                            double ef,
//                            double smearing_width, int smearing_type)
// // calculate smearing correction to retrieve "o-temperature" energy
// {
// }

double gaussian(double x)
// Gaussian
{
    return 2.0 - erfc(x);
}

double fermid(double x)
// Fermi-Dirac
{
    x = -x; // AZ: flip around x-axis so all smearing func.s called the same way
    if (x > FDCUT)
        return 0.0;
    else if (x < -FDCUT)
        return 2.0;
    else
        return 2.0 / (1.0 + exp(x));
}

double delthm(double x)
// Hermite delta expansion
{
    if (x > HMCUT)
        return 2.0;
    else if (x < -HMCUT)
        return 0.0;
    else
        return (2.0 - erfc(x)) + x * exp(-x * x) * ISQRTPI;
}

double spline(double x)
// Gaussian spline
{
    x = -x; // AZ: flip around x-axis so all smearing func.s called the same way
    if (x > 0.0)
        return 2.0 * (EESQH * exp(-(x + ISQRT2) * (x + ISQRT2)));
    else
        return 2.0 * (1.0 - EESQH * exp(-(x - ISQRT2) * (x - ISQRT2)));
}

double poshm(double x)
// Positive Hermite (cold I)
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

double poshm2(double x)
// Positive Hermite (cold II)
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