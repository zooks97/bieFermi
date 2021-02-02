import numpy as np

XACC = 1e-6
JMAX = 50
NMAX = 50
FDCUT = 30
HMCUT = 10
POSHMA = -0.5634

PI = np.pi
EE = np.e
SQRT2 = np.sqrt(2)
SQRTPI = np.sqrt(np.pi)
EESH = np.sqrt(np.e) / 2
SQ2I = np.sqrt(2) / 2
PIEESQQ = np.sqrt(np.e * np.pi) / 4


def efermi(bands, weights, nelecs, smearing_width, smearing_type):
    eigmin = np.min(bands)
    eigmax = np.max(bands)

    # print("n xe1      xe2      f          fmid       f * fmid")

    xe1 = eigmin
    xe2 = eigmax
    xe0 = (xe1 + xe2) / 2
    n = 1

    f = smear(bands, weights, xe1, nelecs, smearing_width, smearing_type)
    fmid = smear(bands, weights, xe2, nelecs, smearing_width, smearing_type)

    # print(f"{n} {xe1:0.6f} {xe2:0.6f} {f:0.6f} {fmid:0.6f} {f * fmid:0.6f}")

    while f * fmid >= 0:
        if n >= NMAX:
            raise Exception(
                "Could not bracket Fermi energy. Is the electronic temperature too small?"
            )

        n += 1
        xe1 = xe0 - n * smearing_width
        xe2 = xe0 + (n - 0.5) * smearing_width

        f = smear(bands, weights, xe1, nelecs, smearing_width, smearing_type)
        fmid = smear(bands, weights, xe2, nelecs, smearing_width,
                     smearing_type)

        # print(
        #     f"{n} {xe1:0.6f} {xe2:0.6f} {f:0.6f} {fmid:0.6f} {f * fmid:0.6f}")

    if f < 0:
        rtbis = xe1
        dx = xe2 - xe1
    else:
        rtbis = xe2
        dx = xe1 - xe2

    j = 0
    while (np.abs(dx) > XACC) and (fmid != 0):
        if j >= JMAX:
            raise Exception("Reached maximum number of bisections.")

        dx *= 0.5
        xmid = rtbis + dx

        fmid = smear(bands, weights, xmid, nelecs, smearing_width,
                     smearing_type)

        if fmid <= 0:
            rtbis = xmid

        j += 1

    ef = rtbis

    return ef


def smear(bands, weights, xe, nelecs, smearing_width, smearing_type):
    z = 0
    for i, kpt in enumerate(bands):
        for j, eigval in enumerate(kpt):
            x = (xe - eigval) / smearing_width
            if smearing_type == 1:
                z += weights[i] * (2 - np.erfc(x))
            elif smearing_type == 2:
                z += weights[i] * fermid(-x)
            elif smearing_type == 3:
                z += weights[i] * delthm(x)
            elif smearing_type == 4:
                z += weights[i] * spline(-x)
            elif smearing_type == 5:
                z += weights[i] * poshm(x)
            elif smearing_type == 6:
                z += weights[i] * poshm2(x)
            else:
                raise Exception(f"Unknown smearing type `{smearing_type}`.")
    return z - nelecs


def fermid(xx):
    if xx > FDCUT:
        return 0
    elif xx < -FDCUT:
        return 2
    else:
        return 2 / (1 + np.exp(xx))


def delthm(xx):
    if xx > HMCUT:
        return 2
    elif xx < -HMCUT:
        return 0
    else:
        return (2 - np.erfc(xx)) + xx * np.exp(-xx * xx) / SQRTPI


def spline(x):
    if x > 0:
        return 2 * (EESH * np.exp(-(x + SQ2I) * -(x + SQ2I)))
    else:
        return 2 * (1.0 - EESH * np.exp(-(x + SQ2I) * -(x + SQ2I)))


def poshm(x):
    if x > HMCUT:
        return 2
    elif x < -HMCUT:
        return 0
    else:
        return (2 - np.erfc(x)) + (-2 * POSHMA * x * x + 2 * x +
                                   POSHMA) * np.exp(-x * x) / SQRTPI / 2


def poshm2(x):
    if x > HMCUT:
        return 2
    elif x < -HMCUT:
        return 0
    else:
        return (2 - np.erfc(x - 1 / SQRT2)
                ) + SQRT2 * np.exp(-x * x + SQRT2 * x - 0.5) / SQRTPI
