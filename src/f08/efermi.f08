module efermif08
    implicit none
    private

    double precision, parameter :: XACC = 1.d-6
    integer, parameter :: NMAX = 100000, JMAX = 10000
    double precision, parameter :: FDCUT = 30.0, HMCUT = 10.0
    double precision, parameter :: POSHMA = -0.5634
    double precision, parameter :: SQRT2 = 1.414213562373095
    double precision, parameter :: ISQRT2 = 0.707106781186548
    double precision, parameter :: ISQRTPI = 0.564189583547756   
    double precision, parameter :: I2SQRTPI = 0.2820947917738781
    double precision, parameter :: HSQRTE = 0.824360635350064   

    public :: efermi, smear, gauss, fermid, delthm, spline, poshm, poshm2

contains

function efermi(nkpt, nbnd, bands, weights, nelec, swidth, stype) result(rtb)
    implicit none
    integer, intent(in) :: nkpt, nbnd, nelec, stype
    double precision, intent(in) :: bands(nkpt, nbnd), weights(nkpt), swidth
    double precision :: x0, x1, x2, dx, xmid, f, fmid, rtb
    integer j, n

    ! Get min, max eigenvalue and set as initial bounds
    x1 = minval(bands)
    x2 = maxval(bands)
    x0 = (x1 + x2) * 0.d5

    ! Calculate initial f, fmid
    f = smear(nkpt, nbnd, bands, weights, x1, nelec, swidth, stype)
    fmid = smear(nkpt, nbnd, bands, weights, x2, nelec, swidth, stype)

    ! Find bounds which bracket the Fermi energy
    do n = 1, NMAX
        if (f * fmid >= 0.d0) then
            x1 = x0 - dble(n) * swidth
            x2 = x0 + (dble(n) - 0.d5) * swidth
            f = smear(nkpt, nbnd, bands, weights, x1, nelec, swidth, stype)
            fmid = smear(nkpt, nbnd, bands, weights, x2, nelec, swidth, stype)
        else
            continue
        end if
    end do
    if (f * fmid >= 0.d0) then
        write (*,*) "Could not bracket the Fermi energy. Smearing too small?"
        rtb = 0.d0
        return
    end if

    ! Set initial Fermi energy guess
    if (f < 0.0) then
        dx = x2 - x1
        rtb = x1
    else
        dx = x1 - x2
        rtb = x2
    end if

    ! Do bisection
    do j = 1, JMAX
        if ((abs(dx) <= XACC) .or. (fmid == 0.d0)) then
            return
        end if
        dx = dx * 0.d5
        xmid = rtb + dx
        fmid = smear(nkpt, nbnd, bands, weights, xmid, nelec, swidth, stype)
        if (fmid <= 0.d0) then
            rtb = xmid
        end if
    end do
    write (*,*) "Reached maximum number of bisections"
    return
end function efermi

function smear(nkpt, nbnd, bands, weights, xe, nelec, swidth, stype) result(z)
    implicit none
    integer, intent(in) :: nkpt, nbnd, nelec, stype
    double precision, intent(in) :: bands(nkpt, nbnd), weights(nkpt), swidth, xe
    double precision :: x, z
    integer i, j

    ! Set up an interface containing `func`, which has the form of a
    ! smearing function
    interface
    function func (z)
        double precision :: func
        double precision, intent (in) :: z
    end function func
    end interface
    ! Create a pointer of type `func`
    procedure (func), pointer :: sfunc => null ()
    ! Select the desired smearing function
    select case (stype)
        case (1)
            sfunc => gauss
        case (2)
            sfunc => fermid
        case (3)
            sfunc => delthm
        case (4)
            sfunc => spline
        case (5)
            sfunc => poshm
        case (6)
            sfunc => poshm2
    end select

    z = 0.0
    do j = 1, nbnd
        do i = 1, nkpt
            x = (xe - bands(i, j)) / swidth
            z = z + weights(i) * sfunc(x)
        end do
    end do
    z = z - dble(nelec)
    return
end function smear

pure function gauss(x) result(y)
    implicit none
    double precision, intent(in) :: x
    double precision :: y
    y = 2.0 - erfc(x)
    return
end function gauss

pure function fermid(x) result(y)
    implicit none
    double precision, intent(in) :: x
    double precision :: y
    if (-x > FDCUT) then
        y = 0.0
    else if (-x < -FDCUT) then
        y = 2.0
    else
        y = 2.0 / (1.0 + exp(-x))
    end if
    return
end function fermid

pure function delthm(x) result(y)
    implicit none
    double precision, intent(in) :: x
    double precision :: y
    if (x > HMCUT) then
        y = 2.0
    else if (x < -HMCUT) then
        y = 0.0
    else
        y = (2.0 - erfc(x)) + x * exp(-x * x) * ISQRTPI
    end if
    return
end function delthm

pure function spline(x) result(y)
    implicit none
    double precision, intent(in) :: x
    double precision :: y
    if (-x > 0.0) then
        y = 2.0 * (HSQRTE * exp(-(x + ISQRT2)**2))
    else
        y = 2.0 * (1.d0 - HSQRTE * exp(-(x - ISQRT2)**2))
    end if
    return
end function spline

pure function poshm(x) result(y)
    implicit none
    double precision, intent(in) :: x
    double precision :: y
    if (x > HMCUT) then
        y = 2.0
    else if (x < -HMCUT) then
        y = 0.0
    else
        y = (2.d0 - erfc(x)) + (-2.0 * POSHMA * x**2 + 2.0 * x + POSHMA) * exp(-x**2) * I2SQRTPI
    end if
    return
end function poshm

pure function poshm2(x) result(y)
    implicit none
    double precision, intent(in) :: x
    double precision :: y
    if (x > HMCUT) then
        y = 2.0
    else if (x < -HMCUT) then
        y = 0.0
    else
        y = (2.0 - erfc(x - ISQRT2)) + SQRT2 * exp(-x**2 + SQRT2 * x - 0.5) * ISQRTPI
    end if
    return
end function poshm2

end module efermif08