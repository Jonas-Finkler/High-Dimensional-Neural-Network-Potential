
module mathUtil
    use precision
    use constants
    implicit none

contains

    ! todo: it would probably make sense to move the gaussian and especially the spherical harmonics stuff into their own directories

    ! avoids over-/undeflow problems when calculating log(sum(exp(x(:))))
    function logSumExp(n, x) result(s)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(:) !implicit shape array to avoid creation of temporary array
        real(dp) :: s, m
        integer :: i

        m = maxval(x)
        s = 0._dp
        do i = 1, n
            s = s + exp(x(i) - m)
        end do
        s = m + log(s)
    end function logSumExp

    function delta(i,j) result(d)
        integer, intent(in) :: i,j
        real(dp) :: d

        if (i==j) then
            d = 1._dp
        else
            d = 0._dp
        end if

    end function delta

    subroutine multGaussians3D(n, mean, sig, x, p)
        integer, intent(in) :: n
        real(dp), intent(in) :: mean(3)
        real(dp), intent(in) :: sig(n)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: p(n)

        p(:) = 1._dp / sqrt((2._dp * PI * sig(:)**2)**3) * exp(-0.5_dp * sum((x - mean)**2) / sig(:)**2)

    end subroutine multGaussians3D

    ! this is a symmetric Gaussian for now
    function gaussian3D(mean, sig, x) result(p)
        real(dp), intent(in) :: mean(3), sig, x(3)
        real(dp) :: p

        p = 1._dp / sqrt((2._dp * PI * sig**2)**3) * exp(-0.5_dp * sum((x - mean)**2) / sig**2)
    end function gaussian3D

    subroutine dGaussian3Ddxyz(mean, sig, x, dxyz)
        real(dp), intent(in) :: mean(3), sig, x(3)
        real(dp), intent(out) :: dxyz(3)

        dxyz(:) = 1._dp / sqrt((2._dp * PI * sig**2)**3) &
            * exp(-0.5_dp * sum((x - mean)**2) / sig**2) &
            * sig**2 * (mean(:) - x(:))

    end subroutine dGaussian3Ddxyz

    function dGaussian3DdSigma(mean, sig, x) result(p)
        real(dp), intent(in) :: mean(3), sig, x(3)
        real(dp) :: p
        real(dp) :: r

        r = sum((x-mean)**2)
        p = &
                1._dp / sqrt((2 * PI * sig**2)**3)  *  r / sig**3 * exp(-0.5_dp * r/ sig**2) &
                -3._dp / sig**4 / (2._dp * PI)**1.5_dp  *  exp(-0.5_dp * r / sig**2)

    end function

    function logGaussian3D(mean, sig, x) result(p)
        real(dp), intent(in) :: mean(3), sig, x(3)
        real(dp) :: p

        p = -0.5 * 3 * log(2 * PI * sig**2) + (-0.5_dp) * sum((x - mean)**2) / sig**2
    end function logGaussian3D


    subroutine sphericalHarmonics(xyz, lmax, ylm)
        real(dp), intent(in) :: xyz(3)
        integer, intent(in) :: lmax
        real(dp), intent(out) :: ylm((lmax+1)**2)
        real(dp) :: ap(1,3)
        real(dp) :: ylmtmp(1,(lmax+1)**2)

        ap(1,:) = xyz
        call sphhar(1,ap,lmax,(lmax+1)**2,ylmtmp)
        ylm = ylmtmp(1,:)
    end subroutine


    function factorial(n) result(r)
        integer, intent(in) :: n
        integer :: r
        integer :: i

        r = 1
        do i=2,n
            r = r * i
        end do

    end function factorial

    function semifactorial(n) result(r)
        integer, intent(in) :: n
        integer :: r
        integer :: i

        r = 1
        do i=n,2,-2
            r = r * i
        end do

    end function semifactorial

end module
