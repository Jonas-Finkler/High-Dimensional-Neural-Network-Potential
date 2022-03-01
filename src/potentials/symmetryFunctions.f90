! Created by JAF on 1/30/22.

module symmetryFunctions
    use precision
    use atomicStructure
    implicit none
    type symFunction
        integer :: name
        integer :: sfEls(2) ! only use first element if name = 2
        real(dp) :: rcut
        real(dp) :: eta
        real(dp) :: lambda
        real(dp) :: zeta
        real(dp) :: Rs
        real(dp) :: sfmin, sfmax, sfmean
    end type symFunction

    type symFunctionList
        integer :: nsf
        type(symFunction), allocatable :: sfs(:)
    end type symFunctionList

    type symFunctionValue
        integer :: nNeis
        real(dp) :: sf
        real(dp), allocatable :: dsf(:,:)
    end type symFunctionValue

    type symFunctionContainer
        integer :: nsf
        type(symFunctionValue), allocatable :: sfs(:)
    end type symFunctionContainer

contains

    ! ats should just contain the neighbours. (for performace, result is still fine if all atoms are included)
    ! keep in mind that one atom can be multiple neighbours in periodic cases or even neighbour itself
    ! keeping track of this however is not the purpose of this routine
    ! It therefore only works for free b.c. or periodic systems, where neighbouring cells are added properly
    ! This is function G2 from Behlers 2011 paper
    subroutine symmetryFunction2(iat, nat, ats, els, rcut, eta, Rs, sfEl, sf, dsf)
        integer, intent(in) :: iat
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3, nat)
        integer, intent(in) :: els(nat)
        real(dp), intent(in) :: rcut
        real(dp), intent(in) :: eta
        real(dp), intent(in) :: Rs
        integer, intent(in) :: sfEl
        real(dp), intent(out) :: sf
        real(dp), intent(out) :: dsf(3, nat)
        integer :: i
        real(dp) :: r(3), d, rdfc(3), tmp1, dfc, fc

        sf = 0._dp
        dsf = 0._dp

        do i=1,nat
            if (i/= iat) then
                if (els(i) == sfEl) then
                    r = ats(:, iat) - ats(:, i)
                    d = sqrt(sum(r**2))
                    call cutoff(d, rcut, fc, dfc)
                    rdfc = dfc * r
                    tmp1 = exp(-eta * (d-Rs)**2)
                    dsf(:, i) = -1._dp * tmp1 * rdfc + 2 * eta * (d-Rs) * r / d * tmp1 * fc
                    dsf(:, iat) = dsf(:,iat) - dsf(:,i)
                    sf = sf + tmp1 * fc
                end if
            end if
        end do


    end subroutine symmetryFunction2

    ! This is function G4 from Behlers 2011 paper
    subroutine symmetryFunction3(iat, nat, ats, els, rcut, eta, lambda, zeta, sfEl, sf, dsf)
        integer, intent(in) :: iat
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3,nat)
        integer, intent(in) :: els(nat)
        real(dp), intent(in) :: rcut
        real(dp), intent(in) :: eta
        real(dp), intent(in) :: lambda
        real(dp), intent(in) :: zeta
        integer, intent(in) :: sfEl(2)
        real(dp), intent(out) :: sf
        real(dp), intent(out) :: dsf(3, nat)
        integer :: j, k
        real(dp) :: dij(3), dik(3), djk(3), rij, rik, rjk
        real(dp) :: fcij, dfcij, fcik, dfcik, fcjk, dfcjk
        real(dp) :: costheta, dcosthetaj(3), dcosthetak(3)
        real(dp) :: fcpart, dfcpartj(3), dfcpartk(3), prefactor, cospart, dcospartj(3), dcospartk(3)
        real(dp) :: exppart, dexppartj(3), dexppartk(3)
        real(dp) :: dsfj(3), dsfk(3)


        sf = 0._dp
        dsf(:,:) = 0._dp

        do j=1,nat
            if (j == iat .or. (els(j) /= sfEl(1) .and. els(j) /= sfEl(2))) then
                cycle
            end if
            dij = ats(:,iat) - ats(:,j)
            rij = sqrt(sum(dij**2))
            call cutoff(rij, rcut, fcij, dfcij)
            do k=j+1,nat
                if (k == iat) then
                    cycle
                end if
                if ((els(j) == sfEl(1) .and. els(k) == sfEl(2)) .or. (els(j) == sfEl(2) .and. els(k) == sfEl(1))) then
                    dik = ats(:,iat) - ats(:,k)
                    rik = sqrt(sum(dik**2))
                    call cutoff(rik, rcut, fcik, dfcik) ! is recomputed a lot, could save all cutoffs between central atom
                    djk = ats(:,j) - ats(:,k)
                    rjk = sqrt(sum(djk**2))
                    call cutoff(rjk, rcut, fcjk, dfcjk)
                    costheta = (rjk**2 - rij**2 - rik**2) / (-2._dp * rij * rik) ! can save this
                    dcosthetaj(:) = -1._dp * dik(:) / (rij * rik) + dij(:) * sum(dij * dik) / (rij**3 * rik)
                    dcosthetak(:) = -1._dp * dij(:) / (rij * rik) + dik(:) * sum(dij * dik) / (rik**3 * rij)
                    fcpart = fcij * fcik * fcjk ! can save this
                    dfcpartj(:) =           fcik * (dfcjk * djk(:) * fcij - dfcij * dij(:) * fcjk)
                    dfcpartk(:) =  -1._dp * fcij * (dfcjk * djk(:) * fcik + dfcik * dik(:) * fcjk)
                    prefactor = 2._dp**(1._dp - zeta)
                    cospart = (1._dp + lambda * costheta)**zeta
                    dcospartj(:) = zeta * (1._dp + lambda * costheta)**(zeta - 1._dp) * lambda * dcosthetaj(:)
                    dcospartk(:) = zeta * (1._dp + lambda * costheta)**(zeta - 1._dp) * lambda * dcosthetak(:)
                    ! todo: might be worth retraining an nnp with eta=0 to save these terms
                    exppart =  exp(-eta * (rij**2 + rik**2 + rjk**2)) ! can save this (eta is almost always the same)
                    dexppartj(:) = exppart * eta * 2._dp * (dij(:) - djk(:))
                    dexppartk(:) = exppart * eta * 2._dp * (dik(:) + djk(:))
                    dsfj(:) = prefactor * (&
                            cospart * exppart * dfcpartj(:) + &
                            cospart * dexppartj(:) * fcpart + &
                            dcospartj(:) * exppart * fcpart)
                    dsfk(:) = prefactor * (&
                            cospart * exppart * dfcpartk(:) + &
                            cospart * dexppartk(:) * fcpart + &
                            dcospartk(:) * exppart * fcpart)
                    dsf(:,j) = dsf(:,j) + dsfj(:)
                    dsf(:,k) = dsf(:,k) + dsfk(:)
                    dsf(:,iat) = dsf(:,iat) - (dsfj(:) + dsfk(:))
                    sf = sf + prefactor * cospart * exppart * fcpart

                end if
            end do
        end do



    end subroutine symmetryFunction3

    ! assumes same cutoff for al sfs
    ! is about three times faster than the routines above
    ! sfs have to be sorted by type (name)
    ! sorting of sfs by eta improves performace
    subroutine calcAllSymmetryFunctions(iat, nat, ats, els, sfList, sfs)
        integer, intent(in) :: iat
        integer, intent(in) :: nat
        real(dp), intent(in) :: ats(3, nat)
        integer, intent(in) :: els(nat)
        type(symFunctionList), intent(in) :: sfList ! sorted by type
        type(symFunctionValue), intent(inout) :: sfs(sfList%nsf) ! they have to be already allocated. This is to avoid a bug with ifort and openmp
        integer :: sf3start
        integer :: isf, j, k
        real(dp) :: dij(3), dik(3), djk(3), rij, rik, rjk
        real(dp) :: fcij, dfcij, fcik, dfcik, fcjk, dfcjk
        real(dp) :: costheta, dcosthetaj(3), dcosthetak(3)
        real(dp) :: fcpart, dfcpartj(3), dfcpartk(3), prefactor, cospart, dcospartj(3), dcospartk(3)
        real(dp) :: exppart, dexppartj(3), dexppartk(3)
        real(dp) :: dsfj(3), dsfk(3)
        real(dp) :: tmp1, tmp2, tmp3
        type(symFunction) :: sf
        real(dp) :: rcut
        real(dp) :: lasteta
        integer :: intzeta
        real(dp) :: fcijs(nat), dfcijs(nat), dijs(3,nat), rijs(nat)

        real(dp) :: t, tt

        sf3start = -1
        rcut = sfList%sfs(1)%rcut

        do isf=1,sfList%nsf
            if (sfList%sfs(isf)%name == 3 .and. sf3start < 0) then
                sf3start = isf
            else
                if (sfList%sfs(isf)%name == 2 .and. sf3start > 0) then
                    print*, 'Error: SFs are not sorted.'
                    stop
                end if
            end if
            !allocate(sfs(isf)%dsf(3,nat))
            sfs(isf)%sf = 0._dp
            sfs(isf)%dsf(:,:) = 0._dp
            if ((sfList%sfs(isf)%rcut - rcut)**2 > 1.e-10_dp) then
                print*, 'Error: expecting same rc for all symfunctions in calcAllSymmetryFunctions'
                stop
            end if
        end do

        do j=1,nat
            if (j == iat) then
                cycle
            end if
            dijs(:,j) = ats(:,iat) - ats(:,j)
            rijs(j) = sqrt(sum(dijs(:,j)**2))
            call cutoff(rijs(j), rcut, fcijs(j), dfcijs(j))
        end do

        do j=1,nat
            if (j == iat) then
                cycle
            end if
            dij = dijs(:, j) !ats(:,iat) - ats(:,j)
            rij = rijs(j) !sqrt(sum(dij**2))
            fcij = fcijs(j)
            dfcij = dfcijs(j)
            !call cutoff(rij, rcut, fcij, dfcij)
            ! 2 body SF
            do isf=1,sf3start-1
                sf = sfList%sfs(isf)
                if (sf%sfEls(1) == els(j)) then
                    exppart = exp(-sf%eta * (rij-sf%Rs)**2)
                    sfs(isf)%dsf(:, j) = -1._dp * exppart * dfcij * dij &
                            + 2._dp * sf%eta * (rij-sf%Rs) * dij / rij * exppart * fcij
                    sfs(isf)%dsf(:, iat) = sfs(isf)%dsf(:,iat) - sfs(isf)%dsf(:,j)
                    sfs(isf)%sf = sfs(isf)%sf + exppart * fcij
                    !print*, 'X2', sfs(isf)%sf
                end if
            end do
            do k=j+1,nat
                if (k == iat) then
                    cycle
                end if
                !dik = ats(:,iat) - ats(:,k)
                !rik = sqrt(sum(dik**2))
                !call cutoff(rik, rcut, fcik, dfcik) ! is recomputed a lot, could save all cutoffs between central atom
                dik = dijs(:, k)
                rik = rijs(k)
                fcik = fcijs(k)
                dfcik = dfcijs(k)

                djk = ats(:,j) - ats(:,k)
                rjk = sqrt(sum(djk**2))
                call cutoff(rjk, rcut, fcjk, dfcjk)
                tmp1 = -1._dp / (rij * rik)
                costheta = 0.5_dp * (rjk**2 - rij**2 - rik**2) * tmp1
                tmp2 = sum(dij * dik)
                dcosthetaj(:) = tmp1 * dik(:) + dij(:) * tmp2 / (rij**3 * rik)
                dcosthetak(:) = tmp1 * dij(:) + dik(:) * tmp2 / (rik**3 * rij)
                fcpart = fcij * fcik * fcjk ! can save this
                dfcpartj(:) =           fcik * (dfcjk * djk(:) * fcij - dfcij * dij(:) * fcjk)
                dfcpartk(:) =  -1._dp * fcij * (dfcjk * djk(:) * fcik + dfcik * dik(:) * fcjk)
                ! trivial case is precomputed
                lasteta = 0._dp
                exppart = 1._dp
                dexppartj = 0._dp
                dexppartk = 0._dp
                ! 3 body sfs
                do isf=sf3start,sfList%nsf
                    sf = sfList%sfs(isf)
                    if ((els(j) == sf%sfEls(1) .and. els(k) == sf%sfEls(2)) &
                            .or. (els(j) == sf%sfEls(2) .and. els(k) == sf%sfEls(1))) then
                        ! exponentiation with reals is expensive
                        intzeta = int(sf%zeta)
                        tmp1 = 1._dp + sf%lambda * costheta
                        if (abs(sf%zeta - intzeta) < 1.e-10_dp) then
                            prefactor = 2._dp**(1 - intzeta)
                            tmp2 = tmp1**(intzeta - 1)
                        else
                            prefactor = 2._dp**(1._dp - sf%zeta)
                            tmp2 = tmp1**(sf%zeta - 1._dp)
                        end if
                        cospart = tmp2 * tmp1
                        tmp3 = tmp2 * sf%lambda * sf%zeta !/ (1._dp + sf%lambda * costheta)
                        dcospartj(:) = tmp3 * dcosthetaj(:)
                        dcospartk(:) = tmp3 * dcosthetak(:)
                        if (abs(sf%eta - lasteta) > 1.e-10_dp) then
                            exppart =  exp(-sf%eta * (rij**2 + rik**2 + rjk**2)) ! can save this (eta is almost always the same)
                            dexppartj(:) = exppart * sf%eta * 2._dp * (dij(:) - djk(:))
                            dexppartk(:) = exppart * sf%eta * 2._dp * (dik(:) + djk(:))
                            lasteta = sf%eta
                        end if
                        ! this does not help... why?
                        !tmp1 = cospart * exppart
                        !tmp2 = cospart * fcpart
                        !tmp3 = exppart * fcpart
                        dsfj(:) = prefactor * (&
                                cospart * exppart * dfcpartj(:) + &
                                cospart * dexppartj(:) * fcpart + &
                                tmp2 * sf%lambda * sf%zeta * dcosthetaj(:) * exppart * fcpart)
                                !dcospartj(:) * exppart * fcpart)
                        dsfk(:) = prefactor * (&
                                cospart * exppart * dfcpartk(:) + &
                                cospart * dexppartk(:) * fcpart + &
                                tmp2 * sf%lambda * sf%zeta * dcosthetak(:) * exppart * fcpart)
                                !dcospartk(:) * exppart * fcpart)
                        sfs(isf)%dsf(:,j) = sfs(isf)%dsf(:,j) + dsfj(:)
                        sfs(isf)%dsf(:,k) = sfs(isf)%dsf(:,k) + dsfk(:)
                        sfs(isf)%dsf(:,iat) = sfs(isf)%dsf(:,iat) - (dsfj(:) + dsfk(:))
                        sfs(isf)%sf = sfs(isf)%sf + prefactor * cospart * exppart * fcpart
                        !print*, 'X3', sfs(isf)%dsf(:, iat)
                    end if
                end do
            end do
        end do

    end subroutine

    pure subroutine cutoff(r, rcut, fc, dfc)
        real(dp), intent(in) :: r
        real(dp), intent(in) :: rcut
        real(dp), intent(out) :: fc
        real(dp), intent(out) :: dfc
        real(dp) :: t, tt

        if (r < rcut) then
            t = tanh(1._dp - r / rcut)
            tt = t**2
            fc = tt * t
            dfc = -3._dp / (rcut * r) * (tt - tt**2)
        else
            fc = 0._dp
            dfc = 0._dp
        end if

    end subroutine cutoff

end module symmetryFunctions
