! Created by JAF on 01.03.22.

module globalHdnnps
    use atomicStructure
    use precision
    use hdnnps

    integer :: nGlobalHdnnps
    type(hdnnp), allocatable :: globalHdnnpHandles(:)

contains

    subroutine initGlobalHdnnpsFromInputFile(fname)
        character(len=*), intent(in) :: fname
        integer :: nnunit
        character(len=256) :: txt
        integer :: i

        open(newunit=nnunit, file=fname, action='read')
        read(nnunit, *) nGlobalHdnnps

        allocate(globalHdnnpHandles(nGlobalHdnnps))
        do i=1,nGlobalHdnnps
            read(nnunit, '(A)') txt
            call parseRuNNerInput(globalHdnnpHandles(i), trim(txt))
        end do
        close(nnunit)

    end subroutine

    subroutine globalHdnnpEnergyAndForces(ats, energy, forces, dEdLat, std)
        implicit none
        type(atStruct), intent(in) :: ats
        real(dp), intent(out) :: forces(3,ats%nat), energy
        real(dp), intent(out), optional :: std, dEdLat(3,3)
        integer :: i
        type(neighborList) :: neiLists(ats%nat)
        type(symFunctionContainer) :: sfs(ats%nat)
        real(dp) :: st, et
        real(dp) :: f(3,ats%nat,nGlobalHdnnps), ee(nGlobalHdnnps), dlat(3,3,nGlobalHdnnps)

        energy = 0._dp
        forces = 0._dp

        !  call cpu_time(st)
        call hdnnpCalcNeighborLists(globalHdnnpHandles(1), ats, neiLists)
        !  call cpu_time(et)
        !  print*, 'neiList', et-st
        !  call cpu_time(st)
        call hdnnpCalcSymfunctions(globalHdnnpHandles(1), ats, neiLists, sfs, .true.)
        !  call cpu_time(et)
        !  print*, 'symfunctions', et-st
        !  call cpu_time(st)
        do i=1,nGlobalHdnnps
            call hdnnpCalcNNs(globalHdnnpHandles(i), ats, neiLists, sfs, ee(i), f(:,:,i), dlat(:,:,i))
        end do
        !  call cpu_time(et)
        !  print*, 'nns', et-st
        energy = sum(ee) / nGlobalHdnnps
        !print*, 'std', sqrt(sum((ee-energy)**2) / nhdnnps)
        forces = sum(f, 3) / nGlobalHdnnps

        if (present(dEdLat)) then
            dEdLat = sum(dlat, 3) / nGlobalHdnnps
        end if

        if (present(std)) then
            std = sqrt(sum((ee - energy)**2) / nGlobalHdnnps)
        end if

    end subroutine globalHdnnpEnergyAndForces

    ! adds a bias to the energy of the form bias*std**2
    ! keeps the MD (or MC) away from regions where training data is missing
    subroutine globalHdnnpEnergyAndForcesWithBias(ats, bias, energy, forces, dEdLat, energyBias)
        implicit none
        type(atStruct), intent(in) :: ats
        real(dp), intent(in) :: bias
        real(dp), intent(out) :: forces(3,ats%nat), energy
        real(dp), intent(out), optional :: dEdLat(3,3), energyBias
        integer :: i
        type(neighborList) :: neiLists(ats%nat)
        type(symFunctionContainer) :: sfs(ats%nat)
        real(dp) :: f(3,ats%nat,nGlobalHdnnps), ee(nGlobalHdnnps), dlat(3,3,nGlobalHdnnps)
        real(dp) :: meanEnergy, meanForces(3, ats%nat), meandEdLat(3,3), std

        energy = 0._dp
        forces = 0._dp

        call hdnnpCalcNeighborLists(globalHdnnpHandles(1), ats, neiLists)

        call hdnnpCalcSymfunctions(globalHdnnpHandles(1), ats, neiLists, sfs, .true.)

        do i=1,nGlobalHdnnps
            call hdnnpCalcNNs(globalHdnnpHandles(i), ats, neiLists, sfs, ee(i), f(:,:,i), dlat(:,:,i))
        end do

        meanenergy = sum(ee) / nGlobalHdnnps
        meanforces = sum(f, 3) / nGlobalHdnnps
        meandEdLat = sum(dlat, 3) / nGlobalHdnnps

        std = sqrt(sum((ee - meanEnergy)**2) / nGlobalHdnnps)
        energy = meanEnergy + bias * std**2
        if (present(energyBias)) then
            energyBias = bias * std**2
        end if

        forces = meanForces
        do i=1,nGlobalHdnnps
            forces = forces + 2._dp * bias / nGlobalHdnnps * (ee(i) - meanEnergy) * (f(:,:,i) - meanForces)
        end do

        if (present(dEdLat)) then
            dEdLat = meandEdLat
            do i=1, nGlobalHdnnps
                dEdlat = dEdLat + 2._dp * bias / nGlobalHdnnps * (ee(i) - meanEnergy) * (dlat(:,:,i) - meandEdLat)
            end do
        end if

    end subroutine globalHdnnpEnergyAndForcesWithBias

end module globalHdnnps