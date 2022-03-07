! Created by JAF on 01.03.22.

module globalHdnnps
    use atomicStructure
    use precision
    use hdnnps

    integer :: nGlobalHdnnps
    type(hdnnp), allocatable :: globalHdnnpHandles(:)

contains

    subroutine initGlobalHdnnpsFromInputFile()
        integer :: nnunit
        character(len=256) :: txt
        integer :: i

        open(newunit=nnunit, file='nnps.txt', action='read')
        read(nnunit, *) nGlobalHdnnps

        allocate(globalHdnnpHandles(nGlobalHdnnps))
        do i=1,nGlobalHdnnps
            print*, txt
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

        dEdLat = 0._dp
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

        dEdLat = sum(dlat, 3) / nGlobalHdnnps

        if (present(std)) then
            std = sqrt(sum((ee - energy)**2) / nGlobalHdnnps)
        end if

    end subroutine globalHdnnpEnergyAndForces

end module globalHdnnps