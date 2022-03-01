! Created by jonas on 1/30/22.

module hdnnps
    use precision
    use atomicStructure
    use symmetryFunctions
    use constants
    use neighborLists
    use neuralNetworks
    use sorting
    implicit none

    type hdnnp
        type(symFunctionList) :: sfLists(maxElemNum)
        type(neuralNetwork) :: nnList(maxElemNum)
        real(dp) :: maxRc
        real(dp) :: emin, emax
        real(dp) :: atomenergies(maxElemNum)
        logical :: scalesf
        logical :: centersf
    end type hdnnp

contains

    subroutine initHdnnp(handle) !
        type(hdnnp), intent(out) :: handle
        handle%atomenergies = 0._dp
        handle%scalesf = .false.
        handle%centersf = .false.

        ! init all the weights, scaling, sfs, etc.
    end subroutine

    subroutine parseRuNNerInput(handle, path)
        type(hdnnp), intent(inout) :: handle
        logical, parameter :: debug = .false.
        character(len=*) :: path
        integer :: unit, iostat
        character(len=256) :: line
        integer :: nlayers
        integer, allocatable :: nnodes(:)
        character(len=256), allocatable :: words(:)
        integer :: nwords, i, iel
        character(len=2) :: sfEls(2), sfMainEl, elstr
        real(dp) :: rcut, eta, rs, lambda, zeta
        real(dp) :: maxRCut
        integer :: sfname, isfMainEl, isfEls(2)
        integer, parameter :: nmaxSFs = 200
        integer :: nEls
        character(len=2), allocatable :: elements(:), tmpelements(:)
        integer, allocatable :: ielements(:), index(:)
        type(symFunction), allocatable :: tmpSFlist(:)
        character(len=1), allocatable :: act(:)
        real(dp) :: atomenergy


        call initHdnnp(handle)

        maxRCut = 0._dp

        open(file=trim(path) // 'input.nn', newunit=unit, iostat=iostat, action='read')
        if (iostat /= 0) then
            print*, path // 'input.nn not found'
            stop
        end if
        read(unit, '(A)', iostat=iostat) line
        do while (iostat == 0)
            !line = adjustl(trim(line))
            call splitString(line, nwords, words)
            if (nwords > 0) then
                if (trim(words(1)) == 'number_of_elements') then
                    read(words(2), *) nEls
                    if (debug) print*, 'nEls', nEls
                    allocate(elements(nEls))
                    allocate(ielements(nEls))
                    allocate(tmpelements(nEls))
                    allocate(index(nEls))
                end if

                if (trim(words(1)) == 'elements') then
                    do i=1,nEls
                        read(words(1+i), *) elements(i)
                        allocate(handle%sfLists(elemSymToNum(elements(i)))%sfs(nmaxSFs))
                        handle%sfLists(elemSymToNum(elements(i)))%nsf = 0
                    end do
                    ielements = elemSymToNum(elements)
                    call mergesort(nEls, 1._dp * ielements, index)
                    tmpElements = elements
                    !print*, index
                    do i=1,nEls
                        elements(i) = tmpelements(index(i))
                    end do
                    if (debug) print*, 'elements  ', elements
                end if

                if (trim(words(1)) == 'atom_energy') then
                    read(words(2), *) elstr
                    read(words(3), *) atomenergy
                    if (debug) print*, 'atom_energy', elstr, atomenergy
                    handle%atomenergies(elemsymtonum(elstr)) = atomenergy
                end if

                if (trim(words(1)) == 'scale_symmetry_functions') then
                    handle%scalesf = .true.
                    if (debug) print*, 'scalesf'
                end if

                if (trim(words(1)) == 'center_symmetry_functions') then
                    handle%centersf = .true.
                    if (debug) print*, 'centersf'
                end if

                if (trim(words(1)) == 'global_hidden_layers_short') then
                    read(words(2), *) nlayers
                    if (debug) print*, 'nlayers', nlayers
                    allocate(nnodes(nlayers + 1))
                    allocate(act(nlayers + 1))
                end if

                if (trim(words(1)) == 'global_nodes_short') then
                    do i=1,nlayers
                        read(words(1+i), *) nnodes(i)
                    end do
                    nnodes(nlayers+1) = 1
                    if (debug) print*, 'nnodes', nnodes
                end if

                if (trim(words(1)) == 'global_activation_short') then
                    do i=1,nlayers+1
                        read(words(1+i), *) act(i)
                    end do
                    nnodes(nlayers+1) = 1
                    if (debug) print*, 'act  ', act
                end if

                if (trim(words(1)) == 'symfunction_short') then
                    read(words(3), *) sfname
                    if (sfname == 2) then
                        read(words(2), *) sfMainEl
                        read(words(4), *) sfEls(1)
                        read(words(5), *) eta
                        read(words(6), *) rs
                        read(words(7), *) rcut
                        isfMainEl = elemSymToNum(sfMainEl)
                        isfEls(1) = elemSymToNum(sfEls(1))
                        if (debug) print*, 'sf2', isfMainEl, isfEls(1), eta, rs, rcut

                        if (rcut > maxRCut) then
                            maxRCut = rcut
                        end if

                        i = handle%sfLists(isfMainEl)%nsf + 1
                        if (i > nmaxSFs) then
                            print*, 'Error: more symmetry functions in input file than nmaxSFs allows'
                            stop
                        end if
                        handle%sfLists(isfMainEl)%nsf = i
                        handle%sfLists(isfMainEl)%sfs(i)%name = 2
                        handle%sfLists(isfMainEl)%sfs(i)%sfEls = isfEls
                        handle%sfLists(isfMainEl)%sfs(i)%eta = eta
                        handle%sfLists(isfMainEl)%sfs(i)%rs = rs
                        handle%sfLists(isfMainEl)%sfs(i)%rcut = rcut

                    end if
                    if (sfname == 3) then
                        read(words(2), *) sfMainEl
                        read(words(4), *) sfEls(1)
                        read(words(5), *) sfEls(2)
                        read(words(6), *) eta
                        read(words(7), *) lambda
                        read(words(8), *) zeta
                        read(words(9), *) rcut
                        isfMainEl = elemSymToNum(sfMainEl)
                        isfEls = elemSymToNum(sfEls)
                        if (debug) print*, 'sf3', isfMainEl, isfEls(:), eta, lambda, zeta, rcut

                        if (rcut > maxRCut) then
                            maxRCut = rcut
                        end if

                        i = handle%sfLists(isfMainEl)%nsf + 1
                        if (i > nmaxSFs) then
                            print*, 'Error: more symmetry functions in input file than nmaxSFs allows'
                            stop
                        end if
                        handle%sfLists(isfMainEl)%nsf = i
                        handle%sfLists(isfMainEl)%sfs(i)%name = 3
                        handle%sfLists(isfMainEl)%sfs(i)%sfEls = isfEls
                        handle%sfLists(isfMainEl)%sfs(i)%eta = eta
                        handle%sfLists(isfMainEl)%sfs(i)%lambda = lambda
                        handle%sfLists(isfMainEl)%sfs(i)%zeta = zeta
                        handle%sfLists(isfMainEl)%sfs(i)%rcut = rcut

                    end if
                end if
            end if
            read(unit, '(A)', iostat=iostat) line
        end do
        close(unit)

        handle%maxRc = maxRCut

        ! resize sf lists and sorting
        do i=1,nEls
            if (allocated(tmpSFlist)) deallocate(tmpSFlist)
            allocate(tmpSFlist(handle%sfLists(elemSymToNum(elements(i)))%nsf))
            tmpSFlist(:) = handle%sfLists(elemSymToNum(elements(i)))%sfs(:handle%sfLists(elemSymToNum(elements(i)))%nsf)
            deallocate(handle%sfLists(elemSymToNum(elements(i)))%sfs)
            allocate(handle%sfLists(elemSymToNum(elements(i)))%sfs(handle%sfLists(elemSymToNum(elements(i)))%nsf))
            handle%sfLists(elemSymToNum(elements(i)))%sfs(:) = tmpSFlist(:)
            call sortSymFunctions(handle%sfLists(elemSymToNum(elements(i)))%nsf, handle%sfLists(elemSymToNum(elements(i)))%sfs)
        end do

        call readScaling(nEls, elemSymToNum(elements), path, handle%sfLists, handle%emin, handle%emax)

        do i=1,nEls
            ! init nns
            iel = elemSymToNum(elements(i))
            call nn_new(&
                handle%nnList(iel), &
                handle%sfLists(iel)%nsf, &
                nLayers + 1, &
                nNodes, &
                act)

            ! read weights
            call readWeights(handle%nnList(iel), iel, path)

        end do



    end subroutine

    subroutine sortSymFunctions(n, sfs)
        integer, intent(in) :: n
        type(symFunction), intent(inout) :: sfs(n)
        character(len=256) :: keys(n)
        integer :: i, index(n)
        type(symFunction) :: tmpsfs(n)

        do i=1,n
            call symfunctionstring(sfs(i), keys(i))
        end do

        call mergesortStr(n, keys, index)

        do i=1,n
            tmpsfs(i) = sfs(index(i))
        end do
        sfs = tmpsfs
    end subroutine

    subroutine symfunctionstring(sf, txt)
        type(symFunction), intent(in) :: sf
        character(len=256), intent(out) :: txt

        if (sf%name == 2) then
            write(txt, '(I2,3F9.5,I3)') sf%name, sf%rcut, sf%eta, sf%Rs, sf%sfEls(1)

        end if
        if (sf%name == 3) then
            ! add something to lambda since negative numbers are not sorted correctly
            write(txt, '(I2,4F9.5,3I3)') sf%name, sf%rcut, sf%eta, sf%zeta, sf%lambda + 10, sf%sfEls(1), sf%sfEls(2)
        end if

    end subroutine

    subroutine readScaling(nels, els, path, sfLists, emin, emax)
        integer, intent(in) :: nels
        character(len=*), intent(in) :: path
        integer, intent(in) :: els(nels)
        type(symFunctionList), intent(inout) :: sfLists(maxElemNum)
        real(dp), intent(out) :: emin, emax
        integer :: i, unit, iostat, isf, riel, risf
        real(dp) :: sfmin, sfmax, sfmean

        open(file=trim(path) // 'scaling.data', newunit=unit, iostat=iostat, action='read')
        if (iostat /= 0) then
            print*, trim(path) // 'scaling.data not found'
            stop
        end if
        do i=1,nels
            do isf=1,sfLists(els(i))%nsf
                read(unit, *) riel, risf, sfmin, sfmax, sfmean
                sfLists(els(i))%sfs(isf)%sfmin = sfmin
                sfLists(els(i))%sfs(isf)%sfmax = sfmax
                sfLists(els(i))%sfs(isf)%sfmean = sfmean
                if (riel /= i .or. risf /= isf) then
                    print*, 'Error: wrong index reading scaling.data'
                    stop
                end if
            end do
        end do

        read(unit, *) emin, emax
        close(unit)


    end subroutine

    subroutine readWeights(nn, iel, path)
        type(neuralNetwork), intent(inout) :: nn
        integer, intent(in) :: iel
        character(len=*), intent(in) :: path
        character(len=256) :: txt
        integer :: unit, iostat
        integer :: ilayer, iin, iout

        write(txt,'(A,I3.3,A)') 'weights.', iel, '.data'
        !print*, trim(txt)
        open(file=trim(path) // txt, newunit=unit, iostat=iostat, action='read')
        if (iostat /= 0) then
            print*, path // trim(txt) // 'not found'
            stop
        end if

        do ilayer=1, nn%nLayers
            do iin = 1,nn%layers(ilayer)%nInp
                do iout = 1,nn%layers(ilayer)%nOut
                    read(unit, *) nn%layers(ilayer)%weights(iout, iin)
                end do
            end do
            do iout = 1,nn%layers(ilayer)%nOut
                read(unit, *) nn%layers(ilayer)%bias(iout)
            end do
        end do
        close(unit)

    end subroutine

    ! todo: split this into part: for the FHMC we dont need to recalc the SF if we evaluate multiple HDNNPs
    subroutine hdnnpEnergyAndForces(handle, ats, energy, force)
        type(hdnnp), intent(inout) :: handle
        type(atStruct), intent(in) :: ats
        real(dp), intent(out) :: energy
        real(dp), intent(out) :: force(3, ats%nat)
        logical, parameter :: calcSFsTogether = .true.
        logical, parameter :: measureTiming = .false.
        type(symFunctionContainer) :: sfs(ats%nat)
        type(neighborList) :: neiLists(ats%nat)
        real(dp) :: st, et

        if (measureTiming) call cpu_time(st)
        call hdnnpCalcNeighborLists(handle, ats, neiLists)
        if (measureTiming) call cpu_time(et)
        if (measureTiming) print*, 'time neighborlist', et-st

        if (measureTiming) call cpu_time(st)
        call hdnnpCalcSymfunctions(handle, ats, neiLists, sfs, calcSfsTogether)
        if (measureTiming) call cpu_time(et)
        if (measureTiming) print*, 'time symfunction', et-st

        if (measureTiming) call cpu_time(st)
        call hdnnpCalcNNs(handle, ats, neiLists, sfs, energy, force)
        if (measureTiming) call cpu_time(et)
        if (measureTiming) print*, 'time nns', et-st


    end subroutine hdnnpEnergyAndForces

    subroutine hdnnpCalcNeighborLists(handle, ats, neiLists)
        type(hdnnp), intent(in) :: handle
        type(atStruct), intent(in) :: ats
        type(neighborList), intent(out) :: neiLists(ats%nat)
        type(atStruct) :: tmpats

        tmpats = ats
        if (tmpats%periodic) then
            call as_moveAtomsIntoCell(tmpats)
        end if
        call buildNeighborList(tmpats, handle%maxRc, neiLists)

    end subroutine

    subroutine hdnnpCalcSymfunctions(handle, ats, neiLists, sfs, calcSfsTogether)
        type(hdnnp), intent(in) :: handle
        type(atStruct), intent(in) :: ats
        type(neighborList), intent(in) :: neiLists(ats%nat)
        type(symFunctionContainer), intent(out) :: sfs(ats%nat)
        logical, intent(in), optional :: calcSfsTogether
        integer :: iat, i, nNeis, isf
        type(symFunction) :: tsf
        integer, allocatable :: neiEls(:)
        !integer :: neiEls(1024)
        logical :: calcTogether

        if (present(calcSfsTogether)) then
            calcTogether = calcSfsTogether
        else
            calcTogether = .true.
        end if

        ! allocating these inside the omp loop or even the calcAllSFs subroutine causes a segfault when ifort is used (gfortran works fine)
        do iat=1,ats%nat
            sfs(iat)%nsf = handle%sfLists(ats%el(iat))%nsf
            allocate(sfs(iat)%sfs(sfs(iat)%nsf))
            do isf=1,handle%sfLists(ats%el(iat))%nsf
              allocate(sfs(iat)%sfs(isf)%dsf(3,neiLists(iat)%nNeis))
            end do
        end do

        !$omp parallel do private(nNeis, i, neiEls, isf, tsf)
        do iat=1,ats%nat
            nNeis = neiLists(iat)%nNeis
            if (allocated(neiEls)) deallocate(neiEls)
            allocate(neiEls(nNeis))
            do i=1,nNeis
                neiEls(i) = ats%el(neiLists(iat)%iNeis(i))
            end do
            !sfs(iat)%nsf = handle%sfLists(ats%el(iat))%nsf
            !allocate(sfs(iat)%sfs(sfs(iat)%nsf))
            if (calcTogether) then
                call calcAllSymmetryFunctions(&
                        1, &
                        nNeis, &
                        neiLists(iat)%neis(:,:nNeis), &
                        neiEls(:), &
                        handle%sfLists(ats%el(iat)), &
                        sfs(iat)%sfs)
            else
                do isf=1,handle%sfLists(ats%el(iat))%nsf
                    tsf = handle%sfLists(ats%el(iat))%sfs(isf)
                    allocate(sfs(iat)%sfs(isf)%dsf(3,nNeis))
                    select case (tsf%name)
                    case (2)
                        call symmetryFunction2(1, &
                                nNeis, &
                                neiLists(iat)%neis(:,:nNeis), &
                                neiEls(:nNeis), &
                                tsf%rcut, &
                                tsf%eta, &
                                tsf%Rs, &
                                tsf%sfEls(1), &
                                sfs(iat)%sfs(isf)%sf, &
                                sfs(iat)%sfs(isf)%dsf)
                    case (3)
                        call symmetryFunction3(1, &
                                nNeis, &
                                neiLists(iat)%neis(:,:nNeis), &
                                neiEls(:nNeis), &
                                tsf%rcut, &
                                tsf%eta, &
                                tsf%lambda, &
                                tsf%zeta, &
                                tsf%sfEls(:), &
                                sfs(iat)%sfs(isf)%sf, &
                                sfs(iat)%sfs(isf)%dsf)
                    case default
                        print*, 'Unknown symmetry function', tsf%name
                        stop
                    end select
                    !print*, 'SF', iat, isf, sfs(iat)%sfs(isf)%dsf
                end do
            end if
            ! scale sfs
            do isf=1,handle%sfLists(ats%el(iat))%nsf
                !print*, 'JO', iat, isf, sfs(iat)%sfs(isf)%sf
                tsf = handle%sfLists(ats%el(iat))%sfs(isf)
                if (handle%centersf) then
                    sfs(iat)%sfs(isf)%sf = sfs(iat)%sfs(isf)%sf - tsf%sfmean
                end if
                if (handle%scalesf) then
                    sfs(iat)%sfs(isf)%sf = sfs(iat)%sfs(isf)%sf / (tsf%sfmax - tsf%sfmin)
                    sfs(iat)%sfs(isf)%dsf = sfs(iat)%sfs(isf)%dsf / (tsf%sfmax - tsf%sfmin)
                end if
            end do
        end do
        !omp end parallel do

    end subroutine hdnnpCalcSymfunctions

    subroutine hdnnpCalcNNs(handle, ats, neiLists, sfs, energy, force)
        type(hdnnp), intent(inout) :: handle
        type(atStruct), intent(in) :: ats
        type(neighborList), intent(in) :: neiLists(ats%nat)
        type(symFunctionContainer), intent(in) :: sfs(ats%nat)
        real(dp), intent(out) :: energy
        real(dp), intent(out) :: force(3, ats%nat)
        integer :: iat, isf, i, nsf
        real(dp), allocatable :: sfValues(:), dEdSf(:)
        real(dp) :: ones(1), tmpEnergy(1)

        energy = 0._dp
        force = 0._dp
        ! todo: implement lapack in nn and use batches over atoms of same element for performace
        do iat=1,ats%nat
            nsf = sfs(iat)%nsf
            if (allocated(sfValues)) deallocate(sfValues)
            allocate(sfValues(nsf))
            do isf=1,nsf
                sfValues(isf) = sfs(iat)%sfs(isf)%sf
            end do
            ! todo: use blas in nn forward and backward routine
            call nn_forward(handle%nnList(ats%el(iat)), sfValues, tmpEnergy)
            energy = energy + tmpEnergy(1) + handle%atomenergies(ats%el(iat))
            ones(1) = 1._dp
            if (allocated(dEdSf)) deallocate(dEdSf)
            allocate(dEdSf(nsf))
            ! todo: add a flag to this routine to not calculate the gradients wrt the weights
            call nn_backward(handle%nnList(ats%el(iat)), ones, dEdSF, .false.)
            ! this here is actually the slowest part
            do isf=1,nsf
                do i=1, neiLists(iat)%nneis
                    force(:,neiLists(iat)%ineis(i)) = force(:,neiLists(iat)%ineis(i)) &
                            - dEdSf(isf) * sfs(iat)%sfs(isf)%dsf(:,i)
                end do
            end do
        end do
    end subroutine


end module hdnnps
