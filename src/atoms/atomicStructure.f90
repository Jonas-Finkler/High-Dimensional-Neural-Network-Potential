! Created by Jonas Finkler on 5/18/20.

! convention: start subroutines acting on atomic structures with 'as_'
! this is done because the many subroutines have to be rewritten to use the new datatype
! backward compatibility has to be kept
! todo: maybe rethink or refactor this in the future
module atomicStructure
    use precision
    use constants
    use util
    use pointsets
    use dynamicArrays

    implicit none

    ! type to store atomic structures
    type atStruct
        integer :: nat
        logical :: periodic = .false.
        real(dp)  :: lat(3,3)
        real(dp), allocatable :: ats(:,:), f(:,:), q(:)
        integer, allocatable :: el(:)
        ! only used with volumetric data
        real(dp), allocatable :: nValenceElectrons(:)
        ! only used with density fitting
        real(dp), allocatable :: rnuc(:)
        real(dp) :: qtot, energy
    end type atStruct

contains

    function as_new(nat, pos, el, nValenceElectrons, rnuc, lat) result(ats)
        integer, intent(in) :: nat
        real(dp), intent(in) :: pos(3, nat)
        integer, intent(in), optional :: el(nat)
        real(dp), intent(in), optional :: nValenceElectrons(nat)
        real(dp), intent(in), optional :: rnuc(nat)
        real(dp), intent(in), optional :: lat(3,3)
        type(atStruct) :: ats

        ats%nat = nat
        allocate(ats%ats(3,nat))
        ats%ats = pos
        if (present(lat)) then
            ats%periodic = .true.
            ats%lat = lat
        else
            ats%periodic = .false.
        end if
        if (present(nValenceElectrons)) then
            allocate(ats%nValenceElectrons(nat))
            ats%nValenceElectrons = nValenceElectrons
        end if
        if (present(rnuc)) then
            allocate(ats%rnuc(nat))
            ats%rnuc = rnuc
        end if
        if (present(el)) then
            allocate(ats%el(nat))
            ats%el = el
        end if

    end function as_new

    subroutine as_center(ats)
        type(atStruct) :: ats

        call center(ats%nat, ats%ats)

    end subroutine

    subroutine as_getCenter(ats, c)
        type(atStruct) :: ats
        real(dp), intent(out) :: c(3)
        c = sum(ats%ats, 2) / ats%nat
    end subroutine

    subroutine as_getCenterOfMass(ats, c)
        type(atStruct) :: ats
        real(dp), intent(out) :: c(3)
        integer :: i
        real(dp) :: mtot, m

        c(:) = 0._dp
        mtot = 0._dp
        do i=1,ats%nat
            m = getAtomicMass(ats%el(i))
            c(:) = c(:) + ats%ats(:,i) * m
            mtot = mtot + m
        end do

        c(:) = c(:) / mtot
    end subroutine

    subroutine as_centerMass(ats)
        type(atStruct) :: ats
        real(dp) :: c(3)

        call as_getCenterOfMass(ats, c)

        ats%ats(1,:) = ats%ats(1,:) - c(1)
        ats%ats(2,:) = ats%ats(2,:) - c(2)
        ats%ats(3,:) = ats%ats(3,:) - c(3)

    end subroutine

    subroutine as_copy(from, to)
        type(atStruct), intent(in) :: from
        type(atStruct), intent(out) :: to

        stop 'Error: as_copy is unnecessary and should no be used. Use = instead.'

        to%nat = from%nat
        if (allocated(to%ats)) then
            deallocate(to%ats)
        end if
        if (allocated(from%ats)) then
            allocate(to%ats(3,to%nat))
            to%ats(:,:) = from%ats(:,:)
        end if

        if (allocated(to%el)) then
            deallocate(to%el)
        end if
        if (allocated(from%el)) then
            allocate(to%el(to%nat))
            to%el(:) = from%el(:)
        end if

        if (allocated(to%nValenceElectrons)) then
            deallocate(to%nValenceElectrons)
        end if
        if (allocated(from%nValenceElectrons)) then
            allocate(to%nValenceElectrons(to%nat))
            to%nValenceElectrons(:) = from%nValenceElectrons(:)
        end if

        to%periodic = from%periodic
        to%lat = from%lat

    end subroutine

    subroutine as_readXYZ(filename, ats)
        character(len = *), intent(in) :: filename
        type(atStruct), intent(out) :: ats
        integer :: i, filestat
        character(len=2) :: el

        open(41, file = filename, iostat = filestat, status = "old", action = "read")
        if(filestat == 0) then
            !write(*,*) "sucessfully opened file"
            read(41, *) ats%nat
            read(41, *)
            allocate(ats%ats(3, ats%nat))
            allocate(ats%el(ats%nat))
            do i = 1, ats%nat
                read(41, *) el, ats%ats(:, i)
                ats%el(i) = elemSymToNum(el)
            end do
            ats%periodic = .false.
            close(41)
        else
            print*, 'Error: Could not open file (' // filename // ')'
            stop
        end if
    end subroutine as_readXYZ

    subroutine as_readFhiaims(filename, ats)
        character(len = *), intent(in) :: filename
        type(atStruct), intent(out) :: ats
        integer :: i, filestat, funit
        character(len=2) :: el
        real(dp) :: xyz(3)
        type(dynamicArray) :: tmpAts, tmpEl
        character(len=200) :: tmpstr, placeholder
        integer :: ilat
        real(dp) :: lat(3,3)
        real(dp), allocatable :: elems(:,:)
        integer :: nat


        nat = 0
        ilat = 1
        lat = 0._dp
        tmpAts = newDynamicArray(3)
        tmpEl = newDynamicArray(1)
        open(newunit=funit, file = filename, iostat = filestat, status = "old", action = "read")
        if(filestat /= 0) then
            stop "could not open file"
        end if

        read(funit, '(A)', iostat=filestat) tmpstr
        do while(filestat == 0)
            if (tmpstr(1:4)=='atom') then
                nat = nat + 1
                read(tmpstr,*) placeholder(:5), xyz(:), el
                call tmpAts%append(xyz)
                call tmpEl%append([real(elemSymToNum(el), dp)])

            elseif(tmpstr(1:14)=='lattice_vector') then
                if (ilat > 3) stop 'ERROR: more than 3 lattice vectors'
                read(tmpstr,*) placeholder, lat(:,ilat)
                ilat = ilat + 1
            end if
            read(funit, '(A)', iostat=filestat) tmpstr
        end do

        allocate(elems(1,nat))
        elems = tmpEl%get()
        ats%nat = nat
        allocate(ats%ats(3,nat))
        allocate(ats%el(nat))
        do i=1,nat
            ats%el(i) = int(elems(1,i))
        end do
        ats%ats = tmpAts%get() * ang2bohr
        if (ilat > 1) then
            ats%periodic = .true.
            ats%lat = lat * ang2bohr
        end if

    end subroutine

    subroutine as_writeXYZ(filename, ats)
        character(len = *), intent(in) :: filename
        type(atStruct), intent(in) :: ats

        integer :: i, filestat

        open(41, file = filename, iostat = filestat, action = "write")
        if(filestat == 0) then
            write(41, *) ats%nat
            write(41, *)
            do i = 1, ats%nat
                ! todo: test format statement
                write(41, '(A2,3E15.6)') elemNumToSym(ats%el(i)), ats%ats(:, i)
            end do
            close(41)
        else
            stop "could not open file"
        end if

    end subroutine as_writeXYZ

    subroutine as_readAscii(filename, ats, isAngstrom)
        character(len=*), intent(in) :: filename
        type(atStruct), intent(out) :: ats
        logical, intent(in), optional :: isAngstrom
        integer :: funit, iostat
        character(len=200) :: tmpstr, angstr
        character(len=2) :: elstr
        integer :: i
        logical :: angs

        angs = .false.
        if (present(isAngstrom)) angs = isAngstrom

        open(newunit=funit, file=filename)
        read(funit, '(A)') tmpstr
        read(tmpstr, *, iostat=iostat) ats%nat, angstr
        if (iostat /= 0) then
            read(tmpstr, *, iostat=iostat) ats%nat
        end if
        ats%periodic = .true.
        if (allocated(ats%ats)) deallocate(ats%ats)
        allocate(ats%ats(3,ats%nat))
        if (allocated(ats%el)) deallocate(ats%el)
        allocate(ats%el(ats%nat))

        ats%lat = 0._dp
        read(funit, *) ats%lat(:1, 1), ats%lat(:2, 2)
        read(funit, *) ats%lat(:3, 3)

        do i=1, ats%nat
            read(funit, *) ats%ats(:,i), elstr
            ats%el(i) = elemSymToNum(elstr)
        end do

        ! check if the units are angstroem
        if (angstr(:9) == 'angstroem' .or. angs) then
            ats%ats = ats%ats * ang2bohr
            ats%lat = ats%lat * ang2bohr
        end if
        close(funit)

    end subroutine

    subroutine as_writeAscii(filename, ats, toAngstrom)
        character(len=*), intent(in) :: filename
        type(atStruct), intent(in) :: ats
        logical, intent(in), optional :: toAngstrom
        real(dp) :: newlat(3,3), newats(3, ats%nat)
        integer :: funit, i
        logical :: angs

        angs = .false.
        if (present(toAngstrom)) then
            angs = toAngstrom
        end if

        newlat(:,:) = ats%lat(:,:)
        newats(:,:) = ats%ats(:,:)
        if (angs) then
            newlat = newlat / ang2bohr
            newats = newats / ang2bohr
        end if
        call alignLatticeWithAxes(ats%nat, newlat, newats)

        open(newunit=funit, file=filename)

        ! put nat in comment line. Not required
        if (angs) then
            write(funit, *) ats%nat, 'angstroem'
        else
            write(funit, *) ats%nat
        end if
        write(funit, *) newlat(:1,1), newlat(:2, 2)
        write(funit, *) newlat(:3, 3)

        do i=1,ats%nat
            write(funit,*) newats(:, i), elemNumToSym(ats%el(i))
        end do

        close(funit)


    end subroutine

    subroutine as_writeDataToUnit(ats, energy, charges, forces, unit)
        type(atStruct), intent(in) :: ats
        real(dp), intent(in) :: energy, charges(ats%nat), forces(3, ats%nat)
        integer, intent(in) :: unit
        integer :: i

        write(unit, *) 'begin'
        if (ats%periodic) then
            write(unit, *) 'lattice ', ats%lat(:,1)
            write(unit, *) 'lattice ', ats%lat(:,2)
            write(unit, *) 'lattice ', ats%lat(:,3)
        end if

        do i=1, ats%nat
            write(unit, *) 'atom ', ats%ats(:,i), elemNumToSym(ats%el(i)), charges(i), 0._dp, forces(:,i)
        end do

        write(unit, *) 'energy ', energy
        write(unit, *) 'charge ', nint(sum(charges))
        write(unit, *) 'end'


    end subroutine


    subroutine as_readData(fname, n, ats)
        character(len=*), intent(in) :: fname
        integer, intent(out) :: n
        type(atStruct), allocatable, intent(out) :: ats(:)
        type(dynamicClassArray) :: atsArr
        type(dynamicArray) :: posArr, qArr, fArr, elArr
        character(len=512) :: line, tmpstr
        character(len=2) :: el(1)
        real(dp) :: pos(3), f(3), q(1), e(1)

        integer :: unit, iostat, ilat, i
        type(atStruct) :: tmpAts, freeAts

        atsArr = newDynamicClassArray()

        open(file=fname, newunit=unit)

        read(unit, '(A)', iostat=iostat) line
        do while (iostat==0)
            line = trim(adjustl(line))
            if (line(:5) == 'begin') then
                tmpAts = freeAts ! this way we dont have to deallocate all the arrays manually
                posArr = newDynamicArray(3)
                qArr = newDynamicArray(1)
                fArr = newDynamicArray(3)
                elArr = newDynamicArray(1)
                ilat = 1
            end if
            if (line(:7) == 'lattice') then
                read(line,*) tmpstr, tmpAts%lat(:,ilat)
                ilat = ilat + 1
            end if
            if (line(:4) == 'atom') then
                read(line,*) tmpstr, pos, el, q, e, f
                call posArr%append(pos)
                call qArr%append(q)
                call fArr%append(f)
                call elArr%append([real(elemSymToNum(el), dp)])
            end if
            if (line(:6) == 'energy') then
                read(line,*) tmpstr, tmpAts%energy
            end if
            if (line(:6) == 'charge') then
                read(line,*) tmpstr, tmpAts%qtot
            end if
            if (line(:3) == 'end') then
                tmpAts%nat = posArr%nData
                allocate(tmpAts%ats(3, tmpAts%nat))
                allocate(tmpAts%f(3, tmpAts%nat))
                allocate(tmpAts%el(tmpAts%nat))
                allocate(tmpAts%q(tmpAts%nat))
                tmpAts%periodic = ilat == 4
                tmpAts%ats = posArr%get()
                tmpAts%el = int(reshape(elArr%get(), shape(tmpAts%el)))
                tmpAts%f = fArr%get()
                tmpAts%q = reshape(qArr%get(), shape(tmpAts%q))
                call atsArr%append(tmpAts)
            end if
            read(unit, '(A)', iostat=iostat) line
        end do
        close(unit)

        n = atsArr%nData
        allocate(ats(n))
        do i=1,n
            select type(sel => atsArr%data(i)%val)
                type is (atStruct)
                    ats(i) = sel
            end select
        end do
    end subroutine


    subroutine as_readPOSCAR(fname, ats)
        character(len=*), intent(in) :: fname
        type(atStruct), intent(out) :: ats
        integer :: unit, i, j, c
        real(dp) :: scale
        character(len=256) :: line
        integer, allocatable :: nels(:)
        integer :: n
        character(len=2), allocatable :: els(:)

        open(file=fname, newunit=unit, action='read')
        read(unit, *)
        read(unit, *) scale
        do i=1,3
            read(unit, *) ats%lat(:,i)
        end do
        ats%lat = ats%lat * ang2bohr
        ats%periodic=.true.

        read(unit, '(A)') line
        n = getnCols(line)
        allocate(els(n))
        allocate(nels(n))
        read(line, *) els
        read(unit, *) nels
        read(unit, '(A)') line
        call lowercase(line)
        if (adjustl(trim(line)) /= 'direct') then
            print*, line
            print*, 'ERROR: POSCAR without DIRECT not implemented'
            stop
        end if

        ats%nat = sum(nels)
        allocate(ats%el(ats%nat))
        allocate(ats%ats(3, ats%nat))

        c = 0
        do i=1,n
            do j=1,nels(i)
                c = c+1
                ats%el(c) = elemSymToNum(els(i))
            end do
        end do

        do i=1,ats%nat
            read(unit, *) ats%ats(:,i)
        end do

        call toAbsoluteCoordinates(ats%nat, ats%lat, ats%ats)

        close(unit)
    end subroutine

    subroutine as_writePOSCAR(fname, ats)
        character(len=*), intent(in) :: fname
        type(atStruct), intent(in) :: ats
        integer :: unit, i, iel
        integer :: nels
        integer, allocatable :: els(:), elcount(:)
        real(dp) :: relats(3, ats%nat)

        open(file=fname, newunit=unit, action='write')

        write(unit, *) ''
        write(unit, *) 1._dp

        do i=1,3
            write(unit, *) ats%lat(:,i) / ang2bohr
        end do

        call getDistinctElements(ats%nat, ats%el, nels, els, elcount)
        write(unit,*) elemNumToSym(els)
        write(unit,*) elcount

        write(unit, *) 'DIRECT'

        relats = ats%ats
        call toRelativeCoordinates(ats%nat, ats%lat, relats)

        do iel=1,nels
            do i=1,ats%nat
                if (ats%el(i) == els(iel)) then
                    write(unit, *) relats(:,i)
                end if
            end do
        end do

        close(unit)
    end subroutine

    subroutine as_moveAtomsIntoCell(ats)
        type(atStruct), intent(inout) :: ats

        call moveAtomsIntoCell(ats%nat, ats%lat, ats%ats)

    end subroutine


    subroutine as_elemCount(ats, nelem)
        type(atStruct), intent(in) :: ats
        integer, intent(out) :: nelem(maxElemNum)
        integer :: i

        nelem(:) = 0

        do i=1,ats%nat
            nelem(ats%el(i)) = nelem(ats%el(i)) + 1
        end do

    end subroutine

    subroutine as_reassign(ats, assignment)
        type(atStruct), intent(inout) :: ats
        integer, intent(in) :: assignment(ats%nat)
        real(dp) :: tmp(3, ats%nat)
        integer :: tmpEl(ats%nat)
        integer :: i

        tmp = ats%ats
        tmpEl = ats%el

        do i=1,ats%nat
            ats%ats(:,i) = tmp(:,assignment(i))
            ats%el(i) = tmpEl(assignment(i))
        end do

    end subroutine as_reassign

    subroutine as_reassignInv(ats, assignment)
        type(atStruct), intent(inout) :: ats
        integer, intent(in) :: assignment(ats%nat)
        real(dp) :: tmp(3, ats%nat)
        integer :: tmpEl(ats%nat)
        integer :: i

        tmp = ats%ats
        tmpEl = ats%el

        do i=1,ats%nat
            ats%ats(:,assignment(i)) = tmp(:,i)
            ats%el(assignment(i)) = tmpEl(i)
        end do

    end subroutine as_reassignInv

    ! determines if two structures are identical up to the atomic coordinates
    function as_compare(atsA, atsB) result(ret)
        type(atStruct) :: atsA, atsB
        logical :: ret
        integer :: nElA(maxElemNum), nElB(maxElemNum)

        ret = .true.

        if (atsA%nat /= atsB%nat) then
            ret = .false.
            return
        end if

        call as_elemCount(atsA, nElA)
        call as_elemCount(atsB, nElB)

        if (.not. arrayCompare(nElA, nElB)) then
            ret = .false.
            return
        end if

    end function as_compare

    ! assuming all atoms are in unit cell (only nearest neighbour replica cells are considered.
    function as_shortestDistance(ats, i, j, c) result(d)
        type(atStruct), intent(in) :: ats
        integer, intent(in) :: i, j
        integer, optional, intent(out) :: c(3)
        real(dp) :: d, dtmp
        integer :: x, y, z
        type(atStruct) :: tmpAts
        d = huge(d)

        tmpAts = ats

        if (ats%periodic) then
            !call as_moveAtomsIntoCell(tmpAts)
            do x = -1,1
                do y = -1,1
                    do z = -1,1
                        dtmp = sqrt(sum((tmpAts%ats(:,i) - tmpAts%ats(:,j) &
                                + matmul(tmpAts%lat, [1._dp * x, 1._dp * y, 1._dp * z]))**2))
                        if (dtmp < d) then
                            d = dtmp
                            if (present(c)) c = [x, y, z]
                        end if
                    end do
                end do
            end do
        else
            d = sqrt(sum((tmpAts%ats(:,i) - tmpAts%ats(:,j))**2))
            if (present(c)) c = 0
        end if
    end function

    function as_getRepetitionIndex(n, i)  result(index)
        integer, intent(in) :: n(3), i(3)
        integer :: index

        index = i(3) + (i(1) - 1) * n(2) * n(3) + (i(2) - 1) * n(3)

    end function

    subroutine as_repeatUnitCell(ats, n)
        type(atStruct), intent(inout) :: ats
        integer, intent(in) :: n(3)
        real(dp), allocatable :: atsTmp(:,:)
        integer, allocatable :: elTmp(:)
        integer :: x, y, z, i, c

        allocate(atsTmp(3,ats%nat * product(n)))
        allocate(elTmp(ats%nat * product(n)))

        c = 0
        do x=1,n(1)
            do y=1,n(2)
                do z=1,n(3)
                    do i=1,ats%nat
                        elTmp(c*ats%nat+i) = ats%el(i)
                        atsTmp(:,c*ats%nat+i) = &
                                ats%ats(:,i) + (x-1)*ats%lat(:,1) + (y-1)*ats%lat(:,2) + (z-1)*ats%lat(:,3)
                    end do
                    c = c + 1
                    ! print*, c, as_getRepetitionIndex(n, [x, y, z]) ! must be equal
                end do
            end do
        end do
        ats%lat(:,1) = ats%lat(:,1) * n(1)
        ats%lat(:,2) = ats%lat(:,2) * n(2)
        ats%lat(:,3) = ats%lat(:,3) * n(3)
        deallocate(ats%ats)
        deallocate(ats%el)
        allocate(ats%ats(3,ats%nat * product(n)))
        allocate(ats%el(ats%nat * product(n)))
        ats%ats = atsTmp
        ats%el = elTmp
        ats%nat = ats%nat * product(n)

    end subroutine as_repeatUnitCell

    subroutine as_translate(ats, t)
        type(atStruct), intent(inout) :: ats
        real(dp), intent(in) :: t(3)
        call translate(ats%nat, ats%ats, t)
    end subroutine as_translate

    subroutine dEdLat2Stress(lat, dEdlat, stress)
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(out) :: stress(3,3)
        real(dp), intent(in) :: dEdlat(3,3)

        stress = -1._dp / det3D(lat) * matmul(dEdlat, transpose(lat))

    end subroutine

    subroutine stress2dEdLat(lat, stress, dEdlat)
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(in) :: stress(3,3)
        real(dp), intent(out) :: dEdlat(3,3)
        real(dp) :: invlat(3,3)

        call inv3DM(lat, invlat)
        dEdlat = -1._dp * det3D(lat) * matmul(stress, transpose(invlat))

    end subroutine

end module atomicStructure