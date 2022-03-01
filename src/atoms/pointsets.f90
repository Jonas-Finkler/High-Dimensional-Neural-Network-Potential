! Written by Jonas A. Finkler

! routiens for atomic structures without element information (pointsets)


module pointsets
    use precision
    use linalg
    use dynamicArrays
    implicit none



contains

    !todo: implement reading elements. -> use as_readXyz

    subroutine readXYZ(filename, nat, ats)
        character(len = *), intent(in) :: filename
        real(dp), intent(out), dimension(:, :), allocatable :: ats
        integer, intent(out):: nat
        integer :: i, filestat
        character(len = 2) :: atType

        open(41, file = filename, iostat = filestat, status = "old", action = "read")
        if(filestat == 0) then
            !write(*,*) "sucessfully opened file"
            read(41, *) nat
            read(41, *)
            if (allocated(ats)) deallocate(ats)
            allocate(ats(3, nat))
            do i = 1, nat
                read(41, *) atType, ats(:, i)
            end do
            close(41)
        else
            stop "could not open file"
        end if
    end subroutine readXYZ

    subroutine nreadXYZ(filename, nat, ats)
        character(len = *), intent(in) :: filename
        integer, intent(in):: nat
        real(dp), intent(out) :: ats(3,nat)
        integer :: i, filestat
        character(len = 2) :: atType

        open(41, file = filename, iostat = filestat, status = "old", action = "read")
        if(filestat == 0) then
            read(41, *)
            read(41, *)
            do i = 1, nat
                read(41, *) atType, ats(:, i)
            end do
            close(41)
        else
            stop "could not open file"
        end if
    end subroutine nreadXYZ

    subroutine writeXYZ(filename, nat, ats)
        character(len = *), intent(in) :: filename
        real(dp), intent(in), dimension(:, :) :: ats
        integer, intent(in) :: nat

        integer :: i, filestat

        open(41, file = filename, iostat = filestat, action = "write")
        if(filestat == 0) then
            write(41, *) nat
            write(41, *)
            do i = 1, nat
                write(41, *) "lj", ats(:, i)
            end do
            close(41)
        else
            stop "could not open file"
        end if

    end subroutine writeXYZ

    subroutine readAscii(filename, nat, lat, ats)
        character(len=*), intent(in) :: filename
        integer, intent(out) :: nat
        real(dp), intent(out) :: lat(3,3)
        real(dp), intent(out), dimension(:,:), allocatable :: ats
        integer :: funit, filestat
        character(len=100) :: tmpstr
        real(dp) :: tmplat(6), tmpats(3)
        type(dynamicArray) :: dats

        dats = newDynamicArray(3)

        open(newunit=funit, file=filename)
        read(funit, *)
        read(funit, *) tmplat(1:3)
        read(funit, *) tmplat(4:6)
        lat(:,:) = 0._dp
        lat(1,  1) = tmplat(1)
        lat(1:2,2) = tmplat(2:3)
        lat(1:3,3) = tmplat(4:6)
        write(*,*) lat

        filestat = 0

        read(funit, *, iostat=filestat) tmpats(:), tmpstr
        do while(filestat == 0)
            call dats%append(tmpats)
            read(funit, *, iostat=filestat) tmpats(:), tmpstr
        end do

        nat = dats%ndata

        allocate(ats(3,nat))
        ats(:,:) = dats%get()

        close(funit)

    end subroutine

    subroutine writeAscii(filename, nat, lat, ats)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: nat
        real(dp), intent(in) :: lat(3,3), ats(3, nat)
        real(dp) :: newlat(3,3), newats(3, nat)

        integer :: funit, i

        newlat(:,:) = lat(:,:)
        newats(:,:) = ats(:,:)
        call alignLatticeWithAxes(nat, newlat, newats)

        open(newunit=funit, file=filename)

        ! put nat in comment line. Not required
        write(funit, *) nat
        write(funit, *) newlat(:1,1), newlat(:2, 2)
        write(funit, *) newlat(:3, 3)

        do i=1,nat
            write(funit,*) newats(:, i), 'Si' ! todo: use formating and elements
        end do

        close(funit)

    end subroutine

    subroutine alignLatticeWithAxes(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: lat(3, 3), ats(3, nat)
        real(dp) :: newlat(3,3), R(3,3)
        integer :: i
        real(dp) :: rotdet

        ! calculate lat = R * newlat
        ! then newlat = R^t * lat
        call qrDecomposition(3, lat, R, newlat)
        rotdet = R(1,1) + R(2,2) + R(3,3)
        if (rotdet < 0._dp) then ! we have an improper rotation -> fix it
            R(:,:) = -1._dp * R(:,:)
            newlat(:,:) = -1._dp * newlat(:,:)
        end if

        R = transpose(R)

        do i=1, nat
            ats(:,i) = matMulVec3(R, ats(:,i))
        end do

        lat(:,:) = newlat(:,:)

    end subroutine

    subroutine toRelativeCoordinates(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: ats(3, nat)
        real(dp), intent(in) :: lat(3,3)
        real(dp) :: invlat(3,3)
        integer :: i

        call inv3DM(lat, invlat)

        do i=1, nat
            ats(:,i) = mAtMulVec3(invlat, ats(:,i))
        end do

    end subroutine

    subroutine toAbsoluteCoordinates(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(inout) :: ats(3,nat)
        integer :: i

        do i=1, nat
            ats(:,i) = matMulVec3(lat, ats(:,i))
        end do

    end subroutine

    subroutine moveAtomsIntoCell(nat, lat, ats)
        integer, intent(in) :: nat
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(inout) :: ats(3, nat)
        integer :: i

        call toRelativeCoordinates(nat, lat, ats)

        do i=1, nat
            ats(1,i) = modulo(ats(1,i), 1._dp)
            ats(2,i) = modulo(ats(2,i), 1._dp)
            ats(3,i) = modulo(ats(3,i), 1._dp)
        end do

        call toAbsoluteCoordinates(nat, lat, ats)

    end subroutine

    ! returns the number of cells in the direction of each lattice vector that are needed to find all neighbors within cutoff
    subroutine getNcells(lattice, cutoff, n)
        implicit none
        real(dp), intent(in) :: lattice(3,3)
        real(dp), intent(in) :: cutoff
        integer, intent(out) :: n(3)
        real(dp) :: axb(3), axc(3), bxc(3)
        real(dp) :: proj(3)


        ! calculate the norm vectors of the 3 planes
        axb = cross(lattice(:,1), lattice(:,2))
        axb = axb / sqrt(sum(axb**2))
        axc = cross(lattice(:,1), lattice(:,3))
        axc = axc / sqrt(sum(axc**2))
        bxc = cross(lattice(:,2), lattice(:,3))
        bxc = bxc / sqrt(sum(bxc**2))

        proj(1) = dot(lattice(:,1), bxc(:))
        proj(2) = dot(lattice(:,2), axc(:))
        proj(3) = dot(lattice(:,3), axb(:))

        n(:) = ceiling(cutoff / abs(proj))

    end subroutine getNcells

    subroutine center(nat, ats)
        integer, intent(in) :: nat
        real(dp), intent(inout), dimension(:, :) :: ats
        real(dp), dimension(3) :: cm
        integer :: i

        cm = sum(ats, 2) / nat

        do i = 1, nat
            ats(:, i) = ats(:, i) - cm(:)
        end do
    end subroutine center

    function abs2(nat, ats)
        integer, intent(in) :: nat
        real(dp) :: abs2
        real(dp), intent(in), dimension(3,nat) :: ats
        abs2 = sum(ats**2)
    end function abs2

    function abs2p(ats)
        real(dp) :: abs2p
        real(dp), intent(in), dimension(3) :: ats
        abs2p = sum(ats**2)
    end function abs2p

    !function abs(ats)
    !    real(dp) :: abs
    !    real(dp), intent(in), dimension(:, :) :: ats
    !    abs = sqrt(sum(ats**2))
    !end function abs

    subroutine translate(nat, atoms,r)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: atoms(3, nat)
        real(dp), intent(in) :: r(3)
        integer :: i

        do i=1,nat
            atoms(:,i) = atoms(:,i) + r(:)
        end do
    end subroutine

    function calcRmsd(nat, atsA, atsB) result(rmsd)
        integer, intent(in) :: nat
        real(dp) :: rmsd
        real(dp), intent(in) :: atsA(3,nat), atsB(3,nat)

        rmsd = sqrt(sum((atsA - atsB)**2 ) / nat)
    end function calcRmsd

    subroutine reassign(nat, ats, assignment)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: ats(3,nat)
        integer, intent(in) :: assignment(nat)
        real(dp) :: tmp(3, nat)
        integer :: i

        tmp = ats
        do i=1,nat
            ats(:,i) = tmp(:,assignment(i))
        end do
    end subroutine reassign

    subroutine reassignInv(nat, ats, assignment)
        integer, intent(in) :: nat
        real(dp), intent(inout) :: ats(3,nat)
        integer, intent(in) :: assignment(nat)
        real(dp) :: tmp(3, nat)
        integer :: i

        tmp = ats
        do i=1,nat
            ats(:,assignment(i)) = tmp(:,i)
        end do
    end subroutine reassignInv

    !points are gaussian distributed. (this is quite expensive (lots of math))
    subroutine randomPointSet(nat, ats)
        use random
        integer, intent(in) :: nat
        real(dp), intent(out) :: ats(3, nat)
        !real(dp) :: PI = 4._dp * datan(1._dp)
        !real(dp):: randN(int((3*size(ats,2)+1)/2),2), gaussians(int((3*size(ats,2)+1)/2)*2)
        integer :: i,j!n, i, np

        do j=1,nat
            do i=1,3
                ats(i,j) = random_normal()
            end do
        end do

        !n = size(ats,2)
        !np = n*3
        !if(mod(np,2) /= 0) np = np + 1

        !call random_number(randN)
        !gaussians(:np/2) = dsqrt(-2._dp * dlog(randN(:,1))) * dcos(2._dp * PI * randN(:,2))
        !gaussians(np/2+1:) = dsqrt(-2._dp * dlog(randN(:,1))) * dsin(2 * PI * randN(:,2))
        !ats = reshape(gaussians(:3*n),[3,n])

    end subroutine randomPointSet

    subroutine decomposePointSet(nat, ps, basis, decomposition)
        integer, intent(in) :: nat
        real(dp), intent(in) :: ps(3,nat), basis(3, nat, 3*nat)
        real(dp), intent(out) :: decomposition(3*nat)
        integer :: i

        do i=1,3*nat
            decomposition(i) = sum(basis(:,:,i)*ps(:,:))
        end do
    end subroutine decomposePointSet

    subroutine composePointSet(nat, ps, basis, decomposition)
        integer, intent(in) :: nat
        real(dp), intent(out) :: ps(3,nat)
        real(dp), intent(in) :: decomposition(3*nat), basis(3, nat, 3*nat)
        integer :: i

        ps = 0._dp
        do i=1,3*nat
            ps(:,:) = ps(:,:) + decomposition(i) * basis(:,:,i)
        end do
    end subroutine composePointSet

end module

