! Created by jonas on 16.03.21.

! Some useful code to parse parameters from a file

module parameterParser
    use precision
    use sorting
    implicit none

    integer, parameter :: paramParserStrLen = 512
    type parameterBuffer
        integer :: n
        character(len=paramParserStrLen), allocatable :: names(:)
        character(len=paramParserStrLen), allocatable :: values(:)
        integer, allocatable :: index(:)
    end type parameterBuffer

    interface getParameter
        module procedure getParameterString
        module procedure getParameterInt
        module procedure getParameterReal
        module procedure getParameterLogical
    end interface getParameter

contains

    ! todo: error message if same prameter is given twice
    subroutine readParameterFile(fname, buffer)
        character(len=*), intent(in) :: fname
        type(parameterBuffer), intent(out) :: buffer
        integer :: unit, fstat, i
        character(len=paramParserStrLen) :: name, value
        character(len=paramParserStrLen), allocatable :: tmpStrs(:)
        character(len=paramParserStrLen*10) :: tmpstr
        integer :: start, split

        logical :: found

        buffer%n = 1
        allocate(buffer%names(buffer%n))
        allocate(buffer%values(buffer%n))

        open(newunit=unit, file=fname)
        i = 0
        do
            read(unit, '(A)', iostat=fstat) tmpstr !name, value
            if (fstat /= 0) exit
            start = index(tmpstr, '#')
            if (start > 1) then
                tmpstr = trim(tmpstr(:start-1))
            end if
            split = index(tmpstr, '=')
            if (split == 0) cycle
            name = adjustl(trim(tmpstr(:split-1)))
            value = adjustl(trim(tmpstr(split+1:)))
            !print*, i, trim(name), ':', trim(value)
            !if (name(:1) == '#') cycle
            i = i + 1
            buffer%names(i) = name
            buffer%values(i) = value
            if (i == buffer%n) then
                allocate(tmpStrs(buffer%n))
                tmpStrs(:) = buffer%names(:)
                deallocate(buffer%names)
                allocate(buffer%names(buffer%n * 2))
                buffer%names(:buffer%n) = tmpStrs
                tmpStrs(:) = buffer%values(:)
                deallocate(buffer%values)
                allocate(buffer%values(buffer%n * 2))
                buffer%values(:buffer%n) = tmpStrs
                deallocate(tmpStrs)
                buffer%n = buffer%n * 2
            end if
        end do
        close(unit)
        buffer%n = i

        allocate(buffer%index(buffer%n))
        call mergesortStr(buffer%n, buffer%names, buffer%index)

    end subroutine readParameterFile

    subroutine printParameters(buffer)
        type(parameterBuffer), intent(in) :: buffer
        integer :: i

        do i=1,buffer%n
            print*, trim(buffer%names(buffer%index(i))), ' = ', trim(buffer%values(buffer%index(i)))
        end do
    end subroutine

    ! binary search in the sorted index
    subroutine getParamIndex(buffer, name, idx, found)
        type(parameterBuffer), intent(in) :: buffer
        character(len=*), intent(in) :: name
        integer, intent(out) :: idx
        logical, intent(out) :: found
        integer :: imax, imin

        found = .false.

        imax = buffer%n
        imin = 1
        idx = buffer%n / 2
        if (buffer%n<2) then
            if (buffer%n>1) then
                idx = 1
                found = trim(name) == trim(buffer%names(1))
            end if
            return
        end if

        do while(imax /= imin)
            if (trim(name) > trim(buffer%names(buffer%index(idx)))) then
                imin = idx + 1
            else
                imax = idx
            end if
            idx = imin + (imax - imin) / 2
            !print*, name, imin, idx, imax
            if (trim(name) == trim(buffer%names(buffer%index(idx)))) then
                found = .true.
                idx = buffer%index(idx)
                return
            end if
        end do

    end subroutine getParamIndex

    subroutine getParameterString(buffer, name, res)
        type(parameterBuffer), intent(in) :: buffer
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: res

        integer :: idx
        logical :: found

        call getParamIndex(buffer, name, idx, found)
        if (.not. found) then
            print*, 'Error: Parameter ', trim(name), ' not found in parameter file'
            stop
        end if

        res = trim(buffer%values(idx))

    end subroutine getParameterString

    subroutine getParameterInt(buffer, name, res)
        type(parameterBuffer), intent(in) :: buffer
        character(len=*), intent(in) :: name
        integer, intent(out) :: res
        integer :: base, exponent
        integer :: epos

        integer :: idx
        logical :: found

        call getParamIndex(buffer, name, idx, found)
        if (.not. found) then
            print*, 'Error: Parameter ', trim(name), ' not found in parameter file'
            stop
        end if

        epos = index(buffer%values(idx), '**')
        if (epos==0) then
            read(buffer%values(idx), *) res
        else
            read(buffer%values(idx)(:epos-1), *) base
            read(buffer%values(idx)(epos+2:), *) exponent
            !print*, 'exp format', base, exponent, base**exponent
            res = base**exponent
        end if

    end subroutine getParameterInt

    subroutine getParameterReal(buffer, name, res)
        type(parameterBuffer), intent(in) :: buffer
        character(len=*), intent(in) :: name
        real(dp), intent(out) :: res

        integer :: idx
        logical :: found

        call getParamIndex(buffer, name, idx, found)
        if (.not. found) then
            print*, 'Error: Parameter ', trim(name), ' not found in parameter file'
            stop
        end if

        read(buffer%values(idx), *) res

    end subroutine getParameterReal

    subroutine getParameterLogical(buffer, name, res)
        type(parameterBuffer), intent(in) :: buffer
        character(len=*), intent(in) :: name
        logical, intent(out) :: res

        integer :: idx
        logical :: found

        call getParamIndex(buffer, name, idx, found)
        if (.not. found) then
            print*, 'Error: Parameter ', trim(name), ' not found in parameter file'
            stop
        end if

        read(buffer%values(idx), *) res

    end subroutine getParameterLogical
end module parameterParser