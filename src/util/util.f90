! Created by  on 5/19/20.

! this is the place to put small utility functions
module util
    use precision
    implicit none

contains

    subroutine getDistinctElements(n, A, nel, els, count)
        integer, intent(in) :: n
        integer, intent(in) :: A(n)
        integer, intent(out) :: nel
        integer, intent(out), allocatable :: els(:)
        integer, intent(out), optional, allocatable :: count(:)

        integer :: i, j
        integer :: tmp(n), tmpCount(n)
        logical :: found

        tmp(:) = -1
        tmpCount(:) = 0
        nel = 0
        do i=1,n
            found = .false.
            do j=1,nel
                if (A(i) == tmp(j)) then
                    found = .true.
                    tmpCount(j) = tmpCount(j) + 1
                    exit
                end if
            end do
            if (.not. found) then
                nel = nel + 1
                tmp(nel) = A(i)
                tmpCount(nel) = 1
            end if
        end do

        if (allocated(els)) deallocate(els)
        allocate(els(nel))
        els(:) = tmp(:nel)

        if (present(count)) then
            if(allocated(count)) deallocate(count)
            allocate(count(nel))
            count(:) = tmpCount(:nel)
        end if

    end subroutine getDistinctElements

    function arrayCompare(a, b)
        logical :: arrayCompare
        integer, intent(in) :: a(:), b(:)
        integer :: i
        arrayCompare = .false.
        if(size(a) /= size(b)) return
        do i = 1, size(a)
            if(a(i)/=b(i)) return
        end do
        arrayCompare = .true.
    end function arrayCompare

    subroutine assert(condition, errorMessage)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: errorMessage

        if(.not. condition) then
            write(*,*) 'assertion failed: ', errorMessage
            stop 
        end if

    end subroutine assert

    function getnCols(txt) result(n)
        character(len=*), intent(in) :: txt
        integer :: n
        logical :: isIn
        integer :: i


        isIn = .false.
        n = 0
        do i=1,len(trim(txt))
            if (isIn .and. txt(i:i) == ' ') then
                n = n + 1
            end if
            isIn = txt(i:i) /= ' '
        end do
        if (isIn) then
            n = n + 1
        end if
    end function getnCols

    ! arr has to contain elements in increasing order
    ! i will be the index of the largest element smaller than x
    ! if one of the elements is equal to x, i will be its index
    subroutine binarySearch(n, arr, x, i)
        integer, intent(in) :: n
        real(dp), intent(in) :: arr(n)
        real(dp), intent(in) :: x
        integer, intent(out) :: i
        integer :: imin, imax, imid

        if (arr(1) > x) then
            i = 0
            return
        end if
        if (arr(n) <= x) then
            i = n
            return
        end if
        imin = 1
        imax = n
        imid = (imax + imin) / 2
        do while (imax - imin > 1)
            if (arr(imid) <= x) then
                imin = imid
            else
                imax = imid
            end if
            imid = (imax + imin) / 2
        end do
        i = imin

    end subroutine binarySearch

    subroutine splitString(txt, n, words)
        character(len=*), intent(in) :: txt
        integer, intent(out) :: n
        character(len=256), allocatable, intent(out) :: words(:)
        integer :: i, pos, istart
        logical :: in

        n = getnCols(txt)
        if (allocated(words)) deallocate(words)
        allocate(words(n))

        in = .false.
        pos = 0
        do i=1,len(trim(txt))
            if (txt(i:i) == ' ' .or. txt(i:i) == '\t') then
                if (in) then
                    words(pos) = txt(istart:i)
                    in = .false.
                end if
            else
                if (.not. in) then
                    pos = pos + 1
                    istart = i
                    in = .true.
                end if
            end if

        end do

        if (in) then
            words(pos) = txt(istart:len(trim(txt)))
        end if

    end subroutine splitString

    elemental subroutine lowercase(word)
        character (len=*) , intent(inout) :: word
        integer                            :: i,ic,nlen
        nlen = len(word)
        do i=1,nlen
           ic = ichar(word(i:i))
           if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
        end do
    end subroutine lowercase


end module util

