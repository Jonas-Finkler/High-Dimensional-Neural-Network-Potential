! Created by  on 28.01.20.

module sorting
    use precision
    implicit none

contains

    subroutine mergeSortInPlace(n, vals)
        integer, intent(in) :: n
        real(dp), intent(inout) :: vals(n)
        integer :: index(n)
        integer :: i
        real(dp) :: tmpVals(n)


        call mergesort(n, vals, index)
        tmpVals = vals
        do i=1,n
            vals(i) = tmpVals(index(i))
        end do
    end subroutine

    ! vals(index(i)) will be the ith larges element
    ! mergesort is stable, order of equal elements is preserved
    subroutine mergesort(n, vals, index)
        integer, intent(in) :: n
        real(dp), intent(in) :: vals(n)
        integer, intent(out) :: index(n)
        integer :: i
        real(dp), allocatable :: v(:)

        allocate(v(n))

        v(:) = vals(:)

        do i=1,n
            index(i) = i
        end do

        call mergesort_rec(n, n, vals, index)

    end subroutine mergesort

    ! maybe not very efficient, because referencing v through index is not cache friendly.
    ! the allocate makes it a bit slower, but allows for much larger datasets
    recursive subroutine mergesort_rec(n, nv, v, index)
        integer, intent(in) :: n, nv
        real(dp), intent(in) :: v(nv)
        integer, intent(inout) :: index(n)
        integer, allocatable :: tmpindex(:)
        integer :: p, itmp, i, a, b

        if (n == 1) return

        ! sort
        if (n == 2) then
            if (v(index(1)) > v(index(2))) then
                itmp = index(1); index(1) = index(2); index(2) = itmp
            end if
            return
        end if

        p = n / 2

        allocate(tmpindex(n))
        tmpindex(:) = index(:)

        call mergesort_rec(p, nv, v, tmpindex(:p))
        call mergesort_rec(n-p, nv, v, tmpindex(p+1:))

        ! merge
        a = 1
        b = p+1

        do i=1,n
            if (a>p) then
                index(i) = tmpindex(b)
                b = b + 1
                cycle
            end if
            if (b>n) then
                index(i) = tmpindex(a)
                a = a + 1
                cycle
            end if
            if(v(tmpindex(b)) < v(tmpindex(a))) then
                index(i) = tmpindex(b)
                b = b + 1
            else
                index(i) = tmpindex(a)
                a = a + 1
            end if
        end do

        deallocate(tmpindex)

    end subroutine mergesort_rec


    ! sorting for string. Basically just a copy paste of the above

    ! vals(index(i)) will be the ith larges element
    ! mergesort is stable, order of equal elements is preserved
    subroutine mergesortStr(n, vals, index)
        integer, intent(in) :: n
        character(len=*), intent(in) :: vals(n)
        integer, intent(out) :: index(n)
        integer :: i
        character(len=len(vals)), allocatable :: v(:)

        allocate(v(n))

        v(:) = vals(:)

        do i=1,n
            index(i) = i
        end do

        call mergesort_rec_str(n, n, vals, index)

    end subroutine mergesortStr

    ! maybe not very efficient, because referencing v through index is not cache friendly.
    ! the allocate makes it a bit slower, but allows for much larger datasets
    recursive subroutine mergesort_rec_str(n, nv, v, index)
        integer, intent(in) :: n, nv
        character(len=*), intent(in) :: v(nv)
        integer, intent(inout) :: index(n)
        integer, allocatable :: tmpindex(:)
        integer :: p, itmp, i, a, b

        if (n == 1) return

        ! sort
        if (n == 2) then
            if (v(index(1)) > v(index(2))) then
                itmp = index(1); index(1) = index(2); index(2) = itmp
            end if
            return
        end if

        p = n / 2

        allocate(tmpindex(n))
        tmpindex(:) = index(:)

        call mergesort_rec_str(p, nv, v, tmpindex(:p))
        call mergesort_rec_str(n-p, nv, v, tmpindex(p+1:))

        ! merge
        a = 1
        b = p+1

        do i=1,n
            if (a>p) then
                index(i) = tmpindex(b)
                b = b + 1
                cycle
            end if
            if (b>n) then
                index(i) = tmpindex(a)
                a = a + 1
                cycle
            end if
            if(v(tmpindex(b)) < v(tmpindex(a))) then
                index(i) = tmpindex(b)
                b = b + 1
            else
                index(i) = tmpindex(a)
                a = a + 1
            end if
        end do

        deallocate(tmpindex)

    end subroutine mergesort_rec_str
end module sorting