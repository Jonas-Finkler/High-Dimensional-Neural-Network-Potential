! Created by  on 28.01.20.

program sortTest
    use precision
    use sorting

    implicit none

    integer, parameter :: n = 10000
    real(dp), allocatable :: v(:)
    integer, allocatable ::  index(:)
    integer :: i
    real(dp) :: st, et

    allocate(v(n))
    allocate(index(n))

    call random_number(v)

    call cpu_time(st)
    call mergesort(n, v, index)
    call cpu_time(et)

    do i=1,n-1
        if (v(index(i)) > v(index(i+1))) write(*,*) 'ERROR', i
    end do

    write(*,*) et-st

end program sortTest