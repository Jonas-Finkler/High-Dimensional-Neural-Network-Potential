! Created by  on 29.01.20.

module dynamicIntArrays

    use precision
    implicit none

    type dynamicIntArray
        integer :: ndata
        integer :: size
        Integer, allocatable :: data(:)
    contains
        procedure :: append => dynamicIntArray_append
        procedure :: deallocate => dynamicIntArray_deallocate
        procedure :: get => dynamicIntArray_get
    end type dynamicIntArray

contains

    function newDynamicIntArray() result(da)
        type(dynamicIntArray) :: da

        da%ndata = 0
        da%size = 1
        allocate(da%data(da%size))

    end function newDynamicIntArray

    subroutine dynamicIntArray_append(da, elem)
        class(dynamicIntArray), intent(inout) :: da
        integer, intent(in) :: elem
        integer, allocatable :: tmp(:)

        if (da%ndata + 1 > da%size) then
            allocate(tmp(da%ndata))
            tmp(:) = da%data(:da%ndata)
            deallocate(da%data)
            da%size = da%size * 2
            allocate(da%data(da%size))
            da%data(:da%ndata) = tmp(:)
            deallocate(tmp)
        end if

        da%ndata = da%ndata + 1
        da%data(da%ndata) = elem

    end subroutine dynamicIntArray_append

    function dynamicIntArray_get(da) result(data)
        class(dynamicIntArray), intent(in) :: da
        integer :: data(da%ndata)

        data(:) = da%data(:da%ndata)

    end function dynamicIntArray_get

    subroutine dynamicIntArray_deallocate(da)
        class(dynamicIntArray), intent(inout) :: da

        deallocate(da%data)

    end subroutine dynamicIntArray_deallocate


end module dynamicIntArrays