! Created by  on 29.01.20.

module dynamicArrays

    use precision
    implicit none

    type dynamicArray
        integer :: dim
        integer :: ndata
        integer :: size
        real(dp), allocatable :: data(:,:)
    contains
        procedure :: append => dynamicArray_append
        procedure :: deallocate => dynamicArray_deallocate
        procedure :: get => dynamicArray_get
    end type dynamicArray

    type dcaContainer
        class(*), allocatable :: val
    end type dcaContainer

    type dynamicClassArray
        integer :: ndata
        integer :: size
        type(dcaContainer), allocatable :: data(:)
    contains
        procedure :: append => dynamicClassArray_append
        procedure :: deallocate => dynamicClassArray_deallocate
       ! procedure :: get => dynamicClassArray_get
    end type dynamicClassArray


contains

    function newDynamicArray(dim) result(da)
        integer, intent(in) :: dim
        type(dynamicArray) :: da

        da%dim = dim
        da%ndata = 0
        da%size = 1
        allocate(da%data(da%dim, da%size))

    end function newDynamicArray

    subroutine dynamicArray_append(da, elem)
        class(dynamicArray), intent(inout) :: da
        real(dp), intent(in) :: elem(da%dim)
        real(dp), allocatable :: tmp(:,:)

        if (da%ndata + 1 > da%size) then
            allocate(tmp(da%dim, da%ndata))
            tmp(:,:) = da%data(:, :da%ndata)
            deallocate(da%data)
            da%size = da%size * 2
            allocate(da%data(da%dim, da%size))
            da%data(:,:da%ndata) = tmp(:,:)
            deallocate(tmp)
        end if

        da%ndata = da%ndata + 1
        da%data(:, da%ndata) = elem

    end subroutine dynamicArray_append

    function dynamicArray_get(da) result(data)
        class(dynamicArray), intent(in) :: da
        real(dp) :: data(da%dim, da%ndata)

        data(:,:) = da%data(:,:da%ndata)

    end function dynamicArray_get

    subroutine dynamicArray_deallocate(da)
        class(dynamicArray), intent(inout) :: da

        deallocate(da%data)

    end subroutine dynamicArray_deallocate

    function newDynamicClassArray() result(da)
        type(dynamicClassArray) :: da

        da%ndata = 0
        da%size = 1
        allocate(da%data(da%size))

    end function newDynamicClassArray

    subroutine dynamicClassArray_append(da, elem)
        class(dynamicClassArray), intent(inout) :: da
        !type(dcaContainer), intent(in) :: elem
        class(*), intent(in) :: elem
        type(dcaContainer), allocatable :: tmp(:)
        type(dcaContainer) :: cont

        cont%val = elem

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
        da%data(da%ndata) = cont

    end subroutine dynamicClassArray_append

    !todo: find a way
   ! subroutine dynamicClassArray_get(da, data)
   !     class(dynamicClassArray), intent(in) :: da
   !     class(*), intent(out) :: data(da%ndata)
   !     integer :: i


   !     do i=1, da%ndata
   !         data(i) = da%data(i)%val

   !     end do

   ! end subroutine dynamicClassArray_get

    subroutine dynamicClassArray_deallocate(da)
        class(dynamicClassArray), intent(inout) :: da

        deallocate(da%data)

    end subroutine dynamicClassArray_deallocate

end module dynamicArrays