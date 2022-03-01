! Created by  on 5/14/20.

module doubleLinkedList
    use precision
    implicit none

    type dll
        integer :: len
        type(dllNode), pointer :: node => null()
    end type dll

    type dllNode
        integer :: x ! content of list
        type(dllNode), pointer :: next => null()
        type(dllNode), pointer :: prev => null()
    end type dllNode

    type dllNodePointer
        type(dllNode), pointer :: p
    end type dllNodePointer

contains

    ! initialize a new list
    subroutine dll_new(list)
        type(dll), intent(out) :: list
        list%len = 0
    end subroutine

    ! insert element into list
    ! will bi inserted at the end (list%node%prev -> new element)
    subroutine dll_insert(list, x)
        type(dll), intent(inout) :: list
        integer, intent(in) :: x
        type(dllNode), pointer :: node => null()

        allocate(node)
        node%x = x

        if (list%len == 0) then
            list%node => node
            node%prev => node
            node%next => node
        else
            node%next => list%node
            node%prev => list%node%prev
            list%node%prev%next => node
            list%node%prev => node
        end if

        list%len = list%len + 1

    end subroutine

    ! cobines two lists by inserting other at the end of list
    ! other is empty after the operation
    subroutine dll_join(list, other)
        type(dll), intent(inout) :: list
        type(dll), intent(inout) :: other

        type(dllNode), pointer :: tmp

        if (other%len == 0) return

        if (list%len == 0) then
            list%len = other%len
            list%node => other%node
            return
        end if

        tmp => other%node%prev

        list%node%prev%next => other%node
        other%node%prev =>list%node%prev

        tmp%next => list%node
        list%node%prev => tmp

        list%len = list%len + other%len

        other%node => null()
        other%len = 0

    end subroutine

    ! removes a node from the list
    ! the node is deallocated during the operation
    ! caution: this does not check, if the node is actually part of the list
    subroutine dll_remove(list, node)
        type(dll), intent(inout) :: list
        type(dllNode), pointer, intent(in) :: node
        type(dllNode), pointer :: tmp

        tmp => node
        if (associated(list%node, tmp)) then
            list%node => list%node%next
        end if

        tmp%next%prev => tmp%prev
        tmp%prev%next => tmp%next

        list%len = list%len - 1

        deallocate(tmp)

    end subroutine

    ! print the list and check if the list is consistent
    subroutine dll_print(list)
        type(dll), intent(in) :: list
        integer :: i
        type(dllNode), pointer :: node

        if (list%len == 0) return
        node => list%node
        do i=1,list%len
            print*, node%x, node%next%prev%x, node%prev%next%x
            if (.not. associated(node%next%prev, node)) stop 'dll inconsistency'
            if (.not. associated(node%prev%next, node)) stop 'dll inconsistency'
            node => node%next
        end do

    end subroutine

    ! deallocate all members of the list
    subroutine dll_deallocate(list)
        type(dll), intent(inout) :: list
        type(dllNode), pointer :: node, delNode

        integer :: i
        node => list%node

        do i=1,list%len
            delNode => node
            node => node%next
            deallocate(delNode)
        end do

        list%node => null()
        list%len = 0

    end subroutine dll_deallocate

end module doubleLinkedList