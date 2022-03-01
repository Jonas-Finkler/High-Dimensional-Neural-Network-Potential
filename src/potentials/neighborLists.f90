! Created by jonas on 1/30/22.

module neighborLists
    use precision
    use atomicStructure
    use pointsets
    implicit none

    type neighborList
        integer :: nNeis
        integer :: maxNeis
        real(dp), allocatable :: neis(:,:)
        integer, allocatable :: ineis(:)

    end type neighborList

contains

    subroutine buildNeighborList(ats, rc, neiLists)
        type(atStruct), intent(in) :: ats
        real(dp), intent(in) :: rc
        type(neighborList), intent(out) :: neiLists(ats%nat)
        integer :: i, ncells(3), j, ix, iy, iz, sy, sz, jstart
        real(dp) :: dlat(3), d(3)

        do i=1,ats%nat
            neiLists(i)%nNeis = 1
            neiLists(i)%maxNeis = 16
            allocate(neiLists(i)%neis(3,neiLists(i)%maxNeis))
            allocate(neiLists(i)%ineis(neiLists(i)%maxNeis))
            ! every atom is its own first neighbor
            neiLists(i)%ineis = -1
            neiLists(i)%neis(:, 1) = ats%ats(:,i)
            neiLists(i)%ineis(1) = i
        end do

        if (ats%periodic) then
            call getNcells(ats%lat, rc, ncells)
        else
            ncells = 0
        end if

        do i=1,ats%nat
            do ix=0,ncells(1)
                sy = -ncells(2)
                if (ix==0) sy = 0
                do iy=sy,ncells(2)
                    sz = -ncells(3)
                    if (ix==0.and.iy==0) sz = 0
                    do iz=sz,ncells(3)
                        dlat(:) = ix * ats%lat(:,1) + iy * ats%lat(:,2) + iz * ats%lat(:,3)
                        !dlat(:) = ix * ats%lat(1,:) + iy * ats%lat(2,:) + iz * ats%lat(3,:)
                        jstart = 1
                        if (ix==0 .and. iy==0 .and. iz==0) jstart=i+1
                        do j=jstart,ats%nat
                            !if (ix/=0 .or. iy/=0 .or. iz/=0 .or. i/=j) then
                            d(:) = ats%ats(:,i) - (ats%ats(:,j) + dlat(:))
                            if (sum(d**2) < rc**2) then
                                if (neiLists(i)%nNeis >= neiLists(i)%maxNeis) then
                                    call resizeNeiList(neiLists(i))
                                end if
                                neiLists(i)%nNeis = neiLists(i)%nNeis + 1
                                neiLists(i)%neis(:,neiLists(i)%nNeis) = ats%ats(:,j) + dlat(:)
                                neiLists(i)%ineis(neiLists(i)%nNeis) = j

                                if (neiLists(j)%nNeis >= neiLists(j)%maxNeis) then
                                    call resizeNeiList(neiLists(j))
                                end if
                                neiLists(j)%nNeis = neiLists(j)%nNeis + 1
                                neiLists(j)%neis(:,neiLists(j)%nNeis) = ats%ats(:,i) - dlat(:)
                                neiLists(j)%ineis(neiLists(j)%nNeis) = i
                            end if
                            !end if
                        end do
                    end do
                end do
            end do
        end do



    end subroutine buildNeighborList

    subroutine resizeNeiList(neiList)
        type(neighborList), intent(inout) :: neiList
        real(dp), allocatable :: tmpNeis(:,:)
        integer, allocatable :: tmpiNeis(:)

        allocate(tmpNeis(3,neiList%nNeis))
        allocate(tmpiNeis(neiList%nNeis))
        tmpNeis(:,:) = neiList%neis(:,:)
        tmpiNeis(:) = neiList%ineis(:)
        deallocate(neiList%neis)
        deallocate(neiList%ineis)
        allocate(neiList%neis(3,neiList%nNeis * 2))
        allocate(neiList%ineis(neiList%nNeis * 2))
        neiList%neis(:,:neiList%nNeis) = tmpNeis(:,:neiList%nNeis)
        neiList%ineis(:neiList%nNeis) = tmpiNeis(:neiList%nNeis)
        neiList%maxNeis = neiList%maxNeis * 2


    end subroutine resizeNeiList

end module neighborLists