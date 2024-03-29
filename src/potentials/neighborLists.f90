! Created by jonas on 1/30/22.

module neighborLists
    use precision
    use atomicStructure
    use pointsets
    use dynamicIntArrays
    implicit none

    type neighborList
        integer :: nNeis
        integer :: maxNeis
        real(dp), allocatable :: neis(:,:)
        integer, allocatable :: ineis(:)

    end type neighborList

contains

    ! brute force approach O(n^2)
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

    ! Puts atoms on grid O(1)
    subroutine buildNeighborListGrid(ats, rc, neiLists)
        type(atStruct), intent(in) :: ats
        real(dp), intent(in) :: rc
        type(neighborList), intent(out) :: neiLists(ats%nat)
        integer :: i, ncells(3), j, ix, iy, iz, sy, sz, jat
        integer :: jx, jy, jz, cx, cy, cz, jj
        real(dp) :: dlat(3), d(3), grid(3,3)
        integer :: ngrid(3), ddl
        type(dynamicIntArray), allocatable :: gridAts(:,:,:)
        type(atStruct) :: tmpats
        real(dp) :: relats(3, ats%nat)
        real(dp) :: latorigin(3)
        real(dp) :: targetCellSize

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

        tmpats = ats
        latorigin = 0._dp
        if (.not. ats%periodic) then ! make artificial lattice
            tmpats%lat = 0._dp
            do i=1,3
                latorigin(i) = minval(tmpats%ats(i,:)) - 1.e-4_dp
                tmpats%ats(i,:) = tmpats%ats(i,:) - latorigin(i)
                tmpats%lat(i,i) = maxval(tmpats%ats(i,:)) + 1.e-4_dp ! add a little to the box to prevent atoms from being directly on the border
            end do
        end if
        call getNcells(tmpats%lat, rc, ncells)

        ! this can probably be tuned a bit. In my tests, 0.6 gave good results
        targetCellSize = rc * 0.6_dp
        do i=1,3
            ngrid(i) = int(sqrt(sum(tmpats%lat(:,i)**2)) / targetCellSize) + 1
            grid(:,i) = tmpats%lat(:,i) / ngrid(i)
        end do

        allocate(gridAts(ngrid(1), ngrid(2), ngrid(3)))
        do ix=1,ngrid(1)
            do iy=1,ngrid(2)
                do iz=1,ngrid(3)
                    gridAts(ix, iy, iz) = newDynamicIntArray()
                end do
            end do
        end do
        call as_moveAtomsIntoCell(tmpats)

        ! Sort atoms into grid
        relats = tmpats%ats
        call toRelativeCoordinates(tmpats%nat, tmpats%lat, relats)
        do i=1,ats%nat
            ix = int(relats(1,i) * ngrid(1)) + 1
            iy = int(relats(2,i) * ngrid(2)) + 1
            iz = int(relats(3,i) * ngrid(3)) + 1
            call gridAts(ix, iy, iz)%append(i)
        end do

        ! how many neighboring grid cells do we need to consider?
        call getNcells(grid, rc, ncells)

        !$omp parallel do private(ix, iy, iz, cx, cy, cz, jx, jy, jz, ddl, dlat, j, jat, d)
        do i=1,tmpats%nat
            ! find atom in the grid
            ix = int(relats(1,i) * ngrid(1)) + 1
            iy = int(relats(2,i) * ngrid(2)) + 1
            iz = int(relats(3,i) * ngrid(3)) + 1
            ! loop over all needed neighbor cells
            loop_jx: do jx=-ncells(1),ncells(1)
                loop_jy: do jy=-ncells(2),ncells(2)
                    loop_jz: do jz=-ncells(3),ncells(3)
                        dlat = 0._dp
                        cx = ix + jx
                        cy = iy + jy
                        cz = iz + jz
                        ! wrap around, when periodic
                        if (cz < 1 .or. cz > ngrid(3)) then
                            if (.not. ats%periodic) then ! dont wrap with free b.c.
                                cycle loop_jz
                            end if
                            ddl = int(floor(1._dp * (cz-1) / ngrid(3)))
                            cz = cz - ngrid(3) * ddl
                            dlat = dlat + ats%lat(:,3) * ddl
                        end if
                        if (cy < 1 .or. cy > ngrid(2)) then
                            if (.not. ats%periodic) then
                                cycle loop_jy
                            end if
                            ddl = int(floor(1._dp * (cy-1) / ngrid(2)))
                            cy = cy - ngrid(2) * ddl
                            dlat = dlat + ats%lat(:,2) * ddl
                        end if
                        if (cx < 1 .or. cx > ngrid(1)) then
                            if (.not. ats%periodic) then
                                cycle loop_jx
                            end if
                            ddl = int(floor(1._dp * (cx-1) / ngrid(1)))
                            cx = cx - ngrid(1) * ddl
                            dlat = dlat + ats%lat(:,1) * ddl
                        end if
                        do j=1,gridAts(cx, cy, cz)%ndata ! now look at all the atoms in the cell
                            jat = gridAts(cx, cy, cz)%data(j)
                            if (i/=jat .or. (.not. (jx==0 .and. jy==0 .and. jz==0))) then ! exclude the atom itself (already in the list)
                                d(:) = tmpats%ats(:,i) - (tmpats%ats(:,jat) + dlat(:))
                                if (sum(d**2) < rc**2) then
                                    if (neiLists(i)%nNeis >= neiLists(i)%maxNeis) then
                                        call resizeNeiList(neiLists(i))
                                    end if
                                    neiLists(i)%nNeis = neiLists(i)%nNeis + 1
                                    neiLists(i)%neis(:,neiLists(i)%nNeis) = tmpats%ats(:,jat) + dlat(:) + latorigin(:)
                                    neiLists(i)%ineis(neiLists(i)%nNeis) = jat
                                end if
                            end if
                        end do
                    end do loop_jz
                end do loop_jy
            end do loop_jx
        end do
        !$omp end parallel do

    end subroutine buildNeighborListGrid
end module neighborLists