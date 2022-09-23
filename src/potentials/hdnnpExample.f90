! Created by jonas on 01.03.22.

program hdnnpExample
    use atomicStructure
    use precision
    use hdnnps
    use constants
    implicit none

    type(hdnnp) :: hdnnpHandle
    type(atStruct) :: ats
    integer :: n, i
    real(dp) :: epot, ekin
    real(dp), allocatable :: forces(:,:), m(:,:), p(:,:), lastForces(:,:)
    real(dp) :: dEdLat(3,3), stress(3,3)

    real(dp), parameter :: stepsize = 0.5_dp * femtoseconds2atomic
    real(dp), parameter :: T = 140._dp ! temperature

    character(len=100) :: txt

    ! read the parameters for the hdnnp
    ! this is only a very simple example hdnnp. It will only work for a methylammonium molecule and is probably not physically correct.
    call parseRuNNerInput(hdnnpHandle, './hdnnp-params/')

    ! read geometry
    ! parsers for other file formats are also available (data, ascii, POSCAR)
    call as_readXYZ('input.xyz', ats)
    !call as_readAscii('input.ascii', ats)
    !call as_repeatUnitCell(ats, [3,3,3])

    ! allocate arrays
    allocate(forces(3, ats%nat))
    allocate(lastForces(3, ats%nat))
    allocate(m(3, ats%nat))
    allocate(p(3, ats%nat))

    ! use a 3*nat array for the mass to make calculations easier
    m(1,:) = getAtomicMass(ats%el) * massUnit2ElectronMass ! atomic units use electron mass = 1
    m(2,:) = m(1,:)
    m(3,:) = m(1,:)
    !m(:,:) = getAtomicMass(1) * massUnit2ElectronMass

    ! Boltzmann distributed velocities
    call randomPointSet(ats%nat, p)
    p = p * sqrt(T * kBoltzmann * m)

    ! remove mean velocity
    call center(ats%nat, p)

   ! ! Derivatives of the energy w.r.t the lattice vectors can also be computed ...
   ! call hdnnpEnergyAndForces(hdnnpHandle, ats, epot, forces, dEdLat)
   ! ! ... and transformed to the stress tensor. Of course this only works for periodic systems (not this example HDNNP).
   ! call dEdLat2Stress(ats%lat, dEdLat, stress)


    ! verlet integration
    call hdnnpEnergyAndForces(hdnnpHandle, ats, epot, forces)

    do i = 0, 1000
        ats%ats = ats%ats + p / m * stepSize + 0.5_dp * forces / m * stepSize**2
        lastForces = forces
        call hdnnpEnergyAndForces(hdnnpHandle, ats, epot, forces)
        p = p + stepSize / 2._dp * (forces + lastForces)

        ekin = sum(p**2 / m) / 2._dp
        print*, 'MD', i, i * stepsize / femtoseconds2atomic, epot, ekin, epot + ekin

        ! save trajectory
        write(txt, '(A,I5.5,A)') 'out/', i, '.xyz'
        call as_writeXYZ(txt, ats)
    end do

end program hdnnpExample