! Created by jonas on 01.03.22.

program hdnnpExample
    use atomicStructure
    use precision
    use hdnnps
    use constants
    use globalHdnnps
    implicit none

    type(atStruct) :: ats
    integer :: n, i
    real(dp) :: epot, ekin
    real(dp), allocatable :: forces(:,:), m(:,:), p(:,:), lastForces(:,:)
    real(dp) :: dEdLat(3,3), stress(3,3), energyBias

    real(dp), parameter :: stepsize = 0.8_dp * femtoseconds2atomic
    real(dp), parameter :: T = 300._dp ! temperature
    real(dp), parameter :: bias = 2._dp

    character(len=100) :: txt


    ! read the parameters for the hdnnps (you can get them here: https://github.com/Jonas-Finkler/MaPbI3-HDNNP)
    ! the file nnps.txt should contain the number of potentials in the first line and then the directories with the parameters on the following lines:
    ! 5
    ! .../MaPbI3-HDNNP/nnp-parameters/MaPbI3-1/
    ! .../MaPbI3-HDNNP/nnp-parameters/MaPbI3-2/
    ! ...
    ! .../MaPbI3-HDNNP/nnp-parameters/MaPbI3-5/
    call initGlobalHdnnpsFromInputFile('nnps.txt')

    ! read geometry
    ! parsers for other file formats are also available (data, xyz, POSCAR)
    call as_readAscii('in.ascii', ats, .true.) ! true, when input file is in angstrom units

    ! make the unit cell larger
    call as_repeatUnitCell(ats, [2, 2, 2])

    ! allocate arrays
    allocate(forces(3, ats%nat))
    allocate(lastForces(3, ats%nat))
    allocate(m(3, ats%nat))
    allocate(p(3, ats%nat))

    ! use a 3*nat array for the mass to make calculations easier
    m(1,:) = getAtomicMass(ats%el) * massUnit2ElectronMass ! atomic units use electron mass = 1
    m(2,:) = m(1,:)
    m(3,:) = m(1,:)


    ! Boltzmann distributed velocities
    call randomPointSet(ats%nat, p)
    p = p * sqrt(T * kBoltzmann * m)

    ! remove mean velocity
    call center(ats%nat, p)


    ! get the energy (mean of all predictions plus a bias proportional to the variance)
    ! The potential should be stabe for simulations up to 400K.
    ! But, the bias prevents the MD from going to regions where training data is missing in higher T simulations
    ! or when large MC moves are used that probe untrained configuration space.
    call globalHdnnpEnergyAndForcesWithBias(ats, bias, epot, forces, dEdLat, energyBias)

    ! verlet integration
    do i = 0, 1000
        ats%ats = ats%ats + p / m * stepSize + 0.5_dp * forces / m * stepSize**2
        lastForces = forces
        call globalHdnnpEnergyAndForcesWithBias(ats, bias, epot, forces, dEdLat, energyBias)
        p = p + stepSize / 2._dp * (forces + lastForces)

        ekin = sum(p**2 / m) / 2._dp
        print*, 'MD', i, i * stepsize / femtoseconds2atomic, epot, ekin, epot + ekin, energyBias

        ! save trajectory
        write(txt, '(A,I5.5,A)') 'out/', i, '.ascii'
        call as_writeAscii(txt, ats, .true.)
    end do

end program hdnnpExample
