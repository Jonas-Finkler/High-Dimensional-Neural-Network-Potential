! Created by jonas on 08.09.20.

module parameters
    use precision
    implicit none

    ! todo: these parameters should maybe be read from an input file?

    real(dp) :: ewaldSummationPrecision = 1.e-7 ! precision of the Ewald summation
    real(dp) :: waterModelCutoff = 10._dp ! cutoff radius of the ewald terms for periodic water


end module parameters