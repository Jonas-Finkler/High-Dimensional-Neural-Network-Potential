module inertia
    use precision
    implicit none
contains

    subroutine inertiaTensor(nat, ats, I)
        integer, intent(in) :: nat
        real(dp), intent(in), dimension(3,nat) :: ats
        real(dp), intent(out), dimension(3,3) :: I
        real(dp) :: ats2sum(3)

        ats2sum = sum(ats**2,2)

        I(1,1) = ats2sum(2) + ats2sum(3)
        I(2,2) = ats2sum(1) + ats2sum(3)
        I(3,3) = ats2sum(1) + ats2sum(2)

        I(1,2) = - sum(ats(1,:)*ats(2,:))
        I(1,3) = - sum(ats(1,:)*ats(3,:))
        I(2,3) = - sum(ats(2,:)*ats(3,:))

        I(2,1) = I(1,2)
        I(3,1) = I(1,3)
        I(3,2) = I(2,3)


    end subroutine inertiaTensor

    function inertiaDeterminant(nat, ats)
        real(dp) :: inertiaDeterminant
        integer, intent(in) :: nat
        real(dp), intent(in), dimension(3,nat) :: ats
        real(dp) :: I(3,3)

        call inertiaTensor(nat, ats,I)

        inertiaDeterminant = &
                + I(1,1) * (I(2,2)*I(3,3) - I(2,3)*I(3,2))&
                - I(1,2) * (I(2,1)*I(3,3) - I(2,3)*I(3,1))&
                + I(1,3) * (I(2,1)*I(3,2) - I(2,2)*I(3,1))

    end function inertiaDeterminant

end module inertia
