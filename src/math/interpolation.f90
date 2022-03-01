! Created by jonas on 4/26/21.

module interpolation
    use precision
    implicit none

contains

    function linearInterpol3D(data, pos) result(res)
        real(dp), intent(in) :: data(2,2,2)
        real(dp), intent(in) :: pos(3)
        real(dp) :: res
        real(dp) :: t1(2,2), t2(2)

        t1(:,:) = pos(1) * data(2,:,:) + (1._dp - pos(1)) * data(1,:,:)
        t2(:) = pos(2) * t1(2,:) + (1._dp - pos(2)) * t1(1,:)
        res = pos(3) * t2(2) + (1._dp - pos(3)) * t2(1)

    end function linearInterpol3D


end module interpolation