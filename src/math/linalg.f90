! Written by JAF

module linalg
    use precision
    use constants

    implicit none

contains


    ! modified (numerically stable) Gram-Schmidt process
    subroutine gramSchmidt(d, n, nu, v, u)
        implicit none
        integer, intent(in) :: d ! dimensionality of the vectors
        integer, intent(in) :: n ! number of vectors that are created
        integer, intent(in) :: nu  ! number of vectors already orthogonal and in the first nu columns of u
        real(dp), intent(in) :: v(d, n) ! linearly independent vectors
        real(dp), intent(inout) :: u(d, n) ! the remaining n-nu cols will be filled with orthonormal vectors orthogonal to the first nu ones alredy present

        integer :: i, j

        do i=nu+1,n
          u(:,i) = v(:,i)
          do j=1,i-1!i-1,1,-1
            u(:,i) = u(:,i) - u(:,j) * sum(u(:,i) * u(:,j)) / sum(u(:,j)**2)
          end do
          u(:,i) = u(:,i) / sqrt(sum(u(:,i)**2))
        end do

    end subroutine

    ! evecs(:,i) is eigenvector to eval(i), only needs lower half of matrix
    subroutine eigensystemDiag(n, M, evals, evecs)
        integer, intent(in) :: n
        real(dp), intent(in) :: M(n, n)
        real(dp), intent(out) :: evals(n), evecs(n, n)
        real(dp) :: A(n, n), W(n), LWORKQUERY(1)
        real(dp), allocatable :: WORK(:)
        integer :: LWORK, INFO, LDA

        A = M
        LDA = n

        !DSYEV does not seperate rotation from translation in the first 6 eigenvectors of the hessian, as DGEEV would. But we dont need this anyways.
        call DSYEV('V', 'L', n, A, LDA, W, LWORKQUERY, -1, INFO)

        LWORK = int(LWORKQUERY(1))
        allocate(WORK(LWORK))
        call DSYEV('V', 'L', n, A, LDA, W, WORK, LWORK, INFO)
        deallocate(WORK)
        if(INFO/=0) stop "DSYEV failed"
        evecs = A
        evals = W
    end subroutine eigensystemDiag

    subroutine eigenvaluesDiag(n, m, evals)
        integer, intent(in) :: n
        real(dp), intent(in) :: m(n, n)
        real(dp), intent(out) :: evals(n)
        real(dp) :: LWORKQUERY(1)
        real(dp), allocatable :: A(:,:), W(:)
        real(dp), allocatable :: WORK(:)
        integer :: LWORK, INFO, LDA

        allocate(A(n,n))
        allocate(W(n))

        A = m
        LDA = n

        !DSYEV does not seperate rotation from translation in the first 6 eigenvectors of the hessian, as DGEEV would. But we dont need this anyways.
        call DSYEV('N', 'L', n, A, LDA, W, LWORKQUERY, -1, INFO)

        LWORK = int(LWORKQUERY(1))
        allocate(WORK(LWORK))
        call DSYEV('N', 'L', n, A, LDA, W, WORK, LWORK, INFO)
        deallocate(WORK)
        if(INFO/=0) stop "DSYEV failed"
        evals = W
    end subroutine eigenvaluesDiag

    subroutine eigenvaluesComplex(n, m, evals)
        integer, intent(in) :: n
        complex(qp), intent(in) :: m(n, n)
        complex(qp), intent(out) :: evals(n)
        complex(qp) :: LWORKQUERY(1)
        complex(qp), allocatable :: A(:,:), W(:)
        complex(qp), allocatable :: WORK(:)
        integer :: LWORK, INFO, LDA
        complex(qp) :: evecsL(1,1), evecsR(1,1)
        real(dp), allocatable :: RWORK(:)

        allocate(RWORK(2*n))
        allocate(A(n,n))
        allocate(W(n))

        A = m
        LDA = n

        call ZGEEV('N', 'N', n, A, n, evals, evecsL, 1, evecsR, 1, LWORKQUERY, -1, RWORK, INFO)

        LWORK = int(LWORKQUERY(1))
        allocate(WORK(LWORK))
        call ZGEEV('N', 'N', n, A, n, evals, evecsL, 1, evecsR, 1, WORK, LWORK, RWORK, INFO)
        deallocate(WORK)
        if(INFO/=0) stop "ZGEEV failed"
    end subroutine eigenvaluesComplex

    subroutine solveEqSystems(neq, n, A, x)
        implicit none
        integer, intent(in) :: neq, n
        real(dp), intent(in) :: A(n,n)
        real(dp), intent(inout) :: x(n, neq)
        integer :: IPIV(n), info
        real(dp), allocatable :: M(:,:)

        allocate(M(n,n))
        M(:,:) = A(:,:)

        ! solve equations
        call dgesv(n, neq, M, n, IPIV, x, n, info)
        if (info /= 0) stop 'Error: could not solve linear equation system'
    end subroutine

    subroutine solveSymEqSystems(neq, n, A, x)
        implicit none
        integer, intent(in) :: neq, n
        real(dp), intent(in) :: A(n,n)
        real(dp), intent(inout) :: x(n, neq)
        integer :: IPIV(n), info, lwork
        real(dp), allocatable :: work(:)
        real(dp) :: nwork(1)
        real(dp), allocatable :: M(:,:)

        allocate(M(n,n))
        M(:,:) = A(:,:)

        ! this will get us nwork, the optimal size of the working array
        call dsysv('L', n, neq, M, n, IPIV, x, n, nwork,   -1, info)
        ! allocate working memory
        LWORK = int(nwork(1))
        allocate(WORK(LWork))
        ! solve equations
        call dsysv('L', n, neq, M, n, IPIV, x, n, WORK, LWORK, info)
        deallocate(WORK)
        if (info /= 0) stop 'Error: could not solve linear equation system'
    end subroutine

    subroutine solveSymEqSystemConjGrad(n, A, b, x)
        implicit none
        real(8), parameter :: tol = 1.d-30
        integer, intent(in) :: n
        real(8), intent(in) :: A(n,n), b(n)
        real(8), intent(inout) :: x(n) ! initial guess in, result out.

        real(8) :: r(n), p(n), alpha, Ap(n), rr, beta
        integer :: i

        i = 0
        r(:) = b(:) - matmul(A, x)
        p(:) = r(:)
        rr = sum(r**2)
        do while (rr > tol)
            !Ap(:) = matmul(A, p)
            call DSYMV('L', n, 1.d0, A, n, p, 1, 0.d0, Ap, 1)
            alpha = rr / sum(p * Ap)
            x(:) = x(:) + alpha * p(:)
            r(:) = r(:) - alpha * Ap(:)
            beta = sum(r**2) / rr
            p(:) = r(:) + beta * p(:)
            rr = sum(r**2)
            i = i + 1
        end do
        !print*, i

    end subroutine

    function vecMulVecT(n, a, b) result(c)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n), b(n)
        real(dp) :: c(n,n)
        integer :: i, j

        do i=1,n
            do j=1,n
                c(i,j) = a(i) * b(j)
            end do
        end do


    end function

    function matMulVec3(mat, vec) result(res)
        implicit none
        real(dp), dimension(3, 3), intent(in) :: mat
        real(dp), dimension(3), intent(in) :: vec
        real(dp), dimension(3) :: res

        res(1) = sum(mat(1, :) * vec(:))
        res(2) = sum(mat(2, :) * vec(:))
        res(3) = sum(mat(3, :) * vec(:))
    end function

    function matMulVecNM(n, m, mat, vec) result(res)
        integer, intent(in) :: n, m
        real(dp), intent(in) :: mat(n,m), vec(m)
        real(dp) :: res(n)
        integer :: i

        do i=1,n
            res(i) = sum(mat(i,:) * vec(:))
        end do

    end function

    function cross(a, b) result(c)
        implicit none
        real(dp), dimension(3), intent(in) :: a, b
        real(dp), dimension(3) :: c

        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function

    function det3D(M) result(det)
        real(dp), intent(in) :: M(3,3)
        real(dp) :: det

        det = M(1,1) * M(2,2) * M(3,3) &
            + M(1,2) * M(2,3) * M(3,1) &
            + M(1,3) * M(2,1) * M(3,2) &
            - M(1,1) * M(2,3) * M(3,2) &
            - M(1,2) * M(2,1) * M(3,3) &
            - M(1,3) * M(2,2) * M(3,1)

    end function det3D

    subroutine inv3D(a1, a2, a3, b1, b2, b3)
        implicit none

        real(dp), intent(in), dimension(3) :: a1, a2, a3
        real(dp), intent(out), dimension(3) :: b1, b2, b3
        real(dp) :: det


        det = a1(1) * a2(2) * a3(3) &
            + a1(2) * a2(3) * a3(1) &
            + a1(3) * a2(1) * a3(2) &
            - a1(1) * a2(3) * a3(2) &
            - a1(2) * a2(1) * a3(3) &
            - a1(3) * a2(2) * a3(1)

        b1 = [a2(2) * a3(3) - a2(3) * a3(2), a1(3) * a3(2) - a1(2) * a3(3), a1(2) * a2(3) - a1(3) * a2(2)]
        b2 = [a2(3) * a3(1) - a2(1) * a3(3), a1(1) * a3(3) - a1(3) * a3(1), a1(3) * a2(1) - a1(1) * a2(3)]
        b3 = [a2(1) * a3(2) - a2(2) * a3(1), a1(2) * a3(1) - a1(1) * a3(2), a1(1) * a2(2) - a1(2) * a2(1)]

        b1 = b1 / det
        b2 = b2 / det
        b3 = b3 / det
    end subroutine

    subroutine inv3DM(A, B)
        implicit none
        real(dp), intent(in), dimension(3, 3) :: A
        real(dp), intent(out), dimension(3, 3) :: B
        real(dp), dimension(3) :: a1, a2, a3
        real(dp), dimension(3) :: b1, b2, b3
        real(dp) :: det

        a1 = A(1, :)
        a2 = A(2, :)
        a3 = A(3, :)

        det = det3D(A)

        b1 = [a2(2) * a3(3) - a2(3) * a3(2), a1(3) * a3(2) - a1(2) * a3(3), a1(2) * a2(3) - a1(3) * a2(2)]
        b2 = [a2(3) * a3(1) - a2(1) * a3(3), a1(1) * a3(3) - a1(3) * a3(1), a1(3) * a2(1) - a1(1) * a2(3)]
        b3 = [a2(1) * a3(2) - a2(2) * a3(1), a1(2) * a3(1) - a1(1) * a3(2), a1(1) * a2(2) - a1(2) * a2(1)]

        b1 = b1 / det
        b2 = b2 / det
        b3 = b3 / det

        B(1, :) = b1
        B(2, :) = b2
        B(3, :) = b3
    end subroutine


    ! only for symmetric matrices so far
    subroutine qrDecomposition(n, A, Q, R)
        integer, intent(in) :: n
        real(dp), intent(in) :: A(n,n)
        real(dp), intent(out) :: Q(n,n), R(n,n)

        real(dp) :: tau(n), workq(1)
        real(dp), allocatable :: work(:)
        integer :: lwork, info, i

        R = A
        call dgeqrf(n, n, R, n, tau, workq, -1, info)
        lwork = int(workq(1))
        allocate(work(int(lwork)))
        call dgeqrf(n, n, R, n, tau, work, lwork, info)
        if (info /= 0) stop 'QR decomposition failed!'
        deallocate(work)

        Q = 0._dp
        do i=1,n
            Q(i,i) = 1._dp
        end do

        call dormqr('L', 'N', n, n, n, R, n, tau, Q, n, workq, -1, info)
        lwork = int(workq(1))
        allocate(work(int(lwork)))
        call dormqr('L', 'N', n, n, n, R, n, tau, Q, n, work, lwork, info)
        if (info /= 0) stop 'QR decomposition failed! (2)'
        deallocate(work)

        ! lapack stores information about Q in the lower half of R. Removing it.
        do i=1,n-1
            R(i+1:n,i) = 0._dp
        end do

    end subroutine QRdecomposition

    function dot(a,b)
        real(dp), intent(in) :: a(:), b(:)
        real(dp) :: dot
        if (size(a) /= size(b)) stop 'dot product uf array with not same length'
        dot = sum(a(:) *  b(:))
    end function

    subroutine recLattice(lat, reclat)
        real(dp), intent(in) :: lat(3,3)
        real(dp), intent(out) :: reclat(3,3)

        call inv3DM(transpose(lat), reclat)
        reclat(:,:) = reclat(:,:) * (2._dp * PI )

    end subroutine recLattice

    function normalizedDot(a, b) result(t)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: t

        t = sum(a*b) / sqrt(sum(a**2) * sum(b**2))

    end function normalizedDot

    function vecAngle(a, b) result(t)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: t

        t = acos(sum(a*b) / sqrt(sum(a**2) * sum(b**2)))

    end function vecAngle

    function vlen(a) result(d)
        real(dp), intent(in) :: a(3)
        real(dp) :: d
        d = sqrt(sum(a**2))
    end function vlen

end module linalg