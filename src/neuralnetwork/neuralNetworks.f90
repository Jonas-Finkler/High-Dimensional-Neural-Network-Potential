

! simple, straight forward implementation of a neural network and backpropagation
! written by Jonas Finkler

! todo: use BLAS
module neuralNetworks
    use precision

    implicit none

    type nnLayer
        integer :: nInp, nOut
        real(dp), allocatable :: weights(:,:)
        real(dp), allocatable :: dWeights(:,:)
        real(dp), allocatable :: bias(:)
        real(dp), allocatable :: dBias(:)
        real(dp), allocatable :: c(:,:)
        real(dp), allocatable :: inp(:,:)
        integer :: batchsize
        character :: activation
    end type nnLayer

    type neuralNetwork
        integer :: nLayers, nInp, nOut
        type(nnLayer), allocatable :: layers(:)
        integer :: nParams
    end type neuralNetwork


contains

    function activate(n, in, act) result(out)
        integer, intent(in) :: n
        real(dp), intent(in) :: in(n)
        character, intent(in) :: act
        real(dp) :: out(n)
        integer :: i

        select case (act)
        case ('l') ! linear
            out = in
        case ('t') ! tanh
            out = tanh(in)
        case ('p') ! softplus
            out = log(1._dp + exp(in))
        case ('r') ! relu
            do i=1,n
                out(i) = max(0._dp, in(i))
            end do
        case default
            stop 'ERROR: Invalid activation function!'
        end select

    end function activate

    ! derivative of activation function
    function deactivate(n, in, act) result(out)
        integer, intent(in) :: n
        real(dp), intent(in) :: in(n)
        character, intent(in) :: act
        real(dp) :: out(n)
        integer :: i

        select case (act)
        case ('l') ! linear
            out = 1._dp
        case ('t') ! tanh
            out = 1._dp - tanh(in)**2
        case ('p') ! softplus
            out = 1._dp / (1._dp + exp(-1._dp * in))
        case ('r') ! relu
            do i=1,n
                if (in(i) > 0._dp) then
                    out(i) = 1._dp
                else
                    out(i) = 0._dp
                end if
            end do
        case default
            stop 'ERROR: Invalid activation function!'
        end select

    end function deactivate

    subroutine layer_zeroGrad(layer)
        type(nnLayer), intent(inout) :: layer

        if (allocated(layer%c)) deallocate(layer%c)
        if (allocated(layer%inp)) deallocate(layer%inp)

    end subroutine layer_zeroGrad

    subroutine layer_forward(layer, batchsize, inp, out)
        type(nnLayer), intent(inout) :: layer
        integer, intent(in) :: batchsize
        real(dp), intent(in) :: inp(layer%nInp, batchsize)
        real(dp), intent(out) :: out(layer%nOut, batchsize)
        integer :: b

        if (allocated(layer%inp) .or. allocated(layer%c)) stop 'ERROR: second forward pass'
        allocate(layer%inp(layer%nInp, batchsize))
        allocate(layer%c(layer%nOut, batchsize))

        layer%inp = inp
        layer%batchsize = batchsize

        do b=1,batchsize
            out(:,b) = matmul(layer%weights, inp(:,b)) + layer%bias
            layer%c(:,b) = out(:,b)
            out(:,b) = activate(layer%nOut, out(:,b), layer%activation)
        end do

    end subroutine layer_forward

    subroutine layer_backward(layer, batchsize, gradOut, gradInp, calcdOutdParams)
        type(nnLayer), intent(inout) :: layer
        integer, intent(in) :: batchsize
        real(dp), intent(in) :: gradOut(layer%nOut, batchsize)
        real(dp), intent(out) :: gradInp(layer%nInp, batchsize)
        logical, intent(in), optional :: calcdOutdParams
        real(dp) :: dAct(layer%nOut)
        integer :: i,j,b
        real(dp) :: tmp(layer%nOut)
        logical :: tmpCalcdOutdParams

        tmpCalcdOutdParams = .true. ! default
        if (present(calcdOutdParams)) then
            tmpCalcdOutdParams = calcdOutdParams
        end if

        if (batchsize /= layer%batchsize) stop 'ERROR: Wrong batchsize in backward pass'

        if (.not. (allocated(layer%c) .and. allocated(layer%inp))) stop 'ERROR: Backward pass before forward pass'

        gradInp = 0._dp
        layer%dWeights = 0._dp
        layer%dBias = 0._dp

        do b=1,batchsize
            dAct = deactivate(layer%nOut, layer%c(:, b), layer%activation)
            ! todo: this here could probably be a matrix multiplication
            tmp(:) = gradOut(:, b) * dAct(:)
            !gradInp(:,b) = matmul(tmp(:), layer%weights(:,:))
            do j=1,layer%nInp
                do i=1,layer%nOut
                    if (tmpCalcdOutdParams) then
                        layer%dWeights(i,j) = layer%dWeights(i,j) + tmp(i) * layer%inp(j,b)
                    end if
                    gradInp(j,b) = gradInp(j,b) + tmp(i) * layer%weights(i,j)
                end do
            end do
            if (tmpCalcdOutdParams) then
                layer%dBias(:) = layer%dBias + gradOut(:,b) * dAct(:)
            end if
        end do

    end subroutine layer_backward

    subroutine layer_new(layer, nInp, nNodes, act)
        type(nnLayer), intent(out) :: layer
        integer, intent(in) :: nInp
        integer, intent(in) :: nNodes
        character, intent(in) :: act

        layer%nInp = nInp
        layer%nOut = nNodes
        layer%activation = act
        allocate(layer%weights(nNodes, nInp))
        allocate(layer%bias(nNodes))
        allocate(layer%dWeights(nNodes, nInp))
        allocate(layer%dBias(nNodes))

    end subroutine

    subroutine nn_new(nn, nInp, nLayers, nNodes, act)
        type(neuralNetwork), intent(out) :: nn
        integer, intent(in) :: nInp
        integer, intent(in) :: nLayers
        integer, intent(in) :: nNodes(nLayers)
        character, intent(in) :: act(nLayers)
        integer :: i, tmp

        nn%nLayers = nLayers
        nn%nInp = nInp
        nn%nOut = nNodes(nLayers)

        allocate(nn%layers(nLayers))

        nn%nParams = 0

        tmp = nInp
        do i=1, nLayers
            call layer_new(nn%layers(i), tmp, nNodes(i), act(i))
            tmp = nNodes(i)
            nn%nParams = nn%nParams + nn%layers(i)%nOut * (nn%layers(i)%nInp + 1)
        end do

        ! todo: find better way to initialize weights
        do i=1, nLayers
            call random_number(nn%layers(i)%weights)
            call random_number(nn%layers(i)%bias)
            nn%layers(i)%weights = nn%layers(i)%weights - 0.5_dp
            nn%layers(i)%bias = nn%layers(i)%bias - 0.5_dp
        end do

    end subroutine nn_new

    subroutine nn_zeroGrad(nn)
        type(neuralNetwork), intent(inout) :: nn
        integer :: i

        do i=1,nn%nLayers
            call layer_zeroGrad(nn%layers(i))
        end do

    end subroutine nn_zeroGrad

    ! WARNING: If several forward passes are performed gradient information is only retained for the last one
    subroutine nn_forward(nn, inp, out)
        type(neuralNetwork), intent(inout) :: nn
        real(dp), intent(in) :: inp(nn%nInp)
        real(dp), intent(out) :: out(nn%nOut)
        real(dp) :: tmpInp(nn%nInp, 1), tmpOut(nn%nOut, 1)

        tmpInp(:,1) = inp

        call nn_forwardBatch(nn, 1, tmpInp, tmpOut)

        out = tmpOut(:,1)

    end subroutine nn_forward

    subroutine nn_forwardBatch(nn, batchsize, inp, out)
        type(neuralNetwork), intent(inout) :: nn
        integer, intent(in) :: batchsize
        real(dp), intent(in) :: inp(nn%nInp, batchsize)
        real(dp), intent(out) :: out(nn%nOut, batchsize)
        integer :: i
        real(dp), allocatable :: tmpInp(:,:), tmpOut(:,:)

        call nn_zeroGrad(nn)

        allocate(tmpOut(nn%nInp, batchsize))
        tmpOut = inp
        do i=1,nn%nLayers
            if (allocated(tmpInp)) then
                deallocate(tmpInp)
            end if
            allocate(tmpInp(nn%layers(i)%nInp, batchsize))
            tmpInp = tmpOut
            if (allocated(tmpOut)) then
                deallocate(tmpOut)
            end if
            allocate(tmpOut(nn%layers(i)%nOut, batchsize))
            call layer_forward(nn%layers(i), batchsize, tmpInp, tmpOut)
        end do
        out = tmpOut

        if (allocated(tmpOut)) then
            deallocate(tmpOut)
        end if
        if (allocated(tmpInp)) then
            deallocate(tmpInp)
        end if

    end subroutine nn_forwardBatch

    subroutine nn_backward(nn, gradOut, gradInp, calcdOutdParams)
        type(neuralNetwork), intent(inout) :: nn
        real(dp), intent(in) :: gradOut(nn%nOut) ! gradient w.r.t the nn output (dE / dOut_i)
        real(dp), intent(out) :: gradInp(nn%nInp) ! gradient w.r.t the nn input (dE / dInp_j = sum_i dOut_i / dInp_j * dE / dOut_i)
        logical, intent(in), optional :: calcdOutdParams

        call nn_backwardBatch(nn, 1, gradOut, gradInp, calcdOutdParams)

    end subroutine

    ! subroutine returns sum_i dOut_i / dInp_j * dE / dOut_i
    subroutine nn_backwardBatch(nn, batchsize, gradOut, gradInp, calcdOutdParams)
        type(neuralNetwork), intent(inout) :: nn
        integer, intent(in) :: batchsize
        real(dp), intent(in) :: gradOut(nn%nOut, batchsize) ! gradient w.r.t the nn output (dE / dOut_i)
        real(dp), intent(out) :: gradInp(nn%nInp, batchsize) ! gradient w.r.t the nn input (dE / dInp_j = sum_i dOut_i / dInp_j * dE / dOut_i)
        logical, intent(in), optional :: calcdOutdParams
        integer :: i
        real(dp), allocatable :: tmpInp(:,:), tmpOut(:,:)


        allocate(tmpInp(nn%nOut, batchsize))
        tmpInp = gradOut
        do i=nn%nLayers,1,-1
            if (allocated(tmpOut)) then
                deallocate(tmpOut)
            end if
            allocate(tmpOut(nn%layers(i)%nOut, batchsize))
            tmpOut = tmpInp
            if (allocated(tmpInp)) then
                deallocate(tmpInp)
            end if
            allocate(tmpInp(nn%layers(i)%nInp, batchsize))
            call layer_backward(nn%layers(i), batchsize, tmpOut, tmpInp, calcdOutdParams)
        end do
        gradInp = tmpInp

        if (allocated(tmpOut)) then
            deallocate(tmpOut)
        end if
        if (allocated(tmpInp)) then
            deallocate(tmpInp)
        end if

    end subroutine nn_backwardBatch


    subroutine nn_getParams(nn, params)
        type(neuralNetwork), intent(inout) :: nn
        real(dp), intent(out) :: params(nn%nParams)
        integer :: i, j
        integer :: c

        c = 1
        do i=1,nn%nLayers
            do j=1,nn%layers(i)%nInp
                params(c:c+nn%layers(i)%nOut-1) = nn%layers(i)%weights(:,j)
                c = c + nn%layers(i)%nOut
            end do
            params(c:c+nn%layers(i)%nOut-1) = nn%layers(i)%bias(:)
            c = c + nn%layers(i)%nOut
        end do

        if (c /= nn%nParams + 1) stop 'ERROR: nParams is wrong!'

    end subroutine nn_getParams

    subroutine nn_setParams(nn, params)
        type(neuralNetwork), intent(inout) :: nn
        real(dp), intent(in) :: params(nn%nParams)
        integer :: i, j
        integer :: c

        c = 1
        do i=1,nn%nLayers
            do j=1,nn%layers(i)%nInp
                nn%layers(i)%weights(:,j) = params(c:c+nn%layers(i)%nOut-1)
                c = c + nn%layers(i)%nOut
            end do
            nn%layers(i)%bias(:) = params(c:c+nn%layers(i)%nOut-1)
            c = c + nn%layers(i)%nOut
        end do

        if (c /= nn%nParams + 1) stop 'ERROR: nParams is wrong!'

    end subroutine nn_setParams

    subroutine nn_getdParams(nn, params)
        type(neuralNetwork), intent(inout) :: nn
        real(dp), intent(out) :: params(nn%nParams)
        integer :: i, j
        integer :: c

        c = 1
        do i=1,nn%nLayers
            do j=1,nn%layers(i)%nInp
                params(c:c+nn%layers(i)%nOut-1) = nn%layers(i)%dWeights(:,j)
                c = c + nn%layers(i)%nOut
            end do
            params(c:c+nn%layers(i)%nOut-1) = nn%layers(i)%dBias(:)
            c = c + nn%layers(i)%nOut
        end do

        if (c /= nn%nParams + 1) stop 'ERROR: nParams is wrong!'

    end subroutine nn_getdParams

end module neuralNetworks