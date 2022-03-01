! Created by  on 23.06.20.

program nnTest
    use precision
    use neuralNetworks
    use random
    use kalmanFilter
    implicit none

    integer, parameter :: nInp = 1
    integer, parameter :: nOut = 1
    integer, parameter :: nLayers = 3
    integer, parameter :: nNodes(nLayers) = [3, 3, 1]
    character, parameter :: act(nLayers) = ['t', 't', 'l']
    real(dp), parameter :: lam0 = 0.9987

    integer, parameter :: nx = 40
    integer, parameter :: nPlot = 1000
    real(dp) :: x(nx), y(nx), yy(nx), err, l(nx), dL(nx,1), dIn(nx,1)
    real(dp) :: xPlot(nPlot), yPlot(nPlot), yyPlot(nPlot)
    integer :: i, j, epoch, order(nx)

    type(neuralNetwork) :: nn
    real(dp) :: inp(nInp), out1(nOut), out2(nOut), dOut(nInp), ones(nOut)
    real(dp) :: dx(nInp)
    real(dp), allocatable :: params(:), dParams(:), dxParams(:), P(:,:)
    real(dp) :: onesBatch(nOut, nx), dOutBatch(nOut, nx)
    real(dp) :: lam

    if (nNodes(nLayers) /= nOut) stop 'wrong nOut'


    ! initialize NN
    call nn_new(nn, nInp, nLayers, nNodes, act)

    print*, 'nParams=', nn%nParams

    ! allocate arrays
    allocate(params(nn%nParams))
    allocate(dParams(nn%nParams))
    allocate(dxParams(nn%nParams))

    ! some tests for the NN derivatives
    if (.false.) then
        call random_number(dx)
        dx = dx / sqrt(sum(dx**2)) * 1.e-4_dp

        ones = 1._dp

        call random_number(inp)

        call nn_forward(nn, inp, out1)

        ! dX / dOut = 1 (ones)
        call nn_backward(nn, ones, dOut)

        inp = inp + dx
        call nn_forward(nn, inp, out2)

        ! check nn derivative with finitie difference
        print*, sum(dOut*dx) / sqrt(sum(dx**2)), '?=', sum(out2 - out1) / sqrt(sum(dx**2))

        call random_number(dxParams)
        dxParams = dxParams / sqrt(sum(dxParams**2)) * 1.e-5_dp

        call nn_forward(nn, inp, out1)
        call nn_backward(nn, ones, dOut)
        call nn_getdParams(nn, dParams)

        call nn_getParams(nn, params)
        params = params + dxParams
        call nn_setParams(nn, params)

        call nn_forward(nn, inp, out2)

        ! check nn derivative of weights with finitie difference
        print*, sum(dParams*dxParams) / sqrt(sum(dParams**2)), '?=', sum(out2 - out1) / sqrt(sum(dParams**2))
        print*, ''

    end if



    ! set up some training data for the NN
    do i=1,nx
        x(i) = real(i, dp) / real(nx, dp) * 2._dp - 1.0_dp
        !x(i) = (real(i, dp)-0.5) / real(nx, dp) * 2._dp - 1.0_dp
        !y(i) = sin(50._dp * x(i)) * x(i)
        y(i) = exp(-abs(x(i)))
        !y(i) = 0; if (x(i) > 0) y(i) = 1._dp
        order(i) = i
    end do
    do i=1,nPlot
        xPlot(i) = real(i, dp) / real(nPlot, dp) * 2._dp - 1.0_dp
        yPlot(i) = exp(-abs(xPlot(i)))
        !yPlot(i) = 0; if (xPlot(i) > 0) yPlot(i) = 1._dp
        !y(i) = sin(50._dp * x(i)) * x(i)
    end do
    call nn_forwardBatch(nn, nx, x, yy)
    open(45, file='ref.txt')
    do j=1,nx
        write(45,*) x(j), y(j), yy(j)
    end do
    close(45)

    ! Kalman P matrix
    allocate(P(nn%nParams, nn%nParams))
    P = 0._dp
    do i=1,nn%nParams
        P(i,i) = 1.e5_dp
    end do

    ! kalman lambda
    lam = 0.99


    do epoch=1,20000
        call shuffle(nx, order)

        do i=1,nx
            inp(1) = x(order(i))
            call nn_forward(nn, inp, out1)
            err = (out1(1) - y(order(i)))!**2
            ! using ones to represent dErr / dout1
            ones(1) = 1._dp ! 2._dp * (out1(1) - y(order(i)))
            call nn_backward(nn, ones, dOut)

            ! using full batch in each update
            !call nn_forwardBatch(nn, nx, x, yy)
            !err = sum((yy - y)**2)
            !onesBatch(1,:) = 2._dp * (yy - y)
            !call nn_backwardBatch(nn, nx, onesBatch, dOutBatch)


            !yy(i) = err**2
            !call nn_forwardBatch(nn, nx, x, yy)
            !l = (yy - y)**2
            !print*, sum(l)
            !dL(:,1) = 2._dp * (yy - y)
            !call nn_backwardBatch(nn, nx, dL, dIn)
            call nn_getdParams(nn, dParams)
            call nn_getParams(nn, params)
            ! steepest descent
            !params = params - dParams * 1.e-3_dp ! * err
            ! kalman filter
            call updatekalman(lam0, lam, err, nn%nParams, params, dParams, P)
            call nn_setParams(nn, params)

        end do

        ! calculate error using the batch mode
        call nn_forwardBatch(nn, nx, x, yy)
        print*, sum((y - yy)**2)


        ! write the result to a file every once in a while for plotting
        if(modulo(epoch,10) == 0) then
            print*, 'prediction'
            call nn_forwardBatch(nn, nPlot, xPlot, yyPlot)
            open(45, file='out.txt')
            do j=1,nPlot
                write(45,*) xPlot(j), yPlot(j), yyPlot(j)
            end do
            close(45)
        end if

    end do



end program nnTest