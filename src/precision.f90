module precision
    use iso_c_binding
    integer, parameter :: dp=c_double
    integer, parameter :: qp=kind((1._dp, 1._dp))
end module precision
