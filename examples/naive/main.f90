program program 
    ! use naive_mod

    ! class(data_type), pointer :: data
    ! type(data_real), target :: data_concrete
    
    ! call data_concrete%init(5,5,5,5, .False.)

    ! data => data_concrete
    ! call data%set_value(1,1,1,1, 5d0)

    ! print*, 'hello world', data%get_value(1,1,1,1)

    real(8), dimension(:), allocatable, target :: A
    real(8), dimension(:), pointer :: B
    
    
    allocate(A(-1:3))
    A = (/56,1,2,3,4/)

    B => A

    print*, A(0)
    print*, B(0)
end program 


