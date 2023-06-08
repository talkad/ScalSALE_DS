program program 
    use naive_mod

    class(data_type), pointer :: data
    type(data_real), target :: data_concrete
    
    call data_concrete%init(5,5,5,5, .False.)

    data => data_concrete
    call data%set_value(1,1,1,1, 5d0)

    print*, 'hello world', data%get_value(1,1,1,1)
end program 


