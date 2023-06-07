program main
    implicit none
    
    interface
        subroutine func(type, N, arr_ptr)
            integer, intent(in) :: type, N
            class(*), dimension(:), pointer :: arr_ptr
        end subroutine
    end interface
    
    integer, dimension(:), pointer :: int_ptr
    real(8), dimension(:), pointer :: real_ptr

    call func(1, 5, int_ptr)
    ! Use int_ptr as the allocated array

    call func(2, 10, real_ptr)
    ! Use real_ptr as the allocated array

    deallocate(int_ptr)
    deallocate(real_ptr)
end program

subroutine func(type, N, arr_ptr)
    integer, intent(in) :: type, N
    class(*), dimension(:), pointer :: arr_ptr
    integer, dimension(:), pointer :: arr_int
    real(8), dimension(:), pointer :: arr_real

    if (type == 1) then
        allocate(arr_int(N))
        arr_ptr => arr_int 
    else
        allocate(arr_real(N))
        arr_ptr => arr_real
    end if
end subroutine