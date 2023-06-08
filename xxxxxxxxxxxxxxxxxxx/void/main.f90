! program main
!     implicit none
    
!     interface
!         subroutine func(type, N, arr_ptr)
!             integer, intent(in) :: type, N
!             class(*), dimension(:), pointer :: arr_ptr
!         end subroutine
!     end interface
    
!     class(*), dimension(:), pointer :: int_ptr
!     class(*), dimension(:), pointer :: real_ptr

!         integer, dimension(:), pointer :: arr_int

!         allocate(arr_int(10))
!         int_ptr => arr_int 

!     ! call func(1, 5, int_ptr)
!     ! call func(2, 10, real_ptr)

!     print*, 'size ', size(int_ptr)
    
!     select type (int_ptr)
!         type is (integer)
!             int_ptr(0) = 100
!             int_ptr(1) = 5
!             print*, "im an integer", int_ptr(0)*int_ptr(1)
!         class default
!             print*, 'im not an integer'
!     end select
!     ! print*, 'aaaa ', int_ptr(0)


!     ! select type (int_ptr)
!     ! type is (integer)
!     !     print*, 'deallocate integer'
!     !     deallocate(int_ptr)
!     ! class default
!     !     print*, 'no deallocation'
!     ! end select


!     print*, 'end of program'
! end program

! subroutine func(type, N, arr_ptr)
!     integer, intent(in) :: type, N
!     class(*), dimension(:), pointer :: arr_ptr
!     integer, dimension(:), pointer :: arr_int
!     real(8), dimension(:), pointer :: arr_real

!     if (type == 1) then
!         allocate(arr_int(N))
!         arr_ptr => arr_int 
!     else
!         allocate(arr_real(N))
!         arr_ptr => arr_real
!     end if
! end subroutine







program main
    use iso_c_binding
    type(c_ptr) :: cptr1, cptr2
    integer, target, allocatable :: array(:)
    integer, target :: scalar
    integer, pointer :: pa(:), ps

    scalar = 50
    allocate(array(0:5))
    array = [1,2,3,4,5,6]
    cptr1 = c_loc(array(0))

    call c_f_pointer(cptr1, pa, shape=[6]) 

    print*, pa(3)
end program



! program main
!     class(*), dimension(:), pointer :: cptr1
!     integer, target, allocatable :: array(:)
!     integer, pointer :: new_ptr(:)

!     allocate(array(5))
!     array = [1,2,3,4,5]

!     cptr1 => array
!     new_ptr => cptr1
!     print*, array(3)
! end program

