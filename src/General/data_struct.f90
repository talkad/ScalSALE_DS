module data_struct_base

    use communication_module, only : communication_t
    use communication_parameters_module, only : communication_parameters_t
    use parallel_parameters_module, only: parallel_parameters_t

    type, abstract :: data_struct_t

        integer, public                                       :: nx   
        integer, public                                       :: ny   
        integer, public                                       :: nz
        integer, public                                       :: nmats
        
        type(communication_t), pointer                   :: communication
        type(communication_parameters_t), pointer        :: communication_parameters
        type (parallel_parameters_t), public, pointer    :: parallel_params  
        real(8), dimension(:), allocatable               :: send_buf
        real(8), dimension(:), allocatable               :: recv_buf
        integer :: request

    contains

        procedure, public :: Ptr_coordinates_1d
        procedure, public :: Ptr_coordinates_3d
        procedure, public :: Ptr_coordinates_4d

        generic, public :: Point_to_data => &
            Ptr_coordinates_1d, &
            Ptr_coordinates_3d, &
            Ptr_coordinates_4d

        procedure(get), public, deferred :: get_item
        procedure(add), public, deferred :: add_item
        procedure(print_struct), public, deferred :: print_data
        procedure(remove), public, deferred :: deallocate_data
        procedure(who), public, deferred :: who_am_i

        procedure, public :: Clean_data
        procedure, public :: Set_communication
        procedure, public :: Exchange_virtual_space_blocking
        procedure, private :: Set_send_buf
        procedure, private :: Get_recv_buf
        procedure, public :: Exchange_virtual_space_nonblocking
        procedure, public :: Exchange_end

        end type data_struct_t


        abstract interface
            function get(this, material_type, i, j, k)
                import data_struct_t
                class(data_struct_t), intent(in) :: this
                integer, intent(in) :: i, j, k, material_type
                real(8) :: get
            end function

            subroutine add(this, material_type, i, j, k, val)
                import data_struct_t
                class(data_struct_t), intent(inout) :: this
                integer, intent(in) :: i, j, k, material_type
                real(8), intent(in) :: val
            end subroutine

            subroutine print_struct(this, file_name)
                import data_struct_t
                class(data_struct_t), intent(inout) :: this
                character(len=*), intent(in)        ::   file_name
            end subroutine

            subroutine remove(this)
                import data_struct_t
                class(data_struct_t), intent(inout) :: this
            end subroutine

            subroutine who(this)
                import data_struct_t
                class(data_struct_t), intent(inout) :: this
            end subroutine
        end interface


        contains

        subroutine Clean_data(this)
            class (data_struct_t), intent(in out) :: this
        end subroutine Clean_data

        
        subroutine Ptr_coordinates_1d(this, ptr)
            class (data_struct_t), intent(in out) :: this
            real(8), dimension(:), pointer, intent(out) :: ptr
        end subroutine Ptr_coordinates_1d


        subroutine Ptr_coordinates_3d(this, ptr)
            class (data_struct_t), intent(in out) :: this
            real(8), dimension(:,:,:), pointer, intent(out) :: ptr
        end subroutine Ptr_coordinates_3d


        subroutine Ptr_coordinates_4d(this, ptr)
            class (data_struct_t), intent(in out) :: this
            real(8), dimension(:,:,:,:), pointer, intent(out) :: ptr
        end subroutine Ptr_coordinates_4d


        subroutine Set_communication(this, comm, comm_params)
            class (data_struct_t), intent(in out)                     :: this
            type(communication_t), pointer            :: comm
            type(communication_parameters_t), pointer :: comm_params
            ! integer, dimension(4) :: vals_shape
            integer :: d1,d2,d3
            
            if (associated(comm)) print*, 'harbu darbu 22222'
            if (associated(this%communication)) print*, 'harbu darbu 11111'
            

            print*, 'harbu darbu'
            this%communication => comm
            print*, 'harbu darbu'
            this%communication_parameters => comm_params
            print*, 'harbu darbu'
            this%parallel_params => this%communication%parallel_params
            if (this%communication%is_parallel .eqv. .true.) then
            print*, 'harbu darbu'
            ! vals_shape = shape(this%values)
            ! d1 = vals_shape(2) - 2
            ! d2 = vals_shape(3) - 2
            ! d3 = vals_shape(4) - 2
    
            d1 = this%nx - 2
            d2 = this%ny - 2
            d3 = this%nz - 2
            print*, 'okokok'
            allocate(this%send_buf(0:this%nmats * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+8)-1))
            allocate(this%recv_buf(0:this%nmats * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+8)-1))
            end if
        end subroutine Set_communication


        subroutine Exchange_virtual_space_blocking(this, ghost_width)
            class (data_struct_t), intent(in out) :: this
            real(8), dimension(:,:,:), allocatable :: send_buf, recv_buf
            integer, dimension(4) :: vals_shape
            integer, optional :: ghost_width
            integer :: ghost_width_local
            integer :: ghost_width_local_x
            integer :: i,j,k
            if (.not. present(ghost_width)) then
            ghost_width_local = 1   
            else
            ghost_width_local = 1   
            end if
            if (this%communication%is_parallel .eqv. .true.) then
            call this%Set_send_buf()
            call this%communication%Send_recv_neighbors_diag (this%communication_parameters, this%send_buf, this%recv_buf)
            call this%Get_recv_buf()
            end if
        end subroutine Exchange_virtual_space_blocking


        subroutine Exchange_virtual_space_nonblocking(this, ghost_width)
            class (data_struct_t), intent(in out) :: this
            integer, optional :: ghost_width
            integer :: ghost_width_local
    
            if (.not. present(ghost_width)) then
            ghost_width_local = 1   
            else
            ghost_width_local = 1   
            end if
            if (this%communication%is_parallel .eqv. .true.) then
            call this%Set_send_buf()
            call this%communication%Send_neighbors_diag (this%communication_parameters,&
                                                            this%send_buf, this%recv_buf, this%request)
    
            end if
        end subroutine Exchange_virtual_space_nonblocking


        subroutine Exchange_end(this)
            class (data_struct_t), intent(in out) :: this

            if (this%communication%is_parallel .eqv. .true.) then
                call this%communication%Wait_recv_neighbors_diag (this%communication_parameters,&
                                                                    this%send_buf, this%recv_buf, this%request)
                call this%Get_recv_buf()
                end if
        
        end subroutine Exchange_end


        subroutine Get_recv_buf(this)
            class (data_struct_t), intent(in out) :: this
        end subroutine Get_recv_buf


        subroutine Set_send_buf(this)
            class (data_struct_t), intent(in out) :: this
        end subroutine Set_send_buf



end module data_struct_base