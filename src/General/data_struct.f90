module data_struct_module

    use communication_module, only : communication_t
    use communication_parameters_module, only : communication_parameters_t
    use parallel_parameters_module, only: parallel_parameters_t

    type, abstract :: data_struct_t

        type(communication_t), pointer                   :: communication
        type(communication_parameters_t), pointer        :: communication_parameters
        type (parallel_parameters_t), public, pointer    :: parallel_params  
        real(8), dimension(:), allocatable               :: send_buf
        real(8), dimension(:), allocatable               :: recv_buf
        integer :: request
    contains

        procedure(Clean_data), public, deferred  :: Clean_data
        procedure(Set_communication), public, deferred  :: Set_communication
        procedure(Exchange_virtual_space_blocking), public, deferred  :: Exchange_virtual_space_blocking
        procedure(Set_send_buf), private, deferred :: Set_send_buf 
        procedure(Get_recv_buf), private, deferred :: Get_recv_buf
        procedure(Exchange_virtual_space_nonblocking), public, deferred  :: Exchange_virtual_space_nonblocking
        procedure(Exchange_end), public, deferred  :: Exchange_end
    end type data_struct_t



    abstract interface

        subroutine Clean_data(this)
            import data_struct_t
            class (data_struct_t), intent(in out) :: this
        end subroutine Clean_data


        subroutine Set_communication(this, comm, comm_params)
            import data_struct_t, communication_t, communication_parameters_t
            class (data_struct_t), intent(in out)                     :: this
            type(communication_t), pointer           :: comm
            type(communication_parameters_t), pointer     :: comm_params
        end subroutine Set_communication


        subroutine Exchange_virtual_space_blocking(this, ghost_width)
            import data_struct_t
            class (data_struct_t), intent(in out) :: this
            integer, optional :: ghost_width
        end subroutine Exchange_virtual_space_blocking


        subroutine Exchange_virtual_space_nonblocking(this, ghost_width)
            import data_struct_t
            class (data_struct_t), intent(in out) :: this
            integer, optional :: ghost_width
        end subroutine Exchange_virtual_space_nonblocking


        subroutine Exchange_end(this)
            import data_struct_t
            class (data_struct_t), intent(in out) :: this
        end subroutine Exchange_end


        subroutine Get_recv_buf(this)
            import data_struct_t
            class (data_struct_t), intent(in out) :: this
        end subroutine Get_recv_buf


        subroutine Set_send_buf(this)
            import data_struct_t
            class (data_struct_t), intent(in out) :: this
        end subroutine Set_send_buf

    end interface


end module data_struct_module