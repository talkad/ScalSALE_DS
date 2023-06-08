module naive

   type, abstract :: data_type

        contains 
            procedure(init_interface), deferred :: init

            procedure(get_value_interface), deferred :: get_value
            procedure(get_send_buf_interface), deferred :: get_send_buf
            procedure(get_recv_buf_interface), deferred :: get_recv_buf

            procedure(set_value_interface), deferred :: set_value
            procedure(set_send_buf_interface), deferred :: set_send_buf
            procedure(set_recv_buf_interface), deferred :: set_recv_buf
   end type   data_type


    abstract interface
        pure function get_value_interface(this, m, i, j, k) result(value)
            import :: data_type
            class(data_type), intent(in) :: this
            integer, intent(in) :: m, i, j, k
            real(8) :: value
        end function get_value_interface

        pure function get_send_buf_interface(this, i) result(value)
            import :: data_type
            class(data_type), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value
        end function get_send_buf_interface

        function get_recv_buf_interface(this, i) result(value)
            import :: data_type
            class(data_type), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value
        end function get_recv_buf_interface

        subroutine set_value_interface(this, m, i, j, k, value)
            import :: data_type
            class(data_type), intent(inout) :: this
            real(8), intent(in) :: value
            integer, intent(in) :: m, i, j, k
        end subroutine set_value_interface

        subroutine set_send_buf_interface(this, i, value)
            import :: data_type
            class(data_type), intent(inout) :: this
            real(8), intent(in) :: value
            integer, intent(in) :: i
        end subroutine set_send_buf_interface

        subroutine set_recv_buf_interface(this, i, value)
            import :: data_type
            class(data_type), intent(inout) :: this
            real(8), intent(in) :: value
            integer, intent(in) :: i
        end subroutine set_recv_buf_interface
    end interface



   type, extends(data_type) :: data_real
      real(8), dimension(:,:,:,:), pointer, public          :: values
      real(8), dimension(:), allocatable     :: send_buf
      real(8), dimension(:), allocatable     :: recv_buf

      contains 
         procedure, public :: get_value => get_value_real
         procedure, public :: get_send_buf => get_send_buf_real
         procedure, public :: get_recv_buf => get_recv_buf_real

         procedure, public :: set_value => set_value_real
         procedure, public :: set_send_buf => set_send_buf_real
         procedure, public :: set_recv_buf => set_recv_buf_real
   end type   data_real

   type, extends(data_type) :: data_int
      integer, dimension(:,:,:,:), pointer, public          :: values
      integer, dimension(:), allocatable     :: send_buf
      integer, dimension(:), allocatable     :: recv_buf

      contains 
         procedure, public :: get_value => get_value_int
         procedure, public :: get_send_buf => get_send_buf_int
         procedure, public :: get_recv_buf => get_recv_buf_int

         procedure, public :: set_value => set_value_int
         procedure, public :: set_send_buf => set_send_buf_int
         procedure, public :: set_recv_buf => set_recv_buf_int
   end type   data_int


    contains

        pure function get_value_real(this, m, i, j, k) result(value)
            class(data_real), intent(in) :: this
            integer, intent(in) :: m, i, j, k
            real(8) :: value

            value = this%values(m,i,j,k)
        end function get_value_real

        pure function get_send_buf_real(this, i) result(value)
            class(data_real), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value

            value = this%send_buf(i)
        end function get_send_buf_real

        function get_recv_buf_real(this, i) result(value)
            class(data_real), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value

            value = this%recv_buf(i)
        end function get_recv_buf_real

        subroutine set_value_real(this, m, i, j, k, value)
            class(data_real), intent(inout) :: this
            integer, intent(in) ::  m, i, j, k
            real(8), intent(in) :: value

            this%values(m,i,j,k) = value
        end subroutine set_value_real

        subroutine set_send_buf_real(this, i, value)
            class(data_real), intent(inout) :: this
            integer, intent(in) ::  i
            real(8), intent(in) :: value

            this%send_buf(i) = value
        end subroutine set_send_buf_real

        subroutine set_recv_buf_real(this, i, value)
            class(data_real), intent(inout) :: this
            integer, intent(in) ::  i
            real(8), intent(in) :: value

            this%recv_buf(i) = value
        end subroutine set_recv_buf_real





      pure function get_value_int(this, m, i, j, k) result(value)
            class(data_int), intent(in) :: this
            integer, intent(in) :: m, i, j, k
            real(8) :: value

            value = real(this%values(m,i,j,k))
        end function get_value_int

        pure function get_send_buf_int(this, i) result(value)
            class(data_int), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value

            value = real(this%send_buf(i))
        end function get_send_buf_int

        function get_recv_buf_int(this, i) result(value)
            class(data_int), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value

            value = real(this%recv_buf(i))
        end function get_recv_buf_int

        subroutine set_value_int(this, m, i, j, k, value)
            class(data_int), intent(inout) :: this
            integer, intent(in) ::  m, i, j, k
            real(8), intent(in) :: value

            this%values(m,i,j,k) = int(value)
        end subroutine set_value_int

        subroutine set_send_buf_int(this, i, value)
            class(data_int), intent(inout) :: this
            integer, intent(in) ::  i
            real(8), intent(in) :: value

            this%send_buf(i) = int(value)
        end subroutine set_send_buf_int

        subroutine set_recv_buf_int(this, i, value)
            class(data_int), intent(inout) :: this
            integer, intent(in) ::  i
            real(8), intent(in) :: value

            this%recv_buf(i) = int(value)
        end subroutine set_recv_buf_int

end module

