
module data_4d_module
   use communication_module, only : communication_t
   use communication_parameters_module, only : communication_parameters_t
   use parallel_parameters_module, only: parallel_parameters_t
   implicit none
   private
   public :: data_4d_t

   type, abstract :: abstract_struct_t

        contains 
            procedure(init_interface), deferred :: init

            procedure(get_value_interface), deferred :: get_value
            procedure(get_send_buf_interface), deferred :: get_send_buf
            procedure(get_recv_buf_interface), deferred :: get_recv_buf

            procedure(set_value_interface), deferred :: set_value
            procedure(set_send_buf_interface), deferred :: set_send_buf
            procedure(set_recv_buf_interface), deferred :: set_recv_buf
   end type   abstract_struct_t


    abstract interface
        subroutine init_interface(this, d1, d2, d3, d4, is_parallel)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(inout) :: this
            integer, intent(in) :: d1, d2, d3, d4
            logical, intent(in) :: is_parallel
        end subroutine init_interface

        pure function get_value_interface(this, m, i, j, k) result(value)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(in) :: this
            integer, intent(in) :: m, i, j, k
            real(8) :: value
        end function get_value_interface

        pure function get_send_buf_interface(this, i) result(value)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value
        end function get_send_buf_interface

        function get_recv_buf_interface(this, i) result(value)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(in) :: this
            integer, intent(in) :: i
            real(8) :: value
        end function get_recv_buf_interface

        subroutine set_value_interface(this, m, i, j, k, value)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(inout) :: this
            real(8), intent(in) :: value
            integer, intent(in) :: m, i, j, k
        end subroutine set_value_interface

        subroutine set_send_buf_interface(this, i, value)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(inout) :: this
            real(8), intent(in) :: value
            integer, intent(in) :: i
        end subroutine set_send_buf_interface

        subroutine set_recv_buf_interface(this, i, value)
            import :: abstract_struct_t
            class(abstract_struct_t), intent(inout) :: this
            real(8), intent(in) :: value
            integer, intent(in) :: i
        end subroutine set_recv_buf_interface
    end interface



   type, extends(abstract_struct_t) :: data_real
      real(8), dimension(:,:,:,:), pointer, public          :: values
      real(8), dimension(:), allocatable     :: send_buf
      real(8), dimension(:), allocatable     :: recv_buf

      contains 
         procedure, public :: init => init_real

         procedure, public :: get_value => get_value_real
         procedure, public :: get_send_buf => get_send_buf_real
         procedure, public :: get_recv_buf => get_recv_buf_real

         procedure, public :: set_value => set_value_real
         procedure, public :: set_send_buf => set_send_buf_real
         procedure, public :: set_recv_buf => set_recv_buf_real
   end type   data_real

   type, extends(abstract_struct_t) :: data_int
      integer, dimension(:,:,:,:), pointer, public          :: values
      integer, dimension(:), allocatable     :: send_buf
      integer, dimension(:), allocatable     :: recv_buf

      contains 
         procedure, public :: init => init_int

         procedure, public :: get_value => get_value_int
         procedure, public :: get_send_buf => get_send_buf_int
         procedure, public :: get_recv_buf => get_recv_buf_int

         procedure, public :: set_value => set_value_int
         procedure, public :: set_send_buf => set_send_buf_int
         procedure, public :: set_recv_buf => set_recv_buf_int
   end type   data_int


      


   type :: data_4d_t
      private

      class(abstract_struct_t), pointer :: storage_struct
      ! real(8), dimension(:,:,:,:), pointer, public          :: values
      ! real(8), dimension(:), allocatable     :: send_buf
      ! real(8), dimension(:), allocatable     :: recv_buf

      integer, public                                       :: nx   
      integer, public                                       :: ny   
      integer, public                                       :: nz
      integer, public                                       :: nmats
      type(communication_t), pointer                        :: communication
      type(communication_parameters_t), pointer :: communication_parameters
      type (parallel_parameters_t)   , public, pointer :: parallel_params      

      integer :: pre_calc(14)

      integer :: request
   contains


      procedure, public :: Point_to_data => Ptr_data

      procedure, public :: Clean_data

      procedure, public :: Set_communication

      procedure, public :: Debug_check_nan

      procedure, public :: Exchange_virtual_space_blocking


      procedure, private ::Set_send_buf

      procedure, private ::Get_recv_buf


      procedure, public :: Exchange_virtual_space_nonblocking

      procedure, public :: Exchange_end

      procedure, public :: debug_print

      procedure, public :: Write_data
      generic :: write(unformatted) => Write_data


      procedure, public :: Read_data
      generic :: read(unformatted) => Read_data

   end type data_4d_t

   public :: Get_copy

   interface data_4d_t
      module procedure Constructor_init_arr
      module procedure Constructor_init_val


   end interface data_4d_t

contains

   type(data_4d_t) function Constructor_init_arr(type, d1, d2, d3, d4, is_parallel)
      implicit none
      integer                  , intent(in) :: d1, d2, d3, d4         
      integer                  , intent(in) :: type         
      logical,                 , intent(in) :: is_parallel
      type(data_real), target :: data_concrete_real
      type(data_int), target :: data_concrete_int

      if (type == 0) then
         call data_concrete_real%init(d1, d2, d3, d4, is_parallel)
         
      else 

      end if

      Constructor_init_arr%nx = d1
      Constructor_init_arr%ny = d2
      Constructor_init_arr%nz = d3
      Constructor_init_arr%nmats = d4
  end function

   type(data_4d_t) function Constructor_init_val(initial_val, d1, d2, d3, d4)
      implicit none
      real(8)           , intent(in) :: initial_val  
      integer           , intent(in) :: d1           
      integer           , intent(in) :: d2           
      integer           , intent(in) :: d3           
      integer                  , intent(in) :: d4

      allocate (Constructor_init_val%values (1:d4, 0:d1, 0:d2, 0:d3))
      Constructor_init_val%values = initial_val
      Constructor_init_val%nx = d1
      Constructor_init_val%ny = d2
      Constructor_init_val%nz = d3
   Constructor_init_val%nmats = d4

   end function


   subroutine Ptr_data (this, ptr)
      class (data_4d_t)                    , intent(in)  :: this
      real(8), dimension(:,:,:,:), pointer, intent(out) :: ptr

      ptr => this%values
   end subroutine Ptr_data

   function Get_copy (this)
      class (data_4d_t)       , intent(in)  :: this
      real(8), dimension(:,:,:,:), pointer :: Get_copy

      Get_copy = this%values
   end function Get_copy

   subroutine Clean_data (this)
      class (data_4d_t), intent(in out) :: this

      deallocate (this%values)
   end subroutine Clean_data


   subroutine Set_communication (this, comm, comm_params)
      class (data_4d_t), intent(in out) :: this
      type(communication_t), pointer            :: comm
      type(communication_parameters_t), pointer :: comm_params
      integer, dimension(4) :: vals_shape
      integer :: d1,d2,d3

      this%communication => comm
      this%communication_parameters => comm_params

      this%parallel_params => this%communication%parallel_params
      if (this%communication%is_parallel .eqv. .true.) then
         vals_shape = shape(this%values)
         d1 = vals_shape(2) - 2
         d2 = vals_shape(3) - 2
         d3 = vals_shape(4) - 2

         this%pre_calc(1) = d1+2
         this%pre_calc(2) = d3+2
         this%pre_calc(3) = d2+2
         this%pre_calc(4) = this%pre_calc(3)*this%pre_calc(2)
         this%pre_calc(5) = 2*this%pre_calc(4)
         this%pre_calc(6) = this%pre_calc(1)*this%pre_calc(2)
         this%pre_calc(7) = 2*this%pre_calc(6)
         this%pre_calc(8) = this%pre_calc(1)*this%pre_calc(3)
         this%pre_calc(9) = 2*this%pre_calc(8)
         this%pre_calc(10) = 4*this%pre_calc(3)+4*this%pre_calc(1)
         this%pre_calc(11) = this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(9)
         this%pre_calc(12) = 4*this%pre_calc(2)
         this%pre_calc(13) = this%pre_calc(11)+this%pre_calc(12)
         this%pre_calc(14) = this%pre_calc(13)+this%pre_calc(10) 


         allocate(this%send_buf(0:this%nmats * (this%pre_calc(14)+8)-1))
         allocate(this%recv_buf(0:this%nmats * (this%pre_calc(14)+8)-1))
      end if
   end subroutine Set_communication


   subroutine Exchange_virtual_space_blocking (this, ghost_width)
      class (data_4d_t), intent(in out) :: this
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




   subroutine Exchange_virtual_space_nonblocking (this, ghost_width)
      class (data_4d_t), intent(in out) :: this
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

   subroutine Exchange_end (this)
      class (data_4d_t), intent(in out) :: this

      if (this%communication%is_parallel .eqv. .true.) then
         call this%communication%Wait_recv_neighbors_diag (this%communication_parameters,&
                                                           this%send_buf, this%recv_buf, this%request)
         call this%Get_recv_buf()
      end if

   end subroutine Exchange_end


   subroutine Get_recv_buf(this)
      class (data_4d_t), intent(in out) :: this
      integer, dimension(4) :: vals_shape
      integer :: d1, d2, d3, x, y, z, m
      m=this%nmats

      vals_shape = shape(this%values)
      d1 = vals_shape(2) - 2
      d2 = vals_shape(3) - 2
      d3 = vals_shape(4) - 2

      x = this%parallel_params%my_coords(1)
      y = this%parallel_params%my_coords(2)
      z = this%parallel_params%my_coords(3)

      !from right
      if (x+1 /= this%parallel_params%npx+1) then
         this%values(1:m, d1+1, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(0 : m * (this%pre_calc(4))-1), (/m, d2+2, d3+2/))
      end if

      !from left
      if (x-1 /= 0) then
         this%values(1:m, 0, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(4)) :&
                     m * (this%pre_calc(5))-1), (/m, d2+2, d3+2/))
      end if

      !from front
      if (y+1 /= this%parallel_params%npy+1) then
         this%values(1:m, 0:d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(5)) :&
                     m * (this%pre_calc(5)+this%pre_calc(6))-1), (/m, this%pre_calc(1), d3+2/))
      end if


      !from rear
      if (y-1 /= 0) then
         this%values(1:m, 0:d1+1, 0, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(5)+this%pre_calc(6)) :&
                     m * (this%pre_calc(5)+this%pre_calc(7))-1), (/m, d1+2, d3+2/))
      end if

      !from up
      if (z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0:d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(5)+this%pre_calc(7)) :&
                     m * (this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(8))-1), (/m, d1+2, d2+2/))
      end if

      !from down
      if (z-1 /= 0) then
         this%values(1:m, 0:d1+1, 0:d2+1, 0) = reshape(this%recv_buf(m * (this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(8)) :&
                     m * (this%pre_calc(11))-1), (/m, d1+2, d2+2/))
      end if

      !from diagonal of i+1, j+1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
         this%values(1:m, d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)) :&
                     m * (this%pre_calc(11)+this%pre_calc(2))-1), (/m, d3+2/))
      end if

      !from diagonal of i+1, j-1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
         this%values(1:m, d1+1, 0, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)+this%pre_calc(2)) :&
                     m * (this%pre_calc(11)+2*this%pre_calc(2))-1), (/m, d3+2/))
      end if

      !from diagonal of i-1, j+1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
         this%values(1:m, 0, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)+2*this%pre_calc(2)) :&
                     m * (this%pre_calc(11)+3*this%pre_calc(2))-1), (/m, d3+2/))
      end if

      !from diagonal of i-1, j-1
      if (x-1 /= 0 .and. y-1 /= 0) then
         this%values(1:m, 0, 0, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)+3*this%pre_calc(2)) :&
                     m * (this%pre_calc(13))-1), (/m, d3+2/))
      end if

      !from diagonal of i+1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)) :&
                     m * (this%pre_calc(13)+this%pre_calc(3))-1), (/m, d2+2/))
      end if

      !from diagonal of i+1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
         this%values(1:m, d1+1, 0:d2+1, 0) = reshape(this%recv_buf(m * (this%pre_calc(13)+this%pre_calc(3)) :&
                     m * (this%pre_calc(13)+2*this%pre_calc(3))-1), (/m, d2+2/))
      end if

      !from diagonal of i-1, k+1
      if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)+2*this%pre_calc(3)) :&
                     m * (this%pre_calc(13)+3*this%pre_calc(3))-1), (/m, d2+2/))
      end if

      !from diagonal of i-1, k-1
      if (x-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, 0, 0:d2+1, 0) = reshape(this%recv_buf(m* (this%pre_calc(13)+3*this%pre_calc(3)) :&
                     m * (this%pre_calc(13)+4*this%pre_calc(3))-1), (/m, d2+2/))
      end if

      !from diagonal of j+1, k+1
      if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0:d1+1, d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)) :&
                     m * (this%pre_calc(13)+4*this%pre_calc(3)+this%pre_calc(1))-1), (/m, d1+2/))
      end if

      !from diagonal of j+1, k-1
      if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(1:m, 0:d1+1, d2+1, 0) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)+this%pre_calc(1)) :&
                     m* (this%pre_calc(13)+4*this%pre_calc(3)+2*this%pre_calc(1))-1), (/m, d1+2/))
      end if

      !from diagonal of j-1, k+1
      if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0:d1+1, 0, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)+2*this%pre_calc(1)) :&
                     m * (this%pre_calc(13)+4*this%pre_calc(3)+3*this%pre_calc(1))-1), (/m, d1+2/))
      end if

      !from diagonal of j-1, k-1
      if (y-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, 0:d1+1, 0, 0) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)+3*this%pre_calc(1)) :&
                     m * (this%pre_calc(14))-1), (/m, d1+2/))
      end if

      !from diagonal of i+1, j+1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, d1+1, d2+1, d3+1) = this%recv_buf(m * (this%pre_calc(14)) :&
                     m * (this%pre_calc(14)+1)-1)
      end if

      !from diagonal of i+1, j+1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(1:m, d1+1, d2+1, 0) = this%recv_buf(m * (this%pre_calc(14)+1) :&
                     m * (this%pre_calc(14)+2)-1)
      end if

      !from diagonal of i+1, j-1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, d1+1, 0, d3+1) = this%recv_buf(m * (this%pre_calc(14)+2) :&
                     m * (this%pre_calc(14)+3)-1)
      end if

      !from diagonal of i+1, j-1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, d1+1, 0, 0) = this%recv_buf(m * (this%pre_calc(14)+3) :&
                     m * (this%pre_calc(14)+4)-1)
      end if

      !from diagonal of i-1, j+1, k+1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0, d2+1, d3+1) = this%recv_buf(m * (this%pre_calc(14)+4) :&
                     m * (this%pre_calc(14)+5)-1)
      end if

      !from diagonal of i-1, j+1, k-1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(1:m, 0, d2+1, 0) = this%recv_buf(m * (this%pre_calc(14)+5) :&
                     m * (this%pre_calc(14)+6)-1)
      end if

      !from diagonal of i-1, j-1, k+1
      if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0, 0, d3+1) = this%recv_buf(m * (this%pre_calc(14)+6) :&
                     m * (this%pre_calc(14)+7)-1)
      end if

      !from diagonal of i-1, j-1, k-1
      if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, 0, 0, 0) = this%recv_buf(m * (this%pre_calc(14)+7) :&
                     m * (this%pre_calc(14)+8)-1)
      end if
   end subroutine Get_recv_buf



   subroutine Set_send_buf(this)
      class (data_4d_t), intent(in out) :: this
      integer, dimension(4) :: vals_shape
      integer :: d1, d2, d3, x, y, z, offset
      integer :: m

      m = this%nmats
      vals_shape = shape(this%values)
      d1 = vals_shape(2) - 2
      d2 = vals_shape(3) - 2
      d3 = vals_shape(4) - 2

      x = this%communication%parallel_params%my_coords(1)
      y = this%communication%parallel_params%my_coords(2)
      z = this%communication%parallel_params%my_coords(3)
      offset = this%communication_parameters%dim_offset

      !to right
      if (x+1 /= this%parallel_params%npx+1) then
         this%send_buf(0 : m * (this%pre_calc(4))-1) = &
                        reshape(this%values(1:m, d1 - offset, 0:d2+1, 0:d3+1), (/m * this%pre_calc(3) * this%pre_calc(2)/))
      end if

      !to left
      if (x-1 /= 0) then
         this%send_buf(m * (this%pre_calc(4)) : m * (this%pre_calc(5))-1) = &
                        reshape(this%values(1:m, 1 + offset, 0:d2+1, 0:d3+1), (/m * this%pre_calc(3) * this%pre_calc(2)/))
      end if

      !to front
      if (y+1 /= this%parallel_params%npy+1) then
         this%send_buf(m * (this%pre_calc(5)) : m * (this%pre_calc(5)+this%pre_calc(6))-1) = &
                        reshape(this%values(1:m, 0:d1+1, d2 - offset, 0:d3+1), (/m * this%pre_calc(1) * this%pre_calc(2)/))
      end if

      !to rear
      if (y-1 /= 0) then
         this%send_buf(m * (this%pre_calc(5)+this%pre_calc(6)) : m * (this%pre_calc(5)+this%pre_calc(7))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 1 + offset, 0:d3+1), (/m * this%pre_calc(1) * this%pre_calc(2)/))
      end if

      !to up
      if (z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(5)+this%pre_calc(7)) : m * (this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(8))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 0:d2+1, d3 - offset), (/m * this%pre_calc(1) * this%pre_calc(3)/))
      end if

      !to down
      if (z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(8)) : &
                        m * (this%pre_calc(11))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 0:d2+1, 1 + offset), (/m * this%pre_calc(1) * this%pre_calc(3)/))
      end if

      !to diagonal of i+1, j+1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
         this%send_buf(m * (this%pre_calc(11)) : &
                        m * (this%pre_calc(11)+this%pre_calc(2))-1) = &
                        reshape(this%values(1:m, d1 - offset, d2 - offset, 0:d3+1), (/m * this%pre_calc(2)/))
      end if

      !to diagonal of i+1, j-1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
         this%send_buf(m * (this%pre_calc(11)+this%pre_calc(2)) : &
                        m * (this%pre_calc(11)+2*this%pre_calc(2))-1) = &
                        reshape(this%values(1:m, d1 - offset, 1 + offset, 0:d3+1), (/m * this%pre_calc(2)/))
      end if

      !to diagonal of i-1, j+1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
         this%send_buf(m * (this%pre_calc(11)+2*this%pre_calc(2)) : &
                        m * (this%pre_calc(11)+3*this%pre_calc(2))-1) = &
                        reshape(this%values(1:m, 1 + offset, d2 - offset, 0:d3+1), (/m * this%pre_calc(2)/))
      end if

      !to diagonal of i-1, j-1
      if (x-1 /= 0 .and. y-1 /= 0) then
         this%send_buf(m * (this%pre_calc(11)+3*this%pre_calc(2)) : &
                        m * (this%pre_calc(13))-1) = &
                        reshape(this%values(1:m, 1+offset, 1+offset, 0:d3+1), (/m * this%pre_calc(2)/))
      end if

      !to diagonal of i+1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(13)) : &
                        m * (this%pre_calc(13)+this%pre_calc(3))-1) = &
                        reshape(this%values(1:m, d1 - offset, 0:d2+1, d3 - offset), (/m * this%pre_calc(3)/))
      end if

      !to diagonal of i+1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(13)+this%pre_calc(3)) : &
                        m * (this%pre_calc(13)+2*this%pre_calc(3))-1) = &
                        reshape(this%values(1:m, d1 - offset, 0:d2+1, 1 + offset), (/m * this%pre_calc(3)/))
      end if

      !to diagonal of i-1, k+1
      if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(13)+2*this%pre_calc(3)) : &
                        m * (this%pre_calc(13)+3*this%pre_calc(3))-1) = &
                        reshape(this%values(1:m, 1 + offset, 0:d2+1, d3 - offset), (/m * this%pre_calc(3)/))
      end if

      !to diagonal of i-1, k-1
      if (x-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(13)+3*this%pre_calc(3)) : &
                        m * (this%pre_calc(13)+4*this%pre_calc(3))-1) = &
                        reshape(this%values(1:m, 1 + offset, 0:d2+1, 1 + offset), (/m * this%pre_calc(3)/))
      end if

      !to diagonal of j+1, k+1
      if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)) : &
                        m* (this%pre_calc(13)+4*this%pre_calc(3)+this%pre_calc(1))-1) = &
                        reshape(this%values(1:m, 0:d1+1, d2 - offset, d3 - offset), (/m * this%pre_calc(1)/))
      end if

      !to diagonal of j+1, k-1
      if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)+this%pre_calc(1)) : &
                        m * (this%pre_calc(13)+4*this%pre_calc(3)+2*this%pre_calc(1))-1) = &
                        reshape(this%values(1:m, 0:d1+1, d2 - offset, 1 + offset), (/m * this%pre_calc(1)/))
      end if

      !to diagonal of j-1, k+1
      if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)+2*this%pre_calc(1)) : &
                        m * (this%pre_calc(13)+4*this%pre_calc(3)+3*this%pre_calc(1))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 1 + offset, d3 - offset), (/m * this%pre_calc(1)/))
      end if

      !to diagonal of j-1, k-1
      if (y-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(13)+4*this%pre_calc(3)+3*this%pre_calc(1)) : &
                        m * (this%pre_calc(14))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 1 + offset, 1+offset), (/m * this%pre_calc(1)/))
      end if


      !to diagonal of i+1, j+1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(14)) :&
                        m * (this%pre_calc(14)+1)-1) = &
                        this%values(1:m, d1- offset, d2-offset, d3-offset)
      end if

      !to diagonal of i+1, j+1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(14)+1) :&
                        m * (this%pre_calc(14)+2)-1) = &
                        this%values(1:m, d1-offset, d2-offset, 1 + offset)
      end if

      !to diagonal of i+1, j-1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(14)+2) :&
                        m * (this%pre_calc(14)+3)-1) = &
                        this%values(1:m, d1-offset, 1+offset, d3-offset)
      end if

      !to diagonal of i+1, j-1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(14)+3) :&
                        m * (this%pre_calc(14)+4)-1) = &
                        this%values(1:m, d1-offset, 1+offset, 1+offset)
      end if

      !to diagonal of i-1, j+1, k+1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(14)+4) :&
                        m * (this%pre_calc(14)+5)-1) = &
                        this%values(1:m, 1+offset, d2-offset, d3-offset)
      end if

      !to diagonal of i-1, j+1, k-1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(14)+5) :&
                        m * (this%pre_calc(14)+6)-1) = &
                        this%values(1:m, 1+offset, d2-offset, 1+offset)
      end if

      !to diagonal of i-1, j-1, k+1
      if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(m * (this%pre_calc(14)+6) :&
                        m * (this%pre_calc(14)+7)-1) = &
                        this%values(1:m, 1+offset, 1+offset, d3-offset)
      end if

      !to diagonal of i-1, j-1, k-1
      if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(m * (this%pre_calc(14)+7) :&
                        m * (this%pre_calc(14)+8)-1) = &
                        this%values(1:m, 1+offset, 1+offset, 1+offset)
      end if
   end subroutine Set_send_buf


   subroutine Debug_check_nan(this, caller)
      implicit none
      class (data_4d_t), intent(in out) :: this
      CHARACTER(*) caller

      integer :: i,j,k
      do i = 0, this%nx
         do j = 0, this%ny
            do k = 0, this%nz
!               if (this%values(i,j,k) /= this%values(i,j,k)) then
                  write(*,*) "NaN in ", caller, ": ",i,j,k
                  return
!               end if
            end do
         end do
      end do
   end subroutine Debug_check_nan

   subroutine debug_print(this, caller, flag)
      implicit none
      class (data_4d_t), intent(in out) :: this
      integer, optional :: flag
      CHARACTER(*) caller
      integer :: width
      integer :: i,j,k
      if (.not. present(flag)) then
         width = 0   
      else
         width = flag   
      end if

      write(69,*) "---- data_t:", caller, "----"
      do i = width, this%nx - width
         do j = width, this%ny - width
            do k = width, this%nz - width
!                  write(69,*) i,j,k,this%values(i,j,k)
            end do
         end do
      end do
   end subroutine debug_print

   subroutine Write_data(this, unit, iostat, iomsg)
      class (data_4d_t), intent(in) :: this
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_data @@@'
#endif

      write(unit, iostat=iostat, iomsg=iomsg) &
         shape(this%values)

      write(unit, iostat=iostat, iomsg=iomsg) &
         this%values, &
         this%nx, &
         this%ny, &
         this%nz

#ifdef DEBUG
      write(*,*) &
         'shape_values', &
         shape(this%values)
      write(*,*) &
         'values', &
         this%values, &
         'nx', &
         this%nx, &
         'ny', &
         this%ny, &
         'nz', &
         this%nz, &
         '###'

      write(*,*) '@@@ end Write_data @@@'
#endif

   end subroutine Write_data

   subroutine Read_data(this, unit, iostat, iomsg)
      class (data_4d_t), intent(in out) :: this
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg
      integer, dimension(3) :: shape_values

#ifdef DEBUG
      write(*,*) "@@@ in Read_data @@@"
#endif

      read(unit, iostat=iostat, iomsg=iomsg) &
         shape_values

      deallocate(this%values)
!      allocate(this%values(0:shape_values(1) - 1, 0:shape_values(2)-1, 0:shape_values(3) - 1))

      read(unit, iostat=iostat, iomsg=iomsg) &
         this%values, &
         this%nx, &
         this%ny, &
         this%nz

#ifdef DEBUG
      write(*,*) &
         'shape_values', &
         shape(this%values)
      write(*,*) &
         'values', &
         this%values, &
         'nx', &
         this%nx, &
         'ny', &
         this%ny, &
         'nz', &
         this%nz, &
         '###'

      write(*,*) "@@@ end Read_data @@@"
#endif

   end subroutine Read_data



     subroutine init_real(this, d1, d2, d3, d4, is_parallel)
            class(data_real), intent(inout) :: this
            integer, intent(in) :: d1, d2, d3, d4
            logical, intent(in) :: is_parallel

            allocate (this%values (1:d4, 0:d1, 0:d2, 0:d3))

            if (is_parallel) then
                ! allocate buffers
            end if
        end subroutine init_real

        pure function get_value_real(this, m, i, j, k) result(value)
            class(data_real), intent(in) :: this
            integer, intent(in) :: m, i, j, k
            real(8) :: value

            valuedata_int = this%values(m,i,j,k)
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



        subroutine init_int(this, d1, d2, d3, d4, is_parallel)
            class(data_int), intent(inout) :: this
            integer, intent(in) :: d1, d2, d3, d4
            logical, intent(in) :: is_parallel

            allocate (this%values (1:d4, 0:d1, 0:d2, 0:d3))

            if (is_parallel) then
                ! allocate buffers
            end if
        end subroutine init_int

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

end module data_4d_module
