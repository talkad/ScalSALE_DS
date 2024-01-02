module block_csr_module
    use data_struct_base


    type :: block_t
        real(8), dimension(:,:,:), pointer :: matrix
    end type block_t


    type, extends(data_struct_t) :: block_csr_t

        real(8), dimension(:,:,:,:), pointer, public         :: values

        type(block_t), dimension(:,:,:,:), allocatable :: grid

        integer :: block_size = 8
        integer :: d1_block, d2_block, d3_block

    contains


        procedure, public :: Ptr_coordinates_4d => Ptr_coordinates_4d_block

        procedure, public :: Clean_data => Clean_data_imp

        procedure, public :: Debug_check_nan

        procedure, private ::Set_send_buf => Set_send_buf_imp

        procedure, private ::Get_recv_buf => Get_recv_buf_imp


        procedure, public :: debug_print

        procedure, public :: Write_data
        generic :: write(unformatted) => Write_data


        procedure, public :: Read_data
        generic :: read(unformatted) => Read_data




        
        procedure :: get_item => get_item
        procedure :: add_item => add_item
        procedure :: print_data => print_data
        procedure :: deallocate_data => deallocate_data
        procedure :: who_am_i => who_am_i
        procedure :: reorder => reorder

        procedure, private :: add_boundary

    end type

    public :: Get_copy


    interface block_csr_t
        module procedure Constructor_init_val
    end interface


    contains
    

    function Constructor_init_val(initial_val, d1, d2, d3, d4)
        implicit none

        type(block_csr_t), pointer :: Constructor_init_val
        real(8)           , intent(in) :: initial_val  
        integer           , intent(in) :: d1, d2, d3, d4  
        integer                        :: k,j,i,m, block_size

        print*, 'init_block_csr'

        allocate(Constructor_init_val)
        Constructor_init_val%communication => null()

        allocate(Constructor_init_val%values(0:1,0:1,0:1,0:1))
        Constructor_init_val%values = 0

        Constructor_init_val%nx = d1
        Constructor_init_val%ny = d2
        Constructor_init_val%nz = d3
        Constructor_init_val%nmats = d4

        block_size = Constructor_init_val%block_size

        Constructor_init_val%d1_block = d1/block_size 
        Constructor_init_val%d2_block = d2/block_size 
        Constructor_init_val%d3_block = d3/block_size 

        allocate(Constructor_init_val%grid(1:d4, 0:Constructor_init_val%d1_block, 0:Constructor_init_val%d2_block, 0:Constructor_init_val%d3_block))  ! m,i,j,k

        do k = 0, Constructor_init_val%d3_block
            do j = 0, Constructor_init_val%d2_block
                do i = 0, Constructor_init_val%d1_block
                    do m = 1, Constructor_init_val%nmats
                        Constructor_init_val%grid(m, i, j ,k)%matrix => null()
                    end do
                end do
            end do
        end do        

        call Constructor_init_val%add_boundary()

    end function


    subroutine Ptr_coordinates_4d_block (this, ptr)
        class (block_csr_t)                    , intent(in out)  :: this
        real(8), dimension(:,:,:,:), pointer, intent(out) :: ptr

        ptr => this%values
    end subroutine Ptr_coordinates_4d_block


    subroutine add_boundary(this)
        implicit none
        class(block_csr_t), intent(inout) :: this
        integer :: i,j,k,m
                
        i = 0
        do k = 0, this%nz
            do j = 0, this%ny
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0, .true.)
                end do
            end do
        end do


        i = this%nx
        do k = 0, this%nz
            do j = 0, this%ny
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0, .true.)
                end do
            end do
        end do


        j = 0
        do k = 0, this%nz
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0, .true.)
                end do
            end do
        end do


        j = this%ny
        do k = 0, this%nz
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0, .true.)
                end do
            end do
        end do


        k = 0
        do j = 0, this%ny
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0, .true.)
                end do
            end do
        end do

        k = this%nz
        do j = 0, this%ny
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0, .true.)
                end do
            end do
        end do

    
    end subroutine add_boundary


    subroutine add_item(this, material_type, i, j, k, val, boundry)
        implicit none
        class(block_csr_t), intent(inout) :: this
        integer, intent(in) :: i, j, k, material_type
        real(8), intent(in) :: val
        integer :: i_new, j_new, k_new
        real(8), dimension(:,:,:), allocatable, target :: matrix
        logical, optional, intent(in) ::  boundry

        i_new = i/this%block_size
        j_new = j/this%block_size
        k_new = k/this%block_size

        if (.not. associated(this%grid(material_type, i_new, j_new ,k_new)%matrix )) then 
            if (.not. present(boundry) .and. val == 0) return
            allocate(this%grid(material_type, i_new, j_new ,k_new)%matrix(0:this%block_size, 0:this%block_size, 0:this%block_size))
            this%grid(material_type, i_new, j_new ,k_new)%matrix(:,:,:) = 0
        end if

        this%grid(material_type, i_new, j_new ,k_new)%matrix(mod(i, this%block_size), mod(j, this%block_size) ,mod(k, this%block_size)) = val
    end subroutine add_item


    function get_item(this, material_type, i, j, k)
        implicit none
        class(block_csr_t), intent(in) :: this
        integer, intent(in) :: i, j, k, material_type
        integer :: i_new, j_new, k_new
        real(8) :: get_item
        get_item = 0d0

        i_new = i/this%block_size
        j_new = j/this%block_size
        k_new = k/this%block_size

        if (associated(this%grid(material_type, i_new, j_new ,k_new)%matrix)) then 
            ! print*, 'block get'
            get_item = this%grid(material_type, i_new, j_new ,k_new)%matrix(mod(i, this%block_size), mod(j, this%block_size) ,mod(k, this%block_size))
        end if

    end function get_item


    subroutine reorder(this, update_mapper)
        implicit none
        class(block_csr_t), intent(inout) :: this
        logical, intent(in) :: update_mapper
    end subroutine reorder

    subroutine print_data(this, file_name)
        class(block_csr_t), intent(inout) :: this
        character(len=*), intent(in)        ::   file_name
        integer :: i_new, j_new, k_new
        type(block_t) :: block

        integer :: i,j,k,m
        integer :: unit

        open (unit=414, file=file_name, status = 'replace')  

        ! print*, this%nz, this%ny, this%nx, this%nmats
        
        do k = 0, this%nz
            do j = 0, this%ny
                do i = 0, this%nx
                    do m = 1, this%nmats

                        write(414,*)  m, i, j, k, this%get_item(m,i,j,k)
                        
                    end do
                end do
            end do
        end do
        
        close (414)
        
        ! integer :: i,j,k
        ! real :: result, total
        ! total = 0

        !  do k = 0, this%nz
        !     do j = 0, this%ny
        !         do i = 0, this%nx
        !             do m = 1, this%nmats

        !                 total = total + 1
        !                 i_new = i/this%block_size
        !                 j_new = j/this%block_size
        !                 k_new = k/this%block_size

        !                 if (associated(this%grid(m, i_new, j_new ,k_new)%matrix))  result = result + 1
                                        
        !             end do
        !         end do
        !     end do
        ! end do

        ! print*, 'block csr ratio', result/total
   end subroutine


   subroutine deallocate_data(this)
        class(block_csr_t), intent(inout) :: this
        integer :: i,j,k,m
        
        do k = 0, this%d3_block
            do j = 0, this%d2_block
                do i = 0, this%d1_block
                    do m = 1, this%nmats

                        if (associated(this%grid(m, i, j ,k)%matrix)) deallocate(this%grid(m, i, j ,k)%matrix)
                        
                    end do
                end do
            end do
        end do
    

      if (allocated(this%grid))   deallocate(this%grid)
   end subroutine deallocate_data



    subroutine who_am_i(this)
      class(block_csr_t), intent(inout) :: this

      print*, 'ITS A ME: BLOCK CSR'

    end subroutine who_am_i

    subroutine Ptr_coordinates_4d_csr(this, ptr)
        ! class(block_csr_t), intent(in out) :: this
        ! type(block_t), dimension(:,:,:,:), pointer, intent(out) :: ptr

        ! ptr => this%grid

        ! TODO: implement
    end subroutine Ptr_coordinates_4d_csr



    function Get_copy (this)
        class(block_csr_t)       , intent(in)  :: this
        real(8), dimension(:,:,:,:), pointer :: Get_copy

        ! Get_copy = this%values
    end function Get_copy

    subroutine Clean_data_imp (this)
        class(block_csr_t), intent(in out) :: this

        ! TODO: implement
    end subroutine Clean_data_imp




    subroutine Get_recv_buf_imp(this)
        class(block_csr_t), intent(in out) :: this
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
        this%values(1:m, d1+1, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(0 : m * ((d2+2)*(d3+2))-1), (/m, d2+2, d3+2/))
        end if

        !from left
        if (x-1 /= 0) then
        this%values(1:m, 0, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(m * ((d2+2)*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2))-1), (/m, d2+2, d3+2/))
        end if

        !from front
        if (y+1 /= this%parallel_params%npy+1) then
        this%values(1:m, 0:d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+(d1+2)*(d3+2))-1), (/m, d1+2, d3+2/))
        end if

        !from rear
        if (y-1 /= 0) then
        this%values(1:m, 0:d1+1, 0, 0:d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+(d1+2)*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2))-1), (/m, d1+2, d3+2/))
        end if

        !from up
        if (z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, 0:d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2))-1), (/m, d1+2, d2+2/))
        end if

        !from down
        if (z-1 /= 0) then
        this%values(1:m, 0:d1+1, 0:d2+1, 0) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2))-1), (/m, d1+2, d2+2/))
        end if

        !from diagonal of i+1, j+1
        if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
        this%values(1:m, d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2))-1), (/m, d3+2/))
        end if

        !from diagonal of i+1, j-1
        if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
        this%values(1:m, d1+1, 0, 0:d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2))-1), (/m, d3+2/))
        end if

        !from diagonal of i-1, j+1
        if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
        this%values(1:m, 0, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2))-1), (/m, d3+2/))
        end if

        !from diagonal of i-1, j-1
        if (x-1 /= 0 .and. y-1 /= 0) then
        this%values(1:m, 0, 0, 0:d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2))-1), (/m, d3+2/))
        end if

        !from diagonal of i+1, k+1
        if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2))-1), (/m, d2+2/))
        end if

        !from diagonal of i+1, k-1
        if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
        this%values(1:m, d1+1, 0:d2+1, 0) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2))-1), (/m, d2+2/))
        end if

        !from diagonal of i-1, k+1
        if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, 0, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2))-1), (/m, d2+2/))
        end if

        !from diagonal of i-1, k-1
        if (x-1 /= 0 .and. z-1 /= 0) then
        this%values(1:m, 0, 0:d2+1, 0) = reshape(this%recv_buf(m* (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2))-1), (/m, d2+2/))
        end if

        !from diagonal of j+1, k+1
        if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, 0:d1+1, d2+1, d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2))-1), (/m, d1+2/))
        end if

        !from diagonal of j+1, k-1
        if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
        this%values(1:m, 0:d1+1, d2+1, 0) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2)) :&
                    m* (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2))-1), (/m, d1+2/))
        end if

        !from diagonal of j-1, k+1
        if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, 0:d1+1, 0, d3+1) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2))-1), (/m, d1+2/))
        end if

        !from diagonal of j-1, k-1
        if (y-1 /= 0 .and. z-1 /= 0) then
        this%values(1:m, 0:d1+1, 0, 0) = reshape(this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2))-1), (/m, d1+2/))
        end if

        !from diagonal of i+1, j+1, k+1
        if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, d1+1, d2+1, d3+1) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+1)-1)
        end if

        !from diagonal of i+1, j+1, k-1
        if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
        this%values(1:m, d1+1, d2+1, 0) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+1) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+2)-1)
        end if

        !from diagonal of i+1, j-1, k+1
        if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, d1+1, 0, d3+1) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+2) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+3)-1)
        end if

        !from diagonal of i+1, j-1, k-1
        if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
        this%values(1:m, d1+1, 0, 0) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+3) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+4)-1)
        end if

        !from diagonal of i-1, j+1, k+1
        if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, 0, d2+1, d3+1) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+4) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+5)-1)
        end if

        !from diagonal of i-1, j+1, k-1
        if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
        this%values(1:m, 0, d2+1, 0) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+5) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+6)-1)
        end if

        !from diagonal of i-1, j-1, k+1
        if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%values(1:m, 0, 0, d3+1) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+6) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7)-1)
        end if

        !from diagonal of i-1, j-1, k-1
        if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
        this%values(1:m, 0, 0, 0) = this%recv_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7) :&
                    m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+8)-1)
        end if
    end subroutine Get_recv_buf_imp


    subroutine Set_send_buf_imp(this)
        class(block_csr_t), intent(in out) :: this
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
        this%send_buf(0 : m * ((d2+2)*(d3+2))-1) = &
                        reshape(this%values(1:m, d1 - offset, 0:d2+1, 0:d3+1), (/m * (d2+2) * (d3+2)/))
        end if

        !to left
        if (x-1 /= 0) then
        this%send_buf(m * ((d2+2)*(d3+2)) : m * (2*(d2+2)*(d3+2))-1) = &
                        reshape(this%values(1:m, 1 + offset, 0:d2+1, 0:d3+1), (/m * (d2+2) * (d3+2)/))
        end if

        !to front
        if (y+1 /= this%parallel_params%npy+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)) : m * (2*(d2+2)*(d3+2)+(d1+2)*(d3+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, d2 - offset, 0:d3+1), (/m * (d1+2) * (d3+2)/))
        end if

        !to rear
        if (y-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+(d1+2)*(d3+2)) : m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 1 + offset, 0:d3+1), (/m * (d1+2) * (d3+2)/))
        end if

        !to up
        if (z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)) : m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 0:d2+1, d3 - offset), (/m * (d1+2) * (d2+2)/))
        end if

        !to down
        if (z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 0:d2+1, 1 + offset), (/m * (d1+2) * (d2+2)/))
        end if

        !to diagonal of i+1, j+1
        if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2))-1) = &
                        reshape(this%values(1:m, d1 - offset, d2 - offset, 0:d3+1), (/m * (d3+2)/))
        end if

        !to diagonal of i+1, j-1
        if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2))-1) = &
                        reshape(this%values(1:m, d1 - offset, 1 + offset, 0:d3+1), (/m * (d3+2)/))
        end if

        !to diagonal of i-1, j+1
        if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2))-1) = &
                        reshape(this%values(1:m, 1 + offset, d2 - offset, 0:d3+1), (/m * (d3+2)/))
        end if

        !to diagonal of i-1, j-1
        if (x-1 /= 0 .and. y-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2))-1) = &
                        reshape(this%values(1:m, 1+offset, 1+offset, 0:d3+1), (/m * (d3+2)/))
        end if

        !to diagonal of i+1, k+1
        if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2))-1) = &
                        reshape(this%values(1:m, d1 - offset, 0:d2+1, d3 - offset), (/m * (d2+2)/))
        end if

        !to diagonal of i+1, k-1
        if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2))-1) = &
                        reshape(this%values(1:m, d1 - offset, 0:d2+1, 1 + offset), (/m * (d2+2)/))
        end if

        !to diagonal of i-1, k+1
        if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2))-1) = &
                        reshape(this%values(1:m, 1 + offset, 0:d2+1, d3 - offset), (/m * (d2+2)/))
        end if

        !to diagonal of i-1, k-1
        if (x-1 /= 0 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2))-1) = &
                        reshape(this%values(1:m, 1 + offset, 0:d2+1, 1 + offset), (/m * (d2+2)/))
        end if

        !to diagonal of j+1, k+1
        if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)) : &
                        m* (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, d2 - offset, d3 - offset), (/m * (d1+2)/))
        end if

        !to diagonal of j+1, k-1
        if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, d2 - offset, 1 + offset), (/m * (d1+2)/))
        end if

        !to diagonal of j-1, k+1
        if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 1 + offset, d3 - offset), (/m * (d1+2)/))
        end if

        !to diagonal of j-1, k-1
        if (y-1 /= 0 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2)) : &
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2))-1) = &
                        reshape(this%values(1:m, 0:d1+1, 1 + offset, 1+offset), (/m * (d1+2)/))
        end if


        !to diagonal of i+1, j+1, k+1
        if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+1)-1) = &
                        this%values(1:m, d1- offset, d2-offset, d3-offset)
        end if

        !to diagonal of i+1, j+1, k-1
        if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+1) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+2)-1) = &
                        this%values(1:m, d1-offset, d2-offset, 1 + offset)
        end if

        !to diagonal of i+1, j-1, k+1
        if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+2) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+3)-1) = &
                        this%values(1:m, d1-offset, 1+offset, d3-offset)
        end if

        !to diagonal of i+1, j-1, k-1
        if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+3) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+4)-1) = &
                        this%values(1:m, d1-offset, 1+offset, 1+offset)
        end if

        !to diagonal of i-1, j+1, k+1
        if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+4) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+5)-1) = &
                        this%values(1:m, 1+offset, d2-offset, d3-offset)
        end if

        !to diagonal of i-1, j+1, k-1
        if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+5) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+6)-1) = &
                        this%values(1:m, 1+offset, d2-offset, 1+offset)
        end if

        !to diagonal of i-1, j-1, k+1
        if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+6) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7)-1) = &
                        this%values(1:m, 1+offset, 1+offset, d3-offset)
        end if

        !to diagonal of i-1, j-1, k-1
        if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
        this%send_buf(m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7) :&
                        m * (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+8)-1) = &
                        this%values(1:m, 1+offset, 1+offset, 1+offset)
        end if
    end subroutine Set_send_buf_imp


    subroutine Debug_check_nan(this, caller)
        implicit none
        class(block_csr_t), intent(in out) :: this
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
        class(block_csr_t), intent(in out) :: this
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
        class(block_csr_t), intent(in) :: this
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

        write(unit, iostat=iostat, iomsg=iomsg) &
        shape(this%values)

        write(unit, iostat=iostat, iomsg=iomsg) &
        this%values, &
        this%nx, &
        this%ny, &
        this%nz


    end subroutine Write_data

    subroutine Read_data(this, unit, iostat, iomsg)
        class(block_csr_t), intent(in out) :: this
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        integer, dimension(3) :: shape_values

        read(unit, iostat=iostat, iomsg=iomsg) &
        shape_values

        deallocate(this%values)
    !      allocate(this%values(0:shape_values(1) - 1, 0:shape_values(2)-1, 0:shape_values(3) - 1))

        read(unit, iostat=iostat, iomsg=iomsg) &
        this%values, &
        this%nx, &
        this%ny, &
        this%nz

    end subroutine Read_data

end module