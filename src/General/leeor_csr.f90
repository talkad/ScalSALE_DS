module leeor_csr_module
    use data_struct_base


    type :: leeor_materials_t
        integer :: material_idx
        real(8) :: material_val
        real(8), dimension(:), pointer :: materials
    end type leeor_materials_t


    type, extends(data_struct_t) :: leeor_csr_t

        real(8), dimension(:,:,:,:), pointer, public         :: values

        type(leeor_materials_t), dimension(:,:,:), allocatable :: grid

    contains


        procedure, public :: Ptr_coordinates_4d => Ptr_coordinates_4d_leeor

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


    interface leeor_csr_t
        module procedure Constructor_init_val
    end interface


    contains
    

    function Constructor_init_val(initial_val, d1, d2, d3, d4)
        implicit none

        type(leeor_csr_t), pointer :: Constructor_init_val
        real(8)           , intent(in) :: initial_val  
        integer           , intent(in) :: d1, d2, d3, d4  
        integer                        :: k,j,i,m

        print*, 'init_leeor_csr'

        allocate(Constructor_init_val)
        Constructor_init_val%communication => null()

        allocate(Constructor_init_val%values(0:1,0:1,0:1,0:1))
        Constructor_init_val%values = 0

        Constructor_init_val%nx = d1
        Constructor_init_val%ny = d2
        Constructor_init_val%nz = d3
        Constructor_init_val%nmats = d4

! 1:d4, 
        allocate(Constructor_init_val%grid(0:d1, 0:d2, 0:d3)) 

        do k = 0, d3
            do j = 0, d2
                do i = 0, d1
                    Constructor_init_val%grid(i, j ,k)%materials => null()
                    Constructor_init_val%grid(i, j ,k)%material_idx = -1
                end do
            end do
        end do        

        call Constructor_init_val%add_boundary()

    end function


    subroutine Ptr_coordinates_4d_leeor(this, ptr)
        class (leeor_csr_t)                    , intent(in out)  :: this
        real(8), dimension(:,:,:,:), pointer, intent(out) :: ptr

        ptr => this%values
    end subroutine Ptr_coordinates_4d_leeor


    subroutine add_boundary(this)
        implicit none
        class(leeor_csr_t), intent(inout) :: this
        integer :: i,j,k,m
                
        i = 0
        do k = 0, this%nz
            do j = 0, this%ny
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0)
                end do
            end do
        end do


        i = this%nx
        do k = 0, this%nz
            do j = 0, this%ny
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0)
                end do
            end do
        end do


        j = 0
        do k = 0, this%nz
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0)
                end do
            end do
        end do


        j = this%ny
        do k = 0, this%nz
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0)
                end do
            end do
        end do


        k = 0
        do j = 0, this%ny
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0)
                end do
            end do
        end do

        k = this%nz
        do j = 0, this%ny
            do i = 0, this%nx
                do m = 1, this%nmats
                    call this%add_item(m,i,j,k, 0d0)
                end do
            end do
        end do

    
    end subroutine add_boundary


    subroutine add_item(this, material_type, i, j, k, val, boundry)
        implicit none
        class(leeor_csr_t), intent(inout) :: this
        integer, intent(in) :: i, j, k, material_type
        real(8), intent(in) :: val
        logical, optional, intent(in) :: boundry

        if (this%grid(i,j,k)%material_idx == -1 .and. val == 0) return 

        if (this%grid(i,j,k)%material_idx == -1) then
            this%grid(i,j,k)%material_idx = material_type
            this%grid(i,j,k)%material_val = val
        else
            if (.not. associated(this%grid(i,j,k)%materials) .and. val==0) return

            if (.not. associated(this%grid(i,j,k)%materials) .and. this%grid(i,j,k)%material_idx == material_type)  then
                this%grid(i,j,k)%material_val = val
            else
                if (.not. associated(this%grid(i,j,k)%materials)) then
                    allocate(this%grid(i,j,k)%materials(1:this%nmats))
                    this%grid(i,j,k)%materials = 0d0
                    this%grid(i,j,k)%materials(this%grid(i,j,k)%material_idx) = this%grid(i,j,k)%material_val
                end if

                this%grid(i,j,k)%materials(material_type) = val
            end if


        end if

    end subroutine add_item


    function get_item(this, material_type, i, j, k)
        implicit none
        class(leeor_csr_t), intent(in) :: this
        integer, intent(in) :: i, j, k, material_type
        real(8) :: get_item
        get_item = 0d0

        if (associated(this%grid(i,j,k)%materials)) then
            get_item = this%grid(i,j,k)%materials(material_type)
        else if (this%grid(i,j,k)%material_idx == material_type) then
            get_item = this%grid(i,j,k)%material_val
        end if

    end function get_item

    subroutine reorder(this, update_mapper)
        implicit none
        class(leeor_csr_t), intent(inout) :: this
        logical, intent(in) :: update_mapper
    end subroutine reorder

    subroutine print_data(this, file_name)
        class(leeor_csr_t), intent(inout) :: this
        character(len=*), intent(in)        ::   file_name
        integer :: i_new, j_new, k_new

        integer :: i,j,k,m
        integer :: unit

        open (unit=414, file=file_name, status = 'replace')  
        
        do k = 0, this%nz
            do j = 0, this%ny
                do i = 0, this%nx
                    do m = 1, this%nmats
                        write(414,*) m, i, j, k, this%get_item(m, i, j, k)                      
                    end do
                end do
            end do
        end do
        
        close (414)
   end subroutine


   subroutine deallocate_data(this)
        class(leeor_csr_t), intent(inout) :: this
        integer :: i,j,k,m
        
        do k = 0, this%nz
            do j = 0, this%ny
                do i = 0, this%nz
                    if (associated(this%grid(i, j ,k)%materials)) deallocate(this%grid(i, j ,k)%materials)    
                end do
            end do
        end do
    

      if (allocated(this%grid))   deallocate(this%grid)
   end subroutine deallocate_data



    subroutine who_am_i(this)
      class(leeor_csr_t), intent(inout) :: this

      print*, 'ITS A ME: LEE-OR CSR'

    end subroutine who_am_i

    subroutine Ptr_coordinates_4d_csr(this, ptr)
        ! class(leeor_csr_t), intent(in out) :: this
        ! type(leeor_csr_t), dimension(:,:,:,:), pointer, intent(out) :: ptr

        ! ptr => this%grid

        ! TODO: implement
    end subroutine Ptr_coordinates_4d_csr



    function Get_copy (this)
        class(leeor_csr_t)       , intent(in)  :: this
        real(8), dimension(:,:,:,:), pointer :: Get_copy

        ! Get_copy = this%values
    end function Get_copy

    subroutine Clean_data_imp (this)
        class(leeor_csr_t), intent(in out) :: this

        ! TODO: implement
    end subroutine Clean_data_imp




    subroutine Get_recv_buf_imp(this)
        class(leeor_csr_t), intent(in out) :: this
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
        class(leeor_csr_t), intent(in out) :: this
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
        class(leeor_csr_t), intent(in out) :: this
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
        class(leeor_csr_t), intent(in out) :: this
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
        class(leeor_csr_t), intent(in) :: this
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
        class(leeor_csr_t), intent(in out) :: this
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