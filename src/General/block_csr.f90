module csr_module
    use indexer_module
    use data_struct_base

    type, extends(data_struct_t) :: csr_t

        real(8), dimension(:,:,:,:), pointer, public         :: values

        !! CSR Arguments 
        real(8), dimension(:), pointer :: nz_values
        integer, dimension(:,:,:,:), pointer :: idx_map        
        integer :: padding_size
        integer :: padding_idx

        real(8) :: ratio = 0.01

    contains

        ! generic, public :: Point_to_data => &
        !                     Ptr_coordinates_1d, &
        !                     Ptr_coordinates_4d

        procedure, public :: Ptr_coordinates_1d => Ptr_coordinates_1d_csr
        ! procedure, public :: Ptr_coordinates_4d => Ptr_coordinates_4d_csr




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
        procedure :: update_struct => update_struct
        procedure, private, nopass :: append_real
        procedure, private :: count_new_cells
        procedure, private :: add_boundary

    end type

    public :: Get_copy


    interface csr_t
        module procedure Constructor_init_val
    end interface


    contains
    

    type(csr_t) function Constructor_init_val(initial_val, d1, d2, d3, d4, idx_map)
        implicit none
        integer, dimension(:,:,:,:), allocatable, target, intent(inout) :: idx_map

        real(8)           , intent(in) :: initial_val  
        integer           , intent(in) :: d1, d2, d3, d4  
        integer :: space_size
        integer :: i, j, k, m
        print*, 'init_csr'
        allocate(Constructor_init_val%values(0:1,0:1,0:1,0:1))
        Constructor_init_val%values = 0


        Constructor_init_val%nx = d1
        Constructor_init_val%ny = d2
        Constructor_init_val%nz = d3
        Constructor_init_val%nmats = d4

        ! print*, 'aaaaaaaaaa', d4, d3, d2, d1

        Constructor_init_val%idx_map => idx_map
        space_size = (d1+1) * (d2+1) * (d3+1) * d4 !* 4

        ! print*, 'nz_size', space_size

        Constructor_init_val%padding_size = space_size*Constructor_init_val%ratio
        Constructor_init_val%padding_idx = 0

        allocate(Constructor_init_val%nz_values(0:space_size+int(space_size*Constructor_init_val%ratio)))
        Constructor_init_val%nz_values = 0d0
        ! print*, 'aaaaaaaaaaaaaaaaaaaa',d4,d1,d2,d3

        call Constructor_init_val%add_boundary()

    end function


    subroutine add_boundary(this)
        implicit none
        class(csr_t), intent(inout) :: this
        type(indexer_t), pointer ::  index_mapper
        integer :: i,j,k,m
        integer :: csr_idx
        
        index_mapper => get_instance()
        ! print*, 'aaaaaaaaaa', shape(this%idx_map)
        ! print*, 'bbbbbbbbbb', this%nmats, this%nx, this%ny, this%nz 
        ! print*, 'aaaaaaaaaaaaaa', this%nx, this%ny, this%nz, this%nmats
        csr_idx = 0

        i = 0
        do k = 0, this%nz
            do j = 0, this%ny
                do m = 1, this%nmats
                    this%idx_map(m,i,j,k) = csr_idx
                    this%nz_values(csr_idx) = 0

                    csr_idx = csr_idx + 1
                end do
            end do
        end do


        i = this%nx
        do k = 0, this%nz
            do j = 0, this%ny
                do m = 1, this%nmats
                    this%idx_map(m,i,j,k) = csr_idx
                    this%nz_values(csr_idx) = 0

                    csr_idx = csr_idx + 1
                end do
            end do
        end do


        j = 0
        do k = 0, this%nz
            do i = 0, this%nx
                do m = 1, this%nmats
                    this%idx_map(m,i,j,k) = csr_idx
                    this%nz_values(csr_idx) = 0

                    csr_idx = csr_idx + 1
                end do
            end do
        end do


        j = this%ny
        do k = 0, this%nz
            do i = 0, this%nx
                do m = 1, this%nmats
                    this%idx_map(m,i,j,k) = csr_idx
                    this%nz_values(csr_idx) = 0

                    csr_idx = csr_idx + 1
                end do
            end do
        end do


        k = 0
        do j = 0, this%ny
            do i = 0, this%nx
                do m = 1, this%nmats
                    this%idx_map(m,i,j,k) = csr_idx
                    this%nz_values(csr_idx) = 0

                    csr_idx = csr_idx + 1
                end do
            end do
        end do

        k = this%nz
        do j = 0, this%ny
            do i = 0, this%nx
                do m = 1, this%nmats
                    this%idx_map(m,i,j,k) = csr_idx
                    this%nz_values(csr_idx) = 0

                    csr_idx = csr_idx + 1
                end do
            end do
        end do

        index_mapper%last_idx = csr_idx
    
    end subroutine add_boundary

    ! subroutine print_mapper(file_name)

    !     type(indexer_t), pointer ::  index_mapper
    !     integer, dimension(:,:,:,:), pointer   ::   mapper
    !     integer :: index, m,i,j,k
    !     character(len=*), intent(in)                      ::   file_name

    !     index_mapper => get_instance()
    !     mapper => index_mapper%mapper

    !     open (unit=414, file=file_name, status = 'replace')  
        
    !     print*, shape(mapper)

    !     do k = 0, size(mapper, dim=4)-1
    !         do j = 0, size(mapper, dim=3)-1
    !             do i = 0, size(mapper, dim=2)-1
    !                 do m = 1, size(mapper, dim=1)
    !                     index = mapper(m,i,j,k) 
    !                     write(414,*) m, i , j , k , index
    !                 end do
    !             end do
    !         end do
    !     end do
        
    !     close (414)

    ! end subroutine print_mapper

    pure subroutine add_item(this, material_type, i, j, k, val)
        implicit none
        class(csr_t), intent(inout) :: this
        integer, intent(in) :: i, j, k, material_type
        real(8), intent(in) :: val
        integer :: idx

        idx = this%idx_map(material_type, i, j, k)
        if (idx > -1) this%nz_values(idx) = val
    end subroutine add_item


    ! Assume the order of the new values is j i m
    subroutine update_struct(this, ms, is, js, ks, vals)
        implicit none
        class(csr_t), intent(inout) :: this
        integer, dimension(:), allocatable, intent(in) :: is, js, ks, ms
        real(8), dimension(:), allocatable, intent(in) :: vals
        integer :: m, i, j, k, materials, nx, ny, nz, insertion_idx, idx, num_vals
        integer :: num_padding, new_cells, num_pads, prev_size, new_size
        real(8), dimension(:), allocatable :: temp
        integer :: current_idx
        logical :: logic_debug
        idx = 0

        ! logic_debug = .False.
        ! num_padding = size(this%nz_values) - this%padding_idx
        ! new_cells = this%count_new_cells(ms,is,js,ks)

        ! if (new_cells > num_padding) then
        !     ! rescalse the size of the values array 
        !     num_pads = (new_cells - num_padding) / this%padding_size + 1 

        !     prev_size = size(this%nz_values, dim=1)

        !     new_size = prev_size + num_pads * this%padding_size - 1
        !     ! write(*,*) 'reallocation', prev_size, '->', new_size

        !     allocate(temp(0:new_size))                          ! enlarge array size by factor of 2
        !     temp(0:new_size) = 0d0
        !     temp(0:prev_size-1) = this%nz_values(0:prev_size-1)    ! copy previous values
        !     call move_alloc(temp, this%nz_values)                  ! temp gets deallocated
        ! end if

        ! materials = size(this%idx_map, dim=1)
        ! nx = size(this%idx_map, dim=2)-1
        ! ny = size(this%idx_map, dim=3)-1
        ! nz = size(this%idx_map, dim=4)-1
        ! num_vals = size(is)-1
        
        ! do i=num_vals, 0, -1
        !     if (is(i)/=-1) then
        !         idx = i
        !         exit
        !     end if
        ! end do
        
        ! insertion_idx = size(this%nz_values)-1
        
        ! do k=nz, 0, -1
        !     do j=ny, 0, -1
        !         do i=nx, 0, -1
        !             do m=materials, 1, -1

        !                 current_idx = this%idx_map(m,i,j,k)

        !                 if (ms(idx)==m .and. is(idx)==i .and. js(idx)==j .and. ks(idx)==k) then
        !                     this%nz_values(insertion_idx) = vals(idx)
        !                     this%idx_map(m,i,j,k) = insertion_idx
        !                     insertion_idx = insertion_idx - 1
        !                     idx = idx - 1
        !                 else if (current_idx > -1) then
        !                     this%nz_values(insertion_idx) = this%nz_values(current_idx)
        !                     this%nz_values(current_idx) = 0
        !                     this%idx_map(m,i,j,k) = insertion_idx
        !                     insertion_idx = insertion_idx - 1
        !                 end if 

        !             end do
        !         end do
        !     end do
        ! end do 

        ! idx = 0

        ! do k=0, nz
        !     do j=0, ny
        !         do i=0, nx
        !             do m=1, materials

        !                 current_idx = this%idx_map(m,i,j,k) 

        !                 if (current_idx > -1) then
        !                     this%nz_values(idx) = this%nz_values(current_idx) 
        !                     this%nz_values(current_idx) = 0
        !                     this%idx_map(m,i,j,k) = idx
        !                     idx = idx + 1
        !                 end if

        !             end do
        !         end do
        !     end do
        ! end do

        ! this%padding_idx = idx

    end subroutine update_struct


    function count_new_cells(this, ms, is, js, ks) result(num)
        implicit none
        class(csr_t), intent(inout) :: this
        integer, dimension(:), allocatable, intent(in) :: is, js, ms, ks
        integer :: num, idx

        num = 0

        do idx=0, size(is)-1
            if (this%idx_map(ms(idx), is(idx), js(idx), ks(idx)) == -1) num=num+1
        end do
    end function count_new_cells


    pure function get_item(this, material_type, i, j, k)
        implicit none
        class(csr_t), intent(in) :: this
        integer, intent(in) :: i, j, k, material_type
        real(8) :: get_item
        integer :: idx
        get_item = 0d0
        idx = this%idx_map(material_type, i, j, k)
        
        if (idx > -1) get_item = this%nz_values(idx)
    end function get_item


    ! helper function - append val to array if it is not full, otherwise enlarge the array and append
    subroutine append_real(array, val, idx)
        implicit none
        real(8), dimension(:), allocatable, intent(inout) :: array
        real(8), intent(in) :: val
        integer, intent(in) :: idx

        integer :: prev_size, new_size
        real(8), dimension(:), allocatable :: temp

        prev_size = size(array, dim=1)

        if (prev_size <= idx) then
            new_size = prev_size*2 - 1

            allocate(temp(0:new_size))                   ! enlarge array size by factor of 2
            temp(0:new_size) = 0d0
            temp(0:prev_size-1) = array(0:prev_size-1)   ! copy previous values
            call move_alloc(temp, array)                 ! temp gets deallocated
        end if

        array(idx) = val
    end subroutine append_real




    subroutine Ptr_coordinates_1d_csr (this, ptr)
        class (csr_t)                    , intent(in out)  :: this
        real(8), dimension(:), pointer, intent(out) :: ptr

        ptr => this%nz_values
    end subroutine Ptr_coordinates_1d_csr


    ! subroutine Ptr_coordinates_4d_csr(this, ptr)
    !     class (csr_t), intent(in out) :: this
    !     real(8), dimension(:,:,:,:), pointer, intent(out) :: ptr
    ! end subroutine Ptr_coordinates_4d_csr



    function Get_copy (this)
        class (csr_t)       , intent(in)  :: this
        real(8), dimension(:,:,:,:), pointer :: Get_copy

        ! Get_copy = this%values
    end function Get_copy

    subroutine Clean_data_imp (this)
        class (csr_t), intent(in out) :: this

        deallocate (this%nz_values)
    end subroutine Clean_data_imp




    subroutine Get_recv_buf_imp(this)
        class (csr_t), intent(in out) :: this
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
        class (csr_t), intent(in out) :: this
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
        class (csr_t), intent(in out) :: this
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
        class (csr_t), intent(in out) :: this
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
        class (csr_t), intent(in) :: this
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
        class (csr_t), intent(in out) :: this
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