module indexer_module

    type :: indexer_t

        type(indexer_t), allocatable, private :: indexer
        integer, dimension(:,:,:,:), allocatable   ::   mapper

        logical :: initiated = .False.
        integer :: m, nx, ny, nz

        contains
            procedure, public :: set_dim
            procedure, pass, private :: constructor    
            procedure, public :: get_instance    
    end type

    contains


    function constructor(m, nx, ny, nz)
        implicit none
        integer, intent(in)            :: m, nx, ny, nz
        type(indexer_t), allocatable :: constructor
        
        allocate(constructor)
        allocate(constructor%mapper(1:m,0:nx,0:ny,0:nz))

    end function constructor


    function get_instance(this)
        implicit none
        type(indexer_t), intent(inout) :: this
        type(indexer_t), pointer :: get_instance

        if (.not. this%initiated) then
            return
        else if (.not. allocated(this%indexer)) then
            this = constructor(this%m, this%nx, this%ny, this%nz)
        else
            get_instance = this%indexer
        end if
    end function get_instance


    subroutine set_dim(this, m, nx, ny, nz)
        implicit none
        type(indexer_t), intent(inout) :: this
        integer, intent(in)            :: m, nx, ny, nz

        this%m = m
        this%nx = nx
        this%ny = ny
        this%nz = nz

        this%initiated = .True.
    end subroutine set_dim


end module indexer_module