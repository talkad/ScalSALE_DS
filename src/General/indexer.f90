module indexer_module

    type :: indexer_t
        integer, dimension(:,:,:,:), allocatable   ::   mapper

        logical :: initiated = .False.
        integer :: m, nx, ny, nz
        integer :: last_idx

        contains
            procedure, nopass, private :: mapper_constructor    
    end type

    type(indexer_t), allocatable, target, private :: indexer
    

    contains


    function mapper_constructor(m, nx, ny, nz)
        implicit none
        integer, intent(in)            :: m, nx, ny, nz
        type(indexer_t), target, allocatable :: mapper_constructor
        
        allocate(mapper_constructor)
        allocate(mapper_constructor%mapper(1:m,0:nx,0:ny,0:nz))

        mapper_constructor%mapper(1:m,0:nx,0:ny,0:nz) = -1
        mapper_constructor%last_idx = 0

    end function mapper_constructor

    
    function get_instance(m, nx, ny, nz)
        implicit none
        type(indexer_t), pointer :: get_instance
        integer, optional, intent(in) :: m, nx, ny, nz

        ! print*, 'get_instance'
        if (.not. allocated(indexer) .and. present(m) .and. present(nx) .and. present(ny) .and. present(nz)) then
            indexer = mapper_constructor(m, nx, ny, nz)
            get_instance => indexer
        else
            get_instance => indexer
        end if

        ! print*, 'mapper shape', shape(indexer%mapper)
    end function get_instance


end module indexer_module


