module indexer_module

    type :: indexer_t
        integer, dimension(:,:,:,:), allocatable   ::   mapper

        logical :: initiated = .False.
        integer :: m, nx, ny, nz

        contains
            procedure, nopass, private :: constructor    
    end type

    type(indexer_t), allocatable, target, private :: indexer
    

    contains


    function constructor(m, nx, ny, nz)
        implicit none
        integer, intent(in)            :: m, nx, ny, nz
        class(indexer_t), allocatable :: constructor
        
        allocate(constructor)
        allocate(constructor%mapper(1:m,0:nx,0:ny,0:nz))

    end function constructor

    
    function get_instance(m, nx, ny, nz)
        implicit none
        class(indexer_t), pointer :: get_instance
        integer, optional, intent(in) :: m, nx, ny, nz

        if (.not. allocated(indexer) .and. present(m) .and. present(nx) .and. present(ny) .and. present(nz)) then
            print*, 'bbbb'
            indexer = constructor(m, nx, ny, nz)
            get_instance => indexer
        else
            print*, 'aaaaa'
            get_instance => indexer
        end if
    end function get_instance


end module indexer_module


