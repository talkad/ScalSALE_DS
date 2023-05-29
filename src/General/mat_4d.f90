module mat4d_module
    use sparse_struct_base_module

    type, extends(sparse_struct_base_t) :: mat4d_t
        real(8), dimension(:,:,:,:), allocatable :: matrix

        contains
            procedure :: get_item => get_item
            procedure :: batch_update => batch_update
    end type
    
    interface mat4d_t
        module procedure mat4d_constructor
    end interface

    contains


    function mat4d_constructor(mats, nx, ny, nz)
        integer, intent(in) :: mats, nx, ny, nz
        type(mat4d_t), pointer :: mat4d_constructor

        allocate(mat4d_constructor)
        allocate(mat4d_constructor%matrix(1:mats,0:nx,0:ny,0:nz))
    end function


    pure function get_item(this, material_type, i, j, k)
        implicit none
        class(mat4d_t), intent(in) :: this
        integer, intent(in) :: material_type, i, j, k
        real(8) :: get_item

        get_item = this%matrix(material_type, i, j, k)
    end function

    subroutine batch_update(this, ms, is, js, ks, vals)
        implicit none
        class(mat4d_t), intent(inout) :: this
        integer, dimension(:), allocatable, intent(in) :: ms, is, js, ks
        real(8), dimension(:), allocatable, intent(in) :: vals
        integer :: idx

        do idx=0, size(is)
            ! assumption cell won't be empty
            ! if (vals(idx) == 0) EXIT
            
            this%matrix(ms(idx), is(idx), js(idx), ks(idx)) = vals(idx)
        end do
    end subroutine


end module

