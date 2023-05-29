module sparse_struct_base_module

    type, abstract :: sparse_struct_base_t

    contains
        procedure(get), public, deferred :: get_item
        procedure(update_data), public, deferred :: batch_update
    end type

    abstract interface
        pure function get(this, material_type, i, j, k)
            import sparse_struct_base_t
            class(sparse_struct_base_t), intent(in) :: this
            integer, intent(in) :: i, j, k, material_type
            real(8) :: get
        end function

        subroutine batch_update(this, ms, is, js, ks, vals)
            import sparse_struct_base_t
            class(sparse_struct_base_t), intent(inout) :: this
            integer, dimension(:), allocatable, intent(in) :: is, js, ks, ms
            real(8), dimension(:), allocatable, intent(in) :: vals
        end subroutine

    end interface


end module
