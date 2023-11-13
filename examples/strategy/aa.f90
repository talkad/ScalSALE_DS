module aa_module

    type, abstract :: values_t
        contains
            procedure (get_values_abstract), public, deferred :: get_values
    end type

    interface
        function get_values_abstract(this)
            import values_t
            class(values_t), intent(inout) :: this
            ! generic :: get_values_abstract
        end function
    end interface

    ! type, extends(values_t) :: values_real_t
    !     real(8), dimension(:,:,:,:), pointer          :: values

    !     contains
    !         procedure, public, pointer :: get_values => get_values_real
    ! end type



    ! type :: data_4d_t
    !     private

    !     class(values_t), pointer :: values
    ! end type

    contains

    ! function get_values_real(this)
    !     class(values_t), pointer :: this
    !     pointer :: get_values_real

    !     ! get_values_real => 
    ! end function


end module