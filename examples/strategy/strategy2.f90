
module strategy_mod
    implicit none


    type, abstract :: data_t
        contains
            procedure(get_values_abstract), deferred :: get_values_abstract
    end type data_t

    interface 
        function get_values_abstract(this)
            import data_t
            class(data_t), intent(in) :: this
            class(*), dimension(:,:,:,:), pointer :: get_values_abstract
        end function get_values_abstract
    end interface

    type, extends(data_t) :: data_real_t
        real(8), dimension(:,:,:,:), pointer :: values

        contains
            procedure :: get_values_abstract => get_values_real
    end type data_real_t

    type, extends(data_t) :: data_integer_t
        integer, dimension(:,:,:,:), pointer :: values

        contains
            procedure :: get_values_abstract => get_values_integer
    end type data_integer_t

contains


    function get_values_real(this)
        class(data_real_t), intent(in) :: this
        class(*), dimension(:,:,:,:), pointer :: get_values_real
        get_values_real => this%values
    end function get_values_real


    function get_values_integer(this)
        class(data_integer_t), intent(in) :: this
        class(*), dimension(:,:,:,:), pointer :: get_values_integer
        get_values_integer => this%values
    end function get_values_integer


end module strategy_mod
