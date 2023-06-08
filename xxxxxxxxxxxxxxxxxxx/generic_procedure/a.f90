module math_operations
    implicit none
    
    interface generic_operation
        module procedure generic_add, generic_multiply
    end interface generic_operation
    
    contains
    
    ! Specific procedure for adding integers
    function generic_add(a, b) result(result)
        integer, intent(in) :: a, b
        integer :: result
        result = a + b
    end function generic_add
    
    ! Specific procedure for multiplying real numbers
    function generic_multiply(a, b) result(result)
        real, intent(in) :: a, b
        real :: result
        result = a * b
    end function generic_multiply
    
end module math_operations

program main
    use math_operations
    implicit none
    
    integer :: x, y, sum
    real :: a, b, product
    
    x = 5
    y = 10
    a = 2.5
    b = 3.5
    
    sum = generic_operation(x, y)  ! Invoking generic_add
    product = generic_operation(a, b)  ! Invoking generic_multiply
    
    print*, "Sum:", sum
    print*, "Product:", product
    
end program main