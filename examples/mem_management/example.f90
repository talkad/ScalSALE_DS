program DynamicArrays
  implicit none

  integer, parameter :: outerSize = 100000
  integer, parameter :: innerSize = 20000
  integer, dimension(:), allocatable :: outerArray(:)
  integer, dimension(:), pointer :: innerArray(:)
  integer :: i, j, sum

  allocate(outerArray(outerSize))

  do i = 1, outerSize
    allocate(innerArray(innerSize))
    innerArray = 1 
    outerArray(i) = transfer(innerArray, outerArray(i)) 
  end do

  sum = 0
  do i = 1, outerSize
    innerArray = transfer(outerArray(i), innerArray) 
    do j = 1, innerSize
      sum = sum + innerArray(j)
    end do
  end do

  print *, "Sum of all values in the inner arrays: ", sum


end program DynamicArrays
