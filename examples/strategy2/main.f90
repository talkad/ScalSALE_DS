
! program main
!     ! use strategy_mod


!     ! class(data_t), allocatable :: data
!     ! allocate(data_real_t :: data)
!     ! values => data%get_values_real()

!     class(*), dimension(:,:,:,:), pointer :: values1
!     real(8), dimension(:,:,:,:), pointer :: values2

!     class(*), dimension(:,:,:,:), pointer :: values11
!     real(8), dimension(:,:,:,:), pointer :: values22

!     allocate(values2(10,10,10,10))
!     allocate(values22(10,10,10,10))


!     values1 => values2
!     values11 => values22

!     values1(1,1,1,1) = values11(1,1,1,1)
!     ! print*, 'aaaa', values2(1,1,1,1)


! end program


integer :: this%pre_calc(14)
vals_shape = shape(this%values)
d1 = vals_shape(2) - 2
d2 = vals_shape(3) - 2
d3 = vals_shape(4) - 2


this%pre_calc(1) = d1+2
this%pre_calc(2) = d3+2
this%pre_calc(3) = this%pre_calc(3) = d2+2
this%pre_calc(4) = this%pre_calc(3)*(this%pre_calc(2))
this%pre_calc(5) = 2*this%pre_calc(4)
this%pre_calc(6) = (this%pre_calc(1))*(this%pre_calc(2))
this%pre_calc(7) = 2*this%pre_calc(6)
this%pre_calc(8) = (this%pre_calc(1))*this%pre_calc(3)
this%pre_calc(9) = 2*this%pre_calc(8)
this%pre_calc(10) = 4*(this%pre_calc(3))+4*(this%pre_calc(1))
this%pre_calc(11) = this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(9)
this%pre_calc(12) = 4*(this%pre_calc(2))
this%pre_calc(13) = this%pre_calc(11)+this%pre_calc(12)
this%pre_calc(14) = this%pre_calc(13)+this%pre_calc(10) 



   subroutine Get_recv_buf(this)
      class (data_4d_t), intent(in out) :: this
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
         this%values(1:m, d1+1, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(0 : m * (this%pre_calc(4))-1), (/m, this%pre_calc(3), this%pre_calc(2)/))
      end if

      !from left
      if (x-1 /= 0) then
         this%values(1:m, 0, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(4)) :&
                     m * (this%pre_calc(5))-1), (/m, this%pre_calc(3), this%pre_calc(2)/))
      end if

      !from front
      if (y+1 /= this%parallel_params%npy+1) then
         this%values(1:m, 0:d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(5)) :&
                     m * (this%pre_calc(5)+this%pre_calc(6)))-1), (/m, this%pre_calc(1), this%pre_calc(2)/))
      end if

      !from rear
      if (y-1 /= 0) then
         this%values(1:m, 0:d1+1, 0, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(5)+this%pre_calc(6)) :&
                     m * (this%pre_calc(5)+this%pre_calc(7))-1), (/m, this%pre_calc(1), this%pre_calc(2)/))
      end if

      !from up
      if (z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0:d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(5)+this%pre_calc(7)) :&
                     m * (this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(8))-1), (/m, this%pre_calc(1), this%pre_calc(3)/))
      end if

      !from down
      if (z-1 /= 0) then
         this%values(1:m, 0:d1+1, 0:d2+1, 0) = reshape(this%recv_buf(m * (this%pre_calc(5)+this%pre_calc(7)+this%pre_calc(8)) :&
                     m * (this%pre_calc(11))-1), (/m, this%pre_calc(1), this%pre_calc(3)/))
      end if

      !from diagonal of i+1, j+1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
         this%values(1:m, d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)) :&
                     m * (this%pre_calc(11)+(this%pre_calc(2)))-1), (/m, this%pre_calc(2)/))
      end if

      !from diagonal of i+1, j-1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
         this%values(1:m, d1+1, 0, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)+(this%pre_calc(2))) :&
                     m * (this%pre_calc(11)+2*(this%pre_calc(2)))-1), (/m, this%pre_calc(2)/))
      end if

      !from diagonal of i-1, j+1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
         this%values(1:m, 0, d2+1, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)+2*(this%pre_calc(2))) :&
                     m * (this%pre_calc(11)+3*(this%pre_calc(2)))-1), (/m, this%pre_calc(2)/))
      end if

      !from diagonal of i-1, j-1
      if (x-1 /= 0 .and. y-1 /= 0) then
         this%values(1:m, 0, 0, 0:d3+1) = reshape(this%recv_buf(m * (this%pre_calc(11)+3*(this%pre_calc(2))) :&
                     m * (this%pre_calc(13))-1), (/m, this%pre_calc(2)/))
      end if

      !from diagonal of i+1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)) :&
                     m * (this%pre_calc(13)+(this%pre_calc(3)))-1), (/m, this%pre_calc(3)/))
      end if

      !from diagonal of i+1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
         this%values(1:m, d1+1, 0:d2+1, 0) = reshape(this%recv_buf(m * (this%pre_calc(13)+(this%pre_calc(3))) :&
                     m * (this%pre_calc(13)+2*(this%pre_calc(3)))-1), (/m, this%pre_calc(3)/))
      end if

      !from diagonal of i-1, k+1
      if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0, 0:d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)+2*(this%pre_calc(3))) :&
                     m * (this%pre_calc(13)+3*(this%pre_calc(3)))-1), (/m, this%pre_calc(3)/))
      end if

      !from diagonal of i-1, k-1
      if (x-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, 0, 0:d2+1, 0) = reshape(this%recv_buf(m* (this%pre_calc(13)+3*(this%pre_calc(3))) :&
                     m * (this%pre_calc(13)+4*(this%pre_calc(3)))-1), (/m, this%pre_calc(3)/))
      end if

      !from diagonal of j+1, k+1
      if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0:d1+1, d2+1, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*(this%pre_calc(3))) :&
                     m * (this%pre_calc(13)+4*(this%pre_calc(3))+(this%pre_calc(1)))-1), (/m, this%pre_calc(1)/))
      end if

      !from diagonal of j+1, k-1
      if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(1:m, 0:d1+1, d2+1, 0) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*(this%pre_calc(3))+(this%pre_calc(1))) :&
                     m* (this%pre_calc(13)+4*(this%pre_calc(3))+2*(this%pre_calc(1)))-1), (/m, this%pre_calc(1)/))
      end if

      !from diagonal of j-1, k+1
      if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0:d1+1, 0, d3+1) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*(this%pre_calc(3))+2*(this%pre_calc(1))) :&
                     m * (this%pre_calc(13)+4*(this%pre_calc(3))+3*(this%pre_calc(1)))-1), (/m, this%pre_calc(1)/))
      end if

      !from diagonal of j-1, k-1
      if (y-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, 0:d1+1, 0, 0) = reshape(this%recv_buf(m * (this%pre_calc(13)+4*(this%pre_calc(3))+3*(this%pre_calc(1))) :&
                     m * (this%pre_calc(14))-1), (/m, this%pre_calc(1)/))
      end if

      !from diagonal of i+1, j+1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, d1+1, d2+1, d3+1) = this%recv_buf(m * (this%pre_calc(14)) :&
                     m * (this%pre_calc(14)+1)-1)
      end if

      !from diagonal of i+1, j+1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(1:m, d1+1, d2+1, 0) = this%recv_buf(m * (this%pre_calc(14)+1) :&
                     m * (this%pre_calc(14)+2)-1)
      end if

      !from diagonal of i+1, j-1, k+1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, d1+1, 0, d3+1) = this%recv_buf(m * (this%pre_calc(14)+2) :&
                     m * (this%pre_calc(14)+3)-1)
      end if

      !from diagonal of i+1, j-1, k-1
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, d1+1, 0, 0) = this%recv_buf(m * (this%pre_calc(14)+3) :&
                     m * (this%pre_calc(14)+4)-1)
      end if

      !from diagonal of i-1, j+1, k+1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0, d2+1, d3+1) = this%recv_buf(m * (this%pre_calc(14)+4) :&
                     m * (this%pre_calc(14)+5)-1)
      end if

      !from diagonal of i-1, j+1, k-1
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(1:m, 0, d2+1, 0) = this%recv_buf(m * (this%pre_calc(14)+5) :&
                     m * (this%pre_calc(14)+6)-1)
      end if

      !from diagonal of i-1, j-1, k+1
      if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(1:m, 0, 0, d3+1) = this%recv_buf(m * (this%pre_calc(14)+6) :&
                     m * (this%pre_calc(14)+7)-1)
      end if

      !from diagonal of i-1, j-1, k-1
      if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%values(1:m, 0, 0, 0) = this%recv_buf(m * (this%pre_calc(14)+7) :&
                     m * (this%pre_calc(14)+8)-1)
      end if
   end subroutine Get_recv_buf