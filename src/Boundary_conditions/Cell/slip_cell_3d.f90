
module slip_cell_3d_module
    use cell_boundary_condition_module, only : cell_boundary_condition_t
    use data_module    , only : data_t
    use quantity_module, only : quantity_t

    implicit none
    private

    type, extends(cell_boundary_condition_t), public :: slip_cell_3d_t
        private
    contains

        procedure, public :: Calculate => Slip_cell_3d_calculate
    end type slip_cell_3d_t

contains
    subroutine Slip_cell_3d_calculate (this, c_quantity, edge_num)
        class (slip_cell_3d_t) , intent (in out)     :: this      
        class(quantity_t)   , intent (in out)     :: c_quantity  
        integer             , intent (in)         :: edge_num  

        real(8), dimension(:,:,:), pointer :: values
        real(8), dimension(:,:,:,:), pointer :: values_4d
        integer :: i, j, k, nxp, nyp, nzp,m, nmats
        integer :: debug1, debug2, debug3

        debug1 = 0
        debug2 = 0
        debug3 = 0

        nxp = c_quantity%d1 + 1
        nyp = c_quantity%d2 + 1
        nzp = c_quantity%d3 + 1

        if (associated(c_quantity%data_4d)) then
            call c_quantity%Point_to_data(values_4d)
            nmats = c_quantity%data_4d%nmats
            select case(edge_num)
                case(1)
                    ! print*, 'ddddddddddd11111111111'
                    call debug(values_4d, 'material_results/split_debug1.txt', nzp, nyp, nxp, nmats)
                    i = 1

                    do k = 0, nzp
                        do j = 0, nyp
                            do m = 1, nmats
                                values_4d(m, i-1, j, k) = values_4d(m, i, j, k)
                                ! print*, m, i-1, j, k, 'val', values_4d(m, i, j, k)

                            end do
                        end do
                    end do 

                    ! print*, 'ddddddddddd22222222222222222'
                    call debug(values_4d, 'material_results/split_debug2.txt', nzp, nyp, nxp, nmats)

                case(2)
                    i = nxp
                    do k = 0, nzp
                        do j = 0, nyp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i-1, j, k)
                            end do
                        end do
                    end do

                case(3)
                    j = 1
                    do k = 0, nzp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j-1, k) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do

                case(4)
                    j = nyp
                    do k = 0, nzp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i, j-1, k)
                            end do
                        end do
                    end do


                case(5)
                    k = 1
                    do j = 0, nyp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k-1) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do

                case(6)
                    k = nzp
                    do j = 0, nyp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i, j, k-1)
                            end do
                        end do
                    end do

            end select
        else
            call c_quantity%Point_to_data(values)

            select case(edge_num)
                case(1)
                    i = 1
                    do k = 0, nzp
                        do j = 0, nyp
                            values(i-1, j, k) = values(i, j, k)
                        end do
                    end do


                case(2)
                    i = nxp
                    do k = 0, nzp
                        do j = 0, nyp
                            values(i, j, k) = values(i-1, j, k)
                        end do
                    end do

                case(3)
                    j = 1
                    do k = 0, nzp
                        do i = 0, nxp
                            values(i, j-1, k) = values(i, j, k)
                        end do
                    end do

                case(4)
                    j = nyp
                    do k = 0, nzp
                        do i = 0, nxp
                            values(i, j, k) = values(i, j-1, k)
                        end do
                    end do

                case(5)
                    k = 1
                    do j = 0, nyp
                        do i = 0, nxp
                            values(i, j, k-1) = values(i, j, k)
                        end do
                    end do

                case(6)
                    k = nzp
                    do j = 0, nyp
                        do i = 0, nxp
                            values(i, j, k) = values(i, j, k-1)
                        end do
                    end do
            end select
        end if
        return
    end subroutine


    subroutine debug(arr, file_name, nzp, nyp, nxp, nmats)
        real(8), dimension(:,:,:,:), pointer, intent(in)   ::   arr
        integer, intent(in)                              ::   nzp, nyp, nxp, nmats
        character(len=*), intent(in)                      ::   file_name

        integer :: i,j,k,m
        integer :: unit
        integer :: total_debug
        total_debug = 0


        open (unit=414, file=file_name,  status = 'replace')  
        
        
        do k = 0, nzp
            do j = 0, nyp
                do i = 0, nxp
                    do m = 1, nmats

                        if (arr(m,i,j,k) == 0)   total_debug = total_debug + 1
                        ! write(414,*) k,j,i,m, arr(m,i,j,k)
                        write(414,*) arr(m,i,j,k)

                    end do
                end do
            end do
        end do
        
        close (414)


    end subroutine debug

end module slip_cell_3d_module
