
module ideal_gas_module
   use equation_of_state_module, only : equation_of_state_t
   use data_module, only : data_t
   use constants_module, only : AVOGADRO, K_BOLTZMAN

   use data_struct_base, only : data_struct_t

   implicit none
   private

   type, extends(equation_of_state_t), public :: ideal_gas_t
      private

   contains

      procedure, public :: Calculate => Calculate_ideal_gas

   end type ideal_gas_t



contains





   subroutine Calculate_ideal_gas (this, pressure, sound_vel, density, sie,&
                                   dp_de, dp_drho, dt_de, dt_drho, gamma_gas, atomic_mass, temperature, mat_num, nrg_or_tmp,&
                                   nx, ny, nz, vof, emf)
      implicit none
      class (ideal_gas_t)                , intent (inout) :: this        
      class(data_struct_t), pointer, intent (inout) :: pressure
      class(data_struct_t), pointer, intent (inout) :: sound_vel
      class(data_struct_t), pointer, intent (inout) :: density
      class(data_struct_t), pointer, intent (inout) :: sie
      class(data_struct_t), pointer, intent (inout) :: temperature
      class(data_struct_t), pointer, intent (inout) :: dp_de
      class(data_struct_t), pointer, intent (inout) :: dp_drho
      class(data_struct_t), pointer, intent (inout) :: dt_de
      class(data_struct_t), pointer, intent (inout) :: dt_drho
      class(data_struct_t), pointer, intent (in) :: vof

      real(8)                            , intent (in)    :: atomic_mass 
      real(8)                            , intent (in)    :: gamma_gas   
      integer                            , intent (in)    :: mat_num     

      integer                            , intent (in) :: nrg_or_tmp
      integer                            , intent (in) :: nx  
      integer                            , intent (in) :: ny  
      integer                            , intent (in) :: nz  
      real(8)                            , intent (in) :: emf 

      real(8) :: gamma1 
      real(8) :: atomic_weight 
      integer :: i, j, k

      gamma1 = gamma_gas - 1d0
      atomic_weight = atomic_mass / AVOGADRO

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx

               if (vof%get_item(mat_num,i,j,k) >= emf) then
                  if (nrg_or_tmp == 1) then
                     call sie%add_item(mat_num, i , j, k, 1d0 / gamma1 * K_BOLTZMAN * temperature%get_item(mat_num,i,j,k) / atomic_weight)
                  else
                     call temperature%add_item(mat_num, i , j, k, atomic_weight * sie%get_item(mat_num,i,j,k) * gamma1 / K_BOLTZMAN)
                  end if
                  call pressure%add_item(mat_num,i,j,k,gamma1 * sie%get_item(mat_num,i,j,k) * density%get_item(mat_num,i,j,k) + 1d-25)
                  call sound_vel%add_item(mat_num, i , j, k, gamma1 * gamma_gas * sie%get_item(mat_num,i,j,k))
                  call dp_de%add_item(mat_num, i , j, k, gamma1 * density%get_item(mat_num,i,j,k))
                  call dp_drho %add_item(mat_num, i , j, k, gamma1 * sie%get_item(mat_num,i,j,k))
                  call dt_de%add_item(mat_num, i , j, k, atomic_weight * gamma1 / K_BOLTZMAN)
                  call dt_drho%add_item(mat_num, i , j, k, 0d0)
                  if (sound_vel%get_item(mat_num,i,j,k) <= 0d0) then
                     call sound_vel%add_item(mat_num, i , j, k, 1d-10)
                     call dp_drho%add_item(mat_num, i , j, k, 1d-20)
                     call pressure%add_item(mat_num, i , j, k, 1d-25)
                  end if
               end if
            end do
         end do
      end do

      end subroutine Calculate_ideal_gas

end module ideal_gas_module
