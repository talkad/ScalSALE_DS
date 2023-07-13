
module ideal_gas_module
   use equation_of_state_module, only : equation_of_state_t
   use data_module, only : data_t
   use constants_module, only : AVOGADRO, K_BOLTZMAN
   use indexer_module
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
      real(8), dimension (:), pointer, intent (inout) :: pressure
      real(8), dimension (:), pointer, intent (inout) :: sound_vel
      real(8), dimension (:), pointer, intent (inout) :: density
      real(8), dimension (:), pointer, intent (inout) :: sie
      real(8), dimension (:), pointer, intent (inout) :: temperature
      real(8), dimension (:), pointer, intent (inout) :: dp_de
      real(8), dimension (:), pointer, intent (inout) :: dp_drho
      real(8), dimension (:), pointer, intent (inout) :: dt_de
      real(8), dimension (:), pointer, intent (inout) :: dt_drho
      real(8)                            , intent (in)    :: atomic_mass 
      real(8)                            , intent (in)    :: gamma_gas   
      integer                            , intent (in)    :: mat_num     

      integer                            , intent (in) :: nrg_or_tmp
      integer                            , intent (in) :: nx  
      integer                            , intent (in) :: ny  
      integer                            , intent (in) :: nz  
      real(8), dimension (:), pointer, intent (in) :: vof
      real(8)                            , intent (in) :: emf 

      real(8) :: gamma1 
      real(8) :: atomic_weight 
      integer :: i, j, k
      type(indexer_t), pointer ::  index_mapper
      integer, dimension(:,:,:,:), pointer   ::   mapper
      integer :: csr_idx

      gamma1 = gamma_gas - 1d0
      atomic_weight = atomic_mass / AVOGADRO

      index_mapper => get_instance()
      mapper => index_mapper%mapper

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx

               csr_idx = mapper(mat_num,i,j,k)
               if (csr_idx == -1) cycle

               if (vof(csr_idx) >= emf) then
                  if (nrg_or_tmp == 1) then
                     sie(csr_idx) = 1d0 / gamma1 * K_BOLTZMAN * temperature(csr_idx) / atomic_weight
                  else
                     temperature(csr_idx) = atomic_weight * sie(csr_idx) * gamma1 / K_BOLTZMAN
                  end if
                  pressure (csr_idx) = gamma1 * sie(csr_idx) * density(csr_idx) + 1d-25
                  sound_vel(csr_idx) = gamma1 * gamma_gas * sie(csr_idx)
                  dp_de    (csr_idx) = gamma1 * density(csr_idx)
                  dp_drho  (csr_idx) = gamma1 * sie(csr_idx)
                  dt_de    (csr_idx) = atomic_weight * gamma1 / K_BOLTZMAN
                  dt_drho  (csr_idx) = 0d0
                  if (sound_vel(csr_idx) <= 0d0) then
                     sound_vel(csr_idx) = 1d-10
                     dp_drho  (csr_idx) = 1d-20
                     pressure (csr_idx) = 1d-25
                  end if
               end if
            end do
         end do
      end do

      end subroutine Calculate_ideal_gas

end module ideal_gas_module
