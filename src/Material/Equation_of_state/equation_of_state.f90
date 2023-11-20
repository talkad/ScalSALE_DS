
module equation_of_state_module
   use data_module, only: data_t
   use data_struct_base, only : data_struct_t

   implicit none
   private
   public :: equation_of_state_t, eos_wrapper_t

   type, abstract :: equation_of_state_t
   contains

      procedure(Calculate_interface), deferred :: Calculate
   end type

   abstract interface
      subroutine Calculate_interface (this, pressure, sound_vel, density, sie,&
                                      dp_de, dp_drho, dt_de, dt_drho, gamma_gas, atomic_mass, temperature, mat_num, nrg_or_tmp,&
                                      nx, ny, nz, vof, emf)
         import :: equation_of_state_t, data_t, data_struct_t
         class (equation_of_state_t)        , intent (inout) :: this        
         class(data_struct_t), pointer, intent (inout) :: pressure
         class(data_struct_t), pointer,  intent (inout) :: sound_vel
         class(data_struct_t), pointer,  intent (inout) :: density
         class(data_struct_t), pointer,  intent (inout) :: sie
         class(data_struct_t), pointer,  intent (inout) :: temperature
         class(data_struct_t), pointer,  intent (inout) :: dp_de
         class(data_struct_t), pointer,  intent (inout) :: dp_drho
         class(data_struct_t), pointer,  intent (inout) :: dt_de
         class(data_struct_t), pointer,  intent (inout) :: dt_drho
         class(data_struct_t), pointer,  intent (in) :: vof

         real(8)                            , intent (in)    :: atomic_mass 
         real(8)                            , intent (in)    :: gamma_gas   
         integer                            , intent (in   ) :: mat_num     

         integer                            , intent (in   ) :: nrg_or_tmp
         integer                            , intent (in   ) :: nx
         integer                            , intent (in   ) :: ny
         integer                            , intent (in   ) :: nz
         real(8)                            , intent (in   ) :: emf



      end subroutine Calculate_interface
   end interface

   type eos_wrapper_t
      class(equation_of_state_t), pointer :: eos
   end type eos_wrapper_t

contains

end module equation_of_state_module
