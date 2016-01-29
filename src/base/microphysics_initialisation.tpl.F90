module microphysics_initialisation
   use microphysics_register, only: register_variable, n_compressible_species, n_incompressible_species
   use microphysics_register, only: reset_register => reset
   use microphysics_register, only: model_constraint, MODEL_CONSTRAINT_ISOBARIC, MODEL_CONSTRAINT_ISOMETRIC

   use microphysics_integration, only: integrate_with_constraint
   use isobaric_integration_helpers, only: integrate_isobaric
   use isometric_integration_helpers, only: integrate_isometric

   implicit none

   public init

contains
   subroutine init(microphysics_implementation, system_constraint)
      character(len=*), intent(in) :: microphysics_implementation, system_constraint
      character(len=32) :: configuration = ''
      character(len=100) :: msg

      call init_with_message(microphysics_implementation, system_constraint, msg)
   end subroutine

   subroutine init_with_message(microphysics_implementation, system_constraint, msg)
      character(len=*), intent(in) :: microphysics_implementation, system_constraint
      character(len=32) :: configuration = ''
      character(len=*), intent(out) :: msg

      namelist /microphysics_config/ configuration

      configuration = microphysics_implementation

      call reset_register()

      call register_variable('water_vapour')
      call register_variable('temperature')
      call register_variable('pressure')

      print *, "Microphysics init (", trim(configuration), ")"

      if (configuration == '') then
         print *, "Please choose a microphysics model to use. The available models are:"
         print *, "  {mphys_names_joined}"
         !{mphys_methods_switch_statements}
      else
         msg = "Error the chosen microphysics implementation wasn't found"
      endif

      if (trim(system_constraint) == 'isobaric') then
         print *, "integrating with isobaric constraint"
         integrate_with_constraint => integrate_isobaric
         model_constraint = MODEL_CONSTRAINT_ISOBARIC
      else if (trim(system_constraint) == 'isometric') then
         print *, "integrating with isometric constraint"
         integrate_with_constraint => integrate_isometric
         model_constraint = MODEL_CONSTRAINT_ISOMETRIC
      else
         msg = "Error `system_constraint` must be either `isobaric` or `isometric`"
      endif

      print *, "Microphysics init complete, tracers:"
      print *, "  incompressible:", n_incompressible_species
      print *, "  compressible:", n_compressible_species

   end subroutine init_with_message

   !{microphysics_module_blocks}

end module microphysics_initialisation
