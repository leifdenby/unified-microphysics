module microphysics_initialisation
   use microphysics_register, only: register_variable, n_compressible_species, n_incompressible_species
   use microphysics_register, only: reset_register => reset

   implicit none

   public init

contains
   subroutine init(microphysics_implementation__in)
      character(len=*), optional, intent(in) :: microphysics_implementation__in
      character(len=32) :: configuration = ''

      namelist /microphysics_config/ configuration

      if (present(microphysics_implementation__in)) then
         configuration = microphysics_implementation__in
      else

         open(1, file='{microphysics_configuration_namelist_filename}',form='formatted',status='old')
         read(1, nml=microphysics_config)
         close(1)
      endif

      call reset_register()

      call register_variable('water_vapour')
      call register_variable('temperature')
      call register_variable('pressure')

      if (configuration == '') then
         print *, "Please choose a microphysics model to use. The available models are:"
         print *, "  {mphys_names_joined}"
         !{mphys_methods_switch_statements}
      else
         print *, "Error the chosen microphysics implementation wasn't found"
         call exit(-1)
      endif

      print *, "Microphysics init complete, tracers:"
      print *, "  incompressible:", n_incompressible_species
      print *, "  compressible:", n_compressible_species

   end subroutine init

   !{microphysics_module_blocks}

end module microphysics_initialisation
