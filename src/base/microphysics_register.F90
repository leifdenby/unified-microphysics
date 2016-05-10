module microphysics_register
   use microphysics_constants, only: kreal
   implicit none
   public n_variables, n_compressible_species, n_incompressible_species

   ! will be equal to the number of hydrometeors requested by the specific implementation, and will be set on init
   integer :: n_variables = 0, n_compressible_species = 0, n_incompressible_species = 0

   integer :: model_constraint = 0
   integer, parameter :: MODEL_CONSTRAINT_ISOBARIC = 1
   integer, parameter :: MODEL_CONSTRAINT_ISOMETRIC = 2

   ! for storing heat capacity of all species
   real(kreal), dimension(:), allocatable :: cp_species
   real(kreal), dimension(:), allocatable :: cv_species
   integer, dimension(:), allocatable :: moments
   ! array which will store which variables of the state vector are species
   ! concentrations
   integer, dimension(:), allocatable :: q_species_flag

   integer :: n_moments__max = -1

   ! variables for indexing into state array
   integer :: idx_cwater=0, idx_rain=0, idx_cice=0, idx_graupel=0, idx_water_vapour=0
   integer :: idx_temp=0, idx_pressure=0

   ! pointers to subroutines used in microphysics calculation
   !procedure(), pointer :: q_flux_function => null()

contains
   subroutine reset()
      n_variables = 0
      n_compressible_species = 0
      n_incompressible_species = 0
      n_moments__max = 0
      idx_cwater = 0
      idx_rain = 0
      idx_cice = 0
      idx_graupel = 0
      idx_water_vapour = 0
      idx_temp = 0
      idx_pressure = 0

      if (allocated(cp_species)) then
         deallocate(cp_species)
         allocate(cp_species(0))
      endif
      if (allocated(cv_species)) then
         deallocate(cv_species)
         allocate(cv_species(0))
      endif
      if (allocated(q_species_flag)) then
         deallocate(q_species_flag)
         allocate(q_species_flag(0))
      endif
   end subroutine reset

   subroutine register_variable(var_name, n_moments)
      use microphysics_constants, only: cp_d, cv_d, cp_v, cv_v
      use microphysics_constants, only: cp_l, cv_l, cp_i, cv_i

      character(len=*), intent(in) :: var_name
      integer, intent(in), optional :: n_moments

      real(kreal) :: cp = 0.0, cv = 0.0
      integer :: species_idx = -1
      integer :: is_spec_conc

      is_spec_conc = 1
      n_variables = n_variables+1
      species_idx = n_variables


      if (trim(var_name) == 'water_vapour') then
         n_compressible_species = n_compressible_species+1
         idx_water_vapour = species_idx
         cp = cp_v
         cv = cv_v
      else if (trim(var_name) == 'cloud_water') then
         idx_cwater = species_idx
         cp = cp_l
         cv = cv_l
         n_incompressible_species = n_incompressible_species+1
      else if (trim(var_name) == 'rain') then
         idx_rain = species_idx
         cp = cp_l
         cv = cv_l
         n_incompressible_species = n_incompressible_species+1
      else if (trim(var_name) == 'cloud_ice') then
         idx_cice = species_idx
         cp = cp_i
         cv = cv_i
         n_incompressible_species = n_incompressible_species+1
      else if (trim(var_name) == 'graupel') then
         idx_graupel = species_idx
         cp = cp_i
         cv = cv_i
         n_incompressible_species = n_incompressible_species+1
      else if (trim(var_name) == 'temperature') then
         idx_temp = species_idx
         cp = 0.0
         cv = 0.0
         is_spec_conc = 0
      else if (trim(var_name) == 'pressure') then
         idx_pressure = species_idx
         cp = 0.0
         cv = 0.0
         is_spec_conc = 0
      else
         print *, "Registration methods not implemented for species ", trim(var_name)
         call exit(-1)
      endif

      call array_append_real(cp_species, cp)
      call array_append_real(cv_species, cv)
      call array_append_int(q_species_flag, is_spec_conc)

      if (present(n_moments)) then
         if (n_moments > n_moments__max) then
            n_moments__max = n_moments
         endif
         call array_append_int(moments, n_moments)
      endif

   end subroutine register_variable

   subroutine array_append_real(arr, item)
      real(kreal), dimension(:), allocatable, intent(inout) :: arr
      real(kreal), intent(in) :: item

      real(kreal), dimension(:), allocatable :: temp_arr
      integer :: n

      if (allocated(arr)) then
         n = size(arr)
         allocate(temp_arr(n))
         temp_arr = arr

         deallocate(arr)
         allocate(arr(n+1))
         arr(1:n) = temp_arr
      else
         n = 0
         allocate(arr(n+1))
      endif

      arr(n+1) = item
   end subroutine array_append_real

   subroutine array_append_int(arr, item)
      integer, dimension(:), allocatable, intent(inout) :: arr
      integer, intent(in) :: item

      integer, dimension(:), allocatable :: temp_arr
      integer :: n

      if (allocated(arr)) then
         n = size(arr)
         allocate(temp_arr(n))
         temp_arr = arr

         deallocate(arr)
         allocate(arr(n+1))
         arr(1:n) = temp_arr
      else
         n = 0
         allocate(arr(n+1))
      endif

      arr(n+1) = item
   end subroutine array_append_int

   !function hydrometeor_is_supported(var_name)
      !logical :: hydrometeor_is_supported
      !character(16), intent(in) :: var_name

      !integer :: j = -1

      !hydrometeor_is_supported = .FALSE.
      !do j=1,size(valid_hydrometeor_names)
          !if (trim(var_name) == trim(valid_hydrometeor_names(j))) then
              !hydrometeor_is_supported = .TRUE.
              !exit
          !endif
      !enddo
   !end function hydrometeor_is_supported

end module microphysics_register
