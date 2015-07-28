module microphysics_register
   implicit none
   public n_species

   ! will be equal to the number of hydrometeors requested by the specific implementation, and will be set on init
   integer :: n_species = -1

   ! for storing heat capacity of all tracers (incompressible and compressible)
   real, dimension(:), allocatable :: cp_all
   real, dimension(:), allocatable :: cv_all
   integer, dimension(:), allocatable :: moments
   integer, dimension(:), allocatable :: compressible_species, incompressible_species

   integer :: n_moments__max = -1

   ! variables for indexing into state array
   integer :: idx_cwater, idx_rain, idx_cice, idx_graupel, idx_dry_air, idx_water_vapour

   ! pointers to subroutines used in microphysics calculation
   procedure(), pointer :: q_flux_function => null()

contains
   subroutine register_species(var_name, n_moments__in, cp__in, cv__in)
      use microphysics_constants, only: cp_l, cv_l, cp_i, cv_i, cp_d, cv_d, cp_v, cv_v
      character(len=*), intent(in) :: var_name
      integer, intent(in), optional :: n_moments__in
      real, intent(in), optional :: cp__in, cv__in

      real :: cp, cv
      integer :: species_idx = -1
      logical :: is_compressible = .false.
      integer :: n_moments = -1

      species_idx = n_species+1

      if (trim(var_name) == 'cloud_water') then
         idx_cwater = species_idx
         cp = cp_l
         cv = cv_l
         is_compressible = .false.
      else if (trim(var_name) == 'rain') then
         idx_rain = species_idx
         cp = cp_l
         cv = cv_l
         is_compressible = .false.
      else if (trim(var_name) == 'cloud_ice') then
         idx_cice = species_idx
         cp = cp_i
         cv = cv_i
         is_compressible = .false.
      else if (trim(var_name) == 'graupel') then
         idx_graupel = species_idx
         cp = cp_i
         cv = cv_i
         is_compressible = .false.
      else if (trim(var_name) == 'dry_air') then
         idx_dry_air = species_idx
         cp = cp_d
         cv = cv_d
         is_compressible = .true.
      else if (trim(var_name) == 'water_vapour') then
         idx_water_vapour = species_idx
         cp = cp_v
         cv = cv_v
         is_compressible = .true.
      else
         print *, "Registration methods not implemented for hydrometeor ", trim(var_name)
      endif

      if (present(n_moments__in)) then
         if (is_compressible) then
            print *, "Error: The variable name provided matches a compressible tracer, these should not have more than one moment"
         else
            n_moments = n_moments__in
         endif
      else
         if (is_compressible) then
            n_moments = 1
         else
            print *, "Error: The variable name provided matches an incompressible tracer, please provided the number of moments"
         endif
      endif

      if (present(cp__in)) then
         cp = cp__in
      endif
      if (present(cv__in)) then
         cv = cv__in
      endif

      if (n_moments > n_moments__max) then
         n_moments__max = n_moments
      endif

      call array_append_real(cp_all, cp)
      call array_append_real(cv_all, cv)
      call array_append_int(moments, n_moments)

      if (is_compressible) then
         call array_append_int(compressible_species, 1)
         call array_append_int(incompressible_species, 0)
      else
         call array_append_int(compressible_species, 0)
         call array_append_int(incompressible_species, 1)
      endif

   end subroutine register_species

   subroutine array_append_real(arr, item)
      real, dimension(:), allocatable, intent(inout) :: arr
      real, intent(in) :: item

      real, dimension(:), allocatable :: temp_arr
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
