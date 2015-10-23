module microphysics_register
   use microphysics_constants, only: kreal
   implicit none
   public n_species

   ! will be equal to the number of hydrometeors requested by the specific implementation, and will be set on init
   integer :: n_species = 0, n_gases = 0, n_solids = 0

   ! for storing heat capacity of all tracers (incompressible and compressible)
   real(kreal), dimension(:), allocatable :: cp_gases
   real(kreal), dimension(:), allocatable :: cv_gases
   real(kreal), dimension(:), allocatable :: cp_solids
   real(kreal), dimension(:), allocatable :: cv_solids
   integer, dimension(:), allocatable :: moments
   !integer, dimension(:), allocatable :: compressible_species, incompressible_species

   integer :: n_moments__max = -1

   ! variables for indexing into state array
   integer :: idx_cwater=0, idx_rain=0, idx_cice=0, idx_graupel=0, idx_dry_air=0, idx_water_vapour=0

   ! pointers to subroutines used in microphysics calculation
   procedure(), pointer :: q_flux_function => null()

contains
   subroutine register_compressible_species(var_name, cp__in, cv__in)
      use microphysics_constants, only: cp_d, cv_d, cp_v, cv_v
      character(len=*), intent(in) :: var_name
      real(kreal), intent(in), optional :: cp__in, cv__in

      real(kreal) :: cp, cv
      integer :: species_idx = -1

      n_gases = n_gases+1
      species_idx = n_gases

      if (trim(var_name) == 'dry_air') then
         idx_dry_air = species_idx
         cp = cp_d
         cv = cv_d
      else if (trim(var_name) == 'water_vapour') then
         idx_water_vapour = species_idx
         cp = cp_v
         cv = cv_v
      else
         print *, "Registration methods not implemented for compressible species ", trim(var_name)
      endif

      if (present(cp__in)) then
         cp = cp__in
      endif
      if (present(cv__in)) then
         cv = cv__in
      endif

      call array_append_real(cp_gases, cp)
      call array_append_real(cv_gases, cv)

      n_species = n_species+1
   end subroutine register_compressible_species

   subroutine register_incompressible_species(var_name, n_moments, cp__in, cv__in)
      use microphysics_constants, only: cp_l, cv_l, cp_i, cv_i
      character(len=*), intent(in) :: var_name
      integer, intent(in) :: n_moments
      real(kreal), intent(in), optional :: cp__in, cv__in

      real(kreal) :: cp, cv
      integer :: species_idx = -1

      n_solids = n_solids+1
      species_idx = n_solids

      if (trim(var_name) == 'cloud_water') then
         idx_cwater = species_idx
         cp = cp_l
         cv = cv_l
      else if (trim(var_name) == 'rain') then
         idx_rain = species_idx
         cp = cp_l
         cv = cv_l
      else if (trim(var_name) == 'cloud_ice') then
         idx_cice = species_idx
         cp = cp_i
         cv = cv_i
      else if (trim(var_name) == 'graupel') then
         idx_graupel = species_idx
         cp = cp_i
         cv = cv_i
      else
         print *, "Registration methods not implemented for hydrometeor ", trim(var_name)
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

      call array_append_real(cp_solids, cp)
      call array_append_real(cv_solids, cv)
      call array_append_int(moments, n_moments)

      n_species = n_species+1
   end subroutine register_incompressible_species

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
