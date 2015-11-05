!> Dummy microphysics implementation which displays which functions must be
!> implemented

module mphys_dummy
   use microphysics_register, only: register_variable
   use microphysics_constants, only: kreal

   implicit none
   public init

   real(kreal), parameter :: radwatc = 1.0

   contains
   subroutine init()
      call register_variable('cloud_water', 1)
      call register_variable('rain', 1)
      call register_variable('cloud_ice', 1)
      call register_variable('graupel', 1)
   end subroutine

   pure function dydt(t, y)
      use microphysics_register, only: idx_cwater, idx_water_vapour
      use microphysics_constants, only: L_cond
      integer, parameter :: n_species = 4

      real(kreal), intent(in) :: t
      real(kreal), dimension(n_species), intent(in) :: y
      real(kreal), dimension(n_species) :: dydt

   end function
end module mphys_dummy
