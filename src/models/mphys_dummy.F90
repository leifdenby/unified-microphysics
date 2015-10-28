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

   subroutine calc_dydt(y, n_variables, dydt)
      use microphysics_register, only: idx_cwater, idx_water_vapour
      use microphysics_constants, only: L_cond

      integer, intent(in) :: n_variables
      real(kreal), dimension(n_variables), intent(in) :: y
      real(kreal), dimension(n_variables), intent(out) :: dydt

   end subroutine
end module mphys_dummy
