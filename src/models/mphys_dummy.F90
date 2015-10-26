!> Dummy microphysics implementation which displays which functions must be
!> implemented

module mphys_dummy
   use microphysics_register, only: register_incompressible_species, n_species, n_gases, n_solids
   use microphysics_constants, only: kreal

   implicit none
   public init

   real(kreal), parameter :: radwatc = 1.0


   contains
   subroutine init()
      call register_incompressible_species('cloud_water', 1)
      call register_incompressible_species('rain', 1)
      call register_incompressible_species('cloud_ice', 1)
      call register_incompressible_species('graupel', 1)
   end subroutine

   subroutine calc_dqdt(q_g, q_tr, temp, pressure, dq_g, dq_tr)
      use microphysics_register, only: n_moments__max, idx_cwater, idx_water_vapour
      use microphysics_constants, only: L_cond

      real(kreal), intent(in) :: temp, pressure
      real(kreal), dimension(n_gases), intent(in) :: q_g
      real(kreal), dimension(n_solids,n_moments__max), intent(in) :: q_tr
      real(kreal), dimension(n_gases), intent(out) :: dq_g
      real(kreal), dimension(n_solids,n_moments__max), intent(out) :: dq_tr

   end subroutine
end module mphys_dummy
