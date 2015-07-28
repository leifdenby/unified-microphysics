!> Simple microphysics implementation which only supports the creation of cloud
!> water from water vapour, no rain or ice-phases are included.
module mphys_no_ice
   use microphysics_register, only: register_species, n_species

   implicit none
   public init

   ! radwatc: Cloud particle radius (radius - water - cloud)
   real, parameter :: radwatc = 1.0


   contains
   subroutine init()
      ! register_species(name, n_moments)
      print *, "no ice init"
      call register_species('cloud_water', 1)
   end subroutine

   subroutine calc_dq(q, dt, temp, pressure, dq)
      use microphysics_register, only: n_moments__max, idx_cwater, idx_dry_air, idx_water_vapour

      real, intent(in) :: dt, temp, pressure
      real, dimension(n_species,n_moments__max), intent(in) :: q
      real, dimension(n_species,n_moments__max), intent(out) :: dq

      real :: dq_cond_evap = 0.0

      dq_cond_evap = dqdt_cond_evap_cloudwater(temp, pressure, q(idx_water_vapour,1), q(idx_cwater,1))*dt

      dq(idx_water_vapour,1) = -dq_cond_evap
      dq(idx_cwater,1) = dq_cond_evap

      contains

         function dqdt_cond_evap_cloudwater(temp, pressure, q_v, q_cw)
            use microphysics_common, only: thermal_conductivity
            use microphysics_common, only: water_vapour_diffusivity, saturation_vapour_pressure
            use microphysics_constants, only: R_v, L_cond, rho_w, R_d

            real :: dqdt_cond_evap_cloudwater
            real, intent(in) :: temp, pressure, q_cw, q_v

            real :: diffk, diffd, satpw
            real :: pfak
            real :: satww
            real :: supsatw
            real :: evap
            real :: rgas
            real :: xnwatc  ! number of cloud droplets per unit gas volume
            real :: cp_mixture, cv_mixture

            ! TODO: Move, this is a constant
            real, parameter :: r1 = 1.0, r4 = 4.0
            real, parameter :: xpi = 4.0/3.0*3.14
            real, parameter :: pi = 3.14

            diffk=thermal_conductivity(temp)
            diffd=water_vapour_diffusivity(temp, pressure)

            satpw=saturation_vapour_pressure(temp)

            !pfak=rgasnew/(R_v*pressure)*gasnew
            !satww=satpw*pfak

            ! qv_sat ~ Rd/Rv * pv_sat/p
            satww = satpw * R_d/R_v*satpw/pressure

            supsatw=q_v/satww-r1

            evap=supsatw/((L_cond/(R_v*temp)-r1)*L_cond/(diffk*temp) + R_v*temp/(diffd*satpw))

            xnwatc=q_cw/(rho_w*xpi*radwatc**3)

            dqdt_cond_evap_cloudwater=r4*pi*xnwatc*radwatc*evap

         end function dqdt_cond_evap_cloudwater

   end subroutine
end module mphys_no_ice
