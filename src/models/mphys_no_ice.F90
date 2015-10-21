!> Simple microphysics implementation which only supports the creation of cloud
!> water from water vapour, no rain or ice-phases are included.
module mphys_no_ice
   use microphysics_register, only: register_incompressible_species, n_species, n_gases, n_solids
   use microphysics_constants, only: kreal

   implicit none
   public init

   ! radwatc: Cloud particle radius (radius - water - cloud)
   real(kreal), parameter :: radwatc = 1.0


   contains
   subroutine init()
      ! register_species(name, n_moments)
      print *, "no ice init"
      call register_incompressible_species('cloud_water', 1)
   end subroutine

   subroutine calc_dqdt(q_g, q_tr, temp, pressure, dq_g, dq_tr)
      use microphysics_register, only: n_moments__max, idx_cwater, idx_dry_air, idx_water_vapour
      use microphysics_constants, only: L_cond

      real(kreal), intent(in) :: temp, pressure
      real(kreal), dimension(n_gases), intent(in) :: q_g
      real(kreal), dimension(n_solids,n_moments__max), intent(in) :: q_tr
      real(kreal), dimension(n_gases), intent(out) :: dq_g
      real(kreal), dimension(n_solids,n_moments__max), intent(out) :: dq_tr

      real(kreal) :: dq_cond_evap = 0.0
      real(kreal) :: dq_cloudwater_formation = 0.0
      real(kreal) :: temp_local = 0.0

      dq_cond_evap = dqdt_cond_evap_cloudwater(temp, pressure, q_g(idx_water_vapour), q_tr(idx_cwater,1))

      temp_local = temp + L_cond*dq_cond_evap

      dq_cloudwater_formation = calc_dq_cloudwater_formation(temp_local, pressure, q_g(idx_water_vapour) - dq_cond_evap)

      dq_g(idx_water_vapour) = -dq_cond_evap -dq_cloudwater_formation
      dq_tr(idx_cwater,1) = dq_cond_evap + dq_cloudwater_formation

      contains
         !> Create cloudwater by consuming excess super saturation
         function calc_dq_cloudwater_formation(temp, pressure, q_v) result(dq_cwater)
            use microphysics_common, only: saturation_vapour_pressure
            use microphysics_constants, only: R_v, R_d

            real(kreal), intent(in) :: temp, pressure, q_v
            real(kreal) :: dq_cwater
            real(kreal) :: satpw = 0.0, satww = 0.0

            satpw=saturation_vapour_pressure(temp)

            satww = R_d/R_v*satpw/pressure

            dq_cwater=min(q_v-satww, 0.0)
            
         end function calc_dq_cloudwater_formation

         function dqdt_cond_evap_cloudwater(temp, pressure, q_v, q_cw)
            use microphysics_common, only: thermal_conductivity
            use microphysics_common, only: water_vapour_diffusivity, saturation_vapour_pressure
            use microphysics_constants, only: R_v, L_cond, rho_w, R_d

            real(kreal) :: dqdt_cond_evap_cloudwater
            real(kreal), intent(in) :: temp, pressure, q_cw, q_v

            real(kreal) :: diffk, diffd, satpw
            real(kreal) :: pfak
            real(kreal) :: satww
            real(kreal) :: supsatw
            real(kreal) :: evap
            real(kreal) :: rgas
            real(kreal) :: xnwatc  ! number of cloud droplets per unit gas volume
            real(kreal) :: cp_mixture, cv_mixture

            ! TODO: Move, this is a constant
            real(kreal), parameter :: r1 = 1.0, r4 = 4.0
            real(kreal), parameter :: xpi = 4.0/3.0*3.14
            real(kreal), parameter :: pi = 3.14

            diffk=thermal_conductivity(temp)
            diffd=water_vapour_diffusivity(temp, pressure)

            satpw=saturation_vapour_pressure(temp)

            !pfak=rgasnew/(R_v*pressure)*gasnew
            !satww=satpw*pfak

            ! qv_sat ~ Rd/Rv * pv_sat/p
            satww = R_d/R_v*satpw/pressure

            supsatw=q_v/satww-r1

            evap=supsatw/((L_cond/(R_v*temp)-r1)*L_cond/(diffk*temp) + R_v*temp/(diffd*satpw))

            xnwatc=q_cw/(rho_w*xpi*radwatc**3)

            dqdt_cond_evap_cloudwater=r4*pi*xnwatc*radwatc*evap

         end function dqdt_cond_evap_cloudwater

   end subroutine
end module mphys_no_ice
