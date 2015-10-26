!> Simple microphysics implementation which only supports the creation of cloud
!> water from water vapour, no rain or ice-phases are included.

module mphys_no_ice
   use microphysics_register, only: register_incompressible_species, n_species, n_gases, n_solids
   use microphysics_constants, only: kreal, nan

   implicit none
   public init, rho_f, dqr_dt__autoconversion

   contains
   subroutine init()
      ! register_species(name, n_moments)
      print *, "no ice init"
      call register_incompressible_species('cloud_water', 1)
      call register_incompressible_species('rain', 1)
   end subroutine

   subroutine calc_dqdt(q_g, q_tr, temp, pressure, dqdt_g, dqdt_tr, dTdt)
      use microphysics_register, only: n_moments__max, idx_cwater, idx_dry_air, idx_water_vapour, idx_rain
      use microphysics_constants, only: L_v => L_cond
      use microphysics_common, only: cp_mixture

      real(kreal), intent(in) :: temp, pressure
      real(kreal), intent(out) :: dTdt
      real(kreal), dimension(n_gases), intent(in) :: q_g
      real(kreal), dimension(n_solids,n_moments__max), intent(in) :: q_tr
      real(kreal), dimension(n_gases), intent(out) :: dqdt_g
      real(kreal), dimension(n_solids,n_moments__max), intent(out) :: dqdt_tr

      real(kreal) :: ql = nan, qv = nan, qg = nan, qr = nan, qd = nan
      real(kreal) :: rho = nan, rho_g = nan
      real(kreal) :: dqrdt_autoconv = nan, dqrdt_accre = nan, dqldt_condevap = nan
      real(kreal) :: cp_m = nan

      cp_m = cp_mixture(q_g, q_tr)

      ! pick out specific concentrations from state vectors
      ql = q_tr(idx_cwater, 1)
      qr = q_tr(idx_rain, 1)
      qv = q_g(idx_water_vapour)
      qd = q_g(idx_dry_air)
      qg = qv + qd

      ! compute gas and mixture density using equation of state
      rho = rho_f(qd, qv, ql, qr, pressure, temp)
      rho_g = rho_f(qd, qv, 0.0_kreal, 0.0_kreal, pressure, temp)

      ! compute time derivatives for each process
      dqrdt_autoconv = dqr_dt__autoconversion(ql, qg, rho_g)
      dqrdt_accre    = dqr_dt__accretion(ql, qg, rho_g, qr)
      dqldt_condevap = dql_dt__condensation_evaporation(rho, rho_g, qv, ql, temp, pressure)

      ! combine to create time derivatives for species
      dqdt_g(idx_water_vapour) = -dqldt_condevap
      dqdt_tr(idx_cwater, 1)   =  dqldt_condevap - dqrdt_autoconv - dqrdt_accre
      dqdt_tr(idx_rain, 1)     =                   dqrdt_autoconv + dqrdt_accre

      dTdt = L_v/cp_m*dqldt_condevap

   end subroutine


   function rho_f(qd, qv, ql, qr, p, temp) result(rho)
      use microphysics_constants, only: R_v, R_d, rho_l => rho_w

      real(kreal), intent(in) :: qd, qv, ql, qr, p, temp
      real(kreal) :: rho, rho_inv = nan

      rho_inv = (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l

      rho = 1.0/rho_inv
   end function rho_f

   function dqr_dt__autoconversion(ql, qg, rho_g)
      real(kreal), intent(in) :: ql, qg, rho_g
      real(kreal) :: dqr_dt__autoconversion

      real(kreal) :: k_c = 1.0e-3, a_c = 5.0e-4

      ! TODO: what happens if qg < 0.0 ?
      dqr_dt__autoconversion = k_c*(ql - qg/rho_g*a_c)
      dqr_dt__autoconversion = max(0.0, dqr_dt__autoconversion)
   end function dqr_dt__autoconversion

   function dqr_dt__accretion(ql, qg, rho_g, qr)
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: ql, qg, rho_g, qr
      real(kreal) :: dqr_dt__accretion

      real(kreal), parameter :: G3p5 = 3.32399614155  ! = Gamma(3.5)
      real(kreal), parameter :: N0r = 1.e7  ! [m^-4]
      real(kreal), parameter :: a_r = 201.0  ! [m^.5 s^-1]
      real(kreal), parameter :: rho0 = 1.12

      real(kreal) :: lambda_r = nan

      lambda_r = (pi*(qg*rho_l)/(qr*rho_g)*N0r)**(1./4.)

      dqr_dt__accretion = pi/4.*N0r*a_r*sqrt(rho0/rho_g)*G3p5*lambda_r**(-3.5)*ql

      dqr_dt__accretion = max(0.0, dqr_dt__accretion)
   end function dqr_dt__accretion

   function dql_dt__condensation_evaporation(rho, rho_g, qv, ql, T, p)
      use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
      use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
      use microphysics_common, only: Ka_f => thermal_conductivity
      use microphysics_common, only: Dv_f => water_vapour_diffusivity
      use microphysics_constants, only: Lv => L_cond, R_v
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: rho, rho_g, qv, ql, T, p
      real(kreal) :: dql_dt__condensation_evaporation

      real(kreal), parameter :: r0 = 0.1e-6_kreal  ! initial cloud droplet radius
      real(kreal), parameter :: r_crit = 5.0e-6_kreal ! critical radius after which cloud-droplet number is increased
      real(kreal), parameter :: N0 = 200*1.0e6_kreal  ! initial cloud droplet number

      real(kreal) :: r_c = nan, Nc = nan
      real(kreal) :: pv_sat = nan, qv_sat = nan, Sw = nan
      real(kreal) :: Dv = nan, Fd = nan
      real(kreal) :: Ka = nan, Fk = nan

      real(kreal) :: r4_3 = 4.0_kreal/3.0_kreal
      real(kreal) :: r1_3 = 1.0_kreal/3.0_kreal

      ! calculate saturation concentration of water vapour
      pv_sat = pv_sat_f(T)
      qv_sat = qv_sat_f(T, p)

      if (ql == 0.0_kreal) then
         r_c = r0
      else
         r_c = (ql*rho/(r4_3*pi*N0*rho_l))**r1_3
      endif

      if (r_c > r_crit) then
         ! if cloud droplet radius with initial number of droplets is larger
         ! than a critial size assume that instead more droplets are made,
         ! all with the critical radius
         Nc = ql*rho/(r4_3*pi*rho_l*r_crit**3.0_kreal)
         r_c = r_crit
      else
         Nc = N0
      endif

      ! condensation evaporation of cloud droplets (given number of droplets
      ! and droplet radius calculated above)
      Sw = qv/qv_sat

      Ka = Ka_f(T)
      Fk = (Lv/(R_v*T) - 1._kreal) * (Lv/(Ka*T))

      Dv = Dv_f(T, p)
      Fd = R_v*T/(pv_sat*Dv)

      ! compute rate of change of condensate from diffusion
      dql_dt__condensation_evaporation = 4.*pi*1./rho*Nc*r_c*(Sw - 1.0)/(Fk + Fd)
   end function dql_dt__condensation_evaporation
end module mphys_no_ice
