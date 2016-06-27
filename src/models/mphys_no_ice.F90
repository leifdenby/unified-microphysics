!> Simple microphysics implementation which only supports the creation of cloud
!> water from water vapour and rain-droplets through auto-conversion and
!> accretion, no ice-phases are included.

module mphys_no_ice
   use microphysics_register, only: register_variable
   use microphysics_constants, only: kreal
   use microphysics_register, only: n_variables

   implicit none

   public init, rho_f, dqr_dt__autoconversion

   contains
   subroutine init(p_disable_rain)
      logical, optional :: p_disable_rain

      call register_variable('cloud_water', 1)
      call register_variable('rain', 1)
   end subroutine

   pure function dydt(t, y, c_m)
      use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain, idx_temp, idx_pressure
      use microphysics_constants, only: L_v => L_cond

      real(kreal), dimension(:), intent(in) :: y
      real(kreal), intent(in) :: t, c_m
      real(kreal) :: dydt(size(y))

      real(kreal) :: ql, qv, qg, qr, qd
      real(kreal) :: rho, rho_g
      real(kreal) :: dqrdt_autoconv, dqrdt_accre, dqrdt_condevap, dqldt_condevap
      real(kreal) :: temp, pressure

      ! OBS: it's important to make sure that return variable is initiated to
      ! zero
      dydt = 0.0

      temp = y(idx_temp)
      pressure = y(idx_pressure)

      ! pick out specific concentrations from state vectors
      ql = y(idx_cwater)
      qr = y(idx_rain)
      qv = y(idx_water_vapour)
      qd = 1.0_kreal - ql - qr - qv
      qg = qv + qd

      ! compute gas and mixture density using equation of state
      rho = rho_f(qd, qv, ql, qr, pressure, temp)
      rho_g = rho_f(qd, qv, 0.0_kreal, 0.0_kreal, pressure, temp)

      ! compute time derivatives for each process
      dqrdt_autoconv = dqr_dt__autoconversion(ql, qg, rho_g)
      dqrdt_accre    = dqr_dt__accretion(ql, qg, rho_g, qr)
      dqldt_condevap = dqr_dt__condensation_evaporation(qv=qv, qr=qr, rho=rho, T=temp, p=pressure)
      dqldt_condevap = dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=temp, p=pressure)

      ! combine to create time derivatives for species
      dydt(idx_water_vapour) = -dqldt_condevap                                - dqrdt_condevap
      dydt(idx_cwater)       =  dqldt_condevap - dqrdt_autoconv - dqrdt_accre
      dydt(idx_rain)         =                   dqrdt_autoconv + dqrdt_accre + dqrdt_condevap

      dydt(idx_temp) = L_v/c_m*dqldt_condevap

   end function


   pure function rho_f(qd, qv, ql, qr, p, temp) result(rho)
      use microphysics_constants, only: R_v, R_d, rho_l => rho_w

      real(kreal), intent(in) :: qd, qv, ql, qr, p, temp
      real(kreal) :: rho, rho_inv

      rho_inv = (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l

      rho = 1.0_kreal/rho_inv
   end function rho_f

   pure function dqr_dt__autoconversion(ql, qg, rho_g)
      real(kreal), intent(in) :: ql, qg, rho_g
      real(kreal) :: dqr_dt__autoconversion

      real(kreal), parameter :: k_c = 1.0e-3, a_c = 5.0e-4

      ! TODO: what happens if ql < qg ?
      dqr_dt__autoconversion = k_c*(ql - qg/rho_g*a_c)
      dqr_dt__autoconversion = max(0.0, dqr_dt__autoconversion)
   end function dqr_dt__autoconversion

   pure function dqr_dt__accretion(ql, qg, rho_g, qr)
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: ql, qg, rho_g, qr
      real(kreal) :: dqr_dt__accretion

      real(kreal), parameter :: G3p5 = 3.32399614155  ! = Gamma(3.5)
      real(kreal), parameter :: N0r = 1.e7  ! [m^-4]
      real(kreal), parameter :: a_r = 201.0  ! [m^.5 s^-1]
      real(kreal), parameter :: rho0 = 1.12

      real(kreal) :: lambda_r

      lambda_r = (pi*(qg*rho_l)/(qr*rho_g)*N0r)**(1./4.)

      dqr_dt__accretion = pi/4.*N0r*a_r*sqrt(rho0/rho_g)*G3p5*lambda_r**(-3.5)*ql

      dqr_dt__accretion = max(0.0, dqr_dt__accretion)
   end function dqr_dt__accretion

   !> Condesation/evaporation of cloud-water droplets
   pure function dql_dt__condensation_evaporation(rho, rho_g, qv, ql, T, p)
      use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
      use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
      use microphysics_common, only: Ka_f => thermal_conductivity
      use microphysics_common, only: Dv_f => water_vapour_diffusivity
      use microphysics_constants, only: Lv => L_cond, R_v
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: rho, rho_g, qv, ql, T, p
      real(kreal) :: dql_dt__condensation_evaporation

      real(kreal), parameter :: r0 = 0.1e-6_kreal  ! initial cloud droplet radius
      real(kreal), parameter :: N0 = 200*1.0e6_kreal  ! initial cloud droplet number

      real(kreal) :: r_c, Nc
      real(kreal) :: pv_sat, qv_sat, Sw
      real(kreal) :: Dv, Fd
      real(kreal) :: Ka, Fk

      real(kreal), parameter :: r4_3 = 4.0_kreal/3.0_kreal
      real(kreal), parameter :: r1_3 = 1.0_kreal/3.0_kreal

      ! calculate saturation concentration of water vapour
      pv_sat = pv_sat_f(T)
      qv_sat = qv_sat_f(T, p)
      Sw = qv/qv_sat

      r_c = (ql*rho/(r4_3*pi*N0*rho_l))**r1_3

      if (Sw < 1.0) then
         if (r_c < r0) then
            ! don't allow evaporation if the droplets are smaller than the
            ! aerosol, there's nothing to evaporate then(!)
            r_c = 0.0
         else
         endif
      else
         r_c = max(r_c, r0)
      endif

      Nc = N0

      ! condensation evaporation of cloud droplets (given number of droplets
      ! and droplet radius calculated above)

      Ka = Ka_f(T)
      Fk = (Lv/(R_v*T) - 1._kreal) * (Lv/(Ka*T))*rho_l

      Dv = Dv_f(T, p)
      Fd = R_v*T/(pv_sat*Dv)*rho_l

      ! compute rate of change of condensate from diffusion
      dql_dt__condensation_evaporation = 4.*pi*rho_l/rho*Nc*r_c*(Sw - 1.0)/(Fk + Fd)

   end function dql_dt__condensation_evaporation

   !> Condensation and evaporation of rain. Similar to cloud-water droplet
   !> condensation/evaporation but includes corrections for "ventilation" and
   !> and droplet-size is assumed to follow a Marshall-Palmer distribution:
   !>
   !>     N(r)dr = N0 exp(-l*r) dr
   pure function dqr_dt__condensation_evaporation(qv, qr, rho, T, p)
      use microphysics_common, only: pv_sat_f => saturation_vapour_pressure
      use microphysics_common, only: qv_sat_f => saturation_vapour_concentration
      use microphysics_common, only: Ka_f => thermal_conductivity
      use microphysics_common, only: Dv_f => water_vapour_diffusivity
      use microphysics_common, only: dyn_visc_f => dynamic_viscosity
      use microphysics_constants, only: Lv => L_cond, R_v
      use microphysics_constants, only: pi, rho_l => rho_w

      real(kreal), intent(in) :: qv, qr, rho, T, p
      real(kreal) :: dqr_dt__condensation_evaporation

      real(kreal), parameter :: G2p75 = 1.608359421985546_kreal ! = Gamma(2.75)
        
      ! droplet-size distribution constant
      real(kreal), parameter :: N0 = 1.0e7 ! [m^-4]

      ! fall-speed coefficient taken from the r > 0.5mm expression for
      ! fall-speed from Herzog '98
      real(kreal), parameter :: a_r = 201.
      ! reference density
      real(kreal), parameter :: rho0 = 1.12

      real(kreal) :: pv_sat, qv_sat, Sw, l
      real(kreal) :: Dv, Fd
      real(kreal) :: Ka, Fk
      real(kreal) :: mu, f

      ! can't do cond/evap without any rain-droplets present
      if (qr == 0.0) then
         dqr_dt__condensation_evaporation = 0.0_kreal
      else
         ! computer super/sub-saturation
         qv_sat = qv_sat_f(T, p)
         Sw = qv/qv_sat

         ! size-distribtion length-scale
         l = (8.*rho_l*pi*N0/(qr*rho))**.25

         ! air condutivity and diffusion effects
         Ka = Ka_f(T)
         Fk = (Lv/(R_v*T) - 1)*Lv/(Ka*T)*rho_l

         pv_sat = pv_sat_f(T)
         Dv = Dv_f(T, p)
         Fd = R_v*T/(pv_sat*Dv)*rho_l

         ! compute the ventilation coefficient `f`
         ! dynamic viscosity
         mu = dyn_visc_f(T=T)

         f = 1. + 0.22*(2.*a_r*rho/mu)**.5*(rho0/rho)**.25*G2p75/(l**.25)

         ! compute rate of change of condensate from diffusion
         dqr_dt__condensation_evaporation = 4*pi*rho_l/rho*N0/l**2.*(Sw - 1.0)/(Fk + Fd)*f
      endif
   end function
end module mphys_no_ice
