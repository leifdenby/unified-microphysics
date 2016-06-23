!> Simple microphysics implementation which only supports the creation of cloud
!> water from water vapour, no ice-phases or precipitation creation is
!> implemented

module mphys_warm_no_rain
   use microphysics_register, only: register_variable
   use microphysics_constants, only: kreal
   use microphysics_register, only: n_variables

   implicit none

   public init, rho_f

   contains
   subroutine init()
      call register_variable('cloud_water', 1)
      ! XXX: shouldn't need to allocate this variable too, but somewhere there's
      ! some code that expects it to be around and causes a crash. I haven't got
      ! time to fix it right now
      call register_variable('rain', 1)
      print *, "Note: rain (autoconversion and accretion) disabled"
   end subroutine

   pure function dydt(t, y, c_m)
      use microphysics_register, only: idx_cwater, idx_water_vapour, idx_temp, idx_pressure
      use microphysics_constants, only: L_v => L_cond

      ! XXX: required for f2py compilation, but fixes the number of vars...
      !integer, parameter :: n_variables = 4

      real(kreal), dimension(n_variables), intent(in) :: y
      real(kreal), intent(in) :: t, c_m
      real(kreal) :: dydt(n_variables)

      real(kreal) :: ql, qv, qg, qr, qd
      real(kreal) :: rho, rho_g
      real(kreal) :: dqrdt_autoconv, dqrdt_accre, dqldt_condevap
      real(kreal) :: temp, pressure

      ! OBS: it's important to make sure that return variable is initiated to
      ! zero
      dydt = 0.0

      temp = y(idx_temp)
      pressure = y(idx_pressure)

      ! pick out specific concentrations from state vectors
      ql = y(idx_cwater)
      qr = 0.0_kreal
      qv = y(idx_water_vapour)
      qd = 1.0_kreal - ql - qr - qv
      qg = qv + qd

      ! compute gas and mixture density using equation of state
      rho = rho_f(qd, qv, ql, qr, pressure, temp)
      rho_g = rho_f(qd, qv, 0.0_kreal, 0.0_kreal, pressure, temp)

      ! compute time derivatives for each process
      dqrdt_autoconv = 0.0
      dqrdt_accre    = 0.0
      dqldt_condevap = dql_dt__condensation_evaporation(rho=rho, rho_g=rho_g, qv=qv, ql=ql, T=temp, p=pressure)

      ! combine to create time derivatives for species
      dydt(idx_water_vapour) = -dqldt_condevap
      dydt(idx_cwater)       =  dqldt_condevap

      dydt(idx_temp) = L_v/c_m*dqldt_condevap

   end function


   pure function rho_f(qd, qv, ql, qr, p, temp) result(rho)
      use microphysics_constants, only: R_v, R_d, rho_l => rho_w

      real(kreal), intent(in) :: qd, qv, ql, qr, p, temp
      real(kreal) :: rho, rho_inv

      rho_inv = (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l

      rho = 1.0_kreal/rho_inv
   end function rho_f

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
end module mphys_warm_no_rain
