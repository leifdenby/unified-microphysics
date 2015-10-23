!
! For ice and liquid water we assume that the thermal expansion coefficient is
! zero so that the heat capacities at constant pressure and volume are the same
!> @todo: find sources for constants
module microphysics_constants
   integer, public, parameter :: real2 = selected_real_kind(3),            &
                                real4 = selected_real_kind(6),            &
                                real8 = selected_real_kind(12),           &
                                kreal = real8
   integer, public, parameter :: int4 = selected_int_kind(7),              &
                                int8 = selected_int_kind(13),             &
                                kint = int4

   ! specific heat capacity (at constant pressure) for dry air [J/kg/K]
   real(kreal), parameter :: cp_d = 1004.64_kreal
   ! specific heat capacity (at constant volume) for dry air [J/kg/K]
   real(kreal), parameter :: cv_d = 717.60_kreal

   ! specific heat capacity (at constant pressure) for water vapour [J/kg/K]
   real(kreal), parameter :: cp_v = 1864._kreal
   ! specific heat capacity (at constant volume) for water vapour [J/kg/K]
   real(kreal), parameter :: cv_v = 1402.55_kreal

   ! specific heat capacity (at constant pressure) for liquid water [J/kg/K]
   real(kreal), parameter :: cp_l = 4183._kreal
   ! specific heat capacity (at constant volume) for liquid water [J/kg/K]
   real(kreal), parameter :: cv_l = cp_l

   ! specific heat capacity (at constant pressure) for ice [J/kg/K]
   real(kreal), parameter :: cp_i = 2103._kreal
   ! specific heat capacity (at constant volume) for ice [J/kg/K]
   real(kreal), parameter :: cv_i = cp_i

   ! specific gas constant for water vapour [J/kg/K]
   real(kreal), parameter :: R_v = cp_v - cv_v
   real(kreal), parameter :: R_d = cp_d - cv_d

   ! latent heat of condensation [J/kg]
   real(kreal), parameter :: L_cond = 2500.8*1.0e3_kreal

   ! density of water [kg/m3]
   real(kreal), parameter :: rho_w = 1.0e3_kreal

   real(kreal), parameter :: T0 = 273.15_kreal
   real(kreal), parameter :: epsmach = 1.0e-15_kreal
   real(kreal), parameter :: ps0 = 101325.0

   real(kreal), parameter :: pi = 3.14

   ! Create a placeholder value that we can assign to variables for initial
   ! value, this way it will be clear if a variable hasn't been assigned
   ! properly
   COMPLEX, parameter :: nan = COMPLEX(0., 1.0)


   ! For computation of saturation vapour pressure with Teten's formula
   ! XXX: I think these are from ECHAM, but I'm not sure (Leif)
   real(kreal), parameter :: p0vs=610.7_kreal
   real(kreal), parameter :: a0_lq=17.25_kreal
   real(kreal), parameter :: a1_lq=36._kreal
   real(kreal), parameter :: a0_ice=21.875_kreal
   real(kreal), parameter :: a1_ice=7.66_kreal

end module microphysics_constants
