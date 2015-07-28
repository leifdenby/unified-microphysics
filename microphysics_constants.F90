module microphysics_constants

   ! specific heat capacity (at constant pressure) for dry air [J/kg/K]
   real, parameter :: cp_d = 1.0
   ! specific heat capacity (at constant volume) for dry air [J/kg/K]
   real, parameter :: cv_d = 1.0

   ! specific heat capacity (at constant pressure) for water vapour [J/kg/K]
   real, parameter :: cp_v = 1.0
   ! specific heat capacity (at constant volume) for water vapour [J/kg/K]
   real, parameter :: cv_v = 1.0

   ! specific heat capacity (at constant pressure) for liquid water [J/kg/K]
   real, parameter :: cp_l = 1.0
   ! specific heat capacity (at constant volume) for liquid water [J/kg/K]
   real, parameter :: cv_l = 1.0

   ! specific heat capacity (at constant pressure) for ice [J/kg/K]
   real, parameter :: cp_i = 1.0
   ! specific heat capacity (at constant volume) for ice [J/kg/K]
   real, parameter :: cv_i = 1.0

   ! specific gas constant for water vapour [J/kg/K]
   real, parameter :: R_v = cp_v - cv_v
   real, parameter :: R_d = cp_d - cv_d

   ! latent heat of condensation [J/g]
   real, parameter :: L_cond = 2500.8

   ! density of water [kg/m3]
   real, parameter :: rho_w = 1.0e3

   real, parameter :: T0 = 273.15
   real, parameter :: epsmach = 1.0e-15
   real, parameter :: ps0 = 101325.0
end module microphysics_constants
