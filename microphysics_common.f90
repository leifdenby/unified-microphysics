module microphysics_common

   ! specific heat capacity (at constant pressure) for water vapour [J/kg/K]
   real, parameter :: cp_v = 1.0
   ! specific heat capacity (at constant volume) for water vapour [J/kg/K]
   real, parameter :: cv_v = 1.0
   ! specific gas constant for water vapour [J/kg/K]
   real, parameter :: R_v = cp_v - cv_v

   ! latent heat of condensation [J/g]
   real, parameter :: L_cond = 2500.8

   ! density of water [kg/m3]
   real, parameter :: rho_w = 1.0e3

   contains
      function thermal_conductivity(temp)
         real :: thermal_conductivity
         real, intent(in) :: temp

         real, parameter :: adiffk=2.4e-2
         real, parameter :: bdiffk=8.e-5

         thermal_conductivity=adiffk+bdiffk*temp    
      end function thermal_conductivity

      function water_vapour_diffusivity(temp, pressure)
         real :: water_vapour_diffusivity
         real, intent(in) :: temp, pressure

         real, parameter :: adiffd=2.11e-5
         real, parameter :: bdiffd=1.94

         water_vapour_diffusivity=adiffd*(temp/T0)**bdiffd*ps0/(pressure)
      end function water_vapour_diffusivity

      function saturation_vapour_pressure(temp)
         real :: saturation_vapour_pressure
         real, intent(in) :: temp

         real, parameter :: awat=17.25
         real, parameter :: bwat=36.
         real, parameter :: t0=273.15
         real, parameter :: asat=610.7

         real :: expon2, satpw

         expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
         saturation_vapour_pressure=asat*exp(expon2)

      end function saturation_vapour_pressure
end module microphysics_common
