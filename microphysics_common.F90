module microphysics_common

   implicit none

   contains
      function thermal_conductivity(temp)
         real :: thermal_conductivity
         real, intent(in) :: temp

         real, parameter :: adiffk=2.4e-2
         real, parameter :: bdiffk=8.e-5

         thermal_conductivity=adiffk+bdiffk*temp    
      end function thermal_conductivity

      function water_vapour_diffusivity(temp, pressure)
         use microphysics_constants, only: T0, ps0

         real :: water_vapour_diffusivity
         real, intent(in) :: temp, pressure

         real, parameter :: adiffd=2.11e-5
         real, parameter :: bdiffd=1.94

         water_vapour_diffusivity=adiffd*(temp/T0)**bdiffd*ps0/(pressure)
      end function water_vapour_diffusivity

      function saturation_vapour_pressure(temp)
         use microphysics_constants, only: epsmach

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

      function cv_mixture(q)
         use microphysics_register, only: n_species, n_moments__max, cv_all
         real, dimension(n_species, n_moments__max), intent(in) :: q
         real :: cv_mixture

         cv_mixture = sum(q(:,1)*cv_all)
      end function cv_mixture

      function cp_mixture(q)
         use microphysics_register, only: n_species, n_moments__max, cp_all
         real, dimension(n_species, n_moments__max), intent(in) :: q
         real :: cp_mixture

         cp_mixture = sum(q(:,1)*cp_all)
      end function cp_mixture

end module microphysics_common
