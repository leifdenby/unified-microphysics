module microphysics_common
   use microphysics_constants, only: kreal

   implicit none

   public thermal_conductivity

   contains
      function thermal_conductivity(temp)
         real(kreal) :: thermal_conductivity
         real(kreal), intent(in) :: temp

         real(kreal), parameter :: adiffk=2.4e-2
         real(kreal), parameter :: bdiffk=8.e-5

         thermal_conductivity=adiffk+bdiffk*temp    
      end function thermal_conductivity

      function water_vapour_diffusivity(temp, pressure)
         use microphysics_constants, only: T0, ps0

         real(kreal) :: water_vapour_diffusivity
         real(kreal), intent(in) :: temp, pressure

         real(kreal), parameter :: adiffd=2.11e-5
         real(kreal), parameter :: bdiffd=1.94

         water_vapour_diffusivity=adiffd*(temp/T0)**bdiffd*ps0/(pressure)
      end function water_vapour_diffusivity

      function saturation_vapour_pressure(temp)
         use microphysics_constants, only: epsmach, T0

         real(kreal) :: saturation_vapour_pressure
         real(kreal), intent(in) :: temp

         real(kreal), parameter :: awat=17.25
         real(kreal), parameter :: bwat=36.
         real(kreal), parameter :: asat=610.7

         real(kreal) :: expon2, satpw

         expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
         saturation_vapour_pressure=asat*exp(expon2)

      end function saturation_vapour_pressure

      function cv_mixture(q_gases, q_solids)
         use microphysics_register, only: n_gases, n_solids, n_moments__max, cv_solids, cv_gases
         real(kreal), dimension(n_solids, n_moments__max), intent(in) :: q_solids
         real(kreal), dimension(n_gases), intent(in) :: q_gases
         real(kreal) :: cv_mixture

         cv_mixture = sum(q_solids(:,1)*cv_solids) + sum(q_gases*cv_gases)
      end function cv_mixture

      function cp_mixture(q_gases, q_solids)
         use microphysics_register, only: n_gases, n_solids, n_moments__max, cp_solids, cp_gases
         real(kreal), dimension(n_solids, n_moments__max), intent(in) :: q_solids
         real(kreal), dimension(n_gases), intent(in) :: q_gases
         real(kreal) :: cp_mixture

         cp_mixture = sum(q_solids(:,1)*cp_solids) + sum(q_gases*cp_gases)
      end function cp_mixture

end module microphysics_common
