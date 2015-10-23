module microphysics_common
   use microphysics_constants, only: kreal

   implicit none

   public thermal_conductivity

   contains
      function thermal_conductivity(temp)
         use microphysics_constants, only: T0, a_K, b_K

         real(kreal) :: thermal_conductivity
         real(kreal), intent(in) :: temp

         thermal_conductivity=a_K*(temp-T0) + b_K
      end function thermal_conductivity

      function water_vapour_diffusivity(temp, pressure)
         use microphysics_constants, only: T0, ps0, a_D, b_D

         real(kreal) :: water_vapour_diffusivity
         real(kreal), intent(in) :: temp, pressure

         ! The tabulated values are in reference to P=10000. Pa
         water_vapour_diffusivity=a_D*(temp/T0)**b_D*10000./(pressure)
      end function water_vapour_diffusivity

      function saturation_vapour_pressure(temp)
         use microphysics_constants, only: epsmach, T0, p0vs, a0_lq, a1_lq, a0_ice, a1_ice

         real(kreal) :: saturation_vapour_pressure
         real(kreal), intent(in) :: temp

         real(kreal) :: expon2, satpw

         if (temp > T0) then
             expon2=a0_lq*(temp-t0)/max(temp+a1_lq,epsmach)
         else
             expon2=a0_ice*(temp-t0)/max(temp+a1_ice,epsmach)
         endif
         saturation_vapour_pressure=p0vs*exp(expon2)
      end function saturation_vapour_pressure

      function saturation_vapour_concentration(temp, p)
         use microphysics_constants, only: epsmach, T0, R_v, R_d

         real(kreal), intent(in) :: temp, p
         real(kreal) :: saturation_vapour_concentration

         real(kreal) :: pv_sat, eps

         pv_sat = saturation_vapour_pressure(temp)
         eps = R_d/R_v
         saturation_vapour_concentration = (eps*pv_sat)/(p-(1.-eps)*pv_sat)
      end function saturation_vapour_concentration


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
