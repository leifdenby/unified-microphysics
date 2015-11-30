module microphysics_common
   use microphysics_constants, only: kreal

   implicit none

   public thermal_conductivity

   contains
      pure function thermal_conductivity(temp)
         use microphysics_constants, only: T0, a_K, b_K

         real(kreal) :: thermal_conductivity
         real(kreal), intent(in) :: temp

         thermal_conductivity=a_K*(temp-T0) + b_K
      end function thermal_conductivity

      pure function water_vapour_diffusivity(temp, pressure)
         use microphysics_constants, only: T0, ps0, a_D, b_D

         real(kreal) :: water_vapour_diffusivity
         real(kreal), intent(in) :: temp, pressure

         ! The tabulated values are in reference to P=10000. Pa
         water_vapour_diffusivity=a_D*(temp/T0)**b_D*10000./(pressure)
      end function water_vapour_diffusivity

      pure function saturation_vapour_pressure(temp)
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

      pure function saturation_vapour_concentration(temp, p)
         use microphysics_constants, only: epsmach, T0, R_v, R_d

         real(kreal), intent(in) :: temp, p
         real(kreal) :: saturation_vapour_concentration

         real(kreal) :: pv_sat, eps

         pv_sat = saturation_vapour_pressure(temp)
         eps = R_d/R_v
         saturation_vapour_concentration = (eps*pv_sat)/(p-(1.-eps)*pv_sat)
      end function saturation_vapour_concentration

      pure function cp_gas(y)
         use microphysics_register, only: n_variables, cp_species, q_species_flag, idx_water_vapour
         use microphysics_constants, only: cp_d, cp_v
         real(kreal), dimension(n_variables), intent(in) :: y
         real(kreal) :: cp_gas, q_d

         q_d = 1.0_kreal - sum(y*q_species_flag)

         cp_gas = cp_d*q_d + cp_v*y(idx_water_vapour)/(q_d + y(idx_water_vapour))
      end function

      pure function cv_gas(y)
         use microphysics_register, only: n_variables, cv_species, q_species_flag, idx_water_vapour
         use microphysics_constants, only: cv_d, cv_v
         real(kreal), dimension(n_variables), intent(in) :: y
         real(kreal) :: cv_gas, q_d

         q_d = 1.0_kreal - sum(y*q_species_flag)

         cv_gas = (cv_d*q_d + cv_v*y(idx_water_vapour))/(q_d + y(idx_water_vapour))
      end function

      pure function cv_mixture(y)
         use microphysics_register, only: n_variables, cv_species, q_species_flag
         use microphysics_constants, only: cv_d
         real(kreal), dimension(n_variables), intent(in) :: y
         real(kreal) :: cv_mixture, q_d

         q_d = 1.0_kreal - sum(y*q_species_flag)

         cv_mixture = cv_d*q_d + sum(y*cv_species)
      end function cv_mixture

      pure function cp_mixture(y)
         use microphysics_register, only: n_variables, cp_species, q_species_flag
         use microphysics_constants, only: cp_d
         real(kreal), dimension(n_variables), intent(in) :: y
         real(kreal) :: cp_mixture, q_d

         q_d = 1.0_kreal - sum(y*q_species_flag)

         cp_mixture = cp_d*q_d + sum(y*cp_species)
      end function cp_mixture

end module microphysics_common
