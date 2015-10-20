!> Simple microphysics implementation which causes instantenous saturation
!> adjustment. All excess water vapour is turned into cloud water and there is not
!> re-evaporation.
module mphys_saturation_adjustment
   use microphysics_register, only: register_incompressible_species, n_species, n_gases, n_solids
   use microphysics_constants, only: kreal, T0

   implicit none
   public init

   ! radwatc: Cloud particle radius (radius - water - cloud)
   real(kreal), parameter :: radwatc = 1.0


   contains
   subroutine init()
      call register_incompressible_species('cloud_water', 1)
   end subroutine

   subroutine calc_dq(q_g, q_tr, dt, temp, pressure, dq_g, dq_tr, dtemp)
      use microphysics_register, only: n_moments__max, idx_cwater, idx_dry_air, idx_water_vapour
      use microphysics_constants, only: L_cond

      real(kreal), intent(in) :: dt, temp, pressure
      real(kreal), dimension(n_gases), intent(in) :: q_g
      real(kreal), dimension(n_solids,n_moments__max), intent(in) :: q_tr

      real(kreal), intent(out) :: dtemp
      real(kreal), dimension(n_gases), intent(out) :: dq_g
      real(kreal), dimension(n_solids,n_moments__max), intent(out) :: dq_tr

      real(kreal) :: dq_cond = 0.0
      real(kreal) :: temp_local = 0.0

      dq_cond = dqdt_cond_cloudwater(temp, pressure, q_g(idx_water_vapour))*dt

      dtemp = L_cond*dq_cond
      dq_g(idx_water_vapour) = -dq_cond
      dq_tr(idx_cwater,1) = dq_cond

      contains
         !> Create cloudwater by consuming excess super saturation
         function dqdt_cond_cloudwater(temp, pressure, q_v) result(dq_cwater)
            use microphysics_common, only: saturation_vapour_pressure
            use microphysics_constants, only: R_v, R_d

            real(kreal), intent(in) :: temp, pressure, q_v
            real(kreal) :: dq_cwater
            real(kreal) :: satpw = 0.0, satww = 0.0

            satpw=saturation_vapour_pressure(temp)

            satww = R_d/R_v*satpw/pressure

            dq_cwater=min(q_v-satww, 0.0)
            
         end function dqdt_cond_cloudwater

   end subroutine

   !> Assuming contant pressure adjust temperature and saturation specific
   !> concentration of water vapor.
   !>
   !> Saturation water vapor pressure is calculated using Tetens (1930) formula.
   !!
   !! @param tem initial temperature [K]
   !! @param q_v initial specific concentration of water vapor [kg/kg]
   !! @param prs pressure [Pa]
   !! @returns tem temperature after adjustment [K]
   !! @return q_v saturation specification concentration of water vapor [kg/kg]
   subroutine moist_adjust( tem, q_v, prs )
      implicit none

      real(kreal), intent(inout) :: tem, q_v
      real(kreal), intent(in)    :: prs

      real(kreal) :: q_sat ! saturation water vapour over water or ice
      real(kreal) :: cor
      real(kreal) :: cond  ! condensate


      q_sat = lua( max( 180._kreal, min( tem, 370._kreal))) / prs

      q_sat = min( 0.5_dp, q_sat )
      cor   = 1._kreal / ( 1._kreal - ( Rv/Rd - 1._kreal ) * q_sat )
      q_sat = q_sat * cor

      cond  = (q_v - q_sat) /                                                &
      (1._kreal + q_sat * cor * lub( max(180._kreal,min( tem ,370._kreal))))

      cond  =  max( cond, 0._kreal )

      ! increase temperature through latent heat of condensation
      tem = tem + luc(tem) * cond

      ! remove condensed water from water vapor
      q_v  = q_v - cond

      if ( cond > 0._kreal ) then

         ! calculate specific humidity from saturation water vapor pressure using Tetens formula
         q_sat = lua( max( 180._kreal, min( tem, 370._kreal))) / prs

         q_sat = min( 0.5_dp, q_sat )
         cor   = 1._kreal / ( 1._kreal - ( Rv/Rd - 1._kreal ) * q_sat )
         q_sat = q_sat * cor

         cond  = (q_v - q_sat) /                                             &
         (1._kreal + q_sat * cor * lub( max(180._kreal,min(tem,370._kreal))) )

         cond  =  max( cond, 0._kreal )

         tem = tem + luc(tem) * cond

         q_v  = q_v - cond

      end if

   end subroutine moist_adjust

   real(kreal) function lua(T)
      use microphysics_constants, only: RD => R_d, RV => R_v
      real(kreal), parameter :: C1ES  = 610.78_kreal
      real(kreal), parameter :: C2ES  = C1ES*RD/RV 
      !> "a" constant for saturation vapor pressure over liquid water in Tetens equation
      real(kreal), parameter :: C3LES = 17.269_kreal 
      !> "b" constant for saturation vapor pressure over liquid in Tetens equation
      real(kreal), parameter :: C4LES = 35.86_kreal
      !> "a" constant for saturation vapor pressure over ice in Tetens equation
      real(kreal), parameter :: C3IES = 21.875_kreal
      !> "b" constant for saturation vapor pressure over ice in Tetens equation
      real(kreal), parameter :: C4IES = 7.66_kreal 

      ! lua(T) = es(T) * Rd/Rv
      ! in echam: tlucua

      real(kreal) :: T
      real(kreal) :: CVM3, CVM4 

      if ( T > T0 ) then
         CVM3 = C3LES
         CVM4 = C4LES
      else
         CVM3 = C3IES
         CVM4 = C4IES
      end if

      lua = c2es * exp( CVM3 * (T - T0)/(T - CVM4)) 

   end function lua


   real(kreal) function lub(T)
      use microphysics_constants, only: aLv => L_v
      !> latent heat of evaporation for water [J/kg]
      real(kreal), parameter :: aLv   = 2.5008e6_kreal      ! latent heat of vaporisation (at 0 C)        [J/kg]
      !> "a" constant for saturation vapor pressure over ice in Tetens equation
      real(kreal), parameter :: C3IES = 21.875_kreal
      !> "b" constant for saturation vapor pressure over ice in Tetens equation
      real(kreal), parameter :: C4IES = 7.66_kreal 
      !> "a" constant for saturation vapor pressure over liquid water in Tetens equation
      real(kreal), parameter :: C3LES = 17.269_kreal 
      !> "b" constant for saturation vapor pressure over liquid in Tetens equation
      real(kreal), parameter :: C4LES = 35.86_kreal

      real(kreal), parameter :: C5LES = C3LES*(T0-C4LES)
      real(kreal), parameter :: C5IES = C3IES*(T0-C4IES)

      real(kreal)       :: T
      real(kreal)       :: CVM4, CVM5

      if ( T > T0 ) then
         CVM4 = C4LES
         CVM5 = C5LES * ALV / CPD
      else
         CVM4 = C4IES
         CVM5 = C5IES * ALS / CPD
      end if

      lub = CVM5 / ( T - CVM4 )**2

   end function lub


   !> Calculate temperature dependent specific temperature change due to phase-transition.
   !> 
   !> If temperature is below freezing latent heat of freezing is used, otherwise
   !> latent heat of condensation is used.
   !! @param T temperature at which to calculate latent heat release [K]
   !! @return specific temperature change [K*kg/kg]
   real(kreal) function luc(T)
      ! latent heat

      implicit none

      real(kreal) :: T

      if ( T > T0 ) then
         ! latent heat of evaporation / dry air spec. heat cap. at contant pressure
         ! J/kg / ( J/kg/K ) = [K*kg/kg]
         luc = ALV / CPD
      else
         luc = ALS / CPD
      end if

   end function luc
end module mphys_saturation_adjustment
