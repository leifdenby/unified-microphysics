module mphys_no_ice
   implicit none
   public n_moments, init

   ! Constants specific to this microphysics implementation

   ! radwatc: Cloud particle radius (radius - water - cloud)
   real, parameter :: radwatc = 1.0

   ! One-moment schemes k = 1, mass
   ! Two moments, k = 0, 1 number concentration and mass
   ! Three moments, k = 0, 1, 2 number concentration, mass and radar reflectivity

   integer :: n_moments(1) = [1]

   contains
   subroutine init(hydrometeor_names)
      character(16), dimension(size(n_moments)), intent(out) :: hydrometeor_names
      hydrometeor_names(1)='cloud_water'
   end subroutine

   subroutine d_hydrometeors(q, dt)
      real, intent(in) :: q, dt

      contains
         function dq_cond_evap_cloudwater(temp, pressure, rgasnew, gasnew, wetnew, watcnew)
            use microphysics_common, only: thermal_conductivity, water_vapour_diffusivity, saturation_vapour_pressure
            use microphysics_common, only: R_v, L_cond, rho_w

            real :: dq_cond_evap_cloudwater
            real, intent(in) :: temp, pressure, rgasnew, gasnew, wetnew, watcnew

            real :: diffk, diffd, satpw
            real :: pfak
            real :: satww
            real :: supsatw
            real :: evap
            real :: xnwatc  ! number of cloud droplets per unit gas volume

            ! TODO: Move, this is a constant
            real, parameter :: r1 = 1.0, r4 = 4.0
            real, parameter :: xpi = 4.0/3.0*3.14
            real, parameter :: pi = 3.14

            diffk=thermal_conductivity(temp)
            diffd=water_vapour_diffusivity(temp, pressure)

            satpw=saturation_vapour_pressure(temp)

            pfak=rgasnew/(R_v*pressure)*gasnew
            satww=satpw*pfak

            supsatw=wetnew/satww-r1

            evap=supsatw/((L_cond/(R_v*temp)-r1)*L_cond/(diffk*temp) + R_v*temp/(diffd*satpw))

            xnwatc=watcnew/(rho_w*xpi*radwatc**3)

            dq_cond_evap_cloudwater=r4*pi*xnwatc*radwatc*evap

         end function dq_cond_evap_cloudwater

   end subroutine
end module mphys_no_ice
