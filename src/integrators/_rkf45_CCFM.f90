module rkf45_CCFM
   integer, parameter :: dp = 8
   integer, max_lev = 900

   contains

   subroutine calc_dy(y, dt, dy, n_variables)
      implicit none

      real(dp)                   :: cl_top               ! top of cloud [m]
      integer                    :: cl_top_lev           ! top of cloud [m]

      real(dp)                   :: err                  ! relative error estimate  [-]
      real(dp)                   :: err_v  (7)

      type(t_cldata)             :: y, yn1, yn2, err_y   ! cloud data
      type(t_cldata)             :: y_sol  (max_lev)     ! temp upper limit
      real(dp)                   :: x_sol  (max_lev)     ! 

      real(dp)                   :: x                    ! height of this level     [m]
      type(t_cldata)             :: dydx (5) 

      real(dp),dimension(level)  :: Hgt                  ! height (grid)            [m]

      !----- for debugging output -----------------------!
      real(dp)                   :: lwp, ice, rain, snow
  
  !--- initialisation --------------------------------------------------------!
  Hgt  = geog / g  
  
  x_sol = 0._dp
  y_sol = t_cldata( 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp, &
                    0._dp,0._dp,0._dp,0._dp,0._dp) 
  
  
  !---------------------------------------------------------------------------!
  !--- initialisation of cloud model - cloud characteristics for cloud base --!
 
  ! environmental values
  t_e  = T_ls (cl_base)
  q_v  = q_ls (cl_base)
  p_e  = pres(cl_base)
  
  dx   = ( geog(cl_base - 1) - geog(cl_base) ) / g
  dpdx = ( pres(cl_base - 1) - pres(cl_base) ) / dx
   
  ! saturation mixing ratio over water [kg/kg] 
  q_s = q_satw(t_e, p_e)
  
  ! density
  rho = p_e / Rd / t_e
 
  y   = t_cldata( 0._dp,0._dp,0._dp,0._dp,0._dp,0._dp, &
                  0._dp,0._dp,0._dp,0._dp,0._dp)

  y% tem = T_base               ! cloud temp_erature         [K]
  y% q_l = 0._dp                ! cloud liquid water         [kg/kg]
  y% q_i = 0._dp                ! cloud ice mixing ratio     [kg/kg]
  y% q_r = 0._dp                ! rain water                 [kg/kg]
  y% q_s = 0._dp                ! snow mixing ratio          [kg/kg]
  y% v_z = w_base               ! vertical velocity          [m/s]
  y% rad = r_base               ! cloud radius               [m]
  y% rho = rho                  ! incloud density            [kg/m^3]
  y% sig = r_base**2*pi/area    ! cloud cover fraction       [-]
  y% q_v = q_base               ! water vapour mixing ratio  [kg/kg]
  
 
  !------------------ initial condition --------------------------------------!
  i        = max_lev            ! max levels at the moment        [-]
  y_sol(i) = y                  ! 
  x        = geog(cl_base) / g  ! current height                  [m]
  cl_top   = x                  ! current cloud top               [m]
  x_sol(i) = x * g              ! geopotential of this cloud step [gpm]
  dx       = 10._dp             ! initial stepsize                [m]
  dx2      = dx / 2._dp
  dx6      = dx / 6._dp
  k        = 0

!======================= end:  debug out =====================================!

  !------------------ cloud model integration --------------------------------!

  dydx(1) = derive(x, y, dx2,  t_e, q_v, p_e, dpdx, cdnc, land)

  cloud_integration: do
 
  k = k+1
     ! Fehlberg trick: dydx(1)  = dydx(5) (when step accepted)
 
     t_e = interpolate(x + dx2, Hgt, T_ls, level)   ! temperature  env. (at height x+dx/2) [K]
     q_v = interpolate(x + dx2, Hgt, q_ls, level)   ! mixing ratio env. (at height x+dx/2) [kg/kg]
     p_e = interpolate(x + dx2, Hgt, Pres, level)   ! pressure     env. (at height x+dx/2) [Pa]
     
     p_n = interpolate(x + dx,  Hgt, Pres, level)   ! pressure     env. (at height x+dx)   [Pa]
     dpdx = ( p_n - p_e ) / dx * 2._dp              ! pressure derivative                  [Pa/m] 

     if ( y% v_z + dx2 * dydx(1)% v_z < 0.1_dp ) exit cloud_integration
 
     dydx(2) = derive(x + dx2, y + dx2 * dydx(1), dx2, t_e, q_v, p_e, dpdx, cdnc, land)

     if ( y% v_z + dx2 * dydx(2)% v_z < 0.1_dp ) exit cloud_integration
     
     dydx(3) = derive(x + dx2, y + dx2 * dydx(2), dx2, t_e, q_v, p_e, dpdx, cdnc, land)
     
     if ( y% v_z + dx  * dydx(3)% v_z < 0.1_dp ) exit cloud_integration
    
     dydx(4) = derive(x + dx,  y + dx  * dydx(3), dx,  t_e, q_v, p_e, dpdx, cdnc, land)   
    
     yn1 = y + ( dydx(1) + 2._dp * ( dydx(2) +  dydx(3) ) + dydx(4)           ) * dx/6._dp

     dydx(5) = derive(x + dx,  yn1,                     dx,       t_e, q_v, p_e, dpdx, cdnc, land)
   
     yn2 = y + ( dydx(1) + 2._dp * ( dydx(2) +  dydx(3) )           + dydx(5) ) * dx6
         
     ! absolute error estimate:
     err_y = ( dydx(4) - dydx(5) ) * dx6
        
     ! init error vector
     err_v = 0._dp
   
     ! relative error estimate:  
                              err_v(1) = err_y%tem / yn2%tem
                              err_v(2) = err_y%rad / yn2%rad    
                              err_v(3) = err_y%v_z / yn2%v_z
     if ( yn2%q_l > 1.e-15 )  err_v(4) = err_y%q_l / yn2%q_l
     if ( yn2%q_i > 1.e-15 )  err_v(5) = err_y%q_i / yn2%q_i
     if ( yn2%q_r > 1.e-15 )  err_v(6) = err_y%q_r / yn2%q_r
     if ( yn2%q_s > 1.e-15 )  err_v(7) = err_y%q_s / yn2%q_s
    
     err = sqrt( dot_product(err_v, err_v) )
   
     ! if error within tolerance accept the solution
     if ( err < tol ) then
       
         i = i - 1   ! level counter for adaptive cloudmodel 

         if ( i < 1 ) then
            print *, 'not enough memory space'
            print *, 'height: ',x,' v_z: ', y%v_z
            exit cloud_integration
         end if 
          
         if ( y%q_l < 0._dp ) then
            exit cloud_integration
         end if 

         x = x + dx        ! height - height of this step [m]
            
         cl_top = x        ! height of cloudtop [m]
         cl_top_lev = i
   
         x_sol(i) = x * g  ! height of this cloud step [gpm]
            
         ! cloud model solution
         y        = yn2    ! or yn1

         ! I think the low cloud water case must be fixed here 
         ! it's pretty resistent in the proper Runge-Kutta integration  
         ! a first (and simple) try
         if ( y%q_l < 1.e-09_dp  .and. &
              y%q_l > 1.e-15_dp  .and. &
              y%tem < 260._dp          )  then
            
            ! a correction for the temperatur should be added
            y% q_i = y%q_i + y%q_l - 0.5e-15_dp
            y% q_l = 0.5e-15_dp
            
            dydx(5)% q_l = 0._dp

         end if 
                          
         y_sol(i) = y

         ! Fehlberg trick:
         dydx(1) = dydx(5)
  
         if ( y% v_z < 0.2_dp )  exit cloud_integration           


     end if  ! if error within tolerance accept the solution
    
     err = max(err, 0.001_dp*tol)
 
     ! new stepsize:
     dx = min( q*dx, sqrt(sqrt((0.8_dp*tol/err))) * dx )
     dx = min( dx, dx_max )
     
     dx2 = dx / 2._dp
     dx6 = dx / 6._dp

  end do cloud_integration


  !============== restriction of cloud data to echam levels ==================!     
  !--- initialisation:
  cl_data(1:level) = t_cldata( 0._dp,0._dp,0._dp,0._dp,0._dp,&
                               0._dp,0._dp,0._dp,0._dp,0._dp,&
                               0._dp)!,0._dp,0._dp,0._dp,0._dp )
  cl_info          = t_clinfo( level,.false.,0._dp,0._dp)
  
  ! set cloud type existence flag, if cloud reaches at least the next level
  if ( cl_top > Hgt(cl_base-1) ) then 
  
     call restrict_cldata( level,                                       &
                           y_sol, x_sol, max_lev, cl_top_lev, cl_top,   & !  in
                           geog, pres, cl_base, area,                   & !  in
                           cl_data, cl_info                             ) ! out 
  end if
  !========= end: restriction of cloud data to echam levels ==================!

!============================= debug output ==================================!
if (.true.) then
#ifdef AUTONOM
    write(45,   '(a8, a9  , a13  , 2a11  , a10  , 5a15)')                     &
              '# Level','Height','Temp cloud','Temp env','velocity',          &
              'Radius','humidity','rainwater',                                & 
              'liquid water','cloud ice','cloud snow'
    do i = cl_info%top, cl_base
       write(45,'(i8, f9.1, f13.3, 2f11.3, f10.1, 5es15.5)')                  &
                i, geog(i)/g, cl_data(i)%tem , T_ls(i), cl_data(i)%v_z,       & 
                cl_data(i)%rad, cl_data(i)%q_v, cl_data(i)%q_r,               &
                cl_data(i)%q_l, cl_data(i)%q_i, cl_data(i)%q_s
    end do

    ! some water information 
    lwp  = 0.
    ice  = 0.
    rain = 0.     
    snow = 0.
     
    ! liquid water path = sum( q_l * rho * dx )  [kg/m^2]
    do i = cl_info%top, cl_base
       if ( cl_data(i)%Tem > 10. ) then
          dx   = ( geog(i-1) - geog(i) ) / g                        
          rho  = pres(i) / RD / cl_data(i)%tem 
          lwp  = lwp  +  cl_data(i)%q_l * rho * dx 
          rain = rain +  cl_data(i)%q_r * rho * dx
          ice  = ice  +  cl_data(i)%q_i * rho * dx 
          snow = snow +  cl_data(i)%q_s * rho * dx  
       end if
    end do
      
    print *,'rain : ', rain, '[kg/m^2]'
    print *,'snow : ', snow, '[kg/m^2]'
    print *,'lwp  : ', lwp , '[kg/m^2]'

!    print *,'ice    water path: ', ice , '[kg/m^2]'
!    print *,'lwp  + iwp       : ', lwp  + ice  , '[kg/m^2]'      
!    print *,'rain + snow      : ', rain + snow , '[kg/m^2]'
!    print *,'total            : ', lwp  + ice + rain + snow
#endif
end if
!============================ end: debug output ==============================!

end subroutine cloudmodel

end module rkf45_CCFM
