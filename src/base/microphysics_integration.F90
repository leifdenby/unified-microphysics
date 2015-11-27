
module microphysics_integration
   !use microphysics_register, only: q_flux_function
   !use rungekuttafehlberg,    only: rkf45
   use mphys_no_ice, only: dydt_mphys => dydt, n_species
   use microphysics_constants, only: abs_tol => integration_abs_tol, rel_tol => integration_rel_tol

   implicit none

   !public calc_dy, calc_dy_with_message
   public integrate_with_message, integrate_isometric, integrate_isobaric


   logical, parameter :: debug = .false.
   real(8), parameter :: dx_min = 1.0e-10
   integer, parameter :: max_steps = 10

   procedure(), pointer :: integrate => null()

   contains
      !!!!!!!!!!! Isometric
      pure function dydt_isometric(t, y) result(dydt)
         use microphysics_common, only: cv_mixture

         real(8), intent(in) :: t, y(5)
         real(8) :: c_m, dydt(5)

         c_m = cv_mixture(y)

         dydt = dydt_mphys(t, y, c_m)
      end function dydt_isometric

      pure function fix_y_isometric(y_old, dy) result(y_new)
         real(8), intent(in) :: y_old(5), dy(5)
         real(8) :: y_new(5)

         y_new = y_old + dy
      end function fix_y_isometric

      subroutine integrate_isometric(y, t0, t_end, msg_out)
         real(8), intent(inout) :: y(5)
         real(8), intent(in) :: t_end, t0
         character(len=100), optional :: msg_out
         character(len=100) :: msg

         integer :: m_total
         msg = " "
         m_total = 0

         call integrate_with_message(dydt_isometric, y, increment_state_isometric, t0, t_end, msg, m_total)

         if (msg(1:1) /= " ") then
            print *, "==============================="
            print *, "integration failed"
            print *, y
            print *, msg
            stop(0)
         endif

         if (present(msg_out)) then
            msg_out = msg
         endif
      end subroutine integrate_isometric

      !> Return a new state changed by increment `dy` keeping constant volume,
      !method is intended to guarantee self-consistency of new state `y_new`
      ! TODO: Implement mass-scaling
      pure function increment_state_isometric(y, dy) result(y_new)
         use microphysics_register, only: idx_pressure

         real(8), intent(in) :: y(5), dy(5)
         real(8) :: y_new(5), rho_old

         rho_old = density_eos(y)
         y_new = y + dy

         ! need to update pressure, since we assume constant volume density is
         ! unchanged, use old density to calculate update pressure
         y_new(idx_pressure) = pressure_eos(y_new, rho_old)

      end function increment_state_isometric
      !!!!!!!!!! end isometric

      !!!!!!!!!! isobaric
      subroutine integrate_isobaric(y, t0, t_end, msg_out)
         real(8), intent(inout) :: y(5)
         real(8), intent(in) :: t_end, t0
         character(len=100), optional :: msg_out

         character(len=100) :: msg
         integer :: m_total

         msg = " "
         m_total = 0

         call integrate_with_message(dydt_isobaric, y, increment_state_isobaric, t0, t_end, msg, m_total)

         if (msg(1:1) /= " ") then
            print *, "==============================="
            print *, "integration failed"
            print *, y
            print *, msg
            stop(0)
         endif

         if (present(msg_out)) then
            msg_out = msg
         endif

      end subroutine integrate_isobaric

      pure function dydt_isobaric(t, y) result(dydt)
         use microphysics_common, only: cp_mixture

         real(8), intent(in) :: t, y(5)
         real(8) :: c_m, dydt(5)

         c_m = cp_mixture(y)

         dydt = dydt_mphys(t, y, c_m)
      end function dydt_isobaric

      !> Return a new state changed by increment `dy`, method is intended to
      !guarantee self-consistency of new state `y_new`
      ! TODO: Implement mass-scaling
      pure function increment_state_isobaric(y, dy) result(y_new)
         real(8), intent(in) :: y(5), dy(5)
         real(8) :: y_new(5)

         y_new = y + dy
      end function increment_state_isobaric
      !!!!!!!!!!! end isobaric

      subroutine integrate_with_message(dydt_f, y, fix_y, t0, t_end, msg, m_total)
         real(8), intent(inout) :: y(5)
         real(8), intent(in) :: t_end, t0
         real(8) :: dt_s, t  ! sub-step timestep
         integer, intent(out) :: m_total

         interface
            pure function dydt_f(x, y)
               real(8), intent(in) :: x, y(5)
               real(8) :: dydt_f(5)
            end function
         end interface

         interface
            pure function fix_y(y_old, dy)
               real(8), intent(in) :: y_old(5), dy(5)
               real(8) :: fix_y(5)
            end function
         end interface

         character(len=100), intent(out) :: msg
         integer :: m
         !f2py raise_python_exception msg

         msg = " "

         ! make initial sub-cycling full timestep, in case we can do that
         dt_s = t_end - t0
         t = t0

         !print *, "Called from ATHAM", t0, t_end, dt_s

         do while (t < t_end)
            !print *, "Step start", t, dt_s
            if (t + dt_s > t_end) then
               dt_s = t_end - t
            endif

            ! evolve y, t and dt_s using the integrator
            ! dt_s will be the time-step that was actually used, so that we can
            ! use that for the next integration step

            if (debug) then
               print *, "=> t", t
            endif

            m = 0
            call rkf34_original(dydt_f, y, fix_y, t, dt_s, msg, m)
            m_total = m_total + m

            if (msg(1:1) /= " ") then
               exit
            endif

         enddo
      end subroutine integrate_with_message


      !> Compute pressure from equation of state.
      ! Because state representation doesn't contain density (it is diagnosed)
      ! we need to supply a density so that pressure can be calculated from the
      ! other state variables
      pure function pressure_eos(y, rho) result(p)
         use microphysics_register, only: idx_temp
         use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
         use microphysics_constants, only: R_v, R_d, rho_l => rho_w

         real(8), intent(in) :: y(5), rho

         real(8) :: temp, qd, qv, ql, qr
         real(8) :: p

         temp = y(idx_temp)

         ql = y(idx_cwater)
         qr = y(idx_rain)
         qv = y(idx_water_vapour)
         qd = 1.0 - ql - qr - qv

         p = temp*(qd*R_d + qv*R_v)/(1./rho - rho_l/ql)

      end function pressure_eos

      !> Calculate density of mixture from equation of state
      pure function density_eos(y) result(rho)
         use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
         use microphysics_register, only: idx_temp, idx_pressure
         use microphysics_constants, only: R_v, R_d, rho_l => rho_w

         real(8), intent(in) :: y(5)

         real(8) :: temp, p, qd, qv, ql, qr
         real(8) :: rho

         p = y(idx_pressure)
         temp = y(idx_temp)

         ql = y(idx_cwater)
         qr = y(idx_rain)
         qv = y(idx_water_vapour)
         qd = 1.0 - ql - qr - qv

         rho = 1.0/( (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l )
      end function density_eos


      recursive subroutine rkf34_original(dydx, y, fix_y, x, dx, msg, m)
         integer, parameter :: n = 5
         character(len=100), intent(out) :: msg
         real(8), intent(inout) :: y(n), dx
         real(8), intent(inout) :: x
         integer, intent(in) :: m


         interface
            pure function dydx(x, y)
               real(8), intent(in) :: x, y(5)
               real(8) :: dydx(5)
            end function
         end interface

         interface
            pure function fix_y(y_old, dy)
               real(8), intent(in) :: y_old(5), dy(5)
               real(8) :: fix_y(5)
            end function
         end interface

         real(8) :: max_abs_err
         real(8) :: dx_min__posdef
         real(8) :: max_rel_err, err

         !f2py raise_python_exception msg

         !real(8), parameter :: &
         !a2=2./7.,   b21=2./7., &
         !a3=7./15.,  b31=77./900.,   b32= 343.900, &
         !a4=35./38., b41=805./1444., b42=-77175./54872., b43=97125./54872., &
         !a5=1.0,     b51=79./490.,   b52= 0.0,           b53=2175./3616.,   b54=2166./9065.
         !real(8) :: &
         !c1_1=79./490.,    c2_1=0.0,  c3_1=2175./36.26,  c4_1=2166./9065.,&
         !c1_2=229./1470.,  c2_2=0.0,  c3_2=1125./1813.,   c4_2=13718./81585., c5_2=1./18.

         real(8), parameter :: &
            a2=0.5, b21=0.5, &
            a3=0.5, b31=0.0,        b32= 0.5, &
            a4=1.0, b41=0.0,        b42= 0.0,           b43=1.0, &
            a5=1.0, b51=1./6.,      b52= 1./3.,         b53=1./3.,   b54= 1./6.
         real(8), parameter :: &
            c1_1=0.0, c2_1=0.0,  c3_1=0.0,  c4_1=0.0,  c5_1=1.,&
            c1_2=1./6., c2_2=1./3.,    c3_2=1./3.,     c4_2=0.,   c5_2=1./6.

         real(8) :: k1(n), k2(n), k3(n), k4(n), k5(n)
         real(8) :: abs_err(n), y_n1(n), y_n2(n), s, rel_err(n)
         real(8) :: dydx0(n), max_total_err(n)

         logical :: done = .false.

         integer :: l

         if (debug) then
            print *, ""
            print *, " -- substep -- m dt", m, dx
            print *, 'y', y
         endif

         done = .false.

         dydx0 = dydx(x, y)

         ! TODO: use a vector for abs error here, so that we can have a
         ! different abs error for say specific concentration and temperature
         !
         ! State components that are already zero do not contribute to the
         ! relative error calculation
         dx_min__posdef = minval(abs(y/dydx0), y /= 0.0)
         if (dx_min__posdef > dx_min) then
            if (dx > dx_min__posdef) then
               print *, "Adjusting integration step down, too big to be pos def"
               print *, "dx dx_min__posdef, m", dx, dx_min__posdef, m
               print *, ""
               s = 0.5*dx_min__posdef/dx
               dx = s*dx
               if (debug) then
                  print *, "Initial timestep risks solution becoming negative, scaling timestep"
                  print *, "by s=", s
                  print *, "dx_new", dx
               endif
            endif
         else
            if (debug) then
               print *, "Timestep required for pos def is very small so we just take a single forward"
               print *, "Euler step before runge-kutta integration"
               print *, "dx dx_min, dx_min__posdef", dx, dx_min, dx_min__posdef
               print *, "y=", y
               print *, "dx__pd=", y/dydx0
               print *, "dy=", dx_min__posdef*dydx0
               print *, ""
            endif
            y = y + dx_min__posdef*dydx0
            x = x + dx_min__posdef

            if (debug) then
               print *, "y_new=", y
            endif

            ! TODO: There's got to be a better way than this, we want to make
            ! sure we get exactly to zero
            where (y < 1.0e-40)
               y = 0.0
            endwhere 
         endif

         if (.not. done) then
            k1 = 0.0
            k2 = 0.0
            k3 = 0.0
            k4 = 0.0
            k5 = 0.0

            k1 = dx*dydx(x,       y)
            k2 = dx*dydx(x+a2*dx, fix_y(y, b21*k1))
            k3 = dx*dydx(x+a3*dx, fix_y(y, b31*k1 + b32*k2))
            k4 = dx*dydx(x+a4*dx, fix_y(y, b41*k1 + b42*k2 + b43*k3))
            k5 = dx*dydx(x+a5*dx, fix_y(y, b51*k1 + b52*k2 + b53*k3 + b54*k4))

            y_n1 = fix_y(y, c1_1*k1 + c2_1*k2 + c3_1*k3 + c4_1*k4 + c5_1*k5)
            y_n2 = fix_y(y, c1_2*k1 + c2_2*k2 + c3_2*k3 + c4_2*k4 + c5_2*k5)

            !abs_err = abs(y_n1 - y_n2)
            abs_err = abs(1./6.*(k4 - k5))

            ! TODO: make abs_tol and rel_tol vectors
            max_total_err = (abs_tol + rel_tol*abs(y))
            s = 0.84*(minval(max_total_err/abs_err, abs_err > 0.0))**0.2

            if (any(y_n2 < 0.0)) then
               !msg = "Solution became negative"
               print *, "=> solution became negative"
               print *, "dx s", dx, s
               print *, "y", y
               print *, "y_n2", y_n2
               print *, "max_tot_err", max_total_err
               print *, "abs_err", abs_err

               if (s > 1.0) then
                  msg = "s was huge"
                  done = .true.
               endif
            endif

            if (debug) then
               print *, "max_tot_err=", max_total_err
               print *, "abs_err    =", abs_err
               !print *, "s_all=", max_total_err/abs_err
               !print *, "s=", s
               !print *, "dx=", dx
               !print *, abs_err < max_total_err

               !s = 0.84*(rel_tol*dx/max_abs_err)**0.25
               print *, ":: scaling by s", s
            endif

            if (any(isnan(y_n1)) .or. any(isnan(y_n2))) then
               if (debug) then
                  print *, "Solution became nan"
                  print *, "s dt", s, dx
                  print *, "y=", y
               endif
               if (s > 1.0) then
                  s = 0.1
               endif
               !msg = "solution became nan"
               !done = .true.
            else
               if (all(abs_err < max_total_err)) then
                  done = .true.
                  y = y_n2
                  x = x + dx
               endif 
            endif

            dx = dx*s

            if (.not. done) then
               if (m > max_steps) then
                  msg = "Didn't converge"
                  print *, "last s=", s
                  y = y_n2
               else
                  !if (s >= 1.0) then
                     !print *, "s=", s
                     !print *, "Warning: Incorrect scaling, timestep is growing."
                     !!s = 0.1
                  !endif
                  call rkf34_original(dydx, y, fix_y, x, dx, msg, m+1)
               endif
            endif
         endif

      end subroutine rkf34_original

      !subroutine calc_dy(y, dt, dy, n_variables)
         !real(8), intent(inout) :: y(n_variables)
         !real(8), intent(in) :: dt
         !real(8), intent(out) :: dy(n_variables)
         !integer, intent(in) :: n_variables
         !real(8) :: t0, tn, y0(n_variables)
         !real(8) :: relerr, abserr
         !integer :: iflag
         !character*80 :: msg

         !abserr = 10.0
         !relerr = 1.0e-3

         !t0 = 0.0
         !tn = t0 + dt
         !y0(:) = y
         !call advance(f, n_variables, y, t0, tn, abserr, relerr, iflag, msg)

         !dy = y - y0
      !end subroutine calc_dy

      !subroutine calc_dy_with_message(y, dt, dy, n_variables, msg)
         !real(8), intent(inout) :: y(n_variables)
         !real(8), intent(in) :: dt
         !real(8), intent(out) :: dy(n_variables)
         !integer, intent(in) :: n_variables
         !character*80, intent(out) :: msg
         !!f2py raise_python_exception msg
         !real(8) :: t0, tn, y0(n_variables)
         !real(8) :: relerr, abserr
         !integer :: iflag

         !abserr = 10.0
         !relerr = 1.0e-3

         !t0 = 0.0
         !tn = t0 + dt
         !y0(:) = y
         !call advance(f, n_variables, y, t0, tn, abserr, relerr, iflag, msg)

         !dy = y - y0
      !end subroutine calc_dy_with_message

      !!> Wrapper for calc_dq function function, for the integrator to use
      !!subroutine f(neqn,t,y,yp)
         !!integer, intent(in) :: neqn
         !!real(8), intent(in) :: t, y(neqn)
         !!real(8), intent(out) :: yp(neqn)

         !!yp(:) = 0.0

         !!! time isn't used
         !!call q_flux_function(y, neqn, yp)
      !!end subroutine f

      !!> Wrapper for rkf45 from `odespy`

      !subroutine advance(f,neqn,y,t,tout,relerr,abserr,iflag, msg)
      !integer neqn,iflag
      !double precision y(neqn),t,tout,relerr,abserr
      !external f

      !integer iwork(5)
      !double precision work(3+6*neqn)
      !logical finished
      !integer triedcounter
      !character*80, intent(out) :: msg

      !triedcounter = 0
      !finished = .false.
      !iflag = 1
      !msg = ''

      !do while (.not. finished)
         !call rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
         !triedcounter = triedcounter + 1
         !print *, t, y
!!     succeed
         !if (iflag .lt. 3) go to 300
!!     Tried more than two times. Return with failure.
         !if (triedcounter .gt. 2) go to 300
!!     Tried again when iflag = 3 or 4.
         !if (iflag .lt. 5) go to 500
!!     Invalid inputs
         !if (iflag .eq. 8) then
            !msg = 'Invalid input parameters'
            !go to 300
         !end if
!!     One-step integrator is altered as suggestion.
         !if (iflag .eq. 7) then
            !msg = 'One-step integrator is altered as suggestion'
            !iflag = -1
            !go to 400
         !end if
!!     Too much accuracy requirement.
         !if (iflag .eq. 6) then
            !msg = 'Too much accuracy requirement.'
            !go to 300
         !end if
!!     atol should be positive to continue.
         !msg = 'atol should be positive to continue.'
         !print *, "iflag", iflag

 !300     finished = .true.
 !400     if (iflag .gt. 2) write (*,*) msg
 !500  continue
      !end do
      !return
      !end
end module microphysics_integration
