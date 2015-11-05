
module microphysics_integration
   !use microphysics_register, only: q_flux_function
   !use rungekuttafehlberg,    only: rkf45
   use mphys_no_ice, only: dydt, n_species

   implicit none

   !public calc_dy, calc_dy_with_message

   contains
      subroutine integrate(y, t0, t_end)
         real(8), intent(inout) :: y(n_species+3)
         real(8), intent(in) :: t_end, t0
         real(8) :: dt_s, t  ! sub-step timestep

         character(len=100) :: msg
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
            call rkf34_original(dydt, y, t, dt_s, msg, 0)

            if (msg(1:1) /= " ") then
               print *, "==============================="
               print *, "integration failed"
               print *, y
               print *, msg
               stop(0)
            endif
         enddo
      end subroutine integrate

      recursive subroutine rkf34_original(dydx, y, x, dx, msg, m)
         integer, parameter :: n = n_species + 3
         character(len=100), intent(out) :: msg
         real(8), intent(inout) :: y(n), dx
         real(8), intent(inout) :: x
         integer, intent(in) :: m

         interface
            function dydx(x, y)
               real(8), intent(in) :: x, y(5)
               real(8) :: dydx(5)
            end function
         end interface

         real(8) :: max_abs_err

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
         real(8) :: abs_err(n), y_n1(n), y_n2(n), s

         ! TODO: move these into an input variable
         real(8) :: abs_tol=1.0e-3, rel_tol=1.0e-2, max_rel_err, err

         logical :: done = .false.

         !print *, "Substep, dx=", dx

         done = .false.

         k1 = 0.0
         k2 = 0.0
         k3 = 0.0
         k4 = 0.0
         k5 = 0.0

         k1 = dx*dydx(x,       y)
         k2 = dx*dydx(x+a2*dx, y+b21*k1)
         k3 = dx*dydx(x+a3*dx, y+b31*k1 + b32*k2)
         k4 = dx*dydx(x+a4*dx, y+b41*k1 + b42*k2 + b43*k3)
         k5 = dx*dydx(x+a5*dx, y+b51*k1 + b52*k2 + b53*k3 + b54*k4)

         y_n1 = y + c1_1*k1 + c2_1*k2 + c3_1*k3 + c4_1*k4 + c5_1*k5
         y_n2 = y + c1_2*k1 + c2_2*k2 + c3_2*k3 + c4_2*k4 + c5_2*k5

         !abs_err = abs(y_n1 - y_n2)
         abs_err = abs(1./6.*(k4 - k5))

         max_abs_err = maxval(abs_err)
         max_rel_err = maxval(abs_err/y_n2)

         s = 0.84*(rel_tol*dx/max_abs_err)**0.25

         if (any(isnan(y_n1)) .or. any(isnan(y_n2))) then
            print *, "Warning: Solution became nan"
         else
            if (max_rel_err < rel_tol) then
               done = .true.
            endif 
         endif

         if (done) then
            y = y_n2
            x = x + dx

            ! if the suggested scaling makes the step larger then lets use it
            !if (s > 1.0) then
               !dx = dx*s
            !endif
         else
            if (dx < 1.0e-32) then
               msg = "step size became very small"
            else if (m > 1000) then
               msg = "Didn't converge"
               print *, "last s=", s
               y = y_n2
            else
               if (s >= 1.0) then
                  print *, "s=", s
                  print *, "Warning: Incorrect scaling, timestep is growing."
                  s = 0.1
               endif
               dx = dx*s
               call rkf34_original(dydx, y, x, dx, msg, m+1)
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
