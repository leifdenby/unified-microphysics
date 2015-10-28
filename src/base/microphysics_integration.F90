
module microphysics_integration
   use microphysics_register, only: q_flux_function
   use rungekuttafehlberg,    only: rkf45

   implicit none

   public calc_dy, calc_dy_with_message

   contains
      subroutine calc_dy(y, dt, dy, n_variables)
         real(8), intent(inout) :: y(n_variables)
         real(8), intent(in) :: dt
         real(8), intent(out) :: dy(n_variables)
         integer, intent(in) :: n_variables
         real(8) :: t0, tn, y0(n_variables)
         real(8) :: relerr, abserr
         integer :: iflag
         character*80 :: msg

         abserr = 10.0
         relerr = 1.0e-3

         t0 = 0.0
         tn = t0 + dt
         y0(:) = y
         call advance(f, n_variables, y, t0, tn, abserr, relerr, iflag, msg)

         dy = y - y0
      end subroutine calc_dy

      subroutine calc_dy_with_message(y, dt, dy, n_variables, msg)
         real(8), intent(inout) :: y(n_variables)
         real(8), intent(in) :: dt
         real(8), intent(out) :: dy(n_variables)
         integer, intent(in) :: n_variables
         character*80, intent(out) :: msg
         !f2py raise_python_exception msg
         real(8) :: t0, tn, y0(n_variables)
         real(8) :: relerr, abserr
         integer :: iflag

         abserr = 10.0
         relerr = 1.0e-3

         t0 = 0.0
         tn = t0 + dt
         y0(:) = y
         call advance(f, n_variables, y, t0, tn, abserr, relerr, iflag, msg)

         dy = y - y0
      end subroutine calc_dy_with_message

      !> Wrapper for calc_dq function function, for the integrator to use
      subroutine f(neqn,t,y,yp)
         integer, intent(in) :: neqn
         real(8), intent(in) :: t, y(neqn)
         real(8), intent(out) :: yp(neqn)

         yp(:) = 0.0

         ! time isn't used
         call q_flux_function(y, neqn, yp)
      end subroutine f

      !> Wrapper for rkf45 from `odespy`

      subroutine advance(f,neqn,y,t,tout,relerr,abserr,iflag, msg)
      integer neqn,iflag
      double precision y(neqn),t,tout,relerr,abserr
      external f

      integer iwork(5)
      double precision work(3+6*neqn)
      logical finished
      integer triedcounter
      character*80, intent(out) :: msg

      triedcounter = 0
      finished = .false.
      iflag = 1
      msg = ''

      do while (.not. finished)
         call rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
         triedcounter = triedcounter + 1
         print *, t, y
!     succeed
         if (iflag .lt. 3) go to 300
!     Tried more than two times. Return with failure.
         if (triedcounter .gt. 2) go to 300
!     Tried again when iflag = 3 or 4.
         if (iflag .lt. 5) go to 500
!     Invalid inputs
         if (iflag .eq. 8) then
            msg = 'Invalid input parameters'
            go to 300
         end if
!     One-step integrator is altered as suggestion.
         if (iflag .eq. 7) then
            msg = 'One-step integrator is altered as suggestion'
            iflag = -1
            go to 400
         end if
!     Too much accuracy requirement.
         if (iflag .eq. 6) then
            msg = 'Too much accuracy requirement.'
            go to 300
         end if
!     atol should be positive to continue.
         msg = 'atol should be positive to continue.'
         print *, "iflag", iflag

 300     finished = .true.
 400     if (iflag .gt. 2) write (*,*) msg
 500  continue
      end do
      return
      end
end module microphysics_integration
