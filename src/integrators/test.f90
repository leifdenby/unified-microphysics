program test
   use rungekuttafehlberg, only: rkf45

   implicit none

   integer, parameter :: neqn = 1

   real(8) :: t, tout, relerr, abserr
   real(8) :: y(neqn)
   integer :: iflag
   integer :: iwork(5)
   double precision work(3+6*neqn)
   real(8), parameter :: pi = 3.141592653589793238462643

   print *, "Hello world!"
   y(1) = 0.5
   print *, y
   t = 0.0
   tout = 10*2.0*pi
   abserr = 1.0e-16

   call advance(f, neqn, y, t, tout, relerr, abserr, iflag)
   print *, y



   contains
      subroutine f(neqn,t,y,yp)
         integer, intent(in) :: neqn
         real(8), intent(in) :: t, y(neqn)
         real(8), intent(out) :: yp(neqn)

         yp(neqn) = sin(t)
      end subroutine f

      subroutine advance(f,neqn,y,t,tout,relerr,abserr,iflag)
      integer neqn,iflag
      double precision y(neqn),t,tout,relerr,abserr
      external f

      integer iwork(5)
      double precision work(3+6*neqn)
      logical finished
      integer triedcounter
      character*80 msg

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

 300     finished = .true.
 400     if (iflag .gt. 2) write (*,*) msg
 500  continue
      end do
      return
      end
end


