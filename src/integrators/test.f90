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

end


