
module microphysics_integration
   use microphysics_register, only: q_flux_function
   use rungekuttafehlberg,    only: rkf45

   contains
      subroutine calc_dq(q_g, q_tr, pressure, temp, dt)
      end subroutine calc_dq

      !> Wrapper for calc_dq function function, for the integrator to use
      subroutine f(neqn,t,y,yp)
         integer, intent(in) :: neqn
         real(8), intent(in) :: t, y(neqn)
         real(8), intent(out) :: yp(neqn)

         yp(neqn) = sin(t)
      end subroutine f
end module microphysics_integration
