!> Wrapper which allows calling of microphysics routines from Python
module microphysics_pylib
   use microphysics_initialisation, only: mp_init => init_with_message

   implicit none
contains
   
   subroutine init(microphysics_name, system_constraint, error_mesg)
      character(len=*), intent(in) :: microphysics_name, system_constraint
      character(len=100), intent(out) :: error_mesg
      !f2py raise_python_exception error_mesg

      call mp_init(microphysics_name, system_constraint, error_mesg)

   end subroutine init

   subroutine dydt(y, pressure, dydt_, n_variables, error_mesg)
      use microphysics_register, only: idx_cwater, idx_water_vapour
      use microphysics_constants, only: L_cond, kreal, kint
      !f2py raise_python_exception error_mesg

      real(kreal), intent(in) :: pressure
      integer, intent(in) :: n_variables
      real(kreal), dimension(n_variables), intent(in) :: y
      real(kreal), dimension(n_variables), intent(out) :: dydt_

      character(len=100), intent(out) :: error_mesg

      if (.not. init_called()) then
         error_mesg = "`init` hasn't been called, please run `init(microphysics_name)` before calling any other functions"
      !else
         !call q_flux_function(y, pressure, dydt_)
      endif
   end subroutine dydt

   subroutine integrate_microphysics(y, t, t_end, n_variables, error_mesg)
      use microphysics_integration, only: integrate
      use microphysics_constants, only: kreal

      integer, intent(in) :: n_variables
      real(kreal), intent(in) :: t, t_end, y(n_variables)
      character(len=100), intent(out) :: error_mesg
      !f2py raise_python_exception error_mesg

      call integrate(y, t, t_end, error_mesg, n_variables)
   end subroutine

   function init_called()
      use microphysics_register, only: n_variables

      logical :: init_called

      if (n_variables == 0) then
         init_called = .false.
      else
         init_called = .true.
      endif
   end function init_called

   function mixture_heat_capacity(y, n_variables) result(c_m)
      use microphysics_common, only: cv_mixture, cp_mixture
      use microphysics_register, only: model_constraint, &
                                       MODEL_CONSTRAINT_ISOBARIC, MODEL_CONSTRAINT_ISOMETRIC
      use microphysics_constants, only: kreal

      integer :: n_variables
      real(kreal), intent(in) :: y(n_variables)
      real(kreal) :: c_m

      if (model_constraint == MODEL_CONSTRAINT_ISOBARIC) then
         c_m = cp_mixture(y)
      else if (model_constraint == MODEL_CONSTRAINT_ISOMETRIC) then
         c_m = cv_mixture(y)
      else
         print *, "`mixture_heat_capacity`: Not implemented"
         call exit(-1)
      endif

   end function

end module microphysics_pylib
