!> Wrapper which allows calling of microphysics routines from Python
module microphysics_pylib
   use microphysics_initialisation, only: mp_init => init
   use microphysics_register, only: q_flux_function, n_variables

   implicit none
contains
   
   subroutine init(microphysics_name)
      character(len=*), intent(in) :: microphysics_name

      call mp_init(microphysics_name)

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
      else
         call q_flux_function(y, pressure, dydt_)
      endif
   end subroutine dydt

   function init_called()
      logical :: init_called

      if (n_variables == 0) then
         init_called = .false.
      else
         init_called = .true.
      endif
   end function init_called

end module microphysics_pylib
