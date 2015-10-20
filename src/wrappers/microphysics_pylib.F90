!> Wrapper which allows calling of microphysics routines from Python
module microphysics_pylib
   use microphysics_initialisation, only: mp_init => init
   use microphysics_register, only: cv_gases, q_flux_function
   use microphysics_common, only: cv_mixture, cp_mixture
   use microphysics_register, only: mp_n_solids => n_solids, mp_n_gases => n_gases, mp_n_moments__max => n_moments__max, n_species

   implicit none
contains
   
   subroutine init(microphysics_name)
      character(len=*), intent(in) :: microphysics_name

      call mp_init(microphysics_name)

   end subroutine init

   subroutine raise_python_exception(mesg)
      character(len=*), intent(in) :: mesg

      print *, "Please complain"
   end subroutine raise_python_exception

   subroutine calc_dq(q_g, q_tr, dt, temp, pressure, dq_g, dq_tr, n_gases, n_solids, n_moments__max)
      use microphysics_register, only: idx_cwater, idx_dry_air, idx_water_vapour
      use microphysics_constants, only: L_cond, kreal, kint

      integer, intent(in) :: n_gases, n_solids, n_moments__max
      real(kreal), intent(in) :: dt, temp, pressure
      real(kreal), dimension(n_gases), intent(in) :: q_g
      real(kreal), dimension(n_solids,n_moments__max), intent(in) :: q_tr
      real(kreal), dimension(n_gases), intent(out) :: dq_g
      real(kreal), dimension(n_solids,n_moments__max), intent(out) :: dq_tr

      if (.not. init_called()) then
         print *, "Error: `init` hasn't been called, please run `init(microphysics_name)` before calling any other function"
         call raise_python_exception('test')
      else if (.not. dimensions_valid(n_gases, n_solids, n_moments__max)) then
         print *, "Error: the argument does not have the correct number of dimensions"
         print *, "   shoud have: n_solids=", mp_n_solids, ", n_gases=", mp_n_gases, ", n_moments__max=", mp_n_moments__max
      else
         call q_flux_function(q_g, q_tr, dt, temp, pressure, dq_g, dq_tr)
      endif
   end subroutine calc_dq

   function dimensions_valid(n_gases, n_solids, n_moments__max) result(ok)
      integer, intent(in) :: n_gases, n_solids, n_moments__max
      logical :: ok

      if (n_gases == mp_n_gases .and. n_solids == mp_n_solids .and. mp_n_moments__max == n_moments__max) then
         ok = .true.
      else
         ok = .false.
      endif

   end function dimensions_valid

   function init_called()
      logical :: init_called

      if (n_species == 0) then
         init_called = .false.
      else
         init_called = .true.
      endif
   end function init_called

end module microphysics_pylib
