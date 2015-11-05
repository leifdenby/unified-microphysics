module mphys_water_vapour
   use microphysics_register, only: register_variable
   use microphysics_constants, only: kreal

   implicit none
   public init

   contains
   subroutine init()
   end subroutine

   pure function dydt(t, y)
      integer, parameter :: n_variables = 2
      real(kreal), intent(in) :: t
      real(kreal), dimension(n_variables), intent(in) :: y
      real(kreal), dimension(n_variables) :: dydt

   end function
end module mphys_water_vapour
