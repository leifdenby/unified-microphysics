module mphys_water_vapour
   use microphysics_register, only: register_variable
   use microphysics_constants, only: kreal

   implicit none
   public init

   contains
   subroutine init()
   end subroutine

   pure function dydt(t, y, c_m)
      real(kreal), intent(in) :: t
      real(kreal), dimension(:), intent(in) :: y
      real(kreal), dimension(size(y)) :: dydt
      real(kreal), intent(in) :: c_m

   end function
end module mphys_water_vapour
