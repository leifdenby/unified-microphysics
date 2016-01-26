
!> Helper functions which constrain the state to be physical during multi-step
!> integration

!> Integration helpers for isometric (constant volume integration)
module isometric_integration_helpers
   use mphys_no_ice, only: dydt_mphys => dydt, n_species
   use integrators, only: integrate_with_message

   implicit none

   public integrate_isometric

   contains
      subroutine integrate_isometric(y, t0, t_end, msg_out, m_total, n)
         integer, intent(in) :: n
         real(8), intent(inout), dimension(n) :: y
         real(8), intent(in) :: t_end, t0
         integer, intent(inout) :: m_total
         character(len=100), intent(inout), optional :: msg_out

         call integrate_with_message(dydt_isometric, y, increment_state_isometric, t0, t_end, msg_out, m_total)
      end subroutine integrate_isometric

      pure function dydt_isometric(t, y) result(dydt)
         use microphysics_common, only: cv_mixture

         real(8), intent(in) :: t
         real(8), intent(in), dimension(:) :: y
         real(8) :: c_m, dydt(size(y))

         c_m = cv_mixture(y)

         dydt = dydt_mphys(t, y, c_m)
      end function dydt_isometric

      pure function fix_y_isometric(y_old, dy) result(y_new)
         real(8), intent(in), dimension(:) :: y_old, dy
         real(8), dimension(size(y_old)) :: y_new

         y_new = y_old + dy
      end function fix_y_isometric

      !> Return a new state changed by increment `dy` keeping constant volume,
      !method is intended to guarantee self-consistency of new state `y_new`
      ! TODO: Implement mass-scaling
      function increment_state_isometric(y, dy) result(y_new)
         use microphysics_register, only: idx_pressure

         real(8), intent(in), dimension(:) :: y, dy
         real(8) :: y_new(size(y)), rho_old

         rho_old = density_eos(y)
         y_new = y + dy

         ! need to update pressure, since we assume constant volume density is
         ! unchanged, use old density to calculate update pressure
         y_new(idx_pressure) = pressure_eos(y_new, rho_old)

      end function increment_state_isometric

      !> Calculate density of mixture from equation of state
      pure function density_eos(y) result(rho)
         use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
         use microphysics_register, only: idx_temp, idx_pressure
         use microphysics_constants, only: R_v, R_d, rho_l => rho_w

         real(8), intent(in), dimension(:) :: y

         real(8) :: temp, p, qd, qv, ql, qr
         real(8) :: rho

         p = y(idx_pressure)
         temp = y(idx_temp)

         ql = y(idx_cwater)
         qr = y(idx_rain)
         qv = y(idx_water_vapour)
         qd = 1.0 - ql - qr - qv

         rho = 1.0/( (qd*R_d + qv*R_v)*temp/p + (ql+qr)/rho_l )
      end function density_eos

      !> Compute pressure from equation of state.
      ! Because state representation doesn't contain density (it is diagnosed)
      ! we need to supply a density so that pressure can be calculated from the
      ! other state variables
      pure function pressure_eos(y, rho) result(p)
         use microphysics_register, only: idx_temp
         use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain
         use microphysics_constants, only: R_v, R_d, rho_l => rho_w

         real(8), intent(in) :: rho
         real(8), intent(in), dimension(:) :: y

         real(8) :: temp, qd, qv, ql, qr
         real(8) :: p

         temp = y(idx_temp)

         ql = y(idx_cwater)
         qr = y(idx_rain)
         qv = y(idx_water_vapour)
         qd = 1.0 - ql - qr - qv

         p = temp*(qd*R_d + qv*R_v)/(1./rho - ql/rho_l)

      end function pressure_eos
end module

module isobaric_integration_helpers
   use mphys_no_ice, only: dydt_mphys => dydt, n_species
   use integrators, only: integrate_with_message

   implicit none

   public integrate_isobaric

   contains
      subroutine integrate_isobaric(y, t0, t_end, msg_out, m_total, n)
         integer, intent(in) :: n
         real(8), intent(inout), dimension(n) :: y
         real(8), intent(in) :: t_end, t0
         integer, intent(inout) :: m_total
         character(len=100), intent(inout), optional :: msg_out

         call integrate_with_message(dydt_isobaric, y, increment_state_isobaric, t0, t_end, msg_out, m_total)
      end subroutine integrate_isobaric

      pure function dydt_isobaric(t, y) result(dydt)
         use microphysics_common, only: cp_mixture

         real(8), intent(in) :: t
         real(8), intent(in), dimension(:) :: y
         real(8) :: c_m, dydt(size(y))

         c_m = cp_mixture(y)

         dydt = dydt_mphys(t, y, c_m)
      end function dydt_isobaric

      !> Return a new state changed by increment `dy`, method is intended to
      !> guarantee self-consistency of new state `y_new`
      !! @TODO: Implement mass-scaling
      pure function increment_state_isobaric(y, dy) result(y_new)
         real(8), intent(in), dimension(:) :: y, dy
         real(8), dimension(size(y)) :: y_new

         y_new = y + dy
      end function increment_state_isobaric
end module


module microphysics_integration
   use microphysics_constants, only: abs_tol => integration_abs_tol, rel_tol => integration_rel_tol
   use microphysics_register, only: n_variables

#ifdef MPI
   use mpi
#endif

   implicit none

   !> Will be assigned to one of the "integration helpers" which gaurantee
   !> either isometric or isobaric integration
   procedure(), pointer :: integrate_with_constraint => null()

   contains
      !> Public subroutine that will be called by ATHAM/python-wrapper etc.
      subroutine integrate(y, t0, t_end, msg_out)
         real(8), intent(inout), dimension(n_variables) :: y
         real(8), intent(in) :: t_end, t0
         character(len=100), optional :: msg_out
         character(len=100) :: msg
         real(8) :: y0(size(y))
         integer :: mpi_rank, ierror
         integer :: m_total

         ! Copy the initial state for future reference
         y0(:) = y

         msg = " "
         m_total = 0

         call integrate_with_constraint(y, t0, t_end, msg, m_total, n_variables)

         if (present(msg_out)) then
            !TODO: when calling from ATHAM this "optional" value is set although
            !it it shouldn't be. Can't use it right now
            !msg_out = msg
         else
            if (msg(1:1) /= " ") then
#ifdef MPI
               call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierror)
#endif
               print *, "==============================="
               print *, "integration failed"
               print *, mpi_rank, ":", y0
               print *, msg
               stop(0)
            endif
         endif
      end subroutine integrate
end module microphysics_integration
