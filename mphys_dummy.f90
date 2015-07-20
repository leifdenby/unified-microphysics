module mphys_dummy
    implicit none
    public n_moments, init

    ! One-moment schemes k = 1, mass
    ! Two moments, k = 0, 1 number concentration and mass
    ! Three moments, k = 0, 1, 2 number concentration, mass and radar reflectivity

    integer :: n_moments(2) = [1, 1]


    contains
    subroutine init(hydrometeor_names)
        character(16), dimension(size(n_moments)), intent(out) :: hydrometeor_names
        print *, "Dummy microphysics init called"

        hydrometeor_names(1)='cloud_water'
        hydrometeor_names(2)="rain"
    end subroutine

    subroutine d_hydrometeors(q, dt)
        real, intent(in) :: q, dt
    end subroutine
end module mphys_dummy
