module microphysics_register
    implicit none
    public init, ntgas, mp_ntrac

    integer :: ntgas = 1
    integer :: mp_ntrac = -1  ! will be equal to the number of hydrometeors requested by the specific implementation, and will be set on init

    character(16), dimension(:), allocatable :: hydrometeor_names(:)
    character(16), dimension(4) :: valid_hydrometeor_names
contains
    subroutine init()
        character(len=32) :: configuration = ''
        namelist /microphysics_config/ configuration
        
        ! This list below is the list of the variables that the unified
        ! microphysics currently supports
        valid_hydrometeor_names(1) = 'cloud_water'
        valid_hydrometeor_names(2) = 'rain'
        valid_hydrometeor_names(3) = 'cloud_ice'
        valid_hydrometeor_names(4) = 'graupel'

        open(1, file='{microphysics_configuration_namelist_filename}',form='formatted',status='old')
        read(1, nml=microphysics_config)
        close(1)

        if (configuration == '') then
            print *, "Please choose a microphysics model to use. The available models are:"
            print *, "  {mphys_names_joined}"
        !{mphys_methods_switch_statements}
        else
            print *, "Error the chosen microphysics implementation wasn't found"
            call exit(-1)
        endif

        call verify_hydrometeor_names()

    end subroutine init

    function hydrometeor_index(hydrometeor_name)
        integer :: hydrometeor_index
        character(32), intent(in) :: hydrometeor_name
        integer :: i = -1

        do i=1,mp_ntrac
            if (index(hydrometeor_names(i), hydrometeor_name) /= 0) then
                hydrometeor_index = i
                exit  ! break statement
            endif
        enddo
    end function hydrometeor_index

    function is_hydrometeor(i, hydrometeor_name)
        logical is_hydrometeor
        integer, intent(in) :: i
        character(32), intent(in) :: hydrometeor_name

        is_hydrometeor = index(hydrometeor_names(i), hydrometeor_name) /= 0
    end function is_hydrometeor

    subroutine allocate_local()
        allocate(hydrometeor_names(mp_ntrac))
    end subroutine allocate_local

    subroutine verify_hydrometeor_names()
        logical :: found_variable = .FALSE.
        integer :: i = -1, j = -1

        ! check that all the variables requested by the microphysics
        ! implementation are valid
        do i=1,mp_ntrac
            found_variable = .FALSE.
            do j=1,size(valid_hydrometeor_names)
                if (trim(hydrometeor_names(i)) == trim(valid_hydrometeor_names(j))) then
                    found_variable = .TRUE.
                endif
            enddo

            if (found_variable .eqv. .FALSE.) then
                print *, "Error: the following hydrometeor variable is not currently implemented in the unified microphysics:"
                print *, hydrometeor_names(i)
                call exit(-1)
            endif
        enddo
    end subroutine verify_hydrometeor_names

    !{microphysics_module_blocks}


    subroutine hydrometeor_fluxes(q, dt, dq)
        real, dimension(mp_ntrac), intent(in) :: q
        real, intent(in) :: dt
        real, dimension(mp_ntrac), intent(out) :: dq

    end subroutine hydrometeor_fluxes
end module microphysics_register
