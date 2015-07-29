program main
    use microphysics_initialisation, only: init
    use microphysics_register, only: cv_gases, n_gases, n_moments__max, q_flux_function
    use microphysics_common, only: cv_mixture, cp_mixture

    real, dimension(:,:), allocatable :: q

    print *, "Testing the unified microphysics package"

    call init()

    !allocate(q(n_species, n_moments__max))

    !cv_all(:) = 0.2
    !q(:,1) = 0.5

    !print *, cv_mixture(q)
    !print *, cp_mixture(q)
end
