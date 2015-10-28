program main
    use microphysics_initialisation, only: init

    print *, "Testing the unified microphysics package"

    call init()

    !allocate(q(n_species, n_moments__max))

    !cv_all(:) = 0.2
    !q(:,1) = 0.5

    !print *, cv_mixture(q)
    !print *, cp_mixture(q)
end
