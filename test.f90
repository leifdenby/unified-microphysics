program test
    use microphysics_register, only: init

    print *, "Testing the unified microphysics package"

    call init()
end program test
