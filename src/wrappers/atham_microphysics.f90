module atham_microphysics
   use microphysics_register, only: mp_ntrac, is_hydrometeor, hydrometeor_names

   use atham_module,   only : nx, ny, nz

   logical :: init_completed = .FALSE.

contains
    subroutine hydrometeor_fluxes(q, dt)
        real, dimension(mp_ntrac), intent(in) :: q
        real, intent(in) :: dt
        real, dimension(mp_ntrac), intent(out) :: dq

         do j=nyv,nyh
            jh=min(nyh,j+1)
            jv=j-1
            do i=nxl,nxr
               ir=min(nxr,i+1)
               il=i-1
               if (iflgs(i,j,k)==1) then
               endif
            enddo
         enddo

    end subroutine hydrometeor_fluxes

   subroutine allocate_tracers()
      use atham_module,   only: ntgas, ntrac, ntpas
      use atham_module,   only: cptgas, cvtgas
      use atham_module,   only: radtrac, rhotrac, cptrac

      integer :: i = -1
      
      ! TODO: introduce check to make sure we're not meant to be allocating
      ! spacing for multi-moment formulations, I don't know how to implement
      ! that yet

      ! TODO: Make sure we correctly offset the number of indecies needed if
      ! we have non-microphysics related tracers

      ntgas=1
      ntrac=mp_trac
      ntpas=4

      !itgas_wet=1
      !itrac_hydro(1)=1
      !itrac_hydro(2)=2
      !itrac_hydro(3)=3
      !itrac_hydro(4)=4


      !alpha_tgas(1)=alpha_tet
      !alpha_trac(1:4)=r1
      !alpha_tpas(1:4)=alpha_tet
      !gas_tpas(1:4)=.true.

      ! Set up tracer properties

      ! TODO: These assignments should eventually include setting the tracer
      ! properties in the unified microphysics module
      cptgas(1)=cphum
      cvtgas(1)=cvhum

      do i=0,mp_trac
         if (is_hydrometeor(i, 'cloud_water')) then
            cptrac(i)=4183._kreal
            rhotrac(i)=1000._kreal
            radtrac(i)=0.00001_kreal
         else if (is_hydrometeor(i, 'rain')) then
            cptrac(i)=4183._kreal
            rhotrac(i)=1000._kreal
            radtrac(i)=0.00001_kreal
         else if (is_hydrometeor(i, 'cloud_ice')) then
            cptrac(i)=2103._kreal
            rhotrac(i)= 917._kreal
            radtrac(i)=0.00001_kreal
         else if (is_hydrometeor(i, 'graupel')) then
            cptrac(i)=2103._kreal
            rhotrac(i)= 700._kreal
            radtrac(i)=0.00001_kreal
         else
            print *, "Error: hydrometeor", hydrometeor_names(i), "is not known to ATHAM"
            call exit(01)
         endif
         
      enddo
   end subroutine

   !subroutine assign_pointers()
     
     !use process_data, only: wetnew, watcnew, watpnew, icenew, granew, &
                             !wetflx, watcflx, watpflx, iceflx, graflx, &
                             !wwatc, wwatp, wice, wgra,                 &
                             !radwatc, radice, rhowat, rhoice, rhogra,  &
                             !cpwat, cpice

     !wetnew  => tgasnew(:,:,:,itgas_wet)

     !watcnew => tracnew(:,:,:,itrac_hydro(1))
     !watpnew => tracnew(:,:,:,itrac_hydro(2))
     !icenew  => tracnew(:,:,:,itrac_hydro(3))
     !granew  => tracnew(:,:,:,itrac_hydro(4))
   
     !wetflx  => tgasflx(:,:,:,itgas_wet)
     !watcflx => tracflx(:,:,:,itrac_hydro(1))
     !watpflx => tracflx(:,:,:,itrac_hydro(2))
     !iceflx  => tracflx(:,:,:,itrac_hydro(3))
     !graflx  => tracflx(:,:,:,itrac_hydro(4))
   
     !wwatc   => wtrac(:,:,:,itrac_hydro(1))
     !wwatp   => wtrac(:,:,:,itrac_hydro(2))
     !wice    => wtrac(:,:,:,itrac_hydro(3))
     !wgra    => wtrac(:,:,:,itrac_hydro(4))
     
     !radwatc => radtrac(itrac_hydro(1))
     !radice  => radtrac(itrac_hydro(3))
     
     !rhowat  => rhotrac(itrac_hydro(1))
     !rhoice  => rhotrac(itrac_hydro(3))
     !rhogra  => rhotrac(itrac_hydro(4))

     !cpwat   => cptrac(itrac_hydro(1))
     !cpice   => cptrac(itrac_hydro(3))
 !end subroutine

end module atham_microphysics
