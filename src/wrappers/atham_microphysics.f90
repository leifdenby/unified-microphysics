!> ATHAM-specific wrapper for the unified-microphysics library
!>
!> TODO:
!> - write ATHAM allocation routines
!> - write meta information routines that set the correct output variables for the netCDF output.
!> - create flag for indicating whether "dry air" specific concentration should be stored in ATHAM's state variables
!>
!> NOTES:
!>
!> - ATHAM internally stores specific concentrations in three arrays:
!>     tgasnew: gaseous tracers
!>     tracnew: incompressible tracers (solids and fluids)
!>     tpasnew: passive tracers


module microphysics
   use precision, only: kreal, kint

   use microphysics_initialisation, only: unified_microphysics_init => init

   implicit none
   private

   public init, enable_output, assign_atham_pointers, evolve_microphysics

   integer :: num_output_files = 0

contains
   subroutine init(microphysics_implementation, ntgas_mp, ntrac_mp)
      use microphysics_register, only: n_compressible_species, n_incompressible_species

      character(len=*), intent(in) :: microphysics_implementation
      integer, intent(out) :: ntgas_mp, ntrac_mp

      ! Steps:
      ! 
      ! 1. Check config file to make sure that microphysics should be used at all
      ! 2. Call init on unified-microphysics module
      ! 3. Setup the right number of tracers. In doing this check that the tracers have not already been allocated
      ! 4. Setup the specific heat capacities

      call unified_microphysics_init(microphysics_implementation, 'isometric')

      ntgas_mp = n_compressible_species
      ntrac_mp = n_incompressible_species

   end subroutine init

   ! ATHAM variables:
   !
   ! `wetnew`: water vapour
   ! `watcnew`: cloud water
   ! `watpnew`: rain water
   ! `icenew`: cloud ice
   ! `granew`: graupel

   subroutine assign_atham_pointers()
      use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain, idx_cice, idx_graupel
      use microphysics_register, only: cp_species, cv_species

      use process_data, only: wetnew, watcnew, watpnew, icenew, granew, &
                              wetflx, watcflx, watpflx, iceflx, graflx
                              !wwatc, wwatp, wice, wgra,                 &
                              !radwatc, radice, rhowat, rhoice, rhogra,  &
                              !cpwat, cpice

      use atham_module, only: tgasnew, tracnew, cptgas, cvtgas, cptrac
      use atham_module, only: tgasflx, tracflx

      integer :: ntgas, ntrac

      ntgas = 0
      ntrac = 0

      print *, "Setting up pointers"

      if (idx_water_vapour /= 0) then
         ntgas = ntgas + 1
         wetnew  => tgasnew(:,:,:,ntgas)
         wetflx  => tgasflx(:,:,:,ntgas)
         print *, "wetnew", ntgas
         cvtgas(ntgas) = cv_species(idx_water_vapour)
         cptgas(ntgas) = cp_species(idx_water_vapour)
      endif

      if (idx_cwater /= 0) then
         print *, "cloud water"
         ntrac = ntrac + 1
         watcnew  => tracnew(:,:,:,ntrac)
         watcflx  => tracflx(:,:,:,ntrac)
         cptrac(ntrac) = cp_species(idx_cwater)
      endif

      if (idx_rain /= 0) then
         print *, "rain"
         ntrac = ntrac + 1
         watpnew  => tracnew(:,:,:,ntrac)
         watpflx  => tracflx(:,:,:,ntrac)
         cptrac(ntrac) = cp_species(idx_rain)
      endif

      if (idx_cice /= 0) then
         print *, "cloud ice"
         ntrac = ntrac + 1
         icenew  => tracnew(:,:,:,ntrac)
         iceflx  => tracflx(:,:,:,ntrac)
         cptrac(ntrac) = cp_species(idx_cice)
      endif

      if (idx_graupel /= 0) then
         print *, "graupel", idx_graupel
         ntrac = ntrac + 1
         granew  => tracnew(:,:,:,ntrac)
         graflx  => tracflx(:,:,:,ntrac)
         cptrac(ntrac) = cp_species(idx_graupel)
      endif

      print *, ntgas, ntrac

      print *, "Done"

   end subroutine assign_atham_pointers

   !> Subroutine for enabling output being written for a particular variable
   !>
   !> In the definitions below:
   !>  -  `var_type` is the variable type on of `gas`, `rac` and `pas` which are
   !>     respectivily compressible (gassous tracers), imcompressible tracers
   !>     (solids or fluids) and passive tracers.
   !>  -  `output_type` is the type of output, either at "movie file" or a "picture file"
   !>
   !> `nt{var_type}_{output_type}` defines the number of outputs of `var_type` that are `output_type` format
   !>
   !> `it{var_type}_{output_type}(n) = m` defines the mapping from defined var
   !>     number `n`, to the indexing into array `t{var_type}new` to get the correct tracer variable
   !>
   !> `var_t{var_type}_{output_type}(n)` defines the variable name in the outputfile for output number `n` of `var_type` in `output_type`
   !>
   !> `des_t{var_type}_{output_type}(n)` is a text description for output number `n` of `var_type` in `output_type`
   !>
   !> @todo Enable output into movie files (by adding an optional argument)
   !> @todo Implement better variable descriptions
   subroutine enable_output(var_name, var_name_output)
      use atham_module,   only: itgas_picture, ntgas_picture, var_tgas_picture, des_tgas_picture
      use atham_module,   only: itrac_picture, ntrac_picture, var_trac_picture, des_trac_picture
      use atham_module,   only: itpas_picture, ntpas_picture, var_tpas_picture, des_tpas_picture

      use atham_module,   only: itgas_movie, ntgas_movie, var_tgas_movie, des_tgas_movie
      use atham_module,   only: itrac_movie, ntrac_movie, var_trac_movie, des_trac_movie
      use atham_module,   only: itpas_movie, ntpas_movie, var_tpas_movie, des_tpas_movie

      use microphysics_register, only: idx_cwater, idx_water_vapour
      use microphysics_register, only: idx_rain, idx_cice, idx_graupel

      character(len=*), intent(in) :: var_name
      character(len=*), intent(in) :: var_name_output
      
      character(len=32) :: var_description
      integer :: picture_output_idx = -1
      integer :: movie_output_idx = -1
      integer :: species_idx = -1
      integer :: var_type = -1  ! 0: compressible, 1: incompressible, 2: passive

      var_description = var_name
      movie_output_idx = -1

      ! TODO: This needs fixing, the indexing is hardcoded, needs to be stored
      ! somewhere

      if (trim(var_name) == 'water_vapour') then
         ntgas_picture = ntgas_picture + 1
         picture_output_idx = ntgas_picture
         ntgas_movie = ntgas_movie + 1
         movie_output_idx = ntgas_movie

         species_idx = 1 ! idx_water_vapour-1
         var_type = 0
      else if (trim(var_name) == 'cloud_water') then
         ntrac_picture = ntrac_picture + 1
         picture_output_idx = ntrac_picture
         ntrac_movie = ntrac_movie + 1
         movie_output_idx = ntrac_movie

         species_idx = 1 ! idx_cwater
         var_type = 1
      else if (trim(var_name) == 'rain') then
         ntrac_picture = ntrac_picture + 1
         picture_output_idx = ntrac_picture
         ntrac_movie = ntrac_movie + 1
         movie_output_idx = ntrac_movie

         species_idx = 2 ! idx_rain
         var_type = 1
      else if (trim(var_name) == 'cloud_ice') then
         ntrac_picture = ntrac_picture + 1
         picture_output_idx = ntrac_picture
         ntrac_movie = ntrac_movie + 1
         movie_output_idx = ntrac_movie

         species_idx = 3 ! idx_cice
         var_type = 1
      else if (trim(var_name) == 'graupel') then
         ntrac_picture = ntrac_picture + 1
         picture_output_idx = ntrac_picture
         ntrac_movie = ntrac_movie + 1
         movie_output_idx = ntrac_movie

         species_idx = 4 ! idx_graupel
         var_type = 1
      else
         print *, "Error: output for variable ", trim(var_name), " no implemented"
         call exit(-2)
      endif

      if (var_type == 0) then
         itgas_picture(picture_output_idx) = species_idx
         var_tgas_picture(picture_output_idx) = var_name_output
         des_tgas_picture(picture_output_idx) = var_description

         if (movie_output_idx /= -1) then
            itgas_movie(movie_output_idx) = species_idx
            var_tgas_movie(movie_output_idx) = var_name_output
            des_tgas_movie(movie_output_idx) = var_description
         endif

      else if (var_type == 1) then
         itrac_picture(picture_output_idx) = species_idx
         var_trac_picture(picture_output_idx) = var_name_output
         des_trac_picture(picture_output_idx) = var_description

         if (movie_output_idx /= -1) then
         itrac_movie(movie_output_idx) = species_idx
         var_trac_movie(movie_output_idx) = var_name_output
         des_trac_movie(movie_output_idx) = var_description
         endif
      else
         print *, "Error: output for variable ", trim(var_name), " no implemented"
         call exit(-2)
      endif
   end subroutine enable_output

   !> Updates microphysics-related state variables in ATHAM using the selected
   !> microphysics implementation.
   !
   ! Steps:
   ! 1. Create representation of state in variables that the microphysics
   !    routines understand
   ! 2. Call microphysics routines
   ! 3. Do mass-fix
   ! 4. Update condensate fluxes
   ! 5. Update temperature using defined heat capacities
   subroutine evolve_microphysics(t, dt)
      use atham_module, only: nx, ny, nz, nxl, nxr, nyv, nyh
      use atham_module,   only: x, y, z, xv, yv

      use process_data, only: wetnew, watcnew, watpnew, icenew, granew
      use process_data, only: wetflx, watcflx, watpflx, iceflx, graflx

      use atham_module, only: ntgas, ntrac, tgasnew, tracnew
      use atham_module, only: pnew, tempnew, p0, tetnew
      use atham_module, only: iflgs
      use atham_module, only: tetflx, pflx

      use microphysics_register, only: n_variables, q_species_flag
      use microphysics_register, only: idx_water_vapour, idx_temp, idx_pressure
      use microphysics_register, only: idx_graupel, idx_cwater, idx_rain, idx_cice
      use microphysics_common, only: cp_gas, cv_gas
      use microphysics_integration, only: integrate_microphysics => integrate
      !use mphys_kessler_old, only: integrate_microphysics => integrate

      ! TODO: Remove these once we know how to calculate d_theta
      use atham_module,   only: ntgas, ntrac, tgasnew, tracnew, cptgas

      !> Rg and cp_g from end of last timestep
      use atham_module,   only: rgasnew, cpgasnew
      use atham_module,   only: cptot
      use phys_constants, only: cpair, ps0, pi
      use phys_constants, only: r0, r1, r1h, epsmin

      real(kreal), intent(in) :: t, dt

      real(kreal), dimension(n_variables) :: y_, dydt, dy, y_2
      real(kreal) :: p = 0.0, d_theta=0.0

      integer :: i=-1, j=-1, k=-1, jh=-1, jv=-1, ir=-1, il=-1

      real(kreal) :: temp_n, pres_n, tet_n, cp_g, R_g

      ! Each microphysics routine serves one purpose:
      ! 
      ! Given the current environmental state which may be super/sub-saturated
      ! calculate a termodynamic state consistent with the present water phases
      !
      ! Given:
      !    gasnew, drynew, cp, cv, tempnew
      !    wetnew,watcnew,watpnew,icenew,granew, 
      ! Compute:
      !    wetflx,watcflx,watpflx,iceflx,graflx

      do k=2,nz
         do j=nyv,nyh
            jh=min(nyh,j+1)
            jv=j-1
            do i=nxl,nxr
               ir=min(nxr,i+1)
               il=i-1
               if (iflgs(i,j,k)==1) then
                  p = p0(k) + pnew(i,j,k)

                  y_(:) = 0.0
                  y_(idx_temp) = tempnew(i,j,k)
                  y_(idx_pressure) = p
                  y_(idx_water_vapour) = wetnew(i,j,k)
                  y_(idx_cwater) = watcnew(i,j,k)

                  if (idx_cice /= 0) then
                     y_(idx_cice) = icenew(i,j,k)
                  endif
                  if (idx_graupel /= 0) then
                     y_(idx_graupel) = granew(i,j,k)
                  endif

                  !print *, ""
                  !print *, '(x, z)=', x(i), z(k)
                  !print *, "--> p T q_v, q_c, q_r", p, y_(idx_temp), y_(idx_water_vapour), y_(idx_cwater), y_(idx_rain)

                  y_2(:) = y_(:)
                  call integrate_microphysics(y_, t, t+dt)  ! modifies y_
                  dy = y_ - y_2

                  ! calculate a new potential temperature which is consistent
                  ! with the updated pressure and temperature
                  pres_n = y_(idx_pressure)
                  temp_n = y_(idx_temp)
                  R_g = cp_gas(y_) - cv_gas(y_)
                  cp_g = cp_gas(y_)
                  tet_n = temp_n*(ps0/pres_n)**(R_g/cp_g)

                  tetflx(i,j,k) = tetflx(i,j,k) + (tet_n - tetnew(i,j,k))

                  pflx(i,j,k) = pflx(i,j,k) + (pres_n - p)

                  wetflx(i,j,k)  = wetflx(i,j,k)  + dy(idx_water_vapour)
                  watcflx(i,j,k) = watcflx(i,j,k) + dy(idx_cwater)
                  watpflx(i,j,k) = watpflx(i,j,k) + dy(idx_rain)
                  if (idx_cice /= 0) then
                     iceflx(i,j,k)  = iceflx(i,j,k)  + dy(idx_cice)
                  endif
                  if (idx_graupel /= 0) then
                     graflx(i,j,k)  = graflx(i,j,k)  + dy(idx_graupel)
                  endif

               endif
            enddo
         enddo
      enddo
   end subroutine evolve_microphysics
end module microphysics
