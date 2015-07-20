!> Wrapper for old microphysics implementation in ATHAM v5.0
!>
!
! -> kessler_init(do_diag, kessler_diag_in): init call
! 1. setup output file for diagnostics
!
! -> kessler(dt): call to compute fluxes from host model
! 1. compute air density
! 2. compute contribution from each process (p1-p19)
! 3. compute fluxes
!
!
!----------------------------------------------------------------!
! considered cloud microphysical processes:                      !
!                                                                !
! p1      accretion of cloud droplets by rain                    !
! p2,p3   accretion of ice particles by rain                     !
! p4      accretion of ice particles by graupel                  !
! p5      accretion of cloud droplets by graupel                 !
! p6      accretion of rain by graupel                           !
!                                                                !

! p7      melting of graupel to rain                             !
! p8      melting ice to cloud water                             !
! p9      statistical freezing of cloud to ice                   !
! p10     statistical freezing of rain to graupel                !
!                                                                !

! p11     autoconversion of cloud droplets to rain               !
! p12     autoconversion of icepartiles in graupel               !
!                                                                !

! p13     cond/evap on cloud droplets                            !
! p14     cond/evap on rain droplets                             !
! p15     cond/evap on melting graupel                           !
!                                                                !

! p16     subl/dep on ice                                        !
! p17     subl/dep on graupel                                    !
!                                                                !

! p18     new formation of cloud droplets                        !
! p19     new formation of ice particles                         !
!----------------------------------------------------------------!

module mphys_kessler_old

    implicit none
    public n_moments, init

    ! One-moment schemes k = 1, mass
    ! Two moments, k = 0, 1 number concentration and mass
    ! Three moments, k = 0, 1, 2 number concentration, mass and radar reflectivity

    integer :: n_moments(2) = [1, 1]

  !> `precision` module imports inlined below
  ! use precision, only: kreal, kint
  integer, public, parameter :: real2 = selected_real_kind(3),            &
                                real4 = selected_real_kind(6),            &
                                real8 = selected_real_kind(12),           &
                                kreal = real8
  integer, public, parameter :: int4 = selected_int_kind(7),              &
                                int8 = selected_int_kind(13),             &
                                kint = int4
   !!

   !> constants

   real, parameter :: r1 = 1.0
   real, parameter :: r1q = 0.25
   real, parameter :: r4 = 4.0
   real, parameter :: r3 = 3.0
   real, parameter :: pi = 3.14
   real, parameter :: xpi =4.0/3.0*3.14
   real, parameter :: epsmach = 1.0e-14

   ! density of water [kg/m3]
   real, parameter :: rhowat = 1.0e3

   ! radius of water droplets
   ! radwatc: Cloud particle radius (radius - water - cloud)
   real, parameter :: radwatc = 1.0
   ! radwatp: Hydrometeor radius (radius - water - precip)
   ! rhoice: Density of ice and ice particles
   ! radice: Radius of ice particles (not, I think, ice hydrometeors which are graupel)

   real, parameter :: heatlat = 1.0

   real, parameter :: cphum = 1.0
   real, parameter :: cvhum = 3.0
   real, parameter :: gasf=cphum-cvhum
   real, parameter :: T0 = 273.15
   real, parameter :: ps0 = 101325.0
   real, parameter :: ronull = 1.225

    contains
    subroutine init(hydrometeor_names)
        character(16), dimension(size(n_moments)), intent(out) :: hydrometeor_names
        print *, "Dummy microphysics init called"

        hydrometeor_names(1)='cloud_water'
        hydrometeor_names(2)="rain"
    end subroutine

    subroutine d_hydrometeors(q, dt, dq)
         !use phys_constants, only: r0, r1, r2, r1h, r1q, epsmin, epsmach
         !use phys_constants, only: cphum, cvhum, cpair, ps0
         !use phys_constants, only: heatlat, heatfrz, heatsub
         !use process_data,   only: wetnew, watcnew, watpnew, icenew, granew, &
                                   !wwatp, wgra, radwatc, radice,             &
                                   !rhowat, rhoice, cpwat
         !use atham_module,   only: ntgas, ntrac, tgasnew, tracnew
         !use atham_module,   only: cptgas, cptrac, cptot
         !use atham_module,   only: tss, ronull, viskos
         !use atham_module,   only: iflgs

        real, dimension(size(n_moments)), intent(in) :: q
        real, intent(in) :: dt
        real, dimension(size(n_moments)), intent(out) :: dq


      contains
         function thermal_conductivity(temp)
            real :: thermal_conductivity
            real, intent(in) :: temp

            real, parameter :: adiffk=2.4e-2_kreal
            real, parameter :: bdiffk=8.e-5_kreal

            thermal_conductivity=adiffk+bdiffk*temp    
         end function thermal_conductivity

         function water_vapour_diffusivity(temp, pressure)
            real :: water_vapour_diffusivity
            real, intent(in) :: temp, pressure

            real, parameter :: adiffd=2.11e-5_kreal
            real, parameter :: bdiffd=1.94_kreal

            water_vapour_diffusivity=adiffd*(temp/T0)**bdiffd*ps0/(pressure)
         end function water_vapour_diffusivity

         function saturation_vapour_pressure(temp)
            real :: saturation_vapour_pressure
            real, intent(in) :: temp

            real, parameter :: awat=17.25_kreal
            real, parameter :: bwat=36._kreal
            real, parameter :: t0=273.15_kreal
            real, parameter :: asat=610.7_kreal

            real :: expon2, satpw

            expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
            saturation_vapour_pressure=asat*exp(expon2)

         end function saturation_vapour_pressure

         subroutine cond_evap(watcnew, temp, pressure, rgasnew, gasnew, wetnew, densair, watpnew, viskos, granew)
         ! watcnew: q_c (cloud water specific concentration)
         ! rgasnew: Rg (specific gas constant for gas mixture)
         ! gasnew: q_g (specifc concentration of gas mixture)
         ! wetnew: q_ ??
         ! densair: density of mixture ??
         ! watpnew: q_r (specific concentration of rain water)
         ! viskos: ??
         ! granew: q_graupel (specific concentration of graupel)

         ! p13     cond/evap on cloud droplets                            !
         ! p14     cond/evap on rain droplets                             !
         ! p15     cond/evap on melting graupel                           !

         real, intent(in) :: watcnew, temp, pressure, rgasnew, gasnew, wetnew, densair, watpnew, viskos, granew
         real :: p13, p14, p15

         real, parameter :: pexp=2.75_kreal

         real, parameter :: alamp=421._kreal
         real, parameter :: ap14=6.28e7_kreal, bp14=2.65e8_kreal

         real, parameter :: alamg=71.8_kreal
         real, parameter :: ap15=1.22e5_kreal, bp15=4.11e5_kreal

         real :: xnwatc 
         real :: evap
         real :: diffk
         real :: diffd
         real :: satpw
         real :: pfak
         real :: satww
         real :: supsatw

         real :: conv
         real :: xlamp
         real :: densfac

         real :: xlamg



         !----------------------------------------------------------!
         ! super-/sub-saturation quantities of water vapor          !
         !----------------------------------------------------------!
         diffk=thermal_conductivity(temp)
         diffd=water_vapour_diffusivity(temp, pressure)
         satpw=saturation_vapour_pressure(temp)

         pfak=rgasnew/(gasf*pressure)*gasnew
         satww=satpw*pfak

         supsatw=wetnew/satww-r1

         evap=supsatw/((heatlat/(gasf*temp)-r1)*heatlat/(diffk*temp) &
               +gasf*temp/(diffd*satpw))


         !----------------------------------------------------------!
         ! condensation and evaporation on droplets                 ! 
         !----------------------------------------------------------!
         ! cond/evap on cloud                                       !
         !----------------------------------------------------------!

         p13=r4*pi*xnwatc*radwatc*evap

         !p13=sign(r1,p13)*min(abs(p13*dt),xswatc/xsall*dwetw)
         !p13=max(p13,-watcnew(i,j,k))
         !----------------------------------------------------------!
         ! cond/evap on rain                                        !
         !----------------------------------------------------------!
         conv=gasnew/densair

         densfac=sqrt(ronull/densair)
         xnwatc=watcnew/(rhowat*xpi*radwatc**3)
         xlamp=(watpnew/conv)**r1q/alamp

         p14=conv*evap*(ap14*xlamp*xlamp                        &
         +bp14*xlamp**pexp*sqrt(densfac/viskos))

         !p14=sign(r1,p14)*min(abs(p14*dt),xswatp/xsall*dwetw)
         !p14=max(p14,-watpnew(i,j,k))
         !----------------------------------------------------------!
         ! cond/evap on melting graupel                             !
         !----------------------------------------------------------!
         xlamg=(granew/conv)**r1q/alamg

         p15=conv*evap*(ap15*xlamg*xlamg                        &
         +bp15*xlamg**pexp*sqrt(densfac/viskos))

         !p15=sign(r1,p15)                                  &
         !*min(abs(p15*dt),xsgra/xsall*dwetw)*itemp
         !p15=max(p15,-granew(i,j,k))

         end subroutine cond_evap

         !subroutine particle_droplet_formation
         !!----------------------------------------------------------!
         !! new formation of cloud droplets                          !
         !!----------------------------------------------------------!
         !! remaining supersaturation water vapor forms iceparitcles !
         !!----------------------------------------------------------!
         !xnice = icenew(i,j,k)/(rhoice*xpi*radice**3)
         !p19(i,j)=max(r0,ap18*conv*exp(bp18*max(r0,t0-temp))-xnice)
         !p19(i,j)=(1-itemp)*min(p19(i,j),dwete)
         !!----------------------------------------------------------!
         !! remaining supersaturation water vapor forms clouddroplets!
         !!----------------------------------------------------------!
         !p18(i,j)=max(r0,dwetw-p19(i,j))
         !end subroutine particle_droplet_formation

         !subroutine sublimation_deposition
         !!----------------------------------------------------------!
         !! sublimation and deposition on particels                  !
         !!----------------------------------------------------------!
         !! subl/dep on ice                                          !
         !!----------------------------------------------------------!
         !p16(i,j)=ap16*pi*xnice*radice*subl
         !p16(i,j)=sign(r1,p16(i,j))                                  &
         !*min(abs(p16(i,j)*dt),xsice*dwete/xssol)*(1-itemp)
         !p16(i,j)=max(p16(i,j),-icenew(i,j,k))
         !!----------------------------------------------------------!
         !! subl/dep on graupel                                      !
         !!----------------------------------------------------------!
         !p17(i,j)=conv*subl*(ap17*xlamg*xlamg                        &
         !+bp17*xlamg**pexp*sqrt(densfac/viskos(k)))
         !p17(i,j)=sign(r1,p17(i,j))                                  &
         !*min(abs(p17(i,j)*dt),xsgra*dwete/xssol)*(1-itemp)
         !p17(i,j)=max(p17(i,j),-granew(i,j,k))
         !end subroutine sublimation_deposition

         !subroutine melting_freezing
         !! p7      melting of graupel to rain                             !
         !! p8      melting ice to cloud water                             !
         !! p9      statistical freezing of cloud to ice                   !
         !! p10     statistical freezing of rain to graupel                !
         !!                                                                !
         !!----------------------------------------------------------!
         !! melting and freezing                                     !
         !!----------------------------------------------------------!
         !! melting of graupel to rain                               !
         !!----------------------------------------------------------!
         !p7(i,j)=-p15(i,j)*heatlat/heatfrz                               &
         !+ap7*conv/heatfrz*diffk*(temp-t0)                          &
         !*(bp7*xlamg*xlamg+cp7*xlamg**pexp*sqrt(densfac/viskos(k))) &
         !+cpwat/heatfrz*(temp-t0)*(p5(i,j)+p6(i,j))
         !p7(i,j)=max(r0,min(p7(i,j)*dt,granew(i,j,k),r1h*dampmel)*itemp)
         !p6(i,j)=p6(i,j)*(1-itemp)
         !dampmel=dampmel-p7(i,j)
         !!----------------------------------------------------------!
         !! melting total ice to cloudwater  at t0                   !
         !!----------------------------------------------------------!
         !p8(i,j)=max(r0,min(dampmel,max(icenew(i,j,k)-p4(i,j),r0))*itemp)
         !!----------------------------------------------------------!
         !! statistical freezing of rain to graupel                  !
         !!----------------------------------------------------------!
         !p10(i,j)=conv*ap10*xlamp**7*(exp(bp10*max(r0,t0-temp))-r1)
         !p10(i,j)=max(r0,min(p10(i,j)*dt,                            &
         !watpnew(i,j,k)-p3(i,j)-p6(i,j),r1h*dampmel))
         !dampmel=dampmel-p10(i,j)
         !!----------------------------------------------------------!
         !! statistical freezing of cloud to ice                     !
         !!----------------------------------------------------------!
         !p9(i,j)=ap9*xnwatc*radwatc**6*(exp(bp9*max(r0,t0-temp))-r1)
         !p9(i,j)=max(r0,min(p9(i,j)*dt,watcnew(i,j,k)-p2(i,j),dampmel))
         !end subroutine melting_freezing

         !subroutine autoconversion
         !! p11     autoconversion of cloud droplets to rain               !
         !! p12     autoconversion of icepartiles in graupel               !

         !!----------------------------------------------------------!
         !! autoconversion                                           !
         !!----------------------------------------------------------!
         !! autoconversion of clouddroplets to rain (kessler)        !
         !!----------------------------------------------------------!
         !p11(i,j)=ap11*max(r0,watcnew(i,j,k)-bp11*conv)
         !p11(i,j)=max(r0,min(p11(i,j)*dt,watcnew(i,j,k)))
         !!----------------------------------------------------------!
         !! autoconversion of icepartiles in graupel (lin et al.)    !
         !!----------------------------------------------------------!
         !p12(i,j)=ap12*exp(bp12*(temp-t0))*max(r0,icenew(i,j,k)-ap12*conv)
         !p12(i,j)=max(r0,min(p12(i,j)*dt,icenew(i,j,k))*(1-itemp))

         !end subroutine autoconversion

         !subroutine accretion
         !! p1      accretion of cloud droplets by rain                    !
         !! p2,p3   accretion of ice particles by rain                     !
         !! p4      accretion of ice particles by graupel                  !
         !! p5      accretion of cloud droplets by graupel                 !
         !! p6      accretion of rain by graupel                           !

         !!----------------------------------------------------------!
         !! accretion                                                !
         !!----------------------------------------------------------!
         !! accretion of cloud droplets by rain                      !
         !!----------------------------------------------------------!
         !zylr=min(azylr*xlamp**arexp,bzylr*xlamp**4,czylr*xlamp**5)
         !p1(i,j)=watcnew(i,j,k)*zylr*densfac
         !p1(i,j)=max(r0,min(p1(i,j)*dt,watcnew(i,j,k)))
         !!----------------------------------------------------------!
         !! accretion of ice particles by rain forming rain    at t>t0
         !!                                    forming graupel at t<t0
         !!----------------------------------------------------------!
         !p2(i,j)=icenew(i,j,k)*zylr*densfac
         !p2(i,j)=max(r0,min(p2(i,j)*dt,icenew(i,j,k))*(1-itemp))
         !!----------------------------------------------------------!
         !! change of rainconcentration at t<t0                      !
         !!----------------------------------------------------------!
         !xniceac=p2(i,j)/(rhoice*xpi*radice**3)
         !p3(i,j)=min(r1,xniceac/max(xnwatp,epsmin))*watpnew(i,j,k)
         !p3(i,j)=max(r0,min(p3(i,j),watpnew(i,j,k),r1h*dampmel))
         !dampmel=dampmel-p3(i,j)
         !!----------------------------------------------------------!
         !! accretion of ice particles by graupel                    !
         !!----------------------------------------------------------!
         !p4(i,j)=exp(0.05*(temp-t0))*icenew(i,j,k)*densfac*apacc*xlamg**arexp
         !p4(i,j)=max(r0,min(p4(i,j)*dt,icenew(i,j,k))*(1-itemp))
         !!----------------------------------------------------------!
         !! accretion of clouddroplets by graupel forming graupel at t<t0
         !!                                       forming rain    at t>t0
         !!----------------------------------------------------------!
         !p5(i,j)=watcnew(i,j,k)*densfac*apacc*xlamg**arexp
         !p5(i,j)=min(p5(i,j)*dt,watcnew(i,j,k)-p1(i,j)*itemp)
         !p5(i,j)=max(r0,p5(i,j)*itemp+min(p5(i,j),r1h*dampmel)*(1-itemp))
         !dampmel=dampmel-p5(i,j)*(1-itemp)
         !!----------------------------------------------------------!
         !! accretion of rain by graupel forming graupel t<t0        !
         !!----------------------------------------------------------!
         !p6(i,j)=1.194e15*abs(wgra(i,j,k)-wwatp(i,j,k))*conv         &
         !*(r5 *xlamp**6*xlamg                                   &
         !+r2 *xlamp**5*xlamg*xlamg                              &
         !+r1h*xlamp**4*xlamg**3)
         !p6(i,j)=min(p6(i,j)*dt,watpnew(i,j,k)-p3(i,j))
         !p6(i,j)=max(r0,p6(i,j)*itemp+min(p6(i,j),r1h*dampmel)*(1-itemp))

         !end subroutine accretion


       !-------------------------------------------------------------------
       !subroutine kessler_processes(k,densair,dt,                          &
                          !p1 ,p2 ,p3 ,p4 ,p5 ,p6 ,p7 ,p8 ,p9,              &
                          !p10,p11,p12,p13,p14,p15,p16,p17,p18,p19)

         !real(kreal), dimension(nx,ny,nz), intent(in) :: densair
         !real(kreal), intent(in) :: dt
         !real(kreal), dimension(nx,ny), intent(out) ::                     &
              !p1 ,p2 ,p3 ,p4 ,p5 ,p6 ,p7 ,p8 ,p9,                          &
              !p10,p11,p12,p13,p14,p15,p16,p17,p18,p19
         !integer(kint), intent(in) :: k
         !!----------------------------------------------------------------!
         !! local variables                                                !
         !!----------------------------------------------------------------!
         !real(kreal) :: gasf,pi,xpi,gasnew,drynew,cp,cv,expon,pfac,temp,   &
              !conv,diffk,diffd,densfac,radwatp,radgra,                     &
              !xnwatc,xnwatp,xnice,xswatc,xswatp,xsice,xsgra,xsall,xssol,   &
              !xniceac,xlamp,xlamg,gas,pfak,expon1,satpe,satwe,supsate,subl,&
              !damp,dwete,expon2,satpw,satww,supsatw,evap,dwetw,zylr,       &
              !dampmel,damptet

         !integer(kint) :: i,j,jh,jv,ir,il,iflg,itemp

         !real(kreal), parameter ::                                         &
              !r3=3._kreal,r4=4._kreal,r5=5._kreal,                         &
              !r1th=0.1_kreal,r1hd=0.01_kreal,                              &
              !arwatp=4.35e-3_kreal,argra=2.57e-2_kreal,                    &
              !alamp=421._kreal,alamg=71.8_kreal,                           &
              !adamp=0.1_kreal,bdamp=0.27_kreal,                            &
              !cdamp=1.884_kreal,ddamp=0.0054_kreal,                        &
              !edamp=0.9_kreal,fdamp=4.395_kreal,gdamp=0.015_kreal,         &
              !satemin=1.e-4_kreal,                                         &
              !ap16=6.4_kreal,                                              &
              !ap17=1.22e5_kreal,bp17=4.11e5_kreal,pexp=2.75_kreal,         &
              !ap14=6.28e7_kreal,bp14=2.65e8_kreal,                         &
              !ap15=1.22e5_kreal,bp15=4.11e5_kreal,                         &
              !azylr=3.71e9_kreal,arexp=3.5_kreal,                          &
              !bzylr=1.88e11_kreal,czylr=5.61e15_kreal,                     &
              !apacc=3.56e6_kreal,                                          &
              !ap7=7.603e4_kreal,bp7=1.6_kreal,cp7=5.36_kreal,              &
              !ap10=1.974e14_kreal,bp10=0.66_kreal,                         &
              !ap9=1.76e6_kreal,bp9=0.66_kreal,                             &
              !ap11=1.e-3_kreal,bp11=5.e-4_kreal,                           &
              !ap12=1.e-3_kreal,bp12=2.5e-2_kreal,                          &
              !ap18=1.e-14_kreal,bp18=0.6_kreal
         !real(kreal), parameter ::                                         &
              !t0=273.15_kreal,  asat=610.7_kreal,                          &
              !aice=22.33_kreal, bice=2._kreal
         !!----------------------------------------------------------------!
         !! start loop over grid points                                    !
         !!----------------------------------------------------------------!
         !gasf=cphum-cvhum
         !pi  =r4*atan(r1)
         !xpi =r4/r3*pi
                  !!----------------------------------------------------------!
                  !! thermodynamic quantities                                 !
                  !!----------------------------------------------------------!
                  !gasnew=r1-sum(tracnew(i,j,k,1:ntrac))
                  !gasnew=max(r1hd,min(gasnew,r1))
                  !drynew=gasnew-sum(tgasnew(i,j,k,1:ntgas))
                  !cp=(cpair*drynew+sum(cptgas(1:ntgas)*tgasnew(i,j,k,1:ntgas))) &
                       !/gasnew
                  !expon=rgasnew(i,j,k)/cp
                  !pfac =(ps0/(p0(k)+pnew(i,j,k)))**expon
                  
                  !cv=cp-rgasnew(i,j,k)
                  !cv=cv*gasnew+sum(cptrac(1:ntrac)*tracnew(i,j,k,1:ntrac))    &
                              !*(r1-gasnew)
!#ifndef OLD_TEMP
                  !temp=tempnew(i,j,k)
!#else
                  !iflg=r1h+sign(r1h,tempnew(i,j,k)-tss(k)/pfac+r3)
                  !temp=tempnew(i,j,k)*iflg                                    &
                    !+max(tempnew(i,j,k),                                      &
                        !(tempnew(i,j,k)+tempnew(ir,j,k)*iflgs(ir,j,k)         &
                                       !+tempnew(il,j,k)*iflgs(il,j,k))        &
                        !/(1+iflgs(ir,j,k)+iflgs(il,j,k)),                     &
                        !(tempnew(i,j,k)+tempnew(i,jh,k)*iflgs(i,jh,k)         &
                                       !+tempnew(i,jv,k)*iflgs(i,jv,k))        &
                        !/(1+iflgs(i,jh,k)+iflgs(i,jv,k)),                     &
                        !(tempnew(ir,j,k)*iflgs(ir,j,k)                        &
                        !+tempnew(il,j,k)*iflgs(il,j,k))                       &
                        !/max(1,iflgs(ir,j,k)+iflgs(il,j,k)),                  &
                        !(tempnew(i,jh,k)*iflgs(i,jh,k)                        &
                        !+tempnew(i,jv,k)*iflgs(i,jv,k))                       &
                        !/max(1,iflgs(i,jh,k)+iflgs(i,jv,k)))                  &
                       !*(1-iflg)
!#endif

                  !conv=gasnew/densair(i,j,k)
                  !diffk=adiffk+bdiffk*temp    
                  !diffd=adiffd*(temp/t0)**bdiffd*ps0/(p0(k)+pnew(i,j,k))
                  !densfac=sqrt(ronull(1)/densair(i,j,k))
                  !itemp=nint(r1h-sign(r1h,t0-temp))
                  !!----------------------------------------------------------!
                  !! radius of rain and graupel at scalarpoints               !
                  !!----------------------------------------------------------!
                  !radwatp=max(radwatc,arwatp*(watpnew(i,j,k)/conv)**r1q)
                  !radgra =max(radice,  argra*(granew(i,j,k) /conv)**r1q)
                  !!----------------------------------------------------------!
                  !! compute particle number per kg total mass                !
                  !!----------------------------------------------------------!
                  !xnwatc=watcnew(i,j,k)/(rhowat*xpi*radwatc**3)
                  !xnwatp=watpnew(i,j,k)/(rhowat*xpi*radwatp**3)
                  !xnice = icenew(i,j,k)/(rhoice*xpi*radice**3)
        
                  !xswatc=watcnew(i,j,k)/radwatc
                  !xswatp=watpnew(i,j,k)/radwatp
                  !xsice = icenew(i,j,k)/radice
                  !xsgra = granew(i,j,k)/radgra
              
                  !xsall =max(epsmin,xswatc+xswatp+xsice+xsgra)
                  !xssol =max(epsmin,xsice+xsgra)
        
                  !xlamp=(watpnew(i,j,k)/conv)**r1q/alamp
                  !xlamg=(granew (i,j,k)/conv)**r1q/alamg
!#ifdef DEBUG
                  !xlamp=max(1.e-30_kreal,xlamp)
                  !xlamg=max(1.e-30_kreal,xlamg)
!#endif

                  !!----------------------------------------------------------!
                  !! super-/sub-saturation quantities of water vapor          !
                  !!----------------------------------------------------------!
                  !gas =rgasnew(i,j,k)
                  !pfak=gas/(gasf*(pnew(i,j,k)+p0(k)))*gasnew
        
                  !expon1=aice*(temp-t0)/max(temp-bice,epsmach)
                  !satpe=asat*exp(expon1)
        
                  !satwe =satpe*pfak
        
                  !supsate =wetnew(i,j,k)/satwe-r1
        
                  !subl=supsate/((heatsub/(gasf*temp)-r1)*heatsub/(diffk*temp) &
                       !+gasf*temp/(diffd*satpe))
                  !!----------------------------------------------------------!
                  !! limits and damping                                       !
                  !!----------------------------------------------------------!
                  !damp=max(adamp,min(bdamp,cdamp-ddamp*temp))
                  !damp=max( damp,min(edamp,fdamp-gdamp*temp))
                  !damptet=r3*cv/heatlat
                  !dwete =abs(wetnew(i,j,k)-satwe)*damp
                  !dwete =min(damptet,dwete)
        
                  !xssol =xssol+satemin*(r1h+sign(r1h,supsate))
                  !!----------------------------------------------------------!
                  !! 1. update of temperature and saturation pressure         !
                  !!----------------------------------------------------------!
                  !! compute temperature change                               !
                  !!----------------------------------------------------------!
                  !temp=heatsub*(p16(i,j)+p17(i,j))/cv+temp
                  !!----------------------------------------------------------!
                  !! thermodynamic quantities                                 !
                  !!----------------------------------------------------------!
                  !expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
                  !satpw=asat*exp(expon2)
                  !satww=satpw*pfak
        
                  !diffk=adiffk+bdiffk*temp    
                  !diffd=adiffd*(temp/t0)**bdiffd*ps0/(p0(k)+pnew(i,j,k))
                  !itemp=nint(r1h-sign(r1h,t0-temp))
                  !!----------------------------------------------------------!
                  !! super-/sub-saturation quantities of water vapor          !
                  !!----------------------------------------------------------!
                  !supsatw=wetnew(i,j,k)/satww-r1
        
                  !evap=supsatw/((heatlat/(gasf*temp)-r1)*heatlat/(diffk*temp) &
                       !+gasf*temp/(diffd*satpw))
                  !!----------------------------------------------------------!
                  !! limits and damping                                       !
                  !!----------------------------------------------------------!
                  !damp=max(adamp,min(bdamp,cdamp-ddamp*temp))
                  !damp=max( damp,min(edamp,fdamp-gdamp*temp))
        
                  !dwetw  =abs(wetnew(i,j,k)-max(r0,p16(i,j)+p17(i,j))-satww)*damp
                  !dwetw  =min(damptet,dwetw)
        
                  !xsall =xsall+r1th*(r1h+sign(r1h,supsatw))
                  !!----------------------------------------------------------!
                  !! 2. update of temperature and saturation pressure         !
                  !!----------------------------------------------------------!
                  !! compute temperature change                               !
                  !!----------------------------------------------------------!
                  !temp=heatlat*(p13(i,j)+p14(i,j)+p15(i,j))/cv+temp
                  !!----------------------------------------------------------!
                  !! thermodynamic quantities                                 !
                  !!----------------------------------------------------------!
                  !expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
                  !satpw=asat*exp(expon2)
                  !satww=satpw*pfak
        
                  !diffk=adiffk+bdiffk*temp    
                  !diffd=adiffd*(temp/t0)**bdiffd*ps0/(p0(k)+pnew(i,j,k))
                  !itemp=nint(r1h-sign(r1h,t0-temp))
                  !!----------------------------------------------------------!
                  !! limits and damping                                       !
                  !!----------------------------------------------------------!
                  !dampmel=cv/heatfrz*abs(t0-temp)
                  !dampmel=dampmel-p6(i,j)*(1-itemp)
                  !!----------------------------------------------------------!
                  !! 3. update of temperature and saturation pressure         !
                  !!----------------------------------------------------------!
                  !! compute temperature change                               !
                  !!----------------------------------------------------------!
                  !temp=heatfrz*(-p7(i,j)-p8(i,j)+p9(i,j)+p10(i,j)             &
                       !+p3(i,j)+p5(i,j)*(1-itemp)+p6(i,j))/cv+temp
                  !!----------------------------------------------------------!
                  !! thermodynamic quantities                                 !
                  !!----------------------------------------------------------!
                  !expon1=aice*(temp-t0)/max(temp-bice,epsmach)
                  !expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
                  !satpe=asat*exp(expon1)
                  !satpw=asat*exp(expon2)
                  !satwe=satpe*pfak
                  !satww=satpw*pfak
        
                  !itemp=nint(r1h-sign(r1h,t0-temp))
                  !!----------------------------------------------------------!
                  !! limits and damping                                       !
                  !!----------------------------------------------------------!
                  !damp=max(adamp,min(bdamp,cdamp-ddamp*temp))
                  !damp=max( damp,min(edamp,fdamp-gdamp*temp))
        
                  !dwetw =damp*max(r0,wetnew(i,j,k)                            &
                        !-max(r0,p13(i,j)+p14(i,j)+p15(i,j)+p16(i,j)+p17(i,j)) &
                        !-satww)
                  !dwetw =min(damptet,dwetw)
        
                  !dwete =damp*max(r0,wetnew(i,j,k)                            &
                        !-max(r0,p13(i,j)+p14(i,j)+p15(i,j)+p16(i,j)+p17(i,j)) &
                        !-satwe)
                  !dwete =min(damptet,dwete)
               !endif
            !enddo
         !enddo
       !end subroutine kessler_processes
       !!-------------------------------------------------------------------
       !subroutine kessler_fluxes(k,densair,p1,p2,p3,p4,p5,p6,p7,p8,p9,     &
                                 !p10,p11,p12,p13,p14,p15,p16,p17,p18,p19)
         !!----------------------------------------------------------------!
         !! computation of concentration changes during dt: fluxes         !
         !!----------------------------------------------------------------!
         !use phys_constants, only: r0, r1, r1h, epsmin
         !use phys_constants, only: cpair, ps0, pi
         !use phys_constants, only: heatlat, heatfrz, heatsub
         !use process_data,   only: wetnew,watcnew,watpnew,icenew,granew,   &
                                   !wetflx,watcflx,watpflx,iceflx,graflx
         !use process_data,   only: rhowat, rhoice, rhogra
         !use atham_module,   only: tetnew, density, cptot
         !use atham_module,   only: ntgas, ntrac, tgasnew, tracnew, cptgas
         !use atham_module,   only: tetflx, pflx
         !use atham_module,   only: x, xv, yv, iflgs
         !use atham_module,   only: tss
         !use atham_module,   only: cylindric_geometry, icenter

         !real(kreal), dimension(nx,ny,nz), intent(in) :: densair
         !real(kreal), dimension(nx,ny), intent(in) ::                      &
              !p1,p2,p3,p4,p5,p6,p7,p8,p9,                                  &
              !p10,p11,p12,p13,p14,p15,p16,p17,p18,p19
         !integer(kint), intent(in) :: k
         !!----------------------------------------------------------------!
         !! local variables                                                !
         !!----------------------------------------------------------------!
         !real(kreal) :: gasnew,drynew,cp,expon,pcom,temp,sumpos,sumneg,dflux, &
                        !wet_flx,watc_flx,watp_flx,ice_flx,gra_flx, dtet,      &
                        !gamma, dmair, dvair
         !real(kreal) :: dE_diag(7,nx,ny)
         !real(kreal) :: dyv, dxv, area, rad, radin
         !integer(kint) :: i,j,ir,il,jv,jh,iflg,itemp
     
         !real(kreal), parameter:: r3=3._kreal,r1hd=0.01_kreal,             &
              !t0=273.15_kreal
         !!----------------------------------------------------------------!
         !! start loop over grid points                                    !
         !!----------------------------------------------------------------!
         !do j=nyv,nyh
            !jh=min(nyh,j+1)
            !jv=j-1
            !do i=nxl,nxr
               !ir=min(nxr,i+1)
               !il=i-1
               !if (iflgs(i,j,k)==1) then
                  !gasnew=r1-sum(tracnew(i,j,k,1:ntrac))
                  !gasnew=max(r1hd,min(gasnew,r1))
                  !drynew=gasnew-sum(tgasnew(i,j,k,1:ntgas))
                  !cp=(cpair*drynew+sum(cptgas(1:ntgas)*tgasnew(i,j,k,1:ntgas))) &
                       !/gasnew
                  !expon=rgasnew(i,j,k)/cp
                  !pcom =((p0(k)+pnew(i,j,k))/ps0)**expon
!#ifndef OLD_TEMP
                  !temp=tempnew(i,j,k)
!#else
                  !iflg=r1h+sign(r1h,tempnew(i,j,k)-tss(k)*pcom+r3)
                  !temp=tempnew(i,j,k)*iflg                                    &
                      !+max(tempnew(i,j,k),                                    &
                          !(tempnew(i,j,k)+tempnew(ir,j,k)*iflgs(ir,j,k)       &
                                         !+tempnew(il,j,k)*iflgs(il,j,k))      &
                          !/(1+iflgs(ir,j,k)+iflgs(il,j,k)),                   &
                          !(tempnew(i,j,k)+tempnew(i,jh,k)*iflgs(i,jh,k)       &
                                         !+tempnew(i,jv,k)*iflgs(i,jv,k))      &
                          !/(1+iflgs(i,jh,k)+iflgs(i,jv,k)),                   &
                          !(tempnew(ir,j,k)*iflgs(ir,j,k)                      &
                          !+tempnew(il,j,k)*iflgs(il,j,k))                     &
                          !/max(1,iflgs(ir,j,k)+iflgs(il,j,k)),                &
                          !(tempnew(i,jh,k)*iflgs(i,jh,k)                      &
                          !+tempnew(i,jv,k)*iflgs(i,jv,k))                     &
                          !/max(1,iflgs(i,jh,k)+iflgs(i,jv,k)))                &
                         !*(1-iflg)
!#endif
                  !!----------------------------------------------------------!
                  !! itemp=1 warm                                             !
                  !! itemp=0 cold                                             !
                  !!----------------------------------------------------------!
                  !itemp=nint(r1h-sign(r1h,t0-temp))
                  !!----------------------------------------------------------!
                  !! concentration change of cloud water                      !
                  !!----------------------------------------------------------!
                  !watc_flx=-p1(i,j)-p5(i,j)+p8(i,j)-p9(i,j)                   &
                                 !-p11(i,j)+p13(i,j)+p18(i,j)
                  !watc_flx=min(watc_flx,wetnew(i,j,k)+icenew(i,j,k)*itemp)
                  !watc_flx=max(-watcnew(i,j,k),watc_flx)
                  !!----------------------------------------------------------!
                  !! concentration change of ice                              !
                  !!----------------------------------------------------------!
                  !ice_flx=-p2(i,j)-p4(i,j)-p8(i,j)                            &
                                !+p9(i,j)-p12(i,j)+p16(i,j)+p19(i,j)
                  !ice_flx=min(ice_flx,                                        &
                       !(wetnew(i,j,k)+watcnew(i,j,k)*(1-itemp)))
                  !ice_flx=max(-icenew(i,j,k),ice_flx)
                  !!----------------------------------------------------------!
                  !! concentration change of rain                             !
                  !!----------------------------------------------------------!
                  !watp_flx=p1(i,j)-p3(i,j)-p6(i,j)+p7(i,j)                    &
                       !-p10(i,j)+p11(i,j)+p14(i,j)+max(r0,p15(i,j))           &
                       !+(p2(i,j)+p5(i,j))*itemp
                  !watp_flx=min(watp_flx,                                      &
                       !wetnew(i,j,k)+watcnew(i,j,k)+granew(i,j,k)*itemp)
                  !watp_flx=max(-watpnew(i,j,k),watp_flx)
                  !!----------------------------------------------------------!
                  !! concentration change of graupel                          !
                  !!----------------------------------------------------------!
                  !gra_flx=p3(i,j)+p4(i,j)+p6(i,j)-p7(i,j)                     &
                       !+p10(i,j)+p12(i,j)+p17(i,j)+min(r0,p15(i,j))           &
                       !+(p2(i,j)+p5(i,j))*(1-itemp)
                  !gra_flx=min(gra_flx,                                        &
                       !(wetnew (i,j,k)+icenew (i,j,k)                         &
                       !+watpnew(i,j,k)+watcnew(i,j,k))*(1-itemp))
                  !gra_flx=max(-granew(i,j,k),gra_flx)
                  !!----------------------------------------------------------!
                  !! concentration change of watervapor                       !
                  !!----------------------------------------------------------!
                  !wet_flx=-p13(i,j)-p14(i,j)-p15(i,j)                         &
                       !-p16(i,j)-p17(i,j)-p18(i,j)-p19(i,j)
                  !wet_flx=max(-wetnew(i,j,k),wet_flx)
                  !!----------------------------------------------------------!
                  !! compute temperature change                               !
                  !!----------------------------------------------------------!
                  !dtet=(heatlat*(p13(i,j)+p14(i,j)+p18(i,j))                  &
                       !+heatsub*(p16(i,j)+p17(i,j)+p19(i,j))                  &
                       !+heatfrz*(-p7(i,j)-p8(i,j)+p9(i,j)+p10(i,j)            &
                                 !+p3(i,j)+p5(i,j)*(1-itemp)+p6(i,j)))         &
                                 !/(pcom*cptot(i,j,k))
                  !tetflx(i,j,k)=tetflx(i,j,k)+dtet
                  !!----------------------------------------------------------!
                  !! massfixer (only for total massfluxes)                    !
                  !!----------------------------------------------------------!
                  !sumpos=max(r0,wet_flx)+max(r0,watc_flx)                     &
                        !+max(r0,watp_flx)+max(r0,ice_flx)                     &
                        !+max(r0,gra_flx)
                  !sumneg=abs(min(r0,wet_flx)+min(r0,watc_flx)                 &
                            !+min(r0,watp_flx)+min(r0,ice_flx)                 &
                            !+min(r0,gra_flx))
                  !dflux=sumpos-sumneg
        
                  !sumpos=max(epsmin,sumpos)
                  !sumneg=max(epsmin,sumneg)

                  !wet_flx=wet_flx-max(dflux,r0)*max(r0,wet_flx)/sumpos        &
                                 !+min(dflux,r0)*min(r0,wet_flx)/sumneg
                  !wetflx(i,j,k)=wetflx(i,j,k)+wet_flx
                  
                  !watc_flx=watc_flx-max(dflux,r0)*max(r0,watc_flx)/sumpos     & 
                                   !+min(dflux,r0)*min(r0,watc_flx)/sumneg
                  !watcflx(i,j,k)=watcflx(i,j,k)+watc_flx

                  !watp_flx=watp_flx-max(dflux,r0)*max(r0,watp_flx)/sumpos     &
                                   !+min(dflux,r0)*min(r0,watp_flx)/sumneg
                  !watpflx(i,j,k)=watpflx(i,j,k)+watp_flx

                  !ice_flx=ice_flx-max(dflux,r0)*max(r0,ice_flx)/sumpos        &
                                 !+min(dflux,r0)*min(r0,ice_flx)/sumneg
                  !iceflx(i,j,k)=iceflx(i,j,k)+ice_flx

                  !gra_flx=gra_flx-max(dflux,r0)*max(r0,gra_flx)/sumpos        &
                                 !+min(dflux,r0)*min(r0,gra_flx)/sumneg
                  !graflx(i,j,k)=graflx(i,j,k)+gra_flx

                  !dmair=wet_flx/gasnew
                  !drvair=((watc_flx+watp_flx)/rhowat+ice_flx/rhoice+gra_flx/rhogra) &
                       !*density(i,j,k)*density(i,j,k)*gasnew/densair(i,j,k)
                  !gamma=cp/(cp-rgasnew(i,j,k))
                  !pflx(i,j,k)=pflx(i,j,k)+gamma*(p0(k)+pnew(i,j,k))*(dmair+dvair)
   !!!$               pflx(i,j,k)=pflx(i,j,k)+gamma*(p0(k)+pnew(i,j,k))*(dtet/tetnew(i,j,k)+dmair+dvair)

               !endif
            !enddo
         !enddo
         !!----------------------------------------------------------------!
         !! latent heat analysis                                           !
         !!----------------------------------------------------------------!
         !if (do_en_diag) then
            !do j=nyv,nyh
               !do i=nxl,nxr
                  !if (iflgs(i,j,k)==1) then
                     !itemp=nint(r1h-sign(r1h,t0-tempnew(i,j,k)))
                     
                     !dE_diag(1,i,j) = (heatlat*(p13(i,j)+p14(i,j)+p18(i,j))       &
                                      !+heatsub*(p16(i,j)+p17(i,j)+p19(i,j))       &
                                      !+heatfrz*(-p7(i,j)-p8(i,j)+p9(i,j)+p10(i,j) &
                                      !+p3(i,j)+p5(i,j)*(1-itemp)+p6(i,j)))*density(i,j,k)

                     !dE_diag(2,i,j) = heatlat*(min(r0,p13(i,j))                   &
                                              !+min(r0,p14(i,j))                   &
                                              !+min(r0,p18(i,j)))*density(i,j,k)
                     !dE_diag(3,i,j) = heatlat*(max(r0,p13(i,j))                   &
                                              !+max(r0,p14(i,j))                   &
                                              !+max(r0,p18(i,j)))*density(i,j,k)
                     !dE_diag(4,i,j) = heatfrz*(-p7(i,j)-p8(i,j))*density(i,j,k)
                     !dE_diag(5,i,j) = heatfrz*(p9(i,j)+p10(i,j)+p3(i,j)           &
                                     !+p5(i,j)*(1-itemp)+p6(i,j))*density(i,j,k)
                     !dE_diag(6,i,j) = heatsub*(min(r0,p16(i,j))                   &
                                              !+min(r0,p17(i,j))                   &
                                              !+min(r0,p19(i,j)))*density(i,j,k)
                     !dE_diag(7,i,j) = heatsub*(max(r0,p16(i,j))                   &
                                              !+max(r0,p17(i,j))                   &
                                              !+max(r0,p19(i,j)))*density(i,j,k)
                  !endif
               !enddo
            !enddo
            !if (cylindric_geometry) then 
               !j = nyv
               !i = icenter
               !if (iflgs(i,j,k)==1) then
                  !rad = xv(i) - x(icenter)
                  !area = pi*rad*rad
                  !en_diag(:,k) = en_diag(:,k) + dE_diag(:,i,j) * area
               !endif
               !do i=icenter+1,nxr   
                  !if (iflgs(i,j,k)==1) then
                     !rad = xv(i)-x(icenter)
                     !radin = xv(i-1)-x(icenter)
                     !area = pi*(rad*rad-radin*radin)
                     !en_diag(:,k) = en_diag(:,k) + dE_diag(:,i,j) * area
                  !endif
               !enddo
            !else
               !do  j=nyh,nyv
                  !dyv=yv(j)-yv(j-1)
                  !do i=nxl,nxr
                     !dxv=xv(i)-xv(i-1)
                     !if (iflgs(i,j,k)==1) then
                        !area=dxv*dyv
                        !en_diag(:,k) = en_diag(:,k) + dE_diag(:,i,j) * area
                     !endif
                  !enddo
               !enddo
            !endif
         !endif

       !end subroutine kessler_fluxes

    end subroutine
end module mphys_kessler_old
