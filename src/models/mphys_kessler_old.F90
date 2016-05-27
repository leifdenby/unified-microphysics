!> Almost verbatim reproduction of original "Kessler microphysics" from ATHAM,
!> previously found in `source/processModules/kessler_microphysics.F90`

module mphys_kessler_old

   implicit none
   public init

   integer, parameter :: real2 = selected_real_kind(3),            &
      real4 = selected_real_kind(6),            &
      real8 = selected_real_kind(12),           &
      kreal = real8
   integer, parameter :: int4 = selected_int_kind(7),              &
      int8 = selected_int_kind(13),             &
      kint = int4

   real(kreal), parameter ::                 &
      epsmach=1e-7_kreal,                  &
      epsmin =1.e-20_kreal              

   contains
      subroutine init()
         use microphysics_register, only: register_variable
         call register_variable('cloud_water', 1)
         call register_variable('rain', 1)
         call register_variable('cloud_ice', 1)
         call register_variable('graupel', 1)
      end subroutine

      !> The old "Kessler microphysics" method in ATHAM performed explicit split
      !> integration, so we expose a subroutine to integrate the microphysics
      !> processes here
      subroutine integrate_isobaric(y, t0, t_end)
         use microphysics_register, only: idx_cwater, idx_water_vapour, idx_rain, idx_temp, idx_pressure
         use microphysics_register, only: idx_cice, idx_graupel

         integer, parameter :: n_variables = 7
         real(8), intent(inout) :: y(n_variables)
         real(8), intent(in) :: t_end, t0

         real(kreal) :: p, temp
         real(kreal) :: q_v, q_l, q_r, q_i, q_gr
         ! state variable names used in ATHAM
         real(kreal) :: wetnew, watcnew, watpnew, icenew, granew
         real(kreal) :: wet_flx, watc_flx, watp_flx, ice_flx, gra_flx, pflx
         real(kreal) :: drynew, tempflx
         real(kreal) :: densair
         real(kreal) :: dt = 0.0

         dt = t_end - t0

         p = y(idx_pressure)
         temp = y(idx_temp)
         q_v = y(idx_water_vapour)
         q_l = y(idx_cwater)
         q_r = y(idx_rain)
         q_i = y(idx_cice)
         q_gr = y(idx_graupel)

         ! Assign to variable names used in original implementation
         watpnew = q_r
         watcnew = q_l
         granew = q_gr
         icenew = q_i
         wetnew = q_v

         !print *, "p", p
         !print *, "temp qv ql qr qi qgr", temp, q_v, q_l, q_r, q_i, q_gr
         !print *, "---->"

         call kessler_processes(p,temp, dt, &
            wetnew, watcnew, watpnew, icenew, granew, &
            wet_flx, watc_flx, watp_flx, ice_flx, gra_flx, tempflx, pflx)

         !print *, "-----"

         y(idx_temp) = temp + tempflx
         y(idx_pressure) = p + pflx
         y(idx_water_vapour) = wetnew + wet_flx
         y(idx_cwater) = watcnew + watc_flx
         y(idx_rain) = watpnew + watp_flx
         y(idx_cice) = icenew + ice_flx
         y(idx_graupel) = granew + gra_flx

         !print *, "dtemp dqv dql dqr dqi dqgr", tempflx, wet_flx, watc_flx, watp_flx, ice_flx, gra_flx

      end subroutine integrate_isobaric


      subroutine kessler_processes(p,temp, dt, &
            wetnew, watcnew, watpnew, icenew, granew, &
            wet_flx, watc_flx, watp_flx, ice_flx, gra_flx, tempflx, pflx)

         real(kreal), intent(inout) :: temp
         real(kreal), intent(in) :: p !< Pressure [Pa]
         real(kreal), intent(in) :: dt
         real(kreal), intent(in) :: wetnew, watcnew, watpnew, icenew, granew
         real(kreal), intent(out) :: wet_flx, watc_flx, ice_flx, watp_flx, gra_flx, tempflx, pflx

         real(kreal), parameter ::                  &
            r0=0._kreal,                          &
            r1=1._kreal,                          &
            r2=2._kreal,                          &
            r1h=0.5_kreal,                        &
            r3h=1.5_kreal,                        &
            r1q=0.25_kreal,                       &
            r3q=0.75_kreal

         real(kreal), parameter ::                 &
            epsmach=1e-7_kreal,                  &
            epsmin =1.e-20_kreal              

         ! Physical constants from `Atham/phys_constants.F90`
         real(kreal), parameter ::                 &
            cpair  =1004.64_kreal,               &   ! Air = dry air
            cvair  = 717.60_kreal,               &
            rgasair=cpair-cvair,                 & 
            cpco2  =844._kreal,                  &   ! CO2
            cvco2  =655._kreal,                  &
            rgasco2=cpco2-cvco2,                 &
            cphum  =1864._kreal,                 &   ! H2O
            cvhum  =1402.5_kreal,                &
            rgashum=cphum-cvhum,                 &
            cpso2  =655.6_kreal,                 &   ! SO2
            cvso2  =512.0_kreal,                 &
            rgasso2=cpso2-cvso2

         real(kreal), parameter :: ps0    =101300._kreal

         real(kreal), parameter ::                 &
            heatlat=2500800._kreal,              &
            heatfrz= 335800._kreal,              &
            heatsub=2836600._kreal

         ! Process constants
         real(kreal), parameter ::                                         &
            r3=3._kreal,r4=4._kreal,r5=5._kreal,                         &
            r1th=0.1_kreal,r1hd=0.01_kreal,                              &
            arwatp=4.35e-3_kreal,argra=2.57e-2_kreal,                    &
            alamp=421._kreal,alamg=71.8_kreal,                           &
            adiffk=2.4e-2_kreal,bdiffk=8.e-5_kreal,                      &
            adiffd=2.11e-5_kreal,bdiffd=1.94_kreal,                      &
            adamp=0.1_kreal,bdamp=0.27_kreal,                            &
            cdamp=1.884_kreal,ddamp=0.0054_kreal,                        &
            edamp=0.9_kreal,fdamp=4.395_kreal,gdamp=0.015_kreal,         &
            satemin=1.e-4_kreal,                                         &
            ap16=6.4_kreal,                                              &
            ap17=1.22e5_kreal,bp17=4.11e5_kreal,pexp=2.75_kreal,         &
            ap14=6.28e7_kreal,bp14=2.65e8_kreal,                         &
            ap15=1.22e5_kreal,bp15=4.11e5_kreal,                         &
            azylr=3.71e9_kreal,arexp=3.5_kreal,                          &
            bzylr=1.88e11_kreal,czylr=5.61e15_kreal,                     &
            apacc=3.56e6_kreal,                                          &
            ap7=7.603e4_kreal,bp7=1.6_kreal,cp7=5.36_kreal,              &
            ap10=1.974e14_kreal,bp10=0.66_kreal,                         &
            ap9=1.76e6_kreal,bp9=0.66_kreal,                             &
            ap11=1.e-3_kreal,bp11=5.e-4_kreal,                           &
            ap12=1.e-3_kreal,bp12=2.5e-2_kreal,                          &
            ap18=1.e-14_kreal,bp18=0.6_kreal

         ! original atham values:
         real(kreal), parameter ::                                         &
            t0=273.15_kreal,  asat=610.7_kreal,                          &
            awat=17.25_kreal, bwat=36._kreal,                            &
            aice=22.33_kreal, bice=2._kreal

         real(kreal), parameter :: gasf=cphum-cvhum
         real(kreal), parameter :: pi  =r4*atan(r1)
         real(kreal), parameter :: xpi =r4/r3*pi

         !----------------------------------------------------------------!
         ! local variables                                                !
         !----------------------------------------------------------------!
         real(kreal) :: gasnew,drynew,cp,cv,expon,   &
            conv,diffk,diffd,densfac,radwatp,radgra,                     &
            xnwatc,xnwatp,xnice,xswatc,xswatp,xsice,xsgra,xsall,xssol,   &
            xniceac,xlamp,xlamg,gas,pfak,expon1,satpe,satwe,supsate,subl,&
            damp,dwete,expon2,satpw,satww,supsatw,evap,dwetw,zylr,       &
            dampmel,damptet

         real(kreal) :: pcom

         real(kreal) :: p1 ,p2 ,p3 ,p4 ,p5 ,p6 ,p7 ,p8 ,p9,                          &
            p10,p11,p12,p13,p14,p15,p16,p17,p18,p19

         integer(kint) :: iflg,itemp


         ! Extra constants from ATHAM
         real(kreal), parameter :: ronull = 1.225_kreal
         real(kreal), parameter :: radwatc = 0.00001_kreal
         real(kreal), parameter :: radice = 0.00001_kreal
         real(kreal), parameter :: rhowat = 1000._kreal
         real(kreal), parameter :: rhoice = 917._kreal
         real(kreal), parameter :: rhogra = 917._kreal
         real(kreal), parameter :: cpwat = 4183._kreal
         ! extra needed for mixture heat capacity, taken from `Configurations/convection.F90`
         real(kreal), parameter :: cp_l = 4183._kreal
         real(kreal), parameter :: cp_i = 2103._kreal

         ! TODO: fix these by actually having a velocity of the hydrometeors
         real(kreal), parameter :: wgra = 0.0_kreal, wwatp = 0.0_kreal
         ! TODO: fix by having non-zero viscosity read in from an external profile
         real(kreal), parameter :: viskos = 0.0_kreal

         ! Variable names used in ATHAM
         real(kreal) :: rgasnew
         real(kreal) :: densair, density
         real(kreal) :: dtemp
         real(kreal) :: cptot
         real(kreal) :: sumpos, sumneg, dflux, dmair, dvair

         ! local variables to easy variable ATHAM variable calculation
         real(kreal) :: cp_g, cv_g, cp_m, gamma

         drynew = 1.0 - wetnew - watcnew - watpnew - icenew - granew

         gasnew = drynew + wetnew

         cp_g = (drynew*cpair + wetnew*cphum)/gasnew
         cv_g = (drynew*cvair + wetnew*cvhum)/gasnew
         cp_m = drynew*cpair + wetnew*cphum + (watcnew + watpnew)*cp_l + (icenew + granew)*cp_i

         rgasnew = cp_g - cv_g
         cp = cp_g
         cv = cv_g
         cptot = cp_m

         densair = p/(rgasnew*temp)
         density = (drynew*rgasair + wetnew*rgashum)*temp/p + (watcnew+watpnew)/rhowat + icenew/rhoice + granew/rhogra


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

         ! reset local variables
         p1 = 0.0
         p2 = 0.0
         p3 = 0.0
         p4 = 0.0
         p5 = 0.0
         p6 = 0.0
         p7 = 0.0
         p8 = 0.0
         p9 = 0.0
         p10 = 0.0
         p11 = 0.0
         p12 = 0.0
         p13 = 0.0
         p14 = 0.0
         p15 = 0.0
         p16 = 0.0
         p17 = 0.0
         p18 = 0.0
         p19 = 0.0
         pcom = 0.0

         watc_flx = 0.0
         wet_flx = 0.0
         watc_flx = 0.0
         ice_flx = 0.0
         watp_flx = 0.0
         gra_flx = 0.0
         dtemp = 0.0
         tempflx = 0.0

         ! Start of old "Kessler" microphysics from ATHAM
         gasnew=max(r1hd,min(gasnew,r1))

         conv=gasnew/densair
         ! NB: Fixed bug in thermal conductivity parameterisation
         !diffk=adiffk+bdiffk*temp
         diffk=adiffk+bdiffk*(temp-t0)
         diffd=adiffd*(temp/t0)**bdiffd*ps0/p
         densfac=sqrt(ronull/densair)
         itemp=nint(r1h-sign(r1h,t0-temp))
         !----------------------------------------------------------!
         ! radius of rain and graupel at scalarpoints               !
         !----------------------------------------------------------!
         radwatp=max(radwatc,arwatp*(watpnew/conv)**r1q)
         radgra =max(radice,  argra*(granew /conv)**r1q)
         !----------------------------------------------------------!
         ! compute particle number per kg total mass                !
         !----------------------------------------------------------!
         xnwatc=watcnew/(rhowat*xpi*radwatc**3)
         xnwatp=watpnew/(rhowat*xpi*radwatp**3)
         xnice = icenew/(rhoice*xpi*radice**3)

         xswatc=watcnew/radwatc
         xswatp=watpnew/radwatp
         xsice = icenew/radice
         xsgra = granew/radgra

         xsall =max(epsmin,xswatc+xswatp+xsice+xsgra)
         xssol =max(epsmin,xsice+xsgra)

         xlamp=(watpnew/conv)**r1q/alamp
         xlamg=(granew /conv)**r1q/alamg

         !----------------------------------------------------------!
         ! super-/sub-saturation quantities of water vapor          !
         !----------------------------------------------------------!
         gas =rgasnew

         ! TODO: I believe there should *not* be /gasnew here
         pfak=gas/(gasf*p)*gasnew

         expon1=aice*(temp-t0)/max(temp-bice,epsmach)
         satpe=asat*exp(expon1)

         satwe =satpe*pfak

         supsate =wetnew/satwe-r1

         subl=supsate/((heatsub/(gasf*temp)-r1)*heatsub/(diffk*temp) &
            +gasf*temp/(diffd*satpe))
         !----------------------------------------------------------!
         ! limits and damping                                       !
         !----------------------------------------------------------!
         damp=max(adamp,min(bdamp,cdamp-ddamp*temp))
         damp=max( damp,min(edamp,fdamp-gdamp*temp))
         damptet=r3*cv/heatlat
         dwete =abs(wetnew-satwe)*damp
         dwete =min(damptet,dwete)

         xssol =xssol+satemin*(r1h+sign(r1h,supsate))
         !----------------------------------------------------------!
         ! sublimation and deposition on particels                  !
         !----------------------------------------------------------!
         ! subl/dep on ice                                          !
         !----------------------------------------------------------!
         p16=ap16*pi*xnice*radice*subl
         p16=sign(r1,p16)                                  &
            *min(abs(p16*dt),xsice*dwete/xssol)*(1-itemp)
         p16=max(p16,-icenew)
         !----------------------------------------------------------!
         ! subl/dep on graupel                                      !
         !----------------------------------------------------------!
         ! TODO: `viskos` should come from an external profile
         p17=conv*subl*(ap17*xlamg*xlamg                        &
            +bp17*xlamg**pexp*sqrt(densfac/viskos))
         p17=sign(r1,p17)                                  &
            *min(abs(p17*dt),xsgra*dwete/xssol)*(1-itemp)
         p17=max(p17,-granew)
         !----------------------------------------------------------!
         ! 1. update of temperature and saturation pressure         !
         !----------------------------------------------------------!
         ! compute temperature change                               !
         !----------------------------------------------------------!
         temp=heatsub*(p16+p17)/cv+temp
         !----------------------------------------------------------!
         ! thermodynamic quantities                                 !
         !----------------------------------------------------------!
         expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
         satpw=asat*exp(expon2)
         satww=satpw*pfak

         ! NB: Fixed bug in thermal conductivity parameterisation
         diffk=adiffk+bdiffk*(temp - t0)
         !diffk=adiffk+bdiffk*temp    
         diffd=adiffd*(temp/t0)**bdiffd*ps0/p
         itemp=nint(r1h-sign(r1h,t0-temp))
         !----------------------------------------------------------!
         ! super-/sub-saturation quantities of water vapor          !
         !----------------------------------------------------------!
         supsatw=wetnew/satww-r1


         evap=supsatw/((heatlat/(gasf*temp)-r1)*heatlat/(diffk*temp) &
            +gasf*temp/(diffd*satpw))
         !----------------------------------------------------------!
         ! limits and damping                                       !
         !----------------------------------------------------------!
         damp=max(adamp,min(bdamp,cdamp-ddamp*temp))
         damp=max( damp,min(edamp,fdamp-gdamp*temp))

         dwetw  =abs(wetnew-max(r0,p16+p17)-satww)*damp
         dwetw  =min(damptet,dwetw)

         xsall =xsall+r1th*(r1h+sign(r1h,supsatw))
         !----------------------------------------------------------!
         ! condensation and evaporation on droplets                 ! 
         !----------------------------------------------------------!
         ! cond/evap on cloud                                       !
         !----------------------------------------------------------!
         p13=r4*pi*xnwatc*radwatc*evap
         p13=sign(r1,p13)*min(abs(p13*dt),xswatc/xsall*dwetw)
         p13=max(p13,-watcnew)
         !----------------------------------------------------------!
         ! cond/evap on rain                                        !
         !----------------------------------------------------------!
         ! TODO: non-zero `viskos` based on height `viskos(k)` previously
         p14=conv*evap*(ap14*xlamp*xlamp                        &
            +bp14*xlamp**pexp*sqrt(densfac/viskos))
         p14=sign(r1,p14)*min(abs(p14*dt),xswatp/xsall*dwetw)
         p14=max(p14,-watpnew)
         !----------------------------------------------------------!
         ! cond/evap on melting graupel                             !
         !----------------------------------------------------------!
         ! TODO: non-zero `viskos` based on height `viskos(k)` previously
         p15=conv*evap*(ap15*xlamg*xlamg                        &
            +bp15*xlamg**pexp*sqrt(densfac/viskos))
         p15=sign(r1,p15)                                  &
            *min(abs(p15*dt),xsgra/xsall*dwetw)*itemp
         p15=max(p15,-granew)

         !----------------------------------------------------------!
         ! 2. update of temperature and saturation pressure         !
         !----------------------------------------------------------!
         ! compute temperature change                               !
         !----------------------------------------------------------!
         temp=heatlat*(p13+p14+p15)/cv+temp
         !----------------------------------------------------------!
         ! thermodynamic quantities                                 !
         !----------------------------------------------------------!
         expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
         satpw=asat*exp(expon2)
         satww=satpw*pfak

         ! NB: Fixed bug in thermal conductivity parameterisation
         diffk=adiffk+bdiffk*(temp-t0)
         !diffk=adiffk+bdiffk*temp    
         diffd=adiffd*(temp/t0)**bdiffd*ps0/p
         itemp=nint(r1h-sign(r1h,t0-temp))
         !----------------------------------------------------------!
         ! limits and damping                                       !
         !----------------------------------------------------------!
         dampmel=cv/heatfrz*abs(t0-temp)
         !----------------------------------------------------------!
         ! accretion                                                !
         !----------------------------------------------------------!
         ! accretion of cloud droplets by rain                      !
         !----------------------------------------------------------!
         zylr=min(azylr*xlamp**arexp,bzylr*xlamp**4,czylr*xlamp**5)
         p1=watcnew*zylr*densfac
         p1=max(r0,min(p1*dt,watcnew))
         !----------------------------------------------------------!
         ! accretion of ice particles by rain forming rain    at t>t0
         !                                    forming graupel at t<t0
         !----------------------------------------------------------!
         p2=icenew*zylr*densfac
         p2=max(r0,min(p2*dt,icenew)*(1-itemp))
         !----------------------------------------------------------!
         ! change of rainconcentration at t<t0                      !
         !----------------------------------------------------------!
         xniceac=p2/(rhoice*xpi*radice**3)
         p3=min(r1,xniceac/max(xnwatp,epsmin))*watpnew
         p3=max(r0,min(p3,watpnew,r1h*dampmel))
         dampmel=dampmel-p3
         !----------------------------------------------------------!
         ! accretion of ice particles by graupel                    !
         !----------------------------------------------------------!
         p4=exp(0.05*(temp-t0))*icenew*densfac*apacc*xlamg**arexp
         p4=max(r0,min(p4*dt,icenew)*(1-itemp))
         !----------------------------------------------------------!
         ! accretion of clouddroplets by graupel forming graupel at t<t0
         !                                       forming rain    at t>t0
         !----------------------------------------------------------!
         p5=watcnew*densfac*apacc*xlamg**arexp
         p5=min(p5*dt,watcnew-p1*itemp)
         p5=max(r0,p5*itemp+min(p5,r1h*dampmel)*(1-itemp))
         dampmel=dampmel-p5*(1-itemp)
         !----------------------------------------------------------!
         ! accretion of rain by graupel forming graupel t<t0        !
         !----------------------------------------------------------!
         p6=1.194e15*abs(wgra-wwatp)*conv         &
            *(r5 *xlamp**6*xlamg                                   &
            +r2 *xlamp**5*xlamg*xlamg                              &
            +r1h*xlamp**4*xlamg**3)
         p6=min(p6*dt,watpnew-p3)
         p6=max(r0,p6*itemp+min(p6,r1h*dampmel)*(1-itemp))
         dampmel=dampmel-p6*(1-itemp)
         !----------------------------------------------------------!
         ! melting and freezing                                     !
         !----------------------------------------------------------!
         ! melting of graupel to rain                               !
         !----------------------------------------------------------!
         ! TODO: non-zero `viskos` based on height `viskos(k)` previously
         p7=-p15*heatlat/heatfrz                               &
            +ap7*conv/heatfrz*diffk*(temp-t0)                          &
            *(bp7*xlamg*xlamg+cp7*xlamg**pexp*sqrt(densfac/viskos)) &
            +cpwat/heatfrz*(temp-t0)*(p5+p6)
         p7=max(r0,min(p7*dt,granew,r1h*dampmel)*itemp)
         p6=p6*(1-itemp)
         dampmel=dampmel-p7
         !----------------------------------------------------------!
         ! melting total ice to cloudwater  at t0                   !
         !----------------------------------------------------------!
         p8=max(r0,min(dampmel,max(icenew-p4,r0))*itemp)
         !----------------------------------------------------------!
         ! statistical freezing of rain to graupel                  !
         !----------------------------------------------------------!
         p10=conv*ap10*xlamp**7*(exp(bp10*max(r0,t0-temp))-r1)
         p10=max(r0,min(p10*dt,                            &
            watpnew-p3-p6,r1h*dampmel))
         dampmel=dampmel-p10
         !----------------------------------------------------------!
         ! statistical freezing of cloud to ice                     !
         !----------------------------------------------------------!
         p9=ap9*xnwatc*radwatc**6*(exp(bp9*max(r0,t0-temp))-r1)
         p9=max(r0,min(p9*dt,watcnew-p2,dampmel))
         !----------------------------------------------------------!
         ! autoconversion                                           !
         !----------------------------------------------------------!
         ! autoconversion of clouddroplets to rain (kessler)        !
         !----------------------------------------------------------!
         p11=ap11*max(r0,watcnew-bp11*conv)
         p11=max(r0,min(p11*dt,watcnew))
         !----------------------------------------------------------!
         ! autoconversion of icepartiles in graupel (lin et al.)    !
         !----------------------------------------------------------!
         p12=ap12*exp(bp12*(temp-t0))*max(r0,icenew-ap12*conv)
         p12=max(r0,min(p12*dt,icenew)*(1-itemp))
         !----------------------------------------------------------!
         ! 3. update of temperature and saturation pressure         !
         !----------------------------------------------------------!
         ! compute temperature change                               !
         !----------------------------------------------------------!
         temp=heatfrz*(-p7-p8+p9+p10             &
            +p3+p5*(1-itemp)+p6)/cv+temp
         !----------------------------------------------------------!
         ! thermodynamic quantities                                 !
         !----------------------------------------------------------!
         expon1=aice*(temp-t0)/max(temp-bice,epsmach)
         expon2=awat*(temp-t0)/max(temp-bwat,epsmach)
         satpe=asat*exp(expon1)
         satpw=asat*exp(expon2)
         satwe=satpe*pfak
         satww=satpw*pfak

         itemp=nint(r1h-sign(r1h,t0-temp))
         !----------------------------------------------------------!
         ! limits and damping                                       !
         !----------------------------------------------------------!
         damp=max(adamp,min(bdamp,cdamp-ddamp*temp))
         damp=max( damp,min(edamp,fdamp-gdamp*temp))

         dwetw =damp*max(r0,wetnew                            &
            -max(r0,p13+p14+p15+p16+p17) &
            -satww)
         dwetw =min(damptet,dwetw)

         dwete =damp*max(r0,wetnew                            &
            -max(r0,p13+p14+p15+p16+p17) &
            -satwe)
         dwete =min(damptet,dwete)
         !----------------------------------------------------------!
         ! new formation of cloud droplets                          !
         !----------------------------------------------------------!
         ! remaining supersaturation water vapor forms iceparitcles !
         !----------------------------------------------------------!
         xnice = icenew/(rhoice*xpi*radice**3)
         p19=max(r0,ap18*conv*exp(bp18*max(r0,t0-temp))-xnice)
         p19=(1-itemp)*min(p19,dwete)
         !----------------------------------------------------------!
         ! remaining supersaturation water vapor forms clouddroplets!
         !----------------------------------------------------------!
         p18=max(r0,dwetw-p19)


         ! ==================================
         ! Start of original `kessler_fluxes`

         expon=rgasnew/cp
         pcom =(p/ps0)**expon

         !----------------------------------------------------------!
         ! itemp=1 warm                                             !
         ! itemp=0 cold                                             !
         !----------------------------------------------------------!
         itemp=nint(r1h-sign(r1h,t0-temp))
         !----------------------------------------------------------!
         ! concentration change of cloud water                      !
         !----------------------------------------------------------!
         watc_flx=-p1-p5+p8-p9                   &
            -p11+p13+p18
         watc_flx=min(watc_flx,wetnew+icenew*itemp)
         watc_flx=max(-watcnew,watc_flx)
         !----------------------------------------------------------!
         ! concentration change of ice                              !
         !----------------------------------------------------------!
         ice_flx=-p2-p4-p8                            &
            +p9-p12+p16+p19
         ice_flx=min(ice_flx,                                        &
            (wetnew+watcnew*(1-itemp)))
         ice_flx=max(-icenew,ice_flx)
         !----------------------------------------------------------!
         ! concentration change of rain                             !
         !----------------------------------------------------------!
         watp_flx=p1-p3-p6+p7                    &
            -p10+p11+p14+max(r0,p15)           &
            +(p2+p5)*itemp
         watp_flx=min(watp_flx,                                      &
            wetnew+watcnew+granew*itemp)
         watp_flx=max(-watpnew,watp_flx)
         !----------------------------------------------------------!
         ! concentration change of graupel                          !
         !----------------------------------------------------------!
         gra_flx=p3+p4+p6-p7                     &
            +p10+p12+p17+min(r0,p15)           &
            +(p2+p5)*(1-itemp)
         gra_flx=min(gra_flx,                                        &
            (wetnew +icenew                          &
            +watpnew+watcnew)*(1-itemp))
         gra_flx=max(-granew,gra_flx)
         !----------------------------------------------------------!
         ! concentration change of watervapor                       !
         !----------------------------------------------------------!
         wet_flx=-p13-p14-p15                         &
            -p16-p17-p18-p19
         wet_flx=max(-wetnew,wet_flx)
         !----------------------------------------------------------!
         ! compute temperature change                               !
         !----------------------------------------------------------!
         dtemp=(heatlat*(p13+p14+p18)                  &
            +heatsub*(p16+p17+p19)                  &
            +heatfrz*(-p7-p8+p9+p10            &
            +p3+p5*(1-itemp)+p6))         &
            /cptot
         tempflx=dtemp

         !dtet=(heatlat*(p13+p14+p18)                  &
            !+heatsub*(p16+p17+p19)                  &
            !+heatfrz*(-p7-p8+p9+p10            &
            !+p3+p5*(1-itemp)+p6))         &
            !/(pcom*cptot)
         !tetflx=tetflx+dtemp

         !----------------------------------------------------------!
         ! massfixer (only for total massfluxes)                    !
         !----------------------------------------------------------!
         sumpos=max(r0,wet_flx)+max(r0,watc_flx)                     &
            +max(r0,watp_flx)+max(r0,ice_flx)                     &
            +max(r0,gra_flx)
         sumneg=abs(min(r0,wet_flx)+min(r0,watc_flx)                 &
            +min(r0,watp_flx)+min(r0,ice_flx)                 &
            +min(r0,gra_flx))
         dflux=sumpos-sumneg

         sumpos=max(epsmin,sumpos)
         sumneg=max(epsmin,sumneg)

         wet_flx=wet_flx-max(dflux,r0)*max(r0,wet_flx)/sumpos        &
            +min(dflux,r0)*min(r0,wet_flx)/sumneg

         watc_flx=watc_flx-max(dflux,r0)*max(r0,watc_flx)/sumpos     & 
            +min(dflux,r0)*min(r0,watc_flx)/sumneg

         watp_flx=watp_flx-max(dflux,r0)*max(r0,watp_flx)/sumpos     &
            +min(dflux,r0)*min(r0,watp_flx)/sumneg

         ice_flx=ice_flx-max(dflux,r0)*max(r0,ice_flx)/sumpos        &
            +min(dflux,r0)*min(r0,ice_flx)/sumneg

         gra_flx=gra_flx-max(dflux,r0)*max(r0,gra_flx)/sumpos        &
            +min(dflux,r0)*min(r0,gra_flx)/sumneg

         ! TODO: work out what to do with the pressure correction calculated here
         dmair=wet_flx/gasnew
         dvair=((watc_flx+watp_flx)/rhowat+ice_flx/rhoice+gra_flx/rhogra) &
         *density*density*gasnew/densair
         gamma=cp/(cp-rgasnew)
         pflx=pflx+gamma*p*(dmair+dvair)
         !$               pflx=pflx+gamma*(p0(k)+pnew)*(dtemp/tetnew+dmair+dvair)

      end subroutine kessler_processes

      function dydt(t, y)
         ! use microphysics_register, only: idx_cwater, idx_water_vapour
         ! use microphysics_constants, only: L_cond
         integer, parameter :: n_species = 8

         real(8), intent(in) :: t
         real(8), dimension(n_species), intent(in) :: y
         real(8), dimension(n_species) :: dydt

         print *, "Not implemented"
         call exit(-1)
      end function
   end module 
