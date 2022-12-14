#include "rundeck_opts.h"
!NOT FOR HYCOM:dts,mo,dxypo, n_abioDIC, num_tracers,trmo and oij passed to the subroutine.
      subroutine obio_carbon(gro,vrbos,kmax,i,j,nstep,kdm,n_co2n,
     &                       DTS,mmo,ddxypo,n_abioDIC,num_tracers,SDIC)
c
c  Computes carbon cycling.  Uses Aumont et al (2002; JGR) for
c  semi-labile DOC (because of basic similarities in model
c  formulation, specifically, presense of grazers and particulate
c  detritus).
c  ncar(1) = DOC semi-labile (uM(C))
c  ncar(2) = DIC  (uM(C))
c

      USE MODEL_COM, only: dtsrc
      USE CONSTANT, only: tf
      USE obio_dim
      USE obio_incom, only : cnratio,rlamdoc,rkdoc1,rkdoc2
     .                      ,rlampoc,uMtomgm3,Pzo,stdslp
     .                      ,excz,resz,remin,excp,resp,bn,cchlratio
     .                      ,mgchltouMC,bf,NCrrat,HvO2,O2thr
      USE obio_forc, only: wind,tirrq
      USE obio_com, only : C_tend,obio_P,P_tend,car
     .                    ,tfac,det,D_tend,tzoo,pnoice,pCO2_ij,pHsfc
     .                    ,temp1d,saln1d,dp1d,rhs,alk1d,trmo_unit_factor
     .                    ,rho1d,dicresp,docbac          !@PL added rho1d,docbac,dicresp
#ifdef TRACERS_bio_O2
     .                    ,o21d
#endif
#ifdef TRACERS_Alkalinity
      use obio_com, only: co3_conc
#endif
#ifdef OBIO_mocsy
      use ofluxes, only: oAPRESS
      use ocean,   only: oLAT_DG,r3d
      use obio_incom, only: npratio
      use obio_com, only: p1d
      use mvars
      use obio_diag, only : oijl=>obio_ijl,ijl_omegaA,ijl_omegaC
#endif


      use obio_com, only: co2flux

      use TimeConstants_mod, only: SECONDS_PER_HOUR, DAYS_PER_YEAR,
     &                             HOURS_PER_DAY
      

      use runtimecontrols_mod, only: constco2, pco2_online
      use dictionary_mod, only: get_param
      use domain_decomp_1d, only: am_i_root

      implicit none

      integer,intent(in) :: kdm,nstep,n_co2n,n_abioDIC,num_tracers
      real, intent(in) :: dts,mmo,ddxypo
      real,intent(inout) :: SDIC

!     real, parameter :: awan=0.337d0/(3.6d5) !piston vel coeff., from
!                                             !Wanninkof 1992, but adjusted
!                                             !by OCMIP, and converted from
!                                             !cm/hr to m/s
      real, parameter :: awan=0.251d0/(3.6d5) !piston vel coeff., OCMIP2016
                                              !converted from
                                              !cm/hr to m/s
      integer :: i,j,k

      integer :: nt,kmax
      real  :: rmmzoo,docexcz,rndep,docdep
      real  :: docdet,dicresz,sumdoc,sumutk,sumres,totgro
      real  :: docexcp(nchl),scco2,scco2arg,wssq,rkwco2
      real  :: Ts,tk,tk100,tk1002,ff,xco2,deltco2,flxmolm3,flxmolm3h
      real  :: gro(kdm,nchl)
      real term,sdic_uM,pCO2_abio,dummy
      real bs
      real, save :: atmco2=-1.

      logical vrbos

#ifdef OBIO_mocsy
      real mocsy_ph(kmax), mocsy_pco2(kmax), mocsy_fco2(kmax), 
     .     mocsy_co2(kmax),mocsy_hco3(kmax), mocsy_co3(kmax), 
     .     mocsy_OmegaA(kmax),mocsy_OmegaC(kmax), mocsy_BetaD(kmax), 
     .     mocsy_rhoSW(kmax),mocsy_p(kmax), mocsy_tempis(kmax),
     .     patm(kmax),mocsylat(kmax)
      character(10) :: optCON, optT, optP, optB, optKf, optK1K2
#endif

      bs = 2.d0*bn

! Detrital, bacterial, and grazing components
      do k = 1,kmax
!change: March 15, 2010
!       cchlratio = bn*cnratio
!       mgchltouMC = cchlratio/uMtomgm3


!@PL ratio of o2/nitrate remineralization
!@PL define delta function for O2
         NCrrat = 1.d0
         HvO2 = 1.d0
#ifdef TRACERS_bio_O2
        if (o21d(k).le.O2thr) then
         NCrrat = 0.4d0
         HvO2 = 0.d0
        else
         NCrrat = 1.d0
         HvO2 = 1.d0
        endif
#endif

        !---------------------------------------------------------------
        !DOC
        rmmzoo = obio_P(k,ntyp)/(Pzo+obio_P(k,ntyp))

        docexcz = excz*rmmzoo*obio_P(k,ntyp)  !zoopl production DOC
#ifdef apply_newnadj
        docexcz = 1.5d0*docexcz
#endif
        term = - docexcz*pnoice(k)
        rhs(k,9,14) = term
        P_tend(k,ntyp) = P_tend(k,ntyp) + term

!change: June 1, 2010
        term = bn*docexcz*pnoice(k)
#ifdef apply_nadj
        term = 1.5d0*term
#endif
        rhs(k,1,9) = term
        P_tend(k,1) = P_tend(k,1) + term

        term = bf*docexcz*pnoice(k)
        rhs(k,4,10) = term
        P_tend(k,4) = P_tend(k,4) + term
!endofchange

        rndep = rlamdoc*(obio_P(k,1)/(rkdoc1 + obio_P(k,1)))
        docdep = car(k,1)/(rkdoc2+car(k,1))
        docbac(k) = NCrrat*tfac(k)*rndep*docdep*car(k,1)   !bacterial loss DOC
        docdet = tfac(k)*rlampoc*det(k,1)        !detrital production DOC

!!!!    term = (docexcz*mgchltouMC
!!!! .                   +  docdet/uMtomgm3-docbac)*pnoice(k)
        term = docexcz*mgchltouMC*pnoice(k)
        rhs(k,13,12) = term
        C_tend(k,1) = C_tend(k,1) + term

        term = docdet/uMtomgm3   *pnoice(k)
        rhs(k,13,13) = term
        C_tend(k,1) = C_tend(k,1) + term


        term = -docbac(k)*pnoice(k)
        rhs(k,13,14) = term
        C_tend(k,1) = C_tend(k,1) + term

        !adjust detritus
        term = - docdet * pnoice(k) !carbon/nitrogen detritus
        rhs(k,10,14) = term
        D_tend(k,1) = D_tend(k,1) + term  

        !equivalent amount of DOC from detritus goes into 
        !nitrate, bypassing DON
        term = docdet/cnratio * pnoice(k)
        rhs(k,1,14) = term
        P_tend(k,1) = P_tend(k,1) + term 

        !---------------------------------------------------------------
        !DIC
        dicresz = tzoo*HvO2*resz*obio_P(k,ntyp) !zoopl production DIC (resp)
        term = - dicresz*pnoice(k)
        rhs(k,ntyp,15) =  term
        P_tend(k,ntyp) = P_tend(k,ntyp) + term

        term = dicresz*mgchltouMC * pnoice(k)
        rhs(k,14,15) = term
        C_tend(k,2) = C_tend(k,2) + term

        term = docbac(k) * pnoice(k)
        rhs(k,14,14) = term
        C_tend(k,2) = C_tend(k,2) + term
     
        term = tfac(k)*NCrrat*remin(1)*det(k,1)/uMtomgm3 * pnoice(k)
        rhs(k,14,10) = term
        C_tend(k,2) = C_tend(k,2) + term

      enddo  !k=1,kmax

! Phytoplankton components related to growth
      do k = 1,kmax

      if (tirrq(k) .gt. 0.d0)then

!change: March 15, 2010
!       cchlratio = bn*cnratio
!       mgchltouMC = cchlratio/uMtomgm3
        sumdoc = 0.0
        sumutk = 0.0
        sumres = 0.0
!change June 1, 2010
        docexcp = 0.0
        dicresp = 0.0
!endofchange

        do nt = 1,nchl
         !!totgro = gro(k,nt)*obio_P(k,nt+nnut)
         totgro = gro(k,nt)
!!!!     if (nt.eq.4) totgro = gro(k,nt)*pnoice(k)  !for cocco (see ptend)

!change June 1, 2010
          docexcp(nt) = excp*totgro   !phyto production DOC
          dicresp(nt) = resp*totgro   !phyto production DIC
!endofchange

!change: March 15, 2010
!dont account for ice here -- it is incorporated in gro/totgro

            term = - docexcp(nt) 
            rhs(k,nt+nnut,14) = term
            P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

            term = - dicresp(nt)
            rhs(k,nt+nnut,15) = term
            P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

!change June 1, 2010
          term =  bn*(docexcp(nt)+dicresp(nt))
          rhs(k,1,15) = rhs(k,1,15)+term    !accumulate
          P_tend(k,1) = P_tend(k,1) + term

          term = bf*(docexcp(nt)+dicresp(nt))
          rhs(k,4,15) = rhs(k,4,15)+term    !accumulate
          P_tend(k,4) = P_tend(k,4) + term
!endofchange

           sumdoc = sumdoc + docexcp(nt)
          sumutk = sumutk + totgro
         sumres = sumres + dicresp(nt)

        enddo !nt

!change June 1, 2010
         term = bs*(docexcp(1)+dicresp(1))
         rhs(k,3,15) = term
         P_tend(k,3) = P_tend(k,3) + term
!endofchange

        term = sumdoc * mgchltouMC              !phyto prod DOC
        rhs(k,13,5) = term        
        C_tend(k,1) = C_tend(k,1) + term

!!!!    term = ((sumres-sumutk)*mgchltouMC)     !phyto prod DIC
        term = sumres*mgchltouMC                !phyto prod DIC by respiration
        rhs(k,14,5) = term
        C_tend(k,2) = C_tend(k,2) + term

        term = -sumutk*mgchltouMC               !sink DIC by phyto growth
        rhs(k,14,6) = term
        C_tend(k,2) = C_tend(k,2) + term

      endif !tirrq>0
      enddo !k=1,kmax

!ifdef OBIO_RUNOFF
!#ifdef obio_burial 
!     if (p1d(kmax) >= 3700.) then  !for water columns that extend deeper of 3700m
!      do k=1,kmax  
!       if (p1d(k)>3700.) then
!                  !kg/m2 -> kg/m2/s
!         !term = - dic_river_sink/54629.d0/dtsrc   !negative = burial
!         !term = - dic_river_sink/5000.d0/dtsrc   !negative = burial
!         term = - dic_river_sink/2500.d0/dtsrc   !negative = burial
!    .             * 1e6/12.d0/dp1d(k)     !uM/s
!         write(*,'(a,4i5,2e12.4)')'obio_burial: '
!    .           ,nstep,i,j,k,dic_river_sink,term
!total number of points below 3700m 54629
!         rhs(k,4,14) = term
!         C_tend(k,2) = C_tend(k,2) + term
!          
!       endif
!      enddo
!     endif
!#endif
!endif

c pCO2
      if (pco2_online) then
        !this ppco2 routine comes from OCMIP. I am not using psurf
      !and thus not compute dtco2 because these are computed in PBL
      !for the case of gasexch and progn. atmco2, atmco2=dummy
        if (atmco2<0.) then    ! uninitialized
          if (constco2) then
            call get_param("atmCO2", atmco2)    
            if (AM_I_ROOT())print*, 'atmCO2=', atmco2
          else
            atmco2=0.
          endif
        endif
        !note that atmco2 is not being used in compute_pco2_online(so ok if pass 0 here)
        call compute_pco2_online(nstep,i,j,atmco2,
     .            temp1d(1),saln1d(1),car(1,2),alk1d(1),
     .            obio_P(1,1),obio_P(1,3),pnoice(1),
     .            pco2_ij,pHsfc,vrbos)

      else

        call ppco2tab(temp1d(1),saln1d(1),car(1,2),alk1d(1),pCO2_ij)
        if (vrbos) then
           write(*,'(a,3i5,6e12.4)')
     .      'carbon: OFFLINE',nstep,i,j,temp1d(1),saln1d(1),
     .                        car(1,2),alk1d(1),pCO2_ij,pHsfc
        endif
      endif

#ifdef OBIO_mocsy
      patm(1:kmax) = oAPRESS(i,j)
      mocsylat(1:kmax) = oLAT_DG(j,1)
      optCON='mol/m3'
      optT='Tinsitu'
      optP='m'
      call vars(mocsy_ph, mocsy_pco2, mocsy_fco2, mocsy_co2, 
     .                     mocsy_hco3, mocsy_co3, mocsy_OmegaA, 
     .                     mocsy_OmegaC, mocsy_BetaD, mocsy_rhoSW, 
     .                     mocsy_p, mocsy_tempis,     ! OUTPUT
     .                     temp1d, saln1d,alk1d(:)*r3d(:,i,j)*1.d-3, !need alk umolC/kg->mol/m3 as other fields
     .                     car(:,2)*1.d-3, obio_P(:,3)*1.d-3,     !need fields in mol/m3, mmol/m3->mol/m3
     .                     obio_P(:,1)*1.d-3*(1.d0/npratio), 
     .                     patm(:),                          !Patm
     .                     p1d(:), mocsylat(:), kmax,            ! INPUT
     .                     optCON,optT,optP)  
!    .                     optB='l10', optK1K2='m10', optKf='dg')

      if (vrbos) then
         write(*,'(a,3i5,7e12.4)') 
     .    'mocsy output:',nstep,i,j,mocsy_pco2(1), pco2_ij,
     .                              mocsy_fco2(1), co2flux,
     .                              mocsy_co2(1),
     .                              mocsy_OmegaA(1),
     .                              mocsy_OmegaC(1) 
      endif
      OIJL(i,j,1:kmax,ijl_omegaA) = OIJL(i,j,1:kmax,ijl_omegaA) 
     .                            + mocsy_OmegaA
      OIJL(i,j,1:kmax,ijl_omegaC) = OIJL(i,j,1:kmax,ijl_omegaC) 
     .                            + mocsy_OmegaC
#endif

c Update DIC for sea-air flux of CO2

!this is for gas exchange + ocean biology
      if (n_co2n>0) then
        k = 1
        term = co2flux               ! mol/m2/s
!    .     * SECONDS_PER_HOUR        ! mol/m2/hr    !comment out to keep in /s   July 2016
     .     /dp1d(k)                  ! mol/m3/hr
     .     * 1000.D0                 !units of uM/s (=mili-mol/m3/s)
                                     !do not mulitply by pnoice here, 
                                     !this is done in SURFACE.f (ptype)
        rhs(k,14,16) = term
        C_tend(k,2) = C_tend(k,2) + term

        if (vrbos) then
          write(*,'(a,3i7,3e12.4)')
     .      'obio_carbon (coupled):',
     .      nstep,i,j,dp1d(1),co2flux,term     !this flux should be mol,co2/m2/s
        endif

      !abiotic DIC tracer
      if (n_abioDIC.ne.0) then
!uses CMIP6 values
        k = 1
        Ts = temp1d(k)
        scco2 = 2116.8 - 136.25*Ts + 4.7353*Ts**2 - 0.092307*Ts**3 + 0.000755*Ts**4
        wssq = wind*wind
        if (scco2.lt.0.) then
          scco2arg=1.d-10
          rkwco2=1.d-10
        else
          scco2arg = (scco2/660.D0)**(-0.5)      !Schmidt number
          rkwco2 = awan*wssq*scco2arg           !transfer coeff (units of m/s)
        endif
        tk = tf+Ts
        tk100 = tk*0.01
        tk1002 = tk100*tk100
        ff = exp(-160.7333 + 215.4152/tk100  +       !solub in mol/kg/picoatm
     .         89.8920*log(tk100) - 1.47759*tk1002 +
     .         saln1d(k) * (0.029941 - 0.027455*tk100 +
     .         0.0053407*tk1002))
        xco2 = atmCO2*1013.D0/stdslp
        sdic_uM=SDIC/trmo_unit_factor(1,14)  !14->DIC
        call compute_pco2_online(nstep,i,j,atmco2,
     .            temp1d(1),saln1d(1),sdic_uM,alk1d(1),
     .            obio_P(1,1),obio_P(1,3),pnoice(1),
     .            pCO2_abio,dummy,vrbos)

        deltco2 = (xco2-pCO2_abio)*ff*rho1d(k)*1d-6 !convert ff mol/m3/uatm
        flxmolm3 = (rkwco2*deltco2/dp1d(k))   !units of mol/m3/s
        term = flxmolm3*1000.D0*pnoice(k)    !units of uM/s (=mili-mol/m^3/s)

          !trmo(i,j,1,n_abioDIC) = trmo(i,j,1,n_abioDIC)
          SDIC = SDIC 
     .                          + term*DTS**1e-6*12.d0     !term is in mili-mol/m3/s -> trmo is in kg,C
     &                          *ddxypo*dp1d(1)
!     .                          * dxypo(j)*dp1d(1)
!    .                          * mo(i,j,1)*dxypo(j)/rho_water 
!    .                          * mo(i,j,1)*dxypo(j)/1024.d0
      endif

      else

!this is for only ocean biology but no gas exchange: 
      !when ocean biology but no CO2 gas exch
      !atmco2 is set to constant
        k = 1
        Ts = temp1d(k)
!       scco2 = 2073.1 - 125.62*Ts + 3.6276*Ts**2 - 0.043219*Ts**3
        !new OCMIP2016 values
        scco2 = 2116.8 - 136.25*Ts + 4.7353*Ts**2 - 0.092307*Ts**3 + 0.000755*Ts**4
        wssq = wind*wind
        if (scco2.lt.0.) then
          scco2arg=1.d-10
          rkwco2=1.d-10
        else
          scco2arg = (scco2/660.D0)**(-0.5)      !Schmidt number
          rkwco2 = awan*wssq*scco2arg           !transfer coeff (units of m/s)
        endif
        tk = tf+Ts
        tk100 = tk*0.01
        tk1002 = tk100*tk100
!       ff = exp(-162.8301 + 218.2968/tk100  +       !solub in mol/kg/picoatm
!    .         90.9241*log(tk100) - 1.47696*tk1002 +
!    .         saln1d(k) * (.025695 - .025225*tk100 +
!    .         0.0049867*tk1002))
        !new OCMIP2016 values
        ff = exp(-160.7333 + 215.4152/tk100  +       !solub in mol/kg/picoatm
     .         89.8920*log(tk100) - 1.47759*tk1002 +
     .         saln1d(k) * (0.029941 - 0.027455*tk100 +
     .         0.0053407*tk1002))


        xco2 = atmCO2*1013.D0/stdslp
        deltco2 = (xco2-pCO2_ij)*ff*rho1d(k)*1d-6 !convert ff mol/m3/uatm
        flxmolm3 = (rkwco2*deltco2/dp1d(k))   !units of mol/m3/s
!       flxmolm3h = flxmolm3*SECONDS_PER_HOUR !units of mol/m3/hr       July 2016
        term = flxmolm3*1000.D0*pnoice(k)    !units of uM/s (=mili-mol/m^3/s)
        rhs(k,14,16) = term
        C_tend(k,2) = C_tend(k,2) + term

      !flux sign is (atmos-ocean)>0, i.e. positive flux is INTO the ocean
        co2flux= rkwco2*(xco2-pCO2_ij)*ff*1.0245D-3*pnoice(k)! air-sea co2 flux
     .            *SECONDS_PER_HOUR                             ! mol/m2/hr
     .            *44.d0*HOURS_PER_DAY*DAYS_PER_YEAR            ! grC/m2/yr
        if (vrbos) then
            write(*,'(a,3i5,11e12.4)')'obio_carbon, fluxdiag:',
     .      nstep,i,j,dp1d(k),Ts,saln1d(k),scco2arg,wssq,rkwco2,
     .      xco2,pCO2_ij,ff,flxmolm3,co2flux
        endif

        if (vrbos) then
          write(6,'(a,3i7,9e12.4)')'obio_carbon(watson):',
     .      nstep,i,j,Ts,scco2arg,wssq,rkwco2,ff,xco2,pCO2_ij,
     .      rkwco2*(xco2-pCO2_ij)*ff*1.0245D-3,term     !this flux should have units mol,co2/m2/s
        endif

      !abiotic DIC tracer
      if (n_abioDIC.ne.0) then
         ! trmo(i,j,1,n_abioDIC) = trmo(i,j,1,n_abioDIC)
          SDIC = SDIC 
     .                          + term*DTS**1e-6*12.d0     !term is in mili-mol/m3/s -> trmo is in kg,C
!    .                          * mo(i,j,1)*dxypo(j)/rho_water 
!     .                          * mo(i,j,1)*dxypo(j)/1024.d0
     &                          *mmo*ddxypo/1024.d0
      endif
      endif

      return
      end subroutine obio_carbon

c ---------------------------------------------------------------------------
      subroutine compute_pco2_online(nstep,i,j,atmco2,
     .            T,S,dic,alk,nitr,sili,pnoice,
     .            pco2,pH,vrbos)

      implicit none

      integer nstep,i,j
      real, intent(in) :: T, S, dic, nitr, sili, pnoice, atmco2
      real, intent(inout):: pco2,pH,alk
      logical vrbos

        call ppco2(T,S,dic,alk,nitr,sili,atmCO2,pCO2,pH)

!note: pco2 is computed as if it is 100% open ocean cell. This is why
!in the flux computation below we need to take into account pnoice
!also in the diagnostics

!     !limits on pco2 ---more work needed
! ppco2 does not handle well the extreme salinity cases, such as
! when ice melts/forms, in river outflows.
        if (S.ge.40. .and. pCO2.lt.100.)pCO2=100.
        if (S.le.31. .and. pCO2.gt.800.)pCO2=800.
        if (pCO2 .lt. 100.) pCO2=100.
!        if (pCO2 .gt.1000.) pCO2=1000. @PL removed Mar 05, 2020

        if(vrbos)then
          write(*,'(a,i8,2i5,9e12.4)')
     .      'carbon: ONLINE',
     .      nstep,i,j,T,S,dic,alk,nitr,sili,pCO2,pH,pnoice
        endif

        return
        end subroutine compute_pco2_online

c ---------------------------------------------------------------------------
      subroutine ppco2tab(T,S,car1D,TA,pco21D)
 
c  Computes pCO2 in the surface layer and delta pCO2 with the
c  atmosphere using OCMIP protocols.  Uses pre-computed look-up
c  table to increase computational speed.
c  Var  min     max     increment
c  T0   -2      37.5    0.5
c  sal  30      39.5    0.5
c  DIC  1800    2450    2
c  TA   2000    2500    2
c

      USE obio_dim, only: ALK_CLIM
      use domain_decomp_1d, only: am_i_root
      use filemanager, only: openunit, closeunit

      implicit none

      real, intent(in) :: T, S, car1D
      real, intent(inout):: TA
      real, intent(out):: pco21D
      integer  :: it0,isal,idic,ita

      real, allocatable, dimension(:,:,:,:), save :: pco2tab
!save from moved ifst parts
      integer, parameter :: it0inc=1,nt0=80/it0inc,nsal=20
      integer, parameter :: idicinc=2,ndic=(650+idicinc)/idicinc
      integer, parameter :: itainc=2,nta=(500+itainc)/itainc
      real, parameter :: tabar=2310.0 !mean total alkalinity uE/kg; OCMIP
      real, parameter :: Sbar=34.836  !global mean annual salinity (area-wt)
      real :: Tsfc,sal,DIC
      integer :: iu_bio, i, j, k, nl

      if (.not.allocated(pco2tab)) then
        allocate(pco2tab(nt0, nsal, ndic, nta))
        call openunit('pco2table',iu_bio)
        if (AM_I_ROOT()) then
          print*, '    '
          print*, 'obio_init, pco2tbl: ',nta,ndic,nsal,nt0
        endif
        do nl=1,nta
          do k=1,ndic
            do j=1,nsal
              do i=1,nt0
                read(iu_bio,'(e12.4)')pco2tab(i,j,k,nl)
              enddo
            enddo
          enddo
        enddo
        call closeunit(iu_bio)
        if (AM_I_ROOT()) then
          print*,'BIO: read pCO2 table: ',
     .        pco2tab(1,1,1,1),pco2tab(50,10,100,100)
          print*, '    '
        endif
      endif
c Get pco2
      pco21D = 0.0
       Tsfc = T
        sal = S
        DIC = car1D     !uM
        if (ALK_CLIM.eq.0) TA = tabar*S/Sbar  !adjust alk for salinity
        it0 = nint((Tsfc+2.0)*2.0)/it0inc + 1
       isal = nint((sal-30.0)*2.0) + 1
       idic = (nint(DIC-1800.0))/idicinc + 1
        ita = (nint(TA-2000.0))/itainc + 1

cdiag  write(*,'(a,4f9.3,4i6)')
cdiag.         'ppco2: ',Tsfc,sal,DIC,TA,it0,isal,idic,ita

       !test within bounds
       if ( it0.gt.nt0)  then
           it0= min( it0,nt0 )
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, it0=',it0,nt0,Tsfc
       endif
       if (isal.gt.nsal) then
           isal= min(isal,nsal)
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, isal=',isal,nsal,sal
       endif
       if (idic.gt.ndic) then
           idic= min(idic,ndic)
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, idic=',idic,ndic,DIC
       endif
       if ( ita.gt.nta) then
            ita = min(ita,nta)
c          write(*,'(a,2i5,f9.3)')
c    .             'correction in ppco2tab, ita =',ita,nta,TA
       endif

       if (it0.lt.  1) ita = max(it0,1)
       if (isal.lt. 1) isal= max(isal,1)
       if (idic.lt. 1) idic= max(idic,1)
       if (ita.lt.  1) ita = max(ita,1)

       pco21D = pco2tab(it0,isal,idic,ita)

      return
      end

c ---------------------------------------------------------------------------
      subroutine ppco2(T,S,car1D,TA,nitr,silic,atmCO2,
     .                 pco2,pHsfc)
c
c  Computes pCO2 in the surface layer and delta pCO2 with the 
c  atmosphere using OCMIP protocols.
c
      USE obio_dim, only: ALK_CLIM

      implicit none

      real, parameter :: phmin=7.5, phmax=8.6 ! min/max pH for iteration
      real*8, parameter :: tabar=2310.0D0    !mean total alkalinity uE/kg; OCMIP
      real*8, parameter :: stdslp=1013.25D0  !standard sea level pressure in mb
      real*8, parameter :: Sbar=34.836D0     !global mean annual salinity (area-wt)

      real*8 dtco2,atmCO2
      real*8 T,S,car1D,TA,nitr,silic,pco2
      real*8 phlo,phhi,ph,atmpres,pHsfc
      real*8 Tsfc,sal,DIC,PO4,Si,dic_in,ta_in,pt_in,sit_in,xco2_in
      real*8 co2star,pco2surf,dpco2

c
c  Constants
      dtco2 = 0.0
      pco2 = 0.0
      phlo = pHmin
      phhi = pHmax
      phlo = 1.D0   !range for pH
      phhi =16.D0   !range for pH
      ph = 8.0
c
c  Use OCMIP subroutines
       Tsfc = T
       sal = S
       DIC = car1D             !uM
       PO4 = nitr*0.1          !uM phosphate, converted from NO3/PO4
c                               ratio from Conkright et al 1994, 
c                               global using 1st 3 standard depths 
c                               (0, 10, and 20m)
        Si= silic              !uM Si
c
c  Convert to units for co2calc
       dic_in = dic*1.0E-3  !uM to mol/m3
       if (ALK_CLIM.eq.0) TA = tabar*S/Sbar  !adjust alk for salinity
       ta_in = ta*1024.5d0*1.0E-6  !uE/kg to E/m3
       pt_in = PO4*1.0E-3   !uM to mol/m3
       sit_in = Si*1.0E-3   !uM to mol/m3
       xco2_in = atmco2
       !!!atmpres = slp/stdslp  --not done here; do it in PBL
c
!      print*,Tsfc,sal,dic_in,ta_in,pt_in,sit_in,
!    *      phlo,phhi,ph,xco2_in

       call co2calc_SWS(Tsfc,sal,dic_in,ta_in,pt_in,sit_in,
     *      phlo,phhi,ph,xco2_in,co2star,pco2surf)
        pco2 = pco2surf
       !!!dtco2 = dco2star
       pHsfc = pH
c
!      print*, 'pco2=',pco2

      return
      end


C
C ---------------------------------------------------------------------
C 
      subroutine co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
     &                  ,phlo,phhi,ph,xco2_in
     &                  ,co2star,pCO2surf)
      USE CONSTANT, only: tf
C
C-------------------------------------------------------------------------
C
C Modified from co2calc.f (RCS version 1.8, OCMIP-2) 
C - by A. Mouchet, 2004:
C
C NOW; "All" constants are given on seawater H scale (hSWS) 
C - HOWEVER, dissociation constants for S and F are on the 'free' H scale 
C            (this is necessary to work with hSWS)
C - this routine corrects for inconsistencies in the previous version.
C
C - Other changes:
C   * use ta_iter_SWS instead of ta_iter_1;
C   * hSWS replaces htotal;
C   * moved upward the calculation of bt, st and ft 
C     (needed in evaluation of kb);
C   * added various comments and references.
C
C
C SUBROUTINE CO2CALC_SWS
C
C PURPOSE
C	Calculate delta co2* from total alkalinity and total CO2 at
C temperature (t), salinity (s) and "atmpres" atmosphere total pressure. 
C
C USAGE
C       call co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
C    &                  ,phlo,phhi,ph,xco2_in,atmpres
C    &                  ,co2star,dco2star,pCO2surf,dpco2)
C
C INPUT
C	dic_in = total inorganic carbon (mol/m^3) 
C                where 1 T = 1 metric ton = 1000 kg
C	ta_in  = total alkalinity (eq/m^3) 
C	pt_in  = inorganic phosphate (mol/m^3) 
C	sit_in = inorganic silicate (mol/m^3) 
C	t      = temperature (degrees C)
C	s      = salinity (PSU)
C	phlo   = lower limit of pH range
C	phhi   = upper limit of pH range
C	xco2_in=atmospheric mole fraction CO2 in dry air (ppmv) 
C	atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
C
C       Note: arguments dic_in, ta_in, pt_in, sit_in, and xco2_in are 
C             used to initialize variables dic, ta, pt, sit, and xco2.
C             * Variables dic, ta, pt, and sit are in the common block 
C               "species".
C             * Variable xco2 is a local variable.
C             * Variables with "_in" suffix have different units 
C               than those without.

C OUTPUT
C	co2star  = CO2*water (mol/m^3)
C	dco2star = delta CO2 (mol/m^3)
c       pco2surf = oceanic pCO2 (ppmv)
c       dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)
C
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
c     - Models carry tracers in mol/m^3 (on a per volume basis)
c     - Conversely, this routine, which was written by observationalists 
c       (C. Sabine and R. Key), passes input arguments in umol/kg  
c       (i.e., on a per mass basis)
c     - I have changed things slightly so that input arguments are in mol/m^3,
c     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in) 
c       should be given in mol/m^3; output arguments "co2star" and "dco2star"  
c       are likewise in mol/m^3.

C FILES and PROGRAMS NEEDED
C	drtsafe
C	ta_iter_SWS
C
C--------------------------------------------------------------------------
C
        implicit none

        real*8 t,s,dic_in,ta_in,pt_in,sit_in,phlo,phhi,
     .         ph,xco2_in,co2star,pCO2surf
     
        real*8 st,ft,ff,x1,x2,xacc,hSWS,hSWS2,bt,scl,s15,s2,sqrtis,
     .         tk,dic,permeg,xco2,tk100,tk1002,sqrts,permil,pt,sit,ta,
     .         dlogtk

        real*8 drtsafe

        real*8 invtk,is,is2
        real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
        common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
        common /species/bt,st,ft,sit,pt,dic,ta
        external ta_iter_SWS
C

c       ---------------------------------------------------------------------
C       Change units from the input of mol/m^3 -> mol/kg:
c       (1 mol/m^3)  x (1 m^3/1024.5 kg)
c       where the ocean's mean surface density is 1024.5 kg/m^3
c       Note: mol/kg are actually what the body of this routine uses 
c       for calculations.  
c       ---------------------------------------------------------------------
        permil = 1.d0 / 1024.5d0
c       To convert input in mol/m^3 -> mol/kg 

!       print*,permil,pt_in,sit_in,ta_in,dic_in

        pt=pt_in*permil
        sit=sit_in*permil
        ta=ta_in*permil
        dic=dic_in*permil
      
!       print*,pt,sit,ta,dic

c       ---------------------------------------------------------------------
C       Change units from uatm to atm. That is, atm is what the body of 
c       this routine uses for calculations.
c       ---------------------------------------------------------------------
        permeg=1.e-6
c       To convert input in uatm -> atm
        xco2=xco2_in*permeg
c       ---------------------------------------------------------------------
C
C*************************************************************************
C Calculate all constants needed to convert between various measured
C carbon species. References for each equation are noted in the code. 
C Once calculated, the constants are
C stored and passed in the common block "const". The original version of this
C code was based on the code by Dickson in Version 2 of "Handbook of Methods
C for the Analysis of the Various Parameters of the Carbon Dioxide System
C in Seawater", DOE, 1994 (SOP No. 3, p25-26). 
C
C Derive simple terms used more than once
C
      tk = tf + t
      tk100 = tk/100.0
      tk1002=tk100*tk100
      invtk=1.0/tk
      dlogtk=log(tk)
      is=19.924*s/(1000.-1.005*s)
      is2=is*is
      sqrtis=sqrt(is)
      s2=s*s
      sqrts=sqrt(s)
      s15=s**1.5
      scl=s/1.80655
C
C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C
C Uppstrom (1974)
      bt = 0.000232 * scl/10.811
C Morris & Riley (1966)
      st = 0.14 * scl/96.062
C Riley (1965)
      ft = 0.000067 * scl/18.9984
C
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
C
      ff = exp(-162.8301 + 218.2968/tk100  +
     & 90.9241*log(tk100) - 1.47696*tk1002 +
     & s * (.025695 - .025225*tk100 + 
     & 0.0049867*tk1002))
C
C K0 from Weiss 1974
C
      k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100) +
     & s * (.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))

C
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]     on hSWS
C
C Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
C (Original reference: Dickson and Millero, DSR, 1987)
C
      k1=10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
     & 0.0118 * s + 0.000116*s2))
C
      k2=10**(-1*(1394.7*invtk + 4.777 - 
     & 0.0184*s + 0.000118*s2))
C
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4] on hSWS
C
C Millero p.670 (1995)
C
      k1p = exp(-4576.752*invtk + 115.540 - 18.453 * dlogtk +
     & (-106.736*invtk + 0.69171) * sqrts +
     & (-0.65643*invtk - 0.01844) * s)
C
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4] on hSWS
C
C Millero p.670 (1995)
C
      k2p = exp(-8814.715*invtk + 172.1033 - 27.927 * dlogtk +
     & (-160.340*invtk + 1.3566) * sqrts +
     & (0.37335*invtk - 0.05778) * s)
C
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4] on hSWS
C
C Millero p.670 (1995)
C
      k3p = exp(-3070.75*invtk - 18.126 + 
     & (17.27039*invtk + 2.81197) *
     & sqrts + (-44.99486*invtk - 0.09984) * s)
C
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
C
C Millero p.671 (1995) using data from Yao and Millero (1995)
C change to (mol/ kg soln)
C
      ksi = exp(-8904.2*invtk + 117.400 - 19.334 * dlogtk +
     & (-458.79*invtk + 3.5913) * sqrtis +
     & (188.74*invtk - 1.5998) * is +
     & (-12.1652*invtk + 0.07871) * is2 +
     & log(1.0-0.001005*s))
C
C------------------------------------------------------------------------
C kw = [H][OH] on hSWS
C
C Millero p.670 (1995) using composite data
C
      kw = exp(-13847.26*invtk + 148.9802 - 23.6521 * dlogtk +
     & (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     & sqrts - 0.01615 * s)
C
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4] on free H scale
C
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C change to (mol/ kg soln)
C
      ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     & (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
     & (35474*invtk - 771.54 + 114.723*dlogtk) * is -
     & 2698*invtk*is**1.5 + 1776*invtk*is2 +
     & log(1.0 - 0.001005*s))
C
C------------------------------------------------------------------------
C kf = [H][F]/[HF] on free H scale
C
C Dickson and Riley (1979)
C change to (mol/ kg soln)
C
      kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     & log(1.0 - 0.001005*s)) 
C
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2] on hSWS
C
C Dickson p.673 (1990)
C change from htotal to hSWS
C
      kb=exp( (-8966.90 - 2890.53*sqrts - 77.942*s +
     & 1.728*s15 - 0.0996*s2)*invtk +
     & (148.0248 + 137.1942*sqrts + 1.62142*s) +
     & (-24.4344 - 25.085*sqrts - 0.2474*s) *
     & dlogtk + 0.053105*sqrts*tk +
     & log((1+(st/ks)+(ft/kf))/(1+(st/ks))) )
C
C*************************************************************************
C
C Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
C The solution converges to err of xacc. The solution must be within
C the range x1 to x2.
C
C If DIC and TA are known then either a root finding or iterative method
C must be used to calculate hSWS. In this case we use the Newton-Raphson
C "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
C error trapping removed).
C
C As currently set, this procedure iterates about 12 times. The x1 and x2
C values set below will accomodate ANY oceanographic values. If an initial
C guess of the pH is known, then the number of iterations can be reduced to
C about 5 by narrowing the gap between x1 and x2. It is recommended that
C the first few time steps be run with x1 and x2 set as below. After that,
C set x1 and x2 to the previous value of the pH +/- ~0.5. The current
C setting of xacc will result in co2star accurate to 3 significant figures
C (xx.y). Making xacc bigger will result in faster convergence also, but this
C is not recommended (xacc of 10**-9 drops precision to 2 significant figures).
C
C Parentheses added around negative exponents (Keith Lindsay)
C
      x1 = 10.0**(-phhi)
      x2 = 10.0**(-phlo)
c	xacc = 1.e-10
      xacc = 1.e-14
      hSWS = drtsafe(ta_iter_SWS,x1,x2,xacc)
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
      hSWS2=hSWS*hSWS
      co2star=dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
      !!!co2starair=xco2*ff*atmpres  !not needed here
      !!!dco2star=co2starair-co2star
      ph=-log10(hSWS)

c
c       ---------------------------------------------------------------
cc      Add two output arguments for storing pCO2surf
cc      Should we be using K0 or ff for the solubility here?
c       ---------------------------------------------------------------
        pCO2surf = co2star / ff
        !!!dpCO2    = pCO2surf - xco2*atmpres
C
C  Convert units of output arguments
c      Note: co2star and dco2star are calculated in mol/kg within this routine 
c      Thus Convert now from mol/kg -> mol/m^3
       co2star  = co2star / permil
       !!!dco2star = dco2star / permil

c      Note: pCO2surf and dpCO2 are calculated in atm above. 
c      Thus convert now to uatm
       pCO2surf = pCO2surf / permeg
       !!!dpCO2    = dpCO2 / permeg
C
      return
      end

C
C  ---------------------------------------------------------------------
C 
        subroutine ta_iter_SWS(x,fn,df)
#ifdef TRACERS_Alkalinity
        use obio_com, only: co3_conc
#endif
        implicit none

        real*8 x,fn,df,b2,db,dic,bt,pt,sit,ta,x3,x2,c,a,a2,da,b,
     .         st,ft,ff,hSWS

        real*8 k12,k12p,k123p
        real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
        common /const/k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi,ff,hSWS
        common /species/bt,st,ft,sit,pt,dic,ta
C
C Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
C - by A. Mouchet, 2004:
C Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
C  1) fixed errors in signs, parenthesis and coefficient c in derivative
C  2) changed from Total to Seawater Scale 
C     * c defined for seawater H scale; 
C     * fn and df adapted to KF on free H scale
C     * comments have been adapted
C

C
C This routine expresses TA as a function of DIC, hSWS and constants.
C It also calculates the derivative of this function with respect to 
C hSWS. It is used in the iterative solution for hSWS. In the call
C "x" is the input value for hSWS, "fn" is the calculated value for TA
C and "df" is the value for dTA/dhSWS
C
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0 + st/ks + ft/kf
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0*x2 + 2.0*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0*x + k1
!!#ifdef TOPAZ_params
#ifdef TRACERS_Alkalinity
!     print*,'ta_iter_SWS: co3_conc',dic,k12,b
      co3_conc = 2.0*dic*k12/b
#endif
!!#endif
C
C	fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
C===========================================================================
C
      fn = k1*x*dic/b +
     &       2.0*dic*k12/b +
     &       bt/(1.0 + x/kb) +
     &       kw/x +
     &       pt*k12p*x/a +
     &       2.0*pt*k123p/a +
     &       sit/(1.0 + x/ksi) -
     &       x/c -
     &       st/(1.0 + ks/(x/c)) -
     &       ft/(1.0 + kf/(x/c)) -
     &       pt*x3/a -
     &       ta
C
C	df = dfn/dx
C
       df = ((k1*dic*b) - k1*x*dic*db)/b2 -
     &       2.0*dic*k12*db/b2 -
     &       bt/kb/(1.0+x/kb)**2. -
     &       kw/x2 +
     &       (pt*k12p*(a - x*da))/a2 -
     &       2.0*pt*k123p*da/a2 -
     &       sit/ksi/(1.0+x/ksi)**2. -
     &       1.0/c -
     &       st *(1.0 + ks/(x/c))**(-2.0) *(ks*c/x2) -
     &       ft*(1.0 + kf/(x/c))**(-2.0) *(kf*c/x2) -
     &       pt*x2*(3.0*a-x*da)/a2
C
      return
      end

C
C ---------------------------------------------------------------------
C 
      REAL*8 FUNCTION DRTSAFE(FUNCD,X1,X2,XACC)
C
C	File taken from Numerical Recipes. Modified  R.M.Key 4/94
C
      implicit none

      integer j,maxit
      real*8 xacc,dx,temp,f,df,dxold,swap,xh,xl,x1,fl,x2,fh

      MAXIT=100
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      IF(FL .LT. 0.0) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE=.5*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(DRTSAFE,F,DF)
      DO 100, J=1,MAXIT
        IF(((DRTSAFE-XH)*DF-F)*((DRTSAFE-XL)*DF-F) .GE. 0.0 .OR.
     &     ABS(2.0*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          DRTSAFE=XL+DX
          IF(XL .EQ. DRTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE
          DRTSAFE=DRTSAFE-DX
          IF(TEMP .EQ. DRTSAFE)RETURN
        END IF
        IF(ABS(DX) .LT. XACC)RETURN
        CALL FUNCD(DRTSAFE,F,DF)
        IF(F .LT. 0.0) THEN
          XL=DRTSAFE
          FL=F
        ELSE
          XH=DRTSAFE
          FH=F
        END IF
  100  CONTINUE
      RETURN
      END
