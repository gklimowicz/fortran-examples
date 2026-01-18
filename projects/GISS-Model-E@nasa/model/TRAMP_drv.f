#include "rundeck_opts.h"
      MODULE AMP_AEROSOL
!@sum Driver for Aerosol Microphysics
!@auth Susanne Bauer
      USE AERO_CONFIG, ONLY: NMODES
      USE AERO_PARAM,  ONLY: NEMIS_SPCS
      USE RESOLUTION,   ONLY: LM
      IMPLICIT NONE
      SAVE

C**************  Latitude-Dependant (allocatable) *******************
      ! Mie lookup tables
      REAL*8, DIMENSION(15,17,23,6)      :: AMP_EXT, AMP_ASY, AMP_SCA   !(15,17,23,6) (RE,IM,size,lambda)
      REAL*8, DIMENSION(15,17,23)        :: AMP_Q55
      REAL*8, DIMENSION(23,26,26,26,6)   :: AMP_EXT_CS,AMP_ASY_CS,AMP_SCA_CS !(23,26,26,26,6) (radius,OC,SO4,H2O,lambda)
      REAL*8, DIMENSION(23,26,26,26)     :: AMP_Q55_CS
      ! 1 Dim arrays for Radiation
      REAL*8, DIMENSION(LM,nmodes)       :: Reff_LEV, NUMB_LEV
      REAL*8, DIMENSION(LM,nmodes)       :: MIX_OC, MIX_SU, MIX_AQ
      COMPLEX*8, DIMENSION(LM,nmodes,6)  :: RindexAMP
      REAL*8, DIMENSION(LM,nmodes,7)     :: dry_Vf_LEV
      ! FALSE : one Radiation call
      ! TRUE  : nmodes Radiation calls
      INTEGER                            :: AMP_RAD_KEY = 1 ! 1=Volume Mixing || 2=Core - Shell || 3=Maxwell Garnett

      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)       :: AQsulfRATE !(i,j,l)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: DIAM       ![m](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: DIAM_dry       ![m](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: AMP_dens   !density(i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: AMP_TR_MM  !molec. mass(i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: NACTV      != 1.0D-30  ![#/m^3](i,j,l,nmodes)
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:)     :: VDDEP_AERO != 1.0D-30  ![m/s](i,j,nmodes,2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)         :: ampPM2p5, ampPM10  ! [kg/kg air]

!-------------------------------------------------------------------------------------------------------------------------
!     The array VDDEP_AERO(X,Y,Z,I,1) contains current values for the dry deposition velocities 
!     for aerosol number concentrations for mode I. 
!     The array VDDEP_AERO(X,Y,Z,I,2) contains current values for the dry deposition velocities 
!     for aerosol mass   concentrations for mode I. 
!     Values in VDDEP_AERO are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------
!     The array NACTV(X,Y,Z,I) contains current values of the number of aerosol particles 
!     activated in clouds for each mode I for use outside of the MATRIX microphysical module.
!     Values in NACTV are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------      
!     1 - BC  2-BCmix 3-OC 4-OCmix 5-SS1  6-SS2 7-D1 8-D2
!-------------------------------------------------------------------------------------------------------------------------      
!-------------------------------------------------------------------------------------------------------------------------
!     The array DIAM(x,y,z) contains current values of some measure of average ambient mode diameter for each mode  
!       for use outside of the MATRIX microphysical module where it is calculated. 
!       Values in DIAM are saved at the top of the subr. MATRIX before microphysical evolution 
!       for the current time step is done. 
!
!     The current measure of particle diameter is the diameter of average mass:
! 
!        DIAM(x,y,z) = [ (6/pi) * (Mi/Ni) * (1/D) ]^(1/3)
!
!     with Mi the total mass concentration (including water) in mode i, Ni the number concentration in mode i, and
!     D a constant ambient particle density, currently set to D = DENSP = 1.4 g/cm^3. 
!-------------------------------------------------------------------------------------------------------------------------      
      END MODULE AMP_AEROSOL

      SUBROUTINE MATRIX_DRV
!@vers 2013/03/27
      USE AmpTracersMetadata_mod, only: AMP_MODES_MAP, AMP_NUMB_MAP,
     *  AMP_AERO_MAP
      USE TRACER_COM
!      USE TRACER_COM, only: n_H2SO4, n_M_ACC_SU, n_M_AKK_SU, n_M_BC1_BC,
!     *  n_M_DD1_DU, n_M_DD2_DU, n_M_OCC_OC, n_M_SSA_SS, n_M_SSC_SS,
!     *  n_NH3, nBiomass,nAircraft, nChemistry, nOther, ntmAMPe, nVolcanic, trm, ntmAMPi 
!#ifdef  TRACERS_SPECIAL_Shindell
!      USE TRACER_COM, only: n_HNO3
!#endif
      use OldTracer_mod, only: trname
      USE TRDIAG_COM, only : taijs=>taijs_loc,taijls=>taijls_loc
     *     ,ijts_AMPp,ijlt_AMPm,ijts_AMPpdf, ijts_AMPe
     *     ,itcon_AMP,itcon_AMPm
      USE AMP_AEROSOL
      USE AEROSOL_SOURCES, only: off_HNO3

      USE RESOLUTION, only : im,jm,lm     ! dimensions
      USE ATM_COM, only : 
     $                      t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
      USE MODEL_COM, only : dtsrc
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE CONSTANT,   only:  lhe,mair,gasc,rgas  
      USE FLUXES, only: tr3Dsource,trsource,trflux1
      USE ATM_COM,   only: pmid,pk,byMA,gz, MA   ! midpoint pressure in hPa (mb)
!                                           and pk is t mess up factor
!                                           byMA  1/Air mass (m^2/kg)
      USE AERO_CONFIG
      USE AERO_INIT
      USE AERO_PARAM, only: IXXX, IYYY, ILAY, NEMIS_SPCS
      USE AERO_SETUP 
      USE PBLCOM,     only: EGCM !(LM,IM,JM) 3-D turbulent kinetic energy [m^2/s^2]
      USE DOMAIN_DECOMP_ATM,only: GRID, getDomainBounds, am_i_root
#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups,
     &                      inc_subdd,find_groups
#endif  /* CACHED_SUBDD */

      IMPLICIT NONE

      REAL(8):: TK,RH,PRES,TSTEP,AQSO4RATE,PM(3)
      REAL(8):: AERO(NAEROBOX)     ! aerosol conc. [ug/m^3] or [#/m^3]
      REAL(8):: GAS(NGASES)        ! gas-phase conc. [ug/m^3]
      REAL(8):: EMIS_MASS(NEMIS_SPCS) ! mass emission rates [ug/m^3]
      REAL(8):: SPCMASS(NMASS_SPCS+2)
      REAL(8):: DT_AERO(NDIAG_AERO,NAEROBOX) !NDIAG_AERO=15
      REAL(8):: yS, yM, ZHEIGHT1,WUP,AVOL 
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      INTEGER:: j,l,i,n,J_0, J_1, I_0, I_1, m,nAMP
C**** functions
      REAL(8):: QSAT

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups),k
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     sddarr3d
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     sddarr2d

#endif  /* CACHED_SUBDD */

      call getDomainBounds(grid, J_STRT =J_0, J_STOP =J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

#ifndef  TRACERS_SPECIAL_Shindell
      CALL READ_OFFHNO3(OFF_HNO3)
#endif

      NACTV(I_0:I_1,J_0:J_1,:,:)      = 0.d0 
      VDDEP_AERO(I_0:I_1,J_0:J_1,:,:) = 0.d0 
      DIAM(I_0:I_1,J_0:J_1,:,:)       = 0.d0
      DIAM_dry(I_0:I_1,J_0:J_1,:,:)   = 0.d0
      AMP_dens(I_0:I_1,J_0:J_1,:,:)   = 0.d0
      AMP_TR_MM(I_0:I_1,J_0:J_1,:,:)  = 0.d0

      DO L=1,LM                            
      DO J=J_0,J_1                          
      DO I=I_0,I_1                 

      IXXX = I
      IYYY = J
      ILAY = L
      DT_AERO(:,:) = 0.d0
      EMIS_MASS(:) = 0.d0
      AERO(:)      = 0.d0
! meteo
      TK = pk(l,i,j)*t(i,j,l)           !should be in [K]
      RH = MIN(1.d0,q(i,j,l)/QSAT(TK,lhe,pmid(l,i,j))) ! rH [0-1]
      PRES= pmid(l,i,j)*100.                  ! pmid in [hPa]
      TSTEP=dtsrc
      ZHEIGHT1 = GZ(i,j,l) /1000./9.81
      WUP = SQRT(.6666667*EGCM(l,i,j))  ! updraft velocity

c avol [m3/gb] mass of air pro m3      
      AVOL = MA(l,i,j)*axyp(i,j)/mair*1000.d0*gasc*tk/pres 
! in-cloud SO4 production rate [ug/m^3/s] ::: AQsulfRATE [kg] 
      AQSO4RATE = AQsulfRATE (i,j,l)* 1.d9  / AVOL /dtsrc
c conversion trm [kg/gb] -> [ug /m^3]
      GAS(GAS_H2SO4) = trm(i,j,l,n_H2SO4)* 1.d9 / AVOL! [ug H2SO4/m^3]
c conversion trm [kg/kg] -> [ug /m^3]
#ifdef  TRACERS_SPECIAL_Shindell
      GAS(GAS_HNO3) = trm(i,j,l,n_HNO3)*1.d9 / AVOL!   [ug HNO3/m^3]
#else
      GAS(GAS_HNO3) = off_HNO3(i,j,l)*1.d9 /AVOL !   [ug HNO3/m^3]
#endif
c conversion trm [kg/gb] -> [ug /m^3]
      GAS(GAS_NH3) = trm(i,j,l,n_NH3)* 1.d9 / AVOL!   [ug NH3 /m^3]
!  [kg/s] -> [ug/m3/s]

       DO n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
c conversion trm [kg/gb] -> AERO [ug/m3]
         if(AMP_NUMB_MAP(nAMP).eq. 0) then
       AERO(AMP_AERO_MAP(nAMP)) =trm(i,j,l,n)*1.d9 / AVOL ! ug/m3
          else

       AERO(AMP_AERO_MAP(nAMP)) =trm(i,j,l,n)/ AVOL       !  #/m3
          endif
       ENDDO

       if (L.eq.1) then     
!      Emis Mass [ug/m3/s] <-- trflux1[kg/m2/s]
#ifdef TRACERS_AMP_M4
      EMIS_MASS(2) =  MAX(trflux1(i,j,n_M_ACC_SU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(3) =  MAX(trflux1(i,j,n_M_BC1_BC)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(4) =  MAX(trflux1(i,j,n_M_OCC_OC)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(5) =  MAX(trflux1(i,j,n_M_DD1_DU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(6) =  MAX(trflux1(i,j,n_M_SSS_SS)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(10)=  MAX(trflux1(i,j,n_M_DD2_DU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
#else
      EMIS_MASS(1) =  MAX(trflux1(i,j,n_M_AKK_SU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(2) =  MAX(trflux1(i,j,n_M_ACC_SU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(3) =  MAX(trflux1(i,j,n_M_BC1_BC)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(4) =  MAX(trflux1(i,j,n_M_OCC_OC)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(5) =  MAX(trflux1(i,j,n_M_DD1_DU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(6) =  MAX(trflux1(i,j,n_M_SSA_SS)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(7) =  MAX(trflux1(i,j,n_M_SSC_SS)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
      EMIS_MASS(10)=  MAX(trflux1(i,j,n_M_DD2_DU)*1.d9 / AVOL,0.d0)
     &        *axyp(i,j)
#endif
        endif
!      Emis Mass [ug/m3/s] <-- trflux1[kg/s]
      EMIS_MASS(1) =  EMIS_MASS(1) + ((tr3Dsource(i,j,l,nVolcanic,n_M_AKK_SU)+
     *                                 tr3Dsource(i,j,l,nBiomass,n_M_AKK_SU)+
     *                                 tr3Dsource(i,j,l,nAircraft,n_M_AKK_SU))
     *                                 *1.d9 / AVOL)
      EMIS_MASS(2) =  EMIS_MASS(2) + ((tr3Dsource(i,j,l,nVolcanic,n_M_ACC_SU)+
     *                                 tr3Dsource(i,j,l,nBiomass,n_M_ACC_SU)+
     *                                 tr3Dsource(i,j,l,nAircraft,n_M_ACC_SU))
     *                                 *1.d9 / AVOL)
      EMIS_MASS(3) =  EMIS_MASS(3) + ((tr3Dsource(i,j,l,nBiomass,n_M_BC1_BC)+
     *                                 tr3Dsource(i,j,l,nAircraft,n_M_BC1_BC))
     *                               *1.d9 / AVOL)
      EMIS_MASS(4) =  EMIS_MASS(4) + ((tr3Dsource(i,j,l,nBiomass,n_M_OCC_OC)+
     *                                 tr3Dsource(i,j,l,nAircraft,n_M_OCC_OC))
     *                               *1.d9 / AVOL)

       CALL SPCMASSES(AERO,GAS,SPCMASS)

!=========
! WARNING: EMIS_MASS is only used to modify number, the mass is already modified in ATURB.
!=========
       CALL MATRIX(AERO,GAS,EMIS_MASS,TSTEP,TK,RH,PRES,AQSO4RATE,WUP,DT_AERO,PM) 
c       CALL SIZE_PDFS(AERO,PDF1,PDF2)
 
       DO n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
          if(AMP_NUMB_MAP(nAMP).eq. 0) then
      tr3Dsource(i,j,l,nChemistry,n) =((AERO(AMP_AERO_MAP(nAMP)) *AVOL *1.d-9)
     *        -trm(i,j,l,n)) /dtsrc 
          else
      tr3Dsource(i,j,l,nChemistry,n) =((AERO(AMP_AERO_MAP(nAMP)) *AVOL)
     *        -trm(i,j,l,n)) /dtsrc
          endif   
       ENDDO

      tr3Dsource(i,j,l,nOther,n_H2SO4) =((GAS(GAS_H2SO4)*AVOL *1.d-9)
     *        -trm(i,j,l,n_H2SO4)) /dtsrc 
      tr3Dsource(i,j,l,nChemistry,n_NH3)   =((GAS(GAS_NH3)*AVOL *1.d-9)
     *        -trm(i,j,l,n_NH3)) /dtsrc

#ifdef  TRACERS_SPECIAL_Shindell
      tr3Dsource(i,j,l,3,n_HNO3)  =((GAS(GAS_HNO3)*AVOL * 1.d-9)
     *        -trm(i,j,l,n_HNO3))/dtsrc
#endif
c       DT_AERO(:,:) = DT_AERO(:,:) * dtsrc !DT_AERO [# or ug/m3/s] , taijs [kg m2/kg(air)], byMA [kg/m2]

      if (l.eq.1) then
c - 2d acc output
c      PM1  [ug/m3] - [kg/kg(air)]
        taijs(i,j,ijts_AMPe(1))=taijs(i,j,ijts_AMPe(1)) + PM(1)*1.d-9*rgas*tk/pres 
c      PM2.5
        taijs(i,j,ijts_AMPe(2))=taijs(i,j,ijts_AMPe(2)) + PM(2)*1.d-9*rgas*tk/pres 
        ampPM2p5(i,j) = PM(2)*1.d-9*rgas*tk/pres 
c      PM10
        taijs(i,j,ijts_AMPe(3))=taijs(i,j,ijts_AMPe(3)) + PM(3)*1.d-9*rgas*tk/pres 
        ampPM10(i,j)  = PM(3)*1.d-9*rgas*tk/pres 
        endif

c Update physical properties per mode
       do n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
c Diagnostic of Processes - Sources and Sincs - timestep included
          if(AMP_NUMB_MAP(nAMP).eq. 0) then  !taijs [kg/s] -> in acc [kg/m2*s]
            do m=1,7
              taijs(i,j,ijts_AMPp(m,n)) =taijs(i,j,ijts_AMPp(m,n)) +(DT_AERO(m+8,AMP_AERO_MAP(nAMP))* AVOL * 1.d-9)
              if (itcon_amp(m,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(m+8,AMP_AERO_MAP(nAMP))* AVOL * 1.d-9),itcon_amp(m,n),n)
            end do
             
          else
            taijs(i,j,ijts_AMPp(1,n)) =taijs(i,j,ijts_AMPp(1,n))+(DT_AERO(2,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(1,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(2,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(1,n),n)
            taijs(i,j,ijts_AMPp(2,n)) =taijs(i,j,ijts_AMPp(2,n))+(DT_AERO(3,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(2,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(3,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(2,n),n)
            taijs(i,j,ijts_AMPp(3,n)) =taijs(i,j,ijts_AMPp(3,n))+(DT_AERO(1,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(3,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(1,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(3,n),n)
            do m=4,7
              taijs(i,j,ijts_AMPp(m,n)) =taijs(i,j,ijts_AMPp(m,n))+(DT_AERO(m,AMP_AERO_MAP(nAMP))* AVOL)
              if (itcon_amp(m,n).gt.0) call inc_diagtcb(i,j,(DT_AERO(m,AMP_AERO_MAP(nAMP))* AVOL),itcon_amp(m,n),n)
            end do

          endif
       select case (trname(n)) !taijs [kg * m2/kg air] -> in acc [kg/kg air]
      CASE('N_AKK_1 ','N_ACC_1 ','N_DD1_1 ','N_DS1_1 ','N_DD2_1 ','N_DS2_1 ','N_SSA_1','N_SSC_1','N_OCC_1 ','N_BC1_1 ',
     *     'N_BC2_1 ','N_BC3_1 ','N_DBC_1 ','N_BOC_1 ','N_BCS_1 ','N_MXX_1 ','N_OCS_1 ')
c - 3d acc output
        taijls(i,j,l,ijlt_AMPm(1,n))=taijls(i,j,l,ijlt_AMPm(1,n)) + DIAM(i,j,l,AMP_MODES_MAP(nAMP))
        taijls(i,j,l,ijlt_AMPm(3,n))=taijls(i,j,l,ijlt_AMPm(3,n)) + DIAM_dry(i,j,l,AMP_MODES_MAP(nAMP))
        taijls(i,j,l,ijlt_AMPm(2,n))=taijls(i,j,l,ijlt_AMPm(2,n)) + (NACTV(i,j,l,AMP_MODES_MAP(nAMP))*AVOL*byMA(l,i,j)/axyp(i,j))


c - 2d PRT Diagnostic
        if (itcon_AMPm(1,n) .gt.0) call inc_diagtcb(i,j,(DIAM(i,j,l,AMP_MODES_MAP(nAMP))*1d6),itcon_AMPm(1,n),n) 
        if (itcon_AMPm(3,n) .gt.0) call inc_diagtcb(i,j,(DIAM_dry(i,j,l,AMP_MODES_MAP(nAMP))*1d6),itcon_AMPm(3,n),n) 
        if (itcon_AMPm(2,n) .gt.0) call inc_diagtcb(i,j,NACTV(i,j,l,AMP_MODES_MAP(nAMP))*AVOL ,itcon_AMPm(2,n),n) 
       end select

      enddo !n
      ENDDO !i
      ENDDO !j
      ENDDO !l


#ifdef CACHED_SUBDD
      call find_groups('taijlh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          do n=ntmAMPi,ntmAMPe
            if (trim(subdd%name(k)) /= 'd'//trim(trname(n))) cycle
            nAMP=n-ntmAMPi+1
            sddarr3d(:,:,:)=diam(:,:,:,amp_modes_map(nAMP))
            call inc_subdd(subdd,k,sddarr3d)
          enddo ! n
        enddo ! k
      enddo ! igrp

! Tracer 2D I-J diags
      call find_groups('taijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case(trim(subdd%name(k)))
         case('ampDustload')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_DD1_DU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DS1_DU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DD2_DU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DS2_DU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DBC_DU),dim=3)
     *                  +sum(trm(:,:,:,n_M_MXX_DU),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 
         case('ampBCload')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_BC1_BC),dim=3)
     *                  +sum(trm(:,:,:,n_M_BC2_BC),dim=3)
     *                  +sum(trm(:,:,:,n_M_BC3_BC),dim=3)
     *                  +sum(trm(:,:,:,n_M_DBC_BC),dim=3)
     *                  +sum(trm(:,:,:,n_M_BOC_BC),dim=3)
     *                  +sum(trm(:,:,:,n_M_BCS_BC),dim=3)
     *                  +sum(trm(:,:,:,n_M_MXX_BC),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 
         case('ampNH4load')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_NH4),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 
         case('ampNO3load')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_NO3),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 
         case('ampOAload')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_OCC_OC),dim=3)
     *                  +sum(trm(:,:,:,n_M_BOC_OC),dim=3)
     *                  +sum(trm(:,:,:,n_M_MXX_OC),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 
         case('ampSO4load')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_AKK_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_ACC_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DD1_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DS1_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DD2_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DS2_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_SSA_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_OCC_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_BC1_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_BC2_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_BC3_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_BOC_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_BCS_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_DBC_SU),dim=3)
     *                  +sum(trm(:,:,:,n_M_MXX_SU),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 
         case('ampSSload')
         sddarr2d(:,:)= (sum(trm(:,:,:,n_M_SSA_SS),dim=3)
     *                  +sum(trm(:,:,:,n_M_SSC_SS),dim=3)
     *                  +sum(trm(:,:,:,n_M_MXX_SS),dim=3))
     *                  *byaxyp(:,:)
         call inc_subdd(subdd,k,sddarr2d) 

       end select
        
        enddo ! k
      enddo ! igrp

#endif  /* CACHED_SUBDD */

      RETURN
      END SUBROUTINE MATRIX_DRV
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE AMPtrdens(i,j,l,n)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the actual density per mode
!----------------------------------------------------------------------------------------------------------------------
      USE OldTracer_mod, only: trpdens
      USE AmpTracersMetadata_mod, only: AMP_MODES_MAP, AMP_trm_nm1,
     *  AMP_trm_nm2
      USE TRACER_COM, only: ntmAMPi, ntm, trm
      USE AMP_AEROSOL, only : AMP_dens, AMP_TR_MM
      USE AERO_CONFIG, ONLY: NMODES

      IMPLICIT NONE
      Integer :: i,j,l,n,x,nAMP
      real*8, dimension(:), allocatable :: trpdens_local

      allocate(trpdens_local(NTM))
      do x=1,NTM
        trpdens_local(x)=trpdens(x)
      enddo
 
      nAMP=n-ntmAMPi+1
      if(AMP_MODES_MAP(nAMP).gt.0)
     &  AMP_dens(i,j,l,AMP_MODES_MAP(nAMP)) = 
     &  sum(trpdens_local(AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP)) * 
     &  trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) 
     & / (sum(trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) + 1.0D-30)
      if (AMP_dens(i,j,l,AMP_MODES_MAP(nAMP)).le.0) 
     &  AMP_dens(i,j,l,AMP_MODES_MAP(nAMP)) = 
     &  trpdens_local(AMP_MODES_MAP(nAMP))

      deallocate(trpdens_local)

      RETURN
      END SUBROUTINE AMPtrdens
c -----------------------------------------------------------------
c -----------------------------------------------------------------
      SUBROUTINE AMPtrmass(i,j,l,n)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the actual molecular mass per mode
!----------------------------------------------------------------------------------------------------------------------
      USE OldTracer_mod, only: tr_mm
      USE AmpTracersMetadata_mod, only: AMP_MODES_MAP, AMP_trm_nm1,
     *  AMP_trm_nm2
      USE TRACER_COM, only: ntmAMPi, ntm, trm
      USE AERO_CONFIG, ONLY: NMODES
      USE AMP_AEROSOL, only : AMP_TR_MM

      IMPLICIT NONE
      Integer :: i,j,l,n,x,nAMP
      real*8, dimension(:), allocatable :: tr_mm_local

      allocate(tr_mm_local(NTM))
      do x=1,NTM 
        tr_mm_local(x)=tr_mm(x)
      enddo

      nAMP=n-ntmAMPi+1
      if(AMP_MODES_MAP(nAMP) > 0 .and.
     &  sum(trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) > 0. )
     &  AMP_TR_MM(i,j,l,AMP_MODES_MAP(nAMP)) = 
     &  sum(tr_mm_local(AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP)) * 
     &  trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP))) 
     & / sum(trm(i,j,l,AMP_trm_nm1(nAMP):AMP_trm_nm2(nAMP)))
      if (AMP_TR_MM(i,j,l,AMP_MODES_MAP(nAMP)).le.0) 
     &  AMP_TR_MM(i,j,l,AMP_MODES_MAP(nAMP)) = 
     &  tr_mm_local(AMP_MODES_MAP(nAMP))

      deallocate(tr_mm_local)

      RETURN
      END SUBROUTINE AMPtrmass
c -----------------------------------------------------------------
      SUBROUTINE SPCMASSES(AERO,GAS,SPCMASS)
!----------------------------------------------------------------------------------------------------------------------
!     Routine to calculate the total mass concentration of each model species:
!     SULF, BCAR, OCAR, DUST, SEAS, NO3, NH4. Aerosol water is not treated. 
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_SETUP, ONLY: SULF_MAP, BCAR_MAP, OCAR_MAP, DUST_MAP, SEAS_MAP
      USE AERO_PARAM
      USE AERO_CONFIG
      IMPLICIT NONE
      REAL(8) :: AERO(NAEROBOX)
      REAL(8) :: GAS(NGASES)    
      REAL(8) :: SPCMASS(NMASS_SPCS+2)
      SPCMASS(1) = SUM( AERO( SULF_MAP(:) ) ) + GAS( GAS_H2SO4 )
      SPCMASS(2) = SUM( AERO( BCAR_MAP(:) ) )
      SPCMASS(3) = SUM( AERO( OCAR_MAP(:) ) )
      SPCMASS(4) = SUM( AERO( DUST_MAP(:) ) )
      SPCMASS(5) = SUM( AERO( SEAS_MAP(:) ) )
      SPCMASS(6) = AERO( MASS_NO3 )           + GAS( GAS_HNO3 )
      SPCMASS(7) = AERO( MASS_NH4 )           + GAS( GAS_NH3  )

 
      RETURN
      END SUBROUTINE SPCMASSES
     
      SUBROUTINE SIZE_PDFS(AERO,PDF1,PDF2)
      USE AERO_PARAM, ONLY: PI6, DENSP, IXXX, IYYY, ILAY
      USE AERO_CONFIG, ONLY: NMODES, NAEROBOX,NBINS
      USE AERO_SETUP, ONLY: SIG0, CONV_DPAM_TO_DGN, NUMB_MAP, MODE_NAME
      USE AMP_AEROSOL, only: DIAM, DIAM_dry
      IMPLICIT NONE

      ! Arguments.
       REAL(8), INTENT(IN) :: AERO(NAEROBOX)! aerosol conc. [ug/m^3] or [#/m^3]

      ! Local variables. 

      INTEGER :: I, N
!      INTEGER, PARAMETER :: NBINS = 30! 200    ! number of bins [1]     defined in config 
      REAL(8) :: DGRID(NBINS)              ! fixed diameter        grid [um]
      REAL(8) :: MGRID(NBINS)              ! fixed mass/particle   grid [ug/particle]
      REAL(8) :: DLOWER(NBINS)             ! lower boundary fixed diameter grid [um]
      REAL(8) :: DUPPER(NBINS)             ! upper boundary fixed diameter grid [um]
      REAL(8) :: NTOT(NMODES)              ! number concentration for each mode [#/m^3]
      REAL(8) :: PDF(NBINS,2,NMODES)       ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: DNDLOGD(NMODES)           ! dN/dlog10(Dp) [ #/m^3]
      REAL(8) :: DMDLOGD(NMODES)           ! dM/dlog10(Dp) [ug/m^3]
      REAL(8) :: RDMIN                     ! reciprocal of DMIN to optimize coagulation [1/um]
      REAL(8) :: RDLOGDSC                  ! reciprocal of log10 of the grid spacing [1]
      REAL(8) :: SCALE, F, SUM1, SUM2      ! scratch variables 
      REAL(8) :: DMINL, DMAXL, DG          ! diameters [um]  
      REAL(8) :: FLN                       ! function for lognormal distribution [1]  
      REAL(8), PARAMETER :: DMIN =  0.001D+00   ! smallest particle diameter of the discrete grid [um]
      REAL(8), PARAMETER :: DMAX = 20.000D+00   ! largest  particle diameter of the discrete grid [um]


        DMAXL = DMAX
        DMINL = DMIN

      SCALE    = ( DMAXL / DMINL )**(1.0D+00/REAL(NBINS-1))
      RDLOGDSC = 1.0D+00 / LOG10( SCALE )
      RDMIN    = 1.0D+00 / DMINL
      DO I=1, NBINS
        DGRID(I)  = DMINL * SCALE**(I-1)                  ! [um]
        DLOWER(I) = DGRID(I) / SCALE**0.5D+00             ! [um]
        DUPPER(I) = DGRID(I) * SCALE**0.5D+00             ! [um]
        MGRID(I)  = 1.0D-06 * DENSP * PI6 * DGRID(I)**3   ! [ug/particle]
        DO N=1, NMODES
          DG = 1.0D+06 * DIAM(IXXX,IYYY,ILAY,N) * CONV_DPAM_TO_DGN(N)   ! convert [m] to [um] and Dbar to Dg
          NTOT(N) = AERO( NUMB_MAP(N) )
          F = NTOT(N) * FLN( DGRID(I), DG, SIG0(N) )
          PDF(I,1,N) = F * ( DUPPER(I) - DLOWER(I) )
          PDF(I,2,N) = PDF(I,1,N) * MGRID(I)
          DNDLOGD(N) = PDF(I,1,N) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
          DNDLOGD(N) = MAX( DNDLOGD(N), 1.0D-30 )
          DMDLOGD(N) = PDF(I,2,N) * RDLOGDSC                     ! [ug/m^3]
          DMDLOGD(N) = MAX( DMDLOGD(N), 1.0D-30 )
        ENDDO
c        WRITE(IUNIT,91) I, DGRID(I), DNDLOGD(:)
c        WRITE(JUNIT,91) I, DGRID(I), DMDLOGD(:)
      ENDDO

        PDF1(:) = 0.0D+00
        PDF2(:) = 0.0D+00
      DO N=1, NMODES
        DO I=1, NBINS
          PDF1(I) = PDF1(I)  + PDF(I,1,N)
          PDF2(I) = PDF2(I)  + PDF(I,2,N)
          SUM1 = SUM1 + PDF(I,1,N)
          SUM2 = SUM2 + PDF(I,2,N)
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE SIZE_PDFS


      REAL(8) FUNCTION FLN(X,XG,SIGMAG)
      REAL(8) :: X      ! particle radius or diameter [any units]
      REAL(8) :: XG     ! geometric mean radius or diameter [any units]
      REAL(8) :: SIGMAG ! geometric standard deviation [monodisperse = 1.0]
      REAL(8), PARAMETER :: SQRTTWOPI = 2.506628275D+00
      FLN = EXP(-0.5D+00*(LOG(X/XG)/LOG(SIGMAG))**2) / (X*LOG(SIGMAG)*SQRTTWOPI)
      RETURN
      END FUNCTION FLN
      subroutine alloc_tracer_amp_com(grid)
!@SUM  To alllocate arrays whose sizes now need to be determined
!@+    at run-time
!@auth Susanne Bauer
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use resolution, only     : im,lm
      use amp_aerosol
      use aero_config, only   : nmodes

      IMPLICIT NONE

      type (dist_grid), intent(in) :: grid
      integer :: ier, J_1H, J_0H, I_1H, I_0H
      logical :: init = .false.

      if(init)return
      init=.true.
    
      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO 

! I,J,L
      allocate(  AQsulfRATE(I_0H:I_1H,J_0H:J_1H,LM)   )
! other dimensions
      allocate(  DIAM(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  ampPM10(I_0H:I_1H,J_0H:J_1H)  )
      allocate(  ampPM2p5(I_0H:I_1H,J_0H:J_1H)  )
      allocate(  DIAM_dry(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  AMP_TR_MM(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  AMP_dens(I_0H:I_1H,J_0H:J_1H,LM,nmodes)  )
      allocate(  NACTV(I_0H:I_1H,J_0H:J_1H,LM,nmodes) )
      allocate(  VDDEP_AERO(I_0H:I_1H,J_0H:J_1H,nmodes,2))

      NACTV    = 1.0D-30
      DIAM     = 1.0D-30
      ampPM10  = 1.0D-30
      ampPM2p5 = 1.0D-30
      DIAM_dry = 1.0D-30

      return
      end subroutine alloc_tracer_amp_com
      
