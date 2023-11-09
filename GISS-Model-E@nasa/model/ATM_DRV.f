#include "rundeck_opts.h"

      subroutine atm_phase1
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      USE Dictionary_mod
      Use RESOLUTION, Only: IM,JM,LM
      Use CONSTANT,   Only: SHA,SHV
      Use GEOM,       Only: AXYP
      Use ATM_COM,    Only: MASUM,MA,MAOLD,PMID,PMIDOLD, PK,T,Q,QCL,QCI,
     &                      MUs,MVs, KEA
      Use DIAG_COM,   Only: AIJ=>AIJ_LOC,
     &                      IJ_dSE_DYN,IJ_dKE_DYN,IJ_dTE_DYN
      USE MODEL_COM
      Use DYNAMICS,   Only: nstep,nidyn,nfiltr,mfiltr,dt
      Use DOMAIN_DECOMP_ATM, Only: GRID, GLOBALSUM
      use domain_decomp_atm, only: writei8_parallel
      USE RANDOM
      USE GETTIME_MOD
      USE DIAG_COM, only : ia_src,ia_d5s,ia_d5d,ia_filt
     &     ,oa,koa
     &     ,MODD5S,NDAa, NDA5d,NDA5s
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use SystemTimers_mod
      use RAD_COM, only : nrad,modrd
      use seaice_com, only : si_atm,si_ocn,iceocn ! temporary until
      use lakes_com, only : icelak                ! melt_si calls
      use fluxes, only : atmocn,atmice            ! are moved

#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: Run,fvstate
#endif

#ifndef CUBED_SPHERE
      USE ATMDYN, only : DYNAM,SDRAG
     &     ,COMPUTE_DYNAM_AIJ_DIAGNOSTICS
#endif

#if defined(TRACERS_ON) || defined(TRACERS_OCEAN)
      USE TRACER_COM, only: mtrace
#endif

#if (defined(TRACERS_ON) || defined(TRACERS_OCEAN))  && defined(TRAC_ADV_CPU)
      USE TRACER_COM, only: mtradv
#endif

#ifdef TRACERS_TOMAS
      USE TRACER_COM, only : NBINS, n_ANUM,n_ASO4,n_AECIL, n_AECOB,
     &     n_AOCIL, n_AOCOB,n_ADUST,n_AH2O,n_ANACL
#endif

#ifdef SCM
      USE SCM_COM , only : nstepSCM
#endif

      implicit none

      INTEGER K,M,MSTART,MNOW,MODD5D,months,ioerr,Ldate,istart
      INTEGER :: MDUM = 0
      REAL*8 start,now, DTIME,TOTALT
      integer :: I,J,L,I_0,I_1,J_0,J_1
      Real*8  :: dSEpKE,MMGLOB
      Real*8,Dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     *                 GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *   SEINIT,KEINIT,SEFINAL,KEFINAL

#ifdef TRACERS_TOMAS
      integer :: n
#endif

      I_0 = GRID%I_STRT; I_1 = GRID%I_STOP
      J_0 = GRID%J_STRT; J_1 = GRID%J_STOP

C**** INITIALIZE TIME PARAMETERS
      NSTEP=(Itime-ItimeI)*NIdyn

C****
C**** INTEGRATE DYNAMIC TERMS (DIAGA AND DIAGB ARE CALLED FROM DYNAM)
C****
      CALL CHECKT ('DYNAM0')
         MODD5D=MOD(Itime-ItimeI,NDA5D)
         IF (MODD5D.EQ.0) IDACC(ia_d5d)=IDACC(ia_d5d)+1
#ifndef SCM
         IF (MODD5D.EQ.0) CALL DIAG5A (2,0)
#endif
         IF (MODD5D.EQ.0) CALL DIAGCA (1)

C**** Save MA and PMID before dynamics for Q advection and clouds
        MAOLD(:,:,:) =   MA(:,:,:)
      PMIDOLD(:,:,:) = PMID(:,:,:)

!**** Initialize total energy (J/m^2) before dynamics
      Call CONSERV_SE (SEINIT)
      Call CONSERV_KE (KEINIT)
      call startTimer('Atm. Dynamics')

#ifdef SCM
      nstepSCM = ITIME-ITIMEI
      write(0,*) 'nstepSCM ',nstepSCM
      if( nstepSCM > 0 ) call update_SCM_inputs
#endif

#if !defined(USE_FVCORE)
!**** Latitude-Longitude B-grid dynamics
      CALL DYNAM()
#endif

#if defined(USE_FVCORE)
!**** FV core dynamics
        IF (MOD(Itime-ItimeI,NDAA).eq.0) CALL DIAGA0
      call Run(fvstate)
#endif

#if defined(USE_FVCORE) && !defined(CUBED_SPHERE)
      CALL SDRAG (DTsrc)
#endif

#if defined(USE_FVCORE) || defined(SCM)
           If (MOD(Itime-ItimeI,NDAA) == 0)  Call DIAGA
#endif

#if defined(USE_FVCORE) && !defined(SCM)
           If (MOD(Itime-ItimeI,NDAA) == 0)  Call DIAGB
#endif

#if defined(USE_FVCORE) && defined(CUBE_SPHERE)
           If (MOD(Itime-ItimeI,NDAA) == 0)  Call EPFLUX
#endif

#ifndef CUBED_SPHERE
      call COMPUTE_DYNAM_AIJ_DIAGNOSTICS(MUs, MVs, DT)
#endif

      call COMPUTE_WSAVE
C**** Scale WM mixing ratios to conserve liquid water
      DO J=J_0,J_1
      DO I=I_0,I_1
!         WM(I,J,:) =  WM(I,J,;) * (MAOLD(:,I,J) / MA(:,I,J))
         QCL(I,J,:) = QCL(I,J,:) * (MAOLD(:,I,J) / MA(:,I,J))
         QCI(I,J,:) = QCI(I,J,:) * (MAOLD(:,I,J) / MA(:,I,J))
      END DO
      END DO
      CALL QDYNAM  ! Advection of Q by integrated fluxes
         CALL TIMER (NOW,MDYN)

#if defined(TRACERS_ON)
      CALL TrDYNAM   ! tracer dynamics
#endif

#if defined(TRACERS_ON) && defined(TRACERS_TOMAS)
!TOMAS- This next section of code ratios the higher order moments of
!       aerosol mass to those of aerosol number so the distributions of
!       aerosol mass and number within a grid cell are consistent
      do n=1,NBINS
         call momentfix(n_ANUM(1)-1+n, n_ASO4(1)-1+n)  !sulfate mass
         call momentfix(n_ANUM(1)-1+n, n_ANACL(1) -1+n)  !na+ mass
         call momentfix(n_ANUM(1)-1+n, n_AECOB(1)-1+n) !hydrophobic EC
         call momentfix(n_ANUM(1)-1+n, n_AECIL(1)-1+n)
         call momentfix(n_ANUM(1)-1+n, n_AOCOB(1)-1+n)
         call momentfix(n_ANUM(1)-1+n, n_AOCIL(1)-1+n)
         call momentfix(n_ANUM(1)-1+n, n_ADUST(1)-1+n)
         call momentfix(n_ANUM(1)-1+n, n_AH2O(1)-1+n)  !water mass
      enddo
#endif

#if defined(TRACERS_ON) && defined(TRAC_ADV_CPU)
         CALL TIMER (NOW,MTRADV)
#endif

#if defined(TRACERS_ON) && !defined(TRAC_ADV_CPU)
         CALL TIMER (NOW,MTRACE)
#endif

#ifndef SCM

! Why is this not in a subroutine?

!**** Adjust thermal energy to conserve total energy by dynamics
!**** Compute global energy generated by dynamics: dSEpKE (J)
!**** Compute global mass: MMGLOB (kg); KEFINAL stores MASUM*AXYP
!**** Subtract dSEpKE as heat proportional to mass
!**** Gaseous specific heat capacity: SHG = SHA + Q*(SHV-SHA)
      Call CONSERV_SE (SEFINAL)
      Call CONSERV_KE (KEFINAL)
      Do J=J_0,J_1  ;  Do I=I_0,I_1
              AIJ(I,J,IJ_dSE_DYN) = AIJ(I,J,IJ_dSE_DYN) + 
     +                  (SEFINAL(I,J) - SEINIT(I,J))
              AIJ(I,J,IJ_dKE_DYN) = AIJ(I,J,IJ_dKE_DYN) + 
     +                  (KEFINAL(I,J) - KEINIT(I,J))
         SEFINAL(I,J) = (SEFINAL(I,J) - SEINIT(I,J) +
     +                  (KEFINAL(I,J) - KEINIT(I,J))) * AXYP(I,J)
         KEFINAL(I,J) = MASUM(I,J) * AXYP(I,J)
         EndDo  ;  EndDo
      Call GLOBALSUM (GRID, SEFINAL, dSEpKE, ALL=.True.)
      Call GLOBALSUM (GRID, KEFINAL, MMGLOB, ALL=.True.)
      dSEpKE = dSEpKE / MMGLOB  !  (J/kg)
      Do J=J_0,J_1  ;  Do I=I_0,I_1
!        T(I,J,:) = T(I,J,:)
!    -            - dSEpKE / (PK(:,I,J)*(SHA + Q(I,J,:)*(SHV-SHA)))
!PRESENTLY VAPORMASS=0
         T(I,J,:) = T(I,J,:) - dSEpKE / (PK(:,I,J)*SHA)
           AIJ(I,J,IJ_dTE_DYN) = AIJ(I,J,IJ_dTE_DYN) + dSEpKE*MASUM(I,J) 
         EndDo  ;  EndDo
#endif
      call stopTimer('Atm. Dynamics')

#ifdef CACHED_SUBDD
      call get_subdd_vinterp_coeffs ! doing after dynamics
#endif

C**** Calculate tropopause level and pressure
      CALL CALC_TROP

#ifndef SCM
!**** Compute dynamic variables for the PBL
      CALL PGRAD_PBL
#endif

C**** calculate zenith angle for current time step
      CALL CALC_ZENITH_ANGLE
         CALL CHECKT ('DYNAM ')
         CALL TIMER (NOW,MSURF)
#ifndef SCM
         IF (MODD5D.EQ.0) CALL DIAG5A (7,NIdyn)
#endif
         IF (MODD5D.EQ.0) CALL DIAGCA (2)
#ifndef SCM
         IF (MOD(Itime,NDAY/2).eq.0) CALL DIAG7A
#endif
C****
C**** INTEGRATE SOURCE TERMS
C****

c calculate KE before atmospheric column physics
         call calc_kea_3d(kea)

#ifdef CUBED_SPHERE
c GWDRAG, SDRAG considered as column physics so that their KE
c dissipation gets included in the KE->PE adjustment
      CALL GWDRAG
      CALL SDRAG (DTsrc)
#endif

         IDACC(ia_src)=IDACC(ia_src)+1
         MODD5S=MOD(Itime-ItimeI,NDA5S)
         atmocn%MODD5S = MODD5S
         IF (MODD5S.EQ.0) IDACC(ia_d5s)=IDACC(ia_d5s)+1
#ifndef SCM
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAG5A (1,0)
#endif
         IF (MODD5S.EQ.0.AND.MODD5D.NE.0) CALL DIAGCA (1)

C**** FIRST CALL MELT_SI SO THAT TOO SMALL ICE FRACTIONS ARE REMOVED
C**** AND ICE FRACTION CAN THEN STAY CONSTANT UNTIL END OF TIMESTEP
! todo: move melt_si(ocean) to the end of the ocean driver, and
! possibly unite melt_si(lakes) with the rest of the lakes calls
      CALL MELT_SI(si_ocn,iceocn,atmocn,atmice)
      CALL MELT_SI(si_atm,icelak,atmocn,atmice)
      call seaice_to_atmgrid(atmice)
         CALL UPDTYPE
         CALL TIMER (NOW,MSURF)

C**** CONDENSATION, SUPER SATURATION AND MOIST CONVECTION
      CALL CONDSE
         CALL CHECKT ('CONDSE')
         CALL TIMER (NOW,MCNDS)
#ifndef SCM
         IF (MODD5S.EQ.0) CALL DIAG5A (9,NIdyn)
#endif
         IF (MODD5S.EQ.0) CALL DIAGCA (3)

C**** RADIATION, SOLAR AND THERMAL
      MODRD=MOD(Itime-ItimeI,NRAD)
      CALL RADIA
         CALL CHECKT ('RADIA ')
         CALL TIMER (NOW,MRAD)
#ifndef SCM
         IF (MODD5S.EQ.0) CALL DIAG5A (11,NIdyn)
#endif
         IF (MODD5S.EQ.0) CALL DIAGCA (4)

#ifdef TRACERS_ON
C**** Calculate non-interactive tracer surface sources and sinks
         call set_tracer_2Dsource
         CALL TIMER (NOW,MTRACE)
C**** Add up the non-interactive tracer surface sources.
      call sum_prescribed_tracer_2Dsources(dtsrc)
      call set_strattroptracer_diag(dtsrc)
#endif

      call atm_phase1_exports

      return
      end subroutine atm_phase1


      subroutine atm_phase1_exports
! Copies fields calculated by the atmosphere into the data structures
! seen by physics of the surface components (ocean, ice, land).
! Some fields are already type-classified (e.g. radiative fluxes).
! A per-type breakdown will soon be applied to other fields as
! appropriate (e.g. fields which vary with height will have
! different values over the ocean and the other portions of a gridbox).
! Some fields have already been stored in atmsrf%xxx and are not
! referenced here in Step 1.  For temporary convenience,
! fields depending on surface pressure are being copied into the
! per-surface-type structures by MAtoPMB because subroutine FILTER
! is currently being called after the main surface physics, but
! before special "daily" surface coding which sometimes requires
! surface pressure - once FILTER is absorbed into DYNAM this hack
! can be eliminated.
      use fluxes, only : atmsrf,asflx4,prec,eprec
#ifdef TRACERS_ON
#ifdef TRACERS_WATER
      use fluxes, only : trprec
#endif
#endif
      use rad_com, only : trhr,fsf,trsurf,cosz1
      use fluxes, only : atmocn
      use rad_com, only : dirvis,fsrdif,dirnir,difnir
      implicit none
      integer :: it

#ifdef SCM
      ! update SST every time step for SCM cases
      call daily_ocean(.false.,atmocn)
#endif

      ! Step 1: copy fields not already stored in atmsrf%xxx
      atmsrf%prec = prec
      atmsrf%eprec = eprec
#ifdef TRACERS_WATER
      atmsrf%trprec = trprec
#endif
      atmsrf%cosz1 = cosz1
      atmsrf%flong = trhr(0,:,:)


      ! Copy from gridbox-mean data structure to per-type structures
      do it=1,4
        asflx4(it)%atm_exports_phase1 = atmsrf%atm_exports_phase1
#ifdef TRACERS_WATER
        asflx4(it)%tratm_exports_phase1 = atmsrf%tratm_exports_phase1
#endif
      enddo

      ! fields already available per type
      do it=1,4
        asflx4(it)%fshort = fsf(it,:,:)
        asflx4(it)%trup_in_rad = trsurf(it,:,:)
      enddo


      ! miscellaneous fields only defined for some types
      if (allocated(atmocn%dirvis)) then
        atmocn % DIRVIS = DIRVIS
        atmocn % DIFVIS = FSRDIF
        atmocn % DIRNIR = DIRNIR
        atmocn % DIFNIR = DIFNIR
      endif

      return
      end subroutine atm_phase1_exports


      subroutine atm_exports_phasesrf
! Copies fields calculated by the atmosphere into the data structures
! seen by the physics of the surface components (ocean, ice, land).
! Per-component exports will be introduced as appropriate (e.g. fields
! which vary with height will have different values over the ocean and
! the other portions of a gridbox).
      use fluxes, only : atmsrf,asflx4
      implicit none

      integer :: it

      call get_atm_layer1

      do it=1,4
        asflx4(it)%atm_exports_phasesrf = atmsrf%atm_exports_phasesrf
#ifdef TRACERS_ON
        asflx4(it)%tratm_exports_phasesrf =
     &       atmsrf%tratm_exports_phasesrf
#endif
      enddo

      return
      end subroutine atm_exports_phasesrf

      subroutine get_atm_layer1
C**** Copies first-layer atm. conditions into the 2D arrays
C**** in the atm-surf. coupling data structure.
      use fluxes, only : atmsrf
      use domain_decomp_atm, only : grid, getDomainBounds
      use atm_com, only : t,q,ualij,valij
      use geom, only : imaxj,byaxyp
#ifdef TRACERS_ON
      use tracer_com, only : ntm,trm
#endif
      implicit none
      integer :: n,i,j,i_0,i_1,j_0,j_1
c
      call getDomainBounds(grid, i_strt=i_0,i_stop=i_1,
     &  j_strt=j_0,j_stop=j_1)
c
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        atmsrf%temp1(i,j) = t(i,j,1)
        atmsrf%q1(i,j) = q(i,j,1)
        atmsrf%u1(i,j) = ualij(1,i,j)
        atmsrf%v1(i,j) = valij(1,i,j)
      enddo
      enddo
#ifdef TRACERS_ON
      do n=1,ntm
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        atmsrf%trm1(n,i,j) = trm(i,j,1,n)*byaxyp(i,j)
      enddo
      enddo
      enddo
#endif
      end subroutine get_atm_layer1


      subroutine atm_phase2
      USE TIMINGS, only : ntimemax,ntimeacc,timing,timestr
      use resolution, only : lm
      USE MODEL_COM
      USE DYNAMICS, only : nidyn,nfiltr,mfiltr
      USE GETTIME_MOD
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE TRACER_COM, only: mtrace
#endif
      USE DIAG_COM, only : kvflxo,oa,koa,ia_filt
     &     ,MODD5S,NDAa, NDA5d,NDA5s,NDA4
#ifndef CACHED_SUBDD
      USE SUBDAILY, only : nsubdd,get_subdd,accSubdd
#endif
#ifndef CUBED_SPHERE
#ifndef SCM
      USE ATMDYN, only : FILTER
#endif
#endif
      USE FLUXES, only : atmocn,atmice
      use TimerPackage_mod, only: startTimer => start
      use TimerPackage_mod, only: stopTimer => stop
      use SystemTimers_mod
      implicit none

      REAL*8 start,now

      call seaice_to_atmgrid(atmice)
      CALL ADVSI_DIAG(atmocn,atmice) ! needed to update qflux model, dummy otherwise
C**** SAVE some noon GMT ice quantities
      IF (MOD(Itime+1,NDAY).ne.0 .and. MOD(Itime+1,NDAY/2).eq.0)
     &        call vflx_OCEAN

C**** IF ATURB is used in rundeck then this is a dummy call
C**** CALCULATE DRY CONVECTION ABOVE PBL
      CALL ATM_DIFFUS (2,LM-1,dtsrc)
         CALL CHECKT ('DRYCNV')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (9)

C**** UPDATE DIAGNOSTIC TYPES
         CALL UPDTYPE
C**** ADD DISSIPATED KE FROM COLUMN PHYSICS CALCULATION BACK AS LOCAL HEAT
      CALL DISSIP ! uses kea calculated before column physics
         CALL CHECKT ('DISSIP')
         CALL TIMER (NOW,MSURF)
         IF (MODD5S.EQ.0) CALL DIAGCA (7)
#ifndef SCM
         IF (MODD5S.EQ.0) CALL DIAG5A (12,NIdyn)
#endif

#if defined(CUBED_SPHERE) || defined(SCM) 
      IDACC(ia_filt)=IDACC(ia_filt)+1 ! prevent /0
#else
C**** SEA LEVEL PRESSURE FILTER
      IF (MFILTR.GT.0.AND.MOD(Itime-ItimeI,NFILTR).EQ.0) THEN
           IDACC(ia_filt)=IDACC(ia_filt)+1
           IF (MODD5S.NE.0) CALL DIAG5A (1,0)
           CALL DIAGCA (1)
           CALL FILTER
           CALL CHECKT ('FILTER')
           CALL TIMER (NOW,MDYN)
           CALL DIAG5A (14,NFILTR*NIdyn)
           CALL DIAGCA (8)
      END IF
#endif
#ifdef TRACERS_ON
#ifdef CUBED_SPHERE
! Reinitialize instantaneous consrv qtys (every timestep since
! DIAGTCA is called every timestep for 3D sources)
      CALL DIAGCA (1) ! was not called w/ SLP filter
#endif
C**** 3D Tracer sources and sinks
C**** Tracer gravitational settling for aerosols
      CALL TRGRAV
C**** Tracer radioactive decay (and possible source)
      CALL TDECAY
C**** Calculate 3D tracers sources and sinks

      call tracer_3Dsource
C**** Accumulate tracer distribution diagnostics
      CALL TRACEA
         CALL TIMER (NOW,MTRACE)
         CALL CHECKT ('T3DSRC')
#endif
      call accum_ma_ia_src
C****
C**** WRITE SUB-DAILY DIAGNOSTICS EVERY NSUBDD hours
C****
#ifdef CACHED_SUBDD
      call accum_subdd_atm
#else
      if (Nsubdd.ne.0) then
        call accSubdd
        if (mod(Itime+1,Nsubdd).eq.0) call get_subdd
      end if
#endif
#ifdef TRACERS_DUST
      call ahourly
#endif

#ifndef SCM
      IF (MOD(Itime+1-ItimeI,NDA4).EQ.0) CALL DIAG4A ! at hr 23 E-history
#endif

      IF (Kvflxo.EQ.0.) OA(:,:,4:KOA)=0. ! to prepare for future saves

      return
      end subroutine atm_phase2


      SUBROUTINE INPUT_atm (istart,istart_fixup,do_IC_fixups,
     &     is_coldstart,KDISK_restart,IRANDI)

C****
C**** THIS SUBROUTINE SETS THE PARAMETERS IN THE C ARRAY, READS IN THE
C**** INITIAL CONDITIONS, AND CALCULATES THE DISTANCE PROJECTION ARRAYS
C****
      USE Dictionary_mod
      USE CONSTANT, only : grav
      USE FLUXES, only : nisurf,atmocn,atmice
      USE RESOLUTION, only : ls1=>ls1_nominal,plbot, MFIX,MFRAC
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only :
     *      irand,idacc ,nday,dtsrc ,iyear1,itime,itimei,itimee
     *     ,mdyn,mcnds,mrad,msurf,mdiag, calendar
#ifndef SCM
      USE DIAG_ZONAL, only : imlon
#endif
      USE RANDOM
      USE DYNAMICS, only : USE_UNR_DRAG
      USE ATM_COM, only : ij_debug
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      USE SOMTQ_COM, only : mz,tmom
      USE ATM_COM, only : ma,t
      USE TRACER_COM,only: MTRACE,daily_z
#ifdef TRACERS_SPECIAL_Shindell
     *     ,mchem
#endif
#ifdef TRAC_ADV_CPU
      USE TRACER_COM,only: MTRADV
#endif
#endif
      USE SOIL_DRV, only: init_LSM,daily_earth
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, AM_I_ROOT
#ifndef CUBED_SPHERE
      USE ATMDYN, only : init_ATMDYN
#endif
#ifdef USE_ESMF
      use ATM_COM, only : atmclock
#endif
#ifdef USE_FVCORE
      USE FV_INTERFACE_MOD, only: fvstate,initialize
#endif
#ifndef SCM
#ifndef CUBED_SPHERE
      use UNRDRAG_COM, only: init_UNRDRAG
#endif
#endif
      IMPLICIT NONE
!@var istart start(1-8)/restart(>8)  option
      integer :: istart,istart_fixup,do_IC_fixups
      LOGICAL :: is_coldstart
      INTEGER :: KDISK_restart
      INTEGER :: IRANDI

      INTEGER I,J,L,K

!@nlparam HOURI,DATEI,MONTHI,YEARI        start of model run
!@nlparam TIMEE,HOURE,DATEE,MONTHE,YEARE,IHOURE   end of model run
!@var  IHRI,IHOURE start and end of run in hours (from 1/1/IYEAR1 hr 0)
c      INTEGER ::   HOURI=0 , DATEI=1, MONTHI=1, YEARI=-1, IHRI=-1,
c     *    TIMEE=-1,HOURE=0 , DATEE=1, MONTHE=1, YEARE=-1, IHOURE=-1

      LOGICAL :: redoGH, iniPBL, inilake, iniSNOW
      INTEGER :: J_0, J_1

#ifdef USE_ESMF
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      write(6,*) 'mpi-zone',J_0,' - ',J_1
#endif

      CALL SET_TIMER("ATMOS. DYNAM",MDYN)
      CALL SET_TIMER("CONDENSATION",MCNDS)
      CALL SET_TIMER("   RADIATION",MRAD)
      CALL SET_TIMER("     SURFACE",MSURF)
      CALL SET_TIMER(" DIAGNOSTICS",MDIAG)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      CALL SET_TIMER("     TRACERS",MTRACE)
#ifdef TRAC_ADV_CPU
      CALL SET_TIMER(" TRACER ADV.",MTRADV)
#endif
#endif
#ifdef TRACERS_SPECIAL_Shindell
      CALL SET_TIMER("   CHEMISTRY",MCHEM)
#endif

      call sync_param( "ij_debug",ij_debug , 2)

C****
C**** Set some documentary parameters in the database
C****
      call set_param("IM",IM,'o')
      call set_param("JM",JM,'o')
      call set_param("LM",LM,'o')
      call set_param("LS1",LS1,'o')
      call set_param("PLBOT",Plbot,LM+1,'o')
      call set_param("MFIX" ,MFIX ,LM  ,'o')
      call set_param("MFRAC",MFRAC,LM  ,'o')

      if(istart.eq.2) then
        call read_aic()
        call aic_part2()
      endif

#ifdef USE_ESMF
      call init_esmf_clock_for_modelE( int(dtsrc), atmclock )
#endif

C****
C**** IRANDI seed for random perturbation of current state (if/=0)
C****        tropospheric temperatures are changed by at most 1 degree C
      IF (ISTART.LT.10 .AND. IRANDI.NE.0) THEN
        CALL RINIT (IRANDI)
        CALL PERTURB_TEMPS
        IF (AM_I_ROOT())
     *       WRITE(6,*) 'Initial conditions were perturbed !!',IRANDI
      END IF

#ifdef TRACERS_ON
      if(istart.le.2) then
         Call COMPUTE_GZ (MA,T,TMOM(MZ,:,:,:), DAILY_Z)
        daily_z = daily_z/grav
      endif
      call initTracerGriddedData(istart<=2)
#endif

C****
      CALL RINIT (IRAND)
c Note on FFT initialization: IMLON is defined by the diag_zonal module,
c not by the resolution module.  IMLON==IM for a latlon grid.
#ifndef SCM
      CALL FFT0 (IMLON)  ! CALL FFT0(IM)
      CALL init_QUS(grid,im,jm,lm)
#endif
#ifndef CUBED_SPHERE
      CALL init_ATMDYN
#endif
      call init_sdrag
C**** Initialize nudging
#ifdef NUDGE_ON
      CALL NUDGE_INIT
#endif
C****
C**** Initialize the gravity wave drag scheme
C****
      CALL init_GWDRAG
      call sync_param( "USE_UNR_DRAG", USE_UNR_DRAG )
#ifndef SCM
#ifndef CUBED_SPHERE
      if (USE_UNR_DRAG==1) CALL init_UNRDRAG(calendar)
#endif
#endif

#ifdef SCM
!      read SCM data and initialize model
       call init_SCM
#endif

      CALL init_CLD(istart)
      CALL init_RAD(istart)
      CALL daily_orbit(.false.)             ! not end_of_day
      CALL daily_ch4ox(.false.)             ! not end_of_day
      CALL daily_RAD(.false.)
      if(istart.eq.2) call read_rad_ic

      call atm_phase1_exports

C**** Initialize lake variables (including river directions)
      if(istart.eq.2) call read_agrice_ic
      iniLAKE = is_coldstart
      CALL init_lakeice(inilake,do_IC_fixups)
c      call stop_model('set_noice_defaults prob that fwater==flake'//
c     &     ' not initialized yet?',255)
      call seaice_to_atmgrid(atmice) ! set gtemp etc.
      CALL init_LAKES(inilake,istart_fixup)
C****
C**** INITIALIZE GROUND HYDROLOGY ARRAYS (INCL. VEGETATION)
C**** Recompute Ground hydrology data if redoGH (new soils data)
C****
      if(istart.eq.2) call read_landsurf_ic
      iniSNOW = is_coldstart       ! extract snow data from first soil layer
      redoGH = .false.
      CALL init_LSM(DTsrc/NIsurf,redoGH,iniSNOW,inilake,ISTART)

      CALL daily_EARTH(.false.)            ! not end_of_day

#ifdef CALCULATE_FLAMMABILITY
      CALL init_flammability
#endif

C**** Initialize land ice (must come after oceans)
      if(istart.eq.2) call read_landice_ic
      CALL init_LI(istart_fixup)

C**** Initialize pbl (and read in file containing roughness length data)
      iniPBL = is_coldstart
#ifndef CUBED_SPHERE /* until a better solution is found */
      if (iniPBL) call recalc_agrid_uv   ! PBL needs A-grid winds
#endif
      CALL init_pbl(iniPBL,istart)

! note there is an all-component call to reset_diag in init_diag
      CALL init_DIAG
      CALL UPDTYPE   ! for atm-grid diags

! init tasks transplanted from main program.  cannot disperse to
! corresponding components w/o changing results.
! DAILY_atmdyn adjusts the global mean mass when ISTART=2, so
! moving it before init of other atm components will change results
      CALL DAILY_atmdyn(.false.)           ! not end_of_day
#ifdef USE_FVCORE
C****
C**** Initialize FV dynamical core (ESMF component) if requested
C**** For restarts/continuations, FV state and import files are
C**** copied to the appropriate names by this procedure, and for
C**** cold starts the required IC files are generated.
C****
      ! See FV_INTERFACE.F90
      Call Initialize(fvstate, istart, kdisk_restart)
#endif

! this daily_OCEAN call belongs in ocean init, but
! ISTART=2 result differs if daily_OCEAN comes before init_pbl,
! since daily_OCEAN replaces the GIC values of gtemp
! with the values from the prescribed SST file
      CALL daily_OCEAN(.false.,atmocn)            ! not end_of_day
! need to do this after prescribed-ice daily_OCEAN changes si frac.
! not needed once daily_OCEAN is moved before atm init.
      call seaice_to_atmgrid(atmice) ! debug
      CALL UPDTYPE

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      CALL daily_tracer(.false.)
#endif
      CALL CHECKT ('INPUT ')

      if (AM_I_ROOT()) then
         WRITE (6,'(A14,4I4)') "IM,JM,LM,LS1=",IM,JM,LM,LS1
         WRITE (6,*) "PLbot=",PLbot
      end if

      return

      end subroutine INPUT_atm

      subroutine alloc_drv_atm()
#ifdef SCM
      use Dictionary_mod, only : sync_param
#endif
c Driver to allocate arrays that become dynamic as a result of
c set-up for MPI implementation
      USE DOMAIN_DECOMP_ATM, ONLY : grid,init_grid
      use MODEL_COM, only: calendar
#ifdef GLINT2
      USE DOMAIN_DECOMP_ATM, ONLY : glint2
      use MpiSupport_mod, only: ROOT_PROCESS
      Use glint2_modele
#endif
      USE RESOLUTION, only : im,jm,lm
#ifndef CUBED_SPHERE
#ifndef SCM
      USE MOMENTS, only : initMoments
#endif
#endif
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      use TRACER_COM, only: initTracerCom, alloc_tracer_com
#ifndef TRACERS_ATM_ONLY
      use ghy_tracers, only: initGhyTracers
#endif
#endif
#ifdef TRACERS_AEROSOLS_SEASALT
      use tracers_seasalt, only: alloc_seasalt_sources
#endif  /* TRACERS_AEROSOLS_SEASALT */
      use geom, only: geom_atm
      IMPLICIT NONE
#ifdef GLINT2
      include 'mpif.h'      ! Needed for GLINT2
#endif

c initialize the atmospheric domain decomposition
c for now, CREATE_CAP is only relevant to the cubed sphere grid
      call init_grid(grid, im, jm, lm, CREATE_CAP=.true.)

      call geom_atm

#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call initTracerCom
#ifndef TRACERS_ATM_ONLY
      call initGhyTracers
#endif
#endif

#ifdef GLINT2
      glint2 = glint2_modele_new('GLINT2', 6, 'm', 1,
     &    im, jm,
     &    grid%i_strt_halo, grid%i_stop_halo,
     &    grid%j_strt_halo, grid%j_stop_halo,
     &    grid%i_strt, grid%i_stop, grid%j_strt, grid%j_stop,
     &    grid%j_strt_skp, grid%j_stop_skp,
     &    MPI_COMM_WORLD, ROOT_PROCESS)
#endif /* GLINT2 */

      call alloc_dynamics(grid)
      call alloc_atm_com(grid)
      call alloc_smomtq(grid)

      call alloc_fluxes !(grid)
      call alloc_clouds_com(grid)
      call alloc_ghy_com(grid)
      call alloc_pbl_com(grid)
      call alloc_diag_com(grid)
      call alloc_diag_loc(grid)
      call alloc_strat_com(grid)
      call alloc_rad_com(grid)
      call alloc_lakes(grid)
      call alloc_lakes_com(grid)
      call alloc_landice_com(grid)

#ifdef CALCULATE_FLAMMABILITY
      call alloc_flammability(grid)
#endif
#ifdef BIOGENIC_EMISSIONS
      call alloc_biogenic_emis(grid)
#endif
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN) || (defined TRACERS_WATER)
      call alloc_tracer_com(grid)
#ifdef TRACERS_DRYDEP
      call alloc_trdrydep(grid)
#endif
#ifdef TRACERS_SPECIAL_Lerner
      call alloc_tracer_special_lerner_com(grid)
      call alloc_linoz_chem_com(grid)
#endif
#if (defined CALCULATE_LIGHTNING) || (defined TRACERS_SPECIAL_Shindell)
      call alloc_lightning(grid)
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call alloc_trchem_shindell_com(grid)
      call alloc_tracer_sources(grid)
#endif
#ifdef TRACERS_AEROSOLS_SEASALT
      call alloc_seasalt_sources(grid)
#endif  /* TRACERS_AEROSOLS_SEASALT */
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      call alloc_aerosol_sources(grid)
#endif
#ifdef TRACERS_AMP
      call alloc_tracer_amp_com(grid)
#endif
#ifdef TRACERS_TOMAS
      call alloc_tracer_tomas_com(grid)
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
      CALL alloc_dust(grid)
#endif
#endif
#ifndef SCM
      call alloc_tracer_adv(grid)
#endif
!!! should be done in init_module_ent
      call alloc_ent_com(grid)
#ifdef TRACERS_ON
      call alloc_trdiag_com
#ifdef TRACERS_SPECIAL_Shindell
      call interpolateAltitude()
#endif
#endif
#ifdef NUDGE_ON
      call alloc_nudge(grid)
#endif

#ifdef SCM
      call alloc_SCM_COM()
#endif
#ifndef CUBED_SPHERE
#ifndef SCM
      call initMoments
#endif
#endif
      end subroutine alloc_drv_atm

      subroutine def_rsf_atmvars(fid)
!@sum  def_rsf_atmvars defines atm prognostic array structure in rsf
!@auth M. Kelley
!@ver  beta
      implicit none
      integer :: fid
      call def_rsf_atm    (fid)
      call def_rsf_lakes  (fid)
      call def_rsf_agrice (fid)
      call def_rsf_icedyn (fid)
      call def_rsf_earth  (fid)
      call def_rsf_soils  (fid)
      call def_rsf_vegetation(fid)
      call def_rsf_veg_related(fid)
      call def_rsf_snow   (fid)
      call def_rsf_landice(fid)
      call def_rsf_bldat  (fid)
      call def_rsf_pbl    (fid)
      call def_rsf_clouds (fid)
#ifdef AUTOTUNE_LIGHTNING
      call def_rsf_lightning(fid)
#endif
      call def_rsf_somtq  (fid)
      call def_rsf_rad    (fid)
#ifdef CALCULATE_FLAMMABILITY
      call def_rsf_flammability(fid)
#endif
#ifdef TRACERS_ON
      call tracerIO(fid, 'define')
#endif
      call def_rsf_subdd  (fid)
      call def_rsf_fluxes (fid)
      return
      end subroutine def_rsf_atmvars

      subroutine new_io_atmvars(fid,iorw)
      use model_com, only: ioread, iowrite
      implicit none
      integer, intent(in) :: fid,iorw
      call new_io_atm    (fid,iorw)
      call new_io_lakes  (fid,iorw)
      call new_io_agrice (fid,iorw)
      call new_io_earth  (fid,iorw)
      call new_io_soils  (fid,iorw)
      call new_io_vegetation  (fid,iorw)
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
      call new_io_veg_related(fid,iorw)
        !call io_ent    (kunit,iaction,ioerr) ! io_vegetation handles ent
      call new_io_snow   (fid,iorw)
      call new_io_landice(fid,iorw)
      call new_io_bldat  (fid,iorw)
      call new_io_pbl    (fid,iorw)
      call new_io_clouds (fid,iorw)
#ifdef AUTOTUNE_LIGHTNING
      call new_io_lightning(fid,iorw)
#endif
      call new_io_somtq  (fid,iorw)
      call new_io_rad    (fid,iorw)
      call new_io_icedyn (fid,iorw)
#ifdef CALCULATE_FLAMMABILITY
      call new_io_flammability(fid,iorw)
#endif
#ifdef TRACERS_ON
      select case (iorw)
      case (ioread)
         call tracerIO(fid, 'read_dist')
      case (iowrite)
         call tracerIO(fid, 'write_dist')
      end select

#endif
      call new_io_subdd  (fid,iorw)
      call new_io_fluxes (fid,iorw)
      return
      end subroutine new_io_atmvars

      subroutine daily_atm(end_of_day)
      use filemanager
      use MODEL_COM, only: nday,itime
      use DYNAMICS, only : nidyn
      USE SOIL_DRV, only: daily_earth
      use diag_com, only : kvflxo,iu_vflxo,oa,koa
      use domain_decomp_atm, only: grid,writei8_parallel
      implicit none
      logical, intent(in) :: end_of_day ! not used yet

#ifndef SCM
      call DIAG5A (1,0)
#endif
      call DIAGCA (1)
      CALL daily_atmdyn(.true.)  ! end_of_day
      CALL daily_orbit(.true.)   ! end_of_day
      CALL daily_ch4ox(.true.)   ! end_of_day
      call daily_RAD(.true.)

      call daily_LAKE
      call daily_EARTH(.true.)  ! end_of_day

      call daily_LI
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call daily_tracer(.true.)
#endif
      call CHECKT ('DAILY ')
#ifndef SCM
      call DIAG5A (16,NDAY*NIdyn)
#endif
      call DIAGCA (10)
      call sys_flush(6)
      call UPDTYPE

C****
C**** WRITE INFORMATION FOR OHT CALCULATION EVERY 24 HOURS
C****
      IF (Kvflxo.NE.0.) THEN
        call writei8_parallel(grid,iu_vflxo,nameunit(iu_vflxo),oa,Itime)
C**** ZERO OUT INTEGRATED QUANTITIES
        OA(:,:,4:KOA)=0.
      END IF

      return
      end subroutine daily_atm

      subroutine finalize_atm
#ifndef CACHED_SUBDD
      USE SUBDAILY, only : close_subdd
#endif
#ifdef USE_FVCORE
      USE MODEL_COM, only : kdisk
      USE FV_INTERFACE_MOD, only: fvstate
      USE FV_INTERFACE_MOD, only: Finalize
#endif
      implicit none
#ifdef USE_FVCORE
         call Finalize(fvstate, kdisk)
#endif

#ifndef CACHED_SUBDD
C**** CLOSE SUBDAILY OUTPUT FILES
      CALL CLOSE_SUBDD
#endif

      return
      end subroutine finalize_atm


      SUBROUTINE CHECKT (SUBR)
!@sum  CHECKT Checks arrays for NaN/INF and reasonablness
!@auth Original Development Team

C**** CHECKT IS TURNED ON BY SETTING QCHECK=.TRUE. IN NAMELIST
C**** REMEMBER TO SET QCHECK BACK TO .FALSE. AFTER THE ERRORS ARE
C**** CORRECTED.
      USE CONSTANT, only : tf
      USE RESOLUTION, only : ls1=>ls1_nominal
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v,t,p,q,qcl,qci,pk
#ifdef BLK_2MOM
#endif
      USE MODEL_COM
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, AM_I_ROOT
      USE soil_drv, only : checke
      IMPLICIT NONE
      INTEGER I,J,L
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H, J_1H, I_0,I_1, I_0H,I_1H, njpol
      INTEGER :: I_0STG,I_1STG,J_0STG,J_1STG
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     *     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
c      I_0STG = grid%I_STRT_STGR
c      I_1STG = grid%I_STOP_STGR
      I_0STG = I_0
      I_1STG = I_1
      J_0STG = grid%J_STRT_STGR
      J_1STG = grid%J_STOP_STGR
      njpol = grid%J_STRT_SKP-grid%J_STRT

      IF (QCHECK) THEN
C**** Check all prog. arrays for Non-numbers
        CALL CHECK3B(U(I_0STG:I_1STG,J_0STG:J_1STG,:),
     &       I_0STG,I_1STG,J_0STG,J_1STG,0,LM,SUBR,'u     ')
        CALL CHECK3B(V(I_0STG:I_1STG,J_0STG:J_1STG,:),
     &       I_0STG,I_1STG,J_0STG,J_1STG,0,LM,SUBR,'v     ')
        CALL CHECK3B(T(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'t     ')
        CALL CHECK3B(Q(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'q     ')
        CALL CHECK3B(P(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &       SUBR,'p     ')
!       CALL CHECK3B(WM(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
!    &       SUBR,'wm    ')
        CALL CHECK3B(QCL(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'qcl   ')
        CALL CHECK3B(QCI(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,LM,
     &       SUBR,'qci   ')
#ifdef BLK_2MOM
#endif

        DO J=J_0,J_1
        DO I=I_0,I_1
          IF (Q(I,J,1).gt.1d-1)print*,SUBR," Q BIG ",i,j,Q(I,J,1:LS1)
          IF (T(I,J,1)*PK(1,I,J)-TF.gt.50.) print*,SUBR," T BIG ",i,j
     *         ,T(I,J,1:LS1)*PK(1:LS1,I,J)-TF
        END DO
        END DO
        DO L=1,LM
        DO J=J_0,J_1
        DO I=I_0,I_1
          IF (Q(I,J,L).lt.0.) then
            print*,"After ",SUBR," Q < 0 ",i,j,Q(I,J,L)
            call stop_model('Q<0 in CHECKT',255)
          END IF
          IF (QCL(I,J,L)+QCI(I,J,L).lt.0.) then
         print*,"After ",SUBR," QCL+QCI < 0 ",i,j,QCL(I,J,L),QCI(I,J,L)
            call stop_model('WM<0 in CHECKT',255)
          END IF
#ifdef BLK_2MOM
#endif
        END DO
        END DO
        END DO
C**** Check PBL arrays
        CALL CHECKPBL(SUBR)
C**** Check Ocean arrays
         If (KOCEAN /= 0)  Call CHECKO (SUBR)
C**** Check Ice arrays
        CALL CHECKI(SUBR)
C**** Check Lake arrays
        CALL CHECKL(SUBR)
C**** Check Earth arrays
        CALL CHECKE(SUBR)
C**** Check Land Ice arrays
        CALL CHECKLI(SUBR)
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
C**** check tracers
        CALL CHECKTR(SUBR)
#endif
      EndIf  !  QCHECK

      RETURN
      END SUBROUTINE CHECKT


      subroutine read_aic
!@sum read_AIC for a cold start, read the atmospheric IC file.
!@+   Three input options are currently recognized
!@+      (1) The input file has already been remapped to the model layering
!@+          and contains the variables traditionally expected by the model
!@+          (winds, temperature, specific humidity, surface pressure).
!@+      (2) Winds, temperature, geopotential height, and RH are available
!@+          on constant-pressure levels.  Surface pressure is obtained via
!@+          the hydrostatic assumption.  Remapping from constant-pressure
!@+          levels to the model layering is peformed, and RH is converted
!@+          to specific humidity. The initial coding for this option was
!@+          imported from init_cond/AIC.D771201.f and generalized considerably.
!@+          If tropopause pressure and temperature are available, they
!@+          are used in the vertical remapping.  The coding for the
!@+          assumed vertical structure of RH above the upper troposphere
!@+          will be made more configurable as AIC files with better-quality
!@+          RH data (than the historical AIC) become available.
!@+      (3) Surface pressure, winds, temperature, humidity on "standard-hybrid"
!@+          layers, where the pressure of the kth input layer at horizontal
!@+          location i,j is equal to hyam(k) + hybm(k)*srfp(i,j) and srfp is
!@+          the surface pressure. Remapping to model layering is performed.
!@+          Condensate species are also (optionally) read on other branches,
!@+          but not here (yet), for reasons including the inability of liquid
!@+          and ice to coexist under the cloud microphysics on this branch.
!@+
!@+   Logic will be added to handle other possible input layerings and
!@+   combinations of available fields.
      use constant, only : grav,rgas,lhe
      use VerticalRes, only : lm
      use atm_com, only : zatmo,psrf=>p,
     &     ualij,valij,uout=>u,vout=>v,
     &     tout=>t,qout=>q
      use atm_com, only : traditional_coldstart_aic
      use fluxes, only : atmsrf
      use pario, only : par_open,par_close,read_data,read_dist_data
     &     ,get_dimlen,get_dimlens,variable_exists
      Use DOMAIN_DECOMP_ATM, Only: GRID, GetDomainBounds,
     &     am_i_root,halo_update_column,hasnorthpole,hassouthpole
      implicit none

      real*8, dimension(:,:,:), allocatable :: uin,vin,tin

      real*8, dimension(:), allocatable :: u,v,t,p,plev
      real*8 :: pe(0:lm)
      real*8, dimension(lm) :: pmid, MAdum,PDSIGdum

      INTEGER :: I,J,L,K,N,KM
      integer :: dlens(7)
      integer :: i_0h,i_1h,j_0h,j_1h
      integer :: i_0,i_1,j_0,j_1
      integer :: j_0stg,j_1stg
      logical :: hybridlayer_input
      integer :: fid

      Call GetDomainBounds(GRID,
     &     I_STRT=I_0,I_STOP=I_1,
     &     J_STRT=J_0,J_STOP=J_1)

      Call GetDomainBounds(GRID,
     &     I_STRT_HALO=I_0H,I_STOP_HALO=I_1H,
     &     J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      Call GetDomainBounds(GRID,
     &     J_STRT_STGR=J_0STG,J_STOP_STGR=J_1STG)

      fid = par_open(grid,'AIC','read')

      traditional_coldstart_aic = 
     &           variable_exists(grid,fid,'u')
     &     .and. variable_exists(grid,fid,'v')
     &     .and. variable_exists(grid,fid,'t')
     &     .and. variable_exists(grid,fid,'p')
     &     .and. variable_exists(grid,fid,'q')


      hybridlayer_input =
     &           variable_exists(grid,fid,'hyam')
     &     .and. variable_exists(grid,fid,'ps_a') ! change this name soon

      if(traditional_coldstart_aic) then ! option 1 in subr. header comments
        call read_dist_data(grid,fid,'u',uout)
        call read_dist_data(grid,fid,'v',vout)
        call read_dist_data(grid,fid,'t',tout)
        call read_dist_data(grid,fid,'p',psrf)
        call read_dist_data(grid,fid,'q',qout)
      elseif(hybridlayer_input) then     ! option 3 in subr. header comments
        call relayer_hybridlayer_input
      else                               ! option 2 in subr. header comments
        call read_pzrh_input
      endif

      call par_close(grid,fid)

      contains

      subroutine relayer_hybridlayer_input  ! option 3 in subr. header comments
!@sum relayer_hybridlayer_input read and relayer atm state on standard hybrid layers
!@+   to whatever the model layering is.   Relayering currently performed by vint_logp
!@+   for historical continuity with hindcasting runs that performed the relayering
!@+   with external codes; this means that the relayering is not conservative.
!@+   Other continuity choices:
!@+   - reading u,v at B-grid locations (and aux B-grid surface pressure for relayering)
!@+   - vertically interpolate log(q) rather than q.
      use resolution, only : im
      implicit none
      real*8, dimension(:,:,:), allocatable :: qin
      real*8, dimension(:,:), allocatable :: psrfb
      real*8, dimension(:), allocatable :: hyam,hybm,xin,xout
      real*8 :: byim

      byim = 1d0/real(im,kind=8)
      km = get_dimlen(grid,fid,'level')

      allocate(hyam(km),hybm(km),p(km),xin(km),xout(lm))
      allocate(
     &     uin(i_0h:i_1h,j_0h:j_1h,km),
     &     vin(i_0h:i_1h,j_0h:j_1h,km),
     &     tin(i_0h:i_1h,j_0h:j_1h,km),
     &     qin(i_0h:i_1h,j_0h:j_1h,km),
     &     psrfb(i_0h:i_1h,j_0h:j_1h)
     &     )

      call read_data(grid,fid,'hyam',hyam,bcast_all=.true.)
      call read_data(grid,fid,'hybm',hybm,bcast_all=.true.)

      call read_dist_data(grid,fid,'ps_a',psrf)
      call read_dist_data(grid,fid,'ps_b',psrfb)
      call read_dist_data(grid,fid,'u',uin)
      call read_dist_data(grid,fid,'v',vin)
      call read_dist_data(grid,fid,'t',tin)
      call read_dist_data(grid,fid,'q',qin)

      hyam  = hyam*.01d0  ! Pa -> hPa
      psrf  = psrf*.01d0  ! Pa -> hPa
      psrfb = psrfb*.01d0 ! Pa -> hPa

      if(hasSouthPole(grid)) then
        psrf(:,j_0) = sum(psrf(:,j_0))*byim
        do k=1,km
          tin(:,j_0,k) = sum(tin(:,j_0,k))*byim
          qin(:,j_0,k) = sum(qin(:,j_0,k))*byim
        enddo
      endif
      if(hasNorthPole(grid)) then
        psrf(:,j_1) = sum(psrf(:,j_1))*byim
        do k=1,km
          tin(:,j_1,k) = sum(tin(:,j_1,k))*byim
          qin(:,j_1,k) = sum(qin(:,j_1,k))*byim
        enddo
      endif

C****
C**** Perform vertical interpolation from input pressures to model pressures.
C****

      do j=j_0,j_1
      do i=i_0,i_1
        call calc_vert_amp(psrf(i,j),lm, madum,pdsigdum,pe,pmid)
        p = hyam + hybm*psrf(i,j)
        ! temperature
        xin = tin(i,j,:)
        call vint_logp(km,lm,p,pmid,xin,xout)
        tout(i,j,:) = xout
        ! humidity
        xin = qin(i,j,:)
        xin = log(max(xin,1d-7)) ! to interpolate log(q)
        call vint_logp(km,lm,p,pmid,xin,xout)
        xout = exp(xout)         ! to interpolate log(q)
        qout(i,j,:) = xout
      enddo
      enddo

      do j=j_0stg,j_1stg
      do i=i_0,i_1
        call calc_vert_amp(psrfb(i,j),lm, madum,pdsigdum,pe,pmid)
        p = hyam + hybm*psrfb(i,j)
        ! u-wind
        xin = uin(i,j,:)
        call vint_logp(km,lm,p,pmid,xin,xout)
        uout(i,j,:) = xout
        ! v-wind
        xin = vin(i,j,:)
        call vint_logp(km,lm,p,pmid,xin,xout)
        vout(i,j,:) = xout
      enddo
      enddo
      call recalc_agrid_uv

      atmsrf%TSAVG = tout(:,:,1) ! improved estimate during PBL init

      end subroutine relayer_hybridlayer_input

      subroutine read_pzrh_input  ! option 2 in subr. header comments
      Use DOMAIN_DECOMP_ATM, Only: globalmax
      use resolution, only : im
      implicit none
      real*8, dimension(:), allocatable :: rh
      real*8, dimension(:,:,:), allocatable :: zin,rh1,rhin
      real*8, dimension(:,:), allocatable :: ptrop,ttrop,zsrf
      INTEGER :: K1,KK,KMZIN,KMTROP,KMRH,KMIN
      real*8, dimension(lm) :: xa,xb
      REAL*8 HSRF,TM,PR,PL,DTDZ,RHTROP,WTDN,BYIM
      real*8 :: max_loc,max_zsrf,max_zin,max_ptrop,max_psrat
      real*8, parameter :: GBYR  = GRAV/RGAS

      REAL*8 QSAT ! external function

      UALIJ = 0.
      VALIJ = 0.

      call get_dimlens(grid,fid,'p',n,dlens)
      if(n.ne.1) call stop_model('read_aic: bad ndims for p',255)
      km = dlens(1)

      call get_dimlens(grid,fid,'z',n,dlens)
      if(n.ne.3) call stop_model('read_aic: bad ndims for z',255)
      kmzin = dlens(3)

      call get_dimlens(grid,fid,'rh',n,dlens)
      if(n.ne.3) call stop_model('read_aic: bad ndims for rh',255)
      kmrh = dlens(3)

      allocate(u(0:km+1),v(0:km+1),t(0:km+1),rh(0:km+1),p(0:km+1))
      allocate(plev(km))
      allocate(
     &     uin(i_0h:i_1h,j_0h:j_1h,km),
     &     vin(i_0h:i_1h,j_0h:j_1h,km),
     &     tin(i_0h:i_1h,j_0h:j_1h,km),
     &     rhin(i_0h:i_1h,j_0h:j_1h,km),
     &     zin(i_0h:i_1h,j_0h:j_1h,kmzin),
     &     rh1(i_0h:i_1h,j_0h:j_1h,kmrh))
      allocate(
     &     zsrf(i_0h:i_1h,j_0h:j_1h),
     &     ptrop(i_0h:i_1h,j_0h:j_1h),
     &     ttrop(i_0h:i_1h,j_0h:j_1h))

      call read_data(grid,fid,'p',plev,bcast_all=.true.)
      call read_dist_data(grid,fid,'u',uin)
      call read_dist_data(grid,fid,'v',vin)
      call read_dist_data(grid,fid,'t',tin)
      call read_dist_data(grid,fid,'rhclim',rhin)
      call read_dist_data(grid,fid,'rh',rh1)
      call read_dist_data(grid,fid,'z',zin)
C****   PTROP = tropopause pressure (mb)
C****   TTROP = tropopause temperature (K)
      ptrop = 0. ! if ptrop not present, it is assumed to be zero.
      call read_dist_data(grid,fid,'ptrop',ptrop)
      call read_dist_data(grid,fid,'ttrop',ttrop)

      do k=1,kmrh
        rhin(:,:,k) = rh1(:,:,k)
      enddo

      zsrf(i_0:i_1,j_0:j_1) = zatmo(i_0:i_1,j_0:j_1)/grav

      max_loc = maxval(zsrf(i_0:i_1,j_0:j_1))
      call globalmax(grid,max_loc,max_zsrf)
      max_loc = maxval(zin(i_0:i_1,j_0:j_1,kmzin))
      call globalmax(grid,max_loc,max_zin)
      if(max_zsrf.gt.max_zin) then
        if(am_i_root()) then
          write(6,*) 'Surface topography exceeds largest input z'
          write(6,*) 'max(zsrf), max(zin) = ',max_zsrf,max_zin
        endif
        call stop_model('error in AIC vert. interp.',255)
      endif

      max_loc = maxval(ptrop(i_0:i_1,j_0:j_1))
      call globalmax(grid,max_loc,max_ptrop)
      if(max_ptrop.ge.plev(1)) then
        if(am_i_root()) then
          write(6,*) 'Tropopause pressure exceeds first pressure level.'
          write(6,*) 'max(PTROP),PLEV(1) ',max_ptrop,PLEV(1)
        endif
        call stop_model('error in AIC vert. interp.',255)
      endif

C****
C**** Perform vertical interpolation of prognostic quantities
C**** from pressure surfaces to model levels.
C****

      DO J=J_0,J_1
      DO I=I_0,I_1

C**** Determine lowest input pressure level above the surface
      HSRF = ZSRF(I,J)
      DO K1=1,KMZIN
        IF(ZIN(I,J,K1).gt.HSRF) exit
      enddo

C**** Determine highest input pressure level below the tropopause
      DO KMTROP=KM,1,-1
        IF(PLEV(KMTROP).gt.PTROP(I,J)) exit
      enddo

      !CALL TSBOOK(JDAY,PTROP(I,J),DLATDG*(J-(.5*(1+JM))),TTRPPC,QTROP)
      RHTROP = 0. ! 100.*QTROP/QSAT(TTRPPC,LHE,PTROP(I,J))

c load pressure, temperature, and rh values into 1D arrays.

      kk = 0
      DO K=K1,KM
        if(k.eq.kmtrop+1) then
c insert tropopause values if available.
          kk = kk + 1
          P(KK) = PTROP(I,J)
          T(KK) = TTROP(I,J)
          RH(KK) = RHTROP
        endif
        kk = kk + 1
        P(kk) = PLEV(K)
        T(kk) = TIN(I,J,K)
        if(k.gt.kmtrop) then
C**** To put Q=0 in stratosphere, set R... = 0. below
          RH(kk) = 0.           ! RHJ(J,K)
        elseif(k.le.kmrh)  then
          RH(kk) = RHIN(I,J,K)
        else
          wtdn = (PLEV(K)-PTROP(I,J))/(PLEV(KMRH)-PTROP(I,J))
          RH(kk) = RHIN(I,J,KMRH)*wtdn + RHTROP*(1d0-wtdn)
        endif
      enddo

      KMIN = KK

C**** Calculate surface pressure
      IF(K1.eq.1) THEN
        DTDZ = 0D0
      ELSE
        DTDZ = (TIN(I,J,K1)-TIN(I,J,K1-1)) / (ZIN(I,J,K1)-ZIN(I,J,K1-1))
      ENDIF
      IF(ABS(DTDZ).lt.1.E-5) THEN  ! isothermal
        P(0) = PLEV(K1)
     &       *EXP(GBYR*(ZIN(I,J,K1)-HSRF)/T(1))
      ELSE
        P(0) = PLEV(K1)
     &       *(1.-DTDZ*(ZIN(I,J,K1)-HSRF)/T(1))**(-GBYR/DTDZ)
      ENDIF

C**** Calculate surface air temperature and surface relative humidity
      IF(K1.gt.1) THEN
        wtdn = (P(0)-P(1))/(PLEV(K1-1)-P(1))
        T(0) = wtdn*TIN(I,J,K1-1) + (1d0-wtdn)*T(1)
        RH(0) = wtdn*RHIN(I,J,K1-1) + (1d0-wtdn)*RH(1)
      else
        T(0) = TIN(I,J,1)
        RH(0) = RHIN(I,J,1)
      endif

C****
C**** Calculate output values
C****
C**** Surface pressure, surface air temperature
C****
      PSRF(I,J)  = P(0)
      atmsrf%TSAVG(I,J) = T(0)

      Call CALC_VERT_AMP (P(0),LM, MAdum,PDSIGdum,PE,PMID)

C****
C**** Remap temperature and RH to model layers
C**** Convert RH to specific humidity
C****
      CALL VNTRP1 (KMIN,P,T, LM,PE,XA)
      CALL VNTRP1 (KMIN,P,RH, LM,PE,XB)
      DO L=1,LM
        TOUT(I,J,L) = XA(L)
        QOUT(I,J,L) = max(3.d-6,.01*XB(L)*QSAT(XA(L),LHE,PMID(L)))
      ENDDO

c A-grid winds.  No insertion of tropopause values.
      U(0) = UIN(I,J,max(1,K1-1))
      V(0) = VIN(I,J,max(1,K1-1))
      P(0) = PLEV(max(1,K1-1)) + 1d-6
      kk = 0
      DO K=K1,KM
        kk = kk + 1
        P(kk)=PLEV(K)
        U(kk)=UIN(I,J,K)
        V(kk)=VIN(I,J,K)
      ENDDO
      KMIN = KK
      CALL VNTRP1 (KMIN,P,U, LM,PE,XA)
      CALL VNTRP1 (KMIN,P,V, LM,PE,XB)
      DO L=1,LM
        UALIJ(L,I,J) = XA(L)
        VALIJ(L,I,J) = XB(L)
      ENDDO

      ENDDO ! I
      ENDDO ! J


      byim = 1d0/real(im,kind=8)
      if(hasSouthPole(grid)) then
        psrf(:,j_0) = sum(psrf(:,j_0))*byim
        do k=1,lm
          tout(:,j_0,k) = sum(tout(:,j_0,k))*byim
          qout(:,j_0,k) = sum(qout(:,j_0,k))*byim
        enddo
      endif
      if(hasNorthPole(grid)) then
        psrf(:,j_1) = sum(psrf(:,j_1))*byim
        do k=1,lm
          tout(:,j_1,k) = sum(tout(:,j_1,k))*byim
          qout(:,j_1,k) = sum(qout(:,j_1,k))*byim
        enddo
      endif

      max_loc = maxval(ptrop(i_0:i_1,j_0:j_1)/psrf(i_0:i_1,j_0:j_1))
      call globalmax(grid,max_loc,max_psrat)
      if(max_psrat.gt..9d0) then
        if(am_i_root()) then
          write(6,*) 'Tropopause pressure too high'
          write(6,*) 'max(ptrop/psrf) = ',max_psrat
        endif
        call stop_model('error in AIC vert. interp.',255)
      endif

#if defined(SCM) || defined(CUBED_SPHERE)
      do l=1,lm
        do j=j_0,j_1
          do i=i_0,i_1
            uout(i,j,l) = ualij(l,i,j)
            vout(i,j,l) = valij(l,i,j)
          enddo
        enddo
      enddo
#else
      ! convert velocities to their native grid/orientation
      ! for the moment, this only applies to the B-grid
      ! dynamics scheme.  The cubed-sphere case is handled elsewhere.
      call halo_update_column(grid,ualij)
      call halo_update_column(grid,valij)
      do l=1,lm
        do j=j_0stg,j_1stg
          do i=i_0,i_1-1 ! im-1
            uout(i,j,l) = .25d0*
     &           (ualij(l,i,j-1)+ualij(l,i+1,j-1)
     &           +ualij(l,i,j  )+ualij(l,i+1,j  ))
            vout(i,j,l) = .25d0*
     &           (valij(l,i,j-1)+valij(l,i+1,j-1)
     &           +valij(l,i,j  )+valij(l,i+1,j  ))
          enddo
          i = i_1 ! im
            uout(i,j,l) = .25d0*
     &           (ualij(l,i,j-1)+ualij(l,  1,j-1)
     &           +ualij(l,i,j  )+ualij(l,  1,j  ))
            vout(i,j,l) = .25d0*
     &           (valij(l,i,j-1)+valij(l,  1,j-1)
     &           +valij(l,i,j  )+valij(l,  1,j  ))
        enddo
      enddo
#endif

      end subroutine read_pzrh_input

      SUBROUTINE VNTRP1 (KM,P,AIN,  LMA,PE,AOUT)
C**** Vertically interpolates a 1-D array
C**** Input:       KM = number of input pressure levels
C****            P(K) = input pressure levels (mb)
C****          AIN(K) = input quantity at level P(K)
C****             LMA = number of vertical layers of output grid
C****           PE(L) = output pressure levels (mb) (edges of layers)
C**** Output: AOUT(L) = output quantity: mean between PE(L-1) & PE(L)
C****
      IMPLICIT NONE
      INTEGER :: KM,LMA
      REAL*8 P(0:KM),AIN(0:KM),    PE(0:LMA),AOUT(LMA)
      INTEGER :: K,K1,L
      REAL*8 :: PDN,ADN,PUP,AUP,PSUM,ASUM
C****
      PDN = PE(0)
      ADN = AIN(0)
      K=1
C**** Ignore input levels below ground level pe(0)=p(0)
      IF(P(1).GT.PE(0)) THEN
         DO K1=2,KM
         K=K1
         IF(P(K).LT.PE(0)) THEN  ! interpolate to ground level
           ADN=AIN(K)+(AIN(K-1)-AIN(K))*(PDN-P(K))/(P(K-1)-P(K))
           GO TO 300
         END IF
         END DO
         STOP 'VNTRP1 - error - should not get here'
      END IF
C**** Integrate - connecting input data by straight lines
  300 DO 330 L=1,LMA
      ASUM = 0.
      PSUM = 0.
      PUP = PE(L)
  310 IF(P(K).le.PUP)  GO TO 320
      PSUM = PSUM + (PDN-P(K))
      ASUM = ASUM + (PDN-P(K))*(ADN+AIN(K))/2.
      PDN  = P(K)
      ADN  = AIN(K)
      K=K+1
      IF(K.LE.KM) GO TO 310
      stop 'VNTRP1 - should not happen'
C****
  320 AUP  = AIN(K) + (ADN-AIN(K))*(PUP-P(K))/(PDN-P(K))
      PSUM = PSUM + (PDN-PUP)
      ASUM = ASUM + (PDN-PUP)*(ADN+AUP)/2.
      AOUT(L) = ASUM/PSUM
      PDN = PUP
  330 ADN = AUP
C****
      RETURN
      END SUBROUTINE VNTRP1

      subroutine vint_logp(kmi,kmo,pi,po,input,output)
!@sum vint_logp linear interpolation in log(p)-space
      implicit none
      integer, intent(in):: kmi,kmo ! number of input,output levels
      real*8, dimension(kmi), intent(in) :: pi ! input pressures
      real*8, dimension(kmo), intent(in) :: po ! output pressures
      real*8, dimension(kmi), intent(in) :: input
      real*8, dimension(kmo), intent(out) :: output
      !real*8, dimension(:), allocatable :: lnpi,lnpo
      real*8 :: wtdn
      integer:: ki, ko, kolin1

      !allocate(lnpi(kmi),lnpo(kmo))
      !lnpi=log(pi)
      !lnpo=log(po)

      do kolin1=1,kmo
        if(po(kolin1).le.pi(1)) exit
      enddo
      do ko=1,kolin1-1
        output(ko)=input(1)
      enddo
      !if(kolin1.gt.1) write(6,*) 'KOLIN1 > 1'

      ki = 1
      do ko=kolin1,kmo
        do while(pi(ki+1).gt.po(ko))
          ki = ki + 1
        enddo
        !wtdn=(lnpo(ko)-lnpi(ki+1))/(lnpi(ki)-lnpi(ki+1))
        wtdn = log(po(ko)/pi(ki+1))/log(pi(ki)/pi(ki+1))
        output(ko)=wtdn*input(ki)+(1d0-wtdn)*input(ki+1)
      enddo
      return
      end subroutine vint_logp

      end subroutine read_aic

#ifdef CACHED_SUBDD
      subroutine accum_subdd_atm
C**** interpolate to pressure levels and accumulate the subdd diagnostics
      USE CONSTANT, only : teeny,lhe,lhm,sha,bygrav
      use subdd_mod, only : lmaxsubdd
      use subdd_mod, only : subdd_type,subdd_groups,subdd_ngroups
      use subdd_mod, only : aijph_l1,aijph_l2
     &      ,subdd_npres,subdd_pk, subdd_pres
      use subdd_mod, only : inc_subdd,find_groups
      use atm_com,    only: u,v,t,q,qcl,qci, pdsig,pmid,pedn,pk,
     &                      ualij,valij, zatmo,gz, wsave, ma,masum
     &                     ,ptropo,ltropo
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use resolution, only : lm,mtop
      USE GEOM, only: imaxj
      use fluxes, only : atmsrf,atmice
      use model_com, only : dtsrc
      implicit none

      INTEGER :: LDN,LUP,I,J,L,k,igrp,ngroups,grpids(subdd_ngroups)
      INTEGER :: J_0, J_1, J_0H, J_1H, I_0,I_1
      type(subdd_type), pointer :: subdd
      REAL*8 QSAT, WTDN,WTUP,qinterp,tinterp,qsat_interp
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     vortl,sddarr
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     sddarr2d
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,subdd_npres)
     &     :: sddcp

      real*8 slp ! function

      CALL GET(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               J_STRT=J_0,        J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      call find_groups('aijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
C
      case ('p_surf')
        sddarr2d = pedn(1,:,:)
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('z_surf')
        call inc_subdd(subdd,k,zatmo)
C
      case ('pblht')
        call inc_subdd(subdd,k,atmsrf%dblavg)
C
      case ('shflx')
        sddarr2d = atmsrf%sensht(:,:)/dtsrc
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('lhflx')
        sddarr2d = atmsrf%latht(:,:)/dtsrc
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('pwv')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = sum(q(i,j,1:LM)*pdsig(1:LM,i,j))*100.*bygrav
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
C
C     East-west humidity flux (vert sum)
C     Calculated on primary (A) grid
      case ('puq')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = sum(ualij(1:LM,i,j)*q(i,j,1:LM)*
     &                    pdsig(1:LM,i,j))*100.*bygrav
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
C
C     North-south humidity flux (vert sum)
C     Calculated on primary (A) grid
      case ('pvq')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = sum(valij(1:LM,i,j)*q(i,j,1:LM)*
     &                    pdsig(1:LM,i,j))*100.*bygrav
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d) 
C
      case ('lwp')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = sum(qcl(i,j,1:LM)*pdsig(1:LM,i,j))*100.*bygrav
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('iwp')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = sum(qci(i,j,1:LM)*pdsig(1:LM,i,j))*100.*bygrav
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('column_fmse')
        ! same as total static energy computed by conserv_se,
        ! which is disabled for SCM
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = sum( 
     &           (sha*t(i,j,1:LM)*pk(1:LM,i,j)+
     &            lhe*q(i,j,1:LM)-lhm*qci(i,j,1:LM))*ma(1:LM,i,j))+
     &           zatmo(i,j)*(masum(i,j)+mtop)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('snowdp')
        call inc_subdd(subdd,k,atmice%snowsave)
C
      case ('slp')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          !ts = t(i,j,1)*pek(1,i,j)
          sddarr2d(i,j) =
     &         slp(pedn(1,i,j),atmsrf%tsavg(i,j),bygrav*zatmo(i,j))
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('ptrop')
        call inc_subdd(subdd,k,ptropo)
C
      case ('ttrop')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = t(i,j,ltropo(i,j))*pk(ltropo(i,j),i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
      end select
      enddo
      enddo

      call find_groups('aijph',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('ucp')
        call inc_subdd(subdd,k,ualij,jdim=3)
      case ('vcp')
        call inc_subdd(subdd,k,valij,jdim=3)
      case ('zcp')
        call inc_subdd(subdd,k,gz)
      case ('tcp')
        do l=1,subdd_npres; do j=j_0,j_1; do i=i_0,imaxj(j)
           ldn = aijph_l1(i,j,l)    ;  lup = aijph_l2(i,j,l)
          wtdn = aijph_l1(i,j,l)-ldn; wtup = aijph_l2(i,j,l)-lup
          sddcp(i,j,l) =
     &         (wtdn*t(i,j,ldn) +wtup*t(i,j,lup))*subdd_pk(l)
        enddo;              enddo;        enddo
        call inc_subdd(subdd,k,sddcp)
      case ('qcp')
        do l=1,subdd_npres; do j=j_0,j_1; do i=i_0,imaxj(j)
           ldn = aijph_l1(i,j,l)    ;  lup = aijph_l2(i,j,l)
          wtdn = aijph_l1(i,j,l)-ldn; wtup = aijph_l2(i,j,l)-lup
          if(wtdn+wtup.gt.0.) then
            sddcp(i,j,l) =
     &       exp(wtdn*log(q(i,j,ldn)+teeny) +wtup*log(q(i,j,lup)+teeny))
          else
            sddcp(i,j,l) = 0.
          endif
        enddo;              enddo;        enddo
        call inc_subdd(subdd,k,sddcp)
      case ('vortcp')
#ifdef SCM
        call stop_model('no SUBDD vorticity in SCM mode',255)
#else
        call get_vorticity(vortl)
        call inc_subdd(subdd,k,vortl)
#endif
      case ('wcp')
        sddarr(:,:,1:lm-1) = wsave
        sddarr(:,:,lm) = 0.
        call inc_subdd(subdd,k,sddarr)
      case ('rhcp')
        do l=1,subdd_npres; do j=j_0,j_1; do i=i_0,imaxj(j)
           ldn = aijph_l1(i,j,l)    ;  lup = aijph_l2(i,j,l)
          wtdn = aijph_l1(i,j,l)-ldn; wtup = aijph_l2(i,j,l)-lup
C**** The if-statement was added for nan-locations near topography, surface. 
C**** "else qinterp=0" should not impact results. 
          if(wtdn+wtup.gt.0.) then
             qinterp =
     &     exp(wtdn*log(q(i,j,ldn)+teeny)+wtup*log(q(i,j,lup)+teeny))
          else
             qinterp = 0.0d0
          endif
          tinterp=(wtdn*t(i,j,ldn) +wtup*t(i,j,lup))*subdd_pk(l)
          qsat_interp = QSAT(tinterp,LHE,subdd_pres(l))
          sddcp(i,j,l) = (qinterp/qsat_interp)*100.0d0
         enddo;              enddo;        enddo
         call inc_subdd(subdd,k,sddcp)
      end select
      enddo
      enddo

C**** cached_subdd on model levels
      call find_groups('aijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('t')
        do l=1,lmaxsubdd; do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j,l) = t(i,j,l)*pk(l,i,j)
        enddo;              enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('th')
        call inc_subdd(subdd,k,t)
      case ('q')
        call inc_subdd(subdd,k,q)
      case ('rh')
        do l=1,lmaxsubdd; do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr(i,j,l) = 
     &    q(i,j,l)/QSAT(t(i,j,l)*pk(l,i,j),LHE,pmid(l,i,j))*100.0d0
        enddo;              enddo;        enddo
        call inc_subdd(subdd,k,sddarr)
      case ('z')
        call inc_subdd(subdd,k,gz)
      case ('p_3d')
        call inc_subdd(subdd,k,pmid,jdim=3)
      case ('u')
        call inc_subdd(subdd,k,ualij,jdim=3)
      case ('v')
        call inc_subdd(subdd,k,valij,jdim=3)
#ifndef CUBED_SPHERE
      case ('ub') ! b-grid wind, EW component
        call inc_subdd(subdd,k,u)
      case ('vb') ! b-grid wind, NS component
        call inc_subdd(subdd,k,v)
#endif
      case ('w')
        sddarr(:,:,1:lm-1) = wsave
        sddarr(:,:,lm) = 0.
        call inc_subdd(subdd,k,sddarr)
      end select
      enddo
      enddo
      return
      end subroutine accum_subdd_atm
#endif
