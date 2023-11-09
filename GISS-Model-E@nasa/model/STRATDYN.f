#include "rundeck_opts.h"

!@sum  STRATDYN stratospheric only routines
!@auth Bob Suozzo/Jean Lerner/Gavin Schmidt/Jeff Jonas/Maxwell Kelley

C**** TO DO:
C****   i) A-grid <-> B-grid  should be done with indexes etc.

      subroutine get_kep(kep)
      implicit none
      integer :: kep
      kep = 21
      end subroutine get_kep

      MODULE STRAT
!@sum  STRAT local stratospheric variables for GW drag etc.
!@auth Bob Suozzo/Jean Lerner
      USE RESOLUTION, only : im,jm,lm
#ifdef CUBED_SPHERE
      use cs2ll_utils, only : uv_derivs_type
#endif
      IMPLICIT NONE
      SAVE
!@dbparam XCDNST parameters for GW drag (in param. database)
      REAL*8, DIMENSION(2) :: XCDNST(2)
!@dbparam Cshear parameter for GW shear drag (in param. database)
!@dbparam CMTN parameter for GW MTN drag (in param. database)
!@dbparam CDEF parameter for GW DEF drag (in param. database)
!@dbparam CMC parameter for GW M. Convective drag (database if LM<=40)
!@dbparam SCVMU parameter to Scale Convective MUs (database if LM>40)
C**** (used to be FMC)
      REAL*8 :: CMTN=.5, CDEF=3, CMC=2d-7, Cshear=1, SCVMU=.188d0
!@dbparam PBREAK p. level above which GW drag acts (in param. database)
      REAL*8 :: PBREAK = 500.   ! default is 500mb
!@dbparam PCONPEN level of penetrating moist conv (in param. database)
      REAL*8 :: PCONPEN = 400.   ! default is 400mb
!@dbparam DEFTHRESH threshold for deformation wave (1/s)
      REAL*8 :: DEFTHRESH = 15d-6  ! default is 15x10^-6 s^-1
!@dbparam PBREAKTOP p. level to force GW breaking in top layer
C**** This should be set to 100. (or something similar) to force
C**** breaking of remaining gravity waves in top layer. Otherwise,
C**** momentum passes through model top.
      REAL*8 :: PBREAKTOP = 0.05d0   ! default is 0.05mb

!@var XLIMIT per timestep limit on mixing and drag
      REAL*8 :: XLIMIT=.1d0  ! default

!@var ZVART,ZVARX,ZVARY,ZWT topographic variance
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ZVART,ZVARX,ZVARY,ZWT
!@var DEFRM deformation field
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DEFRM
!@var LDEF,LDEFM deformation levels
      INTEGER LDEF ! ,LDEFM
!@var LBREAK,LSHR,LD2 levels for various GW drag terms
      INTEGER :: LBREAK,LSHR,LD2 = LM   ! need default for LD2
      INTEGER :: LPCNV     !  200mb level, for Convective Waves  

!@dbparam QGWMTN =1 turns on GW Mountain Wave drag terms
!@dbparam QGWSHR =1 turns on GW Shear drag terms
!@dbparam QGWDEF =1 turns on GW Deformation drag terms
!@dbparam QGWCNV =1 turns on GW Convective drag terms
      INTEGER :: QGWMTN = 1, QGWSHR = 1, QGWDEF = 1, QGWCNV = 1

!@dbparam ang_gwd =1 ang mom. lost by GWDRAG is added in troposphere
      INTEGER :: ang_gwd = 1 ! default: GWDRAG does conserve AM

!@param NM number of gravity wave drag sources
      INTEGER, PARAMETER :: NM=15 ! = mtn, deform, shear + 2x6 conv. waves
!@var EKOFJ wavenumbers as a function of GW source type and latitude J
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: EKOFJ ! (NM,J:)

c
c Inputs/outputs for the column gravity wave drag routine
c
      REAL*8, DIMENSION(LM) ::         ! for each layer:
c pressure, press. thick., temp., pot. temp., density, B-V freq, winds
     &     TL,THL,RHO,BVF,UL,VL, ! input
c GW diffusivity, total GWD wind changes
     &     DL, DUT,DVT,DUSDIF          ! output

      REAL*8, DIMENSION(NM) ::         ! indexed by GW source type
c average wavenumber
     &     EK,                         ! input
c mom. flx at src. level and model top, wave dir. unit vector, phase speed
     &     MU_INC,MU_TOP, UR,VR, CN,   ! output
c area weight
     &     WT                          ! output, but wt(1) is input

c GWD wind change profiles indexed by wave source type
      REAL*8  :: DUGWD(LM,NM)          ! output

      REAL*8 ::
c timestep+reciprocal, coriolis parameter, random number, MC mass flux
     &     DTIME,BYDTIME, CORIOL, RANMTN, AIRXS,                   ! input
c x-/y- topog. variance + weight, surf. geopot., 700 mb deformation
     &     ZVARX_cell,ZVARY_cell,ZWT_cell, ZATMO_cell, DEFRM_cell, ! input
c wind components averaged over MC wave generation levels
     &     USRC, VSRC                                              ! output

      INTEGER ::
c base/top layer indices of moist convection
     &     LMC0,LMC1,                  ! input
c lowest level at which there are gravity waves
     &     LDRAG                       ! output

#ifdef CUBED_SPHERE
      type(uv_derivs_type) :: dfm_type ! used for calculating defrm
      real*8 :: EK_globavg ! no horz variation of EK (quasi-uniform gridsize)
#endif
      EndModule STRAT


      Subroutine GWDCOL (PLE,PL,MAUV)
C 
C     GWDCOL  Column Gravity Wave Drag Routine
C 
C****
C**** LD: incident layer - lowest layer at which wave is allowed to
C****   break - set to LM+1 if wave doesn't exist
C**** MU: momentum flux of wave - conserved if the wave doesn't break
C**** MUB: breaking momentum flux - proportional to -(W-CN)**3
C****    MU is the Greek letter mu, for momentum.
C****    MU and MUB must be negative.
C**** CN: wave phase speed (m/s)
C**** W-CN: (local wind in direction of wave) - (phase speed of wave)
C****    calculated as WMC = (U*UR + V*VR - CN)  must be positive!!!
C**** For propagation, we need MUB < MU < 0. When MUB > 0, the wave
C****    cannot propagate, instead, it breaks and deposits all its
C****    momentum.  The level where MUB changes sign is a critical
C****    level.  W-CN also changes sign.
C**** (UR,VR): unit vector in direction of the wave
C**** There is some trickiness in the way (W-CN), MU, and UR,VR
C**** are gotten for MC waves. (UR,VR) is in the direction of the mean
C**** wind in the source region. But adding +/- 10, 20, and 40 m/s means
C**** sometimes W-CN is negative, unless the direction of (UR,VR) is
C**** changed.  7159-7159.9 fixes this so that (UR,VR) agrees with CN.
C**** NM:
C**** 1   Mountain waves
C**** 2   Shear waves
C**** 3-8 Convective waves
C**** 9   Deformation wave
C**** 10-15 Additional conv. waves.  Todo: all conv. waves in contiguous index range
C****
C 
      Use CONSTANT, Only: GRAV,RGAS,kg2mb 
      Use STRAT
      IMPLICIT NONE 
C 
C     *******  Declaration Section  *******  
C
      REAL*8, PARAMETER :: ERR=1d-20, H0=8000., XFROUD=1.  
      REAL*8, PARAMETER :: ROTK = 1.5, RKBY3= ROTK*ROTK*ROTK   
      Real*8,Intent(In) :: PLE(LM+1),PL(LM),MAUV(LM)
C
C     **  LOCAL  **   
C
      Real*8  :: DP(LM), MUB(LM+1,NM),MU(NM)
      REAL*8 :: MU3HOLD, MU37HOLD 
      REAL*8, DIMENSION(LM) :: UEDGE,VEDGE,BYFACS,WMC,DFM,DFR,DFTL
      INTEGER LD(NM)
      INTEGER L,N,LN,NMX,LD1,LTOP,LMAX   
      REAL*8  U0,V0,W0,BV0,ZVAR,P0,DU,DV,DW,CU,CLDDEP,   
     *        FPLUME,CLDHT,WTX,TEDGE,WMCE,BVEDGE,DFMAX,EXCESS,ALFA,    
     *        XDIFF,DFT,DWT,FDEFRM,WSRC,WCHECK,DUTN,PDN,YDN,FLUXUD,  
     *        FLUXVD,PUP,YUP,DX,DLIMIT,FLUXU,FLUXV,MDN,MUP,MUR
C****
C**** INITIALIZE THE MOMENTUM FLUX FOR VARIOUS WAVES
C****
      DP(:)     = MAUV(:)*kg2mb
      MU(:)         = 0.
      UR(:)         = 0
      LD(:)         = LM+1    ! set drag level = LM+1 (turns off drag)
      WT(:)         = 1.      ! initialize area weight to 1
      MUB(:,:)      = 0.
      DUSDIF(:)     = 0. 
      MU_TOP(:)     = 0. 
      DUGWD(:,:)    = 0. 
CRAD  RDI(:)=960.*960./(TL(:)*TL(:))*EXP(-960./TL(:)) 
C
C**** MOUNTAIN WAVES generate at 1 s.d. above topography ...
      IF (QGWMTN.EQ.1 .AND. RANMTN.LT.0.25d0)  THEN 
        LD(1)=LBREAK
        U0 = (UL(1)*MAUV(1) + UL(2)*MAUV(2)) / (MAUV(1) + MAUV(2))
        V0 = (VL(1)*MAUV(1) + VL(2)*MAUV(2)) / (MAUV(1) + MAUV(2))
        W0=SQRT(U0*U0+V0*V0)
        UR(1)=U0/(W0+ ERR)
        VR(1)=V0/(W0+ ERR)
        BV0 = (BVF(1)*MAUV(1) + BVF(2)*MAUV(2)) / (MAUV(1) + MAUV(2))
        ZVAR=ABS(UR(1))*ZVARX_cell+ABS(VR(1))*ZVARY_cell
        IF(ZVAR.LT.XCDNST(1)*XCDNST(1)) ZVAR=0.
C.... if Froude number (U0/BV0*ZSD) > 1
C.... limit ZSD to be consistent with Froude no. (U0/BV0*ZSD) > 1
        IF (ZVAR.GT.(XFROUD*W0/BV0)**2) ZVAR=(XFROUD*W0/BV0)**2
        P0 = (PL(1)*MAUV(1) + PL(2)*MAUV(2)) / (MAUV(1) + MAUV(2))
        WT(1)=ZWT_cell 
        MU(1)=-CMTN*EK(1)/(H0*ROTK)*P0*BV0*W0*ZVAR
        IF(MU(1)*(UL(LBREAK)*UR(1)+VL(LBREAK)*VR(1)).GE.0.) MU(1)=0.
      END IF
C**** DEFORMATION WAVE (X    d cm-2)
      IF(QGWDEF.EQ.1 .AND. DEFRM_cell.GE.DEFTHRESH !  deform. Threshold
     &               .AND. ZATMO_cell/GRAV.LE.1000.) THEN
        FDEFRM=- CDEF*GRAV/(1000.*DEFTHRESH)*DEFRM_cell 
        DU=UL(LDEF +1)-UL(LDEF )
        DV=VL(LDEF +1)-VL(LDEF )
        DW=SQRT(DU**2        + DV**2       )
        UR(9)=DU       /(DW+ ERR)
        VR(9)=DV       /(DW+ ERR)
        CN(9)=UL(LDEF )*UR(9)+VL(LDEF )*VR(9)
        MU(9)=FDEFRM
        LD(9)=LBREAK
      END IF
C**** WIND SHEAR: USE SHEAR BETWEEN LSHR AND LSHR+1 UNLESS CRIT. LEVEL ABOVE..
      IF (QGWSHR.EQ.1) THEN
        ln=lshr
        lsrc_srch: do ! replaces the goto-style coding commented out below
          l=ln
          cu=.5*(ul(l)+ul(l+1))
          do ln=l+1,ld2-1 ! look for crit. level above this source level
            if((ul(ln)-cu)*(ul(ln+1)-cu).lt.0.) cycle lsrc_srch
          enddo
          exit lsrc_srch  ! no crit. level above this source level
        enddo lsrc_srch
c        LN=LSHR
c 160    L=LN
c        CU=.5*(UL(L)+UL(L+1))
c 170    LN=LN+1
c        IF (LN.GE.LD2) GO TO 172
c        IF ((UL(LN)-CU)*(UL(LN+1)-CU).LT.0.) GO TO 160
c        GO TO 170
c 172    CONTINUE
        LD(2)=L+2
        DU=UL(L+1)-UL(L)
        DV=VL(L+1)-VL(L)
        DW=SQRT(DU*DU+DV*DV)
        UR(2)=DU/(DW+ ERR)
        VR(2)=DV/(DW+ ERR)
        CN(2)=.5*((UL(L+1)+UL(L))*UR(2)+(VL(L+1)+VL(L))*VR(2))
        MU(2)=-CORIOL*PLE(L+1)*DW*DW/(240.*H0*BVF(L+1))
        MU(2)=MU(2)*CSHEAR   
      END IF
C**** MOIST CONVECTIVE MASS FLUX BEGINS TWO LEVELS ABOVE CLOUD...
C**** AMPLITUDE DEPENDS ON |U(SOURCE)-C|.
C**** NOTE:  NM LE 2, NM EQ 4, NM GE 8  ARE ALLOWED FOR MC DRAG
      USRC=0.                   ! needed for diagn
      IF(QGWCNV.EQ.1 .AND. AIRXS.GT.0.) THEN
        IF (PLE(LMC1-1).LE.800.) THEN ! skip MC below 800 mb
        VSRC=0.
        NMX=4
        CLDDEP=PLE(LMC0)-PLE(LMC1)
        FPLUME=AIRXS/CLDDEP
        CLDHT=H0*LOG(PLE(LMC0)/PLE(LMC1))
C
        WTX=FPLUME
#ifdef CUBED_SPHERE /* should be done for all configs */
c empirically compensate for redefinition of MC mass exchange as
c the maximum at any level rather than as the sum over MC layers
c (the sum depends on the number of layers)
        WTX = WTX*10.
        WTX = WTX/(1.+WTX)
#endif
        IF (WTX.GT.1. OR. WTX.LT.0.) THEN
          PRINT *, 'WARNING IN GWDRAG, WTX INCORRECT',WTX,FPLUME,CLDDEP
          IF (WTX.GT.1.) WTX=1.
          IF (WTX.LT.0.) call stop_model(' WTX <0 IN GWDRAG',255)
        END IF
        DO L=LMC0,LMC1-1
          USRC=USRC+UL(L)
          VSRC=VSRC+VL(L)
        END DO
        USRC=USRC/(LMC1-LMC0)
        VSRC=VSRC/(LMC1-LMC0)
        WSRC=SQRT(USRC*USRC+VSRC*VSRC)
        UR(3)=USRC/(WSRC+ ERR)
        VR(3)=VSRC/(WSRC+ ERR)
cc      MU(3)=-EK(3)*CMC*BVF(LMC1-1)*PL(LMC1-1)*CLDHT**2
cc      MU(3)=MU(3)*0.1
cc      MU(4)=MU(3)
        MU3HOLD = -EK(3)*CMC*BVF(LMC1-1)*PL(LMC1-1)*CLDHT**2 
        MU(3) = MU3HOLD * SCVMU  !  default = .188, old = .1806 
        MU(4) = MU3HOLD * SCVMU  !  default = .188, old = .1806 
        CN(3)=WSRC-10.
        CN(4)=WSRC+10.
        UR(4)=UR(3)
        VR(4)=VR(3)
        LD(3)=LPCNV
        LD(4)=LPCNV
        WT(3)=WTX
        WT(4)=WTX
C 
C       Additional Convective Waves 
C 
        UR(3+7) = USRC/(WSRC+ ERR) 
        VR(3+7) =-VSRC/(WSRC+ ERR)  
        MU37HOLD = -EK(3+7)*CMC*BVF(LMC1-1)*PL(LMC1-1)*CLDHT**2 
        MU(3+7) = MU37HOLD * SCVMU  !  default = .188, old = .1806
        MU(4+7) = MU37HOLD * SCVMU  !  default = .188, old = .1806
        CN(3+7) = WSRC-10. 
        CN(4+7) = WSRC+10. 
        UR(4+7) = UR(3+7) 
        VR(4+7) = VR(3+7) 
        LD(3+7)=LPCNV
        LD(4+7)=LPCNV
        WT(3+7)=WTX
        WT(4+7)=WTX
C**** If convection is penetrating (i.e. above PCONPEN) do second set
        IF (PLE(LMC1).LT.PCONPEN .AND. NM.GE.8) THEN
          NMX=8
          DO N=3,NMX
            WT(N)=WTX
            LD(N)=LMC1+1
            WT(N+7)=WTX
            LD(N+7)=LMC1+1
          END DO
          LD(5) = LM+1 
          LD(6) = LM+1 
          LD(5+7) = LM+1 
          LD(6+7) = LM+1 
ccc       CN(5)=WSRC-40.
ccc       CN(6)=WSRC+40.
          CN(7)=WSRC-20.
          CN(8)=WSRC+20.
          CN(7+7)=WSRC-20.
          CN(8+7)=WSRC+20.
          DO N=5,NMX
ccc         MU(N)=MU(3)
            UR(N)=UR(3)
            VR(N)=VR(3)
            UR(N+7)=UR(3+7)
            VR(N+7)=VR(3+7)
          END DO
        ENDIF
C 
        IF (PLE(LMC1).LT.PCONPEN .AND. NM.GE.8) THEN
        MU(5) = 0.0 
        MU(6) = 0.0 
        MU(7) = MU3HOLD * SCVMU  !  default = .188, old = .1806
        MU(8) = MU3HOLD * SCVMU  !  default = .188, old = .1806
        MU(5+7) = 0.0 
        MU(6+7) = 0.0 
        MU(7+7) = MU37HOLD * SCVMU  !  default = .188, old = .1806
        MU(8+7) = MU37HOLD * SCVMU  !  default = .188, old = .1806
        END IF 
C
        WCHECK=UL(LD(3))*UR(3)+VL(LD(3))*VR(3)
        DO N=3,NMX
          IF (WCHECK.GT.CN(N)) CYCLE
          UR(N)=-UR(N)
          VR(N)=-VR(N)
          CN(N)=-CN(N)
        END DO
C 
C       Additional 6 CWs 
C
        WCHECK=UL(LD(3+7))*UR(3+7)+VL(LD(3+7))*VR(3+7)
        DO N=3,NMX
          IF (WCHECK.GT.CN(N+7)) CYCLE
          UR(N+7)=-UR(N+7)
          VR(N+7)=-VR(N+7)
          CN(N+7)=-CN(N+7)
        END DO
C 
        END IF ! skipping shallow convection
      END IF
C****
C**** BREAKING MOMENTUM FLUX AT LAYER EDGES
C****
      DO L=2,LM
        UEDGE(L)=.5*(UL(L-1)+UL(L))
        VEDGE(L)=.5*(VL(L-1)+VL(L))
        TEDGE=.5*(TL(L-1)+TL(L))
        BVEDGE=.5*(BVF(L-1)+BVF(L))
c Rind et al. JAS 45 vol. 3, eqn 10: original coding OK
        BYFACS(L)=-.5*GRAV*PLE(L)/(RGAS*RKBY3*BVEDGE*TEDGE)
c        BYFACS(L)=-.5*GRAV*PLE(L)/(RGAS*SQRT(RKBY3)*BVEDGE*TEDGE)
      END DO
      DO N=1,NM
        DO L=LD(N),LM
          WMCE=UEDGE(L)*UR(N)+VEDGE(L)*VR(N) - CN(N)
          IF (WMCE.GE.0.) THEN
            MUB(L,N)=BYFACS(L)*EK(N)*WMCE**3
            IF (MUB(L,N).GE.-ERR) MUB(L,N)=0.
          END IF
        END DO
      END DO
C**** IF TOP LEVEL IS ABOVE PBREAKTOP FORCE BREAKING GW
      IF (PLE(LM).lt.PBREAKTOP) THEN
C**** DEPOSIT REMAINING MOM. FLUX IN TOP LAYER
        MUB(LM+1,:)=0.
C**** DISTRIBUTE CRIT LEVEL NEAR TOP
        DO N=1,NM
          IF (MUB(LM,N).EQ.0.) MUB(LM,N)=.3d0*MUB(LM-1,N)
        END DO
      END IF
C****
C**** DETERMINE INCIDENT MOMENTUM FLUX
C****
C**** INCIDENT FLUX FOR MTN WAVES
      LN=LD(1)
      IF (MU(1).LT.MUB(LN,1)) MU(1)=MUB(LN,1)
      IF (LN.LE.LM) THEN
        IF (BVF(LN).LT.1.D-5) MU(1)=0.
      END IF
C**** INCIDENT FLUX FOR SHR AND MC WAVES
      DO N=2,NM
        LN=LD(N)
        IF (LN.GT.LM) CYCLE
        IF (MU(N).LT.MUB(LN,N)) MU(N)=MUB(LN,N)
      END DO
C**** LOCATE LDRAG - LOWEST LAYER TO APPLY DRAG
      LDRAG=LM+2    ! better default
      DO N=1,NM
        IF (MU(N).LE.-ERR) LDRAG=MIN(LD(N),LDRAG)
      END DO
C 
      MU_INC(:) = MU(:) 
C****
C**** CALCULATE THE DRAG
C****
      IF (LDRAG.GT.LM) RETURN 
      DO N=1,NM
        IF (LD(N).GT.LM) CYCLE
        LD1=LD(N)
        WMC(LD1-1)=UL(LD1-1)*UR(N)+VL(LD1-1)*VR(N)-CN(N)
        DO L=LD1,LM
          DFM(L)=0.
          DFR(L)=0.
          WMC(L)=UL(L)*UR(N)+VL(L)*VR(N)-CN(N)
          IF(WMC(L).LE..01d0) WMC(L)=MAX(.01d0,.5*(WMC(L-1)+WMC(L)))
          LTOP=L
          IF (MUB(L+1,N).EQ.0.) EXIT
        ENDDO

        DO L=LD1,LTOP
          IF (L.EQ.LM)  MU_TOP(N)=MU(N)*UR(N)*WT(N)
C**** RADIATIVE DAMPING IS OMITTED HERE  
          MUR=MU(N)
CRAD  DTW=DP(L)*BVF(L)/(EK(N)*WMC(L)*WMC(L)*RHO(L)*GRAV+ ERR)
CRAD  WAVEM=1.E3*BVF(L)/(ABS(WMC(L))+ ERR)
CRAD  IF (WAVEM.GT.1.6) WAVEM=1.6d0
CRAD  WM3BY2=WAVEM*SQRT(WAVEM)
CRAD  DTR=86400./(RDI(L)*(AZ(L)+BZ(L)*WAVEM*WAVEM/(CZ(L)+WM3BY2))+ERR)
CRAD  MUR=MU(N)*EXP(-2.*DTW/(DTR+ ERR))
          DFR(L)=MU(N)-MUR
C**** MECHANICAL (TURBULENT) DRAG
          MU(N)=MUR
          IF (MUR.LT.MUB(L+1,N)) MU(N)=MUB(L+1,N) ! saturation drag
          DFM(L)=MUR-MU(N)
        ENDDO
C**** LIMIT THE CONVERGENCE TO XLIMIT*(U-C)
        DO L=LTOP,LD1,-1
          DFT=DFR(L)+DFM(L)
          DFMAX=-XLIMIT*WMC(L)*DP(L)*BYDTIME
          IF (DFT.GE.DFMAX) CYCLE
          EXCESS=DFT-DFMAX
          ALFA=DFM(L)/DFT
          DFR(L)=DFR(L)-EXCESS*(1.-ALFA)
          DFM(L)=DFM(L)-EXCESS*ALFA
          IF (L.GT.LD1) THEN
            DFR(L-1)=DFR(L-1)+EXCESS*(1.-ALFA)
            DFM(L-1)=DFM(L-1)+EXCESS*ALFA
          ENDIF
        END DO
C**** COMBINE MECHANICAL AND RADIATIVE DRAG
        DO L=LD1,LTOP
          DFTL(L)=DFM(L)+DFR(L)
        END DO
C**** CALCULATE DIFFUSION COEFFICIENT AND ADD DRAG TO THE WINDS
C**** (DEFORMATION DIFFUSION IS * XDIFF)
        XDIFF=1. 
        DO L=LD1,LTOP
          DL(L)=DL(L)+ABS(DFM(L)*WMC(L))*XDIFF*WT(N)
          DFT=DFTL(L)
          DWT=DFT*DTIME/DP(L)
          IF (DFT.GT.-ERR) CYCLE
C**** SHEAR AND MOUNTAIN ORIENTATION
          DUTN=DWT*UR(N)*WT(N)
          DUT(L)=DUT(L)+DUTN
          DVT(L)=DVT(L)+DWT*VR(N)
          DUGWD(L,N)=DUTN 
        END DO
      END DO ! end loop over N
C****
C**** MOMENTUM DIFFUSION   (DOUBLED)
C****    (limited to XLIMIT per timestep)
      PDN=PL(LDRAG-1)
      MDN=DP(LDRAG-1)
      YDN=0.
      FLUXUD=0.
      FLUXVD=0.
      DO L=LDRAG,LM
        PUP=PL(L)
        MUP=DP(L)
        YUP=GRAV*GRAV*DTIME*RHO(L)*RHO(L)*DL(L)/(DP(L)*BVF(L)*BVF(L))
        DX=   (YDN+YUP)/(PDN-PUP) ! double diffusion coefficient
        DLIMIT=XLIMIT*MIN(MUP,MDN)
        IF (DX.GT.DLIMIT)  DX=DLIMIT  
        FLUXU=DX/(1.+DX*(MUP+MDN)/(MUP*MDN))*(UL(L)-UL(L-1))
        FLUXV=DX/(1.+DX*(MUP+MDN)/(MUP*MDN))*(VL(L)-VL(L-1))
        DUT(L-1)=DUT(L-1)-(FLUXUD-FLUXU)/MDN
        DVT(L-1)=DVT(L-1)-(FLUXVD-FLUXV)/MDN
        DUSDIF(L-1)=(FLUXU-FLUXUD)/MDN
        FLUXUD=FLUXU
        FLUXVD=FLUXV
        YDN=YUP
        MDN=MUP
        PDN=PUP
      END DO
C**** DIFFUSION IN THE TOP LAYER COMES ONLY FROM BELOW.
      DUT(LM)=DUT(LM)-FLUXUD/MDN
      DVT(LM)=DVT(LM)-FLUXVD/MDN
      DUSDIF(LM) = -FLUXUD/MDN
C  
      RETURN 
      END SUBROUTINE GWDCOL


      SUBROUTINE DFUSEQ(AIRM,DFLX,F,MU,AM,AL,AU,B,LM)
!@sum  DFUSEQ calculate tridiagonal terms
!@auth Bob Suozzo/Jean Lerner
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LM
      REAL*8, INTENT(IN), DIMENSION(LM) :: AIRM
      REAL*8, INTENT(IN), DIMENSION(0:LM+1) :: F
      REAL*8, INTENT(IN), DIMENSION(LM+1) :: DFLX
      REAL*8, INTENT(OUT), DIMENSION(LM) :: AM,AL,AU,B
      REAL*8, INTENT(IN) :: MU
      INTEGER L
      REAL*8 BYAIRM

      DO L=1,LM
        BYAIRM = 1./AIRM(L)
        B(L)=     BYAIRM*(DFLX(L+1)*(F(L+1)-F(L))
     *       - DFLX(L)*(F(L)-F(L-1)))
        AM(L)=1.+ (MU*BYAIRM)*(DFLX(L+1)+DFLX(L))
        if (l.lt.lm) AL(L+1)=-(MU*BYAIRM)*DFLX(L)
        AU(L)=-(MU*BYAIRM)*DFLX(L+1)
      END DO
      RETURN
C****
      END SUBROUTINE DFUSEQ


      SUBROUTINE ALLOC_STRAT_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE STRAT
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: I_0H,I_1H, J_0H,J_1H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               I_STRT_HALO=I_0H, I_STOP_HALO=I_1H)

      ALLOCATE(    DEFRM(I_0H:I_1H,J_0H:J_1H),
     *              EKOFJ(nm,J_0H:J_1H), ! not used by cubed atm
     *             ZVART(I_0H:I_1H,J_0H:J_1H),
     *             ZVARX(I_0H:I_1H,J_0H:J_1H),
     *             ZVARY(I_0H:I_1H,J_0H:J_1H),
     *               ZWT(I_0H:I_1H,J_0H:J_1H),
     *         STAT=IER)

      END SUBROUTINE ALLOC_STRAT_COM


      SUBROUTINE init_GWDRAG
!@sum init_GWDRAG
!@auth Jean Lerner
C**** DO_GWDRAG=true activates the printing of the diagnostics
C**** accumulated in the routines contained herein
      USE Dictionary_mod
      USE CONSTANT, only : twopi,kapa
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : pednl00,pmidl00
      USE DYNAMICS, only : do_gwdrag
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds, AM_I_ROOT
      USE GEOM, only : areag
#ifdef CUBED_SPHERE
      USE STRAT, only : EK_globavg,dfm_type
      use cs2ll_utils, only : init_uv_derivs_type
#else
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE,SOUTH,PACK_DATA
      USE GEOM, only : dxyv,dlat_dg
#endif
      USE STRAT, only : xcdnst, qgwmtn, qgwshr, qgwdef, qgwcnv,lbreak
     *     ,ld2,lshr,ldef,zvarx,zvary,zvart,zwt,nm,ekofj, cmtn,Cshear
     *     ,cdef,cmc,pbreak,pbreaktop,defthresh,pconpen,ang_gwd,LPCNV  
     *     ,SCVMU
      use pario, only : par_open,par_close,read_dist_data
      IMPLICIT NONE
      REAL*8 PLEV,PLEVE,EKS,EK1,EK2,EKX
      REAL*8 :: EK_GLOB(JM)
      INTEGER I,J,L,fid

      INTEGER :: I_0,I_1,J_0,J_1,  I_0H,I_1H,J_0H,J_1H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               I_STRT     =I_0,    I_STOP     =I_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               I_STRT_HALO=I_0H,   I_STOP_HALO=I_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** define flag for optional diagnostics
      DO_GWDRAG = .true.

C**** sync gwdrag parameters from input
      call sync_param( "XCDNST", XCDNST, 2 )
      call sync_param( "CMTN", CMTN)
      call sync_param( "CDEF", CDEF)
      call sync_param( "CMC", CMC)
      call sync_param( "SCVMU", SCVMU)
      call sync_param( "CSHEAR", CSHEAR)
      call sync_param( "PBREAK", PBREAK)
      call sync_param( "PCONPEN", PCONPEN)
      call sync_param( "PBREAKTOP", PBREAKTOP)
      call sync_param( "DEFTHRESH", DEFTHRESH)

C**** sync more gwdrag parameters from input
      call sync_param( "QGWMTN", QGWMTN)
      call sync_param( "QGWSHR", QGWSHR)
      call sync_param( "QGWDEF", QGWDEF)
      call sync_param( "QGWCNV", QGWCNV)
      call sync_param( "ANG_GWD", ANG_GWD)

C**** Calculate levels for deformation etc.
C**** Note: these levels work for the 23 layer model, but may
C**** need testing for other resolutions
      DO L=1,LM
        PLEV  = PMIDL00(L)
        PLEVE = PEDNL00(L)
        IF (PLEV.GE.700) LDEF=L
        IF (PLEVE.GE.PBREAK) LBREAK=L+1
        IF (PLEVE.GE.300.) LSHR=L
c        IF (PLEV.GE.200.) LDEFM=L
        IF (PLEV.GE.0.2d0) LD2=L
        IF (PLEV.GE.200.) LPCNV=L   !  200mb layer for CW 3 & 4  
      END DO
      if (AM_I_ROOT()) then
      WRITE (*,*) ' LEVEL FOR DEFORMATION IS: LDEF,PDEF= ',LDEF,
     *       PMIDL00(LDEF)  ! ,' LDEFM=',LDEFM
      WRITE (*,*) ' LEVELS FOR WIND SHEAR GENERATION: LSHR,LD2= ',LSHR
     *     ,LD2
      WRITE (*,*) ' LBREAK=',LBREAK
      end if
C****
C**** TOPOGRAPHY VARIANCE FOR MOUNTAIN WAVES
C****
      fid = par_open(grid,'ZVAR','read')
      call read_dist_data(grid,fid,'zvart',zvart)
      call read_dist_data(grid,fid,'zvarx',zvarx)
      call read_dist_data(grid,fid,'zvary',zvary)
      call read_dist_data(grid,fid,'zwt',zwt)
      call par_close(grid,fid)

#ifndef CUBED_SPHERE
      CALL HALO_UPDATE(GRID, ZVART, from=SOUTH)
      CALL HALO_UPDATE(GRID, ZVARX, from=SOUTH)
      CALL HALO_UPDATE(GRID, ZVARY, from=SOUTH)
      CALL HALO_UPDATE(GRID,   ZWT, from=SOUTH)
      DO J=J_1S,J_0H,-1
      DO I=1,IM
        ZVART(I,J+1)=ZVART(I,J)
        ZVARX(I,J+1)=ZVARX(I,J)
        ZVARY(I,J+1)=ZVARY(I,J)
        ZWT(I,J+1)=ZWT(I,J)
      END DO
      END DO
#endif

C**** define wave number array EK for GWDRAG
C**** EKX is the mean wave number for wave lengths between a 1x1 degree
C**** box and a model grid box weighted by 1/EK; wave_length=sqrt(area)
#ifdef CUBED_SPHERE
      EK1=TWOPI/SQRT(AREAG/(6*IM*IM)) ! 2pi/grid_box_size
      EK2=TWOPI/SQRT(AREAG/(6*90*90)) ! 2pi/1deg_box_size
      IF(EK1.EQ.EK2) THEN
        EK_globavg = EK1
      ELSE
        EK_globavg=(EK2-EK1)/LOG(EK2/EK1) ! weighted mean
      ENDIF
      call init_uv_derivs_type(grid,dfm_type)
#else
      EKS=0.
      ! full grid for generating the global sum
      DO J=2,JM ! _NOT_ J_0STG,J_1STG
        EK1=TWOPI/SQRT(DXYV(J))                ! 2pi/grid_box_size
        EK2=EK1*SQRT((360./IM)*DLAT_DG)        ! 2pi/1x1deg_box_size
        EKX=EK1
        if(EK2.ne.EK1) EKX=(EK2-EK1)/LOG(EK2/EK1) ! weighted mean
        EKS=EKS+EKX*DXYV(J)
        If ((J >= J_0) .and. (J <= J_1)) THEN
          EKOFJ(1,J)=EKX
          EKOFJ(2,J)=EKX
        END IF
      END DO
      EKS=EKS*IM/AREAG
      EKOFJ(3:NM,J_0:J_1)=EKS
      call PACK_DATA(grid, EKOFJ(1,:), EK_GLOB)
      if (AM_I_ROOT()) then
        WRITE (6,970) (J,EK_GLOB(J),J=2,JM)
        WRITE (6,971) EKS
      end if
  970 FORMAT ('0  J,EK:',9X,1P,7(I4,E12.2)/,9(1X,8(I4,E12.2)/))
  971 FORMAT ('   AVG EK: ',4X,E12.2)
#endif

      END SUBROUTINE init_GWDRAG

#ifndef CUBED_SPHERE /* ends with io_strat */

      Subroutine VDIFF (DT1, UT,VT, U,V,MA, T)
!@sum VDIFF Vertical Diffusion in stratosphere
!@auth Bob Suozzo/Jean Lerner
!**** If GWDRAG is disabled then VDIFF needs halo PEDN,PMID
C****
C**** Vertical diffusion coefficient depends on wind shear.
C**** Uses TRIDIAG for implicit scheme (MU=1) as in diffuse53.
C**** This version only does diffusion for lowest LDIFM layers.
C****
      Use CONSTANT,   Only: RGAS,GRAV,TWOPI,KAPA,SHA,kg2mb
      Use RESOLUTION, Only: IM,JM,LM,LS1=>ls1_nominal
      Use ATM_COM,    Only: PEDN,PMID,PK
      USE DYNAMICS, only : mrch
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE, NORTH, SOUTH
      USE GEOM, only : sini=>siniv,cosi=>cosiv,imaxj,rapvn,rapvs,dxyv
     *     ,kmaxj,idij,idjj,rapj
      USE FLUXES, only : atmsrf
      USE DIAG_COM, only : ajl=>ajl_loc,jl_dudtvdif,byim
      Use STRAT,      Only: DEFRM,ANG_GWD
      USE TRIDIAG_MOD, only :  TRIDIAG
      IMPLICIT NONE

      INTEGER, PARAMETER :: LDIFM=LM
      REAL*8, PARAMETER :: BYRGAS = 1./RGAS
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM+1) ::
     *                                                         VKEDDY
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                             RHO, DUT, DVT, DKE
      REAL*8, DIMENSION(0:LDIFM+1) :: UL,VL,TL,PL,RHOL
      Real*8 :: DP(LM),AM(LM),AL(LM),AU(LM),B(LM),DU(LM),DV(LM),
     *          MAUV(LM)
      REAL*8, DIMENSION(LDIFM+1) :: PLE,RHOE,DPE,DFLX,KMEDGE,KHEDGE
      REAL*8, PARAMETER :: MU=1.
      REAL*8, INTENT(INOUT),
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        U,V,T
      REAL*8, INTENT(IN),
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        UT,VT
      Real*8,Intent(In) :: MA(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      REAL*8, INTENT(IN) :: DT1
      REAL*8 G2DT,PIJ,TPHYS,ANGM,DPT,DUANG
      INTEGER I,J,L,IP1,NDT,N,LMAX
      INTEGER :: J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      REAL*8, DIMENSION(:,:), POINTER :: tsurf,qsurf,usurf,vsurf

      tsurf=>atmsrf%tsavg
      qsurf=>atmsrf%qsavg
      usurf=>atmsrf%usavg
      vsurf=>atmsrf%vsavg

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      G2DT=GRAV*GRAV*DT1
C**** Fill in USURF,VSURF at poles (Shouldn't this be done already?)
      IF (HAVE_SOUTH_POLE) THEN
         DO I=2,IM
            USURF(I,1)=USURF(1,1)*COSI(I)-VSURF(1,1)*SINI(I)
            VSURF(I,1)=VSURF(1,1)*COSI(I)+USURF(1,1)*SINI(I)
         END DO
      ENDIF
      IF (HAVE_NORTH_POLE) THEN
         DO I=2,IM
            USURF(I,JM)=USURF(1,JM)*COSI(I)+VSURF(1,JM)*SINI(I)
            VSURF(I,JM)=VSURF(1,JM)*COSI(I)-USURF(1,JM)*SINI(I)
         END DO
      ENDIF
C**** Calculate RHO(I,J,L)
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        RHO(I,J,L)=   BYRGAS*PMID(L,I,J)/(T(I,J,L)*PK(L,I,J))
      END DO
      END DO
      END DO

C**** Fill in T,RHO at poles (again shouldn't this be done already?)
      IF (HAVE_SOUTH_POLE) THEN
         DO L=1,LM
         DO I=1,IM
            T(I,1,L)=T(1,1,L)
            RHO(I,1,L)=RHO(1,1,L)
         END DO
         END DO
      ENDIF
      IF (HAVE_NORTH_POLE) THEN
         DO L=1,LM
         DO I=1,IM
            T(I,JM,L)=T(1,JM,L)
            RHO(I,JM,L)=RHO(1,JM,L)
         END DO
         END DO
      ENDIF
C**** Get Vertical Diffusion Coefficient for this timestep
      CALL GETVK (U,V,VKEDDY,LDIFM)

      CALL HALO_UPDATE(GRID, USURF , from=SOUTH)
      CALL HALO_UPDATE(GRID, VSURF , from=SOUTH)
      CALL HALO_UPDATE(GRID, TSURF , from=SOUTH)
      CALL HALO_UPDATE(GRID, RHO   , from=SOUTH)
      CALL HALO_UPDATE(GRID, T     , from=SOUTH)

C****
C**** U,V Diffusion
C****
      DUT=0 ; DVT=0.
      DO 300 J=J_0STG,J_1STG
      I=IM
      DO 300 IP1=1,IM
C**** Surface values are used for F(0)
C**** Note area weighting for four point means
      PL(0) = (PEDN(1,I,J-1) + PEDN(1,Ip1,J-1))*RAPVN(J-1) +
     +        (PEDN(1,I,J  ) + PEDN(1,Ip1,J  ))*RAPVS(J)
      UL(0)=(USURF(I  ,J-1) + USURF(IP1,J-1))*RAPVN(J-1) +
     *      (USURF(I  ,J)   + USURF(IP1,J  ))*RAPVS(J)
      VL(0)=(VSURF(I  ,J-1) + VSURF(IP1,J-1))*RAPVN(J-1) +
     *      (VSURF(I  ,J  ) + VSURF(IP1,J  ))*RAPVS(J)
      TPHYS=(TSURF(I  ,J-1) + TSURF(IP1,J-1))*RAPVN(J-1) +
     *      (TSURF(I  ,J)   + TSURF(IP1,J  ))*RAPVS(J)
      TL(0)=TPHYS*PL(0)**KAPA
      RHOL(0)=PL(0)/(RGAS*TPHYS)
      DO L=1,MIN(LDIFM+1,LM)
        UL(L)=U(I,J,L)
        VL(L)=V(I,J,L)
        RHOL(L)=(RHO(I,J-1,L)+RHO(IP1,J-1,L))*RAPVN(J-1)+
     *          (RHO(I,J  ,L)+RHO(IP1,J  ,L))*RAPVS(J)
        MAUV(L) = (  MA(L,I,J-1) +   MA(L,Ip1,J-1))*RAPVN(J-1) +
     +            (  MA(L,I,J  ) +   MA(L,Ip1,J  ))*RAPVS(J)
         PLE(L) = (PEDN(L,I,J-1) + PEDN(L,Ip1,J-1))*RAPVN(J-1) +
     +            (PEDN(L,I,J  ) + PEDN(L,Ip1,J  ))*RAPVS(J)
          PL(L) = (PMID(L,I,J-1) + PMID(L,Ip1,J-1))*RAPVN(J-1) +
     +            (PMID(L,I,J  ) + PMID(L,Ip1,J  ))*RAPVS(J)
      END DO
C**** Edge values at LM+1 don't matter since diffusiv flx=F(L)-F(L-1)=0
      L = LM+1   !!!! A lot relies on LDIFM=LM !!!
        PLE(LM+1) = (PEDN(LM+1,I,J-1) + PEDN(LM+1,Ip1,J-1))*RAPVN(J-1)+
     +              (PEDN(LM+1,I,J  ) + PEDN(LM+1,Ip1,J  ))*RAPVS(J)
         PL(LM+1) = PLE(LM+1)
        UL(L)   =UL(L-1)
        VL(L)   =VL(L-1)
        RHOL(L) =RHOL(L-1)
      DO L=1,LDIFM+1
        KMEDGE(L)=VKEDDY(I,J,L)
        KHEDGE(L)=VKEDDY(I,J,L)
        DPE(L)   =PL(L-1)-PL(L)
      END DO
      DP(:) = PLE(1:LM) - PLE(2:LM+1)
C**** RHOE is obtained from average of T (=P/RHO)
      DO L=1,LDIFM+1
        RHOE(L)=2.*PLE(L)*RHOL(L-1)*RHOL(L)/
     *       (PL(L-1)*RHOL(L)+PL(L)*RHOL(L-1))
      END DO
C**** Calculate diffusive flux and number of timesteps NDT
      DO L=1,LDIFM+1
        DFLX(L)=G2DT*RHOE(L)**2*KMEDGE(L)/DPE(L)
      END DO
      NDT=1

      DO N=1,NDT
C**** dq/dt by diffusion as tridiagonal matrix
         Call DFUSEQ (DP,DFLX,UL(0),MU,AM,AL,AU,B,LDIFM)
        CALL TRIDIAG(Al,Am,AU,B,DU,LDIFM)
         Call DFUSEQ (DP,DFLX,VL(0),MU,AM,AL,AU,B,LDIFM)
        CALL TRIDIAG(Al,Am,AU,B,DV,LDIFM)
C**** Update model winds
        IF (MRCH.GT.0) THEN
          DO L=1,LM
            AJL(J,L,JL_DUDTVDIF) = AJL(J,L,JL_DUDTVDIF) + DU(L)*BYIM
            DUT(I,J,L) = DUT(I,J,L) + DU(L)*MAUV(L)*DXYV(J)
            DVT(I,J,L) = DVT(I,J,L) + DV(L)*MAUV(L)*DXYV(J)
            DKE(I,J,L) = DU(L)*(U(I,J,L)+0.5*DU(L))+
     *                   DV(L)*(V(I,J,L)+0.5*DV(L))
          END DO
        END IF
C**** Save AM change and update U,V
        ANGM = 0.
        DO L=1,LM
          ANGM = ANGM - DU(L)*MAUV(L)
          U(I,J,L) = U(I,J,L) + DU(L)
          V(I,J,L) = V(I,J,L) + DV(L)
        END DO

        if (ang_gwd.gt.0) then  ! add in ang mom
          lmax=ls1-1            ! below ptop
          if (ang_gwd.gt.1) lmax=lm ! over whole column
          DUANG = ANGM * kg2mb / (PLE(1) - PLE(LMAX+1))
          IF (MRCH.GT.0) THEN
            DO L=1,LMAX
              DKE(I,J,L) = DKE(I,J,L) + DUANG*(U(I,J,L)+0.5*DUANG)
              DUT(I,J,L) = DUT(I,J,L) + DUANG*MAUV(L)*DXYV(J)
              AJL(J,L,JL_DUDTVDIF) = AJL(J,L,JL_DUDTVDIF) + DUANG*BYIM
            END DO
          END IF
          DO L=1,LMAX
            U(I,J,L) = U(I,J,L) + DUANG
          END DO
        END IF

      END DO
  300 I=IP1
C**** conservation diagnostic
      IF (MRCH.gt.0) THEN
        CALL DIAGCD (GRID,6,UT,VT,DUT,DVT,DT1)

        call regrid_btoa_3d(dke)
c the optional argument diagIndex will be restored once this
c routine goes back into a module
        call addEnergyAsLocalHeat(DKE, T, PK)!, diagIndex=JL_dTdtsdrg)
      END IF

      RETURN
C****
      END SUBROUTINE VDIFF

      SUBROUTINE GETVK (U,V,VKEDDY,LDIFM)
!@sum GETVK calculate vertical diff. coefficient (use high wind shear)
!@auth Bob Suozzo/Jean Lerner
      USE RESOLUTION, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, HALO_UPDATE,
     *                          NORTH, SOUTH
      USE FLUXES, only : atmsrf
      IMPLICIT NONE
!@var Vert. Diffusion coefficent
      REAL*8, INTENT(OUT),
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM+1) ::
     *                                                         VKEDDY
      REAL*8, INTENT(IN),
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                          U,V
      INTEGER, INTENT(IN) :: LDIFM
      REAL*8, PARAMETER :: XEDDY = 10., DV2MAX = 25.**2
      INTEGER I,J,L,IP1
      REAL*8 US1,VS1,DELV2

      INTEGER :: J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      REAL*8, DIMENSION(:,:), POINTER :: usurf,vsurf

      usurf=>atmsrf%usavg
      vsurf=>atmsrf%vsavg

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      CALL HALO_UPDATE(GRID, USURF , from=SOUTH)
      CALL HALO_UPDATE(GRID, VSURF , from=SOUTH)

C**** Calculate surface winds on velocity grid (rewrite!)
C**** Is this calculated in PBL?
      DO J=J_0STG,J_1STG
        I=IM
        DO IP1=1,IM
          US1=.25*(USURF(I,J-1)+USURF(IP1,J-1)+USURF(I,J)+USURF(IP1,J))
          VS1=.25*(VSURF(I,J-1)+VSURF(IP1,J-1)+VSURF(I,J)+VSURF(IP1,J))
          DELV2   =(U(I,J,1)-US1)**2 + (V(I,J,1)-VS1)**2
          VKEDDY(I,J,1)=0.
          IF (DELV2.GT.DV2MAX) VKEDDY(I,J,1)=XEDDY
          I=IP1
        END DO
      END DO

      DO L=2,MIN(LDIFM+1,LM)
      DO J=J_0STG,J_1STG
      DO I=1,IM
        DELV2   =(U(I,J,L)-U(I,J,L-1))**2 + (V(I,J,L)-V(I,J,L-1))**2
        VKEDDY(I,J,L)=0.
        IF (DELV2.GT.DV2MAX) VKEDDY(I,J,L)=XEDDY
      END DO
      END DO
      END DO
      IF (LDIFM.EQ.LM) THEN
        DO J=J_0STG,J_1STG
          DO I=1,IM
            VKEDDY(I,J,LM+1)=0.
          END DO
        END DO
      ENDIF
      DO J=J_0STG,J_1STG
        DO I=1,IM
          VKEDDY(I,J,1)=0.
          VKEDDY(I,J,2)=0.
        END DO
      END DO
      RETURN
C****
      END SUBROUTINE GETVK


      Subroutine GWDRAG (DT1, UT,VT, U,V,MA, T,SZ, CALC_DEFORM)
!@sum  GWDRAG puts a momentum drag in the stratosphere
!@auth Bob Suozzo/Jean Lerner
!**** If GWDRAG is disabled then VDIFF needs halo PEDN,PMID
C****
C**** GWDRAG is called from DYNAM with arguments:
!**** Input: DT1 = timestep (s)
!****      UT,VT = center velocity of time step (m/s)
!****         MA = mass per unit area at end of timestep (kg/m^2)
!**** Output: U,V,T = wind and Potential Temp. at beginning of timestep
C****
      Use CONSTANT,   Only: GRAV,SHA,KAPA,RGAS,kg2mb
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only : dtsrc
      Use ATM_COM,    Only: ZATMO, PEDN,PMID,PDSIG,PK
      USE DYNAMICS, only : mrch
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, ONLY : getDomainBounds, HALO_UPDATE,
     *                          NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE_COLUMN, am_i_root
      USE CLOUDS_COM,       ONLY : AIRX,LMC
      Use STRAT,      Only: NM, ZVARX,ZVARY,ZWT,DEFRM, QGWCNV,EKOFJ,
     *                      pbreaktop,defthresh,pconpen,ang_gwd,LPCNV,  
     *                      TL,THL,RHO,BVF,DL,DUT,DVT,UL,VL,
     *                      DTIME,BYDTIME,CORIOL,AIRXS,
     &     UR, WT, CN, USRC, DUSDIF, MU_TOP, DUGWD, MU_INC, EK,
     &     LDRAG, LMC0, LMC1,
     &     RANMTN_CELL=>RANMTN,
     &     DEFRM_CELL,ZVARX_CELL,ZVARY_CELL,ZWT_CELL,ZATMO_CELL
      USE GEOM, only : dxyv,bydxyv,fcor,imaxj,rapvn,rapvs
      USE DIAG_COM, only : byim, aij=>aij_loc, ajl=>ajl_loc
     *     ,jl_dudtsdif,jl_sdifcoef,jl_gwfirst
     *     ,ij_gw1,ij_gw2,ij_gw3,ij_gw4,ij_gw5
     *     ,ij_gw6,ij_gw7,ij_gw8,ij_gw9
      USE RANDOM
      IMPLICIT NONE
C 
      REAL*8, PARAMETER :: ERR=1d-20 
C 
      REAL*8, INTENT(INOUT),
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                      T,U,V,SZ
      REAL*8, INTENT(IN),
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        UT,VT
      Real*8,Intent(In) :: MA(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      REAL*8, INTENT(IN) :: DT1
      LOGICAL, INTENT(IN) :: CALC_DEFORM

!**** Local variables
      Real*8 :: PLE(LM+1),PL(LM),MAUV(LM)
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *     DUT3,DVT3,DKE,TLS,THLS,BVS
      REAL*8, DIMENSION(IM,LM) :: UIL,VIL,TLIL,THIL,BVIL
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: DUJL
      INTEGER I,IP1,J,L,N
      Real*8 :: BVFSQ,ANGM,DPT,DUANG,XY,AIRX4,DTHR
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: RANMTN  
      integer lmax_angm(nm),lp10,lp2040,lpshr
C
      INTEGER :: J_1, J_0
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****
      DTHR = DT1/DTsrc  ! = 1/(#calls/phys.call)
      DTIME = DT1
      BYDTIME = 1./DTIME
C
C     Set Random numbers for 1/4 Mountain Drag
C
      CALL BURN_RANDOM(IM*(J_0-1))
      DO J=J_0,J_1
      DO I=1,IM
        RANMTN(I,J) = RANDU(XY)
      END DO
      END DO
      CALL BURN_RANDOM(IM*(JM-J_1))
C****
C**** FILL IN QUANTITIES AT POLES
C****
      IF (HAVE_SOUTH_POLE) THEN
         DO I=2,IM
            AIRX(I,1)=AIRX(1,1)
            LMC(1,I,1)=LMC(1,1,1)
            LMC(2,I,1)=LMC(2,1,1)
         END DO
         DO L=1,LM
         DO I=2,IM
            T(I,1,L)=T(1,1,L)
            SZ(I,1,L)=SZ(1,1,L)
         END DO
         END DO
      ENDIF
      IF (HAVE_NORTH_POLE) THEN
         DO I=2,IM
            AIRX(I,JM)=AIRX(1,JM)
            LMC(1,I,JM)=LMC(1,1,JM)
            LMC(2,I,JM)=LMC(2,1,JM)
         END DO
         DO L=1,LM
         DO I=2,IM
            T(I,JM,L)=T(1,JM,L)
            SZ(I,JM,L)=SZ(1,JM,L)
         END DO
         END DO
      ENDIF
C
      CALL HALO_UPDATE(GRID, SZ    , from=SOUTH)
      CALL HALO_UPDATE(GRID, T     , from=SOUTH)
      CALL HALO_UPDATE(GRID, AIRX  , from=SOUTH)
      CALL HALO_UPDATE_COLUMN(GRID, LMC   , from=SOUTH)
C****
C**** DEFORMATION
C****
      If (CALC_DEFORM)  Call DEFORM (MA,U,V)
      DO L=1,LM
        DUT3(:,:,L)=0. ; DVT3(:,:,L)=0. ; DUJL(:,L)=0.; DKE(:,:,L) = 0.
      END DO
      DO L=1,LM
      DO J=J_0STG,J_1STG
      I=IM
      DO IP1=1,IM
         TLS(I,J,L) =
     *    .25*(PK(L,I,J-1)*T(I,J-1,L)+PK(L,IP1,J-1)*T(IP1,J-1,L)+
     *     PK(L,I,J)*T(I,J,L)+PK(L,IP1,J)*T(IP1,J,L))
         THLS(I,J,L) =
     *    .25*(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))
         BVS(I,J,L) =
     *    .5*(SZ(I,J-1,L)+SZ(IP1,J-1,L)+SZ(I,J,L)+SZ(IP1,J,L))
         I=IP1
       END DO
       END DO
       END DO
C****
C**** BEGINNING OF OUTER LOOP OVER I,J
C****
      DO J=J_0STG,J_1STG
C**** parallel reductions
         UIL(:,:)=U(:,J,:)
         VIL(:,:)=V(:,J,:)
         TLIL(:,:)=TLS(:,J,:)
         THIL(:,:)=THLS(:,J,:)
         BVIL(:,:)=BVS(:,J,:)

      CN(:)=0.
      CORIOL=(ABS(FCOR(J-1))+ABS(FCOR(J)))*BYDXYV(J)
      I=IM
      DO IP1=1,IM
C****
C**** CALCULATE VERTICAL ARRAYS
C****
      Do L=1,LM+1
         PLE(L) = (PEDN(L,I,J-1) + PEDN(L,Ip1,J-1))*RAPVN(J-1) +
     +            (PEDN(L,I,J  ) + PEDN(L,Ip1,J  ))*RAPVS(J)  ;  EndDo
      DO L=1,LM
          PL(L) = (PMID(L,I,J-1) + PMID(L,Ip1,J-1))*RAPVN(J-1) +
     +            (PMID(L,I,J  ) + PMID(L,Ip1,J  ))*RAPVS(J)
        MAUV(L) = (  MA(L,I,J-1) +   MA(L,Ip1,J-1))*RAPVN(J-1) +
     +            (  MA(L,I,J  ) +   MA(L,Ip1,J  ))*RAPVS(J)
      TL(L)=TLIL(I,L)
      THL(L)=THIL(I,L)
      RHO(L)=PL(L)/(RGAS*TL(L))
      BVFSQ = BVIL(I,L)*GRAV*GRAV*RHO(L) / (MAUV(L)*kg2mb * THL(L))
      IF (PL(L).GE..4d0) THEN
        BVF(L)=SQRT(MAX(BVFSQ,1.d-10))
      ELSE
        BVF(L)=SQRT(MAX(BVFSQ,1.d-4))
      END IF
C 
      DL(L)=0.
      DUT(L)=0.
      DVT(L)=0.
      UL(L)=UIL(I,L)   ! U(I,J,L)
      VL(L)=VIL(I,L)   ! V(I,J,L)
      END DO
C**** Levels for angular momentum restoration
      do l=lm,1,-1
        if (pl(l).lt.500.) lp10 = l
        if (pl(l).lt.400.) lp2040 = l
        if (pl(l).lt.200.) lpshr = l
      end do
      lpshr  = lpshr -1
      lp2040 = lp2040-1
      lp10   = lp10  -1
      lmax_angm(1) = lpshr  !Mountain
      lmax_angm(9) = lpshr  !Deformation
      lmax_angm(2) = lpshr  !Shear
      lmax_angm(3) = lp10
      lmax_angm(4) = lp10
      lmax_angm(5) = lp2040
      lmax_angm(6) = lp2040
      lmax_angm(7) = lp2040
      lmax_angm(8) = lp2040
C 
C     Additional CWs 
C 
      lmax_angm(3+7) = lp10
      lmax_angm(4+7) = lp10
      lmax_angm(5+7) = lp2040
      lmax_angm(6+7) = lp2040
      lmax_angm(7+7) = lp2040
      lmax_angm(8+7) = lp2040
C 
      RANMTN_CELL = RANMTN(I,J) 
      ZVARX_CELL  = ZVARX(I,J) 
      ZVARY_CELL  = ZVARY(I,J) 
      ZWT_CELL    = ZWT(I,J) 
      ZATMO_CELL  = ZATMO(I,J) 
      DEFRM_CELL  = DEFRM(I,J) 
C 
      AIRXS = 0.0 
      LMC0  = 0 
      LMC1  = 0 
      IF (QGWCNV.EQ.1) THEN
        AIRX4=AIRX(I,J-1)+AIRX(IP1,J-1)+AIRX(I,J)+AIRX(IP1,J)
        AIRXS= ((AIRX(I,J-1)+AIRX(IP1,J-1))*RAPVN(J-1)
     *       +  (AIRX(I,J  )+AIRX(IP1,J  ))*RAPVS(J))
        AIRXS=AIRXS/DXYV(J)
        IF (AIRX4.GT.0.)  THEN
        LMC0=.5+(LMC(1,I,J-1)*AIRX(I,J-1)+LMC(1,IP1,J-1)*AIRX(IP1,J-1)+
     *       LMC(1,I,J)*AIRX(I,J)+LMC(1,IP1,J)*AIRX(IP1,J))/(AIRX4+ ERR)
        LMC1=.5+(LMC(2,I,J-1)*AIRX(I,J-1)+LMC(2,IP1,J-1)*AIRX(IP1,J-1)+
     *       LMC(2,I,J)*AIRX(I,J)+LMC(2,IP1,J)*AIRX(IP1,J))/(AIRX4+ ERR)
        END IF 
      END IF 
C 
      DO N=1,NM 
         EK(N) = EKOFJ(N,J) 
      END DO 

C 
C     Call the GWD Column 
C 
C**************************************
      Call GWDCOL (PLE,PL,MAUV) 
C**************************************
C 
      IF(LDRAG.LE.LM) THEN
C****
C****     UPDATE DIAGNOSTICS
C****
      IF (MRCH.EQ.2) THEN
         AIJ(I,J,IJ_GW1)=AIJ(I,J,IJ_GW1)+MU_INC(9)*UR(9)*DTHR
         AIJ(I,J,IJ_GW2)=AIJ(I,J,IJ_GW2)+MU_INC(1)*UR(1)*DTHR  *WT(1)
         AIJ(I,J,IJ_GW3)=AIJ(I,J,IJ_GW3)+MU_INC(2)*UR(2)*DTHR
         AIJ(I,J,IJ_GW4)=AIJ(I,J,IJ_GW4)
     &        +(MU_INC(3)*UR(3)+MU_INC(3+7)*UR(3+7))*(DTHR*WT(3))
         AIJ(I,J,IJ_GW5)=AIJ(I,J,IJ_GW5)
     &        +(MU_INC(7)*UR(7)+MU_INC(7+7)*UR(7+7))*(DTHR*WT(7))
         AIJ(I,J,IJ_GW6)=AIJ(I,J,IJ_GW6)+MU_INC(5)*UR(5)*DTHR  *WT(5)
         AIJ(I,J,IJ_GW7)=AIJ(I,J,IJ_GW7)+CN(2)*UR(2)*DTHR
         AIJ(I,J,IJ_GW8)=AIJ(I,J,IJ_GW8)+USRC*DTHR
C 
         DO N=1,NM
            AIJ(I,J,IJ_GW9)=AIJ(I,J,IJ_GW9)+MU_TOP(N)*DTHR
         END DO 
C
         DO N=1,9 ! 9 is without additional convective waves
           DO L=1,LM
             AJL(J,L,N+JL_gwFirst-1)=AJL(J,L,N+JL_gwFirst-1)
     &          +DUGWD(L,N)*BYIM
           END DO
           if(n.ge.3 .and. n.le.8) then ! additional convective waves
             DO L=1,LM
               AJL(J,L,N+JL_gwFirst-1)=AJL(J,L,N+JL_gwFirst-1)
     &          +DUGWD(L,N+7)*BYIM
             END DO
           endif
         END DO 

C
         DO L=LDRAG-1,LM
            AJL(J,L,JL_DUDTSDIF)=AJL(J,L,JL_DUDTSDIF)+DUSDIF(L)*BYIM
         END DO
         DO L=LDRAG,LM
          DUJL(J,L)=DL(L)/(BVF(L)*BVF(L))*DTHR
         END DO
C****
C**** Save KE change and diffusion coefficient on A-grid
C****
         DO L=LDRAG-1,LM
          DKE(I,J,L)=DUT(L)*(0.5*DUT(L)+UIL(I,L)) +
     *               DVT(L)*(0.5*DVT(L)+VIL(I,L))
          DUT3(I,J,L) = DUT(L)*MAUV(L)*DXYV(J)
          DVT3(I,J,L) = DVT(L)*MAUV(L)*DXYV(J)
         END DO
      END IF
C 
C**** Save AM change
      if (ang_gwd.gt.0) then    ! add in ang mom
        DO N=1,NM
           ANGM = - Sum (DUGWD(LDRAG-1:LM,N)*MAUV(LDRAG-1:LM))
           DPT  = Sum (MAUV(1:LMAX_ANGM(N)))
          DUANG = ANGM/DPT
          IF (MRCH.eq.2) THEN
            DO L=1,LMAX_ANGM(N)
              DKE(I,J,L) = DKE(I,J,L)+DUANG*(0.5*DUANG+UIL(I,L)+DUT(L))
              DUT3(I,J,L)=DUT3(I,J,L)+DUANG*MAUV(L)*DXYV(J)
            END DO
          END IF
          DO L=1,LMAX_ANGM(N)
             DUT(L) = DUT(L) + DUANG
          END DO
        END DO
      end if
C
C****
C**** UPDATE THE U AND V WINDS
C****
      DO L=1,LM
        UIL(I,L)=UIL(I,L)+DUT(L)
        VIL(I,L)=VIL(I,L)+DVT(L)
      END DO

      ENDIF ! GW drag occurred

C**** END OF LOOP OVER I
      I=IP1
      END DO
C
         U(:,J,:)=UIL(:,:)
         V(:,J,:)=VIL(:,:)
C**** END OF LOOP OVER J
      END DO
C****

      IF (MRCH.EQ.2) THEN
C**** conservation diagnostic
        CALL DIAGCD (GRID,6,UT,VT,DUT3,DVT3,DTIME)

C**** PUT THE KINETIC ENERGY BACK IN AS HEAT
        call regrid_btoa_3d(dke)
c the optional argument diagIndex will be restored once this
c routine goes back into a module
        call addEnergyAsLocalHeat(DKE,T,PK)!, diagIndex=JL_dTdtsdrg)

        DO L=1,LM
          DO J=J_0,J_1
            AJL(J,L,JL_SDIFCOEF)=AJL(J,L,JL_SDIFCOEF)+ DUJL(J,L)*BYIM
          END DO
        END DO

      END IF
C****
      RETURN
      END SUBROUTINE GWDRAG


      Subroutine DEFORM (MA,U,V)
!@sum  DEFORM calculate defomation terms
!@auth Bob Suozzo/Jean Lerner
C****
C**** Deformation terms  DEFRM1=du/dx-dv/dy   DEFRM2=du/dy+dv/dx
C****
      USE RESOLUTION, only : im,jm,lm
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, HALO_UPDATE,
     *                             NORTH, SOUTH,
     *     haveLatitude
      Use DYNAMICS,   Only: MU,MV
      USE GEOM, only : bydxyv,dxyv,dxv,dyp
      USE STRAT, only : ldef,defrm
      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *     DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: U,V
      Real*8,Intent(In):: MA(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)

      REAL*8, DIMENSION(IM) :: UDXS,DUMS1,DUMS2,DUMN1,DUMN2
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *     DEFRM1,DEFRM2  ! ,DEF1A,DEF2A
!      CHARACTER*80 TITLE
      INTEGER I,J,L,IP1,IM1
      REAL*8 UDXN

      INTEGER :: J_1, J_0, J_0H, J_1H
      INTEGER :: J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****
C**** Deformation terms  DEFRM1=du/dx-dv/dy   DEFRM2=du/dy+dv/dx
C**** For spherical coordinates, we are treating DEFRM1 like DIV
C**** and DEFRM2 like CURL (i.e., like FLUX and CIRCULATION),
C**** except the "V" signs are switched
C****

      CALL HALO_UPDATE(GRID, U  , from=NORTH)
      CALL HALO_UPDATE(GRID, V  , from=NORTH)

      L=LDEF
C**** U-terms
      DO 91 J=J_0,J_1
      IM1=IM
      DO 91 I=1,IM
      DEFRM1(I,J) = MU(I,J,L) - MU(Im1,J,L)
   91 IM1=I
      DO 92 I=1,IM
   92 UDXS(I)=0.
      ! Need to seed the recursion if not southern most PE.
      If (.not. HAVE_SOUTH_POLE) THEN
        J=J_0H
        IM1=IM
        DO I=1,IM
          UDXS(I)=.5*(U(IM1,J+1,L)+U(I,J+1,L))*DXV(J+1)
          IM1=I
        END DO
      END IF

      IM1=IM
      DO 93 J=J_0,J_1S
      DO 93 I=1,IM
      UDXN=.5*(U(IM1,J+1,L)+U(I,J+1,L))*DXV(J+1)
      DEFRM2(I,J ) = UDXN-UDXS(I)
      UDXS(I)=UDXN
   93 IM1=I
      IF (HAVE_NORTH_POLE) THEN
         DO 94 I=1,IM
   94    DEFRM2(I,JM) =     -UDXS(I)
      ENDIF
C**** V-terms
      IM1=IM
      DO 98 I=1,IM
      DO 96 J=J_0,J_1S
   96 DEFRM1(I,J) = DEFRM1(I,J) - (MV(I,J+1,L) - MV(I,J,L))
      If (HAVE_SOUTH_POLE)  DEFRM1(I,1 ) = DEFRM1(I,1 ) - MV(I,2 ,L)
      If (HAVE_NORTH_POLE)  DEFRM1(I,JM) = DEFRM1(I,JM) + MV(I,JM,L)
      DO 97 J=J_0S,J_1S
   97 DEFRM2(I,J ) = DEFRM2(I,J ) +
     *  .5*((V(I,J,L)+V(I,J+1,L))-(V(IM1,J,L)+V(IM1,J+1,L)))*DYP(J)
      IF (HAVE_SOUTH_POLE) DEFRM2(I,1 ) = DEFRM2(I,1 ) +
     *                                   (V(I,2 ,L)-V(IM1,2 ,L))*DYP(1)
      IF (HAVE_NORTH_POLE) DEFRM2(I,JM) = DEFRM2(I,JM) +
     *                                   (V(I,JM,L)-V(IM1,JM,L))*DYP(JM)
   98 IM1=I
C**** Convert to UV-grid
      CALL HALO_UPDATE(grid,DEFRM1,FROM=SOUTH)
      CALL HALO_UPDATE(grid,DEFRM2,FROM=SOUTH)
      I=IM
      DO 99 IP1=1,IM
        DUMS1(I)=DEFRM1(I,J_0STG-1)+DEFRM1(IP1,J_0STG-1)
        DUMS2(I)=DEFRM2(I,J_0STG-1)+DEFRM2(IP1,J_0STG-1)
 99     I=IP1
      DO 110 J=J_0STG,J_1STG
      I=IM
      DO 100 IP1=1,IM
      DUMN1(I)=DEFRM1(I,J)+DEFRM1(IP1,J)
      DUMN2(I)=DEFRM2(I,J)+DEFRM2(IP1,J)
  100 I=IP1
      DO 105 I=1,IM
      DEFRM1(I,J)=DUMS1(I)+DUMN1(I)
  105 DEFRM2(I,J)=DUMS2(I)+DUMN2(I)
      DO 110 I=1,IM
      DUMS1(I)=DUMN1(I)
      DUMS2(I)=DUMN2(I)
  110 CONTINUE

      DO 120 J=J_0STG,J_1STG
      I=IM
      DO 120 IP1=1,IM
      DEFRM1(I,J) = DEFRM1(I,J) / (DXYV(J)*
     *              (MA(L,I,J)+MA(L,Ip1,J)+MA(L,I,J-1)+MA(L,Ip1,J-1)))
      DEFRM2(I,J)=.25*DEFRM2(I,J)*BYDXYV(J)
      DEFRM(I,J)=SQRT(DEFRM1(I,J)**2+DEFRM2(I,J)**2)
  120 I=IP1
C**** Set deformation to zero near the poles
      if (haveLatitude(grid, J=2)) DEFRM(1:IM,2) = 0.
      if (haveLatitude(grid, J=3)) DEFRM(1:IM,3) = 0.

      if (haveLatitude(grid, J=JM-1)) DEFRM(1:IM,JM-1) = 0.
      if (haveLatitude(grid, J=JM)) DEFRM(1:IM,JM) = 0.


C****
C**** Print Deformation terms
C****
CWF   IF (MOD(TAU,8.D0).EQ.0.)  RETURN
CW    D1MAX=0.
CW    D2MAX=0.
CWF   DFMAX=0.
CWF   DO 220 J=4,JM-2
CWF   DO 220 I=1,IM
CW    DEF1A(I,J)=ABS(DEFRM1(I,J))
CW    DEF2A(I,J)=ABS(DEFRM2(I,J))
CW    IF(D1MAX.LT.DEF1A(I,J)) D1MAX=DEF1A(I,J)
CW    IF(D2MAX.LT.DEF2A(I,J)) D2MAX=DEF2A(I,J)
CWF   IF(DFMAX.LT.DEFRM(I,J)) DFMAX=DEFRM(I,J)
  220 CONTINUE
CWF   IHR=NINT(TAU)
CW    SCALE=100./10.**(INT(LOG10(D1MAX)))
CW    TITLE='ABS DEFORMATION TERM 1 = DU/DX - DV/DY (XXXXXXXX S-1)'
CW    WRITE (TITLE(41:48),'(1PE8.0)') SCALE
CCW   CALL MAP0(IM,JM,IHR,TITLE,DEF1A,SCALE,0., 0)
CW    SCALE=100./10.**(INT(LOG10(D2MAX)))
CW    TITLE='ABS DEFORMATION TERM 2 = DU/DY - DV/DX (XXXXXXXX S-1)'
CW    WRITE (TITLE(41:48),'(1PE8.0)') SCALE
CCW   CALL MAP0(IM,JM,IHR,TITLE,DEF2A,SCALE,0., 0)
CWF   SCALE=100./10.**(INT(LOG10(DFMAX)))
CWF   TITLE='RMS OF DEFORMATION TERMS 1 AND 2     (XXXXXXXX S-1)'
CWF   WRITE (TITLE(38:45),'(1PE8.0)') SCALE
CWF   CALL MAP0(IM,JM,IHR,TITLE,DEFRM,SCALE,0., 0)
      RETURN
C****
      END SUBROUTINE DEFORM

#endif /* not cubed sphere */
