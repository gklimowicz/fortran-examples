#include "rundeck_opts.h" 
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

      MODULE LAKES
!@sum  LAKES subroutines for Lakes and Rivers
!@auth Gavin Schmidt/Gary Russell
!@ver  2010/08/04 (based on LB265)
!**** Changes from Model III: MO -> MWL (kg), G0M -> GML (J), GZM -> TLAKE (deg C)
      Use CONSTANT, Only: grav,bygrav,shw,rhow,lhm,shi,teeny,undef
#ifdef TRACERS_WATER
      Use OldTracer_mod, Only: trname
      Use TRACER_COM, Only: NTM
#endif
      IMPLICIT NONE
      SAVE

      Integer,Allocatable,Dimension(:,:) :: &
         KDIREC, &  !@var KDIREC directions for river flow, 0 = no flow, 1-8 = anti-clockwise from top RH corner
         KD911,  &  !@var KD911 emergency directions for river flow, 1-8 = anti-clockwise from top RH corner
         IFLOW,JFLOW, &  !@var IFLOW,JFLOW grid box indexes for downstream direction
         IFL911,JFL911   !@var IFL911,JFL911 grid box indexes for emergency downstream direction
      Real*8,Allocatable,Dimension(:,:) :: &
         dZSG , &  !@var dZSG (m) = continental range of solid ground topography in grid cell
         RATE, &  !@var RATE rate of river flow downslope (fraction)
         DHORZ    !@var DHORZ horizontal distance to downstream box (m)
      Real*8,Parameter :: &
         MINMLD    = 1, &  !@param MINMLD minimum mixed layer depth in lake (m)
         HLAKE_MIN = 1, &  !@param HLAKE_MIN minimum sill depth for lake (m)
         TMAXRHO   = 4, &  !@param TMAXRHO temperature of maximum density (pure water) (C)
         KVLAKE = 1d-5, &  !@param KVLAKE lake diffusion constant at mixed layer depth (m^2/s)
         TFL       = 0, &  !@param TFL freezing temperature for lakes (=0 C)
         FLEADLK   = 0, &  !@param FLEADLK lead fraction for lakes
         byZETA=1/.35d0,&  !@param BYZETA reciprocal of solar rad. extinction depth for lake (1/m)
         AC1LMIN = .1, AC2LMIN = .1  !@param AC1LMIN, AC2LMIN minimum ice thickness for lake ice (kg/m^2) not used
      Real*8 :: &
         RIVER_FAC = 1  !@dbparam river_fac Factor to multiply runoff by to balance sea level
      Integer :: &
         INIT_FLAKE = 1, &   !@dbparam init_flake used to make sure FLAKE is properly initialised when using older restart files
                             !@+   = 0 for no change to lakes
                             !@+   = 1 for a complete reset to initial values
                             !@+   = 2 removal of any excess water that may have accumulated
         VARIABLE_LK = 0, &  !@dbparam variable_lk 1 if lakes are to be variable (temporary variable for development purposes)
         LAKE_ICE_MAX =5, &  !@dbparam  Lake ice exceeding lake_ice_max (m) of water equivalent is dumped into ice berg arrays
         LAKE_RISE_MAX = 100 !@dbparam lake_rise_max amount of lake rise (m) over sill level before spillover into next box
                             !@+   (only if lake covers >95% of box)


      CONTAINS


      SUBROUTINE LKSOURC (I0,J0,ROICE,MLAKE,ELAKE,RUN0,FODT,FIDT,SROX,FSR2, &
#ifdef TRACERS_WATER
                          TRLAKEL,TRUN0,TREVAP,TRO,TRI, &
#endif
                          EVAPO,ENRGFO,ACEFO,ACEFI,ENRGFI)
!@sum  LKSOURC applies fluxes to lake in ice-covered and ice-free areas
!@auth Gary Russell/Gavin Schmidt
      Use MODEL_COM, Only: qcheck
      Use DOMAIN_DECOMP_ATM, Only: WRITE_PARALLEL
      IMPLICIT NONE
      Real*8,Intent(InOut) :: &
         MLAKE(2), &  !@var MLAKE = mass in lake layers (kg)
         ELAKE(2)     !@var ELAKE = heat content in lake layers (J/m^2)
      INTEGER,INTENT(IN)  :: I0,J0
      REAL*8, INTENT(IN)  :: ROICE, EVAPO, RUN0, FODT, FIDT, SROX(2)
      REAL*8, INTENT(OUT) :: ENRGFO, ACEFO, ENRGFI, ACEFI
#ifdef TRACERS_WATER
      REAL*8, INTENT(INOUT), DIMENSION(NTM,2) :: TRLAKEL
      REAL*8, INTENT(IN), DIMENSION(NTM) :: TRUN0,TREVAP
      REAL*8, INTENT(OUT), DIMENSION(NTM) :: TRO,TRI
      REAL*8, DIMENSION(NTM) :: DTR2,TRUNO,TRUNI,TRF1,TRF2,FRAC
#ifdef TRACERS_SPECIAL_O18
      REAL*8 fracls
      INTEGER N
#endif
#endif
      REAL*8,PARAMETER :: emin=-1d-10  !@var emin min energy deficit required before ice forms (J/m^2)
      REAL*8 :: ENRGF1, ACEF1, ENRGF2, ACEF2, FHO, FHI, FH0, FH1, FH2, FSR2, &
                ENRGI, ENRGI2, ENRGO, ENRGO2, RUNO, RUNI, TLK2, DM2, DH2, FRATO, FRATI, E2O, E2I
      Character*300 :: out_line  !@var out_line local variable to hold mixed-type output for parallel I/O

!**** initiallize output
      ENRGFO=0. ; ACEFO=0. ; ACEFI=0. ; ENRGFI=0.

!**** Calculate heat and mass fluxes to lake
      ENRGO = FODT-SROX(1)*FSR2 ! in open water
      ENRGO2=     +SROX(1)*FSR2 ! in open water, second layer
      ENRGI = FIDT-SROX(2)*FSR2 ! under ice
      ENRGI2=     +SROX(2)*FSR2 ! under ice, second layer
      RUNO  =-EVAPO
      RUNI  = RUN0
#ifdef TRACERS_WATER
      TRUNO(:)=-TREVAP(:)
      TRUNI(:)= TRUN0(:)
      FRAC(:)=1.
#ifdef TRACERS_SPECIAL_O18
      FRAC(1:NTM) = fracls(1:NTM)  !  fractionation when freezing
#endif
#endif

!**** Bring up mass from second layer if required/allowed
      IF (MLAKE(1)+RUNO.lt.MINMLD*RHOW.and.MLAKE(2).gt.0) THEN
         DM2 = MIN(MLAKE(2),MINMLD*RHOW-(MLAKE(1)+RUNO))
         DH2 = DM2*(ELAKE(2)+(1.-ROICE)*ENRGO2+ROICE*ENRGI2)/MLAKE(2)
#ifdef TRACERS_WATER
         DTR2(:) = DM2*TRLAKEL(:,2)/MLAKE(2)
#endif
      ELSE
         DM2 = 0
         DH2 = 0
#ifdef TRACERS_WATER
         DTR2(:) = 0
#endif
      END IF

!**** Apply fluxes to 2nd layer
      IF (DM2.lt.MLAKE(2)) THEN
         MLAKE(2) = MLAKE(2) - DM2
         ELAKE(2) = ELAKE(2) - DH2 + (1.-ROICE)*ENRGO2 + ROICE*ENRGI2
#ifdef TRACERS_WATER
         TRLAKEL(:,2) = TRLAKEL(:,2) - DTR2(:)
#endif
      ELSE
         MLAKE(2) = 0
         ELAKE(2) = 0
#ifdef TRACERS_WATER
         TRLAKEL(:,2) = 0
#endif
      END IF

      E2O = 0  ;  E2I = 0

!**** Calculate energy in mixed layer (open ocean)
      IF (ROICE.LT.1d0) THEN
         FHO=ELAKE(1)+ENRGO+DH2-(MLAKE(1)+DM2+RUNO)*TFL*SHW
         IF (FHO.LT.emin) THEN ! FLUXES COOL WATER TO FREEZING, FORM ICE
            ACEFO = FHO/(TFL*(SHI-SHW)-LHM)
            ACEFO = MIN(ACEFO,MAX(MLAKE(1)+DM2+RUNO-MINMLD*RHOW,0d0))
            ENRGFO = ACEFO*(TFL*SHI-LHM)
            E2O=FHO-ENRGFO
         END IF
      END IF

      IF (ROICE.GT.0) THEN
!**** Calculate energy in mixed layer (under ice)
         FHI=ELAKE(1)+DH2+ENRGI-(MLAKE(1)+DM2+RUNI)*TFL*SHW
         IF (FHI.LT.emin) THEN ! FLUXES COOL WATER TO FREEZING, FORM ICE
            ACEFI =FHI/(TFL*(SHI-SHW)-LHM)
            ACEFI =MIN(ACEFI,MAX(MLAKE(1)+DM2+RUNI-MINMLD*RHOW,0d0))
            ENRGFI=ACEFI*(TFL*SHI-LHM)
            E2I=FHI-ENRGFI
         END IF
      END IF
#ifdef TRACERS_WATER
      TRO(:)=ACEFO*FRAC(:)*TRLAKEL(:,1)/MLAKE(1)
      TRI(:)=ACEFI*FRAC(:)*TRLAKEL(:,1)/MLAKE(1)
#endif

!**** Update first layer variables
      MLAKE(1) = MLAKE(1) + DM2 + (1-ROICE)*(RUNO -ACEFO)  + ROICE*(RUNI -ACEFI)
      ELAKE(1) = ELAKE(1) + DH2 + (1-ROICE)*(ENRGO-ENRGFO) + ROICE*(ENRGI-ENRGFI)
#ifdef TRACERS_WATER
      TRLAKEL(:,1)=TRLAKEL(:,1)+DTR2(:)+(1.-ROICE)*(TRUNO(:)-TRO(:))+ ROICE*(TRUNI(:)-TRI(:))
#endif

      ACEF1=0. ; ACEF2=0. ; ENRGF1=0. ; ENRGF2=0.
!**** Take remaining energy and start to freeze second layer
      FH2= ELAKE(1)-MLAKE(1)*TFL*SHW
      IF (FH2.LT.emin) THEN
         IF (MLAKE(2).gt.0) THEN
!**** FH2=-ACEF2*(TLK2-TFL)*SHW+ACEF2*LHM
            TLK2    =ELAKE(2)/(MLAKE(2)*SHW)
            ACEF2   =-FH2/(TLK2*SHW-TFL*SHI+LHM)
            ACEF2   =MIN(ACEF2,MLAKE(2))
            ENRGF2  =ACEF2*(TFL*SHI-LHM)
            ELAKE(1)=ELAKE(1)+ACEF2*TLK2*SHW-ENRGF2
            ELAKE(2)=ELAKE(2)-ACEF2*TLK2*SHW
            MLAKE(2)=MLAKE(2)-ACEF2
         END IF
         FH1=ELAKE(1)-MLAKE(1)*TFL*SHW
         IF (FH1.lt.emin) THEN      ! all layer 2 froze, freeze layer 1
            ACEF1=FH1/(TFL*(SHI-SHW)-LHM)
!**** limit freezing if lake would end up below 50cm depth
!**** this implies ice will be colder than TFL
            IF (MLAKE(1)-ACEF1.lt.0.5d0*RHOW) THEN
               ENRGF1=FH1
               ACEF1=MIN(MLAKE(1)-0.2d0*RHOW,MAX(0.4d0*MLAKE(1)+0.6d0*ACEF1 -0.2d0*RHOW,0d0))
!**** force ice temp > -95 C. Sometimes occurs w. extreme super-cooled water
               IF (ENRGF1.lt.ACEF1*(-95.*SHI-LHM)) ACEF1=ENRGF1/(-95.*SHI -LHM)
               if (qcheck) then
                  WRITE(out_line,*) "Min lake level reached during frez: rsi,hlake,"// &
                                    "hice,tlk,tice",i0,j0,roice,mlake(1)/rhow,ACEF1/RHOW &
                                    ,elake(1)/(mlake(1)*shw),(enrgf1/acef1+lhm)/shi
                  CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
               END IF
               IF (ACEF1.gt.MLAKE(1)) then ! water was too cold - nothing possible
                  Write (6,*) "Lk. freeze impossible: i,j,Tw,Ti,hic,hlk = ",i0,j0,elake(1)/(mlake(1)*shw), &
                              (fh1/acef1+lhm)/shi, acef1/rhow,mlake(1)/rhow
                  call stop_model("Lake water too cold during LKSOURC",255)
               END IF
            ELSE
               ENRGF1=ACEF1*(TFL*SHI-LHM)
            END IF
            MLAKE(1)=MLAKE(1)-ACEF1
            ELAKE(1)=MLAKE(1)*TFL*SHW
         END IF
      END IF
#ifdef TRACERS_WATER
      TRF1(:) = ACEF1*FRAC(:)*TRLAKEL(:,1)/(MLAKE(1)+ACEF1)
      TRLAKEL(:,1)=TRLAKEL(:,1)-TRF1(:)
      IF (MLAKE(2).gt.0) THEN
         TRF2(:) = MIN(ACEF2*FRAC(:)/(MLAKE(2)+ACEF2),1d0)*TRLAKEL(:,2)
         TRLAKEL(:,2)=TRLAKEL(:,2)-TRF2(:)
      ELSE ! possibility of complete freezing (and so no frac)
         TRF2(:) = TRLAKEL(:,2)
         TRLAKEL(:,2) = 0.
      END IF
#endif

!**** combine mass and energy fluxes for output
!**** Note that output fluxes are over open water/ice covered fractions
!**** distribute ice fluxes according to flux amounts
      FRATO = 1d0
      FRATI = 1d0
      IF (E2I+E2O.lt.0) THEN
         FRATO = E2O/(E2I*ROICE+E2O*(1.-ROICE))
         FRATI = E2I/(E2I*ROICE+E2O*(1.-ROICE))
      END IF
      ACEFO =ACEFO + (ACEF1 +ACEF2 )*FRATO
      ACEFI =ACEFI + (ACEF1 +ACEF2 )*FRATI
      ENRGFO=ENRGFO+ (ENRGF1+ENRGF2)*FRATO
      ENRGFI=ENRGFI+ (ENRGF1+ENRGF2)*FRATI
#ifdef TRACERS_WATER
      TRO(:)=TRO(:) + (TRF1(:) + TRF2(:))* FRATO
      TRI(:)=TRI(:) + (TRF1(:) + TRF2(:))* FRATI
#endif
      END SUBROUTINE LKSOURC



      SUBROUTINE LKMIX (MLAKE,ELAKE, &
#ifdef TRACERS_WATER
                        TRLAKEL, &
#endif
                        HLAKE,TKE,ROICE,DTSRC)
!@sum  LKMIX calculates mixing and entrainment in lakes
!@auth Gavin Schmidt
      IMPLICIT NONE
!@var MLAKE,ELAKE mass and energy in lake layers (kg,J /m^2)
      REAL*8, INTENT(INOUT), DIMENSION(2) :: MLAKE,ELAKE
!@var TKE turbulent kinetic energy input at surface of lake (J/m^2)
!@var ROICE ice fraction
      REAL*8, INTENT(IN) :: TKE,ROICE
!@var HLAKE sill depth for lake (m)
      REAL*8, INTENT(IN) :: HLAKE
!@var DTSRC source time step (s)
      REAL*8, INTENT(IN) :: DTSRC
#ifdef TRACERS_WATER
!@var TRLAKEL tracer mass in lake layers (kg/m^2)
      REAL*8, INTENT(INOUT), DIMENSION(NTM,2) :: TRLAKEL
      REAL*8, DIMENSION(NTM) :: DTML,TR1N,TR2N,TRLT
#endif
!@param MAXRHO,RHO0,BFAC freshwater density function approximation
      REAL*8, PARAMETER :: MAXRHO=1d3, RHO0=999.842594d0, BFAC=(MAXRHO-RHO0)/16d0

      REAL*8 TLK1, TLK2, HLT, MLT, DTK, E1N, E2N, ATKE, H1, H2, DRHO, DML, DHML

!**** Only mix if there is a second layer!
      IF (MLAKE(2).gt.0) THEN
        TLK1=ELAKE(1)/(MLAKE(1)*SHW)
        TLK2=ELAKE(2)/(MLAKE(2)*SHW)
        HLT=ELAKE(1)+ELAKE(2)
        MLT=MLAKE(1)+MLAKE(2)
#ifdef TRACERS_WATER
        TRLT(:)=TRLAKEL(:,1)+TRLAKEL(:,2)
#endif
!**** Test for static stability
!**** DRHO=RHO(TLK2)-RHO(TLK1)~=(TLK2-TLK1)*dRHOdT((TLK1+TLK2)/2)
!**** Assumes a parabolic density function going through MAXRHO at
!**** TMAXRHO, and RHO0 at T=0. (reasonable up to about 12 C)
        IF ((TMAXRHO-0.5*(TLK1+TLK2))*(TLK2-TLK1).lt.0) THEN
!**** mix uniformly and set MLD to minimum
          MLAKE(1)=MIN(MLT,MAX(MINMLD*RHOW,MLT-HLAKE*RHOW))
          MLAKE(2)=MLT-MLAKE(1)
          ELAKE(1)=HLT*MLAKE(1)/MLT
          ELAKE(2)=HLT*MLAKE(2)/MLT
#ifdef TRACERS_WATER
          TRLAKEL(:,1)=TRLT(:)*MLAKE(1)/MLT
          TRLAKEL(:,2)=TRLT(:)*MLAKE(2)/MLT
#endif
        ELSE ! not unstable, implicitly diffuse heat + entrain
!**** reduce mixing if there is ice cover
          DTK=2.*KVLAKE*(1.-ROICE)*DTSRC*RHOW**2
          E1N=(ELAKE(1)+DTK*HLT/(MLT*MLAKE(2)))/ (1.+DTK/(MLAKE(1)*MLAKE(2)))
          E2N=(ELAKE(2)+DTK*HLT/(MLT*MLAKE(1)))/ (1.+DTK/(MLAKE(1)*MLAKE(2)))
          ELAKE(1)=E1N
          ELAKE(2)=E2N
#ifdef TRACERS_WATER
!**** diffuse tracers using same KV as for heat?
          TR1N(:)=(TRLAKEL(:,1)+DTK*TRLT(:)/(MLT*MLAKE(2)))/ (1.+DTK/(MLAKE(1)*MLAKE(2)))
          TR2N(:)=(TRLAKEL(:,2)+DTK*TRLT(:)/(MLT*MLAKE(1)))/ (1.+DTK/(MLAKE(1)*MLAKE(2)))
          TRLAKEL(:,1)=TR1N(:)
          TRLAKEL(:,2)=TR2N(:)
#endif
!**** entrain deep water if there is available TKE
!**** take a factor of TKE and calculate change in PE
          IF (TKE.gt.0) THEN
            ATKE=0.2d0*TKE      ! 20% of TKE is available for mixing
            H1=MLAKE(1)/RHOW
            H2=MLAKE(2)/RHOW
            DRHO=(TLK2-TLK1)*2d0*BFAC*(TMAXRHO-0.5*(TLK1+TLK2))
            DML=ATKE*BYGRAV/(DRHO*0.5*H1)
            IF (DML*RHOW.lt.MLAKE(2)) THEN
              DHML=DML*ELAKE(2)/H2
              ELAKE(1)=ELAKE(1)+DHML
              ELAKE(2)=ELAKE(2)-DHML
              MLAKE(1)=MLAKE(1)+DML*RHOW
              MLAKE(2)=MLAKE(2)-DML*RHOW
#ifdef TRACERS_WATER
              DTML(:)=DML*TRLAKEL(:,2)/H2
              TRLAKEL(:,1)=TRLAKEL(:,1)+DTML(:)
              TRLAKEL(:,2)=TRLAKEL(:,2)-DTML(:)
#endif
            ELSE                ! entire second layer is entrained
              MLAKE(1)=MLT
              MLAKE(2)=0
              ELAKE(1)=HLT
              ELAKE(2)=0
#ifdef TRACERS_WATER
              TRLAKEL(:,1)=TRLT(:)
              TRLAKEL(:,2)=0.
#endif
            END IF
          END IF
        END IF
      END IF
      RETURN
!****
      END SUBROUTINE LKMIX


      END MODULE LAKES


      SUBROUTINE ALLOC_LAKES (GRID)
!@SUM  To alllocate arrays whose sizes now need to be determined at run-time
!@auth Raul Garza-Robles
      Use DOMAIN_DECOMP_ATM, Only: DIST_GRID, getDomainBounds
      Use LAKES, Only: dZSG, RATE, DHORZ,KDIREC,IFLOW,JFLOW, KD911,IFL911,JFL911
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      integer :: i_0h,i_1h,j_0h,j_1h

      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO

      ALLOCATE ( KDIREC (I_0H:I_1H,J_0H:J_1H), &
                  KD911 (I_0H:I_1H,J_0H:J_1H), &
                  IFLOW (I_0H:I_1H,J_0H:J_1H), &
                  JFLOW (I_0H:I_1H,J_0H:J_1H), &
                  IFL911(I_0H:I_1H,J_0H:J_1H), &
                  JFL911(I_0H:I_1H,J_0H:J_1H), &
                  dZSG  (I_0H:I_1H,J_0H:J_1H), &
                  RATE  (I_0H:I_1H,J_0H:J_1H), &
                  DHORZ (I_0H:I_1H,J_0H:J_1H))
      END SUBROUTINE ALLOC_LAKES



      Subroutine INIT_LAKES (inilake,istart)
!@sum  init_LAKES initiallises lake variables
!@auth Gary L. Russell / Gavin Schmidt
      Use FILEMANAGER
      Use CONSTANT,   Only: TWOPI,RADIUS,rhow,shw,tf,pi,grav
      Use RESOLUTION, Only: im,jm
      Use MODEL_COM,  Only: dtsrc
      Use ATM_COM,    Only: zatmo, traditional_coldstart_aic
      Use GEOM,       Only: axyp,imaxj,lonlat_to_ij,lon2d_dg,lat2d_dg
      Use FLUXES,     Only: atmocn,atmsrf ,flice,focean,fearth0,flake0,fland
      Use GHY_COM,    Only: fearth
      Use LAKES
      Use LAKES_COM
      Use DIAG_COM,   Only: npts,conpt0,icon_LKM,icon_LKE
      Use Dictionary_mod
      Use PARIO,      Only: par_open, par_close, read_dist_data, Variable_Exists
      Use DOMAIN_DECOMP_ATM, Only: GRID,WRITE_PARALLEL, getDomainBounds,HALO_UPDATE, am_i_root
#ifdef TRACERS_WATER
      Use OldTracer_mod, Only: trw0
#endif
#ifdef SCM
      Use SCM_COM, Only: SCMopt,SCMin
#endif

      IMPLICIT NONE
      Logical,Intent(InOut) :: inilake
      Integer,Intent(In)    :: ISTART
!**** Local variables
      Integer :: FROM, J_0, J_1, J_0H, J_1H, J_0S, J_1S, I_0, I_1, I_0H, I_1H, JMIN_FILL, JMAX_FILL, iu_warn, &
                 I,J,IU,JU,ID,JD,INM, &  !@var I,J,IU,JU,ID,JD loop variables
                 fid,ios, iloop_min,iloop_max,jloop_min,jloop_max, get_dir, &
                 iu_RVR  !@var iu_RVR unit number for river direction file
      Logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE, QCON(NPTS), T=.TRUE., F=.FALSE., &
                 QEXIST, &  !  whether dZSG is provided on atmospheric TOPO file 
                 iwrap  !@var iwrap true if I direction is periodic and has no halo
      Real*8  :: ZSOLID,WT,ZMEAN,ZSQ, dZSGmin,yAREA11, dZSG11 = 10
      Real*8  :: SPMIN,SPMAX,SPEED0,SPEED,DZDH,DZDH1,MLK1,fac,fac1, &
                 horzdist_2pts  !  external function for now
      REAL*8, allocatable, dimension(:,:) :: down_lat_loc ,down_lon_loc,down_lat_911_loc,down_lon_911_loc
      INTEGER, dimension(2) :: ij
      REAL*8, dimension(2) :: ll
      REAL*4, dimension(nrvrmx) :: lat_rvr,lon_rvr
      REAL*8, DIMENSION(:,:), POINTER :: GTEMP,GTEMP2,GTEMPR
      Character*10  :: CONPT(NPTS)
      Character*80  :: fmtstr, TITLEI
      Character*300 :: out_line  !@var out_line local variable to hold mixed-type output for parallel I/O
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:), POINTER :: GTRACER
#endif

      GTEMP => ATMOCN%GTEMP
      GTEMP2 => ATMOCN%GTEMP2
      GTEMPR => ATMOCN%GTEMPR
#ifdef TRACERS_WATER
      GTRACER => ATMOCN%GTRACER
#endif

      call getDomainBounds(GRID, J_STRT = J_0,  J_STOP = J_1, &
                             J_STRT_SKP = J_0S, J_STOP_SKP = J_1S, &
                             J_STRT_HALO= J_0H, J_STOP_HALO= J_1H, &
                           HAVE_SOUTH_POLE = HAVE_SOUTH_POLE, &
                           HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      iwrap = I_0H == I_0

!****
!**** LAKECB  MWL      Mass of water in lake (kg)
!****         GML      Liquid lake enthalpy (J)
!****         TLAKE    Temperature of lake surface (C)
!****         HLAKE    Lake sill depth (m)
!****         ALAKE    Lake surface area (m^2)
!****         dZSG     Range of Solid Ground Topography in cell (m)
!****
!**** FIXDCB  FLAKE0    Original Lake fraction (1)
!****

!**** Get parameters from rundeck
      call sync_param("river_fac",river_fac)
      call sync_param("init_flake",init_flake)
      if (init_flake == 1 .and. istart < 9) inilake = .true.
      call sync_param("variable_lk",variable_lk)
      call sync_param("lake_rise_max",lake_rise_max)
      call sync_param("lake_ice_max" ,lake_ice_max)

!**** Read dZLAKE and dZSG from TOPO input file
      fid = par_open(grid,'TOPO','read')
      call read_dist_data(grid,fid,'hlake',hlake)
      QEXIST = Variable_Exists (GRID,FID,'dZSG')
      If (QEXIST)  Call Read_Dist_Data (GRID,FID,'dZSG' ,dZSG)
      call par_close(grid,fid)

      If (QEXIST)  GoTo 100
!****
!**** If dZSG is not provided on atmospheric TOPO file, it is derived here
!**** ZSOLID = solid ground topography = ZATMO/GRAV - FLAKE0*HLAKE - FOCEAN*HOCEAN = ZATMO/GRAV - (FLAKE0+FOCEAN)*HLAKE
!**** dZSG = range of ZSOLID in grid cell = StanDev(neighboring ZSOLID) * 2
!**** Minimum value of dZSG in a cell is roughly equal to minimum dZSG11 of 1x1 cell times square root of ratio of areas:
!**** Minimum dZSGXY = dZSG11 * (AREAXY / AREA11)^.5 = 10 * [AREAXY / (TWOPI/360)^2]^.5
!****
      If (AM_I_ROOT())  Write (6,*) 'LAKES2: dZSG is derived from ZATMO,FLAKE0,HLAKE'
      yAREA11 = 1 / (TWOPI*RADIUS/360)**2
!**** Lat-Lon south pole cells
      If (Have_South_Pole) Then
         If (FOCEAN(1,1) > 0) &
            Then  ;  dZSG(:,1) = 0
            Else  ;  ZSOLID = ZATMO(1,1)/GRAV - FLAKE0(1,1)*HLAKE(1,1) 
                     WT     = .25*IM
                     ZMEAN  = .25*IM * ZSOLID
                     ZSQ    = .25*IM * ZSOLID**2 
                     Do I=1,IM
                        ZSOLID = ZATMO(I,2)/GRAV - (FLAKE0(I,2)+FOCEAN(I,2))*HLAKE(I,2)
                        WT     = WT    + (1-FOCEAN(I,2))
                        ZMEAN  = ZMEAN + (1-FOCEAN(I,2)) * ZSOLID
                        ZSQ    = ZSQ   + (1-FOCEAN(I,2)) * ZSOLID**2  ;  EndDo
                     ZMEAN = ZMEAN / WT
                     dZSGmin = dZSG11 * (AXYP(1,1)*yAREA11)**.5
                     dZSG(1,1) = Max (Sqrt(ZSQ/WT - ZMEAN**2)*2, dZSGmin)
                     dZSG(:,1) = dZSG(1,1)  ;  EndIf  ;  EndIf
!**** Lat-Lon north pole cells
      If (Have_North_Pole) Then
         If (FOCEAN(1,JM) > 0) &
            Then  ;  dZSG(:,JM) = 0
            Else  ;  ZSOLID = ZATMO(1,JM)/GRAV - FLAKE0(1,JM)*HLAKE(1,JM)
                     WT     = .25*IM
                     ZMEAN  = .25*IM * ZSOLID
                     ZSQ    = .25*IM * ZSOLID**2
                     Do I=1,IM
                        ZSOLID = ZATMO(I,JM-1)/GRAV - (FLAKE0(I,JM-1)+FOCEAN(I,JM-1))*HLAKE(I,JM-1)
                        WT     = WT    + (1-FOCEAN(I,JM-1))
                        ZMEAN  = ZMEAN + (1-FOCEAN(I,JM-1)) * ZSOLID
                        ZSQ    = ZSQ   + (1-FOCEAN(I,JM-1)) * ZSOLID**2  ;  EndDo
                     ZMEAN = ZMEAN / WT
                     dZSGmin = dZSG11 * (AXYP(1,JM)*yAREA11)**.5
                     dZSG(1,JM) = Max (Sqrt(ZSQ/WT - ZMEAN**2)*2, dZSGmin)
                     dZSG(:,JM) = dZSG(1,JM)  ;  EndIf  ;  Endif
!**** Non-polar cells
      Do J=J_0S,J_1S  ;  Do I=I_0,I_1
         If (FOCEAN(I,J) > 0) &
            Then  ;  dZSG(I,J) = 0
            Else  ;  ZSOLID = ZATMO(I,J)/GRAV - FLAKE0(I,J)*HLAKE(I,J)
                     WT     = 1
                     ZMEAN  = ZSOLID
                     ZSQ    = ZSOLID**2
                     ZSOLID = ZATMO(I,J-1)/GRAV - (FLAKE0(I,J-1)+FOCEAN(I,J-1))*HLAKE(I,J-1)
                     WT     = WT    + (1-FOCEAN(I,J-1))                        
                     ZMEAN  = ZMEAN + (1-FOCEAN(I,J-1)) * ZSOLID
                     ZSQ    = ZSQ   + (1-FOCEAN(I,J-1)) * ZSOLID**2
                     ZSOLID = ZATMO(I,J+1)/GRAV - (FLAKE0(I,J+1)+FOCEAN(I,J+1))*HLAKE(I,J+1)
                     WT     = WT    + (1-FOCEAN(I,J+1))                        
                     ZMEAN  = ZMEAN + (1-FOCEAN(I,J+1)) * ZSOLID
                     ZSQ    = ZSQ   + (1-FOCEAN(I,J+1)) * ZSOLID**2
                     If (IWRAP .and. I==1) &
                        Then  ;  ZSOLID = ZATMO(IM,J)/GRAV - (FLAKE0(IM,J)+FOCEAN(IM,J))*HLAKE(IM,J)
                                 WT     = WT    + (1-FOCEAN(IM,J))
                                 ZMEAN  = ZMEAN + (1-FOCEAN(IM,J)) * ZSOLID
                                 ZSQ    = ZSQ   + (1-FOCEAN(IM,J)) * ZSOLID**2
                        Else  ;  ZSOLID = ZATMO(I-1,J)/GRAV - (FLAKE0(I-1,J)+FOCEAN(I-1,J))*HLAKE(I-1,J)
                                 WT     = WT    + (1-FOCEAN(I-1,J))
                                 ZMEAN  = ZMEAN + (1-FOCEAN(I-1,J)) * ZSOLID
                                 ZSQ    = ZSQ   + (1-FOCEAN(I-1,J)) * ZSOLID**2  ;  EndIf
                     If (IWRAP .and. I==IM) &
                        Then  ;  ZSOLID = ZATMO(1,J)/GRAV - (FLAKE0(1,J)+FOCEAN(1,J))*HLAKE(1,J)
                                 WT     = WT    + (1-FOCEAN(1,J))
                                 ZMEAN  = ZMEAN + (1-FOCEAN(1,J)) * ZSOLID
                                 ZSQ    = ZSQ   + (1-FOCEAN(1,J)) * ZSOLID**2
                        Else  ;  ZSOLID = ZATMO(I+1,J)/GRAV - (FLAKE0(I+1,J)+FOCEAN(I+1,J))*HLAKE(I+1,J)
                                 WT     = WT    + (1-FOCEAN(I+1,J))
                                 ZMEAN  = ZMEAN + (1-FOCEAN(I+1,J)) * ZSOLID
                                 ZSQ    = ZSQ   + (1-FOCEAN(I+1,J)) * ZSOLID**2  ;  EndIf
                     ZMEAN = ZMEAN / WT
                     dZSGmin = dZSG11 * (AXYP(I,J)*yAREA11)**.5
                     dZSG(I,J) = Max (Sqrt(ZSQ/WT - ZMEAN**2)*2, dZSGmin)  ;  EndIf  ;  EndDo  ;  EndDo

!**** initialise FLAKE if requested (i.e. from older restart files)
  100 If (INILAKE) Then
         print*,"Initialising FLAKE from TOPO file..."
         FLAKE = FLAKE0
      end if

!**** Ensure that HLAKE is a minimum of 1m for FLAKE>0
      call openunit("warn_lakes",iu_warn)
      DO J=J_0,J_1  ;  DO I=I_0, I_1
         IF (FLAKE0(I,J)+FLAKE(I,J).gt.0 .and. HLAKE(I,J).lt.HLAKE_MIN)  THEN
            write(iu_warn,*) "Warning: Fixing HLAKE",i,j,FLAKE(I,J), FLAKE0(I,J),HLAKE(I,J),"--> ",HLAKE_MIN," m"
            HLAKE(I,J)=HLAKE_MIN
         END IF
      END DO  ;  END DO
      call closeunit(iu_warn)

      IF (INILAKE) THEN
!**** Set lake variables from surface temperature
!**** This is just an estimate for the initiallisation
         if(istart==2) then ! pbl has not been initialized yet
           ! todo: get this temperature via other means
           if(traditional_coldstart_aic)  call read_pbl_tsurf_from_nmcfile
         endif
         DO J=J_0,J_1  ;  DO I=I_0,I_1
            IF (FLAKE(I,J).gt.0) THEN
               TLAKE(I,J) = MAX(0d0,atmsrf%TSAVG(I,J)-TF)
               MWL(I,J)   = RHOW*HLAKE(I,J)*FLAKE(I,J)*AXYP(I,J)
               MLK1       = MINMLD*RHOW*FLAKE(I,J)*AXYP(I,J)
               GML(I,J)   = SHW*(MLK1*TLAKE(I,J) +(MWL(I,J)-MLK1)*MAX(TLAKE(I,J),4d0))
               MLDLK(I,J) = MINMLD
#ifdef TRACERS_WATER
               TRLAKE(:,1,I,J)=MLK1*TRW0()
               TRLAKE(:,2,I,J)=(MWL(I,J)-MLK1)*TRW0()
#endif
            ELSE
               TLAKE(I,J) = 0.
               MWL(I,J)   = 0.
               GML(I,J)   = 0.
               MLDLK(I,J) = MINMLD
#ifdef TRACERS_WATER
               TRLAKE(:,:,I,J)=0.
#endif
            END IF
         END DO  ;  END DO
      END IF

      IF (init_flake.eq.2 .and. istart.le.9) THEN
         print*,"Checking for excess water..."
         DO J=J_0,J_1  ;  DO I=I_0,I_1
            if (FLAKE(I,J) > .949d0*(FLAKE(I,J)+FEARTH(I,J)) .and. &
                MWL(I,J) > (HLAKE(I,J)+LAKE_RISE_MAX)*FLAKE(I,J)*RHOW*AXYP(I,J))  then
               print*,"Adjusting lake level:",i,j," from ",MWL(I,J)/(FLAKE(I,J)*RHOW*AXYP(I,J)), &
                      " to ",HLAKE(I,J) +LAKE_RISE_MAX,MLDLK(I,J)
               fac=(MWL(I,J)-(HLAKE(I,J)+LAKE_RISE_MAX)*FLAKE(I,J)*RHOW*AXYP(I,J))  ! excess mass fractional loss of layer 1 mass
               fac1=fac/(MLDLK(I,J)*FLAKE(I,J)*RHOW*AXYP(I,J))
               if (fac1.lt.1) then
#ifdef TRACERS_WATER
                  TRLAKE(:,1,I,J)=TRLAKE(:,1,I,J)*(1d0-fac1)
#endif
                  MLDLK(I,J)=MLDLK(I,J)*(1-fac1)
                  GML(I,J)=GML(I,J)-fac*SHW*TLAKE(I,J)
                  MWL(I,J)=MWL(I,J)-fac
               else
                   call stop_model('INIT_LAKE: Inconsistent ml depth in lakes',255)
               end if
            END IF
         END DO  ;  END DO
      END IF

!**** Set fixed topographic variables when FOCEAN == 0
!**** Solid ground topography is a linear function of cell area or FLAKE
!**** Lake fills like a triangle from lowest solid ground topography
!**** Vertical cross section of circular lake is then a concave parabola
!**** dZSG  = fixed range of continental solid ground topography (m) in grid cell
!**** ZLBOT = fixed lowest altitude of solid ground topography (m) in grid cell = ZSOLID - .5*dZSG
!**** HLAKE = [ZLTOP - ZLBOT] / 2 = FLAKE * dZSG / 2 = average lake thickness (m)
!**** FLAKE = 2 * HLAKE / dZSG = horizontal lake fraction of cell area
!**** MWL   = RHOW * AXYP * HLAKE * FLAKE = RHOW * AXYP * HLAKE * 2 * HLAKE / dZSG = liquid lake mass (kg)
!**** HLAKE = Sqrt [MWL * dZSG / 2 * RHOW * AXYP]
!**** If HLAKE < HLAKE_MIN, then HLAKEnew = HLAKE_MIN and FLAKEnew = MWL / RHOW * AXYP * HLAKEnew
      DO J=J_0,J_1  ;  DO I=I_0,I_1
         If (FLAKE(I,J) > 0 .and. HLAKE(I,J) < HLAKE_MIN) Then
             FLAKE(I,J) = FLAKE(I,J) * HLAKE(I,J) / HLAKE_MIN
             HLAKE(I,J) = HLAKE_MIN  
         END IF
      END DO  ;  END DO

      CALL PRINTLK("IN")
!**** Set GTEMP arrays for lakes
      IF (ISTART.gt.0) THEN
         DO J=J_0,J_1  ;  DO I=I_0,I_1
            IF (FLAKE(I,J).gt.0) THEN
               GTEMP(I,J)=TLAKE(I,J)
               GTEMPR(I,J) =TLAKE(I,J)+TF
#ifdef SCM
               if (SCMopt%Tskin) then
                  GTEMP(I,J) = SCMin%Tskin - TF
                  GTEMPR(I,J) = SCMin%Tskin
               endif
#endif
               IF (MWL(I,J).gt.(1d-10+MLDLK(I,J))*RHOW*FLAKE(I,J)*AXYP(I,J))  THEN
                  GTEMP2(I,J) = (GML(I,J) - TLAKE(I,J)*SHW*MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J)) / &
                                (SHW * (MWL(I,J) - MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J)))
!**** If starting from a possibly corrupted rsf file, check Tlk2
                  IF (GTEMP2(I,J)>TLAKE(I,J)+1.and.GTEMP2(I,J)>10 .and. istart<9) THEN
                     WRITE(6,*) "Warning: Unphysical Tlk2 fixed",I,J, GTEMP(I,J),GTEMP2(I,J)
                     GTEMP2(I,J)=GTEMP(I,J)  ! set to Tlk1
                     GML(I,J)=TLAKE(I,J)*SHW*MWL(I,J)
                  END IF
               ELSE
                  GTEMP2(I,J)=TLAKE(I,J)
               END IF
#ifdef SCM
               if (SCMopt%Tskin) GTEMP2(I,J) = GTEMP(I,J)
#endif
#ifdef TRACERS_WATER
               GTRACER(:,I,J)=TRLAKE(:,1,I,J)/(MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J))
#endif
               atmocn%MLHC(I,J)= SHW*MLDLK(I,J)*RHOW
            END IF
         END DO  ;  END DO
      END IF

!**** setting river directions
      IFLOW=-99.  ; JFLOW=-99.
      IFL911=-99. ; JFL911=-99.
      KDIREC=0 ; KD911=0
      RATE=0.  ; nrvr=0

      allocate(down_lat_loc(I_0H:I_1H,J_0H:J_1H), &
               down_lon_loc(I_0H:I_1H,J_0H:J_1H), &
               down_lat_911_loc(I_0H:I_1H,J_0H:J_1H), &
               down_lon_911_loc(I_0H:I_1H,J_0H:J_1H) )

!**** Read named river mouth positions
      nrvr = 0
      if (file_exists('NAMERVR')) then
         call openunit("NAMERVR",iu_RVR,.false.,.true.)
         read(iu_RVR,*)
         read(iu_RVR,'(a)') fmtstr
         do
            read(iu_RVR,trim(fmtstr),iostat=ios) namervr(nrvr+1),lat_rvr(nrvr+1),lon_rvr(nrvr+1)
            if(ios.ne.0) exit
            nrvr = nrvr + 1
         enddo
         call closeunit(iu_RVR)

         WRITE (out_line,*) 'Named river file read '
         CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
      endif

!**** Read in down stream lat/lon positions
      if (file_exists('RVR')) then
         fid = par_open(grid,'RVR','read')
         call read_dist_data(grid,fid,'down_lat',down_lat_loc)
         call read_dist_data(grid,fid,'down_lon',down_lon_loc)
         call read_dist_data(grid,fid,'down_lat_911',down_lat_911_loc)
         call read_dist_data(grid,fid,'down_lon_911',down_lon_911_loc)
         call par_close(grid,fid)
      else
         do j=j_0,j_1  ;  do i=i_0,i_1
            down_lon_loc(i,j) = lon2d_dg(i,j)
            down_lon_911_loc(i,j) = lon2d_dg(i,j)
            down_lat_loc(i,j) = lat2d_dg(i,j)
            down_lat_911_loc(i,j) = lat2d_dg(i,j)
         enddo  ;  enddo
      endif
      CALL HALO_UPDATE(GRID, down_lat_loc)
      CALL HALO_UPDATE(GRID, down_lon_loc)
      CALL HALO_UPDATE(GRID, down_lat_911_loc)
      CALL HALO_UPDATE(GRID, down_lon_911_loc)

!****
!**** From each box calculate the downstream river box
!****
      if(have_south_pole) then
         jloop_min=1
      else
         jloop_min=j_0h
      endif
      if(have_north_pole) then
         jloop_max=jm
      else
         jloop_max=j_1h
      endif
      do j=jloop_min,jloop_max
         iloop_min=i_0h
         iloop_max=i_1h
         if(j.lt.1 .or. j.gt.jm) then
! avoid nonexistent halo corners of a cube face.
            iloop_min=max(iloop_min,1)
            iloop_max=min(iloop_max,im)
         endif
         do i=iloop_min,iloop_max
            if (down_lon_loc(i,j).gt.-1000.) then
               ll(1)=down_lon_loc(i,j)
               ll(2)=down_lat_loc(i,j)
               call lonlat_to_ij(ll,ij)
               IFLOW(I,J)=ij(1) ; JFLOW(I,J)=ij(2)
               if(iflow(i,j).ge.i_0h .and. iflow(i,j).le.i_1h .and. jflow(i,j).ge.j_0h .and. jflow(i,j).le.j_1h) &
                  DHORZ(I,J) = horzdist_2pts(i,j,iflow(i,j),jflow(i,j))
            ElseIf (DOWN_LON_LOC(I,J) == -9999)  Then  ;  KDIREC(I,J) = 9
            else  ! if land but no ocean, print warning
               IF ((FEARTH0(I,J)+FLICE(I,J)+FLAKE0(I,J).gt.0) .and. FOCEAN(I,J).le.0 ) THEN
                  WRITE(6,*) "Land box has no river direction I,J: ",I,J,FOCEAN(I,J),FLICE(I,J),FLAKE0(I,J),FEARTH0(I,J)
               ELSE
                  IFLOW(I,J)=i ; JFLOW(I,J)=j     ! local flow
               END IF
               DHORZ(I,J) = horzdist_2pts(i,j,i,j)
            end if
            if (down_lon_911_loc(i,j).gt.-1000.) then
               ll(1)=down_lon_911_loc(i,j)
               ll(2)=down_lat_911_loc(i,j)
               call lonlat_to_ij(ll,ij)
               IFL911(I,J)=ij(1) ; JFL911(I,J)=ij(2)
            endif
!**** do we need get_dir? maybe only need to set KD=0 or >0?
            IF (IFLOW(I,J).gt. -99) KDIREC(I,J)=get_dir(I,J,IFLOW(I,J),JFLOW(I,J),IM,JM)
            IF (IFL911(I,J).gt.-99) KD911(I,J)=get_dir(I,J,IFL911(I,J),JFL911(I,J),IM,JM)
         END DO
      END DO

!**** define river mouths
      do inm=1,nrvr
         ll(1)=lon_rvr(inm)
         ll(2)=lat_rvr(inm)
         call lonlat_to_ij(ll,ij)

         IRVRMTH(INM)=ij(1) ; JRVRMTH(INM)=ij(2)
!        write(*,*) "rvr->",namervr(inm),ij(1),ij(2)

         if (IRVRMTH(INM).ge.I_0H .and. IRVRMTH(INM).le.I_1H .and. JRVRMTH(INM).ge.J_0H .and. JRVRMTH(INM).le.J_1H) THEN
            IF (FOCEAN(IRVRMTH(INM),JRVRMTH(INM)).le.0) &
               WRITE (6,*) "Warning: Named river outlet must be in ocean" &
                           ,INM,IRVRMTH(INM),JRVRMTH(INM),NAMERVR(INM) &
                           ,FOCEAN(IRVRMTH(INM),JRVRMTH(INM)),FLICE(IRVRMTH(INM) &
                           ,JRVRMTH(INM)),FLAKE0(IRVRMTH(INM),JRVRMTH(INM)) &
                           ,FEARTH0(IRVRMTH(INM),JRVRMTH(INM))
         end if
      end do

!****
!**** Calculate river flow RATE (per source time step)
!****
      SPEED0= .35d0  ! m/s
      SPMIN = .15d0  ! m/s
      SPMAX = 5.     ! m/s
      DZDH1 = .00005 ! ratio
      DO JU = J_0,J_1S  ;  DO IU=I_0,IMAXJ(JU)
         If (FOCEAN(IU,JU) < 1 .and. KDIREC(IU,JU) <= 8)  Then
            If (KDIREC(IU,JU) >= 1)  Then
               JD=JFLOW(IU,JU)
               ID=IFLOW(IU,JU)
               DZDH  = (ZATMO(IU,JU)-ZATMO(ID,JD)) / (GRAV*DHORZ(IU,JU))
            ELSE
               DZDH  = ZATMO(IU,JU) / (GRAV*DHORZ(IU,JU))
            END IF
            SPEED = SPEED0*DZDH/DZDH1
            IF(SPEED.lt.SPMIN)  SPEED = SPMIN
            IF(SPEED.gt.SPMAX)  SPEED = SPMAX
            RATE(IU,JU) = DTsrc*SPEED/DHORZ(IU,JU)
         END IF
      END DO  ;  END DO

      do j=j_0,j_1  ;  do i=i_0,imaxj(j)
         if (flake(i,j).gt.0.) then
            DLAKE(I,J)=MWL(I,J)/(RHOW*FLAKE(I,J)*AXYP(I,J))
            GLAKE(I,J)=GML(I,J)/(FLAKE(I,J)*AXYP(I,J))
         else
            DLAKE(I,J)=0.
            GLAKE(I,J)=0.
         endif
      enddo  ;  enddo

!**** assume that at the start GHY is in balance with LAKES
      SVFLAKE = FLAKE

!**** Make sure that constraints are satisfied by defining FLAND/FEARTH
!**** as residual terms.
      DO J=J_0,J_1  ;  DO I=I_0,IMAXJ(J)
!!       FLAND(I,J)=1.-FOCEAN(I,J)  !! already set if FOCEAN>0
         IF (FOCEAN(I,J).le.0) THEN
            FLAND(I,J)=1
            IF (FLAKE(I,J).gt.0) FLAND(I,J)=1.-FLAKE(I,J)
         END IF
         FEARTH(I,J)=FLAND(I,J)-FLICE(I,J) ! Earth fraction
      END DO  ;  END DO
      If (HAVE_SOUTH_POLE) Then
         FLAND(2:IM,1)=FLAND(1,1)
         FEARTH(2:IM,1)=FEARTH(1,1)
      End If
      If (HAVE_NORTH_POLE) Then
         FLAND(2:IM,JM)=FLAND(1,JM)
         FEARTH(2:IM,JM)=FEARTH(1,JM)
      End If

!**** Set conservation diagnostics for Lake mass and energy
      CONPT=CONPT0
      CONPT(4)="PREC+LAT M"
      CONPT(5)="SURFACE"   ; CONPT(8)="RIVERS"
      QCON=(/ F, F, F, T, T, F, F, T, T, F, F/)
      CALL SET_CON (QCON,CONPT,"LAK MASS","(KG/M^2)        ","(10**-9 KG/SM^2)",1d0,1d9,icon_LKM)
      QCON=(/ F, F, F, T, T, F, F, T, T, F, F/)
      CALL SET_CON (QCON,CONPT,"LAK ENRG","(10**3 J/M^2)   ", "(10**-3 W/M^2)  ",1d-3,1d3,icon_LKE)
!****
 910  FORMAT (A72)
 911  FORMAT (72A1)
      END SUBROUTINE init_LAKES



      function horzdist_2pts(i1,j1,i2,j2)
      Use constant, Only: radius
      Use geom, Only: lon2d,sinlat2d,coslat2d,axyp
      implicit none
      real*8 :: horzdist_2pts
      integer :: i1,j1,i2,j2
      real*8 :: x1,y1,z1, x2,y2,z2
      if(i1.eq.i2 .and. j1.eq.j2) then ! within same box
         horzdist_2pts = SQRT(AXYP(I1,J1))
      else                      ! Use great circle distance
         x1 = coslat2d(i1,j1)*cos(lon2d(i1,j1))
         y1 = coslat2d(i1,j1)*sin(lon2d(i1,j1))
         z1 = sinlat2d(i1,j1)
         x2 = coslat2d(i2,j2)*cos(lon2d(i2,j2))
         y2 = coslat2d(i2,j2)*sin(lon2d(i2,j2))
         z2 = sinlat2d(i2,j2)
         horzdist_2pts = radius*acos(x1*x2+y1*y2+z1*z2)
      endif
      end function horzdist_2pts



      integer function get_dir(I,J,ID,JD,IM,JM)
      Use domain_decomp_atm, Only: grid
      Use domain_decomp_1d, Only: hasSouthPole, hasNorthPole

!@sum get_dir derives the locally orientated river direction
      integer, intent(in) :: I,J,ID,JD,IM,JM
      integer ::  DI,DJ

      DI=I-ID
      IF (DI.eq.IM-1) DI=-1
      IF (DI.eq.1-IM) DI=1
      DJ=J-JD
      get_dir=-99
      if (DI.eq.-1 .and. DJ.eq.-1) then
         get_dir=1
      elseif (DI.eq.-1 .and. DJ.eq.0) then
         get_dir=8
      elseif (DI.eq.-1 .and. DJ.eq.1) then
         get_dir=7
      elseif (DI.eq.0 .and. DJ.eq.1) then
         get_dir=6
      elseif (DI.eq.0 .and. DJ.eq.0) then
         get_dir=0
      elseif (DI.eq.0 .and. DJ.eq.-1) then
         get_dir=2
      elseif (DI.eq.1 .and. DJ.eq.-1) then
         get_dir=3
      elseif (DI.eq.1 .and. DJ.eq.0) then
         get_dir=4
      elseif (DI.eq.1 .and. DJ.eq.1) then
         get_dir=5
      end if
      if (hasNorthPole(grid) .and. J.eq.JM) then         ! north pole
         if (DI.eq.0) then
            get_dir=6
         else
            get_dir=8
         end if
      elseif (hasSouthPole(grid) .and. J.eq.1) then      ! south pole
         if (DI.eq.0) then
            get_dir=2
         else
            get_dir=8
         end if
! these next two cases need two latitudes to be on same processor
      elseif (hasNorthPole(grid) .and. J.eq.JM-1) then
         if (JD.eq.JM) get_dir=2
      elseif (hasSouthPole(grid) .and. J.eq.2) then
         if (JD.eq.1) get_dir=6
      end if
      if (get_dir.eq.-99) then
         print*,"get_dir error",i,j,id,jd
         get_dir=0
      end if
      end function get_dir



      SUBROUTINE RIVERF
!@sum  RIVERF transports lake water from each grid box downstream
!@auth Gary Russell/Gavin Schmidt
!@ver  2011/11/03 (based on LB265)

      Use CONSTANT, Only: shw,rhow,teeny,bygrav,tf
      Use RESOLUTION, Only: im,jm
      Use MODEL_COM, Only: dtsrc,itime
      Use ATM_COM, Only: zatmo
      Use DOMAIN_DECOMP_ATM, Only: HALO_UPDATE, GRID,getDomainBounds
      Use domain_decomp_1d, Only: hasSouthPole, hasNorthPole

      Use GEOM, Only: axyp,byaxyp,imaxj
      Use DIAG_COM, Only: aij=>aij_loc,ij_ervr,ij_mrvr,ij_f0oc, &
                           jreg,j_rvrd,j_ervr,ij_fwoc,ij_ervro,ij_mrvro, ij_rvrflo,itlake,itlkice,itocean,itoice
      Use GHY_COM, Only: fearth
      Use FLUXES, Only: atmocn,focean,fland
      Use LAKES, Only: kdirec,rate,iflow,jflow,river_fac, kd911,ifl911,jfl911,lake_rise_max
      Use LAKES_COM, Only: tlake,gml,mwl,mldlk,flake,hlake,dlake,glake
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use TimerPackage_Mod, Only: StartTimer=>Start,StopTimer=>Stop

#ifdef SCM
      Use SCM_COM, Only: SCMopt,SCMin
#endif
#ifdef TRACERS_WATER
      Use TRDIAG_COM, Only: taijn =>taijn_loc , tij_rvr, tij_rvro
      Use LAKES_COM, Only: NTM,TRLAKE
#endif

      IMPLICIT NONE
      INTEGER :: FROM,J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0,I_1,I_0H,I_1H
      logical :: have_pole,have_south_pole,have_north_pole
      INTEGER :: ILOOP_MIN,ILOOP_MAX,JLOOP_MIN,JLOOP_MAX
      INTEGER I,J,IU,JU,ID,JD,JR,ITYPE,KD
      Real*8 MWLSILL, &  !  lake mass (kg) below sill depth
            MWlSILLD, &  !  downstream lake mass (kg) below sill depth
            MLM,DMM,DGM,HLK1,DPE,FLFAC, FLAKEU,FLAKED,BYOAREA
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: FLOW,EFLOW
!@var URATE upstream fractional rate of river flow per time step (only for special case)
      REAL*8 :: URATE = 1d-6  ! roughly 10 day e-folding time

#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: DTM
      REAL*8, DIMENSION(NTM,GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: TRFLOW
#endif

      REAL*8, DIMENSION(:,:), POINTER :: RSI,GTEMP,GTEMPR,FLOWO,EFLOWO
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:), POINTER :: GTRACER,TRFLOWO
#endif

!**** MWL (kg) = Lake water in cell, defined even when FLAKE = 0
!****            such as ice sheets, deserts, and partial ocean cells
!**** GML (J)  = Liquid lake enthalpy
!**** TLAKE(C) = Lake surface temperature
!****
!**** If FOCEAN > 0: FLAKE = 0
!****                (ID,JD) = (IU,JU), allowing MWL to exit to ocean
!**** If FOCEAN = 0: FLAKE is between 0 and 95% of (1-FLICE) of a cell
!****                KDIREC = 0: MWL may not exit except via backwash
!****                KDIREC = 1-8: counter-clockwise start (IU+1,JU+1)
!****                              MWL exits via normal downstream flow
!****                              MWL accepts backwash from KDIREC=0,9
!****                KDIREC = 9: MWL swashes water to adjacent KDIREC=9
!****                            MWL may backwash to upstream cells
!****                            Caspian, Aral, Great Salt, Chad
!****
!**** Check whether emergency directions are needed:
!**** If (KDIRECu==0 and FLAKEu > .949*(FLAKEu+FEARTHu) and
!****     MWLu > RHOW*FLAKEu*AXYPu*(HLAKEu+LAKE_RISE_MAX)):
!****    ID,JD,KD = ID911u,JD911u,KS911u
!****    MWLSILLu = RHOW*FLAKEu*AXYPu*(HLAKEu+LAKE_RISE_MAX)
!**** Otherwise: MWLSILLu = RHOW*FLAKEu*AXYPu*HLAKEu
!****
!**** Backwash River Flow (checked first), DMM (kg) < 0:
!**** Downstream water above upstream sill =
!****   MWLSILLd = RHOW*FLAKEd*AXYPd * [HLAKEd + (ZATMOu-ZATOMd)/GRAV]
!**** DMM = URATE*DTSRC * [FLAKEd*AXYPd*(MWLu-MWLSILLu) -
!****                    - FLAKEu*AXYPu*(MWLd-MWLSILLd)] /
!****                     (FLAKEu*AXYPu + FLAKEd*AXYPd)
!****
!**** Downstream Regular, Emergency and to Ocean Flow:
!**** DMM = (MWLu-MWLSILLu) * RATEu
!****
!**** Swash Water Back and Forth in Internal Sea, KDIRECu = KDIRECd = 9:
!**** Water above sill = MWLSILLu = RHOW*FLAKEu*AXYPu*HLAKEu
!****                    MWLSILLd = RHOW*FLAKEd*AXYPd*HLAKEd
!**** DMM = URATE*DTSRC * [FLAKEd*AXYPd*(MWLu-MWLSILLu) -
!****                    - FLAKEu*AXYPu*(MWLd-MWLSILLd)] /
!****                     (FLAKEu*AXYPu + FLAKEd*AXYPd)

      call startTimer('RIVERF()')
      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP =J_0S, J_STOP_SKP =J_1S, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      have_south_pole=hasSouthPole(grid)
      have_north_pole=hasNorthPole(grid)
      have_pole=have_south_pole .or. have_north_pole
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      FLOWO => ATMOCN%FLOWO
      EFLOWO => ATMOCN%EFLOWO

      FLOW = 0. ; EFLOW = 0.
      FLOWO = 0. ; EFLOWO = 0.

#ifdef TRACERS_WATER
      TRFLOWO => ATMOCN%TRFLOWO
      TRFLOW = 0.
      TRFLOWO = 0.
#endif

      RSI => LAKEICE%RSI
      GTEMP => ATMOCN%GTEMP
      GTEMPR => ATMOCN%GTEMPR
#ifdef TRACERS_WATER
      GTRACER => ATMOCN%GTRACER
#endif

      CALL HALO_UPDATE(GRID, FLAND) ! fixed
      CALL HALO_UPDATE(GRID,FEARTH)
      CALL HALO_UPDATE(GRID, HLAKE)
      CALL HALO_UPDATE(grid, FLAKE)
      CALL HALO_UPDATE(grid,   MWL)
      CALL HALO_UPDATE(grid, MLDLK)
      CALL HALO_UPDATE(grid, TLAKE)
      CALL HALO_UPDATE(grid,  RATE) ! fixed
#ifdef TRACERS_WATER
      CALL HALO_UPDATE(grid,  GTRACER, jdim=3)
      CALL HALO_UPDATE(grid,  TRLAKE(:,1,:,:), jdim=3)
#endif

!**** Calculate fluxes downstream if lake mass is above sill height HLAKE (m)
!**** Also allow flow into ocean fraction of same box if KDIREC=0
!**** SPECIAL CASE: If the downstream box has FLAKE=0.95 and KDIREC=0 (i.e.
!**** no outlet) then the only way to prevent excess water build up is
!**** to allow a back-flux. Take account of mean topography change as
!**** well. This is mainly an issue for the Caspian and Aral Seas.
!**** EMERGENCY CASE: If excess accumulation occurs anyway, Use emergency
!**** river direction (KD911) if level is lake_rise_max m above orig depth.
!**** Loop now includes polar boxes

! note on MPI fixes: since different PEs can influence the downstream
! accumulation of FLOW etc, we loop on the haloed variables to ensure
! that contributions from the halo are included in FLOW/FLOWO etc.
! If downstream box is outside the interior, cycle - this is dealt with on
! a separate PE

      If (HAVE_SOUTH_POLE)  Then  ;  JLOOP_MIN = 1
                            Else  ;  JLOOP_MIN = J_0H  ;  EndIf
      If (HAVE_NORTH_POLE)  Then  ;  JLOOP_MAX = JM
                            Else  ;  JLOOP_MAX = J_1H  ;  EndIf

      DO JU=JLOOP_MIN,JLOOP_MAX
        if(i_0.eq.i_0h) then ! no I halo - latlon grid
          iloop_min=1
          iloop_max=IMAXJ(JU)
        else                 !  Cube-Sphere grid
          iloop_min=i_0h
          iloop_max=i_1h
          if(ju.lt.1 .or. ju.gt.jm) then
! avoid nonexistent SW/NW/SE/NE halo corner of a cubed sphere face.
! instead, mark nonexistent cells with a code in the KDIREC array?
            iloop_min=max(iloop_min,1)
            iloop_max=min(iloop_max,im)
          endif
        endif

        DO IU=ILOOP_MIN,ILOOP_MAX
          If (KDIREC(IU,JU) == 9)  GoTo 400  !  large internal seas
!**** determine whether we have an emergency:
!**** i.e. no outlet, max extent, more than 100m above original height
          If (KDIREC(IU,JU) == 0 .and. &
              FLAKE(IU,JU) > .949d0*(FLAKE(IU,JU)+FEARTH(IU,JU)) .and. &
              MWL(IU,JU) > (HLAKE(IU,JU)+LAKE_RISE_MAX)*FLAKE(IU,JU)*RHOW*AXYP(IU,JU) .and. &
              KD911(IU,JU) > 0)  Then
!**** Use emergency directions
            KD=KD911(IU,JU)
            JD=JFL911(IU,JU)
            ID=IFL911(IU,JU)
            MWLSILL = RHOW*(HLAKE(IU,JU)+lake_rise_max)* FLAKE(IU,JU)*AXYP(IU,JU)
!**** No emergency, Use normal directions
          ELSE
            KD=KDIREC(IU,JU)
            JD=JFLOW(IU,JU)
            ID=IFLOW(IU,JU)
            MWLSILL = RHOW*HLAKE(IU,JU)*FLAKE(IU,JU)*AXYP(IU,JU)
          END IF

          If (KD == 0 .and. FLAND(IU,JU)*FOCEAN(IU,JU) == 0)  Cycle
! only calculate for downstream interior + halo boxes.
! this allows for calcs of halo -> interior & interior-> halo flow
          If (JD > J_1H .or. JD < J_0H .or. ID > I_1H .or. ID < I_0H)  Cycle

!****
!**** Apply possible backwash river flow
!****
          If ((KDIREC(ID,JD) >= 1 .and. KDIREC(ID,JD) <= 8) .or. FLAKE(ID,JD) <= .949d0*(FLAKE(ID,JD)+FEARTH(ID,JD)))  GoTo 200
          MWLSILLD = RHOW * AXYP(ID,JD) * FLAKE(ID,JD) * (HLAKE(ID,JD) + byGRAV*Max(ZATMO(IU,JU)-ZATMO(ID,JD),0d0))
          If (MWL(ID,JD) <= MWLSILLD)  GoTo 200
!**** Backwash river flow is invoked
          If (FLAKE(IU,JU) > 0)  Then
            DMM = URATE*DTSRC * (FLAKE(ID,JD)*AXYP(ID,JD)*(MWL(IU,JU)-MWLSILL) - &
                                 FLAKE(IU,JU)*AXYP(IU,JU)*(MWL(ID,JD)-MWLSILLD)) / &
                  (FLAKE(IU,JU)*AXYP(IU,JU) + FLAKE(ID,JD)*AXYP(ID,JD))
            If (DMM >= 0)  GoTo 200
          Else
            DMM = - (MWL(ID,JD)-MWLSILLD)*URATE*DTSRC  ;  EndIf

          DGM = TLAKE(ID,JD)*DMM*SHW  !  TLAKE always defined
#ifdef TRACERS_WATER
          DTM(:) = DMM*GTRACER(:,ID,JD)
#endif
          GoTo 300

!****
!**** Downstream regular, emergency and to ocean river flow
!****
  200     If (MWL(IU,JU) <= MWLSILL)  Cycle
          DMM = (MWL(IU,JU)-MWLSILL)*RATE(IU,JU)
          If (MWL(IU,JU)-DMM < 1d-6)  DMM = MWL(IU,JU)
          DMM = Min(DMM,.5*RHOW*AXYP(IU,JU)) ! minimise 'flood' events!
          If (FLAKE(IU,JU) > 0)  Then
            MLM = RHOW*MLDLK(IU,JU)*FLAKE(IU,JU)*AXYP(IU,JU)
            If (DMM > .95d0*MLM) Write (0,920) IU,JU,ID,JD, DMM/AXYP(IU,JU), MLM/AXYP(IU,JU), MWL(IU,JU)/AXYP(IU,JU), &
               MWLSILL/AXYP(IU,JU), MLDLK(IU,JU), HLAKE(IU,JU), FLAKE(IU,JU), FLAKE(IU,JU)+FEARTH(IU,JU)
  920       FORMAT ('RIVERF: DMM>.95*MLM     DMM       MLM       MWL     MWLSILL    MLDLK    HLAKE    FLAKE   FL+FE' / &
                    I6,I3.3,' >=',I3,I3.3,1X,6F10.3,2F8.4)
            DMM = Min (DMM,.95d0*MLM)  ;  EndIf
          DGM = TLAKE(IU,JU)*DMM*SHW  !  TLAKE always defined
#ifdef TRACERS_WATER
          If (FLAKE(IU,JU) > 0) &
            Then  ;  DTM(:) = DMM*GTRACER(:,IU,JU)
            Else  ;  DTM(:) = DMM*TRLAKE(:,1,IU,JU)/MWL(IU,JU)  ;  EndIf
#endif

!****
!**** Apply Backwash and Downstream flow to FLOW and EFLOW arrays
!****
  300     FLOW(IU,JU)  =  FLOW(IU,JU) - DMM
          EFLOW(IU,JU) = EFLOW(IU,JU) - DGM
#ifdef TRACERS_WATER
          TRFLOW(:,IU,JU) = TRFLOW(:,IU,JU) - DTM(:)
#endif

!**** diagnostics of outward flow (inward flow saved later)
          AIJ(IU,JU,IJ_MRVRO) = AIJ(IU,JU,IJ_MRVRO) + DMM
          AIJ(IU,JU,IJ_ERVRO) = AIJ(IU,JU,IJ_ERVRO) + DGM
#ifdef TRACERS_WATER
          TAIJN(IU,JU,TIJ_RVRO,:) = TAIJN(IU,JU,TIJ_RVRO,:) + DTM(:)*byAXYP(IU,JU)
#endif

!**** Calculate adjustments for poles
          FLFAC = 1  !  default = no pole adjustment
          If (HAVE_POLE)  Then
            If (JU==1 .or. JU==JM)  FLFAC = IM
            If (JD==1 .or. JD==JM)  FLFAC = 1d0/IM  ;  EndIf

! check to ensure that arrays outside the interior are not updated.
          If (JD < J_0 .or. JD > J_1 .or. ID < I_0 .or. ID > I_1)  Cycle

          If (FOCEAN(ID,JD) == 0)  Then
            DPE = 0  !  DMM*(ZATMO(IU,JU)-ZATMO(ID,JD))
            FLOW(ID,JD)  =  FLOW(ID,JD) +  DMM     *FLFAC
            EFLOW(ID,JD) = EFLOW(ID,JD) + (DGM+DPE)*FLFAC
#ifdef TRACERS_WATER
            TRFLOW(:,ID,JD)=TRFLOW(:,ID,JD) +DTM(:)*FLFAC
#endif

          Else ! Save river mouth flow to for output to oceans
!**** DPE: also add potential energy change to ocean.
!**** Normally ocean is at sea level (Duh!), but in some boxes ZATMO
!**** may not be zero if there is land as well, while in the Caspian,
!**** the ocean level is below zero.
!**** Note: this is diasabled until PE of precip is properly calculated
!**** in atmosphere as well. Otherwise, there is an energy imbalance.
            DPE = 0  !  DMM*(ZATMO(IU,JU)-MIN(0d0,ZATMO(ID,JD)))
!**** possibly adjust mass (not heat) to allow for balancing of sea level
            DMM = RIVER_FAC * DMM
            FLOWO(ID,JD) = FLOWO(ID,JD) +  DMM     *FLFAC
            EFLOWO(ID,JD)=EFLOWO(ID,JD) + (DGM+DPE)*FLFAC
#ifdef TRACERS_WATER
            DTM(:) = RIVER_FAC * DTM(:)
            TRFLOWO(:,ID,JD) = TRFLOWO(:,ID,JD) + DTM(:)*FLFAC
#endif
!**** accumulate river runoff diags (moved from ground)
            Call INC_AJ (ID,JD,ITOCEAN,J_RVRD, DMM*byAXYP(ID,JD)*(1-RSI(ID,JD)))
            Call INC_AJ (ID,JD,ITOCEAN,J_ERVR, (DGM+DPE)*byAXYP(ID,JD)*(1-RSI(ID,JD)))
            Call INC_AJ (ID,JD,ITOICE,J_RVRD, DMM*byAXYP(ID,JD)*RSI(ID,JD))
            Call INC_AJ (ID,JD,ITOICE,J_ERVR, (DGM+DPE)*byAXYP(ID,JD)*RSI(ID,JD))
            AIJ(ID,JD,IJ_F0OC) = AIJ(ID,JD,IJ_F0OC) + (DGM+DPE)*byAXYP(ID,JD)
            AIJ(ID,JD,IJ_FWOC) = AIJ(ID,JD,IJ_FWOC) + DMM*byAXYP(ID,JD)
          EndIf
          JR=JREG(ID,JD)
          Call INC_AREG (ID,JD,JR,J_RVRD,DMM*byAXYP(ID,JD))
          Call INC_AREG (ID,JD,JR,J_ERVR,(DGM+DPE)*byAXYP(ID,JD))
          AIJ(ID,JD,IJ_MRVR) = AIJ(ID,JD,IJ_MRVR) + DMM
          AIJ(ID,JD,IJ_ERVR) = AIJ(ID,JD,IJ_ERVR) + DGM+DPE
#ifdef TRACERS_WATER
          TAIJN(ID,JD,TIJ_RVR,:) = TAIJN(ID,JD,TIJ_RVR,:) + DTM(:)*byAXYP(ID,JD)
#endif
#ifdef TRACERS_OBIO_RIVERS
          AIJ(ID,JD,IJ_rvrflo) = AIJ(ID,JD,IJ_rvrflo) + DMM
#endif

          Cycle

!****
!**** KDIREC=9: Check river flow among 4 adjacent cells in same sea
!**** Coding does not work for cells on opposite sides of IDL
!**** Do not count transport twice inside same processor: KD=4 or 6
!****
  400     Do 440 KD=2,8,6
          If (KD==2) &
             Then  ;  If (IU < I_0 .or. IU > I_1 .or. JU > J_1)  Cycle
                      ID=IU  ;  JD=JU+1
             Else  ;  If (JU < J_0 .or. JU > J_1 .or. IU > I_1)  Cycle
                      ID=IU+1  ;  JD=JU  ;  EndIf

          If (KDIREC(ID,JD) /= 9)  Cycle
          If (FLAKE(IU,JU) + FLAKE(ID,JD) == 0)  Cycle
          FLAKEU   = Max (FLAKE(IU,JU), .01d0)
          FLAKED   = Max (FLAKE(ID,JD), .01d0)
          MWLSILL  = RHOW * HLAKE(IU,JU) * FLAKEU * AXYP(IU,JU)
          MWLSILLD = RHOW * HLAKE(ID,JD) * FLAKED * AXYP(ID,JD)
          DMM = URATE*DTSRC * (FLAKED*AXYP(ID,JD)*(MWL(IU,JU)-MWLSILL) - &
                               FLAKEU*AXYP(IU,JU)*(MWL(ID,JD)-MWLSILLD)) / &
                (FLAKEU*AXYP(IU,JU) + FLAKED*AXYP(ID,JD))
          If (DMM > 0)  GoTo 420

!**** DMM < 0: Move water from grid cell (ID,JD) to cell (IU,JU)
          If (MWL(ID,JD) <= 1*RHOW*FLAKE(ID,JD)*AXYP(ID,JD))  GoTo 440
          If (DMM < 1*RHOW*FLAKE(ID,JD)*AXYP(ID,JD) - MWL(ID,JD)) &
              DMM = 1*RHOW*FLAKE(ID,JD)*AXYP(ID,JD) - MWL(ID,JD)
          DGM = TLAKE(ID,JD)*DMM*SHW
          JR = JREG(IU,JU)
          Call INC_AREG (IU,JU,JR,J_RVRD,-DMM*byAXYP(IU,JU))
          Call INC_AREG (IU,JU,JR,J_ERVR,-DGM*byAXYP(IU,JU))
          AIJ(IU,JU,IJ_MRVR) = AIJ(IU,JU,IJ_MRVR) - DMM
          AIJ(IU,JU,IJ_ERVR) = AIJ(IU,JU,IJ_ERVR) - DGM
          AIJ(ID,JD,IJ_MRVRO)= AIJ(ID,JD,IJ_MRVRO)- DMM
          AIJ(ID,JD,IJ_ERVRO)= AIJ(ID,JD,IJ_ERVRO)- DGM
#ifdef TRACERS_WATER
          DTM(:) = DMM*GTRACER(:,ID,JD)
          TAIJN(IU,JU,TIJ_RVR,:) = TAIJN(IU,JU,TIJ_RVR,:) - DTM(:)*byAXYP(IU,JU)
          TAIJN(ID,JD,TIJ_RVRO,:)= TAIJN(ID,JD,TIJ_RVRO,:) - DTM(:)*byAXYP(ID,JD)
#endif
          GoTo 430

!**** DMM > 0: Move water from grid cell (IU,JU) to cell (ID,JD)
  420     If (MWL(IU,JU) <= 1*RHOW*FLAKE(IU,JU)*AXYP(IU,JU))  GoTo 440
          If (DMM > MWL(IU,JU) - 1*RHOW*FLAKE(IU,JU)*AXYP(IU,JU)) &
              DMM = MWL(IU,JU) - 1*RHOW*FLAKE(IU,JU)*AXYP(IU,JU)
          DGM = TLAKE(IU,JU)*DMM*SHW
          JR = JREG(ID,JD)
          Call INC_AREG (ID,JD,JR,J_RVRD,DMM*byAXYP(ID,JD))
          Call INC_AREG (ID,JD,JR,J_ERVR,DGM*byAXYP(ID,JD))
          AIJ(ID,JD,IJ_MRVR) = AIJ(ID,JD,IJ_MRVR) + DMM
          AIJ(ID,JD,IJ_ERVR) = AIJ(ID,JD,IJ_ERVR) + DGM
          AIJ(IU,JU,IJ_MRVRO)= AIJ(IU,JU,IJ_MRVRO)+ DMM
          AIJ(IU,JU,IJ_ERVRO)= AIJ(IU,JU,IJ_ERVRO)+ DGM
#ifdef TRACERS_WATER
          DTM(:) = DMM*GTRACER(:,IU,JU)
          TAIJN(ID,JD,TIJ_RVR,:) = TAIJN(ID,JD,TIJ_RVR,:) + DTM(:)*byAXYP(ID,JD)
          TAIJN(IU,JU,TIJ_RVRO,:)= TAIJN(IU,JU,TIJ_RVRO,:) + DTM(:)*byAXYP(IU,JU)
#endif

!**** Update transportimg river arrays
  430     FLOW(IU,JU)  =  FLOW(IU,JU) - DMM
          FLOW(ID,JD)  =  FLOW(ID,JD) + DMM
          EFLOW(IU,JU) = EFLOW(IU,JU) - DGM
          EFLOW(ID,JD) = EFLOW(ID,JD) + DGM
#ifdef TRACERS_WATER
          TRFLOW(:,IU,JU) = TRFLOW(:,IU,JU) - DTM(:)
          TRFLOW(:,ID,JD) = TRFLOW(:,ID,JD) + DTM(:)
#endif
  440     Continue

        EndDo  !  End of Do IU= loop
      EndDo    !  End of Do JU= loop

!****
!**** Apply net river flow to continental reservoirs
!****
      DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          IF(FLAND(I,J)+FLAKE(I,J).gt.0.) THEN
            MWL(I,J) = MWL(I,J) +  FLOW(I,J)
            GML(I,J) = GML(I,J) + EFLOW(I,J)
#ifdef TRACERS_WATER
            TRLAKE(:,1,I,J) = TRLAKE(:,1,I,J) + TRFLOW(:,I,J)
#endif

!**** remove pathologically small values
            IF (MWL(I,J).lt.1d-6) THEN
              MWL(I,J)=0.
              GML(I,J)=0.
#ifdef TRACERS_WATER
              TRLAKE(:,1:2,I,J) = 0.
#endif
            END IF
            IF (FLAKE(I,J).gt.0) THEN
              HLK1=(MLDLK(I,J)*RHOW)*TLAKE(I,J)*SHW
              MLDLK(I,J) = MLDLK(I,J)+FLOW(I,J) / (RHOW*FLAKE(I,J)*AXYP(I,J))
              TLAKE(I,J)=(HLK1*FLAKE(I,J)*AXYP(I,J)+EFLOW(I,J)) /(MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J)*SHW)
!**** accumulate some diagnostics
              CALL INC_AJ(I,J,ITLAKE,J_RVRD, FLOW(I,J)*BYAXYP(I,J)*(1.-RSI(I,J)))
              CALL INC_AJ(I,J,ITLAKE,J_ERVR,EFLOW(I,J)*BYAXYP(I,J)*(1.-RSI(I,J)))
              CALL INC_AJ(I,J,ITLKICE,J_RVRD, FLOW(I,J)*BYAXYP(I,J) *RSI(I,J))
              CALL INC_AJ(I,J,ITLKICE,J_ERVR,EFLOW(I,J)*BYAXYP(I,J) *RSI(I,J))
            ELSE
              TLAKE(I,J)=GML(I,J)/(SHW*MWL(I,J)+teeny)
!**** accounting fix to ensure river flow with no lakes is counted
              CALL INC_AJ(I,J,ITLAKE,J_RVRD, FLOW(I,J)*BYAXYP(I,J))
              CALL INC_AJ(I,J,ITLAKE,J_ERVR,EFLOW(I,J)*BYAXYP(I,J))
            END IF
          END IF
        END DO
      END DO

      CALL PRINTLK("RV")
!**** Set GTEMP array for lakes
      DO J=J_0, J_1
        DO I=I_0, I_1
          IF (FLAKE(I,J).gt.0) THEN
            GTEMP(I,J)=TLAKE(I,J)
            GTEMPR(I,J) =TLAKE(I,J)+TF
#ifdef SCM
            if (SCMopt%Tskin) then
              GTEMP(I,J) = SCMin%Tskin - TF
              GTEMPR(I,J) = SCMin%Tskin
            endif
#endif

#ifdef TRACERS_WATER
            GTRACER(:,I,J)=TRLAKE(:,1,I,J)/(MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J))
#endif
            atmocn%MLHC(I,J) = SHW*MLDLK(I,J)*RHOW
          END IF
        END DO
      END DO

      do j=j_0,j_1
      do i=i_0,imaxj(j)
        if(flake(i,j).gt.0.) then
          DLAKE(I,J)=MWL(I,J)/(RHOW*FLAKE(I,J)*AXYP(I,J))
          GLAKE(I,J)=GML(I,J)/(FLAKE(I,J)*AXYP(I,J))
        else
          DLAKE(I,J)=0.
          GLAKE(I,J)=0.
        endif
        if(focean(i,j).gt.0.) then
          byoarea = 1.d0/(axyp(i,j)*focean(i,j))
          flowo(i,j) = flowo(i,j)*byoarea
          eflowo(i,j) = eflowo(i,j)*byoarea
#ifdef TRACERS_WATER
          trflowo(:,i,j) = trflowo(:,i,j)*byoarea
#endif
!        else
!          flowo(i,j) = 0.
!          eflowo(i,j) = 0.
!#ifdef TRACERS_WATER
!          trflowo(:,i,j) = 0.
!#endif
        endif
      enddo
      enddo

      call stopTimer('RIVERF()')
      RETURN
!****
      END SUBROUTINE RIVERF



      SUBROUTINE diag_RIVER
!@sum  diag_RIVER prints out the river outflow for various rivers
!@sum  (now parallel)
!@auth Gavin Schmidt

      Use CONSTANT, Only: rhow,teeny,undef
      Use RESOLUTION, Only: im,jm
      Use MODEL_COM, Only: modelEclock
      Use MODEL_COM, Only: jyear0,amon0,jdate0,jhour0,amon,itime,dtsrc,idacc,itime0,nday, calendar
      Use TimeConstants_mod, Only: INT_MONTHS_PER_YEAR
      Use DOMAIN_DECOMP_ATM, Only: GRID,WRITE_PARALLEL, AM_I_ROOT, getDomainBounds, sumxpe
      Use GEOM, Only: byaxyp
      Use DIAG_COM, Only: aij=>aij_loc,ij_mrvr
#ifdef TRACERS_WATER
      Use OldTracer_mod, Only: trname, trw0, itime_tr0,tr_wd_type,nWATER
      Use TRACER_COM, Only: NTM,n_water
      Use TRDIAG_COM, Only: taijn=>taijn_loc
      Use TRDIAG_COM, Only: tij_rvr,to_per_mil,units_tij,scale_tij
#endif
      Use LAKES_COM, Only: irvrmth,jrvrmth,namervr,nrvr
      Use TimeInterval_mod
      Use Rational_mod

      IMPLICIT NONE
      REAL*8 RVROUT(NRVR), RVROUT_root(NRVR), scalervr, days
      INTEGER INM,I,N,J
      LOGICAL increment
#ifdef TRACERS_WATER
      REAL*8 TRVROUT(NRVR,NTM)
#endif
!@var out_line local variable to hold mixed-type output for parallel I/O
      character(len=300) :: out_line
      integer :: I_0, I_1, J_0, J_1
      integer :: year, hour, date
      type (Rational) :: secondsPerYear

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DAYS=(Itime-Itime0)/REAL(nday,kind=8)
      call modelEclock%get(year=year, hour=hour, date=date)
      WRITE(out_line,900) JYEAR0,AMON0,JDATE0,JHOUR0,YEAR,AMON,DATE,HOUR,ITIME,DAYS
      IF (AM_I_ROOT()) CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
!**** convert kg/(source time step) to km^3/mon
      secondsPerYear = calendar%getMaxDaysInYear() * calendar%getSecondsPerDay()
      SCALERVR = 1d-9*real(secondsPerYear) / (INT_MONTHS_PER_YEAR*RHOW*DTSRC)

      RVROUT(:)=0
#ifdef TRACERS_WATER
      TRVROUT(:,:)=0.
#endif
!**** loop over whole grid
      DO J=J_0,J_1
        DO I=I_0,I_1
          DO INM=1,NRVR
            if (I.eq.IRVRMTH(INM).and. J.eq.JRVRMTH(INM)) THEN
              RVROUT(INM) = SCALERVR*AIJ(I,J,IJ_MRVR)/IDACC(1)
#ifdef TRACERS_WATER
              IF (RVROUT(INM).gt.0)  THEN
                DO N=1,NTM
                  if (to_per_mil(n).gt.0) then
                    if (TAIJN(I,J,TIJ_RVR,N_water).gt.0) then
                      TRVROUT(INM,N) = 1d3*(TAIJN(I,J,TIJ_RVR,N) / (trw0(n)*TAIJN(I,J,TIJ_RVR,N_water))-1.)
                    else
                      TRVROUT(INM,N)=undef
                    endif
                  else
                    TRVROUT(INM,N) = scale_tij (TIJ_RVR,n)*TAIJN(I,J,TIJ_RVR,N)/(AIJ(I,J,IJ_MRVR)*BYAXYP(I,J) +teeny)
                  end if
                END DO
              ELSE
                TRVROUT(INM,:)=undef
              END IF
#endif
            end if
          END DO
        END DO
      END DO

!**** gather diags + print out on root processor
      rvrout_root=0.
      call sumxpe(rvrout, rvrout_root, increment=.true.)

      IF (AM_I_ROOT()) THEN
        DO INM=1,NRVR,6
          WRITE(out_line,901) (NAMERVR(I-1+INM),RVROUT_root(I-1+INM),I=1,MIN(6,NRVR+1-INM))
          CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
        END DO
      END IF

#ifdef TRACERS_WATER
      DO N=1,NTM
        if (itime.ge.itime_tr0(n) .and. tr_wd_TYPE(n).eq.nWater) then
          rvrout_root=0.
          call sumxpe(trvrout(:,N), rvrout_root, increment=.true.)

          IF (AM_I_ROOT()) THEN
            WRITE(out_line,*) "River outflow tracer concentration ",trim(units_tij(tij_rvr,n)),":",TRNAME(N)
            CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
            DO INM=1,NRVR,6
              WRITE(out_line,901) (NAMERVR(I-1+INM),RVROUT_root(I-1+INM),I=1,MIN(6,NRVR+1-INM))
              CALL WRITE_PARALLEL(trim(out_line), UNIT=6)
            END DO
          END IF
        end if
      END DO
#endif

      RETURN
!****
 900  FORMAT ('1* River Outflow (km^3/mon) **  From:',I6,A6,I2,',  Hr',I3,6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X &
              ,'Dif:',F7.2,' Days')
 901  FORMAT (' ',A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3,5X, A8,':',F8.3,5X,A8,':',F8.3,5X,A8,':',F8.3)
      END SUBROUTINE diag_RIVER



      SUBROUTINE CHECKL (SUBR)
!@sum  CHECKL checks whether the lake variables are reasonable.
!@auth Gavin Schmidt/Gary Russell
      Use CONSTANT, Only: rhow
      Use RESOLUTION, Only: im,jm
      Use MODEL_COM, Only: qcheck
      Use FLUXES, Only: focean
      Use DOMAIN_DECOMP_ATM, Only: getDomainBounds, GRID
      Use GEOM, Only: axyp,imaxj
#ifdef TRACERS_WATER
      Use OldTracer_mod, Only: trname, t_qlimit
      Use TRACER_COM, Only: NTM
#endif
      Use LAKES
      Use LAKES_COM
      IMPLICIT NONE
      INTEGER :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0,I_1,I_0H,I_1H,njpol
      INTEGER I,J,N !@var I,J loop variables
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKL
#ifdef TRACERS_WATER
      integer :: imax,jmax
      real*8 relerr,errmax
#endif
      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_HALO=J_0H,J_STOP_HALO=J_1H, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      njpol = grid%J_STRT_SKP-grid%J_STRT

!**** Check for NaN/INF in lake data
      CALL CHECK3B (MWL(I_0:I_1,J_0:J_1)  ,I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'mwl')
      CALL CHECK3B (GML(I_0:I_1,J_0:J_1)  ,I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'gml')
      CALL CHECK3B (MLDLK(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'mld')
      CALL CHECK3B (TLAKE(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'tlk')

      QCHECKL = .FALSE.
      DO J=J_0S,J_1S  ;  DO I=I_0, I_1
         IF(FOCEAN(I,J).eq.0.) THEN
!**** check for negative mass
            IF (MWL(I,J).lt.0 .or. MLDLK(I,J).lt.0) THEN
               WRITE(6,*) 'After ',SUBR,': I,J,TSL,MWL,GML,MLD=', I,J,TLAKE(I,J),MWL(I,J),GML(I,J),MLDLK(I,J)
               QCHECKL = .TRUE.
            END IF
!**** check for reasonable lake surface temps
            IF (TLAKE(I,J).ge.50 .or. TLAKE(I,J).lt.-0.5) THEN
               WRITE(6,*) 'After ',SUBR,': I,J,TSL=',I,J,TLAKE(I,J)
               if (TLAKE(I,J).lt.-5.and.FLAKE(I,J).gt.0) QCHECKL = .TRUE.
            END IF
         END IF
!**** Check total lake mass ( <0.4 m, >20x orig depth)
         IF(FLAKE(I,J).gt.0.) THEN
!!          IF(MWL(I,J).lt.0.4d0*RHOW*AXYP(I,J)*FLAKE(I,J)) THEN
!!             WRITE (6,*) 'After ',SUBR,
!!    *           ': I,J,FLAKE,HLAKE,lake level low=',I,J,FLAKE(I,J),
!!    *           HLAKE(I,J),MWL(I,J)/(RHOW*AXYP(I,J)*FLAKE(I,J))
!!          END IF
!           IF(MWL(I,J).gt.RHOW*MAX(20.*HLAKE(I,J),3d1)*AXYP(I,J)*FLAKE(I,J)
            IF(MWL(I,J).gt.RHOW*(HLAKE(I,J)+lake_rise_max)*AXYP(I,J)*FLAKE(I,J))  THEN
               WRITE (6,*) 'After ',SUBR,': I,J,FLAKE,HLAKE,lake level high=',I,J,FLAKE(I,J), &
                           HLAKE(I,J),MWL(I,J)/(RHOW*AXYP(I,J)*FLAKE(I,J))
            END IF
         END IF
      END DO  ;  END DO

#ifdef TRACERS_WATER
      do n=1,ntm
!**** Check for neg tracers in lake
         if (t_qlimit(n)) then
            do j=J_0,J_1  ;  do i=I_0,imaxj(j)
               if (focean(i,j).eq.0) then
                  if (trlake(n,1,i,j).lt.0 .or. trlake(n,2,i,j).lt.0) then
                     print*,"Neg tracer in lake after ",SUBR,i,j,trname(n),trlake(n,:,i,j)
                     QCHECKL=.TRUE.
                  end if
               end if
            end do  ;  end do
         end if
!**** Check conservation of water tracers in lake
         if (trname(n).eq.'Water') then
            errmax = 0. ; imax=I_0 ; jmax=J_0
            do j=J_0,J_1  ;  do i=I_0,imaxj(j)
               if (focean(i,j).eq.0) then
                  if (flake(i,j).gt.0) then
                     relerr = max (abs(trlake(n,1,i,j)-mldlk(i,j)*rhow*flake(i,j)*axyp(i,j))/trlake(n,1,i,j), &
                                   abs(trlake(n,1,i,j)+trlake(n,2,i,j)-mwl(i,j))/(trlake(n,1,i,j)+trlake(n,2,i,j)))
                  else
                     if ((mwl(i,j).eq.0 .and. trlake(n,1,i,j)+trlake(n,2,i,j) .gt.0) .or. &
                         (mwl(i,j).gt.0 .and. trlake(n,1,i,j)+trlake(n,2,i,j).eq.0))  then
                        print*,"CHECKL ",SUBR,i,j,mwl(i,j),trlake(n,1:2,i,j)
                        relerr=0.
                     else
                        if (mwl(i,j).gt.1d-20) then
                           relerr = abs (trlake(n,1,i,j)+trlake(n,2,i,j)-mwl(i,j)) / (trlake(n,1,i,j)+trlake(n,2,i,j))
                        else
                           if (mwl(i,j).gt.0) print*,"CHECKL2 ",SUBR,i,j,mwl(i,j),trlake(n,1:2,i,j)
                           relerr=0.
                        end if
                     end if
                  end if
                  if (relerr.gt.errmax) then
                     imax=i ; jmax=j ; errmax=relerr
                  end if
               end if
            end do  ;  end do
          print*,"Relative error in lake mass after ",trim(subr),":",imax,jmax,errmax,trlake(n,:,imax,jmax), &
                 mldlk(imax,jmax)*rhow*flake(imax,jmax)*axyp(imax,jmax), &
                 mwl(imax,jmax)-mldlk(imax,jmax)*rhow*flake(imax,jmax)*axyp(imax,jmax)
         end if
      end do
#endif

      IF (QCHECKL) call stop_model ('CHECKL: Lake variables out of bounds',255)
      END SUBROUTINE CHECKL



      Subroutine DAILY_LAKE
!@sum  DAILY_LAKE does lake things at the beginning of every day
!@auth Gavin Schmidt, Gary L. Russell
!@ver  2019/02/08
      Use CONSTANT,    Only: rhow,by3,pi,lhm,shi,shw,teeny,tf
      Use RESOLUTION,  Only: im
      Use LAKES,       Only: minmld, variable_lk, hlake_min, lake_ice_max, dZSG
      Use LAKES_COM,   Only: mwl,flake,mldlk,tlake,gml,svflake,hlake,dlake,glake
      Use SEAICE_COM,  Only: lakeice=>si_atm
      Use SEAICE,      Only: ace1i,xsi,ac2oim
      Use GEOM,        Only: axyp,imaxj,byaxyp
      Use GHY_COM,     Only: fearth
      Use FLUXES,      Only: atmocn,dmwldf,dgml,fland,flice,focean
      Use LANDICE_COM, Only: mdwnimp,edwnimp
      Use DIAG_COM,    Only: j_run,j_erun,jreg,j_implm,J_IMPLH, AIJ=>AIJ_LOC,itlkice,itlake, &
                             IJ_MLKtoGR,IJ_HLKtoGR,IJ_IMPMKI,IJ_IMPHKI
      Use DOMAIN_DECOMP_ATM, Only: getDomainBounds, GRID
#ifdef IRRIGATION_ON
      Use IRRIGMOD, Only: read_irrig
#endif  /* IRRIGATION_ON   */
#ifdef TRACERS_WATER
      Use LAKES_COM,   Only:: trlake,ntm
      Use FLUXES,      Only:: dtrl
      Use LANDICE_COM, Only:: trdwnimp
#endif
#ifdef SCM
      Use SCM_COM, Only: SCMopt,SCMin
#endif
      IMPLICIT NONE

      Integer :: i,j, J_0,J_1,I_0,I_1, jr,itm
      Real*8  :: DLAKEnew, FEARTHold,FLAKEnew,FLAKEold,FLpFE, FMSI2,FMSI3,FMSI4, FHSI2,FHSI3,FHSI4, FRAC,FRACI,FRSAT, F_ENTR, &
                 HLK,HMLT, IMLT, MLDnew,MSInew,MWSAT, M1,M2,M1T1,M2T2, PLAKE,PLKIC, RSIold, SNOWnew,SUMH, TLAKEnew
      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI,GTEMP,GTEMPR
      REAL*8, DIMENSION(:,:,:), POINTER :: HSI
#ifdef TRACERS_WATER
      Real*8 :: hlk2,ftsi2(ntm),ftsi3(ntm),ftsi4(ntm),sumt,dtr(ntm) ,tottr(ntm)
      REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI
      REAL*8, DIMENSION(:,:,:), POINTER :: GTRACER
#endif

      RSI => LAKEICE%RSI
      MSI => LAKEICE%MSI
      HSI => LAKEICE%HSI
      SNOWI  => LAKEICE%SNOWI
      GTEMP  => ATMOCN%GTEMP
      GTEMPR => ATMOCN%GTEMPR
#ifdef TRACERS_WATER
      TRSI => LAKEICE%TRSI
      GTRACER => ATMOCN%GTRACER
#endif

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

#ifdef IRRIGATION_ON
!**** Read potential irrigation daily
      call read_irrig(.true.)
#endif  /* IRRIGATION_ON   */

!**** Set fixed topographic variables when FOCEAN == 0
!**** Solid ground topography is a linear function of cell area or FLAKE
!**** Lake fills like a triangle from lowest solid ground topography
!**** Vertical cross section of circular lake is then a concave parabola
!**** dZSG  = fixed linear range of continental solid ground topography (m) in grid cell
!**** ZLBOT = fixed lowest altitude of solid ground topography (m) in grid cell = ZSOLID - .5*dZSG
!**** ZLTOP = ZLBOT + FLAKE * dZSG = lake top altitude (m); ZLTOP is a linear function of FLAKE or cell area

!**** If DLAKEn < HLAKE_MIN, then DLAKEn = HLAKE_MIN and FLAKEn = MWL / RHOW * AXYP * DLAKEn
!**** If FLAKEn > .953125*(FLAKE+FEARTH), then FLAKEn = .953125*FLpFE and DLAKEn = MWLo/RHOW*AXYP*FLAKEn
!**** If FLAKEn == FLAKEo, then lake neither expands nor contracts 

      SVFLAKE=FLAKE  ! save for ghy purposes
      If (variable_lk .ne. 0) then

      DO J=J_0,J_1  ;  DO I=I_0,IMAXJ(J)
         JR=JREG(I,J)
!        FLpFE = FLAKE(I,J) + FEARTH(I,J)
         FLpFE = 1 - FOCEAN(I,J) - FLICE(I,J)
         IF (FLpFE == 0 .or. FOCEAN(I,J) > 0)  Cycle
!**** Save original fractions
         FLAKEold  = FLAKE(I,J)
         FEARTHold = FEARTH(I,J)
         RSIold    = RSI(I,J)
         PLAKE     = FLAKE(I,J)*(1-RSI(I,J))
         PLKIC     = FLAKE(I,J)*   RSI(I,J)
         DGML(I,J) = 0

!**** DLAKEnew = [ZLTOPn - ZLBOT] / 2 = FLAKEnew * dZSG / 2 = average lake thickness (m) because of triangular shape of ZSG
!**** FLAKEnew = 2 * DLAKEnew / dZSG = horizontal lake fraction of cell area
!**** MWLold   = RHOW * AXYP * DLAKEnew * FLAKEnew = RHOW * AXYP * DLAKEnew * 2 * DLAKEnew / dZSG = liquid lake mass (kg)
!**** DLAKEnew = Sqrt [MWLold * dZSG / 2 * RHOW * AXYP]
         DLAKEnew = Sqrt (MWL(I,J) * dZSG(I,J) / (2 * RHOW * AXYP(I,J)))
!**** Invoke limitations on DLAKEn and FLAKEn
         If (DLAKEnew < HLAKE_MIN)  DLAKEnew = HLAKE_MIN
         FLAKEnew = MWL(I,J) / (RHOW * AXYP(I,J) * DLAKEnew)
         If (FLAKEnew > Min(.953125*FLpFE, FLAKE(I,J)+.046875*FEARTH(I,J)))  Then
             FLAKEnew = Min(.953125*FLpFE, FLAKE(I,J)+.046875*FEARTH(I,J))
             DLAKEnew = MWL(I,J) / (RHOW * AXYP(I,J) * FLAKEnew)  ;  EndIf
         If (FLAKEnew < 1/1024d0)  &  !Then
             FLAKEnew = 0

!****
!**** FLAKEnew == FLAKE == 0
!****
         If (FLAKEnew == 0 .and. FLAKE(I,J)== 0)  Cycle  !  end Do I, Do J

!****
!**** FLAKEnew == FLAKE > 0
!****
         If (FLAKEnew == FLAKE(I,J) .and. FLAKE(I,J) > 0)  Cycle  !  end Do I, Do J

!****
!**** FLAKEnew == 0 and FLAKE > 0
!****
         If (FLAKEnew == 0 .and. FLAKE(I,J) > 0)  Then
!**** remove/do not create lakes that are too small
!**** transfer lake ice mass/energy for accounting purposes; do not add ice mass to river - instead Use implicit array
            IMLT = ACE1I + MSI(I,J) + SNOWI(I,J)
            HMLT = Sum(HSI(:,I,J))
            MDWNIMP(I,J) = MDWNIMP(I,J) + PLKIC*IMLT*AXYP(I,J)
            EDWNIMP(I,J) = EDWNIMP(I,J) + PLKIC*HMLT*AXYP(I,J)
!           MWL(I,J)     = MWL(I,J) + PLKIC*IMLT*AXYP(I,J)
!           GML(I,J)     = GML(I,J) + PLKIC*HMLT*AXYP(I,J)
            RSI(I,J)     = 0
            SNOWI(I,J)   = 0
            HSI(1:2,I,J) = - LHM*XSI(1:2)*ACE1I
            HSI(3:4,I,J) = - LHM*XSI(3:4)*AC2OIM
            MSI(I,J)     = AC2OIM
            TLAKE(I,J)   = GML(I,J) / (SHW*MWL(I,J)+teeny)
            GTEMPR(I,J)  = TF
            MLDLK(I,J)   = MINMLD
            FLAKE(I,J)   = 0
            FEARTH(I,J)  = FLpFE
            FLAND(I,J)   = 1
!**** Accumulate diagnostics
            AIJ(I,J,IJ_IMPMKI) = AIJ(I,J,IJ_IMPMKI) + PLKIC*IMLT
            AIJ(I,J,IJ_IMPHKI) = AIJ(I,J,IJ_IMPHKI) + PLKIC*HMLT
            Call INC_AJ (I,J,ITLKICE,J_IMPLM,PLKIC*IMLT)
            Call INC_AJ (I,J,ITLKICE,J_IMPLH,PLKIC*HMLT)
!           Call INC_AJ (I,J,ITLKICE,J_IMELT,PLKIC*IMLT)
!           Call INC_AJ (I,J,ITLKICE,J_HMELT,PLKIC*HMLT)
!           Call INC_AREG (I,J,JR,J_IMELT,PLKIC*IMLT)
!           Call INC_AREG (I,J,JR,J_HMELT,PLKIC*HMLT)
            Call INC_AREG (I,J,JR,J_IMPLM,PLKIC*IMLT)
            Call INC_AREG (I,J,JR,J_IMPLH,PLKIC*HMLT)
#ifdef TRACERS_WATER
            DO ITM=1,NTM
               TRLAKE(ITM,1,I,J) = SUM(TRLAKE(ITM,:,I,J))
               TRLAKE(ITM,2,I,J) = 0
!              TRLAKE(ITM,1,I,J) = TRLAKE(ITM,1,I,J) + TRLAKE(ITM,2,I,J) + RSI(I,J)*FLAKE(I,J)*SUM(TRSI(ITM,:,I,J))*AXYP(I,J)
               TRDWNIMP(ITM,I,J) = TRDWNIMP(ITM,I,J) + SUM(TRSI(ITM,:,I,J))*PLKIC*AXYP(I,J)
               TRSI(ITM,:,I,J) = 0
            END DO
#endif
#ifdef SCM
            if (SCMopt%Tskin)  GTEMPR(I,J) = SCMin%Tskin
#endif
            GoTo 600  ;  EndIf  !  If [FLAKEnew == 0 .and. FLAKE(I,J) > 0]

!****
!**** FLAKEnew > FLAKE
!****
         If (FLAKEnew > FLAKE(I,J))  Then
!**** If FLAKEnew > FLAKEold, then lake expands and new ground beneath lake must be saturated, DLAKEnew does not change
!**** DMWLDF = water deficit over FEARTH fraction (kg/m^2) = saturated ground water - actual ground water
!**** Recompute FLAKEnew based on old values and DLAKEnew which is not modified
!**** MWLold = FLAKEnew * AXYP * RHOW * DLAKEnew + (FLAKEnew - FLAKEold) * AXYP * DMWLDF
!**** MWLold + FLAKEold * AXYP * DMWLDF = FLAKEnew * AXYP * (RHOW * DLAKEnew + DMWLDF)
!**** FLAKEnew = (MWLold + FLAKEold * AXYP * DMWLDF) / AXYP * (RHOW * DLAKEnew + DMWLDF)
!**** MWSAT  = (FLAKEnew - FLAKEold) * AXYP * DMWLDF (kg) = water used to saturate ground beneath incresed FLAKE
!**** DGML   = GML * MWSAT/MWL (J)                       = energy used to saturate ground beneath incresed FLAKE 
!**** MWLnew = MWLold - MWSAT
!**** RSInew = RSIold * FLAKEold / FLAKEnew = PLKIC / FLAKEnew
            If (DMWLDF(I,J) > 0)  Then
               FLAKEnew  = (MWL(I,J) + FLAKE(I,J) *AXYP(I,J)*DMWLDF(I,J)) / (AXYP(I,J) * (RHOW*DLAKEnew + DMWLDF(I,J)))
               MWSAT     = (FLAKEnew - FLAKE(I,J))*AXYP(I,J)*DMWLDF(I,J)
               FRSAT     = MWSAT / MWL(I,J)
               If (MWL(I,J) <= MWSAT) Then
                  Write (0,*) 'DAILY_LAKE: MWLnew is less than water needed to saturate ground beneath new FLAKE'
                  Write (6,*) 'DAILY_LAKE: MWLnew is less than water needed to saturate ground beneath new FLAKE'
                  Call STOP_MODEL ('DAILY_LAKE: FLAKEnew > FLAKE and not enough water for saturation', 255)  ;  EndIf
               MWL(I,J)  = MWL(I,J) - MWSAT
               DGML(I,J) = GML(I,J) * FRSAT
               GML(I,J)  = GML(I,J) - DGML(I,J)
               MLDLK(I,J)= MLDLK(I,J) * (1-FRSAT)
#ifdef TRACERS_WATER
               DTRL(:,I,J)    = (TRLAKE(:,1,I,J)+TRLAKE(:,2,I,J)) * FRSAT
               TRLAKE(:,:,I,J) = TRLAKE(:,:,I,J) * (1-FRSAT)
#endif
               Call INC_AJ (I,J,ITLAKE, J_RUN ,PLAKE*DMWLDF(I,J)*(FLAKEnew-FLAKE(I,J)))
               Call INC_AJ (I,J,ITLKICE,J_RUN ,PLKIC*DMWLDF(I,J)*(FLAKEnew-FLAKE(I,J)))
               Call INC_AJ (I,J,ITLAKE, J_ERUN,PLAKE*DGML(I,J)*byAXYP(I,J))
               Call INC_AJ (I,J,ITLKICE,J_ERUN,PLKIC*DGML(I,J)*byAXYP(I,J))
               AIJ(I,J,IJ_MLKtoGR) = AIJ(I,J,IJ_MLKtoGR) + DMWLDF(I,J)*(FLAKEnew-FLAKE(I,J))
               AIJ(I,J,IJ_HLKtoGR) = AIJ(I,J,IJ_HLKtoGR) + DGML(I,J)*byAXYP(I,J)  ;  EndIf
            RSI(I,J) = PLKIC / FLAKEnew
            GoTo 500  ;  EndIf  !  FLAKEnew > FLAKE

!****
!**** FLAKEnew < FLAKE
!****
         If (FLAKEnew < FLAKE(I,J))  Then
!**** Conserve lake ice
            If (PLKIC > FLAKEnew)  Then  !  crunch ice up ?
!****          Crunch ice up, RSInew = 1
               FRAC = PLKIC / FLAKEnew
               SNOWnew = SNOWI(I,J)*FRAC
               MSInew  = (MSI(I,J)+ACE1I)*FRAC - ACE1I
               RSI(I,J) = 1
               FMSI3 = ACE1I*(FRAC-1)  !  kg/m^2 flux over new fraction
               FMSI2 = FMSI3*XSI(1)
               FMSI4 = FMSI3*XSI(4)
               FHSI2 = FMSI2*HSI(1,I,J) / (XSI(1)*(ACE1I+SNOWI(I,J)))
               If (FMSI3 < FRAC*XSI(2)*(ACE1I+SNOWI(I,J))) &
                  Then  ;  FHSI3 = FMSI3*HSI(2,I,J) / (XSI(2)*(ACE1I+SNOWI(I,J)))
                  Else  ;  FHSI3 = HSI(2,I,J)*FRAC + & 
                                   (FMSI3-FRAC*XSI(2)*(ACE1I+SNOWI(I,J)))*HSI(1,I,J) / (XSI(1)*(ACE1I+SNOWI(I,J)))  ;  EndIf
               If (FMSI4 < FRAC*XSI(3)*MSI(I,J)) &
                  Then  ;  FHSI4 = FMSI4*HSI(3,I,J) / (XSI(3)*MSI(I,J))
                  Else  ;  FHSI4 = HSI(3,I,J)*FRAC + (FMSI4-FRAC*XSI(3)*MSI(I,J))*FHSI3 / FMSI3  ;  EndIf
               HSI(1,I,J) = HSI(1,I,J)*(ACE1I+SNOWNEW) / (ACE1I+SNOWI(I,J))
               HSI(2,I,J) = HSI(2,I,J)*FRAC + FHSI2 - FHSI3
               HSI(3,I,J) = HSI(3,I,J)*FRAC + FHSI3 - FHSI4
               HSI(4,I,J) = HSI(4,I,J)*FRAC + FHSI4
!**** all tracers --> tracer*FRAC, then adjust layering
#ifdef TRACERS_WATER
               FTSI2(:) = FMSI2*TRSI(:,1,I,J) / (XSI(1)*(ACE1I+SNOWI(I,J)))
               IF (FMSI3.LT.FRAC*XSI(2)*(ACE1I+SNOWI(I,J))) THEN
                  Then  ;  FTSI3(:) = FMSI3*TRSI(:,2,I,J) / (XSI(2)*(ACE1I+SNOWI(I,J)))
                  Else  ;  FTSI3(:) = TRSI(:,2,I,J)*FRAC + &
                                      (FMSI3-FRAC*XSI(2)*(ACE1I+SNOWI(I,J)))*TRSI(:,1,I,J) / (XSI(1)*(ACE1I+SNOWI(I,J)))  ;  EndIf
               IF (FMSI4.LT.FRAC*XSI(3)*MSI(I,J)) THEN
                  Then  ;  FTSI4(:) = FMSI4*TRSI(:,3,I,J) / (XSI(3)*MSI(I,J))
                  Else  ;  FTSI4(:) = TRSI(:,3,I,J)*FRAC + (FMSI4-FRAC*XSI(3)*MSI(I,J))*FTSI3(:) / FMSI3  ;  EndIf
               TRSI(:,1,I,J) = TRSI(:,1,I,J)*(ACE1I+SNOWNEW) / (ACE1I+SNOWI(I,J))
               TRSI(:,2,I,J) = TRSI(:,2,I,J)*FRAC + FTSI2(:) - FTSI3(:)
               TRSI(:,3,I,J) = TRSI(:,3,I,J)*FRAC + FTSI3(:) - FTSI4(:)
               TRSI(:,4,I,J) = TRSI(:,4,I,J)*FRAC + FTSI4(:)
#endif
               MSI(I,J)   = MSInew
               SNOWI(I,J) = SNOWnew
            Else
!****          Do not crunch ice up, instead increase RSI
               RSI(I,J) = PLKIC / FLAKEnew
            EndIf  !  crunch ice up ?
         EndIf  !  If [FLAKEnew < FLAKE]

!****
!**** 0 >< FLAKEnew >< FLAKE
!****
  500    HLK = MWL(I,J) / (RHOW * FLAKEnew * AXYP(I,J))  !  potential new water height
         MLDnew = Min (HLK, Max(MINMLD,HLK-HLAKE(I,J)))
         If (MLDLK(I,J)*FLAKE(I,J) < FLAKEnew*MLDnew) Then
            If (FLAKE(I,J) == 0 .or. HLK <= MLDnew) &
!****          new or shallow lake
               Then  ;  MLDLK(I,J) = MLDnew
#ifdef TRACERS_WATER
                        TOTTR(:) = TRLAKE(:,1,I,J) + TRLAKE(:,2,I,J)
                        TRLAKE(:,2,I,J) = TOTTR(:)*(HLK-MLDLK(I,J)) / HLK
                        TRLAKE(:,1,I,J) = TOTTR(:)*     MLDLK(I,J)  / HLK
#endif
!****          transfer of mass from layer 2 to layer 1. adjust layer-1 properties
               Else  ;  f_entr = (FLAKEnew*MLDnew - MLDLK(I,J)*FLAKE(I,J)) / (FLAKEnew*HLK - MLDLK(I,J)*FLAKE(I,J))
                        M1 = MLDLK(I,J)*RHOW
                        M2 = MAX(MWL(I,J)/(FLAKE(I,J)*AXYP(I,J))-M1,0d0)
                        M1T1 = M1*TLAKE(I,J)
                        M2T2 = GML(I,J)/(SHW*FLAKE(I,J)*AXYP(I,J))-M1T1
                        TLAKEnew = (M1T1+f_entr*M2T2)/(M1+f_entr*M2)
                        TLAKE(I,J) = TLAKEnew
                        MLDLK(I,J) = MLDnew
#ifdef TRACERS_WATER
                        DTR(:)=TRLAKE(:,2,I,J)*f_entr
                        TRLAKE(:,1,I,J)=TRLAKE(:,1,I,J)+DTR(:)
                        TRLAKE(:,2,I,J)=TRLAKE(:,2,I,J)-DTR(:)
#endif
               EndIf  !  If [FLAKE(I,J) == 0 .or. HLK <= MLDnew]
         Else
            MLDLK(I,J) = MLDLK(I,J)*FLAKE(I,J) / FLAKEnew
         EndIf  !  If [MLDLK(I,J)*FLAKE(I,J) < FLAKEnew*MLDnew]
!****    adjust land surface fractions
         FLAKE(I,J)  = FLAKEnew
         FEARTH(I,J) = FLpFE - FLAKEnew
         FLAND(I,J)  = 1 - FLAKE(I,J)

!**** Adjust some radiative fluxes for conservation and restartability
!**** Complications due to ice or water going to earth if lake shrinks
  600    If (FLAKE(I,J) > FLAKEold)  Call RESET_SURF_FLUXES (I,J,4,1, FLAKEold, FLAKE(I,J))
         If (FLAKEold > FLAKE(I,J)) Then ! lake shrinks
! originally some open water
            If (PLAKE > 0)  Call RESET_SURF_FLUXES (I,J,1,4, FEARTHold, FEARTHold+PLAKE-FLAKE(I,J)*(1-RSI(I,J)))
! originally some ice, now different
            If (PLKIC > 0 .and. PLKIC /= FLAKE(I,J)*RSI(I,J)) &  !  PLKICold /= PLKICnew
               Call RESET_SURF_FLUXES (I,J,2,4, FEARTHold+PLAKE-FLAKE(I,J)*(1-RSI(I,J)), FEARTH(I,J))
         EndIf

      EndDo  ;  EndDo  !  Do I, Do J
      EndIf  !  If [VARIABLE_LK /= 0]

!****
      Call PRINTLK("DY")

!**** Set GTEMP array for lakes
      DO J=J_0,J_1  ;  DO I=I_0,IMAXJ(J)
         IF (FLAKE(I,J).gt.0) THEN
            DLAKE(I,J)  = MWL(I,J) / (RHOW*FLAKE(I,J)*AXYP(I,J))
            GLAKE(I,J)  = GML(I,J) / (FLAKE(I,J)*AXYP(I,J))
            GTEMP(I,J)  = TLAKE(I,J)
            GTEMPR(I,J) = TLAKE(I,J)+TF
            atmocn%MLHC(I,J) = SHW*MLDLK(I,J)*RHOW
#ifdef TRACERS_WATER
            GTRACER(:,I,J) = TRLAKE(:,1,I,J) / (MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J))
#endif
#ifdef SCM
            if (SCMopt%Tskin) then
               GTEMP(I,J) = SCMin%Tskin - TF
               GTEMPR(I,J) = SCMin%Tskin
            endif
#endif

!****
!****       Dump lake ice exceeding LAKE_ICE_MAX (m) into ice berg arrays
!****
            If (MSI(I,J) > LAKE_ICE_MAX * RHOW) Then     
               IMLT  = MSI(I,J) - LAKE_ICE_MAX * RHOW
               FRACI = IMLT / MSI(I,J)
               HMLT  = Sum(HSI(3:4,I,J)) * FRACI               
               PLKIC = FLAKE(I,J) * RSI(I,J)
               MDWNIMP(I,J) = MDWNIMP(I,J) + PLKIC*IMLT*AXYP(I,J)
               EDWNIMP(I,J) = EDWNIMP(I,J) + PLKIC*HMLT*AXYP(I,J)
#ifdef TRACERS_WATER
               DO ITM=1,NTM    
                  TRDWNIMP(ITM,I,J) = TRDWNIMP(ITM,I,J) + Sum(TRSI(ITM,3:4,I,J))*PLKIC*AXYP(I,J)*FRACI                   
                  TRSI(ITM,3:4,I,J) = TRSI(ITM,3:4,I,J) * (1-FRACI) 
               END DO    
#endif
!**** save some diags
               AIJ(I,J,IJ_IMPMKI) = AIJ(I,J,IJ_IMPMKI) + PLKIC*IMLT    
               AIJ(I,J,IJ_IMPHKI) = AIJ(I,J,IJ_IMPHKI) + PLKIC*HMLT    
               CALL INC_AJ(I,J,ITLKICE,J_IMPLM,PLKIC*IMLT)    
               CALL INC_AJ(I,J,ITLKICE,J_IMPLH,PLKIC*HMLT)    
!              CALL INC_AJ(I,J,ITLKICE,J_IMELT,PLKIC*IMLT)     
!              CALL INC_AJ(I,J,ITLKICE,J_HMELT,PLKIC*HMLT)     
!**** Accumulate regional diagnostics  
               JR = JREG(I,J)
!              CALL INC_AREG(I,J,JR,J_IMELT,PLKIC*IMLT)     
!              CALL INC_AREG(I,J,JR,J_HMELT,PLKIC*HMLT)     
               CALL INC_AREG(I,J,JR,J_IMPLM,PLKIC*IMLT)    
               CALL INC_AREG(I,J,JR,J_IMPLH,PLKIC*HMLT)    
!****
               MSI(I,J)     = (1-FRACI)*MSI(I,J)  ! = LAKE_ICE_MAX * RHOW
               HSI(3:4,I,J) = (1-FRACI)*HSI(3:4,I,J)
            EndIf

         Else  !  FLAKE = FLAKEnew == 0
            DLAKE(I,J)=0.
            GLAKE(I,J)=0.
         END IF
      END DO  ;  END DO
      EndSubroutine DAILY_LAKE



      SUBROUTINE PRECIP_LK
!@sum  PRECIP_LK driver for applying precipitation/melt to lake fraction
!@auth Gavin Schmidt
      Use CONSTANT, Only: rhow,shw,teeny,tf
      Use RESOLUTION, Only: im,jm
#ifdef SCM
      Use SCM_COM, Only: SCMopt,SCMin
#endif
      Use DOMAIN_DECOMP_ATM, Only: GRID,getDomainBounds
      Use GEOM, Only: imaxj,axyp,byaxyp
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use LAKES_COM, Only: mwl,gml,tlake,mldlk,flake,dlake,glake,icelak
      Use FLUXES, Only: atmocn,atmgla,prec,eprec,flice
#ifdef TRACERS_WATER
      Use LAKES_COM, Only: trlake,ntm
      Use FLUXES,    Only: trprec
#endif
      Use DIAG_COM, Only: aj=>aj_loc,j_run,aij=>aij_loc,ij_lk,itlake,itlkice
      IMPLICIT NONE

      REAL*8 PRCP,ENRGP,PLICE,PLKICE,RUN0,ERUN0,POLAKE,HLK1
      INTEGER :: J_0,J_1,J_0H,J_1H,J_0S,J_1S,I_0H,I_1H,I_0,I_1
      INTEGER I,J,ITYPE
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: TRUN0
#endif

      REAL*8, DIMENSION(:,:), POINTER :: RSI,GTEMP,GTEMP2,GTEMPR, RUNPSI,MELTI,EMELTI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:), POINTER :: GTRACER,TRUNPSI,TRMELTI
#endif

      RSI => LAKEICE%RSI
      GTEMP => ATMOCN%GTEMP
      GTEMP2 => ATMOCN%GTEMP2
      GTEMPR => ATMOCN%GTEMPR
      RUNPSI => ICELAK%RUNPSI
      MELTI => ICELAK%MELTI
      EMELTI => ICELAK%EMELTI
#ifdef TRACERS_WATER
      GTRACER => ATMOCN%GTRACER
      TRUNPSI => ICELAK%TRUNPSI
      TRMELTI => ICELAK%TRMELTI
#endif

      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      CALL PRINTLK("PR")

      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
      IF (FLAKE(I,J)+FLICE(I,J).gt.0) THEN
        POLAKE=(1.-RSI(I,J))*FLAKE(I,J)
        PLKICE=RSI(I,J)*FLAKE(I,J)
        PLICE=FLICE(I,J)
        PRCP=PREC(I,J)
        ENRGP=EPREC(I,J)        ! energy of precipitation

!**** calculate fluxes over whole box
        RUN0 = POLAKE*PRCP + PLKICE* RUNPSI(I,J) + PLICE*atmgla%RUNO(I,J)
        ERUN0=POLAKE*ENRGP ! PLKICE*ERUNPSI(I,J) + PLICE*ERUNOLI(I,J) =0

!**** simelt is given as kg/area
        IF (FLAKE(I,J).gt.0) THEN
          RUN0  =RUN0+ MELTI(I,J)
          ERUN0=ERUN0+EMELTI(I,J)
        END IF

        MWL(I,J) = MWL(I,J) +  RUN0*AXYP(I,J)
        GML(I,J) = GML(I,J) + ERUN0*AXYP(I,J)
#ifdef TRACERS_WATER
        TRUN0(:) = POLAKE*TRPREC(:,I,J) + PLKICE*TRUNPSI(:,I,J) + PLICE *atmgla%TRUNO(:,I,J)
        IF(FLAKE(I,J).gt.0) TRUN0(:)=TRUN0(:)+TRMELTI(:,I,J)
        TRLAKE(:,1,I,J)=TRLAKE(:,1,I,J) + TRUN0(:)*AXYP(I,J)
#endif

        IF (FLAKE(I,J).gt.0) THEN
          HLK1=TLAKE(I,J)*MLDLK(I,J)*RHOW*SHW
          MLDLK(I,J)=MLDLK(I,J) + RUN0/(FLAKE(I,J)*RHOW)
          TLAKE(I,J) = (HLK1*FLAKE(I,J) + ERUN0) / (MLDLK(I,J)*FLAKE(I,J)*RHOW*SHW)
          DLAKE(I,J)=MWL(I,J)/(RHOW*FLAKE(I,J)*AXYP(I,J))
          GLAKE(I,J)=GML(I,J)/(FLAKE(I,J)*AXYP(I,J))
          GTEMP(I,J)=TLAKE(I,J)
          GTEMPR(I,J) =TLAKE(I,J)+TF
#ifdef SCM
          if (SCMopt%Tskin) then
            GTEMP(I,J) = SCMin%Tskin - TF
            GTEMPR(I,J) = SCMin%Tskin
          endif
#endif
          IF (MWL(I,J).gt.(1d-10+MLDLK(I,J))*RHOW*FLAKE(I,J)*AXYP(I,J))  THEN
             GTEMP2(I,J) = (GML(I,J) - TLAKE(I,J)*SHW*MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J)) &
                         / (SHW * (MWL(I,J) - MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J)))
          ELSE
            GTEMP2(I,J)=TLAKE(I,J)
          END IF
#ifdef SCM
          if (SCMopt%Tskin) then
            GTEMP2(I,J) = GTEMP(I,J)
          endif
#endif
#ifdef TRACERS_WATER
          GTRACER(:,I,J) = TRLAKE(:,1,I,J) / (MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J))
#endif
          CALL INC_AJ(I,J,ITLAKE ,J_RUN, -PLICE*atmgla%RUNO(I,J)*(1.-RSI(I,J)))
          CALL INC_AJ(I,J,ITLKICE,J_RUN, -PLICE*atmgla%RUNO(I,J)*RSI(I,J))
        ELSE
          TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+teeny)
          DLAKE(I,J)=0.
          GLAKE(I,J)=0.
!**** accounting fix to ensure runoff with no lakes is counted
!**** no regional diagnostics required
          CALL INC_AJ(I,J,ITLAKE,J_RUN,-PLICE*atmgla%RUNO(I,J))
        END IF

!**** save area diag
        AIJ(I,J,IJ_LK) = AIJ(I,J,IJ_LK) + FLAKE(I,J)
      END IF
      END DO
      END DO
      RETURN
!****
      END SUBROUTINE PRECIP_LK



#ifdef IRRIGATION_ON
      SUBROUTINE IRRIG_LK
!@sum  IRRIG_LK driver for calculating irrigation fluxes from lakes/rivers
!@auth Gavin Schmidt
      Use CONSTANT, Only: rhow,shw,teeny
      Use RESOLUTION, Only: im,jm
      Use DOMAIN_DECOMP_ATM, Only: GRID, getDomainBounds
      Use GEOM, Only: imaxj,axyp
      Use DIAG_COM, Only: itearth,jreg,aij=>aij_loc,ij_mwlir,ij_gmlir,ij_irrgw,ij_irrgwE,j_irgw,j_irgwE
      Use LAKES_COM, Only: mwl,gml,tlake,mldlk,flake
      Use LAKES, Only: minmld,hlake_min
      Use IRRIGMOD, Only: irrigate_extract
      Use FLUXES,Only: fland,irrig_water_act, irrig_energy_act
#ifdef TRACERS_WATER
      Use LAKES_COM, Only: trlake,ntm
      Use FLUXES, Only: irrig_tracer_act
#endif
      Use TimerPackage_mod, Only: startTimer => start
      Use TimerPackage_mod, Only: stopTimer => stop
      IMPLICIT NONE
!**** grid box variables
      REAL*8 M1,M2,E1,E2,DM,DE
      REAL*8 :: MWL_to_irrig,GML_to_irrig,irrig_gw,irrig_gw_energy,irrig_water_actij,irrig_energy_actij
#ifdef TRACERS_WATER
      Real*8 :: TRML_to_irrig(NTM,2),TRML_temp(NTM,2),irrig_tracer_actij(ntm),irrig_gw_tracer(ntm)
#endif
      INTEGER I,J,JR
      INTEGER :: J_0,J_1,J_0S,J_1S,I_0,I_1

      call startTimer('PRECIP_LK()')
      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      CALL PRINTLK("IR")

      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
      JR=JREG(I,J)

!**** Remove mass/energy associated with irrigation
      IF (FLAND(I,J).gt.0) THEN

#ifdef TRACERS_WATER
        TRML_temp(:,:) = TRLAKE(:,:,I,J)
#endif
!****   Compute actual irrigation every timestep
        call irrigate_extract (I,J,MWL(I,J),GML(I,J),MLDLK(I,J),TLAKE(I,J),FLAKE(I,J),HLAKE_MIN,MWL_to_irrig,GML_to_irrig, &
                               irrig_gw,irrig_gw_energy,irrig_water_actij,irrig_energy_actij &
#ifdef TRACERS_WATER
                              ,TRML_temp,TRML_to_irrig,irrig_tracer_actij,irrig_gw_tracer &
#endif
                              )
!**** save fluxes for GHY (m/s), (J/s), (kg/s)
        irrig_water_act(i,j) =irrig_water_actij
        irrig_energy_act(i,j)=irrig_energy_actij
#ifdef TRACERS_WATER
        irrig_tracer_act(:,i,j)=irrig_tracer_actij(:)
#endif

        IF (MWL_to_irrig .gt. 0) THEN
!**** update lake mass/energy
        MWL(I,J) = MWL(I,J) - MWL_to_irrig
        GML(I,J) = GML(I,J) - GML_to_irrig
#ifdef TRACERS_WATER
        TRLAKE(:,:,I,J) = TRLAKE(:,:,I,J) - TRML_to_irrig(:,:)
        IF (MWL(I,J).eq.0) TRLAKE(:,:,I,J)=0.  ! round off issues
#endif

! mixed layer depth and surface temperature adjustments for lakes
        if (FLAKE(I,J).gt.0) THEN
          if (MWL_to_irrig.lt.MLDLK(I,J)*FLAKE(I,J)*AXYP(I,J)*RHOW) then ! layer 1 only
            MLDLK(I,J) = MLDLK(I,J) - MWL_to_irrig / (FLAKE(I,j)*AXYP(I,J)*RHOW)
            M1=MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J) ! kg
            M2=max(MWL(I,J)-M1,0d0)
            if (MLDLK(I,J).LT.MINMLD .and. M2.gt.0) THEN ! bring up from layer 2
              E1=TLAKE(I,J)*SHW*M1
              E2=GML(I,J)-E1
              DM=max(MINMLD*RHOW*FLAKE(I,J)*AXYP(I,J)-M1,0d0) ! kg
              DE=DM*E2/(M2+teeny)
              TLAKE(I,J)=(E1+DE)/((M1+DM)*SHW) ! deg C
#ifdef TRACERS_WATER
              TRLAKE(:,1,I,J) = TRLAKE(:,1,I,J) + DM*TRLAKE(:,2,I,J) / (M2+teeny)
              TRLAKE(:,2,I,J) = TRLAKE(:,2,I,J) - DM*TRLAKE(:,2,I,J) / (M2+teeny)
#endif
              MLDLK(I,J) = MLDLK(I,J) + DM/(FLAKE(I,j)*AXYP(I,J)*RHOW)
            end if
          else ! all layer 1 and some layer 2 gone, relayer
            MLDLK(I,J)=MWL(I,J)/(FLAKE(I,J)*AXYP(I,J)*RHOW)
            TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+teeny)
#ifdef TRACERS_WATER
            TRLAKE(:,1,I,J)=TRLAKE(:,1,I,J)+TRLAKE(:,2,I,J)
            TRLAKE(:,2,I,J)=0.
#endif
          end if
        end if

!****   Compute lake- and irrigation-related diagnostics
        AIJ(I,J,IJ_MWLir)=AIJ(I,J,IJ_MWLir)+MWL_to_irrig
        AIJ(I,J,IJ_GMLir)=AIJ(I,J,IJ_GMLir)+GML_to_irrig

        END IF ! MWL_to_irrig .gt. 0

        AIJ(I,J,IJ_irrgw) =AIJ(I,J,IJ_irrgw) +irrig_gw
        AIJ(I,J,IJ_irrgwE)=AIJ(I,J,IJ_irrgwE)+irrig_gw_energy

        CALL INC_AJ(I,J,itearth, j_irgw , irrig_gw)
        CALL INC_AJ(I,J,itearth, j_irgwE, irrig_gw_energy)

      END IF ! FLAND(I,J).gt.0

      END DO  ! i loop
      END DO  ! j loop

      CALL PRINTLK("I2")

      call stopTimer('PRECIP_LK()')
      RETURN
!****
      END SUBROUTINE IRRIG_LK
#endif



      SUBROUTINE GROUND_LK
!@sum  GROUND_LK driver for applying surface fluxes to lake fraction
!@auth Gavin Schmidt
!@calls
      Use CONSTANT, Only: rhow,shw,teeny,tf
      Use RESOLUTION, Only: im,jm
      Use MODEL_COM, Only: dtsrc
#ifdef SCM
      Use SCM_COM, Only: SCMopt,SCMin
#endif
      Use DOMAIN_DECOMP_ATM, Only: GRID, getDomainBounds

      Use GEOM, Only: imaxj,axyp
      Use FLUXES, Only: atmocn,atmgla,atmlnd,flice,fland
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use DIAG_COM, Only: jreg,j_wtr1,j_wtr2,j_run,j_erun,aij=>aij_loc,ij_mwl,ij_gml,itlake,itlkice,itearth
      Use LAKES_COM, Only: icelak,mwl,gml,tlake,mldlk,flake,hlake
#ifdef TRACERS_WATER
      Use LAKES_COM, Only: trlake,ntm
      Use TRDIAG_COM,Only: taijn=>taijn_loc , tij_lk1,tij_lk2
#endif
      Use LAKES, Only: lkmix,lksourc,byzeta,minmld,hlake_min
      Use GHY_COM, Only: fearth
      Use TimerPackage_mod, Only: startTimer => start
      Use TimerPackage_mod, Only: stopTimer => stop
      IMPLICIT NONE
!**** grid box variables
      REAL*8 ROICE, POLAKE, PLKICE, PEARTH, PLICE
!@var MLAKE,ELAKE mass and energy /m^2 for lake model layers
      REAL*8, DIMENSION(2) :: MLAKE,ELAKE
!**** fluxes
      REAL*8 EVAPO, FIDT, FODT, RUN0, ERUN0, RUNLI, RUNE, ERUNE, HLK1,TLK1,TLK2,TKE,SROX(2),FSR2  ! , U2RHO
!**** output from LKSOURC
      REAL*8 ENRGFO, ACEFO, ACEFI, ENRGFI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(NTM) :: TRUN0,TRO,TRI,TREVAP,TOTTRL
      REAL*8, DIMENSION(NTM,2) :: TRLAKEL
#endif
      INTEGER I,J,JR
      INTEGER :: J_0,J_1,J_0S,J_1S,I_0,I_1

      REAL*8, DIMENSION(:,:), POINTER :: RSI,GTEMP,GTEMP2,GTEMPR, RUNOSI,ERUNOSI,EVAPOR,E0
      REAL*8, DIMENSION(:,:,:), POINTER :: DMSI,DHSI,DSSI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), POINTER :: DTRSI
      REAL*8, DIMENSION(:,:,:), POINTER :: GTRACER,TREVAPOR,TRUNOSI
#ifdef TRACERS_DRYDEP
      REAL*8, DIMENSION(:,:,:), POINTER :: TRDRYDEP
#endif
#endif

      RSI => LAKEICE%RSI
      E0 => ATMOCN%E0
      EVAPOR => ATMOCN%EVAPOR
      GTEMP => ATMOCN%GTEMP
      GTEMP2 => ATMOCN%GTEMP2
      GTEMPR => ATMOCN%GTEMPR
#ifdef TRACERS_WATER
      TREVAPOR => ATMOCN%TREVAPOR
#ifdef TRACERS_DRYDEP
      TRDRYDEP => ATMOCN%TRDRYDEP
#endif
      GTRACER => ATMOCN%GTRACER
#endif
      RUNOSI => ICELAK%RUNOSI
      ERUNOSI => ICELAK%ERUNOSI
      DMSI => ICELAK%DMSI
      DHSI => ICELAK%DHSI
      DSSI => ICELAK%DSSI
#ifdef TRACERS_WATER
      TRUNOSI => ICELAK%TRUNOSI
      DTRSI => ICELAK%DTRSI
#endif

      call startTimer('GROUND_LK()')
      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      CALL PRINTLK("GR")

      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
      JR=JREG(I,J)
      ROICE=RSI(I,J)
      PLKICE=FLAKE(I,J)*ROICE
      POLAKE=FLAKE(I,J)*(1.-ROICE)
!**** Add land ice and surface runoff to lake variables
      IF (FLAND(I,J).gt.0) THEN
        PLICE =FLICE(I,J)
        PEARTH=FEARTH(I,J)
        RUNLI=atmgla%RUNO(I,J)
        RUNE =atmlnd%RUNO(I,J)
        ERUNE=atmlnd%ERUNO(I,J)
!**** calculate flux over whole box
        RUN0 =RUNLI*PLICE + RUNE*PEARTH
        ERUN0=             ERUNE*PEARTH
        MWL(I,J) = MWL(I,J) + RUN0*AXYP(I,J)
        GML(I,J) = GML(I,J) +ERUN0*AXYP(I,J)
#ifdef TRACERS_WATER
        TRLAKE(:,1,I,J) = TRLAKE(:,1,I,J) + (atmgla%TRUNO(:,I,J)*PLICE + atmlnd%TRUNO(:,I,J)*PEARTH)*AXYP(I,J)
#endif

        AIJ(I,J,IJ_MWL)=AIJ(I,J,IJ_MWL)+MWL(I,J)
        AIJ(I,J,IJ_GML)=AIJ(I,J,IJ_GML)+GML(I,J)

        IF (FLAKE(I,J).gt.0) THEN
          HLK1=TLAKE(I,J)*MLDLK(I,J)*RHOW*SHW
          MLDLK(I,J)=MLDLK(I,J) + RUN0/(FLAKE(I,J)*RHOW)
          TLAKE(I,J) = (HLK1*FLAKE(I,J)+ERUN0) / (MLDLK(I,J)*FLAKE(I,J)*RHOW*SHW)
#ifdef TRACERS_WATER
          GTRACER(:,I,J) = TRLAKE(:,1,I,J) / (MLDLK(I,J)*RHOW*FLAKE(I,J)*AXYP(I,J))
#endif
          CALL INC_AJ(I,J,ITLAKE ,J_RUN ,-(RUNE*PEARTH+RUNLI*PLICE)*(1.-RSI(I,J)))
          CALL INC_AJ(I,J,ITLKICE,J_RUN ,-(RUNE*PEARTH+RUNLI*PLICE)*    RSI(I,J))
          CALL INC_AJ(I,J,ITLAKE ,J_ERUN,-ERUNE*PEARTH*(1.-RSI(I,J)))
          CALL INC_AJ(I,J,ITLKICE,J_ERUN,-ERUNE*PEARTH*    RSI(I,J))
        ELSE
          TLAKE(I,J)=GML(I,J)/(MWL(I,J)*SHW+teeny)
!**** accounting fix to ensure runoff with no lakes is counted
!**** no regional diagnostics required
          CALL INC_AJ(I,J,ITLAKE,J_RUN, -(RUNE*PEARTH+RUNLI*PLICE))
          CALL INC_AJ(I,J,ITLAKE,J_ERUN,-ERUNE*PEARTH)
        END IF
      END IF

      IF (FLAKE(I,J).gt.0) THEN
        TLK1 =TLAKE(I,J)
        EVAPO=EVAPOR(I,J)     ! evap/dew over open lake (kg/m^2)
        FODT =E0(I,J)         ! net heat over open lake (J/m^2)
        SROX(1)=atmocn%SOLAR(I,J)      ! solar radiation open lake (J/m^2)
        SROX(2)=icelak%SOLAR(I,J)      ! solar radiation through ice (J/m^2)
        FSR2 =EXP(-MLDLK(I,J)*BYZETA)
!**** get ice-lake fluxes from sea ice routine (over ice fraction)
        RUN0 =RUNOSI(I,J) ! includes ACE2M + basal term
        FIDT =ERUNOSI(I,J)
!**** calculate kg/m^2, J/m^2 from saved variables
        MLAKE(1)=MLDLK(I,J)*RHOW
        MLAKE(2)=MAX(MWL(I,J)/(FLAKE(I,J)*AXYP(I,J))-MLAKE(1),0d0)
        ELAKE(1)=TLK1*SHW*MLAKE(1)
        ELAKE(2)=GML(I,J)/(FLAKE(I,J)*AXYP(I,J))-ELAKE(1)
#ifdef TRACERS_WATER
        TRLAKEL(:,:)=TRLAKE(:,:,I,J)/(FLAKE(I,J)*AXYP(I,J))
        TRUN0(:)=TRUNOSI(:,I,J)
        TREVAP(:)=TREVAPOR(:,I,J)
#ifdef TRACERS_DRYDEP
        TREVAP(:) = TREVAP(:) - trdrydep(:,i,j)
#endif
#endif
        IF (MLAKE(2).lt.1d-10) THEN
          MLAKE(1)=MLAKE(1)+MLAKE(2)
          MLAKE(2)=0.
          ELAKE(1)=ELAKE(1)+ELAKE(2)
          ELAKE(2)=0.
#ifdef TRACERS_WATER
          TRLAKEL(:,1)=TRLAKEL(:,1)+TRLAKEL(:,2)
          TRLAKEL(:,2)=0.
#endif
        END IF

!**** Limit FSR2 in the case of thin second layer
        FSR2=MIN(FSR2,MLAKE(2)/(MLAKE(1)+MLAKE(2)))

!**** Apply fluxes and calculate the amount of frazil ice formation
        CALL LKSOURC (I,J,ROICE,MLAKE,ELAKE,RUN0,FODT,FIDT,SROX,FSR2, &
#ifdef TRACERS_WATER
                      TRLAKEL,TRUN0,TREVAP,TRO,TRI, &
#endif
                      EVAPO,ENRGFO,ACEFO,ACEFI,ENRGFI)

!**** Mixing and entrainment
!**** Calculate turbulent kinetic energy for lake
!       U2rho=(1.-ROICE)*SQRT(DMUA(I,J,1)**2+DMVA(I,J,1)**2)/DTSRC
!       TKE=0.5 * (19.3)^(2/3) * U2rho /rhoair ! (m/s)^2
        TKE=0.  ! 3.6d0*U2rho/rhoair*MLAKE(1)  ! (J/m^2)

        CALL LKMIX (MLAKE,ELAKE, &
#ifdef TRACERS_WATER
                    TRLAKEL, &
#endif
                    HLAKE(I,J),TKE,ROICE,DTSRC)

!**** Resave prognostic variables
        MWL(I,J)  =(MLAKE(1)+MLAKE(2))*(FLAKE(I,J)*AXYP(I,J))
        GML(I,J)  =(ELAKE(1)+ELAKE(2))*(FLAKE(I,J)*AXYP(I,J))
        MLDLK(I,J)= MLAKE(1)/RHOW
        IF (MLAKE(2).eq.0.) MLDLK(I,J)=MIN(MINMLD,MLDLK(I,J))
        TLAKE(I,J)= ELAKE(1)/(SHW*MLAKE(1))
        IF (MLAKE(2).gt.0) THEN
          TLK2    = ELAKE(2)/(SHW*MLAKE(2))
        ELSE
          TLK2    = TLAKE(I,J)
        END IF
#ifdef TRACERS_WATER
        IF (MLAKE(2).eq.0. .and. MLAKE(1)-MLDLK(I,J)*RHOW.gt.1d-10) THEN
          TOTTRL(:)=TRLAKEL(:,1)
          TRLAKEL(:,2)=(MLAKE(1)-MLDLK(I,J)*RHOW)*TRLAKEL(:,1)/MLAKE(1)
          TRLAKEL(:,1)=TOTTRL(:)-TRLAKEL(:,2)
        END IF
        TRLAKE(:,:,I,J)=TRLAKEL(:,:)*(FLAKE(I,J)*AXYP(I,J))
        GTRACER(:,I,J)=TRLAKEL(:,1)/(MLDLK(I,J)*RHOW)
#endif
        GTEMP(I,J)=TLAKE(I,J)
        GTEMP2(I,J)=TLK2       ! diagnostic only
        GTEMPR(I,J) =TLAKE(I,J)+TF
#ifdef SCM
        if (SCMopt%Tskin) then
          GTEMP(I,J) = SCMin%Tskin - TF
          GTEMP2(I,J) = SCMin%Tskin - TF
          GTEMPR(I,J) = SCMin%Tskin
        endif
#endif
!**** Open lake diagnostics
        CALL INC_AJ(I,J, ITLAKE,J_WTR1,MLAKE(1)*POLAKE)
        CALL INC_AJ(I,J, ITLAKE,J_WTR2,MLAKE(2)*POLAKE)
!**** Ice-covered ocean diagnostics
        CALL INC_AJ(I,J, ITLKICE,J_WTR1,MLAKE(1)*PLKICE)
        CALL INC_AJ(I,J, ITLKICE,J_WTR2,MLAKE(2)*PLKICE)
!**** regional diags
        CALL INC_AREG(I,J,JR,J_WTR1,MLAKE(1)*FLAKE(I,J))
        CALL INC_AREG(I,J,JR,J_WTR2,MLAKE(2)*FLAKE(I,J))
#ifdef TRACERS_WATER
!**** tracer diagnostics
        TAIJN(I,J,tij_lk1,:)=TAIJN(I,J,tij_lk1,:)+TRLAKEL(:,1) !*PLKICE?
        TAIJN(I,J,tij_lk2,:)=TAIJN(I,J,tij_lk2,:)+TRLAKEL(:,2) !*PLKICE?
#endif

!**** Store mass and energy fluxes for formation of sea ice
        DMSI(1,I,J)=ACEFO
        DMSI(2,I,J)=ACEFI
        DHSI(1,I,J)=ENRGFO
        DHSI(2,I,J)=ENRGFI
        DSSI(:,I,J)=0.     ! always zero salinity
#ifdef TRACERS_WATER
        DTRSI(:,1,I,J)=TRO(:)
        DTRSI(:,2,I,J)=TRI(:)
#endif
      END IF
      END DO  ! i loop
      END DO  ! j loop

      CALL PRINTLK("G2")

      call stopTimer('GROUND_LK()')
      RETURN
!****
      END SUBROUTINE GROUND_LK



      SUBROUTINE PRINTLK(STR)
!@sum  PRINTLK print out selected diagnostics from specified lakes
!@auth Gavin Schmidt
      Use CONSTANT, Only: lhm,byshi,rhow,shw
      Use MODEL_COM, Only: qcheck
      Use GEOM, Only: axyp
      Use LAKES_COM, Only: tlake,mwl,mldlk,gml,flake
#ifdef TRACERS_WATER
      Use LAKES_COM, Only: trlake
#endif
      Use SEAICE, Only: xsi,ace1i,rhoi
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use DOMAIN_DECOMP_ATM, Only: GRID, getDomainBounds
      IMPLICIT NONE
      CHARACTER*2, INTENT(IN) :: STR
      INTEGER, PARAMETER :: NDIAG=4
      INTEGER I,J,N, J_0, J_1
      INTEGER, DIMENSION(NDIAG) :: IDIAG = (/112, 103, 131, 79/), &
                                   JDIAG = (/ 66,  59,  33, 34/)
      REAL*8 HLK2,TLK2, TSIL(4)

      REAL*8, DIMENSION(:,:), POINTER :: RSI,MSI,SNOWI
      REAL*8, DIMENSION(:,:,:), POINTER :: HSI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI
#endif

      RSI => LAKEICE%RSI
      MSI => LAKEICE%MSI
      HSI => LAKEICE%HSI
      SNOWI => LAKEICE%SNOWI
#ifdef TRACERS_WATER
      TRSI => LAKEICE%TRSI
#endif

      IF (.NOT.QCHECK) RETURN

      call getDomainBounds(grid, J_STRT=J_0,      J_STOP=J_1)

      DO N=1,NDIAG
        I=IDIAG(N)
        J=JDIAG(N)
        if (J.lt. J_0 .or. J.gt. J_1) CYCLE
        IF (FLAKE(I,J).gt.0) THEN
          HLK2 = MWL(I,J)/(RHOW*FLAKE(I,J)*AXYP(I,J)) - MLDLK(I,J)
          IF (HLK2.gt.0) THEN
            TLK2 = (GML(I,J) / (SHW*RHOW*FLAKE(I,J)*AXYP(I,J)) - TLAKE(I,J)*MLDLK(I,J)) / HLK2
          ELSE
            TLK2=0.
          END IF
          TSIL(:)=0.
          IF (RSI(I,J).gt.0) THEN
            TSIL(1:2) = (HSI(1:2,I,J) / (XSI(1:2)*(ACE1I+SNOWI(I,J))) + LHM) * bySHI
            TSIL(3:4) = (HSI(3:4,I,J)/(XSI(3:4)*MSI(I,J))+LHM)*BYSHI
          END IF
          WRITE(99,*) STR,I,J,FLAKE(I,J),TLAKE(I,J),TLK2,MLDLK(I,J),HLK2,RSI(I,J),MSI(I,J)/RHOI,SNOWI(I,J)/RHOW,TSIL(1:4)
#ifdef TRACERS_WATER
          Write (99,*) TRLAKE(1,1:2,I,J),MWL(I,J)
#endif
        ELSE
          WRITE(99,*) STR,I,J,TLAKE(I,J),MWL(I,J)
#ifdef TRACERS_WATER
          Write (99,*) TRLAKE(1,1:2,I,J)
#endif
        END IF
      END DO

      RETURN
      END  SUBROUTINE PRINTLK



      SUBROUTINE conserv_LKM(LKM)
!@sum  conserv_LKM calculates lake mass
!@auth Gary Russell/Gavin Schmidt
      Use RESOLUTION, Only: im,jm
      Use FLUXES, Only: fland
      Use DOMAIN_DECOMP_ATM, Only: GRID, getDomainBounds
      Use GEOM, Only: imaxj,byaxyp
      Use LAKES_COM, Only: mwl,flake
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: LKM
      INTEGER :: I,J
      INTEGER :: J_0,J_1,J_0S,J_1S,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S, &
                            HAVE_SOUTH_POLE = HAVE_SOUTH_POLE, HAVE_NORTH_POLE = HAVE_NORTH_POLE )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!****
!**** LAKE MASS (kg/m^2)
!****
      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
        IF (FLAND(I,J)+FLAKE(I,J).gt.0) THEN
          LKM(I,J)=MWL(I,J)*BYAXYP(I,J)
        ELSE
          LKM(I,J)=0.
        ENDIF
      ENDDO
      ENDDO
      IF (HAVE_SOUTH_POLE) LKM(2:im,1) =LKM(1,1)
      IF (HAVE_NORTH_POLE) LKM(2:im,JM)=LKM(1,JM)
      RETURN
      END SUBROUTINE conserv_LKM



      SUBROUTINE conserv_LKE(LKE)
!@sum  conserv_LKE calculates lake energy
!@auth Gary Russell/Gavin Schmidt
      Use RESOLUTION, Only: im,jm
      Use ATM_COM, Only: zatmo
      Use FLUXES, Only: fland
      Use DOMAIN_DECOMP_ATM, Only: GRID, getDomainBounds
      Use GEOM, Only: imaxj,byaxyp
      Use LAKES_COM, Only: gml,mwl,flake
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: LKE
      INTEGER :: I,J
      INTEGER :: J_0,J_1,J_0S,J_1S,I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S, &
                            HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!****
!**** LAKE ENERGY (J/m^2) (includes potential energy (DISABLED))
!****
        DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          IF (FLAND(I,J)+FLAKE(I,J).gt.0) THEN
            LKE(I,J) = GML(I,J)*byAXYP(I,J)  !  + ZATMO(I,J)*MWL(I,J)
          ELSE
            LKE(I,J)=0.
          ENDIF
        END DO
      END DO
      IF (HAVE_SOUTH_POLE) LKE(2:im,1) =LKE(1,1)
      IF (HAVE_NORTH_POLE) LKE(2:im,JM)=LKE(1,JM)
      RETURN
      END SUBROUTINE conserv_LKE



      subroutine diag_river_prep
      Use constant, Only: rhow
      Use domain_decomp_atm, Only: grid,getDomainBounds,sumxpe
      Use constant, Only: rhow
      Use model_com, Only: dtsrc, calendar
      Use TimeConstants_mod, Only: INT_MONTHS_PER_YEAR
      Use diag_com, Only: aij=>aij_loc,ij_mrvr
      Use lakes_com, Only: irvrmth,jrvrmth,nrvrmx,nrvr,rvrout
      Use Rational_mod
      implicit none
      real*8 rvrout_loc(nrvrmx), scalervr
      integer inm,i,j
      integer :: i_0, i_1, j_0, j_1
      type (Rational) :: secondsPerYear

      if(nrvr.lt.1) return
      call getDomainBounds(grid, j_strt=j_0, j_stop=j_1)
      i_0 = grid%i_strt
      i_1 = grid%i_stop
!**** convert kg/(source time step) to km^3/mon
      secondsPerYear = calendar%getMaxDaysInYear() * calendar%getSecondsPerDay()
      SCALERVR = 1d-9*real(secondsPerYear) / (INT_MONTHS_PER_YEAR*RHOW*DTSRC)

!**** fill in the river discharges in the local domain
      rvrout_loc(:)=0
      do j=j_0,j_1
      do i=i_0,i_1
        do inm=1,nrvr
          if (i.eq.irvrmth(inm).and. j.eq.jrvrmth(inm)) then
            rvrout_loc(inm) = scalervr*aij(i,j,ij_mrvr)
          end if
        end do
      end do
      end do
!**** sum over processors to compose the global table
      call sumxpe(rvrout_loc, rvrout)
      return
      end subroutine diag_river_prep



      SUBROUTINE init_lakeice(iniLAKE,do_IC_fixups)
!@sum  init_ice initialises ice arrays
!@auth Original Development Team
      Use CONSTANT, Only: rhows,omega
      Use MODEL_COM, Only: kocean
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use LAKES_COM, Only: icelak
      Use FLUXES, Only: flake0,atmice
      Use Dictionary_mod
      Use DOMAIN_DECOMP_ATM, Only: GRID, getDomainBounds
      Use GEOM, Only: sinlat2d
      Use DIAG_COM, Only: npts,conpt0,icon_LMSI,icon_LHSI
      IMPLICIT NONE
      LOGICAL :: QCON(NPTS), T=.TRUE. , F=.FALSE. , iniLAKE
      CHARACTER CONPT(NPTS)*10
      INTEGER I,J,do_IC_fixups
      integer :: I_0, I_1, J_0, J_1
!****
!**** Extract useful local domain parameters from "grid"
!****
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

!**** clean up ice fraction/sea ice salinity possibly incorrect in I.C.
      if (do_IC_fixups == 1) then
        DO J=J_0, J_1
        DO I=I_0, I_1
          IF (FLAKE0(I,J).eq.0 .and. lakeice%RSI(i,j).gt.0)  lakeice%RSI(I,J) = 0
          IF (lakeice%RSI(I,J).gt.0 .and. FLAKE0(I,J).gt.0)  lakeice%SSI(:,I,J) = 0
        END DO
        END DO
      end if

      DO J=J_0,J_1
      DO I=I_0,I_1
        icelak%coriol(i,j) = ABS(2.*OMEGA*SINLAT2D(I,J))
      ENDDO
      ENDDO

      IF (KOCEAN.EQ.0.and.iniLAKE) THEN
        ! why should lake ice init depend on kocean,iniocean?
        call set_noice_defaults(lakeice,icelak)
      END IF

      !call seaice_to_atmgrid(atmice) ! set gtemp etc.

!**** Set conservation diagnostics for Lake ice mass and energy
      CONPT=CONPT0
      CONPT(3)="LAT. MELT" ; CONPT(4)="PRECIP"
      CONPT(5)="THERMO"
      CONPT(8)="LK FORM"
      QCON=(/ F, F, T, T, T, F, F, T, T, F, F/)
      CALL SET_CON (QCON,CONPT,"LKICE MS","(KG/M^2)        ", "(10**-9 KG/SM^2)",1d0,1d9,icon_LMSI)
      QCON=(/ F, F, T, T, T, F, F, T, T, F, F/)
      CALL SET_CON (QCON,CONPT,"LKICE EN","(10**6 J/M^2)   ", "(10**-3 W/M^2)  ",1d-6,1d3,icon_LHSI)

      END SUBROUTINE init_lakeice



      SUBROUTINE conserv_LMSI(ICE)
!@sum  conserv_LMSI calculates total amount of snow and ice over lakes
!@auth Gavin Schmidt
      Use RESOLUTION, Only: im,jm
      Use GEOM, Only: imaxj
      Use SEAICE, Only: ace1i
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use LAKES_COM, Only: flake
      Use DOMAIN_DECOMP_ATM, Only: GRID,getDomainBounds
      IMPLICIT NONE
!@var ICE total lake snow and ice mass (kg/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: ICE
      INTEGER I,J

!**** Extract useful domain information from grid
      INTEGER J_0, J_1, I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds (GRID, J_STRT=J_0, J_STOP=J_1, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ICE(I,J) = lakeice%RSI(I,J) * (lakeice%MSI(I,J) + ACE1I + lakeice%SNOWI(I,J)) * FLAKE(I,J)
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) ICE(2:im,1) =ICE(1,1)
      IF (HAVE_NORTH_POLE) ICE(2:im,JM)=ICE(1,JM)
      RETURN
!****
      END SUBROUTINE conserv_LMSI



      SUBROUTINE conserv_LHSI(EICE)
!@sum  conserv_LHSI calculates total ice energy over lakes
!@auth Gavin Schmidt
      Use RESOLUTION, Only: im,jm
      Use GEOM, Only: imaxj
      Use SEAICE_COM, Only: lakeice=>si_atm
      Use LAKES_COM, Only: flake
      Use DOMAIN_DECOMP_ATM, Only: GRID,getDomainBounds
      IMPLICIT NONE
!@var EICE total lake snow and ice energy (J/m^2)
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO, GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: EICE
      INTEGER I,J

!**** Extract useful domain information from grid
      INTEGER J_0, J_1, I_0,I_1
      LOGICAL HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds (GRID, J_STRT=J_0, J_STOP=J_1, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        EICE(I,J)=lakeice%RSI(I,J)*FLAKE(I,J)*SUM(lakeice%HSI(:,I,J))
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) EICE(2:im,1) =EICE(1,1)
      IF (HAVE_NORTH_POLE) EICE(2:im,JM)=EICE(1,JM)
      RETURN
!****
      END SUBROUTINE conserv_LHSI



      SUBROUTINE CHECKI(SUBR)
!@sum  CHECKI Checks whether Ice values are reasonable
!@auth Original Development Team
      Use MODEL_COM
      Use GEOM, Only: imaxj
#ifdef TRACERS_WATER
      Use OldTracer_mod, Only: trname, t_qlimit
      Use TRACER_COM, Only: NTM
#endif
      Use SEAICE, Only: lmi,xsi,ace1i,Ti,Ti2b
      Use SEAICE_COM, Only: x=>si_atm
      Use LAKES_COM, Only: flake
      Use FLUXES
      Use DOMAIN_DECOMP_ATM, Only: GRID
      Use DOMAIN_DECOMP_ATM, Only: getDomainBounds
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
!@var QCHECKI true if errors found in seaice
      LOGICAL QCHECKI
      INTEGER I,J,L
      REAL*8 TICE
#ifdef TRACERS_WATER
      integer :: imax,jmax, n
      real*8 relerr,errmax
#endif

      integer :: J_0, J_1, J_0H, J_1H, I_0, I_1, I_0H, I_1H, njpol
      REAL*8 MSI1,SNOWL(2),MICE(2)
!****
!**** Extract useful local domain parameters from "grid"
!****
      call getDomainBounds (grid, J_STRT=J_0, J_STOP=J_1, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      njpol = grid%J_STRT_SKP-grid%J_STRT

!**** Check for NaN/INF in ice data
      CALL CHECK3B (x%RSI(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'rsi   ')
      CALL CHECK3B (x%MSI(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'msi   ')
      CALL CHECK3C (x%HSI(:,I_0:I_1,J_0:J_1),LMI,I_0,I_1,J_0,J_1,NJPOL, SUBR,'hsi   ')
      CALL CHECK3C (x%SSI(:,I_0:I_1,J_0:J_1),LMI,I_0,I_1,J_0,J_1,NJPOL, SUBR,'ssi   ')
      CALL CHECK3B (x%SNOWI(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1, SUBR,'sni   ')

      QCHECKI = .FALSE.
!**** Check for reasonable values for ice variables
      DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          IF (x%RSI(I,J).lt.0 .or. x%RSI(I,j).gt.1 .or. x%MSI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,RSI,MSI=',I,J,x%RSI(I,J),x%MSI(I,J)
            QCHECKI = .TRUE.
          END IF
          IF ( (FOCEAN(I,J)+FLAKE(I,J))*x%RSI(I,J).gt.0) THEN
          MSI1 = ACE1I + x%SNOWI(I,J)
          IF (ACE1I.gt.XSI(2)*MSI1) THEN ! some ice in first layer
            MICE(1) = ACE1I-XSI(2)*MSI1
            MICE(2) = XSI(2)*MSI1
!           SNOWL(1)= SNOW
            SNOWL(1)= x%SNOWI(I,J)
            SNOWL(2)= 0.
          ELSE  ! some snow in second layer
            MICE(1) = 0.
            MICE(2) = ACE1I
            SNOWL(1)= XSI(1)*MSI1
            SNOWL(2)= XSI(2)*MSI1-ACE1I
          ENDIF
          DO L=1,LMI
            IF (L.EQ.1) THEN
              IF(MICE(1).NE.0.) THEN
                TICE = Ti2b (x%HSI(1,I,J)/(XSI(1)*MSI1), 1d3*x%SSI(L,I,J)/MICE(1), SNOWL(1), MICE(1))
              ELSE
                TICE = Ti(x%HSI(1,I,J)/(XSI(1)*MSI1),0d0)
              ENDIF
            ENDIF
            IF (L.EQ.2)  TICE = Ti2b (x%HSI(2,I,J)/(XSI(2)*MSI1), 1d3*x%SSI(L,I,J)/MICE(2), SNOWL(2), MICE(2))

            IF (L.gt.2) TICE = Ti (x%HSI(L,I,J)/(XSI(L)*x%MSI(I,J)), 1d3*x%SSI(L,I,J)/(XSI(L)*x%MSI(I,J)))
            IF (x%HSI(L,I,J).gt.0.or.TICE.gt.1d-10.or.TICE.lt.-80.) THEN
              WRITE (6,'(3a,3i3,6e12.4/1X,6e12.4)') 'After ',SUBR,': I,J,L,TSI=',I,J,L,TICE,x%RSI(I,J)
            WRITE(6,*) x%HSI(:,I,J),x%MSI(I,J),x%SNOWI(I,J),x%SSI(:,I,J)
              IF (TICE.gt.1d-3.or.TICE.lt.-100.) QCHECKI = .TRUE.
            END IF
            IF (x%SSI(L,I,J).lt.0) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,SSI=',I,J,L,x%SSI(:,I,J),x%MSI(I,J),x%SNOWI(I,J),x%RSI(I,J)
              QCHECKI = .TRUE.
            END IF
           IF (L.gt.2 .and. x%SSI(L,I,J).gt.0.04*XSI(L)*x%MSI(I,J)) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,SSI/MSI=',I,J,L,1d3*x%SSI(:,I,J)/(XSI(L)*x%MSI(I,J)),x%SSI(:,I,J),x%MSI(I,J), &
                          x%SNOWI(I,J),x%RSI(I,J)
              QCHECKI = .TRUE.
            END IF
          END DO
          IF (x%SNOWI(I,J).lt.0) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,SNOWI=',I,J,x%SNOWI(I,J)
            QCHECKI = .TRUE.
          END IF
          IF (x%MSI(I,J).gt.10000) THEN
         WRITE(6,*) 'After ',SUBR,': I,J,MSI=',I,J,x%MSI(I,J),x%RSI(I,J)
!            QCHECKI = .TRUE.
          END IF
          END IF
        END DO
      END DO

#ifdef TRACERS_WATER
      do n=1,ntm
!**** check negative tracer mass
        if (t_qlimit(n)) then
        do j=J_0, J_1
          do i=I_0,imaxj(j)
            if ((focean(i,j)+flake(i,j))*x%rsi(i,j).gt.0) then
              do l=1,lmi
                if (x%trsi(n,l,i,j).lt.0.) then
                  print*,"Neg Tracer in sea ice after ",subr,i,j,l,trname(n),x%trsi(n,l,i,j),x%rsi(i,j), &
                         x%msi(i,j),x%ssi(l,i,j),x%snowi(i,j)
                  QCHECKI=.true.
                end if
              end do
            end if
          end do
        end do
        end if
!**** Check conservation of water tracers in sea ice
        if (trname(n).eq.'Water') then
          errmax = 0. ; imax=I_0 ; jmax=J_0
          do j=J_0, J_1
          do i=I_0,imaxj(j)
            if ((focean(i,j)+flake(i,j))*x%rsi(i,j).gt.0) then
              relerr = max (abs(x%trsi(n,1,i,j)-(x%snowi(i,j)+ace1i)*xsi(1)+x%ssi(1,i,j))/x%trsi(n,1,i,j), &
                            abs(x%trsi(n,2,i,j)-(x%snowi(i,j)+ace1i)*xsi(2)+x%ssi(2,i,j))/x%trsi(n,2,i,j), &
                            abs(x%trsi(n,3,i,j)-x%msi(i,j)*xsi(3)+x%ssi(3,i,j))/x%trsi(n,3,i,j), &
                            abs(x%trsi(n,4,i,j)-x%msi(i,j)*xsi(4)+x%ssi(4,i,j))/x%trsi(n,4,i,j))
              if (relerr.gt.errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            end if
          end do
          end do
          write (*,'(A36,A7,A,2I3,11E24.16)') "Relative error in sea ice mass after",trim(subr),":", &
                imax,jmax,errmax,x%trsi(n,:,imax,jmax), &
                (x%snowi(imax,jmax)+ace1i)*xsi(1)-x%ssi(1,imax,jmax), &
                (x%snowi(imax,jmax)+ace1i)*xsi(2)-x%ssi(2,imax,jmax), &
                x%msi(imax,jmax)*xsi(3:4)-x%ssi(3:4,imax,jmax),x%rsi(imax,jmax),x%msi(imax,jmax)
        end if
      end do
#endif

      IF (QCHECKI)  call stop_model ("CHECKI: Ice variables out of bounds", 255)

      END SUBROUTINE CHECKI
