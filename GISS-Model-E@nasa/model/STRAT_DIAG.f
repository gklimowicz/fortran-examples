!@sum  STRAT_DIAG file for special E-P flux diagnostics from strat model
!@auth B. Suozzo/J/ Lerner

      SUBROUTINE EPFLUX
#ifndef CUBED_SPHERE
     &     (U,V,T,P)
#endif
!@sum  EPFLUX calculates finite difference EP Flux on B-grid
!@auth B. Suozzo/J. Lerner
C**** B-grid
C**** INPUT:
C****     U - Zonal wind (corners) (m s-1)
C****     V - Meridional wind (corners) (m s-1)
C****     W - dp/dt (downward) (P-T grid level edges) (mb m2 s-1)
C****                (L=1 is ground level)
C****     T - Potential Temperature (K)  real T = TH*(P(mb))**K
C****     P - Pressure (mb) of column from 100mb to surface
C**** OUTPUT:
C****     FMY(j,l),FEY(j,l) - North. flux of mean,eddy ang. mom.   m3 s-2
C****     FMZ(j,l),FEZ(j,l) - Vert. flux of mean,eddy ang. mom.    m3 mb s-2
C****     COR(j,l),CORR(j,l) - Coriolis term and transf. coriolis term  m3 s-2
C****     FER1(j,l) - North. flux of error term 1   m2 s-2
C****     ER21(j,l),ER22(j,l) - error term 2 parts 1 & 2  m s-2
C****     FMYR(j,l),FEYR(j,l) - Nor. flux of transf. mean,eddy ang. mom. m3 s-2
C****  FMZR(j,l),FEZR(j,l) - Vert. flx of transf. mean,eddy ang. mom. m3 mb s-2
C****     RX(j,l) - <v'th'>/<dth/dp> (U-V grid level edges)  (mb m s-1)
C****                (L=1 is ground level)
C****
C**** Note: W(,,1) is really PIT, the pressure tendency, but is not used
C****
      USE RESOLUTION, only : lm
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, HALO_UPDATE, SOUTH,
     &     NORTH, HALO_UPDATEj
#ifdef CUBED_SPHERE
      USE RESOLUTION, only : ls1,psfmpt
      USE GCDIAG, only : grid,im=>imlon,jm=>jmlat,byim,jl_dpb,
     &     dxv,rapvn,rapvs,fcor,dxyv,cosv,cosp,dxyp
      USE GCDIAG, only : cs2llinta,cs2llintb
      use cs2ll_utils, only : cs2llint_ij,cs2llint_ijl,cs2llint_lluv_3d
      use atm_com, only : pcs=>p,tcs=>t,ualij,valij
      use domain_decomp_atm, only : grid_cs=>grid
      USE DYNAMICS, only : dsig,conv
      use geom, only : byaxyp
#else
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE RESOLUTION, only : im,jm
      USE GEOM, only : dxv,rapvn,rapvs,fcor,dxyv,cosv,cosp
      Use DYNAMICS,  Only: SD
      USE DIAG_COM, only : byim
#endif
      USE DIAG_COM, only : pl=>plm
      USE GC_COM, only : agc=>agc_loc,kagc,kep
      USE ATM_COM, only : pdsigl00
      IMPLICIT NONE

      Real*8,Dimension(:,:,:),Allocatable :: W
#ifdef CUBED_SPHERE
      real*8, dimension(:,:,:), allocatable :: u,v,t, uijl,vijl,wcs
      real*8, dimension(:,:), allocatable :: p
      integer :: I_0cs,I_1cs,J_0cs,J_1cs, I_0Hcs,I_1Hcs,J_0Hcs,J_1Hcs
#else
      REAL*8, INTENT(INOUT), 
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        U,V,T
      REAL*8, INTENT(IN), 
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: P
#endif
C**** ARRAYS CALCULATED HERE:
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        UV,UW,FD
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *        STB,TI,AX
     *     ,VR,WR,RX,UI,VI,WI
      REAL*8  RXCL(LM),UCL(LM), WXXS(IM),WXXN(IM)

      REAL*8 :: XEP(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM,KEP)

c Use CPP to alias into XEP
#define FMY(j,k)  XEP(j,k, 1)
#define FEY(j,k)  XEP(j,k, 2)
#define FMZ(j,k)  XEP(j,k, 3)
#define FEZ(j,k)  XEP(j,k, 4)
#define FMYR(j,k) XEP(j,k, 5)
#define FEYR(j,k) XEP(j,k, 6)
#define FMZR(j,k) XEP(j,k, 7)
#define FEZR(j,k) XEP(j,k, 8)
#define COR(j,k)  XEP(j,k, 9)
#define CORR(j,k) XEP(j,k,10)
#define FER1(j,k) XEP(j,k,11)
#define ER21(j,k) XEP(j,k,12)
#define ER22(j,k) XEP(j,k,13)
c#define VR(j,k)   XEP(j,k,14)
c#define WR(j,k)   XEP(j,k,15)
c#define RX(j,k)   XEP(j,k,16)
c#define UI(j,k)   XEP(j,k,17)
c#define VI(j,k)   XEP(j,k,18)
c#define WI(j,k)   XEP(j,k,19)
c#define DUT(j,k)  XEP(j,k,20)
c#define DUD(j,k)  XEP(j,k,21)


      INTEGER I,J,L,N,IM1,IP1
      REAL*8 UCB,UCT,UCM,ALPH,RXC

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0H, J_1H, J_0S, J_1S, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)


      Allocate (W(IM,J_0H:J_1H,LM))
#ifdef CUBED_SPHERE
      call getDomainBounds(grid_cs,
     &     I_STRT=I_0cs,I_STOP=I_1cs,
     &     J_STRT=J_0cs,J_STOP=J_1cs,
     &     I_STRT_HALO=I_0Hcs,I_STOP_HALO=I_1Hcs,
     &     J_STRT_HALO=J_0Hcs,J_STOP_HALO=J_1Hcs )
      allocate(
     &         uijl(I_0Hcs:I_1Hcs,J_0Hcs:J_1Hcs,lm),
     &         vijl(I_0Hcs:I_1Hcs,J_0Hcs:J_1Hcs,lm),
     &         wcs(I_0Hcs:I_1Hcs,J_0Hcs:J_1Hcs,lm) )
      allocate(u(im,j_0h:j_1h,lm),v(im,j_0h:j_1h,lm))
      allocate(t(im,j_0h:j_1h,lm))
      allocate(p(im,j_0h:j_1h))
      do l=1,lm
        do j=j_0cs,j_1cs
        do i=i_0cs,i_1cs
          uijl(i,j,l) = ualij(l,i,j)
          vijl(i,j,l) = valij(l,i,j)
          wcs(i,j,l) = conv(i,j,l)*byaxyp(i,j)
        enddo
        enddo
      enddo
      call cs2llint_ij(grid_cs,cs2llinta,pcs,p)
      call cs2llint_ijl(grid_cs,cs2llinta,tcs,t)
      call cs2llint_ijl(grid_cs,cs2llinta,wcs,w)
      call cs2llint_lluv_3d(grid_cs,cs2llintb,uijl,vijl,u,v)

      have_domain: if(grid%have_domain) then

      do l=1,ls1-1
        do j=j_0,j_1
          agc(j,l,jl_dpb) = agc(j,l,jl_dpb) + dsig(l)*sum(p(:,j))
        enddo
      enddo
      do l=ls1,lm
        do j=j_0,j_1
          agc(j,l,jl_dpb) = agc(j,l,jl_dpb) + dsig(l)*im*psfmpt
        enddo
      enddo
      do l=1,lm
        do j=j_0,j_1
          w(:,j,l) = w(:,j,l)*dxyp(j)
        enddo
      enddo

#else
      W(:,:,1) = 0
      W(:,:,2:LM) = SD(:,:,1:LM-1)
#endif


C     Use IDACC(4) for calling frequency

C**** Initialise for this call
      xep = 0

C****
C**** ZONAL AVERAGE QUANTITIES
C****
C**** UI(J,L) ... ZONAL AVERAGE U-WIND (m s-1)
C**** VI(J,L) ... ZONAL AVERAGE V-WIND (m s-1)
C**** TI(J,L) ... ZONAL AVERAGE POTENTIAL TEMPERATURE (K)
C**** WI(J,L) ... ZONAL AVERAGE VERTICAL WIND (mb m2 s-1)

      CALL HALO_UPDATE(grid, U, from = NORTH+SOUTH)
      CALL HALO_UPDATE(grid, V, from = NORTH)
      CALL HALO_UPDATE(grid, W, from = SOUTH)
      CALL HALO_UPDATE(grid, T, from = SOUTH)

      CALL AVGVI (U,UI(:,:))
      CALL AVGVI (V,VI(:,:))
      CALL AVGI (W,WI(:,:))
      CALL AVGI (T,TI(:,:))

C****
C**** STB(J,L) ... DELTA THETA             (PRESSURE EDGES)
C****
      DO L=2,LM
        DO J=J_0,J_1
          STB(J,L) = TI(J,L)-TI(J,L-1)
          IF(STB(J,L).LT.1d-2) THEN
CW         IF(STB(J,L).LT.1d-3)
CW   *     WRITE (*,*) 'STB < .001 AT J,L,STB:',J,L,STB(J,L)
            STB(J,L)=1d-2
          ENDIF
        END DO
      END DO
      STB(:,1)  = 1d-2

      CALL HALO_UPDATEj(grid, VI, from = NORTH)
      CALL HALO_UPDATEj(grid, UI, from = NORTH+SOUTH)
      CALL HALO_UPDATEj(grid, WI, from = SOUTH)
      CALL HALO_UPDATEj(grid, TI, from = SOUTH)
      CALL HALO_UPDATEj(grid, STB, from = SOUTH)

C****
C**** CORRELATIONS OF VARIOUS QUANTITIES
C****
C****
C**** FMY(j,l) .............. (VUdx) NORTHWARD MEAN MOMENTUM TRANSPORT
C****                         (m3 s-1)
      DO L=1,LM
        
        IF (HAVE_SOUTH_POLE) FMY(1,L)=0.
        IF (HAVE_NORTH_POLE) FMY(JM,L)=0.
        DO J=J_0S,J_1S
          FMY(J,L)=.25*(DXV(J)*VI(J,L)+DXV(J+1)*VI(J+1,L))
     *         * (UI(J,L)+UI(J+1,L))
        END DO
C****
C**** FEY(j,l), UV(I,J,L) ... (V'U'dx) NORTHWARD EDDY MOMENTUM TRANSPORT
C****                           (m3 s-2)
        DO I=1,IM
          IF (HAVE_SOUTH_POLE) UV(I,1,L)=0.
          IF (HAVE_NORTH_POLE) UV(I,JM,L)=0.
        END DO
        DO J=J_0S,J_1S
          IM1=IM-1
          I=IM
          DO IP1=1,IM
            UV(I,J,L) = (
     *        DXV(J)*(V(IM1,J,L)+2*V(I,J,L)+V(IP1,J,L)-4*VI(J,L))+
     *        DXV(J+1)*(V(IM1,J+1,L)+2*V(I,J+1,L)+V(IP1,J+1,L)-
     *        4*VI(J+1,L)))* 
     *       (U(I,J,L)+U(I,J+1,L)-(UI(J,L)+UI(J+1,L)))
            IM1=I
            I=IP1
          END DO
        END DO
      END DO
      CALL AVGI (UV,FEY(:,:))
      DO L=1,LM
      DO J=J_0,J_1
        FEY(J,L)=.0625d0*FEY(J,L)
      END DO
      END DO
C****
C**** FMZ(j,l) .............. (WUdA) VERTICAL MEAN MOMENTUM TRANSPORT
C**** FEZ(j,l), UW(I,J,L) ... (U'W'dA) VERTICAL EDDY MOMENTUM TRANSPORT
C****                           (m3 mb s-2)
      UW(:,:,1)=0.
      DO L=2,LM
        WXXS(1:IM)=0.
        IF (HAVE_SOUTH_POLE) THEN
           UW(1:IM,1,L)=0.
        ELSE ! need to seed WXXS from lower neighbor
           I=IM
           DO IP1=1,IM
              WXXS(I) = (W(I,J_0-1,L)+W(IP1,J_0-1,L))-2*WI(J_0-1,L)
              I=IP1
           END DO
        END IF
        DO J=J_0STG,J_1STG
          FMZ(J,L) = (WI(J-1,L)*RAPVN(J-1)+WI(J,L)*RAPVS(J))
     *         * (UI(J,L-1)+UI(J,L))
          I=IM
          DO IP1=1,IM
            WXXN(I)=(W(I,J,L)+W(IP1,J,L))-2*WI(J,L)
            UW(I,J,L) = (WXXS(I)*RAPVN(J-1)+WXXN(I)*RAPVS(J))
     *           * ((U(I,J,L-1)+U(I,J,L))-(UI(J,L-1)+UI(J,L)))
            WXXS(I)=WXXN(I)
            I=IP1
          END DO
        END DO
      END DO
      CALL AVGI (UW,FEZ(:,:))
      FMZ(:,1)=0.
      FEZ(:,1)=0.
      If (HAVE_SOUTH_POLE) THEN
        FMZ(1,2:LM)=0.
        FEZ(1,2:LM)=0.
      End If
      DO L=1,LM
      DO J=J_0,J_1
        FEZ(J,L)=.5*FEZ(J,L)
      END DO
      END DO
C****
C**** COR(J,L) .......... CORIOLIS FORCE
C****                        (m3 s-2)
      DO L=1,LM
        IF (HAVE_SOUTH_POLE) 
     *      COR( 2,L)=.5*(2.*FCOR( 1)+FCOR(   2))*VI(2,L)
        IF (HAVE_NORTH_POLE) 
     *      COR(JM,L)=.5*(2.*FCOR(JM)+FCOR(JM-1))*VI(JM,L)
        DO J=MAX(3,J_0S),J_1S
          COR(J,L)=.5*(FCOR(J-1)+FCOR(J))*VI(J,L)
        END DO
      END DO

cBMP FD calc moved out of next loop for ghosting purposes
      DO L=1,LM
        IM1=IM
        DO I=1,IM
          IF (HAVE_SOUTH_POLE) FD(I,1,L)= 0.
          IF (HAVE_NORTH_POLE) FD(I,JM,L)=0.
          DO J=J_0S,J_1S
            FD(I,J,L) = (DXV(J)-DXV(J+1)) *
     *            .25*(U(IM1,J,L)+U(I,J,L)+U(IM1,J+1,L)+U(I,J+1,L))
          END DO
          IM1=I
        END DO
      END DO

      CALL HALO_UPDATE(grid, FD , from = SOUTH)
cBMP

C****
C**** ERROR TERMS
C****
      DO L=1,LM
        FER1(:,L)=0.
        ER21(:,L)=0.
        ER22(:,L)=0.
C**** ERROR TERM 1 IS DIVERGENCE OF FER1(j,l): (m2 s-2)
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          IF (HAVE_SOUTH_POLE) 
     *        FER1(1,L) =FER1(1,L) +(U(I,2,L)**2   - U(I,2,L)**2)
          IF (HAVE_NORTH_POLE)
     *        FER1(JM,L)=FER1(JM,L)+(U(IM1,JM,L)**2- U(IP1,JM,L)**2)
          DO J=J_0S,J_1S
C? FER1(j,l) does not include COSP, yet ER1 = 1/COSV*(FER1(J,)-FER1(J-1,))
            FER1(J,L)= FER1(J,L) + ((U(IM1,J,L)+U(I,J+1,L))**2
     *           - (U(IP1,J,L)+U(I,J+1,L))**2)
          END DO
          IM1=I
          I=IP1
        END DO
C**** ERROR TERM 2 - THE METRIC TERM AS IN RUN 999
        IM1=IM
        DO I=1,IM
c         IF (HAVE_SOUTH_POLE) FD(I,1,L)= 0.
c         IF (HAVE_NORTH_POLE) FD(I,JM,L)=0.
c         DO J=J_0S,J_1S
c           FD(I,J,L) = (DXV(J)-DXV(J+1)) * 
c    *           .25*(U(IM1,J,L)+U(I,J,L)+U(IM1,J+1,L)+U(I,J+1,L))
c         END DO
          DO J=J_0STG,J_1STG
            ALPH=.25*(FD(I,J-1,L)+FD(I,J,L))
            ER22(J,L)=ER22(J,L)+ALPH*(V(IM1,J,L)+V(I,J,L))
          END DO
          IM1=I
        END DO
CE*** ERROR TERM 2, PART 2  -  ER22(j,l):  (m s-2)
CE    IM1=IM-1
CE    I=IM
CE    DO IP1=1,IM
CE    IF (HAVE_SOUTH_POLE) UXXNS=(U(IM1,2,L)+2*U(I,2,L)+U(IP1,2,L))
CE    IF (HAVE_SOUTH_POLE) XUXXYS=2*UXXNS*DXP(2)
CE    DO J=J_0S,J_1S
CE    UXXNN=(U(IM1,J+1,L)+2*U(I,J+1,L)+U(IP1,J+1,L))
CE    XUXXYN=(UXXNS+UXXNN)*(DXP(J+1)-DXP(J-1))
CE    ER22(J,L) = ER22(J,L)-V(I,J,L)*(XUXXYS+XUXXYN)
CE    XUXXYS=XUXXYN
CE    UXXNS=UXXNN
CE    END DO
CE    IF (HAVE_NORTH_POLE) XUXXYN=-2*UXXNS*DXP(JM-1)
CE    IF (HAVE_NORTH_POLE) ER22(JM,L) = ER22(JM,L)-V(I,JM,L)*(XUXXYS+XUXXYN)
CE    IM1=I
CE    I=IP1
CE    END DO
C**** NORMALIZE ERROR TERMS AND CONVERT TO OUTPUT UNITS
      DO J=J_0,J_1
        FER1(J,L)=1./(48*IM)*FER1(J,L)
      END DO
CE      DO J=J_0STG,J_1STG
CE        ER21(J,L) = -1./(2*DXYV(J)*COSV(J))*
CE   *  (FMY(J-1,L)+FEY(J-1,L)+FMY(J,L)+FEY(J,L))*(COSP(J)-COSP(J-1))
CE        ER22(J,L) = 1./(32*IM*DXYV(J))*ER22(J,L)
CE      END DO
      DO J=J_0STG,J_1STG
        ER22(J,L) = 1./(IM*DXYV(J))*ER22(J,L)
      END DO
      END DO
C****
C****  TRANSFORMED CIRCULATION
C****

C****
C**** RX(J,L) ..... ([V'TH']/[DTH/DP]) TRANSFORMATION GAUGE
C****               (U-wind grid, level edges)  (m mb s-1)

      DO L=2,LM
      DO J=J_0STG,J_1STG
        RX(J,L) = 0.
        I=IM
        DO IP1=1,IM
          RX(J,L) = RX(J,L)+(V(I,J,L)+V(IP1,J,L)-2*VI(J,L))
     *      * (T(I,J-1,L)+T(I,J,L) - (TI(J-1,L)+TI(J,L)))
     *      +        (V(I,J,L-1)+V(IP1,J,L-1)-2*VI(J,L-1))
     *       * (T(I,J-1,L-1)+T(I,J,L-1)-
     *         (TI(J-1,L-1)+TI(J,L-1)))
          I=IP1
        END DO
        RX(J,L) = .25*BYIM * RX(J,L) *
     *       (PL(L)-PL(L-1))/(STB(J-1,L)+STB(J,L))
      END DO
      END DO
      DO J = J_0,J_1
        RX(J,1) = 0
      END DO

      CALL HALO_UPDATEj(grid, RX, FROM=NORTH+SOUTH)

C****
C**** VR(J,L) ........  TRANSFORMED WIND
C**** WR(J,L) ........
C****
      DO L=1,LM-1
      DO J=J_0STG,J_1STG
        VR(J,L)=VI(J,L)+(RX(J,L+1)-RX(J,L))/PDSIGL00(L)
!     *       (PSFMPT*DSIG(L))
      END DO
      END DO
      DO J=J_0STG,J_1STG
        VR(J,LM)=VI(J,LM)-RX(J,LM)/PDSIGL00(LM) !(PSFMPT*DSIG(LM))
      END DO
      DO L=1,LM
      DO J=J_0S,J_1S
        WR(J,L)=WI(J,L)+
     *       (RX(J+1,L)*DXV(J+1)-RX(J,L)*DXV(J))
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) THEN
         DO L=1,LM
            WR(1,L)=WI(1,L) + RX(2,L)*DXV(2)
         END DO
      ENDIF
      IF (HAVE_NORTH_POLE) THEN
         DO L=1,LM
            WR(JM,L)=WI(JM,L) - RX(JM,L)*DXV(JM)
         END DO
      ENDIF

C****
C**** TRANSFORMED MEAN FLUXES
C****
C**** FMYR(J,L)     (m3 s-2)
C**** FMZR(J,L)     (m3 mb s-2)
C**** CORR(J,L)     (m3 s-2)
C****

      CALL HALO_UPDATEj(grid, VR, from = NORTH)
      CALL HALO_UPDATEj(grid, WR, from = SOUTH)
      DO L=1,LM

      DO J=J_0S,J_1S
        FMYR(J,L)=.25*(VR(J,L)*DXV(J)+VR(J+1,L)*DXV(J+1))
     *       *  (UI(J,L)+UI(J+1,L))
      END DO
      END DO
      DO L=2,LM
      DO J=J_0STG,J_1STG
        FMZR(J,L) = (WR(J-1,L)*RAPVN(J-1)+WR(J,L)*RAPVS(J))
     *       *  (UI(J,L-1)+UI(J,L))
      END DO
      END DO
      DO L=1,LM
        IF (HAVE_SOUTH_POLE) 
     *      CORR( 2,L)=.5*(2.*FCOR( 1)+FCOR(   2))*VR(2,L)
        IF (HAVE_NORTH_POLE) 
     *      CORR(JM,L)=.5*(2.*FCOR(JM)+FCOR(JM-1))*VR(JM,L)
        DO J=MAX(3,J_0S),J_1S
          CORR(J,L)=.5*(FCOR(J-1)+FCOR(J))*VR(J,L)
        END DO
      END DO
C****
C**** TRANSFORMED EDDY FLUXES
C****
C**** FEYR(J,L)
C****

      IF (HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          RXCL(L)=RX(2,L)*COSV(2)
          UCL(L)=UI(2,L)*DXV(2)
        END DO
      ELSE
        DO L=1,LM
          RXCL(L)=RX(J_0,L)*COSV(J_0)
          UCL(L)=UI(J_0,L)*DXV(J_0)
        END DO
      END IF

      DO J=J_0S,J_1S
        UCB=(UI(J+1,1)*DXV(J+1)+UCL(1))
        UCM=(UI(J+1,2)*DXV(J+1)+UCL(2))
        FEYR(J,1)=0.
        FEYR(J,LM)=0.
        DO L=2,LM-1
          UCT=(UI(J+1,L+1)*DXV(J+1)+UCL(L+1))
c          FEYR(J,L) = .125d0/(COSP(J)*PSFMPT*DSIG(L))
          FEYR(J,L) = .125d0/(COSP(J)*PDSIGL00(L))
     *         *    ((RX(J+1,L)*COSV(J+1)+RXCL(L))*(UCM-UCB)
     *         +    (RX(J+1,L+1)*COSV(J+1)+RXCL(L+1))*(UCT-UCM))
          UCB=UCM
          UCM=UCT
        END DO
        DO L=1,LM
          RXCL(L)=RX(J+1,L)*COSV(J+1)
          UCL(L)=UI(J+1,L)*DXV(J+1)
        END DO
      END DO

C****
C**** FEZR(j,l)
C****
      If (HAVE_SOUTH_POLE) THEN
         RXCL(:)=0.
         UCL(:)=0.
      Else
         RXCL(:) = RX(J_0-1,:)*COSV(J_0-1)
         UCL(:)  = UI(J_0-1,:)* DXV(J_0-1)
      END If
      
      DO J=J_0S,J_1S
      DO L=2,LM
        RXC = RX(J,L)*COSV(J)
        UCB = (UI(J,L)+UI(J,L-1))*DXV(J)
        FEZR(J,L) = -.5*(FCOR(J-1)+FCOR(J))*RX(J,L)
     *       + .125d0/COSV(J)* ((RXC+RXCL(L))
     *       *   ((UI(J,L)+UI(J,L-1))*DXV(J)-(UCL(L)+UCL(L-1)))
     *       +   (RX(J+1,L)*COSV(J+1)+RXC)
     *       *   ((UI(J+1,L)+UI(J+1,L-1))*DXV(J+1)-UCB))
      END DO
      DO L=1,LM
        RXCL(L)=RX(J,L)*COSV(J)
        UCL(L)=UI(J,L)*DXV(J)
      END DO
      END DO
      IF (HAVE_NORTH_POLE) THEN
         DO L=2,LM
           RXC=RX(JM,L)*COSV(JM)
           UCB=(UI(JM,L)+UI(JM,L-1))*DXV(JM)
           FEZR(JM,L) = -FCOR(JM)*RX(JM,L)
     *          + .125d0/COSV(JM)* (RXC+RXCL(L))
     *          *     ((UI(JM,L)+UI(JM,L-1))*DXV(JM) 
     *          -      (UCL(L)+UCL(L-1)))
          END DO
      ENDIF
C**** Add Eulerian circulation to transformed eddy components
      DO L=1,LM
      DO J=J_0,J_1
        FEYR(J,L)=FEY(J,L)+FEYR(J,L)
        FEZR(J,L)=FEZ(J,L)+FEZR(J,L)
      END DO
      END DO

C****
C**** ACCUMULATE EP FLUXES
C****
      DO N=1,KEP-2
      DO L=1,LM
      DO J=J_0,J_1
        AGC(J,L,KAGC-KEP+N)=AGC(J,L,KAGC-KEP+N)+XEP(J,L,N)
c        AEP(J,L,N)=AEP(J,L,N)+XEP(J,L,N)
      END DO
      END DO
      END DO

      Deallocate (W)
#ifdef CUBED_SPHERE
      endif have_domain
      deallocate(u,v,t,p, uijl,vijl,wcs)
#endif

      RETURN
      END SUBROUTINE EPFLUX
C****
#ifndef CUBED_SPHERE
      subroutine EPFLXI(U)
C****     U - Zonal wind (corners) (m s-1)
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE RESOLUTION, only : im,jm,lm
      USE GC_COM, only : agc=>agc_loc,KAGC
      REAL*8, INTENT(INOUT), 
     *        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     *                                                        U

      integer :: j_0h
      call getDomainBounds(grid, j_strt_halo=j_0h)
      call avgvi (u,agc(j_0h,1,KAGC))
      return
C****
      end subroutine epflxi
#endif


      SUBROUTINE EPFLXP(do_print,DUDS,DMF,DEF,DMFR,DEFR,ER1,ER2)
!@sum  EPFLXP prints out diagnostics of E-P Fluxes
!@auth B. Suozzo/J. Lerner
C****
C**** B-grid
C**** INPUT:
C****     AEP contains the following quantities summed over time:
C****     FMY(j,l),FEY(j,l) - North. flux of mean,eddy ang. mom.   m3 s-2
C****     FMZ(j,l),FEZ(j,l) - Vert. flux of mean,eddy ang. mom.    m3 mb s-2
C****     COR(j,l),CORR(j,l) - coriolis term and transf. coriolis term  m3 s-2
C****     FER1(j,l) - North. flux of error term 1   m2 s-2
C****     ER21(j,l),ER22(j,l) - error term 2 parts 1 & 2  m s-2
C****     FMYR(j,l),FEYR(j,l) - Nor. flux of transf. mean,eddy ang. mom. m3 s-2
C****  FMZR(j,l),FEZR(j,l) - Vert. flx of transf. mean,eddy ang. mom. m3 mb s-2
C**** OUTPUT:
C****        DMF,DEF - Divergence of mean,eddy ang. mom.  m3 s-2
C****        DMFR,DEFR - Div. of transf. mean,eddy ang. mom.  m3 s-2
C****        COR(j,l),CORR(j,l) - same as above   m3 s-2
C****        ER1,ER2 - error terms 1 and 2   m s-2
C****   DUD(j,l),DUR - Delta U by Eulerian and transf. circulation  m s-2
C****
      USE RESOLUTION, only : lm
      USE MODEL_COM, only : dtsrce=>dtsrc,idacc
      USE DYNAMICS, only : dsig
      USE DIAG_COM, only : ndaa,ajl,jl_damdc,jl_dammc
     &     ,jl_dudfmdrg,jl_dumtndrg,jl_dushrdrg,jl_dudtsdif,jl_dudtvdif
     &     ,jl_dumcdrgm10,jl_dumcdrgp10,jl_dumcdrgm20,jl_dumcdrgp20
     &     ,jl_dumcdrgm40,jl_dumcdrgp40
      USE GC_COM, only : kagc,kep,agc
     &     ,lname=>lname_gc,sname=>sname_gc,units=>units_gc,pow=>pow_gc
      USE ATM_COM, only : pmidl00,pdsigl00
#ifdef CUBED_SPHERE
      USE GCDIAG, only : im=>imlon,jm=>jmlat,byim,fim,
     &     dxv,cosv,cosp,dxyp,dyv,bydxyv
#else
      USE DIAG_SERIAL, only : JLMAP
      USE RESOLUTION, only : im,jm
      USE GEOM, only : dxv,cosv,cosp,dyv,bydxyv
      USE DIAG_COM, only : fim,byim
#endif
      USE GCDIAG, only : jl_dpb,jk_dudt_sum1,jk_dudt_meanadv,
     &     jk_dudt_eddycnv,jk_dudt_trnsadv,jk_dudt_epflxdiv,
     &     jk_dudt_fderr1,jk_dudt_fderr2
      IMPLICIT NONE
      logical :: do_print
      REAL*8, DIMENSION(JM,LM) :: ! output arrays
     &     DUDS,DMF,DEF,DMFR,DEFR,ER1,ER2

C**** diagnostic information for print out
C**** this should be in an init_ep routine or something
c      integer, parameter :: njl_out=7
c      character(len=lname_strlen) :: lname(njl_out) = (/
c     *     'DU/DT BY EULER CIRC. + CONVEC + DRAG+DIF+ER2',
c     *     'DU/DT BY MEAN ADVECTION                     ',
c     *     'DU/DT BY EDDY CONVERGENCE                   ',
c     *     'DU/DT BY TRANSFORMED ADVECTION              ',
c     *     'DU/DT BY ELIASSEN-PALM DIVERGENCE           ',
c     *     'DU/DT BY F.D. ERROR TERM 1                  ',
c     *     'DU/DT BY F.D. ERROR TERM 2                  '/)
c      character(len=sname_strlen) :: sname(njl_out) = (/
c     *     'dudt_sum1    ','dudt_meanadv ','dudt_eddycnv '
c     *     ,'dudt_trnsadv ','dudt_epflxdiv','dudt_fderr1  '
c     *     ,'dudt_fderr2  '/)
c      character(len=units_strlen) :: units(njl_out) = 'm/s^2'
c      integer, dimension(njl_out) :: pow = (/ -6,-6,-6,-6,-6,-6,-6 /)

C**** ARRAYS CALCULATED HERE:
      REAL*8 ONES(JM+LM),PMO(LM),DP(LM)
      REAL*8 DXCOSV(JM)
      REAL*8, DIMENSION(JM,LM) :: DUR,BYDPJL,DUD

      REAL*8,DIMENSION(JM,LM,KEP) :: XEP
      REAL*8 DTAEP,BYDT,SCALEP,SCALE1,BYIAEP
      INTEGER I,J,L,N,JL

C**** Initialize constants
      DTAEP = DTsrce*NDAA   ! change of definition of NDAA
      BYIAEP=1./(IDACC(4)+1.D-20)
      BYDT=1./(DTSRCE*IDACC(1)+1.D-20)
      DO J=2,JM
        DXCOSV(J) = DXV(J)*COSV(J)
      END DO

      J = 1
      DO L=1,LM
        DUDS(J,L) = 0.
        DMF(J,L) = 0.
        DEF(J,L) = 0.
        DMFR(J,L) = 0.
        DEFR(J,L) = 0.
        ER1(J,L) = 0.
        ER2(J,L) = 0.
      ENDDO

c GISS-ESMF EXCEPTIONAL CASE: OK AS A CONSTANT ON EACH PE
      DO JL=1,JM+LM
        ONES(JL)=1.
      END DO
      DO L=1,LM
        DP(L)  = PDSIGL00(L)
        PMO(L) = PMIDL00(L)
      END DO
!     do j=2,jm
!       ap=0.25*APJ(J,2)/(FIM*IDACC(4)+teeny)
!       call calc_vert_amp(ap,lm,PL,AML,PDSIGL,PEDNL,PMIDL)
!       BYDPJL(J,1:LM)=1./PDSIGL(1:LM)
!     end do
      DO L=1,LM
      DO J=2,JM
        BYDPJL(J,L)=(FIM*IDACC(4))/(AGC(J,L,JL_DPB)+1.D-20)
      END DO
      END DO
C**** Normalize AEP
C     CALL EPFLXF (U)
      DO N=1,KEP-2
      DO L=1,LM
      DO J=1,JM
        XEP(J,L,N)=AGC(J,L,KAGC-KEP+N)*BYIAEP
c        XEP(J,L,N)=AEP(J,L,N)*BYIAEP
      END DO
      END DO
      END DO
      DO N=KEP-1,KEP
      DO L=1,LM
      DO J=1,JM
        XEP(J,L,N)=AGC(J,L,KAGC-KEP+N)
c        XEP(J,L,N)=AEP(J,L,N)
      END DO
      END DO
      END DO
C****
C**** DMF,DEF by horizontal convergence
C****

      DO L=1,LM
      DO J=2,JM
        DMF(J,L) =(FMY(J-1,L)*COSP(J-1)-FMY(J,L)*COSP(J))/COSV(J)
        DEF(J,L) =(FEY(J-1,L)*COSP(J-1)-FEY(J,L)*COSP(J))/COSV(J)
        DMFR(J,L)=(FMYR(J-1,L)*COSP(J-1)-FMYR(J,L)*COSP(J))/COSV(J)
        DEFR(J,L)=(FEYR(J-1,L)*COSP(J-1)-FEYR(J,L)*COSP(J))/COSV(J)
      END DO
      END DO
C****
C**** Add DMF,DEF by vertical convergence
C****
      DO J=2,JM
        DMF(J,LM)  = DMF(J,LM)  - FMZ(J,LM)/BYDPJL(J,LM)
        DEF(J,LM)  = DEF(J,LM)  - FEZ(J,LM)/BYDPJL(J,LM)
        DMFR(J,LM) = DMFR(J,LM) - FMZR(J,LM)/BYDPJL(J,LM)
        DEFR(J,LM) = DEFR(J,LM) - FEZR(J,LM)/BYDPJL(J,LM)
        DO L=1,LM-1
          DMF(J,L)  = DMF(J,L)+(FMZ(J,L+1)-FMZ(J,L))*BYDPJL(J,L)
          DEF(J,L)  = DEF(J,L)+(FEZ(J,L+1)-FEZ(J,L))*BYDPJL(J,L)
          DMFR(J,L) = DMFR(J,L) + 
     &         (FMZR(J,L+1)-FMZR(J,L))*BYDPJL(J,L)
          DEFR(J,L) = DEFR(J,L) + 
     &         (FEZR(J,L+1)-FEZR(J,L))*BYDPJL(J,L)
        END DO
      END DO
C****
C**** ADD COR(j,l),CORR(j,l)
C****
      DO L=1,LM
        DO J=2,JM
          DMF(J,L)  = DMF(J,L)  + COR(J,L)
          DMFR(J,L) = DMFR(J,L) + CORR(J,L)
        END DO
      END DO
C****
C**** ER1 by horizontal convergence (m s-2)
C****

      DO L=1,LM
      DO J=2,JM
        ER1(J,L) = 1./(DYV(J)*COSV(J))*(FER1(J-1,L)-FER1(J,L))
        ER2(J,L) = ER21(J,L)+ER22(J,L)
      END DO
      END DO
C****
C**** DUD(j,l),DUR total change in Eulerian, Transformed wind  (m s-2)
C****
      DO L=1,LM
      DO J=2,JM
        DUR(J,L)  = BYDXYV(J)*(DMFR(J,L) +DEFR(J,L))
        DUD(J,L)  = BYDXYV(J)*(DMF(J,L)+DEF(J,L))
      END DO
      END DO
C****
C**** Print maps of EP fluxes
C**** note: JLMAP (lname,sname,units,power,Pres,Array,ScalP,ScalJ,ScalL)
C****          prints maps of  (Array * SCALEP * ScaleJ * ScaleL)
C****
#ifndef CUBED_SPHERE
      DO L=1,LM
      DO J=2,JM
        DUDS(J,L)=((AJL(J,L,JL_DUDFMDRG)+AJL(J,L,JL_DUMTNDRG))+
     &             (AJL(J,L,JL_DUSHRDRG)+AJL(J,L,JL_DUMCDRGM10))+
     *       ((AJL(J,L,JL_DUMCDRGP10)+AJL(J,L,JL_DUMCDRGM40))+
     &        (AJL(J,L,JL_DUMCDRGP40)+AJL(J,L,JL_DUMCDRGM20)))+
     *        (AJL(J,L,JL_DUMCDRGP20)+
     &         AJL(J,L,JL_DUDTSDIF)+AJL(J,L,JL_DUDTVDIF)))*FIM
      END DO
      END DO
      SCALE1=1./(FIM*DTSRCE*IDACC(1)+1.D-20)
      DO L=1,LM
      DO J=2,JM
        DUDS(J,L) = DUDS(J,L)+
     &        (AJL(J,L,JL_DAMDC)+AJL(J,L,JL_DAMMC))*BYDPJL(J,L)
      END DO
      END DO
      DO L=1,LM
      DO J=2,JM
        DUDS(J,L) = (DUD(J,L)-ER2(J,L)) + DUDS(J,L)*SCALE1
      END DO
      END DO
#endif

      do l=1,lm
        do j=2,jm
          dmf(j,l) = dmf(j,l)*bydxyv(j)
          def(j,l) = def(j,l)*bydxyv(j)
          dmfr(j,l) = dmfr(j,l)*bydxyv(j)
          defr(j,l) = defr(j,l)*bydxyv(j)
        enddo
      enddo

#ifndef CUBED_SPHERE
      if(do_print) then
      SCALEP=1
CW      CALL WRITJL ('DUDT: EUL+SOURCE',DUDS,SCALEP)
CW     /* CALL WRITJL ('DUDT: ENTIRE GCM',DUT,SCALEP) ! AJK-47 DIAGJK */
      n = jk_dudt_sum1
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,DUDS,SCALEP,ONES,ONES,LM,2,2)
      n = jk_dudt_meanadv
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,DMF,SCALEP,ONES,ONES,LM,2,2)
      n = jk_dudt_eddycnv
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,DEF,SCALEP,ONES,ONES,LM,2,2)
      n = jk_dudt_trnsadv
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,DMFR,SCALEP,ONES,ONES,LM,2,2)
      n = jk_dudt_epflxdiv
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,DEFR,SCALEP,ONES,ONES,LM,2,2)
      n = jk_dudt_fderr1
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,ER1,SCALEP,ONES,ONES,LM,2,2)
      n = jk_dudt_fderr2
      CALL JLMAP (LNAME(n),SNAME(n),UNITS(n),POW(n)
     *     ,PMO,ER2,SCALEP,ONES,ONES,LM,2,2)
      endif ! do_print
#endif /* not CUBED_SPHERE */
C****
CW      DO L=1,LM
CW      DO J=J_0STG,J_1STG
CW      DMF(J,L)=DMF(J,L)*BYDXYV(J)
CW      DEF(J,L)=DEF(J,L)*BYDXYV(J)
CW      DMFR(J,L)=DMFR(J,L)*BYDXYV(J)
CW      DEFR(J,L)=DEFR(J,L)*BYDXYV(J)
CW      END DO
CW      END DO
CW      WRITE (36,'(''DU/DT = 10**-6 M S-2'')')
CW      WRITE (36,'(''TR: TRANSFORMED - EULERIAN'')')
CW      SCALEP = 1.E6
CW      CALL WRITJL ('DUDT: MEAN EULER',DMF,SCALEP)
CW      CALL WRITJL ('DUDT: EDDY EULER',DEF,SCALEP)
CW      /* CALL WRITJL ('DUDT: EULER CIRC',DUD,SCALEP) */
CW      CALL WRITJL ('DUDT: MEAN TRANS',DMFR,SCALEP)
CW      CALL WRITJL ('DUDT: EDDY TRANS',DEFR,SCALEP)
CW      CALL WRITJL ('DUDT: TRANS-EULE',DUR,SCALEP)
      RETURN
      END SUBROUTINE EPFLXP

      SUBROUTINE AVGI (X,XI)
!@sum  AVGI average a 3-dimensional array in the x-direction
!@auth B. Suozzo
      USE RESOLUTION, only : lm
#ifdef CUBED_SPHERE
      USE GCDIAG, only : grid,im=>imlon,BYIM
#else
      USE RESOLUTION, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE GEOM, only : imaxj
      USE DIAG_COM, only : BYIM
#endif
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      IMPLICIT NONE

!@var X input 3-D array
      REAL*8, INTENT(IN), 
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: X
!@var XI output zonally averaged 2-D array
      REAL*8, INTENT(OUT), 
     &        DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: XI
      INTEGER I,J,L,NI
      REAL*8 XXI

      INTEGER :: J_1,J_0

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

      DO L=1,LM
        DO J=J_0,J_1
          XXI=0.
#ifdef CUBED_SPHERE
          NI = IM
#else
          NI = IMAXJ(J)
#endif
          DO I=1,NI
            XXI = XXI + X(I,J,L)
          END DO
          XI(J,L) = XXI/NI
        END DO
      END DO

      RETURN
      END SUBROUTINE AVGI

      SUBROUTINE AVGVI (X,XI)
!@sum  AVGVI average a 3-dimensional array in the x-direction (no pole)
!@auth B. Suozzo
      USE RESOLUTION, only : lm
#ifdef CUBED_SPHERE
      USE GCDIAG, only : grid,im=>imlon,BYIM
#else
      USE RESOLUTION, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DIAG_COM, only : BYIM
#endif
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      IMPLICIT NONE

!@var X input 3-D array
      REAL*8, INTENT(IN), 
     &        DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: X
!@var XI output zonally averaged 2-D array
      REAL*8, INTENT(OUT), 
     &        DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: XI
      INTEGER I,J,L
      REAL*8 XXI

      INTEGER :: J_0STG, J_1STG

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG)
      DO L=1,LM
        DO J=J_0STG,J_1STG
          XXI=0.
          DO I=1,IM
            XXI = XXI + X(I,J,L)
          END DO
          XI(J,L) = XXI*BYIM
        END DO
      END DO
C****
      RETURN
      END SUBROUTINE AVGVI
