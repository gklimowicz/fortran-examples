#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif

      SUBROUTINE STPGF(DTS)
C****
!@sum  STPGF calculates pressure gradient force through strait
!@auth Gary Russell
!@ver  2009/05/22
C****
      USE CONSTANT, only : grav,rrt12=>byrt12
      Use OCEAN,   Only: LMO, DXYPO
      Use STRAITS, Only: NMST,LMST,IST,JST,WIST,MUST,DISTPG,
     *                   LMME,OPRESE,HOCEANE, MOE, G0ME,GZME, S0ME,SZME
      IMPLICIT NONE

      INTEGER I,J,K,L,N
      REAL*8, INTENT(IN) :: DTS
      REAL*8 PE,PHIE,GUP,GDN,SUP,SDN,DP,PUP,PDN,VUP,VDN
      REAL*8 VOLGSP
      Real*8,Dimension(LMO,2) :: PEND,PHI,DH
C****
      DO N=1,NMST
C****
C**** Accelerate the channel speed due to the pressure gradient
C**** force, reduce the channel speed due to friction
C****
      DO K=1,2
      I=IST(N,K)
      J=JST(N,K)
C**** Calculate pressure by integrating from the top down
      PE = OPRESE(K,N)
      DO L=1,LMME(K,N)
         PEND(L,K) = PE + GRAV*MOE(K,N,L)*.5
         PE        = PE + GRAV*MOE(K,N,L)  ;  EndDo
C**** Calculate geopotential by integrating from the bottom up,
C**** also calculate the specific volume
      PHIE = -HOCEANE(K,N)*GRAV
      DO L=LMME(K,N),1,-1
      GUP = (G0ME(K,N,L) - 2*RRT12*GZME(K,N,L)) / (MOE(K,N,L)*DXYPO(J))
      GDN = (G0ME(K,N,L) + 2*RRT12*GZME(K,N,L)) / (MOE(K,N,L)*DXYPO(J))
      SUP = (S0ME(K,N,L) - 2*RRT12*SZME(K,N,L)) / (MOE(K,N,L)*DXYPO(J))
      SDN = (S0ME(K,N,L) + 2*RRT12*SZME(K,N,L)) / (MOE(K,N,L)*DXYPO(J))
      DP  = GRAV*MOE(K,N,L)
      PUP = PEND(L,K) - RRT12*DP
      PDN = PEND(L,K) + RRT12*DP
      VUP = VOLGSP(GUP,SUP,PUP)
      VDN = VOLGSP(GDN,SDN,PDN)
      PHI(L,K) = PHIE + (VUP*(.5-RRT12)+VDN*(.5+RRT12))*.5*DP
       DH(L,K) = MOE(K,N,L)*(VUP+VDN)*.5
C**** Calculate PHI at top of current layer (or bottom of next layer)
       PHIE = PHIE + (VDN+VUP)*.5*DP
      EndDo  !  End of Do L=
      EndDo  !  End of Do K=1,2
C**** Subtract Pressure Gradient Force from the mass flux
      DO L=1,LMST(N)
        MUST(L,N) = MUST(L,N) - .5*DTS*WIST(N)*
     *    (( DH(L,2)+ DH(L,1))*(PEND(L,2)-PEND(L,1)) +
     +     (PHI(L,2)-PHI(L,1))*(MOE(2,N,L)+MOE(1,N,L))) / DISTPG(N)
      EndDo  !  End of Do L=1,LMST(N)
      EndDo  !  End of Do N=1,NMST
C****
      RETURN
      END SUBROUTINE STPGF

      SUBROUTINE STADV(DTS)
C****
!@sum  STADV advects tracers and water mass through the straits
!@auth Gary Russell/Gavin Schmidt
!@ver  2009/05/26
C****
      USE CONSTANT, only : teeny
      USE OCEAN, only : dxypo,bydxypo
      Use STRAITS, Only: NMST,LMST,IST,JST,WIST,MUST,DISTPG,MMST,
     *                   G0MST,GXMST,GZMST, S0MST,SXMST,SZMST,
     *                   MOE, G0ME,GXME,GYME,GZME, S0ME,SXME,SYME,SZME,
     *                   KN2
      USE ODIAG, only : olnst,ln_mflx,ln_gflx,ln_sflx

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
      Use STRAITS, Only: TRMST,TXMST,TZMST, TRME,TXME,TYME,TZME
      Use ODIAG, Only: tlnst
#endif

      IMPLICIT NONE

      INTEGER I1,J1,I2,J2,N,L,ITR,K,KK,NN
      REAL*8, INTENT(IN) :: DTS
      REAL*8 MM1,MM2,AM
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif
C****
      DO N=1,NMST
      I1=IST(N,1)
      J1=JST(N,1)
      I2=IST(N,2)
      J2=JST(N,2)
      DO L=1,LMST(N)
      AM = DTS*MUST(L,N)
      MM1 = MOE(1,N,L)*DXYPO(J1)
      MM2 = MOE(2,N,L)*DXYPO(J2)
      CALL STADVT (N,L,AM,MM1,MM2,G0MST(L,N),GXMST(L,N),GZMST(L,N),
     *             G0ME,GXME,GYME,GZME,OLNST(L,N,LN_GFLX),.False.)
      CALL STADVT (N,L,AM,MM1,MM2,S0MST(L,N),SXMST(L,N),SZMST(L,N),
     *             S0ME,SXME,SYME,SZME,OLNST(L,N,LN_SFLX),.True.)
#ifdef TRACERS_OCEAN
      DO ITR = 1,tracerlist%getsize()
        entry=>tracerlist%at(itr)
        CALL STADVT (N,L,AM,MM1,MM2,TRMST(L,N,ITR),TXMST(L,N,ITR),
     *               TZMST(L,N,ITR),TRME(1,1,1,ITR),TXME(1,1,1,ITR),
     *               TYME(1,1,1,ITR),TZME(1,1,1,ITR),
     *               TLNST(L,N,1,ITR),entry%T_QLIMIT)  ;  EndDo
#endif /* def TRACERS_OCEAN */
      MOE(1,N,L) = MOE(1,N,L) - AM*byDXYPO(J1)
      MOE(2,N,L) = MOE(2,N,L) + AM*byDXYPO(J2)
        OLNST(L,N,LN_MFLX) = OLNST(L,N,LN_MFLX) + AM

C**** Limit heat gradients to 8000 (J/kg) = 2 (C) * SHCW (J/kg*C)
C**** to prevent problems in Red Sea area
      If (Abs(GXME(1,N,L)) > 8000*MOE(1,N,L)*DXYPO(J1))
     *    GXME(1,N,L) = Sign(8000*MOE(1,N,L)*DXYPO(J1),GXME(1,N,L))
      If (Abs(GYME(1,N,L)) > 8000*MOE(1,N,L)*DXYPO(J1))
     *    GYME(1,N,L) = Sign(8000*MOE(1,N,L)*DXYPO(J1),GYME(1,N,L))
      If (Abs(GZME(1,N,L)) > 8000*MOE(1,N,L)*DXYPO(J1))
     *    GZME(1,N,L) = Sign(8000*MOE(1,N,L)*DXYPO(J1),GZME(1,N,L))
      If (Abs(GXME(2,N,L)) > 8000*MOE(2,N,L)*DXYPO(J2))
     *    GXME(2,N,L) = Sign(8000*MOE(2,N,L)*DXYPO(J2),GXME(2,N,L))
      If (Abs(GYME(2,N,L)) > 8000*MOE(2,N,L)*DXYPO(J2))
     *    GYME(2,N,L) = Sign(8000*MOE(2,N,L)*DXYPO(J2),GYME(2,N,L))
      If (Abs(GZME(2,N,L)) > 8000*MOE(2,N,L)*DXYPO(J2))
     *    GZME(2,N,L) = Sign(8000*MOE(2,N,L)*DXYPO(J2),GZME(2,N,L))

      EndDo  !  End of Do L=1,LMST(N)

c Copy updated values to the arrays for a neighboring strait if one exists
      do k=1,2
        if(kn2(1,k,n).gt.0) then
          kk = kn2(1,k,n)
          nn = kn2(2,k,n)
          do l=1,lmst(n)
            moe(kk,nn,l) = moe(k,n,l)
            g0me(kk,nn,l) = g0me(k,n,l)
            gxme(kk,nn,l) = gxme(k,n,l)
            gyme(kk,nn,l) = gyme(k,n,l)
            gzme(kk,nn,l) = gzme(k,n,l)
            s0me(kk,nn,l) = s0me(k,n,l)
            sxme(kk,nn,l) = sxme(k,n,l)
            syme(kk,nn,l) = syme(k,n,l)
            szme(kk,nn,l) = szme(k,n,l)
          enddo
#ifdef TRACERS_OCEAN
          do itr = 1,tracerlist%getsize()
          do l=1,lmst(n)
            trme(kk,nn,l,itr) = trme(k,n,l,itr)
            txme(kk,nn,l,itr) = txme(k,n,l,itr)
            tyme(kk,nn,l,itr) = tyme(k,n,l,itr)
            tzme(kk,nn,l,itr) = tzme(k,n,l,itr)
          enddo
          enddo
#endif
        endif
      enddo

      EndDo  !  End of Do N=1,NMST
      RETURN
      END SUBROUTINE STADV

      SUBROUTINE STADVT(N,L,AM,MM1,MM2,RMST,RXST,RZST,RM,RX,RY,RZ,OLN
     *     ,QLIMIT)
C****
!@sum  STADVT advects tracers through the straits (improved calculation)
!@+    GOES BACK TO OLD CODING FOR STABILITY
!@auth Gary Russell/Gavin Schmidt
!@ver  2009/05/26
C****
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst,lmst,ist,jst,xst,yst,mmst
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: AM,MM1,MM2
      REAL*8, INTENT(INOUT) :: RMST,RXST,RZST
      Real*8,Intent(InOut),Dimension(2,NMST,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT) :: OLN
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8 A1,A2,FM1,FZ1,FM2,FZ2,RXY,X1,Y1,X2,Y2,RXold
      INTEGER I1,I2,J1,J2,N,L
C****
      I1=IST(N,1)
      J1=JST(N,1)
      X1=XST(N,1)
      Y1=YST(N,1)
      I2=IST(N,2)
      J2=JST(N,2)
      X2=XST(N,2)
      Y2=YST(N,2)
      IF(AM.LT.0.)  GO TO 200
C****
C**** Water flux is moving from grid box 1 to grid box 2  (AM > 0)
C****
       A1 = AM/MM1
c     FM1 = A1*(RM(1,N,L) + X1*RX(1,N,L) + Y1*RY(1,N,L))
      FM1 = A1*(RM(1,N,L) + (1-A1)*(X1*RX(1,N,L)+Y1*RY(1,N,L)))
      FZ1 = A1*RZ(1,N,L)
       A2 = AM/MMST(L,N)
      FM2 = A2*(RMST + (1-A2)*RXST)
      FZ2 = A2*RZST
C**** Calculate first moments of tracer mass for grid boxes
c      RX(1,N,L) = RX(1,N,L)*(1 - A1*(.5+1.5*X1*X1))
c      RY(1,N,L) = RY(1,N,L)*(1 - A1*(.5+1.5*Y1*Y1))
c      RXST      = RXST*(1-A2)**3 - 3*FM1 + 3*A2*RMST
c      RXold = RX(2,N,L)
c      RX(2,N,L) = RX(2,N,L) + (RX(2,N,L)*AM*(.5 - 1.5*X2*X2) +
c     *     3*X2*(FM2*MM2 - (RM(2,N,L)+Y2*RY(2,N,L))*AM)) / (MM2+AM)
c      RY(2,N,L) = RY(2,N,L) + (RY(2,N,L)*AM*(.5 - 1.5*Y2*Y2) +
c     *     3*Y2*(FM2*MM2 - (RM(2,N,L)+X2*RXold)*AM)) / (MM2+AM)
      RX(1,N,L) = RX(1,N,L)*(1-A1)*(1-A1*X1*X1)
      RY(1,N,L) = RY(1,N,L)*(1-A1)*(1-A1*Y1*Y1)
      RXST      = RXST*(1-2*A2) - FM1 + FM2
      RX(2,N,L) = RX(2,N,L) + X2*(FM2 - (RM(2,N,L)-X2*RX(2,N,L))
     $     *AM/MM2)
      RY(2,N,L) = RY(2,N,L) + Y2*(FM2 - (RM(2,N,L)-Y2*RY(2,N,L))
     $     *AM/MM2)
      GO TO 300
C****
C**** Water flux is moving from grid box 2 to grid box 1  (AM < 0)
C****
  200  A1 = AM/MMST(L,N)
      FM1 = A1*(RMST - (1+A1)*RXST)
      FZ1 = A1*RZST
       A2 = AM/MM2
c     FM2 = A2*(RM(2,N,L) + RX(2,N,L)*X2 + RY(2,N,L)*Y2)
      FM2 = A2*(RM(2,N,L) + (1+A2)*(X2*RX(2,N,L)+Y2*RY(2,N,L)))
      FZ2 = A2*RZ(2,N,L)
C**** Calculate first moments of tracer mass for grid boxes
c      RXold = RX(1,N,L)
c      RX(1,N,L) = RX(1,N,L) - (RX(1,N,L)*AM*(.5 - 1.5*X1*X1) +
c     *     3*X1*(FM1*MM1 - (RM(1,N,L)+Y1*RY(1,N,L))*AM)) / (MM1-AM)
c      RY(1,N,L) = RY(1,N,L) - (RY(1,N,L)*AM*(.5 - 1.5*Y1*Y1) +
c     *     3*Y1*(FM1*MM1 - (RM(1,N,L)+X1*RXold)*AM)) / (MM1-AM)
c      RXST      = RXST*(1+A1)**3 - 3*FM2 + 3*A1*RMST
c      RX(2,N,L) = RX(2,N,L)*(1 + A2*(.5+1.5*X2*X2))
c      RY(2,N,L) = RY(2,N,L)*(1 + A2*(.5+1.5*Y2*Y2))
      RX(1,N,L) = RX(1,N,L) -
     *            X1*(FM1 - (RM(1,N,L)-X1*RX(1,N,L))*AM/MM1)
      RY(1,N,L) = RY(1,N,L) -
     *            Y1*(FM1 - (RM(1,N,L)-Y1*RY(1,N,L))*AM/MM1)
      RXST      = RXST*(1+2*A1) + FM1 - FM2
      RX(2,N,L) = RX(2,N,L)*(1+A2)*(1+A2*X2*X2)
      RY(2,N,L) = RY(2,N,L)*(1+A2)*(1+A2*Y2*Y2)
C****
C**** Calculate new tracer mass, vertical moment of tracer mass,
C**** and horizontal moment of tracer mass for the straits
C****
  300 RM(1,N,L) = RM(1,N,L) - FM1
      RZ(1,N,L) = RZ(1,N,L) - FZ1
      RM(2,N,L) = RM(2,N,L) + FM2
      RZ(2,N,L) = RZ(2,N,L) + FZ2
      RMST      = RMST + (FM1-FM2)
      RZST      = RZST + (FZ1-FZ2)
       OLN      =  OLN + (FM1+FM2)
C****
      if ( QLIMIT ) then ! limit gradients at ends of strait
        RXY = abs(RX(1,N,L)) + abs(RY(1,N,L))
        if ( RXY > RM(1,N,L) ) then
          RX(1,N,L) = RX(1,N,L)*( RM(1,N,L)/(RXY + tiny(RXY)) )
          RY(1,N,L) = RY(1,N,L)*( RM(1,N,L)/(RXY + tiny(RXY)) )
        endif
        if ( abs(RZ(1,N,L)) > RM(1,N,L) )
     *       RZ(1,N,L) = sign(RM(1,N,L), RZ(1,N,L)+0d0)
        RXY = abs(RX(2,N,L)) + abs(RY(2,N,L))
        if ( RXY > RM(2,N,L) ) then
          RX(2,N,L) = RX(2,N,L)*( RM(2,N,L)/(RXY + tiny(RXY)) )
          RY(2,N,L) = RY(2,N,L)*( RM(2,N,L)/(RXY + tiny(RXY)) )
        endif
        if ( abs(RZ(2,N,L)) > RM(2,N,L) )
     *       RZ(2,N,L) = sign(RM(2,N,L), RZ(2,N,L)+0d0)
        if ( abs(RXST) > RMST ) RXST = sign(RMST, RXST)
        if ( abs(RZST) > RMST ) RZST = sign(RMST, RZST)
      end if
C****
      RETURN
      END SUBROUTINE STADVT

      SUBROUTINE STBDRA
!@sum  STBDRA exerts a bottom drag on the lowest layer and a side drag
!@+    on each layer of each strait. Also reduces along strait gradients
!@auth Gary Russell
C**** MUST = MUST*(1-x)  is approximated by  MUST = MUST/(1+x)
C****
C**** STRACB  MMST  Mass of sea water in strait (kg)
C****         MUST  Mass flow through strait (kg/s)
C****         WIST  Width of strait (m)
C****         DIST  Length of strait (m)
C****
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE OCEAN, only : lmo,dts
      USE STRAITS, only : nmst,lmst,must,mmst,wist,dist,sxmst,gxmst
#ifdef TRACERS_OCEAN
      Use STRAITS, Only: TXMST
#endif
      IMPLICIT NONE
      REAL*8, PARAMETER :: BDRAGX=1, SDRAGX=1d-1 
      REAL*8, PARAMETER :: REDUCE=1./(SECONDS_PER_DAY*20.)
      INTEGER N,L

      DO N=1,NMST
C**** Apply bottom drag
      L=LMST(N)
      MUST(L,N) = MUST(L,N) * MMST(L,N)**2 /
     *  (MMST(L,N)**2 + DTS*BDRAGX*ABS(MUST(L,N))*WIST(N)*DIST(N)**2)
      DO L=1,LMST(N)
C**** Apply side drag
        MUST(L,N) = MUST(L,N) * MMST(L,N)*WIST(N) /
     *       (MMST(L,N)*WIST(N) + DTS*SDRAGX*ABS(MUST(L,N))*DIST(N))
C**** Reduce cross strait tracer gradients (20 day restoring to zero)
        GXMST(L,N)=GXMST(L,N)*(1-REDUCE*DTS)
        SXMST(L,N)=SXMST(L,N)*(1-REDUCE*DTS)
#ifdef TRACERS_OCEAN
        TXMST(L,N,:)=TXMST(L,N,:)*(1-REDUCE*DTS)
#endif
      END DO
      END DO
C****
      RETURN
      END SUBROUTINE STBDRA

      SUBROUTINE init_STRAITS(iniOCEAN)
C****
!@sum  init_STRAITS initializes strait variables
!@auth Gary Russell/Gavin Schmidt
!@ver  2009/05/26
C****
      USE CONSTANT, only: twopi,radius
      Use OCEAN,   Only: IM,JM,LMO, ZE, FJEQ,DLON,DLAT, DXYPO,
     *                   MySparseComm_Type, HOCEAN, LMM, ZMID
      Use STRAITS, Only: NMST,LMST, IST,JST, XST,YST, WIST,DIST,DISTPG,
     *                   MMST,MUST, G0MST,GXMST,GZMST,S0MST,SXMST,SZMST,
!     *                   RSIST,RSIXST, MSIST,SSIST,HSIST,
     *                   MOE, G0ME,GXME,GYME,GZME, S0ME,SXME,SYME,SZME,
     *                   kn2,zst,
     *                   HOCEANe,LMMe
#ifdef OCN_GISS_TURB
     *                   ,otkest
#endif
      use domain_decomp_1d, only: am_i_root, getDomainBounds, 
     *                            getMpiCommunicator, broadcast
      USE OCEANR_DIM, only : grid=>ogrid
      Use SparseCommunicator_mod
#ifdef OCN_GISS_TURB
      USE GISSMIX_COM, Only: otke_init_max,emin,emax
#endif
      IMPLICIT NONE

      INTEGER I,J,L,N,I1,J1,I2,J2, ier, k,kk
      REAL*8 FLAT,DFLON,DFLAT,SLAT,DSLON,DSLAT,TLAT,DTLON
     *     ,DTLAT,G01,GZ1,S01,SZ1,G02,GZ2,S02,SZ2
      LOGICAL, INTENT(IN) :: iniOCEAN

      real*8 rLMMe(2,nmst), rLMM(im,grid%j_strt_halo:grid%j_stop_halo)

      !initialize communication data for gathering and scattering
      !straits ocean data
      integer  :: locLB(2),locUB(2),globLB(2),globUB(2)
      integer, allocatable :: points(:,:)

      allocate(points(2, 2*nmst), stat=ier)
c      points(1,:) = RESHAPE( ist(:,:) , (/2*NMST/) )
c      points(2,:) = RESHAPE( jst(:,:) , (/2*NMST/) )
c The two ends of each strait must be adjacent in the list of points!
      do n=1,nmst
        points(1,2*n-1) = ist(n,1)
        points(2,2*n-1) = jst(n,1)
        points(1,2*n  ) = ist(n,2)
        points(2,2*n  ) = jst(n,2)
      enddo

      call getDomainBounds(grid, J_STRT_HALO = locLB(2), 
     &     J_STOP_HALO = locUB(2))

      locLB(1)=1; locUB(1)=grid%IM_WORLD;
      !locLB(2)=jmin; locUB(2)=jmax;

      !bounds global
      globLB(1)=1; globUB(1)=grid%IM_WORLD;
      globLB(2)=1; globUB(2)=grid%JM_WORLD;

      mySparseComm_type = SparseCommunicator(points,
     &                  locLB  , locUB  , globLB , globUB,
     &                  comm = getMpiCommunicator(grid))
      deallocate(points)

C****
C**** Check whether any straits share I,J endpoints, assuming that
C**** straits sharing I,J endpoints are adjacent in the list and
C**** that only 2 straits can share any I,J endpoint.
C****
      if(am_i_root()) then
        kn2 = 0
        do n=1,nmst-1
          do k=1,2; do kk=1,2
            if(ist(n,k).eq.ist(n+1,kk) .and.
     &         jst(n,k).eq.jst(n+1,kk)) then
              kn2(1:2,k ,n  ) = (/ kk, n+1 /)
              kn2(1:2,kk,n+1) = (/ k , n   /)
            endif
          enddo; enddo
        enddo
      endif

      call gather_straits_pairs(HOCEAN,HOCEANe, 1)
      rLMM(:,:) = LMM(:,:)
      call gather_straits_pairs(rLMM,rLMMe, 1)
      if(am_i_root()) LMMe = rLMMe
      call broadcast(grid,LMMe)

C****
C**** Calculate distance of strait, distance between centers of
C**** ocean grid boxes through strait for pressure gradient force,
C**** mass of water in the strait, depth of strait
C****
cc    DLON = TWOPI/IM
cc    DLAT = NINT(180d0/(JM-1))*TWOPI/360d0
cc    FJEQ = (JM+1)/2d0
      DO N=1,NMST
      IF(ZST(N) > 0D0) THEN ! if depth was specified in input file, use it
        DO L=1,LMO-1
          IF(ZMID(L+1).GT.ZST(N)) EXIT
        ENDDO
        LMST(N) = L
      ENDIF
      LMST(N) = MIN(LMST(N),MINVAL(LMMe(:,N)))
      ZST(N) = ZE(LMST(N))
      DSLON = (IST(N,2)+.5*XST(N,2)-IST(N,1)-.5*XST(N,1))*DLON
      DSLAT = (JST(N,2)+.5*YST(N,2)-JST(N,1)-.5*YST(N,1))*DLAT
       SLAT = (JST(N,2)+.5*YST(N,2)+JST(N,1)+.5*YST(N,1)-2*FJEQ)
     *     *DLAT/2
      DIST(N) = RADIUS*SQRT((DSLON*COS(SLAT))**2 + DSLAT*DSLAT)
      DFLON = .5*XST(N,1)*DLON
      DFLAT = .5*YST(N,1)*DLAT
       FLAT = (JST(N,1) + .25*YST(N,1)-FJEQ)*DLAT
      DTLON = .5*XST(N,2)*DLON
      DTLAT = .5*YST(N,2)*DLAT
       TLAT = (JST(N,2) + .25*YST(N,2)-FJEQ)*DLAT
      DISTPG(N) = DIST(N) +
     +   RADIUS*SQRT((DFLON*COS(FLAT))**2 + DFLAT*DFLAT) +
     +   RADIUS*SQRT((DTLON*COS(TLAT))**2 + DTLAT*DTLAT)
      DO L=1,LMST(N)
        MMST(L,N) = DIST(N)*WIST(N)*(ZE(L)-ZE(L-1))/9.7d-4
      END DO
      DO L=LMST(N)+1,LMO
        MMST(L,N) = 0.
      END DO
      END DO
C****

      IF(iniOCEAN) THEN
C****
C**** Initialize the mass flux, potential heat, and salinity in
C**** each strait
C****
      call gather_ocean_straits ! get values at endpoints of straits

      if(am_I_root()) then
      DO N=1,NMST
      I1=IST(N,1)
      J1=JST(N,1)
      I2=IST(N,2)
      J2=JST(N,2)
      DO L=1,LMST(N)
      MUST(L,N) = 0.
      G01 = (G0ME(1,N,L)+XST(N,1)*GXME(1,N,L)+YST(N,1)*GYME(1,N,L)) /
     /       (MOE(1,N,L)*DXYPO(J1))
      GZ1 =  GZME(1,N,L) / (MOE(1,N,L)*DXYPO(J1))
      G02 = (G0ME(2,N,L)+XST(N,2)*GXME(2,N,L)+YST(N,2)*GYME(2,N,L)) /
     /       (MOE(2,N,L)*DXYPO(J2))
      GZ2 =  GZME(2,N,L) / (MOE(2,N,L)*DXYPO(J2))
      G0MST(L,N) = (G01+G02)*MMST(L,N)*.5
      GXMST(L,N) = (G02-G01)*MMST(L,N)*.5
      GZMST(L,N) = (GZ1+GZ2)*MMST(L,N)*.5
#ifdef OCN_GISS_TURB
      otkest(L,N) = min(max(otke_init_max/(float(l)**2),emin),emax)
#endif
      S01 = (S0ME(1,N,L)+XST(N,1)*SXME(1,N,L)+YST(N,1)*SYME(1,N,L)) /
     /       (MOE(1,N,L)*DXYPO(J1))
      SZ1 =  SZME(1,N,L) / (MOE(1,N,L)*DXYPO(J1))
      S02 = (S0ME(2,N,L)+XST(N,2)*SXME(2,N,L)+YST(N,2)*SYME(2,N,L)) /
     /       (MOE(2,N,L)*DXYPO(J2))
      SZ2 =  SZME(2,N,L) / (MOE(2,N,L)*DXYPO(J2))
      S0MST(L,N) = (S01+S02)*MMST(L,N)*.5
      SXMST(L,N) = (S02-S01)*MMST(L,N)*.5
      SZMST(L,N) = (SZ1+SZ2)*MMST(L,N)*.5
      EndDo  !  End of Do L+1,LMST(N)
      DO L=LMST(N)+1,LMO
       MUST(L,N) = 0.
      G0MST(L,N) = 0.
      GXMST(L,N) = 0.
      GZMST(L,N) = 0.
#ifdef OCN_GISS_TURB
      otkest(L,N) = 0.
#endif
      S0MST(L,N) = 0.
      SXMST(L,N) = 0.
      SZMST(L,N) = 0.
      END DO
      END DO
C**** Initialize sea ice in straits
!      RSIST=0.
!      RSIXST=0.
!      MSIST=0.
!      HSIST=0.
!      SSIST=0.
      end if                    ! root-process only

      call bcast_straits(.true.) ! skip tracers

      END IF
C****
      RETURN
      END SUBROUTINE init_STRAITS

c      SUBROUTINE STADVI
c!@sum  STADVI advects sea ice through the straits
c!@auth Gary Russell/Gavin Schmidt
c      USE OCEAN, only : im,jm,dts,dxypo,bydxypo
c      USE SEAICE, only : xsi,ace1i
c      USE SEAICE_COM, only : lmi,rsi,hsi,msi,snowi,ssi
c      USE STRAITS
c      USE ICEDYN_COM, only : rsix,rsiy
c      USE ODIAG, only : olnst,ln_icfl
c#ifdef TRACERS_WATER
c      Use SEAICE_COM, Only: TRSI,NTM
c#endif
c      IMPLICIT NONE
c
c      REAL*8, SAVE ::  WDST(NMST),BYWDST(NMST)
c      REAL*8 MHS(NTRICE,IM,JM),MHSIST(NTRICE,NMST),MHSL(NTRICE)
c      REAL*8 USTDT,DRSIST,ASI,ASIST,DMHSI,ALPHA,RMEAN
c      INTEGER I,J,K,L,N
c      INTEGER, SAVE :: IFIRST=1

c!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
c#ifndef TRACERS_WATER
c      INTEGER, PARAMETER :: NTRICE=2+2*LMI
c#else
c      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
c      INTEGER ITR
c#endif
c
c      IF(IFIRST.EQ.1) THEN
c        IFIRST = 0
c        DO N=1,NMST
c          WDST(N)    = WIST(N)*DIST(N)
c          BYWDST(N)  = 1 / WDST(N)
c        END DO
c      END IF
cC**** set up local MHS array to contain all advected quantities
cC**** MHS(1:2) = MASS, MHS(3:6) = HEAT, MHS(7:10)=SALT
c      MHS(1,:,:) = ACE1I + SNOWI
c      MHS(2,:,:) = MSI
c      MHSIST(1:2,:) = MSIST(1:2,:)
c      DO L=1,LMI
c        MHS(L+2,:,:) = HSI(L,:,:)
c        MHS(L+2+LMI,:,:) = SSI(L,:,:)
c        MHSIST(L+2,:) = HSIST(L,:)
c        MHSIST(L+2+LMI,:) = SSIST(L,:)
c#ifdef TRACERS_WATER
cC**** add tracers to advected arrays
c        DO ITR=1,NTM
c          MHS(L+2+(1+ITR)*LMI,:,:)=TRSI(ITR,L,:,:)
c          MHSIST(L+2+(1+ITR)*LMI,:)=TRSIST(ITR,L,:)
c        END DO
c#endif
c      END DO
cC****
cC**** Loop over all straits
cC****
c      DO 900 N=1,NMST
c      IF(MUST(1,N).GT.0.) GO TO 500
c      IF(RSIST(N).LE.0.)  GO TO 300
cC****
cC**** MUST < 0: sea ice may be leaving western end of strait
cC****
c      USTDT = USIFAC*DTS*MUST(1,N)*DIST(N)/MMST(1,N)
c      IF(RSIXST(N).LE.0.)  GO TO 120
c      IF(USTDT+2*RSIXST(N) .LT. 0.)  GO TO 110
cC**** RSIXST > 0 and no sea ice reaches western end of strait
c      RSIXST(N) = RSIXST(N) + USTDT
c      GO TO 300
cC**** RSIXST > 0 but some sea ice flows into western ocean box
c  110 DRSIST = -RSIST(N)*(USTDT+2*RSIXST(N))/(DIST(N)-2*RSIXST(N))
c      RSIST(N) =  RSIST(N) - DRSIST
c      RSIXST(N) = .5*USTDT
c      GO TO 200
c  120 IF(USTDT+2*RSIXST(N)+DIST(N) .LE. 0.)  GO TO 130
cC**** RSIXST < 0 and some sea ice flows into western ocean box
c      DRSIST   = -RSIST(N)*USTDT / (DIST(N)+2*RSIXST(N))
c      RSIST(N) =  RSIST(N) - DRSIST
c      RSIXST(N) =  RSIXST(N) + .5*USTDT
c      GO TO 200
cC**** RSIXST < 0 and all sea ice in strait flows into western ocean box
c  130 DRSIST   = RSIST(N)
c      RSIST(N) = 0.
c      RSIXST(N) = 0.
cC****
cC**** Western ocean box receives sea ice from strait
cC****
c  200 I = IST(N,1)
c      J = JST(N,1)
c      ASI   = RSI(I,J)*DXYPO(J)
c      ASIST = DRSIST*WDST(N)
c      DO K=1,NTRICE
c        MHSL(K) = ASI*MHS(K,I,J) + ASIST*MHSIST(K,N)
c      END DO
c        OLNST(2,N,LN_ICFL) = OLNST(2,N,LN_ICFL) - ASIST*(MHSL(1)
c     *     +MHSL(2))
c      IF(ASI+ASIST.GT.DXYPO(J))  GO TO 210
c      RSI(I,J)   = RSI(I,J)  + ASIST*BYDXYPO(J)
c      RSIY(I,J)  = RSIY(I,J) + ASIST*BYDXYPO(J)*2*YST(N,1)
c      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)
c      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)
c      IF(RSI(I,J)-RSIY(I,J).GT.1) RSIY(I,J) =  RSI(I,J) - 1
c      IF(RSI(I,J)+RSIY(I,J).GT.1) RSIY(I,J) =  1 - RSI(I,J)
c      IF(RSI(I,J)-RSIX(I,J).gt.1) RSIX(I,J) =  RSI(I,J)-1
c      IF(RSI(I,J)+RSIX(I,J).gt.1) RSIX(I,J) =  1-RSI(I,J)
c      DO K=1,NTRICE
c        MHS(K,I,J) = MHSL(K) / (ASI+ASIST)
c      END DO
c      GO TO 300
cC**** Sea ice crunches into itself and completely covers grid box
c  210 RSI(I,J)  = 1
c      RSIY(I,J) = 0.
c      MHS(1,I,J) = MHSL(1) / (ASI+ASIST)
c      MHS(2,I,J) = (MHSL(1)+MHSL(2))*BYDXYPO(J) - MHS(1,I,J)
c      DO K=1,(NTRICE-2)/LMI
c        MHS(3+LMI*(K-1),I,J) = MHSL(3+LMI*(K-1)) / (ASI+ASIST)
c        MHS(4+LMI*(K-1),I,J) = MHSL(4+LMI*(K-1)) / (ASI+ASIST)
c        DMHSI = (MHSL(3+LMI*(K-1))+MHSL(4+LMI*(K-1))+MHSL(5+LMI*(K-1))
c     *       +MHSL(6+LMI*(K-1)))*(BYDXYPO(J) -1/(ASI+ASIST))
c        MHS(5+LMI*(K-1),I,J) = MHSL(5+LMI*(K-1)) / (ASI+ASIST) + XSI(3)
c     *       *DMHSI
c        MHS(6+LMI*(K-1),I,J) = MHSL(6+LMI*(K-1)) / (ASI+ASIST) + XSI(4)
c     *       *DMHSI
c      END DO
cC****
cC**** Eastern ocean box may send sea ice into strait
cC****
c  300 I = IST(N,2)
c      J = JST(N,2)
c      IF(RSI(I,J).LE.0.)  GO TO 900
c      USTDT = USIFAC*DTS*MUST(1,N)*DIST(N)/MMST(1,N)
c      ALPHA = -USTDT*WIST(N)*BYDXYPO(J)
c      RMEAN = RSI(I,J) + YST(N,2)*RSIY(I,J)*(1-ALPHA)
c      ASI   = ALPHA*RMEAN*DXYPO(J)
c      RSI(I,J)  = RSI(I,J) - ALPHA*RMEAN
c      RSIY(I,J) = RSIY(I,J)*(1-ALPHA)*(1-ALPHA)
c      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)
c      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)
c      IF(RSI(I,J)-RSIY(I,J).GT.1) RSIY(I,J) =  RSI(I,J) - 1
c      IF(RSI(I,J)+RSIY(I,J).GT.1) RSIY(I,J) =  1 - RSI(I,J)
c      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =  RSI(I,J)
c      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) = -RSI(I,J)
cC****
cC**** Sea ice is entering strait from eastern ocean box
cC****
c  400 DRSIST   = ASI*BYWDST(N)
c      RSIXST(N) = (RSIST(N)*(MHSIST(1,N)+MHSIST(2,N))*RSIXST(N) +
c     *            DRSIST*(MHS(1,I,J)+MHS(2,I,J))*.5*(DIST(N)+USTDT))/
c     *           (RSIST(N)*(MHSIST(1,N)+MHSIST(2,N)) +
c     *             DRSIST*(MHS(1,I,J)+MHS(2,I,J)))
c      DO K=1,NTRICE
c        MHSIST(K,N) = (RSIST(N)*MHSIST(K,N) + DRSIST*MHS(K,I,J)) /
c     *       (RSIST(N) + DRSIST)
c      END DO
c      RSIST(N) = RSIST(N) + DRSIST
c         OLNST(1,N,LN_ICFL) = OLNST(2,N,LN_ICFL) -
c     *     ASI*(MHS(1,I,J)+MHS(2,I,J))
c      GO TO 900
cC****
c  500 IF(RSIST(N).LE.0.)  GO TO 700
cC****
cC**** MUST > 0: sea ice may be leaving eastern end of strait
cC****
c      USTDT = USIFAC*DTS*MUST(1,N)*DIST(N)/MMST(1,N)
c      IF(RSIXST(N).GE.0.)  GO TO 520
c      IF(USTDT+2*RSIXST(N) .GT. 0.)  GO TO 510
cC**** RSIXST < 0 and no sea ice reaches eastern end of strait
c      RSIXST(N) = RSIXST(N) + USTDT
c      GO TO 700
cC**** RSIXST < 0 but some sea ice flows into eastern ocean box
c  510 DRSIST = RSIST(N)*(USTDT+2*RSIXST(N)) / (DIST(N)+2*RSIXST(N))
c      RSIST(N) = RSIST(N) - DRSIST
c      RSIXST(N) = .5*USTDT
c      GO TO 600
c  520 IF(USTDT+2*RSIXST(N)-DIST(N) .GE. 0.)  GO TO 530
cC**** RSIXST > 0 and some sea ice flows into eastern ocean box
c      DRSIST   = RSIST(N)*USTDT / (DIST(N)-2*RSIXST(N))
c      RSIST(N) = RSIST(N) - DRSIST
c      RSIXST(N) = RSIXST(N) + .5*USTDT
c      GO TO 600
cC**** RSIXST > 0 and all sea ice in strait flows into eastern ocean box
c  530 DRSIST   = RSIST(N)
c      RSIST(N) = 0.
c      RSIXST(N) = 0.
cC****
cC**** Eastern ocean box receives sea ice from strait
cC****
c  600 I = IST(N,2)
c      J = JST(N,2)
c      ASI   = RSI(I,J)*DXYPO(J)
c      ASIST = DRSIST*WDST(N)
c      DO K=1,NTRICE
c        MHSL(K) = ASI*MHS(K,I,J) + ASIST*MHSIST(K,N)
c      END DO
c         OLNST(2,N,LN_ICFL) = OLNST(2,N,LN_ICFL) +
c     *     ASIST*(MHSIST(1,N)+MHSIST(2,N))
c      IF(ASI+ASIST.GT.DXYPO(J))  GO TO 610
c      RSI(I,J)   = RSI(I,J)  + ASIST*BYDXYPO(J)
c      RSIY(I,J)  = RSIY(I,J) + ASIST*BYDXYPO(J)*2*YST(N,2)
c      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)
c      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)
c      IF(RSI(I,J)-RSIY(I,J).GT.1) RSIY(I,J) =  RSI(I,J) - 1
c      IF(RSI(I,J)+RSIY(I,J).GT.1) RSIY(I,J) =  1 - RSI(I,J)
c      IF(RSI(I,J)-RSIX(I,J).gt.1) RSIX(I,J) =  RSI(I,J) - 1
c      IF(RSI(I,J)+RSIX(I,J).gt.1) RSIX(I,J) =  1 - RSI(I,J)
c      DO K=1,NTRICE
c        MHS(K,I,J) = MHSL(K) / (ASI+ASIST)
c      END DO
c      GO TO 700
cC**** Sea ice crunches into itself and completely covers grid box
c  610 RSI(I,J)  = 1
c      RSIY(I,J) = 0.
c      MHS(1,I,J) = MHSL(1) / (ASI+ASIST)
c      MHS(2,I,J) = (MHSL(1)+MHSL(2))*BYDXYPO(J) - MHS(1,I,J)
c      DO K=1,(NTRICE-2)/LMI
c        MHS(3+LMI*(K-1),I,J) = MHSL(3+LMI*(K-1)) / (ASI+ASIST)
c        MHS(4+LMI*(K-1),I,J) = MHSL(4+LMI*(K-1)) / (ASI+ASIST)
c        DMHSI = (MHSL(3+LMI*(K-1))+MHSL(4+LMI*(K-1))+MHSL(5+LMI*(K-1))
c     *       +MHSL(6+LMI*(K-1)))*(BYDXYPO(J) -1/(ASI+ASIST))
c        MHS(5+LMI*(K-1),I,J) = MHSL(5+LMI*(K-1)) / (ASI+ASIST) + XSI(3)
c     *       *DMHSI
c        MHS(6+LMI*(K-1),I,J) = MHSL(6+LMI*(K-1)) / (ASI+ASIST) + XSI(4)
c     *       *DMHSI
c      END DO
cC****
cC**** Western ocean box may send sea ice into strait
cC****
c  700 I = IST(N,1)
c      J = JST(N,1)
c      IF(RSI(I,J).LE.0.)  GO TO 900
c      USTDT = USIFAC*DTS*MUST(1,N)*DIST(N)/MMST(1,N)
c      ALPHA = USTDT*WIST(N)*BYDXYPO(J)
c      RMEAN = RSI(I,J) + YST(N,1)*RSIY(I,J)*(1-ALPHA)
c      ASI   = ALPHA*RMEAN*DXYPO(J)
c      RSI(I,J)  = RSI(I,J) - ALPHA*RMEAN
c      RSIY(I,J) = RSIY(I,J)*(1-ALPHA)*(1-ALPHA)
c      IF(RSI(I,J)-RSIY(I,J).LT.0.)  RSIY(I,J) =  RSI(I,J)
c      IF(RSI(I,J)+RSIY(I,J).LT.0.)  RSIY(I,J) = -RSI(I,J)
c      IF(RSI(I,J)-RSIY(I,J).GT.1) RSIY(I,J) =  RSI(I,J) - 1
c      IF(RSI(I,J)+RSIY(I,J).GT.1) RSIY(I,J) =  1 - RSI(I,J)
c      IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =  RSI(I,J)
c      IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) = -RSI(I,J)
cC****
cC**** Sea ice is entering strait from western ocean box
cC****
c  800 DRSIST   = ASI*BYWDST(N)
c      RSIXST(N) = (RSIST(N)*(MHSIST(1,N)+MHSIST(2,N))*RSIXST(N) +
c     *            DRSIST*(MHS(1,I,J)+MHS(2,I,J))*.5*(USTDT-DIST(N)))/
c     *           (RSIST(N)*(MHSIST(1,N)+MHSIST(2,N)) +
c     *              DRSIST*(MHS(1,I,J)+MHS(2,I,J)))
c      DO K=1,NTRICE
c        MHSIST(K,N) = (RSIST(N)*MHSIST(K,N) + DRSIST*MHS(K,I,J)) /
c     *       (RSIST(N) + DRSIST)
c      END DO
c      RSIST(N)   =  RSIST(N) + DRSIST
c         OLNST(1,N,LN_ICFL) = OLNST(2,N,LN_ICFL) +
c     *     ASI*(MHS(1,I,J)+MHS(2,I,J))
cC****
c  900 CONTINUE
cC**** set global variables from local array
c      SNOWI(:,:)= MAX(0d0,MHS(1,:,:) - ACE1I)
c      MSI(:,:)  = MHS(2,:,:)
c      MSIST(1:2,:) = MHSIST(1:2,:)
c      DO L=1,LMI
c        HSI(L,:,:) = MHS(L+2,:,:)
c        SSI(L,:,:) = MHS(L+2+LMI,:,:)
c        HSIST(L,:) = MHSIST(L+2,:)
c        SSIST(L,:) = MHSIST(L+2+LMI,:)
c#ifdef TRACERS_WATER
cC**** add tracers to advected arrays
c        DO ITR=1,NTM
c          TRSI(ITR,L,:,:)=MHS(L+2+(1+ITR)*LMI,:,:)
c          TRSIST(ITR,L,:)=MHSIST(L+2+(1+ITR)*LMI,:)
c        END DO
c#endif
c      END DO
cC****
c      RETURN
c      END

      SUBROUTINE CHECKOST(SUBR)
!@sum  CHECKOST Checks whether Straits are reasonable
!@auth Original Development Team
      USE MODEL_COM, only : qcheck,dtsrc
      USE SEAICE, only : xsi,lmi
      USE STRAITS
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#endif
      IMPLICIT NONE

      REAL*8 relerr,errmax
      INTEGER L,n,ns,nmax,lmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif

C**** Check for NaN/INF in ocean data
      IF (QCHECK) THEN
      CALL CHECK3(G0MST,LMO,NMST,1,SUBR,'g0mst')
      CALL CHECK3(GXMST,LMO,NMST,1,SUBR,'gxmst')
      CALL CHECK3(GZMST,LMO,NMST,1,SUBR,'gzmst')
#ifdef OCN_GISS_TURB
      CALL CHECK3(otkest,LMO,NMST,1,SUBR,'otkest')
#endif
      CALL CHECK3(S0MST,LMO,NMST,1,SUBR,'s0mst')
      CALL CHECK3(SXMST,LMO,NMST,1,SUBR,'sxmst')
      CALL CHECK3(SZMST,LMO,NMST,1,SUBR,'szmst')
      CALL CHECK3(MUST ,LMO,NMST,1,SUBR,'must')
!      CALL CHECK3(MSIST,2,NMST,1,SUBR,'msist')
!      CALL CHECK3(SSIST,LMI,NMST,1,SUBR,'hsist')
!      CALL CHECK3(HSIST,LMI,NMST,1,SUBR,'ssist')
!      CALL CHECK3(RSIST,NMST,1,1,SUBR,'rsist')
!      CALL CHECK3(RSIXST,NMST,1,1,SUBR,'rsxst')
#ifdef TRACERS_OCEAN
      CALL CHECK3(TRMST,LMO,NMST,tracerlist%getsize(),SUBR,'trmst')
      CALL CHECK3(TXMST,LMO,NMST,tracerlist%getsize(),SUBR,'txmst')
      CALL CHECK3(TZMST,LMO,NMST,tracerlist%getsize(),SUBR,'tzmst')
c      CALL CHECK3(TRSIST,NTM_ATM,LMI,NMST,SUBR,'trist')
#endif /* def TRACERS_OCEAN */

      DO N=1,NMST
        DO L=1,LMST(N)
          IF (ABS(DTSrc*MUST(L,N)/MMST(L,N)).gt.0.2) THEN
            print*,"After ",SUBR," MASS FLUX > 20% ",L,N,DTSrc*MUST(L,N)
     *           /MMST(L,N)
          END IF
        END DO
      END DO

#ifdef TRACERS_OCEAN
C**** Check conservation of water tracers in straits
      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        if (entry%trname.eq.'Water') then
          errmax = 0. ; nmax=1 ; lmax=1
          do ns=1,nmst
          do l=1,lmst(ns)
            relerr=max(
     *           abs(trmst(l,ns,n)-mmst(l,ns)+s0mst(l,ns)),
     *           abs(txmst(l,ns,n)+sxmst(l,ns)),
     *           abs(tzmst(l,ns,n)+szmst(l,ns)))/
     *           (mmst(l,ns)-s0mst(l,ns))
            if (relerr.gt.errmax) then
              nmax=ns ; lmax=l ; errmax=relerr
            end if
          end do
          end do
          print*,"Relative error in straits fresh water mass after "
     *         ,subr,":",nmax,lmax,errmax,trmst(lmax,nmax,n),mmst(lmax
     *         ,nmax)-s0mst(lmax,nmax),txmst(lmax,nmax,n),-sxmst(lmax
     *         ,nmax),tzmst(lmax,nmax,n),-szmst(lmax,nmax)

C**** now straits ice (obsolete)
c          errmax = 0. ; nmax=1 ; lmax=1
c          do ns=1,nmst
c          do l=1,lmst(ns)
c              relerr=max(
c     *           abs(trsist(n,1,ns)-msist(1,ns)*xsi(1)+ssist(1,ns))
c     *           /(msist(1,ns)*xsi(1)-ssist(1,ns)),abs(trsist(n,2,ns)
c     *           -msist(1,ns)*xsi(2)+ssist(2,ns))/(msist(1,ns)*xsi(2)
c     *           -ssist(2,ns)),abs(trsist(n,3,ns)-msist(2,ns)*xsi(3)
c     *           +ssist(3,ns))/(msist(2,ns)*xsi(3)-ssist(3,ns))
c     *           ,abs(trsist(n,4,ns)-msist(2,ns)*xsi(4)+ssist(4,ns))
c     *           /(msist(2,ns)*xsi(4)-ssist(4,ns)))
c            if (relerr.gt.errmax) then
c              nmax=ns ; lmax=l ; errmax=relerr
c            end if
c          end do
c          end do
c          print*,"Relative error in straits ice mass after ",subr
c     *         ,":",nmax,lmax,errmax,trsist(n,1:4,nmax),msist(1,nmax)
c     *         *xsi(1:2)-ssist(1:2,nmax),msist(2,nmax)*xsi(3:4)
c     *         -ssist(3:4,nmax)
        end if
      end do
#endif /* def TRACERS_OCEAN */

      END IF
C****
      END SUBROUTINE CHECKOST

      SUBROUTINE BCAST_straits (skip_tracers)
      USE STRAITS
      use domain_decomp_1d, only : broadcast
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE
      logical, intent(in) :: skip_tracers

      CALL broadcast(grid, MUST )
      CALL broadcast(grid, G0MST)
      CALL broadcast(grid, GXMST)
      CALL broadcast(grid, GZMST)
#ifdef OCN_GISS_TURB
      CALL broadcast(grid, otkest)
#endif
      CALL broadcast(grid, S0MST)
      CALL broadcast(grid, SXMST)
      CALL broadcast(grid, SZMST)
!      CALL broadcast(grid, RSIST)
!      CALL broadcast(grid, RSIXST)
!      CALL broadcast(grid, MSIST)
!      CALL broadcast(grid, HSIST)
!      CALL broadcast(grid, SSIST)

      if(skip_tracers) return
!#ifdef TRACERS_WATER
!      CALL broadcast(grid, TRSIST)
!#endif
#ifdef TRACERS_OCEAN
      CALL broadcast(grid, TRMST)
      CALL broadcast(grid, TXMST)
      CALL broadcast(grid, TZMST)
#endif /* def TRACERS_OCEAN */
      return
      end SUBROUTINE BCAST_straits

      subroutine gather_ocean_straits
      USE STRAITS
      USE OCEAN
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist
#endif
      implicit none

      call gather_straits_pairs(MO,MOe,lmo)

      call gather_straits_pairs(G0M, G0Me, lmo)
      call gather_straits_pairs(GXMO,GXMe, lmo)
      call gather_straits_pairs(GYMO,GYMe, lmo)
      call gather_straits_pairs(GZMO,GZMe, lmo)

      call gather_straits_pairs(S0M, S0Me, lmo)
      call gather_straits_pairs(SXMO,SXMe, lmo)
      call gather_straits_pairs(SYMO,SYMe, lmo)
      call gather_straits_pairs(SZMO,SZMe, lmo)

      call gather_straits_pairs(OPRESS,OPRESe, 1)

#ifdef TRACERS_OCEAN
      call gather_straits_pairs(TRMO,TRMe, lmo*tracerlist%getsize())
      call gather_straits_pairs(TXMO,TXMe, lmo*tracerlist%getsize())
      call gather_straits_pairs(TYMO,TYMe, lmo*tracerlist%getsize())
      call gather_straits_pairs(TZMO,TZMe, lmo*tracerlist%getsize())
#endif

      return
      end subroutine gather_ocean_straits

      subroutine gather_straits_pairs(X3D,X_pairs,nl)
c Mini-routine to present straits arrays to sparsecommunicator
c dimensioned as (2*nmst,:) rather than (2,nmst,:)
      USE STRAITS, only : nmst
      USE OCEAN, only : im, mySparseComm_type
      use SparseCommunicator_mod, only : gatherPTS
      use OCEANR_DIM, only : grid=>ogrid
      implicit none
      integer :: nl
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,nl) :: X3D
      real*8, dimension(2*nmst,nl) :: X_pairs
      call gatherPTS(mySparseComm_type, X3D, X_pairs)
      return
      end subroutine gather_straits_pairs

      subroutine scatter_ocean_straits
      USE STRAITS
      USE OCEAN
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist
#endif
      implicit none

      call scatter_straits_pairs(MOe,MO,lmo)

      call scatter_straits_pairs(G0Me, G0M, lmo)
      call scatter_straits_pairs(GXMe,GXMO, lmo)
      call scatter_straits_pairs(GYMe,GYMO, lmo)
      call scatter_straits_pairs(GZMe,GZMO, lmo)

      call scatter_straits_pairs(S0Me, S0M, lmo)
      call scatter_straits_pairs(SXMe,SXMO, lmo)
      call scatter_straits_pairs(SYMe,SYMO, lmo)
      call scatter_straits_pairs(SZMe,SZMO, lmo)

c is this needed? opress not updated by straits routines
c      call scatter_straits_pairs(OPRESe,OPRESS, 1)

#ifdef TRACERS_OCEAN
      call scatter_straits_pairs(TRMe,TRMO, lmo*tracerlist%getsize())
      call scatter_straits_pairs(TXMe,TXMO, lmo*tracerlist%getsize())
      call scatter_straits_pairs(TYMe,TYMO, lmo*tracerlist%getsize())
      call scatter_straits_pairs(TZMe,TZMO, lmo*tracerlist%getsize())
#endif

      return

      end subroutine scatter_ocean_straits

      subroutine scatter_straits_pairs(X_pairs,X3D,nl)
c Mini-routine to present straits arrays to sparsecommunicator
c dimensioned as (2*nmst,:) rather than (2,nmst,:)
      USE STRAITS, only : nmst
      USE OCEAN, only : im, mySparseComm_type
      use SparseCommunicator_mod, only : scatterPTS
      use OCEANR_DIM, only : grid=>ogrid
      implicit none
      integer :: nl
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,nl) :: X3D
      real*8, dimension(2*nmst,nl) :: X_pairs
      call scatterPTS(mySparseComm_type, X_pairs, X3D)
      return
      end subroutine scatter_straits_pairs
