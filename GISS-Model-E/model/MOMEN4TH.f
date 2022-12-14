      module MOMENTS
      implicit none
      private

      public ADVECV, moment_enq_order

      contains
cddd      subroutine init_MOM
cddd!@sum  init_MOM sets an order dependent coefficient for AVRX
cddd      USE constant, only : byrt2
cddd      USE DYNAMICS, only : XAVRX
cddd      XAVRX = byrt2 ! for 4th order scheme, 1 for 2nd order scheme
cddd      CALL AVRX0
cddd      RETURN
cddd      END subroutine init_MOM

      subroutine moment_enq_order(order)
!@sum moment_enq_order returns order of the scheme
      implicit none
      integer, intent(out) :: order
      order = 4
      return
      end subroutine moment_enq_order


      SUBROUTINE ADVECV (PA,UT,VT,PB,U,V,P,DT1)
!@sum  ADVECV Advects momentum (incl. coriolis) using 4-th order scheme
!@auth Original development team
      USE constant, only : by3,by6,by12
      USE RESOLUTION, only : im,jm,lm,ls1,psfmpt
      USE DIAG_COM, only : modd5k
      USE GEOM, only : fcor,dxyp,dxv,ravpn,ravps
      USE DYNAMICS, only : pu,pv,sd,spa,dut,dvt,dsig,mrch
      !USE DIAG, only : diagcd
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, CHECKSUM, GET
      USE DOMAIN_DECOMP_1D, only : haveLatitude, NORTH, SOUTH
      USE FILEMANAGER, only : openunit
      IMPLICIT NONE


      REAL*8, DIMENSION(IM , grid%J_STRT_HALO:grid%J_STOP_HALO , LM) ::
     &        U,V,P,UT,VT
      REAL*8, DIMENSION(IM , grid%J_STRT_HALO:grid%J_STOP_HALO     ) ::
     &        PA,PB,FDU2,GY,FX1,GY1

      REAL*8 FLUXQ(IM),DUMMYS(IM),DUMMYN(IM)

C****  ESMF special case:
C****  FD, FX and FDU1 original global dimensions (IM,0:JM+1) or (IM,JM+1)
      REAL*8 FD(IM,grid%J_STRT-1:grid%J_STOP+1)
      REAL*8 FX(IM,grid%J_STRT-1:grid%J_STOP+1)
      REAL*8 FDU1(IM,grid%J_STRT_HALO:grid%J_STOP+1)

      REAL*8, ALLOCATABLE, SAVE :: FDU(:,:,:)
cgsfc      REAL*8, ALLOCATABLE:: FLUXU2(:,:,:), FLUXV2(:,:,:)

      INTEGER,SAVE :: IFIRST = 1
      INTEGER I,J,IP1,IM1,L,K  !@var I,J,IP1,IM1,L,K  loop variables
      REAL*8 VMASS,RVMASS,ALPH,PDT4,SDU,DT1,DT2,DT3,DT6,DT8,DT24
     *     ,FLUX,FLUXU,FLUXV
!**** ESMF: Work arrays FLUX_U,FLUX_V added for distributed parallelization
      REAL*8, DIMENSION(IM , grid%J_STRT_HALO:grid%J_STOP_HALO      ) :: 
     &                  FLUX_U, FLUX_V
      REAL*8 ASDU(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM-1)

      INTEGER :: I_0, I_1, J_0, J_1
      INTEGER :: J_0S, J_1S, J_0SG, J_1SG, J_0H, J_1H
      INTEGER :: J_START, J_STP   !@var J_START, J_STP exceptional loop indices (esmf)
      INTEGER, SAVE :: ALLOCATED = 1
      INTEGER, SAVE :: iunit_hang

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1, 
     &         J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &         J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,      
     &         J_STRT_STGR=J_0SG, J_STOP_STGR=J_1SG)


cgsfc      allocate (FLUXU(IM, J_0H:J_1H, LM)) 
cgsfc      allocate (FLUXV(IM, J_0H:J_1H, LM)) 

!****  ESMF special case:
!****  FDU original global dimension:  (IM,JM+1,LM)
      if( ALLOCATED == 1 ) then
         ALLOCATED = 0
         allocate ( FDU(IM, J_0H:J_1H, LM) )
      endif
        
C****
         IF(MODD5K.LT.MRCH) CALL DIAG5F (U,V)

      IF (haveLatitude(grid, J=JM)) then
        J_STP = JM + 1
      ELSE
        J_STP = J_1
      END IF

C****
C**** Find scaling factors FDU for const. pressure layers (time indep)
C****

      IF (IFIRST .eq. 1) then
        IFIRST = 0

!       call openunit("hang_debug",iunit_hang)

        DO J = J_0S, J_1S
          DO I=1,IM
            FD(I,J)=PSFMPT*DXYP(J)
          END DO
        END DO

        if (haveLatitude(grid, J=1)) then
          FD(:,0)    = 0.
          FD(:,1)    = 2.*PSFMPT*DXYP(1)
        endif
        if (haveLatitude(grid, J=JM)) then
          FD(:,JM)   = 2.*PSFMPT*DXYP(JM)
          FD(:,JM+1) = 0.
        endif

        call HALO_UPDATE( grid, FD(:,J_0H:J_1H), from=SOUTH )
        I=IM
        DO J = J_0, J_STP
        DO IP1=1,IM
          FDU1(I,J) = .25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
          I=IP1
        END DO
        END DO

        I=IM
        IM1=I-1
        call HALO_UPDATE( grid, FDU1(:,J_0H:J_1H), from=NORTH+SOUTH )
        DO J = J_0SG, J_1SG
        DO IP1=1,IM
         FDU2(I,J)=.25*(FDU1(IP1,J)+FDU1(I,J+1)+FDU1(I,J-1)+FDU1(IM1,J))
         IM1=I
         I=IP1
        END DO
        END DO

        DO J=J_0SG,J_1SG
        DO I=1,IM
          FDU2(I,J)=(5.*FDU1(I,J)-2.*FDU2(I,J))*by3
        END DO
        END DO

        DO L=LS1,LM
        DO J = J_0SG, J_1SG
        DO I=1,IM
          FDU(I,J,L)=FDU2(I,J)*DSIG(L)
        END DO
        END DO
        END DO
      END IF

      DT2=DT1/2.
      DT3=DT1*by3
      DT6=DT1*by6
      DT8=DT1/8.
      DT24=DT1*by12/2.
C****
C**** Scale UT and VT which may then be physically interpreted as
C**** the momentum components
C****
C**** USE FOURTH ORDER ON THE INTERIOR.
C**** CONSIDER FLUXES PU AND PV TO BE ZERO BEYOND POLES.
C**** CONSIDER P TO BE ZERO BEYOND THE POLES
      DO J = J_0S, J_1S
        DO I=1,IM
          FD(I,J)   = PA(I,J)*DXYP(J)
        END DO
      END DO
      if (haveLatitude(grid, J=1)) then
        FD(:,0)     = 0.
        FD(:,1)     = 2.*PA(1,1)*DXYP(1)
      ENDIF
      if (haveLatitude(grid, J=JM)) then
        FD(:,JM)    = 2.*PA(1,JM)*DXYP(JM)
        FD(:,JM+1)  = 0.
      ENDIF

      call HALO_UPDATE( grid, FD(:,J_0H:J_1H), from=SOUTH )
      I=IM
      DO J = J_0, J_STP
        DO IP1=1,IM
          FDU1(I,J) = .25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
          I=IP1
        END DO
      END DO

      call HALO_UPDATE( grid, FDU1(:,J_0H:J_1H), from=NORTH+SOUTH)

      I=IM
      IM1=I-1
      DO J = J_0SG, J_1SG
        DO IP1=1,IM
          FDU2(I,J)=.25*(FDU1(IP1,J)+FDU1(I,J+1)+FDU1(I,J-1)+
     &         FDU1(IM1,J))
          IM1=I
          I=IP1
        END DO
      END DO

      DO J = J_0SG, J_1SG
        DO I=1,IM
          FDU2(I,J)=(5.*FDU1(I,J)-2.*FDU2(I,J))*by3
        END DO
      END DO

      DO L=1,LS1-1
        DO J = J_0SG, J_1SG
          DO I=1,IM
            FDU(I,J,L)=FDU2(I,J)*DSIG(L)
          END DO
        END DO
      END DO

C     DUT=0.
C     DVT=0.
      DO L=1,LM
        If (haveLatitude(grid, J=1)) THEN
          DUT(:,1,L)=0.
          DVT(:,1,L)=0.
        END IF
        DO J = J_0SG, J_1SG
          DO I=1,IM
            UT(I,J,L)=UT(I,J,L)*FDU(I,J,L)
            VT(I,J,L)=VT(I,J,L)*FDU(I,J,L)
            DUT(I,J,L)=0.
            DVT(I,J,L)=0.
          END DO
        END DO
      END DO
C****
C**** BEGINNING OF LAYER LOOP : FIND HORIZONTAL FLUXES
C****

      DO 300 L=1,LM
!       write(iunit_hang,*)'L=',l
        I=IM
        DO 175 J = J_0, J_1
          DO 170 IP1=1,IM
            FX(I,J)=.5*(PU(IP1,J,L)+PU(I,J,L))
  170     I=IP1
  175   CONTINUE
C****
      call HALO_UPDATE( grid, PV(:,:,L), from=NORTH )
      DO 180 J = J_0S, J_1S
      DO 180 I=1,IM
      GY(I,J)=.5*(PV(I,J,L)+PV(I,J+1,L))
  180 CONTINUE
C**** ASSUME FLUXES PU AND PV ARE ZERO BEYOND POLE
      if (haveLatitude(grid, J=1)) then
      DO 185 I=1,IM
      GY(I,1)=.5*PV(I,2,L)
      FX(I,0)=0.
  185 CONTINUE
      ENDIF
      if (haveLatitude(grid, J=JM)) then
      DO 186 I=1,IM
      GY(I,JM)=.5*PV(I,JM,L)
      FX(I,JM+1)=0.
 186  CONTINUE
      ENDIF
C****
      call HALO_UPDATE( grid, FX(:,J_0H:J_1H), from=SOUTH+NORTH )
      I=IM
      IM1=I-1
      DO 192 J = J_0, J_1
      DO 190 IP1=1,IM
      FX1(I,J)=.25*(FX(IP1,J)+FX(IM1,J)+FX(I,J+1)+FX(I,J-1))
      IM1=I
  190 I=IP1
  192 CONTINUE
      call HALO_UPDATE( grid, GY, from=SOUTH+NORTH )
      I=IM
      IM1=I-1
      DO 195 J = J_0S, J_1S
      DO 193 IP1=1,IM
      GY1(I,J)=.25*(GY(IP1,J)+GY(IM1,J)+GY(I,J+1)+GY(I,J-1))
      IM1=I
  193 I=IP1
  195 CONTINUE
C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****
      call HALO_UPDATE( grid, FX(:,J_0H:J_1H), from=SOUTH )
      I=IM
C**** Contribution from the West-East mass flux
      DO 215 J = J_0SG, J_1SG
      DO 210 IP1=1,IM
      FLUX = DT3*(FX(I,J)+FX(I,J-1))
      FLUXU = FLUX*(U(IP1,J,L)+U(I,J,L))
      DUT(IP1,J,L) = DUT(IP1,J,L) + FLUXU
      DUT(I,J,L)   = DUT(I,J,L)   - FLUXU
      FLUXV = FLUX*(V(IP1,J,L)+V(I,J,L))
      DVT(IP1,J,L) = DVT(IP1,J,L) + FLUXV
      DVT(I,J,L)   = DVT(I,J,L)   - FLUXV
  210 I = IP1
  215 CONTINUE

      CALL HALO_UPDATE(grid, GY, from=SOUTH)
      CALL HALO_UPDATE(grid, GY1, from=+SOUTH)
      CALL HALO_UPDATE(grid, FX(:,J_0H:J_1H), from=SOUTH)
      CALL HALO_UPDATE(grid, FX1, from=SOUTH)
      CALL HALO_UPDATE(grid, U(:,:,L), from=NORTH+SOUTH)
      CALL HALO_UPDATE(grid, V(:,:,L), from=NORTH+SOUTH)

      If (.not. haveLatitude(grid, J=1)) THEN
        J=J_0-1
        I = IM
        DO IP1=1,IM
C****   Contribution from the South-North mass flux
          FLUX = DT3*(GY(I,J)+GY(IP1,J))
          FLUXU = FLUX*(U(I,J,L)+U(I,J+1,L))
          DUT(I,J+1,L) = DUT(I,J+1,L) + FLUXU
          FLUXV = FLUX*(V(I,J,L)+V(I,J+1,L))
          DVT(I,J+1,L) = DVT(I,J+1,L) + FLUXV
C****   Contribution from the Southwest-Northeast mass flux
          FLUX = DT6*(FX(I,J)+GY(IP1,J)-FX1(I,J)-GY1(IP1,J))
          FLUXU = FLUX*(U(IP1,J+1,L)+U(I,J,L))
          DUT(IP1,J+1,L) = DUT(IP1,J+1,L) + FLUXU
          FLUXV = FLUX*(V(IP1,J+1,L)+V(I,J,L))
          DVT(IP1,J+1,L) = DVT(IP1,J+1,L) + FLUXV
C****   Contribution from the Southeast-Northwest mass flux
          FLUX = DT6*(-FX(I,J)+GY(IP1,J)+FX1(I,J)-GY1(IP1,J))
          FLUXU = FLUX*(U(I,J+1,L)+U(IP1,J,L))
          DUT(I,J+1,L) = DUT(I,J+1,L) + FLUXU
          FLUXV = FLUX*(V(I,J+1,L)+V(IP1,J,L))
          DVT(I,J+1,L) = DVT(I,J+1,L) + FLUXV
          I=IP1
        END DO
      END IF

      I = IM
      DO 230 J=J_0S,J_1S
      DO 220 IP1=1,IM
C**** Contribution from the South-North mass flux
      FLUX = DT3*(GY(I,J)+GY(IP1,J))
      FLUXU = FLUX*(U(I,J,L)+U(I,J+1,L))
      DUT(I,J+1,L) = DUT(I,J+1,L) + FLUXU
      DUT(I,J,L)   = DUT(I,J,L)   - FLUXU
      FLUXV = FLUX*(V(I,J,L)+V(I,J+1,L))
      DVT(I,J+1,L) = DVT(I,J+1,L) + FLUXV
      DVT(I,J,L)   = DVT(I,J,L)   - FLUXV
C**** Contribution from the Southwest-Northeast mass flux
      FLUX = DT6*(FX(I,J)+GY(IP1,J)-FX1(I,J)-GY1(IP1,J))
      FLUXU = FLUX*(U(IP1,J+1,L)+U(I,J,L))
      DUT(IP1,J+1,L) = DUT(IP1,J+1,L) + FLUXU
      DUT(I,J,L)     = DUT(I,J,L)     - FLUXU
      FLUXV = FLUX*(V(IP1,J+1,L)+V(I,J,L))
      DVT(IP1,J+1,L) = DVT(IP1,J+1,L) + FLUXV
      DVT(I,J,L)     = DVT(I,J,L)     - FLUXV
C**** Contribution from the Southeast-Northwest mass flux
      FLUX = DT6*(-FX(I,J)+GY(IP1,J)+FX1(I,J)-GY1(IP1,J))
      FLUXU = FLUX*(U(I,J+1,L)+U(IP1,J,L))
      DUT(I,J+1,L) = DUT(I,J+1,L) + FLUXU
      DUT(IP1,J,L) = DUT(IP1,J,L) - FLUXU
      FLUXV = FLUX*(V(I,J+1,L)+V(IP1,J,L))
      DVT(I,J+1,L) = DVT(I,J+1,L) + FLUXV
      DVT(IP1,J,L) = DVT(IP1,J,L) - FLUXV
  220 I=IP1
  230 CONTINUE
C****
C****

      call HALO_UPDATE( grid, PU(:,:,L), from=SOUTH )
      I=IM
      IM1=I-1
C**** CONTRIBUTION FROM THE BIG STEP WEST-EAST MASS FLUX
      DO 245 J = J_0SG, J_1SG
      DO 240 IP1=1,IM
      FLUX=-DT24*(PU(I,J,L)+PU(I,J-1,L))
      FLUXU=FLUX*(U(IP1,J,L)+U(IM1,J,L))
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU
      DUT(IM1,J,L)=DUT(IM1,J,L)-FLUXU
      FLUXV=FLUX*(V(IP1,J,L)+V(IM1,J,L))
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV
      DVT(IM1,J,L)=DVT(IM1,J,L)-FLUXV
      IM1=I
  240 I=IP1
  245 CONTINUE
C****
C****
C**** CONTRIBUTION FROM THE BIG STEP SOUTH-NORTH MASS FLUX
C**** ESMF: temporary arrays FLUX_U, FLUX_V needed for correct 
C****       distributed parallel computation.

C**** ESMF: First, make sure all needed halo data is in place
      CALL HALO_UPDATE(grid, U(:,:,L), from=NORTH+SOUTH)
      CALL HALO_UPDATE(grid, V(:,:,L), from=NORTH+SOUTH)

      I=IM
      DO J=J_0S,J_1S
        DO IP1=1,IM
          FLUX=-DT24*(PV(I,J,L)+PV(IP1,J,L))
          FLUX_U(I,J)=FLUX*(U(I,J-1,L)+U(I,J+1,L))
          FLUX_V(I,J)=FLUX*(V(I,J-1,L)+V(I,J+1,L))
          I=IP1
        END DO
      END DO
      CALL HALO_UPDATE(grid, FLUX_U, from=NORTH+SOUTH)
      CALL HALO_UPDATE(grid, FLUX_V, from=NORTH+SOUTH)

C**** ESMF: update latitutudes ahead (j+1) first.
      J_START = max(3, J_0H)
 
      DO J=J_START,J_1S
        DO I=1,IM
          DUT(I,J+1,L)=DUT(I,J+1,L)+FLUX_U(I,J)
          DVT(I,J+1,L)=DVT(I,J+1,L)+FLUX_V(I,J)
        END DO
      END DO

C**** ESMF: update latitudes "behind" (j-1) next.
      J_STP = min(J_1+1,JM-1)

      DO J=J_0S+1,J_STP
        DO I=1,IM
          DUT(I,J-1,L)=DUT(I,J-1,L)-FLUX_U(I,J)
          DVT(I,J-1,L)=DVT(I,J-1,L)-FLUX_V(I,J)
        END DO
      END DO


C****
C**** GISS/ESMF exception above.
C****
  300 CONTINUE
C****
C**** VERTICAL ADVECTION OF MOMENTUM
C****
C     DO 310 L=1,LM-1
C     DO 310 J=2,JM
C     I=IM
C     DO 310 IP1=1,IM
C     SDU=DT2*((SD(I,J-1,L)+SD(IP1,J-1,L))*RAVPN(J-1)+
C    *  (SD(I,J,L)+SD(IP1,J,L))*RAVPS(J))
C     DUT(I,J,L)  =DUT(I,J,L)  +SDU*(U(I,J,L)+U(I,J,L+1))
C     DUT(I,J,L+1)=DUT(I,J,L+1)-SDU*(U(I,J,L)+U(I,J,L+1))
C     DVT(I,J,L)  =DVT(I,J,L)  +SDU*(V(I,J,L)+V(I,J,L+1))
C     DVT(I,J,L+1)=DVT(I,J,L+1)-SDU*(V(I,J,L)+V(I,J,L+1))
C 310 I=IP1
      call HALO_UPDATE( grid, SD, from=SOUTH )
      DO L=1,LM-1
cgsfc      DO J=2,JM
      DO J = J_0SG, J_1SG
         DO I=1,IM-1
            ASDU(I,J,L)=DT2*((SD(I,J-1,L)+
     *                        SD(I+1,J-1,L))
     *      *RAVPN(J-1) +
     *      (SD(I,J,L)+SD(I+1,J,L))*RAVPS(J))
         END DO
         ASDU(IM,J,L)=DT2*((SD(IM,J-1,L)+SD(1,J-1,L))
     *                *RAVPN(J-1)+
     *                (SD(IM,J,L)+SD(1,J,L))*RAVPS(J))
      END DO
      END DO
      L=1
      DO J = J_0SG, J_1SG
        DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
        DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
      END DO
      DO L=2,LM-1
      DO J = J_0SG, J_1SG
         DUT(:,J,L)  =DUT(:,J,L)  -ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  -ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
      END DO
      END DO
      L=LM
      DO J = J_0SG, J_1SG
         DUT(:,J,L)=DUT(:,J,L)-ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DVT(:,J,L)=DVT(:,J,L)-ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
      END DO
C**** CALL DIAGNOSTICS
         IF(MODD5K.LT.MRCH) CALL DIAG5D (4,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAGCD (GRID,1,U,V,DUT,DVT,DT1)
      DO L=1,LM
      DO J = J_0SG, J_1SG
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
        VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO

C****
C**** CORIOLIS FORCE
C****

      call HALO_UPDATE( grid,  P, from=SOUTH )
      DO 430 L=1,LM
      IM1=IM
      DO 405 I=1,IM
C     FD(I,1)=FCOR(1)*2.-.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
C     FD(I,JM)=FCOR(JM)*2.+.5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
C**** Set the Coriolis term to zero at the Poles:
        If (haveLatitude(grid, J=1))
     &       FD(I,1) = -.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
        If (haveLatitude(grid, J=JM))
     &       FD(I,JM)=  .5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
  405 IM1=I
      DO 410 J = J_0S, J_1S
      DO 410 I=1,IM
      FD(I,J)=FCOR(J)+.25*(SPA(IM1,J,L)+SPA(I,J,L))*(DXV(J)-DXV(J+1))
  410 IM1=I
      call HALO_UPDATE( grid, FD(:,J_0H:J_1H), from=SOUTH )
  415 CONTINUE
      DO 425 J = J_0SG, J_1SG
      DO 420 I=1,IM
      PDT4=DT8*(P(I,J-1,L)+P(I,J,L))
      ALPH=PDT4*(FD(I,J)+FD(I,J-1))*DSIG(L)
      DUT(I,J,L)=DUT(I,J,L)+ALPH*V(I,J,L)
      DUT(IM1,J,L)=DUT(IM1,J,L)+ALPH*V(IM1,J,L)
      DVT(I,J,L)=DVT(I,J,L)-ALPH*U(I,J,L)
      DVT(IM1,J,L)=DVT(IM1,J,L)-ALPH*U(IM1,J,L)
  420 IM1=I
  425 CONTINUE
  430 CONTINUE
C**** CALL DIAGNOSTICS, ADD CORIOLIS FORCE INCREMENTS TO UT AND VT
         IF(MODD5K.LT.MRCH) CALL DIAG5D (5,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAGCD (GRID,2,U,V,DUT,DVT,DT1)
      DO L=1,LM
      DO J=J_0SG,J_1SG
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)+DUT(I,J,L)
        VT(I,J,L)=VT(I,J,L)+DVT(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO
C****
C**** UNDO SCALING PERFORMED AT BEGINNING OF ADVECV
C****
      DO J = J_0S, J_1S
      DO I=1,IM
        FD(I,J)   = PB(I,J)*DXYP(J)
      END DO
      END DO

      J_START = J_0
      J_STP   = J_1

      IF (haveLatitude(grid, J=1)) then
        FD(:,0)     = 0.
        FD(:,1)     = 2.*PB(1,1)*DXYP(1)
      ENDIF
      IF (haveLatitude(grid, J=JM)) then
        FD(:,JM)    = 2.*PB(1,JM)*DXYP(JM)
        FD(:,JM+1)  = 0.
        J_STP = JM + 1
      ENDIF

      call HALO_UPDATE( grid, FD(:,J_0H:J_1H), from=SOUTH )
      I=IM
      DO J = J_START, J_STP
      DO IP1=1,IM
        FDU1(I,J) = .25*(FD(I,J)+FD(IP1,J)+FD(I,J-1)+FD(IP1,J-1))
        I=IP1
      END DO
      END DO

      call HALO_UPDATE( grid, FDU1(:,J_0H:J_1H), from=NORTH+SOUTH )
      I=IM
      IM1=I-1
      DO J = J_0SG, J_1SG
      DO IP1=1,IM
        FDU2(I,J)=.25*(FDU1(IP1,J)+FDU1(I,J+1)+FDU1(I,J-1)+FDU1(IM1,J))
        IM1=I
        I=IP1
      END DO
      END DO

      DO J = J_0SG, J_1SG
      DO I=1,IM
        FDU2(I,J)=(5.*FDU1(I,J)-2.*FDU2(I,J))*by3
      END DO
      END DO


      DO L=1,LS1-1
      DO J = J_0SG, J_1SG
      DO I=1,IM
        FDU(I,J,L)=FDU2(I,J)*DSIG(L)
      END DO
      END DO
      END DO


C     DUT=0.
C     DVT=0.
      DO L=1,LM
        If (haveLatitude(grid, J=1)) THEN
          DUT(:,1,L)=0.
          DVT(:,1,L)=0.
        End If
      DO J = J_0SG, J_1SG
      DO I=1,IM
        UT(I,J,L)=UT(I,J,L)/FDU(I,J,L)
        VT(I,J,L)=VT(I,J,L)/FDU(I,J,L)
        DUT(I,J,L)=0.
        DVT(I,J,L)=0.
      END DO
      END DO
      END DO

      RETURN
      END SUBROUTINE ADVECV

      end module MOMENTS
