C*****************************************************************************
C PLEASE KEEP THIS NOTE OF MODEL-DEVELOPMENT HISTORY
C Matrix solve uses Thomas algorithm, 10/1991, Jinlun Zhang  
C Spherical coordinate system, 10/27/93, Jinlun Zhang
C Latest finite differencing scheme for treatment of NP,9/9/1996,Jinlun Zhang
C Alternating direction implicit (ADI) method is used, 10/1998, Jinlun Zhang
C For details about ADI dynamics model, see Zhang and Rothrock, "Modeling
C   Arctic Sea ice with an efficient plastic solution", submitted to JGR, 1999
C*****************************************************************************
      SUBROUTINE RELAX(UICE,VICE,ETA,ZETA,DRAGSU,DRAGSV,DRAGAU,DRAGAV
     &,AMASSU,AMASSV,FORCEX,FORCEY,ERROR,UICEC,VICEC)
      implicit none

c********************************
      USE DOMAIN_DECOMP_ICE, ONLY : grid_ICDYN,AM_I_ROOT,HALO_UPDATE
      USE DOMAIN_DECOMP_1D, ONLY  : hasSouthPole, hasNorthPole
      USE TRIDIAG_MOD, only : TRIDIAG_new, TRIDIAG_cyclic
      IMPLICIT NONE

      REAL*8, DIMENSION(grid_ICDYN%I_STRT_HALO:grid_ICDYN%I_STOP_HALO
     &     ,grid_ICDYN%J_STRT_HALO:
     &     grid_ICDYN%J_STOP_HALO) ::
     &         AU,BU,CU,FXY,FXY1
      REAL*8, DIMENSION(grid_ICDYN%I_STRT_HALO:grid_ICDYN%I_STOP_HALO,
     &     ,grid_ICDYN%J_STRT_HALO:
     &     grid_ICDYN%J_STOP_HALO) ::
     &         AV,BV,CV,VRT,U_tmp
      REAL*8, DIMENSION(grid_ICDYN%I_STRT_HALO:grid_ICDYN%I_STOP_HALO,
     &     ,grid_ICDYN%J_STRT_HALO:
     &     grid_ICDYN%J_STOP_HALO) ::
     &         FXYa,FXY1a,ZETAPETA,ZETAMETA
      real*8, dimension(grid_ICDYN%ISTRT_HALO:grid_ICDYN%ISTOP_HALO,
     &     grid_ICDYN%JSTRT_HALO:grid_ICDYN%JSTOP_HALO) :: 
     &     AA1,AA2,AA3,AA4,AA5,AA6,AA7,AA8,AA9
      real*8 AA10
      
      INTEGER :: J_0,J_1, J_0S,J_1S, J_0H,J_1H
      INTEGER :: I_0,I_1, I_0S,I_1S, I_0H,I_1H

C****
C**** Extract local domain parameters from "grid"
C****
      call getDomainBounds(grid_ICDYN, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H,   J_STOP_HALO=J_1H,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S )
      call getDomainBounds(grid_ICDYN, I_STRT     =I_0,    I_STOP     =I_1,
     &               I_STRT_HALO=I_0H,   I_STOP_HALO=I_1H,
     &               I_STRT_SKP =I_0S,   I_STOP_SKP =I_1S )


      if (hasNorthPole(grid_ICDYN)) then
         NYPOLEU=NY
         NYPOLEV=NY1
         NPOL=1
      else
c*** this needs to be fixed, not working on CS grid nor tripolar grid
         NYPOLEU=NYM1  ! fix this
         NYPOLEV=NY    ! fix this
         NPOL=0
      endif


      DO J=J_0,J_1S
      DO I=I_0,I_1
      FORCEX(I,J)=FORCEX(I,J)*UUM(I,J)
      FORCEY(I,J)=FORCEY(I,J)*VVM(I,J)
      END DO
      END DO

C MUST UPDATE HEFF BEFORE CALLING RELAX
C FIRST SET U(2)=U(1)
      DO J=J_0,J_1S
         DO I=I_0,I_1
C     NOW MAKE SURE BDRY PTS ARE EQUAL TO ZERO
            UICE(I,J,2)=UICE(I,J,1)
            VICE(I,J,2)=VICE(I,J,1)
            UICE(I,J,1)=UICE(I,J,3)*UUM(I,J)
            VICE(I,J,1)=VICE(I,J,3)*VVM(I,J)
         enddo
      enddo

      if (hasNorthPole(grid_ICDYN)) then
         NXLCYC=NX
         LCYC=1
         DO I=1,NX1/2
            UICE(I,NY+1,1)=-UICEC(I+(NX1-2)/2,NY)
            UICE(I,NY+1,3)=-UICEC(I+(NX1-2)/2,NY)
            
            UICEC(I,NY+1)=-UICEC(I+(NX1-2)/2,NY)
            VICEC_NP(I)   =-VICEC(I+(NX1-2)/2,NY1)
         END DO
         
         DO I=NX1/2+1,NX
            UICE(I,NY+1,1)=-UICEC(I-(NX1-2)/2,NY)
            UICE(I,NY+1,3)=-UICEC(I-(NX1-2)/2,NY)
            
            UICEC(I,NY+1)=-UICEC(I-(NX1-2)/2,NY)
            VICEC_NP(I)  =-VICEC(I-(NX1-2)/2,NY1)
         END DO
         
         
         DO J=1,NY1
            UICE(1,J,1)=UICEC(NX,J)
            VICE(1,J,1)=VICEC(NX,J)
            UICE(NX1,J,1)=UICEC(2,J)
            VICE(NX1,J,1)=VICEC(2,J)
            UICE(1,J,3)=UICEC(NX,J)
            VICE(1,J,3)=VICEC(NX,J)
            UICE(NX1,J,3)=UICEC(2,J)
            VICE(NX1,J,3)=VICEC(2,J)
            
            UICEC(1,J)=UICEC(NX,J)
            VICEC(1,J)=VICEC(NX,J)
            UICEC(NX1,J)=UICEC(2,J)
            VICEC(NX1,J)=VICEC(2,J)
         END DO
         
         VICEC_NP(1)=VICEC_NP(NX)
         VICEC_NP(NX1)=VICEC_NP(2)
      else
         NXLCYC=NXM1
         LCYC=0
      endif

C GET UICEC at V location and VICEC at U location
      DO J=J_0,J_1S
      DO I=I_0,I_1
      UICECV(I,J)=0.25*(UICEC(I,J)+UICEC(I+1,J)+UICEC(I+1,J-1)
     &+UICEC(I,J-1))
      IF(J.LT.NY1)
     &VICECU(I,J)=0.25*(VICEC(I,J)+VICEC(I,J+1)+VICEC(I-1,J+1)
     &+VICEC(I-1,J))
      END DO
      END DO

      if (hasNorthPole(grid_ICDYN)) then
         DO I=1,NX1/2
            VICECU(I,NY+1)=0.5*(VICECU(I,NY)-VICECU(I+(NX1-2)/2,NY))
         END DO
         
         DO I=NX1/2+1,NX
            VICECU(I,NY+1)=0.5*(VICECU(I,NY)-VICECU(I-(NX1-2)/2,NY))
         END DO
         
         DO J=1,NY1
            UICECV(1,J)=UICECV(NX,J)
            VICECU(1,J)=VICECU(NX,J)
            UICECV(NX1,J)=UICECV(2,J)
            VICECU(NX1,J)=VICECU(2,J)
         END DO
      endif

1413  FORMAT(2I4,6E12.4)

C FIRST DO UICE
C THE FIRST HALF

C**Update halos 
      CALL HALO_UPDATE(grid_ICDYN, ETA)
      CALL HALO_UPDATE(grid_ICDYN, ZETA)
      CALL HALO_UPDATE(grid_ICDYN, ZETAPETA)
      CALL HALO_UPDATE(grid_ICDYN, ZETAMETA)
      CALL HALO_UPDATE(grid_ICDYN, VICE)
      CALL HALO_UPDATE(grid_ICDYN, KXT)


      DO J=J_0S,J_1S
      DO I=I_0,I_1

      FXY(I,J)=DRAGAU(I,J)*VICECU(I,J)+FORCEX(I,J)
     1+DXCR(I,J)*(ZETAMETA(I,J)*(VICECU(I,J+1)-VICECU(I,J))*DYFR(I,J)
     1-ZETAMETA(I-1,J)*(VICECU(I-1,J+1)-VICECU(I-1,J))*DYFR(I-1,J))
     1
     2+0.25*KYU(I,J)*(ZETA(I-1,J)+ZETA(I,J))*DXCR(I,J)
     2*((VICECU(I,J+1)+VICECU(I,J))-(VICECU(I-1,J+1)+VICECU(I-1,J)))
     2
     3+0.25*(VICECU(I-1,J)+VICECU(I,J)+VICECU(I,J+1)+VICECU(I-1,J+1))
     3*DXCR(I,J)*(KYT(I,J)*ZETAPETA(I,J)-KYT(I-1,J)*ZETAPETA(I-1,J))
     3
     4-3.0*0.25*KXU(I,J)*(ETA(I-1,J)+ETA(I,J))*DYGR(I,J)
     4*((VICECU(I-1,J+1)+VICECU(I,J+1))-(VICECU(I-1,J)+VICECU(I,J)))
     4
     5+3.0*0.25*KYU(I,J)*(ETA(I-1,J)+ETA(I,J))*DXCR(I,J)
     5*((VICECU(I,J)+VICECU(I,J+1))-(VICECU(I-1,J)+VICECU(I-1,J+1)))
     5
     6+0.25*DYGR(I,J)*((ETA(I-1,J+1)+ETA(I,J+1)+ETA(I,J)+ETA(I-1,J))
     6*DXVR(I,J+1)*(VICECU(I,J+1)-VICECU(I-1,J+1))
     6                -(ETA(I-1,J-1)+ETA(I,J-1)+ETA(I,J)+ETA(I-1,J))
     6*DXVR(I,J  )*(VICECU(I,J  )-VICECU(I-1,J  )))
     6
     7-0.25*0.25*DYGR(I,J)
     7*(VICECU(I-1,J)+VICECU(I,J)+VICECU(I,J+1)+VICECU(I-1,J+1))
     7*((KXT(I-1,J)*ETA(I-1,J)+KXT(I,J)*ETA(I,J)+KXT(I,J+1)*ETA(I,J+1)
     7+KXT(I-1,J+1)*ETA(I-1,J+1))-(KXT(I-1,J-1)*ETA(I-1,J-1)+KXT(I,J-1)
     7*ETA(I,J-1)+KXT(I,J)*ETA(I,J)+KXT(I-1,J)*ETA(I-1,J)))

      END DO
      END DO

c      print *, 'FXY for UICE  '
c      print 15, (FXY(20,J),J=1,JMT)

      DO J=J_0S,J_1S
      DO I=I_0,I_1
      AA1(I,J)=DXCR(I,J)*DXFR(I,J)*ZETAPETA(I,J)
      AA2(I,J)=DXCR(I,J)*DXFR(I-1,J)*ZETAPETA(I-1,J)
      AA3(I,J)=0.25*DYGR(I,J)*DYUR(I,J+1)
     &*(ETA(I-1,J)+ETA(I,J)+ETA(I,J+1)+ETA(I-1,J+1))
      AA4(I,J)=0.25*DYGR(I,J)*DYUR(I,J)
     &*(ETA(I-1,J-1)+ETA(I,J-1)+ETA(I,J)+ETA(I-1,J))
      AA5(I,J)=-DXCR(I,J)*(KXT(I,J)*ZETAMETA(I,J)-KXT(I-1,J)
     &*ZETAMETA(I-1,J))
      AA6(I,J)=0.5*KXU(I,J)*(ZETAPETA(I-1,J)+ZETAPETA(I,J))
     &/(DXF(I-1,J)+DXF(I,J))
      AA7(I,J)=0.5*KYU(I,J)*(ETA(I-1,J)+ETA(I,J))/(DYU(I,J+1)+DYU(I,J))
      AA8(I,J)=0.25*DYGR(I,J)*(
     &(KYT(I-1,J)*ETA(I-1,J)+KYT(I,J)*ETA(I,J)+KYT(I,J+1)*ETA(I,J+1)
     &+KYT(I-1,J+1)*ETA(I-1,J+1))-(KYT(I-1,J-1)*ETA(I-1,J-1)+KYT(I,J-1)
     &*ETA(I,J-1)+KYT(I,J)*ETA(I,J)+KYT(I-1,J)*ETA(I-1,J)))
      AA9(I,J)=(KXU(I,J)*KXU(I,J)+KYU(I,J)*KYU(I,J))
     &*(ETA(I-1,J)+ETA(I,J))
      AU(I,J)=(-AA2(I,J)+AA6(I,J))*UUM(I,J)
      BU(I,J)=(AA1(I,J)+AA2(I,J)+AA5(I,J)+AA8(I,J)+AA9(I,J)
     &+AMASSU(I,J)*BYDTS*2.0+DRAGSU(I,J))*UUM(I,J)+(1.0-UUM(I,J))
      CU(I,J)=(-AA1(I,J)-AA6(I,J))*UUM(I,J)
      END DO
      END DO
      DO J=J_0S,J_1
      AU(2,J)=0.0
      CU(NXLCYC,J)=0.0
      CU(2,J)=CU(2,J)/BU(2,J)
      END DO

      DO 1200 J=J_0S,J_1
      DO I=I_0,I_1

      IF(I.EQ.2) THEN
      AA10=(AA2(I,J)-AA6(I,J))*UICEC(I-1,J)*UUM(I,J)*FLOAT(LCYC-0)
      ELSE IF(I.EQ.NXLCYC) THEN
      AA10=(AA1(I,J)+AA6(I,J))*UICEC(I+1,J)*UUM(I,J)*FLOAT(LCYC-0)
      ELSE
      AA10=0.0
      END IF

      URT(I)=AA10+FXY(I,J)-(AA3(I,J)+AA4(I,J))*UICE(I,J,2)
     &+(AA3(I,J)+AA7(I,J))*UICE(I,J+1,2)
     &+(AA4(I,J)-AA7(I,J))*UICE(I,J-1,2)
      URT(I)=(URT(I)+AMASSU(I,J)*BYDTS*UICE(I,J,2)*2.0)*UUM(I,J)
      END DO

ccc---> begin tridiag algo
c      DO I=2,NXLCYC
c      CUU(I)=CU(I,J)
c      END DO
c      URT(2)=URT(2)/BU(2,J)
c      DO I=3,NXLCYC
c      IM=I-1
c      CUU(I)=CUU(I)/(BU(I,J)-AU(I,J)*CUU(IM))
c      URT(I)=(URT(I)-AU(I,J)*URT(IM))/(BU(I,J)-AU(I,J)*CUU(IM))
c      END DO
c      DO I=1,NXLCYC-2
c      J1=NXLCYC-I
c      J2=J1+1
c      URT(J1)=URT(J1)-CUU(J1)*URT(J2)
c      END DO
c      DO I=2,NXLCYC
c      UICE(I,J,1)=URT(I)
c      END DO
#ifdef CUBED_SPHERE_ICE
c***  We need a non cyclic tridiag here that takes care of bndry 
c***  conditions between processors
#else
      CALL TRIDIAG_cyclic(AU(2,J),BU(2,J),CU(2,J),URT(2),UICE(2,J,1),
     &               NXLCYC-1)
#endif
ccc<--- end tridiag
1200  CONTINUE

c      print *, 'UICE 1st  '
c      print 15, (UICE(20,J,1),J=1,JMT)
      DO J=J_0S,J_1
      DO I=I_0,I_1
      UICE(I,J,3)=UICE(I,J,1)
      END DO
      END DO

C NOW THE SECOND HALF
      DO J=J_0S,J_1
      DO I=I_0,I_1
      AV(I,J)=-(AA4(I,J)-AA7(I,J))*UUM(I,J)
      BV(I,J)=(AA3(I,J)+AA4(I,J)+AA5(I,J)+AA8(I,J)+AA9(I,J)
     &+AMASSU(I,J)*BYDTS*2.0+DRAGSU(I,J))*UUM(I,J)+(1.0-UUM(I,J))
      CV(I,J)=-(AA3(I,J)+AA7(I,J))*UUM(I,J)
      END DO
      END DO

      if (hasSouthPole(grid_ICDYN)) then
      DO I=2,NXLCYC
      AV(I,2)=0.0
      CV(I,NYPOLEU)=0.0  ! fix this
c      CV(I,2)=CV(I,2)/BV(I,2)
      END DO
      endif

      DO I=I_0,I_1
      DO J=J_0S,J_1
      IF(J.EQ.NYPOLEU) THEN  ! fix this
      AA10=(AA3(I,J)+AA7(I,J))*UICEC(I,J+1)*UUM(I,J)*FLOAT(NPOL-0) 
      ELSE
      AA10=0.0
      END IF

      FXY1(I,J)=AA10+AMASSU(I,J)*BYDTS*UICE(I,J,1)*2.0
     &-(AA1(I,J)+AA2(I,J))*UICE(I,J,1)
     6+(AA1(I,J)+AA6(I,J))*UICE(I+1,J,1)
     6+(AA2(I,J)-AA6(I,J))*UICE(I-1,J,1)
      
      END DO
      END DO

      DO I=I_0,I_1
      DO J=J_0S,J_1
      VRT(J)=FXY(I,J)+FXY1(I,J)
      VRT(J)=VRT(J)*UUM(I,J)
      END DO
      END DO
ccc---> begin tridiag algo
c      DO J=2,NYPOLEU
c      CVV(J)=CV(I,J)
c      END DO
c      VRT(2)=VRT(2)/BV(I,2)
c      DO J=3,NYPOLEU
c      JM=J-1
c      CVV(J)=CVV(J)/(BV(I,J)-AV(I,J)*CVV(JM))
c      VRT(J)=(VRT(J)-AV(I,J)*VRT(JM))/(BV(I,J)-AV(I,J)*CVV(JM))
c      END DO
c      DO J=1,NYPOLEU-2
c      J1=NYPOLEU-J
c      J2=J1+1
c      VRT(J1)=VRT(J1)-CVV(J1)*VRT(J2)
c      END DO
c      DO J=2,NYPOLEU
c      UICE(I,J,1)=VRT(J)
c      END DO

#ifdef CUBED_SPHERE_ICE
c***   We need a 2d decomposed Thomas tridiag algo here
#else
      CALL TRIDIAG_new(AV, BV, CV, VRT, U_tmp, grid_ICDYN, 2, NYPOLE)
#endif

      DO I=I_0,I_1
      DO J=J_0S,J_1
        UICE(I,J,1)=U_tmp(I,J)     
      END DO
      END DO
c      print *, 'UICE 2st  '
c      print 15, (UICE(20,J,1),J=1,JMT)


C NOW DO VICE
C THE FIRST HALF
#ifdef CUBED_SPHERE_ICE

#else
      CALL HALO_UPDATE(GRID_ICDYN,  ETA,FROM=NORTH)
      CALL HALO_UPDATE(GRID_ICDYN, ZETA,FROM=NORTH)
#endif

      DO I=I_0,I_1
      DO J=J_0S,J_1
      FXY(I,J)=-DRAGAV(I,J)*UICECV(I,J)+FORCEY(I,J)
     1+DYCR(I,J)*(ZETAMETA(I,J)*(UICECV(I+1,J)-UICECV(I,J))*DXFR(I,J)
     1-ZETAMETA(I,J-1)*(UICECV(I+1,J-1)-UICECV(I,J-1))*DXFR(I,J-1))
     1
     2+0.25*KXV(I,J)*(ZETA(I,J-1)+ZETA(I,J))*DYCR(I,J)
     2*((UICECV(I+1,J)+UICECV(I,J))-(UICECV(I+1,J-1)+UICECV(I,J-1)))
     2
     3+0.25*(UICECV(I,J-1)+UICECV(I,J)+UICECV(I+1,J)+UICECV(I+1,J-1))
     3*DYCR(I,J)*(KXT(I,J)*ZETAPETA(I,J)-KXT(I,J-1)*ZETAPETA(I,J-1))
     3
     4-3.0*0.25*KYV(I,J)*(ETA(I,J-1)+ETA(I,J))*DXGR(I,J)
     4*((UICECV(I+1,J-1)+UICECV(I+1,J))-(UICECV(I,J-1)+UICECV(I,J)))
     4
     5+3.0*0.25*KXV(I,J)*(ETA(I,J-1)+ETA(I,J))*DYCR(I,J)
     5*((UICECV(I,J)+UICECV(I+1,J))-(UICECV(I,J-1)+UICECV(I+1,J-1)))
     5
     6+0.25*DXGR(I,J)*((ETA(I+1,J-1)+ETA(I+1,J)+ETA(I,J)+ETA(I,J-1))
     6*DYUR(I+1,J)*(UICECV(I+1,J)-UICECV(I+1,J-1))
     6                -(ETA(I-1,J-1)+ETA(I,J-1)+ETA(I,J)+ETA(I-1,J))
     6*DYUR(I,J  )*(UICECV(I,J  )-UICECV(I,J-1  )))
     6
     7-0.25*0.25*DXGR(I,J)
     7*(UICECV(I,J-1)+UICECV(I,J)+UICECV(I+1,J)+UICECV(I+1,J-1))
     7*((KYT(I,J-1)*ETA(I,J-1)+KYT(I,J)*ETA(I,J)+KYT(I+1,J)*ETA(I+1,J)
     7+KYT(I+1,J-1)*ETA(I+1,J-1))-(KYT(I-1,J-1)*ETA(I-1,J-1)+KYT(I,J-1)
     7*ETA(I,J-1)+KYT(I,J)*ETA(I,J)+KYT(I-1,J)*ETA(I-1,J)))

      END DO
      END DO

c      print *, 'FXY for VICE  '
c      print 15, (FXY(20,J),J=1,JMT)

      DO I=I_0,I_1
      DO J=J_0S,J_1

      AA1(I,J)=DYCR(I,J)*DYFR(I,J)*ZETAPETA(I,J)
      AA2(I,J)=DYCR(I,J)*DYFR(I,J-1)*ZETAPETA(I,J-1)
      AA3(I,J)=0.25*DXGR(I,J)*DXVR(I+1,J)
     &*(ETA(I,J-1)+ETA(I,J)+ETA(I+1,J)+ETA(I+1,J-1))
      AA4(I,J)=0.25*DXGR(I,J)*DXVR(I,J)
     &*(ETA(I-1,J-1)+ETA(I,J-1)+ETA(I,J)+ETA(I-1,J))
      AA5(I,J)=-DYCR(I,J)*(KYT(I,J)*ZETAMETA(I,J)-KYT(I,J-1)
     &*ZETAMETA(I,J-1))
      AA6(I,J)=0.5*KYV(I,J)*(ZETAPETA(I,J-1)+ZETAPETA(I,J))
     &/(DYF(I,J-1)+DXF(I,J))
      AA7(I,J)=0.5*KXV(I,J)*(ETA(I,J-1)+ETA(I,J))/(DXV(I+1,J)+DXV(I,J))
      AA8(I,J)=0.25*DXGR(I,J)*(
     &(KXT(I,J-1)*ETA(I,J-1)+KXT(I,J)*ETA(I,J)+KXT(I+1,J)*ETA(I+1,J)
     &+KXT(I+1,J-1)*ETA(I+1,J-1))-(KXT(I-1,J-1)*ETA(I-1,J-1)+KXT(I,J-1)
     &*ETA(I,J-1)+KXT(I,J)*ETA(I,J)+KXT(I-1,J)*ETA(I-1,J)))
      AA9(I,J)=(KXV(I,J)*KXV(I,J)+KYV(I,J)*KYV(I,J))
     &*(ETA(I,J-1)+ETA(I,J))
      AV(I,J)=(-AA2(I,J)+AA6(I,J))*VVM(I,J)
      BV(I,J)=(AA1(I,J)+AA2(I,J)+AA5(I,J)+AA8(I,J)+AA9(I,J)
     &+AMASSV(I,J)*BYDTS*2.0+DRAGSV(I,J))*VVM(I,J)+(1.0-VVM(I,J))
      CV(I,J)=(-AA1(I,J)-AA6(I,J))*VVM(I,J)
      END DO
      END DO
      DO I=I_0,I_1
      AV(I,2)=0.0
      CV(I,NYPOLEV)=0.0  ! fix this
      CV(I,2)=CV(I,2)/BV(I,2)
      END DO

      DO I=I_0,I_1
      DO J=I_0,I_1
      IF(J.EQ.NYPOLEV) THEN   ! fix this
c      AA10=(AA1(I,J)+AA6(I,J))*VICEC(I,J+1)*VVM(I,J)*FLOAT(NPOL-0)
      AA10=(AA1(I,J)+AA6(I,J))*VICEC_NP(I)*VVM(I,J)*FLOAT(NPOL-0)
      ELSE
      AA10=0.0
      END IF
      VRT(J)=AA10+FXY(I,J)-(AA3(I,J)+AA4(I,J))*VICE(I,J,2)
     &+(AA3(I,J)+AA7(I,J))*VICE(I+1,J,2)
     &+(AA4(I,J)-AA7(I,J))*VICE(I-1,J,2)
      VRT(J)=(VRT(J)+AMASSV(I,J)*BYDTS*VICE(I,J,2)*2.0)*VVM(I,J)
      END DO
      END DO
c*** ---> begin tridiag 
c
c      DO J=2,NYPOLEV
c      CVV(J)=CV(I,J)
c      END DO
c      VRT(2)=VRT(2)/BV(I,2)
c      DO J=3,NYPOLEV
c      JM=J-1
c      CVV(J)=CVV(J)/(BV(I,J)-AV(I,J)*CVV(JM))
c      VRT(J)=(VRT(J)-AV(I,J)*VRT(JM))/(BV(I,J)-AV(I,J)*CVV(JM))
c      END DO
c      DO J=1,NYPOLEV-2
c      J1=NYPOLEV-J
c      J2=J1+1
c      VRT(J1)=VRT(J1)-CVV(J1)*VRT(J2)
c      END DO
c      DO J=2,NYPOLEV
c      VICE(I,J,1)=VRT(J)
c      END DO
#ifdef CUBED_SPHERE_ICE
c***  We need a non cyclic tridiag here that takes care of bndry 
c***  conditions between processors
#else
       CALL TRIDIAG_new(AV, BV, CV, VRT, U_tmp, grid_ICDYN, 2, NYPOLE)
#endif

c      print *, 'VICE 1st  '
c      print 15, (VICE(20,J,1),J=1,JMT)
      DO I=I_0,I_1
      DO J=J_0,J_1
      VICE(I,J,3)=VICE(I,J,1)
      END DO
      END DO

C NOW THE SECOND HALF

      DO J=J_0,J_1
      DO I=I_0,I_1
      AU(I,J)=-(AA4(I,J)-AA7(I,J))*VVM(I,J)
      BU(I,J)=(AA3(I,J)+AA4(I,J)+AA5(I,J)+AA8(I,J)+AA9(I,J)
     &+AMASSV(I,J)*BYDTS*2.0+DRAGSV(I,J))*VVM(I,J)+(1.0-VVM(I,J))
      CU(I,J)=-(AA3(I,J)+AA7(I,J))*VVM(I,J)
      END DO
      END DO
      DO J=J_0,J_1
      AU(2,J)=0.0
      CU(NXLCYC,J)=0.0
      CU(2,J)=CU(2,J)/BU(2,J)
      END DO


      DO J=J_0S,J_1
      DO I=I_0,I_1

      IF(I.EQ.2) THEN
      AA10=(AA4(I,J)-AA7(I,J))*VICEC(I-1,J)*VVM(I,J)*FLOAT(LCYC-0) 
      ELSE IF(I.EQ.NXLCYC) THEN
      AA10=(AA3(I,J)+AA7(I,J))*VICEC(I+1,J)*VVM(I,J)*FLOAT(LCYC-0) 
      ELSE
      AA10=0.0
      END IF

      IF(J.EQ.NY1) THEN
      FXY1(I,J)=AA10+AMASSV(I,J)*BYDTS*VICE(I,J,1)*2.0
     &-(AA1(I,J)+AA2(I,J))*VICE(I,J,1)
     6+(AA1(I,J)+AA6(I,J))*VICEC_NP(I)
     6+(AA2(I,J)-AA6(I,J))*VICE(I,J-1,1)
      ELSE
      FXY1(I,J)=AA10+AMASSV(I,J)*BYDTS*VICE(I,J,1)*2.0
     &-(AA1(I,J)+AA2(I,J))*VICE(I,J,1)
     6+(AA1(I,J)+AA6(I,J))*VICE(I,J+1,1)
     6+(AA2(I,J)-AA6(I,J))*VICE(I,J-1,1)
      END IF
      
      END DO
      END DO

      DO J=J_0,J_1
      DO I=I_0,I_1
      URT(I)=FXY(I,J)+FXY1(I,J)
      URT(I)=URT(I)*VVM(I,J)
      END DO

c*** ---> begin tridiag 
c      DO I=2,NXLCYC
c      CUU(I)=CU(I,J)
c      END DO
c      URT(2)=URT(2)/BU(2,J)
c      DO I=3,NXLCYC
c      IM=I-1
c      CUU(I)=CUU(I)/(BU(I,J)-AU(I,J)*CUU(IM))
c      URT(I)=(URT(I)-AU(I,J)*URT(IM))/(BU(I,J)-AU(I,J)*CUU(IM))
c      END DO
c      DO I=1,NXLCYC-2
c      J1=NXLCYC-I
c      J2=J1+1
c      URT(J1)=URT(J1)-CUU(J1)*URT(J2)
c      END DO
c      DO I=2,NXLCYC
c      VICE(I,J,1)=URT(I)
c      END DO
#ifdef CUBED_SPHERE_ICE

#else
        CALL TRIDIAG_cyclic(AU(2,J),BU(2,J),CU(2,J),URT(2),VICE(2,J,1),
     &               NXLCYC-1)
#endif
      END DO

c      print *, 'VICE 2st  '
c      print 15, (VICE(20,J,1),J=1,JMT)

      if (hasNorthPole(grid_ICDYN)) then
      DO J=1,NY1
      UICE(1,J,1)=UICE(NX,J,1)
      VICE(1,J,1)=VICE(NX,J,1)
      UICE(NX1,J,1)=UICE(2,J,1)
      VICE(NX1,J,1)=VICE(2,J,1)
      END DO
      endif

      DO J=J_0S,J_1
      DO I=I_0,I_1
      UICE(I,J,1)=UICE(I,J,1)*UUM(I,J)
      VICE(I,J,1)=VICE(I,J,1)*VVM(I,J)
c       if(j.gt.NY1-3) 
c     &  write(*,'(2i5,2f4.0,6f8.4)') i,j, uum(i,j),vvm(i,j)
c     &, UICE(I,J,1),VICE(I,J,1), UICEC(I,J),VICEC(I,J)
c     &, UICECV(I,J),VICECU(I,J)
      END DO
      END DO
      call flush(6)

      RETURN
      END

#else
      subroutine relax
      return
      end
c
#endif
