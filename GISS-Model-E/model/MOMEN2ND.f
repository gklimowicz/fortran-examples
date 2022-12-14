      module MOMENTS
      USE RESOLUTION, only : JM
      implicit none
      private

      public ADVECV, moment_enq_order
      public initMoments
      !REAL*8, SAVE :: SMASS(JM)

      contains

      SUBROUTINE initMoments
!@sum Implemented to remove IFIRST in ADVEC
      USE GEOM, only : DXYV
      IMPLICIT NONE
      INTEGER :: J
      !DO J=2,JM
      !  SMASS(J)=PSFMPT*DXYV(J)
      !END DO
      END SUBROUTINE initMoments

      subroutine moment_enq_order(order)
!@sum moment_enq_order returns order of the scheme
      implicit none
      integer, intent(out) :: order
      order = 2
      return
      end subroutine moment_enq_order


      Subroutine ADVECV (DT1, U,V,MMEAN, MBEFOR,UT,VT,MAFTER)
!@sum  ADVECV Advects momentum (incl. coriolis) using mass fluxes
!@auth Original development team
!**** Input: DT1 = time step (s)
!****        U,V = mean horizontal velocity during time step (m/s)
!****     MBEFOR = mass distribution at start of time step (kg/m^2)
!****     MAFTER = mass distribution at end of time step (kg/m^2)
!****      MMEAN = mean mass distribution during time step (kg/m^2)
!****     PU=>MU = eastward mass flux (kg/s)
!****     PV=>MV = northward mass flux (kg/s)
!****     SD=>MW = downward vertical mass flux (kg/s)
!**** Output: UT,VT = horizontal velocity updated by advection (m/s)         
!****
      USE RESOLUTION, only : im,jm,lm
      USE DIAG_COM, only : modd5k
      USE DOMAIN_DECOMP_ATM, only : GRID
      Use DOMAIN_DECOMP_1D,  Only: HALO_UPDATE, NORTH,SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude, getDomainBounds
      USE GEOM, only : fcor,dxyv,dxyn,dxys,dxv,ravpn,ravps
     &     ,sini=>siniv,cosi=>cosiv,acor,polwt
      Use DYNAMICS,   Only: PU=>MU,PV=>MV,SD=>MW, SPA,DUT,DVT
      USE DYNAMICS, only : do_polefix,mrch
c      USE DIAG, only : diagcd
      IMPLICIT NONE
      REAL*8 U(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *       V(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
      REAL*8 UT(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *       VT(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
      Real*8 :: MMEAN(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *         MBEFOR(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     *         MAFTER(LM,IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)

!**** Local variables
      Real*8 :: FD(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      INTEGER I,J,IP1,IM1,L  !@var I,J,IP1,IM1,L  loop variables
      Real*8 :: VMASS,ALPH,PDT4,DT1,DT2,DT4,DT8,DT12,DT24
     *     ,FLUX,FLUXU,FLUXV
      REAL*8 FLUXU_N_S,FLUXV_N_S
      REAL*8 FLUXU_SW_NE,FLUXV_SW_NE
      REAL*8 FLUXU_SE_NW,FLUXV_SE_NW
      REAL*8 :: ASDU(IM,GRID%J_STRT_SKP:GRID%J_STOP_HALO,LM-1)
      Real*8 :: UP(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *          VP(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     *          USVS(IM,LM),VSVS(IM,LM),USVN(IM,LM),VSVN(IM,LM)
c pole fix variables
      integer :: hemi,jpo,jns,jv,jvs,jvn,ipole
      real*8 :: utmp,vtmp,wts
      real*8, dimension(im) :: dmt
C****
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_0STG, J_1STG
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO=J_0H,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      IF(MODD5K.LT.MRCH) CALL DIAG5F (U,V)

      DT2=DT1/2.
      DT4=DT1/4.
      DT8=DT1/8.
      DT12=DT1/12.
      DT24=DT1/24.
C****
C**** SCALE UT AND VT WHICH MAY THEN BE PHYSICALLY INTERPRETED AS
C**** MOMENTUM COMPONENTS
C****
      DO 110 J=J_0S,J_1
      I=IM
      DO 110 IP1=1,IM
      Do L=1,LM
         VMASS = .5*((MBEFOR(L,I,J-1)+MBEFOR(L,Ip1,J-1))*DXYN(J-1) +
     +               (MBEFOR(L,I,J  )+MBEFOR(L,Ip1,J  ))*DXYS(J))
         UT(I,J,L) = UT(I,J,L)*VMASS
         VT(I,J,L) = VT(I,J,L)*VMASS  ;  EndDo
  110 I = Ip1
      DUT(:,J_0S:J_1,:) = 0
      DVT(:,J_0S:J_1,:) = 0
C****
C**** BEGINNING OF LAYER LOOP
C****
!**** Use interpolated velocity values at poles
      Call HALO_UPDATE (GRID, U, From=NORTH)
      Call HALO_UPDATE (GRID, V, From=NORTH)
!     Call HALO_UPDATE (GRID, U, From=NORTH3)
!     Call HALO_UPDATE (GRID, V, From=NORTH3)
      If (J_0STG==2)  Then
         USVS(:,:) = U(:,2,:)
         VSVS(:,:) = V(:,2,:)
         U(:,2,:) = POLWT*U(:,2,:) + (1-POLWT)*U(:,3,:)
         V(:,2,:) = POLWT*V(:,2,:) + (1-POLWT)*V(:,3,:)  ;  EndIf

      Call HALO_UPDATE (GRID, U, From=SOUTH)
      Call HALO_UPDATE (GRID, V, From=SOUTH)
!     Call HALO_UPDATE (GRID, U, From=SOUTHJMm1)
!     Call HALO_UPDATE (GRID, V, From=SOUTHJMm1)
      If (J_1STG==JM)  Then
         USVN(:,:) = U(:,JM,:)
         VSVN(:,:) = V(:,JM,:)
         U(:,JM,:) = POLWT*U(:,JM,:) + (1-POLWT)*U(:,JM-1,:)
         V(:,JM,:) = POLWT*V(:,JM,:) + (1-POLWT)*V(:,JM-1,:)  ;  EndIf
      
      Call HALO_UPDATE (GRID, U, From=NORTH)
      Call HALO_UPDATE (GRID, V, From=NORTH)
      DO 300 L=1,LM

C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****
      I=IM
      DO 230 IP1=1,IM
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX
      DO 210 J=J_0STG,J_1STG
      FLUX=DT12*(PU(IP1,J,L)+PU(IP1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
      FLUXU=FLUX*(U(IP1,J,L)+U(I,J,L))
      DUT(IP1,J,L)=DUT(IP1,J,L)+FLUXU
      DUT(I,J,L)  =DUT(I,J,L)  -FLUXU
      FLUXV=FLUX*(V(IP1,J,L)+V(I,J,L))
      DVT(IP1,J,L)=DVT(IP1,J,L)+FLUXV
  210 DVT(I,J,L)  =DVT(I,J,L)  -FLUXV

      IF (haveLatitude(grid, J=1)) THEN ! No southern contribution off boundary.
        FLUXU_N_S  =0.
        FLUXV_N_S  =0.
        FLUXU_SW_NE=0.
        FLUXV_SW_NE=0.
        FLUXU_SE_NW=0.
        FLUXV_SE_NW=0.
      ELSE ! from lower boundary
        J=J_0-1
        FLUX=DT12*(PV(I,J,L)+PV(IP1,J,L)+PV(I,J+1,L)+PV(IP1,J+1,L))
        FLUXU_N_S=FLUX*(U(I,J,L)+U(I,J+1,L))
        FLUXV_N_S=FLUX*(V(I,J,L)+V(I,J+1,L))
        FLUX=DT24*(PU(IP1,J,L)+PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SW_NE=FLUX*(U(IP1,J+1,L)+U(I,J,L))
        FLUXV_SW_NE=FLUX*(V(IP1,J+1,L)+V(I,J,L))
        FLUX=DT24*(-PU(IP1,J,L)-PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SE_NW=FLUX*(U(I,J+1,L)+U(IP1,J,L))
        FLUXV_SE_NW=FLUX*(V(I,J+1,L)+V(IP1,J,L))
      END IF

      DO J=J_0S,J_1S

        DUT(I,J,L)   = DUT(I,J,L)   + FLUXU_N_S
        DVT(I,J,L)   = DVT(I,J,L)   + FLUXV_N_S
        DUT(IP1,J,L) = DUT(IP1,J,L) + FLUXU_SW_NE
        DVT(IP1,J,L) = DVT(IP1,J,L) + FLUXV_SW_NE
        DUT(I,J,L)   = DUT(I,J,L)   + FLUXU_SE_NW
        DVT(I,J,L)   = DVT(I,J,L)   + FLUXV_SE_NW

C**** CONTRIBUTION FROM THE SOUTH-NORTH MASS FLUX

        FLUX=DT12*(PV(I,J,L)+PV(IP1,J,L)+PV(I,J+1,L)+PV(IP1,J+1,L))
        FLUXU_N_S=FLUX*(U(I,J,L)+U(I,J+1,L))
        FLUXV_N_S=FLUX*(V(I,J,L)+V(I,J+1,L))
      
        DUT(I,J,L)=DUT(I,J,L)-FLUXU_N_S
        DVT(I,J,L)=DVT(I,J,L)-FLUXV_N_S

C**** CONTRIBUTION FROM THE SOUTHWEST-NORTHEAST MASS FLUX
        FLUX=DT24*(PU(IP1,J,L)+PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SW_NE=FLUX*(U(IP1,J+1,L)+U(I,J,L))
        FLUXV_SW_NE=FLUX*(V(IP1,J+1,L)+V(I,J,L))

        DUT(I,J,L) = DUT(I,J,L) - FLUXU_SW_NE
        DVT(I,J,L) = DVT(I,J,L) - FLUXV_SW_NE

C**** CONTRIBUTION FROM THE SOUTHEAST-NORTHWEST MASS FLUX
        FLUX=DT24*(-PU(IP1,J,L)-PU(I,J,L)+PV(IP1,J,L)+PV(IP1,J+1,L))
        FLUXU_SE_NW=FLUX*(U(I,J+1,L)+U(IP1,J,L))
        FLUXV_SE_NW=FLUX*(V(I,J+1,L)+V(IP1,J,L))

        DUT(IP1,J,L) = DUT(IP1,J,L) - FLUXU_SE_NW
        DVT(IP1,J,L) = DVT(IP1,J,L) - FLUXV_SE_NW

      END DO

      IF (HAVE_NORTH_POLE) THEN
        J=JM
        DUT(I,J,L) = DUT(I,J,L) + FLUXU_N_S
        DVT(I,J,L) = DVT(I,J,L) + FLUXV_N_S

        DUT(IP1,J,L) = DUT(IP1,J,L) + FLUXU_SW_NE
        DVT(IP1,J,L) = DVT(IP1,J,L) + FLUXV_SW_NE

        DUT(I,J,L) = DUT(I,J,L) + FLUXU_SE_NW
        DVT(I,J,L) = DVT(I,J,L) + FLUXV_SE_NW
      END IF

  230 I=IP1

  300 CONTINUE

!**** Restore uninterpolated values of U,V at poles   
      If (J_0STG==2)  Then
         U(:,2,:) = USVS(:,:)
         V(:,2,:) = VSVS(:,:)  ;  EndIf
      If (J_1STG==JM)  Then
         U(:,JM,:) = USVN(:,:)
         V(:,JM,:) = VSVN(:,:)  ;  EndIf
      
      if(do_polefix.eq.1) then
c Horizontal advection for the polar row is performed upon
c x-y momentum rather than spherical-coordinate momentum,
c eliminating the need for the metric term (and for the latter
c to be exactly consistent with the advection scheme).
c Southwest-Northeast and Southeast-Northwest corner fluxes
c are discarded since they produce erroneous tendencies at the pole
c in models with a polar half-box.
c Explicit cross-polar advection will be ignored until issues with
c corner fluxes and spherical geometry can be resolved.
       
!**** Polar velocities need values next to pole
!     Call HALO_UPDATE (GRID, U, From=NORTH)
!     Call HALO_UPDATE (GRID, V, From=NORTH)
!     Call HALO_UPDATE (GRID, U, From=NORTH3)
!     Call HALO_UPDATE (GRID, V, From=NORTH3)
      Call HALO_UPDATE (GRID, U, From=SOUTH)
      Call HALO_UPDATE (GRID, V, From=SOUTH)
!     Call HALO_UPDATE (GRID, U, From=SOUTHJMm1)
!     Call HALO_UPDATE (GRID, V, From=SOUTHJMm1)

        do ipole=1,2
          if(J_0STG==2 .and. ipole.eq.1) then
            hemi = -1
            jpo = 1
            jns = jpo + 1
            jv = 2           ! why not staggered grid
            jvs = 2          ! jvs is the southernmost velocity row
            jvn = jvs + 1       ! jvs is the northernmost velocity row
            wts = polwt
          else if(J_1STG==JM .and. ipole.eq.2) then
            hemi = +1
            jpo = JM
            jns = jpo - 1
            jv = JM            ! why not staggered grid
            jvs = jv - 1
            jvn = jvs + 1
            wts = 1.-polwt
          else
            cycle
          endif
c loop over layers
      do l=1,lm
c
c Copy u,v into temporary storage and transform u,v to x-y coordinates
c
      do j=jvs,jvn
         do i=1,im
            UP(i,j,l) = cosi(i)*U(I,J,L) - hemi*sini(i)*V(I,J,L)
            VP(i,j,l) = cosi(i)*V(I,J,L) + hemi*sini(i)*U(I,J,L)
         enddo
      enddo
c
c interpolate polar velocities to the appropriate latitude
c
      UP(:,jv,l) = wts*UP(:,jvs,l) + (1.-wts)*UP(:,jvn,l)
      VP(:,jv,l) = wts*VP(:,jvs,l) + (1.-wts)*VP(:,jvn,l)

c
c Compute advective tendencies of xy momentum in the polar rows
c
      dmt(:) = 0.
      dut(:,jv,l) = 0.
      dvt(:,jv,l) = 0.
c
      i = im
      do ip1=1,im
C**** CONTRIBUTION FROM THE WEST-EAST MASS FLUX
      FLUX=DT8*(PU(IP1,JPO,L)+PU(I,JPO,L)+PU(IP1,JNS,L)+PU(I,JNS,L))
      FLUXU=FLUX*(UP(IP1,JV,L)+UP(I,JV,L))
      DUT(IP1,JV,L)=DUT(IP1,JV,L)+FLUXU
      DUT(I,JV,L)  =DUT(I,JV,L)  -FLUXU
      FLUXV=FLUX*(VP(IP1,JV,L)+VP(I,JV,L))
      DVT(IP1,JV,L)=DVT(IP1,JV,L)+FLUXV
      DVT(I,JV,L)  =DVT(I,JV,L)  -FLUXV
      DMT(IP1)=DMT(IP1)+FLUX+FLUX
      DMT(I)  =DMT(I)  -FLUX-FLUX
C**** CONTRIBUTION FROM THE SOUTH-NORTH MASS FLUX
      FLUX=DT8*(PV(I,JVS,L)+PV(IP1,JVS,L)+PV(I,JVN,L)+PV(IP1,JVN,L))
      FLUX=FLUX*HEMI
      DUT(I,JV,L)=DUT(I,JV,L)+FLUX*(UP(I,JVS,L)+UP(I,JVN,L))
      DVT(I,JV,L)=DVT(I,JV,L)+FLUX*(VP(I,JVS,L)+VP(I,JVN,L))
      DMT(I)=DMT(I)+FLUX+FLUX
      i = ip1
      enddo ! i
c
c correct for the too-large dxyv in the polar row
c and convert dut,dvt from xy to polar coordinates
c
      do i=1,im
         dut(i,jv,l) = dut(i,jv,l) +
     &        (acor-1d0)*(dut(i,jv,l)-dmt(i)*UP(i,jv,l))
         dvt(i,jv,l) = dvt(i,jv,l) +
     &        (acor-1d0)*(dvt(i,jv,l)-dmt(i)*VP(i,jv,l))
         utmp = dut(i,jv,l)
         vtmp = dvt(i,jv,l)
         dut(i,jv,l) = cosi(i)*utmp+hemi*sini(i)*vtmp
         dvt(i,jv,l) = cosi(i)*vtmp-hemi*sini(i)*utmp
      enddo
      enddo ! loop over layers
      enddo ! loop over poles
      endif

C****
C**** VERTICAL ADVECTION OF MOMENTUM
C****
C     DO 310 L=1,LM-1
C     DO 310 J=2,JM
C     I=IM
C     DO 310 IP1=1,IM
C     SDU=DT2*((SD(I,JJ(J-1),L)+SD(IP1,JJ(J-1),L))*RAVPN(J-1)+
C    *  (SD(I,JJ(J),L)+SD(IP1,JJ(J),L))*RAVPS(J))
C     DUT(I,J,L)  =DUT(I,J,L)  +SDU*(U(I,J,L)+U(I,J,L+1))
C     DUT(I,J,L+1)=DUT(I,J,L+1)-SDU*(U(I,J,L)+U(I,J,L+1))
C     DVT(I,J,L)  =DVT(I,J,L)  +SDU*(V(I,J,L)+V(I,J,L+1))
C     DVT(I,J,L+1)=DVT(I,J,L+1)-SDU*(V(I,J,L)+V(I,J,L+1))
C 310 I=IP1
      CALL HALO_UPDATE(GRID,SD,FROM=SOUTH)
      DO L=1,LM-1
      DO J=J_0S,J_1
         DO I=1,IM-1
           ASDU(I,J,L)=DT2*(
     *                (SD(I,J-1,L)+SD(I+1,J-1,L))
     *                 *RAVPN(J-1)+
     *                (SD(I,J,L)+SD(I+1,J,L))
     *                 *RAVPS(J) )
         END DO
           ASDU(IM,J,L)=DT2*(
     *                (SD(IM,J-1,L)+SD(1,J-1,L))
     *                  *RAVPN(J-1)+
     *                (SD(IM,J,L)+SD(1,J,L))
     *                  *RAVPS(J) )
      END DO
      END DO

      L=1
      DO J=J_0S,J_1
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
      END DO
      DO L=2,LM-1
        DO J=J_0S,J_1
         DUT(:,J,L)  =DUT(:,J,L)  -ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DUT(:,J,L)  =DUT(:,J,L)  +ASDU(:,J,L)  *(U(:,J,L)+U(:,J,L+1))
         DVT(:,J,L)  =DVT(:,J,L)  -ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
         DVT(:,J,L)  =DVT(:,J,L)  +ASDU(:,J,L)  *(V(:,J,L)+V(:,J,L+1))
        END DO
      END DO
      L=LM
      DO J=J_0S,J_1
         DUT(:,J,L)=DUT(:,J,L)-ASDU(:,J,L-1)*(U(:,J,L-1)+U(:,J,L))
         DVT(:,J,L)=DVT(:,J,L)-ASDU(:,J,L-1)*(V(:,J,L-1)+V(:,J,L))
      END DO
C**** CALL DIAGNOSTICS
         IF(MODD5K.LT.MRCH) CALL DIAG5D (4,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAGCD (grid,1,U,V,DUT,DVT,DT1)
      DO L=1,LM
      DO J=J_0S,J_1
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
      DO L=1,LM
        IM1=IM
        DO I=1,IM
C         FD(I,1)=FCOR(1)*2.  -.5*(SPA(IM1,1,L)+SPA(I,1,L))*DXV(2)
C         FD(I,JM)=FCOR(JM)*2.+.5*(SPA(IM1,JM,L)+SPA(I,JM,L))*DXV(JM)
C****     Set the Coriolis term to zero at the Poles:
          IF(haveLatitude(grid,J=1))
     *    FD(I,1)= -.5*(SPA(IM1,1,L)+
     *                 SPA(I,1,L))*
     *                 DXV(2)
          IF(haveLatitude(grid,J=JM))
     *    FD(I,JM)=  .5*(SPA(IM1,JM,L)+
     *                 SPA(I,JM,L))*
     *                 DXV(JM)
          IM1=I
        END DO

        DO J=J_0S,J_1S
          IM1=IM
          DO I=1,IM
            FD(I,J)=FCOR(J)+.25*(SPA(IM1,J,L)+SPA(I,J,L))
     *             *(DXV(J)-DXV(J+1))
            IM1=I
          END DO
        END DO
        CALL HALO_UPDATE(GRID,FD,FROM=SOUTH)
        DO J=J_0S,J_1
          IM1=IM
          DO I=1,IM
            PDT4 = DT8*(MMEAN(L,I,J-1)+MMEAN(L,I,J))
            ALPH = PDT4*(FD(I,J)+FD(I,J-1))
            DUT(I,J,L)=DUT(I,J,L)+ALPH*V(I,J,L)
            DUT(IM1,J,L)=DUT(IM1,J,L)+ALPH*V(IM1,J,L)
            DVT(I,J,L)=DVT(I,J,L)-ALPH*U(I,J,L)
            DVT(IM1,J,L)=DVT(IM1,J,L)-ALPH*U(IM1,J,L)
            IM1=I
          END DO
        END DO
      END DO

      if(do_polefix.eq.1) then
c apply the full coriolis force at the pole and ignore the metric term
c which has already been included in advective form
         do ipole=1,2
            if(haveLatitude(grid,J=2) .and. ipole.eq.1) then
               jpo = 1
               jns = jpo + 1
               j = 2
            else if(haveLatitude(grid,J=JM) .and. ipole.eq.2) then
               jpo = JM
               jns = jpo - 1
               j = JM 
            else
               cycle
            endif
            do l=1,lm
               dut(:,j,l) = 0.
               dvt(:,j,l) = 0.
               im1=im
               do i=1,im
                  PDT4 = DT8*(MMEAN(L,i,jpo)+MMEAN(L,i,jns))
                  ALPH = PDT4*(2*fcor(jpo) + fcor(jns))
                  dut(i  ,j,l)=dut(i  ,j,l)+alph*v(i  ,j,l)
                  dut(im1,j,l)=dut(im1,j,l)+alph*v(im1,j,l)
                  dvt(i  ,j,l)=dvt(i  ,j,l)-alph*u(i  ,j,l)
                  dvt(im1,j,l)=dvt(im1,j,l)-alph*u(im1,j,l)
                  im1=i
               enddo
            enddo
         enddo
      endif

C**** CALL DIAGNOSTICS
         IF(MODD5K.LT.MRCH) CALL DIAG5D (5,MRCH,DUT,DVT)
         IF(MRCH.GT.0) CALL DIAGCD (grid,2,U,V,DUT,DVT,DT1)
C****
C**** ADD CORIOLIS FORCE INCREMENTS TO UT AND VT
C**** AND UNDO SCALING PERFORMED AT BEGINNING OF ADVECV
C****
      DO 610 J=J_0S,J_1
      I=IM
      DO 610 IP1=1,IM
      Do L=1,LM
         VMASS = .5*((MAFTER(L,I,J-1)+MAFTER(L,Ip1,J-1))*DXYN(J-1) +
     +               (MAFTER(L,I,J  )+MAFTER(L,Ip1,J  ))*DXYS(J))
         UT(I,J,L) = (UT(I,J,L) + DUT(I,J,L)) / VMASS
         VT(I,J,L) = (VT(I,J,L) + DVT(I,J,L)) / VMASS  ;  EndDo
  610 I = Ip1
      DUT(:,J_0S:J_1,:) = 0
      DVT(:,J_0S:J_1,:) = 0
C
      RETURN
      END SUBROUTINE ADVECV

      end module MOMENTS
