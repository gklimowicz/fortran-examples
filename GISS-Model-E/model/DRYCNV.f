#include "rundeck_opts.h"

      SUBROUTINE ATM_DIFFUS(LBASE_MIN,LBASE_MAX,DTIME)
!@sum  ATM_DIFFUS(DRYCNV) mixes air caused by dry convection.
!@+    this version checks base layers lbase_min to lbase_max.
!@auth Original Development Team
      USE CONSTANT, only : lhe,sha,deltx
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v,q,t,pk,pdsig,plij,pedn
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds
      USE DOMAIN_DECOMP_ATM, only : halo_update, checksum
      USE DOMAIN_DECOMP_ATM, only : halo_update_column, checksum_column
      USE DOMAIN_DECOMP_ATM, only : NORTH, SOUTH
      USE GEOM, only : imaxj, kmaxj, ravj, idij, idjj
      USE QUSDEF, only : nmom,zmoms,xymoms
      USE SOMTQ_COM, only : tmom,qmom
      USE DIAG_COM, only : jl_trbhr,jl_damdc,jl_trbdlht
#ifdef TRACERS_ON
      USE TRACER_COM, only: TRM,TRMOM,NTM
      USE TRDIAG_COM, only: JLNT_TURB
#endif
      USE PBLCOM, only : dclev,w2gcm,w2_l1
      IMPLICIT NONE

      integer, intent(in) :: LBASE_MIN,LBASE_MAX
      real*8, intent(in) :: dtime  ! dummy variable
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) ::
     &        UT,VT
      REAL*8, DIMENSION(LM) :: DP
      INTEGER, DIMENSION(IM) :: IDI,IDJ    !@var ID
      REAL*8, DIMENSION(IM) :: RA !@var
      REAL*8, DIMENSION(IM) :: UMS,VMS !@var
      INTEGER I,J,L,K,IMAX,KMAX,IM1,LMAX,LMIN,N
C
      REAL*8  UKP1(IM,LM), VKP1(IM,LM), UKPJM(IM,LM),VKPJM(IM,LM)

      REAL*8  UKM(4,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM)
      REAL*8  VKM(4,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM)
      INTEGER  LRANG(2,IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
C
      REAL*8, DIMENSION(NMOM) :: TMOMS,QMOMS
      REAL*8 PKMS,QMS,TVMS,THETA,RDP,THM

#ifdef TRACERS_ON
      REAL*8, DIMENSION(NMOM,NTM) :: TRMOMS
      REAL*8, DIMENSION(     NTM) :: TRMS
      REAL*8 SDPL,BYSDPL
#endif


      INTEGER ::  J_1, J_0, I_1, I_0
      INTEGER ::  J_1H, J_0H
      INTEGER ::  J_1S, J_0S
      INTEGER ::  J_1STG, J_0STG
      LOGICAL ::  HAVE_NORTH_POLE, HAVE_SOUTH_POLE

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP


      if(LBASE_MAX.GE.LM) call stop_model('DRYCNV: LBASE_MAX.GE.LM',255)

      ! update w2gcm at 1st GCM layer
      w2gcm=0.d0
      do j=j_0,j_1
      do i=i_0,i_1
         w2gcm(1,i,j)=w2_l1(i,j)
      end do
      end do

C****
C**** Update north halos for arrays U and V
C****
      CALL HALO_UPDATE(grid, U, from=NORTH)
      CALL HALO_UPDATE(grid, V, FROM=NORTH)

C**** LOAD U,V INTO UT,VT.  UT,VT WILL BE FIXED DURING DRY CONVECTION
C****   WHILE U,V WILL BE UPDATED.

      UT=U ; VT=V
C**** OUTSIDE LOOPS OVER J AND I
      JLOOP: DO J=J_0,J_1

      IMAX=IMAXJ(J)
      KMAX=KMAXJ(J)
C****
C**** MAIN LOOP
C****
      IM1=IM
      ILOOP: DO I=I_0,IMAX
         DO K=1,KMAX
            RA(K)=RAVJ(K,J)
            IDI(K)=IDIJ(K,I,J)
            IDJ(K)=IDJJ(K,J)
         END DO
C
         LRANG(1,I,J)=-1
         LRANG(2,I,J)=-2
C
      LMAX=LBASE_MIN-1
      lbase_loop: do while(lmax.lt.lbase_max)
      LMIN=LMAX+1
      LMAX=LMIN
      IF (T(I,J,LMIN)*(1.+Q(I,J,LMIN)*deltx).LE.
     *   T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*deltx)) cycle lbase_loop
C**** MIX HEAT AND MOISTURE THROUGHOUT THE UNSTABLE LAYERS
C**** MIX THROUGH TWO LOWER LAYERS
      DP(LMIN)=PDSIG(LMIN,I,J)
      DP(LMIN+1)=PDSIG(LMIN+1,I,J)
      PKMS=PK(LMIN,I,J)*DP(LMIN)+PK(LMIN+1,I,J)*DP(LMIN+1)
      TVMS=T(I,J,LMIN)*(1.+Q(I,J,LMIN)*deltx)*(PK(LMIN,I,J)*DP(LMIN))
     *    +T(I,J,LMIN+1)*(1.+Q(I,J,LMIN+1)*deltx)
     *                                  *(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMS=Q(I,J,LMIN)*DP(LMIN)+Q(I,J,LMIN+1)*DP(LMIN+1)
C**** sum moments to mix over unstable layers
      TMOMS(XYMOMS) =
     &     TMOM(XYMOMS,I,J,LMIN  )*(PK(LMIN  ,I,J)*DP(LMIN  ))  +
     &     TMOM(XYMOMS,I,J,LMIN+1)*(PK(LMIN+1,I,J)*DP(LMIN+1))
      QMOMS(XYMOMS) =
     &     QMOM(XYMOMS,I,J,LMIN  )*(DP(LMIN  ))  +
     &     QMOM(XYMOMS,I,J,LMIN+1)*(DP(LMIN+1))
#ifdef TRACERS_ON
      TRMS(:) = TRM(I,J,LMIN,:)+TRM(I,J,LMIN+1,:)
      TRMOMS(XYMOMS,:) =
     &     TRMOM(XYMOMS,I,J,LMIN,:)+TRMOM(XYMOMS,I,J,LMIN+1,:)
#endif
      IF (LMIN+1.GE.LM) GO TO 150
      THETA=TVMS/PKMS
C**** MIX THROUGH SUBSEQUENT UNSTABLE LAYERS
      DO L=LMIN+2,LM
        IF (THETA.LT.T(I,J,L)*(1.+Q(I,J,L)*deltx)) GO TO 160
        DP(L)=PDSIG(L,I,J)
        PKMS=PKMS+(PK(L,I,J)*DP(L))
        QMS=QMS+Q(I,J,L)*DP(L)
        TVMS=TVMS+T(I,J,L)*(1.+Q(I,J,L)*deltx)*(PK(L,I,J)*DP(L))
        TMOMS(XYMOMS) = TMOMS(XYMOMS) +
     &       TMOM(XYMOMS,I,J,L)*(PK(L,I,J)*DP(L))
        QMOMS(XYMOMS) = QMOMS(XYMOMS) +
     &       QMOM(XYMOMS,I,J,L)*DP(L)
        THETA=TVMS/PKMS
#ifdef TRACERS_ON
      TRMS(:) = TRMS(:) + TRM(I,J,L,:)
      TRMOMS(XYMOMS,:) = TRMOMS(XYMOMS,:) + TRMOM(XYMOMS,I,J,L,:)
#endif
      END DO
  150 L=LM+1
  160 LMAX=L-1
      RDP=1./(PEDN(LMIN,I,J)-PEDN(LMAX+1,I,J))
      QMS=QMS*RDP
      THM=TVMS/(PKMS*(1.+QMS*deltx))
#ifdef TRACERS_ON
        SDPL = 0.d0
        DO L=LMIN,LMAX
          SDPL = SDPL+DP(L)
        ENDDO
        BYSDPL = 1.D0/SDPL
#endif
      DO L=LMIN,LMAX
      CALL INC_AJL(I,J,L,JL_TRBHR,(THM-T(I,J,L))*PK(L,I,J)*PLIJ(L,I,J))
      CALL INC_AJL(I,J,L,JL_TRBDLHT,(QMS-Q(I,J,L))*PDSIG(L,I,J)*LHE/SHA)
      T(I,J,L)=THM
      TMOM(XYMOMS,I,J,L)=TMOMS(XYMOMS)/PKMS
      TMOM(ZMOMS,I,J,L)=0.
      Q(I,J,L)=QMS
      QMOM(XYMOMS,I,J,L)=QMOMS(XYMOMS)*RDP
      QMOM(ZMOMS,I,J,L)=0.
#ifdef TRACERS_ON
      DO N=1,NTM
        CALL INC_TAJLN(I,J,L,JLNT_TURB,N,TRMS(N)*(DP(L)*BYSDPL)-TRM(I,J
     *       ,L,N))
      END DO
      TRM(I,J,L,:) = TRMS(:)*(DP(L)*BYSDPL)
      TRMOM(XYMOMS,I,J,L,:) = TRMOMS(XYMOMS,:)*(DP(L)*BYSDPL)
      TRMOM(ZMOMS,I,J,L,:) = 0.
#endif
      END DO
C**** MIX MOMENTUM THROUGHOUT UNSTABLE LAYERS
      UMS(1:KMAX)=0.
      VMS(1:KMAX)=0.
      DO L=LMIN,LMAX
         DO K=1,KMAX
            UMS(K)=UMS(K)+UT(IDI(K),IDJ(K),L)*DP(L)
            VMS(K)=VMS(K)+VT(IDI(K),IDJ(K),L)*DP(L)
         ENDDO
      ENDDO
      UMS(1:KMAX)=UMS(1:KMAX)*RDP
      VMS(1:KMAX)=VMS(1:KMAX)*RDP
      LRANG(1,I,J)=LMIN
      LRANG(2,I,J)=LMAX
c     DO L=LMIN,LMAX
c        DO K=1,KMAX
c           U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)
c    &           +(UMS(K)-UT(IDI(K),IDJ(K),L))*RA(K)
c           V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)
c    &           +(VMS(K)-VT(IDI(K),IDJ(K),L))*RA(K)
c           CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,
c    &           (UMS(K)-UT(IDI(K),IDJ(K),L))*PLIJ(L,I,J)*RA(K))
c        ENDDO
c     ENDDO
      DO L=LMIN,LMAX
         IF(J.EQ.1)  THEN
            DO K=1,KMAX
               UKP1(K,L)=(UMS(K)-UT(IDI(K),IDJ(K),L))
               VKP1(K,L)=(VMS(K)-VT(IDI(K),IDJ(K),L))
            END DO
         ELSE IF(J.EQ.JM)  THEN
            DO K=1,KMAX
               UKPJM(K,L)=(UMS(K)-UT(IDI(K),IDJ(K),L))
               VKPJM(K,L)=(VMS(K)-VT(IDI(K),IDJ(K),L))
            END DO
         ELSE
            DO K=1,KMAX
               UKM(K,I,J,L)=(UMS(K)-UT(IDI(K),IDJ(K),L))
               VKM(K,I,J,L)=(VMS(K)-VT(IDI(K),IDJ(K),L))
            END DO
         END IF
      ENDDO
C
      enddo lbase_loop
C**** ACCUMULATE BOUNDARY LAYER DIAGNOSTICS
      if(lbase_min.eq.1) then ! was called from surfce
         DCLEV(I,J)=LMAX
      endif
      IM1=I
      ENDDO ILOOP
      ENDDO JLOOP
C
C     NOW REALLY UPDATE THE MODEL WINDS
C
C***  ...first update halo (J_0-1 values) of UKM,VKM, and PLIJ.
      call halo_update_column(grid, UKM, from=SOUTH)
      call halo_update_column(grid, VKM, from=SOUTH)
      call halo_update_column(grid,PLIJ, from=SOUTH)
      call halo_update_column(grid,LRANG,from=SOUTH)

#ifdef SCM
C
      i=I_TARG
      j=J_TARG
        KMAX=KMAXJ(J)
        DO K=1,KMAX
           IDJ(K)=IDJJ(K,J)
           RA(K) =RAVJ(K,J)
        END DO
c       DO I=1,IM
          LMIN=LRANG(1,I,J)
          LMAX=LRANG(2,I,J)
          DO L=LMIN,LMAX
          DO K=1,KMAX
            IDI(K)=IDIJ(K,I,J)
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)*RA(K)
            CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,UKM(K,I,J,L)*PLIJ(L,I
     *           ,J)*RA(K))
          END DO ; END DO
c       END DO
#else
      IF (HAVE_SOUTH_POLE) then
       J=1
       DO K=1,KMAXJ(J)
         IDI(K)=IDIJ(K,1,J)
         IDJ(K)=IDJJ(K,J)
         RA(K) =RAVJ(K,J)
       END DO
       LMIN=LRANG(1,1,J)
       LMAX=LRANG(2,1,J)
       DO L=LMIN,LMAX
        DO K=1,KMAXJ(J)
          U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKP1(K,L)*RA(K)
          V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKP1(K,L)*RA(K)
          CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,UKP1(K,L)*PLIJ(L,1,J)
     *         *RA(K))
       END DO ; END DO
      ElSE ! need contribution from southern neighbor
        J=J_0-1
        KMAX=KMAXJ(J)
        DO K=1,KMAX
           IDJ(K)=IDJJ(K,J)
           RA(K) =RAVJ(K,J)
        END DO
        DO I=1,IM
          LMIN=LRANG(1,I,J)
          LMAX=LRANG(2,I,J)
          DO L=LMIN,LMAX
          DO K=1,KMAX
            IF (IDJ(K) == J_0) THEN
              IDI(K)=IDIJ(K,I,J)
              U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)*RA(K)
              V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)*RA(K)
              CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,UKM(K,I,J,L)*PLIJ(L
     *             ,I,J)*RA(K))
            END IF
          END DO ; END DO
        END DO
      END IF   !END SOUTH POLE
C
      DO J=J_0S, J_1-1       !J_1S computed below
        KMAX=KMAXJ(J)
        DO K=1,KMAX
           IDJ(K)=IDJJ(K,J)
           RA(K) =RAVJ(K,J)
        END DO
        DO I=1,IM
          LMIN=LRANG(1,I,J)
          LMAX=LRANG(2,I,J)
          DO L=LMIN,LMAX
          DO K=1,KMAX
            IDI(K)=IDIJ(K,I,J)
            U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)*RA(K)
            V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)*RA(K)
            CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,UKM(K,I,J,L)*PLIJ(L,I
     *           ,J)*RA(K))
          END DO ; END DO
        END DO
      END DO
C
      IF (HAVE_NORTH_POLE) THEN
       J=JM
       KMAX=KMAXJ(J)
       DO K=1,KMAX
         IDI(K)=IDIJ(K,1,J)
         IDJ(K)=IDJJ(K,J)
         RA(K) =RAVJ(K,J)
       END DO
       LMIN=LRANG(1,1,J)
       LMAX=LRANG(2,1,J)
       DO L=LMIN,LMAX
        DO K=1,KMAX
          U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKPJM(K,L)*RA(K)
          V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKPJM(K,L)*RA(K)
          CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,UKPJM(K,L)*PLIJ(L,1,J)
     *         *RA(K))
       END DO ; END DO

      ELSE
C**** First half of loop cycle for j=j_1 for internal blocks
        J=J_1

        KMAX=KMAXJ(J)
        DO K=1,2
           IDJ(K)=IDJJ(K,J)
           RA(K) =RAVJ(K,J)
        END DO
        DO I=1,IM
          LMIN=LRANG(1,I,J)
          LMAX=LRANG(2,I,J)
          DO L=LMIN,LMAX
            DO K=1,2
              IDI(K)=IDIJ(K,I,J)
              U(IDI(K),IDJ(K),L)=U(IDI(K),IDJ(K),L)+UKM(K,I,J,L)*RA(K)
              V(IDI(K),IDJ(K),L)=V(IDI(K),IDJ(K),L)+VKM(K,I,J,L)*RA(K)
              CALL INC_AJL(IDI(K),IDJ(K),L,JL_DAMDC,UKM(K,I,J,L)*PLIJ(L
     *             ,I,J)*RA(K))
            END DO
          END DO
        END DO
      ENDIF   !END NORTH POLE
#endif

C***

      RETURN
      END SUBROUTINE ATM_DIFFUS


      subroutine apply_fluxes_to_atm(dt)
!@sum applies earth fluxes to the first layer of the atmosphere
!@vers 2013/03/27
!@auth Original Development Team
      USE MODEL_COM, only : qcheck
#ifdef SCM
     *                      ,I_TARG,J_TARG
#endif
      USE RESOLUTION, only : im,jm
      USE ATM_COM, only : u,v,t,q,byMA,MA,pk
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds
      USE DOMAIN_DECOMP_ATM, only : halo_update,checksum
      USE DOMAIN_DECOMP_ATM, only : halo_update_column,checksum_column
      USE DOMAIN_DECOMP_ATM, only : NORTH, SOUTH
      USE GEOM, only : imaxj,kmaxj,ravj,idij,idjj,siniv,cosiv,axyp
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,trm,trmom,trname,t_qlimit
#ifdef TRACERS_WATER
     *     ,trw0,tr_wd_TYPE,nWATER
#endif
      USE FLUXES, only : trflux1
#endif
      USE FLUXES, only : dth1,dq1,uflux1,vflux1,qflux1
      implicit none
      REAL*8, PARAMETER :: qmin=1.d-12
      integer i,j,k,n
      real*8, intent(in) :: dt
      real*8 hemi,trmin
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &                       usave,vsave
      INTEGER :: J_0,J_1,J_0S,J_1S,J_0STG,J_1STG,I_0,I_1
      LOGICAL :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE       )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=j_0,j_1
        do i=i_0,imaxj(j)
          t(i,j,1) = t(i,j,1) + dth1(i,j)/pk(1,i,j)
          q(i,j,1) = q(i,j,1) + dq1(i,j)
        end do
      end do

#ifdef TRACERS_ON
      do n=1,ntm
        do j=j_0,j_1
          do i=i_0,imaxj(j)
            trm(i,j,1,n) = trm(i,j,1,n) + trflux1(i,j,n)*dt*axyp(i,j)
            trmin=0.d0
#ifdef TRACERS_WATER
            IF(tr_wd_TYPE(n).eq.nWATER) trmin =
     &         qmin*trw0(n)*MA(1,i,j)*axyp(i,j)
#endif
            if (t_qlimit(n).and.trm(i,j,1,n).lt.trmin) then
              if (qcheck) write(99,*) trname(n),I,J,' TR1:',
     *        trm(i,j,1,n),'->',trmin
              trm(i,j,1,n) = trmin
              trmom(:,i,j,1,n)=0.
            end if
          end do
        end do
      end do
#endif
c****
c**** add in surface friction to first layer wind
c****
      usave=u(:,:,1) ; vsave=v(:,:,1)

C   *....update halo (J_0-1 values) of  byMA, uflux1, vflux1.
      Call HALO_UPDATE(grid, uflux1, from=SOUTH)
      Call HALO_UPDATE(grid, vflux1, from=SOUTH)
      Call HALO_UPDATE(grid, byMA(1,:,:), from=SOUTH)


#ifdef SCM
cccc if SCM - update winds ????
      i = I_TARG
      j= J_TARG
      do k=1,2
         u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *          ravj(k,j)*uflux1(i,j)*dt*byMA(1,I,J)
         v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *          ravj(k,j)*vflux1(i,j)*dt*byMA(1,I,J)
      enddo
#else
c**** SOUTH POLE BOX
      if (HAVE_SOUTH_POLE) then
      j=1
        hemi=-1.
        do i=1,imaxj(j)
        do k=1,kmaxj(j)
          u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *     ravj(k,j)*(uflux1(i,j)*cosiv(k)+vflux1(i,j)*siniv(k)*hemi)
     *     *dt*byMA(1,I,J)
          v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *     ravj(k,j)*(vflux1(i,j)*cosiv(k)-uflux1(i,j)*siniv(k)*hemi)
     *     *dt*byMA(1,I,J)
        end do
        end do
      Else
        j = j_0-1

        do i=1,imaxj(j)
        do k=3,4
           If (idjj(k,j) == j_0) Then
              u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *             ravj(k,j)*uflux1(i,j)*dt*byMA(1,I,J)
              v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *             ravj(k,j)*vflux1(i,j)*dt*byMA(1,I,J)
           end if
        end do
        end do

      end if                    !SOUTH POLE

      IF (HAVE_NORTH_POLE) then
        j=jm
        hemi=1.
        do i=1,imaxj(j)
          do k=1,kmaxj(j)
            u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *       ravj(k,j)*(uflux1(i,j)*cosiv(k)+vflux1(i,j)*siniv(k)*hemi)
     *       *dt*byMA(1,I,J)
            v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *       ravj(k,j)*(vflux1(i,j)*cosiv(k)-uflux1(i,j)*siniv(k)*hemi)
     *       *dt*byMA(1,I,J)
          end do
        end do
      End If

c**** non polar boxes
      do j=J_0S,J_1-1
        do i=1,imaxj(j)
        do k=1,kmaxj(j)
          u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*uflux1(i,j)*dt*byMA(1,I,J)
          v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *           ravj(k,j)*vflux1(i,j)*dt*byMA(1,I,J)
        end do
        end do
      end do

C****For distr. parallelization: North-most lattitude of internal blocks.
C   *--->(First half of j=j_1 loop cycle for internal blocks --k=1,2)
      IF (HAVE_NORTH_POLE) then
c***        j=jm
c***        hemi=1.
c***        do i=1,imaxj(j)
c***          do k=1,kmaxj(j)
c***            u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
c***     *       ravj(k,j)*(uflux1(i,j)*cosiv(k)+vflux1(i,j)*siniv(k)*hemi)
c***     *       *dt*byMA(1,I,J)
c***            v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
c***     *       ravj(k,j)*(vflux1(i,j)*cosiv(k)-uflux1(i,j)*siniv(k)*hemi)
c***     *       *dt*byMA(1,I,J)
c***          end do
c***        end do
      Else
          j=j_1
          do i=1,imaxj(j)
            do k=1,2
              u(idij(k,i,j),idjj(k,j),1)=u(idij(k,i,j),idjj(k,j),1) -
     *               ravj(k,j)*uflux1(i,j)*dt*byMA(1,I,J)
              v(idij(k,i,j),idjj(k,j),1)=v(idij(k,i,j),idjj(k,j),1) -
     *               ravj(k,j)*vflux1(i,j)*dt*byMA(1,I,J)
            end do
          end do
        ENDIF   !.not. NORTH POLE
#endif

      return
      end subroutine apply_fluxes_to_atm
