#include "rundeck_opts.h"

      SUBROUTINE DYNAM2
!@vers 2013/04/02
      USE model_com, only : im,jm,lm,ls1,
     &     u,v,t,p,NIdyn,dt,DTsrc,NSTEP,mrch,ndaa
      USE DYNAMICS, only : pu,pv,sd, MUs,MVs,MWs
      USE DYNAMICS, only : gz,phi
      USE SOMTQ_COM, only : tmom,mz
      USE DOMAIN_DECOMP_1D, only : grid, GET, globalsum
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, NORTH, SOUTH
      USE ATMDYN, only : sdrag
      USE GEOM, only : imaxj
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     PT, PX, PSAVE, AM1, AM2, FPEU,FPEV, PPGF
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     UT,VT,TT,TZ,TZT,MMA,UX,VX

      REAL*8 DTFS,DTLF, DAMSUM
      INTEGER I,J,L   !@var I,J,L  loop variables
      INTEGER NS, MODDA

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_HALO = J_0H, J_STOP_HALO = J_1H,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)


      DTFS=DT*2./3.
      DTLF=2.*DT

      call halo_update(grid,p)
      call halo_update(grid,u)
      call halo_update(grid,v)
      call halo_update(grid,t)

      DO L=1,LM
         MUs(:,:,L) = 0.
         MVs(:,:,L) = 0.
         MWs(:,:,L) = 0.
      ENDDO

      DO L=1,LM
      DO J=J_0H,J_1H
      DO I=1,IM
        UX(I,J,L)  = U(I,J,L)
        UT(I,J,L)  = U(I,J,L)
        VX(I,J,L)  = V(I,J,L)
        VT(I,J,L)  = V(I,J,L)
      ENDDO
      ENDDO
      DO J=J_0,J_1
      DO I=1,IM
! copy z-moment of temperature into contiguous memory
         TZ(I,J,L)  = TMOM(MZ,I,J,L)
      ENDDO
      ENDDO
      ENDDO
      call halo_update(grid,tz,from=south)

C**** INITIAL FORWARD STEP, QX = Q + .667*DT*F(Q)
      MRCH=0
      CALL AIR_MASS_FLUX (U,V,P)
      CALL UPDATE_P (P,PX,DTFS)
      CALL UPDATE_UV (P,UX,VX,PX,U,V,P,T,TZ,DTFS)

C**** INITIAL BACKWARD STEP IS ODD, QT = Q + DT*F(QX)
      MRCH=-1
      CALL AIR_MASS_FLUX (UX,VX,PX)
      CALL UPDATE_P (P,PT,DT)
      CALL UPDATE_UV (P,UT,VT,PT,UX,VX,PX,T,TZ,DT)

      do ns=2,nidyn,2

C**** DIAGNOSTICS prep
        MRCH=2
        MODDA=MOD(NSTEP+NS-NIdyn+NDAA*NIdyn+2,NDAA*NIdyn+2)
        IF(MODDA.LT.MRCH) CALL DIAGA0

C**** EVEN LEAP FROG STEP, Q = Q + 2*DT*F(QT)
        MRCH=2
        PSAVE(:,:) = P(:,:)
        CALL AIR_MASS_FLUX (UT,VT,PT)
C**** ADVECT TEMPERATURE and accumulate mass fluxes
        do l=1,lm
          do j=j_0s-1,j_1
            do i=1,imaxj(j)
              tt(i,j,l) = t(i,j,l)
              tzt(i,j,l) = tz(i,j,l)
            enddo
          enddo
        enddo
        call calc_ma (psave,mma)
        CALL HALO_UPDATE(grid,PU, FROM=NORTH)
        CALL AADVT3 (MMA,T,TMOM, TZ, SD,PU,PV, DTLF)
        do l=1,lm
          do j=j_0s-1,j_1
            do i=1,imaxj(j)
              tt(i,j,l) = .5*(tt(i,j,l)+t(i,j,l))
              tzt(i,j,l) = .5*(tzt(i,j,l)+tz(i,j,l))
            enddo
          enddo
        enddo
        CALL UPDATE_P (PSAVE,P,DTLF)
        CALL UPDATE_UV (PSAVE,U,V,P,UT,VT,PT,TT,TZT,DTLF)

C**** DIAGNOSTICS
        IF(MODDA.LT.MRCH) THEN
          CALL CALC_AMPK(LS1-1)
          CALL DIAGA
          CALL DIAGB
          CALL EPFLUX (U,V,T,P)
        ENDIF

C**** ODD LEAP FROG STEP, QT = QT + 2*DT*F(Q)
        if(ns == nidyn) exit ! no need to further update the odd state
        MRCH=-2
        PSAVE(:,:) = PT(:,:)
        CALL AIR_MASS_FLUX (U,V,P)
        CALL UPDATE_P (PSAVE,PT,DTLF)
        CALL UPDATE_UV (PSAVE,UT,VT,PT,U,V,P,T,TZ,DTLF)

      enddo ! end loop over leapfrog steps

      CALL CALC_AMPK(LS1-1)
      GZ=PHI

c apply stratospheric drag
      CALL SDRAG (DTsrc)

c
c call gravity wave drag
c
      mrch = 2
      ut = u; vt = v
      CALL GWDRAG (P,U,V,UT,VT,T,TZ,DTsrc,.true.) ! strat
c      CALL VDIFF (P,U,V,UT,VT,T,DTsrc) ! strat

c      CALL CALC_AMPK(LS1-1)

c apply north-south filter to U and V once per physics timestep
      call conserv_amb_ext(u,am1) ! calculate ang. mom. before filter
      call fltry3(u,1d0) ! 2nd arg could be set using DT_YUfilter
      call fltry3(v,1d0) ! 2nd arg could be set using DT_YVfilter
      call conserv_amb_ext(u,am2) ! calculate ang. mom. after filter
      am2(:,j_0stg:j_1stg) = am1(:,j_0stg:j_1stg)-am2(:,j_0stg:j_1stg)
      if(have_south_pole) am2(:,1) = 0.
      call globalsum(grid,am2,damsum,all=.true.)
      call add_am_as_solidbody_rotation(u,damsum) ! maintain global ang. mom.
c apply east-west filter to U and V once per physics timestep
      CALL FLTRUV2(U,V,1d0)

c slp filter
c      CALL FILTER2

c adjust T/Q moments
      call tmom_topo_adjustments
      call qmom_topo_adjustments

      RETURN
      END SUBROUTINE DYNAM2

      SUBROUTINE AIR_MASS_FLUX (U,V,P)
!@sum  AIR_MASS_FLUX Calculates horizontal/vertical air mass fluxes
!@+    Input: U,V velocities, P pressure
!@+    Output: PIT  pressure tendency (mb m^2/s)
!@+            SD   sigma dot (mb m^2/s)
!@+            PU,PV horizontal mass fluxes (mb m^2/s)
!@+            CONV  horizontal mass convergence (mb m^2/s)
!@+            SPA
!@auth Original development team
      USE model_com, only : im,jm,lm,ls1,psfmpt,byim,dsig
      USE GEOM,  only : dxv=>dxlatv,dyp
      USE DYNAMICS, only : pit,sd,conv,pu,pv,spa
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      IMPLICIT NONE
C**** CONSTANT PRESSURE AT L=LS1 AND ABOVE, PU,PV CONTAIN DSIG
!@var U,V input velocities (m/s)
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,LM) ::
     &     U,V
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     &     P

      INTEGER I,J,L,jmin_pv,jmax_pv
      REAL*8 DXDSIG,DYDSIG
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     &     PATU,PDUM
      REAL*8, DIMENSION(IM,grid%j_strt_halo-1:grid%j_stop_halo) ::
     &     PATV

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      pdum(:,j_1) = p(:,j_1-1)
      call halo_update(grid,pdum,from=south)

C****
C**** BEGINNING OF LAYER LOOP
C****
C****
C**** COMPUTATION OF MASS FLUXES     P,T  PU     PRIMARY GRID ROW
C**** ARAKAWA SCHEME B               PV   U,V    SECONDARY GRID ROW
C****

      do j=max(2,j_0h),j_1s
        do i=1,im-1
          patu(i,j)=(p(i,j)+p(i+1,j))
        enddo
        i=im
          patu(i,j)=(p(i,j)+p(1  ,j))
      enddo
      do j=j_0s,j_1s+1
        do i=1,im
          patv(i,j)=(p(i,j)+p(i,j-1))
        enddo
      enddo
      if(j_0s.gt.2) then
        j = j_0s-1
        do i=1,im
          patv(i,j)=(p(i,j)+pdum(i,j))
        enddo
      endif
      PIT(:,J_0S-1:J_1) = 0.
      do l=1,lm
        if(l.eq.ls1) then
          patu(:,:) = 2.*psfmpt
          patv(:,:) = 2.*psfmpt
        endif
c
c compute pv
c
        do J=max(2,j_0s-1), j_1s+1
          DXDSIG = 0.25D0*DXV(J)*DSIG(L)
          i=1
            PV(I,J,L)=(V(I,J,L)+V(IM ,J,L))*DXDSIG*patv(i,j)
          do i=2,im
            PV(I,J,L)=(V(I,J,L)+V(I-1,J,L))*DXDSIG*patv(i,j)
          enddo
        enddo
c
c compute pu
c
        do j=max(2,j_0h),j_1s
          do i=1,im
            SPA(I,J,L)=U(I,J,L)+U(I,J+1,L)
          enddo
        enddo
        if(l.lt.ls1) then
          if(have_south_pole .or. j_0h.eq.2)
     &         spa(:,   2,l) = spa(:,   2,l)*patu(:,   2)
          if(have_north_pole) spa(:,jm-1,l) = spa(:,jm-1,l)*patu(:,jm-1)
        endif
        CALL AVRX2 (SPA(1,J_0H,L),(/MAX(2,J_0H),J_1S/))
        if(l.lt.ls1) then
          if(have_south_pole .or. j_0h.eq.2)
     &         spa(:,   2,l) = spa(:,   2,l)/patu(:,   2)
          if(have_north_pole) spa(:,jm-1,l) = spa(:,jm-1,l)/patu(:,jm-1)
        endif
        do j=max(2,j_0h),j_1s
          DYDSIG = 0.25D0*DYP(J)*DSIG(L)
          do i=1,im
            PU(I,J,L)=DYDSIG*SPA(I,J,L)*PATU(I,J)
          enddo
        enddo ! j

c compute horizontal mass convergence
        do j=max(2,j_0h),j_1s
          i=1
            CONV(I,J,L)=(PU(IM ,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
            PIT(I,J) = PIT(I,J) + CONV(I,J,L)
          do i=2,im
            CONV(I,J,L)=(PU(I-1,J,L)-PU(I,J,L)+PV(I,J,L)-PV(I,J+1,L))
            PIT(I,J) = PIT(I,J) + CONV(I,J,L)
          enddo
        enddo
        IF(have_south_pole) then
          CONV(:,1 ,L)=-sum(pv(:,2 ,l))*byim
          PIT(:,1) = PIT(:,1) + CONV(:,1,L)
        ENDIF
        If(have_north_pole) then
          CONV(:,JM,L)=+sum(pv(:,jm,l))*byim
          PIT(:,JM) = PIT(:,JM) + CONV(:,JM,L)
        ENDIF

      enddo ! l

C****
C**** END OF HORIZONTAL ADVECTION LAYER LOOP
C****
C**** COMPUTE SD, SIGMA DOT
      DO J=J_0S-1,J_1
      DO I=1,IM
        CONV(I,J,1) = PIT(I,J)
        SD(I,J,LM-1) = CONV(I,J,LM)
      ENDDO
      ENDDO
      DO L=LM-2,LS1-1,-1
      DO J=J_0S-1,J_1
      DO I=1,IM
        SD(I,J,L)=SD(I,J,L+1)+CONV(I,J,L+1)
        CONV(I,J,L+1) = SD(I,J,L)
      ENDDO
      ENDDO
      ENDDO
      DO L=LS1-2,1,-1
      DO J=J_0S-1,J_1
      DO I=1,IM
        SD(I,J,L)=SD(I,J,L+1)+CONV(I,J,L+1)-DSIG(L+1)*PIT(I,J)
        CONV(I,J,L+1) = SD(I,J,L)
      ENDDO
      ENDDO
      ENDDO

c      CALL HALO_UPDATE(GRID,SD, FROM=SOUTH)
c      CALL HALO_UPDATE(grid,PV, FROM=SOUTH)
c      CALL HALO_UPDATE(grid,PU, FROM=NORTH)

      RETURN
      END SUBROUTINE AIR_MASS_FLUX

      SUBROUTINE UPDATE_P (P,PA,DT1)
!@sum  ADVECM Calculates updated column pressures using mass fluxes
!@auth Original development team
      USE MODEL_COM, only : im,jm,lm,ptop,mrch,zatmo,u,v,t,q
      USE GEOM, only : bydxyp,imaxj
      USE DYNAMICS, only : pit
      USE DOMAIN_DECOMP_1D, only : grid, GET, am_i_root
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, GLOBALSUM
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: P(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8, INTENT(OUT) :: PA(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      REAL*8, INTENT(IN) :: DT1
      INTEGER I,J,L  !@var I,J,L  loop variables
      INTEGER IM1 ! @var IM1 = I - 1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      INTEGER :: n_exception, n_exception_all

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1, J_STRT_HALO = J_0H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE )

C**** COMPUTE PA, THE NEW SURFACE PRESSURE
      n_exception = 0
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          PA(I,J)=P(I,J)+(DT1*PIT(I,J)*BYDXYP(J))
          IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
            n_exception = n_exception + 1
          End If
        END DO
      END DO
      Call GLOBALSUM(grid, n_exception, n_exception_all, all=.true.)
      IF (n_exception_all > 0)  Then
        Do J = J_0, J_1
          im1 = im
          DO I = 1, IMAXJ(J)
            IF (PA(I,J)+PTOP.GT.1160. .or. PA(I,J)+PTOP.LT.350.) THEN
              WRITE (6,990) I,J,MRCH,P(I,J),PA(I,J),ZATMO(I,J),DT1,
     &             (U(IM1,J,L),U(I,J,L),U(IM1,J+1,L),U(I,J+1,L),
     &             V(IM1,J,L),V(I,J,L),V(IM1,J+1,L),V(I,J+1,L),
     &             T(I,J,L),L=1,LM)
            END IF
            im1 = i
          END DO
        END DO
        if(am_i_root()) write(6,*) "Pressure diagnostic error"
        call stop_model('ADVECM: Pressure diagnostic error',11)
      END IF

      IF (have_south_pole) PA(2:IM, 1)=PA(1,1)
      IF (have_north_pole) PA(2:IM,JM)=PA(1,JM)
      CALL HALO_UPDATE(grid, PA)

C****
      RETURN
  990 FORMAT (/'0PRESSURE DIAGNOSTIC     I,J,MRCH,P,PA=',3I4,2F10.2/
     &  '     ZATMO=',F10.3,' DT=',F6.1/
     &  '0    U(I-1,J)     U(I,J)   U(I-1,J+1)    U(I,J+1)    V(I-1,J)',
     &   '     V(I,J)   V(I-1,J+1)    V(I,J+1)     T(I,J)    '/
     &  (1X,9F12.3))
      END SUBROUTINE UPDATE_P

      SUBROUTINE UPDATE_UV (PA,UT,VT,PB,U,V,PS,T,SZ,DT1)
!@sum  UPDATE_UV Advects momentum (incl. coriolis) using mass fluxes
!@+    and applies the pressure gradient force.
!@auth Original development team
      use constant, only : radius,pi,twopi,omega2
      USE model_com, only : im,jm,lm,ls1,psfmpt,byim,dsig
      USE DOMAIN_DECOMP_1D, only : GRID,getDomainBounds
      USE GEOM, only : dlon,dlat,dxyv,dxyn,dxys,ravpn,ravps
     &     ,sini=>sinu,cosi=>cosu,coslatv,sinlatv,lonv,latv
      USE DYNAMICS, only : pu,pv,sd,spa
c from pgf:
      USE CONSTANT, only : grav,rgas,kapa,bykapa,bykapap1,bykapap2
      USE MODEL_COM, only : ptop,sige,sig,bydsig,zatmo
      USE DYNAMICS, only : phi
      USE GEOM, only : imaxj,dxv=>dxlatv,dyv
      USE DOMAIN_DECOMP_1D, Only : halo_update

      IMPLICIT NONE
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &     u,v,ut,vt,t,sz
      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     pa,pb,ps
c
      INTEGER I,J,L,JJ
      REAL*8 VMASS,RVMASS,DT1,DT2,DT12,DT24,FLUX,FUUP,FVUP

      real*8, dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     FUNS,FVNS,FUSWNE,FVSWNE,FUNWSE,FVNWSE,FUDN,FVDN,
     &     uxy,vxy,vmassr,vmass2,uaby4

      REAL*8, dimension(im) :: asdu,corfac,dut,dvt,fuew,fvew

      real*8 :: hemi,upol,vpol,dlat_pol,r_pol,xtmp,ytmp
     &     ,dltmp,lonup,latup,lonx,laty,uhxy,vhxy,unew,vnew
     &     ,corfacj,metfacj,utmp,vtmp
      integer :: jvpo,iup,jup,im1

c from pgf:
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)::
     &     phidn,pkdn,pgfx,spv,patu,p
      real*8, dimension(im) :: dvta
      REAL*8 PKE(LS1:LM+1)
      REAL*8 DT4
      REAL*8 PDN,PKPDN,PKPPDN,PUP,PKUP,PKPUP,PKPPUP,DP,P0,X
     &     ,BYDP,pkdnl,alph,TZBYDP,FACTORU,FACTORV,
     &     dphidt,dphidtz,dphimdt,dphimdtz,dspvdt,dspvdtz,dpk,dpkp,dpkpp

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)

      DT2=DT1/2.
      DT12=DT1/12.
      DT24=DT1/24.

      dlat_pol = pi/2.-latv(jm)
      r_pol = radius*dlat_pol

      DO J=J_0STG,J_1STG
        do i=1,im
          FUDN(i,J)  = 0.
          FVDN(i,J)  = 0.
        enddo
        DO I=1,IM-1
          VMASS2(I,J)=.5*((PB(I,J-1)+PB(I+1,J-1))*DXYN(J-1)
     &                 +(PB(I,J)+PB(I+1,J))*DXYS(J))
          VMASSR(I,J)=.5*((PA(I,J-1)+PA(I+1,J-1))*DXYN(J-1)
     &                 +(PA(I,J)+PA(I+1,J))*DXYS(J))/VMASS2(I,J)
        END DO
        I=IM
          VMASS2(I,J)=.5*((PB(I,J-1)+PB(1,J-1))*DXYN(J-1)
     &                 +(PB(I,J)+PB(1,J))*DXYS(J))
          VMASSR(I,J)=.5*((PA(I,J-1)+PA(1,J-1))*DXYN(J-1)
     &                 +(PA(I,J)+PA(1,J))*DXYS(J))/VMASS2(I,J)
      END DO

c from pgf:
      DT4=DT1/4.
      DO L=LS1,LM+1
        PKE(L)=(PSFMPT*SIGE(L)+PTOP)**KAPA
      END DO
      IF (have_south_pole) PGFX(:,1)=0.
      IF (have_north_pole) PGFX(:,JM)=0.
      p(:,:) = ps(:,:)
      DO J=max(J_0H,1),J_1 ! compute halo lats
      DO I=1,IMAXJ(J)
        pkdn(i,j)=(p(i,j)+ptop)**kapa
        phidn(i,j)=zatmo(i,j)
      ENDDO
      ENDDO

C****
C**** BEGINNING OF MAIN LAYER LOOP
C****
      DO L=1,LM

      if(l.eq.ls1) then
        do j=j_0stg,j_1stg
        do i=1,im
          vmassr(i,j)=1d0
          vmass2(i,j)=PSFMPT*DXYV(J)
        enddo
        enddo
      endif

C****
C**** HORIZONTAL ADVECTION OF MOMENTUM
C****

c Simple semi-lagrangian advection at the poles.
c Must be performed before ut,vt are updated at lower latitudes.

c directional conventions:
c in the NH, +y points from IDL to GM, +x from 90 E to 90 W
c in the SH, +y points from GM to IDL, +x from 90 W to 90 E

      do j=j_0stg,j_1stg
        if(j.gt.2 .and. j.lt.jm) cycle
        corfacj = dt1*omega2*sinlatv(j)
        if(latv(j).lt.0.) then
          hemi = -1.
          do jj=2,3     ! convert upwind velocities to xy components
            do i=1,im
              uxy(i,jj) = cosi(i)*ut(i,jj,l)-hemi*sini(i)*vt(i,jj,l)
              vxy(i,jj) = cosi(i)*vt(i,jj,l)+hemi*sini(i)*ut(i,jj,l)
            enddo
          enddo
          upol =  sum(uxy(:,2))*byim
          vpol =  sum(vxy(:,2))*byim
          jvpo = 2
        else
          hemi = +1.
          do jj=jm-1,jm ! convert upwind velocities to xy components
            do i=1,im
              uxy(i,jj) = cosi(i)*ut(i,jj,l)-hemi*sini(i)*vt(i,jj,l)
              vxy(i,jj) = cosi(i)*vt(i,jj,l)+hemi*sini(i)*ut(i,jj,l)
            enddo
          enddo
          upol =  sum(uxy(:,jm))*byim
          vpol =  sum(vxy(:,jm))*byim
          jvpo = jm
        endif
        do i=1,im
          uhxy = cosi(i)*u(i,j,l)-hemi*sini(i)*v(i,j,l)
          vhxy = cosi(i)*v(i,j,l)+hemi*sini(i)*u(i,j,l)
          xtmp =  r_pol*sini(i)     -uhxy*dt1
          ytmp = -r_pol*cosi(i)*hemi-vhxy*dt1
          dltmp = sqrt(xtmp**2+ytmp**2)/radius
          latup = hemi*(pi/2.-dltmp)
          lonup = atan2(hemi*ytmp,xtmp)+pi/2.
          if(lonup.lt.0.) lonup=lonup+twopi
          iup = int(1d0 + lonup/dlon)
          if(iup.gt.1) then
            lonx = (lonup - lonv(iup-1))/dlon
            im1 = iup-1
          else
            lonx = lonup/dlon
            im1 = im
          endif
          if(dltmp.gt.dlat_pol) then ! regular lat-lon interpolation
            jup = int(2d0 + (latup-latv(2))/dlat)
            if(jup.gt.3 .and. jup.lt.jm-1) then
              if(hemi.eq.-1.) then
                jup = 3
              else
                jup = jm-1
              endif
              write(6,*) 'sl advecv: jup out of bounds'
c              call stop_model('sl advecv: jup out of bounds',255)
            endif
            laty = (latup - latv(jup))/dlat
            unew =(1.-lonx)*((1.-laty)*uxy(im1,jup)+laty*uxy(im1,jup+1))
     &           +   (lonx)*((1.-laty)*uxy(iup,jup)+laty*uxy(iup,jup+1))
            vnew =(1.-lonx)*((1.-laty)*vxy(im1,jup)+laty*vxy(im1,jup+1))
     &           +   (lonx)*((1.-laty)*vxy(iup,jup)+laty*vxy(iup,jup+1))
          else                  ! interpolation within the polar cell
            laty = dltmp/dlat_pol
            unew =(1.-lonx)*((1.-laty)*upol+laty*uxy(im1,jvpo))
     &           +   (lonx)*((1.-laty)*upol+laty*uxy(iup,jvpo))
            vnew =(1.-lonx)*((1.-laty)*vpol+laty*vxy(im1,jvpo))
     &           +   (lonx)*((1.-laty)*vpol+laty*vxy(iup,jvpo))
          endif
          ut(i,j,l) = cosi(i)*unew+hemi*sini(i)*vnew
     &         +corfacj*v(i,j,l)
          vt(i,j,l) = cosi(i)*vnew-hemi*sini(i)*unew
     &         -corfacj*u(i,j,l)
        enddo
      enddo ! j

! compute a-grid u.  spa contains a factor of 2
      do j=max(2,j_0s-1),j_1s
        i=1
          uaby4(i,j) = .0625*(spa(im,j,l)+spa(i,j,l))
        do i=2,im
          uaby4(i,j) = .0625*(spa(i-1,j,l)+spa(i,j,l))
        enddo
      enddo

C**** SOUTH-NORTH FLUXES AT PRIMARY LATITUDES
      DO J=max(2,J_0STG-1),min(J_1STG,jm-1) ! compute halo fluxes
        DO I=1,IM-1
          FLUX=DT12*(PV(I,J,L)+PV(I+1,J,L)+PV(I,J+1,L)+PV(I+1,J+1,L))
          FUNS(I,J)=FLUX*(U(I,J,L)+U(I,J+1,L))
          FVNS(I,J)=FLUX*(V(I,J,L)+V(I,J+1,L))
        ENDDO
        I=IM
          FLUX=DT12*(PV(I,J,L)+PV(1  ,J,L)+PV(I,J+1,L)+PV(1  ,J+1,L))
          FUNS(I,J)=FLUX*(U(I,J,L)+U(I,J+1,L))
          FVNS(I,J)=FLUX*(V(I,J,L)+V(I,J+1,L))
C**** SOUTHWEST-NORTHEAST AND SOUTHEAST-NORTHWEST FLUXES
        I=1
          FLUX=DT24*(+PU(IM ,J,L)+PU(I,J,L)+PV(I,J,L)+PV(I,J+1,L))
          FUSWNE(I,J)=FLUX*(U(I,J+1,L)+U(IM ,J,L))
          FVSWNE(I,J)=FLUX*(V(I,J+1,L)+V(IM ,J,L))
          FLUX=DT24*(-PU(IM ,J,L)-PU(I,J,L)+PV(I,J,L)+PV(I,J+1,L))
          FUNWSE(I,J)=FLUX*(U(IM ,J+1,L)+U(I,J,L))
          FVNWSE(I,J)=FLUX*(V(IM ,J+1,L)+V(I,J,L))
        DO I=2,IM
          FLUX=DT24*(+PU(I-1,J,L)+PU(I,J,L)+PV(I,J,L)+PV(I,J+1,L))
          FUSWNE(I,J)=FLUX*(U(I,J+1,L)+U(I-1,J,L))
          FVSWNE(I,J)=FLUX*(V(I,J+1,L)+V(I-1,J,L))
          FLUX=DT24*(-PU(I-1,J,L)-PU(I,J,L)+PV(I,J,L)+PV(I,J+1,L))
          FUNWSE(I,J)=FLUX*(U(I-1,J+1,L)+U(I,J,L))
          FVNWSE(I,J)=FLUX*(V(I-1,J+1,L)+V(I,J,L))
        ENDDO
      ENDDO

C****
C**** VERTICAL AND HORIZONTAL CONVERGENCES OF MOMENTUM, CORIOLIS/METRIC TERM
C****
      if(l.lt.lm) then
        DO J=max(3,J_0STG),min(jm-1,J_1STG)
          corfacj=dt1*omega2*sinlatv(j)
          metfacj=dt1*sinlatv(j)/(radius*coslatv(j))
          I=1
            FLUX=DT12*(PU(IM ,J,L)+PU(IM ,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
            FUEW(I)=FLUX*(U(IM ,J,L)+U(I,J,L))
            FVEW(I)=FLUX*(V(IM ,J,L)+V(I,J,L))
          DO I=2,IM
            FLUX=DT12*(PU(I-1,J,L)+PU(I-1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
            FUEW(I)=FLUX*(U(I-1,J,L)+U(I,J,L))
            FVEW(I)=FLUX*(V(I-1,J,L)+V(I,J,L))
          ENDDO
          DO I=1,IM-1
            ASDU(I)=DT2*(
     &           (SD(I,J-1,L)+SD(I+1,J-1,L))*RAVPN(J-1)+
     &           (SD(I,J  ,L)+SD(I+1,J  ,L))*RAVPS(J  ))
            DUT(I) =
     &           +FUEW(I)-FUEW(I+1) +FUNS(I,J-1)-FUNS(I,J)
     &         +FUSWNE(I,J-1)-FUSWNE(I+1,J) -FUNWSE(I,J)+FUNWSE(I+1,J-1)
            DVT(I) =
     &           +FVEW(I)-FVEW(I+1) +FVNS(I,J-1)-FVNS(I,J)
     &         +FVSWNE(I,J-1)-FVSWNE(I+1,J) -FVNWSE(I,J)+FVNWSE(I+1,J-1)
            corfac(i)=corfacj+metfacj*
     &           (uaby4(i+1,j)+uaby4(i+1,j-1)+uaby4(i,j)+uaby4(i,j-1))
          END DO
          I=IM
          ASDU(I)=DT2*(
     &         (SD(I,J-1,L)+SD(1,J-1,L))*RAVPN(J-1)+
     &         (SD(I,J  ,L)+SD(1,J  ,L))*RAVPS(J  ))
          corfac(i)=corfacj+metfacj*
     &         (uaby4(1,j)+uaby4(1,j-1)+uaby4(i,j)+uaby4(i,j-1))
          DUT(I) =
     &         +FUEW(I)-FUEW(1  ) +FUNS(I,J-1)-FUNS(I,J)
     &         +FUSWNE(I,J-1)-FUSWNE(1  ,J) -FUNWSE(I,J)+FUNWSE(1  ,J-1)
          DVT(I) =
     &         +FVEW(I)-FVEW(1  ) +FVNS(I,J-1)-FVNS(I,J)
     &         +FVSWNE(I,J-1)-FVSWNE(1  ,J) -FVNWSE(I,J)+FVNWSE(1  ,J-1)
          do i=1,im
            FUUP  =ASDU(I)  *(U(I,J,L)+U(I,J,L+1))
            FVUP  =ASDU(I)  *(V(I,J,L)+V(I,J,L+1))
            DUT(I)  =DUT(I)  -FUDN(I,J) + FUUP
            DVT(I)  =DVT(I)  -FVDN(I,J) + FVUP
            FUDN(I,J) = FUUP
            FVDN(I,J) = FVUP
            RVMASS=1./(VMASS2(I,J)*DSIG(L))
            UT(I,J,L)=UT(I,J,L)*VMASSR(I,J)+DUT(I)*RVMASS
     &           +corfac(i)*v(i,j,l)
            VT(I,J,L)=VT(I,J,L)*VMASSR(I,J)+DVT(I)*RVMASS
     &           -corfac(i)*u(i,j,l)
          enddo
        ENDDO
      else
        DO J=max(3,J_0STG),min(jm-1,J_1STG)
          corfacj=dt1*omega2*sinlatv(j)
          metfacj=dt1*sinlatv(j)/(radius*coslatv(j))
          I=1
            FLUX=DT12*(PU(IM ,J,L)+PU(IM ,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
            FUEW(I)=FLUX*(U(IM ,J,L)+U(I,J,L))
            FVEW(I)=FLUX*(V(IM ,J,L)+V(I,J,L))
          DO I=2,IM
            FLUX=DT12*(PU(I-1,J,L)+PU(I-1,J-1,L)+PU(I,J,L)+PU(I,J-1,L))
            FUEW(I)=FLUX*(U(I-1,J,L)+U(I,J,L))
            FVEW(I)=FLUX*(V(I-1,J,L)+V(I,J,L))
          ENDDO
          DO I=1,IM-1
            DUT(I) =
     &           +FUEW(I)-FUEW(I+1) +FUNS(I,J-1)-FUNS(I,J)
     &         +FUSWNE(I,J-1)-FUSWNE(I+1,J) -FUNWSE(I,J)+FUNWSE(I+1,J-1)
            DVT(I) =
     &           +FVEW(I)-FVEW(I+1) +FVNS(I,J-1)-FVNS(I,J)
     &         +FVSWNE(I,J-1)-FVSWNE(I+1,J) -FVNWSE(I,J)+FVNWSE(I+1,J-1)
            corfac(i)=corfacj+metfacj*
     &           (uaby4(i+1,j)+uaby4(i+1,j-1)+uaby4(i,j)+uaby4(i,j-1))
          ENDDO
          I=IM
          corfac(i)=corfacj+metfacj*
     &         (uaby4(1,j)+uaby4(1,j-1)+uaby4(i,j)+uaby4(i,j-1))
          DUT(I) =
     &         +FUEW(I)-FUEW(1  ) +FUNS(I,J-1)-FUNS(I,J)
     &         +FUSWNE(I,J-1)-FUSWNE(1  ,J) -FUNWSE(I,J)+FUNWSE(1  ,J-1)
          DVT(I) =
     &         +FVEW(I)-FVEW(1  ) +FVNS(I,J-1)-FVNS(I,J)
     &         +FVSWNE(I,J-1)-FVSWNE(1  ,J) -FVNWSE(I,J)+FVNWSE(1  ,J-1)
          do i=1,im
            DUT(I)  =DUT(I)  -FUDN(I,J)
            DVT(I)  =DVT(I)  -FVDN(I,J)
            RVMASS=1./(VMASS2(I,J)*DSIG(L))
            UT(I,J,L)=UT(I,J,L)*VMASSR(I,J)+DUT(I)*RVMASS
     &           +corfac(i)*v(i,j,l)
            VT(I,J,L)=VT(I,J,L)*VMASSR(I,J)+DVT(I)*RVMASS
     &           -corfac(i)*u(i,j,l)
          enddo
        ENDDO
      endif

C****
C**** PRESSURE GRADIENT FORCE
C****
        IF(L.GE.LS1) THEN  ! constant-pressure levels: only need phi
          pdn=sige(l)*psfmpt+ptop
          pup=sige(l+1)*psfmpt+ptop
          pkdnl=pke(l)
          pkup=pke(l+1)
          pkpdn=pkdnl*pdn
          dp=dsig(l)*psfmpt
          bydp=1./dp
          p0=sig(l)*psfmpt+ptop
          pkpup=pkup*pup
          dpk = (pkdnl-pkup)*bykapa
          dpkp = (pkpdn-pkpup)*bykapap1
          dpkpp = (pkpdn*pdn-pkpup*pup)*bykapap2
          dphidt = dpk
          dphimdt = bykapa*(pkdnl-dpkp*bydp)
          dphidtz = (dphidt*p0 -dpkp)*2.*bydp
          dphimdtz = (dphimdt*p0 +bykapap1*(bydp*dpkpp-pkpdn))*2.*bydp
          dphidt = dphidt*rgas
          dphimdt = dphimdt*rgas
          dphidtz = dphidtz*rgas
          dphimdtz = dphimdtz*rgas
          do j=max(j_0h,1),j_1 ! compute halo lats
            do i=1,imaxj(j)
              phi(i,j,l)=phidn(i,j)+dphimdt*t(i,j,l)+dphimdtz*sz(i,j,l)
              phidn(i,j)=phidn(i,j)+dphidt*t(i,j,l)+dphidtz*sz(i,j,l)
            enddo
          enddo
          IF (have_south_pole) PHI(2:IM, 1,L)=PHI(1, 1,L)
          IF (have_north_pole) PHI(2:IM,JM,L)=PHI(1,JM,L)
          DO J=Max(2,J_0STG-1),min(jm-1,J_1STG)
            do i=1,im-1
              PGFX(I,J) = PHI(I+1,J,L)-PHI(I,J,L)
            enddo
            I=IM
            PGFX(I,J) = PHI(1,J,L)-PHI(I,J,L)
          END DO
          CALL AVRX2 (PGFX,(/MAX(2,J_0H),MIN(JM-1,J_1H)/))
          DO J=J_0STG,J_1STG
            factoru = -dt1*.5/dxv(j)
            FACTORv = -DT1*.5/DYV(J)
            DO I=1,IM
              UT(I,J,L)=UT(I,J,L)+factoru*(PGFX(I,J)+PGFX(I,J-1))
              dvta(i) = factorv*(PHI(I,J,L)-PHI(I,J-1,L))
            END DO
            do i=1,im-1
              vt(i,j,l) = vt(i,j,l) + (dvta(i)+dvta(i+1))
            enddo
            i=im
            vt(i,j,l) = vt(i,j,l) + (dvta(i)+dvta(1  ))
          END DO
        ELSE                    ! sigma layers
          DO J=max(J_0H,1),J_1 ! compute halo lats
          DO I=1,IMAXJ(J)
            pdn=sige(l)*p(i,j)+ptop
            pup=sige(l+1)*p(i,j)+ptop
            pkpdn=pkdn(i,j)*pdn
            dp=dsig(l)*p(i,j)
            bydp=1./dp
            p0=sig(l)*p(i,j)+ptop
            pkup=pup**kapa
            pkpup=pkup*pup
            dpk = (pkdn(i,j)-pkup)*bykapa
            dpkp = (pkpdn-pkpup)*bykapap1
            dpkpp = (pkpdn*pdn-pkpup*pup)*bykapap2
            dphidt = dpk
            dphimdt = bykapa*(pkdn(i,j)-dpkp*bydp)
            dspvdt = dpkp-ptop*dpk
            dphidtz = dphidt*p0 -dpkp
            dphimdtz = dphimdt*p0 +bykapap1*(bydp*dpkpp-pkpdn)
            dspvdtz = dspvdt*p0+ptop*dpkp-dpkpp
            tzbydp = 2.*sz(i,j,l)*bydp
            spv(i,j) = rgas*bydp*(dspvdt*t(i,j,l)+dspvdtz*tzbydp)
            phi(i,j,l) = phidn(i,j)
     &           +rgas*(dphimdt*t(i,j,l)+dphimdtz*tzbydp)
            phidn(i,j) = phidn(i,j)
     &           +rgas*(dphidt*t(i,j,l)+dphidtz*tzbydp)
            pkdn(i,j) = pkup
          ENDDO
          ENDDO
          IF (have_south_pole) THEN
            SPV(2:IM,1)=SPV(1,1)
            PHI(2:IM,1,L)=PHI(1,1,L)
          ENDIF
          IF (have_north_pole) THEN
            SPV(2:IM,JM)=SPV(1,JM)
            PHI(2:IM,JM,L)=PHI(1,JM,L)
          ENDIF
          DO J=Max(2,J_0STG-1),min(jm-1,J_1STG)
            do i=1,im-1
              patu(i,j) = P(I+1,J)+P(I,J)
              alph = 4.*SPV(I+1,J)*SPV(I,J)/(SPV(I+1,J)+SPV(I,J))
              PGFX(I,J) = patu(i,j)*(PHI(I+1,J,L)-PHI(I,J,L)) +
     &             alph*(P(I+1,J)-P(I,J))
            enddo
            I=IM
            patu(i,j) = P(1,J)+P(I,J)
            alph = 4.*SPV(1,J)*SPV(I,J)/(SPV(1,J)+SPV(I,J))
            PGFX(I,J) = patu(i,j)*(PHI(1,J,L)-PHI(I,J,L)) +
     &           alph*(P(1,J)-P(I,J))
          END DO
          CALL AVRX2 (PGFX,(/MAX(2,J_0H),MIN(JM-1,J_1H)/))
          DO J=Max(2,J_0STG-1),min(jm-1,J_1STG)
            pgfx(:,j) = pgfx(:,j)/patu(:,j)
          END DO
          DO J=J_0STG,J_1STG
            factoru = -dt1*.5/dxv(j)
            FACTORv = -DT1*.5/DYV(J)
            DO I=1,IM
              UT(I,J,L)=UT(I,J,L)+factoru*(PGFX(I,J)+PGFX(I,J-1))
              alph = 4.*SPV(I,J-1)*SPV(I,J)/(SPV(I,J-1)+SPV(I,J))
              dvta(i)=factorv*(PHI(I,J,L)-PHI(I,J-1,L)+
     &             alph*(P(I,J)-P(I,J-1))/(P(I,J)+P(I,J-1)) )
            END DO
            do i=1,im-1
              vt(i,j,l) = vt(i,j,l) + (dvta(i)+dvta(i+1))
            enddo
            i=im
            vt(i,j,l) = vt(i,j,l) + (dvta(i)+dvta(1  ))
          END DO
        ENDIF

        call isotropuv2(ut(1,j_0h,l),vt(1,j_0h,l))

      ENDDO ! end loop over l

c Fill halos of updated u,v
      call halo_update(grid,ut)
      call halo_update(grid,vt)

      RETURN
      END SUBROUTINE UPDATE_UV

      SUBROUTINE AVRX2(X,jrange)
!@sum  AVRX Smoothes zonal mass flux and geopotential near the poles
!@auth Original development team
      USE MODEL_COM, only : im,jm,imh
      USE GEOM, only : dlon,dxp,dyp,bydyp
      !USE DYNAMICS, only : xAVRX
C**** THIS VERSION OF AVRX DOES SO BY TRUNCATING THE FOURIER SERIES.
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      IMPLICIT NONE
      REAL*8, INTENT(INOUT) ::
     &     X(IM,grid%J_STRT_HALO:grid%J_STOP_HALO)
      Integer, Intent(In) :: jrange(2)
      REAL*8, ALLOCATABLE, SAVE  :: DRAT(:)
      REAL*8, SAVE ::  BYSN(IMH)
      REAL*8, DIMENSION(0:IMH) :: AN,BN
CCC   INTEGER, SAVE :: NMIN(grid%J_STRT_HALO:grid%J_STOP_HALO)
      INTEGER, ALLOCATABLE, SAVE :: NMIN(:)
      INTEGER J,N
      LOGICAL, SAVE :: init = .false.
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, J0, J1

c      if ( present(X) ) goto 1000

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_HALO = J_0H, J_STOP_HALO = J_1H,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
C

      IF (.NOT. init) THEN
        init = .true.
C       CALL FFT0(IM)
        j0 = MAX(1,J_0H)
        j1 = MIN(JM,J_1H)
        ALLOCATE(DRAT(j0:j1), NMIN(j0:j1))
        DO N=1,IMH
          BYSN(N)=1d0/SIN(.5*DLON*N)
        END DO
        DO J=j0,j1
          DRAT(J) = DXP(J)*BYDYP(3)
          DO N=IMH,1,-1
            IF(BYSN(N)*DRAT(J) .GT.1.) THEN
              NMIN(J) = N+1
              EXIT
            ENDIF
          END DO
        END DO
      END IF

c      RETURN
C****
!!!      ENTRY AVRX (X)
c 1000 continue
C****

c      If (Present(jrange)) Then
        j0 = jrange(1)
        j1 = jrange(2)
c      Else
c        call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)
c        j0=J_0S
c        j1=J_1S
c      End If

      DO J=j0,j1
        IF (DRAT(J).GT.1) CYCLE
        CALL FFT (X(1,J),AN,BN)
        DO N=NMIN(J),IMH-1
          AN(N)=BYSN(N)*DRAT(J) * AN(N)
          BN(N)=BYSN(N)*DRAT(J) * BN(N)
        END DO
        AN(IMH) = BYSN(IMH)*DRAT(J) * AN(IMH)
        CALL FFTI(AN,BN,X(1,J))
      END DO

      RETURN
      END SUBROUTINE AVRX2

      SUBROUTINE FILTER2
!@sum  FILTER Performs 8-th order shapiro filter in zonal direction
!@auth Original development team
!@calls SHAP1D
      USE CONSTANT, only : kapa,rgas
      USE MODEL_COM, only : im,jm,lm,ls1,t,p,q,wm,zatmo,ptop,byim,sig
      USE SOMTQ_COM, only : tmom,qmom
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      USE ATMDYN, only : shap1d
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     X,Y,POLD,PRAT
      INTEGER I,J,L,NSHAP
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      REAL*8 initialTotalEnergy, finalTotalEnergy
      real*8 getTotalEnergy ! external for now
      real*8, dimension(im) :: rhosrf,pgfx

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

C**** Initialise total energy (J/m^2)
      initialTotalEnergy = getTotalEnergy()

      do j=j_0s,j_1s
        pold(:,j)=p(:,j)        ! save old pressure
        do i=1,im
          rhosrf(i) = ((p(i,j)+ptop)**(1.-kapa))/(rgas*t(i,j,1))
        enddo
        do i=1,im-1
          pgfx(i) = (p(i+1,j)-p(i,j))+
     &         .5*(rhosrf(i+1)+rhosrf(i))*(zatmo(i+1,j)-zatmo(i,j))
        enddo
        i=im
          pgfx(i) = (p(1,j)-p(i,j))+
     &         .5*(rhosrf(1)+rhosrf(i))*(zatmo(1,j)-zatmo(i,j))
        pgfx = pgfx - sum(pgfx)*byim
        x(1,j) = 0.
        do i=2,im
          x(i,j) = x(i-1,j) + pgfx(i-1)
        enddo
        y(:,j) = x(:,j)
      enddo
      if(im.gt.144) then
        nshap = 6
      else
        nshap = 8
      endif
      call shap1d (nshap,x)
      call isotropslp2(x)
      do j=j_0s,j_1s
        p(:,j) = p(:,j) + 1d0*(x(:,j)-y(:,j))
      enddo

C**** Scale mixing ratios (incl moments) to conserve mass/heat
      DO J=J_0S,J_1S
        DO I=1,IM
          PRAT(I,J)=POLD(I,J)/P(I,J)
        END DO
      END DO
      DO L=1,LS1-1
      DO J=J_0S,J_1S
      DO I=1,IM
c adjust pot. temp. to maintain unchanged absolute temp.
         T(I,J,L)= T(I,J,L)*
     &       ((POLD(I,J)*SIG(L)+PTOP)/(P(I,J)*SIG(L)+PTOP))**KAPA
         Q(I,J,L)= Q(I,J,L)*PRAT(I,J)
        WM(I,J,L)=WM(I,J,L)*PRAT(I,J)
        QMOM(:,I,J,L)=QMOM(:,I,J,L)*PRAT(I,J)
      END DO
      END DO
      END DO

      CALL CALC_AMPK(LS1-1)

C**** This fix adjusts thermal energy to conserve total energy TE=KE+PE
      finalTotalEnergy = getTotalEnergy()
      call addEnergyAsDiffuseHeat(finalTotalEnergy - initialTotalEnergy)

      RETURN
      END SUBROUTINE FILTER2

      subroutine fltry3(q3d,strength)
!@sum  fltry3 noise reduction filter for a velocity-type field
!@sum  at secondary latitudes
      use model_com, only : im,jm,lm
      use domain_decomp_1d, only : get,grid,halo_update,north,south
      implicit none
      real*8 :: dt
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: q3d
      real*8 :: strength
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: yn
      real*8 by4ton,yvby4ton
      integer i,j,l,jj,n,nshap,j_0stg,j_1stg,j_1f
      real*8, dimension(im) :: yjm1,yj
      logical :: have_south_pole,have_north_pole

      call getDomainBounds(grid, j_strt_stgr = j_0stg, j_stop_stgr = j_1stg,
     &         have_south_pole = have_south_pole,
     &         have_north_pole = have_north_pole)

      if(im.gt.144) then
        nshap = 6
      else
        nshap = 8
      endif
      by4ton=1./(4.**nshap)
      yvby4ton = min(strength,1d0)*by4ton*((-1)**(nshap))

      if(have_north_pole) then
        j_1f=jm-1
      else
        j_1f=j_1stg
      endif

      do l=1,lm
        do j=j_0stg,j_1stg
          yn(:,j,l)=q3d(:,j,l)
        enddo
      enddo
      do n=1,nshap
        call halo_update(grid, yn)
        do l=1,lm
          if(have_south_pole) then ! pole-crossing conditions
            yjm1(1:im/2)    = -yn(im/2+1:im,2,l)
            yjm1(im/2+1:im) = -yn(1:im/2,2,l)
          else
            yjm1(:)   = yn(:,j_0stg-1,l)
          endif
          do j=j_0stg,j_1f
            do i=1,im
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)+yn(i,j+1,l)
              yjm1(i) = yj(i)
            enddo
          enddo
          if(have_north_pole) then ! pole-crossing conditions
            j=jm
            do i=1,im/2
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)-yn(i+im/2,j,l)
            enddo
            do i=im/2+1,im
              yj(i)   = yn(i,j,l)
              yn(i,j,l) = yjm1(i)-yj(i)-yj(i)-yj(i-im/2)
            enddo
          endif
        enddo                 ! l
      enddo                   ! nshap
      do l=1,lm
        do j=j_0stg,j_1stg
          q3d(:,j,l) = q3d(:,j,l) -yn(:,j,l)*yvby4ton
        enddo
      enddo
      return
      end subroutine fltry3

      SUBROUTINE FLTRUV2(U,V,strength)
!@sum  FLTRUV Filters 2 gridpoint noise from the velocity fields
!@auth Original development team
      USE CONSTANT, only : sha
      USE MODEL_COM, only : im,jm,lm,byim,mrch,dt,t,ang_uv
     *  ,DT_XUfilter,DT_XVfilter,DT_YVfilter,DT_YUfilter
     &  ,do_polefix
      USE DYNAMICS, only : pdsig
      USE GEOM, only : dxyn,dxys
c      USE DIAG, only : diagcd
C**********************************************************************
C**** FILTERING IS DONE IN X-DIRECTION WITH A 8TH ORDER SHAPIRO
C**** FILTER. THE EFFECT OF THE FILTER IS THAT OF DISSIPATION AT
C**** THE SMALLEST SCALES.
C**********************************************************************
      USE DOMAIN_DECOMP_1D, only : grid, GET
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     *     INTENT(INOUT) :: U,V
      REAL*8, INTENT(IN) :: STRENGTH
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     *     DUT,DVT,USAVE,VSAVE
      REAL*8 X(IM),YV(max(2*JM,IM)),DP(IM)
      REAL*8 XUby4toN,XVby4toN,YVby4toN,YUby4toN
      REAL*8 :: DT1=0.
      INTEGER I,J,K,L,N,IP1  !@var I,J,L,N  loop variables
      REAL*8 YV2,YVJ,YVJM1,X1,XI,XIM1
      INTEGER :: NSHAP ! NSHAP MUST BE EVEN
      REAL*8 angm,dpt,D2V,D2U,by4ton
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0STG, J_1STG, J_0S, J_1S, J_0H, J_1H
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &               J_STRT_HALO = J_0H, J_STOP_HALO = J_1H,
     &               J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)
C****

      if(im.gt.144) then
        nshap = 6
      else
        nshap = 8
      endif
      by4ton=1./(4.**nshap)

      USAVE=U ; VSAVE=V
      XUby4toN = by4toN*strength
      XVby4toN = by4toN*strength
C****
C**** Filtering in east-west direction
C****
      DO 350 L=1,LM
C**** Filter U component of velocity
      DO 240 J=J_0STG,J_1STG
      DO 210 I=1,IM
  210 X(I) = U(I,J,L)
      DO 230 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 220 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  220 XIM1 = XI
  230 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 240 I=1,IM
  240 U(I,J,L) = U(I,J,L) - X(I)*XUby4toN
C**** Filter V component of velocity
      DO 340 J=J_0STG,J_1STG
      DO 310 I=1,IM
  310 X(I) = V(I,J,L)
      DO 330 N=1,NSHAP
      X1   = X(1)
      XIM1 = X(IM)
      DO 320 I=1,IM-1
      XI   = X(I)
      X(I) = XIM1-XI-XI+X(I+1)
  320 XIM1 = XI
  330 X(IM)= XIM1-X(IM)-X(IM)+X1
      DO 340 I=1,IM
  340 V(I,J,L) = V(I,J,L) - X(I)*XVby4toN
  350 CONTINUE

C**** Conserve angular momentum along latitudes
c***  The following halo is not needed because PDSIG halo is up to date
c***      CALL HALO_UPDATE_COLUMN(grid, PDSIG, FROM=SOUTH)
      DO L=1,LM
        DO J=J_0STG,J_1STG
          ANGM=0.
          DPT=0.
          I=IM
          DO IP1=1,IM
            DP(I)=0.5*((PDSIG(L,IP1,J-1)+PDSIG(L,I,J-1))*DXYN(J-1)
     *           +(PDSIG(L,IP1,J  )+PDSIG(L,I,J  ))*DXYS(J  ))
            ANGM=ANGM-DP(I)*(U(I,J,L)-USAVE(I,J,L))
            DPT=DPT+DP(I)
            I=IP1
          END DO
          DO I=1,IM
            if (ang_uv.eq.1) U(I,J,L)=U(I,J,L)+ANGM/DPT
            DUT(I,J,L)=(U(I,J,L)-USAVE(I,J,L))*DP(I)
            DVT(I,J,L)=(V(I,J,L)-VSAVE(I,J,L))*DP(I)
          END DO
        END DO
      END DO

C**** Call diagnostics only for even time step
c      IF (MRCH.eq.2) THEN
c        CALL DIAGCD(grid,5,UT,VT,DUT,DVT,DT1)
c      END IF

      RETURN
      END SUBROUTINE FLTRUV2

      subroutine isotropslp2(slp)
      use MODEL_COM, only : im,jm,dtsrc
      USE DOMAIN_DECOMP_1D, Only : GET,grid
      use GEOM, only : cosp,dxp
c      USE DYNAMICS, only : COS_LIMIT
      USE ATMDYN, only : shap1
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: slp
      real*8 :: fac,cos_limit,k
      integer :: ipole,j,jcut,jinc,jp
      integer :: j_0s, j_1s, hemi

      call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      if(im.gt.144) then
        cos_limit=.23d0 ! necessary?
        k=5d3
      else
        cos_limit=.15d0
        k=1d3
      endif

      Do j = j_0s, j_1s
        If(cosp(j) .gt. cos_limit) Cycle
        fac = k*dtsrc/(dxp(j)*dxp(j))
        Call shap1(slp(1,j),im,fac)
      enddo

      return
      end subroutine isotropslp2

      subroutine isotropuv2(u,v)
!@sum  isotropuv isotropizes the velocity field in the near-polar row(s)
!@auth M. Kelley
      USE MODEL_COM, only : im,imh,jm,dt
      USE DOMAIN_DECOMP_1D, Only : GET, grid
      USE GEOM, only : sinlatv,coslatv,dxv=>dxlatv,cosi=>cosu,sini=>sinu
c      USE DYNAMICS, only : COS_LIMIT
      USE ATMDYN, only : shap1
      implicit none
      real*8, parameter :: klo=1d3,khi=1d7
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     *  U, V
      real*8 :: fac,k,cos_limit
      real*8, dimension(im) :: ua,va
      real*8, dimension(0:imh) :: an,bn
      integer :: i,j,hemi
      Integer :: J_0STG, J_1STG

      if(im.gt.144) then
        cos_limit=.23d0 ! necessary?
      else
        cos_limit=.15d0
      endif

      call getDomainBounds(grid, J_STRT_STGR = J_0STG, J_STOP_STGR = J_1STG)

      do j = J_0STG, J_1STG
        If(coslatv(j) .gt. cos_limit) Cycle
        hemi = sign(1.d0,sinlatv(j))

c compute xy velocities
        do i=1,im
          ua(i) = cosi(i)*u(i,j)-hemi*sini(i)*v(i,j)
          va(i) = cosi(i)*v(i,j)+hemi*sini(i)*u(i,j)
        enddo
c filter the xy velocities
        if(j.eq.2 .or. j.eq.jm) then
! really strong filtering right at the pole
          call fft(ua,an,bn)
          an(2:imh) = 0.
          bn(2:imh) = 0.
          call ffti(an,bn,ua)
          call fft(va,an,bn)
          an(2:imh) = 0.
          bn(2:imh) = 0.
          call ffti(an,bn,va)
        else
          k = maxval(abs(u(:,j)))*2.*dt/dxv(j)
          if(k.lt.0.5) then
            k = klo
          else if(k.gt.1d0) then
            k = khi
          else
            k = klo + 2d0*(k-0.5)*(khi-klo)
          endif
          fac = k*dt/(dxv(j)*dxv(j))
          call shap1(ua,im,fac)
          call shap1(va,im,fac)
        endif
c convert xy velocities back to polar coordinates
        do i=1,im
          u(i,j) = cosi(i)*ua(i)+hemi*sini(i)*va(i)
          v(i,j) = cosi(i)*va(i)-hemi*sini(i)*ua(i)
        enddo
      enddo                     ! j

      return
      end subroutine isotropuv2

      SUBROUTINE AADVT3 (MMA,RM,RMOM,RZ,SD,PU,PV,DT)
!@sum  AADVT advection driver
!@vers 2013/03/27
!@auth G. Russell, modified by Maxwell Kelley
c****
c**** AADVT advects tracers using the Quadradic Upstream Scheme.
c****
c**** input:
c****  pu,pv,sd (kg/s) = east-west,north-south,vertical mass fluxes
c****      qlimit = whether moment limitations should be used
C****         DT (s) = time step
c****
c**** input/output:
c****     rm = tracer concentration
c****   rmom = moments of tracer concentration
c****     mma (kg) = fluid mass
c****
      USE DOMAIN_DECOMP_1D, only: grid, get
      USE DOMAIN_DECOMP_1D, only:
     &     HALO_UPDATE,HALO_UPDATE_COLUMN,
     &     NORTH,SOUTH
      USE QUSDEF
      USE QUSCOM, ONLY : IM,JM,LM
      USE DYNAMICS, only : MUs,MVs,MWs
      USE GEOM, only : imaxj
      IMPLICIT NONE

      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) ::
     &                  rm,mma,rz
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM)
     &               :: rmom

      REAL*8, INTENT(IN) :: DT
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm),
     &    intent(in) :: pu,pv
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm-1),
     &    intent(in) :: sd

      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     mflx,mwdn,fdn

      REAL*8, dimension(nmom,im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     fmomdn

      INTEGER :: I,J,L,N
      integer :: jmin_x,jmax_x
      Real*8  :: byMMA

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c halo updates
c      CALL HALO_UPDATE (grid, mma) ! not needed
c      CALL HALO_UPDATE(grid, rm) ! already done by DYNAM
      CALL HALO_UPDATE_COLUMN(grid, rmom)


C****
C**** Advect the tracer using the quadratic upstream scheme
C****
      do j=j_0,j_1
      do i=1,im
        mwdn(i,j) = 0.
        fdn(i,j) = 0.
        fmomdn(:,i,j) = 0.
      enddo
      enddo

      DO L=1,LM

C**** ACCUMULATE MASS FLUXES FOR TRACERS and Q (halo included)
        MUs(:,:,L) = MUs(:,:,L)+PU(:,:,L)*DT
        MVs(:,:,L) = MVs(:,:,L)+PV(:,:,L)*DT
        IF (L.LT.LM) MWs(:,:,L) = MWs(:,:,L)+SD(:,:,L)*DT

C****
C**** convert from concentration to mass units
C****
        DO J=J_0S-1,J_1S+1
        DO I=1,IMAXJ(J)
          RM(I,J,L)=RM(I,J,L)*MMA(I,J,L)
          RMOM(:,I,J,L)=RMOM(:,I,J,L)*MMA(I,J,L)
        enddo
        enddo

        jmin_x = max(2,j_0h); jmax_x = min(jm-1,j_1h)
        do j=jmin_x,jmax_x
          mflx(:,j)=pu(:,j,l)*(.5*dt)
        enddo
        CALL AADVTX3(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &       MMA(1,j_0h,l),MFLX,jmin_x,jmax_x)

c include halo lat to the south
        mflx(:,J_0S-1:J_1S)=pv(:,J_0S:J_1S+1,l)*dt
        if (HAVE_NORTH_POLE) mflx(:,jm)=0.
        CALL AADVTY3(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &       MMA(1,j_0h,l),MFLX)

        if(l.gt.1) then
          do j=j_0,j_1
          do i=1,imaxj(j)
            MFLX(i,j)=SD(i,j,L-1)*(-DT)
c high vertical resolution near steep topography needs courz checks
            if(     mflx(i,j).gt.mma(i,j,l-1)) then
c              write(6,'(a6,3i4,f6.2)')
c     &             'courz ',i,j,l-1,mflx(i,j)/mma(i,j,l-1)
              rmom(zmoms,i,j,l-1) = 0.
            elseif(-mflx(i,j).gt.mma(i,j,l  )) then
c              write(6,'(a6,3i4,f6.2)')
c     &             'courz ',i,j,l-1,mflx(i,j)/mma(i,j,l)
              rmom(zmoms,i,j,l) = 0.
            endif
          enddo
          enddo

          CALL AADVTZ3(RM(1,j_0h,l-1),RMOM(1,1,j_0h,l-1),
     &         MMA(1,j_0h,l-1),MFLX,mwdn,fdn,fmomdn)

          jmin_x = j_0s; jmax_x = j_1s
          do j=jmin_x,jmax_x
            mflx(:,j)=pu(:,j,l-1)*(.5*dt)
          enddo
          CALL AADVTX3(RM(1,j_0h,l-1),RMOM(1,1,j_0h,l-1),
     &         MMA(1,j_0h,l-1),MFLX,jmin_x,jmax_x)

          DO J=J_0,J_1
          DO I=1,IM
            byMMA = 1 / MMA(I,J,L-1)
            RM(I,J,L-1) = RM(I,J,L-1)*byMMA
            RMOM(:,I,J,L-1) = RMOM(:,I,J,L-1)*byMMA
            rz(i,j,l-1)=rmom(mz,i,j,l-1)
          enddo
          enddo

        endif

      ENDDO ! end loop over l

      l=lm
      mflx(:,j_0:j_1)=0.
      CALL AADVTZ3(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &     MMA(1,j_0h,l),MFLX,mwdn,fdn,fmomdn)
      jmin_x = j_0s; jmax_x = j_1s
      do j=jmin_x,jmax_x
        mflx(:,j)=pu(:,j,l)*(.5*dt)
      enddo
      CALL AADVTX3(RM(1,j_0h,l),RMOM(1,1,j_0h,l),
     &     MMA(1,j_0h,l),MFLX,jmin_x,jmax_x)
      DO J=J_0,J_1
      DO I=1,IM
        byMMA = 1 / MMA(I,J,L)
        RM(I,J,L) = RM(I,J,L)*byMMA
        RMOM(:,I,J,L) = RMOM(:,I,J,L)*byMMA
        rz(i,j,l)=rmom(mz,i,j,l)
      enddo
      enddo

c fill halos of updated quantities
      call halo_update(grid,rm)
      call halo_update(grid,rz,from=south)

      RETURN
      END

      subroutine aadvtx3(rm,rmom,mass,mu,jmin,jmax)
!@sum  AADVTX advection driver for x-direction
!@auth Maxwell Kelley
c****
c**** aadvtx advects tracers in the west to east direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mu (kg) = west-east mass flux, positive eastward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use DOMAIN_DECOMP_1D, only : grid, GET
      use QUSDEF
      use QUSCOM, only : im,jm
      implicit none
      integer :: jmin,jmax
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                  rm,mass,mu
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &                  rmom
      integer :: i,ip1,ii,j,ns,nstep
      real*8 :: courmax,bynstep,frac1,fracm,fw,fe,feim,amw,dm2,mold,mnew
     &     ,bymnew
      real*8, dimension(nmom) :: fmomw,fmome,fmomeim
      real*8, dimension(im) :: am,mass_i

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1, J_0S, J_1S
      call getDomainBounds(grid, J_STRT = J_0 , J_STOP=J_1,
     &             J_STRT_SKP=J_0S,J_STOP_SKP=J_1S )

      do j=jmin,jmax
c****
c**** decide how many timesteps to take
c****
      nstep=0
      courmax = 2.
      do while(courmax.gt.1. .and. nstep.lt.20)
        nstep = nstep+1
        bynstep = 1d0/real(nstep,kind=8)
        am(:) = mu(:,j)*bynstep
        mass_i(:)  = mass(:,j)
        courmax = 0.
        do ns=1,nstep
          i = im
          do ip1=1,im
            if(am(i).gt.0.) then
               courmax = max(courmax,+am(i)/mass_i(i))
            else
               courmax = max(courmax,-am(i)/mass_i(ip1))
            endif
            i = ip1
          enddo
          if(ns.lt.nstep) then
             i = im
             do ip1=1,im
                mass_i(ip1) = mass_i(ip1) + (am(i)-am(ip1))
                i = ip1
             enddo
          endif
        enddo
      enddo

c****
c**** loop over timesteps
c****
      do ns=1,nstep

c-----------------------------------------------------------
      ! calculate tracer mass flux f
c-----------------------------------------------------------
c--------------------------------------------------------------------
      ! calculate tracer fluxes of slopes and curvatures
c--------------------------------------------------------------------
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      i=im
      if(am(i).lt.0.) then  ! air mass flux is negative
        ii=1
        frac1=+1.
      else                      ! air mass flux is positive
        ii=i
        frac1=-1.
      endif
      fracm=am(i)/mass(ii,j)
      frac1=fracm+frac1
      fw=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &     (frac1+fracm)*rmom(mxx,ii,j)))
      fmomw(mx)=am(i)*(fracm*fracm*(rmom(mx,ii,j)
     &     -3.*frac1*rmom(mxx,ii,j))-3.*fw)
      fmomw(mxx)=am(i)*(am(i)*fracm**3 *rmom(mxx,ii,j)
     &     -5.*(am(i)*fw+fmomw(mx)))
      ! cross moments
      fmomw(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
      fmomw(mxy) = am(i)*(fracm*fracm*rmom(mxy,ii,j)-3.*fmomw(my))
      fmomw(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
      fmomw(mzx) = am(i)*(fracm*fracm*rmom(mzx,ii,j)-3.*fmomw(mz))
      fmomw(myy) = fracm*rmom(myy,ii,j)
      fmomw(mzz) = fracm*rmom(mzz,ii,j)
      fmomw(myz) = fracm*rmom(myz,ii,j)
      feim = fw
      fmomeim(:) = fmomw(:)
      amw = am(im)

      do i=1,im-1
         if(am(i).lt.0.) then ! air mass flux is negative
            ii=i+1
            frac1=+1.
         else                 ! air mass flux is positive
            ii=i
            frac1=-1.
         endif
         fracm=am(i)/mass(ii,j)
         frac1=fracm+frac1
         fe=fracm*(rm(ii,j)-frac1*(rmom(mx,ii,j)-
     &        (frac1+fracm)*rmom(mxx,ii,j)))
         fmome(mx)=am(i)*(fracm*fracm*(rmom(mx,ii,j)
     &        -3.*frac1*rmom(mxx,ii,j))-3.*fe)
         fmome(mxx)=am(i)*(am(i)*fracm**3 *rmom(mxx,ii,j)
     &        -5.*(am(i)*fe+fmome(mx)))
      ! cross moments
         fmome(my)  = fracm*(rmom(my,ii,j)-frac1*rmom(mxy,ii,j))
         fmome(mxy) = am(i)*(fracm*fracm*rmom(mxy,ii,j)-3.*fmome(my))
         fmome(mz)  = fracm*(rmom(mz,ii,j)-frac1*rmom(mzx,ii,j))
         fmome(mzx) = am(i)*(fracm*fracm*rmom(mzx,ii,j)-3.*fmome(mz))
         fmome(myy) = fracm*rmom(myy,ii,j)
         fmome(mzz) = fracm*rmom(mzz,ii,j)
         fmome(myz) = fracm*rmom(myz,ii,j)

         mold=mass(i,j)
         mnew=mold+amw-am(i)
         bymnew = 1./mnew
         dm2=amw+am(i)
         rm(i,j)=rm(i,j)+fw-fe
      !
         rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm(i,j)
     &     +mold*(fw+fe))+(fmomw(mx)-fmome(mx)))*bymnew
         rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm(i,j)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw-fe)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &        mold*(fmomw(my)+fmome(my))) +
     &        (fmomw(mxy)-fmome(mxy)))*bymnew
         rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
         rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmomw(mz)+fmome(mz))) +
     &        (fmomw(mzx)-fmome(mzx)))*bymnew
      !
         rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
         rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

         mass(i,j) = mnew

         amw = am(i)
         fw = fe
         fmomw(:) = fmome(:)

      enddo ! i

      i = im
      fe = feim
      fmome(:) = fmomeim(:)
      mold=mass(i,j)
      mnew=mold+amw-am(i)
      bymnew = 1./mnew
      dm2=amw+am(i)
      rm(i,j)=rm(i,j)+fw-fe
      !
      rmom(mx,i,j)=(rmom(mx,i,j)*mold-3.*(-dm2*rm(i,j)
     &     +mold*(fw+fe))+(fmomw(mx)-fmome(mx)))*bymnew
      rmom(mxx,i,j) = (rmom(mxx,i,j)*mold*mold
     &     +2.5*rm(i,j)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fw-fe)-fmomw(mx)
     &     -fmome(mx))+dm2*rmom(mx,i,j)*mnew)
     &     +(fmomw(mxx)-fmome(mxx))) * (bymnew*bymnew)
      ! cross moments
      rmom(my,i,j)=rmom(my,i,j)+fmomw(my)-fmome(my)
      rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(my,i,j) +
     &     mold*(fmomw(my)+fmome(my))) +
     &     (fmomw(mxy)-fmome(mxy)))*bymnew
      rmom(mz,i,j)=rmom(mz,i,j)+fmomw(mz)-fmome(mz)
      rmom(mzx,i,j)=(rmom(mzx,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &     mold*(fmomw(mz)+fmome(mz))) +
     &     (fmomw(mzx)-fmome(mzx)))*bymnew
      !
      rmom(myy,i,j)=rmom(myy,i,j)+fmomw(myy)-fmome(myy)
      rmom(mzz,i,j)=rmom(mzz,i,j)+fmomw(mzz)-fmome(mzz)
      rmom(myz,i,j)=rmom(myz,i,j)+fmomw(myz)-fmome(myz)

      mass(i,j) = mnew

      enddo ! ns
      enddo ! j

      return
c****
      end subroutine aadvtx3

      subroutine aadvty3(rm,rmom,mass,mv)
!@sum  AADVTY advection driver for y-direction
!@auth Maxwell Kelley
c****
c**** aadvty advects tracers in the south to north direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg) = north-south mass flux, positive northward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use DOMAIN_DECOMP_1D, only : grid, get
      use QUSDEF
      use QUSCOM, only : im,jm,byim
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &                  rm,mass,mv
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:
     &                          grid%J_STOP_HALO) :: rmom
      integer :: i,j,jj
      REAL*8 :: m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np
      real*8, dimension(im) :: mvj,fs
      real*8, dimension(nmom,im) :: fmoms
      real*8, dimension(nmom) :: fmomn
      real*8 :: frac1,fracm,fn,mold,mnew,bymnew,dm2

c****Get relevant local distributed parameters
      INTEGER J_0,J_1,J_0H,J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0,
     &               J_STOP = J_1,
     &               J_STRT_HALO = J_0H,
     &               J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** scale polar boxes to their full extent
! set horizontal moments to zero at pole
      if (HAVE_SOUTH_POLE) then
        mass(:,1)=mass(1,1)*im
        m_sp = mass(1,1 )
        rm(:,1)=rm(1,1)*im
        rm_sp = rm(1,1 )
        rmom(zomoms,1,1 )=rmom(zomoms,1,1 )*im
        rmom(ihmoms,1,1) = 0.
        do i=2,im
           rmom(:,i,1 )=rmom(:,1,1 )
        enddo
        rzm_sp  = rmom(mz ,1,1 )
        rzzm_sp = rmom(mzz,1,1 )
      end if                       !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm)=mass(1,jm)*im
        m_np = mass(1,jm)
        rm(:,jm)=rm(1,jm)*im
        rm_np = rm(1,jm)
        rmom(zomoms,1,jm)=rmom(zomoms,1,jm)*im
        rmom(ihmoms,1,jm) = 0.
        do i=2,im
           rmom(:,i,jm)=rmom(:,1,jm)
        enddo
        rzm_np  = rmom(mz ,1,jm)
        rzzm_np = rmom(mzz,1,jm)
      end if                       !NORTH POLE

c-----------------------------------------------------------
      ! calculate tracer mass flux f
      ! and fluxes of slopes and curvatures fmom
      ! update tracer mass, moments of tracer mass, air mass distribution
c--------------------------------------------------------------------
      if(have_north_pole) then
        mv(:,jm) = 0.
      endif
      if(have_south_pole) then
        mvj(:) = 0.
        fs(:) = 0.
        fmoms(:,:) = 0.
      else
        j=j_0-1
        do i=1,im
          if(mv(i,j).lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                  ! air mass flux is positive
            jj=j
            frac1=-1.
          endif
          fracm=mv(i,j)/mass(i,jj)
          frac1=fracm+frac1
          fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &         (frac1+fracm)*rmom(myy,i,jj)))
          fmomn(my)=mv(i,j)*(fracm*fracm*(rmom(my,i,jj)
     &         -3.*frac1*rmom(myy,i,jj))-3.*fn)
          fmomn(myy)=mv(i,j)*(mv(i,j)*fracm**3 *rmom(myy,i,jj)
     &         -5.*(mv(i,j)*fn+fmomn(my)))
          fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
          fmomn(myz) = mv(i,j)*
     &         (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
          fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
          fmomn(mxy) = mv(i,j)*
     &         (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
          fmomn(mzz) = fracm*rmom(mzz,i,jj)
          fmomn(mxx) = fracm*rmom(mxx,i,jj)
          fmomn(mzx) = fracm*rmom(mzx,i,jj)
          mvj(i) = mv(i,j)
          fs(i) = fn
          fmoms(:,i) = fmomn(:)
        enddo ! i
      endif
      do j=j_0,j_1
      do i=1,im
         if(mv(i,j).lt.0.) then ! air mass flux is negative
            jj=j+1
            frac1=+1.
          else                   ! air mass flux is positive
            jj=j
            frac1=-1.
         endif
         fracm=mv(i,j)/mass(i,jj)
         frac1=fracm+frac1
         fn=fracm*(rm(i,jj)-frac1*(rmom(my,i,jj)-
     &        (frac1+fracm)*rmom(myy,i,jj)))
         fmomn(my)=mv(i,j)*(fracm*fracm*(rmom(my,i,jj)
     &        -3.*frac1*rmom(myy,i,jj))-3.*fn)
         fmomn(myy)=mv(i,j)*(mv(i,j)*fracm**3 *rmom(myy,i,jj)
     &        -5.*(mv(i,j)*fn+fmomn(my)))
         fmomn(mz)  = fracm*(rmom(mz,i,jj)-frac1*rmom(myz,i,jj))
         fmomn(myz) = mv(i,j)*
     &        (fracm*fracm*rmom(myz,i,jj)-3.*fmomn(mz))
         fmomn(mx)  = fracm*(rmom(mx,i,jj)-frac1*rmom(mxy,i,jj))
         fmomn(mxy) = mv(i,j)*
     &        (fracm*fracm*rmom(mxy,i,jj)-3.*fmomn(mx))
         fmomn(mzz) = fracm*rmom(mzz,i,jj)
         fmomn(mxx) = fracm*rmom(mxx,i,jj)
         fmomn(mzx) = fracm*rmom(mzx,i,jj)

         mold=mass(i,j)
         mnew=mold+mvj(i)-mv(i,j)
         bymnew = 1./mnew
         dm2=mvj(i)+mv(i,j)
         rm(i,j)=rm(i,j)+fs(i)-fn
      !
         rmom(my,i,j)=(rmom(my,i,j)*mold-3.*(-dm2*rm(i,j)
     &     +mold*(fs(i)+fn))+(fmoms(my,i)-fmomn(my)))*bymnew
         rmom(myy,i,j) = (rmom(myy,i,j)*mold*mold
     &     +2.5*rm(i,j)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fs(i)-fn)-fmoms(my,i)
     &     -fmomn(my))+dm2*rmom(my,i,j)*mnew)
     &     +(fmoms(myy,i)-fmomn(myy))) * (bymnew*bymnew)
      ! cross moments
         rmom(mz,i,j)=rmom(mz,i,j)+fmoms(mz,i)-fmomn(mz)
         rmom(myz,i,j)=(rmom(myz,i,j)*mold-3.*(-dm2*rmom(mz,i,j) +
     &        mold*(fmoms(mz,i)+fmomn(mz))) +
     &        (fmoms(myz,i)-fmomn(myz)))*bymnew
         rmom(mx,i,j)=rmom(mx,i,j)+fmoms(mx,i)-fmomn(mx)
         rmom(mxy,i,j)=(rmom(mxy,i,j)*mold-3.*(-dm2*rmom(mx,i,j) +
     &        mold*(fmoms(mx,i)+fmomn(mx))) +
     &        (fmoms(mxy,i)-fmomn(mxy)))*bymnew
      !
         rmom(mzz,i,j)=rmom(mzz,i,j)+fmoms(mzz,i)-fmomn(mzz)
         rmom(mxx,i,j)=rmom(mxx,i,j)+fmoms(mxx,i)-fmomn(mxx)
         rmom(mzx,i,j)=rmom(mzx,i,j)+fmoms(mzx,i)-fmomn(mzx)

         mass(i,j) = mnew

         mvj(i) = mv(i,j)
         fs(i) = fn
         fmoms(:,i) = fmomn(:)

      enddo ! i
      enddo ! j

c**** average and unscale polar boxes
! horizontal moments are zero at pole
      if (HAVE_SOUTH_POLE) then
        mass(1,1 ) = (m_sp + sum(mass(:,1 )-m_sp))*byim
        rm(1,1 ) = (rm_sp + sum(rm(:,1 )-rm_sp))*byim
        rmom(mz ,1,1 ) = (rzm_sp +
     *       sum(rmom(mz ,:,1 )-rzm_sp ))*byim
        rmom(mzz,1,1 ) = (rzzm_sp+
     *       sum(rmom(mzz,:,1 )-rzzm_sp))*byim
        rmom(ihmoms,1,1) = 0
      end if   !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(1,jm) = (m_np + sum(mass(:,jm)-m_np))*byim
        rm(1,jm) = (rm_np + sum(rm(:,jm)-rm_np))*byim
        rmom(mz ,1,jm) = (rzm_np +
     &       sum(rmom(mz ,:,jm)-rzm_np ))*byim
        rmom(mzz,1,jm) = (rzzm_np+
     &       sum(rmom(mzz,:,jm)-rzzm_np))*byim
        rmom(ihmoms,1,jm) = 0.
      end if  !NORTH POLE

      return
      end subroutine aadvty3

      subroutine aadvtz3(rm,rmom,mass,mw,mwdn,fdn,fmomdn)
!@sum  AADVTZ advection driver for z-direction
!@auth Maxwell Kelley
c****
c**** aadvtz advects tracers in the upward vertical direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mw (kg) = vertical mass flux, positive upward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      use DOMAIN_DECOMP_1D, only : grid, GET
      use QUSDEF
      use QUSCOM, only : im,jm
      USE GEOM, only : imaxj
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,2)
     &        :: rm,mass
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo,2)
     &        :: rmom
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo)
     &        :: fmomdn
      REAL*8, dimension(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     &     mw,mwdn,fdn

      real*8, dimension(nmom) :: fmomup
      integer :: i,j,l,ll
      real*8 :: frac1,fracm,fup,mold,mnew,bymnew,dm2

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1
      call getDomainBounds( grid, J_STRT=J_0 , J_STOP=J_1 )

c-----------------------------------------------------------
      ! calculate tracer mass flux f
      ! and fluxes of slopes and curvatures fmom
      ! update tracer mass, moments of tracer mass, air mass distribution
c--------------------------------------------------------------------
      do j=j_0,j_1
      do i=1,imaxj(j)
         if(mw(i,j).lt.0.) then ! air mass flux is negative
            ll=2
            frac1=+1.
          else                   ! air mass flux is positive
            ll=1
            frac1=-1.
         endif
         fracm=mw(i,j)/mass(i,j,ll)
         frac1=fracm+frac1
         fup=fracm*(rm(i,j,ll)-frac1*(rmom(mz,i,j,ll)-
     &        (frac1+fracm)*rmom(mzz,i,j,ll)))
         fmomup(mz)=mw(i,j)*(fracm*fracm*(rmom(mz,i,j,ll)
     &        -3.*frac1*rmom(mzz,i,j,ll))-3.*fup)
         fmomup(mzz)=mw(i,j)*(mw(i,j)*fracm**3 *rmom(mzz,i,j,ll)
     &        -5.*(mw(i,j)*fup+fmomup(mz)))
         fmomup(my)  = fracm*(rmom(my,i,j,ll)-frac1*rmom(myz,i,j,ll))
         fmomup(myz) = mw(i,j)*
     &        (fracm*fracm*rmom(myz,i,j,ll)-3.*fmomup(my))
         fmomup(mx)  = fracm*(rmom(mx,i,j,ll)-frac1*rmom(mzx,i,j,ll))
         fmomup(mzx) = mw(i,j)*
     &        (fracm*fracm*rmom(mzx,i,j,ll)-3.*fmomup(mx))
         fmomup(myy) = fracm*rmom(myy,i,j,ll)
         fmomup(mxx) = fracm*rmom(mxx,i,j,ll)
         fmomup(mxy) = fracm*rmom(mxy,i,j,ll)

         mold=mass(i,j,1)
         mnew=mold+mwdn(i,j)-mw(i,j)
         bymnew = 1./mnew
         dm2=mwdn(i,j)+mw(i,j)
         rm(i,j,1)=rm(i,j,1)+fdn(i,j)-fup
      !
         rmom(mz,i,j,1)=(rmom(mz,i,j,1)*mold-3.*(-dm2*rm(i,j,1)
     &     +mold*(fdn(i,j)+fup))+(fmomdn(mz,i,j)-fmomup(mz)))*bymnew
         rmom(mzz,i,j,1) = (rmom(mzz,i,j,1)*mold*mold
     &     +2.5*rm(i,j,1)*(mold*mold-mnew*mnew-3.*dm2*dm2)
     &     +5.*(mold*(mold*(fdn(i,j)-fup)-fmomdn(mz,i,j)
     &     -fmomup(mz))+dm2*rmom(mz,i,j,1)*mnew)
     &     +(fmomdn(mzz,i,j)-fmomup(mzz))) * (bymnew*bymnew)
      ! cross moments
         rmom(my,i,j,1)=rmom(my,i,j,1)+fmomdn(my,i,j)-fmomup(my)
         rmom(myz,i,j,1)=(rmom(myz,i,j,1)*mold-3.*(-dm2*rmom(my,i,j,1) +
     &        mold*(fmomdn(my,i,j)+fmomup(my))) +
     &        (fmomdn(myz,i,j)-fmomup(myz)))*bymnew
         rmom(mx,i,j,1)=rmom(mx,i,j,1)+fmomdn(mx,i,j)-fmomup(mx)
         rmom(mzx,i,j,1)=(rmom(mzx,i,j,1)*mold-3.*(-dm2*rmom(mx,i,j,1) +
     &        mold*(fmomdn(mx,i,j)+fmomup(mx))) +
     &        (fmomdn(mzx,i,j)-fmomup(mzx)))*bymnew
      !
         rmom(myy,i,j,1)=rmom(myy,i,j,1)+fmomdn(myy,i,j)-fmomup(myy)
         rmom(mxx,i,j,1)=rmom(mxx,i,j,1)+fmomdn(mxx,i,j)-fmomup(mxx)
         rmom(mxy,i,j,1)=rmom(mxy,i,j,1)+fmomdn(mxy,i,j)-fmomup(mxy)

         mass(i,j,1) = mnew

         mwdn(i,j) = mw(i,j)
         fdn(i,j) = fup
         fmomdn(:,i,j) = fmomup(:)

      enddo ! i
      enddo ! j

      if(grid%have_south_pole) then
        j=1
        do i=2,im
          mass(i,j,1) = mass(1,j,1)
          rm(i,j,1) = rm(1,j,1)
c          rmom(:,i,j,1) = rmom(:,1,j,1)
        enddo
      endif
      if(grid%have_north_pole) then
        j=jm
        do i=2,im
          mass(i,j,1) = mass(1,j,1)
          rm(i,j,1) = rm(1,j,1)
c          rmom(:,i,j,1) = rmom(:,1,j,1)
        enddo
      endif

      return
c****
      end subroutine aadvtz3

      SUBROUTINE CALC_MA(p,amp)
      USE MODEL_COM, only : im,jm,lm,ls1,dsig,psf,ptop
      USE GEOM, only : dxyp
      USE DOMAIN_DECOMP_1D, Only : grid, GET
      implicit none
      REAL*8, dimension(IM,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(IM,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: amp
      integer :: j,l
c**** Extract domain decomposition info
      INTEGER :: J_0H, J_1H
      call getDomainBounds(grid, J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

      DO L=1,LM
        IF(L.LT.LS1) THEN
          do j=max(1,J_0H),min(jm,J_1H)
            amp(:,j,l) = p(:,j)*dxyp(j)*dsig(l)
          enddo
        ELSE
          do j=max(1,J_0H),min(jm,J_1H)
            amp(:,j,l) = (psf-ptop)*dxyp(j)*dsig(l)
          enddo
        END IF
      enddo
      return
      end subroutine calc_ma

      subroutine tmom_topo_adjustments
c
c Modifies "horizontal" moments of temperature above steep topographic
c slopes to prevent excessive subsidence warming in gridcells downwind
c of high surface altitudes.
c
      use model_com, only : im,jm,lm,ls1,zatmo,t
      use dynamics, only : MUs,MVs
      use qusdef, only : mx,mxx,my,myy
      use somtq_com, only : tmom
      use domain_decomp_atm, only : grid,getDomainBounds,halo_update
      implicit none
      integer :: i,j,l
      real*8 :: zthresh,te1,te2,te1_sv,te2_sv
      INTEGER :: J_0S,J_1S

C**** define local grid
      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)

      call halo_update(grid,t)    ! already haloed ?
      call halo_update(grid,MVs)  ! already haloed ?

      do j=j_0s,j_1s
      do i=2,im-1 ! skipping wraparound for now

        zthresh = zatmo(i,j) - 12000. ! ~1200 m

c
c north-south
c
        if(zatmo(i,j-1).lt.zthresh .or. zatmo(i,j+1).lt.zthresh) then
          do l=1,ls1-1
            te1_sv = t(i,j,l)-tmom(my,i,j,l)+tmom(myy,i,j,l)
            te2_sv = t(i,j,l)+tmom(my,i,j,l)+tmom(myy,i,j,l)
            te1 = te1_sv
            te2 = te2_sv
            if(zatmo(i,j-1).lt.zthresh .and. MVs(i,j-1,l).lt.0.) then
              te1 = min(te1, .5*t(i,j,l)+.5*t(i,j-1,l) )
c              te1 = min(te1, .25*t(i,j,l)+.75*t(i,j-1,l) )
            endif
            if(zatmo(i,j+1).lt.zthresh .and. MVs(i,j  ,l).gt.0.) then
              te2 = min(te2, .5*t(i,j,l)+.5*t(i,j+1,l) )
c              te2 = min(te2, .25*t(i,j,l)+.75*t(i,j+1,l) )
            endif
            if(te1.lt.te1_sv .or. te2.lt.te2_sv) then
              tmom(my ,i,j,l) = .5*(te2-te1)
              tmom(myy,i,j,l) = .5*(te2+te1)-t(i,j,l)
            endif
          enddo
        endif
c
c east-west
c
        if(zatmo(i-1,j).lt.zthresh .or. zatmo(i+1,j).lt.zthresh) then
          do l=1,ls1-1
            te1_sv = t(i,j,l)-tmom(mx,i,j,l)+tmom(mxx,i,j,l)
            te2_sv = t(i,j,l)+tmom(mx,i,j,l)+tmom(mxx,i,j,l)
            te1 = te1_sv
            te2 = te2_sv
            if(zatmo(i-1,j).lt.zthresh .and. MUs(i-1,j,l).lt.0.) then
              te1 = min(te1, .5*t(i,j,l)+.5*t(i-1,j,l) )
c              te1 = min(te1, .25*t(i,j,l)+.75*t(i-1,j,l) )
            endif
            if(zatmo(i+1,j).lt.zthresh .and. MUs(i  ,j,l).gt.0.) then
              te2 = min(te2, .5*t(i,j,l)+.5*t(i+1,j,l) )
c              te2 = min(te2, .25*t(i,j,l)+.75*t(i+1,j,l) )
            endif
            if(te1.lt.te1_sv .or. te2.lt.te2_sv) then
              tmom(mx ,i,j,l) = .5*(te2-te1)
              tmom(mxx,i,j,l) = .5*(te2+te1)-t(i,j,l)
            endif
          enddo
        endif

      enddo
      enddo

      return
      end subroutine tmom_topo_adjustments
