#include "rundeck_opts.h"

#ifndef SCM
#ifdef CUBED_SPHERE
      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@sum  This version works for a nonorthogonal grid
!@auth M. Kelley
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE GEOM, only : ddx_ci,ddx_cj,ddy_ci,ddy_cj
      USE ATM_COM, only : pmid,pedn,pk, t, zatmo,phi,
     *                    dpdx_by_rho,dpdx_by_rho_0,
     *                    dpdy_by_rho,dpdy_by_rho_0
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, HALO_UPDATE
      IMPLICIT NONE
      REAL*8 :: by_rho1
      real*8 :: dp1i ,dpsi ,dp1j ,dpsj , dg1i ,dg1j ,dgsi ,dgsj
      real*8 :: dp1dx,dpsdx,dp1dy,dpsdy, dg1dx,dg1dy,dgsdx,dgsdy
      INTEGER I,J

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, I_0, I_1
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP


C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the primary grids

      CALL HALO_UPDATE(grid, P)
      CALL HALO_UPDATE(grid, PHI(:,:,1))
c      CALL HALO_UPDATE(grid, ZATMO)

      DO J=J_0,J_1
      DO I=I_0,I_1
        dp1i = pmid(1,i+1,j) - pmid(1,i-1,j)
        dpsi = pedn(1,i+1,j) - pedn(1,i-1,j)
        dg1i = phi(i+1,j,1)-phi(i-1,j,1)
        dgsi = zatmo(i+1,j)-zatmo(i-1,j)
        dp1j = pmid(1,i,j+1) - pmid(1,i,j-1)
        dpsj = pedn(1,i,j+1) - pedn(1,i,j-1)
        dg1j = phi(i,j+1,1)-phi(i,j-1,1)
        dgsj = zatmo(i,j+1)-zatmo(i,j-1)
        dp1dx = dp1i*ddx_ci(i,j) + dp1j*ddx_cj(i,j)
        dp1dy = dp1i*ddy_ci(i,j) + dp1j*ddy_cj(i,j)
        dpsdx = dpsi*ddx_ci(i,j) + dpsj*ddx_cj(i,j)
        dpsdy = dpsi*ddy_ci(i,j) + dpsj*ddy_cj(i,j)
        dg1dx = dg1i*ddx_ci(i,j) + dg1j*ddx_cj(i,j)
        dg1dy = dg1i*ddy_ci(i,j) + dg1j*ddy_cj(i,j)
        dgsdx = dgsi*ddx_ci(i,j) + dgsj*ddx_cj(i,j)
        dgsdy = dgsi*ddy_ci(i,j) + dgsj*ddy_cj(i,j)
        by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(pmid(1,I,J))
        DPDX_BY_RHO  (I,J) = dp1dx*by_rho1 + dg1dx
        DPDX_BY_RHO_0(I,J) = dpsdx*by_rho1 + dgsdx
        DPDY_BY_RHO(  I,J) = dp1dy*by_rho1 + dg1dy
        DPDY_BY_RHO_0(I,J) = dpsdy*by_rho1 + dgsdy
      ENDDO
      ENDDO

      return
      END SUBROUTINE PGRAD_PBL


#else
      SUBROUTINE PGRAD_PBL
!@sum  PGRAD_PBL calculates surface/layer 1 pressure gradients for pbl
!@auth Ye Cheng
C**** As this is written, it must be called after the call to CALC_AMPK
C**** after DYNAM (since it uses pk/pmid). It would be better if it used
C**** SPA and PU directly from the dynamics. (Future work).
      USE CONSTANT, only : rgas
      USE RESOLUTION, only : im,jm
      USE GEOM, only : bydyp,bydxp,cosip,sinip
      USE ATM_COM, only : pmid,pedn,pk, t, zatmo,phi,
     *                    dpdx_by_rho,dpdx_by_rho_0,
     *                    dpdy_by_rho,dpdy_by_rho_0
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only : haveLatitude
      IMPLICIT NONE
      REAL*8 by_rho1,dpx1,dpy1,dpx0,dpy0,hemi,byim
      INTEGER I,J,K,IP1,IM1,J1
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)


      byim = 1d0/real(im,kind=8)

C**** (Pressure gradient)/density at first layer and surface
C**** to be used in the PBL, at the primary grids

      ! for dPdy/rho at non-pole grids
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(grid, ZATMO, FROM=SOUTH+NORTH)

      DO I=1,IM
        DO J=J_0S,J_1S
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDY_BY_RHO(I,J) = (100*(PMID(1,I,J+1)-PMID(1,I,J-1))*by_rho1
     2         +PHI(I,J+1,1)-PHI(I,J-1,1))*BYDYP(J)*.5d0
          DPDY_BY_RHO_0(I,J) =(100*(PEDN(1,I,J+1)-PEDN(1,I,J-1))*by_rho1
     2         +ZATMO(I,J+1)-ZATMO(I,J-1))*BYDYP(J)*.5d0
        END DO
      END DO

      ! for dPdx/rho at non-pole grids

      DO J=J_0S,J_1S
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          by_rho1=(rgas*t(I,J,1)*pk(1,I,J))/(100.*pmid(1,I,J))
          DPDX_BY_RHO(I,J) = (100*(PMID(1,Ip1,J)-PMID(1,Im1,J))*by_rho1
     2         +PHI(IP1,J,1)-PHI(IM1,J,1))*BYDXP(J)*.5d0
          DPDX_BY_RHO_0(I,J) =(100*(PEDN(1,Ip1,J)-PEDN(1,Im1,J))*by_rho1
     2         +ZATMO(IP1,J)-ZATMO(IM1,J))*BYDXP(J)*.5d0
          IM1=I
          I=IP1
        END DO
      END DO

      ! at poles

      IF (haveLatitude(grid, J=1)) THEN
        hemi = -1.; J1 = 2
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,1)  =dpx1*BYIM
        DPDY_BY_RHO(1,1)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,1)=dpx0*BYIM
        DPDY_BY_RHO_0(1,1)=dpy0*BYIM
      END IF

      If (haveLatitude(grid, J=JM)) THEN
          hemi= 1.; J1=JM-1
        dpx1=0. ; dpy1=0.
        dpx0=0. ; dpy0=0.
        DO K=1,IM
          dpx1=dpx1+(DPDX_BY_RHO(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO(K,J1)*SINIP(K))
          dpy1=dpy1+(DPDY_BY_RHO(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO(K,J1)*SINIP(K))
          dpx0=dpx0+(DPDX_BY_RHO_0(K,J1)*COSIP(K)
     2         -hemi*DPDY_BY_RHO_0(K,J1)*SINIP(K))
          dpy0=dpy0+(DPDY_BY_RHO_0(K,J1)*COSIP(K)
     2         +hemi*DPDX_BY_RHO_0(K,J1)*SINIP(K))
        END DO
        DPDX_BY_RHO(1,JM)  =dpx1*BYIM
        DPDY_BY_RHO(1,JM)  =dpy1*BYIM
        DPDX_BY_RHO_0(1,JM)=dpx0*BYIM
        DPDY_BY_RHO_0(1,JM)=dpy0*BYIM
      END IF

      END SUBROUTINE PGRAD_PBL
#endif
#endif


      SUBROUTINE CALC_AMPK(LMAX)
!@sum  CALC_AMPK calculate air mass and pressure arrays
!@vers 2014/04/09
!@auth Jean Lerner/Gavin Schmidt
      USE CONSTANT, only : bygrav,kapa
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : pdsig,pmid,pk,pedn,pek,MA,byMA,MASUM
      USE DOMAIN_DECOMP_ATM, Only : grid, getDomainBounds, HALO_UPDATE
      USE FLUXES, only : atmsrf,asflx4
      IMPLICIT NONE

      INTEGER :: I,J,L,IT  !@var I,J,L  loop variables
      INTEGER, INTENT(IN) :: LMAX !@var LMAX max. level for update
      REAL*8, DIMENSION(LMAX) :: AML,PDSIGL,PMIDL
      REAL*8, DIMENSION(LMAX+1) :: PEDNL
c**** Extract domain decomposition info
      Integer :: J_0H, J_1H, I_0H, I_1H
      Call getDomainBounds (GRID, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

C**** Calculate air mass, layer pressures, P**K, and sqrt(P)
C**** Note that only layers LS1 and below vary as a function of surface
C**** pressure. Routine should be called with LMAX=LM at start, and
C**** subsequentaly with LMAX=LS1-1
C**** Note Air mass is calculated in (kg/m^2)

      Do J=J_0H,J_1H
        DO I=I_0H,I_1H

          Call CALC_VERT_AMP (PEDN(1,I,J),LMAX, AML,PDSIGL,PEDNL,PMIDL)

          DO L=1,MIN(LMAX,LM)
            PDSIG(L,I,J) = PDSIGL(L)
            PMID (L,I,J) = PMIDL (L)
            PEDN (L,I,J) = PEDNL (L)
            MA   (L,I,J) = AML   (L)
            PK   (L,I,J) = PMIDL (L)**KAPA
            PEK  (L,I,J) = PEDNL (L)**KAPA
            byMA (L,I,J) = 1 / MA(L,I,J)
          END DO
          MASUM(I,J) = Sum (MA(:,I,J))
          atmsrf%P1(I,J) = PMID(1,I,J)
          atmsrf%SRFPK(I,J) = PEK(1,I,J)
          atmsrf%AM1(I,J)   =   MA(1,I,J)
          atmsrf%BYAM1(I,J) = byMA(1,I,J)

          IF (LMAX.ge.LM) THEN
            PEDN(LM+1:LMAX+1,I,J) = PEDNL(LM+1:LMAX+1)
            PEK (LM+1:LMAX+1,I,J) = PEDN(LM+1:LMAX+1,I,J)**KAPA
          END IF
        END DO
      END DO

      atmsrf%SRFP(:,:) = PEDN(1,:,:)
      do it=1,4
        asflx4(it)%SRFP = atmsrf%SRFP
      enddo
      RETURN
      END SUBROUTINE CALC_AMPK


      Subroutine MAtoPMB
!@sum  Compute P (mb) arrays from HALOed MA (kg/m^2)
!@vers 2015/05/19
      Use CONSTANT,   Only: KAPA,KG2MB
      Use RESOLUTION, Only: LM, MTOP
#ifndef STDHYB
      Use RESOLUTION, Only: MFIXS
#endif
      Use ATM_COM,    Only: MA,MASUM,byMA, PDSIG,PMID,PEDN,PK,PEK, P
      Use FLUXES,     Only: ATMSRF,ASFLX4
      Use DOMAIN_DECOMP_ATM, Only :GRID, HALO_UPDATE_COLUMN
      Implicit None
      Integer :: L,ITYPE

      Call HALO_UPDATE_COLUMN (GRID, MA)

          MASUM(:,:) = 0
      PEDN(LM+1,:,:) = MTOP*KG2MB
       PEK(LM+1,:,:) = PEDN(LM+1,:,:)**KAPA
      Do L=LM,1,-1
          MASUM(:,:) = MA(L,:,:) + MASUM(:,:)
        PDSIG(L,:,:) = MA(L,:,:)*KG2MB
         PMID(L,:,:) = PEDN(L+1,:,:) + PDSIG(L,:,:)*.5
         PEDN(L,:,:) = PEDN(L+1,:,:) + PDSIG(L,:,:)
           PK(L,:,:) = PMID(L,:,:)**KAPA
          PEK(L,:,:) = PEDN(L,:,:)**KAPA
         byMA(L,:,:) = 1 / MA(L,:,:)  ;  EndDo

          P(:,:) = (MASUM(:,:)
#ifndef STDHYB
     &      -MFIXs
#endif
     &     )*KG2MB

      ATMSRF%  AM1(:,:) =   MA(1,:,:)
      ATMSRF%byAM1(:,:) = byMA(1,:,:)
      ATMSRF%   P1(:,:) = PMID(1,:,:)
      ATMSRF% SRFP(:,:) = PEDN(1,:,:)
      ATMSRF%SRFPK(:,:) =  PEK(1,:,:)

      Do ITYPE=1,4
         ASFLX4(ITYPE)%SRFP(:,:) = ATMSRF%SRFP(:,:)  ;  EndDo
      Return
      EndSubroutine MAtoPMB


      SUBROUTINE CALC_AMP(p,amp)
!@sum  CALC_AMP Calc. AMP: kg air*grav/100, incl. const. pressure strat
!@auth Jean Lerner/Max Kelley
      USE RESOLUTION, only : ls1=>ls1_nominal,plbot
      USE RESOLUTION, only : im,jm,lm
      USE DYNAMICS, only : dsig
      USE GEOM, only : axyp
C****
      USE DOMAIN_DECOMP_ATM, Only : grid, getDomainBounds
      implicit none
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: p
      REAL*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: amp
      integer :: j,l
c**** Extract domain decomposition info
      INTEGER :: I_0, I_1, J_0, J_1
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, J_STRT=J_0,
     &               J_STOP=J_1)

C
      DO L=1,LM
        IF(L.LT.LS1) THEN
ccc   do l=1,ls1-1
          do j=J_0,J_1
            amp(I_0:I_1,j,l) = p(I_0:I_1,j)*axyp(I_0:I_1,j)*dsig(l)
          enddo
ccc   enddo
        ELSE
ccc   do l=ls1,lm
          do j=J_0,J_1
            amp(I_0:I_1,j,l) =
     &           (plbot(1)-plbot(ls1))*axyp(I_0:I_1,j)*dsig(l)
          enddo
        END IF
ccc   enddo
      enddo
C
      return
C****
      end subroutine calc_amp

      SUBROUTINE CALC_TROP
!@sum  CALC_TROP (to calculate tropopause height and layer)
!@auth J. Lerner
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : t
      USE GEOM, only : imaxj
      USE DIAG_COM, only : aij => aij_loc, ij_ptrop, ij_ttrop
      USE ATM_COM, only : pk, pmid, PTROPO, LTROPO
#ifdef etc_subdd
     &   ,TTROPO  
#endif
      USE DOMAIN_DECOMP_ATM, Only : grid, getDomainBounds
      USE CONSTANT, only : planet_name
      IMPLICIT NONE
      INTEGER I,J,L,IERR
      REAL*8, DIMENSION(LM) :: TL
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, I_0,I_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &         HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &         HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if(trim(planet_name).ne.'Earth') then
        ! We cannot call earth-specific routine tropwmo.
        ! Set the "tropopause" height to the top of the model.
        do j=J_0,J_1
        do i=I_0,imaxj(j)
          ltropo(i,j) = lm-1
          ptropo(i,j) = pmid(lm-1,i,j)
        enddo
        enddo
      else

C**** Find WMO Definition of Tropopause to Nearest L
        do j=J_0,J_1        
        do i=I_0,imaxj(j)
          do l=1,lm
            TL(L)=T(I,J,L)*PK(L,I,J)
          end do
          CALL TROPWMO(TL,PMID(1,I,J),PK(1,I,J),PTROPO(I,J),LTROPO(I,J)
     *         ,IERR)
          IF (IERR.gt.0) print*,"TROPWMO error: ",i,j
          AIJ(I,J,IJ_PTROP)=AIJ(I,J,IJ_PTROP)+PTROPO(I,J)
          AIJ(I,J,IJ_TTROP)=AIJ(I,J,IJ_TTROP)+TL(LTROPO(I,J))
#ifdef etc_subdd
          TTROPO(I,J)=TL(LTROPO(I,J))  ! extra subdaily
#endif
        end do
        end do
      endif


      IF (have_south_pole) THEN
#ifdef etc_subdd
        TTROPO(2:IM,1) = TTROPO(1,1)  ! extra subdaily
#endif
        PTROPO(2:IM,1) = PTROPO(1,1)
        LTROPO(2:IM,1) = LTROPO(1,1)
      END IF
      IF (have_north_pole) THEN
#ifdef etc_subdd
        TTROPO(2:IM,1) = TTROPO(1,1)  ! extra subdaily
#endif
        PTROPO(2:IM,JM)= PTROPO(1,JM)
        LTROPO(2:IM,JM)= LTROPO(1,JM)
      END IF

      END SUBROUTINE CALC_TROP


      subroutine tropwmo(ptm1, papm1, pk, ptropo, ltropp,ierr)
!@sum  tropwmo calculates tropopause height according to WMO formula
!@auth D. Nodorp/T. Reichler/C. Land
!@+    GISS Modifications by Jean Lerner/Gavin Schmidt
!@alg  WMO Tropopause Definition
!@+
!@+ From A Temperature Lapse Rate Definition of the Tropopause Based on
!@+ Ozone, J. M. Roe and W. H. Jasperson, 1981
!@+
!@+ In the following discussion the lapse rate is defined as -dT/dz.
!@+
!@+ The main features of the WMO tropopause definition are as follows:
!@+ * The first tropopause (i.e., the conventional tropopause) is
!@+   defined as the lowest level at which the lapse rate decreases to 2
!@+   K/km or less, and the average lapse rate from this level to any
!@+   level within the next higher 2 km does not exceed 2 K/km.
!@+ * If above the first tropopause the average lapse rate between any
!@+   level and all higher levels within 1 km exceed 3 K/km, then a
!@+   second tropopause is defined by the same criterion as under the
!@+   statement above. This tropopause may be either within or above the
!@+   1 km layer.
!@+ * A level otherwise satisfying the definition of tropopause, but
!@+   occuring at an altitude below that of the 500 mb level will not be
!@+   designated a tropopause unless it is the only level satisfying the
!@+   definition and the average lapse rate fails to exceed 3 K/km over
!@+   at least 1 km in any higher layer.
!@+ * (GISS failsafe) Some cases occur when the lapse rate never falls
!@+   below 2 K/km. In such cases the failsafe level is that where the
!@+   lapse rate first falls below 3 K/km. If this still does not work
!@+   (ever?), the level is set to the pressure level below 30mb.
!@+
      USE RESOLUTION, only : klev=>lm
      USE CONSTANT, only : zkappa=>kapa,zzkap=>bykapa,grav,rgas
      implicit none

      real*8, intent(in), dimension(klev) :: ptm1, papm1, pk
      real*8, intent(out) :: ptropo
      integer, intent(out) :: ltropp,ierr
      real*8, dimension(klev) :: zpmk, zpm, za, zb, ztm, zdtdz
!@param zgwmo min lapse rate (* -1) needed for trop. defn. (-K/km)
!@param zgwmo2 GISS failsafe minimum lapse rate (* -1) (-K/km)
!@param zdeltaz distance to check for lapse rate changes (km)
!@param zfaktor factor for caluclating height from pressure (-rgas/grav)
!@param zplimb min pressure at which to define tropopause (mb)
      real*8, parameter :: zgwmo  = -2d-3, zgwmo2=-3d-3,
     *     zdeltaz = 2000.0, zfaktor = -GRAV/RGAS, zplimb=500.
      real*8 zptph, zp2km, zag, zbg, zasum, zaquer, zptf
      integer iplimb,iplimt, jk, jj, kcount, ltset,l
      logical ldtdz
c****
c****  2. Calculate the height of the tropopause
c****  -----------------------------------------
      ltset = -999
      ierr=0
      iplimb=1
c**** set limits based on pressure
      do jk=2,klev-1
        if (papm1(jk-1).gt.600d0) then
          iplimb=jk
        else
          if (papm1(jk).lt.30d0) exit
        end if
      end do
      iplimt=jk
c****
c****  2.1 compute dt/dz
c****  -----------------
c****       ztm  lineare Interpolation in p**kappa
c****     gamma  dt/dp = a * kappa + papm1(jx,jk)**(kappa-1.)

      do jk=iplimb+1,iplimt       ! -1 ?????
        zpmk(jk)=0.5*(pk(jk-1)+pk(jk))

        zpm(jk)=zpmk(jk)**zzkap ! p mitte

        za(jk)=(ptm1(jk-1)-ptm1(jk))/(pk(jk-1)-pk(jk))
        zb(jk) = ptm1(jk)-(za(jk)*pk(jk))

        ztm(jk)=za(jk)*zpmk(jk)+zb(jk) ! T mitte
        zdtdz(jk)=zfaktor*zkappa*za(jk)*zpmk(jk)/ztm(jk)
      end do
c****
c****  2.2 First test: valid dt/dz ?
c****  -----------------------------
c****
      do 1000 jk=iplimb+1,iplimt-1

c**** GISS failsafe test
        if (zdtdz(jk).gt.zgwmo2.and.ltset.ne.1) then
          ltropp=jk
          ltset =1
        end if
c****
        if (zdtdz(jk).gt.zgwmo .and. ! dt/dz > -2K/km
     &       zpm(jk).le.zplimb) then ! zpm not too low
          ltropp = jk
          ltset = 1
c****
c****  2.3 dtdz is valid > something in German
c****  ----------------------------------------
c****    1.lineare in p^kappa (= Dieters neue Methode)

          zag = (zdtdz(jk)-zdtdz(jk+1))/
     &         (zpmk(jk)-zpmk(jk+1)) ! a-gamma
          zbg = zdtdz(jk+1) - zag*zpmk(jk+1) ! b-gamma
          if(((zgwmo-zbg)/zag).lt.0.) then
            zptf=0.
          else
            zptf=1.
          end if
          zptph = zptf*abs((zgwmo-zbg)/zag)**zzkap
          ldtdz=zdtdz(jk+1).lt.zgwmo
          if(.not.ldtdz) zptph=zpm(jk)
c****
c****  2.4 2nd test: dt/dz above 2km must not be lower than -2K/km
c****  -----------------------------------------------------------
c****
          zp2km = zptph + zdeltaz*zpm(jk)
     &         / ztm(jk)*zfaktor ! p at ptph + 2km
          zasum = 0.0           ! zdtdz above
          kcount = 0            ! number of levels above
c****
c****  2.5 Test until pm < p2km
c****  --------------------------
c****
          do jj=jk,iplimt-1
            if(zpm(jj).gt.zptph) cycle ! doesn't happen
            if(zpm(jj).lt.zp2km) goto 2000 ! ptropo valid
            zasum = zasum+zdtdz(jj)
            kcount = kcount+1
            zaquer = zasum/float(kcount) ! dt/dz mean
            if(zaquer.le.zgwmo) goto 1000 ! dt/dz above < 2K/1000
                                          ! discard it
          end do                ! test next level
          goto 2000
        endif
 1000 continue                  ! next level
 2000 continue

      if (ltset.eq.-999) then
        ltropp=iplimt-1  ! default = last level below 30mb
        print*,"In tropwmo ltropp not set, using default: ltropp ="
     *       ,ltropp
        write(6,'(12(I4,5F10.5,/))') (l,ptm1(l),papm1(l),pk(l),zdtdz(l)
     *       ,zpm(l),l=iplimb+1,iplimt-1)
        ierr=1
      end if
      ptropo = papm1(ltropp)
c****
      return
      end subroutine tropwmo

#ifndef CUBED_SPHERE
      module zonalmean_mod
      contains
      subroutine zonalmean_ij2ij(arr,arr_zonal)
c Computes zonal means of arr and stores the result in arr_zonal.
c Lat-lon version.
      use resolution, only : im
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     arr,         ! input
     &     arr_zonal    ! output
      integer :: j
      do j=grid%j_strt,grid%j_stop
        arr_zonal(:,j)=sum(arr(1:imaxj(j),j))/real(imaxj(j),kind=8)
      enddo
      return
      end subroutine zonalmean_ij2ij
      end module zonalmean_mod
#endif

!If running SCM use dummy routines
#ifndef SCM


      function getTotalEnergy() result(totalEnergy)
!@sum  getTotalEnergy returns the sum of kinetic and potential energy.
!@auth Tom Clune (SIVO)
      use GEOM, only: AXYP, AREAG
      use DOMAIN_DECOMP_ATM, only: grid, GLOBALSUM, getDomainBounds
      REAL*8 :: totalEnergy
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     KEIJ,PEIJ,TEIJ
      INTEGER :: I,J
      integer :: I_0, I_1, J_0, J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

      call conserv_PE(PEIJ)
      call conserv_KE(KEIJ)
      DO J=J_0,J_1
      DO I=I_0,I_1
        TEIJ(I,J)= (KEIJ(I,J) + PEIJ(I,J))*AXYP(I,J)/AREAG
      ENDDO
      ENDDO

      CALL GLOBALSUM(grid, TEIJ, totalEnergy, ALL=.true.)

      end function getTotalEnergy

      subroutine addEnergyAsDiffuseHeat(deltaEnergy)
!@sum  addEnergyAsDiffuseHeat adds in energy increase as diffuse heat.
!@auth Tom Clune (SIVO)
      use CONSTANT, only: sha, mb2kg
      USE RESOLUTION, only : psf, pmtop
      USE RESOLUTION, only : lm
      use ATM_COM, only: T,PK
      use DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      implicit none
      real*8, intent(in) :: deltaEnergy

      real*8 :: ediff
      integer :: l
      integer :: I_0, I_1, J_0, J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      ediff = deltaEnergy / ((PSF-PMTOP)*SHA*mb2kg)

      do l=1,lm
        T(I_0:I_1,J_0:J_1,L)=T(I_0:I_1,J_0:J_1,L)
     &       -ediff/PK(L,I_0:I_1,J_0:J_1)
      end do

      end subroutine addEnergyAsDiffuseHeat

      SUBROUTINE DISSIP
!@sum DISSIP adds in dissipated KE (m^2/s^2) as heat locally
!@auth Gavin Schmidt
      USE ATM_COM, only : t,dke,kea,pk
      IMPLICIT NONE
C**** temporarily store latest KE in DKE array
      call calc_kea_3d(dke)
      dke(:,:,:) = dke(:,:,:) - kea(:,:,:)
      call addEnergyAsLocalHeat(DKE, T, PK)

      END SUBROUTINE DISSIP

C***** Add in dissipiated KE as heat locally
      subroutine addEnergyAsLocalHeat(deltaKE, T, PK)!, diagIndex)
!@sum  addEnergyAsLocalHeat adds in dissipated kinetic energy as heat locally.
!@sum  deltaKE is on the A grid (J/kg)
!@auth Tom Clune (SIVO)
      use CONSTANT, only: SHA
      use GEOM, only: IMAXJ
      use RESOLUTION, only: LM
      use DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     deltaKE,T
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo) :: PK
c      integer, optional, intent(in) :: diagIndex

      integer :: i, j, l
      real*8 :: ediff
      integer :: I_0, I_1, J_0, J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%i_strt
      I_1 = grid%i_stop

      DO L=1,LM
      DO J=J_0,J_1
      DO I=I_0,IMAXJ(J)
        ediff = deltaKE(I,J,L) / (SHA*PK(L,I,J))
        T(I,J,L)=T(I,J,L)-ediff
c        if (present(diagIndex)) then
c          CALL INC_AJL(I,J,L,diagIndex,-ediff)
c        end if
      END DO
      END DO
      END DO
      end subroutine addEnergyAsLocalHeat

C**** Calculate 3D vertical velocity (take MWs which has units
C**** mb*m2, needs division by physics time step)
C**** and convert to WSAVE, units of m/s):

      subroutine COMPUTE_WSAVE
      use CONSTANT, only: rgas, bygrav
      use DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      use GEOM, only: byaxyp
      USE RESOLUTION, only : im,jm,lm
      use MODEL_COM, only: DTsrc
      use ATM_COM, only: T
      use ATM_COM, only: wsave, MWs,pk,pedn
      implicit none

      integer :: i, j, l
      integer :: I_0, I_1, J_0, J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do l=1,lm-1
        do j=J_0,J_1
        do i=I_0,I_1
         wsave(i,j,l)=MWs(i,j,l)*byaxyp(i,j)*
     &   rgas*0.5*(T(i,j,l)*pk(l,i,j)+T(i,j,l+1)*
     &   pk(l+1,i,j))*bygrav/(DTsrc*pedn(l+1,i,j))
        end do
        end do
      end do

      end subroutine COMPUTE_WSAVE


      Subroutine COMPUTE_GZ (MAM,S0,SZ, GZ)
!**** Input: MAM = mean mass distribution during time step (kg/m^2)
!****      S0,SZ = potential temperature and vertical gradient (K)
!**** Output: GZ = mean geopotential of layers (m^2/s^2)
      Use CONSTANT,   Only: GRAV,RGAS,KAPA,
     *                      zK=>byKAPA,zKp1=>byKAPAp1,zKp2=>byKAPAp2
      Use RESOLUTION, Only: IM,JM,LM, MTOP
      Use ATM_COM,    Only: ZATMO
      Use GEOM,       Only: IMAXJ
      Use DOMAIN_DECOMP_ATM, Only: GRID, GetDomainBounds
      Implicit None

      Integer :: I,J,L, I1,IN,J1,JN
      Logical :: QSP,QNP  
      Real*8,Dimension(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     *                    GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: MAM
      Real*8,Dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     *                 GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: S0,SZ,GZ
      Real*8 :: HUNDREDTHeKAPA, M,PU,PKU,PKPU,PKPPU,DP,zDP,Y,X,
     *          DGZU(LM),DGZA(LM),PD,PKD,PKPD,PKPPD,GZD

      I1 = GRID%I_STRT  ;  IN = GRID%I_STOP
      J1 = GRID%J_STRT  ;  JN = GRID%J_STOP
      Call GetDomainBounds (GRID, HAVE_SOUTH_POLE = QSP,
     *                            HAVE_NORTH_POLE = QNP)
      HUNDREDTHeKAPA = .01d0**KAPA

      Do J=J1,JN  ;  Do I=I1,IMAXJ(J)
!**** Integrate pressures from the top down
         M   = MTOP
         PU  = M*GRAV
         PKU = PU**KAPA  ;  PKPU = PKU*PU  ;  PKPPU = PKPU*PU
         Do L=LM,1,-1
            DP  = MAM(L,I,J)*GRAV
            zDP = 1 / DP
            Y   = SZ(I,J,L)*2*zDP*HUNDREDTHeKAPA
            X   = S0(I,J,L)*HUNDREDTHeKAPA + Y*(PU+.5*DP)
            PD  = PU + DP
            PKD = PD**KAPA  ;  PKPD = PKD*PD  ;  PKPPD = PKPD*PD   
!           AdM = RGAS*(X*(PKD-PKU)*zK - Y*(PKPD-PKPU)*zKp1)/GRAV
            DGZU(L) = RGAS*(X*(PKD-PKU)*zK - Y*(PKPD-PKPU)*zKp1)
            DGZA(L) = RGAS*(X*(DP*PKD - (PKPD-PKPU)*zKp1)*zK -
     -                Y*(DP*PKPD - (PKPPD-PKPPU)*zKp2)*zKp1)*zDP
!           AdM(I,J,L) = DGZU(L)*byGRAV
!             P(I,J,L) = GRAV*(M + .5*MAM(L,I,J))
            M   = M + MAM(L,I,J)
            PU  = PD
            PKU = PKD  ;  PKPU = PKPD  ;  PKPPU=PKPPD  ;  EndDo
!**** Integrate altitude from the bottom up
         GZD = ZATMO(I,J)
         Do L=1,LM
            GZ(I,J,L) = GZD + DGZA(L)
!      IF(J==1.OR.J==JM) WRITE (6,*) 'COMP_GZ:',J,L,GZD,GZ(I,J,L)
            GZD = GZD + DGZU(L)  ;  EndDo  ;  EndDo  ;  EndDo

      If (QSP)  Then
         Do L=1,LM
            GZ(2:IM,1,L) = GZ(1,1,L)  ;  EndDo  ;  EndIf
      If (QNP)  Then
         Do L=1,LM
            GZ(2:IM,JM,L) = GZ(1,JM,L)  ;  EndDo  ;  EndIf
      Return
      EndSubroutine COMPUTE_GZ


!if running SCM end
#endif

      function nij_before_j0(j0)
#ifdef CUBED_SPHERE
      use resolution, only : im,jm
      use domain_decomp_atm, only : grid
#else
      use geom, only : imaxj
#endif
      implicit none
      integer :: nij_before_j0,j0
#ifdef CUBED_SPHERE
      nij_before_j0 = im*((grid%tile-1)*jm + (j0-1))
#else
      nij_before_j0 = SUM(IMAXJ(1:J0-1))
#endif
      return
      end function nij_before_j0

      function nij_after_j1(j1)
      use resolution, only : im,jm
#ifdef CUBED_SPHERE
      use domain_decomp_atm, only : grid
#else
      use geom, only : imaxj
#endif
      implicit none
      integer :: nij_after_j1,j1
#ifdef CUBED_SPHERE
      nij_after_j1 = im*((6-grid%tile)*jm + (jm-j1))
#else
      nij_after_j1 = SUM(IMAXJ(J1+1:JM))
#endif
      return
      end function nij_after_j1

      function nij_after_i1(i1)
      use resolution, only : im,jm
      implicit none
      integer :: nij_after_i1,i1
#ifdef CUBED_SPHERE
      nij_after_i1 = im-i1
#else
      nij_after_i1 = 0
#endif
      return
      end function nij_after_i1
