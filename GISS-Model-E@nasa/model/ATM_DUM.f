#include "rundeck_opts.h"

      SUBROUTINE conserv_AM(AM)
      USE DOMAIN_DECOMP_ATM, only : GRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: AM
c empty - not needed yet.
      am = 0d0
      return
      end subroutine


      SUBROUTINE conserv_KE(KE)
c calculate column sum of kinetic energy on the A grid, using A-grid winds.
      USE CONSTANT, only : mb2kg
      USE RESOLUTION, only : lm
      USE ATM_COM, only : ua=>ualij,va=>valij,pdsig
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: KE
      integer :: i,j,l,I_0,I_1,J_0,J_1
      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      call recalc_agrid_uv
      do j=j_0,j_1
      do i=i_0,i_1
        ke(i,j) = 0d0
        do l=1,lm
          ke(i,j) = ke(i,j)+(ua(l,i,j)**2+va(l,i,j)**2)*pdsig(l,i,j)
        enddo
        ke(i,j) = .5*ke(i,j)*mb2kg
      enddo
      enddo
      return
      end subroutine conserv_KE

      SUBROUTINE calc_kea_3d(kea)
c calculates .5*(u**2+v**2) on the A grid, using A-grid winds.
      USE RESOLUTION, only : lm
      USE DOMAIN_DECOMP_ATM, only :  GRID,getDomainBounds
      USE ATM_COM, only : ua=>ualij,va=>valij
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) :: KEA
      integer :: i,j,l,I_0,I_1,J_0,J_1
      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      call recalc_agrid_uv
      do j=j_0,j_1
      do i=i_0,i_1
        do l=1,lm
          kea(i,j,l)=.5*(ua(l,i,j)**2+va(l,i,j)**2)
        enddo
      enddo
      enddo
      return
      end subroutine calc_kea_3d

      subroutine recalc_agrid_uv
!@sum recalc_agrid_uv Computes a-grid u,v from native-grid u,v
!@var u   x-component of wind, D grid
!@var v   y-component of wind, D grid
!@var ua  x-component of wind, A grid, latlon orientation, l,i,j order
!@var va  y-component of wind, A grid, latlon orientation, l,i,j order
!@auth M. Kelley
      USE RESOLUTION, only : lm
      USE ATM_COM, only : u,v,ua=>ualij,va=>valij
      USE DOMAIN_DECOMP_ATM, only : grid,getDomainBounds
      use FV_StateMod, only : INTERP_DGRID_TO_AGRID
      implicit none
      real*8, dimension(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,
     &     lm) :: ua_tmp,va_tmp,ud_tmp,vd_tmp
      integer :: i,j,l,I_0,I_1,J_0,J_1
      call getDomainBounds(grid, 
     &      I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
      do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          ud_tmp(i,j,l) = u(i,j,l) ! strip halo cells from u
          vd_tmp(i,j,l) = v(i,j,l) ! strip halo cells from v
        enddo
        enddo
      enddo
      call interp_dgrid_to_agrid(ud_tmp, vd_tmp, ua_tmp, va_tmp,
     &       rotate=.true.)
      do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          ua(l,i,j) = ua_tmp(i,j,l)
          va(l,i,j) = va_tmp(i,j,l)
        enddo
        enddo
      enddo
      return
      end subroutine recalc_agrid_uv

      subroutine replicate_uv_to_agrid(u_r,v_r,k,ursp,vrsp,urnp,vrnp)
c replicate the d-grid u and v surrounding a-grid cells into
c a-grid arrays u_r, v_r containing 2 u and 2 v at each cell.
c
c           ---U_d---            --------- 
c          |    |    |          |         |
c          |  U_r(2) |          |/V_r(1)  |
c          |         |         V_d       V_d 
c          |  U_r(1) |          |  V_r(2)/|
c          |    |    |          |         |
c           ---U_d---            --------- 
c
c
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v ! u,v are u_d and v_d
      USE DOMAIN_DECOMP_ATM, only : getDomainBounds,grid,halo_update
      implicit none
      integer :: k
      real*8, dimension(k,lm,
     &     grid%i_strt_halo:grid%i_stop_halo,
     &     grid%j_strt_halo:grid%j_stop_halo) :: u_r,v_r
c ursp,vrsp,urnp,vrnp are only for compatibility with the latlon configuration
c and are not set here
      real*8, dimension(im,lm) :: ursp,vrsp,urnp,vrnp
      integer :: i,j,l,i_0,i_1,j_0,j_1
      if(k.ne.2) call stop_model('replicate_uv_to_agrid: bad k',255)
c
      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
c for now, using symmetric scalar halo_update plus trivial edge adjustments.
c in future, this could be done better using a vector version, but see
c the discussion below in avg_replicated_duv_to_vgrid()
      call halo_update(grid, u)
      call halo_update(grid, v)
c at the top edge of an odd face, -u is the halo v
      if(mod(grid%tile,2).eq.1 .and. j_1.eq.jm) then
        u(i_0:i_1,jm+1,:) = -v(i_0:i_1,jm+1,:)
      endif
c at the right edge of an even face, -v is the halo u
      if(mod(grid%tile,2).eq.0 .and. i_1.eq.im) then
        v(im+1,j_0:j_1,:) = -u(im+1,j_0:j_1,:)
      endif
      do j=j_0,j_1
      do i=i_0,i_1
      do l=1,lm
        u_r(1,l,i,j) = u(i,j  ,l)
        u_r(2,l,i,j) = u(i,j+1,l)
        v_r(1,l,i,j) = v(i,j  ,l)
        v_r(2,l,i,j) = v(i+1,j,l)
      enddo
      enddo
      enddo
      return
      end subroutine replicate_uv_to_agrid


      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
c average the a-grid tendencies du,dv of replicated d-grid winds to
c d-grid points.  the .5 averaging factor is currently applied by
c the column physics, but this will soon change
c
c     ---------
c    |         |         --------- --------- 
c    |  du(1)  |        |         |         |
c    |    |    |        |         |/dv(1)   |
c     ---U_d---         |        V_d        |  
c    |    |    |        |   dv(2)/|         |
c    |  du(2)  |        |         |         |
c    |         |         --------- --------- 
c     ---------
c
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : u,v ! u,v are u_d and v_d
      USE DOMAIN_DECOMP_ATM, only : getDomainBounds,grid,halo_update
      implicit none
      integer :: k
      real*8, dimension(k,lm,grid%i_strt_halo:grid%i_stop_halo,
     &                       grid%j_strt_halo:grid%j_stop_halo) :: du,dv
c dusp,dvsp,dunp,dvnp are only for compatibility with the latlon configuration
c and are not used here.
      real*8, dimension(im,lm) :: dusp,dvsp,dunp,dvnp
      integer :: i,j,l,i_0,i_1,j_0,j_1
      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
c for now, using symmetric scalar halo_update plus trivial edge adjustments.
c the storage pattern in du,dv prevents use of a vector halo_update here.
      call halo_update(grid, du, jdim=4)
      call halo_update(grid, dv, jdim=4)
c in the left halo of an odd face, dv is -du
      if(mod(grid%tile,2).eq.1 .and. i_0.eq.1) then
        dv(2,:,0,j_0:j_1) = -du(2,:,0,j_0:j_1)
      endif
c in the bottom halo of an even face, du is -dv
      if(mod(grid%tile,2).eq.0 .and. j_0.eq.1) then
        du(2,:,i_0:i_1,0) = -dv(2,:,i_0:i_1,0)
      endif
      do j=j_0,j_1
      do i=i_0,i_1
      do l=1,lm
        u(i,j,l) = u(i,j,l) + !.5*
     &       (du(2,l,i,j-1)+du(1,l,i,j))
        v(i,j,l) = v(i,j,l) + !.5*
     &       (dv(2,l,i-1,j)+dv(1,l,i,j))
      enddo
      enddo
      enddo
      return
      end subroutine avg_replicated_duv_to_vgrid

      subroutine regrid_atov_1d(u_a,v_a,uv1d)
c regrid an a-grid vector u_a,v_a into its d-grid components u,d
c and pack those components into a 1-dimensional array uv1d
c
c           ----------- ----------- 
c          |           |           |      U_d is stored in odd n
c          |           |           |      V_d              even
c          |           |           |
c      uv1d(n+1)   uv1d(n+3)       |
c          |           |           |  
c          |           |           |
c          |           |           |
c           --uv1d(n)-- -uv1d(n+2)- 
c
      USE DOMAIN_DECOMP_ATM, only : grid,getDomainBounds
      use FV_StateMod, only : INTERP_AGRID_TO_DGRID
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo)  :: u_a,v_a
      real*8, dimension((1+grid%i_stop-grid%i_strt)*
     &                  (1+grid%j_stop-grid%j_strt)*2) :: uv1d
      real*8, dimension(:,:), allocatable :: u_d,v_d
      integer :: i,j,n,I_0,I_1,J_0,J_1
      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
c allocate u_d, v_d with the extra row/column for interp routine
      allocate(u_d(i_0:i_1,j_0:j_1+1),v_d(i_0:i_1+1,j_0:j_1))
      call interp_agrid_to_dgrid(
     &     u_a(i_0:i_1,j_0:j_1), v_a(i_0:i_1,j_0:j_1), u_d, v_d)
      n = 0
      do j=j_0,j_1
      do i=i_0,i_1
        n = n + 1
        uv1d(n) = u_d(i,j)
        n = n + 1
        uv1d(n) = v_d(i,j)
      enddo
      enddo
      deallocate(u_d,v_d)
      return
      end subroutine regrid_atov_1d

      subroutine get_nuv(nuv)
c reports the number of u and v in the domain
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: nuv
      nuv = 2*(1+grid%i_stop-grid%i_strt)*(1+grid%j_stop-grid%j_strt)
      return
      end subroutine get_nuv

      subroutine get_uv_of_n(n,uv)
c u is associated with odd n, v with even n
      USE RESOLUTION, only : lm
      USE ATM_COM, only : u,v
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: i,j
      call get_ij_of_n(n,i,j)
      if(mod(n,2).eq.1) then
        uv(1:lm) = u(i,j,1:lm)
      else
        uv(1:lm) = v(i,j,1:lm)
      endif
      return
      end subroutine get_uv_of_n

      subroutine store_uv_of_n(n,uv)
c u is associated with odd n, v with even n
      USE RESOLUTION, only : lm
      USE ATM_COM, only : u,v
      implicit none
      integer :: n
      real*8, dimension(lm) :: uv
      integer :: i,j
      call get_ij_of_n(n,i,j)
      if(mod(n,2).eq.1) then
        u(i,j,1:lm) = uv(1:lm)
      else
        v(i,j,1:lm) = uv(1:lm)
      endif
      return
      end subroutine store_uv_of_n

      subroutine get_vpkey_of_n(n,vpkey)
c return an identifier of the spatial location of velocity index n.
c for the d-grid, every n is a different spatial location, as its
c velocities are not collocated.
      implicit none
      integer :: n,vpkey
      vpkey = n
      return
      end subroutine get_vpkey_of_n

      subroutine get_regrid_info_for_n(n,ilist,jlist,wts,nnbr)
c u is associated with odd n, v with even n
      implicit none
      integer :: n
      integer, dimension(2) :: ilist,jlist
      real*8, dimension(2) :: wts
      integer :: i,j
      integer :: nnbr
      nnbr = 2
      wts = .5
      call get_ij_of_n(n,i,j)
      if(mod(n,2).eq.1) then ! u_d
        ilist = i; jlist = (/ j-1, j /)
      else                   ! v_d
        ilist = (/ i-1, i /); jlist = j
      endif
      return
      end subroutine get_regrid_info_for_n

      subroutine get_ij_of_n(n,i,j)
      use resolution, only : im
      USE DOMAIN_DECOMP_ATM, only : GRID
      implicit none
      integer :: n
      integer :: i,j
      integer :: nv,njm1,ni
      ni = (1+grid%i_stop-grid%i_strt)
      nv = 1+(n-1)/2
      njm1 = (nv-1)/ni
      j = grid%j_strt + njm1
      i = grid%i_strt + nv - njm1*ni - 1
      return
      end subroutine get_ij_of_n


      SUBROUTINE SDRAG(DT1)
!@sum  SDRAG puts a drag on the winds in the top layers of the atmosphere
!@auth Coefficients from the Original Development Team
      USE CONSTANT, only : grav,rgas,sha
      USE RESOLUTION, only : lm
      USE RESOLUTION, only : ls1,psfmpt
      USE DYNAMICS, only : x_sdrag,csdragl,lsdrag
     &     ,ang_sdrag,Wc_Jdrag,wmax,vsdragl,dsig
      USE ATM_COM, only : u,v,p,t,pk,pdsig,pedn,ualij,valij
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, HALO_UPDATE
      IMPLICIT NONE
!@var DT1 time step (s)
      REAL*8, INTENT(IN) :: DT1
!@var LSDRAG lowest level at which SDRAG_lin is applied
C**** SDRAG_CONST is applied in stratosphere and below the SDRAG_LIN
C**** regime (but not above P_CSDRAG)
      real*8 wl,tl,rho,cdn
      integer i,j,l
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     dpdu,dpdv,du,dv
      real*8, dimension(ls1:lm,grid%i_strt_halo:grid%i_stop_halo,
     &                         grid%j_strt_halo:grid%j_stop_halo) :: x
      real*8 xjud,xmid

      integer :: i_0, i_1, j_0, j_1
      call getDomainBounds(grid, i_strt = i_0, i_stop = i_1,
     &               j_strt = j_0, j_stop = j_1)

c
c calculate the slowdown factor X on the A grid
c
c      call recalc_agrid_uv ! not needed

      DO J=J_0,J_1
      DO I=I_0,I_1
        DO L=LS1,LM
          TL=T(I,J,L)*PK(L,I,J)
          RHO=PEDN(L+1,I,J)/(RGAS*TL)
          WL=SQRT(UALIJ(L,I,J)**2 +VALIJ(L,I,J)**2)
          xjud=1.
          if(Wc_JDRAG.gt.0.) xjud=(Wc_JDRAG/(Wc_JDRAG+min(WL,wmax)))**2
C**** WL is restricted to Wmax by adjusting X, if necessary;
C**** the following is equivalent to first reducing (U,V), if necessary,
C**** then finding the drag and applying it to the reduced winds
          IF(L.ge.LSDRAG) then
            CDN=(X_SDRAG(1)+X_SDRAG(2)*min(WL,wmax))*xjud
          else
            CDN=CSDRAGl(l)*xjud
          endif
          X(L,I,J)=DT1*RHO*CDN*min(WL,wmax)*GRAV*VSDRAGL(L)/PDSIG(L,I,J)
          if (wl.gt.wmax) X(L,I,J) = 1. - (1.-X(L,I,J))*wmax/wl
        END DO
      END DO
      END DO

c
c Apply the slowdown factor on the D grid
c
      call halo_update(grid,x,jdim=3)
      dpdu(:,:)=0.
      dpdv(:,:)=0.
      do l=ls1,lm
      do j=j_0,j_1
      do i=i_0,i_1
        xmid = min(.5*(x(l,i,j-1)+x(l,i,j)),1d0)
        dpdu(i,j) = dpdu(i,j) + u(i,j,l)*xmid*dsig(l)
        u(i,j,l) = u(i,j,l)*(1.-xmid)
        xmid = min(.5*(x(l,i-1,j)+x(l,i,j)),1d0)
        dpdv(i,j) = dpdv(i,j) + v(i,j,l)*xmid*dsig(l)
        v(i,j,l) = v(i,j,l)*(1.-xmid)
      enddo
      enddo
      enddo

c
c Conserve the column integrals of U and V by uniformly adding
c the lost stratospheric momentum to the tropospheric layers
c
      if(ang_sdrag.eq.1) then
        do j=j_0,j_1
        do i=i_0,i_1
          du(i,j) = dpdu(i,j)*psfmpt*2./(p(i,j)+p(i,j-1))
          dv(i,j) = dpdv(i,j)*psfmpt*2./(p(i-1,j)+p(i,j))
        enddo
        enddo
        do l=1,ls1-1
        do j=j_0,j_1
        do i=i_0,i_1
          u(i,j,l) = u(i,j,l) + du(i,j)
          v(i,j,l) = v(i,j,l) + dv(i,j)
        enddo
        enddo
        enddo
      endif

      return
      end subroutine sdrag

#ifdef CALC_GWDRAG

      SUBROUTINE GWDRAG
      USE CONSTANT, only : omega,grav,sha,kapa,rgas,radius
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      USE RESOLUTION, only : lm
      USE MODEL_COM, only : dtsrc
      USE ATM_COM, only : ualij,valij,pk,pedn,t,u,v,zatmo
      USE SOMTQ_COM, only : tmom,mz
      USE CLOUDS_COM,       ONLY : AIRX,LMC   
      USE STRAT,            ONLY : GWDCOL, NM, ZVARX,ZVARY,ZWT,DEFRM   
     *     ,QGWCNV,EK_globavg,dfm_type,ldef
     *     ,pbreaktop,defthresh,pconpen,ang_gwd
      USE STRAT, only : PL,PLE,TL,THL,RHO,BVF,DL,DUT,DVT,UL,VL,
     &     DP, DTIME, BYDTIME, CORIOL, AIRXS,
     &     UR, WT, CN, USRC, DUSDIF, MU_TOP, DUGWD,
     &     LDRAG, LMC0, LMC1,
     &     MU_INC, EK,
     &     RANMTN_CELL=>RANMTN,
     &     DEFRM_CELL,ZVARX_CELL,ZVARY_CELL,ZWT_CELL,ZATMO_CELL,
     &     XLIMIT
      USE GEOM, only : axyp,sinlat2d
      USE DIAG_COM, only : aij=>aij_loc
     *     ,jl_gwfirst,jl_dudtsdif,jl_sdifcoef
     *     ,ij_gw1,ij_gw2,ij_gw3,ij_gw4,ij_gw5
     *     ,ij_gw6,ij_gw7,ij_gw8,ij_gw9
      USE RANDOM
      use cs2ll_utils, only : uv_derivs_cs_agrid
      use FV_StateMod, only : INTERP_AGRID_TO_DGRID
      IMPLICIT NONE
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,LM) :: DUA,DVA ! no halo!
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP,   ! extra j for interp
     &                  GRID%J_STRT:GRID%J_STOP+1,LM) :: DUD
      REAL*8, DIMENSION(GRID%I_STRT:GRID%I_STOP+1, ! extra i for interp
     &                  GRID%J_STRT:GRID%J_STOP,LM) :: DVD
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     RANMTN,ULDEF,VLDEF
      REAL*8, DIMENSION(LM) :: AML
      INTEGER I,J,L,N,LTOP
      REAL*8 BVFSQ,ANGM,DPT,DUANG,XY
      integer lmax_angm(nm),lp10,lp2040,lpshr
C
      INTEGER :: I_0,I_1, J_0,J_1
      integer :: nij_before_j0,nij_after_j1,nij_after_i1 ! funcs
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1, J_STRT=J_0,J_STOP=J_1)

      BYDTIME = 1./DTsrc
      DTIME = DTsrc
      EK(:) = EK_globavg ! for now, no horizontal variation
      xlimit = .1d0

c
c recalculate A-grid winds
c
      call recalc_agrid_uv

C****
C**** DEFORMATION
C****
      do j=j_0,j_1
        do i=i_0,i_1
          uldef(i,j) = ualij(ldef,i,j)
          vldef(i,j) = valij(ldef,i,j)
        enddo
      enddo
      call uv_derivs_cs_agrid(grid,dfm_type,uldef,vldef,deform=defrm)
      defrm = defrm/radius ! move division by radius to cs2ll_utils

c
c Set Random numbers for Mountain Drag intermittency
c
      CALL BURN_RANDOM(nij_before_j0(J_0))
      DO J=J_0,J_1
        CALL BURN_RANDOM((I_0-1))
        DO I=I_0,I_1
          RANMTN(I,J) = RANDU(XY)
        ENDDO
        CALL BURN_RANDOM(nij_after_i1(I_1))
      ENDDO
      CALL BURN_RANDOM((nij_after_j1(J_1)))

      DO L=1,LM
        DUA(:,:,L)=0.
        DVA(:,:,L)=0.
      ENDDO

C****
C**** LOOP OVER I,J
C****
      DO J=J_0,J_1
      DO I=I_0,I_1

      CN(:)=0.
      CORIOL=2.*omega*abs(sinlat2d(i,j))

C****
C**** CALCULATE 1D ARRAYS
C****
      Call CALC_VERT_AMP (PEDN(1,I,J),LM, AML,DP,PLE,PL)
      DO L=1,LM
        TL(L)=PK(L,I,J)*T(I,J,L)
        THL(L)=T(I,J,L)
        RHO(L)=PL(L)/(RGAS*TL(L))
        BVFSQ=2.*TMOM(MZ,I,J,L)/(DP(L)*THL(L))*GRAV*GRAV*RHO(L)
        IF (PL(L).GE..4d0) THEN
          BVF(L)=SQRT(MAX(BVFSQ,1.d-10))
        ELSE
          BVF(L)=SQRT(MAX(BVFSQ,1.d-4))
        ENDIF
        DL(L)=0.
        DUT(L)=0.
        DVT(L)=0.
        UL(L)=UALIJ(L,I,J)
        VL(L)=VALIJ(L,I,J)
      ENDDO
C**** Levels for angular momentum restoration
      do l=lm,1,-1
        if (pl(l).lt.500.) lp10 = l
        if (pl(l).lt.400.) lp2040 = l
        if (pl(l).lt.200.) lpshr = l
      enddo
      if (pl(lpshr)-ple(lpshr-1).lt.pl(lpshr)-200.) lpshr=lpshr-1
      if (pl(lp2040)-ple(lp2040-1).lt.pl(lp2040)-400.) lp2040=lp2040-1
      if (pl(lp10)-ple(lp10-1).lt.pl(lp10)-500.) lp10=lp10-1
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
      IF(QGWCNV.EQ.1) THEN
        AIRXS=AIRX(I,J)/AXYP(I,J)
        IF (AIRXS.GT.0.)  THEN
          LMC0=LMC(1,I,J)
          LMC1=LMC(2,I,J)
        ENDIF 
      ENDIF 

c
c     Call the column GWD routine
c 
      CALL GWDCOL

      IF (LDRAG.LE.LM) THEN

C**** AM conservation
      if (ang_gwd.gt.0) then    ! add in ang mom
        DO N=1,NM
          ANGM = 0.
          DO L=LDRAG-1,LM
            ANGM = ANGM - DUGWD(L,N)*DP(L)
          ENDDO
          DPT=0
          DO L=1,LMAX_ANGM(N)
            DPT=DPT+DP(L)
          ENDDO
          DUANG = ANGM/DPT
          DO L=1,LMAX_ANGM(N)
             DUT(L) = DUT(L) + DUANG
          ENDDO
        ENDDO
      endif

C****
C**** Store the U and V increments
C****
      DO L=1,LM ! LDRAG-1,LM?
        DUA(I,J,L)=DUT(L)
        DVA(I,J,L)=DVT(L)
      ENDDO

c
c     UPDATE DIAGNOSTICS
c
      AIJ(I,J,IJ_GW1)=AIJ(I,J,IJ_GW1)+MU_INC(9)*UR(9)
      AIJ(I,J,IJ_GW2)=AIJ(I,J,IJ_GW2)+MU_INC(1)*UR(1)  *WT(1)
      AIJ(I,J,IJ_GW3)=AIJ(I,J,IJ_GW3)+MU_INC(2)*UR(2)
      AIJ(I,J,IJ_GW4)=AIJ(I,J,IJ_GW4)+MU_INC(3)*UR(3)  *WT(3)
      AIJ(I,J,IJ_GW5)=AIJ(I,J,IJ_GW5)+MU_INC(7)*UR(7)  *WT(7)
      AIJ(I,J,IJ_GW6)=AIJ(I,J,IJ_GW6)+MU_INC(5)*UR(5)  *WT(5)
      AIJ(I,J,IJ_GW7)=AIJ(I,J,IJ_GW7)+CN(2)*UR(2)
      AIJ(I,J,IJ_GW8)=AIJ(I,J,IJ_GW8)+USRC
      DO N=1,NM 
        AIJ(I,J,IJ_GW9)=AIJ(I,J,IJ_GW9)+MU_TOP(N) 
      ENDDO 

      DO N=1,NM 
        DO L=1,LM 
          call inc_ajl(i,j,l,N+JL_gwFirst-1,DUGWD(L,N))
        ENDDO 
      ENDDO 

      DO L=LDRAG-1,LM
        call inc_ajl(i,j,l-1,JL_DUDTSDIF,DUSDIF(L))
        call inc_ajl(i,j,l-1,JL_SDIFCOEF,DL(L)/(BVF(L)*BVF(L)))
      ENDDO

      ENDIF ! LDRAG.LE.LM
      ENDDO ! END OF LOOP OVER I
      ENDDO ! END OF LOOP OVER J

c
c Interpolate wind increments to the D grid
c
      call interp_agrid_to_dgrid(dua,dva, dud,dvd)

c
c Update D-grid winds
c
      do l=1,lm
        do j=j_0,j_1
          do i=i_0,i_1
            u(i,j,l) = u(i,j,l) + dud(i,j,l)
            v(i,j,l) = v(i,j,l) + dvd(i,j,l)
          enddo
        enddo
      enddo

C**** conservation diagnostic
c      CALL DIAGCD (GRID,6,UT,VT,DUT3,DVT3,DT1)

      RETURN
      END SUBROUTINE GWDRAG

#endif

      module zonalmean_mod
!@auth M. Kelley
!@sum  Mini-module for longitudinal averaging of a field on a
!@+    domain-decomposed cubed-sphere grid.
      use resolution, only : im
      use domain_decomp_atm, only : grid
      use domain_decomp_1d, only : dist_grid,init_grid,sumxpe,broadcast
      use cs2ll_utils, only : cs2llint_type,
     &     init_cs2llint_type,cs2llint_ij
      implicit none
      save
      private

      public :: zonalmean_ij2ij

!@param imlon,jmlat number of lons/lats for temporary field used for
!@+     longitudinal averages.  imlon is chosen to give roughly equivalent
!@+     east-west resolution, and jmlat is chosen as double the
!@+     equivalent north-south resolution to improve the accuracy of
!@+     final interpolations to the latitudes of cubed-sphere gridcells
      integer, parameter :: imlon=4*im,jmlat=4*im
!@var lons,lats lons/lats of temporary field used for lon avgs
      real*8 :: lons(imlon), lats(jmlat)
!@var have_latlon indicates whether this PE computes any lon avgs
!@+   (true only for min(NPES,jmlat) PEs
      logical :: have_latlon
      real*8 :: dlon,dlat
      real*8 :: byimlon
!@var arr_latlon temporary lat-lon field for lon avgs
!@var jlatwt weights for interpolating lon avgs to the latitudes
!@+   of local gridcells
      real*8, dimension(:,:), allocatable :: arr_latlon,jlatwt

!@var grid_zonal grid info for temporary latlon field
      type(dist_grid) :: grid_zonal

!@var cs2llint_info used during CS -> latlon interpolations
      type(cs2llint_type) :: cs2llint_info

      contains

      subroutine zonalmean_ij2ij(arr,arr_zonal)
!@sum Computes zonal means of arr and stores the result in arr_zonal.
!@+   Uses the cs2llint_ij procedure to obtain a temporary lat-lon
!@+   instance of the field, which is then averaged in longitude and
!@+   interpolated to the latitude of each cubed sphere gridcell.
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     arr,         ! input
     &     arr_zonal    ! output
      integer :: i_0,i_1,j_0,j_1,i,j,jlat
      real*8 :: wtjlat
      real*8, dimension(jmlat) :: arr_loc,arr_jlat

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

c
c check whether grid and interp info need initializing
c
      if(.not.allocated(arr_latlon)) then
        call init_zonalmean_info
      endif

c
c interpolate to obtain latlon array
c
      call cs2llint_ij(grid,cs2llint_info,arr,arr_latlon)

c
c average in longitude and broadcast to all PEs
c
      arr_loc = 0.
      if(have_latlon) then
        do jlat=grid_zonal%js,grid_zonal%je
          arr_loc(jlat) = sum(arr_latlon(:,jlat))*byimlon
        enddo
      endif
c pack_dataj unavailable since domain_decomp_1d refuses NPES > jmlat,
      call sumxpe(arr_loc,arr_jlat) ! so using sumxpe for now
      call broadcast(grid, arr_jlat)

c
c interpolate to the latitudes of CS gridcells
c
      do j=j_0,j_1
        do i=i_0,i_1
          wtjlat = jlatwt(i,j)
          jlat = wtjlat
          wtjlat = wtjlat-jlat
          arr_zonal(i,j) =
     &         wtjlat*arr_jlat(jlat+1)+(1.-wtjlat)*arr_jlat(jlat)
        enddo
      enddo

      return
      end subroutine zonalmean_ij2ij

      subroutine init_zonalmean_info
      use constant, only : pi
      use geom, only : lat2d
      implicit none
      integer :: ilon,jlat,i,j,jn
      integer :: i_0,i_1,j_0,j_1
      integer :: i_0h,i_1h,j_0h,j_1h

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      byimlon = 1d0/real(imlon,kind=8)
      have_latlon = .true.

c
c coordinates of latlon grid on which zonal means are calculated
c
      dlon = 2.*pi*byimlon
      do ilon=1,imlon
        lons(ilon) = -pi + dlon*(real(ilon,kind=8)-.5)
      enddo
      dlat = pi/real(jmlat,kind=8)
      do jlat=1,jmlat
        lats(jlat) = -pi/2. + dlat*(real(jlat,kind=8)-.5)
      enddo

c
c initialize a latlon dist_grid
c
      if(grid%nproc .le. jmlat-2) then
        call init_grid(grid_zonal,imlon,jmlat,1)
      elseif(grid%nproc .ge. jmlat) then
! domain_decomp_1d will refuse this many PEs.  Temporary workaround:
! fill the dist_grid elements needed by init_cs2llint_type, assign
! 1 j per processor, with imaginary lats beyond the north pole
        grid_zonal%gid = grid%gid
        grid_zonal%npx = imlon
        grid_zonal%npy = jmlat
        grid_zonal%js = 1+grid%gid
        grid_zonal%je = grid_zonal%js
        grid_zonal%jsd = grid_zonal%js-1
        grid_zonal%jed = grid_zonal%je+1
        have_latlon = grid_zonal%je <= jmlat
      else
        call stop_model('jmlat prob in init_zonalmean_info',255)
      endif

c
c set up CS -> latlon interp info
c
      allocate(arr_latlon(imlon,grid_zonal%jsd:grid_zonal%jed))
      call init_cs2llint_type(grid,grid_zonal,lons,lats,1,jmlat,
     &     cs2llint_info)

c
c weights for interpolation in latitude
c
      allocate(jlatwt(i_0:i_1,j_0:j_1))
      do j=j_0,j_1
        do i=i_0,i_1
          jn = min(jmlat,1+int((lat2d(i,j)-lats(1))/dlat))
          do jlat=jn,1,-1
            if(lat2d(i,j) .ge. lats(jlat)) exit
          enddo
          jlatwt(i,j) = jlat + (lat2d(i,j)-lats(jlat))/dlat
        enddo
      enddo

      return
      end subroutine init_zonalmean_info

      end module zonalmean_mod

      subroutine compute_cp_vvel(pu,pv,sd,p)
c From horizontal mass fluxes pu/pv, computes vertical mass flux in
c a constant-pressure vertical coordinate whose pressure levels are
c the global means of each layer, and interpolates them to
c terrain-following coordinate surfaces.
c Values for the model constant-pressure layers are copied from sd.
c Outputs are placed in arrays wcp,wcpsig in the DYNAMICS module.
      use RESOLUTION, only: im,lm
      use RESOLUTION, only: ls1,ptop
      USE ATM_COM, only : pednl00
      use DOMAIN_DECOMP_ATM, only: grid, getDomainBounds, halo_update
      USE DYNAMICS, only : wcp,wcpsig,sige,dsig
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     pu,pv
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm-1) ::
     &     sd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: p
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo) ::
     &     pucp,pvcp
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     dpsdt
      real*8 :: patu,patv,wtdn
      integer :: i,j,l,ldn,lup
      real*8, dimension(lm+1) :: pecp,pedge
      real*8, dimension(lm) :: uofl,vofl,bydpcp
      integer, parameter :: lmxmax=2*lm
      integer :: ltop,lmx
      real*8, dimension(lmxmax) :: dpx
      integer, dimension(lmxmax) :: lmod,lcp
      integer :: I_0, I_1, J_0, J_1

      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP
      J_0 = GRID%J_STRT
      J_1 = GRID%J_STOP

      call halo_update(grid,pu)
      call halo_update(grid,pv)
      call halo_update(grid,p)

c at the top edge of an odd face, mv is the halo mu
      if(mod(grid%tile,2).eq.1 .and. j_1.eq.im) then
        pv(i_0:i_1,im+1,:) = pu(i_0:i_1,im+1,:)
      endif
c at the right edge of an even face, mu is the halo mv
      if(mod(grid%tile,2).eq.0 .and. i_1.eq.im) then
        pu(im+1,j_0:j_1,:) = pv(im+1,j_0:j_1,:)
      endif

c redistribute pu,pv to constant-pressure reference layers
      ltop = ls1-1
      pecp(2:lm+1) = pednl00(2:lm+1)
      pecp(1) = 1d30            ! ensure that all column mass is included
      pedge(ls1) = ptop
      bydpcp(:) = 1d0/(pecp(1:lm)-pecp(2:lm+1))

      do j=j_0,j_1
        do i=i_0,i_1+1
          patu = .5*(p(i-1,j)+p(i,j))
          do l=1,ls1-1
            pedge(l) = patu*sige(l)+ptop
            pucp(l,i,j) = 0d0
            uofl(l) = pu(i,j,l)/(patu*dsig(l))
          enddo
          call get_dx_intervals(
     &         pedge,ltop,pecp,ltop,dpx,lmod,lcp,lmx,lmxmax)
          do l=1,lmx
            pucp(lcp(l),i,j)  = pucp(lcp(l),i,j)  + dpx(l)*uofl(lmod(l))
          enddo
        enddo
      enddo
      do j=j_0,j_1+1
        do i=i_0,i_1
          patv = .5*(p(i,j-1)+p(i,j))
          do l=1,ls1-1
            pedge(l) = patv*sige(l)+ptop
            pvcp(l,i,j) = 0d0
            vofl(l) = pv(i,j,l)/(patv*dsig(l))
          enddo
          call get_dx_intervals(
     &         pedge,ltop,pecp,ltop,dpx,lmod,lcp,lmx,lmxmax)
          do l=1,lmx
            pvcp(lcp(l),i,j)  = pvcp(lcp(l),i,j)  + dpx(l)*vofl(lmod(l))
          enddo
        enddo
      enddo

c compute vertical velocity in the constant-pressure coordinate
      do l=ls1-1,lm-1
        wcp(i_0:i_1,j_0:j_1,l) = sd(i_0:i_1,j_0:j_1,l)
      enddo
      wcp(i_0:i_1,j_0:j_1,lm) = 0.
      do l=ls1-1,2,-1
        do j=j_0,j_1
          do i=i_0,i_1
            wcp(i,j,l-1) = wcp(i,j,l)
     &          +(pucp(l,i,j)-pucp(l,i+1,j))
     &          +(pvcp(l,i,j)-pvcp(l,i,j+1))
          enddo
        enddo
      enddo
      l=1
      do j=j_0,j_1
        do i=i_0,i_1
          dpsdt(i,j) = wcp(i,j,l)
     &         +(pucp(l,i,j)-pucp(l,i+1,j))
     &         +(pvcp(l,i,j)-pvcp(l,i,j+1))
        enddo
      enddo

c interpolate constant-pressure vertical velocity to sigma-layer edges
c and subtract the contribution from the tendency of surface pressure
      pecp(1:ls1) = pecp(2:ls1+1)
      do j=j_0,j_1
        do i=i_0,i_1
          lup = 2
          do l=1,ls1-2
            pedge(l) = p(i,j)*sige(l+1)+ptop
            do while(pecp(lup) .gt. pedge(l))
              lup = lup + 1
            enddo
            ldn = lup - 1
            wtdn = (pedge(l)-pecp(lup))*bydpcp(lup)
            wcpsig(i,j,l) =
     &           wtdn*wcp(i,j,ldn)+(1.-wtdn)*wcp(i,j,lup)
     &           -sige(l+1)*dpsdt(i,j)
          enddo
        enddo
      enddo
      do l=ls1-1,lm-1
        wcpsig(i_0:i_1,j_0:j_1,l) = sd(i_0:i_1,j_0:j_1,l)
      enddo
      wcpsig(i_0:i_1,j_0:j_1,lm) = 0.

      return
      end subroutine compute_cp_vvel
