#include "rundeck_opts.h"

      module PBL_DRV
      use resolution, only : im,jm
      implicit none

c     input data:
!@var uocean,vocean ocean/ice velocities for use in drag calulation
      real*8 uocean,vocean
!@var evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
      real*8 :: evap_max,fr_sat
      real*8 :: cdnl(im,jm)

ccc   some variables introduced in later code, first for tracers,
ccc   then globally. hope they are not used ...
      real*8 :: psurf, trhr0


      contains

      subroutine pbl(i,j,itype,ptype)
!@sum  simple replacement for PBL
!@vers 2013/03/26
!@auth I. Aleinov/G. Russell

C    input: ZS1,TGV,TKV,QG_SAT,HEMI,DTSURF,POLE,UOCEAN,VOCEAN
C    output:US,VS,WS,WSM,WSH,TSV,QSRF,PSI,DBL,KMS,KHS,KQS,PPBL
C          ,UG,VG,WG,W2_1

      use constant, only :  rgas,grav,deltx,twopi ! ,teeny omega2,
      use atm_com, only :  t,q,u,v
 !    &     , only : im,jm,lm, t,q,u,v,p,ptop,ls1,psf,itime
      use geom, only : idij,idjj,kmaxj,rapj,cosiv,siniv !,sinp
 !     use dynamics, only : pk !,pedn  ! pmid,
 !    &    ,dpdx_by_rho,dpdy_by_rho,dpdx_by_rho_0,dpdy_by_rho_0
      use SEAICE_COM, only : si_atm
      use seaice, only : ACE1I
      use socpbl, only :  ! npbl=>n ,
     &     zs1,tgv,tkv,qg_sat,hemi,pole    ! rest local
     &     ,us,vs,ws,wsm,wsh,tsv,qsrf,psi,dbl,kms,khs,kqs
     &     ,ug,vg,wg,zgs,w2_1,wint
      use pblcom
ccc for simple model:
      use atm_com, only : MA,pedn,pek,pk
      use QUSDEF, only : mz
      use SOMTQ_COM, only : tmom,qmom
      implicit none

      integer, intent(in) :: i,j  !@var i,j grid point
      integer, intent(in) :: itype  !@var itype surface type
      real*8, intent(in) :: ptype  !@var ptype percent surface type

ccc local vars
ccc model parameters
      real*8, parameter :: S1byG1=.57735d0
      real*8, parameter :: cdma=2.,cdmc=4., cdha=2.,cdhc=4.!, cdqa=2.,cdqc=4.
      real*8, parameter :: wsbyw1=.75d0, ciax=.3d0
ccc input from other parts of the model
      real*8 pg, ps, t1, dtdz, q1, dqdz, u1, v1
ccc temporary variables
!@var tps potential temperature at surface
!@var tks temperature at surface in K
!@var tvs virtual temperature at surface
!@var tvg virtual temperature at ground
      real*8 tps,tks,tvs,tvg, qs,wssq,cdn
      real*8 betas,betag,rigs,cdm,cdh,cia,cosc,sinc

  !    real*8 cm,ch,cq ! should be computed I.A.

      integer k

ccc input vars:
      ! zs1 - not used
      tvg = tgv
      ! tkv - not used
      ! qg_sat - not used


ccc need for simple model: pg, ma1, t1, dtdz, q1, dqdz, u1, v1

!@var pg pressure at ground (Pa)
      pg = pedn(1,i,j) * 100.d0

!@var ps pressure at surface (Pa)
      ps = pg - grav*MA(1,i,j)*.5d0*(1.d0-S1byG1)

      q1 = q(i,j,1)
      dqdz = qmom(mz,i,j,1)
      qs = q1 - dqdz*S1byG1
      !qs = q1
      qs = max(qs,0.d0)

      t1 = t(i,j,1)
      dtdz = tmom(mz,i,j,1)
      tps = t1 - dtdz*S1byG1
      tks = tps*pek(1,i,j)
      tvs = tks*(1.+qs*deltx)

!@var u1 x-velocity on top of pbl
!@var v1 y-velocity on top of pbl
      u1=0.d0 ; v1=0.d0
      ! pole and hemi are determined before pbl is called
      if (pole) then
        do k=1,kmaxj(j)
          u1 = u1 + rapj(k,j)*(u(idij(k,i,j),idjj(k,j),1)*cosiv(k) -
     2                        hemi*v(idij(k,i,j),idjj(k,j),1)*siniv(k))
          v1 = v1 + rapj(k,j)*(v(idij(k,i,j),idjj(k,j),1)*cosiv(k) +
     2                        hemi*u(idij(k,i,j),idjj(k,j),1)*siniv(k))
        end do
      else
        do k=1,kmaxj(j)
          u1 = u1 + u(idij(k,i,j),idjj(k,j),1)*rapj(k,j)
          v1 = v1 + v(idij(k,i,j),idjj(k,j),1)*rapj(k,j)
        end do
      endif

ccc wind speed:
      wssq = (u1*u1 + v1*v1)*wsbyw1*wsbyw1
      ws = sqrt(wssq)

ccc compute cdn
      select case(itype)
      case(1) ! Ocean
        CDN   = .0009d0 + .0000025d0*WSSQ
      case(2) ! Ocean ice
        CDN   = .00095417d0 +.0000005d0*ACE1I+.0000005d0*si_atm%MSI(I,J)
      case(3,4) ! Land ice, Land, snow
        CDN = CDNL(I,J)
      case default
        call stop_model("wrong itype in call to pbl",255)
      
      end select

ccc compute specific volumes
      betas = rgas*tvs/pg  ! ??? shouldn't be PS ?
      betag = rgas*tvg/pg

      if ( betag > betas ) then ! unstable stratification
        rigs = 0.
        cdm  = cdn*cdma
        cdh  = cdn*cdha
        cia  = twopi*.0625/(1.+ws*ciax)
      else ! atmosphere is stable with respect to the ground
        rigs = (pg-ps)*(betas-betag) / (wssq+.1d0)
        cdm  = cdn*cdma / (1.+cdmc*rigs*rigs)
        cdh  = cdn*cdha / (1.+cdhc*rigs*rigs)
        cia  = twopi*(.09375-.03125/(1.+4.*rigs*rigs)) / (1.+ws*ciax)
      endif

      cosc = cos(cia)
      sinc = sin(cia)*hemi
      us   = (u1*cosc - v1*sinc)*wsbyw1
      vs   = (v1*cosc + u1*sinc)*wsbyw1

cccc ***************  output ************

      wsm = sqrt((us-uocean)**2+(vs-vocean)**2)
      wsh = wsm
      tsv = tvs
      qsrf = qs

      ! have no idea what this variable is...
      ! hope this hack is correct I.A.
      wint = wsm
      
      !diag psi,
      !diag dbl,
      !kms,
      !kqs
      !diag ug,vg,wg,
      ! w2_1 used to set w2gcm, but w2gcm is not used
      w2_1 = 666.

ccc hack to remove dependence on these numbers
      khs = 1.d0

      !print *,i,j,itype, t(i,j,1)*pek(1,i,j), qs, tsv/(1.+qs*deltx)

      cmgs(i,j,itype)=cdm
      chgs(i,j,itype)=cdh
      cqgs(i,j,itype)=cdh
      ipbl(i,j,itype)=1  ! ipbl is used in subroutine init_pbl

      !wg    =sqrt(ug*ug+vg*vg)

      !psitop=atan2(vg,ug+teeny)
      !psisrf=atan2(vs,us+teeny)
      !psi   =psisrf-psitop
      ustar_pbl(i,j,itype)=666
      !ts=tsv/(1.+qsrf*deltx)
      wsavg(i,j)=wsavg(i,j)+ws*ptype
      tsavg(i,j)=tsavg(i,j)+tks*ptype
      if(itype.ne.4) qsavg(i,j)=qsavg(i,j)+qsrf*ptype
      usavg(i,j)=usavg(i,j)+us*ptype
      vsavg(i,j)=vsavg(i,j)+vs*ptype
      tauavg(i,j)=tauavg(i,j)+cdm*ws*ws*ptype

      !uflux(I,J)=uflux(I,J)+ufluxs*PTYPE
      !vflux(I,J)=vflux(I,J)+vfluxs*PTYPE
      !tflux(I,J)=tflux(I,J)+tfluxs*PTYPE
      !qflux(I,J)=qflux(I,J)+qfluxs*PTYPE

      tgvAVG(I,J)=tgvAVG(I,J)+tgv*PTYPE
      qgAVG(I,J)=qgAVG(I,J)+qg_sat*PTYPE
      w2_l1(I,J)=w2_l1(I,J)+w2_1*PTYPE

      RETURN
      END SUBROUTINE PBL


      end module PBL_DRV


      subroutine init_pbl(inipbl,istart)
c -------------------------------------------------------------
c These routines include the array ipbl which indicates if the
c  computation for a particular ITYPE was done last time step.
c Sets up the initialization of wind, temperature, and moisture
c  fields in the boundary layer. The initial values of these
c  fields are obtained by solving the static equations of the
c  Level 2 model. This is used when starting from a restart
c  file that does not have this data stored.
c -------------------------------------------------------------
      USE FILEMANAGER
      USE PBLCOM
      use pbl_drv, only : cdnl

      IMPLICIT NONE

!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
      integer :: istart
!@var iu_CDN unit number for roughness length input file
      integer :: iu_CDNL

ccc read Gary's file here
      call openunit("CDNL",iu_CDNL,.TRUE.,.true.)
      call readt (iu_CDNL,0,cdnl,im*jm,cdnl,1)
      call closeunit(iu_CDNL)

      return
      end subroutine init_pbl


      subroutine loadbl
!@sum loadbl initiallise boundary layer calc each surface time step
!@auth Ye Cheng
      USE RESOLUTION, only : jm
      USE GEOM, only : imaxj
      USE PBLCOM, only :
     *     wsavg,tsavg,qsavg,usavg,vsavg,tauavg
     &     ,uflux,vflux,tflux,qflux,tgvavg,qgavg,w2_l1
      IMPLICIT NONE
      integer i,j  !@var i,j,iter,lpbl loop variable

      do j=1,jm
      do i=1,imaxj(j)

C**** initialise some pbl common variables
          WSAVG(I,J)=0.
          TSAVG(I,J)=0.
          QSAVG(I,J)=0.
          USAVG(I,J)=0.
          VSAVG(I,J)=0.
          TAUAVG(I,J)=0.
          TGVAVG(I,J)=0.
          QGAVG(I,J)=0.
          w2_l1(I,J)=0.

          uflux(I,J)=0.
          vflux(I,J)=0.
          tflux(I,J)=0.
          qflux(I,J)=0.

      end do
      end do

      return
      end subroutine loadbl


      SUBROUTINE CHECKPBL(SUBR)
!@sum  CHECKPBL Checks whether PBL data are reasonable
!@auth Original Development Team
      USE RESOLUTION, only : im,jm
      USE PBLCOM, only : wsavg,tsavg,qsavg,dclev,usavg,vsavg,tauavg
     *     ,ustar_pbl,uflux,vflux,tflux,qflux,tgvavg,qgavg,w2_l1
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

C**** Check for NaN/INF in boundary layer data
      CALL CHECK3(wsavg,IM,JM,1,SUBR,'wsavg')
      CALL CHECK3(tsavg,IM,JM,1,SUBR,'tsavg')
      CALL CHECK3(qsavg,IM,JM,1,SUBR,'qsavg')
      CALL CHECK3(dclev,IM,JM,1,SUBR,'dclev')
      CALL CHECK3(usavg,IM,JM,1,SUBR,'usavg')
      CALL CHECK3(vsavg,IM,JM,1,SUBR,'vsavg')
      CALL CHECK3(tauavg,IM,JM,1,SUBR,'tauavg')
      CALL CHECK3(ustar_pbl,IM,JM,4,SUBR,'ustar')

      CALL CHECK3(uflux,IM,JM,1,SUBR,'uflux')
      CALL CHECK3(vflux,IM,JM,1,SUBR,'vflux')
      CALL CHECK3(tflux,IM,JM,1,SUBR,'tflux')
      CALL CHECK3(qflux,IM,JM,1,SUBR,'qflux')

      CALL CHECK3(tgvavg,IM,JM,1,SUBR,'tgvavg')
      CALL CHECK3(qgavg,IM,JM,1,SUBR,'qgavg')
      CALL CHECK3(w2_l1,IM,JM,1,SUBR,'w2_l1')

      END SUBROUTINE CHECKPBL


