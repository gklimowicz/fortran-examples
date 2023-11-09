#include "rundeck_opts.h"

!#define ROUGHL_HACK

      module PBL_DRV
!@sum module PBL_DRV is to compute
!@+  the turbulent transport of momentum, heat and moisture between
!@+  the surface and the middle of the first GCM layer
!@+  to find the values of the various PBL variables at the surface.
!@+  It contains the subroutine PBL.

      use SOCPBL, only : t_pbl_args, xdelt
      use SOCPBL, only : alloc_pbl_args, dealloc_pbl_args
      implicit none

      private

      public t_pbl_args, pbl, xdelt
      public alloc_pbl_args, dealloc_pbl_args
      public dbls0,slope0,dbl_max_stable
      real*8, parameter :: dbls0=10.d0, slope0=0.5d0
      real*8, parameter :: dbl_max_stable=dbls0+slope0*500. ! meters

      contains

      SUBROUTINE PBL(I,J,IHC,ITYPE,PTYPE,pbl_args,atm)
!@sum PBL contains code common for all surface.
!@+  It calculates pbl profiles, for each surface type,
!@+  to find the values of
!@+  the various PBL variables at the surface,
!@+  and accumulates diagnostics and output.
!@+  It is called from within
!@+  the subroutine SURFCE (for itype=1, 2 and 3,
!@+  i.e., surface types ocean,
!@+  seaice and landice respectively, in SURFACE.f),
!@+  and from within the subroutine earth
!@+  (for itype=4, i.e., surface type land, in GHY_DRV.f).
!@+  Dynamic equations for the mean turbulent variables
!@+  are integrated over npbl(=8) sublayers
!@+  between the surface (sublayer 1)
!@+  and the middle of the first GCM layer (sublayer npbl),
!@+  using tridiagonal method.
!@auth Greg. Hartke/Ye Cheng
!@var DDMS downdraft mass flux in kg/(m^2 s), (i,j)
!@var TDN1 downdraft temperature in K, (i,j)
!@var QDN1 downdraft humidity in kg/kg, (i,j)

      USE EXCHANGE_TYPES
      USE CONSTANT, only :  rgas,grav,omega2,deltx,teeny,lhe,lhs
      USE CONSTANT, only : planet_name
      USE GEOM, only : sinlat2d
      USE ATM_COM, only : pk
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
      USE CLOUDS_COM, only : ddm1
      USE CLOUDS_COM, only : DDMS,TDN1,QDN1,DDML
#ifdef TRACERS_ON
      USE TRACER_COM, only : ntm,trdn1
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: trradius, trpdens, tr_mm
#endif
#endif
#ifdef TRACERS_AMP
      USE AmpTracersMetadata_mod, only: AMP_MODES_MAP, AMP_NUMB_MAP
      USE TRACER_COM, only: ntmAMPi, ntmAMPe
      USE AMP_AEROSOL, only : DIAM, AMP_dens,AMP_TR_MM
      USE AERO_SETUP,  only : CONV_DPAM_TO_DGN
#endif
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif

      use SOCPBL, only : npbl=>n, zgs, advanc
      USE PBLCOM

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: IHC  !@var Height class (or 1)
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type
      type (t_pbl_args) :: pbl_args
c     lmonin=-133.d0
      class (atmsrf_xchng_vars) :: atm
      REAL*8 Ts,tmp
      real*8 :: qsat ! external
      real*8, parameter :: S1byG1=.57735d0

#ifdef TRACERS_ON
      integer nx,n
#ifndef TRACERS_TOMAS
     *     ,nAMP
#endif
#endif

c
      REAL*8 ztop,coriol,rhosrf
      REAL*8 qtop,utop,vtop,ufluxs,vfluxs,tfluxs,qfluxs,psitop,psisrf
      INTEGER k
!@var uocean,vocean ocean/ice velocities for use in drag calulation
!@var evap_max maximal evaporation from unsaturated soil
!@var  fr_sat fraction of saturated soil
!@var ZS1    = height of the first model layer (m)
!@var TGV    = virtual potential temperature of the ground (K)
!@+            (if xdelt=0, TGV is the actual temperature)
!@var TKV    = virtual potential temperature of first model layer (K)
!@+            (if xdelt=0, TKV is the actual temperature)
!@var WS     = magnitude of the surface wind (m/s)
!@var PSI    = angular diff. btw geostrophic and surface winds (rads)
!@var WG     = magnitude of the geostrophic wind (m/s)
!@var HEMI   = 1 for northern hemisphere, -1 for southern hemisphere
!@var TG     = bulk ground temperature (K)
!@var ELHX   = latent heat for saturation humidity (J/kg)
!@var dskin  = skin-bulk SST difference (C)
!@VAR QSOL   = solar heating (W/m2)
      real*8 zs1,psi,hemi
!@var POLE   = .TRUE. if at the north or south pole, .FALSE. otherwise
c      logical pole

!**** the following is output from advance (mostly passed through pbl_args)
!@var US     = x component of surface wind, positive eastward (m/s)
!@var VS     = y component of surface wind, positive northward (m/s)
!@var WSGCM  = magnitude of the GCM surface wind - ocean currents (m/s)
!@var WSPDF  = mean surface wind calculated from PDF of wind speed (m/s)
!@var WS     = magn. of GCM surf wind - ocean curr + buoyancy + gust (m/s)
!@var TSV    = virtual potential temperature of the surface (K)
!@+            (if xdelt=0, TSV is the actual temperature)
!@var QS     = surface value of the specific moisture
!@var DBL    = boundary layer height (m)
!@var DBLS   = stable boundary layer height (m)
!@var LDBL   = layer immediately above dbl
!@var LDBLS  = layer immediately above dbls
!@var KMS    = momentum transport coefficient at ZGS (m**2/s)
!@var KHS    = heat transport coefficient at ZGS (m**2/s)
!@var KHQ    = moist transport coefficient at ZGS (m**2/s)
!@var USTAR  = friction speed (square root of momentum flux) (m/s)
!@var CM     = drag coefficient (dimensionless surface momentum flux)
!@var CH     = Stanton number   (dimensionless surface heat flux)
!@var CQ     = Dalton number    (dimensionless surface moisture flux)
!@var z0m   = roughness length for momentum,
!@+           prescribed for itype=3,4 but computed for itype=1,2 (m)
!@var z0h   = roughness length for temperature (m)
!@var z0q   = roughness length for water vapor (m)
!@var UG     = eastward component of the geostrophic wind (m/s)
!@var VG     = northward component of the geostrophic wind (m/s)
!@var MDF    = downdraft mass flux (m/s)
!@var WINT   = integrated surface wind speed over sgs wind distribution
      real*8 :: dbl,kms,kqs,cm,ch,cq,z0m,z0h,z0q,ug,vg,w2_1,mdf
!@var dtdt_gcm temp. tendency from processes other than turbulence (K/s)
      real*8 ::  dpdxr,dpdyr,dpdxr0,dpdyr0,dtdt_gcm
      real*8 ::  mdn  ! ,mup
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
#if defined(TRACERS_ON)
!@var  tr local tracer profile (passive scalars)
      real*8, dimension(npbl,pbl_args%ntx) :: tr
      real*8, dimension(ntm) :: trnradius,trndens,trnmm
      real*8 :: rts,rtsdt
#endif
      real*8 dbls,ustar,lmonin
      integer ldbls,ldbl

      pbl_args%psurf = atm%srfp(i,j)
      pbl_args%QSOL = atm%fshort(i,j)*atm%cosz1(i,j) ! solar heating

      !qtop=q(i,j,1)
      qtop = atm%q1(i,j)
      pbl_args%tkv = (atm%temp1(I,J)*(1.+qtop*xdelt))*atm%srfpk(i,j)
      pbl_args%ZS1=.5d-2*RGAS*pbl_args%TKV*atm%AM1(i,j)/atm%p1(i,j)
      pbl_args%hemi = sign(1d0,atm%lat(i,j))

      utop = atm%u1(i,j) !ua(1,i,j)
      vtop = atm%v1(i,j) !va(1,i,j)

      if(itype < 4) then
        pbl_args%evap_max = 1.
        pbl_args%fr_sat = 1.    ! entire surface is saturated
        pbl_args%qg_aver = pbl_args%qg_sat ! QG_AVER=QG_SAT
#ifdef SCM
        if( SCMopt%Qskin )then ! force skin water vapor mixing ratio
          pbl_args%qg_aver = SCMin%Qskin
        endif
#endif
        if(itype==1) then
          pbl_args%elhx = lhe
        else
          pbl_args%elhx = lhs
        endif
      endif

      if(itype > 2) then
        pbl_args%uocean = 0.
        pbl_args%vocean = 0.
      endif

      pbl_args%trhr0 = atm%flong(i,j)

#ifdef TRACERS_ON
C**** Set up tracers for PBL calculation if required
      do nx=1,pbl_args%ntx
        n=pbl_args%ntix(nx)
C**** Calculate first layer tracer concentration
          pbl_args%trtop(nx)=
     &         atm%trm1(n,i,j)*atm%byam1(i,j)
      end do
      call tracer_lower_bc(i,j,itype,pbl_args,atm)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS) 
      call dust_emission_prep(i,j,itype,pbl_args)
#endif

#endif


ccc extract data needed in driver from the pbl_args structure
      zs1 = pbl_args%zs1
      hemi = pbl_args%hemi
c      pole = pbl_args%pole

#ifdef USE_PBL_E1
      pbl_args%ddml_eq_1=.false.
#else
      pbl_args%ddml_eq_1=DDML(i,j).eq.1
#endif
      ! Redelsperger et al. 2000, eqn(13), J. Climate, 13, 402-421
      ! tprime,qprime are the pertubation of t and q due to gustiness

      ! pick up one of the following two expressions for gusti

      ! for down draft:
      if(pbl_args%ddml_eq_1) then
        mdn=max(DDMS(i,j), -0.07d0)
        pbl_args%gusti=log(1.-600.4d0*mdn-4375.*mdn*mdn)
      else
        pbl_args%gusti=0.
      endif

      ! for up draft:
      ! mup=min(DDMS(i,j), 0.1d0)
      ! pbl_args%gusti=log(1.+386.6d0*mup-1850.*mup*mup)

C        ocean and ocean ice are treated as rough surfaces
C        roughness lengths from Brutsaert for rough surfaces

      IF (ITYPE.GT.2) THEN
        Z0M=ROUGHL(I,J)           ! 30./(10.**ROUGHL(I,J))
      ENDIF
      ztop=zgs+zs1  ! zs1 is calculated before pbl is called
c     IF (pbl_args%TKV.EQ.pbl_args%TGV)
c    &     pbl_args%TGV = 1.0001d0*pbl_args%TGV

      dbl = bldep(i,j)
      ug   = ugeo(i,j) !ua(ldbl,i,j)
      vg   = vgeo(i,j) !va(ldbl,i,j)
      coriol=sinlat2d(i,j)*omega2

      upbl(:)=atm%uabl(:,i,j)
      vpbl(:)=atm%vabl(:,i,j)
      tpbl(:)=atm%tabl(:,i,j)
      qpbl(:)=atm%qabl(:,i,j)
      epbl(1:npbl-1)=atm%eabl(1:npbl-1,i,j)

#ifdef TRACERS_ON
      do nx=1,pbl_args%ntx
        tr(:,nx)=atm%trabl(:,pbl_args%ntix(nx),i,j)
      end do

      do n = 1,ntm
#ifdef TRACERS_DRYDEP
           trnradius(n) = trradius(n)
           trndens(n)   = trpdens(n)
           trnmm(n)     = tr_mm(n)
#endif
#ifdef TRACERS_AMP
       if (n.ge.ntmAMPi.and.n.le.ntmAMPe) then
         nAMP=n-ntmAMPi+1
        if(AMP_MODES_MAP(nAMP).gt.0) then
         if(DIAM(i,j,1,AMP_MODES_MAP(nAMP)).gt.0.) then
          if(AMP_NUMB_MAP(nAMP).eq. 0) then    ! Mass
        trnradius(n)=0.5*DIAM(i,j,1,AMP_MODES_MAP(nAMP))
          else                              ! Number
        trnradius(n)=0.5*DIAM(i,j,1,AMP_MODES_MAP(nAMP))
     +               *CONV_DPAM_TO_DGN(AMP_MODES_MAP(nAMP))
          endif

           call AMPtrdens(i,j,1,n)
           call AMPtrmass(i,j,1,n)

          trndens(n) =AMP_dens(i,j,1,AMP_MODES_MAP(nAMP))
          trnmm(n)   =AMP_TR_MM(i,j,1,AMP_MODES_MAP(nAMP))
        endif   
        endif   
       endif 
#endif  
      enddo
#endif

      cm=atm%cmgs(i,j)
      ch=atm%chgs(i,j)
      cq=atm%cqgs(i,j)
      dpdxr  = DPDX_BY_RHO(i,j)
      dpdyr  = DPDY_BY_RHO(i,j)
      dpdxr0 = DPDX_BY_RHO_0(i,j)
      dpdyr0 = DPDY_BY_RHO_0(i,j)

      mdf = ddm1(i,j)

!!! put some results from above to pbl_args
      pbl_args%dbl = dbl
      pbl_args%ug = ug
      pbl_args%vg = vg
      pbl_args%wg = sqrt(ug*ug+vg*vg)
      pbl_args%cm = cm
      pbl_args%ch = ch
      pbl_args%cq = cq

      ! if ddml_eq_1=.false.,
      ! i.e., either USE_PBL_E1 or DDML(i,j) is not 1,
      ! then tdns,qdns,tprime,qprime are not in use

      if (pbl_args%ddml_eq_1) then
        pbl_args%tdns=TDN1(i,j)*atm%srfpk(i,j)/pk(1,i,j)
        pbl_args%qdns=QDN1(i,j)
#ifdef TRACERS_ON
        do nx=1,pbl_args%ntx
          pbl_args%trdn1(nx)=TRDN1(pbl_args%ntix(nx),i,j)
        end do
#endif
      else
        pbl_args%tdns=0.d0
        pbl_args%qdns=0.d0
#ifdef TRACERS_ON
        pbl_args%trdn1(:)=0.
#endif
      endif

      dtdt_gcm = (pbl_args%tkv - t1_after_aturb(i,j)*
     &     atm%srfpk(i,j))/pbl_args%dtsurf
      call advanc( pbl_args,coriol,utop,vtop,qtop,ztop,mdf
     &     ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &     ,dtdt_gcm,u1_after_aturb(i,j),v1_after_aturb(i,j)
     &     ,i,j,ihc,itype
     &     ,kms,kqs,z0m,z0h,z0q,w2_1,ufluxs,vfluxs,tfluxs,qfluxs
     &     ,upbl,vpbl,tpbl,qpbl,epbl
#if defined(TRACERS_ON)
     &     ,tr,trnradius,trndens,trnmm
#endif
     &     )

      atm%uabl(:,i,j)=upbl(:)
      atm%vabl(:,i,j)=vpbl(:)
      atm%tabl(:,i,j)=tpbl(:)
      atm%qabl(:,i,j)=qpbl(:)
      atm%eabl(1:npbl-1,i,j)=epbl(1:npbl-1)
      rhosrf=100.*pbl_args%psurf/(rgas*pbl_args%tsv)
#ifdef TRACERS_ON
      do nx=1,pbl_args%ntx
        n = pbl_args%ntix(nx)
        atm%trabl(:,n,i,j)=tr(:,nx)
        atm%travg(n,i,j) = pbl_args%trs(nx)
        atm%travg_byvol(n,i,j) = pbl_args%trs(nx)*rhosrf
      end do
#ifdef TRACERS_DRYDEP
      atm%dep_vel(:,i,j)=pbl_args%dep_vel(:)
      atm%gs_vel(:,i,j)=pbl_args%gs_vel(:)
      do nx=1,pbl_args%ntx
        n = pbl_args%ntix(nx)
        rts=rhosrf*pbl_args%trs(nx)
        rtsdt=rts*pbl_args%dtsurf        ! kg*s/m^3
        atm%drydflx(n,i,j)=-rtsdt*
     &       (atm%dep_vel(n,i,j)+atm%gs_vel(n,i,j)) ! kg/m2
      enddo
#ifdef ACCMIP_LIKE_DIAGS
      atm%stomatal_dep_vel(i,j)=pbl_args%stomatal_dep_vel
#endif /* ACCMIP_LIKE_DIAGS */
#endif /* TRACERS_DRYDEP */
#endif /* TRACERS_ON */

#ifdef SCM
      if ( SCMopt%geo .and. SCMopt%ustar ) then
c**** force friction speed and surface drag coefficient
        pbl_args%ustar = SCMin%ustar
        pbl_args%cm = (SCMin%ustar/pbl_args%ws)**2
      endif
#endif

      atm%cmgs(i,j)=pbl_args%cm
      atm%chgs(i,j)=pbl_args%ch
      atm%cqgs(i,j)=pbl_args%cq
      atm%ipbl(i,j)=1  ! ipbl is used in subroutine init_pbl

      psitop=atan2(vg,ug+teeny)
      psisrf=atan2(pbl_args%vs,pbl_args%us+teeny)
      psi   =psisrf-psitop
      atm%ustar_pbl(i,j)=pbl_args%ustar
      atm%lmonin_pbl(i,j)=pbl_args%lmonin
C ******************************************************************
      TS=pbl_args%TSV/(1.+pbl_args%QSRF*xdelt)
      if(planet_name.eq.'Earth') then
      if ( ts.lt.152d0 .or. ts.gt.423d0 ) then
        write(6,*) 'PBL: Ts bad at',i,j,' itype',itype,ts
        if (ts.gt.1d3) call stop_model("PBL: Ts out of range",255)
        if (ts.lt.50d0) call stop_model("PBL: Ts out of range",255)
      end if
      endif
      atm%tsavg(i,j) = ts
      atm%qsavg(i,j) = pbl_args%qsrf
      atm%usavg(i,j) = pbl_args%us
      atm%vsavg(i,j) = pbl_args%vs
      atm%wsavg(i,j) = pbl_args%ws
      atm%rsavg(i,j) = atm%qsavg(i,j)/
     &     qsat(ts,pbl_args%elhx,atm%srfp(i,j))


      atm%TAUAVG(I,J) = pbl_args%CM*pbl_args%WS*pbl_args%WS*rhosrf
      atm%tgvAVG(I,J) = pbl_args%tgv
      atm%qgAVG(I,J) = pbl_args%qg_aver
      atm%gustiwind(i,j) = pbl_args%gusti
      atm%dblavg(i,j) = dbl
      atm%rhoavg(i,j) = rhosrf
      atm%w2_l1(I,J) = w2_1
      atm%ciaavg(i,j)  =  psi
      atm%khsavg(i,j)  =  pbl_args%khs
      atm%wspdf(i,j) = pbl_args%wspdf

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      atm%wsgcm(i,j) = pbl_args%wsgcm
      atm%wsubwd(i,j) = pbl_args%wsubwd
      atm%wsubtke(i,j) = pbl_args%wsubtke
      atm%wsubwm(i,j) = pbl_args%wsubwm
#endif

ccc put drive output data to pbl_args structure
      pbl_args%psi = psi ! maybe also should be moved to ADVANC
                         ! or completely otside of PBL* ?

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      call PBL_adiurn_dust(I,J,ITYPE,PTYPE,pbl_args,atm)
#endif

#ifdef TRACERS_ON
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS) 
      call save_dust_emission_vars(i,j,itype,pbl_args)
#endif
#endif

      RETURN
      END SUBROUTINE PBL

#ifdef TRACERS_ON
      subroutine tracer_lower_bc(i,j,itype,pbl_args,atm)
      use exchange_types
      use geom, only : byaxyp
      use OldTracer_mod, only :
     &     trname,dodrydep,tr_wd_type,nWATER,nPART,nGAS
      use tracer_com, only: n_co2n
      implicit none
      integer, intent(in) :: i,j,itype  !@var itype surface type
      type (t_pbl_args) :: pbl_args
      class (atmsrf_xchng_vars) :: atm
c
      integer :: n,nx
c

c todo: try to merge the itype==4 and itype<4 cases
      if(itype.le.3) then ! ice/water surfaces
C 
C**** Set up b.c. for tracer PBL calculation if required
        do nx=1,pbl_args%ntx
          n=pbl_args%ntix(nx)
C**** set defaults
          pbl_args%trsfac(nx)=0.
          pbl_args%trconstflx(nx)=0.
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
#ifdef TRACERS_WATER
          pbl_args%tr_evap_max(nx)=1.
C**** This distinguishes water from gases or particle
          if ( tr_wd_TYPE(n) == nWATER ) then
c          trgrnd(nx)=atmgla%gtracer(n,i,j)
C**** trsfac and trconstflx are multiplied by cq*ws and QG_SAT in PBL
            pbl_args%trsfac(nx)=1.
            pbl_args%trconstflx(nx)=atm%gtracer(n,i,j)
          else if ( tr_wd_TYPE(n) == nGAS .or.
     &           tr_wd_TYPE(n) == nPART ) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
            if(dodrydep(n)) then
              pbl_args%trsfac(nx)=1. !then multiplied by deposition velocity in PBL
#ifdef TRACERS_WATER
              pbl_args%tr_evap_max(nx)=1.d30
c            trgrnd(nx)=0.
#endif
            end if
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
C**** Now send kg/m^2/s to PBL, and divided by rho there.
#ifndef SKIP_TRACER_SRCS
            pbl_args%trconstflx(nx)=atm%trflux_prescr(n,i,j) ! kg/m^2/s
#endif /*SKIP_TRACER_SRCS*/

#ifdef TRACERS_WATER
          endif
#endif

          if (n==n_co2n) then
! transplanted from SURFACE.f:
            IF (pbl_args%ocean) THEN ! OCEAN
              !trgrnd(nx)=atm%gtracer(n,i,j)
              pbl_args%trsfac(nx)=1.
            END IF
        !need to redo this here because the previous line has changed trconstflx to zero.
        !because we have no sources. is there a better way to do this?
          !call stop_model('why 1/area in the following line?',255)
            pbl_args%trconstflx(nx)=atm%gtracer(n,i,j) * byaxyp(i,j) !kg,co2/kg,air/m2
          endif

        end do
      else
c LAND SURFACE
        do nx=1,pbl_args%ntx
          n=pbl_args%ntix(nx)
C**** set defaults
          pbl_args%trsfac(nx)=0.
          pbl_args%trconstflx(nx)=0.
#ifdef TRACERS_WATER
C**** Set surface boundary conditions for tracers depending on whether
C**** they are water or another type of tracer
C**** The select is used to distinguish water from gases or particle
! select removed because of OMP compiler bug
!        select case (tr_wd_TYPE(n))
!        case (nWATER)
          if (tr_wd_TYPE(n) .eq. nWATER) then
C**** no fractionation from ground (yet)
C**** trsfac and trconstflx are multiplied by cq*ws and QG in PBL
            pbl_args%trsfac(nx)=1.
            pbl_args%trconstflx(nx)=atm%gtracer(n,i,j)
!        case (nGAS, nPART)
          elseif(tr_wd_TYPE(n).eq.nGAS .or. tr_wd_TYPE(n).eq.nPART) then
#endif
C**** For non-water tracers (i.e. if TRACERS_WATER is not set, or there
C**** is a non-soluble tracer mixed in.)
C**** Calculate trsfac (set to zero for const flux)
#ifdef TRACERS_DRYDEP
            if(dodrydep(n)) pbl_args%trsfac(nx) = 1.
          !then multiplied by deposition velocity in PBL
#endif
C**** Calculate trconstflx (m/s * conc) (could be dependent on itype)
            pbl_args%trconstflx(nx)=atm%trflux_prescr(n,i,j) ! kg/m^2/s
#ifdef TRACERS_WATER
!        end select
          end if
#endif
        end do

#ifdef TRACERS_WATER
c**** water tracers are also flux limited
        do nx=1,pbl_args%ntx
          n=pbl_args%ntix(nx)
C       pbl_args%tr_evap_max(nx) = evap_max * trsoil_rat(nx)
          pbl_args%tr_evap_max(nx)=pbl_args%evap_max*atm%gtracer(n,i,j)
#ifdef TRACERS_DRYDEP
          if(dodrydep(n)) pbl_args%tr_evap_max(nx) = 1.d30
#endif
        end do
#endif
      endif ! itype check

      return
      end subroutine tracer_lower_bc
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS) 
      subroutine dust_emission_prep(i,j,itype,pbl_args)
      use constant, only : by3
      use fluxes, only : pprec,pevap
      use trdust_mod, only : hbaij,ricntd
      use clouds_com, only : ddml
      implicit none
      integer, intent(in) :: i,j  !@var i,j grid point
      integer, intent(in) :: itype  !@var itype surface type
      type (t_pbl_args) :: pbl_args
c
      pbl_args%hbaij=hbaij(i,j)
      pbl_args%ricntd=ricntd(i,j)
      pbl_args%pprec=pprec(i,j)
      pbl_args%pevap=pevap(i,j)
c**** fractional area of moist convection * fraction of downdrafts
c**** (=1/3), but only if downdrafts reach the lowest atmospheric
c**** layer. It's only needed to constrain dust emission due to
c**** downdrafts for soil type earth, but it's needed here to calculate
c**** wspdf in PBL.f for the other soil types.
      IF (ddml(i,j) == 1) THEN
        pbl_args%mcfrac=By3
      ELSE
        pbl_args%mcfrac=0.D0
      END IF
      return
      end subroutine dust_emission_prep

      subroutine save_dust_emission_vars(i,j,itype,pbl_args)
      use trdust_mod, only : hbaij,ricntd
      implicit none
      integer, intent(in) :: i,j  !@var i,j grid point
      integer, intent(in) :: itype  !@var itype surface type
      type (t_pbl_args) :: pbl_args
      hbaij(i,j)=pbl_args%hbaij
      ricntd(i,j)=pbl_args%ricntd
      return
      end subroutine save_dust_emission_vars
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      SUBROUTINE PBL_adiurn_dust(I,J,ITYPE,PTYPE,pbl_args,atm)
      use exchange_types
      USE CONSTANT, only :  rgas,grav,deltx,teeny
      use OldTracer_mod, only : dodrydep,trname
      USE TRACER_COM, only : ntm,ntm_dust
      use SOCPBL, only : npbl=>n
      USE PBLCOM
      USE DIAG_COM, only :
     *      adiurn=>adiurn_loc,ndiupt,ndiuvar,ijdd,adiurn_dust
#ifdef USE_HDIURN
     *     ,hdiurn=>hdiurn_loc
#endif
     *     ,idd_wtke,idd_wd,idd_wm,idd_wsgcm,idd_wspdf,idd_wtrsh
#ifdef TRACERS_DUST
     *     ,idd_emis,idd_emis2
     *     ,idd_ws2,idd_ustar,idd_us3,idd_stress,idd_lmon
     *     ,idd_rifl,idd_zpbl1,idd_uabl1,idd_vabl1,idd_uvabl1,idd_tabl1
     *     ,idd_qabl1,idd_zhat1,idd_e1,idd_km1,idd_ri1,idd_grav,idd_turb
#endif

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J  !@var I,J grid point
      INTEGER, INTENT(IN) :: ITYPE  !@var ITYPE surface type
      REAL*8, INTENT(IN) :: PTYPE  !@var PTYPE percent surface type
      type (t_pbl_args) :: pbl_args
      class (atmsrf_xchng_vars) :: atm

#ifdef TRACERS_MINERALS
      INTEGER,PARAMETER :: n_idxd=6
#else
#ifdef TRACERS_DUST
      INTEGER,PARAMETER :: n_idxd=16+6*npbl+4*(npbl-1)
#endif
#endif

      REAL*8 ts,ws,rhosrf
      integer nx,n
      INTEGER :: kr,ii,idxd(n_idxd)
      REAL*8 :: tmp(NDIUVAR)


      IF (adiurn_dust == 1 .and. pbl_args%moddd==0) THEN
C**** QUANTITIES ACCUMULATED HOURLY FOR DIAGDD
        DO KR=1,NDIUPT
          IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
            tmp(:) = 0.
            ii = i
            idxd=(/ idd_wsgcm,idd_wspdf,idd_wtke,idd_wd,idd_wm,idd_wtrsh
#ifdef TRACERS_DUST
     &        ,idd_emis,   idd_emis2,  idd_turb,   idd_grav,   idd_ws2
     &        ,idd_ustar,  idd_us3,    idd_stress, idd_lmon,   idd_rifl,
     &       (idd_zpbl1+ii-1,ii=1,npbl), (idd_uabl1+ii-1,ii=1,npbl),
     *       (idd_vabl1+ii-1,ii=1,npbl), (idd_uvabl1+ii-1,ii=1,npbl),
     *       (idd_tabl1+ii-1,ii=1,npbl), (idd_qabl1+ii-1,ii=1,npbl),
     *       (idd_zhat1+ii-1,ii=1,npbl-1), (idd_e1+ii-1,ii=1,npbl-1),
     *       (idd_km1+ii-1,ii=1,npbl-1), (idd_ri1+ii-1,ii=1,npbl-1)
#endif
     &      /)
            rhosrf=100.*pbl_args%psurf/(rgas*pbl_args%tsv)
            TS=pbl_args%TSV/(1.+pbl_args%QSRF*xdelt)
            ws = pbl_args%ws
            tmp(idd_wsgcm)=pbl_args%wsgcm*ptype
            tmp(idd_wspdf)=pbl_args%wspdf*ptype
            tmp(idd_wtke)=pbl_args%wsubtke*ptype
            tmp(idd_wd)=pbl_args%wsubwd*ptype
            tmp(idd_wm)=pbl_args%wsubwm*ptype
            if(itype==4) then        ! land
              tmp(idd_wtrsh)=pbl_args%wtrsh*ptype
#ifdef TRACERS_DUST
              tmp(idd_emis)=0.D0
              tmp(idd_emis2)=0.D0
              do n=1,Ntm_dust
                tmp(idd_emis)=tmp(idd_emis)
     &               +pbl_args%dust_flux(n)*ptype
                tmp(idd_emis2)=tmp(idd_emis2)
     &               +pbl_args%dust_flux2(n)*ptype
              enddo
#endif
            endif
#ifdef TRACERS_DUST
            tmp(idd_turb)=0.D0
            tmp(idd_grav)=0.D0
#ifdef TRACERS_DRYDEP
            do nx=1,pbl_args%ntx
              n = pbl_args%ntix(nx)
              if (dodrydep( n )) then
                select case(trname( n ))
                case('Clay','Silt1','Silt2','Silt3','Silt4','Silt5')
                  tmp( idd_turb ) = tmp( idd_turb ) + ptype * rhosrf
     &                 * pbl_args%trs(nx) * pbl_args%dep_vel(n)
                  tmp( idd_grav ) = tmp( idd_grav ) + ptype * rhosrf
     &                 * pbl_args%trs(nx) * pbl_args%gs_vel(n)
                end select
              end if
            end do
#endif
            tmp(idd_ws2)=ws*ws*ptype
            tmp(idd_ustar)=pbl_args%ustar*ptype
            tmp(idd_us3)=ptype*pbl_args%ustar**3
            tmp(idd_stress)=rhosrf*pbl_args%cm*(pbl_args%ws**2)*ptype
            tmp(idd_lmon)=pbl_args%lmonin*ptype
            tmp(idd_rifl)=
     &           +ptype*grav*(ts-pbl_args%tg)*pbl_args%zgs/
     &           (ws*ws*pbl_args%tg)
            tmp(idd_zpbl1:idd_zpbl1+npbl-1)=ptype*pbl_args%z(1:npbl)
            tmp(idd_uabl1:idd_uabl1+npbl-1)=
     *           ptype*atm%uabl(1:npbl,i,j)
            tmp(idd_vabl1:idd_vabl1+npbl-1)=
     *           ptype*atm%vabl(1:npbl,i,j)
            tmp(idd_uvabl1:idd_uvabl1+npbl-1)=ptype*sqrt(
     *           atm%uabl(1:npbl,i,j)*atm%uabl(1:npbl,i,j)+
     *           atm%vabl(1:npbl,i,j)*atm%vabl(1:npbl,i,j))
            tmp(idd_tabl1:idd_tabl1+npbl-1)=
     *           ptype*atm%tabl(1:npbl,i,j)
            tmp(idd_qabl1:idd_qabl1+npbl-1)=
     *           ptype*atm%qabl(1:npbl,i,j)
            tmp(idd_zhat1:idd_zhat1+npbl-2)=ptype
     *           *pbl_args%zhat(1:npbl-1)
            tmp(idd_e1:idd_e1+npbl-2)=atm%eabl(1:npbl-1,i,j)*ptype
            tmp(idd_km1:idd_km1+npbl-2)=ptype*pbl_args%km(1:npbl-1)
            tmp(idd_ri1:idd_ri1+npbl-2)=ptype*pbl_args%gh(1:npbl-1)
     *           /(pbl_args%gm(1:npbl-1)+1d-20)
#endif
            ADIURN(idxd(:),kr,pbl_args%ih)=
     &      ADIURN(idxd(:),kr,pbl_args%ih)  +tmp(idxd(:))
#ifdef USE_HDIURN
            HDIURN(idxd(:),kr,pbl_args%ihm)=
     &      HDIURN(idxd(:),kr,pbl_args%ihm) +tmp(idxd(:))
#endif
          END IF
        END DO
      END IF

      RETURN
      END SUBROUTINE PBL_adiurn_dust
#endif

      end module PBL_DRV
      
      subroutine read_pbl_tsurf_from_nmcfile
      USE DOMAIN_DECOMP_ATM, only : GRID
      use fluxes, only : atmsrf
      use pario, only : par_open,par_close,read_dist_data
      implicit none
      integer :: fid
      fid = par_open(grid,'AIC','read')
      call read_dist_data(grid,fid,'tsurf',atmsrf%tsavg)
      call par_close(grid,fid)
      atmsrf%tgvavg(:,:) = atmsrf%tsavg(:,:) ! not used for init. set anyway.
      return
      end subroutine read_pbl_tsurf_from_nmcfile

      subroutine init_pbl(inipbl,istart)
!@sum init_pbl sets up the initialization of wind,
!@+  virtual potential temperature, and specific humidity
!@+  fields in the boundary layer (between
!@+  the surface and the middle
!@+  of the first GCM layer). The initial values of these
!@+  fields are obtained by solving the static equations for these
!@+  fields using the turbulence model of Cheng et al. (2002).
!@+  These initial values are used when starting from a restart
!@+  file that does not have these data stored.
!@+  It is called by subroutine INPUT (in MODELE.f).
c -------------------------------------------------------------
c These routines include the array ipbl which indicates if the
c  computation for a particular ITYPE was done last time step.
c -------------------------------------------------------------
      USE Dictionary_mod
      USE CONSTANT, only : lhe,lhs,tf,omega2,deltx,UNDEF_VAL,rgas,grav
      USE ATM_COM, only : u,v,p,t,q
      USE ATM_COM, only : traditional_coldstart_aic
      USE GEOM, only : imaxj,sinlat2d
#ifdef TRACERS_ON
      use TRACER_COM, only: NTM
#endif
!      USE SOCPBL, only : dpdxr,dpdyr,dpdxr0,dpdyr0

      USE SOCPBL, only : npbl=>n,zgs,inits,XCDpbl,ccoeff0,skin_effect
     &     ,xdelt, maxNTM, calc_wspdf
      USE GHY_COM, only : fearth
      USE PBLCOM
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : WRITET_PARALLEL, getDomainBounds
      USE ATM_COM, only : pmid,pk,pedn,pek,pdsig
     &    ,DPDX_BY_RHO,DPDY_BY_RHO,DPDX_BY_RHO_0,DPDY_BY_RHO_0
     &    ,ua=>ualij,va=>valij
      USE SEAICE_COM, only : si_atm
      USE FLUXES, only : atmocn,atmice,atmgla,atmlnd,flice,fland
     &     ,asflx,atmsrf
#ifdef SCM
      USE SCM_COM, only : SCMopt,SCMin
#endif
      use pario, only : par_open,par_close,read_dist_data

      IMPLICIT NONE
C**** ignore ocean currents for initialisation.
      real*8, parameter :: uocean=0.,vocean=0.
!@var inipbl whether to init prog vars
      logical, intent(in) :: inipbl
!@var istart what kind of (re)start is being done
      integer, intent(in) :: istart

!@var fid unit number for roughness length input file
      integer :: fid
      integer :: ilong  !@var ilong  longitude identifier
      integer :: jlat   !@var jlat  latitude identifier
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,4) ::
     *                                                      tgvdat

      integer :: ipatch,itype4  !@var itype surface type
      integer i,j,k,lpbl !@var i,j,k loop variable
      real*8 pland,pwater,plice,psoil,poice,pocean,
     *     ztop,elhx,coriol,tgrndv,pij,ps,psk,qgrnd
     *     ,utop,vtop,qtop,ttop,zgrnd,cm,ch,cq,ustar,lmonin,tmp
      real*8 qsat
      real*8 ::  dpdxr,dpdyr,dpdxr0,dpdyr0
      real*8, dimension(npbl) :: upbl,vpbl,tpbl,qpbl
      real*8, dimension(npbl-1) :: epbl
      real*8 ug,vg
      real*8, allocatable :: buf(:,:)
      real*8 :: canopy_height, fv
      integer, save :: roughl_from_file = 0

      real*8 :: cdm

      integer :: I_1, I_0, J_1, J_0
      integer :: I_1H, I_0H, J_1H, J_0H

       character*80 :: titrrr
       real*8 rrr(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &            grid%J_STRT_HALO:grid%J_STOP_HALO)


      ! todo: perhaps use presence of LKTAB file to automatically determine this
      call sync_param( 'calc_wspdf', calc_wspdf )
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      calc_wspdf = 1
#endif

      if(calc_wspdf == 1) call init_wspdf_mod

#ifdef TRACERS_ON
       maxNTM = NTM
#endif


        titrrr = "roughness length over land"
        rrr = 0.
        tgvdat = UNDEF_VAL

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,       J_STOP=J_1)

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      if(istart==2) then ! replace with cold vs warm start logic

        if(traditional_coldstart_aic) then
          ! todo: get this via other means
          call read_pbl_tsurf_from_nmcfile
        endif
        CDM=.001d0

        DO J=J_0,J_1
        DO I=I_0,I_1
          pblht(i,j) = (.5*pdsig(1,i,j)/pmid(1,i,j))*
     &         (rgas*t(i,j,1)*pk(1,i,j)/grav)
          atmsrf%usavg(i,j) = ua(1,i,j)
          atmsrf%vsavg(i,j) = va(1,i,j)
          atmsrf%wsavg(i,j) =
     &         sqrt(atmsrf%usavg(i,j)**2 + atmsrf%vsavg(i,j)**2)
C**** SET SURFACE MOMENTUM TRANSFER TAU0
          atmsrf%TAUAVG(I,J)=1.*CDM*atmsrf%WSAVG(I,J)**2  ! air density = 1 kg/m3
C**** Initialize surface friction velocity
          do ipatch=1,size(asflx)
            asflx(ipatch)%USTAR_pbl(I,J)=atmsrf%WSAVG(I,J)*SQRT(CDM)
            asflx(ipatch)%lmonin_pbl(I,J)=100.d0
          enddo
          atmsrf%ustar_pbl(i,j)=0.1d0
          atmsrf%lmonin_pbl(i,j)=100.d0
C**** SET SURFACE SPECIFIC HUMIDITY FROM FIRST LAYER HUMIDITY
          atmsrf%QSAVG(I,J)=Q(I,J,1)
          atmsrf%QGAVG(I,J)=Q(I,J,1)
#ifdef SCM
          if( SCMopt%Qskin )then ! force skin water vapor mixing ratio
            atmsrf%qgavg(i,j) = SCMin%Qskin
          endif
#endif
        ENDDO
        ENDDO
      endif

C things to be done regardless of inipbl

      call sync_param( 'roughl_from_file', roughl_from_file )
!!#if ( ! defined ROUGHL_HACK )
      if ( roughl_from_file .ne. 0 ) then
        allocate ( buf(I_0H:I_1H, J_0H:J_1H) )
        buf = 0.
        fid = par_open(grid,'CDN','read')
        call read_dist_data(grid,fid,'al30rl',buf)
        call par_close(grid,fid)
        roughl(:,:)=30./(10.**buf(:,:))
        deallocate ( buf )
      endif
!!#endif

      call sync_param( 'XCDpbl', XCDpbl )
      call sync_param( 'skin_effect', skin_effect )
#ifdef SCM  
      if( .not. SCMopt%sfcQrad )then
c**** disable skin_effect when ignoring longwave atmospheric heating
        skin_effect=0
      endif
#endif

      do j=J_0,J_1
        do i=I_0,I_1
C**** fix roughness length for ocean ice that turned to land ice
          if (si_atm%snowi(i,j).lt.-1.and.flice(i,j).gt.0)
     &         roughl(i,j)=30./(10.**1.84d0)
          if (fland(i,j).gt.0.and.roughl(i,j) .gt. 29.d0) then
            print*,"Roughness length not defined for i,j",i,j
     *           ,roughl(i,j),fland(i,j),flice(i,j)
            print*,"Setting to .01"
            roughl(i,j)=1d-2
          end if
        end do
      end do

      call ccoeff0
      call getztop(zgs,ztop)

      if(.not.inipbl) return

      do j=J_0,J_1
      do i=I_0,I_1
        pland=fland(i,j)
        pwater=1.-pland
        plice=flice(i,j)
        psoil=fearth(i,j)
        poice=si_atm%rsi(i,j)*pwater
        pocean=pwater-poice
        if (pocean.le.0.) then
          tgvdat(i,j,1)=0.
        else
          tgvdat(i,j,1)=atmocn%gtemp(i,j)+TF
        end if
        if (poice.le.0.) then
          tgvdat(i,j,2)=0.
        else
          tgvdat(i,j,2)=atmice%gtemp(i,j)+TF
        end if
        if (plice.le.0.) then
          tgvdat(i,j,3)=0.
        else
          tgvdat(i,j,3)=atmgla%gtemp(i,j)+TF
        end if
        if (psoil.le.0.) then
          tgvdat(i,j,4)=0.
        else
          tgvdat(i,j,4)=atmlnd%gtemp(i,j)+TF
        end if
      end do
      end do

      do ipatch=1,size(asflx)
        itype4 = asflx(ipatch)%itype4
        if ((itype4.eq.1).or.(itype4.eq.4)) then
          elhx=lhe
        else
          elhx=lhs
        endif

        do j=J_0,J_1
          jlat=j
          do i=I_0,imaxj(j)
            coriol=sinlat2d(i,j)*omega2
            tgrndv=tgvdat(i,j,itype4)
            if (tgrndv.eq.0.) then
              asflx(ipatch)%ipbl(i,j)=0
              go to 200
            endif
            ilong=i
            pij=p(i,j)
            ps=pedn(1,i,j)    !pij+ptop
            psk=pek(1,i,j)    !expbyk(ps)
            qgrnd=qsat(tgrndv,elhx,ps)
#ifdef SCM
            if( SCMopt%Qskin )then ! force skin water vapor mixing ratio
              qgrnd = SCMin%Qskin
            endif
#endif
            utop = ua(1,i,j)
            vtop = va(1,i,j)
            qtop=q(i,j,1)
            ttop=t(i,j,1)*(1.+qtop*xdelt)*psk
            t1_after_aturb(i,j) = ttop/psk
            u1_after_aturb(i,j) = utop
            v1_after_aturb(i,j) = vtop

            zgrnd=.1d0 ! formal initialization
            if (itype4.gt.2) zgrnd=roughl(i,j) !         30./(10.**roughl(i,j))

            if (itype4.gt.2) rrr(i,j) = zgrnd

            dpdxr  = DPDX_BY_RHO(i,j)
            dpdyr  = DPDY_BY_RHO(i,j)
            dpdxr0 = DPDX_BY_RHO_0(i,j)
            dpdyr0 = DPDY_BY_RHO_0(i,j)
#ifdef SCM
            utop = u(i,j,1)
            vtop = v(i,j,1)
            ug = utop
            vg = vtop
#endif
            call inits(tgrndv,qgrnd,zgrnd,zgs,ztop,utop,vtop,
     2                 ttop,qtop,coriol,cm,ch,cq,ustar,lmonin,
     3                 uocean,vocean,ilong,jlat,itype4
     &                 ,dpdxr,dpdyr,dpdxr0,dpdyr0
     &                 ,upbl,vpbl,tpbl,qpbl,epbl,ug,vg)
            asflx(ipatch)%cmgs(i,j)=cm
            asflx(ipatch)%chgs(i,j)=ch
            asflx(ipatch)%cqgs(i,j)=cq

            do lpbl=1,npbl
              asflx(ipatch)%uabl(lpbl,i,j)=upbl(lpbl)
              asflx(ipatch)%vabl(lpbl,i,j)=vpbl(lpbl)
              asflx(ipatch)%tabl(lpbl,i,j)=tpbl(lpbl)
              asflx(ipatch)%qabl(lpbl,i,j)=qpbl(lpbl)
            end do

            do lpbl=1,npbl-1
              asflx(ipatch)%eabl(lpbl,i,j)=epbl(lpbl)
            end do

            asflx(ipatch)%ipbl(i,j)=1
            asflx(ipatch)%ustar_pbl(i,j)=ustar
            asflx(ipatch)%lmonin_pbl(i,j)=lmonin

 200      end do
        end do
      end do

      !write(981) titrrr,rrr
#ifndef CUBED_SPHERE
      call WRITET_PARALLEL(grid,981,"fort.981",rrr,titrrr)
#endif

      return
 1000 format (1x,//,1x,'completed initialization, itype = ',i2,//)
      end subroutine init_pbl

      subroutine loadbl
!@sum loadbl initializes boundary layer calculation each surface time step.
!@+  It checks to see if ice has
!@+  melted or frozen out of one grid box (i,j).
!@+  It is called from subroutine SURFCE (in SURFACE.f).
!@auth Greg. Hartke/Ye Cheng
c ----------------------------------------------------------------------
c             This routine checks to see if ice has
c              melted or frozen out of a grid box.
c
c For ITYPE=1 (ocean; melted ocean ice since last time step):
c  If there was no computation made for ocean at the last time step,
c  this time step may start from ocean ice result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done. Also deals with newly created lake (from land)
c
c For ITYPE=2 (ocean ice; frozen from ocean since last time step):
c  If there was no computation made for ocean ice at the last time step,
c  this time step may start from ocean result. If there was no
c  ocean nor ocean ice computation at the last time step, nothing
c  need be done.
c
c For ITYPE=3 (land ice; frozen on land since last time step):
c  If there was no computation made for land ice at the last time step,
c  this time step may start from land result. If there was no
c  land ice nor land computation at the last time step, nothing
c  need be done.
c
c For ITYPE=4 (land; melted land ice since last time step):
c  If there was no computation made for land at the last time step,
c  this time step may start from land ice result. If there was no
c  land nor land ice computation at the last time step, nothing
c  need be done. Also deal with newly created earth (from lake)
c
c In the current version of the GCM, there is no need to check the
c  land or land ice components of the grid box for ice formation and
c  melting because pland and plice are fixed. The source code to do
c  this is retained and deleted in the update deck in the event this
c  capability is added in future versions of the model.
c ----------------------------------------------------------------------
      USE EXCHANGE_TYPES
      USE MODEL_COM
      USE GEOM, only : imaxj
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE FLUXES, only : asflx,atmocns,atmices,atmlnds
#ifdef GLINT2
#define ATMGLAX atmglas_hp
#else
#define ATMGLAX atmglas
#endif
      USE FLUXES, only : ATMGLAX

      IMPLICIT NONE
      integer i,j,ip  !@var i,j,ip loop variable
      type(atmsrf_xchng_vars), pointer :: ain,aout

      integer :: J_1, J_0, I_1, I_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

! NOTE ON SUBDIVISIONS OF SURFACE TYPES:
! when initializing profiles for a newly existent surface type,
! values are taken from the _first_ instance of the preexisting
! "donor" surface type (index 1).  This follows the coding before
! the introduction of subdivisions of surface types.  If an instance
! of the surface type already exists in the gridbox, new instances
! of that type could be initialized from it.

      do j=J_0,J_1
        do i=I_0,imaxj(j)

c ******* itype=1: Ocean

          do ip=1,ubound(atmocns,1)
            aout => atmocns(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmices(1)%ipbl(i,j).eq.1) then
                ain => atmices(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              elseif (atmlnds(1)%ipbl(i,j).eq.1) then ! initialise from land
                ain => atmlnds(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

c ******* itype=2: Ocean ice

          do ip=1,ubound(atmices,1)
            aout => atmices(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmocns(1)%ipbl(i,j).eq.1) then
                ain => atmocns(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

c ******* itype=3: Land ice
! We use atmglas_hp here because PBL() is run on
! atmglas_hp, not atmglas
          do ip=1,ubound(ATMGLAX,1)
            aout => ATMGLAX(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (atmlnds(1)%ipbl(i,j).eq.1) then
                ain => atmlnds(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

c ******* itype=4: Land

          do ip=1,ubound(atmlnds,1)
            aout => atmlnds(ip)%atmsrf_xchng_vars
            if (aout%ipbl(i,j).eq.0) then
              if (ATMGLAX(1)%ipbl(i,j).eq.1) then
                ain => ATMGLAX(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              elseif (atmocns(1)%ipbl(i,j).eq.1) then
                ain => atmocns(1)%atmsrf_xchng_vars
                call setbl(ain,aout,i,j)
              endif
            endif
          enddo

C**** initialise some pbl common variables

          do ip=1,size(asflx)
            asflx(ip)%ipbl(i,j) = 0 ! - will be set to 1s when pbl is called
          enddo
        end do
      end do

      return
      end subroutine loadbl

      subroutine setbl(ain,aout,i,j)
!@sum setbl initiallise bl from another surface type for one grid box
!@+  It is called from subroutine loadbl.
!@auth Ye Cheng
      USE EXCHANGE_TYPES
      USE PBLCOM, only : npbl
      IMPLICIT NONE
      type (atmsrf_xchng_vars) :: ain,aout
      integer, INTENT(IN) :: i,j
      integer lpbl  !@var lpbl loop variable
      aout%uabl(:,i,j)=ain%uabl(:,i,j)
      aout%vabl(:,i,j)=ain%vabl(:,i,j)
      aout%tabl(:,i,j)=ain%tabl(:,i,j)
      aout%qabl(:,i,j)=ain%qabl(:,i,j)
      aout%eabl(:,i,j)=ain%eabl(:,i,j)
#ifdef TRACERS_ON
      aout%trabl(:,:,i,j)=ain%trabl(:,:,i,j)
#endif
      aout%cmgs(i,j)=ain%cmgs(i,j)
      aout%chgs(i,j)=ain%chgs(i,j)
      aout%cqgs(i,j)=ain%cqgs(i,j)
      aout%ustar_pbl(i,j)=ain%ustar_pbl(i,j)
      aout%lmonin_pbl(i,j)=ain%lmonin_pbl(i,j)
      return
      end subroutine setbl

      subroutine getztop(zgs,ztop)
!@sum  getztop computes the value of ztop which is the height in meters
!@+  of the first GCM layer from the surface.
!@+  This subroutine only needs to be called when the BL fields require
!@+  initialization.
!@+  This form for z1 = zgs + zs1 (in terms of GCM parameters) yields an
!@+  average value for zs1. The quantity theta was computed on the
!@+  assumption of zs1=200 m from the original 9-layer model (actually
!@+  was misconstrued as z1 = 200m when it should have been zs1 = 200m)
!@+  and is then applied to all vertical resolutions.
!@auth Greg. Hartke/Ye Cheng
!@var zgs The height of the surface layer.
!@var ztop The height of the top of the BL simulation domain.
!@+   Corresponds to averaged height of the middle of first model layer.

      USE CONSTANT, only : rgas,grav
      USE RESOLUTION, only : psf
      use ATM_COM, only : pednl00 ! use plbot from res file instead
      IMPLICIT NONE

      REAL*8, INTENT(IN) :: ZGS
      REAL*8, INTENT(OUT) :: ZTOP
      real*8, parameter :: theta=269.0727251d0

      ztop=zgs+0.5d0*(pednl00(1)-pednl00(2))*rgas*theta/(grav*psf)

      return
      end subroutine getztop

      subroutine get_dbl
!@sum
!@+   called from SURFACE.f
      USE FLUXES, only : atmsrf
      USE CONSTANT, only :  rgas,grav,omega,omega2,deltx,teeny
      USE ATM_COM, only : t,q,ua=>ualij,va=>valij
      USE ATM_COM, only : pmid,pk
      use SOCPBL, only : zgs
      USE PBLCOM
      USE GEOM, only : imaxj,sinlat2d
      use PBL_DRV
      use domain_decomp_atm, only : grid
      use PBL_DRV, only : dbls0,slope0,dbl_max_stable
      USE RESOLUTION, only : lm
      implicit none
      integer :: i,j,l,ldbl,ldbls
      real*8 :: ztop,coriol,dbl,ustar,lmonin,tmp,dbls
      real*8 :: zpbl,pl1,tl1,pl,tl,tbar
      real*8 :: thbar ! function

      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)

        ztop = zgs +
     &       .5d-2*RGAS*((atmsrf%temp1(I,J)*(1.+atmsrf%q1(i,j)*xdelt))*
     &       atmsrf%srfpk(i,j))*atmsrf%AM1(i,j)/atmsrf%p1(i,j)
      ! FIND THE PBL HEIGHT IN METERS (DBL) AND THE GCM LAYER IMMEDIATELY
      ! ABOVE (ldbl) AT WHICH TO COMPUTE UG AND VG.
      coriol=sinlat2d(i,j)*omega2
      ldbl=max(int(dclev(i,j)+.5d0),1)
      dbl=max(pblht(i,j),dbls0)
      ustar=atmsrf%ustar_pbl(i,j)
      lmonin=atmsrf%lmonin_pbl(i,j)
      if(lmonin.gt.0.) then
        ! ATMOSPHERE IS STABLE WITH RESPECT TO THE GROUND
        tmp=max(abs(coriol),omega)
        dbls=dbls0+slope0*(abs(lmonin*ustar/tmp))**.5d0
        dbls=max(min(dbls,dbl_max_stable),zgs)
        if(dbls.le.ztop) then
          ldbls=1
        else
          zpbl=ztop
          pl1=pmid(1,i,j)         ! pij*sig(1)+ptop
          tl1=t(i,j,1)*(1.+xdelt*q(i,j,1))*pk(1,i,j)
          do l=2,lm
            pl=pmid(l,i,j)        !pij*sig(l)+ptop
            tl=t(i,j,l)*(1.+xdelt*q(i,j,l))*pk(l,i,j) !virtual,absolute
            tbar=thbar(tl1,tl)
            zpbl=zpbl-(rgas/grav)*tbar*(pl-pl1)/(pl1+pl)*2.
            if (zpbl.ge.dbls) exit
            pl1=pl
            tl1=tl
          end do
          ldbls=l
        endif
        dbl=min(dbl,dbls)
        ldbl=min(ldbl,ldbls)
      endif

      ugeo(i,j) = ua(ldbl,i,j)
      vgeo(i,j) = va(ldbl,i,j)
      bldep(i,j) = dbl

      enddo
      enddo
      return
      end subroutine get_dbl

      SUBROUTINE CHECKPBL(SUBR)
!@sum CHECKPBL checks whether PBL data are reasonable
!@+  (to check if the data contain NaN/INF).
!@+  It is called by subroutine CHECKT, the latter is called
!@+  from PROGRAM GISS_modelE (in MODELE.f).
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE PBLCOM, only : dclev
      USE FLUXES, only : atmsrf
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

      integer :: I_1, I_0, J_1, J_0, njpol
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1,
     *               J_STRT=J_0, J_STOP=J_1)
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in boundary layer data
      CALL CHECK3B(atmsrf%wsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'wsavg')
      CALL CHECK3B(atmsrf%tsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'tsavg')
      CALL CHECK3B(atmsrf%qsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'qsavg')
      CALL CHECK3B(dclev(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'dclev')
      CALL CHECK3B(atmsrf%usavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'usavg')
      CALL CHECK3B(atmsrf%vsavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'vsavg')
      CALL CHECK3B(atmsrf%tauavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'tauavg')

      CALL CHECK3B(atmsrf%tgvavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'tgvavg')
      CALL CHECK3B(atmsrf%qgavg(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'qgavg')
      CALL CHECK3B(atmsrf%w2_l1(I_0:I_1,J_0:J_1),
     &     I_0,I_1,J_0,J_1,NJPOL,1,SUBR,'w2_l1')

      END SUBROUTINE CHECKPBL

C**** Functions and subroutines for sub-gridscale wind distribution calc

      module wspdf_mod
      use PolynomialInterpolator_mod, only: interpolator2D

!@param kim dimension 1 of lookup table for mean surface wind speed integration
!@param kjm dimension 2 of lookup table for mean surface wind speed integration
      INTEGER,PARAMETER :: kim=234,kjm=234
!@var table1 array for lookup table for calculation of mean surface wind speed
!@+          local to each grid box
      REAL*8, DIMENSION(Kim,Kjm) :: table1
!@var x11 index of table1 for GCM surface wind speed from 0 to 50 m/s
!@var x21 index of table1 for sub grid scale velocity scale (sigma)
      REAL*8 :: x11(kim),x21(kjm)
      type(interpolator2D), save :: wsgInterp

      end module wspdf_mod

      subroutine init_wspdf_mod
      use wspdf_mod
      use filemanager
      use domain_decomp_atm, only : grid,broadcast,am_i_root
      implicit none
      integer :: i,j,k,ierr,io_data
      REAL*8 :: zsum
      CHARACTER :: cierr*3,name*256

c**** read in lookup table for calculation of mean surface wind speed from PDF
      IF (am_i_root()) THEN
        CALL openunit('LKTAB1',io_data,.TRUE.,.TRUE.)
        READ (io_data,IOSTAT=ierr) table1
        name=TRIM(nameunit(io_data))
        CALL closeunit(io_data)
        IF (ierr == 0) then
          write(6,*) ' Read from file '//TRIM(name)
        ELSE
          WRITE(cierr,'(I2)') ierr
          write(6,*) ' READ ERROR ON FILE '//TRIM(name)//' rc='//cierr
        END IF
      END IF
      CALL broadcast(grid,ierr)
      IF (ierr.ne.0) CALL stop_model('init_wspdf_mod: READ ERROR',255)
      CALL broadcast(grid,table1)

c**** index of table for sub grid scale velocity (sigma) from 0.0001 to 50 m/s
      zsum=0.D0
      DO j=1,Kjm
        IF (j <= 30) THEN
          zsum=zsum+0.0001d0+FLOAT(j-1)*0.00008d0
          x21(j)=zsum
        ELSE IF (j > 30) THEN
          zsum=zsum-0.055254d0+0.005471d0*FLOAT(j)-
     &         1.938365d-4*FLOAT(j)**2.d0+
     &         3.109634d-6*FLOAT(j)**3.d0-
     &         2.126684d-8*FLOAT(j)**4.d0+
     &         5.128648d-11*FLOAT(j)**5.d0
          x21(j)=zsum
        END IF
      END DO
c**** index of table for GCM surface wind speed from 0.0001 to 50 m/s
      x11(:)=x21(:)
      wsgInterp = interpolator2D(x11,x21,table1)

      end subroutine init_wspdf_mod

      subroutine sig(tke,mdf,dbl,delt,ch,ws,tsv,wtke,wd,wm)
!@sum calculate sub grid scale velocities
!@auth Reha Cakmur/Gavin Schmidt
      use constant, only : by3,grav
      implicit none
!@var tke local turbulent kinetic energy at surface (J/kg)
!@var mdf downdraft mass flux from moist convection (m/s)
!@var dbl boundary layer height (m)
!@var delt difference between surface and ground T (ts-tg) (K)
!@var tsv virtual surface temperature (K)
!@var ch local heat exchange coefficient
!@var ws grid box mean wind speed (m/s)
      real*8, intent(in):: tke,mdf,dbl,delt,tsv,ch,ws
      real*8, intent(out):: wtke,wd,wm

C**** TKE contribution  ~ sqrt( 2/3 * tke)
      wtke=sqrt(2d0*tke*by3)

C**** dry convection/turbulence contribution ~(Qsens*g*H/rho*cp*T)^(1/3)
      if (delt.lt.0d0) then
        wd=(-delt*ch*ws*grav*dbl/tsv)**by3
      else
        wd=0.d0
      endif
C**** moist convection contribution beta/frac_conv = 200.
      wm=200.d0*mdf
      return
      end SUBROUTINE sig

      subroutine integrate_sgswind(sig,wt,wmin,wmax,ws,icase,wint)
!@sum Integrate sgswind distribution for different cases
!@auth Reha Cakmur/Gavin Schmidt
      use constant, only : by3
      implicit none
!@var sig distribution parameter (m/s)
!@var wt threshold velocity (m/s)
!@var ws resolved wind speed (m/s)
!@var wmin,wmax min and max wind speed for integral
      real*8, intent(in):: sig,wt,ws,wmin,wmax
      real*8, intent(out) :: wint
      integer, intent(in) :: icase
      integer, parameter :: nstep = 100
      integer i
      real*8 x,sgsw,bysig2

C**** integrate distribution from wmin to wmax using Simpsons' rule
C**** depending on icase, integral is done on w, w^2 or w^3
C**** Integration maybe more efficient with a log transform (putting
C**** more points near wmin), but that's up to you...
C**** Use approximate value for small sig and unresolved delta function
      if ((wmax-wmin).lt.sig*dble(nstep)) then
        bysig2=1./(sig*sig)
        wint=0
        do i=1,nstep+1
          x=wmin+(wmax-wmin)*(i-1)/dble(nstep)
          if (i.eq.1.or.i.eq.nstep+1) then
            wint=wint+sgsw(x,ws,wt,bysig2,icase)*by3
          elseif (mod(i,2).eq.0) then
            wint=wint+4d0*sgsw(x,ws,wt,bysig2,icase)*by3
          else
            wint=wint+2d0*sgsw(x,ws,wt,bysig2,icase)*by3
          end if
        end do
        wint=wint*(wmax-wmin)/dble(nstep)
      else                      ! approximate delta function
        if (ws.ge.wmin.and.ws.le.wmax) then
          wint = ws**(icase-1)*(ws-wt)
        else
          wint = 0.
        end if
      end if

      return
      end

      real*8 function sgsw(x,ws,wt,bysig2,icase)
!@sum sgsw function to be integrated for sgs wind calc
!@auth Reha Cakmur/Gavin Schmidt
      use constant, only : pi
      use SpecialFunctions_mod, only : besselI0
      implicit none
      real*8, intent(in) :: x,ws,wt,bysig2
      integer, intent(in) :: icase
      real*8 :: besf,exx

      besf=x*ws*bysig2
      if (besf.lt.200d0 ) then
        exx=exp(-0.5d0*(x*x+ws*ws)*bysig2)
        sgsw=(x)**(icase)*(x-wt)*bysig2*besselI0(besf)*exx
      else ! use bessel function expansion for large besf
        sgsw=(x)**(icase)*(x-wt)*bysig2*exp(-0.5d0*(x-ws)**2*bysig2)
     *       /sqrt(besf*2*pi)
      end if

      end function sgsw

      SUBROUTINE get_wspdf(wsubtke,wsubwd,wsubwm,mcfrac,wsgcm,wspdf)
!@sum calculates mean surface wind speed using the integral over the
!@sum probability density function of the wind speed from lookup table
!@auth Reha Cakmur/Jan Perlwitz
      use wspdf_mod, only : kim,kjm,table1,x11,x21
      use wspdf_mod, only : wsgInterp

      implicit none

c Input:
!@var wsubtke velocity scale of sub grid scale turbulence
!@var wsubwd velocity scale of sub grid scale dry convection
!@var wsubwm velocity scale of sub grid scale moist convection
!@var wsgcm GCM surface wind

      REAL(KIND=8),INTENT(IN) :: wsubtke,wsubwd,wsubwm,mcfrac,wsgcm

c Output:
!@var wspdf mean surface wind speed from integration over PDF

      REAL*8,INTENT(OUT) :: wspdf

!@var sigma standard deviation of sub grid fluctuations
      REAL*8 :: sigma,ans,dy

      REAL*8 :: work_wspdf1,work_wspdf2,wsgcm1
      CHARACTER(14) :: fname='WARNING_in_PBL'
      CHARACTER(9) :: subr='get_wspdf'
      CHARACTER(5) :: vname1='wsgcm',vname2='sigma'

      wsgcm1=wsgcm

c     This is the case when wsgcm is very small and we set it
c     equal to one of the smallest values in the table index

      IF (wsgcm1 < 0.0005) wsgcm1=0.0005D0

c     If sigma <= 0.0005:

      wspdf=wsgcm1

      CALL check_upper_limit(wsgcm1,x11(Kim),fname,subr,vname1)
c     There is no moist convection, sigma is composed of TKE and DRY
c     convective velocity scale
      IF (wsubwm == 0.) THEN
        sigma=wsubtke+wsubwd
        IF (sigma > 0.0005) THEN
          CALL check_upper_limit(sigma,x21(Kjm),fname,subr,vname2)
c     Linear Polynomial fit (Default)
          ans = wsgInterp%interpolate2dlin([wsgcm1,sigma])
          wspdf=exp(ans)
        END IF
      ELSE

c     When there is moist convection, the sigma is the combination of
c     all three subgrid scale parameters (i.e. independent or dependent)
c     Takes into account that the moist convective velocity scale acts
c     only over the area with downdrafts (mcfrac).

        work_wspdf1=0.D0
        sigma=wsubtke+wsubwd+wsubwm
        IF (sigma > 0.0005) THEN
          CALL check_upper_limit(sigma,x21(Kjm),fname,subr,vname2)
c     Linear Polynomial fit (Default)
          ans = wsgInterp%interpolate2dlin([wsgcm1,sigma])
          work_wspdf1=exp(ans)*mcfrac
        END IF

        work_wspdf2=0.D0
        sigma=wsubtke+wsubwd
        IF (sigma > 0.0005) THEN
          CALL check_upper_limit(sigma,x21(Kjm),fname,subr,vname2)
c     Linear Polynomial fit (Default)
          ans = wsgInterp%interpolate2dlin([wsgcm1,sigma])
          work_wspdf2=exp(ans)*(1.d0-mcfrac)
        END IF
        wspdf=work_wspdf1+work_wspdf2

      END IF

      RETURN
      END SUBROUTINE get_wspdf

      SUBROUTINE check_upper_limit(var,varm,cwarn,subr,vname)
!@sum checks whether variable var exceeds maximum value varm. If it does
!@+   set it to maximum value and writes warning message to file cwarn
!@auth Jan Perlwitz
!@ver 1.0

      USE filemanager,ONLY : openunit,closeunit
      USE model_com,ONLY : itime

      IMPLICIT NONE

!@var var Variable to be checked and set to maximum value if larger (inout)
!@var varm Maximum value of variable to be checked (in)
!@var cwarn File name for warning message (in)
!@var subr Name of subroutine (in)
!@var vname Name of variable (in)
      REAL(KIND=8),INTENT(IN) :: varm
      CHARACTER(*),INTENT(IN) :: cwarn,subr,vname

      REAL(KIND=8),INTENT(INOUT) :: var

      INTEGER :: iu

      IF (var > varm) THEN
        CALL openunit(cwarn,iu)
        DO
          READ(iu,*,END=1)
          CYCLE
 1        EXIT
        END DO
        WRITE(iu,*) 'Warning in ',subr,':'
        WRITE(iu,*) ' Variable ',vname,' exceeds allowed maximum value.'
        WRITE(iu,*) ' ',vname,' set to maximum value.'
        WRITE(iu,*) ' itime, max of ',vname,', ',vname,':',itime,varm
     &       ,var
        CALL closeunit(iu)
        var=varm        
      END IF

      RETURN
      END SUBROUTINE check_upper_limit

