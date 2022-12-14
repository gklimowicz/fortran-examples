#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_ON
#undef TRACERS_WATER
#endif

!#define ROUGHL_HACK

c******************   TRACERS             ******************************
#ifdef TRACERS_ON
      module ghy_tracers

#ifdef TRACERS_WATER
      use ghy_h, only : ghy_tr_str
#endif
      use sle001, only :
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use sle001,ONLY : aevap
#endif
      use OldTracer_mod, only: itime_tr0, trname, needtrs
      use tracer_com, only: NTM
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use tracer_com, only: Ntm_dust, n_soildust
#endif
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only:dodrydep
#endif
#ifdef TRACERS_WATER
      use OldTracer_mod, only: nWATER, nGAS, nPARt, tr_wd_TYPE
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      USE trdust_mod,ONLY : imDust
#endif
      use trdiag_com, only : taijn=>taijn_loc,
     *   taijs=>taijs_loc,ijts_isrc,jls_isrc,
     *   itcon_surf
#ifdef TRACERS_WATER
     *     ,tij_soil,tij_snow,tre_acc
#endif
#ifdef TRACERS_DRYDEP
     *     ,tij_drydep,tij_gsdep,itcon_dd
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
     &     ,ijts_spec,jls_spec
#endif
#ifdef TRACERS_ON
     &     ,trcsurf,trcSurfByVol
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE trdust_mod,ONLY : nDustEmij,nDustEm2ij,nDustEmjl,nDustEm2jl
     &     ,nDustEv1ij,nDustEv2ij,nDustWthij
     &     ,nDustEv1jl,nDustEv2jl,nDustWthjl
#endif
      use fluxes, only : trflux1
#ifdef TRACERS_WATER
     *     ,trprec
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
     &     ,prec,pprec,pevap,dust_flux_glob,dust_flux2_glob
#endif

#ifdef TRACERS_WATER
      use ghy_com, only : tr_w_ij,tr_wsn_ij !, ntixw
#endif

ccc extra stuff which was present in "earth" by default
#ifdef TRACERS_WATER
      use ghy_com, only : ngm,nlsn
      use constant, only : rhow
#endif
      use geom, only : byaxyp

      implicit none
      private

      public ghy_tracers_set_step
      public ghy_tracers_finish_step
      public ghy_tracers_set_cell
      public ghy_tracers_set_cell_stage2
      public ghy_tracers_save_cell
      public initGhyTracers

      integer, allocatable :: ntix(:)
      integer ntx
      integer, parameter :: itype=4

      contains

      subroutine initGhyTracers
      allocate(ntix(ntm))
      end subroutine initGhyTracers

      subroutine ghy_tracers_set_step(
#ifdef TRACERS_WATER
     &     ghy_tr
#endif
     &     )
!@sum tracers stuff to be called at the beginning of ghy time step
!@+   i.e. before the i,j loop
      use model_com, only : itime
      implicit none
#ifdef TRACERS_WATER
      type (ghy_tr_str) :: ghy_tr
      integer :: ntg, ntixw(ntm)
      logical :: is_water(ntm)
#endif
      integer n

      ntx=0
#ifdef TRACERS_WATER
      ntg=0
#endif
      do n=1,ntm
        if (itime_tr0(n).le.itime .and. needtrs(n)) then
          ntx=ntx+1
          ntix(ntx) = n
#ifdef TRACERS_WATER
          if ( tr_wd_TYPE(n)==nWATER .or. tr_wd_TYPE(n)==nPART ) then
            ntg = ntg + 1
            ntixw(ntg) = n
            if ( tr_wd_TYPE(n)==nWATER ) then
              is_water(ntg) = .true.
            else
              is_water(ntg) = .false.
            endif
          endif
#endif
        end if
      end do

#ifdef TRACERS_WATER
      ghy_tr%ntg = ntg
      allocate( ghy_tr%trpr(ntg), ghy_tr%trdd(ntg), ghy_tr%trirrig(ntg),
     &     ghy_tr%tr_surf(ntg),
     &     ghy_tr%tr_w(ntg,0:ngm,2), ghy_tr%tr_wsn(ntg,nlsn,2) )
      allocate( ghy_tr%atr_evap(ntg),ghy_tr%atr_rnff(ntg),
     &     ghy_tr%atr_g(ntg), ghy_tr%ntixw(ntg), ghy_tr%is_water(ntg) )
      ghy_tr%ntixw(1:ntg) = ntixw(1:ntg)
      ghy_tr%is_water(1:ntg) = is_water(1:ntg)
#endif

      end subroutine ghy_tracers_set_step


      subroutine ghy_tracers_finish_step(
#ifdef TRACERS_WATER
     &     ghy_tr
#endif
     &     )
!@sum tracers stuff to be called at the end of ghy time step
!@+   i.e. after the i,j loop
      implicit none
#ifdef TRACERS_WATER
      type (ghy_tr_str) :: ghy_tr

      deallocate( ghy_tr%trpr, ghy_tr%trdd, ghy_tr%trirrig,
     &     ghy_tr%tr_surf,
     &     ghy_tr%tr_w, ghy_tr%tr_wsn )
      deallocate( ghy_tr%atr_evap,ghy_tr%atr_rnff,
     &     ghy_tr%atr_g, ghy_tr%ntixw, ghy_tr%is_water )
      ghy_tr%ntg = 0
#endif

      end subroutine ghy_tracers_finish_step

      subroutine ghy_tracers_set_cell(i,j,ptype,qg,evap_max,pbl_args
     &  ,qm1
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     &     )
!@sum tracers code to be called before the i,j cell is processed
      use model_com, only : dtsrc
      use pbl_drv, only : t_pbl_args

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE model_com,ONLY : modelEclock
      USE geom,ONLY : axyp
      USE ghy_com,ONLY : wearth,aiearth,wfcs
      use trdust_mod,only : nDustBins,d_dust,ers_data
     &     ,dustSourceFunction,frclay,frsilt,dryhr,vtrsh
     &     ,mineralFractions
#endif
#ifdef TRACERS_WATER
      use fluxes, only : atmlnd
#ifdef IRRIGATION_ON
      use fluxes, only : irrig_tracer_act
#endif
#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: qg,evap_max,ptype
      type (t_pbl_args), intent(inout) :: pbl_args
      real*8, intent(in) :: qm1
      integer n,nx
      integer :: month, dayOfYear

#ifdef TRACERS_WATER
      type (ghy_tr_str) :: ghy_tr
#endif

c**** pass indices of tracer arrays to PBL
      pbl_args%ntix(1:ntm) = ntix(1:ntm)
      pbl_args%ntx = ntx

ccc tracers variables
#ifdef TRACERS_WATER
      do nx=1,ghy_tr%ntg
        n = ghy_tr%ntixw(nx)
        ! prognostic vars
        ghy_tr%tr_w(nx,0,1) = 0.d0
        ghy_tr%tr_w(nx,1:ngm,1) = tr_w_ij(n,1:ngm,1,i,j)
        ghy_tr%tr_w(nx,0:ngm,2) = tr_w_ij(n,0:ngm,2,i,j)
        ghy_tr%tr_wsn(nx,1:nlsn,1:2) = tr_wsn_ij(n,1:nlsn, 1:2, i, j)
        ! flux in
        ghy_tr%trpr(nx) = trprec(n,i,j)/dtsrc ! kg/m^2 s (in precip)
#ifndef GHY_DRYDEP_FIX
#ifdef TRACERS_DRYDEP
        ghy_tr%trdd(nx) = atmlnd%trdrydep(n,i,j)/dtsrc   ! kg/m^2 s (dry dep.)
#else
        ghy_tr%trdd(nx) = 0
#endif
#endif
        ! for non-zero irrigation the following line should be set
        ! to tracer flux due to irrigation
#ifdef IRRIGATION_ON
        ghy_tr%trirrig(nx) = irrig_tracer_act(n,i,j)*rhow/ ptype ! tracers in irrigation kg/m^2 s
#ifdef TRACERS_IRRIGATION_WATER_ONLY
        if ( tr_wd_TYPE(n) .ne. nWATER) then
          ghy_tr%trirrig(nx) = 0.d0
        endif
#endif
#else
        ghy_tr%trirrig(nx) = 0.
#endif
        ! concentration of tracers in atm. water at the surface
        if (qm1.gt.0 .and. tr_wd_TYPE(n)==nWATER) then
          ghy_tr%tr_surf(nx) = atmlnd%trm1(n,i,j)*rhow/qm1 ! kg/m^3
        else
          ghy_tr%tr_surf(nx) = 0.
        end if
      enddo
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      ! todo: move (some of) this to subroutine dust_emission_prep
      call modelEclock%get(month=month, dayOfYear=dayOfYear)
      pbl_args%snow=atmlnd%snowe(i,j)
      pbl_args%wearth=wearth(i,j)
      pbl_args%aiearth=aiearth(i,j)
      pbl_args%wfcs=wfcs(i,j)
      pbl_args%ers_data=ers_data(i,j,month)
      pbl_args%dustSourceFunction = dustSourceFunction(i,j)
      pbl_args%frclay=frclay(i,j)
      pbl_args%frsilt=frsilt(i,j)
      pbl_args%vtrsh=vtrsh(i,j)
      pbl_args%dryhr=dryhr(i,j)
c**** prescribed dust emission
      pbl_args%d_dust( 1:nDustBins ) = d_dust( i, j, 1:nDustBins,
     &     dayOfYear )/ SECONDS_PER_DAY / axyp( i, j ) / ptype
c**** mineral fractions of emitted dust aerosols
      pbl_args%mineralFractions(:)=mineralFractions(i,j,:)
#endif

      end subroutine ghy_tracers_set_cell


      subroutine ghy_tracers_set_cell_stage2(i,j,pbl_args
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     &     )
!@sum tracers code to be called (after pbl) before the i,j cell is processed
      use pbl_drv, only : t_pbl_args
      USE MODEL_COM, only : qcheck
      USE TRACER_COM, only : trm
#ifdef TRACERS_WATER
      use fluxes, only : atmlnd
#endif
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: dodrydep
#endif
      implicit none
      integer, intent(in) :: i,j
      type (t_pbl_args), intent(inout) :: pbl_args
      integer n,nx
      real*8 :: tdryd, td1
#ifdef TRACERS_WATER
      type (ghy_tr_str) :: ghy_tr
#endif

#ifdef TRACERS_WATER
      ghy_tr%trdd(:) = 0.d0
#ifdef TRACERS_DRYDEP

      do nx=1,ghy_tr%ntg
        n = ghy_tr%ntixw(nx)
        if(.not.dodrydep(n)) cycle
        tdryd = atmlnd%drydflx(n,i,j) ! kg/m2
       ! tdd = tdryd             ! kg/m2
        td1 = (atmlnd%trsrfflx(n,i,j)
     &       +atmlnd%trflux_prescr(n,i,j)
     &       )*pbl_args%dtsurf           ! kg/m2
        if(trm(i,j,1,n)*byaxyp(i,j)+(td1+tdryd).lt.0.and.tdryd.lt.0)then
          if (qcheck) write(99,*) "limiting tdryd surface",i,j,n,tdryd
     *         ,trm(i,j,1,n),td1,pbl_args%trs(nx),pbl_args%trtop(nx)
          tdryd= -max(trm(i,j,1,n)*byaxyp(i,j)+td1,0d0)
          !tdryd=tdd
        end if
        ghy_tr%trdd(nx) = -tdryd/pbl_args%dtsurf   ! kg/m^2 s (dry dep.)
      enddo

#endif
#endif
      end subroutine ghy_tracers_set_cell_stage2


      subroutine ghy_tracers_save_cell(i,j,ptype,dtsurf,rhosrf,pbl_args
     &     ,excess_C_delta
#ifdef INTERACTIVE_WETLANDS_CH4
     &     ,ra_temp,ra_sat,ra_gwet
#endif
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     &     )
!@sum tracers code to be called after the i,j cell is processed
      use model_com, only : itime,qcheck
      use fluxes, only : atmlnd,nisurf
      use pbl_drv, only : t_pbl_args
      use ghy_com, only : tearth
      use somtq_com, only : mz
#ifdef TRACERS_COSMO
      USE COSMO_SOURCES, only : BE7D_acc
#endif
      USE TRACER_COM, only: ntm
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only: n_ANUM, n_AH2O, xk, nbins
      USE TOMAS_EMIS 
#endif
 !     use socpbl, only : dtsurf
      use geom, only : axyp
      use sle001, only : nsn,fb,fv
#ifdef TRACERS_GASEXCH_land_CO2
     &     ,agpp,arauto,asoilresp
#endif
#ifdef BIOGENIC_EMISSIONS
      use trdiag_com,ONLY : ijs_isoprene
#endif
#ifdef PS_BVOC
      use trdiag_com,ONLY : ijs_isoprene
      use sle001, only : aipp
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      use constant, only : tf
      use tracer_sources, only : n__temp,n__sat,n__gwet
#endif
cddd#ifdef TRACERS_GASEXCH_land_CO2
cddd      USE FLUXES, only : TRGASEX
cddd#endif
      implicit none
      integer, intent(in) :: i,j
      real*8, intent(in) :: ptype,dtsurf,rhosrf
      type (t_pbl_args), intent(in) :: pbl_args
      real*8 :: byNIsurf
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      integer :: n1
#endif
      integer n,nx
#ifdef TRACERS_DRYDEP
      real*8 tdryd,tdd,td1,rtsdt,rts,depvel,gsvel
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
      real*8 trc_flux
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      real*8, intent(IN) :: ra_temp,ra_sat,ra_gwet
#endif
#ifdef TRACERS_WATER
      type (ghy_tr_str) :: ghy_tr
#endif
      real*8 :: excess_C_delta
#ifdef TRACERS_TOMAS
      INTEGER ss_bin,num_bin,du_bin,num_bin2
      real*8 ss_num(nbins),dust_num(nbins),tot_dust,tot_seasalt
#endif
ccc tracers
      byNIsurf=1.d0/real(NIsurf)

#ifdef TRACERS_GASEXCH_land_CO2
      !!! hack - assume 1 tracer
      n = 1
c     redistribute excess_C in the atmosphere (assume the time scale
c     of about a year)
!!      delta_C = excess_C(i,j)/SECONDS_PER_YEAR*dtsurf
!!      excess_C(i,j) = excess_C(i,j) - delta_C
cddd      TRGASEX(n,4,I,J) =
cddd     &     (arauto+asoilresp-agpp)/dtsurf
      atmlnd%trsrfflx(n,i,j)=atmlnd%trsrfflx(n,i,j)+
     &     (arauto+asoilresp-agpp+excess_C_delta)/dtsurf
     &     *44.d0/12.d0
      taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n))
     &     + (arauto+asoilresp-agpp+excess_C_delta) * axyp(i,j)*ptype
     &     *44.d0/12.d0
#endif

#ifdef TRACERS_WATER
      do nx=1,ghy_tr%ntg
        n = ghy_tr%ntixw(nx)
        tr_w_ij(n,1:ngm,1,i,j) = ghy_tr%tr_w(nx,1:ngm,1)
        tr_w_ij(n,0:ngm,2,i,j) = ghy_tr%tr_w(nx,0:ngm,2)
        tr_wsn_ij(n,1:nlsn, 1:2, i, j) = ghy_tr%tr_wsn(nx,1:nlsn,1:2)
      enddo
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
c     saves precipitation for dust emission calculation at next time step
      pprec(i,j)=prec(i,j)
c     saves evaporation for dust emission calculation at next time step
      pevap(i,j)=aevap
#endif
#ifdef TRACERS_WATER
ccc accumulate tracer evaporation and runoff
      do nx=1,ghy_tr%ntg
        n=ghy_tr%ntixw(nx)
        atmlnd%trevapor(n,i,j) = atmlnd%trevapor(n,i,j)
     &       + ghy_tr%atr_evap(nx)
        atmlnd%truno(n,i,j) = atmlnd%truno(n,i,j) +ghy_tr%atr_rnff(nx)
        atmlnd%gtracer(n,i,j) = ghy_tr%atr_g(nx)
        atmlnd%trsrfflx(n,i,j)=atmlnd%trsrfflx(n,i,j)+
     &       ghy_tr%atr_evap(nx)/dtsurf
      enddo
#endif
#ifdef TRACERS_TOMAS
      ss_bin=0
      num_bin=0
      du_bin=0
      num_bin2=0
#endif
      DO nx=1,ntx
        n=ntix(nx)

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP)  ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
C**** technicallly some of these are ocean emissions, but if
C**** fixed datasets are used, it can happen over land as well.

        select case (trname(n))
        case ('DMS')
          trc_flux=pbl_args%DMS_flux
#if (defined TRACERS_AEROSOLS_SEASALT) || (defined TRACERS_AMP)
        case ('seasalt1', 'M_SSA_SS')
          trc_flux=pbl_args%ss1_flux
        case ('seasalt2', 'M_SSC_SS')
          trc_flux=pbl_args%ss2_flux
        case ('M_SSS_SS')
          trc_flux=(pbl_args%ss1_flux+pbl_args%ss2_flux)
#endif  /* TRACERS_AEROSOLS_SEASALT || TRACERS_AMP */
#ifdef TRACERS_AEROSOLS_OCEAN
        case ('OCocean')
          trc_flux=pbl_args%OCocean_flux
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_AMP
        case ('M_DD1_DU')
          trc_flux=sum(pbl_args%dust_flux(1:2))
        case ('M_DD2_DU')
          trc_flux=sum(pbl_args%dust_flux(3:4))
        case ('M_DDD_DU')
          trc_flux=sum(pbl_args%dust_flux(1:4))
#endif
!TOMAS - I don't think that seasalt emission needs here.. 
#ifdef TRACERS_TOMAS

          CASE ('ADUST_01','ADUST_02','ADUST_03','ADUST_04'
     &          ,'ADUST_05','ADUST_06','ADUST_07','ADUST_08'
     &          ,'ADUST_09','ADUST_10','ADUST_11','ADUST_12'
#ifdef TOMAS_12_3NM 
     *    ,'ADUST_13','ADUST_14','ADUST_15'
#endif
     &         )

          du_bin=du_bin+1

          trc_flux=pbl_args%dust_flux(1)*scalesizeclay(du_bin)
     &         +sum(pbl_args%dust_flux(2:4))*scalesizesilt(du_bin)
          
          dust_num(du_bin)=trc_flux/sqrt(xk(du_bin)*xk(du_bin+1))
       
          case ('ANUM__01','ANUM__02','ANUM__03','ANUM__04',
     &         'ANUM__05','ANUM__06','ANUM__07','ANUM__08',
     &         'ANUM__09','ANUM__10','ANUM__11','ANUM__12'
#ifdef TOMAS_12_3NM 
     *    ,'ANUM__13','ANUM__14','ANUM__15'
#endif
     &         )
           num_bin2=num_bin2+1
           trc_flux=dust_num(num_bin2)
#endif
        case default
          trc_flux=0
        end select

        atmlnd%trsrfflx(n,i,j)=atmlnd%trsrfflx(n,i,j)+trc_flux
#ifndef TRACERS_TOMAS
        if (ijts_isrc(1,n)>0) then
           taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &       trc_flux*axyp(i,j)*ptype*dtsurf
        end if
#else
        if(n.lt.n_ANUM(1).or.n.ge.n_AH2O(1))THEN !for dust and other?
           if (ijts_isrc(1,n)>0) then
              taijs(i,j,ijts_isrc(1,n))=taijs(i,j,ijts_isrc(1,n)) +
     &             trc_flux*axyp(i,j)*ptype*dtsurf
           end if
        elseif(n.ge.n_ANUM(1).and. n.lt.n_AH2O(1))THEN
!ijts_isrc(2,n) for number: DUST number emission
           if (ijts_isrc(2,n)>0) then
              taijs(i,j,ijts_isrc(2,n))=taijs(i,j,ijts_isrc(2,n)) +
     &             trc_flux*axyp(i,j)*ptype*dtsurf
           end if
        endif
#endif

#ifdef TRACERS_AMP
        if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)
#else
#ifdef TRACERS_TOMAS
        if(n.lt.n_ANUM(1).or.n.ge.n_AH2O(1))THEN !for dust and other?
           if(jls_isrc(1,n)>0)  call inc_tajls(i,j,1,jls_isrc(1,n),
     *          trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols?
           if (itcon_surf(1,n).gt.0) call inc_diagtcb(i,j,
     *       trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(1,n),n)

        elseif(n.ge.n_ANUM(1).and. n.lt.n_AH2O(1))THEN
!jls_isrc(2,n) for number: DUST number emission
           if(jls_isrc(2,n)>0)  call inc_tajls(i,j,1,jls_isrc(2,n),
     *          trc_flux*axyp(i,j)*ptype*dtsurf) ! why not for all aerosols?
           if (itcon_surf(5,n).gt.0) call inc_diagtcb(i,j,
     *          trc_flux*axyp(i,j)*ptype*dtsurf,itcon_surf(5,n),n)
        endif

!jls_isrc is only for DU and SS.  
!Unlike, itcon_surf, this does not need a different value for emission type.
#endif
#ifndef TRACERS_TOMAS
        if(jls_isrc(1,n)>0)  call inc_tajls(i,j,1,jls_isrc(1,n),
     *       trc_flux*axyp(i,j)*ptype*dtsurf)   ! why not for all aerosols?
#endif
#endif

#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
ccc dust emission from earth
        select case (trname(n))
        case ('Clay','Silt1','Silt2','Silt3','Silt4','Silt5','ClayIlli'
     &         ,'ClayKaol','ClaySmec','ClayCalc','ClayQuar','ClayFeld'
     &         ,'ClayHema','ClayGyps','ClayIlHe','ClayKaHe','ClaySmHe'
     &         ,'ClayCaHe','ClayQuHe','ClayFeHe','ClayGyHe','Sil1Quar'
     &         ,'Sil1Feld','Sil1Calc','Sil1Hema','Sil1Gyps','Sil1Illi'
     &         ,'Sil1Kaol','Sil1Smec','Sil1QuHe','Sil1FeHe','Sil1CaHe'
     &         ,'Sil1GyHe','Sil1IlHe','Sil1KaHe','Sil1SmHe','Sil2Quar'
     &         ,'Sil2Feld','Sil2Calc','Sil2Hema','Sil2Gyps','Sil2Illi'
     &         ,'Sil2Kaol','Sil2Smec','Sil2QuHe','Sil2FeHe','Sil2CaHe'
     &         ,'Sil2GyHe','Sil2IlHe','Sil2KaHe','Sil2SmHe','Sil3Quar'
     &         ,'Sil3Feld','Sil3Calc','Sil3Hema','Sil3Gyps','Sil3Illi'
     &         ,'Sil3Kaol','Sil3Smec','Sil3QuHe','Sil3FeHe','Sil3CaHe'
     &         ,'Sil3GyHe','Sil3IlHe','Sil3KaHe','Sil3SmHe','Sil4Quar'
     &         ,'Sil4Feld','Sil4Calc','Sil4Hema','Sil4Gyps','Sil4Illi'
     &         ,'Sil4Kaol','Sil4Smec','Sil4QuHe','Sil4FeHe','Sil4CaHe'
     &         ,'Sil4GyHe','Sil4IlHe','Sil4KaHe','Sil4SmHe','Sil5Quar'
     &         ,'Sil5Feld','Sil5Calc','Sil5Hema','Sil5Gyps','Sil5Illi'
     &         ,'Sil5Kaol','Sil5Smec','Sil5QuHe','Sil5FeHe','Sil5CaHe'
     &         ,'Sil5GyHe','Sil5IlHe','Sil5KaHe','Sil5SmHe')
          n1=n-n_soildust+1
          atmlnd%trsrfflx(n,i,j)=atmlnd%trsrfflx(n,i,j)+
     &         pbl_args%dust_flux(n1)
          taijs(i,j,ijts_isrc(nDustEmij,n))
     &         =taijs(i,j,ijts_isrc(nDustEmij,n))
     &         +pbl_args%dust_flux(n1)
     &         *axyp(i,j)*ptype*dtsurf
          if (jls_isrc(nDustEmjl,n)>0) call inc_tajls(i,j,1,jls_isrc(
     &       nDustEmjl,n),pbl_args%dust_flux(n1)*axyp(i,j)*ptype*dtsurf)
          IF ( imDust == 0 .or. imDust >= 3 ) THEN
            taijs(i,j,ijts_isrc(nDustEm2ij,n))
     &           =taijs(i,j,ijts_isrc(nDustEm2ij,n))
     &           +pbl_args%dust_flux2(n1)
     &           *axyp(i,j)*ptype*dtsurf
          if (jls_isrc(nDustEm2jl,n) > 0)
     &           call inc_tajls(i,j,1,jls_isrc(nDustEm2jl,n),
     &           pbl_args%dust_flux2(n1)*axyp(i,j)*ptype*dtsurf)
          END IF

        end select
#endif

#ifdef BIOGENIC_EMISSIONS
        select case (trname(n))
        case ('Isoprene')
          atmlnd%trsrfflx(n,i,j)=atmlnd%trsrfflx(n,i,j)+
     &         pbl_args%emisop
          taijs(i,j,ijs_isoprene)=taijs(i,j,ijs_isoprene)+
     &    pbl_args%emisop*axyp(i,j)*ptype*dtsurf
        end select
#endif
#ifdef PS_BVOC
        select case (trname(n))
        case ('Isoprene')
C Flux in kg/m2/s - put back into /s
          atmlnd%trsrfflx(n,i,j)=atmlnd%trsrfflx(n,i,j)+aipp/dtsurf
c          taijs(i,j,ijs_isoprene)=taijs(i,j,ijs_isoprene)+
c     &    aipp*axyp(i,j)*ptype*dtsurf
          taijs(i,j,ijs_isoprene)=taijs(i,j,ijs_isoprene)+
     &    aipp*axyp(i,j)*ptype

        end select
#endif
      end do

ccc not sure about the code below. hopefully that''s what is meant above
#ifdef TRACERS_WATER
      do nx=1,ghy_tr%ntg
        n=ghy_tr%ntixw(nx)
        taijn(i,j,tij_soil,n)=taijn(i,j,tij_soil,n) + (
     &       fb*(sum( ghy_tr%tr_w(nx,1:ngm,1) ))+
     &       fv*(sum( ghy_tr%tr_w(nx,0:ngm,2) ))
     *       )
        taijn(i,j,tij_snow,n)=taijn(i,j,tij_snow,n) + (
     &       fb*(sum( ghy_tr%tr_wsn(nx,1:nsn(1),1) ))+
     &       fv*(sum( ghy_tr%tr_wsn(nx,1:nsn(2),2) ))
     *       )
      TRE_acc(n,i,j)=TRE_acc(n,i,j)+atmlnd%trevapor(n,i,j)*ptype
      enddo
#endif
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
c ..........
c Accumulates dust events. One diagnostic field for all dust tracers
c ..........
      taijs(i,j,ijts_spec(nDustEv1ij))=taijs(i,j,ijts_spec(nDustEv1ij))
     &     +pbl_args%dust_event1
      call inc_tajls2(i,j,1,jls_spec(nDustEv1jl),pbl_args%dust_event1)
      taijs(i,j,ijts_spec(nDustEv2ij))=taijs(i,j,ijts_spec(nDustEv2ij))
     &     +pbl_args%dust_event2
      call inc_tajls2(i,j,1,jls_spec(nDustEv2jl),pbl_args%dust_event2)
      taijs(i,j,ijts_spec(nDustWthij))=taijs(i,j,ijts_spec(nDustWthij))
     &     +pbl_args%wtrsh*ptype
      call inc_tajls2(i,j,1,jls_spec(nDustWthjl),pbl_args%wtrsh*ptype)
#endif

c     ..........
c     save global variables for subdaily diagnostics
c     ..........

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      DO n=1,Ntm_dust
        dust_flux_glob( i, j, n ) = dust_flux_glob( i, j, n ) +
     &       pbl_args%dust_flux( n ) * ptype / nisurf
        dust_flux2_glob( i, j, n ) = dust_flux2_glob( i, j, n ) +
     &       pbl_args%dust_flux2( n ) * ptype / nisurf
      END DO
#endif

#ifdef INTERACTIVE_WETLANDS_CH4
C**** update running-average of ground temperature:
      call running_average(ra_temp,I,J,dble(nisurf),n__temp)
      call running_average(ra_sat,I,J,dble(nisurf),n__sat)
      call running_average(ra_gwet,I,J,dble(nisurf),n__gwet)
#endif
      end subroutine ghy_tracers_save_cell

      end module ghy_tracers
#endif
c******************   END   TRACERS             ************************


c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************


      module soil_drv
!@sum soil_drv contains variables and routines for the
!@+   TerraE Global Land Model (code for TerraE is located
!@+   in 'model/giss_LSM'
!@auth I. Alienov/F. Abramopolous
      use resolution, only : im,jm
      use socpbl, only : npbl=>n

#ifdef IRRIGATION_ON
      use irrigmod, only : init_irrigmod
#endif  /* IRRIGATION_ON */

      implicit none
      private
      save

      public daily_earth
      public ground_e
      public init_LSM
      public earth
      !public conserv_wtg
      !public conserv_htg
      !public conserv_wtg_1
      public checke
      public snow_cover

!@dbparam snow_cover_coef coefficient for topography variance in
!@+       snow cover parameterisation for albedo
      real*8 :: snow_cover_coef = .15d0
!@dbparam land_CO2_bc_flag type of CO2 BC to be used by Land Surface
!@+   0 - fixed, 1 - transient (from radiation), 2 - interactive
      integer :: land_CO2_bc_flag = 1
!@dbparam land_CO2_bc CO2 concentration (ppm) to be used by Land Surfacw
!@+   in "fixed" case (for land_CO2_bc_flag = 0)
      real*8 :: land_CO2_bc = 280.d0

      integer variable_lk
!@dbparam sbgc_spinup turn on the algorithm for fast soil carbon spinup. 
!@+   You really should know what you are doing to use it! Ask first.
      logical :: sbgc_spinup = .false.

      contains

      subroutine earth (ns,moddsf,moddd)
!@sum EARTH calculates surface fluxes of sensible heat,
!@+   evaporation, thermal radiation, and momentum drag over earth.
!@auth I. Alienov/F. Abramopolous
c****
      use constant, only : grav,rgas,lhe,lhs
     *     ,sha,tf,rhow,deltx,gasc,stbo
      use resolution, only : im,jm
      use model_com, only : modelEclock
      use model_com, only : dtsrc,nday,itime
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin,nstepSCM
#endif
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, AM_I_ROOT
      use geom, only : imaxj,lat2d
      use rad_com, only :
     &      FSRDIR, SRVISSURF, CO2ppm
      !use surf_albedo, only: albvnh   ! added 5/23/03 from RADIATION.f
      !albvnh(9,6,2)=albvnh(sand+8veg,6bands,2hemi) - only need 1st band
      use sle001, only : advnc,evap_limits,
!!!     &    pr,htpr,prs,htprs,
     &     w,snowd,tp,fice,
     &    ashg,alhg, !fv,fb,
     &    aevap,abetad,
     &    aruns,arunu,aeruns,aerunu,
     &    tbcs,tsns,ae0,ijdebug,thets
!     &     ,ghy_counter=>counter
!!!     &    qm1,qs,
!!!     &    pres,rho,ts,ch,srht,trht
!!!     &   ,vs,vs0,tprime,qprime
      use fluxes, only : atmlnd,prec,eprec
     *     ,precss,nisurf, asflx
      use ghy_com, only : snowbv, fearth,
     &     fr_snow_ij,
     *     tearth,tsns_ij,wearth,aiearth,
     &     evap_max_ij, fr_sat_ij, qg_ij, top_dev_ij,
     &     soil_surf_moist, snowd_ij=>snowd
#ifdef CACHED_SUBDD
      use ghy_com, only : gsaveL
      use constant, only : undef
#endif
      use pbl_drv, only : pbl, t_pbl_args, xdelt
      use pbl_drv, only : alloc_pbl_args, dealloc_pbl_args

      use snow_drvm, only : snow_cover_same_as_rad
      use snow_model, only : i_earth,j_earth


#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups
     &     ,inc_subdd,find_groups
#endif

#ifdef TRACERS_ON
      use ghy_tracers, only : ghy_tracers_set_step,ghy_tracers_set_cell,
     &     ghy_tracers_save_cell, ghy_tracers_finish_step,
     &     ghy_tracers_set_cell_stage2
#endif
#ifdef WATER_PROPORTIONAL
      use tracer_com, only : NTM,trm
      use geom, only : axyp
#endif
      use ent_com, only : entcells, excess_C
      !use ent_mod, only : ent_prescribe_vegupdate
      use ent_drv, only : update_vegetation_data
      !use ent_mod, only : ent_cell_print
#ifdef IRRIGATION_ON
      use fluxes, only : irrig_water_act, irrig_energy_act
#ifdef TRACERS_WATER
     *     ,irrig_tracer_act
#endif
#endif
      use ghy_com, only : ngm,imt,nlsn,LS_NFRAC,dz_ij,sl_ij,q_ij,qk_ij
     *     ,top_index_ij,top_dev_ij
     &     ,w_ij,ht_ij,snowbv,nsn_ij,dzsn_ij,wsn_ij
     &     ,hsn_ij,fr_snow_ij,shc_soil_texture
     &     ,Qf_ij
#ifdef TRACERS_ON
!!! need for Ca hack
      use geom, only : byaxyp
      use tracer_com, only : trm
#endif
#ifdef TRACERS_WATER
      use GHY_h, only : ghy_tr_str
#endif
#ifdef TRACERS_GASEXCH_land_CO2
      use tracer_com, only : n_CO2n
#endif
      use TimeConstants_mod, only: SECONDS_PER_YEAR
      implicit none

      integer, intent(in) :: ns,moddsf,moddd
      integer i,j,k,itype,ibv
      real*8 shdt,evhdt,rcdmws,rcdhws
     *     ,cdq,cdm,cdh,elhx,tg,srheat,tg1,ptype,trheat    !,dhgs
     *     ,rhosrf,ma1,tfs,th1,thv1,psk,ps
     *     ,spring,q1,dlwdt,irrig,htirrig
      real*8 fb,fv
!@var rhosrf0 estimated surface air density
      real*8 rhosrf0

      real*8, dimension(:,:), pointer :: trdn,srdn

      real*8 qsat
      integer jr
!@var qg rel. humidity at the ground, defined: total_evap = Cq V (qg-qs)
!@var qg_nsat rel. humidity at non-saturated fraction of soil
      real*8 qg, qg_nsat
c****
c**** fearth    soil covered land fraction (1)
c****
c**** snowe     earth snow amount (kg/m**2)
c**** tearth    earth temperature of first layer (c)
c**** wearth    earth water of first layer (kg/m**2)
c**** aiearth   earth ice of first layer (kg/m**2)
c****
c**** wbare  1-6 water of bare soil layer 1-6 (m)
c**** wvege   0  water of vegetation canopy (m)
c****        1-6 water of vegetated soil layer 1-6 (m)
c**** htbare  0  bare soil layer 0 is unused
c****        1-6 heat content of bare soil layer 1-6 (j m-2)
c**** htvege  0  heat content of vegetation canopy (j m-2)
c****        1-6 heat content of vegetated soil layer 1-6 (j m-2)
c**** snowbv  1  snow depth over bare soil (m)
c****         2  snow depth over vegetated soil (m)
c****

c**** input/output for PBL
      type (t_pbl_args) pbl_args
      real*8 qg_sat,ts,qs

      INTEGER :: ii, ivar, kr
#ifdef TRACERS_DUST
      INTEGER :: n,n1
#endif

      real*8 Ca !@Ca concentration of CO2 at surface (mol/m3)
      real*8 vis_rad, direct_vis_rad, cos_zen_angle
      !integer hemi(1:IM,grid%J_STRT:grid%J_STOP)
      !integer :: JEQUATOR=JM/2
      real*8 :: excess_C_delta
#ifdef WATER_PROPORTIONAL
      real*8, dimension(ntm) :: conc1
      integer :: lpbl,itr
      real*8 :: trconcflx
#endif
#ifdef TRACERS_WATER
      type (ghy_tr_str) :: ghy_tr
#endif

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  ngm                              ) :: sddarr3d
#endif

C****   define local grid
      integer J_0, J_1, J_0H, J_1H, I_0, I_1

!@var end_of_day_flag - hack to pass "end of day" flag to Ent
      logical :: end_of_day_flag

      integer, save :: counter=0
      integer :: dayOfYear

      dayOfYear = modelEclock%getDayOfYear()

      counter = counter + 1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      trdn => atmlnd%flong
      srdn => atmlnd%fshort

    ! dtsurf=dtsrc/nisurf

      spring=-1.
      if((dayOfYear.ge.32).and.(dayOfYear.le.212)) spring=1.

      IF (MOD(Itime,NDAY).eq.0) THEN
        end_of_day_flag = .true.
      else
        end_of_day_flag = .false.
      endif

c****
c**** outside loop over time steps, executed nisurf times every hour
c****
ccc set i,j - independent stuff for tracers
#ifdef TRACERS_ON
      call ghy_tracers_set_step(
#ifdef TRACERS_WATER
     &     ghy_tr
#endif
     &     )
#endif

c****
c**** outside loop over j and i, executed once for each grid point
c****

      !--- at the moment update vegetation every time step
      !hemi(:,max(JEQUATOR+1,J_0):J_1) = 1
      !hemi(:,J_0:min(JEQUATOR,J_1)) = -1
  !    call ent_prescribe_vegupdate(dtsrc/nisurf,entcells,
  !   &     hemi,jday,jyear,.false.,.false.)

       ! moved to daily_earth
!      call update_vegetation_data( entcells,
!     &     im, jm, I_0, I_1, J_0, J_1, jday, jyear )

      !call ent_cell_print(900+counter, entcells)

!!      call update_vegetation_data( entcells,
!!     &     im, jm, I_0, I_1, J_0, J_1, jday, jyear )

      call alloc_pbl_args(pbl_args)

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      pbl_args % moddd = moddd
      pbl_args % ih = 1+modelEclock%getHour()
      pbl_args % ihm = pbl_args%ih+(modelEclock%getDate()-1)*24
#endif

      loop_j: do j=J_0,J_1
  !    hemi=1.
  !    if(j.le.jm/2) hemi=-1.

c**** conditions at the south/north pole
  !    pole= ( j.eq.1 .or. j.eq.jm )

      loop_i: do i=I_0,imaxj(j)

      pbl_args%dtsurf=dtsrc/nisurf
      !pbl_args%hemi = sign(1d0,lat2d(i,j))
c      pbl_args%pole= ( j.eq.1 .or. j.eq.jm )

c****
c**** determine surface conditions
c****
      ptype=fearth(i,j)
      ps=atmlnd%srfp(i,j)
      psk=atmlnd%srfpk(i,j)
      th1=atmlnd%temp1(i,j) !t(i,j,1)
      q1=atmlnd%q1(i,j) !q(i,j,1)
      thv1=th1*(1.+q1*xdelt)
      !pbl_args%tkv=thv1*psk
  !    tkv=thv1*psk

      tfs=tf*ptype
      ma1=atmlnd%am1(i,j)
      !!!qm1=q1*ma1
c     rhosrf=100.*ps/(rgas*tsv)
c     rhosrf=100.*ps/(rgas*tkv)
c****
c**** earth
c****
      if (ptype.le.0.) then
  !      ipbl(4,i,j)=0
        cycle loop_i
      endif
      itype=4

      !!!pr=prec(i,j)/(dtsrc*rhow)
C**** This variable was originally assoicated with super-saturated
C**** large-scale precip, but has not worked for many moons.
C**** If you want to reinstate it, uncomment this calculation.
c      prs=precss(i,j)/(dtsrc*rhow)
      !!!prs=0.
      !!!htpr=eprec(i,j)/dtsrc

ccc tracers variables

      tg1 = tsns_ij(i,j)
      srheat=srdn(i,j)*atmlnd%cosz1(i,j)
      atmlnd%solar(i,j) = atmlnd%solar(i,j) + pbl_args%dtsurf*srheat
c****
c**** boundary layer interaction
c****
      !pbl_args%zs1=.5d-2*rgas*pbl_args%tkv*ma1/atmlnd%p1(i,j)
  !    zs1=.5d-2*rgas*tkv*ma1/pmid(1,i,j)

c**** loop over ground time steps
      tg=tg1+tf
      elhx=lhe
      if(tg1.lt.0.)  elhx=lhs
      pbl_args%qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step
      pbl_args%tg=tg
      pbl_args%elhx=elhx
      !pbl_args%qsol=srheat   ! solar heating
  !    qg_sat=qsat(tg,elhx,ps)  !  replacing with qs from prev step

      qg = qg_ij(i,j)
      ! if ( qg > 999.d0 ) qg = qg_sat
      pbl_args%qg_aver = qg
  !    qg_aver = qg

      pbl_args%tgv=tg*(1.+qg*xdelt)
  !    tgv=tg*(1.+qg*xdelt)

  !    pbl_args%psurf=ps
  !    psurf=ps

c      pbl_args%trhr0 = TRDN(I,J)
  !    trhr0 = TRHR(0,I,J)
      pbl_args%tr4 = atmlnd%gtempr(I,J)**4

      rhosrf0=100.*ps/(rgas*pbl_args%tgv) ! estimated surface density
C**** Obviously there are no ocean currents for earth points, but
C**** variables set for consistency with surfce
  !    pbl_args%uocean=0 ; pbl_args%vocean=0
      pbl_args%ocean = .FALSE.
  !    uocean=0 ; vocean=0

c***********************************************************************
c****
ccc actually PBL needs evap (kg/m^2*s) / rho_air
      pbl_args%evap_max = evap_max_ij(i,j) * 1000.d0 / rhosrf0
      pbl_args%fr_sat = fr_sat_ij(i,j)
  !    evap_max = evap_max_ij(i,j) * 1000.d0 / rhosrf0
  !    fr_sat = fr_sat_ij(i,j)

c**** call tracers stuff
#ifdef TRACERS_ON
      call ghy_tracers_set_cell(i,j,ptype,pbl_args%qg_sat
     &     ,pbl_args%evap_max,pbl_args, q1*ma1
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     &     )
#endif
#ifdef WATER_PROPORTIONAL
#ifdef TRACERS_ON
      call stop_model('GHY: should not get here',255)
#endif
      pbl_args%ntx = 0  ! tracers are activated in PBL, but GHY has none
#endif

      call get_canopy_temperaure(pbl_args%canopy_temperature, i, j)

      call pbl(i,j,1,itype,ptype,pbl_args,atmlnd)

#ifdef GHY_DRYDEP_FIX
      call ghy_tracers_set_cell_stage2(i,j,pbl_args
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     &     )
#endif
c****
      cdm = pbl_args%cm ! cmgs(itype,i,j)
      cdh = pbl_args%ch ! chgs(itype,i,j)
      cdq = pbl_args%cq ! cqgs(itype,i,j)
c***********************************************************************
c**** calculate qs
      qs=pbl_args%qsrf
      ts=pbl_args%tsv/(1.+qs*xdelt)
c**** calculate rhosrf*cdm*ws
      rhosrf=100.*ps/(rgas*pbl_args%tsv)
      rcdmws=cdm*pbl_args%ws*rhosrf
      rcdhws=cdh*pbl_args%ws*rhosrf
      trheat=trdn(i,j)
c***********************************************************************
c****
c  define extra variables to be passed in surfc:
      !!!!pres   = ps
      !veg_pres = ps
      !!!rho    = rhosrf
      !!!vs     = pbl_args%ws
      !!!vs0    = pbl_args%ws0
      !!!tprime = pbl_args%tprime
      !!!qprime = pbl_args%qprime
      !veg_ws = pbl_args%ws
      !!!ch     = cdh
      !veg_ch = cdh
      !!!srht   = srheat
      !veg_srht = srheat
      !!!trht   = trheat
c***********************************************************************
c****
c**** calculate ground fluxes
c     call qsbal
!!! insert htprs here ???

ccc stuff needed for dynamic vegetation
!!! changed units of Ca from mol/m^3 to ppm (also in advnc() )
      select case(land_CO2_bc_flag)
      case(2)
#ifdef TRACERS_GASEXCH_land_CO2
        Ca = trm(i,j,1,n_CO2n)*byaxyp(i,j)/ma1*29.d0/44.d0 *1.d6
!     &       *ps*100.0/gasc/ts
#else
        call stop_model("land_CO2_bc_flag==2 with no C02 tracers",255)
#endif
      case(1)
        Ca = CO2ppm  !*(1.0D-06)*ps*100.0/gasc/ts
      case(0)
        Ca = land_CO2_bc !*(1.0D-06)*ps*100.0/gasc/ts
      end select

!!      vis_rad        = SRVISSURF(i,j)
      vis_rad        = SRVISSURF(i,j)*atmlnd%cosz1(i,j)*.82
      direct_vis_rad = vis_rad * FSRDIR(i,j)
      cos_zen_angle = atmlnd%cosz1(i,j)
!!!!      call ghinij (i,j)
      !call veg_set_cell(i,j)
!i move the following part inside GHY
!i      call ent_get_value( i, j,
!i     &     canopy_holding_capacity=ws_can,
!i     &     canopy_heat_capacity=shc_can,
!i     &     fraction_of_vegetated_soil=fv )
!i      fb = 1.d0 - fv
!!!!      call advnc( entcells(i,j), Ca,
!!!!     &     cos_zen_angle, vis_rad, direct_vis_rad )
!!!!      call evap_limits( vegcell,
!!!!     &     .false., evap_max_ij(i,j), fr_sat_ij(i,j) )

      !call veg_save_cell(i,j)
!!!!      call ghy_save_cell(i,j)
! snow / var lakes problem at this cell (underwater snow in initial data)
!      if (i==47 .and. j==33) then
!        print *,i,j
!      endif
      call get_fb_fv( fb, fv, i, j )
      !write(401,*) "cdh", i,j,cdh
      ijdebug=i*1000+j
      i_earth = i ; j_earth = j ! for snow debug statements
      !ghy_counter = counter

      irrig = 0.d0
      htirrig = 0.d0
#ifdef IRRIGATION_ON
      if (ptype>0.d0)then
         irrig = irrig_water_act(i,j) / ptype
         htirrig = irrig_energy_act(i,j) / ptype
      endif
#endif

      call advnc(
     &     entcells(i,j), Ca,
     &     cos_zen_angle, vis_rad, direct_vis_rad,
     &     Qf_ij(i,j),
     &     w_ij(0:ngm,1:LS_NFRAC,i,j),
     &     ht_ij(0:ngm,1:LS_NFRAC,i,j),
     &     nsn_ij    (1:2, i, j),
     &     dzsn_ij   (1:nlsn, 1:2, i, j),
     &     wsn_ij    (1:nlsn, 1:2, i, j),
     &     hsn_ij    (1:nlsn, 1:2, i, j),
     &     fr_snow_ij(1:2, i, j),
     &     top_index_ij(i, j),
     &     top_dev_ij(i, j),
     &     dz_ij(i,j,1:ngm),
     &     q_ij(i,j,1:imt,1:ngm),
     &     qk_ij(i,j,1:imt,1:ngm),
     &     sl_ij(i,j),
     &     fb,
     &     fv,
     &     prec(i,j)/(dtsrc*rhow),
     &     eprec(i,j)/dtsrc,
     &     precss(i,j)/(dtsrc*rhow),    !0.d0,
     &     0.d0, ! computed inside
     &     irrig,
     &     htirrig,
     &     srheat,
     &     trheat,
     &     pbl_args%tsv/(1.+qs*xdelt),
     &     pbl_args%qsrf,
     &     ps,
     &     rhosrf,
     &     cdh,
     &     q1*ma1,
     &     pbl_args%ws,
     &     pbl_args%ws0,
     &     pbl_args%gusti,
     &     pbl_args%tprime,
     &     pbl_args%qprime,
     &     end_of_day_flag
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     &     )


      call evap_limits(
     &     .false., evap_max_ij(i,j), fr_sat_ij(i,j) )

      !call veg_save_cell(i,j)
      !call ghy_save_cell(i,j)
      if ( fb > 0 ) snowbv(1,i,j)   = snowd(1)
      if ( fv > 0 ) snowbv(2,i,j)   = snowd(2)

      tg1=tsns

c**** computing ground humidity to be used on next time step
      !qg_ij(i,j) = qs  !!! - this seemed to work ok
      !! trying more precise value for qg :  qsat(tg1+tf,elhx,ps)
      qg_sat = qsat(tg1+tf,elhx,ps) ! saturated soil
      qg_nsat = qs              ! non-sat soil, no evap
      if ( rcdhws > 1.d-30 ) then   ! correction to non-sat, due to evap
        qg_nsat = qg_nsat + evap_max_ij(i,j)/(0.001*rcdhws)
      endif
      qg_nsat = min( qg_nsat, qg_sat )
      qg_ij(i,j) = fr_sat_ij(i,j) * qg_sat
     &     + (1.d0 -fr_sat_ij(i,j)) * qg_nsat

ccc Save canopy temperature.
ccc canopy_temp_ij is not used so far... do we need it?
cddd      canopy_temp_ij(i,j) = tp(0,2)  !nyk
cddd
cddd      call get_canopy_temperaure(pbl_args%canopy_temperature, i, j)
cddd
cddd      write(951,*) pbl_args%canopy_temperature,canopy_temp_ij(i,j),i,j


c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      if ( snow_cover_same_as_rad == 0 ) then
        ! recompute snow fraction using different formula
        do ibv=1,2
          call snow_cover(atmlnd%fr_snow_rad(ibv,i,j),
     &         snowbv(ibv,i,j), top_dev_ij(i,j) )
          atmlnd%fr_snow_rad(ibv,i,j) = min (
     &         atmlnd%fr_snow_rad(ibv,i,j), fr_snow_ij(ibv, i, j) )
        enddo
      else
        ! snow fraction same as in snow model
        atmlnd%fr_snow_rad(:,i,j) = fr_snow_ij(:, i, j)
      endif

c**** snowe used in RADIATION
c     snowe(i,j)=1000.*(snowd(1)*fb+snowd(2)*fv)
c workaround for uninitialzed snowd multiply by zero
      atmlnd%snowe(i,j)=1000.*
     &     ( merge( snowd(1)*fb, 0d0, fb > 0 ) +
     &       merge( snowd(2)*fv, 0d0, fv > 0 ) )
      atmlnd%snow(i,j) = atmlnd%snowe(i,j)
      atmlnd%snowfr(i,j) =
     *       ( fb*atmlnd%fr_snow_rad(1,i,j)
     *       + fv*atmlnd%fr_snow_rad(2,i,j) )
      atmlnd%snowdp(i,j) =
     &       ( fb*fr_snow_ij(1,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(1,i,j),1,i,j) )
     &       + fv*fr_snow_ij(2,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(2,i,j),2,i,j) ) )

      do ibv=1,2
        if ( fr_snow_ij(ibv, i, j) > 0.001d0 ) then
          snowd_ij(ibv,i,j) = sum(dzsn_ij(1:nsn_ij(ibv,i,j), ibv, i, j))
        else
          snowd_ij(ibv,i,j) = 0.d0
        endif
      enddo

cddd      if (i==23 .and. j==10) then
cddd        write(755,*) "counter", counter
cddd        write(755,*) "snowe(i,j)", snowe(i,j)
cddd        write(755,*) "snowd(1)",snowd(1)
cddd        write(755,*) "snowd(2)",snowd(2)
cddd        write(755,*) "fb",fb
cddd        write(755,*) "fv",fv
cddd      endif
      !if ( j>23 ) write(995,*) i,j, snowe(i,j)
c**** tearth used only internaly in GHY_DRV
      tearth(i,j) = tbcs
      tsns_ij(i,j) = tsns
c**** wearth+aiearth are used in radiation only
      wearth(i,j)=1000.*( fb*w(1,1)*(1.-fice(1,1)) +
     &     fv*(w(1,2)*(1.-fice(1,2))+w(0,2)*(1.-fice(0,2))) )
      aiearth(i,j)=1000.*( fb*w(1,1)*fice(1,1) +
     &     fv*(w(1,2)*fice(1,2)+w(0,2)*fice(0,2)) )
      atmlnd%gtemp(i,j)=tsns_ij(i,j)
      atmlnd%gtempr(i,j) =tearth(i,j)+tf
      soil_surf_moist(i,j) = 1000.*(fb*w(1,1) + fv*w(1,2))/dz_ij(i,j,1)
      atmlnd%bare_soil_wetness(i,j) =
     &     w(1,1) / ( thets(1,1)*dz_ij(i,j,1) )

#ifdef SCM
      if( SCMopt%Tskin )then
c**** force ground temperature
        atmlnd%gtemp(i,j) = SCMin%Tskin - tf
        atmlnd%gtempr(i,j) = SCMin%Tskin
      endif
      if( SCMopt%sflx )then
C**** impose specified surface heat fluxes
        ashg  = SCMin%shf*pbl_args%dtsurf
        alhg  = SCMin%lhf*pbl_args%dtsurf
        aevap = SCMin%lhf*pbl_args%dtsurf/lhe
      endif
#endif
c**** calculate fluxes using implicit time step for non-ocean points
      atmlnd%uflux1(i,j) = rcdmws*(pbl_args%us) !-uocean)
      atmlnd%vflux1(i,j) = rcdmws*(pbl_args%vs) !-vocean)
c**** accumulate surface fluxes and prognostic and diagnostic quantities
      atmlnd%evapor(i,j)=atmlnd%evapor(i,j)+aevap
      shdt=-ashg
      evhdt=-alhg
C**** calculate correction for different TG in radiation and surface
      dLWDT = pbl_args%dtsurf*
     &     (atmlnd%TRUP_in_rad(I,J) - STBO*(tearth(i,j)+TF)**4)
#ifdef SCM
      if( .not. SCMopt%sfcQrad )then 
c**** zero atmospheric longwave cooling associated with surface
        dLWDT = 0.
      endif
#endif

      atmlnd%dth1(i,j)=-(SHDT+dLWDT)/(sha*ma1)
      atmlnd%dq1(i,j) = aevap/ma1
      atmlnd%sensht(i,j) = atmlnd%sensht(i,j)+SHDT
      atmlnd%latht(i,j) = atmlnd%latht(i,j) + EVHDT
  !    qsavg(i,j)=qsavg(i,j)+qs*ptype

c**** save runoff for addition to lake mass/energy resevoirs
      atmlnd%runo (i,j)= atmlnd%runo (i,j)+ aruns+ arunu
      atmlnd%eruno(i,j)= atmlnd%eruno(i,j)+aeruns+aerunu
c****
      atmlnd%e0(i,j)=atmlnd%e0(i,j) + ae0 - eprec(i,j)
      !e1(i,j,4)=e1(i,j,4)+af1dt

c     redistribute excess_C in the atmosphere (assume the time scale
c     of about a year)
      excess_C_delta = excess_C(i,j)/SECONDS_PER_YEAR*pbl_args%dtsurf
      excess_C(i,j) = excess_C(i,j) - excess_C_delta

      call ghy_diag(i,j,jr,kr,ns,moddsf
     &     ,rcdmws,cdm,cdh,cdq,qg,dlwdt
     &     ,pbl_args, pbl_args%dtsurf
     &     ,excess_C_delta
#ifdef TRACERS_DUST
     &     ,n,n1
#endif
     &     )


c**** update tracers
#ifdef TRACERS_ON
      call ghy_tracers_save_cell(i,j,ptype,pbl_args%dtsurf,rhosrf
     &   ,pbl_args,excess_C_delta
#ifdef INTERACTIVE_WETLANDS_CH4
     & ,tg1,ts-tf,1.d2*abetad/real(nisurf)
#endif
#ifdef TRACERS_WATER
     &     ,ghy_tr
#endif
     & )
#endif


#ifdef WATER_PROPORTIONAL
c calculate fluxes of atmosphere-only water tracers.
c assume 1-way until up/down fluxes of vapor are available
c as a PBL diagnostic.
      conc1(1:ntm) = trm(i,j,1,1:ntm)/
     &     (q1*axyp(i,j)*atmlnd%am1(i,j)+1d-20)
      do itr=1,ntm
        if(aevap.ge.0.) then
          trconcflx = atmlnd%gtracer(itr,i,j)
        else
          trconcflx = conc1(itr)
        endif
        atmlnd%trevapor(itr,i,j) = atmlnd%trevapor(itr,i,j) +
     &       aevap*trconcflx
        atmlnd%trsrfflx(itr,i,j)=atmlnd%trsrfflx(itr,i,j)+
     &       aevap*trconcflx/(pbl_args%dtsurf)
c fill in pbl profile in case it is used to initialize
c another surface type
        do lpbl=1,npbl
          asflx(itype)%trabl(lpbl,itr,i,j)=
     &              conc1(itr)*asflx(itype)%qabl(lpbl,i,j)
        enddo
      enddo ! itr
#endif

      end do loop_i
      end do loop_j


      call soil_spinup_driver

      call dealloc_pbl_args(pbl_args)
      !call dump_ent_C_diags

      ! land water deficit for changing lake fractions
      !!! not working with Ent
!#ifndef uSE_ENT
      call compute_water_deficit
!#endif

#ifdef TRACERS_ON
      call ghy_tracers_finish_step(
#ifdef TRACERS_WATER
     &     ghy_tr
#endif
     &     )
#endif

#ifdef CACHED_SUBDD
C****
C**** Collect high-frequency outputs
C****
      if(ns.eq.nisurf) then ! do only once per physics timestep
C
      call find_groups('aijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
C
      case ('smst')
        where(fearth.gt.0.)
          sddarr2d = soil_surf_moist
        elsewhere
          sddarr2d = 0.
        end where
        call inc_subdd(subdd,k,sddarr2d)
C
      case ('rnft')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          sddarr2d(i,j) = atmlnd%runo(i,j)*fearth(i,j)
        enddo;        enddo
        call inc_subdd(subdd,k,sddarr2d)
      end select
      enddo
      enddo
C
      call find_groups('gijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('GT')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          if(fearth(i,j).gt.0.) then
            sddarr3d(i,j,1:ngm) = gsaveL(i,j,1:ngm,1)
          else
            sddarr3d(i,j,1:ngm) = undef
          end if
        enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
C
      case ('GW')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          if(fearth(i,j).gt.0.) then
            sddarr3d(i,j,1:ngm) = gsaveL(i,j,1:ngm,2)
          else
            sddarr3d(i,j,1:ngm) = undef
          end if
        enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
C
      case ('GI')
        do j=j_0,j_1; do i=i_0,imaxj(j)
          if(fearth(i,j).gt.0.) then
            sddarr3d(i,j,1:ngm) = gsaveL(i,j,1:ngm,3)
          else
            sddarr3d(i,j,1:ngm) = undef
          end if
        enddo;         enddo
        call inc_subdd(subdd,k,sddarr3d)
      end select
      enddo
      enddo

      endif ! once per physics timestep
#endif

      return
      end subroutine earth

      subroutine dump_ent_C_diags
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,READT_PARALLEL
      USE DOMAIN_DECOMP_1D, only : WRITET_PARALLEL
      use ent_mod, only: entcelltype_public, debug_carbon
      use ent_com, only : entcells
      !---
      real*8, dimension(16,im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     total,
     &          C_lab, C_fol, C_sw, C_hw, C_froot, C_croot, C_soil
      real*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     total_sum, C_soil_sum, C_lab_sum, C_fol_sum
      integer, save :: counter = 0
      integer, save :: fc = 1000
      character*80 :: title
      integer :: k, i, j
      integer :: I_1, I_0, J_1, J_0
      integer :: I_1H, I_0H, J_1H, J_0H

      CALL getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     *               J_STRT=J_0,       J_STOP=J_1)

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
 

      if ( mod(counter,48*30) == 0 ) then

        fc = fc + 1
        do j=J_0,J_1
          do i=I_0,I_1
            call debug_carbon(entcells(i,j), total(:,i,j),
     &           C_lab(:,i,j), C_fol(:,i,j), C_sw(:,i,j), C_hw(:,i,j),
     &           C_froot(:,i,j), C_croot(:,i,j), C_soil(:,i,j))
          enddo
        enddo

        do k=1,16
          write(title,*) "total ",k
          call WRITET_PARALLEL(grid,fc,"foo",total(k,:,:),title)
          write(title,*) "C_lab ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_lab(k,:,:),title)
          write(title,*) "C_fol ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_fol(k,:,:),title)
          write(title,*) "C_sw ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_sw(k,:,:),title)
          write(title,*) "C_hw ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_hw(k,:,:),title)
          write(title,*) "C_froot ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_froot(k,:,:),title)
          write(title,*) "C_croot ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_croot(k,:,:),title)
          write(title,*) "C_soil ",k
          call WRITET_PARALLEL(grid,fc,"foo",C_soil(k,:,:),title)
         enddo

         total_sum = 0.d0
         C_soil_sum = 0.d0
         C_lab_sum = 0.d0
         C_fol_sum = 0.d0
         do k=1,16
           total_sum(:,:) = total_sum(:,:) + total(k,:,:)
           C_soil_sum(:,:) = C_soil_sum(:,:) + C_soil(k,:,:)
           C_lab_sum(:,:) = C_lab_sum(:,:) + C_lab(k,:,:)
           C_fol_sum(:,:) = C_fol_sum(:,:) + C_fol(k,:,:)
         enddo

         write(title,*) "total sum"
         call WRITET_PARALLEL(grid,fc,"foo",total_sum(:,:),title)
         write(title,*) "C_soil sum"
         call WRITET_PARALLEL(grid,fc,"foo",C_soil_sum(:,:),title)
         write(title,*) "C_lab sum"
         call WRITET_PARALLEL(grid,fc,"foo",C_lab_sum(:,:),title)
         write(title,*) "C_fol sum"
         call WRITET_PARALLEL(grid,fc,"foo",C_fol_sum(:,:),title)

      endif

      counter = counter + 1

      end subroutine dump_ent_C_diags

c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine ghy_diag(i,j,jr,kr,ns,moddsf
     &     ,rcdmws,cdm,cdh,cdq,qg,dlwdt, pbl_args, dtsurf
     &     ,excess_C_delta
#ifdef TRACERS_DUST
     &     ,n,n1
#endif
     &     )

      use diag_com , only : itearth,aij=>aij_loc
     *     ,tsfrez=>tsfrez_loc,tdiurn,jreg, ij_psoil,ij_vsfr,ij_bsfr
     *     ,ij_rune, ij_arunu, ij_pevap, ij_beta
     *     ,j_erun,j_run,j_type
     *     ,ij_gpp,ij_ipp,ij_rauto,ij_clab,ij_lai,ij_ecvf
     *     ,ij_soilresp,ij_soilCpoolsum,ij_gbetat,ij_gbetpen,ij_gevppen
     *     ,ij_gbsbet,ij_gbsevp,ij_gbst,ij_gvst,ij_gbvswt,ij_gconatm
     *     ,ij_gconcan,ij_gdcevp,ij_gwtbl,ij_gvswet,ij_gwcevp
     *     ,tf_day1,tf_last
     &     ,ij_aflmlt,ij_aeruns,ij_aerunu,ij_fveg
     &     ,ij_htsoil,ij_htsnow,ij_aintrcp
     &     ,ij_evapsn,ij_irrW, ij_irrE, ij_landcarbon
#if (defined HEALY_LM_DIAGS)
     &     ,ij_crops,j_crops,CROPS_DIAG
#endif
      use constant, only : tf,lhe
#ifdef TRACERS_DUST
     &     ,grav
#endif
      use fluxes, only : atmlnd,nisurf
      use model_com, only : modelEclock
      use model_com, only : dtsrc,nday,itime
      use TimeConstants_mod, only: DAYS_PER_YEAR, INT_HOURS_PER_DAY
      use DOMAIN_DECOMP_ATM, only : grid
      use geom, only : axyp,lat2d
      use sle001, only :
     &     tp
     &    ,fv,fb,atrg,ashg,alhg
     &    ,abetad,abetav,abetat
     &    ,abetap,abetab,abeta
     &    ,acna,acnc,agpp,aipp,arauto,aclab,alai,asoilresp,asoilCpoolsum
     &    ,aevap,aevapw,aevapd,aevapb,aevapbs,aevapvs
     &    ,aruns,arunu,aeruns,aerunu,aflmlt,aintercep
     &    ,aepc,aepb,aepp,zw,tbcs,tsns
     &    ,qs,ts,ngr=>n,ht,hsn,fr_snow,nsn
     &    ,tg2av,wtr2av,ace2av
     &    ,tg_L,wtr_L,ace_L
     &    ,airrig,aeirrig,alandC
      use ent_com, only : excess_C
      use ghy_com, only : gdeep, gsaveL, fearth
      USE CLOUDS_COM, only : DDMS

      use pbl_drv, only : t_pbl_args

#ifdef TRACERS_DUST
      USE tracer_com,ONLY : Ntm_dust,n_clay
#ifdef TRACERS_DRYDEP
      use OldTracer_mod, only: dodrydep
#endif
#endif

#ifdef ENT_DEBUG_DIAGS
      use diag_com , only : ij_ent_debug
      use sle001, only : ent_debug_buf
      use ent_debug_mod, only : SIZE_ENT_DEBUG
#endif


      implicit none
      integer, intent(in) :: i,j,ns,moddsf
      real*8, intent(in) :: rcdmws,cdm,cdh,cdq,qg,dlwdt
      type (t_pbl_args) :: pbl_args
      real*8, intent(in) :: dtsurf, excess_C_delta

      real*8 timez
      real*8 trhdt,tg1,shdt,ptype,srheat,srhdt
      real*8 warmer,spring,trheat,evhdt
      integer, parameter :: itype=4
      integer :: kr,jr,ih,ihm,k
#ifdef TRACERS_DUST
      INTEGER :: n,n1
#endif

ccc the following values are returned by PBL
      real*8 us,vs,ws,psi,dbl,khs,ug,vg,wg,gusti,qsrf,tsv,ps,elhx
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
     &     ,wsgcm,wspdf
#endif
      real*8, external :: qsat
      integer :: dayOfYear

      dayOfYear = modelEclock%getDayOfYear()

      us = pbl_args%us
      vs = pbl_args%vs
      ws = pbl_args%ws
      psi = pbl_args%psi
      dbl = pbl_args%dbl
      khs = pbl_args%khs
      ug = pbl_args%ug
      vg = pbl_args%vg
      wg = pbl_args%wg
      gusti = pbl_args%gusti
      qsrf = pbl_args%qsrf ; tsv = pbl_args%tsv ; ps = pbl_args%psurf
      elhx = pbl_args%elhx
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      wsgcm=pbl_args%wsgcm
      wspdf=pbl_args%wspdf
#endif

      timez=dayOfYear+(mod(itime,nday)+(ns-1.)/nisurf)/nday ! -1 ??
      if(dayOfYear.le.31) timez=timez+DAYS_PER_YEAR

      ih=1+modelEclock%getHour()
      ihm = ih+(modelEclock%getDate()-1)*INT_HOURS_PER_DAY

      spring=-1.
      if((dayOfYear.ge.32).and.(dayOfYear.le.212)) spring=1.

      if(lat2d(i,j).lt.0.)  then
         warmer=-spring
       else
         warmer=spring
      end if

      ptype=fearth(i,j)
      jr=jreg(i,j)

      srheat=atmlnd%fshort(i,j)*atmlnd%cosz1(i,j)
      srhdt=srheat*dtsurf
      trheat=atmlnd%flong(i,j)

      !tg1=tbcs
      tg1=tsns
      shdt=-ashg
      evhdt=-alhg

#ifdef ENT_DEBUG_DIAGS
      aij(i,j, ij_ent_debug:ij_ent_debug+SIZE_ENT_DEBUG-1)=
     &     aij(i,j, ij_ent_debug:ij_ent_debug+SIZE_ENT_DEBUG-1)
     &     + ent_debug_buf(:)*ptype
#endif

      aij(i,j,ij_psoil)=aij(i,j,ij_psoil)+ptype/nisurf
      aij(i,j,ij_fveg)=aij(i,j,ij_fveg)+fv/nisurf
      aij(i,j,ij_vsfr)=aij(i,j,ij_vsfr)+fv*ptype/nisurf
      aij(i,j,ij_bsfr)=aij(i,j,ij_bsfr)+fb*ptype/nisurf
      aij(i,j,ij_gbsevp)=aij(i,j,ij_gbsevp)+aevapb*ptype*fb
      aij(i,j,ij_gdcevp)=aij(i,j,ij_gdcevp)+aevapd*ptype*fv
      aij(i,j,ij_gwcevp)=aij(i,j,ij_gwcevp)+aevapw*ptype*fv
      aij(i,j,ij_gbsbet)=aij(i,j,ij_gbsbet)+(abetab/nisurf)*fb*ptype
      aij(i,j,ij_gbetpen)=aij(i,j,ij_gbetpen)+(abetap/nisurf)*ptype
      aij(i,j,ij_gbvswt)=aij(i,j,ij_gbvswt)+(abeta/nisurf)*ptype
      aij(i,j,ij_gconatm)=aij(i,j,ij_gconatm)+(acna/nisurf)*ptype
      aij(i,j,ij_gconcan)=aij(i,j,ij_gconcan)+(acnc/nisurf)*ptype*fv
      aij(i,j,ij_gpp)=aij(i,j,ij_gpp)+agpp*ptype
      aij(i,j,ij_ipp)=aij(i,j,ij_ipp)+aipp*ptype
      aij(i,j,ij_rauto)=aij(i,j,ij_rauto)+arauto*ptype
      aij(i,j,ij_clab)=aij(i,j,ij_clab)+(aclab/nisurf)*ptype
      aij(i,j,ij_lai)=aij(i,j,ij_lai)+(alai/nisurf)*ptype
      aij(i,j,ij_soilresp)=aij(i,j,ij_soilresp)+asoilresp*ptype
      aij(i,j,ij_soilCpoolsum)=aij(i,j,ij_soilCpoolsum)
     &     + (asoilCpoolsum/nisurf)*ptype
      aij(i,j,ij_landcarbon)=aij(i,j,ij_landcarbon)
     &     + ( (alandC+excess_C(i,j) )/nisurf)*ptype
      aij(i,j,ij_ecvf)=aij(i,j,ij_ecvf)+excess_C_delta*ptype
#ifdef HEALY_LM_DIAGS 
      aij(i,j,ij_crops)=aij(i,j,ij_crops)+CROPS_DIAG(i,j)*100.
#endif
      aij(i,j,ij_gvswet)=aij(i,j,ij_gvswet)+(abetav/nisurf)*fv*ptype
      aij(i,j,ij_gbetat)=aij(i,j,ij_gbetat)+(abetat/nisurf)*fv*ptype
      aij(i,j,ij_gevppen)=aij(i,j,ij_gevppen)+aepp*ptype
      aij(i,j,ij_evapsn)=aij(i,j,ij_evapsn)+(aevapbs+aevapvs)*ptype
      if (moddsf.eq.0) then
!      Temperatures of bare soil layers
        do k=1,6
          aij(i,j,ij_gbst+k-1)=aij(i,j,ij_gbst+k-1)+tp(k,1)*ptype*fb
        end do
!      Temperatures of canopy and vegetated soil layers
        do k=0,6
          aij(i,j,ij_gvst+k)=aij(i,j,ij_gvst+k)+tp(k,2)*ptype*fv
        end do
!     Water table depth
        aij(i,j,ij_gwtbl)=aij(i,j,ij_gwtbl)+(fb*zw(1)+fv*zw(2))*ptype
      end if
ccc accumulate total heat storage
      if (moddsf.eq.0) then
        aij(i,j,ij_htsoil)=aij(i,j,ij_htsoil) +
     &       (fb*sum(ht(1:ngr,1)) + fv*sum(ht(0:ngr,2)))*ptype
        aij(i,j,ij_htsnow)=aij(i,j,ij_htsnow)
     &       + (fb*fr_snow(1)*sum(hsn(1:nsn(1),1))
     &         +fv*fr_snow(2)*sum(hsn(1:nsn(2),2)))*ptype
      endif
      trhdt=trheat*dtsurf-atrg
      atmlnd%trheat(i,j) = atmlnd%trheat(i,j) + trhdt
c           for radiation find composite values over earth
c           for diagnostic purposes also compute gdeep 1 2 3
      !!!call retp2 (tg2av,wtr2av,ace2av)
      atmlnd%gtemp2(i,j) = tg2av
      gdeep(i,j,1)=tg2av
      gdeep(i,j,2)=wtr2av
      gdeep(i,j,3)=ace2av
      gsaveL(i,j,:,1)=tg_L(:)
      gsaveL(i,j,:,2)=wtr_L(:)
      gsaveL(i,j,:,3)=ace_L(:)

      aij(i,j,ij_rune)=aij(i,j,ij_rune)+aruns*ptype
      aij(i,j,ij_arunu)=aij(i,j,ij_arunu)+arunu*ptype
      aij(i,j,ij_aeruns)=aij(i,j,ij_aeruns)+aeruns*ptype
      aij(i,j,ij_aerunu)=aij(i,j,ij_aerunu)+aerunu*ptype
      aij(i,j,ij_pevap)=aij(i,j,ij_pevap)+(aepc+aepb)*ptype
      aij(i,j,ij_aflmlt)=aij(i,j,ij_aflmlt)+aflmlt*ptype
      aij(i,j,ij_aintrcp)= aij(i,j,ij_aintrcp)+aintercep*ptype*fv
#ifdef IRRIGATION_ON
      aij(i,j,ij_irrW)=aij(i,j,ij_irrW)+airrig*ptype
      aij(i,j,ij_irrE)=aij(i,j,ij_irrE)+aeirrig*ptype
#else
c      aij(i,j,ij_irrW)=0.d0
c      aij(i,j,ij_irrE)=0.d0
c      aij(i,j,ij_irrW_tot)=0.d0
#endif
      if ( warmer >= 0 ) then
        if(ts.lt.tf) tsfrez(i,j,tf_day1)=timez
        tsfrez(i,j,tf_last)=timez
      else
        if ( tsfrez(i,j,tf_last)+.03 >= timez .and. ts >= tf )
     $       tsfrez(i,j,tf_last)=timez
      endif

      if(tg1.lt.tdiurn(i,j,1)) tdiurn(i,j,1)=tg1
      if(tg1.gt.tdiurn(i,j,2)) tdiurn(i,j,2)=tg1
      if(ts.lt.tdiurn(i,j,3)) tdiurn(i,j,3)=ts
      if(ts.gt.tdiurn(i,j,4)) tdiurn(i,j,4)=ts

c**** quantities accumulated for regions in diagj
      call inc_areg(i,j,jr,j_erun ,(aeruns+aerunu)*ptype)
      call inc_areg(i,j,jr,j_run  ,  (aruns+arunu)*ptype)
#ifdef RJH_LM_CROPS
      call inc_areg(i,j,jr,j_crops,CROPS_DIAG(i,j)*ptype)
#endif

c**** quantities accumulated for latitude-longitude maps in diagij
      aij(i,j,ij_beta)=aij(i,j,ij_beta)+(abetad/nisurf)*fv*ptype
c      if ( moddsf == 0 ) then
chyd       aij(i,j,ij_arunu)=aij(i,j,ij_arunu)
chyd      *  +   (40.6*psoil+.72*(2.*(tss-tfs)-(qsatss-qss)*lhe/sha))
c      endif

c**** quantities accumulated for surface type tables in diagj
      call inc_aj(i,j,itearth,j_erun ,(aeruns+aerunu)*ptype)
      call inc_aj(i,j,itearth,j_run  ,  (aruns+arunu)*ptype)
#ifdef HEALY_LM_DIAGS
      call inc_aj(i,j,itearth,j_crops  ,  CROPS_DIAG(i,j)*ptype)
#endif

c**** quantities accumulated for subdd
      !R_acc(I,J)=R_acc(I,J)+(aruns+arunu)*ptype

      end subroutine ghy_diag

c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************
c***********************************************************************

      subroutine snow_cover( fract_snow, snow_water, top_dev )
!@sum computes snow cover from snow water eq. and topography
!@var fract_snow snow cover fraction (0-1)
!@var snow_water snow water equivalent (m)
!@var top_dev standard deviation of the surface elevation
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use constant, only : teeny
      real*8, intent(out) :: fract_snow
      real*8, intent(in) :: snow_water, top_dev

      ! using formula from the paper by A. Roesch et al
      ! (Climate Dynamics (2001), 17: 933-946)
      fract_snow =
ccc     $     .95d0 * tanh( 100.d0 * snow_water ) *
ccc                               currently using only topography part
     $     sqrt ( 1000.d0 * snow_water /
     $     (1000.d0 * snow_water + teeny + snow_cover_coef * top_dev) )

      end subroutine snow_cover


      subroutine init_LSM(dtsurf,redogh,inisnow,inilake,istart)
      real*8, intent(in) :: dtsurf
      integer, intent(in) :: istart
      logical, intent(in) :: redogh, inisnow, inilake

      call init_gh(dtsurf,istart)


      call init_veg( istart, redogh )
      call init_land_surface(redogh,inisnow,inilake,istart)

#ifdef IRRIGATION_ON
      call init_irrigmod()
#endif

      end subroutine init_LSM


      subroutine init_gh(dtsurf,istart)
c**** modifications needed for split of bare soils into 2 types
      use Dictionary_mod, only : sync_param, get_param
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use fluxes, only : focean
      use diag_com, only : npts,icon_wtg,icon_htg,conpt0
      use sle001, only : hl0, dt, 
     &     minGroundTemperature,  maxGroundTemperature
      use snow_model, only : minSnowTemperature
      use ghy_com
      use snow_drvm, only : snow_cover_coef2=>snow_cover_coef
     &     ,snow_cover_same_as_rad
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      implicit none

      real*8, intent(in) :: dtsurf
      integer, intent(in) :: istart
      integer fid
      logical :: qcon(npts)
      integer i, j
      logical ghy_data_missing
      character conpt(npts)*10
c****
!@dbparam ghy_default_data if == 1 reset all GHY data to defaults
!@+ (do not read it from files)
      integer :: ghy_default_data = 0

C**** define local grid
      integer I_0, I_1, J_0, J_1
      integer I_0H, I_1H, J_0H, J_1H
      logical present_land

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

c**** set conservation diagnostics for ground water mass and energy
      conpt=conpt0
      conpt(4)="EARTH"
      qcon=(/ .false., .false., .false., .true., .false., .false.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,conpt,"GRND WTR","(kg/m^2)        ",
     *     "(10^-9 kg/s/m^2)",1d0,1d9,icon_wtg)
      qcon=(/ .false., .false., .false., .true., .false., .false.
     $     , .false., .false., .true., .false., .false./)
      call set_con(qcon,conpt,"GRND ENG","(10**6 J/m^2)   ",
     *     "(10^-3 W/m^2)   ",1d-6,1d3,icon_htg)

c**** read rundeck parameters
      call sync_param( "snow_cover_coef", snow_cover_coef )
! hack. snow_cover_coef should be moved to snow_drvm
      snow_cover_coef2 = snow_cover_coef
      call sync_param( "snow_cover_same_as_rad", snow_cover_same_as_rad)
      call sync_param( "wsn_max", wsn_max )
      call sync_param( "ghy_default_data", ghy_default_data )
      !call  get_param( "variable_lk", variable_lk )
      !call  get_param( "init_flake", init_flake )
      call sync_param( "land_CO2_bc_flag", land_CO2_bc_flag )
      call sync_param( "land_CO2_bc", land_CO2_bc )

      call sync_param( "minGroundTemperature", minGroundTemperature)
      call sync_param( "maxGroundTemperature", maxGroundTemperature)
      minSnowTemperature = minGroundTemperature

      call sync_param( "sbgc_spinup", sbgc_spinup)

c**** read land surface parameters or use defaults
      if ( ghy_default_data == 0 ) then ! read from files

c**** read soils parameters
        if(file_exists('SOIL')) then
          fid = par_open(grid,'SOIL','read')
          call read_dist_data(grid,fid,'dz',dz_ij)
          call read_dist_data(grid,fid,'q',q_ij)
          call read_dist_data(grid,fid,'qk',qk_ij)
          call read_dist_data(grid,fid,'sl',sl_ij)
          call par_close(grid,fid)
        endif

c**** read topmodel parameters
        if(file_exists('TOP_INDEX')) then
          fid = par_open(grid,'TOP_INDEX','read')
          top_index_ij = 0.     ! in case top_index is absent in the file
          call read_dist_data(grid,fid,'top_index',top_index_ij)
          call read_dist_data(grid,fid,'top_dev',top_dev_ij)
          call par_close(grid,fid)
        endif

      else  ! reset to default data
        if ( istart>0 .and. istart<8 ) then ! reset all
          call reset_gh_to_defaults( .true. )
        else   ! do not reset ghy prognostic variables
          call reset_gh_to_defaults( .false. )
        endif
      endif

c**** time step for ground hydrology
      dt=dtsurf

      ! no need to continue computations for postprocessing
      if (istart.le.0) return

c**** check whether ground hydrology data exist at this point.
      ghy_data_missing = .false.
      do j=J_0,J_1
        do i=I_0,I_1
          present_land = .false.
          if (variable_lk==0) then
            if ( fearth(i,j) > 0.d0 ) present_land = .true.
          else
            if ( focean(i,j) < 1.d0 ) present_land = .true.
          endif
          if ( present_land ) then
            if ( top_index_ij(i,j).eq.-1. ) then
              print *,"No top_index data: i,j=",i,j,top_index_ij(i,j)
              ghy_data_missing = .true.
            end if
            if ( sum(dz_ij(i,j,1:ngm)).eq.0 ) then
              print *, "No soil data: i,j=",i,j,dz_ij(i,j,1:ngm)
              ghy_data_missing = .true.
            endif
            if(allocated(earth_tp)) then
              ! cold-start from temperature and relative wetness.
              ! need to insert appropriate checks here.
            else
              if (w_ij(1,1,i,j)<1.d-10 .and. w_ij(1,2,i,j)<1.d-10) then
                print*,"No gh data in restart file: i,j=",i,j,
     &             w_ij(:,1,i,j),w_ij(:,2,i,j)
                ghy_data_missing = .true.
              endif
            endif
          end if
        enddo
      enddo
      if ( ghy_data_missing ) then
        write(6,*) 'Ground Hydrology data is missing at some pts'
        write(6,*) 'If you have a non-standard land mask, please'
        write(6,*) 'consider using extended GH data and rfs file.'
        call stop_model(
     &       'Ground Hydrology data is missing at some cells',255)
      endif

      call hl0

      return
      end subroutine init_gh


!******************************************************************


      subroutine init_veg( istart, redogh )
      use constant, only : twopi
      use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
      use model_com, only : modelEclock
      use model_com, only : itime,nday
      use ent_drv, only : init_module_ent
      integer, intent(in) :: istart
      logical, intent(in) :: redogh
      !--- local
      integer year, dayOfYear

c**** cosday, sinday should be defined (reset once a day in daily_earth)
      dayOfYear=1+mod(itime/nday,int(EARTH_DAYS_PER_YEAR))

      call modelEclock%get(year=year, dayOfYear=dayOfYear)
      CALL init_module_ent(istart.le.2, dayOfYear, year, istart)

      end subroutine init_veg

!******************************************************************

      subroutine init_land_surface(redogh,inisnow,inilake,istart)
      !use veg_drv, only : spgsn
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use Dictionary_mod, only : sync_param, get_param
      use constant, only : tf, lhe, rhow, shw_kg=>shw
      use ghy_com
      use model_com, only : itime
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin
#endif
      use fluxes, only : atmlnd,focean, flice
      use ent_com, only : entcells
      use ent_mod
      use ent_drv, only : map_ent2giss !YKIM-temp hack
#ifdef TRACERS_WATER
      use OldTracer_mod, only: tr_wd_TYPE, nWATER, itime_tr0, needtrs
      use tracer_com, only : NTM
#endif
      use rad_com, only : snoage
      implicit none
      integer, intent(in) :: istart
      logical, intent(in) :: redogh, inisnow, inilake

      !--- local
      real*8, external :: qsat
      real*8, parameter :: spgsn=.1d0 !@var spgsn specific gravity of snow
      real*8 :: ws_can, shc_can, ht_cap_can, fice_can, aa
      real*8 :: fb, fv
      integer :: reset_canopy_ic=0, reset_snow_ic=0, do_IC_fixups=0
  !!  integer init_flake
      logical present_land
#ifdef TRACERS_WATER
      real*8 trsoil_tot,wsoil_tot,fm,height_can
      integer n
#endif
      integer I_0, I_1, J_0, J_1, i, j, ibv
      logical :: my_inisnow

      my_inisnow = inisnow
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1 )

      call sync_param( "reset_canopy_ic", reset_canopy_ic )
      call sync_param( "reset_snow_ic", reset_snow_ic )
      call  get_param( "variable_lk", variable_lk )
  !!  call  get_param( "init_flake", init_flake )
      if (istart < 9) do_IC_fixups = 1

c**** recompute ground hydrology data if necessary (new soils data)
      if (redogh .or.
           ! new_io_soils checks whether the soil IC file has
           ! IC in intensive units, but has no mechanism to
           ! set the redogh flag.  Check allocation status instead.
     &     allocated(earth_tp)
     &     ) then

        my_inisnow = .false. ! deal with snow here and do not touch it later
        do j=J_0,J_1
          do i=I_0,I_1
            w_ij(:,:,i,j)=0.d0
            ht_ij(:,:,i,j)=0.d0
            !snowbv(:,i,j)=0.d0
            if ( focean(i,j) >= 1.d0 ) cycle
            if ( fearth(i,j) <= 0.d0 .and. variable_lk==0 ) cycle
            call ent_get_exports( entcells(i,j),
     &           canopy_heat_capacity=shc_can,
     &           canopy_max_H2O=ws_can )
            call tp_sat_2_ht_w(
     &           w_ij(:,:,i,j), ht_ij(:,:,i,j),
     &           nsn_ij    (1:2, i, j),
     &           dzsn_ij   (1:nlsn, 1:2, i, j),
     &           wsn_ij    (1:nlsn, 1:2, i, j),
     &           hsn_ij    (1:nlsn, 1:2, i, j),
     &           fr_snow_ij(1:2, i, j),

     &           earth_sat(:,:,i,j), earth_ice(:,:,i,j),
     &           earth_tp(:,:,i,j),  snowbv(:,i,j),
     &           tsnowtop(i,j),

     &           ws_can, shc_can,
     &           q_ij(i,j,:,:), dz_ij(i,j,:)
     &           )

            atmlnd%snowe(i,j) = 1000.*snowbv(1,i,j)
            tearth(i,j) = earth_tp(1,1,i,j)
            atmlnd%bare_soil_wetness(i,j) = earth_sat(1,1,i,j)
            aiearth(i,j) = 0. ! not used?
            snoage(3,i,j) = 0.
            evap_max_ij(i,j) = 1d30 ! trigger initialization of qg_ij
            fr_sat_ij(i,j) = 0.
            qg_ij(i,j) = 0.
            tsns_ij(i,j) = tearth(i,j)
          end do
        end do
        write (*,*) 'ground hydrology data was made from ground data'
        if(allocated(earth_tp)) then
          deallocate(earth_tp, earth_sat, earth_ice, tsnowtop)
        endif
      end if

      if (do_IC_fixups == 1) then
c**** set vegetation temperature to tearth(i,j) if requested
c**** in this case also set canopy water and water tracers to 0
      if ( reset_canopy_ic .ne. 0) then
        do j=J_0,J_1
          do i=I_0,I_1
            present_land = .false.
            if (variable_lk==0) then
              if ( fearth(i,j) > 0.d0 ) present_land = .true.
            else
              if ( focean(i,j) < 1.d0 ) present_land = .true.
            endif

            if( .not. present_land ) cycle

            call get_fb_fv( fb, fv, i, j )
            if ( fv <= 0.d0 ) cycle

            !!! probably will not work
         !!!call stop_model("reset_canopy_ic not implemented for Ent",255)
            call ent_get_exports( entcells(i,j),
     &           canopy_heat_capacity=ht_cap_can )
            ! assume canopy completey dry
            w_ij(0,2,i,j) = 0.d0
            fice_can = 0.d0
            call temperature_to_heat( ht_ij(0,2,i,j),
     &           tsns_ij(i,j), fice_can, w_ij(0,2,i,j), ht_cap_can )
#ifdef TRACERS_WATER
            tr_w_ij(:,0,2,i,j) = 0.d0
#endif
          enddo
        enddo
      endif



c**** if not initialized yet, set evap_max_ij, fr_sat_ij, qg_ij
c**** to something more appropriate
      do j=J_0,J_1
        do i=I_0,I_1
          !if ( fearth(i,j) .le. 0.d0 ) cycle
          if ( focean(i,j) >= 1.d0 ) cycle
          if( evap_max_ij(i,j) < .9d0 ) cycle
          qg_ij(i,j) = qsat(tsns_ij(i,j)+tf,lhe,atmlnd%srfp(i,j))
          fr_sat_ij(i,j) = 0.d0
          evap_max_ij(i,j) = 0.d0
        enddo
      enddo


c**** the following is needed for very old restart files only
c**** (i.e. the files which contained snow as part of 1st soil layer)
      if (my_inisnow) then
        do j=J_0,J_1
          do i=I_0,I_1
            nsn_ij(:,i,j)     = 1
            dzsn_ij(:,:,i,j)  = 0.d0
            wsn_ij(:,:,i,j)   = 0.d0
            hsn_ij(:,:,i,j)   = 0.d0
            fr_snow_ij(:,i,j) = 0.d0
            if ( focean(i,j) >= 1.d0 ) cycle
            if ( fearth(i,j) <= 0.d0 .and. variable_lk==0 ) cycle

           call set_snow1( w_ij (1:ngm,1:LS_NFRAC,i,j) ,
     &           ht_ij(1:ngm,1:LS_NFRAC,i,j) ,
     &           snowbv(1:2,i,j), q_ij(i,j,:,:), dz_ij(i,j,:),
     &           nsn_ij    (1:2, i, j),
     &           dzsn_ij   (1:nlsn, 1:2, i, j),
     &           wsn_ij    (1:nlsn, 1:2, i, j),
     &           hsn_ij    (1:nlsn, 1:2, i, j),
     &           fr_snow_ij(1:2, i, j) )


          end do
        end do
      end if

c**** fix initial conditions for soil water if necessry
        do j=J_0,J_1
          do i=I_0,I_1
            if ( focean(i,j) >= 1.d0 ) cycle
            if ( fearth(i,j) <= 0.d0 .and. variable_lk==0 ) cycle
            call fix_water_ic( w_ij(:,:,i,j),
     &           q_ij(i,j,:,:), dz_ij(i,j,:), i, j )
          enddo
        enddo

c**** fix initial conditions for soil heat if necessry
        do j=J_0,J_1
          do i=I_0,I_1
            if ( focean(i,j) >= 1.d0 ) cycle
            if ( fearth(i,j) <= 0.d0 .and. variable_lk==0 ) cycle
            !call stop_model("fix*_ic not implemented for Ent",255)
            !shc_can = 1.d30 ! i.e. don''t check canopy heat
            call ent_get_exports( entcells(i,j),
     &           canopy_heat_capacity=shc_can )

            call fix_heat_ic(w_ij(:,:,i,j), ht_ij(:,:,i,j),
     &           shc_can,
     &           q_ij(i,j,:,:), dz_ij(i,j,:), i, j )
          enddo
        enddo

c**** remove all land snow from initial conditions
c**** (useful when changing land/vegetation mask)
      if ( reset_snow_ic .ne. 0 ) then
        do j=J_0,J_1
          do i=I_0,I_1
            nsn_ij(:, i, j) = 1
            wsn_ij(:, :, i, j) = 0.d0
            hsn_ij(:, :, i, j) = 0.d0
            dzsn_ij(:, :, i, j) = 0.d0
            fr_snow_ij(:, i, j) = 0.d0
          enddo
        enddo
      endif

!!! hack - remove underwater snow + snow in landice boxes (dealt with separately)
!!! (should not be present in restart file in the first place!)
      do j=J_0,J_1
        do i=I_0,I_1
          if ( fearth(i,j) == 0.d0 ) then
!            if ( maxval(fr_snow_ij(:,i,j)) > 0.d0 )
!      &           print *,"removing snow from ",i,j,
!      &           " : cell under lake or land ice"
            nsn_ij(:, i, j) = 1
            wsn_ij(:, :, i, j) = 0.d0
            hsn_ij(:, :, i, j) = 0.d0
            dzsn_ij(:, :, i, j) = 0.d0
            fr_snow_ij(:, i, j) = 0.d0
          endif
        enddo
      enddo
      end if ! do_IC_fixups

c**** set snow fraction for albedo computation (used by RAD_DRV.f)
      atmlnd%fr_snow_rad(:,:,:) = 0.d0
      do j=J_0,J_1
        do i=I_0,I_1
          if ( fearth(i,j) > 0.d0 ) then
            do ibv=1,2
              call snow_cover(atmlnd%fr_snow_rad(ibv,i,j),
     &             snowbv(ibv,i,j), top_dev_ij(i,j) )
              atmlnd%fr_snow_rad(ibv,i,j) = min (
     &             atmlnd%fr_snow_rad(ibv,i,j), fr_snow_ij(ibv,i,j) )
            enddo
          endif
        enddo
      enddo

      ! initialize underwater fraction for variable lakes
      if ( inilake .and. variable_lk > 0 )
     &     call init_underwater_soil


      ! land water deficit for changing lake fractions
      !!! not working with Ent
!#ifndef uSE_ENT
      call compute_water_deficit
!#endif

c**** set gtemp array
      do j=J_0,J_1
        do i=I_0,I_1
          if (fearth(i,j).gt.0) then
            call get_fb_fv( fb, fv, i, j )
            atmlnd%snow(i,j)=atmlnd%snowe(i,j)
            atmlnd%snowfr(i,j) =
     *           ( fb*atmlnd%fr_snow_rad(1,i,j)
     *           + fv*atmlnd%fr_snow_rad(2,i,j) )
            atmlnd%snowdp(i,j) =
     &       ( fb*fr_snow_ij(1,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(1,i,j),1,i,j) )
     &       + fv*fr_snow_ij(2,i,j)
     &           * sum( dzsn_ij(1:nsn_ij(2,i,j),2,i,j) ) )
            atmlnd%gtemp(i,j)=tsns_ij(i,j)
            atmlnd%gtempr(i,j) =tearth(i,j)+tf
#ifdef SCM
            if( SCMopt%Tskin )then
              atmlnd%gtemp(i,j) = SCMin%Tskin - tf
              atmlnd%gtempr(i,j) = SCMin%Tskin
            endif
#endif
          end if
        end do
      end do

      call set_roughness_length

#ifdef TRACERS_WATER
ccc still not quite correct (assumes fw=1)
      do j=J_0,J_1
        do i=I_0,I_1
          atmlnd%gtracer(:,i,j)=0.  ! default
          !if (fearth(i,j).le.0.d0) cycle
          if (focean(i,j) .ge. 1.d0) cycle
          !fb=afb(i,j) ; fv=1.-fb
          call get_fb_fv( fb, fv, i, j )
      call ent_get_exports( entcells(i,j),
     &     canopy_height=height_can
     &     )
          fm=1.d0-exp(-snowbv(2,i,j)/((height_can*spgsn) + 1d-12))
          if ( fm < 1.d-3 ) fm=0.d0

          call compute_gtracer( ntm, fb, fv, fm, w_ij(:,:,i,j),
     &         fr_snow_ij(:,i,j), wsn_ij(:,:,i,j),
     &         tr_w_ij(:,:,:,i,j), tr_wsn_ij(:,:,:,i,j),
     &         atmlnd%gtracer(:,i,j) )

cddd          wsoil_tot=fb*( w_ij(1,1,i,j)*(1.d0-fr_snow_ij(1,i,j))
cddd     &     + wsn_ij(1,1,i,j)*fr_snow_ij(1,i,j) )
cddd     &     + fv*( w_ij(0,2,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))   !*1.d0
cddd     &     + wsn_ij(1,2,i,j)*fm*fr_snow_ij(2,i,j) )
cddd          do n=1,ntm
cddd            if (itime_tr0(n).gt.itime) cycle
cddd            if ( .not. needtrs(n) ) cycle
cddd            ! should also restrict to TYPE=nWATER ?
cddd            if ( wsoil_tot > 1.d-30 ) then
cddd            gtracer(n,i,j) = (
cddd     &           fb*( tr_w_ij(n,1,1,i,j)*(1.d0-fr_snow_ij(1,i,j))
cddd     &           + tr_wsn_ij(n,1,1,i,j) )         !*fr_snow_ij(1,i,j)
cddd     &           + fv*( tr_w_ij(n,0,2,i,j)*(1.d0-fm*fr_snow_ij(2,i,j))
cddd     &           + tr_wsn_ij(n,1,2,i,j)*fm ) )    !*fr_snow_ij(2,i,j)
cddd     &           /(rhow*wsoil_tot)
cddd            else
cddd              gtracer(n,i,j) = 0.
cddd            end if
cddd          enddo
        end do
      end do
#endif

      end subroutine init_land_surface

!********************************************************************

      subroutine old_gic_2_modele(
     &           w, ht,snowbv,
     &           wearth, aiearth, tearth, snowe )
      real*8 :: w(:,:), ht(:,:),snowbv(:),
     &           wearth, aiearth, tearth, snowe

      call stop_model(
     &     "conversion old_gic_2_modele not supported yet",255)

      end subroutine old_gic_2_modele


      subroutine fix_water_ic( w, q, dz, i, j )
      use model_com, only : qcheck
      use ghy_com, only : ngm
      use sle001, only : get_soil_properties
      !-- out
      real*8, intent(inout) :: w(0:,:)
      !-- in
      real*8, intent(in) :: q(:,:), dz(:)
      integer, intent(in) :: i, j
      !--- local
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2)
      real*8 wmax,wmin
      integer ibv, k

c outer loop over ibv
      do ibv=1,2

        call get_soil_properties( q, dz,
     &     thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )

c initialize soil (w, ht) from earth_*
        do k=1,ngm
          wmax = dz(k)*thets(k,ibv)
          wmin = dz(k)*thetm(k,ibv)
          if ( w(k,ibv) > wmax ) then
            if (qcheck) print *,"fix_water_ic: reducing water, cell ",
     *            i, j, "layer ", k, ibv, w(k,ibv), " -> ", wmax
            if (qcheck) print *, "q= ", q(:,k)
            w(k,ibv) = wmax
          endif
          if ( w(k,ibv) < wmin ) then
            if (qcheck) print *,"fix_water_ic: increasing water, cell ",
     &            i, j, "layer ", k, ibv,
     &           w(k,ibv), " -> ", wmin
            if (qcheck) print *, "q= ", q(:,k)
            w(k,ibv) = wmin
          endif
        enddo
      enddo

      end subroutine fix_water_ic


      subroutine fix_heat_ic(
     &     w, ht,
     &     shc_can,
     &     q, dz, i, j )

      use constant, only : rhow
     &     ,shi_kg=>shi,lhm
      use snow_model, only: snow_redistr, snow_fraction
      use ghy_com, only : ngm
      use sle001, only : get_soil_properties
      !-- out
      real*8, intent(inout) :: ht(0:,:)
      !-- in
      real*8, intent(in) :: w(0:,:)
      real*8, intent(in) :: shc_can
      real*8, intent(in) :: q(:,:), dz(:)
      integer i, j
      !--- local
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2), temp, fice
      real*8 :: t_max=60.d0, t_min=-80.d0
      integer ibv, k

c for canopy:
      call heat_to_temperature(temp,
     &     fice, ht(0,2), w(0,2), shc_can)
      if ( temp > t_max .or. temp < t_min ) then
        print *,"fix_heat_ic: correcting canopy temp", temp,
     &       "cell= ", i,j
        temp = max( temp, t_min )
        temp = min( temp, t_max )
        call temperature_to_heat( ht(0,2),
     &       temp, fice, w(0,2), shc_can )
      endif

c outer loop over ibv
      do ibv=1,2

        call get_soil_properties( q, dz,
     &     thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )

        do k=1,ngm

          call heat_to_temperature(temp,
     &         fice, ht(k,ibv), w(k,ibv), shc(k,ibv))
          if ( temp > t_max .or. temp < t_min ) then
            print *,"fix_heat_ic: correcting soil temp", temp,
     &           "cell= ", i,j, "layer= ", ibv, k
            temp = max( temp, t_min )
            temp = min( temp, t_max )
            call temperature_to_heat( ht(k,ibv),
     &           temp, fice, w(k,ibv), shc(k,ibv) )
          endif

        enddo
      enddo

      end subroutine fix_heat_ic


      subroutine tp_sat_2_ht_w(
     &     w, ht,
     &     nsn, dzsn, wsn, hsn, fr_snow,

     &     earth_sat, earth_ice,
     &     earth_tp, snowd,
     &     tsnowtop,

     &     ws_can, shc_can,
     &     q, dz )

      use constant, only : rhow
     &     ,shi_kg=>shi,lhm
      use snow_model, only: snow_redistr, snow_fraction
      use ghy_com, only : ngm
      use sle001, only : get_soil_properties
      !-- out
      real*8, intent(out) :: w(0:,:), ht(0:,:)
      real*8, intent(out) :: dzsn(:,:), wsn(:,:), hsn(:,:), fr_snow(:)
      integer, intent(out) :: nsn(:)
      !-- in
      real*8, intent(in) :: earth_sat(0:,:), earth_ice(0:,:),
     &      earth_tp(0:,:)
      real*8, intent(in) ::  snowd(:), tsnowtop
      real*8, intent(in) :: ws_can, shc_can
      real*8, intent(in) :: q(:,:), dz(:)
      !--- local
!@var shi heat capacity of pure ice (J/m^3 C)
      real*8, parameter :: shi= shi_kg * rhow
!@var fsn latent heat of melt (J/m^3)
      real*8, parameter :: fsn= lhm * rhow
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2), tsn1(2)
      integer ibv, k

c for canopy:
      w(0,2) = ws_can*earth_sat(0,2)
      call temperature_to_heat( ht(0,2),
     &     earth_tp(0,2), earth_ice(0,2), w(0,2), shc_can )

c outer loop over ibv
      do ibv=1,2

        call get_soil_properties( q, dz,
     &     thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )

c initialize soil (w, ht) from earth_*
        do k=1,ngm
          w(k,ibv) = thets(k,ibv)*dz(k)*earth_sat(k,ibv)
          call temperature_to_heat(ht(k,ibv),
     &         earth_tp(k,ibv), earth_ice(k,ibv), w(k,ibv), shc(k,ibv) )
        enddo

        call reset_snow_to_zero

        if ( snowd(ibv) <= 0.d0 ) cycle

c if there is a snow put it all in the first layer (assume rho_snow = 200)
        dzsn(1,ibv)=snowd(ibv) * 5.d0
        wsn(1,ibv)=snowd(ibv)
        fr_snow(ibv) = 1.d0

c set snow temperature to average of snow top and first soil layer (or 0)
        tsn1(ibv) = min( .5d0*(tsnowtop+earth_tp(1,ibv)), 0.d0 )

c use snow temperature to get the heat of the snow
        call temperature_to_heat(hsn(1,ibv),
     &         tsn1(ibv), 1.d0, wsn(1,ibv), 0.d0)

          ! hack to get identical results ??? (need this ??)
        hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn

        call snow_fraction(dzsn(:,ibv), nsn(ibv), 0.d0, 0.d0,
     &       1.d0, fr_snow(ibv) )

        if ( fr_snow(ibv) == 0.d0 ) then
          call reset_snow_to_zero
          cycle
        endif

        call snow_redistr(dzsn(:,ibv), wsn(:,ibv), hsn(:,ibv),
     &       nsn(ibv), 1.d0/fr_snow(ibv) )

!!! debug : check if  snow is redistributed correctly
        if ( dzsn(1,ibv) > 0.d0 .and. dzsn(1,ibv) < .099d0) then
          call stop_model("set_snow: error in dz",255)
        endif

        if ( dzsn(1,ibv) > 0.d0 ) then
          if (  wsn(1,ibv)/dzsn(1,ibv)  < .1d0)
     &         call stop_model("set_snow: error in dz",255)
        endif

      enddo  ! ibv

      return

      contains

      subroutine reset_snow_to_zero
      ! set one empty layer of som
      nsn(ibv)=1
      dzsn(1,ibv)=0.d0
      wsn(1,ibv)=0.d0
      hsn(1,ibv)=0.d0
      tsn1(ibv)=0.d0
      fr_snow(ibv) = 0.d0
      end subroutine reset_snow_to_zero

      end subroutine tp_sat_2_ht_w


      subroutine set_snow1( w, ht, snowd, q, dz,
     &     nsn, dzsn, wsn, hsn, fr_snow )
!@sum set_snow extracts snow from the first soil layer and initializes
!@+   snow model prognostic variables
!@+   should be called when model restarts from the old restart file
!@+   ( which doesn''t contain new snow model (i.e. 3 layer) data )
c
c input:
c snowd(2) - landsurface snow depth
c w(k,2)   - landsurface water in soil layers
c ht(k,2)  - landsurface heat in soil layers
c fsn      - heat of fusion
c shi      - specific heat of ice
c shc(k,2) - heat capacity of soil layers
c
c output:
c dzsn(lsn,2) - snow layer thicknesses
c wsn(lsn,2)  - snow layer water equivalent depths
c hsn(lsn,2)  - snow layer heat contents
c tsn1(2)     - snow top temperature
c isn(2)      - 0 if no snow, 1 if snow
c nsn(2)      - number of snow layers
c snowd(2)
c w(k,2)
c ht(k,2)
c
c calling sequence:
c
c     assignment of w,ht,snowd
c     call ghinij(i,j,wfc1)
c     call set_snow
c note: only to be called when initializing from landsurface
c       prognostic variables without the snow model.
c
      use constant, only : rhow
     &     ,shi_kg=>shi,lhm
      use snow_model, only: snow_redistr, snow_fraction
      use ghy_com, only : ngm
      use sle001, only : get_soil_properties
      real*8, intent(inout) :: w(:,:), ht(:,:), snowd(:), q(:,:), dz(:)
      integer, intent(out) :: nsn(:)
      real*8, intent(out) :: dzsn(:,:), wsn(:,:), hsn(:,:), fr_snow(:)
      !---
!@var shi heat capacity of pure ice (J/m^3 C)
      real*8, parameter :: shi= shi_kg * rhow
!@var fsn latent heat of melt (J/m^3)
      real*8, parameter :: fsn= lhm * rhow
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2), tsn1(2), fice(2)


      integer ibv

c outer loop over ibv
      do ibv=1,2

        call get_soil_properties( q, dz,
     &     thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )

c initalize all cases to nsn=1
        nsn(ibv)=1


ccc since we don''t know what kind of data we are dealing with,
ccc better check it

        if( snowd(ibv) .gt. w(1,ibv)-dz(1)*thetm(1,ibv)  ) then
          write(99,*) 'snowd corrected: old=', snowd(ibv)
          snowd(ibv) = w(1,ibv)-dz(1)*thetm(1,ibv) - 1.d-10
          write(99,*) '                 new=', snowd(ibv)
          if ( snowd(ibv) .lt. -0.01d0 )
     &         call stop_model('set_snow: neg. snow',255)
          if ( snowd(ibv) .lt. 0.d0 ) snowd(ibv) = 0.d0 ! rounding error
        endif


c if there is no snow, set isn=0.  set snow variables to 0.d0
        if(snowd(ibv).le.0.d0)then
          dzsn(1,ibv)=0.d0
          wsn(1,ibv)=0.d0
          hsn(1,ibv)=0.d0
          tsn1(ibv)=0.d0
          fr_snow(ibv) = 0.d0
        else


c given snow, set isn=1.d0
c!!!        dzsn(1,ibv)=snowd(ibv)/spgsn
c!!!  replacing prev line considering rho_snow = 200
          dzsn(1,ibv)=snowd(ibv) * 5.d0
          wsn(1,ibv)=snowd(ibv)
c!!! actually have to compute fr_snow and modify dzsn ...
          fr_snow(ibv) = 1.d0

c given snow, temperature of first layer can''t be positive.
c the top snow temperature is the temperatre of the first layer.
cddd          if(fsn*w(1,ibv)+ht(1,ibv).lt.0.d0)then
cddd            tsn1(ibv)=(ht(1,ibv)+w(1,ibv)*fsn)/(shc(1,ibv)+w(1,ibv)*shi)
cddd          else
cddd            tsn1(ibv)=0.d0
cddd          endif
cddd          print *,"tsn = ", tsn1(ibv)
          call heat_to_temperature(tsn1(ibv),
     &         fice(ibv), ht(1,ibv), w(1,ibv), shc(1,ibv))
cddd          print *,"      ", tsn1(ibv)
          tsn1(ibv) = min( tsn1(ibv), 0.d0 )

c use snow temperature to get the heat of the snow
cddd          hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
cddd          print *,"hsn = ", hsn(1,ibv)
          call temperature_to_heat(hsn(1,ibv),
     &         tsn1(ibv), 1.d0, wsn(1,ibv), 0.d0)
cddd          print *,"      ", hsn(1,ibv)

          ! hack to get identical results
          hsn(1,ibv)=tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn

c subtract the snow from the landsurface prognositic variables
          w(1,ibv)=w(1,ibv)-wsn(1,ibv)
          ht(1,ibv)=ht(1,ibv)-hsn(1,ibv)

ccc and now limit all the snow to 5cm water equivalent
          if ( snowd(ibv) .gt. 0.05d0 ) then
            snowd(ibv) = 0.05d0
            dzsn(1,ibv)= snowd(ibv) * 5.d0
            wsn(1,ibv)= snowd(ibv)
            hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
          endif

ccc and now limit all the snow to 50cm water equivalent (for debug)
cddd          if ( snowd(ibv) .gt. 0.0005d0 ) then
cddd            snowd(ibv) = 0.5d0
cddd            dzsn(1,ibv)= snowd(ibv) * 5.d0
cddd            wsn(1,ibv)= snowd(ibv)
cddd            hsn(1,ibv)= tsn1(ibv)*wsn(1,ibv)*shi-wsn(1,ibv)*fsn
cddd          endif

ccc redistribute snow over the layers and recompute fr_snow
ccc (to make the data compatible with snow model)
          if ( .not. ( hsn(1,ibv) > 0. .or.  hsn(1,ibv) <= 0. ) )
     &        call stop_model("ERR in init_snow: NaN", 255)

          call snow_fraction(dzsn(:,ibv), nsn(ibv), 0.d0, 0.d0,
     &         1.d0, fr_snow(ibv) )
          if ( fr_snow(ibv) > 0.d0 ) then
            call snow_redistr(dzsn(:,ibv), wsn(:,ibv), hsn(:,ibv),
     &           nsn(ibv), 1.d0/fr_snow(ibv) )
            if ( .not. ( hsn(1,ibv) > 0. .or.  hsn(1,ibv) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)
          else
            snowd(ibv) = 0.d0
            dzsn(1,ibv)=0.d0
            wsn(1,ibv)=0.d0
            hsn(1,ibv)=0.d0
            tsn1(ibv)=0.d0
            fr_snow(ibv) = 0.d0
          endif

        endif
!!! debug : check if  snow is redistributed correctly
        if ( dzsn(1,ibv) > 0.d0 .and. dzsn(1,ibv) < .099d0) then
          call stop_model("set_snow: error in dz",255)
        endif

        if ( dzsn(1,ibv) > 0.d0 ) then
          if (  wsn(1,ibv)/dzsn(1,ibv)  < .1d0)
     &         call stop_model("set_snow: error in dz",255)
        endif


      enddo

            if ( .not. ( hsn(1,1) > 0. .or.  hsn(1,1) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)
            if ( .not. ( hsn(1,2) > 0. .or.  hsn(1,2) <= 0. ) )
     &        call stop_model("ERR in init_snow 2: NaN", 255)



      return
      end subroutine set_snow1



      subroutine reset_gh_to_defaults( reset_prognostic )
      !use model_com, only: vdata
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      use fluxes, only : atmlnd
      use ghy_com
      logical, intent(in) :: reset_prognostic
      integer i,j

C**** define local grid
      integer I_0, I_1, J_0, J_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

ccc ugly, should fix later
      call stop_model(
     &     "reset_gh_to_defaults not implemented yet for Ent",255)
      ! worked without Ent ...
!#else , i.e. #ifndef uSE_ENT
!      call reset_veg_to_defaults( reset_prognostic )
!#endif

      do j=J_0,J_1
      do i=I_0,I_1

      dz_ij(i,j,1:ngm)= (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      q_ij(i,j,1:imt,1:ngm)=
     &     reshape( (/  0.33491623d+00,  0.52958947d+00,
     &     0.13549370d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.32995611d+00,  0.52192056d+00,  0.14812243d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.32145596d+00,
     &     0.48299056d+00,  0.19555295d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.47638881d+00,  0.40400982d+00,
     &     0.11959970d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.99985123d-01,  0.95771909d-01,  0.41175738d-01,
     &     0.00000000d+00,  0.76306665d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.10000000d+01 /), (/imt,ngm/) )
      qk_ij(i,j,1:imt,1:ngm)=
     &     reshape( (/  0.34238762d+00,  0.52882469d+00,
     &     0.12878728d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.32943058d+00,  0.52857041d+00,  0.14199871d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.30698991d+00,
     &     0.52528000d+00,  0.16772974d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.39890009d+00,  0.43742162d+00,
     &     0.16367787d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.46536058d+00,  0.39922065d+00,  0.13541836d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.00000000d+00,  0.00000000d+00,  0.00000000d+00,
     &     0.10000000d+01 /), (/imt,ngm/) )
      sl_ij(i,j)= 0.22695422d+00
      top_index_ij(i,j)= 0.10832934d+02
      top_dev_ij(i,j)= 0.21665636d+03

      if ( .not. reset_prognostic ) cycle

      atmlnd%snowe(i,j)= 0.65458111d-01
      tearth(i,j)= -0.12476520d+00
      tsns_ij(i,j)= -0.12476520d+00
      wearth(i,j)=  0.29203081d+02
      aiearth(i,j)=  0.93720329d-01
      w_ij(:,1,i,j) = (/0.d0,  0.17837750d-01,  0.40924843d-01,
     &     0.77932012d-01,  0.11919649d+00,  0.57237469d-01,
     &     0.10000000d-11 /)
      w_ij(:,2,i,j) = (/  0.10000000d-11,  0.29362259d-01,
     &     0.50065177d-01,  0.82533140d-01,  0.10383620d+00,
     &     0.31552459d-01,  0.10000000d-11 /)
      ht_ij(:,1,i,j)= (/  0.00000000d+00, -0.15487181d+07,
     &     -0.50720067d+07,  0.18917623d+07,  0.77174974d+07,
     &     0.21716040d+08,  0.44723067d+08 /)
      ht_ij(:,2,i,j)= (/ -0.13991376d+05, -0.53165599d+05,
     &     0.65443775d+06,  0.29276050d+07,  0.81455096d+07,
     &     0.21575081d+08,  0.45952255d+08 /)
      snowbv(1:2,i,j)= (/  0.00000000d+00,  0.65458111d-04 /)

      enddo
      enddo

      end subroutine reset_gh_to_defaults



      subroutine checke(subr)
!@sum  checke checks whether arrays are reasonable over earth
!@auth original development team
      use model_com, only : itime
      use geom, only : imaxj
      use constant, only : rhow
      !use veg_com, only : afb
      use fluxes, only : atmlnd
      use ghy_com, only : tearth,wearth,aiearth,w_ij,ht_ij
     *     ,snowbv,ngm,fearth,wsn_ij,fr_snow_ij,nsn_ij,LS_NFRAC,wfcs
#ifdef TRACERS_WATER
     &     ,tr_w_ij,tr_wsn_ij
      use OldTracer_mod, only: trname, t_qlimit
      USE TRACER_COM, only : NTM
#endif
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      implicit none
!@var subr identifies where check was called from
      character*6, intent(in) :: subr

      real*8 x,tgl,wtrl,acel
      integer i,j,imax,jmax,nsb,nsv,n
      real*8, parameter :: EPS=1.d-12
      logical QCHECKL
      real*8 relerr, errmax, fb, fv

C**** define local grid
      integer I_0, I_1, J_0, J_1, njpol

      QCHECKL = .false.
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     *     J_STRT=J_0, J_STOP=J_1)
      njpol = grid%J_STRT_SKP-grid%J_STRT

c**** check for nan/inf in earth data
      call check3c(w_ij(1:ngm,1,I_0:I_1,J_0:J_1) ,ngm  ,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'wb')
      call check3c(w_ij(0:ngm,2,I_0:I_1,J_0:J_1) ,ngm+1,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'wv')
      call check3c(ht_ij(0:ngm,1,I_0:I_1,J_0:J_1),ngm+1,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'hb')
      call check3c(ht_ij(0:ngm,2,I_0:I_1,J_0:J_1),ngm+1,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'hv')
      call check3c(snowbv(1:LS_nfrac,I_0:I_1,J_0:J_1),2    ,
     *                 I_0,I_1,J_0,J_1,njpol,subr,'sn')

c**** check for reasonable temperatures over earth
      x=1.001
      do j=J_0,J_1
        do i=I_0,imaxj(j)
          if (fearth(i,j).gt.0.) then
            tgl=tearth(i,j)
            wtrl=wearth(i,j)
            acel=aiearth(i,j)
            if ((tgl+60.)*(60.-tgl).le.0.) write (6,901) subr,i,j,itime
     *           ,fearth(i,j),'tg1 off',atmlnd%snowe(i,j),tgl,wtrl,acel
            if (wtrl.lt.0..or.acel.lt.0..or.(wtrl+acel).gt.x*wfcs(i
     *           ,j)) write(6,901) subr,i,j,itime,fearth(i,j),'wtr off'
     *           ,atmlnd%snowe(i,j),tgl,wtrl,acel,wfcs(i,j)
          end if
        end do
      end do

      !print *,"checke: w(51,33) ", w_ij(:,:,51,33)
      !print *,"checke: fearth(51,33) ", fearth(51,33),afb(51,33)


c**** check tracers
#ifdef TRACERS_WATER
      do n=1,ntm
        ! check for neg tracers
        if (t_qlimit(n)) then
          do j=J_0, J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              if ( minval( tr_w_ij(n,:,:,i,j)   ) < -1.d15 .or.
     &             minval( tr_wsn_ij(n,:,:,i,j) ) < -1.d15 ) then
                print*,"Neg tracer in earth after ",SUBR,i,j,trname(n)
     &               , "tr_soil= ", tr_w_ij(n,:,:,i,j)
     &               , "TR_SNOW= ", tr_wsn_ij(n,:,:,i,j)
                QCHECKL=.TRUE.
              end if

            end do
          end do
        end if

        ! check if water == water
        if (trname(n) == 'Water') then
          errmax = 0. ; imax=I_0 ; jmax=J_0
          do j=J_0, J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              !fb = afb(i,j)
              !fv = 1.d0 - fb
              call get_fb_fv( fb, fv, i, j )
              relerr= ( fb*sum(abs(
     &             tr_w_ij(n,1:ngm,1,i,j) - w_ij(1:ngm,1,i,j)*rhow))
     &             + fv*sum(abs(
     &             tr_w_ij(n,0:ngm,2,i,j) - w_ij(0:ngm,2,i,j)*rhow)) )
     &             /(  rhow*( fb*sum(w_ij(1:ngm,1,i,j))
     &             + fv*sum(w_ij(0:ngm,2,i,j)) ) + EPS  )
              if (relerr > errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            enddo
          enddo
          !fb = afb(imax,jmax)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, imax, jmax )
          print*,"Relative error in soil after ",trim(subr),":"
     &         ,imax,jmax,errmax
     &         ,( tr_w_ij(n,1:ngm,1,imax,jmax)
     &         - w_ij(1:ngm,1,imax,jmax)*rhow )*fb
     &         ,( tr_w_ij(n,0:ngm,2,imax,jmax)
     &         - w_ij(0:ngm,2,imax,jmax)*rhow )*fv
     &         , rhow*( fb*sum(w_ij(1:ngm,1,imax,jmax))
     &         +        fv*sum(w_ij(0:ngm,2,imax,jmax)) )
cddd     &         ,w_ij(1:ngm,1,imax,jmax)*rhow
cddd     &         ,w_ij(0:ngm,2,imax,jmax)*rhow

          errmax = 0. ; imax=I_0 ; jmax=J_0
          do j=J_0, J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              !fb = afb(i,j)
              !fv = 1.d0 - fb
              call get_fb_fv( fb, fv, i, j )
              nsb = nsn_ij(1,i,j)
              nsv = nsn_ij(2,i,j)
              relerr= ( fb*sum(abs( tr_wsn_ij(n,1:nsb,1,i,j)
     &             - wsn_ij(1:nsb,1,i,j)*rhow*fr_snow_ij(1,i,j) ))
     &             + fv*sum(abs( tr_wsn_ij(n,1:nsv,2,i,j)
     &             - wsn_ij(1:nsv,2,i,j)*rhow*fr_snow_ij(2,i,j))) )
     &             /(rhow*(
     &             fb*sum(wsn_ij(1:nsb,1,i,j))*fr_snow_ij(1,i,j) +
     &             fv*sum(wsn_ij(1:nsv,2,i,j))*fr_snow_ij(2,i,j) ) +EPS)
              if (relerr > errmax) then
                imax=i ; jmax=j ; errmax=relerr
              end if
            enddo
          enddo
          !fb = afb(imax,jmax)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, imax, jmax )
          nsb = nsn_ij(1,imax,jmax)
          nsv = nsn_ij(2,imax,jmax)
          print*,"Relative error in snow after ",trim(subr),":"
     &         ,imax,jmax,errmax
     &         ,( tr_wsn_ij(n,1:nsb,1,imax,jmax)
     &         - wsn_ij(1:nsb,1,imax,jmax)*rhow
     &         *fr_snow_ij(1,imax,jmax) )*fb
     &         ,( tr_wsn_ij(n,1:nsv,2,imax,jmax)
     &         - wsn_ij(1:nsv,2,imax,jmax)*rhow
     &         *fr_snow_ij(2,imax,jmax) )*fv
     &         ,rhow*( fb*sum(wsn_ij(1:nsb,1,imax,jmax))
     &         * fr_snow_ij(1,imax,jmax)
     &         +       fv*sum(wsn_ij(1:nsv,2,imax,jmax))
     &         * fr_snow_ij(2,imax,jmax) )

cddd     &         ,wsn_ij(1:nsn_ij(1,imax,jmax),1,imax,jmax)*rhow
cddd     &         *fr_snow_ij(1,imax,jmax)
cddd     &         ,wsn_ij(1:nsn_ij(1,imax,jmax),2,imax,jmax)*rhow
cddd     &         *fr_snow_ij(2,imax,jmax)



        endif
      enddo
#endif

      IF (QCHECKL)
     &     call stop_model('CHECKL: Earth variables out of bounds',255)

      return
 901  format ('0subr,i,j,i-time,pearth,',a7,2i4,i10,f5.2,1x
     *     ,a/' snw,tg1,wtr1,ice1 (wfc1)',5f12.4)

      end subroutine checke

      subroutine daily_earth(end_of_day)
!@sum  daily_earth performs daily tasks for earth related functions
!@auth original development team
!@calls RDLAI
      use constant, only : rhow,twopi,tf
      use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
      use model_com, only : modelEclock
      use model_com, only : nday
      use fluxes, only : nisurf,focean,atmlnd
      use geom, only : imaxj,lat2d
      use diag_com, only : aij=>aij_loc
     *     ,tdiurn,ij_strngts,ij_dtgdts,ij_tmaxe,ij_tmaxc
     *     ,ij_tdsl,ij_tmnmx,ij_tdcomp, ij_dleaf
      use ghy_com, only : fearth, wsn_max,
     &     q_ij,dz_ij,ngm,w_ij,wfcs
     &     ,aalbveg
      use surf_albedo, only: albvnh, updsur  !nyk
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      !use sle001, only : fb,fv,ws
      use sle001, only : get_soil_properties
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
      use ent_drv, only : update_vegetation_data
      !!use ent_com, only : entcells
      !!use ent_mod, only : ent_get_exports

      implicit none
      real*8 tsavg,wfc1
      real*8 aleafmass, aalbveg0, fvp, sfv  !nyk veg ! , aleafmasslast
      integer i,j,itype
      integer northsouth,iv  !nyk
      logical, intent(in) :: end_of_day
      real*8 fb, fv, ws_can
      real*8 thetm(ngm,2), thets(ngm,2), shc(ngm,2)
      integer ibv
!      integer hemi(1:IM,grid%J_STRT:grid%J_STOP)
!      integer :: JEQUATOR=JM/2
C**** define local grid
      integer I_0, I_1, J_0, J_1
      real*8 ws11,ws12
      integer :: year, dayOfYear

      call modelEclock%get(year=year, dayOfYear=dayOfYear)

C**** Extract useful local domain parameters from "grid"
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if (end_of_day) call accumulate_excess_C(0)

      ! update water and heat in land surface fractions if
      ! lake fraction changed
      if (end_of_day .and. variable_lk > 0) call update_land_fractions

      if (end_of_day .and. wsn_max>0) call remove_extra_snow_to_ocean

!!! testing
!!!      aalbveg(:,:) = 0.08D0
!!!      return

! the following is probably not needed since we update vegetation
! on every time step
 !     if (jday==1) then ! i guess we need to call it only once per year
 !       hemi(:,JEQUATOR+1:J_1) = 1
 !       hemi(:,J_0:JEQUATOR) = -1
 !       call ent_prescribe_vegupdate(entcells,hemi,jday,jyear, .true.)
 !     endif

      if (end_of_day )
     &     call update_vegetation_data( entcells,
     &     im, jm, I_0, I_1, J_0, J_1, dayOfYear, year )
      !if(cond_scheme.eq.2) call updsur (0,jday)
      ! we don''t use cond_scheme==1 any more, so call it always
      call updsur (0,dayOfYear)
c****
      call set_roughness_length

            !albvnh(9,6,2)=albvnh(1+8veg,6bands,2hemi), band 1 is VIS.
      atmlnd%bare_soil_wetness(:,:) = 0.d0 ! make sure that it is initialized
      do j=J_0,J_1
        do i=I_0,I_1
          if(lat2d(i,j).lt.0.) then !nyk added northsouth
            northsouth=1            !southern hemisphere
          else
            northsouth=2            !northern hemisphere
          end if
          wfcs(i,j)=24.
          ! if (fearth(i,j).gt.0.) then
          !if (focean(i,j) < 1.d0 ) then
          if (variable_lk==0) then
            if ( fearth(i,j) <= 0.d0 ) cycle
          else
            if ( focean(i,j) >= 1.d0 ) cycle
          endif

!!!            call ghinij(i,j)
            do ibv=1,2
              call get_soil_properties( q_ij(i,j,:,:), dz_ij(i,j,:),
     &             thets(1:,ibv), thetm(1:,ibv), shc(1:,ibv) )
            enddo

            call get_fb_fv( fb, fv, i, j )

            if ( fv > 0.d0 ) then
              call ent_get_exports(entcells(i,j),canopy_max_H2O=ws_can)
            else
              ws_can = 0.d0
            endif
!!!            wfc1=fb*ws(1,1)+fv*(ws_can+ws(1,2))
cddd            wfc1=fb*thets(1,1)*dz_ij(i,j,1) +
cddd     &           fv*( ws_can + thets(1,2)*dz_ij(i,j,1) )
            ws11 = thets(1,1)*dz_ij(i,j,1)
            ws12 = thets(1,2)*dz_ij(i,j,1)
            wfc1=fb*ws11 +
     &           fv*( ws_can + ws12 )
            wfcs(i,j)=rhow*wfc1 ! canopy part changes
cddd            write(934,*) "wfcs", i,j,wfcs(i,j)

            atmlnd%bare_soil_wetness(i,j) =
     &           w_ij(1,1,i,j) / ( thets(1,1)*dz_ij(i,j,1) )

         !!! this diag belongs to Ent - commenting out
          !end if
        end do
      end do

      if (end_of_day) then
        do j=J_0,J_1
        do i=I_0,imaxj(j)
          tsavg=tdiurn(i,j,5)/(nday*nisurf)
          if(32.+1.8*tsavg.lt.65.)
     *         aij(i,j,ij_strngts)=aij(i,j,ij_strngts)+(33.-1.8*tsavg)
          aij(i,j,ij_dtgdts)=aij(i,j,ij_dtgdts)+18.*((tdiurn(i,j,2)-
     *         tdiurn(i,j,1))/(tdiurn(i,j,4)-tdiurn(i,j,3)+1.d-20)-1.)
          aij(i,j,ij_tdsl)=aij(i,j,ij_tdsl)+
     *         (tdiurn(i,j,4)-tdiurn(i,j,3))*fearth(i,j)
          aij(i,j,ij_tdcomp)=aij(i,j,ij_tdcomp)+
     *         (tdiurn(i,j,6)-tdiurn(i,j,9))
          aij(i,j,ij_tmaxe)=aij(i,j,ij_tmaxe)+
     *         (tdiurn(i,j,4)-tf)*fearth(i,j)
          aij(i,j,ij_tmaxc)=aij(i,j,ij_tmaxc) + (tdiurn(i,j,6)-tf)
          if (tdiurn(i,j,6)-tf.lt.aij(i,j,ij_tmnmx))
     *         aij(i,j,ij_tmnmx)=tdiurn(i,j,6)-tf
        end do
        end do
      end if


#ifdef TRACERS_DRYDEP
      CALL RDLAI ! read leaf area indices for tracer dry deposition
#endif

      if (end_of_day) call accumulate_excess_C(1)

      return
      end subroutine daily_earth

      subroutine ground_e
!@sum  ground_e driver for applying surface fluxes to land fraction
!@auth original development team
      use geom, only : imaxj
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      use ghy_com, only : tearth,wearth,aiearth,w_ij
     *     ,snowbv,gdeep,fearth
      !use veg_com, only : afb
      use diag_com, only : itearth,aij=>aij_loc
     *     ,jreg,j_wtr1,j_ace1,j_wtr2,j_ace2
     *     ,ij_f0e,ij_evape,ij_gwtr
     *     ,ij_gice, ij_gwtr1
     *     ,ij_gbssnd,ij_gvssnd,ij_gbsw,ij_gvsw
      use fluxes, only : atmlnd,eprec
      implicit none
!@var ae0_ep total heat from atm. to land surface withoot eprec (J/m^2)
      real*8 snow,ae0_ep,evap,wtr1,wtr2,ace1,ace2
     *     ,pearth,enrgp,fb,fv
      integer i,j,jr,k

C**** define local grid
      integer J_0, J_1, J_0H, J_1H, I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0      , J_STOP=J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
      do i=I_0,imaxj(j)
      pearth=fearth(i,j)
      jr=jreg(i,j)
      if (pearth.gt.0) then

        snow=atmlnd%snow(i,j)!snowe(i,j)
        !tg1 = tearth(i,j)
        wtr1= wearth(i,j)
        ace1=aiearth(i,j)
        !tg2=gdeep(i,j,1)
        wtr2=gdeep(i,j,2)
        ace2=gdeep(i,j,3)
        ae0_ep=atmlnd%e0(i,j)
        evap=atmlnd%evapor(i,j)
        enrgp=eprec(i,j)      ! including latent heat
        call get_fb_fv( fb, fv, i, j )

c**** accumulate diagnostics
c**** the following is the actual snow cover of the snow model
c        scove = pearth *
c     *       ( afb(i,j)*fr_snow_ij(1,i,j)
c     *       + (1.-afb(i,j))*fr_snow_ij(2,i,j) )
c**** the following computes the snow cover as it is used in RAD_DRV.f
c        scove = pearth * atmlnd%snowfr(i,j)
c     *       ( fb*atmlnd%fr_snow_rad(1,i,j)
c     *       + fv*atmlnd%fr_snow_rad(2,i,j) )

        !if (snowe(i,j).gt.0.) scove=pearth

        call inc_aj(i,j,itearth,j_wtr1,wtr1*pearth)
        call inc_aj(i,j,itearth,j_ace1,ace1*pearth)
        call inc_aj(i,j,itearth,j_wtr2,wtr2*pearth)
        call inc_aj(i,j,itearth,j_ace2,ace2*pearth)

        call inc_areg(i,j,jr,j_wtr1,wtr1*pearth)
        call inc_areg(i,j,jr,j_ace1,ace1*pearth)
        call inc_areg(i,j,jr,j_wtr2,wtr2*pearth)
        call inc_areg(i,j,jr,j_ace2,ace2*pearth)

        aij(i,j,ij_f0e)  =aij(i,j,ij_f0e)  +(ae0_ep+enrgp)*pearth
        aij(i,j,ij_gwtr) =aij(i,j,ij_gwtr)+(wtr1+ace1+wtr2+ace2)*pearth
        aij(i,j,ij_gwtr1) =aij(i,j,ij_gwtr1)+(wtr1+ace1)*pearth
        aij(i,j,ij_gice) =aij(i,j,ij_gice)+(ace1+ace2)*pearth
        aij(i,j,ij_evape)=aij(i,j,ij_evape)+evap*pearth
!     Water in vegetated layers 0 (canopy) - 6 and
!     bare soil layers 1-6
        do k=1,6
          aij(i,j,ij_gbsw+k-1)=aij(i,j,ij_gbsw+k-1)+w_ij(k,1,i,j)
     &         *pearth*fb
        end do
        do k=0,6
          aij(i,j,ij_gvsw+k)=aij(i,j,ij_gvsw+k)+w_ij(k,2,i,j)
     &         *pearth*fv
        end do

        aij(i,j,ij_gbssnd)=aij(i,j,ij_gbssnd)+snowbv(1,i,j)*pearth*fb
        aij(i,j,ij_gvssnd)=aij(i,j,ij_gvssnd)+snowbv(2,i,j)*pearth*fv

      end if
c****
      end do
      end do

      end subroutine ground_e


      subroutine set_roughness_length
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use Dictionary_mod, only : sync_param, get_param
      use PBLCOM, only : roughl
      use fluxes, only : focean, flice
      use ghy_com, only : top_dev_ij
      use ent_com, only : entcells
      use ent_mod
      use ent_drv, only : map_ent2giss !YKIM-temp hack
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      implicit none
      real*8,dimension(N_COVERTYPES) :: fr_cover0, h_cover0
      real*8 :: fr_cover(12), z0_veg
!     original Model II (1983) values (except crops)
!      real*8, parameter :: z0_cover(12) =
!     &     (/0.005d0, 0.01d0, 0.01d0,  0.018d0, 0.32d0, 1.d0, 1.d0,
!     &     2.d0, 0.15d0, 0.005d0, 0.d0, 0.d0 /)
!     these are values which I think are more appropriate (based on
!     info from Monin & Yaglom (~h/7.5) and other literature) -I.A.
      real*8, parameter :: z0_cover(12) =
     &     (/0.005d0, 0.015d0, 0.15d0,  0.5d0, 0.7d0, 1.5d0, 1.5d0,
     &     2.d0, 0.2d0, 0.005d0, 0.d0, 0.d0 /)

      character*80 :: titrrr, titvvv
      real*8 rrr(im,grid%J_STRT_HALO:grid%J_STOP_HALO)
      real*8 vvv(im,grid%J_STRT_HALO:grid%J_STOP_HALO,10)
      real*8 :: roughl_lice = -1.d30
      integer I_0, I_1, J_0, J_1, i, j, k, fid
      titrrr = "roughness length over land"
      rrr = 0.

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     *               I_STRT=I_0, I_STOP=I_1 )

      call sync_param( "roughl_lice", roughl_lice )

      roughl(:,:) = 0.d0
      if(file_exists('ROUGHL')) then
        ! read arbitrary distribution from file
        fid = par_open(grid,'ROUGHL','read')
        call read_dist_data(grid,fid,'roughl',rrr)
        call par_close(grid,fid)
      else
        ! compute non-veg roughness length as in Model II
        do j=J_0,J_1
        do i=I_0,I_1
          if ( focean(i,j) >= 1.d0 ) cycle
          rrr(i,j) = .6d0 * 0.041d0 * top_dev_ij(i,j)**0.71d0
        enddo
        enddo
      endif

      ! make sure that roughness length is always > 0
      do j=J_0,J_1
      do i=I_0,I_1
        if ( focean(i,j) >= 1.d0 ) cycle
        rrr(i,j) = max( rrr(i,j), 0.005d0 ) ! 0.005 is flat desert
      enddo
      enddo

      vvv = 0
      do j=J_0,J_1
        do i=I_0,I_1
            !if ( fearth(i,j) <= 0.d0 ) cycle
          if ( focean(i,j) >= 1.d0 ) cycle
          call ent_get_exports( entcells(i,j),
     &         vegetation_fractions=fr_cover0,
     &         vegetation_heights=h_cover0 )
          call map_ent2giss(fr_cover0,h_cover0,fr_cover) !temp hack: ent pfts->giss veg
          vvv(i,j,:) = fr_cover(1:10)
          z0_veg = (1.d0-flice(i,j))*sum( fr_cover(:)*z0_cover(:) )
     &         + flice(i,j)*0.005d0
          rrr(i,j) = max ( rrr(i,j), z0_veg )
        enddo
      enddo
c**** hack to reset roughl for non-standard land ice fractions
      if ( roughl_lice > 0.d0 ) then
        do j=J_0,J_1
          do i=I_0,I_1
            if ( flice(i,j) > 0.d0 ) rrr(i,j) =
     &           ( rrr(i,j)*(1.d0-flice(i,j)) + roughl_lice*flice(i,j) )
          enddo
        enddo
      endif

      do j=J_0,J_1
        do i=I_0,I_1
          roughl(i,j) = rrr(i,j)
        enddo
      enddo
!!!#endif
      titrrr = " rough len + veg"
      !write(982) titrrr,real(rrr,kind=4)
      !call WRITET_PARALLEL(grid,982,"fort.982",rrr,titrrr)

      do k=1,10
        write(titvvv,*) k
          !write(772) titvvv, real(vvv(:,:,k),kind=4)
        !call WRITET_PARALLEL(grid,773,"fort.773",vvv(:,:,k),titvvv)
      enddo

      end subroutine set_roughness_length


      subroutine accumulate_excess_C(flag)
!@sum accumulate the increment of total carbon stored by LSM
!@auth I. ALeinov
      use constant, only : rhow
      use fluxes, only : focean
      use resolution, only : im,jm
      use geom, only : imaxj,AXYP
      use ghy_com, only : ngm
      use ent_com, only : entcells, excess_C
      use ent_mod, only : ent_get_exports
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      implicit none
      integer, intent(in) :: flag
      real*8, pointer, dimension(:,:), save :: cell_C_old => null()
      real*8 :: cell_C

      integer i,j
      integer :: J_0, J_1 ,I_0,I_1
      integer :: counter = 0

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &     I_STRT=I_0, I_STOP=I_1)

      if (.not. associated(cell_C_old))
     &     allocate( cell_C_old(I_0:I_1, J_0:J_1) )

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          if( focean(i,j) >= 1.d0 ) cycle
          ! call get_fb_fv( fb, fv, i, j )

          call ent_get_exports( entcells(i,j),
     &         C_entcell=cell_C )

          if ( flag == 0 ) then
            cell_C_old(i,j) = cell_C
          else
            excess_C(i,j) = excess_C(i,j)
     &           + cell_C_old(i,j) - cell_C
          endif
        end do
      end do

cddd      write(997,*) 'counter= ', counter
cddd      do j=J_0,J_1
cddd        do i=I_0,I_1
cddd          write(997,*) i,j,excess_C(i,j)
cddd        enddo
cddd      enddo
cddd      counter = counter + 1

      end subroutine accumulate_excess_C


      subroutine soil_spinup_driver
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use geom, only : imaxj
      use ghy_com, only : fearth
      use ent_com, only : entcells
      use ent_mod, only : get_soil_spinup_acc, reset_soil_spinup_acc
     &     , ent_spinup_sbgc_step
      implicit none
      integer, save :: counter=0
ccc   This code should be called right after the land surface time step
ccc   to accumulate the climatology and to eventually apply it to
ccc   spinup the soil carbon. No attempt was made to make this algorithm
ccc   user friendly. You have to manually set up in the definitions below
ccc   the period for which the climatology is accumulated and the number
ccc   of cycles that it will be applied to the soil algorithm.
ccc   Keep in mind that everything should be finished in one run.
ccc   (Climatology is kept in the memory and is not saved.)
ccc   This algorithm assumes assumes that there are 4 land surface
ccc   time steps per hour (hardwired !). If this changes, the algorithm
ccc   needs to be updated.
ccc   Keep in mind that this algorithm resets Ent accumulators which are used
ccc   for this spinup. So it is preferable to call it even if the spinup
ccc   is not performed, to make sure that one does not overflow these
ccc   accumulators.
c**** settings for soil carbon spinup (explicitly assumes 4 setp per hour)
c**** By default everything is done on the fly, so make sure that you 
c**** have enough memory to store the forcings and enough time for: 
c**** Run the model for sbgc_spinup_start-1 steps to get to equilibrium. 
c**** Accumulate forcings for the next sbgc_spinup_hours hours. 
c**** Do the spinup by applying them sbgc_spinup_cycles times. 
c**** Continue the run. 
      logical, parameter :: sbgc_spinup_read_from_file = .false.
      logical, parameter :: sbgc_spinup_dump_to_file = .false.
      
      !integer, parameter :: sbgc_spinup_start = 4*24*365 * 5 +1 ! 1-start now 
      !integer, parameter :: sbgc_spinup_start = 4*24*2 +1 ! wait 2 days 
      integer, parameter :: sbgc_spinup_start = 1 ! 1-start now 
      !integer, parameter :: sbgc_spinup_hours = 24*365*10 ! 10 yrs climatology
      !integer, parameter :: sbgc_spinup_hours = 24*38 ! 1 month climatology
      !integer, parameter :: sbgc_spinup_hours = 24*3 ! 3 days climatology
      integer, parameter :: sbgc_spinup_hours = 24*365 * 9 ! used for C4MIP
      integer, parameter :: sbgc_spinup_steps = sbgc_spinup_hours*4
      !integer, parameter :: sbgc_spinup_cycles = 1000 !how many times to repeat
      integer, parameter :: sbgc_spinup_cycles = 750 ! used for C4MIP
      real*8 :: soil_spinup_acc_buf(7,11)
      real*4, allocatable, save :: soil_spinup_acc_buf_save(:,:,:,:,:)
      logical, save :: first_call = .true.
      integer, save :: sbgc_counter = -1
      integer :: istep, spinup_iter
      
      integer J_0, J_1, I_0, I_1, i, j

      counter = counter + 1
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &     I_STRT=I_0, I_STOP=I_1)

      ! sbgc_spinup - initialize if necessary
      if ( first_call .and. sbgc_spinup ) then
        first_call = .false.
        allocate( soil_spinup_acc_buf_save(
     &       7, 11, I_0:I_1, J_0:J_1, sbgc_spinup_hours ) )
        soil_spinup_acc_buf_save(:,:,:,:,:) = 0.d0
        if ( sbgc_spinup_read_from_file ) then
          do j=J_0,J_1
            read(1000+j) soil_spinup_acc_buf_save(:,:,:,j,:)
          enddo
        endif
      endif

      ! sbgc_spinup - start accumulating the climatology
      if ( sbgc_spinup .and. counter == sbgc_spinup_start ) then
        print *,"starting collecting climatology for sbgc_spinup", J_0
        sbgc_counter = 1
      endif

      if ( sbgc_counter >= 0 ) then
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            if ( fearth(i,j) <= 0.d0 ) cycle
            !print *,"acc ", J_0, sbgc_counter
            call get_soil_spinup_acc(entcells(i,j), soil_spinup_acc_buf)
            soil_spinup_acc_buf_save(:,:,i,j,(sbgc_counter-1)/4+1) =
     &           soil_spinup_acc_buf_save(:,:,i,j,(sbgc_counter-1)/4+1)+
     &           soil_spinup_acc_buf(:,:)
            !call reset_soil_spinup_acc(entcells(i,j))
          enddo
        enddo
      endif

      ! reset accumulator in either case
      do j=J_0,J_1
        do i=I_0,imaxj(j)
          if ( fearth(i,j) <= 0.d0 ) cycle
          call reset_soil_spinup_acc(entcells(i,j))
        enddo
      enddo

      ! sbgc_spinup - do the spinup
      if ( sbgc_counter == sbgc_spinup_steps ) then
        sbgc_counter = -1
        do spinup_iter=1,sbgc_spinup_cycles
          print *,"spinup_iter ", spinup_iter, J_0
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              if ( fearth(i,j) <= 0.d0 ) cycle
              do istep=1,sbgc_spinup_hours
                soil_spinup_acc_buf(:,:) =
     &               soil_spinup_acc_buf_save(:,:,i,j,istep)
                call ent_spinup_sbgc_step(entcells(i,j), 3600.d0,
     &               soil_spinup_acc_buf)
              enddo !istep
            enddo !i
          enddo !j
        enddo
        ! dump forcings to file if requested
        if ( sbgc_spinup_dump_to_file ) then
          do j=J_0,J_1
            write(1000+j) soil_spinup_acc_buf_save(:,:,:,j,:)
          enddo
        endif
      endif

      if ( sbgc_counter >= 0 ) sbgc_counter = sbgc_counter + 1

      end subroutine soil_spinup_driver


      end module soil_drv


      subroutine conserv_wtg(waterg)
!@sum  conserv_wtg calculates ground water incl snow
!@auth Gavin Schmidt
      use constant, only : rhow
      use fluxes, only : focean, flice
      use resolution, only : im,jm
      use geom, only : imaxj,AXYP
      use ghy_com, only : ngm,w_ij,wsn_ij,fr_snow_ij,nsn_ij,fearth
      !use veg_com, only : afb
      use LAKES_COM, only : flake
      use LANDICE_COM,only : MDWNIMP
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      implicit none
!@var waterg ground water (kg/m^2)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO),intent(out)::
     &     waterg

      integer i,j,n
      real*8 wij,fb,fv

C**** define local grid
      integer :: J_0, J_1 ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
      do i=I_0,imaxj(j)
          !if (fearth(i,j).gt.0) then
          if( focean(i,j) + flice(i,j) < 1.d0 ) then
            !fb=afb(i,j)
            !fv=(1.d0-fb)
            call get_fb_fv( fb, fv, i, j )
            wij=fb*sum( w_ij(1:ngm,1,i,j) )
     &       +  fv*sum( w_ij(0:ngm,2,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( wsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( wsn_ij(1:nsn_ij(2,i,j),2,i,j))
            waterg(i,j)=fearth(i,j)*wij*rhow
     &           + flake(i,j)*sum( w_ij(0:ngm,3,i,j) )*rhow
!!! hack to check remove_extra_snow
!!!            waterg(j)=waterg(j)+MDWNIMP(i,j)/axyp(i,j)
          else
            waterg(i,j)=0
          end if
       end do
      end do
      if (HAVE_SOUTH_POLE) waterg(2:im,1) =waterg(1,1)
      if (HAVE_NORTH_POLE) waterg(2:im,jm)=waterg(1,jm)
c****
      end subroutine conserv_wtg


      subroutine conserv_wtg_1(waterg,fearth,flake)
!@sum  conserv_wtg calculates ground water incl snow
!@auth Gavin Schmidt
      use constant, only : rhow
      use resolution, only : im,jm
      use fluxes, only : focean
      use geom, only : imaxj
      use ghy_com, only : ngm,w_ij,wsn_ij,fr_snow_ij,nsn_ij
      !use veg_com, only : afb
      !use LAKES_COM, only : flake
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      implicit none
      real*8,dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                 GRID%J_STRT_HALO:GRID%J_STOP_HALO),
     &     intent(in) :: fearth, flake
!@var waterg ground water (kg/m^2)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO),intent(out)::
     &     waterg

      integer i,j,n
      real*8 wij,fb,fv

C**** define local grid
      integer :: J_0, J_1 ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          !if (fearth(i,j).gt.0) then
          if( focean(i,j) < 1.d0 ) then
            !fb=afb(i,j)
            !fv=(1.d0-fb)
            call get_fb_fv( fb, fv, i, j )
            wij=fb*sum( w_ij(1:ngm,1,i,j) )
     &       +  fv*sum( w_ij(0:ngm,2,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( wsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( wsn_ij(1:nsn_ij(2,i,j),2,i,j))
            waterg(i,j)=fearth(i,j)*wij*rhow
     &           + flake(i,j)*sum( w_ij(0:ngm,3,i,j) )*rhow
          else
            waterg(i,j)=0
          end if
       end do
      end do

      if (HAVE_SOUTH_POLE) waterg(2:im,1) =waterg(1,1)
      if (HAVE_NORTH_POLE) waterg(2:im,jm)=waterg(1,jm)
c****
      end subroutine conserv_wtg_1



      subroutine conserv_htg(heatg)
!@sum  conserv_htg calculates ground energy incl. snow energy
!@auth Gavin Schmidt
      use fluxes, only : focean, flice
      use resolution, only : im,jm
      use geom, only : imaxj, axyp
      use ghy_com, only : ngm,ht_ij,fr_snow_ij,nsn_ij,hsn_ij
     *     ,fearth
      !use veg_com, only : afb
      use LAKES_COM, only : flake
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      implicit none
!@var heatg ground heat (J/m^2)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: heatg

      integer i,j
      real*8 hij,fb,fv

C**** define local grid
      integer J_0, J_1 ,I_0,I_1
      logical :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          heatg(i,j)=0
          !if (fearth(i,j).le.0) cycle
          if ( focean(i,j) + flice(i,j) >= 1.d0 ) cycle
          !fb=afb(i,j)
          !fv=(1.d0-fb)
          call get_fb_fv( fb, fv, i, j )
          hij=fb*sum( ht_ij(1:ngm,1,i,j) )
     &       +  fv*sum( ht_ij(0:ngm,2,i,j) )
     &       +  fb*fr_snow_ij(1,i,j)*sum( hsn_ij(1:nsn_ij(1,i,j),1,i,j))
     &       +  fv*fr_snow_ij(2,i,j)*sum( hsn_ij(1:nsn_ij(2,i,j),2,i,j))
          heatg(i,j)=fearth(i,j)*hij
     &           + flake(i,j)*sum( ht_ij(0:ngm,3,i,j) )
        end do
      end do
      if (HAVE_SOUTH_POLE) heatg(2:im,1) =heatg(1,1)
      if (HAVE_NORTH_POLE) heatg(2:im,jm)=heatg(1,jm)
c****
      end subroutine conserv_htg



      subroutine check_ghy_conservation( flag )
ccc debugging program: cam be put at the beginning and at the end
ccc of the 'surface' to check water conservation
      use constant, only : rhow
      use geom, only : imaxj
      use resolution, only : im,jm
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      use fluxes, only : atmlnd,prec
      use ghy_com, only : ngm,w_ij,ht_ij,snowbv,dz_ij
     *     ,fearth
      !use veg_com, only : afb
      implicit none
      integer flag
      real*8 total_water(im,jm), error_water
      real*8, save :: old_total_water(im,jm)
!      real*8 total_energy(im,jm), error_energy
!      real*8, save :: old_total_energy(im,jm)
      integer i,j,n
      real*8 fb,fv
ccc enrgy check not implemented yet ...

C**** define local grid
      integer J_0, J_1 ,I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do j=J_0,J_1
        do i=I_0,imaxj(j)
          if ( fearth(i,j) <= 0.d0 ) cycle

ccc just checking ...
          do n = 1,ngm
            if ( dz_ij(i,j,n) .le. 0.d0 )
     &           call stop_model('incompatible dz',255)
          enddo

          !fb = afb(i,j)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, i, j )
          total_water(i,j) = fb*sum( w_ij(1:ngm,1,i,j) )
     &         + fv*sum( w_ij(0:ngm,2,i,j) )
     &         + fb*snowbv(1,i,j) + fv*snowbv(2,i,j)
        end do
      end do

      ! call stop_model('just testing...',255)

      if ( flag == 0 ) then
        old_total_water(:,:) = total_water(:,:)
        return
      endif

      do j=J_0,J_1
        do i=I_0,imaxj(j)

          !print *,'fearth = ', i, j, fearth(i,j)

          if ( fearth(i,j) <= 0.d0 ) cycle
          !fb = afb(i,j)
          !fv = 1.d0 - fb
          call get_fb_fv( fb, fv, i, j )
          error_water = ( total_water(i,j) - old_total_water(i,j) )*rhow
     &         - prec(i,j) + atmlnd%evapor(i,j) + atmlnd%runo(i,j)

          !print *, 'err H2O: ', i, j, error_water

  !        if ( abs( error_water ) > 1.d-9 ) print *, 'error'
          if ( abs( error_water ) > 1.d-9 ) call stop_model(  ! was -15
     &         'check_ghy_conservation: water conservation problem',255)

        end do
      end do

      end subroutine check_ghy_conservation


      subroutine compute_water_deficit
      use constant, only : twopi,rhow
      use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
      use ghy_com, only : ngm,imt,LS_NFRAC,dz_ij,q_ij
     &     ,w_ij,fearth
      !use veg_com, only : ala !,afb
      use fluxes, only : focean
      use sle001, only : thm
      use fluxes, only : DMWLDF
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds

      implicit none
cddd      integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1
      integer k,ibv,m
      real*8 :: w_tot(2),w_stor(2)
      real*8 :: w(0:ngm,2),dz(ngm),q(imt,ngm)
cddd      real*8 :: cosday,sinday,alai
      real*8 :: fb,fv

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &     J_STRT=J_0, J_STOP=J_1)

cddd      cosday=cos(twopi/EARTH_DAYS_PER_YEAR*jday)
cddd      sinday=sin(twopi/EARTH_DAYS_PER_YEAR*jday)

      DMWLDF(:,:) = 0.d0

      do j=J_0,J_1
        do i=I_0,I_1

          !if( focean(i,j) >= 1.d0 ) then
          ! this condition should be switched to focean(i,j) >= 1.d0
          ! once all ground arrays are properly initialized for focean(i,j)<1
          if( fearth(i,j) <= 0.d0 ) then
            DMWLDF(i,j) = 0.d0
            cycle
          endif

          w(0:ngm,1:2) = w_ij(0:ngm,1:2,i,j)
          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

          !fb = afb(i,j)
          !fv=1.-fb
          call get_fb_fv( fb, fv, i, j )

          w_stor(:) = 0.d0
          w_tot(:) = 0.d0
          do ibv=1,2
            do k=1,ngm
              do m=1,imt-1
                w_stor(ibv) = w_stor(ibv) + q(m,k)*thm(0,m)*dz(k)
              end do
              w_tot(ibv) = w_tot(ibv) + w(k,ibv)
            end do
          end do

          ! include canopy water here
cddd          alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
cddd          alai=max(alai,1.d0)
cddd          w_stor(2) = w_stor(2) + .0001d0*alai
          w_tot(2) = w_tot(2) + w(0,2)

          ! total water deficit on kg/m^2
          DMWLDF(i,j) = rhow * ( fb*(w_stor(1) - w_tot(1))
     &         + fv*(w_stor(2) - w_tot(2)) )
          DMWLDF(i,j) = max( DMWLDF(i,j), 0.d0 )
        enddo
      enddo


      !print *,"DMWLDF=",DMWLDF

      end subroutine compute_water_deficit


      subroutine init_underwater_soil

!!!! UNFINISHED
      use constant, only : twopi,rhow,shw_kg=>shw
      use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
      use ghy_com, only : ngm,imt,LS_NFRAC,dz_ij,q_ij
     &     ,w_ij,ht_ij,fearth,shc_soil_texture
#ifdef TRACERS_WATER
     &     ,tr_w_ij
      use OldTracer_mod, only: needtrs, itime_tr0
      use TRACER_COM, only : NTM
      use model_com, only : itime
#endif
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
      use sle001, only : thm
      use fluxes, only : focean,DMWLDF
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds

      implicit none
cddd      integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1
      integer k,ibv,m
      real*8 :: w_stor(0:ngm), ht_cap(0:ngm)
      real*8 :: w(0:ngm,2),dz(ngm),q(imt,ngm)
      !!real*8 :: cosday,sinday,alai
      real*8 :: fb,fv
      real*8 shc_layer, aa, tp, tpb, tpv, ficeb, ficev, fice
#ifdef TRACERS_WATER
      integer n
      real*8 wsoil_tot
#endif

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &     J_STRT=J_0, J_STOP=J_1)


cddd      cosday=cos(twopi/EARTH_DAYS_PER_YEAR*jday)
cddd      sinday=sin(twopi/EARTH_DAYS_PER_YEAR*jday)

      do j=J_0,J_1
        do i=I_0,I_1

          if( focean(i,j) >= 1.d0 ) then
            w_ij (0:ngm,3,i,j) = 0.d0
            ht_ij(0:ngm,3,i,j) = 0.d0
            cycle
          endif

          !w(0:ngm,1:2) = w_ij(0:ngm,1:2,i,j)
          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

          call ent_get_exports( entcells(i,j),
     &         canopy_heat_capacity=ht_cap(0) )
          !fb = afb(i,j)
          !fv=1.-fb
          call get_fb_fv( fb, fv, i, j )

          ! compute max water storage and heat capacity
          do k=1,ngm
            w_stor(k) = 0.d0
            do m=1,imt-1
              w_stor(k) = w_stor(k) + q(m,k)*thm(0,m)*dz(k)
            enddo
            shc_layer = 0.d0
            do m=1,imt
              shc_layer = shc_layer + q(m,k)*shc_soil_texture(m)
            enddo
            ht_cap(k) = (dz(k)-w_stor(k)) * shc_layer
          enddo

          ! include canopy water here
          !!alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
          !!alai=max(alai,1.d0)
          !!w_stor(0) = .0001d0*alai*fv
          w_stor(0) = 0.d0
!! set above
          ! we will use as a reference average temperature of the
          ! lowest layer
          call heat_to_temperature( tpb, ficeb,
     &         ht_ij(ngm,1,i,j), w_ij(ngm,1,i,j), ht_cap(ngm) )
          call heat_to_temperature( tpv, ficev,
     &         ht_ij(ngm,2,i,j), w_ij(ngm,2,i,j), ht_cap(ngm) )
          tp   = fb*tpb + fv*tpv
          fice = fb*ficeb + fv*ficev

          ! set underground fraction to tp, fice and saturated water
          !! nothing in canopy layer (ie. temp = 0C)
          w_ij(0,3,i,j) = 0.d0
          ht_ij(0,3,i,j) = 0.d0
          do k=1,ngm
            w_ij(k,3,i,j) = w_stor(k)
            call temperature_to_heat( ht_ij(k,3,i,j),
     &           tp, fice, w_ij(k,3,i,j), ht_cap(k) )
          enddo

#ifdef TRACERS_WATER
      ! set underwater tracers to average land tracers (ignore canopy)
          do k=1,ngm
            wsoil_tot = fb*w_ij(k,1,i,j) + fv*w_ij(k,2,i,j)
            do n=1,ntm
              if (itime_tr0(n).gt.itime) cycle
              if ( .not. needtrs(n) ) cycle
              tr_w_ij(n,k,3,i,j) = 0.d0
              if ( wsoil_tot > 1.d-30 ) then
                tr_w_ij(n,k,3,i,j) = (
     &               fb*tr_w_ij(n,k,1,i,j) + fv*tr_w_ij(n,k,2,i,j)
     &               ) / (wsoil_tot) * w_ij(k,3,i,j) !!! removed /rhow
              end if
            enddo
          enddo
#endif


        enddo
      enddo

      end subroutine init_underwater_soil

      subroutine get_canopy_temperaure(CanTemp, i, j)
!@sum returns canopy temperature for the cell i,j
!@+   returns -1d30 for a cell with no vegetation
      use constant, only : rhow, shw_kg=>shw
      use ghy_com, only : w_ij, ht_ij, fearth
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
      real*8, intent(out) :: CanTemp
      integer, intent(in) :: i, j
      real*8 :: fb, fv, fice, can_ht_cap

      if ( fearth(i,j) <= 0.d0 ) then
        CanTemp = -1.d30
        return
      endif

      call get_fb_fv( fb, fv, i, j )
      if ( fv <= 0.d0 ) then
        CanTemp = -1.d30
        return
      endif

      call ent_get_exports( entcells(i,j),
     &     canopy_heat_capacity=can_ht_cap )
      call heat_to_temperature( CanTemp, fice,
     &     ht_ij(0,2,i,j), w_ij(0,2,i,j), can_ht_cap)

      end subroutine get_canopy_temperaure


      subroutine heat_to_temperature(tp, fice, ht, w, ht_cap)
      use constant, only : rhow, lhm, shw_kg=>shw, shi_kg=>shi
      real*8, intent(out) :: tp, fice
      real*8, intent(in) :: ht, w, ht_cap
      ! volumetric quantities
      real*8, parameter :: lhmv=lhm*rhow, shwv=shw_kg*rhow,
     &     shiv=shi_kg*rhow

      fice = 0.d0
      if( lhmv*w+ht < 0.d0 ) then ! all frozen
        tp = ( ht + w*lhmv )/( ht_cap + w*shiv )
        fice = 1.d0
      else if( ht > 0.d0 ) then ! all melted
        tp = ht /( ht_cap + w*shwv )
      else  ! part frozen
        tp = 0.d0
        if( w > 1d-12 ) fice = -ht /( lhmv*w )
      endif

      end subroutine heat_to_temperature


      subroutine temperature_to_heat(ht, tp, fice, w, ht_cap)
      use constant, only : rhow, lhm, shw_kg=>shw, shi_kg=>shi
      real*8, intent(in) :: tp, fice, w, ht_cap
      real*8, intent(out) :: ht
      ! volumetric quantities
      real*8, parameter :: lhmv=lhm*rhow, shwv=shw_kg*rhow,
     &     shiv=shi_kg*rhow
      real*8 lht_ice

      lht_ice = fice*w*lhmv

      if( tp > 0.d0 ) then
        ht = tp*( ht_cap + w*shwv ) - lht_ice
      else
        ht = tp*( ht_cap + w*shiv ) - lht_ice
      endif

      end subroutine temperature_to_heat


      subroutine update_land_fractions !(jday)

!!!! UNFINISHED
      use constant, only : twopi,rhow,shw_kg=>shw
      use TimeConstants_mod, only: EARTH_DAYS_PER_YEAR
      use ghy_com, only : ngm,imt,dz_ij,q_ij
     &     ,w_ij,ht_ij,fr_snow_ij,fearth
     &     ,snowbv,top_dev_ij
#ifdef TRACERS_WATER
     &     ,tr_w_ij,tr_wsn_ij
      use TRACER_COM, only : NTM
#endif
      !use veg_com, only : ala !,afb
      use LAKES_COM, only : flake, svflake
      use sle001, only : thm
      use fluxes, only : atmlnd,focean,DMWLDF, DGML
#ifdef TRACERS_WATER
     &     ,DTRL
#endif
      use GEOM, only : BYAXYP
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      use soil_drv, only : snow_cover ! conserv_wtg_1
      use snow_drvm, only : snow_cover_same_as_rad

      implicit none
      !! integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1,ibv
      integer k,m
      real*8 :: w_stor(0:ngm)
      real*8 :: dz(ngm),q(imt,ngm)
      !!real*8 :: cosday,sinday,alai
      real*8 :: fb,fv
      real*8 :: dw, dw_soil, dw_lake, dht, dht_soil, dht_lake
      real*8 :: sum_water, ht_per_m3

      real*8 dfrac
#ifdef TRACERS_WATER
      real*8 dtr_soil(ntm),dtr_lake(ntm),dtr(ntm),tr_per_m3(ntm)
#endif
      real*8 tmp_before(0:ngm), tmp_after(0:ngm)
c      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
c     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
c     &     w_before_j, w_after_j
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     fearth_old

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &     J_STRT=J_0, J_STOP=J_1)


      fearth_old = fearth + (flake - svflake)

      ! call conserv_wtg_1(w_before_j,fearth_old,svflake)

      do j=J_0,J_1
        do i=I_0,I_1

          if( focean(i,j) >= 1.d0 ) cycle

          if ( svflake(i,j) == flake(i,j) ) cycle

          call get_fb_fv( fb, fv, i, j )

          if ( flake(i,j) < svflake(i,j) ) then ! lake shrunk
            ! no external fluxes, just part of underwater fraction
            ! is transformed into a land fraction

            ! no changes to underwater fraction
            ! just redistribute added quantities over land fractions


            dfrac = svflake(i,j) - flake(i,j)
            ! make sure old fearth >= 0 (i.e. remove round-off errors)
            dfrac = min(dfrac, fearth(i,j))

            tmp_before = ( ht_ij(0:ngm,3,i,j) )*svflake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)-dfrac) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j)-dfrac)

            do k=1,ngm
              dw  = dfrac*w_ij(k,3,i,j)
              dht = dfrac*ht_ij(k,3,i,j)
              w_ij(k,1:2,i,j) =
     &             (w_ij(k,1:2,i,j)*(fearth(i,j)-dfrac) + dw)
     &             / fearth(i,j)
              ht_ij(k,1:2,i,j) =
     &             (ht_ij(k,1:2,i,j)*(fearth(i,j)-dfrac) + dht)
     &             / fearth(i,j)
#ifdef TRACERS_WATER
              dtr(:) = dfrac*tr_w_ij(:,k,3,i,j)
              do ibv=1,2
                tr_w_ij(:,k,ibv,i,j) =
     &               (tr_w_ij(:,k,ibv,i,j)*(fearth(i,j)-dfrac) + dtr(:))
     &               / fearth(i,j)
              enddo
#endif
            enddo
            !vegetation:
cddd            if ( fv > 0.d0 ) then
cddd              w_ij(0,2,i,j) =
cddd     &             (w_ij(0,2,i,j)*(fearth(i,j)-dfrac) + dw/fv)
cddd     &             / fearth(i,j)
cddd              ht_ij(0,2,i,j) =
cddd     &             (ht_ij(0,2,i,j)*(fearth(i,j)-dfrac) + dht/fv)
cddd     &             / fearth(i,j)
cddd            endif
              w_ij(0,2,i,j) =
     &             (w_ij(0,2,i,j)*(fearth(i,j)-dfrac))
     &             / fearth(i,j)
              ht_ij(0,2,i,j) =
     &             (ht_ij(0,2,i,j)*(fearth(i,j)-dfrac))
     &             / fearth(i,j)
#ifdef TRACERS_WATER
              tr_w_ij(:,0,2,i,j) =
     &             (tr_w_ij(:,0,2,i,j)*(fearth(i,j)-dfrac))
     &             / fearth(i,j)
#endif

            ! change snow fraction
            fr_snow_ij(1:2,i,j) =
     &           fr_snow_ij(1:2,i,j)*(1.d0-dfrac/fearth(i,j))
#ifdef TRACERS_WATER
            ! tr_wsn is spread over entire cell (i.e. *fr_snow)
            tr_wsn_ij(:,:,1:2,i,j) = tr_wsn_ij(:,:,1:2,i,j) *
     &           (1.d0 - dfrac/fearth(i,j))
#endif

            tmp_after = ( ht_ij(0:ngm,3,i,j) )*flake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j))
            !print *,"after", (tmp_after), (tmp_after-tmp_before)

          else if ( flake(i,j) > svflake(i,j) ) then ! lake expanded

            ! no need to change land values
            ! for underwater fraction:

            dz(1:ngm) = dz_ij(i,j,1:ngm)
            q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

            ! compute max water storage and heat capacity
            do k=1,ngm
              w_stor(k) = 0.d0
              do m=1,imt-1
                w_stor(k) = w_stor(k) + q(m,k)*thm(0,m)*dz(k)
              enddo
            enddo
            ! include canopy water here
!!!cddd            alai=ala(1,i,j)+cosday*ala(2,i,j)+sinday*ala(3,i,j)
!!!cddd            alai=max(alai,1.d0)
!!!cddd            w_stor(0) = .0001d0*alai*fv
            ! no uderlake water in canopy
            w_stor(0) = 0.d0
            ! allow any amount of water in upper soil layer
            w_stor(1) = 1.d30

            dfrac = flake(i,j) - svflake(i,j)

            tmp_before = ( ht_ij(0:ngm,3,i,j) )*svflake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)+dfrac) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j)+dfrac)

            sum_water = DMWLDF(i,j)*dfrac/rhow
            if ( sum_water > 1.d-30 ) then
              ht_per_m3 = DGML(i,j)*BYAXYP(I,J)/sum_water
            else
              ht_per_m3 = 0.d0
            endif
#ifdef TRACERS_WATER
            if ( sum_water > 1.d-30 ) then
              tr_per_m3(:) = DTRL(:,i,j)*BYAXYP(I,J)/sum_water
            else
              tr_per_m3(:) = 0.d0
            endif
#endif
            do k=ngm,1,-1  ! do not loop over canopy
              ! dw, dht - total amounts of water and heat added to
              ! underwater fraction
              dw_soil = dfrac*( fb*w_ij(k,1,i,j) + fv*w_ij(k,2,i,j) )
              dw  = min( dfrac*w_stor(k), sum_water + dw_soil )
              dw_lake = dw - dw_soil
              dht_soil = dfrac*( fb*ht_ij(k,1,i,j) + fv*ht_ij(k,2,i,j) )
              dht_lake = dw_lake*ht_per_m3
              dht = dht_soil + dht_lake
              sum_water = sum_water - dw_lake

              w_ij(k,3,i,j) =
     &             (w_ij(k,3,i,j)*svflake(i,j) + dw)/flake(i,j)
              ht_ij(k,3,i,j) =
     &             (ht_ij(k,3,i,j)*svflake(i,j) + dht)/flake(i,j)
#ifdef TRACERS_WATER
              dtr_soil(:) =
     &             dfrac*(fb*tr_w_ij(:,k,1,i,j) + fv*tr_w_ij(:,k,2,i,j))
              dtr_lake(:) = dw_lake*tr_per_m3(:)
              dtr(:) = dtr_soil(:) + dtr_lake(:)
              tr_w_ij(:,k,3,i,j) =
     &             (tr_w_ij(:,k,3,i,j)*svflake(i,j) + dtr(:))/flake(i,j)
#endif
            enddo

            ! dump canopy water into first layer (underwater canopy is 0C)
            dw = dfrac*( fv*w_ij(0,2,i,j) )
            dht = dfrac*( fv*ht_ij(0,2,i,j) )
              w_ij(1,3,i,j) =
     &             w_ij(1,3,i,j) + dw/flake(i,j)
              ht_ij(1,3,i,j) =
     &             ht_ij(1,3,i,j) + dht/flake(i,j)
#ifdef TRACERS_WATER
            dtr(:) = dfrac*( fv*tr_w_ij(:,0,2,i,j) )
            tr_w_ij(:,1,3,i,j) =
     &             tr_w_ij(:,1,3,i,j) + dtr(:)/flake(i,j)
#endif
            tmp_after = ( ht_ij(0:ngm,3,i,j) )*flake(i,j) +
     &           ( ht_ij(0:ngm,1,i,j) )*fb*(fearth(i,j)) +
     &           ( ht_ij(0:ngm,2,i,j) )*fv*(fearth(i,j))

            ! change snow fraction
            if ( fearth(i,j) <= 0.d0 ) then
              print *, "farctions:",i,j,
     &             focean(i,j), fearth(i,j), flake(i,j), svflake(i,j)
              call stop_model("update_land_fractions: fearth<=0",255)
            endif

            fr_snow_ij(1:2,i,j) =
     &           fr_snow_ij(1:2,i,j)*(1.d0+dfrac/fearth(i,j))
#ifdef TRACERS_WATER
            ! tr_wsn is spread over entire cell (i.e. *fr_snow)
            tr_wsn_ij(:,:,1:2,i,j) = tr_wsn_ij(:,:,1:2,i,j) *
     &           (1.d0 + dfrac/fearth(i,j))
#endif


            ! hack to deal with snow in empty fractions (fb, fv)
            if( fb <= 0.d0 )
     &           fr_snow_ij(1,i,j) = min( .95d0, fr_snow_ij(1,i,j) )
            if( fv <= 0.d0 )
     &           fr_snow_ij(2,i,j) = min( .95d0, fr_snow_ij(2,i,j) )

            if ( fr_snow_ij(1,i,j) > 1.d0 .or.
     &           fr_snow_ij(2,i,j) > 1.d0 ) then
              print *,"FR_SNOW_ERROR" ,
     &             i,j,dfrac,fearth(i,j),fb,fv,fr_snow_ij(1:2,i,j)
              call stop_model(
     &             "update_land_fractions: fr_snow_ij > 1",255)
            endif

          endif

c**** Also reset snow fraction for albedo computation
          if ( snow_cover_same_as_rad == 0 ) then
            ! recompute snow fraction using different formula
            do ibv=1,2
               call snow_cover(atmlnd%fr_snow_rad(ibv,i,j),
     &                snowbv(ibv,i,j), top_dev_ij(i,j) )
               atmlnd%fr_snow_rad(ibv,i,j) = min (
     &            atmlnd%fr_snow_rad(ibv,i,j), fr_snow_ij(ibv,i,j) )
            enddo
          else
            ! snow fraction same as in snow model
            atmlnd%fr_snow_rad(:,i,j) = fr_snow_ij(:, i, j)
          endif

          call set_new_ghy_cells_outputs

          ! reset "FLUXES" arrays
          svflake(i,j) = flake(i,j)
          DMWLDF(i,j) = 0.d0
          DGML(i,j) = 0.d0
        enddo
      enddo


      !call conserv_wtg_1(w_after_j,fearth,flake)
      !print *,"UPDATE_LF CONS_WTRG ", w_after_j-w_before_j


      end subroutine update_land_fractions


      subroutine set_new_ghy_cells_outputs
!@sum set output data for newly created earth cells (when lake shrinks)
      use constant, only : rhow,tf,lhe,lhs
      use ghy_com, only : ngm,imt,dz_ij,q_ij
     &     ,w_ij,ht_ij,fr_snow_ij,fearth,qg_ij
     &     ,shc_soil_texture,tearth,wearth,aiearth
     &     ,tsns_ij
#ifdef TRACERS_WATER
     &     ,tr_w_ij
      use OldTracer_mod, only: needtrs,itime_tr0
      use TRACER_COM, only : NTM
      use model_com, only : itime
#endif
#ifdef SCM
      use SCM_COM, only : SCMopt,SCMin
#endif
      use FLUXES, only : atmlnd
      !use veg_com, only : afb
      use sle001, only : thm
      use LAKES_COM, only : flake, svflake
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds

      implicit none
      !---
      real*8, external :: qsat
      real*8 :: EPS=1.d-12
      integer i,j,I_0,I_1,J_0,J_1
      real*8 :: dz(ngm), q(imt,ngm), w_stor(ngm), ht_cap(ngm)
      real*8 tg1, fice, fb, fv, tpb, tpv, ficeb, ficev, shc_layer
      real*8 dfrac, elhx, ps
      integer k, m
#ifdef TRACERS_WATER
      integer n
      real*8 wsoil_tot
#endif

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &     J_STRT=J_0, J_STOP=J_1)

      do j=J_0,J_1
        do i=I_0,I_1

          dfrac = svflake(i,j) - flake(i,j)
          if ( fearth(i,j) < EPS .or. fearth(i,j)-dfrac > EPS ) cycle

          dz(1:ngm) = dz_ij(i,j,1:ngm)
          q(1:imt,1:ngm) = q_ij(i,j,1:imt,1:ngm)

          !fb = afb(i,j)
          !fv=1.-fb
          call get_fb_fv( fb, fv, i, j )

      ! compute max water storage and heat capacity
          do k=1,ngm
            w_stor(k) = 0.d0
            do m=1,imt-1
              w_stor(k) = w_stor(k) + q(m,k)*thm(0,m)*dz(k)
            enddo
            shc_layer = 0.d0
            do m=1,imt
              shc_layer = shc_layer + q(m,k)*shc_soil_texture(m)
            enddo
            ht_cap(k) = (dz(k)-w_stor(k)) * shc_layer
          enddo

          call heat_to_temperature( tpb, ficeb,
     &         ht_ij(1,1,i,j), w_ij(1,1,i,j), ht_cap(ngm) )
          call heat_to_temperature( tpv, ficev,
     &         ht_ij(1,2,i,j), w_ij(1,2,i,j), ht_cap(ngm) )
          tg1   = fb*tpb + fv*tpv
          fice = fb*ficeb + fv*ficev

          elhx=lhe
          if(tg1.lt.0.)  elhx=lhs
          ps=atmlnd%srfp(i,j)
      ! ground humidity to be used on next time step
          qg_ij(i,j) = qsat(tg1+tf,elhx,ps) ! all saturated
      ! snow fraction same as in snow model
          atmlnd%fr_snow_rad(:,i,j) = 0.d0 ! no snow in new cell

c**** snowe used in RADIATION
          atmlnd%snowe(i,j) = 1000.*0.d0
          atmlnd%snow(i,j) = 0.!snowe(i,j)
          atmlnd%snowfr(i,j) = 0d0
          atmlnd%snowdp(i,j) = 0d0
c**** tearth used only internaly in GHY_DRV
          tearth(i,j) = sqrt(sqrt(fb*(tpb+tf)**4 + fv*(tpv+tf)**4)) - tf
          tsns_ij(i,j) = tg1
c**** wearth+aiearth are used in radiation only
          wearth(i,j)=1000.*( fb*w_ij(1,1,i,j)*(1.-ficeb) +
     &         fv*(w_ij(1,2,i,j)*(1.-ficev)) )
          aiearth(i,j)=1000.*( fb*w_ij(1,1,i,j)*ficeb +
     &         fv*w_ij(1,2,i,j)*ficev )
          atmlnd%gtemp(i,j)=tsns_ij(i,j)
          atmlnd%gtempr(i,j) =tearth(i,j)+tf
#ifdef SCM
          if( SCMopt%Tskin )then
            atmlnd%gtemp(i,j) = SCMin%Tskin - tf
            atmlnd%gtempr(i,j) = SCMin%Tskin
          endif
#endif

#ifdef TRACERS_WATER
      ! I use vegetated ground insted of canopy since canopy has 0 H2O
          wsoil_tot = fb*w_ij(1,1,i,j) + fv*w_ij(1,2,i,j)
          do n=1,ntm
            if (itime_tr0(n).gt.itime) cycle
            if ( .not. needtrs(n) ) cycle
            ! should also restrict to TYPE=nWATER ?
            if ( wsoil_tot > 1.d-30 ) then
              atmlnd%gtracer(n,i,j) = (
     &             fb*tr_w_ij(n,1,1,i,j) + fv*tr_w_ij(n,1,2,i,j)
     &             ) / (rhow*wsoil_tot)
            else
              atmlnd%gtracer(n,i,j) = 0.
            end if
         !!! test
            atmlnd%gtracer(n,i,j) = 1.d0
          enddo
#endif
        enddo
      enddo

      end subroutine set_new_ghy_cells_outputs


      subroutine remove_extra_snow_to_ocean ! bad idea ?
      use constant, only : rhow
      use ghy_com, only : nsn_ij, dzsn_ij, wsn_ij, hsn_ij,
     &     fr_snow_ij,fearth,w_ij,ngm,wsn_max
#ifdef TRACERS_WATER
     &     ,tr_wsn_ij
      use TRACER_COM, only : NTM
#endif
      !use veg_com, only : afb
      use LAKES_COM, only : flake
      use LANDICE_COM, only : MDWNIMP, EDWNIMP
#ifdef TRACERS_WATER
     &     ,TRDWNIMP
#endif
      use GEOM, only : AXYP
      Use DIAG_COM, Only: ITEARTH,AIJ=>AIJ_LOC, J_IMPLM,J_IMPLH,
     *                    IJ_IMPMGR,IJ_IMPHGR, JREG
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds

      implicit none
      !! integer, intent(in) :: jday
      !---
      integer i,j,I_0,I_1,J_0,J_1,J_0H,J_1H
      real*8 fbv(2),wsn(3),hsn(3),dzsn(3),fr_snow,wsn_tot,d_wsn,eta
      real*8 dw,dh   ! ,tot0
      integer ibv,nsn,JR

      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1, 
     &     J_STRT=J_0, J_STOP=J_1,
     &     J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)

      do j=J_0,J_1
        do i=I_0,I_1

          if( fearth(i,j) <= 0.d0 ) cycle

          JR=JREG(I,J)

          !fbv(1) = afb(i,j)
          !fbv(2) =1.d0 - fbv(1)
          call get_fb_fv( fbv(1), fbv(2), i, j )

c          tot0=(fbv(1)*sum( w_ij(1:ngm,1,i,j) )
c     &         +  fbv(2)*sum( w_ij(0:ngm,2,i,j) )
c     &         +  fbv(1)*fr_snow_ij(1,i,j)*sum( wsn_ij(1:nsn_ij(1,i,j),1
c     *         ,i,j))+  fbv(2)*fr_snow_ij(2,i,j)*sum( wsn_ij(1:nsn_ij(2
c     *         ,i,j),2,i,j)))*fearth(i,j)*rhow+flake(i,j)*sum(
c     *         w_ij(0:ngm,3,i,j) )*rhow


          do ibv=1,2
            if ( fbv(ibv) <= 0.d0 ) cycle

            nsn = nsn_ij(ibv, i, j)
            wsn = 0 ; hsn = 0 ; dzsn = 0
            wsn(1:nsn) = wsn_ij(1:nsn, ibv, i, j)
            hsn(1:nsn) = hsn_ij(1:nsn, ibv, i, j)
            dzsn(1:nsn) = dzsn_ij(1:nsn, ibv, i, j)
            fr_snow = fr_snow_ij(ibv, i, j)
            wsn_tot = sum( wsn(1:nsn) )
            if ( wsn_tot > WSN_MAX ) then
              ! check if snow structure ok for thick snow
              if ( nsn < 3)
     &             call stop_model("remove_extra_snow: nsn<3",255)
              d_wsn = wsn_tot - WSN_MAX
              ! fraction of snow to be removed:
              eta = d_wsn/(wsn(2)+wsn(3))
              wsn_ij(2:3, ibv, i, j)  = (1.d0-eta)*wsn(2:3)
              hsn_ij(2:3, ibv, i, j)  = (1.d0-eta)*hsn(2:3)
              dzsn_ij(2:3, ibv, i, j) = (1.d0-eta)*dzsn(2:3)
              ! extra water and energy
              dw = eta*(wsn(2)+wsn(3))*fr_snow*fbv(ibv)*rhow
              MDWNIMP(i,j) = MDWNIMP(i,j) + dw*fearth(i,j)*axyp(i,j)
              dh = eta*(hsn(2)+hsn(3))*fr_snow*fbv(ibv)
              EDWNIMP(i,j) = EDWNIMP(i,j) + dh*fearth(i,j)*axyp(i,j)
#ifdef TRACERS_WATER
              TRDWNIMP(1:ntm,i,j) = TRDWNIMP(1:ntm,i,j) +
     &             eta*(
     &             tr_wsn_ij(1:ntm,2,ibv,i,j)+tr_wsn_ij(1:ntm,3,ibv,i,j)
     &             )*fbv(ibv)*fearth(i,j)*axyp(i,j)
              tr_wsn_ij(1:ntm,2:3,ibv,i,j) =
     &             (1.d0-eta)*tr_wsn_ij(1:ntm,2:3,ibv,i,j)
#endif
              CALL INC_AJ(I,J,ITEARTH,J_IMPLH,dh*fearth(i,j))
              CALL INC_AJ(I,J,ITEARTH,J_IMPLM,dw*fearth(i,j))
              CALL INC_AREG(I,J,JR,J_IMPLH,dh*fearth(i,j))
              CALL INC_AREG(I,J,JR,J_IMPLM,dw*fearth(i,j))
              AIJ(I,J,IJ_IMPMGR) = AIJ(I,J,IJ_IMPMGR) + DW*FEARTH(I,J)
              AIJ(I,J,IJ_IMPHGR) = AIJ(I,J,IJ_IMPHGR) + DH*FEARTH(I,J)

              !print *,"remove_extra_snow", i,j,ibv,wsn_tot,eta,dw,dh

            endif
          enddo
c          print
c     *         *,i,j,tot0,(fbv(1)*sum( w_ij(1:ngm,1,i,j) )+  fbv(2)*sum(
c     *         w_ij(0:ngm,2,i,j) )+  fbv(1)*fr_snow_ij(1,i,j)*sum(
c     *         wsn_ij(1:nsn_ij(1,i,j),1,i,j))+  fbv(2)*fr_snow_ij(2,i,j)
c     *         *sum( wsn_ij(1:nsn_ij(2,i,j),2,i,j)))*fearth(i,j)*rhow
c     *         +flake(i,j)*sum(w_ij(0:ngm,3,i,j) )*rhow

        enddo
      enddo

      end subroutine remove_extra_snow_to_ocean


      subroutine get_fb_fv( fb, fv, i, j )
!@sum this is a hack to hyde dependence on Ent/non-Ent vegetation
!@+   in most of the code. It returns fb,fv - fractions of bare and
!@+   vegetated soil
      use ent_com, only : entcells
      use ent_mod, only : ent_get_exports
      implicit none
      real*8, intent(out) :: fb, fv
      integer, intent(in) :: i, j

      call ent_get_exports( entcells(i,j),
     &     fraction_of_vegetated_soil=fv )
      if ( fv > 1.d0 - 1.d-6 ) fv = 1.d0  ! get rid of round-off errors
      if ( fv < 1.d-6 ) fv = 0.d0         ! get rid of round-off errors
      fb = 1.d0 - fv
      end subroutine get_fb_fv

#ifdef CACHED_SUBDD
      subroutine gijlh_defs(arr,nmax,decl_count)
c 3D outputs (model horizontal grid on soil layers).
      use model_com, only : dtsrc,nday
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)

      decl_count = 0

      arr(next()) = info_type_(
     &  sname = 'GT',
     &  lname = 'Soil Temperature Layers 1-6, Land',
     &  units = 'C'
     &     )
c
! This note copied from DIAG.f version:
! 8/13/10: for RELATIVE wetness, edit giss_LSM/GHY.f
! and activate the corresponding lines where wtr_L is set
      arr(next()) = info_type_(
     &  sname = 'GW',
     &  lname = 'Ground Wetness Layers 1-6, Land',
     &  units = 'm'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'GI',
     &  lname = 'Ground Ice Layers 1-6, Land',
     &  units = 'liq. equiv. m'
     &     )

      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine gijlh_defs
#endif /* CACHED_SUBDD */
