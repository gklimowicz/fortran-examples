#include "rundeck_opts.h"

!#define AG2OG_PRECIP_BUNDLE
!#define OG2AG_TOC2SST_BUNDLE
!#define AG2OG_OCEANS_BUNDLE
!#define OG2AG_OCEANS_BUNDLE

!#define BUNDLE_INTERP

#ifdef BUNDLE_INTERP

#ifdef CUBED_SPHERE
      subroutine bundle_interpolation(lstr, remap, copy_np, do_np_avg)
!@auth Igor Aleinov, Denis Gueyffier, Maxwell Kelley 
      USE cs2ll_utils, only : xgridremap_type,xgridremap_lij
      USE ArrayBundle_mod, only : lookup_str, ab_bundle, ab_unbundle
      implicit none
      type (lookup_str) :: lstr
      type (xgridremap_type) :: remap
      integer :: copy_np, do_np_avg
c
      real*8, pointer :: buf_s(:,:,:), buf_d(:,:,:)
      integer :: i,n,im,jm

      call ab_bundle( lstr, buf_s, buf_d )

      if(copy_np>0) then
        jm = copy_np
        do i=2,size(buf_s,2)
          buf_s(:,i,jm) = buf_s(:,1,jm)
        enddo
      endif

      call xgridremap_lij( remap, buf_s, buf_d)

      if(do_np_avg>0) then
        im = size(buf_d,2)
        jm = do_np_avg
        do n=1,size(buf_d,1)
          buf_d(n,:,jm) = sum(buf_d(n,:,jm))/real(im)
        enddo
      endif

      call ab_unbundle( lstr, buf_s, buf_d )

      end subroutine bundle_interpolation
#else
      subroutine bundle_interpolation(lstr, htype, copy_np, do_np_avg)
!@auth Igor Aleinov, Denis Gueyffier, Maxwell Kelley 
      USE hntrp_mod, only : hntrp_type,band_pack_column,hntr8_band_lij
      USE ArrayBundle_mod, only : lookup_str,get_bounds,ab_copy,
     &     ab_bundle,ab_unbundle
      implicit none
      type (lookup_str) :: lstr
      type (hntrp_type) :: htype
      integer :: copy_np, do_np_avg
c
      real*8, pointer :: buf_s(:,:,:), buf_d(:,:,:), buf_band(:,:,:)
      integer :: i,n,im,jm
      integer :: si_0, si_1, sj_0, sj_1,
     &           di_0, di_1, dj_0, dj_1, lm, jmin,jmax

      call get_bounds(lstr, si_0, si_1, sj_0, sj_1,
     &                      di_0, di_1, dj_0, dj_1)

      if (dj_0 .eq. sj_0 .and. dj_1 .eq. sj_1) then   
         call ab_copy(lstr)
         return
      endif

      call ab_bundle( lstr, buf_s, buf_d )

      lm = size(buf_s,1)
      jmin = htype%bpack%jband_strt
      jmax = htype%bpack%jband_stop
      allocate(buf_band(lm,htype%ima,jmin:jmax))
      call band_pack_column(htype%bpack, buf_s, buf_band)
      if(copy_np>0) then
        jm = copy_np
        do i=2,size(buf_band,2)
          buf_band(:,i,jm) = buf_band(:,1,jm)
        enddo
      endif
! unbundle uses not-zero crit. for wts, so init buf_d to 0
! to avoid fpe probs over land points
      buf_d = 0.
      call hntr8_band_lij(buf_band, htype, buf_d)
      deallocate(buf_band)
      if(do_np_avg>0) then
        im = size(buf_d,2)
        jm = do_np_avg
        do n=1,size(buf_d,1)
          buf_d(n,:,jm) = sum(buf_d(n,:,jm))/real(im)
        enddo
      endif
      call ab_unbundle( lstr, buf_s, buf_d )

      end subroutine bundle_interpolation
#endif

      subroutine lon_avg(arr,IM)
!@sum longitudinal average
      implicit none
      real*8 :: arr(1:IM),asum
      integer :: IM

      asum=sum(arr)
      arr(:)=asum/real(IM)

      end subroutine lon_avg

      subroutine copy_pole(arr,IM)
      implicit none
      real*8 :: arr(1:IM)
      integer :: IM
      
      arr(:)=arr(1)

      end subroutine copy_pole

#endif /* BUNDLE_INTERP */

#ifdef AG2OG_PRECIP_BUNDLE
      SUBROUTINE AG2OG_precip(atm,ice)
!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko, Denis Gueyffier

      USE RESOLUTION, only : aIM=>im, aJM=>jm
      USE OCEAN,      only : oIM=>im, oJM=>jm

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE OCEANR_DIM, only : ogrid

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      USE OFLUXES, only : oTRPREC, oTRUNPSI
#endif
#endif
      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI
      use ocean, only : remap_a2o
      USE ArrayBundle_mod
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atm
      type(iceocn_xchng_vars) :: ice
c
      INTEGER N
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT1(:,:),oWEIGHT1(:,:)
      type (lookup_str) :: lstr
      integer :: i,j,l
      integer :: copy_np,do_np_avg

      if(hasNorthPole(agrid)) then
        copy_np = aJM
      else
        copy_np = 0
      endif
      if(hasNorthPole(ogrid)) then
        do_np_avg = oJM
      else
        do_np_avg = 0
      endif

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- need to allocate separate array for different weights
      allocate(aweight1(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight1(oGRID%I_STRT_HALO:oGRID%I_STOP_HALO
     &           ,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      !-- initializing "lstr"
      call ab_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aWEIGHT(:,:) = 1.-ice%RSI(:,:) !!  open ocean fraction
      call ab_add(lstr, aWEIGHT, oWEIGHT, shape(aWEIGHT), 'ij' )

      call ab_add(lstr, atm%PREC, oPREC, shape(atm%PREC), 'ij',
     &     aWEIGHT, oWEIGHT )

      call ab_add(lstr, atm%EPREC, oEPREC, shape(atm%EPREC), 'ij',
     &     aWEIGHT, oWEIGHT)
#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      call ab_add(lstr, atm%TRPREC, oTRPREC, shape(atm%TRPREC), 'lij', 
     &     aWEIGHT, oWEIGHT )
#endif

      aWEIGHT1(:,:) = ice%RSI(:,:)   
!     oWEIGHT1(:,:) = 1.-oWEIGHT(:,:)            ! using the property REGRID(1-RSI)=1-REGRID(RSI)  
      call ab_add(lstr, aWEIGHT1, oWEIGHT1, shape(aWEIGHT1), 'ij')        ! this line should be removed when previous line is uncommented 

      call ab_add(lstr, ice%RUNPSI, oRUNPSI, shape(ice%RUNPSI), 'ij',
     &     aWEIGHT1, oWEIGHT1)

      call ab_add(lstr, ice%SRUNPSI, oSRUNPSI, shape(ice%SRUNPSI), 
     &     'ij', aWEIGHT1, oWEIGHT1)

      call ab_add(lstr, ice%ERUNPSI, oERUNPSI, shape(ice%ERUNPSI), 
     &     'ij', aWEIGHT1, oWEIGHT1)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
      call ab_add(lstr, ice%TRUNPSI, oTRUNPSI, shape(ice%TRUNPSI),'lij',
     &     aWEIGHT1, oWEIGHT1)
#endif

      call ab_add( lstr, ice%RSI, oRSI, shape(ice%RSI), 'ij')

c*   actual interpolation here
      call bundle_interpolation(lstr,remap_A2O,copy_np,do_np_avg)

      deallocate(aweight,oweight)
      deallocate(aweight1,oweight1)

      RETURN
      END SUBROUTINE AG2OG_precip

#else
      SUBROUTINE AG2OG_precip(atm,ice)

!@sum  INT_AG2OG_precip is for interpolation of precipitation
!!       arrays from atmospheric grid to the ocean grid
!@auth Larissa Nazarenko

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist
#endif

      USE OFLUXES, only : ocnatm, ocnice

      USE INT_AG2OG_MOD, only : INT_AG2OG

      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atm
      type(iceocn_xchng_vars) :: ice

      INTEGER N

      REAL*8, allocatable :: aWEIGHT(:,:)
      integer :: ntmocn

#ifdef TRACERS_OCEAN
      ntmocn=tracerlist%getsize()
#endif
      allocate(aweight(atm%I_0H:atm%I_1H,
     &                 atm%J_0H:atm%J_1H))

      aWEIGHT(:,:) = 1.- ice%RSI(:,:)  !!  open ocean fraction

      CALL INT_AG2OG(atm%PREC,ocnatm%PREC, aWEIGHT)

      CALL INT_AG2OG(atm%EPREC,ocnatm%EPREC, aWEIGHT)

      aWEIGHT(:,:) = ice%RSI(:,:)
      CALL INT_AG2OG(ice%RUNPSI,ocnice%RUNPSI, aWEIGHT)

      CALL INT_AG2OG(ice%SRUNPSI,ocnice%SRUNPSI, aWEIGHT)

      CALL INT_AG2OG(ice%ERUNPSI,ocnice%ERUNPSI, aWEIGHT)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(ice%RSI,ocnice%RSI, aWEIGHT)

#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)

      aWEIGHT(:,:) = 1.- ice%RSI(:,:)
      CALL INT_AG2OG(atm%TRPREC,ocnatm%TRPREC, aWEIGHT, NTMocn)

      aWEIGHT(:,:) = ice%RSI(:,:)
      CALL INT_AG2OG(ice%TRUNPSI,ocnice%TRUNPSI, aWEIGHT, NTMocn)
 
#endif

      deallocate(aweight)

      RETURN
      END SUBROUTINE AG2OG_precip
#endif /* ag2og_precip_bundle */

#ifdef OG2AG_TOC2SST_BUNDLE
      SUBROUTINE OG2AG_TOC2SST(atm)
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko, Denis Gueyffier
      USE CONSTANT, only : tf
      USE RESOLUTION, only : aIM=>IM, aJM=>JM
      USE OCEAN, only : oIM=>IM, oJM=>JM
     *                , oFOCEAN_loc=>FOCEAN
     *                , oIMAXJ=>IMAXJ
#ifndef CUBED_SPHERE
      Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      USE DOMAIN_DECOMP_1D, only : hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : oRSI,ocnatm

      USE MODEL_COM, only: nstep=>itime
      use ocean, only : remap_O2A
      USE ArrayBundle_mod
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atm
c
      INTEGER N
      INTEGER IER, I,J,K,L, NT
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S
      INTEGER :: aI_0H,aI_1H, aJ_0H,aJ_1H
      INTEGER :: I_0,I_1, J_0,J_1
      REAL*8 :: UNP,VNP
      REAL*8, ALLOCATABLE :: oTOT_CHLO_loc(:,:)
      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aWEIGHT2(:,:),oWEIGHT2(:,:)

      real*8, dimension(:,:), allocatable :: agtemp,agtemp2,asss,amlhc
      real*8, dimension(:,:,:), allocatable :: opgas_loc,apgas_loc

      type (lookup_str) :: lstr,lstr_uv
      integer :: copy_np,do_np_avg

      if(hasNorthPole(ogrid)) then
        copy_np = oJM
      else
        copy_np = 0
      endif
      if(hasNorthPole(agrid)) then
        do_np_avg = aJM
      else
        do_np_avg = 0
      endif

      aI_0H = atm%i_0h; aI_1H = atm%i_1h
      aJ_0H = atm%j_0h; aJ_1H = atm%j_1h

      allocate(agtemp(aI_0H:aI_1H,aJ_0H:aJ_1H)); agtemp = 0.
      allocate(agtemp2(aI_0H:aI_1H,aJ_0H:aJ_1H)); agtemp2 = 0.
      allocate(asss(aI_0H:aI_1H,aJ_0H:aJ_1H)); asss = 0.
      allocate(amlhc(aI_0H:aI_1H,aJ_0H:aJ_1H)); amlhc = 0.

      !-- need to allocate separate array for different weights
      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &                ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))
      allocate(aweight2(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &                 ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight2(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))


      !-- initializing "lstr"
      call ab_init( lstr, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)
      call ab_init( lstr_uv, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)

      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP

      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)

      nt = ocnatm%gasex_index%getsize()
      ALLOCATE(opgas_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO,
     &                    nt) ,STAT = IER)
      opgas_loc = 0.
      ALLOCATE(apgas_loc(
     &     aGRID%I_STRT_HALO:aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO:aGRID%J_STOP_HALO,nt) ,STAT = IER)
      apgas_loc = 0.

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      oweight(:,oGRID%J_STRT_HALO) = 0.
      oweight(:,oGRID%J_STOP_HALO) = 0.
      call ab_add(lstr, oWEIGHT, aWEIGHT, shape(oWEIGHT),'ij')

      oWEIGHT2(:,:) = oFOCEAN_loc(:,:) * (1.d0 - oRSI(:,:))
      oweight2(:,oGRID%J_STRT_HALO) = 0.
      oweight2(:,oGRID%J_STOP_HALO) = 0.
      call ab_add( lstr, oWEIGHT2, aWEIGHT2, shape(oWEIGHT2),'ij')

      call ab_add(lstr, ocnatm%gtemp, agtemp, shape(ocnatm%gtemp), 
     &     'ij', oWEIGHT, aWEIGHT)
      call ab_add(lstr, ocnatm%gtemp2, agtemp2, shape(ocnatm%gtemp2), 
     &     'ij', oWEIGHT, aWEIGHT)
      call ab_add(lstr, ocnatm%sss, asss, shape(ocnatm%sss), 
     &     'ij', oWEIGHT, aWEIGHT)
      call ab_add(lstr, ocnatm%mlhc, amlhc, shape(ocnatm%mlhc), 
     &     'ij', oWEIGHT, aWEIGHT)
      call ab_add(lstr, ocnatm%ogeoza, atm%ogeoza, shape(ocnatm%ogeoza), 
     &     'ij', oWEIGHT, aWEIGHT)

#ifdef TRACERS_OCEAN
C**** surface tracer concentration

#ifdef TRACERS_WATER
      call ab_add( lstr, ocnatm%gtracer, atm%gtracer,
     &     shape(ocnatm%gtracer), 'lij', oWEIGHT, aWEIGHT)
#endif

      if (ocnatm%chl_defined) then
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation

        DO J=oJ_0,oJ_1
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTOT_CHLO_loc(I,J) = ocnatm%chl(i,j) !tot_chlo(I,J)
            ELSE
              oTOT_CHLO_loc(I,J)=0.
            END IF
          END DO
        END DO

        call ab_add( lstr, oTOT_CHLO_loc, atm%CHL, 
     &     shape(oTOT_CHLO_loc), 'ij', oWEIGHT2, aWEIGHT2) 
        atm%chl_defined=.true.
      endif

!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      if (ocnatm%updated) then
        DO NT = 1,atm%gasex_index%getsize()
          l=atm%gasex_index%at(nt)
          DO J=oJ_0,oJ_1
            DO I=oI_0,oIMAXJ(J)
              IF (oFOCEAN_loc(I,J).gt.0.) THEN
                !pco2 is in uatm, convert to kg,CO2/kg,air
                opgas_loc(I,J,nt) = ocnatm%gtracer(l, I,J)
     .            * atm%vol2mass(nt)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air
              ELSE
                opgas_loc(I,J,nt)=0.
              END IF
            END DO
          END DO
        END DO
        call ab_add( lstr, opgas_loc, apgas_loc, 
     &       shape(opgas_loc), 'ijk', oWEIGHT2, aWEIGHT2) 
      endif
#endif

      call bundle_interpolation(lstr,remap_O2A,copy_np,do_np_avg)

      where(atm%focean.gt.0.)
        atm%gtemp = agtemp
        atm%gtempr = agtemp + tf
      end where
      where(atm%focean.gt.0.) atm%gtemp2 = agtemp2
      atm%sss = 0.
      where(atm%focean.gt.0.) atm%sss = asss
      where(atm%focean.gt.0.) atm%mlhc = amlhc

#ifdef TRACERS_OCEAN
      if (ocnatm%updated .and. nstep.ne.0) then
        do nt=1,atm%gasex_index%getsize()
         !gtracer is first set in TRACER_DRV, then atrac is interpolated
         !here from pco2 in the ocean and later OCNDYN sets gtracer=atrac
         !Therefore in timesetep 0 (cold start) pco2 has not yet been defined
         !and atrac has to be hard coded in here, so that we do not have
         !urealistic tracer flux at the air-sea interface during step0.
          atm%gtracer(atm%gasex_index%at(nt),:,:) = apgas_loc(:,:,nt)
        end do
      endif
#endif


c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      call ab_add(lstr_uv, ocnatm%UOSURF, atm%UOSURF,
     &     shape(ocnatm%UOSURF),'ij')
      call ab_add(lstr_uv, ocnatm%VOSURF, atm%VOSURF,
     &     shape(ocnatm%VOSURF),'ij')

      call bundle_interpolation(lstr_uv,remap_O2A,0,0)

#ifndef CUBED_SPHERE
      if(hasNorthPole(agrid)) then ! latlon atm needs single polar vector
        UNP = SUM(atm%UOSURF(:,aJM)*aCOSI(:))*2/aIM
        VNP = SUM(atm%UOSURF(:,aJM)*aSINI(:))*2/aIM
        atm%UOSURF(1,aJM) = UNP
        atm%VOSURF(1,aJM) = VNP
      endif
#endif

      DEALLOCATE(oTOT_CHLO_loc, opgas_loc)

      deallocate(oweight,aweight,oweight2,aweight2)

      I_0 = atm%I_0
      I_1 = atm%I_1
      J_0 = atm%J_0
      J_1 = atm%J_1

#ifdef TRACERS_WATER
#ifndef TRACERS_OCEAN
      DO J=J_0,J_1
        DO I=I_0,atm%IMAXJ(J)
          IF (atm%FOCEAN(I,J).gt.0.) THEN
            atm%GTRACER(:,I,J)=atm%trw0(:)
          END IF
        END DO
      END DO
#endif
#endif

C**** do poles
      if (atm%HAVE_NORTH_POLE) then
      IF (atm%FOCEAN(1,J_1).gt.0) THEN
        DO I=2,I_1
          atm%GTEMP(I,J_1)=atm%GTEMP(1,J_1)
          atm%GTEMP2(I,J_1)=atm%GTEMP2(1,J_1)
          atm%GTEMPR(I,J_1) =atm%GTEMPR(1,J_1)
          atm%SSS(I,J_1)=atm%SSS(1,J_1)
          atm%MLHC(I,J_1)=atm%MLHC(1,J_1)
          atm%UOSURF(I,J_1) = atm%UOSURF(1,J_1)
          atm%VOSURF(I,J_1) = atm%VOSURF(1,J_1)
          atm%OGEOZA(I,J_1) = atm%OGEOZA(1,J_1)
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
          atm%GTRACER(:,I,J_1)=atm%GTRACER(:,1,J_1)
#endif
        END DO
      END IF
      end if

      if (atm%HAVE_SOUTH_POLE) then
      IF (atm%FOCEAN(1,1).gt.0) THEN
        DO I=2,I_1
          atm%GTEMP(I,1)=atm%GTEMP(1,1)
          atm%GTEMP2(I,1)=atm%GTEMP2(1,1)
          atm%GTEMPR(I,1) =atm%GTEMPR(1,1)
          atm%SSS(I,1)=atm%SSS(1,1)
          atm%MLHC(I,1)=atm%MLHC(1,1)
          atm%UOSURF(I,1) = atm%UOSURF(1,1)
          atm%VOSURF(I,1) = atm%VOSURF(1,1)
          atm%OGEOZA(I,1) = atm%OGEOZA(1,1)
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
          atm%GTRACER(:,I,1)=atm%GTRACER(:,1,1)
#endif
        END DO
      END IF
      end if

      RETURN
      END SUBROUTINE OG2AG_TOC2SST

#else
      SUBROUTINE OG2AG_TOC2SST(atm)
!@sum  OG2AG_oceans: ocean arrays used in the subr. TOC2SST are gathered
!       on the ocean grid, interpolated to the atm. grid, and scattered
!!      on the atm. grid
!@auth Larissa Nazarenko
      USE CONSTANT, only : tf
      USE OCEAN, only : oIM=>IM, oJM=>JM
     *                , oFOCEAN_loc=>FOCEAN
     *                , oIMAXJ=>IMAXJ
      USE OCEANR_DIM, only : ogrid
      USE OFLUXES, only : ocnatm
      USE OFLUXES, only : oRSI
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: tracerlist
#endif
      USE MODEL_COM, only: nstep=>itime
      USE INT_OG2AG_MOD, only : INT_OG2AG

      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atm

      INTEGER IER, I,J, L, NT
      INTEGER :: aI_0H,aI_1H, aJ_0H,aJ_1H, aIM,aJM
      INTEGER oJ_0,oJ_1, oI_0,oI_1, oJ_0S,oJ_1S, oJ_0H,oJ_1H
      integer :: j_0,j_1,i_0,i_1
      REAL*8, allocatable :: oWEIGHT(:,:)
      REAL*8, ALLOCATABLE :: oTOT_CHLO_loc(:,:),opgas_loc(:,:)
      REAL*8, ALLOCATABLE, DIMENSION(:,:)  :: awksp
      real*8 :: unp,vnp

      oI_0 = oGRID%I_STRT
      oI_1 = oGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oJ_0S = oGRID%j_STRT_SKP
      oJ_1S = oGRID%j_STOP_SKP
      oJ_0H = oGRID%j_STRT_HALO
      oJ_1H = oGRID%j_STOP_HALO

      aI_0H = atm%i_0h; aI_1H = atm%i_1h
      aJ_0H = atm%j_0h; aJ_1H = atm%j_1h

      ALLOCATE
     *  (oTOT_CHLO_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     * , STAT = IER)

      allocate(oWEIGHT(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) )


      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(ocnatm%ogeoza, atm%ogeoza, oWEIGHT, .FALSE.)


c Discontinued method for ocean C-grid -> atm A-grid:
c use a variant of INT_OG2AG aware of C-grid staggering.
c May be reinstated later when cubed-sphere INT_OG2AG
c has this capability.
c      oUO1(:,:) = UO(:,:,1)
c      oVO1(:,:) = VO(:,:,1)
c      oWEIGHT(:,:) = 1.d0
c      CALL INT_OG2AG(oUO1,oVO1,aUO1,aVO1, oWEIGHT
c     *             , IVSPO,IVNPO)
c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      oWEIGHT(:,:) = 1.d0

      CALL INT_OG2AG(ocnatm%UOSURF, atm%UOSURF, oWEIGHT, .FALSE.,
     &     AvgPole=.FALSE.)
      CALL INT_OG2AG(ocnatm%VOSURF, atm%VOSURF, oWEIGHT, .FALSE.,
     &     AvgPole=.FALSE.)

#ifndef CUBED_SPHERE
      if(atm%have_north_pole) then ! latlon atm needs single polar vector
        aIM = atm%I_1
        aJM = atm%J_1
        !UNP = SUM(aUO1(:,aJM)*aCOSI(:))*2/aIM
        !VNP = SUM(aUO1(:,aJM)*aSINI(:))*2/aIM
        !aUO1(1,aJM) = UNP
        !aVO1(1,aJM) = VNP
        UNP = SUM(atm%UOSURF(:,aJM)*atm%COSI(:))*2/aIM
        VNP = SUM(atm%UOSURF(:,aJM)*atm%SINI(:))*2/aIM
        atm%UOSURF(1,aJM) = UNP
        atm%VOSURF(1,aJM) = VNP
      endif
#endif

#ifdef TRACERS_OCEAN

#ifdef TRACERS_WATER
C**** surface tracer concentration
      !oWEIGHT(:,:) = MO(:,:,1)*oFOCEAN_loc(:,:)
      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(ocnatm%gtracer,atm%gtracer,oWEIGHT,
     &                         tracerlist%getsize(),atm%focean)
#endif

      if (ocnatm%chl_defined) then
!total ocean chlorophyll. Units are kg,chlorophyll/m3 of seawater
!tot_chlo is defined over all ocean points. Here only use open water
!chorophyll, because that is what is seen by radiation
        DO J=oJ_0,oJ_1
          oWEIGHT(:,J) = oFOCEAN_loc(:,J) * (1.d0 - oRSI(:,J))
          DO I=oI_0,oIMAXJ(J)
            IF (oFOCEAN_loc(I,J).gt.0.) THEN
              oTOT_CHLO_loc(I,J) = ocnatm%chl(i,j) !tot_chlo(I,J)
            ELSE
              oTOT_CHLO_loc(I,J)=0.
            END IF
          END DO
        END DO
        CALL INT_OG2AG(oTOT_CHLO_loc,atm%CHL, oWEIGHT, .FALSE.)
        atm%chl_defined=.true.
      endif

#ifndef STANDALONE_OCEAN
#ifdef TRACERS_GASEXCH_ocean
!partial CO2 pressure in seawater. Units are uatm.
!defined only over open ocean cells, because this is what is
!involved in gas exchage with the atmosphere.
      if (ocnatm%updated) then
        ALLOCATE
     *    (opgas_loc(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO) ,STAT=IER)
        opgas_loc = 0 ! ensure "unused" values at poles don't pollute results.
        do nt=1, atm%gasex_index%getsize()
          DO J=oJ_0,oJ_1
           !oWEIGHT(:,J) = oFOCEAN_loc(:,J)*(1.d0-oRSI(:,J))
            oWEIGHT(:,J) = oFOCEAN_loc(:,J)      !pco2 is calculated at ocean surface; not only open-water
            DO I=oI_0,oIMAXJ(J)
              IF (oFOCEAN_loc(I,J).gt.0.) THEN
                !pco2 is in uatm, convert to kg,CO2/kg,air
                opgas_loc(I,J)=
     &                  ocnatm%gtracer(atm%gasex_index%at(nt),I,J)
     .                * atm%vol2mass(nt)* 1.d-6 ! ppmv (uatm) -> kg,CO2/kg,air
              ELSE
                opgas_loc(I,J)=0.
              END IF
            END DO
          END DO
          CALL INT_OG2AG(opgas_loc,atm%work1, oWEIGHT, .FALSE.)

         !gtracer is first set in TRACER_DRV, then atrac is interpolated
         !here from pco2 in the ocean and later OCNDYN sets gtracer=atrac
         !Therefore in timesetep 0 (cold start) pco2 has not yet been defined
         !and atrac has to be hard coded in here, so that we do not have
         !urealistic tracer flux at the air-sea interface during step0.
          if (nstep.ne.0)
     &          atm%gtracer(atm%gasex_index%at(nt),:,:)=atm%work1
        end do

        deallocate(opgas_loc)
      endif
#endif
#endif
#endif

      DEALLOCATE(oTOT_CHLO_loc)

      I_0 = atm%I_0
      I_1 = atm%I_1
      J_0 = atm%J_0
      J_1 = atm%J_1

#ifdef TRACERS_WATER
#ifndef TRACERS_OCEAN
      DO J=J_0,J_1
        DO I=I_0,atm%IMAXJ(J)
          IF (atm%FOCEAN(I,J).gt.0.) THEN
            atm%GTRACER(:,I,J)=atm%trw0(:)
          END IF
        END DO
      END DO
#endif
#endif

      allocate(awksp(aI_0H:aI_1H, aJ_0H:aJ_1H))
      awksp = 0.

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      atm%SSS(:,:)=0.
      CALL INT_OG2AG(ocnatm%gtemp ,awksp, oWEIGHT, .FALSE.)
      where(atm%focean.gt.0.)
        atm%gtemp = awksp
        atm%gtempr = awksp + tf
      end where
      CALL INT_OG2AG(ocnatm%gtemp2,awksp, oWEIGHT, .FALSE.)
      where(atm%focean.gt.0.) atm%gtemp2 = awksp
      CALL INT_OG2AG(ocnatm%sss   ,awksp, oWEIGHT, .FALSE.)
      where(atm%focean.gt.0.) atm%sss = awksp
      CALL INT_OG2AG(ocnatm%mlhc  ,awksp, oWEIGHT, .FALSE.)
      where(atm%focean.gt.0.) atm%mlhc = awksp

      deallocate(awksp)

      deallocate(oweight)

C**** do poles
      if (atm%HAVE_NORTH_POLE) then
      IF (atm%FOCEAN(1,J_1).gt.0) THEN
        DO I=2,I_1
          atm%GTEMP(I,J_1)=atm%GTEMP(1,J_1)
          atm%GTEMP2(I,J_1)=atm%GTEMP2(1,J_1)
          atm%GTEMPR(I,J_1) =atm%GTEMPR(1,J_1)
          atm%SSS(I,J_1)=atm%SSS(1,J_1)
          atm%MLHC(I,J_1)=atm%MLHC(1,J_1)
          atm%UOSURF(I,J_1) = atm%UOSURF(1,J_1)
          atm%VOSURF(I,J_1) = atm%VOSURF(1,J_1)
          atm%OGEOZA(I,J_1) = atm%OGEOZA(1,J_1)
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
#ifndef STANDALONE_OCEAN
          atm%GTRACER(:,I,J_1)=atm%GTRACER(:,1,J_1)
#endif
#endif
        END DO
      END IF
      end if

      if (atm%HAVE_SOUTH_POLE) then
      IF (atm%FOCEAN(1,1).gt.0) THEN
        DO I=2,I_1
          atm%GTEMP(I,1)=atm%GTEMP(1,1)
          atm%GTEMP2(I,1)=atm%GTEMP2(1,1)
          atm%GTEMPR(I,1) =atm%GTEMPR(1,1)
          atm%SSS(I,1)=atm%SSS(1,1)
          atm%MLHC(I,1)=atm%MLHC(1,1)
          atm%UOSURF(I,1) = atm%UOSURF(1,1)
          atm%VOSURF(I,1) = atm%VOSURF(1,1)
          atm%OGEOZA(I,1) = atm%OGEOZA(1,1)
#if (defined TRACERS_WATER) || (defined TRACERS_OCEAN)
#ifndef STANDALONE_OCEAN
          atm%GTRACER(:,I,1)=atm%GTRACER(:,1,1)
#endif
#endif
        END DO
      END IF
      end if

      RETURN
      END SUBROUTINE OG2AG_TOC2SST

#endif /* OG2AG_TOC2SST_BUNDLE */

#ifdef AG2OG_OCEANS_BUNDLE

      SUBROUTINE AG2OG_oceans(atm,ice)
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko, Denis Gueyffier

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM, oJM=>JM
     &                , oSINI=>SINIC, oCOSI=>COSIC 
#if (defined TRACERS_OCEAN)
      USE OCN_TRACER_COM, only: tracerlist
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid
      use ofluxes, only : ocnatm
      USE OFLUXES, only : oRSI, oSOLARw, oSOLARi, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMUA,oDMVA
     *     , oGMELT, oEGMELT
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
     *     , oTRGMELT
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
#endif
#endif
      Use GEOM,  only : aIMAXJ=>IMAXJ
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
#ifndef CUBED_SPHERE
       Use GEOM,  only : aCOSI=>COSIP,aSINI=>SINIP
#endif

      use ocean, only : remap_A2O
      USE ArrayBundle_mod

      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      use runtimecontrols_mod, only: ocn_cfc, tracers_oceanbiology
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atm
      type(iceocn_xchng_vars) :: ice
c
      REAL*8, DIMENSION(:,:), ALLOCATABLE ::
     &     aRSIwt,oRSIwt,
     &     aROCwt,oROCwt,
     &     aOCNwt,oOCNwt,
     &     aDMUA1tmp,aDMVA1tmp,
     &     atmp03,atmp04,atmp05
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE ::
     &     atmp09
      REAL*8 :: aDMUAnp,aDMVAnp,aDMUAsp,aDMVAsp
      INTEGER I,J,N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER aJ_0H,aJ_1H, aI_0H,aI_1H
      INTEGER oJ_0,oJ_1
      INTEGER oJ_0H,oJ_1H, oI_0H,oI_1H
      REAL*8,
     *     DIMENSION(aGRID%I_STRT:aGRID%I_STOP,
     *     aGRID%J_STRT:aGRID%J_STOP)::
     *     aFact
      REAL*8 :: oDMUA1sp,oDMVA1sp,oDMUA1np,oDMVA1np
      type (lookup_str) :: lstr,lstr_uv
      integer :: copy_np,do_np_avg

      aJ_0H = aGRID%J_STRT_HALO
      aJ_1H = aGRID%J_STOP_HALO
      aI_0H = aGRID%I_STRT_HALO
      aI_1H = aGRID%I_STOP_HALO
      oJ_0H = oGRID%J_STRT_HALO
      oJ_1H = oGRID%J_STOP_HALO
      oI_0H = oGRID%I_STRT_HALO
      oI_1H = oGRID%I_STOP_HALO

      allocate(aRSIwt(aI_0H:aI_1H,aJ_0H:aJ_1H))
      allocate(oRSIwt(oI_0H:oI_1H,oJ_0H:oJ_1H))
      allocate(aROCwt(aI_0H:aI_1H,aJ_0H:aJ_1H))
      allocate(oROCwt(oI_0H:oI_1H,oJ_0H:oJ_1H))
      allocate(aOCNwt(aI_0H:aI_1H,aJ_0H:aJ_1H))
      allocate(oOCNwt(oI_0H:oI_1H,oJ_0H:oJ_1H))

      !-- initializing "lstr"
      call ab_init( lstr, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )
      call ab_init( lstr_uv, aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO,
     &     oGRID%I_STRT_HALO, oGRID%I_STOP_HALO,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO )

      aJ_0 = aGRID%j_STRT
      aJ_1 = aGRID%j_STOP
      aI_0 = aGRID%I_STRT
      aI_1 = aGRID%I_STOP
      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP

      if(hasNorthPole(agrid)) then
        copy_np = aJM
      else
        copy_np = 0
      endif
      if(hasNorthPole(ogrid)) then
        do_np_avg = oJM
      else
        do_np_avg = 0
      endif

      call ab_add( lstr, ice%RSI, oRSI, shape(ice%RSI),'ij')

      allocate(atmp03(aI_0H:aI_1H,aJ_0H:aJ_1H))
      allocate(atmp04(aI_0H:aI_1H,aJ_0H:aJ_1H))
      allocate(atmp05(aI_0H:aI_1H,aJ_0H:aJ_1H))
      atmp03=0.d0
      atmp04=0.d0
      atmp05=0.d0

      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (atm%focean(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(atm%focean(I,J))
            atmp03(I,J) = ice%MELTI(I,J)*aFact(I,J)
            atmp04(I,J) = ice%EMELTI(I,J)*aFact(I,J)
            atmp05(I,J) = ice%SMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO

      call ab_add( lstr, atm%FLOWO, oFLOWO, shape(atm%FLOWO), 'ij')
      call ab_add( lstr, atm%EFLOWO, oEFLOWO, shape(atm%EFLOWO), 'ij')
      call ab_add( lstr, atmp03, oMELTI, shape(atmp03), 'ij')
      call ab_add( lstr, atmp04, oEMELTI, shape(atmp04), 'ij')
      call ab_add( lstr, atmp05, oSMELTI, shape(atmp05), 'ij')
      call ab_add( lstr, atm%GMELT, oGMELT, shape(atm%GMELT), 'ij')
      call ab_add( lstr, atm%EGMELT, oEGMELT, shape(atm%EGMELT), 'ij')
      call ab_add( lstr, ice%APRESS, oAPRESS, shape(ice%APRESS), 'ij')

      aRSIwt(:,:) = ice%RSI(:,:)
      call ab_add( lstr, aRSIwt, oRSIwt, shape(aRSIwt), 'ij')
      aROCwt(:,:) = 1.d0 - ice%RSI(:,:)
      ! why not just set oROCwt = 1-oRSIwt?
      call ab_add( lstr, aROCwt, oROCwt, shape(aROCwt), 'ij')
      aOCNwt(:,:) = atm%focean(:,:)
      call ab_add(lstr, aOCNwt, oOCNwt, shape(aOCNwt), 'ij')


      call ab_add( lstr, ice%RUNOSI, oRUNOSI, shape(ice%RUNOSI),
     &     'ij', aRSIwt, oRSIwt)
      call ab_add( lstr, ice%ERUNOSI, oERUNOSI, shape(ice%ERUNOSI), 
     &     'ij', aRSIwt, oRSIwt)
      call ab_add( lstr, ice%SRUNOSI, oSRUNOSI, shape(ice%SRUNOSI), 
     &     'ij', aRSIwt, oRSIwt)
      call ab_add( lstr, ice%SOLAR, oSOLARi, shape(ice%SOLAR), 
     &     'ij', aRSIwt, oRSIwt)

      call ab_add(lstr, atm%E0, oE0, shape(atm%E0), 'ij',
     &     aROCwt, oROCwt)
      call ab_add(lstr, atm%EVAPOR, oEVAPOR, shape(atm%EVAPOR),
     &     'ij', aROCwt, oROCwt)
      call ab_add(lstr, atm%SOLAR, oSOLARw, shape(atm%SOLAR),
     &     'ij', aROCwt, oROCwt)


#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      allocate(atmp09(tracerlist%getsize(),
     &      aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &     ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      atmp09=0.d0

      call ab_add( lstr, atm%TRFLOWO,oTRFLOWO,shape(atm%TRFLOWO),'lij')


      DO J=aJ_0,aJ_1
        DO I=aI_0,aIMAXJ(J)
          IF (atm%focean(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/atm%focean(I,J)
            atmp09(:,I,J) = ice%TRMELTI(:,I,J)*aFact(I,J)
          END IF
        END DO
      END DO

      call ab_add( lstr, atmp09,oTRMELTI,shape(atmp09),'lij')

      call ab_add( lstr, atm%TRGMELT,oTRGMELT,shape(atm%TRGMELT),'lij')
      call ab_add( lstr, ice%TRUNOSI,oTRUNOSI,shape(ice%TRUNOSI),'lij',
     &     aRSIwt,oRSIwt)
      call ab_add(lstr,atm%TREVAPOR,oTREVAPOR,shape(atm%TREVAPOR),'lij',
     &     aROCwt,oROCwt)

#ifdef TRACERS_DRYDEP
      call ab_add(lstr,atm%TRDRYDEP,oTRDRYDEP,shape(atm%TRDRYDEP),'lij')
#endif

#endif
#endif
#ifdef TRACERS_OCEAN
      call ab_add(lstr,atm%TRGASEX,ocnatm%TRGASEX,
     &     shape(atm%TRGASEX),'lij')
#endif

      if (allocated(atm%dirvis)) then
       call ab_add(lstr,atm%DIRVIS,ocnatm%DIRVIS,shape(atm%DIRVIS),'ij',
     &     aOCNwt,oOCNwt)
       call ab_add(lstr,atm%DIFVIS,ocnatm%DIFVIS,shape(atm%DIFVIS),'ij',
     &     aOCNwt,oOCNwt)
       call ab_add(lstr,atm%DIRNIR,ocnatm%DIRNIR,shape(atm%DIRNIR),'ij',
     &     aOCNwt,oOCNwt)
       call ab_add(lstr,atm%DIFNIR,ocnatm%DIFNIR,shape(atm%DIFNIR),'ij',
     &     aOCNwt,oOCNwt)
      endif

      if (tracers_oceanbiology.or.ocn_cfc) then
        call ab_add( lstr,atm%COSZ1,ocnatm%COSZ1,shape(atm%COSZ1),'ij',
     &     aOCNwt,oOCNwt)
        call ab_add( lstr,atm%WSAVG,ocnatm%WSAVG,shape(atm%WSAVG),'ij',
     &     aOCNwt,oOCNwt)
      endif

c*   actual interpolation here
      call bundle_interpolation(lstr,remap_A2O,copy_np,do_np_avg)
c*

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      deallocate(atmp09)
#endif
#endif


      allocate(aDMUA1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDMVA1tmp(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))

#ifndef CUBED_SPHERE
      if (oIM .eq. aIM .and. oJM .eq. aJM) then
         aDMUA1tmp(:,:) = atm%DMUA(:,:)
         aDMVA1tmp(:,:) = atm%DMVA(:,:)
      else
         
         if (hasNorthPole(agrid)) then
            aDMUAnp = atm%DMUA(1,aJM)
            aDMVAnp = atm%DMVA(1,aJM)
            aDMUA1tmp(:,aJM) = aDMUAnp*aCOSI(:) + aDMVAnp*aSINI(:)
            aDMVA1tmp(:,aJM) = aDMVAnp*aCOSI(:) - aDMUAnp*aSINI(:)
         endif
         if (hasSouthPole(agrid)) then
            aDMUAsp = atm%DMUA(1,1)
            aDMVAsp = atm%DMVA(1,1)
            aDMUA1tmp(:,1) = aDMUAsp*aCOSI(:) - aDMVAsp*aSINI(:)
            aDMVA1tmp(:,1) = aDMVAsp*aCOSI(:) + aDMUAsp*aSINI(:)
         endif

         do j=max(2,aGRID%J_STRT_HALO),min(aJM-1,aGRID%J_STOP_HALO) ! exclude poles
            do i=1,aIMAXJ(J)
               if (atm%focean(i,j).gt.0.) then
                  aDMUA1tmp(i,j) = atm%DMUA(i,j) 
               endif
            enddo
         enddo
         
         do j=max(2,aGRID%J_STRT_HALO),min(aJM-1,aGRID%J_STOP_HALO) ! exclude poles
            do i=1,aIMAXJ(J)
               if (atm%focean(i,j).gt.0.) then
                  aDMVA1tmp(i,j) = atm%DMVA(i,j) 
               endif
            enddo
         enddo
         
      endif
#else
      do j=aGRID%J_STRT_HALO,aGRID%J_STOP_HALO
         do i=aGRID%I_STRT_HALO,aGRID%I_STOP_HALO
            if (atm%focean(i,j).gt.0.) then
               aDMUA1tmp(i,j) = atm%DMUA(i,j) 
            endif
         enddo
      enddo 
      
      do j=aGRID%J_STRT_HALO,aGRID%J_STOP_HALO
         do i=aGRID%I_STRT_HALO,aGRID%I_STOP_HALO
            if (atm%focean(i,j).gt.0.) then
               aDMVA1tmp(i,j) = atm%DMVA(i,j) 
            endif
         enddo
      enddo
#endif /* not CUBED_SPHERE */

      call ab_add(lstr_uv, aDMUA1tmp, oDMUA, shape(aDMUA1tmp), 'ij',
     &     aROCwt, oROCwt)
      call ab_add(lstr_uv, aDMVA1tmp, oDMVA, shape(aDMVA1tmp), 'ij',
     &     aROCwt, oROCwt)
      call bundle_interpolation(lstr_uv,remap_A2O,0,0)

#ifdef CUBED_SPHERE
      if (hasSouthPole(ogrid)) then
         oDMUA1sp =  SUM(oDMUA(:,  1)*oCOSI(:))*2/oIM
         oDMVA1sp = -SUM(oDMUA(:,  1)*oSINI(:))*2/oIM
         oDMUA(1,  1) = oDMUA1sp
         oDMVA(1,  1) = oDMVA1sp
      endif
      if (hasNorthPole(oGrid)) then
         oDMUA1np =  SUM(oDMUA(:,oJM)*oCOSI(:))*2/oIM
         oDMVA1np =  SUM(oDMUA(:,oJM)*oSINI(:))*2/oIM
         oDMUA(1,oJM) = oDMUA1np
         oDMVA(1,oJM) = oDMVA1np
      endif
#else
      if (oIM .ne. aIM .or. oJM .ne. aJM) then
         
         if (hasSouthPole(ogrid)) then
            oDMUA1sp =  SUM(oDMUA(:,  1)*oCOSI(:))*2/oIM
            oDMVA1sp = -SUM(oDMUA(:,  1)*oSINI(:))*2/oIM
            oDMUA(1,  1) = oDMUA1sp
            oDMVA(1,  1) = oDMVA1sp
         endif
         if (hasNorthPole(oGrid)) then
            oDMUA1np =  SUM(oDMUA(:,oJM)*oCOSI(:))*2/oIM
            oDMVA1np =  SUM(oDMUA(:,oJM)*oSINI(:))*2/oIM
            oDMUA(1,oJM) = oDMUA1np
            oDMVA(1,oJM) = oDMVA1np
         endif

      endif
#endif

      ! ocean model still wants ice-masked wind stress
      oDMUA = oDMUA*oROCwt
      oDMVA = oDMVA*oROCwt

      deallocate(aRSIwt,oRSIwt,aROCwt,oROCwt,aOCNwt,oOCNwt)

      deallocate(atmp03,atmp04,atmp05,
     &     aDMUA1tmp,aDMVA1tmp)


      END SUBROUTINE AG2OG_oceans
#else

      SUBROUTINE AG2OG_oceans(atm,ice)
!@sum  AG2OG_oceans: all atmospheric arrays used in the subr. OCEANS are gathered
!       on the atmospheric grid, interpolated to the ocean grid, and scattered
!!      on the ocean grid
!@auth Larissa Nazarenko
      USE MODEL_COM, only: nstep=>itime
#if (defined TRACERS_OCEAN)
      USE OCN_TRACER_COM, only: tracerlist
#endif
      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : ocnatm, ocnice

      USE INT_AG2OG_MOD, only : INT_AG2OG

      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      use runtimecontrols_mod, only: ocn_cfc, tracers_oceanbiology
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atm
      type(iceocn_xchng_vars) :: ice

      INTEGER I,J, N
      INTEGER aJ_0,aJ_1, aI_0,aI_1
      INTEGER aJ_0H,aJ_1H, aI_0H,aI_1H
      INTEGER oJ_0,oJ_1,oI_0,oI_1

      real*8, dimension(:,:), allocatable :: aweight,atmp,aFact
      REAL*8, allocatable :: atmp2(:,:,:)
      integer :: ntmocn

#ifdef TRACERS_OCEAN
      ntmocn=tracerlist%getsize()
#endif

      aJ_0 = atm%J_0
      aJ_1 = atm%J_1
      aI_0 = atm%I_0
      aI_1 = atm%I_1

      aJ_0H = atm%J_0H
      aJ_1H = atm%J_1H
      aI_0H = atm%I_0H
      aI_1H = atm%I_1H

      allocate ( aweight(aI_0H:aI_1H,aJ_0H:aJ_1H) )
      allocate ( atmp(aI_0H:aI_1H,aJ_0H:aJ_1H) )
      allocate ( aFact(aI_0H:aI_1H,aJ_0H:aJ_1H) )

      oJ_0 = oGRID%j_STRT
      oJ_1 = oGRID%j_STOP
      oI_0 = oGRID%i_STRT
      oI_1 = oGRID%i_STOP

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(ice%RSI,ocnice%RSI, aWEIGHT)
      
      CALL INT_AG2OG(atm% FLOWO,ocnatm% FLOWO, aWEIGHT)
      CALL INT_AG2OG(atm%EFLOWO,ocnatm%EFLOWO, aWEIGHT)

      CALL INT_AG2OG(atm% GMELT,ocnatm% GMELT, aWEIGHT)
      CALL INT_AG2OG(atm%EGMELT,ocnatm%EGMELT, aWEIGHT)

      atmp=0.
      DO J=aJ_0,aJ_1
        DO I=aI_0,atm%IMAXJ(J)
          IF (atm%FOCEAN(I,J).gt.0.) THEN
            aFact(I,J) = 1.d0/(atm%FOCEAN(I,J))
            atmp(I,J) = ice%MELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,ocnice%MELTI, aWEIGHT)

      atmp=0.
      DO J=aJ_0,aJ_1
        DO I=aI_0,atm%IMAXJ(J)
          IF (atm%FOCEAN(I,J).gt.0.) THEN
            atmp(I,J) = ice%EMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,ocnice%EMELTI, aWEIGHT)

      atmp=0.
      DO J=aJ_0,aJ_1
        DO I=aI_0,atm%IMAXJ(J)
          IF (atm%FOCEAN(I,J).gt.0.) THEN
            atmp(I,J) = ice%SMELTI(I,J)*aFact(I,J)
          END IF
        END DO
      END DO
      CALL INT_AG2OG(atmp,ocnice%SMELTI, aWEIGHT)

      aWEIGHT(:,:) = ice%RSI(:,:)
      CALL INT_AG2OG(ice%RUNOSI,ocnice%RUNOSI, aWEIGHT)

      CALL INT_AG2OG(ice%ERUNOSI,ocnice%ERUNOSI, aWEIGHT)

      CALL INT_AG2OG(ice%SRUNOSI,ocnice%SRUNOSI, aWEIGHT)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(ice%APRESS,ocnice%APRESS, aWEIGHT)

      aWEIGHT(:,:) = 1.d0 - ice%RSI(:,:)
      CALL INT_AG2OG(atm%E0,ocnatm%E0, aWEIGHT)

      CALL INT_AG2OG(atm%EVAPOR,ocnatm%EVAPOR, aWEIGHT)

      CALL INT_AG2OG(atm%SOLAR,ocnatm%SOLAR, aWEIGHT)

      aWEIGHT(:,:) = ice%RSI(:,:)
      CALL INT_AG2OG(ice%SOLAR,ocnice%SOLAR, aWEIGHT)

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      aWEIGHT(:,:) = 1.d0

      CALL INT_AG2OG(atm%TRFLOWO,ocnatm%TRFLOWO, aWEIGHT, ntmocn)

      allocate ( atmp2(ntmocn,aI_0H:aI_1H,aJ_0H:aJ_1H) )
      atmp2=0.

      DO N=1,ntmocn
        DO J=aJ_0,aJ_1
          DO I=aI_0,atm%IMAXJ(J)
            IF (atm%FOCEAN(I,J).gt.0.) THEN
              aFact(I,J) = 1.d0/(atm%FOCEAN(I,J))
              atmp2(N,I,J) = ice%TRMELTI(N,I,J)*aFact(I,J)
            END IF
          END DO
        END DO
      END DO
      CALL INT_AG2OG(atmp2,ocnice%TRMELTI, aWEIGHT, ntmocn)
      deallocate(atmp2)

      aWEIGHT(:,:) = ice%RSI(:,:)
      CALL INT_AG2OG(ice%TRUNOSI,ocnice%TRUNOSI, aWEIGHT, ntmocn)

      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(atm%TRGMELT,ocnatm%TRGMELT, aWEIGHT, ntmocn)

      aWEIGHT(:,:) = 1.d0 - ice%RSI(:,:)
      CALL INT_AG2OG(atm%TREVAPOR,ocnatm%TREVAPOR, aWEIGHT, ntmocn)

#ifdef TRACERS_DRYDEP
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(atm%TRDRYDEP,ocnatm%TRDRYDEP, aWEIGHT, ntmocn)
#endif
#endif
#endif

#ifdef TRACERS_OCEAN
      aWEIGHT(:,:) = 1.d0
      CALL INT_AG2OG(atm%TRGASEX,ocnatm%TRGASEX, aWEIGHT,
     &     atm%gasex_index%getsize())
#endif

      if (allocated(atm%dirvis)) then
        aWEIGHT(:,:) = atm%FOCEAN(:,:)
        CALL INT_AG2OG(atm%DIRVIS,ocnatm%DIRVIS, aWEIGHT)
        CALL INT_AG2OG(atm%DIFVIS,ocnatm%DIFVIS, aWEIGHT)
        CALL INT_AG2OG(atm%DIRNIR,ocnatm%DIRNIR, aWEIGHT)
        CALL INT_AG2OG(atm%DIFNIR,ocnatm%DIFNIR, aWEIGHT)
      endif

      if (tracers_oceanbiology.or.ocn_cfc) then
        aWEIGHT(:,:) = atm%FOCEAN(:,:)
        CALL INT_AG2OG(atm%COSZ1,ocnatm%COSZ1, aWEIGHT)
        CALL INT_AG2OG(atm%WSAVG,ocnatm%WSAVG, aWEIGHT)
        ocnatm%WSAVG =ocnatm%WSAVG*(1.d0-ocnice%RSI)
      endif

      aWEIGHT(:,:) = 1.d0 - ice%RSI(:,:)      
      CALL INT_AG2OG(atm%DMUA,atm%DMVA,
     &     ocnatm%DMUA,ocnatm%DMVA,
     &     aWEIGHT,atm%FOCEAN,atm%SINI,atm%COSI)
      ! ocean model still wants ice-masked wind stress
      do j=oj_0,oj_1
        ocnatm%DMUA(:,j) = ocnatm%DMUA(:,j)*(1.d0-ocnice%RSI(:,j))
        ocnatm%DMVA(:,j) = ocnatm%DMVA(:,j)*(1.d0-ocnice%RSI(:,j))
      enddo

      deallocate(aweight,atmp,aFact)

      END SUBROUTINE AG2OG_oceans
#endif /*ag2og_oceans_bundle */


#ifdef OG2AG_OCEANS_BUNDLE      
      SUBROUTINE OG2AG_oceans(ice)
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko, Denis Gueyffier

      USE RESOLUTION, only : aIM=>IM, aJM=>JM

      USE OCEAN, only : oIM=>IM,oJM=>JM
     *                , oFOCEAN_loc=>FOCEAN

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: tracerlist
#endif

      USE DOMAIN_DECOMP_ATM, only : agrid=>grid
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole
      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : oDMSI, oDHSI, oDSSI
#ifdef TRACERS_OCEAN
     *     , oDTRSI
#endif

      use ocean, only : remap_O2A
      USE ArrayBundle_mod
      use domain_decomp_1d, only: hasNorthPole, hasSouthPole

      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      IMPLICIT NONE
      type(iceocn_xchng_vars) :: ice

      REAL*8, allocatable :: aWEIGHT(:,:),oWEIGHT(:,:)
      REAL*8, allocatable :: aDMSItmp(:,:,:)
      REAL*8, allocatable :: aDHSItmp(:,:,:)
      REAL*8, allocatable :: aDSSItmp(:,:,:)
      REAL*8, allocatable :: otmp(:,:,:),atmp(:,:,:)
      integer :: i,j,l
      type (lookup_str) :: lstr
      integer :: copy_np,do_np_avg
#ifdef TRACERS_OCEAN
      integer :: ntmocn
#endif

#ifdef TRACERS_OCEAN
      ntmocn=tracerlist%getsize()
#endif
      if(hasNorthPole(ogrid)) then
        copy_np = oJM
      else
        copy_np = 0
      endif
      if(hasNorthPole(agrid)) then
        do_np_avg = aJM
      else
        do_np_avg = 0
      endif

      allocate(aweight(aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(oweight(oIM,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO))

      allocate(aDMSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDHSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))
      allocate(aDSSItmp(2,aGRID%I_STRT_HALO:aGRID%I_STOP_HALO
     &           ,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO))


      !-- initializing "lstr"
      call ab_init( lstr, 1, oIM,
     &     oGRID%J_STRT_HALO, oGRID%J_STOP_HALO,
     &     aGRID%I_STRT_HALO, aGRID%I_STOP_HALO,
     &     aGRID%J_STRT_HALO, aGRID%J_STOP_HALO)


      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      call ab_add(lstr, oWEIGHT, aWEIGHT, shape(oWEIGHT), 'ij')

      call ab_add(lstr, oDMSI, aDMSItmp, shape(oDMSI),
     &     'lij', oWEIGHT, aWEIGHT)
      call ab_add(lstr, oDHSI, aDHSItmp, shape(oDHSI),
     &     'lij', oWEIGHT, aWEIGHT) 
      call ab_add(lstr, oDSSI, aDSSItmp, shape(oDSSI),
     &     'lij', oWEIGHT, aWEIGHT) 

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)
      IF (ntmocn == ice%NTM) THEN
         allocate(otmp(ntmocn*2,ogrid%I_STRT_HALO:ogrid%I_STOP_HALO,
     &        ogrid%J_STRT_HALO:ogrid%J_STOP_HALO))
         allocate(atmp(ntmocn*2,agrid%I_STRT_HALO:agrid%I_STOP_HALO,
     &        agrid%J_STRT_HALO:agrid%J_STOP_HALO))
         otmp = 0.
         DO J=ogrid%J_STRT,oGRID%J_STOP
         DO I=ogrid%I_STRT,oGRID%I_STOP
           IF (oFOCEAN_loc(I,J).gt.0.) THEN
             otmp(:,I,J) = reshape(oDTRSI(:,:,I,J),(/2*ntmocn/))
           ELSE
             otmp(:,I,J) = 0.
           END IF
         ENDDO
         ENDDO
         call ab_add(lstr,otmp,atmp,shape(otmp),'lij',oWEIGHT,aWEIGHT)
      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif

c*   actual interpolation here
      call bundle_interpolation(lstr,remap_O2A,copy_np,do_np_avg)
c*

      do j=agrid%j_strt,agrid%j_stop
      do i=agrid%i_strt,agrid%i_stop
        if(ice%FWATER(i,j) > 0.) then
           ice%DMSI(:,i,j) = aDMSItmp(:,i,j)
           ice%DHSI(:,i,j) = aDHSItmp(:,i,j)
           ice%DSSI(:,i,j) = aDSSItmp(:,i,j)
        endif
      enddo
      enddo

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)
      IF (ntmocn == ice%NTM) THEN
        do j=agrid%j_strt,agrid%j_stop
        do i=agrid%i_strt,agrid%i_stop
          if(ice%FWATER(i,j) > 0.) then
            ice%DTRSI(:,:,i,j) = reshape(atmp(:,i,j),(/ntmocn,2/))
          endif
        enddo
        enddo
        deallocate(otmp,atmp)
      ELSE
      ENDIF
#endif

      deallocate(oWEIGHT, aWEIGHT)
      deallocate(aDMSItmp,aDHSItmp,aDSSItmp)

      RETURN
      END SUBROUTINE OG2AG_oceans

#else
      SUBROUTINE OG2AG_oceans(ice)
!@sum  OG2AG_oceans: ocean arrays for sea ice formation calculated in the
!       subr. OCEANS are gathered on the ocean grid, interpolated to the
!!      atmospheric grid, and scattered on the atmospheric ocean grid
!@auth Larissa Nazarenko

      USE OCEAN, only : oIM=>IM,oJM=>JM
     *                , oFOCEAN_loc=>FOCEAN

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only: tracerlist
#endif

      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : ocnice

      USE INT_OG2AG_MOD, only : INT_OG2AG

      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      IMPLICIT NONE
      type(iceocn_xchng_vars) :: ice

      REAL*8,
     * DIMENSION(oIM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)::
     * oWEIGHT

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: oDTR,aDTR
#endif
      INTEGER :: I,J,aIM
      INTEGER :: aJ_0,aJ_1,aJ_0H,aJ_1H, oJ_0,oJ_1,oJ_0H,oJ_1H
      integer :: ntmocn

#ifdef TRACERS_OCEAN
      ntmocn=tracerlist%getsize()
#endif

      aJ_0 = ice%J_0
      aJ_1 = ice%J_1
      aIM = ice%I_1

      aJ_0H = ice%J_0H
      aJ_1H = ice%J_1H

      call getDomainBounds(ogrid, J_STRT=oJ_0, J_STOP=oJ_1,
     &               J_STRT_HALO=oJ_0H, J_STOP_HALO=oJ_1H)

      oWEIGHT(:,:) = oFOCEAN_loc(:,:)
      CALL INT_OG2AG(ocnice%DMSI,ice%DMSI, oWEIGHT, 2, ice%FWATER)

      CALL INT_OG2AG(ocnice%DHSI,ice%DHSI, oWEIGHT, 2, ice%FWATER)

      CALL INT_OG2AG(ocnice%DSSI,ice%DSSI, oWEIGHT, 2, ice%FWATER)

#if (defined TRACERS_OCEAN) && (defined TRACERS_ON)

      IF (ntmocn == ice%NTM) THEN

        oWEIGHT(:,:) = oFOCEAN_loc(:,:)

        allocate(oDTR(ntmocn*2,oIM,oJ_0H:oJ_1H))
        allocate(aDTR(ntmocn*2,aIM,aJ_0H:aJ_1H))
        do j=oJ_0,oJ_1
        do i=1,oIM
          oDTR(:,i,j) = reshape(ocnice%DTRSI(1:ntmocn,1:2,i,j),
     .                                                 (/2*ntmocn/))
        enddo
        enddo
        CALL INT_OG2AG(oDTR,aDTR, oWEIGHT, 2*ntmocn, ice%FWATER)
        do j=aJ_0,aJ_1
        do i=1,aIM
          ice%DTRSI(:,:,i,j) = reshape(aDTR(:,i,j),(/ntmocn,2/))
        enddo
        enddo
        deallocate(oDTR,aDTR)

      ELSE
C**** if number of ocean and atm. tracers are not the same
C**** do something in here

      END IF
#endif
      RETURN
      END SUBROUTINE OG2AG_oceans

#endif /* OG2AG_OCEANS_BUNDLE */

      MODULE IGOG_regrid_info
!@sum IGOG_info saves instances of hntrp_type for use in
!@+   ice <-> ocean regrids
      use hntrp_mod, only : hntrp_type
      implicit none
      save

      type(hntrp_type) :: hntrp_i2o_u ! ice u C -> ocn u C
      type(hntrp_type) :: hntrp_i2o_v ! ice v C -> ocn v C
      logical :: hntrp_i2o_uv_need_init = .true.

      type(hntrp_type) :: hntrp_o2i_u ! ocn u C -> ice u A
      type(hntrp_type) :: hntrp_o2i_v ! ocn v C -> ice v A
      logical :: hntrp_o2i_uv_need_init = .true.

      END MODULE IGOG_regrid_info


      SUBROUTINE IG2OG_oceans(ice)
!@sum IG2OG_oceans interpolates relevant DYNSI outputs to the
!@+   ocean grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      use ofluxes,    only : ocnice
      USE OCEAN,  only : oIM=>im,oJM=>jm
      USE OCEAN,  only : oDLATM=>DLATM
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK,
     &      hasNorthPole, hasSouthPole
      USE EXCHANGE_TYPES, only : iceocn_xchng_vars
      implicit none
      type(iceocn_xchng_vars) :: ice
      real*8, dimension(:,:), allocatable ::
     &     ones_band,idmui_band,idmvi_band
      integer :: iIM,iJM,jmin,jmax

c assumption: every dynsi PE is an ocean PE
      if(.not. oGRID%have_domain) return

      iIM = ice%grid%im_world
      iJM = ice%grid%jm_world

c      if(iIM.eq.oIM .and. iJM.eq.oJM) then
c        ocnice%DMUI(:,:) = iDMUI(:,:)
c        ocnice%DMVI(:,:) = iDMVI(:,:)
c        return
c      endif

      if(hntrp_i2o_uv_need_init) then
        call Init_Hntrp_Type(hntrp_i2o_u,
     &     ice%GRID, .5d0,ice%DLATM,
     &     oGRID, .5d0,oDLATM,
     &     0.d0)
        call Init_Hntrp_Type(hntrp_i2o_v,
     &     ice%GRID, 0.d0,ice%DLATM,
     &     oGRID, 0.d0,oDLATM,
     &     0.d0,
     &     JMA_4interp=iJM-1,JMB_4interp=oJM-1) ! for secondary lats
        hntrp_i2o_uv_need_init = .false.
      endif

c regrid DMUI from ice C to ocn C
      jmin = hntrp_i2o_u%bpack%jband_strt
      jmax = hntrp_i2o_u%bpack%jband_stop
      ALLOCATE(iDMUI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_u%bpack, ice%DMUI, iDMUI_band)
      if(jmax == iJM) then
        iDMUI_band(:,iJM) = 0.
      endif
      if(iIM.eq.oIM .and. iJM.eq.oJM) then ! no interp necessary
        ocnice%DMUI(:,jmin:jmax) = iDMUI_band(:,jmin:jmax)
      else
        call HNTR8_band (ones_band, iDMUI_band, hntrp_i2o_u,ocnice%DMUI)
      endif
      deallocate(iDMUI_band,ones_band)

c regrid DMVI from ice C to ocn C
      jmin = hntrp_i2o_v%bpack%jband_strt
      jmax = hntrp_i2o_v%bpack%jband_stop
      ALLOCATE(iDMVI_band(iIM,jmin:jmax),ones_band(iIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_i2o_v%bpack, ice%DMVI, iDMVI_band)
      if(jmax == iJM) then
        iDMVI_band(:,iJM) = 0.
      endif
      if(iIM.eq.oIM .and. iJM.eq.oJM) then ! no interp necessary
        ocnice%DMVI(:,jmin:jmax) = iDMVI_band(:,jmin:jmax)
      else
        call HNTR8_band (ones_band, iDMVI_band, hntrp_i2o_v,ocnice%DMVI)
      endif
      deallocate(iDMVI_band,ones_band)

      if(hasNorthPole(oGRID)) then ! INT_AG2OG_Vector2 set them to zero
        ocnice%DMUI(:,oJM) = 0.0
        ocnice%DMVI(:,oJM) = 0.0
      endif

      return
      END SUBROUTINE IG2OG_oceans

      SUBROUTINE OG2IG_uvsurf(ice,atm)
!@sum OG2IG_uvsurf interpolates ocean surface velocity to the
!@+   DYNSI A-grid
!@auth M. Kelley
      use hntrp_mod
      use IGOG_regrid_info
      USE OCEAN,  only     : oIM=>im,oJM=>jm
      USE OCEAN,  only : oDLATM=>DLATM
      USE OCEANR_DIM, only : oGRID
      USE DOMAIN_DECOMP_1D, only : BAND_PACK, 
     &     hasNorthPole, hasSouthPole
      USE OCEAN, only : UO,VO
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      implicit none
      type(iceocn_xchng_vars) :: ice
      type(atmocn_xchng_vars) :: atm
c
      real*8, dimension(:,:), allocatable ::
     &     ones_band,ocnu_band,ocnv_band
      integer :: iIM,iJM,jmin,jmax

c assumption: every dynsi PE is an ocean PE
      if(.not. oGRID%have_domain) return

      iIM = ice%grid%im_world
      iJM = ice%grid%jm_world

#ifndef CUBED_SPHERE
c DYNSI grid == ATM grid, so simply replicate ATM copy
c (Note this requires that TOC2SST has been called first)
      call band_pack(ice%pack_a2i, atm%uosurf, ice%uosurf)
      call band_pack(ice%pack_a2i, atm%vosurf, ice%vosurf)
      return
#endif

c call HNTRP
      if(hntrp_o2i_uv_need_init) then
        call Init_Hntrp_Type(hntrp_o2i_u,
     &       oGRID, .5d0,oDLATM,
     &       ice%GRID, 0.d0,ice%DLATM,
     &       0.d0)
        call Init_Hntrp_Type(hntrp_o2i_v,
     &       oGRID, 0.d0,oDLATM,
     &       ice%GRID, 0.d0,ice%DLATM,
     &       0.d0,
     &       JMA_4interp=oJM-1) ! secondary ocean lats
        hntrp_o2i_uv_need_init = .false.
      endif
c regrid uosurf
      jmin = hntrp_o2i_u%bpack%jband_strt
      jmax = hntrp_o2i_u%bpack%jband_stop
      ALLOCATE(ocnu_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_o2i_u%bpack, uo(:,:,1), ocnu_band)
      call HNTR8_band (ones_band, ocnu_band, hntrp_o2i_u, ice%uosurf)
      deallocate(ocnu_band,ones_band)
c regrid vosurf
      jmin = hntrp_o2i_v%bpack%jband_strt
      jmax = hntrp_o2i_v%bpack%jband_stop
      ALLOCATE(ocnv_band(oIM,jmin:jmax),ones_band(oIM,jmin:jmax))
      ones_band(:,:) = 1d0
      call BAND_PACK (hntrp_o2i_v%bpack, vo(:,:,1), ocnv_band)
      call HNTR8_band (ones_band, ocnv_band, hntrp_o2i_v, ice%vosurf)
      deallocate(ocnv_band,ones_band)
      if(hasNorthPole(ice%grid)) then
c Eventually, uosurf/vosurf will be interpolated directly to
c the ice B grid, which has no polar point.  For now, setting
c A-grid values at the polar point to zero
        ice%uosurf(:,iJM) = 0.
        ice%vosurf(:,iJM) = 0.
      endif
      return
      END SUBROUTINE OG2IG_uvsurf
