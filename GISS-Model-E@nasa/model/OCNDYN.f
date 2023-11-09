#include "rundeck_opts.h"
#ifdef TRACERS_ATM_ONLY
#undef TRACERS_WATER
#endif

#ifndef ODIFF_FIXES_2017
#define ODIFF_FIXES_2017
#endif

!      SUBROUTINE OCEANS_old
!C****
!!@sum  OCEANS integrates ocean source terms and dynamics
!!@auth Gary Russell / Gavin Schmidt
!!@ver  2009/06/30
!C****
!      USE CONSTANT, only : rhows,grav
!      USE MODEL_COM, only : idacc,msurf
!#ifdef TRACERS_OceanBiology
!     * ,nstep=>itime
!#endif
!      USE DIAG_COM, only : modd5s
!      USE OCEANRES, only : NOCEAN,dZO
!      USE OCEAN, only : im,jm,lmo,ndyno,mo,g0m,gxmo,gymo,gzmo,
!     *    s0m,sxmo,symo,szmo,dts,dtofs,dto,dtolf,mdyno,msgso,
!     *    ogeoz,ogeoz_sv,opbot,ze,lmm,imaxj, UO,VONP,IVNP, ! VOSP,IVSP,
!     *    OBottom_drag,OCoastal_drag,focean,
!     *    G0M_glob,GXMO_glob,GYMO_glob,GZMO_glob,
!     *    S0M_glob,SXMO_glob,SYMO_glob,SZMO_glob,
!     *    scatter_ocean, gather_ocean
!#ifdef TRACERS_OCEAN
!     *    ,trmo,txmo,tymo,tzmo
!     *    ,trmo_glob,txmo_glob,tymo_glob,tzmo_glob
!#endif
!      USE OCEAN_DYN, only : mmi,smu,smv,smw
!      USE DOMAIN_DECOMP_1D, only : get, AM_I_ROOT
!      USE OCEANR_DIM, only : grid=>ogrid
!      USE ODIAG, only : oijl=>oijl_loc,oij=>oij_loc,
!     *    ijl_mo,ijl_g0m,ijl_s0m, ijl_gflx, ijl_sflx, ijl_mfw2,
!     *    ijl_mfu,ijl_mfv,ijl_mfw, ijl_ggmfl,ijl_sgmfl,ij_ssh,ij_pb
!#ifdef TRACERS_OCEAN
!     *    ,toijl=>toijl_loc,
!     *     toijl_conc,toijl_tflx,toijl_gmfl
!#endif
!
!#ifdef TRACERS_OCEAN
!      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
!#endif
!
!      IMPLICIT NONE
!      Integer*4 I,J,L,N,NS,NO  ; real*8 now
!      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
!     *       MM0,MM1, UM0,UM1, VM0,VM1
!#ifdef TRACERS_OCEAN
!      type(ocn_tracer_entry), pointer :: entry
!#endif
!
!c**** Extract domain decomposition info
!      INTEGER :: J_0, J_1, J_0H
!      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1, J_STRT_HALO = J_0H)
!
!C***  Get the data from the atmospheric grid to the ocean grid
!      call AG2OG_oceans
!      OGEOZ_SV(:,:)=OGEOZ(:,:)
!
!C***  Interpolate DYNSI outputs to the ocean grid
!C***  (at present, only the ice-ocean stress is used)
!      call IG2OG_oceans
!
!C**** Apply surface fluxes to ocean
!      CALL GROUND_OC
!         CALL CHECKO('GRNDOC')
!
!C**** Apply ice/ocean and air/ocean stress to ocean
!      CALL OSTRES
!         CALL CHECKO('OSTRES')
!         CALL TIMER (NOW,MSURF)
!         IF (MODD5S == 0) CALL DIAGCO (11)
!
!C**** Apply ocean vertical mixing
!      CALL OCONV
!         CALL CHECKO('OCONV ')
!
!C**** Apply bottom and coastal drags
!      if (OBottom_drag  == 1) CALL OBDRAG
!      if (OCoastal_drag == 1) CALL OCOAST
!
!C**** Add ocean biology
!#ifdef TRACERS_OceanBiology
!      call obio_model
!      IF (MODD5S.EQ.0) CALL DIAGCO (5)
!#endif
!
!         CALL TIMER (NOW,MSGSO)
!
!C****
!C**** Integrate Ocean Dynamics
!C****
!      Do 500 NO=1,NOCEAN
!
!C**** initialize summed mass fluxes
!      DO L=1,LMO
!        SMU(:,J_0:J_1,L) = 0
!        SMV(:,J_0:J_1,L) = 0
!        SMW(:,J_0:J_1,L) = 0  ;  EndDo
!
!      CALL OVtoM (MMI,UM0,VM0)
!      CALL OPGF0
!
!      NS = NDYNO / NOCEAN
!C**** Initial Forward step,  QMX = QM0 + DT*F(Q0)
!      Call OFLUX  (NS,MMI,.FALSE.)
!      Call OADVM  (MM1,MMI,DTOFS)
!C     Call STADVM (MM1,MUST,DTOFS,.False.)
!      Call OADVV  (UM1,VM1,UM0,VM0,DTOFS)
!C     Call OTIDEW (NS,DTOFS)
!      Call OPGF   (UM1,VM1,DTOFS)
!C     Call STPGF  (MUST1,MUST,DTOFS)
!C     Call OTIDEV (MM1,UM1,VM1,DTOFS)
!      Call OMtoV  (MM1,UM1,VM1)
!C**** Initial Backward step,  QM1 = QM0 + DT*F(Q0)
!      Call OFLUX  (NS,MMI,.FALSE.)
!      Call OADVM  (MM1,MMI,DTO)
!C     Call STADVM (MM1,MUST,DTO,.False.)
!      Call OADVV  (UM1,VM1,UM0,VM0,DTO)
!C     Call OTIDEW (NS,DTO)
!      Call OPGF   (UM1,VM1,DTO)
!C     Call STPGF  (MUST1,MUST,DTO)
!C     Call OTIDEV (MM1,UM1,VM1,DTO)
!      Call OMtoV  (MM1,UM1,VM1)
!C**** First even leap frog step,  Q2 = Q0 + 2*DT*F(Q1)
!      Call OFLUX  (NS,MMI,.TRUE.)
!      Call OADVM  (MM0,MMI,DTOLF)
!C     Call STADVM (MM0,MUST1,DTOLF,.True.)
!      Call OADVV  (UM0,VM0,UM0,VM0,DTOLF)
!C     Call OTIDEW (NS,DTOLF)
!      Call OPGF   (UM0,VM0,DTOLF)
!C     Call STPGF  (MUST,MUST,DTOLF)
!C     Call OTIDEV (MM0,UM0,VM0,DTOLF)
!      Call OMtoV  (MM0,UM0,VM0)
!      NS=NS-1
!C**** Odd leap frog step,  Q3 = Q1 + 2*DT*F(Q2)
!  420 Continue
!      Call OFLUX  (NS,MM1,.FALSE.)
!      Call OADVM  (MM1,MM1,DTOLF)
!C     Call STADVM (MM1,MUST,DTOLF,.False.)
!      Call OADVV  (UM1,VM1,UM1,VM1,DTOLF)
!C     Call OTIDEW (NS,DTOLF)
!      Call OPGF   (UM1,VM1,DTOLF)
!C     Call STPGF  (MUST1,MUST1,DTOLF)
!C     Call OTIDEV (MM1,UM1,VM1,DTOLF)
!      Call OMtoV  (MM1,UM1,VM1)
!      NS=NS-1
!C**** Even leap frog step,  Q4 = Q2 + 2*DT*F(Q3)
!      Call OFLUX  (NS,MM0,.TRUE.)
!      Call OADVM  (MM0,MM0,DTOLF)
!C     Call STADVM (MM0,MUST1,DTOLF,.True.)
!      Call OADVV  (UM0,VM0,UM0,VM0,DTOLF)
!C     Call OTIDEW (NS,DTOLF)
!      Call OPGF   (UM0,VM0,DTOLF)
!C     Call STPGF  (MUST,MUST,DTOLF)
!C     Call OTIDEV (MM0,UM0,VM0,DTOLF)
!      Call OMtoV  (MM0,UM0,VM0)
!C**** Check for end of leap frog time scheme
!      NS=NS-1
!      IF(NS.GT.1)  GO TO 420
!C     if(j_0 ==  1) UO(IVSP,1 ,:) = VOSP(:) ! not needed if Mod(IM,4)=0
!      if(j_1 == JM) UO(IVNP,JM,:) = VONP(:) ! not needed if Mod(IM,4)=0
!      DO L=1,LMO
!        OIJL(:,:,L,IJL_MFU) = OIJL(:,:,L,IJL_MFU) + SMU(:,:,L)
!        OIJL(:,:,L,IJL_MFV) = OIJL(:,:,L,IJL_MFV) + SMV(:,:,L)
!        OIJL(:,:,L,IJL_MFW) = OIJL(:,:,L,IJL_MFW) + SMW(:,:,L)
!        OIJL(:,:,L,IJL_MFW2)= OIJL(:,:,L,IJL_MFW2)+SMW(:,:,L)*SMW(:,:,L)
!      END DO
!C**** Advection of Potential Enthalpy and Salt
!      CALL OADVT (G0M,GXMO,GYMO,GZMO,DTOLF,.FALSE.
!     *        ,OIJL(1,J_0H,1,IJL_GFLX))
!      CALL OADVT (S0M,SXMO,SYMO,SZMO,DTOLF,.TRUE.
!     *        ,OIJL(1,J_0H,1,IJL_SFLX))
!
!#ifdef TRACERS_OCEAN
!      DO N=1,tracerlist%getsize() 
!        entry=>tracerlist%at(n)
!        CALL OADVT(TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N)
!     *       ,TYMO(1,J_0H,1,N),TZMO(1,J_0H,1,N),DTOLF,entry%t_qlimit
!     *       ,TOIJL(1,J_0H,1,TOIJL_TFLX,N))
!      END DO
!#endif
!
!        CALL CHECKO ('OADVT ')
!        DO L=1,LMO
!          OIJL(:,:,L,IJL_MO)  = OIJL(:,:,L,IJL_MO) +  MO(:,:,L)
!          OIJL(:,:,L,IJL_G0M) = OIJL(:,:,L,IJL_G0M) + G0M(:,:,L)
!          OIJL(:,:,L,IJL_S0M) = OIJL(:,:,L,IJL_S0M) + S0M(:,:,L)
!        END DO
!        DO J=J_0,J_1
!          DO I=1,IMAXJ(J)
!            IF (FOCEAN(I,J).gt.0) THEN
!              OIJ(I,J,IJ_SSH) = OIJ(I,J,IJ_SSH) + OGEOZ(I,J)
!              OIJ(I,J,IJ_PB)  = OIJ(I,J,IJ_PB)  +
!     +             (OPBOT(I,J)-ZE(LMM(I,J))*RHOWS*GRAV)
!            END IF
!c     do L=1,LMO
!c     write(*,'(a,5i5,3e12.4)')'for samar, ocndyn, ze',
!c    .  nstep,i,j,l,lmm(i,j),ze(l),dzo(l),ZE(LMM(I,J))
!c     enddo
!          END DO
!        END DO
!
!#ifdef TRACERS_OCEAN
!        DO N=1,tracerlist%getsize()
!          DO L=1,LMO
!            TOIJL(:,:,L,TOIJL_CONC,N)=TOIJL(:,:,L,TOIJL_CONC,N)
!     *           +TRMO(:,:,L,N)
!          END DO
!        END DO
!#endif
!
!      call gather_ocean_straits()
!
!      IF(AM_I_ROOT()) THEN
!C****
!C**** Acceleration and advection of tracers through ocean straits
!C****
!        CALL STPGF(DTS/NOCEAN)
!        CALL STADV(DTS/NOCEAN)
!          CALL CHECKO_serial ('STADV0')
!        IF (NO .EQ. NOCEAN) THEN
!          CALL STCONV
!          CALL STBDRA
!        END IF
!      END IF
!      call scatter_ocean_straits()
!      call BCAST_straits (.false.)
!        CALL CHECKO ('STADV ')
!
!  500 Continue  !  End of Do-loop NO=1,NOCEAN
!
!        CALL TIMER (NOW,MDYNO)
!        IF (MODD5S == 0) CALL DIAGCO (12)
!
!C**** Apply Wajsowicz horizontal diffusion to UO and VO ocean currents
!      CALL ODIFF(DTS)
!      CALL OABFILx ! binary filter
!      CALL OABFILy ! binary filter
!      CALL CHECKO ('ODIFF0')
!
!C**** Apply GM + Redi tracer fluxes
!      CALL GMKDIF
!      CALL GMFEXP(G0M,GXMO,GYMO,GZMO,.FALSE.,OIJL(1,J_0H,1,IJL_GGMFL))
!      CALL GMFEXP(S0M,SXMO,SYMO,SZMO,.TRUE. ,OIJL(1,J_0H,1,IJL_SGMFL))
!#ifdef TRACERS_OCEAN
!      DO N = 1,tracerlist%getsize()
!        entry=>tracerlist%at(n)
!        CALL GMFEXP(TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N),TYMO(1,J_0H,1,N),
!     *    TZMO(1,J_0H,1,N),entry%t_qlimit,TOIJL(1,J_0H,1,TOIJL_GMFL,N))
!      END DO
!#endif
!      CALL CHECKO ('GMDIFF')
!      CALL TIMER (NOW,MSGSO)
!
!C**** remove STADVI since it is not really consistent with ICEDYN
!c      CALL STADVI
!c        CALL CHECKO ('STADVI')
!
!#ifdef TRACERS_OCEAN
!      CALL OC_TDECAY(DTS)
!#ifdef TRACERS_AGE_OCEAN
!      CALL OCN_TR_AGE(DTS)
!#endif
!#endif
!
!        CALL TIMER (NOW,MSGSO)
!C***  Get the data from the ocean grid to the atmospheric grid
!      CALL TOC2SST
!      call OG2AG_oceans
!
!C***  Interpolate ocean surface velocity to the DYNSI grid
!      call OG2IG_uvsurf
!
!      RETURN
!      END SUBROUTINE OCEANS_old

      SUBROUTINE init_OCEAN(iniOCEAN,istart,atmocn,dynsice)
!@sum init_OCEAN initializes ocean variables
!@auth Original Development Team
      USE FILEMANAGER, only : openunit,closeunit
      USE Dictionary_mod
      USE CONSTANT, only : twopi,radius,by3,grav,rhow
      USE MODEL_COM, only : dtsrc,kocean
      USE OCEAN, only : im,jm,lmo,focean,lmm
     *     ,lmu,lmv,hatmo,hocean,ze,mo,g0m,s0m,zmid,dzoe,bydzoe
     *     ,uo,vo,uod,vod,dxypo,ogeoz,kpl
     *     ,dts,dtolf,dto,dtofs,nocean,mdyno,msgso
     *     ,ndyno,imaxj,ogeoz_sv,bydts,lmo_min,j1o
     *     ,OBottom_drag,OCoastal_drag,OTIDE,oc_salt_mean
      USE OCEAN, only : use_qus,
     *     GXMO,GYMO,GZMO, GXXMO,GYYMO,GZZMO, GXYMO,GYZMO,GZXMO,
     *     SXMO,SYMO,SZMO, SXXMO,SYYMO,SZZMO, SXYMO,SYZMO,SZXMO
      USE OCEAN, only : nbyzmax,
     &     nbyzm,nbyzu,nbyzv,nbyzc,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc
      USE OCEANRES, only : dZO
      USE OCFUNC, only : vgsp,tgsp,hgsp,agsp,bgsp,cgs
      USE SW2OCEAN, only : init_solar
      USE DOMAIN_DECOMP_1D, only : getDomainBounds,halo_update
      use DOMAIN_DECOMP_1D, only: hasSouthPole, hasNorthPole
      use DOMAIN_DECOMP_1D, only: pack_data,broadcast,globalmin
      USE OCEANR_DIM, only : grid=>ogrid
      USE OCEAN, only : remap_a2o,remap_o2a
#ifdef CUBED_SPHERE
      use cs2ll_utils, only : init_xgridremap_type
      use regrid_com, only : xA2O_root,
     &     read_xgrid_file=>init_regrid_root
#else
      use hntrp_mod, only : init_hntrp_type
      USE OCEAN, only : oDLATM=>DLATM
#endif
#ifdef TRACERS_OCEAN
      use model_com, only : nday,nssw,ndisk,itimee
      use ocean, only : ntrtrans,motr
      Use OCEAN, Only: oc_tracer_mean
      Use OCN_TRACER_COM, Only: tracerlist, ocn_tracer_entry
#endif
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
#ifdef OCN_GISS_TURB
      USE GISS_OTURB, only : gissmix_init
#endif
#ifdef OCN_GISS_SM
      USE GISS_SM, only : giss_sm_init
#endif
      use pario, only : par_open,par_close,read_dist_data,
     &     variable_exists
      IMPLICIT NONE
c
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER, INTENT(IN) :: istart
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: dynsice
c
      INTEGER I,J,L,N,fid,iu_OFTAB,IP1,IM1,LMIJ,I1,J1,I2,J2
     *     ,II,JJ,flagij,nt,j1o_loc
      CHARACTER*80 TITLE
      REAL*8 FJEQ,SM,SG0,SGZ,SS0,SSZ
      LOGICAL :: postProc
      logical :: qexist(im)
#ifdef CUBED_SPHERE
      integer, allocatable :: ones(:)
#endif
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_NORTH_POLE

      INTEGER, DIMENSION(IM,JM) :: LMM_glob

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1
     *      ,J_STRT_SKP  = J_0S, J_STOP_SKP  = J_1S
     *      ,J_STRT_HALO  = J_0H, J_STOP_HALO  = J_1H
     *      ,HAVE_NORTH_POLE = HAVE_NORTH_POLE)

!     soon obsolete: postProcessing case
      postProc = .false. ; if(istart < 1) postProc = .true.

      if (LMO_MIN .lt. 2) then
        write (*,*) ' Make minimum number of ocean layers equal to 2'
        stop
      end if

C****
C**** Check that KOCEAN is set correctly
C****
      IF (KOCEAN == 0) THEN
        call stop_model(
     &       "Must have KOCEAN > 0 for interactive ocean runs",255)
      END IF
C****
C**** Select drag options and ocean tides
C****
      call sync_param("OBottom_drag",OBottom_drag)
      call sync_param("OCoastal_drag",OCoastal_drag)
      Call SYNC_PARAM ("OTIDE",OTIDE)

C**** define initial condition options for global mean
      call sync_param("oc_salt_mean",oc_salt_mean)
#ifdef TRACERS_OCEAN
      call sync_param("oc_tracer_mean",oc_tracer_mean,
     &            tracerlist%getsize())
#endif

      call alloc_ofluxes(atmocn)     ! should be moved back to alloc_ocean

C****
C**** set up time steps from atmospheric model
C****
      call sync_param("DTO",DTO)
      call sync_param('nocean',nocean)

      DTS=DTSRC
      BYDTS=1d0/DTS
      NDYNO=2*NINT(.5*DTS/DTO)
      DTO=DTS/NDYNO
      DTOLF=2.*DTO
      DTOFS=2.*DTO*BY3
C**** Set up timing indexes
      CALL SET_TIMER(" OCEAN DYNAM",MDYNO)
      CALL SET_TIMER(" OCEAN PHYS.",MSGSO)
C****
C**** Arrays needed each ocean model run
C****

C**** Calculate ZE
      ZE(0) = 0d0
      DO L = 1,LMO
        ZE(L) = ZE(L-1) + dZO(L)
        ZMID(L) = .5D0*(ZE(L)+ZE(L-1))
      END DO

C**** Calculate DZOE
      DZOE(0) = 0d0
      DO L = 1,LMO-1
        DZOE(L) = .5*(DZO(L)+DZO(L+1))
        BYDZOE(L) = 1D0/DZOE(L)
      END DO

C**** Read in table function for specific volume
      CALL openunit("OFTAB",iu_OFTAB,.TRUE.,.TRUE.)
      READ  (iu_OFTAB) TITLE,VGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,TGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,CGS
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,HGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,AGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      READ  (iu_OFTAB) TITLE,BGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      call closeunit(iu_OFTAB)

c-------------------------------------------------------------------
c Begin ocean-processors-only code region
      ocean_processors_only: if(grid%have_domain) then
c-------------------------------------------------------------------

      CALL OFFT0(IM)
      If (OTIDE > 0)  Call OTIDE0

C**** Calculate J1O = least J with some ocean
      j1o_loc = huge(j1o_loc)
      do j=j_0,j_1
        if(sum(focean(:,j)) > 0d0) then
          j1o_loc = j
          exit
        endif
      enddo
      call globalmin(grid,j1o_loc,j1o)
      write(6,*) "Minimum J with some ocean:",J1O
C**** Fix to 4 for temporary consistency
      J1O=4
      write(6,*) "Fixed J1O:",J1O

C**** Calculate LMM and modify HOCEAN
c      call sync_param("LMO_min",LMO_min)
      do j=max(1,j_0h),min(j_1h,jm)
        do i=1,im
          lmm(i,j) = 0
          if(focean(i,j).le.0.) cycle
          do l=lmo_min,lmo-1
            if(hatmo(i,j)+hocean(i,j) .le. 5d-1*(ze(l)+ze(l+1))) exit
          enddo
          lmm(i,j)=l
          hocean(i,j) = -hatmo(i,j) + ze(l)
        enddo
      enddo
 
C**** Calculate LMU
      do j=max(1,j_0h),min(j_1h,jm)
        i=im
        do ip1=1,im
          lmu(i,j) = min(lmm(i,j),lmm(ip1,j))
          i=ip1
        enddo
      enddo


C**** Calculate LMV
      call pack_data(grid,lmm,lmm_glob)
      call broadcast(grid,lmm_glob)
      do j=max(1,j_0h),min(j_1h,jm-1)
        do i=1,im
          lmv(i,j) = min(lmm_glob(i,j),lmm_glob(i,j+1))
        enddo
      enddo

C****
C**** Tabulate basin start/end indices
C****
      do l=1,lmo
        do j=max(2,j_0h),min(jm-1,j_1h)
          qexist(:) = (l <= lmm(:,j))
          call get_i1i2(qexist,im,nbyzm(j,l),
     &         i1yzm(1,j,l),i2yzm(1,j,l),nbyzmax)
          qexist(:) = (l <= lmu(:,j))
          call get_i1i2(qexist,im,nbyzu(j,l),
     &         i1yzu(1,j,l),i2yzu(1,j,l),nbyzmax)
          qexist(:) = (l <= lmv(:,j))
          call get_i1i2(qexist,im,nbyzv(j,l),
     &         i1yzv(1,j,l),i2yzv(1,j,l),nbyzmax)
        enddo
        do j=max(2,j_0h),min(jm-1,j_1)
          do i=1,im-1
            qexist(i) =
     &           (l <= lmv(i,j)) .or. (l <= lmv(i+1,j)) .or.
     &           (l <= lmu(i,j)) .or. (l <= lmu(i,j+1))
          enddo
          i=im
          qexist(i) =
     &           (l <= lmv(i,j)) .or. (l <= lmv(  1,j)) .or.
     &           (l <= lmu(i,j)) .or. (l <= lmu(i,j+1))
          call get_i1i2(qexist,im,nbyzc(j,l),
     &         i1yzc(1,j,l),i2yzc(1,j,l),nbyzmax)
        enddo
        if(hasSouthPole(grid)) then
          nbyzm(1,l) = 0
          nbyzu(1,l) = 0
          nbyzv(1,l) = 0
          nbyzc(1,l) = 0
        endif
        if(hasNorthPole(grid)) then
          nbyzm(jm,l) = 0
          nbyzu(jm,l) = 0
          nbyzv(jm,l) = 0
          if(l <= lmm(1,jm)) nbyzm(jm,l) = 1
          i1yzm(1,jm,l) = 1
          i2yzm(1,jm,l) = 1
        endif
      enddo

      IF(iniOCEAN) THEN

C**** Initialize a run from ocean initial conditions

      fid = par_open(grid,'OIC','read')

      native_oic: if(variable_exists(grid,fid,'mo')) then
        ! oic with full set of native variables
        call read_dist_data(grid,fid,'mo',mo)
        call read_dist_data(grid,fid,'g' ,g0m)
        call read_dist_data(grid,fid,'gz',gzmo)
        call read_dist_data(grid,fid,'s' ,s0m)
        call read_dist_data(grid,fid,'sz',szmo)

        call halo_update(grid,mo)
        call halo_update(grid,g0m)
        call halo_update(grid,gzmo)
        call halo_update(grid,s0m)
        call halo_update(grid,szmo)

C**** Calculate layer mass from column mass and check for consistency
        DO J=J_0,J_1
        DO I=1,IM
          LMIJ=LMM(I,J)
          DO L=LMIJ+1,LMO
            MO(I,J,L) = 0.
          enddo
C**** if there is a problem try nearest neighbour
          IF((LMM(I,J).GT.0).AND.(ABS(MO(I,J,1)/ZE(1)-1d3).GT.5d1)) THEN
            WRITE (6,931) I,J,LMIJ,MO(I,J,1),ZE(1)
            II=0 ; JJ=0 ; flagij=0
            do j1=0,2
              do i1=0,4
                if(flagij.eq.0) then
                  if(j1.eq.0) jj=j
                  if(j1.eq.1) jj=j-1
                  if(j1.eq.2) jj=j+1

                  if(i1.eq.0) ii=i
                  if(i1.eq.1) ii=i-1
                  if(i1.eq.2) ii=i+1
                  if(i1.eq.3) ii=i-2
                  if(i1.eq.4) ii=i+2
                  if(i1.eq.5) ii=i-3
                  if(i1.eq.6) ii=i+3
                  if(ii.gt.im) ii=ii-im
                  if(ii.lt.1) ii=ii+im
                  if(jj.gt.jm) then
                    jj=jm
                    if(ii.le.im/2) ii=ii+im/2
                    if(ii.gt.im/2) ii=ii-im/2
                  endif
                  if(jj.lt.1) then
                    jj=1
                    if(ii.le.im/2) ii=ii+im/2
                    if(ii.gt.im/2) ii=ii-im/2
                  endif
                  IF ((MO(II,JJ,1).gt.0) .and. (LMM(II,JJ).ge.LMM(I,J)))
     *                 flagij=1
                endif
              enddo
            enddo
            IF (flagij.ne.0) THEN
              MO(I,J,1:LMM(I,J))=MO(II,JJ,1:LMM(I,J))
              G0M(I,J,1:LMM(I,J))=G0M(II,JJ,1:LMM(I,J))
              S0M(I,J,1:LMM(I,J))=S0M(II,JJ,1:LMM(I,J))
              GZMO(I,J,1:LMM(I,J))=GZMO(II,JJ,1:LMM(I,J))
              SZMO(I,J,1:LMM(I,J))=SZMO(II,JJ,1:LMM(I,J))
              WRITE (6,*) "Inconsistency at ",I,J,"fixed from :",II,JJ
            END IF
          END IF
        enddo
        enddo
      else ! not native oic
        call tempsalt_oic(fid,mo,g0m,gzmo,s0m,szmo)
        call halo_update(grid,mo)
        call halo_update(grid,g0m)
        call halo_update(grid,gzmo)
        call halo_update(grid,s0m)
        call halo_update(grid,szmo)
      endif native_oic

      call par_close(grid,fid)

C**** Initialize velocity field and slopes of pot. heat and salinity to
C**** zero
      UO=0
      VO=0
      UOD=0
      VOD=0
      GXMO=0
      GYMO=0
      SXMO=0
      SYMO=0
C**** Define mean value of mass, potential heat, and salinity at poles
      DO 370 L=1,LMO
      if(HAVE_NORTH_POLE) then ! average polar ocean fields
        J=JM
        SM  = 0.
        SG0 = 0.
        SGZ = 0.
        SS0 = 0.
        SSZ = 0.
        DO I=1,IM
          SM  = SM  +   MO(I,J,L)
          SG0 = SG0 + G0M(I,J,L)
          SGZ = SGZ + GZMO(I,J,L)
          SS0 = SS0 + S0M(I,J,L)
          SSZ = SSZ + SZMO(I,J,L)
        end do
        DO I=1,IM
          MO(I,J,L)   = SM /IM
          G0M(I,J,L) = SG0/IM
          GZMO(I,J,L) = SGZ/IM
          S0M(I,J,L) = SS0/IM
          SZMO(I,J,L) = SSZ/IM
        end do
      end if
C**** Define East-West horizontal gradients
      IM1=IM-1
      I=IM
      DO 345 J=J_0S,J_1S
      DO 345 IP1=1,IM
      IF(LMM(I  ,J).LT.L)  GO TO 344
      IF(LMM(IM1,J).GE.L)  GO TO 342
      IF(LMM(IP1,J).LT.L)  GO TO 344
      GXMO(I,J,L) = .5*(G0M(IP1,J,L)-G0M(I,J,L))
      SXMO(I,J,L) = .5*(S0M(IP1,J,L)-S0M(I,J,L))
      GO TO 344
  342 IF(LMM(IP1,J).GE.L)  GO TO 343
      GXMO(I,J,L) = .5*(G0M(I,J,L)-G0M(IM1,J,L))
      SXMO(I,J,L) = .5*(S0M(I,J,L)-S0M(IM1,J,L))
      GO TO 344
  343 GXMO(I,J,L) = .25*(G0M(IP1,J,L)-G0M(IM1,J,L))
      SXMO(I,J,L) = .25*(S0M(IP1,J,L)-S0M(IM1,J,L))
  344 IM1=I
  345 I=IP1
C**** Define North-South horizontal gradients
      DO 354 J=J_0S,J_1S
      DO 354 I=1,IM
      IF(LMM(I,J  ).LT.L)  GO TO 354
      IF(LMM(I,J-1).GE.L)  GO TO 352
      IF(LMM(I,J+1).LT.L)  GO TO 354
      GYMO(I,J,L) = .5*(G0M(I,J+1,L)-G0M(I,J,L))
      SYMO(I,J,L) = .5*(S0M(I,J+1,L)-S0M(I,J,L))
      GO TO 354
  352 IF(LMM(I,J+1).GE.L)  GO TO 353
      GYMO(I,J,L) = .5*(G0M(I,J,L)-G0M(I,J-1,L))
      SYMO(I,J,L) = .5*(S0M(I,J,L)-S0M(I,J-1,L))
      GO TO 354
  353 GYMO(I,J,L) = .25*(G0M(I,J+1,L)-G0M(I,J-1,L))
      SYMO(I,J,L) = .25*(S0M(I,J+1,L)-S0M(I,J-1,L))
  354 CONTINUE
C**** Multiply specific quantities by mass
      DO 360 J=J_0,J_1
      DO 360 I=1,IM
      G0M(I,J,L)  = G0M(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GXMO(I,J,L) = GXMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GYMO(I,J,L) = GYMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      GZMO(I,J,L) = GZMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      S0M(I,J,L)  = S0M(I,J,L)*(MO(I,J,L)*DXYPO(J))
      SXMO(I,J,L) = SXMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
      SYMO(I,J,L) = SYMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
  360 SZMO(I,J,L) = SZMO(I,J,L)*(MO(I,J,L)*DXYPO(J))
  370 CONTINUE
C**** Initiallise geopotential field (needed by KPP)
      OGEOZ = 0.
      OGEOZ_SV = 0.
C**** Initiallise kpl
      kpl=3

      END IF ! end if(iniOCEAN)

C**** Extend ocean data to added layers at bottom if necessary
      if (.not.postProc .and. lmo_min .gt. 1) then
        do j=j_0s,j_1s
        do i=1,im
          if (lmm(i,j) == lmo_min .and. MO(i,j,lmo_min) == 0.) then
            do l=2,lmo_min
              if (MO(i,j,l) == 0.) then
                MO(i,j,l)   = (ZE(L)-ZE(L-1))*RHOW*
     *                     (1.+S0M(i,j,l-1)/(MO(i,j,l-1)*DXYPO(J)))
                G0M(i,j,l)  = G0M(i,j,l-1) *(MO(i,j,l)/MO(i,j,l-1))
                GXMO(i,j,l) = GXMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                GYMO(i,j,l) = GYMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                GZMO(i,j,l) = 0.
                S0M(i,j,l)  = S0M(i,j,l-1) *(MO(i,j,l)/MO(i,j,l-1))
                SXMO(i,j,l) = SXMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                SYMO(i,j,l) = SYMO(i,j,l-1)*(MO(i,j,l)/MO(i,j,l-1))
                SZMO(i,j,l) = 0.
              end if
            end do
          end if
        end do
        end do
      end if

C**** zero out unphysical values (that might have come from a
C**** restart file with different topography)
      if (iniOCEAN) then
        DO J=J_0S,J_1S
          DO I=1,IMAXJ(J)
            VO(I,J,LMV(I,J)+1:LMO)=0.
          END DO
        END DO
        DO J=J_0,J_1
          DO I=1,IMAXJ(J)
            UO(I,J,LMU(I,J)+1:LMO)=0.
            MO(I,J,LMM(I,J)+1:LMO)=0.
          END DO
        END DO
      end if

#ifdef TRACERS_OCEAN
      motr(:,:,:) = mo(:,:,:)
#endif

!      if(istart.eq.2 .and. use_qus.eq.1) then
!        allocate(zero3d(im,j_0h:j_1h,lmo))
!        allocate(marr(im,j_0h:j_1h,lmo))
!        do l=1,lmo
!          do j=j_0,j_1
!            do n=1,nbyzm(j,l)
!              do i=i1yzm(n,j,l),i2yzm(n,j,l)
!                MARR(I,J,L) = MO(I,J,L)*DXYPO(J)
!              enddo
!            enddo
!          enddo
!        enddo
!        zero3d = 0d0
!        gzmo = 0.; gzzmo = 0.
!        szmo = 0.; szzmo = 0.
!        call relax_zmoms(marr,g0m,zero3d,gzmo,gzzmo)
!        call relax_zmoms(marr,s0m,zero3d,szmo,szzmo)
!        deallocate(zero3d,marr)
!      endif

C**** Initialize straits arrays
      call init_STRAITS(iniOCEAN)

C**** Adjust global mean salinity if required (only at first start up)
      if (iniOCEAN .and. oc_salt_mean.ne.-999.)
     *     call adjust_mean_salt

C**** Initialize solar radiation penetration arrays
      call init_solar

C**** Initialize KPP mixing scheme
      call alloc_kpp_com(grid) ! alloc moved here after lsrpd is set
      call kmixinit(ZE)

#ifdef OCN_GISS_TURB
C**** Initialize GISS mixing scheme
      call gissmix_init(iniOCEAN)
#endif
#ifdef OCN_GISS_SM
C**** Initialize GISS SM mixing scheme
      call giss_sm_init(iniOCEAN)
#endif

C***  Initialize ODIFF
      call init_ODIFF(grid)

#ifdef TRACERS_OCEAN
      call tracer_ic_ocean(atmocn)
#endif

c-------------------------------------------------------------------
c End ocean-processors-only code region
      endif ocean_processors_only
c-------------------------------------------------------------------

#ifdef TRACERS_OCEAN
      ! Sanity checks on tracer transport timestep.
      if(mod(nday,ntrtrans).ne.0) then
        call stop_model('mod(nday,ntrtrans).ne.0',255)
      endif
      ! The following requirement avoids the need to save
      ! partial mass flux accumulations in the restart file.
      if(mod(nssw,ntrtrans).ne.0) then
        call stop_model('mod(nssw,ntrtrans).ne.0',255)
      endif

!@PL
!      if (mod(itimee,ntrtrans).ne.0) then
!       call stop_model('houre*2 must be a multiple of ntrtrans',255)
!      endif

      if (mod(ndisk,ntrtrans).ne.0) then
       call stop_model('ndisk must be a multiple of ntrtrans',255)
      endif
#endif

#ifdef CUBED_SPHERE
      call read_xgrid_file(xA2O_root,
     &     atmocn%grid%im_world,atmocn%grid%jm_world,6,
     &     im,                 jm,                 1)
c*** fill in vector full of ones
      allocate(ones(xA2O_root%xgridroot%ncells))
      ones=1
c*** initialize remapping derived types
      call init_xgridremap_type(atmocn%grid,grid,
     &     xA2O_root%xgridroot%ncells,
     &     xA2O_root%xgridroot%ijcub(1,:),
     &     xA2O_root%xgridroot%ijcub(2,:),
     &     xA2O_root%xgridroot%tile,
     &     xA2O_root%xgridroot%ijlatlon(1,:),
     &     xA2O_root%xgridroot%ijlatlon(2,:),
     &     ones,
     &     xA2O_root%xgridroot%xgrid_area,remap_a2o)

      call init_xgridremap_type(grid,atmocn%grid,
     &     xA2O_root%xgridroot%ncells,
     &     xA2O_root%xgridroot%ijlatlon(1,:),
     &     xA2O_root%xgridroot%ijlatlon(2,:),
     &     ones,
     &     xA2O_root%xgridroot%ijcub(1,:),
     &     xA2O_root%xgridroot%ijcub(2,:),
     &     xA2O_root%xgridroot%tile,
     &     xA2O_root%xgridroot%xgrid_area,remap_o2a)

      deallocate(ones)
#else
      call Init_Hntrp_Type(remap_a2o,
     &     atmocn%GRID, 0.d0,atmocn%DLATM,
     &      GRID, 0.d0,oDLATM,
     &     0.d0)
      call Init_Hntrp_Type(remap_o2a,
     &      GRID, 0.d0,oDLATM,
     &     atmocn%GRID, 0.d0,atmocn%DLATM,
     &     0.d0)
#endif

C**** Initialize ocean diagnostics metadata
      call init_ODIAG(atmocn)

C**** Set atmospheric surface variables
      IF (.not.postProc) CALL TOC2SST(atmocn)

C***  Interpolate ocean surface velocity to the DYNSI grid
! is this still necessary?
      IF (.not.postProc) CALL OG2IG_uvsurf(dynsice,atmocn)


      RETURN
C**** Terminate because of improper start up
  820 call stop_model('init_OCEAN: Error reading ocean IC',255)
C****
  931 FORMAT ('0Inconsistency between LMM and M:',3I4,2F10.1)

      RETURN
C****
      END SUBROUTINE init_OCEAN

      subroutine tempsalt_oic(fid,mo,g0,gz,s0,sz)
C**** Create MO/G/S on model levels from input file T/S on model
C**** horizontal grid but arbitrary depths.  Input T/S are assumed
C**** to be at specific depths rather than means over depth intervals.
      use ocean, only : im,lmo,lmo_min
      use ocean, only : focean,lmom=>lmm,hocean
      use ocean, only : dzo,zoe=>ze
      use constant, only : grav,rhow
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      use pario, only : read_dist_data,read_data,
     &     variable_exists,get_dimlens,read_attr
      !use ofluxes, only : oapress
      implicit none
      integer :: fid
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     MO,  !  layer mass (kg/m2)
     &     G0,  !  mean potential specific enthalpy (J/kg)
     &     GZ,  !  vertical gradient of potential enthalpy (J/kg)
     &     S0,  !  mean salinity (psu)
     &     SZ   !  vertical gradient of salinity (psu)
c
      real*8,   parameter ::
     &     RHOO=1035,  ! precise value does not matter
     &     zRT3=1/3**.5d0, zRT12=1/12**.5d0
C****
      integer :: kmoic  ! number of OIC levels
      real*8 :: missing

      real*8, dimension(:,:,:), allocatable ::
     &    TOIC,  !  temperature (C)
     &    SOIC   !  salinity (psu)

C**** Variables for single vertical column
      real*8, dimension(:), allocatable ::
     &     ZOIC,!  input depths (m)
     &     GK,  !  potential specific enthalpy (J/kg)
     &     SK,  !  salinity (psu)
     &     PK,  !  pressure (Pa)
     &     VK   !  specific volume (m^3/kg)

C**** Ocean functions of Pressure (Pa), Temperature (C), and Salinity
      real*8 ::
     *       PHPTS,        !  potential specific enthalpy (J/kg)
     *      VOLPTS,VOLGSP  !  specific volume (m^3/kg)
C****
      integer I,J,K,L,LM, ITER, KMIJ, N
      integer :: j_0,j_1,j_0h,j_1h
      logical :: have_north_pole
      real*8 :: PE,ZE, PUP,PDN, VUP,VDN, dMCOL
      real*8 :: ZLO  !  liquid ocean height (m)
      real*8 :: PSL !  sea level pressure (Pa)

      integer :: ndims,dlens(7),ndgood
      character(len=8) :: vnames(3)

      !call stop_model('is oapress set?',255)

      call getdomainbounds(grid,
     &     j_strt = j_0, j_stop = j_1,
     &     j_strt_halo = j_0h, j_stop_halo = j_1h,
     &     have_north_pole = have_north_pole)

c
c basic sanity checks, reading of metadata
c
      vnames(1) = 'depth'
      vnames(2) = 'temp'
      vnames(3) = 'salt'

      do n=1,3
        if(.not.variable_exists(grid,fid,trim(vnames(n)))) then
          call stop_model(
     &         'tempsalt_oic currently expects OIC variable named '//
     &         trim(vnames(n)),255)
        endif
        call get_dimlens(grid,fid,trim(vnames(n)),ndims,dlens)
        if(trim(vnames(n)).eq.'depth') then
          ndgood = 1
        else
          ndgood = 3
        endif
        if(ndims.ne.ndgood) then
          call stop_model(
     &         'tempsalt_oic: wrong rank for OIC variable '//
     &         trim(vnames(n)),255)
        endif
        if(trim(vnames(n)).eq.'depth') kmoic = dlens(1)
      enddo
      missing = 12345678.
      call read_attr(grid,fid,'temp','missing',i,missing)
      call read_attr(grid,fid,'temp','_FillValue',i,missing)
      if(missing == 12345678.) then
        call stop_model(
     &       'tempsalt_oic: no missing value set for temp',255)
      endif

c
c read input variables
c
      allocate(gk(kmoic),sk(kmoic),pk(kmoic),vk(kmoic))
      allocate(zoic(kmoic))
      call read_data(grid,fid,'depth',zoic,bcast_all=.true.)
      do k=2,kmoic
        if(zoic(k).lt.0. .or. zoic(k).lt.zoic(k-1)) then
          call stop_model(
     &     'tempsalt_oic: zoic not monotonically increasing downward',
     &         255)
        endif
      enddo

      allocate(toic(im,j_0h:j_1h,kmoic),soic(im,j_0h:j_1h,kmoic))
      toic = 0.
      soic = 0.
      call read_dist_data(grid,fid,'temp',toic)
      call read_dist_data(grid,fid,'salt',soic)
      where(soic.ne.missing) soic = soic*1d-3 ! psu -> kg/kg

      g0 = missing  ;  gz = missing
      s0 = missing  ;  sz = missing
      mo = 0.

      do j=j_0,j_1
      do i=1,im
        if (focean(i,j) == 0.) cycle

        LM = LMOM(I,J)

C****
C**** Convert ocean temperature to potential specific enthalpy
C****
        GK(:) = MISSING
        K = 1
        PK(K) = 0
        VK(K) = VOLPTS (PK(K), TOIC(I,J,K), SOIC(I,J,K))
        GK(K) =  PHPTS (PK(K), TOIC(I,J,K), SOIC(I,J,K))
        do k=2,kmoic
          If (TOIC(I,J,K) <= MISSING .or. SOIC(I,J,K) <= MISSING) exit
          kmij = k
          VK(K) = VK(K-1)
          do iter=1,3           ! use 3 iterations of PK, VK to converge
            PK(K) = PK(K-1) + GRAV*(ZOIC(K)-ZOIC(K-1))*2/(VK(K-1)+VK(K))
            VK(K) = VOLPTS (PK(K), TOIC(I,J,K), SOIC(I,J,K))
          enddo
          GK(K) =  PHPTS (PK(K), TOIC(I,J,K), SOIC(I,J,K))
        enddo

C****
C**** Remap G and S to model layers; calculate vertical gradients
C****
        SK(1:KMIJ) = SOIC(I,J,1:KMIJ)
        call VLKtoLZ (KMIJ,LMOM(I,J), ZOIC,ZOE(0:LM), GK(1:KMIJ),
     &       G0(I,J,1:LM),GZ(I,J,1:LM), missing, .true.)
        call VLKtoLZ (KMIJ,LMOM(I,J), ZOIC,ZOE(0:LM), SK(1:KMIJ),
     &       S0(I,J,1:LM),SZ(I,J,1:LM), missing, .true.)
C****
C**** Iteratively solve for MO so that integrated Z matches ZLO
C****

C**** Calculate atmospheric surface pressure (Pa)
        PSL = 101325 ! constant as it was in offline version
        !ZLO = - MSI(I,J)/RHOW ! account for the weight of sea ice
        ZLO = 0.!- OAPRESS(I,J)/RHOW/GRAV ! account for the weight of sea ice and atm
                                          ! why rhow and not rhoo
        MO(I,J,1:LM) = RHOO*dZO(1:LM) !  initial guess
C**** Add heights from each layer from ZSOLID (= - HOCEAN)
        do iter=1,10
          PE = PSL - 101325
          ZE = - HOCEAN(I,J)
          DO L=1,LM
            PUP = PE + MO(I,J,L)*GRAV*(.5-zRT12)
            PDN = PE + MO(I,J,L)*GRAV*(.5+zRT12)
            VUP = VOLGSP (G0(I,J,L)-zRT3*GZ(I,J,L),
     *           (S0(I,J,L)-zRT3*SZ(I,J,L)), PUP)
            VDN = VOLGSP (G0(I,J,L)+zRT3*GZ(I,J,L),
     *           (S0(I,J,L)+zRT3*SZ(I,J,L)), PDN)
            PE  = PE + MO(I,J,L)*GRAV
            ZE  = ZE + MO(I,J,L)*(VUP+VDN)*.5
          enddo
          dMCOL = RHOO*(ZE-ZLO) !  excess mass in column
          !MO(I,J,1:LM) = MO(I,J,1:LM) - dMCOL*MFO(1:LM,LM)
          MO(I,J,1:LM) = MO(I,J,1:LM) - dMCOL*(dZO(1:LM) / ZOE(LM))
          if(abs(ZE-ZLO) < 1d-6) exit
        enddo
      enddo
      enddo

      deallocate(toic,soic)

      end subroutine tempsalt_oic


      subroutine VLKtoLZ (KM,LM, MK,ME, RK, RL,RZ, missing,
     &     fill_downward)
C****
C**** VLKtoLZ assumes a continuous piecewise linear tracer distribution,
C**** defined by input tracer concentrations RK at KM specific points.
C**** MK in the downward vertical mass coordinate.
C**** R(M) = {RK(K-1)*[MK(K)-M] + RK(K)*[M-MK(K-1)]} / [MK(K)-MK(K-1)]
C****               when MK(K-1) < M < MK(K).
C**** R(M) = RK(1)  when M < MK(1).
C**** R(M) = is undefined when MK(KM) < M.
C****
C**** VLKtoLZ integrates this tracer distribution over the LM output
C**** layers defined by their layer edges ME, calculating the tracer
C**** mass RM of each layer and the vertical gradient RZ.
C**** RNEW(M) = RL(L) + RZ(L)*[M-MC(L)]/dM(L) when ME(L-1) < M < ME(L)
C**** where MC(L) = .5*[ME(L-1)+ME(L)] and dM(L) = ME(L)-ME(L-1)
C**** Mean concentration of output layers is RL(L) = RM(L)/dM(L).
C****
C**** If ME(L-1) < MK(KM) < ME(L), then RL(L) and RZ(L) are calculated
C**** from the input profile up to MK(KM); RL(L+1:LM) and RZ(L+1:LM)
C**** for deeper layers are undefined, set to MISSING.
C****
C**** Input:  KM = number of input edges
C****         LM = number of output cells
C****         MK = mass coordinates of input points (kg/m^2)
C****         ME = mass coordinates of output layer edges (kg/m^2)
C****         RK = tracer concentration at input points
C****
C**** Output: RL = mean tracer concentration of each output layer
C****         RZ = vertical gradient of tracer mass of each output layer
C****
C**** Internal: RM = integrated tracer mass of output layers (kg/m^2)
C****           RQ = integrated tracer mass times mass (kg^2/m^4)
C****
      implicit none
      integer :: km,lm
      Real*8 MK(KM),ME(0:LM), RK(KM), RL(LM),RZ(LM), RM(1024),RQ(1024)
      real*8 :: missing
      logical :: fill_downward
      real*8 :: mc
      integer :: k,l,ll
C     If (LM > 1024)  Stop 'LM exceeds internal dimentions in VLKtoLZ'
C****
      RM(1:LM) = 0
      RQ(1:LM) = 0
      K = 1
      L = 1
      MC = .5*(ME(L)+ME(L-1))
C****
C**** Integrate layers with M < MK(1)
C****
      If (ME(0) < MK(1))  GoTo 20
C**** MK(1) <= ME(0), determine K such that MK(K-1) <= ME(0) < MK(K)
   10 If (K == KM)  GoTo 200  ;  K = K+1
      If (MK(K) <= ME(0))  GoTo 10
      GoTo 130  !  MK(K-1) <= ME(0) < MK(K)
C**** ME(0) < MK(1), determine output cell containing MK(1)
   20 If (MK(1) < ME(L))  GoTo 30
C**** ME(L-1) < ME(L) < MK(1), integrate RM from ME(L-1) to ME(L)
      RM(L) = RK(1)*(ME(L)-ME(L-1))
      RQ(L) = 0
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      GoTo 20
C**** ME(L-1) < MK(1) < ME(L), integrate RM from ME(L-1) to MK(1)
   30 RM(L) = RK(1)*(MK(1)-ME(L-1))
      RQ(L) = RK(1)*(MK(1)-ME(L-1))*(.5*(MK(1)+ME(L-1))-MC)
      If (K == KM)  GoTo 220  ;  K = K+1
C****
C**** Integrate layers with MK(1) < M < MK(KM)
C****
  100 If (ME(L) < MK(K))  GoTo 120
C**** ME(L-1) < MK(K-1) < MK(K) < ME(L), integrate from MK(K-1) to MK(K)
      RM(L) = RM(L) + (RK(K)-RK(K-1))*(MK(K)+MK(K-1))/2 +
     +                RK(K-1)*MK(K)-RK(K)*MK(K-1)
      RQ(L) = RQ(L) +
     +  (RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +  (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+MK(K-1))/2 +
     +  (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C**** ME(L-1) < MK(K-1) < ME(L) < MK(K), integrate from MK(K-1) to ME(L)
  120 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(ME(L)+MK(K-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-MK(K-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*MK(K-1)+MK(K-1)*MK(K-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+MK(K-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-MK(K-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
  130 If (MK(K) < ME(L))  GoTo 160
C**** MK(K-1) < ME(L-1) < ME(L) < MK(K), integrate from ME(L-1) to ME(L)
  140 RM(L) = ((RK(K)-RK(K-1))*(ME(L)+ME(L-1))/2 +
     +         (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (ME(L)-ME(L-1)) /
     /        (MK(K)-MK(K-1))
      RQ(L) =
     +  ((RK(K)-RK(K-1))*(ME(L)*ME(L)+ME(L)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(ME(L)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (ME(L)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (L == LM)  GoTo 300  ;  L = L+1  ;  MC = .5*(ME(L)+ME(L-1))
      If (ME(L) < MK(K))  GoTo 140
C**** MK(K-1) < ME(L-1) < MK(K) < ME(L), integrate from ME(L-1) to MK(K)
  160 RM(L) = RM(L) + ((RK(K)-RK(K-1))*(MK(K)+ME(L-1))/2 +
     +                 (RK(K-1)*MK(K)-RK(K)*MK(K-1))) * (MK(K)-ME(L-1))
     /              / (MK(K)-MK(K-1))
      RQ(L) = RQ(L) +
     +  ((RK(K)-RK(K-1))*(MK(K)*MK(K)+MK(K)*ME(L-1)+ME(L-1)*ME(L-1))/3 +
     +   (RK(K-1)*(MK(K)+MC)-RK(K)*(MK(K-1)+MC))*(MK(K)+ME(L-1))/2 +
     +   (RK(K)*MK(K-1)-RK(K-1)*MK(K))*MC) * (MK(K)-ME(L-1)) /
     /  (MK(K)-MK(K-1))
      If (K == KM)  GoTo 220  ;  K = K+1
      GoTo 100
C****
C**** Calculate RL and RZ from RM and RQ when MK(KM) < ME(LM)
C****
C**** MK(KM) <= ME(0)
  200 RL(:) = MISSING
      RZ(:) = MISSING
      Return
C**** ME(L-1) < MK(KM) < ME(L)
  220 Do 230 LL=1,L-1
      RL(LL) =   RM(LL) / (ME(LL)-ME(LL-1))
  230 RZ(LL) = 6*RQ(LL) / (ME(LL)-ME(LL-1))**2
      RL(L)  =   RM(L)  / (MK(KM)-ME(L-1))
      RZ(L)  = 6*(RQ(L) + .5*(ME(L)-MK(KM))*RM(L)) / (MK(KM)-ME(L-1))**2
C**** Vertical gradient is extrapolated half way to .5*[MK(KM)+ME(L)]
      RZ(L)  = RZ(L) * (.5*(MK(KM)+ME(L))-ME(L-1)) / (MK(KM)-ME(L-1))
      if(l.lt.lm) then
        if(fill_downward) then
        ! set output points beyond deepest input
        ! point using deepest interpolated value
          rl(l+1:lm) = rl(l)
          rz(l+1:lm) = 0.
        else
          RL(L+1:LM) = MISSING
          RZ(L+1:LM) = MISSING
        endif
      endif
      Return
C****
C**** Calculate RL and RZ from RM and RQ when ME(LM) < MK(KM)
C****
  300 Do 310 L=1,LM
      RL(L) =   RM(L) / (ME(L)-ME(L-1))
  310 RZ(L) = 6*RQ(L) / (ME(L)-ME(L-1))**2
      return
      end subroutine vlktolz

      subroutine init_odiff(grid)
      use OCEAN, only: BYDXYV, BYDXYPJM, UYPB, UYPA, FSLIP
      USE DOMAIN_DECOMP_1D, ONLY : dist_grid,GETDomainBounds,AM_I_ROOT,
     &     HALO_UPDATE, NORTH, SOUTH
      USE OCEAN, only: KHP,KHV,TANP,TANV,BYDXV,BYDXP,BYDYV,
     *     BYDYP,UXA,UXB,UXC,UYA,UYB,UYC,VXA,VXB,VXC,VYA,VYB,VYC
      use ocean, only: dxpo, dypo, dxvo, dyvo, cospo, cosvo, dxyvo,
     &     rlat, IM, JM, LMO, dlat, dxypo, lmu, lmv, bydxypo
      use constant, only: twopi, omega, radius, rhows
      use dictionary_mod
      implicit none
      type (dist_grid) :: grid

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      integer :: I, IM1, J, L
      integer :: ip1
      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE
      REAL*8 DSV,DSP,VLAT,DT2,DTU,DTV,VX,VY,VT,UT,UX,UY

      REAL*8, DIMENSION(grid%j_strt_halo:grid%j_stop_halo) ::
     *     KYPXP,KXPYV,KYVXV,KXVYP
      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo,2) ::
     *      DUDX,DUDY,DVDX,DVDY

!@var AKHMIN minimum horizontal viscosity (m2/s)
!@dbparam AKHFAC tuning factor for horz viscosity
      real*8 :: AKHMIN,AKHFAC=1d0

      if(is_set_param('AKHFAC')) then
        call get_param('AKHFAC',AKHFAC)
      endif

! resolution-dependent min. viscosity settings
      if(jm==46) then
        akhmin = 1.5d8
      elseif(jm==90) then
        akhmin = 5.d6
      elseif(jm==180) then
        akhmin = 1.d5
      else
        call stop_model('init_odiff: set AKHMIN for your res.',255)
      endif

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)


      DO J=J_0,J_1S
C**** Calculate KH = rho_0 BETA* L_Munk^3 where DX=L_Munk
c      KHP(J)=2d0*RHOWS*OMEGA*COSP(J)*(DXP(J)**3)/RADIUS ! tracer lat
c      KHV(J)=2d0*RHOWS*OMEGA*COSV(J)*(DXV(J)**3)/RADIUS ! (v vel pts)
C**** Calculate KH=rho_0 BETA (sqrt(3) L_Munk/pi)^3, L_Munk=min(DX,DY)
        DSP=MIN(DXPO(J),DYPO(J))*2.*SQRT(3.)/TWOPI  ! tracer lat
        DSV=MIN(DXVO(J),DYVO(J))*2.*SQRT(3.)/TWOPI  ! v vel pts
        KHP(J)=AKHFAC*2d0*RHOWS*OMEGA*COSPO(J)*(DSP**3)/RADIUS ! A-grid
        KHV(J)=AKHFAC*2d0*RHOWS*OMEGA*COSVO(J)*(DSV**3)/RADIUS ! V-grid
        KHP(J)=MAX(KHP(J),AKHFAC*AKHMIN)
        KHV(J)=MAX(KHV(J),AKHFAC*AKHMIN)
        BYDXYV(J)=1D0/DXYVO(J)
        BYDXV(J)=1D0/DXVO(J)
        BYDXP(J)=1D0/DXPO(J)
        BYDYV(J)=1D0/DYVO(J)
        BYDYP(J)=1D0/DYPO(J)
        KYPXP(J)=KHP(J)*DYPO(J)*BYDXP(J)
#ifdef ODIFF_FIXES_2017 /* for hemispheric symmetry */
        KXPYV(J)=KHV(J)*DXVO(J)*BYDYV(J)
        KXVYP(J)=KHP(J)*DXPO(J)*BYDYP(J)
#else
        KXPYV(J)=KHV(J)*DXPO(J)*BYDYV(J)
        KXVYP(J)=KHP(J)*DXVO(J)*BYDYP(J)
#endif
        KYVXV(J)=KHV(J)*DYVO(J)*BYDXV(J)
C**** Discretisation errors need TANP/V to be defined like this
        TANP(J)=TAN(RLAT(J))*TAN(0.5*DLAT)/(RADIUS*0.5*DLAT)
        VLAT = DLAT*(J+0.5-0.5*(1+JM))
        TANV(J)=TAN(VLAT)*SIN(DLAT)/(DLAT*RADIUS)
c       write(*,'(a,2i5,3e12.4)')'for samar, dx,dy:',
c    .      nstep,j,dxpo(j),dypo(j),dxypo(j)
      END DO
      !make halo_update if 2 is not present
      !barrier synchoronization is required before sending the message,
      !j=2 is computed before the message is sent
      if( HAVE_SOUTH_POLE ) then
        KHV(1)=KHV(2)
        BYDXV(1)=1D0/DXVO(1)
        BYDYV(1)=1D0/DYVO(1)
        BYDYP(1)=1D0/DYPO(1)
        TANP(1) = 0
c       write(*,'(a,2i5,3e12.4)')'for samar, dx,dy:',
c    .       nstep,1,dxpo(1),dypo(1),dxypo(1)
      endif
      if( HAVE_NORTH_POLE ) then
        BYDXP(JM)=1D0/DXPO(JM)
        BYDXYPJM=1D0/(DXYPO(JM)*IM)
        TANP(JM) = 0
c       write(*,'(a,2i5,3e12.4)')'for samar, dx,dy:',
c    .       nstep,jm,dxpo(jm),dypo(jm),dxypo(jm)
      endif

      CALL HALO_UPDATE(grid,TANP (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(grid,BYDXP(grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,TANV (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,KHP (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,KXVYP (grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH+NORTH)
!      DYPO has no halo !!
!      CALL HALO_UPDATE(grid,DYPO (grid%j_strt_halo:grid%j_stop_halo) ,
!     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,KXPYV,FROM=SOUTH)
      CALL HALO_UPDATE(grid,KYVXV,FROM=SOUTH)
      CALL HALO_UPDATE(grid,KYPXP,FROM=SOUTH)

C****
C**** Calculate operators fixed in time for U and V equations
C****
      UXA=0. ; UXB=0. ; UXC=0. ; UYA=0. ; UYB=0. ; UYC=0.
      VXA=0. ; VXB=0. ; VXC=0. ; VYA=0. ; VYB=0. ; VYC=0.
      UYPB=0.; UYPA=0.

      DO L=1,LMO
C**** Calculate flux operators
C**** i.e DUDX(1) and DUDX(2) are the coefficients of u1 and u2 for
C**** calculating centered difference K_h du/dx
C**** including metric terms in y derivatives
        DUDX=0.
        DUDY=0.
        DVDX=0.
        DVDY=0.
        DO J=max(2,J_0S-1),J_1S
          I=IM
          DO IP1=1,IM
            IF (L.LE.LMU(IP1,J)) DUDX(IP1,J,1) = KYPXP(J)
            IF (L.LE.LMU(I,J)) THEN
              DUDX(IP1,J,2) = -KYPXP(J)
              IF (L.LE.LMU(I,J+1)) THEN
                DUDY(I,J,1) =  KXPYV(J)*(1. +0.5*TANV(J)*DYVO(J))
                DUDY(I,J,2) = -KXPYV(J)*(1. -0.5*TANV(J)*DYVO(J))
              ELSE
                DUDY(I,J,2) = -(1.-FSLIP)*2d0*KXPYV(J)
              END IF
            ELSE
              IF (L.LE.LMU(I,J+1)) DUDY(I,J,1) = (1.-FSLIP)*2d0*KXPYV(J)
            END IF
            IF (L.LE.LMV(I,J+1)) THEN
#ifdef ODIFF_FIXES_2017 /* for hemispheric symmetry */
            if(j+1.ne.jm) then
              DVDY(I,J+1,1) = KXVYP(J+1)*(1. + 0.5*TANP(J+1)*DYPO(J+1))
            else ! previous line blows up at the NP - decision pending
              DVDY(I,J+1,1) = KXVYP(J)*(1. + 0.5*TANP(J)*DYPO(J))
            endif
#else
            DVDY(I,J+1,1) = KXVYP(J)*(1. + 0.5*TANP(J)*DYPO(J))
#endif
            ENDIF
            IF (L.LE.LMV(I,J)) THEN
#ifdef ODIFF_FIXES_2017 /* for hemispheric symmetry */
              if(j+1.ne.jm) then
              DVDY(I,J+1,2) = -KXVYP(J+1)*(1.-0.5*TANP(J+1)*DYPO(J+1))
              else ! previous line blows up at the NP - decision pending
              DVDY(I,J+1,2) = -KXVYP(J)*(1.-0.5*TANP(J)*DYPO(J))
              endif
#else
              DVDY(I,J+1,2) = -KXVYP(J)*(1.-0.5*TANP(J)*DYPO(J))
#endif
              IF (L.LE.LMV(IP1,J)) THEN
                DVDX(I,J,1) =  KYVXV(J)
                DVDX(I,J,2) = -KYVXV(J)
              ELSE
                DVDX(I,J,2) = -(1.-FSLIP)*2d0*KYVXV(J)
              END IF
            ELSE
              IF (L.LE.LMV(IP1,J)) DVDX(I,J,1) = (1.-FSLIP)*2d0*KYVXV(J)
            END IF
            I=IP1
          END DO
        END DO
#ifdef ODIFF_FIXES_2017 /* for SP as a 1-gridpoint island */
        if( HAVE_SOUTH_POLE ) then
          j = 1
          do i=1,im
            !IF (L.LE.LMU(I,J)) THEN
            !ELSE
            IF(L.LE.LMU(I,J+1)) DUDY(I,J,1) = (1.-FSLIP)*2d0*KXPYV(1)
            !ENDIF
          enddo
        endif
#endif
        CALL HALO_UPDATE(grid, DUDY(:,grid%j_strt_halo:grid%j_stop_halo,
     *     :), FROM=SOUTH)
C****
C**** Combine to form tri-diagonal operators including first metric term
C****
        DO J=J_0S,J_1S
          IM1=IM
          DO I=1,IM
            IF(L.LE.LMU(IM1,J)) THEN
              UXA(IM1,J,L) = -DUDX(IM1,J,2)                  *BYDXYPO(J)
              UXB(IM1,J,L) = (DUDX(I  ,J,2) -DUDX(IM1,J  ,1))*BYDXYPO(J)
              UXC(IM1,J,L) =  DUDX(I  ,J,1)                  *BYDXYPO(J)
              UYA(IM1,J,L) = -DUDY(IM1,J-1,2)                *BYDXYPO(J)
     *                     + 0.5*TANP(J)*KHP(J)*BYDYP(J)
              UYB(IM1,J,L) =(DUDY(IM1,J  ,2)-DUDY(IM1,J-1,1))*BYDXYPO(J)
#ifdef ODIFF_FIXES_2017
     *                     - TANP(J)*TANP(J)*KHP(J)
#else
     *                     + TANP(J)*TANP(J)*KHP(J)
#endif
              UYC(IM1,J,L) = DUDY(IM1,J  ,1)                 *BYDXYPO(J)
     *                     - 0.5*TANP(J)*KHP(J)*BYDYP(J)
            END IF
            IF (L.LE.LMV(I,J)) THEN
              VXA(I,J,L) = -DVDX(IM1,J,2)                 *BYDXYV(J)
              VXB(I,J,L) = (DVDX(I  ,J,2) - DVDX(IM1,J,1))*BYDXYV(J)
              VXC(I,J,L) =  DVDX(I  ,J,1)                 *BYDXYV(J)
              VYA(I,J,L) = -DVDY(I,J  ,2)                 *BYDXYV(J)
     *                   + 0.5*TANV(J)*KHV(J)*BYDYV(J)
              VYB(I,J,L) = (DVDY(I,J+1,2) - DVDY(I  ,J,1))*BYDXYV(J)
#ifdef ODIFF_FIXES_2017
     *                   - TANV(J)*TANV(J)*KHV(J)
#else
     *                   + TANV(J)*TANV(J)*KHV(J)
#endif
              VYC(I,J,L) =  DVDY(I,J+1,1)                 *BYDXYV(J)
     *                   - 0.5*TANV(J)*KHV(J)*BYDYV(J)
            END IF
            IM1=I
          END DO
        END DO
C**** At North Pole
        !following uses jm and jm-1 data, if both are not on same
        !processor (need to check this), halo_update is required.
        IF(HAVE_NORTH_POLE) THEN
          IF(L.LE.LMU(1,JM)) THEN
            DO I=1,IM
              UYPB(L) = UYPB(L) - DUDY(I,JM-1,1)
              UYPA(I,L) = -DUDY(I,JM-1,2)*BYDXYPJM
            END DO
            UYPB(L) = UYPB(L)*BYDXYPJM
          END IF
        END IF
      END DO
!**  IFIRST block ends here
      end subroutine init_odiff

      subroutine get_i1i2(q,im,n,i1,i2,nmax)
c Finds the number of intervals over which q==true, and the
c beginning/ending indices of each interval.
c Wraparound is disabled.
      implicit none
      integer, intent(in) :: im  ! number of points
      logical, dimension(im), intent(in) :: q ! logical mask
      integer, intent(in) :: nmax ! max number of contiguous intervals
      integer, intent(out) :: n ! number of contiguous intervals
      integer, dimension(nmax) :: i1,i2 ! start,end of each interval
      integer :: i,nn
c first find the number of intervals
      n = 0
      do i=1,im-1
        if(q(i) .and. .not.q(i+1)) n = n + 1
      enddo
      if(q(im)) n = n + 1
      if(n.gt.nmax) then
        write(6,*) 'for get_i1i2, increase nmax to ',n
        call stop_model('get_i1i2: n>nmax',255)
      endif
c find the start/end of each interval
      i = 1
      do nn=1,n
        do while(.not.q(i))
          i = i + 1
        enddo
        i1(nn) = i
        do while(q(i))
          i = i + 1
          if(i.gt.im) exit
        enddo
        i2(nn) = i-1
      enddo
      return
      end subroutine get_i1i2

      SUBROUTINE daily_OCEAN(end_of_day,atmocn)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day
      type(atmocn_xchng_vars) :: atmocn

C**** Only do this at end of the day
      IF (end_of_day) THEN

C**** Add glacial melt from Antarctica and Greenland
        CALL GLMELT(SECONDS_PER_DAY)
        if(associated(atmocn%consrv)) then ! not (yet) true for ocean-only
          CALL DIAGCO (10,ATMOCN)
        endif

c uncomment following call to activate tracers at arbitrary times
c#ifdef TRACERS_OCEAN
c      call tracer_ic_ocean(atmocn)
c#endif

C**** set gtemp arrays for ocean
        CALL TOC2SST(atmocn)
      END IF
C****
      RETURN
      END SUBROUTINE daily_OCEAN

      subroutine def_rsf_ocean(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      use ocean
      use straits
      USE OCEANR_DIM, only : grid=>ogrid
#ifdef TRACERS_OCEAN
      Use OCN_TRACER_COM, Only : tracerlist, ocn_tracer_entry
#endif
      use pario, only : defvar
      use domain_decomp_1d, only : getDomainBounds
      implicit none
      integer fid   !@var fid file id
      integer :: n
      integer :: i_0h,i_1h, j_0h,j_1h
      real*8, dimension(:,:), allocatable :: arrdum
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif

      call defvar(grid,fid,mo,'mo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,uo,'uo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,vo,'vo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,uod,'uod(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,vod,'vod(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,g0m,'g0m(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gxmo,'gxmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gymo,'gymo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gzmo,'gzmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,s0m,'s0m(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,sxmo,'sxmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,symo,'symo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,szmo,'szmo(dist_imo,dist_jmo,lmo)')
      if(use_qus==1) then
      call defvar(grid,fid,gxxmo,'gxxmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gyymo,'gyymo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gzzmo,'gzzmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gxymo,'gxymo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gyzmo,'gyzmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gzxmo,'gzxmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,sxxmo,'sxxmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,syymo,'syymo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,szzmo,'szzmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,sxymo,'sxymo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,syzmo,'syzmo(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,szxmo,'szxmo(dist_imo,dist_jmo,lmo)')
      endif
      call defvar(grid,fid,ogeoz,'ogeoz(dist_imo,dist_jmo)')
      call defvar(grid,fid,ogeoz_sv,'ogeoz_sv(dist_imo,dist_jmo)')
      call defvar(grid,fid,kpl,'kpl(dist_imo,dist_jmo)')
      if(nmst.gt.0) then
c straits arrays
      call defvar(grid,fid,must,'must(lmo,nmst)')
      call defvar(grid,fid,g0mst,'g0mst(lmo,nmst)')
      call defvar(grid,fid,gxmst,'gxmst(lmo,nmst)')
      call defvar(grid,fid,gzmst,'gzmst(lmo,nmst)')
      call defvar(grid,fid,s0mst,'s0mst(lmo,nmst)')
      call defvar(grid,fid,sxmst,'sxmst(lmo,nmst)')
      call defvar(grid,fid,szmst,'szmst(lmo,nmst)')
!      call defvar(grid,fid,rsist,'rsist(nmst)')
!      call defvar(grid,fid,rsixst,'rsixst(nmst)')
!      call defvar(grid,fid,msist,'msist(two,nmst)')
!      call defvar(grid,fid,hsist,'hsist(lmi,nmst)')
!      call defvar(grid,fid,ssist,'ssist(lmi,nmst)')
      endif
#ifdef TRACERS_OCEAN
c tracer arrays
      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
        call defvar(grid,fid,trmo(:,:,:,n),
     &       'trmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,txmo(:,:,:,n),
     &       'txmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,tymo(:,:,:,n),
     &       'tymo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,tzmo(:,:,:,n),
     &       'tzmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        if(use_qus==1) then
        call defvar(grid,fid,txxmo(:,:,:,n),
     &       'txxmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,tyymo(:,:,:,n),
     &       'tyymo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,tzzmo(:,:,:,n),
     &       'tzzmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,txymo(:,:,:,n),
     &       'txymo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,tyzmo(:,:,:,n),
     &       'tyzmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        call defvar(grid,fid,tzxmo(:,:,:,n),
     &       'tzxmo_'//trim(entry%trname)//'(dist_imo,dist_jmo,lmo)')
        endif
      enddo
      if(nmst.gt.0) then
c tracer arrays in straits
      call defvar(grid,fid,trmst,'trmst(lmo,nmst,ntmo)')
      call defvar(grid,fid,txmst,'txmst(lmo,nmst,ntmo)')
      call defvar(grid,fid,tzmst,'tzmst(lmo,nmst,ntmo)')
!#ifdef TRACERS_WATER
!      call defvar(grid,fid,trsist,'trsist(ntmo,lmi,nmst)')
!#endif
      endif
#ifdef TRACERS_OceanBiology
      call def_rsf_obio(fid)
#endif
#endif

      call getDomainBounds(grid, i_strt_halo=i_0h,i_stop_halo=i_1h,
     &               j_strt_halo=j_0h,j_stop_halo=j_1h)
      allocate(arrdum(i_0h:i_1h,j_0h:j_1h))
      call defvar(grid, fid, arrdum, 'wliqo(dist_imo,dist_jmo)')
      call defvar(grid, fid, arrdum, 'ekliqo(dist_imo,dist_jmo)')
      call defvar(grid, fid, arrdum, 'epliqo(dist_imo,dist_jmo)')
      call defvar(grid, fid, arrdum, 'sliqo(dist_imo,dist_jmo)')
      deallocate(arrdum)

#ifdef OCN_GISS_TURB
      call def_rsf_gissmix(fid)
#endif
#ifdef OCN_GISS_SM
      call def_rsf_giss_sm(fid)
#endif

      return
      end subroutine def_rsf_ocean

      subroutine new_io_ocean(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use ocean
      use straits
      USE OCEANR_DIM, only : grid=>ogrid
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
#ifdef TRACERS_OCEAN
      Use OCN_TRACER_COM, Only : tracerlist, ocn_tracer_entry
#endif
      use domain_decomp_1d, only : getDomainBounds
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n
      integer :: i_0h,i_1h, j_0h,j_1h
      real*8, dimension(:,:), allocatable :: arrdum
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif

      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'mo',mo)
        call write_dist_data(grid,fid,'uo',uo)
        call write_dist_data(grid,fid,'vo',vo)
        call write_dist_data(grid,fid,'uod',uod)
        call write_dist_data(grid,fid,'vod',vod)
        call write_dist_data(grid,fid,'g0m',g0m)
        call write_dist_data(grid,fid,'gxmo',gxmo)
        call write_dist_data(grid,fid,'gymo',gymo)
        call write_dist_data(grid,fid,'gzmo',gzmo)
        call write_dist_data(grid,fid,'s0m',s0m)
        call write_dist_data(grid,fid,'sxmo',sxmo)
        call write_dist_data(grid,fid,'symo',symo)
        call write_dist_data(grid,fid,'szmo',szmo)
        if(use_qus==1) then
        call write_dist_data(grid,fid,'gxxmo',gxxmo)
        call write_dist_data(grid,fid,'gyymo',gyymo)
        call write_dist_data(grid,fid,'gzzmo',gzzmo)
        call write_dist_data(grid,fid,'gxymo',gxymo)
        call write_dist_data(grid,fid,'gyzmo',gyzmo)
        call write_dist_data(grid,fid,'gzxmo',gzxmo)
        call write_dist_data(grid,fid,'sxxmo',sxxmo)
        call write_dist_data(grid,fid,'syymo',syymo)
        call write_dist_data(grid,fid,'szzmo',szzmo)
        call write_dist_data(grid,fid,'sxymo',sxymo)
        call write_dist_data(grid,fid,'syzmo',syzmo)
        call write_dist_data(grid,fid,'szxmo',szxmo)
        endif
        call write_dist_data(grid,fid,'ogeoz',ogeoz)
        call write_dist_data(grid,fid,'ogeoz_sv',ogeoz_sv)
        call write_dist_data(grid,fid,'kpl',kpl)
        if(nmst.gt.0) then
c straits arrays
        call write_data(grid,fid,'must',must)
        call write_data(grid,fid,'g0mst',g0mst)
        call write_data(grid,fid,'gxmst',gxmst)
        call write_data(grid,fid,'gzmst',gzmst)
        call write_data(grid,fid,'s0mst',s0mst)
        call write_data(grid,fid,'sxmst',sxmst)
        call write_data(grid,fid,'szmst',szmst)
!        call write_data(grid,fid,'rsist',rsist)
!        call write_data(grid,fid,'rsixst',rsixst)
!        call write_data(grid,fid,'msist',msist)
!        call write_data(grid,fid,'hsist',hsist)
!        call write_data(grid,fid,'ssist',ssist)
        endif
#ifdef TRACERS_OCEAN
c tracer arrays
        do n=1,tracerlist%getsize()
          entry=>tracerlist%at(n)
          call write_dist_data(grid,fid,'trmo_'//trim(entry%trname),
     &         trmo(:,:,:,n))
          call write_dist_data(grid,fid,'txmo_'//trim(entry%trname),
     &         txmo(:,:,:,n))
          call write_dist_data(grid,fid,'tymo_'//trim(entry%trname),
     &         tymo(:,:,:,n))
          call write_dist_data(grid,fid,'tzmo_'//trim(entry%trname),
     &         tzmo(:,:,:,n))
          if(use_qus==1) then
          call write_dist_data(grid,fid,'txxmo_'//trim(entry%trname),
     &           txxmo(:,:,:,n))
          call write_dist_data(grid,fid,'tyymo_'//trim(entry%trname),
     &           tyymo(:,:,:,n))
          call write_dist_data(grid,fid,'tzzmo_'//trim(entry%trname),
     &           tzzmo(:,:,:,n))
          call write_dist_data(grid,fid,'txymo_'//trim(entry%trname),
     &           txymo(:,:,:,n))
          call write_dist_data(grid,fid,'tyzmo_'//trim(entry%trname),
     &           tyzmo(:,:,:,n))
          call write_dist_data(grid,fid,'tzxmo_'//trim(entry%trname),
     &           tzxmo(:,:,:,n))
          endif
        enddo
        if(nmst.gt.0) then
c tracer arrays in straits
        call write_data(grid,fid,'trmst',trmst)
        call write_data(grid,fid,'txmst',txmst)
        call write_data(grid,fid,'tzmst',tzmst)
!#ifdef TRACERS_WATER
!        call write_data(grid,fid,'trsist',trsist)
!#endif
        endif
#endif
        call getDomainBounds(grid, i_strt_halo=i_0h,i_stop_halo=i_1h,
     &                 j_strt_halo=j_0h,j_stop_halo=j_1h)
        allocate(arrdum(i_0h:i_1h,j_0h:j_1h))
        if(grid%have_domain) call conserv_OMS(arrdum)
        call write_dist_data(grid,fid,'wliqo',arrdum)
        if(grid%have_domain) call conserv_OKE(arrdum)
        call write_dist_data(grid,fid,'ekliqo',arrdum)
        if(grid%have_domain) call conserv_OCE(arrdum)
        call write_dist_data(grid,fid,'epliqo',arrdum)
        if(grid%have_domain) call conserv_OSL(arrdum)
        call write_dist_data(grid,fid,'sliqo',arrdum)
        deallocate(arrdum)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'mo',mo)
        call read_dist_data(grid,fid,'uo',uo)
        call read_dist_data(grid,fid,'vo',vo)
        call read_dist_data(grid,fid,'uod',uod)
        call read_dist_data(grid,fid,'vod',vod)
        call read_dist_data(grid,fid,'g0m',g0m)
        call read_dist_data(grid,fid,'gxmo',gxmo)
        call read_dist_data(grid,fid,'gymo',gymo)
        call read_dist_data(grid,fid,'gzmo',gzmo)
        call read_dist_data(grid,fid,'s0m',s0m)
        call read_dist_data(grid,fid,'sxmo',sxmo)
        call read_dist_data(grid,fid,'symo',symo)
        call read_dist_data(grid,fid,'szmo',szmo)
        if(use_qus==1) then
        call read_dist_data(grid,fid,'gxxmo',gxxmo)
        call read_dist_data(grid,fid,'gyymo',gyymo)
        call read_dist_data(grid,fid,'gzzmo',gzzmo)
        call read_dist_data(grid,fid,'gxymo',gxymo)
        call read_dist_data(grid,fid,'gyzmo',gyzmo)
        call read_dist_data(grid,fid,'gzxmo',gzxmo)
        call read_dist_data(grid,fid,'sxxmo',sxxmo)
        call read_dist_data(grid,fid,'syymo',syymo)
        call read_dist_data(grid,fid,'szzmo',szzmo)
        call read_dist_data(grid,fid,'sxymo',sxymo)
        call read_dist_data(grid,fid,'syzmo',syzmo)
        call read_dist_data(grid,fid,'szxmo',szxmo)
        endif
        call read_dist_data(grid,fid,'ogeoz',ogeoz)
        call read_dist_data(grid,fid,'ogeoz_sv',ogeoz_sv)
        call read_dist_data(grid,fid,'kpl',kpl)
        if(nmst.gt.0) then
c straits arrays
        call read_data(grid,fid,'must',must,bcast_all=.true.)
        call read_data(grid,fid,'g0mst',g0mst,bcast_all=.true.)
        call read_data(grid,fid,'gxmst',gxmst,bcast_all=.true.)
        call read_data(grid,fid,'gzmst',gzmst,bcast_all=.true.)
        call read_data(grid,fid,'s0mst',s0mst,bcast_all=.true.)
        call read_data(grid,fid,'sxmst',sxmst,bcast_all=.true.)
        call read_data(grid,fid,'szmst',szmst,bcast_all=.true.)
!        call read_data(grid,fid,'rsist',rsist,bcast_all=.true.)
!        call read_data(grid,fid,'rsixst',rsixst,bcast_all=.true.)
!        call read_data(grid,fid,'msist',msist,bcast_all=.true.)
!        call read_data(grid,fid,'hsist',hsist,bcast_all=.true.)
!        call read_data(grid,fid,'ssist',ssist,bcast_all=.true.)
        endif
#ifdef TRACERS_OCEAN
c tracer arrays
        do n=1,tracerlist%getsize()
          entry=>tracerlist%at(n)
          call read_dist_data(grid,fid,'trmo_'//trim(entry%trname),
     &         trmo(:,:,:,n))
          call read_dist_data(grid,fid,'txmo_'//trim(entry%trname),
     &         txmo(:,:,:,n))
          call read_dist_data(grid,fid,'tymo_'//trim(entry%trname),
     &         tymo(:,:,:,n))
          call read_dist_data(grid,fid,'tzmo_'//trim(entry%trname),
     &         tzmo(:,:,:,n))
          if(use_qus==1) then
          call read_dist_data(grid,fid,'txxmo_'//trim(entry%trname),
     &         txxmo(:,:,:,n))
          call read_dist_data(grid,fid,'tyymo_'//trim(entry%trname),
     &         tyymo(:,:,:,n))
          call read_dist_data(grid,fid,'tzzmo_'//trim(entry%trname),
     &         tzzmo(:,:,:,n))
          call read_dist_data(grid,fid,'txymo_'//trim(entry%trname),
     &         txymo(:,:,:,n))
          call read_dist_data(grid,fid,'tyzmo_'//trim(entry%trname),
     &         tyzmo(:,:,:,n))
          call read_dist_data(grid,fid,'tzxmo_'//trim(entry%trname),
     &         tzxmo(:,:,:,n))
          endif
        enddo
        if(nmst.gt.0) then
c tracer arrays in straits
        call read_data(grid,fid,'trmst',trmst,bcast_all=.true.)
        call read_data(grid,fid,'txmst',txmst,bcast_all=.true.)
        call read_data(grid,fid,'tzmst',tzmst,bcast_all=.true.)
!#ifdef TRACERS_WATER
!        call read_data(grid,fid,'trsist',trsist,bcast_all=.true.)
!#endif
        endif
#endif
      end select

#ifdef TRACERS_OceanBiology
      call new_io_obio(fid,iaction)
#endif

#ifdef OCN_GISS_TURB
      call new_io_gissmix(fid,iaction)
#endif
#ifdef OCN_GISS_SM
      call new_io_giss_sm(fid,iaction)
#endif

      return
      end subroutine new_io_ocean

      SUBROUTINE CHECKO_serial(SUBR)
!@sum  CHECKO Checks whether Ocean variables are reasonable (serial vers
!@auth Original Development Team
      USE CONSTANT, only : byrt3,teeny
      USE MODEL_COM, only : qcheck
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#endif
      USE OCEAN, only : im,jm,lmo,dxypo,focean,
     *     imaxj, lmm
     *  ,mo, uo, vo
     &  ,g0m,gxmo,gymo,gzmo
     &  ,s0m,sxmo,symo,szmo
c     *  , mo=>mo_glob
c     *  ,g0m=>g0m_glob, gxmo=>gxmo_glob,gymo=>gymo_glob,gzmo=>gzmo_glob
c     *  ,s0m=>s0m_glob, sxmo=>sxmo_glob,symo=>symo_glob,szmo=>szmo_glob
c     *  ,uo=>uo_glob, vo=>vo_glob
#ifdef TRACERS_OCEAN
     *  ,trmo, txmo, tymo, tzmo
c     *  ,trmo=>trmo_glob, txmo=>txmo_glob, tymo=>tymo_glob,
c     *   tzmo=>tzmo_glob
#endif
      IMPLICIT NONE
      REAL*8 SALIM,GO1,SO1,relerr,errmax,temgs
      LOGICAL QCHECKO
      INTEGER I,J,L,n,imax,jmax,lmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      integer :: njpol
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif

      njpol = 1
C**** Check for NaN/INF in ocean data
      IF (QCHECK) THEN
      CALL CHECK3B(MO  ,1,IM,1,JM,NJPOL,LMO,SUBR,'mo ')
      CALL CHECK3B(G0M ,1,IM,1,JM,NJPOL,LMO,SUBR,'g0m')
      CALL CHECK3B(GXMO,1,IM,1,JM,NJPOL,LMO,SUBR,'gxm')
      CALL CHECK3B(GYMO,1,IM,1,JM,NJPOL,LMO,SUBR,'gym')
      CALL CHECK3B(GZMO,1,IM,1,JM,NJPOL,LMO,SUBR,'gzm')
      CALL CHECK3B(S0M ,1,IM,1,JM,NJPOL,LMO,SUBR,'s0m')
      CALL CHECK3B(SXMO,1,IM,1,JM,NJPOL,LMO,SUBR,'sxm')
      CALL CHECK3B(SYMO,1,IM,1,JM,NJPOL,LMO,SUBR,'sym')
      CALL CHECK3B(SZMO,1,IM,1,JM,NJPOL,LMO,SUBR,'szm')
      CALL CHECK3B(UO  ,1,IM,1,JM,NJPOL,LMO,SUBR,'uo ')
      CALL CHECK3B(VO  ,1,IM,1,JM,NJPOL,LMO,SUBR,'vo ')
#ifdef TRACERS_OCEAN
      CALL CHECK4B(TRMO,1,IM,1,JM,NJPOL,LMO,tracerlist%getsize(),
     &                                                SUBR,'tzm')
#endif

C**** Check for variables out of bounds
      QCHECKO=.FALSE.
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).gt.0.) THEN
C**** Check potential specific enthalpy/salinity
          DO L=1,LMM(I,J)
          GO1 = G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          SO1 = S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          IF(GO1.lt.-10000. .or. GO1.gt.200000.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,GO=',I,J,L,GO1,TEMGS(GO1
     *           ,SO1)
            IF (GO1.lt.-20000. .or. GO1.gt.200000.) QCHECKO=.TRUE.
          END IF
          IF(SO1.gt.0.045 .or. SO1.lt.0.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,SO=',I,J,L,1d3*SO1
            IF ((SO1.gt.0.05 .or. SO1.lt.0.) .and. .not. (I == 47.and.
     *           J == 30)) QCHECKO=.TRUE.
          END IF
          END DO
C**** Check all ocean currents
          DO L = 1,LMO
            IF(ABS(UO(I,J,L)).gt.7. .or. ABS(VO(I,J,L)).gt.5) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,UO,VO=',I,J,L,UO(I,J,L)
     *             ,VO(I,J,L)
              QCHECKO=.TRUE.
            END IF
          END DO
C**** Check first layer ocean mass
          IF(MO(I,J,1).lt.2000. .or. MO(I,J,1).gt.20000.) THEN
!!          IF (I == 47.and.(J == 33.or.J == 34)) GOTO 230 ! not Caspian
            WRITE (6,*) 'After ',SUBR,': I,J,MO=',I,J,MO(I,J,1)
!!          QCHECKO=.TRUE.
          END IF
C**** Check ocean salinity in each eighth box for the first layer
 230      SALIM = .045
!!        IF(JM == 46 .and. I == 47 .and. J == 30) GOTO 240 ! not Persia
C                                                         n Gulf   !.048
          IF(.5*ABS(SXMO(I,J,1))+.5*ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J
     *         ,1)).lt.S0M(I,J,1) .and.(.5*ABS(SXMO(I,J,1))+.5
     *         *ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J,1))+S0M(I,J,1))
     *         /(MO(I,J,1)*DXYPO(J)).lt.SALIM)  GO TO 240
          WRITE (6,*) 'After ',SUBR,': I,J,S0,SX,SY,SZ=',I,J,
     *         1d3*S0M(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3*SXMO(I,J,1)/(MO(I
     *         ,J,1)*DXYPO(J)),1d3*SYMO(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3
     *         *SZMO(I,J,1)/(MO(I,J,1)*DXYPO(J))
!!        QCHECKO=.TRUE.
 240      CONTINUE
        END IF
      END DO
      END DO

#ifdef TRACERS_OCEAN
      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
C**** Check for negative tracers
        if (entry%t_qlimit) then
        do l=1,lmo
        do j=1,jm
        do i=1,imaxj(j)
          if (l.le.lmm(i,j)) then
            if (trmo(i,j,l,n).lt.0) then
              write(6,*) "Neg Tracer in ocean after ",subr,i,j,l,
     *             entry%trname,trmo(i,j,l,n)
              QCHECKO=.true.
            end if
          end if
        end do
        end do
        end do
        end if
C**** Check conservation of water tracers in ocean
        if (entry%trname == 'Water') then
          errmax = 0. ; imax=1 ; jmax=1 ; lmax=1
          do l=1,lmo
          do j=1,jm
          do i=1,imaxj(j)
            if (l.le.lmm(i,j)) then
              relerr=max(
     *             abs(trmo(i,j,l,n)-mo(i,j,l)*dxypo(j)+s0m(i,j,l)),
     *             abs(txmo(i,j,l,n)+sxmo(i,j,l)),
     *             abs(tymo(i,j,l,n)+symo(i,j,l)),
     *             abs(tzmo(i,j,l,n)+szmo(i,j,l)))/
     *             (mo(i,j,l)*dxypo(j)-s0m(i,j,l))
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; lmax=l ; errmax=relerr
            end if
            end if
          end do
          end do
          end do
          print*,"Relative error in ocean fresh water mass after ",subr
     *         ,":",imax,jmax,lmax,errmax,trmo(imax,jmax,lmax,n),mo(imax
     *         ,jmax,lmax)*dxypo(jmax)-s0m(imax,jmax,lmax),txmo(imax
     *         ,jmax,lmax,n),-sxmo(imax,jmax,lmax),tymo(imax,jmax,lmax,n
     *         ),-symo(imax,jmax ,lmax),tzmo(imax,jmax,lmax,n),
     *         -szmo(imax,jmax,lmax)
        end if
      end do
#endif

      IF (QCHECKO)
     &     call stop_model("QCHECKO: Ocean Variables out of bounds",255)

      END IF
C****
      CALL CHECKOST(SUBR)
C****
      END SUBROUTINE CHECKO_serial

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether Ocean variables are reasonable (parallel ve
!@auth Original Development Team
      USE CONSTANT, only : byrt3,teeny
      USE MODEL_COM, only : qcheck
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#endif
      USE OCEAN
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds, AM_I_ROOT
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8 SALIM,GO1,SO1,relerr,errmax,temgs
      LOGICAL QCHECKO
      INTEGER I,J,L,n,imax,jmax,lmax
!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

c**** Extract domain decomposition info
      INTEGER :: J_0S, J_0, J_1, J_0H, J_1H, JM_loc, njpol
      INTEGER :: J_0STG,J_1STG
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif

      call getDomainBounds(grid, J_STRT_SKP=J_0S, 
     *   J_STRT=J_0, J_STOP=J_1,
     *   J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)
      J_0STG = grid%J_STRT_STGR
      J_1STG = grid%J_STOP_STGR
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in ocean data
      IF (QCHECK) THEN
      JM_loc = J_1H - J_0H + 1
      CALL CHECK3B(MO(:,J_0:J_1,:)  ,1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'mo ')
      CALL CHECK3B(G0M(:,J_0:J_1,:) ,1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'g0m')
      CALL CHECK3B(GXMO(:,J_0:J_1,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'gxm')
      CALL CHECK3B(GYMO(:,J_0:J_1,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'gym')
      CALL CHECK3B(GZMO(:,J_0:J_1,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'gzm')
      CALL CHECK3B(S0M(:,J_0:J_1,:) ,1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'s0m')
      CALL CHECK3B(SXMO(:,J_0:J_1,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'sxm')
      CALL CHECK3B(SYMO(:,J_0:J_1,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'sym')
      CALL CHECK3B(SZMO(:,J_0:J_1,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'szm')
      CALL CHECK3B(UO(:,J_0:J_1,:)  ,1,IM,J_0,J_1,NJPOL,LMO,
     &     SUBR,'uo ')
      CALL CHECK3B(VO(:,J_0STG:J_1STG,:),1,IM,J_0STG,J_1STG,0,LMO,
     &     SUBR,'vo ')
#ifdef TRACERS_OCEAN
      CALL CHECK4B(TRMO(:,J_0:J_1,:,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     tracerlist%getsize(), SUBR,'trm')
      CALL CHECK4B(TXMO(:,J_0:J_1,:,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     tracerlist%getsize(), SUBR,'txm')
      CALL CHECK4B(TYMO(:,J_0:J_1,:,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     tracerlist%getsize(), SUBR,'tym')
      CALL CHECK4B(TZMO(:,J_0:J_1,:,:),1,IM,J_0,J_1,NJPOL,LMO,
     &     tracerlist%getsize(), SUBR,'tzm')
#endif

C**** Check for variables out of bounds
      QCHECKO=.FALSE.
      DO J=J_0,J_1
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).gt.0.) THEN
C**** Check potential specific enthalpy/salinity
          DO L=1,LMM(I,J)
          GO1 = G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          SO1 = S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
          IF(GO1.lt.-10000. .or. GO1.gt.200000.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,GO=',I,J,L,GO1,TEMGS(GO1
     *           ,SO1)
            IF (GO1.lt.-20000. .or. GO1.gt.200000.) QCHECKO=.TRUE.
          END IF
          IF(SO1.gt.0.045 .or. SO1.lt.0.) THEN
            WRITE (6,*) 'After ',SUBR,': I,J,L,SO=',I,J,L,1d3*SO1
            IF ((SO1.gt.0.05 .or. SO1.lt.0.) .and. .not. (I == 47.and.
     *           J == 30)) QCHECKO=.TRUE.
          END IF
          END DO
C**** Check all ocean currents
          DO L = 1,LMO
            IF(ABS(UO(I,J,L)).gt.7. .or. ABS(VO(I,J,L)).gt.5) THEN
              WRITE (6,*) 'After ',SUBR,': I,J,L,UO,VO=',I,J,L,UO(I,J,L)
     *             ,VO(I,J,L)
              QCHECKO=.TRUE.
            END IF
          END DO
C**** Check first layer ocean mass
          IF(MO(I,J,1).lt.2000. .or. MO(I,J,1).gt.20000.) THEN
!!          IF (I == 47.and.(J == 33.or.J == 34)) GOTO 230 ! not Caspian
            WRITE (6,*) 'After ',SUBR,': I,J,MO=',I,J,MO(I,J,1)
!!          QCHECKO=.TRUE.
          END IF
C**** Check ocean salinity in each eighth box for the first layer
 230      SALIM = .045
!!        IF(JM == 46 .and. I == 47 .and. J == 30) GOTO 240 ! not Persia
!                                                         n Gulf   !.048
          IF(.5*ABS(SXMO(I,J,1))+.5*ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J
     *         ,1)).lt.S0M(I,J,1) .and.(.5*ABS(SXMO(I,J,1))+.5
     *         *ABS(SYMO(I,J,1))+BYRT3*ABS(SZMO(I,J,1))+S0M(I,J,1))
     *         /(MO(I,J,1)*DXYPO(J)).lt.SALIM)  GO TO 240
          WRITE (6,*) 'After ',SUBR,': I,J,S0,SX,SY,SZ=',I,J,
     *         1d3*S0M(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3*SXMO(I,J,1)/(MO(I
     *         ,J,1)*DXYPO(J)),1d3*SYMO(I,J,1)/(MO(I,J,1)*DXYPO(J)),1d3
     *         *SZMO(I,J,1)/(MO(I,J,1)*DXYPO(J))
!!        QCHECKO=.TRUE.
 240      CONTINUE
        END IF
      END DO
      END DO

#ifdef TRACERS_OCEAN
      do n=1,tracerlist%getsize()
        entry=>tracerlist%at(n)
C**** Check for negative tracers
        if (entry%t_qlimit) then
        do l=1,lmo
        do j=j_0,j_1
        do i=1,imaxj(j)
          if (l.le.lmm(i,j)) then
            if (trmo(i,j,l,n).lt.0) then
              write(6,*) "Neg Tracer in ocean after ",subr,i,j,l,
     *             entry%trname,trmo(i,j,l,n)
              QCHECKO=.true.
            end if
          end if
        end do
        end do
        end do
        end if
C**** Check conservation of water tracers in ocean
        if (entry%trname == 'Water') then
          errmax = 0. ; imax=1 ; jmax=J_0 ; lmax=1
          do l=1,lmo
          do j=j_0,j_1
          do i=1,imaxj(j)
            if (l.le.lmm(i,j)) then
              relerr=max(
     *             abs(trmo(i,j,l,n)-mo(i,j,l)*dxypo(j)+s0m(i,j,l)),
     *             abs(txmo(i,j,l,n)+sxmo(i,j,l)),
     *             abs(tymo(i,j,l,n)+symo(i,j,l)),
     *             abs(tzmo(i,j,l,n)+szmo(i,j,l)))/
     *             (mo(i,j,l)*dxypo(j)-s0m(i,j,l))
            if (relerr.gt.errmax) then
              imax=i ; jmax=j ; lmax=l ; errmax=relerr
            end if
            end if
          end do
          end do
          end do
          print*,"Relative error in ocean fresh water mass after ",subr
     *         ,":",imax,jmax,lmax,errmax,trmo(imax,jmax,lmax,n),mo(imax
     *         ,jmax,lmax)*dxypo(jmax)-s0m(imax,jmax,lmax),txmo(imax
     *         ,jmax,lmax,n),-sxmo(imax,jmax,lmax),tymo(imax,jmax,lmax,n
     *         ),-symo(imax,jmax ,lmax),tzmo(imax,jmax,lmax,n),
     *         -szmo(imax,jmax,lmax)
        end if
      end do
#endif

      IF (QCHECKO)
     &     call stop_model("QCHECKO: Ocean Variables out of bounds",255)

      END IF
C****
      IF (AM_I_ROOT()) CALL CHECKOST(SUBR)
C****
      END SUBROUTINE CHECKO


      Subroutine CONSERV_OMS (OMASS)
C****
!@sum   CONSERV_OMS calculates zonal ocean mass (kg/m^2) on ocean grid
!@auth  Gavin Schmidt/G. Russell/D. Gueyffier
!@ver   2009/04/23
C****
      Use OCEAN,            Only: IMO=>IM,JMO=>JM,oXYP,LMM,MO,imaxj
      Use STRAITS,          Only: NMST,IST,JST,LMST,MMST
      Use OCEANR_DIM,       Only: oGRID
      USE DOMAIN_DECOMP_1D, Only: GETDomainBounds
      Implicit None
!@var aOMASS zonal ocean mass per whole latitude band area (kg/m^2)
      Real*8    :: OMASS(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      Integer*4 :: J_0,J_1,I,J,N
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
      call getDomainBounds(oGRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE  )

      Do J=J_0,J_1
         Do I=1,IMAXJ(J)
            OMASS(I,J) = Sum(MO(I,J,:LMM(I,J)))
         EndDo
      EndDo

      IF(HAVE_SOUTH_POLE) OMASS(2:IMO,1)   = OMASS(1,1)
      IF(HAVE_NORTH_POLE) OMASS(2:IMO,JMO) = OMASS(1,JMO)

C**** Include ocean mass of straits
      Do N=1,NMST
         I = IST(N,1)  ;  J = JST(N,1)
         If (J >= J_0 .and. J <= J_1)  OMASS(I,J) = OMASS(I,J) +
     +      .5 * Sum(MMST(1:LMST(N),N)) / oXYP(I,J)
         I = IST(N,2)  ;  J = JST(N,2)
         If (J >= J_0 .and. J <= J_1)  OMASS(I,J) = OMASS(I,J) +
     +      .5 * Sum(MMST(1:LMST(N),N)) / oXYP(I,J)
      EndDo
      End Subroutine CONSERV_OMS


      Subroutine CONSERV_OSL (OSALT)
C****
!@sum   CONSERV_OSL calculates ocean salt (kg/m^2) on ocean grid
!@auth  Gavin Schmidt/G. Russell/D. Gueyffier
!@ver   2009/04/23
C****
      USE OCEAN, only : IMO=>IM,JMO=>JM,oXYP,LMM, S0M, imaxj
      Use STRAITS,       Only: NMST,IST,JST,LMST, S0MST
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds
      Use OCEANR_DIM,    only: oGRID
      Implicit None
!@var OSALT zonal ocean salt per whole latitude band area (kg/m^2)
      Real*8    :: OSALT(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      Integer*4 :: J_0,J_1,I,J,N
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
      call getDomainBounds(oGRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE  )

      Do J=J_0,J_1
         Do I=1,IMAXJ(J)
            OSALT(I,J) = Sum(S0M(I,J,:LMM(I,J)))/oXYP(I,J)
         EndDo
      EndDo

      IF(HAVE_SOUTH_POLE) OSALT(2:IMO,1)   = OSALT(1,1)
      IF(HAVE_NORTH_POLE) OSALT(2:IMO,JMO) = OSALT(1,JMO)

C**** Include ocean salt of straits
      Do N=1,NMST
         I = IST(N,1)  ;  J = JST(N,1)
         If (J >= J_0 .and. J <= J_1)  OSALT(I,J) = OSALT(I,J) +
     +      .5 * Sum(S0MST(1:LMST(N),N)) / oXYP(I,J)
         I = IST(N,2)  ;  J = JST(N,2)
         If (J >= J_0 .and. J <= J_1)  OSALT(I,J) = OSALT(I,J) +
     +      .5 * Sum(S0MST(1:LMST(N),N)) / oXYP(I,J)
      EndDo
      End Subroutine CONSERV_OSL


      Subroutine CONSERV_OCE (OCEANE)
C****
!@sum   CONSERV_OCE calculates ocean potential enthalpy (J/m^2)
!@auth  Gavin Schmidt/G. Russell/D. Gueyffier
!@ver   2009/04/23
C****
      USE OCEAN, only : IMO=>IM,JMO=>JM,LMM, oXYP,G0M, imaxj
      Use STRAITS,       Only: NMST,IST,JST,LMST, G0MST
      USE DOMAIN_DECOMP_1D, only :  GETDomainBounds
      Use OCEANR_DIM,    Only: oGRID
      Implicit None
!@var aOCEANE zonal ocean potential enthalpy per band area (J/m^2)
      Real*8    :: OCEANE(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      Integer*4 :: J_0,J_1,I,J,N
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
      call getDomainBounds(oGRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE  )

      Do J=J_0,J_1
         Do I=1,IMAXJ(J)
          OCEANE(I,J) = Sum(G0M(I,J,:LMM(I,J)))/oXYP(I,J)
        EndDo
      EndDo

      IF(HAVE_SOUTH_POLE) OCEANE(2:IMO,1)   = OCEANE(1,1)
      IF(HAVE_NORTH_POLE) OCEANE(2:IMO,JMO) = OCEANE(1,JMO)

C**** Include ocean potential enthalpy of straits
      Do N=1,NMST
         I = IST(N,1)  ;  J = JST(N,1)
         If (J >= J_0 .and. J <= J_1)  OCEANE(I,J) = OCEANE(I,J) +
     +      .5 * Sum(G0MST(1:LMST(N),N)) / oXYP(I,J)
         I = IST(N,2)  ;  J = JST(N,2)
         If (J >= J_0 .and. J <= J_1)  OCEANE(I,J) = OCEANE(I,J) +
     +      .5 * Sum(G0MST(1:LMST(N),N)) / oXYP(I,J)
      EndDo
      End Subroutine CONSERV_OCE


      Subroutine CONSERV_OKE (OKE)
C****
!@sum   CONSERV_OKE calculates ocean kinetic energy
!@auth  Gavin Schmidt/G. Russell/D. Gueyffier
!@ver   2009/04/23
C****
      USE OCEAN, only : IMO=>IM,JMO=>JM, IVSPO=>IVSP,IVNPO=>IVNP
     *     , oXYP, MO,UO,VO, LMOM=>LMM
      USE DOMAIN_DECOMP_1D, only : GETDomainBounds, SOUTH, HALO_UPDATE
      Use OCEANR_DIM,    Only: oGRID
      Implicit None
!@var OKE zonal ocean kinetic energy (J/m^2)
      Real*8    :: OKE(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
      Integer*4 :: J_0,J_1,I,J,L,Ip1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
      call getDomainBounds(oGRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE  )
      Call HALO_UPDATE (oGRID, VO, From=South)

      Do J=max(J_0,2),min(J_1,JMO-1)
         do I=2,IMO
            OKE(I,J) = 0.0
            Do L=1,LMOM(I,J)
               OKE(I,J) = OKE(I,J) +
     &               MO(I,J,L)*(UO(I-1,J,L)**2+UO(I,J,L)**2
     &                        + VO(I,J-1,L)**2+VO(I,J,L)**2 )
            EndDo
         EndDo
c*     I=1 here
         OKE(1,J) = 0.0
         Do L=1,LMOM(1,J)
            OKE(1,J) = OKE(1,J) +
     &           MO(1,J,L)*(UO(IMO,J,L)**2+UO(1,J,L)**2
     &                    + VO(1,J-1,L)**2+VO(1,J,L)**2 )
         EndDo
      EndDo

C**** OKE at south pole
      If (HAVE_SOUTH_POLE)  Then
         OKE(1,1) = 0.
         Do L=1,LMOM(1,1)
            OKE(1,1) = OKE(1,1)+  MO(1,1,L) *
     &                (1.5*(UO(IMO,1,L)**2 + UO(IVSPO,1,L)**2)
     &               + Sum(VO(:,1,L)**2)/IMO)
         EndDo
         OKE(2:IMO,1)=OKE(1,1)
      EndIf

C**** OKE at north pole
      If (HAVE_NORTH_POLE)  Then
         OKE(1,JMO) = 0.
         Do L=1,LMOM(1,JMO)
            OKE(1,JMO) = MO(1,JMO,L) *
     &                   (1.5*(UO(IMO,JMO,L)**2 + UO(IVNPO,JMO,L)**2)
     &                   + Sum(VO(:,JMO-1,L)**2)/IMO)
        EndDo
        OKE(2:IMO,JMO)=OKE(1,JMO)
      EndIf
      OKE(:,J_0:J_1) = OKE(:,J_0:J_1)*.25

      End Subroutine CONSERV_OKE

      Subroutine CONSERV_OAM (OAM)
C****
!@sum   CONSERV_OAM calculates ocean angular momentum (kg/s)
!@auth  Gavin Schmidt/G. Russell/D. Gueyffier
!@ver   2009/04/23
C****
      Use CONSTANT,      Only: RADIUS,OMEGA
      USE OCEAN, only : IMO=>IM,JMO=>JM, IVSPO=>IVSP,IVNPO=>IVNP
     *     , oDXYP=>DXYPO, LMOM=>LMM,LMOU=>LMU, MO,UO
     *     , oCOSQ=>COSQ, oCOSM=>COSM

      USE DOMAIN_DECOMP_1D, only : GETDomainBounds
      Use OCEANR_DIM,    Only: oGRID
      Implicit None
      Real*8    :: OAM(IMO,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     *             OMASS,UMLx2
      Integer*4 :: J_0,J_1,I,J,L,Ip1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
      call getDomainBounds(oGRID, J_STRT=J_0, J_STOP=J_1,
     &     HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &     HAVE_NORTH_POLE = HAVE_NORTH_POLE  )

      Do J=max(J_0,2),min(J_1,JMO-1)
         Do I=2,IMO
            OMASS = Sum(MO(I,J,1:LMOM(I,J)))
            UMLx2 = 0.
            Do L=1,LMOU(I,J)
               UMLx2 = MO(I,J,L)*(UO(I-1,J,L)+UO(I,J,L))
            EndDo

            OAM(I,J) = (.5*UMLx2*oCOSM(J)+RADIUS*OMEGA*oCOSQ(J)*OMASS)
     &                 *RADIUS
         EndDo
c     I=1 here
         OMASS = Sum(MO(1,J,1:LMOM(1,J)))
         UMLx2 = 0.
         Do L=1,LMOU(1,J)
            UMLx2 = MO(1,J,L)*(UO(IMO,J,L)+UO(1,J,L))
         EndDo

         OAM(1,J) = (.5*UMLx2*oCOSM(J) + RADIUS*OMEGA*oCOSQ(J)*OMASS)*
     *        RADIUS
      EndDo

      If (HAVE_SOUTH_POLE) then
         OAM(1,1) = RADIUS*RADIUS*OMEGA*oCOSQ(1)*
     *              Sum(MO(1,1,:LMOM(1,1)))
         OAM(2:IMO,1) = OAM(1,1)
      EndIf
      If (HAVE_NORTH_POLE) then
         OAM(1,JMO) = RADIUS*RADIUS*OMEGA*oCOSQ(JMO)*
     *               Sum(MO(1,JMO,:LMOM(1,JMO)))
         OAM(2:IMO,JMO) = OAM(1, JMO)
      EndIf

      End Subroutine CONSERV_OAM

      Subroutine OVtoM (MM,UM,VM)
C****
C**** OVtoM converts density and velocity in concentration units
C**** to mass and momentum in mass units
C**** Input:  MO (kg/m^2), UO (m/s), VO (m/s)
C**** Define: VOSP from UO(IVSP,1), VONP from UO(IVNP,JM)
C**** Output: MM (kg), UM (kg*m/s), VM (kg*m/s)
C****
      Use OCEAN, Only: IM,JM,LMO, IVSP,IVNP, DXYP=>DXYPO, COSU,SINU,
     *                 COSI=>COSIC,SINI=>SINIC, MO,UO,VO, VONP

!      Use DOMAIN_DECOMP_1D, Only: GRID, HALO_UPDATE, NORTH
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Intent(Out),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: MM,UM,VM
      Integer*4 I,J,L, J1,JN,J1P,JNP,JNQ,JNR
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
      JNR = Min(JN+1,JM-1)  !    5      9     JM-1   Halo, Exclude NP
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID, MO, FROM=NORTH)
C****
      DO 50 L=1,LMO
C**** Define VOSP and VONP
C     If (QSP)  VOSP(L) = UO(IVSP,1 ,L)
      If (QNP)  VONP(L) = UO(IVNP,JM,L)
C**** Convert density to mass
      DO 10 J=J1,JNR
   10 MM(:,J,L) = MO(:,J,L)*DXYP(J)
C     If (QSP)  MM(1,1 ,L) = MO(1,1 ,L)*DXYP(1)
      If (QNP)  MM(1,JM,L) = MO(1,JM,L)*DXYP(JM)
C**** Convert U velocity to U momentum
      DO 30 J=J1P,JNP
      DO 20 I=1,IM-1
   20 UM(I ,J,L) = UO(I ,J,L)*(MM(I,J,L)+MM(I+1,J,L))*.5
   30 UM(IM,J,L) = UO(IM,J,L)*(MM(1,J,L)+MM(IM ,J,L))*.5
C**** Convert V velocity to V momentum
      DO 40 J=J1P,JNQ
   40 VM(:,J,L) = VO(:,J,L)*(MM(:,J,L)+MM(:,J+1,L))*.5
C     If (QSP)  VM(:, 1  ,L) = VO(:, 1  ,L)*(MM(:, 2  ,L)+MM(1, 1,L))*.5
      If (QNP)  VM(:,JM-1,L) = VO(:,JM-1,L)*(MM(:,JM-1,L)+MM(1,JM,L))*.5
C**** Convert U and V velocity to U and V momentum at poles
C     If (QSP)  UM(IM,1 ,L) = UO(IM,1 ,L)*MM(1,1 ,L)
      If (QNP)  UM(IM,JM,L) = UO(IM,JM,L)*MM(1,JM,L)
C     If (QSP)  UM(IVSP,1 ,L) =   VOSP(L)*MM(1,1 ,L)
      If (QNP)  UM(IVNP,JM,L) =   VONP(L)*MM(1,JM,L)
C**** Define MO, UO and VO at poles
C     If (QSP)  Then
C       MO(:,1 ,L) = MO(1 ,1 ,L)
C       UO(:,1 ,L) = UO(IM,1 ,L)*COSU(:) - VOSP(L)*SINU(:)
C       VO(:,0 ,L) = VOSP(L)*COSI(:) + UO(IM,1 ,L)*SINI(:)  ;  EndIf
      If (QNP)  Then
        MO(:,JM,L) = MO(1 ,JM,L)
        UO(:,JM,L) = UO(IM,JM,L)*COSU(:) + VONP(L)*SINU(:)
        VO(:,JM,L) = VONP(L)*COSI(:) - UO(IM,JM,L)*SINI(:)  ;  EndIf
   50 Continue
      Return
      EndSubroutine OVtoM

      Subroutine OMtoV (MM,UM,VM)
C****
C**** OMtov converts mass and momentum in mass units
C**** to density and velocity in concentration units
C**** Input:  MM (kg), UM (kg*m/s), VM (kg*m/s)
C**** Output: MO (kg/m^2), UO (m/s), VO (m/s)
C**** Define: VOSP from UM(IVSP,1)/MM, VONP from UM(IVNP,JM)/MM
C****
      Use OCEAN, Only: IM,JM,LMO,IVSP,IVNP,
     *                  LMOM=>LMM, LMOU=>LMU, LMOV=>LMV,
     *                  COSU,SINU, COSI=>COSIC,SINI=>SINIC,
     *                  zDXYP=>BYDXYPO, MO,UO,VO, VONP

!      Use DOMAIN_DECOMP_1D, Only: GRID, HALO_UPDATE, NORTH
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8, !!! Intent(IN), (except for HALO_UPDATEs)
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: MM,UM,VM
      Integer*4 I,J,L, J1,JN,J1P,JNP,JNQ
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID, MM, FROM=NORTH)
C**** Convert mass to density
      Do 10 J=J1P,JNP
      Do 10 I=1,IM
      Do 10 L=1,LMOM(I,J)
   10 MO(I,J,L) = MM(I,J,L)*zDXYP(J)
C     If (QSP)  Then
C       Do L=1,LMOM(1,1)
C  11   MO(:,1,L) = MM(1,1,L)*zDXYP(1)  ;  EndIf
      If (QNP)  Then
        Do 12 L=1,LMOM(1,JM)
   12   MO(:,JM,L) = MM(1,JM,L)*zDXYP(JM)  ;  EndIf
C**** Convert U momentum to U velocity
      Do 30 J=J1P,JNP
      Do 20 I=1,IM-1
      Do 20 L=1,LMOU(I,J)
   20 UO(I,J,L) = UM(I,J,L)*2/(MM(I,J,L)+MM(I+1,J,L))
      Do 30 L=1,LMOU(IM,J)
   30 UO(IM,J,L) = UM(IM,J,L)*2/(MM(IM,J,L)+MM(1,J,L))
C**** Convert V momentum to V velocity
      Do 40 J=J1P,JNQ
      Do 40 I=1,IM
      Do 40 L=1,LMOV(I,J)
   40 VO(I,J,L) = VM(I,J,L)*2/(MM(I,J,L)+MM(I,J+1,L))
C     If (QSP)  Then
C       Do 41 I=1,IM
C       Do 41 L=1,LMOV(I,1)
C  41   VO(I,1,L) = VM(I,1,L)*2/(MM(1,1,L)+MM(I,2,L))  ;  EndIf
      If (QNP)  Then
        Do 42 I=1,IM
        Do 42 L=1,LMOV(I,JM-1)
   42   VO(I,JM-1,L) = VM(I,JM-1,L)*2/(MM(1,JM,L)+MM(I,JM-1,L))  ; EndIf
C**** Convert momentum to velocity at south pole, define UO and VO
C     If (QSP)  Then
C       Do 50 L=1,LMOM(1,1)
C       UO(IM,1,L) = UM(IM,1,L)/MM(1,1,L)
C       VOSP(L)  = UM(IVSP,1,L)/MM(1,1,L)
C       UO(:,1,L) = UO(IM,1,L)*COSU(:) - VOSP(L)*SINU(:)
C       VO(:,0,L) = VOSP(L)*COSI(:) + UO(IM,1,L)*SINI(:)
C  50   Continue  ;  EndIf
C**** Convert momentum to velocity at north pole, define UO and VO
      If (QNP)  Then
        Do 60 L=1,LMOM(1,JM)
        UO(IM,JM,L) = UM(IM,JM,L)/MM(1,JM,L)
        VONP(L)   = UM(IVNP,JM,L)/MM(1,JM,L)
        UO(:,JM,L) = UO(IM,JM,L)*COSU(:) + VONP(L)*SINU(:)
        VO(:,JM,L) = VONP(L)*COSI(:) - UO(IM,JM,L)*SINI(:)
   60   Continue  ;  EndIf
      Return
      EndSubroutine OMtoV

      Subroutine OFLUX (NS,MM,QEVEN)
C****
C**** OFLUX calculates the fluid fluxes
C**** Input:  M (kg/m^2), U (m/s), V (m/s)
C**** Output: MU (kg/s) = DY * U * M
C****         MV (kg/s) = DX * V * M
C****         MW (kg/s) = MW(L-1) + CONV-CONVs*dZO/ZOE + MM-MMs*dZO/ZOE
C****       CONV (kg/s) = MU-MU + MV-MV
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, ZOE=>ZE,dZO, DTO,
     *                 DXYP=>DXYPO,DYP=>DYPO,DXV=>DXVO, MO,UO,VO
      Use OCEAN_DYN, Only: MU,MV,MW, SMU,SMV,SMW, CONV

!      Use DOMAIN_DECOMP_1D, Only: GRID, HALO_UPDATE, NORTH,SOUTH
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH,SOUTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Integer*4,Intent(In) :: NS     !  decrementing leap-frog step
      Logical*4,Intent(In) :: QEVEN  !  whether called from even step
      Real*8   ,Intent(In),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: MM
C****
      Real*8    zNSxDTO, MVS,dMVS(IM),dMVSm, MVN,dMVN(IM),dMVNm,
     *          CONVs,MMs
      Integer*4 I,J,L,LM,IMAX, J1,JN,J1P,J1H,JNP,JNQ
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID, MO, FROM=NORTH+SOUTH)
      Call HALO_UPDATE (GRID, VO, FROM=SOUTH)
      zNSxDTO = 1 / (NS*DTO)
C****
C**** Compute fluid fluxes for the C grid
C****
C**** Smooth the West-East velocity near the poles
      Do 110 L=1,LMO
      Do 110 J=J1P,JNP
  110 MU(:,J,L) = UO(:,J,L)
      Call OPFIL (MU)
C**** Compute MU, the West-East mass flux, at non-polar points
      DO 430 L=1,LMO
      DO 130 J=J1P,JNP
      DO 120 I=1,IM-1
  120 MU(I ,J,L) = .5*DYP(J)*MU(I ,J,L)*(MO(I,J,L)+MO(I+1,J,L))
  130 MU(IM,J,L) = .5*DYP(J)*MU(IM,J,L)*(MO(1,J,L)+MO(IM ,J,L))
C**** Compute MV, the South-North mass flux
      DO 210 J=J1H,JNP
  210 MV(:,J,L) = .5*DXV(J)*VO(:,J,L)*(MO(:,J,L)+MO(:,J+1,L))
C****
C**** Compute MU so that CONV is identical for each polar triangle
C****
C     If (QSP) then
C       MVS = Sum(MV(:,1,L)) / IM
C       dMVS(1) = 0
C       Do 310 I=2,IM
C 310   dMVS(I) = dMVS(I-1) + (MV(I,1,L)-MVS)
C       dMVSm   = Sum(dMVS(:)) / IM
C       MU(:,1,L) = dMVSm - dMVS(:)  ;  EndIf
      If (QNP) then
        MVN = Sum(MV(:,JM-1,L)) / IM
        dMVN(1) = 0
        Do 320 I=2,IM
  320   dMVN(I) = dMVN(I-1) + (MV(I,JM-1,L)-MVN)
        dMVNm   = Sum(dMVN(:)) / IM
        MU(:,JM,L) = dMVN(:) - dMVNm  ;  EndIf
C****
C**** Compute horizontal fluid convergence at non-polar points
C****
      Do 420 J=J1P,JNP
      Do 410 I=2,IM
  410 CONV(I,J,L) = MU(I-1,J,L)-MU(I,J,L) + (MV(I,J-1,L)-MV(I,J,L))
  420 CONV(1,J,L) = MU(IM ,J,L)-MU(1,J,L) + (MV(1,J-1,L)-MV(1,J,L))
C     If (QSP)  CONV(1,1 ,L) = - MVS
      If (QNP)  CONV(1,JM,L) =   MVN
  430 Continue
C****
C**** Compute vertically integrated column convergence and mass
C****
      Do 630 J=J1,JN
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 630 I=1,IMAX
      LM=1
      If (LMOM(I,J) <= 1)  GoTo 620
      LM=LMOM(I,J)
      CONVs = Sum(CONV(I,J,1:LM)) / ZOE(LM)
      MMs   = Sum(  MM(I,J,1:LM)) / ZOE(LM)
C****
C**** Compute MW, the downward fluid flux
C****
      MW(I,J,1) = CONV(I,J,1) - CONVs*dZO(1) +
     +           (  MM(I,J,1) -   MMs*dZO(1)) * zNSxDTO
      Do 610 L=2,LM-1
  610 MW(I,J,L) = CONV(I,J,L) - CONVs*dZO(L) + MW(I,J,L-1) +
     +           (  MM(I,J,L) -   MMs*dZO(L)) * zNSxDTO
  620 MW(I,J,LM:LMO-1) = 0
  630 Continue
C****
C**** Sum mass fluxes to be used for advection of tracers
C****
      If (.not.QEVEN)  Return
      Do 710 L=1,LMO
      SMU(:,J1 :JN ,L) = SMU(:,J1 :JN ,L) + MU(:,J1 :JN ,L)
      SMV(:,J1H:JNP,L) = SMV(:,J1H:JNP,L) + MV(:,J1H:JNP,L)
      If (L==LMO)  GoTo 710
      SMW(:,J1P:JNP,L) = SMW(:,J1P:JNP,L) + MW(:,J1P:JNP,L)
  710 Continue
C     If (QSP)  SMW(1,1 ,:) = SMW(1,1 ,:) + MW(1,1 ,:)
      If (QNP)  SMW(1,JM,:) = SMW(1,JM,:) + MW(1,JM,:)
      Return
      EndSubroutine OFLUX

      Subroutine OPFIL (X)
C****
C**** OPFIL smoothes X in the zonal direction by reducing coefficients
C**** of its Fourier series for high wave numbers near the poles.
C**** The ocean polar filter here works with AVR4X5LD.Z12.gas1 .
C****
      Use CONSTANT, Only: TWOPI
      Use OCEAN, Only: IM,JM,LMO,J1O, DLON, DXP=>DXPO, DYP=>DYPO,
     *  JMPF=>J40S  !  greatest J in SH where polar filter is applied
      Use FILEMANAGER, Only: OPENUNIT, CLOSEUNIT

!      Use DOMAIN_DECOMP_1D, Only: GRID
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
C****
      Real*8,Intent(InOut) ::
     *  X(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO)
      Integer*4,Parameter :: IMz2=IM/2  !  highest possible wave number
      Integer*4,Save ::
     *  JMEX,    !  exclude unfiltered latitudes, store J=JM-1 in JMEX
     *  NBASM,   !  maximum number of ocean basins at any latitude
     *  INDM,    !  number of entries in array REDUCO
     *  NMIN(JM) !  minimum wave number for Fourier smoothing
      Real*8,Save :: SMOOTH(IMz2,JM)  !  multiplies Fourier coefficient
C****
      Integer*2,Save,Allocatable,Dimension(:,:) ::
     *  NBAS    !  number of ocean basins for given layer and latitude
      Integer*2,Save,Allocatable,Dimension(:,:,:) ::
     *  IMINm1, !  western most cell in basin minus 1
     *  IWIDm1  !  number of cells in basin minus 1
      Integer*4,Save,Allocatable,Dimension(:,:) ::
     *  INDEX   !  index to concatenated matricies
      Real*4,Save,Allocatable,Dimension(:) ::
     *  REDUCO  !  concatenation of reduction matricies
      Real*4,Allocatable,Dimension(:) ::
     *  REDUCO_glob  !  concatenation of reduction matricies
C****
      Character*80 TITLE
      Integer*4 I,I1,INDX,IWm2,IWm1, J,JA,JX, K,L, N,NB, IU_AVR
      integer :: indx_min,indx_max
      Real*8 AN(0:IMz2), !  Fourier cosine coefficients
     *       BN(0:IMz2), !  Fourier sine coefficients
     *          Y(IM*2), !  original copy of X that wraps around IDL
     *       DRAT, REDUC
C****
      Integer*4,Save :: IFIRST=1, J1,JN,J1P,JNP
      If (IFIRST == 0)  GoTo 100
      IFIRST = 0
      JMEX = 2*JMPF-1
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
C**** Calculate SMOOTHing factor for longitudes without land cells
      Do 30 J=J1O,JM-1
      DRAT = DXP(J)/DYP(3)
      Do 20 N=IMz2,1,-1
      SMOOTH(N,J) = IMz2*DRAT/N
   20 If (SMOOTH(N,J) >= 1)  GoTo 30
C     N = 0
   30 NMIN(J) = N+1
C**** Read in reduction contribution matrices from disk.  Only keep
C**** the matrices needed for the latitudes on this processor.
      Call OPENUNIT ('AVR',IU_AVR,.True.,.True.)
      Read (IU_AVR) TITLE,NBASM,INDM
      Allocate (IMINm1(NBASM,LMO,J1O:JMEX), NBAS(LMO,J1O:JMEX),
     *          IWIDm1(NBASM,LMO,J1O:JMEX), INDEX(IM,2:JMPF),
     *          REDUCO_glob(INDM))
      Read (IU_AVR) TITLE,NBAS,IMINm1,IWIDm1,INDEX,REDUCO_glob
      Call CLOSEUNIT (IU_AVR)
      Write (6,*) 'Read from unit',IU_AVR,': ',TITLE
      indx_min = indm+1; indx_max = -1
      do J=Max(J1O,J1P),JNP
        JX=J  ;  If(J > JMPF) JX=J+2*JMPF-JM
        JA=J  ;  If(J > JMPF) JA=JM+1-J
        If (JA > JMPF)  cycle   !  skip latitudes J=JMPF+1,JM-JMPF
        do L=1,LMO
          If (IWIDm1(1,L,JX) >= IM)  cycle
          do NB=1,NBAS(L,JX)
            IWm1 = IWIDm1(NB,L,JX)
            indx_min = min(indx_min,INDEX(IWm1,JA)+1)
            indx_max = max(indx_max,INDEX(IWm1,JA)+IWm1**2)
          enddo
        enddo
      enddo
      if(indx_min.le.indx_max) then
        allocate(reduco(indx_min:indx_max))
        reduco(indx_min:indx_max) = reduco_glob(indx_min:indx_max)
      endif
      deallocate(reduco_glob)
 100  CONTINUE
C****
C**** Loop over J and L.  JX = eXclude unfiltered latitudes
C****                     JA = Absolute latitude
C****
      Do 410 J=Max(J1O,J1P),JNP
      JX=J  ;  If(J > JMPF) JX=J+2*JMPF-JM
      JA=J  ;  If(J > JMPF) JA=JM+1-J
      If (JA > JMPF)  GoTo 410  !  skip latitudes J=JMPF+1,JM-JMPF
      Do 400 L=1,LMO
      If (IWIDm1(1,L,JX) >= IM)  GoTo 300
C****
C**** Land cells exist at this latitude and layer, loop over ocean
C**** basins.
C****
      Do 270 NB=1,NBAS(L,JX)
      I1   = IMINm1(NB,L,JX) + 1
      IWm2 = IWIDm1(NB,L,JX) - 1
      INDX = INDEX(IWm2+1,JA)
      If (I1+IWm2 > IM)  GoTo 200
C**** Ocean basin does not wrap around the IDL.
C**** Copy X to temporary array Y and filter X in place.
      Do 110 I=I1,I1+IWm2
  110 Y(I) = X(I,J,L)
      Do 140 I=I1,I1+IWm2
      REDUC = 0
      Do 130 K=I1,I1+IWm2
      INDX=INDX+1
  130 REDUC = REDUC + REDUCO(INDX)*Y(K)
  140 X(I,J,L) = X(I,J,L) - REDUC
      GoTo 270
C**** Ocean basin wraps around the IDL.
C**** Copy X to temporary array Y and filter X in place.
  200 Do 210 I=I1,IM
  210 Y(I) = X(I,J,L)
      Do 220 I=IM+1,I1+IWm2
  220 Y(I) = X(I-IM,J,L)
      Do 240 I=I1,IM
      REDUC = 0
      Do 230 K=I1,I1+IWm2
      INDX=INDX+1
  230 REDUC = REDUC + REDUCO(INDX)*Y(K)
  240 X(I,J,L) = X(I,J,L) - REDUC
      Do 260 I=1,I1+IWm2-IM
      REDUC = 0
      Do 250 K=I1,I1+IWm2
      INDX=INDX+1
  250 REDUC = REDUC + REDUCO(INDX)*Y(K)
  260 X(I,J,L) = X(I,J,L) - REDUC
  270 Continue
      GoTo 400
C****
C**** No land cells at this latitude and layer,
C**** perform standard polar filter
C****
  300 Call OFFT (X(1,J,L),AN,BN)
      Do 310 N=NMIN(J),IMz2-1
      AN(N) = AN(N)*SMOOTH(N,J)
  310 BN(N) = BN(N)*SMOOTH(N,J)
      AN(IMz2) = AN(IMz2)*SMOOTH(IMz2,J)
      Call OFFTI (AN,BN,X(1,J,L))
C****
  400 Continue
  410 Continue
      Return
      EndSubroutine OPFIL

      Subroutine OADVM (MM2,MM0,DT)
C****
C**** OADVM calculates the updated mass
C**** Input:  MM0 (kg), DT (s), CONV (kg/s), MW (kg/s)
C**** Output: MM2 (kg) = MM0 + DT*[CONV + MW(L-1) - MW(L)]
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM
      Use OCEAN_DYN, Only: MW, CONV

!      Use DOMAIN_DECOMP_1D, Only: GRID
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Intent(Out):: MM2(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO)
      Real*8,Intent(In) :: MM0(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO)
     *                    ,DT
      Integer*4 I,J,L,LM,IMAX, J1,JN
C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
C****
C**** Compute the new mass MM2
C****
      Do 20 J=J1,JN
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 20 I=1,IMAX
      If (LMOM(I,J)==0)  GoTo 20
      If (LMOM(I,J)==1)  Then
        MM2(I,J,1) = MM0(I,J,1) + DT*CONV(I,J,1)
        GoTo 20  ;  EndIf
      LM = LMOM(I,J)
      MM2(I,J,1) = MM0(I,J,1) + DT*(CONV(I,J,1)-MW(I,J,1))
      Do 10 L=2,LM-1
   10 MM2(I,J,L) = MM0(I,J,L) + DT*(CONV(I,J,L) + MW(I,J,L-1)-MW(I,J,L))
      MM2(I,J,LM) = MM0(I,J,LM) + DT*(CONV(I,J,LM) + MW(I,J,LM-1))
   20 Continue
      Return
      EndSubroutine OADVM

      Subroutine OADVV (UM2,VM2,UM0,VM0,DT1)
C****
C**** OADVV advects oceanic momentum (with coriolis force)
C**** Input:  MO (kg/m^2), UO (m/s), VO (m/s) = from odd solution
C****         MU (kg/s), MV (kg/s), MW (kg/s) = fluid fluxes
C**** Output: UM2 (kg*m/s) = UM0 + DT*(MU*U-MU*U + MV*U-MV*U + M*CM*V)
C****         VM2 (kg*m/s) = VM0 + DT*(MU*V-MU*V + MV*V-MV*V - M*CM*U)
C****
      Use CONSTANT, Only: RADIUS,OMEGA
      Use OCEAN, Only: IM,JM,LMO, IVSP,IVNP, DXV=>DXVO,
     *                 COSU,SINU, SINxY,TANxY, MO,UO,VO
      Use OCEAN_DYN, Only: MU,MV,MW

!      Use DOMAIN_DECOMP_1D, Only: GRID, HALO_UPDATE, NORTH,SOUTH
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH,SOUTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Intent(Out),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: UM2,VM2
      Real*8,Intent(In),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: UM0,VM0
      Real*8,Intent(In) :: DT1
C****
      Real*8 DT2,DT4,DT6,DT12
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     *       DUM,DVM
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     *       UMVc,UMVw,UMVe
      Real*8 UMU,VMV,VMUc,VMUs,VMUn, FLUX, dUMSP,dVMSP,dUMNP,dVMNP,
     *       USINYJ(IM),UTANUJ(IM)
      Integer*4 I,J,L, J1,JN,J1P,J1H,JNP,JNR
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNR = Min(JN+1,JM-1)  !    5      9     JM-1   Halo, Exclude NP
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
C     Call HALO_UPDATE (GRID, MO, FROM=NORTH)        !  copied in OFLUX
      Call HALO_UPDATE (GRID, UO, FROM=NORTH+SOUTH)
      Call HALO_UPDATE (GRID, VO, FROM=NORTH)
C     Call HALO_UPDATE (GRID, VO, FROM=SOUTH)        !  copied in OFLUX
      Call HALO_UPDATE (GRID, MU, FROM=NORTH)
      Call HALO_UPDATE (GRID, MV, FROM=NORTH)
C     Call HALO_UPDATE (GRID, MV, FROM=SOUTH)        !  recalc in OFLUX
      Call HALO_UPDATE (GRID, MW, FROM=NORTH)
C****
      DT2  = DT1/2
      DT4  = DT1/4
      DT6  = DT1/6
      Dt12 = DT1/12
C****
C**** Horizontal advection of momentum
C****
      Do 480 L=1,LMO
C**** Zero out momentum changes
      DUM(:,:,L) = 0
      DVM(:,:,L) = 0
C**** Contribution of eastward flux to U momentum
      Do 310 J=J1,JN
      UMU = DT4*(UO(1,J,L)+UO(IM,J,L))*(MU(1,J,L)+MU(IM,J,L))
      DUM(1 ,J,L) = DUM(1 ,J,L) + UMU
      DUM(IM,J,L) = DUM(IM,J,L) - UMU
      Do 310 I=2,IM
      UMU = DT4*(UO(I,J,L)+UO(I-1,J,L))*(MU(I,J,L)+MU(I-1,J,L))
      DUM(I  ,J,L) = DUM(I  ,J,L) + UMU
  310 DUM(I-1,J,L) = DUM(I-1,J,L) - UMU
C**** Contribution of northward center to U momentum
      Do 320 J=J1H,JNP
      UMVc(IM,J) = DT6*(UO(IM,J,L)+UO(IM,J+1,L))*(MV(IM,J,L)+MV(1,J,L))
      DUM(IM,J  ,L) = DUM(IM,J  ,L) - UMVc(IM,J)
      DUM(IM,J+1,L) = DUM(IM,J+1,L) + UMVc(IM,J)
      Do 320 I=1,IM-1
      UMVc(I,J) = DT6*(UO(I,J,L)+UO(I,J+1,L))*(MV(I,J,L)+MV(I+1,J,L))
      DUM(I,J  ,L) = DUM(I,J  ,L) - UMVc(I,J)
  320 DUM(I,J+1,L) = DUM(I,J+1,L) + UMVc(I,J)
C**** Contribution of northward corner fluxes to U momentum
      Do 330 J=J1H,JNP
      UMVe(1,J) = DT12*(UO(IM,J,L)+UO(1 ,J+1,L))*MV(1,J,L)
      UMVw(1,J) = DT12*(UO(1 ,J,L)+UO(IM,J+1,L))*MV(1,J,L)
      DUM(1 ,J  ,L) = DUM(1 ,J  ,L) - UMVw(1,J)
      DUM(IM,J  ,L) = DUM(IM,J  ,L) - UMVe(1,J)
      DUM(1 ,J+1,L) = DUM(1 ,J+1,L) + UMVe(1,J)
      DUM(IM,J+1,L) = DUM(IM,J+1,L) + UMVw(1,J)
      Do 330 I=2,IM
      UMVe(I,J) = DT12*(UO(I-1,J,L)+UO(I  ,J+1,L))*MV(I,J,L)
      UMVw(I,J) = DT12*(UO(I  ,J,L)+UO(I-1,J+1,L))*MV(I,J,L)
      DUM(I  ,J  ,L) = DUM(I  ,J  ,L) - UMVw(I,J)
      DUM(I-1,J  ,L) = DUM(I-1,J  ,L) - UMVe(I,J)
      DUM(I  ,J+1,L) = DUM(I  ,J+1,L) + UMVe(I,J)
  330 DUM(I-1,J+1,L) = DUM(I-1,J+1,L) + UMVw(I,J)
C**** Contribution of eastward center flux to V momentum
      Do 340 J=J1,JNP
      VMUc = DT6*(VO(IM,J,L)+VO(1,J,L))*(MU(IM,J,L)+MU(IM,J+1,L))
      DVM(IM,J,L) = DVM(IM,J,L) - VMUc
      DVM(1 ,J,L) = DVM(1 ,J,L) + VMUc
      Do 340 I=1,IM-1
      VMUc = DT6*(VO(I,J,L)+VO(I+1,J,L))*(MU(I,J,L)+MU(I,J+1,L))
      DVM(I  ,J,L) = DVM(I  ,J,L) - VMUc
  340 DVM(I+1,J,L) = DVM(I+1,J,L) + VMUc
C**** Contribution of eastward corner and northward fluxes to V momentum
C     If (QSP)  Then
C       J = 1
C       VMV  = DT4 *(VO(IM,1,L)+V0(IM,0,L))*MV(IM,1,L)
C       VMUn = DT12*(V0(IM,0,L)+VO(1 ,1,L))*MU(IM,1,L)
C       VMUs = DT12*(VO(IM,1,L)+V0(1 ,0,L))*MU(IM,1,L)
C       DVM(IM,1,L) = DVM(IM,1,L) - VMUs + VMV
C       DVM(1 ,1,L) = DVM(1 ,1,L) + VMUn
C       Do 350 I=1,IM-1
C       VMV  = DT4* (VO(I,1,L)+V0(I  ,0,L))*MV(I,1,L)
C       VMUn = DT12*(V0(I,0,L)+VO(I+1,1,L))*MU(I,1,L)
C       VMUs = DT12*(VO(I,1,L)+V0(I+1,0,L))*MU(I,1,L)
C       DVM(I  ,1,L) = DVM(I  ,1,L) - VMUs + VMV
C 350   DVM(I+1,1,L) = DVM(I+1,1,L) + VMUn  ;  EndIf
      Do 360 J=J1P,JNR
      VMV  = DT4 *(VO(IM,J  ,L)+VO(IM,J-1,L))*(MV(IM,J,L)+MV(IM,J-1,L))
      VMUn = DT12*(VO(IM,J-1,L)+VO(1 ,J  ,L))* MU(IM,J,L)
      VMUs = DT12*(VO(IM,J  ,L)+VO(1 ,J-1,L))* MU(IM,J,L)
      DVM(IM,J  ,L) = DVM(IM,J  ,L) - VMUs + VMV
      DVM(IM,J-1,L) = DVM(IM,J-1,L) - VMUn - VMV
      DVM(1 ,J  ,L) = DVM(1 ,J  ,L) + VMUn
      DVM(1 ,J-1,L) = DVM(1 ,J-1,L) + VMUs
      Do 360 I=1,IM-1
      VMV  = DT4* (VO(I,J  ,L)+VO(I  ,J-1,L))*(MV(I,J,L)+MV(I,J-1,L))
      VMUn = DT12*(VO(I,J-1,L)+VO(I+1,J  ,L))* MU(I,J,L)
      VMUs = DT12*(VO(I,J  ,L)+VO(I+1,J-1,L))* MU(I,J,L)
      DVM(I  ,J  ,L) = DVM(I  ,J  ,L) - VMUs + VMV
      DVM(I  ,J-1,L) = DVM(I  ,J-1,L) - VMUn - VMV
      DVM(I+1,J  ,L) = DVM(I+1,J  ,L) + VMUn
  360 DVM(I+1,J-1,L) = DVM(I+1,J-1,L) + VMUs
      If (QNP)  Then
C       J = JM
        VMV  = DT4 *(VO(IM,JM  ,L)+VO(IM,JM-1,L))*MV(IM,JM-1,L)
        VMUn = DT12*(VO(IM,JM-1,L)+VO(1 ,JM  ,L))*MU(IM,JM,L)
        VMUs = DT12*(VO(IM,JM  ,L)+VO(1 ,JM-1,L))*MU(IM,JM,L)
        DVM(IM,JM-1,L) = DVM(IM,JM-1,L) - VMUn - VMV
        DVM(1 ,JM-1,L) = DVM(1 ,JM-1,L) + VMUs
        Do 370 I=1,IM-1
        VMV  = DT4* (VO(I,JM  ,L)+VO(I  ,JM-1,L))*MV(I,JM-1,L)
        VMUn = DT12*(VO(I,JM-1,L)+VO(I+1,JM  ,L))*MU(I,JM,L)
        VMUs = DT12*(VO(I,JM  ,L)+VO(I+1,JM-1,L))*MU(I,JM,L)
        DVM(I  ,JM-1,L) = DVM(I  ,JM-1,L) - VMUn - VMV
  370   DVM(I+1,JM-1,L) = DVM(I+1,JM-1,L) + VMUs  ;  EndIf
C****
C**** Coriolis force and metric term
C****
C**** U component
C     If (QSP)  Then
C       Do 410 I=1,IM-1
C 410   dUM(I,1,L) = dUM(I,1,L) +
C    +    DT2*OMEGA*RADIUS*SINxY(1)*(MV(I,1,L)+MV(I+1,1,L)) +
C    +    TANxY(1)*(UMVw(I,1)+UMVc(I,1)+UMVe(I+1,1))*.5
C       dUM(IM,1,L) = dUM(IM,1,L) +
C    +    DT2*OMEGA*RADIUS*SINxY(1)*(MV(IM,1,L)+MV(1,1,L)) +
C    +    TANxY(1)*(UMVw(IM,1)+UMVc(IM,1)+UMVe(1,1))*.5  !  EndIf
      Do 420 J=J1P,JNP
      dUM(IM,J,L) = dUM(IM,J,L) + DT2*OMEGA*RADIUS*SINxY(J)*
     *    (MV(IM,J-1,L)+MV(1,J-1,L)+MV(IM,J,L)+MV(1,J,L)) +
     +  TANxY(J)*(UMVe(IM,J-1)+UMVc(IM,J-1)+UMVw(1,J-1) +
     +            UMVw(IM,J  )+UMVc(IM,J  )+UMVe(1,J  ))*.5
      Do 420 I=1,IM-1
  420 dUM(I,J,L) = dUM(I,J,L) + DT2*OMEGA*RADIUS*SINxY(J)*
     *    (MV(I,J-1,L)+MV(I+1,J-1,L)+MV(I,J,L)+MV(I+1,J,L)) +
     +  TANxY(J)*(UMVe(I,J-1)+UMVc(I,J-1)+UMVw(I+1,J-1) +
     +            UMVw(I,J  )+UMVc(I,J  )+UMVe(I+1,J  ))*.5
      If (QNP)  Then
        Do 430 I=1,IM-1
  430   dUM(I,JM,L) = dUM(I,JM,L) +
     +    DT2*OMEGA*RADIUS*SINxY(JM)*(MV(I,JM-1,L)+MV(I+1,JM-1,L)) +
     +    TANxY(JM)*(UMVe(I,JM-1)+UMVc(I,JM-1)+UMVw(I+1,JM-1))*.5
        dUM(IM,JM,L) = dUM(IM,JM,L) +
     +    DT2*OMEGA*RADIUS*SINxY(JM)*(MV(IM,JM-1,L)+MV(1,JM-1,L)) +
     +    TANxY(JM)*(UMVe(IM,JM-1)+UMVc(IM,JM-1)+UMVw(1,JM-1))*.5
        EndIf
C**** V component
      Do 450 J=J1,JNP
      Do 440 I=1,IM
      USINYJ(I) =  UO(I,J,L)*SINxY(J) + UO(I,J+1,L)*SINxY(J+1)
  440 UTANUJ(I) = (UO(I,J,L)*TANxY(J) + UO(I,J+1,L)*TANxY(J+1))*
     *            (UO(I,J,L)+UO(I,J+1,L))
      dVM(1,J,L) = dVM(1,J,L) - DT4*DXV(J)*(MO(1,J,L)+MO(1,J+1,L))*
     *  (OMEGA*RADIUS*(USINYJ(IM)+USINYJ(1)) +
     +  .25*(UTANUJ(IM)+UTANUJ(1)) - (TANxY(J)+TANxY(J+1))*
     *    (UO(IM,J,L)-UO(1,J,L))*(UO(IM,J+1,L)-UO(1,J+1,L))/12)
      Do 450 I=2,IM
  450 dVM(I,J,L) = dVM(I,J,L) - DT4*DXV(J)*(MO(I,J,L)+MO(I,J+1,L))*
     *  (OMEGA*RADIUS*(USINYJ(I-1)+USINYJ(I)) +
     +  .25*(UTANUJ(I-1)+UTANUJ(I)) - (TANxY(J)+TANxY(J+1))*
     *    (UO(I-1,J,L)-UO(I,J,L))*(UO(I-1,J+1,L)-UO(I,J+1,L))/12)
  480 Continue
C****
C**** Vertical advection of momentum
C****
C**** U component
      Do 560 J=J1,JN
      If (J==1)   GoTo 520
      If (J==JM)  GoTo 540
      Do 510 L=1,LMO-1
      FLUX = DT4*(MW(IM,J,L)+MW(1,J,L))*(UO(IM,J,L)+UO(IM,J,L+1))
      dUM(IM,J,L)   = dUM(IM,J,L)   - FLUX
      dUM(IM,J,L+1) = dUM(IM,J,L+1) + FLUX
      Do 510 I=1,IM-1
      FLUX = DT4*(MW(I,J,L)+MW(I+1,J,L))*(UO(I,J,L)+UO(I,J,L+1))
      dUM(I,J,L)   = dUM(I,J,L)   - FLUX
  510 dUM(I,J,L+1) = dUM(I,J,L+1) + FLUX
      GoTo 560
  520 Do 530 L=1,LMO-1
      Do 530 I=1,IM
      FLUX = DT2*MW(1,1,L)*(UO(I,1,L)+UO(I,1,L+1))
      dUM(I,1,L)   = dUM(I,1,L)   - FLUX
  530 dUM(I,1,L+1) = dUM(I,1,L+1) + FLUX
      GoTo 560
  540 Do 550 L=1,LMO-1
      Do 550 I=1,IM
      FLUX = DT2*MW(1,JM,L)*(UO(I,JM,L)+UO(I,JM,L+1))
      dUM(I,JM,L)   = dUM(I,JM,L)   - FLUX
  550 dUM(I,JM,L+1) = dUM(I,JM,L+1) + FLUX
  560 Continue
C**** V component
      Do 660 J=J1,JNP
      If (J==1)  GoTo 620
      If (J==JM-1)  GoTo 640
      Do 610 L=1,LMO-1
      Do 610 I=1,IM
      FLUX = DT4*(MW(I,J,L)+MW(I,J+1,L))*(VO(I,J,L)+VO(I,J,L+1))
      dVM(I,J,L)   = dVM(I,J,L)   - FLUX
  610 dVM(I,J,L+1) = dVM(I,J,L+1) + FLUX
      GoTo 660
  620 Do 630 L=1,LMO-1
      Do 630 I=1,IM
      FLUX = DT4*(MW(1,1,L)+MW(I,2,L))*(VO(I,1,L)+VO(I,1,L+1))
      dVM(I,1,L)   = dVM(I,1,L)   - FLUX
  630 dVM(I,1,L+1) = dVM(I,1,L+1) + FLUX
      GoTo 660
  640 Do 650 L=1,LMO-1
      Do 650 I=1,IM
      FLUX = DT4*(MW(I,JM-1,L)+MW(1,JM,L))*(VO(I,JM-1,L)+VO(I,JM-1,L+1))
      dVM(I,JM-1,L)   = dVM(I,JM-1,L)   - FLUX
  650 dVM(I,JM-1,L+1) = dVM(I,JM-1,L+1) + FLUX
  660 Continue
C****
C**** Add changes to momentum
C****
      Do 730 L=1,LMO
C**** Calculate dUM and dVM at poles, then update UM and VM
C     If (QSP)  Then
C       dUMSP =   Sum (dUM(:,1,L)*COSU(:)) * 2/IM
C       dVMSP = - Sum (dUM(:,1,L)*SINU(:)) * 2/IM
C       UM2(IVSP,1,L) = UM0(IVSP,1,L) + dVMSP
C       UM2(IM  ,1,L) = UM0(IM  ,1,L) + dUMSP  ;  EndIf
      If (QNP)  Then
        dUMNP = Sum (dUM(:,JM,L)*COSU(:)) * 2/IM
        dVMNP = Sum (dUM(:,JM,L)*SINU(:)) * 2/IM
        UM2(IM  ,JM,L) = UM0(IM  ,JM,L) + dUMNP
        UM2(IVNP,JM,L) = UM0(IVNP,JM,L) + dVMNP  ;  EndIf
C**** Update UM and VM away from poles
      Do 710 J=J1P,JNP
  710 UM2(:,J,L) = UM0(:,J,L) + dUM(:,J,L)
      Do 720 J=J1,JNP
  720 VM2(:,J,L) = VM0(:,J,L) + dVM(:,J,L)
  730 Continue
      Return
      EndSubroutine OADVV

      Subroutine OPGF (UM,VM,DT1)
C****
C**** OPGF adds the pressure gradient force to the momentum
C**** Input: G0M (J), GZM, S0M (kg), SZM, DT (s), MO (kg/m^2)
C**** Output: UM (kg*m/s) = UM - DT*(DH*D(P)+MO*D(PHI))*DYP
C****         VM (kg*m/s) = VM - DT*(DH*D(P)+MO*D(PHI))*DXV
C****
      Use CONSTANT, Only: GRAV
      Use OCEAN, Only: IM,JM,LMO, IVNP,IVSP,
     *                 LMOM=>LMM,LMOU=>LMU,LMOV=>LMV,
     *                 MO, G0M,GZM=>GZMO, S0M,SZM=>SZMO,
     *                 FOCEAN, mZSOLID=>HOCEAN, DXPGF,DYPGF,
     *                 COSI=>COSIC,SINI=>SINIC, OPRESS,OGEOZ,OPBOT
      Use OCEAN_DYN, Only: VBAR,dH, GUP,GDN, SUP,SDN

!      Use DOMAIN_DECOMP_1D, Only: GRID, HALO_UPDATE, NORTH
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Real*8,Intent(Out),
     *  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) :: UM,VM
      Real*8,Intent(In) :: DT1
      Real*8,External   :: VOLGSP
C****
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     *          dUM,dVM, P,ZG
      Real*8    PUP(LMO),PDN(LMO), PE,ZGE,VUP,VDN,dZGdP,
     *          dUMNP(LMO),dVMNP(LMO), DT2
      Integer*4 I,J,L,IMAX, J1,JN,J1P,JNP,JNH
      Logical*4 QSP,QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
C     QSP = J1==1           !    T      F      F
      QNP = JN==JM          !    F      F      T
C****
      Call HALO_UPDATE (GRID,OPRESS,FROM=NORTH)
C     Call HALO_UPDATE (GRID,  MO, FROM=NORTH)  !  copied in OFLUX
C     Call HALO_UPDATE (GRID, GUP, FROM=NORTH)  !  recalc in OPGF0
C     Call HALO_UPDATE (GRID, GDN, FROM=NORTH)  !  recalc in OPGF0
C     Call HALO_UPDATE (GRID, SUP, FROM=NORTH)  !  recalc in OPGF0
C     Call HALO_UPDATE (GRID, SDN, FROM=NORTH)  !  recalc in OPGF0
      DT2 = DT1/2
C****
C**** Calculate the mass weighted pressure P (Pa),
C**** geopotential ZG (m^2/s^2), and layer thickness dH (m)
C****
      Do 130 J=J1,JNH
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 130 I=1,IMAX
      If (FOCEAN(I,J) == 0)  GoTo 130
C**** Calculate pressure by integrating from the top down
      PE = OPRESS(I,J)
      Do 110 L=1,LMOM(I,J)
      P(I,J,L) = PE     + MO(I,J,L)*GRAV*.5
      PUP(L) = P(I,J,L) - MO(I,J,L)*GRAV*z12eH
      PDN(L) = P(I,J,L) + MO(I,J,L)*GRAV*z12eH
  110 PE     = PE       + MO(I,J,L)*GRAV
C**** save bottom pressure diagnostic
      OPBOT(I,J)=PE
C**** Calculate geopotential by integrating from the bottom up
      ZGE = - mZSOLID(I,J)*GRAV
      Do 120 L=LMOM(I,J),1,-1
      VUP = VOLGSP (GUP(I,J,L),SUP(I,J,L),PUP(L))
      VDN = VOLGSP (GDN(I,J,L),SDN(I,J,L),PDN(L))
      dZGdP = VUP*(.5-z12eH) + VDN*(.5+z12eH)
      VBAR(I,J,L) = (VUP + VDN)*.5
      ZG(I,J,L) = MO(I,J,L)*GRAV*.5*dZGdP + ZGE
      dH(I,J,L) = MO(I,J,L)*VBAR(I,J,L)
  120 ZGE = ZGE + dH(I,J,L)*GRAV
      OGEOZ(I,J) = ZGE
  130 Continue
C**** Copy VBAR and DH to all longitudes at north pole
      If (QNP) Then
        Do 140 L=1,LMOM(1,JM)
        VBAR(2:IM,JM,L) = VBAR(1,JM,L)
  140     DH(2:IM,JM,L) =   DH(1,JM,L)
      EndIf
C****
C**** Calculate smoothed East-West Pressure Gradient Force
C****
      Do 220 J=J1P,JNP
      Do 210 I=1,IM-1
      Do 210 L=1,LMOU(I,J)
  210 dUM(I,J,L) = DT2*DYPGF(J)*
     *             ((dH(I,J,L)+dH(I+1,J,L))*( P(I,J,L)- P(I+1,J,L)) +
     +              (MO(I,J,L)+MO(I+1,J,L))*(ZG(I,J,L)-ZG(I+1,J,L)))
      Do 220 L=1,LMOU(IM,J)
  220 dUM(IM,J,L) = DT2*DYPGF(J)*
     *              ((dH(IM,J,L)+dH(1,J,L))*( P(IM,J,L)- P(1,J,L)) +
     +               (MO(IM,J,L)+MO(1,J,L))*(ZG(IM,J,L)-ZG(1,J,L)))
      Call OPFIL (dUM)
C****
C**** Calculate North-South Pressure Gradient Force
C****
      Do 340 J=J1,JNP
      If (J==JM-1)  GoTo 320
      Do 310 I=1,IM
      Do 310 L=1,LMOV(I,J)
  310 dVM(I,J,L) = DT2*DXPGF(J)*
     *  ((dH(I,J,L)+dH(I,J+1,L))*( P(I,J,L)- P(I,J+1,L)) +
     +   (MO(I,J,L)+MO(I,J+1,L))*(ZG(I,J,L)-ZG(I,J+1,L)))
      GoTo 340
  320 Do 330 I=1,IM
      Do 330 L=1,LMOV(I,JM-1)
  330 dVM(I,JM-1,L) = DT2*DXPGF(JM-1)*
     *  ((dH(I,JM-1,L)+dH(1,JM,L))*( P(I,JM-1,L)- P(1,JM,L)) +
     +   (MO(I,JM-1,L)+MO(1,JM,L))*(ZG(I,JM-1,L)-ZG(1,JM,L)))
  340 Continue
C****
C**** Pressure Gradient Force at north pole
C****
      If (QNP)  Then
        dUMNP(:) = 0  ;  dVMNP(:) = 0
        Do 510 I=1,IM
        Do 510 L=1,LMOV(I,JM-1)
        dUMNP(L) = dUMNP(L) - SINI(I)*dVM(I,JM-1,L)
  510   dVMNP(L) = dVMNP(L) + COSI(I)*dVM(I,JM-1,L)
        Do 520 L=1,LMOM(1,JM)
        dUMNP(L) = dUMNP(L)*4*DXPGF(JM) / (IM*DXPGF(JM-1))
  520   dVMNP(L) = dVMNP(L)*4*DXPGF(JM) / (IM*DXPGF(JM-1))  ;  EndIf
C****
C**** Add pressure gradient force to momentum
C****
C**** Update UM away from poles
      Do 630 J=J1,JNP
      Do 620 I=1,IM
      Do 610 L=1,LMOU(I,J)
  610 UM(I,J,L) = UM(I,J,L) + dUM(I,J,L)
C**** Update VM away from poles
      Do 620 L=1,LMOV(I,J)
  620 VM(I,J,L) = VM(I,J,L) + dVM(I,J,L)
  630 Continue
C**** UM and VM at poles
      If (QNP)  Then
        Do 650 L=1,LMOM(1,JM)
        UM(IM  ,JM,L) = UM(IM  ,JM,L) + dUMNP(L)
  650   UM(IVNP,JM,L) = UM(IVNP,JM,L) + dVMNP(L)  ;  EndIf
      Return
      End Subroutine OPGF

      Subroutine OPGF0
C****
C**** OPGF0 calculates GUP, GDN, SUP and SDN inside each ocean cell,
C**** which will be kept constant during the dynamics of momentum
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM,
     *                 G0M,GZM=>GZMO, S0M,SZM=>SZMO
      Use OCEAN_DYN, Only: MMI, GUP,GDN, SUP,SDN
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Integer*4 I,J,L,IMAX, J1,JN,JNH
C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
C****
      Call HALO_UPDATE (GRID, MMI, FROM=NORTH)
      Call HALO_UPDATE (GRID, G0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, GZM, FROM=NORTH)
      Call HALO_UPDATE (GRID, S0M, FROM=NORTH)
      Call HALO_UPDATE (GRID, SZM, FROM=NORTH)
C****
      Do 10 J=J1,JNH
      IMAX=IM  ;  If(J==1.or.J==JM) IMAX=1
      Do 10 I=1,IMAX
      Do 10 L=1,LMOM(I,J)
      GUP(I,J,L) = (G0M(I,J,L) - 2*z12eH*GZM(I,J,L)) / MMI(I,J,L)
      GDN(I,J,L) = (G0M(I,J,L) + 2*z12eH*GZM(I,J,L)) / MMI(I,J,L)
      SUP(I,J,L) = (S0M(I,J,L) - 2*z12eH*SZM(I,J,L)) / MMI(I,J,L)
   10 SDN(I,J,L) = (S0M(I,J,L) + 2*z12eH*SZM(I,J,L)) / MMI(I,J,L)
      Return
      End Subroutine OPGF0

      SUBROUTINE OADVT (RM,RX,RY,RZ,DT,QLIMIT, OIJL)
!@sum  OADVT advects tracers using the linear upstream scheme.
!@auth Gary Russell
C****
C**** Input:  MB (kg) = mass before advection
C****          DT (s) = time step
C****       MU (kg/s) = west to east mass flux
C****       MV (kg/s) = south to north mass flux
C****       MW (kg/s) = downward vertical mass flux
C****          QLIMIT = whether slope limitations should be used
C**** Output: RM (kg) = tracer mass
C****   RX,RY,RZ (kg) = first moments of tracer mass
C****       OIJL (kg) = diagnostic accumulation of tracer mass flux
C****
      USE OCEAN, only : im,jm,lmo
      USE OCEAN_DYN, only : mb=>mmi,smu,smv,smw

      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),     DIMENSION
     *     (IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: RM,RX,RY,RZ
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO,3) :: OIJL
      REAL*8,
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MA
      INTEGER I,J,L, J_0H
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT

      logical :: HAVE_NORTH_POLE  ! ,HAVE_SOUTH_POLE
         call getDomainBounds(grid, J_STRT_HALO=J_0H)
         call getDomainBounds(grid, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C        call getDomainBounds(grid, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

C****
C**** Load mass after advection from mass before advection
C****
      MA = MB
C****
C**** Advect the tracer using the Slopes Scheme
C****
      CALL OADVTX (RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))
      CALL OADVTY (RM,RX,RY,RZ,MA,SMV,     DT,QLIMIT,OIJL(1,J_0H,1,2))
      CALL OADVTZ (RM,RX,RY,RZ,MA,SMW,     DT,QLIMIT,OIJL(1,J_0H,1,3))
      CALL OADVTX (RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))
C**** Fill in values at the poles
      if (HAVE_NORTH_POLE) then
      DO 20 L=1,LMO
      DO 20 I=1,IM
      RM(I,JM,L) = RM(1,JM,L)
   20 RZ(I,JM,L) = RZ(1,JM,L)
      end if
C     if (HAVE_SOUTH_POLE) then
C     DO 30 L=1,LMO
C     DO 30 I=1,IM
C     RM(I, 1,L) = RM(IM,1,L)
C  30 RZ(I, 1,L) = RZ(IM,1,L)
C     end if

      RETURN
      END SUBROUTINE OADVT

      Subroutine OADVTX (RM,RX,RY,RZ,MM,MU,DT,QLIMIT,OIJL)
C****
!@sum   OADVTX advects tracer in X direction via Linear Upstream Scheme
!@auth  Gary L. Russell
!@ver   2009/01/14
C****
C**** If QLIMIT is true, gradients are limited to prevent mean tracer
C**** from becoming negative; Abs[A(I)] must not exceed 1.
C****
C**** Input: DT (s) = time step
C****     MU (kg/s) = west to east mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                   MM (kg) = ocean mass
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, zDXYP=>BYDXYPO, FOCEAN
      Use OCEANR_DIM, Only: oGRID
      Implicit None
C**** Interface variables
      Logical*4,Intent(In) :: QLIMIT
      Real*8,   Intent(In) :: DT
      Real*8,Intent(InOut),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ, MM,MU, OIJL
C**** Local variables
      Integer*4 I,J,L, Im1,Ip1, J1P,JNP, N,NCOURANT
      Real*8    AM(IM),A(IM),FM(IM),FX(IM),FY(IM),FZ(IM),
     *          MMnew, zCOURANT, RXY
C**** Extract domain decomposition band parameters
      J1P = Max (oGRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (oGRID%J_STOP, JM-1)  !  Exclude north pole
C****
C**** Loop over layers and latitudes
C****
      Do 320 L=1,LMO
      Do 320 J=J1P,JNP
C****
C**** Determine number of Courant iterations NCOURANT
C****
      NCOURANT = 1
      Im1=IM-1
      I = IM
      Do 40 Ip1=1,IM
      If (MU(I,J,L)) 10,30,20
C**** MU < 0, insist that  -AM(I)/N < MM(Ip1)
C****                 and  -AM(I)/N < MM(Ip1) + [AM(I)-AM(Ip1)]*(N-1)/N
   10 MMnew = MM(Ip1,J,L) + DT*(MU(I,J,L)-MU(Ip1,J,L))
      If (MMnew <= 0)  GoTo 801
      If (-DT*MU(I,J,L) > NCOURANT*MM(Ip1,J,L))  Then
        NCOURANT = Ceiling (-DT*MU(I,J,L)/MM(Ip1,J,L))
        Write (6,901) NCOURANT,I,J,L, DT*MU(I,J,L)*zDXYP(J),
     *                MM(Ip1,J,L)*zDXYP(J),MMnew*zDXYP(J)  ;  EndIf
      If (-DT*MU(Ip1,J,L) > NCOURANT*MMnew)  Then
        NCOURANT = Ceiling (-DT*MU(Ip1,J,L)/MMnew)
        Write (6,901) NCOURANT,I,J,L, DT*MU(I,J,L)*zDXYP(J),
     *                MM(Ip1,J,L)*zDXYP(J),MMnew*zDXYP(J)  ;  EndIf
      GoTo 30
C**** MU > 0, insist that  AM(I)/N < MM(I)
C****                 and  AM(I)/N < MM(I) + [AM(Im1)-AM(I)]*(N-1)/N
   20 MMnew = MM(I,J,L) + DT*(MU(Im1,J,L)-MU(I,J,L))
      If (MMnew <= 0)  GoTo 802
      If (DT*MU(I,J,L) > NCOURANT*MM(I,J,L))  Then
        NCOURANT = Ceiling (DT*MU(I,J,L)/MM(I,J,L))
        Write (6,901) NCOURANT,I,J,L, DT*MU(I,J,L)*zDXYP(J),
     *                MM(I,J,L)*zDXYP(J),MMnew*zDXYP(J)  ;  EndIf
      If (DT*MU(Im1,J,L) > NCOURANT*MMnew)  Then
        NCOURANT = Ceiling (DT*MU(Im1,J,L)/MMnew)
        Write (6,901) NCOURANT,I,J,L, DT*MU(I,J,L)*zDXYP(J),
     *                MM(I,J,L)*zDXYP(J),MMnew*zDXYP(J)  ;  EndIf
   30 Im1=I
      I=Ip1
   40 Continue
C**** Loop over Courant iterations
      zCOURANT = 1d0 / NCOURANT
      Do 320 N=1,NCOURANT
C****
C**** Calculate FM (kg), FX (kg**2), FY (kg) and FZ (kg)
C****
      I=IM
      Do 140 Ip1=1,IM
      AM(I) = DT*MU(I,J,L)*zCOURANT
      If (AM(I)) 110,120,130
C**** Ocean mass flux is negative
  110 A(I)  = AM(I) / MM(Ip1,J,L)
      FM(I) = A(I) * (RM(Ip1,J,L) - (1+A(I))*RX(Ip1,J,L))
      FX(I) = AM(I) * (A(I)*A(I)*RX(Ip1,J,L) - 3*FM(I))
      FY(I) = A(I)*RY(Ip1,J,L)
      FZ(I) = A(I)*RZ(Ip1,J,L)
      GoTo 140
C**** Ocean mass flux is zero
  120 A(I)  = 0
      FM(I) = 0
      FX(I) = 0
      FY(I) = 0
      FZ(I) = 0
      GoTo 140
C**** Ocean mass flux is positive
  130 A(I)  = AM(I) / MM(I,J,L)
      FM(I) = A(I) * (RM(I,J,L) + (1-A(I))*RX(I,J,L))
      FX(I) = AM(I) * (A(I)*A(I)*RX(I,J,L) - 3*FM(I))
      FY(I) = A(I)*RY(I,J,L)
      FZ(I) = A(I)*RZ(I,J,L)
  140 I=Ip1
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      If (.not.QLIMIT)  GoTo 300
      Im1=IM
      Do 290 I=1,IM
      If (A(Im1) >= 0)  GoTo 240
C**** Water is leaving through left edge: 2 or 3 divisions
      If (FM(Im1) <= 0)  GoTo 210
C**** Left most division is negative, RML = -FM(I-1) < 0: Case 2 or 4
      RX(I,J,L) = RM(I,J,L) / (1+A(Im1))
      FM(Im1) = 0
      FX(Im1) = AM(Im1)*A(Im1)*A(Im1)*RX(I,J,L)
      If (A(I) <= 0)  GoTo 290
      FM(I) = A(I) * (RM(I,J,L) + (1-A(I))*RX(I,J,L))
      FX(I) = AM(I) * (A(I)*A(I)*RX(I,J,L) - 3*FM(I))
      GoTo 290
C**** Left most division is non-negative, RML = -FM(I-() > 0:
C**** Case 1, 3 or 5
  210 If (A(I) <= 0)  GoTo 230
C**** Water is leaving through right edge: 3 divisions
      If (FM(I) >= 0)  GoTo 290
C**** Right most division is negative, RMR = FM(I) < 0: Case 3 or 5
  220 RX(I,J,L) = -RM(I,J,L) / (1-A(I))
      FM(I) = 0
      FX(I) = AM(I)*A(I)*A(I)*RX(I,J,L)
      FM(Im1) = A(Im1) * (RM(I,J,L) - (1+A(Im1))*RX(I,J,L))
      FX(Im1) = AM(Im1) * (A(Im1)*A(Im1)*RX(I,J,L) - 3*FM(Im1))
      GoTo 290
C**** No water is leaving through right edge: 2 divisions
  230 If (RM(I,J,L)+FM(Im1) >= 0)  GoTo 290
C**** Right most division is negative, RMR = RM(I,J)+FM(I-1) < 0: Case 3
      RX(I,J,L) = RM(I,J,L) / A(Im1)
      FM(Im1) = -RM(I,J,L)
      FX(Im1) = AM(Im1)*(A(Im1)+3)*RM(I,J,L)
      GoTo 290
C**** No water is leaving through the left edge: 1 or 2 divisions
  240 If (A(I) <= 0)  GoTo 290
C**** Water is leaving through right edge: 2 divisions
      If (FM(I) >= 0)  GoTo 250
C**** Right most division is negative, RMR = FM(I) < 0: Case 3
      RX(I,J,L) = -RM(I,J,L) / (1-A(I))
      FM(I) = 0
      FX(I) = AM(I)*A(I)*A(I)*RX(I,J,L)
      GoTo 290
C**** Right most division is non-negative, RMR = FM(I) > 0: Case 1 or 2
  250 If (RM(I,J,L)-FM(I) >= 0)  GoTo 290
C**** Left most division is negative, RML = RM(I,J)-FM(I) < 0: Case 2
      RX(I,J,L) = RM(I,J,L) / A(I)
      FM(I) = RM(I,J,L)
      FX(I) = AM(I)*(A(I)-3)*RM(I,J,L)
C****
  290 Im1=I
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
  300 Im1=IM
      Do 310 I=1,IM
      If (L > LMOM(I,J))  GoTo 310
      RM(I,J,L) = RM(I,J,L) +  (FM(Im1)-FM(I))
      RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FX(Im1)-FX(I)) +
     +    3*((AM(Im1)+AM(I))*RM(I,J,L)-MM(I,J,L)*(FM(Im1)+FM(I))))
     /  / (MM(I,J,L)+AM(Im1)-AM(I))
      RY(I,J,L) = RY(I,J,L) + (FY(Im1)-FY(I))
      RZ(I,J,L) = RZ(I,J,L) + (FZ(Im1)-FZ(I))
      MM(I,J,L) = MM(I,J,L) + (AM(Im1)-AM(I))
        OIJL(I,J,L) = OIJL(I,J,L) + FM(I)
C**** Limit tracer gradients if necessary
      If (QLIMIT)  Then
        If (RM(I,J,L) < 0)  GoTo 831
        RXY = Abs(RX(I,J,L)) + Abs(RY(I,J,L))
        If (RXY > RM(I,J,L))  Then
          RX(I,J,L) = RX(I,J,L)*(RM(I,J,L) / (RXY+Tiny(RXY)))
          RY(I,J,L) = RY(I,J,L)*(RM(I,J,L) / (RXY+Tiny(RXY)))  ;  EndIf
        If (Abs(RZ(I,J,L)) > RM(I,J,L))
     *    RZ(I,J,L) = Sign (RM(I,J,L),RZ(I,J,L))  ;  EndIf
C****
  310 Im1=I     !  End of Do loop over I
  320 Continue  !  End of Do loops over N,J,L
      Return
C****
  801 Write (6,*) 'MMnew < 0 in OADVTX'
      Write (6,*) 'I,J,L,MOold,MOnew,AM(Im1),AM(I) (kg/m^2)=',Ip1,J,L,
     *            MM(Ip1,J,L)*zDXYP(J), MMnew*zDXYP(J),
     *            DT*MU(I,J,L)*zDXYP(J), DT*MU(Ip1,J,L)*zDXYP(J)
      Call STOP_MODEL ('OADVTX',255)
  802 Write (6,*) 'MMnew < 0 in OADVTX'
      Write (6,*) 'I,J,L,MOold,MOnew,AM(Im1),AM(I) (kg/m^2)=',I,J,L,
     *            MM(I,J,L)*zDXYP(J), MMnew*zDXYP(J),
     *            DT*MU(Im1,J,L)*zDXYP(J), DT*MU(I,J,L)*zDXYP(J)
      Call STOP_MODEL ('OADVTX',255)
  831 Write (6,*) 'RM < 0 in OADVTX:',I,J,L,RM(I,J,L)
      Call STOP_MODEL ('OADVTX',255)
C****
  901 Format (/ 'NCOURANT,I,J,L,AM,MMold,MMnew(kg/m^2)=',4I5,3F12.2)
      EndSubroutine OADVTX

      subroutine oadvty (rm,rx,ry,rz,mo,mv,dt,qlimit,oijl)
!@sum  OADVTY advection driver for y-direction
!@auth Gary Russell, modified by T.Clune, R. Ruedy
c****
c**** oadvty advects tracers in the south to north direction using the
c**** linear upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg/s) = north-south mo flux, positive northward
c****        qlimit = whether slope limitations should be used
c****
c**** input/output:
c****     rm     (kg) = tracer mass
c****   rx,ry,rz (kg) = 1st moments of tracer mass
c****     mo     (kg) = ocean mass
c****
!      use DOMAIN_DECOMP_1D, only : grid, get, halo_update
      use DOMAIN_DECOMP_1D, only : getDomainBounds, halo_update
      USE OCEANR_DIM, only : grid=>ogrid

      use DOMAIN_DECOMP_1D, only : halo_update_column
      use DOMAIN_DECOMP_1D, only : NORTH, SOUTH, AM_I_ROOT
      use OCEAN, only : im,jm,lmo,lmm,focean
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &                  rm, rx,ry,rz, mo,mv, oijl
      real*8, intent(in) :: dt
      logical, intent(in) ::  qlimit

      REAL*8, dimension(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &        BM, fm,fx,fy,fz
      integer :: i,j,L,ierr,ICKERR, err_loc(3)
      REAL*8, DIMENSION(LMO) ::
     &     m_np,rm_np,rzm_np ! ,m_sp,rm_sp,rzm_sp

c****Get relevant local distributed parameters
      INTEGER J_0,J_1,J_0H,J_1H,J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0,
     &               J_STOP = J_1, J_STOP_SKP=J_1S,
     &               J_STRT_HALO = J_0H,
     &               J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** loop over layers
      ICKERR=0
      do L=1,lmo

c****   fill in and save polar values
c       if (HAVE_SOUTH_POLE) then
c         m_sp(L) = mo(1,1,L)
c         mo(2:im,1,L) = mo(1,1,L)
c         rm_sp(L) = rm(1,1,L)
c         rm(2:im,1,L) = rm(1,1,L)
c         rzm_sp(L)  = rz(1,1,L)
c         rz(2:im,1,L)  = rz(1,1,L)
c       end if                       !SOUTH POLE

        if (HAVE_NORTH_POLE) then
          m_np(L) = mo(1,jm,L)
          mo(2:im,jm,L) = mo(1,jm,L)
          rm_np(L) = rm(1,jm,l)
          rm(2:im,jm,L) = rm(1,jm,L)
          rzm_np(L)  = rz(1,jm,l)
          rz(2:im,jm,L)  = rz(1,jm,L)
        end if                       !NORTH POLE
c****
c****   convert flux to water mass moving north (if >0)
        do j=j_0,j_1
        do i=1,im
          if(focean(i,j).gt.0. .and. L.le.lmm(i,j)) then
            bm(i,j,L)=dt*mv(i,j,L)
          else
            bm(i,j,L)=0.
            mo(i,j,L)=0.
          end if
        end do
        end do

c****   POLES: set horiz. moments to zero
        IF (HAVE_SOUTH_POLE) THEN
          rx(:,1,L) = 0. ; ry(:,1,L) = 0.
        end if
        IF (HAVE_NORTH_POLE) THEN
          bm(:,jm,L) = 0.
          rx(:,jm,L) = 0. ; ry(:,jm,L) = 0.
        END IF
      end do   ! loop over layers

c****
c**** call 1-d advection routine
c****
        call advec_lin_1D_custom(
     &     rm(1,j_0h,1), rx(1,j_0h,1),ry(1,j_0h,1),rz(1,j_0h,1),
     &     fm(1,j_0h,1), fx(1,j_0h,1),fy(1,j_0h,1),fz(1,j_0h,1),
     &     mo(1,j_0h,1), bm(1,j_0h,1),
     &     qlimit,ierr,err_loc)

        if (ierr.gt.0) then
          write(6,*) "Error in oadvty: i,j,l=",err_loc
          if (ierr == 2) then
ccc         write(0,*) "Error in qlimit: abs(b) > 1"
ccc         call stop_model('Error in qlimit: abs(b) > 1',11)
            ICKERR=ICKERR+1
          endif
        end if
! horizontal moments are zero at pole
c       IF (HAVE_SOUTH_POLE) then
c          rx(:,1, :) = 0  ; ry(:,1, :) = 0
c       end if
        IF (HAVE_NORTH_POLE) then
           rx(:,jm,:) = 0. ; ry(:,jm,:) = 0.
        end if

      do l=1,lmo

c****   average and update polar boxes
c       if (HAVE_SOUTH_POLE) then
c         mo(:,1 ,l) = (m_sp(l) + sum(mo(:,1 ,l)-m_sp(l)))/im
c         rm(:,1 ,l) = (rm_sp(l) + sum(rm(:,1 ,l)-rm_sp(l)))/im
c         rz(:,1 ,l) = (rzm_sp(l) + sum(rz(:,1 ,l)-rzm_sp(l) ))/im
c       end if   !SOUTH POLE

        if (HAVE_NORTH_POLE .and. L.le.lmm(1,jm)) then
          mo(1,jm,l) = m_np(l)   + sum(bm(:,jm-1,l))/im
          rm(1,jm,l) = rm_np(l)  + sum(fm(:,jm-1,l))/im
          rz(1,jm,l) = rzm_np(l) + sum(fz(:,jm-1,l))/im
          if(mo(1,jm,l)<0. .or. (qlimit.and.rm(1,jm,l)<0.)) then
            ICKERR=ICKERR+1
            write(0,*) 'oadvty: mo or salt<0 at North Pole layer',
     *      L,mo(1,jm,L),rm(1,jm,L)
          endif
          mo(2:im,jm,l)=mo(1,jm,l)
          rm(2:im,jm,l)=rm(1,jm,l)
          rz(2:im,jm,l)=rz(1,jm,l)
        end if  !NORTH POLE

      enddo ! end loop over levels
c
c**** sum into oijl
c
      do l=1,lmo ; do j=J_0,J_1S
          oijl(:,j,l)  = oijl(:,j,l) + fm(:,j,l)
      enddo ; enddo
C
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in oadvty',11)
C
      return
c****
      end subroutine oadvty

      subroutine advec_lin_1D_custom(s,sx,sy,sz, f,fx,fy,fz, mass,dm,
     *     qlimit,ierr, err_loc)
!@sum  advec_lin_1d_custom is a parallel variant of adv1d but for a
!@sum  linear upstream scheme.
!@auth T. Clune, R. Ruedy
c--------------------------------------------------------------
c adv1d advects tracers in j-direction using the lin ups scheme
c--------------------------------------------------------------
      use ocean, only: im,jm,lmo,lmm,focean

!      USE DOMAIN_DECOMP_1D, only: grid, GET
      USE DOMAIN_DECOMP_1D, only: GETDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      USE DOMAIN_DECOMP_1D, only: NORTH, SOUTH
      USE DOMAIN_DECOMP_1D, only: HALO_UPDATE, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only: CHECKSUM
      implicit none
      !
!@var s mean tracer amount (kg or J)
!@var sx,sy,sz lin tracer moments (kg or J)
!@var f tracer flux (diagnostic output) (kg or J)
!@var fx,fy,fz tracer moment flux (diagnostic output) (kg or J)
!@var mass mass field (kg)
!@var dm mass flux (kg)
!@var qlimit true if negative tracer is to be avoided
!@var ierr, nerr error codes
      logical, intent(in) :: qlimit
      REAL*8, dimension(im, grid%j_strt_halo:grid%j_stop_halo,lmo)
     *  :: s,sx,sy,sz, f,fx,fy,fz, mass,dm
      integer :: n,np1,nm1,nn,ns
      integer,intent(out) :: ierr,err_loc(3) ! 3 dimensions
      REAL*8 :: fracm,frac1,mnew
      INTEGER :: i
      INTEGER :: l
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE

      ierr=0

      call getDomainBounds(grid, J_STRT = J_0, J_STOP=J_1,
     & J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     & HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

      CALL HALO_UPDATE(grid, mass, FROM=NORTH)
      CALL HALO_UPDATE(grid, s, FROM=NORTH)
      CALL HALO_UPDATE(grid, sx, FROM=NORTH)
      CALL HALO_UPDATE(grid, sy, FROM=NORTH)
      CALL HALO_UPDATE(grid, sz, FROM=NORTH)

      DO l=1,lmo
        Do i=1, im
          Call calc_tracer_mass_flux() ! f from s,sy,mass,dm
        End Do                         ! temp. store dm/mass in fy
      end do

      CALL HALO_UPDATE(grid, fy, FROM=NORTH+SOUTH) ! still holding dm/m
      ! Limit fluxes to maintain positive mean values?
      If (qlimit) Then

        Call HALO_UPDATE(grid, f, FROM=NORTH+SOUTH)
        DO l=1,lmo
          Do i=1, im
            Call apply_limiter() ! adjusting f,sy
          End Do
        end do
        Call HALO_UPDATE(grid, f, FROM=NORTH+SOUTH)
      End If
c--------------------------------------------------------------------
         ! calculate tracer fluxes of slopes fx,fy,fz
c--------------------------------------------------------------------
      DO l=1,lmo
        Do i=1, im
          Call tracer_slopes()
        end do
      end do
c-------------------------------------------------------------------
c update tracer mass, moments of tracer mass, air mass distribution
c-------------------------------------------------------------------
      CALL HALO_UPDATE(grid, f,  FROM=SOUTH)
      CALL HALO_UPDATE(grid, dm, FROM=SOUTH)
      CALL HALO_UPDATE(grid, fx, FROM=SOUTH)
      CALL HALO_UPDATE(grid, fy, FROM=SOUTH)
      CALL HALO_UPDATE(grid, fz, FROM=SOUTH)

      DO l=1,lmo
        DO i = 1, im
          Call update_tracer_mass() ! s also: sx,sy,sz, mass
        end do
      enddo

      return

      Contains

      Integer Function NeighborByFlux(n, dm)
        Integer, Intent(In) :: n
        Real*8, Intent(In) :: dm

        Integer :: nn

        If (dm < 0) Then ! air mass flux is negative
          nn=n+1
        else ! air mass flux is positive
          nn=n
        endif
        NeighborByFlux = nn
      End Function NeighborByFlux

      Function FluxFraction(dm) Result(frac)
        Real*8, Intent(In) :: dm
        Real*8 :: frac

        If (dm < 0 ) Then
          frac = +1.
        Else ! Flux non negative
          frac = -1.
        End If
      End Function FluxFraction

      Function MassFraction(dm, mass) Result (fracm)
        Real*8, Intent(In) :: dm
        Real*8, Intent(In) :: mass
        Real*8 :: fracm

        If (mass > 0.0d0) Then
          fracm = dm / mass
        Else
          fracm = 0.d0
        End If
      End Function MassFraction

      Subroutine calc_tracer_mass_flux()

      Do n = J_0, J_1

        nn = NeighborByFlux(n, dm(i,n,l))
        fracm = MassFraction(dm(i,n,l), mass(i,nn,l))

        frac1 = fracm + FluxFraction(dm(i,n,l))

        f(i,n,l)=fracm*(s(i,nn,l)-frac1*sy(i,nn,l))
      ! temporary storage of fracm in fy, to be used below
        fy(i,n,l)=fracm
      !
      enddo
      End Subroutine calc_tracer_mass_flux

      Subroutine tracer_slopes()

      Do n = J_0, J_1

        nn = NeighborByFlux(n, dm(i,n,l))
      ! retrieving fracm, which was stored in fy
        fracm=fy(i,n,l)
      !
        fy(i,n,l)=
     &     dm(i,n,l)*(fracm*fracm*sy(i,nn,l)-3.*f(i,n,l))
      ! cross moments
        fx(i,n,l)=fracm*sx(i,nn,l)
        fz(i,n,l)=fracm*sz(i,nn,l)
      enddo
      End Subroutine tracer_slopes

      Subroutine apply_limiter
      REAL*8 :: an, anm1, fn, fnm1, sn, syn

c     If (HAVE_SOUTH_POLE) Then
c        n = J_0
c        an = fy(i,n,l) ! reading fracm which was stored in fx
c        anm1 = 0
c        fn = f(i,n,l)
c        fnm1 = 0
c        sn = s(i,n,l)
c        syn = sy(i,n,l)
c        call limitq_lin(anm1,an,fnm1,fn,sn,syn)
c        f(i,n,l) = fn
c        sy(i,n,l) = syn
c     End If

      DO n = J_0S, J_1S+1

         an = fy(i,n,l) ! reading fracm which was stored in fy
         anm1 = fy(i,n-1,l)
         fn = f(i,n,l)
         fnm1 = f(i,n-1,l)
         sn = s(i,n,l)
         syn = sy(i,n,l)
         call limitq_lin(anm1,an,fnm1,fn,sn,syn)
         f(i,n,l) = fn
         f(i,n-1,l) = fnm1
         sy(i,n,l) = syn

      enddo
      End Subroutine apply_limiter

      Subroutine update_tracer_mass() ! non-polar only
      REAL*8 :: sxy

      Do n=J_0S,J_1S

         if(focean(i,n).le.0. .or. l.gt.lmm(i,n)) cycle
         mnew=mass(i,n,l) + dm(i,n-1,l)-dm(i,n,l)

         s(i,n,l)=s(i,n,l) + (f(i,n-1,l)-f(i,n,l))

         sy(i,n,l)=( sy(i,n,l)*mass(i,n,l) + (fy(i,n-1,l)-fy(i,n,l))
     &              + 3.*( (dm(i,n-1,l)+dm(i,n,l))*s(i,n,l)
     &                     -mass(i,n,l)*(f(i,n-1,l)+f(i,n,l)) ) )/mnew
      ! cross moments
         sx(i,n,l) = sx(i,n,l) + (fx(i,n-1,l)-fx(i,n,l))
         sz(i,n,l) = sz(i,n,l) + (fz(i,n-1,l)-fz(i,n,l))
      !
         if (qlimit) then ! limit tracer gradients
           sxy = abs(sx(i,n,l)) + abs(sy(i,n,l))
           if ( sxy > s(i,n,l) ) then
             sx(i,n,l) = sx(i,n,l)*( s(i,n,l)/(sxy + tiny(sxy)) )
             sy(i,n,l) = sy(i,n,l)*( s(i,n,l)/(sxy + tiny(sxy)) )
           end if
           if ( abs(sz(i,n,l)) > s(i,n,l) )
     *       sz(i,n,l) = sign(s(i,n,l),sz(i,n,l)+0.d0)
         end if
c------------------------------------------------------------------
         mass(i,n,l) = mnew
         if(mass(i,n,l).le.0. .or. (qlimit.and.s(i,n,l).lt.0.)) then
            ierr=2
            err_loc=(/ i, n, l /)
            write(0,*) 'oadvty: mo or salt<0 at',err_loc,mnew,s(i,n,l)
            return
         endif
c-----------------------------------------------------------------

      enddo
      End Subroutine Update_Tracer_Mass

      end subroutine advec_lin_1D_custom

      subroutine limitq_lin(anm1,an,fnm1,fn,sn,sx)
!@sum  limitq adjusts moments to maintain non-neg. tracer means/fluxes
!@auth G. Russell, modified by Maxwell Kelley
        implicit none
        REAL*8 :: anm1,an,fnm1,fn,sn,sx
c local variables
        REAL*8 :: sl,sc,sr, frl,frl1, frr,frr1, gamma,g13ab,
     &       fr,fr1, fsign,su,sd
c****
c**** modify the tracer moments so that the tracer mass in each
c**** division is non-negative
c****
c**** no water leaving the box
        if(anm1.ge.0. .and. an.le.0.) return
c**** water is leaving through both the left and right edges
        if(anm1.lt.0. .and. an.gt.0.) then
           sl = -fnm1
           sr = +fn
c**** all divisions are non-negative
           if(sl.ge.0. .and. sr.ge.0.) return
c**** at least one division is negative
           frl = anm1
           frl1 = frl+1.
           frr = an
           frr1 = frr-1.
           if(sl.lt.0.) then           ! leftmost division
              sx = sn/frl1
              sr = frr*(sn-frr1*sx)
              sl = 0.
           else                        ! rightmost division
              sx = sn/frr1
              sl = -frl*(sn-frl1*sx)
              sr = 0.
           endif
           fnm1 = -sl
           fn   = +sr
        else
c**** water is leaving only through one edge
           if(an.gt.0.)  then ! right edge
              fr=an
              sd=fn
              fsign=-1.
           else                  ! left edge
              fr=anm1
              sd=-fnm1
              fsign=1.
           endif
           su = sn-sd
           if(sd.ge.0. .and. su.ge.0.) return
           fr1=fr+fsign
           if(sd.lt.0.)  then
c**** downstream division is negative
              sx = sn/fr1
              su = sn
           else
c**** upstream division is negative
              sx = sn/fr
              su = 0.
           endif
           sd = sn - su
           if(an.gt.0.) then
              fn=sd
           else
              fnm1=-sd
           endif
        endif
        return
      end subroutine limitq_lin

      SUBROUTINE OADVTZ (RM,RX,RY,RZ,MO,MW,DT,QLIMIT,OIJL)
C****
C**** OADVTZ advects tracers in the vertical direction using the
C**** linear upstream scheme.  If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C**** Input: DT (s) = time step
C****     MW (kg/s) = downward vertical mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm,focean

!      use domain_decomp_1d, only : grid, get
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ, OIJL, MO
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO-1) :: MW
      LOGICAL*4, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT
      REAL*8, DIMENSION(0:LMO) :: CM,C,FM,FX,FY,FZ
      INTEGER I,J,L,LMIJ,ICKERR,IMIN,IMAX
      REAL*8 SBMN,SFMN,SFZN,RXY

      INTEGER :: J_0,J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

C****
C**** Loop over latitudes and longitudes
      ICKERR=0
      DO J=J_0,J_1
        IMIN=1
        IMAX=IM
        IF (J == 1) IMIN=IM
        IF (J == JM) IMAX=1
        DO I=IMIN,IMAX
      CM(0) = 0.
       C(0) = 0.
      FM(0) = 0.
      FX(0) = 0.
      FY(0) = 0.
      FZ(0) = 0.
      LMIJ=LMM(I,J)
      IF(LMIJ.LE.1)  GO TO 330
      CM(LMIJ) = 0.
       C(LMIJ) = 0.
      FM(LMIJ) = 0.
      FX(LMIJ) = 0.
      FY(LMIJ) = 0.
      FZ(LMIJ) = 0.
C****
C**** Calculate FM (kg), FX (kg), FY (kg) and FZ (kg**2)
C****
      DO 120 L=1,LMIJ-1
      CM(L) = DT*MW(I,J,L)
      IF(CM(L).LT.0.)  GO TO 110
C**** Ocean mass flux is positive
      C(L)  = CM(L)/MO(I,J,L)
      IF(C(L).GT.1d0)  WRITE (6,*) 'C>1:',I,L,C(L),MO(I,J,L)
      FM(L) = C(L)*(RM(I,J,L)+(1d0-C(L))*RZ(I,J,L))
      FX(L) = C(L)*RX(I,J,L)
      FY(L) = C(L)*RY(I,J,L)
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,J,L)-3d0*FM(L))
      GO TO 120
C**** Ocean mass flux is negative
  110 C(L)  = CM(L)/MO(I,J,L+1)
      IF(C(L).LT.-1d0)  WRITE (6,*) 'C<-1:',I,L,C(L),MO(I,J,L+1)
      FM(L) = C(L)*(RM(I,J,L+1)-(1d0+C(L))*RZ(I,J,L+1))
      FX(L) = C(L)*RX(I,J,L+1)
      FY(L) = C(L)*RY(I,J,L+1)
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,J,L+1)-3d0*FM(L))
  120 CONTINUE
C****
C**** Modify the tracer moments so that the tracer mass in each
C**** division is non-negative
C****
      IF(.NOT.QLIMIT)  GO TO 300
      DO 290 L=1,LMIJ
      IF(C(L-1).GE.0.)  GO TO 240
C**** Water is leaving through the bottom edge: 2 or 3 divisions
      IF(FM(L-1).LE.0.)  GO TO 210
C**** Bottom most division is negative, RMB = -FM(L-1) < 0: Case 2 or 4
      RZ(I,J,L) = RM(I,J,L)/(1d0+C(L-1))
      FM(L-1) = 0.
      FZ(L-1) = CM(L-1)*C(L-1)*C(L-1)*RZ(I,J,L)
      IF(C(L).LE.0.)  GO TO 290
      FM(L) = C(L)*(RM(I,J,L)+(1d0-C(L))*RZ(I,J,L))
      FZ(L) = CM(L)*(C(L)*C(L)*RZ(I,J,L)-3d0*FM(L))
      GO TO 290
C**** Bottom most division is non-negative, RMB = -FM(L-1) > 0:
C**** Case 1, 3 or 5
  210 IF(C(L).LE.0.)  GO TO 230
C**** Water is leaving through the top edge: 3 divisions
      IF(FM(L).GE.0.)  GO TO 290
C**** Top most division is negative, RMT = FM(L) < 0: Case 3 or 5
      RZ(I,J,L) = -RM(I,J,L)/(1d0-C(L))
      FM(L) = 0.
      FZ(L) = CM(L)*C(L)*C(L)*RZ(I,J,L)
      FM(L-1) = C(L-1)*(RM(I,J,L)-(1d0+C(L-1))*RZ(I,J,L))
      FZ(L-1) = CM(L-1)*(C(L-1)*C(L-1)*RZ(I,J,L)-3d0*FM(L-1))
      GO TO 290
C**** No water is leaving through the top edge: 2 divisions
  230 IF(RM(I,J,L)+FM(L-1).GE.0.)  GO TO 290
C**** Top most division is negative, RMT = RM(I,J,L)+FM(L-1) < 0: Case 3
      RZ(I,J,L) = RM(I,J,L)/C(L-1)
      FM(L-1) = -RM(I,J,L)
      FZ(L-1) = CM(L-1)*(C(L-1)+3d0)*RM(I,J,L)
      GO TO 290
C**** No water is leaving through the bottom edge: 1 or 2 divisions
  240 IF(C(L).LE.0.)  GO TO 290
C**** Water is leaving through the top edge: 2 divisions
      IF(FM(L).GE.0.)  GO TO 250
C**** Top most division is negative, RMT = FM(L) < 0: Case 3
      RZ(I,J,L) = -RM(I,J,L)/(1d0-C(L))
      FM(L) = 0.
      FZ(L) = CM(L)*C(L)*C(L)*RZ(I,J,L)
      GO TO 290
C**** Top most division is non-negative, RMT = FM(L) > 0: Case 1 or 2
  250 IF(RM(I,J,L)-FM(L).GE.0.)  GO TO 290
C**** Bottom most division is negative, RMB = RM(I,J,L)-FM(L) < 0: Cas 2
      RZ(I,J,L) = RM(I,J,L)/C(L)
      FM(L) = RM(I,J,L)
      FZ(L) = CM(L)*(C(L)-3d0)*RM(I,J,L)
C****
  290 CONTINUE
C****
C**** Calculate new tracer mass and first moments of tracer mass
C****
  300 DO 310 L=1,LMIJ
      RM(I,J,L) = RM(I,J,L) + (FM(L-1)-FM(L))
      RX(I,J,L) = RX(I,J,L) + (FX(L-1)-FX(L))
      RY(I,J,L) = RY(I,J,L) + (FY(L-1)-FY(L))
      RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZ(L-1)-FZ(L))
     *  + 3d0*((CM(L-1)+CM(L))*RM(I,J,L)-MO(I,J,L)*(FM(L-1)+FM(L))))
     *  / (MO(I,J,L)+CM(L-1)-CM(L))
C****
      if ( QLIMIT ) then ! limit tracer gradients
        RXY = abs(RX(I,J,L)) + abs(RY(I,J,L))
        if ( RXY > RM(I,J,L) ) then
          RX(I,J,L) = RX(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
          RY(I,J,L) = RY(I,J,L)*( RM(I,J,L)/(RXY + tiny(RXY)) )
        end if
        if ( abs(RZ(I,J,L)) > RM(I,J,L) )
     *       RZ(I,J,L) = sign(RM(I,J,L), RZ(I,J,L)+0d0)
      end if
C****
      MO(I,J,L) = MO(I,J,L) +  CM(L-1)-CM(L)
         IF(MO(I,J,L).LE.0.)              ICKERR=ICKERR+1
         IF(QLIMIT.AND.RM(I,J,L).LT.0.)   ICKERR=ICKERR+1
  310 CONTINUE
         DO 320 L=1,LMIJ-1
  320    OIJL(I,J,L) = OIJL(I,J,L) + FM(L)
  330 CONTINUE
      END DO
      END DO

C**** IF NO ERROR HAS OCCURRED - RETURN, ELSE STOP
      IF(ICKERR == 0)  RETURN
      DO J=J_0,J_1
        IMIN=1
        IMAX=IM
        IF (J == 1) IMIN=IM
        IF (J == JM) IMAX=1
        DO I=IMIN,IMAX
        LMIJ=LMM(I,J)
        DO L=1,LMIJ
          IF(FOCEAN(I,J).gt.0 .and. MO(I,J,L).LE.0.)  GO TO 800
          IF(QLIMIT .AND. RM(I,J,L).LT.0.) GO TO 810
        END DO
      END DO
      END DO
      WRITE(6,*) 'ERROR CHECK INCONSISTENCY: OADVTZ ',ICKERR
      call stop_model("OAVDTZ",255)

  800 WRITE (6,*) 'MO<0 in OADVTZ:',I,J,L,MO(I,J,L)
  810 WRITE (6,*) 'RM in OADVTZ:',I,J,L,RM(I,J,L)
c      WRITE (6,*) 'C=',(L,C(L),L=0,LMIJ)
      call stop_model("OADVTZ",255)
      END SUBROUTINE OADVTZ

      Subroutine OBDRAG
!@sum  OBDRAG exerts a drag on the Ocean Model's bottom layer
!@+    define OCN_GISS_TURB and idrag=1 in GISS_OTURB module to include tidal
!@+    enhancement of bottom drag
!@auth Gary Russell and Armando Howard
      Use OCEAN, Only: im,jm,lmo,IVNP,J1O, mo,uo,vo, lmu,lmv, dts,
     *                 COSI=>COSIC,SINI=>SINIC
      use domain_decomp_1d, only : getDomainBounds, halo_update, 
     *                             north, south
      USE OCEANR_DIM, only : grid=>ogrid
#ifdef OCN_GISS_TURB
      USE GISS_OTURB, only : taubx,tauby, ! x,y components of velocity flux (m/s)^2
     &                       rhobot       ! ocean bottom in-situ density (kg/m^3)
     &                      ,idrag        ! idrag=1: no explicit tides; 0: otherwise
#endif

      Implicit None
      REAL*8, PARAMETER :: BDRAGX=1d0, SDRAGX=1d-1
      REAL*8,DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO):: UT,VT
      Integer*4 I,J,L,Ip1,Im1
      REAL*8 WSQ

      INTEGER :: J_0, J_0S,J_1S  ; logical :: have_north_pole
      REAL*8 bdragfac   !density times C_D (u^2 + u_t^2)^1/2 (kg/(m^2 s))
#ifdef OCN_GISS_TURB
      REAL*8 taub       !velocity flux (m/s)^2
      REAL*8 taubbyu    !velocity flux divided by velocity (m/s)
#endif

      call getDomainBounds(grid, J_STRT=J_0, 
     *     J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *     have_north_pole=have_north_pole)
C****
C**** UO = UO*(1-x/y)  is approximated by  UO*y/(y+x)  for stability
C****
      call halo_update (grid, mo, FROM=NORTH)
      call halo_update (grid, uo, FROM=NORTH)
      call halo_update (grid, vo, FROM=SOUTH)
#ifdef OCN_GISS_TURB
      call halo_update (grid, taubx, FROM=NORTH)
      call halo_update (grid, tauby, FROM=SOUTH)
      call halo_update (grid, rhobot, FROM=NORTH)
#endif
C**** Save UO,VO into UT,VT which will be unchanged
      DO 10 L=1,LMO
        UT(:,:,L) = UO(:,:,L)
 10     VT(:,:,L) = VO(:,:,L)
C****
C**** Reduce West-East ocean current
C****
C**** Bottom drag in the interior
      DO 120 J=max(J1O,J_0),J_1S
      I=IM
      DO 110 IP1=1,IM
      IF(LMU(I,J) <= 0)  GO TO 110
      L=LMU(I,J)
      WSQ = UT(I,J,L)*UT(I,J,L) + 1d-20 +
     *  .25*(VT(I,J  ,L)*VT(I,J  ,L) + VT(IP1,J  ,L)*VT(IP1,J  ,L)
     *     + VT(I,J-1,L)*VT(I,J-1,L) + VT(IP1,J-1,L)*VT(IP1,J-1,L))
#ifndef OCN_GISS_TURB
        bdragfac=BDRAGX*SQRT(WSQ)
#else
        if(idrag.eq.0) then
          bdragfac=BDRAGX*SQRT(WSQ)
        else
          taub = (taubx(i,j)*taubx(i,j) +
     *      .25*(tauby(i,j)*tauby(i,j) +
     &           tauby(ip1,j)*tauby(ip1,j)
     *         + tauby(i,j-1)*tauby(i,j-1) +
     &           tauby(ip1,j-1)*tauby(ip1,j-1)))
     &      **0.5d0
          taubbyu=taub/SQRT(WSQ)
          bdragfac=0.5*(rhobot(i,j)+rhobot(ip1,j))*taubbyu
        endif
#endif
        UO(I,J,L) = UO(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *             (MO(I,J,L)+MO(IP1,J,L) + DTS*bdragfac*2d0)
  110 I=IP1
  120 CONTINUE
C**** Bottom drag at the poles
C     IF(LMU(1,1 or JM) <= 0)  GO TO
      if (have_north_pole) then
        L=LMU(1,JM)
        WSQ = UT(IM,JM,L)**2 + UT(IVNP,JM,L)**2 + 1d-20
#ifndef OCN_GISS_TURB
        bdragfac=BDRAGX*SQRT(WSQ)
#else
        if(idrag.eq.0) then
          bdragfac=BDRAGX*SQRT(WSQ)
        else
          taub = SQRT(taubx(im,jm)*taubx(im,jm) +
     &           taubx(ivnp,jm)*taubx(ivnp,jm))
          taubbyu=taub/SQRT(WSQ)
          bdragfac=rhobot(1,jm)*taubbyu
        endif
#endif
          UO(IM,JM,L) = UO(IM,JM,L) * MO(1,JM,L) /
     /                             (MO(1,JM,L) + DTS*bdragfac)
          UO(IVNP,JM,L) = UO(IVNP,JM,L) * MO(1,JM,L) /
     /                             (MO(1,JM,L) + DTS*bdragfac)
      end if
C****
C**** Reduce South-North ocean current
C****
      Do 240 J=max(J1O,J_0),J_1S
      If (J==JM-1)  GoTo 220
C**** Bottom drag away from north pole
      IM1=IM
      DO 210 I=1,IM
      IF(LMV(I,J) <= 0)  GO TO 210
      L=LMV(I,J)
      WSQ = VT(I,J,L)*VT(I,J,L) + 1d-20 +
     *  .25*(UT(IM1,J+1,L)*UT(IM1,J+1,L) + UT(I,J+1,L)*UT(I,J+1,L)
     *     + UT(IM1,J  ,L)*UT(IM1,J  ,L) + UT(I,J  ,L)*UT(I,J  ,L))
#ifndef OCN_GISS_TURB
        bdragfac=BDRAGX*SQRT(WSQ)
#else
        if(idrag.eq.0) then
          bdragfac=BDRAGX*SQRT(WSQ)
        else
          taub = SQRT(tauby(i,j)*tauby(i,j) +
     *      .25*(taubx(im1,j+1)*taubx(im1,j+1) +
     &           taubx(i,j+1)*taubx(i,j+1)
     *         + taubx(im1,j)*taubx(im1,j) +
     &           taubx(i,j)*taubx(i,j)))
          taubbyu=taub/SQRT(WSQ)
          bdragfac=0.5d0*(rhobot(i,j)+rhobot(i,j+1))*taubbyu
        endif
#endif
        VO(I,J,L) = VO(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *             (MO(I,J,L)+MO(I,J+1,L) + DTS*bdragfac*2d0)
  210 IM1=I
      GoTo 240
C**** Bottom drag near north pole
  220 Im1=IM
      Do 230 I=1,IM
      If (LMV(I,JM-1) <= 0)  GoTo 230
      L = LMV(I,JM-1)
      WSQ = VT(I,JM-1,L)*VT(I,JM-1,L) + 1d-20 +
     +   .5*(UT(IM,JM,L)*COSI(I) + UT(IVNP,JM,L)*SINI(I))**2 +
     +  .25*(UT(Im1,JM-1,L)*UT(Im1,JM-1,L) + UT(I,JM-1,L)*UT(I,JM-1,L))
#ifndef OCN_GISS_TURB
        bdragfac=BDRAGX*SQRT(WSQ)
#else
        if(idrag.eq.0) then
          bdragfac=BDRAGX*SQRT(WSQ)
        else
          taub = SQRT(tauby(i,jm-1)*tauby(i,jm-1) +
     +       .5*(taubx(im,jm)*cosi(i) +
     &           taubx(ivnp,jm)*sini(i))**2 +
     +      .25*(taubx(im1,jm-1)*taubx(im1,jm-1) +
     &           taubx(i,jm-1)*taubx(i,jm-1)))
          taubbyu=taub/SQRT(WSQ)
          bdragfac=0.5d0*(rhobot(i,jm-1)+rhobot(1,jm))*taubbyu
        endif
#endif
        VO(I,JM-1,L) = VO(I,JM-1,L) * (MO(I,JM-1,L)+MO(1,JM,L)) /
     *              (MO(I,JM-1,L)+MO(1,JM,L) + DTS*bdragfac*2)
  230 Im1=I
  240 Continue
      RETURN
C****
      END Subroutine OBDRAG

      SUBROUTINE OCOAST
!@sum OCOAST reduces the horizontal perpendicular gradients of tracers
!@sum in coastline ocean grid boxes
!@auth Gary Russell
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE OCEAN, only : im,jm,dts,lmm,gxmo,gymo,sxmo,symo
#ifdef TRACERS_OCEAN
     *     ,txmo,tymo
      Use OCN_TRACER_COM, Only: tracerlist
#endif
!      use domain_decomp_1d, only : grid, get
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      INTEGER I,IM1,IP1,J,LMIN,L,N
      REAL*8 REDUCE

      integer :: J_0S, J_1S

      call getDomainBounds(grid, J_STRT_SKP=J_0S, J_STOP_SKP=J_1S)


      REDUCE = 1d0 - DTS/(SECONDS_PER_DAY*2d1)
C**** Reduce West-East gradient of tracers
      IM1=IM-1
      I=IM
      DO 120 J=J_0S,J_1S
      DO 120 IP1=1,IM
      LMIN = MIN(LMM(IM1,J),LMM(IP1,J)) + 1
      DO 110 L=LMIN,LMM(I,J)
        GXMO(I,J,L) = GXMO(I,J,L)*REDUCE
        SXMO(I,J,L) = SXMO(I,J,L)*REDUCE
#ifdef TRACERS_OCEAN
        DO N = 1,tracerlist%getsize()
          TXMO(I,J,L,N) = TXMO(I,J,L,N) *REDUCE
        END DO
#endif
 110  CONTINUE
      IM1=I
  120 I=IP1
C**** Reduce South-North gradient of tracers
      DO 220 J=J_0S,J_1S
      DO 220 I=1,IM
      LMIN = MIN(LMM(I,J-1),LMM(I,J+1)) + 1
      DO 210 L=LMIN,LMM(I,J)
        GYMO(I,J,L) = GYMO(I,J,L)*REDUCE
        SYMO(I,J,L) = SYMO(I,J,L)*REDUCE
#ifdef TRACERS_OCEAN
        DO N = 1,tracerlist%getsize()
          TYMO(I,J,L,N) = TYMO(I,J,L,N) *REDUCE
        END DO
#endif
 210  CONTINUE
  220 CONTINUE
      RETURN
      END SUBROUTINE OCOAST

      SUBROUTINE OSTRES
!@sum OSTRES applies the atmospheric surface stress over open ocean
!@sum and the sea ice stress to the layer 1 ocean velocities
!@auth Gary Russell

      USE OCEAN, only : IMO=>IM,JMO=>JM
     *     , IVNP, UO,VO, MO,DXYSO,DXYNO,DXYVO
     *     , LMU,LMV, COSIC,SINIC

      USE DOMAIN_DECOMP_1D, only : getDomainBounds, halo_update, north

      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : oDMUA,oDMVA, oDMUI,oDMVI

      IMPLICIT NONE
      INTEGER I,J,IP1

C****
C**** All stress now defined for whole box, not just ocn or ice fraction
C**** FLUXCB  DMUA(1)  U momentum downward into open ocean (kg/m*s)
C****         DMVA(1)  V momentum downward into open ocean (kg/m*s)
C****         DMUA(2,JM,1)  polar atmo. mass slowed to zero (kg/m**2)
C****         DMUI     U momentum downward from sea ice (kg/m*s)
C****         DMVI     V momentum downward from sea ice (kg/m*s)

      integer :: J_0, J_1, J_0S, J_1S  ; logical :: have_north_pole

      call getDomainBounds(ogrid, J_STRT=J_0, J_STOP=J_1,
     *                 J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *                 have_north_pole=have_north_pole)

C****
C**** Surface stress is applied to U component
C****
      I=IMO
      DO J=J_0S,J_1S
      DO IP1=1,IMO
        IF(LMU(I,J).gt.0.)  UO(I,J,1) = UO(I,J,1) +
     *       (oDMUA(I,J) + oDMUA(IP1,J) + 2d0*oDMUI(I,J)) /
     *       (  MO(I,J,1) +   MO(IP1,J,1))
        I=IP1
      END DO
      END DO
      if (have_north_pole) then
        UO(IMO ,JMO,1) = UO(IMO ,JMO,1) + oDMUA(1,JMO)/MO(1,JMO,1)
        UO(IVNP,JMO,1) = UO(IVNP,JMO,1) + oDMVA(1,JMO)/MO(1,JMO,1)
      end if
C****
C**** Surface stress is applied to V component
C****
      call halo_update(ogrid, odmva, from=north)
      call halo_update(ogrid,    mo, from=north)
      DO J=J_0S,min(J_1S,JMO-2)
      DO I=1,IMO
        IF(LMV(I,J).GT.0.)  VO(I,J,1) = VO(I,J,1) +
     *       (oDMVA(I,J  )*DXYNO(J) + oDMVA(I,J+1)*DXYSO(J+1)
     *      + oDMVI(I,J)*DXYVO(J))  !!  2d0*oDMVI(I,J)*DXYVO(J) - error
     * / (MO(I,J,1)*DXYNO(J) + MO(I,J+1,1)*DXYSO(J+1))
      END DO
      END DO
C**** Surface stress is applied to V component at the North Pole
      if (have_north_pole) then
      DO I=1,IMO
        VO(I,JMO-1,1) = VO(I,JMO-1,1) +
     *    (oDMVA(I,JMO-1)*DXYNO(JMO-1)+
     *    (oDMVA(1,JMO)*COSIC(I) - oDMUA(1,JMO)*SINIC(I))*DXYSO(JMO)
     *   + oDMVI(I,JMO-1)*DXYVO(JMO-1)) /
     *  (MO(I,JMO-1,1)*DXYNO(JMO-1) + MO(I,JMO,1)*DXYSO(JMO))
      END DO
      end if

      RETURN
      END SUBROUTINE OSTRES

      SUBROUTINE GROUND_OC
!@sum  GROUND_OC adds vertical fluxes into the ocean
!@auth Gary Russell/Gavin Schmidt
      USE CONSTANT, only : grav
      USE OCEAN, only : IMO=>IM,JMO=>JM, LMO,LMM
     *     , MO,G0M,S0M, FOCEAN, GZMO, IMAXJ,DXYPO,BYDXYPO, OPRESS
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : ogrid
      USE SEAICE, only : Ei,FSSS
      USE OFLUXES, only : oRSI, oSOLARw,oSOLARi, oE0, oEVAPOR
     *     , oRUNOSI, oERUNOSI, oSRUNOSI
     *     , oFLOWO, oEFLOWO, oAPRESS
     *     , oMELTI, oEMELTI, oSMELTI
     *     , oDMSI, oDHSI, oDSSI
     *     , ocnice
      USE ODIAG, only : oij=>oij_loc,ij_srhflx,ij_srwflx,ij_srhflxi
     *     ,ij_srwflxi,ij_srsflxi,ij_ervr,ij_mrvr,ij_ocnfr
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
      Use OCEAN, Only: TRMO
#endif
#ifdef TRACERS_OCEAN
      Use OFLUXES, Only:  oDTRSI
#ifdef TRACERS_WATER
     *     , oTRFLOWO, oTREVAPOR
     *     , oTRUNOSI, oTRMELTI
#endif
#ifdef TRACERS_DRYDEP
     *     , oTRDRYDEP
#endif
!#ifdef TRACERS_GASEXCH_ocean
!     *     , oTRGASEX
!#endif
#endif
      IMPLICIT NONE
      INTEGER I,J,L,n
      REAL*8 DXYPJ,BYDXYPJ,RUNO,RUNI,ERUNO,ERUNI,SROX(2),G0ML(LMO)
     *     ,MO1,SO1,ROICE,DMOO,DMOI,DEOO,DEOI,GZML(LMO),SRUNO,SRUNI,DSOO
     *     ,DSOI,POCEAN,POICE,P0L,S0L,G0L,TF0,SI0,EI0,DM0(LMO),DE0(LMO)
     *     ,DS0(LMO),GF0,GFREZS,TFREZS,TEMGSP,SHCGS,GF00
#ifdef TRACERS_OCEAN
      REAL*8, DIMENSION(tracerlist%getsize()) :: TRUNO,TRUNI,DTROO,
     &              DTROI,TRO1,FRAC
      REAL*8, DIMENSION(tracerlist%getsize(),LMO) :: DTR0
#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif
#endif

      integer ::  J_1, J_0
      logical :: have_south_pole, have_north_pole
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), dimension(:), pointer :: tlist
      type(ocn_tracer_entry), pointer :: entry
#endif

#ifdef TRACERS_OCEAN
      tlist=>tracerlist%getdata()
#endif
      call getDomainBounds(ogrid, J_STOP=J_1, J_STRT=J_0,
     *                have_north_pole=have_north_pole,
     *                have_south_pole=have_south_pole)

#ifdef TRACERS_OCEAN
C**** define tracer behaviour for ice formation
        do n=1,tracerlist%getsize()
#ifdef TRACERS_WATER
#ifdef TRACERS_SPECIAL_O18
          FRAC(n)=fracls(n)
#else
          entry=>tracerlist%at(n)
          FRAC(n)=entry%trw0    ! removal based on default conc.
#endif
#else
          FRAC(n)=0.            ! no removal of non-water tracers
#endif
        end do
#endif
C****
C**** Add surface source of fresh water and heat
C****
      DO J=J_0,J_1
        DXYPJ=DXYPO(J)
        BYDXYPJ=BYDXYPO(J)
      DO I=1,IMAXJ(J)
      IF(FOCEAN(I,J).gt.0.) THEN
        ROICE = oRSI(I,J)
        POCEAN=FOCEAN(I,J)*(1.-ROICE)
        POICE =FOCEAN(I,J)*ROICE
        DXYPJ=DXYPJ*FOCEAN(I,J)     ! adjust areas for completeness
        BYDXYPJ=BYDXYPJ/FOCEAN(I,J) ! no change to results
C**** set mass & energy fluxes (incl. river/sea ice runoff + basal flux)
        RUNO = oFLOWO(I,J) + oMELTI(I,J) - oEVAPOR(I,J)    !  kg/m^2
        RUNI = oFLOWO(I,J) + oMELTI(I,J) + oRUNOSI(I,J)    !  kg/m^2
        ERUNO=oEFLOWO(I,J) + oEMELTI(I,J) + oE0(I,J)       !  J/m^2
        ERUNI=oEFLOWO(I,J) + oEMELTI(I,J) + oERUNOSI(I,J)  !  J/m^2
        SRUNO=oSMELTI(I,J)                                 !  kg/m^2
        SRUNI=oSMELTI(I,J) + oSRUNOSI(I,J)                 !  kg/m^2
        G0ML(:) =  G0M(I,J,:)
        GZML(:) = GZMO(I,J,:)
        SROX(1)=oSOLARw(I,J) ! open water    J/m^2
        SROX(2)=oSOLARi(I,J) ! through ice   J/m^2
        MO1 = MO(I,J,1)
        SO1 = S0M(I,J,1)
#ifdef TRACERS_OCEAN
        TRO1(:) = TRMO(I,J,1,:)
#ifdef TRACERS_WATER
        TRUNO(:)=oTRFLOWO(:,I,J)+oTRMELTI(:,I,J)-oTREVAPOR(:,I,J)
#ifdef TRACERS_DRYDEP
     *       + otrdrydep(:,1,i,j)  !  kg/m^2
#endif
        TRUNI(:)=oTRFLOWO(:,I,J)+ oTRMELTI(:,I,J) + oTRUNOSI(:,I,J)
#else
! default freshwater tracer
        TRUNO(:)=tlist%trw0*(RUNO-SRUNO)
        TRUNI(:)=tlist%trw0*(RUNI-SRUNI)
#endif
#endif

        CALL OSOURC(ROICE,MO1,G0ML,GZML,SO1,DXYPJ,BYDXYPJ,LMM(I,J),RUNO
     *         ,RUNI,ERUNO,ERUNI,SRUNO,SRUNI,SROX,
#ifdef TRACERS_OCEAN
     *         TRO1,TRUNO,TRUNI,DTROO,DTROI,FRAC,
     &         tracerlist%getsize(),
#endif
     *         DMOO,DEOO,DMOI,DEOI,DSOO,DSOI)

C**** update ocean variables
          MO(I,J,1) = MO1
         S0M(I,J,1) = SO1
         G0M(I,J,:) = G0ML(:)
        GZMO(I,J,:) = GZML(:)
c        write(*,*) "GZMO"
#ifdef TRACERS_OCEAN
        TRMO(I,J,1,:) = TRO1(:)
#endif

C**** Do sweep through lower layers for any possible ice formation
C**** Add evenly over open ocean and ice covered areas
        DM0(:)=0 ; DS0(:)=0. ; DE0(:)=0.
#ifdef TRACERS_OCEAN
        DTR0(:,:)=0.
#endif
        P0L=MO(I,J,1)*GRAV
        DO L=2,LMM(I,J)
          G0L=G0M(I,J,L)/(MO(I,J,L)*DXYPJ)
          S0L=S0M(I,J,L)/(MO(I,J,L)*DXYPJ)
          P0L=P0L + MO(I,J,L)*GRAV*.5
          GF00=GFREZS(S0L)
          GF0=GF00-SHCGS(GF00,S0L)*8.19d-8*P0L  ! ~GFREZSP(S0L,P0L)
          IF(G0L.lt.GF0) THEN
            TF0=TFREZS(S0L)-7.53d-8*P0L
            SI0=FSSS*S0L
            EI0=Ei(TF0,SI0*1d3)
            DM0(L)=MO(I,J,L)*(G0L-GF0)/(EI0-GF0)
            DE0(L)=EI0*DM0(L)
            DS0(L)=SI0*DM0(L)
#ifdef TRACERS_OCEAN
            DTR0(:,L)=TRMO(I,J,L,:)*FRAC(:)*
     *               (DM0(L)-DS0(L))/(MO(I,J,L)*DXYPJ-S0M(I,J,L))
            TRMO(I,J,L,:)=TRMO(I,J,L,:)-DTR0(:,L)*DXYPJ
#endif
            MO(I,J,L) = MO(I,J,L)-DM0(L)
            S0M(I,J,L)=S0M(I,J,L)-DS0(L)*DXYPJ
            G0M(I,J,L)=G0M(I,J,L)-DE0(L)*DXYPJ
          END IF
          P0L=P0L + MO(I,J,L)*GRAV*.5
        END DO
C**** Store mass/energy/salt/tracer fluxes for formation of sea ice
c        write(*,*) "store fluxes"
        oDMSI(1,I,J)=DMOO+SUM(DM0)  !  kg/m^2
        oDMSI(2,I,J)=DMOI+SUM(DM0)  !  kg/m^2
        oDHSI(1,I,J)=DEOO+SUM(DE0)  !  J/m^2
        oDHSI(2,I,J)=DEOI+SUM(DE0)  !  J/m^2
        oDSSI(1,I,J)=DSOO+SUM(DS0)  !  kg/m^2
        oDSSI(2,I,J)=DSOI+SUM(DS0)  !  kg/m^2

#ifdef TRACERS_OCEAN
        if(ocnice%ntm == tracerlist%getsize()) then
          oDTRSI(:,1,I,J)=DTROO(:)+SUM(DTR0(:,:),DIM=2)
          oDTRSI(:,2,I,J)=DTROI(:)+SUM(DTR0(:,:),DIM=2)
        endif
#endif

C**** Calculate pressure anomaly at ocean surface (and scale for areas)
C**** Updated using latest sea ice (this ensures that total column mass
C**** is consistent for OGEOZ calculation).
        OPRESS(I,J) = oAPRESS(I,J)+GRAV*(
     *       (1.-oRSI(I,J))*oDMSI(1,I,J) + oRSI(I,J)*oDMSI(2,I,J))

        OIJ(I,J,IJ_OCNFR)  = OIJ(I,J,IJ_OCNFR) + 1.

C**** Set some ocean diagnostics of net fluxes (downward +ve)
C**** This includes atm/oc + si/oc, rivers + icebergs are separate 
        OIJ(I,J,IJ_SRHFLX)  = OIJ(I,J,IJ_SRHFLX)  +     ! net heat
     *       (ERUNO-oDHSI(1,I,J))*POCEAN +
     *       (ERUNI-oDHSI(2,I,J))*POICE  - oEFLOWO(I,J)

        OIJ(I,J,IJ_SRWFLX)  = OIJ(I,J,IJ_SRWFLX)  +     ! net fw
     *       (RUNO-SRUNO-oDMSI(1,I,J)+oDSSI(1,I,J))*POCEAN+
     *       (RUNI-SRUNI-oDMSI(2,I,J)+oDSSI(2,I,J))*POICE - oFLOWO(I,J)

        OIJ(I,J,IJ_SRHFLXI) = OIJ(I,J,IJ_SRHFLXI) + ! ht from ice
     *       oEMELTI(I,J) - oDHSI(1,I,J)*POCEAN + 
     *       (oERUNOSI(I,J)- oDHSI(2,I,J))*POICE

        OIJ(I,J,IJ_SRWFLXI) = OIJ(I,J,IJ_SRWFLXI) + ! fw from ice
     *       oMELTI(I,J) - oSMELTI(I,J) -
     *       (oDMSI(1,I,J)-oDSSI(1,I,J))*POCEAN +
     *       (oRUNOSI(I,J)-oSRUNOSI(I,J)-oDMSI(2,I,J)+oDSSI(2,I,J))
     *       *POICE

        OIJ(I,J,IJ_SRSFLXI) = OIJ(I,J,IJ_SRSFLXI) + ! salt from ice
     *       (SRUNO-oDSSI(1,I,J))*POCEAN +
     *       (SRUNI-oDSSI(2,I,J))*POICE

        OIJ(I,J,IJ_ERVR) = OIJ(I,J,IJ_ERVR) + oEFLOWO(I,J) ! energy from rivers
        OIJ(I,J,IJ_MRVR) = OIJ(I,J,IJ_MRVR) + oFLOWO(I,J) ! fw from rivers

        END IF
      END DO
      END DO
      if(have_south_pole) OPRESS(2:IMO,1)  = OPRESS(1,1)
      if(have_north_pole) OPRESS(2:IMO,JMO) = OPRESS(1,JMO)

      RETURN
      END SUBROUTINE GROUND_OC

      SUBROUTINE OSOURC (ROICE,MO,G0ML,GZML,S0M,DXYPJ,BYDXYPJ,LMIJ,RUNO
     *     ,RUNI,ERUNO,ERUNI,SRUNO,SRUNI,SROX,
#ifdef TRACERS_OCEAN
     *     TROM,TRUNO,TRUNI,DTROO,DTROI,FRAC,
     &     numtracers,    !df: due to bug in gfortran 4.8.3 and earlier
#endif
     *     DMOO,DEOO,DMOI,DEOI,DSOO,DSOI)
!@sum  OSOURC applies fluxes to ocean in ice-covered and ice-free areas
!@auth Gary Russell/Gavin Schmidt

      USE SW2OCEAN, only : lsrpd,fsr,fsrz
      USE SEAICE, only : fsss, Ei

#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist
#endif

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: ROICE,DXYPJ,BYDXYPJ,RUNO,RUNI,ERUNO,ERUNI
     *     ,SROX(2),SRUNO,SRUNI
      INTEGER, INTENT(IN) :: LMIJ
      REAL*8, INTENT(INOUT) :: MO,G0ML(LSRPD),GZML(LSRPD),S0M
      REAL*8, INTENT(OUT) :: DMOO,DMOI,DEOO,DEOI,DSOO,DSOI
      REAL*8 MOO,GOO,GMOO,GMOI,MOI,GOI,SMOO,SMOI,SOO,SOI,GFOO,GFOI,TFOO
     *     ,TFOI,SIOO,SIOI
      REAL*8 GFREZS,TFREZS,TSOL
      INTEGER L,LSR,N

#ifdef TRACERS_OCEAN
      integer, intent(in) :: numtracers
      REAL*8, DIMENSION(numtracers), INTENT(INOUT) :: TROM,FRAC
      REAL*8, DIMENSION(numtracers), INTENT(IN) :: TRUNO,TRUNI
      REAL*8, DIMENSION(numtracers), INTENT(OUT) :: DTROO,DTROI
      REAL*8, DIMENSION(tracerlist%getsize()) :: TMOO,TMOI
#endif

      DMOO=0. ; DEOO=0. ; DMOI=0. ; DEOI=0. ; DSOO=0. ; DSOI=0.

#ifdef TRACERS_OCEAN
      DTROI(:) = 0. ; DTROO(:) = 0.
#endif

      LSR = MIN(LSRPD,LMIJ)
C****
C**** Open Ocean
C****
      MOO  = MO + RUNO
      GMOO = G0ML(1)*BYDXYPJ + ERUNO
      SMOO = S0M*BYDXYPJ + SRUNO
#ifdef TRACERS_OCEAN
      TMOO(:) = TROM(:)*BYDXYPJ+TRUNO(:)
#endif

      IF (ROICE.lt.1d0) THEN
C**** Remove insolation from layer 1 that goes to lower layers
      IF (LSR.gt.1) GMOO = GMOO - SROX(1)*FSR(2)
      GOO  = GMOO/MOO
      SOO  = SMOO/MOO
      GFOO = GFREZS(SOO)
      IF(GOO.lt.GFOO) THEN
C**** Open ocean is below freezing, calculate
C**** DMOO = mass of ocean that freezes over open fraction from
C**** GOO*MOO = GFOO*(MOO-DMOO) + Ei(TFOO,SIOO*1d3)*DMOO
        TFOO = TFREZS(SOO)
        SIOO = FSSS*SOO
        DMOO = MOO*(GOO-GFOO)/(Ei(TFOO,SIOO*1d3)-GFOO)
        DEOO = Ei(TFOO,SIOO*1d3)*DMOO
        DSOO = SIOO*DMOO
#ifdef TRACERS_OCEAN
        DTROO(:) = TMOO(:)*FRAC(:)*(DMOO-DSOO)/(MOO-SMOO)
#endif
      END IF
      END IF
C****
C**** Ocean underneath the ice
C****
      MOI  = MO + RUNI
      GMOI = G0ML(1)*BYDXYPJ + ERUNI
      SMOI = S0M*BYDXYPJ + SRUNI
#ifdef TRACERS_OCEAN
      TMOI(:) = TROM(:)*BYDXYPJ+TRUNI(:)
#endif
      IF(ROICE.gt.0.) THEN
C**** Remove insolation from layer 1 that goes to lower layers
        IF (LSR.gt.1) GMOI = GMOI - SROX(2)*FSR(2)

        GOI  = GMOI/MOI
        SOI  = SMOI/MOI
        GFOI = GFREZS(SOI)
        IF(GOI.LT.GFOI) THEN
C**** Ocean underneath the ice is below freezing, calculate
C**** DMOI = mass of ocean that freezes under sea ice fraction from
C**** GOI*MOI = GFOI*(MOI-DMOI) + Ei(TFOI,SIOI*1d3)*DMOI
          TFOI = TFREZS(SOI)
          SIOI = FSSS*SOI
          DMOI = MOI*(GOI-GFOI)/(Ei(TFOI,SIOI*1d3)-GFOI)
          DEOI = Ei(TFOI,SIOI*1d3)*DMOI
          DSOI = SIOI*DMOI
#ifdef TRACERS_OCEAN
          DTROI(:) = TMOI(:)*FRAC(:)*(DMOI-DSOI)/(MOI-SMOI)
#endif
        END IF
      END IF
C**** Update first layer variables
      MO     =  (MOI-DMOI)*ROICE + (1.-ROICE)*( MOO-DMOO)
      G0ML(1)=((GMOI-DEOI)*ROICE + (1.-ROICE)*(GMOO-DEOO))*DXYPJ
      S0M    =((SMOI-DSOI)*ROICE + (1.-ROICE)*(SMOO-DSOO))*DXYPJ
#ifdef TRACERS_OCEAN
      TROM(:)=((TMOI(:)-DTROI(:))*ROICE + (1.-ROICE)*(TMOO(:)-DTROO(:)))
     *     *DXYPJ
#endif
C**** add insolation to lower layers
      TSOL=(SROX(1)*(1.-ROICE)+SROX(2)*ROICE)*DXYPJ
      DO L=2,LSR-1
        G0ML(L)=G0ML(L)+TSOL*(FSR(L)-FSR(L+1))
        GZML(L)=GZML(L)+TSOL*FSRZ(L)
      END DO
      G0ML(LSR) = G0ML(LSR) + TSOL*FSR (LSR)
      GZML(LSR) = GZML(LSR) + TSOL*FSRZ(LSR)
C****
      RETURN
      END SUBROUTINE OSOURC

      SUBROUTINE PRECIP_OC(atmocn,iceocn)
!@sum  PRECIP_OC driver for applying precipitation to ocean fraction
!@auth Gary Russell/Gavin Schmidt
      USE OCEAN, only : imo=>im,jmo=>jm
     *     , mo,g0m,s0m,focean,imaxj,dxypo
#ifdef TRACERS_OCEAN
     *     , trmo,mosv0
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry
#endif
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE OCEANR_DIM, only : oGRID
      USE OFLUXES, only : oRSI, oPREC, oEPREC
     *     , oRUNPSI, oSRUNPSI, oERUNPSI
#if (defined TRACERS_OCEAN) && (defined TRACERS_WATER)
     *     , oTRPREC, oTRUNPSI
#endif
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
c
      INTEGER I,J
      integer :: J_0, J_1
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), dimension(:), pointer :: tlist
#endif

#ifdef TRACERS_OCEAN
      tlist=>tracerlist%getdata()
#endif
      call getDomainBounds(ogrid, J_STRT=J_0, J_STOP=J_1)

C**** save surface variables before any fluxes are added
      if(ogrid%have_domain) CALL KVINIT

C**** Convert fluxes on atmospheric grid to oceanic grid
      CALL AG2OG_precip(atmocn,iceocn)
C****
      ocean_processors_only: if(ogrid%have_domain) then

#ifdef TRACERS_OCEAN
      ! save 3D mass before all source/sink terms
      if(allocated(mosv0)) mosv0 = mo
#endif

      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          IF(FOCEAN(I,J).gt.0. .and. oPREC(I,J).gt.0.)  THEN
            MO (I,J,1)= MO(I,J,1) + ((1d0-oRSI(I,J))*oPREC(I,J) +
     *           oRSI(I,J)*oRUNPSI(I,J))*FOCEAN(I,J)
            G0M(I,J,1)=G0M(I,J,1)
     &           +((1d0-oRSI(I,J))*(oEPREC(I,J)*dxypo(j))
     &           +oRSI(I,J)*(oERUNPSI(I,J)*dxypo(j)))*FOCEAN(I,J)
            S0M(I,J,1)=S0M(I,J,1)
     &           +oRSI(I,J)*(oSRUNPSI(I,J)*dxypo(j))*FOCEAN(I,J)
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
            TRMO(I,J,1,:)=TRMO(I,J,1,:)
     &           +((1d0-oRSI(I,J))*(oTRPREC(:,I,J)*dxypo(j))
     &           +oRSI(I,J)*(oTRUNPSI(:,I,J)*dxypo(j)))*FOCEAN(I,J)
#else
#ifndef TRACERS_OceanBiology
            TRMO(I,J,1,:)=TRMO(I,J,1,:)+tlist%trw0*((1d0-oRSI(I,J))*
     .         oPREC(I,J)+oRSI(I,J)*oRUNPSI(I,J))*FOCEAN(I,J)*DXYPO(J)
#endif
#endif
#endif
          END IF
        END DO
      END DO

#ifdef STANDALONE_OCEAN
! surface salinity restoration
      call restore_surface_salinity2!(atmocn)
#endif

      endif ocean_processors_only

C**** Convert ocean surface temp to atmospheric SST array
      CALL TOC2SST(atmocn)

      RETURN
      END SUBROUTINE PRECIP_OC

      SUBROUTINE ODIFF (DTDIFF)
C???? ESMF-exception - ODIFF currently works with global arrays
!@sum  ODIFF applies Wajsowicz horizontal viscosity to velocities
!@auth Gavin Schmidt
C****
C**** ODIFF calculates horizontal Wajsowicz viscosity terms in momentum
C**** equations implicitly using ADI method and assumes no slip/free
C**** slip conditions at the side. K_h (m^2/s) may vary spatially
C**** based on Munk length though must remain isotropic.
C**** (If longitudinal variation is wanted just make K arrays K(I,J))
C**** FSLIP = 0 implies no slip conditions, = 1 implies free slip
C**** Mass variation is included
C****
      use OCEAN, only: FSLIP
      USE OCEAN, only : im,jm,lmo,mo,uo,vo,
     *  IVNP,UONP,VONP, COSU,SINU, COSI=>COSIC,SINI=>SINIC,
     *  lmu,lmv,dxpo,dypo,dxvo,dyvo,bydxypo
      USE OCEAN_DYN, only : dh
      USE TRIDIAG_MOD, only : tridiag, tridiag_new
#ifdef ODIFF_TRIDIAG_CYCLIC
      USE TRIDIAG_MOD, only : tridiag_cyclic
#endif
      USE DOMAIN_DECOMP_1D, ONLY : GETDomainBounds, AM_I_ROOT
      USE OCEANR_DIM, only : grid=>ogrid
      USE DOMAIN_DECOMP_1D, ONLY : HALO_UPDATE, NORTH, SOUTH, broadcast
      USE MODEL_COM, only: nstep=>itime
      USE OCEAN, only: BYDXYV, KHP,KHV,TANP,TANV,BYDXV,BYDXP,BYDYV,
     *     BYDYP,UXA,UXB,UXC,UYA,UYB,UYC,VXA,VXB,VXC,VYA,VYB,VYC
      IMPLICIT NONE

      REAL*8, DIMENSION(IM,grid%j_strt_halo:grid%j_stop_halo) ::
     *      FUX,FUY,FVX,FVY,BYMU,BYMV

C**** Local variables

      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::
     *     AU, BU, CU, RU, UU
      REAL*8, ALLOCATABLE, DIMENSION(:,:) ::
     *     AV, BV, CV, RV, UV
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     AU3D, BU3D, CU3D, RU3D, UU3D
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) ::
     *     AV3D, BV3D, CV3D, RV3D, UV3D
      REAL*8, INTENT(IN) :: DTDIFF
      REAL*8 DSV,DSP,VLAT,DT2,DTU,DTV,VX,VY,VT,UT,UX,UY
      INTEGER I,J,L,IP1,IM1,II

!     domain decomposition
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

c**** Extract domain decomposition info
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               J_STRT_HALO = J_0H,   J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C****

C**** End of initialization from first call to ODIFF
C****
C**** Solve diffusion equations semi implicitly
C****
      DT2=DTDIFF*5d-1     ! half time-step
C**** Store North Pole velocity components, they will not be changed
      if(have_north_pole) then
        UONP(:) = UO(IM  ,JM,:)
        VONP(:) = UO(IVNP,JM,:)
      end if

!     also may later need HALO for DXPO(N), DXVO(S)

      CALL HALO_UPDATE(grid,DH,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,MO,
     *                 FROM=NORTH)

      CALL HALO_UPDATE(grid,UO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,UO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,VO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=NORTH)
      CALL HALO_UPDATE(grid,VO (:,grid%j_strt_halo:grid%j_stop_halo,:) ,
     *                 FROM=SOUTH)

C****

      allocate( AU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( CU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( RU(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( UU(IM,grid%j_strt_halo:grid%j_stop_halo) )

      allocate( AV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( BV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( CV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( RV(IM,grid%j_strt_halo:grid%j_stop_halo) )
      allocate( UV(IM,grid%j_strt_halo:grid%j_stop_halo) )

      allocate( AU3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( BU3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( CU3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( RU3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( UU3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )

      allocate( AV3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( BV3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( CV3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( RV3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )
      allocate( UV3D(IM,grid%j_strt_halo:grid%j_stop_halo,lmo) )

      DO L=1,LMO
C**** Calculate rotating polar velocities from UONP and VONP
      if(have_north_pole) then
        UO(:,JM,L) = UONP(L)*COSU(:) + VONP(L)*SINU(:)
        VO(:,JM,L) = VONP(L)*COSI(:) - UONP(L)*SINI(:)
      end if
C**** Save (0.5*) mass reciprical for velocity points
      DO J=J_0S,J_1S
        I=IM
        DO IP1=1,IM
          IF (L.LE.LMU(I,J)) BYMU(I,J) = 1./(MO(I,J,L)+MO(IP1,J,L))
          IF (L.LE.LMV(I,J)) BYMV(I,J) = 1./(MO(I,J,L)+MO(I,J+1,L))
          I=IP1
        END DO
      END DO

#ifdef ODIFF_FIXES_2017
      if(have_north_pole) then
        call polevel(uo(1,j_0h,l),vo(1,j_0h,l),l)
      endif
#endif

      if( HAVE_NORTH_POLE ) then
        IF (L.LE.LMU(1,JM)) BYMU(1,JM) = 1./MO(1,JM,L)
      endif
C**** Calculate Wajsowicz boundary terms
C**** Need dv/dy,tv,dv/dx for u equation, du/dy,tu,du/dx for v equation
      FUX=0             ! flux in U equation at the x_+ boundary
      FUY=0             ! flux in U equation at the y_+ boundary
      FVX=0             ! flux in V equation at the x_+ boundary
      FVY=0             ! flux in V equation at the y_+ boundary
      DO J=J_0, J_1S
        IM1=IM
        DO I=1,IM
          UT=0          ! mean u*tan on x_+ boundary for V equation
          UY=0          ! mean du/dx on y_+ boundary for V equation
          UX=0          ! mean du/dy on x_+ boundary for V equation
          IF (L.LE.LMU(I  ,J+1)) THEN
            UT=     UO(I  ,J+1,L)*TANP(J+1)
            UX=     UO(I  ,J+1,L)
            UY=     UO(I  ,J+1,L)
          END IF
          IF (L.LE.LMU(I  ,J  )) THEN
            UT=UT + UO(I  ,J  ,L)*TANP(J  )
            UY=UY - UO(I  ,J  ,L)
          END IF
          IF (L.LE.LMU(IM1,J+1)) UX=UX-    UO(IM1,J+1,L)
          UT=0.5*UT
          UX=UX*BYDXP(J+1)
          UY=UY*BYDYV(J)
C****
          VT=0          ! mean v*tan on x_+ boundary for U equation
          VX=0          ! mean dv/dx on y_+ boundary for U equation
          VY=0          ! mean dv/dy on x_+ boundary for U equation
          IF (L.LE.LMV(I  ,J  )) THEN
            VT=     VO(I  ,J  ,L)*TANV(J  )
            VX=     VO(I  ,J  ,L)
            VY=     VO(I  ,J  ,L)
          END IF
          IF (J.GT.1) THEN
            IF  (L.LE.LMV(I  ,J-1)) THEN
            VT=VT + VO(I  ,J-1,L)*TANV(J-1)
            VY=VY - VO(I  ,J-1,L)
            END IF
          END IF
          IF (L.LE.LMV(IM1,J  )) VX=VX - VO(IM1,J  ,L)
          VT=0.5*VT
          VY=VY*BYDYP(J)
          VX=VX*BYDXV(J)
C**** Calculate fluxes (including FSLIP condition)
          IF (FSLIP == 1.) THEN
            IF (L.LE.LMV(I,J) .AND. L.LE.LMV(IM1,J))
     *           FUY(IM1,J)=KHV(J)*VX
            IF (L.LE.LMU(I,J) .AND. L.LE.LMU(I,J+1))
     *           FVX(I  ,J)=KHV(J)*(UY + UT)
          ELSE
            FUY(IM1,J)=KHV(J)*VX
            FVX(I  ,J)=KHV(J)*(UY + UT)
          END IF
          FUX(IM1,J)=KHP(J)*(VY + VT)
          IF (J.LT.JM-1) FVY(I,J)=KHP(J+1)*UX
          IM1=I
        END DO
      END DO

      CALL HALO_UPDATE(grid,FUY (:,grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,FVY (:,grid%j_strt_halo:grid%j_stop_halo),
     *                 FROM=SOUTH)

C**** Calculate tridiagonal matrix for first semi-implicit step (in x)
      AU=0. ; BU=0. ; CU=0. ; RU=0.
      AV=0. ; BV=0. ; CV=0. ; RV=0.

      DO J=J_0S,J_1S
        IM1=IM-1
        I=IM
        DO IP1=1,IM
          BU(I,J) = 1d0
          BV(I,J) = 1d0
          IF (L.LE.LMU(I,J)) THEN
            DTU = DT2*(DH(I,J,L)+DH(IP1,J,L))*BYMU(I,J)
            AU(I,J) =         - DTU*UXA(I,J,L)
            BU(I,J) = BU(I,J) - DTU*UXB(I,J,L)
            CU(I,J) =         - DTU*UXC(I,J,L)
            RU(I,J) = UO(I,J,L) + DTU*(UYA(I,J,L)*UO(I,J-1,L)
     *           +UYB(I,J,L)*UO(I,J,L) + UYC(I,J,L)*UO(I,J+1,L))
C**** Add Wajsowicz cross-terms to RU + second metric term
            RU(I,J) = RU(I,J) + DTU*((DYPO(J)*(FUX(IM1,J) - FUX(I,J))
     *           + DXVO(J)*FUY(I,J) - DXVO(J-1)*FUY(I,J-1))*BYDXYPO(J)
     *           - 0.5*(TANV(J-1)*FUY(I,J-1) + TANV(J)*FUY(I,J)))
          END IF
          IF (L.LE.LMV(I,J)) THEN
            DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
            AV(I,J) =         - DTV*VXA(I,J,L)
            BV(I,J) = BV(I,J) - DTV*VXB(I,J,L)
            CV(I,J) =         - DTV*VXC(I,J,L)
            RV(I,J) = VO(I,J,L) + DTV*(VYA(I,J,L)*VO(I,J-1,L)
     *           +VYB(I,J,L)*VO(I,J,L) + VYC(I,J,L)*VO(I,J+1,L))
C**** Add Wajsowicz cross-terms to RV + second metric term
            RV(I,J) = RV(I,J) + DTV*((DYVO(J)*(FVX(I,J) - FVX(IM1,J))
     *           + DXPO(J)*FVY(I,J-1) - DXPO(J+1)*FVY(I,J))*BYDXYV(J)
     *           + 0.5*(TANP(J-1)*FVY(I,J-1) + TANP(J)*FVY(I,J)))
          END IF
          IM1=I
          I=IP1
        END DO
#ifndef ODIFF_TRIDIAG_CYCLIC
C**** Minor complication due to cyclic nature of boundary condition
C**** Make properly tridiagonal by making explicit cyclic terms
        I = 1
        IF (L.LE.LMU(I,J)) THEN
          DTU = DT2*(DH(I,J,L)+DH(I+1,J,L))*BYMU(I,J)
          AU(I,J) = 0.
          RU(I,J)=RU(I,J) + DTU*UXA(I,J,L)*UO(IM,J,L)
        ENDIF
        IF (L.LE.LMV(I,J)) THEN
          AV(I,J) = 0.
          DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
          RV(I,J)=RV(I,J) + DTV*VXA(I,J,L)*VO(IM,J,L)
        ENDIF
        I = IM
        IF (L.LE.LMU(I,J)) THEN
          CU(I,J) = 0.
          DTU = DT2*(DH(I,J,L)+DH(1,J,L))*BYMU(I,J)
          RU(I,J)=RU(I,J) + DTU*UXC(I,J,L)*UO(1,J,L)
        ENDIF
        IF (L.LE.LMV(I,J)) THEN
          CV(I,J) = 0.
          DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
          RV(I,J)=RV(I,J) + DTV*VXC(I,J,L)*VO(1,J,L)
        ENDIF
#endif
      END DO
C**** At North Pole (no metric terms)
c     BU(IIP) = 1d0
c     BV(IIP) = 1d0
c     IF (L.LE.LMU(1,JM)) THEN
c     DTU = DT2*DH(1,JM,L)*BYMU(1,JM)
c       RU(IIP) = 0.
c       DO I=1,IM       ! include Wajsowicz cross-terms at North Pole
c         RU(IIP) = RU(IIP) + DTU*(UYPA(I,L)*UO(I,JM-1,L)
c    *                      - DXVO(JM-1)*FUY(I,JM-1)*BYDXYPJM)
c       END DO
c     END IF

      DO J = J_0S, J_1S
        DO I=1,IM
        if (BU(I,J).EQ.0.) write (6,*) ' Inside ODIFF: BU=0, I,J= ',i,j
        if (BV(I,J).EQ.0.) write (6,*) ' Inside ODIFF: BV=0, I,J= ',i,j
        END DO
      END DO
C**** Call tridiagonal solver
      DO J = J_0S, J_1S
#ifdef ODIFF_TRIDIAG_CYCLIC
        CALL TRIDIAG_cyclic(AU(:,J), BU(:,J), CU(:,J), RU(:,J),
     &       UO(:,J,L), IM)
        CALL TRIDIAG_cyclic(AV(:,J), BV(:,J), CV(:,J), RV(:,J),
     &       VO(:,J,L), IM)
#else
        CALL TRIDIAG(AU(:,J), BU(:,J), CU(:,J), RU(:,J), UO(:,J,L), IM)
        CALL TRIDIAG(AV(:,J), BV(:,J), CV(:,J), RV(:,J), VO(:,J,L), IM)
#endif
      END DO

      END DO ! end loop over layers

C****
C**** Now do semi implicit solution in y
C**** Recalculate rotating polar velocities from UONP and VONP
C     UO(:,JM,L) = UONP(L)*COSU(:) + VONP(L)*SINU(:)
C     VO(:,JM,L) = VONP(L)*COSI(:) - UONP(L)*SINI(:)

      CALL HALO_UPDATE(grid,UO)
      CALL HALO_UPDATE(grid,VO)

      AU3D=0. ; BU3D=0. ; CU3D=0. ; RU3D=0.; UU3D=0
      AV3D=0. ; BV3D=0. ; CV3D=0. ; RV3D=0.; UV3D=0.

      DO L=1,LMO

C**** Save (0.5*) mass reciprical for velocity points
      DO J=J_0S,J_1S
        I=IM
        DO IP1=1,IM
          IF (L.LE.LMU(I,J)) BYMU(I,J) = 1./(MO(I,J,L)+MO(IP1,J,L))
          IF (L.LE.LMV(I,J)) BYMV(I,J) = 1./(MO(I,J,L)+MO(I,J+1,L))
          I=IP1
        END DO
      END DO
      if( HAVE_NORTH_POLE ) then
        IF (L.LE.LMU(1,JM)) BYMU(1,JM) = 1./MO(1,JM,L)
      endif

C**** Calc. cross-term fluxes + second metric term (at half time step)
C**** Need dv/dy,tv,dv/dx for u equation, du/dy,tu,du/dx for v equation
      FUX=0             ! flux in U equation at the x_+ boundary
      FUY=0             ! flux in U equation at the y_+ boundary
      FVX=0             ! flux in V equation at the x_+ boundary
      FVY=0             ! flux in V equation at the y_+ boundary
      DO J=J_0,J_1S
        IM1=IM-1
        DO I=1,IM
          UT=0         ! mean u*tan on x_+ boundary for V equation
          UY=0         ! mean du/dx on y_+ boundary for V equation
          UX=0         ! mean du/dy on x_+ boundary for V equation
          IF (L.LE.LMU(I  ,J+1)) THEN
            UT=     UO(I  ,J+1,L)*TANP(J+1)
            UX=     UO(I  ,J+1,L)
            UY=     UO(I  ,J+1,L)
          END IF
          IF (L.LE.LMU(I  ,J  )) THEN
            UT=UT + UO(I  ,J  ,L)*TANP(J  )
            UY=UY - UO(I  ,J  ,L)
          END IF
          IF (L.LE.LMU(IM1,J+1)) UX=UX-    UO(IM1,J+1,L)
          UT=0.5*UT
          UX=UX*BYDXP(J+1)
          UY=UY*BYDYV(J)
C****
          VT=0         ! mean v*tan on x_+ boundary for U equation
          VX=0         ! mean dv/dx on y_+ boundary for U equation
          VY=0         ! mean dv/dy on x_+ boundary for U equation
          IF (L.LE.LMV(I  ,J  )) THEN
            VT=     VO(I  ,J  ,L)*TANV(J  )
            VX=     VO(I  ,J  ,L)
            VY=     VO(I  ,J  ,L)
          END IF
          IF (J.GT.1) THEN
            IF (L.LE.LMV(I  ,J-1)) THEN
            VT=VT + VO(I  ,J-1,L)*TANV(J-1)
            VY=VY - VO(I  ,J-1,L)
          END IF
          END IF
          IF (L.LE.LMV(IM1,J  )) VX=VX - VO(IM1,J  ,L)
          VT=0.5*VT
          VY=VY*BYDYP(J)
          VX=VX*BYDXV(J)
C**** Calculate fluxes (including FSLIP condition)
          IF (FSLIP == 1.) THEN
            IF (L.LE.LMV(I,J) .AND. L.LE.LMV(IM1,J))
     *           FUY(IM1,J)=KHV(J)* VX
            IF (L.LE.LMU(I,J) .AND. L.LE.LMU(I,J+1))
     *           FVX(I,J)=KHV(J)*(UY + UT)
          ELSE
            FUY(IM1,J)=KHV(J)* VX
            FVX(I  ,J)=KHV(J)*(UY + UT)
          END IF
          FUX(IM1,J)=KHP(J)*(VY + VT)
          IF (J.LT.JM-1) FVY(I,J)=KHP(J+1)*UX
          IM1=I
        END DO
      END DO

      CALL HALO_UPDATE(grid,FUY (:,grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)
      CALL HALO_UPDATE(grid,FVY (:,grid%j_strt_halo:grid%j_stop_halo) ,
     *                 FROM=SOUTH)

C**** Calculate tridiagonal matrix for second semi-implicit step (in y)
C**** Minor complication due to singular nature of polar box

      IM1=IM-1
      I=IM
      DO IP1=1,IM
        DO J=J_0S,J_1S
          BU3D(I,J,L) = 1d0
          BV3D(I,J,L) = 1d0
          IF (L.LE.LMU(I,J)) THEN
            DTU = DT2*(DH(I,J,L)+DH(IP1,J,L))*BYMU(I,J)
            AU3D(I,J,L) =        - DTU*UYA(I,J,L)
            BU3D(I,J,L) = BU3D(I,J,L) - DTU*UYB(I,J,L)
            IF (J.lt.JM-1) CU3D(I,J,L) =     - DTU*UYC(I,J,L)
            RU3D(I,J,L) = UO(I,J,L) + DTU*(UXA(I,J,L)*UO(IM1,J,L)
     *           +UXB(I,J,L)*UO(I,J,L) + UXC(I,J,L)*UO(IP1,J,L))
c**** Make properly tridiagonal by making explicit polar terms
!mkt  UO(1,JM,L) changed to UO(I,JM,L)
            IF (J == JM-1) RU3D(I,J,L)=
     &           RU3D(I,J,L)+DTU*UYC(I,J,L)*UO(I,JM,L)
C**** Add Wajsowicz cross-terms to RU3D + second metric term
            RU3D(I,J,L)=RU3D(I,J,L)+DTU*((DYPO(J)*(FUX(IM1,J)-FUX(I,J))
     *           + DXVO(J)*FUY(I,J) - DXVO(J-1)*FUY(I,J-1))*BYDXYPO(J)
     *           - 0.5*(TANV(J-1)*FUY(I,J-1) + TANV(J)*FUY(I,J)))
          END IF
          IF (L.LE.LMV(I,J)) THEN
            DTV = DT2*(DH(I,J,L)+DH(I,J+1,L))*BYMV(I,J)
            AV3D(I,J,L) =        - DTV*VYA(I,J,L)
            BV3D(I,J,L) = BV3D(I,J,L) - DTV*VYB(I,J,L)
            IF (J.lt.JM-1) CV3D(I,J,L) =     - DTV*VYC(I,J,L)
            RV3D(I,J,L) = VO(I,J,L) + DTV*(VXA(I,J,L)*VO(IM1,J,L)
     *           +VXB(I,J,L)*VO(I,J,L) + VXC(I,J,L)*VO(IP1,J,L))
c**** Make properly tridiagonal by making explicit polar terms
            IF (J == JM-1) RV3D(I,J,L)=
     &           RV3D(I,J,L)+DTV*VYC(I,J,L)*VO(I,JM,L)
C**** Add Wajsowicz cross-terms to RV + second metric term
            RV3D(I,J,L)=RV3D(I,J,L)+DTV*((DYVO(J)*(FVX(I,J) -FVX(IM1,J))
     *           + DXPO(J)*FVY(I,J-1) - DXPO(J+1)*FVY(I,J))*BYDXYV(J)
     *           + 0.5*(TANP(J-1)*FVY(I,J-1) + TANP(J)*FVY(I,J)))
          END IF
        END DO
        IM1=I
        I=IP1
      END DO
C**** At North Pole (do partly explicitly) no metric terms
c     BU3D(IIP) = 1d0
c     BV3D(IIP) = 1d0
c     IF (L.LE.LMU(1,JM)) THEN
c       DTU = DT2*DH(1,JM,L)*BYMU(1,JM)
c       BU3D(IIP) = BU3D(IIP) - DTU*UYPB(L)
c       RU3D(IIP) = UO(1,JM,L)
c       DO I=1,IM       ! include Wajsowicz cross-terms at North Pole
c         RU3D(IIP)= RU3D(IIP) + DTU*(UYPA(I,L)*UO(I,JM-1,L)
c    *         - DXVO(JM-1)*FUY(I,JM-1)*BYDXYPJM)
c       END DO
c     END IF
      END DO ! end loop over layers


C**** Call tridiagonal solver
      CALL TRIDIAG_new(AU3D, BU3D, CU3D, RU3D, UU3D, grid,
     &       J_LOWER=2,J_UPPER=JM-1)
      CALL TRIDIAG_new(AV3D, BV3D, CV3D, RV3D, UV3D, grid,
     &     J_LOWER=2,J_UPPER=JM-1)
      UO(:,J_0S:J_1S,:)=UU3D(:,J_0S:J_1S,:)
      VO(:,J_0S:J_1S,:)=UV3D(:,J_0S:J_1S,:)

      deallocate(AU, BU, CU, RU, UU)
      deallocate(AV, BV, CV, RV, UV)
      deallocate(AU3D, BU3D, CU3D, RU3D, UU3D)
      deallocate(AV3D, BV3D, CV3D, RV3D, UV3D)

C**** Restore unchanged UONP and VONP into prognostic locations in UO
      if(have_north_pole) then
        UO(IM  ,JM,:) = UONP(:)
        UO(IVNP,JM,:) = VONP(:)
      end if
C****
      RETURN
      END SUBROUTINE ODIFF

      SUBROUTINE TOC2SST(atmocn)
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      call get_exports_layer1
      call OG2AG_TOC2SST(atmocn)
      RETURN
      END SUBROUTINE TOC2SST

      subroutine get_exports_layer1
      use ofluxes, only : ocnatm
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getDomainBounds
      use domain_decomp_1d, only : hasSouthPole,hasNorthPole
      use domain_decomp_1d, only : halo_update,south
      use ocean, only : nbyzm,i1yzm,i2yzm
      use ocean, only : mo,g0m,s0m
      use ocean, only : uo,vo
      use ocean, only : ogeoz,ogeoz_sv
      use ocean, only : im,jm,oxyp,sinpo,sinvo
      USE OCEAN, only :
     &                oCOSI=>COSIC,oSINI=>SINIC
     &               ,IVSPO=>IVSP,IVNPO=>IVNP

#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
      use ocean, only : trmo
      use ocn_tracer_com, only: tracerlist
      use ocn_tracer_com, only: ocn_tracer_entry
#endif
#endif
      implicit none
      real*8 temgs,shcgs  ! funcs
      real*8 :: g,s
      real*8 :: awt1,awt2
      integer :: i,j,l,n,nt
      integer :: j_0,j_1,j_0s,j_1s

#if defined (TRACERS_OCEAN) && defined (TRACERS_WATER)
      type(ocn_tracer_entry), pointer :: entry
#endif

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      call getDomainBounds(grid, J_STRT_SKP = J_0S, J_STOP_SKP = J_1S)

      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        l = 1
        g = g0m(i,j,l)/(mo(i,j,l)*oxyp(i,j))
        s = s0m(i,j,l)/(mo(i,j,l)*oxyp(i,j))
        ocnatm%gtemp(i,j) = temgs(g,s)
        ocnatm%sss(i,j) = 1d3*s
        ocnatm%mlhc(i,j) = mo(i,j,1)*shcgs(g,s)
        l = 2
        g = g0m(i,j,l)/(mo(i,j,l)*oxyp(i,j))
        s = s0m(i,j,l)/(mo(i,j,l)*oxyp(i,j))
        ocnatm%gtemp2(i,j) = temgs(g,s) ! layer 2 for GCM diagnostics only
        ocnatm%ogeoza(i,j) = 0.5d0*(ogeoz(i,j)+ogeoz_sv(i,j))
#ifdef TRACERS_OCEAN
#ifdef TRACERS_WATER
C**** surface tracer concentration
        do nt=1,tracerlist%getsize()
          entry=>tracerlist%at(nt)
          if (entry%conc_from_fw) then ! define conc from fresh water
            ocnatm%gtracer(NT,I,J)=TRMO(I,J,1,NT)/
     &           (MO(I,J,1)*OXYP(I,J)-S0M(I,J,1))
          else                  ! define conc from total sea water mass
            ocnatm%gtracer(NT,I,J)=TRMO(I,J,1,NT)/
     &           (MO(I,J,1)*OXYP(I,J))
          endif
        enddo
#endif
#endif

      enddo
      enddo
      enddo

c
c ocean C-grid -> atm A-grid method requiring fewer INT_OG2AG variants:
c ocean C -> ocean A followed by ocean A -> atm A via INT_OG2AG
c
      call halo_update(grid,vo(:,:,1),from=south)

      do j=j_0s,j_1s
c area weights that would have been used by HNTRP for ocean C -> ocean A
        awt1 = (sinpo(j)-sinvo(j-1))/(sinvo(j)-sinvo(j-1))
        awt2 = 1.-awt1
        i=1
          ocnatm%uosurf(i,j) = .5*(UO(i,j,1)+UO(im,j,1))
          ocnatm%vosurf(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        do i=2,im
          ocnatm%uosurf(i,j) = .5*(UO(i,j,1)+UO(i-1,j,1))
          ocnatm%vosurf(i,j) = VO(i,j-1,1)*awt1+VO(i,j,1)*awt2
        enddo
      enddo
      if(hasSouthPole(GRID)) then
        ocnatm%uosurf(:,1) = 0.
        ocnatm%vosurf(:,1) = 0.
      endif
      if(hasNorthPole(grid)) then ! NP U,V from prognostic polar U,V
        ocnatm%uosurf(:,jm) =
     &       UO(im,jm,1)*oCOSI(:) + UO(IVNPO,jm,1)*oSINI(:)
! ocnatm%vosurf currently has no effect when atm is lat-lon
        ocnatm%vosurf(:,jm) =
     &       UO(IVNPO,jm,1)*oCOSI(:) - UO(im,jm,1)*oSINI(:)
      endif

      end subroutine get_exports_layer1

      SUBROUTINE io_oda(kunit,it,iaction,ioerr)
!@sum  io_oda dummy routine for consistency with uncoupled model
!@auth Gavin Schmidt
      RETURN
      END SUBROUTINE io_oda

      SUBROUTINE ADVSI_DIAG(atmocn,atmice)
!@sum ADVSI_DIAG dummy routine for consistency with qflux model
      use exchange_types, only : atmocn_xchng_vars,atmice_xchng_vars
      implicit none
      type(atmocn_xchng_vars) :: atmocn
      type(atmice_xchng_vars) :: atmice
      RETURN
      END SUBROUTINE ADVSI_DIAG

      SUBROUTINE GLMELT(DT)
!@sum  GLMELT adds glacial melt around Greenland and Antarctica to ocean
!@auth Sukeshi Sheth/Gavin Schmidt
!@ver  2010/07/22

      USE MODEL_COM, only : dtsrc
      USE OCEAN, only : IMO=>IM,JMO=>JM, LMM, IMAXJ,DXYPO
     *     , MO, G0M, ZE, FOCEAN, lmo, zmax_glmelt
      USE OFLUXES, only : oGMELT, oEGMELT
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : ogrid
      USE ODIAG, only : oij=>oij_loc, ij_eicb, ij_micb
#ifdef TRACERS_OCEAN
      use ocean, only : ntrtrans,motr
#ifdef TRACERS_WATER
      Use OCEAN,   Only: TRMO
      Use OFLUXES, Only: oTRGMELT
#endif
#endif

      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DT  !@var DT timestep for GLMELT call
      REAL*8 DZ
      INTEGER I,J,L,MAXGL
      integer :: j_0,j_1

      if(.not. ogrid%have_domain) return
      call getDomainBounds(ogrid, J_STRT=j_0, J_STOP=j_1)

      ! find layer index corresponding to zmax_glmelt
      do maxgl=1,lmo-1
        if(ze(maxgl+1) > zmax_glmelt) exit
      enddo

      DO L=1,MAXGL
C**** divide over depth and scale for time step
        DO J=j_0,j_1
          DO I=1,IMAXJ(J)
            If (L <= LMM(I,J) .and. oGMELT(I,J) > 0)  Then
              DZ=DT*(ZE(L)-ZE(L-1))/(DTsrc*ZE(MIN(MAXGL,LMM(I,J))))
              ! todo: remove dxypo from numerator and denominator
              MO(I,J,L) =MO(I,J,L)+(oGMELT(I,J)*dxypo(j))*DZ/
     &             (DXYPO(J)*FOCEAN(I,J))
              G0M(I,J,L)=G0M(I,J,L)+(oEGMELT(I,J)*dxypo(j))*DZ
#ifdef TRACERS_OCEAN
              if(ntrtrans.gt.1) then
                ! add mo increment to motr as well
                motr(i,j,l) = motr(i,j,l) + ogmelt(i,j)*dz
              endif
#ifdef TRACERS_WATER
              TRMO(I,J,L,:)=TRMO(I,J,L,:)+(oTRGMELT(:,I,J)*dxypo(j))*DZ
#endif
#endif
            END IF
          END DO
        END DO
      END DO
C****
      ! diagnostics
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(lmm(i,j) >0 .and. ogmelt(i,j) > 0) then
          OIJ(I,J,IJ_EICB)=OIJ(I,J,IJ_EICB)+(DT/DTsrc)*oEGMELT(I,J) !/FOCEAN(I,J)
          OIJ(I,J,IJ_MICB)=OIJ(I,J,IJ_MICB)+(DT/DTsrc)* oGMELT(I,J) !/FOCEAN(I,J)
        endif
      enddo
      enddo

      RETURN
      END SUBROUTINE GLMELT


      SUBROUTINE ADJUST_MEAN_SALT
!@sum  ADJUST_MEAN_SALT sets the global mean salinity in the ocean
!@auth Gavin Schmidt

      USE CONSTANT, only : grav
      USE OCEAN, only : oxyp,im,jm,lmo,focean,lmm,mo,s0m,sxmo
     *     ,symo,szmo,dxypo,oc_salt_mean,g0m, imaxj
      USE STRAITS, only : s0mst,sxmst,szmst,nmst,lmst,g0mst,mmst,dist
     *     ,wist

      use DOMAIN_DECOMP_1D, only: GLOBALSUM, getDomainBounds, AM_I_ROOT,
     *     broadcast
      USE OCEANR_DIM, only : ogrid

      IMPLICIT NONE
      REAL*8 :: totalSalt,totalMass

      REAL*8, DIMENSION(IM,ogrid%J_STRT_HALO:ogrid%J_STOP_HALO)::OSALT
     *     ,OMASS
      REAL*8 mean_S,frac_inc,T_ORIG,T_NEW,temgsp,shcgs,pres,g,s
      INTEGER I,J,L,N,J_0,J_1

      call getDomainBounds(ogrid, J_STRT = J_0, J_STOP = J_1)

      call conserv_OSL(OSALT)
      call conserv_OMS(OMASS)
      OSALT(:,:)=OSALT(:,:)*oXYP(:,:)
      OMASS(:,:)=OMASS(:,:)*oXYP(:,:)

      CALL GLOBALSUM(ogrid, OSALT, totalSalt, ALL=.true.)
      CALL GLOBALSUM(ogrid, OMASS, totalMass, ALL=.true.)

      if (AM_I_ROOT()) then
        mean_S=1000*totalSalt/totalMass  ! psu
        frac_inc=oc_salt_mean/mean_S
        write(6,*) "Changing ocean salinity: ",mean_S,frac_inc
      end if
      call broadcast(ogrid, frac_inc)

C**** adjust open ocean salinity
      DO J=J_0,J_1
        DO I=1,IMAXJ(J)
          PRES=0
          DO L=1,LMM(I,J)
            PRES=PRES+MO(I,J,L)*GRAV*.5
            G=G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            S=S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            T_ORIG=TEMGSP (G,S,PRES)
            S0M(I,J,L)=frac_inc*S0M(I,J,L)
            if (frac_inc.lt.1.) then
              SXMO(I,J,L)=frac_inc*SXMO(I,J,L)
              SYMO(I,J,L)=frac_inc*SYMO(I,J,L)
              SZMO(I,J,L)=frac_inc*SZMO(I,J,L)
            end if
            S=S0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            T_NEW=TEMGSP (G,S,PRES)
C**** approximately adjust enthalpy to restore temperature
            G0M(I,J,L)=G0M(I,J,L)+(T_ORIG-T_NEW)*SHCGS(G,S)*MO(I,J,L)
     *           *DXYPO(J)
            G=G0M(I,J,L)/(MO(I,J,L)*DXYPO(J))
            IF (L.lt.LMM(I,J)) PRES=PRES+MO(I,J,L)*GRAV*.5
          END DO
        END DO
      END DO

C**** adjust strait salinity
      if (am_I_root()) then
        DO N=1,NMST
          PRES=0
          DO L=1,LMST(N)
            PRES=PRES+MMST(L,N)*GRAV*.5/(DIST(N)*WIST(N))
            G=G0MST(L,N)/MMST(L,N)
            S=S0MST(L,N)/MMST(L,N)
            T_ORIG=TEMGSP (G,S,PRES)
            S0MST(L,N)=frac_inc*S0MST(L,N)
            if (frac_inc.lt.1) then
              SXMST(L,N)=frac_inc*SXMST(L,N)
              SZMST(L,N)=frac_inc*SZMST(L,N)
            end if
            S=S0MST(L,N)/MMST(L,N)
            T_NEW=TEMGSP (G,S,PRES)
C**** approximately adjust enthalpy to restore temperature
            G0MST(L,N)=G0MST(L,N)+(T_ORIG-T_NEW)*SHCGS(G,S)*MMST(L,N)
            G=G0MST(L,N)/MMST(L,N)
            IF (L.lt.LMST(N)) PRES=PRES+MMST(L,N)*GRAV*.5/(DIST(N)
     *           *WIST(N))
          END DO
        END DO
      end if
      CALL broadcast(ogrid, S0MST)
      CALL broadcast(ogrid, SXMST)
      CALL broadcast(ogrid, SZMST)

C**** Check
      call conserv_OSL(OSALT)
      OSALT(:,:)=OSALT(:,:)*oXYP(:,:)
      CALL GLOBALSUM(ogrid, OSALT, totalSalt, ALL=.true.)

      if (AM_I_ROOT()) then
        mean_S=1000*totalSalt/totalMass  ! psu
        write(6,*) "New ocean salinity: ",mean_S,oc_salt_mean
      end if

      RETURN
      END SUBROUTINE ADJUST_MEAN_SALT

      Subroutine OABFILx
C****
C**** OABFILx applies an N-th order Alternating Binomial Filter to the
C**** ocean currents UO and VO in the East-West direction
C****
      USE OCEAN, only : IM,JM,LMO,J1O, UO,VO, LMU,LMV
      USE OCEANR_DIM, only : grid=>ogrid
      USE OCEAN, only : NORDER, OABFUX, OABFVX
      Implicit None
      Real*8, Dimension(IM) :: X,Y
      Integer :: I,J,L,J_0f,J_1f,N
C****

      J_0f = max(J1O,grid%j_strt)
      J_1f = min(grid%j_stop,JM-1)

      Do 500 L=1,LMO
C****
C**** Filter U component of ocean current in east-west direction
C****
      If (OABFUX <= 0)  GoTo 200
      Do 150 J=J_0f,J_1f
      Do 110 I=1,IM
C     If (L > LMU(I,J))  UO(I,J,L) = 0
  110 X(I) = UO(I,J,L)
      Do 140 N=1,NORDER
      Do 120 I=2,IM
  120 Y(I) = X(I) - X(I-1)
      Y(1) = X(1) - X(IM)
      Do 130 I=1,IM-1
  130 If (L <= LMU(I ,J))  X(I)  = Y(I+1) - Y(I)
  140 If (L <= LMU(IM,J))  X(IM) = Y(1)   - Y(IM)
      Do 150 I=1,IM
  150 UO(I,J,L) = UO(I,J,L) - X(I)*OABFUX
C****
C**** Filter the V component of ocean current in east-west direction
C****
  200 If (OABFVX <= 0)  GoTo 500
      Do 250 J=J_0f,J_1f
      Do 210 I=1,IM
C     If (L > LMV(I,J,L))  VO(I,J,L) = 0
      Y(I) = VO(I,J,L)
  210 X(I) = 0
      Do 240 N=1,NORDER
      Do 220 I=1,IM-1
  220 If (L <= LMV(I ,J) .and. L <= LMV(I+1,J))  X(I)  = Y(I+1) - Y(I)
      If (L <= LMV(IM,J) .and. L <= LMV(1  ,J))  X(IM) = Y(1)  - Y(IM)
      Do 230 I=2,IM
  230 Y(I) = X(I) - X(I-1)
  240 Y(1) = X(1) - X(IM)
      Do 250 I=1,IM
  250 VO(I,J,L) = VO(I,J,L) - Y(I)*OABFVX
  500 Continue
      Return
      EndSubroutine OABFILx

      subroutine oabfily
C****
C**** OABFILy applies an 8-th order Alternating Binomial Filter to the
C**** ocean currents UO and VO in the North-South direction
C****
      USE OCEAN, only : IM,JM,LMO,J1O, UO,VO, LMU,LMV
      USE OCEANR_DIM, only : grid=>ogrid
      USE OCEAN, only : NORDER, by4tonv, by4tonu
      use domain_decomp_1d, only : getDomainBounds,halo_update,north,
     &                             south
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) :: x,y
      integer i,j,l,n,j_0,j_1,j_0p,j_1p
      logical :: have_north_pole

      call getDomainBounds(grid, j_strt = j_0, j_stop = j_1,
     &         have_north_pole = have_north_pole)

      j_0p = max(2,j_0)
      j_1p = min(j_1,jm-1)

      if (by4tonv <= 0 .and. by4tonu <= 0) return
C****
C**** Filter the V component of ocean current in north-south direction
C****
      do l=1,lmo
      do j=j_0p,j_1p
        x(:,j,l) = vo(:,j,l)
      enddo
      enddo
      do n=1,norder
        call halo_update(grid, x, from=south)
c compute odd n-s derivative at primary latitudes
        do l=1,lmo
        do j=j_0p,j_1p
        do i=1,im
          y(i,j,l) = x(i,j,l)-x(i,j-1,l)
        enddo
        enddo
        enddo
        if(have_north_pole) then ! pole-crossing conditions
          j = jm
          do l=1,lmo
          do i=1,im/2
            y(i,j,l) = -x(i+im/2,j-1,l)-x(i,j-1,l)
            y(i+im/2,j,l) = y(i,j,l)
          enddo
c          do i=im/2+1,im
c            y(i,j,l) = -x(i-im/2,j-1,l)-x(i,j-1,l)
c          enddo
          enddo
        endif
c compute even n-s derivative at secondary latitudes
        call halo_update(grid, y, from=north)
        do l=1,lmo
          do j=j_0p,j_1p
          do i=1,im
            if(l <= lmv(i,j)) x(i,j,l) = y(i,j+1,l)-y(i,j,l)
          enddo
          enddo
        enddo
      enddo
      do l=1,lmo
      do j=j_0p,j_1p
        vo(:,j,l) = vo(:,j,l) -x(:,j,l)*by4tonv
      enddo
      enddo

C****
C**** Filter the U component of ocean current in north-south direction
C****
      do l=1,lmo
      do j=j_0p,j_1p
      do i=1,im
        y(i,j,l) = uo(i,j,l)
        x(i,j,l) = 0.
      enddo
      enddo
      enddo
      do n=1,norder
        call halo_update(grid, y, from=north)
c compute n-s odd derivative at secondary latitudes
        do l=1,lmo
        do j=j_0p,j_1p
        do i=1,im
          if(l <= lmu(i,j) .and. l <= lmu(i,j+1))
     &         x(i,j,l) = y(i,j+1,l)-y(i,j,l)
        enddo
        enddo
        enddo
c For the filter, assume that the north pole has no ocean
        if(have_north_pole) x(:,jm-1,:) = 0.
c compute even n-s derivative at primary latitudes
        call halo_update(grid, x, from=south)
        do l=1,lmo
          do j=j_0p,j_1p
          do i=1,im
            y(i,j,l) = x(i,j,l)-x(i,j-1,l)
          enddo
          enddo
        enddo
      enddo
      do l=1,lmo
      do j=j_0p,j_1p
        uo(:,j,l) = uo(:,j,l) -y(:,j,l)*by4tonu
      enddo
      enddo

      return
      end subroutine oabfily
