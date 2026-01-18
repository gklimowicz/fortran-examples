#include "rundeck_opts.h"

#if defined(CUBED_SPHERE)
#else
#define USE_ATM_GLOBAL_ARRAYS
#endif

      SUBROUTINE init_OCEAN(iniOCEAN,istart,atmocn,dynsice)
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,broadcast
      USE SEAICE, only : osurf_tilt
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : delt1, salmin
     &  , nstep0, nstep, time0, time, itest, jtest
     &  , iocnmx, brntop, brnbot, ocnmx_factor_s, ocnmx_factor_t
     &  , diapyn, diapyc, jerlv0, thkdff
     &  , bolus_biharm_constant, bolus_laplc_constant
     &  , bolus_laplc_exponential

      USE HYCOM_ARRAYS_GLOB, only: scatter_hycom_arrays
CTNL  USE HYCOM_CPLER, only : agrid, tempr_o2a, fld_o2a
      USE HYCOM_CPLER, only : tempr_o2a, fld_o2a
      USE DOMAIN_DECOMP_ATM, only : agrid=>grid

      USE hycom_arrays_glob_renamer, only : temp_loc,saln_loc
      USE HYCOM_ATM, only : alloc_hycom_atm
      USE Dictionary_mod
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      implicit none

      logical, intent(in) :: iniOCEAN
      integer, intent(in) :: istart
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: dynsice ! not used here

      integer i,j,ia,ja

      atmocn%need_eflow_gl = .true. ! tell atm that eflow_gl is needed

CTNL  agrid => atmocn%grid ! grid used for pack_data/unpack_data in cpler

      aJ_0 = atmocn%J_0
      aJ_1 = atmocn%J_1
      aI_0 = atmocn%I_0
      aI_1 = atmocn%I_1

      aJ_0H = atmocn%J_0H
      aJ_1H = atmocn%J_1H
      aI_0H = atmocn%I_0H
      aI_1H = atmocn%I_1H

      call alloc_hycom_atm(atmocn)

#ifdef CUBED_SPHERE /* should be done for latlon atm also */
C**** Make sure to use geostrophy for ocean tilt term in ice dynamics
C**** (hycom ocean dynamics does not feel the weight of sea ice).
      osurf_tilt = 0
#endif

      call sync_param("itest", itest)
      call sync_param("jtest", jtest)
      call sync_param("iocnmx", iocnmx)
      call sync_param("brntop", brntop)
      call sync_param("brnbot", brnbot)
      call sync_param("ocnmx_factor_s", ocnmx_factor_s)
      call sync_param("ocnmx_factor_t", ocnmx_factor_t)
      call sync_param("diapyn", diapyn)
      call sync_param("diapyc", diapyc)
      call sync_param("jerlv0", jerlv0)
      call sync_param("thkdff", thkdff)
      call sync_param("bolus_biharm_constant",  bolus_biharm_constant)
      call sync_param("bolus_laplc_constant",   bolus_laplc_constant)
      call sync_param("bolus_laplc_exponential",bolus_laplc_exponential)

      if (iocnmx.ge.0.and.iocnmx.le.2 .or. iocnmx.eq.5 .or. iocnmx.eq.6)
     .                                                              then
        call inikpp
      elseif (iocnmx.eq.3 .or. iocnmx.eq.7) then
        call inigis
      else
         stop 'wrong: need to choose one ocean mixing scheme'
      endif
c
      if (AM_I_ROOT()) then ! work on global grids here
        if ((bolus_laplc_constant*bolus_laplc_exponential==1)
     . .or. (bolus_laplc_constant*bolus_biharm_constant==1)
     . .or. (bolus_laplc_exponential*bolus_biharm_constant==1)) then
          print *,' wrong bolus setting: only one can be true'
          stop 'wrong bolus setting: only one can be true'
        elseif (bolus_laplc_constant==0 .and. bolus_laplc_exponential==0
     .    .and. bolus_biharm_constant==0) then
          stop 'wrong bolus setting: one has to be true'
        end if
c

css   if (istart.eq.2 .or. nstep0.eq.0) call geopar
      call inicon
c
c --- increase temp by 2 deg
c     do 21 j=1,jj
c     do 21 l=1,isp(j)
c     do 21 i=ifp(j,l),ilp(j,l)
c       if (latij(i,j,3).lt.-65..and.lonij(i,j,3).le.5.) then !lat[-90:90],lon[0:360]
c         p(i,j,1)=0.
c         do k=1,kk
c         p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
c         temp(i,j,   k)=temp(i,j,   k)+2.*(min(p(i,j,k+1),200.*onem)-
c    .    min(p(i,j,k),200.*onem))/max(dp(i,j,k),onemu)
c         temp(i,j,kk+k)=temp(i,j,kk+k)+2.*(min(p(i,j,k+1),200.*onem)-
c    .    min(p(i,j,k),200.*onem))/max(dp(i,j,k),onemu)
c         enddo
c       endif
c21   continue
c


!!! I guess I had a good reason for commenting this out...
!!! (probably should done in inicon) IA

!!!! the following is already done in inicon
!!      if (istart.gt.2) then               !istart=2 has done this in inirfn
!!      DO J=1,JM
!!      DO I=1,IM
!!        IF (FOCEAN(I,J).gt.0.) THEN
!!          GTEMP(1,1,I,J)=asst(I,J)
!!! GTEMPR ??
!!#ifdef TRACERS_GASEXCH_ocean
!!        do nt=1,ntm
!!        GTRACER(nt,1,I,J)=atrac(I,J,nt)
!!        enddo
!!#endif
!!        END IF
!!      END DO
 !!     END DO
!!      endif

      endif ! AM_I_ROOT

      call scatter_hycom_arrays

!!! hack needed for serial inicon
      CALL broadcast(ogrid, delt1 )

      CALL broadcast(ogrid, salmin )
      CALL broadcast(ogrid, nstep0 )
      CALL broadcast(ogrid, nstep )
      CALL broadcast(ogrid, time0 )
      CALL broadcast(ogrid, time )

c moved here from inicon:
      if (nstep0.eq.0) then     ! starting from Levitus
        call fld_o2a(temp_loc(:,:,1),atmocn%work1,'sst ini')
        call tempr_o2a(temp_loc(:,:,1),atmocn%work2)
        call fld_o2a(saln_loc(:,:,1),atmocn%sss,'sss ini')
c        call fld_o2a(omlhc,mlhc,'mhl ini')
c
        do ja=aJ_0,aJ_1
          do ia=aI_0,aI_1
            if (atmocn%focean(ia,ja).gt.0.) then
              atmocn%gtemp(ia,ja)=atmocn%work1(ia,ja)
              atmocn%gtempr(ia,ja)=atmocn%work2(ia,ja)
            endif
          enddo
        enddo
c     call findmx(ip,temp,ii,ii,jj,'ini sst')
c     call findmx(ip,saln,ii,ii,jj,'ini sss')
      endif

      do ja=aJ_0,aJ_1
        do ia=aI_0,aI_1
          if (atmocn%focean(ia,ja).gt.0.) then
            if (nstep0.eq.0 .and. atmocn%sss(ia,ja).le.1.) then
              write(*,'(a,2i3,3(a,f6.2))')'chk low saln at agcm ',ia,ja,
     &             ' sss=',atmocn%sss(ia,ja),
     &             ' sst=',atmocn%gtemp(ia,ja),
     &             ' focean=',atmocn%focean(ia,ja)
              stop 'wrong sss in agcm'
            endif
          endif
        enddo
      enddo

c
      END SUBROUTINE init_OCEAN
c
      SUBROUTINE DUMMY_OCN
!@sum  DUMMY necessary entry points for non-dynamic/non-deep oceans
!@auth Gavin Schmidt
css   ENTRY ODYNAM
      !! fix later: implicit none

      ENTRY ODIFS
      ENTRY io_ocdiag
      ENTRY new_io_ocdiag
      ENTRY def_rsf_ocdiag
      ENTRY def_meta_ocdiag
      ENTRY write_meta_ocdiag
      ENTRY set_ioptrs_ocnacc_default
      ENTRY set_ioptrs_ocnacc_extended
      ENTRY init_ODEEP
      ENTRY reset_ODIAG
      ENTRY diag_OCEAN
      entry OSTRUC(QTCHNG)
      entry OCLIM(end_of_day)
      entry OSOURC (ROICE,SMSI,TGW,WTRO,OTDT,RUN0,F0DT,F2DT,RVRRUN
     *           ,RVRERUN,EVAPO,EVAPI,TFW,RUN4O,ERUN4O,RUN4I,ERUN4I
     *           ,ENRGFO,ACEFO,ACE2F,ENRGFI)
css   entry daily_OCEAN(end_of_day)
c      entry PRECIP_OC
      entry GROUND_OC
      entry io_oda(kunit,it,iaction,ioerr)
css   entry io_ocean(iu_GIC,ioread,ioerr)
css   entry CHECKO(SUBR)
c
      ENTRY ADVSI_DIAG
!!      entry alloc_ocean
c --- not calling ice dynamics
css      ENTRY DYNSI
css      ENTRY ADVSI
css      ENTRY io_icedyn
css      ENTRY io_icdiag
css      ENTRY init_icedyn
css      ENTRY reset_icdiag
css      ENTRY diag_ICEDYN
c
      entry diag_OCEAN_prep
      RETURN
      END SUBROUTINE DUMMY_OCN

      subroutine precip_oc(atmocn,iceocn)
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
      return
      end subroutine precip_oc

      subroutine diagco(m,atmocn)
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      integer :: m
      type(atmocn_xchng_vars) :: atmocn
      return
      end subroutine diagco


      subroutine def_rsf_ocean(fid)
!@sum  def_rsf_ocean defines ocean array structure in restart files
!@auth M. Kelley
!@ver  beta
      USE HYCOM_DIM, only : grid=>ogrid
      use pario, only : defvar
      USE HYCOM_SCALARS, only : nstep,time,oddev
      USE HYCOM_ARRAYS
      implicit none
      integer fid   !@var fid file id
      integer :: n
      character(len=14) :: str2d
      character(len=18) :: str3d
      character(len=20) :: str3d2
      str2d ='(idm,dist_jdm)'
      str3d ='(idm,dist_jdm,kdm)'
      str3d2='(idm,dist_jdm,kdmx2)'

#if defined(TRACERS_OceanBiology)
      call def_rsf_obio(fid)
#endif

      call defvar(grid,fid,nstep,'nstep')
      call defvar(grid,fid,time,'time')

      call defvar(grid,fid,u,'uo'//str3d2)
      call defvar(grid,fid,v,'vo'//str3d2)
      call defvar(grid,fid,dp,'dp'//str3d2)
      call defvar(grid,fid,temp,'temp'//str3d2)
      call defvar(grid,fid,saln,'saln'//str3d2)
      call defvar(grid,fid,th3d,'th3d'//str3d2)
      call defvar(grid,fid,ubavg,'ubavg(idm,dist_jdm,three)')
      call defvar(grid,fid,vbavg,'vbavg(idm,dist_jdm,three)')
      call defvar(grid,fid,pbavg,'pbavg(idm,dist_jdm,three)')
      call defvar(grid,fid,pbot,'pbot'//str2d)
      call defvar(grid,fid,psikk,'psikk'//str2d)
      call defvar(grid,fid,thkk,'thkk'//str2d)
      call defvar(grid,fid,dpmixl,'dpmixl(idm,dist_jdm,two)')
      call defvar(grid,fid,uflxav,'uflxav'//str3d)
      call defvar(grid,fid,vflxav,'vflxav'//str3d)
      call defvar(grid,fid,ufxavp,'ufxavp'//str3d)
      call defvar(grid,fid,vfxavp,'vfxavp'//str3d)
      call defvar(grid,fid,diaflx,'diaflx'//str3d)
      call defvar(grid,fid,tracer,'tracer(idm,dist_jdm,kdm,ntrcr)')
      call defvar(grid,fid,dpinit,'dpinit'//str3d)
      call defvar(grid,fid,oddev,'oddev')
      call defvar(grid,fid,uav,'uav'//str3d)
      call defvar(grid,fid,vav,'vav'//str3d)
      call defvar(grid,fid,dpuav,'dpuav'//str3d)
      call defvar(grid,fid,dpvav,'dpvav'//str3d)
      call defvar(grid,fid,dpav,'dpav'//str3d)
      call defvar(grid,fid,temav,'temav'//str3d)
      call defvar(grid,fid,salav,'salav'//str3d)
      call defvar(grid,fid,th3av,'th3av'//str3d)
      call defvar(grid,fid,pbavav,'pbavav'//str2d)
      call defvar(grid,fid,sfhtav,'sfhtav'//str2d)
      call defvar(grid,fid,eminpav,'eminpav'//str2d)
      call defvar(grid,fid,surflav,'surflav'//str2d)
      call defvar(grid,fid,salflav,'salflav'//str2d)
      call defvar(grid,fid,brineav,'brineav'//str2d)
      call defvar(grid,fid,tauxav,'tauxav'//str2d)
      call defvar(grid,fid,tauyav,'tauyav'//str2d)
      call defvar(grid,fid,dpmxav,'dpmxav'//str2d)
      call defvar(grid,fid,oiceav,'oiceav'//str2d)

c write:
c        WRITE (kunit,err=10) nstep,time
c     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
c     . ,uflxav,vflxav,diaflx,tracer,dpinit,oddev,uav,vav,dpuav,dpvav
c     . ,dpav,temav,salav,th3av,pbavav,sfhtav,eminpav
c     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
c     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi  ! agcm grid

c read: note it reads in nstep0,time0 instead of nstep,time
c            READ (kunit,err=10) HEADER,nstep0,time0
c     . ,u,v,dp,temp,saln,th3d,ubavg,vbavg,pbavg,pbot,psikk,thkk,dpmixl
c     . ,uflxav,vflxav,diaflx,tracer,dpinit,oddev,uav,vav,dpuav,dpvav
c     . ,dpav,temav,salav,th3av,pbavav,sfhtav,eminpav
c     . ,surflav,salflav,brineav,tauxav,tauyav,dpmxav,oiceav
c     . ,asst,atempr,sss,ogeoza,uosurf,vosurf,dhsi,dmsi,dssi  ! agcm grid

      return
      end subroutine def_rsf_ocean

      subroutine new_io_ocean(fid,iaction)
!@sum  new_io_ocean read/write ocean arrays from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : ioread,iowrite
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use pario, only : write_dist_data,read_dist_data,
     &     write_data,read_data
      USE HYCOM_DIM, only : grid=>ogrid
      USE HYCOM_SCALARS, only : nstep,time,nstep0,time0,baclin,oddev
      USE HYCOM_ARRAYS
      USE HYCOM_ARRAYS_GLOB, only : gather_hycom_arrays   ! for now
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: n
#if defined(TRACERS_OceanBiology)
        call new_io_obio(fid,iaction)
#endif
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_data(grid,fid,'nstep',nstep)
        call write_data(grid,fid,'time',time)
        call write_dist_data(grid,fid,'uo',u)
        call write_dist_data(grid,fid,'vo',v)
        call write_dist_data(grid,fid,'dp',dp)
        call write_dist_data(grid,fid,'temp',temp)
        call write_dist_data(grid,fid,'saln',saln)
        call write_dist_data(grid,fid,'th3d',th3d)
        call write_dist_data(grid,fid,'ubavg',ubavg)
        call write_dist_data(grid,fid,'vbavg',vbavg)
        call write_dist_data(grid,fid,'pbavg',pbavg)
        call write_dist_data(grid,fid,'pbot',pbot)
        call write_dist_data(grid,fid,'psikk',psikk)
        call write_dist_data(grid,fid,'thkk',thkk)
        call write_dist_data(grid,fid,'dpmixl',dpmixl)
        call write_dist_data(grid,fid,'uflxav',uflxav)
        call write_dist_data(grid,fid,'vflxav',vflxav)
        call write_dist_data(grid,fid,'ufxavp',uflxav)
        call write_dist_data(grid,fid,'vfxavp',vflxav)
        call write_dist_data(grid,fid,'diaflx',diaflx)
        call write_dist_data(grid,fid,'tracer',tracer)
        call write_dist_data(grid,fid,'dpinit',dpinit)
        call write_data(grid,fid,'oddev',oddev)
        call write_dist_data(grid,fid,'uav',uav)
        call write_dist_data(grid,fid,'vav',vav)
        call write_dist_data(grid,fid,'dpuav',dpuav)
        call write_dist_data(grid,fid,'dpvav',dpvav)
        call write_dist_data(grid,fid,'dpav',dpav)
        call write_dist_data(grid,fid,'temav',temav)
        call write_dist_data(grid,fid,'salav',salav)
        call write_dist_data(grid,fid,'th3av',th3av)
        call write_dist_data(grid,fid,'pbavav',pbavav)
        call write_dist_data(grid,fid,'sfhtav',sfhtav)
        call write_dist_data(grid,fid,'eminpav',eminpav)
        call write_dist_data(grid,fid,'surflav',surflav)
        call write_dist_data(grid,fid,'salflav',salflav)
        call write_dist_data(grid,fid,'brineav',brineav)
        call write_dist_data(grid,fid,'tauxav',tauxav)
        call write_dist_data(grid,fid,'tauyav',tauyav)
        call write_dist_data(grid,fid,'dpmxav',dpmxav)
        call write_dist_data(grid,fid,'oiceav',oiceav)
      case (ioread)            ! input from restart file
        call read_data(grid,fid,'nstep',nstep0,bcast_all=.true.)
        call read_data(grid,fid,'time',time0,bcast_all=.true.)
        nstep0=time0*SECONDS_PER_DAY/baclin+.0001
        write(*,'(a,i9,f9.0)')
     &       'chk ocean read at nstep/day=',nstep0,time0
        nstep=nstep0
        time=time0
        call read_dist_data(grid,fid,'uo',u)
        call read_dist_data(grid,fid,'vo',v)
        call read_dist_data(grid,fid,'dp',dp)
        call read_dist_data(grid,fid,'temp',temp)
        call read_dist_data(grid,fid,'saln',saln)
        call read_dist_data(grid,fid,'th3d',th3d)
        call read_dist_data(grid,fid,'ubavg',ubavg)
        call read_dist_data(grid,fid,'vbavg',vbavg)
        call read_dist_data(grid,fid,'pbavg',pbavg)
        call read_dist_data(grid,fid,'pbot',pbot)
        call read_dist_data(grid,fid,'psikk',psikk)
        call read_dist_data(grid,fid,'thkk',thkk)
        call read_dist_data(grid,fid,'dpmixl',dpmixl)
        call read_dist_data(grid,fid,'uflxav',uflxav)
        call read_dist_data(grid,fid,'vflxav',vflxav)
        call read_dist_data(grid,fid,'ufxavp',uflxav)
        call read_dist_data(grid,fid,'vfxavp',vflxav)
        call read_dist_data(grid,fid,'diaflx',diaflx)
        call read_dist_data(grid,fid,'tracer',tracer)
        call read_dist_data(grid,fid,'dpinit',dpinit)
        call read_data(grid,fid,'oddev',oddev,bcast_all=.true.)
        call read_dist_data(grid,fid,'uav',uav)
        call read_dist_data(grid,fid,'vav',vav)
        call read_dist_data(grid,fid,'dpuav',dpuav)
        call read_dist_data(grid,fid,'dpvav',dpvav)
        call read_dist_data(grid,fid,'dpav',dpav)
        call read_dist_data(grid,fid,'temav',temav)
        call read_dist_data(grid,fid,'salav',salav)
        call read_dist_data(grid,fid,'th3av',th3av)
        call read_dist_data(grid,fid,'pbavav',pbavav)
        call read_dist_data(grid,fid,'sfhtav',sfhtav)
        call read_dist_data(grid,fid,'eminpav',eminpav)
        call read_dist_data(grid,fid,'surflav',surflav)
        call read_dist_data(grid,fid,'salflav',salflav)
        call read_dist_data(grid,fid,'brineav',brineav)
        call read_dist_data(grid,fid,'tauxav',tauxav)
        call read_dist_data(grid,fid,'tauyav',tauyav)
        call read_dist_data(grid,fid,'dpmxav',dpmxav)
        call read_dist_data(grid,fid,'oiceav',oiceav)
c certain initialization routines still work with global
c arrays, so we have to gather
        call gather_checkpointed_hycom_arrays
      end select

      return
      end subroutine new_io_ocean

      subroutine gather_checkpointed_hycom_arrays
      ! TODO: see which global-domain arrays are
      ! really needed in non-parallelized init routines
      use hycom_arrays_glob
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      call pack_data( ogrid,  u_loc, u )
      call pack_data( ogrid,  v_loc, v )
      call pack_data( ogrid,  dp_loc, dp )
      call pack_data( ogrid,  temp_loc, temp )
      call pack_data( ogrid,  saln_loc, saln )
      call pack_data( ogrid,  th3d_loc, th3d )
      call pack_data( ogrid,  ubavg_loc, ubavg )
      call pack_data( ogrid,  vbavg_loc, vbavg )
      call pack_data( ogrid,  pbavg_loc, pbavg )
      call pack_data( ogrid,  pbot_loc, pbot )
      call pack_data( ogrid,  psikk_loc, psikk )
      call pack_data( ogrid,  thkk_loc, thkk )
      call pack_data( ogrid,  dpmixl_loc, dpmixl )
      call pack_data( ogrid,  uflxav_loc, uflxav )
      call pack_data( ogrid,  vflxav_loc, vflxav )
      call pack_data( ogrid,  ufxavp_loc, ufxavp )
      call pack_data( ogrid,  vfxavp_loc, vfxavp )
      call pack_data( ogrid,  diaflx_loc, diaflx )
      call pack_data( ogrid,  tracer_loc, tracer )
      call pack_data( ogrid,  dpinit_loc, dpinit )
      call pack_data( ogrid,  uav_loc, uav )
      call pack_data( ogrid,  vav_loc, vav )
      call pack_data( ogrid,  dpuav_loc, dpuav )
      call pack_data( ogrid,  dpvav_loc, dpvav )
      call pack_data( ogrid,  dpav_loc, dpav )
      call pack_data( ogrid,  temav_loc, temav )
      call pack_data( ogrid,  salav_loc, salav )
      call pack_data( ogrid,  th3av_loc, th3av )
CTNL  call pack_data( ogrid,  ubavav_loc, ubavav )
CTNL  call pack_data( ogrid,  vbavav_loc, vbavav )
      call pack_data( ogrid,  pbavav_loc, pbavav )
      call pack_data( ogrid,  sfhtav_loc, sfhtav )
      call pack_data( ogrid,  eminpav_loc, eminpav )
      call pack_data( ogrid,  surflav_loc, surflav )
      call pack_data( ogrid,  salflav_loc, salflav )
      call pack_data( ogrid,  brineav_loc, brineav )
      call pack_data( ogrid,  tauxav_loc, tauxav )
      call pack_data( ogrid,  tauyav_loc, tauyav )
      call pack_data( ogrid,  dpmxav_loc, dpmxav )
      call pack_data( ogrid,  oiceav_loc, oiceav )
      return
      end subroutine gather_checkpointed_hycom_arrays

c
      SUBROUTINE CHECKO(SUBR)
#ifdef USE_ATM_GLOBAL_ARRAYS
!@sum  CHECKO Checks whether Ocean are reasonable
!!      USE MODEL_COM, only : im,jm
!!      USE FLUXES, only : gtemp
!!      USE MODEL_COM, only : focean
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only: pack_block,AM_I_ROOT
c      USE HYCOM_ATM, only : gtemp,gtemp_loc
      IMPLICIT NONE
      integer i,j

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR

c      call pack_block( grid,  GTEMP_loc, GTEMP )
c      if (AM_I_ROOT()) then

c      print *,'SUBR=',SUBR
c      write(*,'(10f7.2)') ((gtemp(1,1,i,j),i=1,10),j=15,20)
c      write(*,'(a)') 'focean'
c      write(*,'(10f7.2)') ((focean(i,j),i=1,10),j=15,20)

c      endif ! AM_I_ROOT
      ! no need to sctter since nothing changed
#endif /* USE_ATM_GLOBAL_ARRAYS */
      END SUBROUTINE CHECKO
c

      SUBROUTINE daily_OCEAN(end_of_day,atmocn)
!@sum  daily_OCEAN performs the daily tasks for the ocean module
!@auth Original Development Team
C****
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: end_of_day
      type(atmocn_xchng_vars) :: atmocn
      RETURN
      END SUBROUTINE daily_OCEAN
c
css   REAL*8 FUNCTION TFREZS (SIN)
C****
C**** TFREZS calculates the freezing temperature of sea water as a
C**** function of salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: SIN (1) = salinity (kg NaCl/kg sea water), from .004 to .04
C****
C**** Output: TFREZS (C) = freezing temperature of sea water
C****
css   IMPLICIT NONE
css    REAL*8, INTENT(IN) :: SIN
css      REAL*8 :: A01 = -.0575d0, A02 = -2.154996D-4, A03 =1.710523D-3
css      REAL*8 S,S32
C****
css      S   = SIN*1.D3
css      S32 = S*DSQRT(S)
css      TFREZS = (A01 + A02*S)*S + A03*S32
css      RETURN
css      END
c
      subroutine gather_odiags
C     nothing to gather - ocean prescribed
      implicit none
      return
      end subroutine gather_odiags


      subroutine alloc_ocean
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      USE HYCOM_DIM, only : init_hycom_grid, alloc_hycom_dim
      USE HYCOM_ARRAYS, only : alloc_hycom_arrays

      USE HYCOM_DIM_GLOB, only : alloc_hycom_dim_glob
      USE HYCOM_ARRAYS_GLOB, only : alloc_hycom_arrays_glob

      USE KPRF_ARRAYS, only : alloc_kprf_arrays, alloc_kprf_arrays_local
      USE HYCOM_ATM, only : alloc_hycom_atm

      implicit none

      ! seems like this is ok place to create ocean grid since nobody
      ! uses it before this call...
      call init_hycom_grid

      !call alloc_hycom_atm

      call alloc_hycom_dim
      call alloc_hycom_arrays

      call alloc_hycom_dim_glob
      call alloc_hycom_arrays_glob

      call alloc_kprf_arrays
      call alloc_kprf_arrays_local

#ifdef TRACERS_OceanBiology
      call alloc_obio_com
#endif

      call geopar(.true.)

      !!call reset_hycom_arrays

      !if (AM_I_ROOT()) then
 !!!       call geopar
      !endif


      end subroutine alloc_ocean

#ifdef THIS_PART_IS_NOT_READY
      subroutine reset_hycom_arrays
      USE HYCOM_DIM_GLOB, only : ii,jj,kk
      USE HYCOM_ARRAYS_GLOB
      USE KPRF_ARRAYS, only : sswflx
      USE HYCOM_SCALARS, only : huge
      implicit none

      integer i,j,k,ja,jb,ia
      real :: zero = 0.


      write (*,*) 'laying out arrays in memory ...'
      do 209 j=1,jj
      do 209 i=1,ii
      p(i,j,:)=huge
      pv(i,j,1)=huge
      pbot(i,j)=huge
      ubavg(i,j,:)=huge
      vbavg(i,j,:)=huge
       utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux1(i,j)=huge
      vflux1(i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfx(i,j)=huge
      pgfy(i,j)=huge
      depthu(i,j)=huge
      depthv(i,j)=huge
      tprime(i,j)=huge
c
      srfhgt(i,j)=huge
      dpmixl(i,j,:)=huge
      oice(i,j)=huge
      taux(i,j)=huge
      tauy(i,j)=huge
      oflxa2o(i,j)=huge
      osalt(i,j)=huge
      oemnp(i,j)=huge
      ustar(i,j)=huge
      sswflx(i,j)=huge
c
      pbavav(i,j)=huge
      sfhtav(i,j)=huge
      dpmxav(i,j)=huge
      oiceav(i,j)=huge
      eminpav(i,j)=huge
      surflav(i,j)=huge
      tauxav(i,j)=huge
      tauyav(i,j)=huge
      salflav(i,j)=huge
      brineav(i,j)=huge
c
      u  (i,j,:   )=huge
      v  (i,j,:   )=huge
      uflx(i,j,:)=huge
      vflx(i,j,:)=huge
      ufxcum(i,j,:)=huge
      vfxcum(i,j,:)=huge
      dpinit(i,j,:)=huge
      dpold (i,j,:)=huge
      dp (i,j,:   )=huge
      dpu(i,j,k   )=huge
      dpv(i,j,k   )=huge
      p (i,j,:)=huge
      pu(i,j,:)=huge
      pv(i,j,:)=huge
c
      th3d(i,j,:)=huge
      thstar(i,j,:)=huge
!      do nt=1,ntrcr
        tracer(i,j,:,:)=zero
!      end do
      uav(i,j,:)=huge
      vav(i,j,:)=huge
      dpuav(i,j,:)=huge
      dpvav(i,j,:)=huge
      dpav (i,j,:)=huge
      temav(i,j,:)=huge
      salav(i,j,:)=huge
      th3av(i,j,:)=huge
      uflxav(i,j,:)=huge
      vflxav(i,j,:)=huge
      ufxavp(i,j,:)=huge
      vfxavp(i,j,:)=huge
      diaflx(i,j,:)=huge
 209  continue
c
      do 210 j=1,jj
      !!ja=mod(j-2+jj,jj)+1
      !!do 210 l=1,isq(j)
      !!do 210 i=ifq(j,l),ilq(j,l)
      do i=1,ii
      pbot(i  ,j  )=0.
      !pbot(i-1,j  )=0.
      !pbot(i  ,ja )=0.
      !pbot(i-1,ja )=0.
      p(i  ,j  ,1)=0.
      !p(i-1,j  ,1)=0.
      !p(i  ,ja ,1)=0.
      !p(i-1,ja ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,:   )=0.
      dp(i  ,j  ,k+kk)=0.
      !dp(i-1,j  ,k   )=0.
      !dp(i-1,j  ,k+kk)=0.
      !dp(i  ,ja ,k   )=0.
      !dp(i  ,ja ,k+kk)=0.
      !dp(i-1,ja ,k   )=0.
 !210  !dp(i-1,ja ,k+kk)=0.
 210  continue
c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
      do 156 j=1,jj
      do 156 i=1,ii !ifu(j,l),ilu(j,l)
      pu(i,j,:)=0.
c
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      dpu(i,j,:   )=0.
c
 156  continue
c
      do 158 j=1,jj
      !do 158 l=1,isp(j)
      do 158 i=1,ii !ifp(j,l),ilp(j,l)+1
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      dpu(:,:,:)=0.
      uflx(:,:,:)=0.
      ufxcum(:,:,:)=0.
      u(:,:,:)=0.
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      pv(:,:,:)=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
c
      dpv(:,:,:   )=0.
c
      depthv(:,:)=0.
      vtotn (:,:)=0.
      pgfy  (:,:)=0.
      vbavg(:,:,:)=0.
      vtotm (:,:)=0.
      vflux (:,:)=0.
      vflux2(:,:)=0.
      vflux3(:,:)=0.
      via(:,:)=0.
      vib(:,:)=0.
c
      dpv(:,:,:   )=0.
      vflx(:,:,:)=0.
      vfxcum(:,:,:)=0.
      v(:,:,:   )=0.
      write (*,*) '... array layout completed'
css   endif                    ! end of nstep=0

      end subroutine reset_hycom_arrays
#endif
