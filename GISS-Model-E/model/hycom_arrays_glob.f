#include "rundeck_opts.h"
c   -----------------------------------------------------------------------------
      module hycom_arrays_glob

      !USE HYCOM_DIM_GLOB

      implicit none
      private
      !public

cddd      public ip
cddd      public iu
cddd      public iv
cddd      public iq
cddd      public ifp
cddd      public ilp
cddd      public isp
cddd      public jfp
cddd      public jlp
cddd      public jsp
cddd      public ifq
cddd      public ilq
cddd      public isq
cddd      public jfq
cddd      public jlq
cddd      public jsq
cddd      public ifu
cddd      public ilu
cddd      public isu
cddd      public jfu
cddd      public jlu
cddd      public jsu
cddd      public ifv
cddd      public ilv
cddd      public isv
cddd      public jfv
cddd      public jlv
cddd      public jsv
cddd      public msk


      public scatter_hycom_arrays, gather_hycom_arrays,
     &     alloc_hycom_arrays_glob, hycom_arrays_checksum_glob

      public u
      public v
      public dp
      public dpold
      public dpu
      public dpv
      public p
      public pu
      public pv
      public latij
      public lonij
      public corio
      public potvor
      public temp
      public saln
      public th3d
      public thstar
      public wgtkap
      public psikk
      public thkk
      public dpmixl
      public srfhgt
      public montg
      public defor1
      public defor2
      public ubavg
      public vbavg
      public pbavg
      public ubrhs
      public vbrhs
      public utotm
      public vtotm
      public utotn
      public vtotn
      public uflux
      public vflux
      public uflux1
      public vflux1
      public uflux2
      public vflux2
      public uflux3
      public vflux3
      public uflx
      public vflx
      public bolusu
      public bolusv
      public uav
      public vav
      public dpuav
      public dpvav
      public temav
      public salav
      public th3av
      public dpav
      public pbavav
      public sfhtav
      public uflxav
      public vflxav
      public vfxavp
      public ufxavp
      public diaflx
      public salflav
      public brineav
      public eminpav
      public surflav
      public tauxav
      public tauyav
      public ufxcum
      public vfxcum
      public dpinit
      public dpmxav
      public oiceav
      public util1
      public util2
      public util3
      public util4
      public scpx
      public scpy
      public scux
      public scuy
      public scvx
      public scvy
      public scqx
      public scqy
      public scu2
      public scv2
      public scp2
      public scq2
      public scuxi
      public scvyi
      public scp2i
      public scq2i
      public pgfx
      public pgfy
      public gradx
      public grady
      public depthu
      public depthv
      public pvtrop
      public depths
      public drag
      public glue
      public zone
      public dampu
      public dampv
      public uja
      public ujb
      public via
      public vib
      public pbot
      public tracer
      public diadff
      public tprime
      public sgain
      public surflx
      public salflx
      public sflxcum
      public hflxcum
      public odhsi
      public odmsi
      public omlhc
      public dmfz
      public taux
      public tauy
      public oemnp
      public oicemlt
      public oflxa2o
      public oice
      public ustar
      public ustarb
      public osalt
      public freshw
      public diafor
      public klist
      public ijlist
      public diag1,diag2,diag3,diag4


!!      c o m m o n
      real, allocatable ::
     . u(:,:,:),v(:,:,:)          ! velocity components
     .,dp(:,:,:),dpold(:,:,:)     ! layer thickness
     .,dpu(:,:,:),dpv(:,:,:)      ! layer thickness at u,v points
     .,p(:,:,:)                   ! interface pressure
     .,pu(:,:,:),pv(:,:,:)        ! interface pres. at u,v points
     .,latij(:,:,:),lonij(:,:,:)  ! latitude/longitude
     .,corio(:,:)                 ! coriolis parameter
     .,potvor(:,:)                ! potential vorticity
     .,temp(:,:,:)                ! temperature
     .,saln(:,:,:)                ! salinity
     .,th3d(:,:,:)                ! potential density
     .,thstar(:,:,:)              ! virtual potential density
     .,wgtkap(:,:)                ! scale factor
     .,psikk(:,:)                 ! init.montg.pot. in bottom layer
     .,thkk(:,:)                  ! init.thstar in bottom layer
     .,dpmixl(:,:,:)              ! Kraus-Turner mixed layer depth
     .,srfhgt(:,:)                ! sea surface height
c
!!      real u,v,dp,dpold,dpu,dpv,p,pu,pv,latij,lonij,corio,potvor,
!!     .     temp,saln,th3d,thstar,psikk,thkk,dpmixl,srfhgt
c
!!      c o m m o n
      real, allocatable ::
     . montg(:,:,:)			! montgomery potential
     .,defor1(:,:),defor2(:,:)		! deformation components
     .,ubavg(:,:,:),vbavg(:,:,:)	! barotropic velocity
     .,pbavg(:,:,:)			! barotropic pressure
     .,ubrhs(:,:),vbrhs(:,:)		! rhs of barotropic u,v eqns.
     .,utotm(:,:),vtotm(:,:)		! total (barotrop.+baroclin.)..
     .,utotn(:,:),vtotn(:,:)		! ..velocities at 2 time levels
     .,uflux(:,:),vflux(:,:)		! horizontal mass fluxes
     .,uflux1(:,:),vflux1(:,:)		! more mass fluxes
     .,uflux2(:,:),vflux2(:,:)		! more mass fluxes
     .,uflux3(:,:),vflux3(:,:)		! more mass fluxes
     .,uflx(:,:,:),vflx(:,:,:)		! more mass fluxes
     .,bolusu(:,:,:),bolusv(:,:,:)	! thickness (bolus) fluxes
c
!!      real montg,defor1,defor2,ubavg,vbavg,pbavg,ubrhs,vbrhs,utotm,
!!     .     vtotm,utotn,vtotn,uflux,vflux,uflux1,vflux1,uflux2,vflux2,
!!     .     uflux3,vflux3,uflx,vflx,bolusu,bolusv
c
!!      c o m m o n /timav/                ! fields needed for time-averaging
      real, allocatable ::
     .   uav(:,:,:),  vav(:,:,:)
     .,dpuav(:,:,:),dpvav(:,:,:)
     .,temav(:,:,:),salav(:,:,:)
     .,th3av(:,:,:), dpav(:,:,:)
     .,pbavav(:,:),sfhtav(:,:)
     .,uflxav(:,:,:),vflxav(:,:,:)
     .,ufxavp(:,:,:),vfxavp(:,:,:)
     .,diaflx(:,:,:)                    ! time integral of diapyc.flux
     .,salflav(:,:),brineav(:,:),eminpav(:,:),surflav(:,:)
     .,tauxav(:,:),tauyav(:,:)
     .,ufxcum(:,:,:),vfxcum(:,:,:),dpinit(:,:,:)
     .,dpmxav(:,:),oiceav(:,:)
c
!!      real uav,vav,dpuav,dpvav,temav,salav,th3av,dpav
!!     .    ,pbavav,sfhtav,uflxav,vflxav,diaflx,salflav,brineav,eminpav
!!     .    ,surflav,ufxcum,vfxcum,dpinit
!!     .    ,dpmxav,oiceav
c
!!      c o m m o n
      real, allocatable ::
     . util1(:,:),util2(:,:)		! arrays for temporary storage
     .,util3(:,:),util4(:,:)		! arrays for temporary storage
c
     .,scpx(:,:),scpy(:,:)		! mesh size at p pts in x,y dir.
     .,scux(:,:),scuy(:,:)		! mesh size at u pts in x,y dir.
     .,scvx(:,:),scvy(:,:)		! mesh size at v pts in x,y dir.
     .,scqx(:,:),scqy(:,:)		! mesh size at q pts in x,y dir.
     .,scu2(:,:),scv2(:,:)		! grid box size at u,v pts
     .,scp2(:,:),scq2(:,:)		! grid box size at p,q pts
     .,scuxi(:,:),scvyi(:,:)		! inverses of scux,scvy
     .,scp2i(:,:),scq2i(:,:)		! inverses of scp2,scq2
c
     .,pgfx(:,:),pgfy(:,:)		! horiz. presssure gradient
     .,gradx(:,:),grady(:,:)		! horiz. presssure gradient
     .,depthu(:,:),depthv(:,:)		! bottom pres. at u,v points
     .,pvtrop(:,:)			! pot.vort. of barotropic flow
     .,depths(:,:)			! water depth
     .,drag(:,:)			! bottom drag
     .,glue(:,:)			! regional viscosity enhancement
     .,zone(:,:,:)			! zone of arctic/main/antarctic
     .,dampu(:,:),dampv(:,:)		! coastal wave damping coeff.
c
!!      real util1,util2,util3,util4,scpx,scpy,scux,scuy,scvx,scvy,
!!     .     scqx,scqy,scu2,scv2,scp2,scq2,scuxi,scvyi,scp2i,scq2i,
!!     .     pgfx,pgfy,gradx,grady,depthu,depthv,pvtrop,depths,drag,
!!     .     glue,dampu,dampv
c
!!      c o m m o n
      real, allocatable ::
     . uja(:,:),ujb(:,:)		! velocities at lateral
     .,via(:,:),vib(:,:)		!          neighbor points
     .,pbot(:,:)			! bottom pressure at t=0
     .,tracer(:,:,:,:)			! tracer
     .,diadff(:,:,:)			!
     .,tprime(:,:)			! temp.change due to surflx
     .,sgain(:,:)			! salin.changes from diapyc.mix.
     .,surflx(:,:)			! surface thermal energy flux
     .,salflx(:,:)			! surface salinity flux
     .,sflxcum(:,:)			! accumulated saltflux
     .,hflxcum(:,:)			! accumulated heatflux
c    .,thkice(:,:)			! grid-cell avg. ice thknss (cm)
c    .,covice(:,:)			! ice coverage (rel.units)
c    .,temice(:,:)			! ice surf.temp.
     .,odhsi(:,:)			! heat borrowed from frozen
     .,odmsi(:,:)			! mass corresponding to odhsi
     .,omlhc(:,:)
     .,dmfz(:,:)            ! ice mass due to freezing
c
      integer, allocatable, dimension (:,:) ::
     .  klist               ! k-index of layer below mixl'r
     . ,ijlist              ! global ij index
c
!!    common/int1/klist

c ---  s w i t c h e s    (if set to .true., then...)
c --- diagno      output model fields and diagnostic messages
c --- thermo      use thermodynamic forcing functions
c --- windf       include wind stress in forcing functions
c --- relax       activate lateral boundary nudging
c --- trcout      advect tracer and save results in history/restart file
c --- dotrcr      perform column physics operations on tracer array(s)
c
!!      logical diagno,thermo,windf,relax,trcout,dotrcr
!!      common/swtchs/diagno,thermo,windf,relax,trcout,dotrcr
c
!!      c o m m o n  /frcing/                   !  monthly forcing fields
      real, allocatable ::
     . taux(:,:)                          !  wind stress in x direction
     .,tauy(:,:)                          !  wind stress in y direction
c    .,wndspd(:,:,:)                      !  wind speed (tke source)
c    .,airtmp(:,:,:)                      !  pseudo air temperature
c    .,vapmix(:,:,:)                      !  atmosph. vapor mixing ratio
c    .,oprec(:,:)                         !  precipitation
c    .,oevap(:,:)                         !  evaportation
     .,oemnp(:,:)                         !  e - p
     .,oicemlt(:,:)                       !  ice melt on o-grid
     .,oflxa2o(:,:),oice(:,:)
     .,ustar(:,:)                         ! surface friction velocity
     .,ustarb(:,:)                        ! bottom friction velocity
     .,osalt(:,:)                         ! saltflux from SI(kg/m*m)
     .,freshw(:,:)                        !  river & glacier runoff
     .,diafor(:,:)                        !  imposed diapycnal forcing
     .,diag1(:,:),diag2(:,:),diag3(:,:),diag4(:,:)
c
      contains

      subroutine scatter_hycom_arrays
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: UNPACK_DATA

      !return

      call unpack_data( ogrid,  u, u_loc )
      call unpack_data( ogrid,  v, v_loc )
      call unpack_data( ogrid,  dp, dp_loc )
      call unpack_data( ogrid,  dpold, dpold_loc )
      call unpack_data( ogrid,  dpu, dpu_loc )
      call unpack_data( ogrid,  dpv, dpv_loc )
      call unpack_data( ogrid,  p, p_loc )
      call unpack_data( ogrid,  pu, pu_loc )
      call unpack_data( ogrid,  pv, pv_loc )
      call unpack_data( ogrid,  latij, latij_loc )
      call unpack_data( ogrid,  lonij, lonij_loc )
      call unpack_data( ogrid,  corio, corio_loc )
      call unpack_data( ogrid,  potvor, potvor_loc )
      call unpack_data( ogrid,  temp, temp_loc )
      call unpack_data( ogrid,  saln, saln_loc )
      call unpack_data( ogrid,  th3d, th3d_loc )
      call unpack_data( ogrid,  thstar, thstar_loc )
      call unpack_data( ogrid,  wgtkap, wgtkap_loc )
      call unpack_data( ogrid,  psikk, psikk_loc )
      call unpack_data( ogrid,  thkk, thkk_loc )
      call unpack_data( ogrid,  dpmixl, dpmixl_loc )
      call unpack_data( ogrid,  srfhgt, srfhgt_loc )
      call unpack_data( ogrid,  montg, montg_loc )
      call unpack_data( ogrid,  defor1, defor1_loc )
      call unpack_data( ogrid,  defor2, defor2_loc )
      call unpack_data( ogrid,  ubavg, ubavg_loc )
      call unpack_data( ogrid,  vbavg, vbavg_loc )
      call unpack_data( ogrid,  pbavg, pbavg_loc )
      call unpack_data( ogrid,  ubrhs, ubrhs_loc )
      call unpack_data( ogrid,  vbrhs, vbrhs_loc )
      call unpack_data( ogrid,  utotm, utotm_loc )
      call unpack_data( ogrid,  vtotm, vtotm_loc )
      call unpack_data( ogrid,  utotn, utotn_loc )
      call unpack_data( ogrid,  vtotn, vtotn_loc )
      call unpack_data( ogrid,  uflux, uflux_loc )
      call unpack_data( ogrid,  vflux, vflux_loc )
      call unpack_data( ogrid,  uflux1, uflux1_loc )
      call unpack_data( ogrid,  vflux1, vflux1_loc )
      call unpack_data( ogrid,  uflux2, uflux2_loc )
      call unpack_data( ogrid,  vflux2, vflux2_loc )
      call unpack_data( ogrid,  uflux3, uflux3_loc )
      call unpack_data( ogrid,  vflux3, vflux3_loc )
      call unpack_data( ogrid,  uflx, uflx_loc )
      call unpack_data( ogrid,  vflx, vflx_loc )
      call unpack_data( ogrid,  bolusu, bolusu_loc )
      call unpack_data( ogrid,  bolusv, bolusv_loc )
      call unpack_data( ogrid,  uav, uav_loc )
      call unpack_data( ogrid,  vav, vav_loc )
      call unpack_data( ogrid,  dpuav, dpuav_loc )
      call unpack_data( ogrid,  dpvav, dpvav_loc )
      call unpack_data( ogrid,  temav, temav_loc )
      call unpack_data( ogrid,  salav, salav_loc )
      call unpack_data( ogrid,  th3av, th3av_loc )
      call unpack_data( ogrid,  dpav, dpav_loc )
      call unpack_data( ogrid,  pbavav, pbavav_loc )
      call unpack_data( ogrid,  sfhtav, sfhtav_loc )
      call unpack_data( ogrid,  uflxav, uflxav_loc )
      call unpack_data( ogrid,  vflxav, vflxav_loc )
      call unpack_data( ogrid,  ufxavp, ufxavp_loc )
      call unpack_data( ogrid,  vfxavp, vfxavp_loc )
      call unpack_data( ogrid,  diaflx, diaflx_loc )
      call unpack_data( ogrid,  salflav, salflav_loc )
      call unpack_data( ogrid,  brineav, brineav_loc )
      call unpack_data( ogrid,  eminpav, eminpav_loc )
      call unpack_data( ogrid,  surflav, surflav_loc )
      call unpack_data( ogrid,  tauxav, tauxav_loc )
      call unpack_data( ogrid,  tauyav, tauyav_loc )
      call unpack_data( ogrid,  ufxcum, ufxcum_loc )
      call unpack_data( ogrid,  vfxcum, vfxcum_loc )
      call unpack_data( ogrid,  dpinit, dpinit_loc )
      call unpack_data( ogrid,  dpmxav, dpmxav_loc )
      call unpack_data( ogrid,  oiceav, oiceav_loc )
      call unpack_data( ogrid,  util1, util1_loc )
      call unpack_data( ogrid,  util2, util2_loc )
      call unpack_data( ogrid,  util3, util3_loc )
      call unpack_data( ogrid,  util4, util4_loc )
      call unpack_data( ogrid,  scpx, scpx_loc )
      call unpack_data( ogrid,  scpy, scpy_loc )
      call unpack_data( ogrid,  scux, scux_loc )
      call unpack_data( ogrid,  scuy, scuy_loc )
      call unpack_data( ogrid,  scvx, scvx_loc )
      call unpack_data( ogrid,  scvy, scvy_loc )
      call unpack_data( ogrid,  scqx, scqx_loc )
      call unpack_data( ogrid,  scqy, scqy_loc )
      call unpack_data( ogrid,  scu2, scu2_loc )
      call unpack_data( ogrid,  scv2, scv2_loc )
      call unpack_data( ogrid,  scp2, scp2_loc )
      call unpack_data( ogrid,  scq2, scq2_loc )
      call unpack_data( ogrid,  scuxi, scuxi_loc )
      call unpack_data( ogrid,  scvyi, scvyi_loc )
      call unpack_data( ogrid,  scp2i, scp2i_loc )
      call unpack_data( ogrid,  scq2i, scq2i_loc )
      call unpack_data( ogrid,  pgfx, pgfx_loc )
      call unpack_data( ogrid,  pgfy, pgfy_loc )
      call unpack_data( ogrid,  gradx, gradx_loc )
      call unpack_data( ogrid,  grady, grady_loc )
      call unpack_data( ogrid,  depthu, depthu_loc )
      call unpack_data( ogrid,  depthv, depthv_loc )
      call unpack_data( ogrid,  pvtrop, pvtrop_loc )
      call unpack_data( ogrid,  depths, depths_loc )
      call unpack_data( ogrid,  drag, drag_loc )
      call unpack_data( ogrid,  glue, glue_loc )
      call unpack_data( ogrid,  zone, zone_loc )
      call unpack_data( ogrid,  dampu, dampu_loc )
      call unpack_data( ogrid,  dampv, dampv_loc )
      call unpack_data( ogrid,  uja, uja_loc )
      call unpack_data( ogrid,  ujb, ujb_loc )
      call unpack_data( ogrid,  via, via_loc )
      call unpack_data( ogrid,  vib, vib_loc )
      call unpack_data( ogrid,  pbot, pbot_loc )
      call unpack_data( ogrid,  tracer, tracer_loc )
      call unpack_data( ogrid,  diadff, diadff_loc )
      call unpack_data( ogrid,  tprime, tprime_loc )
      !!!call unpack_data( ogrid,  sgain, sgain_loc )
      call unpack_data( ogrid,  surflx, surflx_loc )
      call unpack_data( ogrid,  salflx, salflx_loc )
      call unpack_data( ogrid,  sflxcum, sflxcum_loc )
      call unpack_data( ogrid,  hflxcum, hflxcum_loc )
      call unpack_data( ogrid,  odhsi, odhsi_loc )
      call unpack_data( ogrid,  odmsi, odmsi_loc )
      call unpack_data( ogrid,  omlhc, omlhc_loc )
      call unpack_data( ogrid,  dmfz, dmfz_loc )
      call unpack_data( ogrid,  taux, taux_loc )
      call unpack_data( ogrid,  tauy, tauy_loc )
      call unpack_data( ogrid,  oemnp, oemnp_loc )
      call unpack_data( ogrid,  oicemlt, oicemlt_loc )
      call unpack_data( ogrid,  oflxa2o, oflxa2o_loc )
      call unpack_data( ogrid,  oice, oice_loc )
      call unpack_data( ogrid,  ustar, ustar_loc )
      call unpack_data( ogrid,  ustarb, ustarb_loc )
      call unpack_data( ogrid,  osalt, osalt_loc )
      call unpack_data( ogrid,  freshw, freshw_loc )
      call unpack_data( ogrid,  diafor, diafor_loc )
      call unpack_data( ogrid,  klist, klist_loc )
      call unpack_data( ogrid,  ijlist,ijlist_loc )
      call unpack_data( ogrid,  diag1,diag1_loc )
      call unpack_data( ogrid,  diag2,diag2_loc )
      call unpack_data( ogrid,  diag3,diag3_loc )
      call unpack_data( ogrid,  diag4,diag4_loc )

      end subroutine scatter_hycom_arrays


      subroutine gather_hycom_arrays
      use hycom_arrays_glob_renamer
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA

      !return

      call pack_data( ogrid,  u_loc, u )
      call pack_data( ogrid,  v_loc, v )
      call pack_data( ogrid,  dp_loc, dp )
      call pack_data( ogrid,  dpold_loc, dpold )
      call pack_data( ogrid,  dpu_loc, dpu )
      call pack_data( ogrid,  dpv_loc, dpv )
      call pack_data( ogrid,  p_loc, p )
      call pack_data( ogrid,  pu_loc, pu )
      call pack_data( ogrid,  pv_loc, pv )
      call pack_data( ogrid,  latij_loc, latij )
      call pack_data( ogrid,  lonij_loc, lonij )
      call pack_data( ogrid,  corio_loc, corio )
      call pack_data( ogrid,  potvor_loc, potvor )
      call pack_data( ogrid,  temp_loc, temp )
      call pack_data( ogrid,  saln_loc, saln )
      call pack_data( ogrid,  th3d_loc, th3d )
      call pack_data( ogrid,  thstar_loc, thstar )
      call pack_data( ogrid,  wgtkap_loc, wgtkap )
      call pack_data( ogrid,  psikk_loc, psikk )
      call pack_data( ogrid,  thkk_loc, thkk )
      call pack_data( ogrid,  dpmixl_loc, dpmixl )
      call pack_data( ogrid,  srfhgt_loc, srfhgt )
      call pack_data( ogrid,  montg_loc, montg )
      call pack_data( ogrid,  defor1_loc, defor1 )
      call pack_data( ogrid,  defor2_loc, defor2 )
      call pack_data( ogrid,  ubavg_loc, ubavg )
      call pack_data( ogrid,  vbavg_loc, vbavg )
      call pack_data( ogrid,  pbavg_loc, pbavg )
      call pack_data( ogrid,  ubrhs_loc, ubrhs )
      call pack_data( ogrid,  vbrhs_loc, vbrhs )
      call pack_data( ogrid,  utotm_loc, utotm )
      call pack_data( ogrid,  vtotm_loc, vtotm )
      call pack_data( ogrid,  utotn_loc, utotn )
      call pack_data( ogrid,  vtotn_loc, vtotn )
      call pack_data( ogrid,  uflux_loc, uflux )
      call pack_data( ogrid,  vflux_loc, vflux )
      call pack_data( ogrid,  uflux1_loc, uflux1 )
      call pack_data( ogrid,  vflux1_loc, vflux1 )
      call pack_data( ogrid,  uflux2_loc, uflux2 )
      call pack_data( ogrid,  vflux2_loc, vflux2 )
      call pack_data( ogrid,  uflux3_loc, uflux3 )
      call pack_data( ogrid,  vflux3_loc, vflux3 )
      call pack_data( ogrid,  uflx_loc, uflx )
      call pack_data( ogrid,  vflx_loc, vflx )
      call pack_data( ogrid,  bolusu_loc, bolusu )
      call pack_data( ogrid,  bolusv_loc, bolusv )
      call pack_data( ogrid,  uav_loc, uav )
      call pack_data( ogrid,  vav_loc, vav )
      call pack_data( ogrid,  dpuav_loc, dpuav )
      call pack_data( ogrid,  dpvav_loc, dpvav )
      call pack_data( ogrid,  temav_loc, temav )
      call pack_data( ogrid,  salav_loc, salav )
      call pack_data( ogrid,  th3av_loc, th3av )
      call pack_data( ogrid,  dpav_loc, dpav )
      call pack_data( ogrid,  pbavav_loc, pbavav )
      call pack_data( ogrid,  sfhtav_loc, sfhtav )
      call pack_data( ogrid,  uflxav_loc, uflxav )
      call pack_data( ogrid,  vflxav_loc, vflxav )
      call pack_data( ogrid,  ufxavp_loc, ufxavp )
      call pack_data( ogrid,  vfxavp_loc, vfxavp )
      call pack_data( ogrid,  diaflx_loc, diaflx )
      call pack_data( ogrid,  salflav_loc, salflav )
      call pack_data( ogrid,  brineav_loc, brineav )
      call pack_data( ogrid,  eminpav_loc, eminpav )
      call pack_data( ogrid,  surflav_loc, surflav )
      call pack_data( ogrid,  tauxav_loc, tauxav )
      call pack_data( ogrid,  tauyav_loc, tauyav )
      call pack_data( ogrid,  ufxcum_loc, ufxcum )
      call pack_data( ogrid,  vfxcum_loc, vfxcum )
      call pack_data( ogrid,  dpinit_loc, dpinit )
      call pack_data( ogrid,  dpmxav_loc, dpmxav )
      call pack_data( ogrid,  oiceav_loc, oiceav )
      call pack_data( ogrid,  util1_loc, util1 )
      call pack_data( ogrid,  util2_loc, util2 )
      call pack_data( ogrid,  util3_loc, util3 )
      call pack_data( ogrid,  util4_loc, util4 )
      call pack_data( ogrid,  scpx_loc, scpx )
      call pack_data( ogrid,  scpy_loc, scpy )
      call pack_data( ogrid,  scux_loc, scux )
      call pack_data( ogrid,  scuy_loc, scuy )
      call pack_data( ogrid,  scvx_loc, scvx )
      call pack_data( ogrid,  scvy_loc, scvy )
      call pack_data( ogrid,  scqx_loc, scqx )
      call pack_data( ogrid,  scqy_loc, scqy )
      call pack_data( ogrid,  scu2_loc, scu2 )
      call pack_data( ogrid,  scv2_loc, scv2 )
      call pack_data( ogrid,  scp2_loc, scp2 )
      call pack_data( ogrid,  scq2_loc, scq2 )
      call pack_data( ogrid,  scuxi_loc, scuxi )
      call pack_data( ogrid,  scvyi_loc, scvyi )
      call pack_data( ogrid,  scp2i_loc, scp2i )
      call pack_data( ogrid,  scq2i_loc, scq2i )
      call pack_data( ogrid,  pgfx_loc, pgfx )
      call pack_data( ogrid,  pgfy_loc, pgfy )
      call pack_data( ogrid,  gradx_loc, gradx )
      call pack_data( ogrid,  grady_loc, grady )
      call pack_data( ogrid,  depthu_loc, depthu )
      call pack_data( ogrid,  depthv_loc, depthv )
      call pack_data( ogrid,  pvtrop_loc, pvtrop )
      call pack_data( ogrid,  depths_loc, depths )
      call pack_data( ogrid,  drag_loc, drag )
      call pack_data( ogrid,  glue_loc, glue )
      call pack_data( ogrid,  zone_loc, zone )
      call pack_data( ogrid,  dampu_loc, dampu )
      call pack_data( ogrid,  dampv_loc, dampv )
      call pack_data( ogrid,  uja_loc, uja )
      call pack_data( ogrid,  ujb_loc, ujb )
      call pack_data( ogrid,  via_loc, via )
      call pack_data( ogrid,  vib_loc, vib )
      call pack_data( ogrid,  pbot_loc, pbot )
      call pack_data( ogrid,  tracer_loc, tracer )
      call pack_data( ogrid,  diadff_loc, diadff)
      call pack_data( ogrid,  tprime_loc, tprime )
      !!!call pack_data( ogrid,  sgain_loc, sgain )
      call pack_data( ogrid,  surflx_loc, surflx )
      call pack_data( ogrid,  salflx_loc, salflx )
      call pack_data( ogrid,  sflxcum_loc, sflxcum )
      call pack_data( ogrid,  hflxcum_loc, hflxcum )
      call pack_data( ogrid,  odhsi_loc, odhsi )
      call pack_data( ogrid,  odmsi_loc, odmsi )
      call pack_data( ogrid,  omlhc_loc, omlhc )
      call pack_data( ogrid,  dmfz_loc, dmfz )
      call pack_data( ogrid,  taux_loc, taux )
      call pack_data( ogrid,  tauy_loc, tauy )
      call pack_data( ogrid,  oemnp_loc, oemnp )
      call pack_data( ogrid,  oicemlt_loc, oicemlt )
      call pack_data( ogrid,  oflxa2o_loc, oflxa2o )
      call pack_data( ogrid,  oice_loc, oice )
      call pack_data( ogrid,  ustar_loc, ustar )
      call pack_data( ogrid,  ustarb_loc, ustarb )
      call pack_data( ogrid,  osalt_loc, osalt )
      call pack_data( ogrid,  freshw_loc, freshw )
      call pack_data( ogrid,  diafor_loc, diafor )
      call pack_data( ogrid,  klist_loc, klist )
      call pack_data( ogrid,  ijlist_loc, ijlist )
      call pack_data( ogrid,  diag1_loc, diag1 )
      call pack_data( ogrid,  diag2_loc, diag2 )
      call pack_data( ogrid,  diag3_loc, diag3 )
      call pack_data( ogrid,  diag4_loc, diag4 )

      end subroutine gather_hycom_arrays


      subroutine alloc_hycom_arrays_glob
      USE HYCOM_DIM, only :
     &     idm_full=>idm,jdm_full=>jdm,kdm_full=>kdm,ntrcr_full=>ntrcr
      use domain_decomp_1d, only : am_i_root
      integer :: idm,jdm,kdm,ntrcr

      if(am_i_root()) then
        idm = idm_full
        jdm = jdm_full
        kdm = kdm_full
        ntrcr = ntrcr_full
      else
        idm = 1
        jdm = 1
        kdm = 1
        ntrcr = 1
      endif

      allocate(
     &     depths(idm_full,jdm_full)
     &     )
      depths = 0

      allocate(
     . u(idm,jdm,2*kdm),v(idm,jdm,2*kdm)
     .,dp(idm,jdm,2*kdm),dpold(idm,jdm,kdm)
     .,dpu(idm,jdm,2*kdm),dpv(idm,jdm,2*kdm)
     .,p(idm,jdm,kdm+1)
     .,pu(idm,jdm,kdm+1),pv(idm,jdm,kdm+1)
     .,latij(idm,jdm,4),lonij(idm,jdm,4)
     .,corio(idm,jdm)
     .,potvor(idm,jdm)
     .,temp(idm,jdm,2*kdm)
     .,saln(idm,jdm,2*kdm)
     .,th3d(idm,jdm,2*kdm)
     .,thstar(idm,jdm,2*kdm)
     .,wgtkap(idm,jdm)
     .,psikk(idm,jdm)
     .,thkk(idm,jdm)
     .,dpmixl(idm,jdm,2)
     .,srfhgt(idm,jdm) )
c
      allocate(
     . montg(idm,jdm,kdm)
     .,defor1(idm,jdm),defor2(idm,jdm)
     .,ubavg(idm,jdm,3),vbavg(idm,jdm,3)
     .,pbavg(idm,jdm,3)
     .,ubrhs(idm,jdm),vbrhs(idm,jdm)
     .,utotm(idm,jdm),vtotm(idm,jdm)
     .,utotn(idm,jdm),vtotn(idm,jdm)
     .,uflux(idm,jdm),vflux(idm,jdm)
     .,uflux1(idm,jdm),vflux1(idm,jdm)
     .,uflux2(idm,jdm),vflux2(idm,jdm)
     .,uflux3(idm,jdm),vflux3(idm,jdm)
     .,uflx(idm,jdm,kdm),vflx(idm,jdm,kdm)
     .,bolusu(idm,jdm,kdm),bolusv(idm,jdm,kdm) )
c
      allocate(
     .   uav(idm,jdm,kdm),  vav(idm,jdm,kdm)
     .,dpuav(idm,jdm,kdm),dpvav(idm,jdm,kdm)
     .,temav(idm,jdm,kdm),salav(idm,jdm,kdm)
     .,th3av(idm,jdm,kdm), dpav(idm,jdm,kdm)
     .,pbavav(idm,jdm),sfhtav(idm,jdm)
     .,uflxav(idm,jdm,kdm),vflxav(idm,jdm,kdm)
     .,ufxavp(idm,jdm,kdm),vfxavp(idm,jdm,kdm)
     .,diaflx(idm,jdm,kdm)
     .,salflav(idm,jdm),brineav(idm,jdm),eminpav(idm,jdm)
     .,surflav(idm,jdm),tauxav(idm,jdm),tauyav(idm,jdm)
     .,ufxcum(idm,jdm,kdm),vfxcum(idm,jdm,kdm),dpinit(idm,jdm,kdm)
     .,dpmxav(idm,jdm),oiceav(idm,jdm)
     .)
c
      allocate(
     . util1(idm,jdm),util2(idm,jdm)
     .,util3(idm,jdm),util4(idm,jdm)
c
     .,scpx(idm,jdm),scpy(idm,jdm)
     .,scux(idm,jdm),scuy(idm,jdm)
     .,scvx(idm,jdm),scvy(idm,jdm)
     .,scqx(idm,jdm),scqy(idm,jdm)
     .,scu2(idm,jdm),scv2(idm,jdm)
     .,scp2(idm,jdm),scq2(idm,jdm)
     .,scuxi(idm,jdm),scvyi(idm,jdm)
     .,scp2i(idm,jdm),scq2i(idm,jdm)
c
     .,pgfx(idm,jdm),pgfy(idm,jdm)
     .,gradx(idm,jdm),grady(idm,jdm)
     .,depthu(idm,jdm),depthv(idm,jdm)
     .,pvtrop(idm,jdm)
     .,drag(idm,jdm)
     .,glue(idm,jdm)
     .,zone(idm,jdm,3)
     .,dampu(idm,jdm),dampv(idm,jdm) )
c
       allocate(
     . uja(idm,jdm),ujb(idm,jdm)
     .,via(idm,jdm),vib(idm,jdm)
     .,pbot(idm,jdm)
     .,tracer(idm,jdm,kdm,ntrcr)
     .,diadff(idm,jdm,kdm)
     .,tprime(idm,jdm)
     .,sgain(idm,kdm)
     .,surflx(idm,jdm)
     .,salflx(idm,jdm)
     .,sflxcum(idm,jdm)
     .,hflxcum(idm,jdm)
c    .,thkice(idm,jdm)
c    .,covice(idm,jdm)
c    .,temice(idm,jdm)
     .,odhsi(idm,jdm)
     .,odmsi(idm,jdm)
     .,omlhc(idm,jdm)
     .,dmfz(idm,jdm) )
c
      allocate( klist(idm,jdm),ijlist(idm,jdm) )
c
      allocate(
     . taux(idm,jdm)
     .,tauy(idm,jdm)
c    .,wndspd(idm,jdm,4)
c    .,airtmp(idm,jdm,4)
c    .,vapmix(idm,jdm,4)
c    .,oprec(idm,jdm)
c    .,oevap(idm,jdm)
     .,oemnp(idm,jdm)
     .,oicemlt(idm,jdm)
     .,oflxa2o(idm,jdm),oice(idm,jdm)
     .,ustar(idm,jdm)
     .,ustarb(idm,jdm)
     .,osalt(idm,jdm)
     .,freshw(idm,jdm)
     .,diag1(idm,jdm),diag2(idm,jdm),diag3(idm,jdm),diag4(idm,jdm)
     .,diafor(idm,jdm) )
c

      !!return

      u = 0
      v = 0
      dp = 0
      dpold = 0
      dpu = 0
      dpv = 0
      p = 0
      pu = 0
      pv = 0
      latij = 0
      lonij = 0
      corio = 0
      potvor = 0
      temp = 0
      saln = 0
      th3d = 0
      thstar = 0
      wgtkap = 0
      psikk = 0
      thkk = 0
      dpmixl = 1.0     ! TNL: avoid NaN on the first step
      srfhgt = 0
      montg = 0
      defor1 = 0
      defor2 = 0
      ubavg = 0
      vbavg = 0
      pbavg = 0
      ubrhs = 0
      vbrhs = 0
      utotm = 0
      vtotm = 0
      utotn = 0
      vtotn = 0
      uflux = 0
      vflux = 0
      uflux1 = 0
      vflux1 = 0
      uflux2 = 0
      vflux2 = 0
      uflux3 = 0
      vflux3 = 0
      uflx = 0
      vflx = 0
      bolusu = 0
      bolusv = 0
      uav = 0
      vav = 0
      dpuav = 0
      dpvav = 0
      temav = 0
      salav = 0
      th3av = 0
      dpav = 0
      pbavav = 0
      sfhtav = 0
      uflxav = 0
      vflxav = 0
      ufxavp = 0
      vfxavp = 0
      diaflx = 0
      salflav = 0
      brineav = 0
      eminpav = 0
      surflav = 0
      tauxav = 0
      tauyav = 0
      ufxcum = 0
      vfxcum = 0
      dpinit = 0
      dpmxav = 0
      oiceav = 0
      util1 = 0
      util2 = 0
      util3 = 0
      util4 = 0
      scpx = 0
      scpy = 0
      scux = 0
      scuy = 0
      scvx = 0
      scvy = 0
      scqx = 0
      scqy = 0
      scu2 = 0
      scv2 = 0
      scp2 = 0
      scq2 = 0
      scuxi = 0
      scvyi = 0
      scp2i = 0
      scq2i = 0
      pgfx = 0
      pgfy = 0
      gradx = 0
      grady = 0
      depthu = 0
      depthv = 0
      pvtrop = 0
      drag = 0
      glue = 0
      zone = 0
      dampu = 0
      dampv = 0
      uja = 0
      ujb = 0
      via = 0
      vib = 0
      pbot = 0
      tracer = 0
      diadff = 0
      tprime = 0
      sgain = 0
      surflx = 0
      salflx = 0
      sflxcum = 0
      hflxcum = 0
      odhsi = 0
      odmsi = 0
      omlhc = 0
      dmfz = 0
      taux = 0
      tauy = 0
      oemnp = 0
      oicemlt = 0
      oflxa2o = 0
      oice = 0
      ustar = 0
      ustarb = 0
      osalt = 0
      freshw = 0
      diafor = 0
      klist = 0
      ijlist = 0
      diag1=0.; diag2=0.; diag3=0.; diag4=0.

      end subroutine alloc_hycom_arrays_glob


      subroutine hycom_arrays_checksum_glob

      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(u(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(v(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dp(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpold(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpu(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpv(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(p(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pu(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pv(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(latij(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(lonij(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(corio(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(potvor(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(temp(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(saln(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(th3d(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(thstar(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(wgtkap(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(psikk(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(thkk(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpmixl(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(srfhgt(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(montg(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(defor1(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(defor2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ubavg(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vbavg(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pbavg(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ubrhs(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vbrhs(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(utotm(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vtotm(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(utotn(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vtotn(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uflux(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vflux(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uflux1(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vflux1(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uflux2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vflux2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uflux3(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vflux3(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uflx(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vflx(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(bolusu(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(bolusv(:,:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpuav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpvav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(temav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(salav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(th3av(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pbavav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(sfhtav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uflxav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vflxav(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ufxavp(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vfxavp(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(diaflx(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(salflav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(brineav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(eminpav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(surflav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(tauxav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(tauyav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ufxcum(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vfxcum(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpinit(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dpmxav(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(oiceav(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(util1(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(util2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(util3(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(util4(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scpx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scpy(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scux(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scuy(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scvx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scvy(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scqx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scqy(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scu2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scv2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scp2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scq2(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scuxi(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scvyi(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scp2i(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(scq2i(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pgfx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pgfy(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(gradx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(grady(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(depthu(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(depthv(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pvtrop(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(depths(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(drag(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(glue(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(zone(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dampu(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dampv(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(uja(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ujb(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(via(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(vib(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(pbot(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(tracer(:,:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(diadff(:,:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(tprime(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(sgain(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(surflx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(salflx(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(sflxcum(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(hflxcum(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(odhsi(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(odmsi(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(omlhc(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(dmfz(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(klist(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ijlist(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(taux(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(tauy(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(oemnp(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(oicemlt(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(oflxa2o(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(oice(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ustar(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(ustarb(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(osalt(:,:))
c
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(freshw(:,:))
      write(801,*) 'hycom_arrays_glob.f ',__LINE__,sum(diafor(:,:))

      end subroutine hycom_arrays_checksum_glob

      end module hycom_arrays_glob
c
c> Revision history:
c>
c> July 1997 - eliminated 3-D arrays -uold,vold- (used in time smoothing)
c-----------------------------------------------------------------------------
