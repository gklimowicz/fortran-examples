      module kprf_arrays
      USE HYCOM_DIM, only : idm, jdm, kdm
      USE HYCOM_DIM, only : I_0H,I_1H,J_0H,J_1H

      implicit none

      private

      public alloc_kprf_arrays
      public alloc_kprf_arrays_local
      public gather_kprf_arrays
      public scatter_kprf_arrays

      public
     & klstsv,        ! copy of k-index for momtum
     & jerlov,         ! jerlov water type 1-5
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb,           ! saln in lyr. containing mlb
     & hekman,        ! ekman layer thickness
     & hmonob,        ! monin-obukhov length
     & dpbl,          ! turbulent boundary layer depth
     & dpbbl,         ! bottom turbulent boundary layer depth
     & dpmold,        ! mixed layer depth at old time step
     & tmix,          ! mixed layer temperature
     & smix,          ! mixed layer salinity
     & thmix,         ! mixed layer potential density 
     & umix,  vmix,    ! mixed layer velocity
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac,         ! fract. of penetr. red light
     & akpar,          ! photosynthetically available radiation coefficent
     & zgrid,          !  grid levels in meters
     & vcty,           !  vert. viscosity coefficient
     & difs,           !  vert. scalar diffusivity
     & dift,           !  vert. temperature diffusivity
     & ghats,          !  nonlocal transport
     & buoflx,        ! surface buoyancy flux
     & bhtflx,        ! surface buoyancy flux from heat
     & mixflx,        ! mixed layer thermal energy flux
     & sswflx,         ! short-wave rad'n flux at surface
     & lookup,
     & irimax,
     & nb,
     & ribtbl,
     & ridb,
     & slq2b
     &,dri
     &,smb
     &,shb
     &,ssb
     &,back_ra_r
     &,sisamax
     &,ra_rmax
     &,c_y_r0
     &,sm_r1
     &,sh_r1
     &,ss_r1
     &,slq2_r1
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn

     &,klstsv_loc, jerlov_loc
     &,t1sav_loc,s1sav_loc,tmlb_loc,smlb_loc
     &,hekman_loc,hmonob_loc,dpbl_loc,dpbbl_loc
     &,dpmold_loc,tmix_loc,smix_loc,thmix_loc,umix_loc
     &,vmix_loc
     &,akpar_loc
     &,zgrid_loc,vcty_loc,difs_loc,dift_loc,ghats_loc
     &,buoflx_loc,bhtflx_loc,mixflx_loc,sswflx_loc 

!! end public

      integer, allocatable, dimension (:,:) ::
     & klstsv,        ! copy of k-index for momtum
     & jerlov         ! jerlov water type 1-5

!!      common /kpp_integ_arrays/
!!     .  klstsv,jerlov

      real, allocatable, dimension (:,:,:) :: 
     & t1sav,         ! upper sublayer temperature
     & s1sav,         ! upper sublayer salinity
     & tmlb,          ! temp in lyr. containing mlb.
     & smlb           ! saln in lyr. containing mlb

      real, allocatable, dimension (:,:) ::
     & hekman,        ! ekman layer thickness
     & hmonob,        ! monin-obukhov length
     & dpbl,          ! turbulent boundary layer depth
     & dpbbl,         ! bottom turbulent boundary layer depth
     & dpmold,        ! mixed layer depth at old time step
     & tmix,          ! mixed layer temperature
     & smix,          ! mixed layer salinity
     & thmix,         ! mixed layer potential density 
     & umix,  vmix    ! mixed layer velocity

      real, dimension (5) ::
     & betard,        ! red  extinction coefficient
     & betabl,        ! blue extinction coefficient
     & redfac         ! fract. of penetr. red light

      real, allocatable, dimension (:,:,:) ::
     & akpar          ! photosynthetically available radiation coefficent

c --- kpp variables
      real, allocatable, dimension (:,:,:) ::
     & zgrid          !  grid levels in meters
     &,vcty           !  vert. viscosity coefficient
     &,difs           !  vert. scalar diffusivity
     &,dift           !  vert. temperature diffusivity
     &,ghats          !  nonlocal transport

      real, allocatable, dimension (:,:) ::
     & buoflx,        ! surface buoyancy flux
     & bhtflx,        ! surface buoyancy flux from heat
     & mixflx,        ! mixed layer thermal energy flux
     & sswflx         ! short-wave rad'n flux at surface

!local to processors
      integer, allocatable, dimension (:,:) ::
     & klstsv_loc,        ! copy of k-index for momtum
     & jerlov_loc         ! jerlov water type 1-5

      real, allocatable, dimension (:,:,:) ::
     & t1sav_loc,         ! upper sublayer temperature
     & s1sav_loc,         ! upper sublayer salinity
     & tmlb_loc,          ! temp in lyr. containing mlb.
     & smlb_loc           ! saln in lyr. containing mlb

      real, allocatable, dimension (:,:) ::
     & hekman_loc,        ! ekman layer thickness
     & hmonob_loc,        ! monin-obukhov length
     & dpbl_loc,          ! turbulent boundary layer depth
     & dpbbl_loc,         ! bottom turbulent boundary layer depth
     & dpmold_loc,        ! mixed layer depth at old time step
     & tmix_loc,          ! mixed layer temperature
     & smix_loc,          ! mixed layer salinity
     & thmix_loc,         ! mixed layer potential density
     & umix_loc,  vmix_loc    ! mixed layer velocity

      real, allocatable, dimension (:,:,:) ::
     & akpar_loc          ! photosynthetically available radiation coefficent

c --- kpp variables
      real, allocatable, dimension (:,:,:) ::
     & zgrid_loc          !  grid levels in meters
     &,vcty_loc           !  vert. viscosity coefficient
     &,difs_loc           !  vert. scalar diffusivity
     &,dift_loc           !  vert. temperature diffusivity
     &,ghats_loc          !  nonlocal transport

      real, allocatable, dimension (:,:) ::
     & buoflx_loc,        ! surface buoyancy flux
     & bhtflx_loc,        ! surface buoyancy flux from heat
     & mixflx_loc,        ! mixed layer thermal energy flux
     & sswflx_loc         ! short-wave rad'n flux at surface

!!      common /kpp_real_arrays/
!!     .  t1sav,s1sav,tmlb,smlb,hekman,hmonob,dpbl,dpbbl,
!!     .  dpmold,tmix,smix,thmix,umix,vmix,betard,betabl,redfac,akpar,
!!     .  zgrid,vcty,difs,dift,ghats,buoflx,bhtflx,mixflx,sswflx
c
      integer,parameter :: lookup=762
      integer, save ::
     & irimax(-lookup:lookup)
     &,nb

!!      common/giss_int2/
!!     &        irimax,nb
!!      save  /giss_int2/
c
      real, save ::
     & ribtbl(-lookup:lookup)
     &,ridb(  -lookup:lookup)
     &,slq2b( -lookup:lookup,-lookup:lookup)
     &,dri
     &,smb(   -lookup:lookup,-lookup:lookup)
     &,shb(   -lookup:lookup,-lookup:lookup)
     &,ssb(   -lookup:lookup,-lookup:lookup)
     &,back_ra_r(-39:117)
     &,sisamax(  -39:117)
     &,ra_rmax(  -39:117)
     &,c_y_r0(   -39:117)
     &,sm_r1(    -39:117)
     &,sh_r1(    -39:117)
     &,ss_r1(    -39:117)
     &,slq2_r1(  -39:117)
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn

!!      common/giss_real2/
!!     &        ribtbl,ridb,slq2b,dri,smb,shb,ssb,back_ra_r,
!!     &        sisamax,ra_rmax,c_y_r0,sm_r1,sh_r1,ss_r1,slq2_r1,
!!     &        b1,visc_cbu_limit,diff_cbt_limit,
!!     &        theta_rcrp,theta_rcrn
!!      save  /giss_real2/

      contains

!!! switch to this when doing parallelization
      subroutine alloc_kprf_arrays_local

      allocate(
     & klstsv_loc(I_0H:I_1H,J_0H:J_1H),
     & jerlov_loc(I_0H:I_1H,J_0H:J_1H) )

      allocate(
     & t1sav_loc(I_0H:I_1H,J_0H:J_1H,2),
     & s1sav_loc(I_0H:I_1H,J_0H:J_1H,2),
     & tmlb_loc(I_0H:I_1H,J_0H:J_1H,2),
     & smlb_loc(I_0H:I_1H,J_0H:J_1H,2) )

      allocate(
     & hekman_loc(I_0H:I_1H,J_0H:J_1H),
     & hmonob_loc(I_0H:I_1H,J_0H:J_1H),
     & dpbl_loc(I_0H:I_1H,J_0H:J_1H),
     & dpbbl_loc(I_0H:I_1H,J_0H:J_1H),
     & dpmold_loc(I_0H:I_1H,J_0H:J_1H),
     & tmix_loc(I_0H:I_1H,J_0H:J_1H),
     & smix_loc(I_0H:I_1H,J_0H:J_1H),
     & thmix_loc(I_0H:I_1H,J_0H:J_1H),
     & umix_loc(I_0H:I_1H,J_0H:J_1H),  vmix_loc(I_0H:I_1H,J_0H:J_1H) )

      allocate(
     & akpar_loc(I_0H:I_1H,J_0H:J_1H,12) )

      allocate(
     & zgrid_loc(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,vcty_loc(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,difs_loc(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,dift_loc(I_0H:I_1H,J_0H:J_1H,kdm+1)
     &,ghats_loc(I_0H:I_1H,J_0H:J_1H,kdm+1) )

      allocate(
     & buoflx_loc(I_0H:I_1H,J_0H:J_1H),
     & bhtflx_loc(I_0H:I_1H,J_0H:J_1H),
     & mixflx_loc(I_0H:I_1H,J_0H:J_1H),
     & sswflx_loc(I_0H:I_1H,J_0H:J_1H) )

      !initialize to 0
      klstsv_loc=0; jerlov_loc=0;
      t1sav_loc=0; s1sav_loc=0; tmlb_loc=0; smlb_loc=0;
      hekman_loc=0; hmonob_loc=0; dpbl_loc=0; dpbbl_loc=0;
      dpmold_loc=0; tmix_loc=0; 
      smix_loc=0; thmix_loc=0; umix_loc=0; vmix_loc=0;
      akpar_loc=0;
      zgrid_loc=0; vcty_loc=0; difs_loc=0; dift_loc=0; ghats_loc=0;
      buoflx_loc=0; bhtflx_loc=0; mixflx_loc=0; sswflx_loc=0;

      end subroutine alloc_kprf_arrays_local

      subroutine alloc_kprf_arrays

      allocate(
     & klstsv(IDM,JDM),
     & jerlov(IDM,JDM) )

      allocate(
     & t1sav(IDM,JDM,2),
     & s1sav(IDM,JDM,2),
     & tmlb(IDM,JDM,2),
     & smlb(IDM,JDM,2) )

      allocate(
     & hekman(IDM,JDM),
     & hmonob(IDM,JDM),
     & dpbl(IDM,JDM),
     & dpbbl(IDM,JDM),
     & dpmold(IDM,JDM),
     & tmix(IDM,JDM),
     & smix(IDM,JDM),
     & thmix(IDM,JDM),
     & umix(IDM,JDM),  vmix(IDM,JDM) )

      allocate(
     & akpar(IDM,JDM,12) )

      allocate(
     & zgrid(IDM,JDM,kdm+1)
     &,vcty(IDM,JDM,kdm+1)
     &,difs(IDM,JDM,kdm+1)
     &,dift(IDM,JDM,kdm+1)
     &,ghats(IDM,JDM,kdm+1) )

      allocate(
     & buoflx(IDM,JDM),
     & bhtflx(IDM,JDM),
     & mixflx(IDM,JDM),
     & sswflx(IDM,JDM) )

      !initialize to 0
      klstsv=0; jerlov=0;
      t1sav=0; s1sav=0; tmlb=0; smlb=0;
      hekman=0; hmonob=0; dpbl=0; dpbbl=0; dpmold=0; tmix=0; 
      smix=0; thmix=0; umix=0; vmix=0;
      akpar=0;
      zgrid=0; vcty=0; difs=0; dift=0; ghats=0;
      buoflx=0; bhtflx=0; mixflx=0; sswflx=0;

      end subroutine alloc_kprf_arrays

      subroutine gather_kprf_arrays
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA

      call pack_data( ogrid,  klstsv_loc, klstsv )
      call pack_data( ogrid,  jerlov_loc, jerlov )
      call pack_data( ogrid,  t1sav_loc,  t1sav )
      call pack_data( ogrid,  s1sav_loc,  s1sav )
      call pack_data( ogrid,  tmlb_loc,   tmlb )
      call pack_data( ogrid,  smlb_loc,   smlb )
      call pack_data( ogrid,  hekman_loc, hekman )
      call pack_data( ogrid,  hmonob_loc, hmonob )
      call pack_data( ogrid,  dpbl_loc,   dpbl )
      call pack_data( ogrid,  dpbbl_loc,  dpbbl )
      call pack_data( ogrid,  dpmold_loc, dpmold )
      call pack_data( ogrid,  tmix_loc,   tmix )
      call pack_data( ogrid,  smix_loc,   smix )
      call pack_data( ogrid,  thmix_loc,  thmix )
      call pack_data( ogrid,  umix_loc,   umix )
      call pack_data( ogrid,  vmix_loc,   vmix )
      call pack_data( ogrid,  akpar_loc,  akpar )
      call pack_data( ogrid,  zgrid_loc,  zgrid )
      call pack_data( ogrid,  vcty_loc,   vcty )
      call pack_data( ogrid,  difs_loc,   difs )
      call pack_data( ogrid,  dift_loc,   dift )
      call pack_data( ogrid,  ghats_loc,  ghats )
      call pack_data( ogrid,  buoflx_loc, buoflx )
      call pack_data( ogrid,  bhtflx_loc, bhtflx )

      call pack_data( ogrid,  mixflx_loc, mixflx )
      call pack_data( ogrid,  sswflx_loc, sswflx )

      end subroutine gather_kprf_arrays

      subroutine scatter_kprf_arrays
      USE HYCOM_DIM, only : ogrid
      USE DOMAIN_DECOMP_1D, ONLY: UNPACK_DATA

      call unpack_data( ogrid,  klstsv, klstsv_loc )
      call unpack_data( ogrid,  jerlov, jerlov_loc )
      call unpack_data( ogrid,  t1sav,  t1sav_loc )
      call unpack_data( ogrid,  s1sav,  s1sav_loc )
      call unpack_data( ogrid,  tmlb,   tmlb_loc )
      call unpack_data( ogrid,  smlb,   smlb_loc )
      call unpack_data( ogrid,  hekman, hekman_loc )
      call unpack_data( ogrid,  hmonob, hmonob_loc )
      call unpack_data( ogrid,  dpbl,   dpbl_loc )
      call unpack_data( ogrid,  dpbbl,  dpbbl_loc )
      call unpack_data( ogrid,  dpmold, dpmold_loc )
      call unpack_data( ogrid,  tmix,   tmix_loc )
      call unpack_data( ogrid,  smix,   smix_loc )
      call unpack_data( ogrid,  thmix,  thmix_loc )
      call unpack_data( ogrid,  umix,   umix_loc )
      call unpack_data( ogrid,  vmix,   vmix_loc )
      call unpack_data( ogrid,  akpar,  akpar_loc )
      call unpack_data( ogrid,  zgrid,  zgrid_loc )
      call unpack_data( ogrid,  vcty,   vcty_loc )
      call unpack_data( ogrid,  difs,   difs_loc )
      call unpack_data( ogrid,  dift,   dift_loc )
      call unpack_data( ogrid,  ghats,  ghats_loc )
      call unpack_data( ogrid,  buoflx, buoflx_loc )
      call unpack_data( ogrid,  bhtflx, bhtflx_loc )

      call unpack_data( ogrid,  mixflx, mixflx_loc )
      call unpack_data( ogrid,  sswflx, sswflx_loc )

      end subroutine scatter_kprf_arrays

      end module kprf_arrays
