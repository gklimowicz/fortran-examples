
c**** COMMENT program when linked with MODELE
      program testregrid
cc      call test_zonal_loop()    ! tests parallel zonal mean code, using rolled loops
cc      call test_zonal_unrolled() ! tests parallel zonal mean code, using unrolled loops
      call test_regrid()        ! tests parallel regridding code
cc      call offregrid()    ! tests offline regridding for input file
cc      call test_or()
cc      call test_onlineread()
cc      call test_rr2()
cc        call exact_regrid()
cc      call add_dd()
      end program testregrid
c****


      subroutine test_or()
      USE regrid_com, only:is,ie,js,je
      use gatscat_mod

      implicit none
      integer :: npes,ndims
      real*8, dimension(:,:), allocatable :: P

c***
      call init_grid_temp()

c***  Initialize gather&scatter
      call domain_decomp_init
      call gatscat_init
c***  Initialize exchange grid for zonal means
      call init_xgrid_loop()

      open(UNIT=10, FILE="AIC.RES_F20.D771201", 
     *     FORM="UNFORMATTED", 
     *     STATUS="OLD",CONVERT='BIG_ENDIAN')

      allocate(P(is:ie,js:je))
      
      call readt_regrid_parallel(10,"AIC",0,P,1)
      write(6,*) P(:,:)

      deallocate(P)

      end subroutine test_or
c****


      subroutine test_onlineread()
      USE regrid_com, only:ic,jc
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop
     &     ,fv_domain_decomp=>domain_decomp
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain, tile, npes_x, npes_y
      use fv_grid_utils_mod, only : cosa_s,sina_s
     &     ,grid_utils_init
c     &     ,sw_corner,se_corner,ne_corner,nw_corner
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type, dx_const, dy_const
      use fv_control_mod, only : uniform_ppm,c2l_ord
      use gs_domain_decomp, ng=>halo_width
      use gatscat_mod

      implicit none
      integer :: i,j,npz
      integer :: npes,ndims
      integer :: commID
      integer, dimension(:), allocatable :: pelist

      real*8, dimension(:,:), allocatable :: P

      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline'
      logical :: non_ortho
      
      write(6,*) "UNIT TEST : readt_regrid_parallel"

c***  Temporarily use instanciation of grid through fv_grid_tools_mod's init_grid()

      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
c      write(6,*) 'commID ',commID
      npx = ic; npy = npx ! 1x1 resolution
      ng = 3 ! number of ghost zones required

      call fv_domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)
c      call mp_stop()

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

c***  Initialize gather&scatter
      call domain_decomp_init
      write(6,*) "AFT. dpmain_decomp_init"
      call gatscat_init
      write(6,*) "AFT. gatscat_init"


c***  Initialize exchange grid for zonal means
      call init_xgrid_loop()

      open(UNIT=10, FILE="AIC.RES_F20.D771201", 
     *     FORM="UNFORMATTED", 
     *     STATUS="OLD",CONVERT='BIG_ENDIAN')

      allocate(P(is:ie,js:je))
      
      call readt_regrid_parallel(10,"AIC",0,P,1)
      write(6,*) ">>>>>>>>>ARRIJ GID>>>>>>>>>>",gid
      write(6,*) P(:,:)

      call mpp_exit()
      deallocate(pelist)
      deallocate(P)

      end subroutine test_onlineread
c****

c****
      subroutine test_zonal_unrolled()
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop,domain_decomp
      use fv_grid_utils_mod, only : grid_utils_init
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type
      use fv_control_mod, only : uniform_ppm, c2l_ord

      implicit none
      include 'netcdf.inc'
      integer :: npx,npy,npes,ng,ndims,rootpe
      integer :: commID
      integer, dimension(:), allocatable :: pelist

      integer :: l,npz
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline', ofi
      logical :: non_ortho
      integer, parameter :: nedge=4 ! number of edges on each face
      integer :: nij_latlon,nc2,loctile
      real*8 :: gsum

      real*8, allocatable :: tcub_loc(:,:)
      real*8 :: tcub(ic,jc)
      real*8 :: area_latband(jlatm),zonal_mean(jlatm)
      integer :: itile,fid,vid,srt(3),cnt(3),status,ikey,jlat,i,j
      character*200 :: infi,infile
      character*1 :: istr 

c
c     Initialize grid and domain decomposition
c
      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
      write(6,*) 'rootpe= ',mpp_root_pe()
c      npx = 90; npy = npx ! 1x1 resolution
      npx = ic; npy = npx 
      
      ng = 3 ! number of ghost zones required
      call domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

c      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

      rootpe=mpp_root_pe()

c
c     Initialize exchange grid and extract information
c
      call init_xgrid_unrolled()

      write(6,*) "HERE"

ccc    remove these lines - this will be initialize somewhere else
      allocate(tcub_loc(isd:ied,jsd:jed))

      infi='atmos_daily.tile'
      loctile=gid/dom_per_tile
      loctile=loctile+1
      write(istr,'(i1)') loctile
      infile=trim(infi)//istr//'.nc'
      write(*,*) infile
      status = nf_open(trim(infile),nf_nowrite,fid)
      write(*,*) "status=",status
      status = nf_inq_varid(fid,'t_surf',vid)
      if (status .ne. NF_NOERR) write(*,*) "ERROR"
      srt = (/ 1, 1, 1 /)
      cnt = (/ ic, jc, 1 /)
      status = nf_get_vara_double(fid,vid,srt,cnt,tcub(1,1))
      status = nf_close(fid)
      tcub_loc(is:ie,js:je)=tcub(is:ie,js:je)

c      tcub_loc(:,:)=1  !gid
ccc    end remove

      call zonalmean_cs_unrolled(tcub_loc(:,:),
     *       zonal_mean(:),
     *       area_latband(:))

      write(6,*) "THERE"

c     should be pre-computed
      call mpp_sum(area_latband,jlatm)

c     line below should be moved inside zaonal mean routine
      call mpp_sum(zonal_mean,jlatm)
         

      if (gid .eq. rootpe) then
      do jlat=1,jlatm
         if (area_latband(jlat) .gt. 1.d-14) then
            zonal_mean(jlat)=zonal_mean(jlat)
     *           /area_latband(jlat)
         endif
         write(*,*) "zonal_mean unrolled=",zonal_mean(jlat)
      enddo
      endif

      call mpp_exit()
      deallocate(pelist)
      deallocate(tcub_loc)
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine test_zonal_unrolled
c****


      subroutine test_zonal_loop()
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop,domain_decomp
      use fv_mp_mod, only : gid, domain
      use fv_grid_utils_mod, only : grid_utils_init
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type
      use fv_control_mod, only : uniform_ppm, c2l_ord

      implicit none
      include 'netcdf.inc'
      integer :: npx,npy,npes,ng,ndims,rootpe
      integer :: commID
      integer, dimension(:), allocatable :: pelist

      integer :: l,npz
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline', ofi
      logical :: non_ortho
      integer, parameter :: nedge=4 ! number of edges on each face
      integer :: nij_latlon,nc2,loctile
      real*8 :: gsum

      real*8, allocatable :: tcub_loc(:,:)
      real*8 :: tcub(ic,jc)
      real*8 :: area_latband(jlatm),zonal_mean(jlatm)
      integer :: itile,fid,vid,srt(3),cnt(3),status,ikey,jlat,i,j
      character*200 :: infi,infile
      character*1 :: istr 

c
c     Initialize grid and domain decomposition
c
      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
      write(6,*) 'rootpe= ',mpp_root_pe()
c      npx = 90; npy = npx ! 1x1 resolution
      npx = ic; npy = npx 
      
      ng = 3 ! number of ghost zones required
      call domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

c      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

      rootpe=mpp_root_pe()

c
c     Initialize exchange grid and extract information
c
      call init_xgrid_loop()

      write(6,*) "HERE"

ccc    remove these lines - this will be initialize somewhere else
      allocate(tcub_loc(isd:ied,jsd:jed))

      infi='atmos_daily.tile'
      loctile=gid/dom_per_tile
      loctile=loctile+1
      write(istr,'(i1)') loctile
      infile=trim(infi)//istr//'.nc'
      write(*,*) infile
      status = nf_open(trim(infile),nf_nowrite,fid)
      write(*,*) "status=",status
      status = nf_inq_varid(fid,'t_surf',vid)
      if (status .ne. NF_NOERR) write(*,*) "ERROR"
      srt = (/ 1, 1, 1 /)
      cnt = (/ ic, jc, 1 /)
      status = nf_get_vara_double(fid,vid,srt,cnt,tcub(1,1))
      status = nf_close(fid)
      tcub_loc(is:ie,js:je)=tcub(is:ie,js:je)

c      tcub_loc(:,:)=1  !gid
ccc    end remove


c     zonal mean inside rolled loop

      area_latband(:)=0.0
      zonal_mean(:)=0.0

      do i=is,ie
         do j=js,je
            call zonalmean_cs_loop(tcub_loc(i,j),i,j,
     *           zonal_mean(:),area_latband(:))
            enddo
         enddo

c     sum over all domains (not reducing rank)

      call mpp_sum(area_latband,jlatm)
      call mpp_sum(zonal_mean,jlatm)


      if (gid .eq. rootpe) then
      do jlat=1,jlatm
         if (area_latband(jlat) .gt. 1.d-14) then
c        line beloz should be moved in the weights inside zonal mean routine
            zonal_mean(jlat)=zonal_mean(jlat)
     *           /area_latband(jlat)
         endif
         write(*,*) "zonal_mean rolled=",zonal_mean(jlat)
      enddo
      endif

      call mpp_exit()
      deallocate(pelist)
      deallocate(tcub_loc)
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine test_zonal_loop
c****


      subroutine test_regrid()
      use regrid_com
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop,domain_decomp
      use fv_mp_mod, only : domain
      use fv_grid_utils_mod, only : grid_utils_init
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type
      use fv_control_mod, only : uniform_ppm, c2l_ord

      implicit none
      include 'netcdf.inc'

      real*8 :: tlatlon(ilonm,jlatm),alatlon(ilonm,jlatm)
      integer :: i,j,n,itile,rootpe,ierr
      integer :: nfaces,npx,npy,npes,ng,ndims
      integer :: commID, status, fid, vid
      integer, dimension(:), allocatable :: pelist
      real*8, dimension(:,:), allocatable :: tcub_loc

      integer :: l,npz
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline', ofi
      logical :: non_ortho
      integer, parameter :: nedge=4 ! number of edges on each face
      integer :: nij_latlon,nc2
      real*8 :: gsum


      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)

      npx = ic; npy = npx 
      
      ng = 3 ! number of ghost zones required
      call domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

c read weights and initialize target array to zero
      call init_regrid()

      allocate(tcub_loc(isd:ied,jsd:jed))
      tcub_loc(:,:)=gid

      call parallel_regrid_cs2ll(tcub_loc,tlatlon,alatlon)   

      if (gid .eq. 0) then
      ofi='tstout.nc'
      status = nf_open(trim(ofi),nf_write,fid)
      if (status .ne. NF_NOERR) write(*,*) NF_STRERROR(status)
      status = nf_inq_varid(fid,'lwup_sfc',vid)
      write(*,*) NF_STRERROR(status)
      status = nf_put_var_double(fid,vid,tlatlon)
      write(*,*) NF_STRERROR(status)
     
      status = nf_close(fid)

      call mpp_exit()

      endif

      deallocate(pelist)
      deallocate(tcub_loc)
      deallocate(xgrid_area)
      deallocate(ijcub)
      deallocate(ijlatlon)
      deallocate(tile)

      end subroutine test_regrid
c*


      subroutine test_rr2()
      use regrid_com, only : ilonm,jlatm,ic,jc
      use mpp_mod
      use mpp_domains_mod
      use fv_mp_mod, only : mp_start,mp_stop
     &     ,fv_domain_decomp=>domain_decomp
      use fv_mp_mod, only : is, ie, js, je, isd, ied, jsd, jed
      use fv_mp_mod, only : gid, domain, tile, npes_x, npes_y
      use fv_grid_utils_mod, only : cosa_s,sina_s
     &     ,grid_utils_init
c     &     ,sw_corner,se_corner,ne_corner,nw_corner
      use fv_arrays_mod, only: fv_atmos_type
      use fv_grid_tools_mod,  only: init_grid, cosa, sina, area, area_c,
     &     dx, dy, dxa, dya, dxc, dyc, grid_type, dx_const, dy_const
      use fv_control_mod, only : uniform_ppm,c2l_ord
      use gs_domain_decomp, ng=>halo_width
      use gatscat_mod

      implicit none
      integer :: i,j,npz
      integer :: npes,ndims
      integer :: commID
      integer, dimension(:), allocatable :: pelist

c pack
      real*8, dimension(:,:), allocatable :: arrij
      character*80 :: TITLE

      real*4 :: tllr4(ilonm,jlatm)  ! latlon real*4 data read from input file
      real*4 :: X(ilonm,jlatm)
      real*8 :: tdatall(ilonm,jlatm)  ! latlon real*8 data 
      real*8 :: tcubglob(ic,jc,6)  ! global array 
      integer :: n,ipos,ierr,nskip
      type(fv_atmos_type) :: atm
      character*80 :: grid_name = 'Gnomonic'
      character*120:: grid_file = 'Inline'
      logical :: non_ortho

      call mpp_init(MPP_VERBOSE)
c code copied from fv_init:
      npes = mpp_npes()
      allocate(pelist(npes))
      call mpp_get_current_pelist( pelist, commID=commID )
      call mp_start(commID)
c      write(6,*) 'commID ',commID
      npx = ic; npy = npx ! 1x1 resolution
      ng = 3 ! number of ghost zones required

      call fv_domain_decomp(npx+1,npy+1,ntiles,ng,grid_type)
c      call mp_stop()

      ndims = 2
      npz = 5
      call init_grid(atm,grid_name,grid_file,
     &     npx+1, npy+1, npz, ndims, ntiles, ng)

      non_ortho=.true.
      call grid_utils_init(Atm, npx+1, npy+1, npz, Atm%grid, Atm%agrid,
     &     area, area_c, cosa, sina, dx, dy, dxa, dya, non_ortho,
     &     uniform_ppm, grid_type, c2l_ord)

      call domain_decomp_init
      call gatscat_init

c      write(6,*) 'indices', is, ie, js, je, isd, ied, jsd, jed

c***  Initialize exchange grid for zonal means
      call init_xgrid_loop()

      write(6,*) "AFT. init xgrid"

      open(UNIT=10, FILE="AIC.RES_F20.D771201", 
     *     FORM="UNFORMATTED", 
     *     STATUS="OLD",CONVERT='BIG_ENDIAN')

      ipos=1
      nskip=0

      if (gid .eq. 0) then
         do n=1,ipos-1
            read(UNIT=10,IOSTAT=ierr)
         enddo
         read(UNIT=10,IOSTAT=ierr) TITLE, (X,n=1,nskip), tllr4
c     convert from real*4 to real*8
         tdatall=tllr4
     
c
c     Regrid from lat-lon to cubbed sphere, form global array
c

         call root_regrid_ll2cs(tdatall,tcubglob)

      endif

      allocate(arrij   (is:ie,js:je))
      
      call unpack_data(tcubglob,arrij)
      write(6,*) ">>>>>>>>>ARRIJ GID>>>>>>>>>>",gid
      write(6,*) arrij

      call mpp_exit()

      end subroutine test_rr2
c****


