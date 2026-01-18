#include "rundeck_opts.h"
!@auth M. Kelley
!@sum  cubed-atmosphere versions of atm-ocean regridding routines

      module hycom_cpler
      use cs2ll_utils, only : xgridremap_type,cs2llint_type,
     &     init_xgridremap_type,xgridremap_ij,
     &     init_cs2llint_type2,cs2llint_lluv
      USE HYCOM_DIM_GLOB, only : iio,jjo
      USE HYCOM_DIM, only : agrid,ogrid
     &    ,aI_0,aI_1, aJ_0,aJ_1, aI_0H,aI_1H, aJ_0H,aJ_1H
     &    , J_0,  J_1,  J_0H,  J_1H
     &    ,isp,ifp,ilp,ip
      USE CONSTANT, only : tf
      implicit none
      save
      private

      public ssto2a,veca2o,flxa2o,veco2a,tempro2a,cpl_wgt
      public coso_glob,sino_glob

!@var remap_flxa2o,remap_ssto2a atm->ocn,ocn->atm regrid info
      type(xgridremap_type) ::     ! uses area weighting
     &     remap_flxa2o            ! atm A -> ocn A
     &    ,remap_ssto2a            ! ocn A -> atm A
!@var remap_veca2o atm->ocn regrid info
      type(cs2llint_type) ::       ! bilinear interpolation
     &     remap_veca2o            ! atm A -> ocn A

      real*8, dimension(:,:), allocatable :: coso,sino
      real*8, dimension(:,:), allocatable :: coso_glob,sino_glob

      contains

      subroutine cpl_wgt
      USE CONSTANT, ONLY : RADIAN
      USE DOMAIN_DECOMP_1D, only : am_i_root,unpack_data,broadcast
      USE HYCOM_SCALARS, only : flnmcoso
      USE HYCOM_DIM_GLOB, only :
     &     isp_glob=>isp,ifp_glob=>ifp,ilp_glob=>ilp
      USE HYCOM_ARRAYS_GLOB, only : lonij,latij,depths
      use filemanager, only : findunit
      implicit none
      real*8, dimension(:,:), allocatable :: lon_hycom,lat_hycom
      integer :: i,j,l,iz,jz,iu1,m,n
      integer :: ncells
      integer, dimension(:), allocatable :: ones,tile,tile2
      integer, dimension(:,:), allocatable :: ij_cubed,ij_hycom
      integer, dimension(:), allocatable :: ic2,jc2,io2,jo2
      real*8, dimension(:), allocatable :: xgrid_area,xgrid_area2
c netcdf stuff
      include 'netcdf.inc'
      integer :: status,vid,fid,did

c read rotation coeffs between hycom gridlines and geographic north
      allocate(coso(iio,j_0h:j_1h),sino(iio,j_0h:j_1h))
      if(am_i_root()) then
        allocate(coso_glob(iio,jjo),sino_glob(iio,jjo))
        call findunit(iu1)
        open(iu1,file=flnmcoso,form='unformatted',status='old')
        read(iu1) iz,jz,coso_glob,sino_glob
        close(iu1)
        if (iz.ne.iio .or. jz.ne.jjo) then
          write(6,*) ' iz,jz=',iz,jz
          stop '(wrong iz/jz in cososino.8bin)'
        endif
      endif
      call unpack_data(ogrid,coso_glob,coso)
      call unpack_data(ogrid,sino_glob,sino)

c read areas/indices for area-weighted regridding and initialize
c the corresponding data structures
      status = nf_open('aoremap',nf_nowrite,fid)
      status = nf_inq_dimid(fid,'ncells',did)
      status = nf_inq_dimlen(fid,did,ncells)
      allocate(xgrid_area(ncells))
      allocate(ij_cubed(2,ncells))
      allocate(ij_hycom(2,ncells))
      allocate(tile(ncells),ones(ncells))
      allocate(tile2(ncells),ic2(ncells),jc2(ncells))
      allocate(io2(ncells),jo2(ncells),xgrid_area2(ncells))
      status = nf_inq_varid(fid,'xgrid_area',vid)
      status = nf_get_var_double(fid,vid,xgrid_area)
      status = nf_inq_varid(fid,'tile1',vid)
      status = nf_get_var_int(fid,vid,tile)
      status = nf_inq_varid(fid,'tile1_cell',vid)
      status = nf_get_var_int(fid,vid,ij_cubed)
      status = nf_inq_varid(fid,'tile2_cell',vid)
      status = nf_get_var_int(fid,vid,ij_hycom)
      status = nf_close(fid)

      m = 0
      do n=1,ncells
        i = ij_hycom(1,n); j = ij_hycom(2,n)
        if(depths(i,j).eq.0.) cycle
        m = m + 1
        ic2(m) = ij_cubed(1,n)
        jc2(m) = ij_cubed(2,n)
        tile2(m) = tile(n)
        xgrid_area2(m) = xgrid_area(n)
        io2(m) = i
        jo2(m) = j
      enddo

      ones=1
      call init_xgridremap_type(agrid,ogrid,
c     &     ncells,
c     &     ij_cubed(1,:),ij_cubed(2,:),tile,
c     &     ij_hycom(1,:),ij_hycom(2,:),ones,
c     &     xgrid_area,
     &     m,
     &     ic2(1:m),jc2(1:m),tile2(1:m),
     &     io2(1:m),jo2(1:m),ones(1:m),
     &     xgrid_area2(1:m),
     &     remap_flxa2o)

      call init_xgridremap_type(ogrid,agrid,
     &     m,
     &     io2(1:m),jo2(1:m),ones(1:m),
     &     ic2(1:m),jc2(1:m),tile2(1:m),
     &     xgrid_area2(1:m),
     &     remap_ssto2a)
      deallocate(xgrid_area,ij_cubed,ij_hycom,tile,ones)
      deallocate(xgrid_area2,ic2,jc2,io2,jo2,tile2)

c Initialize data structure used for agcm->ogcm bilinear interpolation.
c Convert hycom lons/lats to radians, change lon domain from -pi:+pi,
c and set lat at land gridpoints to a missing value
      allocate(lon_hycom(iio,jjo),lat_hycom(iio,jjo))
      if(am_i_root()) then
        lat_hycom = -1d30
        lon_hycom = -1d30
        do j=1,jjo
          do l=1,isp_glob(j)
            do i=ifp_glob(j,l),ilp_glob(j,l)
              lon_hycom(i,j) = lonij(i,j,3)
              lat_hycom(i,j) = latij(i,j,3)
              if(lon_hycom(i,j) .gt. 180.)
     &             lon_hycom(i,j) = lon_hycom(i,j)-360.
              lon_hycom(i,j) = lon_hycom(i,j)*radian
              lat_hycom(i,j) = lat_hycom(i,j)*radian
            enddo
          enddo
        enddo
      endif
      call broadcast(ogrid,lon_hycom)
      call broadcast(ogrid,lat_hycom)
      call init_cs2llint_type2(agrid,ogrid,lon_hycom,lat_hycom,
     &     remap_veca2o,setup_rot_pol=.true.) ! setup_rot_pol for vectors
      deallocate(lon_hycom,lat_hycom)

      return
      end subroutine cpl_wgt

      subroutine ssto2a(fldo,flda)
c --- mapping sst from ogcm A grid to agcm A grid
c     input: fldo, output: flda
c
      real*8 fldo(iio,J_0H:J_1H),flda(aI_0H:aI_1H,aJ_0H:aJ_1H)
      call xgridremap_ij(remap_ssto2a, fldo, flda)
      return
      end subroutine ssto2a

      subroutine tempro2a(tco,tka)
c --- mapping sqrt(sqrt(temp**4)) from ogcm A grid to agcm A grid
c --- input: tco in deg C; outout: tka in deg K
c
      real*8 tco(iio,J_0H:J_1H),tka(aI_0H:aI_1H,aJ_0H:aJ_1H)
      real*8 tko4(iio,J_0H:J_1H)
      tko4(:,J_0:J_1) = (tco(:,J_0:J_1)+tf)**4
      call xgridremap_ij(remap_ssto2a, tko4, tka)
      tka(aI_0:aI_1,aJ_0:aJ_1) = tka(aI_0:aI_1,aJ_0:aJ_1)**.25d0
      return
      end subroutine tempro2a

      subroutine flxa2o(flda,fldo)
c --- mapping flx from agcm A grid to ogcm A grid
c     input: flda, output: fldo
c
      real*8 flda(aI_0H:aI_1H,aJ_0H:aJ_1H),fldo(iio,J_0H:J_1H)
      call xgridremap_ij(remap_flxa2o, flda, fldo)
      return
      end subroutine flxa2o
      
      subroutine veca2o(tauxa,tauya,tauxo,tauyo)
c --- mapping vector from agcm A grid to ogcm A grid
c --- input  tauxa/tauya (N/m2): E_/N_ward on agcm A grid
c --- output tauxo/tauyo (N/m2): +i_/+j_ward on ogcm A grid
c ---                            (S_/E_ward in Mercator domain)
c
      implicit none
      real*8, dimension(aI_0H:aI_1H,aJ_0H:aJ_1H) :: tauxa,tauya
      real*8, dimension(iio        , J_0H: J_1H) :: tauxo,tauyo
      integer :: i,j,l
      real*8 :: sward,eward
      call cs2llint_lluv(agrid,remap_veca2o,tauxa,tauya,tauxo,tauyo)
      do j=j_0,j_1
      do l=1,isp(j)
      do i=ifp(j,l),ilp(j,l)
        eward = +tauxo(i,j)
        sward = -tauyo(i,j)
        tauxo(i,j)= sward*coso(i,j)+eward*sino(i,j)
        tauyo(i,j)= eward*coso(i,j)-sward*sino(i,j)
      enddo
      enddo
      enddo
      return
      end subroutine veca2o

      subroutine veco2a(tauxo,tauyo,tauxa,tauya)
c --- mapping vector from ogcm C grid to agcm A grid
c --- input  tauxo/tauyo (N/m2): +i_/+j_ward on ogcm C grid
c ---                            (S_/E_ward in Mercator domain)
c --- output tauxa/tauya (N/m2): E_/N_ward on agcm A grid
c
      use domain_decomp_1d, only : halo_update
      implicit none
      real*8, dimension(iio        , J_0H: J_1H) :: tauxo,tauyo
      real*8, dimension(aI_0H:aI_1H,aJ_0H:aJ_1H) :: tauxa,tauya
      real*8, dimension(iio        , J_0H: J_1H) :: nward,eward
      integer :: i,j,l,jb
      real*8 :: sine
c --- rotate taux/tauy to latlon orientation at ocean A grid
      call halo_update(ogrid,tauyo)
      nward=0.
      eward=0.
      do j=j_0,j_1
      jb=j+1
      do l=1,isp(j)           
      do i=ifp(j,l),ilp(j,l)
        if (ip(i,j).eq.1) then
          sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j) ! why not == 1?
          nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .               -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
          eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .               +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
        endif
      enddo
      enddo
      enddo
      call xgridremap_ij(remap_ssto2a, eward, tauxa)
      call xgridremap_ij(remap_ssto2a, nward, tauya)
      return
      end subroutine veco2a

      end module hycom_cpler
