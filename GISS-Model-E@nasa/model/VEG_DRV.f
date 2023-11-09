#include "rundeck_opts.h"


      module veg_drv
!@sum veg_drv contains variables and routines for vegetation driver
!@auth I. Alienov, N. Kiang, Y. Kim

      implicit none
      private
      save

      public get_vdata, get_cropdata, get_laimaxdata, get_soil_C_total

      contains

      subroutine get_vdata(year, vdata, vegnames)
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds!, AM_I_ROOT
      use pario, only: par_open,par_close,read_dist_data,variable_exists
      use filemanager, only : file_exists
      !use vegetation, only : cond_scheme,vegCO2X_off
      !use veg_com
      !use model_com, only : jyear,focean
      !use ghy_com, only : fearth
      use ent_mod, only: N_COVERTYPES,COVER_SAND !ykim - use ent_const to accomodate GISS and Ent PFTs.
      use timestream_mod, only : init_stream,read_stream,timestream      implicit none
      integer, intent(in) :: year
      real*8, intent(out) :: vdata(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      character(len=*) :: vegnames(:)
      !---
      INTEGER :: J_1, J_0, J_1H, J_0H, I_1H, I_0H, I_1, I_0
      integer :: i, j, k, fid
      real*8 :: s
      !character(len=32) :: vegnames(N_COVERTYPES-N_OTHER)
      logical, save :: init = .false.
      integer :: day
      type(timestream), save :: VEGstream(N_COVERTYPES)
      logical, save :: have_veg_file = .false.
      logical, save :: is_timestream = .false.

      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

c**** read land surface parameters or use defaults
      ! vegnames should be set elsewhere.  Do this once the model
      ! actually uses Ent vegetation types as input.
cddd      vegnames = (/
cddd     &     'brightsoil     ','tundra         ','grass          ',
cddd     &     'shrub_and_grass','tree_and_grass ','deciduous      ',
cddd     &     'evergreen      ','rainforest     ','cultivation    ',
cddd     &     'darksoil       '
cddd     &     /)

      day = 1 ! to pass a required argument

      if (.not. init) then
        init = .true.
        have_veg_file = file_exists('VEG')
        if (have_veg_file) then
          fid = par_open(grid,'VEG','read')
          is_timestream = variable_exists(grid,fid,'time')
          call par_close(grid,fid)
          if (is_timestream) then
            if (year<0)
     &           call stop_model("get_vdata: year<0 in timestream",255)
            do k=1,size(vegnames)
              write(6,*) 'GET_VDATA initializing:', trim(vegnames(k))
              call init_stream(grid,VEGstream(k),'VEG',
     &             trim(vegnames(k)),
     &             0d0,1d30,'none',year,day)
            enddo
          endif ! is timestream
        endif ! VEG file exists
      endif

      vdata(:,:,:) = 0.d0 ! set to zero everything not present in the file
      if ( have_veg_file .and. is_timestream ) then
        if(year<0)call stop_model("get_vdata: year<0 in timestream",255)
        do k=1,size(vegnames)
          write(6,*) 'GET_VDATA reading:', trim(vegnames(k)),year
          call read_stream(grid,VEGstream(k),year,day,vdata(:,:,k))
        enddo
      else if ( have_veg_file ) then
        fid = par_open(grid,'VEG','read')
        do k=1,size(vegnames)
          call read_dist_data(grid,fid,trim(vegnames(k)),vdata(:,:,k))
        enddo
        call par_close(grid,fid)
      else
        vdata(:,:,COVER_SAND) = 1.d0 ! all bare soil if no input data available
      endif

c**** zero-out vdata(11) until it is properly read in
!!! already taken care of above
!      do k=N_COVERTYPES-N_OTHER+1, N_COVERTYPES
!        vdata(:,:,k) = 0.
!      end do

      ! make sure that veg fractions are reasonable
      do j=J_0,J_1
        do i=I_0,I_1
          !print *, i,j
          !print *, vdata(i,j,:) 
          do k=1,N_COVERTYPES
            ! get rid of unreasonably small fractions
            if ( vdata(i,j,k) < 1.d-4 ) vdata(i,j,k) = 0.d0
          enddo
          s = sum( vdata(i,j,:) )
          if ( s > .9d0 ) then
            vdata(i,j,:) = vdata(i,j,:)/s
          else if ( s < .1d0 ) then
            print *, "missing veg data at ",i,j,"assume bare soil"
            vdata(i,j,: ) = 0.d0
            vdata(i,j,COVER_SAND) = 1.d0
          else
            print *,i,j,s
            print *, vdata(i,j,:) 
            call stop_model("Incorrect data in VEG file",255)
          endif
        enddo
      enddo

      end subroutine get_vdata

cddd      module cropdata_mod
cddd      use timestream_mod, only : timestream
cddd      implicit none
cddd!@var CROPstream interface for reading and time-interpolating the crop file
cddd!@+   See usage notes in timestream_mod
cddd      type(timestream) :: CROPstream
cddd      logical :: have_crops_file
cddd      end module cropdata_mod
      subroutine get_cropdata(year, cropdata)
!@sum get_cropdata reads timeseries file for crop fraction and
!@+   interpolates to requested year.

      use domain_decomp_atm, only : grid
      use timestream_mod, only : init_stream,read_stream,timestream
      !use cropdata_mod
      use filemanager, only : file_exists
      implicit none
      integer, intent(in) :: year
      real*8 :: cropdata(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO)
c
      logical, save :: init = .false.
      integer :: day
      type(timestream), save :: CROPstream
      logical, save :: have_crops_file

      day = 1 ! to pass a required argument

      if (.not. init) then
        init = .true.
        have_crops_file = file_exists('CROPS')
        if(have_crops_file) then
          call init_stream(grid,CROPstream,'CROPS','crops',
     &         0d0,1d30,'none',year,day)
        endif
      endif

      cropdata = 0.d0 ! make sure halo is intialized
      if(have_crops_file) then
        call read_stream(grid,CROPstream,year,day,cropdata)
        !Set min to zero, since no land mask yet -nyk 1/22/08
        cropdata = max(0d0, cropdata)
      endif

      end subroutine get_cropdata



      subroutine get_laimaxdata(laimaxdata)
      use DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds!, AM_I_ROOT
      use DOMAIN_DECOMP_ATM, only : READT_PARALLEL
      use filemanager
      !use vegetation, only : cond_scheme,vegCO2X_off
      !use veg_com
      !use model_com, only : jyear,focean
      !use ghy_com, only : fearth
      use ent_mod, only: N_COVERTYPES, N_OTHER,COVER_SAND !ykim - use ent_const to accomodate GISS and Ent PFTs.
      implicit none
      real*8, intent(out) :: laimaxdata(
     &     grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,N_COVERTYPES)
      !---
      INTEGER :: J_1, J_0, J_1H, J_0H, I_1H, I_0H, I_1, I_0
      integer :: i, j, k, iu_laimax
      real*8 :: s

      CALL getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

c**** read land surface parameters or use defaults
      call openunit("LAImax",iu_LAIMAX,.true.,.true.)
      do k=1,N_COVERTYPES-N_OTHER
        CALL READT_PARALLEL
     *    (grid,iu_LAIMAX,NAMEUNIT(iu_LAIMAX),laimaxdata(:,:,K),1)
      end do
c**** zero-out laimaxdata(11) until it is properly read in
      do k=N_COVERTYPES-N_OTHER+1, N_COVERTYPES
        laimaxdata(:,:,k) = 0.
      end do
      call closeunit(iu_LAIMAX)

      end subroutine get_laimaxdata



      subroutine get_soil_C_total(ncasa, soil_C_total)
      use DOMAIN_DECOMP_ATM, only : GRID
      use pario, only : par_open,par_close,read_dist_data
      use filemanager, only : file_exists
      implicit none
      integer, intent(in) :: ncasa
      real*8 :: !,intent(out) ::
     &     soil_C_total(ncasa,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO)
      !---
      integer :: fid
      if(file_exists('SOILCARB_global')) then
        fid = par_open(grid,'SOILCARB_global','read')
        call read_dist_data(grid,fid,'soil_C_total',soil_C_total,jdim=3)
        call par_close(grid,fid)
      else
        soil_C_total = 0.
      endif
      end subroutine get_soil_C_total

      end module veg_drv
