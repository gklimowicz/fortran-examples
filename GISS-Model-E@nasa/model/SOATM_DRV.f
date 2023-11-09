#include "rundeck_opts.h"

      module core_data
      use TimeConstants_mod, only: DAYS_PER_YEAR
      implicit none
      save

      real*8, dimension(:,:), allocatable :: runoff

      real*8, dimension(:,:,:), allocatable ::
     &     swdn0,lwdn0,prec0,srfsal0,rsi,
     &     swdn1,lwdn1,prec1,srfsal1,
     &     swdn2,lwdn2,prec2,srfsal2

      real*8, dimension(:,:,:), allocatable ::
     &     psl0,ts0,qs0,us0,vs0,
     &     psl1,ts1,qs1,us1,vs1,
     &     psl2,ts2,qs2,us2,vs2

!@dbparam sss_restore_dt timescale (days) for surf salinity
!@+       relaxation back to observations
      real*8 :: sss_restore_dt=2.*DAYS_PER_YEAR
!@dbparam sss_restore_dtice timescale (days) for surf salinity
!@+       relaxation back to observations in the presence of sea ice
      real*8 :: sss_restore_dtice=30.

!@dbparam interannual_forcing==1 if using interannual atmospheriric
!@+       forcing files
      integer :: interannual_forcing=0
!@dbparam iaf_year_[start,end] first and last years of forcing history.
!@var iaf_year current forcing year
!@+        == jyear modulo 1+iaf_year_end-iaf_year_start
      integer :: iaf_year_start=1948, iaf_year_end=2007
      integer :: iaf_year_sv=-999, iaf_year

      end module core_data

      subroutine alloc_core_data
      use domain_decomp_atm, only : grid
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR
      use core_data
      use Dictionary_mod, only : sync_param
      implicit none
      integer :: j_0h,j_1h, i_0h,i_1h

      call sync_param("sss_restore_dt",sss_restore_dt)
      call sync_param("sss_restore_dtice",sss_restore_dtice)
      call sync_param("interannual_forcing",interannual_forcing)
      if(interannual_forcing==1) then
        call sync_param("iaf_year_start",iaf_year_start)
        call sync_param("iaf_year_end",iaf_year_end)
      endif

      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      allocate(swdn0(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))
      allocate(lwdn0(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))
      allocate(prec0(i_0h:i_1h,j_0h:j_1h,12))
      allocate(srfsal0(i_0h:i_1h,j_0h:j_1h,12))
      allocate(rsi(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))

      allocate(swdn1(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))
      allocate(lwdn1(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))
      allocate(prec1(i_0h:i_1h,j_0h:j_1h,12))
      allocate(srfsal1(i_0h:i_1h,j_0h:j_1h,12))

      allocate(swdn2(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))
      allocate(lwdn2(i_0h:i_1h,j_0h:j_1h,INT_DAYS_PER_YEAR))
      allocate(prec2(i_0h:i_1h,j_0h:j_1h,12))
      allocate(srfsal2(i_0h:i_1h,j_0h:j_1h,12))

      allocate(runoff(i_0h:i_1h,j_0h:j_1h))

      allocate(psl0(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(ts0(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(qs0(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(us0(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(vs0(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))

      allocate(psl1(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(ts1(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(qs1(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(us1(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(vs1(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))

      allocate(psl2(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(ts2(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(qs2(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(us2(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))
      allocate(vs2(i_0h:i_1h,j_0h:j_1h,4*INT_DAYS_PER_YEAR))

      return
      end subroutine alloc_core_data

      subroutine read_core_data
      use domain_decomp_atm, only : grid
      use fluxes, only : atmocn
      use core_data
      use pario, only : par_open,par_close,read_dist_data
      integer :: fid
      integer :: i,j,n
      real*8, dimension(12) :: f12,a12,b12,c12
      real*8, dimension(365) :: f365,a365,b365,c365
      real*8, dimension(1460) :: f1460,a1460,b1460,c1460
c read river discharge
      fid = par_open(grid,'RUNOFF','read')
      call read_dist_data(grid, fid, 'Foxx_o_roff', runoff)
      call par_close(grid,fid)
c read sea surface salinity
      fid = par_open(grid,'SSS','read')
      call read_dist_data(grid, fid, 'SALT', srfsal0)
      call par_close(grid,fid)
c read sea ice fraction
      fid = par_open(grid,'RSI','read')
      call read_dist_data(grid, fid, 'rsi', rsi)
      call par_close(grid,fid)
c
c create parabolic coeffs for each time interval
c
      do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          if(atmocn%focean(i,j).le.0.) cycle
c
          n = 12
          f12 = srfsal0(i,j,:)
          call coeffs1d_pos(f12,a12,b12,c12,n,n)
          srfsal0(i,j,:) = a12
          srfsal1(i,j,:) = b12
          srfsal2(i,j,:) = c12
c
        enddo
      enddo

c
c read normal-year files if appropriate
c
      if(interannual_forcing==0) then
c read precip
      fid = par_open(grid,'PREC','read')
      call read_dist_data(grid, fid, 'PREC_RAW', prec0)
      call par_close(grid,fid)
c read radiation
      fid = par_open(grid,'RAD','read')
      call read_dist_data(grid, fid, 'SWDN', swdn0)
      call read_dist_data(grid, fid, 'LWDN', lwdn0)
      call par_close(grid,fid)
c read SLP and 10-m u,v,t,q
      fid = par_open(grid,'SLP','read')
      call read_dist_data(grid, fid, 'SLP', psl0)
      call par_close(grid,fid)
      fid = par_open(grid,'T10','read')
      call read_dist_data(grid, fid, 'T_10', ts0)
      call par_close(grid,fid)
      fid = par_open(grid,'Q10','read')
      call read_dist_data(grid, fid, 'Q_10', qs0)
      call par_close(grid,fid)
      fid = par_open(grid,'U10','read')
      call read_dist_data(grid, fid, 'U_10', us0)
      call par_close(grid,fid)
      fid = par_open(grid,'V10','read')
      call read_dist_data(grid, fid, 'V_10', vs0)
      call par_close(grid,fid)
c
c create parabolic coeffs for each time interval
c
      do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          if(atmocn%focean(i,j).le.0.) cycle
c
          n = 12
          f12 = prec0(i,j,:)
          call coeffs1d_pos(f12,a12,b12,c12,n,n)
          prec0(i,j,:) = a12
          prec1(i,j,:) = b12
          prec2(i,j,:) = c12
c
          n = 365
          f365 = lwdn0(i,j,:)
          call coeffs1d_pos(f365,a365,b365,c365,n,n)
          lwdn0(i,j,:) = a365
          lwdn1(i,j,:) = b365
          lwdn2(i,j,:) = c365
c
          n = 365
          f365 = swdn0(i,j,:)
          call coeffs1d_pos(f365,a365,b365,c365,n,n)
          swdn0(i,j,:) = a365
          swdn1(i,j,:) = b365
          swdn2(i,j,:) = c365
c
          n = 1460
          f1460 = psl0(i,j,:)
          call coeffs1d_pos(f1460,a1460,b1460,c1460,n,n)
          psl0(i,j,:) = a1460
          psl1(i,j,:) = b1460
          psl2(i,j,:) = c1460
c
          n = 1460
          f1460 = ts0(i,j,:)
          call coeffs1d_pos(f1460,a1460,b1460,c1460,n,n)
          ts0(i,j,:) = a1460
          ts1(i,j,:) = b1460
          ts2(i,j,:) = c1460
c
          n = 1460
          f1460 = qs0(i,j,:)
          call coeffs1d_pos(f1460,a1460,b1460,c1460,n,n)
          qs0(i,j,:) = a1460
          qs1(i,j,:) = b1460
          qs2(i,j,:) = c1460
c
          n = 1460
          f1460 = us0(i,j,:)
          call coeffs1d(f1460,a1460,b1460,c1460,n,n)
          us0(i,j,:) = a1460
          us1(i,j,:) = b1460
          us2(i,j,:) = c1460
c
          n = 1460
          f1460 = vs0(i,j,:)
          call coeffs1d(f1460,a1460,b1460,c1460,n,n)
          vs0(i,j,:) = a1460
          vs1(i,j,:) = b1460
          vs2(i,j,:) = c1460
c
        enddo
      enddo
      endif ! not interannual forcing
      return
      end subroutine read_core_data

      subroutine read_core_interannual_data
c assumptions regarding input files:
c   1. forcing history stored with separate files for each year
c   2. files are padded with four extra timesteps (two at the
c      beginning of the year, two at the end) to facilitate
c      time interpolation.
      use domain_decomp_atm, only : grid
      use fluxes, only : atmocn
      use core_data
      use pario, only : par_open,par_close,read_dist_data
      integer :: fid
      integer :: i,j,n,n2
      integer :: i_0h,i_1h,j_0h,j_1h
      real*8, dimension(16) :: f16,a16,b16,c16
      real*8, dimension(369) :: f369,a369,b369,c369
      real*8, dimension(1464) :: f1464,a1464,b1464,c1464
c
      real*8, dimension(:,:,:), allocatable ::
     &     swdn0_,lwdn0_,prec0_,psl0_,ts0_,qs0_,us0_,vs0_
c
      character(len=4) :: cyr

c
      write(cyr,'(i4)') iaf_year
c
      i_0h = grid%i_strt_halo
      i_1h = grid%i_stop_halo
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      allocate(prec0_(i_0h:i_1h,j_0h:j_1h,16))
      allocate(swdn0_(i_0h:i_1h,j_0h:j_1h,369))
      allocate(lwdn0_(i_0h:i_1h,j_0h:j_1h,369))
      allocate(psl0_(i_0h:i_1h,j_0h:j_1h,1464))
      allocate(ts0_(i_0h:i_1h,j_0h:j_1h,1464))
      allocate(qs0_(i_0h:i_1h,j_0h:j_1h,1464))
      allocate(us0_(i_0h:i_1h,j_0h:j_1h,1464))
      allocate(vs0_(i_0h:i_1h,j_0h:j_1h,1464))

c read precip
      fid = par_open(grid,'IAF/precip.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'PREC_RAW', prec0_)
      call par_close(grid,fid)
c read radiation
      fid = par_open(grid,'IAF/rad.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'SWDN', swdn0_)
      call read_dist_data(grid, fid, 'LWDN', lwdn0_)
      call par_close(grid,fid)
c read SLP and 10-m u,v,t,q
      fid = par_open(grid,'IAF/slp.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'SLP', psl0_)
      call par_close(grid,fid)
      fid = par_open(grid,'IAF/t.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'T_10', ts0_)
      call par_close(grid,fid)
      fid = par_open(grid,'IAF/q.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'Q_10', qs0_)
      call par_close(grid,fid)
      fid = par_open(grid,'IAF/u.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'U_10', us0_)
      call par_close(grid,fid)
      fid = par_open(grid,'IAF/v.'//cyr//'.nc','read')
      call read_dist_data(grid, fid, 'V_10', vs0_)
      call par_close(grid,fid)

c
c create parabolic coeffs for each time interval
c
      do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,grid%i_stop
          if(atmocn%focean(i,j).le.0.) cycle
c
          n = 16; n2 = 12
          f16 = prec0_(i,j,:)
          call coeffs1d_pos(f16,a16,b16,c16,n,n2)
          prec0(i,j,:) = a16(3:14)
          prec1(i,j,:) = b16(3:14)
          prec2(i,j,:) = c16(3:14)
c
          n = 369; n2 = 365
          f369 = lwdn0_(i,j,:)
          call coeffs1d_pos(f369,a369,b369,c369,n,n2)
          lwdn0(i,j,:) = a369(3:367)
          lwdn1(i,j,:) = b369(3:367)
          lwdn2(i,j,:) = c369(3:367)
c
          n = 369; n2 = 365
          f369 = swdn0_(i,j,:)
          call coeffs1d_pos(f369,a369,b369,c369,n,n2)
          swdn0(i,j,:) = a369(3:367)
          swdn1(i,j,:) = b369(3:367)
          swdn2(i,j,:) = c369(3:367)
c
          n = 1464; n2 = 1460
          f1464 = psl0_(i,j,:)
          call coeffs1d_pos(f1464,a1464,b1464,c1464,n,n2)
          psl0(i,j,:) = a1464(3:1462)
          psl1(i,j,:) = b1464(3:1462)
          psl2(i,j,:) = c1464(3:1462)
c
          n = 1464; n2 = 1460
          f1464 = ts0_(i,j,:)
          call coeffs1d(f1464,a1464,b1464,c1464,n,n2)
          ts0(i,j,:) = a1464(3:1462)
          ts1(i,j,:) = b1464(3:1462)
          ts2(i,j,:) = c1464(3:1462)
c
          n = 1464; n2 = 1460
          f1464 = qs0_(i,j,:)
          call coeffs1d_pos(f1464,a1464,b1464,c1464,n,n2)
          qs0(i,j,:) = a1464(3:1462)
          qs1(i,j,:) = b1464(3:1462)
          qs2(i,j,:) = c1464(3:1462)
c
          n = 1464; n2 = 1460
          f1464 = us0_(i,j,:)
          call coeffs1d(f1464,a1464,b1464,c1464,n,n2)
          us0(i,j,:) = a1464(3:1462)
          us1(i,j,:) = b1464(3:1462)
          us2(i,j,:) = c1464(3:1462)
c
          n = 1464; n2 = 1460
          f1464 = vs0_(i,j,:)
          call coeffs1d(f1464,a1464,b1464,c1464,n,n2)
          vs0(i,j,:) = a1464(3:1462)
          vs1(i,j,:) = b1464(3:1462)
          vs2(i,j,:) = c1464(3:1462)
c
        enddo
      enddo

      deallocate(swdn0_,lwdn0_,prec0_,psl0_,ts0_,qs0_,us0_,vs0_)

      return
      end subroutine read_core_interannual_data

      subroutine coeffs1d(fm,c0,c1,c2,n,ndx)
      implicit none
      integer :: n,ndx
      real*8, parameter :: twoby3=2d0/3d0
      real*8, parameter :: c712=7d0/12d0,c112=1d0/12d0
      real*8, dimension(n) :: fm,c0,c1,c2
      integer :: i
      real*8 :: xl,xr,fl,fr,dx
      real*8, dimension(:), allocatable :: fm_
      allocate(fm_(0:n+2))
      fm_(0) = fm(n)
      fm_(1:n) = fm(1:n)
      fm_(n+1:n+2) = fm(1:2)
      dx = 1d0/real(ndx,kind=8)
      xl = -dx*real((n-ndx)/2,kind=8)
      xr = xl
      i = n
      fl = c712*(fm_(i)+fm_(i+1))-c112*(fm_(i-1)+fm_(i+2))
      do i=1,n
        xr = xr + dx
        fr = c712*(fm_(i)+fm_(i+1))-c112*(fm_(i-1)+fm_(i+2))
        c2(i) = (fr+fl-2.*fm(i))/(xl*xl+xr*xr-twoby3*(xr**3-xl**3)/dx)
        c1(i) = (fr-fl)/dx-c2(i)*(xr+xl)
        c0(i) = fl-xl*(c1(i)+c2(i)*xl)
        fl = fr
        xl = xr
      enddo
      deallocate(fm_)
      return
      end subroutine coeffs1d

      subroutine coeffs1d_pos(fm,c0,c1,c2,n,ndx)
      implicit none
      integer :: n,ndx
      real*8, parameter :: twoby3=2d0/3d0,by3=1d0/3d0
      real*8, parameter :: c712=7d0/12d0,c112=1d0/12d0
c      real*8, parameter :: c712=6d0/12d0,c112=0d0/12d0
      real*8, dimension(n) :: fm,c0,c1,c2
      integer :: i
      real*8 :: xl,xr,fl,fr,dx,zmin,fmin,fl_,fr_
      real*8, dimension(:), allocatable :: fm_
      real*8 :: b0,b1,b2,momrat
      allocate(fm_(0:n+2))
      fm_(0) = fm(n)
      fm_(1:n) = fm(1:n)
      fm_(n+1:n+2) = fm(1:2)
      dx = 1d0/real(ndx,kind=8)
      xl = -dx*real((n-ndx)/2,kind=8)
      xr = xl
      i = n
      fl = c712*(fm_(i)+fm_(i+1))-c112*(fm_(i-1)+fm_(i+2))
      do i=1,n
        xr = xr + dx
        fr = c712*(fm_(i)+fm_(i+1))-c112*(fm_(i-1)+fm_(i+2))
        fl_ = max(fl,0d0)
        fr_ = max(fr,0d0)
        if(fm_(i).le.0.) then
          fl_ = 0.
          fr_ = 0.
        endif
        if(fm_(i-1).le.0.) fl_ = 0.
        if(fm_(i+1).le.0.) fr_ = 0.
        b2 = .75*(fr_+fl_-2.*fm(i))
        if(b2.ne.0d0) then
          b1 = .5*(fr_-fl_)
          zmin = -.5d0*b1/b2
          fmin = (fm(i)-b2*by3) + .5d0*b1*zmin ! b0 in parens
          if(abs(zmin).le.1d0 .and. fmin.lt.0.) then
            momrat = (fm(i)/(by3+zmin*zmin))/b2
            fl_ = fm(i)+momrat*(fl_-fm(i))
            fr_ = fm(i)+momrat*(fr_-fm(i))
          endif
        endif
        if(fm_(i).le.0.) then
          fl_ = 0.
          fr_ = 0.
        endif
        if(fm_(i-1).le.0.) fl_ = max(fl_,0.)
        if(fm_(i+1).le.0.) fr_ = max(fr_,0.)
        c2(i) = (fr_+fl_-2.*fm(i))/(xl*xl+xr*xr-twoby3*(xr**3-xl**3)/dx)
        c1(i) = (fr_-fl_)/dx-c2(i)*(xr+xl)
        c0(i) = fl_-xl*(c1(i)+c2(i)*xl)
        fl = fr
        xl = xr
      enddo
      deallocate(fm_)
      return
      end subroutine coeffs1d_pos

      subroutine get_ocean_forcings
      use model_com, only : dtsrc,nday,itime,modelEclock
      use constant, only : lhm,tf
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      use domain_decomp_atm, only : grid
      use geom, only : axyp
      use fluxes, only : atmocn,atmice
      use seaice_com, only : si_ocn
      use core_data
      use domain_decomp_1d, only : globalsum,hasNorthPole
      implicit none
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: evx,prx
      real*8 :: oevg,oprg,ebyp,t,rsiloc,sssresfac,sssresfac_ice
      integer :: i,j,l,n,j6hr,itmod,itperyr,jmm
      integer :: i_0,i_1, j_0,j_1
      integer :: jyear,jday

      jyear = modelEclock%getYear()
      jday  = modelEclock%getDayOfYear()

      if(interannual_forcing==1) then
        if(jyear < iaf_year_start) then
          write(6,*) 'please set model year >= iaf_year_start'
          call stop_model('model year < iaf_year_start',255)
        endif
        iaf_year = iaf_year_start +
     &       mod(jyear-iaf_year_start,1+iaf_year_end-iaf_year_start)
        if(iaf_year.ne.iaf_year_sv) then
          call read_core_interannual_data        
          iaf_year_sv = iaf_year
        endif
      endif

      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      J_0 = grid%J_STRT
      J_1 = grid%J_STOP
      sssresfac     = exp(-dtsrc/(sss_restore_dt   *SECONDS_PER_DAY))
      sssresfac_ice = exp(-dtsrc/(sss_restore_dtice*SECONDS_PER_DAY))
      itperyr = nday*INT_DAYS_PER_YEAR
      itmod = mod(itime,itperyr)
      j6hr = 1+itmod*4/nday
      t = (real(itmod,kind=8)+.5d0)/real(itperyr,kind=8)
      jmm = 1+int(12.*t)
      do j=j_0,j_1
      do i=i_0,i_1
        if(atmocn%focean(i,j).eq.0.) cycle
        atmocn%PREC(i,j) = dtsrc*(
     &       prec0(i,j,jmm)+t*(prec1(i,j,jmm)+t*prec2(i,j,jmm))
     &       )
        atmocn%FLOWO(i,j) = runoff(i,j)*dtsrc
        atmocn%EFLOWO(i,j) = 0. ! for now
        atmocn%fshort(i,j) = (
     &       swdn0(i,j,jday)+t*(swdn1(i,j,jday)+t*swdn2(i,j,jday)) )
        atmocn%flong(i,j) = (
     &       lwdn0(i,j,jday)+t*(lwdn1(i,j,jday)+t*lwdn2(i,j,jday)) )
        atmocn%USAVG(i,j) = (
     &       us0(i,j,j6hr)+t*(us1(i,j,j6hr)+t*us2(i,j,j6hr)) )
        atmocn%VSAVG(i,j) = (
     &       vs0(i,j,j6hr)+t*(vs1(i,j,j6hr)+t*vs2(i,j,j6hr)) )
        atmocn%TSAVG(i,j) = (
     &       ts0(i,j,j6hr)+t*(ts1(i,j,j6hr)+t*ts2(i,j,j6hr)) )
        atmocn%QSAVG(i,j) = (
     &       qs0(i,j,j6hr)+t*(qs1(i,j,j6hr)+t*qs2(i,j,j6hr)) )
        atmocn%SRFP(i,j) = .01d0*(
     &       psl0(i,j,j6hr)+t*(psl1(i,j,j6hr)+t*psl2(i,j,j6hr)) )
        if(atmocn%TSAVG(i,j).gt.tf) then
          atmocn%EPREC(i,j) = 0.
        else
          atmocn%EPREC(i,j) = -lhm*atmocn%PREC(i,j)
        endif
        atmocn%sssobs(i,j) = .001d0*(
     &       srfsal0(i,j,jmm)+t*(srfsal1(i,j,jmm)+t*srfsal2(i,j,jmm))
     &       )
        !atmocn%rsiobs(i,j) = si_ocn%rsi(i,j)
        atmocn%rsiobs(i,j) = rsi(i,j,jday)
        rsiloc = atmocn%rsiobs(i,j)
c        if(rsi.gt.0.) then
c          atmocn%sssresfac(i,j) = sssresfac_ice
c        else
c          atmocn%sssresfac(i,j) = sssresfac
c        endif
        atmocn%sssresfac(i,j) =
     &       rsiloc*sssresfac_ice+(1.-rsiloc)*sssresfac
      enddo
      enddo

c scale prec,runoff so that global evap = P+R
      do j=j_0,j_1
      do i=i_0,i_1
        prx(i,j) = 0.
        evx(i,j) = 0.
        if(atmocn%focean(i,j).eq.0.) cycle
        prx(i,j) = axyp(i,j)*(atmocn%prec(i,j)+atmocn%flowo(i,j))
        rsiloc = si_ocn%rsi(i,j)
        evx(i,j) = axyp(i,j)*((1.-rsiloc)*atmocn%evapor(i,j)
     &                           +rsiloc *atmice%evapor(i,j))
      enddo
      enddo
      if(hasNorthPole(grid)) then
        prx(2:i_1,j_1)=prx(1,j_1)
        evx(2:i_1,j_1)=evx(1,j_1)
      endif
      call globalsum(grid, prx, oprg, all=.true.)
      call globalsum(grid, evx, oevg, all=.true.)
      ebyp = oevg/oprg
      if(grid%gid.eq.0) write(6,*) 'oevg/oprg ',ebyp
      do j=j_0,j_1
      do i=i_0,i_1
        if(atmocn%focean(i,j).eq.0.) cycle
        atmocn%prec(i,j) = atmocn%prec(i,j)*ebyp
        atmocn%flowo(i,j) = atmocn%flowo(i,j)*ebyp
      enddo
      enddo
      return
      end subroutine get_ocean_forcings

      MODULE SURF_ALBEDO
!@sum SURF_ALBEDO contains parameters/variables needed for albedo calc
!@auth A. Lacis/V. Oinas (modifications by I. Aleinov/G. Schmidt)
      implicit none
      save
      private

      public get_surf_albedo

      real*8, dimension(6), parameter :: fr6band=
     &    (/ 0.53d0, 0.06d0, 0.20d0, 0.07d0, 0.08d0, 0.06d0 /)

!@var SRFOAM look up table for ocean foam as a function of wind speed
      REAL*8, PARAMETER :: SRFOAM(25) = (/
     *     0.000,0.000,0.000,0.000,0.001,0.002,0.003,0.005,0.007,0.010,
     *     0.014,0.019,0.025,0.032,0.041,0.051,0.063,0.077,0.094,0.112,
     *     0.138,0.164,0.191,0.218,0.246/)

C

!@var ASHZOI,ANHZOI hemisph.Ice Albedo half-max depth (m) (orig.version)
      REAL*8 :: ASHZOI=.1d0, ANHZOI=.1d0             ! tuning parameters
!@var DMOICE  masking depth for snow on sea ice           (orig.version)
      REAL*8 :: DMOICE = 10.                         ! tuning parameter
!@var AVSCAT,ANSCAT,AVFOAM,ANFOAM for ocean albedo calc
      REAL*8 ::                                      ! tuning parameters
     *     AVSCAT=.0156d0, ANSCAT=0d0, AVFOAM=.2197d0, ANFOAM=.1514d0

!@var ASNALB snow albedo for old snow
!@var AOIALB seaice albedo                            (original version)
      REAL*8, parameter ::
C                        VIS  NIR1  NIR2  NIR3  NIR4  NIR5    NIR
     *     ASNALB(7)=(/.60d0,.55d0,.55d0,.30d0,.10d0,.05d0, .35d0/),
     *     AOIALB(7)=(/.55d0,.50d0,.45d0,.25d0,.10d0,.05d0, .30d0/)

C**** variables that control snow aging calculation (over land)
!@var AGEXPF exponent in snowage calculation depends on hemi/surf type
!@var ALBDIF difference in albedo as function of snowage
      REAL*8 ::
     *     AGEXPF(3,2) = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0,  0.2d0 /), (/3,2/) ),
     *     ALBDIF(3,2) = RESHAPE( (/
C          SH EA   SH OC   SH LI   NH EA   NH OC   NH LI
     *     0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0, 0.35d0/), (/3,2/) )

C**** parameters used for Schramm sea ice albedo scheme (Hansen)
!@var AOImin,AOImax           range for seaice albedo
!@var ASNwet,ASNdry           wet,dry snow albedo over sea ice
!@var AMPmin                  mininimal melt pond albedo
      REAL*8 ::
C                         VIS   NIR1   NIR2   NIR3   NIR4   NIR5
     *     AOImin(6)=(/ .05d0, .05d0, .05d0, .050d0, .05d0, .03d0/),
     *     AOImax(6)=(/ .62d0, .42d0, .30d0, .120d0, .05d0, .03d0/),
     *     ASNwet(6)=(/ .85d0, .75d0, .50d0, .175d0, .03d0, .01d0/),
     *     ASNdry(6)=(/ .90d0, .85d0, .65d0, .450d0, .10d0, .10d0/),
     *     AMPmin(6)=(/ .10d0, .05d0, .05d0, .050d0, .05d0, .03d0/)

!@var GZSNOW asymmetry parameter for snow over sea ice
!@+   from Wiscombe and Warren (1980) JAS
!@+   Note this is used for ice + melt-ponds as well.
      REAL*8, PARAMETER :: GZSNOW(7,2) = RESHAPE( (/
C       VIS     NIR1    NIR2     NIR3     NIR4     NIR5    NIRT
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0,! SH
     * 0.95d0, 0.94d0, 0.905d0, 0.896d0, 0.894d0, 0.89d0, 0.91d0 ! NH
     *     /), (/7,2/) )

      REAL*8, public :: GTAU(51,11,143),TGDATA(122,13)

      contains

      SUBROUTINE GET_SURF_ALBEDO(
     i     COSZ,
     i     HEMI,
     i     AGESN,POCEAN,POICE,
     i     ZOICE,FMP,ZSNWOI,ZMP,
     i     SNOWOI,WMAG,dalbsn,
     i     FLAGS,
     o     OCNALB,SIALB
     &     )
!@sum GETSUR computes surface albedo for each grid box
!@auth A. Lacis/V. Oinas (modifications by I. Aleinov/G. Schmidt)


      implicit none

!********* start  in/out *****************************
C**** config data
      REAL*8, parameter :: SNOAGE_FAC_MAX=.5d0
      INTEGER, parameter :: KSIALB=0,KZSNOW=1
C**** inputs
      REAL*8 COSZ
      INTEGER HEMI
      REAL*8 AGESN,POCEAN,POICE,
     *     ZOICE,FMP,ZSNWOI,ZMP,
     *     SNOWOI,WMAG,
     &     dalbsn
      LOGICAL*4 :: FLAGS
C**** outputs
      REAL*8 OCNALB,SIALB
!********* end in/out *****************************
      INTEGER K,J,L,JH,IWM,JWM
      REAL*8 ASNAGE,FSNAGE,AVSCUM,ANSCUM,WMJ,WMI,FRFOAM
     *   ,KKZSNO,FDZICE,AHMZOI

C**** local variables for albedos for each surface type
      REAL*8 XOCVIS,XOCNIR,EXPSNO

C**** variables used for sea ice albedo calculation (6 bands, Hansen)
      REAL*8, dimension(6) :: BOIVN,BSNVN,XOIVN,XSNVN,almp6,alsf6
      REAL*8 :: patchy,snagfac

C     -----------------------------------------------------------------
C     Ocean Albedo Dependence on Zenith Angle and Wind Speed
C
      REAL*8 XVH2O, WMAG1, X
      XVH2O(WMAG1,X)=.021D0+X*X*(.0421D0+X*(.1283D0+X*(-.04D0+X*(3.117D0
     +              /(5.679D0+WMAG1)+X*.025D0/(.3333D0+WMAG1)))))
C     -----------------------------------------------------------------
C
C-----------------------------------------------------------------------
C     Select albedo computation using KSIALB
C     KSIALB= 0  Schramm oi.alb, Antarc/Greenl alb=.8 (J.Hansen)
C     KSIALB= 1  6-band original albedo - no 'fixups' (Andy Lacis)
C       else     Schramm oi.alb but no land ice 'fixups'


C           Get Albedo, Thermal Flux, Flux Derivative for each Surf Type
C           ------------------------------------------------------------

      JH=HEMI
      KKZSNO=KZSNOW
      IF(COSZ < 0.001) KKZSNO=0
C
      EXPSNO=1.D0
C
      AVSCUM=0.D0
      ANSCUM=0.D0
C
      BOIVN=0. ; BSNVN=0.
      XOIVN=0. ; XSNVN=0.

C
C                                             --------------------------
C                                             Ocean Albedo Specification
C                                             --------------------------
C
      IF(POCEAN > 0.) THEN
        X=0.5D0+(0.5D0-COSZ)
        XOCVIS=XVH2O(WMAG,X)+AVSCAT+AVSCUM
        XOCNIR=XVH2O(WMAG,X)+ANSCAT+ANSCUM
C
        IWM=WMAG
        IF(IWM < 1) IWM=1
        IF(IWM > 24) IWM=24
        JWM=IWM+1
        WMJ=WMAG-IWM
        WMI=1.D0-WMJ
        FRFOAM=WMI*SRFOAM(IWM)+WMJ*SRFOAM(JWM)
C
        XOCVIS=XOCVIS*(1.D0-FRFOAM)+FRFOAM*AVFOAM
        XOCNIR=XOCNIR*(1.D0-FRFOAM)+FRFOAM*ANFOAM
        OCNALB = FR6BAND(1)*XOCVIS + (1.-FR6BAND(1))*XOCNIR
      ENDIF



C
C                                         ------------------------------
C                                         Ocean Ice Albedo Specification
C                                         ------------------------------
      IF(POICE > 0.) THEN
      IF(KSIALB==1) THEN
C****                      original version (6-band)
        EXPSNO=EXP(-SNOWOI/DMOICE)
C**** Set snow albedo over sea ice
        ASNAGE=ALBDIF(2,JH)*EXP(-AGEXPF(2,JH)*AGESN)
        DO L=1,6
          FSNAGE=1.D0
          IF(L > 2) FSNAGE=2.0D0/L
          BSNVN(L)=ASNALB(L)+ASNAGE*FSNAGE
        END DO

C**** set ice albedo
        AHMZOI=ASHZOI
        IF(HEMI==2) AHMZOI=ANHZOI
        FDZICE=ZOICE/(ZOICE+AHMZOI)   ! ZOICE = ice depth (m)
        DO L=1,6
          BOIVN(L)=FDZICE*AOIALB(L)*EXPSNO+BSNVN(L)*(1.D0-EXPSNO)
C**** Set zenith angle dependence if required
          IF (KKZSNO > 0) THEN
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,JH),XOIVN(L))
          ELSE
            XOIVN(L)=BOIVN(L)
          END IF
        END DO
C**** end of original version (6-band)

      else        ! KSIALB not 1
C**** Schramm/J. Hansen's sea ice albedo formulas (6 spectral bands)
C**** Bare ice:
        BOIVN(1:6) = aoimax(1:6)
        if(ZOICE < 1.)then       ! ZOICE: ice depth (at least Z1I=.1m)
          BOIVN(1:4)=aoimin(1:4)+(aoimax(1:4)-aoimin(1:4))*sqrt(ZOICE)
        endif
C**** Snow:    patchy: snow_cover_fraction (<1 if snow depth < .1m)
        patchy = 0.       !  max(0.d0 , min(1.d0, 10.d0*ZSNWOI) )
        if(ZSNWOI > 0.)then
          if(ZSNWOI >= 0.1d0)then      ! snow deeper than .1m
            patchy=1d0
          else
            patchy=ZSNWOI/0.1d0
          endif
          if(FLAGS)then         ! wet snow
            alsf6(1:6)=asnwet(1:6)
          else                  ! dry snow
            alsf6(1:6)=asndry(1:6)
          endif
C       snow aging based on Loth and Graf (1998)
C       Dry, Wet(thick), Wet(thin) snow decreases by
C       0.006,  0.015 and 0.071 per day, respectively (for mean)
C       assume decrease for each band is proportional
          if (FLAGS) then
            if (ZSNWOI > 0.25) then
              snagfac = 0.015d0/0.7d0 * AGESN
            else
              snagfac = 0.071d0/0.7d0 * AGESN
            end if
          else
            snagfac = 0.006d0/0.82d0  * AGESN
          end if
C       make sure snow albedo doesn't get too low!
          snagfac=min(SNOAGE_FAC_MAX,snagfac)
          alsf6(1:6)=alsf6(1:6)*(1.-snagfac)
C       combine bare ice and snow albedos
          BOIVN(1:6)=BOIVN(1:6)*(1.-patchy)+alsf6(1:6)*patchy
        endif
C**** Melt ponds:
        almp6(1:6)=ampmin(1:6)    ! used if melt pond deeper than .5m
        if(ZMP > 0. .and. ZMP < 0.5)then
          almp6(1:6)=almp6(1:6)+(aoimax(1:6)-almp6(1:6))*(1.-2.*ZMP)**2
        end if
c**** combined sea ice albedo
        BOIVN(1:6)=BOIVN(1:6)*(1.-FMP)+almp6(1:6)*FMP
C**** set zenith angle dependence
        IF (KKZSNO > 0) THEN     ! for all surface types
          DO L=1,6
            CALL RXSNOW(BOIVN(L),COSZ,GZSNOW(L,JH),XOIVN(L))
          END DO
        ELSE
          XOIVN(1:6)=BOIVN(1:6)
        END IF
        EXPSNO=1.-patchy

c**** Reduce the Ocean Ice albedo by dalbsn
        DO L=1,2
          BOIVN(L) = max(0.d0,BOIVN(L)+dalbsn/L)
          XOIVN(L) = max(0.d0,XOIVN(L)+dalbsn/L)
        END DO
      end if  !  KSIALB not 1:  Schramm/Hansen
      SIALB = sum(fr6band*xoivn)
C
      ENDIF

      RETURN
      END SUBROUTINE GET_SURF_ALBEDO

      END MODULE SURF_ALBEDO

      SUBROUTINE RXSNOW(RBSNO,XCOSZ,GGSNO,RXSNO)
!@sum RXSNOW calculate zenith angle dependence for snow/ice albedo
!@auth A. Lacis (modified by G. Schmidt)
  !    USE RADPAR, only : gtsalb,sgpgxg
      IMPLICIT NONE
!@var RBSNO diffuse albedo
      REAL*8, INTENT(IN) :: RBSNO
!@var XCOSZ zenith angle
      REAL*8, INTENT(IN) :: XCOSZ
!@var GGSNO Asymmetry parameter for snow
      REAL*8, INTENT(IN) :: GGSNO
!@var RXSNO direct albedo
      REAL*8, INTENT(OUT) :: RXSNO
      INTEGER NDBLS,NN
      REAL*8 XXG,XXT,GGSN,RBSN,FRTOP,TAU,TAUSN,GPFF,PR,PT,DBLS,SECZ,XANB
     *     ,XANX,TANB,TANX,RASB,RASX,BNORM,XNORM,RARB,RARX,XATB,DENOM,DB
     *     ,DX,UB,UX,DRBRAT,RBBOUT

      IF(RBSNO < 0.05D0) THEN
        RXSNO=RBSNO
        RETURN
      ENDIF
      XXG=0.D0
      XXT=0.D0
      GGSN=GGSNO
      IF(GGSNO > 0.9D0) GGSN=0.9D0
      RBSN=RBSNO
      FRTOP=1.D0
      IF(RBSNO > 0.5D0) THEN
        RBSN=0.5D0
        FRTOP=((1.D0-RBSNO)/0.5D0)**2
      ENDIF

      CALL GTSALB(XXG,XXT,RBBOUT,RBSN,GGSN,TAUSN,2)
      CALL SGPGXG(XCOSZ,TAUSN,GGSN,GPFF)
      PR=1.D0-GPFF
      PT=1.D0+GPFF
      DBLS=10.D0+1.44269D0*LOG(TAUSN)
      NDBLS=DBLS
      TAU=TAUSN/2**NDBLS
C     Set optically thin limit values of R,T,X using PI0 renormalization
C     ------------------------------------------------------------------
C
      SECZ=1.D0/XCOSZ
      XANB=EXP(-TAU-TAU)
      XANX=EXP(-TAU*SECZ)
      TANB=PT*XANB
      XXT=(SECZ-2.D0)*TAU
      TANX=PT*SECZ
     +    *(.5D0+XXT*(.25D0+XXT*(.0833333D0+XXT*(.0208333D0+XXT))))*XANX
      RASB=PR*(1.D0-TAU*(2.D0-2.66667D0*TAU*(1.D0-TAU)))
      XXT=(SECZ+2.D0)*TAU
      RASX=PR*SECZ
     +    *(.5D0-XXT*(.25D0-XXT*(.0833333D0-XXT*(.0208333D0-XXT))))
      BNORM=(1.D0-XANB)/(RASB+TANB)
      XNORM=(1.D0-XANX)/(RASX+TANX)
      RASB=RASB*BNORM
      RASX=RASX*XNORM
      TANB=TANB*BNORM
      TANX=TANX*XNORM
      DO NN=1,NDBLS
        RARB=RASB*RASB
        RARX=XANX*RASX
        XATB=XANB+TANB
        DENOM=1.D0-RARB
        DB=(TANB+XANB*RARB)/DENOM
        DX=(TANX+RARX*RASB)/DENOM
        UB=RASB*(XANB+DB)
        UX=RARX+RASB*DX
        RASB=RASB+XATB*UB
        RASX=RASX+XATB*UX
        TANB=XANB*TANB+XATB*DB
        TANX=XANX*TANX+XATB*DX
        XANB=XANB*XANB
        XANX=XANX*XANX
      END DO
      DRBRAT=RASX/RBSN-1.D0
      RXSNO=RBSNO*(1.D0+DRBRAT*FRTOP)
      RETURN
      END SUBROUTINE RXSNOW


      SUBROUTINE SETGTS
CCC   SUBROUTINE GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)
      USE SURF_ALBEDO, only: TGDATA
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: GIN,TAUIN,RBBIN,EGIN
      INTEGER, INTENT(IN) :: KGTAUR
      REAL*8, INTENT(OUT) :: RBBOUT,TAUOUT

      REAL*8, save :: TAUGSA(1001,14),SALBTG(768,14),TAUTGS(768)
     &     ,TAUTGD(122)

      REAL*8 FFKG(4,3),RBBK(3)
      REAL*8, PARAMETER, DIMENSION(14) :: GVALUE = (/.0,.25,.45,.50,.55,
     *     .60,.65,.70,.75,.80,.85,.90,.95,1./)
      REAL*8 CWM,CWE,TIJ,RBBI,RBB,BTAU,CUSPWM,CUSPWE,G,TAU,EG,DELTAU,TI
     *     ,WTJ,WTI,GI,WGI,WGJ,F1,F2,F3,F4,F21,F32,F43,F3221,F4332,A,B,C
     *     ,D,XF,FFCUSP,XEXM,CUSPWT,FFLINR,RB2,RB3,TBB,TB2,TB3,XG,XM
     *     ,XP,RBBB,RI,WRJ,WRI,EI,WEI,WEJ,DELALB,X1,X2,X3,X4,XX,BB,DTAU
      INTEGER I,J,K,KTERPL,IT,JT,IG,JG,ITERPL,IGM,JGP,KG,IR,JR,IE,JE,IEM
     *     ,JEP

      DO 100 I=1,122
      TAUTGD(I)=(I-1)*0.1D0
      IF(I > 24) TAUTGD(I)=(I-24)*0.2D0+2.2D0
      IF(I > 48) TAUTGD(I)=(I-48)*0.5D0+7.0D0
      IF(I > 72) TAUTGD(I)=(I-72)+19.0D0
      IF(I > 96) TAUTGD(I)=(I-96)*5.0D0+40.0D0
      IF(I > 112) TAUTGD(I)=(I-112)*100.0D0+100.0D0
      IF(I == 121) TAUTGD(I)=9999.99D0
      IF(I == 122) TAUTGD(I)=12000.0D0
  100 CONTINUE

      DO 110 I=1,768
      IF(I < 602) TAUTGS(I)=(I-1)*0.05D0
      IF(I > 601) TAUTGS(I)=(I-601)*0.50D0+30.0D0
      IF(I > 741) TAUTGS(I)=(I-741)*50.0D0+100.D0
      IF(I > 758) TAUTGS(I)=(I-758)*1000.D0
  110 CONTINUE

      DO 130 J=1,13
      DO 120 I=1,768
      CWM=0.5
      CWE=0.5
      IF(I > 759) CWM=0.0
      IF(I > 759) CWE=0.0
      TIJ=TAUTGS(I)
      CALL SPLINE(TAUTGD,TGDATA(1,J),122,TIJ,RBBI,CWM,CWE,0)
      SALBTG(I,J)=RBBI
  120 CONTINUE
  130 CONTINUE
      DO 150 J=1,13
      DO 140 I=2,1000
      RBB=(I-1)*0.001D0
      CWM=0.5
      CWE=0.5
      CALL SPLINE(SALBTG(1,J),TAUTGS,768,RBB,BTAU,CWM,CWE,0)
      TAUGSA(I,J)=BTAU
  140 CONTINUE
  150 CONTINUE
      SALBTG(1,:) = 0                                ! 1:14
      TAUGSA(1,:) = 0
      TAUGSA(1001,:) = 10000

      SALBTG(:,14) = SALBTG(:,13)*2 - SALBTG(:,12)   ! 1:768
      TAUGSA(:,14) = TAUGSA(:,13)*2 - TAUGSA(:,12)   ! 1:1001
      RETURN

      ENTRY GTSALB(GIN,TAUIN,RBBOUT,RBBIN,EGIN,TAUOUT,KGTAUR)

      KTERPL=0
      CUSPWM=0.5
      CUSPWE=0.5

      G=GIN
      TAU=TAUIN
      RBB=RBBIN
      EG=EGIN

      RBBOUT=0.0
      TAUOUT=0.0
C                                           ---------------------------
C                                           OPTICAL DEPTH INTERPOLATION
C                                           0.05 ON (0.00 < TAU < 30.0)
C                                           0.50 ON (30.0 < TAU < 100.)
C                                           50.0 ON (100. < TAU < 1000)
C                                           ---------------------------

      IF(KGTAUR == 2) GO TO 300

  200 CONTINUE

      DELTAU=0.05D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT > 599) GO TO 210
      WTJ=TI-IT
      WTI=1.D0-WTJ
      IT=IT+1
      GO TO 240
  210 CONTINUE
      DELTAU=0.50D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT > 199) GO TO 220
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+541
      GO TO 240
  220 CONTINUE
      DELTAU=50.0D0
      TI=TAU/DELTAU
      IT=TI
      IF(IT > 19) GO TO 230
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+649
      GO TO 240
  230 CONTINUE
      DELTAU=1000.0D0
      TI=TAU/DELTAU
      IT=TI
      WTJ=TI-IT
      WTI=1.0-WTJ
      IT=IT+758
  240 CONTINUE
      JT=IT+1

C                                    ---------------------------------
C                                    ASYMMETRY PARAMETER INTERPOLATION
C                                    0.05 CUBIC SPLINE (0.5 < G < 0.9)
C                                    0.25 QUADRATIC ON (0.0 < G < 0.5)
C                                    LINEAR EXTRAP FOR (.95 < G < 1.0)
C                                    ---------------------------------

      GI=G*20.D0
      IF(GI > 10.0) GO TO 250
      IG=2
      JG=3
      ITERPL=1
      GO TO 260

  250 CONTINUE
      ITERPL=4
      IG=GI
      WGJ=GI-IG
      WGI=1.D0-WGJ
      IG=IG-6
      IF(IG > 12) THEN
      ITERPL=2
      IG=12
      ENDIF
      JG=IG+1

  260 CONTINUE

      IGM=IG-1
      JGP=JG+1

      K=0
      DO 270 KG=IGM,JGP
      K=K+1
      F1=SALBTG(IT-1,KG)
      F2=SALBTG(IT  ,KG)
      F3=SALBTG(JT  ,KG)
      F4=SALBTG(JT+1,KG)
      IF(IT == 1) F1=-F3
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WTJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XEXM = (1-2*XF)**2  ! = XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      FFKG(K,1)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      FFKG(K,2)=F2
      FFKG(K,3)=F3
  270 CONTINUE

      IF(ITERPL < 4) GO TO 290

      DO 280 K=1,3
      F1=FFKG(1,K)
      F2=FFKG(2,K)
      F3=FFKG(3,K)
      F4=FFKG(4,K)
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WGJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XEXM = (1-2*XF)**2  ! = XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      RBBK(K)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
  280 CONTINUE
      RBB=RBBK(1)
      RB2=RBBK(2)
      RB3=RBBK(3)
      TBB=TAU
      TB2=TAUTGS(IT)
      TB3=TAUTGS(JT)
      IF(KGTAUR == 1) RETURN
      IF(KTERPL == 1) GO TO 400
      GO TO 300

  290 CONTINUE
      XG=G*2.D0-0.5D0
      IF(ITERPL == 2) XG=G*10.D0-9.D0
      XM=1.D0-XG-XG
      XP=1.D0+XG+XG
      RBB=XM*XP*FFKG(ITERPL+1,1)-XG*XM*FFKG(ITERPL,1)+XG*XP*FFKG(4,1)
      RB2=XM*XP*FFKG(ITERPL+1,2)-XG*XM*FFKG(ITERPL,2)+XG*XP*FFKG(4,2)
      RB3=XM*XP*FFKG(ITERPL+1,3)-XG*XM*FFKG(ITERPL,3)+XG*XP*FFKG(4,3)

      IF(KGTAUR == 1) RETURN
      IF(KTERPL == 1) GO TO 400

  300 CONTINUE
      RBBB=RBB

      RI=RBB*1000.D0
      IR=RI
      WRJ=RI-IR
      WRI=1.D0-WRJ
      IR=IR+1
      JR=IR+1

      EI=EG*20.D0
      IF(EI > 10.0) GO TO 310
      IE=2
      JE=3
      ITERPL=1
      GO TO 320
  310 CONTINUE

      ITERPL=4
      IE=EI
      WEJ=EI-IE
      WEI=1.D0-WEJ
      IE=IE-6
      IF(IE > 12) THEN
      ITERPL=2
      IE=12
      ENDIF
      JE=IE+1
  320 CONTINUE

      DELALB=0.001D0
      IEM=IE-1
      JEP=JE+1
      K=0
      DO 330 KG=IEM,JEP
      K=K+1
      F1=TAUGSA(IR-1,KG)
      F2=TAUGSA(IR  ,KG)
      F3=TAUGSA(JR  ,KG)
      F4=TAUGSA(JR+1,KG)
      IF(IR == 1) F1=-F3
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WRJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XEXM = (1-2*XF)**2  ! = XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      FFKG(K,1)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      FFKG(K,2)=F2
      FFKG(K,3)=F3
  330 CONTINUE
      X1=GVALUE(IE-1)
      X2=GVALUE(IE  )
      X3=GVALUE(JE  )
      X4=GVALUE(JE+1)
      XX=WEJ
      IF(ITERPL < 4) GO TO 350

      DO 340 K=1,3
      F1=FFKG(1,K)
      F2=FFKG(2,K)
      F3=FFKG(3,K)
      F4=FFKG(4,K)
      F21=(F2-F1)
      F32=(F3-F2)
      F43=(F4-F3)
      F3221=(F32+F21)*0.5D0
      F4332=(F43+F32)*0.5D0
      A=F2
      B=F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)
      XF=WEJ
      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XEXM = (1-2*XF)**2  ! = XE**2
      CUSPWT=(1.0-XEXM)*CUSPWM+XEXM*CUSPWE
      FFLINR=A+XF*F32
      RBBK(K)=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
  340 CONTINUE
      TBB=RBBK(1)
      TB2=RBBK(2)
      TB3=RBBK(3)

      IF(KTERPL == 1) GO TO 400
      GO TO 390

  350 CONTINUE
      XG=EG*2.D0-0.5D0
      IF(ITERPL == 2) XG=G*10.D0-9.D0
      XM=1.D0-XG-XG
      XP=1.D0+XG+XG
      TBB=XM*XP*FFKG(ITERPL+1,1)-XG*XM*FFKG(ITERPL,1)+XG*XP*FFKG(4,1)
      TB2=XM*XP*FFKG(ITERPL+1,2)-XG*XM*FFKG(ITERPL,2)+XG*XP*FFKG(4,2)
      TB3=XM*XP*FFKG(ITERPL+1,3)-XG*XM*FFKG(ITERPL,3)+XG*XP*FFKG(4,3)
      IF(KTERPL == 1) GO TO 400
  390 CONTINUE
      KTERPL=1
      TAU=TBB
      G=EGIN
      GO TO 200
  400 CONTINUE
      IF(ABS(WTI*WTJ) < 0.1D0) DTAU=(RBBB-RB2)/(RB3-RB2)
      IF(ABS(WTI*WTJ) >= 0.1D0) THEN
      C=(RB3-RBB)/WTI-(RBB-RB2)/WTJ
      B=(RBB-RB2)/WTJ-WTJ*C
      A=RB2
      BB=B*B+4.D0*C*(RBBB-A)
      IF(BB > 0.D0) DTAU=(SQRT(BB)-B)/(C+C)
      ENDIF
      TAUOUT=(IT-1+DTAU)*DELTAU
      RBBOUT=RBBB

      RETURN
CCC   END SUBROUTINE GTSALB
      END SUBROUTINE SETGTS


      SUBROUTINE SGPGXG(XMU,TAU,G,GG)
      use SURF_ALBEDO, only : gtau
      IMPLICIT NONE
C     ----------------------------------------------------------------
C     COSBAR ADJUSTMENT TO REPRODUCE THE SOLAR ZENITH ANGLE DEPENDENCE
C     FOR AEROSOL ALBEDO FOR OPTICAL THICKNESSES [0.0 < TAU < 10000.0]
C     ----------------------------------------------------------------
      REAL*8, INTENT(IN) :: XMU,TAU,G
      REAL*8, INTENT(OUT) :: GG
      REAL*8 XI,WXI,WXJ,GI,WGI,WGJ,TI,WTJ,WTI
      INTEGER IX,JX,IG,JG,IT,IT0,JT
C                          -------------------------------------------
C                          XMU (COSZ) SOLAR ZENITH ANGLE INTERPOLATION
C                          DATA INTERVAL:  0.02  ON  [0.0 < XMU < 1.0]
C                          -------------------------------------------

      XI=XMU*50.D0+0.999999D0  ! >1 since XMU=COSZ>.001
      IX=XI
      JX=IX+1
      WXJ=XI-IX
      WXI=1.D0-WXJ

C                                      -------------------------------
C                                      COSBAR DEPENDENCE INTERPOLATION
C                                         0.10 ON [0.0 < COSBAR < 1.0]
C                                      -------------------------------

      GI=G*10.D0
      IG=GI
      WGJ=GI-IG
      WGI=1.D0-WGJ
      IG=IG+1
      JG=IG+1

C                            -----------------------------------------
C                               AEROSOL TAU INTERPOLATION INTERVALS
C                            -----------------------------------------
C                            dTau      1      1  (Lin Int)  61     62
C                            0.10 ON [0.00 , 0.00 < TAU < 6.00 , 6.10]
C                                     63     64            92     93
C                            0.50 ON [5.50 , 6.00 < TAU < 20.0 , 20.5]
C                                     94     95           111    112
C                            5.00 ON [15.0 , 20.0 < TAU < 100. , 105.]
C                                     113    114          132    133
C                            50.0 ON [50.0 , 100. < TAU < 1000 , 1050]
C                                            134          143
C                            1000 ON [     , 1000 < TAU < 10000,     ]
C                            -----------------------------------------

      IF(TAU < 6.D0) THEN
      TI=TAU*10.D0+1.
      IT=TI
      WTJ=TI-IT
      IT0=0

      ELSEIF(TAU < 20.D0) THEN
      TI=(TAU-6.D0)*2.00D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=62

      ELSEIF(TAU < 100.D0) THEN
      TI=(TAU-20.D0)*0.20D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=93

      ELSEIF(TAU < 1000.D0) THEN
      TI=(TAU-100.D0)*0.02D0+2.D0
      IT=TI
      WTJ=TI-IT
      IT0=112

      ELSE
      TI=TAU*0.001D0+1.D-6
      IT=TI
      WTJ=TI-IT
      IF(IT > 9) IT=9
      IT0=133
      ENDIF

      WTI=1.D0-WTJ
      IT=IT+IT0
      JT=IT+1
      GG=WGI*(WTI*(WXI*GTAU(IX,IG,IT)+WXJ*GTAU(JX,IG,IT))
     +      + WTJ*(WXI*GTAU(IX,IG,JT)+WXJ*GTAU(JX,IG,JT)))
     +  +WGJ*(WTI*(WXI*GTAU(IX,JG,IT)+WXJ*GTAU(JX,JG,IT))
     +      + WTJ*(WXI*GTAU(IX,JG,JT)+WXJ*GTAU(JX,JG,JT)))

      RETURN

      RETURN
      END SUBROUTINE SGPGXG

      SUBROUTINE SPLINE(X,F,NXF,XX,FF,CUSPWM,CUSPWE,KXTRAP)
      IMPLICIT NONE

      INTEGER, intent(in)  :: NXF,              KXTRAP
      REAL*8 , intent(in)  :: X(NXF),F(NXF),XX, CUSPWM,CUSPWE
      REAL*8 , intent(out) :: FF

C---------------------------------------------------------------------
C
C    SPLINE locates XX between points (F2,X2)(F3,X3) on 4-point spread
C       and returns 4-point Cubic Spline interpolated value FF = F(XX)
C
C    Quadratic Derivatives of Spline are continuous at (F2,X2),(F3,X3)
C    (X-Coordinate may be specified in increasing or decreasing order)
C
C---------------------------------------------------------------------
C
C    Custom Control Parameters:  CUSPWM,CUSPWE,KXTRAP
C------------------------------
C
C    In cases where data points are unevenly spaced and/or data points
C    exhibit abrupt changes in value, Spline Interpolation may produce
C    undesirable bulging of interpolated values. In more extreme cases
C    Linear Interpolation may be less problematic to use.
C
C    Interpolation can be weighted between: Cubic Spline and Linear by
C    adjusting weights CUSPWM and CUSPWE to values between 1.0 and 0.0
C
C    CUSPWM = Cubic Spline Weight at the (X2-X3) Interval Mid-point
C    CUSPWE = Cubic Spline Weight at the (X2-X3) Interval End-points
C
C    For example, with:
C
C    CUSPWM=1.0,CUSPWE=1.0  FF returns Cubic Spline interpolated value
C    CUSPWM=0.0,CUSPWE=0.0  FF returns   Linearly   interpolated value
C
C---------------------------------------------------------------------
C
C     Extrapolation for XX outside of defined interval:  X(1)<->X(NXF)
C
C               KXTRAP = 0    No Extrapolation  (i.e., sets F(XX)=0.0)
C                        1    Fixed Extrapolation (F(XX) = edge value)
C                        2    Linear Extrapolation using 2 edge points
C
C---------------------------------------------------------------------

      REAL*8 x1,x2,x3,x4, x21,x32,x43,x31,x42, betw,FFCUSP,FFLINR,CUSPWT
      REAL*8 f1,f2,f3,f4, f21,f32,f43,f3221,f4332,  a,b,c,d, xf,xe,xexm
      INTEGER K

      K=2
      X2=X(K)
      X3=X(NXF-1)
      BETW=(XX-X2)*(X3-XX)
      IF(BETW <= 0.D0) GO TO 120

  100 CONTINUE
      K=K+1
      X3=X(K)
      BETW=(XX-X2)*(X3-XX)
      IF(BETW >= 0.D0) GO TO 110
      X2=X3
      GO TO 100

  110 CONTINUE
      F3=F(K)
      F4=F(K+1)
      X4=X(K+1)
      F2=F(K-1)
      X2=X(K-1)
      F1=F(K-2)
      X1=X(K-2)
      X21=X2-X1
      X31=X3-X1
      X32=X3-X2
      X43=X4-X3
      X42=X4-X2
      F21=(F2-F1)/(X21*X21)
      F32=(F3-F2)/(X32*X32)
      F43=(F4-F3)/(X43*X43)
      F3221=(F32+F21)/X31*X21
      F4332=(F43+F32)/X42*X43
      A=F2
      B=X32*F3221
      C=3.D0*F32-F3221-F3221-F4332
      D=(F3221+F4332-F32-F32)/X32
      XF=XX-X2

C                             FFCUSP= Cubic Spline Interpolation Result
C                             -----------------------------------------

      FFCUSP=A+XF*(B+XF*(C+XF*D))
      XE=(X3+X2-XX-XX)/X32
      IF(XE < 0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE

C                                   FFLINR= Linear Interpolation Result
C                                   -----------------------------------
      FFLINR=A+XF*F32*X32
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

C                Edge Point Interval Interpolation and/or Extrapolation
C                ------------------------------------------------------
  120 CONTINUE
      BETW=(X2-XX)*(X3-X2)
      IF(BETW < 0.D0) GO TO 140

C                          X(1),X(2)  Edge Point Interval Interpolation
C                          --------------------------------------------
      X1=X(1)
      F1=F(1)
      F2=F(2)
      X21=X2-X1
      F21=(F2-F1)/X21
      XF=XX-X1
      BETW=(X2-XX)*XF
      IF(BETW < 0.D0) GO TO 130
      F3=F(3)
      X3=X(3)
      X32=X3-X2
      X31=X3-X1
      C=((F3-F2)/X32-F21)/X31
      B=F21-X21*C
      A=F1
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F21
      XE=1.D0-2.D0*XF/X21
      IF(XE < 0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

  130 CONTINUE
C                  Extrapolation for XX Outside of Interval X(1) - X(2)
C                  ----------------------------------------------------
C                  IF(KXTRAP == 0)  (No Extrapolation:  sets F(XX)=0.0)
C                  IF(KXTRAP == 1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP == 2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP == 0) FF=0.D0
      IF(KXTRAP == 1) FF=F1
      IF(KXTRAP == 2) FF=F1+XF*F21
      GO TO 160

  140 CONTINUE
C                    X(NXF-1),X(NXF)  Edge Point Interval Interpolation
C                    --------------------------------------------------
      F3=F(NXF)
      X3=X(NXF)
      F2=F(NXF-1)
      X2=X(NXF-1)
      X32=X3-X2
      F32=(F3-F2)/X32
      XF=XX-X3
      BETW=(X2-XX)*(XX-X3)
      IF(BETW < 0.D0) GO TO 150
      F1=F(NXF-2)
      X1=X(NXF-2)
      X21=X2-X1
      X31=X3-X1
      F21=(F2-F1)/X21
      XF=XX-X2

C                    3-Point Quadratic Interpolation for Edge Intervals
C                    --------------------------------------------------
C
C      (Edge Option)     ----------------------------------------------
C                        For Linear Interpolation within Edge Intervals
C                        between X(1),X(2), and between X(NXF-1),X(NXF)
C                        set the value of coefficient C below, to C=0.0
C                        ----------------------------------------------

      C=(F32-F21)/X31
      B=F21+X21*C
      A=F2
      FFCUSP=A+XF*(B+XF*C)
      FFLINR=A+XF*F32
      XE=1.D0-2.D0*XF/X32
      IF(XE < 0.D0) XE=-XE
      XEXM=XE**2
      CUSPWT=(1.D0-XEXM)*CUSPWM+XEXM*CUSPWE
      FF=FFCUSP*CUSPWT+FFLINR*(1.D0-CUSPWT)
      GO TO 160

  150 CONTINUE
C              Extrapolation for X Outside of Interval  X(NXF-1)-X(NXF)
C              --------------------------------------------------------
C                  IF(KXTRAP == 0)  (No Extrapolation:  sets F(XX)=0.0)
C                  IF(KXTRAP == 1)  (Extrapolation at Fixed Edge Value)
C                  IF(KXTRAP == 2)  (2 Edge Point Linear Extrapolation)

      IF(KXTRAP == 0) FF=0.D0
      IF(KXTRAP == 1) FF=F3
      IF(KXTRAP == 2) FF=F3+XF*(F3-F2)/(X3-X2)

  160 CONTINUE
      RETURN
      END SUBROUTINE SPLINE

      MODULE RAD_COM
!@sum  RAD_COM Model radiation arrays and parameters
!@auth Original Development Team
!@ver  1.0
      IMPLICIT NONE
      SAVE

C**** DEFAULT ORBITAL PARAMETERS FOR EARTH
C**** Note PMIP runs had specified values that do not necesarily
C**** coincide with those used as the default, or the output of ORBPAR.
C****                    OMEGT          OBLIQ        ECCEN
C**** DEFAULT (2000 AD): 282.9          23.44        0.0167
C**** PMIP CONTROL:      282.04         23.446       0.016724
C**** PMIP 6kyr BP:      180.87         24.105       0.018682
C**** PMIP LGM (21k):    294.42         22.949       0.018994
!@param OMEGT_def precession angle (degrees from vernal equinox)
      real*8, parameter :: omegt_def = 282.9d0
!@param OBLIQ_def obliquity angle  (degrees)
      real*8, parameter :: obliq_def = 23.44d0
!@param ECCN_def eccentricity
      real*8, parameter :: eccn_def  = .0167d0
!@var OMEGT,OBLIQ,ECCN actual orbital parameters used
      real*8 OMEGT,OBLIQ,ECCN

C**** Database parameters to control orbital parameter calculation
C**** Note: setting calc_orb_par with paleo_orb_yr=-50 (i.e. year 2000)
C**** does not produce exactly the same as the default values.
!@dbparam calc_orb_par_year = PALEO YEAR (BP) to calculate orbital parameters
!@+       on first timestep of each year. PALEO YEAR incremented yearly from
!@+       the value of JYEAR so set IYEAR appropriately (1950 usually).
      integer :: calc_orb_par_year = 0
!@dbparam calc_orb_par = 1 to calc orbital parameters
      integer :: calc_orb_par = 0
!@dbparam paleo_orb_yr is paleo year (BP) for orbital calc
      real*8 :: paleo_orb_yr = -50.  ! (i.e. 2000AD)
!@dbparam calc_orb_par_sp = 1 to directly specify orbital parameters
      integer :: calc_orb_par_sp = 0
!@dbparam paleo_orb_par :: directly specifies orbital parameters
      real*8, dimension(3) :: paleo_orb_par = (/ eccn_def, obliq_def,
     *     omegt_def /)

!@var COSZ1 Mean Solar Zenith angle for curr. physics(not rad) time step
      !REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ1
!@var COSZ_day Mean Solar Zenith angle for current day
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: COSZ_day,DUSK

!@dbparam S0X solar constant multiplication factor
      real*8 :: S0 = 1367.
      REAL*8 :: S0X = 1.
!@dbparam S0_yr,S0_day obs.date of solar constant (if 0: time var)
      INTEGER :: S0_yr = 1951 , S0_day = 182

!@var RSDIST,SIND,COSD orbit related variables computed once a day
      REAL*8 :: RSDIST,SIND,COSD

      END MODULE RAD_COM

      SUBROUTINE ALLOC_RAD_COM
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, ONLY : getDomainBounds
      USE RAD_COM, ONLY : COSZ_day, DUSK !,COSZ1
      IMPLICIT NONE
      INTEGER :: I_0H, I_1H, J_0H, J_1H
      call getDomainBounds(grid,
     &               I_STRT_HALO=I_0H, I_STOP_HALO=I_1H,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      ALLOCATE(
!     &     COSZ1(I_0H:I_1H, J_0H:J_1H),
     &     COSZ_day(I_0H:I_1H,J_0H:J_1H),
     &     DUSK(I_0H:I_1H,J_0H:J_1H)
     &     )
      RETURN
      END SUBROUTINE ALLOC_RAD_COM


      SUBROUTINE ORBIT (DOBLIQ,ECCEN,DOMEGVP,VEDAY,EDPY, DAY,
     *                  SDIST,SIND,COSD,SUNLON,SUNLAT,EQTIME)
C****
C**** ORBIT receives orbital parameters and time of year, and returns
C**** distance from Sun, declination angle, and Sun's overhead position.
C**** Reference for following caculations is:  V.M.Blanco and
C**** S.W.McCuskey, 1961, "Basic Physics of the Solar System", pages
C**** 135 - 151.  Existence of Moon and heavenly bodies other than
C**** Earth and Sun are ignored.  Earth is assumed to be spherical.
C****
C**** Program author: Gary L. Russell 2004/11/16
C**** Angles, longitude and latitude are measured in radians.
C****
C**** Input: ECCEN  = eccentricity of the orbital ellipse
C****        OBLIQ  = latitude of Tropic of Cancer
C****        OMEGVP = longitude of perihelion (sometimes Pi is added) =
C****               = spatial angle from vernal equinox to perihelion
C****                 with Sun as angle vertex
C****        DAY    = days measured since 2000 January 1, hour 0
C****
C****        EDPY  = Earth days per year
C****                tropical year = 365.2425 (Gregorgian Calendar)
C****                tropical year = 365      (Generic Year)
C****        VEDAY = Vernal equinox
C****                79.0 (Generic year Mar 21 hour 0)
C****                79.5 (Generic year Mar 21 hour 12 - PMIP standard)
C****                79.3125d0 for days from 2000 January 1, hour 0 till vernal
C****                     equinox of year 2000 = 31 + 29 + 19 + 7.5/24
C****
C**** Intermediate quantities:
C****    BSEMI = semi minor axis in units of semi major axis
C****   PERIHE = perihelion in days since 2000 January 1, hour 0
C****            in its annual revolution about Sun
C****       TA = true anomaly = spatial angle from perihelion to
C****            current location with Sun as angle vertex
C****       EA = eccentric anomaly = spatial angle measured along
C****            eccentric circle (that circumscribes Earth's orbit)
C****            from perihelion to point above (or below) Earth's
C****            absisca (where absisca is directed from center of
C****            eccentric circle to perihelion)
C****       MA = mean anomaly = temporal angle from perihelion to
C****            current time in units of 2*Pi per tropical year
C****   TAofVE = TA(VE) = true anomaly of vernal equinox = - OMEGVP
C****   EAofVE = EA(VE) = eccentric anomaly of vernal equinox
C****   MAofVE = MA(VE) = mean anomaly of vernal equinox
C****   SLNORO = longitude of Sun in Earth's nonrotating reference frame
C****   VEQLON = longitude of Greenwich Meridion in Earth's nonrotating
C****            reference frame at vernal equinox
C****   ROTATE = change in longitude in Earth's nonrotating reference
C****            frame from point's location on vernal equinox to its
C****            current location where point is fixed on rotating Earth
C****   SLMEAN = longitude of fictitious mean Sun in Earth's rotating
C****            reference frame (normal longitude and latitude)
C****
C**** Output: SIND = sine of declination angle = sin(SUNLAT)
C****         COSD = cosine of the declination angle = cos(SUNLAT)
C****       SUNDIS = distance to Sun in units of semi major axis
C****       SUNLON = longitude of point on Earth directly beneath Sun
C****       SUNLAT = latitude of point on Earth directly beneath Sun
C****       EQTIME = Equation of Time =
C****              = longitude of fictitious mean Sun minus SUNLON
C****
C**** From the above reference:
C**** (4-54): [1 - ECCEN*cos(EA)]*[1 + ECCEN*cos(TA)] = (1 - ECCEN^2)
C**** (4-55): tan(TA/2) = sqrt[(1+ECCEN)/(1-ECCEN)]*tan(EA/2)
C**** Yield:  tan(EA) = sin(TA)*sqrt(1-ECCEN^2) / [cos(TA) + ECCEN]
C****    or:  tan(TA) = sin(EA)*sqrt(1-ECCEN^2) / [cos(EA) - ECCEN]
C****
      USE CONSTANT, only : twopi,pi,radian
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: DOBLIQ,ECCEN,DOMEGVP,DAY,VEDAY,EDPY
      REAL*8, INTENT(OUT) :: SIND,COSD,SDIST,SUNLON,SUNLAT,EQTIME

      REAL*8 MA,OMEGVP,OBLIQ,EA,DEA,BSEMI
     *     ,TAofVE,EAofVE,MAofVE,SUNDIS,TA,SUNX,SUNY,SLNORO
     *     ,VEQLON,ROTATE,SLMEAN
c      REAL*8, PARAMETER :: EDAYzY=365.2425d0, VE2000=79.3125d0
c      REAL*8, PARAMETER :: EDAYzY=365d0, VE2000=79d0  ! original parameters
      REAL*8  EDAYzY,VE2000
C****
      VE2000=VEDAY
      EDAYzY=EDPY
      OMEGVP=DOMEGVP*radian
      OBLIQ=DOBLIQ*radian
C**** Determine EAofVE from geometry: tan(EA) = b*sin(TA) / [e+cos(TA)]
C**** Determine MAofVE from Kepler's equation: MA = EA - e*sin(EA)
C**** Determine MA knowing time from vernal equinox to current day
C****
      BSEMI  = SQRT (1 - ECCEN*ECCEN)
      TAofVE = - OMEGVP
      EAofVE = ATAN2 (BSEMI*SIN(TAofVE), ECCEN+COS(TAofVE))
      MAofVE = EAofVE - ECCEN*SIN(EAofVE)
C     PERIHE = VE2000 - MAofVE*EDAYzY/TWOPI
      MA     = MODULO (TWOPI*(DAY-VE2000)/EDAYzY + MAofVE, TWOPI)
C****
C**** Numerically invert Kepler's equation: MA = EA - e*sin(EA)
C****
      EA  = MA + ECCEN*(SIN(MA) + ECCEN*SIN(2*MA)/2)
   10 dEA = (MA - EA + ECCEN*SIN(EA)) / (1 - ECCEN*COS(EA))
      EA  = EA + dEA
      IF(ABS(dEA).gt.1d-10)  GO TO 10
C****
C**** Calculate distance to Sun and true anomaly
C****
      SUNDIS = 1 - ECCEN*COS(EA)
      TA     = ATAN2 (BSEMI*SIN(EA), COS(EA)-ECCEN)
      SDIST  = SUNDIS*SUNDIS   ! added for compatiblity
C****
C**** Change reference frame to be nonrotating reference frame, angles
C**** fixed according to stars, with Earth at center and positive x
C**** axis be ray from Earth to Sun were Earth at vernal equinox, and
C**** x-y plane be Earth's equatorial plane.  Distance from current Sun
C**** to this x axis is SUNDIS sin(TA-TAofVE).  At vernal equinox, Sun
C**** is located at (SUNDIS,0,0).  At other times, Sun is located at:
C****
C**** SUN = (SUNDIS cos(TA-TAofVE),
C****        SUNDIS sin(TA-TAofVE) cos(OBLIQ),
C****        SUNDIS sin(TA-TAofVE) sin(OBLIQ))
C****
      SIND   = SIN(TA-TAofVE) * SIN(OBLIQ)
      COSD   = SQRT (1 - SIND*SIND)
      SUNX   = COS(TA-TAofVE)
      SUNY   = SIN(TA-TAofVE) * COS(OBLIQ)
      SLNORO = ATAN2 (SUNY,SUNX)
C****
C**** Determine Sun location in Earth's rotating reference frame
C**** (normal longitude and latitude)
C****
      VEQLON = TWOPI*VE2000 - PI + MAofVE - TAofVE  !  modulo 2*Pi
      ROTATE = TWOPI*(DAY-VE2000)*(EDAYzY+1)/EDAYzY
      SUNLON = MODULO (SLNORO-ROTATE-VEQLON, TWOPI)
      IF(SUNLON.gt.PI)  SUNLON = SUNLON - TWOPI
      SUNLAT = ASIN (SIN(TA-TAofVE)*SIN(OBLIQ))
C****
C**** Determine longitude of fictitious mean Sun
C**** Calculate Equation of Time
C****
      SLMEAN = PI - TWOPI*(DAY-FLOOR(DAY))
      EQTIME = MODULO (SLMEAN-SUNLON, TWOPI)
      IF(EQTIME.gt.PI)  EQTIME = EQTIME - TWOPI
C****
      RETURN
      END SUBROUTINE ORBIT

      SUBROUTINE CALC_ZENITH_ANGLE
!@sum calculate zenith angle for current time step
!@auth Gavin Schmidt (from RADIA)
      USE CONSTANT, only : twopi
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE MODEL_COM, only : itime,nday,dtsrc
      !USE RAD_COM, only : cosz1
      USE RAD_COSZ0, only : coszt
      USE FLUXES, only : atmocn
      IMPLICIT NONE
      INTEGER JTIME
      REAL*8 ROT1,ROT2

      JTIME=MOD(ITIME,NDAY)
      ROT1=(TWOPI*JTIME)/NDAY
      ROT2=ROT1+TWOPI*DTsrc/SECONDS_PER_DAY
      CALL COSZT (ROT1,ROT2,atmocn%COSZ1)

      END SUBROUTINE CALC_ZENITH_ANGLE

      SUBROUTINE init_RAD
      USE Dictionary_mod
      USE MODEL_COM, only : dtsrc,iyear1,modelEclock
      USE DOMAIN_DECOMP_1D, only : write_parallel, am_i_root
      USE RAD_COM, only : s0x,s0_yr,s0_day
     *     ,obliq,eccn,omegt,obliq_def,eccn_def,omegt_def
     *     ,calc_orb_par,paleo_orb_yr
     *     ,calc_orb_par_sp,paleo_orb_par,calc_orb_par_year
      use RAD_COSZ0, only : cosz_init
      use SURF_ALBEDO, only : gtau,tgdata
      USE FILEMANAGER, only : openunit,closeunit
      IMPLICIT NONE
      integer :: iu
      REAL*8 pyear
      integer :: jyear

      character(len=300) :: out_line
      character*6 :: skip

C**** sync radiation parameters from input
      call sync_param( "S0X", S0X )
      call sync_param( "S0_yr", S0_yr )
      call sync_param( "S0_day", S0_day )

C**** Set orbital parameters appropriately

      jyear = modelEclock%getYear()

      if (calc_orb_par_year.ne.0) then ! calculate from paleo-year
        ! 0 BP is defined as 1950CE
        pyear = 1950.+JYEAR-IYEAR1-calc_orb_par_year
        write(out_line,*)
     *  "   Calculating Orbital Params for year : ",
     *  pyear,"     (CE);"
        call orbpar(pyear,eccn, obliq, omegt)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*) 'calc_orb_par_year =',calc_orb_par_year
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
     *  "   Calculating Orbital Params for year : ",
     *  pyear,"     (CE);"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*) " Orbital Parameters Calculated:"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f8.7,a,f8.7,a)') "   Eccentricity: ",eccn,
     *       " (default = ",eccn_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f9.6,a,f9.6,a)') "   Obliquity (degs): ",
     *       obliq,
     *       " (default = ",obliq_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f7.3,a,f7.3,a)')
     *       "   Precession (degs from ve): ",
     *       omegt," (default = ",omegt_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
      elseif (calc_orb_par.eq.1) then ! calculate from paleo-year
        pyear=1950.-paleo_orb_yr ! since 0 BP is defined as 1950CE
        call orbpar(pyear, eccn, obliq, omegt)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*) " Orbital Parameters Calculated:"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f8.0,a,f8.0,a)') "   Paleo-year: ",pyear,"
     *       (CE);", paleo_orb_yr," (BP)"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f8.7,a,f8.7,a)') "   Eccentricity: ",eccn,
     *       " (default = ",eccn_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f9.6,a,f9.6,a)') "   Obliquity (degs): ",
     *       obliq,
     *       " (default = ",obliq_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f7.3,a,f7.3,a)')
     *       "   Precession (degs from ve): ",
     *       omegt," (default = ",omegt_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
      elseif(calc_orb_par_sp.eq.1) then
        omegt=paleo_orb_par(3)
        obliq=paleo_orb_par(2)
        eccn=paleo_orb_par(1)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*) " Orbital Parameters Specified:"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f8.0,a,f8.0,a)') "   Paleo-year: ",pyear,"
     *       (CE);", paleo_orb_yr," (BP)"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f8.7,a,f8.7,a)') "   Eccentricity: ",eccn,
     *       " (default = ",eccn_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f9.6,a,f9.6,a)') "   Obliquity (degs): ",
     *       obliq,
     *       " (default = ",obliq_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f7.3,a,f7.3,a)')
     *       "   Precession (degs from ve): ",
     *       omegt," (default = ",omegt_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
      else  ! set from defaults (defined in CONSTANT module)
        omegt=omegt_def
        obliq=obliq_def
        eccn=eccn_def
      end if

      call cosz_init

c
c FOR SNOW ALBEDO ON SEA ICE:
c
C-----------------------------------------------------------------------
CR(1) Reads GTAU Asymmetry Parameter Conversion Table used within SGPGXG
C
C       (SGPGXG does Multiple Scattering Parameterization used in SOLAR)
C       ----------------------------------------------------------------

      call openunit('GTAU',iu,.true.,.true.)
      READ (IU) GTAU,TGDATA
      CALL SETGTS
      call closeunit(iu)

      return
      END SUBROUTINE init_RAD

      subroutine alloc_drv_atm
      use domain_decomp_atm
      use seaice_com, only : alloc_icestate_type,si_atm
      implicit none
      call alloc_ocean
      call alloc_geom
      call alloc_adiag
      call alloc_fluxes
      call alloc_rad_com
      call alloc_core_data
      call alloc_icestate_type(grid,si_atm,'OCEAN')
      return
      end subroutine alloc_drv_atm

      SUBROUTINE INPUT_atm (istart,istart_fixup,do_IC_fixups,
     &     is_coldstart,KDISK_restart,IRANDI)
      implicit none
!@var istart start(1-8)/restart(>8)  option
      integer :: istart,istart_fixup,do_IC_fixups
      LOGICAL :: is_coldstart
      INTEGER :: KDISK_restart
      INTEGER :: IRANDI
c
      call init_RAD
      call read_core_data
      call daily_atm(.false.)
      return
      end subroutine input_atm

      SUBROUTINE DAILY_atm(end_of_day)
!@sum  DAILY performs daily tasks at end-of-day and maybe at (re)starts
!@auth Original Development Team
!@ver  1.0
      USE MODEL_COM, only : iyear1,modelEclock
      USE RAD_COM, only : RSDIST,COSD,SIND,COSZ_day,DUSK,
     *     omegt,obliq,eccn,omegt_def,obliq_def,eccn_def,
     *     calc_orb_par_year
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT, WRITE_PARALLEL
      use RAD_COSZ0, only : daily_cosz
      IMPLICIT NONE
      REAL*8 :: SUNLON,SUNLAT,LAM,EDPY,VEDAY,PYEAR
      INTEGER i,j,l,iy,it
      character(len=300) :: out_line
      LOGICAL, INTENT(IN) :: end_of_day
      integer :: jday,jyear

      jyear = modelEclock%getYear()
      jday  = modelEclock%getDayOfYear()

C**** Tasks to be done at end of day and at each start or restart

C**** CALCULATE SOLAR ANGLES AND ORBIT POSITION
C**** This is for noon (GMT) for new day.

C**** The orbital calculation will need to vary depending on the kind
C**** of calendar adopted (i.e. a generic 365 day year, or a transient
C**** calendar including leap years etc.).  For transient calendars the
C**** JDAY passed to orbit needs to be adjusted to represent the number
C**** of days from Jan 1 2000AD.
c      EDPY=365.2425d0, VEDAY=79.3125d0  ! YR 2000AD
c      JDAY => JDAY + 365 * (JYEAR-2000) + appropriate number of leaps
C**** Default calculation (no leap, VE=Mar 21 hr 0)
c      EDPY=365d0 ; VEDAY=79d0           ! Generic year
C**** PMIP calculation (no leap, VE=Mar 21 hr 12)
      EDPY=365d0 ; VEDAY=79.5d0           ! Generic year
C**** Set orbital parameters appropriately
      if (calc_orb_par_year.ne.0.and.JDAY.eq.1) then ! calculate from paleo-year
        pyear = 1950.+JYEAR-IYEAR1-calc_orb_par_year
        ! 0 BP is defined as 1950CE
        call orbpar(pyear, eccn, obliq, omegt)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*) " Orbital Parameters Calculated:"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
     *  "   Calculating Orbital Params for year : ",
     *  pyear,"     (CE);", ' JYEAR, IYEAR1,calc_orb_par_year=',
     *  JYEAR,IYEAR1,calc_orb_par_year
        call write_parallel(trim(out_line),unit=6)
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f8.7,a,f8.7,a)') "   Eccentricity: ",eccn,
     *       " (default = ",eccn_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f9.6,a,f9.6,a)') "   Obliquity (degs): ",
     *       obliq,
     *       " (default = ",obliq_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,'(a,f7.3,a,f7.3,a)')
     *       "   Precession (degs from ve): ",
     *       omegt," (default = ",omegt_def,")"
        call write_parallel(trim(out_line),unit=6)
        write(out_line,*)
        call write_parallel(trim(out_line),unit=6)
      end if
      CALL ORBIT (OBLIQ,ECCN,OMEGT,VEDAY,EDPY,REAL(JDAY,KIND=8)-.5
     *     ,RSDIST,SIND,COSD,SUNLON,SUNLAT,LAM)
      call daily_cosz(sind,cosd,cosz_day,dusk)

      RETURN
      END SUBROUTINE DAILY_atm

      subroutine atm_phase1
      use model_com, only : idacc
      use mdiag_com, only : ia_cpl
      use seaice_com, only : si_ocn,iceocn ! temporary until melt_si calls
      use fluxes, only : atmocn,atmice     ! are moved to ocean driver
      implicit none

C**** calculate solar zenith angle for current time step
      CALL CALC_ZENITH_ANGLE

      IDACC(ia_cpl)=IDACC(ia_cpl)+1

C***  Get the forcing data
      call get_ocean_forcings

C**** FIRST CALL MELT_SI SO THAT TOO SMALL ICE FRACTIONS ARE REMOVED
C**** AND ICE FRACTION CAN THEN STAY CONSTANT UNTIL END OF TIMESTEP
      CALL MELT_SI(si_ocn,iceocn,atmocn,atmice)

      return
      end subroutine atm_phase1

      subroutine atm_phase2
      end subroutine atm_phase2

      subroutine def_rsf_acc(fid,r4_on_disk)
      use model_com, only : idacc
      use mdiag_com, only : monacc
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use diag_com, only : grid_ioptr,aij_ioptr
      implicit none
      integer fid   !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4
      call defvar(grid,fid,idacc,'monacc(twelve)')
      call defvar(grid,fid,idacc,'idacc(nsampl)')
      call defvar(grid_ioptr,fid,aij_ioptr,'aij(im,dist_jm,kaij)',
     &     r4_on_disk=r4_on_disk)
      return
      end subroutine def_rsf_acc

      subroutine new_io_acc(fid,iaction)
      use model_com, only : ioread,iowrite,iowrite_single,idacc
      use mdiag_com, only : monacc
      use domain_decomp_atm, only : grid
      use pario, only : write_data,read_data,
     &     write_dist_data,read_dist_data
      use diag_com, only : grid_ioptr,aij_ioptr
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite,iowrite_single) ! output to restart or acc file
        call write_data(grid,fid,'monacc',monacc)
        call write_data(grid,fid,'idacc',idacc)
        call write_dist_data(grid_ioptr,fid,'aij',aij_ioptr)
      case (ioread)            ! input from restart or acc file
        call read_data(grid,fid,'monacc',monacc,bcast_all=.true.)
        call read_data(grid,fid,'idacc',idacc,bcast_all=.true.)
        call read_dist_data(grid_ioptr,fid,'aij',aij_ioptr)
      end select
      return
      end subroutine new_io_acc

      subroutine def_rsf_atmvars(fid)
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use fluxes, only : atmocn,atmice
      implicit none
      integer :: fid
      character(len=17) :: ijstr
      ijstr='(dist_im,dist_jm)'
      call defvar(grid,fid,atmocn%cosz1,'cosz'//ijstr)
      call defvar(grid,fid,atmocn%cosz1,'cosz_day'//ijstr)
      call defvar(grid,fid,atmocn%cosz1,'lon2d'//ijstr)
      call defvar(grid,fid,atmocn%cosz1,'lat2d'//ijstr)
      call defvar(grid,fid,atmocn%cosz1,'focean'//ijstr)
      call defvar(grid,fid,atmocn%gtemp,'asst'//ijstr)
      call defvar(grid,fid,atmocn%gtempr,'atempr'//ijstr)
      call defvar(grid,fid,atmocn%sss,'sss'//ijstr)
      call defvar(grid,fid,atmocn%ogeoza,'ogeoza'//ijstr)
      call defvar(grid,fid,atmocn%uosurf,'uosurf'//ijstr)
      call defvar(grid,fid,atmocn%vosurf,'vosurf'//ijstr)
      call defvar(grid,fid,atmocn%mlhc,    'mlhc'//ijstr)
      call defvar(grid,fid,atmocn%evapor,'oevapor'//ijstr)
      call defvar(grid,fid,atmice%evapor,'ievapor'//ijstr)
#ifndef STANDALONE_HYCOM
      call def_rsf_icedyn (fid)   ! move this!!!!
#endif
      return
      end subroutine def_rsf_atmvars

      subroutine new_io_atmvars(fid,iorw)
      use model_com, only : iowrite,ioread
      use domain_decomp_atm, only : grid
      use pario, only : write_dist_data,read_dist_data
      use fluxes, only : atmocn,atmice
      use rad_com, only : cosz_day
      use geom, only : lon2d,lat2d
      implicit none
      integer, intent(in) :: fid,iorw
      select case (iorw)
      case (iowrite)            ! output to restart file
        CALL CALC_ZENITH_ANGLE
        call write_dist_data(grid, fid, 'cosz', atmocn%cosz1)
        call write_dist_data(grid, fid, 'cosz_day', cosz_day)
        call write_dist_data(grid, fid, 'lon2d', lon2d)
        call write_dist_data(grid, fid, 'lat2d', lat2d)
        call write_dist_data(grid, fid, 'focean', atmocn%focean)
        call write_dist_data(grid, fid, 'asst',atmocn%gtemp)
        call write_dist_data(grid, fid, 'atempr',atmocn%gtempr)
        call write_dist_data(grid, fid, 'sss',atmocn%sss)
        call write_dist_data(grid, fid, 'ogeoza',atmocn%ogeoza)
        call write_dist_data(grid, fid, 'uosurf',atmocn%uosurf)
        call write_dist_data(grid, fid, 'vosurf',atmocn%vosurf)
        call write_dist_data(grid, fid, 'mlhc',atmocn%mlhc)
        call write_dist_data(grid, fid, 'oevapor',atmocn%evapor)
        call write_dist_data(grid, fid, 'ievapor',atmice%evapor)
      case (ioread)             ! input from restart file
        call read_dist_data(grid, fid, 'asst',atmocn%gtemp)
        call read_dist_data(grid, fid, 'atempr',atmocn%gtempr)
        call read_dist_data(grid, fid, 'sss',atmocn%sss)
        call read_dist_data(grid, fid, 'ogeoza',atmocn%ogeoza)
        call read_dist_data(grid, fid, 'uosurf',atmocn%uosurf)
        call read_dist_data(grid, fid, 'vosurf',atmocn%vosurf)
        call read_dist_data(grid, fid, 'mlhc',atmocn%mlhc)
        call read_dist_data(grid, fid, 'oevapor',atmocn%evapor)
        call read_dist_data(grid, fid, 'ievapor',atmice%evapor)
      end select
#ifndef STANDALONE_HYCOM
      call new_io_icedyn (fid,iorw)  ! move this!!!
#endif
      return
      end subroutine new_io_atmvars

      subroutine daily_diag(newmonth)
      implicit none
      logical, intent(in) :: newmonth
      return
      end subroutine daily_diag

      subroutine reset_adiag(idum)
      use diag_com, only : aij
      implicit none
      integer :: idum ! not used
      call reset_mdiag  ! move this call to higher level!!!
      aij = 0.
      call reset_icdiag ! move this!!!!
      return
      end subroutine reset_adiag

      subroutine def_meta_atmacc(fid)
      use domain_decomp_atm, only : grid
      use diag_com
      use pario, only : defvar,write_attr
      use cdl_mod, only : defvar_cdl
      use geom, only : lon2d_dg,lat2d_dg
      implicit none
      integer :: fid
      real*8 :: r8dum
      call defvar(grid,fid,r8dum,'time')
      call write_attr(grid,fid,'time','reduction','avg')
      call defvar(grid_nodup,fid,lon2d_dg,'lon(im,dist_jm)')
      call defvar(grid_nodup,fid,lat2d_dg,'lat(im,dist_jm)')
      call write_attr(grid,fid,'aij','reduction','sum')
      call write_attr(grid,fid,'aij','split_dim',3)
      call defvar(grid,fid,ia_aij,'ia_aij(kaij)')
      call defvar(grid,fid,scale_aij,'scale_aij(kaij)')
      call defvar(grid,fid,denom_aij,'denom_aij(kaij)')
      call defvar(grid,fid,sname_aij,'sname_aij(sname_strlen,kaij)')
      call defvar_cdl(grid,fid,cdl_aij,'cdl_aij(cdl_strlen,kcdl_aij)')
      return
      end subroutine def_meta_atmacc

      subroutine write_meta_atmacc(fid)
      use domain_decomp_atm, only : grid
      use diag_com
      use pario, only : write_data,write_dist_data
      use cdl_mod, only : write_cdl
      use geom, only : lon2d_dg,lat2d_dg
      use model_com, only : jyear0,jmon0
      implicit none
      integer :: fid
      real*8 :: r8time
      r8time = real(jyear0,kind=8)+(real(jmon0,kind=8)-.5)/12d0
      call write_data(grid,fid,'time',r8time)
      call write_dist_data(grid_nodup,fid,'lon',lon2d_dg)
      call write_dist_data(grid_nodup,fid,'lat',lat2d_dg)
      call write_data(grid,fid,'ia_aij',ia_aij)
      call write_data(grid,fid,'scale_aij',scale_aij)
      call write_data(grid,fid,'denom_aij',denom_aij)
      call write_data(grid,fid,'sname_aij',sname_aij)
      call write_cdl(grid,fid,'cdl_aij',cdl_aij)
      return
      end subroutine write_meta_atmacc

      subroutine calc_derived_acc_atm
      use diag_com
      implicit none
      integer :: k
#ifdef STANDALONE_HYCOM
      do k=1,kaij
        call panam_nodup(aij(:,:,k),aij_nodup(:,:,k))
      enddo
#endif
      return
      end subroutine calc_derived_acc_atm

      subroutine set_ioptrs_atmacc_default
      use diag_com
      use domain_decomp_atm, only : grid
      implicit none
      aij_ioptr => aij
      grid_ioptr => grid
      return
      end subroutine set_ioptrs_atmacc_default

      subroutine set_ioptrs_atmacc_extended
      use diag_com
      implicit none
#ifdef STANDALONE_HYCOM
      aij_ioptr => aij_nodup
#else
      aij_ioptr => aij
#endif
      grid_ioptr => grid_nodup
      return
      end subroutine set_ioptrs_atmacc_extended

      subroutine def_rsf_longacc
      end subroutine def_rsf_longacc

      subroutine new_io_longacc
      end subroutine new_io_longacc

      subroutine print_diags
      end subroutine print_diags

      subroutine finalize_atm
      end subroutine finalize_atm

      SUBROUTINE SURFACE
      USE CONSTANT, only : rgas,lhe,lhs,sha,tf,shv,shi,stbo,deltx,teeny,
     &     grav,rhow
      USE MODEL_COM, only : dtsrc
      use fluxes, only : atmocn,atmice
#ifdef CONSTANT_OCEAN_ALBEDO
      use fluxes, only : ocean_albedo_const
#endif
      use geom, only : axyp,lat2d
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : getDomainBounds,
     &     GLOBALSUM,hasNorthPole
      USE SEAICE, only : xsi,ace1i,rhoi,byrls,solar_ice_frac
     *     ,tfrez,dEidTi,alami,rhos
      USE SEAICE_COM, only : si_ocn
      USE RAD_COM, only : s0,cosz_day
      USE SURF_ALBEDO, ONLY : GET_SURF_ALBEDO
      USE DIAG_COM, only : aij,
     &     ij_roc,ij_rsi,ij_ocnsw,ij_sisw,ij_ocnalb,ij_sialb,ij_msi,
     &     ij_ocheat,ij_sst,ij_sss,ij_foc,ij_taus,ij_tauus,ij_tauvs
      IMPLICIT NONE

      INTEGER I,J,ITYPE,ITER,HEMI
      REAL*8 POICE,POCEAN,PS
     *     ,ELHX,MSI1,MSI2,dF1dTG,HCG1,HCG2,EVHDT,F1DT
     *     ,CD,CM,CH,CQ,EVHEAT,F0,F1,DSHDTG,DQGDTG
     *     ,DEVDTG,DTRDTG,DF0DTG,DFDTG,DTG,dSNdTG
     *     ,dT2,EVAP,F0DT,PTYPE,TG1,TG2,SRHEAT,SNOW
     *     ,SHDT,TRHDT,RHOSRF,RCDMWS,RCDHWS,RCDQWS
     *     ,SHEAT,TRHEAT,T2DEN,T2CON,T2MUL,F2
     *     ,FSRI(2),zoice,zsnow,fmp,zmp,dalbsn,OCNALB,SIALB
      REAL*8 QSAT,DQSATDT,TR4

      real*8 tg,qg,uocean,vocean,ts,tsv,qs,us,vs,ws,taux,tauy

c
c  subroutine ncar_ocean_fluxes (u_del, t, ts, q, qs, z, avail, &
c                              cd, ch, cq, ustar, bstar       )
      real*8, parameter :: vonkarm = .4d0 ! von karman constant
      real*8, parameter :: z=10. ! level at which air t,q,u,v prescribed
      real*8 :: cd_n10, cq_n10, ch_n10, cd_n10_rt ! neutral 10m drag coefficients
      real*8 :: zeta, x2, x, psi_m, psi_h         ! stability parameters
      real*8 :: w10, ustar, bstar, tstar, qstar, z0, xx, stab, cd_rt
      integer, parameter :: n_itts = 2

C****
      INTEGER :: J_0, J_1, I_0,I_1

      real*8, dimension(:,:), pointer :: rsi,msi,snowi,pond_melt
     &     ,sss,focean,snoage
      real*8, dimension(:,:,:), pointer :: ssi
      logical, dimension(:,:), pointer :: flag_dsws

      rsi => si_ocn%rsi
      msi => si_ocn%msi
      snowi => si_ocn%snowi
      ssi => si_ocn%ssi
      pond_melt => si_ocn%pond_melt
      flag_dsws => si_ocn%flag_dsws
      snoage => atmice%snoage !si_ocn%snoage
      focean => atmocn%focean
      sss => atmocn%sss

      call getDomainBounds(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      call seaice_to_atmgrid(atmice)

C**** Zero out fluxes
      atmocn%E0=0. ; atmocn%EVAPOR=0. ; atmocn%DMUA=0. ; atmocn%DMVA=0.
      atmice%E0=0. ; atmice%EVAPOR=0. ; atmice%DMUA=0. ; atmice%DMVA=0.
      atmice%E1=0.
      atmocn%SOLAR=0.; atmice%SOLAR=0.

      DO J=J_0,J_1
      DO I=I_0,atmocn%IMAXJ(J)
      IF(FOCEAN(I,J).LE.0.) CYCLE

      ! using ocean-seaice composite values
      us = atmocn%usavg(i,j)
      vs = atmocn%vsavg(i,j)
      ts = atmocn%tsavg(i,j)
      qs = atmocn%qsavg(i,j)
      ps = atmocn%srfp(i,j)*100.

      ws = sqrt(us*us+vs*vs)
      atmocn%wsavg(i,j) = ws

C****
C**** DETERMINE SURFACE CONDITIONS
C****
      POICE=RSI(I,J)
      POCEAN=1.-POICE

      if(poice.gt.0.) then
        IF(LAT2D(I,J).lt.0.) THEN
          HEMI = 1
        ELSE
          HEMI = 2
        ENDIF
        zoice = (ace1i+msi(i,j))/rhoi
        zsnow = snowi(i,j)/rhos
        fmp = min(1.6d0*sqrt(pond_melt(i,j)/rhow),1d0)
        zmp = min(0.8d0*fmp,0.9d0*zoice)
        dalbsn = 0.
      endif
      CALL GET_SURF_ALBEDO(
     i     atmocn%COSZ1(I,J),
     i     HEMI,
     i     SNOAGE(I,J),POCEAN,POICE,
     i     ZOICE,fmp,zsnow,zmp,
     i     SNOWI(I,J),WS,dalbsn,
     i     flag_dsws(i,j),
     o     OCNALB,SIALB
     &     )

#if (defined CONSTANT_OCEAN_ALBEDO)
      OCNALB = ocean_albedo_const
#elif (defined CORE2_OCEAN_ALBEDO)
      OCNALB = .069d0 - .011d0*cos(2d0*lat2d(i,j))
#endif

      SRHEAT = atmocn%FSHORT(I,J)*atmocn%COSZ1(I,J)/
     &     (COSZ_DAY(I,J)+teeny)

      aij(i,j,ij_foc) = aij(i,j,ij_foc) + 1.
      aij(i,j,ij_sst) = aij(i,j,ij_sst) + atmocn%GTEMP(I,J)
      aij(i,j,ij_sss) = aij(i,j,ij_sss) + SSS(I,J)

      DO ITYPE=1,2

      if ( ITYPE == 1 ) then
        PTYPE=POCEAN
      else
        PTYPE=POICE
      endif
      if(ptype.le.0.) cycle


      if ( ITYPE == 1 ) then
C****
C**** OPEN OCEAN
C****
        TG1 = atmocn%GTEMP(I,J)
        TG2 = atmocn%GTEMP2(I,J)      ! diagnostic only
        uocean = atmocn%uosurf(i,j)
        vocean = atmocn%vosurf(i,j)
        aij(i,j,ij_roc) = aij(i,j,ij_roc) + ptype
        aij(i,j,ij_ocnsw) = aij(i,j,ij_ocnsw) + ptype*srheat
        aij(i,j,ij_ocnalb) = aij(i,j,ij_ocnalb) + ptype*srheat*ocnalb
        SRHEAT = (1.-OCNALB)*SRHEAT
        atmocn%SOLAR(I,J)=DTSRC*SRHEAT
        ELHX=LHE
c**** sanity check (to prevent rare anomalies that will be dealt with by
C**** addice next time)
        TG1=max(TG1,tfrez(sss(i,j)))

      else
C****
C**** OCEAN ICE
C****
        uocean = atmice%uisurf(i,j)
        vocean = atmice%visurf(i,j)
        TG1=atmice%GTEMP(I,J)
        TG2=atmice%GTEMP2(I,J)
        SNOW=SNOWI(I,J)
        MSI1=SNOW+ACE1I         ! snow and first layer ice mass (kg/m^2)
        MSI2=MSI(I,J)           ! second (physical) layer ice mass (kg/m^2)
C**** determine heat capacity etc for top ice layers
        dF1dTG = 2./(ACE1I/(RHOI*alami(TG1,1d3*((SSI(1,I,J)+SSI(2,I,J))
     *       /ACE1I)))+SNOW*BYRLS)
        HCG1 = dEidTi(TG1,1d3*(SSI(1,I,J)/(XSI(1)*MSI1)))*XSI(1)*MSI1
        HCG2 = dEidTi(TG2,1d3*(SSI(2,I,J)/(XSI(2)*MSI1)))*XSI(2)*MSI1
        aij(i,j,ij_rsi) = aij(i,j,ij_rsi) + ptype
        aij(i,j,ij_sisw) = aij(i,j,ij_sisw) + ptype*srheat
        aij(i,j,ij_sialb) = aij(i,j,ij_sialb) + ptype*srheat*sialb
        aij(i,j,ij_msi) = aij(i,j,ij_msi) + ptype*(ace1i+msi2)
        SRHEAT = (1.-SIALB)*SRHEAT
        atmice%SOLAR(I,J)=DTSRC*SRHEAT
C**** fraction of solar radiation leaving layer 1 and 2
        IF (SRHEAT.gt.0) THEN   ! only bother if there is sun
          call solar_ice_frac(SNOW,MSI2,FLAG_DSWS(I,J),FSRI,2)
        ELSE
          FSRI(1:2) = 0
        END IF
        ELHX=LHS

      endif

C****
C**** BOUNDARY LAYER INTERACTION
C****
      SHDT=0.
      EVHDT=0.
      TRHDT=0.
      F1DT=0.

      tg = tg1+tf
      QG = QSAT(TG,ELHX,.01d0*PS)
      IF (ITYPE.eq.1) QG=0.98d0*QG

      tsv = ts*(1+deltx*qs)
      ws = max(ws, 0.5)  ! 0.5 m/s floor on wind (undocumented NCAR)
      w10 = ws              ! first guess 10m wind
    
      cd_n10 = (2.7/w10+0.142+0.0764*w10)/1e3  ! L-Y eqn. 6a
      cd_n10_rt = sqrt(cd_n10)
      cq_n10 =  34.6 *cd_n10_rt/1e3            ! L-Y eqn. 6b
      stab = 0.5d0 + sign(0.5d0,ts-tg)
      ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3  ! L-Y eqn. 6c
  
      cd = cd_n10     ! first guess for exchange coeff's at z
      ch = ch_n10
      cq = cq_n10
      do iter=1,n_itts      ! Monin-Obukhov iteration
        cd_rt = sqrt(cd)
        ustar = cd_rt*ws                             ! L-Y eqn. 7a
        tstar    = (ch/cd_rt)*(ts-tg)                ! L-Y eqn. 7b
        qstar    = (cq/cd_rt)*(qs-qg)                ! L-Y eqn. 7c
        bstar = grav*(tstar/tsv+qstar/(qs+1/deltx))
        zeta     = vonkarm*bstar*z/(ustar*ustar)    ! L-Y eqn. 8a
        zeta     = sign( min(abs(zeta),10.0), zeta )  ! undocumented NCAR
        x2 = sqrt(abs(1-16*zeta))                     ! L-Y eqn. 8b
        x2 = max(x2, 1.0)                             ! undocumented NCAR
        x = sqrt(x2)
    
        if (zeta > 0) then
          psi_m = -5*zeta                             ! L-Y eqn. 8c
          psi_h = -5*zeta                             ! L-Y eqn. 8c
        else
          psi_m = log((1+2*x+x2)*(1+x2)/8)-2*(atan(x)-atan(1.0)) ! L-Y eqn. 8d
          psi_h = 2*log((1+x2)/2)                                ! L-Y eqn. 8e
        end if
    
        w10 = ws/(1+cd_n10_rt*(log(z/10)-psi_m)/vonkarm)   ! L-Y eqn. 9
        cd_n10 = (2.7/w10+0.142+0.0764*w10)/1e3           ! L-Y eqn. 6a again
        cd_n10_rt = sqrt(cd_n10)
        cq_n10 = 34.6*cd_n10_rt/1e3                       ! L-Y eqn. 6b again
        stab = 0.5d0 + sign(0.5d0,zeta)
        ch_n10 = (18.0*stab+32.7*(1-stab))*cd_n10_rt/1e3  ! L-Y eqn. 6c again
        z0 = 10*exp(-vonkarm/cd_n10_rt)                   ! diagnostic
    
        xx = (log(z/10)-psi_m)/vonkarm
        cd = cd_n10/(1+cd_n10_rt*xx)**2                   ! L-Y 10a
        xx = (log(z/10)-psi_h)/vonkarm
!!$        ch = ch_n10/(1+ch_n10*xx/cd_n10_rt)**2                !       b (bug)
!!$        cq = cq_n10/(1+cq_n10*xx/cd_n10_rt)**2                !       c (bug)
        ! 10b, 10c (corrected code aug2007)
        ch = ch_n10/(1+ch_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10) ! 10b 
        cq = cq_n10/(1+cq_n10*xx/cd_n10_rt)*sqrt(cd/cd_n10) ! 10c
      end do
      cm = cd

C**** CALCULATE RHOSRF*CM*WS AND RHOSRF*CH*WS
      RHOSRF=PS/(RGAS*TG)

      RCDMWS=CM*WS*RHOSRF
      RCDHWS=CH*WS*RHOSRF
      RCDQWS=CQ*WS*RHOSRF
C**** CALCULATE FLUXES OF SENSIBLE HEAT, LATENT HEAT, THERMAL
C****   RADIATION, AND CONDUCTION HEAT (WATTS/M**2) (positive down)
      SHEAT=SHA*(RCDHWS*(TS-TG))
      EVHEAT=(LHE+TG1*SHV)*(RCDQWS*(QS-QG))
      TR4=(TG1+TF)**4
      TRHEAT = atmocn%FLONG(I,J)-STBO*TR4

C**** CASE (1) ! FLUXES USING EXPLICIT TIME STEP FOR OCEAN POINTS
      if ( ITYPE == 1) then
        SHDT = DTSRC*SHEAT
        EVHDT=DTSRC*EVHEAT              ! latent heat flux
        TRHDT=DTSRC*TRHEAT

C**** CASE (2) ! FLUXES USING IMPLICIT TIME STEP FOR ICE POINTS
      else if ( ITYPE == 2 ) then

! heat flux on first/second/third layers (W/m^2)
        F1 = (TG1-TG2)*dF1dTG + SRHEAT*FSRI(1)
        F2 = SRHEAT*FSRI(2)
        EVHEAT=LHE*(RCDQWS*(QS-QG)) ! why is this different to above?
        F0=SRHEAT+TRHEAT+SHEAT+EVHEAT
        dSNdTG=-RCDHWS*SHA
        dQGdTG = QG*DQSATDT(TG,ELHX) ! d(QG)/dTG
        dEVdTG = -dQGdTG*LHE*RCDQWS ! d(EVHEAT)/dTG
        dTRdTG = -4*STBO*sqrt(sqrt(TR4))**3 ! d(TRHEAT)/dTG
        dF0dTG = dSNdTG+dEVdTG+dTRdTG ! d(F0)/dTG

        T2DEN = HCG2+DTSRC*dF1dTG
        T2CON = DTSRC*(F1-F2)/T2DEN
        T2MUL = DTSRC*dF1dTG/T2DEN

        DFDTG=DF0DTG-dF1dTG
        DTG=(F0-F1)*DTSRC/(HCG1-DTSRC*DFDTG)

        IF (TG1+dTG .GT. 0.) dTG = -TG1
        dT2 = T2CON+T2MUL*dTG

        SHDT  = DTSRC*(SHEAT +dTG*dSNdTG) ! sensible
        EVHDT = DTSRC*(EVHEAT+dTG*dEVdTG) ! latent
        TRHDT = DTSRC*(TRHEAT+dTG*dTRdTG) ! thermal flux (J/m^2)
        F1DT = DTSRC*(F1+(dTG*dF1dTG-dT2*dF1dTG))
        TG1 = TG1+dTG          ! first layer sea ice temperature (degC)
        TG2 = TG2+dT2          ! second layer sea ice temperature (degC)

      endif

C**** CALCULATE EVAPORATION
      EVAP =-EVHDT/((LHE+TG1*SHV))

C**** ACCUMULATE SURFACE FLUXES AND PROGNOSTIC AND DIAGNOSTIC QUANTITIES
      F0DT=DTSRC*SRHEAT+TRHDT+SHDT+EVHDT

      if(itype.eq.1) then
        atmocn%E0(I,J)=F0DT
        atmocn%EVAPOR(I,J)=EVAP
        atmocn%DMUA(I,J)=DTSRC*RCDMWS*(US-UOCEAN)
        atmocn%DMVA(I,J)=DTSRC*RCDMWS*(VS-VOCEAN)
        aij(i,j,ij_ocheat) = aij(i,j,ij_ocheat) + ptype*F0DT
      else
        atmice%E0(I,J)=F0DT
        atmice%E1(I,J)=F1DT
        atmice%EVAPOR(I,J)=EVAP
        atmice%DMUA(I,J)=DTSRC*RCDMWS*(US-UOCEAN)
        atmice%DMVA(I,J)=DTSRC*RCDMWS*(VS-VOCEAN)
      endif

      ENDDO   ! end of itype loop
      ENDDO   ! end of I loop
      ENDDO   ! end of J loop

      do j=j_0,j_1
      do i=i_0,atmocn%imaxj(j)
        if(focean(i,j).le.0.) cycle
        poice = rsi(i,j)
        pocean = 1d0-poice
        taux = (atmocn%dmua(i,j)*pocean+atmice%dmua(i,j)*poice)/dtsrc
        tauy = (atmocn%dmva(i,j)*pocean+atmice%dmva(i,j)*poice)/dtsrc
        aij(i,j,ij_tauus) = aij(i,j,ij_tauus) + taux
        aij(i,j,ij_tauvs) = aij(i,j,ij_tauvs) + tauy
        aij(i,j,ij_taus) = aij(i,j,ij_taus) + sqrt(taux**2 + tauy**2)
      enddo
      enddo

      RETURN
      END SUBROUTINE SURFACE

#ifdef STANDALONE_HYCOM
      ! hycom is responsible for its own salinity restoring
#else
      subroutine restore_surface_salinity!(atmocn)
      !use exchange_types, only : atmocn_xchng_vars
      use fluxes, only : atmocn
      use diag_com, only : aij,ij_salres
      use ocean, only : imo=>im,jmo=>jm,imaxj
     &     ,mo,s0m,focean,dxypo
      use domain_decomp_1d, only : getDomainBounds,globalsum,am_i_root
      use oceanr_dim, only : ogrid
      implicit none
      !type(atmocn_xchng_vars) :: atmocn
c
      integer :: i,j, j_0, j_1
      logical :: have_south_pole,have_north_pole
      real*8 :: restore,sal,ds,sumneg,sumpos,negfac,posfac
      real*8, dimension(imo,ogrid%j_strt_halo:ogrid%j_stop_halo) ::
     &     dsneg,dspos


      call getDomainBounds(ogrid, j_strt=j_0, j_stop=j_1,
     &     have_south_pole=have_south_pole,
     &     have_north_pole=have_north_pole)
      dsneg(:,:) = 0.
      dspos(:,:) = 0.
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(focean(i,j).le.0.) cycle
        restore = atmocn%sssresfac(i,j)
        sal = s0m(i,j,1)/(dxypo(j)*mo(i,j,1))
        sal = sal*restore + atmocn%sssobs(i,j)*(1.-restore)
        ds = sal*(dxypo(j)*mo(i,j,1)) - s0m(i,j,1)
        if(ds.gt.0.) then
          dspos(i,j) = ds
        else
          dsneg(i,j) = ds
        endif
      enddo
      enddo
      if(have_south_pole) then
        dsneg(2:imo,1) = dsneg(1,1)
        dspos(2:imo,1) = dspos(1,1)
      endif
      if(have_north_pole) then
        dsneg(2:imo,jmo) = dsneg(1,jmo)
        dspos(2:imo,jmo) = dspos(1,jmo)
      endif
      call globalsum(ogrid,dsneg,sumneg,all=.true.)
      call globalsum(ogrid,dspos,sumpos,all=.true.)
      if(abs(sumneg).gt.sumpos) then
        posfac = 1.
        negfac = -sumpos/sumneg
      else
        negfac = 1.
        posfac = -sumneg/sumpos
      endif
      if(am_i_root()) write(6,*) 'negfac,posfac ',negfac,posfac
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(focean(i,j).le.0.) cycle
        ds = dsneg(i,j)*negfac+dspos(i,j)*posfac
        s0m(i,j,1) = s0m(i,j,1) + ds
        aij(i,j,ij_salres) = aij(i,j,ij_salres) +ds/dxypo(j)
      enddo
      enddo
      return
      end subroutine restore_surface_salinity

      subroutine restore_surface_salinity2!(atmocn)
      !use exchange_types, only : atmocn_xchng_vars
      use fluxes, only : atmocn
      use diag_com, only : aij,ij_salres
      use ocean, only : imo=>im,jmo=>jm,imaxj,lmm
     &     ,mo,s0m,focean,dxypo
      use domain_decomp_1d, only : getDomainBounds,globalsum,am_i_root
      use oceanr_dim, only : ogrid
      implicit none
      !type(atmocn_xchng_vars) :: atmocn
c
      integer :: i,j,l, j_0, j_1
      logical :: have_south_pole,have_north_pole
      real*8 :: restore,sal,salrat,sums,sumds
      real*8, dimension(imo,ogrid%j_strt_halo:ogrid%j_stop_halo) ::
     &     ds,ssum

      call getDomainBounds(ogrid, j_strt=j_0, j_stop=j_1,
     &     have_south_pole=have_south_pole,
     &     have_north_pole=have_north_pole)
      ds(:,:) = 0.
      ssum(:,:) = 0.
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(focean(i,j).le.0.) cycle
        restore = atmocn%sssresfac(i,j)
        sal = s0m(i,j,1)/(dxypo(j)*mo(i,j,1))
        sal = sal*restore + atmocn%sssobs(i,j)*(1.-restore)
        ds(i,j) = sal*(dxypo(j)*mo(i,j,1)) - s0m(i,j,1)
        ssum(i,j) = sum(s0m(i,j,1:lmm(i,j)))
      enddo
      enddo
      if(have_south_pole) then
        ds(2:imo,1) = ds(1,1)
        ssum(2:imo,1) = ssum(1,1)
      endif
      if(have_north_pole) then
        ds(2:imo,jmo) = ds(1,jmo)
        ssum(2:imo,jmo) = ssum(1,jmo)
      endif
      call globalsum(ogrid,ds,sumds,all=.true.)
      call globalsum(ogrid,ssum,sums,all=.true.)
      salrat = sums/(sums+sumds)
      if(am_i_root()) write(6,*) 'salrat ',salrat-1.
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(focean(i,j).le.0.) cycle
        s0m(i,j,1) = s0m(i,j,1) + ds(i,j)
        ssum(i,j) = ssum(i,j) + ds(i,j)
        do l=1,lmm(i,j)
          s0m(i,j,l) = s0m(i,j,l)*salrat
        enddo
        aij(i,j,ij_salres) = aij(i,j,ij_salres)
     &       +(ds(i,j)+ssum(i,j)*(salrat-1.))/dxypo(j)
      enddo
      enddo
      return
      end subroutine restore_surface_salinity2
#endif
