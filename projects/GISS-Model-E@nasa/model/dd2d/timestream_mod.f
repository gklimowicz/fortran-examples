#include "rundeck_opts.h"

!@sum  module timestream_mod is a package to read and time-interpolate
!@+    time-varying gridded data from suitably structured input files.
!@auth M. Kelley
!@ver  beta.  Loose ends here and there.
!@usage
!@+    Notes:
!@+
!@+    - Instantiation of a timestream object via init_stream
!@+    - Obtaining values via read_stream (or get_by_index)
!@+    - Time-registration conventions for datafiles
!@+    - Miscellaneous
!@+
!@+    Instantiation of a timestream object via init_stream
!@+
!@+        Use of this object (fortran derived type) as a file handle
!@+        facilitates the caching of datafile values and other
!@+        information required for efficient time interpolation
!@+        and disk access, and also provides a pass-able container
!@+        for the set of parameters defining the time interpolation:
!@+        time interpolation methods, min/max valid values, imposed
!@+        periodicity, etc.
!@+
!@+        call init_stream(grid,tstream,fbase,vname,qmin,qmax,method,
!@+                         jyear,jday,nskip,msk,cyclic)
!@+                        grid: instance of dist_grid corresponding to
!@+                              the grid on which data will be requested
!@+                     tstream: timestream object to be initialized
!@+             character fbase: path to input data (file, directory, or
!@+                              symbolic link thereto)
!@+             character vname: netcdf variable to be read from fbase
!@+            character method: desired intra-annual time interpolation method
!@+                              (valid values: 'ppm', 'linm2m', 'none')
!@+            real*8 qmin,qmax: min/max valid values for interpolated output
!@+          integer jyear,jday: date info for pre-fetching data needed
!@+                              for first time interpolations
!@+               integer nskip: optional, unused, only applies to giss-fmt files
!@+                  real*8 msk: optional 0/1 mask array defining the set of
!@+                              gridpoints that read_stream should fill
!@+              logical cyclic: optional, whether to impose annual periodicity
!@+                              (only matters for multi-year input data)
!@+        
!@+    Obtaining values via read_stream (or get_by_index)
!@+        
!@+        read_stream() is the primary access method.
!@+
!@+        call read_stream(grid,tstream,jyear,jday,arr)
!@+                        grid: instance of dist_grid
!@+                     tstream: timestream object
!@+          integer jyear,jday: date info
!@+                  real*8 arr: output array (2D i,j or 3D i,j,k or 4D i,j,k,l)
!@+
!@+        Time interpolation may be performed to obtain the output array,
!@+        and partial reads from the datafile occur as necessary.  The
!@+        sequence of dates for successive read_stream calls need not be
!@+        increasing or monotonic, but extra I/O costs are incurred when
!@+        going backward or when jumping far ahead in time.
!@+        Time interpolation is performed on two timescales:
!@+        inter-annual (when needed) and intra-annual.  The former permits
!@+        the use of "climatologies" defined over multi-year periods,
!@+        and the latter permits creation of daily time resolution
!@+        values from longer-period data.  The inter-annual case
!@+        currently only employs linear interpolation, while the
!@+        intra-annual methods include the range of possibilities
!@+        accepted by init_stream.
!@+
!@+        In some circumstances, a calling program may need to refer
!@+        directly to sequences of elements of the datafile. Extant
!@+        examples of this access in Model E include references
!@+        to a full annual cycle of monthly values for some year.
!@+        While calls to read_dist_data(...,record[1]=N) are
!@+        the most straightforward way to access the values, such
!@+        accesses may appear in conjunction with calls to init_stream
!@+        or read_stream.  Since the timestream object currently
!@+        caches a full year of values, the get_by_index()
!@+        procedure was introduced to provide access
!@+        to those data with no additional I/O cost.
!@+        
!@+        call get_by_index(grid,tstream,ind,arr)
!@+                  grid: instance of dist_grid
!@+               tstream: timestream object
!@+           integer ind: index of desired input-file timestep (1-based)
!@+                        for the last year passed to init/read_stream
!@+            real*8 arr: the output array (2D i,j or 3D i,j,k or 4D i,j,k,l)
!@+        
!@+    Time-registration conventions for datafiles
!@+
!@+        Time-dependent variables in datafiles must possess the unlimited
!@+        dimension, though this dimension does not have to be named "time".
!@+        The type and origin of the time axis in the datafile(s) need
!@+        not be known in advance by the calling program, as it is either
!@+        declared in the input metadata or implicitly specified, as per
!@+        the following logic applied during init_stream().  Case 2b
!@+        is likely to be the most common, though splitting timeseries
!@+        into multiple files is sometimes advantageous (case 1).
!@+        
!@+        (1) The path is determined to be a directory. It is searched
!@+            for filenames of the form YYYY.nc, each assumed to contain
!@+            one full year of data, where YYYY is the 4-digit
!@+            calendar year (e.g. 1950).  The time resolution of the data
!@+            is inferred from the size of the unlimited dimension
!@+            in the file corresponding to the current year, or the next
!@+            existing one in the dataset.  A size of 1 denotes annual,
!@+            12 monthly, 73 pentads, 365/366 daily, and so forth.
!@+        (2) The path is a regular file.  Either
!@+            (a) It lacks the time axis information defined below in (b), in
!@+                which case it is assumed to correspond to one year of data
!@+                whose time resolution is determined as in (1) above
!@+            (b) It possesses a time coordinate, i.e. a 1-dimensional
!@+                variable whose sole dimension is the unlimited dimension
!@+                and whose name is the name of the unlimited dimension.
!@+                The "units" attribute of this variable must match one of
!@+                the following patterns (case-sensitive):
!@+                   "years since YYYY"       for annual data
!@+                   "months since YYYY-MM"   for monthly data
!@+                   "pentads since YYYY-PP"  for pentad data
!@+                   "days since YYYY-MM-DD"  for daily data
!@+                where YYYY denotes a 4-digit year, MM or PP a 2-digit month
!@+                or pentad, and DD a 2-digit day of the month.
!@+                Units of hours and minutes
!@+                will be added to the match list in the future as necessary.
!@+                Currently, only the integer part of the time axis values is
!@+                considered. Excepting the first/last, each year of data must
!@+                be complete. Therefore t(n+1)-t(n)-1 must be a non-negative
!@+                integer multiple of the number of steps per year, where t(n)
!@+                is the time value for step n.  If the first year is incomplete,
!@+                it must immediately precede the second year.
!@+
!@+    Miscellaneous
!@+
!@+        getname_firstfile (see its comments) will return the first file
!@+        in a YYYY.nc sequence (time-registration convention #1 above).
!@+
!@+        There is no requirement that all years exist in a timeseries,
!@+        or that the years be evenly spaced.
!@+
!@+        When performed, intra-annual time interpolation is to daily
!@+        resolution.  Interpolation to sub-daily resolution has not
!@+        yet been implemented.  It will be straightforward to do so
!@+        once it is required in Model E for nudged runs, high-frequency
!@+        inputs to single-component "standalone" runs, etc.
!@+        
!@+        The rank of input data can be either 2 or 3 (2 horizontal
!@+        dimensions, one optional 3rd dimension which must the the
!@+        last).  Support for rank-4 arrays may be added in the future.
!@+        
!@+        Files are assumed to be in netcdf format.  An early version of
!@+        this module supported GISS format for 2D data, but that effort
!@+        was postponed pending evaluation of actual interest in the option.
!@+
!@+        The horizontal grid of the input data must match that of the
!@+        calling program.  On-the-fly horizontal remapping may be added
!@+        in the future.
!@+
!@+        Calling programs expose two explicit assumptions about the names
!@+        of input data:
!@+           (1) the path (file or directory name, or symlink thereto)
!@+           (2) the name of the netcdf variable
!@+        If the netcdf variable name is not found in the file, the
!@+        list of global attributes in the file is queried for the
!@+        existence of an attribute name XYZname, where XYZ is the
!@+        netcdf variable name passed from the calling program.  The
!@+        value of XYZname is then assumed to be the variable name.
!@+        This allows users to maintain the original variable names
!@+        in input files, simply adding a global attribute to those
!@+        files as necessary.
!@+        Example: the calling program requests variable "sst", but
!@+          this variable is named "temp" in the datafile.  Add a
!@+          global attribute sstname="temp" to the datafile to use this
!@+          feature.
!@+
!@+        A datafile may contain multiple time-dependent variables,
!@+        but each variable to be read from the file requires its own
!@+        instance of a timestream object. The reading of a given
!@+        variable from a netcdf file is formally independent of the
!@+        presence of other variables in that file.
!@+        
!@+        If the intra-annual time interpolation uses PPM, end-of-period
!@+        values may be read from the datafile rather than calculated
!@+        on the fly from period means.  This option is automatically
!@+        activated for a timestream having pathname XYZ and variable
!@+        name xyz if there also exists a file or directory XYZ_eom
!@+        containing a netcdf variable xyz_eom. XYZ_eom may be a symbolic
!@+        link to the same data as XYZ.
!@+        
!@+        Datafiles are not kept open between calls to read_stream().
!@+

      module timestream_mod
      IMPLICIT NONE
      SAVE
      PRIVATE

!@var giss_fmt,netcdf_fmt format constants
      integer, parameter ::
     &     netcdf_fmt=1, gissclim_fmt=2, fbsa_fmt=3

      integer, parameter :: linm2m=1, ppm=2, nointerp=3

      public :: timestream
      public :: init_stream,read_stream,get_by_index
     &     ,reset_stream_properties,getname_firstfile

      interface read_stream
        module procedure read_stream_2d
        module procedure read_stream_3d
        module procedure read_stream_4d
      end interface read_stream

      interface get_by_index
        module procedure get_by_index_2d
        module procedure get_by_index_3d
        module procedure get_by_index_4d
      end interface get_by_index

      integer, parameter :: firstcall_int=-9999,reset_int=999999

!@type timestream derived type serving as a file handle and container
!@+    for time interpolation parameters and cached input data
      type timestream
!@var ndims number of dimensions of outputs (excluding time)
         integer :: ndims=2
!@var dlens lengths of non-horizontal dimensions (excluding time)
         integer :: dlens(2)=(/1,1/)
!@var lm product(dlens)
         integer :: lm=1
!@var m2r number of time steps in one year of data
         integer :: m2r=1
!@var fmt flag taking a value equal to giss_fmt or netcdf_fmt.
         integer :: fmt=-9999
!@var tinterp_method type of intra-annual time interpolation
         integer :: tinterp_method
!@var eom_from_file flag indicating whether netcdf input file
!@+   contains end-of-month values (automatically determined in read_netcdf)
         logical :: eom_from_file=.false.
!@var cyclic impose that the first year of data is used subsequently
         logical :: cyclic=.false.
!@var daily_data flag indicating whether input file contains daily data
!@var annual_data flag indicating whether input file contains annual data
!@var monthly_data flag indicating whether input file contains monthly data
!@var pentad_data flag indicating whether input file contains pentad data
         logical :: daily_data=.false.
         logical :: annual_data=.false.
         logical :: monthly_data=.false.
         logical :: pentad_data=.false.
!@var multiple_files whether multi-year data are stored one file per year
         logical :: multiple_files=.false.
!@var roff time index corresponding to January of first year, minus 1
         integer :: roff=0
!@var fbase path to data
         character(len=32) :: fbase
!@var firstfile path to first file containing data
!@+    (==fbase if not one file per year)
         character(len=32) :: firstfile
!@var vname netcdf name for qty in the input file
         character(len=32) :: vname
!@var cycl flag indicating whether data repeat each year 0:no 1:yes
!@+   (automatically determined from the input files)
         integer :: cycl
!@var qmin,qmax minimum/maximum valid values for interpolated output
         real*8 :: qmin,qmax
!@var fid file ID
         integer :: fid
!@var nskip unused for netcdf datafiles
         integer :: nskip=0
!@var year_sv year corresponding to current contents of qty
         integer :: year_sv=firstcall_int
!@var year_sv2 year corresponding to current contents of qty2
!@+   (if no interannual interpolation, year_sv2=year_sv)
         integer :: year_sv2=firstcall_int
!@var npad padding size at beginning/end of year
         integer :: npad=0
!@var nfileyrs number of years of available data on disk
         integer :: nfileyrs=0
!@var fileyrs list of years of available data
         integer, dimension(:), allocatable :: fileyrs

!@var msk 0/1 mask array defining the set of gridpoints to fill
!@+   during time interpolation
         real*8, dimension(:,:), allocatable :: msk
!@var qty holds one year of data
!@+   The unpadded size of its time dimension is m2r.  If intra-annual
!@+   time interpolation is performed, the size of the time dimension
!@+   is padded by 2 at the beginning and end (qty allocated over the
!@+   range -1:m2r+2)
         real*8, dimension(:,:,:,:), allocatable :: qty
!@var qty1, qty2 when interannual interpolation is required, these
!@+   arrays hold the data for the years bracketing the current year.
!@+   They are unpadded
         real*8, dimension(:,:,:,:), allocatable :: qty1,qty2
!@var eom holds one year of end-of-month values needed for
!@+   parabolic interpolation (13 months incl. prev. Dec.)
         real*8, dimension(:,:,:,:), allocatable :: eom
      end type timestream

C**** (Simplified) Calendar Related Terms
!@param JDperY,JMperY    number of days,months per year
!@var   JDendOfM(0:12)   last Julian day in month
!@var   JDmidOfM(0:13)   middle Julian day in month
! These calendar arrays have been copied from MODEL_COM to allow
! this module to exist in a component subdirectory.  Bad practice
! but not a big deal yet.
      integer, PARAMETER :: JDPERY = 365, JMPERY = 12
      integer :: JDendOfM(0:JMPERY) = (
     *     /0,31,59,90,120,151,181,212,243,273,304,334,365/)
      integer :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,45,75,106,136,167,197,228,259,289,320,350,381/)

      contains

      subroutine init_stream(grid,tstream,fbase,vname,qmin,qmax,method,
     &     jyear,jday,nskip,msk,cyclic
     &     )
!@sum init_stream initializes an instance of the timestream object.
!@+   see usage notes at beginning of module
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      character(len=*) fbase,vname,method
      real*8 :: qmin,qmax
      integer :: jyear,jday
      integer, optional :: nskip
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo), optional ::
     &     msk
      logical, optional :: cyclic
c
      integer :: jmon
      integer :: i_0,i_1,j_0,j_1
c
      tstream%fbase = fbase
      tstream%vname = vname
      tstream%qmin = qmin
      tstream%qmax = qmax

      select case (method)
      case ('ppm')
        tstream%tinterp_method = ppm
      case ('linm2m')
        tstream%tinterp_method = linm2m
      case ('none')
        tstream%tinterp_method = nointerp
      case default
        call stop_model(
     &       'init_stream: unrecognized time interp method',255)
      end select

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      allocate(tstream%msk(i_0:i_1,j_0:j_1))
      if(present(msk)) then
        tstream%msk = msk(i_0:i_1,j_0:j_1)
      else
        tstream%msk = 1d0
      endif

      call check_format(tstream)

      ! duplicate of read_stream. maybe pass jday to read_year instead
      jmon=1
      do while(jday.gt.JDendOfM(jmon))
        jmon=jmon+1
      enddo

      if(present(cyclic)) tstream%cyclic = cyclic

      call read_year(grid,tstream,jyear,jmon)

      return
      end subroutine init_stream

      subroutine reset_stream_properties(grid,tstream,cyclic)
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      logical, optional :: cyclic
c
      if(present(cyclic)) then
        tstream%cyclic = cyclic
        tstream%year_sv = reset_int
      endif
c
      return
      end subroutine reset_stream_properties

      subroutine getname_firstfile(tstream,fname)
!@sum getname_firstfile
!@+   Many timeseries contain additional metadata that must be queried
!@+   during initialization.   This is often done by independently opening
!@+   and closing the path passed to init_stream().   For time-registration
!@+   convention #1 (YYYY.nc grouping), this method requires obtaining
!@+   the name of an existing YYYY.nc file in the input directory.   The
!@+   routine getname_firstfile can be used to obtain the first such file
!@+   (smallest value of YYYY) for a given timestream object.   If the
!@+   path for that object is a file instead of a directory, the name of
!@+   that file is returned.
      implicit none
      type(timestream) :: tstream
!@var fname returned filename
      character(len=*) :: fname
      fname = tstream%firstfile
      end subroutine getname_firstfile

      subroutine get_by_index_2d(grid,tstream,ind,arr)
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: ind
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: arr
      if(ind.lt.1 .or. ind.gt.size(tstream%qty,4)) then
        call stop_model('get_by_index: bad_index',255)
      endif
      arr(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop) =
     &     tstream%qty(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop
     &     ,1,ind)
      return
      end subroutine get_by_index_2d

      subroutine get_by_index_3d(grid,tstream,ind,arr)
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: ind
      real*8, dimension(grid%i_strt_halo:,
     &                  grid%j_strt_halo:,:) :: arr
      if(ind.lt.1 .or. ind.gt.size(tstream%qty,4)) then
        call stop_model('get_by_index: bad_index',255)
      endif
      arr(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,:) =
     & tstream%qty(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop
     &     ,:,ind)
      return
      end subroutine get_by_index_3d

      subroutine get_by_index_4d(grid,tstream,ind,arr)
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: ind,k,l,kl
      real*8, dimension(grid%i_strt_halo:,
     &                  grid%j_strt_halo:,:,:) :: arr
      if(ind.lt.1 .or. ind.gt.size(tstream%qty,4)) then
        call stop_model('get_by_index: bad_index',255)
      endif
      kl = 0
      do l=1,tstream%dlens(2)
      do k=1,tstream%dlens(1)
        kl = kl + 1
        arr(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,k,l) =
     &      tstream%qty(grid%i_strt:grid%i_stop,grid%j_strt:grid%j_stop,
     &       kl,ind)
      enddo
      enddo
      return
      end subroutine get_by_index_4d

      subroutine read_year(grid,tstream,jyear,jmon)
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: jyear,jmon
!!       if(tstream%fmt.lt.0) call check_format(tstream)
!!       if(tstream%fmt.eq.netcdf_fmt) then
        call read_stream_netcdf(grid,tstream,jyear,jmon)
!!       else
!!         call read_stream_gissfmt(grid,tstream,jyear,jmon)
!!       endif
      end subroutine read_year

      subroutine check_format(tstream,nrecs)
!!       ! determine whether files are GISS- or netcdf-format
!!       USE FILEMANAGER
      implicit none
      type(timestream) :: tstream
      integer, optional :: nrecs
!!       integer, allocatable :: reclens(:)
!!       integer :: nrecs_
!!       if(present(nrecs)) then
!!         nrecs_ = nrecs
!!       else
!!         nrecs_ = 99 ! reasonable default
!!       endif
!!       allocate(reclens(nrecs_))
!!       call get_recordlengths(tstream%fbase,nrecs_,reclens)
!!       deallocate(reclens)
!!       if(nrecs_.le.0) then
        tstream%fmt = netcdf_fmt
!!       else
!!         if(nrecs_.eq.12+tstream%nskip) then
!!           tstream%fmt = gissclim_fmt
!!         else
!!           tstream%fmt = fbsa_fmt
!!         endif
!!       endif
!!       if(present(nrecs)) nrecs = nrecs_
      if(present(nrecs)) nrecs = 0
      end subroutine check_format

      subroutine get_wts(tstream,jyear,jj,wtl,wtr)
      implicit none
      type(timestream) :: tstream
      integer :: jyear,jj
      real*8 :: wtl,wtr
      do jj=1,tstream%nfileyrs
        if(tstream%fileyrs(jj).ge.jyear) exit
      enddo
      if(jj.eq.1) then
        wtl = 0d0
      elseif(jj.gt.tstream%nfileyrs) then
        wtl = 1d0
      else
        wtl = real(tstream%fileyrs(jj)-jyear,kind=8)/
     &        real(tstream%fileyrs(jj)-tstream%fileyrs(jj-1),kind=8)
      endif
      wtr = 1d0-wtl
      return
      end subroutine get_wts

      function have_prev_yr(tstream,jyear)
      implicit none
      logical :: have_prev_yr
      type(timestream) :: tstream
      integer :: jyear
      integer :: jj
      have_prev_yr = .false.
      do jj=1,tstream%nfileyrs
        if(tstream%fileyrs(jj).eq.jyear) then
          if(tstream%fileyrs(jj-1).eq.jyear-1) have_prev_yr=.true.
          exit
        endif
      enddo
      end function have_prev_yr

      function have_curr_yr(tstream,jyear)
      implicit none
      logical :: have_curr_yr
      type(timestream) :: tstream
      integer :: jyear
      integer :: jj
      have_curr_yr = .false.
      do jj=1,tstream%nfileyrs
        if(tstream%fileyrs(jj).eq.jyear) then
          have_curr_yr=.true.
          exit
        endif
      enddo
      end function have_curr_yr

      function have_next_yr(tstream,jyear)
      implicit none
      logical :: have_next_yr
      type(timestream) :: tstream
      integer :: jyear
      integer :: jj
      have_next_yr = .false.
      do jj=1,tstream%nfileyrs
        if(tstream%fileyrs(jj).eq.jyear) then
          if(tstream%fileyrs(jj+1).eq.jyear+1) have_next_yr=.true.
          exit
        endif
      enddo
      end function have_next_yr

      subroutine check_metadata(grid,tstream,jyear)
      use dd2d_utils, only : dist_grid
      use pario, only : par_open, par_close
     &     ,read_data, read_attr
     &     ,get_record_dimlen,get_dimlens,get_record_dimname
     &     ,variable_exists
!!       use param, only : sync_param
!#ifdef IN_MODELE
      use SystemTools, only : stLinkStatus,stFileList
!#endif
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: jyear
c
      integer :: fid,mon,lm
      character(len=200) :: fname,fname_eom,attname
      INTEGER :: J_0,J_1, I_0,I_1, M1,M2, npad
      INTEGER :: I,J,IDUM,RDIMLEN
      logical :: exists,monthly_data,daily_data,annual_data,pentad_data
c
      integer :: jyr,nfileyrs,jj,yr_ind,roff
      integer, dimension(:), allocatable :: fileyrs
      real*8, dimension(:), allocatable :: taxis
      real*8 :: wtl,wtr
      logical :: cyclic,multiple_yrs,multiple_files
      integer :: ndims,dlens(7)
c
      integer :: indx,yr0,mn0,pn0,dy0,tm0,yrx,rdimlen1
      integer :: linkstatus
      character(len=32) :: dname,tunits,tunits_off
      integer, parameter :: max_fname_len=128
      character(len=max_fname_len), allocatable :: flist(:)
      character(len=max_fname_len) :: thisline
      character(len=4) :: c4
      integer :: lsiter,nfiles,ifile,ios

c
      cyclic = tstream%cyclic

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      lm = tstream%lm

      nfileyrs = 0

!#ifdef IN_MODELE
      call stLinkStatus(tstream%fbase, linkstatus)
      if(linkstatus==2) then  ! this is a directory. todo: no hard-coded retcodes
        allocate(flist(1000)) ! 1000 files maximum
        call stFileList(tstream%fbase,flist,nfiles)

        ! determine available years from YYYY.nc file presence
        do lsiter=1,2
        if(lsiter.eq.2) allocate(fileyrs(nfileyrs))
        nfileyrs = 0
        do ifile=1,nfiles
          thisline = adjustl(flist(ifile))
          if(len_trim(thisline).ne.7) cycle
          if(thisline(5:7).ne.'.nc') cycle
          c4 = thisline(1:4)
          read(c4,*,iostat=ios) jyr
          if(ios.ne.0) cycle
          if(jyr.lt.0) cycle
          nfileyrs = nfileyrs + 1
          if(lsiter.eq.2) fileyrs(nfileyrs) = jyr
        enddo
        enddo

        if(nfileyrs.eq.0) then
          if(grid%am_i_globalroot) write(6,*)
     &         'read_stream: empty directory '//trim(tstream%fbase)
          call stop_model(
     &         'read_stream: empty input directory',255)
        endif

        ! Apparently there is no guarantee that stFileList will report
        ! YYYY.nc files in ascending year order, so sort post-hoc
        if(nfileyrs.gt.1) then
          call mergesort(nfileyrs,fileyrs)
        endif
        write(c4,'(i4)') fileyrs(1)
        tstream%firstfile = trim(tstream%fbase)//'/'//c4//'.nc'
      else
        tstream%firstfile = tstream%fbase
      endif ! fbase is directory or not
!#endif IN_MODELE

      multiple_yrs = nfileyrs.gt.0
      multiple_files = multiple_yrs
      if(multiple_yrs) then
          !tstream%cycl=0
          !cyclic=.false.
        allocate(tstream%fileyrs(0:nfileyrs+1))
        tstream%nfileyrs = nfileyrs
        tstream%fileyrs(1:nfileyrs) = fileyrs(1:nfileyrs)
        tstream%fileyrs(0) = fileyrs(1)
        tstream%fileyrs(nfileyrs+1) = fileyrs(nfileyrs)
        ! why not just files fileyrs(1) and skip get_wts call
        call get_wts(tstream,jyear,yr_ind,wtl,wtr)
        jyr = tstream%fileyrs(yr_ind)
        tstream%multiple_files = .true.
        call make_fname(tstream,jyr,fname,fname_eom)
      else
        daily_data = .false.
        annual_data = .false.
        pentad_data = .false.
        monthly_data = .true. ! default?
          !tstream%cycl = 1 ! may be overridden
          !cyclic = .true.  ! may be overridden
        fname = tstream%fbase
        fname_eom = trim(tstream%fbase)//'_eom'
          ! begin new
          ! check whether a time axis variable exists, and its units
        fid = par_open(grid,trim(fname),'read')
        rdimlen = get_record_dimlen(grid,fid)
        dname = 'notaname'
        call get_record_dimname(grid,fid,dname)
        if(variable_exists(grid,fid,dname)) then
          tunits = 'notaname'
          call read_attr(grid,fid,trim(dname),'units',idum,tunits)
          indx = index(tunits,' since ')
          if(indx.gt.1) then
            tunits = adjustl(tunits)
            indx = index(tunits,' since ')
            tunits_off = tunits(indx+7:len_trim(tunits))
            tunits = tunits(1:indx-1)
            do i=1,len_trim(tunits_off)
              if(tunits_off(i:i).eq.'-') tunits_off(i:i)=' '
              if(iachar(tunits_off(i:i)).eq.0) tunits_off(i:i)=' '
            enddo
            select case (trim(tunits))
            case ('years')
              read(tunits_off,'(i4)') yr0
              rdimlen1 = 1
            case ('months')
              read(tunits_off,'(i4,1x,i2)') yr0,mn0
              tm0 = mn0
              rdimlen1 = 12
            case ('pentads')
              read(tunits_off,'(i4,1x,i2)') yr0,pn0
              tm0 = pn0
              rdimlen1 = 73
            case ('days')
              read(tunits_off,'(i4,1x,i2,1x,i2)') yr0,mn0,dy0
              tm0 = dy0+jdendofm(mn0-1)
              rdimlen1 = 365
            case default
              call stop_model('unrecognized time axis units',255)
            end select
            multiple_yrs = rdimlen.gt.rdimlen1
            if(multiple_yrs) then
              nfileyrs = rdimlen/rdimlen1
              daily_data = rdimlen1 == 365
              annual_data = rdimlen1 == 1
              monthly_data = rdimlen1 == 12
              pentad_data = rdimlen1 == 73
              !deallocate(fileyrs)
              allocate(fileyrs(nfileyrs))
              allocate(taxis(rdimlen))
              call read_data(grid,fid,trim(dname),taxis,
     &             bcast_all=.true.)
              if(monthly_data .or. daily_data .or. pentad_data) then
                ! sanity-check time axis
                do jj=1,rdimlen-1
                  i = taxis(jj+1)-taxis(jj)
                  if(mod(i+rdimlen1-1,rdimlen1).ne.0)
     &                 call stop_model('read_stream: jump',255)
                enddo
                ! find roff = -1 + record corresponding to January of first year
                yrx = yr0 + (tm0-1+int(taxis(1)))/rdimlen1
                fileyrs(1) = yrx
                roff = 0
                do jj=2,rdimlen1
                  yrx = yr0 + (tm0-1+int(taxis(jj)))/rdimlen1
                  if(yrx .ne. fileyrs(1)) then
                    if(yrx.ne.fileyrs(1)+1) then
                      call stop_model('read_stream: jump',255)
                    endif
                    !fileyrs(1) = yrx
                    roff = jj-1 -rdimlen1
                    exit
                  endif
                enddo
                i = 1
                do jj=roff+1+rdimlen1,rdimlen,rdimlen1
                  i = i + 1
                  fileyrs(i) = yr0 + (tm0-1+int(taxis(jj)))/rdimlen1
                enddo
                tstream%roff = roff
              else
                fileyrs = yr0 + int(taxis(1:rdimlen))
              endif
              allocate(tstream%fileyrs(0:nfileyrs+1))
              tstream%nfileyrs = nfileyrs
              tstream%fileyrs(1:nfileyrs) = fileyrs(1:nfileyrs)
              tstream%fileyrs(0) = fileyrs(1)
              tstream%fileyrs(nfileyrs+1) = fileyrs(nfileyrs)
              !call get_wts(tstream,jyear,yr_ind,wtl,wtr)
              !jyr = tstream%fileyrs(yr_ind)
            endif
          endif
        endif
        call par_close(grid,fid)
        if(.not. multiple_yrs) then
          tstream%cycl = 1
          cyclic = .true.
        endif
          ! debug print begin
          !if(multiple_yrs) then
          !write(6,*) trim(tunits)
          !write(6,*) trim(tunits_off)
          !write(6,*) 'yr0 ',yr0
          !if(monthly_data .or. daily_data) write(6,*) 'tm0 ',tm0
          !if(daily_data) write(6,*) 'dy0 ',dy0
          !write(6,*) 'rdimlen ',rdimlen
          !write(6,*) 'fileyrs ',fileyrs
          !stop 'here'
          !endif
          ! debug print end
          ! end new
      endif
      !deallocate(fileyrs)
      fid = par_open(grid,trim(fname),'read')
      attname = trim(tstream%vname)//'name'
      call read_attr(grid,fid,'global',trim(attname),idum,
     &     tstream%vname)
      rdimlen = get_record_dimlen(grid,fid,
     &     checkvar=trim(tstream%vname))
      call get_dimlens(grid,fid,trim(tstream%vname),ndims,dlens)
      if(ndims.lt.3 .or. ndims.gt.5) then
        call stop_model('read_stream: bad dimension count'//
     &         'for variable '//trim(tstream%vname),255)
      endif
      call par_close(grid,fid)
      tstream%dlens(:) = 1
      if(ndims.ge.4) tstream%dlens(1) = dlens(3)
      if(ndims.ge.5) tstream%dlens(2) = dlens(4)
      lm = product(tstream%dlens)
      tstream%lm = lm
      tstream%ndims = ndims-1 ! not including time dimension

      if(multiple_files .or. .not.multiple_yrs) then
        if(  (rdimlen.gt.12 .and. rdimlen.lt.73) .or.
     &       (rdimlen.gt.73 .and. rdimlen.lt.365) ) then
          if(grid%am_i_globalroot) write(6,*)
     &         'read_netcdf: bad record dimension length'
          call stop_model('read_netcdf',255)
        endif
        daily_data = rdimlen == 365
        annual_data = rdimlen == 1
        monthly_data = rdimlen == 12
        pentad_data = rdimlen == 73
      endif

      !monthly_data = .not. (daily_data .or. annual_data)
      tstream%monthly_data = monthly_data
      tstream%daily_data = daily_data
      tstream%annual_data = annual_data
      tstream%pentad_data = pentad_data
      tstream%multiple_files = multiple_files

!!       call sync_param( trim(tstream%fbase)//'_cycl', tstream%cycl )
      tstream%cyclic = cyclic

      if(daily_data) then
        tstream%m2r = 365
        tstream%tinterp_method = nointerp
      elseif(annual_data) then
        tstream%m2r = 1
      elseif(pentad_data) then
        tstream%m2r = 73
      elseif(monthly_data) then
        tstream%m2r = 12
        inquire(file=trim(fname_eom), exist=tstream%eom_from_file)
        tstream%eom_from_file = tstream%eom_from_file .and.
     &       allocated(tstream%eom)
      endif
      if(tstream%tinterp_method.eq.ppm .and.
     &     (annual_data.or.pentad_data)) then
        tstream%tinterp_method = linm2m
        if(grid%am_i_globalroot) write(6,*)
     &       'disabling ppm and using linm2m for non-monthly file '//
     &       trim(tstream%fbase)
      endif
      if(tstream%tinterp_method.eq.ppm) then
        if(tstream%eom_from_file) then
          npad = 0
        else
          npad = 2
        endif
      elseif(tstream%tinterp_method.eq.linm2m) then
        npad = 1
      else
        npad = 0
      endif
      if(annual_data .and. cyclic) then
        npad = 0
        tstream%tinterp_method = nointerp
      endif
      tstream%npad = npad
      m1 = 1-npad; m2 = tstream%m2r+npad
      allocate(tstream%qty(I_0:I_1,J_0:J_1,LM,M1:M2))
      tstream%qty = 0.

      if(tstream%tinterp_method.eq.ppm) then
        allocate(tstream%eom(I_0:I_1,J_0:J_1,LM,0:12))
      endif

      end subroutine check_metadata

      subroutine mergesort(n,arr)
      implicit none
      integer :: n
      integer, dimension(n) :: arr
      integer, dimension(:), allocatable :: arr1
      integer :: i,ii,di,nn,ipass,npass,i1,i2,i1max,i2max
      allocate(arr1(n))
      di = 1
      npass = 0
      do while(di.lt.n)
        npass = npass + 1
        di = 2*di
      enddo
      di = 1
      do ipass=1,npass
        di = di*2
        nn = 0
        i1 = 1
        do i=1,n,di
          i2 = min(i1+di/2,n+1)
          i1max = i2-1
          i2max = min(i+di-1,n)
          do ii=1,min(di,n-i+1)
            nn = nn + 1
            if(i1.gt.i1max) then
              arr1(nn) = arr(i2)
              i2 = i2 + 1
            elseif(i2.gt.i2max) then
              arr1(nn) = arr(i1)
              i1 = i1 + 1
            elseif(arr(i1).lt.arr(i2)) then
              arr1(nn) = arr(i1)
              i1 = i1 + 1
            else
              arr1(nn) = arr(i2)
              i2 = i2 + 1
            endif
          enddo
          i1 = i2max + 1
        enddo
        arr(:) = arr1(:)
      enddo
      do i=1,n-1
        if(arr(i).gt.arr(i+1)) then
          stop 'bad mergesort'
        endif
      enddo
      end subroutine mergesort

      subroutine read_stream_netcdf(grid,tstream,jyear,jmon)
!@sum read_stream_netcdf reads one year of qty
      use dd2d_utils, only : dist_grid
      use pario, only : par_open, par_close, read_dist_data, read_data
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: jyear,jmon
c
      integer :: fid,mon,mon1,mon2,monoff,lm,ndims,m2r
      character(len=200) :: fnames(0:3),fnames_eom(0:3),fname_dum
      INTEGER :: J_0,J_1, I_0,I_1
      INTEGER :: I,J,K, npad
      logical :: firstcall,year_reset
     &     ,monthly_data,daily_data,annual_data,pentad_data,need_ends
     &     ,read_prev,read_prev_full,snglread,yearly_varying
c
      integer :: jyr,jj,yr_ind,roff
      real*8 :: wtl,wtr,wtlp,wtrp,wtln,wtrn
      logical :: cyclic,multiple_yrs,multiple_files,do_yr_interp
     &     ,do_yrm1_interp,do_yrp1_interp,continuous,early_in_year
      integer :: record0
      integer :: m1,m2,record1
      real*8, dimension(:,:,:,:,:), allocatable :: qty4d
c

      if(tstream%year_sv == jyear) return
      firstcall = tstream%year_sv == firstcall_int
      if(tstream%cyclic) then
        if(.not.firstcall .and. tstream%year_sv .ne. reset_int) return
      endif

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      fnames(:) = ''
      fnames_eom(:) = ''

      if(firstcall) then
        call check_metadata(grid,tstream,jyear)
      endif

      npad = tstream%npad
      cyclic = tstream%cyclic
      daily_data = tstream%daily_data
      annual_data = tstream%annual_data
      monthly_data = tstream%monthly_data
      pentad_data = tstream%pentad_data
      multiple_yrs = tstream%nfileyrs.gt.0
      multiple_files = tstream%multiple_files
      roff = tstream%roff

      need_ends = npad.gt.0
      yearly_varying = multiple_yrs .and. .not. cyclic
      continuous = yearly_varying .and. need_ends

      lm = tstream%lm
      ndims = tstream%ndims
      m2r = tstream%m2r

      ! check if jyear has gone either (1) backward or (2) forward by more than 1 year
      year_reset = firstcall
      if(.not.firstcall .and. tstream%nfileyrs.gt.0 .and.
     &     (tstream%year_sv .lt. jyear-1 .or. jyear.lt.tstream%year_sv))
     &     then
        !ny: call get_wts(tstream,tstream%year_sv,jj,    wtl,wtr)
        !ny: call get_wts(tstream,jyear,          yr_ind,wtl,wtr)
        ! if interp interval has changed, reset necessary elements of tstream
        year_reset = .true. !ny: yr_ind.lt.jj .or. yr_ind.gt.jj+1
        if(year_reset) then
          tstream%year_sv2 = -9999 ! reset_int?
        endif
      endif

      tstream%year_sv = jyear

      if(continuous) then
        ! interannually varying conditions: time interp for the
        ! beginning of this year needs data from the end of the prev year
        do i=1-npad,0
          tstream%qty(:,:,:,i) = tstream%qty(:,:,:,m2r+i)
        enddo
      endif

      if(multiple_yrs) then
        call get_wts(tstream,jyear,yr_ind,wtl,wtr)
      else
        yr_ind = 1 ! not used
      endif

      if(monthly_data) then
        early_in_year = jmon.lt.3
      elseif(pentad_data) then
        early_in_year = jmon.lt.2
      elseif(annual_data) then
        early_in_year = jmon.lt.8 ! should use day instead
      else
        early_in_year = .false.
      endif

      ! determine whether to read data from a prior year
      do_yrm1_interp = .false.
      if(year_reset .and. multiple_yrs) then
        read_prev = .false.
        if(.not.have_curr_yr(tstream,jyear)) then
          read_prev = .true.
          read_prev_full = .true.
          do_yrm1_interp = continuous .and. early_in_year
        elseif(need_ends .and. early_in_year) then
          if(have_prev_yr(tstream,jyear)) then
            if(yearly_varying) then
              read_prev = .true.
              read_prev_full = .false.
            endif
          else
            read_prev = .true.
            read_prev_full = .true. ! necessary s.t. qty1 set?
            if(jyear.eq.tstream%fileyrs(1)) then
              ! we can only get here having part of the immediately prior year
              read_prev_full = .false.
            endif
            do_yrm1_interp = yearly_varying .and. read_prev_full
          endif
        endif
        if(read_prev) then
          if(read_prev_full) then
            k = 1
          else
            k = 0
          endif
          jyr = tstream%fileyrs(yr_ind-1)
          call make_fname(tstream,jyr,fnames(k),fnames_eom(k))
        endif
      endif

      if(multiple_yrs) then
        if(jyear.gt.tstream%year_sv2 .and.
     &     jyear.le.tstream%fileyrs(tstream%nfileyrs)) then
          tstream%year_sv2 = tstream%fileyrs(yr_ind)
          call make_fname(tstream,tstream%fileyrs(yr_ind),fnames(2),
     &         fnames_eom(2))
        else
          fnames(2) = ''
          fnames_eom(2) = ''
        endif
      else
        fnames(2) = tstream%fbase
        fnames_eom(2) = trim(tstream%fbase)//'_eom'
      endif

      do_yr_interp = multiple_yrs .and. jyear.lt.tstream%year_sv2

      !do_yrp1_interp = do_yr_interp .and. continuous
      do_yrp1_interp = continuous .and. (do_yr_interp .or. year_reset)

      if(continuous .and. jyear.eq.tstream%year_sv2) then
        ! flag that we need to read first part of next available year
        call make_fname(tstream,tstream%fileyrs(yr_ind+1),fnames(3))
        do_yrp1_interp = .not.have_next_yr(tstream,jyear)
      endif

      if(multiple_yrs .and. jyear.eq.tstream%year_sv2 .and.
     &     len_trim(fnames(2)).eq.0) then
        ! only happens when reaching the endpoint of an interp interval
        tstream%qty(:,:,:,1:m2r) = tstream%qty2
      endif

      !
      ! read data
      ! k = 0     last npad periods of year-1
      !     1     entirety of previous available year
      !     2     entirety of this or next available year
      !     3     first npad periods of next available year
      ! k=0,1     only execute on startup/reset

      do k=0,3
        if(len_trim(fnames(k)).eq.0) cycle

        mon1=1; mon2=m2r
        if(k.eq.0) mon1=m2r-npad+1
        if(k.eq.3) mon2=npad
        monoff = 0
        if(k.eq.0) monoff = -m2r
        if(k.eq.3) monoff = +m2r

        jj = yr_ind
        if(k.le.1) jj = max(1,jj-1)
        if(k.eq.3) jj = min(jj+1,tstream%nfileyrs)
        record0 = m2r*(jj-1)+roff

        if(record0.lt.0) then ! first year has only partial data
          if(k.ne.2) call stop_model(
     &         'read_stream_netcdf: should not happen',255)
          mon1 = mon1 - record0
        endif

        snglread = multiple_files .or. .not. multiple_yrs
        if(snglread) record0 = 0
        snglread = snglread .and. (k.eq.1 .or. k.eq.2)

        if(snglread) then
          m1 = 1; m2 = m2r
          record1 = 1
        else
          m1 = mon1+monoff; m2 = mon2+monoff
          record1 = mon1+record0
        endif

        fid = par_open(grid,trim(fnames(k)),'read')
        if(ndims.eq.2) then
          call read_dist_data(grid,fid,trim(tstream%vname),
     &         tstream%qty(:,:,1,m1:m2),record1=record1)
        elseif(ndims.eq.3) then
          call read_dist_data(grid,fid,trim(tstream%vname),
     &         tstream%qty(:,:,:,m1:m2),record1=record1)
        else ! infrequent case - recreate extra dims for read_dist_data
          allocate(qty4d(i_0:i_1,j_0:j_1,
     &         tstream%dlens(1),tstream%dlens(2),m1:m2))
          call read_dist_data(grid,fid,trim(tstream%vname),
     &         qty4d,record1=record1)
          tstream%qty(:,:,:,m1:m2) =
     &         reshape(qty4d,shape(tstream%qty(:,:,:,m1:m2)))
          deallocate(qty4d)
        endif
        call par_close(grid,fid)

        if(k.eq.1) then ! startup
          call check_alloc(tstream)
          tstream%qty1=tstream%qty(:,:,:,1:m2r)
          tstream%qty2=tstream%qty(:,:,:,1:m2r)
        endif
        if(k.eq.2 .and. do_yr_interp) then
          tstream%qty1 = tstream%qty2
          tstream%qty2 = tstream%qty(:,:,:,1:m2r)
        endif
        if(k.eq.3 .and. do_yrp1_interp) then
          call check_alloc(tstream)
          do i=1,npad
!deferred copy  tstream%qty1(:,:,:,i) = tstream%qty(:,:,:,i)
            tstream%qty2(:,:,:,i) = tstream%qty(:,:,:,m2r+i)
          enddo
        endif
      enddo

      ! do year interp

      if(do_yrm1_interp) then
        if(have_curr_yr(tstream,jyear-1)) then
          wtln = 1d0; wtrn = 0d0
        else
          call get_wts(tstream,jyear-1,jj,wtln,wtrn)
        endif
        do i=1-npad,0
          tstream%qty(:,:,:,i) =
     &         wtln*tstream%qty1(:,:,:,m2r+i)
     &        +wtrn*tstream%qty(:,:,:,m2r+i)
        enddo
      endif

      ! deferred copy
      if(do_yrp1_interp .and. len_trim(fnames(3)).ne.0) then
        do i=1,npad
          tstream%qty1(:,:,:,i) = tstream%qty(:,:,:,i)
        enddo
      endif

      if(do_yr_interp) then
        tstream%qty(:,:,:,1:m2r) = wtl*tstream%qty1 +wtr*tstream%qty2
      endif

      if(do_yrp1_interp) then
        call get_wts(tstream,jyear+1,jj,wtlp,wtrp)
        do i=1,npad
          tstream%qty(:,:,:,m2r+i) =
     &         wtlp*tstream%qty1(:,:,:,i)
     &        +wtrp*tstream%qty2(:,:,:,i)
        enddo
      endif

      if(need_ends .and. cyclic .and. year_reset) then
        do i=1-npad,0
          tstream%qty(:,:,:,i) = tstream%qty(:,:,:,m2r+i)
        enddo
        do i=1,npad
          tstream%qty(:,:,:,m2r+i) = tstream%qty(:,:,:,i)
        enddo
      endif

      ! prep for next call
      if(jyear.eq.tstream%year_sv2 .and. yearly_varying .and.
     &   .not.have_next_yr(tstream,jyear) ) then
        call check_alloc(tstream)
        tstream%qty1 = tstream%qty(:,:,:,1:m2r)
        tstream%qty2 = tstream%qty(:,:,:,1:m2r)
      endif

      ! read or create EOM values for PPM
      if(tstream%tinterp_method.eq.ppm .and. monthly_data) then

      if(tstream%eom_from_file) then

        if(yearly_varying) then
        ! interannually varying conditions: time interp for the
        ! beginning of this year needs data from the end of the prev year
          tstream%eom(:,:,:,0) = tstream%eom(:,:,:,12)
        endif

        if(year_reset .and. yearly_varying .and. jmon.eq.1) then
          jyr = tstream%fileyrs(yr_ind-1)
          call make_fname(tstream,jyr,fname_dum,fnames_eom(0))
        endif

        do k=0,2,2
          if(len_trim(fnames_eom(k)).eq.0) cycle

          mon1=1; mon2=12
          if(k.eq.0) mon1=12
          monoff = 0
          if(k.eq.0) monoff = -12

          jj = yr_ind
          if(k.eq.0) jj = max(1,jj-1)
          record0 = 12*(jj-1)+roff

          if(record0.lt.0) then ! first year has only partial data
            if(k.ne.2) call stop_model(
     &           'read_stream_netcdf: should not happen',255)
            mon1 = mon1 - record0
          endif
          
          snglread = multiple_files .or. .not. multiple_yrs
          if(snglread) record0 = 0
          snglread = snglread .and. (k.eq.1 .or. k.eq.2)

          fid = par_open(grid,trim(fnames_eom(k)),'read')
          if(snglread) then
            if(ndims.eq.2) then
              call read_dist_data(grid,fid,trim(tstream%vname)//'_eom',
     &             tstream%eom(:,:,1,1:12))
            elseif(ndims.eq.3) then
              call read_dist_data(grid,fid,trim(tstream%vname)//'_eom',
     &             tstream%eom(:,:,:,1:12))
            else
              call stop_model('implement this case',255)
            endif
          else
            do mon=mon1,mon2
              if(ndims.eq.2) then
               call read_dist_data(grid,fid,trim(tstream%vname)//'_eom',
     &               tstream%eom(:,:,1,mon+monoff),record=mon+record0)
              elseif(ndims.eq.3) then
               call read_dist_data(grid,fid,trim(tstream%vname)//'_eom',
     &               tstream%eom(:,:,:,mon+monoff),record=mon+record0)
              else
                call stop_model('implement this case',255)
              endif
            enddo
          endif
          call par_close(grid,fid)
        enddo

        if(cyclic .and. firstcall)
     &       tstream%eom(:,:,:,0) = tstream%eom(:,:,:,12)

      else ! calculate EOM values
        if(firstcall .or. .not.cyclic) then
          do k=1,lm
            call edginterp(tstream%qty(:,:,k,:),tstream%eom(:,:,k,:),
     &           tstream%qmin,tstream%qmax)
          enddo
          if(trim(tstream%fbase).eq.'OSST') then ! hack checking file name
            do mon=0,12
            do j=j_0,j_1
            do i=i_0,i_1
              if(tstream%eom(i,j,1,mon).lt.-1.8d0)
     &             tstream%eom(i,j,1,mon)=-1.8d0
              if(tstream%qty(i,j,1,mon).le.-1.79d0)
     &             tstream%eom(i,j,1,mon)=-1.8d0
              if(tstream%qty(i,j,1,mon+1).lt.-1.79d0)
     &             tstream%eom(i,j,1,mon)=-1.8d0
            enddo
            enddo
            enddo
          endif
        endif
      endif ! reading versuing calculating eom

      endif ! ppm, monthly

      return
      contains
      subroutine edginterp(x,xe,xmin,xmax)
      real*8 :: x(i_0:i_1,j_0:j_1,-1:14)
      real*8 :: xe(i_0:i_1,j_0:j_1,0:12)
      real*8 :: xmin,xmax
      real*8, parameter :: by12=1d0/12d0,sevby12=7d0*by12
      integer :: i,j,me
      do me=0,12
      do j=j_0,j_1
      do i=i_0,i_1
        xe(i,j,me) = -by12*(x(i,j,me-1)+x(i,j,me+2))
     &            +sevby12*(x(i,j,me  )+x(i,j,me+1))
        xe(i,j,me) = min(max(xmin,xe(i,j,me)),xmax)
        if(minval(x(i,j,me:me+1)).le.xmin) xe(i,j,me)=xmin
        if(maxval(x(i,j,me:me+1)).ge.xmax) xe(i,j,me)=xmax
      enddo
      enddo
      enddo
      end subroutine edginterp
      subroutine check_alloc(tstream)
      type(timestream) :: tstream
      if(.not.allocated(tstream%qty1)) then
        allocate(
     &       tstream%qty1(I_0:I_1,J_0:J_1,tstream%LM,tstream%m2r),
     &       tstream%qty2(I_0:I_1,J_0:J_1,tstream%LM,tstream%m2r)
     &       )
        tstream%qty1 = 0.
        tstream%qty2 = 0.
      endif
      end subroutine check_alloc

      end subroutine read_stream_netcdf

      subroutine make_fname(tstream,jyr,fname,fname_eom)
      type(timestream) :: tstream
      character(len=200) :: fname
      character(len=200), optional :: fname_eom
      integer :: jyr
      character(len=4) :: year_string
      if(tstream%multiple_files) then
        write(year_string,'(i4)') jyr
        fname = trim(tstream%fbase)//'/'//year_string//'.nc'
        if(present(fname_eom)) fname_eom =
     &       trim(tstream%fbase)//'_eom/'//year_string//'.nc'
      else
        fname = tstream%fbase
        if(present(fname_eom)) fname_eom = trim(fname)//'_eom'
      endif
      end subroutine make_fname

      subroutine read_stream_2d(grid,tstream,jyear,jday,arr,tlim)
!@sum read_stream_2d a wrapper to call read_stream_3d for 2d outputs
!@+   by adding an extra dimension of size 1.
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: jyear,jday
      real*8, dimension(grid%i_strt_halo:,grid%j_strt_halo:) :: arr
      real*8, dimension(2,grid%i_strt_halo:grid%i_stop_halo,
     &                    grid%j_strt_halo:grid%j_stop_halo),
     &     optional :: tlim
      real*8, dimension(:,:,:), allocatable :: arr3d
      integer :: i,j
      allocate(arr3d(grid%i_strt_halo:grid%i_stop_halo,
     &               grid%j_strt_halo:grid%j_stop_halo,1))
      if(present(tlim)) then
        call read_stream_3d(grid,tstream,jyear,jday,arr3d,tlim)
      else
        call read_stream_3d(grid,tstream,jyear,jday,arr3d)
      endif
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,grid%i_stop
        if(tstream%msk(i,j).eq.0d0) cycle
        arr(i,j) = arr3d(i,j,1)
      enddo
      enddo
      deallocate(arr3d)
      end subroutine read_stream_2d

      subroutine read_stream_4d(grid,tstream,jyear,jday,arr)!,tlim)
!@sum read_stream_4d a wrapper to call read_stream_3d for 4d outputs
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: jyear,jday
      real*8, dimension(grid%i_strt_halo:,grid%j_strt_halo:,:,:) :: arr
!      real*8, dimension(2,grid%i_strt_halo:grid%i_stop_halo,
!     &                    grid%j_strt_halo:grid%j_stop_halo,:,:),
!     &     optional :: tlim
      real*8, dimension(:,:,:), allocatable :: arr3d
      integer :: k,l,kl
      integer :: j_0,j_1, i_0,i_1
      allocate(arr3d(grid%i_strt_halo:grid%i_stop_halo,
     &               grid%j_strt_halo:grid%j_stop_halo,tstream%lm))
!      if(present(tlim)) then
!        call read_stream_3d(grid,tstream,jyear,jday,arr3d,tlim)
!      else
        call read_stream_3d(grid,tstream,jyear,jday,arr3d)
!      endif
      !arr = reshape(arr3d,shape(arr))
      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop
      kl = 0
      do l=1,tstream%dlens(2)
      do k=1,tstream%dlens(1)
        kl = kl + 1
        arr(i_0:i_1,j_0:j_1,k,l) = arr3d(i_0:i_1,j_0:j_1,kl)
      enddo
      enddo
      deallocate(arr3d)
      end subroutine read_stream_4d

      subroutine read_stream_3d(grid,tstream,jyear,jday,arr,tlim)
!@sum read_stream_3d implementation of read_stream interface for 3D outputs.
!@+   See documentation at beginning of module for the read_stream interface.
!@+   The reads from disk and inter-annual time interpolations are performed
!@+   by read_year. This routine performs the intra-annual time interpolation.
      use dd2d_utils, only : dist_grid
      implicit none
      type(dist_grid) :: grid
      type(timestream) :: tstream
      integer :: jyear,jday
      real*8, dimension(grid%i_strt_halo:,
     &                  grid%j_strt_halo:,:) :: arr
!@var tlim optional argument used for PPM interpolation with limits.
!@+   tlim(1) = t0  tlim(2) = t1
!@+   see ppm_frac,ppm_tlim code blocks below for definition of t0 and t1
      real*8, dimension(2,grid%i_strt_halo:grid%i_stop_halo,
     &                    grid%j_strt_halo:grid%j_stop_halo),
     &     optional :: tlim
      real*8, parameter :: by12=1d0/12d0, teeny=1d-30
      integer i,j,imon,jmon,jdate,tinterp_method
      real*8 time,frac
      integer :: j_0,j_1, i_0,i_1

      real*8 :: a,b,c,e0,e1,csq,t0,t1
      integer, parameter :: ppm_frac = -99*ppm,ppm_tlim=-999*ppm,
     &     ppm_nonneg = -9999*ppm
      integer :: l,lm

      arr = 0.

      lm = tstream%lm

      jmon=1
      do while (jday.gt.jdendofm(jmon))
        jmon=jmon+1
      end do
      jdate=jday-jdendofm(jmon-1)

      call read_year(grid,tstream,jyear,jmon)

      i_0 = grid%i_strt
      i_1 = grid%i_stop
      j_0 = grid%j_strt
      j_1 = grid%j_stop

      tinterp_method = tstream%tinterp_method

      if(tinterp_method.eq.ppm) then
        ! This is not yet coded to work for input data other than monthly
        time=(jdate-.5)/(jdendofm(jmon)-jdendofm(jmon-1))-.5 ! -.5<time<.5
        if(trim(tstream%vname).eq.'ZSI') ! temp hack checking var name
     &       tinterp_method = ppm_tlim
        if(tstream%qmin.eq.0d0) then
          if(tstream%qmax.eq.1d0) then
            tinterp_method = ppm_frac
          elseif(tstream%qmax.ge.1d20) then
            tinterp_method = ppm_nonneg
          else
            call stop_model('unrecognized qmax for ppm',255)
          endif
        endif
      endif

      select case (tinterp_method)

      case(nointerp)
! no time interpolation required
        if(tstream%daily_data) then
          imon = jday
        elseif(tstream%annual_data) then
          imon = 1
        elseif(tstream%pentad_data) then
          imon = 1+(jday-1)/5
        elseif(tstream%monthly_data) then
          imon = jmon
        endif

        do l=1,lm
          arr(i_0:i_1,j_0:j_1,l) = tstream%qty(i_0:i_1,j_0:j_1,l,imon)
        enddo

      case(linm2m)
! linear interpolation between period midpoints

        if(tstream%daily_data) call stop_model
     &       ('linm2m currently not coded for daily data',255)

        if(tstream%monthly_data) then
          if(jday.le.jdmidofm(jmon)) then
            imon = jmon
          else
            imon = jmon + 1
          endif
          frac = real(jdmidofm(imon)-jday,kind=8)/
     &               (jdmidofm(imon)-jdmidofm(imon-1))
        elseif(tstream%annual_data) then
          if(jday.le.183) then
            imon = 1
            frac = real(183-jday,kind=8)/365d0
          else
            imon = 2
            frac = real(548-jday,kind=8)/365d0
          endif
        elseif(tstream%pentad_data) then
          imon = 1+(jday+1)/5
          frac = real(3+(imon-1)*5-jday,kind=8)/5d0
        endif
        do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          if(tstream%msk(i,j).eq.0d0) cycle
          arr(i,j,l) =
     &        frac*tstream%qty(i,j,l,imon-1)
     &  +(1.-frac)*tstream%qty(i,j,l,imon)
        enddo
        enddo
        enddo

      case(ppm)

! Piecewise Parabolic Method.  A parabola is constructed for each period
! whose mean, initial, and final values are a, e0, and e1.
        do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          if(tstream%msk(i,j).eq.0d0) cycle
          a = tstream%qty(i,j,l,jmon)    ! period mean
          e0 = tstream%eom(i,j,l,jmon-1) ! value at beginning of period
          e1 = tstream%eom(i,j,l,jmon)   ! value at end of period
          b = e1-e0                      ! mean of first time derivative
          c = 3.*(e1+e0) - 6.*a          ! second time derivative (curvature)
          arr(i,j,l) = a+b*time+c*(time**2-by12)
        enddo
        enddo
        enddo

      case(ppm_nonneg)

! Piecewise Parabolic Method limiting the interpolant to be non-negative.
! See previous code block for pure-PPM details.
! If the minimum of the unadjusted parabola is less than zero, the
! form of the fit is instead taken as piecewise linear over three time
! intervals: (1) time<t0 (2) t0<time<t1 (3) time>t1
! The constant value in the second interval is 0 if the parabola undershot.
! The values of t0 and t1 can be determined from the additional criteria
! that the mean over the three intervals is a, the initial value of the
! first interval is e0, and the final value of the last interval is e1.
! The length of the first and/or last intervals may be zero.
!
        do l=1,lm
        do j=j_0,j_1
        do i=i_0,i_1
          if(tstream%msk(i,j).eq.0d0) cycle
          a = tstream%qty(i,j,l,jmon)
          if(a.le.0.) then
            !call stop_model('read_stream: bad monthly mean',255)
            arr(i,j,l) = a ! keep constant value
            cycle               
          endif
          e0 = tstream%eom(i,j,l,jmon-1)
          e1 = tstream%eom(i,j,l,jmon)
          b=e1-e0
          c=3.*(e1+e0) - 6.*a
          arr(i,j,l)=a+b*time+c*(time**2-by12) ! default: pure PPM
          if(abs(c) .gt. abs(b)) then ! but check if linear fit is needed
            csq=c*(a*c - .25*b**2 - c**2*by12)
            if(csq.lt.0.) then        ! quadratic fit at apex < 0
              b = .5*(e0**2 + e1**2) / a
              if(e0-b*(time+.5) .gt. 0.)  then
                arr(i,j,l) = e0 - b*(time+.5) !  time < t0
              elseif(e1-b*(.5-time) .gt. 0.)  then
                arr(i,j,l) = e1 - b*(.5-time) !  t1 < time
              else
                arr(i,j,l) = 0.               !  t0 < time < t1
              end if
            end if
          end if
          !if(arr(i,j,l).lt.0.) then
          !  write(6,*) 'negative output',i,j,l,jday,jmon,arr(i,j,l)
          !endif
        enddo
        enddo
        enddo

      case(ppm_frac)

! Piecewise Parabolic Method limiting the interpolant to the range 0-1.
! See previous code block for pure-PPM details.
! If the min/max of the unadjusted parabola is outside the 0-1 range, the
! form of the fit is instead taken as piecewise linear over three time
! intervals: (1) time<t0 (2) t0<time<t1 (3) time>t1
! The constant value in the second interval is 1 if the parabola overshot,
! and 0 if it undershot.
! The values of t0 and t1 can be determined from the additional criteria
! that the mean over the three intervals is a, the initial value of the
! first interval is e0, and the final value of the last interval is e1.
! The length of the first and/or last intervals may be zero.
!
        do j=j_0,j_1
        do i=i_0,i_1
          if(tstream%msk(i,j).eq.0d0) cycle
          a = tstream%qty(i,j,1,jmon)
          if(present(tlim)) tlim(:,i,j) = (/ -1d30, 1d30 /)
          if(a.le.0. .or. a.ge.1.) then
            !call stop_model('read_stream: bad monthly mean',255)
            arr(i,j,1) = a ! keep constant value
            cycle               
          endif
          e0 = tstream%eom(i,j,1,jmon-1)
          e1 = tstream%eom(i,j,1,jmon)
          b=e1-e0
          c=3.*(e1+e0) - 6.*a
          arr(i,j,1)=a+b*time+c*(time**2-by12) ! default: pure PPM
          if(abs(c) .gt. abs(b)) then ! but check if linear fit is needed
            csq=c*(a*c - .25*b**2 - c**2*by12)
            if(csq.lt.0.) then        ! quadratic fit at apex < 0
              b = .5*(e0**2 + e1**2) / a
              if(present(tlim)) then
                tlim(1,i,j) = e0/b - .5d0
                tlim(2,i,j) = .5d0 - e1/b
              endif
              if(e0-b*(time+.5) .gt. 0.)  then
                arr(i,j,1) = e0 - b*(time+.5) !  time < t0
              elseif(e1-b*(.5-time) .gt. 0.)  then
                arr(i,j,1) = e1 - b*(.5-time) !  t1 < time
              else
                arr(i,j,1) = 0.               !  t0 < time < t1
              end if
            elseif(csq.gt.c**2)  then ! quadratic fit at apex > 1
              b = .5*((e0-1.)**2 + (e1-1.)**2) / (a-1.)
              if(e0-b*(time+.5) .lt. 1.)  then
                arr(i,j,1) = e0 - b*(time+.5) !  time < t0
              elseif(e1-b*(.5-time) .lt. 1.)  then
                arr(i,j,1) = e1 - b*(.5-time) !  t1 < time
              else
                arr(i,j,1) = 1.               !  t0 < time < t1
              end if
            end if
          end if
        enddo
        enddo

      case(ppm_tlim)

! Piecewise Parabolic Method for a quantity constrained to be
! zero over the time interval t0<time<t1, where t0 and t1 were
! computed during a previous call to this routine using the
! ppm_frac code block listed above. If the ppm_frac call kept
! its unadjusted parabola, the unadjusted parabola is used
! here also (t0 and t1 both lie outside the interval [-1:1]).
! Otherwise, the fit is piecewise linear over three intervals
! as explained for the ppm_frac case.

        do j=j_0,j_1
        do i=i_0,i_1
          if(tstream%msk(i,j).eq.0d0) cycle
          e0 = tstream%eom(i,j,1,jmon-1)
          e1 = tstream%eom(i,j,1,jmon)
          t0 = tlim(1,i,j)
          t1 = tlim(2,i,j)
          arr(i,j,1) = 0.
          if(t0.lt.-1. .and. t1.gt.1.) then ! unadjusted parabola
            a = tstream%qty(i,j,1,jmon)
            b = e1-e0
            c = 3.*(e1+e0) - 6.*a
            arr(i,j,1)=a+b*time+c*(time**2-by12)
          elseif(time.lt.t0) then
            arr(i,j,1)=e0*(1d0-(time+.5d0)/(t0+.5d0+teeny))
          elseif(time.gt.t1) then
            arr(i,j,1)=e1*(1d0-(.5d0-time)/(.5d0-t1+teeny))
          else
            arr(i,j,1)=0.
          endif
        enddo
        enddo

      end select ! tinterp_method

      return
      end subroutine read_stream_3d

      end module timestream_mod

!!       subroutine read_stream_gissfmt(grid,tstream,jyear,jmon)
!!       use model_com, only : iyear1,ItimeE,Nday
!!       use dd2d_utils, only : dist_grid
!!       use pario_fbsa, only :
!!      &    READ_PARALLEL,BACKSPACE_PARALLEL,MREAD_PARALLEL,READT_PARALLEL
!!      &   ,SKIP_PARALLEL
!!       use filemanager, only : openunit,closeunit
!!       use param, only : sync_param
!!       implicit none
!!       type(dist_grid) :: grid
!!       type(timestream) :: tstream
!!       integer :: jyear,jmon
!! c
!!       integer :: mon,lstmon,nmon_read,iskip
!!       INTEGER :: I_0H,I_1H,J_0H,J_1H
!!       INTEGER :: I,J,M,M1
!! !@var TEMP_LOCAL stores AOST+EOST1 or ARSI+ERST1
!!       REAL*8 :: TEMP_LOCAL(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
!!      &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO,2)
!!       integer :: end_year,end_month,end_day,end_date,end_hour
!!       character(len=4) :: c4
!!       logical :: firstcall
!! 
!!       if(tstream%year_sv == jyear) return
!!       firstcall = tstream%year_sv < 0
!!       if(tstream%cycl.eq.1 .and. .not.firstcall) return
!!       tstream%year_sv = jyear
!! 
!!       ! For convenience:
!!       ! Get the end month/year info for the simulation.
!!       ! It is specified in the rundeck, but only ItimeE
!!       ! is globally visible.
!!       call getdte(ItimeE,Nday,IYear1,
!!      &       end_year,end_month,end_day,
!!      &       end_date,end_hour,c4)
!! 
!!       i_0h = grid%i_strt_halo
!!       i_1h = grid%i_stop_halo
!!       j_0h = grid%j_strt_halo
!!       j_1h = grid%j_stop_halo
!! 
!!       ! determine the number of months of this year that need to be read
!!       if(jyear.lt.end_year) then
!!         nmon_read = 12
!!       else
!!         nmon_read = end_month
!!       endif
!! 
!!       if(firstcall) then ! first call
!!         nmon_read = nmon_read - jmon + 1
!!         allocate(tstream%qty(I_0H:I_1H,J_0H:J_1H,1,1:12))
!!         if(tstream%fmt.eq.gissclim_fmt) then
!!           tstream%cycl = 1
!!         else
!!           tstream%cycl = 0
!!         endif
!! 
!!         call openunit(trim(tstream%fbase),tstream%fid,.true.,.true.)
!! 
!! C**** Skip any leading records
!!         do iskip=1,tstream%nskip
!!           call skip_parallel(tstream%fid)
!!         enddo
!! 
!!         if(tstream%cycl.eq.1) then
!!           ! read climatology
!!           do mon=1,12
!!             CALL READT_PARALLEL(grid,tstream%fid,trim(tstream%fbase),
!!      &           TEMP_LOCAL,1)
!!             tstream%qty(:,:,1,mon)= TEMP_LOCAL(:,:,1)
!!             tstream%eom(:,:,mon)  = TEMP_LOCAL(:,:,2)
!!           enddo
!!           tstream%eom(:,:,0) = tstream%eom(:,:,12)
!!           call closeunit(tstream%fid)
!!         elseif(tstream%cycl.eq.0) then
!!           if (grid%am_i_globalroot) then
!!             write(6,*) '********************************************'
!!             write(6,*) '* Make sure that IYEAR1 is consistent with *'
!!             write(6,*) '*    the data file '//trim(tstream%fbase)
!!             write(6,*) '********************************************'
!!             write(6,*) 'IYEAR1=',IYEAR1
!!           end if
!!           ! advance to the previous month
!!           LSTMON=JMON-1+(JYEAR-IYEAR1)*JMperY
!!           m = lstmon-1
!!           do while(m.lt.lstmon)
!!             call READ_PARALLEL(grid, M, tstream%fid)
!!           enddo
!!           CALL BACKSPACE_PARALLEL( tstream%fid )
!!           ! read the EOM from the previous month
!!           CALL MREAD_PARALLEL(GRID,tstream%fid,trim(tstream%fbase), m,
!!      &         TEMP_LOCAL)
!!           if(jmon.eq.1) then
!!             tstream%eom(:,:,12)  = TEMP_LOCAL(:,:,2)
!!           else
!!             tstream%eom(:,:,jmon-1)  = TEMP_LOCAL(:,:,2)
!!           endif
!!         endif
!! 
!!         call sync_param( trim(tstream%fbase)//'_cycl', tstream%cycl )
!! 
!!       endif
!! 
!!       if(tstream%cycl.eq.0) then
!!         ! time interp for the beginning of next year needs data from the end of this year
!!         tstream%eom(:,:,0) = tstream%eom(:,:,12)
!!         ! read the next year of data
!!         do mon=jmon,jmon+nmon_read-1
!!           CALL MREAD_PARALLEL(GRID,tstream%fid,trim(tstream%fbase),M,
!!      &         TEMP_LOCAL)
!!           tstream%qty(:,:,1,mon)= TEMP_LOCAL(:,:,1)
!!           tstream%eom(:,:,mon)  = TEMP_LOCAL(:,:,2)
!!         enddo
!!       endif
!! 
!!       return
!!       end subroutine read_stream_gissfmt
