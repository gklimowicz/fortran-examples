!@sum  POUT_netcdf default output routines for netcdf formats
!@auth M. Kelley

C****
C**** If other formats are desired, please replace these routines
C**** with ones appropriate for your chosen output format, but with the
C**** same interface
C****
C**** Note: it would be nice to amalgamate IL and JL, but that will
C**** have to wait.

      module ncout
!@sum  ncout handles the writing of output fields in netcdf-format
!@auth M. Kelley
      use MODEL_COM, only : xlabel,lrunid
      use MDIAG_COM, only : acc_period,
     &     sname_strlen,units_strlen,lname_strlen
      implicit none

      include 'netcdf.inc'

      private
      public ::
     &     outfile,jlfile,ilfile
     &    ,open_out,def_dim_out,set_dim_out,close_out
     &    ,status_out,varid_out,out_fid
     &    ,ndims_out,dimids_out,dimlens_out,var_name
     &    ,units,long_name,missing,missing8,real_att_name,real_att
     &    ,disk_dtype,prog_dtype,lat_dg
     &    ,iu_ij,im,jm,lm,lm_req,im_data,iu_ijl
     &    ,iu_ijk,iu_il,iu_j,iu_jc,iu_jl,iu_isccp,iu_diurn,iu_hdiurn
     &    ,iu_wp,def_missing,srt,cnt,write_whole_array
     &    ,sname_strlen,units_strlen,lname_strlen

!@var iu_ij,iu_jl,iu_il,iu_j !  units for selected diag. output
      integer iu_ij,iu_ijk,iu_il,iu_j,iu_jl,iu_isccp,iu_diurn,iu_hdiurn
     &     ,iu_jc,iu_wp,iu_ijl
!@var im,jm,lm,lm_req local dimensions set in open_* routines
!@var im_data inner dimension of arrays passed to pout_*
!     it will differ from im because of array padding to hold means, etc.
      integer :: im,jm,lm,lm_req,im_data
!@var JMMAX maximum conceivable JM
      INTEGER, PARAMETER :: JMMAX=200
!@var LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JMMAX,2) :: LAT_DG

      character(len=80) :: outfile,jlfile,ilfile
!@var status_out return code from netcdf calls
!@var out_fid ID-number of output file
      integer :: status_out,out_fid
!@var ndims_out number of dimensions of current output field
      integer :: ndims_out
!@var varid_out variable ID-number of current output field
      integer :: varid_out
!@var dimids_out dimension ID-numbers (1:ndims_out) of output field
!@var dimlens_out dimension lengths (1:ndims_out) of output field
      integer, dimension(7) :: dimids_out,dimlens_out
!@var srt start indices for output from current output field
!@var cnt number of elements to output from current output field
      integer, dimension(7) :: srt,cnt
!@var var_name string containing name of current output field
      character(len=sname_strlen) :: var_name='' ! is reset to '' after each write
!@var units string containing units of current output field
      character(len=units_strlen) :: units='' ! is reset to '' after each write
!@var long_name description of current output field
      character(len=lname_strlen) :: long_name='' ! is set to '' after each write

!@param missing value to substitute for missing data
      real, parameter :: missing=-1.e30
      real*8, parameter :: missing8=-1.d30
!@var def_missing flag whether to define a missing data attribute
      logical :: def_missing=.false.
!@var write_whole_array flag whether to write the entire array
      logical :: write_whole_array=.false.

      integer, parameter :: real_att_max_size=100
!@var real_att_name name of real attribute to write to output file
      character(len=30) :: real_att_name=''
!@var real_att value(s) of real attribute to write to output file
      real, dimension(real_att_max_size) :: real_att = missing

c netcdf library will convert prog_dtype to disk_dtype
c these are the defaults for converting GCM real*8 to real*4 on disk
!@var disk_dtype data type to write to output file
      integer :: disk_dtype = nf_real ! is set to this after each write
!@var prog_dtype data type of array being written to disk
      integer :: prog_dtype = nf_double ! set to this after each write

      contains

      subroutine open_out
      status_out = nf_create (trim(outfile), nf_clobber, out_fid)
      if(status_out.ne.nf_noerr) then
         write(6,*) 'cannot create: '//trim(outfile)
         call stop_model('stopped in POUT_netcdf.f',255)
      endif
c-----------------------------------------------------------------------
c write global header
c-----------------------------------------------------------------------
      status_out=nf_put_att_text(out_fid,nf_global,'run_name',
     &     lrunid,xlabel)
      status_out=nf_put_att_text(out_fid,nf_global,'acc_period',
     &     len_trim(acc_period),acc_period)
      return
      end subroutine open_out

      subroutine def_dim_out(dim_name,dimlen)
      character(len=20) :: dim_name
      integer :: dimlen
      integer :: tmp_id
      status_out = nf_def_dim(out_fid,trim(dim_name),dimlen,tmp_id)
      return
      end subroutine def_dim_out

      subroutine set_dim_out(dim_name,dim_num)
      character(len=20) :: dim_name
      integer :: dim_num
      integer :: idim,status,ndims_fid
      character(len=20) :: dname
      if(dim_num.lt.1 .or. dim_num.gt.ndims_out)
     &     call stop_model('set_dim_out: invalid dim #',255)
      status = nf_inq_ndims(out_fid,ndims_fid)
      if(ndims_fid.eq.0) call stop_model(
     &     'set_dim_out: no dims defined',255)
      dimids_out(dim_num) = -1
      do idim=1,ndims_fid
        status = nf_inq_dimname(out_fid,idim,dname)
        if(trim(dname).eq.trim(dim_name)) then
          dimids_out(dim_num) = idim
          status = nf_inq_dimlen(out_fid,idim,dimlens_out(dim_num))
          exit
        endif
      enddo
      if(dimids_out(dim_num).eq.-1) then
         print*,dim_name
         call stop_model(
     &        'set_dim_out: invalid dim_name',255)
      endif
      return
      end subroutine set_dim_out

      subroutine close_out
      status_out = nf_close(out_fid)
      return
      end subroutine close_out

      end module ncout

      subroutine wrtgattc(att_name,att,n)
      use ncout
      implicit none
      include 'netcdf.inc'
      integer :: n
      character(len=30) :: att_name
      character(len=n) :: att
c-----------------------------------------------------------------------
c write global attribute of type character
c-----------------------------------------------------------------------
      status_out=nf_put_att_text(out_fid,nf_global,trim(att_name),
     &     len_trim(att),att)
      return
      end subroutine wrtgattc

      subroutine defarr
      use ncout
      implicit none
      include 'netcdf.inc'
      status_out = nf_redef(out_fid)
c check to make sure that variable isn't already defined in output file
      if(nf_inq_varid(out_fid,trim(var_name),varid_out)
     &                 .eq.nf_noerr) then
         write(6,*) 'defarr: variable already defined: ',trim(var_name)
         call stop_model('stopped in POUT_netcdf.f',255)
      endif
      status_out=nf_def_var(out_fid,trim(var_name),disk_dtype,
     &     ndims_out,dimids_out,varid_out)
      if(len_trim(units).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'units',len_trim(units),units)
      if(len_trim(long_name).gt.0) status_out =
     &     nf_put_att_text(out_fid,varid_out,
     &     'long_name',len_trim(long_name),long_name)
      if(len_trim(real_att_name).gt.0) status_out =
     &     nf_put_att_real(out_fid,varid_out,trim(real_att_name),
     &     nf_float,count(real_att.ne.missing),real_att)
c define missing value attribute if requested
      if(def_missing) status_out = nf_put_att_real(out_fid,varid_out,
     &     'missing_value',nf_float,1,missing)
      status_out = nf_enddef(out_fid)
      return
      end subroutine defarr

      subroutine wrtdarr(var)
      use ncout
      implicit none
      REAL*8 :: var(1)
      write_whole_array = .true.
      call wrtdarrn(var,1)
      return
      end subroutine wrtdarr
      subroutine wrtrarr(var)
      use ncout
      implicit none
      REAL*4 :: var(1)
      write_whole_array = .true.
      call wrtrarrn(var,1)
      return
      end subroutine wrtrarr
      subroutine wrtiarr(var)
      use ncout
      implicit none
      INTEGER :: var(1)
      write_whole_array = .true.
      call wrtiarrn(var,1)
      return
      end subroutine wrtiarr
      subroutine wrtcarr(var)
      use ncout
      implicit none
      character(len=*) :: var
      write_whole_array = .true.
      call wrtcarrn(var,1)
      return
      end subroutine wrtcarr

      subroutine setup_arrn(outer_pos)
      use ncout
      implicit none
      include 'netcdf.inc'
      integer :: outer_pos
      integer :: n
      if(ndims_out.lt.2 .and. .not.write_whole_array)
     &     call stop_model('setup_arrn: invalid ndims_out',255)
      if(outer_pos.le.0 .or.
     &     outer_pos.gt.dimlens_out(ndims_out)) then
         write(6,*) 'setup_arrn: invalid outer_pos'
         call stop_model('stopped in POUT_netcdf.f',255)
      endif
      do n=1,ndims_out
         srt(n) = 1
         cnt(n) = dimlens_out(n)
      enddo
      if(.not.write_whole_array) then
         srt(ndims_out) = outer_pos
         cnt(ndims_out) = 1
      endif
c define variable if outer_pos=1
      if(outer_pos.eq.1) then
         call defarr
      else
c check to make sure that var_name is already defined in output file
         if(nf_inq_varid(out_fid,trim(var_name),varid_out)
     &        .ne.nf_noerr) then
            write(6,*) 'setup_arrn: variable not yet defined: ',
     &           trim(var_name)
            if(.not.write_whole_array) write(6,*)
     &           'first call to wrtXarrn must have outer_pos=1'
            call stop_model('stopped in POUT_netcdf.f',255)
         endif
      endif
c restore defaults (was moved here from defarr)
      var_name=''
      units=''
      long_name=''
      real_att_name=''
      real_att=missing
      disk_dtype = nf_real
      prog_dtype = nf_double
      def_missing = .false.
      write_whole_array = .false.
      return
      end subroutine setup_arrn

      subroutine wrtdarrn(var,outer_pos)
c Writes the outer_pos-th element of the outermost dimension of an array.
c When outer_pos=1, defines the array as well; hence the first call
c for a particular array MUST have outer_pos=1.
c Does not currently have the capability to write any real_atts.
      use ncout
      implicit none
      include 'netcdf.inc'
      REAL*8 :: var(1)
      integer :: outer_pos
      integer :: n
      if(outer_pos.eq.1) then ! check for missing values
         prog_dtype=nf_double
         if(write_whole_array) then
            n = product(dimlens_out(1:ndims_out))
            if(any(var(1:n).eq.missing8)) def_missing=.true.
         else ! we can only check var(:,...,:,1) at this point,
c so just assume that there will be some missing values
            def_missing=.true.
         endif
      endif
      call setup_arrn(outer_pos)
      status_out = nf_put_vara_double(out_fid,varid_out,srt,cnt,var)
      return
      end subroutine wrtdarrn
      subroutine wrtrarrn(var,outer_pos)
c like wrtdarrn but for real*4 data
      use ncout
      implicit none
      include 'netcdf.inc'
      REAL*4 :: var(1)
      integer :: outer_pos
      integer :: n
      if(outer_pos.eq.1) then ! check for missing values
         prog_dtype=nf_real
         if(write_whole_array) then
            n = product(dimlens_out(1:ndims_out))
            if(any(var(1:n).eq.missing)) def_missing=.true.
         else ! we can only check var(:,...,:,1) at this point,
c so just assume that there will be some missing values
            def_missing=.true.
         endif
      endif
      call setup_arrn(outer_pos)
      status_out = nf_put_vara_real(out_fid,varid_out,srt,cnt,var)
      return
      end subroutine wrtrarrn
      subroutine wrtiarrn(var,outer_pos)
c like wrtdarrn but for integer data
      use ncout
      implicit none
      include 'netcdf.inc'
      integer :: var(1)
      integer :: outer_pos
      prog_dtype=nf_int
      call setup_arrn(outer_pos)
      status_out = nf_put_vara_int(out_fid,varid_out,srt,cnt,var)
      return
      end subroutine wrtiarrn
      subroutine wrtcarrn(var,outer_pos)
c like wrtdarrn but for character data
      use ncout
      implicit none
      include 'netcdf.inc'
      character(len=*) :: var
      integer :: outer_pos
      disk_dtype=nf_char
      prog_dtype=nf_char
      call setup_arrn(outer_pos)
      status_out = nf_put_vara_text(out_fid,varid_out,srt,cnt,var)
      return
      end subroutine wrtcarrn


      subroutine open_ij(filename,im_gcm,jm_gcm)
!@sum  OPEN_IJ opens the lat-lon binary output file
!@auth M. Kelley
      USE GEOM, only : lon_dg,lat_dg,dxyp,dxv,dyp,dxyv
      USE NCOUT, only : iu_ij,im,jm,set_dim_out,def_dim_out,units
     *     ,outfile,out_fid,ndims_out,open_out,var_name,long_name
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm
!
      character(len=30) :: dim_name

      outfile = trim(filename)//".nc"

! define output file
      call open_out
      iu_ij = out_fid

C**** set dimensions
      im=im_gcm
      jm=jm_gcm

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='lonb'; call def_dim_out(dim_name,im)
      dim_name='latb'; call def_dim_out(dim_name,jm-1)

      ndims_out = 1
c primary grid
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtdarr(lon_dg(1,1))
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtdarr(lat_dg(1,1))
      units='m2'; long_name='primary gridbox area'
      var_name='area';call wrtdarr(dxyp)
      units='m'
      long_name='north-south grid length'
      var_name='dy';call wrtdarr(dyp)
c secondary grid
      dim_name='lonb'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='lonb'; call wrtdarr(lon_dg(1,2))
      dim_name='latb'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latb'; call wrtdarr(lat_dg(2,2))
      units='m2'; long_name='secondary gridbox area'
      var_name='areab';call wrtdarr(dxyv(2))
      units='m'
      long_name='east-west grid length at secondary latitudes'
      var_name='dxb';call wrtdarr(dxv(2))

      return
      end subroutine open_ij

      subroutine close_ij
!@sum  CLOSE_IJ closes the lat-lon binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_ij
      call close_out

      return
      end subroutine close_ij

      subroutine POUT_IJ(TITLE,SNAME,LNAME,UNITS_IN,XIJ,XJ,XSUM,
     &     IGRID,JGRID)
!@sum  POUT_IJ output lat-lon binary records
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var XIJ lat/lon output field
      REAL*8, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field
      REAL*8, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field
      REAL*8, INTENT(IN) :: XSUM
!@var IGRID,JGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IGRID,JGRID

      character(len=30) :: lon_name,lat_name

      out_fid = iu_ij

! (re)set shape of output array
      ndims_out = 2
      if(igrid.eq.1) then
         lon_name='longitude'
      else if(igrid.eq.2) then
         lon_name='lonb'
      else
         call stop_model('pout_ij: unrecognized lon grid',255)
      endif
      if(jgrid.eq.1) then
         lat_name='latitude'
      else if(jgrid.eq.2) then
         lat_name='latb'
      else
         call stop_model('pout_ij: unrecognized lat grid',255)
      endif
      call set_dim_out(lon_name,1)
      call set_dim_out(lat_name,2)

      var_name=sname
      long_name=lname
      units=units_in
      real_att_name='glb_mean'
      real_att(1)=xsum
      if(jgrid.eq.1) then
         call wrtdarr(xij(1,1))
      else
         call wrtdarr(xij(1,2)) ! ignore first latitude row
      endif
      return
      end

      subroutine open_wp(filename,jm_gcm,lm_gcm,lm_req_gcm,lat_dg_gcm)
!@sum  OPEN_WP opens the wave power binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var LM_GCM,JM_GCM,lm_req_gcm dimensions for jl output
      INTEGER, INTENT(IN) :: lm_gcm,jm_gcm,lm_req_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM,2) :: lat_dg_gcm
!
      character(len=30) :: att_name
      integer, parameter :: lenatt=300
      character(len=lenatt) :: att_str
      character(len=1), parameter :: nl=achar(10)

      outfile = trim(filename)//".nc"
      call open_out
      iu_wp = out_fid

      att_name='note'
      att_str=nl//
     & 'Wave power diagnostics are not yet available in netcdf'//nl//
     & 'format.  If you are interested in this capability, please'//nl//
     & 'contact the Model E development team.'
      call wrtgattc(att_name,att_str,lenatt)

      return
      end subroutine open_wp

      subroutine close_wp
!@sum  CLOSE_WP closes the wave power binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_wp
      call close_out

      return
      end subroutine close_wp

      subroutine POUT_WP(TITLE,LNAME,SNAME,UNITS_IN,
     &     J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_WP output wave power in binary format
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var J1,KLMAX variables required by pout_jl
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field, dimensioned to accomodate pout_jl expectations
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
!@var PM the "vertical" coordinate for pout_jl
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM
!@var CX,CY x,y coordinate names
      CHARACTER*16, INTENT(IN) :: CX,CY
      return
      end subroutine pout_wp

      subroutine open_jl(filename,jm_gcm,lm_gcm,lm_req_gcm,lat_dg_gcm)
!@sum  OPEN_JL opens the lat-height binary output file
!@auth M. Kelley
      USE DIAG_COM, only : plm,ple,ple_dn,pmb,kgz_max,zoc,zoc1
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      character(len=30) :: dim_name
!@var LM_GCM,JM_GCM,lm_req_gcm dimensions for jl output
      INTEGER, INTENT(IN) :: lm_gcm,jm_gcm,lm_req_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM,2) :: lat_dg_gcm

      outfile = trim(filename)//".nc"
      call open_out
      iu_jl = out_fid
      jlfile = outfile

C**** set dimensions
      jm=jm_gcm
      lm=lm_gcm
      lm_req=lm_req_gcm
      lat_dg(1:JM,:)=lat_dg_gcm(1:JM,:)

      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='latb'; call def_dim_out(dim_name,jm-1)

      if(index(filename,'.ojl').gt.0) then
        dim_name='odepth'; call def_dim_out(dim_name,lm)
        dim_name='odepth1'; call def_dim_out(dim_name,lm+1)
      else
        dim_name='p'; call def_dim_out(dim_name,lm)
        dim_name='ple'; call def_dim_out(dim_name,lm-1)
        if(lm_req.gt.0) then
          dim_name='prqt'; call def_dim_out(dim_name,lm+lm_req)
        endif
        dim_name='ple_up'; call def_dim_out(dim_name,lm)
        dim_name='ple_dn'; call def_dim_out(dim_name,lm)
        dim_name='pgz'; call def_dim_out(dim_name,kgz_max)
      endif

! put lat,ht into output file
      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude'; call wrtdarr(lat_dg(1,1))
      dim_name='latb'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latb'; call wrtdarr(lat_dg(2,2))

      if(index(filename,'.ojl').gt.0) then
        dim_name='odepth'; call set_dim_out(dim_name,1)
        units='m'
        var_name='odepth'; call wrtdarr(zoc)
        dim_name='odepth1'; call set_dim_out(dim_name,1)
        units='m'
        var_name='odepth1'; call wrtdarr(zoc1)
      else
        dim_name='p'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='p'; call wrtdarr(plm)
        dim_name='ple'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='ple'; call wrtdarr(ple)
        if(lm_req.gt.0) then
          dim_name='prqt'; call set_dim_out(dim_name,1)
          units='mb'
          var_name='prqt'; call wrtdarr(plm)
        endif
        dim_name='ple_up'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='ple_up'; call wrtdarr(ple)
        dim_name='ple_dn'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='ple_dn'; call wrtdarr(ple_dn)
        dim_name='pgz'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='pgz'; call wrtdarr(pmb(1:kgz_max))
      endif

      return
      end subroutine open_jl

      subroutine close_jl
!@sum  CLOSE_JL closes the lat-height binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_jl
      call close_out

      return
      end subroutine close_jl

      subroutine POUT_JL(TITLE,LNAME,SNAME,UNITS_IN,
     &     J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_JL output lat-height binary records
!@auth M. Kelley
      USE RESOLUTION, only : LS1
      USE DIAG_COM, only : plm,ple,ple_dn,kgz_max,zoc,zoc1
      USE NCOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var KLMAX max level to output
!@var J1 minimum j value to output (needed for secondary grid fields)
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field
!@+       (J1:JM,1:KLMAX) is field
!@+       (JM+1:JM+3,1:KLMAX) are global/NH/SH average over L
!@+       (J1:JM+3,LM+LM_REQ+1) are averages over J
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
      REAL*8, DIMENSION(JM,LM+LM_REQ) :: XJL0
      REAL*8, DIMENSION(JM-1,LM+LM_REQ) :: XJL0B
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM

      character(len=30) :: dim_name

      CHARACTER*16, INTENT(IN) :: CX,CY
      INTEGER J,L
      logical :: is_ocn

      is_ocn = index(jlfile,'.ojl').gt.0

      out_fid = iu_jl

! (re)set default shape of output arrays
      ndims_out = 2

      if(j1.eq.1) then
         dim_name='latitude'
      else if(j1.eq.2) then
         dim_name='latb'
      else
         call stop_model('pout_jl: unrecognized latitude grid',255)
      endif
      call set_dim_out(dim_name,1)

      if(is_ocn .and. all(pm(1:klmax).eq.zoc(1:klmax))) then
         dim_name='odepth'
      else if(is_ocn .and. all(pm(1:klmax).eq.zoc1(1:klmax))) then
         dim_name='odepth1'
      else if(lm_req.gt.0 .and. klmax.eq.lm+lm_req) then
         dim_name='prqt'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.plm(1:lm))) then
         dim_name='p'
      else if(klmax.eq.ls1-1 .and.all(pm(1:ls1-1).eq.plm(1:ls1-1))) then
! some arrays only defined in troposphere, extend to stratosphere anyway
         dim_name='p'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple(1:lm))) then
         dim_name='ple_up'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple_dn(1:lm))) then
         dim_name='ple_dn'
      else if(klmax.eq.lm-1 .and. all(pm(1:lm-1).eq.ple(1:lm-1))) then
         dim_name='ple'
      else if(klmax.eq.kgz_max) then
         dim_name='pgz'
      else if(klmax.eq.1) then
 ! sometimes this routine is called to write out 1-dimensional arrays
         ndims_out = 1
      else
         write(6,*) 'klmax =',klmax,title,pm(1:klmax),ple(1:klmax)
         call stop_model('pout_jl: unrecognized vertical grid',255)
      endif
      if(ndims_out.eq.2) call set_dim_out(dim_name,2)

      var_name=sname
      long_name=lname
      units=units_in

      if(.not. is_ocn) then
        real_att_name='g-nh-sh_sums-means'
        real_att(1:3)=XJL(JM+1:JM+3,LM+LM_REQ+1)
      endif

      if(j1.eq.1) then
         xjl0(1:jm,1:lm+lm_req) = xjl(1:jm,1:lm+lm_req)
         call wrtdarr(xjl0)
      else
         xjl0b(1:jm-1,1:lm+lm_req) = xjl(2:jm,1:lm+lm_req)
         call wrtdarr(xjl0b)
      endif

      return
      end

      subroutine open_il(filename,im_gcm,lm_gcm,lm_req_gcm)
!@sum  OPEN_IL opens the lon-height binary output file
!@auth M. Kelley
      USE GEOM, only : lon_dg
      USE DIAG_COM, only : plm,ple,ple_dn,zoc,zoc1
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      character(len=30) :: dim_name
!@var IM_GCM,LM_GCM,lm_req_gcm dimensions for il output
      INTEGER, INTENT(IN) :: im_gcm,lm_gcm,lm_req_gcm

      outfile = trim(filename)//".nc"
      call open_out
      iu_il = out_fid
      ilfile = outfile

      im=im_gcm
      lm=lm_gcm
      lm_req=lm_req_gcm

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='lonb'; call def_dim_out(dim_name,im)

      if(index(filename,'.oil').gt.0) then
        dim_name='odepth'; call def_dim_out(dim_name,lm)
        dim_name='odepth1'; call def_dim_out(dim_name,lm+1)
      else
        dim_name='p'; call def_dim_out(dim_name,lm)
        dim_name='prqt'; call def_dim_out(dim_name,lm+lm_req)
        dim_name='ple_up'; call def_dim_out(dim_name,lm)
        dim_name='ple_dn'; call def_dim_out(dim_name,lm)
        dim_name='ple_int'; call def_dim_out(dim_name,lm-1)
      endif

! put lon,ht into output file
      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude'; call wrtdarr(lon_dg(1,1))
      dim_name='lonb'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='lonb'; call wrtdarr(lon_dg(1,2))

      if(index(filename,'.oil').gt.0) then
        dim_name='odepth'; call set_dim_out(dim_name,1)
        units='m'
        var_name='odepth'; call wrtdarr(zoc)
        dim_name='odepth1'; call set_dim_out(dim_name,1)
        units='m'
        var_name='odepth1'; call wrtdarr(zoc1)
      else
        dim_name='p'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='p'; call wrtdarr(plm)
        dim_name='prqt'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='prqt'; call wrtdarr(plm)
        dim_name='ple_up'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='ple_up'; call wrtdarr(ple)
        dim_name='ple_dn'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='ple_dn'; call wrtdarr(ple_dn)
        dim_name='ple_int'; call set_dim_out(dim_name,1)
        units='mb'
        var_name='ple_int'; call wrtdarr(ple)
      endif

      return
      end subroutine open_il

      subroutine close_il
!@sum  CLOSE_IL closes the lon-height binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_il
      call close_out

      return
      end subroutine close_il

      subroutine POUT_IL(TITLE,sname,lname,unit,I1,ISHIFT,KLMAX,XIL
     *     ,PM,CX,CY,ASUM,GSUM,ZONAL)
!@sum  POUT_IL output lon-height binary records
!@auth M. Kelley
      USE GEOM, only : lon_dg
      USE DIAG_COM, only : plm,ple,ple_dn,zoc,zoc1
      USE NCOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var ISHIFT flag for secondary grid
!@var I1 coordinate index associated with first long. (for wrap-around)
      INTEGER, INTENT(IN) :: KLMAX,ISHIFT,I1
!@var XIL output field
      REAL*8, DIMENSION(IM,KLMAX+1), INTENT(INOUT) :: XIL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(KLMAX), INTENT(IN) :: PM
!@var ASUM vertical mean/sum
      REAL*8, DIMENSION(IM), INTENT(IN) :: ASUM
!@var GSUM total sum/mean
      REAL*8, INTENT(IN) :: GSUM
!@var ZONAL zonal sum/mean
      REAL*8, DIMENSION(KLMAX), INTENT(IN) :: ZONAL

      CHARACTER*16, INTENT(IN) :: CX,CY
      INTEGER I,L
      REAL*8 XTEMP(IM,KLMAX)
      character(len=30) :: dim_name
      character(len=sname_strlen), intent(in) :: sname
      character(len=units_strlen), intent(in) :: unit
      character(len=lname_strlen), intent(in) :: lname
      logical :: is_ocn

      is_ocn = index(ilfile,'.oil').gt.0

      out_fid = iu_il

! (re)set shape of output arrays
      ndims_out = 2

C**** Note that coordinates cannot be wrapped around as a function I1,
C**** therefore shift array back to standard order.
      if (i1.gt.1) then
        xtemp(:,1:klmax)=xil(:,1:klmax)
        xil(1:i1-1,1:klmax) = xtemp(im-i1+2:im,1:klmax)
        xil(i1:im,1:klmax) = xtemp(1:im-i1+1,1:klmax)
      end if

      if(ishift.eq.1) then
         dim_name='longitude'
      else if(ishift.eq.2) then
         dim_name='lonb'
      else
         call stop_model('pout_il: unrecognized longitude grid',255)
      endif
      call set_dim_out(dim_name,1)

      if(is_ocn .and. all(pm(1:klmax).eq.zoc(1:klmax))) then
         dim_name='odepth'
      else if(is_ocn .and. all(pm(1:klmax).eq.zoc1(1:klmax))) then
         dim_name='odepth1'
      elseif(klmax.eq.lm+lm_req) then
         dim_name='prqt'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.plm(1:lm))) then
         dim_name='p'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple(1:lm))) then
         dim_name='ple_up'
      else if(klmax.eq.lm .and. all(pm(1:lm).eq.ple_dn(1:lm))) then
         dim_name='ple_dn'
      else if(klmax.eq.lm-1 .and. all(pm(1:lm-1).eq.ple(1:lm-1))) then
         dim_name='ple_int'
      else
         write(6,*) 'klmax =',klmax,title,pm(1:klmax)
         call stop_model('pout_il: unrecognized vertical grid',255)
      endif
      call set_dim_out(dim_name,2)

      var_name=sname
      long_name=lname
      units=unit

      if(.not. is_ocn) then
        real_att_name='mean'
        real_att(1)=gsum
      endif

      call wrtdarr(xil)

      return
      end

      subroutine open_j(filename,ntypes,jm_gcm,lat_dg_gcm)
!@sum  OPEN_J opens the latitudinal binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var ntypes number of surface types to be output
      integer, intent(in) :: ntypes
!@var JM_GCM dimensions for j output
      INTEGER, INTENT(IN) :: jm_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM) :: lat_dg_gcm

      character(len=30) :: dim_name

      outfile = trim(filename)//".nc"
      call open_out
      iu_j = out_fid

C**** set dimensions
      jm=jm_gcm
      lat_dg(1:JM,1)=lat_dg_gcm(1:JM)

      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='stype'; call def_dim_out(dim_name,ntypes)
      dim_name='stype_clen'; call def_dim_out(dim_name,16)

      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtdarr(lat_dg(1,1))

      return
      end subroutine open_j

      subroutine close_j
!@sum  CLOSE_J closes the latitudinal binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_j
      call close_out

      return
      end subroutine close_j

      subroutine POUT_J(TITLE,SNAME,LNAME,UNITS_IN,BUDG,KMAX,TERRAIN,
     *     iotype)
!@sum  POUT_J output zonal budget file
!@auth M. Kelley
      USE DIAG_COM, only : KAJ
      USE NCOUT
      IMPLICIT NONE
c
      CHARACTER*16, DIMENSION(KAJ),INTENT(INOUT) :: TITLE
!@var LNAME,SNAME,UNITS_IN information strings for netcdf
      CHARACTER(len=lname_strlen), DIMENSION(KAJ),INTENT(IN) :: LNAME
      CHARACTER(len=sname_strlen), DIMENSION(KAJ),INTENT(IN) :: SNAME
      CHARACTER(len=units_strlen), DIMENSION(KAJ),INTENT(IN) :: UNITS_IN
      CHARACTER*16, INTENT(IN) :: TERRAIN
      REAL*8, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
      INTEGER, INTENT(IN) :: KMAX,iotype
      INTEGER K,N,J

      character(len=30) :: dim_name

      out_fid = iu_j

! write out surface type name
      ndims_out = 2
      dim_name='stype_clen'; call set_dim_out(dim_name,1)
      dim_name='stype'; call set_dim_out(dim_name,2)
      var_name='stype_names';call wrtcarrn(terrain,iotype)


! (re)set shape of output array
      ndims_out = 2
      dim_name='latitude'; call set_dim_out(dim_name,1)
      dim_name='stype'; call set_dim_out(dim_name,2)

      DO K=1,KMAX
        var_name=sname(k)
        long_name=lname(k)
        units=units_in(k)
        call wrtdarrn(budg(1,k),iotype)
      END DO

      return
      end

      subroutine open_ijk(filename,im_gcm,jm_gcm,lm_gcm)
!@sum  OPEN_IJK opens the lat-lon-height binary output file
!@auth M. Kelley
      USE GEOM, only : lon_dg,lat_dg
      USE DIAG_COM, only : plm,ple
      USE NCOUT, only : im,jm,lm,iu_ijk,set_dim_out,def_dim_out,out_fid
     *     ,outfile,units,ndims_out,open_out,var_name
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM,LM_GCM dimensions for ijk output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm,lm_gcm
!
      character(len=30) :: dim_name

      outfile = trim(filename)//".nc"

! define output file
      call open_out
      iu_ijk = out_fid

C**** set dimensions
      im=im_gcm
      jm=jm_gcm
      lm=lm_gcm

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='lonb'; call def_dim_out(dim_name,im)
      dim_name='latb'; call def_dim_out(dim_name,jm-1)
      dim_name='p'; call def_dim_out(dim_name,lm)
      dim_name='ple'; call def_dim_out(dim_name,lm)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtdarr(lon_dg(1,1))
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtdarr(lat_dg(1,1))
      dim_name='lonb'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='lonb';call wrtdarr(lon_dg(1,2))
      dim_name='latb'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latb';call wrtdarr(lat_dg(2,2))
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='p'; call wrtdarr(plm)
      dim_name='ple'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple'; call wrtdarr(ple)
      return
      end subroutine open_ijk

      subroutine close_ijk
!@sum  CLOSE_IJK closes the lat-lon-height binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_ijk
      call close_out

      return
      end subroutine close_ijk

      subroutine POUT_IJK(TITLE,SNAME,LNAME,UNITS_IN,XIJK,XJK,XK
     &     ,IJKGRID)
!@sum  POUT_IJK output lat-lon-height binary output file
!@auth M. Kelley
      USE NCOUT
      USE DIAG_COM, only : igride,jgride,kgride
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(INOUT) :: XIJK
!@var XJK lat sum/mean of output field
      REAL*8, DIMENSION(JM,LM), INTENT(INOUT) :: XJK
!@var XK global sum/mean of output field
      REAL*8, DIMENSION(LM), INTENT(IN) :: XK
!@var IJKGRID contains IGRID,JGRID,KGRID in successive bits
      INTEGER, INTENT(IN) :: IJKGRID
!@var IGRID,JGRID,KGRID = 1 for centers, 2 for edges
      INTEGER :: IGRID,JGRID,KGRID,IJKG

      character(len=30) :: dim_name

      integer :: j,l,jpack,lpack

      out_fid = iu_ijk

      ijkg = ijkgrid
      kgrid = ijkg/kgride
      ijkg = ijkg - kgride*kgrid
      jgrid = ijkg/jgride
      ijkg = ijkg - jgride*jgrid
      igrid = ijkg

      if (jgrid.eq.1) then
c pack first-j-empty :,jm,lm array to memory-contiguous :,jm-1,lm array
      lpack = 1
      jpack = 0
      do l=1,lm
      do j=2,jm
         jpack = jpack + 1
         if(jpack.gt.jm) then
            jpack = 1
            lpack = lpack + 1
         endif
         xijk(:,jpack,lpack) = xijk(:,j,l)
         xjk(jpack,lpack) = xjk(j,l)
      enddo
      enddo
      end if

! (re)set shape of output array
      ndims_out = 3

      if (igrid.eq.0) then
        dim_name = 'longitude'
      else
        dim_name = 'lonb'
      end if
      call set_dim_out(dim_name,1)

      if (jgrid.eq.0) then
        dim_name = 'latitude'
      else
        dim_name = 'latb'
      end if
      call set_dim_out(dim_name,2)

      if(kgrid.eq.0) then
        dim_name = 'p'
      else
        dim_name = 'ple'
      endif
      call set_dim_out(dim_name,3)

      var_name=sname
      long_name=lname
      units=units_in

      call wrtdarr(xijk)

      if (jgrid.eq.1) then
c unpack memory-contiguous :,jm-1,lm memory to first-j-empty :,jm,lm array
      lpack = 1 + ((jm-1)*lm)/jm
      jpack = (jm-1)*lm - (lpack-1)*jm
      do l=lm,1,-1
      do j=2,jm
         xijk(:,j,l) = xijk(:,jpack,lpack)
         xjk(j,l) = xjk(jpack,lpack)
         jpack = jpack - 1
         if(jpack.lt.1) then
            jpack = jm
            lpack = lpack - 1
         endif
      enddo
      enddo
      end if

      return
      end

      subroutine open_ijl(filename,im_gcm,jm_gcm,lm_gcm,
     &     kaijl,name_ijl,lname_ijl,units_ijl,lgrid_ijl
     &     )
!@sum  OPEN_IJL opens the lat-lon-layer binary output file
!@auth M. Kelley
      USE GEOM, only : lon_dg_gcm=>lon_dg,lat_dg_gcm=>lat_dg
      USE NCOUT
      USE DIAG_COM, only :
     &     plm,ple,ctr_ml,edg_ml,ctr_cp,edg_cp
      USE MDIAG_COM, only :
     &     sname_strlen,lname_strlen,units_strlen
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM,LM_GCM dimensions for ijl output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm,lm_gcm
!@var kaijl number of predefined output fields
      integer, intent(in) :: kaijl
      CHARACTER(len=lname_strlen), DIMENSION(kaijl) :: lname_ijl
      CHARACTER(len=sname_strlen), DIMENSION(kaijl) :: name_ijl
      CHARACTER(len=units_strlen), DIMENSION(kaijl) :: units_ijl
      integer, dimension(kaijl) :: lgrid_ijl
      INTEGER :: k,l
      REAL*8 :: lev(lm_gcm)
      include 'netcdf.inc'
!
      character(len=30) :: dim_name

      outfile = trim(filename)//".nc"

! define output file
      call open_out
      iu_ijl = out_fid

c
c Set dimensions.  All output fields are assumed to be on the
c primary grid.
c
      im=im_gcm
      jm=jm_gcm
      lm=lm_gcm

      dim_name='longitude'; call def_dim_out(dim_name,im)
      dim_name='latitude'; call def_dim_out(dim_name,jm)
      dim_name='level'; call def_dim_out(dim_name,lm)
      dim_name='p'; call def_dim_out(dim_name,lm)
      dim_name='ple'; call def_dim_out(dim_name,lm)

      ndims_out = 1
      dim_name='longitude'; call set_dim_out(dim_name,1)
      units='degrees_east'
      var_name='longitude';call wrtdarr(lon_dg_gcm(1,1))
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtdarr(lat_dg_gcm(1,1))
      dim_name='level'; call set_dim_out(dim_name,1)
      units=' '
      do l=1,lm
        lev(l)=l
      end do
      var_name='level'; call wrtdarr(lev)
      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='p'; call wrtdarr(plm)
      dim_name='ple'; call set_dim_out(dim_name,1)
      units='mb'
      var_name='ple'; call wrtdarr(ple)

c
c predefine the output fields
c
      ndims_out = 3
      dim_name='longitude'; call set_dim_out(dim_name,1)
      dim_name='latitude'; call set_dim_out(dim_name,2)
      disk_dtype = nf_real
      do k=1,kaijl
        if(trim(lname_ijl(k)).eq.'no output') cycle
        if(lgrid_ijl(k).eq.ctr_ml .or. lgrid_ijl(k).eq.edg_ml) then
          dim_name = 'level'
          def_missing = .false.
        elseif(lgrid_ijl(k).eq.ctr_cp) then
          dim_name = 'p'
          def_missing = .true.
        else
          dim_name = 'ple'
          def_missing = .true.
        endif
        call set_dim_out(dim_name,3)

        var_name=name_ijl(k)
        long_name=lname_ijl(k)
        units=units_ijl(k)

        call defarr
      enddo
c restore some defaults
      def_missing = .false.

      return
      end subroutine open_ijl

      subroutine close_ijl
!@sum  CLOSE_IJL closes the lat-lon-layer binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_ijl
      call close_out

      return
      end subroutine close_ijl

      subroutine POUT_IJL(TITLE,SNAME,LNAME,UNITS_IN,XIJL,XJL,XL,IJGRID)
!@sum  POUT_IJL output lat-lon-layer binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(INOUT) :: XIJL
!@var XJK lat sum/mean of output field
      REAL*8, DIMENSION(JM,LM), INTENT(INOUT) :: XJL
!@var XK global sum/mean of output field
      REAL*8, DIMENSION(LM), INTENT(IN) :: XL
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID

      character(len=30) :: dim_name
      include 'netcdf.inc'
      integer :: status

      out_fid = iu_ijl

      var_name = sname
      status = nf_inq_varid(out_fid,trim(var_name),varid_out)
      status = nf_put_var_double(out_fid,varid_out,xijl)

      return
      end

      subroutine open_isccp(filename,ntau,npres,nisccp)
!@sum  OPEN_ISCCP opens the binary output file of ISCCP histograms
!@auth M. Kelley
      USE DIAG_COM, only : isccp_press,isccp_taum,isccp_lat
      USE NCOUT, only : im,jm,lm,iu_isccp,set_dim_out,def_dim_out,
     &     out_fid,outfile,units,long_name,ndims_out,open_out
     &     ,var_name
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var ntau,npres,nisccp dimensions for isccp output
      INTEGER, INTENT(IN) :: ntau,npres,nisccp
!
      character(len=30) :: dim_name

      outfile = trim(filename)//".nc"

! define output file
      call open_out
      iu_isccp = out_fid

C**** set dimensions
      im=ntau
      jm=npres
      lm=nisccp

      dim_name='tau'; call def_dim_out(dim_name,im)
      dim_name='p'; call def_dim_out(dim_name,jm)
      dim_name='lat'; call def_dim_out(dim_name,lm)

      ndims_out = 1

      dim_name='tau'; call set_dim_out(dim_name,1)
      units='1'
      long_name='mid-point of optical depth for each tau category'
      var_name='tau';call wrtdarr(isccp_taum)

      dim_name='p'; call set_dim_out(dim_name,1)
      units='mb'
      long_name='midpoint pressure for each pressure category'
      var_name='p';call wrtiarr(isccp_press)

      dim_name='lat'; call set_dim_out(dim_name,1)
      units='degrees_north'
      long_name='midpoint latitude for each latitude category'
      var_name='lat'; call wrtdarr(isccp_lat)

      return
      end subroutine open_isccp

      subroutine close_isccp
!@sum  CLOSE_ISCCP closes the binary output file of ISCCP histograms
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_isccp
      call close_out

      return
      end subroutine close_isccp

      subroutine POUT_ISCCP(TITLE,SNAME,LNAME,UNITS_IN,XIJK,TAUM,PRES)
!@sum  POUT_ISCCP outputs tau-height-lat binary output file of ISCCP histograms
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var XIJK tau/height/lat output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(IN) :: XIJK
      REAL*8, INTENT(IN) :: taum(IM), pres(JM)

      character(len=30) :: dim_name

      out_fid = iu_isccp

! (re)set shape of output array
      ndims_out = 3

      dim_name = 'tau'
      call set_dim_out(dim_name,1)
      dim_name = 'p'
      call set_dim_out(dim_name,2)
      dim_name = 'lat'
      call set_dim_out(dim_name,3)

      var_name=sname
      long_name=lname
      units=units_in

      call wrtdarr(xijk)

      return
      end subroutine pout_isccp

      subroutine open_diurn(filename,hr_in_period,NDIUVAR_gcm,kr1,kr2)
!@sum  OPEN_DIURN opens the average diurnal cycle netcdf output file
!@auth M. Kelley
      USE DIAG_COM, only : namdd,ijdd
      USE NCOUT, only : im,jm,lm,im_data,iu_diurn,set_dim_out
     &   ,def_dim_out,out_fid,outfile,units,long_name,ndims_out,open_out
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: kr1,kr2
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: hr_in_period,NDIUVAR_gcm
!
      character(len=3) :: istr,jstr
      character(len=30) :: dim_name,att_name
      character(len=100) :: loc_info
      integer :: n

      outfile = trim(filename)//".nc"

! define output file
      call open_out
      iu_diurn = out_fid

C**** set dimensions
      im=hr_in_period
      im_data=hr_in_period+1 ! default; change before calling pout_diurn
      jm=NDIUVAR_gcm

      dim_name='hour'; call def_dim_out(dim_name,im)

! write the names and i,j indices of diagnostic locations
      loc_info=''
      do n=kr1,kr2
         write(istr,'(i3)') ijdd(1,n)
         istr = adjustl(istr)
         write(jstr,'(i3)') ijdd(2,n)
         jstr = adjustl(jstr)
         loc_info = trim(loc_info)//' '//namdd(n)//'='//
     &        trim(istr)//','//trim(jstr)
      enddo
      att_name='locations'
      call wrtgattc(att_name,loc_info,100)

c      ndims_out = 1
c      dim_name='hour'; call set_dim_out(dim_name,1)
c      units='1'
c      var_name='hour';call wrtdarr(hours)

      return
      end subroutine open_diurn

      subroutine close_diurn
!@sum  CLOSE_DIURN closes the average diurnal cycle netcdf output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_diurn
      call close_out

      return
      end subroutine close_diurn

      subroutine POUT_diurn(SNAME_IN,NAME_IN,UNITS_IN,
     &     FHOUR,NAMDD,IJDD1,IJDD2,HR_IN_PERIOD,kp)
!@sum  POUT_diurn outputs the average diurnal_cycle netcdf output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var SNAME,NAME,UNITS short name, long name, and units strings
      CHARACTER*16, DIMENSION(jm) :: UNITS_IN,NAME_IN,SNAME_IN
      CHARACTER*4, INTENT(IN) :: NAMDD ! name of current output region
      INTEGER, INTENT(IN) :: HR_IN_PERIOD,KP,IJDD1,IJDD2
      REAL*8, DIMENSION(im_data,jm), INTENT(IN) :: FHOUR

      character(len=30) :: dim_name
      integer :: k

      out_fid = iu_diurn

! (re)set shape of output array
      ndims_out = 1

      dim_name = 'hour'
      call set_dim_out(dim_name,1)

      do k=1,kp
         var_name = namdd//'_'//trim(sname_in(k))
         long_name=trim(name_in(k))
         units=trim(units_in(k))
         call wrtdarr(fhour(1,k))
      enddo

      return
      end subroutine POUT_diurn

      subroutine open_hdiurn(filename,hr_in_month,NDIUVAR_gcm,kr1,kr2)
!@sum  OPEN_HDIURN opens the hour-by-hour history netcdf output file
!@auth M. Kelley
! note: hdiurn output routines simply call diurn output routines;
! to avoid file unit conflicts, open_hdiurn cannot be called if
! any files are currently open through open_diurn
      use NCOUT
      use MODEL_COM, only : idacc ! get # of hours in this month
      use DIAG_COM, only : ia_12hr,hr_in_day
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: kr1,kr2
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: hr_in_month,NDIUVAR_gcm
      character(len=30) :: att_name
      integer, parameter :: lenatt=300
      character(len=lenatt) :: att_str
      character(len=1), parameter :: nl=achar(10)
      integer :: nhrs
!
      nhrs = min(hr_in_month,hr_in_day*(idacc(ia_12hr)/2))
      call open_diurn(filename,nhrs,NDIUVAR_gcm,kr1,kr2)
      im_data = hr_in_month ! input arrays always dimensioned the same

      att_name='note'
      att_str=nl//
     & 'Variables calculated in the radiation routine may have'//nl//
     & 'zero values at the beginning of the month if that routine'//nl//
     & 'is not called every hour.  Fill in with the last value'//nl//
     & 'from the previous month.'
      call wrtgattc(att_name,att_str,lenatt)

      return
      end subroutine open_hdiurn

      subroutine close_hdiurn
!@sum  CLOSE_HDIURN closes the hour-by-hour history netcdf output file
!@auth M. Kelley
      IMPLICIT NONE

      call close_diurn
      return
      end subroutine close_hdiurn

      subroutine POUT_hdiurn(SNAME_IN,NAME_IN,UNITS_IN,
     &     FHOUR,NAMDD,IJDD1,IJDD2,HR_IN_PERIOD,kp)
!@sum  POUT_hdiurn output hour-by-hour history netcdf output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var SNAME,NAME,UNITS short name, long name, and units strings
      CHARACTER*16, DIMENSION(jm) :: UNITS_IN,NAME_IN,SNAME_IN
      CHARACTER*4, INTENT(IN) :: NAMDD ! name of current output region
      INTEGER, INTENT(IN) :: HR_IN_PERIOD,KP,IJDD1,IJDD2
      REAL*8, DIMENSION(im_data,jm), INTENT(IN) :: FHOUR

      call pout_diurn(SNAME_IN,NAME_IN,UNITS_IN,
     &     FHOUR,NAMDD,IJDD1,IJDD2,HR_IN_PERIOD,kp)

      return
      end subroutine POUT_hdiurn

      subroutine open_jc(filename,jm_gcm,lat_dg_gcm)
!@sum  OPEN_jc opens the conservation quantity binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var JM_GCM dimensions for j output
      INTEGER, INTENT(IN) :: jm_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM) :: lat_dg_gcm

      character(len=30) :: dim_name

      outfile = trim(filename)//".nc"
      call open_out
      iu_jc = out_fid

C**** set dimensions
      jm=jm_gcm
      lat_dg(1:JM,1)=lat_dg_gcm(1:JM)

      dim_name='latitude'; call def_dim_out(dim_name,jm)

      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)
      units='degrees_north'
      var_name='latitude';call wrtdarr(lat_dg(1,1))

      return
      end subroutine open_jc

      subroutine close_jc
!@sum  CLOSE_JC closes the conservation quantity binary output file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE

      out_fid = iu_jc
      call close_out

      return
      end subroutine close_jc

      subroutine POUT_jc(TITLE,SNAME,LNAME,UNITS_IN,cnslat,KMAX)
!@sum  POUT_JC output zonal conservation binary file
!@auth M. Kelley
      USE NCOUT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KMAX
      CHARACTER*38, DIMENSION(kmax),INTENT(INOUT) :: TITLE
!@var LNAME,SNAME,UNITS dummy strings
      CHARACTER(len=lname_strlen), DIMENSION(kmax),INTENT(IN) :: LNAME
      CHARACTER(len=sname_strlen), DIMENSION(kmax),INTENT(IN) :: SNAME
      CHARACTER(len=units_strlen), DIMENSION(kmax),INTENT(IN) ::UNITS_IN
      REAL*8, DIMENSION(JM+3,kmax), INTENT(IN) :: cnslat
      INTEGER :: K
      character(len=30) :: dim_name

      out_fid = iu_jc

! (re)set shape of output array
      ndims_out = 1
      dim_name='latitude'; call set_dim_out(dim_name,1)

      DO K=1,KMAX
        var_name=sname(k)
        long_name=lname(k)
        units=units_in(k)
        real_att_name='sh-nh-g_sums-means'
        real_att(1:3)=cnslat(jm+1:jm+3,k)
        call wrtdarr(cnslat(1,k))
      END DO

      return
      end subroutine pout_jc
