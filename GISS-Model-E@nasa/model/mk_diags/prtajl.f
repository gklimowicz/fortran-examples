      subroutine prtajl(fid,progargs)
!@sum prtajl prints the fields in an input file whose
!@+   metadata mark them as having been created from
!@+   a modelE AJL-type array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:), allocatable :: lat_dg,vmean,plm,ple,pgz,pm
      real*4, dimension(:), allocatable ::
     &     lat_budg_dg,lat_agc_dg,lat2_agc_dg
      real*4 :: p_radonly(3)
      real*4, dimension(:,:), allocatable :: xjl,xjl_hemis,
     &     xjl_radonly,xjl_radonly_hemis
      character(len=30) :: units
      character(len=80) :: lname,title,outfile
      character(len=40) :: vname,vname_hemis,vname_vmean
      character(len=8) :: tpow
      character(len=4) :: dash='----',blank4
      character(len=20) :: latname
      real*4 :: prtfac,fglob,fnh,fsh
      integer :: j,l,jm,jm_budg,jm_agc,lm,kgz,km,inc,lstr,prtpow,linect,
     &     nargs,k1,k2,lunit,prtpow_vmean
      integer :: lats_per_zone,j1,j2,zone,nzones,lats_this_zone
      integer :: minj,maxj
      logical :: do_giss,all_lats,has_radonly,has_vmean
      real*4, parameter :: missing=-1.e30
      integer :: vid

c
c Various string formats
c
      character(len=80), parameter :: fmtlat =
     &     "('  P(MB)   MEAN G      NH      SH  ',32I4)"
      integer :: status,varid,varid_hemis,varid_vmean,nvars,dimids(2),
     &     plm_dimid,ple_dimid,pgz_dimid,
     &     lat_agc_dimid,lat2_agc_dimid,lat_budg_dimid,
     &     varid_radonly,varid_radonly_hemis
      character(len=132) :: xlabel
      character(len=100) :: fromto

      nargs = iargc()

c
c Look at command-line arguments to see whether GISS-convention
c fortran-format binary output was requested
c
      k1 = index(progargs,'gissfmt=')
      do_giss = k1.gt.0
      if(do_giss) then
        do k2=k1+8,len_trim(progargs)
          if(progargs(k2:k2).eq.' ') exit
        enddo
        outfile=progargs(k1+8:k2-1)
        lunit = 10
        open(lunit,file=trim(outfile),form='unformatted',
     &       convert='big_endian')
      endif

c
c Check whether printout of all lats was requested
c
      all_lats = index(progargs,'all_lats').gt.0

c
c get run ID, time/date info, number of latitudes and levels
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)

      latname = 'lat_budg'
      if(nf_inq_dimid(fid,trim(latname),lat_budg_dimid)==nf_noerr) then
        call get_dimsize(fid,trim(latname),jm_budg)
      else
        jm_budg = 0
        lat_budg_dimid = -99
      endif

      latname = 'lat'
      if(nf_inq_dimid(fid,trim(latname),lat_agc_dimid)==nf_noerr) then
        call get_dimsize(fid,trim(latname),jm_agc)
        status = nf_inq_dimid(fid,'lat2',lat2_agc_dimid)
      else
        jm_agc = 0
        lat_agc_dimid = -99
        lat2_agc_dimid = -99
      endif

      call get_dimsize(fid,'plm',lm)

c
c allocate workspace
c
      allocate(lat_budg_dg(jm_budg))
      allocate(lat_agc_dg(jm_agc))
      allocate(lat2_agc_dg(jm_agc))

      jm = max(jm_agc,jm_budg)
      allocate(lat_dg(jm),vmean(jm+3))
      allocate(xjl_hemis(3,lm))
      allocate(xjl_radonly(jm_budg,lm),xjl_radonly_hemis(3,lm))
      allocate(plm(lm),ple(lm),pm(lm))

      allocate(xjl(1,1)) ! just in case allocation status undefined

c
c read geometry
c
      if(lat_agc_dimid.gt.0) then ! primary agc lats
        call get_var_real(fid,'lat',lat_agc_dg)
      endif
      if(lat2_agc_dimid.gt.0) then ! secondary agc lats
        call get_var_real(fid,'lat2',lat2_agc_dg)
      endif
      if(lat_budg_dimid.gt.0) then ! budg lats
        call get_var_real(fid,'lat_budg',lat_budg_dg)
      endif

      call get_var_real(fid,'plm',plm)
      call get_var_real(fid,'ple',ple)
      status = nf_inq_dimid(fid,'plm',plm_dimid)
      status = nf_inq_dimid(fid,'ple',ple_dimid)
      status = nf_inq_dimid(fid,'pgz',pgz_dimid)
      if(status.eq.nf_noerr) then
        call get_dimsize(fid,'pgz',kgz)
        allocate(pgz(kgz))
        call get_var_real(fid,'pgz',pgz)
      endif
      if(nf_inq_varid(fid,'p_radonly',vid).eq.nf_noerr) then
        status = nf_get_var_real(fid,vid,p_radonly)
      endif

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c Loop over quantities.  An array xyz is printed if there also
c exists an array xyz_hemis containing hemispheric/global averages
c for that quantity.
c
      linect=65

      do varid_hemis=1,nvars
        status = nf_inq_varname(fid,varid_hemis,vname_hemis)
        lstr = len_trim(vname_hemis)
        if(index(vname_hemis,'_hemis').eq.0) cycle
        if(vname_hemis(lstr-5:lstr).ne.'_hemis') cycle
        vname = vname_hemis(1:lstr-6)
        vname_vmean = trim(vname)//'_vmean'
        if(index(trim(vname),'_radonly').gt.0) cycle ! print on top of reg. domain
        has_radonly = nf_noerr.eq.
     &       nf_inq_varid(fid,trim(vname)//'_radonly',varid_radonly)
        if(has_radonly) then
          status = nf_get_var_real(fid,varid_radonly,xjl_radonly)
          status = nf_inq_varid(fid,trim(vname)//'_radonly_hemis',
     &         varid_radonly_hemis)
          status = nf_get_var_real(fid,varid_radonly_hemis,
     &         xjl_radonly_hemis)
        endif
        status = nf_inq_varid(fid,trim(vname),varid)
        has_vmean = nf_noerr.eq.
     &       nf_inq_varid(fid,trim(vname_vmean),varid_vmean)
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        if(trim(units).eq.'unused') cycle
        lname = ''
        status = nf_get_att_text(fid,varid,'long_name',lname)
        prtpow = 0
        status = nf_get_att_int(fid,varid,'prtpow',prtpow)
        prtpow_vmean = 0
        status = nf_get_att_int(fid,varid,'prtpow_vmean',prtpow_vmean)


c
c retrieve horizontal and vertical coordinate info
c
        status = nf_inq_vardimid(fid,varid,dimids)
        if(dimids(1).eq.lat_agc_dimid) then
          jm = jm_agc
          lat_dg(1:jm) = lat_agc_dg(1:jm)
          minj = 1; maxj = jm
        elseif(dimids(1).eq.lat2_agc_dimid) then
          jm = jm_agc
          lat_dg(1:jm) = lat2_agc_dg(1:jm)
          if(lat_dg(1).le.-90.) then
            minj = 2; maxj = jm
          else
            minj = 1; maxj = jm-1
          endif
        elseif(dimids(1).eq.lat_budg_dimid) then
          jm = jm_budg
          lat_dg(1:jm) = lat_budg_dg(1:jm)
          minj = 1; maxj = jm
        else
          stop 'how did we get here'
        endif
        km = lm
        if(dimids(2).eq.plm_dimid) then
          pm(:) = plm(:)
        elseif(dimids(2).eq.ple_dimid) then
          pm(:) = ple(:)
        elseif(dimids(2).eq.pgz_dimid) then
          km = kgz
          pm(1:kgz) = pgz(1:kgz)
        endif


        xjl_hemis = missing
        status = nf_get_var_real(fid,varid_hemis,xjl_hemis)
        if(has_vmean) then
          status = nf_get_var_real(fid,varid_vmean,vmean)
        else
          vmean = missing
        endif
        if(allocated(xjl)) deallocate(xjl)
        allocate(xjl(jm,lm))
        status = nf_get_var_real(fid,varid,xjl)
        if(any(xjl.eq.nf_fill_real)) then
          write(6,*) 'undefined output: ',trim(vname)
          write(6,*) 'run agcstat first'
          cycle
        endif

        nzones = 1
        lats_per_zone = jm
        if(all_lats) then
          do while(lats_per_zone .gt. 32)
            nzones = nzones + 1
            lats_per_zone = jm/nzones
          enddo
          if(lats_per_zone*nzones.ne.jm) stop 'factoring error'
          inc=1
        else
          inc=1+(jm-1)/24
        endif

c
c form title string and rescale fields for ASCII output
c
        blank4 = '    '
        if(prtpow.ne.0) then
          prtfac = 10.**(-prtpow)
          where(xjl.ne.missing) xjl = xjl*prtfac
          where(xjl_hemis.ne.missing) xjl_hemis = xjl_hemis*prtfac
          write (tpow, '(i3)') prtpow
          tpow='10**'//trim(adjustl(tpow))
          units = trim(tpow)//' '//trim(units)
          if(has_radonly) then
            xjl_radonly = xjl_radonly*prtfac
            xjl_radonly_hemis = xjl_radonly_hemis*prtfac
          endif
          where(vmean.ne.missing) vmean = vmean*prtfac
          if(prtpow_vmean.eq.prtpow-1 .and. prtpow_vmean.ne.0) then
            ! could be made general, but GCM only uses a factor of 10
            where(vmean.ne.missing) vmean = vmean*10.
            blank4 = '.1* '
          endif
        endif
        title = trim(lname)//' ('//trim(units)//')'

c
c write binary output
c
        if(do_giss) then
          write(lunit) title,jm,km,1,1,
     &         xjl,
     &         lat_dg,pm(1:km),
     &         real(1.,kind=4),real(1.,kind=4),
     &         'LATITUDE        ',
     &         'PRESSURE (MB)   ',
     &         '                ','                ','NASAGISS',
     &         vmean,
     &         xjl_hemis
        endif

c
c print table
c
        where(xjl.eq.missing) xjl=0.
        where(xjl_hemis.eq.missing) xjl_hemis=0.
        where(vmean.eq.missing) vmean=0.
        linect = linect + km + 7
        if(linect.gt.60) then
          write(6,'(a)') '1'//xlabel
          write(6,'(a)') ' '//fromto
          linect = km+8
        endif
        do zone=1,nzones
          j2 = min(maxj, jm -(zone-1)*lats_per_zone)
          j1 = max(minj, j2 -lats_per_zone +1)
          lats_this_zone = 1+(j2-j1)/inc
          write(6,901) title
          call prtdashes(lats_this_zone)
          write(6,fmtlat) int(lat_dg(j2:j1:-inc))
          call prtdashes(lats_this_zone)
          if(has_radonly) then
            do l=3,1,-1
              fsh  = xjl_radonly_hemis(1,l)
              fnh  = xjl_radonly_hemis(2,l)
              fglob= xjl_radonly_hemis(3,l)
              write(6,902) p_radonly(l),fglob,fnh,fsh,
     &             (nint(xjl_radonly(j,l)),j=j2,j1,-inc)
            enddo
          endif
          do l=km,1,-1
            fsh  = xjl_hemis(1,l)
            fnh  = xjl_hemis(2,l)
            fglob= xjl_hemis(3,l)
            write(6,902) pm(l),fglob,fnh,fsh,
     &           (nint(xjl(j,l)),j=j2,j1,-inc)
          enddo
          call prtdashes(lats_this_zone)
          fsh  =vmean(jm+1)
          fnh  =vmean(jm+2)
          fglob=vmean(jm+3)
          write(6,903) blank4,fglob,fnh,fsh,
     &         (nint(vmean(j)),j=j2,j1,-inc)
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(lat_dg,vmean,xjl,xjl_hemis)
      deallocate(lat_budg_dg,lat_agc_dg,lat2_agc_dg)
      deallocate(plm,ple,pm)
      if(allocated(pgz)) deallocate(pgz)

c
c close binary output file
c
      if(do_giss) close(lunit)

      return
  901 FORMAT ('0',30X,A64)
  902 FORMAT (1X,F8.3,3F8.1,1X,32I4)
  903 FORMAT (1X,A6,2X,3F8.1,1X,32I4)
      end subroutine prtajl

      subroutine prtdashes(n)
      implicit none
      integer :: n,nn
      write(6,'(a2)',advance='NO') '  '
      do nn=1,8+n-1
        write(6,'(a4)',advance='NO') '----'
      enddo
      write(6,'(a4)') '----'
      return
      end subroutine prtdashes
