      subroutine prtareg(fid,progargs)
!@sum prtareg prints the fields in an input file whose
!@+   metadata mark them as having been created from
!@+   a modelE AREG array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:), allocatable :: areg
      character(len=4), dimension(2,23) :: namreg
      character(len=16) :: stitle
      character(len=30) :: fmt
      integer :: nreg
      logical :: has_time
      integer :: status,varid,nvars,nreg_dimid,time_dimid,dimids(7)
      character(len=132) :: xlabel
      character(len=100) :: fromto

c
c Various string formats
c
      character(len=80), parameter ::
     &    fmtreg = "('0',16X,23(1X,A4)/17X,23(1X,A4)/1X,131('-'))"

c
c get run ID, time/date info, number of regions
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'nreg',nreg)
      status = nf_inq_dimid(fid,'nreg',nreg_dimid)
      status = nf_inq_unlimdim(fid,time_dimid)
      has_time = status==nf_noerr

c
c allocate workspace
c
      allocate(areg(nreg))

c
c read region names
c
      call get_var_text(fid,'namreg',namreg)

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c Loop over quantities.  An array xyz is printed if its
c sole dimension name is nreg.
c
      write(6,'(a)') '1'//xlabel
      write(6,'(a16,a100)') '   (REGIONS)    ',fromto
      write(6,fmt=fmtreg) reshape((/namreg(1,1:23),namreg(2,1:23)/),
     &     (/23*2/) )
      do varid=1,nvars
        dimids(1:2) = -1
        status = nf_inq_vardimid(fid,varid,dimids)
        if(dimids(1).ne.nreg_dimid) then
          if(has_time) then
            if(dimids(2).ne.time_dimid) cycle
          elseif(dimids(2).ne.-1) then
            cycle
          endif
        endif
        stitle = ''
        status = nf_get_att_text(fid,varid,'stitle',stitle)
        if(status.ne.nf_noerr) cycle
        if(trim(stitle).eq.'no output') cycle
        fmt = ''
        status = nf_get_att_text(fid,varid,'fmt',fmt)
        if(trim(fmt).eq.'not computed') cycle
        status = nf_get_var_real(fid,varid,areg)
        where(areg.eq.-1.e30) areg=0.
        if(index(fmt,'23I').gt.0) then ! integer format
          write(6,trim(fmt)) stitle,nint(areg(1:23))
        else
          write(6,trim(fmt)) stitle,areg(1:23)
        endif
      enddo
      write(6,fmt=fmtreg) reshape((/namreg(1,1:23),namreg(2,1:23)/),
     &     (/23*2/) )
      write(6,905)

c
c deallocate workspace
c
      deallocate(areg)

      return
  905 FORMAT (1X,131('-'))
      end subroutine prtareg
