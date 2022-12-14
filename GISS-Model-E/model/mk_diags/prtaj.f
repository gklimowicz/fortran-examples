      subroutine prtaj(fid,progargs)
!@sum prtaj prints the fields in an input file whose
!@+   metadata mark them as having been created
!@+   from a modelE AJ array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:), allocatable :: lat_dg,xj
      real*4, dimension(3) :: xj_hemis
      character(len=16), dimension(:), allocatable :: terrain
      character(len=16) :: stitle
      character(len=30) :: fmt
      character(len=40) :: vname,vname_hemis
      real*4 :: fglob,fnh,fsh
      integer :: j,jm,inc,nt,ntype,lstr
      integer :: srt(7),cnt(7)
c
c Various string formats
c
      character(len=80), parameter ::
     &     fmtlat = "('0',131('-')/20X,'G      NH     SH   ',24I4)"

      integer :: status,varid,varid_hemis,nvars
      character(len=132) :: xlabel
      character(len=100) :: fromto

c
c get run ID, time/date info, number of latitudes and number of surface types
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'lat_budg',jm)
      call get_dimsize(fid,'ntype',ntype)

c
c allocate workspace
c
      allocate(lat_dg(jm),xj(jm))
      allocate(terrain(ntype))

c
c read geometry and surface type names
c
      call get_var_real(fid,'lat_budg',lat_dg)
      call get_var_text(fid,'terrain',terrain)

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c Loop over surface types and quantities.  An array xyz
c is printed if there also exists an array xyz_hemis
c containing hemispheric/global averages for that quantity.
c
      inc=1+(jm-1)/24

      do nt=1,ntype
        write(6,'(a)') '1'//xlabel
        write(6,'(a12,a16,a3,a100)')
     &       '0** BUDGETS ',terrain(nt),'** ',fromto
        write(6,fmtlat) nint(lat_dg(jm:inc:-inc))
        write(6,905)
        do varid_hemis=1,nvars
          status = nf_inq_varname(fid,varid_hemis,vname_hemis)
          lstr = len_trim(vname_hemis)
          if(vname_hemis(lstr-5:lstr).ne.'_hemis') cycle
          vname = vname_hemis(1:lstr-6)
          status = nf_inq_varid(fid,trim(vname),varid)
          stitle = ''
          status = nf_get_att_text(fid,varid,'stitle',stitle)
          if(status.ne.nf_noerr) cycle
          if(trim(stitle).eq.'no output') cycle
          fmt = ''
          status = nf_get_att_text(fid,varid,'fmt',fmt)
          if(status.ne.nf_noerr) cycle
          srt(1:3) = (/ 1,nt,1 /)
          cnt(1:3) = (/ 3,1,1 /)
          status = nf_get_vara_real(fid,varid_hemis,srt,cnt,xj_hemis)
          cnt(1:3) = (/ jm,1,1 /)
          status = nf_get_vara_real(fid,varid,srt,cnt,xj)
          where(xj.eq.-1.e30) xj=0.
          where(xj_hemis.eq.-1.e30) xj_hemis=0.
          fsh  = xj_hemis(1)
          fnh  = xj_hemis(2)
          fglob= xj_hemis(3)
          if(index(fmt,'24I').gt.0) then ! integer format
            write(6,trim(fmt)) stitle,fglob,fnh,fsh,
     &           (nint(xj(j)),j=jm,inc,-inc)
          else
            write(6,trim(fmt)) stitle,fglob,fnh,fsh,
     &           (xj(j),j=jm,inc,-inc)
          endif
        enddo
        write(6,fmtlat) nint(lat_dg(jm:inc:-inc))
        write(6,905)
      enddo

c
c deallocate workspace
c
      deallocate(lat_dg,xj,terrain)

      return
  905 FORMAT (1X,131('-'))
      end subroutine prtaj
