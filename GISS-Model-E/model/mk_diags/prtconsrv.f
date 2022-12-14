      subroutine prtconsrv(fid,progargs)
!@sum prtconsrv prints the fields in an input file whose
!@+   metadata mark them as having been created from
!@+   a modelE CONSRV-type array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:), allocatable :: lat_dg,dxyp,xj
      integer, dimension(:), allocatable :: marea
      logical, dimension(:), allocatable :: is_inst_qty,is_sumchg_qty
      integer, dimension(:), allocatable :: hemis_ids,var_ids,n1,n2
      real*4, dimension(3) :: xj_hemis
      character(len=38) :: title
      character(len=40) :: vname,vname_hemis
      real*4 :: fglob,fhem,aglob
      integer :: j,jm,inc,lstr,jhemi,jp1,jpm,jeq,n,nqty,iqty,nconsrv
      character*4, parameter :: hemis(2) = (/' SH ',' NH '/),
     *     dash = ('----')

      integer :: status,varid,varid_hemis,nvars
      character(len=132) :: xlabel
      character(len=100) :: fromto

c
c get run ID, time/date info, number of latitudes
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'lat_budg',jm)

c
c allocate workspace
c
      allocate(lat_dg(jm),xj(jm),dxyp(jm),marea(jm))

c
c read geometry
c
      call get_var_real(fid,'lat_budg',lat_dg)
      call get_var_real(fid,'area_budg',dxyp)
      marea = nint(dxyp*1.e-10)
      aglob = sum(dxyp)*1.e-10

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c scan for conservation quantities
c
      allocate(hemis_ids(nvars),var_ids(nvars))
      allocate(is_inst_qty(nvars),is_sumchg_qty(nvars))
      is_inst_qty(:) = .false.
      n = 0
      do varid_hemis=1,nvars
        status = nf_inq_varname(fid,varid_hemis,vname_hemis)
        lstr = len_trim(vname_hemis)
        if(vname_hemis(lstr-5:lstr).ne.'_hemis') cycle
        vname = vname_hemis(1:lstr-6)
        status = nf_inq_varid(fid,trim(vname),varid)
        title = ''
        status = nf_get_att_text(fid,varid,'long_name',title)
        n = n + 1
        hemis_ids(n) = varid_hemis
        var_ids(n) = varid
        if(index(title,'INSTANT').gt.0) is_inst_qty(n) = .true.
        if(index(title,'SUM OF CHANGES').gt.0) is_sumchg_qty(n) = .true.
      enddo
      nconsrv = n

c
c sanity checks
c
      if(count(is_inst_qty) .ne. count(is_sumchg_qty)) then
        write(6,*) 'parse error: the number of INSTANT fields does '//
     &       'not match the number of'
        write(6,*) 'SUM OF CHANGES fields'
        stop
      endif
      do n=2,nconsrv
        if(is_inst_qty(n) .and. .not. is_sumchg_qty(n-1)) then
          write(6,*) 'The INSTANT field for species N must follow '//
     &         'the SUM OF CHANGES field'
          write(6,*) 'for species N-1.  '//
     &         'Have fields been reordered after scaleacc?'
          stop
        endif
      enddo

      nqty = count(is_inst_qty)
      allocate(n1(nqty),n2(nqty))
      nqty = 0
      do n=1,nconsrv
        if(is_inst_qty(n)) then
          nqty = nqty + 1
          n1(nqty) = n
        endif
        if(is_sumchg_qty(n)) n2(nqty) = n
      enddo

c
c Loop over quantities and hemispheres
c
      inc=1+(jm-1)/24
      jeq = 1+jm/2
      do iqty=1,nqty
        write(6,'(a)') '1'//trim(xlabel)
        write(6,'(a31,a100)') '0Conservation Quantities       ',fromto
        do jhemi=2,1,-1
          jp1=1+(jhemi-1)*(jeq-1)
          jpm=jhemi*(jeq-1)
          write(6,903) (dash,j=jp1,jpm,inc)
          write(6,904) hemis(jhemi),(nint(lat_dg(j)),j=jpm,jp1,-inc)
          write(6,903) (dash,j=jp1,jpm,inc)
          do n=n1(iqty),n2(iqty)
            varid_hemis = hemis_ids(n)
            varid = var_ids(n)
            title = ''
            status = nf_get_att_text(fid,varid,'long_name',title)
            status = nf_get_var_real(fid,varid_hemis,xj_hemis)
            status = nf_get_var_real(fid,varid,xj)
            where(xj.eq.-1.e30) xj=0.
            where(xj_hemis.eq.-1.e30) xj_hemis=0.
            fhem = xj_hemis(jhemi)
            fglob= xj_hemis(3)
            write(6,905) title,fglob,fhem,(nint(xj(j)),j=jpm,jp1,-inc)
          enddo
          write (6,906) aglob,.5*aglob,(marea(j),j=jpm,jp1,-inc)
        enddo
        write(6,*) ''
      enddo
c
c deallocate workspace
c
      deallocate(lat_dg,xj,dxyp,marea)
      deallocate(hemis_ids,var_ids,is_inst_qty,is_sumchg_qty,n1,n2)

      return
  901 FORMAT ('1',A)
  903 FORMAT (1X,28('--'),13(A4,'--'))
  904 FORMAT (41X,'GLOBAL',A7,2X,13I6)
  905 FORMAT (A38,2F9.2,1X,13I6)
  906 FORMAT ('0AREA (10**10 m^2)',20X,2F9.1,1X,13I6)
      end subroutine prtconsrv
