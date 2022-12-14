      subroutine prtolnst(fid,progargs)
!@sum prtolnst prints the fields in an input file whose
!@+   metadata mark them as having been created from
!@+   a modelE OLNST array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      integer, dimension(:), allocatable :: lmst
      character(len=20), dimension(:), allocatable :: name_st
      real*4, dimension(:,:), allocatable :: as
      character(len=80) :: lname,units,title
      integer :: status,varid
      character(len=132) :: xlabel
      character(len=100) :: fromto
      integer :: nmst,lmo,nvars,n,l,lmax,line_length
      integer :: dimids(7),did_zoc,did_zoce,did_nmst
      logical :: print_this
      real*4 :: vsum
      character(len=160) :: dashes,header

      do n=1,len(dashes)
        dashes(n:n) = '-'
      enddo
      header=' Strait               Sum/Mean'//
     &     '     1     2     3     4     5     6     7     8'//
     &     '     9    10    11    12    13    14    15    16'

c
c get run ID, time/date info, etc.
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'zoc',lmo)
      call get_dimsize(fid,'nmst',nmst)

c
c allocate workspace
c
      allocate(as(lmo,nmst),name_st(nmst),lmst(nmst))

c
c read olnst metadata
c
      call get_var_int(fid,'lmst',lmst)
      call get_var_text(fid,'strait_name',name_st)
      line_length = 30 + 6*maxval(lmst)

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c Loop over locations and quantities.  Quantities with dimensions
c zoc,nmst are printed.
c
      status = nf_inq_dimid(fid,'zoc',did_zoc)
      did_zoce = did_zoc
      status = nf_inq_dimid(fid,'zoce',did_zoce)
      status = nf_inq_dimid(fid,'nmst',did_nmst)
      write(6,'(a)') '1'//xlabel
      write(6,'(a)') ' '//fromto
      do varid=1,nvars
        print_this = .true.
        dimids(1:2) = -1
        status = nf_inq_vardimid(fid,varid,dimids)
        if(dimids(2).ne.did_nmst) print_this=.false.
        if(dimids(1).ne.did_zoc .and. dimids(1).ne.did_zoce)
     &       print_this=.false.
        if(.not.print_this) cycle
        status = nf_get_var_real(fid,varid,as)
        lname = ''
        status = nf_get_att_text(fid,varid,'long_name',lname)
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        write(6,*)
        write(6,'(a)') ' '//trim(lname)//' ('//trim(units)//')'
        write(6,'(a)') dashes(1:line_length)
        write(6,'(a)') header(1:line_length)
        write(6,'(a)') dashes(1:line_length)
        do n=1,nmst
          lmax = lmst(n)
          vsum = sum(as(1:lmax,n))
          if(index(lname,'Trans').eq.0) vsum=vsum/lmax
          write(6,'(1x,a20,f8.1,1x,16i6)') name_st(n),vsum,
     &         nint(as(1:lmax,n))
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(as,name_st,lmst)

      return
      end subroutine prtolnst
