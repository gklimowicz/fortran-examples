      subroutine prtadiurn(fid,progargs)
!@sum prtadiurn prints the fields in an input file whose
!@+   metadata mark them as having been created from
!@+   a modelE ADIURN array.
!@+   See conventions.txt for additional info.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(25) :: xh
      integer, dimension(:,:), allocatable :: ijdd
      character(len=80) :: lname
      character(len=4), dimension(:), allocatable :: namdd
      integer :: status,varid
      character(len=132) :: xlabel
      character(len=100) :: fromto
      integer :: ndiupt,kr,kq,ih,nvars
      integer :: dimids(7),dimids_fit(2)

c
c get run ID, time/date info, etc.
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'ndiupt',ndiupt)

c
c allocate workspace
c
      allocate(namdd(ndiupt),ijdd(2,ndiupt))

c
c read adiurn metadata
c
      call get_var_int(fid,'ijdd',ijdd)
      call get_var_text(fid,'namdd',namdd)

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

c
c Loop over locations and quantities.  Quantities with dimensions
c ndiupt,hour are printed.
c
      status = nf_inq_dimid(fid,'ndiupt',dimids_fit(1))
      status = nf_inq_dimid(fid,'hour'  ,dimids_fit(2))
      do kr=1,ndiupt
        write(6,'(a)') '1'//xlabel
        write(6,'(a)') ' '//fromto
        write(6,903) namdd(kr),ijdd(1,kr),ijdd(2,kr),(ih,ih=1,24)
        kq = 0
        do varid=1,nvars
          dimids(1:2) = -1
          status = nf_inq_vardimid(fid,varid,dimids)
          if(any(dimids(1:2).ne.dimids_fit)) cycle
          kq = kq+1
          if(mod(kq-1,5).eq.0) write(6,*)
          status = nf_get_vara_real(fid,varid,(/kr,1/),(/1,24/),xh)
          lname = ''
          status = nf_get_att_text(fid,varid,'long_name',lname)
          where(xh.eq.-1.e30) xh=0.
          xh(25) = sum(xh(1:24))/24.
          write(6,'(A8,25I5)') lname,nint(xh)
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(namdd,ijdd)

      return
  903 FORMAT ('0',A4,I2,',',I2,' ',I2,23I5,'  AVE')
      end subroutine prtadiurn
