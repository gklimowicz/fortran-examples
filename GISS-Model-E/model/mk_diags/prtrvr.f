      subroutine prtrvr(fid,progargs)
!@sum prtrvr prints diagnostics for named rivers
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:), allocatable :: rvr
      character(len=8), dimension(:), allocatable :: namervr
      integer :: nrvr,n,nn
      character(len=22) :: rvrstr
      character(len=132) :: lineout
      integer :: status
      character(len=132) :: xlabel
      character(len=100) :: fromto
      integer :: idacc(12)
c
c get run ID, time/date info, number of named rivers
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'nrvr',nrvr)

c
c allocate workspace
c
      allocate(rvr(nrvr),namervr(nrvr))

c
c read river names and outflows
c
      call get_var_real(fid,'rvr',rvr)
      call get_var_text(fid,'namervr',namervr)
      call get_var_int(fid,'idacc',idacc)
      rvr = rvr/idacc(1)

c
c print the table
c
      write(6,'(a132)') xlabel
      write(6,'(a32,a100)') '1** River Outflow (km^3/mon) ** ',fromto
      lineout=''
      nn=1
      do n=1,nrvr
        write(rvrstr,'(a9,f8.3,5x)') namervr(n)//':',rvr(n)
        lineout(nn:nn+21)=rvrstr
        nn=nn+22
        if(mod(n,6).eq.0) then
          write(6,'(a132)') lineout
          lineout=''
          nn=1
        endif
      enddo
      if(mod(nrvr,6).ne.0) write(6,'(a132)') lineout

c
c deallocate workspace
c
      deallocate(rvr,namervr)

      return

      end subroutine prtrvr
