      subroutine prtisccp(fid,progargs)
!@sum prtisccp prints ISCCP histograms
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:,:,:), allocatable :: aisccp
      real*4, dimension(:), allocatable ::
     &     wisccp,isccp_tau,isccp_late
      integer, dimension(:), allocatable :: isccp_press
      integer :: npres,ntau,nisccp
      integer :: status
      character(len=132) :: xlabel
      character(len=100) :: fromto
      character(len=80) :: title
      character(len=3) :: degs(100)
      integer :: n,ipress
c
c get run ID, time/date info, dimension sizes
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'nisccp',nisccp)
      call get_dimsize(fid,'ntau',ntau)
      call get_dimsize(fid,'npres',npres)

c
c allocate workspace
c
      allocate(aisccp(ntau,npres,nisccp),wisccp(nisccp))
      allocate(isccp_press(npres),isccp_tau(ntau),isccp_late(nisccp+1))

c
c read data
c
      call get_var_real(fid,'aisccp',aisccp)
      call get_var_real(fid,'wisccp',wisccp)
      call get_var_real(fid,'isccp_tau',isccp_tau) ! not used in titles
      call get_var_int(fid,'isccp_press',isccp_press)
      call get_var_real(fid,'isccp_late',isccp_late)

c
c print the tables
c
      do n=1,nisccp+1
        write(degs(n),'(i2)') int(abs(isccp_late(n)))
        if(isccp_late(n).gt.0.) then
          degs(n)(3:3)='N'
        else
          degs(n)(3:3)='S'
        endif
      enddo
      write(6,'(a132)') xlabel
      write(6,'(a100)') fromto
      do n=nisccp,1,-1
        title='ISCCP CLOUD FREQUENCY (NTAU,NPRES) % '//
     &       degs(n)//'-'//degs(n+1)
        write(6,100) title
        do ipress=1,npres
          write(6,101) isccp_press(ipress),
     &         100.*aisccp(2:ntau,ipress,n)/wisccp(n)
        enddo
        write(6,*)
      enddo

c
c deallocate workspace
c
      deallocate(aisccp,wisccp,isccp_press,isccp_tau,isccp_late)

      return
 100  FORMAT (1X,A80/1X,72('-')/3X,
     *     'PRESS\TAU    0.  1.3  3.6  9.4  23   60   > ')
 101  FORMAT (5X,I3,7X,6F5.1)
      end subroutine prtisccp
