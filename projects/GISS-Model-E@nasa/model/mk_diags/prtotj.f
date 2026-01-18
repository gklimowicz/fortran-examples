      subroutine prtotj(fid,progargs)
!@sum prtotj prints ocean northward transports of
!@+   mass, heat, and salt for the Atl, Pac, Ind,
!@+   and global oceans.
!@auth M. Kelley
      implicit none
      include 'netcdf.inc'
      integer :: fid                 ! input file ID
      character(len=160) :: progargs ! options string
      real*4, dimension(:,:), allocatable :: otj
      real*4, dimension(:), allocatable :: lats
      character(len=80) :: lname,units
      integer :: status,varid
      character(len=132) :: xlabel
      character(len=100) :: fromto
      integer :: jmo,j,kq,kc,kb
      character(len=4), dimension(3), parameter ::
     &     qtys=(/'Mass', 'Heat', 'Salt' /)
      character(len=3), dimension(3), parameter ::
     &     circs=(/ 'moc', 'gmf', 'gyr' /)
      character(len=8), dimension(4) :: basin=
     *     (/'Atlantic','Pacific ','Indian  ','Global  ' /)
      character(len=80) :: vname,title

c
c get run ID, time/date info, etc.
c
      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)
      call get_dimsize(fid,'lato2',jmo)

c
c allocate workspace
c
      allocate(otj(jmo,4),lats(jmo))

c
c read otj metadata
c
      call get_var_real(fid,'lato2',lats)

      write(6,'(a)') '1'//xlabel
      write(6,'(a)') ' '//fromto

c
c print total transports by latitude
c
      do kq=1,3
        do kb=1,4
          vname = 'nt_'//qtys(kq)//'_'//basin(kb)(1:3)
          call get_var_real(fid,trim(vname),otj(1,kb))
        enddo
        status = nf_inq_varid(fid,trim(vname),varid)
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        title = 'Northward Transport of '//qtys(kq)//
     &       ' ('//trim(units)//')'
        write(6,*)
        write(6,'(a)') ' '//title
        write(6,'(a)') ' Lat   Atlantic   Pacific    Indian     Global'
        do j=1,jmo
          if(abs(lats(j)).gt.89.99) cycle
          write(6,'(f6.0,4f10.3)') lats(j),otj(j,:)
        enddo
      enddo

c
c print breakdown of transports by circ components
c
      do kq=2,3
        do kb=1,4
          vname = 'nt_'//qtys(kq)//'_'//basin(kb)(1:3)
          call get_var_real(fid,trim(vname),otj(1,1))
          do kc=1,3
            vname = 'nt_'//qtys(kq)//'_'//basin(kb)(1:3)//
     &           '_'//circs(kc)
            call get_var_real(fid,trim(vname),otj(1,kc+1))
          enddo
          status = nf_inq_varid(fid,trim(vname),varid)
          units = ''
          status = nf_get_att_text(fid,varid,'units',units)
          title = 'Northward Transport of '//qtys(kq)//' in '//
     &         trim(basin(kb))//' Ocean '//' ('//trim(units)//')'
          write(6,*)
          write(6,'(a)') ' '//title
          write(6,'(a)')
     &         ' Lat      Total    Overturn   GM_flux   Hor_Gyre'
          do j=1,jmo
            if(abs(lats(j)).gt.89.99) cycle
            write(6,'(f6.0,4f10.3)') lats(j),otj(j,:)
          enddo
        enddo
      enddo

c
c deallocate workspace
c
      deallocate(otj,lats)

      return

      end subroutine prtotj
