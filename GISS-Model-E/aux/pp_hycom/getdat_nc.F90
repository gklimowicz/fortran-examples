      subroutine getdat_nc(flnm,day0,succes)
!
! --- read hybrid model fields (binary) and extract portion of global fields.
! --- then  t r a n s f o r m   to   i s o p y c n i c   fields
!
      use netcdf
      use hycom_dimen
      use hycom_arrays, only: srfhgt,dpmixl,thkice,covice,scp2,    &
          depths,u,v,dp,p,temp,saln,th3d,tracer,uflx,vflx,diaflx , &
          thmix,tmix,smix,umix,vmix,alloc_hycom_arrays
      use const_proc,only: iorign,jorign,grvty,rho,dz,spval,timav, &
          cnvert,itest,jtest,onem

!      use constants,only: cnvert,timav,g,rho,iorign,jorign,spval, &
!                         latij,scvx,scuy,itest,jtest,onem
!      use variables
!      use fields

      implicit none
      character,intent(IN) :: flnm*(*)
      real    ,intent(OUT) :: day0
      logical ,intent(OUT) :: succes
      integer ni,nt,nrec,ncid,id_fld
      logical,parameter :: vrbos=.true.
      real theta(kdm),pnew(idm,jdm,kdm+1)
      real*4 real4(ifull,jfull)
      real real3d(ifull,jfull,kdm),real2d(ifull,jfull),realk4(ifull,jfull,4)
      logical :: check

      integer :: start2d(3)   ! starting indices for the variable in the netcdf file
      integer :: kount2d(3)   ! lengths for the variable in the netcdf file
      integer :: start3d(4)   ! starting indices for the variable in the netcdf file
      integer :: kount3d(4)   ! lengths for the variable in the netcdf file
      real*4  day4
!
      call errhandl (nf90_open (path=trim(flnm), mode=NF90_NOWRITE,ncid=ncid))
      call errhandl (nf90_inq_varid (ncid, 'time', id_fld))
      call errhandl (nf90_get_var  (ncid, id_fld, day4))
      day0=day4
      print *,' day0=',day0
!
      print *,'acquiring theta...'
      call readcdf_1d(ncid,'theta',theta,kdm)
      write(*,'(a,50f6.2)') ' theta=',theta
      write(*,'(a,7i4)') 'iorign,jorign,ifull,jfull,idm,jdm,kdm=',iorign,jorign,ifull,jfull,idm,jdm,kdm
!
      start2d = (/1,1,1/)
      kount2d = (/ifull,jfull,1/)
!
! depths is read in through gtdpth
      if (day0.le.31) then
        print *,'acquiring topo...'
        call readcdf_2d(ncid,'topo',real2d,ifull,jfull,start2d,kount2d,check)
        real4(:,:)=real2d(:,:)
!       call flipp(real4)	! have to be real*4
          do i=1,idm
          do j=1,jdm
            if (abs(real4(i,j)-depths(i,j)) .gt. 0.1) then
              write(*,'(a,2i4,2f8.0)')'mismatched depth at ij,=',i,j,depths(i,j),real4(i,j)
              stop 'read-in depth does not match the archived depth'
            end if
          end do
          end do
          write(*,*)  'chk prescribed depth file matches with archived'
!       call extrct(real4,iorign,jorign,depths)	! input real*4, output: real
      end if
!
      print *,'acquiring scp2...'
      call readcdf_2d(ncid,'scp2',real2d,ifull,jfull,start2d,kount2d,check)
        real4(:,:)=real2d(:,:)
!       call flipp(real4)
        call extrct(real4,iorign,jorign,scp2)
!
      print *,'acquiring srfht...'
      if (timav) then
        call readcdf_2d(ncid,'srfhtav',real2d,ifull,jfull,start2d,kount2d,check)
      else
        call readcdf_2d(ncid,'srfht',real2d,ifull,jfull,start2d,kount2d,check)
      end if
      real4(:,:)=real2d(:,:)
!     call flipp(real4)
      call extrct(real4,iorign,jorign,srfhgt)
!
      print *,'acquiring covice...'
      if (timav) then
        call readcdf_2d(ncid,'coviceav',real2d,ifull,jfull,start2d,kount2d,check)
      else
        call readcdf_2d(ncid,'covice',real2d,ifull,jfull,start2d,kount2d,check)
      end if
      real4(:,:)=real2d(:,:)
!     call flipp(real4)
      call extrct(real4,iorign,jorign,covice)

!     call readcdf_2d(ncid,'eminp',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,eminp)
!     call readcdf_2d(ncid,'srflx',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,srflx)
!     call readcdf_2d(ncid,'radfx',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,radfl)
!     call readcdf_2d(ncid,'snsib',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,snsib)
!     call readcdf_2d(ncid,'latnt',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,latnt)
!     call readcdf_2d(ncid,'zmixl',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,dpmix)
!     call readcdf_2d(ncid,'tice',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,temice)
!     call readcdf_2d(ncid,'zice',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,thkice)
!
! lat/lon are read in through meshsz
!     start3d = (/1,1,1,1/)
!     kount3d = (/ifull,jfull,4,1/)
!
!     call readcdf_3d(ncid,'lat',realk4,ifull,jfull,4,start3d,kount3d,check)
!     do k=1,4
!        real4(:,:)=realk4(:,:,k)
!     call readcdf_2d(ncid,'eminp',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,eminp)
!     call readcdf_2d(ncid,'srflx',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,srflx)
!     call readcdf_2d(ncid,'radfx',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,radfl)
!     call readcdf_2d(ncid,'snsib',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,snsib)
!     call readcdf_2d(ncid,'latnt',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,latnt)
!     call readcdf_2d(ncid,'zmixl',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,dpmix)
!     call readcdf_2d(ncid,'tice',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,temice)
!     call readcdf_2d(ncid,'zice',real2d,ifull,jfull,start2d,kount2d)
!       real4(:,:)=real2d(:,:)
!       call flipp(real4)
!       call extrct(real4,iorign,jorign,thkice)
!
! lat/lon are read in through meshsz
!     start3d = (/1,1,1,1/)
!     kount3d = (/ifull,jfull,4,1/)
!
!     call readcdf_3d(ncid,'lat',realk4,ifull,jfull,4,start3d,kount3d,check)
!     do k=1,4
!        real4(:,:)=realk4(:,:,k)
!        call flipp(real4)
!        call extrct(real4,iorign,jorign,latij(1,1,k))
!     end do
!
!     call readcdf_3d(ncid,'lon',realk4,ifull,jfull,4,start3d,kount3d,check)
!     do k=1,4
!        real4(:,:)=realk4(:,:,k)
!        call flipp(real4)
!        call extrct(real4,iorign,jorign,lonij(1,1,k))
!     end do
!
      print *,'acquiring temp...'
      start3d = (/1,1,1,1/)
      kount3d = (/ifull,jfull,kdm,1/)
      if (timav) then
        call readcdf_3d(ncid,'tempav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      else
        call readcdf_3d(ncid,'temp',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      end if
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipp(real4)
         call extrct(real4,iorign,jorign,temp(1,1,k))
      end do
!
      print *,'acquiring saln...'
      if (timav) then
        call readcdf_3d(ncid,'salnav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      else
        call readcdf_3d(ncid,'saln',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      end if
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipp(real4)
         call extrct(real4,iorign,jorign,saln(1,1,k))
      end do
!
      print *,'acquiring th3d...'
      if (timav) then
        call readcdf_3d(ncid,'th3dav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      else
        call readcdf_3d(ncid,'th3d',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      end if
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipp(real4)
         call extrct(real4,iorign,jorign,th3d(1,1,k))
      end do
!
      print *,'acquiring thik...'
      if (timav) then
        call readcdf_3d(ncid,'thikav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      else
        call readcdf_3d(ncid,'thik',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      end if
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipp(real4)
         call extrct(real4,iorign,jorign,dp(1,1,k))
      end do
!
      print *,'acquiring u......'
      if (timav) then
        call readcdf_3d(ncid,'uav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      else
        call readcdf_3d(ncid,'utotal',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      end if
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipu(real4)
         call extrct(real4,iorign,jorign,u(1,1,k))
      end do
!
      print *,'acquiring v...'
      if (timav) then
        call readcdf_3d(ncid,'vav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      else
        call readcdf_3d(ncid,'vtotal',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      end if
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipv(real4)
         call extrct(real4,iorign,jorign,v(1,1,k))
      end do
!
      print *,'acquiring uflxav...'
      call readcdf_3d(ncid,'uflxav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipu(real4)
         call extrct(real4,iorign,jorign,uflx(1,1,k))
      end do
!
      print *,'acquiring vflxav...'
      call readcdf_3d(ncid,'vflxav',real3d,ifull,jfull,kdm,start3d,kount3d,check)
      do k=1,kdm
         real4(:,:)=real3d(:,:,k)
!        call flipv(real4)
         call extrct(real4,iorign,jorign,vflx(1,1,k))
      end do
!
!     call readcdf_3d(ncid,'diaflx',real3d,ifull,jfull,kdm,start3d,kount3d,check)

      if (ntrcr.ge.1) then
        print *,'acquiring trc1...'
        call readcdf_3d(ncid,'trc1',real3d,ifull,jfull,kdm,start3d,kount3d,check)
        do k=1,kdm
          real4(:,:)=real3d(:,:,k)
!         call flipv(real4)
          call extrct(real4,iorign,jorign,tracer(1,1,k,1))
        end do
      end if

      if (ntrcr.ge.2) then
        print *,'acquiring trc2...'
        call readcdf_3d(ncid,'trc2',real3d,ifull,jfull,kdm,start3d,kount3d,check)
        do k=1,kdm
          real4(:,:)=real3d(:,:,k)
!         call flipv(real4)
          call extrct(real4,iorign,jorign,tracer(1,1,k,2))
        end do
      end if

      if (ntrcr.ge.3) then
        print *,'acquiring trc3...'
        call readcdf_3d(ncid,'trc3',real3d,ifull,jfull,kdm,start3d,kount3d,check)
        do k=1,kdm
          real4(:,:)=real3d(:,:,k)
!         call flipv(real4)
          call extrct(real4,iorign,jorign,tracer(1,1,k,3))
        end do
      end if

      if (ntrcr.ge.4) then
        print *,'acquiring trc4...'
        call readcdf_3d(ncid,'trc4',real3d,ifull,jfull,kdm,start3d,kount3d,check)
        do k=1,kdm
          real4(:,:)=real3d(:,:,k)
!         call flipv(real4)
          call extrct(real4,iorign,jorign,tracer(1,1,k,4))
        end do
      end if

      if (ntrcr.ge.5) stop '(ntrcr cannot > 5)'

!$OMP PARALLEL DO
      do 20 j=1,jdm
!
      do 21 i=1,idm
      p(i,j,1)=0.
      if (ip(i,j).gt.0) then
        srfhgt(i,j)=srfhgt(i,j)*100.*ip(i,j)  !  convert from m to cm
!       eminp(i,j)=eminp(i,j)*365.*86400.     !  eminp m/s -> m/yr
      else
        srfhgt(i,j)=spval
!       dpmix(i,j)=spval
!       sssobs(i,j)=spval
!       rsiobs(i,j)=spval
      end if
 21   continue
!
      do 20 k=1,kdm
      do 20 i=1,idm
!     uflx(i,j,k)=uflx(i,j,k)	! in Sv
!     vflx(i,j,k)=vflx(i,j,k)	! in Sv
      u(i,j,k)=u(i,j,k)*100./dz*iu(i,j)	! convert from m/s to cm/s
      v(i,j,k)=v(i,j,k)*100./dz*iv(i,j)	! convert from m/s to cm/s
      if (ip(i,j).lt.0) then
        dp(i,j,k)=spval
        temp(i,j,k)=spval
        saln(i,j,k)=spval
        th3d(i,j,k)=spval
      end if
 20   continue
!$OMP END PARALLEL DO
!
      call findmx(ip,dp(1,1,2),idm,idm-1,jdm-1,'lyr  2 thk (m)')
      call findmx(ip,dp(1,1,8),idm,idm-1,jdm-1,'lyr  8 thk (m)')
!
!$OMP PARALLEL DO
      do 39 j=1,jdm
!
! --- save mixed layer fields for future use
      do 38 i=1,idm
      thmix(i,j)=th3d(i,j,1)
      tmix(i,j)=temp(i,j,1)
      smix(i,j)=saln(i,j,1)
      umix(i,j)=u(i,j,1)
      vmix(i,j)=v(i,j,1)
      p(i,j,2)=dp(i,j,1)
 38   pnew(i,j,1)=0.
!
      do 39 k=2,kdm
      do 39 i=1,idm
      th3d(i,j,k)=max(th3d(i,j,k-1),th3d(i,j,k))
 39   p(i,j,k+1)=p(i,j,k)+dp(i,j,k)
!$OMP END PARALLEL DO
!
      if (day0.le.31) then
        print *,'shown below: bottom pressure day0=',day0
        call zebra(p(1,1,kdm+1),idm,idm-1,jdm-1)
        write (*,'(a)') 'shown below: ice cover (%)'
        call zebra(covice,idm,idm-1,jdm-1)
      end if
!
! --- transform hybrid fields to isopycnic coord.
!
      if (cnvert) then
        write (*,*) 'now transforming input fields to isopyc.coord.'
        call restep(u,v,temp,saln,tracer,th3d,p,u(1,1,kdm+1),		&
          v(1,1,kdm+1),temp(1,1,kdm+1),saln(1,1,kdm+1),			&
          tracer(1,1,kdm+1,1),th3d(1,1,kdm+1),pnew,theta,kdm,kdm)
! --- converted fields are saved in kdm+k & dpi; raw fields are in kdm & dp
        dp(:,:,:)=0.
!$OMP PARALLEL DO
        do 36 j=1,jdm
        do 36 k=1,kdm
        do 36 i=1,idm
        u(i,j,k)=u(i,j,kdm+k)
        v(i,j,k)=v(i,j,kdm+k)
!       temp(i,j,k)=temp(i,j,kdm+k)
!       saln(i,j,k)=saln(i,j,kdm+k)
!       th3d(i,j,k)=th3d(i,j,kdm+k)
        if (ntrcr.gt.0) tracer(i,j,k,:)=tracer(i,j,kdm+k,:)
!       p(i,j,k+1)=pnew(i,j,k+1)
        if (ip(i,j).gt.0) dp(i,j,k)=pnew(i,j,k+1)-pnew(i,j,k)
 36     continue
!$OMP END PARALLEL DO
      end if                            ! cnvert
!
!     write (*,'(a)') 'shown below: mixed layer density'
!     call zebra(thmix,idm,idm-1,jdm-1)
!     write (*,'(a)') 'shown below: sea surface height (cm)'
!     call zebra(srfhgt,idm,idm-1,jdm-1)
!     write (*,'(a)') 'shown below: SST'
!     call zebra(temp(1,1,2),idm,idm-1,jdm-1)
!     write (*,'(a)') 'shown below: layer 1 thickness (m)'
!     call zebra(dp,idm,idm-1,jdm-1)
!
      succes=.true.
      return
!
 6    succes=.false.
      return
      end

      subroutine errhandl (ret)
      use netcdf
      implicit none

      integer, intent(in) :: ret

      if (ret /= NF90_NOERR) then
         write(6,*) nf90_strerror (ret)
         stop 999
      end if
      return
      end subroutine errhandl

      subroutine readcdf_1d(ncid,varname,field,kdm)
!
! --- extract a 1-D field (real) from a previously opened netcdf file
!
      use netcdf
      implicit none
      integer        ,intent(IN)      :: ncid,kdm
      character*(*),intent(IN)      :: varname
      real           ,intent(OUT)     :: field(kdm)
      integer :: id_fld

      call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
      call errhandl (nf90_get_var (ncid, id_fld, field))
      return
      end subroutine readcdf_1d

      subroutine readcdf_2d(ncid,varname,field,idm,jdm,start,count,check)
!
! --- extract a 2-D field (real) from a previously opened netcdf file
!
      use netcdf
      implicit none
      integer        ,intent(IN)      :: ncid,idm,jdm,start(3),count(3)
      character*(*),intent(IN)      :: varname
      logical,optional,intent(IN)  :: check
      real           ,intent(OUT)     :: field(idm,jdm)
      integer :: iz,jz,tz,id_idm,id_jdm,id_tdm,id_fld

      if (present(check)) then
       if (check) then
! --- verify that dimensions of archived field agree with specified dimensions
        call errhandl (nf90_inq_dimid (ncid, 'idm', id_idm))
        call errhandl (nf90_inq_dimid (ncid, 'jdm', id_jdm))
        call errhandl (nf90_inq_dimid (ncid, 'tdm', id_tdm))

        call errhandl (nf90_inquire_dimension (ncid, id_idm, len=iz))
        call errhandl (nf90_inquire_dimension (ncid, id_jdm, len=jz))
        call errhandl (nf90_inquire_dimension (ncid, id_tdm, len=tz))

        if (iz.ne.idm .or. jz.ne.jdm) then
          print '(2(a,2i5))','(readcdf) array dimensions are',iz,jz,          &
            ', not the epected',idm,jdm
          stop '(readcdf error)'
        end if

        if (start(3)+count(3)-1.gt.tz) then
          print '(2(a,i5))','(readcdf) requested record',                         &
            start(3)+count(3)-1,'  does not exist. records in file:',tz
          stop '(readcdf error)'
        end if
       end if
      end if

      call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
      call errhandl (nf90_get_var (ncid, id_fld, field, start=start, count=count))
      return
      end subroutine readcdf_2d

      subroutine readcdf_3d(ncid,varname,field,idm,jdm,kdm,start,count,check)
!
! --- extract a 3-D field (real) from a previously opened netcdf file
!
      use netcdf
      implicit none
      integer        ,intent(IN)      :: ncid,idm,jdm,kdm,start(4),count(4)
      character*(*),intent(IN)      :: varname
      logical,optional,intent(IN)  :: check
      real           ,intent(OUT)     :: field(idm,jdm,kdm)
      integer :: iz,jz,kz,tz,id_idm,id_jdm,id_kdm,id_tdm,id_fld

      if (present(check)) then
       if (check) then
! --- verify that dimensions of archived field agree with specified dimensions
        call errhandl (nf90_inq_dimid (ncid, 'idm', id_idm))
        call errhandl (nf90_inq_dimid (ncid, 'jdm', id_jdm))
        call errhandl (nf90_inq_dimid (ncid, 'kdm', id_kdm))
        call errhandl (nf90_inq_dimid (ncid, 'tdm',id_tdm))

        call errhandl (nf90_inquire_dimension (ncid, id_idm, len=iz))
        call errhandl (nf90_inquire_dimension (ncid, id_jdm, len=jz))
        call errhandl (nf90_inquire_dimension (ncid, id_kdm, len=kz))
        call errhandl (nf90_inquire_dimension (ncid, id_tdm, len=tz))

        if (iz.ne.idm .or. jz.ne.jdm .or. kz.ne.kdm) then
          print '(2(a,3i5))','(readcdf) array dimensions are',iz,jz,kz,          &
            ', not the epected',idm,jdm,kdm
          stop '(readcdf error)'
        end if

        if (start(4)+count(4)-1.gt.tz) then
          print '(2(a,i5))','(readcdf) requested record',                         &
            start(4)+count(4)-1,'  does not exist. records in file:',tz
          stop '(readcdf error)'
        end if
       end if
      end if

      call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
      call errhandl (nf90_get_var (ncid, id_fld, field, start=start, count=count))
      return
      end subroutine readcdf_3d
