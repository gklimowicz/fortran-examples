      program latlonz3d
c     
c     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c     (1) read in monthly output files (*.out) from hycom
c     average them
c     (2) convert selected (11) fields to lat/lon/z grid of 1x1x33
c     in giss or netcdf format
c     <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
c     
      use const_proc
      use hycom_arrays, only: depths,srfhgt,dpmixl,covice,u,v,dp,p
     .     ,temp,saln,th3d,tracer,alloc_hycom_arrays
      use hycom_dimen
      use hycom_o2a

      implicit none
c     
      real :: day0,day1
      integer, allocatable :: im(:,:)
c     
      integer mo,mon1,i70,i45,ieq,status,ia,k00,ny,m,mm,nt
      integer :: nrec
      real*8  :: area,avgbot
      character flnmin*128,ttl*80,ttl1*80,ttl2*80
     .     ,flnmann*80,flnmout*128
      logical :: succes
      character*26 ttlt(ntrcr)
      character*26,dimension(8),parameter::
     .     ttl0=(/'sea surface height (cm)   '
     .     ,'mixed layer depth (m)     '
     .     ,'sea ice coverage [0:1]    '
     .     ,'eastward  velocity (cm/s) '
     .     ,'northward velocity (cm/s) '
     .     ,'temperature (deg C)       '
     .     ,'salinity (psu)            '
     .     ,'density (sigma2)          '/)
      real :: avg2d(iia,jja),avg3d(iia,jja,k33)
     .     ,sshij(iia,jja),dpmixij(iia,jja),iceij(iia,jja)
     .     ,tz(iia,jja,k33),sz(iia,jja,k33),uz(iia,jja,k33)
     .     ,vz(iia,jja,k33),rz(iia,jja,k33),trz(iia,jja,k33,ntrcr)
     .     ,a2d(iia,jja,12),a3d(iia,jja,k33,12)
     .     ,worka(iia,jja),worko(idm,jdm),depthij(iia,jja)
     .     ,srfhgt_av(idm,jdm),dpmixl_av(idm,jdm)
     .     ,oice_av(idm,jdm),         p_av(idm,jdm,kdm+1)
     .     ,   u_av(idm,jdm,kdm),     v_av(idm,jdm,kdm)
     .     ,temp_av(idm,jdm,kdm),  saln_av(idm,jdm,kdm)
     .     ,th3d_av(idm,jdm,kdm),tracer_av(idm,jdm,kdm,ntrcr)
      real*4 :: real4a(iia,jja)
c     
      character(len=1024) :: cmnd_arg   
      character(len=  80) :: cmnd_title 
      integer :: num_args, cmnd_len, cmnd_status, iarg
      real :: factor
      logical :: need_netcdf = .false.
      logical :: need_giss   = .false.

c     
      call get_command_argument (0, cmnd_arg, cmnd_len, cmnd_status)
      if (cmnd_status /= 0) then
         write (*,*) 'Getting command name failed with cmnd_status = ',
     .        cmnd_status
         stop
      end if
      write (unit=*,fmt='(/a,3x,a)') 'command name =', 
     .     cmnd_arg (1:cmnd_len)

      num_args = command_argument_count ()
      write (unit=*,fmt='(a,i3)') 'number of command arguments = ', 
     .     num_args

c-----------------------------------------------------------------------------
      iarg=1
      call get_command_argument (iarg, cmnd_arg, cmnd_len, cmnd_status)
      if (cmnd_status /= 0) then
         write (unit=*,fmt='(a,i2,a)') 
     .        'get_command_argument failed: cmnd_status = ', 
     .        cmnd_status,  ' arg = ', iarg
         stop
      end if

      flnmout=trim(cmnd_arg)
      write (unit=*,fmt='(/a,i2,2a)') 'command arg=', iarg, 
     .     ' (name output file) = ', trim(flnmout)

      need_netcdf = is_suffix_nc( trim(flnmout) ) 

      if( need_netcdf ) then
         write (unit=*,fmt='(/a)') 
     .        "OUTPUT file is NETCDF format !!!"
      else
         need_giss = .true.
         open(40,file=trim(flnmout),form='unformatted',status='unknown')
         write (unit=*,fmt='(a)') 
     .        'OUTPUT file is open for writing GISS format '
      end if

c-----------------------------------------------------------------------------
      iarg=2
      call get_command_argument (iarg, cmnd_arg, cmnd_len, cmnd_status)
      if (cmnd_status /= 0) then
         write (*,*) 'get_command_argument failed: cmnd_status = ', 
     .        cmnd_status,  ' arg = ', iarg
         stop
      end if

      cmnd_title=trim(cmnd_arg)
      write (unit=*,fmt='(/a,i2,3a)') 'command arg=', iarg, 
     .     ' (title output file) ="',trim(cmnd_title),'"'
c     
      write(*,'(/a,i2)') 'number of tracers =',ntrcr
      nrec=8+ntrcr

      do nt=1,ntrcr
         write(ttlt(nt),'(a,i2.2)') 'tracer No.',nt
      enddo

      write(ttl1,'(i3,a,i3,3x)')  iia,'x',jja
      write(ttl2,'(2(i3,a),i2)') iia,'x',jja,'x',k33
c     Check the length of output title. Should be lesss 80
c     write(ttl,'(a,3x,a,3x,a)') ttl0(1),trim(ttl1),trim(cmnd_title)
      if( len(ttl0(1))+3+len(trim(ttl2))+3+ 
     .    len(trim(cmnd_title)) > 80 ) then
         write(*,*) "The length of the title at the output file "
         write(*,*) "    exceed 80 characters"
         write(*,*) "You should decrease length title from comand line!"
         write(*,*) "Safe length less than 38 caracters !!!"
C     write(*,*) " STOP at file=",__FILE__," line=", __LINE__  
      end if

      call alloc_hycom_arrays
      call alloc_hycom_dimen

      allocate (im(idm,jdm),stat=status)
c     --- determine do-loop limits for u,v,p,q points
      call gtdpth(depths,im)
      call bigrid(depths)
      if (diag) then
         print *,' depth at i,j=',itest,jtest
         do i=itest-5,itest+5
            write(*,'(11f6.0)') (depths(i,j),j=jtest-5,jtest+5)
         enddo
      endif
c     
      call o2a_wgt
      call o2a_sfc(depths,depthij)
      if (diag) then
         i=iatest
         j=jatest
         print *,'topo after o2a',i,j,depths(itest,jtest),depthij(i,j)
      endif
c     
      kij=0
      do 30 i=1,iia
      do 30 j=1,jja
         if (depthij(i,j) <= 0.) goto 32
            do 31 k=2,k33
               if(z33(k-1) <= depthij(i,j) .and. 
     &            z33(k)   >  depthij(i,j)) then
                  kij(i,j)=k-1
                  go to 32
               elseif (z33(k33) <= depthij(i,j)) then
                  kij(i,j)=k33
                  go to 32
               endif
 31         continue
 32      continue
 30   continue

      if (diag) then
         write(*,'(21f5.0)') ((depthij(i,j),i=iatest-10,iatest+10)
     .           ,j=jatest-10,jatest+10)
         print *,' kij main'
         write(*,'(21i5)') ((kij(i,j),i=iatest-10,iatest+10)
     .           ,j=jatest-10,jatest+10)
         i=iatest
         j=jatest
         write(*,'(a,3i4,f7.1)') 'chk kij=',i,j,kij(i,j),depthij(i,j)
      end if

      srfhgt_av=0.
      dpmixl_av=0.
      oice_av=0.
      p_av=0.
      u_av=0.
      v_av=0.
      temp_av=0.
      saln_av=0.
      th3d_av=0.
      tracer_av=0.
c     
      factor = 1./ (num_args-2)
      do iarg = 3, num_args
         call get_command_argument (iarg, cmnd_arg, cmnd_len, 
     &        cmnd_status)
         if (cmnd_status /= 0) then
            write (*,*) 'get_command_argument failed: cmnd_status = ', 
     .           cmnd_status,  ' arg = ', iarg
            stop
         end if
         flnmin=trim(cmnd_arg)
         write (unit=*,fmt='(/a,i3.3,2a)') 'command arg=', iarg, 
     .        ' (name input(reading)  file) = ', trim(flnmin)

c     --- read archive data
         timav=.true.
         cnvert=.false.
         print *,'read data from ',flnmin
         call getdat(flnmin,day0,day1,succes)
c     
         do  i=1,idm
            do  j=1,jdm
            srfhgt_av(i,j)=srfhgt_av(i,j)+srfhgt(i,j)  *factor
            dpmixl_av(i,j)=dpmixl_av(i,j)+dpmixl(i,j,1)*factor
            oice_av(i,j)=  oice_av(i,j)+covice(i,j)  *factor
            do k=1,kdm
               p_av(i,j,k+1)= p_av(i,j,k+1)+ p(i,j,k+1)*factor
               u_av(i,j,k)=   u_av(i,j,k)+   u(i,j,k)*factor
               v_av(i,j,k)=   v_av(i,j,k)+   v(i,j,k)*factor
               temp_av(i,j,k)=temp_av(i,j,k)+temp(i,j,k)*factor
               saln_av(i,j,k)=saln_av(i,j,k)+saln(i,j,k)*factor
               th3d_av(i,j,k)=th3d_av(i,j,k)+th3d(i,j,k)*factor

               do nt=1,ntrcr
                  tracer_av(i,j,k,nt)=tracer_av(i,j,k,nt)+ 
     &                             tracer   (i,j,k,nt)*factor
               enddo      ! loop nt

            enddo      ! loop k

            end do        ! loop j 
         end do           ! loop i
      end do              ! loop arguments 3, .
c     
      call o2a_sfc(srfhgt_av,sshij)
      if (diag) then
         i=iatest
         j=jatest
         print *,'ssh after o2a ',i,j,srfhgt_av(itest,jtest),sshij(i,j)
      endif

      call o2a_sfc(dpmixl_av,dpmixij)
      if (diag) then
         i=iatest
         j=jatest
         print *,'ml aft o2a ',i,j,dpmixl_av(itest,jtest),dpmixij(i,j)
c     call prtmsk(kij,dpmixij,worka,iia,iia,jja,0.,1.,'dpmxl_ij')
      endif

      call o2a_sfc(oice_av,iceij)
      call o2a_3dvec(p_av,u_av,v_av,uz,vz)
      call o2a_3d(p_av,temp_av,tz)
      call o2a_3d(p_av,saln_av,sz)
      call o2a_3d(p_av,th3d_av,rz)

      do nt=1,ntrcr
         call o2a_3d(p_av,tracer_av(1,1,1,nt),trz(1,1,1,nt))
      enddo

CCC   TEST start
c     do i = 1, iia
c        do j = 1, jja
c           do k = 1, k33
c           do nt=1,ntrcr
c              trz(i,j,k,nt) = i + j + k
c           end do
c           tz(i,j,k) = i+j+k
c           end do
c        end do
c     end do
CCC   TEST end

c     --- Shift the first meridian to DATE line ( -179.5 degrees_east )
      sshij = cshift(sshij,iia/2,1)
      dpmixij = cshift(dpmixij,iia/2,1)
      iceij = cshift(iceij,iia/2,1)
      uz = cshift(uz,iia/2,1)
      vz = cshift(vz,iia/2,1)
      tz = cshift(tz,iia/2,1)
      sz = cshift(sz,iia/2,1)
      rz = cshift(rz,iia/2,1)

      do nt=1,ntrcr
         trz(:,:,:,nt) = cshift(trz(:,:,:,nt),iia/2,1)
      enddo

      if (diag) then
c     call prtmsk(kij,temp,worka,iia,iia,jja,0.,1.,'sst')
         i=iatest
         j=jatest
         print '(2i4,f7.1,3x,4(a,9x),a)',i,j,depthij(i,j),
     .        ' t ',' s ',' u ',' v ','trc'
         do k=1,k33
            write(*,'(i2,f7.0,5f12.4)')k,z33(k),tz(i,j,k),sz(i,j,k)
     .           ,uz(i,j,k),vz(i,j,k),trz(i,j,k,1)
         enddo

         i=iatest
         j=jatest
         print *,'temp  at (i-1,j),(i,j),(i+1),(i,j-1),(i,j+1)',i,j
         do k=1,k33
            write(*,'(i2,5f12.4)') k,tz(i-1,j,k),tz(i,j,k),tz(i+1,j,k)
     .           ,tz(i,j-1,k),tz(i,j+1,k)
         enddo

         do k=1,k33
            print *,'  final  temp ', iatest,jatest,k
            do j=jatest+3,jatest-3,-1
            write(*,'(11f7.2)') (tz(i,j,k),i=iatest-3,iatest+3)
            enddo
         enddo
      endif            !  diag

      if( need_giss ) then
         write(ttl,'(a,3x,a,3x,a)') ttl0(1),trim(ttl1),trim(cmnd_title)
         write(40) ttl, sshij
         write(ttl,'(a,3x,a,3x,a)') ttl0(2),trim(ttl1),trim(cmnd_title)
         write(40) ttl, dpmixij
         write(ttl,'(a,3x,a,3x,a)') ttl0(3),trim(ttl1),trim(cmnd_title)
         write(40) ttl, iceij
c     
         write(ttl,'(a,3x,a,3x,a)') ttl0(4),trim(ttl2),trim(cmnd_title)
         write(40) ttl, uz
         write(ttl,'(a,3x,a,3x,a)') ttl0(5),trim(ttl2),trim(cmnd_title)
         write(40) ttl, vz
         write(ttl,'(a,3x,a,3x,a)') ttl0(6),trim(ttl2),trim(cmnd_title)
         write(40) ttl, tz
         write(ttl,'(a,3x,a,3x,a)') ttl0(7),trim(ttl2),trim(cmnd_title)
         write(40) ttl, sz
         write(ttl,'(a,3x,a,3x,a)') ttl0(8),trim(ttl2),trim(cmnd_title)
         write(40) ttl, rz

         do nt=1,ntrcr
            write(ttl,'(a,3x,a,3x,a)') ttlt(nt), 
     &           trim(ttl2),trim(cmnd_title)
            write(40) ttl,trz(:,:,:,nt)
         enddo

         close(40)
      end if 

      if( need_netcdf ) call write_netcdf
c     
C     write(unit=*,fmt='(/,3a)') " +++ END PROGRAM ", __FILE__," +++ "

      CONTAINS 
      
      subroutine write_netcdf

      implicit none
      include 'netcdf.inc'

      integer, parameter :: nvar=12
      integer :: rc,fid,vid(nvar),nd(nvar),dimids(3)
      real*4 :: lons(360),lats(180)
      real*4, parameter :: missing=-999.
      character(len=20), dimension(nvar) :: units,sname
      character(len=80), dimension(nvar) :: lname
c     c
c     c write netcdf format
c     c
      nd(1:3) = 1
      sname(1)='lon'; units(1)='degrees_east'; lname(1)='longitude'
      sname(2)='lat'; units(2)='degrees_north'; lname(2)='latitude'
      sname(3)='z'; units(3)='m'; lname(3)='depth'

      nd(4:6) = 2
      sname(4)='SSH'; units(4)='cm'; lname(4)='sea surface height'
      sname(5)='Zmix'; units(5)='m'; lname(5)='mixed layer depth'
      sname(6)='Icefr'; units(6)='0:1'; lname(6)='sea ice coverage'

      nd(7:12) = 3
      sname(7)='U'; units(7)='cm/s'; lname(7)='eastward velocity'
      sname(8)='V'; units(8)='cm/s'; lname(8)='northward velocity'
      sname(9)='Temp'; units(9)='C'; lname(9)='temperature'
      sname(10)='Sal'; units(10)='psu'; lname(10)='salinity'
      sname(11)='Dens'; units(11)='kg/m3-1000';
      lname(11)='density (sigma2)'
      sname(12)='Trac1'; units(12)='tr1_unit'; lname(12)='Tracer1'
c     
      do m=1,iia
         lons(m) = -180. + (m-.5)*(360./real(iia))
      enddo
      do m=1,jja
         lats(m) = -90. + (m-.5)*(180./real(jja))
      enddo

      rc = nf_create(trim(flnmout),nf_clobber,fid)
      rc = NF_PUT_ATT_TEXT (fid, NF_GLOBAL, 'title',
     &        len(trim(cmnd_title)), trim(cmnd_title))

      rc = nf_def_dim(fid,'lon',iia,dimids(1))
      rc = nf_def_dim(fid,'lat',jja,dimids(2))
      rc = nf_def_dim(fid,'z',k33,dimids(3))
      rc = nf_def_var(fid,'lon',nf_float,1,dimids(1),vid(1))
      rc = nf_def_var(fid,'lat',nf_float,1,dimids(2),vid(2))
      rc = nf_def_var(fid,'z',nf_float,1,dimids(3),vid(3))

      do m=4,12
         rc = nf_def_var(fid,trim(sname(m)),nf_float,nd(m),
     &        dimids,vid(m))
         rc = nf_put_att_real(fid,vid(m),'missing_value',
     &        nf_float,1,missing)
      enddo

      do m=1,12
         rc = nf_put_att_text(fid,vid(m),'units',
     &        len_trim(units(m)),trim(units(m)))
         rc = nf_put_att_text(fid,vid(m),'long_name',
     &        len_trim(lname(m)),trim(lname(m)))
      enddo

      rc = nf_enddef(fid)

      rc = nf_put_var_real(fid,vid(1),lons)
      rc = nf_put_var_real(fid,vid(2),lats)
      rc = nf_put_var_real(fid,vid(3),z33)
      rc = nf_put_var_real(fid,vid(4),sshij)
      rc = nf_put_var_real(fid,vid(5),dpmixij)
      rc = nf_put_var_real(fid,vid(6),iceij)
      rc = nf_put_var_real(fid,vid(7),uz)
      rc = nf_put_var_real(fid,vid(8),vz)
      rc = nf_put_var_real(fid,vid(9),tz)
      rc = nf_put_var_real(fid,vid(10),sz)
      rc = nf_put_var_real(fid,vid(11),rz)

      rc = nf_put_var_real(fid,vid(12),trz(:,:,:,1))

      rc = nf_close(fid)
      end subroutine write_netcdf

      logical function is_suffix_nc( file_name ) result(answer)
C     Returns true if operation is successful, false  otherwise.
      character(len=*), intent(in) :: file_name
      integer :: i
      answer = .false.
C  look for last period
      i = index(file_name, ".", back=.true.) 
C  if no period was found
      if( i==0 ) go to 5 
c  after . should be ONLY two characters
      if( len(trim(file_name)) /= i+2 ) go to 5
      if( file_name(i:i+2) == ".nc" ) answer = .true.
C  normal exit
      return
C error trap
    5 answer = .false.
      return
      end function is_suffix_nc

      end program latlonz3d
