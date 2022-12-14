      subroutine regrid_input(grid)
!@sum Input file regridding routines. Uses x_2gridsroot derived type
!@    to define source and target grids.
!@    Can be called anywhere from modelE.f after dist_grid has been initialized
!@    Will soon become a standalone program
!@auth Denis Gueyffier
      use regrid_com
      use GEOM, only : geom_cs
      implicit none
      type (dist_grid), intent(in) :: grid
      type (x_2gridsroot) :: xll2cs
      integer :: imsource,jmsource,ntilessource,imtarget,
     &     jmtarget,ntilestarget

c      imsource=72
c      jmsource=46
c      imsource=288
c      jmsource=180

c     for RVR F(144x90)->CS32
c      imsource=144
c      jmsource=90

c     for RVR 1x1->CS90
      imsource=360
      jmsource=180

c     for STN
c      imsource=720
c      jmsource=360
      ntilessource=1

      imtarget=grid%IM_WORLD
      jmtarget=grid%JM_WORLD
      ntilestarget=grid%ntiles

c      call init_regrid_root(xll2cs,imsource,jmsource,ntilessource,
c     &     imtarget, jmtarget,ntilestarget)


c***   regrid* = conservative regridding
c***   interp* = bilinear interpolation
      write(*,*) "IN REGRID INPUT"

      if (AM_I_ROOT()) then
c         call regridTOPO(xll2cs)
c         call regridOSST(xll2cs)
c         call regridSICE(xll2cs)
c         call regridCDN(xll2cs)
c         call regridVEG(xll2cs)
c         call regridCROPS(xll2cs)
c         call regridTOPINDEX(xll2cs)
c         call regridSOIL(xll2cs)
c         call regridGLMELT(xll2cs)
c         call regridVEGFRAC(xll2cs)
c         call regridLAI(xll2cs)
      endif
c     call regridGIC(xll2cs,grid)
c     call regridAIC(xll2cs,grid)

c  Workflow for computation of river directions on cubed sphere
c  1) call regridRDSCAL(grid,720,360,...) in regridinput.f. It regrids scalar distance to the ocean 
c     from latlon grid to cube sphere. STN-30p.bin -> STN_CS32 \ STN_CS60 \ STN_CS90
c  2) run program contained in aux/RDto0.CS32.F (..RDtoO.CS90.F), which converts e.g. STN_CS32 to RDtoO.CS32
c  3) river directions can be plotted using aux/RDVCS32.F (..RDVCS90.F)
c  4) RDtoO.CS32 is mapped from (i,j,tile) -> (lat,lon) and the river mouths are added using RDijk2ll_CS 
c     in regridinput. The ascii files mouthnames_DtoO_CS32 and mouthij_DtoO_CS32 contain the mouth names and 
c     (i,j,tile) indices.
c  The resulting file is RDdistocean_CS32.bin

      call geom_cs
c      call interpRD(grid,imsource,jmsource,imtarget,jmtarget)
c      call regridRDSCAL(grid,imsource,jmsource,ntilessource)
      call RDijk2ll_CS(grid,imtarget,jmtarget)

      end subroutine regrid_input
c*


#ifdef CUBED_SPHERE
      subroutine interpRD(grid,ims,jms,imt,jmt)
!@sum  Interpolates river directions from latlon grid to cubed sphere 
!@+    using the following approach:
!@+         1- The latlon river directions vector field is bilinearly 
!@+           interpolated to the cubed sphere. The resulting vector 
!@+           field is still expressed in the polar coordinates basis
!@+         2- The vector field is then projected on a XY panel coordinate 
!@+           system where cube sphere edges form perpendicular lines. 
!@+           This allows us to avoid computing angles in cubed sphere metric.
!@+         3- We compute "quantized" river directions by following vector 
!@+           and picking corresponding rectangular cubed sphere cell 
!@auth Denis Gueyffier
      use regrid_com
      use DOMAIN_DECOMP_ATM, only: halo_update,pack_data
      use geom, only : ll2csxy_vec, lon2d_dg,lat2d_dg
      implicit none
      include 'mpif.h'
      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: twopi = 2d0*pi           !@param twopi 2*pi
      real*8,parameter :: radian = pi/180d0        !@param radian pi/180
      real*8, parameter :: g=0.615479708670387d0 ! g=.5*acos(1/3)
      integer, parameter :: nrvrmx=42
      type (dist_grid), intent(in) :: grid
      integer, intent(in) :: ims,jms,imt,jmt
      character*80 :: titlei,title2,title,name,namecs
      real*4, dimension(ims,jms) :: down_lat,down_lon,down_lat_911
     *     ,down_lon_911
      real*8,allocatable,dimension(:,:) :: ull,vll
      real*4, dimension(nrvrmx) :: lat_rvr,lon_rvr
      character*8, dimension(nrvrmx) :: namervr
      real*8, allocatable :: ucs_loc(:,:),vcs_loc(:,:)
      real*8 :: eps, meps, norm, lon, lat, dlon_dg, dlat_dg,fjeq
      real*8 :: x,y,VX,VY,slope,slope_diag,Xhalf,Yhalf,Xone,Yone
      integer iu_RVR, iu_RVRCS
      integer :: nrvr, i,j, nocomp, quadrant,nij,ierr,
     &     i0,i1,j0,j1,i0h,i1h,j0h,j1h,ii,jj,ti,imj,jmi,nijt,imouth
      integer, allocatable, dimension(:,:) :: kdirec
      real*8, allocatable, dimension(:,:) :: idown,jdown,kdown
      integer*4, dimension(imt,jmt,6) :: idown_g,jdown_g,kdown_g
      real*8,dimension(imt,jmt,6) :: idown_glob,jdown_glob,kdown_glob
     &     ,lon_glob,lat_glob

      allocate(kdirec(grid%I_STRT:grid%I_STOP,
     &     grid%J_STRT:grid%J_STOP))
      allocate(
     &     ull(ims,jms),vll(ims,jms),
     &     ucs_loc(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     vcs_loc(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     idown(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     jdown(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO),
     &     kdown(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO) )

      iu_RVR=20
      iu_RVRCS=21
c      name="RD_modelE_Fa.RVR.bin"
c      namecs="RD_CS32.bin"
      name="RD_modelE_O.RVR.bin"
      namecs="RD_CS90.bin"

      write(*,*) name

      imouth = 1  ! if 0 then input file doesn't contain river mouths

      if (am_i_root()) then
         open( iu_RVR, FILE=name,FORM='unformatted', STATUS='old')
         read(iu_RVR) titlei
         write(*,*) titlei
         if (imouth .eq. 1) then
            read(iu_RVR) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *           ,lon_rvr(1:nrvr)
         else
            read(iu_RVR) title2
         endif

         write(*,*) title2
         read(iu_RVR) title,down_lat
         write(*,*) title
         read(iu_RVR) title,down_lon
         write(*,*) title
         read(iu_RVR) title,down_lat_911
         write(*,*) title
         read(iu_RVR) title,down_lon_911
         write(*,*) title
         close(iu_RVR)
         
         write(110,*) "down_lat",down_lat
         
c     
c     create vector field on latlon grid
c     
         eps=1.d-6
         meps=-1.d-6
         dlat_dg=180./real(jms) ! even spacing (default)
         IF (JMS.eq.46) dlat_dg=180./REAL(jms-1) ! 1/2 box at pole for 4x5
         dlon_dg = 360./real(ims)
         fjeq=0.5*(1+jms)
         
        do j=1,jms
            if (j .ne. 1 .and. j .ne. jms) then
               lat = dlat_dg*(j - fjeq)
            elseif (j .eq. 1) then
               lat = -90.
            elseif (j .eq. jms) then
               lat = 90.
            endif
            do i=1,ims
               lon = -180.+(i-0.5)*dlon_dg
               if (i .eq. 1) lon = -180.+0.5*dlon_dg
               if (abs(down_lon(i,j)) .gt. 1./eps) then
                  ull(i,j)=0.d0
               else
                  ull(i,j)=(down_lon(i,j)- lon)*cos(lat*radian)
               endif
               if (abs(down_lat(i,j)) .gt. 1./eps) then
                  vll(i,j)=0.d0
               else
                  vll(i,j)=down_lat(i,j) - lat 
               endif
               
c*    Test cases of bilinear interpolation: uncomment lines below
c     For normal use, keep lines commented 
c     test1: checking if we interpolate exactly linear and constant fields
c     ull(i,j)=(2*i+3*j) ; vll(i,j) = (i+j)
c     ull(i,j)=1. ; vll(i,j)=-1.
c     test2: checking boundary conditions using periodic functions
c     ull(i,j)=cos(lon*radian) ; vll(i,j)=sin(lon*radian)
c     write(60,*) lon,lat,vll(i,j)
c*    end test cases               

               write(70,200) lon,lat,ull(i,j),vll(i,j)
 200              format(4(1X,f8.3))

        if (lat .ge. 37. .and. abs(down_lon(i,j)) .lt. 1./eps
     &                 .and. abs(down_lat(i,j)) .lt. 1./eps) then
            write(100,200) 
     &         cos(lat*radian)*cos(lon*radian),
     &         cos(lat*radian)*sin(lon*radian),
     &         cos(down_lat(i,j)*radian)*cos(down_lon(i,j)*radian)
     &        -cos(lat*radian)*cos(lon*radian),
     &         cos(down_lat(i,j)*radian)*sin(down_lon(i,j)*radian)
     &        -cos(lat*radian)*sin(lon*radian)
         endif
            enddo
         enddo
      end if

c     broadcast latlon field (ull,vll) to all PES
      nij=ims*jms
      call MPI_BCAST( ull, nij, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 
      call MPI_BCAST( vll, nij, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr ) 


c     interpolate vectors field from latlon points to cubed-sphere points
c     Note that the returned vector field V_cs = (ucs_loc, vcs_loc) is expressed 
c     in the latlon coordinate system
      call bilin_ll2cs_vec(grid,ull,vll,ucs_loc,vcs_loc,ims,jms)

      
      call halo_update(grid,ucs_loc)
      call halo_update(grid,vcs_loc)

      ti=grid%tile

      i0=grid%i_strt
      i0h=grid%i_strt_halo
      j0=grid%j_strt
      j0h=grid%j_strt_halo

      i1=grid%i_stop
      i1h=grid%i_stop_halo
      j1=grid%j_stop
      j1h=grid%j_stop_halo

      write(*,400) "i0h i0 i1 i1h j0h j0 j1 j1h",i0h,i0,i1,i1h,
     &     j0h,j0,j1,j1h
 400   format (A,8(I3))

c*    gather and broadcast lat/lon coordinates of CS cell centers
c*    for plotting purpose only
      call pack_data(grid,lon2d_dg,lon_glob)
      call pack_data(grid,lat2d_dg,lat_glob)
      nijt=imt*jmt*6
      call MPI_BCAST(lat_glob,nijt,MPI_DOUBLE_PRECISION,
     *     0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lon_glob,nijt,MPI_DOUBLE_PRECISION,
     *     0,MPI_COMM_WORLD,ierr )
 
      do j=j0,j1
         do i=i0,i1
            
c     write(20+ti,200) lon2d_dg(i,j),lat2d_dg(i,j),
c     &        vcs_loc(i,j)
            write(60+ti,200) lon2d_dg(i,j),lat2d_dg(i,j)
     &           ,2.*ucs_loc(i,j),2.*vcs_loc(i,j)
            
            nocomp=0
            
c     x,y coordinates of current CS cell center
            x = -1.d0 + 2d0*(i-0.5)/imt
            y = -1.d0 + 2d0*(j-0.5)/imt
            
c     compute components of V_cs in XY plane V_cs = (VX,VY)
            if ( sqrt(ucs_loc(i,j)*ucs_loc(i,j)
     &           +vcs_loc(i,j)*vcs_loc(i,j)) .lt. 1.e-4) then
               kdirec(i,j) = 0
               nocomp =1
               
c               write(*,*) "SMALL (ZERO) UCS VCS"
            else
               call ll2csxy_vec(x,y,ti,ucs_loc(i,j),
     &              vcs_loc(i,j),VX,VY)
               
c     write(*,100) "ucs vcs VX VY=",ucs_loc(i,j),vcs_loc(i,j),
c     &           VX,VY
c     100        format(A,4(1X,F8.2))
               
c     we assume V lies in first quadrant, VX > 0 && VY > 0. 
c     All other cases are deduced by symmetry
               if (VX .gt. eps .and. VY .gt. eps) then
                  quadrant = 1
               elseif ( VX .lt. meps .and. VY .gt. eps) then
                  quadrant=2
                  VX=-VX   
               elseif ( VX .lt. meps .and. VY .lt. meps) then
                  quadrant=3
                  VX=-VX   
                  VY=-VY   
               elseif ( VX .gt. eps .and. VY .lt. meps) then
                  quadrant=4
                  VY=-VY   
               elseif (abs(VX) .le. eps .and. abs(VY) .gt. eps) then
                  if (VY .gt. eps) then
                     kdirec(i,j)=2
                     nocomp=1
                  elseif (VY .lt. meps) then
                     kdirec(i,j)=6
                     nocomp=1
                  endif
               elseif (abs(VY) .le. eps .and. abs(VX) .gt. eps) then
                  if (VX .gt. eps) then
                     kdirec(i,j)=8
                     nocomp=1
                  elseif (VX .lt. meps) then
                     kdirec(i,j)=4
                     nocomp=1
                  endif
               elseif (abs(VX) .le. eps .and. abs(VY) .lt. eps) then
                  kdirec(i,j)=0
                  nocomp=1
               endif
            
            endif    
           
            if (nocomp .ne. 1) then
c     
c     compare slope VY/VX with (Yj+1-Yj+1/2)/(Xi+1-Xi+1/2) to decide which is the prefered direction
c     
               slope=VY/VX
c     write(*,*) "slope=",slope
               
               xone = -1.d0 + 2d0*i/imt ! (xone,yone) upper right corner of cell
               yone = -1.d0 + 2d0*j/imt
               Xone = tan(g*xone)*sqrt(2d0)
               Yone = tan(g*yone)*sqrt(2d0)
               Xhalf = tan(g*x)*sqrt(2d0)
               Yhalf = tan(g*y)*sqrt(2d0)
               slope_diag=(Yone-Yhalf)/(Xone-Xhalf)
               
               if (slope/slope_diag .lt. 1.d0) then
                  if (slope/slope_diag .lt. 0.5d0) then
                     if (quadrant .eq. 1) then
                        kdirec(i,j)=8
                     elseif (quadrant .eq. 2) then 
                        kdirec(i,j)=4
                     elseif (quadrant .eq. 3) then
                        kdirec(i,j)=4
                     elseif (quadrant .eq. 4) then
                        kdirec(i,j)=8
                     endif
                  else
                     if (quadrant .eq. 1) then 
                        kdirec(i,j)=1
                     elseif (quadrant .eq. 2) then
                        kdirec(i,j)=3
                     elseif (quadrant .eq. 3) then
                        kdirec(i,j)=5
                     elseif (quadrant .eq. 4) then
                        kdirec(i,j)=7
                     endif
                  endif
               else   
                  if (slope_diag/slope .lt. 0.5) then 
                     if (quadrant .eq. 1) then 
                        kdirec(i,j)=2
                     elseif (quadrant .eq. 2) then
                        kdirec(i,j)=2
                     elseif (quadrant .eq. 3) then
                        kdirec(i,j)=6
                     elseif (quadrant .eq. 4) then
                        kdirec(i,j)=6
                     endif
                  else
                     if (quadrant .eq. 1) then 
                        kdirec(i,j)=1
                     elseif (quadrant .eq. 2) then
                        kdirec(i,j)=3
                     elseif (quadrant .eq. 3) then
                        kdirec(i,j)=5
                     elseif (quadrant .eq. 4) then
                        kdirec(i,j)=7
                     endif
                  endif
               endif
               
            endif
            
c            write(*,*) "KDIREC(i,j)=",kdirec(i,j)
            
c     
c     The following is for Gary Russell's plotting program
c     
            KDOWN(i,j) = real(ti)

            SELECT CASE (KDIREC(i,j))
         CASE (0)
            IDOWN(i,j) = real(i)
            JDOWN(i,j) = real(j)

         CASE (1)
            IDOWN(i,j) = real(i+1)
            JDOWN(i,j) = real(j+1)

         CASE (2)
            IDOWN(i,j) = real(i)
            JDOWN(i,j) = real(j+1)

         CASE (3)
            IDOWN(i,j) = real(i-1)
            JDOWN(i,j) = real(j+1)

         CASE (4)
            IDOWN(i,j) = real(i-1)
            JDOWN(i,j) = real(j)

         CASE (5)
            IDOWN(i,j) = real(i-1)
            JDOWN(i,j) = real(j-1)

         CASE (6)
            IDOWN(i,j) = real(i)
            JDOWN(i,j) = real(j-1)

         CASE (7)
            IDOWN(i,j) = real(i+1)
            JDOWN(i,j) = real(j-1)

         CASE (8)
            IDOWN(i,j) = real(i+1)
            JDOWN(i,j) = real(j)

         END SELECT
         
      enddo
      enddo

 300  format(A,I3,I3,I3,A,F8.3,F8.3,F8.3)

c     output direction @ i= 53, j=51, tile = 3
      if (ti .eq. 6) write(*,*) "kdirec(53,51)=",kdirec(53,51)
      
      do j=j0,j1
         do i=i0,i1
            imj=i1-j
            jmi=j1-i
            if (idown(i,j) .eq. i1h) then
               ii=i1-j
               if (ti.eq.1) then
                   kdown(i,j)=2
                   idown(i,j)=1
                   jdown(i,j)=j
               endif
               if (ti.eq.2) then
                   kdown(i,j)=4
                   idown(i,j)=imj
                   jdown(i,j)=1
               endif
               if (ti.eq.3) then
                   kdown(i,j)=4
                   idown(i,j)=1
                   jdown(i,j)=j
               endif
               if (ti.eq.4) then
                   kdown(i,j)=6
                   idown(i,j)=imj
                   jdown(i,j)=1
               endif
               if (ti.eq.5) then
                   kdown(i,j)=6
                   idown(i,j)=1
                   jdown(i,j)=j 
               endif
               if (ti.eq.6) then
                   kdown(i,j)=2
                   idown(i,j)=imj
                   jdown(i,j)=1 
               endif         
               write(*,300) "case I+ ->i j tile =",i,j,ti,
     &          "// id jd kd",idown(i,j),jdown(i,j),kdown(i,j)
            elseif (idown(i,j) .eq. i0h .and. (
     &              kdirec(i,j).eq.3 .or.
     &              kdirec(i,j).eq.4 .or. 
     &              kdirec(i,j).eq.5      )) then
               if (ti.eq.1) then 
                   kdown(i,j)=5
                   idown(i,j)=imj
                   jdown(i,j)=j1
               endif
               if (ti.eq.2) then 
                   kdown(i,j)=1
                   idown(i,j)=i1
                   jdown(i,j)=j
               endif
               if (ti.eq.3) then
                   kdown(i,j)=1
                   idown(i,j)=imj
                   jdown(i,j)=j1
               endif
               if (ti.eq.4) then
                   kdown(i,j)=3
                   idown(i,j)=i1
                   jdown(i,j)=j
               endif
               if (ti.eq.5) then
                   kdown(i,j)=3
                   idown(i,j)=imj
                   jdown(i,j)=j1  
               endif
               if (ti.eq.6) then
                   kdown(i,j)=5
                   idown(i,j)=i1
                   jdown(i,j)=j
               endif       
               write(*,300) "case I- ->i j tile =",i,j,ti,
     &          "// id jd kd",idown(i,j),jdown(i,j),kdown(i,j)
             endif

            if (jdown(i,j) .eq. j1h) then
               if (ti.eq.1) then 
                   kdown(i,j)=3
                   idown(i,j)=1
                   jdown(i,j)=jmi
               endif
               if (ti.eq.2) then
                   kdown(i,j)=3
                   idown(i,j)=i
                   jdown(i,j)=1
               endif
               if (ti.eq.3) then 
                   kdown(i,j)=5
                   idown(i,j)=1
                   jdown(i,j)=jmi
               endif
               if (ti.eq.4) then
                   kdown(i,j)=5
                   idown(i,j)=i
                   jdown(i,j)=1
               endif
               if (ti.eq.5) then
                   kdown(i,j)=1
                   idown(i,j)=1
                   jdown(i,j)=jmi
               endif  
               if (ti.eq.6) then
                   kdown(i,j)=1
                   idown(i,j)=i
                   jdown(i,j)=1
               endif        
               write(*,300) "case J+ ->i j tile =",i,j,ti,
     &          "// id jd kd",idown(i,j),jdown(i,j),kdown(i,j)

             elseif (jdown(i,j) .eq. j0h .and. (
     &              kdirec(i,j).eq.5 .or.
     &              kdirec(i,j).eq.6 .or. 
     &              kdirec(i,j).eq.7      )) then
               if (ti.eq.1) then
                   kdown(i,j)=6
                   idown(i,j)=i
                   jdown(i,j)=j1
               endif
               if (ti.eq.2) then
                   kdown(i,j)=6
                   idown(i,j)=i1
                   jdown(i,j)=jmi
               endif
               if (ti.eq.3) then
                   kdown(i,j)=2
                   idown(i,j)=i
                   jdown(i,j)=j1
               endif
               if (ti.eq.4) then
                   kdown(i,j)=2
                   idown(i,j)=i1
                   jdown(i,j)=jmi
               endif
               if (ti.eq.5) then
                   kdown(i,j)=4
                   idown(i,j)=i
                   jdown(i,j)=j1 
               endif 
               if (ti.eq.6) then
                   kdown(i,j)=4
                   idown(i,j)=i1
                   jdown(i,j)=jmi
               endif
               write(*,300) "case J- ->i j tile =",i,j,ti,
     &          "// id jd kd",idown(i,j),jdown(i,j),kdown(i,j)
             endif
c* latlon projection
            write(70+ti,200) lon2d_dg(i,j),lat2d_dg(i,j),
     &           lon_glob(idown(i,j),jdown(i,j),kdown(i,j))
     &           -lon2d_dg(i,j),
     &           lat_glob(idown(i,j),jdown(i,j),kdown(i,j))
     &           -lat2d_dg(i,j)

c*    polar-spherical projection
            write(80+ti,200) 
     &         cos(lat2d_dg(i,j)*radian)*cos(lon2d_dg(i,j)*radian),
     &         cos(lat2d_dg(i,j)*radian)*sin(lon2d_dg(i,j)*radian),
     &         cos(lat_glob(idown(i,j),jdown(i,j),kdown(i,j))*radian)
     &        *cos(lon_glob(idown(i,j),jdown(i,j),kdown(i,j))*radian )
     &        -cos(lat2d_dg(i,j)*radian)*cos(lon2d_dg(i,j)*radian),
     &         cos(lat_glob(idown(i,j),jdown(i,j),kdown(i,j))*radian)
     &        *sin(lon_glob(idown(i,j),jdown(i,j),kdown(i,j))*radian )
     &        -cos(lat2d_dg(i,j)*radian)*sin(lon2d_dg(i,j)*radian)      
         enddo
      enddo

c     gather all river directions before writing to output file
      call pack_data(grid,idown,idown_glob)
      call pack_data(grid,jdown,jdown_glob)
      call pack_data(grid,kdown,kdown_glob)
      
      
      if (am_i_root()) then
         idown_g=idown_glob
         jdown_g=jdown_glob
         kdown_g=kdown_glob
c     write(iu_RVRCS) titlei
         title="downstream idown, jdown ,kdown"
         open( iu_RVRCS, FILE=namecs,FORM='unformatted', 
     &        STATUS='unknown')
         write(iu_RVRCS) title,idown_g,jdown_g,kdown_g
         close(iu_RVRCS)
      endif
      
c     end output for Gary Russell's plotting program
      
      deallocate(ull,vll,kdirec,ucs_loc,vcs_loc,idown,jdown,kdown)

      end subroutine interpRD
c*
#endif


      subroutine RDijk2ll_CS(grid,imt,jmt)
!@sum On CS grid, convert i,j,k (with k=cube face index from 1..6)
!@+   to absolute lat-lon coordinates
!@auth Denis Gueyffier
      use geom, only : lon2d_dg,lat2d_dg
      use regrid_com
      use domain_decomp_atm, only : pack_data 
      implicit none
      type (dist_grid), intent(in) :: grid
      integer, intent(in) :: imt,jmt
      integer, parameter :: nrvr = 41 ! # of river mouths
      character*80 :: name,nameout,title,
     &     title1,title2,title3,title4,title5,title6
      real*8,parameter :: undef=-1.d30  !missing value
      real*4, dimension(imt,jmt,6) :: down_lat,down_lon
      real*4, dimension(imt,jmt,6) :: down_lat_911,down_lon_911
      real*8, dimension(imt,jmt,6) :: lon_glob,lat_glob
      integer, dimension(imt,jmt,6) :: idown,jdown,kdown
      integer :: iu_RD,iu_TOPO,iu_MNAME,iu_MIJ,i,j,k
      LOGICAL, dimension(imt,jmt) :: NODIR
      real*4 :: FOCEAN(imt,jmt,6)
      character*8,dimension(nrvr) :: namervr
      character*2,dimension(nrvr) :: mouthI,mouthJ
      character*1,dimension(nrvr) :: mouthK
      integer,dimension(nrvr) :: imouthI,imouthJ,imouthK
      real*4,dimension(nrvr) :: lat_rvr,lon_rvr
      iu_TOPO=30

      if (am_i_root()) then
c*    Read ocean fraction
      name="Z_CS32"
      open(iu_TOPO,FILE=name,FORM='unformatted', STATUS='old')
      read(iu_TOPO) title,FOCEAN
      close(iu_TOPO)

c*    Read names of river mouths
      iu_MNAME=20
      name="mouthnames_DtoO_CS32"
      open(iu_MNAME,FILE=name,FORM='formatted', STATUS='old')
      READ (iu_MNAME,'(A8)') (namervr(I),I=1,nrvr) !Read mouths names
      write(*,*) namervr
      close(iu_MNAME)

c*    Read i,j,k coordinates of river mouths (k = index of cube face)
      name="mouthij_DtoO_CS32"
      open(iu_MIJ,FILE=name,FORM='formatted', STATUS='old')
      READ (iu_MIJ,'(A2,1X,A2,1X,A1)') (      
     &       mouthI(I),mouthJ(I),mouthK(I), I=1,nrvr)
      close(iu_MIJ)

c*     conversion char to int
      do i=1,nrvr
          read(mouthI(i),'(I2)') imouthI(i)
          read(mouthJ(i),'(I2)') imouthJ(i)
          read(mouthK(i),'(I1)') imouthK(i)
          write(*,*) imouthI(i),imouthJ(i),imouthK(i)
      enddo
    
      nameout="RDdistocean_CS32.bin"
    
      iu_RD=20
      name="RDtoO.CS32"
      write(*,*) name,imt,jmt

c*    read i,j,k coordinates of downstream cells

         open( iu_RD, FILE=name,FORM='unformatted',
     &        STATUS='unknown')
         read(iu_RD) title,idown,jdown,kdown
         close(iu_RD)
      write(*,*) "read RDijk2ll_CS"
      endif

c*    form global arrays for mapping (i,j,k) -> (lon,lat)
      call pack_data(grid,lon2d_dg,lon_glob)
      call pack_data(grid,lat2d_dg,lat_glob)

      if (am_i_root()) then

c*    convert i,j,k coordinates of river mouths to absolute lat-lon coordinates    
      do k=1,nrvr
            lat_rvr(k)=lat_glob(imouthI(k),imouthJ(k)
     &            ,imouthK(k))
            lon_rvr(k)=lon_glob(imouthI(k),imouthJ(k)
     &            ,imouthK(k))
      enddo

c*    convert i,j,k coordinates of downstream cells to absolute lat-lon coordinates    
      do k=1,6
        do j=1,jmt
          do i=1,imt
            if (FOCEAN(i,j,k) .lt. 1.e-6) then
              down_lat(i,j,k)=lat_glob(idown(i,j,k),jdown(i,j,k)
     &            ,kdown(i,j,k))
              down_lon(i,j,k)=lon_glob(idown(i,j,k),jdown(i,j,k)
     &            ,kdown(i,j,k))

c*    dummy emergency directions
              down_lat_911(i,j,k)=down_lat(i,j,k)
              down_lon_911(i,j,k)=down_lon(i,j,k)
            write(130+k,200) lon_glob(i,j,k),lat_glob(i,j,k),
     &           down_lon(i,j,k)-lon_glob(i,j,k),
     &           down_lat(i,j,k)-lat_glob(i,j,k)
            else
              down_lat(i,j,k)=undef
              down_lon(i,j,k)=undef
              down_lat_911(i,j,k)=undef
              down_lon_911(i,j,k)=undef
            write(130+k,200) lon_glob(i,j,k),lat_glob(i,j,k),
     &           0.,0.
            end if

          enddo
        enddo
      enddo

c*    output everything 
      open(iu_RD,FILE=nameout,FORM='unformatted',
     &        STATUS='unknown')
      title1="river directions from dist. to ocean, CS32, 03/10/09"
      title2="Named River Mouths:"
      title3="Latitude of downstream river direction box"
      title4="Longitude of downstream river direction box"
      title5="Latitude of emergency downstream river direction box"
      title6="Longitude of emergency downstream river direction box"

      write(iu_RD) title1   
      write(iu_RD) title2,nrvr,namervr(1:nrvr),lat_rvr(1:nrvr)
     *     ,lon_rvr(1:nrvr)
      write(iu_RD) title3,down_lat
      write(iu_RD) title4,down_lon
      write(iu_RD) title5,down_lat_911
      write(iu_RD) title6,down_lon_911
      close(iu_RD)

      write(*,*) "wrote RD file"

      endif

 200  format(4(1X,f8.3))

      end subroutine RDijk2ll_CS
c*


      subroutine regridTOPO(x2grids)
c
      use regrid_com
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE(nrecmax),name,ncfile
      type (x_2gridsroot), intent(inout) :: x2grids
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      real*8, allocatable :: tsource(:,:,:,:)
      real*4, allocatable :: ttargr4(:,:,:,:)
      real*4, allocatable :: tsourc4(:,:,:,:)
      integer :: iu_TOPO,i,j,k,irec,iunit,imt,jmt,ntt,ims,jms,nts,
     &     maxrec,status,vid,fid,ir

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget
      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource


      write(*,*) "imt jmt ntt",imt,jmt,ntt

      allocate( ttargglob(imt,jmt,ntt,nrecmax),
     &     ttargr4(imt,jmt,ntt,nrecmax),
     &     ones(imt,jmt,ntt), 
     &     tsource(ims,jms,nts,nrecmax),
     &     tsourc4(ims,jms,nts,nrecmax)  )
      
      iu_TOPO=20
      name="Z288X180N"
c      name="Z72X46N.cor4_nocasp"
c      name="Z144X90N_nocasp"

      write(*,*) name
      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='old')

c      write(*,*) "iu_TOPO",iu_TOPO
c      call read_regrid_4D_1R(x2grids,iu_TOPO,TITLE,ttargglob,maxrec)

      tsource(:,:,:,:)=0.
      ttargglob(:,:,:,:)=0.

      write(*,*) "ims, jms, nts=",ims,jms,nts
            
      irec=1

      do
         read(unit=iu_TOPO,END=30) TITLE(irec), tsourc4(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= tsourc4(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      write(*,*) "maxrec",maxrec
      
      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),
     &        ttargglob(:,:,:,ir))
        write(*,*) "TITLE",TITLE(ir)
      enddo

c
c     CONSISTENCY CHECKS: 1) FOCEAN+FLAKE+FGRND+FGICE=1 
c                         2) IF FOCEAN(i,j) > 0 set FGRND=FGRND+FLAKE, FLAKE=0

      ones(:,:,:)=ttargglob(:,:,:,1)
     &     +ttargglob(:,:,:,2)
     &     +ttargglob(:,:,:,3)
     &     +ttargglob(:,:,:,4)
      
      if (any( abs(ones(1:imt,1:jmt,1:ntt) - 1.0) .gt. 1.d-6 ) ) then 
         write(*,*) "WARNING FOCEAN+FLAKE+FGRND+FGICE=1 BROKEN"
      endif

      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if ( (ttargglob(i,j,k,1) .gt. 0.) .and.
     &              (ttargglob(i,j,k,2) .gt. 0. ) ) then
                  write(*,*) "FGRND=FGRND+FLAKE, FLAKE=0"
                  ttargglob(i,j,k,3)=ttargglob(i,j,k,2)
     &                 +ttargglob(i,j,k,3)
                  ttargglob(i,j,k,2)=0.0
               endif
            enddo
         enddo
      enddo

      close(iu_TOPO)

c      name="Z_CS48"
      name="Z_CS90"

      write(*,*) name

      open( iu_TOPO, FILE=name,FORM='unformatted', STATUS='unknown')

      ttargr4=ttargglob

      do ir=1,maxrec
         write(unit=iu_TOPO) TITLE(ir), ttargr4(:,:,:,ir)
         write(*,*) TITLE(ir), ttargr4(:,:,:,ir)
      enddo

      close(iu_TOPO)


c      ncfile="topo6tiles.nc"

c      status = nf_open(trim(ncfile),nf_write,fid)
c      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN FILE"
c      status = nf_inq_varid(fid,'zatmo',vid)
c      write(*,*) NF_STRERROR(status)
c      status = nf_put_var_double(fid,vid,ttargglob(:,:,:,1))
c      write(*,*) "STATUS",NF_STRERROR(status),"<<"
c      status = nf_close(fid)

c      ncfile="topo.nc"

c      status = nf_open(trim(ncfile),nf_write,fid)
c      if (status .ne. NF_NOERR) write(*,*) "UNABLE TO OPEN FILE"
c      status = nf_inq_varid(fid,'zatmo',vid)
c      write(*,*) NF_STRERROR(status)
c      status = nf_put_var_double(fid,vid,ttargglob(:,:,2,1))
c      write(*,*) "STATUS",NF_STRERROR(status),"<<"
c      status = nf_close(fid)

      deallocate(ttargglob,ones,ttargr4,tsource)
   
      end subroutine regridTOPO 
c*



      subroutine regridRDSCAL(grid,ims,jms,nts)
c     
c     regridding scalar quantities used to derive river directions 
c     using simple sampling scheme (non conservative) 
c
      use regrid_com
      use geom, only : lon2d_dg,lat2d_dg
      use domain_decomp_atm, only : pack_data 
      implicit none
      include 'mpif.h'
      type (dist_grid) :: grid
      character*80 :: TITLE1,TITLE2,TITLE3,TITLE4,TITLE5,name
      real*8,allocatable :: ttargglob(:,:,:)
      real*8,allocatable :: tsource(:,:,:)
      real*4,allocatable :: ttargupst4(:,:,:)
      real*4,allocatable :: ttargdist4(:,:,:)
      real*4,allocatable :: tsourc4(:,:,:)
      integer, allocatable :: bassId(:,:,:),tbassId(:,:,:)
      integer:: iu_RDSCAL,i,j,k,iunit,imt,jmt,ntt,ims,jms,nts,
     &     status,nijt,ierr,ilon,jlat
      real*8, allocatable :: lon_glob(:,:,:),lat_glob(:,:,:)
      real*8 :: dlat_dg,dlon_dg

      imt=grid%IM_WORLD
      jmt=grid%JM_WORLD
      ntt=grid%ntiles

      write(*,*) "imt jmt ntt ims jms nts",imt,jmt,ntt,ims,jms,nts

      allocate( 
     &     ttargglob(imt,jmt,ntt),
     &     ttargupst4(imt,jmt,ntt),
     &     ttargdist4(imt,jmt,ntt),
     &     tbassId(imt,jmt,ntt),
     &     tsourc4(ims,jms,nts),
     &     tsource(ims,jms,nts),
     &     bassId(ims,jms,nts),
     &     lon_glob(imt,jmt,ntt),
     &     lat_glob(imt,jmt,ntt)
     &     )
      
      write(*,*) "alloc"

c*    gather and broadcast lat/lon coordinates of CS cell centers
      call pack_data(grid,lon2d_dg,lon_glob)
      call pack_data(grid,lat2d_dg,lat_glob)
      nijt=imt*jmt*6
      write(*,*) "aft pack"

      call MPI_BCAST( lat_glob, nijt, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr )
      call MPI_BCAST( lon_glob, nijt, MPI_DOUBLE_PRECISION,
     *     0, MPI_COMM_WORLD, ierr )
      write(*,*) "aft bcast"

      dlon_dg = 360./dble(ims)
      dlat_dg=180./real(jms) ! even spacing (default)
      IF (jms.eq.46) dlat_dg=180./REAL(jms-1) ! 1/2 box at pole for 4x5
c*

      if (am_i_root()) then
      iu_RDSCAL=20
      name="STN-30p.bin"

      write(*,*) name
      open( iu_RDSCAL, FILE=name,FORM='unformatted', STATUS='old')

      write(*,*) "ims, jms, nts=",ims,jms,nts
      read(iu_RDSCAL) TITLE1
      write(*,*) TITLE1
      read(iu_RDSCAL) TITLE2
      write(*,*) TITLE2
      read(iu_RDSCAL) TITLE3,bassId(:,:,:)
      do k=1,ntt
      do j=1,jmt
      do i=1,imt
      ilon=ims/2+1+int(real((lon_glob(i,j,k)+.01)/dlon_dg))
      jlat=jms/2+1+int(real((lat_glob(i,j,k)+.01)/dlat_dg))
      tbassId(i,j,k)=bassId(ilon,jlat,1)
      enddo
      enddo
      enddo
      write(*,*) TITLE3
      read(iu_RDSCAL) TITLE4,tsourc4(:,:,:)
      tsource(:,:,:)=tsourc4(:,:,:)
      do k=1,ntt
      do j=1,jmt
      do i=1,imt
      ilon=ims/2+1+int(real((lon_glob(i,j,k)+.01)/dlon_dg))
      jlat=jms/2+1+int(real((lat_glob(i,j,k)+.01)/dlat_dg))
      ttargdist4(i,j,k)=tsource(ilon,jlat,1)
      enddo
      enddo
      enddo
      write(*,*) TITLE4
      read(iu_RDSCAL) TITLE5,tsourc4(:,:,:)
      tsource(:,:,:)=tsourc4(:,:,:)
      do k=1,ntt
      do j=1,jmt
      do i=1,imt
      ilon = ims/2+1+int(real((lon_glob(i,j,k)+.01)/dlon_dg))
      jlat = jms/2+1+int(real((lat_glob(i,j,k)+.01)/dlat_dg))
      ttargupst4(i,j,k)=tsource(ilon,jlat,1)
      enddo
      enddo
      enddo
      write(*,*) TITLE5
      close (iu_RDSCAL)

c      name="STN_C90"
      name="STN_CS32"

      write(*,*) name
      open( iu_RDSCAL, FILE=name,FORM='unformatted', STATUS='unknown')
      write(iu_RDSCAL) TITLE1
      write(iu_RDSCAL) TITLE2
      write(iu_RDSCAL) TITLE3,tbassId(:,:,:)
      write(iu_RDSCAL) TITLE4,ttargdist4(:,:,:)
      write(iu_RDSCAL) TITLE5,ttargupst4(:,:,:)
      close(iu_RDSCAL)
      write(*,*) "bassid=",bassId  
      endif

      deallocate(ttargglob,tsource,ttargupst4,ttargdist4,bassId,
     &             tbassId,lon_glob,lat_glob)

      end subroutine regridRDSCAL
c*



      subroutine regridOSST(x2grids)
c     for 1x1 resolution :  Gary on Athena clima1/OBS/AMIP/1x1
c     Also SST1x1_HadISST from Hadley (Jeff Jonas) - check if this  
c     has been interpolated from lower resolution 

      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 :: TITLE,name
      real*4 OSTmean(x2grids%imsource,x2grids%jmsource),
     &     OSTend(x2grids%imsource,x2grids%jmsource)
      integer iu_OSST,i,j,k

      iu_OSST=20
      name="OST4X5.B.1876-85avg.Hadl1.1"
      if (am_i_root()) then
      open(iu_OSST, FILE=name,FORM='unformatted', STATUS='old')
      endif

      call read_regrid_write_4D_2R(x2grids,name,iu_OSST)

      if (am_i_root()) close(iu_OSST)

      end subroutine regridOSST
c


      subroutine regridSICE(x2grids)
c     for 1x1 resolution: Gary on Athena clima1/OBS/AMIP/1x1
c     Also ICE_1x1_HadISST from Hadley (Jeff Jonas) - check if this  
c     has been interpolated from lower resolution 
c     /u/cmrun/SICE4X5.B.1993-2002avg.Hadl1.1

      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name,outunformat, TITLE2(nrecmax)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tin(:,:,:),tout(:,:,:)
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
      integer iu_SICE,i,j,k,ims,jms,nts,imt,jmt,ntt,iuout,
     &    maxrec,ir

      iu_SICE=19
      name="SICE4X5.B.1876-85avg.Hadl1.1"
      write(*,*) name

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),
     &     tin(ims,jms,nts),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:)=0.0
 
      if (am_i_root()) then
      open (iu_SICE, FILE=name,FORM='unformatted', STATUS='old')
      read(unit=iu_SICE) TITLE, tin
      tsource=tin
      endif

      call root_regrid(x2grids,tsource(:,:,:),ttargglob)
      tout(:,:,:)=ttargglob(:,:,:)
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20

      if (am_i_root()) then
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      write(unit=iuout) TITLE,tout(:,:,:)
      write(*,*) TITLE
      endif

      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt),tbig(imt,jmt,2,ntt) )

      tsource1(:,:,:,:)=0.0
      tsource2(:,:,:,:)=0.0

      call read_recs_2R(tsource1,tsource2,iu_SICE,TITLE2,
     &     maxrec,ims,jms,nts)

      write(*,*) "maxrec",maxrec

      do ir=1,maxrec
         call root_regrid(x2grids,tsource1(:,:,:,ir),ttargglob1)
         tout1(:,:,:)=ttargglob1(:,:,:)
         call root_regrid(x2grids,tsource2(:,:,:,ir),ttargglob2)
         tout2(:,:,:)=ttargglob2(:,:,:)
         tbig(:,:,1,:)=tout1(:,:,:)
         tbig(:,:,2,:)=tout2(:,:,:)
         if (am_i_root()) then
         write(unit=iuout) TITLE2(ir),tbig
c         write(unit=iuout) TITLE2(ir),tout1(:,:,:),tout2(:,:,:)
         write(*,*) "TITLE",TITLE2(ir)
         endif
      enddo

      close(iuout)

      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2,
     &     tin,tout,tbig)

      if (am_i_root()) close(iu_SICE)
      
      end subroutine regridSICE



      subroutine regridCDN(x2grids)
c     for 1x1 resolution: Jeff Jonas uses CDN=AL30RL360X180N.rep

      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_CDN

      iu_CDN=20
      name="CD4X500S.ext"
      if (am_i_root()) then
      open (iu_CDN, FILE=name,FORM='unformatted', STATUS='old')
      endif

      call read_regrid_write_4D_1R(x2grids,name,iu_CDN)
      if (am_i_root()) close(iu_CDN)
      
      end subroutine regridCDN
c*




      subroutine regridVEG(x2grids)
c     for 1x1 resolution: Jeff uses VEG=V360X180_no_crops.rep
c     It is identical to V144X90_no_crops.ext (144X90 data 
c     was just transfered to 360X180 grid without any change)
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_VEG
	
      iu_VEG=20
      name="V72X46.1.cor2_no_crops.ext"

      if (am_i_root()) then
      open(iu_VEG,FILE=name,FORM='unformatted', STATUS='old')
c     &     convert='big_endian')
      endif

      call read_regrid_write_veg(x2grids,name,iu_VEG)
           
      if (am_i_root()) close(iu_VEG)
      
      end subroutine regridVEG
c*


      subroutine read_regrid_write_veg(x2grids,name,iuin)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:),data(:,:,:)
c      real*4 :: vadata(11,4,3)
      character*80 :: TITLE
      character*80 :: outunformat
      integer :: ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt

      allocate (tsource(ims,jms,nts),data(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),tout(imt,jmt,ntt))
c     arrsum(imt,jmt,ntt)
     
      
      outunformat=trim(name)//".CS"
      write(*,*) outunformat

      iuout=20

      if (am_i_root()) then
         open( iuout, FILE=outunformat,
     &        FORM='unformatted', STATUS="UNKNOWN")
      endif
      
      do ir=1,10
         if (am_i_root()) then
            read(unit=iuin) TITLE,data
            write(*,*) "TITLE, ir",TITLE,ir
            tsource= data
         endif
         call root_regrid(x2grids,tsource,ttargglob)
c     arrsum(:,:,:)=arrsum(:,:,:)+ttargglob(:,:,:)
         tout=ttargglob
         if (am_i_root()) then
            write(unit=iuout) TITLE,tout
         endif
      enddo

      if (am_i_root()) close(iuout) 

      deallocate(tsource,ttargglob,tout,data)

      end subroutine read_regrid_write_veg
c*


      subroutine regridCROPS(x2grids)
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(inout) :: x2grids
      character*80 name
      integer iu_CROPS

      iu_CROPS=20
      name="CROPS_72X46N.cor4.ext"
      if (am_i_root()) then
         open(iu_CROPS,FILE=name,FORM='unformatted', STATUS='old')
      endif

      call read_regrid_write_4D_1R(x2grids,name,iu_CROPS)

      if (am_i_root()) close(iu_CROPS)

      end subroutine regridCROPS
c*



      subroutine regridTOPINDEX(x2grids)
c     top_index_360x180.ij is not "extended". It has "-1" in ocean 
c     cells (meaning missing data), so you should be careful not to 
c     include those cells in the regridding computations. 
c     Alternatively you can use top_index_360x180.ij.rep which is 
c     direct transfer from top_index_144x90.ij.ext 

      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_TOP_INDEX

      iu_TOP_INDEX=20
      name="top_index_72x46.ij.ext"
      if (am_i_root()) then
         open(iu_TOP_INDEX,FILE=name,FORM='unformatted', STATUS='old')
      endif
      call read_regrid_write_4D_1R(x2grids,name,iu_TOP_INDEX)
      
      if (am_i_root()) close(iu_TOP_INDEX)
      
      end subroutine regridTOPINDEX
c*


      
      subroutine regridSOIL(x2grids)
c
c     for 1x1 resolution: Jeff uses SOIL=S360X180_0098M.rep
c
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      real*4, allocatable :: dz(:,:,:),ftext(:,:,:,:),
     &     ftextk(:,:,:,:),sl(:,:)
      real*4, allocatable :: dzout(:,:,:,:),ftextout(:,:,:,:,:),
     &     ftextkout(:,:,:,:,:),slout(:,:,:),bigarrout(:,:,:,:)
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLE,name,outunformat
      integer iu_SOIL,iuout,ims,jms,nts,imt,jmt,ntt,i,j,k,l,inds

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

      iu_SOIL=20
      name="S4X50093.ext"
      if (am_i_root()) then
      open(iu_SOIL,FILE=name,FORM='unformatted', STATUS='old')
      endif

      allocate (dz(ims,jms,6),ftext(ims,jms,6,5),
     &     ftextk(ims,jms,6,5),sl(ims,jms),
     &     dzout(imt,jmt,6,ntt),ftextout(imt,jmt,6,5,ntt),
     &     ftextkout(imt,jmt,6,5,ntt),
     &     slout(imt,jmt,ntt),bigarrout(imt,jmt,67,ntt) )

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt) )

      if (am_i_root()) then
         read(iu_SOIL) dz,ftext,ftextk,sl

         close(iu_SOIL)
      endif
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat

      write(*,*) dz
      iuout=20
      if (am_i_root()) then
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      endif

      do k=1,6
         do j=1,jms
            do i=1,ims
               tsource(i,j,1)=dz(i,j,k)
            enddo
         enddo
         call root_regrid(x2grids,tsource,ttargglob)
         write(*,*) "k=",k
         dzout(:,:,k,:)=ttargglob(:,:,:)
      enddo

      do k=1,6
         do l=1,5
            tsource(:,:,1)=ftext(:,:,k,l)
            call root_regrid(x2grids,tsource,ttargglob)
            ftextout(:,:,k,l,:)=ttargglob(:,:,:)
         enddo
      enddo


      do k=1,6
         do l=1,5
            tsource(:,:,1)=ftextk(:,:,k,l)
            call root_regrid(x2grids,tsource,ttargglob)
            ftextkout(:,:,k,l,:)=ttargglob(:,:,:)
         enddo
      enddo

      tsource(:,:,1)=sl(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      slout(:,:,:)=ttargglob(:,:,:)
      
      do k=1,6
         bigarrout(:,:,k,:)=dzout(:,:,k,:)
      enddo

      inds=6
      do k=1,6
         do l=1,5
            bigarrout(:,:,inds+k+6*(l-1),:)=ftextout(:,:,k,l,:)
         enddo
      enddo

      inds=inds+30
      do k=1,6
         do l=1,5
            bigarrout(:,:,inds+k+6*(l-1),:)=ftextkout(:,:,k,l,:)
         enddo
      enddo

      inds=inds+30
      bigarrout(:,:,inds+1,:)=slout(:,:,:)

      if (am_i_root()) then
         write(unit=iuout) bigarrout
         
         close(iuout)
      endif

      end subroutine regridSOIL
c*
      
      

      subroutine regridGLMELT(x2grids)
      use regrid_com
      implicit none
      type (x_2gridsroot), intent(in) :: x2grids
      character*80 TITLE,name
      integer iu_GLMELT

      name="GLMELT_4X5.OCN.bin"
      write(*,*) name

      iu_GLMELT=20
      
      if (am_i_root()) then
         open (iu_GLMELT,FILE=name,FORM='unformatted', STATUS='old')
      endif
   
      call read_regrid_write_GLMELT(x2grids,name,iu_GLMELT)

      if (am_i_root()) close(iu_GLMELT)
      
      end subroutine regridGLMELT
c*



      subroutine read_regrid_write_GLMELT(x2grids,name,iuin)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,itile,i,j

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:,:)=0.0

      if (am_i_root()) then
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)

      write(*,*) "maxrec",maxrec

      outunformat=trim(name)//".CS"

      write(*,*) outunformat

      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      endif

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         do itile=1,ntt
           do i=1,imt
             do j=1,jmt
             if (ttargglob(i,j,itile) .le. 0.) then
               tout(i,j,itile)=0.
             else
               tout(i,j,itile)=1.
             endif
            enddo
           enddo
          enddo
          if (am_i_root()) then
             write(unit=iuout) TITLE(ir),tout(:,:,:)
             write(*,*) "TITLE",TITLE(ir)
          endif
      enddo

      if (am_i_root()) close(iuout)

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_GLMELT
c*



      subroutine regridVEGFRAC(x2grids)
c
c     regriding vegetation fraction used by tracers code
c
      use regrid_com
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE(nrecmax),name,oname,TITLEFILE
      type (x_2gridsroot), intent(inout) :: x2grids
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      real*4, allocatable :: ttargr4(:,:,:,:)
      integer :: iu_VEGFRAC,iu_VEGFRACCS,i,j,k,irec,iunit,imt,jmt,
     &     ntt,maxrec,status,vid,fid,ir

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "imt jmt ntt",imt,jmt,ntt

      allocate( ttargglob(imt,jmt,ntt,nrecmax),
     &     ttargr4(imt,jmt,ntt,nrecmax),
     &     ones(imt,jmt,ntt) )

      ttargglob(:,:,:,:) = 0.0


      if (am_i_root()) then
      iu_VEGFRAC=20
      name="vegtype.global4X5.bin"
      write(*,*) name
      open( iu_VEGFRAC, FILE=name,FORM='unformatted', STATUS='old')
      endif

      call read_regrid_4D_1R(x2grids,iu_VEGFRAC,TITLE,ttargglob,maxrec)

c***  
c***  CONSISTENCY CHECKS: 1) sum(fractions)=1 
c***  

      ones=sum(ttargglob(1:imt,1:jmt,1:ntt,1:maxrec),4)
      
      if (any( abs(ones(1:imt,1:jmt,1:ntt) - 1.0) .gt. 1.d-3 ) ) then 
         write(*,*) "WARNING FRACTIONS DO NOT ADD UP TO 1 
     &        EVERYWHERE ON CUBED SPHERE"
         write(*,*) "ones=",ones
      endif

      if (am_i_root()) close(iu_VEGFRAC)

      oname=trim(name)//".CS"

      if (am_i_root()) then
         write(*,*) oname
         
         open( iu_VEGFRACCS, FILE=oname,FORM='unformatted', 
     &        STATUS='unknown')
      endif

      ttargr4=ttargglob

      if (am_i_root()) then
         do ir=1,maxrec
            write(unit=iu_VEGFRACCS) TITLE(ir), ttargr4(:,:,:,ir)
         enddo
         
         close(iu_VEGFRACCS)
      endif

      deallocate(ttargglob,ones,ttargr4)
   
      end subroutine regridVEGFRAC
c*




      subroutine regridLAI(x2grids)
c
c     regriding LAI used by tracers code
c
      use regrid_com
      implicit none
      include 'netcdf.inc'
      character*80 :: TITLE(nrecmax),name,oname,TITLEFILE
      type (x_2gridsroot), intent(inout) :: x2grids
      real*8, allocatable :: ttargglob(:,:,:,:),ones(:,:,:)
      real*4, allocatable :: ttargr4(:,:,:,:)
      integer :: iu_LAI,iu_LAICS,i,j,k,irec,iunit,imt,jmt,
     &     ntt,maxrec,status,vid,fid,ir,imonth
      character(len=2) :: c2month 
      character(len=1) :: c1month 

      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "imt jmt ntt",imt,jmt,ntt

      iu_LAI=20
      
      do imonth=1,12
         if (imonth .lt. 10) then
            write(c1month,'(i1)') imonth
            c2month="0"//c1month
         else
            write(c2month,'(i2)') imonth
         endif
         name='lai'//c2month//'.global.bin'
         write(*,*) name

         if (am_i_root()) then               
            open( iu_LAI, FILE=name,FORM='unformatted', STATUS='old')
         endif

         call read_regrid_write_4D_1R(x2grids,name,iu_LAI)
         
         if (am_i_root()) close(iu_LAI)
      enddo
   
      end subroutine regridLAI
c*



      subroutine regridGIC(x2grids,dd2d)
c
c     for 1x1 resolution: Jeff uses GIC=GIC.360X180.DEC01.1.rep
c
      use DOMAIN_DECOMP_1D,only : am_i_root
      use regrid_com
      use dd2d_utils
      use pario, only : defvar,write_data

      implicit none
      include 'netcdf.inc'

      type (x_2gridsroot), intent(in) :: x2grids
      type (dist_grid), intent(in) :: dd2d

c*    read
      real*8, allocatable :: Tocn(:,:,:),MixLD(:,:)        ! OCN01
      real*8, allocatable :: F(:,:),H(:,:,:),snw(:,:),msi(:,:),
     &     ssi(:,:,:),pond_melt(:,:)
      logical, allocatable :: flag_dsws(:,:)               ! SICE02
      real*8, allocatable :: snowe(:,:),Te(:,:),WTRe(:,:),ICEe(:,:),
     &     SNOage(:,:,:),evmax(:,:),fsat(:,:),gq(:,:)      ! EARTH01
      real*8, allocatable :: W(:,:,:,:),HT(:,:,:,:),SNWbv(:,:,:) ! SOILS03
      real*8, allocatable :: SNOW(:,:),T(:,:,:),
     &     MDWN(:,:),EDWN(:,:)
      real*8 ::  ACCPDA, ACCPDG,  EACCPDA, EACCPDG           !GLAI

c*    write
      real*8, allocatable :: Tocn_out(:,:,:,:),MixLD_out(:,:,:)   ! OCN01
      real*8, allocatable :: F_out(:,:,:),H_out(:,:,:,:),
     &     snw_out(:,:,:),msi_out(:,:,:),ssi_out(:,:,:,:),
     &     pond_melt_out(:,:,:), flag_dsws_out(:,:,:)                ! SICE02
      real*8, allocatable :: snowe_out(:,:,:),Te_out(:,:,:),
     &     WTRe_out(:,:,:), ICEe_out(:,:,:),SNOage_out(:,:,:,:),
     &     evmax_out(:,:,:),fsat_out(:,:,:),gq_out(:,:,:)         ! EARTH01
      real*8, allocatable :: W_out(:,:,:,:,:),
     &     HT_out(:,:,:,:,:),SNWbv_out(:,:,:,:)                   ! SOILS02
      real*8, allocatable :: SNOW_out(:,:,:),T_out(:,:,:,:),
     &     MDWN_out(:,:,:),EDWN_out(:,:,:)        !GLAI
      real*8, allocatable :: tsource(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      character*80 TITLEOCN01,TITLESICE02,TITLEEARTH01,TITLESOILS03,
     &     TITLEGLAIC01,name,outnc
      integer :: iu_GIC,iuout,ims,jms,nts,imt,jmt,ntt
      integer :: i,j,k,l,m,fid,status,ntiles,im,jm,d2,d3
     

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt",ims,jms,nts,imt,jmt,ntt

c      name="GIC.E046D3M20A.1DEC1955.ext"
      name="GIC.360X180.DEC01.1.rep"
c      name="GIC.144X90.DEC01.1.ext"

      iu_GIC=20

      if (am_i_root()) then
      open(iu_GIC,FILE=name,FORM='unformatted', STATUS='old')
      endif

      allocate (Tocn(3,ims,jms),MixLD(ims,jms),
     &     F(ims,jms),H(4,ims,jms),snw(ims,jms),msi(ims,jms),
     &     ssi(4,ims,jms),pond_melt(ims,jms),flag_dsws(ims,jms),
     &     snowe(ims,jms),Te(ims,jms),WTRe(ims,jms),ICEe(ims,jms), 
     &     SNOage(3,ims,jms),evmax(ims,jms),fsat(ims,jms),gq(ims,jms),
     &     W(7,3,ims,jms),HT(7,3,ims,jms),
     &     SNWbv(2,ims,jms),
     &     SNOW(ims,jms),T(2,ims,jms),
     &     MDWN(ims,jms),EDWN(ims,jms))

      allocate (Tocn_out(3,imt,jmt,ntt),MixLD_out(imt,jmt,ntt),
     &     F_out(imt,jmt,ntt),H_out(4,imt,jmt,ntt),
     &     snw_out(imt,jmt,ntt),msi_out(imt,jmt,ntt),
     &     ssi_out(4,imt,jmt,ntt),pond_melt_out(imt,jmt,ntt),
     &     flag_dsws_out(imt,jmt,ntt),
     &     snowe_out(imt,jmt,ntt),Te_out(imt,jmt,ntt),
     &     WTRe_out(imt,jmt,ntt),ICEe_out(imt,jmt,ntt), 
     &     SNOage_out(3,imt,jmt,ntt),evmax_out(imt,jmt,ntt),
     &     fsat_out(imt,jmt,ntt),gq_out(imt,jmt,ntt),
     &     W_out(7,3,imt,jmt,ntt),
     &     HT_out(7,3,imt,jmt,ntt),
     &     SNWbv_out(3,imt,jmt,ntt),
     &     SNOW_out(imt,jmt,ntt),T_out(2,imt,jmt,ntt),
     &     MDWN_out(imt,jmt,ntt),EDWN_out(imt,jmt,ntt))

      allocate (tsource(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt) )

      if (am_i_root()) then
      read(iu_GIC) TITLEOCN01, Tocn,MixLD
      write(*,*) TITLEOCN01
      read(iu_GIC) TITLESICE02, F,H,snw,msi,ssi,pond_melt,flag_dsws
      write(*,*) TITLESICE02
      read(iu_GIC) TITLEEARTH01, snowe,Te,WTRe, ICEe, SNOage,evmax,
     &     fsat,gq
      write(*,*) TITLEEARTH01
      read(iu_GIC) TITLESOILS03, W,HT,SNWbv
      write(*,*) TITLESOILS03
      read(iu_GIC) TITLEGLAIC01, SNOW,T,MDWN,EDWN,
     &     ACCPDA,ACCPDG,EACCPDA,EACCPDG
      write(*,*) TITLEGLAIC01

      close(iu_GIC)
      endif

            
      do k=1,3
         tsource(:,:,1)=Tocn(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         Tocn_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

     
      tsource(:,:,1)=MixLD(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      MixLD_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=F(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      F_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         tsource(:,:,1)=H(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         H_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=snw(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      snw_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=msi(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      msi_out(:,:,:)=ttargglob(:,:,:)

      do k=1,4
         tsource(:,:,1)=ssi(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         ssi_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=pond_melt(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      pond_melt_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1) = 0.d0

      do j=1,jms
         do i=1,ims
            if (flag_dsws(i,j) .eq. .true.) then
               tsource(i,j,1)=1.d0
            else
               tsource(i,j,1)=0.d0
            endif
         enddo
      enddo

      call root_regrid(x2grids,tsource,ttargglob)

c***  Compatibility
      do k=1,ntt
         do j=1,jmt
            do i=1,imt
               if (ttargglob(i,j,k) .ge. 0.5d0) then
                  flag_dsws_out(i,j,k)=1.
               else
                  flag_dsws_out(i,j,k)=0.
               endif
            enddo
         enddo
      enddo

      tsource(:,:,1)=snowe(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      snowe_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=Te(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      Te_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=WTRe(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      WTRe_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=ICEe(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      ICEe_out(:,:,:)=ttargglob(:,:,:)

      do k=1,3
         tsource(:,:,1)=SNOage(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         SNOage_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=evmax(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      evmax_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=fsat(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      fsat_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=gq(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      gq_out(:,:,:)=ttargglob(:,:,:)
      
      
      do k=1,7
         do m=1,3
            tsource(:,:,1)=W(k,m,:,:)
            call root_regrid(x2grids,tsource,ttargglob)
            W_out(k,m,:,:,:)=ttargglob(:,:,:)
         enddo
      enddo

      do k=1,7
         do m=1,3
            tsource(:,:,1)=HT(k,m,:,:)
            call root_regrid(x2grids,tsource,ttargglob)
            HT_out(k,m,:,:,:)=ttargglob(:,:,:)
         enddo
      enddo

      do k=1,2
         tsource(:,:,1)=SNWbv(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         SNWbv_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      SNWbv_out(3,:,:,:) = 0.

      tsource(:,:,1)=SNOW(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      SNOW_out(:,:,:)=ttargglob(:,:,:)

      do k=1,2
         tsource(:,:,1)=T(k,:,:)
         call root_regrid(x2grids,tsource,ttargglob)
         T_out(k,:,:,:)=ttargglob(:,:,:)
      enddo

      tsource(:,:,1)=MDWN(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      MDWN_out(:,:,:)=ttargglob(:,:,:)

      tsource(:,:,1)=EDWN(:,:)
      call root_regrid(x2grids,tsource,ttargglob)
      EDWN_out(:,:,:)=ttargglob(:,:,:)
      

c***  Write Netcdf file
#ifdef TRACERS_WATER
      write(*,*) "STOP TRACERS WATER NOT IMPLEMENTED IN regridinput"
      stop
#endif      
      

      deallocate (Tocn,MixLD,F,H,snw,msi,
     &     ssi,pond_melt,flag_dsws,
     &     snowe,Te,WTRe,ICEe, 
     &     SNOage,evmax,fsat,gq,
     &     W,HT,SNWbv,SNOW,T,MDWN,EDWN )

      deallocate (tsource,ttargglob)     
      
      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc
      
      if (am_i_root()) then
         status = nf_create(outnc,nf_clobber,fid)
         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      endif

    
c***  Define OCN variables
      call defvar(dd2d,fid,Tocn_out,'tocean(d3,im,jm,tile)')
      call defvar(dd2d,fid,MixLD_out,'z1o(im,jm,tile)')
c***  Define SICE variables
      call defvar(dd2d,fid,F_out,'rsi(im,jm,tile)')
      call defvar(dd2d,fid,H_out,'hsi(lmi,im,jm,tile)')
      call defvar(dd2d,fid,snw_out,'snowi(im,jm,tile)')
      call defvar(dd2d,fid,msi_out,'msi(im,jm,tile)')
      call defvar(dd2d,fid,ssi_out,'ssi(lmi,im,jm,tile)')
      call defvar(dd2d,fid,pond_melt_out,
     &     'pond_melt(im,jm,tile)')
      call defvar(dd2d,fid,flag_dsws_out,
     &     'flag_dsws(im,jm,tile)')
c***  Define EARTH variables
      call defvar(dd2d,fid,snowe_out,'snowe(im,jm,tile)')
      call defvar(dd2d,fid,Te_out,'tearth(im,jm,tile)')
      call defvar(dd2d,fid,WTRe_out,'wearth(im,jm,tile)')
      call defvar(dd2d,fid,ICEe_out,'aiearth(im,jm,tile)')
      call defvar(dd2d,fid,SNOage_out,'snoage(d3,im,jm,tile)')
      call defvar(dd2d,fid,evmax_out,
     &     'evap_max_ij(im,jm,tile)')
      call defvar(dd2d,fid,fsat_out,'fr_sat_ij(im,jm,tile)')
      call defvar(dd2d,fid,gq_out,'qg_ij(im,jm,tile)')
c***  Define SOIL variables     ! this is the old SOIL02 version, implement the new version
      call defvar(dd2d,fid,W_out,
     &     'w_ij(zero_to_ngm,ls_nfrac,im,jm,tile)')
      call defvar(dd2d,fid,HT_out,
     &     'ht_ij(zero_to_ngm,ls_nfrac,im,jm,tile)')
      call defvar(dd2d,fid,SNWbv_out,
     &     'snowbv(ls_nfrac,im,jm,tile)')
c***  Define GLAIC variables
      call defvar(dd2d,fid,SNOW_out,'snowli(im,jm,tile)')
      call defvar(dd2d,fid,T_out,'tlandi(d2,im,jm,tile)')
      call defvar(dd2d,fid,MDWN_out,'mdwnimp(im,jm,tile)')
      call defvar(dd2d,fid,EDWN_out,'edwnimp(im,jm,tile)')
      call defvar(dd2d,fid,ACCPDA,'accpda')
      call defvar(dd2d,fid,ACCPDG,'accpdg')
      call defvar(dd2d,fid,EACCPDA,'eaccpda')
      call defvar(dd2d,fid,EACCPDG,'eaccpdg')
      if (am_i_root()) then
         status = nf_enddef(fid)
         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      endif

c***  Write OCN variables
      call write_data(dd2d,fid,'tocean',Tocn_out)
      call write_data(dd2d,fid,'z1o',MixLD_out)
c***  Write SICE variables
      call write_data(dd2d,fid,'rsi',F_out)
      call write_data(dd2d,fid,'hsi',H_out)
      call write_data(dd2d,fid,'snowi',snw_out)
      call write_data(dd2d,fid,'msi',msi_out)
      call write_data(dd2d,fid,'ssi',ssi_out)
      call write_data(dd2d,fid,'pond_melt',pond_melt_out)
      call write_data(dd2d,fid,'flag_dsws',flag_dsws_out)
c***  Write EARTH variables
      call write_data(dd2d,fid,'snowe',snowe_out)
      call write_data(dd2d,fid,'tearth',Te_out)
      call write_data(dd2d,fid,'wearth',WTRe_out)
      call write_data(dd2d,fid,'aiearth',ICEe_out)
      call write_data(dd2d,fid,'snoage',SNOage_out)
      call write_data(dd2d,fid,'evap_max_ij',evmax_out)
      call write_data(dd2d,fid,'fr_sat_ij',fsat_out)
      call write_data(dd2d,fid,'qg_ij',gq_out)
c***  Write SOIL variables     ! this is the old SOIL02 version, implement the new version
      call write_data(dd2d,fid,'w_ij',W_out)
      call write_data(dd2d,fid,'ht_ij',HT_out)
      call write_data(dd2d,fid,'snowbv',SNWbv_out)
c***  Write GLAIC variables
      call write_data(dd2d,fid,'snowli',SNOW_out)
      call write_data(dd2d,fid,'tlandi',T_out)
      call write_data(dd2d,fid,'mdwnimp',MDWN_out)
      call write_data(dd2d,fid,'edwnimp',EDWN_out)
      call write_data(dd2d,fid,'accpda',ACCPDA)
      call write_data(dd2d,fid,'accpdg',ACCPDG)
      call write_data(dd2d,fid,'eaccpda',EACCPDA)
      call write_data(dd2d,fid,'eaccpdg',EACCPDG)


      deallocate (Tocn_out,MixLD_out,F_out,H_out,
     &     snw_out,msi_out,ssi_out,pond_melt_out,
     &     flag_dsws_out,snowe_out,Te_out,
     &     WTRe_out,ICEe_out,SNOage_out,evmax_out,
     &     fsat_out,gq_out,W_out,
     &     HT_out,SNWbv_out,
     &     SNOW_out,T_out,MDWN_out,EDWN_out )

      
      if(am_i_root()) status = nf_close(fid)
      
      write(*,*) "end regrid GIC"

      end subroutine regridGIC
c*



      subroutine regridAIC(x2grids,dd2d)
c
c     for 1x1 resolution : Jeff uses AIC=AIC.RES_X40.D771201N.rep
c
      use regrid_com
      use pario, only : defvar,write_data
      use dd2d_utils
      implicit none
      include 'netcdf.inc'
      type(x_2gridsroot), intent(in) :: x2grids
      type (dist_grid), intent(in) :: dd2d
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:)
      real*4, allocatable :: ts4(:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),tcopy(:,:,:),
     &     tcopy2(:,:,:),tcopy3(:,:,:),tcopy4(:,:,:)
      character*80 :: TITLE,outunformat,outnc,name
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,
     &     status,fid,vid
      integer iu_AIC
     

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "ims,jms,nts,imt,jmt,ntt r8",
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts),ts4(ims,jms,nts),
     &     ttargglob(imt,jmt,ntt),tcopy(imt,jmt,ntt),
     &     tcopy2(imt,jmt,ntt), 
     &     tcopy3(imt,jmt,ntt), 
     &     tcopy4(imt,jmt,ntt))

      tsource(:,:,:)=0.0

      iu_AIC=20

      if (am_i_root()) then
      name="AIC.RES_M20A.D771201"
      open(iu_AIC,FILE=name,FORM='unformatted', STATUS='old')

      outunformat=trim(name)//".CS"
      write(*,*) outunformat

      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      irec=1

      do
         read(unit=iu_AIC,END=30) TITLE, ts4
         tsource=ts4
         write(*,*) "TITLE, irec",TITLE,irec
         call root_regrid(x2grids,tsource,ttargglob)

         if (irec .eq. 1) tcopy=ttargglob
         if (irec .eq. 82) tcopy2=ttargglob
         if (irec .eq. 2) tcopy3=ttargglob
         if (irec .eq. 22) tcopy4=ttargglob

         write(unit=iuout) TITLE,real(ttargglob,KIND=4)
c          write(unit=iuout) TITLE,ttargglob
         irec=irec+1
      enddo
   
 30   continue

      maxrec=irec-1

      write(*,*) "maxrec",maxrec

      close(iuout) 
      close(iu_AIC)

      endif

      write(*,*) "here w r4"

      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc

      if (am_i_root()) then
         write(*,*) "TCOPY=",tcopy
         status = nf_create(outnc,nf_clobber,fid)
         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      endif

      call defvar(dd2d,fid,tcopy,'press(im,jm,tile)')
      call defvar(dd2d,fid,tcopy2,'surf_temp(im,jm,tile)')
      call defvar(dd2d,fid,tcopy3,'u974(im,jm,tile)')
      call defvar(dd2d,fid,tcopy4,'v974(im,jm,tile)')

      if (am_i_root()) then
         status = nf_enddef(fid)
         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      endif

      call write_data(dd2d,fid,'press',tcopy)
      call write_data(dd2d,fid,'surf_temp',tcopy2)
      call write_data(dd2d,fid,'u974',tcopy3)
      call write_data(dd2d,fid,'v974',tcopy4)

      if(am_i_root()) status = nf_close(fid)

      deallocate(tsource,ts4,ttargglob,tcopy,tcopy2,
     &     tcopy3,tcopy4)
 
      end subroutine regridAIC
c*


      subroutine read_recs_1R(tsource,iuin,TITLE,maxrec,im1,jm1,ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*4 :: data(im1,jm1,ntl,nrecmax)
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      write(*,*) "iuin",iuin

      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), data(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= data(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)


      end subroutine read_recs_1R
c*

      subroutine read_recs_1R_r4_r8(tsource,iuin,TITLE,maxrec,
     *     im1,jm1,ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      real*4 :: ts4(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      if (am_i_root()) then

      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), ts4(:,:,:,irec)
         tsource(:,:,:,irec)=ts4(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      endif

      end subroutine read_recs_1R_r4_r8
c*


      subroutine read_recs_1R_r8(tsource,iuin,TITLE,maxrec,im1,jm1,
     *     ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*8, intent(inout) :: tsource(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      if (am_i_root()) then
      
      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), tsource(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      endif

      end subroutine read_recs_1R_r8
c*


      subroutine read_recs_2R(tsource1,tsource2,iuin,
     &     TITLE,maxrec,im1,jm1,ntl)
      use regrid_com
      implicit none
      integer i,j,k,irec
      integer, intent(in) :: im1,jm1,ntl
      real*4 :: data1(im1,jm1,ntl,nrecmax),
     &     data2(im1,jm1,ntl,nrecmax)
      real*8, intent(inout) :: tsource1(im1,jm1,ntl,nrecmax),
     &     tsource2(im1,jm1,ntl,nrecmax)
      integer, intent(in) :: iuin
      character*80, intent(inout) :: TITLE(nrecmax)
      integer, intent(out) :: maxrec

      if (am_i_root()) then
      
      write(*,*) "iuin",iuin
      irec=1

      do
         read(unit=iuin,END=30) TITLE(irec), data1(:,:,:,irec),
     &        data2(:,:,:,irec)

         tsource1(:,:,:,irec)= data1(:,:,:,irec)
         tsource2(:,:,:,irec)= data2(:,:,:,irec)
         irec=irec+1
      enddo

 30   continue

      maxrec=irec-1

      close(iuin)

      endif

      end subroutine read_recs_2R
c*


      subroutine read_regrid_write_4D_1R(x2grids,name,iuin)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:)
      real*4, allocatable :: tout(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),
     &     tout(imt,jmt,ntt))
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat

      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")


      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         tout(:,:,:)=ttargglob(:,:,:)
         write(unit=iuout) TITLE(ir),tout(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
      enddo

      if (am_i_root()) close(iuout) 

      deallocate(tsource,ttargglob,tout)

      end subroutine read_regrid_write_4D_1R
c*


      subroutine read_regrid_write_4D_1R_r8(dd2d,x2grids,name,iuin)
      use regrid_com
      use pario, only : defvar,write_data
      use dd2d_utils
      implicit none
      include 'netcdf.inc'
      type(x_2gridsroot), intent(in) :: x2grids
      type (dist_grid), intent(in) :: dd2d
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),tcopy(:,:,:)
      character*80 :: TITLE(nrecmax),outunformat,outnc
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt,
     &     status,fid,vid

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt r8",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),tcopy(imt,jmt,ntt) )
      tsource(:,:,:,:)=0.0
      
      call read_recs_1R_r4_r8(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      if (am_i_root()) then
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      endif

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         if (ir .eq. 1) tcopy=ttargglob
         if (am_i_root()) then
         write(unit=iuout) TITLE(ir),ttargglob(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
         endif
      enddo

      if (am_i_root()) close(iuout) 

      write(*,*) "here w r8"

      outnc=trim(name)//"-CS.nc"
      write(*,*) outnc

      if (am_i_root()) then
         write(*,*) "TCOPY=",tcopy
         status = nf_create(outnc,nf_clobber,fid)
         if (status .ne. NF_NOERR) write(*,*) "UNABLE TO CREATE FILE"
      endif

      call defvar(dd2d,fid,tcopy,'press(im,jm,tile)')

      if (am_i_root()) then
         status = nf_enddef(fid)
         if (status .ne. NF_NOERR) write(*,*) "Problem with enddef"
      endif

      call write_data(dd2d,fid,'press',tcopy)

      deallocate(tsource,ttargglob,tcopy)

      end subroutine read_regrid_write_4D_1R_r8
c*



      subroutine read_regrid_write_4D_1R_rmax(x2grids,name,iuin,rmax)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      character*80, intent(in) :: name
      integer, intent(in) :: iuin
      integer, intent(in):: rmax
      integer :: iuout
      real*8, allocatable :: tsource(:,:,:,:)
      real*8, allocatable :: ttargglob(:,:,:),arrsum(:,:,:)
      real*4, allocatable :: tout(:,:,:),data(:,:,:,:)
      character*80, allocatable :: TITLE(:)
      character*80 :: outunformat
      integer :: maxrec,irec,ir,ims,jms,nts,imt,jmt,ntt

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource(ims,jms,nts,nrecmax),
     &     ttargglob(imt,jmt,ntt),arrsum(ims,jmt,ntt),
     &     tout(imt,jmt,ntt),
     &     data(ims,jms,nts,rmax),
     &     TITLE(rmax))
      tsource(:,:,:,:)=0.0
      data(:,:,:,:)=0.0
      
      if (am_i_root()) then
      do irec=1,rmax
         read(unit=iuin) TITLE(irec), data(:,:,:,irec)
         write(*,*) "TITLE, irec",TITLE(irec),irec
         tsource(:,:,:,irec)= data(:,:,:,irec)
      enddo
      maxrec=rmax
            
      write(*,*) "maxrec",maxrec
      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")

      endif

      arrsum(:,:,:)=0.

      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),ttargglob)
         arrsum(:,:,:)=arrsum(:,:,:)+ttargglob(:,:,:)
         tout(:,:,:)=ttargglob(:,:,:)
         if (am_i_root()) then
         write(unit=iuout) TITLE(ir),tout(:,:,:)
         write(*,*) "TITLE",TITLE(ir)
         endif

      enddo

      write(*,*) "SUM ARRAY=",arrsum

      if (am_i_root()) close(iuout) 

      deallocate(tsource,ttargglob,arrsum,tout,data,TITLE)

      end subroutine read_regrid_write_4D_1R_rmax
c*



      subroutine read_regrid_4D_1R(x2grids,iuin,TITLE,
     &     ttargglob,maxrec)
      use regrid_com
      implicit none
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      integer, intent(inout) :: maxrec
      real*8, allocatable :: tsource(:,:,:,:)
      real*8 :: ttargglob(x2grids%imtarget,x2grids%jmtarget,
     &     x2grids%ntilestarget,nrecmax)
      character*80, intent(inout) :: TITLE(nrecmax)
      integer :: ir,ims,jms,nts

      
      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource

      allocate (tsource(ims,jms,nts,nrecmax))
      tsource(:,:,:,:)=0.
      ttargglob(:,:,:,:)=0.

      write(*,*) "ims, jms, nts=",ims,jms,nts
            
      call read_recs_1R(tsource,iuin,TITLE,
     &        maxrec,ims,jms,nts)
      
      write(*,*) "maxrec",maxrec
      
      do ir=1,maxrec
         call root_regrid(x2grids,tsource(:,:,:,ir),
     &        ttargglob(:,:,:,ir))
c        write(*,*) "TITLE",TITLE(ir)
      enddo

      deallocate(tsource)

      end subroutine read_regrid_4D_1R
c*



      subroutine read_regrid_write_4D_2R(x2grids,name,iuin)
      use regrid_com
      type(x_2gridsroot), intent(in) :: x2grids
      integer, intent(in) :: iuin
      character*80 name
      integer :: iuout
      real*8, allocatable :: tsource1(:,:,:,:),tsource2(:,:,:,:)
      real*8, allocatable :: ttargglob1(:,:,:),ttargglob2(:,:,:)
      real*4, allocatable :: tout1(:,:,:),tout2(:,:,:),tbig(:,:,:,:)
      character*80 :: TITLE(nrecmax),outunformat
      integer :: maxrec

      ims=x2grids%imsource
      jms=x2grids%jmsource
      nts=x2grids%ntilessource
      imt=x2grids%imtarget
      jmt=x2grids%jmtarget
      ntt=x2grids%ntilestarget

      write(*,*) "iuin ims,jms,nts,imt,jmt,ntt 4D",iuin,
     &     ims,jms,nts,imt,jmt,ntt
      allocate (tsource1(ims,jms,nts,nrecmax),
     &     tsource2(ims,jms,nts,nrecmax),
     &     ttargglob1(imt,jmt,ntt),
     &     ttargglob2(imt,jmt,ntt),
     &     tout1(imt,jmt,ntt),
     &     tout2(imt,jmt,ntt),
     &     tbig(imt,jmt,2,ntt) )

      tsource1(:,:,:,:)=0.0
      tsource2(:,:,:,:)=0.0
      
      call read_recs_2R(tsource1,tsource2,iuin,TITLE,
     &     maxrec,ims,jms,nts)
      
c      write(*,*) "TITLE RECS",TITLE(:)
      write(*,*) "maxrec",maxrec

      if (am_i_root()) then      
      outunformat=trim(name)//".CS"
      
      write(*,*) outunformat
      iuout=20
      open( iuout, FILE=outunformat,
     &     FORM='unformatted', STATUS="UNKNOWN")
      endif

      do ir=1,maxrec
         call root_regrid(x2grids,tsource1(:,:,:,ir),ttargglob1)
         tout1(:,:,:)=ttargglob1(:,:,:)
         call root_regrid(x2grids,tsource2(:,:,:,ir),ttargglob2)
         tout2(:,:,:)=ttargglob2(:,:,:)
         tbig(:,:,1,:)=tout1
         tbig(:,:,2,:)=tout2
         if (am_i_root()) then
         write(unit=iuout) TITLE(ir),tbig
         write(*,*) "TITLE",TITLE(ir)
         endif
      enddo
      
      if (am_i_root()) close(iuout) 

      deallocate(tsource1,tsource2,ttargglob1,ttargglob2,tout1,tout2,
     &     tbig)

      end subroutine read_regrid_write_4D_2R
