C****  
C**** OCNGISS_SM.f GISS ocean sub-mesoscale vertical mixing scheme 2013/02/13
C****
#include "rundeck_opts.h"

      MODULE GISS_SM_COM
!@sum  GISS_SM_COM holds variables related to GISS sub-mesoscale scheme
!@auth Cheng/Kelley/Howard/Leboissetier/Etienne
      USE OCEAN, only : im,jm,lmo
      USE SW2OCEAN, only : lsrpd
      IMPLICIT NONE
      SAVE
!@var au_sm,av_sm store time-averaged velocity components
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: au_sm,av_sm
!@var rx,ry 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: rx,ry
!@var gx,gy 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: gx,gy
!@var sx,sy 
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: sx,sy
!@var rx_sm,ry_sm store time-averaged drhodz,drhody
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: rx_sm,ry_sm
!@var gx_sm,gy_sm store time-averaged dGdz,dGdy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: gx_sm,gy_sm
!@var sx_sm,sy_sm store time-averaged dSdz,dSdy
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: sx_sm,sy_sm
!@var fvb_sm store submesoscale vertical buoyancy flux
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: fvb_sm
!@var p3d stores the 3-d pressure
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: p3d
!@var rho3d stores the 3-d density
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:):: rho3d

      END MODULE GISS_SM_COM


      MODULE GISS_SM
!@sum GISS_SM contains variables/routines for GISS sub-mesoscale scheme
!@ref Canuto et al. 2013 
!@auth Cheng/Kelley/Howard/Leboissetier/Etienne
!@date Feb 2013

      USE GISS_SM_COM
      USE CONSTANT, only : omega,by3,grav,pi

      implicit none

      CONTAINS

      subroutine giss_sm_init(iniOCEAN)
!@sum creates initializes variables for giss sm mixing
!@+   gissmix_init is called from routine init_OCEAN in OCNDYN.f
!@auth Cheng/Kelley/Howard/Leboissetier/Etienne
!@date Feb 2013

      USE FILEMANAGER
      USE DOMAIN_DECOMP_1D, only : READT_PARALLEL
      USE OCEANR_DIM, only : ogrid

      implicit none
      ! in:
      logical, intent(in) :: iniOCEAN
      ! local:
      
      fvb_sm=0.d0
      p3d=0.d0
      rho3d=0.d0

      end subroutine giss_sm_init

      END MODULE GISS_SM


      subroutine giss_sm_mix( 
      ! in:
     &    ze,zg,db,u,v,rx,ry,gx,gy,sx,sy,grav,fc,lmix,lmij
      ! out:
     &   ,fvg,fvs,fvb,p4uv) 

!@sum giss sub-mesoscale model for ocean
!@ref Canuto et al. 2013
!@auth Cheng/Kelley/Howard/Leboissetier/Etienne
!@date Feb 2013

      USE GISS_SM, only : lmo
      USE CONSTANT, only : rt2
      !@var grav acceleration of gravity m/s^2

      implicit none

      ! in:
      real*8 ze(0:lmo)   !@var ze vertical layer-edge depth (m), > 0
      real*8 zg(0:lmo+1) !@var zg vertical grid depth (m), < 0
      real*8 db(lmo)     !@var db -grav/rho*d(rho) (m/s^2) ! vertical
      real*8 u(lmo)      !@var u x-component of vel at mid-layer
      real*8 v(lmo)      !@var v x-component of vel at mid-layer
      real*8 rx(lmo)    !@var rx x-directive of rho
      real*8 ry(lmo)    !@var ry y-directive of rho
      real*8 gx(lmo)    !@var gx x-directive of g
      real*8 gy(lmo)    !@var gy y-directive of g
      real*8 sx(lmo)    !@var sx x-directive of s
      real*8 sy(lmo)    !@var sy y-directive of s
      real*8 grav       !@var grav gravit. accel. on earth (m/s**2)
      real*8 fc         !@var fc Coriolis parameter=2*omega*sin(lat) (1/s)
      integer lmix      !@var lmix edge layer the SM to mix
      integer lmij      !@var lmij layer number in the column
      intent (in) ze,zg,db,u,v,rx,ry,gx,gy,sx,sy,grav,fc,lmix,lmij
      ! out:
      real*8 fvg(0:lmo+1) !@var fvg non-local enthalpy flux ( )
      real*8 fvs(0:lmo+1) !@var fvs non-local salinity flux (m/s)
      real*8 fvb(0:lmo+1) !@var fvb non-local buoyancy flux ( )
      real*8 p4uv(lmo,2)   !@var fuv SM new terms for U,V eqns 
      intent (out) fvg,fvs,fvb,p4uv

      ! local:
      integer l
      integer lmixe ! at the layer edge
      real*8, parameter :: ce=.07d0,gc=5.d0/21.d0,rho0=1d3
      real*8 zeta
      real*8 z(lmo)   ! -ze
      real*8 ui(lmo)  ! ui(l) integral of u from ze(0) to ze(i)
      real*8 vi(lmo)  ! vi(l) integral of v from ze(0) to ze(i)
      real*8 dudz(lmo),dvdz(lmo)
      real*8 zu(lmo),zv(lmo),zui(lmo),zvi(lmo)
      real*8 kx,ky    ! x,y-components of SM diffusivity kappa
      real*8 bx,by    ! x,y-components of buoyancy 
      real*8 Rb,Ro,ome,ome2,afc,sfc,A,A1,A2,A20,At
      real*8 sgmt,AtFc
      real*8 zul,zvl,b0,um,vm
      real*8 h,h2,byh,byhg ! h is interface depth closest to dbl (m)
      real*8 zl,bydz,dudzl,dvdzl,tmp
      real*8 zum,zvm,r3um,r3vm
      real*8 fc2,dbl,bv2l,bxl,byl,ril
      logical, parameter :: GISS=.true.! if not, use Fox-Kemper model

      fvg=0.d0; fvs=0.d0; fvb=0.d0; p4uv=0.d0;
      lmixe=min(lmix,lmij-1)
      h=ze(lmixe)
      byh=1.d0/h
      byhg=1.d0/(-zg(lmix))
      h2=h*h
      b0=-grav/rho0
      fc2=fc*fc
      afc=abs(fc)
      sfc=fc/afc

      if(.not.GISS) then ! Fox-Kemper SM model
        Rb=5.d0
        A20=4*ce*h2/afc*Rb
      else ! GISS model
        Rb=1.d0
        Ro=5.d0
        ome=Ro/rt2
        ome2=ome*ome
        A=ome2/(1+ome2)
        At=ome2/(1+ome2)
        AtFc=At*fc
        A1=A*sfc/ome
        A20=A*h2*Rb/(2*ome*afc)
        ! integrals of u,v over ze from 0 to ze which equal
        ! the minus of integrals of u,v over z from 0 to z
        ! u,v at layer middle, ze (>=0) at layer edge
        ui(1)=ze(1)*u(1)
        vi(1)=ze(1)*v(1)
        do l=2,lmixe
          ui(l)=ui(l-1)+(ze(l)-ze(l-1))*u(l)
          vi(l)=vi(l-1)+(ze(l)-ze(l-1))*v(l)
        end do
        um=ui(lmixe)*byh
        vm=vi(lmixe)*byh
        ! integrals of u,v over ze from 0 to h which equal
        ! the integrals of u,v over z from -h to 0
      endif
      do l=1,lmixe
        if(db(l).gt.0.d0) then
          ! calc Fv at layer edge, l from 1 to lmixe
          bx=.5d0*b0*(rx(l)+rx(l+1))
          by=.5d0*b0*(ry(l)+ry(l+1))
          zl=-ze(l)
          zeta=zl*byh  ! < 0
          if(.not.GISS) then ! Fox-Kemper SM model
            A2=A20*zeta*(1+zeta)*(1+gc*(1+2*zeta)**2)
            kx=A2*bx
            ky=A2*by
          else ! GISS SM model
            A2=A20*zeta*(1+zeta)
            zul=-ui(l)-zl*um
            zvl=-vi(l)-zl*vm
            kx=A*zul+A1*zvl+A2*bx
            ky=A*zvl-A1*zul+A2*by
          endif
          fvg(l)=-.5d0*(kx*(gx(l)+gx(l+1))+ky*(gy(l)+gy(l+1)))
          fvs(l)=-.5d0*(kx*(sx(l)+sx(l+1))+ky*(sy(l)+sy(l+1)))
          fvb(l)=-(kx*bx+ky*by)
          fvg(l)=max(fvg(l),0.d0)
          fvs(l)=max(fvs(l),0.d0)
          fvb(l)=max(fvb(l),0.d0)
        endif
      end do
      if(GISS) then
        ! calc p4uv at layer middle, l from 1 to lmix
        do l=2,lmix-1
          dudz(l)=(u(l-1)-u(l+1))/(zg(l-1)-zg(l+1))
          dvdz(l)=(v(l-1)-v(l+1))/(zg(l-1)-zg(l+1))
        end do
        dudz(1)=(u(1)-u(2))/(zg(1)-zg(2))
        dvdz(1)=(v(1)-v(2))/(zg(1)-zg(2))
        dudz(lmix)=(u(lmix-1)-u(lmix))/(zg(lmix-1)-zg(lmix))
        dvdz(lmix)=(v(lmix-1)-v(lmix))/(zg(lmix-1)-zg(lmix))
        ui(1)=-zg(1)*u(1)
        vi(1)=-zg(1)*v(1)
        do l=2,lmix
          tmp=(zg(l-1)-zg(l))*.5d0
          ui(l)=ui(l-1)+tmp*(u(l-1)+u(l))
          vi(l)=vi(l-1)+tmp*(v(l-1)+v(l))
        end do
        um=ui(lmix)*byhg
        vm=vi(lmix)*byhg
        do l=1,lmix
          if(l.eq.1) then
            dbl=db(1)
          elseif(l.eq.lmij) then
            dbl=db(lmij-1)
          else
            dbl=.5d0*(db(l-1)+db(l))
          endif
          if(dbl.gt.0.d0) then
            ! for dU/dt,dV/dt equations, GISS SM model:
            p4uv(l,1)=-AtFc*(v(l)-vm)            ! R12 u-component
            p4uv(l,2)= AtFc*(u(l)-um)            ! R12 v-component
          endif
        end do
      endif
      return
      end subroutine giss_sm_mix


#ifdef OCN_GISS_SM /*alloc_giss_sm_com,def_rsf_giss_sm,new_io_giss_sm*/

      SUBROUTINE alloc_giss_sm_com(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time
!@+    this subroutine is called near the end of
!@+    subroutine alloc_ocean, the latter is in OCEAN_COM.f
!@+    remember to put #ifdef OCN_GISS_SM around the call 
!@auth Ruedy/Kelley/Cheng
!@date Sep 21,2012

      USE DOMAIN_DECOMP_1D, only : dist_grid,getDomainBounds
      USE GISS_SM_COM

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, L
      INTEGER :: IER

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      ALLOCATE(   au_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   av_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   rx(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   ry(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   gx(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   gy(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   sx(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   sy(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   rx_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   ry_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   gx_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   gy_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   sx_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   sy_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   fvb_sm(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   p3d(IM,J_0H:J_1H,LMO), STAT = IER)
      ALLOCATE(   rho3d(IM,J_0H:J_1H,LMO), STAT = IER)
      au_sm=0.d0; av_sm=0.d0;
      rx=0; ry=0; gx=0; gy=0.; sx=0; sy=0; 
      rx_sm=0; ry_sm=0; gx_sm=0; gy_sm=0; sx_sm=0; sy_sm=0.;
      fvb_sm=0.d0; p3d=0.; rho3d=0.

      END SUBROUTINE alloc_giss_sm_com

      subroutine def_rsf_giss_sm(fid)
!@sum
!@+    this subroutine is called at the end of
!@+    subroutine def_rsf_ocean, the latter is in OCNDYN.f
!@+    remember to put #ifdef OCN_GISS_SM around the call 
!@auth Ruedy/Kelley/Cheng
!@date Sep 21,2012
      use pario, only : defvar
      use oceanr_dim, only : grid=>ogrid
      use giss_sm_com, only : au_sm,av_sm
     &                       ,rx_sm,ry_sm,gx_sm,gy_sm,sx_sm,sy_sm

      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,au_sm,'au_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,av_sm,'av_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,rx_sm,'rx_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,ry_sm,'ry_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gx_sm,'gx_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,gy_sm,'gy_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,sx_sm,'sx_sm(dist_imo,dist_jmo,lmo)')
      call defvar(grid,fid,sy_sm,'sy_sm(dist_imo,dist_jmo,lmo)')
      return
      end subroutine def_rsf_giss_sm

      subroutine new_io_giss_sm(fid,iaction)
!@sum
!@+    this subroutine is called at the end of
!@+    subroutine new_io_ocean, the latter is in OCNDYN.f
!@+    remember to put #ifdef OCN_GISS_SM around the call 
!@auth Ruedy/Kelley/Cheng
!@date Sep 21,2012
      use model_com, only : ioread,iowrite
      use pario, only : read_dist_data,write_dist_data
     &                 ,read_data,write_data
      use oceanr_dim, only : grid=>ogrid
      use giss_sm_com, only : au_sm,av_sm
     &                       ,rx_sm,ry_sm,gx_sm,gy_sm,sx_sm,sy_sm
   
      implicit none
      integer fid     !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'au_sm',au_sm)
        call write_dist_data(grid,fid,'av_sm',av_sm)
        call write_dist_data(grid,fid,'rx_sm',rx_sm)
        call write_dist_data(grid,fid,'ry_sm',ry_sm)
        call write_dist_data(grid,fid,'gx_sm',gx_sm)
        call write_dist_data(grid,fid,'gy_sm',gy_sm)
        call write_dist_data(grid,fid,'sx_sm',sx_sm)
        call write_dist_data(grid,fid,'sy_sm',sy_sm)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'au_sm',au_sm)
        call read_dist_data(grid,fid,'av_sm',av_sm)
        call read_dist_data(grid,fid,'rx_sm',rx_sm)
        call read_dist_data(grid,fid,'ry_sm',ry_sm)
        call read_dist_data(grid,fid,'gx_sm',gx_sm)
        call read_dist_data(grid,fid,'gy_sm',gy_sm)
        call read_dist_data(grid,fid,'sx_sm',sx_sm)
        call read_dist_data(grid,fid,'sy_sm',sy_sm)
      end select
      return
      end subroutine new_io_giss_sm

#endif /*alloc_giss_sm_com,def_rsf_giss_sm,new_io_giss_sm*/
