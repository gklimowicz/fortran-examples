C**** QUSEM12 E001M12 SOMTQ QUSB261AM12
C**** QUSBM9=QUSB140M9 with correction to second order moment calc.
C**** Changes for constant pressure above LS1 + REAL*8
C**** QUS   is Russell quadratic upstream scheme for temperature
C**** and water vapor advection, with limits applied to water vapor.
C**** Changes for constant pressure above LS1
C**** FQU,FQV for additional diagnostics
C**** Routines included: AADVT, AADVTX, AADVTY, AADVTZ

      MODULE QUSCOM
!@sum  QUSCOM contains gcm-specific advection parameters/workspace
!@auth Maxwell Kelley
      USE QUSDEF
      IMPLICIT NONE
      SAVE
      INTEGER :: IM,JM,LM
      INTEGER :: XSTRIDE,YSTRIDE,ZSTRIDE
      REAL*8 :: BYIM
C**** AIR MASS FLUXES
      REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: MFLX
C**** WORKSPACE FOR AADVTX
cc    REAL*8, DIMENSION(:), ALLOCATABLE :: AM,F_I
cc    REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMOM_I
C**** WORKSPACE FOR AADVTY
cc    REAL*8, DIMENSION(:), ALLOCATABLE :: BM,F_J
cc    REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMOM_J
C**** WORKSPACE FOR AADVTZ
cc    REAL*8, DIMENSION(:), ALLOCATABLE :: CM,F_L
cc    REAL*8, DIMENSION(:,:), ALLOCATABLE :: FMOM_L
      END MODULE QUSCOM

      SUBROUTINE init_QUS(grd_dum,IM_GCM,JM_GCM,LM_GCM)
!@sum  init_QUS sets gcm-specific advection parameters/workspace
!@auth Maxwell Kelley
      USE DOMAIN_DECOMP_1D, only : DIST_GRID, getDomainBounds
      use QUSCOM
      USE Dictionary_mod
      INTEGER, INTENT(IN) :: IM_GCM,JM_GCM,LM_GCM
      TYPE (DIST_GRID), INTENT(IN) :: grd_dum

C****
C**** Extract local domain parameters from "grd_dum"
C****
      INTEGER J_0,J_1,J_0H,J_1H
      call getDomainBounds(grd_dum, J_STRT=J_0,       J_STOP=J_1,
     &                  J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)

C**** SET RESOLUTION
      IM = IM_GCM
      JM = JM_GCM
      LM = LM_GCM
      BYIM = 1.D0/REAL(IM,KIND=8)
      XSTRIDE = 1
      YSTRIDE = IM
      ZSTRIDE = IM*(J_1H-J_0H+1)
C**** ALLOCATE SPACE FOR AIR MASS FLUXES
      ALLOCATE(MFLX(IM,J_0H:J_1H,LM))
C**** ALLOCATE WORKSPACE FOR AADVTX
cc    ALLOCATE(AM(IM),F_I(IM),FMOM_I(NMOM,IM))
C**** ALLOCATE WORKSPACE FOR AADVTY
cc    ALLOCATE(BM(JM),F_J(J_0H:J_1H),FMOM_J(NMOM,J_0H:J_1H))
C**** ALLOCATE WORKSPACE FOR AADVTZ
cc    ALLOCATE(CM(LM),F_L(LM),FMOM_L(NMOM,LM))

      call sync_param("prather_limits",prather_limits)

      RETURN
      END SUBROUTINE init_QUS

      Subroutine AADVT (DT, MM,RM,RMOM, QLIMIT, FQU,FQV)
!@sum  AADVT advection driver
!@auth G. Russell, modified by Maxwell Kelley
c****
c**** AADVT advects tracers using the Quadradic Upstream Scheme.
c****
c**** input:
c****  MU,MV,MW (kg/s) = east-west,north-south,vertical mass fluxes
c****      qlimit = whether moment limitations should be used
C****         DT (s) = time step
c****
c**** input/output:
c****     rm = tracer concentration
c****   rmom = moments of tracer concentration
c****     MM (kg) = fluid mass
c****
      Use DYNAMICS, Only: MU,MV,MW
      USE DOMAIN_DECOMP_ATM, only: grid, getDomainBounds
      USE DOMAIN_DECOMP_1D, only: HALO_UPDATE, NORTH,SOUTH
      USE QUSDEF
      USE QUSCOM, ONLY : IM,JM,LM, MFLX
      IMPLICIT NONE

      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: 
     &                  MM,RM
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) 
     &               :: rmom
      REAL*8, INTENT(IN) :: DT
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO), 
     & intent(inout) :: fqu,fqv

      INTEGER :: I,J,L,N
      REAL*8 :: BYMA
c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &               J_STRT_SKP  = J_0S,   J_STOP_SKP  = J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

C**** Initialise diagnostics
      FQU=0.  ; FQV=0.

C**** Fill in values at the poles
C**** SOUTH POLE:
      if (HAVE_SOUTH_POLE) then
        DO L=1,LM
          DO I=2,IM
            RM(I,1 ,L) =   RM(1,1 ,L)
            DO N=1,NMOM
              RMOM(N,I,1 ,L) =  RMOM(N,1,1 ,L)
            enddo
          enddo
        enddo
      end if       !SOUTH POLE

c**** NORTH POLE:
      if (HAVE_NORTH_POLE) then
        DO L=1,LM
          DO I=2,IM
            RM(I,JM,L) =   RM(1,JM,L)
            DO N=1,NMOM
              RMOM(N,I,JM,L) =  RMOM(N,1,JM,L)
            enddo
          enddo
        enddo
      end if    ! NORTH POLE
C****
C**** convert from concentration to mass units
C****
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
         RM    (I,J,L) = RM    (I,J,L)*MM(I,J,L)
         RMOM(:,I,J,L) = RMOM(:,I,J,L)*MM(I,J,L)
      enddo
      enddo
      enddo
C****
C**** Advect the tracer using the quadratic upstream scheme
C****
CC    mflx(:,:,:)=MU(:,:,:)*(.5*dt)
       DO L=1,LM
          MFLX(:,:,L) = MU(:,:,L)*(.5*DT)
       ENDDO
      Call AADVTX (RM,RMOM,MM,MFLX,QLIMIT,FQU)

CC    mflx(:,1:jm-1,:)=MV(:,2:jm,:)*dt
CC    mflx(:,jm,:)=0.
!     Call HALO_UPDATE (GRID, MV, From=NORTH)   Haloed in AFLUX
       DO L=1,LM
          MFLX(:,J_0:J_1S,L) = MV(:,J_0+1:J_1S+1,L)*DT
       ENDDO
       if (HAVE_NORTH_POLE) mflx(:,jm,:)=0.

      Call AADVTY (RM,RMOM,MM,MFLX,QLIMIT,FQV)
CC    mflx(:,:,1:lm-1)=MW(:,:,1:lm-1)*(-dt)
CC    mflx(:,:,lm)=0.
      DO L=1,LM
         IF(L.NE.LM)  THEN
            MFLX(:,:,L) = MW(:,:,L)*(-DT)
         ELSE
            MFLX(:,:,L)=0.
         END IF
      ENDDO
      Call AADVTZ (RM,RMOM,MM,MFLX,QLIMIT)
CC    mflx(:,:,:)=MU(:,:,:)*(.5*dt)
       DO L=1,LM
          MFLX(:,:,L) = MU(:,:,L)*(.5*DT)
       ENDDO
      Call AADVTX (RM,RMOM,MM,MFLX,QLIMIT,FQU)
C****
C**** convert from mass to concentration units
C****
      DO L=1,LM
      DO J=J_0,J_1
      DO I=1,IM
         byMA = 1 / MM(I,J,L)
         RM(I,J,L)=RM(I,J,L)*BYMA
         RMOM(:,I,J,L)=RMOM(:,I,J,L)*BYMA
      enddo
      enddo
      enddo
      RETURN
      END

      subroutine aadvtx(rm,rmom,mass,mu,qlimit,fqu)
!@sum  AADVTX advection driver for x-direction
!@auth Maxwell Kelley
c****
c**** aadvtx advects tracers in the west to east direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mu (kg) = west-east mass flux, positive eastward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      USE DOMAIN_DECOMP_ATM, only: grid
      use DOMAIN_DECOMP_1D, only : getDomainBounds
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, xstride,am,f_i,fmom_i
      use QUSCOM, only : im,jm,lm, xstride
      implicit none
      REAL*8, dimension(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: 
     &                  rm,mass,mu,hfqu
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &                  rmom
      logical ::  qlimit
      REAL*8, INTENT(OUT), 
     &        DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) :: FQU
      REAL*8  AM(IM), F_I(IM), FMOM_I(NMOM,IM)
     &     ,MASS_I(IM), COURMAX, BYNSTEP
      integer :: i,ip1,j,l,ierr,nerr,ICKERR,ns,nstep

c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1, J_0S, J_1S
      call getDomainBounds(grid, J_STRT = J_0 , J_STOP=J_1,
     &             J_STRT_SKP=J_0S,J_STOP_SKP=J_1S )
c**** loop over layers and latitudes
      ICKERR=0
      do l=1,lm
      do j=J_0S,J_1S
c****
c**** decide how many timesteps to take
c****
      nstep=0
      courmax = 2.
      do while(courmax.gt.1. .and. nstep.lt.20)
        nstep = nstep+1
        bynstep = 1d0/real(nstep,kind=8)
        am(:) = mu(:,j,l)*bynstep
        mass_i(:)  = mass(:,j,l)
        courmax = 0.
        do ns=1,nstep
          i = im
          do ip1=1,im
            if(am(i).gt.0.) then
               courmax = max(courmax,+am(i)/mass_i(i))
            else
               courmax = max(courmax,-am(i)/mass_i(ip1))
            endif
            i = ip1
          enddo
          if(ns.lt.nstep) then
             i = im
             do ip1=1,im
                mass_i(ip1) = mass_i(ip1) + (am(i)-am(ip1))
                i = ip1
             enddo
          endif
        enddo
      enddo
      if(courmax.gt.1.) then
         write(6,*) 'aadvtx: j,l,courmax=',j,l,courmax
         ICKERR=ICKERR+1
      endif

c      am(:) = mu(:,j,l)*bynstep ! am already set
      hfqu(:,j,l)  = 0.
c****
c**** loop over timesteps
c****
      do ns=1,nstep
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(1,j,l),rmom(1,1,j,l), f_i,fmom_i, mass(1,j,l),
     &        am, im, qlimit,xstride,xdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtx: i,j,l=",nerr,j,l
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(a) > 1"
CCC       call stop_model('Error in qlimit: abs(a) > 1',11)
          ICKERR=ICKERR+1
        end if
      end if
c****
c**** store tracer flux in fqu array
c****
CCC   fqu(:,j)  = fqu(:,j) + f_i(:)
      hfqu(:,j,l)  = hfqu(:,j,l) + f_i(:)
      enddo ! ns
      enddo ! j
      enddo ! l
c
c     now sum into fqu
c
      do j=J_0S,J_1S
      do l=1,lm
         fqu(:,j)  = fqu(:,j) + hfqu(:,j,l)
      enddo ! j
      enddo ! l
C
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in aadvtx',11)
C
      return
c****
      end subroutine aadvtx

      subroutine aadvty(rm,rmom,mass,mv,qlimit,fqv)
!@sum  AADVTY advection driver for y-direction
!@auth Maxwell Kelley
c****
c**** aadvty advects tracers in the south to north direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg) = north-south mass flux, positive northward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      USE DOMAIN_DECOMP_ATM, only: grid
      use DOMAIN_DECOMP_1D, only : getDomainBounds, halo_update
      use DOMAIN_DECOMP_1D, only : halo_update_column
      use DOMAIN_DECOMP_1D, only : NORTH, SOUTH, AM_I_ROOT
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, ystride,bm,f_j,fmom_j, byim
      use QUSCOM, only : im,jm,lm, ystride,               byim
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) :: 
     &                  rm,mass,mv
      REAL*8, dimension(NMOM,IM,grid%J_STRT_HALO:
     &                          grid%J_STOP_HALO,LM) :: rmom
      logical ::  qlimit
      REAL*8, intent(out), dimension(im,grid%J_STRT_HALO:
     &                                  grid%J_STOP_HALO) :: fqv
      REAL*8      BM(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     &                F_J(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM),
     &        FMOM_J(NMOM,IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM)
      integer :: i,j,l,ierr,ICKERR, err_loc(3)
      REAL*8, DIMENSION(LM) ::
     &     m_sp,m_np,rm_sp,rm_np,rzm_sp,rzm_np,rzzm_sp,rzzm_np


c****Get relevant local distributed parameters
      INTEGER J_0,J_1,J_0H,J_1H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
      call getDomainBounds(grid, J_STRT = J_0,
     &               J_STOP = J_1,
     &               J_STRT_HALO = J_0H,
     &               J_STOP_HALO = J_1H,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)

c**** loop over layers
      ICKERR=0
      fqv = 0
      do l=1,lm

c**** scale polar boxes to their full extent
      if (HAVE_SOUTH_POLE) then
        mass(:,1,l)=mass(:,1,l)*im
        m_sp(l) = mass(1,1 ,l)
        rm(:,1,l)=rm(:,1,l)*im
        rm_sp(l) = rm(1,1 ,l)
        do i=1,im
           rmom(:,i,1 ,l)=rmom(:,i,1 ,l)*im
        enddo
        rzm_sp(l)  = rmom(mz ,1,1 ,l)
        rzzm_sp(l) = rmom(mzz,1,1 ,l)
      end if                       !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm,l)=mass(:,jm,l)*im
        m_np(l) = mass(1,jm,l)
        rm(:,jm,l)=rm(:,jm,l)*im
        rm_np(l) = rm(1,jm,l)
        do i=1,im
           rmom(:,i,jm,l)=rmom(:,i,jm,l)*im
        enddo
        rzm_np(l)  = rmom(mz ,1,jm,l)
        rzzm_np(l) = rmom(mzz,1,jm,l)
      end if                       !NORTH POLE

c***c**** loop over longitudes
c***      do i=1,im
c****
c**** load 1-dimensional arrays
c****
      bm   (:,:,l) = mv(:,:,l) !/nstep

c**** POLES
      IF (HAVE_SOUTH_POLE)  rmom(ihmoms,:,1,l ) = 0. ! horizontal moments are zero at pole
      IF (HAVE_NORTH_POLE) THEN
        bm(:,jm,l)= 0.
        rmom(ihmoms,:,jm,l) = 0.
      END IF
      end do

c****
c**** call 1-d advection routine
c****
        call advection_1D_custom( rm(1,j_0h,1), rmom(1,1,j_0h,1), 
     &     f_j(1,j_0h,1),fmom_j(1,1,j_0h,1), mass(1,j_0h,1),
     &     bm(1,j_0h,1),j_1-j_0+1,LM,(/ (.true.,l=1,lm) /),
     &     qlimit,ystride,ydir,ierr,err_loc)

        if (ierr.gt.0) then
          write(6,*) "Error in aadvty: i,j,l=",err_loc
          if (ierr.eq.2) then
            write(0,*) "Error in qlimit: abs(b) > 1"
ccc         call stop_model('Error in qlimit: abs(b) > 1',11)
            ICKERR=ICKERR+1
          endif
        end if
! horizontal moments are zero at pole
        IF (HAVE_SOUTH_POLE) rmom(ihmoms,:,1, :) = 0
        IF (HAVE_NORTH_POLE) rmom(ihmoms,:,jm,:) = 0.


      do l=1,lm

c**** average and unscale polar boxes
      if (HAVE_SOUTH_POLE) then
        mass(:,1 ,l) = (m_sp(l) + sum(mass(:,1 ,l)-m_sp(l)))*byim
        rm(:,1 ,l) = (rm_sp(l) + sum(rm(:,1 ,l)-rm_sp(l)))*byim
        rmom(mz ,:,1 ,l) = (rzm_sp(l) + 
     *       sum(rmom(mz ,:,1 ,l)-rzm_sp(l) ))*byim
        rmom(mzz,:,1 ,l) = (rzzm_sp(l)+ 
     *       sum(rmom(mzz,:,1 ,l)-rzzm_sp(l)))*byim
      end if   !SOUTH POLE

      if (HAVE_NORTH_POLE) then
        mass(:,jm,l) = (m_np(l) + sum(mass(:,jm,l)-m_np(l)))*byim
        rm(:,jm,l) = (rm_np(l) + sum(rm(:,jm,l)-rm_np(l)))*byim
        rmom(mz ,:,jm,l) = (rzm_np(l) + 
     &       sum(rmom(mz ,:,jm,l)-rzm_np(l) ))*byim
        rmom(mzz,:,jm,l) = (rzzm_np(l)+ 
     &       sum(rmom(mzz,:,jm,l)-rzzm_np(l)))*byim
      end if  !NORTH POLE

      enddo ! end loop over levels
c
c     sum into fqv
c
      do j=J_0,J_1
        do l=1,lm
          fqv(:,j)  = fqv(:,j) + f_j(:,j,l)
        enddo                   ! l
      enddo                     ! j
      if (HAVE_NORTH_POLE) fqv(:,jm) = 0. ! not really needed
C
      IF(ICKERR.GT.0)  CALL stop_model('Stopped in aadvty',11)
C
      return
c****
      end subroutine aadvty


      subroutine aadvtz(rm,rmom,mass,mw,qlimit)
!@sum  AADVTZ advection driver for z-direction
!@auth Maxwell Kelley
c****
c**** aadvtz advects tracers in the upward vertical direction using the
c**** quadratic upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mw (kg) = vertical mass flux, positive upward
c****      qlimit = whether moment limitations should be used
c****
c**** input/output:
c****     rm (kg) = tracer mass
c****   rmom (kg) = moments of tracer mass
c****   mass (kg) = fluid mass
c****
      USE DOMAIN_DECOMP_ATM, only: grid
      use DOMAIN_DECOMP_1D, only : getDomainBounds
      use QUSDEF
ccc   use QUSCOM, only : im,jm,lm, zstride,cm,f_l,fmom_l
      use QUSCOM, only : im,jm,lm, zstride
      implicit none
      REAL*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lm) 
     &        :: rm,mass,mw
      REAL*8, dimension(NMOM,IM,grid%j_strt_halo:grid%j_stop_halo,LM) 
     &        :: rmom
      logical ::  qlimit
      REAL*8  CM(LM),F_L(LM),FMOM_L(NMOM,LM),MASS_L(LM)
      real*8 bynstep,courmax
      integer :: i,j,l,ierr,nerr,ICKERR,ns,nstep
c**** Get useful local parameters for domain decomposition
      integer :: J_0, J_1
      call getDomainBounds( grid, J_STRT=J_0 , J_STOP=J_1 )
c**** loop over latitudes and longitudes
      ICKERR=0
      do j=J_0,J_1
      do i=1,im
c****
c**** decide how many timesteps to take
c****
      nstep=0
      courmax = 2.
      do while(courmax.gt.1. .and. nstep.lt.20)
        nstep = nstep+1
        bynstep = 1d0/real(nstep,kind=8)
        cm(:) = mw(i,j,:)*bynstep
        cm(lm) = 0.
        mass_l(:)  = mass(i,j,:)
        courmax = 0.
        do ns=1,nstep
          do l=1,lm-1
            if(cm(l).gt.0.) then
               courmax = max(courmax,+cm(l)/mass_l(l))
            else
               courmax = max(courmax,-cm(l)/mass_l(l+1))
            endif
          enddo
          if(ns.lt.nstep) then
            do l=1,lm-1
               mass_l(l  ) = mass_l(l  ) - cm(l)
               mass_l(l+1) = mass_l(l+1) + cm(l)
            enddo
          endif
        enddo
      enddo
      if(courmax.gt.1.) then
         write(6,*) 'aadvtz: i,j,courmax=',i,j,courmax
         ICKERR=ICKERR+1
      endif

c      cm(:) = mw(i,j,:)*bynstep ! cm already set
c      cm(lm)= 0.

c****
c**** loop over timesteps
c****
      do ns=1,nstep
c****
c**** call 1-d advection routine
c****
      call adv1d(rm(i,j,1),rmom(1,i,j,1),f_l,fmom_l,mass(i,j,1),
     &        cm,lm,qlimit,zstride,zdir,ierr,nerr)
      if (ierr.gt.0) then
        write(6,*) "Error in aadvtz: i,j,l=",i,j,nerr
        if (ierr.eq.2) then
          write(0,*) "Error in qlimit: abs(c) > 1"
ccc       call stop_model('Error in qlimit: abs(c) > 1',11)
          ICKERR=ICKERR+1
        endif
      end if
      enddo ! ns
      enddo ! i
      enddo ! j
C
      IF(ICKERR.GT.0) call stop_model('Stopped in aadvtz',11)
      return
c****
      end subroutine aadvtz
