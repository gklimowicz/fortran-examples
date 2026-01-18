#include "rundeck_opts.h"
      MODULE NUDGE_COM
!@sum  NUDGE_COM contains all the nudging related variables
!@auth Susanne Bauer/Gavin Schmidt
!@ver
      USE RESOLUTION, only : im,jm,lm
      IMPLICIT NONE
      SAVE
!@param nlevnc vertical levels of input data
#ifdef MERRA_NUDGING
      INTEGER, PARAMETER :: nlevnc =42 ! MERRA
#else
      INTEGER, PARAMETER :: nlevnc =17 ! NCEP (default)
#endif
!@var  U1, V1 NCEP wind at prior ncep timestep (m/s)
!@var  U2, V2 NCEP wind at the following ncep timestep (m/s)
      REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: u1,v1,u2,v2
      REAL*4, ALLOCATABLE, DIMENSION(:,:,:) :: UN1,VN1,UN2,VN2
      REAL*4, DIMENSION(nlevnc) :: pl
      REAL*8, DIMENSION(IM,JM,LM) :: u18,v18,u28,v28
      REAL*8, DIMENSION(IM,JM,LM) :: UN18,VN18,UN28,VN28
      REAL*8, DIMENSION(nlevnc) :: pl8
!@var netcdf integer
      INTEGER :: ncidu,ncidv,uid,vid,plid,tmid
!@var step_rea current second set timestamp
      INTEGER :: step_rea
!@var nts_max max number of time steps in file
      INTEGER :: nts_max
!@var nts_day total number of time steps in one day
      INTEGER :: nts_day
!@var nts_len length of time step in hours
      !variable is "real" to allow for time steps shorter than one hour
      REAL*8  :: nts_len
!@var tau nudging time interpolation
      REAL*8 :: tau
!@param  anudgeu anudgev relaxation constant (1/s)
      REAL*8 :: anudgeu = 0.001d0, anudgev = 0.001d0

      END MODULE NUDGE_COM

c******************************************************************

      SUBROUTINE ALLOC_NUDGE(grid)
!@sum allocate nudging related variables
      USE DOMAIN_DECOMP_1D, ONLY: DIST_GRID
      USE NUDGE_COM

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H
      INTEGER :: IER

      J_0H = grid%J_STRT_HALO
      J_1H = grid%J_STOP_HALO
      ALLOCATE ( u1(IM,J_0H:J_1H,LM),
     $           v1(IM,J_0H:J_1H,LM),
     $           u2(IM,J_0H:J_1H,LM),
     $           v2(IM,J_0H:J_1H,LM),
     $   STAT = IER)
      ALLOCATE ( un1(IM,J_0H:J_1H,nlevnc),
     $           vn1(IM,J_0H:J_1H,nlevnc),
     $           un2(IM,J_0H:J_1H,nlevnc),
     $           vn2(IM,J_0H:J_1H,nlevnc),
     $   STAT = IER)

      RETURN
      END SUBROUTINE ALLOC_NUDGE

c******************************************************************

      SUBROUTINE NUDGE_INIT
!@sum  Initialization for Nudging - called once at beginning of run
!@auth Susanne Bauer/Gavin Schmidt
!@ver
      USE DOMAIN_DECOMP_ATM, only: grid
      USE DOMAIN_DECOMP_1D, only: am_i_root, broadcast
      USE RESOLUTION, only : im,jm,lm
      use model_com, only: modelEclock
      USE MODEL_COM, only : itime,nday
      USE NUDGE_COM
      USE Dictionary_mod
      IMPLICIT NONE
      include 'netcdf.inc'
      character(len=4) :: nstr1,nstr2

C**** Rundeck parameters:
      call sync_param( "ANUDGEU", ANUDGEU )
      call sync_param( "ANUDGEV", ANUDGEV )

C**** initiallise all netcdf parameters etc.

C**** determine current model year
      write(nstr1,'(I0)') modelEclock%getYear()

C**** always need to open at least one file
      if(am_i_root()) then
        call open_nudge_file_init(nstr1)
      end if
C**** broadcast pressure levels just once (since they don't change)
      call broadcast(grid,pl8)
      pl(1:nlevnc)=sngl(pl8(1:nlevnc))

C**** broadcast timestep variables just once (since they shouldn't change either)
      call broadcast(grid,nts_max)
      call broadcast(grid,nts_day)
      call broadcast(grid,nts_len)

c**** determine current nudging file time step
      step_rea = INT((((modelEclock%getDayOfYear() - 1) * 24) +
     &  modelEclock%getHour())/nts_len) + 1  

c**** print time info to PRT file
      if (am_i_root())
     &  print*,'READING REANALYSIS INIT ',modelEclock%getHour(), 
     &    modelEclock%getDayOfYear(), step_rea

C**** read first set of nudged winds
      if (am_i_root()) print*,"nudge init0",step_rea
      call read_nudge_file(un1,vn1,step_rea)

C**** read in second set of nudged winds
      step_rea=step_rea+1
C**** if near end of year, may need to open another file
      if (step_rea.eq.nts_max+1) then
        step_rea=1
        write(nstr2,'(I0)') modelEclock%getYear() + 1
        if (am_i_root()) then
          call close_nudge_file(nstr1)
          call  open_nudge_file(nstr2)
        endif
      end if
      if (am_i_root()) print*,"nudge init1",step_rea
      call read_nudge_file(un2,vn2,step_rea)

      return
      end subroutine nudge_init

c******************************************************************

      SUBROUTINE NUDGE_PREP
!@sum  Nudging of the horizontal wind components to reanalysis
!@+    nudge_prep is called from within dynamic time step
!@auth Susanne Bauer/Gavin Schmidt
!@ver
      USE DOMAIN_DECOMP_1D, only: am_i_root
      USE RESOLUTION, only: im,jm,lm
      use model_com, only: modelEclock
      USE MODEL_COM, only: itime,nday
      USE NUDGE_COM
      IMPLICIT NONE
      include 'netcdf.inc'
      integer step_rea_1
      character(len=4) :: nstr1,nstr2

C**** Check whether nudged winds need to be updated
        step_rea_1 = INT( (((modelEclock%getDayOfYear() - 1) * 24) +
     *  modelEclock%getHour())/nts_len) + 2 

      if (step_rea.lt.step_rea_1 .and. (mod(itime,nday/nts_day).eq.0.)) 
     *     then                 ! they do

C**** move existing second set to first set
        vn1(:,:,:) = vn2(:,:,:)
        un1(:,:,:) = un2(:,:,:)
        step_rea=step_rea_1
C**** check whether new file is needed

        if (step_rea.eq.nts_max+1) then
          step_rea = 1
          write(nstr1,'(I0)') modelEclock%getYear()
          write(nstr2,'(I0)') modelEclock%getYear() + 1
          if (am_i_root()) then
            call close_nudge_file(nstr1)
            call  open_nudge_file(nstr2)
          endif
        end if
C**** read new second set of winds
        call read_nudge_file(un2,vn2,step_rea)
      endif

C**** set time interpolation for this dynamic time step
      tau  = mod(itime,nday/nts_day)/float(nday/nts_day) 

C**** vertical interpolation
      call vinterana2mod(un1,nlevnc,pl,u1)
      call vinterana2mod(vn1,nlevnc,pl,v1)
      call vinterana2mod(un2,nlevnc,pl,u2)
      call vinterana2mod(vn2,nlevnc,pl,v2)

      RETURN
      END SUBROUTINE NUDGE_PREP

c******************************************************************

      SUBROUTINE NUDGE(UGCM,VGCM,DTSTEP)
!@sum  Nudging of the horizontal wind components to reanalysis data sets
!@auth Susanne Bauer
!@ver
      USE RESOLUTION, only : im,jm,lm,plbot
      USE DOMAIN_DECOMP_ATM, only : grid
      USE NUDGE_COM, only : u1,v1,u2,v2,tau,anudgeu,anudgev,pl,nlevnc
      IMPLICIT NONE
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     UGCM, VGCM

      INTEGER i,j,l
      REAL*8  alphau,alphav,a,dtstep

      INTEGER :: J_0SG, J_1SG
C****
C**** Extract useful local domain parameters from "grid"
C****
      J_0SG = grid%J_STRT_STGR
      J_1SG = grid%J_STOP_STGR

C**** do the nudging in an implicit way:
C**** i.e. d(x)/dt = nu*(x0-x)
C****      (x^+ - x^-)/dt = nu ( x0 - x^+)
C****       x^+ = (x^- + nu*dt x0)/(1+nu*dt)
C****  i.e.
C****  x_gcm_new = ( x_gcm_old + alpha * x_reanalys ) / (1+alpha)
C****            where alpha = nu*dt 
C****
      alphau=dtstep * anudgeu   
      alphav=dtstep * anudgev   

      do l=1,lm
        if (plbot(l+1).lt.pl(nlevnc)) cycle
        do j=J_0SG,J_1SG        ! Please pay attention j starts at 2
          do i=1,im
            a=(1.-tau)*u1(i,j,l)+tau*u2(i,j,l) !time interpolation
            ugcm(i,j,l) = (ugcm(i,j,l)+ (a * alphau))/ (1+alphau) !nudging
            a=(1.-tau)*v1(i,j,l)+tau*v2(i,j,l) !time interpolation
            vgcm(i,j,l) = (vgcm(i,j,l)+ (a * alphav))/ (1+alphav) !nudging
          enddo
        enddo
      enddo

      RETURN
      END SUBROUTINE NUDGE

c******************************************************************

      subroutine vinterana2mod(varo,lmo,po,varn)
!@sum  vertical interpolation to gcm grid
!@auth Susanne Bauer
!@ver
      USE RESOLUTION, only : im,jm,lm
      USE ATM_COM, only : PMID ! Pressure at mid point of box (mb)
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, SOUTH
      USE CONSTANT, only : undef
      IMPLICIT NONE

      INTEGER lmo               ! vertical dimensions of input
      INTEGER  i,j,lo,ln,ll,llowest

      real po(lmo)              ! pressure levels in millibars of the input
      REAL*8 scratch(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lmo)
      REAL varo(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lmo) !  Variable on the old grid (input)
      REAL varn(im,grid%J_STRT_HALO:grid%J_STOP_HALO,lm) !  Variable on the new grid (output)
      real dp1,dp2
      integer J_0SG,J_1SG

      J_0SG = grid%J_STRT_STGR
      J_1SG = grid%J_STOP_STGR
      scratch = dble(varo)
      CALL HALO_UPDATE(grid,scratch,FROM=SOUTH)
      varo = sngl(scratch)
      do j= J_0SG, J_1SG        ! Please pay attention j_gcm starts at 2
        do i=1,im

          llowest=1 ! will be 1 for NCEP default

#ifdef MERRA_NUDGING
! because MERRA data has missing values "below ground", so to speak, 
! we want to make sure not to use those, so find the lowest level
! that has actual wind data: 
          loop_ll: do ll=1,lmo
            if(varo(i,j-1,ll).ne.sngl(undef)) then
              llowest=ll
              exit loop_ll
            end if
          end do loop_ll
#endif

          do ln=1,lm
            
            lo=llowest
            if (pmid(ln,i,j).ge.po(llowest))then
              varn(i,j,ln) =  varo(i,j-1,llowest)
            else if (pmid(ln,i,j).le.po(lmo)) then
              varn(i,j,ln) =  varo(i,j-1,lmo)
            else
 10           if ( (pmid(ln,i,j).le.po(lo)).and.
     &             (pmid(ln,i,j).gt.po(lo+1)) )then
                dp1 = (-pmid(ln,i,j)+po(lo))
                dp2 = (pmid(ln,i,j)-po(lo+1))

#ifdef MERRA_NUDGING
                ! Just in case. Can be removed if never happens:
                if(varo(i,j-1,lo)==sngl(undef).or.
     &          varo(i,j-1,lo+1)==sngl(undef))    
     &          call stop_model('undefined wind in vinterana2mod',13)
#endif
                varn(i,j,ln)=((varo(i,j-1,lo) * dp1)
     &               + (varo(i,j-1,lo+1) *dp2)) / (dp1 + dp2)
              else
                lo=lo+1
                if(lo==lmo)call stop_model(
     &          'lo reached lmo in NUDGING.',13)
                if (lo.le.lmo) goto 10
              end if
            endif
            
          enddo
        enddo
      enddo
      
      return
      end  subroutine vinterana2mod

c******************************************************************

      subroutine open_nudge_file_init(nstr)
!@sum open a new nudging file, and read in all necessary variables
!@+   and dimensions, including those that do not change with time
!@+   (should only be called by root)
      USE NUDGE_COM
      implicit none
      include 'netcdf.inc'
      character(len=4) :: nstr
      integer status
      !variables for time-step calculations
      !---------
      character(len=4) :: tname
      !---------
      
      print*, 'IN NUDGE: OPEN NF FILES','  {u,v}'//trim(nstr)//'.nc'
      status=NF_OPEN('u'//trim(nstr)//'.nc',NCNOWRIT,ncidu)
      if(status /= nf_noerr)call nudgeStop('opening U file',status)
      status=NF_OPEN('v'//trim(nstr)//'.nc',NCNOWRIT,ncidv)
      if(status /= nf_noerr)call nudgeStop('opening V file',status)
      
      status=NF_INQ_VARID(ncidu,'level',plid)
      if(status /= nf_noerr)call nudgeStop('finding pl var',status)
      status=NF_INQ_VARID(ncidu,'uwnd',uid)
      if(status /= nf_noerr)call nudgeStop('finding uwnd var',status)
      status=NF_INQ_VARID(ncidv,'vwnd',vid)
      if(status /= nf_noerr)call nudgeStop('finding vwnd var',status)

c**** get levels which don't change as a function of time
      status=NF_GET_VARA_REAL(ncidu,plid,1,nlevnc,pl)
      if(status /= nf_noerr)call nudgeStop('reading pl var',status)
      pl8=0.d0
      pl8(1:nlevnc)=dble(pl(1:nlevnc))
      pl(1:nlevnc)=sngl(pl8(1:nlevnc))

c**** check to see if time dimension is present in file
      status=NF_INQ_DIMID(ncidu,'time',tmid)
      if(status /= nf_noerr)call nudgeStop('reading time dim',status)

c**** determine total/max number of time steps from time dimension
      nts_max = 0
      status=NF_INQ_DIM(ncidu,tmid,tname,nts_max)
      if(status /= nf_noerr)
     &  call nudgeStop('reading time length',status)   

c**** determine if file represents a whole year
      if(mod(nts_max,365) /= 0) 
     &  call nudgeStop('Nudge time not divisible by 365:',nts_max)

c**** determine number of time steps per day
      nts_day = nts_max/365

c**** determine number of hours per time step
      nts_len = 24/nts_day

      return
      end subroutine open_nudge_file_init

c******************************************************************

      subroutine open_nudge_file(nstr)
!@sum open a new nudging file (should only be called by root)
      USE NUDGE_COM
      implicit none
      include 'netcdf.inc'
      character(len=4) :: nstr
      integer status

      print*, 'IN NUDGE: OPEN NF FILES','  {u,v}'//trim(nstr)//'.nc'
      status=NF_OPEN('u'//trim(nstr)//'.nc',NCNOWRIT,ncidu)
      if(status /= nf_noerr)call nudgeStop('opening U file',status)
      status=NF_OPEN('v'//trim(nstr)//'.nc',NCNOWRIT,ncidv)
      if(status /= nf_noerr)call nudgeStop('opening V file',status)

      status=NF_INQ_VARID(ncidu,'uwnd',uid)
      if(status /= nf_noerr)call nudgeStop('finding uwnd var',status)
      status=NF_INQ_VARID(ncidv,'vwnd',vid)
      if(status /= nf_noerr)call nudgeStop('finding vwnd var',status)

      return
      end subroutine open_nudge_file

      subroutine close_nudge_file(nstr)
!@sum close old nudging file (should only be called by root)
      USE NUDGE_COM
      implicit none
      include 'netcdf.inc'
      character(len=4) :: nstr
      integer status

      status=NF_CLOSE(ncidu)
      if(status /= nf_noerr)call nudgeStop('closing u file',status)
      status=NF_CLOSE(ncidv)
      if(status /= nf_noerr)call nudgeStop('closing v file',status)
      
      return
      end subroutine close_nudge_file

      subroutine read_nudge_file(un,vn,timestep)
!@sum return velocities for specific time step
      USE RESOLUTION, only : im,jm
      USE DOMAIN_DECOMP_ATM, only : grid
      USE DOMAIN_DECOMP_1D, only : unpack_data, am_i_root
      USE NUDGE_COM
      implicit none
      include 'netcdf.inc'
      real*4 un(IM,grid%j_strt_halo:grid%j_stop_halo,nlevnc),
     *       vn(IM,grid%j_strt_halo:grid%j_stop_halo,nlevnc)
      real*8 scratch(IM,grid%j_strt_halo:grid%j_stop_halo,nlevnc)
      real*4 scratch_glob(im,jm,nlevnc)
      integer, parameter, dimension(4) :: count= (/im, jm-1, nlevnc, 1/)
      integer start(4), timestep
      integer status

      start= (/ 1, 1, 1, timestep /)
c -----------------------------------------------------------------
c   read u, v
c -----------------------------------------------------------------
      if(am_i_root()) then
        scratch_glob=0.
        status=NF_GET_VARA_REAL(ncidu,uid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
        if(status /= nf_noerr)call nudgeStop('reading u var',status)
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      un=sngl(scratch)

      if(am_i_root()) then
        status=NF_GET_VARA_REAL(ncidv,vid,start,count,
     &                        scratch_glob(:,1:JM-1,:))
        if(status /= nf_noerr)call nudgeStop('reading v var',status)
      endif
      call unpack_data(grid, dble(scratch_glob), scratch)
      vn=sngl(scratch)

      return
      end subroutine read_nudge_file


      subroutine nudgeStop(activityString,status)
!@sum error handling for the netCDF/Fortran interface calls  
      USE DOMAIN_DECOMP_1D, only : am_i_root
      implicit none
      include 'netcdf.inc'
      character(len=*) :: activityString
      integer status

      if (am_i_root())
     & print*, 'WIND NUDGING: while model was '//trim(activityString)
     &//', encountered the NF error: ',trim(trim(nf_strerror(status)))

      call stop_model('Nudging error: see PRT for detailed message.',13)

      end subroutine nudgeStop
