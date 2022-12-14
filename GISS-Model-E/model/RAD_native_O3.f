#include "rundeck_opts.h"

      module RAD_native_O3
!@sum  Declaring a few new variables necessary for reading ozone
!@+ data for radition code at native GCM horizonal resolution.
!@auth Greg Faluvegi

      implicit none

      SAVE

      real*8, allocatable, dimension(:,:,:,:) :: O3YEAR,delta_O3_max_min
      real*8, allocatable, dimension(:,:,:,:) :: O3ICMA,O3JCMA
      real*8, allocatable, dimension(:,:,:)::O3JDAY_native,O3JREF_native

      end module RAD_native_O3


      subroutine alloc_RAD_native_O3(grid)
!@SUM  alllocate RAD_native_O3 arrays for current grid
!@auth Greg Faluvegi
      use domain_decomp_atm, only : dist_grid, getDomainBounds
      use RADPAR, only : NLO3
      use RAD_native_O3
      implicit none

      type (dist_grid), intent(in) :: grid
      integer :: J_1H, J_0H, I_0H, I_1H

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &                 I_STRT_HALO=I_0H, I_STOP_HALO=I_1H )

      allocate(           O3YEAR(I_0H:I_1H,J_0H:J_1H,NLO3,0:12) )
      allocate( delta_O3_max_min(I_0H:I_1H,J_0H:J_1H,NLO3,0:12) )
      allocate(           O3ICMA(I_0H:I_1H,J_0H:J_1H,NLO3,12) )
      allocate(           O3JCMA(I_0H:I_1H,J_0H:J_1H,NLO3,12) )
      allocate(    O3JDAY_native(NLO3,I_0H:I_1H,J_0H:J_1H) )
      allocate(    O3JREF_native(NLO3,I_0H:I_1H,J_0H:J_1H) )
      
      return
      end subroutine alloc_RAD_native_O3


      subroutine native_UPDO3D(JYEARO,JJDAYO)
!@sum A simplified version of UPDO3D to work at GCM resolution,
!@+ more input files, and with a much simpler time interpolation:
!@+ that is, no DSO3 in the stratosphere of CH4 dependence in
!@+ the troposphere [other than that inherent in the offline files].
!@auth G. Faluvegi (directly from A. Lacis/R. Ruedy)

      use filemanager, only : openunit,closeunit,nameunit
      use domain_decomp_atm, only: grid, am_i_root, readt_parallel
      use domain_decomp_atm, only: getDomainBounds
      use Dictionary_mod
      use RADPAR, only: O3yr_max,plb0,plbO3,S00WM2,RATLS0,NLO3
      use RAD_native_O3
      implicit none

! Routine expects up to 12 offline ozone files on 49 pressure levels
! and native GCM horizonal resolution (IM x JM) and no trend file.

      integer, parameter ::
     & NFO3x=12,     !  max. number of decadal ozone files used
     & IYIO3=1850,   !  year of earliest ozone file (not trend file)
     & IYEO3=2000    !  year of last ozone file (not trend file)

      real*4,dimension(NLO3):: delta_O3_now
      character*80 :: title
      logical :: qexist
!@dbparam use_sol_Ox_cycle if =1, a cycle of ozone is appled to
!@+ o3year, as a function of the solar constant cycle.
!@var add_sol is [S00WM2(now)-1/2(S00WM2min+S00WM2max)]/
!@+ [S00WM2max-S00WM2min] so that O3(altered) = O3(default) +
!@+ add_sol*delta_O3_max_min
      integer :: use_sol_Ox_cycle = 0,  NFO3 = NFO3X
      save use_sol_Ox_cycle
      real*8 :: add_sol, S0min, S0max
      character*40, dimension(NFO3X) :: DDFILE
      integer, dimension(NFO3X) :: IYEAR =
     &(/1850,1870,1890,1910,1930,1940,1950,1960,1970,1980,1990,2000/)
      integer :: ifile,idfile,I,J,L,M,N,IY,JY,JYEARX,JYEARO,JJDAYO,MI,MJ
      real*8 :: WTI,WTJ,WTMJ,WTMI,XMI

      integer, SAVE :: IYR=0, JYRNOW=0, IYRDEC=0, IFIRST=1, JYR
      save nfo3,iyear,ddfile,ifile,idfile,S0min,S0max

      INTEGER :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1)

C**** Deal with out-of-range years (incl. starts before 1850)
      if(abs(jyearo) < abs(jyrnow)) jyearo=-jyrnow ! keep cycling
      if(abs(jyearo) > O3YR_max) jyearo=-O3YR_max ! cycle thru O3YR_max

      if(IFIRST==1) then  ! --------------------- first time through ---

      call sync_param("use_sol_Ox_cycle",use_sol_Ox_cycle)

      if(use_sol_Ox_cycle==1)then
        call stop_model('Greg: address delta_O3 distributed read',255)
! I need to turn this into a parallel read, and also
! we need a native-grid file for this...
!       call openunit ("delta_O3",idfile,.true.,.true.)
!       read(idfile)S0min,S0max
!       do m=1,12; do L=1,NLO3
!         read(idfile)title,delta_O3_max_min(:,:,L,M)
!       enddo    ; enddo
!       delta_O3_max_min(:,:,:,0)=delta_O3_max_min(:,:,:,12)
!       call closeunit (idfile)
      endif

      if(plbo3(1) < plb0(1)) plbo3(1)=plb0(1)                  ! ??
C**** Find O3 data files and fill array IYEAR from title(1:4)
      nfo3=0
      if(am_i_root()) write(6,'(/a)') ' List of O3 files used:'
      loop_files: do n=1,nfo3X ! generic names O3file_01,...
        ddfile(n)=' '
        write(ddfile(n),'(a7,i2.2)') 'O3file_',n
        inquire(file=trim(ddfile(n)),exist=qexist)
        if(.not.qexist) exit loop_files
        call openunit(ddfile(n),ifile,.true.,.true.)
        read(ifile) title
        call closeunit(ifile)
        if (am_i_root())
     &  write(6,'(a,a)') ' read O3 file, record 1: ',title
        read(title(1:4),*) IYEAR(n)
        nfo3=nfo3+1
      end do loop_files

      if(nfo3==1) JYEARO=-IYEAR(1)
      if(nfo3==0) call stop_model('native_updo3d: no Ozone files',255)

C**** Prior to first year of data, cycle through first year of data
      if(abs(jyearo) < IYEAR(1)) jyearo=-IYEAR(1)

      IY=0
      if(JYEARO < 0) then  ! 1 year of O3 data is used in a cyclical way
        do n=1,nfo3        ! check for file matching date
           if(IYEAR(n)==-JYEARO) IY=n
        end do
      end if

      if(IY > 0) then
        call openunit(ddfile(IY),ifile,.true.,.true.)
        do M=1,12 ; do L=1,NLO3
          call readt_parallel
     &    (grid,ifile,nameunit(ifile),O3YEAR(:,:,L,M),1)
        enddo ; enddo
        O3YEAR(:,:,:,0)=O3YEAR(:,:,:,12)
        call closeunit(ifile)
        JYRNOW=-JYEARO  ! insures that O3YEAR is no longer changed
      end if

#ifdef TRACERS_SPECIAL_Shindell
! read the 3D field for O3 RCOMPX reference calls
      call openunit ('Ox_ref',ifile,.true.,.true.)
      do L=1,NLO3
        call readt_parallel
     &  (grid,ifile,nameunit(ifile),O3JREF_native(L,:,:),1)
      enddo
      call closeunit(ifile)
#endif
      IFIRST=0
      endif           ! ---------------end of first time through ---


C     To time input data READs, JYEARX is set ahead of JYEARO by 15 days
C     ------------------------------------------------------------------
      if(JYEARO < 0) then     !        ... except in cyclical case
        JYEARX=-JYEARO        ! Use fixed-year decadal climatology
      else
        JYEARX=min(JYEARO+(JJDAYO+15)/366,IYEO3+1) ! +1 for continuity
      end if                                       !         at Dec 15

! Next section skipped of O3YEAR already ready for O3JDAY-defining
! (including cyclical case past O3YR_max):
      IF(JYEARX.ne.JYRNOW .and. JYRNOW<=O3YR_max)THEN

C****
C**** Get 13 months of O3 data O3YEAR starting with the leading December
C****
      do jy=1,nfo3                  ! nfo3 is at least 2, if we get here
        if(iyear(jy) > JYEARx) go to 100
      end do
      jy=nfo3
  100 if(jy <= 1) jy=2
      iy=jy-1
      IYR=IYEAR(IY)
      JYR=IYEAR(JY)

C**** Get first decadal file
      call openunit(ddfile(IY),ifile,.true.,.true.)               ! IYR
      do M=1,12 ; do L=1,NLO3
        call readt_parallel
     &  (grid,ifile,nameunit(ifile),O3ICMA(:,:,L,M),1)
      enddo ; enddo
      call closeunit(ifile)

C**** Define December 
!     Read and use prior decadal file to define prior year December
!     (only when starting up with JYEARO=1870,...,1990 AND non-cyclical)
      if(JYEARX == IYR.and.IYRDEC.ne.JYEARX-1 .and. IY > 1 .and.
     & JYEARO > 0.and.JYEARX <= O3YR_max) then 

        call openunit(ddfile(IY-1),ifile,.true.,.true.)          ! KYR
        do M=1,12 ; do L=1,NLO3
          call readt_parallel
     &    (grid,ifile,nameunit(ifile),O3JCMA(:,:,L,M),1)
        enddo ; enddo
        call closeunit(ifile)

        ! get simple linear interpolation coeffs:
        call native_O3_WTS(IYIO3, IYR, IYEAR(IY-1), JYEARX-1,! in
     &                     WTI, WTJ)                         ! out

        do J=J_0,J_1 ; do L=1,NLO3 ; do I=I_0,I_1     
          O3YEAR(I,J,L,0)=
     &    max(0.d0, WTI*O3ICMA(I,J,L,12)+WTJ*O3JCMA(I,J,L,12) )
        enddo; enddo; enddo

        IYRDEC=JYEARX   ! Set flag to indicate December data is current
      endif 

C**** Get next decadal file
      call openunit(ddfile(JY),ifile,.true.,.true.)               ! JYR
      do M=1,12 ; do L=1,NLO3
        call readt_parallel
     &  (grid,ifile,nameunit(ifile),O3JCMA(:,:,L,M),1)
      enddo ; enddo
      call closeunit (ifile)

      if(JYEARX.ne.IYRDEC) then ! we are not done with prior December

        if(JYEARX==IYRDEC+1) then  ! copy data from M=12 -> M=0
          O3YEAR(:,:,:,0)=O3YEAR(:,:,:,12) ! DEC from prior year
          if(JYEARX > O3YR_max) JYEARX=O3YR_max
          IYRDEC=JYEARX  ! Set flag to indicate December data is current
        else if(JYEARO > IYEAR(1)) then
C         Interpolate prior December from the decadal files - start-up
          ! get simple linear interpolation coeffs:
          call native_O3_WTS(IYIO3, IYR, JYR, JYEARX-1, ! in
     &                       WTI, WTJ)                  ! out

          do J=J_0,J_1 ; do L=1,NLO3 ; do I=I_0,I_1     
            O3YEAR(I,J,L,0)=
     &      max(0.d0, WTI*O3ICMA(I,J,L,12)+WTJ*O3JCMA(I,J,L,12) )
          enddo; enddo; enddo

          IYRDEC=JYEARX  ! Set flag to indicate December data is current
        endif 
      
      end if ! now we are done with prior December


C            Fill in a full year of O3 data by interpolation
C            -----------------------------------------------
      if(JYEARX > O3YR_max) JYEARX=O3YR_max
      ! get simple linear interpolation coeffs:
      call native_O3_WTS(IYIO3, IYR, JYR, JYEARX, ! in
     &                     WTI, WTJ)              ! out

      do M=1,12 ; do J=J_0,J_1 ; do L=1,NLO3 ; do I=I_0,I_1 
        O3YEAR(I,J,L,M)=
     &  max(0.d0, WTI*O3ICMA(I,J,L,M)+WTJ*O3JCMA(I,J,L,M) )
      enddo ; enddo ; enddo ; enddo

      if(JYEARX.ne.IYRDEC) then               ! cyclical start-up case
        O3YEAR(:,:,:,0)=O3YEAR(:,:,:,12)      ! DEC from current year
        IYRDEC=JYEARX  ! Set flag to indicate December data is current
      endif   
      JYRNOW=JYEARX

      ENDIF ! O3YEAR defining was necessary 

C****
C**** O3JDAY_native is interp daily from O3YEAR seasonal data via JJDAYO
C****
C     the formula below yields M near the middle of month M
      XMI=(JJDAYO+JJDAYO+31-(JJDAYO+15)/61+(JJDAYO+14)/61)/61.D0
      MI=XMI
      WTMJ=XMI-MI       !   Intra-year interpolation is linear in JJDAYO
      WTMI=1.D0-WTMJ
      IF(MI > 11) MI=0
      MJ=MI+1
      if(use_sol_Ox_cycle==1)then
        add_sol = (S00WM2*RATLS0-0.5d0*(S0min+S0max))/(S0max-S0min)
        write(6,661)JJDAYO,S00WM2*RATLS0,S0min,S0max,add_sol
      endif
  661 format('JJDAYO,S00WM2*RATLS0,S0min,S0max,frac=',I4,3F9.2,F7.3)
      do J=J_0,J_1 ; do  I=I_0,I_1 
        O3JDAY_native(:,I,J)=WTMI*O3YEAR(I,J,:,MI)+WTMJ*O3YEAR(I,J,:,MJ)
        if(use_sol_Ox_cycle==1) then
          delta_O3_now(:) = WTMI*delta_O3_max_min(I,J,:,MI) +
     &                      WTMJ*delta_O3_max_min(I,J,:,MJ)
          O3JDAY_native(:,I,J) = O3JDAY_native(:,I,J) + 
     &                           add_sol*delta_O3_now(:)
        endif
      enddo ; enddo
      return
      end subroutine native_UPDO3D


      subroutine native_O3_wts (IYI, IY1, IY2, IYX,        !  input
     *                          WTI, WTJ )                 ! output
!@sum this just does simplest possible weighting, so
!@+ hardly even based on O3_WTS.  
!@auth G. Faluvegi
      implicit none
      integer :: IYI     ! first year                       - input
      integer :: IY1,IY2 ! 2 distinct years with O3 data    - input
      integer :: IYX     ! current year                     - input
      real*8 :: WTI,WTJ  ! weights                          - output
      real*8 :: dyear

      dyear=IY2-IY1
      WTI=min(1.d0,max(0.d0,(IY2-IYX)/dyear))
      WTJ=1.d0-WTI

      return
      end subroutine native_O3_WTS


