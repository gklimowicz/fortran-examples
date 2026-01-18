#include "rundeck_opts.h" 
      MODULE AEROSOL_SOURCES
!@sum repository for Koch aerosol sources, features, etc.
!@auth Dorothy Koch
!@ subroutines in this file include:
!@ alloc_aerosol_sources
!@ get_O3_offline
!@ read_mon3Dsources
!@ READ_OFFHNO3
!@ read_DMS_sources
!@ aerosol_gas_chem
!@ SCALERAD
!@ GET_SULFATE
!@ GET_BC_DALBEDO
!@ GRAINS
!@ read_seawifs_chla
      IMPLICIT NONE
      SAVE
      INTEGER, PARAMETER :: ndmssrc  = 1
!@var DMSinput           DMS ocean source (kg/s/m2)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: DMSinput ! DMSinput(im,jm,12)
#ifndef TRACERS_AEROSOLS_SOA
!@var OCT_src    OC Terpene source (kg/s/box)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: OCT_src !(im,jm,12)
#endif  /* TRACERS_AEROSOLS_SOA */
!@var SO2_src_3D SO2 volcanic sources (and biomass) (kg/s)
      INTEGER :: nso2src_3d=0,iso2volcano=0,iso2volcanoexpl=0
      INTEGER :: iso2exvolc=0
      real*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: SO2_src_3D !(im,jm,lm,nso2src_3d)
!@var H2O_src_3D H2O volcanic sources (kg kg-1 s-1)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: H2O_src_3D !(im,jm,lm)
!@var PBLH boundary layer height
!@var MDF is the mass of the downdraft flux
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: 
     *   oh,dho2,perj,tno3,o3_offline  !im,jm,lm
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: ohr,dho2r,perjr,
     *   tno3r,ohsr  !im,jm,lm,12   DMK jmon
      integer :: JmonthCache = -1
      real*8, allocatable, dimension(:,:,:) :: 
     &     ohrCache, dho2rCache, perjrCache, tno3rCache
      
#ifdef BC_ALB
      real*8, ALLOCATABLE, DIMENSION(:,:) :: snosiz
#endif  /* BC_ALB */
#ifdef TRACERS_RADON
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: rn_src
#endif
!var off_HNO3 off-line HNO3 field, used for nitrate and AMP when gas phase chemistry turned off
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:)     ::  off_HNO3
#ifdef TRACERS_AEROSOLS_VBS
!@var VBSemifact factor that distributes organic aerosols in volatility bins
      real*8, allocatable, dimension(:) :: VBSemifact
#endif /* TRACERS_AEROSOLS_VBS */

!@dbparam tune_DMS Multiplication factor for DMS emissions
      real*8 :: tune_DMS=1.

      END MODULE AEROSOL_SOURCES

      SUBROUTINE alloc_aerosol_sources(grid)
!@auth D. Koch
      use domain_decomp_atm, only: dist_grid, getDomainBounds
      use TRACER_COM, only: NTM
      use TRACER_COM, only: ex_volc_num
      use AEROSOL_SOURCES, only: DMSinput,
#ifndef TRACERS_AEROSOLS_SOA
     * OCT_src,
#endif  /* TRACERS_AEROSOLS_SOA */
     * nso2src_3d,SO2_src_3D,iso2volcano,iso2volcanoexpl,H2O_src_3d,
     * iso2exvolc,
     * ohr,dho2r,perjr, tno3r, 
     * ohrCache, dho2rCache, perjrCache, tno3rCache,
     * oh,dho2,perj,tno3,ohsr
     * ,o3_offline
     * ,off_HNO3
#ifdef TRACERS_RADON
     * ,rn_src
#endif
#ifdef TRACERS_AEROSOLS_VBS
     * ,VBSemifact
      use TRACERS_VBS, only: vbs_tr
#endif
#ifdef BC_ALB
      use AEROSOL_SOURCES, only: snosiz
#endif  /* BC_ALB */
      use filemanager, only: file_exists

      use RESOLUTION, only: lm
      
      IMPLICIT NONE
      type (dist_grid), intent(in) :: grid
      integer ::  J_1H, J_0H, I_0H, I_1H
      integer :: IER
      logical :: init = .false.

      if(init)return
      init=.true.

      call getDomainBounds( grid , J_STRT_HALO=J_0H, J_STOP_HALO=J_1H )
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      allocate( DMSinput(I_0H:I_1H,J_0H:J_1H,12) ,STAT=IER)
#ifndef TRACERS_AEROSOLS_SOA
      allocate( OCT_src(I_0H:I_1H,J_0H:J_1H,12) ,STAT=IER)
#endif  /* TRACERS_AEROSOLS_SOA */
      if (file_exists('SO2_VOLCANO')) then
        nso2src_3d=nso2src_3d+1
        iso2volcano=nso2src_3d
      endif
      if (file_exists('SO2_VOLCANO_EXPL')) then
        nso2src_3d=nso2src_3d+1
        iso2volcanoexpl=nso2src_3d
      endif
      if (ex_volc_num>0) then
        nso2src_3d=nso2src_3d+1
        iso2exvolc=nso2src_3d
      endif
      allocate( SO2_src_3D(I_0H:I_1H,J_0H:J_1H,lm,nso2src_3d),STAT=IER )
      allocate( H2O_src_3D(I_0H:I_1H,J_0H:J_1H,lm),STAT=IER )
      allocate( oh(I_0H:I_1H,J_0H:J_1H,lm),dho2(I_0H:I_1H,J_0H:J_1H,lm),
     * perj(I_0H:I_1H,J_0H:J_1H,lm),tno3(I_0H:I_1H,J_0H:J_1H,lm)
     * ,o3_offline(I_0H:I_1H,J_0H:J_1H,lm),STAT=IER )
      allocate( ohr(I_0H:I_1H,J_0H:J_1H,lm),
     * dho2r(I_0H:I_1H,J_0H:J_1H,lm),
     * perjr(I_0H:I_1H,J_0H:J_1H,lm),tno3r(I_0H:I_1H,J_0H:J_1H,lm),
     * ohsr(I_0H:I_1H,J_0H:J_1H,lm),STAT=IER )
      allocate( ohrCache(I_0H:I_1H,J_0H:J_1H,lm),
     * dho2rCache(I_0H:I_1H,J_0H:J_1H,lm),
     * perjrCache(I_0H:I_1H,J_0H:J_1H,lm),
     *     tno3rCache(I_0H:I_1H,J_0H:J_1H,lm))
#ifdef BC_ALB
      allocate( snosiz(I_0H:I_1H,J_0H:J_1H) ,STAT=IER)
#endif  /* BC_ALB */
#ifdef TRACERS_RADON
      allocate( rn_src(I_0H:I_1H,J_0H:J_1H,12) ,STAT=IER)
#endif
c off line 
      allocate(  off_HNO3(I_0H:I_1H,J_0H:J_1H,LM)     )
#ifdef TRACERS_AEROSOLS_VBS
      allocate(VBSemifact(vbs_tr%nbins))
#endif

      return
      end SUBROUTINE alloc_aerosol_sources      
      SUBROUTINE get_O3_offline
!@sum read in ozone fields for aqueous oxidation
c
C**** GLOBAL parameters and variables:
C
      use resolution, only: lm
      use model_com, only: modelEclock
      use filemanager, only: openunit,closeunit
      use aerosol_sources, only: o3_offline
      use domain_decomp_atm, only: grid, getDomainBounds, write_parallel
C     
      implicit none

C**** Local parameters and variables and arguments:
c
!@var nmons: number of monthly input files
      integer, parameter :: nmons=1,levo3=23
      integer, dimension(nmons) :: mon_units,imon
      integer :: i,j,k,l
      integer :: jdlast=0
      logical :: ifirst=.true.
      character*80 title
      character(len=300) :: out_line
      character*10 :: mon_files(nmons) = (/'O3_FIELD'/)
      logical :: mon_bins(nmons)=(/.true./) ! binary file?
      real*8 :: frac
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,levo3,1):: src
      real*8, allocatable, dimension(:,:,:,:) :: tlca, tlcb
      save jdlast,mon_units,imon,ifirst,tlca,tlcb
      INTEGER :: J_1, J_0, J_0H, J_1H, I_0H, I_1H

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

      if (ifirst) then
        call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
        I_0H = grid%I_STRT_HALO
        I_1H = grid%I_STOP_HALO

        allocate(tlca(i_0H:i_1H,j_0H:j_1H,levo3,nmons),
     &           tlcb(i_0H:i_1H,j_0H:j_1H,levo3,nmons))
        ifirst=.false.
      endif
      k=1
      call openunit(mon_files(k),mon_units(k),mon_bins(k),.true.)
      call read_mon3Dsources(levo3,mon_units(k),jdlast,
     & tlca(:,:,:,k),tlcb(:,:,:,k),src(:,:,:,k),frac,imon(k))     
      call closeunit(mon_units(k))
      
      jdlast = modelEclock%getDayOfYear()

      if(levo3 /= LM) then ! might be ok, but you should check
        write(out_line,*)
     &  'make sure levels are correct in get_O3_offline'
        call write_parallel(trim(out_line))
        call stop_model('check on get_O3_offline',255)
      endif
      
      o3_offline(:,J_0:J_1,:)=src(:,J_0:J_1,:,k)

      write(out_line,*)
     &'offline ozone interpolated to current day',frac
      call write_parallel(trim(out_line))
      
      return
      END SUBROUTINE get_O3_offline
      
      
      SUBROUTINE read_mon3Dsources(Ldim,iu,jdlast,tlca,tlcb,data1,
     & frac,imon)
! we need to combine this with the one that is used in TRACERS_
! SPECIAL_Shindell, (called read_monthly_3Dsources) and then move it
! into more general tracer code...
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,Ldim,nm),tlcb(im,jm,Ldim,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others / Greg Faluvegi
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only: idofm=>JDmidOfM
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR,INT_MONTHS_PER_YEAR
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,READT_PARALLEL,
     &     REWIND_PARALLEL,write_parallel
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer Ldim,L,imon,iu,jdlast
      character(len=300) :: out_line
      real*8 :: frac
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::A2D,B2D,dummy
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::tlca,tlcb,data1

      integer :: J_0, J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
C
      if (jdlast == 0) then   ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1                ! imon=January
        if (modelEclock%getDayOfYear() <= 16)  then ! DAYOFYEAR in Jan 1-15, first month is Dec
          do L=1,Ldim*11
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          CALL REWIND_PARALLEL( iu )
        else              ! DAYOFYEAR is in Jan 16 to Dec 16, get first month
  120     imon=imon+1
          if (modelEclock%getDayOfYear() > idofm(imon) .AND. 
     *      imon <= INT_MONTHS_PER_YEAR) go to 120
          do L=1,Ldim*(imon-2)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),dummy,1)
          end do
          DO L=1,Ldim
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),A2D,1)
            tlca(:,J_0:J_1,L)=A2D(:,J_0:J_1)
          END DO
          if (imon == 13)  CALL REWIND_PARALLEL( iu )
        end if
      else                         ! Do we need to read in second month?
        if (modelEclock%getDayOfYear() /= jdlast+1) then ! Check that data is read in daily
          if (modelEclock%getDayOfYear() /= 1 .OR. 
     &      jdlast /= INT_DAYS_PER_YEAR) then
            write(out_line,*)'Bad day values in read_monthly_3Dsources'
     &      //': JDAY,JDLAST=',modelEclock%getDayOfYear(),JDLAST
            call write_parallel(trim(out_line),crit=.true.)
            call stop_model('Bad values in read_monthly_3Dsources',255)
          end if
          imon=imon-INT_MONTHS_PER_YEAR             ! New year
          go to 130
        end if
        if (modelEclock%getDayOfYear() <= idofm(imon)) go to 130
        imon=imon+1                ! read in new month of data
        if (imon == 13) then
          CALL REWIND_PARALLEL( iu  )
        else
          do L=1,Ldim*(imon-1)
            CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),dummy,1)
          end do
        endif
      end if
      DO L=1,Ldim
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),B2D,1)
        tlcb(:,J_0:J_1,L)=B2D(:,J_0:J_1)
      END DO
 130  continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-modelEclock%getDayOfYear())
     &     / (idofm(imon)-idofm(imon-1))
      data1(:,J_0:J_1,:) =
     &     tlca(:,J_0:J_1,:)*frac + tlcb(:,J_0:J_1,:)*(1.-frac)
      return
      end SUBROUTINE read_mon3Dsources

      SUBROUTINE READ_OFFHNO3(OUT)
      use resolution, only: lm
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only : JDendOFM
      USE DOMAIN_DECOMP_ATM, only : grid,am_i_root
      IMPLICIT NONE
      include 'netcdf.inc'
      REAL*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),intent(out) :: OUT
!@param  nlevnc vertical levels of off-line data  - 4x5 model=23, 2x2.5 model = 40
      INTEGER, PARAMETER :: nlevnc =40
      REAL*4, DIMENSION(GRID%I_STRT:GRID%I_STOP,
     &                  GRID%J_STRT:GRID%J_STOP,nlevnc) ::
     &     IN1_nohalo, IN2_nohalo
      REAL*8, DIMENSION(:,:,:), pointer, save :: IN1, IN2
!@var netcdf integer
      INTEGER :: ncid,id
      INTEGER, save :: step_rea=0, first_call=1
!@var time interpoltation
      REAL*8 :: tau
      integer start(4),count(4),status,l
      integer :: i_0,i_1,j_0,j_1
c -----------------------------------------------------------------
c   Initialisation of the files to be read
c ----------------------------------------------------------------     

      if (first_call==1) then
        first_call=0
        allocate( IN1(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *       ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc) )
        allocate( IN2(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *       ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,nlevnc) )
      endif
      if (step_rea.ne.modelEclock%getMonth()) then 
        step_rea = modelEclock%getMonth()
        if ( am_i_root() ) then
        print*,'READING HNO3 OFFLINE ',modelEclock%getMonth(), step_rea
        endif
c -----------------------------------------------------------------
c   Opening of the files to be read
c -----------------------------------------------------------------
        status=NF_OPEN('OFFLINE_HNO3.nc',NCNOWRIT,ncid)
        status=NF_INQ_VARID(ncid,'FELD',id)
C------------------------------------------------------------------
c -----------------------------------------------------------------
c   read
c   this is still latlon-specific.
c   will call read_dist_data for cubed sphere compatibility
c -----------------------------------------------------------------
        i_0 = grid%i_strt
        i_1 = grid%i_stop
        j_0 = grid%j_strt
        j_1 = grid%j_stop
        start(1)=i_0
        start(2)=j_0
        start(3)=1
        start(4)=step_rea
        count(1)=1+(i_1-i_0)
        count(2)=1+(j_1-j_0)
        count(3)=nlevnc
        count(4)=1

        status=NF_GET_VARA_REAL(ncid,id,start,count,IN1_nohalo)
        start(4)=step_rea+1
        if (start(4).gt.12) start(4)=1
        status=NF_GET_VARA_REAL(ncid,id,start,count,IN2_nohalo)

        status=NF_CLOSE(ncid)

        IN1(I_0:I_1,J_0:J_1,:) = IN1_nohalo(I_0:I_1,J_0:J_1,:)
        IN2(I_0:I_1,J_0:J_1,:) = IN2_nohalo(I_0:I_1,J_0:J_1,:)

      endif
C-----------------------------------------------------------------
      tau = 
     & (modelEclock%getDate()-.5)/(JDendOFM(modelEclock%getMonth()) - 
     &     JDendOFM(modelEclock%getMonth()-1))
         do l=1,lm
         OUT(:,:,l) = (1.-tau)*IN1(:,:,l)+tau*IN2(:,:,l)  
         enddo
c -----------------------------------------------------------------
      RETURN
      END SUBROUTINE READ_OFFHNO3
c -----------------------------------------------------------------

      SUBROUTINE read_DMS_sources(swind,itype,i,j,DMS_flux) !!! T
!@sum generates DMS ocean source
!@auth Koch
c Monthly DMS ocean concentration sources are read in and combined
c  with wind and ocean temperature functions to get DMS air surface
c  concentrations
c want kg DMS/m2/s
      use TimeConstants_mod, only: SECONDS_PER_DAY
      USE GEOM, only: axyp
      use OldTracer_mod, only: tr_mm
      USE TRACER_COM, only: n_DMS
      use model_com, only: modelEclock
      USE AEROSOL_SOURCES, only: tune_DMS
      USE AEROSOL_SOURCES, only: DMSinput
#ifdef old_DMS_emis
      USE FLUXES, only: GTEMP
#endif
      implicit none
      integer jread
      REAL*8 akw,erate,SCH,SCHR !!! T,Tc
#ifdef old_DMS_emis
      real*8 Tc ! YHL - FOR another DMS source
#endif
      real*8, PARAMETER :: E1=0.17d0
      real*8, PARAMETER :: E2=2.85d0
      real*8, PARAMETER :: E3=0.612d0
      real*8, PARAMETER :: E4=5.9d0
      real*8, PARAMETER :: E5=26.79d0
      real*8, PARAMETER :: E6=0.612d0
      real*8, PARAMETER :: SCHT=600.d0
      real*8, INTENT(OUT) :: DMS_flux
      real*8, INTENT(IN) :: swind
      integer, INTENT(IN) :: itype,i,j

      DMS_flux=0.d0
        erate=0.d0
        if (itype.eq.1) then
c       if (lm.lt.40) then 
#ifndef old_DMS_emis
c Nightingale et al
        akw = 0.23d0*swind*swind + 0.1d0 * swind
        akw = akw * 0.24d0
        erate=akw*DMSinput(i,j,modelEclock%getMonth())*1.d-9*62.d0 !*tr_mm(nt)
     *       /SECONDS_PER_DAY*tune_DMS
#endif

#ifdef old_DMS_emis
c YUNHA - Liss and Merlivat (1986) code is from GISS GCM II-prime. 
c Liss and Merlivat (1986), use for > lm=40 to moderate DMS flux

       Tc=GTEMP(1,1,I,J) ! YUNHA GTEMP is already Celcius. 

       SCH=2674.d0-147.12d0*Tc+3.726d0*Tc*Tc-0.038d0*Tc*Tc*Tc
       IF(Tc.gt.47.) print*,'BAD_TEMPERATURE_DMS',i,j,Tc,SCH
       SCHR=SCHT/SCH
       if (swind.lt.3.6) then
        akw=0.041*(SCHR)**(2.d0/3.d0)*swind
       else if (swind.lt.13.) then
         akw=(0.68*SWIND - 2.31)*DSQRT(SCHR)
       else
       akw=(1.42*SWIND - 11.8)*DSQRT(SCHR)
       endif  !swind
       erate=akw*DMSinput(i,j,modelEclock%month())*1.d-9*62.d0/
     *      SECONDS_PER_DAY*tune_DMS     !not sure of units

#endif
c       if (lm.ge.40) erate=erate/5.d0   !I think there was an error in input files
c       else
c Liss and Merlivat (1986), use for > lm=40 to moderate DMS flux
c       SCH=2674.d0-147.12d0*Tc+3.726d0*Tc*Tc-0.038d0*Tc*Tc*Tc
c       SCHR=SCHT/SCH
c       if (swind.lt.3.6) then
c       akw=E1*(SCHR)**(2.d0/3.d0)*swind
c       else if (swind.lt.13.) then
c       akw=E2*DSQRT(SCHR)*(swind-3.6d0)
c    *     +E3*(SCHR)**(2.d0/3.d0)
c       else
c       akw=E4*(swind-13.d0)*DSQRT(SCHR)+E5*(swind-3.6d0)*
c    *      DSQRT(SCHR)+E6*(SCHR)**(2.d0/3.d0)
c       endif  !swind
c       erate=akw*DMSinput(i,j,jmon)*1.d-9/sday*tune_DMS !not sure of units
c       endif ! lm
        endif !itype
        DMS_flux=erate          ! units are kg/m2/s
c
      return
      end SUBROUTINE read_DMS_sources

      SUBROUTINE aerosol_gas_chem
!@sum aerosol gas phase chemistry
!@vers 2013/03/27
!@auth Dorothy Koch
      use OldTracer_mod, only: trname, tr_mm
      use TRACER_COM, only: ntm, oh_live, no3_live, trm
      use TRACER_COM, only: coupled_chem, n_BCIA, n_BCII, n_DMS,n_H2O2_s
      use TRACER_COM, only: rsulf1, rsulf2, rsulf3, rsulf4
      use TRACER_COM, only: n_MSA, N_OCII, n_OX, n_SO2, n_OCIA
      use TRACER_COM, only: n_SO4, n_SO4_d1, n_SO4_d2, n_SO4_d3
#ifdef TRACERS_AEROSOLS_VBS
      use TRACER_COM, only: n_BCB, n_isopp1a, n_isopp2a, n_apinp1a,
     &                      n_apinp2a, n_NH4, n_NO3p
#endif  /* TRACERS_AEROSOLS_VBS */
      use TRACER_COM, only: nChemistry, nChemLoss, nOther
#if (defined TRACERS_HETCHEM) || (defined TRACERS_NITRATE)
      use TRACER_COM, only: rxts1, rxts2, rxts3
#endif
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AECOB, n_AECIL, n_AOCOB, n_AOCIL
      use TRACER_COM, only: n_AECIL, n_H2SO4, nbins
#endif
#ifdef TRACERS_AMP
      use TRACER_COM, only: n_H2SO4
#endif
      USE TRDIAG_COM, only : 
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
     *     jls_OHconk,jls_HO2con,jls_NO3,jls_phot
#endif
#ifdef TRACERS_SPECIAL_Shindell
     &     ,jls_OHcon
#endif
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT, getDomainBounds 
      USE DOMAIN_DECOMP_ATM, only: DREAD8_PARALLEL,DREAD_PARALLEL
      USE DOMAIN_DECOMP_ATM, only : GRID, write_parallel
      use resolution, only: im,lm
      use atm_com, only : t,q
      use model_com, only: modelEclock
      USE MODEL_COM, only: dtsrc
      USE ATM_COM, only: pmid,MA,pk,LTROPO,byMA
      USE PBLCOM, only : dclev
      USE GEOM, only: axyp,imaxj,BYAXYP
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      USE AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,
     & dho2,perj,tno3,ohsr,o3_offline, JmonthCache,
     &      ohrCache, dho2rCache, perjrCache, tno3rCache
      USE CONSTANT, only : mair
      use TimeConstants_mod, only: SECONDS_PER_DAY
#ifdef TRACERS_TOMAS
      USE TOMAS_AEROSOL, only : h2so4_chem
#endif
#ifdef TRACERS_AEROSOLS_VBS
      use CONSTANT, only : gasc
      use TRACERS_VBS, only: vbs_tracers, vbs_conditions, 
     &                       vbs_calc, vbs_tr
#endif /* TRACERS_AEROSOLS_VBS */
c Aerosol chemistry
      implicit none
      logical :: ifirst=.true.
      real*8 ppres,te,tt,mm,dmm,ohmc,r1,d1,r2,d2,ttno3,r3,d3,
     * ddno3,dddms,ddno3a,fmom,dtt
      real*8 rk4,ek4,r4,d4
      real*8 r6,d6,ek9,ek9t,ch2o,eh2o,dho2mc,dho2kg,eeee,xk9,
     * r5,d5,dmssink,bdy
#ifdef TRACERS_HETCHEM
     *       ,d41,d42,d43,o3mc,rsulfo3
#endif
      real*8 bciage,ociage
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: ohsr_in
      integer i,j,l,n,iuc,iun,itau,ichemi,itt,
     * ittime,isp,iix,jjx,llx,ii,jj,ll,iuc2,it,najl,j_0,j_1,
     * j_0s,j_1s,mmm,J_0H,J_1H,I_0,I_1
#ifdef TRACERS_SPECIAL_Shindell
!@var maxl chosen tropopause 0=LTROPO(I,J), 1=LS1-1
#endif
#ifdef TRACERS_AEROSOLS_VBS
      type(vbs_tracers) :: vbs_tr_old ! concentrations, ug m-3
      type(vbs_conditions) :: vbs_cond ! current box conditions (meteo+chem)
!@var kg2ugm3 factor to convert kilograms gridbox-1 to ug m-3
      real*8 :: kg2ugm3
#endif /* TRACERS_AEROSOLS_VBS */
      integer maxl,nrecs_skip
      logical :: newMonth
      save ifirst
#ifdef TRACERS_TOMAS
      REAL*8 TAU_hydro
      integer k
#endif

      call getDomainBounds(grid, J_STRT=J_0,J_STOP=J_1,
     *    J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,J_STRT_SKP=J_0S,
     * J_STOP_SKP=J_1S)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** initialise source arrays
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_DMS)=0. ! DMS chem sink
#ifndef TRACERS_AMP
#ifndef TRACERS_TOMAS
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_MSA)=0. ! MSA chem source
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_SO4)=0. ! SO4 chem source
#endif
#endif
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_SO2)=0. ! SO2 chem source
        tr3Dsource(:,j_0:j_1,:,nChemloss,n_SO2)=0. ! SO2 chem sink
        if(n_H2O2_s>0) tr3Dsource(:,j_0:j_1,:,nChemistry,n_H2O2_s)=0. ! H2O2 chem source
        if(n_H2O2_s>0) tr3Dsource(:,j_0:j_1,:,nChemLoss,n_H2O2_s)=0. ! H2O2 chem sink
#ifdef TRACERS_AMP
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_H2SO4)=0. ! H2SO4 chem source
#endif
#ifdef TRACERS_TOMAS
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_H2SO4)=0. ! H2SO4 chem source
        H2SO4_chem(:,j_0:j_1,:)=0.0
        do k=1,nbins
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AECOB(k))=0.
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AECIL(k))=0.
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AOCOB(k))=0.
           tr3Dsource(:,j_0:j_1,:,nChemistry,n_AOCIL(k))=0.
        enddo
#endif
#ifdef TRACERS_HETCHEM
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_SO4_d1) =0. ! SO4 on dust
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_SO4_d2) =0. ! SO4 on dust
        tr3Dsource(:,j_0:j_1,:,nChemistry,n_SO4_d3) =0. ! SO4 on dust
#endif
        if (n_BCII.gt.0) then
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_BCII)=0. ! BCII sink
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_BCIA)=0. ! BCIA source
        end if
        if (n_OCII.gt.0) then
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_OCII)=0. ! OCII sink
          tr3Dsource(:,j_0:j_1,:,nChemistry,n_OCIA)=0. ! OCIA source
        end if
#ifdef TRACERS_AEROSOLS_VBS
        tr3Dsource(:,j_0:j_1,:,nChemistry,vbs_tr%igas)=0.
        tr3Dsource(:,j_0:j_1,:,nChemloss,vbs_tr%igas)=0.
        tr3Dsource(:,j_0:j_1,:,nOther,vbs_tr%igas)=0.
        tr3Dsource(:,j_0:j_1,:,nChemistry,vbs_tr%iaer)=0.
#endif

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
C Coupled mode: use on-line radical concentrations
      if (coupled_chem.eq.1) then
        oh(:,j_0:j_1,:)=oh_live(:,j_0:j_1,:)
        tno3(:,j_0:j_1,:)=no3_live(:,j_0:j_1,:)
c Set h2o2_s =0 and use on-line h2o2 from chemistry
        if(n_H2O2_s>0) trm(:,j_0:j_1,:,n_h2o2_s)=0.0
      endif

      if (coupled_chem.eq.0) then
c Use this for chem inputs from B4360C0M23, from Drew
c      if (ifirst) then
        newMonth = jMonthCache /= modelEclock%getMonth()
        if (newMonth) then
          jMonthCache = modelEclock%getMonth()
          call openunit('AER_CHEM',iuc,.true.)
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),ohrCache,
     &         recs_to_skip=5*(modelEclock%getMonth()-1)+1)    ! 5 recs/month + ichemi for this month
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),dho2rCache)
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),perjrCache)
          call DREAD8_PARALLEL(grid,iuc,nameunit(iuc),tno3rCache)
          call closeunit(iuc)
        end if
        ohr   = ohrCache  
        dho2r = dho2rCache
        perjr = perjrCache
        tno3r = tno3rCache
        if (im.eq.72) then
        call openunit('AER_OH_STRAT',iuc2,.true.)
        nrecs_skip=lm*(modelEclock%getMonth()-1) ! skip all the preceding months
        do ll=1,lm
          call DREAD_PARALLEL(grid,iuc2,nameunit(iuc2),ohsr_in,
     &       recs_to_skip=nrecs_skip)
          ohsr(:,:,ll)=ohsr_in(:,:)*1.D5
          nrecs_skip=0 ! do not skip any more records
        enddo
        call closeunit(iuc2)

c skip poles because there was a bug in the input file over the pole
        do j=j_0s,j_1s   
        do i=i_0,i_1
          maxl=ltropo(i,j)
        do l=maxl,lm
          ohr(i,j,l)=ohsr(i,j,l)
        end do
        end do
        end do
cdmk turning these off because I do not have 2x2.5
c       ifirst=.false.
        else    !need to scale inputs (10^5 mol/cm3)
        ohr(:,:,:)=ohr(:,:,:)*1.D5
        tno3r(:,:,:)=tno3r(:,:,:)*1.D5
        dho2r(:,:,:)=dho2r(:,:,:)*1.D7
        perjr(:,:,:)=perjr(:,:,:)*1.D2
        endif   !im.eq.72
c I have to read in every timestep unless I can find a better way
c
c impose diurnal variability
        CALL SCALERAD
c       write(6,*) ' RRR OXID2 ',ohr(10,45,1),
c    *   oh(10,45,1),dho2r(3,45,1),dho2(3,45,1)
       endif   !coupled_chem.eq.0

C Calculation of gas phase reaction rates
C Now called from tracer_3Dsource
c#ifndef  TRACERS_SPECIAL_Shindell  
c        CALL GET_SULF_GAS_RATES
c#endif

#ifdef TRACERS_HETCHEM
c calculation of heterogeneous reaction rates: SO2 on dust 
      CALL SULFDUST
c calculation of heterogeneous reaction rates: SO2 on seasalt
c      CALL SULFSEAS 
c     if (COUPLED_CHEM.ne.1) then
c     CALL GET_O3_OFFLINE
c     endif
#endif
#endif
      dtt=dtsrc
      !efold time of 1 days
      bciage=(1.d0-exp(-dtsrc/(1.0d0*SECONDS_PER_DAY)))/dtsrc 
      !efold time of 1.6 days
      ociage=(1.d0-exp(-dtsrc/(1.6d0*SECONDS_PER_DAY)))/dtsrc
C**** THIS LOOP SHOULD BE PARALLELISED
      do 20 l=1,lm
      do 21 j=j_0,j_1
      do 22 i=i_0,imaxj(j)
C Initialise       
        bdy = dclev(i,j) 
      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm = MA(l,i,j)*axyp(i,j)
      tt = 1.d0/te

c DMM is number density of air in molecules/cm3
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !oh is alread in units of molecules/cm3

      do 23 n=1,NTM ! ===== THIS IS CHEMISTRY OF Koch AEROSOLS =====

        select case (trname(n))
c    Aging of industrial carbonaceous aerosols 
        case ('BCII')
          tr3Dsource(i,j,l,nChemistry,n)=-bciage*trm(i,j,l,n)
          tr3Dsource(i,j,l,nChemistry,n_BCIA)=bciage*trm(i,j,l,n)

#ifdef TRACERS_AEROSOLS_VBS
        case ('vbsAm2') ! This handles all VBS tracers
          kg2ugm3=1.d9*(1.d2*pmid(l,i,j))*mair/
     &            (MA(l,i,j)*axyp(i,j)*gasc*te)
          vbs_cond%dt=dtsrc
          vbs_cond%OH=ohmc
          vbs_cond%temp=te
          vbs_cond%nvoa=(trm(i,j,l,n_BCII)
     &                  +trm(i,j,l,n_BCIA)
     &                  +trm(i,j,l,n_BCB)
#ifdef TRACERS_AEROSOLS_SOA
     &                  +trm(i,j,l,n_isopp1a)
     &                  +trm(i,j,l,n_isopp2a)
     &                  +trm(i,j,l,n_apinp1a)
     &                  +trm(i,j,l,n_apinp2a)
#endif /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
     &                  +trm(i,j,l,n_ococean)
#endif  /* TRACERS_AEROSOLS_OCEAN */
     &                  +trm(i,j,l,n_msa)
     &                  +trm(i,j,l,n_so4)
#ifdef TRACERS_NITRATE
     &                  +trm(i,j,l,n_nh4)
     &                  +trm(i,j,l,n_no3p)
#endif
     &                  )*kg2ugm3
          vbs_tr_old%gas=trm(i,j,l,vbs_tr%igas)*kg2ugm3
          vbs_tr_old%aer=trm(i,j,l,vbs_tr%iaer)*kg2ugm3

          call vbs_calc(vbs_tr_old,vbs_cond)

          tr3Dsource(i,j,l,nChemistry,vbs_tr%igas)=
     &      vbs_tr%chem_prod/kg2ugm3/vbs_cond%dt
          tr3Dsource(i,j,l,nChemloss,vbs_tr%igas)=
     &      vbs_tr%chem_loss/kg2ugm3/vbs_cond%dt
          tr3Dsource(i,j,l,nOther,vbs_tr%igas)=
     &      -vbs_tr%partition/kg2ugm3/vbs_cond%dt ! partitioning
          tr3Dsource(i,j,l,nChemistry,vbs_tr%iaer)=
     &      vbs_tr%partition/kg2ugm3/vbs_cond%dt
!     &      (vbs_tr%gas-vbs_tr_old%gas)/kg2ugm3/vbs_cond%dt
!      if (sum(vbs_tr_old%gas)+sum(vbs_tr_old%aer) /= 0.) then
!        print '(a,3e)','KOSTAS gas',
!     &                 sum(vbs_tr_old%gas),
!     &                 sum(vbs_tr%gas),
!     &                 sum(vbs_tr_old%gas)+sum(vbs_tr_old%aer)
!        print '(a,3e)','KOSTAS aer',
!     &                 sum(vbs_tr_old%aer),
!     &                 sum(vbs_tr%aer),
!     &                 sum(vbs_tr%gas)+sum(vbs_tr%aer)
!        print '(a,3e)','KOSTAS bud',
!     &                 sum(vbs_tr%chem_prod),
!     &                 sum(vbs_tr%chem_loss),
!     &                 sum(vbs_tr%partition)
!      endif
#else
        case ('OCII')
          tr3Dsource(i,j,l,nChemistry,n)=-ociage*trm(i,j,l,n)
          tr3Dsource(i,j,l,nChemistry,n_OCIA)=ociage*trm(i,j,l,n)
#endif /* TRACERS_AEROSOLS_VBS */

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
        case ('DMS')
C***1.DMS + OH -> 0.75SO2 + 0.25MSA
C***2.DMS + OH -> SO2
C***3.DMS + NO3 -> HNO3 + SO2

          r1=rsulf1(i,j,l)*ohmc 
          d1 = exp(-r1*dtsrc)
          r2=rsulf2(i,j,l)*ohmc
          d2 = exp(-r2*dtsrc)

c     NO3 is in mixing ratio: convert to molecules/cm3
c - not necessary for Shindell source
          if (l.gt.bdy) then
            ttno3=0.d0
            go to 87
          endif
          ttno3 = tno3(i,j,l)   !*6.02d20*ppres/(.082056d0*te)
 87       r3=rsulf3(i,j,l)*ttno3
          d3= exp(-r3*dtsrc)
          ddno3=r3*trm(i,j,l,n)/tr_mm(n)*1000.d0*dtsrc
          dddms=trm(i,j,l,n)/tr_mm(n)*1000.d0
          if (ddno3.gt.dddms) ddno3=dddms

          ddno3=ddno3*0.9
C DMS losses: eqns 1, 2 ,3

          tr3Dsource(i,j,l,nChemistry,n) = trm(i,j,l,n)*(d1*d2-1.)/dtsrc

          dmssink=ddno3*tr_mm(n)/1000.d0

          if (dmssink.gt.trm(i,j,l,n)+
     *        tr3Dsource(i,j,l,nChemistry,n)*dtsrc)
     *         dmssink=trm(i,j,l,n)+
     *                 tr3Dsource(i,j,l,nChemistry,n)*dtsrc
          tr3Dsource(i,j,l,nChemistry,n)=tr3Dsource(i,j,l,nChemistry,n)-
     *                                   dmssink/dtsrc
          
        case ('MSA')
C MSA gain: eqn 1

          tr3Dsource(i,j,l,nChemistry,n) = 0.25d0*Tr_mm(n)/Tr_mm(n_dms)
     *         *trm(i,j,l,n_dms)*(1.d0 -D1)*SQRT(D2)/dtsrc
          
        case ('SO2')
c SO2 production from DMS
          tr3Dsource(i,j,l,nChemistry,n) = (0.75*tr_mm(n)/tr_mm(n_dms)
     *         *trm(i,j,l
     *         ,n_dms)*(1.d0 - d1)*sqrt(d2)+ tr_mm(n)/tr_mm(n_dms)*trm(i
     *         ,j,l,n_dms)*(1.d0 - d2)*sqrt(d1)+dmssink*tr_mm(n)
     *         /tr_mm(n_dms))/dtsrc
#ifdef TRACERS_TOMAS 
! EC/OC aging 
        case ('AECIL_01')
           TAU_hydro=1.5D0*SECONDS_PER_DAY !24.D0*3600.D0 !1.5 day 

           DO K=1,nbins
              tr3Dsource(i,j,l,nChemistry,n_AECIL(K))=
     &             trm(i,j,l,n_AECOB(K))*
     &             (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc 

              tr3Dsource(i,j,l,nChemistry,n_AECOB(K))=
     &             -trm(i,j,l,n_AECOB(K))*
     &            (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc 

!              IF(am_i_root())
!         print*,'ECOB aging',k,n_AECOB(K),(1.D0-EXP(-dtsrc/TAU_hydro))
!     &  , trm(i,j,l,n_AECOB(K))
           ENDDO
           
        case ('AOCIL_01')
           TAU_hydro=1.5D0*SECONDS_PER_DAY !24.D0*3600.D0 !1.5 day 

           DO K=1,nbins
              tr3Dsource(i,j,l,nChemistry,n_AOCIL(K))
     &             =trm(i,j,l,n_AOCOB(K))*
     &             (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc
!     &             4.3D-6
              tr3Dsource(i,j,l,nChemistry,n_AOCOB(K))
     &             =-trm(i,j,l,n_AOCOB(K))*
     &            (1.D0-EXP(-dtsrc/TAU_hydro))/dtsrc 
!     &             4.3D-6
!              IF(am_i_root())
!         print*,'OCOB aging',k,n_AOCOB(K),(1.D0-EXP(-dtsrc/TAU_hydro))
!     &  , trm(i,j,l,n_AOCOB(K))

           ENDDO   

#endif
#ifndef TRACERS_TOMAS                   
          najl = jls_NO3
          call inc_tajls2(i,j,l,najl,ttno3)
#endif
#endif
        end select
        
 23   CONTINUE ! ===== END OF CHEMISTRY OF Koch AEROSOLS ====
 22   CONTINUE
 21   CONTINUE
 20   CONTINUE

      do 30 l=1,lm
      do 31 j=j_0,j_1
      do 32 i=i_0,imaxj(j)

      ppres=pmid(l,i,j)*9.869d-4 !in atm
      te=pk(l,i,j)*t(i,j,l)
      mm = MA(l,i,j)*axyp(i,j)
      tt = 1.d0/te
      dmm=ppres/(.082d0*te)*6.02d20
      ohmc = oh(i,j,l)          !OH is already in units of molecules/cm3
#ifdef TRACERS_HETCHEM
      if (COUPLED_CHEM.ne.1) then
      o3mc=o3_offline(i,j,l)*dmm*(28.0D0/48.0D0)*BYAXYP(I,J)*byMA(L,I,J)
      else
        o3mc=trm(i,j,l,n_Ox)*dmm*(28.0D0/48.0D0)*BYAXYP(I,J)*byMA(L,I,J)
      endif
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS)
      do 33 n=1,NTM
        select case (trname(n))
        case ('SO2')
c oxidation of SO2 to make SO4: SO2 + OH -> H2SO4

          r4=rsulf4(i,j,l)*ohmc
          d4 = exp(-r4*dtsrc)

          IF (d4.GE.1.) d4=0.99999d0
#ifdef TRACERS_HETCHEM
      rsulfo3 = 4.39d11*exp(-4131/te)+( 2.56d3*exp(-966/te)) * 10.d5 !assuming pH=5
      rsulfo3 = exp(-rsulfo3*o3mc *dtsrc) !O3 oxidation Maahs '83
       d41 = exp(-rxts1(i,j,l)*dtsrc)     
       d42 = exp(-rxts2(i,j,l)*dtsrc)     
       d43 = exp(-rxts3(i,j,l)*dtsrc)     
       tr3Dsource(i,j,l,nChemloss,n) = (-trm(i,j,l,n)*(1.d0-d41)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d4)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d42)/dtsrc)
     .                       + ( -trm(i,j,l,n)*(1.d0-d43)/dtsrc)
#else
       tr3Dsource(i,j,l,nChemloss,n) = -trm(i,j,l,n)*(1.d0-d4)/dtsrc 
#ifdef TRACERS_AMP
       tr3Dsource(i,j,l,nChemistry,n_H2SO4)=trm(i,j,l,n)*(1.d0-d4)/dtsrc
     &      *tr_mm(n_H2SO4)/tr_mm(n)
#endif  /* TRACERS_AMP */
#ifdef TRACERS_TOMAS
       H2SO4_chem(i,j,l)=trm(i,j,l,n)*(1.d0-d4)/dtsrc 
     &      *tr_mm(n_H2SO4)/tr_mm(n) 
       tr3Dsource(i,j,l,nChemistry,n_H2SO4)=H2SO4_chem(i,j,l)
#endif  /* TRACERS_TOMAS */
#endif  /* TRACERS_HETCHEM */
c diagnostics to save oxidant fields
c No need to accumulate Shindell version here because it
c   is done elsewhere
c#ifdef TRACERS_SPECIAL_Shindell
c         najl = jls_OHcon
c#else
          najl = jls_OHconk
c#endif
          if (najl > 0) call inc_tajls2(i,j,l,najl,oh(i,j,l))
          najl = jls_HO2con
          if (najl > 0) call inc_tajls2(i,j,l,najl,dho2(i,j,l))

#ifdef TRACERS_HETCHEM
       case ('SO4_d1')
c sulfate production from SO2 on mineral dust aerosol due to O3 oxidation

       tr3Dsource(i,j,l,nChemistry,n)=tr3Dsource(i,j,l,nChemistry,n)
     *         +tr_mm(n)/tr_mm(n_so2)
     *         *(1.d0-d41)*trm(i,j,l,n_so2)            !  SO2
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc
       case ('SO4_d2')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,nChemistry,n) = tr3Dsource(i,j,l,nChemistry,n)
     *         +tr_mm(n)/tr_mm(n_so2)*(1.d0-d42)*trm(i,j,l,n_so2)
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc
       case ('SO4_d3')
c sulfate production from SO2 on mineral dust aerosol

       tr3Dsource(i,j,l,nChemistry,n) = tr3Dsource(i,j,l,nChemistry,n)
     *         +tr_mm(n)/tr_mm(n_so2)*(1.d0-d43)*trm(i,j,l,n_so2)
     *         * (1.d0-rsulfo3)                              !+ O3
     *           /dtsrc

#endif
        case('SO4')
C SO4 production
          tr3Dsource(i,j,l,nChemistry,n) = 
     *         tr3Dsource(i,j,l,nChemistry,n)+tr_mm(n)
     *         /tr_mm(n_so2)*trm(i,j,l,n_so2)*(1.d0 -d4)/dtsrc
        case('H2O2_s')

          if (coupled_chem.ne.1) then

c hydrogen peroxide formation and destruction:
C***5.H2O2 +hv -> 2OH
C***6.H2O2 + OH -> H2O + HO2
C***7.H2O2 + SO2 -> H2O + SO3 (in-cloud, in CB)
C***9.HO2 + HO2 ->H2O2 + O2
C     HO2 + HO2 + M ->
C     HO2 + HO2 + H2O ->
C     HO2 + HO2 + H2O + M ->

          r6 = 2.9d-12 * exp(-160.d0*tt)*ohmc
          d6 = exp(-r6*dtsrc)
          ek9 = 2.2d-13*exp(600.d0*tt)
          ek9t = 1.9d-20*dmm*0.78d0*exp(980.d0*tt)*1.d-13
          ch2o = q(i,j,l)*6.02d20*28.97d0/18.d0*ppres/(.082d0*te)
          eh2o = 1.+1.4d-21*exp(2200.d0*tt)*ch2o
          dho2mc = dho2(i,j,l)  !/mm*1.292/.033*6.02e17
          dho2kg = dho2(i,j,l)*mm*te*.082056d0/(ppres*28.97d0*6.02d20)
          eeee = eh2o*(ek9+ek9t)*dtt*dho2mc
          xk9 = dho2kg*eeee
c         if (i.eq.2.and.l.eq.1.and.j.eq.46) write(6,*) 
c    *    'RRR CHEM DEBUG ',i,j,xk9,dho2kg,eeee,dho2mc
c    *    ,eh2o,ek9,ek9t,dtt
c         if (i.eq.72.and.l.eq.1.and.j.le.46) write(6,*) 
c    *    'RRR CHEM DEBUG ',i,j,xk9,dho2kg,eeee,dho2mc
c H2O2 production: eqn 9
         
          tr3Dsource(i,j,l,nChemistry,n) = tr_mm(n)*xk9/dtsrc
c        if (i.eq.10.and.j.eq.45.and.l.eq.1) then
c        write(6,*) 'RRR OXID H2O2',xk9,dho2kg,eeee
c         endif
c H2O2 losses:5 and 6
          r5 = perj(i,j,l)
          d5 = exp(-r5*dtsrc)

          tr3Dsource(i,j,l,nChemLoss,n)=(trm(i,j,l,n))*(d5*d6-1.d0)
     *         /dtsrc
          
          najl = jls_phot
          if (najl > 0) call inc_tajls(i,j,l,najl,perj(i,j,l))
          endif
        end select

 33   CONTINUE
#endif
 32   CONTINUE
 31   CONTINUE
 30   CONTINUE

      RETURN
      END SUBROUTINE aerosol_gas_chem

      SUBROUTINE SCALERAD
      use constant, only : pi
      use resolution, only: lm
      use AEROSOL_SOURCES, only: ohr,dho2r,perjr,tno3r,oh,dho2,perj,tno3
      USE DOMAIN_DECOMP_ATM, only:GRID, getDomainBounds
      use RAD_COM, only: cosz1,cosz_day,sunset
      implicit none
      real*8, parameter ::
     &     night_frac_min=.01d0 ! minimum night_frac for tno3 scaling
      real*8 stfac,night_frac
      integer i,j,l,i_0,i_1,j_0,j_1
c      real*8, DIMENSION(JM) :: tczen
c      integer, DIMENSION(JM) :: nradn

      call getDomainBounds(grid, 
     &     I_STRT=I_0,I_STOP=I_1, J_STRT=J_0,J_STOP=J_1)

c      nradn(:)=0
c      tczen(:)=0.d0
c      do 100 j = j_0,j_1
c      do 100 i = 1, im
c      if (cosz1(i,j).gt.0.) then
c      tczen(j)=tczen(j)+cosz1(i,j)
c      else
c      nradn(j)=nradn(j)+1
c      endif
c 100  continue

      do j = j_0,j_1
      do i = i_0,i_1
c Get NO3 only if dark, weighted by number of dark hours
        night_frac = 1.-sunset(i,j)/pi
c        night_frac = real(nradn(j))/real(im)
        if (cosz1(i,j).le.0.and.night_frac.gt.night_frac_min) then
          tno3(i,j,:)=tno3r(i,j,:)/night_frac !DMK jmon
        else
          tno3(i,j,:)=0.d0
        endif
        if (cosz1(i,j).gt.0.) then
c          stfac=cosz1(i,j)/tczen(j)*real(IM)
          stfac=cosz1(i,j)/cosz_day(i,j)
          oh(i,j,:)=ohr(i,j,:)*stfac
          perj(i,j,:)=perjr(i,j,:)*stfac
          dho2(i,j,:)=dho2r(i,j,:)*stfac
        end if
      end do
      end do

c        if (I.EQ.1.AND.L.EQ.1) write(6,*)'NO3R',TAU,J,NRADN(I,J)
c            if (l.eq.1.and.j.eq.23.and.i.eq.10) write(6,*)
c    *     'RRR SCALE ',stfac,cosz1(i,j),tczen(j),oh(i,j,l),ohr(i,j,l)
      RETURN
      END SUBROUTINE SCALERAD



      SUBROUTINE GET_SULFATE(pl,temp_in,fcloud,
     *  wa_vol,wmxtr,sulfin,sulfinom,sulfinc,sulfout,tr_left,
     *  tmg,tmd,airm,lhx,dt_sulf,fcld0)

!@sum  GET_SULFATE calculates formation of sulfate from SO2 and H2O2
!@+    within or below convective or large-scale clouds. Gas
!@+    condensation uses Henry's Law if not freezing.
!@auth Dorothy Koch

!**** GLOBAL parameters and variables:
      USE CONSTANT, only: bygasc, MAIR,teeny,mb2kg,gasc,lhe
      use OldTracer_mod, only: trname, mass2vol, tr_mm
      use OldTracer_mod, only: tr_RKD, tr_DHD
      USE TRACER_COM, only: n_H2O2_s,n_SO2
     *     ,NTM
     *     ,lm,n_SO4,n_H2O2,coupled_chem
      use tracer_com, only: aqchem_count,aqchem_list
      USE CLOUDS, only: NTX,DXYPIJ
      USE MODEL_COM, only: dtsrc

      IMPLICIT NONE

!**** Local parameters and variables and arguments:

!@var sulfin amount of precursor used to make product from the gas phase (kg)
!@var sulfinom sulfin=sulfinom*tmg, for avoiding divisions sulfinom=sulfin/tmg
!@+   (dimensionless)
!@var sulfinc amount of precursor used to make product from the condensate (kg)
!@var sulfout total amount of product generated (kg)
!@var tr_left is the amount of precursor left after product is made
!@+   and is now available to condense
!@+   This is a very strange variable, probably wrong!!!
      real*8, dimension(aqchem_count), intent(out) :: sulfin,sulfinom,
     &                                        sulfinc,sulfout,tr_left
!@var fcloud cloud fraction available for tracer condensation. fcloud=fplume
!@+   for convective clouds, and fcloud=fcld for large-scale clouds
!@var fcld0 updated cloud fraction, given the current state of large-scale
!@+   clouds. fcld0=0.d0 for convective clouds.
!@var lhx latent heat of evaporation or sublimation (J/Kg). When equal to lhe
!@+   the cloud is in the ice phase.
!@var finc XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      real*8, intent(in) :: fcloud,fcld0,lhx
      real*8 :: finc
!@var airm layer pressure depth (mb). Multiply by mb2kg to convert to air mass
!@+   per m2, based on the hydrostatic pressure equation:
!@+   pressure (Pa=kg/m/s2) = height (m) * density (kg/m3) * g (m/s2)
!@var amass airmass in kg, calculated by airm*mb2kg*dxypij
      real*8, intent(in) :: airm
      real*8 :: amass

!@var pl pressure at current altitude (mbar)
!@var press pressure at current altitude (Pa)
!@var temp temperature to be used (K), always greater than 230K
!@var tfac exponent factor for temperature dependence calculations (mol/J).
!@+   tfac = (1/temp - 1/298)/R; R=8.31451 J/mol/K
!@var clwc cloud liquid water content (volume water/volume air)
!@var temp_in temperature at current altitude (K)
!@var wmxtr cloud water mixing ratio (kg water/kg air)
      real*8 :: press, temp, tfac, clwc
      real*8, intent(in) :: pl,temp_in,wmxtr

!@param k1so2dissoc0 first dissociation rate of dissolved SO2 at 298K.
!@+     SO2.H2O <--> H+ + HSO3-
!@param dh1so2dissoc enthalpy of dissociation for the first dissociation rate
!@+     of dissolved SO2
!@var k1so2dissoc first dissociation rate of dissolved SO2 at current temp.
      real*8, parameter :: k1so2dissoc0=1.3d-2 ! M (molar)
      real*8, parameter :: dh1so2dissoc=-1.6736d4 !J/mol
      real*8 :: k1so2dissoc ! M (molar)
!@param Hplus concentration of H+ (molar) for pH=4.5
!@var henry modified henry constant of the current tracer at current
!@+   conditions, taking into account the current pH, if needed (moles/J).
!@+   multiply with convert_HSTAR to convert to moles/liter/atm
      real*8, parameter :: Hplus=10.d0**(-4.5d0)
      real*8 :: henry
!@param kso2h2o20 reaction rate of SO2 + H2O2 at 298K. SO2 + H2O2 --> SO3 + H2O
!@param dhso2h2o2 enthalpy of reaction of SO2 + H2O2
!@var kso2h2o2 reaction rate of SO2 + H2O2 at current temp.
      real*8, parameter :: kso2h2o20=6.357d14    !1/(M*M*s)
      real*8, parameter :: dhso2h2o2=3.95d4 !J/mol
      real*8 :: kso2h2o2

!@var ix index of current species
!@var is index of current species in ntx array
!@var ih index of current species in ntx array
!@var isx index of SO2 species in aqchem_list array
!@var ihx index of H2O2 species in aqchem_list array
      integer :: ix,is,ih,isx,ihx

!@var tmg amount of tracer in the gas phase in the cloudy area (kg).
!@var tmd amount of tracer in the aqueous phase (kg).
!@var tmgmol amount of gas phase tracer in cloudy area (moles)
!@var tmdmol amount of gas phase tracer in cloudy area (moles)
!@var tmgrate is the new concentration of species that resulted from the
!@+   dissolution of its gas phase precursor species, on top of what was there
!@+   from the previous timestep. The units are M/kg, meaning molarity
!@+   produced per kilogram reacted.
!@var tmdrate is the new concentration of species that resulted from the
!@+   already dissolved precursor species, on top of what was there
!@+   from the previous timestep. The units are M/kg, meaning molarity
!@+   produced per kilogram reacted.
      real*8, dimension(ntx), intent(in) :: tmg,tmd
      real*8, dimension(aqchem_count) :: tmgmol,tmdmol,tmgrate,tmdrate

!@var wa_vol cloud water volume (liters)
!@var dso4g amount of sulfate produced from the gas phase (moles/kg/kg)
!@var dso4d amount of sulfate produced from the condensate phase (moles/kg/kg)
!@var dso4gt amount of sulfate produced from the gas phase (moles)
!@var dso4dt amount of sulfate produced from the condensate phase (moles)
      real*8,  intent(in) :: wa_vol
      real*8 :: dso4g,dso4d,dso4gt,dso4dt

!@var n index for tracer number loop
      integer :: n

!@var dt_sulf accumulated diagnostic of sulfate chemistry changes
      real*8, dimension(ntx), intent(inout) :: dt_sulf

      sulfin(:)=0.d0
      sulfinom(:)=0.d0
      sulfinc(:)=0.d0
      sulfout(:)=0.d0
      tr_left(:)=1.d0

! if no water clouds or no clouds at all, do nothing
      if (lhx.ne.lhe.or.fcloud.lt.teeny.or.wmxtr.le.teeny) return

! calculate the fraction of tracer mass that becomes condensate
      finc=(fcloud-fcld0)/fcloud
      if (finc.lt.0.d0) finc=0.d0

! calculate some variables for later
      amass=airm*mb2kg*dxypij ! kg
      press = pl*1.d2 ! Pa
! comment from Dorothy Koch:
! calls to this subroutine are sometimes made at stages of the cloud scheme
! at which some but not all tendencies have been applied to temp_in, so we
! impose a lower limit (liquid water is very unlikely to exist below 230 K)
      temp = max(temp_in, 230.d0) ! K
      tfac = (1.d0/temp - 1.d0/298.d0)*bygasc  ! mol/J
! cloud liquid water content
      clwc=wmxtr*mair*press/temp*bygasc/1.d6/fcloud ! volume water/volume air

      k1so2dissoc=k1so2dissoc0*exp(-dh1so2dissoc*tfac) ! SO2.H2O <--> H+ + HSO3-
      kso2h2o2=kso2h2o20*exp(-dhso2h2o2/(gasc*temp)) ! SO2 + H2O2 --> SO3 + H2O

! First allow for formation of sulfate from SO2 and H2O2. Then remaining
! gases may be allowed to dissolve (amount given by tr_left)
! H2O2 + SO2 -> H2O + SO3 -> H2SO4
      do n=1,aqchem_count
        ix=aqchem_list(n)

        select case (trname(ix))
        case('SO2', 'H2O2', 'H2O2_s')

! save some per-tracer values needed later
          select case (trname(ix))
          case('SO2')
            is=ix
            isx=n
          case('H2O2','H2O2_s')
            ih=ix
            ihx=n
            select case (trname(ix))
            case('H2O2')
              if (coupled_chem.eq.0) goto 400
            case('H2O2_s')
              if (coupled_chem.eq.1) goto 400
            end select
          end select

! initial amount of species in the gas and aqueous phases
          tmgmol(n)=1.d3*tmg(ix)/tr_mm(ix) ! gas-phase, in moles
          tmdmol(n)=tmd(ix)*1.d3/tr_mm(ix) ! aqueous phase, in moles

! henry coefficient
          henry=tr_rkd(ix)*exp(-tr_dhd(ix)*tfac) ! moles/J

! partial pressure of gas x henry's law coefficient
          tmgrate(n)=mass2vol(ix)*1.d-3*press/amass*henry

! modified Henry's Law coefficient assuming pH of 4.5
          select case (trname(ix))
          case('SO2')
            henry=henry*(1.d0+ k1so2dissoc/Hplus)
          end select

! rate of production from gaseous and dissolved precursors.
          tmgrate(n)=tmgrate(n)/(1.d0+henry*clwc*gasc*temp) ! dimless henry
          tmdrate(n)=mass2vol(ix)*press/amass*bygasc/temp*1.d-3/clwc

 400      continue
        end select
      enddo

! do not calculate dso4g if there is not enough gas phase to react
      if (tmg(ih).lt.teeny.or.tmg(is).lt.teeny) then
        dso4g=0.d0
        dso4gt=0.d0
        go to 21
      endif

! production from the gas phase, moles/kg/kg
      dso4g=kso2h2o2*k1so2dissoc*tmgrate(ihx)*tmgrate(isx)*dtsrc*wa_vol
      dso4g=dso4g*finc ! increase production based on current cloud water volume
      dso4gt=dso4g*tmg(ih)*tmg(is) ! moles

! can't be more than the moles we started with
      if (dso4gt.gt.tmgmol(isx)) then ! so2
        dso4g=tmgmol(isx)/(tmg(ih)*tmg(is))
        dso4gt=tmgmol(isx)
      endif
      if (dso4gt.gt.tmgmol(ihx)) then ! h2o2
        dso4g=tmgmol(ihx)/(tmg(ih)*tmg(is))
        dso4gt=tmgmol(ihx)
      endif
 21   continue

! do not calculate dso4d if there is not enough dissolved phase to react
      if (tmd(ih).lt.teeny.or.tmd(is).lt.teeny) then
        dso4d=0.d0
        dso4dt=0.d0
        go to 22
      endif

! production from the already-dissolved aqueous phase, moles/kg/kg
      dso4d=kso2h2o2*k1so2dissoc*tmdrate(ihx)*tmdrate(isx)*dtsrc*wa_vol
      dso4dt=dso4d*tmd(ih)*tmd(is) ! moles

! can't be more than the moles we started with
      if (dso4dt.gt.tmdmol(isx)) then
        dso4d=tmdmol(isx)/(tmd(ih)*tmd(is))
        dso4dt=dso4d*tmd(ih)*tmd(is)
      endif
      if (dso4dt.gt.tmdmol(ihx)) then
        dso4d=tmdmol(ihx)/(tmd(ih)*tmd(is))
        dso4dt=dso4d*tmd(ih)*tmd(is)
      endif
 22   continue

! save final concentrations and diagnostics
      do n=1,aqchem_count
        ix=aqchem_list(n)
        select case (trname(ix))
        case('SO4','M_ACC_SU','ASO4__01')
          sulfout(n)=tr_mm(ix)/1.d3*(dso4gt+dso4dt) ! kg
          dt_sulf(ix) = dt_sulf(ix) + sulfout(n)

        case('SO2','H2O2','H2O2_s')
          select case (trname(ix))
          case('SO2')
            sulfinom(n)=-dso4g*tmg(ih)*tr_mm(ix)/1.d3 ! dimensionless
            sulfinc(n)=-dso4d*tmd(ih)*tr_mm(ix)/1.d3*tmd(ix) ! kg
          case('H2O2','H2O2_s')
            select case (trname(ix))
            case('H2O2')
              if (coupled_chem.eq.0) goto 401
            case('H2O2_s')
              if (coupled_chem.eq.1) goto 401
            end select
            sulfinom(n)=-dso4g*tmg(is)*tr_mm(ix)/1.d3 ! dimensionless
            sulfinc(n)=-dso4d*tmd(is)*tr_mm(ix)/1.d3*tmd(ix) ! kg
          end select
          sulfin(n)=sulfinom(n)*tmg(ix) ! kg

          sulfinom(n)=max(-1.d0,sulfinom(n))
          sulfin(n)=max(-tmg(ix),sulfin(n))
          sulfinc(n)=max(-tmd(ix),sulfinc(n))
          tr_left(n)=0.d0
          if (fcloud.gt.abs(sulfinom(n))) then
            tr_left(n)=fcloud+sulfinom(n)
          endif
 401      continue
          dt_sulf(ix)=dt_sulf(ix)+sulfin(n)+sulfinc(n)

        end select
      enddo

      END SUBROUTINE GET_SULFATE

#ifdef BC_ALB
      SUBROUTINE GET_BC_DALBEDO(i,j,bc_dalb,snow_present)
!@sum Calculates change to albedo of snow on ice and snow on land due
!@+     to BC within the snow.
!@+     Parameterization based on Warren and Wiscombe (1980) (21 inputs)
!@+     or actually on Flanner et al. fig 2 r_e=500 (14 inputs)
!@+     or Warren and Wiscombe (1985) (18 input, old vs new, then I
!@+     continue linearly from 19-29 off the plot)
!@auth Dorothy Koch, modified by Kostas Tsigaridis

c gtracer(n,i,j) is tracer concentration in snow on sea ice?

!@param rhow density of pure water [kg m-3]
      USE CONSTANT, only: rhow
!@var tr_wsn_ij tracer amount in snow over land (multiplied by fr_snow) [kg m-2]
c wsn_ij(nsl,2,i,j)
      USE GHY_COM, only: tr_wsn_ij, wsn_ij
!@var si_atm%snowi snow amount on sea ice [kg m-2]
      USE SEAICE_COM, only : si_atm
#ifdef TRACERS_AEROSOLS_Koch
      use TRACER_COM, only: n_BCB,n_BCII,n_BCIA
#endif
#ifdef TRACERS_AMP
      use TRACER_COM, only:n_M_BC1_BC,n_M_BC2_BC,n_M_BC3_BC,n_M_DBC_BC
     *  ,n_M_BOC_BC,n_M_BCS_BC,n_M_MXX_BC
#endif
#ifdef TRACERS_TOMAS
      use TRACER_COM, only: n_AECIL,n_AECOB,nbins
#endif
      USE FLUXES, only: atmice
      IMPLICIT NONE
c Warren and Wiscombe 1985 includes age dependence
      real*8, parameter :: bc(29)=(/1.d0,2.d0,3.d0,4.d0,5.d0,
     * 6.d0,7.d0,8.d0,9.d0,10.d0,20.d0,30.d0,40.d0,50.d0,60.d0,
     * 70.d0,80.d0,90.d0,100.d0,110.d0,120.d0,130.d0,140.d0,
     * 150.d0,160.d0,170.d0,180.d0,190.d0,200.d0/)
      real*8, parameter :: daln(29)=(/0.d0,0.1d0,0.1d0,0.2d0,
     * 0.2d0,0.2d0,0.2d0,0.3d0,0.3d0,0.4d0,0.7d0,0.9d0,1.1d0,
     * 1.3d0,1.5d0,1.6d0,1.8d0,2.d0,2.2d0,2.4d0,2.6d0,2.8d0,
     * 3.d0,3.2d0,3.4d0,3.6d0,3.8d0,4.d0,4.2d0/)
      real*8, parameter :: dalo(29)=(/0.1d0,0.2d0,0.4d0,0.5d0,
     * 0.6d0,0.7d0,0.8d0,0.9d0,1.d0,1.d0,2.d0,2.6d0,3.2d0,
     * 3.8d0,4.3d0,4.8d0,5.2d0,5.5d0,5.9d0,6.3d0,6.7d0,7.1d0,
     * 7.5d0,7.9d0,8.3d0,8.7d0,9.1d0,9.5d0,9.9d0/)
c Flanner et al
c     real*8, parameter :: bc(14)=(/25.d0,50.d0,100.d0,150.d0,
c    * 200.d0,
c    *250.d0,300.d0,400.d0,500.d0,600.d0,700.d0,800.d0,900.d0,
c    * 1000.d0/)
c     real*8, parameter :: dal(14)=(/1.d0,2.d0,3.d0,4.d0,5.d0,
c    * 6.d0,7.d0,8.d0,9.d0,10.d0,11.d0,11.5d0,12.d0,12.5d0/)
c Warren and Wiscomb 1980
c     real*8, parameter :: bc(21)=(/0.05d0,0.075d0,0.1d0,0.2d0,
c    * 0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.d0,2.d0,
c    * 3.d0,4.d0,5.d0,6.d0,7.d0,8.d0,9.d0,10.d0/)
c     real*8, parameter :: dal(21)=(/2.d0,3.d0,4.d0,6.d0,8.d0,
c    * 9.d0,10.d0,11.d0,12.d0,13.d0,14.d0,16.d0,18.d0,20.d0,
c    * 22.d0,24.d0,26.d0,28.d0,30.d0,32.d0,34.d0/)

!@var bcsnowb BC amount in snow over bare soil [kg m-2]
!@var bcsnowv BC amount in snow over vegetation [kg m-2]
!@var sconb BC concentration in snow over bare soil [kg kg-1]
!@var sconv BC concentration in snow over bare soil [kg kg-1]
!@var scon BC concentration in snow over land [ppm by mass]
!@var icon BC concentration in snow over sea ice [ppm by mass]
      real*8 ::  bcsnowb,bcsnowv,sconb,sconv,scon,icon
!@var fb fraction of land with bare soil (1.-fv)
!@var fv fraction of land with vegetation (1.-fb)
!@var bcc BC concentration (=max(scon,icon)) to be used for albedo calculations
!@var rads snow grain size determined in GRAINS
      real*8 :: fv,fb,bcc,rads
      INTEGER n,ib
      INTEGER, INTENT(IN) :: i,j
      REAL*8, INTENT(OUT) :: bc_dalb
      logical, intent(out) :: snow_present
#ifdef TRACERS_AEROSOLS_Koch
      integer, parameter :: nspBC=3
#endif
#ifdef TRACERS_AMP
      integer, parameter :: nspBC=7
#endif
#ifdef TRACERS_TOMAS
      integer, parameter :: nspBC=nbins+nbins
#endif
      integer, dimension(nspBC) :: spBC

! define indices of BC tracers
#ifdef TRACERS_AEROSOLS_Koch
      spBC(1)=n_BCII
      spBC(2)=n_BCIA
      spBC(3)=n_BCB
#endif
#ifdef TRACERS_AMP
      spBC(1)=n_M_BC1_BC
      spBC(2)=n_M_BC2_BC
      spBC(3)=n_M_BC3_BC
      spBC(4)=n_M_DBC_BC
      spBC(5)=n_M_BOC_BC
      spBC(6)=n_M_BCS_BC
      spBC(7)=n_M_MXX_BC
#endif
#ifdef TRACERS_TOMAS
      do n=1,nbins
         spBC(n)=n_AECOB(n)
         spBC(n+nbins)=n_AECIL(n)
      enddo
#endif

! initialize
      bcsnowb=0.d0
      bcsnowv=0.d0
      sconb=0.d0
      sconv=0.d0
      scon=0.d0
      icon=0.d0
      bc_dalb=0.d0
      snow_present = .false.

! get bare soil and vegetation fractions (fb+fv=1.)
      call get_fb_fv( fb, fv, i, j )

! calculate BC concentration in snow layer 1 over bare soil
      if (wsn_ij(1,1,i,j).gt.0.d0) then
        snow_present = .true. ! should this ignore trace amounts of snow?
        do n=1,nspBC
          bcsnowb=bcsnowb+tr_wsn_ij(spBC(n),1,1,i,j)
        enddo
        sconb=bcsnowb/wsn_ij(1,1,i,j)/rhow
      endif

! calculate BC concentration in snow layer 1 over vegetation
      if (wsn_ij(1,2,i,j).gt.0.d0) then
        snow_present = .true. ! should this ignore trace amounts of snow?
        do n=1,nspBC
          bcsnowv=bcsnowv+tr_wsn_ij(spBC(n),1,2,i,j)
        enddo
        sconv=bcsnowv/wsn_ij(1,2,i,j)/rhow
      endif

! calculate mean BC concentration in total snow in layer 1
      scon=(fb*sconb+fv*sconv)*1.d9

! calculate BC concentration in snow over sea ice
      if (si_atm%snowi(i,j).gt.0.d0) then
        snow_present = .true. ! should this ignore trace amounts of snow?
        do n=1,nspBC
          icon=icon+atmice%gtracer(spBC(n),i,j)*1.d9
        enddo
      endif

! use the maximum BC concentration between snow over land and over sea ice
      bcc=max(icon,scon)

! calculate snow grain size
      call GRAINS(i,j,rads)

! calculate BC albedo effect
      do ib=1,28
        if (bcc.gt.bc(ib).and.bcc.lt.bc(ib+1)) then
          bc_dalb=-(daln(ib)
     *       +(rads-100.d0)/900.d0*(dalo(ib)-daln(ib)))/100.d0
          exit
        endif
      enddo
      if (bcc.ge.bc(29)) bc_dalb=-(daln(29)
     *    + (rads-100.d0)/900.d0*(dalo(29)-daln(29)))/100.d0
c     if (bc_dalb.ne.0.) write(6,*) 'alb_write',i,j,bc_dalb,bcc,rads

      END SUBROUTINE GET_BC_DALBEDO

      SUBROUTINE GRAINS(i,j,rads)
!@sum Estimates snow grain size (microns) based on air temperature
!@+     and snow age. From Susan Marshall's PhD thesis
!@auth Dorothy Koch

      USE CONSTANT, only: pi,gasc,tf
      USE FLUXES, only: atmsrf
      use TimeConstants_mod, only: DAYS_PER_YEAR
      USE RAD_COM, only: snoage
      USE AEROSOL_SOURCES, only: snosiz
      IMPLICIT none
      REAL*8 E,A,age,r0,radmm,ert,
     * tfac,area,delrad
      INTEGER, INTENT(IN) :: i,j
      REAL*8, INTENT(OUT) :: rads
      DATA E,A /26020.d0, 29100.d0/
c tsavg(i,j) surface air temperature
c snoage(k,i,j) k=1 ocean ice, k=2 land ice, k=3 land
c   (do these differ within a gridbox?) in days
c RN1 radius from previous timestep
c RADS snow grain radius
c
c Find the age of snow, I assume the age does not
c  vary within the gridbox so just take the max?
       age=DMAX1(snoage(1,i,j),snoage(2,i,j),snoage(3,i,j))
c Use Temperature to check if melting or non-melting snow
       IF (atmsrf%tsavg(i,j).le.tf) then
c Non-melting snow; distinguish between initial or
c  secondary growth rate
        IF (age.lt.13.5) then
c initial growth
         r0=50.
         radmm=r0+ (0.008d0+age)
         rads=radmm*1000.d0
         snosiz(i,j)=rads+r0
        ELSE
c secondary growth
         r0=150.d0
         ert = dexp(-E/(GASC*atmsrf%TSAVG(I,J)))
         tfac = a*ert
         area = (TFAC/DAYS_PER_YEAR) * (AGE-12.5d0)
         radmm=dsqrt(area/pi)
         delrad = radmm*1000.d0
         rads=delrad + r0
         snosiz(i,j)=rads
        ENDIF
       ELSE
c melting snow
        radmm=dsqrt(0.05d0)
        rads=snosiz(i,j)+(radmm*100.d0)
        snosiz(i,j)=rads
       ENDIF
       rads=DMIN1(rads,1000.d0)
       rads=DMAX1(rads,100.d0)
      RETURN
      END SUBROUTINE GRAINS
#endif  /* BC_ALB */
