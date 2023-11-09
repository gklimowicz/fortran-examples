      program vertflux
!@sum vertflux reads the ocean data saved from a GCM run,
!@+   calculates the ocean vertical flux (correcting for any small
!@+   imbalance) and writes the data onto a disk file.
!@+   RSI is interpolated in time from the climatological
!@+   ocean file for the run
!@+   Now produces qfluxes w and w/o ice dynamics
C****
C**** Input files : VFLX = Vert flux data from model run
C****               SICE = seaice file
C****               TOPO = topography
C**** Output files: XCORR = XCORR
C****               SRCOR = SRCOR-line for rundeck
C****               SNOW  = SNOW depth information
C****
!AOO use statements added for domain_decomp and dynamics to pull in
!AOO dynamically allocated arrays:  part 1 of 3
      use domain_decomp_1d, only : init_app, init_grid,grid, finish_app
      use domain_decomp_atm, only : readt_parallel
!AOO end of part 1 of 3
!!    use model_com, only : ioread
      use model_com, only : im,jm,lm
      USE STATIC_OCEAN
      USE DIAG_COM, only : koa
      USE SEAICE_COM, only : rsi,snowi
      USE SEAICE, only : ace1i
      USE FLUXES, only : sss,fwsim
      USE GEOM
      USE FILEMANAGER
      implicit none
      integer first_month, first_year, last_month, last_year, years
      integer months, monthe, year, month, itime1,
     *     last_day, kday, i, j, iu_TOPO, iu_VFLX, iu_XCORR, iu_SRCOR
     *     ,iu_SNOW,iok
      REAL*8 AVFX(IM,JM,2),GSR,GVFXSR(2),SYEAR,SYEARS
      real*8 VFSR, VF(2), XCORR(2)
      REAL*8 ASR(im,jm),AVFXSR(im,jm,2),oa(im,jm,koa)
      REAL*4 OAS(IM,JM,8), month_day(12)
      CHARACTER*80 TITLE(4),TITLE0, RunID, file_name
      character*4 month_name(12), tmonth, tyear
      character*2 tday
      data month_name /'JAN ','FEB ','MAR ','APR ','MAY ','JUN ',
     *                 'JUL ','AUG ','SEP ','OCT ','NOV ','DEC '/
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/
C****
      DATA SYEAR/31536000./
      DATA TITLE/'Uncorrected annual solar radiation (W/m**2)',
     2           'Vertical flux excluding solar radiation (W/m**2)',
     3           'Uncorrected annual vertical flux (W/m**2)',
     4           'Corrected annual vertical flux (W/m**2)'/
C****
C**** Contents of OA(I,J,K) saved from the GCM
C****
C****    K= 1  SNOWOI (instantaneous at noon GMT)
C****       2  FWSIM  (instantaneous at noon GMT)
C****       3  HSIT   (instantaneous at noon GMT)
C****       4  ENRGP  (integrated over the day)
C****       5  SRHDT  (for ocean, integrated over the day)
C****       6  TRHDT  (for ocean, integrated over the day)
C****       7  SHDT   (for ocean, integrated over the day)
C****       8  EVHDT  (for ocean, integrated over the day)
C****       9  TRHDT  (for ocean ice, integrated over the day)
C****      10  SHDT   (for ocean ice, integrated over the day)
C****      11  EVHDT  (for ocean ice, integrated over the day)
C****      12  SRHDT  (for ocean ice, integrated over the day)
C****
C**** Extra array needed for dealing with advected ice
C****      13  HCHSI  (HORIZ CONV SEA ICE ENRG, INTEGRATED OVER THE DAY)
C****
      off_line = .true. ! skip unneeded parts in OCLIM
!AOO calls to init routines for dynamically allocated arrays:part 2 of 3
      call init_app()
      call init_grid(grid,im,jm,lm)
      call alloc_drv()
!AOO end of part 2 of 3
      call getarg(1,RunID )
      call getarg(2,title0)
      read(title0,*) first_month
      call getarg(3,title0)
      read(title0,*) first_year
      call getarg(4,title0)
      read(title0,*) last_month
      call getarg(5,title0)
      read(title0,*) last_year
      years = last_year-first_year+1
      SYEARS = SYEAR*years
      call openunit("SICE",iu_SICE,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_SICE,"SICE",DM,1)
C****
C**** Calculate spherical geometry
C****
      call GEOM_B
C****
C**** Read in FOCEAN - ocean fraction
C****
      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_TOPO,"TOPO",focean,1) ! Ocean fraction
      call closeunit(iu_TOPO)
C****
C**** Zero out the vertical flux and its components
C****
      ASR = 0.
      AVFXSR = 0.
C**** initialise land array
      fland = 1.-focean
C****
      call openunit("OSST",iu_OSST,.true.,.true.)
CCC   call openunit("OCNML",iu_OCNML,.true.,.true.) ! not needed
      call openunit("SNOW",iu_SNOW,.true.,.false.)
      write(6,*) "Outputting Snow amount to unit ", iu_SNOW
C*
C**** define sea surface salinity (needed for OCLIM)
      sss(:,:)=sss0

C**** Loop over days of the year and read the data for each day
C****
      DO year=first_year,last_year
        months = 1
        monthe = 12
        if (year .eq. first_year) months = first_month
        if (year .eq. last_year)  monthe = last_month
        do month = months, monthe
          tmonth = month_name(month)
          write (tyear, '(i4.4)') year
          file_name='VFLXO'//trim(tmonth)//tyear
          call openunit(file_name,iu_VFLX,.true.,.true.)
          last_day = month_day(month)
          do kday = 1,last_day
            call READi8 (iu_VFLX, itime,OA,itime1,IM*JM*koa,iok)
            if(iok.gt.0) go to 555
            jmon = month
            jdate = kday

            kocean = 0
            CALL OCLIM (.true.)
C****
C**** Accumulate ASR, the daily flux of solar radiation,
C****     and AVFXSR, the vertical flux excluding solar radiation
C****
            DO J=1,JM
              DO I=1,IM
                IF(FOCEAN(I,J).LE.0.) CYCLE
                VFSR=(1.-RSI(I,J))*OA(I,J,5)+RSI(I,J)*OA(I,J,12)
                VF(1) = OA(I,J,4)
     *               + (1.-RSI(I,J))*(OA(I,J,6)+OA(I,J,7)+OA(I,J,8))
     *               + RSI(I,J)*(OA(I,J,9)+OA(I,J,10)+OA(I,J,11))
     *               + OA(I,J,13)
                VF(2) = OA(I,J,4)    ! no ice dyn term
     *               + (1.-RSI(I,J))*(OA(I,J,6)+OA(I,J,7)+OA(I,J,8))
     *               + RSI(I,J)*(OA(I,J,9)+OA(I,J,10)+OA(I,J,11))
                ASR(I,J)      = ASR(I,J)      + VFSR
                AVFXSR(I,J,:) = AVFXSR(I,J,:) + VF(:)
              END DO
            END DO
C*
C**** WRITE OUT OA(I,J,1) FOR THE LAST YEAR
C*
            IF (year .EQ. LAST_YEAR) THEN
              write (tday, '(i2.2)') kday
              TITLE0='SNOW AMOUNT ON SEA ICE, '//tyear//', '//tmonth//
     *             ', '//tday//', (kg/m^2)'
C*
              DO J = 1,JM
                DO I = 1,IM
                  SNOWI(I,J) = OA(I,J,1)
                END DO
              END DO
C*
              WRITE(iu_SNOW) TITLE0,SNOWI
            END IF
          end do
          call closeunit(iu_VFLX)
        end do
      end do
C*
      call closeunit(iu_SNOW)
C****
C**** DETERMINE THE SOLAR RADIATION CORRECTION FACTOR, XCORR, SO THAT
C**** THE CORRECTED GLOBAL VERTICAL FLUX IS EQUAL TO ZERO
C****
      GSR      = ASR(1,JM)     *FOCEAN(1,JM)*IM*DXYP(JM)
      GVFXSR(:)= AVFXSR(1,JM,:)*FOCEAN(1,JM)*IM*DXYP(JM)
      DO J=2,JM-1
      DO I=1,IM
        GSR      = GSR       + ASR(I,J)     *FOCEAN(I,J)*DXYP(J)
        GVFXSR(:)= GVFXSR(:) + AVFXSR(I,J,:)*FOCEAN(I,J)*DXYP(J)
      END DO
      END DO
      GSR    = GSR/years
      GVFXSR(:) = GVFXSR(:)/years
      XCORR(:) = -GVFXSR(:)/GSR
      WRITE(6,*)  'GSR,GVFXSR,XCORR=',GSR,GVFXSR,XCORR
      call openunit("XCORR",iu_XCORR,.true.,.false.)
      call openunit("SRCOR",iu_SRCOR,.false.,.false.)
      WRITE(iu_SRCOR,*) '  KOCEAN=1,SRCOR=',XCORR(1),', (with icedyn)'
      WRITE(iu_SRCOR,*) '  KOCEAN=1,SRCOR=',XCORR(2),', (no icedyn)'
      WRITE(iu_XCORR) XCORR
      call closeunit(iu_XCORR)
      call closeunit(iu_SRCOR)
C****
      DO J=1,JM
      DO I=1,IM
        AVFX(I,J,:) = (XCORR(:)*ASR(I,J)+AVFXSR(I,J,:))/SYEARS
      END DO
      END DO
C****
      do i=2,im
        AVFX(i,JM,:) = AVFX(1,JM,:)
      end do
C****
C**** Produce maps of ASR, AVFXSR, ASR+AVFXSR, and XCORR*ASR+AVFXSR
C****
      DO J=1,JM
      DO I=1,IM
        OAS(I,J,5) = FOCEAN(I,J)
        OAS(I,J,1) = ASR(I,J)/SYEARS
        OAS(I,J,2) = AVFXSR(I,J,1)/SYEARS
        OAS(I,J,3) = (ASR(I,J)+AVFXSR(I,J,1))/SYEARS
        OAS(I,J,4) = AVFX(I,J,1)
        OAS(I,J,6) = AVFXSR(I,J,2)/SYEARS
        OAS(I,J,7) = (ASR(I,J)+AVFXSR(I,J,2))/SYEARS
        OAS(I,J,8) = AVFX(I,J,2)
      END DO
      END DO
C**** fix diagnostic polar boxes
      DO I=2,IM
        OAS(I, 1,:) = OAS(1, 1,:)
        OAS(I,JM,:) = OAS(1,JM,:)
      END DO
      CALL MAP1 (IM,JM,0,TITLE(1),OAS(1,1,1),OAS(1,1,5),1.,0.,0)
      CALL MAP1 (IM,JM,0,TITLE(2),OAS(1,1,2),OAS(1,1,5),1.,0.,0)
      CALL MAP1 (IM,JM,0,TITLE(3),OAS(1,1,3),OAS(1,1,5),1.,0.,0)
      CALL MAP1 (IM,JM,0,TITLE(4),OAS(1,1,4),OAS(1,1,5),1.,0.,0)
C**** no ice dyn
      do i=2,4
        TITLE(i)=trim(TITLE(i))//" no IceD"
      end do
      CALL MAP1 (IM,JM,0,TITLE(2),OAS(1,1,6),OAS(1,1,5),1.,0.,0)
      CALL MAP1 (IM,JM,0,TITLE(3),OAS(1,1,7),OAS(1,1,5),1.,0.,0)
      CALL MAP1 (IM,JM,0,TITLE(4),OAS(1,1,8),OAS(1,1,5),1.,0.,0)
C****
!AOO not sure if this is needed, but just in case ...  part 3 of 3
      call finish_app()
!AOO end of part 3 of 3
      STOP
 555  write (0,*) ' Error: Premature end of file ',file_name
      call exit_rc (11)
      END

