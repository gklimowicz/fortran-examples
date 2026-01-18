!**** WATERG.F90     Model E Global Water budget     2018/05/09
!**** Fortran Compile into Executable:  FCE90NETCDF WATERG.F90

!**** Name            IJ_Diag       SourceFile        Subroutine       LName
!**** H20_from_CH4    IJ_H2OCH4     RAD_DRV.f         DAILY_CH4OX      WATER DERIVED FROM CH4 OXIDATION IN STRATOSPHERE
!**** prec            IJ_PREC       CLOUDS2_DRV.F90   CONDSE           PRECIPITATION
!**** pr_grnd         IJ_PRECGR     CLOUDS2_DRV.F90   CONDSE           PRECIPITATION OVER EARTH
!**** pr_lndice       IJ_PRECLI     CLOUDS2_DRV.F90   CONDSE           PRECIPITATION OVER LAND ICE
!**** pr_oice         IJ_PRECSI     CLOUDS2_DRV.F90   CONDSE           PRECIPITATION OVER SEA ICE
!**** pr_oocn         IJ_PRECOO     CLOUDS2_DRV.F90   CONDSE           PRECIPITATION OVER OPEN OCEAN
!**** evap            IJ_EVAP       SURFACE.f         SURFACE_DIAG1    EVAPORATION
!**** evap_land       IJ_EVAPE      GHY_DRV.f         GROUND_E         SOIL EVAPORATION
!**** evap_oice       IJ_EVAPI      SURFACE.f         SURFACE_DIAG1a   OCEAN ICE EVAPORATION
!**** evap_ocn        IJ_EVAPO      SURFACE.f         SURFACE_DAIG1a   OPEN WATER EVAPORATION
!**** evap_lndice     IJ_EVAPLI     SURFACE.f         SURFACE_DIAG1a   LAND ICE EVAPORATION
!**** runoff_soil     IJ_RUNE       GHY_DRV.f         GHY_DIAG         GROUND RUNOFF OVER SOIL
!**** topmlt_oice     IJ_SITOPMLT   SEAICE_DRV.f      SI_DIAGS         SEA ICE SURFACE MELT RATE
!**** runoff_lndice   IJ_RUNLI      SURFACE.f         SURFACE_DIAG3    SURFACE RUNOFF OVER LAND ICE
!**** runoff_ugrnd    IJ_ARUNU      GHY_DRV.f         GHY_DIAG         UNDERGROUND RUNOFF OVER SOIL
!**** irrig_gw        IJ_IRRGW      LAKES.f           IRRIG_LK         IRRIGATION WATER FROM EXTERNAL SOURCE (GRNDWATER)
!**** mwl_irrigate    IJ_MWLIR      LAKES.f           IRROG_LK         MASS OF LAKE/RIVER WATER USED FOR IRRIGATION
!**** MRVR            IJ_MRVR       LAKES.f           RIVERF           Mass Inflow by Rivers
!**** msnflood        IJ_MSNFLOOD   SEAICE_DRV.f      SI_DIAGS         ICE MASS FROZEN by SNOW FLOOD
!**** grlat_oice      IJ_SIGRLT     SEAICE_DRV.f      SI_DIAGS         SEA ICE LATERAL GROWTH RATE
!**** botmlt_oice     IJ_SIBOTMLT   SEAICE_DRV.f      SI_DIAGS         SEA ICE BASAL MELT RATE    
!**** grcong_oice     IJ_SIGRCG     SEAICE_DRV.f      SI_DIAGS         SEA ICE CONGEALATION GROWTH RATE
!**** grfraz_oice     IJ_SIGRFR     SEAICE_DRV.f      SI_DIAGS         SEA ICE FRAZIL GROWTH RATE
!**** mlktogr         IJ_MLKtoGR    LAKES.f           DAILY_LAKE       MASS of EXPANDING LAKE SATURATES GROUND
!**** impm_gr         IJ_IMPMGR     GHY_DRV.f         REMOVE_EXTRA..   IMPLICIT MASS FLUX over GROUND
!**** impm_ki         IJ_IMPMKI     LAKES.f           DAILY_LAKE       IMPLICIT MASS FLUX over LAKE ICE
!**** impm_lanice     IJ_IMPMLI     SURFACE.f         SURFACE_DIAG3    IMPLICIT MASS FLUX over LAND ICE
!**** MICB

      Implicit None
      Integer,Parameter :: IMA=144,JMA= 90, IMO=288,JMO=180, &                                                                            
         DAYSzM(13) = (/ 31,28,31, 30,31,30, 31,31,30, 31,30,31, 365 /)
      Real*8,Parameter :: RADIUS=6371000, SDAY=86400
      Character*3,Parameter :: &
         CMONTH(13) = (/ 'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC', 'ANN' /)
!****
      Integer :: I,J, NARGS, MRSF0,MRSFM,MACC, YEAR0,YEARM, COLON1,COLON2,LENARG,SLASH, IOSTAT,IDAIJ,IDRSF,IDQ
      Logical :: QEXIST
      Real*4  :: aFOCEAN(IMA,JMA), FGICE(IMA,JMA),oFOCEAN(IMO,JMO), &
                  FEARTH(IMA,JMA), FLAKE(IMA,JMA), FWICE(IMA,JMA), &
                  OXIDAT(IMA,JMA), &
                  PRECAT(IMA,JMA),PRECGR(IMA,JMA),PRECGI(IMA,JMA),PRECWI(IMA,JMA),PRECWL(IMA,JMA), &
                  RUNPGR(IMA,JMA),RUNPGI(IMA,JMA),RUNPWI(IMA,JMA), &
                  EVAPAT(IMA,JMA),EVAPGR(IMA,JMA),EVAPGI(IMA,JMA),EVAPWI(IMA,JMA),EVAPWL(IMA,JMA), &
                  RUNSGR(IMA,JMA),RUNSGI(IMA,JMA),RUNSWT(IMA,JMA), &
                  RUNUGR(IMA,JMA), &
                  AQFRGR(IMA,JMA), &
                  IRIGGR(IMA,JMA), &
                  RVRFLO(IMA,JMA), &
                  MLATWT(IMA,JMA), &
                  SNFLWT(IMA,JMA), &
                  MBOTWT(IMA,JMA), &
                  FREZWT(IMA,JMA), &
                  FRAZWT(IMA,JMA), &
                  LKtoGR(IMA,JMA), &
                  FORBGR(IMA,JMA),FORBLI(IMA,JMA),FORBGI(IMA,JMA), &
                  MELBBI(IMA,JMA)
      Real*8  :: aDXYP(JMA),oDXYP(JMO),AREAG, aXYP(IMA,JMA),oXYP(IMO,JMO), &
                 AGR(IMA,JMA),ALI(IMA,JMA),ALL(IMA,JMA),ALK(IMA,JMA), &
                 AGI(IMA,JMA),ASI(IMA,JMA),AOL(IMA,JMA),AOC(IMA,JMA), &
                 RES0(9),RESM(9),ERRO(9), DAYS, &
                 TRAN(9),OXID(9),PREC(9),RUNP(9),EVAP(9),RUNS(9),RUNU(9),AQFR(9),IRIG(9),RVRF(9), &
                 SNFL(9),MLAT(9),MBOT(9),FREZ(9),FRAZ(9),LKGR(9),FORB(9),MELB(9)
      Character*4   :: YRSF0,YRSFM
      Character*9   :: YRMO0,YRMOM,YACC
      Character*10  :: DATE0,DATEM
      Character*80  :: TITLE
      Character*256 :: OPDIR,IFDIR, COMMND, RUNARG, RSFFILE0,RSFFILEM,AIJFILE,ACCFILE, AIJSUB,RSFSUB, RUN,PATH, &
                       FILEIN1,FILEIN2,FILEIN3,FILEIN4, RSFPATH0,RSFPATHM,AIJPATH,ACCPATH,ACCDIR, ZFILEA,ZFILEO
      Include 'netcdf.inc'

      Call GetEnv ('IFDIR',IFDIR)  ;  If (IFDIR .eq. ' ') IFDIR = '/discover/nobackup/projects/giss/prod_input_files'
      Call GetEnv ('OPDIR',OPDIR)  ;  If (OPDIR .eq. ' ') OPDIR = '/discover/nobackup/projects/giss/prod_runs'
      NARGS = IArgC()
      If (NARGS < 2)  GoTo 800

!**** Determine time period from 2-nd (and 3-rd) command line arguments
!**** Determine initial rsf date
      Call GetArg (2,YRMO0)
      YRSF0 = YRMO0(1:4)
      MRSF0 = 1
      If (Len_Trim(YRMO0) > 5)  Read (YRMO0(6:7),*) MRSF0
      Write (DATE0,902) YRSF0,MRSF0
      If (NARGS >= 3)  GoTo 10
!**** Determine final rsf date when 2-nd time is not provided
      If (Len_Trim(YRMO0) > 5) &
         Then  ;  DAYS  = DAYSzM(MRSF0)  !  individual month
                  YRSFM = YRSF0
                  MRSFM = MRSF0 + 1
                  If (MRSF0 == 12)  Then
                     Read  (YRSF0,*) YEAR0
                     Write (YRSFM,903) YEAR0+1
                     MRSFM = 1  ;  EndIf
         Else  ;  DAYS = 365  !  individual year
                  Read  (YRSF0,*) YEAR0
                  Write (YRSFM,903) YEAR0+1
                  MRSFM = 1  ;  EndIf
      GoTo 20
!**** Determine final rsf date when 2-nd time is provided
   10 Call GetArg (3,YRMOM)
      If (Len_Trim(YRMO0) > 4 .or. Len_Trim(YRMOM) > 4)  GoTo 801
      Read  (YRSF0,*) YEAR0
      Read  (YRMOM,*) YEARM  ;  YEARM = YEARM+1
      Write (YRSFM,903) YEARM
      DAYS  = 365 * (YEARM-YEAR0)
      MRSFM = 1
   20 Write (DATEM,902) YRSFM,MRSFM
!**** Determine aij dates
      YACC = YRSF0
      If (NARGS > 2)  YACC = YRSF0 // '-' // YRMOM
      MACC = 13
      If (Len_Trim(YRMO0) > 5)  MACC = MRSF0

!**** Determine RUN, directory RUNDIR, and subdirectories AIJSUB and RSFSUB
      Call GetArg (1,RUNARG)
      COLON1 = Index (RUNARG, ':')
      COLON2 = Index (RUNARG, ':', BACK=.True.)
      LENARG = Len_Trim (RUNARG)
!**** Determine RSF subdirectory
      RSFSUB = ' '
      If (COLON1 < COLON2 .and. COLON2 < LENARG)  RSFSUB = RUNARG(COLON2+1:LENARG) // '/'
!**** Determine AIJ subdirectory
      AIJSUB = ' '
      If (0 < COLON1 .and. COLON1 < LENARG) Then
         If (COLON1+1 < COLON2) &
            Then  ;  AIJSUB = RUNARG(COLON1+1:COLON2-1) // '/'
            Else  ;  AIJSUB = RUNARG(COLON1+1:LENARG) // '/'  ;  EndIf  ;  EndIf
!**** Determine RUN
      If (COLON1 == 0)  COLON1 = LENARG + 1 
      SLASH = Index (RUNARG(1:COLON1-1), '/', BACK=.True.)
      RUN   = RUNARG(SLASH+1:COLON1-1)
      PATH  = RUNARG(1:COLON1-1)
!**** Determine file names
      RSFFILE0 = '1' // CMONTH(MRSF0) // YRSF0 // '.rsf' // Trim(RUN) // '.nc'
      RSFFILEM = '1' // CMONTH(MRSFM) // YRSFM // '.rsf' // Trim(RUN) // '.nc'
      AIJFILE  =    CMONTH(MACC) // Trim(YACC) // '.aij' // Trim(RUN) // '.nc'
      ACCFILE  =    CMONTH(MACC) // Trim(YACC) // '.acc' // Trim(RUN) // '.nc'
!**** Determine RSF and AIJ files when PresetWorkingDirectory is the RUNDIR
      FILEIN1 = Trim(RSFSUB) // RSFFILE0
      Inquire (File=FILEIN1, Exist=QEXIST)
      If (QEXIST) Then
         RSFPATH0 = Trim(RSFSUB) // RSFFILE0
         RSFPATHM = Trim(RSFSUB) // RSFFILEM
         AIJPATH  = Trim(AIJSUB) // AIJFILE
         ACCPATH  = Trim(AIJSUB) // ACCFILE
         ACCDIR   =      AIJSUB
         GoTo 50  ;  EndIf
!**** Determine RSF and AIJ files when RUNDIR = RUN and RUN is a subdirectory of PresentWorkingDirectory
      FILEIN2 = Trim(RUN) // '/' // FILEIN1
      Inquire (File=FILEIN2, Exist=QEXIST)
      If (QEXIST) Then
         RSFPATH0 = Trim(RUN) // '/' // Trim(RSFSUB) // RSFFILE0
         RSFPATHM = Trim(RUN) // '/' // Trim(RSFSUB) // RSFFILEM
         AIJPATH  = Trim(RUN) // '/' // Trim(AIJSUB) // AIJFILE
         ACCPATH  = Trim(RUN) // '/' // Trim(AIJSUB) // ACCFILE
         ACCDIR   = Trim(RUN) // '/' //      AIJSUB
         GoTo 50  ;  EndIf
!**** Determine RSF and AIJ files when RUN is symbolicly lined as subdirectory of $OPDIR
      FILEIN3 = Trim(OPDIR) // '/' // FILEIN2
      Inquire (File=FILEIN3, Exist=QEXIST)
      If (QEXIST) Then
         RSFPATH0 = Trim(OPDIR) // '/' // Trim(RUN) // '/' // Trim(RSFSUB) // RSFFILE0
         RSFPATHM = Trim(OPDIR) // '/' // Trim(RUN) // '/' // Trim(RSFSUB) // RSFFILEM
         AIJPATH  = Trim(OPDIR) // '/' // Trim(RUN) // '/' // Trim(AIJSUB) // AIJFILE 
         ACCPATH  = Trim(OPDIR) // '/' // Trim(RUN) // '/' // Trim(AIJSUB) // ACCFILE 
         ACCDIR   = Trim(OPDIR) // '/' // Trim(RUN) // '/' //      AIJSUB
         GoTo 50  ;  EndIf
!**** Determine RSF and AIJ files when RUNARG contains the entire path, RUNARG(1:COLON1-1)
      FILEIN4 = Trim(PATH) // '/' // FILEIN1
      Inquire (File=FILEIN4, Exist=QEXIST)
      If (QEXIST) Then
         RSFPATH0 = Trim(PATH) // '/' // Trim(RSFSUB) // RSFFILE0
         RSFPATHM = Trim(PATH) // '/' // Trim(RSFSUB) // RSFFILEM
         AIJPATH  = Trim(PATH) // '/' // Trim(AIJSUB) // AIJFILE
         ACCPATH  = Trim(PATH) // '/' // Trim(AIJSUB) // ACCFILE
         ACCDIR   = Trim(PATH) // '/' //      AIJSUB
         GoTo 50  ;  EndIf
      GoTo 803

!**** Open NetCDF file for aij diagnostics
!**** Calculate areas for various reservoirs
   50 Inquire (File=AIJPATH, Exist=QEXIST)
      If (.not.QEXIST) Then
         Inquire (File=ACCPATH, Exist=QEXIST)
         If (.not.QEXIST)  GoTo 805
!**** If .aij file does not exist, but .acc file does, create .aij using scaleacc
	 If (ACCDIR .ne. ' ') &
            Then  ;  COMMND = 'cd ' // Trim(ACCDIR) // ' ; scaleacc ' // Trim(ACCFILE) // ' aij'
            Else  ;  COMMND = 'scaleacc ' // Trim(ACCFILE) // ' aij'  ;  EndIf
         Write (0,900) 'Executing:  ' // Trim(COMMND)
         Call System (COMMND)  ;  EndIf
!**** Read atmosphere areas and fractions for .aij file
      IOSTAT = NF_OPEN (AIJPATH, 0, IDAIJ)  ;  If(IOSTAT /= 0) GoTo 806
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'ocnfr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, aFOCEAN)  !  (%)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'soilfr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FEARTH)  !  (%)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'lakefr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FLAKE)  !  (%)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'oicefr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FWICE)  !  (%)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'landicefr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FGICE)  !  (%)

!**** Read input file Z2HX2fromZ1QX1N obtaining aFOCEAN and FLICE
      ZFILEA = Trim(IFDIR) // '/Z2HX2fromZ1QX1N.BS1'
      Open (11, File=ZFILEA, Form='Unformatted', Status='Old', Err=807)
      Read (11) TITLE,aFOCEAN
      Read (11)
      Read (11)
      Read (11) TITLE,FGICE
      Close(11)
!**** Read input file Z1QX1N obtaining oFOCEAN
      ZFILEO = Trim(IFDIR) // '/Z1QX1N.BS1'
      Open (11, File=ZFILEO, Form='Unformatted', Status='Old', Err=808)
      Read (11) TITLE,oFOCEAN
      Close(11)

!**** Determine ocean grid cell area, multiply polar values times IM
      Call SGEOM (IMA,JMA,RADIUS, aDXYP,AREAG)
      Call SGEOM (IMO,JMO,RADIUS, oDXYP,AREAG)
      Do 60 J=1,JMA
   60 aXYP(:,J) = aDXYP(J)
      Do 70 J=1,JMO
   70 oXYP(:,J) = oDXYP(J)
!**** Note that A??(2:IMA,1) = A??(2:IMA,JMA) = 0  and  A??(1,1) & A??(1,JMA) contain IMA*[area of polar triangle cell]
      aXYP(1,1)   = IMA*aDXYP(1)    ;  aXYP(2:IMA,1)   = 0
      aXYP(1,JMA) = IMA*aDXYP(JMA)  ;  aXYP(2:IMA,JMA) = 0
      oXYP(1,1)   = IMO*oDXYP(1)    ;  oXYP(2:IMO,1)   = 0
      oXYP(1,JMO) = IMO*oDXYP(JMO)  ;  oXYP(2:IMO,JMO) = 0

!**** Compute atmospheric area fractions (m^2)
!**** Note that A??(2:IMA,1) = A??(2:IMA,JMA) = 0  and  A??(1,1) & A??(1,JMA) contain IMA*[area of polar triangle cell]
      Do 80 J=1,JMA  ;  Do 80 I=1,IMA
      AGR(I,J) = aXYP(I,J) * FEARTH(I,J) / 100
      AGI(I,J) = aXYP(I,J) *  FGICE(I,J)
      If (aFOCEAN(I,J) == 0)  Then
         ALI(I,J) = aXYP(I,J) *  FWICE(I,J) / 100
         ALL(I,J) = aXYP(I,J) * (FLAKE(I,J) - FWICE(I,J)) / 100
         ALK(I,J) = aXYP(I,J) *  FLAKE(I,J) / 100
         ASI(I,J) = 0
         AOL(I,J) = 0
         AOC(I,J) = 0
      Else
         ASI(I,J) = aXYP(I,J) *    FWICE(I,J) / 100
         AOL(I,J) = aXYP(I,J) * (aFOCEAN(I,J) - FWICE(I,J) / 100)
         AOC(I,J) = aXYP(I,J) *  aFOCEAN(I,J)
         ALI(I,J) = 0
         ALL(I,J) = 0
         ALK(I,J) = 0  ;  EndIf
      If (Abs(FLAKE(I,J) + FGICE(I,J)*100 + FEARTH(I,J) + aFOCEAN(I,J)*100 - 100) < 1e-4)  GoTo 80
      Write (0,908) I,J, aFOCEAN(I,J)*100, FEARTH(I,J), FGICE(I,J)*100, FLAKE(I,J), &
                           FLAKE(I,J) + FGICE(I,J)*100 + FEARTH(I,J) + aFOCEAN(I,J)*100
   80 Continue

!**** Open and write titles to output file WATERG.TXT
      Open  (3, File='WATERG.TXT')
      Write (3,900) 'ZFILEA = ' // Trim(ZFILEA)
      Write (3,900) 'ZFILEO = ' // Trim(ZFILEO)
      Write (3,910) Trim(RUN)

!**** Open NetCdF rsf file for YEAR0
!**** Calculate and write reservoir data for YEAR0
      IOSTAT = NF_OPEN (RSFPATH0, 0, IDRSF)  ;  If(IOSTAT /= 0) GoTo 810
      Call WATERRES (IDRSF,aXYP,oXYP,AREAG, RES0)
      IOSTAT = NF_CLOSE (IDRSF)
      Write (3,911) DATE0,RES0
      Write (3,911)
!****
!**** Horizontal transport
!****

!****
!**** Methane Oxidation
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'H2O_from_CH4', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, OXIDAT)  !  (10^-6 mm/day)
      OXID(1) = Sum (aXYP(:,:)*OXIDAT(:,:)) * 1d-6  ;  OXID(2:8) = 0
      OXID(9) = OXID(1)
      OXID(:) = OXID(:) * DAYS / AREAG
      Write (3,930) 'CH4 Oxidat',OXID(1),OXID(9)
!****
!**** Precipitation
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'prec', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, PRECAT)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'pr_grnd', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, PRECGR)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'pr_lndice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, PRECGI)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'pr_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, PRECWI)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'pr_oocn', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, PRECWL)  !  (mm/day)
      PREC(1) = - Sum (aXYP(:,:)*PRECAT(:,:))
      PREC(2) = Sum (AGR(:,:)*PRECGR(:,:))
      PREC(3) = Sum (ALI(:,:)*PRECWI(:,:))
      PREC(4) = Sum (ALL(:,:)*PRECWL(:,:))
      PREC(5) = Sum (AGI(:,:)*PRECGI(:,:))
      PREC(6) = 0
      PREC(7) = Sum (ASI(:,:)*PRECWI(:,:))
      PREC(8) = Sum (AOL(:,:)*PRECWL(:,:))
      PREC(9) = Sum (PREC(1:8))
      PREC(:) = PREC(:) * DAYS / AREAG
      Write (3,931) 'Precipitat',PREC(1:5),PREC(7:9)
!****
!**** Surface Runoff from Precipitation
!****
!     IOSTAT = NF_INQ_VARID    (IDAIJ, 'runp_gr', IDQ)
!     IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNPGR)  !  (mm/day)
!     IOSTAT = NF_INQ_VARID    (IDAIJ, 'runp_gi', IDQ)
!     IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNPGI)  !  (mm/day)
!     IOSTAT = NF_INQ_VARID    (IDAIJ, 'runp_wi', IDQ)
!     IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNPWI)  !  (mm/day)
!     RUNP(:) = 0
!     RUNP(2) = - Sum (AGR(:,:)*RUNPGR(:,:))
!     RUNP(3) = - Sum (ALI(:,:)*RUNPWI(:,:))
!     RUNP(5) = - Sum (AGI(:,:)*RUNPGI(:,:))
!     RUNP(7) = - Sum (ASI(:,:)*RUNPWI(:,:))
!     RUNP(4) = - RUNP(2) - RUNP(3) - RUNP(5)
!     RUNP(8) = - RUNP(7)
!     RUNP(9) = Sum (PREC(1:8))
!     RUNP(:) = PREC(:) * DAYS / AREAG
!     Write (3,932) 'PrecRunoff',RUNP(2:5),RUNP(7:9)
!****
!**** Evaporation
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'evap', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, EVAPAT)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'evap_land', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, EVAPGR)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'evap_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, EVAPWI)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'evap_ocn', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, EVAPWL)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'evap_lndice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, EVAPGI)  !  (mm/day)
      EVAP(1) =  Sum (aXYP(:,:)*EVAPAT(:,:))
      EVAP(2) = - Sum (AGR(:,:)*EVAPGR(:,:))
      EVAP(3) = - Sum (ALI(:,:)*EVAPWI(:,:))
      EVAP(4) = - Sum (ALL(:,:)*EVAPWL(:,:))
      EVAP(5) = - Sum (AGI(:,:)*EVAPGI(:,:))
      EVAP(6) = 0
      EVAP(7) = - Sum (ASI(:,:)*EVAPWI(:,:))
      EVAP(8) = - Sum (AOL(:,:)*EVAPWL(:,:))
      EVAP(9) = Sum (EVAP(1:8))
      EVAP(:) = EVAP(:) * DAYS / AREAG
      Write (3,931) 'Evaporatio',EVAP(1:5),EVAP(7:9)
!****
!**** Runoff from Surface Heating
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'runoff_soil', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNSGR)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'topmlt_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNSWT)  !  (mm/day)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'runoff_lndice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNSGI)  !  (mm/day)
      RUNS(:) = 0
      RUNS(2) = - Sum (AGR(:,:)*RUNSGR(:,:))
      RUNS(3) = - Sum (ALK(:,:)*RUNSWT(:,:)) * SDAY
      RUNS(5) = - Sum (AGI(:,:)*RUNSGI(:,:))
      RUNS(4) = - RUNS(2) - RUNS(3) - RUNS(5)
      RUNS(7) = - Sum (AOC(:,:)*RUNSWT(:,:)) * SDAY
      RUNS(8) = - RUNS(7)
      RUNS(9) = Sum (RUNS(1:8))
      RUNS(:) = RUNS(:) * DAYS / AREAG
      Write (3,932) 'SurfRunoff',RUNS(2:5),RUNS(7:9)
!****
!**** Underground Runoff
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'runoff_ugrnd', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RUNUGR)  !  (mm/day)
      RUNU(:) = 0
      RUNU(2) = - Sum (AGR(:,:)*RUNUGR(:,:))
      RUNU(4) = - RUNU(2)
      RUNU(9) = Sum (RUNU(1:8))
      RUNU(:) = RUNU(:) * DAYS / AREAG
      Write (3,935) 'UndergrRun',RUNU(2),RUNU(4),RUNU(9)
!****
!**** Aquifers
!****
      AQFRGR(:,:) = 0
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'irrig_gw', IDQ)
      If (IOSTAT /= 0) GoTo 370
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, AQFRGR)  !  (mm/day)
      AQFR(:) = 0
      AQFR(2) = Sum (aXYP(:,:)*AQFRGR(:,:))
      AQFR(9) = Sum (AQFR(1:8))
      AQFR(:) = AQFR(:) * DAYS / AREAG
      Write (3,936) 'Aquifers  ',AQFR(2),AQFR(9)
!****
!**** Lake irrigation = Total irrigation - aquifers
!****
  370 IOSTAT = NF_INQ_VARID    (IDAIJ, 'mwl_irrigate', IDQ)
      If (IOSTAT /= 0) GoTo 380
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, IRIGGR)  !  (kg/s)
      IRIG(:) = 0
      IRIG(2) = Sum (IRIGGR(:,:))
      IRIG(4) = - IRIG(2)
      IRIG(9) = Sum (IRIG(1:8))
      IRIG(:) = IRIG(:) * DAYS *SDAY / AREAG
      Write (3,935) 'Irrigation',IRIG(2),IRIG(4),IRIG(9)
!****
!**** River Flow: assumes no ocean near SP, no land near NP
!****
  380 IOSTAT = NF_INQ_VARID    (IDAIJ, 'MRVR', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, RVRFLO)  !  (10^5 kg/s)
      RVRF(:) = 0
      Do 381 J=2,JMA-1
      Do 381 I=1,IMA
  381 If (aFOCEAN(I,J) > 0)  RVRF(8) = RVRF(8) + RVRFLO(I,J)
      RVRF(4) = - RVRF(8)
      RVRF(9) = Sum (RVRF(1:8))
      RVRF(:) = RVRF(:) * 1d5 * DAYS * SDAY / AREAG
      Write (3,938) 'River Flow',RVRF(4),RVRF(8:9)
!****
!**** Ice freexing by snow flooding
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'msnflood', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, SNFLWT)  !  (kg/s*m^2)
      SNFL(:) = 0
      SNFL(3) = Sum (ALK(:,:)*SNFLWT(:,:))
      SNFL(4) = - SNFL(3)
      SNFL(7) = Sum (AOC(:,:)*SNFLWT(:,:))
      SNFL(8) = - SNFL(7)
      SNFL(9) = Sum (SNFL(1:8))
      SNFL(:) = SNFL(:) * DAYS * SDAY / AREAG
      Write (3,939) 'Snow Flood',SNFL(3:4),SNFL(7:9)
!****
!**** Lateral Melting of Water Ice
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'grlat_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, MLATWT)  !  (kg/s*m^2)
      MLAT(:) = 0
      MLAT(3) = Sum (ALK(:,:)*MLATWT(:,:))
      MLAT(4) = - MLAT(3)
      MLAT(7) = Sum (AOC(:,:)*MLATWT(:,:))
      MLAT(8) = - MLAT(7)
      MLAT(9) = Sum (MLAT(1:8))
      MLAT(:) = MLAT(:) * DAYS * SDAY / AREAG
      Write (3,939) 'LateraMelt',MLAT(3:4),MLAT(7:9)
!****
!**** Bottom Melting of Water Ice
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'botmlt_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, MBOTWT)  !  (kg/s*m^2)
      MBOT(:) = 0
      MBOT(3) = - Sum (ALK(:,:)*MBOTWT(:,:))
      MBOT(4) = - MBOT(3)
      MBOT(7) = - Sum (AOC(:,:)*MBOTWT(:,:))
      MBOT(8) = - MBOT(7)
      MBOT(9) = Sum (MBOT(1:8))
      MBOT(:) = MBOT(:) * DAYS * SDAY / AREAG
      Write (3,939) 'BottomMelt',MBOT(3:4),MBOT(7:9)
!****
!**** Bottom Freezing of Water Ice
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'grcong_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FREZWT)  !  (kg/s*m^2)
      FREZ(:) = 0
      FREZ(3) = Sum (ALK(:,:)*FREZWT(:,:))
      FREZ(4) = - FREZ(3)
      FREZ(7) = Sum (AOC(:,:)*FREZWT(:,:))
      FREZ(8) = - FREZ(7)
      FREZ(9) = Sum (FREZ(1:8))
      FREZ(:) = FREZ(:) * DAYS * SDAY / AREAG
      Write (3,939) 'BottomFrez',FREZ(3:4),FREZ(7:9)
!****
!**** Frazil Freezing of Water Ice
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'grfraz_oice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FRAZWT)  !  (kg/s*m^2)
      FRAZ(:) = 0
      FRAZ(3) = Sum (ALK(:,:)*FRAZWT(:,:))
      FRAZ(4) = - FRAZ(3)
      FRAZ(7) = Sum (AOC(:,:)*FRAZWT(:,:))
      FRAZ(8) = - FRAZ(7)
      FRAZ(9) = Sum (FRAZ(1:8))
      FRAZ(:) = FRAZ(:) * DAYS * SDAY / AREAG
      Write (3,939) 'FrazilFrez',FRAZ(3:4),FRAZ(7:9)
!****
!**** Expanding Lake to Ground
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'mlktogr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, LKtoGR)  !  (kg/s*m^2)
      LKGR(:) = 0
      LKGR(2) = Sum (aXYP(:,:)*LKtoGR(:,:))
      LKGR(4) = - LKGR(2)
      LKGR(9) = Sum (LKGR(1:8))
      LKGR(:) = LKGR(:) * DAYS * SDAY / AREAG
      Write (3,944) 'ExpandLake',LKGR(2),LKGR(4),LKGR(9)
!****
!**** Formation of Berg Ice
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'impm_gr', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FORBGR)  !  (kg/s*m^2)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'impm_ki', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FORBLI)  !  (kg/s*m^2)
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'impm_lndice', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, FORBGI)  !  (kg/s*m^2)
      FORB(:) = 0
      FORB(2) = - Sum (AGR(:,:)*FORBGR(:,:))
      FORB(3) = - Sum (ALI(:,:)*FORBLI(:,:))
      FORB(5) = - Sum (AGI(:,:)*FORBGI(:,:))
      FORB(6) = - FORB(2) - FORB(3) - FORB(5)
      FORB(9) = Sum (FORB(1:8))
      FORB(:) = FORB(:) * DAYS * SDAY / AREAG
      Write (3,945) 'FormBergIc',FORB(2:3),FORB(5:6),FORB(9)
!****
!**** Melting of Berg Ice
!****
      IOSTAT = NF_INQ_VARID    (IDAIJ, 'MICB', IDQ)
      IOSTAT = NF_GET_VAR_REAL (IDAIJ, IDQ, MELBBI)  !  (10^5 kg/s)
      MELB(:) = 0
      MELB(6) = - Sum (MELBBI(:,:))
      MELB(8) = - MELB(6)
      MELB(9) = Sum (MELB(1:8))
      MELB(:) = MELB(:) * 1d5 * DAYS * SDAY / AREAG
      Write (3,961) 'MeltBergIc',MELB(6),MELB(8:9)
      IOSTAT = NF_CLOSE (IDAIJ)
!****
!**** Open NetCdF rsf file for YEAR0
!**** Calculate and write reservoir data for YEARM
!****
      IOSTAT = NF_OPEN (RSFPATHM, 0, IDRSF)  ;  If(IOSTAT /= 0) GoTo 880
      Call WATERRES (IDRSF,aXYP,oXYP,AREAG, RESM)
      IOSTAT = NF_CLOSE (IDRSF)
      Write (3,911)
      Write (3,911) DATEM,RESM
!****
!**** Calculate error: ERRO = RESM - RES0 - Sum(FLUX)
!****
      ERRO(:) = RESM(:) - RES0(:) - OXID(:) - PREC(:) - RUNP(:) - EVAP(:) - RUNS(:) - RUNU(:) - AQFR(:) - IRIG(:) - RVRF(:) &
                        - SNFL(:) - MLAT(:) - MBOT(:) - FREZ(:) - FRAZ(:) - LKGR(:) - FORB(:) - MELB(:)
      Write (3,911)
      Write (3,911) 'Error     ',ERRO(1:9)
!****
      Call System ('cat WATERG.TXT')
      GoTo 999
!****
  800 Write (6,900) '              Run  Year/M             Dates Used'     
      Write (6,900) '              ===========             =========='     
      Write (6,900) 'Usage: WATERG RUN  1850/1             JAN1850'
      Write (6,900) '       WATERG RUN:AIJ  1850           ANN1850'           
      Write (6,900) '       WATERG RUN::RSF  1850 1859     ANN1850-1859'                  
      Write (6,900) '       WATERG RUN:AIJ:RSF  1850/7     JUL1850'       
      Write (6,900) '       WATERG path/RUN  1850          ANN1850'                                
      Write (6,900) ' '
      Write (6,900) '       RUN directory is RUN, pwd/RUN, $OPDIR/RUN, or path/RUN       2018/05/10'
      Write (6,900) '       AIJ is a subdirectory of RUN directory where .aij (or .acc) file exists'
      Write (6,900) '       RSF is a subdirectory of RUN directory where both .rsf files exists'
      Write (6,900) '       WATERG requires .aij file and .rsf files at both ends of time period' 
      Write (6,900) '       If .aij file does not exist but .acc does, .aij is created using scaleacc'
      Write (6,900) '       WATERG does not execute sumfiles to create annual and multi-year .aij files'
      Write (6,900) '       Stored output file in pwd: WATERG.TXT'
      Write (6,900) '       $IFDIR = ' // Trim(IFDIR)
      Write (6,900) '       $OPDIR = ' // Trim(OPDIR)     
      GoTo 999
  801 Write (0,*) 'Unacceptable command line times: ',YRMO0,YRMOM   ;  Stop 801
  803 Write (0,900) 'Unable to find initial .rsf file.  Tried:'
      Write (0,900) '   ' // Trim(FILEIN1)
      Write (0,900) '   ' // Trim(FILEIN2)
      Write (0,900) '   ' // Trim(FILEIN3)
      Write (0,900) '   ' // Trim(FILEIN4)  ;  Stop 803
  805 Write (0,900) 'Unable to locate .aij file: ' // Trim(AIJPATH)        
      Write (0,900) 'Unable to locate .acc file: ' // Trim(ACCPATH)        ;  Stop 805        
  806 Write (0,900) 'Error opening .aij file: '    // Trim(AIJPATH)        ;  Stop 806
  807 Write (0,900) 'Error opening atmosphere Z file: ' // Trim(ZFILEA)    ;  Stop 807
  808 Write (0,900) 'Error opening ocean Z file: ' // Trim(ZFILEO)         ;  Stop 808
  810 Write (0,900) 'Error opening initial .rsf file: ' // Trim(RSFPATH0)  ;  Stop 810
  880 Write (0,900) 'Error opening final .rsf file: ' // Trim(RSFPATHM)    ;  Stop 880
!****
  900 Format (A)
  902 Format (A4,'/',I2.2,'/01')
  903 Format (I4.4)
  908 Format ('Fractions: I,J,OC,GR,GI,LK,TO=',2I4,5F12.6)
  910 Format (/ A,'      SALT WATER MASS (kg/global area)' // &
      '                Atmos    Ground   LakeIce   LiqLake   LandIce   IceBerg    SeaIce      LiqOcean         Total' / &
      '                =====    ======   =======   =======   =======   =======    ======      ========         =====')
  911 Format (A10,F11.4,6F10.4,2F14.4)
  930 Format (A10,F11.4,74X,F14.4)
  931 Format (A10,F11.4,4F10.4,10X,F10.4,2F14.4)
  932 Format (A10,11X,4F10.4,10X,F10.4,2F14.4)
  935 Format (A10,11X,F10.4,10X,F10.4,44X,F14.4)
  936 Format (A10,11X,F10.4,64X,F14.4)
  938 Format (A10,31X,F10.4,30X,2F14.4)
  939 Format (A10,21X,2F10.4,20X,F10.4,2F14.4)
  944 Format (A10,11X,F10.4,10X,F10.4,44X,F14.4)
  945 Format (A10,11X,2F10.4,10X,2F10.4,24X,F14.4)
  961 Format (A10,51X,F10.4,10X,2F14.4)
  999 End



      Subroutine WATERRES (IDRSF,aXYP,oXYP,AREAG, RES)
!****
!**** Input:    IDRSF = ID of NetCDF input file
!****      aXYP,oXYP = grid cell area
!****          AREAG = global area
!**** Output: RES(9) = water (kg/global area) of 9 reservoirs
!****
      Implicit  None
      Integer,Parameter :: IMA=144,JMA= 90, IMO=288,JMO=180
      Integer,External  :: NF_INQ_VARID, NF_GET_VAR_DOUBLE
      Integer :: IDRSF
      Real*8  :: aXYP(IMA,JMA),oXYP(IMO,JMO),AREAG, RES(9)
!****
      Integer :: J,IMAX, IOSTAT,IDQ
      Real*8  :: WATMO(IMA,JMA),WGRND(IMA,JMA),WLAKI(IMA,JMA),WLIQL(IMA,JMA), &
                 WLANI(IMA,JMA),WICEB(IMA,JMA),WSEAI(IMA,JMA),WLIQO(IMO,JMO)
!****
      RES(:) = 0
!**** Calculate atmospheric water vapor and cloud water: RES(1)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'watmo', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WATMO)
      Do 10 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   10 RES(1) = RES(1) + Sum(WATMO(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate ground water: RES(2)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wgrnd', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WGRND)
      Do 20 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   20 RES(2) = RES(2) + Sum(WGRND(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate lake ice mass: RES(3)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wlaki', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WLAKI)
      Do 30 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   30 RES(3) = RES(3) + Sum(WLAKI(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate liquid lake mass: RES(4)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wliql', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WLIQL)
      Do 40 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   40 RES(4) = RES(4) + Sum(WLIQL(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate land ice: RES(5)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wlani', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WLANI)
      Do 50 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   50 RES(5) = RES(5) + Sum(WLANI(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate ice berg mass: RES(6)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wiceb', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WICEB)
      Do 60 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   60 RES(6) = RES(6) + Sum(WICEB(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate sea ice mass: RES(7)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wseai', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WSEAI)
      Do 70 J=1,JMA
      IMAX=IMA  ;  If(J==1.or.J==JMA) IMAX=1
   70 RES(7) = RES(7) + Sum(WSEAI(1:IMAX,J)*aXYP(1:IMAX,J))
!**** Calculate liquid ocean mass: RES(8)
      IOSTAT = NF_INQ_VARID      (IDRSF, 'wliqo', IDQ)
      IOSTAT = NF_GET_VAR_DOUBLE (IDRSF, IDQ, WLIQO)
      Do 80 J=1,JMO
      IMAX=IMO  ;  If(J==1.or.J==JMO) IMAX=1
   80 RES(8) = RES(8) + Sum(WLIQO(1:IMAX,J)*oXYP(1:IMAX,J))
!**** Divide mass by global area producing kg/m^2
      RES(1:8) = RES(1:8) / AREAG
      RES(9) = Sum (RES(1:8))
      Return
      EndSubroutine WATERRES



      Subroutine SGEOM (IM,JM,RADIUS, DXYP,AREAG)
!****
!**** SGEOM calculates the Spherical GEOMetry for a regular Lat-Lon grid
!****
      Implicit None
      Real*8,Parameter :: TWOPI = 6.283185307179586477d0
      Integer :: IM,JM, J
      Real*8  :: RADIUS, DXYP(*),AREAG, DLON,DLAT,FJEQ,VLATN,VLATS,VSINN,VSINS
!****
      If (JM /= 90 .and. JM /= 180)  GoTo 800
      DLON = TWOPI / IM
      DLAT = TWOPI*.5 / JM
      FJEQ = .5*(1+JM)
      Do 10 J=1,JM
      VLATN = DLAT*(J+.5-FJEQ)  ;  If(J==JM) VLATN =  TWOPI/4
      VLATS = DLAT*(J-.5-FJEQ)  ;  If(J==1)  VLATS = -TWOPI/4
      VSINN = Sin (VLATN)
      VSINS = Sin (VLATS)
   10 DXYP(J) = RADIUS*RADIUS*DLON*(VSINN-VSINS)
      AREAG = 2*TWOPI*RADIUS*RADIUS
      Return
!****
  800 Write (0,*) 'SGEOM: DXYP may not be correct for give JM =',JM
      Stop 800
      EndSubroutine SGEOM
