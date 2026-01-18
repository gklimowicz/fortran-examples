!@sum RMS finds RMS differences from observations and stores results
!@auth Gary Russell/Gavin Schmidt/Reto Ruedy
!@ver 1.0
C**** This program works with a template, the text file I.RMS .
C**** On each line, 8 words and a title have to be listed in the
C**** following order (width of field irrelevant)
C**** word  1        2      3  4  5      6   7      8
C****   Field  ObsSrce Domain M1 M2 Weight Fac Offset    Title
C**** They have to be chosen such that
C****   FieldIMXJM.ObsSrce = name of the obs.file for that line in OBS/
C****   Domain             = region over which observations are usable
C****   M1 and M2          = record numbers corr. to Jan/DJF and Jul/JJA
C****   Weight             = weighting to get a summary score
C****   Fac*obs.data+Off   = converts obs->model data units
C****   Title              = title in the model's lat-lon_diag-files
C****                        (not needed if same as previous line)
C****
C**** Command line arguments: [ IMxJM ] RunID  Jan_file [ Jul_file ]
C****                                          model_output_files
C****               defaults: IMxJM=72x46  Jul_file=Jan_file w/ JAN->JUL
C****
C**** Link (-s) the appropriate topography file to ./TOPO
C**** Link (-s) $CMRUN/OBS/RMSx to ./OBS ; cp $CMRUN/OBS/I.RMSx I.RMS
C****
C**** Now includes the Arcsin Mielke skill score
      IMPLICIT NONE
!@var LOBS number of RMS calculations    limited to LMAX
!@var NDIAG number of target diagnostics limited to LMAX
      INTEGER, PARAMETER :: LMAX=99

!**** Data read from I.RMS file are stored in OBS and TITLE
      CHARACTER OBS(LMAX)*99, TITLE(0:LMAX)*45, TITLE1*80, TITLE2*80

!**** Command line arguments -> grid,  RUN,    MONFILE
      character      line*99,  grid*7, RUN*20, MONFILE(2)*80, date*8
      CHARACTER*7, dimension(2) ::     MONTH =(/'January','   July'/)
      CHARACTER*6 :: BLANK ='      '

      REAL*4 TOTAL(20,2),TOTALA(20,2)
      REAL*4,  ALLOCATABLE :: DXYP(:),weights(:,:),DATA1(:,:),DATA2(:,:)
      REAL*4,  ALLOCATABLE ::          fearth(:,:),DOBS1(:,:),DOBS2(:,:)
      integer, ALLOCATABLE :: imaxj(:)

      CHARACTER*140, DIMENSION(0:LMAX,2) :: RMS, AMS
      DATA RMS /' ',LMAX*' ', ' ',LMAX*' '/
      DATA AMS /' ',LMAX*' ', ' ',LMAX*' '/
!@var NREC array of positions of targets in diagnostic file
      INTEGER, DIMENSION(LMAX) :: NREC
C****
      INTEGER IARGC,NARGS,M,N,L,NRUN,IREC,K,KM,KMIN,KMAX,NMOD
      INTEGER im,jm,ndiag,lobs,Lb,Lb1,len,len1,lname,narg1
      REAL*4 WRMS,RMSDIF,RRMS,XRMS,AMSCORE,XAMS
C****
      NARGS = IARGC()
      IF(NARGS.lt.2)  GO TO 800

!**** Determine the grid size and get runID
      grid='72X46'        ; narg1=2
      call getarg(1,run)

      if (index('123456789',run(1:1)) > 0) then
        grid=run ; n=index(grid,'x')
        if (n>0) grid=grid(1:n-1)//'X'//grid(n+1:7)
        call getarg(2,run) ; narg1=3
        IF(NARGS < 3) go to 800
      end if

      n=index(grid,'X')  ; read(grid(  1:n-1),*) im
      nrun=len_trim(run) ; read(grid(n+1:  7),*) jm

      allocate (  data1(im,jm),   data2(im,jm), dxyp(jm), imaxj(jm) )
      allocate (  DOBS1(im,jm),   DOBS2(im,jm) )
      allocate ( fearth(im,jm), weights(im,jm) )

      CALL GEOM(im,jm,dxyp,imaxj)

!**** Open observation template sheet and collect entries in TITLE/OBS
      OPEN(1,FILE='I.RMS',STATUS="OLD",FORM="FORMATTED",err=960)
      read(1,*) ; read(1,*)      ! skip 2 header lines
      L=0 ; Ndiag=0 ; title(0)=' '
      do while (L < LMAX)
         read(1,'(a)',end=5) line
         if (line(1:1) == '!') cycle       ! ignore comment lines
         if (line(1:1) == ' ') then        ! ignore leading blank lines
           if (L > 0) OBS(L)(46:46)='Y'    ! mark end of group of lines
           cycle
         end if

         call LocWrd (9, line, 99, Lb,len)        ! find loc of title
         if ( len.lt.1 .or. line(Lb:Lb+44).ne.title(n) ) then
           Ndiag=Ndiag+1 ; TITLE(Ndiag)=line(Lb:Lb+44)    ! new title
         end if

         L=L+1 ;   OBS(L)=' '
         write(OBS(L)(19:21),'(i3)') Ndiag
         call LocWrd (1,line,99, Lb,len) ; OBS(L)( 1: 5)=line(Lb:Lb+4)
         Lb1=Lb ; len1=len
         call LocWrd (2,line,99, Lb,len) ; OBS(L)( 7:10)=line(Lb:Lb+3)
         OBS(L)(63:99)=     ! standardized file name of obs. data file ?
     *        line(Lb1:Lb1+len1-1)//trim(grid)//'.'//line(Lb:Lb+len-1)
         call LocWrd (3,line,99, Lb,len) ; OBS(L)(12:17)=line(Lb:Lb+5)
         call LocWrd (4,line,99, Lb,len) ; OBS(L)(41:42)=line(Lb:Lb+1)
         call LocWrd (5,line,99, Lb,len) ; OBS(L)(43:44)=line(Lb:Lb+1)
         call LocWrd (6,line,99, Lb,len) ; OBS(L)(47:49)=line(Lb:Lb+2)
         call LocWrd (7,line,99, Lb,len) ; OBS(L)(51:54)=line(Lb:Lb+3)
         call LocWrd (8,line,99, Lb,len) ; OBS(L)(56:61)=line(Lb:Lb+5)
      end do
    5 Lobs  = L

      close (1)

C**** Read in FEARTH
C****
      open(1,file="TOPO",status="old",form="unformatted",err=950)
      read(1) ; read(1) ; read(1) title1,FEARTH
      close (1)

C**** Use first file to find target diagnostics
C****
      CALL GETARG (narg1,MONFILE(1))      ! January file
      if (nargs > narg1) then
        CALL GETARG (narg1+1,MONFILE(2))  ! July file
      else
        MONFILE(2)='JUL'//MONFILE(1)(4:80)
      end if

      M=1     ! in case MONFILE(1) is not found and err=970 applies
      OPEN(3,FILE=MONFILE(1),STATUS="OLD",FORM="UNFORMATTED",err=970)
      DO N=1,NDIAG
        REWIND (3)
        IREC=1
 10     READ(3,END=20) TITLE1
        IF ( TITLE1(1:len_trim(TITLE(N))) .eq. trim(TITLE(N)) ) THEN
          NREC(N)=IREC
          CYCLE
        ELSE
          IREC=IREC+1
          GOTO 10
        END IF
 20     WRITE (0,*) ' Target diagnostic not found: ',N,TITLE(N)
        WRITE (0,*) ' RMS computation will not be complete.'
        NREC(N)=-1
      END DO
      CLOSE (3)

C**** Calculate new RMS/AMS statistics and fill into RMS/AMS array
C****
!     Here we allow for adding a column to an old score sheet, a
!     feature that is not (yet) implemented above: the new column goes
!     into positions KMIN:KMIN+6 except that "RUN" may need more room
!     (future columns may wipe out the last characters of "RUN")

      KMIN = max(   1, index(RMS(1,1),'       ') )
      KMAX = min( 140, KMIN + max(6,nrun) )
      DO M=1,2
        RMS(0,M)(KMAX-NRUN+1:KMAX) = RUN(1:NRUN)
        AMS(0,M)(KMAX-NRUN+1:KMAX) = RUN(1:NRUN)

        OPEN(3,FILE=MONFILE(M),STATUS="OLD",FORM="UNFORMATTED",err=970)

C**** loop over calculations
        DO L=1,Lobs
C**** find appropriate target
          READ (OBS(L)(19:21),*) NMOD    ! target index
          REWIND(3)
          IF (NREC(NMOD).gt.-1) THEN
            DO N=1,NREC(NMOD)-1
              READ(3)
            END DO
            IF(OBS(L)(1:3).eq.'UVS')  THEN
              READ(3) TITLE1,DATA1
              READ(3) TITLE2,DATA2
            ELSE
              READ(3) TITLE1,DATA1
            END IF

            CALL SCORE(RRMS,AMSCORE,M,OBS(L),DATA1,DATA2,
     *        DOBS1,DOBS2,weights,im,jm,dxyp,imaxj,fearth)
            if (abs(RRMS) < 10000.) then
              WRITE (RMS(L,M)(KMIN:KMIN+6),921) RRMS
            else
              RMS(L,M)(KMIN:KMIN+6) = "Infinit"
            end if
            WRITE (AMS(L,M)(KMIN:KMIN+6),921) AMSCORE
          ELSE
            RMS(L,M)(KMIN:KMIN+6) = "*******"
            AMS(L,M)(KMIN:KMIN+6) = "*******"
          END IF
        END DO
        CLOSE (3)
      END DO
C****
C**** Calculate weighted RMS (printed as final line)
C****
      DO M=1,2
        KM = (KMIN+6)/7
        DO K=1,KM
          TOTAL(K,M) = 0.
          TOTALA(K,M) = 0.
          N = 0
          DO L=1,Lobs
            IF (RMS(L,M)(K*7-6:K*7-6) .ne. "I" .and.
     *          RMS(L,M)(K*7-6:K*7-6) .ne. "*") THEN
              READ (RMS(L,M)(K*7-6:K*7),921) XRMS
              READ (AMS(L,M)(K*7-6:K*7),921) XAMS
              READ (OBS(L)(47:49),*) WRMS
              TOTAL(K,M) = TOTAL(K,M) + XRMS*WRMS
              TOTALA(K,M) = TOTALA(K,M) + XAMS
              N = N + 1
            END IF
          END DO
          TOTALA(K,M)= TOTALA(K,M)/REAL(N)
        END DO
      END DO
C****
C**** Write RMS file to database
C****
      OPEN (3,FILE='RMS_'//trim(RUN),STATUS="UNKNOWN"
     *     ,FORM="FORMATTED")
      OPEN (4,FILE='AMS_'//trim(RUN),STATUS="UNKNOWN"
     *     ,FORM="FORMATTED")
      DO M=1,2
        CALL DATE_AND_TIME (DATE)
        WRITE (3,940) 'RMS    RMS differences for various runs    ' //
     *       DATE(1:4) // '/' // DATE(5:6) // '/' // DATE(7:8)
        WRITE (4,940) 'AMS  Arcsin Mielke scores for various runs ' //
     *       DATE(1:4) // '/' // DATE(5:6) // '/' // DATE(7:8)
        WRITE (3,940)
        WRITE (3,940) MONTH(M)
        WRITE (3,940)
        WRITE (3,901) 'Diag  Obsr Domain',RMS(0,M)(1:KMAX)
        WRITE (3,941) ('   ----',K=1,(KMIN+6)/7)
        WRITE (4,940)
        WRITE (4,940) MONTH(M)
        WRITE (4,940)
        WRITE (4,901) 'Diag  Obsr Domain',RMS(0,M)(1:KMAX)
        WRITE (4,941) ('   ----',K=1,(KMIN+6)/7)
        DO L=1,Lobs
          WRITE (3,901) OBS(L)(1:18),RMS(L,M)(1:KMIN+6)
          IF(OBS(L)(46:46).eq.'Y')  WRITE (3,901)
          WRITE (4,901) OBS(L)(1:18),AMS(L,M)(1:KMIN+6)
          IF(OBS(L)(46:46).eq.'Y')  WRITE (4,901)
        END DO
        WRITE (3,942) ( '  -----',K=1,(KMIN+6)/7)
        WRITE (3,943) (NINT(TOTAL(K,M)),K=1,(KMIN+6)/7)
        IF(M.eq.1)  WRITE (3,*)
        WRITE (4,942) ( '  -----',K=1,(KMIN+6)/7)
        WRITE (4,945) (TOTALA(K,M),K=1,(KMIN+6)/7)
        IF(M.eq.1)  WRITE (4,*)
      END DO
      CLOSE (3)
      CLOSE (4)
      STOP
C****
  800 WRITE (0,*)
     *' Usage: RMS [IMxJM] run jan_diag_file [jul_diag_file] 2006/10/09'
      WRITE (0,*)
     *'        Calculates RMS and Arcsin Mielke statistics for JAN/JUL'
      WRITE (0,*)
      WRITE (0,*) ' Example:  RMS 72x46 E001 JAN1950-1955.ijE001 '
      WRITE (0,*)
      WRITE (0,*) ' Do: ln -s appropriate_topo_file TOPO'
      WRITE (0,*) '     ln -s base/OBS/I.RMS? I.RMS'
      WRITE (0,*) '     ln -s base/OBS/RMS? OBS '
      WRITE (0,*)
     *' in your current directory before executing the command RMS '
      WRITE (0,*)
     *' All data file records are of the form: title*80,a(im,jm),..'
      WRITE (0,*)
      GO TO 999
C****
  901 FORMAT (A17,1X,A)
  921 FORMAT (F7.2)
  940 FORMAT (A)
  941 FORMAT ('----  ---- ------',1X,20A7)
  942 FORMAT (18X,20A7)
  943 FORMAT (18X,20I7)
  945 FORMAT (18X,20F7.2)
  950 WRITE (0,*) ' Error opening topography file TOPO'
      STOP 950
  960 WRITE (0,*) ' Error opening I.RMS'
      STOP 960
  970 WRITE (0,*) ' Error opening model diag file ',MONFILE(M)
      STOP 970
  999 END

      SUBROUTINE SCORE(RMSDIF,AMSCORE,M,OBS,UM,VM,
     *   uo,vo,weights,im,jm,dxyp,imaxj,fearth)
!@sum SCORE calculates the RMS difference and Arcsin Mielke score
!@+   between a model data record and an observation data record
!@auth Gary Russell/Gavin Schmidt
C**** Input: M = 1 for January, 2 for July
C****      UM,VM data arrays
      IMPLICIT NONE
      integer im,jm,imaxj(jm)
      real*4 dxyp(jm)
      real*4 :: TWOPI = 6.2831853
      REAL*4, PARAMETER :: SKIPOBS=-999999., SKIP=-1.e30
      CHARACTER OBS*99, FILEIN*80, TITLE*80
      REAL*4, DIMENSION(IM,JM) :: UM,VM, UO,VO, weights, fearth
      REAL*4 FAC,OFFSET,DLAT
      REAL*4 W,Q,M2,O2,MBAR,OBAR,MDIFF,VARM,VARO,MSE,RMSDIF,AMSCORE
      INTEGER N,NOBS,I,J,J60,M
C**** factors and offsets to transform model output to commensurate
C**** units
      READ(OBS(56:61),*) OFFSET
      READ(OBS(51:54),*) FAC
C****
C**** Open model and observation data files
C****
      FILEIN = 'OBS/'//OBS(63:99)
      OPEN (2,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD',ERR=801)
C****
C**** Skip through unused records
C****
      READ (OBS(39+2*M:40+2*M),*) NOBS
      DO N=1,NOBS-1
        READ (2)
      END DO
C****
C**** Determine weights
C****
      IF(OBS(12:17).eq.'NHGrnd' .or. OBS( 1: 3).eq.'UVS')  THEN
        DO I=1,IM
          WEIGHTS(I,JM/2+1:JM) = DXYP(JM/2+1:JM)*FEARTH(I,JM/2+1:JM)
          WEIGHTS(I,1:JM/2) = 0.
        END DO
      END IF
      IF(OBS(12:17).eq.'Glob60')  THEN
        DLAT = .5*NINT(360./JM)
        J60  = 1 + NINT(JM/2 - 60./DLAT)
        WEIGHTS(:,1:J60-1) = 0.
        WEIGHTS(:,2+JM-J60:JM) = 0.
        DO I=1,IM
          WEIGHTS(I,J60:1+JM-J60) = DXYP(J60:1+JM-J60)
        END DO
      END IF
      IF(OBS(12:17).eq.'Global')  THEN
        DO I=1,IM
          WEIGHTS(I,:)=DXYP(:)
        END DO
      END IF
C****
C**** Calculate global scores
C****
      IF (OBS( 1: 3).ne.'UVS') THEN
        READ (2,ERR=810) TITLE,UO
        W = 0. ; Q = 0.
        MBAR = 0. ; OBAR = 0.
        M2 = 0. ; O2 = 0.
        DO J=1,JM
          DO I=1,IMAXJ(J)
            IF(UM(I,J).le.SKIP .or. UO(I,J).eq.SKIPOBS .or.
     *           WEIGHTS(I,J).eq.0.) CYCLE
            W = W + WEIGHTS(I,J)
            Q = Q + WEIGHTS(I,J)*(UM(I,J)*FAC+OFFSET-UO(I,J))**2
            MBAR = MBAR + WEIGHTS(I,J)*(UM(I,J)*FAC+OFFSET)
            OBAR = OBAR + WEIGHTS(I,J)* UO(I,J)
            M2 = M2 + WEIGHTS(I,J)*(UM(I,J)*FAC+OFFSET)**2
            O2 = O2 + WEIGHTS(I,J)*UO(I,J)**2
          END DO
        END DO
      ELSE
C****
C**** Calculate vector RMS
C****
        READ (2,ERR=810) TITLE,UO,VO
        W = 0. ; Q = 0.
        MBAR = 0. ; OBAR = 0.
        M2 = 0. ; O2 = 0.
        DO J=1,JM
          DO I=1,IM
            IF(UM(I,J).le.SKIP .or. UO(I,J).eq.SKIPOBS .or.
     *           VM(I,J).le.SKIP .or. VO(I,J).eq.SKIPOBS .or.
     *           WEIGHTS(I,J).eq.0) CYCLE
            W = W + WEIGHTS(I,J)
            Q = Q + WEIGHTS(I,J)*((UM(I,J)*FAC+OFFSET-UO(I,J))**2
     *           +(VM(I,J)*FAC+OFFSET-VO(I,J))**2)
            MBAR = MBAR + WEIGHTS(I,J)*((UM(I,J)+VM(I,J))*FAC+OFFSET)
            OBAR = OBAR + WEIGHTS(I,J)*(UO(I,J)+VO(I,J))
            M2 = M2 + WEIGHTS(I,J)*((UM(I,J)*FAC+OFFSET)**2+
     *           (VM(I,J)*FAC+OFFSET)**2)
            O2 = O2 + WEIGHTS(I,J)*(UO(I,J)**2 + VO(I,J)**2)
          END DO
        END DO
      END IF
C****
C**** Final calculation of RMS and AMSCORE and close datafiles
C****
      MSE = Q/(W+1.e-20)
      RMSDIF = SQRT(MSE)
      VARM = M2/(W+1.e-20) - (MBAR/(W+1.e-20))**2
      VARO = O2/(W+1.e-20) - (OBAR/(W+1.e-20))**2
      MDIFF = ((MBAR-OBAR)/(W+1.e-20))**2
      AMSCORE = (4./TWOPI)*ASIN(1.- MSE/(VARO + VARM + MDIFF))
      CLOSE (2)
      RETURN
C****
  801 WRITE (0,*) ' Error opening datafile: ',FILEIN(1:40)
      STOP 801
  810 WRITE (0,*) ' Error reading datafile: ',OBS
      STOP 810
      END subroutine score

      SUBROUTINE GEOM (im,jm,dxyp,imaxj)
      implicit none
      integer j,im,jm,imaxj(jm)
      real*4 DXYP(JM)
      real*4 TWOPI,DLAT,FJEQ,SINS,SINN
C****
C**** Calculates spherical geometry
C****
      TWOPI = 6.2831853
      DLAT  = TWOPI*NINT(360./(JM-1))/720.
      FJEQ  = (1+JM)/2.
C**** Geometric parameters centered at primary latitudes
      DO J=2,JM-1
        SINS    = SIN(DLAT*(J-.5-FJEQ))
        SINN    = SIN(DLAT*(J+.5-FJEQ))
        DXYP(J) = (SINN-SINS)
        imaxj(j)=IM
      end do
      DXYP(JM)= (1.-SINN)
      DXYP(1) = DXYP(JM)
      imaxj(1) = 1 ; imaxj(jm) = 1

      RETURN
      END subroutine geom

      subroutine LocWrd (n,str,lenx,Lbeg,len)
!@sum LocWrd finds location of n-th word in string: str(Lbeg:Lbeg+len-1)
!@auth Reto Ruedy
      implicit none
      integer         , intent(in) :: n,lenx
      character*(lenx), intent(in) :: str
      integer         , intent(out) :: Lbeg,len
      integer :: nw

      len=0  ! flag if string has less than n words

      Lbeg=1
      do while (str(Lbeg:Lbeg) == ' ')             ! skip leading blanks
         Lbeg=Lbeg+1 ; if (Lbeg == lenx) return
      end do

      nw=1
      do while (nw < n)                    ! find beginning of next word
         Lbeg=Lbeg+1 ; if (Lbeg == lenx) return
         if (str(Lbeg:Lbeg) == ' ') then          ! found end of word nw
            Lbeg=Lbeg+1 ; if (Lbeg == lenx) return
            do while (str(Lbeg:Lbeg) == ' ')
               Lbeg=Lbeg+1 ; if (Lbeg == lenx) return
            end do
            nw=nw+1
         end if
      end do

      len=1                               ! find end of word n
      do while (str(Lbeg+len:Lbeg+len) .ne. ' ')
         if (Lbeg+len-1 == lenx) return
         len=len+1
      end do

      return
      end subroutine LocWrd
