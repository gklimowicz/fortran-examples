      PROGRAM AVR_modelE
C**** AVR calculations 144x90 as well
C**** AVR72X46.S     Calculate AVRX coefficients     97/11/20
C****
C**** Compile: FCG AVR72X46.S
C****
C**** AVR72X46 calculates the AVRX coefficients for ocean polar filter.
C**** The size of the output file depends on NSEGM and INDM.  Start
C**** with large values.  The .PRT file will show what values are
C**** needed.  Rerun this program with those values for NSEGM and INDM
C**** in order to minimize the size of the output file.  Also insert
C**** those values in the PARAMETER statement of subroutine OPFIL.
C****
C**** LMO   ZE1   ZERAT   Topog     NSEGM   INDM   File Made
C**** ---   ---   -----   -----     -----   ----   ---------
C****  13    12     1.5   Z4X510S     7    92700   AVR4X5LD
C****  13    12     1.5   Z4X511S     7    90372   AVR4X5LD.Z11
C****  13    12     1.5   Z4X512S     7    90365   AVR4X5LD.Z12
C****  13    12     1.5   Z72X46N     7    90365   AVR72X46.L13
C****  13    12     1.5 Z144X90N_nocasp 10 941701  AVR144X90_nocasp.L13.modelE
C****  13    12     1.5 Z72X46N_9KY  6    90066  AVR72X46.L13.9KY.modelE
C**** 13    12     1.5 Z72X46N_9KY_nobs  6    90066  AVR72X46.L13a.9KY_nobs.modelE lmo_min=1
C****  13    12     1.5 Z72X46N_9KY_nobs  6    90066  AVR72X46.L13b.9KY_nobs.modelE lmo_min=2
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IM=72,JM=46,LMO=13,J1O=4,JFMAX=13,NSEGM=15,INDM=168130
     *     ,JXMAX=2*JFMAX-1,JNOF=JM-2*JFMAX,ZE1=12.,ZERAT=1.5, TWOPI = 6
     *     .283185307179586477d0)
C**** JFMAX = last latitude at which polar filter is applied (40 deg)
C**** JXMAX = maximum J index when NH filtered J's are added to SH
C**** JNOF  = number of interior latitudes that are not filtered
      INTEGER*2 NSEG(LMO,J1O:JXMAX),IMIN(NSEGM,LMO,J1O:JXMAX),
     *    ILEN(NSEGM,LMO,J1O:JXMAX),IMAX(NSEGM,LMO,J1O:JXMAX)
      INTEGER*4 LMM(IM,JM),INDEX(IM,2:JFMAX),LENNOF(2:JFMAX+1),iargc
      REAL*4 FOCEAN(IM,JM),ZATMO(IM,JM),HOCEAN(IM,JM), REDUCO(INDM)
      REAL*8 ZE(0:LMO), REDUC(IM,IM)
      CHARACTER*1 QSEG(IM,2:JFMAX)
      CHARACTER*80 TITLE
      character*100 filein,fileout
      integer ier
      integer*2, allocatable, dimension(:,:,:) :: imin_A,ilen_A
      integer*4, allocatable, dimension(:) :: reduco_A
      COMMON /GEOMCB/ DXYP(JM),DXP(JM),DYP(JM),DXV(JM)
C****
C**** Calculate maximum segment length below which polar filter
C**** is not applied:  DXP(J)*(LEN+1) > DYP(3)*LEN
C****
      if((iargc().eq.1).or.(iargc().eq.2)) then
        CALL GETARG(1,filein)
        fileout="AVR"//trim(filein)
        if(iargc().eq.2) CALL GETARG(2,fileout)
      else
        print*,"AVR72X46 calculates the AVRX coefficients for "
        print*,"ocean polar filter."
        print*,"AVR_modelE <topo filein> opt:<avr fileout>"
        print*,"Contains info to do AVR144X92 but req's recompiling"
        go to 999
      endif

      QSEG(:,:)=' '
      CALL SGEOM
      DO 10 J=2,JFMAX+1
   10 LENNOF(J) = DXP(J)/(DYP(3)-DXP(J))
      WRITE (6,901) (J,J=2,JFMAX+1),LENNOF
C****
C**** Calculate LMM
C****
C**** Calculate ZE
      DO 110 L=0,LMO
  110 ZE(L) = ZE1*(ZERAT**L-1.)/(ZERAT-1.)
C**** Read in FOCEAN, ZATMO and HOCEAN
c      FILEIN = 'Z72X46N_gas.1_nocasp_21k_ice5g'
c      FILEIN = 'Z72X46N_9KY_nobs'
c      FILEIN = '/u/cmrun/Z144X90N_nocasp.1'
      OPEN  (11,FILE=FILEIN,FORM='UNFORMATTED',STATUS='OLD')
      READ  (11) TITLE,FOCEAN
      WRITE (6,*) 'FOCEAN read: ',TITLE
      READ  (11)
      READ  (11)
      READ  (11)
      READ  (11) TITLE,ZATMO
      WRITE (6,*) 'ZATMO read: ',TITLE
      READ  (11) TITLE,HOCEAN
      WRITE (6,*) 'HOCEAN read: ',TITLE
      CLOSE (11)
C**** Calculate LMM
      LMO_min= 2  ! can vary!
      DO 170 I=1,IM*JM
      L=0
      IF(FOCEAN(I,1).le.0.)  GO TO 160
      DO 150 L=LMO_min,LMO-1
  150 IF(HOCEAN(I,1) .le. .5*(ZE(L)+ZE(L+1)))  GO TO 160
C     L=LMO
  160 LMM(I,1)=L
  170 continue
C****
C**** Determine contiguous interfaces of ocean boxes at each layer
C**** and latitude.
C**** NSEG(L,JX)   = number of groups of contiguous interfaces
C**** IMIN(N,L,JX) = zeroeth I of the N-th group of interfaces
C**** IMAX(N,L,JX) = last I of the N-th group of interfaces
C****
      NSEGN=0
      DO 300 JX=J1O,JXMAX
      J =JX
      JA=JX
      IF(JX.gt.JFMAX)  J =JX+JNOF
      IF(JX.gt.JFMAX)  JA=2*JFMAX+1-JX
      DO 300 L=1,LMO
      N=0
C**** Locate the first longitude whose box is land
      IL1=0
      IF(LMM(IM,J).lt.L)  GO TO 220
      DO 210 IL1=1,IM-1
      IF(LMM(IL1,J).lt.L)  GO TO 220
  210 CONTINUE
C**** The grid boxes are ocean for all longitudes
      N=1
      IMIN(N,L,JX) = 0
      IMAX(N,L,JX) = IM
      ILEN(N,L,JX) = IM
      GO TO 290
C**** IL1 is a longitude whose box is land, find the first
C**** subsequent longitude whose box is ocean
  220 DO 230 IW1=IL1+1,IM
      IF(LMM(IW1,J).ge.L)  GO TO 240
  230 CONTINUE
C**** There are no more ocean boxes
      GO TO 290
C**** IW1 is the first longitude whose box is ocean, find the last
C**** longitude whose box is ocean
  240 DO 250 IL1=IW1+1,IM
      IF(LMM(IL1,J).lt.L)  GO TO 270
  250 CONTINUE
C**** The last longitude IM is an ocean box, find the first
C**** longitude starting from 0 whose box is land
      DO 260 IL1=1,IM
      IF(LMM(IL1,J).lt.L)  GO TO 280
  260 CONTINUE
      STOP 260
C**** IL1-1 is the last longitude whose box is ocean
  270 LENGTH = IL1-IW1-1
      IF(LENGTH.le.0)  GO TO 220
      IF(LENGTH.le.LENNOF(JA))  QSEG(LENGTH,JA) = 'N'
      IF(LENGTH.le.LENNOF(JA))  GO TO 220
      N=N+1
      IMIN(N,L,JX) = IW1-1
      IMAX(N,L,JX) = IL1-2
      ILEN(N,L,JX) = LENGTH
      QSEG(LENGTH,JA) = 'R'
      GO TO 220
C**** IL1-1 is the last longitude whose box is ocean, but this is
C**** the last ocean segment because IL1 crossed IM
  280 LENGTH = IM+IL1-IW1-1
      IF(LENGTH.le.0)  GO TO 290
      IF(LENGTH.le.LENNOF(JA))  QSEG(LENGTH,JA) = 'N'
      IF(LENGTH.le.LENNOF(JA))  GO TO 290
      N=N+1
      IMIN(N,L,JX) = IW1-1
      IMAX(N,L,JX) = IM+IL1-2
      ILEN(N,L,JX) = LENGTH
      QSEG(LENGTH,JA) = 'R'
  290 NSEG(L,JX)   = N
      WRITE (6,929) J,L,N,(IMIN(K,L,JX),IMAX(K,L,JX),K=1,N)
      IF(N.gt.NSEGN)  NSEGN=N
  300 CONTINUE
      WRITE (6,*) 'NSEGN,NSEGM =',NSEGN,NSEGM
      WRITE (0,*) 'From AVR: NSEGN,NSEGM =',NSEGN,NSEGM
      WRITE (6,*) (I,(QSEG(I,JA),JA=2,JFMAX),I=1,IM)
C****
C**** Loop over absolute latitudes and segment lengths.  Determine
C**** starting location in the large reduction matrix, calculate
C**** entries in the local reduction matrix, and load them into
C**** the large matrix.
C**** INDEX(ILEN,JA) = starting location of large reduction matrix
C**** REDUCO(IND)   = large reduction matrix
C****
      IND=0
      DO 480 JA=2,JFMAX
      DO 480 IL=1,IM-1
      IF(QSEG(IL,JA).ne.'R')  GO TO 480
C****
C**** Calculate the local reduction matrix for given IL and JA.
C**** A(I) = A(I) - sum REDUC(K,I)*A(K)  for K = 1,IL
C****
      KM = 2*(IL+1)
C**** Zero out local reduction factor matrix
      DO 410 I=1,IL
      DO 410 K=1,IL
  410 REDUC(K,I) = 0.
C**** Loop over wave numbers from IL down to 1.  Calculate the
C**** reduction for each wave number and the grid point
C**** contribution of it from all other grid points.
      DO 440 N=IL,1,-1
      REDUCN = 1. - DXP(JA)*(IL+1)/(DYP(3)*N)
      IF(REDUCN.le.0.)  GO TO 440
      DO 420 I=1,IL
      DO 420 K=1,IL
  420 REDUC(K,I) = REDUC(K,I) +
     +  SIN(TWOPI*N*K/KM)*(SIN(TWOPI*N*I/KM)*REDUCN*4./KM)
C     DO 430 I=(IL+3)/2,IL
C     DO 430 K=1,IL
C 430 REDUC(K,I) = REDUC(IL+1-K,IL+1-I)
  440 CONTINUE
      IF(IL.gt.5)  GO TO 460
      WRITE (6,944) JA,IL
      DO 450 I=1,IL
  450 WRITE (6,945) (REDUC(K,I),K=1,IL)
  460 CONTINUE
C**** Load the reduction contribution matrix into the large
C**** coefficient array
      INDEX(IL,JA)=IND
      DO 470 I=1,IL
      DO 470 K=1,IL
      IND=IND+1
  470 REDUCO(IND) = REDUC(K,I)
  480 continue
      WRITE (6,*) 'IND,INDM =',IND,INDM
      WRITE (0,*) 'From AVR: IND,INDM =',IND,INDM
C****
C**** Write locating parameters and reduction matrix to disk
C****
c      OPEN  (2,FILE=fileout,FORM='UNFORMATTED')
c      write(title,*) "AVR parameters: NSEGM=",NSEGM,",INDM=",INDM
c      write(2) title,NSEGM,INDM      
c      TITLE = 'Polar Filter Reduction Matrix ' //
c     *        'Z72X46N_9KY_nobs, 13L, Lmin=2 '
c      WRITE (2) TITLE,NSEG,IMIN,ILEN,INDEX,REDUCO
c      CLOSE (2)

      allocate(imin_A(nsegn,lmo,j1o:jxmax),stat=ier)
      allocate(ilen_A(nsegn,lmo,j1o:jxmax),stat=ier)
      allocate(reduco_A(ind),stat=ier)
 
      OPEN  (2,FILE=fileout,FORM='UNFORMATTED')
      write(title,*) "AVR parameters: NSEGn=",NSEGn,",IND=",IND
      write(2) title,NSEGn,IND
      TITLE = 'Polar Filter Reduction Matrix ' //
     *        trim(filein)// ', 13L, Lmin=2 '
      WRITE (2) TITLE,NSEG,IMIN_A,ILEN_A,INDEX,REDUCO_A
      CLOSE (2)
     
      WRITE (6,*) 'Output written on unit 2:  ',trim(fileout)
C****
  901 FORMAT ('0',9X,13I4/' LENNOF = ',13I4/)
  929 FORMAT (' J,L,NSEG=',3I3,12(I5,I3))
  930 FORMAT ('0Existence of segment'/6X,'23456789ABCD'/
     *  (I3,'  X',12A1,'X'))
  944 FORMAT ('0J,IL=',2I3,5X,'Reduction contribution matrix')
  945 FORMAT (1X,13F10.4)
 999  continue

      END program AVR_modelE

      SUBROUTINE SGEOM
C****
C**** SGEOM calculates the spherical geometry for the C grid
C****
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (IM=72,JM=46)
!      PARAMETER (IM=144, JM=90)
      COMMON /GEOMCB/ DXYP(JM),DXP(JM),DYP(JM),DXV(JM)
C****
C**** Calculates the spherical geometry for the C grid
C****
      TWOPI = 6.283185307179586477D0
      DLON  = TWOPI/IM
!      DLAT  = .5*TWOPI/REAL(JM)  ! full polar box
      DLAT   = .5*TWOPI/(JM-1)    ! 1/2 polar box
      RADIUS= 6375000.
C**** Geometric parameters centered at secondary latitudes
      FJEQ = .5*(1+JM)
      DO 110 J=1,JM-1
      SINS = SIN(DLAT*(J   -FJEQ))
      SINN = SIN(DLAT*(J+1.-FJEQ))
  110 DXV(J)  = DLON*RADIUS*(SINN-SINS)/DLAT
C**** Geometric parameters centered at primary latitudes
      DO 130 J=2,JM-1
      SINS    = SIN(DLAT*(J-.5-FJEQ))
      SINN    = SIN(DLAT*(J+.5-FJEQ))
      DXYP(J) = DLON*RADIUS*RADIUS*(SINN-SINS)
      DXP(J)  = .5*(DXV(J-1)+DXV(J))
  130 DYP(J)  = DLAT*RADIUS
      DXYP(JM)= DLON*RADIUS*RADIUS*(1.-SINN)
      DXYP(1) = DXYP(JM)
      DXP(1)  = .5*DXV(1)
      DXP(JM) = .5*DXV(JM-1)
      DYP(1)  = .5*DLAT*RADIUS
      DYP(JM) = .5*DLAT*RADIUS
      RETURN
      END subroutine sgeom

