C
C Program to take model O3 model output and prepare files for ozone input
C to the radiation code, based upon Andy's routines: TROPO3_SHINDELL,
C REPART, etc. Expected input files are of the normal GISS IJ binary
C format, with records like:
C   1: Ox L  1 (10^-7 V/V air)               JAN 2051      E013TdsS15M23     
C   2: Ox L  2 (10^-7 V/V air)               JAN 2051      E013TdsS15M23 
C ...
C  23: Ox L 23 (10^-7 V/V air)               JAN 2051      E013TdsS15M23     
C  24: Ox L  1 (10^-7 V/V air)               FEB 2051      E013TdsS15M23     
C ...
C 276: Ox L 23 (10^-7 V/V air)               DEC 2051      E013TdsS15M23 
C
C User may need to change a few parameters, E.g.: 
C ltop -- (if not 23 layer model output). If you change this you should 
C         also change DATPLB to have the nominal pressure levels of your
C         GCM output.
C lgcmtop --  (the level among the 49 rad code levels above which user is
c         providing no input). See below.
C convert_power -- See below
C
C If you are not using the "medium" 4.x5. degree resolution for your
C model output, you should first put your data on that grid, because
C the radiation code uses that grid internally.
C
C Ox input to this program are assumed to be in 1E-7 pppv (i.e. V/V)
C But if you have the model output in say 1E-10 pppm (i.e. kg/kg) simply 
C change the parameter "convert_power" and uncomment the factor of
C *mwAIR/mwO3 in the definition of the CONVERT variable.
C Output of this program are in 1/1000 DU (atm*cm???).
C
C When the radiation levels get too high and you have no
C chemistry output, program pastes in values from the file:
C /raid4/jan2004_o3_shindelltrop_72x46x49x12_1990
C (that's /raid4 availbale on athena at GISS. Let me know if you 
C can't get the file gfaluvegi@giss.nasa.gov)
C This should be fine very high up, but be careful if you are
C not inputting enough Ox data. (and change parameter lgcmtop please).
C You could use a differnt file by changing O3_fn_sup.
C
C Program will prompt user for filename of O3 input, and for a 
C representative year of this data (just for labelling purposes).
C
C auth= g.faluvegi based on a.lacis. date=7.7.06
C
  
      IMPLICIT REAL*8(A-H,O-Z) 

      integer, parameter :: im=72, jm=46, ltop=23, nmon=12, lrad=49,
     & lgcmtop=48 ! <-- the highest rad level with avail. chem data
      real*8, parameter :: convert_power=1.d-7,scale_height=7.991D+05
     &,mwO3=48.D0,mwAIR=28.97D0
      real*8, parameter, dimension(lrad+1) :: PLBMOD = (
     &  /984.,934.,854.,720.,550.,390.,285.,210.,150.,125.,100.,
     &  80.,60.,55.,50.,45.,40.,35.,30.,25.,20.,15.,10.0,7.00,
     &  5.000,4.000,3.000,2.000,1.500,1.000,0.700,0.500,0.400,
     &  0.300,0.200,0.150,0.100,0.070,0.050,0.040,0.030,0.020,
     &  0.015,0.010,0.007,0.005,0.004,0.003,0.001,0.0000001/ )
      real*8, parameter, dimension(ltop+1) :: DATPLB =
     & (/984.,960.,929.,884.,819.,710.,570.,425.,314.,245.
     & ,192.,150.,117.,86.2,56.2,31.6,17.8,10.0, 4.6, 1.5
     & , 0.5,0.145,0.0312,0.00206/)
      real*4, dimension(im,jm) :: A, B, XY
      real*4, dimension(im,jm,ltop,nmon) :: O3ICMA
      real*4, dimension(im,jm,lrad,nmon) :: O3YEAR, paste
      real*8, dimension(lgcmtop) :: O3RADLEV
      real*8, dimension(ltop) :: DATO3
      character*80 :: title,title2,O3_fn_in,O3_fn_out,O3_fn_sup
      character*3, parameter, dimension(nmon+1) :: CHMON = (
     &/'JAN','FEB','MAR','APR','MAY','JUN'
     &,'JUL','AUG','SEP','OCT','NOV','DEC','ANN'/ )
      character*2 :: CHLAY
      character*4 :: CYR
      character*5 :: CHPLB5,CHPLT5
          
      O3_fn_sup='jan2004_o3_shindelltrop_72x46x49x12_1990'
      write(6,*)'Enter filename of tracer output'
      read(5,*) O3_fn_in
      O3_fn_out='o3_shindell_72x46x49x12_from_file_'//trim(o3_fn_in)
      write(6,*)'Enter year this file represents (YYYY)'
      read(5,*) CYR

      CONST=scale_height*convert_power/PLBMOD(1) !! *mwAIR/mwO3

C Read model Ox data and the file for pasting values above those data:
C Also convert units:

      OPEN(12,FILE=trim(O3_fn_sup),FORM='UNFORMATTED',STATUS='OLD') 
      OPEN(11,FILE=trim(O3_fn_in),FORM='UNFORMATTED',STATUS='OLD') 
      do M=1,nmon
        do L=1,ltop      
          FACTOR=(DATPLB(L)-DATPLB(L+1))*CONST
          READ(11) TITLE,A
          O3ICMA(:,:,L,M)=A(:,:)*FACTOR 
        enddo
        do L=1,lrad
          READ(12) TITLE,B
          paste(:,:,L,M)=B(:,:)
        enddo
      enddo 
      CLOSE(11)
      CLOSE(12)

C Call the radiation code subroutine for re-binning data from
C model to radiation levels. Paste above model levels:

      do M=1,nmon 
        do j=1,jm
          do i=1,im   
            DATO3(1:ltop)=O3ICMA(I,J,1:ltop,M) 
            CALL REPART(DATO3,DATPLB,ltop+1,O3RADLEV,
     &      PLBMOD(1:lgcmtop),lgcmtop)
            do L=1,lgcmtop-1
              O3YEAR(I,J,L,M)=O3RADLEV(L)
            enddo
            do L=lgcmtop,lrad
              O3YEAR(I,J,L,M)=paste(I,J,L,M)
            enddo
          enddo
        enddo
      enddo

       
C Write results to file format to be read by radiation code:

      OPEN(25,FILE=O3_fn_out,FORM='UNFORMATTED',STATUS='NEW')
      do m=1,nmon 
        title2(1:4)=CYR
        title2(5:5)=' '
        title2(6:8)=CHMON(m)
        title2(10:12)=' L='
        do L=1,lrad ! was lgcmtop   
          WRITE(CHLAY,'(I2)') L 
          IF(PLBMOD(L).GE.1.0) THEN
            WRITE(CHPLB5,'(F5.1)') PLBMOD(L)
            WRITE(CHPLT5,'(F5.1)') PLBMOD(L+1)  
          ELSE
            WRITE(CHPLB5,'(F5.3)') PLBMOD(L)
            WRITE(CHPLT5,'(F5.3)') PLBMOD(L+1)  
          ENDIF
          title2(13:14)=CHLAY
          title2(15:17)=' p='
          title2(18:22)=CHPLB5
          title2(23:23)='-' 
          title2(24:28)=CHPLT5
          title2(34:80)=' from file '//trim(o3_fn_in)
          do j=1,jm
            do i=1,im 
              XY(I,J)=O3YEAR(I,J,L,M)   
              IF(XY(I,J) < 0.) WRITE(6,66) JYEAR,M,L,J,I,XY(I,J)
   66         FORMAT(/5I5,F10.5)
            enddo 
          enddo   
          WRITE(25) title2,XY   
        enddo
      enddo
      CLOSE(25)
      
      END


      SUBROUTINE REPART(FXL,XLB,NXB,GYL,YLB,NYB)
      IMPLICIT NONE

C     ------------------------------------------------------------------
C
C     REPART   Repartitions FXL (a histogram-tyoe distribution function)
C              where XLB depicts the NXB partitions that define FXL data
C              FXL is assumed to be constant between XLB(N) AND XLB(N+1)
C
C              GYL(N) is the new histogram distribution function defined
C              by the (input) NYB partitions YLB(N).  The YLB partitions
C              can differ from XLB in number and/or spacing, or in upper
C              and lower limits.  XLB and YLB coordinates are assumed to
C              be linear both increasing or decreasing in the same sense
C
C     RETURNS  GYL as a histogram-type distribution function, with NYB-1
C              values assumed to be constant between YLB(N) and YLB(N+1)
C
C       NOTE:  The column amount of a vertically distributed quantity is
C              conserved, within the Repartition Interval that is common
C              to to XLB and YLB (layer bottom edge and top edge) limits
C
C     ------------------------------------------------------------------
      INTEGER, INTENT(IN) :: NXB,NYB
      REAL*8, INTENT(IN) :: FXL(NXB-1),XLB(NXB),YLB(NYB)
      REAL*8, INTENT(OUT) :: GYL(NYB-1)
      INTEGER NXF,NYG
      REAL*8 SUMG,SUMY,XA,YA,XB,YB,XAYA,PART
      INTEGER I,J

      NXF=NXB-1
      NYG=NYB-1
      SUMG=0.D0
      DO 50 I=1,NYG
      GYL(I)=0.D0
   50 CONTINUE
      SUMY=0.D0
      I=1
      XA=XLB(I)
      J=1
      YA=YLB(J)
      XB=XLB(I+1)
      IF(XB < XA) GO TO 200
  100 CONTINUE
      YB=YLB(J+1)
      IF(YB > XA) GO TO 110
      GYL(J)=0.D0
      J=J+1
      IF(J > NYG) GO TO 160
      YA=YB
      GO TO 100
  110 CONTINUE
      XB=XLB(I+1)
      IF(XB > YA) GO TO 120
      I=I+1
      IF(I > NXF) GO TO 160
      XA=XB
      GO TO 110
  120 CONTINUE
      XAYA=XA
      IF(YA > XA) XAYA=YA
      IF(YB > XB) GO TO 130
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG
      J=J+1
      IF(J > NYG) GO TO 160
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 120
  130 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I > NXF) GO TO 140
      XA=XB
      XB=XLB(I+1)
      GO TO 120
  140 CONTINUE
      GYL(J)=SUMG
  150 CONTINUE
      J=J+1
      IF(J > NYG) GO TO 160
      GYL(J)=0.D0 
      GO TO 150
  160 CONTINUE
      GO TO 300
      
  200 CONTINUE
      YB=YLB(J+1)
      IF(YB < XA) GO TO 210
      GYL(J)=0.D0 
      J=J+1
      IF(J > NYG) GO TO 260
      YA=YB
      GO TO 200
  210 CONTINUE
      XB=XLB(I+1)
      IF(XB < YA) GO TO 220
      I=I+1 
      IF(I > NXF) GO TO 260
      XA=XB
      GO TO 210
  220 CONTINUE
      XAYA=XA
      IF(YA < XA) XAYA=YA
      IF(YB < XB) GO TO 230
      PART=(YB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      GYL(J)=SUMG
      J=J+1
      IF(J > NYG) GO TO 260
      SUMG=0.D0
      SUMY=0.D0
      YA=YB
      YB=YLB(J+1)
      GO TO 220
  230 CONTINUE
      PART=(XB-XAYA)/(XB-XA)
      SUMG=SUMG+PART*FXL(I)
      SUMY=SUMY+PART
      I=I+1
      IF(I > NXF) GO TO 240
      XA=XB
      XB=XLB(I+1)
      GO TO 220
  240 CONTINUE
      GYL(J)=SUMG
  250 CONTINUE
      J=J+1
      IF(J > NYG) GO TO 260
      GYL(J)=0.D0
      GO TO 250
  260 CONTINUE

  300 CONTINUE
      RETURN
      END SUBROUTINE REPART

