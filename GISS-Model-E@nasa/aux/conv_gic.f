      program convert_gic_files
C**** convert GIC files from model II' (B399) format to modelE format
C**** compile with: gmake conv_gic.o
C**** f90 -o conv_gic conv_gic.o ../model/ *.o -O2 -64 -mips4 \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE CONSTANT, only : lhm,shi
      USE MODEL_COM, only : im,jm,lm,iowrite,focean
      USE GHY_COM, only : snowe,tearth,wearth,aiearth,snoage,wbare,wvege
     *     ,htbare,htvege,snowbv,ngm,evap_max_ij,fr_sat_ij,qg_ij
      use veg_com, only : cint,qfol,cnc_ij
      USE STATIC_OCEAN, only : tocean,z1o
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi,pond_melt,flag_dsws
      USE SEAICE, only : ace1i,xsi,ac2oim,ssi0
      USE LANDICE_COM, only : tlandi,snowli
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60
      INTEGER IARGC,iu_GIC,I,J,L,N,ioerr,kunit,iaction
      REAL*8 MSI1,X
      INTEGER*4 len

      IF (IARGC().lt.2) THEN
        PRINT*,"Convert GIC files from old format to new"
        PRINT*,"conv_gic filename output_file"
        STOP
      END IF

      CALL GETARG(1,infile)
      CALL GETARG(2,outfile)

      iu_GIC=9

      print*,trim(infile)

      OPEN(iu_GIC,FILE=trim(infile),FORM="UNFORMATTED",STATUS="OLD")

C**** Note that old GDATA file only has GDATA(1:14)!!

      READ(iu_GIC,END=810) SNOWI,SNOWE,
     *     ((HSI(1,I,J),I=1,IM),J=1,JM),TEARTH,WEARTH,AIEARTH,
     *     ((HSI(2,I,J),I=1,IM),J=1,JM),((X,I=1,IM),J=1,JM),
     *     (((SNOAGE(L,I,J),I=1,IM),J=1,JM),L=1,3),SNOWLI,
     *     (((TLANDI(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     (((HSI(L,I,J),I=1,IM),J=1,JM),L=3,4),
     *     (((wbare(L,I,J),I=1,IM),J=1,JM),L=1,NGM),
     *     (((wvege(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((htbare(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((htvege(L,I,J),I=1,IM),J=1,JM),L=0,NGM),
     *     (((snowbv(L,I,J),I=1,IM),J=1,JM),L=1,2),
     *     ((TOCEAN(1,I,J),I=1,IM),J=1,JM),RSI
      CLOSE (iu_GIC)


C**** define sea ice defaults: Use AC2OIM instead of MSI.
      DO J=1,JM
        DO I=1,IM
          MSI1=SNOWI(I,J)+ACE1I
          MSI(I,J) = AC2OIM
          RSI(I,J) = MIN(MAX(0d0,RSI(I,J)),1d0)
          TOCEAN(2:3,I,J) = TOCEAN(1,I,J)
          POND_MELT(I,J)=0.
          FLAG_DSWS(I,J)=.FALSE.
          IF (FOCEAN(I,J).gt.0) THEN
            SSI(3:4,I,J)=SSI0 * XSI(3:4)*MSI(I,J)
            IF (ACE1I*XSI(1).gt.SNOWI(I,J)*XSI(2)) THEN
              SSI(1,I,J)=SSI0 * (ACE1I-(ACE1I+SNOWI(I,J))* XSI(2))
              SSI(2,I,J)=SSI0 * (ACE1I+SNOWI(I,J))* XSI(2)
            ELSE
              SSI(1,I,J)=0.
              SSI(2,I,J)=SSI0 * ACE1I
            END IF
          ELSE
            SSI(1:4,I,J) = 0
          END IF
          HSI(1:2,I,J) = (SHI*MIN(0d0,HSI(1:2,I,J))-LHM)*XSI(1:2)*MSI1
     *         +LHM*SSI(1:2,I,J)
          HSI(3:4,I,J) = (SHI*MIN(0d0,HSI(3:4,I,J))-LHM)*XSI(3:4)*AC2OIM
     *         +LHM*SSI(3:4,I,J) 
        END DO
      END DO

      print*,TOCEAN(1,JM-1,2),RSI(45,2),SNOWI(45,2)

C**** set default values for evaporation limiting arrays
      evap_max_ij(:,:) = 1.d0
      fr_sat_ij(:,:) = 1.d0
      qg_ij(:,:) = 0.d0

C**** set default foliage values
ccc   these have been moved out of GIC
c      Qfol=3.D-6
c      Cint=0.0127D0
c      cnc_ij=0.d0

      OPEN(iu_GIC,FILE=trim(outfile),
     *     FORM="UNFORMATTED",STATUS="UNKNOWN")

      kunit=iu_GIC
      iaction=iowrite
      ioerr=-1
      call io_ocean  (kunit,iaction,ioerr)
      call io_seaice (kunit,iaction,ioerr)
      call io_earth  (kunit,iaction,ioerr)
      call io_soils  (kunit,iaction,ioerr)
      call io_landice(kunit,iaction,ioerr)
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT="
     *     ,kunit
      close (kunit)

      print*,"New GIC file written out to ",trim(outfile)
      stop
 800  print*,"Error reading in file"
      stop
 810  print*,"Premature End of file"
      stop
      end
