!@sum ODEEP contains routines used for Qflux mixed layer with deep diff.
!@auth G. Schmidt/G. Russell
      MODULE ODEEP_COM
!@sum  ODEEP_COM defines the variables for deep diffusing Qflux model
!@auth Gavin Schmidt/Gary Russell
      USE RESOLUTION, only : im,jm
      IMPLICIT NONE
      SAVE

!@param LMOM number of layers for deep ocean diffusion
c      INTEGER, PARAMETER :: LMOM = 9    ! good for 1000m
      INTEGER, PARAMETER :: LMOM = 12    ! good for 5015m
!@var TG3M Monthly accumulation of temperatures at base of mixed layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TG3M ! (IM,JM,12)
!@var RTGO Temperature anomaly in thermocline
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: RTGO ! (LMOM,IM,JM)
      REAL*8 RTGO_diag(LMOM,IM,JM) ! used in diagnostic print routine
!@var STG3 accumulated temperature at base of mixed layer
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: STG3 ! (IM,JM)
!@var DTG3 accumulated temperature diff. from initial monthly values
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: DTG3 ! (IM,JM)
!@var EDO ocean vertical diffusion (m^2/s) - fixed, not distributed
      REAL*8, DIMENSION(IM,JM) :: EDO
!@var DZ thermocline layer thickness (m)
      REAL*8, DIMENSION(LMOM) :: DZ
!@var DZO,BYDZO distance between centres in thermocline layer (m)
      REAL*8, DIMENSION(LMOM-1) :: DZO,BYDZO

      END MODULE ODEEP_COM

      SUBROUTINE init_ODEEP(iniOCEAN)
!@sum  init_ODEEP initialise deep ocean arrays
!@auth G. Schmidt
      USE FILEMANAGER, only : openunit,closeunit
      USE RESOLUTION, only : im,jm
      USE ODEEP_COM, only : tg3m,stg3,dtg3,rtgo,dz,dzo,bydzo,edo,lmom
      USE STATIC_OCEAN, only : tocean
      USE DOMAIN_DECOMP_ATM, only : GRID, am_I_root, unpack_data

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: iniOCEAN
      INTEGER :: iu_tg3m,iu_EDDY,L
      CHARACTER*80 TITLE
      real*8, allocatable :: tg3m_glob(:,:,:)

!@param FAC ratio of adjacent deep ocean layers
C**** NOTE: For LMOM is 9 this value gives a total depth of 1000m
C**** For any different number of layers, the effective depth is given
C**** by the equation Z=10*(1-x^(LMOM-1))/(1-x).
C**** In particular, for LMOM=12, the total depth is 5015m
      REAL*8, PARAMETER :: FAC=1.705357255658901d0

C**** READ IN EDDY DIFFUSIVITY AT BASE OF MIXED LAYER
      CALL openunit("EDDY",iu_EDDY,.TRUE.,.TRUE.)
      CALL READT (iu_EDDY,0,IM*JM,EDO,1)
      call closeunit(iu_EDDY)

C**** DEFINE THE VERTICAL LAYERING EVERYWHERE EXCEPT LAYER 1 THICKNESS
      DZ(2)=10.
      DZO(1)=0.5*DZ(2)          ! 10./SQRT(FAC)
      BYDZO(1)=1./DZO(1)
      DO L=2,LMOM-1
        DZ(L+1)=DZ(L)*FAC
        DZO(L)=0.5*(DZ(L+1)+DZ(L)) !DZO(L-1)*FAC
        BYDZO(L)=1./DZO(L)
      END DO

C**** read in initial conditions and climatology for the temperatures at
C**** base of mixed layer, initialise deep arrays for start of run
      if (iniOCEAN) then
        stg3=0. ; dtg3=0. ; rtgo=0.

        call openunit("TG3M",iu_tg3m,.true.,.true.)
        allocate (tg3m_glob(im,jm,12))
        if(am_I_root()) then
          READ(iu_tg3m) TITLE,TG3M_glob
          WRITE(6,*) "Read from TG3M",TITLE
        end if
        call unpack_data(grid,tg3m_glob,tg3m)
        deallocate (tg3m_glob)
        call closeunit (iu_tg3m)

      end if
      return
C****
      end subroutine init_odeep

      SUBROUTINE io_ocean(kunit,iaction,ioerr)
!@sum  io_ocean reads and writes ocean arrays to file
!@auth Gavin Schmidt
      USE MODEL_COM, only : ioread,iowrite,irsficno,lhead
      USE STATIC_OCEAN
      USE ODEEP_COM
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : am_I_root,pack_data,unpack_data
      USE DOMAIN_DECOMP_1D, only : pack_column, unpack_column
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCNDEEP01"
      real*8, allocatable :: Tocean_glob(:,:,:),z1o_glob(:,:)
      real*8, allocatable, dimension(:,:) :: STG3_glob,DTG3_glob
      real*8, allocatable, dimension(:,:,:) :: RTGO_glob,TG3M_glob

      WRITE(MODULE_HEADER(lhead+1:80),'(a45,i2,a)') 'R8 To(3,ijm)'//
     *     ', dim(ijm): MixLD,Stg3,dtg3,rtgo(',lmom,',ijm),'//
     *     'tg3m(12,ijm)'

      allocate (Tocean_glob(3,im,jm),z1o_glob(im,jm))
      allocate (STG3_glob(im,jm),DTG3_glob(im,jm))
      allocate (RTGO_glob(lmom,im,jm),TG3M_glob(im,jm,12))

      SELECT CASE (IACTION)
      CASE (:IOWRITE)            ! output to standard restart file
!       WRITE (kunit,err=10) MODULE_HEADER,TOCEAN,Z1O,STG3,DTG3,RTGO,
!    *     TG3M
        call pack_column(grid,tocean,tocean_glob)
        call pack_data(grid,z1o , z1o_glob)
        call pack_data(grid,STG3,STG3_glob)
        call pack_data(grid,DTG3,DTG3_glob)
        call pack_column(grid,RTGO  ,  RTGO_glob)
        call pack_data(grid,tg3m,tg3m_glob)
        if(am_I_root()) WRITE (kunit,err=10) MODULE_HEADER,
     *    TOCEAN_glob,Z1O_glob,STG3_glob,DTG3_glob,RTGO_glob,TG3M_glob
        
      CASE (IOREAD:)            ! input from restart file
        SELECT CASE (IACTION)
        CASE (IRSFICNO)         ! no ocean initial conditions
C**** Note this is for reading in a full rsf file, but from a qflux run
C**** We do not check HEADER here because it will be wrong. The other
C**** data MUST be initialised by setting iniOCEAN=.TRUE. in init_ODEEP.
          if(am_I_root()) READ (kunit) HEADER,TOCEAN_glob,Z1O_glob
          call unpack_column(grid,tocean_glob,tocean)
          call unpack_data(grid,z1o_glob,z1o)
        CASE DEFAULT            ! restart file
          if(am_I_root()) then 
            READ (kunit,err=10) HEADER,TOCEAN_glob,Z1O_glob
     *             ,STG3_glob,DTG3_glob,RTGO_glob,TG3M_glob
            IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
            PRINT*,"Discrepancy in module version ",HEADER,MODULE_HEADER
            GO TO 10
            END IF
          end if
          call unpack_column(grid,tocean_glob,tocean)
          call unpack_data(grid,z1o_glob,z1o)
          call unpack_data(grid,STG3_glob,STG3)
          call unpack_data(grid,DTG3_glob,DTG3)
          call unpack_column(grid,RTGO_glob,RTGO)
          call unpack_data(grid,TG3M_glob,TG3M)
        END SELECT
      END SELECT
       
      deallocate(tocean_glob,z1o_glob,STG3_glob,DTG3_glob)
      deallocate(RTGO_glob,TG3M_glob)

      RETURN
 10   IOERR=1
      RETURN
C****
      END SUBROUTINE io_ocean

      SUBROUTINE io_ocdiag(kunit,it,iaction,ioerr)
!@sum  io_ocdiag reads and writes ocean diagnostic arrays to file
!@auth Gavin Schmidt
      USE MODEL_COM, only : ioread,iowrite,irsfic,irerun,iowrite_single
     *     ,ioread_single,lhead
      USE RESOLUTION, only : im,jm
      USE ODEEP_COM
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : am_I_root
      USE DOMAIN_DECOMP_1D, only : pack_column, unpack_column

      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var HEADER Character string label for individual records
      CHARACTER*80 :: HEADER, MODULE_HEADER = "OCDIAGDEEP01"
!@var it input/ouput value of hour
      INTEGER, INTENT(INOUT) :: it
!@var RTGO4 dummy variable for reading in single precision
      REAL*4 RTGO4(LMOM,IM,JM)
      real*8, allocatable :: rtgo_glob(:,:,:),rtgo_loc(:,:,:)
      
      INTEGER I,J

C**** no output required for rsf files. Only acc files
      write(MODULE_HEADER(lhead+1:80),'(a,i2,a)')
     *     'R4 RTGO(',lmom,'im,jm)'

      allocate(rtgo_glob(LMOM,IM,JM),
     *     rtgo_loc(LMOM,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                   grid%J_STRT_HALO:grid%J_STOP_HALO))
      SELECT CASE (IACTION)
      CASE (IOWRITE_SINGLE)     ! output to acc file
        call pack_column(grid,rtgo,rtgo_glob)
        if(am_I_root())
     *    WRITE (kunit,err=10) MODULE_HEADER,REAL(RTGO_glob,KIND=4)     
           
      CASE (ioread_single)    ! read in from acc file
        if(am_I_root()) then
           READ (kunit,err=10) HEADER,RTGO4
           rtgo_glob = RTGO4
           IF (HEADER(1:LHEAD).NE.MODULE_HEADER(1:LHEAD)) THEN
             PRINT*,"Discrepancy in module version ",HEADER
     *             ,MODULE_HEADER
             GO TO 10
           END IF
        end if
        call unpack_column(grid,rtgo_glob,rtgo_loc)
C**** sum RTGO over input files
        RTGO=RTGO+RTGO_loc
      END SELECT
 
      deallocate(rtgo_glob,rtgo_loc)

      RETURN

 10   IOERR=1
      deallocate(rtgo_glob,rtgo_loc)

      RETURN
C****
      END SUBROUTINE io_ocdiag

      SUBROUTINE reset_odiag(isum)
!@sum reset_odiag zeros out ocean diagnostics if needed
!@auth G. Schmidt
      USE ODEEP_COM, only : rtgo
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: isum

C**** there is a confusion of definition for rtgo
C****   i) for the model run, it is the actual temperature anomaly
C****  ii) for post-processing, it is the accumulatated anomaly
C**** Thus it is only initiallised here for case ii).
      if (isum.eq.1) rtgo=0

      return
      end subroutine reset_odiag

      SUBROUTINE conserv_OCE(OCEANE)
!@sum  conserv_OCE calculates ocean energy for Qflux ocean
!@auth Gavin Schmidt
      USE CONSTANT, only : shw,rhows
      USE RESOLUTION, only : im,jm
      USE FLUXES, only : focean
      USE GEOM, only : imaxj
      USE STATIC_OCEAN, only : tocean,z1o,z12o
      USE ODEEP_COM, only : dz,rtgo,lmom
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      IMPLICIT NONE
!@var OCEANE ocean energy (J/M^2)
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) :: OCEANE
      INTEGER I,J,L
      integer :: J_0, J_1, I_0,I_1
      logical :: HAVE_NORTH_POLE, HAVE_SOUTH_POLE

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &          HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &          HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      DO J=J_0, J_1
      DO I=I_0,IMAXJ(J)
        IF (FOCEAN(I,J).gt.0) THEN
          OCEANE(I,J)=(TOCEAN(1,I,J)*Z1O(I,J)
     *         +TOCEAN(2,I,J)*(Z12O(I,J)-Z1O(I,J)))*SHW*RHOWS
          DO L=2,LMOM
            OCEANE(I,J)=OCEANE(I,J)+(RTGO(L,I,J)*DZ(L)*SHW*RHOWS)
          END DO
        ELSE
          OCEANE(I,J)=0
        END IF
      END DO
      END DO
      IF (HAVE_SOUTH_POLE) OCEANE(2:im,1) =OCEANE(1,1)
      IF (HAVE_NORTH_POLE) OCEANE(2:im,JM)=OCEANE(1,JM)
C****
      END SUBROUTINE conserv_OCE

      SUBROUTINE ODIFS
!@sum  ODFIS calculates heat diffusion at the base of the mixed layer
!@+    compares that to the control run's temperature, calls odffus,
!@+    and reduces the upper ocean temperatures by the amount of heat
!@+    that is diffused into the thermocline
!@auth Gary Russell/G. Schmidt
!@calls ODFFUS
      USE FILEMANAGER
      USE CONSTANT, only : tf
      use TimeConstants_mod, only: SECONDS_PER_DAY, DAYS_PER_YEAR
      USE GEOM, only : imaxj
      USE ODEEP_COM, only : tg3m,rtgo,stg3,dtg3,edo,dz,dzo,bydzo,lmom
      USE SEAICE_COM, only : si_ocn
      USE DIAG_COM, only : aj,j_ftherm,itocean,itoice
      USE FLUXES, only : atmocn,focean
      USE STATIC_OCEAN, only : z12o,tocean
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      IMPLICIT NONE
      REAL*8, PARAMETER :: PERDAY=1./DAYS_PER_YEAR
!@param ALPHA degree of implicitness (1 fully implicit,0 fully explicit)
      REAL*8, PARAMETER :: ALPHA=.5d0
      REAL*8 :: ADTG3
      INTEGER I,J,L
      integer :: J_0, J_1, I_0,I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C****
C**** ACCUMULATE OCEAN TEMPERATURE AT MAXIMUM MIXED LAYER
C****
      DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          STG3(I,J)=STG3(I,J)+TOCEAN(3,I,J)
        END DO
      END DO
C****
C**** AT THE END OF EACH MONTH, UPDATE THE OCEAN TEMPERATURE
C**** DIFFERENCE AND REPLACE THE MONTHLY SUMMED TEMPERATURE
C****
      IF(JDATE.EQ.1) THEN
      DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          DTG3(I,J)=DTG3(I,J)+(STG3(I,J)-TG3M(I,J,JMON))
          TG3M(I,J,JMON)=STG3(I,J)
          STG3(I,J)=0.
        END DO
      END DO
      END IF
C****
C**** DIFFUSE THE OCEAN TEMPERATURE DIFFERENCE OF THE UPPER LAYERS
C**** INTO THE THERMOCLINE AND REDUCE THE UPPER TEMPERATURES BY THE
C**** HEAT THAT IS DIFFUSED DOWNWARD
C****
      DO J=J_0, J_1
        DO I=I_0,IMAXJ(J)
          IF(FOCEAN(I,J).GT.0.) THEN

            ADTG3=DTG3(I,J)*PERDAY
            RTGO(1,I,J)=ADTG3
C**** Set first layer thickness
            DZ(1)=Z12O(I,J)

            CALL ODFFUS (SECONDS_PER_DAY,ALPHA,EDO(I,J),DZ,BYDZO,
     &                   RTGO(1,I,J),LMOM)

            DO L=1,3
              TOCEAN(L,I,J)=TOCEAN(L,I,J)+(RTGO(1,I,J)-ADTG3)
            END DO
            AJ(J,J_FTHERM,ITOCEAN)=AJ(J,J_FTHERM,ITOCEAN)-(RTGO(1,I,J)
     *           -ADTG3)*Z12O(I,J)*FOCEAN(I,J)*(1.-si_ocn%RSI(I,J))
            AJ(J,J_FTHERM,ITOICE )=AJ(J,J_FTHERM,ITOICE )-(RTGO(1,I,J)
     *           -ADTG3)*Z12O(I,J)*FOCEAN(I,J)*si_ocn%RSI(I,J)
            atmocn%GTEMP(I,J) = TOCEAN(1,I,J)
            atmocn%GTEMP2(I,J) = TOCEAN(2,I,J)
            atmocn%GTEMPR(I,J)    = TOCEAN(1,I,J)+TF
          END IF
        END DO
      END DO

      RETURN
      END SUBROUTINE ODIFS

      SUBROUTINE ODFFUS (DT,ALPHA,ED,DZ,BYDZO,R,LMIJ)
!@sum  ODFFUS calculates the vertical mixing of a tracer
!@auth Gavin Schmidt/Gary Russell
!@calls TRIDIAG
      USE TRIDIAG_MOD, only : tridiag
      IMPLICIT NONE
!@var LMIJ IS THE NUMBER OF VERTICAL LAYERS
      INTEGER, INTENT(IN) :: LMIJ
!@var ED diffusion coefficient between adjacent layers (m**2/s)
!@var ALPHA determines the time scheme (0 explicit,1 fully implicit)
!@var DT time step (s)
      REAL*8, INTENT(IN) :: ED,ALPHA,DT
!@var DZ the depth of the layers (m)
!@var BYDZO is the inverse of depth between layer centers (1/m)
      REAL*8, INTENT(IN) :: DZ(LMIJ),BYDZO(LMIJ-1)
!@var R tracer concentration
      REAL*8, INTENT(INOUT) :: R(LMIJ)

      REAL*8 AM(LMIJ),BM(LMIJ),CM(LMIJ),DM(LMIJ)
      INTEGER L
C**** SET UP TRIDIAGONAL MATRIX ENTRIES AND RIGHT HAND SIDE
      AM(1)=0
      BM(1)=DZ(1)+ALPHA*DT*ED*BYDZO(1)
      CM(1)=     -ALPHA*DT*ED*BYDZO(1)
      DM(1)=DZ(1)*R(1)-(1.-ALPHA)*DT*ED*(R(1)-R(2))*BYDZO(1)

      DO L=2,LMIJ-1
        AM(L)=     -ALPHA*DT* ED*BYDZO(L-1)
        BM(L)=DZ(L)+ALPHA*DT*(ED*BYDZO(L-1)+ED*BYDZO(L))
        CM(L)=     -ALPHA*DT*               ED*BYDZO(L)
        DM(L)=DZ(L)*R(L)+(1.-ALPHA)*DT*(ED*(R(L-1)-R(L))*BYDZO(L-1)
     *                                 -ED*(R(L)-R(L+1))*BYDZO(L))
      END DO

      AM(LMIJ)=        -ALPHA*DT*ED*BYDZO(LMIJ-1)
      BM(LMIJ)=DZ(LMIJ)+ALPHA*DT*ED*BYDZO(LMIJ-1)
      CM(LMIJ)=0.
      DM(LMIJ)=DZ(LMIJ)*R(LMIJ)+(1.-ALPHA)*DT*ED*
     *         (R(LMIJ-1)-R(LMIJ))*BYDZO(LMIJ-1)

      CALL TRIDIAG(AM,BM,CM,DM,R,LMIJ)

      RETURN
      END SUBROUTINE ODFFUS

      SUBROUTINE CHECKO(SUBR)
!@sum  CHECKO Checks whether deep ocean values are reasonable
!@auth Original Development Team
      USE RESOLUTION, only : im,jm
      USE FLUXES, only : focean
      USE ODEEP_COM, only : lmom,stg3,dtg3,tg3m,rtgo
      USE STATIC_OCEAN, only : tocean
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds
      IMPLICIT NONE

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      LOGICAL QCHECKO
      INTEGER I,J,J_0,J_1,J_0H,J_1H,I_0,I_1,njpol

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     *     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      njpol = grid%J_STRT_SKP-grid%J_STRT

C**** Check for NaN/INF in ocean data
      CALL CHECK3C(TOCEAN(:,I_0:I_1,J_0:J_1),3,I_0,I_1,J_0,J_1,NJPOL,
     &     SUBR,'toc')
      CALL CHECK3B(DTG3(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1,
     &     SUBR,'dtg3')
      CALL CHECK3B(TG3M(I_0:I_1,J_0:J_1,:),I_0,I_1,J_0,J_1,NJPOL,12,
     &     SUBR,'tg3m')
      CALL CHECK3B(STG3(I_0:I_1,J_0:J_1),I_0,I_1,J_0,J_1,NJPOL,1 ,
     &     SUBR,'stg3')
      CALL CHECK3C(RTGO(:,I_0:I_1,J_0:J_1),LMOM,I_0,I_1,J_0,J_1,NJPOL,
     &     SUBR,'rtgo')

      QCHECKO = .FALSE.
C**** Check for reasonable values for ocean variables
      DO J=J_0, J_1
        DO I=I_0,I_1
          IF (focean(i,j)>0 .and.
     &        (TOCEAN(1,I,J).lt.-2. .or. TOCEAN(1,I,J).gt.50.)) THEN
            WRITE(6,*) 'After ',SUBR,': I,J,TOCEAN=',I,J,TOCEAN(1:3,I,J)
            QCHECKO = .TRUE.
          END IF
       END DO
      END DO
      IF (QCHECKO)
     *     call stop_model("CHECKO: Ocean variables out of bounds",255)

      END SUBROUTINE CHECKO

      SUBROUTINE diag_OCEAN_prep
      RETURN
      END SUBROUTINE diag_OCEAN_prep

      SUBROUTINE diag_OCEAN
!@sum  diag_OCEAN prints out diagnostics for ocean
!@$    ESMF: It should only be called from a serial region.
!@$          It is NOT parallelized.
!@auth Gavin Schmidt
      USE RESOLUTION, only : jm
      USE MODEL_COM, only : lrunid,xlabel,idacc
      USE FLUXES, only : focean
      USE GEOM, only : imaxj,lat_dg
      USE ODEEP_COM, only : lmom,rtgo=>rtgo_diag,dz
      USE DIAG_COM, only : qdiag,zoc
      USE MDIAG_COM, only : acc_period
     &     ,sname_strlen,units_strlen,lname_strlen
      USE DIAG_SERIAL, only : JLMAP
      IMPLICIT NONE
      CHARACTER(len=lname_strlen) :: LNAME
      CHARACTER(len=sname_strlen) :: SNAME
      CHARACTER(len=units_strlen) :: UNITS

      INTEGER I,J,L
      REAL*8 ATGO(JM,LMOM),SCALED,ONES(JM),focnj

C**** OPEN PLOTTABLE OUTPUT FILE IF DESIRED
      IF(QDIAG) call open_jl(trim(acc_period)//'.o'//XLABEL(1:LRUNID),
     *  jm,lmom,0,lat_dg)
      LNAME="Zonally averaged deep ocean temperature anomaly"
      SNAME="tgo_deep_anom"
      UNITS="DEGREES C"
C**** calculate zonal average
      DO L=1,LMOM
        DO J=1,JM
          ATGO(J,L)=0.
          focnj=0.
          DO I=1,IMAXJ(J)
            ATGO(J,L)=ATGO(J,L)+RTGO(L,I,J)*focean(i,j)
            focnj=focnj+focean(i,j)
          END DO
          if(focnj.gt.0.) ATGO(J,L)=ATGO(J,L)/focnj
        END DO
      END DO
C**** depths are calculated from base of the mixed layer
      ZOC(1)=0.
      DO L=2,LMOM
        ZOC(L)=ZOC(L-1)+DZ(L)
      END DO
      SCALED=1./IDACC(12)
      ONES(1:JM)=1.
C**** Print out a depth/latitude plot of the deep ocean temp anomaly
      CALL 
     *  JLMAP(LNAME,SNAME,UNITS,-1,ZOC,ATGO,SCALED,ONES,ONES,LMOM,2,1)
C****
      if(qdiag) call close_jl

      RETURN
      END SUBROUTINE diag_OCEAN

      SUBROUTINE ALLOC_ODEEP(grid)
      USE ODEEP_COM, only  : lmom,TG3M,RTGO,sTG3,dTG3
      USE DOMAIN_DECOMP_ATM, only : DIST_GRID,getDomainBounds
      IMPLICIT NONE
      INTEGER :: J_0H,J_1H,IER,I_0H,I_1H
      TYPE (DIST_GRID), INTENT(IN) :: grid

      call getDomainBounds(GRID,J_STRT_HALO=J_0H,J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      ALLOCATE(TG3M(I_0H:I_1H,J_0H:J_1H,12),
     &         RTGO(LMOM,I_0H:I_1H,J_0H:J_1H),
     &         sTG3(I_0H:I_1H,J_0H:J_1H),
     &         dTG3(I_0H:I_1H,J_0H:J_1H),
     &    STAT=IER)

      END SUBROUTINE ALLOC_ODEEP

      SUBROUTINE gather_odiags ()
!@sum  collect the local acc-arrays into global arrays
!@+    run-time
!@auth Reto Ruedy

      USE ODEEP_COM, only  : RTGO,RTGO_diag
      use domain_decomp_atm, only : grid, pack_column

      IMPLICIT NONE
      
      call pack_column (grid, RTGO, RTGO_diag)

      END SUBROUTINE gather_odiags

