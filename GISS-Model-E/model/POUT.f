#include "rundeck_opts.h"

!@sum  POUT default output routines for standard GISS formats
!@auth Gavin Schmidt

C****
C**** If other formats are desired, please replace these routines
C**** with ones appropriate for your chosen output format, but with the
C**** same interface
C****
C**** Note: it would be nice to amalgamate IL and JL, but that will
C**** have to wait.


      module gissout
!@sum gissout contains variables for outputting GISS format binaries
!@auth G. Schmidt
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      implicit none
!@var iu_ij,iu_jl,iu_il,iu_j,iu_diurn,iu_hdiurn,iu_jc !units for selected 
!@+   diag. output
      integer iu_ij,iu_ijk,iu_il,iu_j,iu_jl,iu_diurn,iu_hdiurn,iu_isccp
     *     ,iu_ijl,iu_jc
!@var im,jm,lm,lm_req local dimensions set in open_* routines
      integer :: im,jm,lm,lm_req,ndiuvar
!@var JMMAX maximum conceivable JM
      INTEGER, PARAMETER :: JMMAX=200
!@var LAT_DG latitude of mid points of primary and sec. grid boxs (deg)
      REAL*8, DIMENSION(JMMAX,2) :: LAT_DG

      end module gissout

      subroutine open_ij(filename,im_gcm,jm_gcm)
!@sum  OPEN_IJ opens the lat-lon binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm

C**** open unit for ij
      call openunit(filename,iu_ij,.true.,.false.)

C**** set dimensions
      im=im_gcm
      jm=jm_gcm

      return
      end subroutine open_ij

      subroutine close_ij
!@sum  CLOSE_IJ closes the lat-lon binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ij)
      return
      end subroutine close_ij

      subroutine POUT_IJ(TITLE,SNAME,LNAME,UNITS,XIJ,XJ,XSUM,
     &     IGRID,JGRID)
!@sum  POUT_IJ output lat-lon binary records
!@auth Gavin Schmidt
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS
!@var XIJ lat/lon output field
      REAL*8, DIMENSION(IM,JM), INTENT(IN) :: XIJ
!@var XJ lat sum/mean of output field
      REAL*8, DIMENSION(JM), INTENT(IN) :: XJ
!@var XSUM global sum/mean of output field
      REAL*8, INTENT(IN) :: XSUM
!@var IGRID,JGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IGRID,JGRID

      WRITE(iu_ij) TITLE,REAL(XIJ,KIND=4),REAL(XJ,KIND=4),
     *             REAL(XSUM,KIND=4)
      return
      end subroutine POUT_IJ

      subroutine open_wp(filename,jm_gcm,lm_gcm,lm_req_gcm,lat_dg_gcm)
!@sum  OPEN_WP uses OPEN_JL to open the wave power binary output file
!@auth M. Kelley
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var LM_GCM,JM_GCM,lm_req_gcm dimensions for jl output
      INTEGER, INTENT(IN) :: lm_gcm,jm_gcm,lm_req_gcm
!@var lat_dg_gcm the "horizontal" coordinate for pout_jl
      REAL*8, INTENT(IN), DIMENSION(JM_GCM,2) :: lat_dg_gcm

      call open_jl(filename,jm_gcm,lm_gcm,lm_req_gcm,lat_dg_gcm)

      return
      end subroutine open_wp

      subroutine close_wp
!@sum  CLOSE_WP uses CLOSE_JL to close the wave power binary output file
!@auth M. Kelley
      IMPLICIT NONE
      call close_jl
      return
      end subroutine close_wp

      subroutine POUT_WP(TITLE,LNAME,SNAME,UNITS_IN,
     &     J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_WP calls POUT_JL to output wave power binary records
!@auth M. Kelley
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var J1,KLMAX variables required by pout_jl
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field, dimensioned to accomodate pout_jl expectations
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
!@var PM the "vertical" coordinate for pout_jl
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM
!@var CX,CY x,y coordinate names
      CHARACTER*16, INTENT(IN) :: CX,CY
      call POUT_JL(TITLE,LNAME,SNAME,UNITS_IN,
     &     J1,KLMAX,XJL,PM,CX,CY)
      return
      end subroutine pout_wp

      subroutine open_jl(filename,jm_gcm,lm_gcm,lm_req_gcm,lat_dg_gcm)
!@sum  OPEN_JL opens the lat-height binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var LM_GCM,JM_GCM,lm_req_gcm dimensions for jl output
      INTEGER, INTENT(IN) :: lm_gcm,jm_gcm,lm_req_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM,2) :: lat_dg_gcm

      call openunit(filename,iu_jl,.true.,.false.)

C**** set dimensions
      jm=jm_gcm
      lm=lm_gcm
      lm_req=lm_req_gcm
      lat_dg(1:JM,:)=lat_dg_gcm(1:JM,:)

      return
      end subroutine open_jl

      subroutine close_jl
!@sum  CLOSE_JL closes the lat-height binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_jl)
      return
      end subroutine close_jl

      subroutine POUT_JL(TITLE,LNAME,SNAME,UNITS_IN,
     &     J1,KLMAX,XJL,PM,CX,CY)
!@sum  POUT_JL output lat-height binary records
!@auth Gavin Schmidt
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS_IN
!@var KLMAX max level to output
!@var J1 minimum j value to output (needed for secondary grid fields)
      INTEGER, INTENT(IN) :: KLMAX,J1
!@var XJL output field
!@+       (J1:JM,1:KLMAX) is field
!@+       (JM+1:JM+3,1:KLMAX) are global/NH/SH average over L
!@+       (J1:JM+3,LM+LM_REQ+1) are averages over J
      REAL*8, DIMENSION(JM+3,LM+LM_REQ+1), INTENT(IN) :: XJL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(KLMAX), INTENT(IN) :: PM

      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = '                '
      REAL*8 XCOOR(JM)
      INTEGER J,L,JXMAX

      JXMAX = JM-J1+1
      XCOOR(1:JXMAX) = LAT_DG(J1:JM,J1)

      WRITE (iu_jl) TITLE,JXMAX,KLMAX,1,1,
     *     ((REAL(XJL(J1+J-1,L),KIND=4),J=1,JXMAX),L=1,KLMAX)
     *     ,(REAL(XCOOR(J),KIND=4),J=1,JXMAX)
     *     ,(REAL(PM(L),KIND=4),L=1,KLMAX)
     *     ,REAL(1.,KIND=4),REAL(1.,KIND=4)
     *     ,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(REAL(XJL(J,LM+LM_REQ+1),KIND=4),J=J1,JM+3)
     *     ,((REAL(XJL(J,L),KIND=4),J=JM+1,JM+3),L=1,KLMAX)

      return
      end

      subroutine open_il(filename,im_gcm,lm_gcm,lm_req_gcm)
!@sum  OPEN_IL opens the lon-height binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,LM_GCM,lm_req_gcm dimensions for il output
      INTEGER, INTENT(IN) :: im_gcm,lm_gcm,lm_req_gcm

      call openunit(filename,iu_il,.true.,.false.)

C**** set units
      im=im_gcm
      lm=lm_gcm
      lm_req=lm_req_gcm

      return
      end subroutine open_il

      subroutine close_il
!@sum  CLOSE_IL closes the lon-height binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_il)
      return
      end subroutine close_il

      subroutine POUT_IL(TITLE,sname,lname,unit,I1,ISHIFT,KLMAX,XIL
     *     ,PM,CX,CY,ASUM,GSUM,ZONAL)
!@sum  POUT_IL output lon-height binary records
!@auth Gavin Schmidt
      USE GISSOUT
#ifndef CUBED_SPHERE
      USE GEOM, only : lon_dg
#endif
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, INTENT(IN) :: TITLE*80
!@var KLMAX max level to output
!@var ISHIFT flag for secondary grid
!@var I1 coordinate index associated with first long. (for wrap-around)
      INTEGER, INTENT(IN) :: KLMAX,ISHIFT,I1
!@var XIL output field
      REAL*8, DIMENSION(IM,LM+LM_REQ+1), INTENT(IN) :: XIL
!@var PM pressure levels (MB)
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: PM
!@var ASUM vertical mean/sum
      REAL*8, DIMENSION(IM), INTENT(IN) :: ASUM
!@var GSUM total sum/mean
      REAL*8, INTENT(IN) :: GSUM
!@var ZONAL zonal sum/mean
      REAL*8, DIMENSION(LM+LM_REQ), INTENT(IN) :: ZONAL

      character(len=sname_strlen), intent(in) :: sname
      character(len=units_strlen), intent(in) :: unit
      character(len=lname_strlen), intent(in) :: lname
      CHARACTER*16, INTENT(IN) :: CX,CY
      CHARACTER*16, PARAMETER :: CBLANK = ' '
      REAL*8 XCOOR(IM)
      INTEGER I,L

C**** Allow for the possibility of wrap-around arrays
#ifndef CUBED_SPHERE
      XCOOR(1:IM-I1+1) = LON_DG(I1:IM,ISHIFT)
      IF (I1.gt.1) XCOOR(IM-I1+2:IM) = LON_DG(1:I1-1,ISHIFT)
#endif

      WRITE (iu_il) TITLE,IM,KLMAX,1,1,
     *     ((REAL(XIL(I,L),KIND=4),I=1,IM),L=1,KLMAX)
     *     ,(REAL(XCOOR(I),KIND=4),I=1,IM)
     *     ,(REAL(PM(L),KIND=4),L=1,KLMAX)
     *     ,REAL(0.,KIND=4),REAL(0.,KIND=4)
     *     ,CX,CY,CBLANK,CBLANK,'NASAGISS'
     *     ,(REAL(ASUM(I),KIND=4),I=1,IM),REAL(GSUM,KIND=4)
     *     ,(REAL(ZONAL(L),KIND=4),L=1,KLMAX)

      return
      end subroutine POUT_IL

      subroutine close_j
!@sum  CLOSE_J closes the latitudinal budget-page ascii output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_j)
      return
      end subroutine close_j

      subroutine open_j(filename,ntypes,jm_gcm,lat_dg_gcm)
!@sum  OPEN_J opens the latitudinal budget-page ascii output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var ntypes number of surface types to be output
      integer, intent(in) :: ntypes
!@var JM_GCM dimensions for j output
      INTEGER, INTENT(IN) :: jm_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM) :: lat_dg_gcm

      call openunit(filename,iu_j,.false.,.false.)

C**** set dimensions
      jm=jm_gcm
      lat_dg(1:JM,1)=lat_dg_gcm(1:JM)

      return
      end subroutine open_j

      subroutine POUT_J(TITLE,SNAME,LNAME,UNITS,BUDG,KMAX,TERRAIN,
     *     iotype)
!@sum  POUT_J output zonal budget ascii file (aplot format)
!@auth Gavin Schmidt
      USE GISSOUT
      USE DIAG_COM, only : KAJ
      IMPLICIT NONE
      CHARACTER*16, DIMENSION(KAJ),INTENT(INOUT) :: TITLE
      CHARACTER*16 :: NEWTIT
!@var LNAME,SNAME,UNITS dummy strings
      CHARACTER(len=lname_strlen), DIMENSION(KAJ),INTENT(IN) :: LNAME
      CHARACTER(len=sname_strlen), DIMENSION(KAJ),INTENT(IN) :: SNAME
      CHARACTER(len=units_strlen), DIMENSION(KAJ),INTENT(IN) :: UNITS
      CHARACTER*16, INTENT(IN) :: TERRAIN
      REAL*8, DIMENSION(JM+3,KAJ), INTENT(IN) :: BUDG
      INTEGER, INTENT(IN) :: KMAX,iotype
      INTEGER K,N,J,n1

C**** Convert spaces in TITLE to underscore
C**** Try simply removing spaces for compactness
      DO K=1,KMAX
        newtit=' '
        n1=1
        do n=2,len_trim(title(K))    ! skip control character
          if (title(K)(n:n).ne.' ') then
            n1=n1+1
            newtit(n1:n1)=title(K)(n:n)
          end if
        end do
        title(K)=newtit
      END DO

      WRITE(iu_j,*) "Zonal Budgets for surface type ",TERRAIN
      WRITE(iu_j,*) "Latitude"
      WRITE(iu_j,*) "Zonal Average"
      WRITE(iu_j,'(A4,100A)') "Lat",(TRIM(TITLE(K)(1:14)),K=1,KMAX)

      DO J=1,JM
        WRITE(iu_j,'(I4,100(1X,F8.3))') NINT(LAT_DG(J,1)),
     *       (BUDG(J,K),K=1,KMAX)
      END DO
      WRITE(iu_j,*)
C**** output hemispheric and global means
      WRITE(iu_j,'(A4,100F8.3)') "NH",(BUDG(JM+1,K),K=1,KMAX)
      WRITE(iu_j,'(A4,100F8.3)') "SH",(BUDG(JM+2,K),K=1,KMAX)
      WRITE(iu_j,'(A4,100F8.3)') "GLOB",(BUDG(JM+3,K),K=1,KMAX)
      WRITE(iu_j,*)

      return
      end

      subroutine open_ijk(filename,im_gcm,jm_gcm,lm_gcm)
!@sum  OPEN_IJK opens the lat-lon-height binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM,LM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm,lm_gcm

      call openunit(filename,iu_ijk,.true.,.false.)

C**** set dimensions
      im=im_gcm
      jm=jm_gcm
      lm=lm_gcm

      return
      end subroutine open_ijk

      subroutine close_ijk
!@sum  CLOSE_IJK closes the lat-lon-height binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ijk)
      return
      end subroutine close_ijk

      subroutine POUT_IJK(TITLE,SNAME,LNAME,UNITS,XIJK,XJK,XK,IJGRID)
!@sum  POUT_IJK outputs lat-lon-height binary records
!@auth M. Kelley
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(IN) :: XIJK
!@var XJK lat sum/mean of output field
      REAL*8, DIMENSION(JM,LM), INTENT(IN) :: XJK
!@var XK global sum/mean of output field
      REAL*8, DIMENSION(LM), INTENT(IN) :: XK
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID
      INTEGER :: K

      DO K=1,LM
         WRITE(iu_ijk) TITLE(K), REAL(XIJK(:,:,K),KIND=4),
     &     REAL(XJK(:,K),KIND=4), REAL(XK(K),KIND=4)
      ENDDO
      return
      end subroutine POUT_IJK

      subroutine open_ijl(filename,im_gcm,jm_gcm,lm_gcm,
     &     kaijl,name_ijl,lname_ijl,units_ijl,lgrid_ijl)
!@sum  OPEN_IJK opens the lat-lon-layer binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM,LM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: im_gcm,jm_gcm,lm_gcm
! the following are only for plug compatibility with the netcdf version
      integer :: kaijl
      CHARACTER(len=lname_strlen), DIMENSION(kaijl) :: lname_ijl
      CHARACTER(len=sname_strlen), DIMENSION(kaijl) :: name_ijl
      CHARACTER(len=units_strlen), DIMENSION(kaijl) :: units_ijl
      integer, dimension(kaijl) :: lgrid_ijl

      call openunit(filename,iu_ijl,.true.,.false.)

C**** set dimensions
      im=im_gcm
      jm=jm_gcm
      lm=lm_gcm

      return
      end subroutine open_ijl

      subroutine close_ijl
!@sum  CLOSE_IJL closes the lat-lon-layer binary output file
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_ijl)
      return
      end subroutine close_ijl

      subroutine POUT_IJL(TITLE,SNAME,LNAME,UNITS,XIJL,XJL,XL,IJGRID)
!@sum  POUT_IJL outputs lat-lon-layer binary records
!@auth M. Kelley
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(IN) :: XIJL
!@var XJK lat sum/mean of output field
      REAL*8, DIMENSION(JM,LM), INTENT(IN) :: XJL
!@var XK global sum/mean of output field
      REAL*8, DIMENSION(LM), INTENT(IN) :: XL
!@var IJGRID = 1 for primary lat-lon grid, 2 for secondary lat-lon grid
      INTEGER, INTENT(IN) :: IJGRID
      INTEGER :: K

      DO K=1,LM
         WRITE(iu_ijl) TITLE(K), REAL(XIJL(:,:,K),KIND=4),
     &     REAL(XJL(:,K),KIND=4), REAL(XL(K),KIND=4)
      ENDDO
      return
      end subroutine POUT_IJL

      subroutine open_isccp(filename,ntau,npres,nisccp)
!@sum  OPEN_ISCCP opens the binary output file of ISCCP histograms
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var IM_GCM,JM_GCM,LM_GCM dimensions for ij output
      INTEGER, INTENT(IN) :: ntau,npres,nisccp

      call openunit(filename,iu_isccp,.true.,.false.)

C**** set dimensions
      im=ntau
      jm=npres
      lm=nisccp

      return
      end subroutine open_isccp

      subroutine close_isccp
!@sum  CLOSE_ISCCP closes the binary output file of ISCCP histograms
!@auth M. Kelley
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_isccp)
      return
      end subroutine close_isccp

      subroutine POUT_ISCCP(TITLE,SNAME,LNAME,UNITS,XIJK,TAUM,PRES)
!@sum  POUT_ISCCP outputs tau-height-lat binary output file of ISCCP histograms
!@auth M. Kelley
      USE GISSOUT
      IMPLICIT NONE
!@var TITLE 80 byte title including description and averaging period
      CHARACTER, DIMENSION(LM), INTENT(IN) :: TITLE*80
!@var SNAME short name of field
      CHARACTER(len=sname_strlen), INTENT(IN) :: SNAME
!@var LNAME long name of field
      CHARACTER(len=lname_strlen), INTENT(IN) :: LNAME
!@var UNITS units of field
      CHARACTER(len=units_strlen), INTENT(IN) :: UNITS
!@var XIJK lat/lon/height output field
      REAL*8, DIMENSION(IM,JM,LM), INTENT(IN) :: XIJK
      REAL*8, INTENT(IN) :: taum(im), pres(jm)
      INTEGER :: K
      CHARACTER*16, PARAMETER ::
     &     CX = 'OPTICAL DEPTH   ',
     &     CY = 'PRESSURE        ',
     &     CBLANK = '                '
      do k=1,lm
c restored "GISS 4D format"
        write(iu_isccp) title(k),im,jm,1,1
     &       ,real(xijk(:,:,k),kind=4)
     &       ,real(taum(1:im),kind=4)
     &       ,real(pres(1:jm),kind=4)
     &       ,real(1.,kind=4),real(1.,kind=4)
     &       ,cx,cy,cblank,cblank,'NASAGISS'

      enddo
      return
      end subroutine POUT_ISCCP

      subroutine close_diurn
!@sum  CLOSE_DIURN closes the hourly diurnal_cycle ascii output file
!@auth J. Lerner     
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_diurn)
      return
      end subroutine close_diurn

      subroutine open_diurn(filename,hr_in_day,NDIUVAR_gcm,kr1,kr2)
!@sum  OPEN_DIURN opens the hourly diurnal_cycle ascii output file
!@auth J. Lerner     
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: kr1,kr2 ! dummy variables
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: hr_in_day,NDIUVAR_gcm
      call openunit(filename,iu_diurn,.false.,.false.)
      NDIUVAR = NDIUVAR_gcm
      return
      end subroutine open_diurn

      subroutine POUT_diurn(SNAME,NAME,UNITS,FHOUR,NAMDD,IJDD1,IJDD2,
     &     HR_IN_DAY,kp)
!@sum  POUT_diurn output hourly diurnal_cycle ascii file (aplot format)
!@auth J. Lerner     
      USE GISSOUT
      IMPLICIT NONE
!@var NAME,UNITS dummy strings
      CHARACTER*16, DIMENSION(NDIUVAR) :: UNITS,NAME,SNAME
      CHARACTER*4, INTENT(IN) :: NAMDD !names of boxes
      INTEGER, INTENT(IN) :: HR_IN_DAY,KP,IJDD1,IJDD2
      REAL*8, DIMENSION(HR_IN_DAY+1,NDIUVAR), INTENT(IN) :: FHOUR
      INTEGER K,N,I

C**** Convert spaces in TITLE to underscore
C**** Try simply removing spaces for compactness
      DO K=1,kp
        do n=2,len_trim(name(K))    ! skip leading blank
          if (name(k)(n:n).eq.' ') name(k)(n:n)='_'
        end do
        do n=2,len_trim(units(K))    
          if (units(k)(n:n).eq.' ') units(k)(n:n)='_'
        end do
      END DO

      WRITE(iu_diurn,*) "Hourly Means for Region ",NAMDD,' at (',
     &    IJDD1,',',IJDD2,')'
      WRITE(iu_diurn,*) "Hour"
      WRITE(iu_diurn,*) "Hourly Mean"
      WRITE(iu_diurn,'(A4,100A)') "Hour",(NAME(K),K=1,kp)

      DO I=1,HR_IN_DAY
        WRITE(iu_diurn,'(I4,100(1X,F8.3))') I,(FHOUR(i,K),K=1,kp)
      END DO
      WRITE(iu_diurn,*)
C**** output daily mean
      WRITE(iu_diurn,'(A4,100F8.3)') "AVE",(FHOUR(HR_IN_DAY+1,K),K=1,kp)
      WRITE(iu_diurn,*)
      WRITE(iu_diurn,'(A4,100A)') "Units",(UNITS(K),K=1,kp)
      WRITE(iu_diurn,*)
      WRITE(iu_diurn,*)

      return
      end subroutine POUT_diurn


      subroutine close_hdiurn
!@sum  CLOSE_HDIURN closes the hourly diurnal_cycle ascii output file
!@auth J. Lerner
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_hdiurn)
      return
      end subroutine close_hdiurn

      subroutine open_hdiurn(filename,hr_in_month,NDIUVAR_gcm,kr1,kr2)
!@sum  OPEN_HDIURN opens the hourly diurnal_cycle ascii output file
!@auth J. Lerner
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: kr1,kr2 ! dummy variables
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: hr_in_month,NDIUVAR_gcm
      call openunit(filename,iu_hdiurn,.false.,.false.)
      NDIUVAR = NDIUVAR_gcm
      return
      end subroutine open_hdiurn

      subroutine POUT_hdiurn(SNAME,NAME,UNITS,FHOUR,NAMDD,IJDD1,IJDD2,
     &     HR_IN_period,kp)
!@sum  POUT_hdiurn output hourly diurnal_cycle ascii file (aplot format)
!@auth J. Lerner
      USE GISSOUT
      IMPLICIT NONE
!@var NAME,UNITS dummy strings
      CHARACTER*16, DIMENSION(NDIUVAR) :: UNITS,NAME,SNAME
      CHARACTER*4, INTENT(IN) :: NAMDD !names of boxes
      INTEGER, INTENT(IN) :: HR_IN_period,KP,IJDD1,IJDD2
      REAL*8, DIMENSION(HR_IN_period,NDIUVAR), INTENT(IN) :: FHOUR
      INTEGER K,N,I

C**** Convert spaces in TITLE to underscore
C**** Try simply removing spaces for compactness
      DO K=1,kp
        do n=2,len_trim(name(K))    ! skip leading blank
          if (name(k)(n:n).eq.' ') name(k)(n:n)='_'
        end do
        do n=2,len_trim(units(K))
          if (units(k)(n:n).eq.' ') units(k)(n:n)='_'
        end do
      END DO

      WRITE(iu_hdiurn,*) "Hourly Value for Region ",NAMDD,' at (',
     &    IJDD1,',',IJDD2,')'
      WRITE(iu_hdiurn,*) "Hour"
      WRITE(iu_hdiurn,*) "Hourly Value"
      WRITE(iu_hdiurn,'(A4,100A)') "Hour",(NAME(K),K=1,kp)

      DO I=1,HR_IN_period
        WRITE(iu_hdiurn,'(I4,100(1X,F9.3))') I,(FHOUR(i,K),K=1,kp)
      END DO
      WRITE(iu_hdiurn,*)

      return
      end subroutine POUT_hdiurn

      subroutine close_jc
!@sum  CLOSE_jc closes the conservation quantity ascii output file
!@auth J. Lerner
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
      call closeunit(iu_jc)
      return
      end subroutine close_jc

      subroutine open_jc(filename,jm_gcm,lat_dg_gcm)
!@sum  OPEN_jc opens the conservation quantity ascii output file
!@auth J. Lerner
      USE GISSOUT
      USE FILEMANAGER
      IMPLICIT NONE
!@var FILENAME output file name
      CHARACTER*(*), INTENT(IN) :: filename
!@var JM_GCM dimensions for j output
      INTEGER, INTENT(IN) :: jm_gcm
!@var lat_dg_gcm latitude of mid points of grid boxs (deg)
      REAL*8, INTENT(IN), DIMENSION(JM_GCM) :: lat_dg_gcm

      call openunit(filename,iu_jc,.false.,.false.)

C**** set dimensions
      jm=jm_gcm
      lat_dg(1:JM,1)=lat_dg_gcm(1:JM)

      return
      end subroutine open_jc

      subroutine POUT_jc(TITLE,SNAME,LNAME,UNITS,cnslat,KMAX)
!@sum  POUT_JC output zonal conservation ascii file (aplot format)
!@auth J. Lerner
      USE GISSOUT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KMAX
      CHARACTER*38, DIMENSION(kmax),INTENT(INOUT) :: TITLE
!@var LNAME,SNAME,UNITS dummy strings
      CHARACTER(len=lname_strlen), DIMENSION(kmax),INTENT(IN) :: LNAME
      CHARACTER(len=sname_strlen), DIMENSION(kmax),INTENT(IN) :: SNAME
      CHARACTER(len=units_strlen), DIMENSION(kmax),INTENT(IN) :: UNITS
      REAL*8, DIMENSION(JM+3,kmax), INTENT(IN) :: cnslat
      INTEGER K,N,J,KK,K1,K2

C**** Convert spaces in TITLE to underscore
      DO K=1,KMAX
        title(K)(1:1)=' '
        do n=2,len_trim(title(K))    ! skip control character
          if (title(k)(n:n).eq.' ') title(k)(n:n)='_'
        end do
      END DO

      DO KK=1,KMAX/50+1   ! split into 50 parameter chunks
      K1=1+(KK-1)*50
      K2=MIN(KMAX,KK*50)
      WRITE(iu_jc,*) "Zonal conservation diagnostics"
      WRITE(iu_jc,*) "Latitude"
      WRITE(iu_jc,*) "Zonal Average"
      WRITE(iu_jc,'(A,100A)') "Lat",(TRIM(TITLE(K)(1:38)),K=K1,K2)

      DO J=1,JM
        WRITE(iu_jc,'(F6.1,50(1X,1pE11.4))') LAT_DG(J,1),
     *       (cnslat(J,K),K=K1,K2)
      END DO
      WRITE(iu_jc,*)
C**** output hemispheric and global means
      WRITE(iu_jc,'(A4,50(1X,1pE11.4))')"NH",(cnslat(JM+2,K),K=K1,K2)
      WRITE(iu_jc,'(A4,50(1X,1pE11.4))')"SH",(cnslat(JM+1,K),K=K1,K2)
      WRITE(iu_jc,'(A4,50(1X,1pE11.4))')"GLOB",(cnslat(JM+3,K),K=K1,K2)
      WRITE(iu_jc,*)
      END DO

      return
      end subroutine pout_jc


