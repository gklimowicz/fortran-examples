      MODULE AERO_SUBS


!@sum     This module contains various aerosol microphysical routines.
!@auth    Susanne Bauer/Doug Wright
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_PARAM
      USE AERO_CONFIG
      IMPLICIT NONE

      CONTAINS


      SUBROUTINE MASSADJ(AERO,GAS,SPCMASS1,SPCMASS2,EMIS_MASS,AQSO4RATE,TSTEP)
!----------------------------------------------------------------------------------------------------------------------
!     This routine rescales all aerosol and gas-phase species to enforce
!     mass conservation to machine precision.
!----------------------------------------------------------------------------------------------------------------------
      USE AERO_SETUP, ONLY: SULF_MAP, BCAR_MAP, OCAR_MAP, DUST_MAP, SEAS_MAP
      IMPLICIT NONE

      ! Arguments.
 
      REAL(8), INTENT(INOUT) :: AERO(NAEROBOX)         ! aerosol conc. [ug/m^3] or [#/m^3]
      REAL(8), INTENT(INOUT) :: GAS(NGASES)            ! gas-phase conc. [ug/m^3]
      REAL(8), INTENT(IN)    :: SPCMASS1(NMASS_SPCS+2) ! initial total mass spc. conc. [ug/m^3]
      REAL(8), INTENT(INOUT) :: SPCMASS2(NMASS_SPCS+2) ! final   total mass spc. conc. [ug/m^3]
      REAL(8), INTENT(IN)    :: EMIS_MASS(NEMIS_SPCS)  ! mass emission rates [ug/m^3/s]
      REAL(8), INTENT(IN)    :: AQSO4RATE              ! in-cloud SO4 production rate [ug/m^3/s]
      REAL(8), INTENT(IN)    :: TSTEP                  ! model physics time step [s]

      ! Local variables.

      INTEGER       :: I
      REAL(8)       :: SCALE(NMASS_SPCS+2)             ! scale factor for mass adjustment
      REAL(8), SAVE :: SCALEMAX = 1.0D-80
      REAL(8), SAVE :: SCALEMIN = 1.0D+80
      
      !----------------------------------------------------------------------------------------------------------------
      ! Get the precise mass conc. that should exist at the end of the time
      ! step, divided by the actual mass conc. at the end of the time step.
      !----------------------------------------------------------------------------------------------------------------
      SPCMASS2(:) = SPCMASS2(:) + TINYDENOM 
      SCALE(1) = ( SPCMASS1(1) + ( AQSO4RATE + EMIS_MASS(1) + EMIS_MASS(2)  ) * TSTEP ) / SPCMASS2(1) 
      SCALE(2) = ( SPCMASS1(2) + (             EMIS_MASS(3) + EMIS_MASS(8)  ) * TSTEP ) / SPCMASS2(2) 
      SCALE(3) = ( SPCMASS1(3) + (             EMIS_MASS(4) + EMIS_MASS(9)  ) * TSTEP ) / SPCMASS2(3) 
      SCALE(4) = ( SPCMASS1(4) + (             EMIS_MASS(5) + EMIS_MASS(10) ) * TSTEP ) / SPCMASS2(4) 
      SCALE(5) = ( SPCMASS1(5) + (             EMIS_MASS(6) + EMIS_MASS(7)  ) * TSTEP ) / SPCMASS2(5) 
      SCALE(6) = ( SPCMASS1(6)                                                        ) / SPCMASS2(6) 
      SCALE(7) = ( SPCMASS1(7)                                                        ) / SPCMASS2(7) 

      ! WRITE(*,'(7F14.9)') SCALE(:)
      ! WRITE(*,'(7E14.6)') SPCMASS1(6), SPCMASS2(6), SPCMASS1(7), SPCMASS2(7)
      !----------------------------------------------------------------------------------------------------------------

      AERO( SULF_MAP(:) ) = AERO( SULF_MAP(:) ) * SCALE(1)
      AERO( BCAR_MAP(:) ) = AERO( BCAR_MAP(:) ) * SCALE(2)
      AERO( OCAR_MAP(:) ) = AERO( OCAR_MAP(:) ) * SCALE(3)
      AERO( DUST_MAP(:) ) = AERO( DUST_MAP(:) ) * SCALE(4)
      AERO( SEAS_MAP(:) ) = AERO( SEAS_MAP(:) ) * SCALE(5)
      AERO( MASS_NO3    ) = AERO( MASS_NO3    ) * SCALE(6)
      AERO( MASS_NH4    ) = AERO( MASS_NH4    ) * SCALE(7)
      GAS ( GAS_H2SO4   ) = GAS ( GAS_H2SO4   ) * SCALE(1)
      GAS ( GAS_HNO3    ) = GAS ( GAS_HNO3    ) * SCALE(6)
      GAS ( GAS_NH3     ) = GAS ( GAS_NH3     ) * SCALE(7)
       
      
!----------------------------------------------------------------------------------------------------------------------
!     Track the maximum and minimum scale factors required.
!----------------------------------------------------------------------------------------------------------------------
      IF( WRITE_LOG ) THEN  
        WRITE(31,90000) SPCMASS1(:)
        WRITE(31,90000) SPCMASS2(:)
        WRITE(31,90000) SCALE(:)
        WRITE(31,*) '  '
        WRITE(32,90000) SCALE(:)
        DO I=1, NMASS_SPCS+2
          IF(     SCALE(I) .GT. SCALEMAX ) THEN
            SCALEMAX = SCALE(I)
            WRITE(33,90001) SCALE(I), SCALEMAX, SCALEMIN, I, SPCMASS1(I), SPCMASS2(I)
          ELSEIF( SCALE(I) .LT. SCALEMIN ) THEN
            SCALEMIN = SCALE(I)
            WRITE(33,90001) SCALE(I), SCALEMAX, SCALEMIN, I, SPCMASS1(I), SPCMASS2(I)
          ENDIF
        ENDDO
      ENDIF

90000 FORMAT(7D15.6)
90001 FORMAT(3D15.6,I6,2D15.6)
      RETURN
      END SUBROUTINE MASSADJ



      
      SUBROUTINE SIZE_PDFS(AERO,PDF1,PDF2)
      USE AERO_PARAM, ONLY: PI6, DENSP, IXXX, IYYY, ILAY
      USE AERO_CONFIG, ONLY: NMODES, NAEROBOX,NBINS
      USE AERO_SETUP, ONLY: SIG0, CONV_DPAM_TO_DGN, NUMB_MAP, MODE_NAME
      USE AERO_DIAM
      IMPLICIT NONE

      ! Arguments.
       REAL(8), INTENT(IN) :: AERO(NAEROBOX)! aerosol conc. [ug/m^3] or [#/m^3]

      ! Local variables. 

      INTEGER :: I, N
!      INTEGER, PARAMETER :: NBINS = 30! 200    ! number of bins [1]     defined in config 
      REAL(8) :: DGRID(NBINS)              ! fixed diameter        grid [um]
      REAL(8) :: MGRID(NBINS)              ! fixed mass/particle   grid [ug/particle]
      REAL(8) :: DLOWER(NBINS)             ! lower boundary fixed diameter grid [um]
      REAL(8) :: DUPPER(NBINS)             ! upper boundary fixed diameter grid [um]
      REAL(8) :: NTOT(NMODES)              ! number concentration for each mode [#/m^3]
      REAL(8) :: PDF(NBINS,2,NMODES)       ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF1(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: PDF2(NBINS)               ! number or mass conc. at each grid point [#/m^3] or [ug/m^3]       
      REAL(8) :: DNDLOGD(NMODES)           ! dN/dlog10(Dp) [ #/m^3]
      REAL(8) :: DMDLOGD(NMODES)           ! dM/dlog10(Dp) [ug/m^3]
      REAL(8) :: RDMIN                     ! reciprocal of DMIN to optimize coagulation [1/um]
      REAL(8) :: RDLOGDSC                  ! reciprocal of log10 of the grid spacing [1]
      REAL(8) :: SCALE, F, SUM1, SUM2      ! scratch variables 
      REAL(8) :: DMINL, DMAXL, DG          ! diameters [um]  
      REAL(8) :: FLN                       ! function for lognormal distribution [1]  
      REAL(8), PARAMETER :: DMIN =  0.001D+00   ! smallest particle diameter of the discrete grid [um]
      REAL(8), PARAMETER :: DMAX = 20.000D+00   ! largest  particle diameter of the discrete grid [um]


        DMAXL = DMAX
        DMINL = DMIN

      SCALE    = ( DMAXL / DMINL )**(1.0D+00/REAL(NBINS-1))
      RDLOGDSC = 1.0D+00 / LOG10( SCALE )
      RDMIN    = 1.0D+00 / DMINL
      DO I=1, NBINS
        DGRID(I)  = DMINL * SCALE**(I-1)                  ! [um]
        DLOWER(I) = DGRID(I) / SCALE**0.5D+00             ! [um]
        DUPPER(I) = DGRID(I) * SCALE**0.5D+00             ! [um]
        MGRID(I)  = 1.0D-06 * DENSP * PI6 * DGRID(I)**3   ! [ug/particle]
        DO N=1, NMODES
          DG = 1.0D+06 * DIAM(IXXX,IYYY,ILAY,N) * CONV_DPAM_TO_DGN(N)   ! convert [m] to [um] and Dbar to Dg
          NTOT(N) = AERO( NUMB_MAP(N) )
          F = NTOT(N) * FLN( DGRID(I), DG, SIG0(N) )
          PDF(I,1,N) = F * ( DUPPER(I) - DLOWER(I) )
          PDF(I,2,N) = PDF(I,1,N) * MGRID(I)
          DNDLOGD(N) = PDF(I,1,N) * RDLOGDSC * 1.0D-06           ! convert from [#/m^3] to [#/cm^3]
          DNDLOGD(N) = MAX( DNDLOGD(N), 1.0D-30 )
          DMDLOGD(N) = PDF(I,2,N) * RDLOGDSC                     ! [ug/m^3]
          DMDLOGD(N) = MAX( DMDLOGD(N), 1.0D-30 )
        ENDDO
c        WRITE(IUNIT,91) I, DGRID(I), DNDLOGD(:)
c        WRITE(JUNIT,91) I, DGRID(I), DMDLOGD(:)
      ENDDO

        PDF1(:) = 0.0D+00
        PDF2(:) = 0.0D+00
      DO N=1, NMODES
        DO I=1, NBINS
          PDF1(I) = PDF1(I)  + PDF(I,1,N)
          PDF2(I) = PDF2(I)  + PDF(I,2,N)
          SUM1 = SUM1 + PDF(I,1,N)
          SUM2 = SUM2 + PDF(I,2,N)
        ENDDO
      ENDDO

      RETURN 
      END SUBROUTINE SIZE_PDFS


      REAL(8) FUNCTION FLN(X,XG,SIGMAG)
      REAL(8) :: X      ! particle radius or diameter [any units]
      REAL(8) :: XG     ! geometric mean radius or diameter [any units]
      REAL(8) :: SIGMAG ! geometric standard deviation [monodisperse = 1.0]
      REAL(8), PARAMETER :: SQRTTWOPI = 2.506628275D+00
      FLN = EXP(-0.5D+00*(LOG(X/XG)/LOG(SIGMAG))**2) / (X*LOG(SIGMAG)*SQRTTWOPI)
      RETURN
      END FUNCTION FLN


      REAL(8) FUNCTION GETXNUM(NI,NJ,DGNI,DGNJ,XLSGI,XLSGJ)
!---------------------------------------------------------------------------------------------------------------------
!
! GETXNUM = ln( Dij / Dgi ) / ( sqrt(2) * ln(Sgi) ), where
!
!      Dij is the diameter of intersection,
!      Dgi is the median diameter of the smaller size mode, and
!      Sgi is the geometric standard deviation of smaller mode.
!
! A quadratic equation is solved to obtain GETXNUM, following the method of Press et al. 1992.
!  
! REFERENCES:
!
!  1. Binkowski, F.S. and S.J. Roselle, Models-3 Community Multiscale Air Quality (CMAQ) 
!     model aerosol component 1: Model Description.  J. Geophys. Res., Vol 108, No D6, 4183
!     doi:10.1029/2001JD001409, 2003.
!  2. Press, W.H., S.A. Teukolsky, W.T. Vetterling, and B.P. Flannery, Numerical Recipes in 
!     Fortran 77 - 2nd Edition. Cambridge University Press, 1992.
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE
      
      ! Arguments.
      
      REAL(8) :: NI         ! Aitken       mode number concentration [#/m^3] 
      REAL(8) :: NJ         ! accumulation mode number concentration [#/m^3] 
      REAL(8) :: DGNI       ! Aitken       mode geo. mean diameter [um] 
      REAL(8) :: DGNJ       ! accumulation mode geo. mean diameter [um]
      REAL(8) :: XLSGI      ! Aitken       mode ln(geo. std. dev.) [1]
      REAL(8) :: XLSGJ      ! accumulation mode ln(geo. std. dev.) [1]

      ! Local variables. 

      REAL(8) :: AA, BB, CC, DISC, QQ, ALFA, L, YJI
      REAL(8), PARAMETER :: SQRT2 = 1.414213562D+00

      ALFA = XLSGI / XLSGJ
      YJI = LOG( DGNJ / DGNI ) / ( SQRT2 * XLSGI )
      L = LOG( ALFA * NJ / NI)

      ! Calculate quadratic equation coefficients & discriminant.
      
      AA = 1.0D+00 - ALFA * ALFA
      BB = 2.0D+00 * YJI * ALFA * ALFA
      CC = L - YJI * YJI * ALFA * ALFA
      DISC = BB * BB - 4.0D+00 * AA * CC

      ! If roots are imaginary, return a negative GETAF value so that no IMTR takes place.
      
      IF( DISC .LT. 0.0D+00 ) THEN
        GETXNUM = - 5.0D+00         ! ERROR IN INTERSECTION
        RETURN
      ENDIF
      
      ! Equation 5.6.4 of Press et al. 1992.
      
      QQ = -0.5D+00 * ( BB + SIGN( 1.0D+00, BB ) * SQRT(DISC) )

      ! Return solution of the quadratic equation that corresponds to a
      ! diameter of intersection lying between the median diameters of the 2 modes.
      
      GETXNUM = CC / QQ       ! See Equation 5.6.5 of Press et al.
      
      ! WRITE(*,*)'GETXNUM = ', GETXNUM
      RETURN
      END FUNCTION GETXNUM

      END MODULE AERO_SUBS
 
