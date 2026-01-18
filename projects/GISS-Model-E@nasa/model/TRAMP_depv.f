      MODULE AERO_DEPV
      USE AERO_PARAM,  ONLY: NLAYS, IXXX, IYYY, ILAY
      USE AERO_CONFIG, ONLY: NMODES
      USE AMP_AEROSOL, ONLY: VDDEP_AERO
!-------------------------------------------------------------------------------------------------------------------------
!     The array VDDEP_AERO(X,Y,Z,I,1) contains current values for the dry deposition velocities 
!     for aerosol number concentrations for mode I. 
!     The array VDDEP_AERO(X,Y,Z,I,2) contains current values for the dry deposition velocities 
!     for aerosol mass   concentrations for mode I. 
!     Values in VDDEP_AERO are saved in subr. MATRIX at each time step. 
!-------------------------------------------------------------------------------------------------------------------------
      
      CONTAINS


      SUBROUTINE GET_AERO_DEPV(N,TK,RHOA,XLM,AMU,WSTAR,USTAR,RA,DGN_DDEP,XLS_DDEP,DEN_DDEP)
!----------------------------------------------------------------------------------------------------------------------
!     Calculate deposition velocity for Aitken, accumulation, and coarse modes.
!     Reference: Binkowski F. S., and U. Shankar, The regional particulate model 
!     1. Model description and preliminary results. J. Geophys. Res., 100, D12, 26191-26209, 1995.
!
!     12-18-06, DLW: Derived from the CMAQ routine aero_depv.f originally coded by F. S. Binkowski. :: 
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments. 

      INTEGER :: N                    ! number of modes [1]
      REAL(8) :: TK                   ! air temperature [K]
      REAL(8) :: RHOA                 ! air density  [kg/m^3]
      REAL(8) :: XLM                  ! atmospheric mean free path [m]
      REAL(8) :: AMU                  ! atmospheric dynamic viscosity [kg/(m s)]
      REAL(8) :: WSTAR                ! convective velocity scale [m/s]
      REAL(8) :: USTAR                ! friction velocity [m/s]
      REAL(8) :: RA                   ! aerodynamic resistance [s/m]
      REAL(8) :: DGN_DDEP(N)          ! geo. mean diameter    for each mode [um] 
      REAL(8) :: XLS_DDEP(N)          ! ln(geo. std. dev.)    for each mode [1]
      REAL(8) :: DEN_DDEP(N)          ! avg. particle density for each mode [g/cm^3]

      ! Local variables. 

      INTEGER :: I      
      REAL(8) :: DGN_M                ! geo. mean diameter    [m] 
      REAL(8) :: DEN_KGM3             ! avg. particle density [kg/m^3]
      REAL(8) :: VDEP(2)              ! deposition velocities [m/s]


      DO I=1, N
        DGN_M    = DGN_DDEP(I) * 1.0D-06    ! convert from [um] to [m]
        DEN_KGM3 = DEN_DDEP(I) * 1.0D+03    ! convert from [g/cm^3] to [kg/m^3]
        CALL GETDEP_V( TK, RHOA, XLM, AMU, WSTAR, USTAR, RA, DGN_M, XLS_DDEP(I), DEN_KGM3, VDEP )
!       VDEP(:) = MIN( VDEP(:), 10.0D+00 )  ! cap at 10 [m/s] = 1000 [cm/s]; should have no effect 
        VDDEP_AERO(IXXX,IYYY,I,1) = VDEP(1) ! for deposition of number [m/s]
        VDDEP_AERO(IXXX,IYYY,I,2) = VDEP(2) ! for deposition of mass   [m/s]
      ENDDO

      RETURN
      END SUBROUTINE GET_AERO_DEPV


      SUBROUTINE GETDEP_V( BLKTA, BLKDENS, XLM, AMU, BLKWSTAR, BLKUSTAR, BLKRA, DGACC, XXLSGAC, PDENSAC, VDEP )
!----------------------------------------------------------------------------------------------------------------------
!     Calculate deposition velocity for Aitken, accumulation, and coarse modes.
!     Reference: Binkowski F. S., and U. Shankar, The regional particulate model 
!     1. Model description and preliminary results. J. Geophys. Res., 100, D12, 26191-26209, 1995.
!
!     12-18-06, DLW: Derived from the CMAQ routine aero_depv.f originally coded by F. S. Binkowski. :: 
!----------------------------------------------------------------------------------------------------------------------
      IMPLICIT NONE

      ! Arguments. 

      REAL(8) :: BLKTA    ! air temperature [K]
      REAL(8) :: BLKDENS  ! air density  [kg/m^3]
      REAL(8) :: XLM      ! atmospheric mean free path [m]
      REAL(8) :: AMU      ! atmospheric dynamic viscosity [kg/(m s)]
      REAL(8) :: BLKWSTAR ! convective velocity scale [m/s]
      REAL(8) :: BLKUSTAR ! friction velocity [m/s]
      REAL(8) :: BLKRA    ! aerodynamic resistance [s/m]
      REAL(8) :: DGACC    ! geo. mean diamter [m]
      REAL(8) :: XXLSGAC  ! ln(geo. std. dev.) [1]
      REAL(8) :: PDENSAC  ! average particle density [kg/m^3]
      REAL(8) :: VDEP(2)  ! deposition  velocity [ m/s ]

      ! Local variables.

      REAL(8) :: KNACC    ! Modal Knudsen [1]
      REAL(8) :: DCHAT0A  ! Modal particle diffusivity for number [m^2/s]
      REAL(8) :: DCHAT3A  ! Modal particle diffusivity for mass   [m^2/s]
      REAL(8) :: VGHAT0A  ! Modal sedimentation velocity for number [m/s]
      REAL(8) :: VGHAT3A  ! Modal sedimentation velocity for mass   [m/s]
      REAL(8) :: DCONST1, DCONST1A
      REAL(8) :: DCONST2, DCONST3A
      REAL(8) :: SC0A     ! Schmidt numbers for number
      REAL(8) :: SC3A     ! Schmidt numbers for 3rd moment
      REAL(8) :: ST0A     ! Stokes numbers for number
      REAL(8) :: ST3A     ! Stokes numbers for 3rd moment
      REAL(8) :: RD0A     ! canopy resistance for number
      REAL(8) :: RD3A     ! canopy resisteance for 3rd moment
      REAL(8) :: UTSCALE  ! scratch function of USTAR and WSTAR
      REAL(8) :: NU       ! kinematic viscosity [ m**2 s**-1 ]
      REAL(8) :: USTFAC   ! scratch function of USTAR, NU, and GRAV

      REAL(8), PARAMETER :: BHAT     = 1.246D+00             ! Constant from Cunningham slip correction
      REAL(8), PARAMETER :: PI       = 3.141592653589793D+00        
      REAL(8), PARAMETER :: PI6      = PI / 6.0D+00
      REAL(8), PARAMETER :: THREEPI  = 3.0D+00 * PI
      REAL(8), PARAMETER :: ONE3     = 1.0D+00 / 3.0D+00
      REAL(8), PARAMETER :: TWO3     = 2.0D+00 / 3.0D+00
      REAL(8), PARAMETER :: AVO      = 6.0221367D+23         ! Avogadro's Constant  [1/mol]
      REAL(8), PARAMETER :: RGASUNIV = 8.314510D+00          ! universal gas const  [J/mol/K]
      REAL(8), PARAMETER :: BOLTZ    = RGASUNIV / AVO        ! Boltzmann's Constant [J/K]
      !----------------------------------------------------------------------------------------------------------------
      ! Value is the mean of polar and equatorial values. CRC Handbook (76th Ed) page 14-6. (FSB)
      !----------------------------------------------------------------------------------------------------------------
      REAL(8), PARAMETER :: GRAV        = 9.80622D+00        ! mean gravitational accel [m/s^2]
      REAL(8), PARAMETER :: DGACC_MAX   = 1.0D-06            ! 1.0 um, min. value for elimination of impaction term [m]
      REAL(8), PARAMETER :: XXLSGAC_MAX = 0.6931472D+00      ! ln(2),  min. value for elimination of impaction term [1]

      ! Scratch variables for standard deviations.

      REAL(8) :: L2SGAC           
      REAL(8) :: EAC1             
      REAL(8) :: ESAC04    
      REAL(8) :: ESAC08    
      REAL(8) :: ESAC16    
      REAL(8) :: ESAC20     
      REAL(8) :: ESAC28       
      REAL(8) :: ESAC32    
      REAL(8) :: ESAC64    


      KNACC  = 2.0D+00 * XLM / DGACC
      L2SGAC = XXLSGAC * XXLSGAC
      EAC1   = EXP( 0.125D+00 * L2SGAC )
      ESAC04 = EAC1**4
      ESAC08 = ESAC04 * ESAC04
      ESAC16 = ESAC08 * ESAC08
      ESAC20 = ESAC16 * ESAC04
      ESAC28 = ESAC20 * ESAC08
      ESAC32 = ESAC16 * ESAC16
      ESAC64 = ESAC32 * ESAC32

      DCONST1  = BOLTZ * BLKTA / ( THREEPI * AMU )
      DCONST1A = DCONST1 / DGACC
      DCONST2  = GRAV / ( 18.0D+00 * AMU )
      DCONST3A = DCONST2 * PDENSAC * DGACC*DGACC

      DCHAT0A = DCONST1A * ( ESAC04  + BHAT * KNACC * ESAC16  )
      DCHAT3A = DCONST1A * ( ( 1.0D+00 / ESAC20 ) + BHAT * KNACC / ESAC32 )
      VGHAT0A = DCONST3A * ( ESAC16  + BHAT * KNACC * ESAC04  )
      VGHAT3A = DCONST3A * ( ESAC64  + BHAT * KNACC * ESAC28  )

      NU      = AMU / BLKDENS
      USTFAC  = BLKUSTAR * BLKUSTAR / ( GRAV * NU )
      UTSCALE = BLKUSTAR + 0.24D+00 * BLKWSTAR * BLKWSTAR / BLKUSTAR

      SC0A = NU / DCHAT0A
      ST0A = MAX ( VGHAT0A * USTFAC, 0.01D+00 )
      IF( DGACC .LT. DGACC_MAX ) THEN                        ! Not a coarse mode ...
        RD0A = 1.0D+00 / ( UTSCALE * ( SC0A**( -TWO3 ) + 10.0**( -3.0D+00 / ST0A ) ) )
      ELSE
        RD0A = 1.0D+00 / ( UTSCALE * ( SC0A**( -TWO3 ) ) )   ! Eliminate impaction term for coarse modes as in CMAQ.
      ENDIF
      VDEP(1) = VGHAT0A + 1.0D+00 / ( BLKRA + RD0A + RD0A * BLKRA * VGHAT0A )  ! For deposition of number.

      SC3A = NU / DCHAT3A
      ST3A = MAX( VGHAT3A * USTFAC , 0.01D+00 )
      IF( DGACC .LT. DGACC_MAX  ) THEN                       ! Not a coarse mode ...
        RD3A = 1.0D+00 / ( UTSCALE * ( SC3A**( -TWO3 ) + 10.0**( -3.0D+00 / ST3A ) ) )
      ELSE
        RD3A = 1.0D+00 / ( UTSCALE * ( SC3A**( -TWO3 ) ) )   ! Eliminate impaction term for coarse modes as in CMAQ.
      ENDIF
      VDEP(2) = VGHAT3A + 1.0D+00 / ( BLKRA + RD3A + RD3A * BLKRA * VGHAT3A )  ! For deposition of mass.

      RETURN
      END SUBROUTINE GETDEP_V
      
      
      END MODULE AERO_DEPV
