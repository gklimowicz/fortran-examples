#include "rundeck_opts.h"

!@sum TOMAS_act : this subroutine computes clouds microphysical 
!@+    properties based on Nenes and Seinfeld (2003) and Fountoukis and 
!@+    Nenes (2005). 
!@+   NOTE from Yunha Lee : ModelE-style in-code documentation is only 
!@+    done for the main subroutine CALCND. For the other subroutines, 
!@+    original in-code comments from Georgia Tech are kept.   
          
C=======================================================================
C
!@sum   SUBROUTINE CALCNd
!@auth  ATHANASIOS NENES
!@+     GISS GCM (TOMAS microphys) interface
!@+     This subroutine calls the Nenes and Seinfeld (2003) sectional
!@+     activation parameterization, with the water vapor condensation
!@+     modifications of Fountoukis and Nenes (2005).
!
!@+   Modified by YUNHA LEE (JUNE 2011)
!@+   ALL variables defined with IMPLICIT are re-defined accordingly.
!     Common block is removed as well (TOMAS_ACTV module is used instead)
!
C
C=======================================================================
C
      MODULE TOMAS_ACTV
      
C======================================================================C
!@+ *** INCLUDE FILE 'PARAMETR.INC'
!@+ *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!@+     AND VARIABLES. 
!@+
!@auth  ATHANASIOS NENES
!@+     Modified by YUNHA LEE (JUNE 2011) 
!@+     PARAMETR.INC is moved to MODULE TOMAS_ACTV 
C
C======================================================================C
!YUNHA      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      INTEGER NSMX,NpGauss, n,NSEC,ISEC
      PARAMETER (NSMX = 2000, NpGauss=10, n=400)
      INTEGER ::  INTEGMOD,I,J
      REAL*8 :: NPART(NSMX+1)
      LOGICAL  :: CRIT1, CRIT2, CRIT3, DPMXINF
      INTEGER :: IDPMF,MAXIT
      INTEGER :: ISMIN,ISLW1, ISLW2,NITER

      REAL*8 :: WPARCEL, TEMP, PRES
     &    ,SC(NSMX), AKOH(NSMX), ACCOM
     &   ,SSRAT, DPMX, ALFA, BET1, BET2 
     &   ,AMW, AMA, DENW, CPAIR, DHV, AKA
     &   ,DV,  PSAT, DAIR, SURT
     &   ,DPNMX, WPDBG(NpGauss), PDDBG(NpGauss)
     &   ,NADBG(NpGauss), SMDBG(NpGauss)
     &   ,EPS, XGS(NpGauss), WGSi(NpGauss)
     &   ,ZERO, GREAT, SQ2PI
   

      DATA AMW     /18d-3/                   ! Water molecular weight
      DATA AMA     /29d-3/                   ! Air molecular weight
C
      DATA ACCOM /0.06/                     ! Default accommodation coef !it is different to one above -cyhl
C
      DATA MAXIT   /60/                      ! Max iterations for solution
      DATA EPS     /5d-4/                    ! Convergence criterion  ! CYHL-different to weichun
      DATA INTEGMOD/0/   ! Integral mode (0=Normal,1=Dp>>Dc, 2=Dp~Dc)
      DATA IDPMF   /0/   ! Kin.limit.fit (0=Avg fit,1=Singl-mode, 2=Tri-modal)

      DATA ZERO    /0d0/
      DATA GREAT   /1D30/
      DATA SQ2PI   /2.5066282746d0/

      REAL*8  GRV(26,26)
  
C
C *** Data For Row
C
      DATA GRV/
     &0.13E-19,0.79E-19,0.25E-18,0.60E-18,0.13E-17,0.27E-17,0.55E-17,
     &0.11E-16,0.22E-16,0.42E-16,0.81E-16,0.16E-15,0.31E-15,0.59E-15,
     &0.11E-14,0.22E-14,0.43E-14,0.84E-14,0.16E-13,0.32E-13,0.62E-13,
     &0.12E-12,0.24E-12,0.46E-12,0.88E-12,0.17E-11,0.79E-19,0.27E-19,
     &0.16E-18,0.50E-18,0.12E-17,0.26E-17,0.55E-17,0.11E-16,0.22E-16,
     &0.43E-16,0.84E-16,0.16E-15,0.32E-15,0.61E-15,0.12E-14,0.23E-14,
     &0.44E-14,0.86E-14,0.17E-13,0.33E-13,0.63E-13,0.12E-12,0.24E-12,
     &0.46E-12,0.89E-12,0.17E-11,0.25E-18,0.16E-18,0.54E-19,0.32E-18,
     &0.99E-18,0.24E-17,0.53E-17,0.11E-16,0.22E-16,0.44E-16,0.86E-16,
     &0.17E-15,0.32E-15,0.63E-15,0.12E-14,0.24E-14,0.46E-14,0.88E-14,
     &0.17E-13,0.33E-13,0.65E-13,0.13E-12,0.24E-12,0.47E-12,0.91E-12,
     &0.17E-11,0.60E-18,0.50E-18,0.32E-18,0.11E-18,0.63E-18,0.62E-17,
     &0.18E-16,0.44E-16,0.95E-16,0.20E-15,0.40E-15,0.78E-15,0.15E-14,
     &0.30E-14,0.58E-14,0.11E-13,0.22E-13,0.42E-13,0.81E-13,0.16E-12,
     &0.30E-12,0.59E-12,0.11E-11,0.22E-11,0.42E-11,0.81E-11,0.13E-17,
     &0.12E-17,0.99E-18,0.63E-18,0.44E-17,0.32E-16,0.11E-15,0.28E-15,
     &0.63E-15,0.13E-14,0.27E-14,0.53E-14,0.11E-13,0.20E-13,0.40E-13,
     &0.77E-13,0.15E-12,0.29E-12,0.56E-12,0.11E-11,0.21E-11,0.40E-11,
     &0.78E-11,0.15E-10,0.29E-10,0.55E-10,0.27E-17,0.26E-17,0.24E-17,
     &0.62E-17,0.32E-16,0.21E-16,0.14E-15,0.46E-15,0.11E-14,0.25E-14,
     &0.52E-14,0.10E-13,0.21E-13,0.41E-13,0.79E-13,0.15E-12,0.30E-12,
     &0.58E-12,0.11E-11,0.22E-11,0.42E-11,0.80E-11,0.15E-10,0.30E-10,
     &0.57E-10,0.11E-09,0.55E-17,0.55E-17,0.53E-17,0.18E-16,0.11E-15,
     &0.14E-15,0.66E-16,0.41E-15,0.13E-14,0.33E-14,0.72E-14,0.15E-13,
     &0.30E-13,0.60E-13,0.12E-12,0.23E-12,0.44E-12,0.85E-12,0.17E-11,
     &0.32E-11,0.61E-11,0.12E-10,0.23E-10,0.44E-10,0.84E-10,0.16E-09,
     &0.11E-16,0.11E-16,0.11E-16,0.44E-16,0.28E-15,0.46E-15,0.41E-15,
     &0.17E-15,0.11E-14,0.34E-14,0.83E-14,0.18E-13,0.38E-13,0.77E-13,
     &0.15E-12,0.30E-12,0.58E-12,0.11E-11,0.22E-11,0.42E-11,0.81E-11,
     &0.15E-10,0.30E-10,0.57E-10,0.11E-09,0.21E-09,0.22E-16,0.22E-16,
     &0.22E-16,0.95E-16,0.63E-15,0.11E-14,0.13E-14,0.11E-14,0.42E-15,
     &0.25E-14,0.81E-14,0.20E-13,0.44E-13,0.90E-13,0.18E-12,0.36E-12,
     &0.71E-12,0.14E-11,0.27E-11,0.51E-11,0.99E-11,0.19E-10,0.36E-10,
     &0.70E-10,0.13E-09,0.25E-09,0.42E-16,0.43E-16,0.44E-16,0.20E-15,
     &0.13E-14,0.25E-14,0.33E-14,0.34E-14,0.25E-14,0.97E-15,0.59E-14,
     &0.19E-13,0.45E-13,0.99E-13,0.21E-12,0.41E-12,0.82E-12,0.16E-11,
     &0.31E-11,0.60E-11,0.12E-10,0.22E-10,0.43E-10,0.82E-10,0.16E-09,
     &0.30E-09,0.81E-16,0.84E-16,0.86E-16,0.40E-15,0.27E-14,0.52E-14,
     &0.72E-14,0.83E-14,0.81E-14,0.59E-14,0.22E-14,0.13E-13,0.41E-13,
     &0.10E-12,0.22E-12,0.46E-12,0.92E-12,0.18E-11,0.36E-11,0.69E-11,
     &0.13E-10,0.26E-10,0.49E-10,0.94E-10,0.18E-09,0.34E-09,0.16E-15,
     &0.16E-15,0.17E-15,0.78E-15,0.53E-14,0.10E-13,0.15E-13,0.18E-13,
     &0.20E-13,0.19E-13,0.13E-13,0.48E-14,0.29E-13,0.90E-13,0.22E-12,
     &0.48E-12,0.99E-12,0.20E-11,0.39E-11,0.77E-11,0.15E-10,0.29E-10,
     &0.55E-10,0.11E-09,0.20E-09,0.38E-09,0.31E-15,0.32E-15,0.32E-15,
     &0.15E-14,0.11E-13,0.21E-13,0.30E-13,0.38E-13,0.44E-13,0.45E-13,
     &0.41E-13,0.29E-13,0.10E-13,0.61E-13,0.19E-12,0.47E-12,0.10E-11,
     &0.21E-11,0.43E-11,0.84E-11,0.16E-10,0.32E-10,0.61E-10,0.12E-09,
     &0.22E-09,0.42E-09,0.59E-15,0.61E-15,0.63E-15,0.30E-14,0.20E-13,
     &0.41E-13,0.60E-13,0.77E-13,0.90E-13,0.99E-13,0.10E-12,0.90E-13,
     &0.61E-13,0.22E-13,0.13E-12,0.41E-12,0.98E-12,0.22E-11,0.45E-11,
     &0.90E-11,0.18E-10,0.34E-10,0.66E-10,0.13E-09,0.24E-09,0.46E-09,
     &0.11E-14,0.12E-14,0.12E-14,0.58E-14,0.40E-13,0.79E-13,0.12E-12,
     &0.15E-12,0.18E-12,0.21E-12,0.22E-12,0.22E-12,0.19E-12,0.13E-12,
     &0.45E-13,0.27E-12,0.85E-12,0.21E-11,0.45E-11,0.93E-11,0.19E-10,
     &0.37E-10,0.71E-10,0.14E-09,0.26E-09,0.49E-09,0.22E-14,0.23E-14,
     &0.24E-14,0.11E-13,0.77E-13,0.15E-12,0.23E-12,0.30E-12,0.36E-12,
     &0.41E-12,0.46E-12,0.48E-12,0.47E-12,0.41E-12,0.27E-12,0.95E-13,
     &0.56E-12,0.18E-11,0.43E-11,0.93E-11,0.19E-10,0.39E-10,0.76E-10,
     &0.15E-09,0.28E-09,0.53E-09,0.43E-14,0.44E-14,0.46E-14,0.22E-13,
     &0.15E-12,0.30E-12,0.44E-12,0.58E-12,0.71E-12,0.82E-12,0.92E-12,
     &0.99E-12,0.10E-11,0.98E-12,0.85E-12,0.56E-12,0.20E-12,0.12E-11,
     &0.36E-11,0.88E-11,0.19E-10,0.40E-10,0.79E-10,0.15E-09,0.30E-09,
     &0.56E-09,0.84E-14,0.86E-14,0.88E-14,0.42E-13,0.29E-12,0.58E-12,
     &0.85E-12,0.11E-11,0.14E-11,0.16E-11,0.18E-11,0.20E-11,0.21E-11,
     &0.22E-11,0.21E-11,0.18E-11,0.12E-11,0.40E-12,0.24E-11,0.75E-11,
     &0.18E-10,0.39E-10,0.80E-10,0.16E-09,0.31E-09,0.59E-09,0.16E-13,
     &0.17E-13,0.17E-13,0.81E-13,0.56E-12,0.11E-11,0.17E-11,0.22E-11,
     &0.27E-11,0.31E-11,0.36E-11,0.39E-11,0.43E-11,0.45E-11,0.45E-11,
     &0.43E-11,0.36E-11,0.24E-11,0.82E-12,0.49E-11,0.15E-10,0.37E-10,
     &0.79E-10,0.16E-09,0.32E-09,0.62E-09,0.32E-13,0.33E-13,0.33E-13,
     &0.16E-12,0.11E-11,0.22E-11,0.32E-11,0.42E-11,0.51E-11,0.60E-11,
     &0.69E-11,0.77E-11,0.84E-11,0.90E-11,0.93E-11,0.93E-11,0.88E-11,
     &0.75E-11,0.49E-11,0.17E-11,0.98E-11,0.31E-10,0.73E-10,0.16E-09,
     &0.32E-09,0.63E-09,0.62E-13,0.63E-13,0.65E-13,0.30E-12,0.21E-11,
     &0.42E-11,0.61E-11,0.81E-11,0.99E-11,0.12E-10,0.13E-10,0.15E-10,
     &0.16E-10,0.18E-10,0.19E-10,0.19E-10,0.19E-10,0.18E-10,0.15E-10,
     &0.98E-11,0.34E-11,0.20E-10,0.61E-10,0.15E-09,0.31E-09,0.63E-09,
     &0.12E-12,0.12E-12,0.13E-12,0.59E-12,0.40E-11,0.80E-11,0.12E-10,
     &0.15E-10,0.19E-10,0.22E-10,0.26E-10,0.29E-10,0.32E-10,0.34E-10,
     &0.37E-10,0.39E-10,0.40E-10,0.39E-10,0.37E-10,0.31E-10,0.20E-10,
     &0.67E-11,0.39E-10,0.12E-09,0.29E-09,0.61E-09,0.24E-12,0.24E-12,
     &0.24E-12,0.11E-11,0.78E-11,0.15E-10,0.23E-10,0.30E-10,0.36E-10,
     &0.43E-10,0.49E-10,0.55E-10,0.61E-10,0.66E-10,0.71E-10,0.76E-10,
     &0.79E-10,0.80E-10,0.79E-10,0.73E-10,0.61E-10,0.39E-10,0.13E-10,
     &0.77E-10,0.24E-09,0.56E-09,0.46E-12,0.46E-12,0.47E-12,0.22E-11,
     &0.15E-10,0.30E-10,0.44E-10,0.57E-10,0.70E-10,0.82E-10,0.94E-10,
     &0.11E-09,0.12E-09,0.13E-09,0.14E-09,0.15E-09,0.15E-09,0.16E-09,
     &0.16E-09,0.16E-09,0.15E-09,0.12E-09,0.77E-10,0.26E-10,0.15E-09,
     &0.46E-09,0.88E-12,0.89E-12,0.91E-12,0.42E-11,0.29E-10,0.57E-10,
     &0.84E-10,0.11E-09,0.13E-09,0.16E-09,0.18E-09,0.20E-09,0.22E-09,
     &0.24E-09,0.26E-09,0.28E-09,0.30E-09,0.31E-09,0.32E-09,0.32E-09,
     &0.31E-09,0.29E-09,0.24E-09,0.15E-09,0.50E-10,0.29E-09,0.17E-11,
     &0.17E-11,0.17E-11,0.81E-11,0.55E-10,0.11E-09,0.16E-09,0.21E-09,
     &0.25E-09,0.30E-09,0.34E-09,0.38E-09,0.42E-09,0.46E-09,0.49E-09,
     &0.53E-09,0.56E-09,0.59E-09,0.62E-09,0.63E-09,0.63E-09,0.61E-09,
     &0.56E-09,0.46E-09,0.29E-09,0.97E-10
     & /


      REAL*8 TUROCEAN(26,26)
  
C
C *** Data For Row
C
      DATA TUROCEAN/
     &0.64E-16,0.34E-15,0.84E-15,0.16E-14,0.28E-14,0.46E-14,0.79E-14,
     &0.97E-14,0.14E-13,0.20E-13,0.34E-13,0.46E-13,0.50E-13,0.68E-13,
     &0.96E-13,0.13E-12,0.19E-12,0.29E-12,0.43E-12,0.62E-12,0.89E-12,
     &0.13E-11,0.18E-11,0.27E-11,0.72E-11,0.12E-10,0.34E-15,0.96E-16,
     &0.54E-15,0.13E-14,0.26E-14,0.45E-14,0.72E-14,0.13E-13,0.15E-13,
     &0.22E-13,0.32E-13,0.55E-13,0.72E-13,0.79E-13,0.11E-12,0.16E-12,
     &0.22E-12,0.33E-12,0.50E-12,0.73E-12,0.10E-11,0.15E-11,0.21E-11,
     &0.31E-11,0.94E-11,0.18E-10,0.84E-15,0.54E-15,0.14E-15,0.86E-15,
     &0.21E-14,0.41E-14,0.72E-14,0.12E-13,0.20E-13,0.25E-13,0.34E-13,
     &0.50E-13,0.84E-13,0.11E-12,0.13E-12,0.18E-12,0.25E-12,0.37E-12,
     &0.56E-12,0.86E-12,0.12E-11,0.16E-11,0.24E-11,0.36E-11,0.12E-10,
     &0.24E-10,0.16E-14,0.13E-14,0.86E-15,0.22E-15,0.14E-14,0.33E-14,
     &0.65E-14,0.11E-13,0.18E-13,0.32E-13,0.39E-13,0.55E-13,0.80E-13,
     &0.12E-12,0.16E-12,0.21E-12,0.30E-12,0.42E-12,0.62E-12,0.97E-12,
     &0.13E-11,0.17E-11,0.26E-11,0.41E-11,0.15E-10,0.32E-10,0.28E-14,
     &0.26E-14,0.21E-14,0.14E-14,0.33E-15,0.22E-14,0.53E-14,0.10E-13,
     &0.18E-13,0.29E-13,0.51E-13,0.62E-13,0.87E-13,0.13E-12,0.19E-12,
     &0.26E-12,0.35E-12,0.49E-12,0.68E-12,0.10E-11,0.15E-11,0.18E-11,
     &0.29E-11,0.54E-11,0.19E-10,0.41E-10,0.46E-14,0.45E-14,0.41E-14,
     &0.33E-14,0.22E-14,0.51E-15,0.34E-14,0.85E-14,0.17E-13,0.29E-13,
     &0.46E-13,0.81E-13,0.99E-13,0.14E-12,0.21E-12,0.31E-12,0.41E-12,
     &0.59E-12,0.81E-12,0.11E-11,0.18E-11,0.22E-11,0.33E-11,0.69E-11,
     &0.30E-10,0.52E-10,0.79E-14,0.72E-14,0.72E-14,0.65E-14,0.53E-14,
     &0.34E-14,0.78E-15,0.55E-14,0.13E-13,0.26E-13,0.46E-13,0.74E-13,
     &0.12E-12,0.16E-12,0.23E-12,0.37E-12,0.50E-12,0.65E-12,0.97E-12,
     &0.13E-11,0.21E-11,0.33E-11,0.45E-11,0.86E-11,0.44E-10,0.86E-10,
     &0.97E-14,0.13E-13,0.12E-13,0.11E-13,0.10E-13,0.85E-14,0.55E-14,
     &0.12E-14,0.87E-14,0.21E-13,0.42E-13,0.73E-13,0.12E-12,0.18E-12,
     &0.26E-12,0.42E-12,0.66E-12,0.83E-12,0.11E-11,0.16E-11,0.28E-11,
     &0.46E-11,0.74E-11,0.12E-10,0.61E-10,0.13E-09,0.14E-13,0.15E-13,
     &0.20E-13,0.18E-13,0.18E-13,0.17E-13,0.13E-13,0.87E-14,0.18E-14,
     &0.14E-13,0.34E-13,0.67E-13,0.12E-12,0.19E-12,0.28E-12,0.48E-12,
     &0.74E-12,0.12E-11,0.15E-11,0.20E-11,0.37E-11,0.68E-11,0.11E-10,
     &0.20E-10,0.83E-10,0.18E-09,0.20E-13,0.22E-13,0.25E-13,0.32E-13,
     &0.29E-13,0.29E-13,0.26E-13,0.21E-13,0.14E-13,0.29E-14,0.22E-13,
     &0.55E-13,0.11E-12,0.19E-12,0.32E-12,0.47E-12,0.86E-12,0.13E-11,
     &0.20E-11,0.27E-11,0.44E-11,0.94E-11,0.16E-10,0.28E-10,0.11E-09,
     &0.25E-09,0.34E-13,0.32E-13,0.34E-13,0.39E-13,0.51E-13,0.46E-13,
     &0.46E-13,0.42E-13,0.34E-13,0.22E-13,0.45E-14,0.35E-13,0.88E-13,
     &0.17E-12,0.32E-12,0.57E-12,0.80E-12,0.15E-11,0.23E-11,0.34E-11,
     &0.58E-11,0.10E-10,0.21E-10,0.38E-10,0.15E-09,0.32E-09,0.46E-13,
     &0.55E-13,0.50E-13,0.55E-13,0.62E-13,0.81E-13,0.74E-13,0.73E-13,
     &0.67E-13,0.55E-13,0.35E-13,0.70E-14,0.56E-13,0.14E-12,0.28E-12,
     &0.59E-12,0.10E-11,0.14E-11,0.28E-11,0.41E-11,0.84E-11,0.14E-10,
     &0.29E-10,0.50E-10,0.18E-09,0.41E-09,0.50E-13,0.72E-13,0.84E-13,
     &0.80E-13,0.87E-13,0.99E-13,0.12E-12,0.12E-12,0.12E-12,0.11E-12,
     &0.88E-13,0.56E-13,0.11E-13,0.86E-13,0.23E-12,0.47E-12,0.11E-11,
     &0.19E-11,0.27E-11,0.50E-11,0.11E-10,0.23E-10,0.41E-10,0.78E-10,
     &0.21E-09,0.48E-09,0.68E-13,0.79E-13,0.11E-12,0.12E-12,0.13E-12,
     &0.14E-12,0.16E-12,0.18E-12,0.19E-12,0.19E-12,0.17E-12,0.14E-12,
     &0.86E-13,0.17E-13,0.13E-12,0.40E-12,0.80E-12,0.20E-11,0.35E-11,
     &0.54E-11,0.14E-10,0.34E-10,0.59E-10,0.11E-09,0.25E-09,0.56E-09,
     &0.96E-13,0.11E-12,0.13E-12,0.16E-12,0.19E-12,0.21E-12,0.23E-12,
     &0.26E-12,0.28E-12,0.32E-12,0.32E-12,0.28E-12,0.23E-12,0.13E-12,
     &0.27E-13,0.22E-12,0.69E-12,0.14E-11,0.35E-11,0.64E-11,0.15E-10,
     &0.43E-10,0.86E-10,0.14E-09,0.30E-09,0.63E-09,0.13E-12,0.16E-12,
     &0.18E-12,0.21E-12,0.26E-12,0.31E-12,0.37E-12,0.42E-12,0.48E-12,
     &0.47E-12,0.57E-12,0.59E-12,0.47E-12,0.40E-12,0.22E-12,0.43E-13,
     &0.35E-12,0.12E-11,0.27E-11,0.63E-11,0.18E-10,0.44E-10,0.11E-09,
     &0.20E-09,0.36E-09,0.70E-09,0.19E-12,0.22E-12,0.25E-12,0.30E-12,
     &0.35E-12,0.41E-12,0.50E-12,0.66E-12,0.74E-12,0.86E-12,0.80E-12,
     &0.10E-11,0.11E-11,0.80E-12,0.69E-12,0.35E-12,0.68E-13,0.57E-12,
     &0.21E-11,0.52E-11,0.17E-10,0.51E-10,0.12E-09,0.26E-09,0.43E-09,
     &0.78E-09,0.29E-12,0.33E-12,0.37E-12,0.42E-12,0.49E-12,0.59E-12,
     &0.65E-12,0.83E-12,0.12E-11,0.13E-11,0.15E-11,0.14E-11,0.19E-11,
     &0.20E-11,0.14E-11,0.12E-11,0.57E-12,0.11E-12,0.98E-12,0.36E-11,
     &0.14E-10,0.50E-10,0.13E-09,0.29E-09,0.50E-09,0.87E-09,0.43E-12,
     &0.50E-12,0.56E-12,0.62E-12,0.68E-12,0.81E-12,0.97E-12,0.11E-11,
     &0.15E-11,0.20E-11,0.23E-11,0.28E-11,0.27E-11,0.35E-11,0.35E-11,
     &0.27E-11,0.21E-11,0.98E-12,0.18E-12,0.17E-11,0.11E-10,0.41E-10,
     &0.13E-09,0.31E-09,0.54E-09,0.93E-09,0.62E-12,0.73E-12,0.86E-12,
     &0.97E-12,0.10E-11,0.11E-11,0.13E-11,0.16E-11,0.20E-11,0.27E-11,
     &0.34E-11,0.41E-11,0.50E-11,0.54E-11,0.64E-11,0.63E-11,0.52E-11,
     &0.36E-11,0.17E-11,0.29E-12,0.46E-11,0.33E-10,0.11E-09,0.30E-09,
     &0.55E-09,0.96E-09,0.89E-12,0.10E-11,0.12E-11,0.13E-11,0.15E-11,
     &0.18E-11,0.21E-11,0.28E-11,0.37E-11,0.44E-11,0.58E-11,0.84E-11,
     &0.11E-10,0.14E-10,0.15E-10,0.18E-10,0.17E-10,0.14E-10,0.11E-10,
     &0.46E-11,0.80E-12,0.13E-10,0.85E-10,0.25E-09,0.52E-09,0.96E-09,
     &0.13E-11,0.15E-11,0.16E-11,0.17E-11,0.18E-11,0.22E-11,0.33E-11,
     &0.46E-11,0.68E-11,0.94E-11,0.10E-10,0.14E-10,0.23E-10,0.34E-10,
     &0.43E-10,0.44E-10,0.51E-10,0.50E-10,0.41E-10,0.33E-10,0.13E-10,
     &0.23E-11,0.43E-10,0.20E-09,0.46E-09,0.90E-09,0.18E-11,0.21E-11,
     &0.24E-11,0.26E-11,0.29E-11,0.33E-11,0.45E-11,0.74E-11,0.11E-10,
     &0.16E-10,0.21E-10,0.29E-10,0.41E-10,0.59E-10,0.86E-10,0.11E-09,
     &0.12E-09,0.13E-09,0.13E-09,0.11E-09,0.85E-10,0.43E-10,0.94E-11,
     &0.12E-09,0.37E-09,0.81E-09,0.27E-11,0.31E-11,0.36E-11,0.41E-11,
     &0.54E-11,0.69E-11,0.86E-11,0.12E-10,0.20E-10,0.28E-10,0.38E-10,
     &0.50E-10,0.78E-10,0.11E-09,0.14E-09,0.20E-09,0.26E-09,0.29E-09,
     &0.31E-09,0.30E-09,0.25E-09,0.20E-09,0.12E-09,0.29E-10,0.24E-09,
     &0.68E-09,0.72E-11,0.94E-11,0.12E-10,0.15E-10,0.19E-10,0.30E-10,
     &0.44E-10,0.61E-10,0.83E-10,0.11E-09,0.15E-09,0.18E-09,0.21E-09,
     &0.25E-09,0.30E-09,0.36E-09,0.43E-09,0.50E-09,0.54E-09,0.55E-09,
     &0.52E-09,0.46E-09,0.37E-09,0.24E-09,0.64E-10,0.46E-09,0.12E-10,
     &0.18E-10,0.24E-10,0.32E-10,0.41E-10,0.52E-10,0.86E-10,0.13E-09,
     &0.18E-09,0.25E-09,0.32E-09,0.41E-09,0.48E-09,0.56E-09,0.63E-09,
     &0.70E-09,0.78E-09,0.87E-09,0.93E-09,0.96E-09,0.96E-09,0.90E-09,
     &0.81E-09,0.68E-09,0.46E-09,0.13E-09
     & /

      REAL*8  TUREARTH(26,26)
  
C
C *** Data For Row
C
      DATA TUREARTH/
     &0.14E-15,0.43E-15,0.97E-15,0.19E-14,0.32E-14,0.52E-14,0.91E-14,
     &0.11E-13,0.16E-13,0.24E-13,0.41E-13,0.57E-13,0.63E-13,0.87E-13,
     &0.12E-12,0.18E-12,0.26E-12,0.40E-12,0.61E-12,0.91E-12,0.13E-11,
     &0.20E-11,0.28E-11,0.40E-11,0.11E-10,0.17E-10,0.43E-15,0.20E-15,
     &0.67E-15,0.15E-14,0.30E-14,0.52E-14,0.84E-14,0.15E-13,0.18E-13,
     &0.26E-13,0.38E-13,0.67E-13,0.90E-13,0.10E-12,0.14E-12,0.21E-12,
     &0.30E-12,0.46E-12,0.71E-12,0.11E-11,0.15E-11,0.22E-11,0.32E-11,
     &0.47E-11,0.14E-10,0.25E-10,0.97E-15,0.67E-15,0.30E-15,0.11E-14,
     &0.24E-14,0.47E-14,0.83E-14,0.13E-13,0.24E-13,0.29E-13,0.42E-13,
     &0.62E-13,0.11E-12,0.14E-12,0.17E-12,0.24E-12,0.35E-12,0.51E-12,
     &0.80E-12,0.13E-11,0.17E-11,0.24E-11,0.36E-11,0.54E-11,0.18E-10,
     &0.35E-10,0.19E-14,0.15E-14,0.11E-14,0.44E-15,0.17E-14,0.39E-14,
     &0.76E-14,0.13E-13,0.22E-13,0.38E-13,0.47E-13,0.67E-13,0.10E-12,
     &0.16E-12,0.21E-12,0.28E-12,0.41E-12,0.58E-12,0.88E-12,0.14E-11,
     &0.20E-11,0.25E-11,0.40E-11,0.62E-11,0.22E-10,0.46E-10,0.32E-14,
     &0.30E-14,0.24E-14,0.17E-14,0.65E-15,0.26E-14,0.63E-14,0.12E-13,
     &0.22E-13,0.35E-13,0.62E-13,0.77E-13,0.11E-12,0.16E-12,0.25E-12,
     &0.34E-12,0.48E-12,0.69E-12,0.98E-12,0.15E-11,0.23E-11,0.27E-11,
     &0.44E-11,0.80E-11,0.28E-10,0.58E-10,0.52E-14,0.52E-14,0.47E-14,
     &0.39E-14,0.26E-14,0.97E-15,0.42E-14,0.10E-13,0.20E-13,0.35E-13,
     &0.57E-13,0.10E-12,0.12E-12,0.18E-12,0.28E-12,0.41E-12,0.56E-12,
     &0.82E-12,0.12E-11,0.16E-11,0.27E-11,0.34E-11,0.49E-11,0.10E-10,
     &0.44E-10,0.75E-10,0.91E-14,0.84E-14,0.83E-14,0.76E-14,0.63E-14,
     &0.42E-14,0.15E-14,0.67E-14,0.16E-13,0.32E-13,0.56E-13,0.92E-13,
     &0.16E-12,0.20E-12,0.31E-12,0.50E-12,0.69E-12,0.92E-12,0.14E-11,
     &0.20E-11,0.32E-11,0.50E-11,0.69E-11,0.13E-10,0.65E-10,0.12E-09,
     &0.11E-13,0.15E-13,0.13E-13,0.13E-13,0.12E-13,0.10E-13,0.67E-14,
     &0.22E-14,0.11E-13,0.26E-13,0.52E-13,0.91E-13,0.15E-12,0.24E-12,
     &0.35E-12,0.56E-12,0.90E-12,0.12E-11,0.17E-11,0.24E-11,0.42E-11,
     &0.70E-11,0.11E-10,0.19E-10,0.90E-10,0.19E-09,0.16E-13,0.18E-13,
     &0.24E-13,0.22E-13,0.22E-13,0.20E-13,0.16E-13,0.11E-13,0.34E-14,
     &0.17E-13,0.43E-13,0.84E-13,0.15E-12,0.24E-12,0.37E-12,0.64E-12,
     &0.10E-11,0.16E-11,0.21E-11,0.30E-11,0.55E-11,0.10E-10,0.17E-10,
     &0.30E-10,0.12E-09,0.26E-09,0.24E-13,0.26E-13,0.29E-13,0.38E-13,
     &0.35E-13,0.35E-13,0.32E-13,0.26E-13,0.17E-13,0.52E-14,0.28E-13,
     &0.70E-13,0.14E-12,0.24E-12,0.42E-12,0.64E-12,0.12E-11,0.19E-11,
     &0.29E-11,0.39E-11,0.66E-11,0.14E-10,0.24E-10,0.43E-10,0.17E-09,
     &0.36E-09,0.41E-13,0.38E-13,0.42E-13,0.47E-13,0.62E-13,0.57E-13,
     &0.56E-13,0.52E-13,0.43E-13,0.28E-13,0.82E-14,0.46E-13,0.11E-12,
     &0.23E-12,0.43E-12,0.79E-12,0.11E-11,0.22E-11,0.34E-11,0.51E-11,
     &0.87E-11,0.16E-10,0.33E-10,0.58E-10,0.22E-09,0.47E-09,0.57E-13,
     &0.67E-13,0.62E-13,0.67E-13,0.77E-13,0.10E-12,0.92E-13,0.91E-13,
     &0.84E-13,0.70E-13,0.46E-13,0.13E-13,0.75E-13,0.19E-12,0.38E-12,
     &0.81E-12,0.15E-11,0.19E-11,0.41E-11,0.62E-11,0.13E-10,0.21E-10,
     &0.44E-10,0.76E-10,0.26E-09,0.60E-09,0.63E-13,0.90E-13,0.11E-12,
     &0.10E-12,0.11E-12,0.12E-12,0.16E-12,0.15E-12,0.15E-12,0.14E-12,
     &0.11E-12,0.75E-13,0.21E-13,0.12E-12,0.32E-12,0.66E-12,0.15E-11,
     &0.27E-11,0.40E-11,0.76E-11,0.17E-10,0.35E-10,0.63E-10,0.12E-09,
     &0.31E-09,0.71E-09,0.87E-13,0.10E-12,0.14E-12,0.16E-12,0.16E-12,
     &0.18E-12,0.20E-12,0.24E-12,0.24E-12,0.24E-12,0.23E-12,0.19E-12,
     &0.12E-12,0.33E-13,0.19E-12,0.57E-12,0.12E-11,0.29E-11,0.52E-11,
     &0.82E-11,0.22E-10,0.52E-10,0.92E-10,0.17E-09,0.38E-09,0.82E-09,
     &0.12E-12,0.14E-12,0.17E-12,0.21E-12,0.25E-12,0.28E-12,0.31E-12,
     &0.35E-12,0.37E-12,0.42E-12,0.43E-12,0.38E-12,0.32E-12,0.19E-12,
     &0.54E-13,0.33E-12,0.10E-11,0.21E-11,0.54E-11,0.99E-11,0.24E-10,
     &0.67E-10,0.14E-09,0.22E-09,0.46E-09,0.92E-09,0.18E-12,0.21E-12,
     &0.24E-12,0.28E-12,0.34E-12,0.41E-12,0.50E-12,0.56E-12,0.64E-12,
     &0.64E-12,0.79E-12,0.81E-12,0.66E-12,0.57E-12,0.33E-12,0.90E-13,
     &0.56E-12,0.19E-11,0.42E-11,0.10E-10,0.28E-10,0.70E-10,0.17E-09,
     &0.32E-09,0.54E-09,0.10E-08,0.26E-12,0.30E-12,0.35E-12,0.41E-12,
     &0.48E-12,0.56E-12,0.69E-12,0.90E-12,0.10E-11,0.12E-11,0.11E-11,
     &0.15E-11,0.15E-11,0.12E-11,0.10E-11,0.56E-12,0.15E-12,0.96E-12,
     &0.34E-11,0.85E-11,0.28E-10,0.83E-10,0.19E-09,0.41E-09,0.67E-09,
     &0.12E-08,0.40E-12,0.46E-12,0.51E-12,0.58E-12,0.69E-12,0.82E-12,
     &0.92E-12,0.12E-11,0.16E-11,0.19E-11,0.22E-11,0.19E-11,0.27E-11,
     &0.29E-11,0.21E-11,0.19E-11,0.96E-12,0.26E-12,0.18E-11,0.62E-11,
     &0.24E-10,0.84E-10,0.22E-09,0.47E-09,0.78E-09,0.13E-08,0.61E-12,
     &0.71E-12,0.80E-12,0.88E-12,0.98E-12,0.12E-11,0.14E-11,0.17E-11,
     &0.21E-11,0.29E-11,0.34E-11,0.41E-11,0.40E-11,0.52E-11,0.54E-11,
     &0.42E-11,0.34E-11,0.18E-11,0.46E-12,0.32E-11,0.19E-10,0.71E-10,
     &0.22E-09,0.51E-09,0.86E-09,0.14E-08,0.91E-12,0.11E-11,0.13E-11,
     &0.14E-11,0.15E-11,0.16E-11,0.20E-11,0.24E-11,0.30E-11,0.39E-11,
     &0.51E-11,0.62E-11,0.76E-11,0.82E-11,0.99E-11,0.10E-10,0.85E-11,
     &0.62E-11,0.32E-11,0.81E-12,0.93E-11,0.61E-10,0.19E-09,0.51E-09,
     &0.90E-09,0.15E-08,0.13E-11,0.15E-11,0.17E-11,0.20E-11,0.23E-11,
     &0.27E-11,0.32E-11,0.42E-11,0.55E-11,0.66E-11,0.87E-11,0.13E-10,
     &0.17E-10,0.22E-10,0.24E-10,0.28E-10,0.28E-10,0.24E-10,0.19E-10,
     &0.93E-11,0.26E-11,0.29E-10,0.16E-09,0.46E-09,0.89E-09,0.15E-08,
     &0.20E-11,0.22E-11,0.24E-11,0.25E-11,0.27E-11,0.34E-11,0.50E-11,
     &0.70E-11,0.10E-10,0.14E-10,0.16E-10,0.21E-10,0.35E-10,0.52E-10,
     &0.67E-10,0.70E-10,0.83E-10,0.84E-10,0.71E-10,0.61E-10,0.29E-10,
     &0.84E-11,0.97E-10,0.39E-09,0.82E-09,0.15E-08,0.28E-11,0.32E-11,
     &0.36E-11,0.40E-11,0.44E-11,0.49E-11,0.69E-11,0.11E-10,0.17E-10,
     &0.24E-10,0.33E-10,0.44E-10,0.63E-10,0.92E-10,0.14E-09,0.17E-09,
     &0.19E-09,0.22E-09,0.22E-09,0.19E-09,0.16E-09,0.97E-10,0.38E-10,
     &0.28E-09,0.73E-09,0.14E-08,0.40E-11,0.47E-11,0.54E-11,0.62E-11,
     &0.80E-11,0.10E-10,0.13E-10,0.19E-10,0.30E-10,0.43E-10,0.58E-10,
     &0.76E-10,0.12E-09,0.17E-09,0.22E-09,0.32E-09,0.41E-09,0.47E-09,
     &0.51E-09,0.51E-09,0.46E-09,0.39E-09,0.28E-09,0.12E-09,0.58E-09,
     &0.13E-08,0.11E-10,0.14E-10,0.18E-10,0.22E-10,0.28E-10,0.44E-10,
     &0.65E-10,0.90E-10,0.12E-09,0.17E-09,0.22E-09,0.26E-09,0.31E-09,
     &0.38E-09,0.46E-09,0.54E-09,0.67E-09,0.78E-09,0.86E-09,0.90E-09,
     &0.89E-09,0.82E-09,0.73E-09,0.58E-09,0.30E-09,0.11E-08,0.17E-10,
     &0.25E-10,0.35E-10,0.46E-10,0.58E-10,0.75E-10,0.12E-09,0.19E-09,
     &0.26E-09,0.36E-09,0.47E-09,0.60E-09,0.71E-09,0.82E-09,0.92E-09,
     &0.10E-08,0.12E-08,0.13E-08,0.14E-08,0.15E-08,0.15E-08,0.15E-08,
     &0.14E-08,0.13E-08,0.11E-08,0.66E-09
     & /

      end module TOMAS_ACTV

      
      SUBROUTINE CALCNd (TPARC,PPARC,TPi,MLi,NSECi,WPARC,NACT,SMAX
     &                   ,RHOSI,GCMLWMR,EPSILON,AUTO,DIFFLWMR,DIFFEPS,
     &     pland)
     
C
cYHL      INCLUDE 'calcNdNS.inc'
      USE TOMAS_ACTV
      USE TRACER_COM,only : NBINS ! nbins=nseci
      IMPLICIT NONE
      integer nseci
      REAL*8 TPi(nseci),MLi(nseci),TPARC,PPARC,NACT,WPARC,
     &                 SMAX,
     &                 GCMLWMR, EPSILON,AUTO(6),RHOSI,DIFFLWMR,DIFFEPS
      real*8           PLAND
      LOGICAL EX
C
      ACCOM  = 0.06d0                                  ! Accommod.coeff.

      CALL CCNSPEC (TPARC,PPARC,TPi,MLi,NSECi)        ! CCN props
c
      CALL PDFACTIV (WPARC,0d0,TPARC,PPARC,NACT,SMAX) ! Calc Nd; single W
C
!      print*,'PDFACTV',NACT,SMAX,WPARC
      NACT = MAX(NACT,40d6)    ! Restrict Minismum droplet number
      CALL GROWTH (WPARC,SMAX,NACT,RHOSI,GCMLWMR,
     &             EPSILON,AUTO,DIFFLWMR,DIFFEPS,pland)
      
!      print*,'GROWTH',NACT,SMAX,WPARC
c      print*,'NACT',i,j,l,WPARC, NACT
C
C *** END OF SUBROUTINE CALCND ****************************************
C
      END



C=======================================================================
C
C *** SUBROUTINE CCNSPEC
C *** THIS SUBROUTINE CALCULATES THE CCN SPECTRUM OF THE AEROSOL USING
C     THE APPROPRIATE FORM OF KOHLER THEORY     
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE CCNSPEC (TPARC,PPARC,TP,ML,NDATA)
C
cyhl      INCLUDE 'calcNdNS.inc'
      USE TOMAS_ACTV
      USE CONSTANT, only : PI, grav,gasc
      IMPLICIT NONE
      REAL*8 TP(NDATA), ML(NDATA), SCMED(NSMX)
      REAL SFT
      INTEGER NDATA
      REAL*8 AKH,BKH,PAR2,SCT,ANT,TPARC,PPARC
C
C *** Setup common block variables; input to parameterization
C
      TEMP = TPARC
      PRES = PPARC
      NSEC = NDATA
C
C *** Calculate thermophysical properties
C
      CALL PROPS
      SURT = SFT(SNGL(TEMP))         ! Water Surface Tension (J m-2)
C
C *** Loop over each mode; number and critical supersaturation
C
      AKH  = 4D0*AMW*SURT/GASC/TEMP/DENW   ! Kohler param
      DO I=1,NSEC
         BKH      = ML(I)/MAX(TP(I),1d-10)     ! Moles*vhf per particle
         BKH      = BKH/(PI/6.0)/(DENW/AMW)
         PAR2     = SQRT(4d0*AKH**3/27d0/BKH)
ccc         SC(I)    = EXP(PAR2) - 1D0
ccc         SCMED(I)    = PAR2
         SC(I)    = PAR2
         NPART(I) = MAX(TP(I),1d-10)      ! Min concentration: 1d-10 particles per m3
         AKOH(I)  = AKH
!      PRINT*,'CCNSPEC',I,NSEC,AKH,BKH,SC(I)
      ENDDO
C
C *** Sort with increasing Sc
C
      DO I=1,NSEC-1
         DO J=I,NSEC
            IF (SC(I).GT.SC(J)) THEN
               SCT      = SC(I)
               ANT      = NPART(I)
ccc               AKH      = AKOH(I)   ! Not needed 'cause s.t. = const
               SC(I)    = SC(J)
               NPART(I) = NPART(J)
ccc               AKOH(I)  = AKOH(J)   ! Not needed 'cause s.t. = const
               SC(J)    = SCT
               NPART(J) = ANT
ccc               AKOH(J)  = AKH      ! Not needed 'cause s.t. = const
            ENDIF
         ENDDO
      ENDDO
      NPART(NSEC+1) = 0d0
C
C *** Calculate Sc boundaries
C
      DO I=1,NSEC
         SCMED(I) = SC(I)
      ENDDO
      DO I=2,NSEC-1
         SC(I) = SQRT(SCMED(I)*SCMED(I+1))
      ENDDO
      SC(1)    = SC(2)*(SCMED(1)/SCMED(2))
      SC(NSEC) = SC(NSEC-1)*(SCMED(NSEC)/SCMED(NSEC-1))
C
C *** END OF SUBROUTINE CCNSPEC ****************************************
C
      RETURN
      END




C=======================================================================
C
C *** SUBROUTINE PDFACTIV
C *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING 
C     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
C     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
C     and Nenes (in preparation). THIS ROUTINE CALCULATES FOR A PDF OF
C     UPDRAFT VELOCITIES.
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE PDFACTIV (WPARC,SIGW,TPARC,PPARC,NACT,SMAX)
C
cyhl      INCLUDE 'calcNdNS.inc'
      USE TOMAS_ACTV
      IMPLICIT NONE
      REAL*8 TPART, NACT, NACTI, SMAX, SIGW, TPARC, PPARC,
     * WPARC,PLIMT,PROBI,WHI,WLO,SCAL,WPI,SMAXI,
     * X1,Y1,X2,Y2,X3,Y3
      REAL             PDF
C
C *** Setup common block variables; input to parameterization
C
      TEMP = TPARC
      PRES = PPARC
C
C *** Case where updraft is very small
C
      IF (WPARC.LE.1d-6) THEN
         SMAX  = 1d-5
         NACT  = 0d0
         ISEC  = 1
         DPNMX = GREAT
         RETURN
      ENDIF
C
C *** Single updraft case
C
      IF (SIGW.LT.1d-10) THEN
         CALL ACTIVATE (WPARC,NACT,SMAX)
         WPDBG(1) = WParc                      ! Save debug info
         PDDBG(1) = 1.0
         NADBG(1) = NACT
         SMDBG(1) = SMAX
C
C *** PDF of updrafts
C
      ELSE           
         NACT  = ZERO
         SMAX  = ZERO
         PLIMT = 1e-3     ! Probability of High Updraft limit
         PROBI = SQRT(-2.0*LOG(PLIMT*SIGW*SQ2PI))
         WHI   = WPARC + SIGW*PROBI             ! Upper updrft limit
         WLO   = 0.05  ! WPARC - SIGW*PROBI     ! Low updrft limit
         SCAL  = 0.5*(WHI-WLO)                  ! Scaling for updrafts
         DO I=1,Npgauss                 
            WPI  = WLO + SCAL*(1.0-XGS(i))      ! Updraft
            CALL ACTIVATE (WPI,NACTI,SMAXI)     ! # of drops
            PDF  = (1.0/SQ2PI/SIGW)*EXP(-0.5*((WPI-WPARC)/SIGW)**2) ! Prob. of updrafts        ! Probability of updraft        ! Probability of updraft
            NACT = NACT + WGSi(i)*(PDF*NACTI)    ! Integral for drops
            SMAX = SMAX + WGSi(i)*(PDF*SMAXI)    ! Integral for Smax
            WPDBG(I) = WPI                      ! Save debug info
            PDDBG(I) = PDF
            NADBG(I) = NACTI
            SMDBG(I) = SMAXI
            IF (PDF.LT.PLIMT) GOTO 100
         ENDDO
 100     NACT = NACT*SCAL                       ! Scale Integrals
         SMAX = SMAX*SCAL
      ENDIF
C

      DPNMX = DPMX                              ! Debug info (cyhl - not in wei-chun's code)
      IF (DPMXINF) DPNMX = GREAT ! (cyhl - not in wei-chun's code)
C
      RETURN
C
C *** END OF SUBROUTINE PDFACTIV ****************************************
C
      END



C=======================================================================
C
C *** SUBROUTINE WACTIV
C *** THIS SUBROUTINE CALCULATES THE UPDRAFT NECESSARY TO ACHIEVE A DROP
C     CONCENTRATION. 
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE WACTIV (NACT, WPARC, SMAX)
C
cyhl      INCLUDE 'calcNdNS.inc'
      USE TOMAS_ACTV
      IMPLICIT NONE
      REAL*8 TPART, NACT, NACT1, NACT2, NACT3, WPARC, SMAX
     &     ,X1,Y1,X2,Y2,X3,Y3
C
C *** Initialization
C
      CALL PROPS           !  Thermophysical properties
C
C *** INITIAL VALUES FOR BISECTION **************************************
C
      X1   = 1e-3          ! Low value of updraft
      CALL ACTIVATE (X1, NACT1, SMAX)
      Y1   = NACT1/NACT - 1d0
C
      X2   = 20           ! High value of updraft
      CALL ACTIVATE (X2, NACT2, SMAX)
      Y2   = NACT2/NACT - 1d0
C
C *** PERFORM BISECTION *************************************************
C
20    DO 30 I=1,MAXIT
         X3   = 0.5*(X1+X2)
         CALL ACTIVATE (X3, NACT3, SMAX)
         Y3   = NACT3/NACT - 1d0
C
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
C
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
         NITER = I
C
30    CONTINUE
C
C *** CONVERGED ; RETURN ************************************************
C
40    X3    = 0.5*(X1+X2)
      CALL ACTIVATE (X3, NACT3, SMAX)
      Y3    = NACT3/NACT - 1d0
      WPARC = X3
C
      RETURN
C
C *** END OF SUBROUTINE ACTIVATE ****************************************
C
      END


C=======================================================================
C
C *** SUBROUTINE ACTIVATE
C *** THIS SUBROUTINE CALCULATES THE CCN ACTIVATION FRACTION ACCORDING 
C     TO THE Nenes and Seinfeld (2003) PARAMETERIZATION, WITH
C     MODIFICATION FOR NON-CONTUNUUM EFFECTS AS PROPOSED BY Fountoukis
C     and Nenes (in preparation).
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE ACTIVATE (WPARC, NACT, SMAX)
C
cyhl      INCLUDE 'calcNdNS.inc'
      USE CONSTANT, only : grav,gasc,pi
      USE TOMAS_ACTV
      IMPLICIT NONE
      REAL*8 TPART, NACT, WPARC, SMAX,BETA,CF1,CF2,
     *     X1,Y1,X2,Y2,X3,Y3,SINTEG1,SINTEG2
C
C *** Initialization
C
      CALL PROPS           !  Thermophysical properties
      WPARCEL = WPARC      !  Save input to common block
C
C *** Parameterization "coefficients"
C
      ALFA = GRAV*AMW*DHV/CPAIR/GASC/TEMP/TEMP - GRAV*AMA/GASC/TEMP
C
      BET1 = PRES*AMA/PSAT/AMW + AMW*DHV*DHV/CPAIR/GASC/TEMP/TEMP
      BET2 = GASC*TEMP*DENW/PSAT/DV/AMW/4d0 + 
     &       DHV*DENW/4d0/AKA/TEMP*(DHV*AMW/GASC/TEMP - 1d0)
      BETA = 0.5d0*PI*BET1*DENW/BET2/ALFA/WPARC/DAIR
C
      CF1  = 1.d0/SQRT(BET2*ALFA*WPARC)   ! For objective function
      CF2  = 2d0/3d0

Cyhl
      IF(NSEC.EQ.15) NSEC=NSEC-3 !for FAST-TOMAS 12bin+3nm(3bins) 
cyhl      IF(NSEC.EQ.36) NSEC=NSEC-6

C *** INITIAL VALUES FOR BISECTION **************************************
C
      X1   = SC(1)     ! MIN cloud supersaturation = first section bound
      CALL SINTEGRAL(X1,NACT,SINTEG1,SINTEG2) !NACT ==> NDRPL for wei-chun's code
      Y1   = (SINTEG1*CF1+SINTEG2*CF2)*BETA*X1 - 1d0
C
      X2   = SC(NSEC)   ! MAX cloud supersaturation = last section bound
      CALL SINTEGRAL(X2,NACT,SINTEG1,SINTEG2)
      Y2   = (SINTEG1*CF1+SINTEG2*CF2)*BETA*X2 - 1d0
C
      IF (SIGN(1.d0,Y1).LT.ZERO .AND. SIGN(1.d0,Y2).LT.ZERO) THEN  ! Double-check bounds
         NACT = NPART(NSEC)                    ! Condensation always .lt. supply --> all activate    
         SMAX = SC(NSEC)
         RETURN 
      ELSEIF (SIGN(1.d0,Y1).GT.ZERO .AND. SIGN(1.d0,Y2).GT.ZERO) THEN
         NACT = NPART(2)                    ! Condensation always .gt. supply --> almost none activate    
         SMAX = SC(2)
         RETURN 
      ENDIF
C
C *** PERFORM BISECTION *************************************************
C
20    DO 30 I=1,MAXIT
         X3   = 0.5*(X1+X2)
         CALL SINTEGRAL(X3,NACT,SINTEG1,SINTEG2)
         Y3 = (SINTEG1*CF1+SINTEG2*CF2)*BETA*X3 - 1d0
C
         IF (SIGN(1.d0,Y1)*SIGN(1.d0,Y3) .LE. ZERO) THEN  ! (Y1*Y3 .LE. ZERO)
            Y2    = Y3
            X2    = X3
         ELSE
            Y1    = Y3
            X1    = X3
         ENDIF
C
         IF (ABS(X2-X1) .LE. EPS*X1) GOTO 40
         NITER = I
C
30    CONTINUE
C
C *** CONVERGED ; RETURN ************************************************
C
40    X3   = 0.5*(X1+X2)
      CALL SINTEGRAL(X3,NACT,SINTEG1,SINTEG2)
      Y3   = (SINTEG1*CF1+SINTEG2*CF2)*BETA*X3 - 1d0
      SMAX = X3
C
      RETURN
C
C *** END OF SUBROUTINE ACTIVATE ****************************************
C
      END


C======================================================================C
C *** SUBROUTINE GROWTH
C *** THIS SUBROUTINE CALCULATES THE CONDENSATION GROWTH OF DROPLETS
C     AFTER THE SMAX
C
C=====================================================================      
      SUBROUTINE GROWTH (WPARC,SMAX,NACT,
     &                   RHOSI,GCMLWMR,EPSILON,AUTO,DIFFLWMR,DIFFEPS,
     &                   pland)
      USE TOMAS_ACTV
      IMPLICIT NONE
      INTEGER NMDM,ISW
      PARAMETER (NMDM=3)     ! MAX # OF LOGNORMAL MODES.
c      PARAMETER        (NSEC=50,NSEC3=150)

      REAL*8 WPARC, NACT,SMAX,EPSFIT,
     &                 TPI(NMDM), DPGI(NMDM),  SIGI(NMDM),
c$$$     &                 NPART1(NSEC),DPART1(NSEC),
c$$$     &                 NPART2(NSEC),DPART2(NSEC),
c$$$     &                 NPART3(NSEC),DPART3(NSEC),
c$$$     &                 NPART(NSEC3),DPARTS(NSEC3),
c$$$     &                 SC(NSEC3),DPSMAX(NSEC3),
     &                 DPMEAN,SIGMA,NISUM,EPSILON,LWC,
     &                 DELTAT,
     &                 DENA,SOLD,SNEW,LWMR,DWLDT,DSDT,
     &                 P6,P4,KKAUTO,MC,BH94,PGISS,AUTOK,
     &                 DPGROW,NGROW,EPSILONSMAX,RHOSI,
     &                 DPFIT(20),NFIT(20),GCMLWMR,DIFFLWMR,DIFFEPS,
     &     AUTO(6),DIFFND              ! all autoconversion I need -Yunha Lee
      real*8           PLAND
C
C *** SETUP CONSTANTS
C
      EPSILON = 0.44

      ISW=4
      IF (GCMLWMR .EQ. 0 ) THEN
      AUTO(1) = 0
      AUTO(2) = 0
      AUTO(3) = 0
      AUTO(4) = 0
      AUTO(5) = 0
      AUTO(6) = 0
C      AUTO = 0
      ELSE
      GCMLWMR = MAX(GCMLWMR,1D-7)

cyhl this is for SCE          
      CALL FIT (EPSILON,GCMLWMR,NACT,RHOSI,DPFIT,NFIT)
      CALL CLDPROP (DPFIT,NFIT,20,RHOSI,NISUM,DPMEAN,SIGMA,EPSFIT,LWMR)
      DIFFLWMR = (GCMLWMR-LWMR)/GCMLWMR
      DIFFND = NACT - NISUM
      DIFFEPS = (EPSILON-EPSFIT)/EPSILON

      CALL SCE (DPFIT,NFIT,20,ISW,AUTOK,pland)
cyhl up here. 

CYHL this is to autoconversion parameterization. 
      DENA  = 1.2D0
      EPSILON = 0.44
      LWC=GCMLWMR*DENA
      CALL AUTOCON (LWC,NACT,EPSILON,P6,P4,KKAUTO,MC,BH94,PGISS)
CYHL      AUTO = KKAUTO

      AUTO(1)=P6
      AUTO(2)=KKAUTO
      AUTO(3)=MC
      AUTO(4)=BH94
      AUTO(5)=PGISS
      AUTO(6)=AUTOK

cyhl up here. 

c       PRINT*,'AUTOCON',AUTO(1),AUTO(4),AUTO(6)
       if(ISNAN(AUTO(6)))THEN
          print*,'wrong NS Dpfit',Dpfit,'Nfit',Nfit,'Nsum',Nisum,
     &         'dpmean',dpmean, 'sigma',sigma,'epsfit',epsfit,
     &         'lwmr',lwmr,'gcmlwmr',GCMLWMR,'nact',nact 
          call stop_model('Wrong AUTO in TOMAS_act',255)
       ENDIF
      ENDIF


      
c       PRINT*,AUTO,NACT,555
C      PRINT*,AUTO,133
C      DO 70 I =1, M
C      DPOLD(I) = DPGROW(I)
C70    CONTINUE

C
C***  TIME ITERATION FOR SUPERSATURATION AND DROPLET GROWTH
C

C      DELTAT = 1
C      J = 1
C      LWC = 0D0

C      DO WHILE ((GCMLWC-LWC) .GT. 1D-5)
      
C      DWLDT = 0
      
C      DO 80 I =1, M
      
C      DWLDT = DWLDT + NGROW(I)*DPOLD(I)

C80    CONTINUE

C      DWLDT = PI/2D0*DENW/DENA/BET2*SOLD*DWLDT
C      DSDT = ALFA*WPARCEL -BET1*DWLDT
C      SNEW = SOLD + DSDT*DELTAT
      
C      SNEW = MAX( 0D0, SNEW)

C      DO 90 I =1, M

C      DPNEW(I) = SQRT(DPOLD(I)**2 + 2/BET2*(SOLD*DELTAT + 0.5*
C     &DSDT*DELTAT**2))
C90    CONTINUE

C      CALL CLDPROP (DPNEW,NGROW,M,NISUM,DPMEAN,SIGMA,EPSILON,LWMR)
C      LWC=LWMR*DENA

C      CALL AUTOCON (LWC,NISUM,EPSILON,P6,P4,KKAUTO,MC,BH94,PGISS)

C      SOLD = SNEW

C     DO 110 I = 1, M
C      DPOLD(I) = DPNEW(I)
C110   CONTINUE

C      J = J + 1
C      ENDDO !TIME LOOP
C      ENDDO !WHILE
      RETURN
C
C *** END OF SUBROUTINE GROWTH ****************************************
C
      END

C======================================================================C
C *** SUBROUTINE CLDPROP
C *** THIS SUBROUTINE CALCULATES THE  DPMEAN, SIGMA AT EACH CLOUD HEIGHT
C
C
C
C *** WRITTEN BY WEI-CHUN HSIEH
C
C *** Small change by YUNHA LEE, 07.18.2011 
C     (DENA TO RHOSI) Because  DENA is actually RHOSI. 
C======================================================================C
!      SUBROUTINE CLDPROP (DPNEW,NPART,NSEC3,DENA,NISUM,DPMEAN,SIGMA,
!     &           	EPSILON,LWMR)
 
      SUBROUTINE CLDPROP (DPNEW,NPART,NSEC3,RHOSI,NISUM,DPMEAN,SIGMA,
     &           	EPSILON,LWMR)
     
      USE CONSTANT, only: gasc,pi,grav
      REAL*8 NIDP3,DPNEW(NSEC3),NPART(NSEC3),LWMR,
     &                 NITDPSUM,NISUM,DPMEAN,SIGSUM,SIGMA,EPSILON
     *     ,DENW,RHOSI !RHOSI is added by YUNHA LEE 
      INTEGER I,J

      DENW  = 1D3
C     DENA  = 1.2D0 

      NITDPSUM = 0D0
      NISUM = 0D0
      NIDP3 = 0D0

      DO I = 1, NSEC3
      NITDPSUM = NITDPSUM + DPNEW(I)*NPART(I)
      NISUM = NISUM + NPART(I)
      NIDP3 = NIDP3 + NPART(I)*DPNEW(I)**3D0
      ENDDO

      IF(NISUM.EQ.0.) PRINT*,'ZERO NISUM',NPART

      DPMEAN = NITDPSUM/NISUM
CYHL      LWMR = DENW/DENA*PI/6D0*NIDP3
      LWMR = DENW/RHOSI*PI/6D0*NIDP3
      
      SIGSUM = 0D0
      DO I = 1, NSEC3
      IF ( NPART(I) .GT. 0D0 ) THEN
      SIGSUM = SIGSUM + (DPNEW(I)-DPMEAN)**2*NPART(I)
      ENDIF
      ENDDO

      SIGMA = SQRT(SIGSUM/NISUM)

      EPSILON = SIGMA/DPMEAN

      RETURN
      END


C======================================================================C
C *** SUBROUTINE AUTOCON
C *** THIS SUBROUTINE CALCULATES THE AUTOCONVERSION RATES
C
C
C
C *** WRITTEN BY WEI-CHUN HSIEH
C
C======================================================================C
      SUBROUTINE AUTOCON (LWC,ND,EPSILON,P6,P4,KKAUTO,MC,BH94,PGISS)
      USE CONSTANT, only : pi
      REAL*8 LWC,LWMR,ND,EPSILON,P6,P4,KKAUTO,MC,BH94,PGISS,
     &                 K1,K2,BETA6,ALPHA6,BETA4,ALPHA4,ALPHAMC,K,DENW
     
      DENW  = 1D3
      DENA  = 1.2D0
      
      DENW=1000D0
      K1=1.19D8
      K2=1.9D17
      BETA6=((1D0+3D0*EPSILON**2D0)*(1D0+4D0*EPSILON**2D0)*
     & (1D0+5D0*EPSILON**2D0)/
     & (1D0+EPSILON**2D0)/(1D0+2D0*EPSILON**2D0))**(1D0/6D0)
      ALPHA6=(3D0/4D0/PI/DENW)**2*K2*BETA6**6D0*(LWC/ND)**(2D0/3D0)
      P6 = ALPHA6*ND**(-1D0/3D0)*LWC**(7D0/3D0)

      BETA4=(1D0+3D0*EPSILON**2D0)**(1D0/4D0)/((1D0+2D0*EPSILON**2D0)*
     & (1D0+EPSILON**2D0))**(1D0/12D0)
      ALPHA4=PI*K1*(3D0/4D0/PI/DENW)**(4D0/3D0)*0.55D0*BETA4**4D0
      P4=ALPHA4*ND**(-1D0/3D0)*LWC**(7D0/3D0)

      LWMR=LWC/1.2D0
      KKAUTO=1350D0*LWMR**(2.47D0)*(ND*1D-6)**(-1.79D0)

      ALPHAMC = PI*K1*(3D0/4D0/PI/DENW)**(4D0/3D0)*0.55D0 !0.55 EMC
      MC = ALPHAMC*ND**(-1D0/3D0)*(LWC)**(7D0/3D0)

      K=EPSILON**(-2D0)
      BH94=6D25*(K-1D0)**(-1.7D0)*(ND*1D-6)**(-3.3D0)*(LWC*1D-3)**
     & (4.7D0)
      BH94=BH94*1D3 !%CONVERT CGS TO MKS UNIT
      
C GISS OCEAN
      PGISS = 1D-4*LWMR*(1D0-EXP(-((LWC/0.0005D0)**4D0)))

C GISS LAND
C      PGISS = 1D-4*LWMR*(1D0-EXP(-((LWC/0.001D0)**4D0)))

      RETURN
      END

C======================================================================C
C *** SUBROUTINE FIT
C *** THIS SUBROUTINE GIVE FITTED GAMMA DISTRIBUTION
C
C
C
C *** WRITTEN BY WEI-CHUN HSIEH
C
C======================================================================C
C      SUBROUTINE FIT(EPSILON,LWMR,ND,DPART,NPART)
      SUBROUTINE FIT(EPSILON,LWMR,ND,RHOSI,DBIN,NBIN)
      USE CONSTANT, only: gasc,pi,grav
      PARAMETER (NSEC=40)
      PARAMETER (N=40)
      INTEGER I,J
      REAL*8 DENW,DENA,EMIN

      DATA DENW    /1000D0/
      DATA DENA    /1.2D0/
      DATA EMIN /1.E-9/
      REAL*8 K,THETA,N0,LWC,LWMR,ND,EPSILON,DMX,DMN,RHOSI,
     &                 DLND,DAMH(NSEC+1),DPART(NSEC),NPART(NSEC),NTOT,
     &                 GA,GK,GKP3,GIN,GIM,GIP,E(N),R(N),D(N),NBOTT(N),
     &                 EB(N+1),RB(N+1),DB(N+1),NBIN(20),DBIN(20),NTOT2,
     &                 GINSMALL,GINLARGE
      DMX = 50D-6
      DMN = 0.4D-6

      LWC = LWMR*RHOSI
      K = EPSILON**(-2D0)

      CALL GAMMA(K,GA)
      GK=GA
      CALL GAMMA(K+3,GA)
      GKP3=GA

      THETA = (6D0/PI*LWC/DENW/ND*GK/GKP3)**(1D0/3D0)
      N0 = ND/GK/THETA**K
      
      SCAL=2
      AX=2.D0**(1.0/SCAL)
C MASS AND RADIUS GRID
      E(1)=EMIN*0.5*(AX+1.)
      R(1)=1000.*DEXP(DLOG(3.*E(1)/(4.*PI))/3.)
      D(1)=R(1)*2.E-6

      EB(1)=EMIN
      RB(1)=1000.*DEXP(DLOG(3.*EB(1)/(4.*PI))/3.)
      DO I=2,N
         E(I)=AX*E(I-1)
         R(I)=1000.*DEXP(DLOG(3.*E(I)/(4.*PI))/3.)
         D(I)=R(I)*2.D-6
         EB(I)=AX*EB(I-1)
         RB(I)=1000.*DEXP(DLOG(3.*EB(I)/(4.*PI))/3.)
         DB(I)=RB(I)*2.D-6
         IF (R(I) .GE. 2. .AND. R(I) .LE. 40.) THEN
         ENDIF
      ENDDO
         EMAX=AX*EB(N)
         RB(N+1)=1000.*DEXP(DLOG(3.*EMAX/(4.*PI))/3.)
         DB(N+1)=RB(N+1)*2.D-6
      
      
      DO I=1,N
      NBOTT(I)=0.
      ENDDO
      
      NTOT = 0
      DO I=11,30
      CALL INCOG(K,DB(I)/THETA,GIN,GIM,GIP)
      GINSMALL = GIN
      CALL INCOG(K,DB(I+1)/THETA,GIN,GIM,GIP)
      GINLARGE = GIN
      NBOTT(I) = N0*THETA**K*(GINLARGE-GINSMALL)
      NTOT = NTOT+NBOTT(I)
cyhl Added by YUNHA Lee 
c  DEC 2011 - NTOT is zero because GINLARGE and GINSMALL are same.  
c  However, DB(I) and DB(I+1) are different. 
c  When does gamma function give the same values for different DB(I)/THETA?
      if(I.EQ.11.AND.GINLARGE.EQ.GINSMALL)THEN
         PRINT*,'NTOT=0 BUT SET TO 100.'
         NTOT=100.
      ENDIF
cyhl
      ENDDO
      
      NTOT2 = 0 
      M = 1
      DO I =11,30
      DBIN(M)=D(I)
      NBIN(M)=NBOTT(I)/NTOT*ND
      NTOT2 = NTOT2+NBIN(M)
      M = M + 1
      ENDDO

CYHL      if(NTOT2.EQ.0.)PRINT*,'NTOT2',ND,NTOT,THETA,N0,
CYHL     & GINLARGE,GINSMALL,EPSILON
C      CALL INCOG(A,X,GIN,GIM,GIP)

      RETURN
      END


C======================================================================C
C *** SUBROUTINE SCE
C *** THIS SUBROUTINE CALCULATES SCE AUTOCONVERSION RATES
C
C
C
C *** WRITTEN BY WEI-CHUN HSIEH
C
C======================================================================C
      SUBROUTINE SCE (DPNEW,NPAR,M,ISW,AUTOK,pland)
      USE TOMAS_ACTV,ONLY : GRV,TUROCEAN,TUREARTH
      USE CONSTANT, only : PI,gasc,grav

      IMPLICIT NONE
      INTEGER NH,M
      PARAMETER (NH=400)
      INTEGER ISW,I,J
      CHARACTER*80 ISWSTR
      REAL*8 DPNEW(M),NPAR(M),M0,AUTOK,KERNEL,
     &                 MASS,R,CK,CKT,CK0,D0,pland,DENW
C      DIMENSION MASS(NH),R(NH),CK(NH,NH),CKT(NH,NH),CK0(NH,NH)
      DIMENSION MASS(M),R(M),CK(M,M),CKT(M,M),CK0(M,M)
C     &          CKISW0(N,N),CKISW4(N,N)
C      COMMON /BOTT/ CKISW0(N,N),CKISW4(N,N)
C      DIMENSION MASS(100),R(100),CK(100,100),CKT(100,100),CK0(100,100)

      DENW  = 1D3
      
      D0=40D-6
C%D0=50E-6;
      M0=PI/6D0*D0**3D0*DENW
      
      DO I = 1,M
      MASS(I) = PI/6D0*DENW*DPNEW(I)**3D0
      R(I) =DPNEW(I)/2D0*1D6 !% R IN UNIT MICRO
      ENDDO


      AUTOK=0

CYHL      IF (ITYP .EQ. 1) THEN
      DO I=1,M
      DO J=1,M
      CK(I,J)=(1.d0-pland)*TUROCEAN(I,J) + pland*TUREARTH(I,J)   ! Assemble diagnostics
      ENDDO
      ENDDO
CYHL      ELSEIF (ITYP .EQ. 2) THEN
CYHL      DO I=1,M
CYHL      DO J=1,M
CYHL      CK(I,J)=TUREARTH(I,J)
CYHL      ENDDO
CYHL      ENDDO
CYHL      ENDIF 

      DO I = 1,M
      DO J = I,M
        IF ((MASS(I)+MASS(J)) .GT. M0) THEN

         AUTOK = AUTOK +  CK(I,J)*MASS(I)*NPAR(J)*NPAR(I)

C         ICOUNT = ICOUNT + 1

        ENDIF

C        JCOUNT = JCOUNT + 1
      ENDDO
      ENDDO
      
      RETURN
      END


C======================================================================C
C *** SUBROUTINE KERNEL
C *** THIS SUBROUTINE CALCULATES LONG'S KERNEL (POLYNOMIAL FORM)
C
C
C
C *** WRITTEN BY WEI-CHUN HSIEH
C
C======================================================================C
      SUBROUTINE LONG(D1,D2,KERNEL)
C      FUNCTION RATE = KERNAL(D1,D2)
      REAL*8 D1,D2,K1,K2,KERNEL
      K2 = 2.59D15
      K1 = 3.03D3

C      D1,D2 IN METERS

C%D1 = MAX(D1,D2);
C%D2 = MIN(D1,D2);

      IF (D1 .LT. 100D-6) THEN
          KERNEL = K2*(D1**6D0+D2**6D0)
      ELSE
          KERNEL = K1*(D1**3D0+D2**3D0)
      ENDIF
      
      RETURN
      END

C=======================================================================
C
C *** SUBROUTINE SINTEGRAL
C *** THIS SUBROUTINE CALCULATES THE CONDENSATION INTEGRALS, ACCORDING
C     TO THE POPULATION SPLITTING ALGORITHM OF Nenes and Seinfeld (2003)
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE SINTEGRAL (SPAR,NACT,INTEGMAT,INTEGDC)
C
      USE TOMAS_ACTV
      IMPLICIT NONE
      REAL*8 NACT,NACT1,NACT2,NACT3,NACT4,SPART, SPAR,
     &                 INTEG1,INTEG2,INTEG3,INTEG4, INTEGMAT,INTEGDC
     * ,AKOHI,SCMIN,DESCR,SSPLT1,SSPLT2,DIF,CURDIF
      REAL             SFT
C
C *** Calculate ISEC that corresponds to SPAR (i.e. Smax)
C
cyhl this should be changed for different lowest size boundary. 

      IF(NSEC.EQ.15) NSEC=NSEC-3 
cyhl 

      IF (SC(1).GT.SPAR) THEN
         ISEC = 1
      ELSE IF (SC(NSEC).LT. SPAR) THEN
         ISEC = NSEC + 1
      ELSE
         DO I=2,NSEC
            ISEC = I
            IF (SC(I-1).LE.SPAR .AND. SPAR.LT.SC(I)) GOTO 10
         ENDDO
      ENDIF
C
C *** Calculate population boundaries; assign default values
C
  10  ISLW1 = ISEC        ! Initial values
      ISLW2 = ISEC
      ISMIN = ISEC
C
C *** DPMX is size below which CCN are kinetically limited.
C
C * To calculate Spart, assume that surface tension = water
C
      SURT   = SFT(SNGL(TEMP))       ! Surface Tension for water (J m-2)
      AKOHI  = 4D0*18d-3*SURT/8.314d0/TEMP/1d3            ! Kohler param
C
      IF (IDPMF.EQ.0) THEN            ! Regression fit
         DPMX   = 1d-7*SPAR**(-0.6176)   ! Average fit
      ELSEIF (IDPMF.EQ.1) THEN
         DPMX   = 3d-7*SPAR**(-0.5056)   ! Single mode
      ELSE
         DPMX   = 7d-8*SPAR**(-0.6415)   ! Tri-modal fit
      ENDIF
      SCMIN  = (2d0/3d0)*AKOHI/DPMX   ! Minimum Scritical
C
C *** Here is where the criterion with the descriminant is put. When it
C     is < 0, then set CRIT2 = .TRUE. Otherwise, set the two values of
C     SSPLT and continue.
   
      DESCR  = 1d0 - (16d0/9d0)*ALFA*WPARCEL*BET2*(AKOHI/SPAR**2)**2
      IF (DESCR.LE.0d0) THEN
         CRIT2  = .TRUE.             ! SSPLT1,SSPLT2 do not exist
      ELSE
         CRIT2  = .FALSE.
         SSPLT1 = 0.5d0*(1d0-SQRT(DESCR)) ! min root of both
         SSPLT2 = 0.5d0*(1d0+SQRT(DESCR)) ! max root of both
         SSPLT1 = SQRT(SSPLT1)*SPAR       ! Multiply ratios with Smax
         SSPLT2 = SQRT(SSPLT2)*SPAR
      ENDIF
C
C *** Now scan for the appropriate indices in the CCN spectrum.
C
      IF (.NOT.CRIT2) THEN           ! ISLW1,ISLW2 exist: find them.
         DIF = 1d30                                     ! Find ISLW1
         IF (ISEC.GT.1) THEN
            DO I=MIN(ISEC,NSEC),2,-1
               CURDIF = SC(I)-SSPLT1
               IF(CURDIF.LT.0d0) GOTO 20
               IF(CURDIF .LT. DIF) THEN
                 ISLW1 = I
                 DIF  = CURDIF
               ENDIF
            ENDDO
         ENDIF
C
  20     DIF = 1d30                                     ! Find ISLW2
         IF (ISEC.GT.1) THEN
            DO I=MIN(ISEC,NSEC),2,-1
               CURDIF = SC(I)-SSPLT2
               IF(CURDIF.LT.0d0) GOTO 25
               IF(CURDIF .LT. DIF) THEN
                 ISLW2 = I
                 DIF  = CURDIF
               ENDIF
            ENDDO
         ENDIF
C
      ENDIF
C
  25  DIF = 1d30                       ! Determine ISMIN
      IF (ISEC.GT.1) THEN
         DO I=MIN(ISEC,NSEC),1,-1
            CURDIF = SC(I)-SCMIN
            IF(CURDIF.LT.0d0) GOTO 30
            IF(CURDIF .LT. DIF) THEN
              ISMIN = I
              DIF   = CURDIF
            ENDIF
         ENDDO
      ENDIF
C
C *** Now, assign the population flags their appropriate values
C
C      CRIT1  = .TRUE. when ISMIN >= ISLW1
C      CRIT2  = .TRUE. when ISLW1  = ISLW2, or ISLWs nonexistent
C      CRIT3  = .TRUE. when ISLW2  = ISEC
C
  30  IF (INTEGMOD.EQ.0) THEN       ! Normal integration mode
         IF (ISEC.NE.1) THEN        ! IF ISEC is more than 1
            CRIT2 =  CRIT2 .OR. ISLW1.EQ.ISLW2  ! Second criterion
            IF (.NOT.CRIT2) THEN
               ISMIN = ISLW1-1      !Point of difference from CASE09
               CRIT3 = ISLW2.EQ.ISEC            ! Third criterion
            ELSE
               CRIT3 =.FALSE.
               ISLW1 = ISMIN+1                  ! To make integration right
               ISLW2 = ISMIN+1       
            ENDIF
            CRIT1 = ISLW1.EQ.ISMIN+1            ! First criterion
         ELSE                       ! ISEC = 1, integrate with Dp >> Dc
            CRIT1 = .TRUE.
            CRIT2 = .TRUE.
            CRIT3 = .TRUE.
            ISMIN = 1
         ENDIF
C
      ELSEIF (INTEGMOD.EQ.1) THEN   ! Dp >> Dc integration mode for all CCN
         CRIT1 = .TRUE.
         CRIT2 = .FALSE.
         CRIT3 = .TRUE.
         ISMIN = 1
         ISLW1 = 2
         ISLW2 = ISEC
      ELSEIF (INTEGMOD.EQ.2) THEN   ! Dp ~ Dc integration mode for all CCN
         CRIT1 = .TRUE.
         CRIT2 = .TRUE.
         CRIT3 = .FALSE.
         ISMIN = 1
         ISLW2 = 2
         ISLW1 = 2
      ELSE
         STOP 'Bad Integral mode specified. Abort'
      ENDIF   
C
C *** Calculate integrals depending on the populations existing.
C
      CALL SINTEG1 (SPAR,1,ISMIN,NACT1,INTEG1) ! Inertial CCN
C
      IF (CRIT1) THEN              ! Large CCN 
         NACT2  = 0d0
         INTEG2 = 0d0
      ELSE
         CALL SINTEG2 (SPAR,ISMIN+1,ISLW1,NACT2,INTEG2)
      ENDIF
C
      IF (CRIT2) THEN              ! Mature drops
         NACT3  = 0d0
         INTEG3 = 0d0
      ELSE
         CALL SINTEG1 (SPAR,ISLW1+1,ISLW2,NACT3,INTEG3)
      ENDIF
C
      IF (CRIT3) THEN              ! High Sc (recently activated) CCN
         NACT4  = 0d0
         INTEG4 = 0d0
      ELSE
         CALL SINTEG2 (SPAR,ISLW2+1,ISEC,NACT4,INTEG4)
      ENDIF
C
      NACT     = NACT1 + NACT2 + NACT3 + NACT4
      INTEGMAT = INTEG1 + INTEG3   ! Integral from "mature drops"
      INTEGDC  = INTEG2 + INTEG4   ! Integral from "recently activated"
C
      RETURN
C
C *** END OF SUBROUTINE SINTEGRAL ***************************************
C
      END



C=======================================================================
C
C *** SUBROUTINE SINTEG1
C *** This routine computes the condensation rate with the assumption
C     that the growth term is very large compared to Dc, or Dc is too
C     large to be attained during the droplet growth.
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE SINTEG1 (SPAR,IST,IEN,NACT,INTEG)
C
      USE TOMAS_ACTV
      IMPLICIT NONE
      INTEGER IEN,IST
      REAL*8 NACT, SPART, INTEG, SPAR,SJM1,SJ,DELS,BJM1,
     & BJ
C
C *** Calculate integral
C
cyhl this should be changed for different lowest size boundary. 

      IF(NSEC.EQ.15) NSEC=NSEC-3 
cyhl 
      INTEG = 0d0
      DO I=IST,IEN
         IF (I.EQ.1) THEN         ! S(j-1)
            SJM1 = 0d0
         ELSE
            SJM1 = SC(I-1)
         ENDIF
         SJ = SC(I)
C
         IF (I.EQ.1) THEN         ! S(j)-S(j-1)
            DELS = SC(1)
         ELSEIF (I.EQ.NSEC+1) THEN
            DELS = (SPAR-SC(NSEC))
         ELSE
            DELS = (SC(I)-SC(I-1))
         ENDIF
C
         BJM1=SJM1*SQRT(SPAR*SPAR-SJM1*SJM1) + 
     &             SPAR*SPAR*DASIN(MIN(SJM1,SPAR)/SPAR)
         BJ  =SJ  *SQRT(MAX(SPAR*SPAR-SJ*SJ,+0d0)) + 
     &             SPAR*SPAR*DASIN(MIN(SPAR,SJ)/SPAR)
C
Cyhl  Added by Yunha Lee - TO AVOID INTEG IS TO "NaN". 
         IF(I.GT.1)THEN
           if(SC(I).EQ.SC(I-1))THEN
             BJ=BJM1
             DELS=1
           ENDIF                ! SET BY YHL 05/20/2009
         ENDIF

         INTEG = INTEG + 0.5d0*(BJ-BJM1)*NPART(I)/DELS
      ENDDO
C
C *** Calculate # of drops
C
      IF (SC(IEN).GT.SPAR) THEN   ! # of drops, last section, if Sparc < SJ
         NACT = NPART(IEN)*(SPAR-SJM1)/DELS
      ELSE
         NACT = NPART(IEN)
      ENDIF
C
      DO I=IEN-1,IST,-1
         NACT = NACT + NPART(I)
      ENDDO
C
      RETURN
C
C *** END OF SUBROUTINE SINTEG1 *****************************************
C
      END



C=======================================================================
C
C *** SUBROUTINE SINTEG2
C *** This routine computes the condensation rate with the assumption 
C     the droplet size is of order Dc.
C
C *** WRITTEN BY ATHANASIOS NENES
C
C=======================================================================
C
      SUBROUTINE SINTEG2 (SPAR,IST,IEN,NACT,INTEG)
C
      USE TOMAS_ACTV
      IMPLICIT NONE
      INTEGER IEN,IST
      REAL*8 NACT, SPART, INTEG, SPAR,SJM1,SJ,DELS,BJM1,
     & BJ
C
      IF (IST.EQ.1) STOP 'Bad IST in SINTEG2'
C
C *** Calculate integral
C
cyhl this should be changed for different lowest size boundary. 

      IF(NSEC.EQ.15) NSEC=NSEC-3 
cyhl 

      INTEG = 0d0
      SJM1  = SPAR
      DELS  = 1d0
C
      DO I=IST,IEN
         SJM1 = SC(I-1)             ! S(j-1)
         IF (I.EQ.IEN) THEN         ! S(j)
            SJ = SPAR
         ELSE
            SJ = SC(I)
         ENDIF
         IF (I.EQ.NSEC+1) THEN      ! S(j) - S(j-1)
            DELS = (SPAR-SC(NSEC))
         ELSE
            DELS = (SC(I)-SC(I-1))
         ENDIF

Cyhl     TO AVOID INTEG IS TO "NaN". 
         if(I.GT.1)THEN
           if(SC(I).EQ.SC(I-1))THEN
             BJ=BJM1
             DELS=1
           ENDIF                ! SET BY YHL 05/20/2009
         ENDIF

         INTEG = INTEG + AKOH(I)*LOG(SJ/SJM1)*NPART(I)/DELS
      ENDDO
C
C *** Calculate # of drops
C                                                      
      IF (SC(IEN).GT.SPAR) THEN     ! # of drops, last section, if Sparc < SJ
         NACT = NPART(IEN)*(SPAR-SJM1)/DELS
      ELSE
         NACT = NPART(IEN)
      ENDIF
C
      DO I=IEN-1,IST,-1
         NACT = NACT + NPART(I)
      ENDDO
C
      RETURN
C
C *** END OF SUBROUTINE SINTEG2 *****************************************
C
      END


C======================================================================C
C *** SUBROUTINE PROPS
C *** THIS SUBROUTINE CALCULATES THE THERMOPHYSICAL PROPERTIES
C
C *** WRITTEN BY ATHANASIOS NENES
C
C======================================================================C
      SUBROUTINE PROPS
      USE TOMAS_ACTV
      USE CONSTANT, only: gasc,pi,grav,lhe
      IMPLICIT NONE
      REAL*8 DBIG,DLOW,COEF,PRESA

      REAL  VPRES, SFT
      real( kind=8 ) :: wv_psat ! external function, jan perlwitz, Nov 2017
C
      DENW  = 1D3                           ! WATER DENSITY
      DHV   = 2.25D6                        ! WATER ENTHALPY OF VAPORIZATION
      CPAIR = 1.0061D3                      ! AIR CP
      PRESA = PRES/1.013D5                  ! PRESSURE (PA)
      DAIR  = PRES*AMA/GASC/TEMP            ! AIR DENSITY
C
      AKA   = (4.39+0.071*TEMP)*1D-3        ! AIR THERMAL CONDUCTIVITY
C
      DV    = (0.211D0/PRESA)*(TEMP/273D0)**1.94
      DV    = DV*1D-4                       ! WATER VAPOR DIFFUSIVITY IN AIR
      DBIG  = 5.0D-6
      DLOW  = 0.207683*((ACCOM)**(-0.33048))
      DLOW  = DLOW*1D-6
C
C DV AVERAGE
C
      COEF  = ((2*PI*AMW/(GASC*TEMP))**0.5)
C
      DV    = (DV/(DBIG-DLOW))*((DBIG-DLOW)-(2*DV/ACCOM)*COEF*
     &        (DLOG((DBIG+(2*DV/ACCOM)*COEF)/(DLOW+(2*DV/ACCOM)*
     &        COEF))))                      ! NON-CONTINUUM EFFECTS
C
C      PSAT  = VPRES(SNGL(TEMP))*(1E5/1.0D3) ! SATURATION VAPOR PRESSURE
      psat = wv_psat( temp, lhe ) * 100.d0 ! [Pa], jan perlwitz, Nov 2017
C
      SURT  = SFT(SNGL(TEMP))       ! SURFACE TENSION FOR WATER (J M-2)
C
      RETURN
C
C *** END OF SUBROUTINE PROPS *******************************************
C
      END



C======================================================================C
C *** FUNCTION VPRES
C *** THIS FUNCTION CALCULATES SATURATED WATER VAPOUR PRESSURE AS A
C     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -50 AND
C     50 C.
C
C ======================== ARGUMENTS / USAGE ==========================C
C  INPUT:
C     [T]
C     REAL VARIABLE.
C     AMBIENT TEMPERATURE EXPRESSED IN KELVIN.
C
C  OUTPUT:
C     [VPRES]
C     REAL VARIABLE.
C     SATURATED VAPOR PRESSURE EXPRESSED IN MBAR.
C
C======================================================================C
      REAL FUNCTION VPRES (T)
      REAL A(0:6), T
      DATA A/6.107799610E+0, 4.436518521E-1, 1.428945805E-2,
     &       2.650648471E-4, 3.031240396E-6, 2.034080948E-8,
     &       6.136820929E-11/
C
C CALCULATE POLYNOMIAL (WITHOUT EXPONENTIATION).
C
      TTEMP = T-273
      VPRES = A(6)*TTEMP
      DO I=5,1,-1
         VPRES = (VPRES + A(I))*TTEMP
      ENDDO
      VPRES = VPRES + A(0)
C
C END OF FUNCTION VPRES
C
      RETURN
      END



C======================================================================C
C *** FUNCTION SFT
C *** THIS FUNCTION CALCULATES WATER SURFACE TENSION AS A
C     FUNCTION OF TEMPERATURE. VALID FOR TEMPERATURES BETWEEN -40 AND
C     40 C.
C
C ======================== ARGUMENTS / USAGE ==========================C
C  INPUT:
C     [T]
C     REAL VARIABLE.
C     AMBIENT TEMPERATURE EXPRESSED IN KELVIN.
C
C  OUTPUT:
C     [SFT]
C     REAL VARIABLE.
C     SURFACE TENSION EXPRESSED IN J M-2.
C
C======================================================================C
      REAL FUNCTION SFT (T)
      REAL T
C
      TPARS = T-273
      SFT   = 0.0761-1.55E-4*TPARS
C
      RETURN
      END


C ***********************************************************************
C
      SUBROUTINE GAULEG (X,W,N)
C
C CALCULATION OF POINTS AND WEIGHTS FOR N POINT GAUSS INTEGRATION
C ***********************************************************************
      REAL*8 X,W,EPS,X1,X2
      DIMENSION X(N), W(N)
      PARAMETER (EPS=1.E-6)
      PARAMETER (X1=-1.0, X2=1.0)
C
C CALCULATION
C
      M=(N+1)/2D0
      XM=0.5D0*(X2+X1)
      XL=0.5D0*(X2-X1)
      DO 12 I=1,M
        Z=COS(3.141592654D0*(I-.25D0)/(N+.5D0))
1       CONTINUE
          P1=1.D0
          P2=0.D0
          DO 11 J=1,N
            P3=P2
            P2=P1
            P1=((2.D0*J-1.)*Z*P2-(J-1.D0)*P3)/J
11        CONTINUE
          PP=N*(Z*P1-P2)/(Z*Z-1.D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS)GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.D0*XL/((1.D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
12    CONTINUE
      RETURN
      END



C
C======================================================================C
      FUNCTION GAMMLN(XX)
C
C======================================================================C
      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)
      RETURN
      END


      
C======================================================================C
C      MATHEMATICAL ROUTINES
C
C======================================================================C
      SUBROUTINE INCOG(A,X,GIN,GIM,GIP)
C
C       ==================================================
C       PURPOSE: COMPUTE THE INCOMPLETE GAMMA FUNCTION
C                R(A,X), ?(A,X) AND P(A,X)
C       INPUT :  A   --- PARAMETER ( A ? 170 )
C                X   --- ARGUMENT
C       OUTPUT:  GIN --- R(A,X)
C                GIM --- ?(A,X)
C                GIP --- P(A,X)
C       ROUTINE CALLED: GAMMA FOR COMPUTING ?(X)
C       ==================================================C
!        IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT NONE
      REAL*8 XAM,X,A, GIN,GA,GIM,GIP,S,R,T0
      INTEGER K

        XAM=-X+A*DLOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'A AND/OR X TOO LARGE'
           STOP
        ENDIF
        IF (X.EQ.0.0) THEN
           GIN=0.0
           CALL GAMMA(A,GA)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           CALL GAMMA(A,GA)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           CALL GAMMA(A,GA)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
        END


        SUBROUTINE GAMMA(X,GA)
C
C       =================================================C       PURPOSE: COMPUTE GAMMA FUNCTION ?(X)
C       INPUT :  X  --- ARGUMENT OF ?(X)
C                       ( X IS NOT EQUAL TO 0,-1,-2,???)
C       OUTPUT:  GA --- ?(X)
C       =================================================C
!        IMPLICIT REAL*8 (A-H,O-Z)

        USE CONSTANT, only : PI
        IMPLICIT NONE
        REAL*8 X,GA,Z,R,GR,G
        INTEGER K,M1,M

        DIMENSION G(26)
!cyhl        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


C=======================================================================
C
C *** SUBROUTINE getCCN
C *** WRITTEN BY ATHANASIOS NENES
C *** GISS GCM (TOMAS microphys) interface
C *** This subroutine gets the aerosol microphysics properties from TOMAS.
C
C=======================================================================
C
      SUBROUTINE getCCN(TRM,BOXVi,Tot_MLi,
     &     Tt,TPi,MLi,NMxi,NSECi)
C
      use OldTracer_mod, only: tr_mm
      USE TRACER_COM, only : ntm
     &     ,n_ASO4,n_ANACL,n_AECIL,
     &     n_AOCIL,n_ANUM,nbins
      implicit none 
      INTEGER KK,nseci,k,NMxi
      real*8, dimension(ntm) :: trm
      REAL*8 TPi(NMxi), MLi(NMxi), SLFi, BOXVi, Mso4, Mna,
     &     Mocil,Tt,VhfSo4, VhfNacl,Vhfocil, Tot_MLi



      NSECi = nbins       !  Get number of sections
      SLFi  = 0d0         !  Total sulfate (ug m-3)
      Tt    = 0d0         !  Total particles (per m3)
C
      DO k=1,NSECi
         TPi(k) = TRM(n_ANUM(1)-1+k)   ! Number of particles (total)
         TPi(k) = TPi(k)/BOXVi                      ! Concentration (per m3)
         if (NSECi.eq.15.and.k.gt.3)THEN
           Tt=Tt + TPi(k)       ! Total concentration (per m3)
C
         ELSE !should be changed to lower size bin. 
           Tt=Tt+Tpi(k)
!           PRINT*,'getCCN in',k,TRM(n_ANUM(1)-1+k),Tt,BOXVi
         endif

         Mso4   = TRM(n_ASO4(1)-1+k) ! vhf is 2.5 - ammonium bisulfate
         SLFi   = SLFi + Mso4                       ! accumulate sulfate
         Mso4   = Mso4/(tr_mm(n_ASO4(1))*1e-3)/BOXVi   ! moles of SO4
         VhfSo4 = 2.5                               ! vhf is 2.5 - ammonium bisulfate
C
         Mna    = TRM(n_ANACL(1)-1+k)     ! vhf is 2.0
         Mna    = Mna/(tr_mm(n_ANACL(1))*1e-3)/BOXVi     ! moles of NaCl
         VhfNacl= 2                                 ! vhf is 2 - sodium chloride

         Mocil    = TRM(n_AOCIL(1)-1+k)     ! vhf is 2.0
         Mocil    = Mocil/(tr_mm(n_AOCIL(1))*1e-3)/BOXVi     ! moles of OCIL
         Vhfocil= 1.43                               ! vhf is 1 - hydrophilic carbon 
cyhl  Vhfocil is followed by what Jeff used in a look-up table.  I am not sure this is good value. 

C
         MLi(k) = VhfSo4*Mso4 + VhfNacl*Mna+Vhfocil*Mocil         ! Total moles in solution, per m3
         MLi(k) = MAX(MLi(k),1d-20)                 ! Avoid numerical exception errors

         if (NSECi.eq.15.and.k.gt.3) then
           Tot_MLi=Tot_MLi+MLi(k) !total solute molar concentration for B-L relationship. YHL (2008)
         elseif(nseci.ne.15)then
           Tot_MLi=Tot_MLi+MLi(k)
         endif

      ENDDO
C
      SLFi = 1d6*SLFi*1000./BOXVi                   ! Sulfate concentration, ug m-3
C
      RETURN
      END
