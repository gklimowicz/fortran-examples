#include "rundeck_opts.h"

c -----------------------------------------------------------------

      module AMP_Utilities_mod

      contains

      character*2 function aerosolkind(tracerName) 
      character(len=*), intent(in) :: tracerName
      aerosolkind = tracerName(7:8)
      end function aerosolkind

      end module AMP_Utilities_mod

      SUBROUTINE SETAMP(EXT,SCT,GCB,TAB)
!@sum Calculation of extinction, asymmetry and scattering for AMP Aerosols
!@sum Calculation of absorption in the longwave
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer 
      USE domain_decomp_atm,ONLY: am_i_root

      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA,
     +                       AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS,
     +                       Reff_LEV, NUMB_LEV, RindexAMP, AMP_Q55, dry_Vf_LEV,
     +                       MIX_OC, MIX_SU, MIX_AQ, AMP_RAD_KEY
      USE AERO_CONFIG, only: NMODES
      USE AERO_PARAM,  only: DG_AKK,DG_ACC,DG_DD1,DG_DS1,DG_DD2, 
     +                       DG_DS2,DG_SSA,DG_SSC,DG_OCC,DG_BC1,
     +                       DG_BC2,DG_BC3,DG_DBC,DG_BOC,DG_BCS,DG_MXX

      USE RESOLUTION,  only: lm
      USE MODEL_COM,   only: itime,itimeI
      USE TRACER_COM,  only: TRM
      USE RADPAR,      only: aesqex,aesqsc,aesqcb,FSTOPX,FTTOPX,AMP_TAB_SPEC

      IMPLICIT NONE
      INTEGER, save:: Ifirstrad = 1
      ! Arguments: Optical Parameters dimension(lm,wavelength)
      REAL(8), INTENT(OUT) :: EXT(LM,6)       ! Extinction, SW
      REAL(8), INTENT(OUT) :: SCT(LM,6)       ! Single Scattering Albedo, SW
      REAL(8), INTENT(OUT) :: GCB(LM,6)       ! Asymmetry Factor, SW
      REAL(8), INTENT(OUT) :: TAB(LM,33)      ! Thermal absorption Cross section, LW
      REAL(8), DIMENSION(LM,NMODES) :: TTAUSV

      ! Local
      
      INTEGER l,n,w,MA,MB,MC,MD,NA,NS
      REAL*8 sizebins(23), Mie_IM(17), Mie_RE(15), HELP, AMP_TAB(33) 
      REAL*8 CORE_CLASS(nmodes), SHELL_CLASS(nmodes),Reff_mode(nmodes),Vf(6),CS_Mix(26)
      REAL*8 a,b,AMPEXT,AMPSCA,AMPASY
      DATA sizebins/0.002, 0.005,0.01,0.05,0.08,0.1,0.13,0.17,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2,1.5,2.,3.,5.,10./
      DATA CS_Mix/0.,0.04,0.08,0.12,0.16,0.2,0.24,0.28
     +          ,0.32,0.36,0.4,0.44,0.48,0.52,0.56,0.6
     +          ,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96,1.0/
      DATA Mie_RE/1.25,1.3,1.35,1.4,1.45,1.5,1.55,1.6,1.65,1.7,1.75,1.8,1.85,1.9,1.9/
      DATA Mie_IM/0.0,0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5,1.0/
c                        AKK  ACC  DD1  DS1  DD2  DS2  SSA  SSC  OCC  BC1  BC2  BC3  DBC  BOC  BCS  MXX
c                        1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16
      DATA CORE_CLASS   /1,   1,   6,   6,   6,   6,   2,   2,   4,   5,   5,   5,   6,   4,   5,   6/
      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0/
c      DATA SHELL_CLASS  /0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   1,   1,   2/
      DATA REFF_mode / 0.026D+00, 0.075D+00, 1.160D+00, 2.000D+00, 1.260D+00,
     +                 2.00D+00 , 0.12D+00 , 2.D+00   , 0.075D+00, 0.050D+00,  
     +                 0.100D+00, 0.100D+00, 0.330D+00, 0.100D+00, 0.070D+00, 0.100D+00/    
  
      CHARACTER*3 :: MODE_NAME(nmodes)=(/'AKK','ACC','DD1','DS1','DD2',
     +                             'DS2','SSA','SSC','OCC','BC1','BC2',
     +                                   'BC3','DBC','BOC','BCS','MXX'/)
c                   NA1= SO4  NA2=SS  NA3=NO3 NA4=OC NA5=BC NA6=DU

      EXT(:,:)    = 0.d0
      SCT(:,:)    = 0.d0
      GCB(:,:)    = 0.d0
      TAB(:,:)    = 0.d0
      TTAUSV(:,:) = 0.d0
C Longwave Pre calculate TAB: ---------------------------------------------------------------------------------------------

      if ( Ifirstrad==1 ) then
      Ifirstrad = 0
      DO n = 1,nmodes
      NA = CORE_CLASS(n)
      NS = SHELL_CLASS(n)
      Vf(:)=0.d0 
      CALL GET_LW(NA,NS,Reff_mode(n),AMP_TAB,Vf)
      AMP_TAB_SPEC(:,n)=AMP_TAB(:)
      enddo
      endif

      if (itime.ne.itimeI) then 
          IF (AMP_RAD_KEY == 1 .or. AMP_RAD_KEY ==3) THEN
         
c Shortwave: ---------------------------------------------------------------------------------------------    

      DO l = 1,lm
      DO n = 1,nmodes

         w = 6    ! aot at 550
         do MD = 1,23
            if (Reff_LEV(l,n) .le. sizebins(md)) goto 100            
         enddo
 100      continue  
          if (MD.gt.1) then
            MD = min(23,MD)
          b = sizebins(md) - Reff_LEV(l,n)
          a = Reff_LEV(l,n) - sizebins(md-1)
          endif
c---- INTERNAL MIXTURE ---------------------------------------------        
            do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 200
            enddo
 200        continue
            do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 201
            enddo
 201        continue
            MA = min(15,MA)              
            MB = min(17,MB)
            if (MD.gt.1) then
          TTAUSV(l,n) = NUMB_LEV(l,n) * (a/(b+a)* AMP_Q55(MA,MB,MD) + b/(a+b) * AMP_Q55(MA,MB,MD-1))
            else
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55(MA,MB,MD)
            endif
C----------------------------------------------------------------------
      DO w = 1,6  !wavelength
c---- INTERNAL MIXTURE ---------------------------------------------        
         do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 401
         enddo
 401     continue
         do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 402
         enddo
 402     continue
          MA = min(15,MA)
          MB = min(17,MB)

          if (MD.gt.1) then
          AMPEXT = (a/(b+a)* AMP_EXT(MA,MB,MD,w)  + b/(a+b) * AMP_EXT(MA,MB,MD-1,w))
          AMPSCA = (a/(b+a)* AMP_SCA(MA,MB,MD,w)  + b/(a+b) * AMP_SCA(MA,MB,MD-1,w))
          AMPASY = (a/(b+a)* AMP_ASY(MA,MB,MD,w)  + b/(a+b) * AMP_ASY(MA,MB,MD-1,w))
          else
          AMPEXT = AMP_EXT(MA,MB,MD,w) 
          AMPSCA = AMP_SCA(MA,MB,MD,w) 
          AMPASY = AMP_ASY(MA,MB,MD,w) 
          endif
c--------------------------------------------------------------------
          EXT(l,w) = EXT(l,w) + ( AMPEXT * TTAUSV(l,n) * FSTOPX(n))
          HELP     = ((GCB(l,w) * SCT(l,w) ) + (AMPASY * AMPSCA * TTAUSV(l,n)* FSTOPX(n))) 
          SCT(l,w) = SCT(l,w) +  (AMPSCA * TTAUSV(l,n) * FSTOPX(n))
          GCB(l,w) = HELP / (SCT(l,w)+ 1.D-10)
 
          aesqex(l,w,n)= AMPEXT * TTAUSV(l,n) 
          aesqsc(l,w,n)= AMPSCA * TTAUSV(l,n) 
          aesqcb(l,w,n)= AMPASY * aesqsc(l,w,n)

      ENDDO   ! wave
      ENDDO   ! modes
      ENDDO   ! level

         ENDIF       ! AMP_RAD_KEY=1or3

c --------------------------------------------------------------------------------------------------------    
c --------------------------------------------------------------------------------------------------------    

         IF (AMP_RAD_KEY == 2) THEN
c Shortwave: ---------------------------------------------------------------------------------------------    

      DO l = 1,lm
      DO n = 1,nmodes

         w = 6    ! aot at 550
         do MD = 1,23
            if (Reff_LEV(l,n) .le. sizebins(md)) goto 500            
         enddo
 500      continue  
          if (MD.gt.1) then
            MD = min(23,MD)
          b = sizebins(md) - Reff_LEV(l,n)
          a = Reff_LEV(l,n) - sizebins(md-1)
          endif

       select case (MODE_NAME(n))
c---- INTERNAL MIXTURE ---------------------------------------------        
       case ('AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','DBC','MXX')
   
            do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 600
            enddo
 600        continue
            do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 601
            enddo
 601        continue
            MA = min(15,MA)              
            MB = min(17,MB)
            if (MD.gt.1) then
          TTAUSV(l,n) = NUMB_LEV(l,n) * (a/(b+a)* AMP_Q55(MA,MB,MD) + b/(a+b) * AMP_Q55(MA,MB,MD-1))
            else
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55(MA,MB,MD)
            endif

C------ CORE SHELL -------------------------------------------------
       case ('BC1','BC2','BC3','BOC','BCS')
   
         do MA = 1,26
            if (MIX_OC(l,n) .le. CS_MIX(MA)) goto 701
         enddo
 701     continue
         do MB = 1,26
            if (MIX_SU(l,n) .le. CS_MIX(MB)) goto 702
         enddo
 702     continue
         do MC = 1,26
            if (MIX_AQ(l,n) .le. CS_MIX(MC)) goto 703
         enddo
 703     continue
         
          MA = min(26,MA-1)
          MB = min(26,MB-1)
          MC = min(26,MC-1)
          MA = max(1,MA-1)
          MB = max(1,MB-1)
          MC = max(1,MC-1)

         if (MD.gt.1) then
          TTAUSV(l,n) = NUMB_LEV(l,n) * (a/(b+a)* AMP_Q55_CS(MD,MA,MB,MC) + b/(a+b) * AMP_Q55_CS(MD-1,MA,MB,MC))
         else
          TTAUSV(l,n) = NUMB_LEV(l,n) * AMP_Q55_CS(MD,MA,MB,MC)
         endif
       end select
C----------------------------------------------------------------------
      DO w = 1,6  !wavelength
c---- INTERNAL MIXTURE ---------------------------------------------        

       select case (MODE_NAME(n))

       case ('AKK','ACC','DD1','DS1','DD2','DS2','SSA','SSC','OCC','DBC','MXX')
   
         do MA = 1,15
            if ( real    (RindexAMP(l,n,w)) .le. Mie_RE(MA)) goto 801
         enddo
 801     continue
         do MB = 1,17
            if ( aimag   (RindexAMP(l,n,w)) .le. Mie_IM(MB)) goto 802
         enddo
 802     continue
          MA = min(15,MA)
          MB = min(17,MB)

          if (MD.gt.1) then
          AMPEXT = (a/(b+a)* AMP_EXT(MA,MB,MD,w)  + b/(a+b) * AMP_EXT(MA,MB,MD-1,w))
          AMPSCA = (a/(b+a)* AMP_SCA(MA,MB,MD,w)  + b/(a+b) * AMP_SCA(MA,MB,MD-1,w))
          AMPASY = (a/(b+a)* AMP_ASY(MA,MB,MD,w)  + b/(a+b) * AMP_ASY(MA,MB,MD-1,w))
          else
          AMPEXT = AMP_EXT(MA,MB,MD,w) 
          AMPSCA = AMP_SCA(MA,MB,MD,w) 
          AMPASY = AMP_ASY(MA,MB,MD,w) 
          endif
C------ CORE SHELL -------------------------------------------------
       case ('BC1','BC2','BC3','BOC','BCS')
   
         do MA = 1,26
            if (MIX_OC(l,n) .le. CS_MIX(MA)) goto 805
         enddo
 805     continue
         do MB = 1,26
            if (MIX_SU(l,n) .le. CS_MIX(MB)) goto 806
         enddo
 806     continue
         do MC = 1,26
            if (MIX_AQ(l,n) .le. CS_MIX(MC)) goto 807
         enddo
 807     continue

          MA = min(26,MA-1)
          MB = min(26,MB-1)
          MC = min(26,MC-1)
          MA = max(1,MA-1)
          MB = max(1,MB-1)
          MC = max(1,MC-1)
                    
          if (MD.gt.1) then
          AMPEXT = (a/(b+a)* AMP_EXT_CS(MD,MA,MB,MC,w)  + b/(a+b) * AMP_EXT_CS(MD-1,MA,MB,MC,w))
          AMPSCA = (a/(b+a)* AMP_SCA_CS(MD,MA,MB,MC,w)  + b/(a+b) * AMP_SCA_CS(MD-1,MA,MB,MC,w))
          AMPASY = (a/(b+a)* AMP_ASY_CS(MD,MA,MB,MC,w)  + b/(a+b) * AMP_ASY_CS(MD-1,MA,MB,MC,w))
          else
          AMPEXT = AMP_EXT_CS(MD,MA,MB,MC,w)
          AMPSCA = AMP_SCA_CS(MD,MA,MB,MC,w)
          AMPASY = AMP_ASY_CS(MD,MA,MB,MC,w)
          endif
        end select
c--------------------------------------------------------------------
          EXT(l,w) = EXT(l,w) + ( AMPEXT * TTAUSV(l,n) * FSTOPX(n))
          HELP     = ((GCB(l,w) * SCT(l,w) ) + (AMPASY * AMPSCA * TTAUSV(l,n)* FSTOPX(n))) 
          SCT(l,w) = SCT(l,w) +  (AMPSCA * TTAUSV(l,n) * FSTOPX(n))
          GCB(l,w) = HELP / (SCT(l,w)+ 1.D-10)
 
          aesqex(l,w,n)= AMPEXT * TTAUSV(l,n) 
          aesqsc(l,w,n)= AMPSCA * TTAUSV(l,n) 
          aesqcb(l,w,n)= AMPASY * aesqsc(l,w,n)

      ENDDO   ! wave
      ENDDO   ! modes
      ENDDO   ! level

        ENDIF     ! AMP_RAD_KEY = 2


C Longwave: ---------------------------------------------------------------------------------------------

      DO l = 1,lm
      DO n = 1,nmodes
         TAB(l,:) = TAB(l,:) + (AMP_TAB_SPEC(:,n) *  TTAUSV(l,n) * FTTOPX(n))
      ENDDO   ! modes
      ENDDO   ! level
      endif
  
      RETURN
      END SUBROUTINE SETAMP
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE SETAMP_LEV(i,j,l)
!@sum Calulates effective Radius and Refractive Index for Mixed Aerosols
!@sum Puts AMP Aerosols in 1 dimension CALLED in RADIA
!@auth Susanne Bauer

      USE AMP_AEROSOL, only: DIAM,Reff_LEV, NUMB_LEV, RindexAMP,
     +  dry_Vf_LEV,MIX_OC,MIX_SU,MIX_AQ,AMP_dens, AMP_RAD_KEY
      USE AmpTracersMetadata_mod,  only: AMP_AERO_MAP, AMP_NUMB_MAP,
     +  AMP_MODES_MAP
      USE TRACER_COM,  only: TRM, ntmAMPi,ntmAMPe
      USE AERO_CONFIG, only: NMODES
      USE AERO_SETUP,  only: SIG0, CONV_DPAM_TO_DGN   !(nmodes * npoints) lognormal parameters for each mode
      USE GEOM,        only: BYDXYP ! inverse area of gridbox [m-2]

      USE AERO_ACTV, only: DENS_SULF, DENS_DUST,DENS_SEAS, DENS_BCAR, DENS_OCAR
      use OldTracer_mod, only: trname
      use AMP_utilities_mod, only: aerosolkind
      IMPLICIT NONE

      ! Arguments: 
      INTEGER, INTENT(IN) :: i,j,l

      ! Local
      INTEGER n,w,s,nAMP
      REAL*8,     DIMENSION(nmodes,7) :: VolFrac, VMass=0
      REAL*8                          :: H2O, NO3 
      REAL(8), PARAMETER :: TINYNUMER = 1.0D-30 
      COMPLEX*8, DIMENSION(6,7)      :: Ri
c     Variables for Maxwell Garnett:
      REAL*8                           :: V_bc, V_host
      COMPLEX*8                        :: M_mg, M_bc, M_host
      CHARACTER*3 :: MODE_NAME(nmodes)=(/'AKK','ACC','DD1','DS1','DD2',
     +                             'DS2','SSA','SSC','OCC','BC1','BC2',
     +                                   'BC3','DBC','BOC','BCS','MXX'/)
c Andies data incl Solar weighting - integral over 6 radiation band
      DATA Ri/(1.46099,    0.0764233)  ,(1.48313,  0.000516502),    !Su
     &        (1.49719,  1.98240e-05)  ,(1.50793,  1.64469e-06),
     &        (1.52000,  1.00000e-07)  ,(1.52815,  1.00000e-07),

c     &        (1.80056,     0.605467)  ,(1.68622,     0.583112),    !Bc
c     &        (1.63586,     0.551897)  ,(1.59646,     0.515333), 
c     &        (1.57466,     0.484662)  ,(1.56485,     0.487992),
c only 550nm values
c     &        (1.85,     0.71)  ,(1.85,     0.71),    !Bc
c     &        (1.85,     0.71)  ,(1.85,     0.71), 
c     &        (1.85,     0.71)  ,(1.85,     0.71),
cBond + Berstroem, all wavelength
     &        (2.15,     1.05)  ,(1.98,     0.79),    !Bc
     &        (1.92,     0.75)  ,(1.86,     0.73), 
     &        (1.85,     0.71)  ,(1.85,     0.71),

     &        (1.46099,    0.0761930)  ,(1.48313,   0.00470000),    !Oc
     &        (1.49719,   0.00470000)  ,(1.50805,   0.00480693),
     &        (1.52000,   0.00540000)  ,(1.52775,    0.0144927), 

     &        (1.47978,    0.0211233)  ,(1.50719,   0.00584169),    !Du
     &        (1.51608,   0.00378434)  ,(1.52998,   0.00178703),
     &        (1.54000,  0.000800000)  ,(1.56448,   0.00221463),

     &        (1.46390,   0.00571719)  ,(1.45000,      0.00000),    !Ss
     &        (1.45000,      0.00000)  ,(1.45000,      0.00000),
     &        (1.45000,      0.00000)  ,(1.45000,      0.00000),
 
     &        (1.46099,  0.0764233)    ,(1.48313,  0.000516502),    !No3
     &        (1.49719,  1.98240e-05)  ,(1.50793,  1.64469e-06),
     &        (1.52000,  1.00000e-07)  ,(1.52815,  1.00000e-07),

     &        (1.26304,    0.0774872)  ,(1.31148,  0.000347758),    !H2O
     &        (1.32283,  0.000115835)  ,(1.32774,  3.67435e-06),
     &        (1.33059,  1.58222e-07)  ,(1.33447,  3.91074e-08)/


       ! + Effective Radius [um] per Mode = geometric mass mean radius
       DO n=1,nmodes
         Reff_LEV(l,n) = DIAM(i,j,l,n)*CONV_DPAM_TO_DGN(n)*exp(5.*(sig0(n)**(-2))/2.)* 0.5e6
       ENDDO

       ! + Mass and Number Concentration
       DO n=ntmAMPi,ntmAMPe
         nAMP=n-ntmAMPi+1
           if(trname(n) .eq.'M_NO3') NO3 =trm(i,j,l,n)
           if(trname(n) .eq.'M_H2O') H2O =trm(i,j,l,n)
           if(AMP_NUMB_MAP(nAMP).eq. 0) then  ! Volume fraction
             if (trname(n).ne.'M_NO3'.and.trname(n).ne.'M_H2O'.and.
     &           trname(n).ne.'M_NH4') then
               select case (aerosolKind(trname(n)))
               case ('SU')
                  VMass(AMP_MODES_MAP(nAMP),1) =trm(i,j,l,n)/DENS_SULF
               case ('BC')
                  VMass(AMP_MODES_MAP(nAMP),2) =trm(i,j,l,n)/DENS_BCAR
               case ('OC')
                  VMass(AMP_MODES_MAP(nAMP),3) =trm(i,j,l,n)/DENS_OCAR
               case ('DU')
                  VMass(AMP_MODES_MAP(nAMP),4) =trm(i,j,l,n)/DENS_DUST
               case ('SS')
                  VMass(AMP_MODES_MAP(nAMP),5) =trm(i,j,l,n)/DENS_SEAS
               end select
             endif
           else                           ! Number
!          [ - ]                        [#/gb]         [m-2]      
           NUMB_LEV(l,AMP_NUMB_MAP(nAMP)) =trm(i,j,l,n) * bydxyp(j)
          endif
       ENDDO

       ! + Volume Fraction
       DO n=1,nmodes  ![#/m2]         pi/4     [m2]
        NUMB_LEV(l,n) = NUMB_LEV(l,n)* 0.7853 * DIAM(i,j,l,n)**2
        ! NO3   
        VMass(n,6) = VMass(n,1) / Sum(VMass(:,1)) * NO3/ 1720.
        ! H2O
        VMass(n,7) = VMass(n,1) /(Sum(VMass(:,1)) + TINYNUMER)  * H2O /1000.
       ENDDO

      DO s=1,7  ! loop over species 
        DO n=1,nmodes           ! loop over modes
          Volfrac(n,s) = VMass(n,s) / (Sum(VMass(n,:)) + TINYNUMER)
          dry_Vf_LEV(l,n,s) = VMass(n,s) / (Sum(VMass(n,1:6)) + TINYNUMER)
      ! Core Shell Composition
          if (n.eq.14) then     ! BOC
            MIX_OC(l,n) = VMass(n,3) / (VMass(n,1) + VMass(n,2) + VMass(n,3) + VMass(n,7) + TINYNUMER)
            MIX_SU(l,n) = VMass(n,1) / (VMass(n,1) + VMass(n,2) + VMass(n,3) + VMass(n,7) + TINYNUMER)
            MIX_AQ(l,n) = VMass(n,7) / (VMass(n,1) + VMass(n,2) + VMass(n,3) + VMass(n,7) + TINYNUMER)
          endif
          if (n.eq.15) then     ! BCS
            MIX_OC(l,n) = 0.d0
            MIX_SU(l,n) = VMass(n,1) / (VMass(n,1) + VMass(n,2) + VMass(n,7) + TINYNUMER)
            MIX_AQ(l,n) = VMass(n,7) / (VMass(n,1) + VMass(n,2) + VMass(n,7) + TINYNUMER)
          endif
          if (n.ge.10.and.n.le.12) then ! BC123
            MIX_OC(l,n) = 0.d0
            MIX_SU(l,n) = VMass(n,1) / (VMass(n,1) + VMass(n,2) + VMass(n,7) + TINYNUMER)
            MIX_AQ(l,n) = VMass(n,7) / (VMass(n,1) + VMass(n,2) + VMass(n,7) + TINYNUMER)
          endif
        ENDDO
      ENDDO
 
      ! + Refractive Index of Aerosol mix per mode and wavelength
      
      RindexAMP(l,:,:) = 0.d0
      DO s=1,7                  ! loop over species 
        DO w=1,6                ! loop over wavelength
          DO n=1,nmodes         ! loop over modes
            RindexAMP(l,n,w) = RindexAMP(l,n,w) + ( Volfrac(n,s) * Ri(w,s))
          ENDDO
        ENDDO
      ENDDO

      if (AMP_RAD_KEY == 3) then      ! - - - Maxwell Garnett Mixing Rule
        DO w=1,6                ! loop over wavelength
          DO n=1,nmodes         ! loop over modes
       select case (MODE_NAME(n))
       case ('BC1','BC2','BC3','BOC','BCS')
             M_bc   = Ri(w,2)
             M_host = ( Volfrac(n,1) * Ri(w,1))
             DO s=3,7                  ! loop over species other that BC
             M_host = M_host + ( Volfrac(n,s) * Ri(w,s))
             ENDDO
             V_bc   = Volfrac(n,2)
             V_host = Volfrac(n,1)+Volfrac(n,3)+Volfrac(n,4)+Volfrac(n,5)+Volfrac(n,6)+Volfrac(n,7)
             M_mg = M_host**2  * (M_bc**2 + 2.d0 * M_host**2 + 2.d0 * V_bc * (M_bc    - M_host   ) ) 
     +                         / (M_bc**2 + 2.d0 * M_host**2 -        V_host*(M_bc**2 - M_host**2) )

             RindexAMP(l,n,w) = SQRT( M_mg)
       end select      
          ENDDO
        ENDDO
      endif

      RETURN
      END SUBROUTINE SETAMP_LEV
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE SETUP_RAD
!@sum Initialization for Radiation incl. Aerosol Microphysics
!@auth Susanne Bauer

      USE AMP_AEROSOL, only: AMP_EXT, AMP_ASY, AMP_SCA, AMP_Q55,
     +           AMP_EXT_CS, AMP_ASY_CS, AMP_SCA_CS, AMP_Q55_CS  
	  
	  IMPLICIT NONE
      include 'netcdf.inc'
      integer start(4),count(4),count3(3),status
      integer start2(5),count2(5),count32(4)
      integer ncid, id1, id2, id3, id4,ncid2

      real*4, DIMENSION(15,17,23,6) ::  ASY,SCA,EXT
      real*4, DIMENSION(15,17,23) ::    QEX
      real*4, DIMENSION(23,26,26,26,6) ::  CS_ASY,CS_SCA,CS_EXT
      real*4, DIMENSION(23,26,26,26) ::    CS_QEX
c -----------------------------------------------------------------
c   Opening of the files to be read: MIE TABLES
c -----------------------------------------------------------------
          status=NF_OPEN('AMP_MIE_TABLES',NCNOWRIT,ncid)
          status=NF_INQ_VARID(ncid,'ASYM',id1)
          status=NF_INQ_VARID(ncid,'QEXT',id2)
          status=NF_INQ_VARID(ncid,'QSCT',id3)
          status=NF_INQ_VARID(ncid,'Q55E',id4)
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start(1)=1
          start(2)=1
          start(3)=1
          start(4)=1

          count(1)=15
          count(2)=17
          count(3)=23
          count(4)=6
          count3(1)=15
          count3(2)=17
          count3(3)=23
 

          status=NF_GET_VARA_REAL(ncid,id1,start,count,ASY)
          status=NF_GET_VARA_REAL(ncid,id2,start,count,EXT)
          status=NF_GET_VARA_REAL(ncid,id3,start,count,SCA)
          status=NF_GET_VARA_REAL(ncid,id4,start,count3,QEX)

          status=NF_CLOSE('AMP_MIE_TABLES',NCNOWRIT,ncid)

          AMP_ASY = ASY 
          AMP_EXT = EXT 
          AMP_SCA = SCA 
          AMP_Q55 = QEX
c -----------------------------------------------------------------
c   Opening of the files to be read: Core Shell Mie Tables
c   Core is BC, Shell Material is OC, SO4 and H2O
c -----------------------------------------------------------------
          status=NF_OPEN('AMP_CORESHELL_TABLES',NCNOWRIT,ncid2)
          status=NF_INQ_VARID(ncid2,'CS_ASYM',id1)
          status=NF_INQ_VARID(ncid2,'CS_QEXT',id2)
          status=NF_INQ_VARID(ncid2,'CS_QSCT',id3)
          status=NF_INQ_VARID(ncid2,'CS_Q55E',id4)
c -----------------------------------------------------------------
c   read
c -----------------------------------------------------------------
          start2(1)=1
          start2(2)=1
          start2(3)=1
          start2(4)=1
          start2(5)=1

          count2(1)=23
          count2(2)=26
          count2(3)=26 
          count2(4)=26
          count2(5)=6
          count32(1)=23
          count32(2)=26
          count32(3)=26
          count32(4)=26
 

          status=NF_GET_VARA_REAL(ncid2,id1,start2,count2,CS_ASY)
          status=NF_GET_VARA_REAL(ncid2,id2,start2,count2,CS_EXT)
          status=NF_GET_VARA_REAL(ncid2,id3,start2,count2,CS_SCA)
          status=NF_GET_VARA_REAL(ncid2,id4,start2,count32,CS_QEX)

          status=NF_CLOSE('AMP_CORESHELL_TABLES',NCNOWRIT,ncid2)

          AMP_ASY_CS = CS_ASY 
          AMP_EXT_CS = CS_EXT 
          AMP_SCA_CS = CS_SCA 
          AMP_Q55_CS = CS_QEX
      RETURN
      END SUBROUTINE SETUP_RAD
c -----------------------------------------------------------------

c -----------------------------------------------------------------
      SUBROUTINE GET_LW(NA,NS,AREFF,TQAB,Vf)
!@sum Calculation of LW absorption for AMP aerosols
!@sum Called in SETAER / RCOMPX
!@auth Susanne Bauer

      USE RADPAR, only: TRUQEX, TRSQEX, TRDQEX, TRUQSC, TRSQSC, TRDQSC
     *                  , REFU22, REFS25, REFD25
      INTEGER, intent(IN) :: NA,NS
      REAL*8,  intent(in) :: areff,Vf(6)
      REAL*8   TQEX(33),TQSC(33),TQAB(33),TQEX_S(33),TQSC_S(33)
      REAL*8   QXAERN(25),QSAERN(25)
      REAL*8   wts,wta
      INTEGER  n0,k,n,nn

c CORE  
       IF(NA==0) THEN
         TQAB(:)= 0d0
       ENDIF
                                      !                               1   2   3   4
      IF(NA > 0 .and. NA < 5) THEN    !    NA : Aerosol compositions SO4,SEA,NO3,OC
        N0=0
        IF(NA==2) N0=22
        IF(NA==3) N0=44
        IF(NA==4) N0=88
        DO K=1,33
        DO N=1,22
        NN=N0+N
        QXAERN(N)=TRUQEX(K,NN)
        QSAERN(N)=TRUQSC(K,NN)
        ENDDO
        CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
        ENDDO
      ENDIF
                              
      IF(NA==5) THEN                  !   NA : Aerosol compositions BC
        DO K=1,33
        QXAERN(:)=TRSQEX(K,:)    ! 1:25
        QSAERN(:)=TRSQSC(K,:)    ! 1:25
        CALL SPLINE(REFS25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFS25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
        ENDDO
      ENDIF

                                      !                              6
      IF(NA==6) THEN                  !   NA : Aerosol composition DST
       DO K=1,33
        QXAERN(:)=TRDQEX(K,:)    ! 1:25
        QSAERN(:)=TRDQSC(K,:)    ! 1:25
        CALL SPLINE(REFD25,QXAERN,25,AREFF,TQEX(K),1.D0,1.D0,1)
        CALL SPLINE(REFD25,QSAERN,25,AREFF,TQSC(K),1.D0,1.D0,1)
        TQAB(K)=TQEX(K)-TQSC(K)
      ENDDO
      ENDIF

c SHELL
c         IF(NS > 0 .and. NS < 5) THEN    !    NS : Aerosol compositions SO4,SEA,NO3,OC
c         N0=0
c         IF(NS==2) N0=22
c         IF(NS==3) N0=44
c         IF(NS==4) N0=88


c         DO K=1,33
c         DO N=1,22
c         NN=N0+N
c         IF (NS==1) WTS=Vf(1)                      ! <- shell fraction of aerosol composition
c         IF (NS==2) WTS=Vf(5)                      ! <- shell fraction of aerosol composition
c         WTA=1.D0-WTS
c         QXAERN(N)=TRUQEX(K,NN)
c         QSAERN(N)=TRUQSC(K,NN)
c         ENDDO
c         CALL SPLINE(REFU22,QXAERN,22,AREFF,TQEX_S(K),1.D0,1.D0,1)
c         CALL SPLINE(REFU22,QSAERN,22,AREFF,TQSC_S(K),1.D0,1.D0,1)
c         TQAB(K)=(TQEX(K)*WTA + TQEX_S(K)*WTS)-(TQSC(K)*WTA + TQSC_S(K)*WTS)
c         ENDDO

c      ENDIF

      RETURN
      END SUBROUTINE GET_LW
c -----------------------------------------------------------------

