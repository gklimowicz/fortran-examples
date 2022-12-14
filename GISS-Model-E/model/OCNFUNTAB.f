C****   
C**** FUNTABLE.OCN   Table functions of ocean parameters   10/09/91
C****
      MODULE OCFUNC
!@sum  OCFUNC contains the ocean function lookup tables
!@auth Gary Russell/Gavin Schmidt      
      IMPLICIT NONE
      SAVE
C**** All look up tables take specific pot. enthalpy (J/kg), 
C**** salinity (ppt) and pressure (MPa = 10^6 Pa) as arguments
!@var VGSP specific volume (m^3/kg)
!@var TGSP potential temperature (C)
!@var HGSP specific enthalpy of seawater (J/kg)
!@var AGSP thermal expansion coefficient (kg/m**3/C)
!@var BGSP saline expansion coefficient (kg/m**3/PSU)
!@var CGS specific heat capacity of sea water (J/kg*C)
      REAL*8, DIMENSION(-2:40,0:40,0:39) :: VGSP,TGSP,HGSP,AGSP,BGSP
      REAL*8, DIMENSION(-2:40,0:40) :: CGS
      END MODULE OCFUNC

      REAL*8 FUNCTION VOLGSP (G,S,P)
C****
C**** VOLGSP returns a linearly interpolated specific volume from
C**** an input table that depends on potential specific enthalpy,
C**** salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: VOLGSP (m**3/kg) = specific volume of sea water
C****
      USE OCFUNC, only : v=>vgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      VOLGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*V(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*V(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*V(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*V(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION VOLGS (G,S)
C****
C**** VOLGS returns a linearly interpolated specific volume from
C**** an input table that depends on potential specific enthalpy
C**** and salinity.  Pressure is assumed to be the normal surface
C**** ocean reference pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: VOLGS (m**3/kg) = specific volume of sea water
C****
      USE OCFUNC, only : v=>vgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S
      REAL*8 GG,SS
      INTEGER IG,JS
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      VOLGS = (JS-SS+1)*((IG-GG+1)*V(IG  ,JS  ,0)
     *                 + (GG-IG  )*V(IG+1,JS  ,0))
     *      + (SS-JS  )*((IG-GG+1)*V(IG  ,JS+1,0)
     *                 + (GG-IG  )*V(IG+1,JS+1,0))
      RETURN
      END

      REAL*8 FUNCTION TEMGS (G,S)
C****
C**** TEMGS returns a linearly interpolated temperature from an
C**** input table that depends on potential specific enthalpy and
C**** salinity.  Pressure is assumed to be the normal surface
C**** ocean reference pressure. 
C**** 
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: TEMGS (C) = temperature of sea water
C****
      USE OCFUNC, only : t=>tgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S
      REAL*8 GG,SS
      INTEGER IG,JS
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      TEMGS = (JS-SS+1)*((IG-GG+1)*T(IG  ,JS  ,0)
     *                 + (GG-IG  )*T(IG+1,JS  ,0))
     *      + (SS-JS  )*((IG-GG+1)*T(IG  ,JS+1,0)
     *                 + (GG-IG  )*T(IG+1,JS+1,0))
      RETURN
      END

      REAL*8 FUNCTION SHCGS (G,S)
C****
C**** SHCGS returns a linearly interpolated specific heat capacity
C**** from an input table that depends on potential specific enthalpy
C**** and salinity.  Pressure is assumed to be the normal surface
C**** ocean reference pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: SHCGS (J/kg*C) = specific heat capacity of sea water
C****
      USE OCFUNC, only : c=>cgs
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      SHCGS = (JS-SS+1)*((IG-GG+1)*C(IG  ,JS  )
     *                 + (GG-IG  )*C(IG+1,JS  ))
     *      + (SS-JS  )*((IG-GG+1)*C(IG  ,JS+1)
     *                 + (GG-IG  )*C(IG+1,JS+1))
      RETURN
      END

      REAL*8 FUNCTION GFREZS (S)
C****
C**** GFREZS returns a linearly interpolated freezing point of
C**** potential specific enthalpy that depends on salinity.  Pressure
C**** is assumed to be the normal surface ocean reference pressure.
C****
C**** Input: S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: GFREZS (J/kg) = freezing point of potential specific
C****                         enthalpy
C****
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: S
      REAL*8 F(0:40)
      DATA F /                                                0.000,
     * -232.482d0,  -461.136d0,  -688.288d0,  -914.774d0, -1141.078d0,
     *-1367.530d0, -1594.374d0, -1821.798d0, -2049.957d0, -2278.978d0,
     *-2508.971d0, -2740.030d0, -2972.237d0, -3205.667d0, -3440.386d0,
     *-3676.454d0, -3913.926d0, -4152.851d0, -4393.287d0, -4635.259d0,
     *-4878.815d0, -5123.994d0, -5370.830d0, -5619.358d0, -5869.609d0,
     *-6121.615d0, -6375.402d0, -6631.001d0, -6888.436d0, -7147.733d0,
     *-7408.917d0, -7672.010d0, -7937.037d0, -8204.019d0, -8472.976d0,
     *-8743.931d0, -9016.917d0, -9291.927d0, -9568.992d0, -9848.131d0/
      REAL*8 SS
      INTEGER JS
C****
      SS = S*1000.
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
C****
      GFREZS = (JS-SS+1)*F(JS) + (SS-JS)*F(JS+1)
      RETURN
      END

      REAL*8 FUNCTION TFREZS (SIN)
C****
C**** TFREZS calculates the freezing temperature of sea water as a
C**** function of salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: SIN (1) = salinity (kg NaCl/kg sea water), from .004 to .04
C****
C**** Output: TFREZS (C) = freezing temperature of sea water
C****
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: SIN
      REAL*8 :: A01 = -.0575d0, A02 = -2.154996D-4, A03 =1.710523D-3
      REAL*8 S,S32
C****
      S   = SIN*1.D3
      S32 = S*DSQRT(S)
      TFREZS = (A01 + A02*S)*S + A03*S32
      RETURN
      END

      REAL*8 FUNCTION HETGSP (G,S,P)
C****
C**** HETGSP returns a linearly interpolated specific enthalpy from
C**** an input table that depends on potential specific enthalpy,
C**** salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: HETGSP (J/kg) = specific enthalpy of sea water
C****
      USE OCFUNC, only : h=>hgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      HETGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*H(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*H(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*H(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*H(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*H(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*H(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*H(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*H(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION ALPHAGSP (G,S,P)
C****
C**** ALPHAGSP returns a linearly interpolated thermal expansion
C**** coefficient from an input table that depends on potential
C**** specific enthalpy, salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: ALPHAGSP (kg/m**3/C) = thermal expansion coefficient 
C****
      USE OCFUNC, only : a=>agsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      ALPHAGSP=(KP-PP+1)*((JS-SS+1)*((IG-GG+1)*A(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*A(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*A(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*A(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*A(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*A(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*A(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*A(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION BETAGSP (G,S,P)
C****
C**** BETAGSP returns a linearly interpolated saline expansion
C**** coefficient from an input table that depends on potential 
C**** specific enthalpy, salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: BETAGSP (kg/m**3/PSU) = saline expansion coefficient 
C****
      USE OCFUNC, only : b=>bgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      BETAGSP =(KP-PP+1)*((JS-SS+1)*((IG-GG+1)*B(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*B(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*B(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*B(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*B(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*B(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*B(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*B(IG+1,JS+1,KP+1)))
      RETURN
      END

      REAL*8 FUNCTION TEMGSP (G,S,P)
C****
C**** TEMGSP returns a linearly interpolated in situ temperature
C**** from an input table that depends on potential 
C**** specific enthalpy, salinity, and pressure.
C****
C**** Input: G (J/kg) = potential specific enthalpy,
C****                   from -8000 to 160000
C****           S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****          P (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 78.E6
C****
C**** Output: TEMGSP (C) = in situ temperature
C****
      USE OCFUNC, only : t=>tgsp
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: G,S,P
      REAL*8 GG,SS,PP
      INTEGER IG,JS,KP
C****
      GG = G*2.5d-4   ! /4000.
      SS = S*1000.
      PP = P*5d-7     ! /2.D6
      IG = INT(GG+2.) - 2
      IF(IG.LT.-2)  IG = -2
      IF(IG.GE.40)  IG = 39
      JS = SS
C     IF(JS.LT. 0)  JS =  0
      IF(JS.GE.40)  JS = 39
      KP = PP
      IF(KP.LT. 0)  KP =  0
      IF(KP.GE.39)  KP = 38
C****
      TEMGSP = (KP-PP+1)*((JS-SS+1)*((IG-GG+1)*T(IG  ,JS  ,KP  )
     *                             + (GG-IG  )*T(IG+1,JS  ,KP  ))
     *                  + (SS-JS  )*((IG-GG+1)*T(IG  ,JS+1,KP  )
     *                             + (GG-IG  )*T(IG+1,JS+1,KP  )))
     *       + (PP-KP  )*((JS-SS+1)*((IG-GG+1)*T(IG  ,JS  ,KP+1)
     *                             + (GG-IG  )*T(IG+1,JS  ,KP+1))
     *                  + (SS-JS  )*((IG-GG+1)*T(IG  ,JS+1,KP+1)
     *                             + (GG-IG  )*T(IG+1,JS+1,KP+1)))
      RETURN
      END

! All functions below this line were imported from the offline OIC
! program. Todo: remove implicit typing and gotos.

      Real*8 Function VOLPTS (PIN,T,SIN)
C****
C**** VOLPTS calculates the specific volume of sea water as a
C**** function of pressure, temperature and salinity.
C**** The reference for this function is:
C**** "Tenth report of the joint panel on oceanographic tables and
C**** standards", Sidney, British Columbia, Canada, 1-5 September
C**** 1980, sponsored by Unesco, ICES, SCOR, IAPSO.
C**** Also see:
C**** N.P. Fofonoff, 1985.  Physical Properties of Seawater: A New
C**** Salinity Scale and Equation of State for Seawater.  Journal
C**** of Geophysical Research, volume 90, pp 3332-3342.
C****
C**** Input: PIN (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 1.E8
C****          T (C)  = temperature, from -2 to 40
C****        SIN (1)  = salinity (kg NaCl/kg sea water), from 0 to .042
C****
C**** Output: VOLPTS (m^3/kg) = specific volume of sea water
C****
      Implicit None
      Real*8,Intent(In) :: PIN,T,SIN
      Real*8 :: A0,A1,A2,A3,A4,A5, B0,B1,B2,B3,B4, C0,C1,C2, D0,
     *          E0,E1,E2,E3,E4, F0,F1,F2,F3, G0,G1,G2, H0,H1,H2,H3,
     *          I0,I1,I2, J0, K0,K1,K2, M0,M1,M2,
     *          P,S,S32,KW,AW,BW,KO,A,B,K,DENSTW,DENST0
      Data A0,A1,A2,A3,A4,A5 /999.842594, 6.793952D-2,
     *  -9.095290D-3, 1.001685D-4, -1.120083D-6, 6.536332D-9/
      Data B0,B1,B2,B3,B4 /8.24493D-1, -4.0899D-3, 7.6438D-5,
     *  -8.2467D-7, 5.3875D-9/
      Data C0,C1,C2 /-5.72466D-3, 1.0227D-4, -1.6546D-6/
      Data D0 /4.8314D-4/
      Data E0,E1,E2,E3,E4 /19652.21, 148.4206, -2.327105,
     *  1.360477D-2, -5.155288D-5/
      Data F0,F1,F2,F3 /54.6746, -.603459, 1.09987D-2, -6.1670D-5/
      Data G0,G1,G2 /7.944D-2, 1.6483D-2, -5.3009D-4/
      Data H0,H1,H2,H3 /3.239908, 1.43713D-3, 1.16092D-4, -5.77905D-7/
      Data I0,I1,I2 /2.2838D-3, -1.0981D-5, -1.6078D-6/
      Data J0 /1.91075D-4/
      Data K0,K1,K2 /8.50935D-5, -6.12293D-6, 5.2787D-8/
      Data M0,M1,M2 /-9.9348D-7, 2.0816D-8, 9.1697D-10/
C****
      P  = PIN*1.D-5
      S  = SIN*1.D3
      S32= S*Sqrt(S)
      KW = E0+(E1+(E2+(E3+E4*T)*T)*T)*T
      AW = H0+(H1+(H2+H3*T)*T)*T
      BW = K0+(K1+K2*T)*T
      KO = KW + (F0+(F1+(F2+F3*T)*T)*T)*S + (G0+(G1+G2*T)*T)*S32
      A  = AW + (I0+(I1+I2*T)*T)*S + J0*S32
      B  = BW + (M0+(M1+M2*T)*T)*S
      K  = KO + A*P + B*P**2
      DENSTW = A0+(A1+(A2+(A3+(A4+A5*T)*T)*T)*T)*T
      DENST0 = DENSTW + (B0+(B1+(B2+(B3+B4*T)*T)*T)*T)*S
     *                + (C0+(C1+C2*T)*T)*S32
     *                +  D0*S**2
      VOLPTS = (1.-P/K)/DENST0
      Return
      End


      Real*8 Function SHCPTS (PIN,T,SIN)
C****
C**** SHCPTS calculates the specific heat capacity of sea water as
C**** a function of pressure, temperature and salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C**** Also see:
C**** N.P. Fofonoff, 1985.  Physical Properties of Seawater: A New
C**** Salinity Scale and Equation of State for Seawater.  Journal
C**** of Geophysical Research, volume 90, pp 3332-3342.
C****
C**** Input: PIN (Pa) = pressure above normal atmospheric pressure,
C****                   from 0 to 1.E8
C****          T (C)  = temperature, from 0 to 35
C****        SIN (1)  = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: SHCPTS (J/kg*C) = specific heat capacity of sea water
C****                           with standard deviation error of
C****                           .636 (J/C*kg)
C****
      Implicit None
      Real*8,Intent(In) :: PIN,T,SIN
      Real*8 :: A000,A001,A002, A010,A011,A012, A020,A021,A022, A030,
     *          A040, A100,A101,A102, A110,A111,A112, A120,A121,A122,
     *          A130,A131, A140,A141,
     *          A200,A201,A202, A210,A211, A220,A221, A230,A231, A240,
     *          A300,A301, A310,A311,A312, A320,A321, A330,
     *          P,S,S32
      Data A000/ 4217.4    /, A001/-7.643575  /, A002/  .1770383 /,
     *     A010/-3.720283  /, A011/  .1072763 /, A012/-4.07718D-3/,
     *     A020/  .1412855 /, A021/-1.38385D-3/, A022/ 5.148D-5  /,
     *     A030/-2.654387D-3/,
     *     A040/ 2.093236D-5/,
     *     A100/-4.9592D-1 /, A101/ 4.9247D-3 /, A102/-1.2331D-4 /,
     *     A110/ 1.45747D-2/, A111/-1.28315D-4/, A112/-1.517D-6  /,
     *     A120/-3.13885D-4/, A121/ 9.802D-7  /, A122/ 3.122D-8  /,
     *     A130/ 2.0357D-6 /, A131/ 2.5941D-8 /,
     *     A140/ 1.7168D-8 /, A141/-2.9179D-10/,
     *     A200/ 2.4931D-4 /, A201/-2.9558D-6 /, A202/ 9.971D-8  /,
     *     A210/-1.08645D-5/, A211/ 1.17054D-7/,
     *     A220/ 2.87533D-7/, A221/-2.3905D-9 /,
     *     A230/-4.0027D-9 /, A231/ 1.8448D-11/,
     *     A240/ 2.2956D-11/,
     *     A300/-5.422D-8  /, A301/ 5.540D-10 /,
     *     A310/ 2.6380D-9 /, A311/-1.7682D-11/, A312/-1.4300D-12/,
     *     A320/-6.5637D-11/, A321/ 3.513D-13 /,
     *     A330/ 6.136D-13 /
C****
      P  = PIN*1.D-5
      S  = SIN*1.D3
      S32= S*Sqrt(S)
      SHCPTS = A000+(A010+(A020+(A030+A040*T)*T)*T)*T
     *  +   S*(A001+(A011+ A021                 *T)*T)
     *  + S32*(A002+(A012+ A022                 *T)*T)
     *  +     (A100+(A110+(A120+(A130+A140*T)*T)*T)*T
     *  +   S*(A101+(A111+(A121+(A131+A141*T)*T)*T)*T)
     *  + S32*(A102+(A112+ A122                 *T)*T)
     *  +     (A200+(A210+(A220+(A230+A240*T)*T)*T)*T
     *  +   S*(A201+(A211+(A221+ A231        *T)*T)*T)
     *  + S32* A202
     *  +     (A300+(A310+(A320+ A330        *T)*T)*T
     *  +   S*(A301+(A311+ A321                 *T)*T)
     *  + S32*       A312                          *T )*P)*P)*P
      Return
      End


      Real*8 Function ATGPTS (PIN,T,SIN)
C****
C**** ATGPTS calculates the adiabatic lapse rate of sea water as
C**** a function of pressure, temperature and salinity.
C**** The reference for this function is:
C**** N.P. Fofonoff and R.C. Millard Jr., 1983.  Algorithms for
C**** Computation of Fundamental Properties of Seawater.  UNESCO
C**** Technical Papers in Marine Science, volume 44.
C****
C**** Input: PIN (Pa) = pressure above normal atmospheric pressure
C****          T (C)  = temperature
C****        SIN (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: ATGPTS (C/Pa) = adiabatic lapse tate of sea water, at
C****                         S = .035, error < .006 (C) when used
C****                         to calculate potential temperature
C****
      Implicit None
      Real*8,Intent(In) :: PIN,T,SIN
      Real*8 :: A000,A010,A020,A030, A001,A011, A100,A110,A120,A130,
     *          A101,A111, A200,A210,A220,
     *          P,S
      Data A000/ 3.5803D-5 /, A010/ 8.5258D-6 /, A020/-6.8360D-8 /,
     *                        A030/ 6.6228D-10/,
     *     A001/ 1.8932D-6 /, A011/-4.2393D-8 /,
     *     A100/ 1.8741D-8 /, A110/-6.7795D-10/, A120/ 8.7330D-12/,
     *                        A130/-5.4481D-14/,
     *     A101/-1.1351D-10/, A111/ 2.7759D-12/,
     *     A200/-4.6206D-13/, A210/ 1.8676D-14/, A220/-2.1687D-16/
C****
      P = PIN*1.D-4
      S = SIN*1.D3 - 35.
      ATGPTS = A000+(A010+(A020+A030*T)*T)*T
     *    + S*(A001+ A011*T)
     *    +   (A100+(A110+(A120+A130*T)*T)*T
     *    + S*(A101+ A111*T)
     *    +   (A200+(A210+ A220        *T)*T)*P)*P
      ATGPTS = ATGPTS*1.D-4
      Return
      End


      Real*8 Function PTPTS (P,T,S)
C****
C**** PTPTS calculates the potential temperature of sea water as
C**** a function of pressure, temperature and salinity.
C**** At pressures above 0, PTPTS solves the differential equation:
C**** dT/dP = ATG = TK/SHC * dSVOL/dT .
C**** At pressure = 0, PTPTS = T.
C****
C**** Input: P (Pa) = pressure above normal atmospheric pressure
C****        T (C)  = temperature
C****        S (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: PTPTS (C) = potential temperature of sea water,
C****                     with maximum error of .004 (C) ?
C****
      Implicit None
      Real*8,Intent(In) :: P,T,S
      Real*8,External   :: ATGPTS
      Real*8  :: DP,T0,P0,T1
      Integer :: N,NM
!****
      NM = 1 + Abs(P)/2.D6
      DP = P/NM
      T0 = T
      DO 10 N=NM,1,-1
      P0 = N*DP
      T1 = T0 - .5*DP*ATGPTS(P0-.25*DP,T0,S)
   10 T0 = T0 -    DP*ATGPTS(P0-.50*DP,T1,S)
      PTPTS = T0
      Return
      End


      Real*8 Function DELHTS (T,SIN)
C****
C**** DELHTS calculates the change of specific enthalpy of sea
C**** water as salinity changes from 0 to an input value, as a
C**** function of temperature at atmospheric pressure.
C**** The reference for this function is:
C**** Frank J. Millero and Wing H. Leung, 1976.  The Thermodynamics
C**** of Seawater at One Atmosphere.  American Journal of Science,
C**** volume 276.
C****
C**** Input: T (C) = temperature, from -2 to 40
C****        S (1) = salinity (kg NaCl/kg sea water), from 0 to .04
C****
C**** Output: DELHTS (J/kg) = change of specific heat of sea water
C****
      Implicit None
      Real*8,Intent(In) :: T,SIN
      Real*8 :: A01,A03,A02, A11,A13,A12, A21,A23,A22, A31,A33,A32,
     *          S
      Data A01/ 3.4086D-3/, A03/ 7.9350D-4/, A02/-4.7989D-4/,
     *     A11/-6.3798D-5/, A13/ 1.0760D-4/, A12/ 6.3787D-6/,
     *     A21/ 1.3877D-6/, A23/-6.3923D-7/, A22/-1.1647D-7/,
     *     A31/-1.0512D-8/, A33/ 8.60D-9  /, A32/ 5.717D-10/
C****
      S = SIN*1.D3
      DELHTS = (A01+(A11+(A21+A31*T)*T)*T
     *       + (A03+(A13+(A23+A33*T)*T)*T)*Sqrt(S)
     *       + (A02+(A12+(A22+A32*T)*T)*T)*S)*S*1.D3
      Return
      End


      Real*8 Function HETPTS (P,T,S)
C****
C**** HETPTS calculates the specific enthalpy of sea water as a
C**** function of pressure, temperature and salinity.
C**** At pressures above 0, HETPTS solves the differential equation:
C**** dH/dP = V .
C**** At pressure = 0 and salinities above 0,
C**** H(0,T,S) = H(0,T,0) + DELHTS(T,S)
C**** At pressure = 0 and salinity = 0, HETPTS solves the differential
C**** equation:  dH/dT = C .
C****
C**** Input: P (Pa) = pressure above normal atmospheric pressure
C****        T (C)  = temperature
C****        S (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: HETPTS (J/kg) = specific enthalpy of sea water
C****
      Implicit None
      Real*8,Intent(In) :: P,T,S
      Real*8,External   :: ATGPTS,VOLPTS,DELHTS,SHCPTS
      Real*8  :: H0,T0,DP,P1,TX,T1,DT
      Integer :: N,NM
C****
C**** Calculate H(P,T,S) - H(0,T,S) by integrating  dH/dP = V  at
C**** constant entropy and salinity
C****
      H0 = 0.
      T0 = T
      If (P==0)  GoTo 20
      NM = 1 + Abs(P)/2.D6
      DP = P/NM
C**** For first step, integrate T down 1/2 DP, and H down full DP
      P1 = DP*(NM-.5)
      TX = T0 - .25*DP*ATGPTS(P1+.375*DP,T0,S)
      T1 = T0 - .50*DP*ATGPTS(P1+.250*DP,TX,S)
      H0 = H0 +     DP*VOLPTS(P1,T1,S)
C**** For subsequent steps, integrate T and H down full DP
      Do 10 N=NM-1,1,-1
      P1 = DP*(N-.5)
      T0 = T1 - .5*DP*ATGPTS(P1+.75*DP,T1,S)
      T1 = T1 -    DP*ATGPTS(P1+.50*DP,T0,S)
   10 H0 = H0 +    DP*VOLPTS(P1,T1,S)
C**** For last step, integrate T down 1/2 DP
      TX = T1 - .25*DP*ATGPTS(.375*DP,T1,S)
      T0 = T1 - .50*DP*ATGPTS(.250*DP,TX,S)
C****
C**** Calculate H(P,T,S) - H(0,T,0) = H(P,T,S) - H(0,T,S) + DELHTS(T,S)
C****
   20 H0 = H0 + DELHTS(T0,S)
C****
C**** Calculate H(0,T,0) - H(0,0,0) by integrating  dH/dT = C  at
C**** constant pressure and salinity
C****
      NM = 1 + Abs(T0)
      DT = T0/NM
      Do 30 N=0,NM-1
      T1 = DT*(N+.5)
   30 H0 = H0 + DT*SHCPTS(0.D0,T1,0.D0)
C****
      HETPTS = H0
      Return
      End


      Real*8 Function PHPTS (P,T,S)
C****
C**** PHPTS calculates the potential specific enthalpy of sea water
C**** as a function of pressure, temperature and salinity.
C****
C**** Input: P (Pa) = pressure above normal atmospheric pressure
C****        T (C)  = temperature
C****        S (1)  = salinity (kg NaCl/kg sea water)
C****
C**** Output: PHPTS (J/kg) = potential specific enthalpy of sea water
C****
      Implicit None
      Real*8,Intent(In) :: P,T,S
      Real*8,External   :: PTPTS,HETPTS
      Real*8 :: A
      A = PTPTS(P,T,S)
      PHPTS = HETPTS(0.D0,A,S)
      Return
      End
