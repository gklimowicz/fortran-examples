      Subroutine WATER (SUBR)
!@sum  WATER writes global salt water reservoirs to unit 6
!@auth Gary L. Russell
!@ver  2018/06/06

!* Salt Water Budget (kg) in Model E

!* Variables defined on atmosphere grid are not preceeded by 'a'
!* Variables defined on ocean grid are preceeded by 'o'

!* The section "Without Fluxes" indicates the prognostic variables that
!* contribute to the salt water mass (kg) for each reservoir

!* Without Fluxes (kg):
!* AtmosQ  = CONSERV_WM   = ATMO = Q*MA*aXYP + WM*MA*aXYP
!* Ground  = CONSERV_WTG  = GRNW = ground formula
!* LakeIce = CONSERV_LMSI = LAKI = (SNOWI+ACE1I+MSI)*RSI*FLAKE*aXYP
!* LiqLake = CONSERV_LKM  = LIQL = MWL
!* LandIce = CONSERV_MLI  = LANI = (SNOWLI+ACE1LI+ACE2LI)*FLICE*aXYP
!* IceBerg = CONSERV_MICB = ICEB = MDWNIMP
!* SeaIce  = CONSERV_OMSI = SEAI = (SNOWI+ACE1I+MSI)*RSI*aFOCEAN*aXYP
!* LiqOcen = CONSERV_OMS  = LIQO = oMO*oXYP

!* After many subroutines, water mass has been removed from some reservoirs but not yet been added to the next reservoirs;
!* water mass is temporarily stored in a "flux" array.
!* After each subroutine, reservoirs which do not have unadded fluxes are not mentioned; for other reservoirs flux arrays are shown.
!* Other than 6 character SUBR designators below, no fluxes are needed for water mass accounting. 

!* MELTIO : ATM_DRV.f : Call MELT_SI(Ocn) : IceOcn%MELTI removed from SeaIce and attributed to LiqOcen
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* MELTIK : ATM_DRV.f : Call MELT_SI(Lak) : IceLak%MELTI removed from LakeIce and attributed to LiqLake
!* LiqLake = LIQL + IceLak%MELTI*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* CONDSE : ATM_DRV.f : Call CONDSE : PREC removed from AtmosQ and attributed to LandIce,Ground,Lakeice,LiqLake,SeaIce,LiqOcen
!* Ground  = GRNW + PREC*FEARTH*aXYP
!* LakeIce = LAKI + PREC*SI_Atm%RSI*FLAKE*aXYP
!* LiqLake = LIQL + IceLak%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*FLAKE*aXYP
!* LandIce = LANI + PREC*FLICE*aXYP
!* SeaIce  = SEAI + PREC*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*aFOCEAN*aXYP

!* PRECIO : MODELE.f : Call PRECIP_SI(Ocn) : PREC added to SeaIce; IceOcn%RUNPSI removed from SeaIce and attributed to LiqOcen  
!* Ground  = GRNW + PREC*FEARTH*aXYP
!* LakeIce = LAKI + PREC*SI_Atm%RSI*FLAKE*aXYP
!* LiqLake = LIQL + IceLak%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*FLAKE*aXYP
!* LandIce = LANI + PREC*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*aFOCEAN*aXYP + IceOcn%RUNPSI*SI_Atm%RSI*AtmIce%FWATER*aXYP

!* PRECOC : MODELE.f : Call PRECIP_OC : PREC and IceOcn%RUNPSI added to LiqOcen
!* Ground  = GRNW + PREC*FEARTH*aXYP
!* LakeIce = LAKI + PREC*SI_Atm%RSI*FLAKE*aXYP
!* LiqLake = LIQL + IceLak%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*FLAKE*aXYP
!* LandIce = LANI + PREC*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* PRECIK : SURFACE.f : Call PRECIP_SI(Lak) : PREC added to LakeIce; IceLak%RUNPSI removed from LakeIce and attributed to LiqLake  
!* Ground  = GRNW + PREC*FEARTH*aXYP
!* LiqLake = LIQL + IceLak%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*FLAKE*aXYP + IceLak%RUNPSI*SI_Atm%RSI*FLAKE*aXYP
!* LandIce = LANI + PREC*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* PRECGI : SURFACE.f : Call PRECIP_LI : PREC added to LandIce; AtmGla%RUNO removed from LandIce and attributed to LiqLake
!*                                       AtmGla%IMPLM removed from LandIce and attributed to IceBerg
!* Ground  = GRNW + PREC*FEARTH*aXYP
!* LiqLake = LIQL + IceLak%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*FLAKE*aXYP + IceLak%RUNPSI*SI_Atm%RSI*FLAKE*aXYP +AtmGla%RUNO*FLICE*aXYP
!* IceBerg = ICEB + AtmGla%IMPLM*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* IRRIG : SURFACE.f : Call IRRIG_LK : IRRIG_WATER_ACT from LiqLake and aquifer water are attributed to Ground  
!* Ground  = GRNW + PREC*FEARTH*aXYP + IRRIG_WATER_ACT
!* LiqLake = LIQL + IceLak%MELTI*aXYP + PREC*(1-SI_Atm%RSI)*FLAKE*aXYP + IceLak%RUNPSI*SI_Atm%RSI*FLAKE*aXYP +AtmGla%RUNO*FLICE*aXYP
!* IceBerg = ICEB + AtmGla%IMPLM*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* PRECLK : SURFACE.f : Call PRECIP_LK : IceLak%MELTI, PREC, IceLak%RUNPSI and RUNOLI added to LiqLake 
!* Ground  = GRNW + PREC*FEARTH*aXYP + IRRIG_WATER_ACT
!* IceBerg = ICEB + AtmGla%IMPLM*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* EARTH1 : SURFACE.f : Call EARTH : PREC added to Ground; IRRIG_WATER_ACT added to Ground; 
!*                                   AtmLnd%EVAPOR removed from Ground and attributed to AtmosQ;
!*                                   AtmLnd%RUNO removed from Ground and attributed to LiqLake 
!* AtmosQ  = ATMO + AtmLnd%EVAPOR*FEARTH*aXYP
!* LiqLake = LIQL + AtmLnd%RUNO*FEARTH*aXYP
!* IceBerg = ICEB + AtmGla%IMPLM*FLICE*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP

!* ATMDIF : ATM_DRV.f : Call ATM_DIFFUS : AtmLnd%EVAPOR added to AtmosQ from Ground;
!*                                        EVAPOR added to AtmosQ and attributed from LakeIce,LiqLake,LandIce,SeaIce,LiqOcen
!* LakeIce = LAKI - AtmIce%EVAPOR*SI_Atm%RSI*FLAKE*aXYP
!* LiqLake = LIQL + AtmLnd%RUNO*FEARTH*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*FLAKE*aXYP
!* LandIce = LANI - AtmGlas(1)%EVAPOR*FLICE*aXYP
!* IceBerg = ICEB + AtmGla%IMPLM*FLICE*aXYP
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP

!* GRNDGI : SURFACE.f : Call GROUND_LI : AtmGlas(1)%EVAPOR removed from LandIce;
!*                                       AtmGla%RUNO removed from LandIce and attributed to LiqLake
!*                                       additional AtmGla%IMPLM removed from LandIce and attributed to IceBerg 
!* LakeIce = LAKI - AtmIce%EVAPOR*SI_Atm%RSI*FLAKE*aXYP
!* LiqLake = LIQL + AtmLnd%RUNO*FEARTH*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*FLAKE*aXYP + AtmGla%RUNO*FLICE*aXYP
!* IceBerg = ICEB + AtmGla%IMPLM*FLICE*aXYP
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP

!* SDIAG3 : SURFACE.f : Call SURFACE_DIAG3 : AtmGla%IMPLM added to IceBerg
!* LakeIce = LAKI - AtmIce%EVAPOR*SI_Atm%RSI*FLAKE*aXYP
!* LiqLake = LIQL + AtmLnd%RUNO*FEARTH*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*FLAKE*aXYP + AtmGla%RUNO*FLICE*aXYP
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP

!* GRNDIK : SURFACE.f : Call GROUND_SI(Lak) : AtmIce%EVAPOR removed from LakeIce;
!*                                            IceLak%RUNOSI removed from LakeIce and attributed to LiqLake
!* LiqLake = LIQL + AtmLnd%RUNO*FEARTH*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*FLAKE*aXYP + AtmGla%RUNO*FLICE*aXYP +
!*                + IceLak%RUNOSI*SI_Atm%RSI*FLAKE*aXYP
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP

!* GRNDLK : SURFACE.f : Call GROUND_LK : AtmOcn%EVAPOR removed from LiqLake;
!*                                       AtmLnd%RUNO, AtmGla%RUNO, and IceLak%RUNOSI added to LiqLake;
!*                                       IceOcn%DMSI removed from LiqLake and attributed to LakeIce 
!* LakeIce = LAKI + IceOcn%DMSI(1)*(1-SI_Atm%RSI)*FLAKE*aXYP + IceOcn%DMSI(2)*SI_Atm%RSI*FLAKE*aXYP
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP

!* RIVERF : SURFACE.f : Call RIVERF : AtmOcn%FLOWO removed from LiqLake and attributed to LiqOcen
!* LakeIce = LAKI + IceOcn%DMSI(1)*(1-SI_Atm%RSI)*FLAKE*aXYP + IceOcn%DMSI(2)*SI_Atm%RSI*FLAKE*aXYP
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP + AtmOcn%FLOWO*aFOCEAN*aXYP

!* FORMIK : SURFACE.f : Call FORM_SI(Lak) : IceOcn%DMSI added to LakeIce
!* SeaIce  = SEAI - AtmIce%EVAPOR*SI_Atm%RSI*aFOCEAN*aXYP
!* LiqOcen = LIQO + IceOcn%MELTI - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP + AtmOcn%FLOWO*aFOCEAN*aXYP

!* GRNDIO : OCN_DRV.f : Call GROUND_SI(Ocn) : AtmIce%EVAPOR removed from SeaIce;
!*                                            IceOcn%RUNOSI removed from SeaIce and attributed to LiqOcen
!* LiqOcen = LIQO + IceOcn%MELTI*aXYP - AtmOcn%EVAPOR*(1-SI_Atm%RSI)*aFOCEAN*aXYP + AtmOcn%FLOWO*aFOCEAN*aXYP +
!*                + IceOcn%RUNOSI*SI_Atm%RSI*aFOCEAN*aXYP*aXYP

!* GRNDOC : OCNDYN2.f : Call GROUND_OC : IceOcn%MELTI, AtmOcn%FLOWO, and IceOcn%RUNOSI added to LiqOcen;
!*                                       AtmOcn%EVAPOR removed from LiqOcen; oDMSI removed from LiqOcen and attributed to SeaIce
!* SeaIce  = SEAI + oDMSI(1)*(1-oRSI)*oFOCEAN*oXYP + oDMSI(2)*oRSI*oFOCEAN*oXYP =

!* FORMIO : OCN_DRV.f : Call FORM_SI(Ocn) : oDMSI added to SeaIce

!**** Input: SUBR = 6 character string that labels ouput line and WATER recognizes flux arrays are to be included  

      Use CONSTANT,    Only: RHOW
      Use RESOLUTION,  Only: IMA=>IM, JMA=>JM
      Use MODEL_COM,   Only: DTSRC,ITIME,NDAY,IYEAR1,AMON
      Use GEOM,        Only: aXYP,AREAG, aIMAXJ=>IMAXJ
      Use LANDICE,     Only: ACCPDA,ACCPDG
      Use GHY_COM,     Only: FEARTH
      Use LAKES_COM,   Only: FLAKE, IceLak
      Use SEAICE_COM,  Only: SI_ATM,IceOcn
      Use FLUXES,      Only: AtmGlas,AtmGla,AtmOcn,AtmIce,AtmLnd, PREC, IRRIG_WATER_ACT, FLICE, aFOCEAN=>FOCEAN
      Use OFLUXES,     Only: oRSI,oDMSI
      Use OCEAN,       Only: IMO=>IM, JMO=>JM, oXYP, oIMAXJ=>IMAXJ, oFOCEAN=>FOCEAN
      Use DOMAIN_DECOMP_ATM, Only: aGrid=>Grid,AM_I_ROOT,GLOBALSUM
      Use OCEANR_DIM,        Only: oGrid
      Use TimeConstants_Mod, Only: EDperY=>Earth_Days_Per_Year

      Implicit  None
      Character*6,Intent(In) :: SUBR  !@var SUBR identifies after which subroutine WATER was called

!**** Local Variables
      Integer :: I,J,L,N, I1A,INA,J1A,JNA, I1O,INO,J1O,JNO, JYEAR,JMON,JDAY,JDATE,JHOUR
      Logical :: QMELTIO,QMELTIK,QCONDSE,QPRECIO, QPRECOC,QPRECIK,QPRECGI,QIRRIG,  QPRECLK,QEARTH1,QATMDIF,QGRNDGI, &
                 QSDIAG3,QGRNDIK,QGRNDLK, QRIVERF,QFORMIK,QGRNDIO,QGRNDOC, QFORMIO,  QLL,QROOT
      Real*8  :: A(aGrid%I_STRT_HALO:aGrid%I_STOP_HALO,aGrid%J_STRT_HALO:aGrid%J_STOP_HALO),AGLOBE, &
                 O(oGrid%I_STRT_HALO:oGrid%I_STOP_HALO,oGrid%J_STRT_HALO:oGrid%J_STOP_HALO),OGLOBE, &
                 ATMO,GRNW,IRRI,LAKI,LIQL,LANI,ICEB,SEAI,LIQO,TOTL

      QMELTIO = SUBR.eq.'MELTIO'  ;  QMELTIK = SUBR.eq.'MELTIK'  ;  QCONDSE = SUBR.eq.'CONDSE'  ;  QPRECIO = SUBR.eq.'PRECIO'
      QPRECOC = SUBR.eq.'PRECOC'  ;  QPRECIK = SUBR.eq.'PRECIK'  ;  QPRECGI = SUBR.eq.'PRECGI'  ;  QIRRIG  = SUBR.eq.'IRRIG '
      QPRECLK = SUBR.eq.'PRECLK'  ;  QEARTH1 = SUBR.eq.'EARTH1'  ;  QATMDIF = SUBR.eq.'ATMDIF'  ;  QGRNDGI = SUBR.eq.'GRNDGI'
      QSDIAG3 = SUBR.eq.'SDIAG3'  ;  QGRNDIK = SUBR.eq.'GRNDIK'  ;  QGRNDLK = SUBR.eq.'GRNDLK'  ;  QRIVERF = SUBR.eq.'RIVERF'
      QFORMIK = SUBR.eq.'FORMIK'  ;  QGRNDIO = SUBR.eq.'GRNDIO'  ;  QGRNDOC = SUBR.eq.'GRNDOC'  ;  QFORMIO = SUBR.eq.'FORMIO'
      QROOT   = AM_I_ROOT()

      Call GetDte (ITIME,NDAY,IYEAR1, JYEAR,JMON,JDAY,JDATE,JHOUR, AMON)

!**** Extract domain decomposition band parameters
      I1A = aGrid%I_STRT  ;  INA = aGrid%I_STOP  ;  J1A = aGrid%J_STRT  ;  JNA = aGrid%J_STOP
      I1O = oGrid%I_STRT  ;  INO = oGrid%I_STOP  ;  J1O = oGrid%J_STRT  ;  JNO = oGrid%J_STOP
      QLL = aGrid%I_STRT_HALO == I1A

!****
!**** Calculate atmospheric water mass ATMO
!****
      Call CONSERV_WM (A)
      Do 110 J=J1A,JNA
      Do 110 I=I1A,INA
  110 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,ATMO)
!**** Include AtmLnd%EVAPOR after EARTH1
      If (QEARTH1)  Then
         Do 121 J=J1A,JNA
         Do 120 I=I1A,aIMAXJ(J)
  120    A(I,J) = AtmLnd%EVAPOR(I,J) * FEARTH(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  121    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  ATMO = ATMO + AGLOBE  ;  EndIf

!****
!**** Calculate ground water GRNW
!****
      Call CONSERV_WTG (A)
      Do 210 J=J1A,JNA
      Do 210 I=I1A,INA
  210 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,GRNW)
!**** Include PREC after CONDSE,PRECIO,PRECOC,PRECIK,PRECGI,IRRIG,PRECLK
      If (QCONDSE .or. QPRECIO .or. QPRECOC .or. QPRECIK .or. QPRECGI .or. QIRRIG .or. QPRECLK)  Then
         Do 221 J=J1A,JNA
         Do 220 I=I1A,aIMAXJ(J)
  220    A(I,J) = PREC(I,J) * FEARTH(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  221    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  GRNW = GRNW + AGLOBE  ;  EndIf
!**** Include IRRIG_WATER_ACT (m/s) after IRRIG,PRECLK
      If (QIRRIG .or. QPRECLK)  Then
         Do 231 J=J1A,JNA
         Do 230 I=I1A,aIMAXJ(J)
  230    A(I,J) = IRRIG_WATER_ACT(I,J) * DTSRC * RHOW * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  231    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  GRNW = GRNW + AGLOBE  ;  EndIf

!****
!**** Calculate lake ice mass LAKI
!****
      Call CONSERV_LMSI (A)
      Do 310 J=J1A,JNA
      Do 310 I=I1A,INA
  310 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,LAKI)
!**** Include PREC after CONDSE,PRECIO,PRECOC
      If (QCONDSE .or. QPRECIO .or. QPRECOC)  Then
         Do 321 J=J1A,JNA
         Do 320 I=I1A,aIMAXJ(J)
  320    A(I,J) = PREC(I,J) * SI_Atm%RSI(I,J) * FLAKE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  321    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LAKI = LAKI + AGLOBE  ;  EndIf
!**** Subtract AtmIce%EVAPOR after ATMDIF,GRNDGI,SDIAG3
      If (QATMDIF .or. QGRNDGI .or. QSDIAG3)  Then
         Do 351 J=J1A,JNA
         Do 350 I=I1A,aIMAXJ(J)
  350    A(I,J) = AtmIce%EVAPOR(I,J) * SI_Atm%RSI(I,J) * FLAKE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  351    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LAKI = LAKI - AGLOBE  ;  EndIf
!**** Include IceOcn%DMSI after GRNDLK,RIVERF
      If (QGRNDLK .or. QRIVERF)  Then
         Do 361 J=J1A,JNA
         Do 360 I=I1A,aIMAXJ(J)
  360    A(I,J) = (IceLak%DMSI(1,I,J) * (1 - SI_Atm%RSI(I,J)) + &
                   IceLak%DMSI(2,I,J) *      SI_Atm%RSI(I,J)) * FLAKE(I,J) * aXYP(I,J)
         If (J==1 .or. J==JMA)  A(2:IMA,J) = A(1,J)
  361    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LAKI = LAKI + AGLOBE  ;  EndIf

!****
!**** Calculate liquid lake mass LIQL
!****
      Call CONSERV_LKM (A)
      Do 410 J=J1A,JNA
      Do 410 I=I1A,INA
  410 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,LIQL)
!**** Include IceLak%MELTI after MELTIK,CONDSE,PRECIO,PRECOC,PRECIK,PRECGI,IRRIG
      If (QMELTIK .or. QCONDSE .or. QPRECIO .or. QPRECOC .or. QPRECIK .or. QPRECGI .or. QIRRIG)  Then
         Do 421 J=J1A,JNA
         Do 420 I=I1A,aIMAXJ(J)
         A(I,J) = 0
  420    If (IceLak%FWATER(I,J) > 0)  A(I,J) = IceLak%MELTI(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  421    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf
!**** Include PREC after CONDSE,PRECIO,PRECOC,PRECIK,PRECGI,IRRIG
      If (QCONDSE .or. QPRECIO .or. QPRECOC .or. QPRECIK .or. QPRECGI .or. QIRRIG)  Then
         Do 431 J=J1A,JNA
         Do 430 I=I1A,aIMAXJ(J)
  430    A(I,J) = PREC(I,J) * (1 - SI_Atm%RSI(I,J)) * FLAKE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  431    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf
!**** Include IceLak%RUNPSI after PRECIK,PRECGI,IRRIG
      If (QPRECIK .or. QPRECGI .or. QIRRIG)  Then
         Do 441 J=J1A,JNA
         Do 440 I=I1A,aIMAXJ(J)
  440    A(I,J) = IceLak%RUNPSI(I,J) * SI_Atm%RSI(I,J) * FLAKE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  441    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf
!**** Include AtmGla%RUNO after PRECGI,IRRIG
      If (QPRECGI .or. QIRRIG)  Then
         Do 451 J=J1A,JNA
         Do 450 I=I1A,aIMAXJ(J)
  450    A(I,J) = AtmGla%RUNO(I,J) * FLICE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  451    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf
!**** Include AtmLnd%RUNO after EARTH1,ATMDIF,GRNDGI,SDIAG3,GRNDIK
      If (QEARTH1 .or. QATMDIF .or. QGRNDGI .or. QSDIAG3 .or. QGRNDIK)  Then
         Do 461 J=J1A,JNA
         Do 460 I=I1A,aIMAXJ(J)
  460    A(I,J) = AtmLnd%RUNO(I,J) * FEARTH(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  461    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf
!**** Subtract AtmOcn%EVAPOR after ATMDIF,GRNDGI,SDIAG3,GRNDIK
      If (QATMDIF .or. QGRNDGI .or. QSDIAG3 .or. QGRNDIK)  Then
         Do 471 J=J1A,JNA
         Do 470 I=I1A,aIMAXJ(J)
  470    A(I,J) = AtmOcn%EVAPOR(I,J) * (1 - SI_Atm%RSI(I,J)) * FLAKE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  471    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL - AGLOBE  ;  EndIf
!**** Include AtmGla%RUNO after GRNDGI,SDIAG3,GRNDIK
      If (QGRNDGI .or. QSDIAG3 .or. QGRNDIK)  Then
         Do 481 J=J1A,JNA
         Do 480 I=I1A,aIMAXJ(J)
  480    A(I,J) = AtmGla%RUNO(I,J) * FLICE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  481    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf
!**** Include IceLak%RUNOSI after GRNDIK
      If (QGRNDIK)  Then
         Do 491 J=J1A,JNA
         Do 490 I=I1A,aIMAXJ(J)
  490    A(I,J) = IceLak%RUNOSI(I,J) * SI_Atm%RSI(I,J) * FLAKE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  491    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQL = LIQL + AGLOBE  ;  EndIf

!****
!**** Calculate land ice LANI
!****
      Call CONSERV_MLI (A)
      Do 510 J=J1A,JNA
      Do 510 I=I1A,INA
  510 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,LANI)
!**** Include PREC after CONDSE,PRECIO,PRECOC,PRECIK
      If (QCONDSE .or. QPRECIO .or. QPRECOC .or. QPRECIK)  Then
         Do 521 J=J1A,JNA
         Do 520 I=I1A,aIMAXJ(J)
  520    A(I,J) = PREC(I,J) * FLICE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  521    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LANI = LANI + AGLOBE  ;  EndIf
!**** Subtract AtmGlas(1)%EVAPOR after ATMDIF
      If (QATMDIF)  Then
         Do 531 J=J1A,JNA
         Do 530 I=I1A,aIMAXJ(J)
  530    A(I,J) = AtmGlas(1)%EVAPOR(I,J) * FLICE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  531    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LANI = LANI - AGLOBE  ;  EndIf

!****
!**** Calculate ice berg mass ICEB
!****
      Call CONSERV_MICB (A)
      Do 610 J=J1A,JNA
      Do 610 I=I1A,INA
  610 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,ICEB)
!**** Include IMPLM after PRECGI,IRRIG,PRECLK,EARTH1,ATMDIF,GRNDGI
      If (QPRECGI .or. QIRRIG .or. QPRECLK .or. QEARTH1 .or. QATMDIF .or. QGRNDGI)  Then
         Do 621 J=J1A,JNA
         Do 620 I=I1A,aIMAXJ(J)
  620    A(I,J) = AtmGla%IMPLM(I,J) * FLICE(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  621    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LANI = LANI + AGLOBE  ;  EndIf
!**** Subtract (ACCPDA+ACCPDG)/EDperY after DAILYO,DAILYI
!*    If (QROOT .and. (SUBR=='DAILYO' .or. SUBR=='DAILYI'))  ICEB = ICEB - (ACCPDA+ACCPDG)/EDperY

!****
!**** Calculate sea ice mass SEAI
!****
      Call CONSERV_OMSI (A)
      Do 710 J=J1A,JNA
      Do 710 I=I1A,INA
  710 A(I,J) = A(I,J) * aXYP(I,J)
      Call GLOBALSUM (aGrid,A,SEAI)
!**** Include PREC after CONDSE
      If (QCONDSE)  Then
         Do 721 J=J1A,JNA
         Do 720 I=I1A,aIMAXJ(J)
  720    A(I,J) = PREC(I,J) * SI_Atm%RSI(I,J) * aFOCEAN(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  721    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  SEAI = SEAI + AGLOBE  ;  EndIf
!**** Subtract AtmIce%EVAPOR after ATMDIF,GRNDGI,SDIAG3,GRNDIK,GRNDLK,RIVERF,FORMIK
      If (QATMDIF .or. QGRNDGI .or. QSDIAG3 .or. QGRNDIK .or. QGRNDLK .or. QRIVERF .or. QFORMIK)  Then
         Do 731 J=J1A,JNA
         Do 730 I=I1A,aIMAXJ(J)
  730    A(I,J) = AtmIce%EVAPOR(I,J) * SI_Atm%RSI(I,J) * aFOCEAN(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  731    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  SEAI = SEAI - AGLOBE  ;  EndIf
!**** Include oDMSI after GRNDOC
      If (QGRNDOC)  Then
         Do 741 J=J1O,JNO
         Do 740 I=I1O,oIMAXJ(J)
  740    O(I,J) = (oDMSI(1,I,J) * (1 - oRSI(I,J)) + oDMSI(2,I,J) * oRSI(I,J)) * oFOCEAN(I,J) * oXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMO))  O(2:IMO,J) = O(1,J)
  741    Continue
         Call GLOBALSUM (oGrid,O,AGLOBE)
         If (QROOT)  SEAI = SEAI + AGLOBE  ;  EndIf

!****
!**** Calculate liquid ocean mass LIQO
!****
      Call CONSERV_OMS (O)
      Do 810 J=J1O,JNO
      Do 810 I=1,IMO
  810 O(I,J) = O(I,J) * oXYP(I,J)
      Call GLOBALSUM (oGrid,O,LIQO)
!**** Include IceOcn%MELTI after MELTIO,MELTIK,CONDSE,PRECIO,PRECOC,PRECIK,PRECGI,IRRIG ,PRECLK,
!****                            EARTH1,ATMDIF,GRNDGI,SDIAG3,GRNDIK,GRNDLK,RIVERF,FORMIK,GRNDIO
      If (QMELTIO .or. QMELTIK .or. QCONDSE .or. QPRECIO .or. QPRECOC .or. QPRECIK .or. QPRECGI .or. QIRRIG  .or. QPRECLK .or. &
          QEARTH1 .or. QATMDIF .or. QGRNDGI .or. QSDIAG3 .or. QGRNDIK .or. QGRNDLK .or. QRIVERF .or. QFORMIK .or. QGRNDIO)  Then
         Do 821 J=J1A,JNA
         Do 820 I=I1A,aIMAXJ(J)
         A(I,J) = 0
  820    If (IceOcn%FWATER(I,J) > 0)  A(I,J) = IceOcn%MELTI(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  821    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQO = LIQO + AGLOBE  ;  EndIf
!**** Include PREC after CONDSE,PRECIO
      If (QCONDSE .or. QPRECIO)  Then
         Do 831 J=J1A,JNA
         Do 830 I=I1A,aIMAXJ(J)
  830    A(I,J) = PREC(I,J) * (1 - SI_Atm%RSI(I,J)) * aFOCEAN(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  831    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQO = LIQO + AGLOBE  ;  EndIf
!**** Include IceOcn%RUNPSI after PRECIO
      If (QPRECIO)  Then
         Do 841 J=J1A,JNA
         Do 840 I=I1A,aIMAXJ(J)
  840    A(I,J) = IceOcn%RUNPSI(I,J) * SI_Atm%RSI(I,J) * aFOCEAN(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  841    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQO = LIQO + AGLOBE  ;  EndIf
!**** Subtract AtmOcn%EVAPOR after ATMDIF,GRNDGI,SDIAG3,GRNDIK,GRNDLK,RIVERF,FORMIK,GRNDIO
      If (QATMDIF .or. QGRNDGI .or. QSDIAG3 .or. QGRNDIK .or. QGRNDLK .or. QRIVERF .or. QFORMIK .or. QGRNDIO)  Then
         Do 851 J=J1A,JNA
         Do 850 I=I1A,aIMAXJ(J)
  850    A(I,J) = AtmOcn%EVAPOR(I,J) * (1 - SI_Atm%RSI(I,J)) * aFOCEAN(I,J) *aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  851    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQO = LIQO - AGLOBE  ;  EndIf
!**** Include IceOcn%RUNOSI after GRNDIO
      If (QGRNDIO)  Then
         Do 861 J=J1A,JNA
         Do 860 I=I1A,aIMAXJ(J)
  860    A(I,J) = IceOcn%RUNOSI(I,J) * SI_Atm%RSI(I,J) * aFOCEAN(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  861    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQO = LIQO + AGLOBE  ;  EndIf
!**** Include AtmOcn%FLOWO after RIVERF,FORMIK,GRNDIO
      If (QRIVERF .or. QFORMIK .or. QGRNDIO)  Then
         Do 871 J=J1A,JNA
         Do 870 I=I1A,aIMAXJ(J)
  870    A(I,J) = AtmOcn%FLOWO(I,J) * aFOCEAN(I,J) * aXYP(I,J)
         If (QLL .and. (J==1 .or. J==JMA))  A(2:IMA,J) = A(1,J)
  871    Continue
         Call GLOBALSUM (aGrid,A,AGLOBE)
         If (QROOT)  LIQO = LIQO + AGLOBE  ;  EndIf

!**** Write global data to unit 6
      If (.not.QROOT)  GoTo 900
      ATMO = ATMO/AREAG
      GRNW = GRNW/AREAG
      LAKI = LAKI/AREAG
      LIQL = LIQL/AREAG
      LANI = LANI/AREAG
      ICEB = ICEB/AREAG
      SEAI = SEAI/AREAG
      LIQO = LIQO/AREAG
      TOTL = ATMO + GRNW + LAKI + LIQL + LANI + ICEB + SEAI + LIQO

      If (SUBR .eq. 'MAIN  ')  Write (6,988)
      Write (6,989)
      Write (6,989) JMON,JDATE,JHOUR, SUBR, ATMO,GRNW,LAKI,LIQL,LANI,ICEB,SEAI,LIQO,TOTL

  900 Return
!****
  988 Format ('WATER:' // &
              'WATER:    Date   (kg/m^2)   AtmosQ      Ground   LakeIce     LiqLake    LandIce    IceBerg', &
              '     SeaIce     LiquidOcean       Summation')
  989 Format ('WATER:',I4.2,2('/',I2.2), A8, F11.6,F12.6,F10.6,F12.6,3F11.6,2F16.6)
      EndSubroutine WATER


!!!   Calls to subroutine WATER are added to following source codes (based on master branch of 2018/03/07):

!!!   D:/discover/nobackup/grussell/LLC28/model>> diff MODELE.f MODELE_WATER.f  
!!!   293a294,295
!!!   >            Call WATER ('MAIN  ')      
!!!   303a306
!!!   >            Call WATER ('DAILY ')
!!!   321a325
!!!   >            Call WATER ('PRECIO')
!!!   322a327
!!!   >            Call WATER ('PRECOC')

!!!   D:/discover/nobackup/grussell/LLC28/model>> diff SURFACE.f SURFACE_WATER.f
!!!   260c260
!!!   >            Call WATER ('PRECIK')
!!!   263a264
!!!   >            Call WATER ('PRECGI')
!!!   269a271
!!!   >            Call WATER ('PRECGI')
!!!   289a292
!!!   >            Call WATER ('IRRIG ')
!!!   291a295
!!!   >            Call WATER ('PRECLK')
!!!   856a861
!!!   >            If (NS==1)  Call WATER ('EARTH1')
!!!   1049a1055
!!!   >            Call WATER ('ATMDIF')
!!!   1072a1079
!!!   >            Call WATER ('GRNDGI')
!!!   1081a1089
!!!   >            Call WATER ('GRNDGI')
!!!   1101a1110
!!!   >            Call WATER ('SDIAG3')
!!!   1107a1117
!!!   >            Call WATER ('GRNDIK')
!!!   1109a1120
!!!   >            Call WATER ('GRNDLK')
!!!   1112a1124
!!!   >            Call WATER ('RIVERF')
!!!   1113a1126
!!!   >            Call WATER ('FORMIK')

!!!   D:/discover/nobackup/grussell/LLC28/model>> diff ATM_DRV.f ATM_DRV_WATER.f
!!!   257a258
!!!   >            Call WATER ('MELTIO')
!!!   258a260
!!!   >            Call WATER ('MELTIK')
!!!   264a267
!!!   >            Call WATER ('CONDSE')

!!!   D:/discover/nobackup/grussell/LLC28/model>> diff OCN_DRV.f OCN_DRV_WATER.f
!!!   32a33
!!!   >            Call WATER ('GRNDIO')
!!!   46a48
!!!   >            Call WATER ('FORMIO')

!!!   D:/discover/nobackup/grussell/LLC28/model>> diff OCNDYN2.f OCNDYN2_WATER.f
!!!   86a87
!!!   >            Call WATER ('GRNDOC')

