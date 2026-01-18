#include "rundeck_opts.h"

      Subroutine CARBON (SUBR)
!@sum  CARBON writes global carbon reservoirs to unit 6
!@auth Original Development Team
!@ver  2010/11/24

C**** Output: Line written to unit 6 of carbon reservoirs (10^-6 kg/m^2)
C****   Summed carbon reservoirs (kg) are:
C****   DIAT = TRMO(5)  + TRMST(5)
C****   CHLO = TRMO(6)  + TRMST(6)
C****   CYAN = TRMO(7)  + TRMST(7)
C****   COCC = TRMO(8)  + TRMST(8)
C****   HERB = TRMO(9)  + TRMST(9)
C****   NDET = TRMO(10) + TRMST(10)
C****    DOC = TRMO(13) + TRMST(13)
C****    DIC = TRMO(14) + TRMST(14)
C****   TOTL = Sum (TRMO + TRMST)
C****
C**** Input: SUBR = 6 character string which labels ouput line

      Use MODEL_COM, Only: DTSRC, modeleclock
      USE MODEL_COM, only: nstep=>itime
      Use GEOM,      Only: AREAG, aXYP
      Use OCEAN,     Only: LMO, oXYP, TRMO,
     *                     oFOCEAN=>FOCEAN
      Use STRAITS,   Only: TRMST
      Use FLUXES,    Only: atmocn
      Use OFLUXES,   Only: ocnatm
      Use DOMAIN_DECOMP_ATM, Only: aGRID=>GRID, AM_I_ROOT,GLOBALSUM
      Use OCEANR_DIM,       Only: oGRID

      Implicit  None
!@var SUBR identifies where CHECK was called from
      Character*6,Intent(In) :: SUBR

C**** Local Variables
      Integer*4,Save :: IFIRST=1
      Integer*4 L,LMAX, J1A,JNA, J1O,JNO
      Real*8,Save :: LAST = 0
      Real*8 :: DIAT(LMO),CHLO(LMO),CYAN(LMO),COCC(LMO),
     *     HERB(LMO),NDET(LMO), DOC(LMO), DIC(LMO),TOTL(LMO),
     *     A(agrid%im_world,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO), aFLUX,
     *     O(ogrid%im_world,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), oFLUX
      integer :: idx_co2,ii,jj

#ifdef OBIO_QUIET_MODE
      return
#endif

      idx_co2=atmocn%gasex_index%getindex(atmocn%n_co2n)
      if (idx_co2<1) return
C**** Extract domain decomposition band parameters
      J1A = aGRID%J_STRT
      JNA = aGRID%J_STOP
      J1O = oGRID%J_STRT
      JNO = oGRID%J_STOP

C**** Calculate DIAT
      Do 10 L=1,LMO
      O(:,:) = TRMO(:,:,L,5)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,DIAT(L))
      If (AM_I_ROOT())  DIAT(L) = DIAT(L) + Sum(TRMST(L,:,5))
   10 Continue

C**** Calculate CHLO
      Do 20 L=1,LMO
      O(:,:) = TRMO(:,:,L,6)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,CHLO(L))
      If (AM_I_ROOT())  CHLO(L) = CHLO(L) + Sum(TRMST(L,:,6))
   20 Continue

C**** Calculate CYAN
      Do 30 L=1,LMO
      O(:,:) = TRMO(:,:,L,7)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,CYAN(L))
      If (AM_I_ROOT())  CYAN(L) = CYAN(L) + Sum(TRMST(L,:,7))
   30 Continue

C**** Calculate COCC
      Do 40 L=1,LMO
      O(:,:) = TRMO(:,:,L,8)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,COCC(L))
      If (AM_I_ROOT())  COCC(L) = COCC(L) + Sum(TRMST(L,:,8))
   40 Continue

C**** Calculate HERB
      Do 50 L=1,LMO
      O(:,:) = TRMO(:,:,L,9)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,HERB(L))
      If (AM_I_ROOT())  HERB(L) = HERB(L) + Sum(TRMST(L,:,9))
   50 Continue

C**** Calculate NDET
      Do 60 L=1,LMO
      O(:,:) = TRMO(:,:,L,10)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,NDET(L))
      If (AM_I_ROOT())  NDET(L) = NDET(L) + Sum(TRMST(L,:,10))
   60 Continue

C**** Calculate DOC
      Do 70 L=1,LMO
      O(:,:) = TRMO(:,:,L,13)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,DOC(L))
      If (AM_I_ROOT())  DOC(L) = DOC(L) + Sum(TRMST(L,:,13))
   70 Continue

C**** Calculate DIC
      Do 80 L=1,LMO
      O(:,:) = TRMO(:,:,L,14)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,DIC(L))
      If (AM_I_ROOT())  DIC(L) = DIC(L) + Sum(TRMST(L,:,14))
   80 Continue
      If (SUBR == 'SURFCE')  Then
         A(:,:)=atmocn%trgasex(idx_co2,:,:)*atmocn%focean(:,:)*aXYP(:,:)
         If (J1A==1)    A(2:agrid%im_world,1)   = A(1,1)
         If (JNA==agrid%jm_world)
     &      A(2:agrid%im_world,agrid%jm_world) = A(1,agrid%jm_world)
         Call GLOBALSUM (aGRID,A,aFLUX)
         If (AM_I_ROOT())
     *      DIC(1) = DIC(1) + aFLUX * DTSRC * 12.*1d-3
         EndIf
      If (SUBR == 'AG2OG_' .or. SUBR == 'GRNDOC' .or. SUBR == 'OCONV ')
     .  Then
         O(:,J1O:JNO) = ocnatm%trgasex(idx_co2,:,J1O:JNO) 
     .                * oFOCEAN(:,J1O:JNO) * oXYP(:,J1O:JNO)
         If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
         If (JNO==ogrid%jm_world)
     &      O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
         Call GLOBALSUM (oGRID,O,oFLUX)
         If (AM_I_ROOT())
     *      DIC(1) = DIC(1) + oFLUX * DTSRC * 12.*1d-3
         EndIf

C**** Write global data to unit 6
      If (.not.AM_I_ROOT())  GoTo 900
      DIAT(:) = DIAT(:) * 1d6 / AREAG
      CHLO(:) = CHLO(:) * 1d6 / AREAG
      CYAN(:) = CYAN(:) * 1d6 / AREAG
      COCC(:) = COCC(:) * 1d6 / AREAG
      HERB(:) = HERB(:) * 1d6 / AREAG
      NDET(:) = NDET(:) * 1d6 / AREAG
       DOC(:) =  DOC(:) * 1d6 / AREAG
       DIC(:) =  DIC(:) * 1d6 / AREAG
      TOTL(:) = DIAT(:) + CHLO(:) + CYAN(:) + COCC(:) +
     +          HERB(:) + NDET(:) +  DOC(:) +  DIC(:)

      If (IFIRST == 1)  Then
         Write (6,909)
         IFIRST = 0  ;  Endif
C     Write (6,910)
C     LMAX = LMO  ;  If(SUBR=='SURFCE' .or. SUBR=='AG2OG_')  LMAX = 1
C     Do 100 L=1,LMAX
C 100 Write (6,910) SUBR,modeleclock%getyear(),modeleclock%getmonth(),
C    *              modeleclock%getdate(),modeleclock%gethour(),
C    *              DIAT(L),CHLO(L),CYAN(L),COCC(L),
C    *              HERB(L),NDET(L), DOC(L), DIC(L),TOTL(L)
      Write (6,910) SUBR,modeleclock%getyear(),modeleclock%getmonth(),
     *              modeleclock%getdate(),modeleclock%gethour(),
     *              Sum(DIAT(:)),Sum(CHLO(:)),Sum(CYAN(:)),
     *              Sum(COCC(:)),Sum(HERB(:)),Sum(NDET(:)),
     *              Sum( DOC(:)),Sum( DIC(:)),Sum(TOTL(:)),
     *              Sum(TOTL(:))-LAST
      LAST = Sum (TOTL(:))

  900 Return
C****
  909 Format ('CARBON  (10^-6 kg/m^2)           DIAT     CHLO ',
     *        '    CYAN     COCC     HERB       NDET         DOC  ',
     *        '         DIC            Total    Dif ')
  910 Format ('CARBON:',A7,I6,3('/',I2.2),5F9.3,F11.3,F13.3,2F14.3
     .         ,F11.6)
      End Subroutine CARBON

      Subroutine NITR (SUBR)
!@sum  NITR   writes global nitrate reservoirs to unit 6
!@auth Original Development Team
!@ver  2010/11/24

C**** Output: Line written to unit 6 of nitrate reservoirs (10^-6 kg/m^2)
C****   Summed nitrate reservoirs (kg) are:
C****   DIAT = TRMO(5)  + TRMST(5)
C****   CHLO = TRMO(6)  + TRMST(6)
C****   CYAN = TRMO(7)  + TRMST(7)
C****   COCC = TRMO(8)  + TRMST(8)
C****   HERB = TRMO(9)  + TRMST(9)
C****   NDET = TRMO(10) + TRMST(10)
C****    DOC = TRMO(13) + TRMST(13)
C****    DIC = TRMO(14) + TRMST(14)
C****    NIT = TRMO(1)  + TRMST(1)
C****   AMMO = TRMO(2)  + TRMST(2)
C****   TOTL = Sum (TRMO + TRMST)
C****
C**** Input: SUBR = 6 character string which labels ouput line

      Use MODEL_COM, Only: modeleclock
      Use GEOM,      Only: AREAG
      Use OCEAN,     Only: LMO, TRMO
      Use STRAITS,   Only: TRMST
      Use DOMAIN_DECOMP_ATM, Only: aGRID=>GRID, AM_I_ROOT,GLOBALSUM
      Use OCEANR_DIM,       Only: oGRID

      Implicit  None
!@var SUBR identifies where CHECK was called from
      Character*6,Intent(In) :: SUBR

C**** Local Variables
      Integer*4,Save :: IFIRST=1
      Integer*4 L,LMAX, J1A,JNA, J1O,JNO
      Real*8,Save :: LAST = 0
      Real*8 :: DIAT(LMO),CHLO(LMO),CYAN(LMO),COCC(LMO),DIC(LMO),
     *     HERB(LMO),NDET(LMO), DOC(LMO), NIT(LMO),TOTL(LMO),
     *     A(agrid%im_world,aGRID%J_STRT_HALO:aGRID%J_STOP_HALO), aFLUX,
     *     O(ogrid%im_world,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), oFLUX,
     *     AMMO(LMO)

#ifdef OBIO_QUIET_MODE
      return
#endif

C**** Extract domain decomposition band parameters
      J1A = aGRID%J_STRT
      JNA = aGRID%J_STOP
      J1O = oGRID%J_STRT
      JNO = oGRID%J_STOP

C**** Calculate DIAT
      Do 10 L=1,LMO
      O(:,:) = TRMO(:,:,L,5)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &       O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,DIAT(L))
      If (AM_I_ROOT())  DIAT(L) = DIAT(L) + Sum(TRMST(L,:,5))
   10 Continue

C**** Calculate CHLO
      Do 20 L=1,LMO
      O(:,:) = TRMO(:,:,L,6)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,CHLO(L))
      If (AM_I_ROOT())  CHLO(L) = CHLO(L) + Sum(TRMST(L,:,6))
   20 Continue

C**** Calculate CYAN
      Do 30 L=1,LMO
      O(:,:) = TRMO(:,:,L,7)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,CYAN(L))
      If (AM_I_ROOT())  CYAN(L) = CYAN(L) + Sum(TRMST(L,:,7))
   30 Continue

C**** Calculate COCC
      Do 40 L=1,LMO
      O(:,:) = TRMO(:,:,L,8)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &       O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,COCC(L))
      If (AM_I_ROOT())  COCC(L) = COCC(L) + Sum(TRMST(L,:,8))
   40 Continue

C**** Calculate HERB
      Do 50 L=1,LMO
      O(:,:) = TRMO(:,:,L,9)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,HERB(L))
      If (AM_I_ROOT())  HERB(L) = HERB(L) + Sum(TRMST(L,:,9))
   50 Continue

C**** Calculate NDET
      Do 60 L=1,LMO
      O(:,:) = TRMO(:,:,L,10)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &     O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,NDET(L))
      If (AM_I_ROOT())  NDET(L) = NDET(L) + Sum(TRMST(L,:,10))
   60 Continue

C**** Calculate DOC
      Do 70 L=1,LMO
      O(:,:) = TRMO(:,:,L,13)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,DOC(L))
      If (AM_I_ROOT())  DOC(L) = DOC(L) + Sum(TRMST(L,:,13))
   70 Continue

C**** Calculate DIC
      Do 75 L=1,LMO
      O(:,:) = TRMO(:,:,L,14)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,DIC(L))
      If (AM_I_ROOT())  DIC(L) = DIC(L) + Sum(TRMST(L,:,14))
   75 Continue

C**** Calculate NIT
      Do 80 L=1,LMO
      O(:,:) = TRMO(:,:,L,1)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,NIT(L))
      If (AM_I_ROOT()) NIT(L) = NIT(L) + Sum(TRMST(L,:,1))
   80 Continue

C**** Calculate AMMO
      Do 85 L=1,LMO
      O(:,:) = TRMO(:,:,L,2)
      If (J1O==1)    O(2:ogrid%im_world,1)   = O(1,1)
      If (JNO==ogrid%jm_world)
     &    O(2:ogrid%im_world,ogrid%jm_world) = O(1,ogrid%jm_world)
      Call GLOBALSUM (oGRID,O,AMMO(L))
      If (AM_I_ROOT()) AMMO(L) = AMMO(L) + Sum(TRMST(L,:,2))
   85 Continue

C**** Write global data to unit 6
      If (.not.AM_I_ROOT())  GoTo 900
      DIAT(:) = DIAT(:)/(106/16) * 1d6 / AREAG
      CHLO(:) = CHLO(:)/(106/16) * 1d6 / AREAG
      CYAN(:) = CYAN(:)/(106/16) * 1d6 / AREAG
      COCC(:) = COCC(:)/(106/16) * 1d6 / AREAG
      HERB(:) = HERB(:)/(106/16) * 1d6 / AREAG
      NDET(:) = NDET(:)/(106/16) * 1d6 / AREAG
       DOC(:) =  DOC(:)/(106/16) * 1d6 / AREAG
       NIT(:) = NIT(:) * 1d6 / AREAG
      AMMO(:) =AMMO(:) * 1d6 / AREAG
      TOTL(:) = DIAT(:) + CHLO(:) + CYAN(:) + COCC(:) +
     +          HERB(:) + NDET(:) +  DOC(:) + NIT(:) + AMMO(:)

      If (IFIRST == 1)  Then
         Write (6,909)
         IFIRST = 0  ;  Endif
C     Write (6,910)
C     LMAX = LMO  ;  If(SUBR=='SURFCE' .or. SUBR=='AG2OG_')  LMAX = 1
C     Do 100 L=1,LMAX
C 100 Write (6,910) SUBR,modeleclock%getyear(),modeleclock%getmonth(),
C    *              modeleclock%getdate(),modeleclock%gethour(),
C    *              DIAT(L),CHLO(L),CYAN(L),COCC(L),
C    *              HERB(L),NDET(L), DOC(L), NIT(L),TOTL(L)
      Write (6,910) SUBR,modeleclock%getyear(),modeleclock%getmonth(),
     *              modeleclock%getdate(),modeleclock%gethour(),
     *              Sum(DIAT(:)),Sum(CHLO(:)),Sum(CYAN(:)),
     *              Sum(COCC(:)),Sum(HERB(:)),Sum(NDET(:)),
     *              Sum( DOC(:)),Sum(AMMO(:)),Sum(NIT(:)),Sum(TOTL(:)),
     *              Sum(TOTL(:))-LAST
      LAST = Sum (TOTL(:))

  900 Return
C****
  909 Format ('NITRATE  (10^-6 kg/m^2)           DIAT     CHLO ',
     *        '    CYAN     COCC     HERB       NDET         DOC  ',
     *        '       AMMO        NITR            Total    Dif ')
  910 Format ('NITRATE:',A7,I6,3('/',I2.2),5F9.3,F11.3,F13.3,3F14.3
     .         ,F11.6)
      End Subroutine NITR
