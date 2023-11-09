      Subroutine OWALL
!@sum  OWALL exerts a drag on Ocean Model's velocities based on vertical walls that the flow encounters
!@auth Gary L. Russell
!@ver  2015/03/13

      Use OCEAN,      Only: IM,JM,LMO,IVSP,IVNP, LMOM=>LMM,LMOU=>LMU,LMOV=>LMV, COSU,SINU, DTS, OWALLU,OWALLV, &
                            UO,VO,UOD,VOD, GXMO,GYMO,GXXMO,GYYMO,GXYMO,GYZMO,GZXMO, SXMO,SYMO,SXXMO,SYYMO,SXYMO,SYZMO,SZXMO, &
                            OCOASTAL_DRAG,VCOAST,USE_QUS
      Use OCEANR_DIM, Only: OGRID
      Use DOMAIN_DECOMP_1D, Only: AM_I_ROOT
#ifdef    TRACERS_OCEAN
          Use OCEAN,  Only: NTM, TXMO,TYMO,TXXMO,TYYMO,TXYMO,TYZMO,TZXMO
#endif
      Implicit None

!**** Local variables
      Integer :: I,J,L,N, J1,JN,J1P,JNP
      Real*8  :: WMAG, UOIM,UOIV

!**** Extract domain decomposition
      J1 = OGRID%J_STRT  ;  JN = OGRID%J_STOP  ;  J1P = Max(J1,2)  ;  JNP = Min(JN,JM-1)

!**** Apply wall drag to UO and VOD components away from poles
!**** UO = UO - UO*DTS*WALL*|WO| =
!****    = UO * (1 - DTS*WALL*|WO|) ~
!****    ~ UO / (1 + DTS*WALL*|WO|)
      Do 140 J=J1P,JNP
      Do 120 I=1,IM-1
      Do 110 L=LMOU(I,J),1,-1
      If (OWALLU(I,J,L) + OWALLU(I+1,J,L) <= 0)  GoTo 120
      WMAG  = Sqrt (UO(I,J,L)**2 + VOD(I,J,L)**2)
       UO(I,J,L) =  UO(I,J,L) / (1 + DTS * .5*(OWALLU(I,J,L) + OWALLU(I+1,J,L)) * WMAG)
  110 VOD(I,J,L) = VOD(I,J,L) / (1 + DTS * .5*(OWALLV(I,J,L) + OWALLV(I+1,J,L)) * WMAG)
  120 Continue
      Do 130 L=LMOU(IM,J),1,-1
      If (OWALLU(IM,J,L) + OWALLU(1,J,L) <= 0)  GoTo 140
      WMAG  = Sqrt (UO(IM,J,L)**2 + VOD(IM,J,L)**2)
       UO(IM,J,L) =  UO(IM,J,L) / (1 + DTS * .5*(OWALLU(IM,J,L) + OWALLU(1,J,L)) * WMAG)
  130 VOD(IM,J,L) = VOD(IM,J,L) / (1 + DTS * .5*(OWALLV(IM,J,L) + OWALLV(1,J,L)) * WMAG)
  140 Continue

!**** Apply wall drag to VO and UOD components
      Do 220 J=J1,JNP  ;  Do 220 I=1,IM
      Do 210 L=LMOV(I,J),1,-1
      If (OWALLV(I,J,L) + OWALLV(I,J+1,L) <= 0)  GoTo 220
      WMAG  = Sqrt (VO(I,J,L)**2 + UOD(I,J,L)**2)
       VO(I,J,L) =  VO(I,J,L) / (1 + DTS * .5*(OWALLV(I,J,L) + OWALLV(I,J+1,L)) * WMAG)
  210 UOD(I,J,L) = UOD(I,J,L) / (1 + DTS * .5*(OWALLU(I,J,L) + OWALLU(I,J+1,L)) * WMAG)
  220 Continue

      If (J1 /= 1)  GoTo 400
!**** Apply wall drag to UO(IM) and UO(IVSP) component at south pole
      Do 310 L=LMOM(1,1),1,-1
      If (OWALLU(IM,1,L) <= 0)  GoTo 400
      WMAG = Sqrt (UO(IM,1,L)**2 + UO(IVSP,1,L)**2)
      UOIM = UO(IM  ,1,L) / (1 + DTS*OWALLU(IM  ,1,L)*WMAG)
      UOIV = UO(IVSP,1,L) / (1 + DTS*OWALLV(IVSP,1,L)*WMAG)
  310 UO(:,1,L) = UOIM*COSU(:) - UOIV*SINU(:)

  400 If (JN /= JM)  GoTo 500
!**** Apply wall drag to UO(IM) and UO(IVSP) component at south pole
      Do 410 L=LMOM(1,JM),1,-1
      If (OWALLU(IM,JM,L) <= 0)  GoTo 500
      WMAG = Sqrt (UO(IM,JM,L)**2 + UO(IVNP,JM,L)**2)
      UOIM = UO(IM  ,JM,L) / (1 + DTS*OWALLU(IM  ,JM,L)*WMAG)
      UOIV = UO(IVNP,JM,L) / (1 + DTS*OWALLV(IVNP,JM,L)*WMAG)
  410 UO(:,JM,L) = UOIM*COSU(:) + UOIV*SINU(:)

!**** Exert drag on Ocean Model's strait velocities based on vertical walls that the flow encounters
!**** UST (m/s) = strait velocity = MUST*DIST/MMST
!****
!****  UST =  UST -  UST*DTS*WALL*|UST| =
!**** MUST = MUST - MUST*DTS*WALL*|MUST|*DIST/MMST =
!****      = MUST * (1 - DTS*WALL*|MUST|*DIST/MMST) ~
!****      ~ MUST / (1 + DTS*WALL*|MUST|*DIST/MMST) =
!****      = MUST*MMST / (MMST + DTS*WALL*|MUST|*DIST) =
! 500 If (AM_I_ROOT())  Then
!        Do 510 N=1,NMST  ;  Do 510 L=1,LMST(N)
! 510    MUST(L,N) = MUST(L,N)*MMST(L,N) / (MMST(L,N) + DTS*OWALST(L,N)*Abs(MUST(L,N))*DIST(N))  ;  EndIf

  500 If (OCOASTAL_DRAG /= 2 .or. VCOAST == 0)  Return
!**** Reduce horizontal gradients of heat and salt encountering walls
      Do 610 J=J1P,JNP  ;  Do 610 I=1,IM  ;  Do 610 L=1,LMOM(I,J)
      GXMO(I,J,L) =  GXMO(I,J,L) * (1 - DTS*VCOAST*OWALLU(I,J,L))
      SXMO(I,J,L) =  SXMO(I,J,L) * (1 - DTS*VCOAST*OWALLU(I,J,L))   
      GYMO(I,J,L) =  GYMO(I,J,L) * (1 - DTS*VCOAST*OWALLV(I,J,L))      
  610 SYMO(I,J,L) =  SYMO(I,J,L) * (1 - DTS*VCOAST*OWALLV(I,J,L))

      If (USE_QUS == 1)  Then
         Do 620 J=J1P,JNP  ;  Do 620 I=1,IM  ;  Do 620 L=1,LMOM(I,J)
         GXXMO(I,J,L) = GXXMO(I,J,L) * (1 - DTS*VCOAST*OWALLU(I,J,L))      
         GZXMO(I,J,L) = GZXMO(I,J,L) * (1 - DTS*VCOAST*OWALLU(I,J,L))      
         SXXMO(I,J,L) = SXXMO(I,J,L) * (1 - DTS*VCOAST*OWALLU(I,J,L))      
         SZXMO(I,J,L) = SZXMO(I,J,L) * (1 - DTS*VCOAST*OWALLU(I,J,L))      
         GYYMO(I,J,L) = GYYMO(I,J,L) * (1 - DTS*VCOAST*OWALLV(I,J,L))      
         GYZMO(I,J,L) = GYZMO(I,J,L) * (1 - DTS*VCOAST*OWALLV(I,J,L))      
         SYYMO(I,J,L) = SYYMO(I,J,L) * (1 - DTS*VCOAST*OWALLV(I,J,L))      
         SYZMO(I,J,L) = SYZMO(I,J,L) * (1 - DTS*VCOAST*OWALLV(I,J,L))      
         GXYMO(I,J,L) = GXYMO(I,J,L) * (1 - DTS*VCOAST*(OWALLU(I,J,L)+OWALLV(I,J,L))*.5)      
  620    SXYMO(I,J,L) = SXYMO(I,J,L) * (1 - DTS*VCOAST*(OWALLU(I,J,L)+OWALLV(I,J,L))*.5)  ;  EndIf

!**** Reduce horizontal gradients of tracers encountering walls
#ifdef   TRACERS_OCEAN
      Do 630 N=1,NTM  ;  Do 630 J=J1P,JNP  ;  Do 630 I=1,IM  ;  Do 630 L=1,LMOM(I,J)
      TXMO(I,J,L,N) =  TXMO(I,J,L,N) * (1 - DTS*VCOAST*OWALLU(I,J,L))
  630 TYMO(I,J,L,N) =  TYMO(I,J,L,N) * (1 - DTS*VCOAST*OWALLV(I,J,L))

      If (USE_QUS == 1)  Then
         Do 640 N=1,NTM  ;  Do 640 J=J1P,JNP  ;  Do 640 I=1,IM  ;  Do 640 L=1,LMOM(I,J)
         TXXMO(I,J,L,N) = TXXMO(I,J,L,N) * (1 - DTS*VCOAST*OWALLU(I,J,L))
         TZXMO(I,J,L,N) = TZXMO(I,J,L,N) * (1 - DTS*VCOAST*OWALLU(I,J,L))
         TYYMO(I,J,L,N) = TYYMO(I,J,L,N) * (1 - DTS*VCOAST*OWALLV(I,J,L))
         TYZMO(I,J,L,N) = TYZMO(I,J,L,N) * (1 - DTS*VCOAST*OWALLV(I,J,L))
  640    TXYMO(I,J,L,N) = TXYMO(I,J,L,N) * (1 - DTS*VCOAST*(OWALLU(I,J,L)+OWALLV(I,J,L))*.5)  ;  EndIf
#endif
      Return
      EndSubroutine OWALL
