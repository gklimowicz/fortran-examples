#include "rundeck_opts.h"

!#define LUS_VERT_ADV

c Will add more documentation if this version becomes the modelE default.

      SUBROUTINE OCEANS(atmocn,iceocn,dynsice)
C****
      USE CONSTANT, only : rhows,grav
      USE MODEL_COM, only : msurf,itime,itimei,DTSRC
      USE OCEAN, only : im,jm,lmo,ndyno,nocean,mo,g0m,s0m,
     *    dts,dtofs,dto,dtolf,mdyno,msgso,dxypo,
     *    ogeoz,ogeoz_sv,opbot,ze,lmm,imaxj, UO,VO,VONP,IVNP, ! VOSP,IVSP,
     *    OBottom_drag,OCoastal_drag,OTIDE,uod,vod,lmu,lmv
      USE OCEAN, only : use_qus,
     *     GXMO,GYMO,GZMO, GXXMO,GYYMO,GZZMO, GXYMO,GYZMO,GZXMO,
     *     SXMO,SYMO,SZMO, SXXMO,SYYMO,SZZMO, SXYMO,SYZMO,SZXMO
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv, i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      USE OCEAN_DYN, only : mmi,smu,smv,smw
      USE DOMAIN_DECOMP_1D, only : getDomainBounds, AM_I_ROOT, 
     &     halo_update, south,north, hasSouthPole, hasNorthPole
      USE OCEANR_DIM, only : grid=>ogrid
      USE ODIAG, only : oijl=>oijl_loc,oij=>oij_loc,
     *    ijl_mo,ijl_g0m,ijl_s0m,  ijl_gflx, ijl_sflx, ijl_mfw2,
     *    ijl_mfu,ijl_mfv,ijl_mfw, ijl_ggmfl,ijl_sgmfl,ij_ssh,ij_pb,
     *    IJ_dEPO_Dyn
      USE OFLUXES, only : ocnatm
#ifdef TRACERS_OCEAN
      USE OCN_TRACER_COM, only : tracerlist, ocn_tracer_entry,n_age
     &         ,n_vent,n_gasx,n_wms1,n_wms2,n_wms3,n_ocfc11,n_ocfc12
     &         ,n_sf6,n_abioDIC
      USE OCEAN, only : trmo,
     &     txmo,tymo,tzmo,txxmo,tyymo,tzzmo,txymo,tyzmo,tzxmo
      Use ODIAG, Only: toijl=>toijl_loc,
     *               toijl_conc,toijl_tflx,toijl_gmfl
      use ocean, only : do_tracer_trans,ntrtrans,
     &     asmu,asmv,asmw,motr,mosv0
#endif
      USE EXCHANGE_TYPES, only : atmocn_xchng_vars,iceocn_xchng_vars
      IMPLICIT NONE
      type(atmocn_xchng_vars) :: atmocn
      type(iceocn_xchng_vars) :: iceocn
      type(iceocn_xchng_vars) :: dynsice
c
      Integer*4 I,J,L,N,NS,NST,NO,NEVEN ; real*8 now
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     &     MO1,MO2, UO1,UO2,UOD1,UOD2, VO1,VO2,VOD1,VOD2
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     OPBOT1,OPBOT2, G0INIT,G0FINAL
      real*8 :: relfac,dt_odiff,TIME
      real*8 :: dtdum,byno,mrat_st

      INTEGER it,jt

c**** Extract domain decomposition info
      INTEGER :: J_0, J_1, J_0H,J_1H, J_0S,J_1S
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif

      ocnatm%updated=.true.
      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1,
     &     J_STRT_SKP = J_0S, J_STOP_SKP = J_1S,
     &     J_STRT_HALO = J_0H, J_STOP_HALO = J_1H)

      byno = 1d0/nocean

C***  Get the data from the atmospheric grid to the ocean grid
      call AG2OG_oceans(atmocn,iceocn)
#ifdef TRACERS_GASEXCH_ocean_CO2
         Call CARBON ('AG2OG_')
         Call NITR ('AG2OG_')
#endif

C***  Interpolate DYNSI outputs to the ocean grid
C***  (at present, only the ice-ocean stress is used)
      call IG2OG_oceans(dynsice)

c-------------------------------------------------------------------
c Begin ocean-processors-only code region
      ocean_processors_only: if(grid%have_domain) then
c-------------------------------------------------------------------

      OGEOZ_SV(:,:)=OGEOZ(:,:)

C**** Apply surface fluxes to ocean
      CALL GROUND_OC
         CALL CHECKO('GRNDOC')
#ifdef TRACERS_GASEXCH_ocean_CO2
         Call CARBON ('GRNDOC')
         Call NITR ('GRNDOC')
#endif

#ifdef TRACERS_OCEAN
      ! set or increment the seawater mass array corresponding to tracer state
      if(ntrtrans.eq.1) then
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          motr(i,j,l) = mo(i,j,l)
        enddo
        enddo
        enddo
        enddo
      else !if(ntrtrans.gt.1) then
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          motr(i,j,l) = motr(i,j,l) + (mo(i,j,l)-mosv0(i,j,l))
        enddo
        enddo
        enddo
        enddo
      endif
#endif

C**** Apply ice/ocean and air/ocean stress to ocean
      CALL OSTRES2
         CALL CHECKO('OSTRES')
         CALL TIMER (NOW,MSURF)
         IF (ATMOCN%MODD5S == 0) CALL DIAGCO (11,atmocn)

C**** Apply ocean vertical mixing
      CALL OCONV
         CALL CHECKO('OCONV ')
#ifdef TRACERS_GASEXCH_ocean_CO2
         Call CARBON ('OCONV ')
         Call NITR ('OCONV ')
#endif

C**** Apply bottom and coastal drags
      if (OBottom_drag  == 1) CALL OBDRAG2
      if (OCoastal_drag == 1) CALL OCOAST

C**** Add ocean biology
#ifdef TRACERS_OceanBiology
      call ocnstate_derived
      call obio_model(0)
#ifdef TRACERS_GASEXCH_ocean_CO2
         Call CARBON ('OBIO_M')
         Call NITR ('OBIO_M')
#endif
#ifndef STANDALONE_OCEAN
      IF (ATMOCN%MODD5S == 0) CALL DIAGCO (13,atmocn)
#endif
#endif

         CALL TIMER (NOW,MSGSO)

c relax UOD,VOD toward 4-pt avgs of UO,VO
      call halo_update(grid,uo,from=north)
      call halo_update(grid,vo,from=south)
      relfac = .005d0
      do l=1,lmo
        if(hasNorthPole(grid)) then
          j = jm-1
          call polevel(uo(1,j_0h,l),vo(1,j_0h,l),l)
          i=1
          if(l.le.lmv(i,j)) then
            uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &           uo(im,j,l)+uo(i,j,l)+2.*uo(i,j+1,l)
     &           )
          endif
          do n=1,nbyzv(j,l)
            do i=max(2,i1yzv(n,j,l)),i2yzv(n,j,l)
              uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &             uo(i-1,j,l)+uo(i,j,l)+2.*uo(i,j+1,l)
     &             )
            enddo
          enddo
        endif
        do j=j_0s,min(jm-2,j_1s)
          i=1
          if(l.le.lmv(i,j)) then
            uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &           uo(im,j,l)+uo(i,j,l)+uo(im,j+1,l)+uo(i,j+1,l)
     &           )
          endif
          do n=1,nbyzv(j,l)
            do i=max(2,i1yzv(n,j,l)),i2yzv(n,j,l)
              uod(i,j,l) = (1.-relfac)*uod(i,j,l) + relfac*.25*(
     &             uo(i-1,j,l)+uo(i,j,l)+uo(i-1,j+1,l)+uo(i,j+1,l)
     &             )
            enddo
          enddo
        enddo
        do j=j_0s,j_1s
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
              vod(i,j,l) = (1.-relfac)*vod(i,j,l) + relfac*.25*(
     &             vo(i,j-1,l)+vo(i+1,j-1,l)+vo(i,j,l)+vo(i+1,j,l)
     &             )
            enddo
          enddo
          i=im
          if(l.le.lmu(i,j)) then
            vod(i,j,l) = (1.-relfac)*vod(i,j,l) + relfac*.25*(
     &           vo(i,j-1,l)+vo(1,j-1,l)+vo(i,j,l)+vo(1,j,l)
     &           )
          endif
        enddo
      enddo

      DO L=1,LMO
        MO1(:,:,L) = 0
        UO1(:,:,L) = 0
        VO1(:,:,L) = 0
        UOD1(:,:,L) = 0
        VOD1(:,:,L) = 0
        MO2(:,:,L) = 0
        UO2(:,:,L) = 0
        VO2(:,:,L) = 0
        UOD2(:,:,L) = 0
        VOD2(:,:,L) = 0
      EndDo

C****
C**** Integrate Ocean Dynamics
C****
      Do NO=1,NOCEAN

c check which of these are necessary
      call halo_update(grid,mo)
      call halo_update(grid,uo)
      call halo_update(grid,vo)

      CALL ODHORZ0

      do l=1,lmo
        SMU(:,J_0:J_1,L) = 0
        SMV(:,J_0H:J_1,L) = 0
c        SMW(:,J_0:J_1H,L) = 0 ! not summed
        do j=max(1,j_0h),min(jm,j_1h)
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              mo1(i,j,l) = mo(i,j,l)
              mo2(i,j,l) = mo(i,j,l)
            enddo
          enddo
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              uo1(i,j,l) = uo(i,j,l)
              uo2(i,j,l) = uo(i,j,l)
              vod1(i,j,l) = vod(i,j,l)
              vod2(i,j,l) = vod(i,j,l)
            enddo
          enddo
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              vo1(i,j,l) = vo(i,j,l)
              vo2(i,j,l) = vo(i,j,l)
              uod1(i,j,l) = uod(i,j,l)
              uod2(i,j,l) = uod(i,j,l)
            enddo
          enddo
        enddo
        if(hasNorthPole(grid)) then
          mo1(:,jm,l) = mo(:,jm,l)
          mo2(:,jm,l) = mo(:,jm,l)
          uo1(:,jm,l) = uo(:,jm,l)
          uo2(:,jm,l) = uo(:,jm,l)
          vo1(:,jm,l) = vo(:,jm,l)
          vo2(:,jm,l) = vo(:,jm,l)
        endif
      enddo
      do j=max(1,j_0h),min(jm,j_1h)
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            opbot1(i,j) = opbot(i,j)
            opbot2(i,j) = opbot(i,j)
          enddo
        enddo
      enddo

c
c advance the short-timestep horizontal dynamics
c
c initialize the odd state
      If (OTIDE > 0)  Then
         TIME = DTSRC*ITIME + .5*DTOFS  !  seconds since 2000/01/01/00
     +        + DTO*((NO-1)*NDYNO/NOCEAN)
         Call OTIDEW (TIME)
         Call OTIDEV (.False.,DTOFS, UO2,VO2,UOD2,VOD2)  ;  EndIf
      Call ODHORZ(MO ,UO ,VO ,UOD ,VOD ,OPBOT ,
     &            MO2,UO2,VO2,UOD2,VOD2,OPBOT2, DTOFS,.false.)
      If (OTIDE > 0)  Then
         TIME = DTSRC*ITIME + .5*DTO
     +        + DTO*((NO-1)*NDYNO/NOCEAN)
         Call OTIDEW (TIME)
         Call OTIDEV (.False.,DTO, UO1,VO1,UOD1,VOD1)  ;  EndIf
      Call ODHORZ(MO2,UO2,VO2,UOD2,VOD2,OPBOT2,
     &            MO1,UO1,VO1,UOD1,VOD1,OPBOT1, DTO,.false.)
c loop over the leapfrog steps
      neven = NDYNO / (2*NOCEAN)
      do n=1,neven
c update the even state         
        If (OTIDE > 0)  Then
           TIME = DTSRC*ITIME + .5*DTOLF
     +          + DTO*((NO-1)*NDYNO/NOCEAN + 2*N-2)
           Call OTIDEW (TIME)
           Call OTIDEV (.True.,DTOLF, UO,VO,UOD,VOD)  ;  EndIf
        Call ODHORZ(MO1,UO1,VO1,UOD1,VOD1,OPBOT1,
     &              MO ,UO ,VO ,UOD ,VOD ,OPBOT , DTOLF,.true.)
        if(n == neven) exit ! no need to further update the odd state
c update the odd state
        If (OTIDE > 0)  Then
           TIME = DTSRC*ITIME + .5*DTOLF
     +          + DTO*((NO-1)*NDYNO/NOCEAN + 2*N-1)
           Call OTIDEW (TIME)
           Call OTIDEV (.False.,DTOLF, UO1,VO1,UOD1,VOD1)  ;  EndIf
        Call ODHORZ(MO ,UO ,VO ,UOD ,VOD ,OPBOT ,
     &              MO1,UO1,VO1,UOD1,VOD1,OPBOT1, DTOLF,.false.)
      enddo

c is this still needed?
      if(j_1 == JM) UO(IVNP,JM,:) = VONP(:) ! not needed if Mod(IM,4)=0

c
c long-timestep vertical redistribution of mass (and momentum)
c
      Call OFLUXV

c
c long-timestep advection of potential enthalpy, salt, and tracers
c
      dtdum = DTOLF
      Call CONSERV_OCE (G0INIT)
      if(use_qus==1) then
        CALL OADVT3 (MO1,G0M,
     &     GXMO,GYMO,GZMO, GXXMO,GYYMO,GZZMO, GXYMO,GYZMO,GZXMO,
     &     DTDUM,.FALSE.,OIJL(1,J_0H,1,IJL_GFLX))
        CALL OADVT3 (MO1,S0M,
     &     SXMO,SYMO,SZMO, SXXMO,SYYMO,SZZMO, SXYMO,SYZMO,SZXMO,
     &     DTDUM,.TRUE.,OIJL(1,J_0H,1,IJL_SFLX))
      else
        CALL OADVT2 (MO1,G0M,GXMO,GYMO,GZMO,DTDUM,.FALSE.
     *        ,OIJL(1,J_0H,1,IJL_GFLX))
        CALL OADVT2 (MO1,S0M,SXMO,SYMO,SZMO,DTDUM,.TRUE.
     *        ,OIJL(1,J_0H,1,IJL_SFLX))
      endif
      Call CONSERV_OCE (G0FINAL)

c
c diagnostics
c
      DO L=1,LMO
        do j=j_0s,j_1s
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              OIJL(I,J,L,IJL_MFU) = OIJL(I,J,L,IJL_MFU) +
     &             SMU(I,J,L)*DTOLF
            enddo
          enddo
        enddo
        do j=j_0,j_1
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              OIJL(I,J,L,IJL_MFV) = OIJL(I,J,L,IJL_MFV) +
     &             SMV(I,J,L)*DTOLF
            enddo
          enddo
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              OIJL(I,J,L,IJL_MFW) = OIJL(I,J,L,IJL_MFW) +
     &             SMW(I,J,L)*DTOLF
              OIJL(I,J,L,IJL_MFW2)= OIJL(I,J,L,IJL_MFW2)+
     &             byno*(SMW(I,J,L)/neven)**2
              OIJL(I,J,L,IJL_MO)  = OIJL(I,J,L,IJL_MO) +  MO(I,J,L)*byno
              OIJL(I,J,L,IJL_G0M) = OIJL(I,J,L,IJL_G0M) +G0M(I,J,L)*byno
              OIJL(I,J,L,IJL_S0M) = OIJL(I,J,L,IJL_S0M) +S0M(I,J,L)*byno
            enddo
          enddo
        enddo
      ENDDO

      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            OIJ(I,J,IJ_SSH) = OIJ(I,J,IJ_SSH) + OGEOZ(I,J)*byno
            OIJ(I,J,IJ_PB)  = OIJ(I,J,IJ_PB)  +
     &           (OPBOT(I,J)-ZE(LMM(I,J))*RHOWS*GRAV)*byno
            OIJ(I,J,IJ_dEPO_Dyn) = OIJ(I,J,IJ_dEPO_Dyn) +
     +           (G0FINAL(I,J) - G0INIT(I,J))
          enddo
        enddo
      enddo


#ifdef TRACERS_OCEAN
      if(ntrtrans.gt.1) then
        ! accumulate mass fluxes
        asmu = asmu + smu
        asmv = asmv + smv
        asmw = asmw + smw
      else
        if(use_qus==1) then
          DO N=1,tracerlist%getsize()
            entry=>tracerlist%at(n)
            CALL OADVT3(MO1,TRMO(1,J_0H,1,N),
     &       TXMO (1,J_0H,1,N),TYMO (1,J_0H,1,N),TZMO (1,J_0H,1,N),
     &       TXXMO(1,J_0H,1,N),TYYMO(1,J_0H,1,N),TZZMO(1,J_0H,1,N),
     &       TXYMO(1,J_0H,1,N),TYZMO(1,J_0H,1,N),TZXMO(1,J_0H,1,N),
     &       dtdum,entry%t_qlimit,TOIJL(1,J_0H,1,TOIJL_TFLX,N))
          ENDDO
        else
          DO N=1,tracerlist%getsize()
            entry=>tracerlist%at(n)
            CALL OADVT2(MO1,TRMO(1,J_0H,1,N),TXMO(1,J_0H,1,N)
     *           ,TYMO(1,J_0H,1,N),TZMO(1,J_0H,1,N),dtdum,entry%t_qlimit
     *           ,TOIJL(1,J_0H,1,TOIJL_TFLX,N))
          ENDDO
        endif
      endif

c
c tracer diagnostics
c
c todo for ntrtrans>1 case: sample toijl right after long-step transport
c rather than here?   Or scale here by mo/motr.
        DO N=1,tracerlist%getsize()
          DO L=1,LMO
            TOIJL(:,:,L,TOIJL_CONC,N)=TOIJL(:,:,L,TOIJL_CONC,N)
     *           +TRMO(:,:,L,N)*byno
          END DO
        END DO
#endif

        CALL CHECKO ('OADVT ')

      if(use_qus==1) then
        ! Save seawater mass pre-straits to adjust second-order
        ! moments (SOMs) post-straits.  This is easier than
        ! propagating the SOMs through the straits code.
        ! The need for adjustment primarily arises from the
        ! fact that the moments have extensive units.
        ! First-order moments are already propagated through
        ! the straits code.
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          mo1(i,j,l) = mo(i,j,l)
        enddo
        enddo
        enddo
        enddo
      endif
      call gather_ocean_straits()

      IF(AM_I_ROOT()) THEN
C****
C**** Acceleration and advection of tracers through ocean straits
C****
          CALL STPGF(DTS/NOCEAN)
          CALL STADV(DTS/NOCEAN)
          CALL CHECKO_serial ('STADV0')
        IF (NO .EQ. NOCEAN) THEN
          CALL STCONV
          CALL STBDRA
        END IF
      END IF
      call scatter_ocean_straits()
      call BCAST_straits (.false.)

      if(use_qus==1) then
        ! Adjust second-order moments post-straits in upstream cells.
        ! Seawater emerging from a strait (downstream) is assumed
        ! to have zero SOMs for now.
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          if(mo1(i,j,l) .le. mo(i,j,l)) cycle
          mrat_st = mo(i,j,l)/mo1(i,j,l)
          gxxmo(i,j,l) = gxxmo(i,j,l)*mrat_st
          gyymo(i,j,l) = gyymo(i,j,l)*mrat_st
          gzzmo(i,j,l) = gzzmo(i,j,l)*mrat_st
          gxymo(i,j,l) = gxymo(i,j,l)*mrat_st
          gyzmo(i,j,l) = gyzmo(i,j,l)*mrat_st
          gzxmo(i,j,l) = gzxmo(i,j,l)*mrat_st
          sxxmo(i,j,l) = sxxmo(i,j,l)*mrat_st
          syymo(i,j,l) = syymo(i,j,l)*mrat_st
          szzmo(i,j,l) = szzmo(i,j,l)*mrat_st
          sxymo(i,j,l) = sxymo(i,j,l)*mrat_st
          syzmo(i,j,l) = syzmo(i,j,l)*mrat_st
          szxmo(i,j,l) = szxmo(i,j,l)*mrat_st
#ifdef TRACERS_OCEAN
          txxmo(i,j,l,:) = txxmo(i,j,l,:)*mrat_st
          tyymo(i,j,l,:) = tyymo(i,j,l,:)*mrat_st
          tzzmo(i,j,l,:) = tzzmo(i,j,l,:)*mrat_st
          txymo(i,j,l,:) = txymo(i,j,l,:)*mrat_st
          tyzmo(i,j,l,:) = tyzmo(i,j,l,:)*mrat_st
          tzxmo(i,j,l,:) = tzxmo(i,j,l,:)*mrat_st
#endif
        enddo
        enddo
        enddo
        enddo
      endif

#ifdef TRACERS_OCEAN
      if(ntrtrans.gt.1) then
        ! increment the seawater mass array corresponding to tracer state
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          motr(i,j,l) = motr(i,j,l) + (mo(i,j,l)-mo1(i,j,l))
        enddo
        enddo
        enddo
        enddo
      endif
#endif

        CALL CHECKO ('STADV ')

      ENDDO  !  End of Do-loop NO=1,NOCEAN

        CALL TIMER (NOW,MDYNO)
        IF (ATMOCN%MODD5S == 0) CALL DIAGCO (12,atmocn)

c
c recalculate vbar etc. for ocean physics
c
      CALL ODHORZ0

      call ocnstate_derived

C**** Apply Wajsowicz horizontal diffusion to UO and VO ocean currents
C**** every 3 hours
      dt_odiff = 3.*3600.
      if(mod(itime,int(dt_odiff/dts)).eq.0) then
        CALL ODIFF(dt_odiff)
      endif
c      CALL OABFILx ! binomial filter
c      CALL OABFILy ! binomial filter
c      CALL CHECKO ('ODIFF0')


#ifdef TRACERS_OCEAN
      do_tracer_trans = mod(1+itime-itimei,ntrtrans).eq.0
#endif

C****
C**** Mesoscale tracer transports
C****
      call ocnmeso_drv

#ifdef TRACERS_OCEAN

      ! set motr for source-term routines below that refer to it
      if(ntrtrans.eq.1) then
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          motr(i,j,l) = mo(i,j,l)
        enddo
        enddo
        enddo
        enddo

C****
C**** Resolved-flow tracer transports
C****
      elseif(do_tracer_trans) then
        ! This is a tracer transport timestep.  Note that in contrast
        ! to the ntrtrans==1 case, the resolved advection is
        ! performed after the mesoscale transport.

        ! First copy accumulated mass fluxes into arrays used by OADVT
        ! and reset accumulators.
        smu = asmu
        smv = asmv
        smw = asmw
        call halo_update(grid,smv,from=south)
        asmu = 0.
        asmv = 0.
        asmw = 0.
        do l=1,lmo
        do j=j_0,j_1
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              if(motr(i,j,l).lt.0d0) then
                call stop_model('motr < 0 aft meso',255)
              endif
              mmi(i,j,l) = motr(i,j,l)*dxypo(j)
            enddo
          enddo
        enddo
        enddo

        if(use_qus.ne.1) then
          call stop_model('OADVT2 not ready for ntrtrans>1',255)
        endif

        ! Perform the advection.
        DO N=1,tracerlist%getsize()
          entry=>tracerlist%at(n)
          CALL OADVT3(MO1,TRMO(1,J_0H,1,N),
     &       TXMO (1,J_0H,1,N),TYMO (1,J_0H,1,N),TZMO (1,J_0H,1,N),
     &       TXXMO(1,J_0H,1,N),TYYMO(1,J_0H,1,N),TZZMO(1,J_0H,1,N),
     &       TXYMO(1,J_0H,1,N),TYZMO(1,J_0H,1,N),TZXMO(1,J_0H,1,N),
     &       dtdum,entry%t_qlimit,TOIJL(1,J_0H,1,TOIJL_TFLX,N))
        ENDDO

        ! Save post-advection seawater mass corresponding to tracer state.
        ! To prevent buildup of roundoff error, then reset motr to mo.
        do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          motr(i,j,l) = mo1(i,j,l)/dxypo(j)
          ! sanity check
          !if(abs(motr(i,j,l)-mo(i,j,l)) .gt. 1d-11*mo(i,j,l)) then
          !  write(6,*) 'motr error ',motr(i,j,l),mo(i,j,l)
          !endif
          motr(i,j,l) = mo(i,j,l)
        enddo
        enddo
        enddo
        enddo
      endif
#endif

#ifdef TRACERS_GASEXCH_ocean_CO2
         Call CARBON ('OCNMESO')
         Call NITR ('OCNMESO')
#endif
      CALL CHECKO ('GMDIFF')

#ifdef TRACERS_OCEAN
      CALL OC_TDECAY(DTS)
      if (n_age.gt.0) CALL OCN_TR_AGE(DTS)
      if (n_vent.gt.0) CALL OCN_TR_VENT(DTS)
      if (n_gasx.gt.0) CALL OCN_TR_GASX(DTS)
      if ( n_wms1+n_wms2+n_wms3 .gt. 0)
           ! if-test here is not really needed, since
           ! tests of individual n_wms>0 are within the routine
     &     CALL OCN_TR_WaterMass(DTS)

      if (n_ocfc11.gt.0) CALL OCN_TR_CFC(DTS,11)   
      if (n_ocfc12.gt.0) CALL OCN_TR_CFC(DTS,12)   
      if (n_sf6.gt.0) CALL OCN_TR_CFC(DTS,6)   
!     if (n_abioDIC.gt.0) is defined in obio_carbon.f
#endif

      CALL TIMER (NOW,MSGSO)

c-------------------------------------------------------------------
c End ocean-processors-only code region
c-------------------------------------------------------------------

      else ! no ocean domain, call the timers anyway
        CALL TIMER (NOW,MSURF)
        CALL TIMER (NOW,MSGSO)
        CALL TIMER (NOW,MDYNO)
        CALL TIMER (NOW,MSGSO)
      endif ocean_processors_only

C***  Get the data from the ocean grid to the atmospheric grid
      CALL TOC2SST(atmocn)
      call OG2AG_oceans(iceocn)

C***  Interpolate ocean surface velocity to the DYNSI grid
      call OG2IG_uvsurf(dynsice,atmocn)

      RETURN
      END SUBROUTINE OCEANS

      Subroutine OFLUXV
      use constant, only : grav
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM,LMOU=>LMU,LMOV=>LMV,
     &     ZOE=>ZE,dZO, DTS, DXYPO, BYDXYPO, OPRESS, dtolf
      Use OCEAN, only : opbot,mo,uo,vo
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      Use OCEAN_DYN, Only: SMW
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH,SOUTH
      USE OCEANR_DIM, only : grid=>ogrid
      Implicit None
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: msum
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mb,mtmp,mwtmp
      Real*8 :: mwfac,mfinal
      Integer*4 I,J,L,LM, J1,JN,J1P,J1H,JNP,JNQ,JNH, N
      Logical QNP
C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNQ = Min(JN,JM-2)    !    4      8     JM-2   Exclude NP,NP-1
      JNH = Min(JN+1,JM)
      QNP = JN==JM          !    F      F      T
C****

      do l=1,lmo
        do j=j1p,jn
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              mb(i,j,l) = mo(i,j,l)
            enddo
          enddo
        enddo
        if(qnp) then
          j = jm
          do i=2,im
            mb(i,j,l) = mo(i,j,l)
          enddo
        endif
      enddo

      do j=j1p,jn !jnh
        mwfac = dxypo(j)/dtolf
        do n=1,nbyzm(j,2)
          do i=i1yzm(n,j,2),i2yzm(n,j,2)
            msum(i,j) = (opbot(i,j)-opress(i,j))/grav
            lm=lmom(i,j)
            mfinal = msum(i,j)*dzo(1)/zoe(lm)
c            mw(i,j) = mwfac*(mo(i,j,1)-mfinal)
            smw(i,j,1) = mwfac*(mo(i,j,1)-mfinal)
            mo(i,j,1) = mfinal
          enddo
        enddo
      enddo

      do l=2,lmo-1

      do j=j1p,jn !jnh
        mwfac = dxypo(j)/dtolf
        do n=1,nbyzm(j,l+1)
          do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
            lm=lmom(i,j)
            mfinal = msum(i,j)*dzo(l)/zoe(lm)
c            mw(i,j) = mw(i,j) + mwfac*(mo(i,j,l)-mfinal)
            smw(i,j,l) = smw(i,j,l-1) + mwfac*(mo(i,j,l)-mfinal)
            mo(i,j,l) = mfinal
          enddo
        enddo
      enddo

      enddo ! l

c
c update bottom layer mass
c
      do j=j1p,jn !jnh
        do n=1,nbyzm(j,2)
          do i=i1yzm(n,j,2),i2yzm(n,j,2)
            lm=lmom(i,j)
            mfinal = msum(i,j)*dzo(lm)/zoe(lm)
            mo(i,j,lm) = mfinal
          enddo
        enddo
      enddo

      if(qnp)  then ! fill pole
        j = jm
        do l=1,lmom(1,j)
          do i=2,im
            mo(i,j,l) = mo(1,j,l)
            smw(i,j,l) = smw(1,j,l)
          enddo
        enddo
      endif

c
c for now, vertical advection of U,V uses the simplest upstream scheme
c
      call halo_update(grid,smw,from=north)
      call halo_update(grid,mb,from=north)
      do l=1,lmo
        do j=j1p,jnp
          mwfac = bydxypo(j)
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),min(i2yzu(n,j,l),im-1)
              mtmp(i,j,l) = .5*(mb(i,j,l)+mb(i+1,j,l))
              mwtmp(i,j,l) = .5*(smw(i,j,l)+smw(i+1,j,l))*mwfac
            enddo
            i = im
            if(l <= lmou(i,j)) then
              mtmp(i,j,l) = .5*(mb(i,j,l)+mb(1,j,l))
              mwtmp(i,j,l) = .5*(smw(i,j,l)+smw(1,j,l))*mwfac
            endif
          enddo
        enddo
      enddo
      do j=j1p,jnp
        do n=1,nbyzu(j,1)
          do i=i1yzu(n,j,1),i2yzu(n,j,1)
            mwtmp(i,j,lmou(i,j)) = 0.
          enddo
        enddo
      enddo
      call oadvuz(uo,mtmp,mwtmp,dtolf,j1p,jnp,nbyzu,i1yzu,i2yzu)
      do l=1,lmo
        do j=j1p,jnp
          do n=1,nbyzv(j,l)
            do i=i1yzv(n,j,l),i2yzv(n,j,l)
              mtmp(i,j,l) = .5*(mb(i,j,l)+mb(i,j+1,l))
              mwtmp(i,j,l) =
     &             .5*(smw(i,j,l)*bydxypo(j)+smw(i,j+1,l)*bydxypo(j+1))
            enddo
          enddo
        enddo
      enddo
      do j=j1p,jnp
        do n=1,nbyzv(j,1)
          do i=i1yzv(n,j,1),i2yzv(n,j,1)
            mwtmp(i,j,lmov(i,j)) = 0.
          enddo
        enddo
      enddo
      call oadvuz(vo,mtmp,mwtmp,dtolf,j1p,jnp,nbyzv,i1yzv,i2yzv)

      Return
      End Subroutine OFLUXV

      module opfil2_coeffs
      use constant
      use ocean, only : im,lmo
      implicit none
      save
! abbreviation: MBF(S) = matrix-based filtering (segment)
!@param hwid_max maximum half-width of MBF stencil
      integer, parameter :: hwid_max=15
!@var nbas_fft number of segments over which to apply fft-based filter
!@var nbas_fil number of segments over which to apply matrix-based filter
      integer :: nbas_fft,nbas_fil
!@var jfft the nth value of jfft is the j index of the nth fft filter segment
!@var jfil the nth value of jfil is the j index of the nth MBFS
!@var i1fil,i2fil the nth value of i1fil,i2fil is the left/right longitude index
!@+   of the nth MBFS
      integer, dimension(:), allocatable ::
     &     jfft,jfil,i1fil,i2fil,indx_fil
!@var n1fft,n2fft the lth value of n1fft,n2fft is the first/last segment index for
!@+   fft-based filtering at layer l
!@var n1fil,n2fil the lth value of n1fil,n2fil is the first/last segment index for
!@+   MBF at layer l
      integer, dimension(lmo) :: n1fft,n2fft,n1fil,n2fil
!@var reduco is the concatenation of the filtering matrices needed for the local
!@+   domain (including halo rows)
      real*8, dimension(:), allocatable :: reduco

! the following are used for fft-based filtering
      integer, parameter :: IMz2=IM/2  !  highest possible wave number
      integer, allocatable :: NMIN(:) !  minimum wave number for Fourier smoothing
      real*8, allocatable  :: SMOOTH(:,:)  !  multiplies Fourier coefficient

      contains
      subroutine calc_opfil2_coeffs
      use oceanr_dim, only : grid=>ogrid
      use ocean, only : jm,lmo,dxpo,dypo,lmu
      use ocean, only : nbyzmax,nbyzu,i1yzu,i2yzu
      implicit none
C****
C**** The ocean dynamics time step is chosen as the largest convenient
C**** time step such that gravity waves of length 2*dY do not cause
C**** the numerical solution of the momentum equation to diverge.
C**** Shorter gravity waves that are resolved at high latitudes in the
C**** zonal direction could cause the numerical solution to diverge.
C**** To prevent this, certain fields are spectrally analyzed and
C**** coefficients of waves shorter than 2*dY are multiplied by the
C**** wave length divided by 2*dY.
C****
C**** For an ocean basin that is IWID grid cells wide, the fields of
C**** interest are defined on IWID-1 grid cell edges.  The fields are
C**** spectrally analyzed on 2*IWID periodic values which include 0's
C**** on each coast and IWID-1 reflected values of opposite sign.
C**** The cosine spectral coefficients are all zero, but the sine
C**** coefficients for wave numbers 1 to IWID-1 are generally nonzero.
C**** The shortest wave (wave number IWID-1) in a basin that is IWID
C**** cells wide has length 2*IWID*dX/(IWID-1).
C****
C**** IWMIN is calculated such that the length of the shortest wave in this
C**** minimum basin is 2*IWMIN*DXP(J)/(IWMIN-1) which must be less than 2*DYP(3)
C****
C**** Since filtered field values are fixed linear combinations of
C**** the IWID-1 unfiltered values, the matrix coefficents of linear
C**** combinations are calculated once in this program, but are used
C**** repetively in the ocean polar filter subroutine.  The matrix
C**** coefficients depend upon IWID and latitude (dX).
C****
C**** The polar filter is not applied at the poles (J = 1 or JM) nor
C**** at latitudes for which the aspect ratio dY/dX < 1.
C****

      integer :: i,j,k,l,n,nn,ja,indx,iwide,km,iw,ie,iwmin,iter,ibas
      real*8 reduc(im-1,im-1)   !  single polar filter reduction matrix
      real*8 sintab(im-1)
      real*8 :: reducn,drat
      integer, dimension(im-1) :: k1,k2
      integer, allocatable :: indx_sv(:,:)
      integer :: nbas,i1bas(nbyzmax),i2bas(nbyzmax)
      integer :: js0,js1,jtmp,js0a,js1a
C****

      js0 = max(2   ,grid%j_strt_halo)
      js1 = min(jm-1,grid%j_stop_halo)

C**** Calculate factors for fft-based filtering
      allocate(nmin(js0:js1))
      allocate(smooth(imz2,js0:js1))
      do j=js0,js1
        drat = dxpo(j)/dypo(3)
        do n=imz2,1,-1
          smooth(n,j) = imz2*drat/n
          if (smooth(n,j) >= 1d0) exit
        enddo
        nmin(j) = n+1
      enddo

C**** Calculate factors for matrix-based filtering
      n1fft(:) = huge(n1fft(1))
      n2fft(:) = -1
      n1fil(:) = huge(n1fil(1))
      n2fil(:) = -1

      js0a = min(js0, jm+1-js0) !  absolute latitude index
      js1a = min(js1, jm+1-js1) !  absolute latitude index
      if(js0a > js1a) then
        jtmp = js0a
        js0a = js1a
        js1a = jtmp
      endif
      if(js0 <= jm/2 .and. js1 > jm/2) js1a = jm/2
      allocate(indx_sv(im-1,js0a:js1a))

      do iter=1,2 ! two passes: first is for counting array sizes
      indx_sv(:,:) = -1
      nbas_fft = 0
      nbas_fil = 0
      indx = 0
      do l=1,lmo
      do j=js0,js1
        if(dxpo(j) >= dypo(3)) cycle
        ja = min(j, jm+1-j)  !  absolute latitude index

        iwmin = ceiling (dypo(3) / (dypo(3)-dxpo(ja)))

        if(nbyzu(j,l)  .eq.1 .and.
     &     i1yzu(1,j,l).eq.1 .and.
     &     i2yzu(1,j,l).eq.im ) then ! all ocean, use fft-based filtering
          n = nbas_fft + 1
          nbas_fft = n
          if(iter.eq.1) cycle ! first pass for counting only
          jfft(n) = j
          n1fft(l) = min(n1fft(l),n)
          n2fft(l) = max(n2fft(l),n)
        else                         ! some land,  use matrix-based filtering
          if(l <= lmu(im,j) .and. i1yzu(1,j,l).eq.1) then ! IDL crossing
            nbas = nbyzu(j,l)-1
            if(nbas.gt.1) then
              i1bas(1:nbas-1) = i1yzu(2:nbas,j,l)
              i2bas(1:nbas-1) = i2yzu(2:nbas,j,l)
            endif
            i1bas(nbas) = i1yzu(nbyzu(j,l),j,l)
            i2bas(nbas) = i2yzu(1,j,l) + im
          else
            nbas = nbyzu(j,l)
            i1bas(:) = i1yzu(:,j,l)
            i2bas(:) = i2yzu(:,j,l)
          endif
          n = 0
          do ibas=1,nbas ! loop over basins
            iw = i1bas(ibas)
            ie = i2bas(ibas)
            iwide = ie-iw+2
            if (iwide < iwmin) cycle
            n = nbas_fil+1
            nbas_fil = n
            if(indx_sv(iwide,ja).lt.0) then
              indx_sv(iwide,ja) = indx
              do i=1,iwide-1
                k1(i) = max(      1,i-hwid_max)
                k2(i) = min(iwide-1,i+hwid_max)
              enddo
              if(iter.eq.1) then ! first pass for counting only
                do i=1,iwide-1
                  indx=indx+1+k2(i)-k1(i)
                enddo
              else
                km = 2*iwide
                reduc(:,:) = 0
                do nn=iwide-1,1,-1
                  reducn = (1 - dxpo(ja)*iwide/(dypo(3)*nn))*4/km
                  if (reducn <= 0) cycle
                  do i=1,iwide-1
                    sintab(i) = sin(twopi*nn*i/km)
                  enddo
                  do i=1,iwide-1
                  do k=k1(i),k2(i)
                    reduc(k,i) = reduc(k,i) +
     &                   sintab(k)*(sintab(i)*reducn)
                  enddo
                  enddo
                enddo
                do i=1,iwide-1
                do k=k1(i),k2(i)
                  indx=indx+1
                  reduco(indx) = reduc(k,i)
                enddo
                enddo
              endif
            endif
            if(iter.eq.1) cycle ! first pass for counting only
            jfil(n) = j
            i1fil(n) = iw
            i2fil(n) = ie
            indx_fil(n) = indx_sv(iwide,ja)
            n1fil(l) = min(n1fil(l),n)
            n2fil(l) = max(n2fil(l),n)
          enddo ! ibas
        endif
      enddo ! j
      enddo ! l
      if(iter.eq.1) then ! allocate once we know the counts
        allocate(jfft(nbas_fft))
        allocate(jfil(nbas_fil))
        allocate(i1fil(nbas_fil),i2fil(nbas_fil),indx_fil(nbas_fil))
        allocate(reduco(indx))
      endif
      enddo ! iter

      deallocate(indx_sv)

      end subroutine calc_opfil2_coeffs

      end module opfil2_coeffs

      subroutine opfil2(x,l,jmin,jmax)
C****
C**** OPFIL smoothes X in the zonal direction by reducing coefficients
C**** of its Fourier series for high wave numbers near the poles.
C****
      use ocean, only: im,jm
      use oceanr_dim, only : grid=>ogrid
      use opfil2_coeffs
      implicit none
C****
      real*8,intent(inout) :: x(im,grid%j_strt_halo:grid%j_stop_halo)
      integer*4, intent(in) :: l,jmin,jmax
      integer :: i,i1,i2,j,k,k1,k2,indx,n,nn
      real*8 AN(0:IMz2), !  Fourier cosine coefficients
     *       BN(0:IMz2), !  Fourier sine coefficients
     *          Y(IM*2), !  original copy of X that wraps around IDL
     *       REDUC
C****
      integer, save :: ifirst=1
      integer :: j1,jn,j1p,jnp

      if(ifirst == 1) then
        ifirst = 0
        call calc_opfil2_coeffs
      endif

C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP

C****
C**** FFT-based filter for all-ocean latitudes
C****
      do nn=n1fft(l),n2fft(l)
        j = jfft(nn)
        if(j < jmin .or. j > jmax) cycle
        call offt (x(1,j),an,bn)
        do n=nmin(j),imz2-1
          an(n) = an(n)*smooth(n,j)
          bn(n) = bn(n)*smooth(n,j)
        enddo
        an(imz2) = an(imz2)*smooth(imz2,j)
        call offti (an,bn,x(1,j))
      enddo

C****
C**** Per-basin filter
C****
      do n=n1fil(l),n2fil(l)
        j = jfil(n)
        if(j < jmin .or. j > jmax) cycle
        indx = indx_fil(n)
        i1 = i1fil(n)
        i2 = i2fil(n)
        if(i2 <= im)  then
C**** Ocean basin does not wrap around the IDL.
C**** Copy X to temporary array Y and filter X in place.
          do i=i1,i2
            y(i) = x(i,j)
          enddo
          do i=i1,i2
            reduc = 0
            k1 = max(i1,i-hwid_max)
            k2 = min(i2,i+hwid_max)
            do k=k1,k2
              indx=indx+1
              reduc = reduc + reduco(indx)*y(k)
            enddo
            x(i,j) = x(i,j) - reduc
          enddo
        else
C**** Ocean basin wraps around the IDL.
C**** Copy X to temporary array Y and filter X in place.
          do i=i1,im
            y(i) = x(i,j)
          enddo
          do i=im+1,i2
            y(i) = x(i-im,j)
          enddo
          do i=i1,im
            reduc = 0
            k1 = max(i1,i-hwid_max)
            k2 = min(i2,i+hwid_max)
            do k=k1,k2
              indx=indx+1
              reduc = reduc + reduco(indx)*y(k)
            enddo
            x(i,j) = x(i,j) - reduc
          enddo
          do i=1,i2-im
            reduc = 0
            k1 = max(i1,im+i-hwid_max)
            k2 = min(i2,im+i+hwid_max)
            do k=k1,k2
              indx=indx+1
              reduc = reduc + reduco(indx)*y(k)
            enddo
            x(i,j) = x(i,j) - reduc
          enddo
        endif  ! wraparound or not        
      enddo

      return
      end subroutine opfil2

      Subroutine ODHORZ(MOH,UOH,VOH,UODH,VODH,OPBOTH,
     &                  MO ,UO ,VO ,UOD ,VOD ,OPBOT, DT,qeven)
      Use CONSTANT, Only: GRAV,omega,UNDEF_VAL
      Use OCEAN, Only: IM,JM,LMO,
     *                 LMOM=>LMM,LMOU=>LMU,LMOV=>LMV,
     *                 mZSOLID=>HOCEAN,OGEOZ,
     &     SINVO,SINPO,DXPO,DYPO,DXYVO,DXVO,DYVO,DXYPO
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv,nbyzc,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv, i1yzc,i2yzc
      Use OCEAN_DYN, Only: SMU,SMV,VBAR,dZGdP,MMI
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH, SOUTH
     &     ,GLOBALMAX
      USE OCEANR_DIM, only : grid=>ogrid
      Implicit None
      Real*8, Intent(In),
     &     Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     &     MOH,UOH,VOH,UODH,VODH
      Real*8, Intent(In),
     &     Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OPBOTH
      Real*8, Intent(Inout),
     &  Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     &     MO,UO,VO,UOD,VOD
      Real*8, Intent(Inout),
     &     Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) :: OPBOT
      Real*8,Intent(In) :: DT
      Logical, Intent(In) :: qeven
c
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     ZG,P,PDN,DH,USMOOTH,PGFX,PGFY,MU,MV,UA,VA,VORT,KE
      Real*8 corofj,mufac,mvfac,bydx,bydy,mmid,xeven,
     &     convij,convfac,dp,pgf4pt,pgfac,uq,vq,uasmooth
      Integer*4 I,J,L, J1,JN,J1P,JNP,J1H,JNH, j1a, jminpfu,jmaxpfu, N
      integer :: smallmo_loc,smallmo
      Logical QNP

C****
C**** Extract domain decomposition band parameters
C**** ASSUME THAT J1 NEVER EQUALS 2 AND JN NEVER EQUALS JM-1
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      J1P = Max(J1,2)       !    2      5     JM-3   Exclude SP
      JNP = Min(JN,JM-1)    !    4      8     JM-1   Exclude NP
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
      QNP = JN==JM          !    F      F      T
      JMINpfu = max(2,J1-1)
      JMAXpfu = min(JNH,JM-1)
      j1h = max(1,grid%j_strt_halo)
      j1a = grid%j_strt_halo

      if(qeven) then
        xeven = 1.
      else
        xeven = 0.
      endif

      ZG(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      P(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      PDN(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      DH(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      USMOOTH(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      UA(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      VA(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL
      KE(:,[GRID%J_STRT_HALO,GRID%J_STOP_HALO])=UNDEF_VAL

c
c initialize pressure and geopotential at the ocean bottom
c
      do j=j1h,jnh
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            pdn(i,j) = opboth(i,j)
            OGEOZ(I,J) = - mZSOLID(I,J)*GRAV
          enddo
        enddo
      enddo

! zero only once with bottom-up looping. see whether these are
! all necessary
      usmooth(:,:) = 0.
      mu(:,:) = 0.
      mv(:,:) = 0.
      pgfx(:,:) = 0.
      pgfy(:,:) = 0.
      vort(:,:) = 0.

      smallmo_loc = 0
c
c loop over layers, starting at the ocean bottom
c
      do l=lmo,1,-1

C**** Apply polar filter to West-East velocity
      if(any(nbyzu(jminpfu:jmaxpfu,l).gt.0)) then
        Do J=jminpfu,jmaxpfu
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              usmooth(I,J) = UOH(I,J,L)
            enddo
          enddo
        enddo
        Call OPFIL2(usmooth,L,jminpfu,jmaxpfu)
      endif
      if(qnp) then
        usmooth(:,jm) = uoh(:,jm,l)
      endif

c
c calculate pressure, geopotential, and kinetic energy
c
      do j=j1h,jnh
        i = 1
        if(l <= lmom(i,j)) then
          dp = moh(i,j,l)*grav
          dh(i,j) = moh(i,j,l)*vbar(i,j,l)
          p(i,j) = pdn(i,j) - .5*dp
          zg(i,j) = ogeoz(i,j) + dp*.5*dzgdp(i,j,l)
          pdn(i,j) = pdn(i,j) - dp
          ogeoz(i,j) = ogeoz(i,j) + dh(i,j)*grav
        endif
        do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            dp = moh(i,j,l)*grav
            dh(i,j) = moh(i,j,l)*vbar(i,j,l)
            p(i,j) = pdn(i,j) - .5*dp
            zg(i,j) = ogeoz(i,j) + dp*.5*dzgdp(i,j,l)
            pdn(i,j) = pdn(i,j) - dp
            ogeoz(i,j) = ogeoz(i,j) + dh(i,j)*grav
          enddo
        enddo
      enddo
      do j=j1,jnh
        i = 1
        if(l <= lmom(i,j)) then
          uasmooth = .5*(usmooth(im,j)+usmooth(i,j))
          ua(i,j) = .5*(uoh(im,j,l)+uoh(i,j,l))
          va(i,j) = .5*(voh(i,j-1,l)+voh(i,j,l))
          ke(i,j) = .5*(ua(i,j)*uasmooth + va(i,j)**2)
        endif
        do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            uasmooth = .5*(usmooth(i-1,j)+usmooth(i,j))
            ua(i,j) = .5*(uoh(i-1,j,l)+uoh(i,j,l))
            va(i,j) = .5*(voh(i,j-1,l)+voh(i,j,l))
            ke(i,j) = .5*(ua(i,j)*uasmooth + va(i,j)**2)
          enddo
        enddo
      enddo

c
c fill pole
c
      if(qnp) then
        if(l <= lmom(1,jm)) then
          j = jm
          do i=2,im
            dh(i,j) = dh(1,j)
            p(i,j) = p(1,j)
            zg(i,j) = zg(1,j)
            ke(i,j) = ke(1,j)
            ua(i,j) = .5*(usmooth(i-1,j)+usmooth(i,j))
          enddo
        endif
      endif

c
c store, and smooth, the east-west pressure gradient
c
      DO J=J1P,JMAXpfu
        mufac = .5*dypo(j)
        bydx = 1d0/dxpo(j)
        do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
            mmid = MOH(I,J,L)+MOH(I+1,J,L)
            MU(I,J) = mufac*usmooth(i,j)*mmid
            SMU(I,J,L) = SMU(I,J,L) + MU(I,J)*xeven
            pgfx(i,j) =
     &           (ZG(I,J)-ZG(I+1,J))+(P(I,J)-P(I+1,J))*
     &           (dH(I,J)+dH(I+1,J))/mmid
            pgfx(i,j) = pgfx(i,j)*bydx
          enddo
        enddo
        I=IM
        if(l <= lmou(i,j)) then
          mmid = MOH(I,J,L)+MOH(1,J,L)
          MU(I,J) = mufac*usmooth(i,j)*mmid
          SMU(I,J,L) = SMU(I,J,L) + MU(I,J)*xeven
          pgfx(i,j) =
     &         (ZG(I,J)-ZG(1,J))+(P(I,J)-P(1,J))*
     &         (dH(I,J)+dH(1,J))/mmid
          pgfx(i,j) = pgfx(i,j)*bydx
        endif
      ENDDO
      if(any(nbyzu(j1p:JMAXpfu,l).gt.0)) then
        Call OPFIL2(pgfx,l,J1P,JMAXpfu)
      endif

      do j=j1p-1,jnp
        bydy = 1d0/dyvo(j)
        if(j == jm-1) bydy = bydy*2./3.
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            mmid = (moh(i,j,l)+moh(i,j+1,l))
            pgfy(i,j) = (ZG(I,J)-ZG(I,J+1))+(P(I,J)-P(I,J+1))*
     &             (dH(I,J)+dH(I,J+1))/mmid
            pgfy(i,j) = pgfy(i,j)*bydy
          enddo
        enddo
      enddo

c
c calculate vorticity at cell corners
c
      do j=j1p-1,jnp  ! merge with pgfy loop
        do n=1,nbyzc(j,l)
          do i=i1yzc(n,j,l),min(im-1,i2yzc(n,j,l))
            vort(i,j) = (dxpo(j)*usmooth(i,j)-dxpo(j+1)*usmooth(i,j+1)
     &                  +dyvo(j)*(voh(i+1,j,l)-voh(i,j,l)))/dxyvo(j)
          enddo
        enddo
        i=im
        vort(i,j) = (dxpo(j)*usmooth(i,j)-dxpo(j+1)*usmooth(i,j+1)
     &              +dyvo(j)*(voh(1,j,l)-voh(i,j,l)))/dxyvo(j)
      enddo

c
c update UO, VOD
c
      Do J=J1P,JNP
        bydx = 1d0/dxpo(j)
        corofj = 2.*omega*sinpo(j)
        do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),min(im-1,i2yzu(n,j,l))
            vq = .25*(va(i,j)+va(i+1,j))*(vort(i,j-1)+vort(i,j))
            uo(i,j,l) = uo(i,j,l) + dt*(
     &           pgfx(i,j) +(ke(i,j)-ke(i+1,j))*bydx
     &           +vodh(i,j,l)*corofj + vq)
            pgf4pt=.25*(pgfy(i,j)+pgfy(i+1,j)+pgfy(i,j-1)+pgfy(i+1,j-1))
            vod(i,j,l) = vod(i,j,l) + dt*(pgf4pt - uoh(i,j,l)*corofj)
          enddo
        enddo
        i = im
        if(l <= lmou(i,j)) then
          vq = .25*(va(i,j)+va(1,j))*(vort(i,j-1)+vort(i,j))
          uo(i,j,l) = uo(i,j,l) + dt*(
     &         pgfx(i,j) +(ke(i,j)-ke(1,j))*bydx
     &         +vodh(i,j,l)*corofj + vq)
          pgf4pt = .25*(pgfy(i,j)+pgfy(1,j)+pgfy(i,j-1)+pgfy(1,j-1))
          vod(i,j,l) = vod(i,j,l) + dt*(pgf4pt - uoh(i,j,l)*corofj)
        endif
      enddo

c
c update VO, UOD
c

      do j=j1h,jnp
        bydy = 1d0/dyvo(j)
        if(j == jm-1) bydy = bydy*2./3.
        corofj = 2.*omega*sinvo(j)
        mvfac = .5*dxvo(j)
        pgfac = .25
        if(j == jm-1) pgfac = .5
        i = 1
        if(l <= lmov(i,j)) then
          mmid = (moh(i,j,l)+moh(i,j+1,l))
          mv(i,j) = mvfac*voh(i,j,l)*mmid
          smv(i,j,l) = smv(i,j,l) + mv(i,j)*xeven
          uq = .25*(ua(i,j)+ua(i,j+1))*(vort(im,j)+vort(i,j))
          vo(i,j,l) = vo(i,j,l) + dt*(
     &         pgfy(i,j) +(ke(i,j)-ke(i,j+1))*bydy
     &         -uodh(i,j,l)*corofj -uq)
          pgf4pt = pgfac*
     &         (pgfx(im,j)+pgfx(i,j)+pgfx(im,j+1)+pgfx(i,j+1))
          uod(i,j,l) = uod(i,j,l) + dt*(pgf4pt + voh(i,j,l)*corofj)
        endif
        do n=1,nbyzv(j,l)
          do i=max(2,i1yzv(n,j,l)),i2yzv(n,j,l)
            mmid = (moh(i,j,l)+moh(i,j+1,l))
            mv(i,j) = mvfac*voh(i,j,l)*mmid
            smv(i,j,l) = smv(i,j,l) + mv(i,j)*xeven
            uq = .25*(ua(i,j)+ua(i,j+1))*(vort(i-1,j)+vort(i,j))
            vo(i,j,l) = vo(i,j,l) + dt*(
     &           pgfy(i,j) +(ke(i,j)-ke(i,j+1))*bydy
     &           -uodh(i,j,l)*corofj -uq)
            pgf4pt = pgfac*
     &           (pgfx(i-1,j)+pgfx(i,j)+pgfx(i-1,j+1)+pgfx(i,j+1))
            uod(i,j,l) = uod(i,j,l) + dt*(pgf4pt + voh(i,j,l)*corofj)
          enddo
        enddo
      enddo

c
c update polar velocities
c
      call polevel(uo(1,j1a,l),vo(1,j1a,l),l)


c
c update MO
c
      Do J=J1P,JNP
        convfac = dt/dxypo(j)
        i = 1
        if(l <= lmom(i,j)) then
          convij = convfac*(MU(IM ,J)-MU(I,J) + MV(I,J-1)-MV(I,J))
          mo(i,j,l) = mo(i,j,l) + convij
          opbot(i,j) = opbot(i,j) + convij*grav
          if(mo(i,j,l)*dxypo(j) < 0.75d0*mmi(i,j,l)) then
            write(6,*) 'small mo',i,j,l,
     &           mo(i,j,l)*dxypo(j)/mmi(i,j,l),lmom(i,j)
            smallmo_loc = 1
          endif
        endif
        do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            convij = convfac*(MU(I-1,J)-MU(I,J) +MV(I,J-1)-MV(I,J))
            mo(i,j,l) = mo(i,j,l) + convij
            opbot(i,j) = opbot(i,j) + convij*grav
            if(mo(i,j,l)*dxypo(j) < 0.75d0*mmi(i,j,l)) then
              write(6,*) 'small mo',i,j,l,
     &             mo(i,j,l)*dxypo(j)/mmi(i,j,l),lmom(i,j)
              smallmo_loc = 1
            endif
          enddo
        enddo
      enddo
      If (QNP)  then
        j = jm
        convij = dt*sum(mv(:,jm-1))/(im*dxypo(jm))
        mo(1,j,l) = mo(1,j,l) + convij
        opbot(1,j) = opbot(1,j) + convij*grav
        do i=2,im
          mo(i,j,l) = mo(1,j,l)
        enddo
      endif

      enddo ! end loop over layers

      CALL GLOBALMAX(grid, smallmo_loc, smallmo)
      if(smallmo .gt. 0) call stop_model('small mo',255)

c emergency vertical regrid
c      Call OFLUXV(opbot,mo,qeven)

c fill the halos of just-updated quantities
      call halo_update(grid,mo)
      call halo_update(grid,opbot)
      call halo_update(grid,uo)
      call halo_update(grid,vo)

      Return
      End Subroutine ODHORZ

      subroutine polevel(u,v,l)
      Use OCEAN, Only: IM,JM,COSIC,SINIC, SINU,COSU
      USE OCEAN, only :
     &     nbyzm,nbyzu,nbyzv,
     &     i1yzm,i2yzm, i1yzu,i2yzu, i1yzv,i2yzv
      USE OCEANR_DIM, only : grid=>ogrid
      use DOMAIN_DECOMP_1D, only: hasNorthPole
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) :: u,v
      integer, intent(in) :: l
      integer :: i,j,n
      real*8 :: unp,vnp
c
c calculate U and V at the north pole
c
      if(hasNorthPole(grid)) then
        unp = 0.
        vnp = 0.
        j = jm-1
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            unp = unp - sinic(i)*v(i,j)
            vnp = vnp + cosic(i)*v(i,j)
          enddo
        enddo
        unp = unp*2/im
        vnp = vnp*2/im
c        vonp(l) = vnp
        do i=1,im
          u(i,jm) = unp*cosu(i)  + vnp*sinu(i)
          v(i,jm) = vnp*cosic(i) - unp*sinic(i)
        enddo
      endif
      return
      end subroutine polevel

      subroutine ocnstate_derived
!@sum ocnstate_derived define intensive-units thermodynamic quantities
!@+   from extensive-units state variables
      use ocean, only : g3d,t3d,s3d,p3d,rho=>r3d,vbar=>v3d
      use constant, only: grav
      use ocean, only : G0M,GZM=>GZMO, S0M,SZM=>SZMO, OPRESS, FOCEAN, MO
      use ocean, only : im,jm,lmo,dxypo
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds,halo_update_column
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      implicit none

c
      Real*8, External :: VOLGSP,temgsp

      integer :: i,j,l,n

      INTEGER :: J_0, J_1
      LOGICAL :: HAVE_NORTH_POLE,HAVE_SOUTH_POLE

      Real*8, Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Real*8 :: gup,gdn,sup,sdn,pm,vup,vdn,bym

      real*8, dimension(0:lmo) :: pe

c**** Extract domain decomposition info

      call getdomainbounds(grid, j_strt = j_0, j_stop = j_1,
     &               have_north_pole = have_north_pole,
     &               have_south_pole = have_south_pole)


      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        pe(0) = opress(i,j)
        do l=1,lmm(i,j)
          bym = 1d0/(mo(i,j,l)*dxypo(j))
          pe(l) = pe(l-1) + mo(i,j,l)*grav
          g3d(l,i,j) = g0m(i,j,l)*bym
          s3d(l,i,j) = s0m(i,j,l)*bym
          p3d(l,i,j) = .5*(pe(l)+pe(l-1))
          pm = p3d(l,i,j)
C**** In-situ temperature
          t3d(l,i,j) = temgsp(g3d(l,i,j),s3d(l,i,j),pm)
C**** Specific volume ref to mid-point pressure
          gup = (g0m(i,j,l)-2*z12eh*gzm(i,j,l))*bym
          gdn = (g0m(i,j,l)+2*z12eh*gzm(i,j,l))*bym
          sup = (s0m(i,j,l)-2*z12eh*szm(i,j,l))*bym
          sdn = (s0m(i,j,l)+2*z12eh*szm(i,j,l))*bym
          sup = max(0d0,sup)
          sdn = max(0d0,sdn)
          vup = volgsp(gup,sup,pm)
          vdn = volgsp(gdn,sdn,pm)
          vbar(l,i,j) = (vup + vdn)*.5
C**** In-situ density = 1/specific volume
          rho(l,i,j)  = 1d0/vbar(l,i,j)
        enddo
      enddo
      enddo
      enddo

C**** Copy to all longitudes at poles

      if(have_north_pole) then
        do l=1,lmm(1,jm)
          rho(l,2:im,jm) = rho(l,1,jm)
          vbar(l,2:im,jm) = vbar(l,1,jm)
          g3d(l,2:im,jm) = g3d(l,1,jm)
          s3d(l,2:im,jm) = s3d(l,1,jm)
          p3d(l,2:im,jm) = p3d(l,1,jm)
        enddo
      endif
      if(have_south_pole) then
        do l=1,lmm(1,1)
          rho(l,2:im,1) = rho(l,1,1)
          vbar(l,2:im,1) = vbar(l,1,1)
          g3d(l,2:im,1) = g3d(l,1,1)
          s3d(l,2:im,1) = s3d(l,1,1)
          p3d(l,2:im,1) = p3d(l,1,1)
        enddo
      endif

      call halo_update_column(grid,vbar)
      call halo_update_column(grid,rho)
      call halo_update_column(grid,g3d)
      call halo_update_column(grid,s3d)
      call halo_update_column(grid,p3d)

      end subroutine ocnstate_derived

      Subroutine ODHORZ0
C****
C**** OPGF0 calculates GUP, GDN, SUP and SDN inside each ocean cell,
C**** which will be kept constant during the dynamics of momentum
C****
      use constant, only : grav
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, MO,UO,VO,OPBOT,OPRESS,
     *                 G0M,GZM=>GZMO, S0M,SZM=>SZMO, DXYPO
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      Use OCEAN_DYN, Only: MMI, GUP,GDN, SUP,SDN, VBAR, dZGdP
      use ocean_dyn, only : dh3d=>dh ! for interface with other ocean routines
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, NORTH
      USE OCEANR_DIM, only : grid=>ogrid

      Implicit None
      Real*8,Parameter :: z12eH=.28867513d0  !  z12eH = 1/SQRT(12)
      Integer*4 I,J,L, J1,J1H,JN,J1A,JNH, N
      Logical :: QNP
      Real*8,External   :: VOLGSP
      Real*8,Dimension(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LMO) ::
     *          P
      Real*8 :: PUP,PDN, VUP,VDN, smean
C****
C**** Extract domain decomposition band parameters
C****                          Band1  Band2  BandM
      J1  = GRID%J_STRT     !    1      5     JM-3   Band minimum
      J1H = Max(J1-1,1)     !    1      4     JM-4   Halo minimum
      JN  = GRID%J_STOP     !    4      8     JM     Band maximum
      JNH = Min(JN+1,JM)    !    5      9     JM     Halo maximum
      QNP = JN==JM          !    F      F      T
      j1a = grid%j_strt_halo

C****
      Call HALO_UPDATE (GRID, G0M)
      Call HALO_UPDATE (GRID, GZM)
      Call HALO_UPDATE (GRID, S0M)
      Call HALO_UPDATE (GRID, SZM)

C**** Calculate pressure by integrating from the top down
      Call HALO_UPDATE (GRID,OPRESS)
      l=1
      do j=j1h,jnh
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            opbot(i,j) = OPRESS(I,J)
          enddo
        enddo
      enddo
      do l=1,lmo
      do j=j1h,jnh
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            P(I,J,L) = opbot(i,j)  + MO(I,J,L)*GRAV*.5
            opbot(i,j) = opbot(i,j) + MO(I,J,L)*GRAV
          enddo
        enddo
      enddo
      enddo

C****
      do l=lmo,1,-1
        Do J=J1H,JNH
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              MMI(I,J,L) = MO(I,J,L)*DXYPO(J)
              GUP(I,J,L) = (G0M(I,J,L) - 2*z12eH*GZM(I,J,L))/MMI(I,J,L)
              GDN(I,J,L) = (G0M(I,J,L) + 2*z12eH*GZM(I,J,L))/MMI(I,J,L)
              SUP(I,J,L) = (S0M(I,J,L) - 2*z12eH*SZM(I,J,L))/MMI(I,J,L)
              SDN(I,J,L) = (S0M(I,J,L) + 2*z12eH*SZM(I,J,L))/MMI(I,J,L)
              smean = s0m(i,j,l)/mmi(i,j,l)
              SUP(I,J,L) = MAX(SUP(I,J,L), .5d0*smean)
              SDN(I,J,L) = MAX(SDN(I,J,L), .5d0*smean)
              PUP = P(I,J,L) - MO(I,J,L)*GRAV*z12eH
              PDN = P(I,J,L) + MO(I,J,L)*GRAV*z12eH
              VUP = VOLGSP (GUP(I,J,L),SUP(I,J,L),PUP)
              VDN = VOLGSP (GDN(I,J,L),SDN(I,J,L),PDN)
              dZGdP(I,J,L) = VUP*(.5-z12eH) + VDN*(.5+z12eH)
              VBAR(I,J,L) = (VUP + VDN)*.5
              DH3D(I,J,L) = MO(I,J,L)*VBAR(I,J,L)
            enddo
          enddo
        enddo

C**** Copy to all longitudes at north pole
        if(QNP) then
          if(l <= lmom(1,jm)) then
            DH3D(2:IM,JM,L) = DH3D(1,JM,L)
            VBAR(2:IM,JM,L) = VBAR(1,JM,L)
            dZGdP(2:IM,JM,L) = dZGdP(1,JM,L)
            MO(2:IM,JM,L) = MO(1,JM,L)
          endif
        endif

c initialize polar velocities
        call polevel(uo(1,j1a,l),vo(1,j1a,l),l)

      enddo

      Return
      End Subroutine ODHORZ0

      SUBROUTINE OADVT2 (MA,RM,RX,RY,RZ,DT,QLIMIT, OIJL)
!@sum  OADVT advects tracers using the linear upstream scheme.
C****
C**** Input:  MB (kg) = mass before advection
C****          DT (s) = time step
C****       MU (kg/s) = west to east mass flux
C****       MV (kg/s) = south to north mass flux
C****       MW (kg/s) = downward vertical mass flux
C****          QLIMIT = whether slope limitations should be used
C**** Output: RM (kg) = tracer mass
C****   RX,RY,RZ (kg) = first moments of tracer mass
C****       OIJL (kg) = diagnostic accumulation of tracer mass flux
C****
      USE OCEAN, only : im,jm,lmo
      USE OCEAN_DYN, only : mb=>mmi,smu,smv,smw

!      use domain_decomp_1d, only : grid, get
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),     DIMENSION
     *     (IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MA,RM,RX,RY,RZ
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO,3) :: OIJL
      INTEGER I,J,L, J_0H
      LOGICAL, INTENT(IN) :: QLIMIT
      REAL*8, INTENT(IN) :: DT

      logical :: HAVE_NORTH_POLE  ! ,HAVE_SOUTH_POLE
         call getDomainBounds(grid, J_STRT_HALO=J_0H)
         call getDomainBounds(grid, HAVE_NORTH_POLE=HAVE_NORTH_POLE)
C        call getDomainBounds(grid, HAVE_SOUTH_POLE=HAVE_SOUTH_POLE)

C****
C**** Load mass after advection from mass before advection
C****
      MA = MB
C****
C**** Advect the tracer using the Slopes Scheme
C****
      CALL OADVTX2(RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))
      CALL OADVTY2(RM,RX,RY,RZ,MA,SMV,     DT,QLIMIT,OIJL(1,J_0H,1,2))
      CALL OADVTZ2(RM,RX,RY,RZ,MA,SMW,     DT,QLIMIT,OIJL(1,J_0H,1,3))
      CALL OADVTX2(RM,RX,RY,RZ,MA,SMU,5d-1*DT,QLIMIT,OIJL(1,J_0H,1,1))

C**** Fill in values at the poles
      if (HAVE_NORTH_POLE) then
      DO 20 L=1,LMO
      DO 20 I=1,IM
      RM(I,JM,L) = RM(1,JM,L)
   20 RZ(I,JM,L) = RZ(1,JM,L)
      end if
C     if (HAVE_SOUTH_POLE) then
C     DO 30 L=1,LMO
C     DO 30 I=1,IM
C     RM(I, 1,L) = RM(IM,1,L)
C  30 RZ(I, 1,L) = RZ(IM,1,L)
C     end if

      RETURN
      END SUBROUTINE OADVT2

      Subroutine OADVTX2(RM,RX,RY,RZ,MM,MU,DT,QLIMIT,OIJL)
C****
!@sum   OADVTX advects tracer in X direction via Linear Upstream Scheme
C****
C**** If QLIMIT is true, gradients are limited to prevent mean tracer
C**** from becoming negative; Abs[A(I)] must not exceed 1.
C****
C**** Input: DT (s) = time step
C****     MU (kg/s) = west to east mass flux
C****        QLIMIT = whether slope limitations should be used
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                   MM (kg) = ocean mass
C****
      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOU=>LMU
      USE OCEAN, only : nbyzu,i1yzu,i2yzu, nbyzm,i1yzm,i2yzm
      Use OCEANR_DIM, Only: oGRID
      Implicit None
C**** Interface variables
      Logical,Intent(In) :: QLIMIT
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  MU
      Real*8,Intent(InOut),
     *  Dimension(IM, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ, MM,OIJL
C**** Local variables
      Integer :: I,J,L, J1P,JNP, N,NC,NCOURANT
      real*8, dimension(im) :: mudt
      Real*8 :: AM,A,FM,FX,FY,FZ,MMnew, zCOURANT,mcheck,courmax,rxlimit
      Real*8 :: AMim1,FMim1,FXim1,FYim1,FZim1,
     &          AM_im,FM_im,FX_im,FY_im,FZ_im

C**** Extract domain decomposition band parameters
      J1P = Max (oGRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (oGRID%J_STOP, JM-1)  !  Exclude north pole

      if(qlimit) then
        rxlimit = 1d0
      else
        rxlimit = 0d0
      endif

C****
C**** Loop over layers and latitudes
C****
      Do L=1,LMO
      Do J=J1P,JNP
        if(nbyzu(j,l) == 0) cycle


c
c adjust mass fluxes so that courant numbers are never greater than 1
c
        courmax = 0.
        mudt(1:2) = mu(1:2,j,l)*dt
        mudt(im) = mu(im,j,l)*dt
        i=1
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(im)-mudt(i))
          else
            mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        do n=1,nbyzu(j,l)
          i = i1yzu(n,j,l)
          if(i.gt.1) mudt(i-1) = 0. ! zero mudt at west boundary
          mudt(i) = mu(i,j,l)*dt
          do i=max(2,i1yzu(n,j,l)),min(i2yzu(n,j,l),im-1)
            mudt(i+1) = mu(i+1,j,l)*dt
            if(mudt(i).ge.0.) then
              mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
            else
              mcheck = -mm(i+1,j,l)-min(0d0,mudt(i)-mudt(i+1))
            endif
            courmax = max(courmax,mudt(i)/mcheck)
          enddo
        enddo
        i = im
        if(l <= lmou(i,j)) then
          if(mudt(i).ge.0.) then
            mcheck = mm(i,j,l)+min(0d0,mudt(i-1)-mudt(i))
          else
            mcheck = -mm(1,j,l)-min(0d0,mudt(i)-mudt(1))
          endif
          courmax = max(courmax,mudt(i)/mcheck)
        endif
        if(courmax.gt.1.) then
          ncourant = 1+int(courmax)
          write(6,*) 'ncourant ',j,l,ncourant
          zcourant = 1d0/ncourant
          do n=1,nbyzu(j,l)
            do i=i1yzu(n,j,l),i2yzu(n,j,l)
              mudt(i) = mudt(i)*zcourant
            enddo
          enddo
        else
          ncourant = 1
        endif

        do nc=1,ncourant

c
c calculate fluxes at the dateline
c
          i = im
          if(l <= lmou(i,j)) then
            AM = mudt(I)
            if(am.ge.0.) then   ! mass flux is positive or zero
              A  = AM / MM(I,J,L)
              rx(i,j,l) = rx(i,j,l)-rxlimit*sign(min(0d0,
     &             rm(i,j,l)-abs(rx(i,j,l))),rx(i,j,l))
              FM = A * (RM(I,J,L) + (1-A)*RX(I,J,L))
              FX = AM * (A*A*RX(I,J,L) - 3*FM)
              FY = A*RY(I,J,L)
              FZ = A*RZ(I,J,L)
            else                ! mass flux is negative
              A  = AM / MM(1,J,L)
              rx(1,j,l) = rx(1,j,l)-rxlimit*sign(min(0d0,
     &             rm(1,j,l)-abs(rx(1,j,l))),rx(1,j,l))
              FM = A * (RM(1,J,L) - (1+A)*RX(1,J,L))
              FX = AM * (A*A*RX(1,J,L) - 3*FM)
              FY = A*RY(1,J,L)
              FZ = A*RZ(1,J,L)
            endif
          else
            am = 0.
            fm = 0.
            fx = 0.
            fy = 0.
            fz = 0.
          endif
          amim1 = am
          fmim1 = fm
          fxim1 = fx
          fyim1 = fy
          fzim1 = fz
          am_im = am
          fm_im = fm
          fx_im = fx
          fy_im = fy
          fz_im = fz

c
c loop over basins
c
          do n=1,nbyzm(j,l)
            if(i1yzm(n,j,l)>1 .and. i1yzm(n,j,l)==i2yzm(n,j,l)) cycle
            do i=i1yzm(n,j,l),min(i2yzm(n,j,l),im-1)
              AM = mudt(I)
              if(am.ge.0.) then ! mass flux is positive or zero
                A  = AM / MM(I,J,L)
                rx(i,j,l) = rx(i,j,l)-rxlimit*sign(min(0d0,
     &               rm(i,j,l)-abs(rx(i,j,l))),rx(i,j,l))
                FM = A * (RM(I,J,L) + (1-A)*RX(I,J,L))
                FX = AM * (A*A*RX(I,J,L) - 3*FM)
                FY = A*RY(I,J,L)
                FZ = A*RZ(I,J,L)
              else              ! mass flux is negative
                A  = AM / MM(I+1,J,L)
                rx(i+1,j,l) = rx(i+1,j,l)-rxlimit*sign(min(0d0,
     &               rm(i+1,j,l)-abs(rx(i+1,j,l))),rx(i+1,j,l))
                FM = A * (RM(I+1,J,L) - (1+A)*RX(I+1,J,L))
                FX = AM * (A*A*RX(I+1,J,L) - 3*FM)
                FY = A*RY(I+1,J,L)
                FZ = A*RZ(I+1,J,L)
              endif
              MMnew = MM(I,J,L) + (AMim1-AM)
              RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
              RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &             3d0*((AMim1+AM)*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))
     &             / MMnew
              RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
              RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
              MM(I,J,L) = MMnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
              amim1 = am
              fmim1 = fm
              fxim1 = fx
              fyim1 = fy
              fzim1 = fz
            enddo
          enddo

c
c update at dateline
c
          i = im
          if(l <= lmom(i,j)) then
            am = am_im
            fm = fm_im
            fx = fx_im
            fy = fy_im
            fz = fz_im
            MMnew = MM(I,J,L) + (AMim1-AM)
            RM(I,J,L) = RM(I,J,L) +  (FMim1-FM)
            RX(I,J,L) = (RX(I,J,L)*MM(I,J,L) + (FXim1-FX) +
     &           3d0*((AMim1+AM)*RM(I,J,L)-MM(I,J,L)*(FMim1+FM)))
     &           / MMnew
            RY(I,J,L) = RY(I,J,L) + (FYim1-FY)
            RZ(I,J,L) = RZ(I,J,L) + (FZim1-FZ)
            MM(I,J,L) = MMnew
            OIJL(I,J,L) = OIJL(I,J,L) + FM
          endif
        enddo ! nc
      enddo ! j
      enddo ! l

      Return
      EndSubroutine OADVTX2

      Subroutine OADVTY2(RM,RX,RY,RZ,MO,MV,DT,QLIMIT,OIJL)
C****
!@sum   OADVTY advects tracer in Y direction via Linear Upstream Scheme
c****
c**** oadvty advects tracers in the south to north direction using the
c**** linear upstream scheme.  if qlimit is true, the moments are
c**** limited to prevent the mean tracer from becoming negative.
c****
c**** input:
c****     mv (kg/s) = north-south mo flux, positive northward
c****        qlimit = whether slope limitations should be used
c****
c**** input/output:
c****     rm     (kg) = tracer mass
c****   rx,ry,rz (kg) = 1st moments of tracer mass
c****     mo     (kg) = ocean mass
c****
C****

      Use OCEAN, Only: IM,JM,LMO, LMOM=>LMM, LMOV=>LMV
      USE OCEAN, only : nbyzm,i1yzm,i2yzm, nbyzv,i1yzv,i2yzv
      Use OCEANR_DIM, Only: grid=>oGRID
      use DOMAIN_DECOMP_1D, only : halo_update, 
     &     hasSouthPole, hasNorthPole
      Implicit None
C**** Interface variables
      Logical,Intent(In) :: QLIMIT
      Real*8,   Intent(In) :: DT
      Real*8,Intent(In),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  MV
      Real*8,Intent(InOut),
     *  Dimension(IM, GRID%J_STRT_HALO:GRID%J_STOP_HALO, LMO) ::
     *  RM,RX,RY,RZ, MO,OIJL
C**** Local variables
      logical :: qnp
      Integer :: I,J,L, J1H,J1P,JNP,JN, N, imin,imax
      Real*8 :: BM,B,FM,FX,FY,FZ,Mnew,rylimit
      real*8, dimension(im) :: BMjm1,FMjm1,FXjm1,FYjm1,FZjm1

C**** Extract domain decomposition band parameters
      J1H = Max (GRID%J_STRT-1, 1)   !  Exclude south pole
      J1P = Max (GRID%J_STRT, 2)     !  Exclude south pole
      JNP = Min (GRID%J_STOP, JM-1)  !  Exclude north pole
      JN = GRID%J_STOP
      QNP = hasNorthPole(grid)

      if(qlimit) then
        rylimit = 1d0
      else
        rylimit = 0d0
      endif

      call halo_update(grid,mo)
      call halo_update(grid,rm)
      call halo_update(grid,rx)
      call halo_update(grid,ry)
      call halo_update(grid,rz)

      do l=1,lmo

        do i=1,im
          bmjm1(i) = 0.
          fmjm1(i) = 0.
          fxjm1(i) = 0.
          fyjm1(i) = 0.
          fzjm1(i) = 0.
        enddo

        if(qnp) then
          j = jm
          do i=2,im
            mo(i,j,l) = mo(1,j,l)
            rm(i,j,l) = rm(1,j,l)
            rz(i,j,l) = rz(1,j,l)
          enddo
          do i=1,im
            rx(i,j,l) = 0.
            ry(i,j,l) = 0.
          enddo
        endif

        j = j1h
        do n=1,nbyzv(j,l)
          do i=i1yzv(n,j,l),i2yzv(n,j,l)
            bm = mv(i,j,l)*dt
            if(bm.ge.0.) then   ! mass flux is positive or zero
              B  = BM / MO(I,J,L)
              ry(i,j,l) = ry(i,j,l)-rylimit*sign(min(0d0,
     &             rm(i,j,l)-abs(ry(i,j,l))),ry(i,j,l))
              FM = B * (RM(I,J,L) + (1-B)*RY(I,J,L))
              FY = BM * (B*B*RY(I,J,L) - 3*FM)
              FX = B*RX(I,J,L)
              FZ = B*RZ(I,J,L)
            else                ! mass flux is negative
              B  = BM / MO(i,j+1,L)
              ry(i,j+1,l) = ry(i,j+1,l)-rylimit*sign(min(0d0,
     &             rm(i,j+1,l)-abs(ry(i,j+1,l))),ry(i,j+1,l))
              FM = B * (RM(i,j+1,L) - (1+B)*RY(i,j+1,L))
              FY = BM * (B*B*RY(i,j+1,L) - 3*FM)
              FX = B*RX(i,j+1,L)
              FZ = B*RZ(i,j+1,L)
            endif
            bmjm1(i) = bm
            fmjm1(i) = fm
            fxjm1(i) = fx
            fyjm1(i) = fy
            fzjm1(i) = fz
          enddo
        enddo

        do j=j1p,jn
          do n=1,nbyzm(j,l)
            imin=i1yzm(n,j,l)
            imax=i2yzm(n,j,l)
            if(j == jm) imax=im
            do i=imin,imax
              bm = mv(i,j,l)*dt
              if(bm.ge.0.) then ! mass flux is positive or zero
                B  = BM / MO(I,J,L)
                ry(i,j,l) = ry(i,j,l)-rylimit*sign(min(0d0,
     &               rm(i,j,l)-abs(ry(i,j,l))),ry(i,j,l))
                FM = B * (RM(I,J,L) + (1-B)*RY(I,J,L))
                FY = BM * (B*B*RY(I,J,L) - 3*FM)
                FX = B*RX(I,J,L)
                FZ = B*RZ(I,J,L)
              else              ! mass flux is negative
                B  = BM / MO(i,j+1,L)
                ry(i,j+1,l) = ry(i,j+1,l)-rylimit*sign(min(0d0,
     &               rm(i,j+1,l)-abs(ry(i,j+1,l))),ry(i,j+1,l))
                FM = B * (RM(i,j+1,L) - (1+B)*RY(i,j+1,L))
                FY = BM * (B*B*RY(i,j+1,L) - 3*FM)
                FX = B*RX(i,j+1,L)
                FZ = B*RZ(i,j+1,L)
              endif
              Mnew = MO(I,J,L) + (BMjm1(i)-BM)
              RM(I,J,L) = RM(I,J,L) +  (FMjm1(i)-FM)
              RY(I,J,L) = (RY(I,J,L)*MO(I,J,L) + (FYjm1(i)-FY) + 3d0*
     &             ((BMjm1(i)+BM)*RM(I,J,L)-MO(I,J,L)*(FMjm1(i)+FM)))
     &             / Mnew
              RX(I,J,L) = RX(I,J,L) + (FXjm1(i)-FX)
              RZ(I,J,L) = RZ(I,J,L) + (FZjm1(i)-FZ)
              MO(I,J,L) = Mnew
              OIJL(I,J,L) = OIJL(I,J,L) + FM
              bmjm1(i) = bm
              fmjm1(i) = fm
              fxjm1(i) = fx
              fyjm1(i) = fy
              fzjm1(i) = fz
            enddo
          enddo
        enddo

c
c average the pole
c
        if(qnp) then
          if(l<=lmom(1,jm)) then
            j = jm
            mo(:,j,l) = sum(mo(:,j,l))/im
            rm(:,j,l) = sum(rm(:,j,l))/im
            rz(:,j,l) = sum(rz(:,j,l))/im
            do i=1,im
              rx(i,j,l) = 0.
              ry(i,j,l) = 0.
            enddo
          endif
        endif

      enddo ! l

      return
      end subroutine oadvty2

      SUBROUTINE OADVTZ2(RM,RX,RY,RZ,MO,MW,DT,QLIMIT,OIJL)
C****
C**** OADVTZ advects tracers in the vertical direction using the
C**** linear upstream scheme.  If QLIMIT is true, the gradients are
C**** limited to prevent the mean tracer from becoming negative.
C****
C****
C**** Input: DT (s) = time step
C****     MW (kg/s) = downward vertical mass flux
C**** Input and Output: RM (kg) = tracer mass
C****             RX,RY,RZ (kg) = first moments of tracer mass
C****                    M (kg) = ocean mass
C****
      USE OCEAN, only : im,jm,lmo,lmm
      USE OCEAN, only : nbyzm,i1yzm,i2yzm
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid

      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     *  RM,RX,RY,RZ, OIJL, MO
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MW
      REAL*8, INTENT(IN) :: DT
      LOGICAL, INTENT(IN) :: QLIMIT
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     CMUP,FMUP,FXUP,FYUP,FZUP,FMUP_ctr
      REAL*8 :: CM,C,FM,FX,FY,FZ,MNEW,rzlim,no_rzlim
      real*8 :: fm_ctr,r_edge,wtdn,rdn,rup,rm_lus
c Applied when qlimit is true, edgmax is the maximum allowed ratio of
c r_edge to gridbox-mean r.  edgmax = 2 can be used if there is never
c flow out of the bottom and top edges simultaneously.
      real*8, parameter :: edgmax=1.5
      INTEGER :: I,J,L,N,LDN

      INTEGER :: J_0,J_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

      if(qlimit) then
        rzlim = 1d0
      else
        rzlim = 0d0
      endif
      no_rzlim = 1d0-rzlim

      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            cmup(i,j) = 0.
            fmup(i,j) = 0.
            fxup(i,j) = 0.
            fyup(i,j) = 0.
            fzup(i,j) = 0.
            fmup_ctr(i,j) = 0.
          enddo
        enddo
      enddo

      do l=1,lmo
        do j=j_0,j_1
          do n=1,nbyzm(j,l)
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              CM = DT*MW(I,J,L)
              ldn = min(l+1,lmm(i,j))
              wtdn = mo(i,j,l)/(mo(i,j,l)+mo(i,j,ldn)) ! use dzo instead?
              rdn = rm(i,j,ldn)/mo(i,j,ldn)
              rup = rm(i,j,l  )/mo(i,j,l  )
              R_edge = wtdn*rdn+(1.-wtdn)*rup
              if(cm.ge.0.) then ! mass flux is downward or zero
                C  = CM/MO(I,J,L)
                IF(C.GT.1d0)  WRITE (6,*) 'C>1:',I,J,L,C,MO(I,J,L)
                rz(i,j,l) = rz(i,j,l)-rzlim*sign(min(0d0,
     &               rm(i,j,l)-abs(rz(i,j,l))),rz(i,j,l))
                FM = C*(RM(I,J,L)+(1d0-C)*RZ(I,J,L))
                FX = C*RX(I,J,L)
                FY = C*RY(I,J,L)
                FZ = CM*(C*C*RZ(I,J,L)-3d0*FM)
                R_edge = R_edge*no_rzlim+rzlim*min(R_edge,edgmax*rup)
                fm_ctr = cm*(c*rup + (1.-c)*r_edge)
              else              ! mass flux is upward
                C  = CM/MO(I,J,L+1)
                IF(C.LT.-1d0)  WRITE (6,*) 'C<-1:',I,J,L,C,MO(I,J,L+1)
                rz(i,j,l+1) = rz(i,j,l+1)-rzlim*sign(min(0d0,
     &               rm(i,j,l+1)-abs(rz(i,j,l+1))),rz(i,j,l+1))
                FM = C*(RM(I,J,L+1)-(1d0+C)*RZ(I,J,L+1))
                FX = C*RX(I,J,L+1)
                FY = C*RY(I,J,L+1)
                FZ = CM*(C*C*RZ(I,J,L+1)-3d0*FM)
                R_edge = R_edge*no_rzlim+rzlim*min(R_edge,edgmax*rdn)
                fm_ctr = cm*(-c*rdn + (1.+c)*r_edge)
              endif
#ifdef LUS_VERT_ADV
              fm_ctr = fm
#endif
              RM_lus = RM(I,J,L) + (FMUP(I,J)-FM)
              RM(I,J,L) = RM(I,J,L) + (FMUP_ctr(I,J)-FM_ctr)
              mnew = MO(I,J,L) + CMUP(I,J)-CM
              RZ(I,J,L) = (RZ(I,J,L)*MO(I,J,L) + (FZUP(I,J)-FZ) +3d0*
     &             ((CMUP(I,J)+CM)*RM_lus-MO(I,J,L)*(FMUP(I,J)+FM)))
     &             / mnew
              MO(I,J,L) = mnew
              cmup(i,j) = cm
              fmup(i,j) = fm
              fmup_ctr(i,j) = fm_ctr
              fzup(i,j) = fz
              RX(I,J,L) = RX(I,J,L) + (FXUP(I,J)-FX)
              fxup(i,j) = fx
              RY(I,J,L) = RY(I,J,L) + (FYUP(I,J)-FY)
              fyup(i,j) = fy
              OIJL(I,J,L) = OIJL(I,J,L) + FM_ctr
            enddo
          enddo
        enddo
      enddo

      return
      END SUBROUTINE OADVTZ2

      SUBROUTINE OADVUZ(R,M,MW,DT,jmin,jmax,nbyz,i1yz,i2yz)
c simplest upstream scheme, for vertical direction
      USE OCEAN, only : im,jm,lmo
      USE OCEAN, only :  nbyzmax
      use domain_decomp_1d, only : getDomainBounds
      USE OCEANR_DIM, only : grid=>ogrid
      IMPLICIT NONE
      REAL*8, INTENT(INOUT),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: R,M
      REAL*8, INTENT(IN),
     *  DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: MW
      REAL*8, INTENT(IN) :: DT
      INTEGER, INTENT(IN) :: jmin,jmax
      INTEGER, INTENT(IN),
     &     DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) :: nbyz
      INTEGER, INTENT(IN),
     &     DIMENSION(nbyzmax,grid%J_STRT_HALO:grid%J_STOP_HALO,LMO) ::
     &     i1yz,i2yz
c
      REAL*8, DIMENSION(IM,grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     CMUP,FMUP
      REAL*8 :: CM,FM,MNEW
      INTEGER :: I,J,L,N
      do j=jmin,jmax
        do n=1,nbyz(j,1)
          do i=i1yz(n,j,1),i2yz(n,j,1)
            cmup(i,j) = 0.
            fmup(i,j) = 0.
          enddo
        enddo
      enddo
      do l=1,lmo
        do j=jmin,jmax
          do n=1,nbyz(j,l)
            do i=i1yz(n,j,l),i2yz(n,j,l)
              CM = DT*MW(I,J,L)
              if(cm.ge.0.) then ! mass flux is downward or zero
                FM = CM*R(I,J,L)
              else              ! mass flux is upward
                FM = CM*R(I,J,L+1)
              endif
              mnew = M(I,J,L) + CMUP(I,J)-CM
              R(I,J,L) = (R(I,J,L)*M(I,J,L) + (FMUP(I,J)-FM))/mnew
              cmup(i,j) = cm
              fmup(i,j) = fm
              M(I,J,L) = mnew
            enddo
          enddo
        enddo
      enddo
      return
      END SUBROUTINE OADVUZ

      SUBROUTINE OSTRES2
!@sum OSTRES applies the atmospheric surface stress over open ocean
!@sum and the sea ice stress to the layer 1 ocean velocities
!@auth Gary Russell

      USE OCEAN, only : IMO=>IM,JMO=>JM
     *     , IVNP, UO,VO, UOD,VOD, MO,DXYSO,DXYNO,DXYVO
     *     , LMU,LMV, COSIC,SINIC

      USE DOMAIN_DECOMP_1D, only : getDomainBounds, halo_update, north,
     *     south

      USE OCEANR_DIM, only : ogrid

      USE OFLUXES, only : oDMUA,oDMVA, oDMUI,oDMVI

      IMPLICIT NONE
      INTEGER I,J,IM1,IP1

C****
C**** All stress now defined for whole box, not just ocn or ice fraction
C**** FLUXCB  DMUA(1)  U momentum downward into open ocean (kg/m*s)
C****         DMVA(1)  V momentum downward into open ocean (kg/m*s)
C****         DMUA(2,JM,1)  polar atmo. mass slowed to zero (kg/m**2)
C****         DMUI     U momentum downward from sea ice (kg/m*s)
C****         DMVI     V momentum downward from sea ice (kg/m*s)

      integer :: J_0, J_1, J_0S, J_1S  ; logical :: have_north_pole

      call getDomainBounds(ogrid, J_STRT=J_0, J_STOP=J_1,
     *                 J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     *                 have_north_pole=have_north_pole)

C****
C**** Surface stress is applied to U component
C****
      DO J=J_0S,J_1S
      I=IMO
      DO IP1=1,IMO
        IF(LMU(I,J).gt.0.)  UO(I,J,1) = UO(I,J,1) +
     *       (oDMUA(I,J) + oDMUA(IP1,J) + 2d0*oDMUI(I,J)) /
     *       (  MO(I,J,1) +   MO(IP1,J,1))
        I=IP1
      END DO
      END DO
      if (have_north_pole) then
        UO(IMO ,JMO,1) = UO(IMO ,JMO,1) + oDMUA(1,JMO)/MO(1,JMO,1)
        UO(IVNP,JMO,1) = UO(IVNP,JMO,1) + oDMVA(1,JMO)/MO(1,JMO,1)
      end if
      call halo_update(ogrid, odmvi, from=south)
      DO J=J_0S,J_1S
        I=IMO
        DO IP1=1,IMO
          IF(LMU(I,J).gt.0.) THEN
            VOD(I,J,1) = VOD(I,J,1) + (
     *           oDMVA(I,J)+oDMVA(IP1,J)
     *        +.5*(oDMVI(I,J-1)+oDMVI(IP1,J-1)+oDMVI(I,J)+oDMVI(IP1,J))
     *           )/(MO(I,J,1)+MO(IP1,J,1))
          ENDIF
          I=IP1
        END DO
      END DO
C****
C**** Surface stress is applied to V component
C****
      call halo_update(ogrid, odmva, from=north)
      call halo_update(ogrid,    mo, from=north)
      DO J=J_0S,min(J_1S,JMO-2)
      DO I=1,IMO
        IF(LMV(I,J).GT.0.)  VO(I,J,1) = VO(I,J,1) +
     *       (oDMVA(I,J)*DXYNO(J) + oDMVA(I,J+1)*DXYSO(J+1)
     *      + oDMVI(I,J)*DXYVO(J))  !!  2d0*oDMVI(I,J)*DXYVO(J) - error
     * / (MO(I,J,1)*DXYNO(J) + MO(I,J+1,1)*DXYSO(J+1))
      END DO
      END DO
C**** Surface stress is applied to V component at the North Pole
      if (have_north_pole) then
      DO I=1,IMO
        VO(I,JMO-1,1) = VO(I,JMO-1,1) +
     *    (oDMVA(I,JMO-1)*DXYNO(JMO-1)+
     *    (oDMVA(1,JMO  )*COSIC(I) - oDMUA(1,JMO)*SINIC(I))*DXYSO(JMO)
     *   + oDMVI(I,JMO-1)*DXYVO(JMO-1)) /
     *  (MO(I,JMO-1,1)*DXYNO(JMO-1) + MO(I,JMO,1)*DXYSO(JMO))
      END DO
      end if
      call halo_update(ogrid, odmua, from=north)
      call halo_update(ogrid, odmui, from=north)
      DO J=J_0S,min(J_1S,JMO-2)
        IM1=IMO
        DO I=1,IMO
          IF(LMV(I,J).GT.0.) THEN
            UOD(I,J,1) = UOD(I,J,1) + (
     *           oDMUA(I,J)+oDMUA(I,J+1)
     *     +.5*(oDMUI(IM1,J)+oDMUI(I,J)+oDMUI(IM1,J+1)+oDMUI(I,J+1))
     *           )/(MO(I,J,1)+MO(I,J+1,1))
          ENDIF
          IM1=I
        END DO
      END DO
      RETURN
      END SUBROUTINE OSTRES2

      Subroutine OBDRAG2
!@sum  OBDRAG exerts a drag on the Ocean Model's bottom layer
!@+    define OCN_GISS_TURB and idrag=1 in GISS_OTURB module to include tidal
!@+    enhancement of bottom drag
!@auth Gary Russell and Armando Howard
!@ver  2010/01/08, 2012/03/16
      Use OCEAN, Only: IM,JM,LMO,IVNP,J1O, MO,UO,VO,UOD,VOD, LMU,LMV,
     *                 DTS, COSI=>COSIC,SINI=>SINIC
      Use DOMAIN_DECOMP_1D, Only: HALO_UPDATE, SOUTH,NORTH
      Use OCEANR_DIM,       Only: oGRID
#ifdef OCN_GISS_TURB
      USE GISS_OTURB, only : taubx,tauby, ! x,y components of velocity flux (m/s)^2
     &                       rhobot       ! ocean bottom in-situ density (kg/m^3)
     &                      ,idrag        ! idrag=1: no explicit tides; 0: otherwise
#endif

      Implicit None
      REAL*8, PARAMETER :: BDRAGX=1d0,  ! kg/m^3
     &                     SDRAGX=1d-1
      Integer*4 I,IP1,J,L,J1,JN,JNP
      Real*8    WSQ
      REAL*8 bdragfac   !density times C_D (u^2 + u_t^2)^1/2 (kg/(m^2 s))
#ifdef OCN_GISS_TURB
      REAL*8 taubbyu    !magnitude of velocity flux divided by magnitude of velocity (m/s)
#endif

C**** Define decomposition band parameters
      J1 = oGRID%J_STRT
      JN = oGRID%J_STOP
      JNP = JN  ;  If(JN==JM) JNP = JM-1
      
      Call HALO_UPDATE (oGRID, MO, From=NORTH)
#ifdef OCN_GISS_TURB
      Call HALO_UPDATE (oGRID, rhobot, From=NORTH)
#endif
C****
C**** Reduce ocean current at east edges of cells
C**** UO = UO*(1-x/y)  is approximated by  UO*y/(y+x)  for stability
C**** 
      Do J=Max(J1O,J1),JNP
        I=IM
        DO IP1=1,IM
          IF(LMU(I,J) > 0) THEN
            L=LMU(I,J)
            WSQ = UO(I,J,L)**2 + VOD(I,J,L)**2 + 1d-20
#ifndef OCN_GISS_TURB
            bdragfac=BDRAGX*SQRT(WSQ)
#else
            if(idrag.eq.0) then
              bdragfac=BDRAGX*SQRT(WSQ)
            else
              taubbyu=SQRT((taubx(i,j)**2+tauby(i,j)**2)/WSQ)
              bdragfac=0.5*(rhobot(i,j)+rhobot(ip1,j))*taubbyu
            endif
#endif
            UO(I,J,L) = UO(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *           (MO(I,J,L)+MO(IP1,J,L) + DTS*bdragfac*2d0)
            VOD(I,J,L) = VOD(I,J,L) * (MO(I,J,L)+MO(IP1,J,L)) /
     *           (MO(I,J,L)+MO(IP1,J,L) + DTS*bdragfac*2d0)
          ENDIF
          I=IP1
        ENDDO
      ENDDO
C****
C**** Reduce ocean current at north edges of cells
C****
      Do J=Max(J1O,J1),JNP
        DO I=1,IM
          IF(LMV(I,J) > 0) THEN
            L=LMV(I,J)   
            WSQ = VO(I,J,L)**2 + UOD(I,J,L)**2 + 1d-20
#ifndef OCN_GISS_TURB
            bdragfac=BDRAGX*SQRT(WSQ)
#else
            if(idrag.eq.0) then
              bdragfac=BDRAGX*SQRT(WSQ)
            else
              taubbyu=SQRT((taubx(i,j)**2+tauby(i,j)**2)/WSQ)
              bdragfac=0.5*(rhobot(i,j)+rhobot(i,j+1))*taubbyu
            endif
#endif
            VO(I,J,L) = VO(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *           (MO(I,J,L)+MO(I,J+1,L) + DTS*bdragfac*2d0)
            UOD(I,J,L) = UOD(I,J,L) * (MO(I,J,L)+MO(I,J+1,L)) /
     *           (MO(I,J,L)+MO(I,J+1,L) + DTS*bdragfac*2d0)
          ENDIF
        ENDDO
      ENDDO
      RETURN
C****
      END Subroutine OBDRAG2
