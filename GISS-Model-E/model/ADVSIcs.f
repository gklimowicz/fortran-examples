#include "rundeck_opts.h"

      SUBROUTINE ADVSI(atmice)
!@sum  ADVSI advects sea ice
!@+    Currently set up to advect ice on AGCM grid (i.e. usidt/vsidt are
!@+    on the AGCM grid, and RSI/MSI/HSI etc. are unchanged)
!@+    At some point this will change (USIDT/VSIDT on ice grid, and RSI
!@+    etc. will need to be interpolated back and forth).
!@auth Gary Russell/Gavin Schmidt
!@auth rewrite for cubed sphere: M. Kelley
c NOTE: CURRENTLY ASSUMING THAT THERE IS NO TRANSPORT OF ICE TO/FROM
c EQUATORIAL CUBE FACES.  WILL UPGRADE AS NEEDED.
      USE CONSTANT, only : byshi,lhm,grav
      USE MODEL_COM, only : kocean,dts=>dtsrc
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, 
     &     HALO_UPDATE
      USE GEOM, only : axyp,byaxyp,
     &     dlxsina,dlysina, ull2ucs,vll2ucs, ull2vcs,vll2vcs
      USE ICEDYN_COM, only : foa,byfoa
      USE SEAICE, only : ace1i,xsi
      USE SEAICE_COM, only : si_ocn,lmi
#ifdef TRACERS_WATER
     *     ,ntm
#endif
      USE ICEDYN, only : grid_icdyn,usi,vsi
      USE ICEDYN_COM, only : i2a_uc,i2a_vc,UVLLATUC,UVLLATVC,CONNECT
      use cs2ll_utils, only : ll2csint_lij
      USE EXCHANGE_TYPES, only : atmice_xchng_vars
      IMPLICIT NONE
      type(atmice_xchng_vars) :: atmice
!
!@var NTRICE max. number of tracers to be advected (mass/heat/salt+)
#ifndef TRACERS_WATER
      INTEGER, PARAMETER :: NTRICE=2+2*LMI
#else
      INTEGER, PARAMETER :: NTRICE=2+(2+NTM)*LMI
      INTEGER ITR
      REAL*8 TRSNOW(NTM), TRICE(NTM)
#endif
      INTEGER I,J,L,K,IM
      REAL*8 DMHSI,ASI,YRSI,XRSI,FRSI,SICE,COUR,FAO,CNEW,
     &     ull,vll
C****
C**** USI, VSI  latlon-oriented U,V components of
C****           sea ice velocity (m/s)
C****
C**** FAW    flux of surface water area (m^2) = USIDT*DYP
C**** FASI   flux of sea ice area (m^2) = USIDT*DYP*RSIedge
C**** FMSI   flux of sea ice mass (kg) or heat (J) or salt (kg)
C****
      REAL*8, DIMENSION(ntrice,grid%I_STRT_HALO:grid%I_STOP_HALO) ::
     &     FMSI, FMSI_jm1
      REAL*8, DIMENSION(ntrice) :: AMSI, FMSI_im1
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO) ::
     &     FASI, FXSI, FYSI, FAW,
     &     FASI_jm1, FXSI_jm1, FYSI_jm1, FAW_jm1
      REAL*8 ::
     &     FASI_im1, FXSI_im1, FYSI_im1, FAW_im1

!@var MHS mass/heat/salt content of sea ice
      REAL*8 MHS(NTRICE,grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO)

      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO) ::
     &     UDYDT,VDXDT

      REAL*8, DIMENSION(2,grid_icdyn%im_world,
     &         grid_icdyn%J_STRT_HALO:grid_icdyn%J_STOP_HALO) ::
     &     uvll

      REAL*8, DIMENSION(:,:), POINTER ::
     &     RSI,MSI,SNOWI,FOCEAN,RSIX,RSIY,RSISAVE
      REAL*8, DIMENSION(:,:,:), POINTER :: HSI,SSI
#ifdef TRACERS_WATER
      REAL*8, DIMENSION(:,:,:,:), POINTER :: TRSI
#endif

      INTEGER I_0,I_1,J_0,J_1, I_0Y,I_1Y
!@var coastfac: A proportionality factor to compute the component
!@+   of advective velocity which limits ice buildup along
!@+   coastlines.  (At some gridcells, negative feedbacks on
!@+   ice production are not able to assert themselves when
!@+   sea ice does not reside on the ocean grid.)
      real*8 :: coastfac

C**** Get grid parameters
      IM = grid%im_world
      call getDomainBounds(grid, 
     *     I_STRT=I_0, I_STOP=I_1, J_STRT=J_0, J_STOP=J_1)

      focean => atmice%focean
      rsi => si_ocn%rsi
      rsix => si_ocn%rsix
      rsiy => si_ocn%rsiy
      rsisave => si_ocn%rsisave
      msi => si_ocn%msi
      hsi => si_ocn%hsi
      ssi => si_ocn%ssi
      snowi => si_ocn%snowi
#ifdef tracers_water
      trsi => si_ocn%trsi
#endif

      atmice%musi(:,:)=0
      atmice%husi(:,:)=0
      atmice%susi(:,:)=0
      atmice%mvsi(:,:)=0
      atmice%hvsi(:,:)=0
      atmice%svsi(:,:)=0
#ifdef tracers_water
      atmice%tusi(:,:,:)=0
      atmice%tvsi(:,:,:)=0
#endif

      i_0y = max(1 ,i_0-1)
      i_1y = min(im,i_1+1)

C**** Reduce ice concentration gradients if ice amounts decreased
      DO J=J_0,J_1
      DO I=I_0,I_1
        IF (RSI(I,J).gt.1d-4) THEN
          IF (RSISAVE(I,J).gt.RSI(I,J)) THEN ! reduce gradients
            FRSI=(RSISAVE(I,J)-RSI(I,J))/RSISAVE(I,J)
            RSIX(I,J)=RSIX(I,J)*(1.-FRSI)
            RSIY(I,J)=RSIY(I,J)*(1.-FRSI)
          END IF
          IF(RSI(I,J)-RSIX(I,J).lt.0.)  RSIX(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIX(I,J).lt.0.)  RSIX(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIX(I,J).gt.1d0) RSIX(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIX(I,J).gt.1d0) RSIX(I,J) =1d0-RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).lt.0.)  RSIY(I,J) =    RSI(I,J)
          IF(RSI(I,J)+RSIY(I,J).lt.0.)  RSIY(I,J) =   -RSI(I,J)
          IF(RSI(I,J)-RSIY(I,J).gt.1d0) RSIY(I,J) =    RSI(I,J)-1d0
          IF(RSI(I,J)+RSIY(I,J).gt.1d0) RSIY(I,J) =1d0-RSI(I,J)
        ELSE
          RSIX(I,J) = 0.  ; RSIY(I,J) = 0.
        END IF
C**** update RSISAVE for diagnostics
        RSISAVE(I,J)=RSI(I,J)
C**** set up local MHS array to contain all advected quantities
C**** MHS(1:2) = MASS, MHS(3:2+LMI) = HEAT, MHS(3+LMI:2+2*LMI)=SALT
C**** Currently this is on atmospheric grid
        MHS(1,I,J) = ACE1I + SNOWI(I,J)
        MHS(2,I,J) = MSI(I,J)
        DO L=1,LMI
          MHS(L+2,I,J) = HSI(L,I,J)
          MHS(L+2+LMI,I,J) = SSI(L,I,J)
        END DO
#ifdef TRACERS_WATER
C**** add tracers to advected arrays
        DO ITR=1,NTM
          IF (SNOWI(I,J)*XSI(2).gt.XSI(1)*ACE1I) THEN ! layer 1:all snow
            SICE=SSI(1,I,J)+SSI(2,I,J)
            TRSNOW(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J)*MAX(1.
     *           -(ACE1I-SICE)/(XSI(2)*(ACE1I+SNOWI(I,J))-SICE),0d0)
          ELSE                  ! first layer is snow and some ice
            TRSNOW(ITR) = TRSI(ITR,1,I,J)*MIN(SNOWI(I,J)/(XSI(1)*(ACE1I
     *           +SNOWI(I,J))-SSI(1,I,J)),1d0)
          END IF
          TRICE(ITR) = TRSI(ITR,1,I,J) + TRSI(ITR,2,I,J) - TRSNOW(ITR)
          MHS(1+2+(1+ITR)*LMI,I,J)=TRSNOW(ITR)
          MHS(2+2+(1+ITR)*LMI,I,J)=TRICE(ITR)
          DO L=3,LMI
            MHS(L+2+(1+ITR)*LMI,I,J)=TRSI(ITR,L,I,J)
          END DO
        END DO
#endif
      ENDDO ! i
      ENDDO ! j


      CALL HALO_UPDATE(grid, MSI)

C****
C**** Interpolate to obtain latlon-oriented ice velocities at
C**** cell edges, and transform to CS orientation.  ll2csint_lij
C**** fills any halo cells in its outputs.
C****
      do j=grid_icdyn%j_strt,grid_icdyn%j_stop
        do i=1,grid_icdyn%im_world
          uvll(1,i,j) = usi(i,j)*dts
          uvll(2,i,j) = vsi(i,j)*dts
        enddo
      enddo
      call ll2csint_lij(grid_icdyn,i2a_uc,uvll,uvllatuc,
     &     is_ll_vector=.true.)
      call ll2csint_lij(grid_icdyn,i2a_vc,uvll,uvllatvc,
     &     is_ll_vector=.true.)

      coastfac =
     &          1d-3 ! convert kg/m2 ice mass to ice thickness
     &         *1d-1 ! 10 cm/s speed for 1 m thickness difference over ~100 km
     &         *(real(im,kind=8)/90d0) ! scale by gridlength

      do j=j_0-1,j_1
      do i=i_0y,i_1y
        ull = uvllatvc(1,i,j+1)
        vll = uvllatvc(2,i,j+1)
        vdxdt(i,j) = (ull*ull2vcs(i,j+1)+vll*vll2vcs(i,j+1))
        if(connect(i,j)*connect(i,j+1).eq.0) then
          vdxdt(i,j)=0.
        elseif(connect(i,j)+connect(i,j+1).lt.30) then
          ! 
          if(focean(i,j).lt.focean(i,j+1)) then
            vdxdt(i,j) = vdxdt(i,j)
     &           + dts*min(+10d0,max(0d0,msi(i,j)-msi(i,j+1))*coastfac)
          else
            vdxdt(i,j) = vdxdt(i,j)
     &           + dts*max(-10d0,min(0d0,msi(i,j)-msi(i,j+1))*coastfac)
          endif
        endif
        vdxdt(i,j) = dlxsina(i,j+1)*vdxdt(i,j)
      enddo
      enddo
      do j=j_0,j_1
      do i=i_0-1,i_1
        ull = uvllatuc(1,i+1,j)
        vll = uvllatuc(2,i+1,j)
        udydt(i,j) = (ull*ull2ucs(i+1,j)+vll*vll2ucs(i+1,j))
        if(connect(i,j)*connect(i+1,j).eq.0) then
          udydt(i,j)=0.
        elseif(connect(i,j)+connect(i+1,j).lt.30) then
          if(focean(i,j).lt.focean(i+1,j)) then
            udydt(i,j) = udydt(i,j)
     &           + dts*min(+10d0,max(0d0,msi(i,j)-msi(i+1,j))*coastfac)
          else
            udydt(i,j) = udydt(i,j)
     &           + dts*max(-10d0,min(0d0,msi(i,j)-msi(i+1,j))*coastfac)
          endif
        endif
        udydt(i,j) = dlysina(i+1,j)*udydt(i,j)
      enddo
      enddo

c for now, no transport across cube edges
      if(i_0 .eq.  1) udydt(0 ,:) = 0.
      if(i_1 .eq. im) udydt(im,:) = 0.
      if(j_0 .eq.  1) vdxdt(:, 0) = 0.
      if(j_1 .eq. im) vdxdt(:,im) = 0.

C****
C**** Update halos of transported quantities
C****
      CALL HALO_UPDATE(grid, RSI)
      CALL HALO_UPDATE(grid, RSIX)
      CALL HALO_UPDATE(grid, RSIY)
      CALL HALO_UPDATE(grid, MHS, jdim=3)

C****
C**** Transport in the Y direction
C****
      j = j_0-1
      do i=i_0y,i_1y ! updating i-halo cells for subsequent x-sweep
        faw(i) = vdxdt(i,j)  ! should compute vdxdt here
        if(faw(i).le.0.) then
c**** sea ice velocity is southward at grid box edge
          cour = faw(i)*byaxyp(i,j+1)
          fao = faw(i)*focean(i,j+1)
          fasi(i)=fao*(rsi(i,j+1)-(1d0+cour)*rsiy(i,j+1))
          fxsi(i)=fao*rsix(i,j+1)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j+1)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j+1)
        else
c**** sea ice velocity is northward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsiy(i,j))
          fxsi(i)=fao*rsix(i,j)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        end if
        faw_jm1(i)  = faw(i)
        fasi_jm1(i) = fasi(i)
        fxsi_jm1(i) = fxsi(i)
        fysi_jm1(i) = fysi(i)
        fmsi_jm1(1:ntrice,i) = fmsi(1:ntrice,i)
      enddo

      do j=j_0,j_1
      do i=i_0y,i_1y ! updating i-halo cells for subsequent x-sweep

c _jm1 qtys are already zero. no need to re-zero.
        if(vdxdt(i,j-1).eq.0. .and. vdxdt(i,j).eq.0.) cycle

        faw(i) = vdxdt(i,j)  ! should compute vdxdt here
        if(faw(i).le.0.) then
c**** sea ice velocity is southward at grid box edge
          cour = faw(i)*byaxyp(i,j+1)
          fao = faw(i)*focean(i,j+1)
          fasi(i)=fao*(rsi(i,j+1)-(1d0+cour)*rsiy(i,j+1))
          fxsi(i)=fao*rsix(i,j+1)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j+1)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j+1)
        else
c**** sea ice velocity is northward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsiy(i,j))
          fxsi(i)=fao*rsix(i,j)
          fysi(i)=faw(i)*(cour*fao*rsiy(i,j)-3d0*fasi(i))
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        end if
c accumulate transports
        atmice%mvsi(i,j)=sum(fmsi(1:2,i))
        atmice%hvsi(i,j)=sum(fmsi(3:2+lmi,i))
        atmice%svsi(i,j)=sum(fmsi(3+lmi:2+2*lmi,i))
#ifdef TRACERS_WATER
        do itr=1,ntm
          atmice%tvsi(i,j,itr)=
     &         sum(fmsi(3+(1+itr)*lmi:2+(2+itr)*lmi,i))
        enddo
#endif
        if(faw_jm1(i).gt.0. .or. faw(i).lt.0.) then
! when there is inflow, use general-case formulas
          asi = rsi(i,j)*foa(i,j)
          do k=1,ntrice
            amsi(k) = asi*mhs(k,i,j) + (fmsi_jm1(k,i)-fmsi(k,i))
          enddo
          asi = asi + (fasi_jm1(i)-fasi(i))
          if(asi.le.foa(i,j)) then
            yrsi = (rsiy(i,j)*axyp(i,j)*foa(i,j)+
     &           (fysi_jm1(i)-fysi(i)) + 3d0*((faw_jm1(i)+
     &           faw(i))*asi-axyp(i,j)*(fasi_jm1(i)+fasi(i))))
     &           / (axyp(i,j) + (faw_jm1(i)-faw(i)))
            rsi(i,j)  = asi*byfoa(i,j)
            rsiy(i,j) = yrsi*byfoa(i,j)
            rsix(i,j) = rsix(i,j) + (fxsi_jm1(i)-fxsi(i))*byfoa(i,j)
            if(asi.gt.0) mhs(1:ntrice,i,j) = amsi(1:ntrice)/asi
          else
c**** sea ice crunches into itself and completely covers grid box
            rsi(i,j)   = 1d0
            rsix(i,j)  = 0.
            rsiy(i,j)  = 0.
            mhs(1,i,j) = amsi(1)/asi
            mhs(2,i,j) =(amsi(1)+amsi(2))*byfoa(i,j) - mhs(1,i,j)
            do k=1,(ntrice-2)/lmi
              mhs(3+lmi*(k-1),i,j) = amsi(3+lmi*(k-1)) / asi
              mhs(4+lmi*(k-1),i,j) = amsi(4+lmi*(k-1)) / asi
              dmhsi = (amsi(3+lmi*(k-1))+amsi(4+lmi*(k-1))
     &             +amsi(5+lmi*(k-1))
     &             +amsi(6+lmi*(k-1)))*(byfoa(i,j) -1d0 / asi )
              mhs(5+lmi*(k-1),i,j) = amsi(5+lmi*(k-1)) / asi +
     &             xsi(3)*dmhsi
              mhs(6+lmi*(k-1),i,j) = amsi(6+lmi*(k-1)) / asi +
     &             xsi(4)*dmhsi
            end do

          endif
        else
! when there is only outflow, use simpler formulas.
! why is mhs not updated here.
          rsi(i,j)  =  rsi(i,j) + (fasi_jm1(i)-fasi(i))*byfoa(i,j)
          cnew = 1d0+(faw_jm1(i)*focean(i,j-1)
     &               -faw(i)*focean(i,j))*byfoa(i,j)
          rsix(i,j) = rsix(i,j)*cnew
          rsiy(i,j) = rsiy(i,j)*cnew**2
        endif

        faw_jm1(i)  = faw(i)
        fasi_jm1(i) = fasi(i)
        fxsi_jm1(i) = fxsi(i)
        fysi_jm1(i) = fysi(i)
        fmsi_jm1(1:ntrice,i) = fmsi(1:ntrice,i)

C**** Limit RSIX and RSIY so that sea ice is positive at the edges.
        rsi(i,j) = max(0d0,rsi(i,j))
        if(rsi(i,j)-rsix(i,j).lt.0.)  rsix(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsix(i,j).lt.0.)  rsix(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsix(i,j).gt.1d0) rsix(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsix(i,j).gt.1d0) rsix(i,j) =1d0-rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).lt.0.)  rsiy(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsiy(i,j).lt.0.)  rsiy(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).gt.1d0) rsiy(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsiy(i,j).gt.1d0) rsiy(i,j) =1d0-rsi(i,j)

      enddo ! i
      enddo ! j

C****
C**** Transport in the X direction
C****

      do j=j_0,j_1

      i=i_0-1
      faw(i) = udydt(i,j)  ! should compute udydt here
      if(faw(i).le.0.) then
c**** sea ice velocity is westward at grid box edge
        cour = faw(i)*byaxyp(i+1,j)
        fao = faw(i)*focean(i+1,j)
        fasi(i)=fao*(rsi(i+1,j)-(1d0+cour)*rsix(i+1,j))
        fxsi(i)=faw(i)*(cour*fao*rsix(i+1,j)-3d0*fasi(i))
        fysi(i)=fao*rsiy(i+1,j)
        fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i+1,j)
      else
c**** sea ice velocity is eastward at grid box edge
        cour = faw(i)*byaxyp(i,j)
        fao = faw(i)*focean(i,j)
        fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsix(i,j))
        fxsi(i)=faw(i)*(cour*fao*rsix(i,j)-3d0*fasi(i))
        fysi(i)=fao*rsiy(i,j)
        fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
      endif
      faw_im1 = faw(i)
      fasi_im1 = fasi(i)
      fxsi_im1 = fxsi(i)
      fysi_im1 = fysi(i)
      fmsi_im1(1:ntrice) = fmsi(1:ntrice,i)

      do i=i_0,i_1

c _im1 qtys are already zero. no need to re-zero.
        if(udydt(i-1,j).eq.0. .and. udydt(i,j).eq.0.) cycle

        faw(i) = udydt(i,j) ! should compute udydt here
        if(faw(i).le.0.) then
c**** sea ice velocity is westward at grid box edge
          cour = faw(i)*byaxyp(i+1,j)
          fao = faw(i)*focean(i+1,j)
          fasi(i)=fao*(rsi(i+1,j)-(1d0+cour)*rsix(i+1,j))
          fxsi(i)=faw(i)*(cour*fao*rsix(i+1,j)-3d0*fasi(i))
          fysi(i)=fao*rsiy(i+1,j)
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i+1,j)
        else
c**** sea ice velocity is eastward at grid box edge
          cour = faw(i)*byaxyp(i,j)
          fao = faw(i)*focean(i,j)
          fasi(i)=fao*(rsi(i,j)+(1d0-cour)*rsix(i,j))
          fxsi(i)=faw(i)*(cour*fao*rsix(i,j)-3d0*fasi(i))
          fysi(i)=fao*rsiy(i,j)
          fmsi(1:ntrice,i) = fasi(i)*mhs(1:ntrice,i,j)
        endif
c accumulate transports
        atmice%musi(i,j)=sum(fmsi(1:2,i))
        atmice%husi(i,j)=sum(fmsi(3:2+lmi,i))
        atmice%susi(i,j)=sum(fmsi(3+lmi:2+2*lmi,i))
#ifdef TRACERS_WATER
        do itr=1,ntm
          atmice%tusi(i,j,itr)=
     &         sum(fmsi(3+(1+itr)*lmi:2+(2+itr)*lmi,i))
        enddo
#endif
        if(faw_im1.gt.0. .or. faw(i).lt.0.) then
! when there is inflow, use general-case formulas
          asi = rsi(i,j)*foa(i,j)
          do k=1,ntrice
            amsi(k) = asi*mhs(k,i,j) + (fmsi_im1(k)-fmsi(k,i))
          enddo
          asi = asi + (fasi_im1-fasi(i))
          if(asi.le.foa(i,j)) then
            xrsi = (rsix(i,j)*axyp(i,j)*foa(i,j)+
     &        (fxsi_im1-fxsi(i)) + 3d0*((faw_im1+faw(i))*asi-
     &           axyp(i,j)*(fasi_im1+fasi(i))))
     &           / (axyp(i,j) + (faw_im1-faw(i)))
            rsi(i,j)  = asi*byfoa(i,j)
            rsix(i,j) = xrsi*byfoa(i,j)
            rsiy(i,j) = rsiy(i,j) + (fysi_im1-fysi(i))*byfoa(i,j)
            if (asi.gt.0) mhs(1:ntrice,i,j) = amsi(1:ntrice)/asi
          else
c**** sea ice crunches into itself and completely covers grid box
            rsi(i,j)   = 1d0
            rsix(i,j)  = 0.
            rsiy(i,j)  = 0.
            mhs(1,i,j) = amsi(1)/asi
            mhs(2,i,j) =(amsi(1)+amsi(2))*byfoa(i,j) - mhs(1,i,j)
            do k=1,(ntrice-2)/lmi
              mhs(3+lmi*(k-1),i,j) = amsi(3+lmi*(k-1)) / asi
              mhs(4+lmi*(k-1),i,j) = amsi(4+lmi*(k-1)) / asi
              dmhsi = (amsi(3+lmi*(k-1))+amsi(4+lmi*(k-1))
     &             +amsi(5+lmi*(k-1))
     &             +amsi(6+lmi*(k-1)))*(byfoa(i,j) -1d0/ asi)
              mhs(5+lmi*(k-1),i,j) = amsi(5+lmi*(k-1)) / asi +
     &             xsi(3)*dmhsi
              mhs(6+lmi*(k-1),i,j) = amsi(6+lmi*(k-1)) / asi +
     &             xsi(4)*dmhsi
            enddo
          endif
        else
! when there is only outflow, use simpler formulas
! why is mhs not updated here.
          rsi(i,j)  =  rsi(i,j) + (fasi_im1-fasi(i))*byfoa(i,j)
          cnew = 1d0+(faw_im1*focean(i-1,j)
     &               -faw(i )*focean(i  ,j))*byfoa(i,j)
          rsix(i,j) = rsix(i,j)*cnew**2
          rsiy(i,j) = rsiy(i,j)*cnew
        endif
        faw_im1  = faw(i)
        fasi_im1 = fasi(i)
        fxsi_im1 = fxsi(i)
        fysi_im1 = fysi(i)
        fmsi_im1(1:ntrice) = fmsi(1:ntrice,i)
        
C**** Limit RSIX and RSIY so that sea ice is positive at the edges.
c why is it necessary to do this after the advection?
        rsi(i,j) = max(0d0,rsi(i,j))
        if(rsi(i,j)-rsix(i,j).lt.0.)  rsix(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsix(i,j).lt.0.)  rsix(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsix(i,j).gt.1d0) rsix(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsix(i,j).gt.1d0) rsix(i,j) =1d0-rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).lt.0.)  rsiy(i,j) =    rsi(i,j)
        if(rsi(i,j)+rsiy(i,j).lt.0.)  rsiy(i,j) =   -rsi(i,j)
        if(rsi(i,j)-rsiy(i,j).gt.1d0) rsiy(i,j) =    rsi(i,j)-1d0
        if(rsi(i,j)+rsiy(i,j).gt.1d0) rsiy(i,j) =1d0-rsi(i,j)

      enddo ! i
      enddo ! j

      IF (KOCEAN.ge.1) THEN ! full ocean calculation, adjust sea ice
C**** set global variables from local array
C**** Currently on atmospheric grid, so no interpolation necessary
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
C**** Fresh water sea ice mass convergence (needed for qflux model)
            atmice%MSICNV(I,J) =
     &             RSI(I,J)*(MHS(1,I,J)+MHS(2,I,J)-SUM(MHS(3
     &           +LMI:2*LMI+2,I,J))) - RSISAVE(I,J)*(ACE1I+SNOWI(I,J)
     &           +MSI(I,J)-SUM(SSI(1:LMI,I,J)))
C**** sea ice prognostic variables
            SNOWI(I,J)= MAX(0d0,MHS(1,I,J) - ACE1I)
            MSI(I,J)  = MHS(2,I,J)
            DO L=1,LMI
              HSI(L,I,J) = MHS(L+2,I,J)
            END DO
C**** ensure that salinity is only associated with ice
            SICE=MHS(1+2+LMI,I,J)+MHS(2+2+LMI,I,J)
            IF (SNOWI(I,J).gt.XSI(2)*(ACE1I+SNOWI(I,J))) THEN
              SSI(1,I,J)=0.
            ELSE
              SSI(1,I,J)=SICE*(XSI(1)*ACE1I-XSI(2)*SNOWI(I,J))/ACE1I
            END IF
            SSI(2,I,J)=SICE-SSI(1,I,J)
C**** correction of heat energy to compensate for salinity fix
            HSI(1,I,J)=HSI(1,I,J)-(MHS(1+2+LMI,I,J)-SSI(1,I,J))*LHM
            HSI(2,I,J)=HSI(2,I,J)+(MHS(1+2+LMI,I,J)-SSI(1,I,J))*LHM
            DO L=3,LMI
               SSI(L,I,J) = MHS(L+2+LMI,I,J)
            END DO
#ifdef TRACERS_WATER
C**** reconstruct tracer arrays
            DO ITR=1,NTM
              IF (ACE1I.gt.XSI(2)*(SNOWI(I,J)+ACE1I)) THEN
                TRSI(ITR,1,I,J)= MHS(2+2+(1+ITR)*LMI,I,J) *(ACE1I
     *               -XSI(2)*(SNOWI(I,J)+ACE1I))/ACE1I +MHS(1+2+(1+ITR)
     *               *LMI,I,J)
              ELSE
                TRSI(ITR,1,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)*XSI(1)*(ACE1I
     *               +SNOWI(I,J))/SNOWI(I,J)
              END IF
              TRSI(ITR,2,I,J)= MHS(1+2+(1+ITR)*LMI,I,J)+MHS(2+2+(1+ITR)
     *             *LMI,I,J)-TRSI(ITR,1,I,J)
              DO L=3,LMI
                TRSI(ITR,L,I,J)=MHS(L+2+(1+ITR)*LMI,I,J)
              END DO
            END DO
#endif
            atmice%FWSIM(I,J)=RSI(I,J)*(ACE1I+SNOWI(I,J)+MSI(I,J)-
     *           SUM(SSI(1:LMI,I,J)))
            END IF
          END DO
        END DO
      ELSE          ! fixed SST case, save implied heat convergence
        DO J=J_0,J_1
          DO I=I_0,I_1
            IF (FOCEAN(I,J).gt.0) THEN
              atmice%HSICNV(I,J)=(RSI(I,J)*SUM(MHS(3:2+LMI,I,J))
     *             -RSISAVE(I,J)*SUM(HSI(1:LMI,I,J)))
C**** reset sea ice concentration
              RSI(I,J)=RSISAVE(I,J)
            END IF
          END DO
        END DO
      END IF
C****
      RETURN
      END SUBROUTINE ADVSI
