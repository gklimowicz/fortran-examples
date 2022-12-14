#include "rundeck_opts.h"

      module ATMDYN
      implicit none

      contains

      SUBROUTINE init_ATMDYN
      return
      end SUBROUTINE init_ATMDYN

      SUBROUTINE DYNAM
      Use RESOLUTION, Only: IM,JM,LM
      USE SOMTQ_COM,  only: tmom,mz
      USE ATM_COM,    only: MA,t,q,PMID,PEDN,MUs,MVs,MWs
      USE DOMAIN_DECOMP_ATM, only : grid

      Real*8 :: TZ(IM,grid%J_STRT_HALO:grid%J_STOP_HALO,LM)
      INTEGER L

      do L=1,LM
         MUs(:,:,L) = 0.
         MVs(:,:,L) = 0.
         MWs(:,:,L) = 0.
      ENDDO

      CALL CALC_AMPK(LM)

      call SCM_FORCN
  
      CALL tq_zmom_init(T,Q,PMID,PEDN)

      DO L=1,LM
         TZ(:,:,L)  = TMOM(MZ,:,:,L)
      ENDDO

      Call PGF_SCM (T,TZ,MA)

      return
      END SUBROUTINE DYNAM


      SUBROUTINE SCM_FORCN
c     apply large-scale forcings to T, Q, U, V

      USE MODEL_COM,  only: DTSRC
      USE ATM_COM,    only: P,T,Q,PK,U,V,PMID
      USE RESOLUTION, only: LM
      USE CONSTANT,   only: KAPA, OMEGA, GRAV, RGAS
      USE GEOM, only : sinlat2d
      USE SCM_COM,    only: SCMopt,SCMin
#ifdef CACHED_SUBDD
      USE DOMAIN_DECOMP_ATM, only : grid,getDomainBounds
      use subdd_mod, only : subdd_groups,subdd_ngroups,subdd_type
     &     ,inc_subdd,find_groups
#endif

      IMPLICIT NONE

      real*8, dimension(LM) :: Tabs(LM),
     &                         SCM_ver_u_adv,SCM_ver_v_adv,
     &                         SCM_nudge_T,SCM_nudge_Q,
     &                         SCM_force_T,SCM_force_Q,
     &                         SCM_force_U,SCM_force_V
      real*8 f_cor

      INTEGER L
#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      integer :: i,i_0,i_1,j,j_0,j_1,k
      type(subdd_type), pointer :: subdd
      REAL*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &                  sddarr3d
      CALL getDomainBounds(grid,
     &     I_STRT=I_0,I_STOP=I_1,J_STRT=J_0,J_STOP=J_1)
#endif

c     operate on absolute temperature
      do L = 1,LM
        Tabs(L) = T(1,1,L)*PK(L,1,1)
      enddo

      SCM_force_T = 0.
      SCM_force_Q = 0.
      SCM_force_U = 0.
      SCM_force_V = 0.
      SCM_nudge_T = 0.
      SCM_nudge_Q = 0.

      do L = 1,LM

c       fix winds if specified and no Coriolis acceleration

        if( SCMopt%wind .and. .not. SCMopt%geo )then
          U(1,1,L) = SCMin%U(L)
          V(1,1,L) = SCMin%V(L)
        endif

c       large-scale forcings

        if( ( SCMopt%omega .and. .not. SCMopt%ls_v ) .or. SCMopt%w )then
c       *** apply omega defined at layer bottom to upwind gradient

          SCMin%SadvV(L) = 0.
          SCMin%QadvV(L) = 0.

          if ( L < LM ) then ! omega assumed zero at top of layer LM
            if ( SCMin%Omega(L+1) > 0. ) then ! upwind gradient above
              SCMin%SadvV(L) = -SCMin%Omega(L+1)*
     &           (T(1,1,L+1)-T(1,1,L))*PK(L,1,1)/
     &           (PMID(L+1,1,1)-PMID(L,1,1))
              SCMin%QadvV(L) = -SCMin%Omega(L+1)*
     &           (Q(1,1,L+1)-Q(1,1,L))/
     &           (PMID(L+1,1,1)-PMID(L,1,1))
            endif
          endif

          if ( L > 1 ) then ! no atmospheric gradient through surface
            if ( SCMin%Omega(L) < 0. ) then ! upwind gradient below
              SCMin%SadvV(L) = -SCMin%Omega(L)*
     &           (T(1,1,L)-T(1,1,L-1))*PK(L,1,1)/
     &           (PMID(L,1,1)-PMID(L-1,1,1))
              SCMin%QadvV(L) = -SCMin%Omega(L)*
     &           (Q(1,1,L)-Q(1,1,L-1))/
     &           (PMID(L,1,1)-PMID(L-1,1,1))
            endif
          endif

        else if( .not. SCMopt%ls_v )then
c       *** otherwise no LS vertical flux divergence if not specified
          SCMin%SadvV(L) = 0.
          SCMin%QadvV(L) = 0.
        endif

        if( .not. SCMopt%ls_h )then
c       *** no horizontal forcings
          SCMin%TadvH(L) = 0.
          SCMin%QadvH(L) = 0.
        endif

        SCM_force_T(L) = (SCMin%TadvH(L)+SCMin%SadvV(L))*DTSRC
        SCM_force_Q(L) = (SCMin%QadvH(L)+SCMin%QadvV(L))*DTSRC

        Tabs(L) = Tabs(L) + SCM_force_T(L)
        Q(1,1,L) = Q(1,1,L) + SCM_force_Q(L)

        if( Q(1,1,L) < 0. )then
          SCM_force_Q(L) = -Q(1,1,L)
          Q(1,1,L) = 0.0
        endif

c       apply nudging terms

        if( SCMopt%nudge )then
c       *** calculate nudging toward observed profile

          SCM_nudge_T(L) = (SCMin%T(L)-Tabs(L))/SCMopt%tau*DTSRC
          SCM_nudge_Q(L) = (SCMin%Q(L)-Q(1,1,L))/SCMopt%tau*DTSRC

          if( SCMopt%Fnudge )then
            SCM_nudge_T(L) = SCM_nudge_T(L)*SCMin%Fnudge(L)
            SCM_nudge_Q(L) = SCM_nudge_Q(L)*SCMin%Fnudge(L)
          endif

          Tabs(L) = Tabs(L) + SCM_nudge_T(L)
          Q(1,1,L) = Q(1,1,L) + SCM_nudge_Q(L)

        endif

      enddo

c     *** apply changes to actual prognostic variable (potential temperature)
      do L = 1,LM
        T(1,1,L) = Tabs(L)/PK(L,1,1)
      enddo

c     when Coriolis forcing used (computed from geostrophic winds),
c     also possibly apply vertical advection to horizontal winds

      if ( SCMopt%geo ) then

        f_cor = 2.*omega*sinlat2d(1,1)

        SCM_ver_u_adv = 0.
        SCM_ver_v_adv = 0.

        if( .not. SCMopt%ls_h_UV )then
c       *** no horizontal wind forcings
          SCMin%UadvH(:) = 0.
          SCMin%VadvH(:) = 0.
        endif

        do L = 1,LM

          if( SCMopt%VadvHwind )then
c         *** apply omega defined at layer bottom to upwind gradient

            if ( L < LM ) then ! omega assumed zero at top of layer LM
              if ( SCMin%Omega(L+1) > 0. ) then ! upwind gradient above
                 SCM_ver_u_adv(L) = -SCMin%Omega(L+1)*
     &              (U(1,1,L+1)-U(1,1,L))/
     &              (PMID(L+1,1,1)-PMID(L,1,1))
                 SCM_ver_v_adv(L) = -SCMin%Omega(L+1)*
     &              (V(1,1,L+1)-V(1,1,L))/
     &              (PMID(L+1,1,1)-PMID(L,1,1))
              endif
            endif

            if ( L > 1 ) then ! no atmospheric gradient through surface
              if ( SCMin%Omega(L) < 0. ) then ! upwind gradient below
                 SCM_ver_u_adv(L) = -SCMin%Omega(L)*
     &              (U(1,1,L)-U(1,1,L-1))/
     &              (PMID(L,1,1)-PMID(L-1,1,1))
                 SCM_ver_v_adv(L) = -SCMin%Omega(L)*
     &              (V(1,1,L)-V(1,1,L-1))/
     &              (PMID(L,1,1)-PMID(L-1,1,1))
              endif
            endif
          endif

          SCM_force_U(L) = (SCM_ver_u_adv(L)+SCMin%UadvH(L))*DTSRC
          SCM_force_V(L) = (SCM_ver_v_adv(L)+SCMin%VadvH(L))*DTSRC

c         apply combined forcings to horizontal winds
          U(1,1,L) = U(1,1,L) + SCM_force_U(L) +
     &        f_cor*(V(1,1,L)-SCMin%Vg(L))*DTSRC
          V(1,1,L) = V(1,1,L) + SCM_force_V(L) -
     &        f_cor*(U(1,1,L)-SCMin%Ug(L))*DTSRC

        enddo      ! L = 1,LM
      endif        ! use geostrophic winds for Coriolis forcing

#ifdef CACHED_SUBDD
C****
C**** Collect some high-frequency outputs
C****
      call find_groups('fijlh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      do k=1,subdd%ndiags
      select case (subdd%name(k))
      case ('dq_ls')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_force_Q(l)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dth_ls')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_force_T(l)/PK(l,1,1)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('du_ls')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_force_U(l)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dv_ls')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_force_V(l)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dq_nudge')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_nudge_Q(l)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      case ('dth_nudge')
        do j=j_0,j_1; do i=i_0,i_1; do l=1,lm
          sddarr3d(i,j,l) = SCM_nudge_T(l)/PK(l,1,1)
        enddo;        enddo;        enddo
        call inc_subdd(subdd,k,sddarr3d)
      end select
      enddo
      enddo
#endif

      return
      END SUBROUTINE SCM_FORCN 

  
      SUBROUTINE SDRAG(DT1)
      REAL*8, INTENT(IN) :: DT1 
      return
      END SUBROUTINE SDRAG 


      Subroutine PGF_SCM (S0,SZ,MAM)
!@SCM-version   Computes geopotential consistent with PGF 
!**** Input: MAM = mean mass distribution during time step (kg/m^2)
!****      S0,SZ = potential temperature and vertical gradient (K)
!**** Output: GZ = PHI (m^2/s^2) = atmospheric geopotential

!**** R (J/kg*C) = gas constant = 287 for dry air
!**** K          = exponent of exner function = R/SHA
!**** M (kg/m^2) = vertical coordinate = air mass above the level
!**** DM(kg/m^2) = layer mass difference = MAM
!**** P (Pa)     = pressure = M*GRAV
!**** DP(Pa)     = layer pressure difference = PD - PU
!**** A (m^3/kg) = specific volume = R*T / P
!**** S (K)      = potential temperature = S0 - SZ*2*(M-M0)/(MD-MU) =
!****            = S0 - SZ*2*(P-P0)/(PD-PU) = S0 - SZ*2*(P-P0)/DP =
!****            = S0+SZ*2*P0/DP - SZ*2*P/DP
!**** T (K)      = temperature = S * P(mb)^K = S*.01^K * P^K =
!****            = [(S0+SZ*2*P0/DP)*.01^K - P*(SZ*2/DP)*.01^K]*P^K =
!****            = (X - P*Y)*P^K = X*P^K - Y*P^(K+1)

!**** Integral of A*dM from MU to MD (from top to bottom of layer)
!**** Int[A*dM] = Int[R*T*dP/P*G] = R*Int{[X*P^(K-1) - Y*P^K]*dP}/G =
!**** = R*{X*P^K/K - Y*P^(K+1)/(K+1)}/G from PU to PD =
!**** = R*{X*(PD^K-PU^K)/K - Y*[PD^(K+1)-PU^(K+1)]/(K+1)}/G

!**** Compute DGZ thickness everwhere in a layer from layer bottom
!**** G*dZ = - A*dP = - (R*T/P)*dP = - R*[X*P^(K-1) - Y*P^K]*dP 
!**** DGZ = - Int{R*[X*P^(K-1) - Y*P^K]*dP} from PD to P =
!****     = - R*{X*(P^K-PD^K)/K - Y*[P^(K+1)-PD^(K+1)]/(K+1)}
!****     = R*{X*(PD^K-P^K)/K - Y*[PD^(K+1)-P^(K+1)]/(K+1)}

!**** DGZup = R*{X*(PD^K-PU^K)/K - Y*[PD^(K+1)-PU^(K+1)]/(K+1)}
!**** Int[A*dM] = DGZup/G

!**** Compute mass weighted average value of DGZ in a layer
!**** DGZave = Int{DGZ*dP}/DP from PD to PU =
!**** = R*Int({X*(PD^K-P^K)/K - Y*[PD^(K+1)-P^(K+1)]/(K+1)}*dP)/DP =
!**** = R*{X*[P*PD^K - P^(K+1)/(K+1)]/K -
!****    - Y*[P*PD^(K+1) - P^(K+2)/(K+2)]/(K+1)}/DP =
!**** = R*(X*{(PD-PU)*PD^K - [PD^(K+1)-PU^(K+1)]/(K+1)}/K -
!****    - Y*{(PD-PU)*PD^(K+1) - [PD^(K+2)-PU^(K+2)]/(K+2)}/(K+1))/DP =
!**** = R*(X*{DP*PD^K - [PD^(K+1)-PU^(K+1)]/(K+1)}/K -
!****    - Y*{DP*PD^(K+1) - [PD^(K+2)-PU^(K+2)]/(K+2)}/(K+1))/DP =

!**** GZave(L) = GZATMO + Sum[DGZup(1:L-1)] + DGZave(L)  

      Use CONSTANT,   Only: GRAV,RGAS,KAPA,byGRAV,
     *                      zK=>byKAPA,zKp1=>byKAPAp1,zKp2=>byKAPAp2
      Use RESOLUTION, Only: IM,JM,LM, MTOP
      Use ATM_COM,    Only: ZATMO, GZ,PHI
      Implicit None
      Real*8,Dimension(1,1,LM) :: S0,SZ
      Real*8,Dimension(LM,1,1) :: MAM
!**** Local variables
      Real*8  :: DGZU(LM),DGZA(LM),
     *           M,PU,PKU,PKPU,PKPPU,DP,zDP,X,Y,PD,PKD,PKPD,PKPPD,GZD,
     *           HUNDREDTHeKAPA
      Integer :: I,J,L

      HUNDREDTHeKAPA = .01d0**KAPA
      I=1  ;  J=1
!**** Integrate pressures from the top down
      M   = MTOP
      PU  = M*GRAV
      PKU = PU**KAPA  ;  PKPU = PKU*PU  ;  PKPPU = PKPU*PU       
      Do L=LM,1,-1
         DP  = MAM(L,I,J)*GRAV
         zDP = 1 / DP
         Y   = SZ(I,J,L)*2*zDP*HUNDREDTHeKAPA
         X   = S0(I,J,L)*HUNDREDTHeKAPA + Y*(PU+.5*DP)
         PD  = PU + DP
         PKD = PD**KAPA  ;  PKPD = PKD*PD  ;  PKPPD = PKPD*PD
!        AdM = RGAS*(X*(PKD-PKU)*zK - Y*(PKPD-PKPU)*zKp1)/GRAV
         DGZU(L) = RGAS*(X*(PKD-PKU)*zK - Y*(PKPD-PKPU)*zKp1)
         DGZA(L) = RGAS*(X*(DP*PKD - (PKPD-PKPU)*zKp1)*zK -
     -                   Y*(DP*PKPD - (PKPPD-PKPPU)*zKp2)*zKp1)*zDP
         M   = M + MAM(L,I,J)
         PU  = PD
         PKU = PKD  ;  PKPU = PKPD  ;  PKPPU=PKPPD  ;  EndDo     
!**** Integrate altitude from the bottom up
      GZD = ZATMO(I,J)
      Do L=1,LM
         GZ(I,J,L) = GZD + DGZA(L)
         GZD = GZD + DGZU(L)  ;  EndDo

      PHI(:,:,:) = GZ(:,:,:)
      Return
      End Subroutine PGF_SCM


C**** Dummy routines

      SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS( MUs,MVs,dt)
!@sum COMPUTE_DYNAM_AIJ_DIAGNOSTICS Dummy
      use DOMAIN_DECOMP_ATM, only: grid

      real*8, intent(in) :: MUs(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: MVs(:,grid%J_STRT_HALO:,:)
      real*8, intent(in) :: dt

      return
      END SUBROUTINE COMPUTE_DYNAM_AIJ_DIAGNOSTICS

      end module ATMDYN

      SUBROUTINE conserv_KE(RKE)
!@sum  conserv_KE calculates A-grid column-sum atmospheric kinetic energy,
!@sum  multiplied by cell area
!@auth Gary Russell/Gavin Schmidt
      IMPLICIT NONE

      REAL*8, DIMENSION(1,1) :: RKE

      RKE = 0.
      !call stop_model('calculate a-grid value instead',255)

      RETURN
C****
      END SUBROUTINE conserv_KE

      SUBROUTINE calc_kea_3d(kea)
!@sum  calc_kea_3d calculates square of wind speed on the A grid
      USE RESOLUTION, only: lm
      USE ATM_COM,    only: u,v
      IMPLICIT NONE
      REAL*8, DIMENSION(1,1,LM) :: KEA

      RETURN

      END SUBROUTINE calc_kea_3d

      subroutine recalc_agrid_uv
      USE ATM_COM,    only: u,v
      USE ATM_COM,    only: ua=>ualij,va=>valij
      implicit none

      ua(:,1,1)=u(1,1,:)
      va(:,1,1)=v(1,1,:)

      return
      end subroutine recalc_agrid_uv

      subroutine replicate_uv_to_agrid(ur,vr,k,ursp,vrsp,urnp,vrnp)
      USE RESOLUTION, only: lm
      USE ATM_COM,    only: u,v
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,1,1) :: UR,VR
      real*8, dimension(1,lm) :: ursp,vrsp,urnp,vrnp ! not used
      integer :: l
      if(k.ne.1)
     &     call stop_model('incorrect k in replicate_uv_to_agrid',255)
      do l=1,lm
        ur(1,l,1,1) = u(1,1,l)
        vr(1,l,1,1) = v(1,1,l)
      enddo ! l
      return
      end subroutine replicate_uv_to_agrid

      subroutine avg_replicated_duv_to_vgrid(du,dv,k,
     &     dusp,dvsp,dunp,dvnp)
      USE RESOLUTION, only: lm
      USE ATM_COM,    only: u,v
      implicit none
      integer :: k
      REAL*8, DIMENSION(k,LM,1,1) :: DU,DV
      real*8, dimension(1,lm) :: dusp,dvsp,dunp,dvnp ! not used
      integer :: l

      if(k.ne.1) call stop_model(
     &     'incorrect k in avg_replicated_duv_to_vgrid',255)

      do l=1,lm
        u(1,1,l)=u(1,1,l)+du(1,l,1,1)
        v(1,1,l)=v(1,1,l)+dv(1,l,1,1)
      enddo ! l

      return
      end subroutine avg_replicated_duv_to_vgrid

      SUBROUTINE QDYNAM
      return
      END SUBROUTINE QDYNAM

#ifdef CACHED_SUBDD
      subroutine fijlh_defs(arr,nmax,decl_count)
c
c 3D outputs
c
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      use model_com, only: dtsrc
      use constant, only: kapa
      use TimeConstants_mod, only: SECONDS_PER_DAY
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c
      decl_count = 0
c
      arr(next()) = info_type_(
     &  sname = 'dq_ls',
     &  lname = 'moisture tendency from large-scale forcings',
     &  units = 'kg/kg/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dth_ls',
     &  lname = 'theta tendency from large-scale forcings',
     &  units = 'K/day',
     &  scale = 1000.**kapa/dtsrc*SECONDS_PER_DAY
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'du_ls',
     &  lname = 'zonal wind tendency from large-scale forcings',
     &  units = 'm/s/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dv_ls',
     &  lname = 'meridional wind tendency from large-scale forcings',
     &  units = 'm/s/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dq_nudge',
     &  lname = 'moisture tendency from nudging',
     &  units = 'kg/kg/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dth_nudge',
     &  lname = 'theta tendency from nudging',
     &  units = 'K/day',
     &  scale = 1000.**kapa/dtsrc*SECONDS_PER_DAY
     &     )
c
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine fijlh_defs
#endif
