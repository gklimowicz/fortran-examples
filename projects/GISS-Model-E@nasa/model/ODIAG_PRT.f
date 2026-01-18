#include "rundeck_opts.h"

      subroutine diag_ocean_prep
      use oceanr_dim, only : grid=>ogrid
      implicit none
      if(.not. grid%have_domain) return
      call basin_prep
      call oijl_prep
      return
      end subroutine diag_ocean_prep

      SUBROUTINE diag_OCEAN
      END SUBROUTINE diag_OCEAN

      SUBROUTINE STRMIJ (MFU,FAC,OLNST,FACST,SF)
!@sum STRMIJ calculates the latitude by longitude stream function for a
!@+   given quantity.
C****
C**** Input:
!@var MFU  = west-east and south-north tracer fluxes (kg/s)
!@var   FAC = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      Use OCEAN, Only: IM,JM,LMO, oDLAT_DG,oFJEQ=>FJEQ
      USE STRAITS, only : nmst,lmst
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN), DIMENSION(IM,JM) :: MFU
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(IM,JM) :: SF
      Integer*4 :: I,J,L, NSUM
      REAL*8 TSUM
C****
C**** Integrate up from South Pole
C****
      SF=0
      DO J=2,JM-1
        DO I=1,IM
          SF(I,J) = SF(I,J-1) + MFU(I,J)*FAC
        END DO
C****
C**** Add strait flow from to the Stream Function
C****
        CALL STRMIJ_STRAITS(J,SF,OLNST,FACST)
      EndDo  !  End of Do J=2,JM-1

C****
C**** Recalibrate SF to be 0 over middle of North America
C**** Include UO cells whose centers reside in 110-90 W, 32-40 N
C****
      NSUM = 0
      TSUM = 0
      Do J = Ceiling(oFJEQ+32/oDLAT_DG), Floor(oFJEQ+40/oDLAT_DG)
      Do I = Ceiling(70*IM/360.), Floor(90*IM/360.)
        TSUM = TSUM + SF(I,J)
        NSUM = NSUM + 1  ;  EndDo  ;  EndDo
      SF(:,:) = SF(:,:) - TSUM/NSUM
      RETURN
      END

      subroutine oijl_prep
c
c Convert oijl accumulations into the desired units
c
      use ocean, only : im,jm,lmo,lmm,imaxj,focean,dxypo,dxvo,dypo,dts
      use ocean, only : nbyzm,i1yzm,i2yzm
      use ocnmeso_com, only : use_tdmix
      use odiag, only : koijl,oijl_out,oijl=>oijl_loc,ijl_area
     &     ,igrid_oijl,jgrid_oijl,sname_oijl
     &     ,ijl_mo,ijl_mou,ijl_mov,ijl_g0m,ijl_s0m,ijl_ptm,ijl_pdm
     &     ,ijl_mfu,ijl_mfv,ijl_mfw,ijl_mfw2,ijl_ggmfl,ijl_sgmfl
#ifdef TDMIX_AUX_DIAGS
     &     ,ijl_gsymmf,ijl_ssymmf
#endif
     &     ,ijl_wgfl,ijl_wsfl,ijl_kvm,ijl_kvg,ijl_kvx,ijl_gflx,ijl_sflx
     &     ,ijl_mfub,ijl_mfvb,ijl_mfwb,ijl_isdm,ijl_pdm2
     &     ,oij=>oij_loc,ij_sf,olnst,ln_mflx
#ifdef OCN_GISS_TURB
     &     ,ijl_ri,ijl_rrho,ijl_bv2,ijl_otke,ijl_kvs,ijl_kvc,ijl_buoy
#endif
#ifdef OCN_GISS_SM
     &     ,ijl_fvb
#endif
#ifdef OCEAN_TENDENCY_DIAGS
     &     ,olnst,ln_gflx,ln_sflx
#endif
      use odiag, only : ia_oijl
#ifdef TRACERS_OCEAN
      use odiag, only :
     &     ktoijlx,toijl_out,divbya_toijl,kn_toijl,toijl_loc,toijl_conc
     &    ,toijl_tflx,toijl_gmfl
      USE OCN_TRACER_COM, only : n_Water, tracerlist, ocn_tracer_entry
#endif
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : am_i_root,halo_update,south
     &     ,hassouthpole,hasnorthpole
     &     ,pack_data,unpack_data ! for horz stream function
      use mdiag_com, only : ia_cpl
      use model_com, only : idacc
      use constant, only : grav
      use kpp_com, only : use_tdiss
#ifdef OCEAN_TENDENCY_DIAGS
      use straits, only: ist,jst,lmst
      use domain_decomp_1d, only : broadcast
#endif
      use straits, only : nmst
      implicit none
      integer i,j,l,k,kk,n
      real*8 mass,gos,sos,temgs,volgs,volgsp,fac,facst,dpr
      real*8 wdenom,xedge,wtdn,wtup,massdn,massup
      integer :: j_0,j_1,j_0s,j_1s
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     mfu,pres
      real*8, dimension(:,:), allocatable :: mfu_glob,sf_glob
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mfub,mfvb
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,0:lmo) ::
     &     mfwb,dmfwb
      real*8 :: sncor
      integer :: ib
#ifdef TRACERS_OCEAN
      type(ocn_tracer_entry), pointer :: entry
#endif
#ifdef OCEAN_TENDENCY_DIAGS
      integer, parameter :: max_num_tends=16 ! maximum number of tendency outputs
      integer :: flxind(max_num_tends)
      character(len=32) :: tendname(max_num_tends)
      real*8 :: sign_end
#endif
      real*8 dumarr(1,1)

      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0s = grid%j_strt_skp
      j_1s = grid%j_stop_skp

      oijl_out(:,:,:,:) = 0.

c
c Cell-centered quantities. Some conversions to per square meter
c
      pres = 0. ! need to add atm + sea ice press, but this is a small
                ! effect for typical uses of in-situ density diag
      do l=1,lmo
      do j=j_0,j_1
      do i=1,imaxj(j)
        mass = oijl(i,j,l,ijl_mo)*dxypo(j)
        if(focean(i,j).le..5 .or. mass.le.0.) cycle
        oijl_out(i,j,l,ijl_mo) = mass
        oijl_out(i,j,l,ijl_g0m) = oijl(i,j,l,ijl_g0m)
        oijl_out(i,j,l,ijl_s0m) = oijl(i,j,l,ijl_s0m)
c
c compute potential temperature, potential density, in-situ density
c
        dpr = grav*oijl(i,j,l,ijl_mo)/idacc(ia_oijl(ijl_mo))
        pres(i,j) = pres(i,j) + .5d0*dpr
        gos = oijl(i,j,l,ijl_g0m) / mass
        sos = oijl(i,j,l,ijl_s0m) / mass
        oijl_out(i,j,l,ijl_ptm) = mass*temgs(gos,sos)
        oijl_out(i,j,l,ijl_pdm) = mass*(1d0/volgs(gos,sos)-1000d0) !sigma0
        oijl_out(i,j,l,ijl_pdm2)= mass*(1d0/volgsp(gos,sos,2d7)-1d3) !sigma2
        oijl_out(i,j,l,ijl_isdm) =
     &       mass*(1d0/volgsp(gos,sos,pres(i,j))-1000d0)
        pres(i,j) = pres(i,j) + .5d0*dpr
      enddo
      enddo
      enddo


      dmfwb = 0.
      mfub = oijl(:,:,:,ijl_mfub)
      mfvb = oijl(:,:,:,ijl_mfvb)
      mfwb = 0.
      call halo_update(grid,mfvb,from=south)
      if(use_tdmix==1) then
        ! Calculate dmfwb, the bolus-induced component of the mass
        ! flux used in the remapping (vertical advective) tracer flux.
        ! It is used below for post-hoc repartitioning of accumulated
        ! advective fluxes into resolved- and bolus-velocity components.
        ! See notes in TDMIX.
        mfwb(:,:,1:lmo) = oijl(:,:,1:lmo,ijl_mfwb)
        do l=1,lmo-1
        do j=j_0s,j_1s
          i=1
          if(l.lt.lmm(i,j)) then
            dmfwb(i,j,l) = dmfwb(i,j,l-1) + (
     &         (mfub(im,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
     &         )!/dxypo(j)
          endif
          do n=1,nbyzm(j,l+1)
          do i=max(2,i1yzm(n,j,l+1)),i2yzm(n,j,l+1)
            dmfwb(i,j,l) = dmfwb(i,j,l-1) + (
     &         (mfub(i-1,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
     &         )!/dxypo(j)
          enddo
          enddo
        enddo
        enddo
      else
c****
c**** derive bolus vertical mass flux from bolus horizontal mass fluxes
c**** for skew-GM
        do l=1,lmo-1
        do j=j_0s,j_1s
          i=1
          if(l.lt.lmm(i,j)) then
            mfwb(i,j,l) = mfwb(i,j,l-1) + (
     &         (mfub(im,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &         )
          endif
          do n=1,nbyzm(j,l+1)
          do i=max(2,i1yzm(n,j,l+1)),i2yzm(n,j,l+1)
            mfwb(i,j,l) = mfwb(i,j,l-1) + (
     &         (mfub(i-1,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &         )
          enddo
          enddo
        enddo
        enddo
        oijl(:,:,:,ijl_mfwb) = mfwb(:,:,1:lmo)
      endif
c
c Vertical fluxes.  Some conversions to per square meter
c
      do l=1,lmo-1
      do j=j_0,j_1
      do n=1,nbyzm(j,l+1)
      do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
        oijl_out(i,j,l,ijl_area) = idacc(ia_cpl)*dxypo(j)
          
        oijl_out(i,j,l,ijl_mfw) = oijl(i,j,l,ijl_mfw)
     &       -dmfwb(i,j,l) ! subtract bolus-induced part of remap flux
        oijl_out(i,j,l,ijl_mfwb) = oijl(i,j,l,ijl_mfwb)
     &       +dmfwb(i,j,l) !      add bolus-induced part of remap flux


        massup = oijl(i,j,l  ,ijl_mo)*dxypo(j)
        massdn = oijl(i,j,l+1,ijl_mo)*dxypo(j)
        wtdn = massup/(massup+massdn)
        wtup = 1d0-wtdn

        ! centered approximation of layer edge g for repartitioning
        xedge =
     &       wtup* (oijl(i,j,l  ,ijl_g0m) / massup)
     &      +wtdn* (oijl(i,j,l+1,ijl_g0m) / massdn)
        oijl_out(i,j,l,ijl_gflx+2) = oijl(i,j,l,ijl_gflx+2)
     &       -dmfwb(i,j,l)*xedge ! subtract bolus-induced part of remap flux
        oijl_out(i,j,l,ijl_ggmfl+2) = oijl(i,j,l,ijl_ggmfl+2)
     &       +dmfwb(i,j,l)*xedge !      add bolus-induced part of remap flux

        ! centered approximation of layer edge s for repartitioning
        xedge =
     &       wtup* (oijl(i,j,l  ,ijl_s0m) / massup)
     &      +wtdn* (oijl(i,j,l+1,ijl_s0m) / massdn)
        oijl_out(i,j,l,ijl_sflx+2) = oijl(i,j,l,ijl_sflx+2)
     &       -dmfwb(i,j,l)*xedge ! subtract bolus-induced part of remap flux
        oijl_out(i,j,l,ijl_sgmfl+2) = oijl(i,j,l,ijl_sgmfl+2)
     &       +dmfwb(i,j,l)*xedge !      add bolus-induced part of remap flux

        oijl_out(i,j,l,ijl_mfw2) = oijl(i,j,l,ijl_mfw2)/dxypo(j)
#ifdef TDMIX_AUX_DIAGS
        oijl_out(i,j,l,ijl_gsymmf+2) = oijl(i,j,l,ijl_gsymmf+2)
        oijl_out(i,j,l,ijl_ssymmf+2) = oijl(i,j,l,ijl_ssymmf+2)
#endif

        oijl_out(i,j,l,ijl_wgfl) = oijl(i,j,l,ijl_wgfl)
        oijl_out(i,j,l,ijl_wsfl) = oijl(i,j,l,ijl_wsfl)
        oijl_out(i,j,l,ijl_kvm) = oijl(i,j,l,ijl_kvm)*dxypo(j)
        oijl_out(i,j,l,ijl_kvg) = oijl(i,j,l,ijl_kvg)*dxypo(j)
        if(use_tdiss==1) then
          oijl_out(i,j,l,ijl_kvx) = oijl(i,j,l,ijl_kvx)*dxypo(j)
        endif
#ifdef OCN_GISS_TURB
        oijl_out(i,j,l,ijl_kvs) = oijl(i,j,l,ijl_kvs)*dxypo(j)
        oijl_out(i,j,l,ijl_kvc) = oijl(i,j,l,ijl_kvc)*dxypo(j)
        oijl_out(i,j,l,ijl_ri) = oijl(i,j,l,ijl_ri)*dxypo(j)
        oijl_out(i,j,l,ijl_rrho) = oijl(i,j,l,ijl_rrho)*dxypo(j)
        oijl_out(i,j,l,ijl_bv2) = oijl(i,j,l,ijl_bv2)*dxypo(j)
        oijl_out(i,j,l,ijl_buoy) = oijl(i,j,l,ijl_buoy)*dxypo(j)
        oijl_out(i,j,l,ijl_otke) = oijl(i,j,l,ijl_otke)*dxypo(j)
#endif
#ifdef OCN_GISS_SM
        oijl_out(i,j,l,ijl_fvb) = oijl(i,j,l,ijl_fvb)*dxypo(j)
#endif
      enddo
      enddo
      enddo
      enddo

c
c Horizontal fluxes.  Some conversions to per meter
c
      do k=1,koijl
        call halo_update(grid,oijl(:,:,:,k),from=south)
      enddo
      do l=1,lmo
        do j=j_0s,j_1s
        do i=1,im
          oijl_out(i,j,l,ijl_mfu) = oijl(i,j,l,ijl_mfu)
          oijl_out(i,j,l,ijl_mfub) = oijl(i,j,l,ijl_mfub)
          oijl_out(i,j,l,ijl_gflx) = oijl(i,j,l,ijl_gflx)
          oijl_out(i,j,l,ijl_sflx) = oijl(i,j,l,ijl_sflx)
          oijl_out(i,j,l,ijl_ggmfl) = oijl(i,j,l,ijl_ggmfl)
          oijl_out(i,j,l,ijl_sgmfl) = oijl(i,j,l,ijl_sgmfl)
#ifdef TDMIX_AUX_DIAGS
          oijl_out(i,j,l,ijl_gsymmf) = oijl(i,j,l,ijl_gsymmf)
          oijl_out(i,j,l,ijl_ssymmf) = oijl(i,j,l,ijl_ssymmf)
#endif
        enddo
        do i=1,im-1
          oijl_out(i,j,l,ijl_mou) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(i+1,j,l,ijl_mo))*dypo(j)
        enddo
        i=im
          oijl_out(i,j,l,ijl_mou) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(1,j,l,ijl_mo))*dypo(j)
        enddo ! j
        do j=max(2,j_0),j_1
        do i=1,im
          oijl_out(i,j,l,ijl_mfv) = oijl(i,j-1,l,ijl_mfv)
          oijl_out(i,j,l,ijl_mfvb) = oijl(i,j-1,l,ijl_mfvb)
          oijl_out(i,j,l,ijl_mov) =
     &         .5*(oijl(i,j,l,ijl_mo)+oijl(i,j-1,l,ijl_mo))*dxvo(j-1)
          oijl_out(i,j,l,ijl_gflx+1) = oijl(i,j-1,l,ijl_gflx+1)
          oijl_out(i,j,l,ijl_sflx+1) = oijl(i,j-1,l,ijl_sflx+1)
          oijl_out(i,j,l,ijl_ggmfl+1) = oijl(i,j-1,l,ijl_ggmfl+1)
          oijl_out(i,j,l,ijl_sgmfl+1) = oijl(i,j-1,l,ijl_sgmfl+1)
#ifdef TDMIX_AUX_DIAGS
          oijl_out(i,j,l,ijl_gsymmf+1) = oijl(i,j-1,l,ijl_gsymmf+1)
          oijl_out(i,j,l,ijl_ssymmf+1) = oijl(i,j-1,l,ijl_ssymmf+1)
#endif
        enddo
        enddo ! j
      enddo

#ifdef OCEAN_TENDENCY_DIAGS
      ! compute convergences of 3D resolved/parameterized fluxes
      mfwb(:,:,0) = 0.

      k = 0
      k = k+1; flxind(k) = ijl_gflx; tendname(k) = 'g_advtend'
      k = k+1; flxind(k) = ijl_sflx; tendname(k) = 's_advtend'
      k = k+1; flxind(k) = ijl_ggmfl; tendname(k) = 'g_mesotend'
      k = k+1; flxind(k) = ijl_sgmfl; tendname(k) = 's_mesotend'
#ifdef TDMIX_AUX_DIAGS
      k = k+1; flxind(k) = ijl_gsymmf; tendname(k) = 'g_mesotend_sym'
      k = k+1; flxind(k) = ijl_ssymmf; tendname(k) = 's_mesotend_sym'
#endif
      kk = k
      ! currently, convergences are in assumed positions in oijl (xflux+3)
      do k=1,kk ! so sanity-check the names first
        if(trim(sname_oijl(flxind(k)+3)).ne.trim(tendname(k))) then
          call stop_model(
     &         'oijl_prep: OCEAN_TENDENCY_DIAGS name mismatch',255)
        endif
      enddo
      do k=1,kk
        mfwb(:,:,1:lmo) = oijl_out(:,:,1:lmo,flxind(k)+2)
        call do_fluxconv_3d(
     &       oijl(:,:,:,flxind(k)+0),
     &       oijl(:,:,:,flxind(k)+1),
     &       mfwb,
     &       oijl_out(:,:,:,flxind(k)+3) )
      enddo

      ! Deal with the straits
      call broadcast(grid, olnst)
      do n=1,nmst
        sign_end = -.5 ! .5 is the scale factor for olnst
        do k=1,2
          i=ist(n,k)
          j=jst(n,k)
          if(j.ge.j_0 .and. j.le.j_1) then
            do l=1,lmst(n)
              oijl_out(i,j,l,ijl_gflx+3) = oijl_out(i,j,l,ijl_gflx+3)
     &             + olnst(l,n,ln_gflx)*sign_end
              oijl_out(i,j,l,ijl_sflx+3) = oijl_out(i,j,l,ijl_sflx+3)
     &             + olnst(l,n,ln_sflx)*sign_end
            enddo
          endif
          sign_end = -sign_end
        enddo
      enddo
#endif

      ! Fill poles
      do k=1,koijl
        if(jgrid_oijl(k).ne.1 .or. igrid_oijl(k).ne.1) cycle
        do l=1,lmo
          !if(hassouthpole(grid)) then
          !  j = j_0
          !  oijl_out(2:im,j,l,k) = oijl_out(1,j,l,k)
          !endif
          if(hasnorthpole(grid)) then
            j = j_1
            oijl_out(2:im,j,l,k) = oijl_out(1,j,l,k)
          endif
        enddo
      enddo

C****
C**** Calculate Horizontal Mass Stream Function
C****
      do j=j_0s,j_1s
      do i=1,im
        mfu(i,j) = 0.
        do l=1,lmm(i,j)
          mfu(i,j) = mfu(i,j) + oijl(i,j,l,ijl_mfu)
        enddo
      enddo
      enddo
      if(am_i_root()) then
        allocate(mfu_glob(im,jm),sf_glob(im,jm))
      else
        allocate(mfu_glob(1,1),sf_glob(1,1))
      endif
      call pack_data(grid,mfu,mfu_glob)
      if(am_i_root()) then
        FAC   = -1d-9/dts
        FACST = -1d-9/dts
        if(nmst.gt.0) then
          CALL STRMIJ(MFU_GLOB,FAC,OLNST(1,1,LN_MFLX),FACST,SF_GLOB)
        else
          CALL STRMIJ(MFU_GLOB,FAC,DUMARR,FACST,SF_GLOB)
        endif
      endif
      call unpack_data(grid,sf_glob,oij(:,:,ij_sf))
      deallocate(mfu_glob,sf_glob)

#ifdef TRACERS_OCEAN
C****
C**** Tracers
C****
      toijl_out(:,:,:,:) = 0.
      kk = 1
      toijl_out(:,:,:,kk) = oijl_out(:,:,:,ijl_mo)
      do kk=2,ktoijlx
        k = kn_toijl(1,kk)
        n = kn_toijl(2,kk)
        entry=>tracerlist%at(n)
        if(k.le.0 .or. n.le.0) cycle
        toijl_out(:,:,:,kk) = toijl_loc(:,:,:,k,n)
        if(divbya_toijl(kk)) then
          do l=1,lmo; do j=j_0,j_1
            toijl_out(:,j,l,kk) = toijl_out(:,j,l,kk)/dxypo(j)
          enddo; enddo
        endif
        if(entry%to_per_mil>0 .and. n.ne.n_Water) then
          toijl_out(:,:,:,kk) = 1d3*(toijl_out(:,:,:,kk)/entry%trw0
     &         -toijl_loc(:,:,:,TOIJL_conc,n_water))
        endif
        if(use_tdmix==1 .and.
     &       (k.eq.toijl_gmfl+2 .or. k.eq.toijl_tflx+2)) then
          if(k.eq.toijl_gmfl+2) then
            sncor = +1d0
          else
            sncor = -1d0
          endif
          do l=1,lmo-1
          do j=j_0,j_1
          do ib=1,nbyzm(j,l+1)
          do i=i1yzm(ib,j,l+1),i2yzm(ib,j,l+1)
            massup = oijl(i,j,l  ,ijl_mo)*dxypo(j)
            massdn = oijl(i,j,l+1,ijl_mo)*dxypo(j)
            wtdn = massup/(massup+massdn)
            wtup = 1d0-wtdn
            ! centered approximation of layer edge tracer for repartitioning
            xedge =
     &           wtup* (toijl_loc(i,j,l  ,toijl_conc,n) / massup)
     &          +wtdn* (toijl_loc(i,j,l+1,toijl_conc,n) / massdn)
            toijl_out(i,j,l,kk) = toijl_out(i,j,l,kk) 
                 ! add or subtract bolus-induced part of remap flux
     &           +(dmfwb(i,j,l)/dxypo(j))*xedge*sncor
          enddo
          enddo
          enddo
          enddo
        endif
      enddo
#endif

      return

      contains

      subroutine do_fluxconv_3d(mfub,mfvb,mfwb,conv_out)
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     mfub,mfvb,conv_out
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,0:lmo) ::
     &     mfwb

      real*8 :: byim

      !call halo_update(grid,mfvb,from=south)
      do l=1,lmo
        do j=grid%j_strt_skp,grid%j_stop_skp
          i=1
          if(l.le.lmm(i,j)) then
            conv_out(i,j,l) =
     &         (mfub(im,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
          endif
          do n=1,nbyzm(j,l)
          do i=max(2,i1yzm(n,j,l)),i2yzm(n,j,l)
            conv_out(i,j,l) =
     &         (mfub(i-1,j,l)-mfub(i,j,l))
     &        +(mfvb(i,j-1,l)-mfvb(i,j,l))
     &        +(mfwb(i,j,l-1)-mfwb(i,j,l))
          enddo
          enddo
        enddo
      enddo
      if(hasnorthpole(grid)) then
        j = grid%j_stop
        byim = 1d0/real(im,kind=8)
        do l=1,lmm(1,j)
          conv_out(:,j,l) = sum(mfvb(:,j-1,l))*byim
     &         +(mfwb(1,j,l-1)-mfwb(1,j,l))
        enddo
      endif
      end subroutine do_fluxconv_3d

      end subroutine oijl_prep

      Subroutine STRMIJ_STRAITS (J,SF,OLNST,FACST)
C****
C**** Add strait flow to IxJ Strean Function
C****
C**** Input:   J = ocean model latitude index
C****      OLNST = ocean strait mass flux (kg/s)
C****      FACST = global scaling factor for ocean straits
C**** Output: SF = IxJ stream function (kg/s)
C****
      Use OCEAN,   Only: IM,JM,LMO
      Use STRAITS, Only: NMST,LMST, IST,JST
      Implicit None

      Integer*4,Intent(In) :: J
      Real*8,Intent(InOut) :: SF(IM,JM)
      Real*8,Intent(In)    :: OLNST(LMO,NMST),FACST

C**** Local variables
      Integer*4 N,LM, I1,J1, I2,J2

      Do 40 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      LM = LMST(N)
      If (J2 - J1) 10,20,30

C**** JST(N,2) < JST(N,1)
   10 If (J==J1 .or. J==J2)  Then
         SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST*.5/(J1-J2)
         GoTo 40  ;  EndIf
      If (J2 < J .and. J < J1)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST/(J1-J2)
      GoTo 40

C**** JST(N,2) = JST(N,1)
   20 If (J==J1)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) + Sum(OLNST(1:LM,N))*FACST
      GoTo 40

C**** JST(N,2) > JST(N,1)
   30 If (J==J1 .or. J==J2)  Then
         SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST*.5/(J2-J1)
         GoTo 40  ;  EndIf
      If (J1 < J .and. J < J2)
     *   SF(I1:I2-1,J) = SF(I1:I2-1,J) +
     *                   Sum(OLNST(1:LM,N))*FACST/(J2-J1)
   40 Continue
      Return
      EndSubroutine STRMIJ_STRAITS
