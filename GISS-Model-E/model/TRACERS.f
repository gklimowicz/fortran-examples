#include "rundeck_opts.h"

#ifdef TRACERS_ON
!@sum  TRACERS: generic tracer routines used for all tracers
!@+    Routines included:
!@+      Generic diags: set_generic_tracer_diags
!@+      Apply previously set sources: apply_tracer_sources
!@+      Radioactive Decay: tdecay
!@+      Gravitaional Settling: trgrav
!@+      Check routine: checktr
!@auth Jean Lerner/Gavin Schmidt
      MODULE apply3d
!@sum apply3d is used simply so that I can get optional arguments
!@+   to work. If anyone can some up with something neater, let me know.
      use TracerSource_mod, only: TracerSource3D
      use Tracer_mod, only: Tracer
      use OldTracer_mod, only: trname
      USE TRACER_COM, only : NTM,trm,trmom,alter_sources,tracers
      USE CONSTANT, only : teeny
      USE RESOLUTION, only: lm
      USE MODEL_COM, only : dtsrc
      USE GEOM, only : imaxj,byaxyp,lat2d_dg,lon2d_dg
      USE QUSDEF, only: nmom
#ifndef SKIP_TRACER_SRCS
      USE FLUXES, only : tr3Dsource
#endif
      USE TRDIAG_COM, only : jls_3Dsource,itcon_3Dsrc
     *     ,ijts_3Dsource,taijs=>taijs_loc
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds, am_i_root
      use EmissionRegion_mod, only: numRegions, regions

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE apply_tracer_3Dsource( ns , n , momlog )
      USE CONSTANT, only : UNDEF_VAL
!@sum apply_tracer_3Dsource adds 3D sources to tracers
!@auth Jean Lerner/Gavin Schmidt
!@var MOM true (default) if moments are to be modified
      logical, optional, intent(in) :: momlog
      integer, intent(in) :: n,ns
      real*8 fred(grid%i_strt:grid%i_stop),
     &     dtrm(grid%i_strt_halo:grid%i_stop_halo,
     &          grid%j_strt_halo:grid%j_stop_halo,lm),
     &     dtrml(lm),eps

      logical :: domom
      integer najl,i,j,l,naij,kreg,nsect,nn
      INTEGER :: J_0, J_1, I_0, I_1

      type (TracerSource3D), pointer :: source
      class (Tracer), pointer :: pTracer

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Ensure that this is a valid tracer and source
      if (n.eq.0 .or. ns.eq.0) then
         return
      end if
#ifndef SKIP_TRACER_SRCS
C**** parse options
      domom=.true.
      if( present(momlog) ) then
        domom=momlog
      end if
C**** This is tracer independent coding designed to work for all
C**** 3D sources.
C**** Modify tracer amount, moments, and diagnostics
      najl = jls_3Dsource(ns,n)
      naij = ijts_3Dsource(ns,n)

      eps = tiny(trm(i_0,j_0,1,n))
      fred = UNDEF_VAL
      dtrm = UNDEF_VAL
      do l=1,lm
      do j=j_0,j_1
        do i=i_0,imaxj(j)
          dtrm(i,j,l) = tr3Dsource(i,j,l,ns,n)*dtsrc
C**** calculate fractional loss and update tracer mass

          fred(i) = max(0.d0,
     &         1.+min(0.d0,dtrm(i,j,l))/(trm(i,j,l,n)+eps))

          trm(i,j,l,n) = trm(i,j,l,n)+dtrm(i,j,l)
          if(fred(i).le.1d-16) trm(i,j,l,n) = 0.
        end do
        if(domom .and. any(fred.lt.1.)) then
          do i=i_0,imaxj(j)
            trmom(:,i,j,l,n) = trmom(:,i,j,l,n)*fred(i)
          enddo
        endif
      end do
      if (naij.gt.0) then
        do j=j_0,j_1
          do i=i_0,imaxj(j)
            taijs(i,j,naij) = taijs(i,j,naij) + dtrm(i,j,l)
          end do
        end do
      end if
      end do ! l

      if(jls_3Dsource(ns,n) > 0) then
        do j=j_0,j_1
          do i=i_0,imaxj(j)
            dtrml(:) = dtrm(i,j,:)
            call inc_tajls_column(i,j,1,lm,lm,najl,dtrml)
          enddo
        enddo
      endif
#endif

      if (itcon_3Dsrc(ns,n).gt.0)
     *  call DIAGTCA(itcon_3Dsrc(ns,n),n)

C****
      RETURN
      END SUBROUTINE apply_tracer_3Dsource

      end module apply3d

      SUBROUTINE set_generic_tracer_diags
!@sum set_generic_tracer_diags init trace gas attributes and diagnostics
!@auth J. Lerner
!@calls sync_param
      use OldTracer_mod, only: trName, tr_wd_TYPE, mass2vol
      use OldTracer_mod, only: dowetdep, dodrydep, trradius, ntrocn
      use OldTracer_mod, only: ntm_power, nwater, src_dist_index
      USE CONSTANT, only: mair
      USE MODEL_COM, only: dtsrc
      USE FLUXES, only : nisurf,atmice
      USE DIAG_COM, only: ia_src,ia_12hr,ir_log2,ir_0_71
      USE TRACER_COM, only: ntm, n_Water, nPart
      USE TRDIAG_COM
      use Dictionary_mod
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      implicit none
      integer :: l,k,n,nd
      character*10, DIMENSION(NTM) :: CMR,CMRWT
      logical :: T=.TRUE. , F=.FALSE.
      character*50 :: unit_string

      atmice%taijn => taijn_loc
#ifdef TRACERS_WATER
      allocate(atmice%do_accum(ntm))
      do n=1,ntm
        atmice%do_accum(n) =
     &       (tr_wd_TYPE(n).eq.nWater .or. tr_wd_TYPE(n).eq.nPART)
      enddo
#endif

C**** Get factor to convert from mass to volume mr if necessary
      do n=1,ntm
        if (to_volume_MixRat(n) .eq.1) then
          MMR_to_VMR(n) = mass2vol(n)
          cmr(n) = 'V/V air'
        else
          MMR_to_VMR(n) = 1.d0
          cmr(n) = 'kg/kg air'
        endif
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          cmrwt(n) = 'per mil'
          cmr(n) = 'per mil'     ! this overrides to_volume_ratio
          MMR_to_VMR(n) = 1.
        else
          cmrwt(n) = 'kg/m^2'
        end if
#endif
      end do
C****
C**** TAJLN(J,L,KQ,N)  (SUM OVER LONGITUDE AND TIME OF)
C****
C**** jlq_power Exponent associated with a physical process
C****      (for printing tracers).   (ntm_power+jlq_power)
C**** jls_power Exponent associated with a source/sink (for printing)
C**** Defaults for JLN
      scale_jlq(:) = 1./DTsrc
      jgrid_jlq(:) = 1
      ia_jlq(:) = ia_src

C**** Tracer concentration
      do n=1,ntm
        k = 1        ! <<<<< Be sure to do this
        jlnt_conc = k
        sname_jln(k,n) = trim(trname(n))//'_CONCENTRATION'
        lname_jln(k,n) = trim(trname(n))//' CONCENTRATION'
        jlq_power(k) = 0.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),cmr(n))
        scale_jlq(k) = 1.d0
        scale_jln(n) = MMR_to_VMR(n)
#ifdef TRACERS_WATER
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif
C**** Tracer mass
        k = k + 1
        jlnt_mass = k
        sname_jln(k,n) = trim(trname(n))//'_MASS'
        lname_jln(k,n) = trim(trname(n))//' MASS'
        jlq_power(k) = 4

#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) ||\
    (defined TRACERS_TOMAS) || (defined TRACERS_AEROSOLS_SEASALT)
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)+13
     *       ,'kg')
#else
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)
     *       ,'kg/m^2')
#endif
        scale_jlq(k) = 1.d0
#ifdef TRACERS_WATER
C****   TRACER CONCENTRATION IN CLOUD WATER
        k = k + 1
        jlnt_cldh2o = k
        sname_jln(k,n) = trim(trname(n))//'_WM_CONC'
        lname_jln(k,n) = trim(trname(n))//' CLOUD WATER CONCENTRATION'
        jlq_power(k) = 4
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k)
     *       ,'kg/kg water')
        scale_jlq(k) = 1.d0
        if (to_per_mil(n).gt.0) units_jln(k,n) = unit_string(0,cmr(n))
#endif
C**** Physical processes affecting tracers
C****   F (TOTAL NORTHWARD TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_nt_tot = k
        sname_jln(k,n) = 'tr_nt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL NORTHWARD TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11
        jgrid_jlq(k) = 2
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   STM/SM (MEAN MERIDIONAL N.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_nt_mm = k
        sname_jln(k,n) = 'tr_nt_mm_'//trname(n)
        lname_jln(k,n) = 'NORTHWARD TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10
        jgrid_jlq(k) = 2
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   F (TOTAL VERTICAL TRANSPORT OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_tot = k
        sname_jln(k,n) = 'tr_vt_tot_'//trname(n)
        lname_jln(k,n) = 'TOTAL VERTICAL TRANSPORT OF '//
     &     trim(trname(n))//' MASS'
        jlq_power(k) = 11
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   STM/SM (MEAN MERIDIONAL V.T. OF TRACER MASS)  (kg)
        k = k + 1
        jlnt_vt_mm = k
        sname_jln(k,n) = 'tr_vt_mm_'//trname(n)
        lname_jln(k,n) = 'VERTICAL TRANS. OF '//
     &     trim(trname(n))//' MASS BY MERIDIONAL CIRC.'
        jlq_power(k) = 10.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY MOIST CONVEC)(kg)
        k = k + 1
        jlnt_mc = k
        sname_jln(k,n) = 'tr_mc_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY MOIST CONVECTION'
        jlq_power(k) = 10
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY Large-scale CONDENSE)  (kg)
        k = k + 1
        jlnt_lscond = k
        sname_jln(k,n) = 'tr_lscond'//trname(n)
        lname_jln(k,n) ='CHANGE OF '//
     &     trim(trname(n))//' MASS BY LARGE-SCALE CONDENSE'
        jlq_power(k) = 10.
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')
C****   TMBAR-TM (CHANGE OF TRACER MASS BY DRY CONVEC)  (kg)
        k = k + 1
        jlnt_turb = k
        sname_jln(k,n) = 'tr_turb_'//trname(n)
        lname_jln(k,n) = 'CHANGE OF '//
     &     trim(trname(n))//' MASS BY TURBULENCE/DRY CONVECTION'
        jlq_power(k) = 10
        units_jln(k,n) = unit_string(ntm_power(n)+jlq_power(k),'kg/s')


        if (k.gt. ktajl) then
           if (AM_I_ROOT()) write (6,*)
     &          'tjl_defs: Increase ktajl=',ktajl,' to at least ',k
           call stop_model('ktajl too small',255)
        end if
      end do

C**** CONTENTS OF TAIJLN(I,J,LM,N)  (SUM OVER TIME OF)
C****        TML (M*M * KG TRACER/KG AIR)
C**** Set defaults that are true for all tracers and layers
      do n=1,ntm
        ijtm_power(n) = ntm_power(n)+4 !n for integrated mass
        ijtc_power(n) = ntm_power(n)+1 !n for concentration
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) ijtc_power(n) = 0
#endif
      end do
C**** Tracer concentrations (TAIJLN)
      do n=1,ntm
        write(sname_ijt(n),'(a)') trim(TRNAME(n))
        write(lname_ijt(n),'(a)') trim(TRNAME(n))
        if (to_conc(n).eq.1) then   ! diag in kg/m3
          units_ijt(n) = unit_string(ijtc_power(n),'kg/m3')
          scale_ijt(n) = 10.**(-ijtc_power(n))
        else ! mixing ratio
          units_ijt(n) = unit_string(ijtc_power(n),cmr(n))
          scale_ijt(n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
        end if
      end do

C**** AIJN
C****     1  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER/KG AIR)
C****     2  TRS (SURFACE TRACER CONC.) (M*M * KG TRACER/KG AIR)
C****     3  TM (SUM OVER ALL LAYERS) (M*M * KG TRACER)

      do n=1,ntm
C**** Summation of mass over all layers
      k = 0        ! <<<<< Be sure to do this
      if (src_dist_index(n)<=1) then
        k = k+1
        tij_mass = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Total_Mass'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Total Mass'
        units_tij(k,n) = unit_string(ijtm_power(n),'kg/m^2')
        scale_tij(k,n) = 10.**(-ijtm_power(n))
C**** Average concentration over layers
        k = k+1
        tij_conc = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_Average'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' Average'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n) = MMR_to_VMR(n)*10.**(-ijtc_power(n))
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        endif
#endif
C**** Surface concentration
        k = k+1
        tij_surf = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//' At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),cmr(n))
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
#ifdef TRACERS_WATER
        if (to_per_mil(n) .eq.1) then
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        endif
#endif
C**** Surface concentration by volume (units kg/m^3)
        k = k+1
        tij_surfbv = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *     '_byVol_At_Surface'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' By Vol At Surface'
        units_tij(k,n) = unit_string(ijtc_power(n),'kg/m^3')
        scale_tij(k,n)=MMR_to_VMR(n)*10.**(-ijtc_power(n))/
     *                 REAL(NIsurf,KIND=8)
      endif ! if (src_dist_index(n)<=1) then
C**** Tropopause flux Diagnostics
        k = k+1
        tij_strop = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *     '_StratTropflux'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Flux Tropopause'
        units_tij(k,n) = unit_string(ijtc_power(n),'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n))
#ifdef TRACERS_WATER
C**** the following diagnostics are set assuming that the particular
C**** tracer exists in water.
C**** Tracers in precipitation (=Wet deposition)
      k = k+1
      tij_prec = k
      if (dowetdep(n)) then
        if (to_per_mil(n) .eq.1) then
          write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_prec'
          write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *         ' in Precip'
          units_tij(k,n)=cmrwt(n)
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_wet_dep'
          write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *         ' Wet Deposition'
          units_tij(k,n)=unit_string(ijtc_power(n)-5,trim(cmrwt(n))
     *         //'/s')
          scale_tij(k,n)=10.**(-ijtc_power(n)+5)/dtsrc
        end if
      end if
C**** Tracers in evaporation
      if (src_dist_index(n)<=1) then
        k = k+1
        tij_evap = k
        if (tr_wd_type(n).eq.nWater) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_in_evap'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' in Evaporation'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(ijtc_power(n),cmrwt(n))
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          units_tij(k,n)=unit_string(ijtc_power(n),trim(cmrwt(n))//'/s')
          scale_tij(k,n)=10.**(-ijtc_power(n))/dtsrc
        end if
#ifdef TRACERS_SPECIAL_O18
C**** Tracers at sea surface
        k = k+1
        tij_owiso = k
        write(sname_tij(k,n),'(a,i2)')
     &        trim(TRNAME(n))//'_Sea_Surface'
        write(lname_tij(k,n),'(a,i2)')
     &        trim(TRNAME(n))//' at Sea Surface'
        !Convert to permil if specified:
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=cmrwt(n)
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        else
          units_tij(k,n)='kg/kg fresh water'
          !Scale quantity to account for extra surface flux time steps:
          scale_tij(k,n)=1.d0/REAL(NIsurf,KIND=8)
          !Set special denominator name in order to use ocean fraction:
          dname_tij(k,n) = 'ocnfrac'
        endif
#endif
        endif !if (tr_wd_type(n).eq.nWater)
C**** Tracers in river runoff (two versions - for inflow and outflow)
        k = k+1
        tij_rvr = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_rvr'
        lname_tij(k,n) = trim(TRNAME(n))//' in River Inflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
        k = k+1
        tij_rvro = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_rvro'
        lname_tij(k,n) = trim(TRNAME(n))//' in River Outflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers in iceberg runoff 
        k = k+1
        tij_icb = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_icb'
        lname_tij(k,n) = trim(TRNAME(n))//' in Iceberg Inflow'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers in sea ice
        k = k+1
        atmice%tij_seaice = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_ice'
        lname_tij(k,n) = trim(TRNAME(n))//' in Sea Ice'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,cmrwt(n))
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
          dname_tij(k,n) = 'oicefr'
c        denom_tij(k,n)=n_Water ! if kg/kg units for non-water-isotopes
        end if
C**** Tracers conc. in ground component (ie. water or ice surfaces)
        k = k+1
        tij_grnd = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_at_Grnd'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' at Ground'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
          denom_tij(k,n)=n_Water
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)/REAL(NIsurf,KIND=8)
        end if
C**** Tracers conc. in lakes (layer 1)
        k = k+1
        tij_lk1 = k
        sname_tij(k,n) = trim(TRNAME(n))//'_Lake1'
        lname_tij(k,n) = trim(TRNAME(n))//' Lakes layer 1'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers conc. in lakes (layer 2)
        k = k+1
        tij_lk2 = k
        sname_tij(k,n) = trim(TRNAME(n))//'_Lake2'
        lname_tij(k,n) = trim(TRNAME(n))//' Lakes layer 2'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers conc. in soil water
        k = k+1
        tij_soil = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_Soil'
        lname_tij(k,n) = trim(TRNAME(n))//' Soil Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracers conc. in land snow water
        k = k+1
        tij_snow = k
        sname_tij(k,n) = trim(TRNAME(n))//'_in_Snow'
        lname_tij(k,n) = trim(TRNAME(n))//' Land Snow Water'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)+3,'kg/kg wat')
          scale_tij(k,n)=10.**(-ijtc_power(n)-3)
        end if
        denom_tij(k,n)=n_Water
C**** Tracer ice-ocean flux
        k = k+1
        atmice%tij_icocflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_ic_oc_flx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Ice-Ocean Flux'
        if (to_per_mil(n) .eq.1) then
          units_tij(k,n)=unit_string(0,cmrwt(n))
          denom_tij(k,n)=n_Water
          scale_tij(k,n)=1.
        else
          units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
          scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
        end if
C**** Tracers integrated E-W atmospheric flux
        k = k+1
        tij_uflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_uflx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' E-W Atmos Flux'
        units_tij(k,n)=unit_string(ijtc_power(n)+10,'kg/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)-10)/DTsrc
C**** Tracers integrated N-S atmospheric flux
        k = k+1
        tij_vflx = k
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_vflx'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' N-S Atmos Flux'
        units_tij(k,n)=unit_string(ijtc_power(n)+10,'kg/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)-10)/DTsrc
C**** Tracers integrated E-W sea ice flux
        k = k+1
        atmice%tij_tusi = k
        sname_tij(k,n) = trim(TRNAME(n))//'_tusi'
        lname_tij(k,n) = trim(TRNAME(n))//' E-W Ice Flux'
        units_tij(k,n) = unit_string(ntrocn(n),'kg/s')
        scale_tij(k,n) = (10.**(-ntrocn(n)))/DTsrc
C**** Tracers integrated N-S sea ice flux
        k = k+1
        atmice%tij_tvsi = k
        sname_tij(k,n) = trim(TRNAME(n))//'_tvsi'
        lname_tij(k,n) = trim(TRNAME(n))//' N-S Ice Flux'
        units_tij(k,n) = unit_string(ntrocn(n),'kg/s')
        scale_tij(k,n) = (10.**(-ntrocn(n)))/DTsrc
      endif ! if (src_dist_index(n)<=1) 
#endif /* TRACERS_WATER */
#ifdef TRACERS_DRYDEP
C**** Tracers dry deposition flux.
      k = k+1
      tij_drydep = k
      if (dodrydep(n)) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_dry_dep'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Dry Deposition'
        units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
      end if
      k = k+1
      tij_gsdep = k
      if (trradius(n).gt.0) then
        write(sname_tij(k,n),'(a,i2)') trim(TRNAME(n))//'_gs_dep'
        write(lname_tij(k,n),'(a,i2)') trim(TRNAME(n))//
     *       ' Gravitational Settling'
        units_tij(k,n)=unit_string(ijtc_power(n)-5,'kg/m^2/s')
        scale_tij(k,n)=10.**(-ijtc_power(n)+5)/DTsrc
      end if
#endif

      if (k .gt. ktaij) then
        if (AM_I_ROOT()) write (6,*)
     &   'tij_defs: Increase ktaij=',ktaij,' to at least ',k
        call stop_model('ktaij too small',255)
      end if

      end do !ntm

c
c Collect denominator short names for later use
c
      do n=1,ntm
      do k=1,ktaij
        nd = denom_tij(k,n)
        if(nd.gt.0) dname_tij(k,n) = sname_tij(k,nd)
      enddo
      enddo

      RETURN
      END SUBROUTINE set_generic_tracer_diags

      SUBROUTINE sum_prescribed_tracer_2Dsources(dtstep)
!@sum apply_tracer_2Dsource adds surface sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE GEOM, only : imaxj,byaxyp
      USE QUSDEF, only : mz,mzz
      USE TRACER_COM, only : NTM,ntsurfsrc
#ifdef TRACERS_TOMAS
     &     ,n_ASO4,n_ANACL,n_AECOB,n_AOCOB,n_ADUST,n_SO2
#endif
#ifndef SKIP_TRACER_SRCS
      USE FLUXES, only : trsource
#endif
      USE FLUXES, only : trflux1,atmsrf
      USE TRDIAG_COM, only : taijs=>taijs_loc
      USE TRDIAG_COM, only : ijts_source,jls_source,itcon_surf
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,ns,naij,najl,j,i
      REAL*8, DIMENSION(grid%I_STRT_HALO:grid%I_STOP_HALO
     *     ,grid%J_STRT_HALO:grid%J_STOP_HALO) :: dtracer

      INTEGER :: J_0, J_1, I_0, I_1
      INTEGER ntsurf !same as ntsurfsrc

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      
C**** This is tracer independent coding designed to work for all
C**** surface sources.
C**** Note that tracer flux is added to first layer either implicitly
C**** in ATURB or explicitly in 'apply_fluxes_to_atm' call in SURFACE.

      do n=1,ntm
        trflux1(:,:,n) = 0.
        ntsurf = ntsurfsrc(n) 
#ifdef TRACERS_TOMAS

! Overwrite with first bin 
       if(n.ge.n_ASO4(1).and.n.lt.n_ANACL(1))
     .               ntsurf=ntsurfsrc(n_SO2) !so4
       if(n.ge.n_AECOB(1).and.n.lt.n_AOCOB(1))
     .               ntsurf=ntsurfsrc(n_AECOB(1)) !ecob
       if(n.ge.n_AOCOB(1).and.n.lt.n_ADUST(1))
     .               ntsurf=ntsurfsrc(n_AOCOB(1)) !ocob + ocil
#endif
#ifndef SKIP_TRACER_SRCS
        do ns=1,ntsurf
C**** diagnostics
          naij = ijts_source(ns,n)
          IF (naij > 0) THEN
          taijs(:,:,naij) = taijs(:,:,naij) + trsource(:,:,ns,n)*dtstep
          ENDIF
          najl = jls_source(ns,n)
          IF (najl > 0) THEN
            DO J=J_0,J_1
              DO I=I_0,imaxj(j)
                call inc_tajls(i,j,1,najl,trsource(i,j,ns,n)*dtstep)
              END DO
            END  DO
          END IF
          DO J=J_0,J_1
            do i=i_0,imaxj(j)
              dtracer(i,j)=trsource(i,j,ns,n)*dtstep
            end do
          end do
          if (itcon_surf(ns,n).gt.0)
     *         call DIAGTCB(dtracer,itcon_surf(ns,n),n)
C**** trflux1 is total flux into first layer
          trflux1(:,:,n) = trflux1(:,:,n)+trsource(:,:,ns,n)*byaxyp(:,:)
        end do
#endif
        atmsrf%trflux_prescr(n,:,:) = trflux1(:,:,n)
      end do

#ifdef TRACERS_TOMAS
#ifdef ALT_EMISS_COAG
! The subgridcoag_drv_2d call (which adjusts trflux_prescr) has been
! moved from SURFACE in order to avoid "double-counting" trflux_prescr.
! trflux_prescr currently affects the interactive surface fluxes, but
! subgridcoag_drv_2d does not account for this.   Application of the
! full coagulation increment to trflux_prescr BEFORE the interactive
! surface fluxes is most consistent with the current model structure.
! Other possible routes not taken:
! (1) Inclusion of coagulation tendency terms within the interactive
!     surface flux calculation (complicated).
! (2) Modification of subgridcoag_drv_2d to only see the part of
!     trflux_prescr not consumed by downward interactive fluxes
!     (may not reflect the original intent).
C**** Apply subgrid coagulation for freshly emitted particles.
      call subgridcoag_drv_2D(dtstep)
#endif
#endif

      RETURN
      END SUBROUTINE sum_prescribed_tracer_2Dsources

      SUBROUTINE set_strattroptracer_diag(dtstep)
!@sum safe Tracer Fluxes at the Tropopause
!@auth Susanne Bauer
      USE GEOM, only : imaxj,byaxyp
      USE TRACER_COM, only : NTM
      USE TRDIAG_COM, only : taijn=>taijn_loc
      USE TRDIAG_COM, only : TSCF3D=>tscf3d_loc
      USE TRDIAG_COM, only : tij_strop
      USE ATM_COM,    only : LTROPO
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds

      IMPLICIT NONE
      INTEGER n,j,i
      REAL*8, INTENT(IN) :: dtstep
      INTEGER :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

#ifndef SKIP_TRACER_DIAGS
      do j=j_0,j_1
      do i=i_0,imaxj(j)
      do n=1,ntm
      taijn(i,j,tij_strop,n) = taijn(i,j,tij_strop,n)
     &         + TSCF3D(i,j,ltropo(i,j),n)*byaxyp(i,j)/dtstep

      end do ! tracer n
      enddo
      enddo
#endif /*SKIP_TRACER_DIAGS*/


      RETURN
      END SUBROUTINE set_strattroptracer_diag

      SUBROUTINE apply_tracer_2Dsource(dtstep)
!@sum apply_tracer_2Dsource adds surface sources to tracers
!@auth Jean Lerner/Gavin Schmidt
      USE GEOM, only : imaxj,axyp
      USE QUSDEF, only : mz,mzz
      USE TRACER_COM, only : NTM,trm,trmom
      USE FLUXES, only : trflux1,atmsrf
      USE DOMAIN_DECOMP_ATM, ONLY : GRID, getDomainBounds
      use oldtracer_mod, only: src_dist_index
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: dtstep
      INTEGER n,i,j
      REAL*8 ftr1,dewflux,tinyReal8

      INTEGER :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      
      tinyReal8=tiny(dewflux)

C**** This is tracer independent coding designed to work for all
C**** surface sources.
C**** Note that tracer flux is added to first layer either implicitly
C**** in ATURB or explicitly in 'apply_fluxes_to_atm' call in SURFACE.

      do n=1,ntm

       do j=j_0,j_1
C**** modify vertical moments (only from non-interactive sources)
c this is disabled until vertical moments are also modified
c during the vertical transport between layer 1 and the others
c        trmom( mz,:,j,1,n) = trmom( mz,:,j,1,n)-1.5*trflux_prescr(:,j,n)
c     *       *dtstep*axyp(:,j)
c        trmom(mzz,:,j,1,n) = trmom(mzz,:,j,1,n)+0.5*trflux_prescr(:,j,n)
c     *       *dtstep*axyp(:,j)

C**** Add prescribed and interactive sources
c        trflux1(:,j,n) = trflux1(:,j,n)+atmsrf%trsrfflx(n,:,j)
        trflux1(:,j,n) =
     &        atmsrf%trflux_prescr(n,:,j)+atmsrf%trsrfflx(n,:,j)
       end do

C**** Technically speaking the vertical moments should be modified here
C**** as well. But for consistency with water vapour we only modify
C**** moments for dew.
       if (src_dist_index(n)==0) then
        do j=J_0,J_1
          do i=i_0,imaxj(j)
            dewflux=-atmsrf%trsrfflx(n,i,j)*axyp(i,j)*dtstep
            ! The previous criteria were: 
            ! if(atmsrf%trsrfflx(n,i,j).lt.0 .and. trm(i,j,1,n).gt.0)then
            ! The new "dewflux > tinyReal8" takes care of the first of those,
            ! plus a saftey margin so we don't divide by an exceedingly small
            ! number. The max() in the denominator of ftr1 already prevents
            ! a negative ftr1, but I leave in the trm criterion just so we
            ! don't alter moments when trm is negative (as before).
            if ( dewflux > tinyReal8 .and. trm(i,j,1,n) > 0.) then
              ftr1=dewflux/max(dewflux,trm(i,j,1,n))
              trmom(:,i,j,1,n)=trmom(:,i,j,1,n)*(1.-ftr1)
            end if
          end do
        end do
       endif
      end do
C****
      RETURN
      END SUBROUTINE apply_tracer_2Dsource

      SUBROUTINE TDECAY
!@sum TDECAY decays radioactive tracers every source time step
!@auth Gavin Schmidt/Jean Lerner
      USE RESOLUTION, only: im,jm,lm
      USE MODEL_COM, only : itime,dtsrc
#ifndef SKIP_TRACER_SRCS
      USE FLUXES, only : tr3Dsource
#endif
      USE GEOM, only : imaxj
      use OldTracer_mod, only: itime_tr0, trname, trdecay
      use TRACER_COM, only: nChemistry
      USE TRACER_COM, only : NTM
     &     ,trm,trmom,n_Pb210, n_Rn222
#ifdef TRACERS_WATER
     *     ,trwm
      USE SEAICE_COM, only : si_atm,si_ocn
#ifndef TRACERS_ATM_ONLY
      USE LAKES_COM, only : trlake
      USE LANDICE_COM, only : trlndi,trsnowli
      USE GHY_COM, only : tr_w_ij,tr_wsn_ij
#endif
#endif
      USE TRDIAG_COM, only : jls_decay,itcon_decay
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      IMPLICIT NONE
      real*8, dimension(ntm) :: expdec
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &                  grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: told
      logical, save :: ifirst=.true.
      integer n,najl,j,l,i
      integer :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      if (ifirst) then
        expdec = 1.
        ifirst = .false.
      end if

      do n=1,ntm
        if (trdecay(n).gt.0. .and. itime.ge.itime_tr0(n)) then
          expdec(n)=exp(-trdecay(n)*dtsrc)
C**** Atmospheric decay
          told(:,:,:)=trm(:,:,:,n)

#ifdef TRACERS_WATER
     *               +trwm(:,:,:,n)
          trwm(:,:,:,n)   = expdec(n)*trwm(:,:,:,n)
#endif
#ifndef SKIP_TRACER_SRCS
          if (trname(n) .eq. "Rn222" .and. n_Pb210.gt.0) then
            tr3Dsource(:,:,:,nChemistry,n_Pb210)=
     *        trm(:,:,:,n)*(1-expdec(n))*210./222./dtsrc
          end if
#endif

          trm(:,:,:,n)    = expdec(n)*trm(:,:,:,n)
          trmom(:,:,:,:,n)= expdec(n)*trmom(:,:,:,:,n)

#ifdef TRACERS_WATER
C**** Note that ocean tracers are dealt with by separate ocean code.
C**** Decay sea ice tracers
#ifndef TRACERS_ATM_ONLY
          si_atm%trsi(n,:,:,:)   = expdec(n)*si_atm%trsi(n,:,:,:)
#endif
          if(si_atm%grid%im_world .ne. si_ocn%grid%im_world) then
            call stop_model(
     &           'TDECAY: tracers in sea ice are no longer on the '//
     &           'atm. grid - please move the next line',255)
          endif
#ifndef TRACERS_ATM_ONLY
          si_ocn%trsi(n,:,:,:)   = expdec(n)*si_ocn%trsi(n,:,:,:)
C**** ...lake tracers
          trlake(n,:,:,:) = expdec(n)*trlake(n,:,:,:)
C**** ...land surface tracers
          tr_w_ij(n,:,:,:,:) = expdec(n)*tr_w_ij(n,:,:,:,:)
          tr_wsn_ij(n,:,:,:,:)= expdec(n)*tr_wsn_ij(n,:,:,:,:)
          trsnowli(n,:,:,:) = expdec(n)*trsnowli(n,:,:,:)
          trlndi(n,:,:,:)   = expdec(n)*trlndi(n,:,:,:)
#endif
#endif
C**** atmospheric diagnostics
          najl = jls_decay(n)
          do l=1,lm
          do j=J_0,J_1
            do i=I_0,imaxj(j)
              call inc_tajls(i,j,l,najl,trm(i,j,l,n)
#ifdef TRACERS_WATER
     *             +trwm(i,j,l,n)
#endif
     *             -told(i,j,l))
            enddo
          enddo
          enddo

          call DIAGTCA(itcon_decay(n),n)
        end if
      end do
C****
      return
      end subroutine tdecay


      SUBROUTINE TRGRAV
!@sum TRGRAV gravitationally settles particular tracers
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : grav,deltx,lhe,rgas,visc_air
      USE RESOLUTION, only: im,jm,lm
      USE MODEL_COM, only : itime,dtsrc
      USE ATM_COM, only : t,q
      USE GEOM, only : imaxj,byaxyp
      USE SOMTQ_COM, only : mz,mzz,mzx,myz,zmoms
      USE ATM_COM, only : gz,pmid,pk
      use OldTracer_mod, only: trradius, itime_tr0, trname, trpdens
      USE TRACER_COM, only : NTM,trm,trmom
#ifdef TRACERS_AMP
      USE TRACER_COM, only : ntmAMPi,ntmAMPe
      USE AmpTracersMetadata_mod, only: AMP_MODES_MAP, AMP_NUMB_MAP
      USE AMP_AEROSOL, only : DIAM, AMP_dens
      USE AERO_SETUP,  only : CONV_DPAM_TO_DGN
#endif
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only : nbins,n_ASO4,xk
      USE CONSTANT,   only : pi 
#endif
      USE TRDIAG_COM, only : jls_grav
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      IMPLICIT NONE
      real*8 :: stokevdt,press,fgrfluxd,qsat,vgs,tr_radius,tr_dens,temp
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,lm) :: told,airden,visc,rh
     *     ,gbygz
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO) :: fluxd, fluxu
      integer n,najl,i,j,l
#ifndef TRACERS_TOMAS
     &     ,nAMP
#endif
      integer :: J_0, J_1, I_0, I_1
      logical :: hydrate
#ifdef TRACERS_TOMAS
      integer binnum,k
!@var vs : gravitational settling velocity at each bin (m s-1)
      real*8, dimension(grid%I_STRT_HALO:grid%I_STOP_HALO,
     &     grid%J_STRT_HALO:grid%J_STOP_HALO,lm,NBINS) :: vs !gravitational settling velocity (m s-1)
!@var Dp_gr : particle diameter (m)
      real*8 Dp_gr(nbins)         
      real*8 density_gr(nbins)  !density (kg/m3) of current size bin
      real*8 mp                 !particle mass (kg)
      real*8 mu                 !air viscosity (kg/m s)
#endif

      call getDomainBounds(grid, J_STRT = J_0, J_STOP = J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Calculate some tracer independent arrays      
C**** air density + relative humidity (wrt water) + air viscosity
      do l=1,lm
        do j=J_0,J_1
          do i=I_0,imaxj(j)
            press=pmid(l,i,j)
            temp=pk(l,i,j)*t(i,j,l)
            airden(i,j,l)=100.d0*press/(rgas*temp*(1.+q(i,j,l)*deltx))
            rh(i,j,l)=q(i,j,l)/qsat(temp,lhe,press)
            visc(i,j,l)=visc_air(temp)
            if (l.eq.1) then
              gbygz(i,j,l)=0.
            else
              gbygz(i,j,l)=grav/(gz(i,j,l)-gz(i,j,l-1))
            end if
          end do
        end do
      end do

C**** Gravitational settling
      do n=1,ntm
        if (trradius(n).gt.0. .and. itime.ge.itime_tr0(n)) then
C**** need to hydrate the sea salt before determining settling
          hydrate = (trname(n).eq.'seasalt1'.or.trname(n).eq.'seasalt2')

          fluxd=0.
          do l=lm,1,-1          ! loop down
            do j=J_0,J_1
              do i=I_0,imaxj(j)

C*** save original tracer mass
                told(i,j,l)=trm(i,j,l,n)
C**** set incoming flux from previous level
                fluxu(i,j)=fluxd(i,j)

C**** set particle properties
                tr_dens = trpdens(n)
                tr_radius = trradius(n)

#ifdef TRACERS_AMP
       if (n.ge.ntmAMPi.and.n.le.ntmAMPe) then
         nAMP=n-ntmAMPi+1
        if(AMP_MODES_MAP(nAMP).gt.0) then
        if(DIAM(i,j,l,AMP_MODES_MAP(nAMP)).gt.0.) then
        if(AMP_NUMB_MAP(nAMP).eq. 0) then    ! Mass
        tr_radius=0.5*DIAM(i,j,l,AMP_MODES_MAP(nAMP))
        else                              ! Number
        tr_radius=0.5*DIAM(i,j,l,AMP_MODES_MAP(nAMP))
     +            *CONV_DPAM_TO_DGN(AMP_MODES_MAP(nAMP))
        endif

        call AMPtrdens(i,j,l,n)
        tr_dens =AMP_dens(i,j,l,AMP_MODES_MAP(nAMP))
        endif   
        endif   
       endif 
#endif  

#ifndef TRACERS_TOMAS
C**** calculate stokes velocity (including possible hydration effects
C**** and slip correction factor)
                stokevdt=dtsrc*vgs(airden(i,j,l),rh(i,j,l),tr_radius
     *               ,tr_dens,visc(i,j,l),hydrate)
#else 
       
       if(n.lt.n_ASO4(1))then
!     no size resolved aerosol tracer (e.g. NH4)
          stokevdt=dtsrc*vgs(airden(i,j,l),rh(i,j,l),tr_radius
     *         ,tr_dens,visc(i,j,l),hydrate)
          
       elseif(n.ge.n_ASO4(1)) then
         if(n.eq.n_ASO4(1))THEN
C     02/20/2012 - TOMAS trgrav is modified to be able to reproduce the model output
           call dep_getdp(i,j,l,Dp_gr,density_gr) 
           do k=1,nbins
C     APR 2015 - FIX vs with slip correction factor (use vgs now)
             
cyhl              vs(I,J,L,k)=density_gr(k)*(Dp_gr(k)**2)*grav
cyhl     *             /18.d0/visc(i,j,l) 

             vs(I,J,L,k) = 
     *            vgs(airden(i,j,l),0.d0,Dp_gr(k)/2.,density_gr(k),
     *            visc(i,j,l),hydrate) 
           enddo
         endif
          binnum=mod(N-n_ASO4(1)+1,NBINS)
          if (binnum.eq.0) binnum=NBINS
          stokevdt=dtsrc*vs(i,j,l,binnum) !grav. settling velocity for TOMAS model
       endif !size-resolved aerosols

#endif
C**** Calculate height differences using geopotential
C**** Next line causes problems in high vertical resolution models. Limit it for now:
C****           fgrfluxd=stokevdt*gbygz(i,j,l) 
                fgrfluxd=min(stokevdt*gbygz(i,j,l),1.d0) 
                fluxd(i,j) = trm(i,j,l,n)*fgrfluxd ! total flux down
                trm(i,j,l,n) = trm(i,j,l,n)*(1.-fgrfluxd)+fluxu(i,j)
                if (1.-fgrfluxd.le.1d-16) trm(i,j,l,n) = fluxu(i,j)
                trmom(zmoms,i,j,l,n)=trmom(zmoms,i,j,l,n)*(1.-fgrfluxd)
              end do
            end do
          end do
          najl = jls_grav(n)
          IF (najl > 0) THEN
            do l=1,lm
              do j=J_0,J_1
                do i=I_0,imaxj(j)
                  call inc_tajls(i,j,l,najl,trm(i,j,l,n)-told(i,j,l))
                enddo
              enddo
            enddo
          END IF
        end if
      end do

C****

      return
      end subroutine trgrav

      REAL*8 FUNCTION vgs(airden,rh1,tr_radius,tr_dens,visc,hydrate)
!@sum vgs returns settling velocity for tracers (m/s)
!@auth Gavin Schmidt/Reha Cakmur
      USE CONSTANT, only : by3,pi,gasc,avog,rt2,deltx
     *     ,mair,grav
      IMPLICIT NONE
      real*8, intent(in) ::  airden,rh1,tr_radius,tr_dens,visc
      logical, intent(in) :: hydrate
      real*8  wmf,frpath
      real*8, parameter :: dair=3.65d-10 !m diameter of air molecule
C**** coefficients for slip factor calculation
      real*8, parameter :: s1=1.247d0, s2=0.4d0, s3=1.1d0
C**** coefficients for hydration radius if required 
      real*8, parameter :: c1=0.7674d0, c2=3.079d0, c3=2.573d-11,
     *     c4=-1.424d0
      real*8 dens, rad, rh

C**** Determine if hydration of aerosols is required
      if (hydrate) then
        rh=max(0.01d0,min(rh1,0.99d0))
C**** hydrated radius (Gerber, 1988)
        rad=(c1*tr_radius**(c2)/(c3*tr_radius**(c4)-log10(rh))
     *       + tr_radius**3)**by3
C**** hydrated density
        dens=((rad**3 - tr_radius**3)*1000.d0
     *       + tr_radius**3*tr_dens)/rad**3
      else                      ! take dry values
        dens = tr_dens
        rad = tr_radius 
      end if

C**** calculate stokes velocity
      vgs=2.*grav*dens*rad**2/(9.*visc)

C**** slip correction factor
c wmf is the additional velocity if the particle size is small compared
c   to the mean free path of the air; important in the stratosphere
      frpath=1d-3*mair/(pi*rt2*avog*airden*(dair)**2.)
      wmf=frpath/tr_radius*(s1+s2*exp(-s3*tr_radius/frpath))
      vgs=(1.d0+wmf)*vgs
C****
      return
      end function vgs

      SUBROUTINE read_monthly_sources(iu,jdlast,tlca,tlcb,data,
     *  frac,imon)
!@sum Read in monthly sources and interpolate to current day
!@+   Calling routine must have the lines:
!@+      real*8 tlca(im,jm,nm),tlcb(im,jm,nm)
!@+      integer imon(nm)   ! nm=number of files that will be read
!@+      data jdlast /0/
!@+      save jdlast,tlca,tlcb,imon
!@+   Input: iu, the fileUnit#; jdlast
!@+   Output: interpolated data array + two monthly data arrays
!@auth Jean Lerner and others
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_ATM, only : getDomainBounds, AM_I_ROOT
      USE DOMAIN_DECOMP_ATM, only : READT_PARALLEL, REWIND_PARALLEL
      USE RESOLUTION, only: im,jm
      use model_com, only: modelEclock
      USE JulianCalendar_mod, only: idofm=>JDmidOfM
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR, 
     &  INT_MONTHS_PER_YEAR
      implicit none
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     tlca,tlcb,data
      real*8 :: frac
      integer :: imon
      integer :: iu,jdlast

      integer :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

      if (jdlast.EQ.0) then ! NEED TO READ IN FIRST MONTH OF DATA
        imon=1          ! imon=January
        if (modelEclock%getDayOfYear().le.16)  then ! JDAY in Jan 1-15, first month is Dec
          CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),tlca,12)
          CALL REWIND_PARALLEL( iu )
        else            ! JDAY is in Jan 17 to Dec 16, get first month
  120     imon=imon+1
          if (modelEclock%getDayOfYear().gt.idofm(imon) 
     &        .AND. imon.le.INT_MONTHS_PER_YEAR)
     &        go to 120
          CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),tlca,imon-1)
          if (imon.eq.13)  CALL REWIND_PARALLEL( iu )
        end if
      else              ! Do we need to read in second month?
        if (modelEclock%getDayOfYear().ne.jdlast+1) then ! Check that data is read in daily
          if (modelEclock%getDayOfYear().ne.1 .OR. 
     *      jdlast.ne.INT_DAYS_PER_YEAR) then
            if (AM_I_ROOT()) write(6,*)
     *      'Incorrect values in Tracer Source:JDAY,JDLAST=',
     *      modelEclock%getDayOfYear(),JDLAST
            call stop_model('stopped in TRACERS.f',255)
          end if
          imon=imon-INT_MONTHS_PER_YEAR  ! New year
          go to 130
        end if
        if (modelEclock%getDayOfYear().le.idofm(imon)) go to 130
        imon=imon+1     ! read in new month of data
        tlca(I_0:I_1,J_0:J_1) = tlcb(I_0:I_1,J_0:J_1)
        if (imon.eq.13) CALL REWIND_PARALLEL( iu  )
      end if
      CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),tlcb,1)
  130 continue
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-modelEclock%getDayOfYear())
     & / (idofm(imon)-idofm(imon-1))
      data(I_0:I_1,J_0:J_1) = tlca(I_0:I_1,J_0:J_1)*frac + 
     & tlcb(I_0:I_1,J_0:J_1)*(1.-frac)
      return
      end subroutine read_monthly_sources
#endif  /* TRACERS_ON */

      subroutine checktr(subr)
!@sum  CHECKTR Checks whether atmos tracer variables are reasonable
!@vers 2013/03/26
!@auth Gavin Schmidt
#ifdef TRACERS_ON
      USE CONSTANT, only : teeny
      USE RESOLUTION, only: im,jm,lm
      USE ATM_COM, only : q,qcl,qci
      USE GEOM, only : axyp,imaxj
      USE SOMTQ_COM, only : qmom
      USE ATM_COM, only : MA
      USE FLUXES, only : atmocn,atmice,atmgla,atmlnd
      use OldTracer_mod, only: trname, t_qlimit
      USE TRACER_COM, only: ntm, trmom, trm, nmom
#ifdef TRACERS_WATER
      USE TRACER_COM, only: trwm
#endif
      USE DOMAIN_DECOMP_ATM, ONLY: GRID, getDomainBounds, AM_I_ROOT
      IMPLICIT NONE
      LOGICAL QCHECKT
      INTEGER I,J,L,N,m, imax,jmax,lmax
      REAL*8 relerr, errmax,errsc,tmax,amax,qmax,wmax,twmax,qmomax(nmom)
     *     ,tmomax(nmom)
#ifdef TRACERS_WATER
      real*8 :: qc
#endif

!@var SUBR identifies where CHECK was called from
      CHARACTER*6, INTENT(IN) :: SUBR
      INTEGER :: J_0, J_1, nj, I_0,I_1
      call getDomainBounds(GRID, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      nj = J_1 - J_0 + 1

      !CALL CHECK4(gtracer(1,1,1,J_0),NTM,4,IM,nJ,SUBR,'GTRACE')
      CALL CHECK3(atmocn%gtracer(1,1,J_0),NTM,IM,nJ,SUBR,'GTRACO')
      CALL CHECK3(atmice%gtracer(1,1,J_0),NTM,IM,nJ,SUBR,'GTRACI')
      CALL CHECK3(atmgla%gtracer(1,1,J_0),NTM,IM,nJ,SUBR,'GTRACG')
      CALL CHECK3(atmlnd%gtracer(1,1,J_0),NTM,IM,nJ,SUBR,'GTRACE')
      do n=1,NTM
        CALL CHECK4(trmom(:,:,J_0:J_1,:,n),NMOM,IM,nJ,LM,SUBR,
     *       'X'//trname(n))
        CALL CHECK3(trm(:,J_0:J_1,:,n),IM,nJ,LM,SUBR,trname(n))
#ifdef TRACERS_WATER
        CALL CHECK3(trwm(:,J_0:J_1,:,n),IM,nJ,LM,SUBR,'QCL'//trname(n))
#endif

C**** check for negative tracer amounts (if t_qlimit is set)
        if (t_qlimit(n)) then
          QCHECKT=.false.
          do l=1,lm
          do j=j_0,j_1
          do i=i_0,imaxj(j)
            if (trm(i,j,l,n).lt.0) then
              if (AM_I_ROOT())
     *         write(6,*) "Negative mass for ",trname(n),i,j,l,trm(i,j,l
     *             ,n)," after ",SUBR,"."
              QCHECKT=.true.
            end if
          end do
          end do
          end do
          if (QCHECKT)
     &         call stop_model("CHECKTR: Negative tracer amount",255)
        end if

C**** check whether air mass is conserved

        if (trname(n).eq.'Air') then
          errmax = 0. ; lmax=1 ; imax=I_0 ; jmax=J_0; tmax=0. ; amax=0.
          do l=1,lm
          do j=j_0,j_1
          do i=i_0,imaxj(j)
            relerr=abs(trm(i,j,l,n)-ma(l,i,j)*axyp(i,j))/
     *           (ma(l,i,j)*axyp(i,j))
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
              tmax=trm(i,j,l,n) ; amax=ma(l,i,j)*axyp(i,j)
            end if
          end do
          end do
          end do
          print*,"Relative error in air mass after ",trim(subr),":",imax
     *         ,jmax,lmax,errmax,tmax,amax
        end if

#ifdef TRACERS_WATER
        if (trname(n).eq.'Water') then
          errmax = 0. ; lmax=1 ; imax=I_0 ; jmax=J_0
          tmax=0. ; twmax=0. ; qmax=0. ; wmax=0.
          tmomax = 0. ; qmomax = 0. 
          do l=1,lm
          do j=j_0,j_1
          do i=i_0,imaxj(j)
            errsc=(q(i,j,l)+sum(abs(qmom(:,i,j,l))))*ma(l,i,j)*axyp(i,j)
            if (errsc.eq.0.) errsc=1.
            relerr=abs(trm(i,j,l,n)-q(i,j,l)*ma(l,i,j)*axyp(i,j))/errsc
            qc = qcl(i,j,l)+qci(i,j,l) ! add liquid and ice to compare to trwm
            if (qc.gt.0 .and. trwm(i,j,l,n).gt.1.) relerr
     *           =max(relerr,(trwm(i,j,l,n)-qc*ma(l,i,j)*
     *           axyp(i,j))/(qc*ma(l,i,j)*axyp(i,j)))
            if ((qc.eq.0 .and.trwm(i,j,l,n).gt.1) .or. 
     *           (qc.gt.teeny .and.trwm(i,j,l,n).eq.0))
     *           print*,"Condensate mismatch: ",subr,i,j,l,
     &           trwm(i,j,l,n),qc*ma(l,i,j)*axyp(i,j)
            do m=1,nmom
              relerr=max(relerr,(trmom(m,i,j,l,n)-qmom(m,i,j,l)*ma(l,i,j
     *             )*axyp(i,j))/errsc)
            end do
            if (relerr.gt.errmax) then
              lmax=l ; imax=i ; jmax=j ; errmax=relerr
              tmax=trm(i,j,l,n) ; qmax=q(i,j,l)*ma(l,i,j)*axyp(i,j)
              twmax=trwm(i,j,l,n) ; wmax=qc*ma(l,i,j)*axyp(i,j)
              tmomax(:)=trmom(:,i,j,l,n)
              qmomax(:)=qmom(:,i,j,l)*ma(l,i,j)*axyp(i,j)
            end if
          end do
          end do
          end do
          print*,"Relative error in water mass after ",trim(subr),":"
     *         ,imax,jmax,lmax,errmax,tmax,qmax,twmax,wmax,(tmomax(m)
     *         ,qmomax(m),m=1,nmom) 
        end if
#endif
      end do
#endif
      return
      end subroutine checktr

      subroutine setup_emis_sectors
!@sum setup_emis_sectors reads from the rundeck the 
!@+ geographic regions and sectors associated with tracer
!@+ emissions and saves names.
!@auth Greg Faluvegi

      use TRACER_COM, only : n_max_sect,
     & n_max_reg,alter_sources,ef_REG_IJ,
     & ef_fact,num_sectors,sect_name
      USE DOMAIN_DECOMP_ATM, only: GRID,getDomainBounds
      use DOMAIN_DECOMP_ATM, only: AM_I_ROOT,writet_parallel
      USE GEOM, only: lat2d_dg, lon2d_dg, imaxj
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      use Dictionary_mod, only : sync_param
      use EmissionRegion_mod, only: initializeEmissionsRegions
      use EmissionRegion_mod, only: regions, numRegions

      implicit none

      integer :: i,j,n,iu
      character*80 :: title
      character*2 :: fnum
      character*124 :: sectors_are

      sectors_are=' '
      call sync_param("sectors_are",sectors_are)

      call initializeEmissionsRegions()

! see how many sectors there are, save names in array:
      num_sectors=0
      i=1
      do while(i < len(sectors_are))
        j=index(sectors_are(i:len(sectors_are))," ")
        if (j > 1) then
          num_sectors=num_sectors+1
          i=i+j
        else
          i=i+1
        end if
      enddo
      if (num_sectors > n_max_sect) call stop_model
     &("n_max_sect must be increased",255)
      if(num_sectors > 0 ) read(sectors_are,*)
     & sect_name(1:num_sectors)

      end subroutine setup_emis_sectors

      subroutine setup_emis_sectors_regions
!@sum Reads the factors associated with each sector and
!@+ region. Output IJ map of regions.
!@auth Greg Faluvegi

      use TRACER_COM, only : n_max_sect,
     & n_max_reg,alter_sources,ef_REG_IJ,
     & ef_fact,num_sectors,sect_name
      USE DOMAIN_DECOMP_ATM, only: GRID,getDomainBounds
      use DOMAIN_DECOMP_ATM, only: AM_I_ROOT,writet_parallel
      USE GEOM, only: lat2d_dg, lon2d_dg, imaxj
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      use Dictionary_mod, only : sync_param
      use EmissionRegion_mod, only: initializeEmissionsRegions
      use EmissionRegion_mod, only: regions, numRegions

      implicit none

      integer :: i,j,n,iu
      character*80 :: title
      character*2 :: fnum
      character*124 :: sectors_are

      INTEGER J_0, J_1, I_0, I_1
      call getDomainBounds(grid, J_STRT=J_0,J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      sectors_are=' '
      call sync_param("sectors_are",sectors_are)

      call initializeEmissionsRegions()

! read the actual emission altering factors:
      do n=1,num_sectors
         write(fnum,'(I2.2)') n
        ef_fact(n,:)=1.d0
        call sync_param("SECT_"//fnum,ef_fact(n,:),numRegions)
      enddo

! check if there is any altering requested:
      alter_sources = .false.
      do i=1,num_sectors; do j=1,numRegions
        if(ef_fact(i,j)/=1.0 .and. ef_fact(i,j)/=-1.d30)
     &  alter_sources = .true.
      enddo ; enddo

! write out regions to GISS-format IJ file, each restart:
      if(alter_sources)then
        ef_REG_IJ(I_0:I_1,J_0:J_1)=0.d0
        do j=J_0,J_1 
          do i=i_0,imaxj(j)
            do n=1,numRegions
              if (regions(n)%hasLatLon(lat2d_dg(i,j), lon2d_dg(i,j)))
     &             ef_REG_IJ(i,j)= max(ef_REG_IJ(i,j),dble(n))
            enddo
          enddo
        enddo
        title='Regions defined in rundeck for altering tracer sources'
        call openunit('EF_REG',iu,.true.)
        call WRITET_PARALLEL(grid,iu,nameunit(iu),ef_REG_IJ,title)
        call closeunit(iu)
      endif

      end subroutine setup_emis_sectors_regions

      subroutine tracerIO(fid, action)
!@sum tracerIO() provides a generic interface for IO actions
!@+   on the full list of tracers.
!@auth T. Clune
      use ParallelIo_mod
      use pario, only : read_data,defvar,write_data
      use domain_decomp_atm, only : grid
      USE Dictionary_mod
      USE TRACER_COM, only: ntm, TRmom, TRM, coupled_chem
      USE TRACER_COM, only: ntm, nmom, no3_live, oh_live
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM, only: yNO3,pHOx,pNOx,pOx,yCH3O2,yC2O3,
     &yROR,yXO2,yAldehyde,yXO2N,yRXPAR,ss,ydms,yso2,sulfate,pNO3
#ifdef TRACERS_dCO
     &,ydC217O3,ydC218O3,yd13C2O3
     &,yd13CXPAR
     &,yd17OROR,yd18OROR,yd13CROR
     &,yd17Oald,yd18Oald,yd13Cald
     &,ydCH317O2,ydCH318O2,yd13CH3O2
#endif  /* TRACERS_dCO */
     &,SF3,SF2,pClOx,pClx,pOClOx,pBrOx,yCl2,yCl2O2
     &,mostRecentNonZeroAlbedo
#ifdef INTERACTIVE_WETLANDS_CH4 
      use TRACER_SOURCES, only: day_ncep,DRA_ch4,sum_ncep,PRS_ch4,
     & HRA_ch4,iday_ncep,i0_ncep,iHch4,iDch4,i0ch4,first_ncep,first_mod,
     & avg_model,avg_ncep
#endif
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef BC_ALB
      USE AEROSOL_SOURCES, only : snosiz
#endif  /* BC_ALB */
      use trdiag_com, only: trcSurfMixR_acc,trcSurfByVol_acc
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      USE fluxes,ONLY : pprec,pevap
      USE trdust_mod,ONLY : hbaij,ricntd
      use trdust_drv, only: def_rsf_trdust
      use trdust_drv, only: new_io_trdust
#endif
#ifdef TRACERS_WATER
      USE TRACER_COM, only: trwm
#endif
#if (defined CUBED_SPHERE) || (defined TRACERS_VOLCEXP)
      USE TRACER_COM, only: daily_z
#endif
      use OldTracer_mod, only: trName
      use model_com, only : ioread,iowrite

#ifdef TRACERS_AMP
      use amp_aerosol, only : diam,nactv
#endif

      implicit none

      integer, intent(in) :: fid
      character(len=*), intent(in) :: action

      type (ParallelIo) :: handle
      character(len=:), allocatable :: ijldims
#ifdef TRACERS_SPECIAL_Shindell
      character(len=:), allocatable :: ijcdims
#endif
      integer :: n

      ijldims='(dist_im,dist_jm,lm)' 
#ifdef TRACERS_SPECIAL_Shindell
      ijcdims='(dist_im,dist_jm,topLevelOfChemistry)'
#endif
      handle = ParallelIo(grid, fid)

      do n=1,NTM
        call doVar(handle,action,trm(:,:,:,n),
     &        'trm_'//trim(trname(n))//ijldims)
        call doVar(handle,action,trmom(:,:,:,:,n),
     &       'trmom_'//trim(trname(n))//'(nmom,dist_im,dist_jm,lm)',
     &       jdim=3)
#ifdef TRACERS_WATER
        call doVar(handle,action,trwm(:,:,:,n),
     &       'trwm_'//trim(trname(n))//ijldims)
#endif
      enddo

#if (defined CUBED_SPHERE) || (defined TRACERS_VOLCEXP)
c daily_z is currently only needed for CS
      call doVar(handle,action,daily_z,'daily_z'//ijldims)
#endif

#ifdef TRACERS_SPECIAL_Shindell       

      handle = ParallelIo(grid, fid, 'TRACERS_SPECIAL_Shindell')

      call doVar(handle,action,ss,
     & 'ss(n_rj,topLevelOfChemistry,dist_im,dist_jm)',jdim=4)
      call doVar(handle,action,yNO3,'yNO3'//ijcdims)
      call doVar(handle,action,pHOx,'pHOx'//ijcdims)
      call doVar(handle,action,pNOx,'pNOx'//ijcdims)
      call doVar(handle,action,pNO3,'pNO3'//ijcdims)
      call doVar(handle,action,pOx ,'pOx'//ijcdims)
      call doVar(handle,action,yCH3O2,'yCH3O2'//ijcdims)
#ifdef TRACERS_dCO
      call doVar(handle,action,ydCH317O2,'ydCH317O2'//ijcdims)
      call doVar(handle,action,ydCH318O2,'ydCH318O2'//ijcdims)
      call doVar(handle,action,yd13CH3O2,'yd13CH3O2'//ijcdims)
#endif  /* TRACERS_dCO */
      call doVar(handle,action,yC2O3,'yC2O3'//ijcdims)
#ifdef TRACERS_dCO
      call doVar(handle,action,ydC217O3,'ydC217O3'//ijcdims)
      call doVar(handle,action,ydC218O3,'ydC218O3'//ijcdims)
      call doVar(handle,action,yd13C2O3,'yd13C2O3'//ijcdims)
#endif  /* TRACERS_dCO */
      call doVar(handle,action,yROR,'yROR'//ijcdims)
#ifdef TRACERS_dCO
      call doVar(handle,action,yd17OROR,'yd17OROR'//ijcdims)
      call doVar(handle,action,yd18OROR,'yd18OROR'//ijcdims)
      call doVar(handle,action,yd13CROR,'yd13CROR'//ijcdims)
#endif  /* TRACERS_dCO */
      call doVar(handle,action,yXO2,'yXO2'//ijcdims)
      call doVar(handle,action,yXO2N,'yXO2N'//ijcdims)
      call doVar(handle,action,yAldehyde,'yAldehyde'//ijcdims)
#ifdef TRACERS_dCO
      call doVar(handle,action,yd17Oald,'yd17Oald'//ijcdims)
      call doVar(handle,action,yd18Oald,'yd18Oald'//ijcdims)
      call doVar(handle,action,yd13Cald,'yd13Cald'//ijcdims)
#endif  /* TRACERS_dCO */
      call doVar(handle,action,yRXPAR,'yRXPAR'//ijcdims)
#ifdef TRACERS_dCO
      call doVar(handle,action,yd13CXPAR,'yd13CXPAR'//ijcdims)
#endif  /* TRACERS_dCO */
      call doVar(handle,action,ydms,'ydms'//ijcdims)
      call doVar(handle,action,ySO2,'ySO2'//ijcdims)
      call doVar(handle,action,sulfate,'sulfate'//ijldims) ! stays ijldims
      if(trim(action) == 'read_dist') then
           ! read_dist is a badly chosen synonym for read
        if(is_set_param("coupled_chem"))
     &       call get_param( "coupled_chem", coupled_chem )
      endif
      if(coupled_chem == 1) then
        call doVar(handle,action,oh_live,'oh_live'//ijldims)   ! stays ijldims
        call doVar(handle,action,no3_live,'no3_live'//ijldims) ! stays ijldims
      endif
      call doVar(handle,action,SF3,'SF3'//ijcdims)
      call doVar(handle,action,SF2,'SF2'//ijcdims)
      call doVar(handle,action,pClOx,'pClOx'//ijcdims)
      call doVar(handle,action,pClx,'pClx'//ijcdims)
      call doVar(handle,action,pOClOx,'pOClOx'//ijcdims)
      call doVar(handle,action,pBrOx,'pBrOx'//ijcdims)
      call doVar(handle,action,yCl2,'yCl2'//ijcdims)
      call doVar(handle,action,yCl2O2,'yCl2O2'//ijcdims)

#ifdef INTERACTIVE_WETLANDS_CH4 
      handle = ParallelIo(grid, fid, 'INTERACTIVE_WETLANDS_CH4')

      call doVar(handle,action,day_ncep,
     &     'day_ncep(dist_im,dist_jm,max_days,nra_ncep)')
      call doVar(handle,action,dra_ch4,
     &     'dra_ch4(dist_im,dist_jm,max_days,nra_ch4)')
      call doVar(handle,action,sum_ncep,
     &     'sum_ncep(dist_im,dist_jm,nra_ncep)')
      call doVar(handle,action,prs_ch4,
     &     'prs_ch4(dist_im,dist_jm,nra_ch4)')
      call doVar(handle,action,HRA_ch4,
     &     'HRA_ch4(dist_im,dist_jm,maxHR_ch4,nra_ch4)')
      call doVar(handle,action,i0ch4,
     &     'i0ch4(dist_im,dist_jm,nra_ch4)')
      call doVar(handle,action,iDch4,
     &     'iDch4(dist_im,dist_jm,nra_ch4)')
      call doVar(handle,action,iHch4,
     &     'iHch4(dist_im,dist_jm,nra_ch4)')
      call doVar(handle,action,first_mod,
     &     'first_mod(dist_im,dist_jm,nra_ch4)')
      call doVar(handle,action,avg_model,
     &     'avg_model(dist_im,dist_jm,nra_ch4)')
      call doVar(handle,action,avg_ncep, 
     &     'avg_ncep(dist_im,dist_jm,nra_ncep)')
      select case (action)
        case('define')
          call defvar(grid,fid,iday_ncep,'iday_ncep(nra_ncep)')
          call defvar(grid,fid,i0_ncep,'i0_ncep(nra_ncep)')
          call defvar(grid,fid,first_ncep,'first_ncep(nra_ncep)')
        case('read_dist')
          call read_data(grid,fid,'iday_ncep',iday_ncep,
     &       bcast_all=.true.)
          call read_data(grid,fid,'i0_ncep',i0_ncep,
     &       bcast_all=.true.)
          call read_data(grid,fid,'first_ncep',first_ncep,
     &       bcast_all=.true.)
        case ('write_dist')
          call write_data(grid,fid,'iday_ncep',iday_ncep)
          call write_data(grid,fid,'i0_ncep',i0_ncep)
          call write_data(grid,fid,'first_ncep',first_ncep)
      end select
!        call doVar(handle,action,iday_ncep,'iday_ncep(nra_ncep)')
!        call doVar(handle,action,i0_ncep,'i0_ncep(nra_ncep)')
!        call doVar(handle,action,first_ncep,'first_ncep(nra_ncep)')
#endif /* INTERACTIVE_WETLANDS_CH4 */
#endif /* TRACERS_SPECIAL_Shindell */

      call doVar(handle,action,trcSurfMixR_acc
     &     ,'trcSurfMixR_acc(dist_im,dist_jm,Ntm)')
      call doVar(handle,action,trcSurfByVol_acc
     &     ,'trcSurfByVol_acc(dist_im,dist_jm,Ntm)')

#ifdef TRACERS_SPECIAL_Shindell
      handle = ParallelIo(grid, fid,'TRACERS_SPECIAL_Shindell')
      call doVar(handle,action,mostRecentNonZeroAlbedo,
     & 'mostRecentNonZeroAlbedo(dist_im,dist_jm)')
#endif

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      handle = ParallelIo(grid, fid,
     &   'TRACERS_DUST||TRACERS_MINERALS')
      call doVar(handle,action,hbaij,'hbaij(dist_im,dist_jm)')
      call doVar(handle,action,ricntd,'ricntd(dist_im,dist_jm)')
      call doVar(handle,action,pprec,'pprec(dist_im,dist_jm)')
      call doVar(handle,action,pevap,'pevap(dist_im,dist_jm)')

      select case (action)
      case ('define')
         call def_rsf_trdust(fid)
      case ('read_dist')
         call new_io_trdust(fid,ioread)
      case ('write_dist')
         call new_io_trdust(fid,iowrite)
      end select

#endif

#ifdef BC_ALB
      call doVar(handle,action,snosiz,'snosiz(dist_im,dist_jm)')
#endif  /* BC_ALB */

#ifdef TRACERS_AMP
      ! restartability hack until matrix code refactored to
      ! re-diagnose these qtys on demand
      call doVar(handle,action,diam,
     &     'amp_diam(dist_im,dist_jm,lm,nmodes)')
      call doVar(handle,action,nactv,
     &     'amp_nactv(dist_im,dist_jm,lm,nmodes)')
#endif

      return
      end subroutine tracerIo
     

#ifdef CACHED_SUBDD
      subroutine tijh_defs(arr,nmax,decl_count)
! 2D tracer outputs (model horizontal grid).
! Each tracer output must be declared separately (no bundling).
      use model_com, only : dtsrc,nday
      use subdd_mod, only : info_type, sched_rad, reduc_max
      use OldTracer_mod, only: trname
#ifdef TRACERS_WATER 
      use OldTracer_mod, only : nWater, tr_wd_type
#endif
      use tracer_com, only : ntm
      use trdiag_com, only : to_volume_MixRat, save_dry_aod
      use radpar, only: nraero_aod=>NTRACE
      use rad_com, only: ntrix_aod,nraero_rf,ntrix_rf,diag_fc
      use RunTimeControls_mod, only: tracers_amp, tracers_tomas
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
! types of aods to be saved
! The name will be any combination of {,TRNAME}{as,cs,dry}{,a}aod
      character(len=10), dimension(3) :: ssky=(/'as ','cs ','dry'/),
     &                                lsky=(/'All-sky  ','Clear-sky',
     &                                       'Dry aeros'/)
      character(len=10), dimension(2) :: sabs=(/' ','a'/),
     &                                labs=(/'          ','absorption'/)
      character(len=10), dimension(2) :: sfrc=(/'swf','lwf'/),
     &                                lfrc=(/'shortwave','longwave '/)
      character(len=10) :: spcname
! types of PM/tracer surface amounts to be saved
! The name will be any combination of PM{2p5,10}{l1,s}{m,c}
!                            I.E. {species}{location}{units}
! and similar format for any tracer: trname(){l1,s}{m,c}.
! In practice did not include the l1c (L=1 cocentration) case
      character(len=20), dimension(2) :: 
     &   ssiz=(/'2p5','10 '/), lsiz=(/'PM2.5','PM10 '/),
     &   sloc=(/'l1','s '/),   lloc=(/'L=1    ','Surface'/),
     &   sunt=(/'m','c'/),
     &   lunt=(/'Mass Mixing Ratio','Concentration    '/),
     &   uunt=(/'kg species / kg air','kg m-3             '/)
      character*80 :: unitString,unitString2
      integer :: s,a,n,f,u,l,p

      decl_count = 0

! Optical Depths

      do s=1,size(ssky)
      if (ssky(s).eq.'dry' .and. save_dry_aod.eq.0) cycle
      do a=1,size(sabs)
      do n=1,nraero_aod+1 ! +1 for total
        if (n<=nraero_aod) then
          spcname = trim(trname(ntrix_aod(n)))
        else
          spcname = ''
        endif
        arr(next()) = info_type_(
     &    sname = trim(spcname)//trim(ssky(s))//trim(sabs(a))//'aod',
     &    lname = trim(spcname)//' '//trim(lsky(s))//' '//
     &            trim(labs(a))//' aerosol optical depth',
     &    units = '-',
     &    sched = sched_rad
     &       )
      enddo ! n
      enddo ! a
      enddo ! s

! Forcing

      do f=1,size(sfrc)
      do n=1,nraero_rf
        if (diag_fc==2) then
          spcname = trim(trname(ntrix_rf(n)))
        else if (diag_fc==1) then
          if (tracers_amp) then
            spcname='AMP'
          elseif (tracers_tomas) then
            spcname='TOMAS'
          else
            spcname='OMA'
          endif
        endif
        arr(next()) = info_type_(
     &    sname = trim(sfrc(f))//'_'//trim(spcname),
     &    lname = trim(spcname)//' '//trim(lfrc(f))//' forcing',
     &    units = 'W m-2',
     &    sched = sched_rad
     &       )
      enddo ! n
      enddo ! f

! Surface Tracer Amount

      do n=1,ntm
        ! L=1 and surface mixing ratios:
        u=1
        if (to_volume_MixRat(n) == 1) then
          unitString='mole species / mole air'
          unitString2='Volume Mixing Ratio'
        else
          unitString=trim(uunt(u))
          unitString2=trim(lunt(u))
        endif 
        do l=1,size(sloc) 
          arr(next()) = info_type_(
     &    sname = trim(trname(n))//trim(sloc(l))//trim(sunt(u)),
     &    lname=
     &    trim(trname(n))//' '//trim(lloc(l))//' '//trim(unitString2),
     &    units = trim(unitString)
     &    )
        end do
        ! surface concentrations:
        u=2
        l=2
        arr(next()) = info_type_(
     &  sname = trim(trname(n))//trim(sloc(l))//trim(sunt(u)),
     &  lname=trim(trname(n))//' '//trim(lloc(l))//' '//trim(lunt(u)),
     &  units = trim(uunt(u))
     &  )
      end do ! ntm

! Tracer Load (column mass)

      do n=1,ntm
        arr(next()) = info_type_(
     &  sname = trim(trname(n))//'load',
     &  lname = trim(trname(n))//' Column Mass',
     &  units = 'kg m-2'
     &  )
      end do ! ntm

#ifdef TRACERS_AMP
      arr(next()) = info_type_(
     &  sname = 'ampBCload',
     &  lname = 'BC Column Mass',
     &  units = 'kg m-2'
     &  )
      arr(next()) = info_type_(
     &  sname = 'ampDustload',
     &  lname = 'Dust Column Mass',
     &  units = 'kg m-2'
     &  )
      arr(next()) = info_type_(
     &  sname = 'ampNH4load',
     &  lname = 'NH4 Column Mass',
     &  units = 'kg m-2'
     &  )
      arr(next()) = info_type_(
     &  sname = 'ampNO3load',
     &  lname = 'NO3 Column Mass',
     &  units = 'kg m-2'
     &  )
      arr(next()) = info_type_(
     &  sname = 'ampOAload',
     &  lname = 'OA Column Mass',
     &  units = 'kg m-2'
     &  )
      arr(next()) = info_type_(
     &  sname = 'ampSO4load',
     &  lname = 'SO4 Column Mass',
     &  units = 'kg m-2'
     &  )
      arr(next()) = info_type_(
     &  sname = 'ampSSload',
     &  lname = 'SS Column Mass',
     &  units = 'kg m-2'
     &  )
#endif

! Surface Particulate Matter Amount

      do p=1,size(ssiz)
        ! L=1 and surface mixing ratios (always mass):
        u=1
        do l=1,size(sloc)
          arr(next()) = info_type_(
     &    sname = 'PM'//trim(ssiz(p))//trim(sloc(l))//trim(sunt(u)),
     &    lname = 
     &     trim(lsiz(p))//' '//trim(lloc(l))//' '//trim(lunt(u)),
     &    units = trim(uunt(u))
     &    )
        end do 
        ! surface concentrations:
        u=2
        l=2
        arr(next()) = info_type_(
     &  sname = 'PM'//trim(ssiz(p))//trim(sloc(l))//trim(sunt(u)),
     &  lname =
     &   trim(lsiz(p))//' '//trim(lloc(l))//' '//trim(lunt(u)),
     &  units = trim(uunt(u))
     &  )
      end do ! p (PM size)

#ifdef TRACERS_SPECIAL_Shindell
      arr(next()) = info_type_(
     &  sname = 'MRNO2l1', ! because not a tracer
     &  lname = 'L=1 NO2 mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRNOl1', ! because not a tracer
     &  lname = 'L=1 NO mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRO3l1max', ! because not a tracer
     &  lname = 'Maximum Daily L=1 O3 mixing ratio',
     &  units = 'mole species / mole air',
     &  reduc = reduc_max
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'O3col', ! not "load", to contrast with tracers
     &  lname = 'O3 Column Mass',
     &  units = 'kg m-2'
     &  )
#endif

#ifdef TRACERS_WATER
! Water tracer/isotope precipitation 
      do n=1,ntm
        if (tr_wd_type(n).eq.nWater) then !Is it a water (isotope) tracer?
          !Set 'subdd' variable/object meta-data:
          arr(next()) = info_type_(
     &      sname = trim(trname(n))//'_in_prec',
     &      lname = trim(trname(n))//' in Precip',
     &      units = 'kg/m^2/s',
     &      scale = 1./dtsrc !kg/m2 -> kg/m2/s
     &      )
        end if
      end do !ntm

! Water tracer/isotope evaporation
      do n=1,ntm
        if (tr_wd_type(n).eq.nWater) then !Is it a water (isotope) tracer?
          !Set 'subdd' variable/object meta-data:
          arr(next()) = info_type_(
     &      sname = trim(trname(n))//'_in_evap',
     &      lname = trim(trname(n))//' in Evap',
     &      units = 'kg/m^2/s',
     &      scale = 1./dtsrc !kg/m2 -> kg/m2/s
     &      )
        end if
      end do !ntm
#endif

      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine tijh_defs

      subroutine tijlh_defs(arr,nmax,decl_count)
! 3D tracer outputs (model horizontal grid and layers).
! Each tracer output must be declared separately (no bundling).
      use model_com, only : dtsrc,nday
      use subdd_mod, only : info_type, sched_rad
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      use tracer_com, only : ntm
      use OldTracer_mod, only: trname
      use radpar, only: nraero_aod=>NTRACE
      use rad_com, only: ntrix_aod
      use trdiag_com, only : to_volume_MixRat, save_dry_aod
      implicit none
      integer :: nmax,decl_count
      integer :: n
      character*80 :: unitString
      type(info_type) :: arr(nmax)
! types of aods to be saved
! The name will be any combination of {,TRNAME}{as,cs,dry}{,a}aod3d
      character(len=10), dimension(3) :: ssky=(/'as ','cs ','dry'/),
     &                                lsky=(/'All-sky  ','Clear-sky',
     &                                       'Dry aeros'/)
      character(len=10), dimension(2) :: sabs=(/' ','a'/),
     &                               labs=(/'          ','absorption'/),
     &                               lcoef=(/'extinction','absorption'/)
      character(len=10) :: spcname
      integer :: s,a

      decl_count = 0

      ! First, diagnostics available for all tracers:
      do n=1,ntm

        ! 3D mixing ratios (SUBDD string is just tracer name):
        if (to_volume_MixRat(n) == 1) then
          unitString='mole species / mole air'
        else
          unitString='kg species / kg air'
        endif
        arr(next()) = info_type_(
     &    sname = trim(trname(n)),
     &    lname = trim(trname(n))//' mixing ratio',
     &    units = trim(unitString)
     &    )

#ifdef TRACERS_AMP
! AMP aerosol diameters
        spcname=trim(trname(n))
        if ((spcname(1:2) == 'N_') .and. (spcname(6:7) == '_1')) then
          arr(next()) = info_type_(
     &      sname = 'd'//trim(trname(n)),
     &      lname = trim(trname(n))//' mass mean diameter',
     &      units = 'm'
     &         )
        endif
#endif  /* TRACERS_AMP */

      end do ! tracers loop

! 3d AOD

      do s=1,size(ssky)
      if (ssky(s).eq.'dry' .and. save_dry_aod.eq.0) cycle
      do a=1,size(sabs)
      do n=1,nraero_aod+1 ! +1 for total
        if (n<=nraero_aod) then
          spcname = trim(trname(ntrix_aod(n)))
        else
          spcname = ''
        endif
        arr(next()) = info_type_(
     &    sname = trim(spcname)//trim(ssky(s))//trim(sabs(a))//'aod3d',
     &    lname = trim(spcname)//' '//trim(lsky(s))//' '//
     &            trim(labs(a))//' aerosol optical depth',
     &    units = '-',
     &    sched = sched_rad
     &       )
        arr(next()) = info_type_(
     &    sname = trim(spcname)//trim(ssky(s))//trim(sabs(a))//
     &    'bcoef3d',
     &    lname = trim(spcname)//' '//trim(lsky(s))//' '//
     &            trim(lcoef(a))//' coefficient',
     &    units = 'm-1',
     &    sched = sched_rad
     &       )
      enddo ! n
      enddo ! a
      enddo ! s

      ! Other tracer diags on model levels:

#ifdef TRACERS_SPECIAL_Shindell
      arr(next()) = info_type_(
     &  sname = 'MRNO2', ! because not a tracer
     &  lname = 'NO2 mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRNO', ! because not a tracer
     &  lname = 'NO mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRO3', ! because not a tracer
     &  lname = 'O3 mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'OH_conc', ! because not a tracer
     &  lname = 'OH concentration',
     &  units = 'molecules cm-3'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'HO2_conc', ! because not a tracer
     &  lname = 'HO2 concentration',
     &  units = 'molecules cm-3'
     &  )

      arr(next()) = info_type_(
     &  sname = 'JO1D', ! because not a tracer
     &  lname = 'O3-->O1D+O2 photolysis rate',
     &  units = 's-1'
     &  )

      arr(next()) = info_type_(
     &  sname = 'JNO2', ! because not a tracer
     &  lname = 'NO2-->NO+O photolysis rate',
     &  units = 's-1'
     &  )
#endif /* TRACERS_SPECIAL_Shindell */

      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine tijlh_defs


      subroutine tijph_defs(arr,nmax,decl_count)
! 3D tracer outputs (model horizontal grid and constant pressure levels)
! Each tracer output must be declared separately (no bundling).
      use model_com, only : dtsrc,nday
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      use tracer_com, only : ntm
      use OldTracer_mod, only: trname
      use trdiag_com, only : to_volume_MixRat
      implicit none
      integer :: nmax,decl_count
      integer :: n
      character*80 :: unitString
      type(info_type) :: arr(nmax)

      decl_count = 0

      ! First, diagnostics available for all tracers:
      do n=1,ntm

        ! 3D mixing ratios (SUBDD string is just tracer name with cp
        ! appended):
        if (to_volume_MixRat(n) .eq.1) then
          unitString='mole species / mole air'
        else
          unitString='kg species / kg air'
        endif
        arr(next()) = info_type_(
     &    sname = trim(trname(n))//'cp',
     &    lname = trim(trname(n))//' mixing ratio',
     &    units = trim(unitString)
     &    )

      end do ! tracers loop

      ! Other tracer diags on constant pressure levels:

#ifdef TRACERS_SPECIAL_Shindell
      arr(next()) = info_type_(
     &  sname = 'MRNO2cp', ! because not a tracer
     &  lname = 'NO2 mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRNOcp', ! because not a tracer
     &  lname = 'NO mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'MRO3cp', ! because not a tracer
     &  lname = 'O3 mixing ratio',
     &  units = 'mole species / mole air'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'OH_conccp', ! because not a tracer
     &  lname = 'OH concentration',
     &  units = 'molecules cm-3'
     &  )
C
      arr(next()) = info_type_(
     &  sname = 'HO2_conccp', ! because not a tracer
     &  lname = 'HO2 concentration',
     &  units = 'molecules cm-3'
     &  )
#endif /* TRACERS_SPECIAL_Shindell */

      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine tijph_defs


      subroutine accumCachedTracerSUBDDs

      use domain_decomp_atm, only : grid
      USE resolution, only: LM
      use geom, only : byaxyp
      use atm_com, only    : byma
      use tracer_com, only : ntm,trm,mass2vol
      use OldTracer_mod, only: trname, pm10fact, pm2p5fact
      use trdiag_com, only : to_volume_MixRat,trcsurf,trcSurfByVol
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups
     &     ,inc_subdd,find_groups, LmaxSUBDD
#ifdef TRACERS_AMP
      use AMP_AEROSOL, only: ampPM2p5, ampPM10
#endif
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd
      integer :: L, n, k
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: sddarr2d
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  LM                               ) :: sddarr3d
      real*8 :: convert

      ! Tracer 3D diags on model levels
      call find_groups('taijlh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          ntm_loop: do n=1,ntm
            ! tracer 3D mixing ratios (SUBDD names are just tracer name):
            if(trim(trname(n)).eq.trim(subdd%name(k))) then
              if (to_volume_MixRat(n) == 1) then
                convert=mass2vol(n)
              else
                convert=1.d0
              endif
              do L=1,LmaxSUBDD
                sddarr3d(:,:,L) = 
     &          trm(:,:,L,n)*convert*byaxyp(:,:)*byma(L,:,:)
              end do
              call inc_subdd(subdd,k,sddarr3d)
              exit ntm_loop
            end if
          end do ntm_loop
        enddo ! k
      enddo ! igroup

      ! Tracer 3D diags on constant pressure levels
      call find_groups('taijph',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          ntm_loop2: do n=1,ntm
            ! tracer 3D mixing ratios (SUBDD names are tracer name with cp appended):
            if(trim(trname(n))//'cp'.eq.trim(subdd%name(k))) then
              if (to_volume_MixRat(n) == 1) then
                convert=mass2vol(n)
              else
                convert=1.d0
              endif
              do L=1,LM ! not LmaxSUBDD in case pressure requested above that
                sddarr3d(:,:,L) =
     &          trm(:,:,L,n)*convert*byaxyp(:,:)*byma(L,:,:)
              end do
              call inc_subdd(subdd,k,sddarr3d)
              exit ntm_loop2
            end if
          end do ntm_loop2
        enddo ! k
      enddo ! igroup

      ! Tracer 2D I-J diags
      call find_groups('taijh',grpids,ngroups)
      do igrp=1,ngroups
      subdd => subdd_groups(grpids(igrp))
      diag_loop: do k=1,subdd%ndiags
        ntm_loop3: do n=1,ntm

          ! tracer surface mixing ratios:
          if(trim(trname(n))//'sm'.eq.trim(subdd%name(k))) then
            if (to_volume_MixRat(n) == 1) then
              sddarr2d(:,:)=trcsurf(:,:,n)*mass2vol(n)
            else
              sddarr2d(:,:)=trcsurf(:,:,n)
            endif
            call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop
          end if

          ! tracer surface concentrations:
          if(trim(trname(n))//'sc'.eq.trim(subdd%name(k))) then
            sddarr2d(:,:)=trcSurfByVol(:,:,n)
            call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop
          end if

          ! tracer L=1 mixing ratios:
          if(trim(trname(n))//'l1m'.eq.trim(subdd%name(k))) then
            if (to_volume_MixRat(n) == 1) then
              sddarr2d(:,:)=
     &          trm(:,:,1,n)*mass2vol(n)*byaxyp(:,:)*byma(1,:,:)
            else
              sddarr2d(:,:)=trm(:,:,1,n)*byaxyp(:,:)*byma(1,:,:)
            endif
            call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop
          end if

          ! tracer column load:
          if(trim(trname(n))//'load'.eq.trim(subdd%name(k))) then
            sddarr2d(:,:)=sum(trm(:,:,:,n),dim=3)*byaxyp(:,:)
            call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop
          end if

        enddo ntm_loop3

! Particulate matter to be treated differently for mass-based 
! aerosols or not:
#ifdef TRACERS_TOMAS
        select case(trim(subdd%name(k)))
        case('PM2p5sm','PM2p5l1m','PM2p5sc',
     &    'PM10sm','PM10l1m','PM10sc')
          call tomas_pm_subdd_accum(subdd,k,trim(subdd%name(k)))
          cycle diag_loop
        end select
#elif (defined TRACERS_AMP)
        select case(trim(subdd%name(k)))
        ! L=1 PM2.5 mass mixing ratio:
         case('PM2p5l1m')
         sddarr2d(:,:)= ampPM2p5(:,:)     ! kg/kg air
         call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop

        ! L=1 PM10 mass mixing ratio:
         case('PM10l1m')
         sddarr2d(:,:)= ampPM10(:,:)      ! kg/kg air
         call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop

        end select
#else
        select case(trim(subdd%name(k)))
        ! surface PM2.5 mass mixing ratio:
        case('PM2p5sm')
          sddarr2d(:,:)=0.d0
          do n=1,ntm
            if(pm2p5fact(n)/=0.)
     &      sddarr2d(:,:)=sddarr2d(:,:)+pm2p5fact(n)*trcsurf(:,:,n)
          end do
          call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop 

        ! L=1 PM2.5 mass mixing ratio:
        case('PM2p5l1m')
          sddarr2d(:,:)=0.d0
          do n=1,ntm
            if(pm2p5fact(n)/=0.)
     &      sddarr2d(:,:)=sddarr2d(:,:)+pm2p5fact(n)*
     &            trm(:,:,1,n)*byaxyp(:,:)*byma(1,:,:)
          end do
          call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop

        ! surface PM2.5 concentration:
        case('PM2p5sc')
          sddarr2d(:,:)=0.d0
          do n=1,ntm
            if(pm2p5fact(n)/=0.)
     &      sddarr2d(:,:)=sddarr2d(:,:)+pm2p5fact(n)*trcSurfByVol(:,:,n)
          end do
          call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop

        ! surface PM10 mass mixing ratio:
        case('PM10sm')
          sddarr2d(:,:)=0.d0
          do n=1,ntm
            if(pm10fact(n)/=0.)
     &      sddarr2d(:,:)=sddarr2d(:,:)+pm10fact(n)*trcsurf(:,:,n)
          end do
          call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop  

        ! L=1 PM10 mass mixing ratio:
        case('PM10l1m')
          sddarr2d(:,:)=0.d0
          do n=1,ntm
            if(pm10fact(n)/=0.)
     &      sddarr2d(:,:)=sddarr2d(:,:)+pm10fact(n)*
     &            trm(:,:,1,n)*byaxyp(:,:)*byma(1,:,:)
          end do
          call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop

        ! surface PM10 concentration:
        case('PM10sc')
          sddarr2d(:,:)=0.d0
          do n=1,ntm
            if(pm10fact(n)/=0.)
     &      sddarr2d(:,:)=sddarr2d(:,:)+pm10fact(n)*trcSurfByVol(:,:,n)
          end do
          call inc_subdd(subdd,k,sddarr2d) ; cycle diag_loop

        end select
#endif /* -not- TRACERS_TOMAS, TRACERS_AMP sections */

      enddo diag_loop
      enddo ! igroup

      end subroutine accumCachedTracerSUBDDs
#endif /* CACHED_SUBDD */

#if (defined TRACERS_SPECIAL_Shindell) || (defined TRACERS_AEROSOLS_Koch) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS) ||\
    (defined TRACERS_GASEXCH_GCC)

      subroutine get_aircraft_tracer
     & (nTracer,fileName,year,xday,phi,need_read,AIRCstream)
!@sum  get_aircraft_tracer to define the 3D source of tracers from aircraft
!@auth Drew Shindell? / Greg Faluvegi / Jean Learner
      use RESOLUTION, only : im,jm,lm
      use model_com, only: itime, master_yr
      use domain_decomp_atm, only: GRID,getDomainBounds,write_parallel
      use constant, only: bygrav
      use filemanager, only: openunit,closeunit,is_fbsa
      use fluxes, only: tr3Dsource
      use geom, only: axyp
      use OldTracer_mod, only: itime_tr0, trname
      use OldTracer_mod, only: set_first_aircraft, first_aircraft
      use TRACER_COM, only: ntm_chem_beg,ntm_chem_end,nAircraft
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) || \
    (defined TRACERS_TOMAS)
      use TRACER_COM, only: aer_int_yr
      use TRACER_COM, only: SO2_int_yr
      use TRACER_COM, only: NH3_int_yr
      use TRACER_COM, only: BC_int_yr
      use TRACER_COM, only: OC_int_yr
#endif
      use Dictionary_mod, only: is_set_param, get_param
      use RAD_COM, only: o3_yr
      use timestream_mod, only : read_stream, timestream, init_stream

      IMPLICIT NONE
 
!@param Laircr the number of layers of aircraft data read from file
      INTEGER, PARAMETER :: Laircr=25
!@param aircraft_Tyr1, aircraft_Tyr2 the starting and ending years
!@+     for transient tracer aircraft emissions (= means non transient)
      integer :: aircraft_Tyr1=0,aircraft_Tyr2=0
!@var airtracer 3D source of tracer from aircraft (on model levels)
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM)
     &     :: airtracer
!@var fileName the name of the aircraft source file for this tracer
      character(len=*), intent(IN) :: fileName
!@var nTracer the index of the tracer in current call in ntm arrays
!@+   for example n_NOx or n_M_BC1_BC
      integer, intent(IN) :: year,xday,nTracer
      integer :: xyear
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM),
     &     intent(IN) :: phi
      logical, intent(IN) :: need_read

      integer :: fileUnit 
      integer L,i,j,k,LL
      interface
        real*8 function get_src_fact(n,OA_not_OC)
          integer, intent(in) :: n
          logical, intent(in), optional :: OA_not_OC
        end function get_src_fact
      end interface

!@var src holds the tracer source returned from actual reading routine
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Laircr):: src
!@var zmod approx. geometric height at model layer(m), phi/grav
      real*8, dimension(LM) :: zmod
!@var zairL heights of CMIP5,CMIP6 aircraft emissions (km)
      real*4, parameter, dimension(Laircr) :: zairL = ! alt in km:
     & (/0.305, 0.915, 1.525, 2.135, 2.745, 3.355, 3.965, 4.575, 5.185,
     & 5.795, 6.405, 7.015, 7.625, 8.235001, 8.845, 9.455001, 10.065,
     & 10.675, 11.285, 11.895, 12.505, 13.115, 13.725, 14.335, 14.945/)
      integer :: J_1, J_0, I_0, I_1, do_ppm
      logical :: trans_emis=.false.,isItFbsa=.true.
      integer :: yr1=0, yr2=0, copy_master_yr, cyclic_yr
      ! note the AIRCstream is passed before init_stream is called for it
      ! (and in the case of fbsa files init_stream is never called for it.)
      type(timestream) :: AIRCstream
 
! Aircraft tracer source input is monthly, on 25 levels.
! Read it in here and interpolated each day.

      ! for fortran binary sequential access files, transient emissions/
      ! start/end years are determined from rundeck parameters:
      isItFbsa=is_fbsa(fileName)
      if (isItFbsa) then
        if (is_set_param("aircraft_Tyr1")) then
          call get_param("aircraft_Tyr1",aircraft_Tyr1)
        else
          call stop_model("Must provide aircraft_Tyr1 via rundeck",255)
        end if
        if (is_set_param("aircraft_Tyr2")) then
          call get_param("aircraft_Tyr2",aircraft_Tyr2)
        else
          call stop_model("Must provide aircraft_Tyr2 via rundeck",255)
        end if
        if (aircraft_Tyr1==aircraft_Tyr2) then
          trans_emis=.false.; yr1=0; yr2=0
        else
          trans_emis=.true.; yr1=aircraft_Tyr1; yr2=aircraft_Tyr2
        end if
      end if

      if (nTracer == 0) then
        call stop_model("nTracer undefined in get_aircraft_tracer",255)
      end if

      if (itime < itime_tr0(nTracer)) goto 999 ! returns w/o doing source

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)

      ! Determine year of emissions to use and (for nc emissions)
      ! whether the timestream should be cyclic or not.
      ! (Actual year 'year' has been passed in. Allow override of this
      ! if say, {master,o3,aer_int}_yr non-zero):
      call get_param('master_yr',copy_master_yr,default=0)
      cyclic_yr=copy_master_yr
#ifdef TRACERS_SPECIAL_Shindell
      if ((nTracer>=ntm_chem_beg).and.(nTracer<=ntm_chem_end)) then
        call get_param('o3_yr',cyclic_yr,default=copy_master_yr)
        select case (trname(nTracer))
        case ('NOx')
          call get_param('NOx_yr',cyclic_yr,default=cyclic_yr)
        case ('CO')
          call get_param('CO_yr',cyclic_yr,default=cyclic_yr)
        case ('Alkenes', 'Paraffin')
          call get_param('VOC_yr',cyclic_yr,default=cyclic_yr)
        end select
      else
#endif
#if (defined TRACERS_AEROSOLS_Koch) || (defined TRACERS_AMP) || \
    (defined TRACERS_TOMAS)
        call get_param('aer_int_yr',cyclic_yr,default=copy_master_yr)
        select case (trname(nTracer))
        case ('SO2', 'SO4', 'M_ACC_SU', 'M_AKK_SU', 'ASO4__01')
          call get_param('SO2_int_yr',cyclic_yr,default=cyclic_yr)
        case ('NH3')
          call get_param('NH3_int_yr',cyclic_yr,default=cyclic_yr)
        case ('BCII', 'BCB', 'M_BC1_BC', 'M_BOC_BC', 'AECOB_01')
          call get_param('BC_int_yr',cyclic_yr,default=cyclic_yr)
        case ('OCII', 'OCB', 'M_OCC_OC', 'M_BOC_OC', 'AOCOB_01',
     &        'vbsAm2', 'vbsAm1', 'vbsAz', 'vbsAp1', 'vbsAp2',
     &        'vbsAp3', 'vbsAp4', 'vbsAp5', 'vbsAp6')
          call get_param('OC_int_yr',cyclic_yr,default=cyclic_yr)
        end select
#endif
#ifdef TRACERS_SPECIAL_Shindell
      end if
#endif
      cyclic_yr=abs(cyclic_yr)
      xyear=year
      if (cyclic_yr > 0) xyear=cyclic_yr

      if (isItFbsa) then
        ! for old giss binary files, skip execessive reading by disallowing
        ! emissions before year 1900 (if transient emissions requested):
        if (trans_emis .and. xyear < 1900) goto 999
      end if

      if (need_read) then

        ! Monthly sources are interpolated to the current day
        ! Units are kg m-2 s-1, so no conversion is necessary:

        if (isItFbsa) then
          call openunit(fileName,fileUnit,.true.)
          call read_monthly_3Dsources(Laircr,fileUnit,
     &         src,trans_emis,yr1,yr2,xyear,xday)
          call closeunit(fileUnit)
        else
          if(first_aircraft(nTracer)) then
            call set_first_aircraft(nTracer, .false.)
            call get_param('nc_emis_use_ppm_interp',do_ppm,default=1)
            if (do_ppm==1) then
              call init_stream(grid,AIRCstream,fileName,
     &        trim(trname(nTracer)), 0d0, 1d30, 'ppm',
     &        xyear, xday, cyclic = (cyclic_yr > 0) )
            else
              call init_stream(grid,AIRCstream,fileName,
     &        trim(trname(nTracer)), 0d0, 1d30, 'linm2m',
     &        xyear, xday, cyclic = (cyclic_yr > 0) )
            end if
          end if
          call read_stream(grid,AIRCstream,xyear,xday,src)
        end if

! Place aircraft sources onto model levels:
        airtracer = 0.d0
        do j=J_0,J_1
          do i=I_0,I_1
            zmod(:)=phi(i,j,:)*bygrav*1.d-3 ! km
            do LL=1,Laircr
              if (src(i,j,LL) > 0.) then
                loop_L: do L=1,LM
                  if (zairL(LL) <= zmod(L)) then
                    airtracer(i,j,L) = airtracer(i,j,L) +
     &                                 src(i,j,LL)*axyp(i,j)
                    exit loop_L
                  end if
                  if(L==LM)call stop_model("aircraft lev. problem",255)
                end do loop_L
              end if ! is there a source?
            end do ! LL aircraft levels
          end do ! I
        end do ! J

      end if ! read was needed

      tr3Dsource(I_0:I_1,J_0:J_1,:,nAircraft,nTracer) =
     & airtracer(I_0:I_1,J_0:J_1,:)*get_src_fact(nTracer)

999   continue
      return
      end subroutine get_aircraft_tracer
 

      SUBROUTINE read_monthly_3Dsources
     & (Ldim,iu,data1,trans_emis,yr1,yr2,xyear,xday)
!@sum Read in monthly sources and interpolate to current day
!@auth Jean Lerner and others / Greg Faluvegi
      USE RESOLUTION, only : im,jm
      USE JulianCalendar_mod, only: idofm=>JDmidOfM
      USE FILEMANAGER, only : NAMEUNIT
      USE DOMAIN_DECOMP_ATM, only : GRID,getDomainBounds,READT_PARALLEL
     &     ,REWIND_PARALLEL
     &     ,write_parallel,backspace_parallel,am_i_root
      implicit none
!@var Ldim how many vertical levels in the read-in file?
!@var L dummy vertical loop variable
      integer :: Ldim,L,imon,iu,ipos,k,nn,k2,kstep=10
      character(len=300) :: out_line
      real*8 :: frac, alpha
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::A2D,B2D,dummy
      real*8, DIMENSION(GRID%I_STRT_HALO:GRID%I_STOP_HALO
     *     ,GRID%J_STRT_HALO:GRID%J_STOP_HALO,Ldim) ::tlca,tlcb,data1
     *     ,sfc_a,sfc_b
      logical, intent(in):: trans_emis
      integer, intent(in):: yr1,yr2,xyear,xday
     
      integer :: J_0, J_1, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)     
      call getDomainBounds(grid, I_STRT=I_0, I_STOP=I_1)     

C No doubt this code can be combined/compressed, but I am going to
C do the transient and non-transient cases separately for the moment:

! -------------- non-transient emissions ----------------------------!
      if(.not.trans_emis) then
C
      imon=1                ! imon=January
      if (xday <= 16)  then ! DAY in Jan 1-15, first month is Dec
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo  
        call rewind_parallel(iu)
      else              ! DAY is in Jan 16 to Dec 16, get first month
        do while(xday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        if(imon/=2)then ! avoids advancing records at start of file
          if(am_i_root())write(6,*) 'Not using this first record:'
          call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-2))
        end if
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo   
        if(imon==13) call rewind_parallel(iu)
      end if
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(I_0:I_1,J_0:J_1,L)=B2D(I_0:I_1,J_0:J_1)
      enddo 
c**** Interpolate two months of data to current day
      frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
      data1(I_0:I_1,J_0:J_1,:) =
     & tlca(I_0:I_1,J_0:J_1,:)*frac + tlcb(I_0:I_1,J_0:J_1,:)*(1.-frac)
      write(out_line,*) '3D source monthly factor=',frac
      call write_parallel(trim(out_line))

! --------------- transient emissions -------------------------------!
      else
        ! 3D source files as of now have no meta-data so assume
        ! that transient time slices are decadal:
        kstep=10
        ipos=1
        k2=yr1
        alpha=0.d0 ! before start year, use start year value
        if(xyear>yr2.or.(xyear==yr2.and.xday>=183))then
          alpha=1.d0 ! after end year, use end year value
          ipos=(yr2-yr1)/kstep
          k2=yr2-kstep
        endif
        do k=yr1,yr2-kstep,kstep
          if(xyear>k .or. (xyear==k.and.xday>=183)) then
            if(xyear<k+kstep .or. (xyear==k+kstep.and.xday<183))then
              ipos=1+(k-yr1)/kstep ! (integer artithmatic)
              alpha=real(xyear-k)/real(kstep)
              k2=k
              exit
            endif
          endif
        enddo
!
! read the two necessary months from the first decade:
!
      imon=1                ! imon=January
      if (xday <= 16)  then ! DAY in Jan 1-15, first month is Dec
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! DAY is in Jan 16 to Dec 16, get first month
        do while(xday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        if(imon/=2 .or. ipos/=1)then ! avoids advancing records at start of file
          if(am_i_root())write(6,*) 'Not using this first record:' 
          call readt_parallel
     &    (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        end if
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCC write(6,*) 'Not using this first record:'
CCCCC call readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(I_0:I_1,J_0:J_1,L)=B2D(I_0:I_1,J_0:J_1)
      enddo
      frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
      sfc_a(I_0:I_1,J_0:J_1,:) =
     & tlca(I_0:I_1,J_0:J_1,:)*frac + tlcb(I_0:I_1,J_0:J_1,:)*(1.-frac)
      call rewind_parallel( iu )

      ipos=ipos+1
      imon=1                ! imon=January
      if (xday <= 16)  then ! DAY in Jan 1-15, first month is Dec
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*11)
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        do nn=1,12*Ldim; call backspace_parallel(iu); enddo
      else              ! DAY is in Jan 16 to Dec 16, get first month
        do while(xday > idofm(imon) .AND. imon <= 12)
          imon=imon+1
        enddo
        if(am_i_root())write(6,*) 'Not using this first record:'
        call readt_parallel
     &  (grid,iu,nameunit(iu),dummy,(ipos-1)*12*Ldim+Ldim*(imon-2))
        do L=1,Ldim
          call readt_parallel(grid,iu,nameunit(iu),A2D,1)
          tlca(I_0:I_1,J_0:J_1,L)=A2D(I_0:I_1,J_0:J_1)
        enddo
        if(imon==13)then
          do nn=1,12*Ldim; call backspace_parallel(iu); enddo
        endif
      end if
CCCCCCwrite(6,*) 'Not using this first record:'
CCCCCCcall readt_parallel(grid,iu,nameunit(iu),dummy,Ldim*(imon-1))
      do L=1,Ldim
        call readt_parallel(grid,iu,nameunit(iu),B2D,1)
        tlcb(I_0:I_1,J_0:J_1,L)=B2D(I_0:I_1,J_0:J_1)
      enddo
      frac = float(idofm(imon)-xday)/(idofm(imon)-idofm(imon-1))
      sfc_b(I_0:I_1,J_0:J_1,:) =
     & tlca(I_0:I_1,J_0:J_1,:)*frac + tlcb(I_0:I_1,J_0:J_1,:)*(1.-frac)

! now interpolate between the two time periods:

      data1(I_0:I_1,J_0:J_1,:) = sfc_a(I_0:I_1,J_0:J_1,:)*(1.d0-alpha) 
     & + sfc_b(I_0:I_1,J_0:J_1,:)*alpha

      write(out_line,*) '3D source at',
     &100.d0*alpha,' % of period this day ',k2,' to this day ',k2+kstep,
     &' and monthly fraction= ',frac 
      call write_parallel(trim(out_line))

      endif ! transient or not

      return
      end subroutine read_monthly_3Dsources

#endif /* defined TRACERS_SPECIAL_Shindell or Koch/AMP/TOMAS aerosols */
