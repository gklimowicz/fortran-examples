#include "rundeck_opts.h"

      module gcdiag
      implicit none

      save

!@var JK_xxx, JL_xxx names for AGC indices
!@+   JL/JK refer to model versus constant-pressure levels
      INTEGER ::
     &     JK_dpa ,JK_dpb ,JK_temp, JK_theta ,JK_u
     &    ,JK_v ,JK_eddke ,JK_totke ,JK_eddntsh
     &    ,JK_totntdse ,JK_eddntgeo ,JK_eddntlh
     &    ,JK_totntlh ,JK_totntke ,JK_eddntmom
     &    ,JK_totntmom ,JK_p2kedpgf ,JK_dpsqr ,JK_nptsavg
     &    ,JK_vvel ,JK_eddvtdse ,JK_totvtdse ,JK_eddvtlh
     &    ,JK_totvtlh ,JK_vtgeoeddy ,JK_barekegen ,JK_potvort
     &    ,JK_vtpv ,JK_vtpveddy ,JK_nptsavg1 ,JK_totvtke
     &    ,JK_vtameddy ,JK_totvtam ,JK_sheth ,JK_dudtmadv
     &    ,JK_dtdtmadv ,JK_dudttem ,JK_dtdttem ,JK_epflxncp
     &    ,JK_epflxvcp ,JK_uinst ,JK_totdudt ,JK_tinst
     &    ,JK_totdtdt ,JK_eddvtpt, JK_rhoedz
     &    ,JK_psi, JK_dpb2
      INTEGER :: jk_dudt_econv,
     *  jk_dudt_epdiv,jk_stdev_dp,
     *  jk_dtempdt_econv,jl_phi_amp_wave1,jl_phi_phase_wave1,
     *  jk_vt_se_eddy,jk_tot_vt_se,jk_psi_tem,jk_epflx_v,
     *  jk_nt_eqgpv,jk_dyn_conv_eddy_geop,jk_nt_sheat_e,
     *  jk_dyn_conv_dse,jk_seke,
     *  jk_nt_dse_se,jk_nt_dse_e,
     *  jk_nt_see,jk_tot_nt_se,
     *  jk_nt_am_stand_eddy,
     *  jk_we_flx_nor,jk_we_flx_div,jk_refr_ind_wave1,
     *  jk_del_qgpv,jk_nt_lh_se,jk_wstar,jk_vstar
      INTEGER ::
     &     jl_totntlh,jl_zmfntlh,jl_totvtlh,jl_zmfvtlh,jl_ape,
     &     jl_epflxv,jl_epflxn,jl_40,jl_47,jl_zmfntmom,jl_totntmom,
     &     jl_nt_lh_e,jl_vt_lh_e,jl_dpb,jl_dpsig,jl_dpsigv,jl_epflx_div

      integer :: ! outputs from STRAT_DIAG
     &     jk_dudt_sum1,jk_dudt_meanadv,jk_dudt_eddycnv,
     &     jk_dudt_trnsadv,jk_dudt_epflxdiv,
     &     jk_dudt_fderr1,jk_dudt_fderr2

!@var IJK_xxx AIJK diagnostic names
      INTEGER ::
     &     IJK_UB, IJK_VB, IJK_DSE, IJK_DPB, IJK_TB, IJK_W, IJK_PF,
     &     IJK_UV, IJK_VQ, IJK_VT, IJK_UU, IJK_VV, IJK_TT, IJK_PHI,
     &     IJK_BAREKEGEN


!@var force_gc_vmean a mechanism to force vertical averaging for qtys
!@+   at layer edges or that otherwise lack layer weighting info.
      logical, dimension(:), allocatable :: force_gc_vmean

      end module gcdiag


      subroutine gc_defs
      use CONSTANT, only : twopi,rgas,lhe,bygrav,sha
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : qcheck,dtsrc
      use DYNAMICS, only : dt,do_gwdrag
      use GCDIAG
      use DIAG_COM
      use GC_COM
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use GEOM, only : lat_dg
      use cdl_mod
      use MDIAG_COM, only : make_timeaxis
      use diag_zonal, only : xwon
      implicit none
      integer :: k,kk
      character(len=10) :: ystr,zstr,powstr
      type(cdl_type) :: cdl_dum
      logical :: set_miss
c
      allocate(force_gc_vmean(kagcx))
      force_gc_vmean = .false.
      do k=1,kagcx
         write(sname_gc(k),'(a3,i3.3)') 'AGC',k
         lname_gc(k) = 'unused'
         units_gc(k) = 'unused'
         pow_gc(k) = 0
         scale_gc(k) = 1.
         ia_gc(k) = ia_dga
         denom_gc(k) = 0
         jgrid_gc(k) = 1
         lgrid_gc(k) = ctr_cp
      enddo
c
      k=0
c
      k=k+1
      jk_dpa = k             !'AJK01'
      sname_gc(k) = 'dp_cp1' !   DP=MIN(PM(K),PS)-MIN(PM(K+1),PS)
      lname_gc(k) =  'PRESSURE DIFFERENCES (CP,PT)' ! DP (PT GRID)
      units_gc(k) = 'mb'
      scale_gc(k) = byim
      jgrid_gc(k) = 1
c
      k=k+1
      jk_dpb = k
      sname_gc(k) = 'dp_cp2' !'AJK02'
      lname_gc(k) = 'PRESSURE DIFFERENCES (CP,UV)' ! DP (UV GRID)
      units_gc(k) = 'mb'
      scale_gc(k) = byim
      jgrid_gc(k) = 2
c
      k=k+1
      jl_dpb = k ! only used as a denominator
      sname_gc(k) = 'jl_dpb'
c      lname_gc(k) = 'MASS AT SECONDARY LATITUDES' ! in DIAGB
c      units_gc(k) = 'mb'
      scale_gc(k) = 1.
      jgrid_gc(k) = 2
c
      k=k+1
      jl_dpsig = k ! only used as a denominator
      sname_gc(k) = 'dp_dpsig'
c      lname_gc(k) = 'PRESSURE DIFFERENCES (SIG,PT)'
c      units_gc(k) = 'mb'
      scale_gc(k) = 1
      jgrid_gc(k) = 1
c
      k=k+1
      jl_dpsigv = k ! only used as a denominator
      sname_gc(k) = 'dp_dpsigv'
c      lname_gc(k) = 'PRESSURE DIFFERENCES (SIG,UV)'
c      units_gc(k) = 'mb'
      scale_gc(k) = 1
      jgrid_gc(k) = 2
c
      k=k+1
      jk_rhoedz = k ! only used as a denominator
      sname_gc(k) = 'rhoedz'
      lname_gc(k) = 'LAYER EDGE DENSITY * G * DZ'
      units_gc(k) = 'mb'
      scale_gc(k) = 1
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_nptsavg1 = k
      sname_gc(k) = 'npts_avg1' !'AJK35'
      lname_gc(k) = 'NUMBER OF GRIDPOINTS IN AVERAGE (PT GRID)'
      units_gc(k) = '1'
      jgrid_gc(k) = 1
c
      k=k+1
      jk_nptsavg = k
      sname_gc(k) = 'npts_avg' !'AJK24'
      lname_gc(k) = 'NUMBER OF GRIDPOINTS INCLUDED IN AVERAGE (CP)'
      units_gc(k) = '1'
      scale_gc(k) = XWON !TWOPI/(DLON*FIM)
      jgrid_gc(k) = 2
c
      k=k+1
      jk_dpsqr = k
      sname_gc(k) = 'dp_sqr' !'AJK23'
      lname_gc(k) = 'SQUARE OF PRESSURE THICK. FOR EXISTING UV POINTS'
      units_gc(k) = 'mb*mb'
      jgrid_gc(k) = 2
      denom_gc(k) = jk_nptsavg
c
      k=k+1
      jk_dpb2 = k
      sname_gc(k) = 'dp_cp_existing'
      lname_gc(k) = 'PRESSURE THICK. FOR EXISTING CP,UV POINTS'
      units_gc(k) = 'mb'
      jgrid_gc(k) = 2
      denom_gc(k) = jk_nptsavg
c
      k=k+1
      jk_temp = k
      sname_gc(k) = 'temp' !'AJK03'
      lname_gc(k) = 'TEMPERATURE' !'(TX-TFREZ)*DP'
      units_gc(k) = 'C'
      scale_gc(k) = 1.
      jgrid_gc(k) = 1
      denom_gc(k) = jk_dpa
      pow_gc_vmean(k) = -1
c
      k=k+1
      jk_theta = k
      sname_gc(k) = 'pot_temp' !'AJK06'
      lname_gc(k) = 'POTENTIAL TEMPERATURE' !'TH*DP'
      units_gc(k) = 'K'
      scale_gc(k) = p1000k
      jgrid_gc(k) = 1
      denom_gc(k) = jk_dpa
c
      k=k+1
      jk_potvort = k
      sname_gc(k) = 'pot_vort' !'AJK32'
      lname_gc(k) = 'POTENTIAL VORTICITY (CP)'
      units_gc(k) = 'K/(mb*s)'
      pow_gc(k) = -6
      scale_gc(k) = byim*p1000k
      jgrid_gc(k) = 1
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_u = k
      sname_gc(k) = 'u' !'AJK08'
      lname_gc(k) = 'ZONAL WIND (U COMPONENT)' !'U*DP4  (UV GRID)'
      units_gc(k) = 'm/s' !'100 PA*m/s'
      pow_gc(k) = -1
      scale_gc(k) = 1.
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_v = k
      sname_gc(k) = 'v' !'AJK09'
      lname_gc(k) = 'MERIDIONAL WIND (V COMPONENT)' !'V*DP4  (UV GRID)'
      units_gc(k) = 'm/s' !'100 PA*m/s'
      pow_gc(k) = -2
      pow_gc_vmean(k) = -3
      scale_gc(k) = 1.
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_eddke = k
      sname_gc(k) = 'eddy_ke' !'AJK10'
      lname_gc(k) = 'EDDY KINETIC ENERGY'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .5
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_totke = k
      sname_gc(k) = 'tot_ke' !'AJK11'
      lname_gc(k) = 'TOTAL KINETIC ENERGY'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .5
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_eddntsh = k
      sname_gc(k) = 'nt_sheat_eddy'
      lname_gc(k) = 'NORTH. TRANS. OF SENSIBLE HEAT BY EDDIES'
      units_gc(k) = 'W/mb'
      scale_gc(k) = SHA*XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 11
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_eddntgeo = k
      sname_gc(k) = 'edd_nt_geo' !'AJK14'
      lname_gc(k) = 'NORTH. TRANS. GEOPOT. BY EDDIES'
      units_gc(k) = 'W/mb'
      scale_gc(k) = SHA*XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 11
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_totntdse = k
      sname_gc(k) = 'tot_nt_dse'
      lname_gc(k) = 'TOTAL NORTH. TRANSPORT OF DRY STATIC ENERGY'
      units_gc(k) = 'W/mb'
      scale_gc(k) = XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 12
      pow_gc_vmean(k) = 11
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_eddntlh = k
      sname_gc(k) = 'nt_lh_e'
      lname_gc(k) = 'NORTHWARD TRANSPORT OF LATENT HEAT BY EDDIES'
      units_gc(k) = 'W/mb'
      pow_gc(k) = 10
      scale_gc(k) = lhe*XWON*FIM*1d2*BYGRAV
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_totntlh = k
      sname_gc(k) = 'tot_nt_lh' !'AJK17'
      lname_gc(k) = 'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT'
      units_gc(k) = 'W/mb'
      pow_gc(k) = 10
      scale_gc(k) = LHE*XWON*FIM*100.*BYGRAV
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_totntke = k
      sname_gc(k) = 'tot_nt_ke' !'AJK19'
      lname_gc(k) = 'TOTAL NORTHWARD TRANSPORT OF KINETIC ENERGY'
      units_gc(k) = 'W/mb'
      pow_gc(k) = 9
      scale_gc(k) = .5*XWON*FIM*1d2*BYGRAV
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_eddntmom = k
      sname_gc(k) = 'nt_u_eddy' !'AJK20'
      lname_gc(k) = 'NORTH. TRANS. ZONAL MOM. BY EDDIES'
      units_gc(k) = 'm^2/s^2'
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_totntmom = k
      sname_gc(k) = 'tot_nt_u'
      lname_gc(k) = 'TOTAL NORTH. TRANS. ZONAL MOM.' ! DIAGJK adds uearth term
      units_gc(k) = 'm^2/s^2'
      pow_gc(k) = 1
      jgrid_gc(k) = 2
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_p2kedpgf = k
      sname_gc(k) = 'p2k_eddy_pgf' !'AJK22'
      lname_gc(k) = 'P-K BY EDDY PRESSURE GRADIENT FORCE'
      units_gc(k) = 'W/(m^2*mb)'
      scale_gc(k) = .5*1d2*BYGRAV/DT
      jgrid_gc(k) = 2
      pow_gc(k) = -4
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_vvel = k
      sname_gc(k) = 'vvel' !'AJK25'
      lname_gc(k) = 'VERTICAL VELOCITY (POSITIVE UPWARD)'
      units_gc(k) = 'mb/s'
      pow_gc(k) = -5
      scale_gc(k) = -byim
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_eddvtdse = k
      sname_gc(k) = 'vt_dse_e'
      lname_gc(k) = 'VERT. TRANS. OF DRY STATIC ENERGY BY EDDIES (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      pow_gc(k) = -1
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_totvtdse = k
      sname_gc(k) = 'tot_vt_dse' !'AJK27'
      lname_gc(k) = 'TOTAL LGE SCALE VERT. TRANS. OF DRY STAT ENRG (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      pow_gc(k) = 1
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_eddvtlh = k
      sname_gc(k) = 'vt_lh_eddy'
      lname_gc(k) = 'VERTICAL TRANSPORT OF LATENT HEAT BY EDDIES (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM*LHE
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_totvtlh = k
      sname_gc(k) = 'tot_vt_lh' !'AJK29'
      lname_gc(k) = 'TOTAL LGE SCALE VERT. TRANS. OF LATENT HEAT (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM*LHE
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_vtgeoeddy = k
      sname_gc(k) = 'vt_geopot_eddy' !'AJK30'
      lname_gc(k) = 'VERT. TRANS. OF GEOPOT. ENERGY BY EDDIES (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      pow_gc(k) = -2
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_barekegen = k
      sname_gc(k) = 'baroc_eddy_ke_gen' !'AJK31'
      lname_gc(k) = 'BAROCLINIC EDDY KINETIC ENERGY GEN.'
      units_gc(k) = 'W/(m^2*mb)'
      scale_gc(k) = .5*RGAS*1d2*BYGRAV
      pow_gc(k) = -4
      jgrid_gc(k) = 1
      lgrid_gc(k) = ctr_cp ! according to print routine
      denom_gc(k) = jk_dpb
c
      k=k+1
      jk_vtpv = k
      sname_gc(k) = 'vt_pv' !'AJK33'
      lname_gc(k) = 'VERT. TRANSPORT OF POTENTIAL VORTICITY (CP)'
      units_gc(k) = 'K/s^2'
      pow_gc(k) = -9
      scale_gc(k) = -P1000K
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_vtpveddy = k
      sname_gc(k) = 'vt_pv_eddy' !'AJK34'
      lname_gc(k) = 'VERT. TRANS. OF POT. VORT. BY EDDIES (CP)'
      units_gc(k) = 'K/S^2'
      pow_gc(k) = -9
      scale_gc(k) = -P1000K
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_totvtke = k
      sname_gc(k) = 'tot_vt_ke' !'AJK36'
      lname_gc(k) ='TOTAL LGE SCALE VERT. TRANS. OF KINETIC ENRG (CP)'
      units_gc(k) = 'W/m^2'
      pow_gc(k) = -1
      scale_gc(k) = -.5*100.*BYGRAV*BYIM
      jgrid_gc(k) = 2
      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_vtameddy = k
      sname_gc(k) = 'vt_u_eddy' !'AJK37'
      lname_gc(k) = 'EDDY VERTICAL ZONAL MOM. FLUX (CP)'
      units_gc(k) = 'm^2/s^2'
      pow_gc(k) = -3
      scale_gc(k) = -100.*BYGRAV*BYIM
      jgrid_gc(k) = 2
      lgrid_gc(k) = edg_cp
      denom_gc(k) = jk_rhoedz
c
      k=k+1
      jk_totvtam = k
      sname_gc(k) = 'tot_vt_u' !'AJK38'
      lname_gc(k) = 'TOTAL VERTICAL ZONAL MOM. FLUX (CP)'
      units_gc(k) = 'm^2/s^2'
      pow_gc(k) = -2
      scale_gc(k) = -100.*BYGRAV*BYIM
      jgrid_gc(k) = 2
      lgrid_gc(k) = edg_cp
      denom_gc(k) = jk_rhoedz
c
      k=k+1
      jk_sheth = k
      sname_gc(k) = 'sheth' !'AJK39'
c      lname_gc(k) = 'unknown'
c      units_gc(k) = 'unknown'
c
      k=k+1
      jk_dudtmadv = k
      sname_gc(k) = 'dudt_mean_advec' !'AJK40'
      lname_gc(k) = 'DU/DT BY MEAN ADVECTION (CP)'
      units_gc(k) = 'm/s^2'
      pow_gc(k) = -6
      scale_gc(k) = 1.
      jgrid_gc(k) = 2
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_dtdtmadv = k
      sname_gc(k) = 'dtempdt_mean_advec' !'AJK41'
      lname_gc(k) = 'DTEMP/DT BY MEAN ADVECTION (CP)'
      units_gc(k) = 'K/DAY'
      pow_gc(k) = -1
      scale_gc(k) = SECONDS_PER_DAY
      jgrid_gc(k) = 1
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_dudttem = k
      sname_gc(k) = 'dudt_advec_tem' !'AJK42'
      lname_gc(k) = 'DU/DT BY TRANSFORMED ADVECTION (CP)'
      units_gc(k) = 'm/s^2'
      pow_gc(k) = -6
      scale_gc(k) = 1.
      jgrid_gc(k) = 2
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_dtdttem = k
      sname_gc(k) = 'dtempdt_advec_tem' !'AJK43'
      lname_gc(k) = 'DTEMP/DT BY TRANSFORMED ADVECTION (CP)'
      units_gc(k) = 'K/DAY'
      pow_gc(k) = -1
      scale_gc(k) = SECONDS_PER_DAY
      jgrid_gc(k) = 1
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_epflxncp = k
      sname_gc(k) = 'epflx_north_cp' !'AJK44'
c      lname_gc(k) = 'NORTHWARD COMP. OF ELIASSEN-PALM FLUX (CP)'
c      units_gc(k) = 'unknown'
c      jgrid_gc(k) = 2
c
      k=k+1
      jk_epflxvcp = k
      sname_gc(k) = 'epflx_vert_cp' !'AJK45'
c      lname_gc(k) = 'VERTICAL COMP. OF ELIASSEN-PALM FLUX (CP)'
c      units_gc(k) = 'unknown'
c      jgrid_gc(k) = 1
c      lgrid_gc(k) = edg_cp
c
      k=k+1
      jk_uinst = k
      sname_gc(k) = 'u_inst' !'AJK46'
      lname_gc(k) = 'INSTANTANEOUS ZONAL AVERAGE OF ZONAL WIND'
      units_gc(k) = 'm/s'
      jgrid_gc(k) = 1
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_totdudt = k
      sname_gc(k) = 'tot_dudt' !'AJK47'
      lname_gc(k) = 'DU/DT   TOTAL CHANGE (CP)'
      units_gc(k) = 'm/s^2'
      pow_gc(k) = -6
      jgrid_gc(k) = 2
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_tinst = k
      sname_gc(k) = 't_inst' !'AJK48'
      lname_gc(k) = 'INSTANTANEOUS ZONAL AVERAGE OF TEMPERATURE'
      units_gc(k) = 'K'
      jgrid_gc(k) = 1
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_totdtdt = k
      sname_gc(k) = 'dtempdt' !'AJK49'
      lname_gc(k) = 'DTEMP/DT   TOTAL CHANGE (CP)'
      units_gc(k) = 'K/DAY'
      pow_gc(k) = -1
      jgrid_gc(k) = 1
      scale_gc(k) = SECONDS_PER_DAY
      force_gc_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jk_eddvtpt = k
      sname_gc(k) = 'edd_vt_pt'
c      lname_gc(k) = 'EDDY VERTICAL TRANSPORT OF POT. TEMP.'
c      units_gc(k) = 'unknown'
c      jgrid_gc(k) = 1
c      lgrid_gc(k) = edg_cp
c
      k=k+1
      jl_ape = k
      sname_gc(k) = 'avail_pe' !
      lname_gc(k) = 'AVAILABLE POTENTIAL ENERGY'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .5*RGAS
      ia_gc(k) = ia_d5s
      jgrid_gc(k) = 1
      denom_gc(k) = jl_dpsig
c
      k=k+1
      jl_epflxv = k
      sname_gc(k) = 'epflx_vert' !
      lname_gc(k) = 'VERTICAL ELIASSEN-PALM FLUX'
      units_gc(k) = 'm^2/s^2'
      pow_gc(k) = -3
      scale_gc(k) = .5*100.*BYGRAV*BYIM
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_cp
      denom_gc(k) = jk_rhoedz ! should be jl_rhoe instead
c
      k=k+1
      jl_epflxn = k
      sname_gc(k) = 'epflx_north' !
      lname_gc(k) = 'NORTHWARD ELIASSEN-PALM FLUX'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = 1.
      jgrid_gc(k) = 2
      denom_gc(k) = jl_dpsigv ! why not jl_dpb
c
      k=k+1
      jl_zmfntmom = k
      sname_gc(k) = 'zmf_nt_mom' !
c      lname_gc(k) = 'NORTH TRANS ZON. MOM. BY ZON. MEAN FLOW'
c      units_gc(k) = 'unknown'
c
      k=k+1
      jl_totntmom = k
      sname_gc(k) = 'tot_nt_mom' !
c      lname_gc(k) = 'TOTAL NORTH TRANS ZON. MOM.'
c      units_gc(k) = 'unknown'
c
      k=k+1
      jl_zmfntlh = k
      sname_gc(k) = 'jl_zmf_nt_lh'    ! used in DIAGJK but not printed
      lname_gc(k) = 'MEAN MERIDIONAL NORTHWARD TRANS. OF LATENT HEAT'
      units_gc(k) = 'W/mb'
      pow_gc(k) = 9
      scale_gc(k) = 100.*bygrav*LHE*XWON*fim/DTsrc
      ia_gc(k) = ia_src
      jgrid_gc(k) = 2
      denom_gc(k) = jl_dpb
c
      k=k+1
      jl_totntlh = k
      sname_gc(k) = 'jl_tot_nt_lh'
      lname_gc(k) = 'TOTAL NORTHWARD TRANSPORT OF LATENT HEAT (QDYN)'
      units_gc(k) = 'W/mb'
      pow_gc(k) = 10
      scale_gc(k) = 100.*bygrav*LHE*XWON*fim/DTsrc
      ia_gc(k) = ia_src
      jgrid_gc(k) = 2
      denom_gc(k) = jl_dpb
c
      k=k+1
      jl_zmfvtlh = k     ! used in DIAGJK but not printed
      sname_gc(k) = 'jl_zmf_vt_lh'
      lname_gc(k) = 'MEAN MERIDIONAL VERTICAL TRANS. OF LATENT HEAT'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = 100.*BYGRAV*LHE*XWON*byim/DTsrc
      ia_gc(k) = ia_src
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_ml
c
      k=k+1
      jl_totvtlh = k
      sname_gc(k) = 'jl_tot_vt_lh'
      lname_gc(k) = 'TOTAL VERTICAL TRANSPORT OF LATENT HEAT'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = 100.*BYGRAV*LHE*XWON*byim/DTsrc
      ia_gc(k) = ia_src
      jgrid_gc(k) = 1
      lgrid_gc(k) = edg_ml
c
      k=k+1
      jl_47 = k
      sname_gc(k) = 'AJL47' !V-V*  =D((V-VI)*(T-TI)/DTHDP)/DP
      lname_gc(k) = 'unknown'
      units_gc(k) = 'unknown'

c
c derived JKs
c
      k=k+1
      jk_psi = k
      sname_gc(k) = 'psi_cp'
      lname_gc(k) = 'STREAM FUNCTION (CP)'
      units_gc(k) = 'kg/s'
      scale_gc(k) = 100.*BYGRAV
      pow_gc(k) = 9
      jgrid_gc(k) = 2
      lgrid_gc(k) = edg_cp
      k = k + 1
      jl_epflx_div = k                        ; jgrid_gc(k) = 1
      sname_gc(k) = 'epflx_div'
      lname_gc(k) = 'DIVERGENCE OF THE ELIASSEN-PALM FLUX'
      units_gc(k) = 'm/s^2'
      scale_gc(k) = 1.
      pow_gc(k) = -6
      denom_gc(k) = jl_dpsig
      k = k + 1
      jl_nt_lh_e = k
      sname_gc(k) = 'nt_lh_eddy'        ; jgrid_gc(k) = 2
      lname_gc(k) = 'N. TRANSPORT OF LATENT HEAT BY EDDIES (QDYN)'
      units_gc(k) = 'W/mb'
      scale_gc(k) = 100.*bygrav*xwon*lhe*fim/DTsrc
      pow_gc(k) = 10
      ia_gc(k) = ia_src
      denom_gc(k) = jl_dpb
      k = k + 1
      jl_vt_lh_e = k
      sname_gc(k) = 'vt_lh_eddy1'        ; jgrid_gc(k) = 1
      lname_gc(k) = 'V. TRANSPORT OF LATENT HEAT BY EDDIES'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = 100.*bygrav*xwon*lhe*byim/DTsrc
      pow_gc(k) = 0
      ia_gc(k) = ia_src
      lgrid_gc(k) = edg_ml
c
      k = k + 1
      jk_dudt_econv = k                       ; jgrid_gc(k) = jgrid_u
      sname_gc(k) = 'dudt_eddy_conv'
      lname_gc(k) = 'DU/DT BY EDDY CONVERGENCE (CP)'
      units_gc(k) = '10**-6 m/s^2'
      scale_gc(k)= 1.D6
      denom_gc(k) = jk_dpa
      k = k + 1
      jk_dudt_epdiv = k                       ; jgrid_gc(k) = jgrid_u
      sname_gc(k) = 'dudt_epdiv'
      lname_gc(k) = 'DU/DT BY ELIASSEN-PALM DIVERGENCE (CP)'
      units_gc(k) = 'm/s^2'
      pow_gc(k) = -6
      scale_gc(k) = 1.
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dtempdt_econv = k                    ; jgrid_gc(k) = 1
      sname_gc(k) = 'dtempdt_eddy_conv'
      lname_gc(k) = 'DTEMP/DT BY EDDY CONVERGENCE (CP)'
      units_gc(k) = 'K/DAY'
      scale_gc(k) = SECONDS_PER_DAY
      pow_gc(k) = -1
      denom_gc(k) = jk_dpa
      k = k + 1
      jk_vt_se_eddy = k                       ; jgrid_gc(k) = 1
      sname_gc(k) = 'vt_se_eddy'
      lname_gc(k) ='VERTICAL TRANSPORT OF STATIC ENERGY BY EDDIES (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      ia_gc(k) = ia_dga
      lgrid_gc(k) = edg_cp
      k = k + 1
      jk_tot_vt_se = k                        ; jgrid_gc(k) = 1
      sname_gc(k) = 'tot_vt_se'
      lname_gc(k) =
     &    'TOTAL LARGE SCALE VERT. TRANS. OF STATIC ENRG (CP)'
      units_gc(k) = 'W/m^2'
      scale_gc(k) = -100.*BYGRAV*BYIM
      pow_gc(k) = 1
      ia_gc(k) = ia_dga
      lgrid_gc(k) = edg_cp
      k = k + 1
      jk_psi_tem = k                          ; jgrid_gc(k) = 2
      sname_gc(k) = 'psi_tem'
      lname_gc(k) = 'TRANSFORMED STREAM FUNCTION (CP)'
      units_gc(k) = 'kg/s'
      scale_gc(k) = 100.*BYGRAV*XWON
      ia_gc(k) = ia_dga
      pow_gc(k) = 9
      lgrid_gc(k) = edg_cp
c      k = k + 1
c      jk_epflx_v = k                          ; jgrid_gc(k) = 1
c      sname_gc(k) = 'epflx_vert_cp'
c      lname_gc(k) = 'VERTICAL ELIASSEN-PALM FLUX (CP)'
c      units_gc(k) = 'm^2/s^2'
c      scale_gc(k) = -100.*BYGRAV*BYIM
c      ia_gc(k) = ia_dga
c      pow_gc(k) = -2
      k = k + 1
      jk_nt_eqgpv = k                         ; jgrid_gc(k) = 1
      sname_gc(k) = 'nt_eddy_qgpv'
      lname_gc(k) = 'NORTH. TRANS. OF EDDY Q-G POT. VORTICITY'
      units_gc(k) = '10**-6 m/s^2'
      scale_gc(k) = 1.D6
      denom_gc(k) = jk_dpa
      k = k + 1
      jk_dyn_conv_eddy_geop = k               ; jgrid_gc(k) = 1
      sname_gc(k) = 'dyn_conv_eddy_geop'
      lname_gc(k) = 'DYNAMIC CONVERGENCE OF EDDY GEOPOTENTIAL'
      units_gc(k) = 'W/(m^2*mb)'
      scale_gc(k) = 1d2*BYGRAV
      pow_gc(k) = -4
      denom_gc(k) = jk_dpa
      k = k + 1
      jk_dyn_conv_dse = k                     ; jgrid_gc(k) = 1
      sname_gc(k) = 'dyn_conv_dse'
      lname_gc(k) = 'DYNAMIC CONVERGENCE OF DRY STATIC ENERGY'
      units_gc(k) = 'W/(m^2*mb)'
      scale_gc(k) = 1d2*BYGRAV
      pow_gc(k) = -2
      denom_gc(k) = jk_dpa
      k = k + 1
      jk_nt_dse_e = k                         ; jgrid_gc(k) = 2
      sname_gc(k) = 'nt_dse_eddy'
      lname_gc(k) = 'NORTH. TRANS. OF DRY STATIC ENERGY BY EDDIES'
      units_gc(k) = 'W/mb'
      scale_gc(k) = XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 11
      denom_gc(k) = jk_dpb
      k = k + 1
      jk_nt_see = k                           ; jgrid_gc(k) = 2
      sname_gc(k) = 'nt_se_eddy'
      lname_gc(k) = 'NORTH.TRANSPORT OF STATIC ENERGY BY EDDIES'
      units_gc(k) = 'W/mb'
      scale_gc(k) = XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 11
      denom_gc(k) = jk_dpb
      k = k + 1
      jk_tot_nt_se = k                        ; jgrid_gc(k) = 2
      sname_gc(k) = 'tot_nt_se'
      lname_gc(k) = 'TOTAL NORTHWARD TRANSPORT OF STATIC ENERGY'
      units_gc(k) = 'W/mb'
      scale_gc(k) = XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 12
      pow_gc_vmean(k) = 11
      denom_gc(k) = jk_dpb
      k = k + 1
      jk_we_flx_nor = k                       ; jgrid_gc(k) = 2
      sname_gc(k) = 'we_flx_nor'
      lname_gc(k) = 'NORTHWARD WAVE ENERGY FLUX'
c      units_gc(k) = '10**11 JOULES/METER/UNIT SIGMA'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .25
      k = k + 1
      jk_we_flx_div = k                       ; jgrid_gc(k) = 1
      sname_gc(k) = 'we_flx_div'
      lname_gc(k) = 'DIVERGENCE OF THE WAVE ENERGY FLUX'
      units_gc(k) = 'm/s^2'
      scale_gc(k) = 1
      pow_gc(k) = -6
      k = k + 1
      jk_del_qgpv = k                         ; jgrid_gc(k) = 2
      sname_gc(k) = 'del_qgpv'
      lname_gc(k) = 'Q-G POT. VORTICITY CHANGE OVER LATITUDES'
      units_gc(k) = '1/(m*s)'
      scale_gc(k) = 1.
      pow_gc(k) = -12
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_wstar = k                            ; jgrid_gc(k) = 1
      sname_gc(k) = 'wstar'
      lname_gc(k) = 'W*    RESIDUAL VERTICAL VELOCITY'
      units_gc(k) = 'mb/s'
      pow_gc(k) = -5
      scale_gc(k) = 1
      lgrid_gc(k) = edg_cp
      k = k + 1
      jk_vstar = k                            ; jgrid_gc(k) = 2
      sname_gc(k) = 'vstar'
      lname_gc(k) = 'V* = V - D(V''TH''/DTHDP)/DP'
      units_gc(k) = 'm/s'
      pow_gc(k) = -2
      denom_gc(k) = jk_dpb

      if(kep.gt.0) then ! outputs from STRAT_DIAG
      k = k + 1
      jk_dudt_sum1 = k
      sname_gc(k) = 'dudt_sum1'
      lname_gc(k) = 'DU/DT BY EULER CIRC. + CONVEC + DRAG+DIF+ER2'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dudt_meanadv = k
      sname_gc(k) = 'dudt_meanadv'
      lname_gc(k) = 'DU/DT BY MEAN ADVECTION'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dudt_eddycnv = k
      sname_gc(k) = 'dudt_eddycnv'
      lname_gc(k) = 'DU/DT BY EDDY CONVERGENCE'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dudt_trnsadv = k
      sname_gc(k) = 'dudt_trnsadv'
      lname_gc(k) = 'DU/DT BY TRANSFORMED ADVECTION'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dudt_epflxdiv = k
      sname_gc(k) = 'dudt_epflxdiv'
      lname_gc(k) = 'DU/DT BY ELIASSEN-PALM DIVERGENCE'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dudt_fderr1 = k
      sname_gc(k) = 'dudt_fderr1'
      lname_gc(k) = 'DU/DT BY F.D. ERROR TERM 1'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_dudt_fderr2 = k
      sname_gc(k) = 'dudt_fderr2'
      lname_gc(k) = 'DU/DT BY F.D. ERROR TERM 2'
      units_gc(k) = 'm/s^2'
      jgrid_gc(k) = 2
      pow_gc(k) = -6
      force_gc_vmean(k) = .true. ! b/c no denom
      endif ! kep.gt.0

      if(k.gt.kagc) then
        if(am_i_root()) then
          write(6,*) 'gc_defs: Increase kagc=',kagc,' to at least ',k
        endif
        call stop_model( 'kagc too small', 255 )
      end if

c
c more derived JK-arrays
c
      k = KAGC

      k = k + 1
      jk_seke = k                             ; jgrid_gc(k) = jgrid_ke
      sname_gc(k) = 'stand_eddy_ke'
      lname_gc(k) = 'STANDING EDDY KINETIC ENERGY'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = .5
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_nt_dse_se = k                        ; jgrid_gc(k) = 2
      sname_gc(k) = 'nt_dse_stand_eddy'
      lname_gc(k) = 'NOR. TRANS. OF DRY STAT. ENERGY BY STAND. EDDIES'
      units_gc(k) = 'W/mb'
      scale_gc(k) = XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 11
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_nt_lh_se = k
      sname_gc(k) = 'nt_lh_stand_eddy'        ; jgrid_gc(k) = 2
      lname_gc(k) = 'N. TRANSPORT OF LATENT HEAT BY STAND. EDDIES'
      units_gc(k) = 'W/mb'
      scale_gc(k) = lhe*XWON*FIM*1d2*BYGRAV
      pow_gc(k) = 9
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_nt_am_stand_eddy = k                 ; jgrid_gc(k) = 2
      sname_gc(k) = 'nt_u_stand_eddy'
      lname_gc(k) = 'NORTH. TRANS. ZONAL MOM. BY STAND. EDDIES'
      units_gc(k) = 'm^2/s^2'
      scale_gc(k) = 1.
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jk_refr_ind_wave1 = k  !!!!! Refraction Inicies must be in order
      jgrid_gc(k) = 2
      sname_gc(k) = 'refr_ind_wave1'
      lname_gc(k) = 'REFRACTION INDEX FOR WAVE NUMBER 1'
      units_gc(k) = '10**-8 m^-2'
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jgrid_gc(k) = 2
      sname_gc(k) = 'refr_ind_wave2'
      lname_gc(k) = 'REFRACTION INDEX FOR WAVE NUMBER 2'
      units_gc(k) = '10**-8 m^-2'
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jgrid_gc(k) = 2
      sname_gc(k) = 'refr_ind_wave3'
      lname_gc(k) = 'REFRACTION INDEX FOR WAVE NUMBER 3'
      units_gc(k) = '10**-8 m^-2'
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jgrid_gc(k) = 2
      sname_gc(k) = 'refr_ind_wave6'
      lname_gc(k) = 'REFRACTION INDEX FOR WAVE NUMBER 6'
      units_gc(k) = '10**-8 m^-2'
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jgrid_gc(k) = 2
      sname_gc(k) = 'refr_ind_wave9'
      lname_gc(k) = 'REFRACTION INDEX FOR WAVE NUMBER 9'
      units_gc(k) = '10**-8 m^-2'
      force_gc_vmean(k) = .true. ! b/c no denom
      k = k + 1
      jl_phi_amp_wave1 = k                    ; jgrid_gc(k) = 1
      sname_gc(k) = 'phi_amp_wave1'
      lname_gc(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units_gc(k) = 'METERS'
      k = k + 1
      sname_gc(k) = 'phi_amp_wave2'           ; jgrid_gc(k) = 1
      lname_gc(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units_gc(k) = 'METERS'
      k = k + 1
      sname_gc(k) = 'phi_amp_wave3'           ; jgrid_gc(k) = 1
      lname_gc(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units_gc(k) = 'METERS'
      k = k + 1
      sname_gc(k) = 'phi_amp_wave4'           ; jgrid_gc(k) = 1
      lname_gc(k) ='AMPLITUDE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units_gc(k) = 'METERS'
      k = k + 1
      jl_phi_phase_wave1 = k                  ; jgrid_gc(k) = 1
      sname_gc(k) = 'phi_phase_wave1'
      lname_gc(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 1'
      units_gc(k) = 'DEG WEST LONG'
      k = k + 1
      sname_gc(k) = 'phi_phase_wave2'         ; jgrid_gc(k) = 1
      lname_gc(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 2'
      units_gc(k) = 'DEG WEST LONG'
      k = k + 1
      sname_gc(k) = 'phi_phase_wave3'         ; jgrid_gc(k) = 1
      lname_gc(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 3'
      units_gc(k) = 'DEG WEST LONG'
      k = k + 1
      sname_gc(k) = 'phi_phase_wave4'         ; jgrid_gc(k) = 1
      lname_gc(k) = 'PHASE OF GEOPOTENTIAL HEIGHT FOR WAVE NUMBER 4'
      units_gc(k) = 'DEG WEST LONG'
      k = k + 1
      jk_stdev_dp = k                         ; jgrid_gc(k) = 2
      sname_gc(k) = 'stdev_dp'
      lname_gc(k) = 'STANDARD DEVIATION OF PRESSURE DIFFERENCES'
      units_gc(k) = 'MB'
      scale_gc(k) = 1.
c Check the count
      if (k .gt. KAGCx) then
        write (6,*) 'Increase KAGCx=',KAGCx,' to at least ',k
        call stop_model('JK_TITLES: KAGCx too small',255)
      end if

c Print out long names of diagnostic fields if QCHECK=.TRUE.
      if (AM_I_ROOT()) then
         write (6,*) 'Number of AGC diagnostics defined: kagcmax=',k
         if(qcheck) then
           do kk=1,k
             write (6,'(i4,'':'',a)') kk,trim(lname_gc(kk))
           end do
         endif
      end if

      do k=1,kagcx
        if(denom_gc(k).ne.0) cycle ! already set
        if(lgrid_gc(k).eq.edg_cp) force_gc_vmean(k)=.true.
        if(lgrid_gc(k).eq.edg_ml) force_gc_vmean(k)=.true.
      enddo

      lat_gc(:) = lat_dg(:,1)
      lat_gc2(:) = lat_dg(:,2)

c
c Declare the dimensions and metadata of AGC output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      call init_cdl_type('cdl_gc',cdl_gc)
      call add_coord(cdl_gc,'lat',jmlat,units='degrees_north',
     &     coordvalues=lat_gc)
      call add_coord(cdl_gc,'lat2',jmlat,units='degrees_north',
     &     coordvalues=lat_gc2)
      call add_coord(cdl_gc,'pgz',kgz,units='mb',
     &     coordvalues=pmb(1:kgz))
      call add_dim(cdl_gc,'shnhgm',3)
      call add_dim(cdl_gc,'lat_plus3',jmlat+3)
      call add_dim(cdl_gc,'lat2_plus3',jmlat+3)
      cdl_dum = cdl_gc
      call merge_cdl(cdl_dum,cdl_heights,cdl_gc)

      do k=1,kagcx
        if(trim(units_gc(k)).eq.'unused') cycle
        if(jgrid_gc(k).eq.1) then
          ystr='lat'
        else
          ystr='lat2'
        endif
        if(sname_gc(k)(1:7).eq.'phi_amp' .or.
     &     sname_gc(k)(1:9).eq.'phi_phase') then
          zstr='(pgz,'
        elseif(lgrid_gc(k).eq.ctr_ml .or. lgrid_gc(k).eq.ctr_cp) then
          zstr='(plm,'
        else
          zstr='(ple,'
        endif
        set_miss = denom_gc(k).ne.0
        call add_var(cdl_gc,
     &       'float '//trim(sname_gc(k))//trim(zstr)//trim(ystr)//') ;',
     &       units=trim(units_gc(k)),
     &       long_name=trim(lname_gc(k)),
     &       auxvar_string=
     &          'float '//trim(sname_gc(k))//'_hemis'//
     &          trim(zstr)//'shnhgm) ;',
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis
     &       )
        if(pow_gc(k).ne.0) then
          write(powstr,'(i3)') pow_gc(k)
          call add_varline(cdl_gc,
     &         trim(sname_gc(k))//':prtpow = '//trim(powstr)//' ;')
        endif
        if(pow_gc_vmean(k).ne.0) then
          write(powstr,'(i3)') pow_gc_vmean(k)
          call add_varline(cdl_gc,
     &        trim(sname_gc(k))//':prtpow_vmean = '//trim(powstr)//' ;')
        endif
        if(denom_gc(k).gt.0 .or. force_gc_vmean(k)) then
          if(make_timeaxis) then ! hacky logic
            call add_var(cdl_gc,'float '//trim(sname_gc(k))//
     &           '_vmean(time,'//trim(ystr)//'_plus3) ;')
          else
            call add_var(cdl_gc,'float '//trim(sname_gc(k))//
     &           '_vmean('//trim(ystr)//'_plus3) ;')
          endif
        endif
      enddo


! exceptions to forced vertical averaging for layer-edge quantities
      force_gc_vmean(jk_rhoedz) = .false. ! this is a denom: want sum instead

      return
      end subroutine gc_defs

      subroutine ijk_defs
      use CONSTANT, only : bygrav,tf,sha,rgas
      use RESOLUTION, only : im,jm,lm
      use MODEL_COM, only : qcheck,dtsrc
      USE GCDIAG
      use DIAG_COM, only : ia_dga,kaijk,name_ijk,scale_ijk,
     &     lname_ijk,units_ijk,jgrid_ijk,off_ijk,denom_ijk,ia_ijk,
     &     igridc,igride,jgridc,jgride,kgridc,kgride,ijkgridc
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      USE GEOM, only : lon_dg,lat_dg
      use MDIAG_COM, only : make_timeaxis
      use DIAG_COM, only : cdl_ijk,cdl_heights
      use cdl_mod
      implicit none
      integer :: k,kk
      type(cdl_type) :: cdl_dum
c igrid,jgrid,kgrid = 1 for centers, 2 for edges
      integer :: igrid,jgrid,kgrid,ijkg
      character(len=8) :: xstr,ystr,zstr
      logical :: set_miss
c
      do k=1,kaijk
         write(name_ijk(k),'(a4,i3.3)') 'AIJK',k
         lname_ijk(k) = 'unused'
         units_ijk(k) = 'unused'
         scale_ijk(k) = 1.
         jgrid_ijk(k) = igride+jgride+kgridc
         off_ijk(k)   = 0.
         denom_ijk(k) = 0
         ia_ijk(k) = ia_dga
      enddo
c
      k=0
c
      k=k+1
      IJK_DPB=k
      name_ijk(k) = 'dpb' !'DPB'
      lname_ijk(k) = 'DELTA-P'!           b-grid'
      units_ijk(k) = '100 PA'
c
      k=k+1
      IJK_UB=k
      name_ijk(k) = 'ub' !'UDPB'
      lname_ijk(k) = 'U-WIND'!            x delta p, b-grid'
      units_ijk(k) = 'm/s'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_VB=k
      name_ijk(k) = 'vb' !'VDPB'
      lname_ijk(k) = 'V-WIND'!            x delta p, b-grid'
      units_ijk(k) = 'm/s'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_DSE=k
      name_ijk(k) = 'dse' !'DSEDPB'
      lname_ijk(k) = 'DRY STAT. ENERGY'!  x delta p x 4, b-grid'
      units_ijk(k) = 'm^2/s^2'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_TB=k
      name_ijk(k) = 'tb' !'TDPB'
      lname_ijk(k) = 'TEMPERATURE'!       x delta p x 4, b-grid'
      units_ijk(k) = 'C'
      off_ijk(k)   = -TF
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_W=k
      name_ijk(k) = 'w' !'WDPB'
      lname_ijk(k) = 'OMEGA'!  a-grid
      units_ijk(k) = 'Pa/s'
      scale_ijk(k) = 100.
      jgrid_ijk(k) = igridc+jgridc+kgridc ! eventually use kgride instead
c
      k=k+1
      IJK_PF=k
      name_ijk(k) = 'pf'
      lname_ijk(k) = 'FRAC ABOVE SURF'
      units_ijk(k) = '%'
      scale_ijk(k) = 100.0
c
      k=k+1
      IJK_UV=k
      name_ijk(k) = 'UxV'
      lname_ijk(k) = 'U-WIND*V-WIND'
      units_ijk(k) = 'm**2/s**2'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_VQ=k
      name_ijk(k) = 'VxQ'
      lname_ijk(k) = 'V-WIND*SPECIFIC HUMIDITY'
      units_ijk(k) = '(m/s)*kg/kg'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_VT=k
      name_ijk(k) = 'VxT'
      lname_ijk(k) = 'V-WIND*TEMPERATURE'
      units_ijk(k) = '(m/s)*K'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_UU=k
      name_ijk(k) = 'UxU'
      lname_ijk(k) = 'U-WIND*U-WIND'
      units_ijk(k) = 'm**2/s**2'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_VV=k
      name_ijk(k) = 'VxV'
      lname_ijk(k) = 'V-WIND*V-WIND'
      units_ijk(k) = 'm**2/s**2'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_TT=k
      name_ijk(k) = 'TxT'
      lname_ijk(k) = 'TEMPERATURE*TEMPERATURE'
      units_ijk(k) = 'K**2'
      denom_ijk(k) = IJK_DPB
c
      k=k+1
      IJK_BAREKEGEN = k ! JK version is average of this over longitude.
      name_ijk(k) = 'baroc_eddy_ke_gen' ! Caveat Interpretor.
      lname_ijk(k) = 'BAROCLINIC EDDY KINETIC ENERGY GEN.'
      units_ijk(k) = 'W/m^2'
      scale_ijk(k) = .5*RGAS*1d2*BYGRAV
      jgrid_ijk(k) = ijkgridc
!      denom_ijk(k) = IJK_DP ! for res-independent W/(mb*m^2) units
c
c composite outputs
c
      k = k + 1
      IJK_PHI = k
      name_ijk(k) = 'z'
      lname_ijk(k) = 'HEIGHT'
      units_ijk(k) = 'm'
      scale_ijk(k) = BYGRAV
      denom_ijk(k) = IJK_DPB

      if (AM_I_ROOT()) then
         if (k .gt. kaijk) then
            write (6,*) 'ijk_defs: Increase kaijk=',kaijk,' to ',k
            call stop_model( 'kaijk too small', 255 )
         end if
         write (6,*) 'Number of AIJK diagnostics defined: kaijkmax=',k
         if(qcheck) then
           do kk=1,k
             write (6,'(i4,'':'',a)') kk,trim(lname_ijk(kk))
           end do
         endif
      end if

c
c Declare the dimensions and metadata of AIJK output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      call init_cdl_type('cdl_ijk',cdl_ijk)
      call add_coord(cdl_ijk,'lon',im,units='degrees_east',
     &     coordvalues=lon_dg(:,1))
      call add_coord(cdl_ijk,'lon2',im,units='degrees_east',
     &     coordvalues=lon_dg(:,2))
      call add_coord(cdl_ijk,'lat',jm,units='degrees_north',
     &     coordvalues=lat_dg(:,1))
      call add_coord(cdl_ijk,'lat2',jm,units='degrees_north',
     &     coordvalues=lat_dg(:,2))
      cdl_dum = cdl_ijk
      call merge_cdl(cdl_dum,cdl_heights,cdl_ijk)

      do k=1,kaijk
        if(trim(units_ijk(k)).eq.'unused') cycle
        ijkg = jgrid_ijk(k) ! contains igrid,jgrid,kgrid in successive bits
        kgrid = ijkg/kgride
        ijkg = ijkg - kgride*kgrid
        jgrid = ijkg/jgride
        ijkg = ijkg - jgride*jgrid
        igrid = ijkg
        if(igrid.eq.0) then
          xstr='lon) ;'
        else
          xstr='lon2) ;'
        endif
        if(jgrid.eq.0) then
          ystr='lat,'
        else
          ystr='lat2,'
        endif
        if(kgrid.eq.0) then
          zstr='(plm,'
        else
          zstr='(ple,'
        endif
        set_miss = denom_ijk(k).ne.0
        call add_var(cdl_ijk,
     &       'float '//trim(name_ijk(k))//
     &       trim(zstr)//trim(ystr)//trim(xstr),
     &       units=trim(units_ijk(k)),
     &       long_name=trim(lname_ijk(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
      enddo

      return
      end subroutine ijk_defs


      SUBROUTINE DIAGB
!@sum DIAGB calculate constant pressure diagnostics from within DYNAM
C****
C**** CONTENTS OF AGC(J,K,N)  (SUM OVER LONGITUDE AND TIME OF)
C****   See gc_defs for contents
C****
C**** CONTENTS OF AIJK(I,J,K,N)   (SUM OVER TIME OF)
C****   See ijks_defs for contents
C****
      Use CONSTANT,   Only: kg2mb,lhe,omega,sha,tf,teeny, radius
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only : idacc,mdyn,mdiag
      Use ATM_COM,    Only: PHI,MA,PDSIG,PMID,PEDN,U,V,T,Q,QCL,QCI
      USE GEOM, only : bydxyp,bydxyv,rapvs,rapvn,
     &     COSV,DXV,DXYN,DXYP,DXYS,DXYV,DYP,DYV,FCOR,IMAXJ
      USE DIAG_COM, only : imh,fim,byim,ia_dga,ndaa
     &    ,aijk=>aijk_loc, ! adiurn,hdiurn
     &     ndiupt,hr_in_day
     *     ,idd_w,ijdd
     &     ,aij=>aij_loc,ij_puq,ij_pvq,ij_dsev
      USE GCDIAG
      USE GC_COM, only : jeq,
     &    agc=>agc_loc,speca,nspher,klayer,nwav_dag
      Use DYNAMICS,  Only: CONV,dut,dvt,SD
      USE DIAG_LOC, only : w,tx,pm,pl,pmo,plo
     &     ,ldna,lupa
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : getDomainBounds,CHECKSUM,HALO_UPDATE
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATEj, HALO_UPDATE_COLUMN
      USE DOMAIN_DECOMP_1D, only : SOUTH, NORTH, GLOBALSUM
      USE DOMAIN_DECOMP_1D, only : SUMXPE, broadcast, AM_I_ROOT
      USE PRECISION_MOD
      USE GETTIME_MOD
      IMPLICIT NONE

      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE,KE_jsum
      REAL*8, DIMENSION
     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part
      REAL*8, DIMENSION(IM,GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &     ZX,STB,UDX
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &     STJK,DPJK,UJK,VJK,WJK,TJK,
     &     PSIJK,UP,TY,PSIP,WTJK,UVJK,WUJK
      REAL*8, DIMENSION(IM) :: X1,X1tmp,WPA2_of_lon
      REAL*8, DIMENSION(LM) :: SHETH,DPM,DTH,PDSIGL
      REAL*8, DIMENSION(LM+1) :: PEDNL

      INTEGER ::
     &     I,IH,IHM,IM1,INCH,INCHM,IP1,IZERO,J,J45N,
     &     JHEMI,K,KDN,KR,KS,KS1,KSPHER,KUP,KX,L,
     &     LUP,N,NM,KM

      REAL*8 ::
     &     begin, BYDP,BYFIM,DP,DPDN,DP4,
     &     DPE,DPI,DPK,DPSQI,DPUP,DPUV,DUTI,DUTK,DVTI,DVTK,FIMI,
     &     MAUV,
     &     PAI,PAK,PDN,PMK,PQ4I,PQ4K,PQV4I,PS4I,
     &     PS4K,PSIY,PSV4I,PT4I,PT4K,PTK,PTV4I,PUI,PUK,PUP,
     &     PUVI,PV2,PV2I,PVI,PVK,PWWI,PWWVI,PY,PZ4I,PZ4K,
     &     PZV4I,QK,QKI, SDL,SDLm1,
     &     SMALL,SQRTDP,THK,THKI,THPI,TK,TKI,TPI,
     &     UDUTI,    UEARTH,UK,UY,VDVTI,VK,VSTAR,W2,W2I,W4,
     &     WI,WNP,WPA2I,WPV4I,WQI,WSP,WSTAR,WTHI,
     &     WTI,WUP,WZI,ZK,ZKI,
     &     FIMIofK(LM),UKIofK(LM),W4IofK(LM),WKE4IofK(LM),WU4IofK(LM)
     &     ,AMRHT,AMRHQ,AMUV,AMVQ,AMVT,AMUU,AMVV,AMTT

c local vars for transplanted DIAGA calculations
      Real*8 :: dudp,dthdp,umn,thmn,fphi,dudx,pvthp,CsLMN,CsL(IM),
     &          CsLPU,upe,vpe, p4i,p4,pu4i,pv4i,puv4i,t4,z4,sp2
      integer :: ldn
      real*8, dimension(im) :: thsec

      REAL*8, PARAMETER :: BIG=1.E20
      REAL*8 :: QSAT
      REAL*8 :: pm_ge_ps(im,grid%j_strt_halo:grid%j_stop_halo,lm)
      INTEGER :: J_0, J_1, J_0S, J_1S, J_0STG, J_1STG, J_0H
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

      CALL GETTIME(BEGIN)

      call getDomainBounds(grid, J_STRT=J_0,         J_STOP=J_1,
     &               J_STRT_SKP=J_0S,    J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               J_STRT_HALO=J_0H,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)

      call halo_update(grid, tx)

C****
C**** EASTWARD TRANSPORTS (TRANSPLANTED FROM DIAGA)
C****
      CALL HALO_UPDATE(grid, U, FROM=NORTH)
      DO L=1,LM
      DO J=J_0S,J_1S
      I=IM
      DO IP1=1,IM
        AIJ(I,J,IJ_PUQ) = AIJ(I,J,IJ_PUQ)+(PDSIG(L,I,J)+PDSIG(L,IP1,J))*
     *       (U(I,J,L)+U(I,J+1,L))*(Q(I,J,L)+Q(IP1,J,L))*.125
        I=IP1
      END DO
      END DO
      END DO
C****
C**** NORTHWARD TRANSPORTS (TRANSPLANTED FROM DIAGA)
C****
      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
      CALL HALO_UPDATE(grid, Q, FROM=SOUTH)

      DO J=J_0STG,J_1STG
         Do L=1,LM
            AGC(J,L,JL_DPB) = AGC(J,L,JL_DPB) +
     +         .5 * Sum(PDSIG(L,:,J-1) + PDSIG(L,:,J))  ;  EndDo
      DO L=1,LM
        PU4I=0.
        PV4I=0.
        PUV4I = 0  ;  P4I = 0 
        I=IM
        DO IP1=1,IM
          P4=PDSIG(L,I,J-1)+PDSIG(L,IP1,J-1)+PDSIG(L,I,J)+PDSIG(L,IP1,J)
          P4I = P4I + P4
          PU4I=PU4I+P4*U(I,J,L)
          PV4I=PV4I+P4*V(I,J,L)
          PUV4I=PUV4I+P4*U(I,J,L)*V(I,J,L)
          T4=TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L)
          Z4=PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L)
          AIJ(I,J,IJ_DSEV)=AIJ(I,J,IJ_DSEV)+P4*(SHA*T4+Z4)*V(I,J,L)
     *       *DXV(J)*0.0625d0
          SP2 = PDSIG(L,IP1,J-1) + PDSIG(L,IP1,J)
          AIJ(IP1,J,IJ_PVQ)=AIJ(IP1,J,IJ_PVQ)+.125*SP2
     *       *(V(I,J,L)+V(IP1,J,L))*(Q(IP1,J-1,L)+Q(IP1,J,L))
          I=IP1
        END DO
        AGC(J,L,JL_ZMFNTMOM)=AGC(J,L,JL_ZMFNTMOM)+.25*PU4I*PV4I/P4I
        AGC(J,L,JL_TOTNTMOM)=AGC(J,L,JL_TOTNTMOM)+.25*PUV4I
      END DO
      END DO

C****
C**** ELIASSEN PALM FLUX ON MODEL LAYERS (TRANSPLANTED FROM DIAGA)
C****
C**** NORTHWARD COMPONENT
      CALL HALO_UPDATE(grid, T, FROM=SOUTH)

      DO 868 J=J_0STG,J_1STG
      DO 868 L=1,LM
      DUDP=0.
      DTHDP=0.
      UMN=0.
      THMN=0.
      LDN=LDNA(L)
      LUP=LUPA(L)
      I=IM
      DO 864 IP1=1,IM
      DUDP=DUDP+U(I,J,LUP)-U(I,J,LDN)
      DTHDP=DTHDP+T(I,J,LUP)+T(I,J-1,LUP)-T(I,J,LDN)-T(I,J-1,LDN)
      UMN=UMN+U(I,J,L)
      THMN=THMN+T(I,J,L)+T(I,J-1,L)
      THSEC(I)=T(I,J,L)+T(IP1,J,L)+T(I,J-1,L)+T(IP1,J-1,L)
  864 I=IP1
      UMN=UMN*BYIM
      THMN=2.*THMN/FIM
      FPHI=0.
      SMALL=.0002d0*FIM*T(1,J,L)
c      IF (DTHDP.LT.SMALL) WRITE (6,999) J,L,DTHDP,SMALL
      IF (DTHDP.LT.SMALL) DTHDP=SMALL
      I=IM
      DO 866 Ip1=1,IM
      MAUV = (MA(L,I,J  ) + MA(L,Ip1,J  )) * RAPVS(J)+
     *       (MA(L,I,J-1) + MA(L,Ip1,J-1)) * RAPVN(J-1)
      FPHI = FPHI + MAUV*V(I,J,L)*(.5*(THSEC(I)-THMN)*DUDP/DTHDP
     *   -U(I,J,L)+UMN)
  866 I=Ip1
  868 AGC(J,L,JL_EPFLXN)=AGC(J,L,JL_EPFLXN)+FPHI

C**** VERTICAL COMPONENT
      CALL HALO_UPDATE(grid, V, FROM=NORTH)

      DO 878 J=J_0S,J_1S
      CsLMN = 0  ;  CsL(:) = 0
      DO 878 L=1,LM-1
      CsLMN = Sum(CONV(:,J,L)) / IM + CsLMN
      DTHDP = Sum(T(:,J,L+1)-T(:,J,L))
      THMN  = Sum(T(:,J,L+1)+T(:,J,L)) / IM
      DUDX  = DXV(J+1)*Sum(U(:,J+1,L)+U(:,J+1,L+1))
     -      - DXV(J)  *Sum(U(:,J  ,L)+U(:,J  ,L+1))
      SMALL=.0001d0*FIM*T(1,J,L+1)
      IF (DTHDP.LT.SMALL) DTHDP=SMALL
      PVTHP = 0  ;  CsLPU = 0
      IM1=IM
      DO 874 I=1,IM
      CsL(I) = CsL(I) + CONV(I,J,L)
      UPE = U(IM1,J,L  )+U(IM1,J+1,L  )+U(I,J,L  )+U(I,J+1,L  )+
     *      U(IM1,J,L+1)+U(IM1,J+1,L+1)+U(I,J,L+1)+U(I,J+1,L+1)
      VPE = V(IM1,J,L  )+V(IM1,J+1,L  )+V(I,J,L  )+V(I,J+1,L  )+
     *      V(IM1,J,L+1)+V(IM1,J+1,L+1)+V(I,J,L+1)+V(I,J+1,L+1)
      PVTHP = PVTHP + VPE*MA(L,I,J)*(T(I,J,L)+T(I,J,L+1)-THMN)
      CsLPU = CsLPU - UPE*(CsL(I)-CsLMN)
  874 IM1=I
      AGC(J,L,JL_EPFLXV) = AGC(J,L,JL_EPFLXV) + .25*byDXYP(J)*
     &     ((.5*FIM*FCOR(J)-.25*DUDX)*PVTHP/DTHDP + CsLPU)*kg2mb
  878 CONTINUE

c
c END OF CALCULATIONS MOVED FROM DIAGA
c

      pm_ge_ps(:,:,:)=-1.
C****
C**** INTERNAL QUANTITIES T,TH
C****
      KM=LM
      DO 170 J=J_0,J_1
      DO 170 K=1,KM
      DPI=0.
      TPI=0.
      THPI=0.
      FIMI=0.
      DO 160 I=1,IMAXJ(J)
C**** FIND L=L(K) AND LUP=L(K+1) S.T. P(LUP).GT.P(K+1)
      If (PM(K+1) >= PEDN(1,I,J))  GoTo 160
      L=1
      PDN = PEDN(1,I,J)
      If (PM(K) >= PEDN(1,I,J))  GoTo 120
      PDN=PM(K)
  110 If (PM(K) > PEDN(L+1,I,J))  GoTo 120
      L=L+1
      GO TO 110
  120 LUP=L
  130 If (PM(K+1) >= PEDN(LUP+1,I,J) .or. LUP == LM)  GoTo 140
      LUP=LUP+1
      GO TO 130
  140 CONTINUE
C**** ACCUMULATE HERE
      DPI=DPI+PDN-PM(K+1)
      FIMI=FIMI+1.
  150 PUP=PEDN(L+1,I,J)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      TPI=TPI+(TX(I,J,L)-TF)*DP
      THPI=THPI+T(I,J,L)*DP
      IF (L.EQ.LUP) GO TO 160
      L=L+1
      PDN = PEDN(L,I,J)
      GO TO 150
  160 CONTINUE
      AGC(J,K,JK_NPTSAVG1)=AGC(J,K,JK_NPTSAVG1)+FIMI
      AGC(J,K,JK_DPA)=AGC(J,K,JK_DPA)+DPI
      AGC(J,K,JK_TEMP)=AGC(J,K,JK_TEMP)+TPI
      AGC(J,K,JK_THETA)=AGC(J,K,JK_THETA)+THPI
         TJK(J,K)=THPI/(DPI+teeny)
         IF (IDACC(ia_dga).EQ.1) AGC(J,K,JK_TINST)=TJK(J,K)
         AGC(J,K,JK_TOTDTDT)=TJK(J,K)-AGC(J,K,JK_TINST)
  170 CONTINUE
C****
C**** CALCULATE STABILITY AT ODD LEVELS ON PU GRID
C****
      DO 230 J=J_0,J_1
      I=IMAXJ(J)
      DO 230 IP1=1,IMAXJ(J)
      PLO(:) = .5*(PMID(:,I,J) + PMID(:,Ip1,J))
       PL(:) = .5*(PEDN(:,I,J) + PEDN(:,Ip1,J))
      DO 180 L=1,LM-1
      DTH(L)=(T(I,J,L)+T(IP1,J,L)-T(I,J,L+1)-T(IP1,J,L+1))/
     *  (2.*(PLO(L)-PLO(L+1)))
  180 CONTINUE
      DO 220 K=1,KM
      STB(I,J,K)=0.
      IF (PM(K+1).GE.PL(1)) GO TO 220
      PMK=PMO(K)
      If (PM(K) > PL(1))  PMK = .5*(PL(1) + PM(K+1))
      L=2
      IF (PMK.GE.PL(2)) GO TO 210
  190 LUP=L+1
      IF (L.EQ.LM) GO TO 210
      IF (PMK.GE.PL(LUP)) GO TO 200
      L=LUP
      GO TO 190
  200 DPUP=PMK-PL(LUP)
      DPDN=PL(L)-PMK
      STB(I,J,K)=(DTH(L-1)*DPUP+DTH(L)*DPDN)/(DPUP+DPDN+teeny)
      GO TO 220
C**** SPECIAL CASES,  L=2, L=LM
  210 STB(I,J,K)=DTH(L-1)
  220 CONTINUE
  230 I=IP1
C**** CALCULATE STJK; THE MEAN STATIC STABILITY
      DO 260 J=J_0,J_1
      DO 260 K=1,KM
      STJK(J,K)=0.
      DPJK(J,K)=0.
      I=IMAXJ(J)
      DO 250 IP1=1,IMAXJ(J)
      If (PM(K+1) > .5*(PEDN(1,I,J) + PEDN(1,Ip1,J)))  GoTo 250
      STJK(J,K)=STJK(J,K)+STB(I,J,K)
      DPJK(J,K)=DPJK(J,K)+1.
  250 I=IP1
      STJK(J,K)=STJK(J,K)/(DPJK(J,K)+teeny)
      SMALL=.0001
      IF (ABS(STJK(J,K)).LT.SMALL) STJK(J,K)=-SMALL
  260 CONTINUE
C****
C**** CONSTANT PRESSURE DIAGNOSTICS:  FLUX, ENERGY, ANGULAR MOMENTUM
C****
      ZX(:,:,:)=0.
      IF(HAVE_SOUTH_POLE) THEN
        UDX(:,1,:)=0.
      ENDIF

C Needs to check later if all these halo calls are necessary.
C DIAGA may have contained relevant halo calls
C and since DIAGB is called immediately after DIAGA
C there may not be a need for these calls if
C the concerned arrays have not been updated
C from the previous halo call.
c      CALL HALO_UPDATE(grid, TX, FROM=SOUTH)
c      CALL HALO_UPDATE(grid, PHI, FROM=SOUTH)
c      CALL HALO_UPDATE(grid, Q, FROM=SOUTH)
c      CALL HALO_UPDATE(grid, T, FROM=SOUTH)
      CALL HALO_UPDATEj(grid, STJK, FROM=SOUTH)
c***      DO L=1,LM
c***         CALL HALO_UPDATE(grid, STJK(:,L), FROM=SOUTH)
c***      END DO

      DO 390 J=J_0STG,J_1STG
      I=IM
      DO 280 IP1=1,IM
      UDX(I,J,1:KM) = 0
      DO L=1,LM
        PDSIGL(L) = (PDSIG(L,I,J  ) + PDSIG(L,IP1,J  ))*RAPVS(J) +
     +              (PDSIG(L,I,J-1) + PDSIG(L,IP1,J-1))*RAPVN(J-1)
        DUT(I,J,L)=DUT(I,J,L)/(PDSIGL(L)*DXYV(J))
        DVT(I,J,L)=DVT(I,J,L)/(PDSIGL(L)*DXYV(J))
      END DO
  280 I=IP1
      DO 350 K=1,KM
      DPI=0.
      DPSQI=0.
      FIMI=0.
      PUI=0.
      PVI=0.
      PWWI=0.
      PT4I=0.
      PTV4I=0.
      PZ4I=0.
      PZV4I=0.
      PQ4I=0.
      PQV4I=0.
      PWWVI=0.
      PUVI=0.
      DVTI=0.
      VDVTI=0.
      DUTI=0.
      UDUTI=0.
      PS4I=0.
      PSV4I=0.
      I=IM
      DO 340 IP1=1,IM
      DO 286 L=1,LM+1
  286 PL(L) = (PEDN(L,I,J  ) + PEDN(L,Ip1,J  )) * RAPVS(J)+
     +        (PEDN(L,I,J-1) + PEDN(L,Ip1,J-1)) * RAPVN(J-1)
      If (PM(K+1) >= PL(1)) Then
        pm_ge_ps(i,j,k) = 1.
        UDX(I,J,K)=BIG
      ELSE
        L=1
        PDN = PL(1)
        If (PM(K) >= PL(1))  GoTo 300
        PDN=PM(K)
  290   IF (PM(K).GT.PL(L+1)) GO TO 300
        L=L+1
        If (L <= LM)  GoTo 290
  300   LUP=L
  310   IF (PM(K+1).GE.PL(LUP+1)) GO TO 320
        LUP=LUP+1
        If (LUP <= LM)  GoTo 310
  320   CONTINUE
        DPK=PDN-PM(K+1)
        PUK=0.
        PVK=0.
        PT4K=0.
        PZ4K=0.
        PQ4K=0.
        DUTK=0.
        DVTK=0.
        PS4K=0.
C**** FOR AMIP
        AMVQ=0.
        AMVT=0.
        AMUU=0.
        AMVV=0.
        AMUV=0.
        AMTT=0.
C**** END AMIP
C**** INTERPOLATE HERE
  330 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      DP4=.25*DP
      PUK=PUK+DP*U(I,J,L)
      PVK=PVK+DP*V(I,J,L)
      PT4K=PT4K+(TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L))*DP4
      PZ4K=PZ4K+(PHI(I,J-1,L)+PHI(IP1,J-1,L)+PHI(I,J,L)+PHI(IP1,J,L))
     &     *DP4
      PQ4K=PQ4K+(Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L))*DP4
      DUTK=DUTK+DP*DUT(I,J,L)
      DVTK=DVTK+DP*DVT(I,J,L)
      PS4K=PS4K+(T(I,J-1,L)+T(IP1,J-1,L)+T(I,J,L)+T(IP1,J,L))*DP4
C**** FOR AMIP 2
      AMRHT=.25*(TX(I,J-1,L)+TX(IP1,J-1,L)+TX(I,J,L)+TX(IP1,J,L))
      AMRHQ=.25*(Q(I,J-1,L)+Q(IP1,J-1,L)+Q(I,J,L)+Q(IP1,J,L))
      AMVQ=AMVQ+DP*V(I,J,L)*AMRHQ
      AMVT=AMVT+DP*V(I,J,L)*AMRHT
      AMUU=AMUU+DP*U(I,J,L)*U(I,J,L)
      AMVV=AMVV+DP*V(I,J,L)*V(I,J,L)
      AMUV=AMUV+DP*U(I,J,L)*V(I,J,L)
      AMTT=AMTT+DP*AMRHT*AMRHT
C**** END AMIP
      IF (LUP.EQ.L) GO TO 332
      L=L+1
      PDN=PL(L)
      If (L <= LM)  GoTo 330
C**** ACCUMULATE HERE
  332 FIMI=FIMI+1.
      DPI=DPI+DPK
      DPSQI=DPSQI+DPK*DPK
      IF (DPK.LT.teeny) DPK=teeny
      BYDP=1./DPK
      PUI=PUI+PUK
      PVI=PVI+PVK
      PWWI=PWWI+BYDP*(PUK*PUK+PVK*PVK)
      PWWVI=PWWVI+BYDP*BYDP*(PUK*PUK+PVK*PVK)*PVK
      PUVI=PUVI+BYDP*PUK*PVK
      PT4I=PT4I+PT4K
      PTV4I=PTV4I+BYDP*PT4K*PVK
      PZ4I=PZ4I+PZ4K
      PZV4I=PZV4I+BYDP*PZ4K*PVK
      PQ4I=PQ4I+PQ4K
      PQV4I=PQV4I+BYDP*PQ4K*PVK
      DVTI=DVTI+DVTK
      VDVTI=VDVTI+BYDP*PVK*DVTK
      DUTI=DUTI+DUTK
      UDUTI=UDUTI+BYDP*PUK*DUTK
!!    IF(SKIPSE.EQ.1.) GO TO 334
      AIJK(I,J,K,IJK_UB) =AIJK(I,J,K,IJK_UB) +PUK
      AIJK(I,J,K,IJK_VB) =AIJK(I,J,K,IJK_VB) +PVK
      AIJK(I,J,K,IJK_DSE)=AIJK(I,J,K,IJK_DSE)+SHA*PT4K+PZ4K
      AIJK(I,J,K,IJK_DPB)=AIJK(I,J,K,IJK_DPB)+DPK
      AIJK(I,J,K,IJK_TB) =AIJK(I,J,K,IJK_TB) +PT4K
      AIJK(I,J,K,IJK_PF)  =AIJK(I,J,K,IJK_PF)+1.
C     *  *  *  FOR AMIP 2  *  *  *
      AIJK(I,J,K,IJK_UV)=AIJK(I,J,K,IJK_UV)+AMUV
      AIJK(I,J,K,IJK_VQ)=AIJK(I,J,K,IJK_VQ)+AMVQ
      AIJK(I,J,K,IJK_VT)=AIJK(I,J,K,IJK_VT)+AMVT
      AIJK(I,J,K,IJK_UU)=AIJK(I,J,K,IJK_UU)+AMUU
      AIJK(I,J,K,IJK_VV)=AIJK(I,J,K,IJK_VV)+AMVV
      AIJK(I,J,K,IJK_TT)=AIJK(I,J,K,IJK_TT)+AMTT
C**** END AMIP
C**** EDDY TRANSPORT OF THETA;  VORTICITY
  334   PS4I=PS4I+PS4K
        PSV4I=PSV4I+BYDP*PVK*PS4K
        UDX(I,J,K)=BYDP*PUK*DXV(J)
!ESMF   IF (UDX(I,J-1,K).LT.BIG) ZX(I,J-1,K)=UDX(I,J,K)-UDX(I,J-1,K)
!ESMF   IF (UDX(I,J-1,K).GE.BIG) ZX(I,J-1,K)=0.
!ESMF   IF (ZX(I,J-1,K).GE.BIG) ZX(I,J-1,K)=0.
      END IF                            !---> (PM(K+1).GE.PS)

  340 I=IP1     !-->END I Loop (IP1) from IM to IM-1 (1-IM).
      DPM(K)=DPI/(FIMI+teeny)
      DPJK(J,K)=DPI
      AGC(J,K,JK_DPB)=AGC(J,K,JK_DPB)+DPI
      AGC(J,K,JK_DPSQR)=AGC(J,K,JK_DPSQR)+DPSQI
      AGC(J,K,JK_NPTSAVG)=AGC(J,K,JK_NPTSAVG)+FIMI
      IF (DPI.LT.teeny) DPI=teeny
      AGC(J,K,JK_U)=AGC(J,K,JK_U)+PUI
      AGC(J,K,JK_V)=AGC(J,K,JK_V)+PVI
      AGC(J,K,JK_EDDKE)=AGC(J,K,JK_EDDKE)+PWWI-(PUI*PUI+PVI*PVI)/DPI
      AGC(J,K,JK_TOTKE)=AGC(J,K,JK_TOTKE)+PWWI
      AGC(J,K,JK_EDDNTSH)=AGC(J,K,JK_EDDNTSH)+PTV4I-PT4I*PVI/DPI
      AGC(J,K,JK_TOTNTDSE)=AGC(J,K,JK_TOTNTDSE)+SHA*PTV4I+PZV4I
      AGC(J,K,JK_EDDNTGEO)=AGC(J,K,JK_EDDNTGEO)+PZV4I-PZ4I*PVI/DPI
      AGC(J,K,JK_EDDNTLH)=AGC(J,K,JK_EDDNTLH)+PQV4I-PQ4I*PVI/DPI
      AGC(J,K,JK_TOTNTLH)=AGC(J,K,JK_TOTNTLH)+PQV4I
      AGC(J,K,JK_TOTNTKE)=AGC(J,K,JK_TOTNTKE)+PWWVI
      AGC(J,K,JK_EDDNTMOM)=AGC(J,K,JK_EDDNTMOM)+PUVI-PUI*PVI/DPI
      AGC(J,K,JK_TOTNTMOM)=AGC(J,K,JK_TOTNTMOM)+PUVI
      AGC(J,K,JK_P2KEDPGF)=AGC(J,K,JK_P2KEDPGF)+VDVTI+UDUTI-
     *   (PUI*DUTI+PVI*DVTI)/DPI
      SHETH(K)=(PSV4I-PS4I*PVI/DPI)*DXYV(J)/(STJK(J-1,K)*DXYN(J-1)+
     *   STJK(J,K)*DXYS(J))
         UJK(J,K)=PUI/DPI
         VJK(J,K)=PVI/DPI
         PSIJK(J,K)=SHETH(K)/DPI
         UVJK(J,K)=(PUVI-PUI*PVI/DPI)/DPI
         IF (IDACC(ia_dga).EQ.1) AGC(J,K,JK_UINST)=UJK(J,K)
         AGC(J,K,JK_TOTDUDT)=UJK(J,K)-AGC(J,K,JK_UINST)
  350 AGC(J,K,JK_SHETH)=AGC(J,K,JK_SHETH)+SHETH(K)
  390 CONTINUE

C**** ZX for distributed parallelization
c****
      CALL HALO_UPDATE( grid, UDX, from=NORTH )
      CALL HALO_UPDATE( grid, pm_ge_ps, from=NORTH)

      DO J=J_0,J_1S
        DO K=1,KM
          DO I=1,IM
            if (pm_ge_ps(i,j+1,k) < 0) then
            IF (UDX(I,J,K).LT.BIG ) ZX(I,J,K)=-UDX(I,J,K)+UDX(I,J+1,K)
            IF (UDX(I,J,K).GE.BIG)  ZX(I,J,K)=0.
            IF (ZX(I,J,K).GE.BIG)   ZX(I,J,K)=0
            end if
          END DO
        END DO
      END DO
C****
C**** alternate vertical mass flux diagnostic (from SD)
C****
!**** Linearly interpolate SD to constant pressure coordinates between
!**** model layer edges, assumes SD = 0 at top and bottom of column.
!**** For LM-1 model edges: PM(2:LM),PEDN(2:LM),SD(1:LM-1),W(1:LM-1)  
      Do J=J_0,J_1  ;  Do I=1,IM
      W(I,J,KM) = 0
      K = KM-1  ;  L = LM  ;  SDL = 0  ;  SDLm1 = SD(I,J,LM-1)
  410 If (PM(K+1) <= PEDN(L,I,J))
     *   Then  ;  W(I,J,K) = ((PEDN(L,I,J)-PM(K+1))*SDL +
     +                        (PM(K+1)-PEDN(L+1,I,J))*SDLm1) /
     /                       ((PEDN(L,I,J)-PEDN(L+1,I,J)) * DXYP(J))
                  AIJK(I,J,K,IJK_W) = AIJK(I,J,K,IJK_W) + W(I,J,K)
                  If (K==1) Cycle
                  K = K-1  ;  GoTo 410
         Else  ;  L = L-1  ;  SDL = SDLm1
                  If (L > 1)
     *               Then  ;  SDLm1 = SD(I,J,L-1)  ;  GoTo 410  ;  EndIf
                  If (L == 1)
     *               Then  ;  SDLm1 = 0  ;  GoTo 410
                     Else  ;  W(I,J,1:K) = 0  ;  Cycle  ;  EndIf
         EndIf  ;  EndDo  ;  EndDo

C**** ACCUMULATE ALL VERTICAL WINDS
!!    DO 558 J=J_0,J_1
!!    DO 558 I=1,IM
!!    DO KR=1,NDIUPT
!!       IF(I.EQ.IJDD(1,KR).AND.J.EQ.IJDD(2,KR)) THEN
!!*** Warning:     This diagnostic has 3 flaws   (?)
!!***          1 - It assumes that DTsrc=1hr, (DTsrc=3600.)
!!***          2 - since DTdaa-Ndaa*DTsrc=2*DTdyn rather than 0,
!!***              some hours are skipped once in a while
!!***          3 - Some of the first Ndaa hours are skipped at the
!!***              beginning of a month and overcounted at the end;
!!***              this happens to balance out, if and only if
!!***              mod(days_in_month,ndaa)=0  (i.e. February if Ndaa=7)
!!***          In addition, IHM occasionally is out-of-bounds.
!!          IH=JHOUR+1
!!          IHM = IH+(JDATE-1)*24
!!          DO INCH=1,NDAA
!!            IF(IH.GT.HR_IN_DAY) IH=IH-HR_IN_DAY
!!            ADIURN(IDD_W,KR,IH)=ADIURN(IDD_W,KR,IH)+1.E5*W(I,J,3)
!!            HDIURN(IDD_W,KR,IHM)=HDIURN(IDD_W,KR,IHM)+1.E5*W(I,J,3)
!!            IH=IH+1
!!            IHM=IHM+1
!!          END DO
!!       END IF
!!    END DO
!!558 CONTINUE

      DO 565 J=J_0,J_1
      DO 565 K=1,KM
      WI=0.
      DO I=1,IMAXJ(J)
        WI=WI+W(I,J,K)
      END DO
  565 AGC(J,K,JK_VVEL)=AGC(J,K,JK_VVEL)+WI
C****
C**** ACCUMULATE T,Z,Q VERTICAL TRANSPORTS
C****
      DO 610 J=J_0,J_1
      DO 610 K=2,KM
      WI=0.
      TKI=0.
      QKI=0.
      ZKI=0.
      WTI=0.
      WQI=0.
      WZI=0.
         THKI=0.
         WTHI=0.
      FIMI=0.
      DO 600 I=1,IMAXJ(J)
      DO 569 L=1,LM
  569 PLO(L) = PMID(L,I,J)
      If (PM(K) >= PEDN(1,I,J))  GoTo 600
      L=1
      IF (PM(K).GE.PLO(1)) GO TO 580
  570 LUP=L+1
      IF (L.EQ.LM) GO TO 580
      IF (PM(K).GE.PLO(LUP)) GO TO 575
      L=LUP
      GO TO 570
  575 DPUP=PM(K)-PLO(LUP)
      DPDN=PLO(L)-PM(K)
      BYDP=1./(DPDN+DPUP)
      TK=BYDP*(TX(I,J,L)*DPUP+TX(I,J,LUP)*DPDN)
      QK=Q(I,J,L)*Q(I,J,LUP)/(BYDP*(Q(I,J,L)*DPDN+Q(I,J,LUP)*DPUP)+
     *  teeny)
      ZK=BYDP*(PHI(I,J,L)*DPUP+PHI(I,J,LUP)*DPDN)
         THK=BYDP*(T(I,J,L)*DPUP+T(I,J,LUP)*DPDN)
      GO TO 590
C**** SPECIAL CASES;  L=1, L=LM
  580 TK=TX(I,J,L)
      QK=Q(I,J,L)
      ZK=PHI(I,J,L)
         THK=T(I,J,L)
C**** MERIDIONAL AVERAGING
  590 WI=WI+W(I,J,K)
      TKI=TKI+TK
      QKI=QKI+QK
      ZKI=ZKI+ZK
      WTI=WTI+W(I,J,K)*TK
      WQI=WQI+W(I,J,K)*QK
      WZI=WZI+W(I,J,K)*ZK
         THKI=THKI+THK
         WTHI=WTHI+W(I,J,K)*THK
      FIMI=FIMI+1.
  600 CONTINUE
      BYFIM=teeny
      IF (FIMI.GT.teeny) BYFIM=1./FIMI
      AGC(J,K-1,JK_TOTVTDSE)=AGC(J,K-1,JK_TOTVTDSE)+SHA*WTI+WZI
      AGC(J,K-1,JK_TOTVTLH)=AGC(J,K-1,JK_TOTVTLH)+WQI
C     AGC(J,K-1,JK_BAREKEGEN)=AGC(J,K-1,JK_BAREKEGEN)+WTI-BYFIM*WI*TKI
      WJK(J,K)=BYFIM*WI
      WTJK(J,K)=BYFIM*(WTHI-BYFIM*WI*THKI)
      if(j>1 .and. j<jm) then
        AGC(J,K-1,JK_EDDVTPT)=AGC(J,K-1,JK_EDDVTPT)+WTJK(J,K)
        AGC(J,K-1,JK_EDDVTDSE)=AGC(J,K-1,JK_EDDVTDSE)+SHA*WTI+WZI
     &       -BYFIM*(SHA*TKI+ZKI)*WI
        AGC(J,K-1,JK_EDDVTLH)=AGC(J,K-1,JK_EDDVTLH)+WQI-BYFIM*QKI*WI
        AGC(J,K-1,JK_VTGEOEDDY)=AGC(J,K-1,JK_VTGEOEDDY)+WZI-BYFIM*WI*ZKI
      endif
  610 CONTINUE
C****
C**** BAROCLINIC EDDY KINETIC ENERGY GENERATION
C****
      DO J=J_0S,J_1S
      DO K=1,KM
      FIMI=0.
      W2I=0.
      PAI=0.
      WPA2I=0.
      DO 626 I=1,IMAXJ(J)
      DO 611 L=1,LM
  611 PL(L) = PEDN(L,I,J)
      If (PM(K+1) >= PEDN(1,I,J))  GoTo 626
      L=1
      PDN = PEDN(1,I,J)
      If (PM(K) >= PEDN(1,I,J))  GoTo 614
      PDN=PM(K)
  612 IF (PM(K).GT.PL(L+1)) GO TO 614
      L=L+1
      If (L <= LM)  GoTo 612
  614 LUP=L
  616 IF (PM(K+1).GE.PL(LUP+1)) GO TO 618
      LUP=LUP+1
      If (LUP <= LM)  GoTo 616
  618 CONTINUE
      PTK=0.
C**** INTERPOLATE HERE
  620 PUP=PL(L+1)
      IF (LUP.EQ.L) PUP=PM(K+1)
      DP=PDN-PUP
      PTK=PTK+DP*TX(I,J,L)
      IF (LUP.EQ.L) GO TO 622
      L=L+1
      PDN=PL(L)
      If (L <= LM)  GoTo 620
C**** ACCUMULATE HERE
  622 FIMI=FIMI+1.
      WUP=0.
      IF (K.LT.KM) WUP=W(I,J,K+1)
      W2I=W2I+W(I,J,K)+WUP
      PY=PMO(K)
      If (PM(K) >= PL(1))  PY = .5*(PL(1) + PM(K+1))
      PAK=PTK/PY
      PAI=PAI+PAK
      WPA2I=WPA2I+(W(I,J,K)+WUP)*PAK
      WPA2_of_lon(I)=(W(I,J,K)+WUP)*PAK
  626 CONTINUE
      AGC(J,K,JK_BAREKEGEN)=AGC(J,K,JK_BAREKEGEN)-
     &     (WPA2I-W2I*PAI/(FIMI+teeny))
      DO I=1,IMAXJ(J)
        If (PM(K+1) < PEDN(1,I,J))  Then
          AIJK(I,J,K,IJK_BAREKEGEN)=AIJK(I,J,K,IJK_BAREKEGEN)-
     &         (WPA2_of_lon(I)-W2I*PAI/(FIMI*FIMI+teeny))
        ENDIF
      ENDDO
      ENDDO ! K
      ENDDO ! J
C****
C**** ACCUMULATE UV VERTICAL TRANSPORTS
C****
C**** RESCALE POLAR WINDS
      DO 640 K=1,KM
        IF(HAVE_SOUTH_POLE) THEN
          WSP=W(1,1,K)/FIM
          DO I=1,IM
            W(I,1,K)=WSP
          ENDDO
        ENDIF
        IF(HAVE_NORTH_POLE) THEN
          WNP=W(1,JM,K)/FIM
          DO I=1,IM
            W(I,JM,K)=WNP
          ENDDO
        ENDIF
  640 CONTINUE

C P already halo'ed; no need      CALL CHECKSUM(grid, P, __LINE__, __FILE__)
C P already halo'ed; no need     CALL HALO_UPDATE(grid, P, FROM=SOUTH)
      CALL HALO_UPDATE(grid, W, FROM=SOUTH)

      DO 710 J=J_0STG,J_1STG
      UEARTH=RADIUS*OMEGA*COSV(J)
        W4IofK(:) = 0
        UKIofK(:) = 0
       WU4IofK(:) = 0
      WKE4IofK(:) = 0
       FIMIofK(:) = 0
      I=IM
      DO 700 IP1=1,IM
      DO 660 L=1,LM
  660 PLO(L) = .25*(PMID(L,I,J-1) + PMID(L,Ip1,J-1) +
     +              PMID(L,I,J)   + PMID(L,Ip1,J))
      Do 695 K=2,KM
      IF (PM(K) >= .25*(PEDN(1,I,J-1) + PEDN(1,Ip1,J-1) +
     +                  PEDN(1,I,J)   + PEDN(1,Ip1,J)))  GoTo 695
      L=1
      IF (PM(K).GE.PLO(1)) GO TO 680
  670 LUP=L+1
      IF (L.EQ.LM) GO TO 680
      IF (PM(K).GE.PLO(LUP)) GO TO 675
      L=LUP
      GO TO 670
  675 DPUP=PM(K)-PLO(LUP)
      DPDN=PLO(L)-PM(K)
      BYDP=1./(DPDN+DPUP)
      UK=BYDP*(U(I,J,L)*DPUP+U(I,J,LUP)*DPDN)
      VK=BYDP*(V(I,J,L)*DPUP+V(I,J,LUP)*DPDN)
      GO TO 690
C**** SPECIAL CASES;  L=1,L=LM
  680 UK=U(I,J,L)
      VK=V(I,J,L)
C**** MERIDIONAL AVERAGING
  690 W4=.25*(W(I,J-1,K)+W(IP1,J-1,K)+W(I,J,K)+W(IP1,J,K))
        W4IofK(K) =   W4Iofk(K) + W4
        UKIofK(K) =   UKIofk(K) + UK
       WU4IofK(K) =  WU4Iofk(K) + W4*UK
      WKE4IofK(K) = WKE4Iofk(K) + W4*(UK*UK + VK*VK)
       FIMIofK(K) =  FIMIofk(K) + 1
  695 Continue
  700 I=IP1
      Do 710 K=2,KM
      byFIM = 1 / (FIMIofK(K)+teeny)
      WUJK(J,K) = (WU4IofK(K) - W4IofK(K)*UKIofK(K)*byFIM)*byFIM
      AGC(J,K-1,JK_TOTVTKE)  = AGC(J,K-1,JK_TOTVTKE)  + WKE4IofK(K)
      AGC(J,K-1,JK_VTAMEDDY) = AGC(J,K-1,JK_VTAMEDDY) +
     +                         WU4IofK(K) - byFIM*W4IofK(K)*UKIofK(K)
  710 AGC(J,K-1,JK_TOTVTAM)  = AGC(J,K-1,JK_TOTVTAM)  + WU4IofK(K)  !  + W4I*UEARTH
C****
C**** POTENTIAL VORTICITY AND VERTICAL TRANSPORT OF POT. VORT.
C****
      DO 760 J=J_0S,J_1S
      JHEMI=1
      IF (J.LT.1+JM/2) JHEMI=-1
      DO 730 K=1,KM
      PVI=0.
      DO 720 I=1,IM
      DUT(I,J,K)=JHEMI*STB(I,J,K)*(ZX(I,J,K)-FCOR(J))*BYDXYP(J)
  720 PVI=PVI+DUT(I,J,K)
  730 AGC(J,K,JK_POTVORT)=AGC(J,K,JK_POTVORT)+PVI
      DO 760 K=2,KM
      W2I=0.
      PV2I=0.
      WPV4I=0.
      FIMI=0.
      I=IM
      DO 740 IP1=1,IM
      If (PM(K) >= .5*(PEDN(1,I,J)+PEDN(1,Ip1,J)))  GoTo 740
      W2=.5*(W(I,J,K)+W(IP1,J,K))
      W2I=W2I+W2
      PV2=.5*(DUT(I,J,K-1)+DUT(I,J,K))
      PV2I=PV2I+PV2
      WPV4I=WPV4I+W2*PV2
      FIMI=FIMI+1.
  740 I=IP1
      AGC(J,K-1,JK_VTPV)=AGC(J,K-1,JK_VTPV)+WPV4I
  760 AGC(J,K-1,JK_VTPVEDDY)=AGC(J,K-1,JK_VTPVEDDY)+
     &     WPV4I-W2I*PV2I/(FIMI+teeny)
C****
C**** SPECIAL MEAN/EDDY DIAGNOSTICS ARE CALCULATED
C****
      DO 770 J=J_0STG,J_1STG
      DO 765 K=2,KM
      DPE=PMO(K)-PMO(K-1)
      UP(J,K)=(UJK(J,K)-UJK(J,K-1))/DPE
  765 PSIP(J,K)=(PSIJK(J,K)-PSIJK(J,K-1))/DPE
      UP(J,1)=UP(J,2)
      PSIP(J,1)=PSIP(J,2)
  770 CONTINUE
      DO 780 K=1,KM
      KUP=K+1
      IF (K.EQ.KM) KUP=KM
      KDN=K-1
      IF (K.EQ.1) KDN=1

      CALL HALO_UPDATEj(grid, TJK, FROM=SOUTH)
c***      DO L=1,LM
c***         CALL HALO_UPDATE(grid, TJK(:,L), FROM=SOUTH)
c***      END DO

      DO 780 J=J_0STG,J_1STG
      TY(J,K)=(TJK(J,K)-TJK(J-1,K))/DYV(J)
C**** E-P FLUX NORTHWARD COMPONENT
      AGC(J,K,JK_EPFLXNCP)=AGC(J,K,JK_EPFLXNCP)+
     &     PSIJK(J,K)*(UJK(J,KUP)-UJK(J,KDN))/
     *  (PMO(KUP)-PMO(KDN))-UVJK(J,K)
  780 CONTINUE

      CALL HALO_UPDATEj(grid, PSIJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, UJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, VJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, WJK, FROM=NORTH)
      CALL HALO_UPDATEj(grid, UP, FROM=NORTH)
      CALL HALO_UPDATEj(grid, TY, FROM=NORTH)
      CALL HALO_UPDATEj(grid, PSIP, FROM=NORTH)

c***      DO L=1,LM
c***         CALL HALO_UPDATE(grid, PSIJK(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, UJK(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, VJK(:,L), FROM=NORTH)
c***         If (L > 1) THEN
c***           CALL HALO_UPDATE(grid, WJK(:,L), FROM=NORTH)
c***         END IF
c***         CALL HALO_UPDATE(grid, UP(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, TY(:,L), FROM=NORTH)
c***         CALL HALO_UPDATE(grid, PSIP(:,L), FROM=NORTH)
c***      END DO

      DO 800 J=J_0S,J_1S
      DO 800 K=2,KM-1
      UY=(UJK(J+1,K)*DXV(J+1)-UJK(J,K)*DXV(J)-FCOR(J))/DXYP(J)
      PSIY=(PSIJK(J+1,K)*DXV(J+1)-PSIJK(J,K)*DXV(J))/DXYP(J)
C**** ZONAL MEAN MOMENTUM EQUATION   (MEAN ADVECTION)
      AGC(J,K,JK_DUDTMADV)=AGC(J,K,JK_DUDTMADV)-
     &     .5*UY*(VJK(J,K)+VJK(J+1,K))-
     *  .25*((UP(J+1,K+1)+UP(J,K+1))*WJK(J,K+1)+(UP(J+1,K)+UP(J,K))*
     *   WJK(J,K))
C**** ZONAL MEAN HEAT EQUATION   (MEAN ADVECTION)
      AGC(J,K,JK_DTDTMADV)=AGC(J,K,JK_DTDTMADV)-
     &     .5*(TY(J,K)*VJK(J,K)+TY(J+1,K)*VJK(J+1,K))
     *  -.5*STJK(J,K)*(WJK(J,K+1)+WJK(J,K))
C**** LAGRANGIAN MEAN MOMENTUM EQUATION  (MEAN ADVECTION)
      VSTAR=.5*(VJK(J,K)+VJK(J+1,K)-.5*(PSIP(J,K)+PSIP(J,K+1)
     *  +PSIP(J+1,K)+PSIP(J+1,K+1)))
      WSTAR=.5*(WJK(J,K)+WJK(J,K+1))+PSIY
      AGC(J,K,JK_DUDTTEM)=AGC(J,K,JK_DUDTTEM)-
     &     UY*VSTAR-.25*(UP(J,K)+UP(J+1,K)+
     *  UP(J,K+1)+UP(J+1,K+1))*WSTAR
      AGC(J,K,JK_DTDTTEM)=AGC(J,K,JK_DTDTTEM)-
     &     .5*(TY(J+1,K)+TY(J,K))*VSTAR-
     *  STJK(J,K)*WSTAR
C**** VERTICAL E-P FLUX
      AGC(J,K-1,JK_EPFLXVCP)=AGC(J,K-1,JK_EPFLXVCP)-
     &     WUJK(J,K)-.5*PSIJK(J,K)*UY
      AGC(J,K,JK_EPFLXVCP)=AGC(J,K,JK_EPFLXVCP)-.5*PSIJK(J,K)*UY
  800 CONTINUE
C****
C**** SPECTRAL ANALYSIS OF KINETIC ENERGIES AT CONSTANT PRESSURE
C****
      IZERO=0
      NM=1+IM/2
      J45N=2+.75*(JM-1.)
c      KS1=LS1
C**** TOTAL THE KINETIC ENERGIES
      KE(:,:)=0.
      KE_part(:,:,:)=0.

C P already halo'ed; no need      CALL CHECKSUM(grid, P, __LINE__, __FILE__)
C P already halo'ed; no need     CALL HALO_UPDATE(grid, P, FROM=SOUTH)

      DO J=J_0STG,J_1STG
        DO K=1,KM
          KSPHER=KLAYER(K)
          IF (J.GT.JEQ) KSPHER=KSPHER+1
          DO KX=IZERO,LM,LM
            DO I=1,IM
              Ip1=I+1  ;  If(I==IM) Ip1=1
              DPUV=0.
              DO 2025 L=1,LM
              PLO(L) = (PMID(L,I,J  ) + PMID(L,Ip1,J  ))*RAPVS(J) +
     *                 (PMID(L,I,J-1) + PMID(L,Ip1,J-1))*RAPVN(J-1)
 2025         PL(L)  = (PEDN(L,I,J  ) + PEDN(L,Ip1,J  ))*RAPVS(J) +
     *                 (PEDN(L,I,J-1) + PEDN(L,Ip1,J-1))*RAPVN(J-1)
              IF (PM(K+1).GE.PLO(1)) GO TO 2090           ! really ?? not PL?
              L=1
              PDN = PL(1)
              IF (PM(K).GE.PLO(1)) GO TO 2040             ! really ?? not PL?
              PDN=PM(K)
 2030         IF (PM(K).GT.PL(L+1)) GO TO 2040
              L=L+1
              If (L <= LM)  GoTo 2030
 2040         LUP=L
 2050         IF (PM(K+1).GE.PL(LUP+1)) GO TO 2060
              LUP=LUP+1
              If (LUP <= LM)  GoTo 2050
 2060         CONTINUE
C**** ACCUMULATE HERE
              SQRTDP=SQRT(PDN-PM(K+1))
 2070         PUP=PL(L+1)
              IF (LUP.EQ.L) PUP=PM(K+1)
              DP=PDN-PUP
              IF(KX.EQ.IZERO) DPUV=DPUV+DP*U(I,J,L)
              IF(KX.EQ.LM)    DPUV=DPUV+DP*V(I,J,L)
              IF (LUP.EQ.L) GO TO 2080
              L=L+1
              PDN=PL(L)
              If (L <= LM)  GoTo 2070
 2080         IF (SQRTDP.EQ.0.) SQRTDP=teeny
              DPUV=DPUV/SQRTDP
 2090         X1(I)=DPUV
            ENDDO
            CALL FFTE (X1,X1tmp)
            X1=X1tmp
            IF (J.NE.JEQ) THEN
              DO N=1,NM
                KE_part(N,J,KSPHER)=KE_part(N,J,KSPHER)+X1(N)*DXYV(J)
              ENDDO
              IF (J.EQ.J45N) THEN
                DO N=1,NM
                  KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                  X1(N)*DXYV(J)
                ENDDO
              ENDIF
            ELSE
              DO N=1,NM
                KE_part(N,J,KSPHER+2)=KE_part(N,J,KSPHER+2)+
     &                                X1(N)*DXYV(J)
                KE_part(N,J,KSPHER  )=KE_part(N,J,KSPHER)+
     &                                .5D0*X1(N)*DXYV(J)
                KE_part(N,J,KSPHER+1)=KE_part(N,J,KSPHER+1)+
     &                                .5D0*X1(N)*DXYV(J)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO

c      CALL GLOBALSUM(grid, KE_part, KE, ALL=.true.) ! uses transposes
      KE_jsum = sum(KE_part(:,J_0:J_1,:),2)
      CALL SUMXPE(KE_jsum, KE)
! changing the PE count causes roundoff diffs when using sumxpe, so
      if(am_i_root()) call reduce_precision(ke,1d-9) ! remove trailing bits

      call broadcast(grid, KE)

      DO 2150 KS=1,NSPHER
      DO 2150 N=1,NM
 2150 SPECA(N,18,KS)=SPECA(N,18,KS)+KE(N,KS)
C**** ACCUMULATE TIME USED IN DIAGA
      CALL TIMEOUT(BEGIN,MDIAG,MDYN)
      RETURN
      END SUBROUTINE DIAGB


      SUBROUTINE DIAG5A (M5,NDT)
C****
C**** THIS DIAGNOSTICS ROUTINE PRODUCES A SPECTRAL ANALYSIS OF KINETIC
C**** AND AVAILABLE POTENTIAL ENERGIES AND THEIR TRANSFER RATES BY
C**** VARIOUS ATMOSPHERIC PROCESSES.
C****
C**** THE PARAMETER M INDICATES WHAT IS STORED IN SPECA(N,M,KSPHER),
C**** IT ALSO INDICATES WHEN DIAG5A IS BEING CALLED.
C**** M=1  MEAN STANDING KINETIC ENERGY            BEFORE SOURCES
C****   2  MEAN KINETIC ENERGY                     BEFORE DYNAMICS
C****   3  MEAN POTENTIAL ENERGY
C****   4  CONVERSION OF K.E. BY ADVECTION         AFTER ADVECTION
C****   5  CONVERSION OF K.E. BY CORIOLIS FORCE    AFTER CORIOLIS TERM
C****   6  CONVERSION FROM P.E. INTO K.E.          AFTER PRESS GRAD FORC
C****   7  CHANGE OF K.E. BY DYNAMICS              AFTER DYNAMICS
C****   8  CHANGE OF P.E. BY DYNAMICS
C****   9  CHANGE OF K.E. BY CONDENSATION          AFTER CONDENSATION
C****  10  CHANGE OF P.E. BY CONDENSATION
C****  11  CHANGE OF P.E. BY RADIATION             AFTER RADIATION
C****  12  CHANGE OF K.E. BY SURFACE               AFTER SURFACE
C****  13  CHANGE OF P.E. BY SURFACE
C****  14  CHANGE OF K.E. BY FILTER                AFTER FILTER
C****  15  CHANGE OF P.E. BY FILTER
C****  16  CHANGE OF K.E. BY DAILY                 AFTER DAILY
C****  17  CHANGE OF P.E. BY DAILY
C****  18  UNUSED
C****  19  LAST KINETIC ENERGY
C****  20  LAST POTENTIAL ENERGY
C****
      USE CONSTANT, only : sha, UNDEF_VAL
      Use RESOLUTION, Only: IM,JM,LM
      USE MODEL_COM, only : IDACC,MDIAG
      Use ATM_COM,    Only: PMID,PEDN,PDSIG,PK,T,U,V,ZATMO
      USE GEOM, only : AREAG,DXYN,DXYP,DXYS,imaxj
      USE GC_COM, only : speca,atpe,nspher,kspeca,klayer,agc=>agc_loc
     &     ,jeq
      USE DIAG_COM, only : SQRTM
      USE DIAG_COM, only : imh,fim,byim
      USE GCDIAG, only : jl_ape
      USE DIAG_LOC, only : lupa,ldna
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : getDomainBounds
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, AM_I_ROOT
      USE DOMAIN_DECOMP_1D, only : GLOBALSUM, SOUTH, WRITE_PARALLEL
      USE DOMAIN_DECOMP_1D, only : SUMXPE, broadcast
      USE PRECISION_MOD
      IMPLICIT NONE
      INTEGER :: M5,NDT

!**** Local variables
      REAL*8, DIMENSION(IM) :: X, Xtmp
      REAL*8, DIMENSION(IMH+1,NSPHER) :: KE,KE_jsum,APE
      REAL*8, DIMENSION
     &  (IMH+1,GRID%J_STRT_HALO:GRID%J_STOP_HALO,NSPHER) :: KE_part
      REAL*8, DIMENSION(IMH+1,4,LM) :: VAR,VAR_jsum
      REAL*8, DIMENSION(IMH+1,4,LM,GRID%J_STRT_HALO:GRID%J_STOP_HALO)
     &     :: VAR_part
      REAL*8, DIMENSION(2) :: TPE
      REAL*8               :: TPE_sum
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO) :: TPE_psum
CMoved to DAGCOM so it could be declared allocatable      REAL*8, SAVE, DIMENSION(IM,JM) :: SQRTM

      REAL*8, DIMENSION(LM) :: THGM,GMEAN

      INTEGER, PARAMETER :: IZERO=0

      INTEGER, DIMENSION(KSPECA), PARAMETER ::
     &     MTPEOF=(/0,0,1,0,0,0,0,2,0,3,  4,0,5,0,6,0,7,0,0,8/)

      INTEGER :: I,IJL2,IP1,J,J45N,JH,JHEMI,JP,K,KS,KSPHER,L,LDN,
     &     LUP,MAPE,MKE,MTPE,N,NM

      REAL*8 :: SUMI,SUMT,NOW
      REAL*8, DIMENSION(grid%J_STRT_HALO:grid%J_STOP_HALO,LM) ::
     &     GMEAN_part
      REAL*8, DIMENSION(GRID%J_STRT_HALO:GRID%J_STOP_HALO,LM) ::
     &        THGM_part

      INTEGER :: J_0S,J_1S,J_0STG,J_1STG,J_0,J_1
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE

#ifdef SCM
      return
#endif

      TPE_psum([grid%J_STRT_HALO,grid%J_STOP_HALO]) = UNDEF_VAL
      TPE = 0.0d0 

      call getDomainBounds(GRID, J_STRT_SKP=J_0S   , J_STOP_SKP=J_1S,
     &               J_STRT_STGR=J_0STG, J_STOP_STGR=J_1STG,
     &               HAVE_SOUTH_POLE=HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE=HAVE_NORTH_POLE)
      J_0=GRID%J_STRT
      J_1=GRID%J_STOP

      NM=1+IM/2
      J45N=2.+.75*(JM-1.)
      IJL2=IM*JM*LM*2

      MKE=M5
      MAPE=M5
C****
C**** Note: KSPHER has been re-arranged from previous models to better
C****       deal with optionally resolved stratosphere. The higher
C****       values are only used if the model top is high enough.
C****
C**** KSPHER=1 SOUTHERN TROPOSPHERE         2 NORTHERN TROPOSPHERE
C****        3 EQUATORIAL TROPOSPHERE       4 45 DEG NORTH TROPOSPHERE
C****
C****        5 SOUTHERN LOW STRATOSPHERE    6 NORTHERN LOW STRATOSPHERE
C****        7 EQUATORIAL LOW STRATOSPHERE  8 45 DEG NORTH LOW STRATOSPH
C****
C****        9 SOUTHERN MID STRATOSPHERE   10 NORTHERN MID STRATOSPHERE
C****       11 EQUATORIAL MID STRATOSPHERE 12 45 DEG NORTH MID STRATOSPH
C****
C****       13 SOUTHERN UPP STRATOSPHERE   14 NORTHERN UPP STRATOSPHERE
C****       15 EQUATORIAL UPP STRATOSPHERE 16 45 DEG NORTH UPP STRATOSPH
C****
      GO TO (200,200,810,810,810,  810,200,810,205,810,
     *       296,205,810,205,810,  205,810,810,810,810),M5

C***  810 WRITE (6,910) M5
  810 CALL WRITE_PARALLEL(M5, UNIT=6, format=
     & "('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.  M5=',I5)")
C****  910 FORMAT ('0INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.  M5=',I5)
      call stop_model('INCORRECT VALUE OF M5 WHEN CALLING DIAG5A.',255)
C**** MASS FOR KINETIC ENERGY
  200 CONTINUE

C****
  205 CONTINUE

      MAPE=MKE+1
      KE(:,:)=0.
      KE_part(:,:,:)=0.
C**** CURRENT KINETIC ENERGY
      DO L=1,LM
        DO J=J_0STG,J_1STG
          I=IM
          Do Ip1=1,IM
             SQRTM(I,J) =
     =          Sqrt(.5*((PDSIG(L,I,J  )+PDSIG(L,Ip1,J  ))*DXYS(J) +
     +                   (PDSIG(L,I,J-1)+PDSIG(L,Ip1,J-1))*DXYN(J-1)))
             I=Ip1  ;  EndDo
          IF (J <= JEQ) THEN
            KSPHER=KLAYER(L)
          ELSE
            KSPHER=KLAYER(L)+1
          END IF
          DO K = IZERO,LM,LM
            IF(K.EQ.IZERO)X(1:IM)=U(1:IM,J,L)*SQRTM(1:IM,J)
            IF(K.EQ.LM)   X(1:IM)=V(1:IM,J,L)*SQRTM(1:IM,J)
            CALL FFTE (X,Xtmp)
            X=Xtmp
            IF (J.EQ.JEQ) THEN
              DO N=1,NM
                KE_part(N,J,KSPHER+2) = KE_part(N,J,KSPHER+2) + X(N)
                KE_part(N,J,KSPHER  ) = KE_part(N,J,KSPHER  ) + X(N)*.5
                KE_part(N,J,KSPHER+1) = KE_part(N,J,KSPHER+1) + X(N)*.5
              ENDDO
cgsfc              IF(K.EQ.LM)KSPHER=KSPHER+1
            ELSE
              DO N=1,NM
                KE_part(N,J,KSPHER) = KE_part(N,J,KSPHER) + X(N)
              ENDDO
              IF (J.EQ.J45N) THEN
                DO N=1,NM
                  KE_part(N,J,KSPHER+2) = KE_part(N,J,KSPHER+2) + X(N)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO

c      CALL GLOBALSUM(grid, KE_part, KE, ALL=.true.) ! uses transposes
      KE_jsum = sum(KE_part(:,J_0:J_1,:),2)
      CALL SUMXPE(KE_jsum, KE)
! changing the PE count causes roundoff diffs when using sumxpe, so
      if(am_i_root()) call reduce_precision(ke,1d-9) ! remove trailing bits
      call broadcast(grid, KE)

      IF (NDT /= 0) THEN
C**** TRANSFER RATES AS DIFFERENCES OF KINETIC ENERGY
        DO KS=1,NSPHER
          DO N=1,NM
           SPECA(N,MKE,KS)=SPECA(N,MKE,KS)+(KE(N,KS)-SPECA(N,19,KS))/NDT
          END DO
        END DO
      END IF
      DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,19,KS)=KE(N,KS)
        END DO
      END DO

C****
C**** POTENTIAL ENERGY
C****
  296 CONTINUE


C****
C**** AVAILABLE POTENTIAL ENERGY
C****
C**** Calculate global means for each layer of
C**** pot. temp (thgm) and static stability (gmean)
C****
      DO L=1,LM
        LDN=LDNA(L)
        LUP=LUPA(L)
        DO J=J_0,J_1
          GMEAN_part(J,L)=0.
          THGM_part(J,L)=0.
          DO I=1,IMAXJ(J)
            THGM_part(J,L) = THGM_part(J,L) +T(I,J,L)*Sqrt(PDSIG(L,I,J))
            GMEAN_part(J,L)=GMEAN_part(J,L)+
     +                           PMID(L,I,J)*(T(I,J,LUP)-T(I,J,LDN)) /
     /                         (PK(L,I,J)*(PMID(LDN,I,J)-PMID(LUP,I,J)))
          ENDDO
          GMEAN_part(J,L)=GMEAN_part(J,L)*DXYP(J)
          THGM_part(J,L)=THGM_part(J,L)*DXYP(J)
        ENDDO
        IF(HAVE_SOUTH_POLE) THEN
          THGM_part(1,L)=THGM_part(1,L)*FIM
          GMEAN_part(1,L)=GMEAN_part(1,L)*FIM
        ENDIF
        IF(HAVE_NORTH_POLE) THEN
          THGM_part(JM,L)=THGM_part(JM,L)*FIM
          GMEAN_part(JM,L)=GMEAN_part(JM,L)*FIM
        ENDIF
      ENDDO

      CALL GLOBALSUM(grid,THGM_part(:,1:LM),THGM(1:LM),ALL=.TRUE.)
      THGM=THGM/AREAG
      CALL GLOBALSUM(grid,GMEAN_part(:,1:LM),GMEAN(1:LM),ALL=.TRUE.)

      GMEAN(1:LM) = AREAG / (GMEAN(1:LM) + 1d-20)

      APE(:,:)=0.

C**** SPECTRAL ANALYSIS OF AVAILABLE POTENTIAL ENERGY

      DO L = 1, LM
        LDN=LDNA(L)
        LUP=LUPA(L)

        DO J=J_0,J_1
          VAR_part(:,:,L,J)=0
          IF (J < JEQ) THEN
            JHEMI = 1
          ELSE
            JHEMI = 2
          END IF

          DO I=1,IM
             X(I) = T(I,J,L)*Sqrt(PDSIG(L,I,J)) - THGM(L)
          END DO

          if(m5.eq.7) then
            AGC(J,L,JL_APE)=AGC(J,L,JL_APE)+SUM(X*X)*GMEAN(L)*BYIM*BYIM
          endif

          IF(J.EQ.1 .or. J.EQ.JM) THEN
            VAR_part(1,JHEMI,L,J)=.5*(X(1)**2)*DXYP(1)*FIM
            cycle ! skip poles for spectral analysis
          END IF

          CALL FFTE (X,Xtmp)
          X=Xtmp
          VAR_part(1:NM,JHEMI,L,J)=X(1:NM)*DXYP(J)
          IF (J == JEQ-1) THEN
            VAR_part(1:NM,3,L,J)=X(1:NM)*DXYP(J)
          ELSEIF (J == J45N-1) THEN
            VAR_part(1:NM,4,L,J)=X(1:NM)*DXYP(J)
          END IF
        END DO ! J
      END DO ! L

c      CALL GLOBALSUM(grid, VAR_part, VAR) ! not parallelized
      VAR_jsum = sum(VAR_part(:,:,:,J_0:J_1),4)
      CALL SUMXPE(VAR_jsum, VAR)

      IF (AM_I_ROOT()) THEN
        DO L = 1, LM
          KS=KLAYER(L)
          DO JHEMI=1,4
            DO N=1,NM
              APE(N,KS)=APE(N,KS)+VAR(N,JHEMI,L)*GMEAN(L)
            END DO
            KS=KS+1
          END DO
        END DO
! changing the PE count causes roundoff diffs when using sumxpe, so
        call reduce_precision(ape,1d-9) ! remove trailing bits
      END IF

C**** CURRENT TOTAL POTENTIAL ENERGY
 450  CONTINUE

      IF (HAVE_SOUTH_POLE) THEN
        J=1
        SUMT = Sum(T(1,J,:)*PK(:,1,J)*PDSIG(:,1,J))
        TPE_psum(J) = FIM*DXYP(J)*(ZATMO(1,J)*PEDN(1,1,J) + SUMT*SHA)
      END IF
      IF (HAVE_NORTH_POLE) THEN
        J=JM
        SUMT = Sum(T(1,J,:)*PK(:,1,J)*PDSIG(:,1,J))
        TPE_psum(J) = FIM*DXYP(J)*(ZATMO(1,J)*PEDN(1,1,J) + SUMT*SHA)
      END IF
      DO J=J_0S, J_1S
        SUMI=0
        DO I=1,IM
          SUMT = Sum(T(I,J,:)*PK(:,I,J)*PDSIG(:,I,J))
          SUMI = SUMI + ZATMO(I,J)*PEDN(1,I,J) + SUMT*SHA
        END DO
        TPE_psum(J) = SUMI*DXYP(J)
      END DO
      CALL GLOBALSUM(grid, TPE_psum, TPE_sum, TPE)

      IF (NDT /= 0) THEN
        MTPE=MTPEOF(MAPE)

C**** TRANSFER RATES AS DIFFERENCES FOR POTENTIAL ENERGY
        DO KS=1,NSPHER
          DO N=1,NM
        SPECA(N,MAPE,KS)=SPECA(N,MAPE,KS)+(APE(N,KS)-SPECA(N,20,KS))/NDT
          END DO
        END DO

        ATPE(MTPE,1)=ATPE(MTPE,1)+(TPE(1)-ATPE(8,1))/NDT
        ATPE(MTPE,2)=ATPE(MTPE,2)+(TPE(2)-ATPE(8,2))/NDT
      END IF
      DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,20,KS)=APE(N,KS)
        END DO
      END DO
      ATPE(8,1)=TPE(1)
      ATPE(8,2)=TPE(2)

      IF (M5.EQ.2) THEN
C**** ACCUMULATE MEAN KINETIC ENERGY AND MEAN POTENTIAL ENERGY
        DO KS=1,NSPHER
        DO N=1,NM
          SPECA(N,2,KS)=SPECA(N,2,KS)+KE(N,KS)
          SPECA(N,3,KS)=SPECA(N,3,KS)+APE(N,KS)
        END DO
        END DO
        ATPE(1,1)=ATPE(1,1)+TPE(1)
        ATPE(1,2)=ATPE(1,2)+TPE(2)
      END IF

      CALL TIMER (NOW,MDIAG)
      RETURN
      END SUBROUTINE DIAG5A


      SUBROUTINE DIAG7A
C****
C**** THIS ROUTINE ACCUMULATES A TIME SEQUENCE FOR SELECTED
C**** QUANTITIES AND FROM THAT PRINTS A TABLE OF WAVE FREQUENCIES.
C****
      USE CONSTANT, only : grav,bygrav
      USE RESOLUTION, only : im,jm,lm
      USE MODEL_COM, only : IDACC,MDIAG
      USE ATM_COM, only : PMID,U,V,PHI
      USE GC_COM, only : nwav_dag,wave,max12hr_sequ,kwp,re_and_im
     &     ,jeq,j50n
      USE DIAG_COM, only : ia_12hr,imh
      USE DIAG_LOC, only : ldex
      USE DOMAIN_DECOMP_ATM, only : GRID
      USE DOMAIN_DECOMP_1D, only : getDomainBounds,SUMXPE,AM_I_ROOT
      IMPLICIT NONE

      REAL*8, DIMENSION(0:IMH) :: AN,BN
      INTEGER, PARAMETER :: KM=6,KQMAX=12
      INTEGER :: NMAX=nwav_dag
      REAL*8, DIMENSION(IM,KM) :: HTRD,HTRD_loc
      REAL*8, DIMENSION(IM,LM) :: UEQ,UEQ_loc,VEQ,VEQ_loc
      REAL*8, DIMENSION(KM), PARAMETER ::
     &     PMB=(/922.,700.,500.,300.,100.,10./),
     &     GHT=(/500.,2600.,5100.,8500.,15400.,30000./)
      REAL*8 :: PIJ50N,PL,PLM1,SLOPE,NOW
      INTEGER I,IDACC9,K,KQ,L,N
      INTEGER :: J_0, J_1

#ifdef SCM
      return
#endif

      call getDomainBounds(GRID,J_STRT=J_0,J_STOP=J_1)

      IDACC9=IDACC(ia_12hr)+1
      IDACC(ia_12hr)=IDACC9
      IF (IDACC9.GT.Max12HR_sequ) RETURN

      IF(J_0 <= JEQ .and. JEQ <= J_1) THEN
        UEQ_loc=U(:,JEQ,:)
        VEQ_loc=V(:,JEQ,:)
      ELSE
        UEQ_loc=0d0
        VEQ_loc=0d0
      ENDIF

      CALL SUMXPE(UEQ_loc, UEQ)
      CALL SUMXPE(VEQ_loc, VEQ)

      IF(AM_I_ROOT()) THEN
        DO KQ=1,3
          CALL FFT(UEQ(1,LDEX(KQ)),AN,BN)
          DO N=1,NMAX
            WAVE(1,IDACC9,N,2*KQ-1)=AN(N)
            WAVE(2,IDACC9,N,2*KQ-1)=BN(N)
          ENDDO
          CALL FFT(VEQ(1,LDEX(KQ)),AN,BN)
          DO N=1,NMAX
            WAVE(1,IDACC9,N,2*KQ)=AN(N)
            WAVE(2,IDACC9,N,2*KQ)=BN(N)
          ENDDO
        ENDDO
      ENDIF

      IF(J_0 <= J50N .and. J50N <= J_1) THEN
        DO 150 I=1,IM
          K=1
          L=1
          PL = PMID(1,I,J50N)
 130      L=L+1
          PLM1=PL
          PL = PMID(L,I,J50N)
          IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
C**** ASSUME THAT PHI IS LINEAR IN LOG P
          SLOPE=(PHI(I,J50N,L-1)-PHI(I,J50N,L))/LOG(PLM1/PL)
 140      HTRD_loc(I,K)=
     &         (PHI(I,J50N,L)+SLOPE*LOG(PMB(K)/PL))*BYGRAV-GHT(K)
          IF (K.GE.KM) GO TO 150
          K=K+1
          IF (PMB(K).LT.PL.AND.L.LT.LM) GO TO 130
          GO TO 140
 150    CONTINUE
      ELSE
        HTRD_loc(:,:)=0d0
      ENDIF

      CALL SUMXPE(HTRD_loc, HTRD)

      IF(AM_I_ROOT()) THEN
        DO KQ=7,KQMAX
          CALL FFT(HTRD(1,KQ-6),AN,BN)
          DO N=1,NMAX
            WAVE(1,IDACC9,N,KQ)=AN(N)
            WAVE(2,IDACC9,N,KQ)=BN(N)
          END DO
        END DO
      ENDIF

      CALL TIMER (NOW,MDIAG)
      RETURN
      END SUBROUTINE DIAG7A


      subroutine speca_prep
!@sum speca_prep apply scale factors to SPECA and ATPE accumulations
!@+   to promote model independence of offline ASCII-table printer routines.
!@+   Accumulations are also rescaled by the number of physics timesteps
!@+   to avoid passing detailed idacc-info ( ia_d5[dfs], ia_filt, ia_12hr )
      use constant, only : grav,rgas,teeny
      use dynamics, only : dt
      use model_com, only : idacc
      use gc_com, only : speca,speca_out,kspeca
     &     ,atpe,atpe_out,ktpe
      use diag_com, only :
     &     ia_src,ia_d5s,ia_filt,ia_12hr,ia_d5f,ia_d5d,ia_dga,ia_inst
      implicit none
      integer :: m,hemis
      real*8, dimension(kspeca) :: scalet
      integer, dimension(ktpe), parameter ::
     &     mapeof=(/3,8,10,11,13,15,17,20/)


      ! todo: remove powers of 10 from the scale factors
      IF (IDACC(ia_inst).LT.1) IDACC(ia_inst)=1
      SCALET(1)=100.D-17/(GRAV*IDACC(ia_dga)+teeny)
      SCALET(19)=100.D-17/(GRAV*IDACC(ia_inst))
      SCALET(20)=SCALET(19)*RGAS
      SCALET(2)=SCALET(19)*IDACC(ia_inst)/(IDACC(ia_d5d)+teeny)
      SCALET(3)=SCALET(2)*RGAS
      SCALET(4)=100.D-12/(GRAV*DT*IDACC(ia_d5f)+teeny)
      SCALET(5)=SCALET(4)
      SCALET(6)=SCALET(4)
      SCALET(7)=100.D-12/(GRAV*DT*(IDACC(ia_d5d)+teeny))
      SCALET(8)=SCALET(7)*RGAS
      SCALET(9)=100.D-12/(GRAV*DT*(IDACC(ia_d5s)+teeny))
      SCALET(10)=SCALET(9)*RGAS
      SCALET(11)=SCALET(10)
      SCALET(12)=SCALET(9)
      SCALET(13)=SCALET(10)
      SCALET(14)=100.D-12/(GRAV*DT*(IDACC(ia_filt)+teeny))
      SCALET(15)=SCALET(14)*RGAS
      SCALET(16)=100.D-12/(GRAV*DT*(.5*IDACC(ia_12hr)+teeny))
      SCALET(17)=SCALET(16)*RGAS
      SCALET(18)=100.D-17/(GRAV*IDACC(ia_dga)+teeny)

      do m=1,kspeca
        speca_out(:,m,:) = speca(:,m,:) * scalet(m) * idacc(ia_src)
      enddo

      do hemis=1,2
        atpe_out(:,hemis) = atpe(:,hemis)*(scalet(mapeof(:))/rgas)
     &       * idacc(ia_src)
      enddo

      end subroutine speca_prep


      subroutine diaggc_prep
c Calculate derived GC outputs
      use constant, only : kapa,lhe,sha,radius,omega,rgas,tf,teeny
      use resolution, only : im,jm,lm
      use model_com, only : idacc,dtsrc
      use dynamics, only : dt,dsig
      use domain_decomp_atm, only : am_i_root
      use diag_com, only : kdiag, ia_dga,ndaa,
     &     plm, pme=>ple_dn, ple, byim
      use gc_com, only : kagc,jgrid_gc,hemis_gc,vmean_gc,
     &     agc_in=>agc, agc=>agc_out, kep, kagcx, denom_gc
      use gcdiag
      use geom, only : dxyp,dxyv,dxv,bydxyp,cosv,cosp,imaxj,fcor
      implicit none
      integer :: j,j1,j2,k,l,jg,ind,n,ldn,lup,kk
      real*8, dimension(jm,2) :: wtj
      integer, dimension(20) :: inds
      real*8, dimension(jm) :: bydapo,dxcosv
      real*8, dimension(lm) :: pkm
      real*8, dimension(jm,lm) :: temp,vx,sheth,shethe
      real*8 :: wtup,tedge,rhoe,delpe
      REAL*8, DIMENSION(JM,LM) :: ! outputs of EPFLXP
     &     DUDS,DMF,DEF,DMFR,DEFR,ER1,ER2

      if(.not.am_i_root()) return

      bydapo(:)=bydxyp(:)
      bydapo(1)=bydapo(1)*byim
      bydapo(jm)=bydapo(jm)*byim
      dxcosv(:) = dxv(:)*cosv(:)

      pkm(1:lm) = plm(1:lm)**kapa

      agc(:,:,:) = agc_in(:,:,:)

c stream function
      do j=1,jm
        agc(j,lm,jk_psi) = 0d0
        do l=lm-1,1,-1
          agc(j,l,jk_psi)=agc(j,l+1,jk_psi)-agc(j,l+1,jk_v)*dxv(j)
        enddo
      enddo

c
c conversions of northward transports from per unit length to
c integrals around each latitude circle
c
      inds(1:6) = (/
     &     jk_eddntsh, jk_eddntgeo, jk_totntdse,
     &     jk_eddntlh, jk_totntlh,  jk_totntke /)
      do k=1,6
        ind = inds(k)
        do j=2,jm
          agc(j,:,ind) = agc(j,:,ind)*dxv(j)
        enddo
      enddo

c
c conversions to per unit area
c
      inds(1:2) = (/ JL_totvtlh, JL_zmfvtlh /)
      do k=1,2
        ind = inds(k)
        do l=1,lm
          agc(:,l,ind) = agc(:,l,ind)*bydapo(:)
        enddo
      enddo

c
c multiplications by p**kapa
c
      inds(1:2) = (/ jk_dtdtmadv, jk_dtdttem /)
      do k=1,2
        ind = inds(k)
        do l=1,lm
          agc(:,l,ind) = agc(:,l,ind)*pkm(l)
        enddo
      enddo

c
c Denominators which are zonally integrated p(:,:)-ptop (or psf-ptop)
c
      do j=1,jm
        n = jl_dpsig
        agc(j,:,n) = sum(agc(j,:,jk_dpa))/imaxj(j)
        n = jl_dpsigv
        agc(j,:,n) = sum(agc(j,:,jl_dpb))
      enddo

c
c Transformed Eulerian Mean diagnostics
c
C****
C**** SHETH = V'TH'/DTH/DP
C**** SHETHE = 4*SHETH AT INTERFACES x DXV
C**** VX = V-V* IS D/DP(SHETH)
C****
      sheth = agc(:,:,jk_sheth)/(agc(:,:,jk_dpb)+teeny)

      shethe = 0.
      DO J=2,JM
        SHETHE(J,LM)=SHETH(J,LM)
        DO L=1,LM-1
          WTUP=AGC(J,L,JK_DPB)/(AGC(J,L+1,JK_DPB)+AGC(J,L,JK_DPB)+teeny)
          SHETHE(J,L)=(SHETH(J,L)*(1.-WTUP)+SHETH(J,L+1)*WTUP)
        ENDDO
        SHETHE(J,:) = SHETHE(J,:)*DXV(J)*IM*IDACC(ia_dga)
      ENDDO
      n = jk_psi_tem
      agc(:,:,n) = agc(:,:,jk_psi) + shethe

      vx = 0.
      DO J=2,JM
        LDN=1
        DO L=1,LM
          LUP=L+1
          IF (L.EQ.LM) LUP=LM
          VX(J,L)=0.
          IF (AGC(J,L,JK_DPB).gt.0.) then
            IF (AGC(J,LDN,JK_DPB).EQ.0.) LDN=LDN+1
            VX(J,L)=AGC(J,L,JK_DPB)*(
     &           AGC(J,LUP,JK_SHETH)/AGC(J,LUP,JK_DPB)-
     &           AGC(J,LDN,JK_SHETH)/AGC(J,LDN,JK_DPB))/
     &           (PLE(LUP)-PLE(LDN) +.5*
     &           (AGC(J,LUP,JK_DPB)/AGC(J,LUP,JK_NPTSAVG)-
     &            AGC(J,LDN,JK_DPB)/AGC(J,LDN,JK_NPTSAVG)))
          endif
          LDN=L
        ENDDO
      ENDDO
      n = jk_Vstar
      agc(:,:,n) = agc(:,:,jk_v) - vx

      n = jk_Wstar
      agc(:,:,n) = 0.
      DO L=2,LM-1
        DO J=2,JM-1
          AGC(J,L-1,N)=-(DXV(jm/2)/(IM*DXYV(J)))*
     &         (AGC(J+1,L,JK_PSI_TEM)/DXV(J+1)
     &         -AGC(J  ,L,JK_PSI_TEM)/DXV(J  ))
        ENDDO
      ENDDO

c
c NOR. TRANSPORT OF QUASI-GEOSTROPHIC POT. VORTICITY BY EDDIES
c
      n = jk_nt_eqgpv
      agc(:,:,n) = 0.
      DO L=1,LM
      DO J=2,JM-1
        AGC(J,L,N) = AGC(J,L,JK_DPA)*BYDXYP(J)*(
     &   (AGC(J,  L,JK_EDDNTMOM)*DXCOSV(J  )/(AGC(J,  L,JK_DPB)+teeny)-
     &    AGC(J+1,L,JK_EDDNTMOM)*DXCOSV(J+1)/(AGC(J+1,L,JK_DPB)+teeny)
     &       )/COSP(J)
     &      + FCOR(J)*(VX(J,L)+VX(J+1,L))/
     &     (AGC(J,L,JK_DPB)+AGC(J+1,L,JK_DPB)+teeny) )
      ENDDO
      ENDDO

c
c some eddy northward transports
c
      n = jk_nt_dse_e
      agc(:,:,n) = sha*agc(:,:,jk_eddntsh)+agc(:,:,jk_eddntgeo)
      n = jk_nt_see
      agc(:,:,n) = agc(:,:,jk_nt_dse_e)+lhe*agc(:,:,jk_eddntlh)
      n = jk_tot_nt_se
      agc(:,:,n) = agc(:,:,jk_totntdse)+lhe*agc(:,:,jk_totntlh)
      n = jl_nt_lh_e
      agc(:,:,n) = agc(:,:,jl_totntlh)-agc(:,:,jl_zmfntlh)
      n = jk_totntmom
      do l=1,lm
        do j=2,jm
          agc(j,l,n) = agc(j,l,n) + radius*omega*cosv(j)*agc(j,l,jk_v)
        enddo
      enddo

c
c some eddy vertical transports
c
      n = jl_vt_lh_e
      agc(:,:,n)=agc(:,:,Jl_totvtlh)-agc(:,:,Jl_zmfvtlh)
      n = jk_vt_se_eddy
      agc(:,:,n) = agc(:,:,jk_eddvtdse)+lhe*agc(:,:,jk_eddvtlh)
      n = jk_tot_vt_se
      agc(:,:,n) = agc(:,:,jk_totvtdse)+lhe*agc(:,:,jk_totvtlh)

c
c tendency terms in the zonal momentum equation
c
      IF (KEP.gt.0) THEN  ! detailed breakdown from STRAT_DIAG package
        CALL EPFLXP(.false.,DUDS,DMF,DEF,DMFR,DEFR,ER1,ER2)
        agc(:,:,jk_dudt_sum1)     = duds*idacc(ia_dga)
        agc(:,:,jk_dudt_meanadv)  = dmf *idacc(ia_dga)
        agc(:,:,jk_dudt_eddycnv)  = def *idacc(ia_dga)
        agc(:,:,jk_dudt_trnsadv)  = dmfr*idacc(ia_dga)
        agc(:,:,jk_dudt_epflxdiv) = defr*idacc(ia_dga)
        agc(:,:,jk_dudt_fderr1)   = er1 *idacc(ia_dga)
        agc(:,:,jk_dudt_fderr2)   = er2 *idacc(ia_dga)
      ELSE                ! tendencies from approximate eddy fluxes
        agc(:,:,jk_dudt_econv) = 0.
        agc(:,:,jk_dudt_epdiv) = 0.
        do l=2,lm-1
        do j=2,jm-1
          if(agc(j,l,jk_dpa).gt.0.) then
            n = jk_dudt_econv
            agc(j,l,n) =
     &           (AGC(J  ,L,JK_EDDNTMOM)*DXV(J)
     &           -AGC(J+1,L,JK_EDDNTMOM)*DXV(J+1))/DXYP(J)
     &           +.5*(  ! note using DPA as denom, not avg of DPB
     &           AGC(J  ,L,JK_VTAMEDDY)-AGC(J  ,l-1,JK_VTAMEDDY)
     &          +AGC(J+1,L,JK_VTAMEDDY)-AGC(J+1,l-1,JK_VTAMEDDY))
          endif
          n = jk_dudt_epdiv
          agc(j,l,n) =
     &         ( AGC(J+1,L,JK_EPFLXNCP)*DXCOSV(J+1)-
     &           AGC(J  ,L,JK_EPFLXNCP)*DXCOSV(J) )/(DXYP(J)*COSP(J))
     &         +.5*(
     &           AGC(J  ,l-1,JK_EPFLXVCP)-AGC(J  ,L,JK_EPFLXVCP)
     &          +AGC(J+1,l-1,JK_EPFLXVCP)-AGC(J+1,L,JK_EPFLXVCP))/
     &         (PLE(l-1)-PLE(L))
        enddo
        enddo
      ENDIF

c
c rate of change of temperature and zonal wind
c
      if(IDACC(ia_dga).GT.1) then
        do l=1,lm
          n = JK_TOTDTDT
          agc(:,l,n) = agc(:,l,n)*pkm(l)*idacc(ia_dga)/
     &         ((IDACC(ia_dga)-1)*(DTsrc*NDAA+DT+DT))
          n = JK_TOTDUDT
          agc(:,l,n) = agc(:,l,n)*idacc(ia_dga)/
     &         ((IDACC(ia_dga)-1)*(DTsrc*NDAA+DT+DT))
        enddo
      endif

c
c flux divergences/convergences
c
      DO L=1,LM
C     agc(1,l,jk_dyn_conv_dse)=-AGC(2,l,JK_TOTNTDSE)
        agc(1,l,jk_dyn_conv_dse)=0.
        agc(jm,l,jk_dyn_conv_dse)=AGC(JM,l,JK_TOTNTDSE)*BYDAPO(JM)
C     agc(1,l,jk_dyn_conv_eddy_geop)=-(AGC(2,l,JK_EDDNTGEO))
        agc(1,l,jk_dyn_conv_eddy_geop)=0.
        agc(jm,l,jk_dyn_conv_eddy_geop)=AGC(JM,l,JK_EDDNTGEO)*BYDAPO(JM)
        DO J=2,JM-1
          agc(j,l,jk_dyn_conv_dse)=
     &         (AGC(J,l,JK_TOTNTDSE)-AGC(J+1,l,JK_TOTNTDSE))*BYDXYP(J)
          agc(j,l,jk_dyn_conv_eddy_geop)=
     &         (AGC(J,l,JK_EDDNTGEO)-AGC(J+1,l,JK_EDDNTGEO))*BYDXYP(J)
        enddo
      enddo
      DO L=1,LM-1
      DO J=1,JM
        agc(j,l,jk_dyn_conv_dse)=
     &       agc(j,l,jk_dyn_conv_dse)    +AGC(J,l,JK_TOTVTDSE)
        agc(j,l+1,jk_dyn_conv_dse)=
     &       agc(j,l+1,jk_dyn_conv_dse)  -AGC(J,l,JK_TOTVTDSE)
        agc(j,l,jk_dyn_conv_eddy_geop)=
     &       agc(j,l,jk_dyn_conv_eddy_geop)  +AGC(J,l,JK_VTGEOEDDY)
        agc(j,l+1,jk_dyn_conv_eddy_geop)=
     &       agc(j,l+1,jk_dyn_conv_eddy_geop)-AGC(J,l,JK_VTGEOEDDY)
      ENDDO
      ENDDO

      n = jk_dtempdt_econv
      agc(:,:,n) = 0.
      do l=2,lm-1
      do j=2,jm-1
        if (AGC(J,L,JK_DPA).le.0.) cycle
        AGC(J,L,N)=
     &       (AGC(J,L,JK_EDDNTSH)-AGC(J+1,L,JK_EDDNTSH))/DXYP(J)
     &      +(AGC(J,L,JK_EDDVTPT)-AGC(J,L-1,JK_EDDVTPT))*PKM(L)*IM
      enddo
      enddo

      n = jl_epflx_div
      agc(:,:,n) = 0.
      DO J=2,JM-1
        DO L=1,LM
          AGC(J,L,N) =
     &          agc(j,l,jl_dpsig)*bydapo(j)*
     &         (AGC(J+1,L,JL_EPFLXN)*DXV(J+1)/AGC(J+1,L,JL_DPSIGV)
     &         -AGC(J  ,L,JL_EPFLXN)*DXV(J  )/AGC(J  ,L,JL_DPSIGV))
        ENDDO
        agc(j, 1,n) = agc(j,1,n) +.5*byim*agc(j,1,jl_epflxv)/dsig(1)
        do l=2,lm
          agc(j,l,n) = agc(j,l,n) +.5*byim*
     &         (agc(j,l,jl_epflxv)-agc(j,l-1,jl_epflxv))/dsig(l)
        enddo
      ENDDO

C****
C**** Calculate a density field on tracer grid, edge pressure
C****     (not quite ok if K-1 is underground?)
c NOTE: THIS IS AT LAYER BOTTOM, NOT LAYER TOP.

      kk = 0
      do k=1,kagcx
        if(denom_gc(k)==jk_rhoedz) then
          kk = kk + 1
          inds(kk) = k
        endif
      enddo

      temp(:,:) = tf + agc(:,:,jk_temp)/(agc(:,:,jk_dpa)+teeny)
      do l=2,lm
        do j=1,jm
          wtup = (plm(l-1)-pme(l))/(plm(l-1)-plm(l))
          tedge = (1.-wtup)*temp(j,l-1) + wtup*temp(j,l)
          rhoe = 100.*pme(l)/(rgas*tedge)
          delpe = wtup*(pme(l)-ple(l))
     &     + (1.-wtup)*(pme(l-1)-pme(l))
          agc(j,l,jk_rhoedz) = idacc(ia_dga)*delpe
          agc(j,l,inds(1:kk)) = agc(j,l,inds(1:kk))*(delpe/rhoe)
        enddo
      enddo
      l = 1
      do j=1,jm
        rhoe = 100.*pme(l)/(rgas*temp(j,l))
        delpe = .5d0*(pme(l)-pme(l+1))
        agc(j,l,jk_rhoedz) = idacc(ia_dga)*delpe
        agc(j,l,inds(1:kk)) = agc(j,l,inds(1:kk))*(delpe/rhoe)
      enddo

c
c copy pressure thickness to conditionally averaged output instance
c
      agc(:,:,jk_dpb2) = agc(:,:,jk_dpb)

c
c compute hemispheric/global means and vertical sums
c
      wtj(:,1) = 2.*dxyp(:)/sum(dxyp)
      wtj(2:jm,2) = 2.*dxyv(2:jm)/sum(dxyp)
      wtj(jm/2+1,2) = wtj(jm/2+1,2)*.5
      do k=1,kagc
        jg = jgrid_gc(k)
        do l=1,lm
          j1 = jg; j2 = jm/2+jg-1
          hemis_gc(1,l,k) = sum(agc(j1:j2,l,k)*wtj(j1:j2,jg))
          j1 = jm/2+1; j2 = jm
          hemis_gc(2,l,k) = sum(agc(j1:j2,l,k)*wtj(j1:j2,jg))
          hemis_gc(3,l,k) = .5*(hemis_gc(1,l,k)+hemis_gc(2,l,k))
        enddo
        if(force_gc_vmean(k)) then
          ! Note dsig is the wrong weight for layer-edge qtys.
          ! It is used for the time being to match the behavior of JLMAP.
          do j=jg,jm
            vmean_gc(j,1,k) = sum(dsig(:)*agc(j,:,k))/sum(dsig)
          enddo
          do j=jm+1,jm+3
            vmean_gc(j,1,k) = sum(dsig(:)*hemis_gc(j-jm,:,k))/sum(dsig)
          enddo
        else
          vmean_gc(jg:jm,1,k) = sum(agc(jg:jm,:,k),dim=2)
          vmean_gc(jm+1:jm+3,1,k) = sum(hemis_gc(:,:,k),dim=2)
        endif
      enddo
      return
      end subroutine diaggc_prep

      subroutine calc_derived_aijk
!@sum Calculate derived AIJK diagnostics prior to printing
!@auth Group
      use constant, only : sha
      use resolution, only : im,lm
      use diag_com, only : aijk => aijk_loc
      use gcdiag, only : ijk_dse,ijk_tb,ijk_phi
      use domain_decomp_atm, only : grid
      implicit none
      integer :: i,j,l,j_0stg,j_1stg

      j_0stg = grid%j_strt_stgr
      j_1stg = grid%j_stop_stgr

      do l=1,lm
        do j=j_0stg,j_1stg
          do i=1,im
            aijk(i,j,l,ijk_phi) =
     &           aijk(i,j,l,ijk_dse) - sha*aijk(i,j,l,ijk_tb)
          enddo
        enddo
      enddo

      return
      end subroutine calc_derived_aijk
