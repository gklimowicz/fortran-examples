#include "rundeck_opts.h"
      subroutine def_acc
c-----------------------------------------------------------------------
c define acc names, units, etc
c-----------------------------------------------------------------------
      use DIAG_COM
      implicit none
      integer :: k
      call uvgrid_defs
      call tsf_defs
      call j_defs
      name_reg=name_j
      do k=1,kaj ! to avoid naming conflicts, put a prefix
         name_j(k) = 'J_'//trim(name_j(k))
         name_reg(k) = 'reg_'//trim(name_reg(k))
      enddo
      call jl_defs
      call sjl_defs
      call ij_defs
c      call il_defs
#ifndef SCM
      call wave_defs
      call gc_defs
#endif
      call ijl_defs
#ifndef SCM
      call ijk_defs
#endif
      call diurn_defs
      call ijhc_defs

#ifndef STANDALONE_OCEAN
#if (defined TRACERS_ON) || (defined TRACERS_OCEAN)
      call InitTracerDiagMetadata()
#endif
#endif
      return
      end subroutine def_acc

      subroutine uvgrid_defs
!@sum uvgrid_defs defines the latitudinal placement of U-wind and
!@+   KE diagnostics
!@auth M. Kelley
      use DIAG_COM
c      use GEOM, only : jg_u, jg_ke
#ifdef CUBED_SPHERE
      jgrid_u = 1
      jgrid_ke = 1
#else
      jgrid_u = 2 !jg_u
      jgrid_ke = 2 !jg_ke
#endif
c#ifdef AGRID_DIAG
c      jgrid_u = 1
c      jgrid_ke = 1
c#endif
      return
      end subroutine uvgrid_defs

      subroutine j_defs
!@sum j_defs definitions for j_xx zonal budget diagnostics
!@+   diags are printed out in the order they are defined
!@auth G. Schmidt/M. Kelley
      use CONSTANT, only : grav,shw,rgas,omega,bygrav,gamd
     &     ,radian,radius
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : dtsrc,kocean,qcheck
      use DIAG_COM
      use DIAG_COM_RAD
      use DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use MDIAG_COM, only : make_timeaxis
      use cdl_mod
      USE FLUXES, only : atmice
      implicit none
      character(len=30), parameter ::
     &     fmt906='(A16,3F7.2,2X,24F4.1)'
     &    ,fmt907='(A16,3F7.2,2X,24I4)'
     &    ,fmt911='(A16,3F7.3,2X,24F4.1)'
     &    ,fmt912='(A16,3F7.3,2X,24I4)'
     &    ,fmt909='(A16,1X,23I5)'
     &    ,fmt910='(A16,1X,23F5.1)'
     &    ,fmtnone='not computed'
      integer :: k,kk
      character(len=30) :: sname
      character(len=8) :: namreg_1word(23)
      logical :: set_miss
c
      do k=1,kaj
         write(name_j(k),'(a2,i3.3)') 'AJ',k
         lname_j(k) = 'unused'
         units_j(k) = 'unused'
         stitle_j(k)= 'no output'
         scale_j(k) = 1.
         ia_j(k)    = 1.
         fmt_j(k) = fmt907
         fmt_reg(k) = fmt909
         iden_j(k)  = 0
      enddo
c
      k=0
c
      k=k+1
      J_SRINCP0= k ! SRINCP0 (W/M**2)                              2 RD
      name_j(k) = 'inc_sw'
      lname_j(k) = 'SOLAR RADIATION INCIDENT ON PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' INC SW (W/m^2)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k = k + 1
      J_ALBP0  = k ! J_SRINCP0-J_SRNFP0 will be calculated after acc
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'plan_alb'
      lname_j(k) = ' TOTAL PLANETARY ALBEDO'
      units_j(k) = '%'
      stitle_j(k)= ' PLANETARY ALBDO'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_PLAVIS = k ! PLAVIS*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'plan_alb_vis'
      lname_j(k) = 'PLANETARY ALBEDO IN VISUAL'
      units_j(k) = '%'
      stitle_j(k)= ' PLAN ALB VISUAL'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_PLANIR = k ! PLANIR*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'plan_alb_nir'
      lname_j(k) = 'PLANETARY ALBEDO IN NEAR IR'
      units_j(k) = '%'
      stitle_j(k)= ' PLAN ALB NEARIR'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_ALBG  = k ! J_SRINCG-J_SRNFG will be calculated after acc
      iden_j(k)  = J_SRINCG ! J_SRINCG not defined yet.  fixup below.
      name_j(k) = 'surf_alb'
      lname_j(k) = 'GROUND ALBEDO'
      units_j(k) = '%'
      stitle_j(k)= ' SURFACE G ALBDO'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_ALBVIS = k ! ALBVIS*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'surf_alb_vis'
      lname_j(k) = 'GROUND ALBEDO IN VISUAL'
      units_j(k) = '%'
      stitle_j(k)= ' SURF ALB VISUAL'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_ALBNIR = k ! ALBNIR*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'surf_alb_nir'
      lname_j(k) = 'GROUND ALBEDO IN NEAR IR'
      units_j(k) = '%'
      stitle_j(k)= ' SURF ALB NEARIR'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_SRRVIS = k ! SRRVIS*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'atm_alb_vis'
      lname_j(k) = 'ATMOSPHERIC ALBEDO IN VISUAL'
      units_j(k) = '%'
      stitle_j(k)= '0ATMO ALB VISUAL'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_SRRNIR = k ! SRRNIR*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'atm_alb_nir'
      lname_j(k) = 'ATMOSPHERIC ALBEDO IN NEAR IR'
      units_j(k) = '%'
      stitle_j(k)= ' ATMO ALB NEARIR'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_SRAVIS = k ! SRAVIS*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'atm_abs_vis'
      lname_j(k) = 'ATMOSPHERIC ABSORPTION IN VISUAL'
      units_j(k) = 'W/m**2'
      stitle_j(k)= ' ATMO ABS VISUAL'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k = k + 1
      J_SRANIR = k ! SRANIR*S0*COSZ (W/m**2)                       2 RD
      iden_j(k)  = J_SRINCP0
      name_j(k) = 'atm_abs_nir'
      lname_j(k) = 'ATMOSPHERIC ABSORPTION IN NEAR IR'
      units_j(k) = 'W/m**2'
      stitle_j(k)= ' ATMO ABS NEARIR'
      scale_j(k) = 100.
      fmt_j(k) = fmt912
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRNFP0=  k ! SRNFP0 (W/m**2)                               2 RD
      name_j(k) = 'sw_abs_p0'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0SW ABS BELOW P0'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRNFP1=  k ! SRNFP1 (W/m**2)                               2 RD
      name_j(k) = 'sw_abs_p1'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BELOW PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS BELOW P1'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRABS =  k ! SRABSATM=AJ(SRNFP0)-AJ(SRNFG) (W/m**2)        2 D1
      name_j(k) = 'sw_abs_atm'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY ATMOSPHERE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS BY ATMOS'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRINCG=  k ! SRINCG (W/m**2)                               2 RD
      name_j(k) = 'sw_inc_z0'
      lname_j(k) = 'SOLAR RADIATION INCIDENT ON GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW INC ON Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_SRNFG =  k ! SRNFG (W/m**2)                                2 RD
      name_j(k) = 'sw_abs_z0'
      lname_j(k) = 'SOLAR RADIATION ABSORBED BY GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' SW ABS AT Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TRNFP0=  k ! TRNFP0 (W/m**2)2 D1
      name_j(k) = 'net_lw_p0'
      lname_j(k) = 'THERMAL RADIATION EMITTED BY PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0NET LW AT P0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TRNFP1=  k ! TRNFP1=AJ(ALBNIR) D1
      name_j(k) = 'net_lw_p1'
      lname_j(k) = 'NET THERMAL RADIATION AT PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET LW AT P1   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TRHDT =  k ! TRHDT (J/m**2)                                1 SF
      name_j(k) = 'net_lw_z0'
      lname_j(k) = 'NET THERMAL RADIATION AT GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET LW AT Z0   '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_TRINCG= k ! TRINCG (W/m**2)                               2 RD
      name_j(k) = 'lw_inc_z0'
      lname_j(k) = 'THERMAL RADIATION INCIDENT ON GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' LW INC ON Z0   '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_LWCORR= k ! LWUP(RAD)-LW(SURF)   (W/m2)
      name_j(k) = 'LWCORR'
      lname_j(k) = 'LW RAD/SURF CORRECTION'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' LW TO L1'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_BRTEMP= k ! BTEMPW-TF                                     2 RD
      name_j(k) = 'btemp_window'
      lname_j(k) = 'BRIGHTNESS TEMP THROUGH WINDOW REGION'
      units_j(k) = 'C'
      stitle_j(k)= ' LW WINDOW BTEMP'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_RNFP0 = k ! RNFP0=AJ(SRNFP0)+AJ(TRNFP0) (W/m**2)          2 D1
      name_j(k) = 'net_rad_p0'
      lname_j(k) = 'NET RADIATION OF PLANET'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0NET RAD AT P0  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLRTOA= k ! clear sky radiative forcing top of atmosphere
      name_j(k) = 'net_clr_toa'
      lname_j(k) = 'NET CLEAR SKY RADIATION AT P0'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET CLR RAD  P0'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_RNFP1 = k ! RNFP1=AJ(SRNFP1)+AJ(TRNFP1) (W/m**2)          2 D1
      name_j(k) = 'net_rad_p1'
      lname_j(k) = 'NET RADIATION BELOW PTOP'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET RAD AT P1  '
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLRTRP= k ! clear sky radiative forcing tropopause
      name_j(k) = 'net_clr_trp'
      lname_j(k) = 'NET CLEAR-SKY RADIATION AT TROPOPAUSE (WMO)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET CLR RAD TRP'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_TOTTRP= k ! total radiative forcing tropopause
      name_j(k) = 'net_tot_trp'
      lname_j(k) = 'NET RADIATION AT TROPOPAUSE (WMO)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET RAD (TROPP)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_RHDT  = k ! RHDT=A1BYA2*AJ(SRNFG)*DTS+AJ(TRHDT)(J/m^2)    1 D1
      name_j(k) = 'net_rad_z0'
      lname_j(k) = 'NET RADIATION ABSORBED BY GROUND'
      units_j(k) = 'J/m^2'
      stitle_j(k)= ' NET RAD AT Z0  '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_SHDT  = k ! SHEATDT (J/m**2)                              1 SF
      name_j(k) = 'snsht_flx'
      lname_j(k) = 'SENSIBLE HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0SENSBL HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_EVHDT = k ! EVHDT (J/m**2)                                1 SF
      name_j(k) = 'evht_flx'
      lname_j(k) = 'LATENT HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' EVAPOR HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_EPRCP = k ! ENERGP (J/m**2)                               1 CN
      name_j(k) = 'prec_ht_flx'
      lname_j(k) = 'PRECIPITATION HEAT FLUX INTO THE GROUND'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' PRECIP HEAT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt912
c
      k=k+1
      J_HZ0   = k ! HEATZ0=RHDT+SHDT+EVHDT+EPRCP (J/m**2)   1 D1
      name_j(k) = 'nt_ht_z0'
      lname_j(k) = 'NET HEATING AT GROUND SURFACE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT AT Z0 '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_ERVR = k ! Energy in river discharge                     1 GP
      name_j(k) = 'ht_rvr_disch'
      lname_j(k) = 'HEAT IN RIVER DISCHARGE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0HT RVR DISCH   '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_HZ1   = k ! Net heating at ocean surface                  1 D1
      name_j(k) = 'net_ht_z1'
      lname_j(k) = 'NET HEAT AT SURFACE (INCL RIVERS)'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HEAT Z0+RVR'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_ERUN = k ! ERUNOFF (J/m**2)                                1 GP
      name_j(k) = 'ht_runoff'
      lname_j(k) = 'HEAT RUNOFF'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0HEAT RUNOFF '
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      atmice%J_HMELT = k ! net amount of energy associated with ice melt/form
      name_j(k) = 'ht_ice_melt'
      lname_j(k) = 'NET HEAT OF ICE MELT/FORMATION'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HT ICE MLT/FORM'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt912
      fmt_reg(k) = fmt910
c
C**** Note this is used for ice in fixed SST runs, but for ocean in
C**** qflux runs. Over land, it is always used for landice changes.
      k=k+1
      J_IMPLH = k               !                                 1 GP
      atmice%J_IMPLH = k
      name_j(k) = 'impl_ht'
      lname_j(k) = 'DOWNWARD IMPLICIT HEAT FLUX AT ICE BASE/OCN ML'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' DWN IMPL HT FLX'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
      fmt_reg(k) = fmt910
c
      k=k+1
      J_HZ2   = k !
      name_j(k) = 'net_ht_hz2'
      lname_j(k) = 'NET HEAT CONVERGENCE INTO GROUND TYPE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' NET HT CNV GRND'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
c
      k=k+1
      J_PRCPSS= k ! PRCPSS (100 PA)                               1 CN
      name_j(k) = 'ssprec'
      lname_j(k) = 'SUPER SATURATION PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= '1SS PRECIP(MM/D)'
      scale_j(k) = 100.*SECONDS_PER_DAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_PRCPMC= k ! PRCPMC (100 PA)                               1 CN
      name_j(k) = 'mcprec'
      lname_j(k) = 'MOIST CONVECTIVE PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' MC PRECIP(MM/D)'
      scale_j(k) = 100.*SECONDS_PER_DAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_PRCP  = k ! PRCP=AJ(PRCPSS)+AJ(PRCPMC) (100 PA)           1 D1
      name_j(k) = 'prec'
      lname_j(k) = 'PRECIPITATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' PRECIP (MM/DAY)'
      scale_j(k) = 100.*SECONDS_PER_DAY/(DTsrc*GRAV)
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_EVAP  = k ! EVAP (KG/m**2)                                1 GD
      name_j(k) = 'evap'
      lname_j(k) = 'EVAPORATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' EVAPOR (MM/DAY)'
      scale_j(k) = SECONDS_PER_DAY/DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_IRGW = k               !                                 1 GP
      name_j(k) = 'irrig_external'
      lname_j(k) = 'IRRIGATION WATER FROM EXTERNAL SOURCE (GRNDWATER)'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' IRR ADD(MM/DAY)'
      scale_j(k) = 1000.d0 * SECONDS_PER_DAY
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_IRGWE = k               !                                 1 GP
      name_j(k) = 'irrig_externalE'
      lname_j(k) = 'IRRIGATION HEAT FROM EXTERNAL SOURCE (GRNDWATER)'
      units_j(k) = 'W/m2'
      stitle_j(k)= ' IRRE AD(W/m2)'
      scale_j(k) = 1.d0
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_RUN  = k ! RUNOFF (KG/m**2)                                1 GP
      name_j(k) = 'wat_runoff'
      lname_j(k) = 'WATER RUNOFF AT GROUND SURFACE'
      units_j(k) = 'mm/day'
      stitle_j(k)= '0WATER RUNOFF'
      scale_j(k) = SECONDS_PER_DAY/DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_RVRD  = k ! RIVER DISCHARGE                           1 GP
      name_j(k) = 'river_discharge'
      lname_j(k) = 'RIVER DISCHARGE'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' RVR DISCH(MM/D)'
      scale_j(k) = SECONDS_PER_DAY/DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      atmice%J_IMELT = k ! net amount ice melt/formation in oc/lk  1 GP
      name_j(k) = 'ice_melt'
      lname_j(k) = 'NET ICE MELTING/FORMATION'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' ICE MELT/FORM  '
      scale_j(k) = SECONDS_PER_DAY/DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
C**** Note this is used for ice in fixed SST runs, but for ocean in
C**** qflux runs. Over land, it is always used for landice changes.
      k=k+1
      J_IMPLM = k               !                                 1 GP
      atmice%J_IMPLM = k
      name_j(k) = 'impl_m_flux'
      lname_j(k) =
     *     'DOWNWARD IMPLICIT FRESHWATER FLUX AT ICE BASE/OCN ML'
      units_j(k) = 'mm/day'
      stitle_j(k)= ' DWN IMPL WT FLX'
      scale_j(k) = SECONDS_PER_DAY/DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_H2OCH4 = k               !                                 1 GP
      name_j(k) = 'h2o_from_ch4'
      lname_j(k) = 'WATER DERIVED FROM CH4 OXIDATION IN STRATOSPHERE'
      units_j(k) = '10^-6 mm/day'
      stitle_j(k)= ' H2O BY CH4(x1M)'
      scale_j(k) = 2d6
      ia_j(k) = ia_12hr
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      atmice%J_SMELT= k ! salt flux associated with ice melt/formation   1 GD
      name_j(k) = 's_ice_melt'
      lname_j(k) = 'SALT IN ICE MELT/FORMATION'
      units_j(k) = '10^-3 kg/m^2/day'
      stitle_j(k)= '0SALT MELT (x1K)'
      scale_j(k) = 1000.*SECONDS_PER_DAY/DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
c
      k=k+1
      J_TG1   = k ! TG1 (K-TF)                                    1 GD
      name_j(k) = 'tg1'
      lname_j(k) = 'TEMPERATURE OF GROUND LAYER 1'
      units_j(k) = '.1 C'
      stitle_j(k)= '0TG1 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c
      k=k+1
      J_TG2   = k ! TG2 (K-TF)                                    1 GD
      name_j(k) = 'tg2'
      lname_j(k) = 'TEMPERATURE OF GROUND LAYER 2'
      units_j(k) = '.1 C'
      stitle_j(k)= ' TG2 (.1 C)     '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c
      k=k+1
      J_TSRF  = k ! TS (K-TF)                                     3 SF
      name_j(k) = 'tsurf'
      lname_j(k) = 'SURFACE AIR TEMPERATURE'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T SURF (.1 C)  '
      scale_j(k) = 10.
      ia_j(k) = ia_srf
c
      k=k+1
      J_TX    = k ! TX (K-TF)  (INTEGRAL OVER ATMOSPHERE OF)      4 DA
      name_j(k) = 'tair'
      lname_j(k) = 'AIR TEMPERATURE'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T AIR (.1 C)   '
      scale_j(k) = 10.
      ia_j(k) = ia_dga
c
      k=k+1
      J_TX1   = k ! TX1 (K-TF)                                    4 DA
      name_j(k) = 't1'
      lname_j(k) = 'TEMPERATURE OF AIR LAYER 1'
      units_j(k) = '.1 C'
      stitle_j(k)= ' T1 (.1 C)      '
      scale_j(k) = 10.
      ia_j(k) = ia_dga
c
      k=k+1
      J_DTDJS = k ! T(J+1)-T(J-1)  (SUM OVER STRATOSPHERE OF)     4 DA
      name_j(k) = 'dtdlat_strat'
      lname_j(k) = 'STRATO TEMP CHANGE PER DEGREE LATITUDE'
      units_j(k) = 'deg C/deg lat'
      stitle_j(k)= '0DT/DLAT(STRAT) '
      scale_j(k) = 100.*radius*radian
      ia_j(k) = ia_dga
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_DTDJT = k ! T(J+1)-T(J-1)  (SUM OVER TROPOSPHERE OF)      4 DA
      name_j(k) = 'dtdlat_trop'
      lname_j(k) = 'TROPO TEMP CHANGE PER DEGREE LATITUDE'
      units_j(k) = 'deg C/deg lat'
      stitle_j(k)= ' DT/DLAT(TROPO) '
      scale_j(k) = 100.*radius*radian
      ia_j(k) = ia_dga
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_DTSGST= k ! DTH/DPHI  (STRATOSPHERE)                      4 DA
      name_j(k) = 'sstab_strat'
      lname_j(k) = 'STRATOSPHERIC STATIC STABILITY'
      units_j(k) = 'C/km'
      stitle_j(k)= '0STAT STB(STRAT)'
      scale_j(k) = 1.D3*GRAV*P1000K
      ia_j(k) = ia_dga
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_DTDGTR= k ! DTH/DPHI  (TROPOSPHERE)                       4 DA
      name_j(k) = 'sstab_trop'
      lname_j(k) = 'TROPOSPHERIC STATIC STABILITY'
      units_j(k) = 'C/km'
      stitle_j(k)= ' STAT STB(TROPO)'
      scale_j(k) = 1.D3*GRAV*P1000K
      ia_j(k) = ia_dga
      fmt_j(k) = fmt906
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_RICST = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (STRATOSPHERE)  4 DA
      name_j(k) = 'rich_num_strat'
      lname_j(k) = 'STRATOSPHERIC RICHARDSON NUMBER'
      units_j(k) = '1'
      stitle_j(k)= '0RICH NUM(STRAT)'
      scale_j(k) = 1.
      ia_j(k) = ia_dga
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_RICTR = k ! .0625*DTH*DLNP/(DU*DU+DV*DV)  (TROPOSPHERE)   4 DA
      name_j(k) = 'rich_num_trop'
      lname_j(k) = 'TROPOSPHERIC RICHARDSON NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' RICH NUM(TROPO)'
      scale_j(k) = 1.
      ia_j(k) = ia_dga
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_ROSST = k ! UMAX/(FCOR)  (STRATOSPHERE)              4 DA
      name_j(k) = 'ross_num_strat'
      lname_j(k) = 'STRATOSPHERIC ROSSBY NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' ROSS NUM(STRAT)'
      scale_j(k) = 1./(L_ROSSBY_NUMBER)
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_ROSTR = k ! UMAX/(FCOR)  (TROPOSPHERE)               4 DA
      name_j(k) = 'ross_num_trop'
      lname_j(k) = 'TROPOSPHERIC ROSSBY NUMBER'
      units_j(k) = '1'
      stitle_j(k)= ' ROSS NUM(TROPO)'
      scale_j(k) = 1./(L_ROSSBY_NUMBER)
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_LSTR  = k ! SQRT(DPHI*DLOGTH)/SINJ  (STRATOSPHERE)           4 DA
      name_j(k) = 'ross_radius_strat'
      lname_j(k) = 'ROSSBY RADIUS IN THE STRATOSPHERE'
      units_j(k) = '10**6 m'
      stitle_j(k)= ' L(STRAT)(10**6)'
      scale_j(k) = 1d-6!/(2.*OMEGA)
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_LTRO  = k ! SQRT(DPHI*DLOGTH)/SINJ  (TROPOSPHERE)            4 DA
      name_j(k) = 'ross_radius_trop'
      lname_j(k) = 'ROSSBY RADIUS IN THE TROPOSPHERE'
      units_j(k) = '10**6 m'
      stitle_j(k)=  ' L(TROP) (10**6)'
      scale_j(k) = 1d-6!/(2.*OMEGA)
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_GAM   = k ! GAM  (K/m)  (*SIG(TROPOSPHERE)/GRAV)          4 DA
      name_j(k) = 'lapse_rate'
      lname_j(k) = 'MEAN LAPSE RATE'
      units_j(k) = 'K/km'
      stitle_j(k)= '0GAM(K/KM)      '
      scale_j(k) = 1d3*GRAV
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_GAMM  = k ! GAMM  (K-S**2/m**2)  (SIG(TROPOSPHERE)/GAMD)  4 DA
      name_j(k) = 'lapse_rate_m'
      lname_j(k) = 'MOIST ADIABATIC LAPSE RATE'
      units_j(k) = 'K/km'
      stitle_j(k)= ' GAMM(K/KM)     '
      scale_j(k) = 1d3*GAMD
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_GAMC  = k ! GAMC  (K/m)                                   4 DA
      name_j(k) = 'lapse_rate_c'
      lname_j(k) = 'GAMC'
      units_j(k) = 'K/Km'
      stitle_j(k)= ' GAMC(K/KM)     '
      scale_j(k) = 1d3
      ia_j(k) = ia_dga
      fmt_j(k) = fmt911
      fmt_reg(k) = fmtnone
c
      k=k+1
      J_PCLDSS= k ! PCLDSS (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'sscld'
      lname_j(k) = 'SUPER SATURATION CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= '0TOT SUP SAT CLD'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c
      k=k+1
      J_PCLDMC= k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE)       2 RD
      name_j(k) = 'mccld'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= ' TOT MST CNV CLD'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c
      k=k+1
      J_PCLD  = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)         2 RD
      name_j(k) = 'totcld'
      lname_j(k) = 'TOTAL CLOUD COVER'
      units_j(k) = '%'
      stitle_j(k)= ' TOTAL CLD COVER'
      scale_j(k) = 100.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CLDDEP= k ! PBOTMC-PTOPMC (100 PA)                        2 RD
      iden_j(k)  = J_PCLDMC
      name_j(k) = 'mc_clddp'
      lname_j(k) = 'MOIST CONVECTIVE CLOUD DEPTH'
      units_j(k) = 'mb'
      stitle_j(k)= ' MC CLD DPTH(MB)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
      fmt_j(k) = fmt907
c
      k=k+1
      J_QP    = k ! Q*P (100 PA)  (INTEGRAL OVER ATMOSPHERE OF)   4 DA
      name_j(k) = 'atmh2o'
      lname_j(k) = 'WATER CONTENT OF ATMOSPHERE'
      units_j(k) = 'mm'
      stitle_j(k)= ' H2O OF ATM (MM)'
      scale_j(k) = 100.*BYGRAV
      ia_j(k) = ia_dga
      fmt_reg(k) = fmt910
c
      k=k+1
      J_WTR1  = k ! WTR1 (KG/m**2)                                1 GD
      name_j(k) = 'wat_g1'
      lname_j(k) = 'WATER IN GROUND LAYER 1'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= '0WATER IN G1    '
      scale_j(k) = 1.
      ia_j(k) = ia_src
      fmt_reg(k) = fmt910
c
      k=k+1
      J_ACE1  = k ! ACE1 (KG/m**2)                                1 GD
      atmice%J_ACE1 = k
      name_j(k) = 'ice_g1'
      lname_j(k) = 'ICE IN GROUND LAYER 1'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' ICE IN G1      '
      scale_j(k) = 1.
      ia_j(k) = ia_src
      fmt_reg(k) = fmt910
c
      k=k+1
      J_WTR2  = k ! WTR2 (KG/m**2)                                1 GD
      name_j(k) = 'wat_g2'
      lname_j(k) = 'WATER IN GROUND LAYER 2'
      units_j(k) = '10^2 kg/m^2'
      stitle_j(k)= ' WATER G2 x0.01 '
      scale_j(k) = 1d-2
      ia_j(k) = ia_src
c
      k=k+1
      J_ACE2  = k ! ACE2 (KG/m**2)                                1 GD
      atmice%J_ACE2 = k
      name_j(k) = 'ice_g2'
      lname_j(k) = 'ICE IN GROUND LAYER 2'
      units_j(k) = '10^2 kg/m^2'
      stitle_j(k)= ' ICE G2   x0.01 '
      scale_j(k) = 1d-2
      ia_j(k) = ia_src
c
      k=k+1
      J_SNOW  = k ! SNOW (KG/m**2)                                1 GD
      atmice%J_SNOW = k
      name_j(k) = 'snowdp'
      lname_j(k) = 'SNOW DEPTH'
      units_j(k) = 'kg/m^2'
      stitle_j(k)= ' SNOW DEPTH     '
      scale_j(k) = 1.
      ia_j(k) = ia_src
      fmt_reg(k) = fmt910
c
      k=k+1
      J_RSNOW = k ! PSNOW (1)                                     4 DA
      atmice%J_RSNOW = k
      name_j(k) = 'snow_cover'
      lname_j(k) = 'SNOW COVER'
      units_j(k) = '%'
      stitle_j(k)= ' SNOW COVER     '
      scale_j(k) = 100.
      ia_j(k) = ia_src
      fmt_reg(k) = fmt910
c
      k=k+1
      J_RSI   = k ! RSI (1)                                       1 GD
      atmice%J_RSI = k
      name_j(k) = 'ocn_lak_ice_frac'
      lname_j(k) = 'OCEAN/LAKE ICE COVER'
      units_j(k) = '%'
      stitle_j(k)= ' OC/LK ICE COVER'
      scale_j(k) = 100.
      ia_j(k) = ia_src
      fmt_reg(k) = fmt910
c
      IF (KOCEAN.gt.0) THEN ! only for non-fixed SST runs
      k=k+1
      J_OHT   = k ! OCEAN TRANSPORT                               1 GD
      name_j(k) = 'ocn_ht_trans'
      lname_j(k) = 'CONVERGED OCEAN HEAT TRANSPORT'
      units_j(k) = 'W/m^2'
      stitle_j(k)= '0OCEAN TRNS CONV'
      scale_j(k) = 1./DTSRC
      ia_j(k) = ia_src
      fmt_j(k) = fmt912
c
      k=k+1
      J_FTHERM= k ! ENERGY DIFFUSION INTO THERMOCLINE (W/m**2) dly odeep
      name_j(k) = 'ht_thermocline'
      lname_j(k) = 'ENERGY DIFFUSION INTO THE THERMOCLINE'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' HT INTO THRMOCL'
      scale_j(k) = 2d3*SHW/SECONDS_PER_DAY
      ia_j(k) = ia_12hr
      fmt_j(k) = fmt911
      fmt_reg(k) = fmt910
      END IF
c
      k=k+1
      J_TYPE  = k ! PTYPE                                         1 GD
      name_j(k) = 'surf_type_frac'
      lname_j(k) = 'SURF TYPE FRACT'
      units_j(k) = '%'
      stitle_j(k)= ' SURF TYPE FRACT'
      scale_j(k) = 100.
      ia_j(k) = ia_srf
      fmt_reg(k) = fmtnone

c None of the following will be printed out:
      k=k+1
      J_HSURF = k ! TRNFP0-TRNFG (W/m**2)                         2 RD
      name_j(k) = 'HSURF'
      lname_j(k) = 'SURFACE THERMAL HEATING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_HATM  = k ! TRNFP1-TRNFG (W/m**2)                         2 RD
      name_j(k) = 'HATM'
      lname_j(k) = 'ATMOSPHERIC THERMAL HEATING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= 'no output'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
#ifdef HEALY_LM_DIAGS
c
      k=k+1
      J_VTAU  = k ! Volcanic TAU*-20 (W/m^2)
      name_j(k) = 'vol_tau_forc'
      lname_j(k) = 'VOLCANIC OPTICAL DEPTH FORCING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' VOLC (.1 W/m^2)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_GHG  = k ! Greenhouse forcing (W/m^2)
      name_j(k) = 'ghg_forc'
      lname_j(k) = 'GREENHOUSE GAS EFFECTIVE FORCING'
      units_j(k) = 'W/m^2'
      stitle_j(k)= ' GHG Fe(.1 W/m^2)'
      scale_j(k) = 1.
      ia_j(k) = ia_rad
c
      k=k+1
      J_CROPS  = k ! Crop cover (%)
      name_j(k) = 'crop_cover'
      lname_j(k) = 'CROP COVER'
      units_j(k) = '%'
      stitle_j(k)= ' CROP COVER '
      scale_j(k) = 1.
      ia_j(k) = ia_src
#endif

c fixups for cases when denominator indices were not yet defined.
c change specifification of denominators to use name_j instead.
      iden_j(J_ALBG)  = J_SRINCG

c
      if (k .gt. kaj) then
        if(am_i_root()) then
          write (6,*) 'j_defs: Increase kaj=',kaj,' to at least ',k
        endif
        call stop_model( 'kaj too small', 255 )
      endif
      if(am_i_root()) then
        write (6,*) 'Number of AJ diagnostics defined: kajmax=',k
        if(qcheck) then
          do kk=1,k
            write (6,'(i4,'':'',a)') kk,trim(lname_j(kk))
          end do
        endif
      endif

      iden_reg(:) = iden_j(:)

c set denominator to surface type fraction for qtys that are
c not ratios
      do k=1,kaj
        if(k.eq.J_TYPE) cycle
        if(iden_j(k).gt.0) cycle
        iden_j(k)  = J_TYPE
      enddo

c
c Declare the dimensions and metadata of AJ output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      call init_cdl_type('cdl_latbudg',cdl_latbudg)
      call add_coord(cdl_latbudg,'lat_budg',jm_budg,
     &     units='degrees_north',coordvalues=lat_budg)
#ifndef SCM
      call add_dim(cdl_latbudg,'shnhgm',3)
      call add_dim(cdl_latbudg,'lat_budg_plus3',jm_budg+3)
#endif
      call add_var(cdl_latbudg,'float area_budg(lat_budg) ;',
     &       units='m^2')
      call add_vardata(cdl_latbudg,'area_budg',dxyp_budg)

      cdl_j = cdl_latbudg  ! invoke a copy method later
      call add_dim(cdl_j,'ntype',ntype_out)
      call add_dim(cdl_j,'stype_strlen',16)
      call add_var(cdl_j,'char terrain(ntype,stype_strlen) ;')
      call add_vardata(cdl_j,'terrain',terrain)
      do k=1,kaj
        if(trim(stitle_j(k)).eq.'no output') cycle
        sname = 'J_'//trim(name_j(k))
        set_miss = iden_j(k).ne.0
        call add_var(cdl_j,
     &       'float '//trim(sname)//'(ntype,lat_budg) ;',
     &       long_name=trim(lname_j(k)),
     &       auxvar_string=
     &         'float '//trim(sname)//'_hemis(ntype,shnhgm) ;',
     &       units=trim(units_j(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
        call add_varline(cdl_j,
     &         trim(sname)//':fmt = "'//trim(fmt_j(k))//'" ;')
        call add_varline(cdl_j,
     &       trim(sname)//':stitle = "'//trim(stitle_j(k))//'" ;')
      enddo

c
c Declare the dimensions and metadata of AREG output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      call init_cdl_type('cdl_reg',cdl_reg)
      call add_dim(cdl_reg,'nreg',nreg)
      call add_dim(cdl_reg,'twenty_three',23)
      call add_dim(cdl_reg,'eight',8)
      call add_var(cdl_reg,'char namreg(twenty_three,eight) ;')
      do k=1,23!nreg
        namreg_1word(k) = namreg(1,k)//namreg(2,k)
      enddo
      call add_vardata(cdl_reg,'namreg',namreg_1word)
      do k=1,kaj
        if(trim(stitle_j(k)).eq.'no output') cycle
        if(trim(fmt_reg(k)).eq.'not computed') cycle
        sname = 'reg_'//trim(name_j(k))
        set_miss = iden_reg(k).ne.0
        call add_var(cdl_reg,
     &       'float '//trim(sname)//'(nreg) ;',
     &       long_name=trim(lname_j(k)),
     &       units=trim(units_j(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
        call add_varline(cdl_reg,
     &       trim(sname)//':fmt = "'//trim(fmt_reg(k))//'" ;')
        call add_varline(cdl_reg,
     &       trim(sname)//':stitle = "'//trim(stitle_j(k))//'" ;')
      enddo

c
c Declare the dimensions and metadata of CONSRV output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      cdl_consrv = cdl_latbudg  ! invoke a copy method later
      do k=1,kcmx
        sname = trim(name_consrv(k))
        call add_var(cdl_consrv,
     &       'float '//trim(sname)//'(lat_budg) ;',
     &       long_name=trim(title_con(k)))
#ifndef SCM
        call add_var(cdl_consrv,
     &       'float '//trim(sname)//'_hemis(shnhgm) ;')
#endif
      enddo

      return
      end subroutine j_defs


      subroutine ij_defs
      use constant
      use MODEL_COM
      use TimeConstants_mod, only : SECONDS_PER_DAY
      use DIAG_COM
      use DIAG_COM_RAD
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use fluxes, only : nisurf,atmice
#ifdef TRACERS_ON
      use TRACER_COM, only: gasex_index
      use oldtracer_mod, only: trname
#endif
      use cdl_mod
      use MDIAG_COM, only : make_timeaxis,sname_strlen,lname_strlen
      use geom
      use dynamics, only : do_gwdrag,ido_gwdrag
      use rad_com, only: nradfrc,TAero_aod_diag,cloud_rad_forc
      use rad_com, only: aer_rad_forc
#ifdef DETAILED_FIRE_OUTPUT
      use flammability_com, only: ij_flamV
      use ent_const, only : N_COVERTYPES
#endif
      USE SOCPBL, only : calc_wspdf
      use dictionary_mod
      implicit none
      integer :: i,k,kk,k1,l,n,ngx,nq
      character(len=16) :: ijstr
      real*8 x_dummy(im)
#ifdef ENT_DEBUG_DIAGS
      integer ent_k1, ent_k2
      character*3 :: ent_s1, ent_s2
#endif
      logical :: set_miss
! The following local variables are used in the definition of groups of
! 2D outputs collected into output fields having a third dimension.
! The groupings are applied if fuse_groups is true.
!@dbparam fuse_groups whether to output groups of related
!@+         diagnostics in arrays of rank+1
      logical :: fuse_groups=.false.
!@var index1 index of the first acc element of a group of outputs
!@var dim3info_index for a group using an already-defined third dimension,
!@+      this is index of the acc element containing information
!@+      about the third dimension
      integer, dimension(:), allocatable :: index1,dim3info_index
!@var coord3 the coordinate axes corresponding to the third dimensions.
      real*8, dimension(:,:), allocatable :: coord3 ! (1:ncoordvalues,k)
!@var name3 variable names used for grouped output
!@var dim3name the names of the third dimensions for grouped output
!@var dim3units units of coordinate axes of the third dimensions (if existing)
      character(len=sname_strlen), dimension(:), allocatable ::
     &     name3,dim3name,dim3units
!@var lname3 long names of grouped output variable
      character(len=lname_strlen), dimension(:), allocatable ::
     &     lname3
!@var ij_xxx_diminfo indices saved to set dim3info_index
      integer :: ij_cp_diminfo,ij_ghy_diminfo,ij_aer_diminfo
      character(len=80) :: varstr,varstrll,long_name,auxvar_string
      character(len=3), dimension(8), parameter :: aer_in_rad_name=
     &  (/'SO4','SEA','ANT','OCX','BCI','BCB','DUS','VOL'/)
      character(len=1) :: aer_band
      integer :: kr
c
      do k=1,kaij
         write(name_ij(k),'(a3,i3.3)') 'AIJ',k
         lname_ij(k) = 'unused'
         units_ij(k) = 'unused'
         ia_ij(k) = ia_src
         scale_ij(k) = 1.
         igrid_ij(k) = 1
         jgrid_ij(k) = 1
         ir_ij(k) = ir_pct
         denom_ij(k) = 0
      enddo

      do k=1,kaijmm
         write(name_ijmm(k),'(a6,i3.3)') 'AIJmm',k
         lname_ijmm(k) = 'unused'
         units_ijmm(k) = 'unused'
         scale_ijmm(k) = 1.
      enddo

      ! NB: fuse_groups only used for aij category at the moment,
      ! so no special aij-specific substring in the flag name yet
      if(is_set_param('fuse_groups')) then
        call get_param('fuse_groups',fuse_groups)
      endif
      allocate(index1(kaij),name3(kaij),dim3name(kaij),coord3(100,kaij))
      allocate(dim3info_index(kaij),lname3(kaij),dim3units(kaij))
      do k=1,kaij
        index1(k) = 0
        dim3info_index(k) = 0
        coord3(:,k) = -1d30
        dim3units(k) = ''
      enddo

c
      k=0
C**** AIJ diagnostic names:
C**** NAME     NO.    DESCRIPTION   (SCALE)*IDACC  LOCATION
C**** ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c
!**** Horizontal Fractions
      k=k+1 !
      IJ_POCEAN = k ! not accumulated
      lname_ij(k) = 'OCEAN FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'ocnfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_RSOI = k ! POICE (1)            1 GD
      atmice%IJ_RSOI = k
      lname_ij(k) = 'OCEAN/LAKE ICE COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'oicefr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_POPOCN = k
      lname_ij(k) = 'OPEN OCEAN FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'opocnfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_LK = k ! PLAKE                                       4 DA
      lname_ij(k) = 'LAKE FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'lakefr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_POPWAT = k
      lname_ij(k) = 'OPEN WATER FRACTION'    ! (focean+flake)*(1-rsi)
      units_ij(k) = '%'
      name_ij(k) = 'opwatfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PWATER = k
      lname_ij(k) = 'WATER FRACTION'    ! (focean+flake)
      units_ij(k) = '%'
      name_ij(k) = 'pwatfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PSOIL = k
      lname_ij(k) = 'SOIL FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'soilfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_BSFR = k
      lname_ij(k) = 'BARE SOIL FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'bsfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_VSFR = k
      lname_ij(k) = 'VEGETATION FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'vsfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1
      IJ_FVEG = k
      name_ij(k) = 'fr_veg' !
      lname_ij(k) = 'FRACTION OF VEGETATED SOIL'
      units_ij(k) = '%'
      scale_ij(k) = 100.
      ia_ij(k) = ia_src
      !ir_ij(k) = ir_0_3550
c
      k=k+1 !
      IJ_LI = k ! PLICE                                       4 DA
      lname_ij(k) = 'LAND ICE FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'landicefr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      IJ_PMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'PRESSURE FREQUENCY at ' // Trim(PMNAME(L))//'mb'  
         units_ij(k) = '%'
         name_ij(k) = 'p_freq_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 100.
      EndDo
c
      k=k+1 !
      IJ_RSNW = k ! PSNOW (1)            1 GD
      atmice%IJ_RSNW = k
      lname_ij(k) = 'SNOW COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'snowfr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_RSIT = k ! POICE+PLICE+(IF SNOW)PEARTH               4 DA
      atmice%IJ_RSIT = k
      lname_ij(k) = 'SNOW AND ICE COVERAGE'
      units_ij(k) = '%'
      name_ij(k) = 'snowicefr'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1
      IJ_FRMP = k
      lname_ij(k) = 'SEA ICE MELT POND FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'FRMP'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
      denom_ij(k) = IJ_RSOI

      k=k+1 !
      IJ_CLDCV = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'TOTAL CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'pcldt'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PMCCLD = k ! PCLDMC (1)  (COMPOSITE OVER ATMOSPHERE) 2 RD
      lname_ij(k) = 'CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'pmccld'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_TCLDI = k
      lname_ij(k) = 'FRACTION OF TIME FOR ISCCP CLOUD'
      units_ij(k) = '%'
      name_ij(k) = 'pcldt_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_CLDCV1 = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'TAU>1 CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'pcldt_tau1'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.

      k=k+1 !
      IJ_MCCVTP = k ! PCLDMC (1)  (PICK UP FROM MSTCNV)
      lname_ij(k) = 'MC TOP CLOUD COVER'
      units_ij(k) = ''
      name_ij(k) = 'mccvtp'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1 !
      IJ_MCCVBS = k ! PCLDMC (1)  (PICK UP FROM MSTCNV)
      lname_ij(k) = 'MC BASE CLOUD COVER'
      units_ij(k) = ''
      name_ij(k) = 'mccvbs'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1 !
      IJ_CLRSKY = k ! not accumulated
      lname_ij(k) = 'CLEAR SKY FRACTION'
      units_ij(k) = '%'
      name_ij(k) = 'clrsky'
      ia_ij(k) = ia_ij(IJ_CLDCV)
      scale_ij(k) = 100.
c
!**** Water Mass
      IJ_QPMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'SPECIFIC HUMIDITY at ' // Trim(PMNAME(L)) //'mb'
         units_ij(k) = 'g/kg'
         name_ij(k) = 'q_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 1d3
         ir_ij(k) = ir_0_18
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_QPMB1
         coord3(L,IJ_QPMB1) = PMB(L)
      EndDo
      name3(IJ_QPMB1) = 'qcp'
      lname3(IJ_QPMB1) = 'SPECIFIC HUMIDITY'
      dim3name(IJ_QPMB1) = 'pcp'
      dim3units(IJ_QPMB1) = 'mb'
      ij_cp_diminfo = IJ_QPMB1
c
      k=k+1 !
      IJ_QS   = k ! QS                                (NO PRT)  3 SF
      lname_ij(k) = 'SURFACE AIR SPECIFIC HUMIDITY'
      units_ij(k) = '10^-4 g/g'
      name_ij(k) = 'qsurf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.d4
      ir_ij(k) = ir_0_180
c
      IJ_RHPMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'RELATIVE HUMIDITY at ' // Trim(PMNAME(L)) //'mb'
         units_ij(k) = '%'
         name_ij(k) = 'rh_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 1d2
         ir_ij(k) = ir_pct
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_RHPMB1
         coord3(L,IJ_RHPMB1) = PMB(L) ! let RH coord be distinct for other const-press
      EndDo
      name3(IJ_RHPMB1) = 'rhcp'
      lname3(IJ_RHPMB1) = 'RELATIVE HUMIDITY'
      dim3name(IJ_RHPMB1) = 'prh'
      dim3units(IJ_RHPMB1) = 'mb'
c
      k=k+1 !
      IJ_RH1 = k !
      lname_ij(k) = 'LAYER 1 RELATIVE HUMIDITY'
      units_ij(k) = '%'
      name_ij(k) = 'rh_layer1'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1d2
      ir_ij(k) = ir_pct
c
      k=k+1 !
      IJ_RHs  = k ! RHs                               (NO PRT)  3 SF
      lname_ij(k) = 'SURFACE AIR RELATIVE HUMIDITY'
      units_ij(k) = '%'
      name_ij(k) = 'RHsurf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.d2
c
      k=k+1
      IJ_WMSUM = k ! CLOUD LIQUID WATER PATH (kg/m**2)             1 CL
      lname_ij(k) = 'CLOUD LIQUID WATER PATH'
      units_ij(k) = '.1 kg/m^2'
      name_ij(k) = 'clwp'
      ia_ij(k) = ia_src
      scale_ij(k) = 10.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_QM = k ! ATMOSPHERIC WATER VAPOUR COLUMN (kg/m**2)             1 CL
      lname_ij(k) = 'ATMOSPHERIC WATER VAPOUR COLUMN'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'qatm'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_CLDW = k ! CLOUD CONDENSED WATER COLUMN               1 CL
      lname_ij(k) = 'CLOUD CONDENSED WATER COLUMN'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'cldw'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1 !
      IJ_CLDI = k ! CLOUD CONDENSED ICE COLUMN                 1 CL
      lname_ij(k) = 'CLOUD CONDENSED ICE COLUMN'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'cldi'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1
      IJ_LWPrad = k ! LIQUID WATER PATH SEEN BY RADIATION (kg/m**2), includes precip
      lname_ij(k) = 'LIQUID WATER PATH SEEN BY RADIATION'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'LWPrad'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
c
      k=k+1
      IJ_IWPrad = k ! ICE WATER PATH SEEN BY RADIATION (kg/m**2), includes precip
      lname_ij(k) = 'ICE WATER PATH SEEN BY RADIATION'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'IWPrad'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
c
      k=k+1
      ij_tclssct = k ! cumulative layer-average stratiform liquid cloud opacity (not output)
      ia_ij(k) = ia_src
c
      k=k+1
      ij_rclssct = k
      lname_ij(k) = 'Large-scale cld droplet r_eff near cld top'
      units_ij(k) = 'micron'
      name_ij(k) = 'sscldtop_reffcl'
      ia_ij(k) = ia_src
      denom_ij(k) = ij_tclssct
c
#ifdef CLD_AER_CDNC
      k=k+1
      ij_nclssct = k
      lname_ij(k) = 'Large-scale cld droplet num conc near cld top'
      units_ij(k) = '#/cm^3'
      name_ij(k) = 'sscldtop_ncl'
      ia_ij(k) = ia_src
      denom_ij(k) = ij_tclssct
c
#endif
#ifdef TRACERS_AMP
#ifdef BLK_2MOM
      k=k+1
      ij_ccnssct = k
      lname_ij(k) = 'CCN activated near cld top'
      units_ij(k) = '#/cm^3'
      name_ij(k) = 'sscldtop_ccn'
      ia_ij(k) = ia_src
      denom_ij(k) = ij_tclssct
      scale_ij(k) = 1e-6
c
#endif
#endif
      k=k+1 !
      IJ_SNOW = k ! SNOW (KG/m**2)       1 GD
      atmice%IJ_SNOW = k
      lname_ij(k) = 'SNOW DEPTH'    ! 'SNOW MASS'
      units_ij(k) = 'mm H2O'
      name_ij(k) = 'snowdp'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_26_150
c
      k=k+1
      IJ_GBSSND = k
      name_ij(k) = 'bs_snowdp' !
      lname_ij(k) = 'SNOW DEPTH OVER BARE SOIL'
      units_ij(k) = 'mm H2O'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_3550
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      IJ_GVSSND = k
      name_ij(k) = 'vs_snowdp' !
      lname_ij(k) = 'SNOW DEPTH OVER VEG SOIL'
      units_ij(k) = 'mm H2O'
      scale_ij(k) = 1000.
      ia_ij(k) = ia_src
      ir_ij(k) = ir_0_3550
      denom_ij(k) = IJ_VSFR
c
      k=k+1 !
      IJ_ZSNOW = k ! snow thickness over all surface types
      atmice%IJ_ZSNOW = k
      lname_ij(k) = 'SNOW THICKNESS'
      units_ij(k) = 'm'
      name_ij(k) = 'zsnow'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_4
c
      k=k+1
      atmice%IJ_MLTP = k
      lname_ij(k) = 'SEA ICE MELT POND MASS'
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'MLTP'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = IJ_RSOI

      k=k+1 !
      IJ_BETA = k ! BETA (1)             1 GD
      lname_ij(k) = 'GROUND WETNESS (VEG ROOTS)'
      units_ij(k) = '%'
      name_ij(k) = 'beta'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_VSFR
c
      k=k+1 !
      IJ_GWTR = k ! WATER1+WATER2+ICE1+ICE2 (EARTH POINTS ONLY) 1 GD
      lname_ij(k) = 'TOTAL EARTH WATER' ! includes ice
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'gwtr'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_GICE = k ! ICE1+ICE2 (EARTH POINTS ONLY) 1 GD
      lname_ij(k) = 'TOTAL EARTH ICE' ! includes ice
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'gice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_GWTR1 = k ! WATER1+ICE1 (EARTH POINTS ONLY) 1 GD
      lname_ij(k) = 'WATER IN UPPER EARTH LAYER' ! includes ice
      units_ij(k) = 'kg/m^2'
      name_ij(k) = 'gwtr1'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_GBSW = k
      name_ij(k) = 'bs_wlay1' !
      lname_ij(k) = 'LAYER 1 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_wlay2' !
      lname_ij(k) = 'LAYER 2 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_wlay3' !
      lname_ij(k) = 'LAYER 3 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_wlay4' !
      lname_ij(k) = 'LAYER 4 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_wlay5' !
      lname_ij(k) = 'LAYER 5 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_wlay6' !
      lname_ij(k) = 'LAYER 6 BARE SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      IJ_GBSBET = k
      name_ij(k) = 'bs_beta' !
      lname_ij(k) = 'BARE SOIL WETNESS, BETA'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      IJ_GVSW = k
      name_ij(k) = 'vs_wcan' !
      lname_ij(k) = 'VEGETATION CANOPY SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_wlay1' !
      lname_ij(k) = 'LAYER 1 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_wlay2' !
      lname_ij(k) = 'LAYER 2 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_wlay3' !
      lname_ij(k) = 'LAYER 3 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_wlay4' !
      lname_ij(k) = 'LAYER 4 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_wlay5' !
      lname_ij(k) = 'LAYER 5 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_wlay6' !
      lname_ij(k) = 'LAYER 6 VEGETATED SOIL WATER'
      units_ij(k) = 'mm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      IJ_GVSWET = k
      name_ij(k) = 'vs_wetness' !
      lname_ij(k) = 'VEGETATED SOIL WETNESS'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      IJ_GBVSWT = k
      name_ij(k) = 'bvs_wet' !
      lname_ij(k) = 'BARE & VEGETATED SOIL WETNESS'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_GWTBL = k
      name_ij(k) = 'wtbl_depth' !
      lname_ij(k) = 'AVERAGE WATER TABLE DEPTH'
      units_ij(k) = 'm'
      ia_ij(k) = ia_srf
      scale_ij(k) = -1.
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_MWL = k
      lname_ij(k) = 'MASS OF LAKE AND RIVER WATER'
      units_ij(k) = '10^10 kg'
      name_ij(k) = 'mwl'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      IJ_GBETPEN = k
      name_ij(k) = 'beta_pen' !
      lname_ij(k) = 'PENMAN SOIL WETNESS, BETA'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_DTGDTS = k ! 18*(DEL(TG)/DEL(TS)-1),DEL=diurn_range     dly_ea
      lname_ij(k) = 'PLANT WATER STRESS'
      units_ij(k) = '1'
      name_ij(k) = 'plant_wstress'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.*30.
      ir_ij(k) = ir_m190_530
c
      k=k+1 !
      IJ_MCCON   = k !
      lname_ij(k) = 'MOIST CONV COUNT'
      units_ij(k) = '1'
      name_ij(k) = 'mccon'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_3550
c
!**** Salt
      k=k+1 !
      IJ_SSS = k               !      3 SF
      lname_ij(k) = 'SEA SURFACE SALINITY'    ! layer 1
      units_ij(k) = 'psu'
      name_ij(k) = 'sss'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      denom_ij(k) = IJ_POCEAN
c
      k=k+1
      atmice%IJ_SSI1 = k
      lname_ij(k) = 'SEA ICE SALINITY (MASS LAYER 1)'
      units_ij(k) = 'psu'
      name_ij(k) = 'SSI1'
      ia_ij(k) = ia_src
      scale_ij(k) = 1d3
      denom_ij(k) = IJ_RSOI

      k=k+1
      atmice%IJ_SSI2 = k
      lname_ij(k) = 'SEA ICE SALINITY (MASS LAYER 2)'
      units_ij(k) = 'psu'
      name_ij(k) = 'SSI2'
      ia_ij(k) = ia_src
      scale_ij(k) = 1d3
      denom_ij(k) = IJ_RSOI

      k=k+1 !
      atmice%IJ_STIO = k ! NET SALT AT ICE-OCEAN INTERFACE
      lname_ij(k) = 'NET ICE-OCEAN SALT'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'netst_icoc'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_POCEAN
c
!**** Pressure
      k=k+1 !
      IJ_CLDTPPR = k ! P-CLOUD TOP   (100 PA)                  2 RD
      lname_ij(k) = 'CLOUD TOP PRESSURE'
      units_ij(k) = 'mb'
      name_ij(k) = 'cldtpp'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
      denom_ij(k) = IJ_CLDCV
c
c     near cloud top P:     P at level down to which cloud opt.depth = 1
      k=k+1 !
      IJ_CLDT1P  = k ! P-CLOUD TOP   (100 PA)                  2 RD
      lname_ij(k) = 'CLOUD TAU=1 PRESSURE'
      units_ij(k) = 'mb'
      name_ij(k) = 'cldtpp_tau1'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
      denom_ij(k) = IJ_CLDCV1
c
      k=k+1 !
      IJ_MCCLDTP = k ! P-MC CLOUD TOP   (100 PA)
      lname_ij(k) = 'CONVECTIVE CLOUD TOP PRESSURE'
      units_ij(k) = 'mb'
      name_ij(k) = 'mccldtp'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
      denom_ij(k) = IJ_MCCVTP
c
      k=k+1 !
      IJ_MCCLDBS = k ! P-MC CLOUD BASE  (100 PA)
      lname_ij(k) = 'CONVECTIVE CLOUD BASE PRESSURE'
      units_ij(k) = 'mb'
      name_ij(k) = 'mccldbs'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
      denom_ij(k) = IJ_MCCVBS
c
      k=k+1 !
      IJ_PRES = k ! PIJ (100 PA)  (NO PRINTOUT)  4 DA
      lname_ij(k) = 'SURFACE PRESSURE'
      units_ij(k) = 'mb'
      name_ij(k) = 'prsurf'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
c
      k=k+1 !
      IJ_PRESQ = k ! PIJ (100 PA)  (NO PRINTOUT)  4 DA
      lname_ij(k) = 'SURFACE PRESSURE (INCL. Q)'
      units_ij(k) = 'mb'
      name_ij(k) = 'prsurfq'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
c
      k=k+1 !
      IJ_SLP  = k ! PSL (100 PA-1000)  (USING TS)             4 DA
      lname_ij(k) = 'SEA LEVEL PRESSURE'
      units_ij(k) = 'mb-1000'
      name_ij(k) = 'slp'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_SLPQ  = k ! PSL (100 PA-1000)  (USING TS)             4 DA
      lname_ij(k) = 'SEA LEVEL PRESSURE (INCL. Q)'
      units_ij(k) = 'mb-1000'
      name_ij(k) = 'slpq'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      K = K+1
      IJ_PVS = K
       NAME_IJ(K) = 'PVS'
      LNAME_IJ(K) = 'SURFACE VAPOR PRESSURE'
      UNITS_IJ(K) = 'mb'
      ia_ij(k) = IA_IJ(IJ_QS  )

!**** Temperature
      k=k+1 !
      IJ_BTMPW = k ! BTEMPW-TF (K-TF)                         2 RD
      lname_ij(k) = 'BRIGHTNESS TEMP THRU WNDW' ! window region
      units_ij(k) = 'C'
      name_ij(k) = 'btemp_window'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1
      IJ_PTROP = k
      lname_ij(k) = 'TROPOPAUSE PRESSURE (WMO)'
      units_ij(k) = 'mb'
      name_ij(k) = 'ptrop'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1
      IJ_TTROP = k
      lname_ij(k) = 'TROPOPAUSE TEMPERATURE (WMO)'
      units_ij(k) = 'K'
      name_ij(k) = 'ttrop'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1 !
      IJ_CLDTPT = k !
      lname_ij(k) = 'CLOUD TOP TEMPERATURE'
      units_ij(k) = 'C'
      name_ij(k) = 'cldtpt'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_CLDCV

c     near cloud top T:     T at level down to which cloud opt.depth = 1
      k=k+1 !
      IJ_CLDT1T = k !
      lname_ij(k) = 'CLOUD TAU=1 TEMPERATURE'
      units_ij(k) = 'C'
      name_ij(k) = 'cldtpt_tau1'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_CLDCV1
c
      k=k+1 !
      IJ_CTPI = k
      lname_ij(k) = 'CLOUD TOP PRESSURE (ISCCP)'
      units_ij(k) = 'mb'
      name_ij(k) = 'cldtpp_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_1775
      denom_ij(k) = IJ_TCLDI
c
      IJ_TPMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'TEMPERATURE at ' // Trim(PMNAME(L)) // 'mb'
         units_ij(k) = 'C'
         name_ij(k) = 't_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 1
         ir_ij(k) = ir_m80_28
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_TPMB1
      EndDo
      name3(IJ_TPMB1) = 'tcp'
      lname3(IJ_TPMB1) = 'TEMPERATURE'
      dim3info_index(IJ_TPMB1) = ij_cp_diminfo
c
      k=k+1 !
      IJ_TS   = k ! TS (K-TF)                                 3 SF
      lname_ij(k) = 'SURFACE AIR TEMPERATURE'
      units_ij(k) = 'C'
      name_ij(k) = 'tsurf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_TSLI = k ! SURF AIR TEMP OVER LAND ICE  (C)  NISURF*1 SF
      lname_ij(k) = 'SURF AIR TEMP OVER LAND ICE'
      units_ij(k) = 'C'
      name_ij(k) = 'tsurf_lndice'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_TG1  = k ! TG1 (K-TF)                                1 GD
      lname_ij(k) = 'GROUND TEMPERATURE'
      units_ij(k) = 'C'
      name_ij(k) = 'tgrnd'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1
      atmice%IJ_TSICE = k
      lname_ij(k) = 'SEA ICE SURFACE TEMPERATURE'
      units_ij(k) = 'K'
      name_ij(k) = 'ts_oice'
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_RSOI
c
      k=k+1
      atmice%IJ_TSI = k
      lname_ij(k) = 'SEA ICE TEMPERATURE (MASS LAYER 2)'
      units_ij(k) = 'C'
      name_ij(k) = 'TEMPSI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = IJ_RSOI

      k=k+1
      IJ_GVST = k
      name_ij(k) = 'can_temp' !
      lname_ij(k) = 'CANOPY TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_tlay1' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_tlay2' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 2 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_tlay3' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 3 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_tlay4' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 4 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_tlay5' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 5 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      name_ij(k) = 'vs_tlay6' !
      lname_ij(k) = 'VEGETATED SOIL LAYER 6 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      IJ_GBST = k
      name_ij(k) = 'bs_tlay1' !
      lname_ij(k) = 'BARE SOIL LAYER 1 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_tlay2' !
      lname_ij(k) = 'BARE SOIL LAYER 2 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_tlay3' !
      lname_ij(k) = 'BARE SOIL LAYER 3 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_tlay4' !
      lname_ij(k) = 'BARE SOIL LAYER 4 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_tlay5' !
      lname_ij(k) = 'BARE SOIL LAYER 5 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      name_ij(k) = 'bs_tlay6' !
      lname_ij(k) = 'BARE SOIL LAYER 6 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_BSFR
c
      k=k+1 !
      IJ_TGO  = k               !      3 SF
      lname_ij(k) = 'SEA SURFACE TEMPERATURE'    ! layer 1
      units_ij(k) = 'C'
      name_ij(k) = 'sst'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      IJ_TOC2 = k ! TGO2= TOCEAN(2)  (C)
      lname_ij(k) = 'OCEAN TEMPERATURE BELOW MIXED LAYER' ! lyr 2
      units_ij(k) = 'C'
      name_ij(k) = 'TOC2'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      ir_ij(k) = ir_m9_26
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      IJ_TGO2 = k ! TGO12= TOCEAN(3) (C)
      lname_ij(k) = 'OCEAN TEMPERATURE AT ANN-MAX MIXED-LAYER' ! layer 3
      units_ij(k) = 'C'
      name_ij(k) = 'TGO2'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      ir_ij(k) = ir_m9_26
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      IJ_TMNMX  = k ! MIN(DIURNAL MAX OF COMPOSITE TS)      12 MN
      lname_ij(k) = 'SURFC AIR TEMPERATURE: LOWEST DIURNAL HIGH'
      units_ij(k) = 'C'
      name_ij(k) = 'TMNMX'
      ia_ij(k) = ia_inst
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
c
      k=k+1 !
      IJ_TMAXE = k ! MAX TS OVER EARTH FOR CURRENT DAY (C)       dly_ea
      lname_ij(k) = 'SURFACE AIR TEMPERATURE: DIURNAL HIGH/SOIL'
      units_ij(k) = 'C'
      name_ij(k) = 'TMAXE'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_TMAXC = k ! MAX composite TS FOR CURRENT DAY (C)        dly_ea
      lname_ij(k) = 'SURFACE AIR TEMPERATURE: DIURNAL HIGH' ! composite
      units_ij(k) = 'C'
      name_ij(k) = 'TMAXC'
      ia_ij(k) = ia_12hr   ! really ia_24hr
      scale_ij(k) = 2.*1.  ! really 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = 0
c
      k=k+1 !
      IJ_TDSL = k ! DIURNAL DELTA TS (K) OVER SOIL (NO PRT)      dly_ea
      lname_ij(k) = 'DIURNAL SURF AIR TEMP RANGE OVER SOIL'
      units_ij(k) = 'K'
      name_ij(k) = 'dtdiurn_soil'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      ir_ij(k) = ir_0_18
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_TDCOMP = k
      lname_ij(k) = 'DIURNAL SURF AIR TEMP RANGE' ! composite
      units_ij(k) = 'C'
      name_ij(k) = 'dtdiurn'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.
      ir_ij(k) = ir_0_18
c
!**** Energy
      k=k+1
      IJ_HTSNOW = k
      lname_ij(k) = 'TOTAL LAND SNOW HEAT STORAGE'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'snow_heat' !
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_HTSOIL = k
      lname_ij(k) = 'TOTAL SOIL HEAT STORAGE'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'soil_heat' !
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_GML = k
      lname_ij(k) = 'ENTHALPY OF LAKE AND RIVER WATER'
      units_ij(k) = '10^15 J'
      name_ij(k) = 'gml'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-15
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      atmice%IJ_SIHC = k   ! note this is referenced to water at 0 deg C
      lname_ij(k) = 'SEA ICE HEAT CONTENT'
      units_ij(k) = 'J/m^2'
      name_ij(k) = 'hc_oice'
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_RSOI
c
      k=k+1 !
      IJ_STRNGTS = k ! MAX(0,65F-TS_daily_avg in F)              dly_ea
      lname_ij(k) = 'MONTHLY HEATING' ! monthly heating need ?
      units_ij(k) = 'degF days'
      name_ij(k) = 'heat_deg_days'
      ia_ij(k) = ia_12hr
      scale_ij(k) = 2.*30.
      ir_ij(k) = ir_0_3550
c
!**** Geopotential Height
      IJ_ZPMB1 = k+1     
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'HEIGHT at ' // Trim(PMNAME(L)) // 'mb'
         units_ij(k) = 'm'
         name_ij(k) = 'z_' // PMNAME(L)     
         ia_ij(k) = ia_dga
         scale_ij(k) = 1       
         ir_ij(k) = ir_m190_530
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_ZPMB1
      EndDo
      name3(IJ_ZPMB1) = 'zcp'
      lname3(IJ_ZPMB1) = 'HEIGHT'
      dim3info_index(IJ_ZPMB1) = ij_cp_diminfo
c
      k=k+1 !
      IJ_PBLHT   = k !
      lname_ij(k) = 'PBL HEIGHT'
      units_ij(k) = 'm'
      name_ij(k) = 'pblht'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_3550
c
      k=k+1 !
      IJ_SSH = k               !      3 SF
      lname_ij(k) = 'SEA SURFACE HEIGHT'
      units_ij(k) = 'm'
      name_ij(k) = 'ssh'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      atmice%IJ_MSI = k ! ACE2OI+ACE1I= (MSI2+MSI1)*POICE/RHOI (m)   1 GD
      lname_ij(k) = 'OCEAN ICE THICKNESS'
      units_ij(k) = 'm'
      name_ij(k) = 'ZSI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./RHOI
      ir_ij(k) = ir_0_4
      denom_ij(k) = IJ_RSOI

      k=k+1 !
      atmice%IJ_SITF = k
      lname_ij(k) = 'OCEAN ICE TIME FRACTION'
      units_ij(k) = '1.0'
      name_ij(k) = 'sitimefrac'
      scale_ij(k) = 1.

      k=k+1 !
      atmice%IJ_SIMASS = k
      lname_ij(k) = 'OCEAN ICE MASS PER AREA'
      units_ij(k) = 'kg/m2'
      name_ij(k) = 'simass'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.

      k=k+1 !
      atmice%IJ_SIVOL = k
      lname_ij(k) = 'OCEAN ICE VOLUME PER AREA'
      units_ij(k) = 'm'
      name_ij(k) = 'sivol'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./RHOI
c
!**** Velocity and Momentum
      k=k+1 !
      IJ_UJET = k ! UJET (M/S)                                4 DA
      lname_ij(k) = 'U COMPONENT OF JET WINDS'
      units_ij(k) = 'm/s'
      name_ij(k) = 'ujet'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      igrid_ij(k) = 1 ! now using a-grid winds
      jgrid_ij(k) = 1
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_VJET = k ! VJET (M/S)                                4 DA
      lname_ij(k) = 'V COMPONENT OF JET WINDS'
      units_ij(k) = 'm/s'
      name_ij(k) = 'vjet'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      igrid_ij(k) = 1 ! now using a-grid winds
      jgrid_ij(k) = 1
      ir_ij(k) = ir_m38_106
!
      IJ_UPMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'EASTWARD VELOCITY at ' // Trim(PMNAME(L)) //'mb'
         units_ij(k) = 'm/s'
         name_ij(k) = 'u_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 1       
         ir_ij(k) = ir_m38_106
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_UPMB1
      EndDo
      name3(IJ_UPMB1) = 'ucp'
      lname3(IJ_UPMB1) = 'EASTWARD VELOCITY'
      dim3info_index(IJ_UPMB1) = ij_cp_diminfo
!
      IJ_VPMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k) = 'NORTHWARD VELOCITY at ' // Trim(PMNAME(L))//'mb'
         units_ij(k) = 'm/s'
         name_ij(k) = 'v_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 1       
         ir_ij(k) = ir_m38_106
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_VPMB1
      EndDo
      name3(IJ_VPMB1) = 'vcp'
      lname3(IJ_VPMB1) = 'NORTHWARD VELOCITY'
      dim3info_index(IJ_VPMB1) = ij_cp_diminfo
c
      k=k+1 !
      IJ_US   = k ! US (M/S)                                  3 SF
      lname_ij(k) = 'U COMPONENT OF SURFACE AIR WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'usurf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_VS   = k ! VS (M/S)                                  3 SF
      lname_ij(k) = 'V COMPONENT OF SURFACE AIR WIND'
      units_ij(k) = 'm/s'
      name_ij(k) = 'vsurf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_WS   = k ! SURFACE WIND SPEED (M/S)                  3 SF
      lname_ij(k) = 'SURFACE WIND SPEED'
      units_ij(k) = 'm/s'
      name_ij(k) = 'wsurf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_GUSTI   = k !
      lname_ij(k) = 'GUSTI WIND'
      units_ij(k) = 'M/S'
      name_ij(k) = 'gusti'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_3550
c
      k=k+1 !
      IJ_TAUUS = k ! TAUUS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'U COMPON OF MOMENTUM SRF DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'tauus'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !
      IJ_TAUVS = k ! TAUVS (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'V COMPON OF MOMENTUM SRF DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'tauvs'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_m2650_950
c
      k=k+1 !
      IJ_TAUS = k ! TAUS  (MOM. SURF. DRAG) (kg/m**2) (NO PRT)  3 SF
      lname_ij(k) = 'MAG OF MOMENTUM SURFACE DRAG'
      units_ij(k) = 'g/m*s^2'
      name_ij(k) = 'tausmag'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1000.
      ir_ij(k) = ir_0_710
c
c     k=k+1 !
c     IJ_WLM  = k ! WIND SPEED IN TOP LAYER (M/S) before SDRAG 1 SD
c     lname_ij(k) = 'WIND SPEED IN TOP LAYER'
c     units_ij(k) = 'm/s'
c     name_ij(k) = 'WLM'
c     ia_ij(k) = ia_src
c     scale_ij(k) = 1.
c     igrid_ij(k) = 2
c     jgrid_ij(k) = 2
c     ir_ij(k) = ir_0_26_150
c
      IJ_OMEGAPMB1 = k+1
      Do L=1,KGZ
         k=k+1
         lname_ij(k)='DOWNWARD PRESSURE FLUX at '//Trim(PMNAME(L))//'mb'
         units_ij(k) = 'Pa/s'
         name_ij(k) = 'omega_' // PMNAME(L)
         ia_ij(k) = ia_dga
         scale_ij(k) = 1       
!        ir_ij(k) = ir_m38_106
         denom_ij(k) = IJ_PMB1 + L - 1
         index1(k) = IJ_OMEGAPMB1
      EndDo
      name3(IJ_OMEGAPMB1) = 'omegacp'
      lname3(IJ_OMEGAPMB1) = 'DOWNWARD PRESSURE FLUX'
      dim3info_index(IJ_OMEGAPMB1) = ij_cp_diminfo

!**** Vertical Mass Fluxes
      k=k+1
      IJ_H2OCH4 = k  !  1 GP
      lname_ij(k) = 'WATER DERIVED FROM CH4 OXIDATION IN STRATOSPHERE'
      units_ij(k) = '10^-6 mm/day'
      name_ij(k)  = 'H2O_from_CH4'
      ia_ij(k)    = ia_12hr  !  accumulated daily, 2* in scale
      scale_ij(k) = 2d6      
!
      k=k+1 !
      IJ_PREC = k ! PREC (mm/day)       1 CN
      lname_ij(k) = 'PRECIPITATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'prec'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_PRECGR = k ! PREC OVER EARTH (mm/day)       1 CN
      lname_ij(k) = 'PRECIPITATION OVER EARTH'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'pr_grnd'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_AINTRCP = k ! aintercep
      lname_ij(k) = 'PRECIPITATION INTERCEPTED BY CANOPY'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'prec_int_canopy'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_VSFR
c
      k=k+1 !
      IJ_PRECLI = k ! PREC OVER LAND ICE (mm/day)       1 CN
      lname_ij(k) = 'PRECIPITATION OVER LAND ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'pr_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_PRECOO = k ! PREC OVER OPEN OCEAN (mm/day)       1 CN
      lname_ij(k) = 'PRECIPITATION OVER OPEN OCEAN'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'pr_oocn'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      denom_ij(k) = IJ_POPWAT
c
      k=k+1 !
      IJ_PRECSI = k ! PREC OVER SEA ICE (mm/day)       1 CN
      lname_ij(k) = 'PRECIPITATION OVER SEA ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'pr_oice'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      denom_ij(k) = IJ_RSOI
c
      k=k+1 !
      IJ_PRECMC = k ! PREC MC (mm/day)       1 CN
      lname_ij(k) = 'CONVECTIVE PRECIPITATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'prec_mc'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_EVAP = k ! EVAP (mm/day)       1 SF
      lname_ij(k) = 'EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'evap'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_PEVAP = k ! POTENTIAL EVAPORATION (KG/m**2)         1 EA
      lname_ij(k) = 'POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'pot_evap'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_26_150
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_EVAPLI = k ! EVAP OVER LAND ICE  (KG/m**2)          1 GD
      lname_ij(k) = 'LAND ICE EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'evap_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_EVAPO = k ! EVAP*PWATER  (KG/m**2)                  1 GD
      lname_ij(k) = 'OPEN WATER EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'evap_ocn'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_POPWAT
c
      k=k+1 !
      IJ_EVAPI = k ! EVAP*POICE  (KG/m**2)                   1 GD
      lname_ij(k) = 'OCEAN ICE EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'evap_oice'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_RSOI
c
      k=k+1 !
      IJ_EVAPE = k ! EVAP OVER EARTH  (KG/m**2)              1 GD
      lname_ij(k) = 'SOIL EVAPORATION'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'evap_land'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
c     iw built-in
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_GBSEVP = k
      name_ij(k) = 'bs_evap' !
      lname_ij(k) = 'BARE SOIL EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_BSFR
c
      k=k+1
      IJ_GDCEVP = k
      name_ij(k) = 'drycan_evap' !
      lname_ij(k) = 'DRY CANOPY EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      IJ_GWCEVP = k
      name_ij(k) = 'wetcan_evap' !
      lname_ij(k) = 'WET CANOPY EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_VSFR
c
      k=k+1
      IJ_GEVPPEN = k
      name_ij(k) = 'pev_pen' !
      lname_ij(k) = 'PENMAN POTENTIAL EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_26_150
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_EVAPSN = k
      name_ij(k) = 'snow_evap' !
      lname_ij(k) = 'LAND SNOW EVAPORATION'
      units_ij(k) = 'mm/day'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_SNWF = k ! SNOW FALL  (KG/m**2)                     1 PR
      lname_ij(k) = 'SNOW FALL (H2O EQUIV)'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'snowfall'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_0_3_15
c
      k=k+1 !
      IJ_AFLMLT = k ! snow melt water  (KG/m**2 /s)                1 PG
      lname_ij(k) = 'SNOW MELT FLUX'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'snow_melt'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_RUNE = k ! RUN1 OVER EARTH  (KG/m**2)                1 PG
      lname_ij(k) = 'GROUND RUNOFF OVER SOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'runoff_soil'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_RUNLI = k ! RUN1 OVER LAND ICE  (KG/m**2) (NO PRT)    1 PG
      lname_ij(k) = 'SURFACE RUNOFF OVER LAND ICE'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'runoff_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_ARUNU = k ! ARUNU                                      1 EA
      lname_ij(k) = 'UNDERGROUND RUNOFF OVER SOIL'
      units_ij(k) = 'mm/day'
      name_ij(k) = 'runoff_ugrnd'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      atmice%IJ_SISNWF = k
      lname_ij(k) = 'SEA ICE SNOWFALL RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'prsn_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_RSOI
c
      k=k+1
      atmice%IJ_SIGRFR = k
      lname_ij(k) = 'SEA ICE FRAZIL GROWTH RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'grfraz_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_PWATER
c
      k=k+1
      atmice%IJ_SIGRCG = k
      lname_ij(k) = 'SEA ICE CONGELATION GROWTH RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'grcong_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_PWATER
c
      k=k+1
      atmice%IJ_SIGRLT = k   ! i.e. negative of lateral melt
      lname_ij(k) = 'SEA ICE LATERAL GROWTH RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'grlat_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_PWATER
c
      k=k+1
      atmice%IJ_SNTOSI = k   ! includes snow to ice and seawater to ice terms
      lname_ij(k) = 'SNOW ICE FORMATION RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'snotoice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_PWATER
c
      k=k+1
      atmice%IJ_SITOPMLT = k
      lname_ij(k) = 'SEA ICE SURFACE MELT RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'topmlt_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_PWATER
c
      k=k+1
      atmice%IJ_SIBOTMLT = k
      lname_ij(k) = 'SEA ICE BASAL MELT RATE'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'botmlt_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_PWATER
c
      K = K+1
      atmice%IJ_MSNFLOOD = K
      LNAME_IJ(K) = 'ICE MASS FROZEN by SNOW FLOOD'
      UNITS_IJ(K) = 'kg/s*m^2'
      NAME_IJ(K)  = 'msnflood'
      SCALE_IJ(K) = 1 / DTSRC
      IA_IJ(K)    = IA_SRC
      DENOM_IJ(K) = IJ_PWATER
c
      k=k+1 !
      IJ_FWOC = k ! NET FRESH WATER AT Z0 (INCL RIVERS + ICEBERGS)
      lname_ij(k) = 'NET FRESH WATER INTO OCEAN'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'netfw_osurf'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      atmice%IJ_FWIO = k ! NET FRESH WATER AT ICE-OCEAN INTERFACE
      lname_ij(k) = 'NET ICE-OCEAN FRESH WATER'
      units_ij(k) = 'kg/m^2/s'
      name_ij(k) = 'netfw_icoc'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      IJ_IMPMLI = k ! IMPLICIT MASS FLUX over LAND ICE (kg/s*m^2)           1 SF
      lname_ij(k) = 'IMPLICIT MASS FLUX over LAND ICE'
      units_ij(k) = 'kg/s*m^2'
      name_ij(k) = 'impm_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_IMPMGR = k ! IMPLICIT MASS FLUX over GROUND (kg/s*m^2)           1 SF
      lname_ij(k) = 'IMPLICIT MASS FLUX over GROUND'
      units_ij(k) = 'kg/s*m^2'
      name_ij(k) = 'impm_gr'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_IMPMKI = k ! IMPLICIT MASS FLUX over LAKE ICE (kg/s*m^2)           1 SF
      lname_ij(k) = 'IMPLICIT MASS FLUX over LAKE ICE'
      units_ij(k) = 'kg/s*m^2'
      name_ij(k) = 'impm_ki'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_RSOI
c
      IF (KOCEAN.ne.1) THEN
        k=k+1
        atmice%IJ_SMFX = k
        lname_ij(k) = 'SEA ICE IMPLICIT MASS FLUX'
        units_ij(k) = 'kg/m^2'
        name_ij(k) = 'SIMSFX'
        ia_ij(k) = ia_12hr
        scale_ij(k) = 2.
      END IF
c
      K = K+1
      IJ_MLKtoGR  = K
      LNAME_IJ(K) = 'MASS of EXPANDING LAKE SATURATES GROUND'
      UNITS_IJ(K) = 'kg/s*m^2'
      NAME_IJ(K)  = 'mlktogr'
      IA_IJ(K)    = IA_12HR
      SCALE_IJ(K) = 2 / SECONDS_PER_DAY
c
!**** Vertical Energy Fluxes
      K = K+1
      IJ_dSE_Dyn = K  !  ATM_PHASE1
      LNAME_IJ(K) = 'Static Energy Generated by Dynamics'
      UNITS_IJ(K) = 'W/m^2'
       NAME_IJ(K) = 'dSE_Dyn'
      SCALE_IJ(K) = 1 / DTsrc
         IA_IJ(K) = IA_SRC

      K = K+1
      IJ_dKE_Dyn = K  !  ATM_PHASE1
      LNAME_IJ(K) = 'Kinetic Energy Generated by Dynamics'
      UNITS_IJ(K) = 'W/m^2'
       NAME_IJ(K) = 'dKE_Dyn'
      SCALE_IJ(K) = 1 / DTsrc
         IA_IJ(K) = IA_SRC

      K = K+1
      IJ_dTE_Dyn = K  !  ATM_PHASE1
      LNAME_IJ(K) = 'Uniform Total Energy by Dynamics, Drag, Filter'
      UNITS_IJ(K) = 'W/m^2'
       NAME_IJ(K) = 'dTE_Dyn'
      SCALE_IJ(K) = 1 / DTsrc
         IA_IJ(K) = IA_SRC

      K = K+1
      ATMICE%IJ_dHSI_Dyn = K  !  SI_DIAGS
      LNAME_IJ(K) = 'Heat Content of Sea Ice by Advection'
      UNITS_IJ(K) = 'W/m^2'
       NAME_IJ(K) = 'dHSI_Dyn'
      SCALE_IJ(K) = 1 / DTsrc
         IA_IJ(K) = IA_SRC

      k=k+1
      IJ_SHDT = k ! SHDT (J/m**2)        1 SF
      lname_ij(k) = 'SENSIBLE HEAT FLUX'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'sensht'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m265_95
c
      k=k+1 !
      IJ_SHDTLI = k ! SHDT OVER LAND ICE  (J/m**2)           1 SF
      lname_ij(k) = 'SENS HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'sensht_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m38_106
      denom_ij(k) = IJ_LI
c
      k=k+1
      IJ_SISH = k
      lname_ij(k) = 'SEA ICE SENSIBLE HEAT FLUX'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'sh_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_RSOI
c
      K = K+1
      IJ_HWV = K
       NAME_IJ(K) = 'HWV'
      LNAME_IJ(K) = 'LATENT HEAT FLUX'
      UNITS_IJ(K) = 'W/m^2'
      ia_ij(k) = IA_IJ(IJ_EVAP)
      scale_ij(k) = 2500000 / DTsrc

      k=k+1 !
      IJ_EVHDT = k ! EVHDT OVER LAND ICE  (J/m**2)           1 SF
      lname_ij(k) = 'LATENT HEAT FLUX OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'latht_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m38_106
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_NETH = k    ! SRHDT+TRHDT+SHDT+EVHDT+ENRGP (J/m**2)   1 SC
      lname_ij(k) = 'NET HEATING AT GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'netht_grnd'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_AERUNS = k ! energy of surf runoff                      1 EA
      lname_ij(k) = 'HEAT OF SURFACE RUNOFF OVER SOIL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'e_runoff_surf'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_AERUNU = k ! energy of underground runoff              1 EA
      lname_ij(k) = 'HEAT OF UNDERGROUND RUNOFF OVER SOIL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'e_runoff_ugrnd'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0/DTsrc
      ir_ij(k) = ir_m1_3
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      atmice%IJ_F0OI = k ! F0DT*POICE, NET HEAT AT Z0  (J/m**2)     1 GD
      lname_ij(k) = 'NET HEAT INTO OCEAN ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'netht_oice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
c     iw built-in
      ir_ij(k) = ir_m530_190
      denom_ij(k) = IJ_RSOI
c
      k=k+1 !
      IJ_F0LI = k ! F0DT, NET HEAT AT Z0 OVER LAND ICE  (J/m**2) 1 GD
      lname_ij(k) = 'NET HEAT INTO LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'netht_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m530_190
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_F0E  = k ! F0DT, NET HEAT AT Z0 OVER EARTH  (J/m**2) 1 GD
      lname_ij(k) = 'NET HEAT INTO LAND SURFACE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'netht_land'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m530_190
      denom_ij(k) = IJ_PSOIL
c
      K = K+1
      atmice%IJ_HSNFLOOD = K
      LNAME_IJ(K) = 'ICE HEAT FROZEN by SNOW FLOOD'
      UNITS_IJ(K) = 'W/m^2'
      NAME_IJ(K)  = 'hsnflood'
      SCALE_IJ(K) = 1 / DTSRC
      IA_IJ(K)    = IA_SRC
      DENOM_IJ(K) = IJ_PWATER
c
      k=k+1 !
      IJ_F0OC = k ! NET HEAT INTO OCEAN (INCL. RIVERS + ICEBERGS)
      lname_ij(k) = 'NET HEAT INTO OCEAN'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'netht_osurf'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
c     iw built-in
      ir_ij(k) = ir_m530_190
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      atmice%IJ_HTIO = k ! NET HEAT AT ICE-OCEAN INTERFACE
      lname_ij(k) = 'NET ICE-OCEAN HEAT'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'netht_icoc'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_POCEAN
c
      k=k+1 !
      IJ_F1LI = k ! F1DT OVER LAND ICE  (J/m**2)             1 PG
      lname_ij(k) = 'CONDUCTION AT LYR 1 BOTTOM OVER LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'F1LI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m530_190
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_IMPHLI = k ! IMPLICIT HEAT FLUX over LAND ICE (W/m^2)           1 SF
      lname_ij(k) = 'IMPLICIT HEAT FLUX over LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'imph_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_LI
c
      k=k+1 !
      IJ_IMPHGR = k ! IMPLICIT HEAT FLUX over GROUND (W/m^2)           1 SF
      lname_ij(k) = 'IMPLICIT HEAT FLUX over GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'imph_gr'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !
      IJ_IMPHKI = k ! IMPLICIT HEAT FLUX over LAKE ICE (W/m^2)           1 SF
      lname_ij(k) = 'IMPLICIT HEAT FLUX over LAKE ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'imph_ki'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      denom_ij(k) = IJ_RSOI
c
      K = K+1
      IJ_HLKtoGR  = K
      LNAME_IJ(K) = 'HEAT of EXPANDING LAKE SATURATES GROUND'
      UNITS_IJ(K) = 'W/m^2'
      NAME_IJ(K)  = 'hlktogr'
      IA_IJ(K)    = IA_12HR
      SCALE_IJ(K) = 2 / SECONDS_PER_DAY
c
!**** Radiation
      k=k+1 !
      IJ_SRINCP0 = k ! SRINCP0 (W/m**2)                        3 SR
      lname_ij(k) = 'INCIDENT SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'incsw_toa'
      ia_ij(k) = ia_src
      scale_ij(k) = 1
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_SRNFP0 = k ! SRNFP0 (W/m**2)                         2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srnf_toa'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_CLR_SRUPTOA = k ! SRUFP0*CLRSKY (W/m**2)            2 RD
      lname_ij(k) = 'CLR SKY OUT SOLAR RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'swup_toa_clrsky'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_CLRSKY
c
      k=k+1 !
      IJ_SRNTP = k   ! SRNTP (W/m**2)                          2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, TROPO'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srn_tropo'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_CLR_SRNTP = k   ! SRNTP_CLR (W/m**2)                   2 RD
      lname_ij(k) = 'NET CLEAR-SKY SOLAR RADIATION, TROPO'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srnclr_tropo'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m190_530
      denom_ij(k) = IJ_CLRSKY
c
      k=k+1 !
      IJ_SRINCG = k ! SRINCG (W/m**2)                         2 RD
      lname_ij(k) = 'INCIDENT SOLAR RADIATION, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'incsw_grnd'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_CLR_SRINCG = k ! SRINCG*CLRSKY (W/m**2)            2 RD
      lname_ij(k) = 'CLR SKY INCIDENT SOLAR RADIATION, SRF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'incsw_grnd_clrsky'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_CLRSKY
c
      k=k+1
      IJ_SISWD = k
      lname_ij(k) = 'SEA ICE DOWNWARD SHORTWAVE RADIATION'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'swd_oice'
      scale_ij(k) = 1.
      ia_ij(k) = ia_rad
      denom_ij(k) = IJ_RSOI
c
      k=k+1
      IJ_SISWU = k
      lname_ij(k) = 'SEA ICE UPWARD SHORTWAVE RADIATION'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'swu_oice'
      scale_ij(k) = 1.
      ia_ij(k) = ia_rad
      denom_ij(k) = IJ_RSOI
c
      k=k+1 !
      IJ_SRNFG = k   ! SRNFG (W/m**2)                          2 RD
      lname_ij(k) = 'NET SOLAR RADIATION, SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srnf_grnd'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_CLR_SRNFG = k ! SRNFG*CLRSKY (W/m**2)            2 RD
      lname_ij(k) = 'CLR SKY NET SOLAR RADIATION, SRF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'swn_grnd_clrsky'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_CLRSKY
c
      k=k+1 !
      IJ_SRREF = k ! PLAVIS*S0*COSZ (W/m**2)                  2 RD
      lname_ij(k) = 'REFLECTED SOLAR RADIATION IN VISUAL'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srrefvis_toa'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_SRVISSURF = k ! SRVISSURF (W/m**2)                         2 RD
      lname_ij(k) = 'TOTAL VISIBLE SOLAR RADIATION AT SURFACE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srvissurf'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_SRVDIR = k ! FSRDIR*SRVISSURF (W/m**2)                  2 RD
      lname_ij(k) = 'DIRECT VISIBLE SOLAR RADIATION AT SURFACE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srvdir'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1
      IJ_SRVIS = k ! ALBVIS*S0*COSZ (W/m**2)                  2 RD
      lname_ij(k) = 'REFLECTED SOLAR RADIATION IN VISUAL AT SURF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srrefvis_grnd'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_TRNFP0 = k ! TRNFP0 (W/m**2)                         2 RS
      lname_ij(k) = 'NET THERMAL RADIATION, TOA'   ! >0 if down !
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trnf_toa'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_CLR_TRUPTOA = k ! TRUFP0*CLRSKY (W/m**2)            2 RD
      lname_ij(k) = 'CLR SKY OUT THERMAL RADIATION, TOA'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trup_toa_clrsky'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_CLRSKY
c                          non-negligible clouds (opt.depth>1)
      k=k+1 !
      IJ_TRNTP = k   ! TRNTP (W/m**2)                          2 RD
      lname_ij(k) = 'NET THERMAL RADIATION, TROPO'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trn_tropo'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
c
      k=k+1 !
      IJ_CLR_TRNTP = k   ! TRNTP_CLR (W/m**2)                   2 RD
      lname_ij(k) = 'NET CLEAR-SKY THERMAL RADIATION, TROPO'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trnclr_tropo'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_CLRSKY
c
      k=k+1 !
      IJ_TRSDN = k ! (W/m**2)                         2 RS
      lname_ij(k) = 'THERMAL RADIATION DOWN, SURF'   ! >0 if down !
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trdn_surf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_CLR_TRDNG = k ! TRDNG*CLRSKY (W/m**2)            2 RD
      lname_ij(k) = 'CLR SKY THERMAL RADIATION DOWN, SRF'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trdn_grnd_clrsky'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_710
      denom_ij(k) = IJ_CLRSKY
c
      k=k+1 !
      IJ_TRSUP = k ! (W/m**2)                         2 RS
      lname_ij(k) = 'THERMAL RADIATION UP, SURF'   ! >0 if up !
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trup_surf'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_m530_190
c
      K = K+1
      IJ_RTSE = K
       NAME_IJ(K) = 'RTSE'
      LNAME_IJ(K) = 'NET THERMAL RADIATION from SURFACE'    ! >0 if up !
      UNITS_IJ(K) = 'W/m^2'
      ia_ij(k) = IA_IJ(IJ_TRSUP)

      k=k+1
      IJ_TRHDT = k ! TRHDT OVER LAND ICE  (J/m**2)           1 SF
      lname_ij(k) = 'NET THERMAL RADIATION INTO LAND ICE'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'trht_lndice'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m38_106
      denom_ij(k) = IJ_LI
c
      k=k+1
      IJ_SILWD = k
      lname_ij(k) = 'SEA ICE DOWNWARD LONGWAVE RADIATION'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'lwd_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_RSOI
c
      k=k+1
      IJ_SILWU = k
      lname_ij(k) = 'SEA ICE UPWARD LONGWAVE RADIATION'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'lwu_oice'
      scale_ij(k)=1/DTsrc
      ia_ij(k) = ia_src
      denom_ij(k) = IJ_RSOI
c
      k=k+1 !
      IJ_RNFP1 = k ! RNFP1 (W/m**2)
      lname_ij(k) = 'NET RADIATION, P1'   ! >0 if down !
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'rnf_p1'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      ir_ij(k) = ir_m530_190
c
      k=k+1 !
      IJ_SRTR = k ! SRHDT+TRHDT (J/m**2)                    1 RD/SF
      lname_ij(k) = 'NET RADIATION AT GROUND'
      units_ij(k) = 'W/m^2'
      name_ij(k) = 'srtrnf_grnd'
      ia_ij(k) = ia_src
      scale_ij(k) = 1./DTsrc
      ir_ij(k) = ir_m95_265
c
      if (nradfrc>0) then
        k=k+1 !
        IJ_SWCRF = k ! SW cloud radiative forcing (W/m**2)      2 RD
        lname_ij(k) = 'SW CLOUD RADIATIVE FORCING, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swcrf_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m265_95
c
        k=k+1 !
        IJ_LWCRF = k ! LW cloud radiative forcing (W/m**2)      2 RD
        lname_ij(k) = 'LW CLOUD RADIATIVE FORCING, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwcrf_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265

        if (TAero_aod_diag > 0) then
          do n=1,8
            do kr=1,6
              if (TAero_aod_diag == 2 .and. kr /= 6) cycle ! only save band6
              write(aer_band,'(i1)')kr

              k=k+1
              ij_nintaerext(kr,n)=k
              name_ij(k)='NINText_'//
     &                    aer_in_rad_name(n)//'_band'//aer_band
              lname_ij(k)='Extinction for '//
     &                    aer_in_rad_name(n)//' band'//aer_band
              units_ij(k)='1'
              ia_ij(k)=ia_rad_frc
              scale_ij(k)=1.

              k=k+1
              ij_nintaersca(kr,n)=k
              name_ij(k)='NINTsca_'//
     &                    aer_in_rad_name(n)//'_band'//aer_band
              lname_ij(k)='Scattering for '//
     &                    aer_in_rad_name(n)//' band'//aer_band
              units_ij(k)='1'
              ia_ij(k)=ia_rad_frc
              scale_ij(k)=1.

              k=k+1
              ij_nintaerasy(kr,n)=k
              name_ij(k)='NINTasy_'//
     &                    aer_in_rad_name(n)//'_band'//aer_band
              lname_ij(k)='Asymmetry parameter for '//
     &                    aer_in_rad_name(n)//' band'//aer_band
              units_ij(k)='1'
              ia_ij(k)=ia_rad_frc
              scale_ij(k)=1.
            enddo ! kr
          enddo ! n
        endif

        k=k+1 !
        IJ_SW_CS_noA = k ! SW aerosol free radiative forcing (W/m**2)      2 RD
        lname_ij(k) = 'SW Clear Sky Aerosol FREE RADIATIVE FORCING, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'sw_cs_noa_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.

        k=k+1 !
        IJ_LW_CS_noA = k ! LW aerosol free radiative forcing (W/m**2)      2 RD
        lname_ij(k) = 'LW Clear Sky Aerosol FREE RADIATIVE FORCING, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lw_cs_noa_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.

        k=k+1 !
        IJ_SW_AS_noA = k ! SW aerosol free radiative forcing (W/m**2)      2 RD
        lname_ij(k) = 'SW All Sky Aerosol FREE RADIATIVE FORCING, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'sw_as_noa_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.

        k=k+1 !
        IJ_LW_AS_noA = k ! LW aerosol free radiative forcing (W/m**2)      2 RD
        lname_ij(k) = 'LW All Sky Aerosol FREE RADIATIVE FORCING, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lw_as_noa_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
c

         if (cloud_rad_forc.eq.2) then
        k=k+1 !
        IJ_SWCRF2 = k ! SW cloud radiative forcing (W/m**2) without aerosols and Ozone
        lname_ij(k) = 'SW CLOUD RF NO AER NO OX, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swcrf_toa2'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m265_95
c
        k=k+1 !
        IJ_LWCRF2 = k ! LW cloud radiative forcing (W/m**2) without aerosols and Ozone
        lname_ij(k) = 'LW CLOUD RF NO AER NO Ox, TOA'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwcrf_toa2'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265
         end if

c
        if (aer_rad_forc > 0) then
          IJ_SWAERRF = k+1   ! TOA SW aerosol rad forcing (W/m**2)
          DO N=1,8
            k=k+1
            lname_ij(k) = 'SW AER RADIATIVE FORCING, TOA N='//char(N+48)
            units_ij(k) = 'W/m^2'
            name_ij(k) = 'swaerrf_toa_'//char(n+48)
            ia_ij(k) = ia_rad_frc
            scale_ij(k) = 1.
            ir_ij(k) = ir_m95_265
          END DO
c
          IJ_LWAERRF = k+1   ! TOA LW aerosol rad forcing (W/m**2)
          DO N=1,8
            k=k+1
            lname_ij(k) = 'LW AER RADIATIVE FORCING, TOA N='//char(N+48)
            units_ij(k) = 'W/m^2'
            name_ij(k) = 'lwaerrf_toa_'//char(n+48)
            ia_ij(k) = ia_rad_frc
            scale_ij(k) = 1.
            ir_ij(k) = ir_m95_265
          END DO
c
          IJ_SWAERSRF = k+1   ! Surf SW aerosol rad forcing (W/m**2)
          DO N=1,8
            k=k+1
           lname_ij(k) = 'SW AER RADIATIVE FORCING, SURF N='//char(N+48)
            units_ij(k) = 'W/m^2'
            name_ij(k) = 'swaerrf_surf_'//char(n+48)
            ia_ij(k) = ia_rad_frc
            scale_ij(k) = 1.
            ir_ij(k) = ir_m95_265
          END DO
c
          IJ_LWAERSRF = k+1   ! Surf LW aerosol rad forcing (W/m**2)
          DO N=1,8
            k=k+1
           lname_ij(k) = 'LW AER RADIATIVE FORCING, SURF N='//char(N+48)
            units_ij(k) = 'W/m^2'
            name_ij(k) = 'lwaerrf_surf_'//char(n+48)
            ia_ij(k) = ia_rad_frc
            scale_ij(k) = 1.
            ir_ij(k) = ir_m95_265
          END DO
c
          IJ_SWAERABS = k+1   ! Atm. abs. by aerosol (W/m**2)
          DO N=1,8
            k=k+1
            lname_ij(k) = 'SW AER ATMOS. ABSORPTION N='//char(N+48)
            units_ij(k) = 'W/m^2'
            name_ij(k) = 'swaeraa_toa_'//char(n+48)
            ia_ij(k) = ia_rad_frc
            scale_ij(k) = 1.
            ir_ij(k) = ir_m95_265
          END DO
c
          IJ_LWAERABS = k+1   ! Atm. abs. by aerosol (W/m**2)
          DO N=1,8
            k=k+1
            lname_ij(k) = 'LW AER ATMOS. ABSORPTION N='//char(N+48)
            units_ij(k) = 'W/m^2'
            name_ij(k) = 'lwaeraa_toa_'//char(n+48)
            ia_ij(k) = ia_rad_frc
            scale_ij(k) = 1.
            ir_ij(k) = ir_m95_265
          END DO
c
          k=k+1
          IJ_SWAERRFNT = k   ! NET TOA SW aerosol rad forcing (W/m**2)
          lname_ij(k) = 'SW AER RADIATIVE FORCING, TOA NET'
          units_ij(k) = 'W/m^2'
          name_ij(k) = 'swaerrf_toa_net'
          ia_ij(k) = ia_rad_frc
          scale_ij(k) = 1.
          ir_ij(k) = ir_m95_265
c
          k=k+1
          IJ_LWAERRFNT = k   ! NET TOA LW aerosol rad forcing (W/m**2)
          lname_ij(k) = 'LW AER RADIATIVE FORCING, TOA NET'
          units_ij(k) = 'W/m^2'
          name_ij(k) = 'lwaerrf_toa_net'
          ia_ij(k) = ia_rad_frc
          scale_ij(k) = 1.
          ir_ij(k) = ir_m95_265
c
          k=k+1
          IJ_SWAERSRFNT = k   ! NET Surf SW aerosol rad forcing (W/m**2)
          lname_ij(k) = 'SW AER RADIATIVE FORCING, SURF NET'
          units_ij(k) = 'W/m^2'
          name_ij(k) = 'swaerrf_surf_net'
          ia_ij(k) = ia_rad_frc
          scale_ij(k) = 1.
          ir_ij(k) = ir_m95_265
c
          k=k+1
          IJ_LWAERSRFNT = k   ! NET Surf LW aerosol rad forcing (W/m**2)
          lname_ij(k) = 'LW AER RADIATIVE FORCING, SURF NET'
          units_ij(k) = 'W/m^2'
          name_ij(k) = 'lwaerrf_surf_net'
          ia_ij(k) = ia_rad_frc
          scale_ij(k) = 1.
          ir_ij(k) = ir_m95_265
c
          k=k+1               ! unused ????
          IJ_SWAERABSNT = k   ! NET Atm. abs. by aerosol (W/m**2)
          lname_ij(k) = 'SW AER ATMOS. ABSORPTION NET'
          units_ij(k) = 'W/m^2'
          name_ij(k) = 'swaeraa_toa_net'
          ia_ij(k) = ia_rad
          scale_ij(k) = 1.
          ir_ij(k) = ir_m95_265
c
          k=k+1                     ! unused ????
          IJ_LWAERABSNT = k         ! NET Atm. abs. by aerosol (W/m**2)
          lname_ij(k) = 'LW AER ATMOS. ABSORPTION NET'
          units_ij(k) = 'W/m^2'
          name_ij(k) = 'lwaeraa_toa_net'
          ia_ij(k) = ia_rad
          scale_ij(k) = 1.
          ir_ij(k) = ir_m95_265
        endif
c
        k=k+1 !
        IJ_SWDCLS = k ! SW clear-sky down radiation surf (W/m**2) 2 RD
        lname_ij(k) = 'SW CLR-SKY DOWNWARD RADIATION, SURFACE METHOD 2'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swdcls'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265
c
        k=k+1 !
        IJ_SWNCLS = k ! SW clear-sky net radiation surf (W/m**2) 2 RD
        lname_ij(k) = 'SW CLR-SKY NET RADIATION, SURFACE METHOD 2'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swncls'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265
c
        k=k+1 !
        IJ_LWDCLS = k ! LW clear-sky down radiation surf (W/m**2) 2 RD
        lname_ij(k) = 'LW CLR-SKY DOWNWARD RADIATION SURFACE METHOD 2'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwdcls'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265
c
        k=k+1 !
        IJ_SWNCLT = k ! SW clear-sky net radiation TOA (W/m**2) 2 RD
        lname_ij(k) = 'SW CLR-SKY NET RADIATION TOA METHOD 2'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swnclt'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265
c
        k=k+1 !
        IJ_LWNCLT = k ! LW clear-sky net radiation TOA (W/m**2) 2 RD
        lname_ij(k) = 'LW CLR-SKY NET RADIATION TOA METHOD 2'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwnclt'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
        ir_ij(k) = ir_m95_265
      endif
c
!**** Horizontal Mass Fluxes
      k=k+1 !
      IJ_FMU  = k ! EAST-WEST MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'EAST-WEST MASS FLUX'
      units_ij(k) = '10^10 kg/s'
      name_ij(k) = 'fmu'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10*100.*BYGRAV/DTsrc
      igrid_ij(k) = 2
      jgrid_ij(k) = 1
      ir_ij(k) = ir_m38_106
c
      k=k+1 !
      IJ_FMV  = k ! NORTH-SOUTH MASS FLUX (KG/S) 100./GRAV/3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH MASS FLUX'
      units_ij(k) = '10^10 kg/s'
      name_ij(k) = 'fmv'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10*100.*BYGRAV/DTsrc
      igrid_ij(k) = 1
      jgrid_ij(k) = 2
c
      k=k+1 !
      IJ_PUQ  = k ! P*U*Q (VERT. INTEGRATED) (100 PA*M/S) 4 DA
      lname_ij(k) = 'EAST-WEST HUMIDITY FLUX (VERT SUM)'
      units_ij(k) = 'mb*m/s'
      name_ij(k) = 'puq'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      igrid_ij(k) = 2
      jgrid_ij(k) = 1
      ir_ij(k) = ir_m45_130
c
      k=k+1 !
      IJ_PVQ  = k ! P*V*Q (VERT. INTEGRATED) (100 PA*M/S) 4 DA
      lname_ij(k) = 'NORTH-SOUTH HUMIDITY FLUX (VERT SUM)'
      units_ij(k) = 'mb*m/s'
      name_ij(k) = 'pvq'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.
      igrid_ij(k) = 1
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m45_130
c
      k=k+1 !
      IJ_MRVR = k ! Mass Inflow by Rivers (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'Mass Inflow by Rivers'
      units_ij(k) = '10^5 kg/s'
      name_ij(k) = 'MRVR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-5/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_MRVRO = k ! Mass Outflow by Rivers (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'Mass Outflow by Rivers'
      units_ij(k) = '10^5 kg/s'
      name_ij(k) = 'MRVRO'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-5/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_MICB = k ! Mass Inflow by Icebergs (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'Mass Inflow by Icebergs'
      units_ij(k) = '10^5 kg/s'
      name_ij(k) = 'MICB'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-5/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1
      atmice%IJ_MUSI=k
      lname_ij(k)="Sea ice EW mass flux"
      name_ij(k)="musi"
      units_ij(k)="10^7 kg/s"
      scale_ij(k)=1d-7/dtsrc
      igrid_ij(k)=2
c
      k=k+1
      atmice%IJ_MVSI=k
      lname_ij(k)="Sea ice NS mass flux"
      name_ij(k)="mvsi"
      units_ij(k)="10^7 kg/s"
      scale_ij(k)=1d-7/dtsrc       
      jgrid_ij(k)=2
c
      k=k+1       
      atmice%IJ_SUSI=k
      lname_ij(k)="Sea ice EW salt flux"
      name_ij(k)="susi"  
      units_ij(k)="10^3 kg/s"
      scale_ij(k)=1d-3/dtsrc
      igrid_ij(k)=2
c      
      k=k+1
      atmice%IJ_SVSI=k
      lname_ij(k)="Sea ice NS salt flux"
      name_ij(k)="svsi"
      units_ij(k)="10^3 kg/s"
      scale_ij(k)=1d-3/dtsrc
      jgrid_ij(k)=2
c
!**** Horizontal Energy Fluxes
      k=k+1 !
      IJ_DSEV  = k
       ! P4*(SHA*T4+Z4)*V1*DSIG*DXV (100 W*M/S**2) (UV GRID) 4 DA
      lname_ij(k) = 'TOTAL NT DRY STAT ENRGY' !  NT: NORTHWARD TRANSPORT
      units_ij(k) = '10^14 W'
      name_ij(k) = 'nt_dse'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1.d-14*100.*BYGRAV
      igrid_ij(k) = 2
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m95_265
c
      k=k+1 !
      IJ_FGZU = k ! EAST-WEST GEOPOTENTIAL FLUX (W)    /3600.*1 DY
      lname_ij(k) = 'EAST-WEST GEOPOTENTIAL FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'fgzu'
      ia_ij(k) = ia_src
      SCALE_IJ(K) = 1d-10 / DTSRC
      igrid_ij(k) = 2
      jgrid_ij(k) = 1
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_FGZV = k ! NORTH-SOUTH GEOPOTENTIAL FLUX (W)  /3600.*1 DY
      lname_ij(k) = 'NORTH-SOUTH GEOPOTENTIAL FLUX'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'fgzv'
      ia_ij(k) = ia_src
      SCALE_IJ(K) = 1d-10 / DTSRC
      igrid_ij(k) = 1
      jgrid_ij(k) = 2
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !     NOTE: INFLOW IS DEFINED AS MASS ENTERING A BOX
      IJ_ERVR = k ! Energy Inflow by Rivers (10**10 W) E-10/DTS*1 RV
      lname_ij(k) = 'Energy Inflow by Rivers'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'ERVR'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !    NOTE: OUTFLOW IS DEFINED AS MASS LEAVING A BOX
      IJ_ERVRO = k ! Energy Outflow by Rivers (10**10 W) E-10/DTS*1 RV
      lname_ij(k) = 'Energy Outflow by Rivers'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'ERVRO'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1 !
      IJ_EICB = k ! Energy Inflow by Icebergs (10**10 W) E-10/DTS*1 RV
      lname_ij(k) = 'Energy Inflow by Icebergs'
      units_ij(k) = '10^10 W'
      name_ij(k) = 'EICB'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-10/DTsrc
      ir_ij(k) = ir_m1325_475
c
      k=k+1
      atmice%IJ_HUSI=k
      lname_ij(k)="Sea ice EW heat flux"
      name_ij(k)="husi"
      units_ij(k)="10^12 W"
      scale_ij(k)=1d-12/dtsrc
      igrid_ij(k)=2
c
      k=k+1
      atmice%IJ_HVSI=k
      lname_ij(k)="Sea ice NS heat flux"
      name_ij(k)="hvsi"
      units_ij(k)="10^12 W"
      scale_ij(k)=1d-12/dtsrc
      jgrid_ij(k)=2
c
!**** Stability
      k=k+1 !
      IJ_DTDP = k ! DTHETA/DPHI (K S**2/m**2) IN TROPOSPHERE  4 DA
      lname_ij(k) = 'TROP STATIC STABILITY'
      units_ij(k) = 'C/km'
      name_ij(k) = 'dtdz_tropo'
      ia_ij(k) = ia_dga
      scale_ij(k) = 1000.*GRAV*P1000K
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_SSTABX = k ! PEAK DTHETA/DPHI (K S**2/m**2) IN PBL   1 CL
      lname_ij(k) = 'PEAK STATIC STABILITY IN PBL'
      units_ij(k) = 'C/km'
      name_ij(k) = 'dtdzmax_pbl'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.*GRAV*P1000K
      ir_ij(k) = ir_0_18
c
!**** Clouds
      k=k+1 !
      IJ_PCLDL = k ! PCLD(LOW) (1)                            2 RD
      lname_ij(k) = 'LOW LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'pcldl'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PCLDM = k ! PCLD(MID) (1)                            2 RD
      lname_ij(k) = 'MIDDLE LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'pcldm'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PCLDH = k ! PCLD(HIGH) (1)                           2 RD
      lname_ij(k) = 'HIGH LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'pcldh'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PCLDL_SS = k
      lname_ij(k) = 'LARGE-SCALE LOW LEVEL CLOUDINESS'
      units_ij(k) = '%'
      name_ij(k) = 'pcldl_ss'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PSCLD = k ! SHALLOW CONVECTIVE CLOUD COVER  (1)     1 CL
      lname_ij(k) = 'SHALLOW CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'pscld'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_PDCLD = k ! DEEP CONVECTIVE CLOUD COVER     (1)     1 CL
      lname_ij(k) = 'DEEP CONVECTIVE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'pdcld'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_DCNVFRQ = k ! DEEP CONVECTIVE CLOUD OCCURRENCE (1)   1 CL
      lname_ij(k) = 'DEEP CONVECTIVE CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'dcnvfrq'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_SCNVFRQ = k ! SHALLOW CONVECTIVE CLOUD OCCURRENCE (1) 1 CL
      lname_ij(k) = 'SHALLOW CONV CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'scnvfrq'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_CNVFRQ = k
      lname_ij(k) = 'CONV CLOUD FREQUENCY'
      units_ij(k) = '%'
      name_ij(k) = 'cnvfrq'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
C**** Note these two diagnostics are NOT the total cloud cover (that is got by
C**** summing the low+ mid+high diagnostics). Instead, this is the
C**** fraction of time that a cloud appears in the grid box (which may
C**** well cover less than 100% of the box) and the amount of time that
C**** clouds could be detected. TCLDI is needed for weighting the cloud
C**** top pressure and optical depth, and SCLDI is needed for weighting
C**** frequency diags.
c
      k=k+1 !
      IJ_SCLDI = k
      lname_ij(k) = 'FRACTION OF SUNLIT FOR ISCCP CLOUD'
      units_ij(k) = '%'
      name_ij(k) = 'pclds_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_LCLDI = k
      lname_ij(k) = 'LOW LEVEL CLOUDINESS (ISCCP)'
      units_ij(k) = '%'
      name_ij(k) = 'pcldl_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_SCLDI
c
      k=k+1 !
      IJ_MCLDI = k
      lname_ij(k) =
     *     'MIDDLE LEVEL CLOUDINESS (ISCCP)'
      units_ij(k) = '%'
      name_ij(k) = 'pcldm_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_SCLDI
c
      k=k+1 !
      IJ_HCLDI = k
      lname_ij(k) = 'HIGH LEVEL CLOUDINESS (ISCCP)'
      units_ij(k) = '%'
      name_ij(k) = 'pcldh_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_SCLDI
c
      k=k+1 !
      IJ_TAUI = k
      lname_ij(k) = 'CLOUD OPTICAL DEPTH (ISCCP)'
      units_ij(k) = ''
      name_ij(k) = 'optd_isccp'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = IJ_TCLDI

      k=k+1 !
      IJ_WTRCLD = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'WATER CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'wtrcld'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_ICECLD = k ! PCLD (1)  (COMPOSITE OVER ATMOSPHERE)   2 RD
      lname_ij(k) = 'ICE CLOUD COVER'
      units_ij(k) = '%'
      name_ij(k) = 'icecld'
      ia_ij(k) = ia_rad
      scale_ij(k) = 100.
c
      k=k+1 !
      IJ_OPTDW = k
      lname_ij(k) = 'WATER CLOUD OPTICAL DEPTH'
      units_ij(k) = ''
      name_ij(k) = 'optdw'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      denom_ij(k) = IJ_WTRCLD
c
      k=k+1 !
      IJ_OPTDI = k
      lname_ij(k) = 'ICE CLOUD OPTICAL DEPTH'
      units_ij(k) = ''
      name_ij(k) = 'optdi'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      denom_ij(k) = IJ_ICECLD
c
!**** Vegetation and Ground Hydrology variables
      k=k+1
      IJ_GCONATM = k
      name_ij(k) = 'cond_atm' !
      lname_ij(k) = 'CONDUCTANCE OF ATMOSPHERE'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      ir_ij(k) = ir_0_71
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_GCONCAN = k
      name_ij(k) = 'cond_can' !
      lname_ij(k) = 'CONDUCTANCE OF CANOPY'
      units_ij(k) = '.01 m/s'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      ir_ij(k) = ir_0_71
      denom_ij(k) = IJ_VSFR
c
      k=k+1 ! nyk 4/25/03
      IJ_GPP = k    !kg[C]/m2/s
      lname_ij(k) = 'GROSS PRIMARY PRODUCTIVITY'
      !from kg[C]/m2/s, typical range to 30 gC/m2/year
      units_ij(k) = 'g[C]/m2/day'
      name_ij(k) = 'gpp'
      !Scale for mg/m2/day
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY*1000./DTsrc    !scale from kg/s to g/day
      denom_ij(k) = IJ_PSOIL
c     iw  built-in
C NADINE
      k=k+1 !
      IJ_IPP = k    !kg[C]/m2/s
      lname_ij(k) = 'ISOPRENE EMISSION FROM VEG'
      !from kg[C]/m2/s, typical range to 30 gC/m2/year
      units_ij(k) = 'g[C]/m2/day'
      name_ij(k) = 'ipp'
      !Scale for mg/m2/day
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY*1000./DTsrc    !scale from kg/s to g/day
c      iw_ij(k) = iw_soil     !Weight over land
c     iw  built-in
c
      k=k+1 ! nyk 1/10/08
      IJ_RAUTO = k    !kg[C]/m2/s original units
      lname_ij(k) = 'AUTOTROPHIC RESPIRATION'
      units_ij(k) = 'g[C]/m2/day'
      name_ij(k) = 'rauto'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY*1000./DTsrc    !scale from kg/s to g/day
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 ! nyk 1/10/08
      IJ_CLAB = k    !kg[C]/m2
      lname_ij(k) = 'PLANT LABILE CARBON'
      units_ij(k) = 'kg[C]/m2'
      name_ij(k) = 'C_lab'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !PK 1/11/08
      IJ_SOILRESP = k    !kg[C]/m2/s in Ent
      lname_ij(k) = 'SOIL RESPIRATION'
      units_ij(k) = 'g[C]/m2/day'
      name_ij(k) = 'soilresp'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY*1000./DTsrc    !scale from kg/m2/s to g/m2/d
      denom_ij(k) = IJ_PSOIL
c
      k=k+1 !PK 1/11/08
      IJ_SOILCPOOLSUM = k    !g[C]/m2 in Ent
      lname_ij(k) = 'SOIL ORGANIC CARBON POOL'
      units_ij(k) = 'kg[C]/m2'
      name_ij(k) = 'soilCpool'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d-3    !scale from g/m2 to kg/m2
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_LANDCARBON = k
      lname_ij(k) = 'TOTAL LAND ORGANIC CARBON (incl. excess)'
      units_ij(k) = 'kg[C]/m2'
      name_ij(k) = 'landCtot'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0   !the scale is already in kg/m2
      denom_ij(k) = IJ_PSOIL
c
      k=k+1
      IJ_ECVF = k
      lname_ij(k) = 'EXCESS C FLUX DUE TO VEG FRACTIONS CHANGE'
      !from kg[C]/m2/s, typical range to 30 gC/m2/year
      units_ij(k) = 'g[C]/m2/day'
      name_ij(k) = 'ecvf'
      !Scale for mg/m2/day
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY*1000./DTsrc    !scale from kg/s to g/day
      denom_ij(k) = IJ_PSOIL
c
#ifdef ENT_DEBUG_DIAGS
      ij_ent_debug = k+1
      do ent_k1=1,16+11
      do ent_k2=1,16
        write(ent_s1,'(i3.3)') ent_k1
        write(ent_s2,'(i3.3)') ent_k2
      k=k+1 ! nyk 1/10/08
      !IJ_RAUTO = k    !kg[C]/m2/s original units
      lname_ij(k) = 'Ent diag '//ent_s1//ent_s2
      units_ij(k) = 'g[C]/m2/day'
      name_ij(k) = 'ra'//ent_s1//ent_s2
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY*1000./DTsrc    !scale from kg/s to g/day
      denom_ij(k) = IJ_PSOIL
      enddo
      enddo
#endif
c
      k=k+1 ! nyk 5/12/03
      IJ_DLEAF = k    !kg[C]/m2, IJ_DLEAF is accumulated daily.
!      lname_ij(k) = 'LEAF MASS CHANGE'
!      units_ij(k) = 'mg[C]/m2/day'
!      name_ij(k) = 'DLEAF'
!      ia_ij(k) = ia_12hr      !Accumulated once daily, put 2.* in scale
!      scale_ij(k) = 2.*1000000.   !Scale for daily, scale from kg to mg
      lname_ij(k) = 'ANNUAL LEAF MAX GROWTH'
      units_ij(k) = 'g[C]/m2/yr'
      name_ij(k) = 'dleaf'
      ia_ij(k) = ia_inst      !Accumulat instantaneous value for max-min
      scale_ij(k) = 1000.    !Scale from kg to g
      denom_ij(k) = IJ_PSOIL
c     iw  built-in
c
      k=k+1 !YKIM 10/20/09
      IJ_LAI = k    !m2/m2 in Ent
      lname_ij(k) = 'LEAF AREA INDEX'
      units_ij(k) = 'm2/m2'
      name_ij(k) = 'LAI'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1
      IJ_GBETAT = k
      name_ij(k) = 'beta_trans' !
      lname_ij(k) = 'TRANSPIRATION EFFICIENCY, BETAT'
      units_ij(k) = '%'
      ia_ij(k) = ia_src
      scale_ij(k) = 100.
      denom_ij(k) = IJ_VSFR
c
#ifdef CLD_AER_CDNC
      k=k+1
      IJ_dzwm = k
      lname_ij(k) = 'Warm Moist Cnv Cld DZ'
      units_ij(k) = 'm'
      name_ij(k) = 'dzwm'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1
      IJ_dzim = k
      lname_ij(k) = 'Cold Moist Cnv cld DZ'
      units_ij(k) = 'm'
      name_ij(k) = 'dzim'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1
      IJ_dzws = k
      lname_ij(k) = 'Warm Large-scale cld DZ'
      units_ij(k) = 'm'
      name_ij(k) = 'dzws'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1
      IJ_dzis = k
      lname_ij(k) = 'Cold Large-scale cld DZ'
      units_ij(k) = 'm'
      name_ij(k) = 'dzis'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
c
      k=k+1
      IJ_3dNWM = k
      lname_ij(k) = '2D Warm Moist Cnv CDNC '
      units_ij(k) = 'cm^-3'
      name_ij(k) = 'Nwm3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzwm
c
      k=k+1
      IJ_3dNIM = k
      lname_ij(k) = '2D Cold Moist Cnv CDNC '
      units_ij(k) = 'cm^-3'
      name_ij(k) = 'Nim3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzim
c
      k=k+1
      IJ_3dRWM = k
      lname_ij(k) = '2D Warm Moist Conv Reff '
      units_ij(k) = 'um'
      name_ij(k) = 'Rwm3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzwm
c
      k=k+1
      IJ_3dRIM = k
      lname_ij(k) = '2D Cold Moist Conv Reff '
      units_ij(k) = 'um'
      name_ij(k) = 'Rim3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzim
c
      k=k+1
      IJ_3dLWM = k
      lname_ij(k) = '2D Warm Moist Conv LWC  '
      units_ij(k) = 'g m-3'
      name_ij(k) = 'Lwm3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzwm
c
      k=k+1
      IJ_3dLIM = k
      lname_ij(k) = '2D Cold Moist Conv LWC  '
      units_ij(k) = 'g m-3'
      name_ij(k) = 'Lim3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzim
c
      k=k+1
      IJ_3dNWS = k
      lname_ij(k) = '2D Warm Large-scale CDNC '
      units_ij(k) = 'cm^-3'
      name_ij(k) = 'Nws3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzws
c
      k=k+1
      IJ_3dNIS = k
      lname_ij(k) = '2D Cold Large-scale CDNC '
      units_ij(k) = 'cm^-3'
      name_ij(k) = 'Nis3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzis
c
      k=k+1
      IJ_3dRWS = k
      lname_ij(k) = '2D Warm Large-scale Reff '
      units_ij(k) = 'um'
      name_ij(k) = 'Rws3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzws
c
      k=k+1
      IJ_3dRIS = k
      lname_ij(k) = '2D Cold Large-scale Reff '
      units_ij(k) = 'um'
      name_ij(k) = 'Ris3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzis
c
      k=k+1
      IJ_3dLWS = k
      lname_ij(k) = '2D Warm Large-scale LWC '
      units_ij(k) = 'g m-3'
      name_ij(k) = 'Lws3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzws
c
      k=k+1
      IJ_3dLIS = k
      lname_ij(k) = '2D Cold Large-scale LWC '
      units_ij(k) = 'g m-3'
      name_ij(k) = 'Lis3d'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.
      denom_ij(k) = ij_dzis
#endif
c
#ifdef CHL_DIAGNOSTIC
      k=k+1
      IJ_CHL = k
      lname_ij(k) = 'Total Chlorophyll'
      units_ij(k) = 'mg/m^3'
      scale_ij(k) = 1.
      name_ij(k) = 'chl'
      ia_ij(k) = ia_rad
      scale_ij(k) = 1.
      denom_ij(k) = IJ_POCEAN
#endif /* #ifdef CHL_DIAGNOSTIC */
c

#ifdef TRACERS_ON
      allocate(ij_kw(gasex_index%getsize()))
      allocate(ij_alpha(gasex_index%getsize()))
      allocate(ij_gasx(gasex_index%getsize()))
      do ngx=1, gasex_index%getsize()
        n=gasex_index%at(ngx)
        k=k+1
        IJ_Kw(ngx) = k
        lname_ij(k) = 'Transfer Velocity '//trname(n)
        units_ij(k) = 'm/s'
        scale_ij(k) = 1.
        name_ij(k) = 'Kw_gas_'//trname(n)
        ia_ij(k) = ia_srf
        scale_ij(k) = 1.
        denom_ij(k) = IJ_POCEAN

        k=k+1
        IJ_alpha(ngx) = k
        lname_ij(k) = 'Solubility '//trname(n)
        units_ij(k) = 'mol/m3/uatm'
        scale_ij(k) = 1.
        name_ij(k) = 'alpha_gas_'//trname(n)
        ia_ij(k) = ia_srf
        scale_ij(k) = 1.
        denom_ij(k) = IJ_POCEAN

        k=k+1
        IJ_gasx(ngx) = k
        lname_ij(k) = 'Gas Exchange Flux '//trname(n)
        units_ij(k) = 'mol/m2/yr'
        scale_ij(k) = 1.
        name_ij(k) = 'trgasex_'//trname(n)
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
        denom_ij(k) = IJ_POCEAN
      end do
#endif

#ifdef TRACERS_OBIO_RIVERS
      k=k+1
      IJ_rvrflo = k  ! Mass Outflow by Rivers without gmelt (10**5 kg/s)  E-5/DTS*1 RV
      lname_ij(k) = 'Mass Outflow by Rivers without gmelt'
      units_ij(k) = '10^5 kg/s'
!      scale_ij(k) = 1.
      scale_ij(k) = 1.d-5/DTsrc
      name_ij(k) = 'FLOWO'
      ia_ij(k) = ia_src
!      denom_ij(k) = IJ_POCEAN

#endif
c
#ifdef CLD_AER_CDNC
c
      k=k+1 !
      IJ_WISUM = k ! ICE WATER PATH (kg/m**2)             1 CL
      lname_ij(k) = 'ICE WATER PATH'
      units_ij(k) = '.1 kg/m^2'
      name_ij(k) = 'iwp'
      ia_ij(k) = ia_src
      scale_ij(k) = 10.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_WMCLWP = k ! MC LIQUID WATER PATH (kg/m**2)             1 CL
      lname_ij(k) = 'MC LIQUID WATER PATH'
      units_ij(k) = '.1 kg/m^2'
      name_ij(k) = 'mclwp'
      ia_ij(k) = ia_src
      scale_ij(k) = 10.
      ir_ij(k) = ir_0_18
c
      k=k+1 !
      IJ_WMCTWP = k ! MC Total WATER PATH (kg/m**2)             1 CL
      lname_ij(k) = 'MC TOTAL WATER PATH'
      units_ij(k) = '.1 kg/m^2'
      name_ij(k) = 'mctwp'
      ia_ij(k) = ia_src
      scale_ij(k) = 10.
      ir_ij(k) = ir_0_18
#endif
c
c     k=k+1 !
cfree IJ_EMTMOM = k ! INCIDENT MTN EAST MOM. FLUX (MB-M/S**2)  1 SD
cfree lname_ij(k) = 'INCIDENT MTN EAST MOMENTUM FLUX'
cfree units_ij(k) = 'mb m/s^2'
cfree name_ij(k) = 'EMTMOM'
cfree ia_ij(k) = ia_src
cfree scale_ij(k) = 1./DTsrc
c
c     k=k+1 !
cfree IJ_SMTMOM = k ! INCIDENT MTN SOUTH MOM. FLUX (MB-M/S**2) 1 SD
cfree lname_ij(k) = 'INCIDENT MTN SOUTH MOMENTUM FLUX'
cfree units_ij(k) = 'mb m/s^2'
cfree name_ij(k) = 'SMTMOM'
cfree ia_ij(k) = ia_src
cfree scale_ij(k) = 1./DTsrc
c
      k=k+1 !
      IJ_LKICE = k
      lname_ij(k) = 'LAKE ICE WEIGHTING' ! for lake freeze/thaw diags
      units_ij(k) = ' '
      name_ij(k) = 'LKICEWT'
      ia_ij(k) = ia_inst
c
      k=k+1 !
      IJ_LKON = k
      lname_ij(k) = 'LAST ICE-FREE DAY (SH-58,NH-242)'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'lkonday'
      ia_ij(k) = ia_inst
      denom_ij(k) = IJ_LKICE
c
      k=k+1 !
      IJ_LKOFF = k
      lname_ij(k) = 'LAST ICED-UP DAY (SH-58,NH-242)'
      units_ij(k) = 'JULIAN DAY'
      name_ij(k) = 'lkoffday'
      ia_ij(k) = ia_inst
      denom_ij(k) = IJ_LKICE
c
#ifdef IRRIGATION_ON
      k=k+1
      IJ_IRRW_TOT = k
      lname_ij(k) = 'POTENTIAL IRRIGATION'
      units_ij(k) = 'mm/d'
      name_ij(k) = 'irrig_w_tot'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0 ! per total grid area; =IJ_PSOIL for per earth
c
      k=k+1
      IJ_IRRW = k
      lname_ij(k) = 'TOTAL (ACTUAL) IRRIGATION WATER'
      units_ij(k) = 'mm/d'
      name_ij(k) = 'irrig_w'
      ia_ij(k) = ia_src
      scale_ij(k) = SECONDS_PER_DAY/DTsrc
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      IJ_IRRE = k
      lname_ij(k) = 'HEAT OF IRRIGATION'
      units_ij(k) = 'W/m2'
      name_ij(k) = 'irrig_e'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0/DTsrc
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      IJ_GMLir = k
      lname_ij(k) = 'ENTHALPY OF LAKE/RIVER LOST TO IRRIGATION'
      units_ij(k) = 'J'
      name_ij(k) = 'gml_irrigate'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      IJ_irrgw = k
      lname_ij(k) = 'IRRIGATION WATER FROM EXTERNAL SOURCE (GRNDWATER)'
      units_ij(k) = 'mm/d'
      name_ij(k) = 'irrig_gw'
      ia_ij(k) = ia_src
      scale_ij(k) = 1000.d0*SECONDS_PER_DAY
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      IJ_irrgwE = k
      lname_ij(k) = 'HEAT OF EXTERNALLY ADDED IRRIGATION'
      units_ij(k) = 'W/m2'
      name_ij(k) = 'irrig_gwE'
      ia_ij(k) = ia_src
      scale_ij(k) = 1.d0
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
c
      k=k+1
      IJ_MWLir = k
      lname_ij(k) = 'MASS OF LAKE/RIVER WATER USED FOR IRRIGATION'
      units_ij(k) = 'kg/s'
      name_ij(k) = 'mwl_irrigate'
      ia_ij(k) = ia_src
      scale_ij(k) = 1 / DTSRC
      !ir_ij(k) = ir_m1_3
      denom_ij(k) = 0
#endif
c
c Gravity Wave diagnostics
      iDO_GWDRAG = 0
      if (DO_GWDRAG) then
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW1 = k
      name_ij(k) = 'ij_def_drag_mom_flux'
      lname_ij(k) = 'DEFORM. DRAG MOM FLUX'
      units_ij(k) = '.1 N/m**2'
      scale_ij(k) = 10.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW2 = k
      name_ij(k) = 'ij_mtn_wave_mom_flux'
      lname_ij(k) = 'MTN WAVE MOM. FLUX'
      units_ij(k) = '.1 N/m**2'  ! dynes/cm^2
      scale_ij(k) = 10.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW3 = k
      name_ij(k) = 'ij_shr_wave_mom_flux'
      lname_ij(k) = 'SHEAR WAVE MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW4 = k
      name_ij(k) = 'ij_mc_c_m10r_mom_flux'
      lname_ij(k) = 'MC C=-10R MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW5 = k
      name_ij(k) = 'ij_mc_c_m20r_mom_flux'
      lname_ij(k) = 'MC C=-20R MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW6 = k
      name_ij(k) = 'ij_mc_c_m40r_mom_flux'
      lname_ij(k) = 'MC C=-40R MOM. FLUX'
      units_ij(k) = '.001 N/m**2'
      scale_ij(k) = 1000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW7 = k
      name_ij(k) = 'ij_phase_speed_of_wind_shear'
      lname_ij(k) = 'PHASE SPEED OF SHEAR WAVE'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m45_130
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW8 = k
      name_ij(k) = 'ij_source_speed_of_mc'
      lname_ij(k) = 'MC SOURCE WIND SPEED'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m45_130
c
      k=k+1                ! ij diags from gwdrag calculations
      IJ_GW9 = k
      name_ij(k) = 'ij_exit_tot_mom_flux'
      lname_ij(k) = 'EXIT TOT. MOM. FLUX'
      units_ij(k) = '.0001 N/m**2'
      scale_ij(k) = 10000.*100.*BYGRAV
      ia_ij(k) = ia_src
      ir_ij(k) = ir_m1_3
      iDO_GWDRAG = k-IJ_GW1+1
      END IF
c
      k=k+1 !
      IJ_DSKIN   = k !
      lname_ij(k) = 'SKIN TEMPERATURE OFFSET (OCEAN)'
      units_ij(k) = '0.1 C'
      name_ij(k) = 'dskin'
      ia_ij(k) = ia_srf
      scale_ij(k) = 10.
      ir_ij(k) = ir_m9_26
c
      k=k+1 !
      IJ_DSKINSNOW   = k !
      lname_ij(k) = 'SKIN TEMPERATURE OFFSET (SNOW ON OC/LK ICE)'
      units_ij(k) = '0.1 C'
      name_ij(k) = 'dskinsnow'
      ia_ij(k) = ia_srf
      scale_ij(k) = 10.
      ir_ij(k) = ir_m9_26
c
      if (calc_wspdf == 1) then
        k=k+1
        ij_wspdf = k
        lname_ij(k) = 'PDF MEAN SURFACE WIND SPEED'
        name_ij(k) = 'wspdf'  
        units_ij(k) = 'm/s'
        ia_ij(k) = ia_srf
      endif
                                    
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      k=k+1
      ij_wsgcm = k
      lname_ij(k) = 'GCM SURFACE WIND SPEED'
      name_ij(k) = 'wsgcm'
      units_ij(k) = 'm/s'
      ia_ij(k)= ia_srf
      k=k+1
      IJ_wdry = k
      lname_ij(k) = 'DRY CONVECTIVE VELOCITY SCALE'
      name_ij(k) = 'wsubwd'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_srf
      k=k+1
      IJ_wtke = k
      lname_ij(k) = 'TKE VELOCITY SCALE'
      name_ij(k) = 'wsubtke'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_srf
      k=k+1
      IJ_wmoist = k
      lname_ij(k) = 'MOIST CONVECTIVE VELOCITY SCALE'
      name_ij(k) = 'wsubwm'
      units_ij(k) = 'm/s'
      ia_ij(k) = ia_srf
#endif
#ifdef CALCULATE_FLAMMABILITY
      k=k+1        !rh for flammability 
        ij_flam_rh = k
        lname_ij(k) = 'rh for FLAMMABILITY'
        units_ij(k) = 'none'
        name_ij(k) = 'FLAMM_rh'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.e0
      k=k+1        ! prec for flammability 
        ij_flam_prec = k
        lname_ij(k) = 'prec for FLAMMABILITY'
        units_ij(k) = 'mm/day'
        name_ij(k) = 'FLAMM_prec'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.e0
      k=k+1        ! tsurf for flammability 
        ij_flam_tsurf = k
        lname_ij(k) = 'tsurf for FLAMMABILITY'
        units_ij(k) = 'K'
        name_ij(k) = 'FLAMM_tsurf'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.e0
      k=k+1        ! flammability of vegetation
        ij_flam = k
        lname_ij(k) = 'VEGETATION FLAMMABILITY'
        units_ij(k) = 'none'
        name_ij(k) = 'FLAMM'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.e0
#if (defined DYNAMIC_BIOMASS_BURNING)&&(defined DETAILED_FIRE_OUTPUT)
      k=k+1   ! vegetation fractions used in fire model
        ij_flamV = k
        lname_ij(k) = 'FLAMMABLE VEG FRACTION'
        units_ij(k) = 'fraction of whole grid'
        name_ij(k) = 'FVFRAC'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
#endif /* DYNAMIC_BIOMASS_BURNING && DETAILED_FIRE_OUTPUT */
      k=k+1        ! vegetation density for fire model purposes
        ij_fvden = k
        lname_ij(k) = 'FIRE MODEL VEGETATION DENSITY'
        units_ij(k) = 'none'
        name_ij(k) = 'FVDEN'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
#ifdef ANTHROPOGENIC_FIRE_MODEL
      k=k+1        ! frac dynamic BB emis from non suppression
        ij_nsuppress = k
        lname_ij(k) = 'FRAC OF DYN BB EMIS DUE TO NONSUPPRESS'
        units_ij(k) = 'none'
        name_ij(k) = 'f_nsuppress'
        ia_ij(k) = ia_inst
        scale_ij(k) = 1.
      k=k+1        ! frac dynamic biomass burning emis from lightning 
        ij_cgign = k
        lname_ij(k) = 'FRAC OF DYN BB EMIS DUE TO CG LIGT IGN. ONLY'
        units_ij(k) = 'm-2 s-1'
        name_ij(k) = 'f_ignCG'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
      k=k+1        ! frac dynamic biomass burning emis from humans
        ij_humanign = k
        lname_ij(k) = 'FRAC OF DYN BB EMIS DUE TO HUMAN IGN. ONLY'
        units_ij(k) = 'm-2 s-1'
        name_ij(k) = 'f_ignHUMAN'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
#endif /* ANTHROPOGENIC_FIRE_MODEL */
      k=k+1        ! The Fire Count (no need to save for ubiquitous
        ij_fireC = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'FIRE COUNT FOR DYN BIOBURN USING ANTHRO MODEL'
        units_ij(k) = 'm-2 s-1'
        name_ij(k) = 'fireCount'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
      k=k+1        !gridbox RH used for burnt area calculation 
        ij_barh1 = k 
        lname_ij(k) = 'RH for BA calc'
        units_ij(k) = ''
        name_ij(k) = 'BA_RH'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
      k=k+1        !gridbox wsurf used for burnt area calculation 
        ij_bawsurf = k 
        lname_ij(k) = 'wsurf for BA calc'
        units_ij(k) = 'm/s'
        name_ij(k) = 'BA_wsurf'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
      k=k+1        !gridbox burnt area of tree surface type 
        ij_ba_tree = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'TREE BA FOR DYN BIOBURN USING ANTHRO MODEL'
        units_ij(k) = 'm2'
        name_ij(k) = 'BA_tree'
        ia_ij(k) = ia_inst
        scale_ij(k) = 1.
      k=k+1        !gridbox burnt area of shrub surface type 
        ij_ba_shrub = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'SHRUB BA FOR DYN BIOBURN USING ANTHRO MODEL'
        units_ij(k) = 'm2'
        name_ij(k) = 'BA_shrub'
        ia_ij(k) = ia_inst
        scale_ij(k) = 1.
      k=k+1        !gridbox burnt area of grass surface type 
        ij_ba_grass = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'GRASS BA FOR DYN BIOBURN USING ANTHRO MODEL'
        units_ij(k) = 'm2'
        name_ij(k) = 'BA_grass'
        ia_ij(k) = ia_inst
        scale_ij(k) = 1.
      k=k+1        !gridbox fire spread area of tree surface type 
        ij_a_tree = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'FIRE SPREAD AREA FOR TREES'
        units_ij(k) = 'm2'
        name_ij(k) = 'A_tree'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
      k=k+1        !gridbox fire spread area of shrub surface type 
        ij_a_shrub = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'FIRE SPREAD AREA FOR SHRUBS'
        units_ij(k) = 'm2'
        name_ij(k) = 'A_shrub'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
      k=k+1        !gridbox fire spread area of grass surface type 
        ij_a_grass = k ! case, since it is constant factor times flammability)
        lname_ij(k) = 'FIRE SPREAD AREA FOR GRASS'
        units_ij(k) = 'm2'
        name_ij(k) = 'A_grass'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.
#endif /* CALCULATE_FLAMMABILITY */
#if(defined CALCULATE_LIGHTNING)||(defined TRACERS_SPECIAL_Shindell)
      k=k+1        ! lightning flash rate
        ij_flash = k
        lname_ij(k) = 'LIGHTNING FLASH DENSITY'
        units_ij(k) = '1.e-10 flashes/m2/s'
        name_ij(k) = 'FLASH'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.e10/DTsrc
      k=k+1        ! lightning cloud-to-ground flash rate
        ij_CtoG = k
        lname_ij(k) = 'LIGHTNING CLOUD TO GROUND FLASH DENSITY'
        units_ij(k) = '1.e-10 flashes/m2/s'
        name_ij(k) = 'CtoG'
        ia_ij(k) = ia_src
        scale_ij(k) = 1.e10/DTsrc
#endif
#ifdef ACCMIP_LIKE_DIAGS
      k=k+1
        ij_fcghg(1,1) = k
        lname_ij(k) = 'SW TOA CH4 RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swch4_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(2,1) = k
        lname_ij(k) = 'LW TOA CH4 RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwch4_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(1,2) = k
        lname_ij(k) = 'SW TOA N2O RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swn2o_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(2,2) = k
        lname_ij(k) = 'LW TOA N2O RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwn2o_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(1,3) = k
        lname_ij(k) = 'SW TOA CFC11 RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swc11_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(2,3) = k
        lname_ij(k) = 'LW TOA CFC11 RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwc11_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(1,4) = k
        lname_ij(k) = 'SW TOA CFC12 RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'swc12_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
      k=k+1
        ij_fcghg(2,4) = k
        lname_ij(k) = 'LW TOA CFC12 RADIATIVE FORCING'
        units_ij(k) = 'W/m^2'
        name_ij(k) = 'lwc12_toa'
        ia_ij(k) = ia_rad_frc
        scale_ij(k) = 1.
#endif /* ACCMIP_LIKE_DIAGS */

c
c the following are not accumulated
c

      k = k + 1

      ij_topo = k
      name_ij(k) = 'topog'
      lname_ij(k) = 'TOPOGRAPHY'
      units_ij(k) = 'METERS'
      ir_ij(k) = ir_0_3550
      ia_ij(k) = ia_inst

      k = k + 1
      ij_fland = k
      name_ij(k) = 'frac_land'
      lname_ij(k) = 'LAND COVERAGE'
      units_ij(k) = '%'
      scale_ij(k) = 100.

      k = k + 1
      ij_jet = k
      name_ij(k) = 'jet_speed'
      lname_ij(k) = 'JET SPEED'
      units_ij(k) = 'm/s'
      scale_ij(k) = scale_ij(ij_ujet)
      ia_ij(k) = ia_ij(ij_ujet)
      ir_ij(k) = ir_0_71

      k = k + 1
      ij_wsmn = k
      name_ij(k) = 'rt_usmn2_vsmn2'
      lname_ij(k) = 'SURF WIND SPEED FROM Uav,Vav'
      units_ij(k) = 'm/s'
      scale_ij(k) = scale_ij(ij_us)
      ia_ij(k) = ia_ij(ij_us)
      ir_ij(k) = ir_0_18

      k = k + 1
      ij_jetdir = k
      name_ij(k) = 'jet_dir'
      lname_ij(k) = 'JET DIRECTION'
      units_ij(k) = 'CW NOR'
      scale_ij(k) = 360./twopi
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_angl

      k = k + 1
      ij_wsdir = k
      name_ij(k) = 'srf_wind_dir'
      lname_ij(k) = 'SURFACE WIND DIRECTION'
      units_ij(k) = 'CW NOR'
      scale_ij(k) = 360./twopi
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_angl

      k = k + 1
      ij_netrdp = k
      name_ij(k) = 'net_rad_planet'
      lname_ij(k) = 'NET RAD. OF PLANET'
      units_ij(k) = 'W/m^2'
      ir_ij(k) = ir_m530_190
      ia_ij(k) = ia_ij(ij_trnfp0)

      k = k + 1
      ij_albp = k
      name_ij(k) = 'plan_alb'
      lname_ij(k) = 'PLANETARY ALBEDO'
      units_ij(k) = '%'
      scale_ij(k) = 100.
      ia_ij(k) = ia_ij(ij_srincp0)
      denom_ij(k) = ij_srincp0

      k = k + 1
      ij_albg = k
      name_ij(k) = 'grnd_alb'
      lname_ij(k) = 'GROUND ALBEDO'
      units_ij(k) = '%'
      scale_ij(k) = 100.
      ia_ij(k) = ia_ij(ij_srincg)
      denom_ij(k) = ij_srincg

      k = k + 1
      ij_albv = k
      name_ij(k) = 'vis_alb'
      lname_ij(k) = 'VISUAL ALBEDO'
      units_ij(k) = '%'
      scale_ij(k) = 100.
      ia_ij(k) = ia_ij(ij_srref)
      denom_ij(k) = ij_srincp0

      k = k + 1
      ij_albgv = k
      name_ij(k) = 'grnd_alb_vis'
      lname_ij(k) = 'GROUND ALBEDO IN VISUAL RANGE'
      units_ij(k) = '%'
      scale_ij(k) = 100.
      ia_ij(k) = ia_ij(ij_srvis)
      denom_ij(k) = ij_srincp0

      k = k + 1
      ij_ntdsese = k
      name_ij(k) = 'stand_eddy_nt_dse'
      lname_ij(k) = 'NT DRY STAT ENR BY ST ED' ! NORTHWD TRANSP
      units_ij(k) = 'E14 WT'

      k = k + 1
      ij_ntdsete = k
      name_ij(k) = 'trans_eddy_nt_dse'
      lname_ij(k) = 'NT DRY STAT ENR BY TR ED' ! NORTHWD TRANSP
      units_ij(k) = 'E14 WT'

      k = k + 1
      ij_grow = k
      name_ij(k) = 'grow_seas'
      lname_ij(k) = 'GROWING SEASON'
      units_ij(k) = 'days'
      ir_ij(k) = ir_0_180
      ia_ij(k) = ia_inst

C**** Also include MSU radiation diagnostics here

      k=k+1 !
      ij_msutlt = k
      name_ij(k) = 'Tmsu-TLT'
      lname_ij(k) = 'MSU-TLT TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_m80_28

      k = k + 1
      ij_msutmt = k
      name_ij(k) = 'Tmsu_TMT'
      lname_ij(k) = 'MSU-TMT TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_m80_28

      k = k + 1
      ij_msutls = k
      name_ij(k) = 'Tmsu_TLS'
      lname_ij(k) = 'MSU-TLS TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_m80_28

C**** Add in Stratospheric Sounding Units (3 channels)
      k = k + 1
      ij_ssu1 = k
      name_ij(k) = 'Tssu_ch1'
      lname_ij(k) = 'SSU-Ch 1 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_m80_28

      k = k + 1
      ij_ssu2 = k
      name_ij(k) = 'Tssu_ch2'
      lname_ij(k) = 'SSU-Ch 2 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_m80_28

      k = k + 1
      ij_ssu3 = k
      name_ij(k) = 'Tssu_ch3'
      lname_ij(k) = 'SSU-Ch 3 TEMPERATURE'
      units_ij(k) = 'C'
      ia_ij(k) = ia_inst
      ir_ij(k) = ir_m80_28

C**** Also include Land-Ocean Temperature Index
      k = k + 1
      ij_LOTI = k
      name_ij(k) = 'L-O_TI'
      lname_ij(k) = 'Land-Ocean TEMPERATURE Index'
      units_ij(k) = 'C'
      ia_ij(k) = ia_srf
      ir_ij(k) = ir_m80_28

C****
      k = k + 1
      ij_Tatm = k
      name_ij(k) = 'Tatm'
      lname_ij(k) = 'ATMOSPHERIC TEMPERATURE'
      units_ij(k) = 'C'
      denom_ij(k) = IJ_PRES
      ia_ij(k) = ia_dga ! IA_IJL(IJK_TX)

      k=k+1
      IJ_TMINC = k ! MIN composite TS FOR CURRENT DAY (C)
      lname_ij(k) = 'SURFACE AIR TEMPERATURE: DIURNAL LOW' ! composite
      units_ij(k) = 'C'
      name_ij(k) = 'TMINC'
      ia_ij(k) = ia_12hr   ! really ia_24hr
      scale_ij(k) = 2.*1.  ! really 1.
      ir_ij(k) = ir_m80_28
      denom_ij(k) = 0

#ifdef HEALY_LM_DIAGS
      k=k+1 !
      IJ_CROPS  = k ! CROPS (%)                                1 GD
      lname_ij(k) = 'CROP COVER'
      units_ij(k) = '%'
      name_ij(k) = 'crops'
      ia_ij(k) = ia_srf
      scale_ij(k) = 1.
      ir_ij(k) = ir_pct
#endif


      if (AM_I_ROOT()) then
         if (k .gt. kaij) then
            write (6,*) 'ij_defs: Increase kaij=',kaij,' to at least ',k
            call stop_model( 'kaij too small', 255 )
         end if

         write (6,*) 'Number of AIJ diagnostics defined: kaijmax=',k
         if(qcheck) then
           do kk=1,k
             write (6,'(i4,'':'',a)') kk,trim(lname_ij(kk))
           end do
         endif
      end if

c
c Define min/max diagnostics
c
      k=0
c
      k=k+1
      IJ_TSURFmin = k
      lname_ijmm(k) = 'Minimum surface air temperature'
      units_ijmm(k) = 'C'
      name_ijmm(k) = 'tsurf_min'
      scale_ijmm(k) = -1.  ! flip max -> min
c
      k=k+1
      IJ_TSURFmax = k
      lname_ijmm(k) = 'Maximum surface air temperature'
      units_ijmm(k) = 'C'
      name_ijmm(k) = 'tsurf_max'
      scale_ijmm(k) = 1.


      if (AM_I_ROOT()) then
         if (k .gt. kaijmm) then
            write (6,*) 'ij_defs: Increase kaijmm=',kaijmm,
     &          ' to at least ',k
            call stop_model( 'kaijmm too small', 255 )
         end if
      end if


c
c Declare the dimensions and metadata of AIJ output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      call init_cdl_type('cdl_aij',cdl_ij_template)
#ifndef SCM
      call add_dim(cdl_ij_template,'shnhgm',3)
#endif
#ifdef CUBED_SPHERE
      ijstr='tile,y,x);'
      do i=1,im
        x_dummy(i) = -1d0 + 2d0*(dble(i)-.5d0)/im
      enddo
      call add_coord(cdl_ij_template,'x',im,
     &     long_name='nondimensional cube coordinate',
     &     coordvalues=x_dummy)
      call add_coord(cdl_ij_template,'y',im,
     &     long_name='nondimensional cube coordinate',
     &     coordvalues=x_dummy)
      call add_dim(cdl_ij_template,'tile',6)
      call add_dim(cdl_ij_template,'nv',4)
      call add_var(cdl_ij_template,'float lon('//trim(ijstr),
     &     units='degrees_east')
      call add_varline(cdl_ij_template,'lon:bounds = "lonbds" ;')
      call add_var(cdl_ij_template,'float lat('//trim(ijstr),
     &     units='degrees_north')
      call add_varline(cdl_ij_template,'lat:bounds = "latbds" ;')
      call add_var(cdl_ij_template,'float lonbds(tile,y,x,nv) ;',
     &     units='degrees_east')
      call add_var(cdl_ij_template,'float latbds(tile,y,x,nv) ;',
     &     units='degrees_north')

      call init_cdl_type('cdl_aij_latlon',cdl_ij_latlon_template)
      call add_dim(cdl_ij_latlon_template,'shnhgm',3)
      call add_coord(cdl_ij_latlon_template,'lon',1,
     &     units='degrees_east')
      call add_coord(cdl_ij_latlon_template,'lat',1,
     &     units='degrees_north')
      call add_var(cdl_ij_latlon_template,'float axyp(lat,lon) ;',
     &     units='m^2',long_name='gridcell area')
      cdl_ij_latlon = cdl_ij_latlon_template ! invoke a copy method later
#else
      ijstr='lat,lon);'
      call add_coord(cdl_ij_template,'lon',im,units='degrees_east',
     &     coordvalues=lon_dg(:,1))
      call add_coord(cdl_ij_template,'lat',jm,units='degrees_north',
     &     coordvalues=lat_dg(:,1))
#endif
      call add_var(cdl_ij_template,'float axyp('//trim(ijstr),
     &     units='m^2',long_name='gridcell area')

      cdl_ij = cdl_ij_template ! invoke a copy method later
      cdl_ijmm = cdl_ij_template

      ! erase the grouping information if not applying it
      if(.not.fuse_groups) index1(:) = 0

      do k=1,kaij
        if(trim(units_ij(k)).eq.'unused') cycle
        if(index1(k).eq.0) then
          set_miss = denom_ij(k).ne.0
          varstr='float '//trim(name_ij(k))//'('//trim(ijstr)
          varstrll='float '//trim(name_ij(k))//'(lat,lon);'
          long_name=lname_ij(k)
          auxvar_string=
     &         'float '//trim(name_ij(k))//'_hemis(shnhgm);'
        elseif(index1(k).eq.k) then
          name_ij(k) = name3(k)
          nq = count(index1==k)
          if(count(coord3(:,k).ne.-1d30).eq.nq) then
            if(len_trim(dim3units(k)).gt.0) then
              call add_coord(cdl_ij,trim(dim3name(k)),nq,
     &             coordvalues=coord3(1:nq,k),units=trim(dim3units(k)))
#ifdef CUBED_SPHERE
              call add_coord(cdl_ij_latlon,trim(dim3name(k)),nq,
     &             coordvalues=coord3(1:nq,k),units=trim(dim3units(k)))
#endif
            else
              call add_coord(cdl_ij,trim(dim3name(k)),nq,
     &             coordvalues=coord3(1:nq,k))
#ifdef CUBED_SPHERE
              call add_coord(cdl_ij_latlon,trim(dim3name(k)),nq,
     &             coordvalues=coord3(1:nq,k))
#endif
            endif
          elseif(dim3info_index(k).eq.0) then
            call add_dim(cdl_ij,trim(dim3name(k)),nq)
#ifdef CUBED_SPHERE
            call add_dim(cdl_ij_latlon,trim(dim3name(k)),nq)
#endif
          else
            dim3name(k) = dim3name(dim3info_index(k))
          endif
          set_miss = any(denom_ij(k:k+nq-1).ne.0)
          varstr='float '//trim(name3(k))//'('//trim(dim3name(k))//
     &             ','//trim(ijstr)
          varstrll='float '//trim(name3(k))//'('//trim(dim3name(k))//
     &             ',lat,lon);'
          long_name=lname3(k)
          auxvar_string=
     &         'float '//trim(name3(k))//'_hemis('//
     &                trim(dim3name(k))//',shnhgm);'
        else
          cycle
        endif
        call add_var(cdl_ij,
     &       trim(varstr),
     &       units=trim(units_ij(k)),
     &       long_name=trim(long_name),
     &       auxvar_string=trim(auxvar_string),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#ifdef CUBED_SPHERE
        call add_var(cdl_ij_latlon,
     &       trim(varstrll),
     &       units=trim(units_ij(k)),
     &       long_name=trim(long_name),
     &       auxvar_string=trim(auxvar_string),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#endif
      enddo

      do k=1,kaijmm
        if(trim(units_ijmm(k)).eq.'unused') cycle
        call add_var(cdl_ijmm,
     &       'float '//trim(name_ijmm(k))//'('//trim(ijstr),
     &       units=trim(units_ijmm(k)),
     &       long_name=trim(lname_ijmm(k)))
      enddo

      deallocate(index1,name3,dim3name,coord3)
      deallocate(dim3info_index,lname3,dim3units)

      return
      end subroutine ij_defs


      subroutine jl_defs
      use CONSTANT, only : grav,twopi,sha,rgas,bygrav,radius,lhe
     &     ,bymrat
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : dtsrc,qcheck
      use DYNAMICS, only : nidyn,do_gwdrag
      USE ATM_COM, only : pmidl00,lm_req
      use DIAG_COM
      use DIAG_COM_RAD
#ifdef CUBED_SPHERE
      use GCDIAG, only : fim
#endif
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use cdl_mod
      use MDIAG_COM, only : make_timeaxis
      implicit none
      integer :: l,k,kk
      character(len=10) :: zstr,powstr
      real*8, dimension(lm) :: one_to_lm
      logical :: set_miss
c
      do k=1,kajl
         write(sname_jl(k),'(a3,i3.3)') 'AJL',k
         lname_jl(k) = 'unused'
         units_jl(k) = 'unused'
         scale_jl(k) = 1.
         pow_jl(k) = 0
         jgrid_jl(k) = 1
         lgrid_jl(k) = ctr_ml ! default for most qtys
         denom_jl(k) = 0 ! jl_dpa NOT set as default
         ia_jl(k) = ia_src
      enddo
c
      k=0
c
      k=k+1
      jl_dpa = k
      sname_jl(k) = 'jl_dpa'
      lname_jl(k) = 'MASS AT PRIMARY GRID CELLS' ! not printed
      units_jl(k) = 'mb'
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_dpasrc = k
      sname_jl(k) = 'jl_dpasrc'
c      lname_jl(k) = 'MASS AT PRIMARY GRID CELLS (SRC TIME)' ! not printed
c      units_jl(k) = 'kg/m2'
c
      k=k+1
      jl_mcmflx = k
      sname_jl(k) = 'mc_mflx' !                  'FMX(MC)'
      lname_jl(k) = 'VERTICAL MASS EXCHANGE FROM MOIST CONVECTION'
      units_jl(k) = 'mb/s'
      scale_jl(k) = 1./(DTsrc)
      pow_jl(k) = -5
      lgrid_jl(k) = edg_ml
c
      k=k+1
      jl_mcdflx = k
      sname_jl(k) = 'mc_dflx' !                  'DDMX(MC)'
      lname_jl(k) = 'DOWNDRAFT MASS FLUX FROM MOIST CONVECTION'
      units_jl(k) = 'mb/s'
      scale_jl(k) = 1./(DTsrc)
      pow_jl(k) = -5
      lgrid_jl(k) = edg_ml
c
      k=k+1
      jl_srhr = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'srad_heat' !
      lname_jl(k) = 'SOLAR RADIATION HEATING RATE' !'SRHR'
      units_jl(k) = 'K/DAY' !'W/m^2'
      pow_jl(k) = -2
      scale_jl(k) = 1.D-2*GRAV*SECONDS_PER_DAY/SHA
      ia_jl(k) = ia_rad
c
      k=k+1
      !jl_srhr_radonly = k
      sname_jl(k) = 'srad_heat_radonly' !
      lname_jl(k) = 'SOLAR RADIATION HEATING IN RADIATION-ONLY LAYERS'
      units_jl(k) = 'K/DAY'
      pow_jl(k) = -2
      scale_jl(k) = 1.D-2*GRAV*SECONDS_PER_DAY/SHA
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_trcr = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'trad_cool' !
      lname_jl(k) = 'THERMAL RADIATION COOLING RATE' !'TRHR'
      units_jl(k) = 'K/DAY' !'W/m^2'
      pow_jl(k) = -2
      scale_jl(k) = -1.D-2*GRAV*SECONDS_PER_DAY/SHA
      ia_jl(k) = ia_rad
c
      k=k+1
      !jl_trcr_radonly = k
      sname_jl(k) = 'trad_cool_radonly' !
      lname_jl(k) = 'THERMAL RADIATION COOLING IN RADIATION-ONLY LAYERS'
      units_jl(k) = 'K/DAY'
      pow_jl(k) = -2
      scale_jl(k) = -1.D-2*GRAV*SECONDS_PER_DAY/SHA
      ia_jl(k) = ia_rad
c
      k = k + 1
      jl_rad_cool = k ! not accumulated
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'rad_cool'
      lname_jl(k) = 'TOTAL RADIATION COOLING RATE'
      units_jl(k) = 'W/(m^2*mb)'
      scale_jl(k) = -1.
      ia_jl(k) = ia_rad
      pow_jl(k) = -2
c
      k = k + 1
      !jl_rad_cool_radonly = k ! not accumulated
      sname_jl(k) = 'rad_cool_radonly'
      lname_jl(k) = 'TOTAL COOLING RATE, RADIATION-ONLY LAYERS'
      units_jl(k) = 'W/(m^2*mb)'
      scale_jl(k) = -1.
      ia_jl(k) = ia_rad
      pow_jl(k) = -2
c
      k=k+1
      jl_sshr = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'lscond_heat' !
      lname_jl(k) = 'HEATING BY LARGE SCALE CONDENSATION' !'DTX(SS)*DP'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_trbhr = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'turb_heat' !
      lname_jl(k) = 'HEATING BY TURBULENCE' !'DT(DC)*DP'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_dtdyn = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'dtempdt_dynamics' !
      lname_jl(k) = 'DTEMP/DT BY DYNAMICS'
      units_jl(k) = 'K/DAY'
      pow_jl(k) = -1
      scale_jl(k) = SECONDS_PER_DAY*NIDYN/(2*dtsrc) ! 1/dt_lf(days)
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_totcld = k
      sname_jl(k) = 'totcld' !
      lname_jl(k) = 'TOTAL CLOUD COVER' !'PCLD (TOTAL)'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_sscld = k
      sname_jl(k) = 'sscld' !
      lname_jl(k) = 'SUPER SATURATION CLOUD COVER' !'PCLD (SS)'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_mccld = k
      sname_jl(k) = 'mccld' !
      lname_jl(k) = 'MOIST CONVECTIVE CLOUD COVER' !'PCLD (MC)'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dtdtsdrg = k
      sname_jl(k) = 'dtempdt_sdrag' !
      lname_jl(k) = 'DTEMP/DT BY STRATOSPHERIC DRAG'
      units_jl(k) = 'K/DAY'
      pow_jl(k) = -1
      scale_jl(k) = SECONDS_PER_DAY/(IM*DTsrc)
      ia_jl(k) = ia_src
c
      k=k+1
      jl_damdc = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'dudt_dc' !
      lname_jl(k) = 'CHANGE OF U-WIND BY TURBULENCE'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
c
      k=k+1
      jl_dammc = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'dudt_mc' !          'DU(MC)*DP'
      lname_jl(k) = 'CHANGE OF U-WIND BY MOIST CONV'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
c
      k=k+1
      jl_epacwt = k
      sname_jl(k) = 'epac_wt' !
      lname_jl(k) = 'EAST PACIFIC WEIGHT'
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_uepac = k
      denom_jl(k) = jl_epacwt
      sname_jl(k) = 'u_epac' !
      lname_jl(k) = 'U WIND AVERAGED OVER EAST PACIFIC'
      units_jl(k) = '10^-1 m/s'
      scale_jl(k) = 10.
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_vepac = k
      denom_jl(k) = jl_epacwt
      sname_jl(k) = 'v_epac' !
      lname_jl(k) = 'V WIND AVERAGED OVER EAST PACIFIC'
      units_jl(k) = '10^-1 m/s'
      scale_jl(k) = 10.
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_wepac = k
      denom_jl(k) = jl_epacwt ! not quite right, but ok
      sname_jl(k) = 'vvel_epac' !
      lname_jl(k) = 'VERTICAL VELOCITY FOR EAST PACIFIC'
      units_jl(k) = '10**-5 m/s'
      scale_jl(k) = -1.D5*RGAS/(GRAV)
      ia_jl(k) = ia_dga
      lgrid_jl(k) = edg_ml
c
      k=k+1
      jl_wpacwt = k
      sname_jl(k) = 'wpac_wt' !
      lname_jl(k) = 'WEST PACIFIC WEIGHT'
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_uwpac = k
      denom_jl(k) = jl_wpacwt
      sname_jl(k) = 'u_wpac' !
      lname_jl(k) = 'U WIND AVERAGED OVER WEST PACIFIC'
      units_jl(k) = '10^-1 m/s'
      scale_jl(k) = 10.
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_vwpac = k
      denom_jl(k) = jl_wpacwt
      sname_jl(k) = 'v_wpac' !
      lname_jl(k) = 'V WIND AVERAGED OVER WEST PACIFIC'
      units_jl(k) = 'm/s'
      pow_jl(k) = -1
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_wwpac = k
      denom_jl(k) = jl_wpacwt ! not quite right, but ok
      sname_jl(k) = 'vvel_wpac' !
      lname_jl(k) = 'VERTICAL VELOCITY FOR WEST PACIFIC'
      units_jl(k) = '10**-5 m/s'
      scale_jl(k) = -1.D5*RGAS/(GRAV)
      ia_jl(k) = ia_dga
      lgrid_jl(k) = edg_ml
c
      k=k+1
      jl_dudtsdrg = k
      sname_jl(k) = 'dudt_sdrag' !
      lname_jl(k) = 'DU/DT BY SDRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./(DTsrc)
      jgrid_jl(k) = jgrid_u
c
      k=k+1
      jl_mcldht = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'moist_lat' !
      lname_jl(k) = 'CHANGE OF LATENT HEAT BY CTEI'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_trbke = k
      sname_jl(k) = 'tke' !
      lname_jl(k) = 'TURBULENT KINETIC ENERGY'
      units_jl(k) = 'm^2/s^2'
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_trbdlht = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'turb_lat' !
      lname_jl(k) = 'CHANGE OF LATENT HEAT BY TURBULENCE'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_mcheat = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'tot_ht_mc' !
      lname_jl(k) = 'TOTAL HEATING BY MOIST CONVECTION (Q1-QR)'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_mcdeep = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'tot_ht_deepmc' !
      lname_jl(k) = 'TOTAL HEATING BY DEEP MOIST CONVECTION (Q1-QR)'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_mcshlw = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'tot_ht_shlwmc' !
      lname_jl(k) = 'TOTAL HEATING BY SHALLOW MOIST CONVECTION (Q1-QR)'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_mcdry = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'tot_dry_mc' !
      lname_jl(k) = 'TOTAL DRYING BY MOIST CONVECTION (Q2)'
      units_jl(k) = 'W/(m^2*mb)'
      pow_jl(k) = -2
      scale_jl(k) = 100.*BYGRAV*SHA/DTsrc
c
      k=k+1
      jl_rhe= k
      sname_jl(k) = 'rhe'
      lname_jl(k) = 'EFFECTIVE RELATIVE HUMIDITY' ! from cloud scheme
      units_jl(k) = '%'
      scale_jl(k) = 100.
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_cldmc= k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'cldmc'   ! no output
      lname_jl(k) = 'MOIST CONVECTIVE CLOUD FRACTION' !from cloud scheme
      units_jl(k) = '%'
      scale_jl(k) = 100.
c
      k=k+1
      jl_cldss= k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'cldss'  ! no output
      lname_jl(k) = 'LARGE-SCALE CLOUD FRACTION' ! from cloud scheme
      units_jl(k) = '%'
      scale_jl(k) = 100.
c
      k=k+1
      jl_csizmc= k
      denom_jl(k) = jl_cldmc
      sname_jl(k) = 'csizmc'
      lname_jl(k) = 'MOIST CONVECTIVE EFFECTIVE CLOUD PARTICLE SIZE'
      units_jl(k) = 'micron'
c
      k=k+1
      jl_csizss= k
      denom_jl(k) = jl_cldss
      sname_jl(k) = 'csizss'
      lname_jl(k) = 'LARGE-SCALE EFFECTIVE CLOUD PARTICLE SIZE'
      units_jl(k) = 'micron'

#ifdef CLD_AER_CDNC
c ! Menon added diag for CDNC
c
      k=k+1
      jl_cnumwm= k
      denom_jl(k) = jl_cldmc
      sname_jl(k) = 'cnumwm'
      lname_jl(k) = 'WARM MOIST CONVECTIVE CLOUD DROPLET NUMBER'
      units_jl(k) = 'cm^-3'
c
      k=k+1
      jl_cnumws= k
      denom_jl(k) = jl_cldss
      sname_jl(k) = 'cnumws'
      lname_jl(k) = 'WARM LARGE-SCALE CLOUD DROPLET NUMBER'
      units_jl(k) = 'cm^-3'
c
      k=k+1
      jl_cnumim= k
      denom_jl(k) = jl_cldmc
      sname_jl(k) = 'cnumim'
      lname_jl(k) = 'COLD MOIST CONVECTIVE CLOUD DROPLET NUMBER'
      units_jl(k) = 'cm^-3'
c
      k=k+1
      jl_cnumis= k
      denom_jl(k) = jl_cldss
      sname_jl(k) = 'cnumis'
      lname_jl(k) = 'COLD LARGE-SCALE CLOUD DROPLET NUMBER'
      units_jl(k) = 'cm^-3'
#endif
c
      k=k+1
      jl_wcld = k
      sname_jl(k) = 'wcld' !
      lname_jl(k) = 'WATER CLOUD COVER' !'PCLD'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_icld = k
      sname_jl(k) = 'icld' !
      lname_jl(k) = 'ICE CLOUD COVER' !'PCLD'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_wcldwt = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'wcldwt' ! only used as weight.
      lname_jl(k) = 'WATER CLOUD COVER' !'PCLD*DP'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_icldwt = k
      denom_jl(k) = jl_dpa
      sname_jl(k) = 'icldwt' ! only used as weight.
      lname_jl(k) = 'ICE CLOUD COVER' !'PCLD*DP'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_wcod = k
      denom_jl(k) = jl_wcldwt
      sname_jl(k) = 'wcod' !
      lname_jl(k) = 'WATER CLOUD OPTICAL DEPTH' ! od*wcldcv*dp
      units_jl(k) = '/1000mb'
      scale_jl(k) = 1000.
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_icod = k
      denom_jl(k) = jl_icldwt
      sname_jl(k) = 'icod' !
      lname_jl(k) = 'ICE CLOUD OPTICAL DEPTH'   ! od*icldcv*dp
      units_jl(k) = '/1000mb'
      scale_jl(k) = 1000.
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_wcsiz= k
      denom_jl(k) = jl_wcod
      sname_jl(k) = 'wcsiz'
      lname_jl(k) = 'EFFECTIVE WATER CLOUD PARTICLE SIZE' ! SIZ*OPT.DPTH
      units_jl(k) = 'micron'
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_icsiz= k
      denom_jl(k) = jl_icod
      sname_jl(k) = 'icsiz'
      lname_jl(k) = 'EFFECTIVE ICE CLOUD PARTICLE SIZE' ! SIZ*OPT.DPTH
      units_jl(k) = 'micron'
      ia_jl(k) = ia_rad
c
      k=k+1
      jl_dwasrc = k
      sname_jl(k) = 'jl_dwasrc'
      lname_jl(k) = 'CONDENSATE MASS AT PRIMARY GRID CELLS'
      units_jl(k) = 'kg/m2'
      scale_jl(k) = 1. ! not printed

c
      k=k+1
      jk_dpwt = k
      sname_jl(k) = 'dp_cp'
      lname_jl(k) =  'PRESSURE DIFFERENCES (CP,PT)' ! DP (PT GRID)
      units_jl(k) = 'mb'
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_tx = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'tx' !'AJK03'
      lname_jl(k) = 'TEMPERATURE' !'(TX-TFREZ)*DP'
      units_jl(k) = 'C'
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      !jk_tx_rqt = k
      sname_jl(k) = 'tx_radonly' !'AJK03'
      lname_jl(k) = 'TEMPERATURE IN RADIATION-ONLY LAYERS'
      units_jl(k) = 'C'
c      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_hght = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'height' !'AJK04'
      lname_jl(k) = 'HEIGHT' !'PHI*DP'
      units_jl(k) = 'm'
      pow_jl(k) = 2
      scale_jl(k) = BYGRAV
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_q = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'q' !'AJK05'
      lname_jl(k) = 'SPECIFIC HUMIDITY' !'Q*DP'
      units_jl(k) = 'g/kg'
      pow_jl(k) = -2
      scale_jl(k) = 1.e3
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_rh = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'rh' !'AJK07'
      lname_jl(k) = 'RELATIVE HUMIDITY' !'RH*DP'
      units_jl(k) = '%'
      scale_jl(k) = 100.
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_cldh2o = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'cldh2o' !'AJK51'
      lname_jl(k) = 'TOTAL CLOUD WATER CONTENT'
      units_jl(k) = 'kg/kg'
      pow_jl(k) = -6
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_cldwtr = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'cldwtr'
      lname_jl(k) = 'CLOUD LIQUID WATER CONTENT'
      units_jl(k) = 'kg/kg'
      pow_jl(k) = -6
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jk_cldice = k
      denom_jl(k) = jk_dpwt
      sname_jl(k) = 'cldice'
      lname_jl(k) = 'CLOUD ICE WATER CONTENT'
      units_jl(k) = 'kg/kg'
      pow_jl(k) = -6
      lgrid_jl(k) = ctr_cp
      ia_jl(k) = ia_dga
c
      k=k+1
      jl_mchr = k
c      denom_jl(k) = jl_xxx
      sname_jl(k) = 'AJL13'
      lname_jl(k) = 'DT(MC)*P  DRY HEATING'
      units_jl(k) = '100 PA*K'
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_mchphas = k
c      denom_jl(k) = jl_xxx
      sname_jl(k) = 'AJL50'
      lname_jl(k) = 'DT(MC)*P  CHANGE OF PHASE'
      units_jl(k) = '100 PA*K'
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_mcdtotw = k
c      denom_jl(k) = jl_xxx
      sname_jl(k) = 'mc_del_tot_wat' !'CLHE*DQ(MC BEFORE COND)*P'
      lname_jl(k) = 'CHANGE IN TOTAL WATER BY MOIST CONV'
      units_jl(k) = '100 PA*K'
      force_jl_vmean(k) = .true. ! b/c no denom
c
      if (DO_GWDRAG) then

      k=k+1
      jl_gwFirst = k   ! The next consecutive 9 are Gravity Wave Diags
      jl_dumtndrg = k
      sname_jl(k) = 'dudt_mtndrg' !
      lname_jl(k) = 'DU/DT BY STRAT MTN DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      pow_jl_vmean(k) = -7
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dushrdrg = k
      sname_jl(k) = 'dudt_shrdrg'
      lname_jl(k) = 'DU/DT BY STRAT SHR DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dumcdrgm10 = k
      sname_jl(k) = 'dudt_mcdrgm10' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-10'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dumcdrgp10 = k
      sname_jl(k) = 'dudt_mcdrgp10' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+10'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dumcdrgm40 = k
      sname_jl(k) = 'dudt_mcdrgm40' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-40'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dumcdrgp40 = k
      sname_jl(k) = 'dudt_mcdrgp40' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+40'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dumcdrgm20 = k
      sname_jl(k) = 'dudt_mcdrgm20' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=-20'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dumcdrgp20 = k
      sname_jl(k) = 'dudt_mcdrgp20' !
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+20'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c Last of the Gravity Wave JL's
      k=k+1
      jl_dudfmdrg = k
      sname_jl(k) = 'dudt_dfmdrg' !
      lname_jl(k) = 'DU/DT BY STRAT DEFORM DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      ia_jl(k) = ia_src
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom

C**** Some extra GWDRAG related diags
c
      k=k+1
      jl_sdifcoef = k
      sname_jl(k) = 'strat_diff_coeff' !
      lname_jl(k) = 'STRAT. DIFFUSION COEFF'
      units_jl(k) = 'm^2/s'
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dudtsdif = k
      sname_jl(k) = 'dudt_sdiff' !     ! gwdrag
      lname_jl(k) = 'DU/DT BY GRAVITY WAVE MOMENTUM DIFFUSION'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k=k+1
      jl_dudtvdif = k
      sname_jl(k) = 'dudt_vdiff' !     ! vdiff
      lname_jl(k) = 'DU/DT BY VERTICAL DIFFUSION'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom

c combination GWDRAG diags
      k = k + 1
      jl_mcdrgpm10 = k
      sname_jl(k) = 'dudt_mcdrgpm10'
      lname_jl(k) = 'DU/DT BY STRAT. MC DRAG  C=+/-10R'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k = k + 1
      jl_mcdrgpm40 = k
      sname_jl(k) = 'dudt_mcdrgpm40' !AJL24+25
      lname_jl(k) = 'DU/DT BY STRAT. MC DRAG  C=+/-40R'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k = k + 1
      jl_mcdrgpm20 = k
      sname_jl(k) = 'dudt_mcdrgpm20' !AJL26+27
      lname_jl(k) = 'DU/DT BY STRAT MC DRAG C=+/-20R'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom
c
      k = k + 1
      jl_sumdrg = k
      sname_jl(k) = 'dudt_sumdrg' !AJL(18+20-27)
      lname_jl(k) = 'ZONAL WIND CHANGE BY MTN+DEFORM+SHR+MC DRAG'
      units_jl(k) = 'm/s^2'
      pow_jl(k) = -6
      scale_jl(k) = 1./DTsrc
      jgrid_jl(k) = jgrid_u
      force_jl_vmean(k) = .true. ! b/c no denom

      end if

      if (AM_I_ROOT()) then
         if (k .gt. kajl) then
            write (6,*) 'jl_defs: Increase kajl=',kajl,' to at least ',k
            call stop_model( 'kajl too small', 255 )
         end if

         write (6,*) 'Number of AJL diagnostics defined: kajlmax=',k
         if(qcheck) then
           do kk=1,k
             write (6,'(i4,'':'',a)') kk,trim(lname_jl(kk))
           end do
         endif
      end if

      do k=1,kajl
        if(denom_jl(k).ne.0) cycle ! already set
        if(lgrid_jl(k).eq.edg_ml) force_jl_vmean(k) = .true.
      enddo

c
c Declare the dimensions and metadata of AJL output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      call add_coord(cdl_heights,'plm',lm,units='mb',
     &     coordvalues=plm(1:lm))
      call add_coord(cdl_heights,'ple',lm,units='mb',
     &     coordvalues=ple)
      do l=1,lm
        one_to_lm(l) = l
      enddo
      call add_coord(cdl_heights,'level',lm,
     &     coordvalues=one_to_lm)
      call merge_cdl(cdl_latbudg,cdl_heights,cdl_jl_template)

      cdl_jl = cdl_jl_template ! invoke a copy method later

      call add_coord(cdl_jl,'p_radonly',lm_req,units='mb',
     &     coordvalues=pmidl00(lm+1:lm+lm_req))

      do k=1,kajl
        if(trim(units_jl(k)).eq.'unused') cycle
c        call get_zstr(lgrid_jl(k),zstr)
        if(lgrid_jl(k).eq.ctr_ml .or. lgrid_jl(k).eq.ctr_cp) then
          zstr='plm'
        else
          zstr='ple'
        endif
        set_miss = denom_jl(k).ne.0
        call add_var(cdl_jl, 'float '//trim(sname_jl(k))//'('//
     &       trim(zstr)//',lat_budg) ;',
     &       units=trim(units_jl(k)),
     &       long_name=trim(lname_jl(k)),
     &       auxvar_string=
     &         'float '//trim(sname_jl(k))//'_hemis('//
     &          trim(zstr)//',shnhgm) ;',
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
        if(pow_jl(k).ne.0) then
          write(powstr,'(i2)') pow_jl(k)
          call add_varline(cdl_jl,
     &         trim(sname_jl(k))//':prtpow = '//trim(powstr)//' ;')
        endif
        if(pow_jl_vmean(k).ne.0) then
          write(powstr,'(i2)') pow_jl_vmean(k)
          call add_varline(cdl_jl,
     &        trim(sname_jl(k))//':prtpow_vmean = '//trim(powstr)//' ;')
        endif
#ifndef SCM
        if(denom_jl(k).gt.0 .or. force_jl_vmean(k)) then
          if(make_timeaxis) then ! hacky logic
            call add_var(cdl_jl, 'float '//trim(sname_jl(k))//
     &           '_vmean(time,lat_budg_plus3) ;')
          else
            call add_var(cdl_jl, 'float '//trim(sname_jl(k))//
     &           '_vmean(lat_budg_plus3) ;')
          endif
        endif
#endif
      enddo

      return
      end subroutine jl_defs

      subroutine sjl_defs
      use CONSTANT, only : grav,sha,bygrav
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : qcheck
      use DIAG_COM
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      implicit none
      integer :: k,kk
c
      do k=1,kasjl
         write(name_sjl(k),'(a4,i3.3)') 'ASJL',k
         lname_sjl(k) = 'unused'
         units_sjl(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_sjl(k) = 'ASJL01'
      lname_sjl(k) = 'TX'
      units_sjl(k) = 'C'
      scale_sjl(k) = 1.
      ia_sjl(k) = ia_dga
c
      k=k+1
      name_sjl(k) = 'ASJL02'
      lname_sjl(k) = 'HEIGHT'
      units_sjl(k) = '100 M'
      scale_sjl(k) = .01*BYGRAV
      ia_sjl(k) = ia_dga
c
      k=k+1
      name_sjl(k) = 'srad_heat' !'ASJL03'
      lname_sjl(k) = 'SOLAR RADIATION HEATING RATE' !'SRHR'
      units_sjl(k) = '10**-2 K/DAY' !'W/m^2'
      scale_sjl(k) = 100.D-2*GRAV*SECONDS_PER_DAY/SHA
      ia_sjl(k) = ia_rad
c
      k=k+1
      name_sjl(k) = 'trad_cool' !'ASJL04'
      lname_sjl(k) = 'THERMAL RADIATION COOLING RATE' !'TRHR'
      units_sjl(k) = '10**-2 K/DAY' !'W/m^2'
      scale_sjl(k) = -100.D-2*GRAV*SECONDS_PER_DAY/SHA
      ia_sjl(k) = ia_rad
c
      k=k+1
      name_sjl(k) = 'rad_cool' ! not accumulated
      lname_sjl(k) = 'TOTAL RADIATION COOLING RATE'
      units_sjl(k) = 'W/(m^2*mb)'
      scale_sjl(k) = -100.
      ia_sjl(k) = ia_rad
c
      if (AM_I_ROOT()) then
         write (6,*) 'Number of ASJL diagnostics defined: kasjlmax=',k
         if(.not.qcheck) return
         do kk=1,k
            write (6,'(i4,'':'',a)') kk,trim(lname_sjl(kk))
         end do
      end if
         return
      end subroutine sjl_defs

      subroutine ijl_defs
      use CONSTANT, only : bygrav,sha,rgas
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : dtsrc
      use DIAG_COM
      use DIAG_COM_rad
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use cdl_mod
      use MDIAG_COM, only : make_timeaxis
      implicit none
      integer :: k,kk
      character(len=16) :: zstr,hstr,tstr
      logical :: set_miss
c
      do k=1,kaijl
         write(name_ijl(k),'(a4,i3.3)') 'AIJL',k
         lname_ijl(k) = 'no output'
         units_ijl(k) = 'unused'
         scale_ijl(k) = 1.
         denom_ijl(k) = 0
         ia_ijl(k) = ia_src
         lgrid_ijl(k) = ctr_ml
         jgrid_ijl(k) = 1
      enddo

c
      k=0
c
      k=k+1
      IJL_DP = k  ! not output - this field serves only as a weight
      ia_ijl(k) = ia_dga
c
      k=k+1
      IJK_DP = k  ! not output - this field serves only as a weight
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
#ifdef AIJL_CP_TRANSPORTS /* for vertical integrals.  output always? */
      name_ijl(k) = 'dpcp'
      lname_ijl(k) = 'pressure thickness'
      units_ijl(k) = 'mb'
#endif
c
      k=k+1
      IJL_U = k   ! e-w wind on model layers
      ia_ijl(k) = ia_dga
#ifdef CUBED_SPHERE
! latlon B-grid config only wants this in DIAGIL
      name_ijl(k) = 'u'
      lname_ijl(k) = 'east-west velocity'
      units_ijl(k) = 'm/s'
#endif
c
      k=k+1
      IJL_V = k   ! n-s wind on model layers
      ia_ijl(k) = ia_dga
#ifdef CUBED_SPHERE
! latlon B-grid config only wants this in DIAGIL
      name_ijl(k) = 'v'
      lname_ijl(k) = 'north-south velocity'
      units_ijl(k) = 'm/s'
#endif
c
      k=k+1
      IJL_W = k   ! vertical velocity
      ia_ijl(k) = ia_dga
#ifdef CUBED_SPHERE
! latlon B-grid config only wants this in DIAGIL
      name_ijl(k) = 'w'
      lname_ijl(k) = 'vertical velocity'
      units_ijl(k) = 'Pa/s'
      scale_ijl(k) = 100./DTsrc
      lgrid_ijl(k) = edg_cp ! constant pressure levels
#endif
c
      k=k+1
      IJK_TX = k
      name_ijl(k) = 'temp'
      lname_ijl(k) = 'TEMPERATURE'
      units_ijl(k) = 'C'
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1
      IJK_Q = k
      name_ijl(k) = 'q' !'QDP'
      lname_ijl(k) = 'SPECIFIC HUMIDITY'
      units_ijl(k) = 'kg/kg'
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1
      IJK_RH = k
      name_ijl(k) = 'rh'
      lname_ijl(k) = 'RELATIVE HUMIDITY'
      units_ijl(k) = '%'
      scale_ijl(k) = 100.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp
c
      k=k+1
      IJL_RC = k   ! no 3D output yet - presently used only in DIAGIL
      ia_ijl(k) = ia_rad
      denom_ijl(k) = IJL_DP
c
      k=k+1
      IJL_MC = k   ! no 3D output yet - presently used only in DIAGIL
      ia_ijl(k) = ia_src
      denom_ijl(k) = IJL_DP
c
      k=k+1
      IJL_CF=k
      name_ijl(k) = 'cf'
      lname_ijl(k) = 'CLOUD FRACTION'
      units_ijl(k) = '%'
      scale_ijl(k) = 100.
      ia_ijl(k) = ia_rad
c
      k=k+1        ! moist convective air mass flux (model layers)
      IJL_MCamFX   = k
      name_ijl(k)  = 'MCamFX'
      lname_ijl(k) = 'MOIST CONVECTIVE AIR MASS FLUX'
      units_ijl(k) = '1e-4 kg/m2/s'
      scale_ijl(k) = 1.e4*100.*BYGRAV/DTsrc
      denom_ijl(k) = 0
      ia_ijl(k)    = ia_src
      lgrid_ijl(k) = edg_ml
c
      k=k+1        ! mass fraction of cloud liquid water (model layers)
      IJL_cldwtr   = k
      name_ijl(k)  = 'wtrcld'
      lname_ijl(k) = 'Cloud Liquid Water Mass Mixing Ratio'
      units_ijl(k) = 'kg/kg'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJL_DP
      ia_ijl(k)    = ia_src
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! mass fraction of cloud ice (model layers)
      IJL_cldice   = k
      name_ijl(k)  = 'icecld'
      lname_ijl(k) = 'Cloud Ice Water Mass Mixing Ratio'
      units_ijl(k) = 'kg/kg'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJL_DP
      ia_ijl(k)    = ia_src
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! mass fraction of cloud+precip liquid seen by radiation (model layers)
      IJL_QLrad    = k
      name_ijl(k)  = 'QLrad'
      lname_ijl(k) = 'Liquid Water Mass Mixing Ratio Seen by Radiation'
      units_ijl(k) = 'kg/kg'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJL_DP
      ia_ijl(k)    = ia_rad
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! mass fraction of cloud+precip ice seen by radiation (model layers)
      IJL_QIrad    = k
      name_ijl(k)  = 'QIrad'
      lname_ijl(k) = 'Ice Water Mass Mixing Ratio Seen by Radiation'
      units_ijl(k) = 'kg/kg'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJL_DP
      ia_ijl(k)    = ia_rad
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! opacity of cloud water seen by radiation (model layers)
      IJL_wtrtau   = k
      name_ijl(k)  = 'wtrtau'
      lname_ijl(k) = 'Cloud Water Opacity Seen by Radiation'
      units_ijl(k) = '-'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_rad
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! opacity of cloud ice seen by radiation (model layers)
      IJL_icetau   = k
      name_ijl(k)  = 'icetau'
      lname_ijl(k) = 'Cloud Ice Opacity Seen by Radiation'
      units_ijl(k) = '-'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_rad
      lgrid_ijl(k) = ctr_ml
c
CC    3D drying and latent heating profiles
      if (lh_diags.eq.1) then

      k=k+1
      IJL_LLH=k
      name_ijl(k) = 'LLH'
      lname_ijl(k) = 'Heating by Large Scale Condensation'
      units_ijl(k) = 'C/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
      denom_ijl(k) = IJL_DP
c
      k=k+1
      IJL_MCTLH=k
      name_ijl(k) = 'CTLH'
      lname_ijl(k) = 'Heating by Moist Convection'
      units_ijl(k) = 'C/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
      denom_ijl(k) = IJL_DP
c
      k=k+1
      IJL_MCDLH=k
      name_ijl(k) = 'CDLH'
      lname_ijl(k) = 'Heating by Deep Convection'
      units_ijl(k) = 'C/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
      denom_ijl(k) = IJL_DP
c
      k=k+1
      IJL_MCSLH=k
      name_ijl(k) = 'CSLH'
      lname_ijl(k) = 'Heating by Shallow Convection'
      units_ijl(k) = 'C/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
      denom_ijl(k) = IJL_DP
c
      k=k+1
      IJL_LDRY=k
      name_ijl(k) = 'LSCDRY'
      lname_ijl(k) = 'Drying by Large Scale Condensation (Q2)'
      units_ijl(k) = 'kg/kg/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
c      
      k=k+1
      IJL_TMCDRY=k
      name_ijl(k) = 'TMCDRY'
      lname_ijl(k) = 'Drying by Moist Convection (Q2)'
      units_ijl(k) = 'kg/kg/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
c      
      k=k+1
      IJL_DMCDRY=k
      name_ijl(k) = 'DMCDRY'
      lname_ijl(k) = 'Drying by Deep Convection (Q2)'
      units_ijl(k) = 'kg/kg/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src
c
      k=k+1
      IJL_SMCDRY=k
      name_ijl(k) = 'SMCDRY'
      lname_ijl(k) = 'Drying by Shallow Convection (Q2)'
      units_ijl(k) = 'kg/kg/d'
      scale_ijl(k) = SECONDS_PER_DAY/DTsrc
      ia_ijl(k) = ia_src

      endif ! lh_diags==1

c
#ifdef CLD_AER_CDNC
      k=k+1
      IJL_CFWM=k
      name_ijl(k) = 'cfwm'
      lname_ijl(k) = 'Warm C Cloud Fraction'
      units_ijl(k) = '1'
      ia_ijl(k) = ia_src
c
      k=k+1
      IJL_CFWS=k
      name_ijl(k) = 'cfws'
      lname_ijl(k) = 'Warm S Cloud Fraction'
      units_ijl(k) = '1'
      ia_ijl(k) = ia_src
c
      k=k+1
      IJL_CFIM=k
      name_ijl(k) = 'cfim'
      lname_ijl(k) = 'Cold C Cloud Fraction'
      units_ijl(k) = '1'
      ia_ijl(k) = ia_src
c
      k=k+1
      IJL_CFIS=k
      name_ijl(k) = 'cfis'
      lname_ijl(k) = 'Cold S Cloud Fraction'
      units_ijl(k) = '1'
      ia_ijl(k) = ia_src
c
      k=k+1
      IJL_REWM=k
      name_ijl(k) = 'rewm'
      lname_ijl(k) = 'Warm C Reff'
      units_ijl(k) = 'um'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfwm
c
      k=k+1
      IJL_REWS=k
      name_ijl(k) = 'rews'
      lname_ijl(k) = 'Warm S Reff'
      units_ijl(k) = 'um'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfws
c
      k=k+1
      IJL_CDWM=k
      name_ijl(k) = 'cdwm'
      lname_ijl(k) = 'Warm C CDNC'
      units_ijl(k) = 'cm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfwm
c
      k=k+1
      IJL_CDWS=k
      name_ijl(k) = 'cdws'
      lname_ijl(k) = 'Warm S CDNC'
      units_ijl(k) = 'cm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfws

#ifdef TRACERS_TOMAS

      k=k+1
      IJL_CDTOMAS=k
      name_ijl(k) = 'CDNC_TOMAS'
      lname_ijl(k) = 'Warm S CDNC TOMAS'
      units_ijl(k) = 'cm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfws

#endif
c
      k=k+1
      IJL_CWWM=k
      name_ijl(k) = 'cwwm'
      lname_ijl(k) = 'Warm C LWC'
      units_ijl(k) = 'gm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfwm
c
      k=k+1
      IJL_CWWS=k
      name_ijl(k) = 'cwws'
      lname_ijl(k) = 'Warm S LWC'
      units_ijl(k) = 'gm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfws
c
      k=k+1
      IJL_REIM=k
      name_ijl(k) = 'reim'
      lname_ijl(k) = 'Cold C Reff'
      units_ijl(k) = 'um'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfim
c
      k=k+1
      IJL_REIS=k
      name_ijl(k) = 'reis'
      lname_ijl(k) = 'Cold S Reff'
      units_ijl(k) = 'um'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfis
c
      k=k+1
      IJL_CDIM=k
      name_ijl(k) = 'cdim'
      lname_ijl(k) = 'Cold C CDNC'
      units_ijl(k) = 'cm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfim
c
      k=k+1
      IJL_CDIS=k
      name_ijl(k) = 'cdis'
      lname_ijl(k) = 'Cold S CDNC'
      units_ijl(k) = 'cm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfis
c
      k=k+1
      IJL_CWIM=k
      name_ijl(k) = 'cwim'
      lname_ijl(k) = 'Cold C LWC'
      units_ijl(k) = 'gm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfim
c
      k=k+1
      IJL_CWIS=k
      name_ijl(k) = 'cwis'
      lname_ijl(k) = 'Cold S LWC'
      units_ijl(k) = 'gm-3'
      ia_ijl(k) = ia_src
      denom_ijl(k) = ijl_cfis
#endif /* CLD_AER_CDNC */
c
      k=k+1        ! temperature (model layers)
      ijl_tempL    = k
      name_ijl(k)  = 'TempL'
      lname_ijl(k) = 'Layer Temperature'
      units_ijl(k) = 'K'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_dga
      lgrid_ijl(k) = ctr_ml
#ifdef ACCMIP_LIKE_DIAGS
c
      k=k+1        ! grid box geometric thickness (model layers)
      ijl_gridh    = k
      name_ijl(k)  = 'GridH'
      lname_ijl(k) = 'Grid Box Geom Thickness'
      units_ijl(k) = 'm'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_dga
      lgrid_ijl(k) = ctr_ml
#endif
c
      k=k+1        ! specific humidity (model layers)
      ijl_husl     = k
      name_ijl(k)  = 'SpHuL'
      lname_ijl(k) = 'Specific Humidity'
      units_ijl(k) = 'kg/kg'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_dga
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! grid box geometric height (model layers)
      ijl_zL       = k
      name_ijl(k)  = 'z'
      lname_ijl(k) = 'Grid Box Geom Height'
      units_ijl(k) = 'm'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_dga
      lgrid_ijl(k) = ctr_ml
c
      k=k+1        ! Air mass on model layers
      ijl_airmass  = k
      name_ijl(k)  = 'airmass'
      lname_ijl(k) = 'Air Mass'
      units_ijl(k) = 'kg/m2/layer'
      scale_ijl(k) = 1.
      ia_ijl(k)    = ia_src ! to match taijl one
      lgrid_ijl(k) = ctr_ml
#ifdef AIJL_CP_TRANSPORTS
c
      k=k+1        ! u on constant-pressure layers
      ijk_ucp      = k
      name_ijl(k)  = 'ucp'
      lname_ijl(k) = 'Eastward velocity'
      units_ijl(k) = 'm/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! v on constant-pressure layers
      ijk_vcp      = k
      name_ijl(k)  = 'vcp'
      lname_ijl(k) = 'Northward velocity'
      units_ijl(k) = 'm/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! u*t on constant-pressure layers
      ijk_utcp     = k
      name_ijl(k)  = 'utcp'
      lname_ijl(k) = 'Eastward temperature flux'
      units_ijl(k) = 'K * m/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! v*t on constant-pressure layers
      ijk_vtcp     = k
      name_ijl(k)  = 'vtcp'
      lname_ijl(k) = 'Northward temperature flux'
      units_ijl(k) = 'K * m/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! u*q on constant-pressure layers
      ijk_uqcp     = k
      name_ijl(k)  = 'uqcp'
      lname_ijl(k) = 'Eastward humidity flux'
      units_ijl(k) = 'kg/kg * m/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! v*q on constant-pressure layers
      ijk_vqcp     = k
      name_ijl(k)  = 'vqcp'
      lname_ijl(k) = 'Northward humidity flux'
      units_ijl(k) = 'kg/kg * m/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! u*phi on constant-pressure layers
      ijk_uphicp   = k
      name_ijl(k)  = 'uphicp'
      lname_ijl(k) = 'Eastward geopotential flux'
      units_ijl(k) = 'm2/s2 * m/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
      k=k+1        ! v*phi on constant-pressure layers
      ijk_vphicp   = k
      name_ijl(k)  = 'vphicp'
      lname_ijl(k) = 'Northward geopotential flux'
      units_ijl(k) = 'm2/s2 * m/s'
      scale_ijl(k) = 1.
      denom_ijl(k) = IJK_DP
      ia_ijl(k) = ia_dga
      lgrid_ijl(k) = ctr_cp ! constant pressure levels
c
#endif
c
      if (k .gt. kaijl) then
        if(am_i_root())
     &       write (6,*) 'ijl_defs: Increase kaijl=',kaijl,' to ',k
        call stop_model( 'kaijl too small', 255 )
      end if
      if(AM_I_ROOT())
     &     write (6,*) 'Number of AIJL diagnostics defined: kaijlmax=',k

c
c Declare the dimensions and metadata of AIJL output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).
c
      call merge_cdl(cdl_ij_template,cdl_heights,cdl_ijl_template)
      cdl_ijl = cdl_ijl_template ! invoke a copy method later

#ifdef CUBED_SPHERE
      call merge_cdl(cdl_ij_latlon_template,cdl_heights,
     &     cdl_ijl_latlon_template)
      cdl_ijl_latlon = cdl_ijl_latlon_template ! invoke a copy method later
      tstr='(tile,'
      hstr=',y,x) ;'
#else
      tstr='('
      hstr=',lat,lon) ;'
#endif
      do k=1,kaijl
        if(trim(units_ijl(k)).eq.'unused') cycle
        call get_zstr(lgrid_ijl(k),zstr)
        set_miss = denom_ijl(k).ne.0
        call add_var(cdl_ijl,
     &       'float '//trim(name_ijl(k))//
     &       trim(tstr)//trim(zstr)//trim(hstr),
     &       units=trim(units_ijl(k)),
     &       long_name=trim(lname_ijl(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#ifdef CUBED_SPHERE
        call add_var(cdl_ijl_latlon, 'float '//
     &       trim(name_ijl(k))//'('//trim(zstr)//',lat,lon);',
     &       units=trim(units_ijl(k)),
     &       long_name=trim(lname_ijl(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
#endif
      enddo

      return
      end subroutine ijl_defs

#ifndef SCM
      subroutine wave_defs
      use GC_COM
      use MODEL_COM, only : qcheck
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      implicit none
      integer :: k,kk
c
      do k=1,kwp
         write(name_wave(k),'(a5,i3.3)') 'WAVEP',k
         lname_wave(k) = 'unused'
         units_wave(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_wave(k) = 'U850_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V850_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'U300_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V300_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'U050_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'V050_EQ'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z922_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z700_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z500_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z300_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z100_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      k=k+1
      name_wave(k) = 'Z010_50N'
      lname_wave(k) = 'unknown'
      units_wave(k) = 'unknown'
c
      if (AM_I_ROOT()) then
         write (6,*) 'Number of Wave diagnostics defined: kwavemax=',k
         if(.not.qcheck) return
         do kk=1,k
            write (6,'(i4,'':'',a)') kk,trim(name_wave(kk))
         end do
      end if
      return
      end subroutine wave_defs
#endif

      subroutine tsf_defs
      use DIAG_COM
      implicit none
      integer :: k
c
      do k=1,ktsf
         write(name_tsf(k),'(a6,i3.3)') 'TSFREZ',k
         lname_tsf(k) = 'unused'
         units_tsf(k) = 'unused'
      enddo
c
      k=0
c
      k=k+1
      name_tsf(k) = 'TSFREZ1'
      lname_tsf(k) = 'FIRST DAY OF GROWING SEASON'
      units_tsf(k) = 'JULIAN DAY'
      tf_day1 = k
c
      k=k+1
      name_tsf(k) = 'TSFREZ2'
      lname_tsf(k) = 'LAST DAY OF GROWING SEASON'
      units_tsf(k) = 'JULIAN DAY'
      tf_last = k
c
      k=k+1
      name_tsf(k) = 'LKICEON'
      lname_tsf(k) = 'LAST DAY OF ICE-FREE LAKE'
      units_tsf(k) = 'JULIAN DAY'
      tf_lkon = k
c
      k=k+1
      name_tsf(k) = 'LKICEOFF'
      lname_tsf(k) = 'LAST DAY OF ICED-UP LAKE'
      units_tsf(k) = 'JULIAN DAY'
      tf_lkoff = k
c
      return
      end subroutine tsf_defs

      subroutine diurn_defs
!@sum  diurn_defs definitions for diurnal diagnostic accumulated arrays
!@auth G. Schmidt
      use CONSTANT, only : sha,rgas,twopi,grav
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use MODEL_COM, only : dtsrc,qcheck
      use FLUXES, only : nisurf
      use DIAG_COM
      use DIAG_COM_RAD
      use SOCPBL, only : npbl=>n
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use cdl_mod
#ifdef TRACERS_AMP
      use AERO_CONFIG, ONLY: nmodes
#endif
      implicit none
      integer :: k,kk,l, lmax_dd0=5 ! why?
      character*2 lst(lm)
      real*8 :: dummy_hrs(hr_in_month)
      logical :: set_miss

C**** define levels strings
      do l=1,lm
        if (l.lt.10) write(lst(l)(1:2),'(I1,1X)') l
        if (l.ge.10) write(lst(l)(1:2),'(I2)') l
      end do
c
      do k=1,NDIUVAR
         write(name_dd(k),'(a5,i3.3)') 'DIURN',k
         lname_dd(k) = 'unused'
         units_dd(k) = 'unused'
         scale_dd(k) = 1.
         denom_dd(k) = 0
      enddo
c
      k=0
c
      k=k+1
      IDD_ISW=k
      name_dd(k)='INCSW'
      units_dd(k)='W/m^2'
      scale_dd(k)=1.
      lname_dd(k)=' INC SW RADIATION'
c
      k=k+1
      IDD_PALB=k
      name_dd(k)='PLALB'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' P ALBD '
c
      k=k+1
      IDD_GALB=k
      name_dd(k)='GRDALB'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' G ALBD '
c
      k=k+1
      IDD_ABSA=k   ! absorbed solar radiation only
      name_dd(k)='ABSATM'
      units_dd(k)='W/m^2'
      scale_dd(k)=1.
      lname_dd(k)=' ABS ATM'
c
      k=k+1
      IDD_ECND=k
      name_dd(k)='ENRGCND'
      units_dd(k)='W/m^2'
      scale_dd(k)=100.*SHA/(GRAV*DTsrc)
      lname_dd(k)=' E CNDS '
c
      k=k+1
      IDD_SPR=k
      name_dd(k)='SRFPRS'
      units_dd(k)='mb'
      scale_dd(k)=1.
      lname_dd(k)=' SRF PRS'
c
        IDD_PT5=k+1
        do l=lmax_dd0,1,-1
          k=k+1
          name_dd(k)='POTT'//lst(l)
          units_dd(k)='K'
          scale_dd(k)=1.
          lname_dd(k)=' PT '//lst(l)
        end do
c
      k=k+1
      IDD_TS=k
      name_dd(k)='TSURF'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' TS     '
c
      k=k+1
      IDD_TG1=k
      name_dd(k)='TGRND'
      units_dd(k)='K'
      scale_dd(k)=1.
      lname_dd(k)=' TG1    '
c
        IDD_Q5=k+1
        do l=lmax_dd0,1,-1
          k=k+1
          name_dd(k)='Q'//lst(l)
          units_dd(k)='1d-5 kg/kg'
          scale_dd(k)=1d5
          lname_dd(k)=' Q '//lst(l)
        end do
c
      k=k+1
      IDD_QS=k
      name_dd(k)='QSURF'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' QS     '
c
      k=k+1
      IDD_QG=k
      name_dd(k)='QGRND'
      units_dd(k)='1d-5 kg/kg'
      scale_dd(k)=1d5
      lname_dd(k)=' QG     '
c
      k=k+1
      IDD_SWG=k
      name_dd(k)='SWGRND'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' SW ON G'
c
      k=k+1
      IDD_LWG=k
      name_dd(k)='LWGRND'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' LW AT G'
c
      k=k+1
      IDD_SH=k
      name_dd(k)='SENSHT'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' SNSB HT'
c
      k=k+1
      IDD_LH=k
      name_dd(k)='LATHT'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' LAT HT '
c
      k=k+1
      IDD_HZ0=k
      name_dd(k)='NETHT'
      units_dd(k)='W/m^2'
      scale_dd(k)=NISURF/DTsrc
      lname_dd(k)=' HEAT Z0'
c
      k=k+1
      IDD_UG=k
      name_dd(k)='UGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' UG*10  '
c
      k=k+1
      IDD_VG=k
      name_dd(k)='VGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' VG*10  '
c
      k=k+1
      IDD_WG=k
      name_dd(k)='WGEOS'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WG*10  '
c
      k=k+1
      IDD_US=k
      name_dd(k)='USURF'
      units_dd(k)='0.1 m/s'
      scale_dd(k)=10.
      lname_dd(k)=' US*10  '
c
      k=k+1
      IDD_VS=k
      name_dd(k)='VSURF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' VS*10  '
c
      k=k+1
      IDD_WS=k
      name_dd(k)='WSURF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WS*10  '
c
      k=k+1
      IDD_CIA=k
      name_dd(k)='CRISOANG'
      units_dd(k)=''
      scale_dd(k)=360./TWOPI
      lname_dd(k)=' ALPHA0 '
c
      k=k+1
      IDD_RIS=k
      name_dd(k)='RIS1'
      units_dd(k)='0.01'
      scale_dd(k)=100.
      lname_dd(k)=' RIS1*E2'
c
      k=k+1
      IDD_RIG=k
      name_dd(k)='RIGS'
      units_dd(k)='0.01'
      scale_dd(k)=100.
      lname_dd(k)=' RIGS*E2'
c
      k=k+1
      IDD_CM=k
      name_dd(k)='CM'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CM*E4  '
c
      k=k+1
      IDD_CH=k
      name_dd(k)='CH'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CH*E4  '
c
      k=k+1
      IDD_CQ=k
      name_dd(k)='CQ'
      units_dd(k)='1d-4'
      scale_dd(k)=1d4
      lname_dd(k)=' CQ*E4  '
c
      k=k+1
      IDD_EDS=k
      name_dd(k)='EDS'
      units_dd(k)='0.1'
      scale_dd(k)=10.
      lname_dd(k)=' EDS1*10'
c
      k=k+1
      IDD_DBL=k
      name_dd(k)='DBL'
      units_dd(k)='m'
      scale_dd(k)=1.
      lname_dd(k)=' DBL    '
c
      k=k+1
      IDD_DCF=k
      name_dd(k)='DCFREQ'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' DC FREQ'
c
      k=k+1
      IDD_LDC=k
      denom_dd(k) = IDD_DCF
      name_dd(k)='LDC'
      units_dd(k)='0.1 L'
      scale_dd(k)=10.
      lname_dd(k)=' LDC*10 '
c
      k=k+1
      IDD_PR=k
      name_dd(k)='PREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SECONDS_PER_DAY/(DTsrc*GRAV)
      lname_dd(k)=' PRC*100'  ! check scale
c
      k=k+1
      IDD_EV=k
      name_dd(k)='EVAP'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*SECONDS_PER_DAY*NISURF/DTsrc
      lname_dd(k)=' EVP*100'
c
      k=k+1
      IDD_DMC=k
      name_dd(k)='DEEPMC'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' DEEP MC'
c
      k=k+1
      IDD_SMC=k
      name_dd(k)='SHLWMC'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' SHLW MC'
c
        IDD_CL7=k+1
        do l=7,1,-1
          k=k+1
          name_dd(k)='CLD'//lst(l)
          units_dd(k)='%'
          scale_dd(k)=100.
          lname_dd(k)=' CLD '//lst(l)
        end do
c
!!    k=k+1
!!    IDD_W=k
!!    name_dd(k)='VERTVEL'
!!    units_dd(k)='1d-5 m/s'
!!    scale_dd(k)=1.
!!    lname_dd(k)=' W TO-5 '
c
      k=k+1
      IDD_CCV=k
      name_dd(k)='CLDCOV'
      units_dd(k)='%'
      scale_dd(k)=100.
      lname_dd(k)=' C COVER'
c
      k=k+1
      IDD_SSP=k
      name_dd(k)='SSPREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SECONDS_PER_DAY/(DTsrc*GRAV)
      lname_dd(k)=' SSP*100'
c
      k=k+1
      IDD_MCP=k
      name_dd(k)='MCPREC'
      units_dd(k)='0.01 mm/day'
      scale_dd(k)=100.*100.*SECONDS_PER_DAY/(DTsrc*GRAV)
      lname_dd(k)=' MCP*100'

c
      k=k+1
      IDD_aot=k
      name_dd(k)='AOT'
      units_dd(k)=' '
      scale_dd(k)= 1.
      lname_dd(k)=' AOT'
c
      k=k+1
      IDD_aot2=k
      name_dd(k)='QSCAT'
      units_dd(k)=' '
      scale_dd(k)= 1.
      lname_dd(k)='qscat'

#ifdef TRACERS_AMP
c
      k=k+1
      IDD_lwp=k
      name_dd(k)='LWP'
      units_dd(k)='kg/m2'
      scale_dd(k)= 1.
      lname_dd(k)=' LWP'

      k=k+1
      IDD_so2=k
      name_dd(k)='SO2'
      units_dd(k)='ug/m3'
      scale_dd(k)= 1.
      lname_dd(k)=' SO2'

c Mode Diagnostics
       idd_diam=k+1
        do l=1,nmodes
          k=k+1
          name_dd(k)='DIAM_M'//lst(l)
          units_dd(k)='m'
          scale_dd(k)=1.
          lname_dd(k)=' DIAM_M'//lst(l)
        end do

       idd_numb=k+1
        do l=1,nmodes
          k=k+1
          name_dd(k)='N_M'//lst(l)
          units_dd(k)='#/m3'
          scale_dd(k)=1.
          lname_dd(k)=' N_M'//lst(l)
        end do
c Mass Tracer
       idd_mass=k+1
        do l=1,38
          k=k+1
          name_dd(k)='MASS_M'//lst(l)
          units_dd(k)='ug/m3'
          scale_dd(k)=1.
          lname_dd(k)=' MASS_M'//lst(l)
        end do
c Column Diagnostics
       idd_ncL=k+1
        do l=1,lmax_dd0
          k=k+1
          name_dd(k)='TotNumb_L'//lst(l)
          units_dd(k)='#/cm3'
          scale_dd(k)=1.
          lname_dd(k)=' TotNumb_L'//lst(l)
        end do

       idd_ccn=k+1
        do l=1,lmax_dd0
          k=k+1
          name_dd(k)='CCN_L'//lst(l)
          units_dd(k)='#/cm3'
          scale_dd(k)=1.
          lname_dd(k)=' CCN_L'//lst(l)
        end do

       idd_cdnc=k+1
        do l=1,lmax_dd0
          k=k+1
          name_dd(k)='CDNC_L'//lst(l)
          units_dd(k)='#/cm3'
          scale_dd(k)=1.
          lname_dd(k)=' CDNC_L'//lst(l)
        end do

       idd_lwc=k+1
        do l=1,lmax_dd0
          k=k+1
          name_dd(k)='LWC_L'//lst(l)
          units_dd(k)='kg/m3'
          scale_dd(k)=1.
          lname_dd(k)=' LWC_L'//lst(l)
        end do

       idd_pres=k+1
        do l=1,lmax_dd0
          k=k+1
          name_dd(k)='Pres_L'//lst(l)
          units_dd(k)='hPa'
          scale_dd(k)=1.
          lname_dd(k)=' PRES_L'//lst(l)
        end do     
#endif

      k=k+1
      IDD_WSPDF=k
      name_dd(k)='WSPDF'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WSPDF*10'
c

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS)
      IF (adiurn_dust == 1) THEN
c
      k=k+1
      IDD_WTKE=k
      name_dd(k)='WTKE'
      units_dd(k)='m/s'
      scale_dd(k)=1.
      lname_dd(k)=' WTKE'
c
      k=k+1
      IDD_WD=k
      name_dd(k)='WD'
      units_dd(k)='m/s'
      scale_dd(k)=1.
      lname_dd(k)=' WD'
c
      k=k+1
      IDD_WM=k
      name_dd(k)='WM'
      units_dd(k)='m/s'
      scale_dd(k)=1.
      lname_dd(k)=' WM'
c
      k=k+1
      IDD_WSGCM=k
      name_dd(k)='WSGCM'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WSGCM*10'
c
      k=k+1
      IDD_WTRSH=k
      name_dd(k)='WTRSH'
      units_dd(k)='m/s'
      scale_dd(k)=10.
      lname_dd(k)=' WTRSH*10'
c
      END IF
#endif
#ifdef TRACERS_DUST
      IF (adiurn_dust == 1) THEN
c
        IDD_U1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='U_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' U_L'//lst(l)
        end do
c
        IDD_V1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='V_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' V_L'//lst(l)
        end do
c
        IDD_UV1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='UV_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' UV_L'//lst(l)
        end do
c
        IDD_T1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='T_L'//lst(l)
          units_dd(k)='K'
          scale_dd(k)=1.
          lname_dd(k)=' T_L'//lst(l)
        end do
c
        IDD_QQ1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='Q_L'//lst(l)
          units_dd(k)='kg/kg'
          scale_dd(k)=1.
          lname_dd(k)=' Q_L'//lst(l)
        end do
c
        IDD_P1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='P_L'//lst(l)
          units_dd(k)='100*mb'
          scale_dd(k)=100.
          lname_dd(k)=' P_L'//lst(l)
        end do
c
        IDD_W1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='W_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' W_L'//lst(l)
        end do
c
        IDD_PHI1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='PHI_L'//lst(l)
          units_dd(k)='m'
          scale_dd(k)=1.
          lname_dd(k)=' PHI_L'//lst(l)
        end do
c
        IDD_LOAD1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='LOAD_L'//lst(l)
          units_dd(k)='10^-5 kg/m**2'
          scale_dd(k)=1.e5
          lname_dd(k)=' LOAD_L'//lst(l)
        end do
c
        IDD_CONC1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='CONC_L'//lst(l)
          units_dd(k)='10^-8 kg/kg air'
          scale_dd(k)=1.e8
          lname_dd(k)=' LOAD_L'//lst(l)
        end do
c
      k=k+1
      IDD_EMIS=k
      name_dd(k)='EMIS'
      units_dd(k)='10^-13 kg/m^2/s'
      scale_dd(k)=1.e13
      lname_dd(k)=' EMIS'
c
      k=k+1
      IDD_EMIS2=k
      name_dd(k)='EMIS2'
      units_dd(k)='10^-13 kg/m^2/s'
      scale_dd(k)=1.e13
      lname_dd(k)=' EMIS2'
c
      k=k+1
      IDD_WET=k
      name_dd(k)='DEPWET'
      units_dd(k)='10^-13 kg/m^2/s'
      scale_dd(k)=1.e13
      lname_dd(k)=' DEPOWET'
c
      k=k+1
      IDD_GRAV=k
      name_dd(k)='DEPGRAV'
      units_dd(k)='10^-13 kg/m^2/s'
      scale_dd(k)=1.e13
      lname_dd(k)=' DEPOGRAV'
c
      k=k+1
      IDD_TURB=k
      name_dd(k)='DEPTURB'
      units_dd(k)='10^-13 kg/m^2/s'
      scale_dd(k)=1.e13
      lname_dd(k)=' DEPOTURB'
c
        IDD_TAU1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='TAU_L'//lst(l)
          units_dd(k)='1.'
          scale_dd(k)=1.
          lname_dd(k)=' TAU_L'//lst(l)
        end do
c
        IDD_TAU_CS1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='TAU_CS_L'//lst(l)
          units_dd(k)='1.'
          scale_dd(k)=1.
          lname_dd(k)=' TAU_CS_L'//lst(l)
        end do
c
        IDD_SR1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='SRNFLB_L'//lst(l)
          units_dd(k)='W/m**2'
          scale_dd(k)=1.
          lname_dd(k)=' SRNFLB_L'//lst(l)
        end do
c
        IDD_TR1=k+1
        do l=1,lmax_dd2
          k=k+1
          name_dd(k)='TRNFLB_L'//lst(l)
          units_dd(k)='W/m**2'
          scale_dd(k)=1.
          lname_dd(k)=' TRNFLB_L'//lst(l)
        end do
c
      k=k+1
      IDD_WS2=k
      name_dd(k)='WS2'
      units_dd(k)='m^2/s^2'
      scale_dd(k)=1.
      lname_dd(k)=' WS2'
c
      k=k+1
      IDD_USTAR=k
      name_dd(k)='USTAR'
      units_dd(k)='m/s'
      scale_dd(k)=1.
      lname_dd(k)=' USTAR'
c
      k=k+1
      IDD_US3=k
      name_dd(k)='USTAR3'
      units_dd(k)='m^3/s^3'
      scale_dd(k)=1.
      lname_dd(k)=' USTAR3'
c
      k=k+1
      IDD_STRESS=k
      name_dd(k)='WSTRESS'
      units_dd(k)='Nm^-2'
      scale_dd(k)=1.
      lname_dd(k)=' WSTRESS'
c
      k=k+1
      IDD_LMON=k
      name_dd(k)='LMONIN'
      units_dd(k)='m'
      scale_dd(k)=1.
      lname_dd(k)=' LMONIN'
c
      k=k+1
      IDD_RIFL=k
      name_dd(k)='RI_FLUX'
      units_dd(k)='1.'
      scale_dd(k)=1.
      lname_dd(k)=' RI_FLUX'
c
        IDD_ZPBL1=k+1
        do l=1,npbl
          k=k+1
          name_dd(k)='ZPBL_L'//lst(l)
          units_dd(k)='m'
          scale_dd(k)=1.
          lname_dd(k)=' ZPBL_L'//lst(l)
        end do
c
        IDD_UABL1=k+1
        do l=1,npbl
          k=k+1
          name_dd(k)='UABL_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' UABL_L'//lst(l)
        end do
c
        IDD_VABL1=k+1
        do l=1,npbl
          k=k+1
          name_dd(k)='VABL_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' VABL_L'//lst(l)
        end do
c
        IDD_UVABL1=k+1
        do l=1,npbl
          k=k+1
          name_dd(k)='UVABL_L'//lst(l)
          units_dd(k)='m/s'
          scale_dd(k)=1.
          lname_dd(k)=' UVABL_L'//lst(l)
        end do
c
        IDD_TABL1=k+1
        do l=1,npbl
          k=k+1
          name_dd(k)='TABL_L'//lst(l)
          units_dd(k)='K'
          scale_dd(k)=1.
          lname_dd(k)=' TABL_L'//lst(l)
        end do
c
        IDD_QABL1=k+1
        do l=1,npbl
          k=k+1
          name_dd(k)='QABL_L'//lst(l)
          units_dd(k)='kg/kg'
          scale_dd(k)=1.
          lname_dd(k)=' QABL_L'//lst(l)
        end do
c
        IDD_ZHAT1=k+1
        do l=1,npbl-1
          k=k+1
          name_dd(k)='ZHAT_L'//lst(l)
          units_dd(k)='m'
          scale_dd(k)=1.
          lname_dd(k)=' ZHAT_L'//lst(l)
        end do
c
        IDD_E1=k+1
        do l=1,npbl-1
          k=k+1
          name_dd(k)='TKE_L'//lst(l)
          units_dd(k)='m^2/s^2'
          scale_dd(k)=1.
          lname_dd(k)=' TKE_L'//lst(l)
        end do
c
        IDD_KM1=k+1
        do l=1,npbl-1
          k=k+1
          name_dd(k)='KM_L'//lst(l)
          units_dd(k)='m^2/s'
          scale_dd(k)=1.
          lname_dd(k)=' KM_L'//lst(l)
        end do
c
        IDD_RI1=k+1
        do l=1,npbl-1
          k=k+1
          name_dd(k)='RI_L'//lst(l)
          units_dd(k)='1.'
          scale_dd(k)=1.
          lname_dd(k)=' RI_L'//lst(l)
        end do
c
      END IF
#endif

      if (AM_I_ROOT()) then
         if (k .gt. Ndiuvar) then
            write (6,*) 'idd_defs: Increase Ndiuvar=',Ndiuvar,
     &           ' to at least ',k
            call stop_model( 'Ndiuvar too small', 255 )
         end if

         write (6,*) 'Number of Diurn diagnostics defined: kaddmax=',k
         if(qcheck) then
           do kk=1,k
             write (6,'(i4,'':'',a)') kk,trim(lname_dd(kk))
           end do
         endif
      end if

c
c Declare the dimensions and metadata of ADIURN output fields using
c netcdf CDL notation.  The C convention for dimension ordering
c must be used (reversed wrt Fortran).  Information needed for
c printing ASCII tables of the output is stored here as well.
c
      call init_cdl_type('cdl_dd',cdl_dd)
      call add_dim(cdl_dd,'ndiupt',ndiupt)
      call add_dim(cdl_dd,'namdd_strlen',4)
      call add_dim(cdl_dd,'two',2)
      call add_var(cdl_dd,'char namdd(ndiupt,namdd_strlen) ;')
      call add_vardata(cdl_dd,'namdd',namdd)
      call add_var(cdl_dd,'int ijdd(ndiupt,two) ;')
      call add_vardata(cdl_dd,'ijdd',
     &     reshape(ijdd,(/size(ijdd)/)) )
      do k=1,ndiuvar
        if(trim(lname_dd(k)).eq.'unused') cycle
        set_miss = denom_dd(k).ne.0
        call add_var(cdl_dd,
     &       'float '//trim(name_dd(k))//'(hour,ndiupt) ;',
     &       long_name=trim(lname_dd(k)),
     &       units=trim(units_dd(k)),
     &       set_miss=set_miss)
      enddo

      do k=1,hr_in_month
        dummy_hrs(k) = k
      enddo

#ifdef USE_HDIURN
c Declare the dimensions and metadata of HDIURN output fields
      cdl_hd = cdl_dd
      call add_coord(cdl_hd,'hour',hr_in_month,
     &     coordvalues=dummy_hrs(1:hr_in_month))
#endif

      call add_coord(cdl_dd,'hour',hr_in_day,
     &     coordvalues=dummy_hrs(1:hr_in_day))

      return
      end subroutine diurn_defs

      subroutine get_zstr(lgrid,zstr)
      use diag_com, only : ctr_ml,edg_ml,ctr_cp,edg_cp
      implicit none
      integer, intent(in) :: lgrid
      character(len=*), intent(out) :: zstr
      zstr=''
      select case(lgrid)
      case (ctr_ml,edg_ml)
        zstr='level'
      case(ctr_cp)
        zstr='plm'
      case(edg_cp)
        zstr='ple'
      end select
      return
      end subroutine get_zstr
