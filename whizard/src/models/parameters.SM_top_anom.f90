! parameters.SM_top_anom.f90
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module parameters_sm_top_anom
  use kinds
  use constants
  use sm_physics !NODEP!
  implicit none
  private

  real(default), dimension(27), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs, ig, unit, half

  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn
  real(default), public :: ttop, tbot, tch, ttau, tw
  real(default), public :: ltop, lbot, lc, ltau, lw
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, ghmm, &
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww, &
       gvl_qbub_n, gvl_qw, gvl_qw_u, gvl_qw_d, &
       gsl_dttr, gsr_dttr, gsl_dttl, gsl_dbtl, &
       c_quqd1_1, c_quqd1_2, c_quqd8_1, c_quqd8_2
  real(default), public :: vev, lambda, gi_flag, norm_flag, norm, &
       n_tvaa, n_vlrz, n_tvaz, n_vlrw, n_tlrw, n_tvag, n_sph
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn, &
       tvaa, tvaabb, vlrz, vlrcz, vlruz, tvaz, tcvaz, tuvaz, tvazbb, tcvaa, tuvaa, &
       vlrw, tlrw, tvag, tcvag, tuvag, sph, &
       gvlr_qbub, gvlr_qbub_u, gvlr_qbub_d, gvlr_qbub_e, &
       gvlr_qgug, gslr_dbtr
  integer, public :: fun_flag
  logical, public :: bz=.false., bw=.false., ba=.false.

  public :: import_from_whizard, model_update_alpha_s, &
       gmom, gtva_tta, gtva_tca, gtva_tua, gvlr_ttz, gvlr_tcz, gvlr_tuz, &
       gtva_ttz, gtva_tcz, gtva_tuz, gvlr_btw, gvlr_tbw, &
       gtlr_btw, gtrl_tbw, gtlr_btwz, gtrl_tbwz, gtlr_btwa, gtrl_tbwa, &
       gtva_ttww, gtva_bba, gtva_bbz, gtva_bbww, &
       gtva_ttg, gtva_ttgg, gtva_tcg, gtva_tug, gtva_tcgg, gtva_tugg, gsp_tth

contains

  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(73), intent(in) :: par_array
    integer, intent(in) :: scheme
    type :: parameter_set
       real(default) :: gf
       real(default) :: mZ
       real(default) :: mW
       real(default) :: mH
       real(default) :: alphas
       real(default) :: me
       real(default) :: mmu
       real(default) :: mtau
       real(default) :: ms
       real(default) :: mc
       real(default) :: mb
       real(default) :: mtop
       real(default) :: wtop
       real(default) :: wZ
       real(default) :: wW
       real(default) :: wH
       real(default) :: khgaz
       real(default) :: khgaga
       real(default) :: khgg
       real(default) :: xi0
       real(default) :: xipm
       real(default) :: tvA
       real(default) :: taA
       real(default) :: tcvA
       real(default) :: tcaA
       real(default) :: tuvA
       real(default) :: tuaA
       real(default) :: vlZ
       real(default) :: vrZ
       real(default) :: vlcZ
       real(default) :: vrcZ
       real(default) :: vluZ
       real(default) :: vruZ
       real(default) :: tvZ
       real(default) :: taZ
       real(default) :: tcvZ
       real(default) :: tcaZ
       real(default) :: tuvZ
       real(default) :: tuaZ
       real(default) :: vlWRe
       real(default) :: vlWIm
       real(default) :: vrWRe
       real(default) :: vrWIm
       real(default) :: tlWRe
       real(default) :: tlWIm
       real(default) :: trWRe
       real(default) :: trWIm
       real(default) :: tvG
       real(default) :: taG
       real(default) :: tcvG
       real(default) :: tcaG
       real(default) :: tuvG
       real(default) :: tuaG
       real(default) :: sH
       real(default) :: pH
       real(default) :: lam
       real(default) :: fun
       real(default) :: nrm
       real(default) :: gi
       real(default) :: re_CqW
       real(default) :: re_Cquqd1_1
       real(default) :: im_Cquqd1_1
       real(default) :: re_Cquqd1_2
       real(default) :: im_Cquqd1_2
       real(default) :: re_Cquqd8_1
       real(default) :: im_Cquqd8_1
       real(default) :: re_Cquqd8_2
       real(default) :: im_Cquqd8_2
       real(default) :: Rt
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    !!! This corresponds to 1/alpha = 137.03598949333
    real(default), parameter :: &
         alpha = 1.0_default/137.03598949333_default
    e_em = sqrt(4.0_default * PI * alpha)
    par%gf     = par_array(1)
    par%mZ     = par_array(2)
    par%mW     = par_array(3)
    par%mH     = par_array(4)
    par%alphas = par_array(5)
    par%me     = par_array(6)
    par%mmu    = par_array(7)
    par%mtau   = par_array(8)
    par%ms     = par_array(9)
    par%mc     = par_array(10)
    par%mb     = par_array(11)
    par%mtop   = par_array(12)
    par%wtop   = par_array(13)
    par%wZ     = par_array(14)
    par%wW     = par_array(15)
    par%wH     = par_array(16)
    par%khgaz  = par_array(17)
    par%khgaga = par_array(18)
    par%khgg   = par_array(19)
    par%xi0    = par_array(20)
    par%xipm   = par_array(21)
    par%tvA    = par_array(22)
    par%taA    = par_array(23)
    par%tcvA   = par_array(24)
    par%tcaA   = par_array(25)
    par%tuvA   = par_array(26)
    par%tuaA   = par_array(27)
    par%vlZ    = par_array(28)
    par%vrZ    = par_array(29)
    par%vlcZ   = par_array(30)
    par%vrcZ   = par_array(31)
    par%vluZ   = par_array(32)
    par%vruZ   = par_array(33)
    par%tvZ    = par_array(34)
    par%taZ    = par_array(35)
    par%tcvZ   = par_array(36)
    par%tcaZ   = par_array(37)
    par%tuvZ   = par_array(38)
    par%tuaZ   = par_array(39)
    par%vlWRe  = par_array(40)
    par%vlWIm  = par_array(41)
    par%vrWRe  = par_array(42)
    par%vrWIm  = par_array(43)
    par%tlWRe  = par_array(44)
    par%tlWIm  = par_array(45)
    par%trWRe  = par_array(46)
    par%trWIm  = par_array(47)
    par%tvG    = par_array(48)
    par%taG    = par_array(49)
    par%tcvG   = par_array(50)
    par%tcaG   = par_array(51)
    par%tuvG   = par_array(52)
    par%tuaG   = par_array(53)
    par%sH     = par_array(54)
    par%pH     = par_array(55)
    par%lam    = par_array(56)
    par%fun    = par_array(57)
    par%nrm    = par_array(58)
    par%gi     = par_array(59)
    par%re_CqW = par_array(60)
    par%re_Cquqd1_1 = par_array(61)
    par%im_Cquqd1_1 = par_array(62)
    par%re_Cquqd1_2 = par_array(63)
    par%im_Cquqd1_2 = par_array(64)
    par%re_Cquqd8_1 = par_array(65)
    par%im_Cquqd8_1 = par_array(66)
    par%re_Cquqd8_2 = par_array(67)
    par%im_Cquqd8_2 = par_array(68)
    par%Rt     = par_array(69)
    par%v      = par_array(70)
    par%cw     = par_array(71)
    par%sw     = par_array(72)
    par%ee     = par_array(73)
    mass(1:27) = 0
    width(1:27) = 0
    mass(3) = par%ms
    mass(4) = par%mc
    mass(5) = par%mb
    mass(6) = par%mtop
    width(6) = par%wtop
    mass(11) = par%me
    mass(13) = par%mmu
    mass(15) = par%mtau
    mass(23) = par%mZ
    width(23) = par%wZ
    mass(24) = par%mW
    width(24) = par%wW
    mass(25) = par%mH
    width(25) = par%wH
    mass(26) =  par%xi0 * mass(23)
    width(26) =  0
    mass(27) =  par%xipm * mass(24)
    width(27) =  0
    ttop = 4.0_default * mass(6)**2 / mass(25)**2
    tbot = 4.0_default * mass(5)**2 / mass(25)**2
    tch  = 4.0_default * mass(4)**2 / mass(25)**2
    ttau = 4.0_default * mass(15)**2 / mass(25)**2
    tw   = 4.0_default * mass(24)**2 / mass(25)**2
    ltop = 4.0_default * mass(6)**2 / mass(23)**2
    lbot = 4.0_default * mass(5)**2 / mass(23)**2
    lc   = 4.0_default * mass(4)**2 / mass(23)**2
    ltau = 4.0_default * mass(15)**2 / mass(23)**2
    lw   = 4.0_default * mass(24)**2 / mass(23)**2
    vev = par%v
    e = par%ee
    sinthw = par%sw
    sin2thw = par%sw**2
    costhw = par%cw
    tanthw = sinthw/costhw
    qelep = - 1
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
    ig = cmplx (0.0_default, 1.0_default, kind=default) * g
    unit = 1.0_default
    half = 0.5_default
    gcc = - g / 2 / sqrt (2.0_default)
    gncneu(1) = - g / 2 / costhw * ( + 0.5_default)
    gnclep(1) = - g / 2 / costhw * ( - 0.5_default - 2 * qelep * sin2thw)
    gncup(1)  = - g / 2 / costhw * ( + 0.5_default - 2 * qeup  * sin2thw)
    gncdwn(1) = - g / 2 / costhw * ( - 0.5_default - 2 * qedwn * sin2thw)
    gncneu(2) = - g / 2 / costhw * ( + 0.5_default)
    gnclep(2) = - g / 2 / costhw * ( - 0.5_default)
    gncup(2)  = - g / 2 / costhw * ( + 0.5_default)
    gncdwn(2) = - g / 2 / costhw * ( - 0.5_default)
    qlep = - e * qelep
    qup = - e * qeup
    qdwn = - e * qedwn
    qw = e
    iqw = imago * qw
    gzww = g * costhw
    igzww = imago * gzww
    gwww = g
    igwww = imago * gwww
    gw4 = gwww**2
    gzzww = gzww**2
    gazww = gzww * qw
    gaaww = qw**2
    ghww = mass(24) * g
    ghhww = g**2 / 2.0_default
    ghzz = mass(23) * g / costhw
    ghhzz = g**2 / 2.0_default / costhw**2
    ghtt = - mass(6) / vev
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    ghmm = - mass(13) / vev
    gh3 = - 3 * mass(25)**2 / vev
    gh4 = - 3 * mass(25)**2 / vev**2
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs
    fun_flag = nint(par%fun)

    lambda = par%lam
    norm_flag = par%nrm
    norm = vev / lambda**2
    if ( norm_flag > 0. ) then
      n_tvaa  = e / mass(6)
      n_vlrz  = g / 2 / costhw / 2
      n_tvaz  = 1.0_default / vev
      n_vlrw  = g / sqrt(2.0_default) / 2
      n_tlrw  = g / sqrt(2.0_default) / mass(24) / 2
      n_tvag  = gs / mass(6)
      n_sph   = sqrt(0.5_default)
    else
      n_tvaa = e * norm
      n_vlrz = mass(23) * norm / 2
      n_tvaz = norm
      n_vlrw = sqrt(2.0_default) * mass(24) * norm / 2
      n_tlrw = sqrt(2.0_default) * norm / 2
      n_tvag = gs * norm
      n_sph  = sqrt(0.5_default) * vev * norm
    end if

    tvaa(1)  = n_tvaa * par%tvA
    tvaa(2)  = n_tvaa * par%taA * imago
    tcvaa(1) = n_tvaa * par%tcvA
    tcvaa(2) = n_tvaa * par%tcaA * imago
    tuvaa(1) = n_tvaa * par%tuvA
    tuvaa(2) = n_tvaa * par%tuaA * imago
    vlrz(1)  = n_vlrz * par%vlZ
    vlrz(2)  = n_vlrz * par%vrZ
    vlrcz(1) = n_vlrz * par%vlcZ
    vlrcz(2) = n_vlrz * par%vrcZ
    vlruz(1) = n_vlrz * par%vluZ
    vlruz(2) = n_vlrz * par%vruZ
    tvaz(1)  = n_tvaz * par%tvZ
    tvaz(2)  = n_tvaz * par%taZ * imago
    tcvaz(1) = n_tvaz * par%tcvZ
    tcvaz(2) = n_tvaz * par%tcaZ * imago
    tuvaz(1) = n_tvaz * par%tuvZ
    tuvaz(2) = n_tvaz * par%tuaZ * imago
    vlrw(1)  = n_vlrw * ( par%vlWRe + par%vlWIm * imago )
    vlrw(2)  = n_vlrw * ( par%vrWRe + par%vrWIm * imago )
    tlrw(1)  = n_tlrw * ( par%tlWRe + par%tlWIm * imago )
    tlrw(2)  = n_tlrw * ( par%trWRe + par%trWIm * imago )
    tvag(1)  = n_tvag * par%tvG
    tvag(2)  = n_tvag * par%taG * imago
    tcvag(1) = n_tvag * par%tcvG
    tcvag(2) = n_tvag * par%tcaG * imago
    tuvag(1) = n_tvag * par%tuvG
    tuvag(2) = n_tvag * par%tuaG * imago
    sph(1)   = n_sph  * par%sH
    sph(2)   = n_sph  * par%pH  * imago
    tvaabb(1) = 0.0_default
    tvaabb(2) = 0.0_default
    tvazbb(1) = 0.0_default
    tvazbb(2) = 0.0_default

    gi_flag = par%gi
    if ( gi_flag > 0. ) then

      if ( abs(par%vlWRe) > 0. ) then
        print *, "WARNING: gauge invariance and vanishing anomalous"
        print *, "  bbZ vector coupling implies a relation"
        print *, "           vl_ttZ ~ vl_tbW_Re ."
        print *, "  Inferring vl_ttZ from vl_tbW_Re and IGNORING any"
        print *, "  inconsistent values set!"

        vlrz(1) = real(vlrw(1)) * sqrt(2.0_default) / costhw

        print *, "  vl_ttZ = ", real(vlrz(1))/n_vlrz
      end if

      if ( ( abs(par%tlWRe) > 0. ).or.( abs(par%tlWIm) > 0. ) ) then
        print *, "WARNING: anomalous tbW tensor couplings are related to anomalous"
        print *, "  bbZ and bbA tensor couplings by gauge invariance:"
        print *, "  Inferring bottom couplings from tl_tbW"

        tvaabb(1) = real(tlrw(1))        * sinthw / sqrt(2.0_default)
        tvaabb(2) = aimag(tlrw(1))*imago * sinthw / sqrt(2.0_default)
        tvazbb(1) = real(tlrw(1))        * costhw / sqrt(2.0_default)
        tvazbb(2) = aimag(tlrw(1))*imago * costhw / sqrt(2.0_default)

        print *, "  tv_bbA = ", real(tvaabb(1))/n_tvaa
        print *, "  ta_bbA = ", aimag(tvaabb(2))/n_tvaa
        print *, "  tv_bbZ = ", real(tvazbb(1))/n_tvaz
        print *, "  ta_bbZ = ", aimag(tvazbb(2))/n_tvaz
      end if

      if ( ( abs(par%tvZ) > 0. ).or.( abs(par%taZ) > 0. ) ) then
        bz = .true.
      end if
      if ( ( abs(par%trWRe) > 0. ).or.( abs(par%trWIm) > 0. ) ) then
        bw = .true.
      end if
      if ( ( abs(par%tvA) > 0. ).or.( abs(par%taA)> 0. ) ) then
        ba = .true.
      end if

      if ( bz.or.bw.or.ba ) then
        print *, "WARNING: top anomalous tensor couplings to W, A and Z"
        print *, "  are interrelated by gauge invariance:"
        print *, "  Inferring Z couplings from W/A couplings according to"
        print *, "  the relation in the model file and IGNORING any inconsistent"
        print *, "  values set manually! (Exception: only tX_ttZ != 0:"
        print *, "  tr_tbW ~ tv_ttZ + i*ta_ttZ and tX_ttA = 0)"

        if ( ( bw.and.bz ).and..not.ba ) then
          tvaa(1) = ( real(tlrw(2))        - sqrt(2.0_default)*costhw*tvaz(1) ) / ( 2.0_default*sinthw )
          tvaa(2) = ( aimag(tlrw(2))*imago - sqrt(2.0_default)*costhw*tvaz(2) ) / ( 2.0_default*sinthw )
        else if ( bz.and..not.bw ) then
          tlrw(2) = sqrt(2.0_default)*costhw*( tvaz(1) + tvaz(2) ) + 2.0_default*sinthw*( tvaa(1) + tvaa(2) )
        else
          tvaz(1) = ( real(tlrw(2))        - 2.0_default*sinthw*tvaa(1) ) / ( sqrt(2.0_default)*costhw )
          tvaz(2) = ( aimag(tlrw(2))*imago - 2.0_default*sinthw*tvaa(2) ) / ( sqrt(2.0_default)*costhw )
        end if

        print *, "  tv_ttA = ", real(tvaa(1))/n_tvaa
        print *, "  ta_ttA = ", aimag(tvaa(2))/n_tvaa
        print *, "  tv_ttZ = ", real(tvaz(1))/n_tvaz
        print *, "  ta_ttZ = ", aimag(tvaz(2))/n_tvaz
        print *, "  tr_tbW = ", real(tlrw(2))/n_tlrw, " + ", aimag(tlrw(2))/n_tlrw, "I"
      end if
    end if

    !!! Contact interactions: don't forget the i^2 = -1 factors from the internal handling!
    ! gvlr_qgug(1)   = - gs / sqrt(2) / 2 / (2*lambda**2) * par%re_CqG  ! *i^2
    ! gvlr_qgug(2)   = - gs / sqrt(2) / 2 / (2*lambda**2) * par%re_CuG  ! *i^2
    ! gvlr_qbub(1)   = - g  / 2 * tanthw / (2*lambda**2) * par%re_CqB  ! *i^2
    ! gvlr_qbub(2)   = - g  / 2 *tanthw / (2*lambda**2) * par%re_CuB  ! *i^2
    gvlr_qgug(1)   =   gs**2 / 48 * (par%Rt/1000.)**2  ! *i^2
    gvlr_qgug(2)   =   gs**2 / 48 * (par%Rt/1000.)**2  ! *i^2
    gvlr_qbub(1)   =   0
    gvlr_qbub(2)   =   0
    gvlr_qbub_u(1) =   1.0_default / 3.0_default / 2
    gvlr_qbub_u(2) =   4.0_default / 3.0_default / 2
    gvlr_qbub_d(1) =   1.0_default / 3.0_default / 2
    gvlr_qbub_d(2) = - 2.0_default / 3.0_default / 2
    gvlr_qbub_e(1) = - 1.0_default / 2
    gvlr_qbub_e(2) = - 2.0_default / 2
    gvl_qbub_n     = - 1.0_default / 2
    !!! C_qW normalization
    ! gvl_qw         = - g / 2 / (lambda**2) * par%re_CqW  ! *i^2
    !!! C_qq3 normalization
    ! gvl_qw         = - 4.0_default / 2 / (lambda**2) * par%re_CqW  ! *i^2
    !!! v3_4f normalization
    gvl_qw         = - 1.0_default / 2 / (lambda**2) * par%re_CqW  ! *i^2
    gvl_qw_u       =   0.5_default / 2
    gvl_qw_d       = - 0.5_default / 2
    ! gsl_dttr       = - mass(6) / 2 / (sqrt(2.0_default)*vev*lambda**2) * 2.0_default * par%re_CDu         ! *i^2
    ! gsr_dttr       = - mass(6) / 2 / (sqrt(2.0_default)*vev*lambda**2) * (par%re_CDd + par%im_CDd*imago)  ! *i^2
    ! gsl_dttl       = - mass(6) / 2 / (sqrt(2.0_default)*vev*lambda**2) * (par%re_CDd - par%im_CDd*imago)  ! *i^2
    ! gslr_dbtr(1)   = - mass(6) / 2 / (sqrt(2.0_default)*vev*lambda**2) * 2.0_default * par%re_CDu         ! *i^2
    ! gslr_dbtr(2)   =   mass(6) / 2 / (sqrt(2.0_default)*vev*lambda**2) * (par%re_CDd + par%im_CDd*imago)  ! *i^2
    ! gsl_dbtl       =   mass(6) / 2 / (sqrt(2.0_default)*vev*lambda**2) * (par%re_CDd - par%im_CDd*imago)  ! *i^2
    gsl_dttr       =   0
    gsr_dttr       =   0
    gsl_dttl       =   0
    gslr_dbtr(1)   =   0
    gslr_dbtr(2)   =   0
    gsl_dbtl       =   0
    c_quqd1_1      =   1.0_default / 2 / (lambda**2) * (par%re_Cquqd1_1 + par%im_Cquqd1_1*imago)  ! *i^2
    c_quqd1_2      =   1.0_default / 2 / (lambda**2) * (par%re_Cquqd1_2 + par%im_Cquqd1_2*imago)  ! *i^2
    !!! Color flow basis, divide by (sqrt(2))**2
    c_quqd8_1      =   1.0_default / 4 / (lambda**2) * (par%re_Cquqd8_1 + par%im_Cquqd8_1*imago)  ! *i^2
    c_quqd8_2      =   1.0_default / 4 / (lambda**2) * (par%re_Cquqd8_2 + par%im_Cquqd8_2*imago)  ! *i^2

  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs
  end subroutine model_update_alpha_s

  pure function gmom (k2, i, coeff, lam) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    complex(default), dimension(2), intent(in) :: coeff
    real(default), intent(in) :: lam
    select case (fun_flag)
      case (0)
        c = coeff(i)
      case (1)
        c = coeff(i) / (1.0_default + k2/lam**2)
      case (2)
        c = coeff(i) / (1.0_default + k2/lam**2)**2
      case (3)
        c = coeff(i) * exp(-k2/lam**2)
    end select
  end function gmom

  pure function gtva_tta (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tvaa, lambda)
  end function gtva_tta

  pure function gtva_tca (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tcvaa, lambda)
  end function gtva_tca

  pure function gtva_tua (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tuvaa, lambda)
  end function gtva_tua

  pure function gtva_bba (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tvaabb, lambda)
  end function gtva_bba

  pure function gvlr_ttz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, vlrz, lambda)
  end function gvlr_ttz

  pure function gvlr_tcz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, vlrcz, lambda)
  end function gvlr_tcz

  pure function gvlr_tuz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, vlruz, lambda)
  end function gvlr_tuz

  pure function gtva_ttz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tvaz, lambda)
  end function gtva_ttz

  pure function gtva_tcz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tcvaz, lambda)
  end function gtva_tcz

  pure function gtva_tuz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tuvaz, lambda)
  end function gtva_tuz

  pure function gtva_bbz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tvazbb, lambda)
  end function gtva_bbz

  pure function gvlr_btw (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, vlrw, lambda)
  end function gvlr_btw

  pure function gvlr_tbw (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, conjg(vlrw), lambda)
  end function gvlr_tbw

  pure function gtlr_btw (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tlrw, lambda)
  end function gtlr_btw

  pure function gtrl_tbw (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, conjg(tlrw), lambda)
  end function gtrl_tbw

  pure function gtlr_btwa (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - sinthw * gtlr_btw(k2, i)
  end function gtlr_btwa

  pure function gtrl_tbwa (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - sinthw * gtrl_tbw(k2, i)
  end function gtrl_tbwa

  pure function gtlr_btwz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - costhw * gtlr_btw(k2, i)
  end function gtlr_btwz

  pure function gtrl_tbwz (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - costhw * gtrl_tbw(k2, i)
  end function gtrl_tbwz

  pure function gtva_ttww (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! additional factor 2 wrt. ward relation to cancel explicit 1/2 in "gtlr"
    select case (i)
      case (1)
        c = - sqrt(2.0_default) * real(gtlr_btw(k2, 2))
      case (2)
        c = - sqrt(2.0_default) * aimag(gtlr_btw(k2, 2)) * imago
    end select
  end function gtva_ttww

  pure function gtva_bbww (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! additional factor 2 wrt. ward relation to cancel explicit 1/2 in "gtlr"
    select case (i)
      case (1)
        c = - sqrt(2.0_default) * real(gtlr_btw(k2, 1))
      case (2)
        c =   sqrt(2.0_default) * aimag(gtlr_btw(k2, 1)) * imago
    end select
  end function gtva_bbww

  pure function gtva_ttg (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tvag, lambda)
  end function gtva_ttg

  pure function gtva_ttgg (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - gtva_ttg(k2, i)
  end function gtva_ttgg

  pure function gtva_tcg (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tcvag, lambda)
  end function gtva_tcg

  pure function gtva_tug (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, tuvag, lambda)
  end function gtva_tug

  pure function gtva_tcgg (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - gtva_tcg(k2, i)
  end function gtva_tcgg

  pure function gtva_tugg (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    !!! don't touch this relative factor: fixed by ward identity!
    c = - gtva_tug(k2, i)
  end function gtva_tugg

  pure function gsp_tth (k2, i) result (c)
    complex(default) :: c
    real(default), intent(in) :: k2
    integer, intent(in) :: i
    c = - gmom (k2, i, sph, lambda)
  end function gsp_tth

end module parameters_sm_top_anom
