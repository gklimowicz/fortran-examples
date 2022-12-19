! parameters.SM_top_anom.f90 --
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
       tvaa, tvaabb, vlrz, vlrcz, vlruz, tvaz, tcvaz, tvazbb, tcvaa, tuvaa, &
       tuvaz, vlrw, tlrw, tvag, tcvag, tuvag, sph, &
       gvlr_qbub, gvlr_qbub_u, gvlr_qbub_d, gvlr_qbub_e, &
       gvlr_qgug, gslr_dbtr
  integer, public :: fun_flag
  logical, public :: bz=.false., bw=.false., ba=.false.

  public :: init_parameters, model_update_alpha_s, &
       gmom, gtva_tta, gtva_tca, gtva_tua, gtva_tuz, gvlr_tuz, &
       gvlr_ttz, gvlr_tcz, gtva_ttz, gtva_tcz, gvlr_btw, gvlr_tbw, &
       gtlr_btw, gtrl_tbw, gtlr_btwz, gtrl_tbwz, gtlr_btwa, gtrl_tbwa, &
       gtva_ttww, gtva_bba, gtva_bbz, gtva_bbww, &
       gtva_ttg, gtva_ttgg, gtva_tcg, gtva_tug, gtva_tcgg, gtva_tugg, gsp_tth

  real(default), parameter :: &
          GF = 1.16639E-5_default   ! Fermi constant
  !!! This corresponds to 1/alpha = 137.03598949333
  real(default), parameter :: &
          alpha = 1.0_default/137.03598949333_default
  complex(default), parameter :: &
          alphas =  0.1178_default  ! Strong coupling constant (Z point)

contains

  subroutine init_parameters

    tvaa(1)     = 1
    tvaa(2)     = 1 * imago
    tcvaa(1)    = 1
    tcvaa(2)    = 1 * imago
    tuvaa(1)    = 1
    tuvaa(2)    = 1 * imago
    tuvaz(1)    = 1
    tuvaz(2)    = 1 * imago
    vlrz(1)     = 1
    vlrz(2)     = 1
    vlrcz(1)    = 1
    vlrcz(2)    = 1 * imago
    vlruz(1)    = 1
    vlruz(2)    = 1 * imago
    tvaz(1)     = 1
    tvaz(2)     = 1 * imago
    tcvaz(1)    = 1
    tcvaz(2)    = 1 * imago
    vlrw(1)     = 1 + 1 * imago
    vlrw(2)     = 1 + 1 * imago
    tlrw(1)     = 1 + 1 * imago
    tlrw(2)     = 1 + 1 * imago
    tvag(1)     = 1
    tvag(2)     = 1 * imago
    tcvag(1)    = 1
    tcvag(2)    = 1 * imago
    tuvag(1)    = 1
    tuvag(2)    = 1 * imago
    sph(1)      = 0
    sph(2)      = 0
    lambda      = 2000
    fun_flag    = 0
!     norm_flag   = 1
    gi_flag     = 1

    mass(1:27)  = 0
    width(1:27) = 0
!     mass(3)     = 0.095_default        ! s-quark mass 
!     mass(4)     = 1.2_default          ! c-quark mass
     mass(5)     = 4.2_default          ! b-quark mass
     mass(6)     = 173.1_default        ! t-quark mass
!     width(6)    = 1.523_default        ! t-quark width
!     mass(11)    = 0.000510997_default  ! electron mass
!     mass(13)    = 0.105658389_default  ! muon mass
!     mass(15)    = 1.77705_default      ! tau-lepton mas
    mass(23)    = 91.1882_default      ! Z-boson mass
    width(23)   = 2.443_default        ! Z-boson width
    mass(24)    = 80.419_default       ! W-boson mass
    width(24)   = 2.049_default        ! W-boson width
    mass(25)    = 200._default          ! Higgs mass
    width(25)   = 1.419_default        ! Higgs width

    ttop = 4.0_default * mass(6)**2 / mass(25)**2
    tbot = 4.0_default * mass(5)**2 / mass(25)**2
    tch  = 4.0_default * mass(4)**2 / mass(25)**2
    ttau = 4.0_default * mass(15)**2 / mass(25)**2
    tw   = 4.0_default * mass(24)**2 / mass(25)**2  
!     ltop = 4.0_default * mass(6)**2 / mass(23)**2
!     lbot = 4.0_default * mass(5)**2 / mass(23)**2  
!     lc   = 4.0_default * mass(4)**2 / mass(23)**2
!     ltau = 4.0_default * mass(15)**2 / mass(23)**2
!     lw   = 4.0_default * mass(24)**2 / mass(23)**2

    e_em = sqrt(4.0_default * PI * alpha)
    vev = 1 / sqrt (sqrt (2.0_default) * GF)  ! v (Higgs vev)
    ! costhw = mass(24) / mass(23)  ! cos(theta-W)
    costhw = 0.881901_default
    sinthw = sqrt (1.0_default-costhw**2)  ! sin(theta-W)
    sin2thw = sinthw**2
    tanthw = sinthw/costhw
    ! e = 2.0_default * sinthw * mass(24) / vev  ! em-coupling (GF scheme)
    e = 0.349196_default

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
    iqw = imago*qw
    gzww = g * costhw
    igzww = imago*gzww
    gwww = g
    igwww = imago*gwww
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
    gs = sqrt(2.0_default*PI*alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs

    norm = vev / lambda**2
    n_tvaa = e * norm
    n_vlrz = mass(23) * norm / 2
    n_tvaz = norm
    n_vlrw = sqrt(2.0_default) * mass(24) * norm / 2
    n_tlrw = sqrt(2.0_default) * norm / 2
    ! Assuming that imag(gs).eq.0 ...
    n_tvag = real (gs, kind=default) * norm
    n_sph  = sqrt(0.5_default) * vev * norm

    tvaabb(1) = 0.0_default
    tvaabb(2) = 0.0_default
    tvazbb(1) = 0.0_default
    tvazbb(2) = 0.0_default

    if ( gi_flag > 0. ) then

      if ( abs(real(vlrw(1))) > 0. ) then
        print *, "WARNING: gauge invariance and vanishing anomalous"
        print *, "  bbZ vector coupling implies a relation"
        print *, "           vl_ttZ ~ vl_tbW_Re ."
        print *, "  Inferring vl_ttZ from vl_tbW_Re and IGNORING any"
        print *, "  inconsistent values set!"

        vlrz(1) = real(vlrw(1)) * sqrt(2.0_default) / costhw

        print *, "  vl_ttZ = ", real(vlrz(1))/n_vlrz
      end if

      if ( ( abs(real(tlrw(1))) > 0. ).or.( abs(aimag(tlrw(1))) > 0. ) ) then
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

      if ( ( abs(real(tvaz(1))) > 0. ).or.( abs(aimag(tvaz(2))) > 0. ) ) then
        bz = .true.
      end if
      if ( ( abs(real(tlrw(2))) > 0. ).or.( abs(real(tlrw(2))) > 0. ) ) then
        bw = .true.
      end if
      if ( ( abs(real(tvaa(1))) > 0. ).or.( abs(aimag(tvaa(2)))> 0. ) ) then
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

    gvlr_qgug(1)   =   0
    gvlr_qgug(2)   =   0
    gvlr_qbub(1)   =   0
    gvlr_qbub(2)   =   0
    gvlr_qbub_u(1) =   1.0_default / 3.0_default / 2
    gvlr_qbub_u(2) =   4.0_default / 3.0_default / 2
    gvlr_qbub_d(1) =   1.0_default / 3.0_default / 2
    gvlr_qbub_d(2) = - 2.0_default / 3.0_default / 2
    gvlr_qbub_e(1) = - 1.0_default / 2
    gvlr_qbub_e(2) = - 2.0_default / 2
    gvl_qbub_n     = - 1.0_default / 2
    gvl_qw         =   0
    gvl_qw_u       =   0.5_default / 2
    gvl_qw_d       = - 0.5_default / 2
    gsl_dttr       =   0
    gsr_dttr       =   0
    gsl_dttl       =   0
    gslr_dbtr(1)   =   0
    gslr_dbtr(2)   =   0
    gsl_dbtl       =   0
    c_quqd1_1      =   0
    c_quqd1_2      =   0
    c_quqd8_1      =   0
    c_quqd8_2      =   0

  end subroutine init_parameters

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
      case default
        c = 0 ! hoping we never get here ...
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
      case default
        c = 0 ! hoping we never get here ...
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
      case default
        c = 0 ! hoping we never get here ...
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
