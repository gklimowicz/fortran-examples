! parameters.MSSM.omega.f90
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
module parameters_mssm_grav
  use kinds
  use constants
  implicit none
  private
  public :: import_from_whizard, model_update_alpha_s
  real(kind=default), dimension(70), save, public :: mass = 0, width = 0
  real(kind=default), parameter, public :: GeV = 1.0_default
  real(kind=default), parameter, public :: MeV = GeV / 1000
  real(kind=default), parameter, public :: keV = MeV / 1000
  real(kind=default), parameter, public :: TeV = GeV * 1000
  real(kind=default), save, public :: &
       alpha = 1.0_default / 137.0359895_default, &
       sin2thw = 0.23124_default
  integer, save, public :: &
       sign1 = +1, sign2 = +1, sign3 = +1, sign4 = +1
  real(kind=default), save, public :: &
       sigch1 = +1, sigch2 = +1 
  complex(kind=default), save, private :: vev
  real(kind=default), public, save :: sind = 0._default, & 
    cosd = 1._default, sinckm12 = 0.223_default, &
    sinckm13 = 0.004_default, sinckm23 = 0.04_default, &
    tana = 30._default, tanb = 30._default, as = 0._default
  real(kind=default), public, save :: dummy1, cos2am2b, sin2am2b, sinamb, & 
    sinapb, cosamb, cosapb, cos4be, sin4be, sin4al, sin2al, sin2be, cos2al, &
    cos2be, cosbe, sinbe, cosal, sinal, costhw, sinthw 
  real(kind=default), public, save :: q_lep, q_up, q_down 
  complex(kind=default), public, save :: gcc, qchar, qdwn, qup, qlep, & 
        gz, g, e, gs
  complex(kind=default), save, public :: xia = 1, xi0 = 1, xipm = 1
  complex(kind=default), dimension(2), public, save :: gncdwn
  complex(kind=default), dimension(2), public, save :: gncup
  complex(kind=default), dimension(2), public, save :: gnclep
  complex(kind=default), dimension(2), public, save :: gncneu
  complex(kind=default), public, save :: g_yuk_ch1_sn1_2_c, & 
    g_yuk_ch1_sn1_2, g_yuk_ch1_sn1_1_c, g_yuk_ch1_sn1_1, g_yuk_ch2_sn1_2_c, &
    g_yuk_ch2_sn1_2, g_yuk_ch2_sn1_1_c, g_yuk_ch2_sn1_1
  complex(kind=default), public, save :: g_yuk_ch2_su1_1_2_c, &
    g_yuk_ch2_su1_1_2, g_yuk_ch2_sd1_1_2_c, g_yuk_ch2_sd1_1_2, &
    g_yuk_ch1_su1_1_2_c, g_yuk_ch1_su1_1_2, g_yuk_ch1_sd1_1_2_c, &
    g_yuk_ch1_sd1_1_2, g_yuk_ch2_su1_1_1_c, g_yuk_ch2_su1_1_1, &
    g_yuk_ch2_sd1_1_1_c, g_yuk_ch2_sd1_1_1, g_yuk_ch1_su1_1_1_c, & 
    g_yuk_ch1_su1_1_1, g_yuk_ch1_sd1_1_1_c, g_yuk_ch1_sd1_1_1, &
    g_yuk_ch2_su1_2_2_c, g_yuk_ch2_su1_2_2, g_yuk_ch2_sd1_2_2_c, &
    g_yuk_ch2_sd1_2_2, g_yuk_ch1_su1_2_2_c, g_yuk_ch1_su1_2_2, & 
    g_yuk_ch1_sd1_2_2_c, g_yuk_ch1_sd1_2_2, g_yuk_ch2_su1_2_1_c, &
    g_yuk_ch2_su1_2_1, g_yuk_ch2_sd1_2_1_c, g_yuk_ch2_sd1_2_1, &
    g_yuk_ch1_su1_2_1_c, g_yuk_ch1_su1_2_1, g_yuk_ch1_sd1_2_1_c, &
    g_yuk_ch1_sd1_2_1
  complex(kind=default), public, save :: g_yuk_n4_sn1_3_c, g_yuk_n4_sn1_3, &
    g_yuk_n4_sn1_2_c, g_yuk_n4_sn1_2, g_yuk_n4_sn1_1_c, g_yuk_n4_sn1_1, &
    g_yuk_n3_sn1_3_c, g_yuk_n3_sn1_3, g_yuk_n3_sn1_2_c, g_yuk_n3_sn1_2, &
    g_yuk_n3_sn1_1_c, g_yuk_n3_sn1_1, g_yuk_n2_sn1_3_c, g_yuk_n2_sn1_3, &
    g_yuk_n2_sn1_2_c, g_yuk_n2_sn1_2, g_yuk_n2_sn1_1_c, g_yuk_n2_sn1_1, &
    g_yuk_n1_sn1_3_c, g_yuk_n1_sn1_3, g_yuk_n1_sn1_2_c, g_yuk_n1_sn1_2, &
    g_yuk_n1_sn1_1_c, g_yuk_n1_sn1_1, g_yuk_ch2_sl2_3_c, g_yuk_ch2_sl2_3, &
    g_yuk_ch2_sl1_3_c, g_yuk_ch2_sl1_3, g_yuk_ch2_sl1_2_c, g_yuk_ch2_sl1_2, &
    g_yuk_ch2_sl1_1_c, g_yuk_ch2_sl1_1, g_yuk_ch1_sl2_3_c, g_yuk_ch1_sl2_3, &
    g_yuk_ch1_sl1_3_c, g_yuk_ch1_sl1_3, g_yuk_ch1_sl1_2_c, g_yuk_ch1_sl1_2, &
    g_yuk_ch1_sl1_1_c, g_yuk_ch1_sl1_1, ghsu2sd2_3_3_c, ghsu2sd2_3_3, &
    ghsu2sd1_3_3_c, ghsu2sd1_3_3, ghsu1sd2_3_3_c, ghsu1sd2_3_3, ghsu1sd1_3_3_c, &
    ghsu1sd1_3_3, ghsu2sd2_3_2_c, ghsu2sd2_3_2, ghsu2sd1_3_2_c, ghsu2sd1_3_2, &
    ghsu1sd2_3_2_c, ghsu1sd2_3_2, ghsu1sd1_3_2_c, ghsu1sd1_3_2, ghsu2sd2_3_1_c, &
    ghsu2sd2_3_1, ghsu2sd1_3_1_c, ghsu2sd1_3_1, ghsu1sd2_3_1_c
  complex(kind=default), public, save :: ghsu1sd2_3_1, ghsu1sd1_3_1_c, &
    ghsu1sd1_3_1, ghsu2sd2_2_3_c, ghsu2sd2_2_3, ghsu2sd1_2_3_c, ghsu2sd1_2_3, &
    ghsu1sd2_2_3_c, ghsu1sd2_2_3, ghsu1sd1_2_3_c, ghsu1sd1_2_3, ghsu1sd1_2_2_c, &
    ghsu1sd1_2_2, ghsu1sd1_2_1_c, ghsu1sd1_2_1, ghsu2sd2_1_3_c, ghsu2sd2_1_3, & 
    ghsu2sd1_1_3_c, ghsu2sd1_1_3, ghsu1sd2_1_3_c, ghsu1sd2_1_3, ghsu1sd1_1_3_c, & 
    ghsu1sd1_1_3, ghsu1sd1_1_2_c, ghsu1sd1_1_2, ghsu1sd1_1_1_c, ghsu1sd1_1_1, & 
    gh2sn1sn1_3, gh1sn1sn1_3, ghsnsl2_3_c, ghsnsl2_3, ghsnsl1_3_c, ghsnsl1_3, &
    gh2sd2sd2_3, gh2su2su2_3, gh2sl2sl2_3, gh1sd2sd2_3, gh1su2su2_3, gh1sl2sl2_3
  complex(kind=default), public, save :: g_yuk_n4_sd2_2_c, g_yuk_n4_sd2_2, &
    g_yuk_n4_su2_2_c, g_yuk_n4_su2_2, g_yuk_n4_sl2_2_c, g_yuk_n4_sl2_2, &
    g_yuk_n3_sd2_2_c, g_yuk_n3_sd2_2, g_yuk_n3_su2_2_c, g_yuk_n3_su2_2, &
    g_yuk_n3_sl2_2_c, g_yuk_n3_sl2_2, g_yuk_n2_sd2_2_c, g_yuk_n2_sd2_2, &
    g_yuk_n2_su2_2_c, g_yuk_n2_su2_2, g_yuk_n2_sl2_2_c, g_yuk_n2_sl2_2, &
    g_yuk_n1_sd2_2_c, g_yuk_n1_sd2_2, g_yuk_n1_su2_2_c, g_yuk_n1_su2_2, &
    g_yuk_n1_sl2_2_c, g_yuk_n1_sl2_2, g_yuk_n4_sd1_2_c, g_yuk_n4_sd1_2, &
    g_yuk_n4_su1_2_c, g_yuk_n4_su1_2, g_yuk_n4_sl1_2_c, g_yuk_n4_sl1_2, &
    g_yuk_n3_sd1_2_c, g_yuk_n3_sd1_2, g_yuk_n3_su1_2_c, g_yuk_n3_su1_2, &
    g_yuk_n3_sl1_2_c, g_yuk_n3_sl1_2, g_yuk_n2_sd1_2_c, g_yuk_n2_sd1_2, &
    g_yuk_n2_su1_2_c, g_yuk_n2_su1_2, g_yuk_n2_sl1_2_c, g_yuk_n2_sl1_2, &
    g_yuk_n1_sd1_2_c, g_yuk_n1_sd1_2, g_yuk_n1_su1_2_c, g_yuk_n1_su1_2, &
    g_yuk_n1_sl1_2_c, g_yuk_n1_sl1_2, g_yuk_n4_sd2_1_c, g_yuk_n4_sd2_1, &
    g_yuk_n4_su2_1_c, g_yuk_n4_su2_1, g_yuk_n4_sl2_1_c, g_yuk_n4_sl2_1, &
    g_yuk_n3_sd2_1_c, g_yuk_n3_sd2_1, g_yuk_n3_su2_1_c, g_yuk_n3_su2_1, &
    g_yuk_n3_sl2_1_c, g_yuk_n3_sl2_1, g_yuk_n2_sd2_1_c, g_yuk_n2_sd2_1, &
    g_yuk_n2_su2_1_c, g_yuk_n2_su2_1, g_yuk_n2_sl2_1_c, g_yuk_n2_sl2_1, &
    g_yuk_n1_sd2_1_c, g_yuk_n1_sd2_1, g_yuk_n1_su2_1_c, g_yuk_n1_su2_1, &
    g_yuk_n1_sl2_1_c, g_yuk_n1_sl2_1, g_yuk_n4_sd1_1_c, g_yuk_n4_sd1_1, &
    g_yuk_n4_su1_1_c, g_yuk_n4_su1_1, g_yuk_n4_sl1_1_c, g_yuk_n4_sl1_1, &
    g_yuk_n3_sd1_1_c, g_yuk_n3_sd1_1, g_yuk_n3_su1_1_c, g_yuk_n3_su1_1, &
    g_yuk_n3_sl1_1_c, g_yuk_n3_sl1_1, g_yuk_n2_sd1_1_c, g_yuk_n2_sd1_1, &
    g_yuk_n2_su1_1_c, g_yuk_n2_su1_1, g_yuk_n2_sl1_1_c, g_yuk_n2_sl1_1, &
    g_yuk_n1_sd1_1_c, g_yuk_n1_sd1_1, g_yuk_n1_su1_1_c, g_yuk_n1_su1_1, &
    g_yuk_n1_sl1_1_c, g_yuk_n1_sl1_1
  complex(kind=default), public, save :: gh2sd2sd1_3, gh2su2su1_3, &
    gh2sl2sl1_3, gh1sd2sd1_3, gh1su2su1_3, gh1sl2sl1_3, &
    gh2sd1sd2_3, gh2su1su2_3, gh2sl1sl2_3, gh1sd1sd2_3, &
    gh1su1su2_3, gh1sl1sl2_3, gh2sd1sd1_3, &
    gh2su1su1_3, gh2sl1sl1_3, gh1sd1sd1_3, gh1su1su1_3, gh1sl1sl1_3, &
    gh2sn1sn1_2, gh1sn1sn1_2, ghsnsl1_2_c, ghsnsl1_2, &
    gh2sd2sd2_2, gh2su2su2_2, gh2sl2sl2_2, &
    gh1sd2sd2_2, gh1su2su2_2, gh1sl2sl2_2, &
    gh2sd1sd1_2, gh2su1su1_2, gh2sl1sl1_2, gh1sd1sd1_2, &
    gh1su1su1_2, gh1sl1sl1_2, gh2sn1sn1_1, gh1sn1sn1_1
  !!! complex(kind=default), public, save :: ghsnsl2_1, ghsnsl2_1_c & 
  !!!   ghsnsl2_2_c, ghsnsl2_2, 
  complex(kind=default), public, save :: ghsnsl1_1_c, ghsnsl1_1, &
    gh2sd2sd2_1, gh2su2su2_1, gh2sl2sl2_1, &
    gh1sd2sd2_1, gh1su2su2_1, gh1sl2sl2_1, &
    gh2sd1sd1_1, gh2su1su1_1, gh2sl1sl1_1, gh1sd1sd1_1, &
    gh1su1su1_1, gh1sl1sl1_1 
  complex(kind=default), public, save :: gasl2sl2_3, gasl2sl1_3, &
    gasl1sl2_3, gasl1sl1_3 !!! , gasl2sl2_2, gasl2sl1_2, gasl1sl2_2, &
    !!! gasl1sl1_2, gasl2sl2_1, gasl2sl1_1, gasl1sl2_1, gasl1sl1_1
  complex(kind=default), public, save :: gasu2su2_3, gasu2su1_3, &
    gasu1su2_3, gasu1su1_3 !!! , gasu2su2_2, gasu2su1_2, gasu1su2_2, &
    !!! gasu1su1_2, gasu2su2_1, gasu2su1_1, gasu1su2_1, gasu1su1_1
  complex(kind=default), public, save :: gasd2sd2_3, gasd2sd1_3, &
    gasd1sd2_3, gasd1sd1_3 !!! , gasd2sd2_2, gasd2sd1_2, gasd1sd2_2, & 
    !!! gasd1sd1_2, gasd2sd2_1, gasd2sd1_1, gasd1sd2_1, gasd1sd1_1 
  complex(kind=default), public, save :: g_h43_321susd, g_h43_312susd, &
    g_h43_322susd, g_h43_311susd, g_h43_221susd, g_h43_212susd, g_h43_222susd, &
    g_h43_211susd, g_h43_121susd
  complex(kind=default), public, save :: g_h43_112susd, g_h43_122susd, &
    g_h43_111susd, g_h42_321susd, g_h42_312susd, g_h42_322susd, g_h42_311susd, &
    g_h42_211susd, g_h42_111susd, g_h41_321susd, g_h41_312susd, g_h41_322susd, & 
    g_h41_311susd, g_h41_211susd, g_h41_111susd, &
    g_h4312slsn, g_h4311slsn, g_h3321slsl, g_h3312slsl, g_h2321slsl, &
    g_h2312slsl, g_h2322slsl, g_h2311slsl, g_h2311snsn, g_h1321slsl, &
    g_h1312slsl, g_h1322slsl, g_h1311slsl, g_h1311snsn, g_h3321sdsd, &
    g_h3312sdsd, g_h3321susu, g_h3312susu, g_h2321sdsd, g_h2312sdsd, &
    g_h2322sdsd, g_h2311sdsd, g_h2321susu, g_h2312susu, g_h2322susu, &
    g_h2311susu, g_h1321sdsd, g_h1312sdsd, g_h1322sdsd, g_h1311sdsd, &
    g_h1321susu, g_h1312susu, g_h1322susu, g_h1311susu, g_h4211slsn, &
    g_h2222slsl, g_h2211slsl
  complex(kind=default), public, save :: g_h2211snsn, &
    g_h1222slsl, g_h1211slsl, g_h1211snsn, g_h2222sdsd, g_h2211sdsd, & 
    g_h2222susu, g_h2211susu, g_h1222sdsd, g_h1211sdsd, g_h1222susu, & 
    g_h1211susu, g_h4111slsn, g_h2122slsl, g_h2111slsl, g_h2111snsn, &
    g_h1122slsl, g_h1111slsl, g_h1111snsn, g_h2122sdsd, g_h2111sdsd, & 
    g_h2122susu, g_h2111susu, g_h1122sdsd, g_h1111sdsd, &
    g_h1122susu, g_h1111susu, gnzn_4_4, gnzn_3_3, gnzn_2_2, &
    gnzn_1_1, rnch_42, lnch_42, rnc_42
  complex(kind=default), public, save :: gcicih1_1_1, gcicih1_2_2, &
    gcicih1_3_3, gcicih1_4_4, gcicih2_1_1, gcicih2_2_2, gcicih2_3_3, &
    gcicih2_4_4, gcicia_1_1, gcicia_2_2, gcicia_3_3, gcicia_4_4 
  !!! complex(kind=default), public, save :: g_h3112susu, g_h3121susu, &
  !!!  g_h3112sdsd, g_h3121sdsd, g_h3112slsl, g_h3121slsl, g_h3212susu, & 
  !!!  g_h3221susu, g_h3212sdsd, g_h3221sdsd, g_h3212slsl, g_h3221slsl, &
  !!! complex(kind=default), public, save :: g_h4112slsn, g_h4212slsn, &
  complex(kind=default), public, save :: lnc_42, rnch_41, &
    lnch_41, rnc_41, lnc_41, rnch_32, lnch_32, rnc_32, lnc_32, rnch_31, & 
    lnch_31, rnc_31, lnc_31, rnch_22, lnch_22, rnc_22, lnc_22, rnch_21, &
    lnch_21, rnc_21, lnc_21, rnch_12, lnch_12, rnc_12, lnc_12, rnch_11, & 
    lnch_11, rnc_11, lnc_11, rcn_24, lcn_24, rcn_23, lcn_23, rcn_22, &
    lcn_22, rcn_21, lcn_21 
  complex(kind=default), public, save :: gch1c_1_1, gch1c_2_2, & 
    gch2c_1_1, gch2c_2_2, gcac_1_1, gcac_2_2            
  complex(kind=default), public, save :: rcn_14, lcn_14, &
    rcn_13, lcn_13, rcn_12, lcn_12, rcn_11, lcn_11, ap_22, vp_22, &
    ap_21, vp_21, ap_12, vp_12, ap_11, vp_11, pnna_44, snna_44, & 
    pnnh2_44, snnh2_44, pnnh1_44
  complex(kind=default), public, save :: snnh1_44, axial0_44, vector0_44, &
    pnna_34, snna_34, pnnh2_34, snnh2_34, pnnh1_34, &
    snnh1_34, axial0_34, vector0_34, pnna_33, snna_33, &
    pnnh2_33, snnh2_33, pnnh1_33, snnh1_33, axial0_33, vector0_33 
  complex(kind=default), public, save :: pnna_24, & 
    snna_24, pnnh2_24, snnh2_24, pnnh1_24, snnh1_24, &
    axial0_24, vector0_24, pnna_23, snna_23, pnnh2_23, snnh2_23, &
    pnnh1_23, snnh1_23, axial0_23, vector0_23, pnna_22, snna_22, &
    pnnh2_22, snnh2_22, pnnh1_22, snnh1_22, axial0_22, vector0_22, &
    pnna_14, snna_14, pnnh2_14, snnh2_14, pnnh1_14, snnh1_14, &
    axial0_14, vector0_14, pnna_13, snna_13, pnnh2_13, snnh2_13, &
    pnnh1_13, snnh1_13, axial0_13, vector0_13, pnna_12, snna_12, &
    pnnh2_12
  complex(kind=default), public, save :: snnh2_12, pnnh1_12, snnh1_12, &
    axial0_12, vector0_12, pnna_11, snna_11, pnnh2_11, snnh2_11, &
    pnnh1_11, snnh1_11, axial0_11, vector0_11, gglwsu2sd1_3_3_c, gglwsu1sd2_3_3_c, &
    gglwsu2sd2_3_3_c, gglwsu1sd1_3_3_c, gglwsu2sd1_3_3, gglwsu1sd2_3_3, &
    gglwsu2sd2_3_3, gglwsu1sd1_3_3, gglwsu2sd1_3_2_c, gglwsu1sd2_3_2_c, &
    gglwsu2sd2_3_2_c, gglwsu1sd1_3_2_c, gglwsu2sd1_3_2, gglwsu1sd2_3_2, &
    gglwsu2sd2_3_2, gglwsu1sd1_3_2, gglwsu2sd1_3_1_c, gglwsu1sd2_3_1_c, &
    gglwsu2sd2_3_1_c, gglwsu1sd1_3_1_c, gglwsu2sd1_3_1, gglwsu1sd2_3_1, &
    gglwsu2sd2_3_1, gglwsu1sd1_3_1, gglwsu2sd1_2_3_c, gglwsu1sd2_2_3_c, &
    gglwsu2sd2_2_3_c, gglwsu1sd1_2_3_c, gglwsu2sd1_2_3, gglwsu1sd2_2_3, &
    gglwsu2sd2_2_3, gglwsu1sd1_2_3, gglwsu2sd1_2_2_c, gglwsu1sd2_2_2_c, &
    gglwsu2sd2_2_2_c, gglwsu1sd1_2_2_c, gglwsu2sd1_2_2, gglwsu1sd2_2_2, &
    gglwsu2sd2_2_2, gglwsu1sd1_2_2, gglwsu2sd1_2_1_c, gglwsu1sd2_2_1_c, &
    gglwsu2sd2_2_1_c, gglwsu1sd1_2_1_c, gglwsu2sd1_2_1, gglwsu1sd2_2_1, &
    gglwsu2sd2_2_1, gglwsu1sd1_2_1, gglwsu2sd1_1_3_c, gglwsu1sd2_1_3_c, &
    gglwsu2sd2_1_3_c, gglwsu1sd1_1_3_c, gglwsu2sd1_1_3, gglwsu1sd2_1_3
  complex(kind=default), public, save :: gglwsu2sd2_1_3, gglwsu1sd1_1_3, &
    gglwsu2sd1_1_2_c, gglwsu1sd2_1_2_c, gglwsu2sd2_1_2_c, gglwsu1sd1_1_2_c, &
    gglwsu2sd1_1_2, gglwsu1sd2_1_2, gglwsu2sd2_1_2, gglwsu1sd1_1_2, &
    gglwsu2sd1_1_1_c, gglwsu1sd2_1_1_c, gglwsu2sd2_1_1_c, gglwsu1sd1_1_1_c, &
    gglwsu2sd1_1_1, gglwsu1sd2_1_1, gglwsu2sd2_1_1, gglwsu1sd1_1_1, mix_sd322, &
    mix_sd321, mix_sd312, mix_sd311, mix_sd222, mix_sd221, mix_sd212, &
    mix_sd211, mix_sd122, mix_sd121, mix_sd112, mix_sd111, mix_su322, &
    mix_su321, mix_su312, mix_su311, mix_su222, mix_su221, mix_su212, &
    mix_su211, mix_su122, mix_su121, mix_su112, mix_su111, mix_sl322, &
    mix_sl321, mix_sl312, mix_sl311, mix_sl222, mix_sl221, mix_sl212, &
    mix_sl211, mix_sl122, mix_sl121, mix_sl112, mix_sl111, gglsd2sd1_3, &
    gglsd1sd2_3, gglsd2sd2_3, gglsd1sd1_3, gglsu2su1_3, gglsu1su2_3, &
    gglsu2su2_3, gglsu1su1_3, gglsd2sd1_2, gglsd1sd2_2, gglsd2sd2_2, &
    gglsd1sd1_2, gglsu2su1_2, gglsu1su2_2, gglsu2su2_2
  complex(kind=default), public, save :: gglsu1su1_2, gglsd2sd1_1, &
    gglsd1sd2_1, gglsd2sd2_1, gglsd1sd1_1, gglsu2su1_1, gglsu1su2_1, &
    gglsu2su2_1, gglsu1su1_1, gglpsqsq, gglglsqsq, gzwpsu2sd1_3_3_c, &
    gzwpsu1sd2_3_3_c, gzwpsu2sd2_3_3_c, gzwpsu1sd1_3_3_c, gzwpsu2sd1_3_3, &
    gzwpsu1sd2_3_3, gzwpsu2sd2_3_3, gzwpsu1sd1_3_3, gpwpsu2sd1_3_3_c, &
    gpwpsu1sd2_3_3_c, gpwpsu2sd2_3_3_c, gpwpsu1sd1_3_3_c, gpwpsu2sd1_3_3, &
    gpwpsu1sd2_3_3, gpwpsu2sd2_3_3, gpwpsu1sd1_3_3, gzwpsu2sd1_3_2_c, &
    gzwpsu1sd2_3_2_c, gzwpsu2sd2_3_2_c, gzwpsu1sd1_3_2_c, gzwpsu2sd1_3_2, &
    gzwpsu1sd2_3_2, gzwpsu2sd2_3_2, gzwpsu1sd1_3_2, gpwpsu2sd1_3_2_c, &
    gpwpsu1sd2_3_2_c, gpwpsu2sd2_3_2_c, gpwpsu1sd1_3_2_c, gpwpsu2sd1_3_2, &
    gpwpsu1sd2_3_2, gpwpsu2sd2_3_2, gpwpsu1sd1_3_2, gzwpsu2sd1_3_1_c, &
    gzwpsu1sd2_3_1_c, gzwpsu2sd2_3_1_c, gzwpsu1sd1_3_1_c, gzwpsu2sd1_3_1, &
    gzwpsu1sd2_3_1, gzwpsu2sd2_3_1, gzwpsu1sd1_3_1, gpwpsu2sd1_3_1_c, &
    gpwpsu1sd2_3_1_c, gpwpsu2sd2_3_1_c, gpwpsu1sd1_3_1_c, gpwpsu2sd1_3_1, &
    gpwpsu1sd2_3_1, gpwpsu2sd2_3_1, gpwpsu1sd1_3_1, gzwpsu2sd1_2_3_c, &
    gzwpsu1sd2_2_3_c, gzwpsu2sd2_2_3_c, gzwpsu1sd1_2_3_c, gzwpsu2sd1_2_3, &
    gzwpsu1sd2_2_3, gzwpsu2sd2_2_3, gzwpsu1sd1_2_3, gpwpsu2sd1_2_3_c, &
    gpwpsu1sd2_2_3_c
  complex(kind=default), public, save :: gpwpsu2sd2_2_3_c, gpwpsu1sd1_2_3_c, &
    gpwpsu2sd1_2_3, gpwpsu1sd2_2_3, gpwpsu2sd2_2_3, gpwpsu1sd1_2_3, &
    gzwpsu2sd1_2_2_c, gzwpsu1sd2_2_2_c, gzwpsu2sd2_2_2_c, gzwpsu1sd1_2_2_c, &
    gzwpsu2sd1_2_2, gzwpsu1sd2_2_2, gzwpsu2sd2_2_2, gzwpsu1sd1_2_2, &
    gpwpsu2sd1_2_2_c, gpwpsu1sd2_2_2_c, gpwpsu2sd2_2_2_c, gpwpsu1sd1_2_2_c, &
    gpwpsu2sd1_2_2, gpwpsu1sd2_2_2, gpwpsu2sd2_2_2, gpwpsu1sd1_2_2, &
    gzwpsu2sd1_2_1_c, gzwpsu1sd2_2_1_c, gzwpsu2sd2_2_1_c, gzwpsu1sd1_2_1_c, &
    gzwpsu2sd1_2_1, gzwpsu1sd2_2_1, gzwpsu2sd2_2_1, gzwpsu1sd1_2_1, &
    gpwpsu2sd1_2_1_c, gpwpsu1sd2_2_1_c, gpwpsu2sd2_2_1_c, gpwpsu1sd1_2_1_c, &
    gpwpsu2sd1_2_1, gpwpsu1sd2_2_1, gpwpsu2sd2_2_1, gpwpsu1sd1_2_1, &
    gzwpsu2sd1_1_3_c, gzwpsu1sd2_1_3_c, gzwpsu2sd2_1_3_c, gzwpsu1sd1_1_3_c, &
    gzwpsu2sd1_1_3, gzwpsu1sd2_1_3, gzwpsu2sd2_1_3, gzwpsu1sd1_1_3, &
    gpwpsu2sd1_1_3_c, gpwpsu1sd2_1_3_c, gpwpsu2sd2_1_3_c, gpwpsu1sd1_1_3_c, &
    gpwpsu2sd1_1_3, gpwpsu1sd2_1_3, gpwpsu2sd2_1_3, gpwpsu1sd1_1_3, &
    gzwpsu2sd1_1_2_c, gzwpsu1sd2_1_2_c, gzwpsu2sd2_1_2_c, gzwpsu1sd1_1_2_c, &
    gzwpsu2sd1_1_2, gzwpsu1sd2_1_2, gzwpsu2sd2_1_2, gzwpsu1sd1_1_2, &
    gpwpsu2sd1_1_2_c, gpwpsu1sd2_1_2_c, gpwpsu2sd2_1_2_c, gpwpsu1sd1_1_2_c, &
    gpwpsu2sd1_1_2, gpwpsu1sd2_1_2, gpwpsu2sd2_1_2
  complex(kind=default), public, save :: gpwpsu1sd1_1_2, gzwpsu2sd1_1_1_c, &
    gzwpsu1sd2_1_1_c, gzwpsu2sd2_1_1_c, gzwpsu1sd1_1_1_c, gzwpsu2sd1_1_1, &
    gzwpsu1sd2_1_1, gzwpsu2sd2_1_1, gzwpsu1sd1_1_1, gpwpsu2sd1_1_1_c, &
    gpwpsu1sd2_1_1_c, gpwpsu2sd2_1_1_c, gpwpsu1sd1_1_1_c, gpwpsu2sd1_1_1, &
    gpwpsu1sd2_1_1, gpwpsu2sd2_1_1, gpwpsu1sd1_1_1, gwzsl2sn_3_c, gwzsl1sn_3_c, &
    gwzsl2sn_3, gwzsl1sn_3, gpwsl2sn_3_c, gpwsl1sn_3_c, gpwsl2sn_3, gpwsl1sn_3, &
    gwwsd2sd1_3, gwwsd1sd2_3, gwwsd2sd2_3, gwwsd1sd1_3, gwwsu2su1_3, &
    gwwsu1su2_3, gwwsu2su2_3, gwwsu1su1_3, gwwsn1sn1_3, gwwsl2sl1_3, &
    gwwsl1sl2_3, gwwsl2sl2_3, gwwsl1sl1_3, gzpsd2sd1_3, gzpsd1sd2_3, &
    gzpsd2sd2_3, gzpsd1sd1_3, gzpsu2su1_3, gzpsu1su2_3, gzpsu2su2_3, &
    gzpsu1su1_3, gzpsl2sl1_3, gzpsl1sl2_3, gzpsl2sl2_3, gzpsl1sl1_3, &
    gzzsd2sd1_3, gzzsd1sd2_3, gzzsd2sd2_3, gzzsd1sd1_3, gzzsu2su1_3, &
    gzzsu1su2_3, gzzsu2su2_3, gzzsu1su1_3, gzzsn1sn1_3, gzzsl2sl1_3, &
    gzzsl1sl2_3, gzzsl2sl2_3, gzzsl1sl1_3, gwzsl2sn_2_c, gwzsl1sn_2_c, &
    gwzsl2sn_2, gwzsl1sn_2, gpwsl2sn_2_c, gpwsl1sn_2_c
  complex(kind=default), public, save :: gpwsl2sn_2, gpwsl1sn_2, &
    gwwsd2sd1_2, gwwsd1sd2_2, gwwsd2sd2_2, gwwsd1sd1_2, gwwsu2su1_2, &
    gwwsu1su2_2, gwwsu2su2_2, gwwsu1su1_2, gwwsn1sn1_2, gwwsl2sl1_2, &
    gwwsl1sl2_2, gwwsl2sl2_2, gwwsl1sl1_2, gzpsd2sd1_2, gzpsd1sd2_2, &
    gzpsd2sd2_2, gzpsd1sd1_2, gzpsu2su1_2, gzpsu1su2_2, gzpsu2su2_2, &
    gzpsu1su1_2, gzpsl2sl1_2, gzpsl1sl2_2, gzpsl2sl2_2, gzpsl1sl1_2, &
    gzzsd2sd1_2, gzzsd1sd2_2, gzzsd2sd2_2, gzzsd1sd1_2, gzzsu2su1_2, &
    gzzsu1su2_2, gzzsu2su2_2, gzzsu1su1_2, gzzsn1sn1_2, gzzsl2sl1_2, &
    gzzsl1sl2_2, gzzsl2sl2_2, gzzsl1sl1_2, gwzsl2sn_1_c, gwzsl1sn_1_c, &
    gwzsl2sn_1, gwzsl1sn_1, gpwsl2sn_1_c, gpwsl1sn_1_c, gpwsl2sn_1, gpwsl1sn_1, &
    gwwsd2sd1_1, gwwsd1sd2_1, gwwsd2sd2_1, gwwsd1sd1_1, gwwsu2su1_1, &
    gwwsu1su2_1, gwwsu2su2_1, gwwsu1su1_1, gwwsn1sn1_1, gwwsl2sl1_1, &
    gwwsl1sl2_1, gwwsl2sl2_1, gwwsl1sl1_1, gzpsd2sd1_1, gzpsd1sd2_1, &
    gzpsd2sd2_1, gzpsd1sd1_1, gzpsu2su1_1, gzpsu1su2_1, gzpsu2su2_1, &
    gzpsu1su1_1
  complex(kind=default), public, save :: gzpsl2sl1_1, gzpsl1sl2_1, &
    gzpsl2sl2_1, gzpsl1sl1_1, gzzsd2sd1_1, gzzsd1sd2_1, gzzsd2sd2_1, &
    gzzsd1sd1_1, gzzsu2su1_1, gzzsu1su2_1, gzzsu2su2_1, gzzsu1su1_1, &
    gzzsn1sn1_1, gzzsl2sl1_1, gzzsl1sl2_1, gzzsl2sl2_1, gzzsl1sl1_1, gppsdsd, &
    gppsusu, gppslsl, gsl2_3snw_c, gsl1_3snw_c, gsl2_3snw, gsl1_3snw, &
    gsd2zsd1_3, gsd1zsd2_3, gsd2zsd2_3, gsd1zsd1_3, gsu2zsu1_3, gsu1zsu2_3, &
    gsu2zsu2_3, gsu1zsu1_3, gsn1zsn1_3, gsl2zsl1_3, gsl1zsl2_3, gsl2zsl2_3, &
    gsl1zsl1_3, gsl2_2snw_c, gsl1_2snw_c, gsl2_2snw, gsl1_2snw, gsd2zsd1_2, &
    gsd1zsd2_2, gsd2zsd2_2, gsd1zsd1_2, gsu2zsu1_2, gsu1zsu2_2, gsu2zsu2_2, &
    gsu1zsu1_2, gsn1zsn1_2, gsl2zsl1_2, gsl1zsl2_2, gsl2zsl2_2, gsl1zsl1_2, &
    gsl2_1snw_c, gsl1_1snw_c, gsl2_1snw, gsl1_1snw, gsd2zsd1_1, gsd1zsd2_1, &
    gsd2zsd2_1, gsd1zsd1_1, gsu2zsu1_1, gsu1zsu2_1, gsu2zsu2_1, gsu1zsu1_1, &
    gsn1zsn1_1, gsl2zsl1_1, gsl1zsl2_1
  complex(kind=default), public, save :: gsl2zsl2_1, gsl1zsl1_1,  &
    gs2ws1_3_3_c, gs1ws2_3_3_c, gs2ws2_3_3_c, gs1ws1_3_3_c, gs2ws1_3_3, &
    gs1ws2_3_3, gs2ws2_3_3, gs1ws1_3_3, gs2ws1_3_2_c, gs1ws2_3_2_c, &
    gs2ws2_3_2_c, gs1ws1_3_2_c, gs2ws1_3_2, gs1ws2_3_2, gs2ws2_3_2, &
    gs1ws1_3_2, gs2ws1_3_1_c, gs1ws2_3_1_c, gs2ws2_3_1_c, gs1ws1_3_1_c, &
    gs2ws1_3_1, gs1ws2_3_1, gs2ws2_3_1, gs1ws1_3_1, gs2ws1_2_3_c, &
    gs1ws2_2_3_c, gs2ws2_2_3_c, gs1ws1_2_3_c, gs2ws1_2_3, gs1ws2_2_3, &
    gs2ws2_2_3, gs1ws1_2_3, gs2ws1_2_2_c, gs1ws2_2_2_c, gs2ws2_2_2_c, &
    gs1ws1_2_2_c, gs2ws1_2_2, gs1ws2_2_2, gs2ws2_2_2, gs1ws1_2_2, &
    gs2ws1_2_1_c, gs1ws2_2_1_c, gs2ws2_2_1_c, gs1ws1_2_1_c, gs2ws1_2_1, &
    gs1ws2_2_1, gs2ws2_2_1, gs1ws1_2_1, gs2ws1_1_3_c, gs1ws2_1_3_c, &
    gs2ws2_1_3_c, gs1ws1_1_3_c, gs2ws1_1_3, gs1ws2_1_3, gs2ws2_1_3, &
    gs1ws1_1_3, gs2ws1_1_2_c, gs1ws2_1_2_c, gs2ws2_1_2_c, gs1ws1_1_2_c, &
    gs2ws1_1_2, gs1ws2_1_2, gs2ws2_1_2, gs1ws1_1_2, gs2ws1_1_1_c, &
    gs1ws2_1_1_c
  complex(kind=default), public, save :: gs2ws2_1_1_c, gs1ws1_1_1_c, &
    gs2ws1_1_1, gs1ws2_1_1, gs2ws2_1_1, gs1ws1_1_1, g_yuk15_3, g_yuk14_3, &
    g_yuk13_3, g_yuk12_3, g_yuk11_3, g_yuk10_3, g_yuk9_3, g_yuk8_3, &
    g_yuk7_3, g_yuk6_3, g_yuk15_2, g_yuk14_2, g_yuk13_2, g_yuk12_2, &
    g_yuk11_2, g_yuk10_2, g_yuk9_2, g_yuk8_2, g_yuk7_2, g_yuk6_2, g_yuk15_1, &
    g_yuk14_1, g_yuk13_1, g_yuk12_1, g_yuk11_1, g_yuk10_1, g_yuk9_1, &
    g_yuk8_1, g_yuk7_1, g_yuk6_1, ghhww, gh2h2ww, gh1h1ww, gaaww, ghh2wp, &
    ghawp, ghawz, gh2az, gh1az, ghaw, ghh1wp, ghh2wz, ghh1wz, ghphmpz, & 
    ghphmpp, ghphmzz, gh2h2zz, gh1h1zz, gaazz, ghhp, ghhz, gh2zz, gh1zz, &
    ghh2w, ghh1w, gh2ww, gh1ww, gh4_11, gh4_10, gh4_9, gh4_8, gh4_7, gh4_6, &
    gh4_5, gh4_4
  complex(kind=default), public, save :: gh4_3, gh4_2, gh4_1, gh3_8, &
    gh3_7, gh3_6, gh3_5, gh3_4, gh3_3, gh3_2, gh3_1, mu, ad_3, au_3, al_3, &
    ad_2, au_2, al_2, ad_1, au_1, al_1, mv_22, mv_21, mv_12, mv_11, mu_22, &
    mu_21, mu_12, mu_11, mn_44, mn_43, mn_42, mn_41, mn_34, mn_33, mn_32, &
    mn_31, mn_24, mn_23, mn_22, mn_21, mn_14, mn_13, mn_12, mn_11
  !!! complex(kind=default), public, save :: sinthsu3, &
  !!!   sinthsu2, sinthsu1, sinthsd3, sinthsd2, sinthsd1, sinthsl3, sinthsl2, &
  !!!   sinthsl1, costhsu3, costhsu2, costhsu1, costhsd3, costhsd2, costhsd1, &
  !!!   costhsl3, costhsl2, costhsl1 
  complex(kind=default), public, save :: eta1, eta2, eta3, eta4
  complex(kind=default), public, save :: eidelta, cosckm23, cosckm13, &
    cosckm12, vckm_33, vckm_32, vckm_31, vckm_23, vckm_22, vckm_21, vckm_13, &
    vckm_12, vckm_11, gpzww, gppww, gzzww, gw4, igwww, igzww, iqw, igs, &
    gssq
  complex(kind=default), public, save :: gccq_3_3_c, gccq_3_3, & 
    gccq_3_2_c, gccq_3_2, gccq_3_1_c, gccq_3_1, gccq_2_3_c, gccq_2_3, & 
    gccq_2_2_c, gccq_2_2, gccq_2_1_c, gccq_2_1, gccq_1_3_c, gccq_1_3, &
    gccq_1_2_c, gccq_1_2, gccq_1_1_c, gccq_1_1
  complex(kind=default), dimension(2), public, save :: g_yuk_gsd2_3_c, &
    g_yuk_gsd2_3, g_yuk_gsu2_3_c, g_yuk_gsu2_3, g_yuk_gsd1_3_c, &
    g_yuk_gsd1_3, g_yuk_gsu1_3_c, g_yuk_gsu1_3, g_yuk_n4_sd2_3_c, &
    g_yuk_n4_sd2_3, g_yuk_n4_su2_3_c, g_yuk_n4_su2_3, g_yuk_n4_sl2_3_c, &
    g_yuk_n4_sl2_3, g_yuk_n3_sd2_3_c, g_yuk_n3_sd2_3, g_yuk_n3_su2_3_c, &
    g_yuk_n3_su2_3, g_yuk_n3_sl2_3_c, g_yuk_n3_sl2_3, g_yuk_n2_sd2_3_c, &
    g_yuk_n2_sd2_3, g_yuk_n2_su2_3_c, g_yuk_n2_su2_3, g_yuk_n2_sl2_3_c, &
    g_yuk_n2_sl2_3, g_yuk_n1_sd2_3_c, g_yuk_n1_sd2_3, g_yuk_n1_su2_3_c, &
    g_yuk_n1_su2_3, g_yuk_n1_sl2_3_c, g_yuk_n1_sl2_3, g_yuk_n4_sd1_3_c, &
    g_yuk_n4_sd1_3, g_yuk_n4_su1_3_c, g_yuk_n4_su1_3, g_yuk_n4_sl1_3_c, &
    g_yuk_n4_sl1_3, g_yuk_n3_sd1_3_c, g_yuk_n3_sd1_3, g_yuk_n3_su1_3_c, &
    g_yuk_n3_su1_3, g_yuk_n3_sl1_3_c, g_yuk_n3_sl1_3, g_yuk_n2_sd1_3_c, &
    g_yuk_n2_sd1_3, g_yuk_n2_su1_3_c, g_yuk_n2_su1_3, g_yuk_n2_sl1_3_c, &
    g_yuk_n2_sl1_3, g_yuk_n1_sd1_3_c, g_yuk_n1_sd1_3, g_yuk_n1_su1_3_c, &
    g_yuk_n1_su1_3, g_yuk_n1_sl1_3_c, g_yuk_n1_sl1_3
  complex(kind=default), dimension(2), public, save :: g_yuk_ch2_su2_3_3_c, &
    g_yuk_ch2_su2_3_3, g_yuk_ch2_sd2_3_3_c, g_yuk_ch2_sd2_3_3, &
    g_yuk_ch2_su1_3_3_c, g_yuk_ch2_su1_3_3, g_yuk_ch2_sd1_3_3_c, &
    g_yuk_ch2_sd1_3_3, g_yuk_ch1_su2_3_3_c, g_yuk_ch1_su2_3_3, &
    g_yuk_ch1_sd2_3_3_c, g_yuk_ch1_sd2_3_3, g_yuk_ch1_su1_3_3_c, &
    g_yuk_ch1_su1_3_3, g_yuk_ch1_sd1_3_3_c, g_yuk_ch1_sd1_3_3, &
    g_yuk_ch2_su2_3_2_c, g_yuk_ch2_su2_3_2, g_yuk_ch2_sd2_3_2_c, &
    g_yuk_ch2_sd2_3_2, g_yuk_ch2_su1_3_2_c, g_yuk_ch2_su1_3_2, &
    g_yuk_ch2_sd1_3_2_c, g_yuk_ch2_sd1_3_2, g_yuk_ch1_su2_3_2_c, &
    g_yuk_ch1_su2_3_2, g_yuk_ch1_sd2_3_2_c, g_yuk_ch1_sd2_3_2, &
    g_yuk_ch1_su1_3_2_c, g_yuk_ch1_su1_3_2, g_yuk_ch1_sd1_3_2_c, &
    g_yuk_ch1_sd1_3_2, g_yuk_ch2_su2_3_1_c, g_yuk_ch2_su2_3_1, &
    g_yuk_ch2_sd2_3_1_c, g_yuk_ch2_sd2_3_1, g_yuk_ch2_su1_3_1_c, &
    g_yuk_ch2_su1_3_1, g_yuk_ch2_sd1_3_1_c, g_yuk_ch2_sd1_3_1, &
    g_yuk_ch1_su2_3_1_c, g_yuk_ch1_su2_3_1, g_yuk_ch1_sd2_3_1_c, &
    g_yuk_ch1_sd2_3_1, g_yuk_ch1_su1_3_1_c, g_yuk_ch1_su1_3_1, &
    g_yuk_ch1_sd1_3_1_c, g_yuk_ch1_sd1_3_1, g_yuk_ch2_su2_2_3_c, &
    g_yuk_ch2_su2_2_3, g_yuk_ch2_sd2_2_3_c, g_yuk_ch2_sd2_2_3, &
    g_yuk_ch2_su1_2_3_c, g_yuk_ch2_su1_2_3, g_yuk_ch2_sd1_2_3_c, &
    g_yuk_ch2_sd1_2_3, g_yuk_ch1_su2_2_3_c, g_yuk_ch1_su2_2_3, &
    g_yuk_ch1_sd2_2_3_c, g_yuk_ch1_sd2_2_3, g_yuk_ch1_su1_2_3_c, &
    g_yuk_ch1_su1_2_3, g_yuk_ch1_sd1_2_3_c, g_yuk_ch1_sd1_2_3, &
    g_yuk_ch2_su2_1_3_c, g_yuk_ch2_su2_1_3, g_yuk_ch2_sd2_1_3_c, &
    g_yuk_ch2_sd2_1_3, g_yuk_ch2_su1_1_3_c, g_yuk_ch2_su1_1_3, &
    g_yuk_ch2_sd1_1_3_c, g_yuk_ch2_sd1_1_3, g_yuk_ch1_su2_1_3_c, &
    g_yuk_ch1_su2_1_3, g_yuk_ch1_sd2_1_3_c, g_yuk_ch1_sd2_1_3, &
    g_yuk_ch1_su1_1_3_c, g_yuk_ch1_su1_1_3, g_yuk_ch1_sd1_1_3_c, &
    g_yuk_ch1_sd1_1_3, g_yuk_ch2_sn1_3_c, g_yuk_ch2_sn1_3, &
    g_yuk_ch1_sn1_3_c, g_yuk_ch1_sn1_3
  complex(kind=default), dimension(2), public, save :: gcac_2_1, &
    gch2c_2_1, gch1c_2_1, gcac_1_2, gch2c_1_2, gch1c_1_2, gcicia_3_4, &
    gcicih2_3_4, gcicih1_3_4, gcicia_2_4, gcicih2_2_4, gcicih1_2_4, &
    gcicia_2_3, gcicih2_2_3, gcicih1_2_3, gcicia_1_4, gcicih2_1_4, &
    gcicih1_1_4, gcicia_1_3, gcicih2_1_3, gcicih1_1_3, gcicia_1_2, &
    gcicih2_1_2, gcicih1_1_2, g_chn_4_2, gcwn_2_4, g_chn_3_2, gcwn_2_3, &
    g_chn_2_2, gcwn_2_2, g_chn_1_2, gcwn_2_1, g_chn_4_1, gcwn_1_4, g_chn_3_1, &
    gcwn_1_3, g_chn_2_1, gcwn_1_2, g_chn_1_1, gcwn_1_1, g_nhc_4_2, gnwc_4_2, &
    g_nhc_4_1, gnwc_4_1, g_nhc_3_2, gnwc_3_2, g_nhc_3_1, gnwc_3_1, g_nhc_2_2, &
    gnwc_2_2, g_nhc_2_1, gnwc_2_1, g_nhc_1_2, gnwc_1_2, g_nhc_1_1, gnwc_1_1, &
    gczc_2_2, gczc_2_1, gczc_1_2, gczc_1_1, gnzn_3_4, gnzn_2_4, gnzn_2_3, &
    gnzn_1_4, gnzn_1_3, gnzn_1_2, g_yuk2_3_3, g_yuk2_3_2, g_yuk2_3_1, &
    g_yuk2_2_3, g_yuk2_1_3, g_yuk1_3_3, g_yuk1_3_2, g_yuk1_3_1, g_yuk1_2_3, &
    g_yuk1_1_3
  !!! Specific gravitino parameters
  real(kind=default), public, save :: m_planck
  real(kind=default), public, save :: gzll, gzlr, gzul, gzur, &
        gzdl, gzdr
  complex(kind=default), public, save :: &
      ggrav, ggrneu1, ggrneu2, ggrneu3, ggrneu4, ggravn
  complex(kind=default), dimension(2), public, save :: &
      ggrch1, ggrch2, ggrch1c, ggrch2c, ggrzneu1, ggrzneu2, ggrzneu3, &
      ggrzneu4, ggraneu1, ggraneu2, ggraneu3, ggraneu4, ggr4neu1, &  
      ggr4neu2, ggr4neu3, ggr4neu4, ggr4ach1, ggr4ach2, ggr4ach1c, &
      ggr4ach2c, ggr4zch1, ggr4zch2, ggr4zch1c, ggr4zch2c
  complex(kind=default), dimension(2), public, save :: &
      ggravu11, ggravu12, ggravu11c, ggravu12c, ggravu21, ggravu22, ggravu21c, &
      ggravu22c, ggravu31, ggravu32, ggravu31c, ggravu32c, &
      ggravd11, ggravd12, ggravd11c, ggravd12c, ggravd21, ggravd22, ggravd21c, &
      ggravd22c, ggravd31, ggravd32, ggravd31c, ggravd32c, &
      ggravl11, ggravl12, ggravl11c, ggravl12c, ggravl21, ggravl22, ggravl21c, &
      ggravl22c, ggravl31, ggravl32, ggravl31c, ggravl32c, &
      ggrhch1, ggrhch2, ggrhch1c, ggrhch2c, ggrh1neu1, ggrh1neu2, &      
      ggrh1neu3, ggrh1neu4, ggrh2neu1, ggrh2neu2, ggrh2neu3, &
      ggrh2neu4, ggrh3neu1, ggrh3neu2, ggrh3neu3, ggrh3neu4
  complex(kind=default), dimension(2), public, save :: &
       ggr4asl11, ggr4asl12, ggr4asl21, ggr4asl22, &
       ggr4asl31, ggr4asl32, ggr4asl11c, ggr4asl12c, ggr4asl21c, &
       ggr4asl22c, ggr4asl31c, ggr4asl32c, ggr4asu11, ggr4asu12, &       
       ggr4asu21, ggr4asu22, ggr4asu31, ggr4asu32, ggr4asu11c, ggr4asu12c, &
       ggr4asu21c, ggr4asu22c, ggr4asu31c, ggr4asu32c, ggr4asd11, &
       ggr4asd12, ggr4asd21, ggr4asd22, ggr4asd31, ggr4asd32, ggr4asd11c, &
       ggr4asd12c, ggr4asd21c, ggr4asd22c, ggr4asd31c, ggr4asd32c, &
       ggr4zsl11, ggr4zsl12, ggr4zsl21, ggr4zsl22, ggr4zsl31, ggr4zsl32, &
       ggr4zsl11c, ggr4zsl12c, ggr4zsl21c, ggr4zsl22c, ggr4zsl31c, &
       ggr4zsl32c, ggr4zsu11, ggr4zsu12, ggr4zsu21, ggr4zsu22, ggr4zsu31, & 
       ggr4zsu32, ggr4zsu11c, ggr4zsu12c, ggr4zsu21c, ggr4zsu22c, ggr4zsu31c, &
       ggr4zsu32c, ggr4zsd11, ggr4zsd12, ggr4zsd21, ggr4zsd22, ggr4zsd31, &
       ggr4zsd32, ggr4zsd11c, ggr4zsd12c, ggr4zsd21c, ggr4zsd22c, ggr4zsd31c, &       
       ggr4zsd32c
  complex(kind=default), public, save :: &
       ggr4zsn, ggr4zsnc, ggr4wsn, ggr4wsnc, ggr4wsl11, ggr4wsl12, ggr4wsl21, &
       ggr4wsl22, ggr4wsl31, ggr4wsl32, ggr4wsl11c, ggr4wsl12c, ggr4wsl21c, &
       ggr4wsl22c, ggr4wsl31c, ggr4wsl32c, ggr4wsu11, ggr4wsu12, ggr4wsu21, &
       ggr4wsu22, ggr4wsu31, ggr4wsu32, ggr4wsu11c, ggr4wsu12c, ggr4wsu21c, &
       ggr4wsu22c, ggr4wsu31c, ggr4wsu32c, ggr4wsd11, ggr4wsd12, ggr4wsd21, &
       ggr4wsd22, ggr4wsd31, ggr4wsd32, ggr4wsd11c, ggr4wsd12c, ggr4wsd21c, &
       ggr4wsd22c, ggr4wsd31c, ggr4wsd32c
  complex(kind=default), dimension(2), public, save :: &
       ggr4glsu11, ggr4glsu12, ggr4glsu21, ggr4glsu22, ggr4glsu31, ggr4glsu32, &
       ggr4glsu11c, ggr4glsu12c, ggr4glsu21c, ggr4glsu22c, ggr4glsu31c, &
       ggr4glsu32c, ggr4glsd11, ggr4glsd12, ggr4glsd21, ggr4glsd22, ggr4glsd31, &
       ggr4glsd32, ggr4glsd11c, ggr4glsd12c, ggr4glsd21c, ggr4glsd22c, ggr4glsd31c, &
       ggr4glsd32c, ggr4zh1_1, ggr4zh1_2, ggr4zh1_3, ggr4zh1_4, ggr4zh2_1,  &
       ggr4zh2_2, ggr4zh2_3, ggr4zh2_4, ggr4zh3_1, ggr4zh3_2, ggr4zh3_3, ggr4zh3_4, &  
       ggr4wh_1, ggr4wh_2, ggr4wh_1c, ggr4wh_2c, ggr4ha1, ggr4ha2, ggr4ha1c, &
       ggr4ha2c, ggr4hz1, ggr4hz2, ggr4hz1c, ggr4hz2c  

contains

  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(138), intent(in) :: par_array
    integer, intent(in) :: scheme
    type :: parameter_set 
       real(default) :: gf
       real(default) :: mZ
       real(default) :: wZ
       real(default) :: mW
       real(default) :: wW
       real(default) :: me
       real(default) :: mmu
       real(default) :: mtau
       real(default) :: ms
       real(default) :: mc
       real(default) :: mb
       real(default) :: mtop
       real(default) :: wtop
       real(default) :: alphas
       real(default) :: mtype
       real(default) :: m_zero
       real(default) :: m_half
       real(default) :: A0
       real(default) :: tanb
       real(default) :: sgn_mu
       real(default) :: lambda
       real(default) :: m_mes
       real(default) :: n5
       real(default) :: c_grav
       real(default) :: m_grav
       real(default) :: ae_33
       real(default) :: au_33
       real(default) :: ad_33
       real(default) :: mh
       real(default) :: wh
       real(default) :: mhh
       real(default) :: mha
       real(default) :: mhpm
       real(default) :: whh
       real(default) :: whpm
       real(default) :: wha
       real(default) :: al_h
       real(default) :: mu_h
       real(default) :: tanb_h
       real(default) :: msu1
       real(default) :: msd1
       real(default) :: msc1
       real(default) :: mss1
       real(default) :: mstop1
       real(default) :: msb1
       real(default) :: msu2
       real(default) :: msd2
       real(default) :: msc2
       real(default) :: mss2
       real(default) :: mstop2
       real(default) :: msb2
       real(default) :: mse1
       real(default) :: msne
       real(default) :: msmu1
       real(default) :: msnmu
       real(default) :: mstau1
       real(default) :: msntau
       real(default) :: mse2
       real(default) :: msmu2
       real(default) :: mstau2
       real(default) :: mgg
       real(default) :: mgr
       real(default) :: mch1
       real(default) :: mch2
       real(default) :: mneu1
       real(default) :: mneu2
       real(default) :: mneu3
       real(default) :: mneu4
       real(default) :: wsu1
       real(default) :: wsd1
       real(default) :: wsc1
       real(default) :: wss1
       real(default) :: wstop1
       real(default) :: wsb1
       real(default) :: wsu2
       real(default) :: wsd2
       real(default) :: wsc2
       real(default) :: wss2
       real(default) :: wstop2
       real(default) :: wsb2
       real(default) :: wse1
       real(default) :: wsne
       real(default) :: wsmu1
       real(default) :: wsnmu
       real(default) :: wstau1
       real(default) :: wsntau
       real(default) :: wse2
       real(default) :: wsmu2
       real(default) :: wstau2
       real(default) :: wgg
       real(default) :: wch1
       real(default) :: wch2
       real(default) :: wneu1
       real(default) :: wneu2
       real(default) :: wneu3
       real(default) :: wneu4
       real(default) :: wgr
       real(default) :: mt_11
       real(default) :: mt_12
       real(default) :: mt_21
       real(default) :: mt_22
       real(default) :: mb_11
       real(default) :: mb_12
       real(default) :: mb_21
       real(default) :: mb_22
       real(default) :: ml_11
       real(default) :: ml_12
       real(default) :: ml_21
       real(default) :: ml_22
       real(default) :: mn_11
       real(default) :: mn_12
       real(default) :: mn_13
       real(default) :: mn_14
       real(default) :: mn_21
       real(default) :: mn_22
       real(default) :: mn_23
       real(default) :: mn_24
       real(default) :: mn_31
       real(default) :: mn_32
       real(default) :: mn_33
       real(default) :: mn_34
       real(default) :: mn_41
       real(default) :: mn_42
       real(default) :: mn_43
       real(default) :: mn_44
       real(default) :: mu_11
       real(default) :: mu_12
       real(default) :: mu_21
       real(default) :: mu_22
       real(default) :: mv_11
       real(default) :: mv_12
       real(default) :: mv_21
       real(default) :: mv_22
       real(default) :: m_pl
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    real(kind=default) :: sinthw, costhw, qelep, qeup, qedwn, v
    par%gf     = par_array(1)
    par%mZ     = par_array(2)
    par%wZ     = par_array(3)
    par%mW     = par_array(4)
    par%wW     = par_array(5)
    par%me     = par_array(6)
    par%mmu    = par_array(7)
    par%mtau   = par_array(8)
    par%ms     = par_array(9)
    par%mc     = par_array(10)
    par%mb     = par_array(11)
    par%mtop   = par_array(12)
    par%wtop   = par_array(13)
    par%alphas = par_array(14)
    par%mtype  = par_array(15)
    par%m_zero = par_array(16)
    par%m_half = par_array(17)
    par%A0     = par_array(18)
    par%tanb   = par_array(19)
    par%sgn_mu = par_array(20)
    par%lambda = par_array(21)
    par%m_mes  = par_array(22)
    par%n5     = par_array(23)
    par%c_grav = par_array(24)
    par%m_grav = par_array(25)
    par%ae_33  = par_array(26)
    par%au_33  = par_array(27)
    par%ad_33  = par_array(28)
    par%mh     = par_array(29)
    par%wh     = par_array(30)
    par%mhh    = par_array(31)
    par%mha    = par_array(32)
    par%mhpm   = par_array(33)
    par%whh    = par_array(34)
    par%whpm   = par_array(35)
    par%wha    = par_array(36)
    par%al_h   = par_array(37)
    par%mu_h   = par_array(38)
    par%tanb_h = par_array(39)
    par%msu1   = par_array(40)
    par%msd1   = par_array(41)
    par%msc1   = par_array(42)
    par%mss1   = par_array(43)
    par%mstop1 = par_array(44)
    par%msb1   = par_array(45)
    par%msu2   = par_array(46)
    par%msd2   = par_array(47)
    par%msc2   = par_array(48)
    par%mss2   = par_array(49)
    par%mstop2 = par_array(50)
    par%msb2   = par_array(51)
    par%mse1   = par_array(52)
    par%msne   = par_array(53)
    par%msmu1  = par_array(54)
    par%msnmu  = par_array(55)
    par%mstau1 = par_array(56)
    par%msntau = par_array(57)
    par%mse2   = par_array(58)
    par%msmu2  = par_array(59)
    par%mstau2 = par_array(60)
    par%mgg    = par_array(61)
    par%mgr    = par_array(62)
    par%mch1   = par_array(63)
    par%mch2   = par_array(64)
    par%mneu1  = par_array(65)
    par%mneu2  = par_array(66)
    par%mneu3  = par_array(67)
    par%mneu4  = par_array(68)
    par%wsu1   = par_array(69)
    par%wsd1   = par_array(70)
    par%wsc1   = par_array(71)
    par%wss1   = par_array(72)
    par%wstop1 = par_array(73)
    par%wsb1   = par_array(74)
    par%wsu2   = par_array(75)
    par%wsd2   = par_array(76)
    par%wsc2   = par_array(77)
    par%wss2   = par_array(78)
    par%wstop2 = par_array(79)
    par%wsb2   = par_array(80)
    par%wse1   = par_array(81)
    par%wsne   = par_array(82)
    par%wsmu1  = par_array(83)
    par%wsnmu  = par_array(84)
    par%wstau1 = par_array(85)
    par%wsntau = par_array(86)
    par%wse2   = par_array(87)
    par%wsmu2  = par_array(88)
    par%wstau2 = par_array(89)
    par%wgg    = par_array(90)
    par%wch1   = par_array(91)
    par%wch2   = par_array(92)
    par%wneu1  = par_array(93)
    par%wneu2  = par_array(94)
    par%wneu3  = par_array(95)
    par%wneu4  = par_array(96)
    par%wgr    = par_array(97)
    par%mt_11  = par_array(98)
    par%mt_12  = par_array(99)
    par%mt_21  = par_array(100)
    par%mt_22  = par_array(101)
    par%mb_11  = par_array(102)
    par%mb_12  = par_array(103)
    par%mb_21  = par_array(104)
    par%mb_22  = par_array(105)
    par%ml_11  = par_array(106)
    par%ml_12  = par_array(107)
    par%ml_21  = par_array(108)
    par%ml_22  = par_array(109)
    par%mn_11  = par_array(110)
    par%mn_12  = par_array(111)
    par%mn_13  = par_array(112)
    par%mn_14  = par_array(113)
    par%mn_21  = par_array(114)
    par%mn_22  = par_array(115)
    par%mn_23  = par_array(116)
    par%mn_24  = par_array(117)
    par%mn_31  = par_array(118)
    par%mn_32  = par_array(119)
    par%mn_33  = par_array(120)
    par%mn_34  = par_array(121)
    par%mn_41  = par_array(122)
    par%mn_42  = par_array(123)
    par%mn_43  = par_array(124)
    par%mn_44  = par_array(125)
    par%mu_11  = par_array(126)
    par%mu_12  = par_array(127)
    par%mu_21  = par_array(128)
    par%mu_22  = par_array(129)
    par%mv_11  = par_array(130)
    par%mv_12  = par_array(131)
    par%mv_21  = par_array(132)
    par%mv_22  = par_array(133)
    par%m_pl   = par_array(134)
    par%v      = par_array(135)
    par%cw     = par_array(136)
    par%sw     = par_array(137)
    par%ee     = par_array(138)
    mass(1:70) = 0
    width(1:70) = 0
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
    mass(25) = par%mh
    width(25) = par%wh
    mass(26) =  xi0 * mass(23)
    width(26) =  0
    mass(27) =  xipm * mass(24)
    width(27) =  0
    mass(35) = par%mHH
    width(35) = par%wHH
    mass(36) = par%mHA
    width(36) = par%wHA
    mass(37) = par%mHpm
    width(37) = par%wHpm
    mass(39) = par%mgr
    width(39) = par%wgr
    mass(41) = par%msd1
    width(41) = par%wsd1
    mass(42) = par%msu1
    width(42) = par%wsu1
    mass(43) = par%mss1
    width(43) = par%wss1
    mass(44) = par%msc1
    width(44) = par%wsc1
    mass(45) = par%msb1
    width(45) = par%wsb1
    mass(46) = par%mstop1
    width(46) = par%wstop1
    mass(47) = par%msd2
    width(47) = par%wsd2
    mass(48) = par%msu2
    width(48) = par%wsu2
    mass(49) = par%mss2
    width(49) = par%wss2
    mass(50) = par%msc2
    width(50) = par%wsc2
    mass(51) = par%msb2
    width(51) = par%wsb2
    mass(52) = par%mstop2
    width(52) = par%wstop2
    mass(53) = par%mse1
    width(53) = par%wse1
    mass(54) = par%msne
    width(54) = par%wsne
    mass(55) = par%msmu1
    width(55) = par%wsmu1
    mass(56) = par%msnmu
    width(56) = par%wsnmu
    mass(57) = par%mstau1
    width(57) = par%wstau1
    mass(58) = par%msntau
    width(58) = par%wsntau
    mass(59) = par%mse2
    width(59) = par%wse2
    mass(61) = par%msmu2
    width(61) = par%wsmu2
    mass(63) = par%mstau2
    width(63) = par%wstau2
    mass(64) = par%mgg
    width(64) = par%wgg
    mass(65) = abs(par%mneu1)
    width(65) = par%wneu1
    mass(66) = abs(par%mneu2)
    width(66) = par%wneu2
    mass(67) = abs(par%mneu3)
    width(67) = par%wneu3
    mass(68) = abs(par%mneu4)
    width(68) = par%wneu4
    mass(69) = abs(par%mch1)
    width(69) = par%wch1
    mass(70) = abs(par%mch2)
    width(70) = par%wch2
    sigch1   = sign (1._default, par%mch1)
    sigch2   = sign (1._default, par%mch2)
    sign1    = sign (1, int(par%mneu1))
    sign2    = sign (1, int(par%mneu2))
    sign3    = sign (1, int(par%mneu3))
    sign4    = sign (1, int(par%mneu4))
    vckm_11  = 1
    vckm_12  = 0
    vckm_13  = 0
    vckm_21  = 0
    vckm_22  = 1
    vckm_23  = 0
    vckm_31  = 0
    vckm_32  = 0
    vckm_33  = 1
    v = 2 * par%mW * par%sw / par%ee
    e = par%ee
    !!! This should not be in the color flow basis !!!
    as = par%alphas
    tanb = par%tanb_h
    tana = tan(par%al_h)
    select case (sign1)
       case (1)
          eta1 = (1.0_default,0.0_default)
       case (-1)
          eta1 = (0.0_default,1.0_default)
       case default 
          print *, 'sign1', sign1
          stop "parameters_MSSM: No definite sign neutralino1"
    end select
    select case (sign2)
       case (1)
          eta2 = (1.0_default,0.0_default)
       case (-1)
          eta2 = (0.0_default,1.0_default)
       case default 
          print *, 'sign2', sign2
          stop "parameters_MSSM: No definite sign neutralino2"
    end select
    select case (sign3)
       case (1)
          eta3 = (1.0_default,0.0_default)
       case (-1)
          eta3 = (0.0_default,1.0_default)
       case default 
          print *, 'sign3', sign3
          stop "parameters_MSSM: No definite sign neutralino3"
    end select
    select case (sign4)
       case (1)
          eta4 = (1.0_default,0.0_default)
       case (-1)
          eta4 = (0.0_default,1.0_default)
       case default 
          print *, 'sign4', sign4
          stop "parameters_MSSM: No definite sign neutralino4"
    end select
    sinthw = par%sw
    sin2thw = sinthw**2
    costhw = par%cw
    qelep = - 1.0_default
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
!  real(kind=default), public, save :: gs = 0_default, &
!    sind = 0_default, cosd = 1_default, sinckm12 = 0.223_default, &
!    sinckm13 = 0.004_default, sinckm23 = 0.04_default, &
!    tana = 30_default, tanb = 30_default 
!  complex(kind=default), public, save :: eidelta, cosckm23, cosckm13, &
!    cosckm12, vckm_33, vckm_32, vckm_31, vckm_23, vckm_22, vckm_21, vckm_13, &
!    vckm_12, vckm_11, gpzww, gppww, gzzww, gw4, igwww, igzww, iqw, igs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! The Planck mass as an input parameter
    m_planck = par%m_pl
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    call setup_parameters1
    call setup_parameters2
    call setup_parameters3
    call setup_parameters4
    call setup_parameters5
    call setup_parameters6
    call setup_parameters7
    call setup_parameters8
    call setup_parameters9
    call setup_parameters10
    call setup_parameters11
    call setup_parameters12
    call setup_parameters13
    call setup_parameters14
    call setup_parameters15
    call setup_parameters16
    call setup_parameters17
contains
 subroutine setup_parameters1 ()
    g = (e / sinthw)
    gz = (g / costhw)
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default * PI * par%alphas)
    igs = (imago * gs)
    q_lep  = (- 1.0_default) 
    q_up   = (2.0_default / 3.0_default)
    q_down = (- 1.0_default / 3.0_default)    
    qlep = - e * qelep   !!! This is the negative particle charge !!! 
    qup = - e * qeup     !!! This is the negative particle charge !!! 
    qdwn = - e * qedwn   !!! This is the negative particle charge !!! 
    qchar = ( - e)       !!! This is the negative particle charge !!! 
    ! qlep = ((-1.0_default) * e)
    ! qup = ((2.0_default / 3.0_default) * e)
    ! qdwn = (((-1.0_default) / 3.0_default) * e)
    gcc = (g / (2.0_default * sqrt (2.0_default)))
    gssq = (gs / sqrt (2.0_default))
    iqw = imago * e
    igzww = imago * g * costhw
    gw4 = (g**2)
    gzzww = ((g**2) * (costhw**2))
    gppww = (e**2)
    gpzww = (e * g * costhw) 
    sinal = sin (par%al_h)
    cosal = cos (par%al_h)
    sinbe = (tanb / sqrt ((1.0_default + (tanb**2))))
    cosbe = (1.0_default / sqrt ((1.0_default + (tanb**2))))
    eidelta = (cosd + (imago * sind))
    cos2be = ((cosbe**2) - (sinbe**2))
    cos2al = ((cosal**2) - (sinal**2))
    sin2be = (2.0_default * cosbe * sinbe)
    sin2al = (2.0_default * cosal * sinal)
    sin4al = (2.0_default * cos2al * sin2al)
    sin4be = (2.0_default * cos2be * sin2be)
    cos4be = ((cos2be**2) - (sin2be**2))
    cosapb = ((cosal * cosbe) - (sinal * sinbe))
    cosamb = ((cosal * cosbe) + (sinal * sinbe))
    sinapb = ((cosal * sinbe) + (sinal * cosbe))
    sinamb = ((sinal * cosbe) - (sinbe * cosal))
    sin2am2b = (2.0_default * sinamb * cosamb)
    cos2am2b = ((cosamb**2) - (sinamb**2))
    !!! costhsl1 = 1.0_default   !!! These are obsolete
    !!! costhsl2 = 1.0_default   !!! These are obsolete
    !!! costhsl3 = 1.0_default   !!! These are obsolete
    !!! costhsd1 = 1.0_default   !!! These are obsolete
    !!! costhsd2 = 1.0_default   !!! These are obsolete
    !!! costhsd3 = 1.0_default   !!! These are obsolete
    !!! costhsu1 = 1.0_default   !!! These are obsolete
    !!! costhsu2 = 1.0_default   !!! These are obsolete
    !!! costhsu3 = 1.0_default   !!! These are obsolete
    !!! sinthsl1 = 0.0_default   !!! These are obsolete
    !!! sinthsl2 = 0.0_default   !!! These are obsolete
    !!! sinthsl3 = 0.0_default   !!! These are obsolete
    !!! sinthsd1 = 0.0_default   !!! These are obsolete
    !!! sinthsd2 = 0.0_default   !!! These are obsolete
    !!! sinthsd3 = 0.0_default   !!! These are obsolete
    !!! sinthsu1 = 0.0_default   !!! These are obsolete
    !!! sinthsu2 = 0.0_default   !!! These are obsolete
    !!! sinthsu3 = 0.0_default   !!! These are obsolete
    mn_11 = eta1 * par%mn_11
    mn_12 = eta1 * par%mn_12
    mn_13 = eta1 * par%mn_13
    mn_14 = eta1 * par%mn_14
    mn_21 = eta2 * par%mn_21
    mn_22 = eta2 * par%mn_22
    mn_23 = eta2 * par%mn_23
    mn_24 = eta2 * par%mn_24
    mn_31 = eta3 * par%mn_31
    mn_32 = eta3 * par%mn_32
    mn_33 = eta3 * par%mn_33
    mn_34 = eta3 * par%mn_34
    mn_41 = eta4 * par%mn_41
    mn_42 = eta4 * par%mn_42
    mn_43 = eta4 * par%mn_43
    mn_44 = eta4 * par%mn_44
    !!! Checked by JR !!! 
    mu_11 = par%mu_11                !!! Rotat. matrix containing phi_R
    mu_12 = par%mu_12                !!! Rotat. matrix containing phi_R
    mu_21 = par%mu_21                !!! Rotat. matrix containing phi_R
    mu_22 = par%mu_22                !!! Rotat. matrix containing phi_R
    mv_11 = sigch1 * par%mv_11       !!! Rotat. matrix containing phi_L
    mv_12 = sigch1 * par%mv_12       !!! Rotat. matrix containing phi_L
    mv_21 = sigch2 * par%mv_21       !!! Rotat. matrix containing phi_L
    mv_22 = sigch2 * par%mv_22       !!! Rotat. matrix containing phi_L
    al_1 = 0
    au_1 = 0
    ad_1 = 0
    al_2 = 0
    au_2 = 0
    ad_2 = 0
    al_3 = par%Ae_33
    au_3 = par%Au_33
    ad_3 = par%Ad_33
    mu = par%mu_h
    mix_sl111 = 1.0_default
    mix_sl112 = 0.0_default
    mix_sl122 = 1.0_default
    mix_sl121 = 0.0_default
    mix_sl211 = 1.0_default
    mix_sl212 = 0.0_default
    mix_sl222 = 1.0_default
    mix_sl221 = 0.0_default
    mix_su111 = 1.0_default
    mix_su112 = 0.0_default
    mix_su122 = 1.0_default
    mix_su121 = 0.0_default
    mix_su211 = 1.0_default
    mix_su212 = 0.0_default
    mix_su222 = 1.0_default
    mix_su221 = 0.0_default
    mix_sd111 = 1.0_default
    mix_sd112 = 0.0_default
    mix_sd122 = 1.0_default
    mix_sd121 = 0.0_default
    mix_sd211 = 1.0_default
    mix_sd212 = 0.0_default
    mix_sd222 = 1.0_default
    mix_sd221 = 0.0_default
    !!! Checked by JR !!! 
    mix_sl311 = par%ml_11
    mix_sl312 = par%ml_12
    mix_sl321 = par%ml_21
    mix_sl322 = par%ml_22
    mix_su311 = par%mt_11
    mix_su312 = par%mt_12
    mix_su321 = par%mt_21
    mix_su322 = par%mt_22
    mix_sd311 = par%mb_11
    mix_sd312 = par%mb_12
    mix_sd321 = par%mb_21
    mix_sd322 = par%mb_22
    gh3_1 = ((mass(23) * (gz / 2.0_default) * cos2be * cosapb) -  &
      (mass(24) * g * cosamb))
    gh3_2 = ((mass(24) * g * sinamb) - ( &
      (gz / 2.0_default) * mass(23) * cos2be * sinapb))
    gh3_3 = ((gz / 2.0_default) * mass(23) * ( &
      (2.0_default * sin2al * cosapb) + (cos2al * sinapb)))
    gh3_4 = ( - ( &
      (3.0_default / 2.0_default) * gz * mass(23) * cos2al * cosapb))
    gh3_5 = ( - ( &
      (3.0_default / 2.0_default) * gz * mass(23) * cos2al * sinapb))
    gh3_6 = ((gz / 2.0_default) * mass(23) * ((cos2al * cosapb) -  &
      (2.0_default * sin2al * sinapb)))
    gh3_7 = ((gz / 2.0_default) * mass(23) * cos2be * cosapb)
    gh3_8 = ( - ((gz / 2.0_default) * mass(23) * cos2be * sinapb))
    gh4_1 = ( - (((gz**2) / 2.0_default) * (cos2be**2)))
 end subroutine setup_parameters1
 subroutine setup_parameters2 ()
    gh4_2 = ((((gz**2) / 4.0_default) * cos2al * cos2be) - (( &
      (g**2) / 2.0_default) * (cosamb**2)))
    gh4_3 = ( - ((((gz**2) / 4.0_default) * cos2al * cos2be) + (( &
      (g**2) / 2.0_default) * (sinamb**2))))
    gh4_4 = ((((g**2) / 2.0_default) * cosamb * sinamb) - (( &
      (gz**2) / 4.0_default) * sin2al * cos2be))
    gh4_5 = ( - (((gz**2) / 4.0_default) * (cos2be**2)))
    gh4_6 = ( - ((3.0_default / 4.0_default) * (gz**2) * (cos2al**2)))
    gh4_7 = (((gz**2) / 4.0_default) * (1.0_default - (3.0_default *  &
      (sin2al**2))))
    gh4_8 = ( - ((3.0_default / 8.0_default) * (gz**2) * sin4al))
    gh4_9 = (((gz**2) / 4.0_default) * cos2al * cos2be)
    gh4_10 = ( - (((gz**2) / 4.0_default) * sin2al * cos2be))
    gh4_11 = ( - ((3.0_default / 4.0_default) * (gz**2) * (cos2be**2)))
    ghaw = ( - (imago * (g / 2.0_default)))
    gh1az = (imago *  &
      (gz / 2.0_default) * cosamb)
    gh2az = (imago *  &
      (gz / 2.0_default) * sinamb)
    gh1ww = ( - (g * mass(24) * sinamb))
    gh2ww = (g * mass(24) * cosamb)
    ghh1w = ((g / 2.0_default) * cosamb)
    ghh2w = ((g / 2.0_default) * sinamb)
    gh1zz = ( - (gz * mass(23) * sinamb))
    gh2zz = (gz * mass(23) * cosamb)
    ghhz = ((gz / 2.0_default) * (1.0_default -  &
      (2.0_default * sin2thw)))
    ghhp = e
    gaazz = ((gz**2) / 2.0_default)
    gh1h1zz = gaazz
    gh2h2zz = gaazz
    ghphmzz = (gaazz * (((2.0_default * (costhw**2)) - 1.0_default)**2))
    ghphmpp = (2.0_default * (e**2))
    ghphmpz = (e * gz * ((2.0_default * (costhw**2)) - 1.0_default))
    ghh1wz = ( - ( &
      (1.0_default / 2.0_default) * g * gz * sin2thw * cosamb))
    ghh2wz = ( - ( &
      (1.0_default / 2.0_default) * g * gz * sin2thw * sinamb))
    ghh1wp = (e * (g / 2.0_default) * cosamb)
    ghh2wp = (e * (g / 2.0_default) * sinamb)
    gaaww = ((g**2) / 2.0_default)
    gh1h1ww = gaaww
    gh2h2ww = gaaww
    ghhww = gaaww
    ghawz = (imago * g * gz *  &
      (1.0_default / 2.0_default) * sin2thw)
    ghawp = ( - (imago * e * g *  &
      (1.0_default / 2.0_default)))
    g_yuk6_1 = (gcc * (mass(11) / mass(24)) * tanb)
    g_yuk7_1 = ( - ((g / 2.0_default) * (mass(11) / mass(24)) *  &
      (cosal / cosbe)))
    g_yuk8_1 = ((g / 2.0_default) * (mass(11) / mass(24)) * (sinal / cosbe))
    g_yuk9_1 = (imago *  &
      (g / 2.0_default) * (mass(11) / mass(24)) * tanb)
    g_yuk10_1 = ( - ((g / 2.0_default) * (mass(2) / mass(24)) *  &
      (sinal / sinbe)))
    g_yuk11_1 = ( - ((g / 2.0_default) * (mass(2) / mass(24)) *  &
      (cosal / sinbe)))
    g_yuk12_1 = (imago *  &
      (g / 2.0_default) * (mass(2) / mass(24)) * (1.0_default / tanb))
    g_yuk13_1 = ( - ((g / 2.0_default) * (mass(1) / mass(24)) *  &
      (cosal / cosbe)))
    g_yuk14_1 = ((g / 2.0_default) * (mass(1) / mass(24)) * (sinal / cosbe))
    g_yuk15_1 = (imago *  &
      (g / 2.0_default) * (mass(1) / mass(24)) * tanb)
    g_yuk6_2 = (gcc * (mass(13) / mass(24)) * tanb)
    g_yuk7_2 = ( - ((g / 2.0_default) * (mass(13) / mass(24)) *  &
      (cosal / cosbe)))
    g_yuk8_2 = ((g / 2.0_default) * (mass(13) / mass(24)) * (sinal / cosbe))
    g_yuk9_2 = (imago *  &
      (g / 2.0_default) * (mass(13) / mass(24)) * tanb)
    g_yuk10_2 = ( - ((g / 2.0_default) * (mass(4) / mass(24)) *  &
      (sinal / sinbe)))
    g_yuk11_2 = ( - ((g / 2.0_default) * (mass(4) / mass(24)) *  &
      (cosal / sinbe)))
    g_yuk12_2 = (imago *  &
      (g / 2.0_default) * (mass(4) / mass(24)) * (1.0_default / tanb))
    g_yuk13_2 = ( - ((g / 2.0_default) * (mass(3) / mass(24)) *  &
      (cosal / cosbe)))
    g_yuk14_2 = ((g / 2.0_default) * (mass(3) / mass(24)) * (sinal / cosbe))
    g_yuk15_2 = (imago *  &
      (g / 2.0_default) * (mass(3) / mass(24)) * tanb)
    g_yuk6_3 = (gcc * (mass(15) / mass(24)) * tanb)
    g_yuk7_3 = ( - ((g / 2.0_default) * (mass(15) / mass(24)) *  &
      (cosal / cosbe)))
    g_yuk8_3 = ((g / 2.0_default) * (mass(15) / mass(24)) * (sinal / cosbe))
    g_yuk9_3 = (imago *  &
      (g / 2.0_default) * (mass(15) / mass(24)) * tanb)
    g_yuk10_3 = ( - ((g / 2.0_default) * (mass(6) / mass(24)) *  &
      (sinal / sinbe)))
    g_yuk11_3 = ( - ((g / 2.0_default) * (mass(6) / mass(24)) *  &
      (cosal / sinbe)))
    g_yuk12_3 = (imago *  &
      (g / 2.0_default) * (mass(6) / mass(24)) * (1.0_default / tanb))
    g_yuk13_3 = ( - ((g / 2.0_default) * (mass(5) / mass(24)) *  &
      (cosal / cosbe)))
    g_yuk14_3 = ((g / 2.0_default) * (mass(5) / mass(24)) * (sinal / cosbe))
    g_yuk15_3 = (imago *  &
      (g / 2.0_default) * (mass(5) / mass(24)) * tanb)
    gccq_1_1 = (gcc * vckm_11)
    gccq_1_1_c = (gcc * conjg (vckm_11))
    gccq_1_2 = (gcc * vckm_12)
    gccq_1_2_c = (gcc * conjg (vckm_12))
    gccq_1_3 = (gcc * vckm_13)
    gccq_1_3_c = (gcc * conjg (vckm_13))
    gccq_2_1 = (gcc * vckm_21)
    gccq_2_1_c = (gcc * conjg (vckm_21))
    gccq_2_2 = (gcc * vckm_22)
    gccq_2_2_c = (gcc * conjg (vckm_22))
    gccq_2_3 = (gcc * vckm_23)
    gccq_2_3_c = (gcc * conjg (vckm_23))
    gccq_3_1 = (gcc * vckm_31)
    gccq_3_1_c = (gcc * conjg (vckm_31))
    gccq_3_2 = (gcc * vckm_32)
    gccq_3_2_c = (gcc * conjg (vckm_32))
    gccq_3_3 = (gcc * vckm_33)
    gccq_3_3_c = (gcc * conjg (vckm_33))
    gs1ws1_1_1 = ( - (gcc * 2.0_default * vckm_11 *  &
      conjg (mix_su111) * mix_sd111))
    gs2ws2_1_1 = ( - (gcc * 2.0_default * vckm_11 *  &
      conjg (mix_su121) * mix_sd121))
    gs1ws2_1_1 = ( - (gcc * 2.0_default * vckm_11 *  &
      conjg (mix_su111) * mix_sd121))
    gs2ws1_1_1 = ( - (gcc * 2.0_default * vckm_11 *  &
      conjg (mix_su121) * mix_sd111))
    gs1ws1_1_1_c = conjg (gs1ws1_1_1) 
    gs2ws2_1_1_c = conjg (gs2ws2_1_1)
    gs1ws2_1_1_c = conjg (gs1ws2_1_1)
    gs2ws1_1_1_c = conjg (gs2ws1_1_1)
    gs1ws1_1_2 = ( - (gcc * 2.0_default * vckm_12 *  &
      conjg (mix_su111) * mix_sd211))
    gs2ws2_1_2 = ( - (gcc * 2.0_default * vckm_12 *  &
      conjg (mix_su121) * mix_sd221))
    gs1ws2_1_2 = ( - (gcc * 2.0_default * vckm_12 *  &
      conjg (mix_su111) * mix_sd221))
    gs2ws1_1_2 = ( - (gcc * 2.0_default * vckm_12 *  &
      conjg (mix_su121) * mix_sd211))
    gs1ws1_1_2_c = conjg (gs1ws1_1_2)
    gs2ws2_1_2_c = conjg (gs2ws2_1_2)
    gs1ws2_1_2_c = conjg (gs1ws2_1_2)
    gs2ws1_1_2_c = conjg (gs2ws1_1_2)
    gs1ws1_1_3 = ( - (gcc * 2.0_default * vckm_13 *  &
      conjg (mix_su111) * mix_sd311))
    gs2ws2_1_3 = ( - (gcc * 2.0_default * vckm_13 *  &
      conjg (mix_su121) * mix_sd321))
    gs1ws2_1_3 = ( - (gcc * 2.0_default * vckm_13 *  &
      conjg (mix_su111) * mix_sd321))
    gs2ws1_1_3 = ( - (gcc * 2.0_default * vckm_13 *  &
      conjg (mix_su121) * mix_sd311))
    gs1ws1_1_3_c = conjg (gs1ws1_1_3) 
    gs2ws2_1_3_c = conjg (gs2ws2_1_3)
    gs1ws2_1_3_c = conjg (gs1ws2_1_3)
    gs2ws1_1_3_c = conjg (gs2ws1_1_3)
    gs1ws1_2_1 = ( - (gcc * 2.0_default * vckm_21 *  &
      conjg (mix_su211) * mix_sd111))
    gs2ws2_2_1 = ( - (gcc * 2.0_default * vckm_21 *  &
      conjg (mix_su221) * mix_sd121))
    gs1ws2_2_1 = ( - (gcc * 2.0_default * vckm_21 *  &
      conjg (mix_su211) * mix_sd121))
    gs2ws1_2_1 = ( - (gcc * 2.0_default * vckm_21 *  &
      conjg (mix_su221) * mix_sd111))
    gs1ws1_2_1_c = conjg (gs1ws1_2_1)  
    gs2ws2_2_1_c = conjg (gs2ws2_2_1)
    gs1ws2_2_1_c = conjg (gs1ws2_2_1)
    gs2ws1_2_1_c = conjg (gs2ws1_2_1)
    gs1ws1_2_2 = ( - (gcc * 2.0_default * vckm_22 *  &
      conjg (mix_su211) * mix_sd211))
    gs2ws2_2_2 = ( - (gcc * 2.0_default * vckm_22 *  &
      conjg (mix_su221) * mix_sd221))
    gs1ws2_2_2 = ( - (gcc * 2.0_default * vckm_22 *  &
      conjg (mix_su211) * mix_sd221))
    gs2ws1_2_2 = ( - (gcc * 2.0_default * vckm_22 *  &
      conjg (mix_su221) * mix_sd211))
    gs1ws1_2_2_c = conjg (gs1ws1_2_2)
    gs2ws2_2_2_c = conjg (gs2ws2_2_2)
    gs1ws2_2_2_c = conjg (gs1ws2_2_2)
    gs2ws1_2_2_c = conjg (gs2ws1_2_2)
    gs1ws1_2_3 = ( - (gcc * 2.0_default * vckm_23 *  &
      conjg (mix_su211) * mix_sd311))
    gs2ws2_2_3 = ( - (gcc * 2.0_default * vckm_23 *  &
      conjg (mix_su221) * mix_sd321))
    gs1ws2_2_3 = ( - (gcc * 2.0_default * vckm_23 *  &
      conjg (mix_su211) * mix_sd321))
    gs2ws1_2_3 = ( - (gcc * 2.0_default * vckm_23 *  &
      conjg (mix_su221) * mix_sd311))
    gs1ws1_2_3_c = conjg (gs1ws1_2_3)
    gs2ws2_2_3_c = conjg (gs2ws2_2_3)
    gs1ws2_2_3_c = conjg (gs1ws2_2_3)
    gs2ws1_2_3_c = conjg (gs2ws1_2_3)
    gs1ws1_3_1 = ( - (gcc * 2.0_default * vckm_31 *  &
      conjg (mix_su311) * mix_sd111))
    gs2ws2_3_1 = ( - (gcc * 2.0_default * vckm_31 *  &
      conjg (mix_su321) * mix_sd121))
 end subroutine setup_parameters2
 subroutine setup_parameters3 ()
    gs1ws2_3_1 = ( - (gcc * 2.0_default * vckm_31 *  &
      conjg (mix_su311) * mix_sd121))
    gs2ws1_3_1 = ( - (gcc * 2.0_default * vckm_31 *  &
      conjg (mix_su321) * mix_sd111))
    gs1ws1_3_1_c = conjg (gs1ws1_3_1) 
    gs2ws2_3_1_c = conjg (gs2ws2_3_1)
    gs1ws2_3_1_c = conjg (gs1ws2_3_1)
    gs2ws1_3_1_c = conjg (gs2ws1_3_1)
    gs1ws1_3_2 = ( - (gcc * 2.0_default * vckm_32 *  &
      conjg (mix_su311) * mix_sd211))
    gs2ws2_3_2 = ( - (gcc * 2.0_default * vckm_32 *  &
      conjg (mix_su321) * mix_sd221))
    gs1ws2_3_2 = ( - (gcc * 2.0_default * vckm_32 *  &
      conjg (mix_su311) * mix_sd221))
    gs2ws1_3_2 = ( - (gcc * 2.0_default * vckm_32 *  &
      conjg (mix_su321) * mix_sd211))
    gs1ws1_3_2_c = conjg (gs1ws1_3_2) 
    gs2ws2_3_2_c = conjg (gs2ws2_3_2)
    gs1ws2_3_2_c = conjg (gs1ws2_3_2)
    gs2ws1_3_2_c = conjg (gs2ws1_3_2)
    gs1ws1_3_3 = ( - (gcc * 2.0_default * vckm_33 *  &
      conjg (mix_su311) * mix_sd311))
    gs2ws2_3_3 = ( - (gcc * 2.0_default * vckm_33 *  &
      conjg (mix_su321) * mix_sd321))
    gs1ws2_3_3 = ( - (gcc * 2.0_default * vckm_33 *  &
      conjg (mix_su311) * mix_sd321))
    gs2ws1_3_3 = ( - (gcc * 2.0_default * vckm_33 *  &
      conjg (mix_su321) * mix_sd311))
    gs1ws1_3_3_c = conjg (gs1ws1_3_3)
    gs2ws2_3_3_c = conjg (gs2ws2_3_3)
    gs1ws2_3_3_c = conjg (gs1ws2_3_3)
    gs2ws1_3_3_c = conjg (gs2ws1_3_3)
    gsl1zsl1_1 = ((gz / 2.0_default) * ((2.0_default * sin2thw) -  &
      (mix_sl111 * conjg (mix_sl111))))
    gsl2zsl2_1 = ((gz / 2.0_default) * ((2.0_default * sin2thw) -  &
      (mix_sl121 * conjg (mix_sl121))))
    gsl1zsl2_1 = ((( - gz) / 2.0_default) * conjg (mix_sl111) * mix_sl121)
    gsl2zsl1_1 = conjg (gsl1zsl2_1)
    gsn1zsn1_1 = (gz / 2.0_default)
    gsu1zsu1_1 = ((gz / 2.0_default) * ((mix_su111 * conjg (mix_su111)) - ( &
      (4.0_default / 3.0_default) * sin2thw)))
    gsu2zsu2_1 = ((gz / 2.0_default) * ((mix_su121 * conjg (mix_su121)) - ( &
      (4.0_default / 3.0_default) * sin2thw)))
    gsu1zsu2_1 = ((gz / 2.0_default) * conjg (mix_su111) * mix_su121)
    gsu2zsu1_1 = conjg (gsu1zsu2_1)
    gsd1zsd1_1 = ((gz / 2.0_default) * (( &
      (2.0_default / 3.0_default) * sin2thw) - (mix_sd111 *  &
      conjg (mix_sd111))))
    gsd2zsd2_1 = ((gz / 2.0_default) * (( &
      (2.0_default / 3.0_default) * sin2thw) - (mix_sd121 *  &
      conjg (mix_sd121))))
    gsd1zsd2_1 = ((( - gz) / 2.0_default) * conjg (mix_sd111) * mix_sd121)
    gsd2zsd1_1 = conjg (gsd1zsd2_1)
    gsl1_1snw = (gcc * 2.0_default * mix_sl111)
    gsl2_1snw = (gcc * 2.0_default * mix_sl121)
    gsl1_1snw_c = (gcc * 2.0_default * conjg (mix_sl111))
    gsl2_1snw_c = (gcc * 2.0_default * conjg (mix_sl121))
    gsl1zsl1_2 = ((gz / 2.0_default) * ((2.0_default * sin2thw) -  &
      (mix_sl211 * conjg (mix_sl211))))
    gsl2zsl2_2 = ((gz / 2.0_default) * ((2.0_default * sin2thw) -  &
      (mix_sl221 * conjg (mix_sl221))))
    gsl1zsl2_2 = ((( - gz) / 2.0_default) * conjg (mix_sl211) * mix_sl221)
    gsl2zsl1_2 = conjg (gsl1zsl2_2)
    gsn1zsn1_2 = (gz / 2.0_default)
    gsu1zsu1_2 = ((gz / 2.0_default) * ((mix_su211 * conjg (mix_su211)) - ( &
      (4.0_default / 3.0_default) * sin2thw)))
    gsu2zsu2_2 = ((gz / 2.0_default) * ((mix_su221 * conjg (mix_su221)) - ( &
      (4.0_default / 3.0_default) * sin2thw)))
    gsu1zsu2_2 = ((gz / 2.0_default) * conjg (mix_su211) * mix_su221)
    gsu2zsu1_2 = conjg (gsu1zsu2_2)
    gsd1zsd1_2 = ((gz / 2.0_default) * (( &
      (2.0_default / 3.0_default) * sin2thw) - (mix_sd211 *  &
      conjg (mix_sd211))))
    gsd2zsd2_2 = ((gz / 2.0_default) * (( &
      (2.0_default / 3.0_default) * sin2thw) - (mix_sd221 *  &
      conjg (mix_sd221))))
    gsd1zsd2_2 = ((( - gz) / 2.0_default) * conjg (mix_sd211) * mix_sd221)
    gsd2zsd1_2 = conjg (gsd1zsd2_2)
    gsl1_2snw = (gcc * 2.0_default * mix_sl211)
    gsl2_2snw = (gcc * 2.0_default * mix_sl221)
    gsl1_2snw_c = (gcc * 2.0_default * conjg (mix_sl211))
    gsl2_2snw_c = (gcc * 2.0_default * conjg (mix_sl221))
    gsl1zsl1_3 = ((gz / 2.0_default) * ((2.0_default * sin2thw) -  &
      (mix_sl311 * conjg (mix_sl311))))
    gsl2zsl2_3 = ((gz / 2.0_default) * ((2.0_default * sin2thw) -  &
      (mix_sl321 * conjg (mix_sl321))))
    gsl1zsl2_3 = ((( - gz) / 2.0_default) * conjg (mix_sl311) * mix_sl321)
    gsl2zsl1_3 = conjg (gsl1zsl2_3)
    gsn1zsn1_3 = (gz / 2.0_default)
    gsu1zsu1_3 = ((gz / 2.0_default) * ((mix_su311 * conjg (mix_su311)) - ( &
      (4.0_default / 3.0_default) * sin2thw)))
    gsu2zsu2_3 = ((gz / 2.0_default) * ((mix_su321 * conjg (mix_su321)) - ( &
      (4.0_default / 3.0_default) * sin2thw)))
    gsu1zsu2_3 = ((gz / 2.0_default) * conjg (mix_su311) * mix_su321)
    gsu2zsu1_3 = conjg (gsu1zsu2_3)
    gsd1zsd1_3 = ((gz / 2.0_default) * (( &
      (2.0_default / 3.0_default) * sin2thw) - (mix_sd311 *  &
      conjg (mix_sd311))))
    gsd2zsd2_3 = ((gz / 2.0_default) * (( &
      (2.0_default / 3.0_default) * sin2thw) - (mix_sd321 *  &
      conjg (mix_sd321))))
    gsd1zsd2_3 = ((( - gz) / 2.0_default) * conjg (mix_sd311) * mix_sd321)
    gsd2zsd1_3 = conjg (gsd1zsd2_3)
    gsl1_3snw = (gcc * 2.0_default * mix_sl311)
    gsl2_3snw = (gcc * 2.0_default * mix_sl321)
    gsl1_3snw_c = (gcc * 2.0_default * conjg (mix_sl311))
    gsl2_3snw_c = (gcc * 2.0_default * conjg (mix_sl321))
    gppslsl = (2.0_default * (e**2))
    gppsusu = ((8.0_default / 9.0_default) * (e**2))
    gppsdsd = ((2.0_default / 9.0_default) * (e**2))
    gzzsl1sl1_1 = (((gz**2) / 2.0_default) * (((1.0_default -  &
      (4.0_default * sin2thw)) * (mix_sl111 * conjg (mix_sl111))) +  &
      (4.0_default * (sin2thw**2))))
    gzzsl2sl2_1 = (((gz**2) / 2.0_default) * (((1.0_default -  &
      (4.0_default * sin2thw)) * (mix_sl121 * conjg (mix_sl121))) +  &
      (4.0_default * (sin2thw**2))))
    gzzsl1sl2_1 = (((gz**2) / 2.0_default) * (1.0_default -  &
      (4.0_default * sin2thw)) * mix_sl111 * conjg (mix_sl121))
    gzzsl2sl1_1 = conjg(gzzsl1sl2_1)
    gzzsn1sn1_1 = ((gz**2) / 2.0_default)
    gzzsu1su1_1 = (((gz**2) / 2.0_default) * (((1.0_default - ( &
      (8.0_default / 3.0_default) * sin2thw)) * (mix_su111 *  &
      conjg (mix_su111))) + ((sin2thw**2) *  &
      (16.0_default / 9.0_default))))
    gzzsu2su2_1 = (((gz**2) / 2.0_default) * (((1.0_default - ( &
      (8.0_default / 3.0_default) * sin2thw)) * (mix_su121 *  &
      conjg (mix_su121))) + ((sin2thw**2) *  &
      (16.0_default / 9.0_default))))
    gzzsu1su2_1 = (((gz**2) / 2.0_default) * (1.0_default - (sin2thw *  &
      (8.0_default / 3.0_default))) * mix_su111 * conjg (mix_su121))
    gzzsu2su1_1 = conjg(gzzsu1su2_1)
    gzzsd1sd1_1 = (((gz**2) / 2.0_default) * (((1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * (mix_sd111 * conjg (mix_sd111))) +  &
      ((sin2thw**2) * (4.0_default / 9.0_default))))
    gzzsd2sd2_1 = (((gz**2) / 2.0_default) * (((1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * (mix_sd121 * conjg (mix_sd121))) +  &
      ((sin2thw**2) * (4.0_default / 9.0_default))))
    gzzsd1sd2_1 = (((gz**2) / 2.0_default) * (1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * mix_sd111 * conjg (mix_sd121))
    gzzsd2sd1_1 = conjg(gzzsd1sd2_1)
    gzpsl1sl1_1 = (e * gz * ((mix_sl111 * conjg (mix_sl111)) -  &
      (2.0_default * sin2thw)))
    gzpsl2sl2_1 = (e * gz * ((mix_sl121 * conjg (mix_sl121)) -  &
      (2.0_default * sin2thw)))
    gzpsl1sl2_1 = (e * gz * mix_sl111 * conjg (mix_sl121))
    gzpsl2sl1_1 = (e * gz * mix_sl121 * conjg (mix_sl111))
    gzpsu1su1_1 = (e * gz * (2.0_default / 3.0_default) * ((mix_su111 *  &
      conjg (mix_su111)) - (sin2thw * (4.0_default / 3.0_default))))
    gzpsu2su2_1 = (e * gz * (2.0_default / 3.0_default) * ((mix_su121 *  &
      conjg (mix_su121)) - (sin2thw * (4.0_default / 3.0_default))))
    gzpsu1su2_1 = (e * gz * (2.0_default / 3.0_default) * mix_su111 *  &
      conjg (mix_su121))
    gzpsu2su1_1 = (e * gz * (2.0_default / 3.0_default) * mix_su121 *  &
      conjg (mix_su111))
    gzpsd1sd1_1 = (e * gz * (1.0_default / 3.0_default) * ((mix_sd111 *  &
      conjg (mix_sd111)) - (sin2thw * (2.0_default / 3.0_default))))
    gzpsd2sd2_1 = (e * gz * (1.0_default / 3.0_default) * ((mix_sd121 *  &
      conjg (mix_sd121)) - (sin2thw * (2.0_default / 3.0_default))))
    gzpsd1sd2_1 = (e * gz * (1.0_default / 3.0_default) * mix_sd111 *  &
      conjg (mix_sd121))
    gzpsd2sd1_1 = (e * gz * (1.0_default / 3.0_default) * mix_sd121 *  &
      conjg (mix_sd111))
    gwwsl1sl1_1 = (((g**2) / 2.0_default) * (mix_sl111 * conjg (mix_sl111)))
    gwwsl2sl2_1 = (((g**2) / 2.0_default) * (mix_sl121 * conjg (mix_sl121)))
    gwwsl1sl2_1 = (((g**2) / 2.0_default) * mix_sl111 * conjg (mix_sl121))
    gwwsl2sl1_1 = (((g**2) / 2.0_default) * mix_sl121 * conjg (mix_sl111))
    gwwsn1sn1_1 = ((g**2) / 2.0_default)
    gwwsu1su1_1 = (((g**2) / 2.0_default) * (mix_su111 * conjg (mix_su111)))
    gwwsu2su2_1 = (((g**2) / 2.0_default) * (mix_su121 * conjg (mix_su121)))
    gwwsu1su2_1 = (((g**2) / 2.0_default) * mix_su111 * conjg (mix_su121))
    gwwsu2su1_1 = (((g**2) / 2.0_default) * mix_su121 * conjg (mix_su111))
    gwwsd1sd1_1 = (((g**2) / 2.0_default) * (mix_sd111 * conjg (mix_sd111)))
    gwwsd2sd2_1 = (((g**2) / 2.0_default) * (mix_sd121 * conjg (mix_sd121)))
    gwwsd1sd2_1 = (((g**2) / 2.0_default) * mix_sd111 * conjg (mix_sd121))
    gwwsd2sd1_1 = (((g**2) / 2.0_default) * mix_sd121 * conjg (mix_sd111))
    gpwsl1sn_1 = ( - (e * 2.0_default * gcc * mix_sl111))
    gpwsl2sn_1 = ( - (e * 2.0_default * gcc * mix_sl121))
    gpwsl1sn_1_c = ( - (e * 2.0_default * gcc * conjg (mix_sl111)))
    gpwsl2sn_1_c = ( - (e * 2.0_default * gcc * conjg (mix_sl121)))
    gwzsl1sn_1 = (gcc * gz * 2.0_default * sin2thw * mix_sl111)
 end subroutine setup_parameters3
 subroutine setup_parameters4 ()
    gwzsl2sn_1 = (gcc * gz * 2.0_default * sin2thw * mix_sl121)
    gwzsl1sn_1_c = (gcc * gz * 2.0_default * sin2thw * conjg (mix_sl111))
    gwzsl2sn_1_c = (gcc * gz * 2.0_default * sin2thw * conjg (mix_sl121))
    gzzsl1sl1_2 = (((gz**2) / 2.0_default) * (((1.0_default -  &
      (4.0_default * sin2thw)) * (mix_sl211 * conjg (mix_sl211))) +  &
      (4.0_default * (sin2thw**2))))
    gzzsl2sl2_2 = (((gz**2) / 2.0_default) * (((1.0_default -  &
      (4.0_default * sin2thw)) * (mix_sl221 * conjg (mix_sl221))) +  &
      (4.0_default * (sin2thw**2))))
    gzzsl1sl2_2 = (((gz**2) / 2.0_default) * (1.0_default -  &
      (4.0_default * sin2thw)) * mix_sl211 * conjg (mix_sl221))
    gzzsl2sl1_2 = conjg(gzzsl1sl2_2)
    gzzsn1sn1_2 = ((gz**2) / 2.0_default)
    gzzsu1su1_2 = (((gz**2) / 2.0_default) * (((1.0_default - ( &
      (8.0_default / 3.0_default) * sin2thw)) * (mix_su211 *  &
      conjg (mix_su211))) + ((sin2thw**2) *  &
      (16.0_default / 9.0_default))))
    gzzsu2su2_2 = (((gz**2) / 2.0_default) * (((1.0_default - ( &
      (8.0_default / 3.0_default) * sin2thw)) * (mix_su221 *  &
      conjg (mix_su221))) + ((sin2thw**2) *  &
      (16.0_default / 9.0_default))))
    gzzsu1su2_2 = (((gz**2) / 2.0_default) * (1.0_default - (sin2thw *  &
      (8.0_default / 3.0_default))) * mix_su211 * conjg (mix_su221))
    gzzsu2su1_2 = conjg(gzzsu1su2_2)
    gzzsd1sd1_2 = (((gz**2) / 2.0_default) * (((1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * (mix_sd211 * conjg (mix_sd211))) +  &
      ((sin2thw**2) * (4.0_default / 9.0_default))))
    gzzsd2sd2_2 = (((gz**2) / 2.0_default) * (((1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * (mix_sd221 * conjg (mix_sd221))) +  &
      ((sin2thw**2) * (4.0_default / 9.0_default))))
    gzzsd1sd2_2 = (((gz**2) / 2.0_default) * (1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * mix_sd211 * conjg (mix_sd221))
    gzzsd2sd1_2 = conjg(gzzsd1sd2_2)
    gzpsl1sl1_2 = (e * gz * ((mix_sl211 * conjg (mix_sl211)) -  &
      (2.0_default * sin2thw)))
    gzpsl2sl2_2 = (e * gz * ((mix_sl221 * conjg (mix_sl221)) -  &
      (2.0_default * sin2thw)))
    gzpsl1sl2_2 = (e * gz * mix_sl211 * conjg (mix_sl221))
    gzpsl2sl1_2 = (e * gz * mix_sl221 * conjg (mix_sl211))
    gzpsu1su1_2 = (e * gz * (2.0_default / 3.0_default) * ((mix_su211 *  &
      conjg (mix_su211)) - (sin2thw * (4.0_default / 3.0_default))))
    gzpsu2su2_2 = (e * gz * (2.0_default / 3.0_default) * ((mix_su221 *  &
      conjg (mix_su221)) - (sin2thw * (4.0_default / 3.0_default))))
    gzpsu1su2_2 = (e * gz * (2.0_default / 3.0_default) * mix_su211 *  &
      conjg (mix_su221))
    gzpsu2su1_2 = (e * gz * (2.0_default / 3.0_default) * mix_su221 *  &
      conjg (mix_su211))
    gzpsd1sd1_2 = (e * gz * (1.0_default / 3.0_default) * ((mix_sd211 *  &
      conjg (mix_sd211)) - (sin2thw * (2.0_default / 3.0_default))))
    gzpsd2sd2_2 = (e * gz * (1.0_default / 3.0_default) * ((mix_sd221 *  &
      conjg (mix_sd221)) - (sin2thw * (2.0_default / 3.0_default))))
    gzpsd1sd2_2 = (e * gz * (1.0_default / 3.0_default) * mix_sd211 *  &
      conjg (mix_sd221))
    gzpsd2sd1_2 = (e * gz * (1.0_default / 3.0_default) * mix_sd221 *  &
      conjg (mix_sd211))
    gwwsl1sl1_2 = (((g**2) / 2.0_default) * (mix_sl211 * conjg (mix_sl211)))
    gwwsl2sl2_2 = (((g**2) / 2.0_default) * (mix_sl221 * conjg (mix_sl221)))
    gwwsl1sl2_2 = (((g**2) / 2.0_default) * mix_sl211 * conjg (mix_sl221))
    gwwsl2sl1_2 = (((g**2) / 2.0_default) * mix_sl221 * conjg (mix_sl211))
    gwwsn1sn1_2 = ((g**2) / 2.0_default)
    gwwsu1su1_2 = (((g**2) / 2.0_default) * (mix_su211 * conjg (mix_su211)))
    gwwsu2su2_2 = (((g**2) / 2.0_default) * (mix_su221 * conjg (mix_su221)))
    gwwsu1su2_2 = (((g**2) / 2.0_default) * mix_su211 * conjg (mix_su221))
    gwwsu2su1_2 = (((g**2) / 2.0_default) * mix_su221 * conjg (mix_su211))
    gwwsd1sd1_2 = (((g**2) / 2.0_default) * (mix_sd211 * conjg (mix_sd211)))
    gwwsd2sd2_2 = (((g**2) / 2.0_default) * (mix_sd221 * conjg (mix_sd221)))
    gwwsd1sd2_2 = (((g**2) / 2.0_default) * mix_sd211 * conjg (mix_sd221))
    gwwsd2sd1_2 = (((g**2) / 2.0_default) * mix_sd221 * conjg (mix_sd211))
    gpwsl1sn_2 = ( - (e * 2.0_default * gcc * mix_sl211))
    gpwsl2sn_2 = ( - (e * 2.0_default * gcc * mix_sl221))
    gpwsl1sn_2_c = ( - (e * 2.0_default * gcc * conjg (mix_sl211)))
    gpwsl2sn_2_c = ( - (e * 2.0_default * gcc * conjg (mix_sl221)))
    gwzsl1sn_2 = (gcc * gz * 2.0_default * sin2thw * mix_sl211)
    gwzsl2sn_2 = (gcc * gz * 2.0_default * sin2thw * mix_sl221)
    gwzsl1sn_2_c = (gcc * gz * 2.0_default * sin2thw * conjg (mix_sl211))
    gwzsl2sn_2_c = (gcc * gz * 2.0_default * sin2thw * conjg (mix_sl221))
    gzzsl1sl1_3 = (((gz**2) / 2.0_default) * (((1.0_default -  &
      (4.0_default * sin2thw)) * (mix_sl311 * conjg (mix_sl311))) +  &
      (4.0_default * (sin2thw**2))))
    gzzsl2sl2_3 = (((gz**2) / 2.0_default) * (((1.0_default -  &
      (4.0_default * sin2thw)) * (mix_sl321 * conjg (mix_sl321))) +  &
      (4.0_default * (sin2thw**2))))
    gzzsl1sl2_3 = (((gz**2) / 2.0_default) * (1.0_default -  &
      (4.0_default * sin2thw)) * mix_sl311 * conjg (mix_sl321))
    gzzsl2sl1_3 = conjg(gzzsl1sl2_3)
    gzzsn1sn1_3 = ((gz**2) / 2.0_default)
    gzzsu1su1_3 = (((gz**2) / 2.0_default) * (((1.0_default - ( &
      (8.0_default / 3.0_default) * sin2thw)) * (mix_su311 *  &
      conjg (mix_su311))) + ((sin2thw**2) *  &
      (16.0_default / 9.0_default))))
    gzzsu2su2_3 = (((gz**2) / 2.0_default) * (((1.0_default - ( &
      (8.0_default / 3.0_default) * sin2thw)) * (mix_su321 *  &
      conjg (mix_su321))) + ((sin2thw**2) *  &
      (16.0_default / 9.0_default))))
    gzzsu1su2_3 = (((gz**2) / 2.0_default) * (1.0_default - (sin2thw *  &
      (8.0_default / 3.0_default))) * mix_su311 * conjg (mix_su321))
    gzzsu2su1_3 = conjg(gzzsu1su2_3)
    gzzsd1sd1_3 = (((gz**2) / 2.0_default) * (((1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * (mix_sd311 * conjg (mix_sd311))) +  &
      ((sin2thw**2) * (4.0_default / 9.0_default))))
    gzzsd2sd2_3 = (((gz**2) / 2.0_default) * (((1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * (mix_sd321 * conjg (mix_sd321))) +  &
      ((sin2thw**2) * (4.0_default / 9.0_default))))
    gzzsd1sd2_3 = (((gz**2) / 2.0_default) * (1.0_default - (sin2thw *  &
      (4.0_default / 3.0_default))) * mix_sd311 * conjg (mix_sd321))
    gzzsd2sd1_3 = conjg(gzzsd1sd2_3)
    gzpsl1sl1_3 = (e * gz * ((mix_sl311 * conjg (mix_sl311)) -  &
      (2.0_default * sin2thw)))
    gzpsl2sl2_3 = (e * gz * ((mix_sl321 * conjg (mix_sl321)) -  &
      (2.0_default * sin2thw)))
    gzpsl1sl2_3 = (e * gz * mix_sl311 * conjg (mix_sl321))
    gzpsl2sl1_3 = (e * gz * mix_sl321 * conjg (mix_sl311))
    gzpsu1su1_3 = (e * gz * (2.0_default / 3.0_default) * ((mix_su311 *  &
      conjg (mix_su311)) - (sin2thw * (4.0_default / 3.0_default))))
    gzpsu2su2_3 = (e * gz * (2.0_default / 3.0_default) * ((mix_su321 *  &
      conjg (mix_su321)) - (sin2thw * (4.0_default / 3.0_default))))
    gzpsu1su2_3 = (e * gz * (2.0_default / 3.0_default) * mix_su311 *  &
      conjg (mix_su321))
    gzpsu2su1_3 = (e * gz * (2.0_default / 3.0_default) * mix_su321 *  &
      conjg (mix_su311))
    gzpsd1sd1_3 = (e * gz * (1.0_default / 3.0_default) * ((mix_sd311 *  &
      conjg (mix_sd311)) - (sin2thw * (2.0_default / 3.0_default))))
    gzpsd2sd2_3 = (e * gz * (1.0_default / 3.0_default) * ((mix_sd321 *  &
      conjg (mix_sd321)) - (sin2thw * (2.0_default / 3.0_default))))
    gzpsd1sd2_3 = (e * gz * (1.0_default / 3.0_default) * mix_sd311 *  &
      conjg (mix_sd321))
    gzpsd2sd1_3 = (e * gz * (1.0_default / 3.0_default) * mix_sd321 *  &
      conjg (mix_sd311))
    gwwsl1sl1_3 = (((g**2) / 2.0_default) * (mix_sl311 * conjg (mix_sl311)))
    gwwsl2sl2_3 = (((g**2) / 2.0_default) * (mix_sl321 * conjg (mix_sl321)))
    gwwsl1sl2_3 = (((g**2) / 2.0_default) * mix_sl311 * conjg (mix_sl321))
    gwwsl2sl1_3 = (((g**2) / 2.0_default) * mix_sl321 * conjg (mix_sl311))
    gwwsn1sn1_3 = ((g**2) / 2.0_default)
    gwwsu1su1_3 = (((g**2) / 2.0_default) * (mix_su311 * conjg (mix_su311)))
    gwwsu2su2_3 = (((g**2) / 2.0_default) * (mix_su321 * conjg (mix_su321)))
    gwwsu1su2_3 = (((g**2) / 2.0_default) * mix_su311 * conjg (mix_su321))
    gwwsu2su1_3 = (((g**2) / 2.0_default) * mix_su321 * conjg (mix_su311))
    gwwsd1sd1_3 = (((g**2) / 2.0_default) * (mix_sd311 * conjg (mix_sd311)))
    gwwsd2sd2_3 = (((g**2) / 2.0_default) * (mix_sd321 * conjg (mix_sd321)))
    gwwsd1sd2_3 = (((g**2) / 2.0_default) * mix_sd311 * conjg (mix_sd321))
    gwwsd2sd1_3 = (((g**2) / 2.0_default) * mix_sd321 * conjg (mix_sd311))
    gpwsl1sn_3 = ( - (e * 2.0_default * gcc * mix_sl311))
    gpwsl2sn_3 = ( - (e * 2.0_default * gcc * mix_sl321))
    gpwsl1sn_3_c = ( - (e * 2.0_default * gcc * conjg (mix_sl311)))
    gpwsl2sn_3_c = ( - (e * 2.0_default * gcc * conjg (mix_sl321)))
    gwzsl1sn_3 = (gcc * gz * 2.0_default * sin2thw * mix_sl311)
    gwzsl2sn_3 = (gcc * gz * 2.0_default * sin2thw * mix_sl321)
    gwzsl1sn_3_c = (gcc * gz * 2.0_default * sin2thw * conjg (mix_sl311))
    gwzsl2sn_3_c = (gcc * gz * 2.0_default * sin2thw * conjg (mix_sl321))
    gpwpsu1sd1_1_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_11 *  &
      conjg (mix_su111) * mix_sd111)
    gpwpsu2sd2_1_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_11 *  &
      conjg (mix_su121) * mix_sd121)
    gpwpsu1sd2_1_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_11 *  &
      conjg (mix_su111) * mix_sd121)
    gpwpsu2sd1_1_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_11 *  &
      conjg (mix_su121) * mix_sd111)
    gpwpsu1sd1_1_1_c = conjg (gpwpsu1sd1_1_1) 
    gpwpsu2sd2_1_1_c = conjg (gpwpsu2sd2_1_1)
    gpwpsu1sd2_1_1_c = conjg (gpwpsu1sd2_1_1)
    gpwpsu2sd1_1_1_c = conjg (gpwpsu2sd1_1_1)
    gzwpsu1sd1_1_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_11 *  &
      conjg (mix_su111) * mix_sd111))
    gzwpsu2sd2_1_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_11 *  &
      conjg (mix_su121) * mix_sd121))
    gzwpsu1sd2_1_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_11 *  &
      conjg (mix_su111) * mix_sd121))
    gzwpsu2sd1_1_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_11 *  &
      conjg (mix_su121) * mix_sd111))
    gzwpsu1sd1_1_1_c = conjg (gzwpsu1sd1_1_1) 
    gzwpsu2sd2_1_1_c = conjg (gzwpsu2sd2_1_1)
    gzwpsu1sd2_1_1_c = conjg (gzwpsu1sd2_1_1)
    gzwpsu2sd1_1_1_c = conjg (gzwpsu2sd1_1_1)
    gpwpsu1sd1_1_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_12 *  &
      conjg (mix_su111) * mix_sd211)
    gpwpsu2sd2_1_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_12 *  &
      conjg (mix_su121) * mix_sd221)
    gpwpsu1sd2_1_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_12 *  &
      conjg (mix_su111) * mix_sd221)
    gpwpsu2sd1_1_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_12 *  &
      conjg (mix_su121) * mix_sd211)
    gpwpsu1sd1_1_2_c = conjg (gpwpsu1sd1_1_2) 
    gpwpsu2sd2_1_2_c = conjg (gpwpsu2sd2_1_2)
    gpwpsu1sd2_1_2_c = conjg (gpwpsu1sd2_1_2)
    gpwpsu2sd1_1_2_c = conjg (gpwpsu2sd1_1_2)
    gzwpsu1sd1_1_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_12 *  &
      conjg (mix_su111) * mix_sd211))
 end subroutine setup_parameters4
 subroutine setup_parameters5 ()
    gzwpsu2sd2_1_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_12 *  &
      conjg (mix_su121) * mix_sd221))
    gzwpsu1sd2_1_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_12 *  &
      conjg (mix_su111) * mix_sd221))
    gzwpsu2sd1_1_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_12 *  &
      conjg (mix_su121) * mix_sd211))
    gzwpsu1sd1_1_2_c = conjg (gzwpsu1sd1_1_2) 
    gzwpsu2sd2_1_2_c = conjg (gzwpsu2sd2_1_2)
    gzwpsu1sd2_1_2_c = conjg (gzwpsu1sd2_1_2)
    gzwpsu2sd1_1_2_c = conjg (gzwpsu2sd1_1_2)
    gpwpsu1sd1_1_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_13 *  &
      conjg (mix_su111) * mix_sd311)
    gpwpsu2sd2_1_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_13 *  &
      conjg (mix_su121) * mix_sd321)
    gpwpsu1sd2_1_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_13 *  &
      conjg (mix_su111) * mix_sd321)
    gpwpsu2sd1_1_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_13 *  &
      conjg (mix_su121) * mix_sd311)
    gpwpsu1sd1_1_3_c = conjg (gpwpsu1sd1_1_3) 
    gpwpsu2sd2_1_3_c = conjg (gpwpsu2sd2_1_3)
    gpwpsu1sd2_1_3_c = conjg (gpwpsu1sd2_1_3)
    gpwpsu2sd1_1_3_c = conjg (gpwpsu2sd1_1_3)
    gzwpsu1sd1_1_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_13 *  &
      conjg (mix_su111) * mix_sd311))
    gzwpsu2sd2_1_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_13 *  &
      conjg (mix_su121) * mix_sd321))
    gzwpsu1sd2_1_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_13 *  &
      conjg (mix_su111) * mix_sd321))
    gzwpsu2sd1_1_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_13 *  &
      conjg (mix_su121) * mix_sd311))
    gzwpsu1sd1_1_3_c = conjg (gzwpsu1sd1_1_3) 
    gzwpsu2sd2_1_3_c = conjg (gzwpsu2sd2_1_3)
    gzwpsu1sd2_1_3_c = conjg (gzwpsu1sd2_1_3)
    gzwpsu2sd1_1_3_c = conjg (gzwpsu2sd1_1_3)
    gpwpsu1sd1_2_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_21 *  &
      conjg (mix_su211) * mix_sd111)
    gpwpsu2sd2_2_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_21 *  &
      conjg (mix_su221) * mix_sd121)
    gpwpsu1sd2_2_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_21 *  &
      conjg (mix_su211) * mix_sd121)
    gpwpsu2sd1_2_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_21 *  &
      conjg (mix_su221) * mix_sd111)
    gpwpsu1sd1_2_1_c = conjg (gpwpsu1sd1_2_1) 
    gpwpsu2sd2_2_1_c = conjg (gpwpsu2sd2_2_1)
    gpwpsu1sd2_2_1_c = conjg (gpwpsu1sd2_2_1)
    gpwpsu2sd1_2_1_c = conjg (gpwpsu2sd1_2_1)
    gzwpsu1sd1_2_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_21 *  &
      conjg (mix_su211) * mix_sd111))
    gzwpsu2sd2_2_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_21 *  &
      conjg (mix_su221) * mix_sd121))
    gzwpsu1sd2_2_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_21 *  &
      conjg (mix_su211) * mix_sd121))
    gzwpsu2sd1_2_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_21 *  &
      conjg (mix_su221) * mix_sd111))
    gzwpsu1sd1_2_1_c = conjg (gzwpsu1sd1_2_1) 
    gzwpsu2sd2_2_1_c = conjg (gzwpsu2sd2_2_1)
    gzwpsu1sd2_2_1_c = conjg (gzwpsu1sd2_2_1)
    gzwpsu2sd1_2_1_c = conjg (gzwpsu2sd1_2_1)
    gpwpsu1sd1_2_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_22 *  &
      conjg (mix_su211) * mix_sd211)
    gpwpsu2sd2_2_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_22 *  &
      conjg (mix_su221) * mix_sd221)
    gpwpsu1sd2_2_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_22 *  &
      conjg (mix_su211) * mix_sd221)
    gpwpsu2sd1_2_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_22 *  &
      conjg (mix_su221) * mix_sd211)
    gpwpsu1sd1_2_2_c = conjg (gpwpsu1sd1_2_2)
    gpwpsu2sd2_2_2_c = conjg (gpwpsu2sd2_2_2)
    gpwpsu1sd2_2_2_c = conjg (gpwpsu1sd2_2_2)
    gpwpsu2sd1_2_2_c = conjg (gpwpsu2sd1_2_2)
    gzwpsu1sd1_2_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_22 *  &
      conjg (mix_su211) * mix_sd211))
    gzwpsu2sd2_2_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_22 *  &
      conjg (mix_su221) * mix_sd221))
    gzwpsu1sd2_2_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_22 *  &
      conjg (mix_su211) * mix_sd221))
    gzwpsu2sd1_2_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_22 *  &
      conjg (mix_su221) * mix_sd211))
    gzwpsu1sd1_2_2_c = conjg (gzwpsu1sd1_2_2) 
    gzwpsu2sd2_2_2_c = conjg (gzwpsu2sd2_2_2)
    gzwpsu1sd2_2_2_c = conjg (gzwpsu1sd2_2_2)
    gzwpsu2sd1_2_2_c = conjg (gzwpsu2sd1_2_2)
    gpwpsu1sd1_2_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_23 *  &
      conjg (mix_su211) * mix_sd311)
    gpwpsu2sd2_2_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_23 *  &
      conjg (mix_su221) * mix_sd321)
    gpwpsu1sd2_2_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_23 *  &
      conjg (mix_su211) * mix_sd321)
    gpwpsu2sd1_2_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_23 *  &
      conjg (mix_su221) * mix_sd311)
    gpwpsu1sd1_2_3_c = conjg (gpwpsu1sd1_2_3) 
    gpwpsu2sd2_2_3_c = conjg (gpwpsu2sd2_2_3)
    gpwpsu1sd2_2_3_c = conjg (gpwpsu1sd2_2_3)
    gpwpsu2sd1_2_3_c = conjg (gpwpsu2sd1_2_3)
    gzwpsu1sd1_2_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_23 *  &
      conjg (mix_su211) * mix_sd311))
    gzwpsu2sd2_2_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_23 *  &
      conjg (mix_su221) * mix_sd321))
    gzwpsu1sd2_2_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_23 *  &
      conjg (mix_su211) * mix_sd321))
    gzwpsu2sd1_2_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_23 *  &
      conjg (mix_su221) * mix_sd311))
    gzwpsu1sd1_2_3_c = conjg (gzwpsu1sd1_2_3) 
    gzwpsu2sd2_2_3_c = conjg (gzwpsu2sd2_2_3)
    gzwpsu1sd2_2_3_c = conjg (gzwpsu1sd2_2_3)
    gzwpsu2sd1_2_3_c = conjg (gzwpsu2sd1_2_3)
    gpwpsu1sd1_3_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_31 *  &
      conjg (mix_su311) * mix_sd111)
    gpwpsu2sd2_3_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_31 *  &
      conjg (mix_su321) * mix_sd121)
    gpwpsu1sd2_3_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_31 *  &
      conjg (mix_su311) * mix_sd121)
    gpwpsu2sd1_3_1 = (e * gcc * (2.0_default / 3.0_default) * vckm_31 *  &
      conjg (mix_su321) * mix_sd111)
    gpwpsu1sd1_3_1_c = conjg (gpwpsu1sd1_3_1) 
    gpwpsu2sd2_3_1_c = conjg (gpwpsu2sd2_3_1)
    gpwpsu1sd2_3_1_c = conjg (gpwpsu1sd2_3_1)
    gpwpsu2sd1_3_1_c = conjg (gpwpsu2sd1_3_1)
    gzwpsu1sd1_3_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_31 *  &
      conjg (mix_su311) * mix_sd111))
    gzwpsu2sd2_3_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_31 *  &
      conjg (mix_su321) * mix_sd121))
    gzwpsu1sd2_3_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_31 *  &
      conjg (mix_su311) * mix_sd121))
    gzwpsu2sd1_3_1 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_31 *  &
      conjg (mix_su321) * mix_sd111))
    gzwpsu1sd1_3_1_c = conjg (gzwpsu1sd1_3_1) 
    gzwpsu2sd2_3_1_c = conjg (gzwpsu2sd2_3_1)
    gzwpsu1sd2_3_1_c = conjg (gzwpsu1sd2_3_1)
    gzwpsu2sd1_3_1_c = conjg (gzwpsu2sd1_3_1)
    gpwpsu1sd1_3_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_32 *  &
      conjg (mix_su311) * mix_sd211)
    gpwpsu2sd2_3_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_32 *  &
      conjg (mix_su321) * mix_sd221)
    gpwpsu1sd2_3_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_32 *  &
      conjg (mix_su311) * mix_sd221)
    gpwpsu2sd1_3_2 = (e * gcc * (2.0_default / 3.0_default) * vckm_32 *  &
      conjg (mix_su321) * mix_sd211)
    gpwpsu1sd1_3_2_c = conjg (gpwpsu1sd1_3_2) 
    gpwpsu2sd2_3_2_c = conjg (gpwpsu2sd2_3_2)
    gpwpsu1sd2_3_2_c = conjg (gpwpsu1sd2_3_2)
    gpwpsu2sd1_3_2_c = conjg (gpwpsu2sd1_3_2)
    gzwpsu1sd1_3_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_32 *  &
      conjg (mix_su311) * mix_sd211))
    gzwpsu2sd2_3_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_32 *  &
      conjg (mix_su321) * mix_sd221))
    gzwpsu1sd2_3_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_32 *  &
      conjg (mix_su311) * mix_sd221))
    gzwpsu2sd1_3_2 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_32 *  &
      conjg (mix_su321) * mix_sd211))
    gzwpsu1sd1_3_2_c = conjg (gzwpsu1sd1_3_2) 
    gzwpsu2sd2_3_2_c = conjg (gzwpsu2sd2_3_2)
    gzwpsu1sd2_3_2_c = conjg (gzwpsu1sd2_3_2)
    gzwpsu2sd1_3_2_c = conjg (gzwpsu2sd1_3_2)
    gpwpsu1sd1_3_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_33 *  &
      conjg (mix_su311) * mix_sd311)
    gpwpsu2sd2_3_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_33 *  &
      conjg (mix_su321) * mix_sd321)
    gpwpsu1sd2_3_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_33 *  &
      conjg (mix_su311) * mix_sd321)
    gpwpsu2sd1_3_3 = (e * gcc * (2.0_default / 3.0_default) * vckm_33 *  &
      conjg (mix_su321) * mix_sd311)
    gpwpsu1sd1_3_3_c = conjg (gpwpsu1sd1_3_3) 
    gpwpsu2sd2_3_3_c = conjg (gpwpsu2sd2_3_3)
    gpwpsu1sd2_3_3_c = conjg (gpwpsu1sd2_3_3)
    gpwpsu2sd1_3_3_c = conjg (gpwpsu2sd1_3_3)
    gzwpsu1sd1_3_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_33 *  &
      conjg (mix_su311) * mix_sd311))
    gzwpsu2sd2_3_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_33 *  &
      conjg (mix_su321) * mix_sd321))
    gzwpsu1sd2_3_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_33 *  &
      conjg (mix_su311) * mix_sd321))
    gzwpsu2sd1_3_3 = ( - (gcc * gz *  &
      (2.0_default / 3.0_default) * sin2thw * vckm_33 *  &
      conjg (mix_su321) * mix_sd311))
    gzwpsu1sd1_3_3_c = conjg (gzwpsu1sd1_3_3)
    gzwpsu2sd2_3_3_c = conjg (gzwpsu2sd2_3_3)
    gzwpsu1sd2_3_3_c = conjg (gzwpsu1sd2_3_3)
    gzwpsu2sd1_3_3_c = conjg (gzwpsu2sd1_3_3)
    gglglsqsq = (gs**2)
 end subroutine setup_parameters5
 subroutine setup_parameters6 ()
    gglpsqsq = 2.0_default * e * gs / 3.0_default
    gglsu1su1_1 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su111 * conjg (mix_su111))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu2su2_1 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su121 * conjg (mix_su121))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu1su2_1 = (gz * gs * (1.0_default / 2.0_default) * mix_su111 *  &
      conjg (mix_su121))
    gglsu2su1_1 = (gz * gs * (1.0_default / 2.0_default) * mix_su121 *  &
      conjg (mix_su111))
    gglsd1sd1_1 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd111 * conjg (mix_sd111))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd2sd2_1 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd121 * conjg (mix_sd121))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd1sd2_1 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd111 * conjg (mix_sd121)))
    gglsd2sd1_1 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd121 * conjg (mix_sd111)))
    gglsu1su1_2 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su211 * conjg (mix_su211))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu2su2_2 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su221 * conjg (mix_su221))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu1su2_2 = (gz * gs * (1.0_default / 2.0_default) * mix_su211 *  &
      conjg (mix_su221))
    gglsu2su1_2 = (gz * gs * (1.0_default / 2.0_default) * mix_su221 *  &
      conjg (mix_su211))
    gglsd1sd1_2 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd211 * conjg (mix_sd211))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd2sd2_2 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd221 * conjg (mix_sd221))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd1sd2_2 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd211 * conjg (mix_sd221)))
    gglsd2sd1_2 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd221 * conjg (mix_sd211)))
    gglsu1su1_3 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su311 * conjg (mix_su311))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu2su2_3 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su321 * conjg (mix_su321))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu1su2_3 = (gz * gs * (1.0_default / 2.0_default) * mix_su311 *  &
      conjg (mix_su321))
    gglsu2su1_3 = (gz * gs * (1.0_default / 2.0_default) * mix_su321 *  &
      conjg (mix_su311))
    gglsd1sd1_3 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd311 * conjg (mix_sd311))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd2sd2_3 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd321 * conjg (mix_sd321))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd1sd2_3 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd311 * conjg (mix_sd321)))
    gglsd2sd1_3 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd321 * conjg (mix_sd311)))
    gglwsu1sd1_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su111) * mix_sd111)
    gglwsu2sd2_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su121) * mix_sd121)
    gglwsu1sd2_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su111) * mix_sd121)
    gglwsu2sd1_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su121) * mix_sd111)
    gglwsu1sd1_1_1_c = conjg (gglwsu1sd1_1_1) 
    gglwsu2sd2_1_1_c = conjg (gglwsu2sd2_1_1)
    gglwsu1sd2_1_1_c = conjg (gglwsu1sd2_1_1)
    gglwsu2sd1_1_1_c = conjg (gglwsu2sd1_1_1)
    gglwsu1sd1_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su111) * mix_sd211)
    gglwsu2sd2_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su121) * mix_sd221)
    gglwsu1sd2_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su111) * mix_sd221)
    gglwsu2sd1_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su121) * mix_sd211)
    gglwsu1sd1_1_2_c = conjg (gglwsu1sd1_1_2)
    gglwsu2sd2_1_2_c = conjg (gglwsu2sd2_1_2)
    gglwsu1sd2_1_2_c = conjg (gglwsu1sd2_1_2)
    gglwsu2sd1_1_2_c = conjg (gglwsu2sd1_1_2)
    gglwsu1sd1_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su111) * mix_sd311)
    gglwsu2sd2_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su121) * mix_sd321)
    gglwsu1sd2_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su111) * mix_sd321)
    gglwsu2sd1_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su121) * mix_sd311)
    gglwsu1sd1_1_3_c = conjg (gglwsu1sd1_1_3)  
    gglwsu2sd2_1_3_c = conjg (gglwsu2sd2_1_3)
    gglwsu1sd2_1_3_c = conjg (gglwsu1sd2_1_3)
    gglwsu2sd1_1_3_c = conjg (gglwsu2sd1_1_3)
    gglwsu1sd1_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su211) * mix_sd111)
    gglwsu2sd2_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su221) * mix_sd121)
    gglwsu1sd2_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su211) * mix_sd121)
    gglwsu2sd1_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su221) * mix_sd111)
    gglwsu1sd1_2_1_c = conjg (gglwsu1sd1_2_1) 
    gglwsu2sd2_2_1_c = conjg (gglwsu2sd2_2_1)
    gglwsu1sd2_2_1_c = conjg (gglwsu1sd2_2_1)
    gglwsu2sd1_2_1_c = conjg (gglwsu2sd1_2_1)
    gglwsu1sd1_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su211) * mix_sd211)
    gglwsu2sd2_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su221) * mix_sd221)
    gglwsu1sd2_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su211) * mix_sd221)
    gglwsu2sd1_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su221) * mix_sd211)
    gglwsu1sd1_2_2_c = conjg (gglwsu1sd1_2_2) 
    gglwsu2sd2_2_2_c = conjg (gglwsu2sd2_2_2)
    gglwsu1sd2_2_2_c = conjg (gglwsu1sd2_2_2)
    gglwsu2sd1_2_2_c = conjg (gglwsu2sd1_2_2)
    gglwsu1sd1_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su211) * mix_sd311)
    gglwsu2sd2_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su221) * mix_sd321)
    gglwsu1sd2_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su211) * mix_sd321)
    gglwsu2sd1_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su221) * mix_sd311)
    gglwsu1sd1_2_3_c = conjg (gglwsu1sd1_2_3)
    gglwsu2sd2_2_3_c = conjg (gglwsu2sd2_2_3)
    gglwsu1sd2_2_3_c = conjg (gglwsu1sd2_2_3)
    gglwsu2sd1_2_3_c = conjg (gglwsu2sd1_2_3)
    gglwsu1sd1_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su311) * mix_sd111)
    gglwsu2sd2_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su321) * mix_sd121)
    gglwsu1sd2_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su311) * mix_sd121)
    gglwsu2sd1_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su321) * mix_sd111)
    gglwsu1sd1_3_1_c = conjg (gglwsu1sd1_3_1)
    gglwsu2sd2_3_1_c = conjg (gglwsu2sd2_3_1)
    gglwsu1sd2_3_1_c = conjg (gglwsu1sd2_3_1)
    gglwsu2sd1_3_1_c = conjg (gglwsu2sd1_3_1)
    gglwsu1sd1_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su311) * mix_sd211)
    gglwsu2sd2_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su321) * mix_sd221)
    gglwsu1sd2_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su311) * mix_sd221)
 end subroutine setup_parameters6
 subroutine setup_parameters7 ()
    gglwsu2sd1_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su321) * mix_sd211)
    gglwsu1sd1_3_2_c = conjg (gglwsu1sd1_3_2)
    gglwsu2sd2_3_2_c = conjg (gglwsu2sd2_3_2)
    gglwsu1sd2_3_2_c = conjg (gglwsu1sd2_3_2)
    gglwsu2sd1_3_2_c = conjg (gglwsu2sd1_3_2)
    gglwsu1sd1_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su311) * mix_sd311)
    gglwsu2sd2_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su321) * mix_sd321)
    gglwsu1sd2_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su311) * mix_sd321)
    gglwsu2sd1_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su321) * mix_sd311)
    gglwsu1sd1_3_3_c = conjg (gglwsu1sd1_3_3) 
    gglwsu2sd2_3_3_c = conjg (gglwsu2sd2_3_3)
    gglwsu1sd2_3_3_c = conjg (gglwsu1sd2_3_3)
    gglwsu2sd1_3_3_c = conjg (gglwsu2sd1_3_3)
    axial0_11 = real ((mn_14 * conjg (mn_14)) - (mn_13 * conjg (mn_13))) &
      / 2.0_default
    snnh1_11 = 2.0_default * ( - real ((mn_12 - ( &
      (sinthw / costhw) * mn_11)) * ((sinal * mn_13) + (cosal * mn_14)))) 
    snnh2_11 = 2.0_default * real ((mn_12 - ((sinthw / costhw) * mn_11)) * &
      ((cosal * mn_13) - (sinal * mn_14))) 
    pnna_11 = 2.0_default * cmplx (0.0_default, real ((mn_12 - (mn_11 * (sinthw / costhw))) * ( &
      (mn_13 * sinbe) - (mn_14 * cosbe))),kind=default) 
    vector0_12 = cmplx (0.0_default, aimag ((mn_14 * conjg (mn_24)) -  &
      (mn_13 * conjg (mn_23))), kind=default) / 2.0_default
    axial0_12 = real ((mn_14 * conjg (mn_24)) - (mn_13 * conjg (mn_23))) &
      / 2.0_default
    snnh1_12 = ( - real (((mn_12 - ((sinthw / costhw) * mn_11)) * ((sinal &
       * mn_23) + (cosal * mn_24))) + ((mn_22 - ((sinthw / costhw) * mn_21)) &
      * ((sinal * mn_13) + (cosal * mn_14))))) 
    pnnh1_12 = ( - cmplx (0.0_default, aimag (((mn_12 - ( &
      (sinthw / costhw) * mn_11)) * ((sinal * mn_23) + (cosal * mn_24))) + ( &
      (mn_22 - ((sinthw / costhw) * mn_21)) * ((sinal * mn_13) +  &
      (cosal * mn_14)))), kind=default))
    snnh2_12 = real (((mn_12 - ((sinthw / costhw) * mn_11)) * ((cosal * & 
      mn_23) - (sinal * mn_24))) + ((mn_22 - ((sinthw / costhw) * mn_21)) * ( &
      (cosal * mn_13) - (sinal * mn_14)))) 
    pnnh2_12 = cmplx (0.0_default, aimag (((mn_12 - ((sinthw / costhw) & 
      * mn_11)) * ((cosal * mn_23) - (sinal * mn_24))) + ((mn_22 - ( &
      (sinthw / costhw) * mn_21)) * ((cosal * mn_13) -  &
      (sinal * mn_14)))), kind=default)
    snna_12 = - aimag (((mn_12 - (mn_11 *  &
      (sinthw / costhw))) * ((mn_23 * sinbe) - (mn_24 * cosbe))) + ((mn_22 -  &
      (mn_21 * (sinthw / costhw))) * ((mn_13 * sinbe) -  &
      (mn_14 * cosbe))))
    pnna_12 = cmplx (0.0_default, real (((mn_12 - (mn_11 * (sinthw / costhw))) * ((mn_23 * sinbe) &
      - (mn_24 * cosbe))) + ((mn_22 - (mn_21 * (sinthw / costhw))) * ( &
      (mn_13 * sinbe) - (mn_14 * cosbe)))),kind=default)
    vector0_13 = cmplx (0.0_default, aimag ((mn_14 * conjg (mn_34)) - (mn_13 * &
      conjg (mn_33))), kind=default) / 2.0_default 
    axial0_13 = real ((mn_14 * conjg (mn_34)) - (mn_13 * conjg (mn_33))) &
      / 2.0_default
    snnh1_13 = ( - real (((mn_12 - ((sinthw / costhw) * mn_11)) * ((sinal * &
       mn_33) + (cosal * mn_34))) + ((mn_32 - ((sinthw / costhw) * mn_31)) * ( &
      (sinal * mn_13) + (cosal * mn_14))))) 
    pnnh1_13 = ( - cmplx (0.0_default, aimag (((mn_12 - ( &
      (sinthw / costhw) * mn_11)) * ((sinal * mn_33) + (cosal * mn_34))) + ( &
      (mn_32 - ((sinthw / costhw) * mn_31)) * ((sinal * mn_13) +  &
      (cosal * mn_14)))), kind=default)) 
    snnh2_13 = real (((mn_12 - ((sinthw / costhw) * mn_11)) * ((cosal * & 
      mn_33) - (sinal * mn_34))) + ((mn_32 - ((sinthw / costhw) * mn_31)) * ( &
      (cosal * mn_13) - (sinal * mn_14)))) 
    pnnh2_13 = cmplx (0.0_default, aimag (((mn_12 - ((sinthw / costhw) * & 
      mn_11)) * ((cosal * mn_33) - (sinal * mn_34))) + ((mn_32 - ( &
      (sinthw / costhw) * mn_31)) * ((cosal * mn_13) -  &
      (sinal * mn_14)))), kind=default)
    snna_13 = - aimag (((mn_12 - (mn_11 *  &
      (sinthw / costhw))) * ((mn_33 * sinbe) - (mn_34 * cosbe))) + ((mn_32 -  &
      (mn_31 * (sinthw / costhw))) * ((mn_13 * sinbe) -  &
      (mn_14 * cosbe))))
    pnna_13 = cmplx (0.0_default, real (((mn_12 - (mn_11 * (sinthw / costhw))) * ((mn_33 * sinbe) &
      -  (mn_34 * cosbe))) + ((mn_32 - (mn_31 * (sinthw / costhw))) * ( &
      (mn_13 * sinbe) - (mn_14 * cosbe)))),kind=default)
    vector0_14 = cmplx (0.0_default, aimag ((mn_14 * conjg (mn_44)) - (mn_13 * & 
      conjg (mn_43))), kind=default) / 2.0_default
    axial0_14 = real ((mn_14 * conjg (mn_44)) - (mn_13 * conjg (mn_43))) &
      / 2.0_default
    snnh1_14 = ( - real (((mn_12 - ((sinthw / costhw) * mn_11)) * ((sinal * & 
       mn_43) + (cosal * mn_44))) + ((mn_42 - ((sinthw / costhw) * mn_41)) * ( &
      (sinal * mn_13) + (cosal * mn_14))))) 
    pnnh1_14 = ( - cmplx (0.0_default, aimag (((mn_12 - ( &
      (sinthw / costhw) * mn_11)) * ((sinal * mn_43) + (cosal * mn_44))) + ( &
      (mn_42 - ((sinthw / costhw) * mn_41)) * ((sinal * mn_13) +  &
      (cosal * mn_14)))), kind=default))
    snnh2_14 = real (((mn_12 - ((sinthw / costhw) * mn_11)) * ((cosal * & 
      mn_43) - (sinal * mn_44))) + ((mn_42 - ((sinthw / costhw) * mn_41)) * ( &
      (cosal * mn_13) - (sinal * mn_14)))) 
    pnnh2_14 = cmplx (0.0_default, aimag (((mn_12 - ((sinthw / costhw) * &
      mn_11)) * ((cosal * mn_43) - (sinal * mn_44))) + ((mn_42 - ( &
      (sinthw / costhw) * mn_41)) * ((cosal * mn_13) -  &
      (sinal * mn_14)))), kind=default) 
    snna_14 = - aimag (((mn_12 - (mn_11 *  &
      (sinthw / costhw))) * ((mn_43 * sinbe) - (mn_44 * cosbe))) + ((mn_42 -  &
      (mn_41 * (sinthw / costhw))) * ((mn_13 * sinbe) -  &
      (mn_14 * cosbe))))
    pnna_14 = cmplx (0.0_default, real (((mn_12 - (mn_11 * (sinthw / costhw))) * ((mn_43 * sinbe) &
      -  (mn_44 * cosbe))) + ((mn_42 - (mn_41 * (sinthw / costhw))) * ( &
      (mn_13 * sinbe) - (mn_14 * cosbe)))),kind=default)
    axial0_22 = real ((mn_24 * conjg (mn_24)) - (mn_23 * conjg (mn_23))) &
      / 2.0_default
    snnh1_22 = 2.0_default * ( - real ((mn_22 - ( &
      (sinthw / costhw) * mn_21)) * ((sinal * mn_23) + (cosal * mn_24))))
    snnh2_22 = 2.0_default * real ((mn_22 - ((sinthw / costhw) * mn_21)) & 
      * ((cosal * mn_23) - (sinal * mn_24))) 
    pnna_22 = 2.0_default * cmplx (0.0_default, real ((mn_22 - (mn_21 * (sinthw / costhw))) * ( &
      (mn_23 * sinbe) - (mn_24 * cosbe))),kind=default) 
    vector0_23 = cmplx (0.0_default, aimag ((mn_24 * conjg (mn_34)) - (mn_23 * & 
      conjg (mn_33))), kind=default) / 2.0_default
    axial0_23 = real ((mn_24 * conjg (mn_34)) - (mn_23 * conjg (mn_33))) &
      / 2.0_default
    snnh1_23 = ( - real (((mn_22 - ((sinthw / costhw) * mn_21)) * ((sinal * &
      mn_33) + (cosal * mn_34))) + ((mn_32 - ((sinthw / costhw) * mn_31)) * ( &
      (sinal * mn_23) + (cosal * mn_24))))) 
    pnnh1_23 = ( - cmplx (0.0_default, aimag (((mn_22 - ( &
      (sinthw / costhw) * mn_21)) * ((sinal * mn_33) + (cosal * mn_34))) + ( &
      (mn_32 - ((sinthw / costhw) * mn_31)) * ((sinal * mn_23) +  &
      (cosal * mn_24)))), kind=default))
    snnh2_23 = real (((mn_22 - ((sinthw / costhw) * mn_21)) * ((cosal * & 
      mn_33) - (sinal * mn_34))) + ((mn_32 - ((sinthw / costhw) * mn_31)) * ( &
      (cosal * mn_23) - (sinal * mn_24)))) 
    pnnh2_23 = cmplx (0.0_default, aimag (((mn_22 - ((sinthw / costhw) * &
      mn_21)) * ((cosal * mn_33) - (sinal * mn_34))) + ((mn_32 - ( &
      (sinthw / costhw) * mn_31)) * ((cosal * mn_23) -  &
      (sinal * mn_24)))), kind=default) 
    snna_23 = - aimag (((mn_22 - (mn_21 *  &
      (sinthw / costhw))) * ((mn_33 * sinbe) - (mn_34 * cosbe))) + ((mn_32 -  &
      (mn_31 * (sinthw / costhw))) * ((mn_23 * sinbe) -  &
      (mn_24 * cosbe))))
    pnna_23 = cmplx (0.0_default, real (((mn_22 - (mn_21 * (sinthw / costhw))) * ((mn_33 * sinbe) &
      - (mn_34 * cosbe))) + ((mn_32 - (mn_31 * (sinthw / costhw))) * ( &
      (mn_23 * sinbe) - (mn_24 * cosbe)))),kind=default) 
    vector0_24 = cmplx (0.0_default, aimag ((mn_24 * conjg (mn_44)) - (mn_23 * & 
      conjg (mn_43))), kind=default) / 2.0_default
    axial0_24 = real ((mn_24 * conjg (mn_44)) - (mn_23 * conjg (mn_43))) &
      / 2.0_default
    snnh1_24 = - real (((mn_22 - ((sinthw / costhw) * mn_21)) * ((sinal * & 
      mn_43) + (cosal * mn_44))) + ((mn_42 - ((sinthw / costhw) * mn_41)) * ( &
      (sinal * mn_23) + (cosal * mn_24))))
    pnnh1_24 = ( - cmplx (0.0_default, aimag (((mn_22 - ( &
      (sinthw / costhw) * mn_21)) * ((sinal * mn_43) + (cosal * mn_44))) + ( &
      (mn_42 - ((sinthw / costhw) * mn_41)) * ((sinal * mn_23) +  &
      (cosal * mn_24)))), kind=default))
    snnh2_24 = real (((mn_22 - ((sinthw / costhw) * mn_21)) * ((cosal * & 
      mn_43) - (sinal * mn_44))) + ((mn_42 - ((sinthw / costhw) * mn_41)) * ( &
      (cosal * mn_23) - (sinal * mn_24)))) 
    pnnh2_24 = cmplx (0.0_default, aimag (((mn_22 - ((sinthw / costhw) * &
      mn_21)) * ((cosal * mn_43) - (sinal * mn_44))) + ((mn_42 - ( &
      (sinthw / costhw) * mn_41)) * ((cosal * mn_23) -  &
      (sinal * mn_24)))), kind=default) 
    snna_24 = - aimag (((mn_22 - (mn_21 *  &
      (sinthw / costhw))) * ((mn_43 * sinbe) - (mn_44 * cosbe))) + ((mn_42 -  &
      (mn_41 * (sinthw / costhw))) * ((mn_23 * sinbe) -  &
      (mn_24 * cosbe))))
    pnna_24 = cmplx (0.0_default, real (((mn_22 - (mn_21 * (sinthw / costhw))) * ((mn_43 * sinbe) &
      - (mn_44 * cosbe))) + ((mn_42 - (mn_41 * (sinthw / costhw))) * ( &
      (mn_23 * sinbe) - (mn_24 * cosbe)))),kind=default)
    axial0_33 = real ((mn_34 * conjg (mn_34)) - (mn_33 * conjg (mn_33))) &
      / 2.0_default
    snnh1_33 = 2.0_default * ( - real ((mn_32 - ( &
      (sinthw / costhw) * mn_31)) * ((sinal * mn_33) + (cosal * mn_34)))) 
    snnh2_33 = 2.0_default * real ((mn_32 - ((sinthw / costhw) * mn_31)) & 
      * ((cosal * mn_33) - (sinal * mn_34))) 
    pnna_33 = 2.0_default * cmplx (0.0_default, real ((mn_32 - (mn_31 * (sinthw / costhw))) * ( &
      (mn_33 * sinbe) - (mn_34 * cosbe))),kind=default) 
 end subroutine setup_parameters7
 subroutine setup_parameters8 ()
    vector0_34 = cmplx (0.0_default, aimag ((mn_34 * conjg (mn_44)) - (mn_33 * & 
      conjg (mn_43))), kind=default) / 2.0_default
    axial0_34 = real ((mn_34 * conjg (mn_44)) - (mn_33 * conjg (mn_43))) &
      / 2.0_default
    snnh1_34 = ( - real (((mn_32 - ((sinthw / costhw) * mn_31)) * ((sinal * &
      mn_43) + (cosal * mn_44))) + ((mn_42 - ((sinthw / costhw) * mn_41)) * ( &
      (sinal * mn_33) + (cosal * mn_34))))) 
    pnnh1_34 = ( - cmplx (0.0_default, aimag (((mn_32 - ( &
      (sinthw / costhw) * mn_31)) * ((sinal * mn_43) + (cosal * mn_44))) + ( &
      (mn_42 - ((sinthw / costhw) * mn_41)) * ((sinal * mn_33) +  &
      (cosal * mn_34)))), kind=default))
    snnh2_34 = real (((mn_32 - ((sinthw / costhw) * mn_31)) * ((cosal * & 
      mn_43) - (sinal * mn_44))) + ((mn_42 - ((sinthw / costhw) * mn_41)) * ( &
      (cosal * mn_33) - (sinal * mn_34)))) 
    pnnh2_34 = cmplx (0.0_default, aimag (((mn_32 - ((sinthw / costhw) * & 
      mn_31)) * ((cosal * mn_43) - (sinal * mn_44))) + ((mn_42 - ( &
      (sinthw / costhw) * mn_41)) * ((cosal * mn_33) -  &
      (sinal * mn_34)))), kind=default) 
    snna_34 = - aimag (((mn_32 - (mn_31 *  &
      (sinthw / costhw))) * ((mn_43 * sinbe) - (mn_44 * cosbe))) + ((mn_42 -  &
      (mn_41 * (sinthw / costhw))) * ((mn_33 * sinbe) -  &
      (mn_34 * cosbe))))
    pnna_34 = cmplx (0.0_default, real (((mn_32 - (mn_31 * (sinthw / costhw))) * ((mn_43 * sinbe) &
      -  (mn_44 * cosbe))) + ((mn_42 - (mn_41 * (sinthw / costhw))) * ( &
      (mn_33 * sinbe) - (mn_34 * cosbe)))),kind=default)
    axial0_44 = real ((mn_44 * conjg (mn_44)) - (mn_43 * conjg (mn_43))) &
      / 2.0_default
    snnh1_44 = 2.0_default * ( - real ((mn_42 - ( &
      (sinthw / costhw) * mn_41)) * ((sinal * mn_43) + (cosal * mn_44))))
    snnh2_44 = 2.0_default * real ((mn_42 - ((sinthw / costhw) * mn_41)) & 
      * ((cosal * mn_43) - (sinal * mn_44))) 
    pnna_44 = 2.0_default * cmplx (0.0_default, real ((mn_42 - (mn_41 * (sinthw / costhw))) * ( &
      (mn_43 * sinbe) - (mn_44 * cosbe))),kind=default)
    vp_11 = ((((1.0_default -  &
      (2.0_default * sin2thw)) / 4.0_default) * ((mv_12 * conjg (mv_12)) &
       + (conjg (mu_12) * mu_12))) + (((costhw**2) / 2.0_default) * ( &
      (mv_11 * conjg (mv_11)) + (conjg (mu_11) * mu_11))))
    ap_11 = ((((1.0_default -  &
      (2.0_default * sin2thw)) / 4.0_default) * ((mv_12 *  &
      conjg (mv_12)) - (conjg (mu_12) * mu_12))) + (( &
      (costhw**2) / 2.0_default) * ((mv_11 * conjg (mv_11)) - ( &
      conjg (mu_11) * mu_11))))
    vp_12 = ((((1.0_default -  &
      (2.0_default * sin2thw)) / 4.0_default) * ((mv_12 * conjg (mv_22)) &
       + (conjg (mu_12) * mu_22))) + (((costhw**2) / 2.0_default) * ( &
      (mv_11 * conjg (mv_21)) + (conjg (mu_11) * mu_21))))
    ap_12 = ((((1.0_default -  &
      (2.0_default * sin2thw)) / 4.0_default) * ((mv_12 *  &
      conjg (mv_22)) - (conjg (mu_12) * mu_22))) + (( &
      (costhw**2) / 2.0_default) * ((mv_11 * conjg (mv_21)) - ( &
      conjg (mu_11) * mu_21))))
    vp_21 = conjg (vp_12)
    ap_21 = conjg (ap_12)
    vp_22 = ((((1.0_default -  &
      (2.0_default * sin2thw)) / 4.0_default) * ((mv_22 * conjg (mv_22)) &
       + (conjg (mu_22) * mu_22))) + (((costhw**2) / 2.0_default) * ( &
      (mv_21 * conjg (mv_21)) + (conjg (mu_21) * mu_21))))
    ap_22 = ((((1.0_default -  &
      (2.0_default * sin2thw)) / 4.0_default) * ((mv_22 *  &
      conjg (mv_22)) - (conjg (mu_22) * mu_22))) + (( &
      (costhw**2) / 2.0_default) * ((mv_21 * conjg (mv_21)) - ( &
      conjg (mu_21) * mu_21))))
    lcn_11 = ((conjg (mn_12) * mv_11 * sqrt (2.0_default)) - ( &
      conjg (mn_14) * mv_12))
    rcn_11 = ((mn_12 * conjg (mu_11) * sqrt (2.0_default)) + (mn_13 *  &
      conjg (mu_12)))
    lnch_11 = (cosbe * ((conjg (mn_14) * conjg (mv_11)) + ((conjg (mv_12) /  &
      sqrt (2.0_default)) * (conjg (mn_12) + ((sinthw / costhw) *  &
      conjg (mn_11))))))
    rnch_11 = (sinbe * ((mn_13 * mu_11) - ((mu_12 / sqrt (2.0_default)) *  &
      (mn_12 + ((sinthw / costhw) * mn_11)))))
    lcn_12 = ((conjg (mn_22) * mv_11 * sqrt (2.0_default)) - ( &
       conjg (mn_24) * mv_12))
    rcn_12 = ((mn_22 * conjg (mu_11) * sqrt (2.0_default)) + (mn_23 *  &
      conjg (mu_12)))
    lnch_21 = (cosbe * ((conjg (mn_24) * conjg (mv_11)) + ((conjg (mv_12) /  &
      sqrt (2.0_default)) * (conjg (mn_22) + ((sinthw / costhw) *  &
      conjg (mn_21))))))
    rnch_21 = (sinbe * ((mn_23 * mu_11) - ((mu_12 / sqrt (2.0_default)) *  &
      (mn_22 + ((sinthw / costhw) * mn_21)))))
    lcn_13 = ((conjg (mn_32) * mv_11 * sqrt (2.0_default)) - ( &
       conjg (mn_34) * mv_12))
    rcn_13 = ((mn_32 * conjg (mu_11) * sqrt (2.0_default)) + (mn_33 *  &
      conjg (mu_12)))
    lnch_31 = (cosbe * ((conjg (mn_34) * conjg (mv_11)) + ((conjg (mv_12) /  &
      sqrt (2.0_default)) * (conjg (mn_32) + ((sinthw / costhw) *  &
      conjg (mn_31))))))
    rnch_31 = (sinbe * ((mn_33 * mu_11) - ((mu_12 / sqrt (2.0_default)) *  &
      (mn_32 + ((sinthw / costhw) * mn_31)))))
    lcn_14 = ((conjg (mn_42) * mv_11 * sqrt (2.0_default)) - ( &
       conjg (mn_44) * mv_12))
    rcn_14 = ((mn_42 * conjg (mu_11) * sqrt (2.0_default)) + (mn_43 *  &
      conjg (mu_12)))
    lnch_41 = (cosbe * ((conjg (mn_44) * conjg (mv_11)) + ((conjg (mv_12) /  &
      sqrt (2.0_default)) * (conjg (mn_42) + ((sinthw / costhw) *  &
      conjg (mn_41))))))
    rnch_41 = (sinbe * ((mn_43 * mu_11) - ((mu_12 / sqrt (2.0_default)) *  &
      (mn_42 + ((sinthw / costhw) * mn_41)))))
    lcn_21 = ((conjg (mn_12) * mv_21 * sqrt (2.0_default)) - ( &
       conjg (mn_14) * mv_22))
    rcn_21 = ((mn_12 * conjg (mu_21) * sqrt (2.0_default)) + (mn_13 *  &
      conjg (mu_22)))
    lnch_12 = (cosbe * ((conjg (mn_14) * conjg (mv_21)) + ((conjg (mv_22) /  &
      sqrt (2.0_default)) * (conjg (mn_12) + ((sinthw / costhw) *  &
      conjg (mn_11))))))
 end subroutine setup_parameters8
subroutine setup_parameters9 ()
    rnch_12 = (sinbe * ((mn_13 * mu_21) - ((mu_22 / sqrt (2.0_default)) *  &
      (mn_12 + ((sinthw / costhw) * mn_11)))))
    lcn_22 = ((conjg (mn_22) * mv_21 * sqrt (2.0_default)) - ( &
       conjg (mn_24) * mv_22))
    rcn_22 = ((mn_22 * conjg (mu_21) * sqrt (2.0_default)) + (mn_23 *  &
      conjg (mu_22)))
    lnch_22 = (cosbe * ((conjg (mn_24) * conjg (mv_21)) + ((conjg (mv_22) /  &
      sqrt (2.0_default)) * (conjg (mn_22) + ((sinthw / costhw) *  &
      conjg (mn_21))))))
    rnch_22 = (sinbe * ((mn_23 * mu_21) - ((mu_22 / sqrt (2.0_default)) *  &
      (mn_22 + ((sinthw / costhw) * mn_21)))))
    lcn_23 = ((conjg (mn_32) * mv_21 * sqrt (2.0_default)) - ( &
       conjg (mn_34) * mv_22))
    rcn_23 = ((mn_32 * conjg (mu_21) * sqrt (2.0_default)) + (mn_33 *  &
      conjg (mu_22)))
    lnch_32 = (cosbe * ((conjg (mn_34) * conjg (mv_21)) + ((conjg (mv_22) /  &
      sqrt (2.0_default)) * (conjg (mn_32) + ((sinthw / costhw) *  &
      conjg (mn_31))))))
    rnch_32 = (sinbe * ((mn_33 * mu_21) - ((mu_22 / sqrt (2.0_default)) *  &
      (mn_32 + ((sinthw / costhw) * mn_31)))))
    lcn_24 = ((conjg (mn_42) * mv_21 * sqrt (2.0_default)) - ( &
       conjg (mn_44) * mv_22))
    rcn_24 = ((mn_42 * conjg (mu_21) * sqrt (2.0_default)) + (mn_43 *  &
      conjg (mu_22)))
    lnch_42 = (cosbe * ((conjg (mn_44) * conjg (mv_21)) + ((conjg (mv_22) /  &
      sqrt (2.0_default)) * (conjg (mn_42) + ((sinthw / costhw) *  &
      conjg (mn_41))))))
    rnch_42 = (sinbe * ((mn_43 * mu_21) - ((mu_22 / sqrt (2.0_default)) *  &
      (mn_42 + ((sinthw / costhw) * mn_41)))))
    lnc_11 = conjg (lcn_11)
    rnc_11 = conjg (rcn_11)
    lnc_12 = conjg (lcn_21)
    rnc_12 = conjg (rcn_21)
    lnc_21 = conjg (lcn_12)
    rnc_21 = conjg (rcn_12)
    lnc_22 = conjg (lcn_22)
    rnc_22 = conjg (rcn_22)
    lnc_31 = conjg (lcn_13)
    rnc_31 = conjg (rcn_13)
    lnc_32 = conjg (lcn_23)
    rnc_32 = conjg (rcn_23)
    lnc_41 = conjg (lcn_14)
    rnc_41 = conjg (rcn_14)
    lnc_42 = conjg (lcn_24)
    rnc_42 = conjg (rcn_24)
    gnzn_1_1 = (gz * axial0_11)
    gnzn_2_2 = (gz * axial0_22)
    gnzn_3_3 = (gz * axial0_33)
    gnzn_4_4 = (gz * axial0_44)
    dummy1 = ( - gs)
    !!! JR check 01.04.2005
    g_h1111susu = (gz * mass(23) * ((1.0_default / 2.0_default) -  &
      (sin2thw * q_up)) * sinapb)
    g_h1122susu = (gz * mass(23) *  q_up * sinapb * sin2thw)
    g_h1111sdsd = (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * q_down)) * sinapb)
    g_h1122sdsd = (gz * mass(23) * q_down * sinapb * sin2thw)
    g_h2111susu = ( - (gz * mass(23) * ((1.0_default / 2.0_default) -  &
      (sin2thw * q_up)) * cosapb))
    g_h2122susu = ( - gz * mass(23) * q_up * cosapb * sin2thw) 
    g_h2111sdsd = ( - (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * q_down)) * cosapb)) 
    g_h2122sdsd = ( - gz * mass(23) * q_down * cosapb * sin2thw)
    !!! g_h3112susu = - (imago * ((g * mass(2) * (( &
    !!!   conjg (au_1) * cosbe) + (mu * sinbe))) /  &
    !!!   (2.0_default * mass(24) * sinbe)))
    !!! g_h3121susu = conjg (g_h3112susu)
    !!! g_h3112sdsd = - (imago * ((g * mass(1) * (( &
    !!!   conjg (ad_1) * tanb) + mu)) / (2.0_default * mass(24))))
    !!! g_h3121sdsd = conjg (g_h3112sdsd)
    g_h1111snsn = (gz * mass(23) * (1.0_default / 2.0_default) * sinapb) 
    g_h1111slsl = (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * ( - 1.0_default))) * sinapb) 
    g_h1122slsl = (gz * mass(23) * ( - 1.0_default) * sinapb * sin2thw) 
    g_h2111snsn = ( - (gz * mass(23) * (1.0_default / 2.0_default)) * cosapb)
    g_h2111slsl = ( - (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * ( - 1.0_default))) * cosapb))
    g_h2122slsl = ( - (gz * mass(23) * ( - 1.0_default) * cosapb * sin2thw)) 
    !!! g_h3112slsl = - (imago * ((g * mass(11) * (( &
    !!!   conjg (al_1) * tanb) + mu)) / (2.0_default * mass(24))))
    !!! g_h3121slsl = conjg (g_h3112slsl)
    g_h4111slsn = ((- g / (sqrt (2.0_default) * mass(24))) * (mass(24)**2) * sin2be)
    !!! g_h4112slsn = (sqrt (2.0_default) * ((g * mass(11) * ((conjg ( &
    !!!   al_1) * sinbe) + (mu * cosbe))) /  &
    !!!   (2.0_default * mass(24) * cosbe)))
    g_h1211susu = g_h1111susu
    g_h1222susu = g_h1122susu
    g_h1211sdsd = g_h1111sdsd
    g_h1222sdsd = g_h1122sdsd 
    g_h2211susu = g_h2111susu
    g_h2222susu = g_h2122susu
    g_h2211sdsd = g_h2111sdsd
    g_h2222sdsd = g_h2122sdsd    
    !!! g_h1211susu = (gz * mass(23) * ((1.0_default / 2.0_default) -  &
    !!!   (sin2thw * q_up)) * sinapb)
    !!! g_h1222susu = (gz * mass(23) * q_up * sinapb * sin2thw)
    !!! g_h1211sdsd = (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
    !!!   (sin2thw * q_down)) * sinapb)
    !!! g_h1222sdsd = (gz * mass(23) * q_down * sinapb * sin2thw)
    !!! g_h2211susu = ( - (gz * mass(23) * ((1.0_default / 2.0_default) -  &
    !!!   (sin2thw * q_up)) * cosapb))
    !!! g_h2222susu = ( - gz * mass(23) * q_up * cosapb * sin2thw)
end subroutine setup_parameters9
 subroutine setup_parameters10 ()
    !!1 g_h2211sdsd = ( - (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
    !!1   (sin2thw * q_down)) * cosapb))
    !!1 g_h2222sdsd = ( - gz * mass(23) * q_down * cosapb * sin2thw) 
    !!! g_h3212susu = - (imago * ((g * mass(4) * (( &
    !!!   conjg (au_2) * cosbe) + (mu * sinbe))) /  &
    !!!   (2.0_default * mass(24) * sinbe)))
    !!! g_h3221susu = conjg (g_h3212susu)
    !!! g_h3212sdsd = - (imago * ((g * mass(3) * (( &
    !!!   conjg (ad_2) * sinbe) + (mu * cosbe))) /  &
    !!!   (2.0_default * mass(24) * cosbe)))
    !!! g_h3221sdsd = conjg (g_h3212sdsd)
    g_h1211snsn = g_h1111snsn
    g_h1211slsl = g_h1111slsl
    g_h1222slsl = g_h1122slsl
    g_h2211snsn = g_h2111snsn
    g_h2211slsl = g_h2111slsl
    g_h2222slsl = g_h2122slsl
    !!! g_h1211snsn = (gz * mass(23) * (1.0_default / 2.0_default))
    !!! g_h1211slsl = (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
    !!!   (sin2thw * ( - 1.0_default))) * sinapb)
    !!! g_h1222slsl = (gz * mass(23) * ( - 1.0_default) * sinapb * sin2thw)
    !!! g_h2211snsn = ( - (gz * mass(23) * (1.0_default / 2.0_default)))
    !!! g_h2211slsl = ( - (gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
    !!!   (sin2thw * ( - 1.0_default))) * cosapb))
    !!! g_h2222slsl = ( - (gz * mass(23) * ( - 1.0_default) * cosapb * sin2thw))
    !!! g_h3212slsl = - (imago * ((g * mass(13) * (( &
    !!!   conjg (al_2) * sinbe) + (mu * cosbe))) /  &
    !!!   (2.0_default * mass(24) * cosbe)))
    !!! g_h3221slsl = conjg (g_h3212slsl)
    g_h4211slsn = ((- g / (sqrt (2.0_default) * mass(24))) * (mass(24)**2) * sin2be)
    !!! g_h4212slsn = (sqrt (2.0_default) * ((g * mass(13) * ((conjg ( &
    !!!   al_2) * sinbe) + (mu * cosbe))) /  &
    !!!   (2.0_default * mass(24) * cosbe)))
    g_h1311susu = ((gz * mass(23) * ((1.0_default / 2.0_default) -  &
      (sin2thw * q_up)) * sinapb) - ((g *  &
      (mass(6)**2) * cosal) / (mass(24) * sinbe)))
    g_h1322susu = ((gz * mass(23) * q_up * sinapb * sin2thw) - ((g *  &
      (mass(6)**2) * cosal) / (mass(24) * sinbe)))
    g_h1312susu = - ((g * mass(6) * ((conjg (au_3) * cosal) + (  &
      mu * sinal))) / (2.0_default * mass(24) * sinbe))
    g_h1321susu = conjg (g_h1312susu)
    g_h1311sdsd = ((gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * q_down)) * sinapb) + ((g *  &
      (mass(5)**2) * sinal) / (mass(24) * cosbe)))
    g_h1322sdsd = ((gz * mass(23) * q_down * sinapb * sin2thw) + ((g *  &
      (mass(5)**2) * sinal) / (mass(24) * cosbe)))
    g_h1312sdsd = ((g * mass(5) * ((conjg (ad_3) * sinal) + ( &
       mu * cosal))) / (2.0_default * mass(24) * cosbe))
    g_h1321sdsd = conjg (g_h1312sdsd)
    g_h2311susu = ( - ((gz * mass(23) * ((1.0_default / 2.0_default) -  &
      (sin2thw * q_up)) * cosapb) + ((g *  &
      (mass(6)**2) * sinal) / (mass(24) * sinbe))))
    g_h2322susu = ( - ((gz * mass(23) * q_up * cosapb * sin2thw) + ((g *  &
      (mass(6)**2) * sinal) / (mass(24) * sinbe))))
    g_h2312susu = ((g * mass(6) * ((conjg (- au_3) * sinal) + ( &
       mu * cosal))) / (2.0_default * mass(24) * sinbe))
    g_h2321susu = conjg (g_h2312susu)
    g_h2311sdsd = ( - ((gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * q_down)) * cosapb) + ((g *  &
      (mass(5)**2) * cosal) / (mass(24) * cosbe))))
    g_h2322sdsd = ( - ((gz * mass(23) * q_down * cosapb * sin2thw) + ((g *  &
      (mass(5)**2) * cosal) / (mass(24) * cosbe))))
    g_h2312sdsd = ((g * mass(5) * ((conjg (- ad_3) * cosal) + ( &
       mu * sinal))) / (2.0_default * mass(24) * cosbe))
    g_h2321sdsd = conjg (g_h2312sdsd)
    g_h3312susu = - (imago * ((g * mass(6) * (( &
      conjg (au_3) * cosbe) + (mu * sinbe))) /  &
      (2.0_default * mass(24) * sinbe)))
    g_h3321susu = conjg (g_h3312susu)
    g_h3312sdsd = - (imago * ((g * mass(5) * (( &
      conjg (ad_3) * sinbe) + (mu * cosbe))) /  &
      (2.0_default * mass(24) * cosbe)))
    g_h3321sdsd = conjg (g_h3312sdsd)
    g_h1311snsn = (gz * mass(23) * (1.0_default / 2.0_default) * sinapb)
    g_h1311slsl = ((gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * ( - 1.0_default))) * sinapb) + ((g *  &
      (mass(15)**2) * sinal) / (mass(24) * cosbe)))
    g_h1322slsl = ((gz * mass(23) * ( - 1.0_default) * sinapb * sin2thw) + ( &
      (g * (mass(15)**2) * sinal) / (mass(24) * cosbe)))
    g_h1312slsl = ((g * mass(15) * ((conjg (al_3) * sinal) + ( &
       mu * cosal))) / (2.0_default * mass(24) * cosbe))
    g_h1321slsl = conjg (g_h1312slsl)
    g_h2311snsn = ( - (gz * mass(23) * (1.0_default / 2.0_default) * cosapb))
    g_h2311slsl = ( - ((gz * mass(23) * (( - (1.0_default / 2.0_default)) -  &
      (sin2thw * ( - 1.0_default))) * cosapb) + ((g *  &
      (mass(15)**2) * cosal) / (mass(24) * cosbe))))
    g_h2322slsl = ( - ((gz * mass(23) * ( - 1.0_default) * cosapb * sin2thw) + ( &
      (g * (mass(15)**2) * cosal) / (mass(24) * cosbe))))
    g_h2312slsl = ((g * mass(15) * ((conjg (- al_3) * cosal) + ( &
       mu * sinal))) / (2.0_default * mass(24) * cosbe))
    g_h2321slsl = conjg (g_h2312slsl)
    g_h3312slsl = - (imago * ((g * mass(15) * (( &
      conjg (al_3) * sinbe) + (mu * cosbe))) /  &
      (2.0_default * mass(24) * cosbe)))
    g_h3321slsl = conjg (g_h3312slsl)
    g_h4311slsn = ((g / (sqrt (2.0_default) * mass(24))) * (( &
      (mass(15)**2) * tanb) - ((mass(24)**2) * sin2be)))
    g_h4312slsn = (sqrt (2.0_default) * ((g * mass(15) * ((conjg ( &
      al_3) * sinbe) + (mu * cosbe))) /  &
      (2.0_default * mass(24) * cosbe)))
    g_h41_111susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_11 * ( - ( &
      (mass(24)**2) * sin2be)))
    g_h41_211susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_12 * ( - ( &
      (mass(24)**2) * sin2be)))
    g_h41_311susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_13 * (( - ( &
      (mass(24)**2) * sin2be)) + (((mass(5)**2) * tanb) + ((mass(2)**2) / tanb))))
    g_h41_322susd = ((sqrt (2.0_default) * g * mass(2) * mass(5) * vckm_13) /  &
      (mass(24) * sin2be))
    g_h41_312susd = (((g * mass(5)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_13 * (mu + ( &
      conjg (ad_3) * tanb)))
    g_h41_321susd = (((g * mass(2)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_13 * (conjg (mu) +  &
      (au_1 / tanb)))
    g_h42_111susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_21 * ( - ( &
      (mass(24)**2) * sin2be)))
    g_h42_211susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_22 * ( - ( &
      (mass(24)**2) * sin2be)))
    g_h42_311susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_23 * (( - ( &
      (mass(24)**2) * sin2be)) + (((mass(5)**2) * tanb) + ((mass(4)**2) / tanb))))
    g_h42_322susd = ((sqrt (2.0_default) * g * mass(4) * mass(5) * vckm_23) /  &
      (mass(24) * sin2be))
    g_h42_312susd = (((g * mass(5)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_23 * (mu + ( &
      conjg (ad_3) * tanb)))
    g_h42_321susd = (((g * mass(4)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_23 * (conjg (mu) +  &
      (au_2 / tanb)))
    g_h43_111susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_31 * (( - ( &
      (mass(24)**2) * sin2be)) + (((mass(1)**2) * tanb) + ((mass(6)**2) / tanb))))
    g_h43_122susd = ((sqrt (2.0_default) * g * mass(6) * mass(1) * vckm_31) /  &
      (mass(24) * sin2be))
    g_h43_112susd = (((g * mass(1)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_31 * (mu + ( &
      conjg (ad_1) * tanb)))
    g_h43_121susd = (((g * mass(6)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_31 * (conjg (mu) +  &
      (au_3 / tanb)))
    g_h43_211susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_32 * (( - ( &
      (mass(24)**2) * sin2be)) + (((mass(3)**2) * tanb) + ((mass(6)**2) / tanb))))
    g_h43_222susd = ((sqrt (2.0_default) * g * mass(6) * mass(3) * vckm_32) /  &
      (mass(24) * sin2be))
    g_h43_212susd = (((g * mass(3)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_32 * (mu + ( &
      conjg (ad_2) * tanb)))
    g_h43_221susd = (((g * mass(6)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_32 * (conjg (mu) +  &
      (au_3 / tanb)))
    g_h43_311susd = ((g / (sqrt (2.0_default) * mass(24))) * vckm_33 * (( - ( &
      (mass(24)**2) * sin2be)) + (((mass(5)**2) * tanb) + ((mass(6)**2) / tanb))))
    g_h43_322susd = ((sqrt (2.0_default) * g * mass(6) * mass(5) * vckm_33) /  &
      (mass(24) * sin2be))
    g_h43_312susd = (((g * mass(5)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_33 * (mu + ( &
      conjg (ad_3) * tanb)))
    g_h43_321susd = (((g * mass(6)) / ( &
      sqrt (2.0_default) * mass(24))) * vckm_33 * (conjg (mu) +  &
      (au_3 / tanb)))
 end subroutine setup_parameters10
 subroutine setup_parameters11 ()
    gh1sl1sl1_1 = g_h1111slsl
    gh1su1su1_1 = g_h1111susu
    gh1sd1sd1_1 = g_h1111sdsd
    gh2sl1sl1_1 = g_h2111slsl
    gh2su1su1_1 = g_h2111susu
    gh2sd1sd1_1 = g_h2111sdsd
    !!! gasl1sl1_1 = ((conjg (mix_sl111) * mix_sl112 * g_h3112slsl) + ( &
    !!!   conjg (mix_sl112) * mix_sl111 * g_h3121slsl))
    !!! gasu1su1_1 = ((conjg (mix_su111) * mix_su112 * g_h3112susu) + ( &
    !!!   conjg (mix_su112) * mix_su111 * g_h3121susu))
    !!! gasd1sd1_1 = ((conjg (mix_sd111) * mix_sd112 * g_h3112sdsd) + ( &
    !!!   conjg (mix_sd112) * mix_sd111 * g_h3121sdsd))
    !!! gasl1sl2_1 = ((conjg (mix_sl111) * mix_sl122 * g_h3112slsl) + ( &
    !!!   conjg (mix_sl112) * mix_sl121 * g_h3121slsl))
    !!! gasu1su2_1 = ((conjg (mix_su111) * mix_su122 * g_h3112susu) + ( &
    !!!   conjg (mix_su112) * mix_su121 * g_h3121susu))
    !!! gasd1sd2_1 = ((conjg (mix_sd111) * mix_sd122 * g_h3112sdsd) + ( &
    !!!   conjg (mix_sd112) * mix_sd121 * g_h3121sdsd))
    !!! gasl2sl1_1 = ((conjg (mix_sl121) * mix_sl112 * g_h3112slsl) + ( &
    !!!   conjg (mix_sl122) * mix_sl111 * g_h3121slsl))
    !!! gasu2su1_1 = ((conjg (mix_su121) * mix_su112 * g_h3112susu) + ( &
    !!!   conjg (mix_su122) * mix_su111 * g_h3121susu))
    !!! gasd2sd1_1 = ((conjg (mix_sd121) * mix_sd112 * g_h3112sdsd) + ( &
    !!!   conjg (mix_sd122) * mix_sd111 * g_h3121sdsd))
    gh1sl2sl2_1 = g_h1122slsl
    gh1su2su2_1 = g_h1122susu
    gh1sd2sd2_1 = g_h1122sdsd
    gh2sl2sl2_1 = g_h2122slsl
    gh2su2su2_1 = g_h2122susu
    gh2sd2sd2_1 = g_h2122sdsd
    !!! gasl2sl2_1 = ((conjg (mix_sl121) * mix_sl122 * g_h3112slsl) + ( &
    !!!   conjg (mix_sl122) * mix_sl121 * g_h3121slsl))
    !!! gasu2su2_1 = ((conjg (mix_su121) * mix_su122 * g_h3112susu) + ( &
    !!!   conjg (mix_su122) * mix_su121 * g_h3121susu))
    !!! gasd2sd2_1 = ((conjg (mix_sd121) * mix_sd122 * g_h3112sdsd) + ( &
    !!!   conjg (mix_sd122) * mix_sd121 * g_h3121sdsd))
    !!! ghsnsl1_1 = g_h4111slsn
    !!! ghsnsl1_1 = ((conjg (mix_sl111) * g_h4111slsn) + ( &
    !!!   conjg (mix_sl112) * g_h4112slsn))
    !!! ghsnsl1_1_c = conjg (ghsnsl1_1)
    ghsnsl1_1 = g_h4111slsn
    !!! ghsnsl2_1 = ((conjg (mix_sl121) * g_h4111slsn) + ( &
    !!!   conjg (mix_sl122) * g_h4112slsn))
    ghsnsl1_1_c = conjg (ghsnsl1_1)
    gh1sn1sn1_1 = g_h1111snsn
    gh2sn1sn1_1 = g_h2111snsn
    gh1sl1sl1_2 = g_h1211slsl
    gh1su1su1_2 = g_h1211susu
    gh1sd1sd1_2 = g_h1211sdsd
    gh2sl1sl1_2 = g_h2211slsl
    gh2su1su1_2 = g_h2211susu
    gh2sd1sd1_2 = g_h2211sdsd
    !!! gasl1sl1_2 = ((conjg (mix_sl211) * mix_sl212 * g_h3212slsl) + ( &
    !!!   conjg (mix_sl212) * mix_sl211 * g_h3221slsl))
    !!! gasu1su1_2 = ((conjg (mix_su211) * mix_su212 * g_h3212susu) + ( &
    !!!   conjg (mix_su212) * mix_su211 * g_h3221susu))
    !!! gasd1sd1_2 = ((conjg (mix_sd211) * mix_sd212 * g_h3212sdsd) + ( &
    !!!   conjg (mix_sd212) * mix_sd211 * g_h3221sdsd))
    !!! gasl1sl2_2 = ((conjg (mix_sl211) * mix_sl222 * g_h3212slsl) + ( &
    !!!   conjg (mix_sl212) * mix_sl221 * g_h3221slsl))
    !!! gasu1su2_2 = ((conjg (mix_su211) * mix_su222 * g_h3212susu) + ( &
    !!!   conjg (mix_su212) * mix_su221 * g_h3221susu))
    !!! gasd1sd2_2 = ((conjg (mix_sd211) * mix_sd222 * g_h3212sdsd) + ( &
    !!!   conjg (mix_sd212) * mix_sd221 * g_h3221sdsd))
    !!! gasl2sl1_2 = ((conjg (mix_sl221) * mix_sl212 * g_h3212slsl) + ( &
    !!!   conjg (mix_sl222) * mix_sl211 * g_h3221slsl))
    !!! gasu2su1_2 = ((conjg (mix_su221) * mix_su212 * g_h3212susu) + ( &
    !!!   conjg (mix_su222) * mix_su211 * g_h3221susu))
    !!! gasd2sd1_2 = ((conjg (mix_sd221) * mix_sd212 * g_h3212sdsd) + ( &
    !!!   conjg (mix_sd222) * mix_sd211 * g_h3221sdsd))
    gh1sl2sl2_2 = g_h1222slsl
    gh1su2su2_2 = g_h1222susu
    gh1sd2sd2_2 = g_h1222sdsd
    gh2sl2sl2_2 = g_h2222slsl
    gh2su2su2_2 = g_h2222susu
    gh2sd2sd2_2 = g_h2222sdsd
    !!! gasl2sl2_2 = ((conjg (mix_sl221) * mix_sl222 * g_h3212slsl) + ( &
    !!!   conjg (mix_sl222) * mix_sl221 * g_h3221slsl))
    !!! gasu2su2_2 = ((conjg (mix_su221) * mix_su222 * g_h3212susu) + ( &
    !!!   conjg (mix_su222) * mix_su221 * g_h3221susu))
    !!! gasd2sd2_2 = ((conjg (mix_sd221) * mix_sd222 * g_h3212sdsd) + ( &
    !!!   conjg (mix_sd222) * mix_sd221 * g_h3221sdsd))
    ghsnsl1_2 = g_h4211slsn
    !!! ghsnsl1_2 = ((conjg (mix_sl211) * g_h4211slsn) + ( &
    !!!   conjg (mix_sl212) * g_h4212slsn))
    ghsnsl1_2_c = conjg (ghsnsl1_2)
    !!! ghsnsl2_2 = g_h4211slsn
    !!! ghsnsl2_2 = ((conjg (mix_sl221) * g_h4211slsn) + ( &
    !!!   conjg (mix_sl222) * g_h4212slsn))
    !!! ghsnsl2_2_c = conjg (ghsnsl2_2)
    gh1sn1sn1_2 = g_h1211snsn
    gh2sn1sn1_2 = g_h2211snsn
    gh1sl1sl1_3 = ((conjg (mix_sl311) * mix_sl311 * g_h1311slsl) + ( &
      conjg (mix_sl312) * mix_sl312 * g_h1322slsl) + ( &
      conjg (mix_sl311) * mix_sl312 * g_h1312slsl) + ( &
      conjg (mix_sl312) * mix_sl311 * g_h1321slsl))
    gh1su1su1_3 = ((conjg (mix_su311) * mix_su311 * g_h1311susu) + ( &
      conjg (mix_su312) * mix_su312 * g_h1322susu) + ( &
      conjg (mix_su311) * mix_su312 * g_h1312susu) + ( &
      conjg (mix_su312) * mix_su311 * g_h1321susu))
    gh1sd1sd1_3 = ((conjg (mix_sd311) * mix_sd311 * g_h1311sdsd) + ( &
      conjg (mix_sd312) * mix_sd312 * g_h1322sdsd) + ( &
      conjg (mix_sd311) * mix_sd312 * g_h1312sdsd) + ( &
      conjg (mix_sd312) * mix_sd311 * g_h1321sdsd))
    gh2sl1sl1_3 = ((conjg (mix_sl311) * mix_sl311 * g_h2311slsl) + ( &
      conjg (mix_sl312) * mix_sl312 * g_h2322slsl) + ( &
      conjg (mix_sl311) * mix_sl312 * g_h2312slsl) + ( &
      conjg (mix_sl312) * mix_sl311 * g_h2321slsl))
    gh2su1su1_3 = ((conjg (mix_su311) * mix_su311 * g_h2311susu) + ( &
      conjg (mix_su312) * mix_su312 * g_h2322susu) + ( &
      conjg (mix_su311) * mix_su312 * g_h2312susu) + ( &
      conjg (mix_su312) * mix_su311 * g_h2321susu))
    gh2sd1sd1_3 = ((conjg (mix_sd311) * mix_sd311 * g_h2311sdsd) + ( &
      conjg (mix_sd312) * mix_sd312 * g_h2322sdsd) + ( &
      conjg (mix_sd311) * mix_sd312 * g_h2312sdsd) + ( &
      conjg (mix_sd312) * mix_sd311 * g_h2321sdsd))
    gasl1sl1_3 = ((conjg (mix_sl311) * mix_sl312 * g_h3312slsl) + ( &
      conjg (mix_sl312) * mix_sl311 * g_h3321slsl))
    gasu1su1_3 = ((conjg (mix_su311) * mix_su312 * g_h3312susu) + ( &
      conjg (mix_su312) * mix_su311 * g_h3321susu))
    gasd1sd1_3 = ((conjg (mix_sd311) * mix_sd312 * g_h3312sdsd) + ( &
      conjg (mix_sd312) * mix_sd311 * g_h3321sdsd))
    gh1sl1sl2_3 = ((conjg (mix_sl311) * mix_sl321 * g_h1311slsl) + ( &
      conjg (mix_sl312) * mix_sl322 * g_h1322slsl) + ( &
      conjg (mix_sl311) * mix_sl322 * g_h1312slsl) + ( &
      conjg (mix_sl312) * mix_sl321 * g_h1321slsl))
    gh1su1su2_3 = ((conjg (mix_su311) * mix_su321 * g_h1311susu) + ( &
      conjg (mix_su312) * mix_su322 * g_h1322susu) + ( &
      conjg (mix_su311) * mix_su322 * g_h1312susu) + ( &
      conjg (mix_su312) * mix_su321 * g_h1321susu))
    gh1sd1sd2_3 = ((conjg (mix_sd311) * mix_sd321 * g_h1311sdsd) + ( &
      conjg (mix_sd312) * mix_sd322 * g_h1322sdsd) + ( &
      conjg (mix_sd311) * mix_sd322 * g_h1312sdsd) + ( &
      conjg (mix_sd312) * mix_sd321 * g_h1321sdsd))
    gh2sl1sl2_3 = ((conjg (mix_sl311) * mix_sl321 * g_h2311slsl) + ( &
      conjg (mix_sl312) * mix_sl322 * g_h2322slsl) + ( &
      conjg (mix_sl311) * mix_sl322 * g_h2312slsl) + ( &
      conjg (mix_sl312) * mix_sl321 * g_h2321slsl))
    gh2su1su2_3 = ((conjg (mix_su311) * mix_su321 * g_h2311susu) + ( &
      conjg (mix_su312) * mix_su322 * g_h2322susu) + ( &
      conjg (mix_su311) * mix_su322 * g_h2312susu) + ( &
      conjg (mix_su312) * mix_su321 * g_h2321susu))
    gh2sd1sd2_3 = ((conjg (mix_sd311) * mix_sd321 * g_h2311sdsd) + ( &
      conjg (mix_sd312) * mix_sd322 * g_h2322sdsd) + ( &
      conjg (mix_sd311) * mix_sd322 * g_h2312sdsd) + ( &
      conjg (mix_sd312) * mix_sd321 * g_h2321sdsd))
    gasl1sl2_3 = ((conjg (mix_sl311) * mix_sl322 * g_h3312slsl) + ( &
      conjg (mix_sl312) * mix_sl321 * g_h3321slsl))
    gasu1su2_3 = ((conjg (mix_su311) * mix_su322 * g_h3312susu) + ( &
      conjg (mix_su312) * mix_su321 * g_h3321susu))
    gasd1sd2_3 = ((conjg (mix_sd311) * mix_sd322 * g_h3312sdsd) + ( &
      conjg (mix_sd312) * mix_sd321 * g_h3321sdsd))
    gh1sl2sl1_3 = ((conjg (mix_sl321) * mix_sl311 * g_h1311slsl) + ( &
      conjg (mix_sl322) * mix_sl312 * g_h1322slsl) + ( &
      conjg (mix_sl321) * mix_sl312 * g_h1312slsl) + ( &
      conjg (mix_sl322) * mix_sl311 * g_h1321slsl))
    gh1su2su1_3 = ((conjg (mix_su321) * mix_su311 * g_h1311susu) + ( &
      conjg (mix_su322) * mix_su312 * g_h1322susu) + ( &
      conjg (mix_su321) * mix_su312 * g_h1312susu) + ( &
      conjg (mix_su322) * mix_su311 * g_h1321susu))
    gh1sd2sd1_3 = ((conjg (mix_sd321) * mix_sd311 * g_h1311sdsd) + ( &
      conjg (mix_sd322) * mix_sd312 * g_h1322sdsd) + ( &
      conjg (mix_sd321) * mix_sd312 * g_h1312sdsd) + ( &
      conjg (mix_sd322) * mix_sd311 * g_h1321sdsd))
    gh2sl2sl1_3 = ((conjg (mix_sl321) * mix_sl311 * g_h2311slsl) + ( &
      conjg (mix_sl322) * mix_sl312 * g_h2322slsl) + ( &
      conjg (mix_sl321) * mix_sl312 * g_h2312slsl) + ( &
      conjg (mix_sl322) * mix_sl311 * g_h2321slsl))
    gh2su2su1_3 = ((conjg (mix_su321) * mix_su311 * g_h2311susu) + ( &
      conjg (mix_su322) * mix_su312 * g_h2322susu) + ( &
      conjg (mix_su321) * mix_su312 * g_h2312susu) + ( &
      conjg (mix_su322) * mix_su311 * g_h2321susu))
    gh2sd2sd1_3 = ((conjg (mix_sd321) * mix_sd311 * g_h2311sdsd) + ( &
      conjg (mix_sd322) * mix_sd312 * g_h2322sdsd) + ( &
      conjg (mix_sd321) * mix_sd312 * g_h2312sdsd) + ( &
      conjg (mix_sd322) * mix_sd311 * g_h2321sdsd))
    gasl2sl1_3 = ((conjg (mix_sl321) * mix_sl312 * g_h3312slsl) + ( &
      conjg (mix_sl322) * mix_sl311 * g_h3321slsl))
 end subroutine setup_parameters11
subroutine setup_parameters12 ()
    gasu2su1_3 = ((conjg (mix_su321) * mix_su312 * g_h3312susu) + ( &
      conjg (mix_su322) * mix_su311 * g_h3321susu))
    gasd2sd1_3 = ((conjg (mix_sd321) * mix_sd312 * g_h3312sdsd) + ( &
      conjg (mix_sd322) * mix_sd311 * g_h3321sdsd))
    gh1sl2sl2_3 = ((conjg (mix_sl321) * mix_sl321 * g_h1311slsl) + ( &
      conjg (mix_sl322) * mix_sl322 * g_h1322slsl) + ( &
      conjg (mix_sl321) * mix_sl322 * g_h1312slsl) + ( &
      conjg (mix_sl322) * mix_sl321 * g_h1321slsl))
    gh1su2su2_3 = ((conjg (mix_su321) * mix_su321 * g_h1311susu) + ( &
      conjg (mix_su322) * mix_su322 * g_h1322susu) + ( &
      conjg (mix_su321) * mix_su322 * g_h1312susu) + ( &
      conjg (mix_su322) * mix_su321 * g_h1321susu))
    gh1sd2sd2_3 = ((conjg (mix_sd321) * mix_sd321 * g_h1311sdsd) + ( &
      conjg (mix_sd322) * mix_sd322 * g_h1322sdsd) + ( &
      conjg (mix_sd321) * mix_sd322 * g_h1312sdsd) + ( &
      conjg (mix_sd322) * mix_sd321 * g_h1321sdsd))
    gh2sl2sl2_3 = ((conjg (mix_sl321) * mix_sl321 * g_h2311slsl) + ( &
      conjg (mix_sl322) * mix_sl322 * g_h2322slsl) + ( &
      conjg (mix_sl321) * mix_sl322 * g_h2312slsl) + ( &
      conjg (mix_sl322) * mix_sl321 * g_h2321slsl))
    gh2su2su2_3 = ((conjg (mix_su321) * mix_su321 * g_h2311susu) + ( &
      conjg (mix_su322) * mix_su322 * g_h2322susu) + ( &
      conjg (mix_su321) * mix_su322 * g_h2312susu) + ( &
      conjg (mix_su322) * mix_su321 * g_h2321susu))
    gh2sd2sd2_3 = ((conjg (mix_sd321) * mix_sd321 * g_h2311sdsd) + ( &
      conjg (mix_sd322) * mix_sd322 * g_h2322sdsd) + ( &
      conjg (mix_sd321) * mix_sd322 * g_h2312sdsd) + ( &
      conjg (mix_sd322) * mix_sd321 * g_h2321sdsd))
    gasl2sl2_3 = ((conjg (mix_sl321) * mix_sl322 * g_h3312slsl) + ( &
      conjg (mix_sl322) * mix_sl321 * g_h3321slsl))
    gasu2su2_3 = ((conjg (mix_su321) * mix_su322 * g_h3312susu) + ( &
      conjg (mix_su322) * mix_su321 * g_h3321susu))
    gasd2sd2_3 = ((conjg (mix_sd321) * mix_sd322 * g_h3312sdsd) + ( &
      conjg (mix_sd322) * mix_sd321 * g_h3321sdsd))
    ghsnsl1_3 = ((conjg (mix_sl311) * g_h4311slsn) + ( &
      conjg (mix_sl312) * g_h4312slsn))
    ghsnsl1_3_c = conjg (ghsnsl1_3)
    ghsnsl2_3 = ((conjg (mix_sl321) * g_h4311slsn) + ( &
      conjg (mix_sl322) * g_h4312slsn))
    ghsnsl2_3_c = conjg (ghsnsl2_3)
    gh1sn1sn1_3 = g_h1311snsn
    gh2sn1sn1_3 = g_h2311snsn
    ghsu1sd1_1_1 = g_h41_111susd
    ghsu1sd1_1_1_c = conjg (ghsu1sd1_1_1)
    ghsu1sd1_1_2 = g_h41_211susd
    ghsu1sd1_1_2_c = conjg (ghsu1sd1_1_2)
    ghsu1sd1_1_3 = ((conjg (mix_su111) * mix_sd311 * g_h41_311susd) + ( &
      conjg (mix_su112) * mix_sd312 * g_h41_322susd) + ( &
      conjg (mix_su111) * mix_sd312 * g_h41_312susd) + ( &
      conjg (mix_su112) * mix_sd311 * g_h41_321susd))
    ghsu1sd1_1_3_c = conjg (ghsu1sd1_1_3)
    ghsu1sd2_1_3 = ((conjg (mix_su111) * mix_sd321 * g_h41_311susd) + ( &
      conjg (mix_su112) * mix_sd322 * g_h41_322susd) + ( &
      conjg (mix_su111) * mix_sd322 * g_h41_312susd) + ( &
      conjg (mix_su112) * mix_sd321 * g_h41_321susd))
    ghsu1sd2_1_3_c = conjg (ghsu1sd2_1_3)
    ghsu2sd1_1_3 = ((conjg (mix_su121) * mix_sd311 * g_h41_311susd) + ( &
      conjg (mix_su122) * mix_sd312 * g_h41_322susd) + ( &
      conjg (mix_su121) * mix_sd312 * g_h41_312susd) + ( &
      conjg (mix_su122) * mix_sd311 * g_h41_321susd))
    ghsu2sd1_1_3_c = conjg (ghsu2sd1_1_3)
    ghsu2sd2_1_3 = ((conjg (mix_su121) * mix_sd321 * g_h41_311susd) + ( &
      conjg (mix_su122) * mix_sd322 * g_h41_322susd) + ( &
      conjg (mix_su121) * mix_sd322 * g_h41_312susd) + ( &
      conjg (mix_su122) * mix_sd321 * g_h41_321susd))
    ghsu2sd2_1_3_c = conjg (ghsu2sd2_1_3)
    ghsu1sd1_2_1 = g_h42_111susd
    ghsu1sd1_2_1_c = conjg (ghsu1sd1_2_1)
    ghsu1sd1_2_2 = g_h42_211susd
    ghsu1sd1_2_2_c = conjg (ghsu1sd1_2_2)
    ghsu1sd1_2_3 = ((conjg (mix_su211) * mix_sd311 * g_h42_311susd) + ( &
      conjg (mix_su212) * mix_sd312 * g_h42_322susd) + ( &
      conjg (mix_su211) * mix_sd312 * g_h42_312susd) + ( &
      conjg (mix_su212) * mix_sd311 * g_h42_321susd))
    ghsu1sd1_2_3_c = conjg (ghsu1sd1_2_3)
    ghsu1sd2_2_3 = ((conjg (mix_su211) * mix_sd321 * g_h42_311susd) + ( &
      conjg (mix_su212) * mix_sd322 * g_h42_322susd) + ( &
      conjg (mix_su211) * mix_sd322 * g_h42_312susd) + ( &
      conjg (mix_su212) * mix_sd321 * g_h42_321susd))
    ghsu1sd2_2_3_c = conjg (ghsu1sd2_2_3)
    ghsu2sd1_2_3 = ((conjg (mix_su221) * mix_sd311 * g_h42_311susd) + ( &
      conjg (mix_su222) * mix_sd312 * g_h42_322susd) + ( &
      conjg (mix_su221) * mix_sd312 * g_h42_312susd) + ( &
      conjg (mix_su222) * mix_sd311 * g_h42_321susd))
    ghsu2sd1_2_3_c = conjg (ghsu2sd1_2_3)
    ghsu2sd2_2_3 = ((conjg (mix_su221) * mix_sd321 * g_h42_311susd) + ( &
      conjg (mix_su222) * mix_sd322 * g_h42_322susd) + ( &
      conjg (mix_su221) * mix_sd322 * g_h42_312susd) + ( &
      conjg (mix_su222) * mix_sd321 * g_h42_321susd))
    ghsu2sd2_2_3_c = conjg (ghsu2sd2_2_3)
    ghsu1sd1_3_1 = ((conjg (mix_su311) * mix_sd111 * g_h43_111susd) + ( &
      conjg (mix_su312) * mix_sd112 * g_h43_122susd) + ( &
      conjg (mix_su311) * mix_sd112 * g_h43_112susd) + ( &
      conjg (mix_su312) * mix_sd111 * g_h43_121susd))
    ghsu1sd1_3_1_c = conjg (ghsu1sd1_3_1)
    ghsu1sd2_3_1 = ((conjg (mix_su311) * mix_sd121 * g_h43_111susd) + ( &
      conjg (mix_su312) * mix_sd122 * g_h43_122susd) + ( &
      conjg (mix_su311) * mix_sd122 * g_h43_112susd) + ( &
      conjg (mix_su312) * mix_sd121 * g_h43_121susd))
    ghsu1sd2_3_1_c = conjg (ghsu1sd2_3_1)
    ghsu2sd1_3_1 = ((conjg (mix_su321) * mix_sd111 * g_h43_111susd) + ( &
      conjg (mix_su322) * mix_sd112 * g_h43_122susd) + ( &
      conjg (mix_su321) * mix_sd112 * g_h43_112susd) + ( &
      conjg (mix_su322) * mix_sd111 * g_h43_121susd))
    ghsu2sd1_3_1_c = conjg (ghsu2sd1_3_1)
    ghsu2sd2_3_1 = ((conjg (mix_su321) * mix_sd121 * g_h43_111susd) + ( &
      conjg (mix_su322) * mix_sd122 * g_h43_122susd) + ( &
      conjg (mix_su321) * mix_sd122 * g_h43_112susd) + ( &
      conjg (mix_su322) * mix_sd121 * g_h43_121susd))
    ghsu2sd2_3_1_c = conjg (ghsu2sd2_3_1)
    ghsu1sd1_3_2 = ((conjg (mix_su311) * mix_sd211 * g_h43_211susd) + ( &
      conjg (mix_su312) * mix_sd212 * g_h43_222susd) + ( &
      conjg (mix_su311) * mix_sd212 * g_h43_212susd) + ( &
      conjg (mix_su312) * mix_sd211 * g_h43_221susd))
    ghsu1sd1_3_2_c = conjg (ghsu1sd1_3_2)
    ghsu1sd2_3_2 = ((conjg (mix_su311) * mix_sd221 * g_h43_211susd) + ( &
      conjg (mix_su312) * mix_sd222 * g_h43_222susd) + ( &
      conjg (mix_su311) * mix_sd222 * g_h43_212susd) + ( &
      conjg (mix_su312) * mix_sd221 * g_h43_221susd))
    ghsu1sd2_3_2_c = conjg (ghsu1sd2_3_2)
    ghsu2sd1_3_2 = ((conjg (mix_su321) * mix_sd211 * g_h43_211susd) + ( &
      conjg (mix_su322) * mix_sd212 * g_h43_222susd) + ( &
      conjg (mix_su321) * mix_sd212 * g_h43_212susd) + ( &
      conjg (mix_su322) * mix_sd211 * g_h43_221susd))
    ghsu2sd1_3_2_c = conjg (ghsu2sd1_3_2)
    ghsu2sd2_3_2 = ((conjg (mix_su321) * mix_sd221 * g_h43_211susd) + ( &
      conjg (mix_su322) * mix_sd222 * g_h43_222susd) + ( &
      conjg (mix_su321) * mix_sd222 * g_h43_212susd) + ( &
      conjg (mix_su322) * mix_sd221 * g_h43_221susd))
    ghsu2sd2_3_2_c = conjg (ghsu2sd2_3_2)
    ghsu1sd1_3_3 = ((conjg (mix_su311) * mix_sd311 * g_h43_311susd) + ( &
      conjg (mix_su312) * mix_sd312 * g_h43_322susd) + ( &
      conjg (mix_su311) * mix_sd312 * g_h43_312susd) + ( &
      conjg (mix_su312) * mix_sd311 * g_h43_321susd))
    ghsu1sd1_3_3_c = conjg (ghsu1sd1_3_3)
    ghsu1sd2_3_3 = ((conjg (mix_su311) * mix_sd321 * g_h43_311susd) + ( &
      conjg (mix_su312) * mix_sd322 * g_h43_322susd) + ( &
      conjg (mix_su311) * mix_sd322 * g_h43_312susd) + ( &
      conjg (mix_su312) * mix_sd321 * g_h43_321susd))
    ghsu1sd2_3_3_c = conjg (ghsu1sd2_3_3)
    ghsu2sd1_3_3 = ((conjg (mix_su321) * mix_sd311 * g_h43_311susd) + ( &
      conjg (mix_su322) * mix_sd312 * g_h43_322susd) + ( &
      conjg (mix_su321) * mix_sd312 * g_h43_312susd) + ( &
      conjg (mix_su322) * mix_sd311 * g_h43_321susd))
    ghsu2sd1_3_3_c = conjg (ghsu2sd1_3_3)
    ghsu2sd2_3_3 = ((conjg (mix_su321) * mix_sd321 * g_h43_311susd) + ( &
      conjg (mix_su322) * mix_sd322 * g_h43_322susd) + ( &
      conjg (mix_su321) * mix_sd322 * g_h43_312susd) + ( &
      conjg (mix_su322) * mix_sd321 * g_h43_321susd))
    ghsu2sd2_3_3_c = conjg (ghsu2sd2_3_3)
    g_yuk_ch1_sl1_1_c = ((( - g) / 2.0_default) * mu_11)
    g_yuk_ch1_sl1_1 = conjg (g_yuk_ch1_sl1_1_c)
    g_yuk_ch1_sl1_2_c = ((( - g) / 2.0_default) * mu_11)
    g_yuk_ch1_sl1_2 = conjg (g_yuk_ch1_sl1_2_c)
    g_yuk_ch1_sl1_3_c = ((((( - g) / 2.0_default) * mu_11) *  &
      conjg (mix_sl311)) + (((gcc * mass(15) * mu_12) / (mass(24) * cosbe)) * &
      conjg (mix_sl312)))
    g_yuk_ch1_sl1_3 = conjg (g_yuk_ch1_sl1_3_c)
    g_yuk_ch1_sl2_3_c = ((((( - g) / 2.0_default) * mu_11) *  &
      conjg (mix_sl321)) + (((gcc * mass(15) * mu_12) / (mass(24) * cosbe)) * &
      conjg (mix_sl322)))
    g_yuk_ch1_sl2_3 = conjg (g_yuk_ch1_sl2_3_c)
    g_yuk_ch1_sn1_1_c = ( - ((g / 2.0_default) * mv_11))
    g_yuk_ch1_sn1_1 = conjg (g_yuk_ch1_sn1_1_c)
    g_yuk_ch1_sn1_2_c = ( - ((g / 2.0_default) * mv_11))
    g_yuk_ch1_sn1_2 = conjg (g_yuk_ch1_sn1_2_c)
    g_yuk_ch2_sl1_1_c = ((( - g) / 2.0_default) * mu_21)
    g_yuk_ch2_sl1_1 = conjg (g_yuk_ch2_sl1_1_c)
    g_yuk_ch2_sl1_2_c = ((( - g) / 2.0_default) * mu_21)
    g_yuk_ch2_sl1_2 = conjg (g_yuk_ch2_sl1_2_c)
    g_yuk_ch2_sl1_3_c = ((((( - g) / 2.0_default) * mu_21) *  &
      conjg (mix_sl311)) + (((gcc * mass(15) * mu_22) / (mass(24) * cosbe)) * &
      conjg (mix_sl312)))
    g_yuk_ch2_sl1_3 = conjg (g_yuk_ch2_sl1_3_c)
    g_yuk_ch2_sl2_3_c = ((((( - g) / 2.0_default) * mu_21) *  &
      conjg (mix_sl321)) + (((gcc * mass(15) * mu_22) / (mass(24) * cosbe)) * &
      conjg (mix_sl322)))
    g_yuk_ch2_sl2_3 = conjg (g_yuk_ch2_sl2_3_c)
    g_yuk_ch2_sl2_3_c = conjg (g_yuk_ch2_sl2_3)
    g_yuk_ch2_sn1_1_c = ( - ((g / 2.0_default) * mv_21))
    g_yuk_ch2_sn1_1 = conjg (g_yuk_ch2_sn1_1_c)
    g_yuk_ch2_sn1_2_c = ( - ((g / 2.0_default) * mv_21))
    g_yuk_ch2_sn1_2 = conjg (g_yuk_ch2_sn1_2_c)
    g_yuk_ch1_sd1_1_1 = ( - ((g / 2.0_default) * conjg (mu_11) * vckm_11))
    g_yuk_ch1_sd1_1_1_c = conjg (g_yuk_ch1_sd1_1_1)
    g_yuk_ch1_su1_1_1 = ( - ((g / 2.0_default) * conjg (mv_11) * vckm_11))
    g_yuk_ch1_su1_1_1_c = conjg (g_yuk_ch1_su1_1_1)
    g_yuk_ch1_sd1_1_2 = ( - ((g / 2.0_default) * conjg (mu_11) * vckm_12))
    g_yuk_ch1_sd1_1_2_c = conjg (g_yuk_ch1_sd1_1_2)
    g_yuk_ch1_su1_1_2 = ( - ((g / 2.0_default) * conjg (mv_11) * vckm_12))
    g_yuk_ch1_su1_1_2_c = conjg (g_yuk_ch1_su1_1_2)
    g_yuk_ch1_sd1_2_1 = ( - ((g / 2.0_default) * conjg (mu_11) * vckm_21))
    g_yuk_ch1_sd1_2_1_c = conjg (g_yuk_ch1_sd1_2_1)
    g_yuk_ch1_su1_2_1 = ( - ((g / 2.0_default) * conjg (mv_11) * vckm_21))
    g_yuk_ch1_su1_2_1_c = conjg (g_yuk_ch1_su1_2_1)
    g_yuk_ch1_sd1_2_2 = ( - ((g / 2.0_default) * conjg (mu_11) * vckm_22))
    g_yuk_ch1_sd1_2_2_c = conjg (g_yuk_ch1_sd1_2_2)
    g_yuk_ch1_su1_2_2 = ( - ((g / 2.0_default) * conjg (mv_11) * vckm_22))
    g_yuk_ch1_su1_2_2_c = conjg (g_yuk_ch1_su1_2_2)
    g_yuk_ch2_sd1_1_1 = ( - ((g / 2.0_default) * conjg (mu_21) * vckm_11))
    g_yuk_ch2_sd1_1_1_c = conjg (g_yuk_ch2_sd1_1_1)
    g_yuk_ch2_su1_1_1 = ( - ((g / 2.0_default) * conjg (mv_21) * vckm_11))
    g_yuk_ch2_su1_1_1_c = conjg (g_yuk_ch2_su1_1_1)
    g_yuk_ch2_sd1_1_2 = ( - ((g / 2.0_default) * conjg (mu_21) * vckm_12))
    g_yuk_ch2_sd1_1_2_c = conjg (g_yuk_ch2_sd1_1_2)
    g_yuk_ch2_su1_1_2 = ( - ((g / 2.0_default) * conjg (mv_21) * vckm_12))
    g_yuk_ch2_su1_1_2_c = conjg (g_yuk_ch2_su1_1_2)
    g_yuk_ch2_sd1_2_1 = ( - ((g / 2.0_default) * conjg (mu_21) * vckm_21))
    g_yuk_ch2_sd1_2_1_c = conjg (g_yuk_ch2_sd1_2_1)
    g_yuk_ch2_su1_2_1 = ( - ((g / 2.0_default) * conjg (mv_21) * vckm_21))
    g_yuk_ch2_su1_2_1_c = conjg (g_yuk_ch2_su1_2_1)
    g_yuk_ch2_sd1_2_2 = ( - ((g / 2.0_default) * conjg (mu_21) * vckm_22))
    g_yuk_ch2_sd1_2_2_c = conjg (g_yuk_ch2_sd1_2_2)
    g_yuk_ch2_su1_2_2 = ( - ((g / 2.0_default) * conjg (mv_21) * vckm_22))
    g_yuk_ch2_su1_2_2_c = conjg (g_yuk_ch2_su1_2_2)
   g_yuk_n1_sn1_1 = (gcc * ((mn_11 * (sinthw / costhw)) - mn_12))
   g_yuk_n1_sn1_1_c = conjg (g_yuk_n1_sn1_1)
   g_yuk_n1_sn1_2 = (gcc * ((mn_11 * (sinthw / costhw)) - mn_12))
   g_yuk_n1_sn1_2_c = conjg (g_yuk_n1_sn1_2)
   g_yuk_n1_sn1_3 = (gcc * ((mn_11 * (sinthw / costhw)) - mn_12))
   g_yuk_n1_sn1_3_c = conjg (g_yuk_n1_sn1_3)
    g_yuk_n2_sn1_1 = (gcc * ((mn_21 * (sinthw / costhw)) - mn_22))
    g_yuk_n2_sn1_1_c = conjg (g_yuk_n2_sn1_1)
    g_yuk_n2_sn1_2 = (gcc * ((mn_21 * (sinthw / costhw)) - mn_22))
    g_yuk_n2_sn1_2_c = conjg (g_yuk_n2_sn1_2)
    g_yuk_n2_sn1_3 = (gcc * ((mn_21 * (sinthw / costhw)) - mn_22))
    g_yuk_n2_sn1_3_c = conjg (g_yuk_n2_sn1_3)
    g_yuk_n3_sn1_1 = (gcc * ((mn_31 * (sinthw / costhw)) - mn_32))
    g_yuk_n3_sn1_1_c = conjg (g_yuk_n3_sn1_1)
    g_yuk_n3_sn1_2 = (gcc * ((mn_31 * (sinthw / costhw)) - mn_32))
    g_yuk_n3_sn1_2_c = conjg (g_yuk_n3_sn1_2)
    g_yuk_n3_sn1_3 = (gcc * ((mn_31 * (sinthw / costhw)) - mn_32))
    g_yuk_n3_sn1_3_c = conjg (g_yuk_n3_sn1_3)
    g_yuk_n4_sn1_1 = (gcc * ((mn_41 * (sinthw / costhw)) - mn_42))
    g_yuk_n4_sn1_1_c = conjg (g_yuk_n4_sn1_1)
    g_yuk_n4_sn1_2 = (gcc * ((mn_41 * (sinthw / costhw)) - mn_42))
    g_yuk_n4_sn1_2_c = conjg (g_yuk_n4_sn1_2)
    g_yuk_n4_sn1_3 = (gcc * ((mn_41 * (sinthw / costhw)) - mn_42))
    g_yuk_n4_sn1_3_c = conjg (g_yuk_n4_sn1_3)
    g_yuk_n1_sl1_1 = (gcc * (mn_12 + ((sinthw * mn_11) / costhw)))
    g_yuk_n1_sl1_1_c = conjg (g_yuk_n1_sl1_1)
    g_yuk_n1_sl2_1 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_11)) / costhw)
    g_yuk_n1_sl2_1_c = conjg (g_yuk_n1_sl2_1)
    g_yuk_n1_su1_1 = (( - gcc) * (mn_12 + ((sinthw * mn_11) /  &
      (3.0_default * costhw))))
    g_yuk_n1_su1_1_c = conjg (g_yuk_n1_su1_1)
    g_yuk_n1_su2_1 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_11)) / costhw)
    g_yuk_n1_su2_1_c = conjg (g_yuk_n1_su2_1)
    g_yuk_n1_sd1_1 = (gcc * (mn_12 - ((sinthw * mn_11) /  &
      (costhw * 3.0_default))))
    g_yuk_n1_sd1_1_c = conjg (g_yuk_n1_sd1_1)
    g_yuk_n1_sd2_1 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_11)) / costhw)
    g_yuk_n1_sd2_1_c = conjg (g_yuk_n1_sd2_1)
    g_yuk_n2_sl1_1 = (gcc * (mn_22 + ((sinthw * mn_21) / costhw)))
    g_yuk_n2_sl1_1_c = conjg (g_yuk_n2_sl1_1)
    g_yuk_n2_sl2_1 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_21)) / costhw)
    g_yuk_n2_sl2_1_c = conjg (g_yuk_n2_sl2_1)
    g_yuk_n2_su1_1 = (( - gcc) * (mn_22 + ((sinthw * mn_21) /  &
      (3.0_default * costhw))))
    g_yuk_n2_su1_1_c = conjg (g_yuk_n2_su1_1)
    g_yuk_n2_su2_1 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_21)) / costhw)
    g_yuk_n2_su2_1_c = conjg (g_yuk_n2_su2_1)
    g_yuk_n2_sd1_1 = (gcc * (mn_22 - ((sinthw * mn_21) /  &
      (costhw * 3.0_default))))
    g_yuk_n2_sd1_1_c = conjg (g_yuk_n2_sd1_1)
    g_yuk_n2_sd2_1 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_21)) / costhw)
    g_yuk_n2_sd2_1_c = conjg (g_yuk_n2_sd2_1)
    g_yuk_n3_sl1_1 = (gcc * (mn_32 + ((sinthw * mn_31) / costhw)))
    g_yuk_n3_sl1_1_c = conjg (g_yuk_n3_sl1_1)
    g_yuk_n3_sl2_1 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_31)) / costhw)
    g_yuk_n3_sl2_1_c = conjg (g_yuk_n3_sl2_1)
    g_yuk_n3_su1_1 = (( - gcc) * (mn_32 + ((sinthw * mn_31) /  &
      (3.0_default * costhw))))
    g_yuk_n3_su1_1_c = conjg (g_yuk_n3_su1_1)
    g_yuk_n3_su2_1 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_31)) / costhw)
    g_yuk_n3_su2_1_c = conjg (g_yuk_n3_su2_1)
    g_yuk_n3_sd1_1 = (gcc * (mn_32 - ((sinthw * mn_31) /  &
      (costhw * 3.0_default))))
    g_yuk_n3_sd1_1_c = conjg (g_yuk_n3_sd1_1)
    g_yuk_n3_sd2_1 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_31)) / costhw)
    g_yuk_n3_sd2_1_c = conjg (g_yuk_n3_sd2_1)
    g_yuk_n4_sl1_1 = (gcc * (mn_42 + ((sinthw * mn_41) / costhw)))
    g_yuk_n4_sl1_1_c = conjg (g_yuk_n4_sl1_1)
    g_yuk_n4_sl2_1 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_41)) / costhw)
    g_yuk_n4_sl2_1_c = conjg (g_yuk_n4_sl2_1)
    g_yuk_n4_su1_1 = (( - gcc) * (mn_42 + ((sinthw * mn_41) /  &
      (3.0_default * costhw))))
    g_yuk_n4_su1_1_c = conjg (g_yuk_n4_su1_1)
    g_yuk_n4_su2_1 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_41)) / costhw)
    g_yuk_n4_su2_1_c = conjg (g_yuk_n4_su2_1)
    g_yuk_n4_sd1_1 = (gcc * (mn_42 - ((sinthw * mn_41) /  &
      (costhw * 3.0_default))))
    g_yuk_n4_sd1_1_c = conjg (g_yuk_n4_sd1_1)
    g_yuk_n4_sd2_1 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_41)) / costhw)
    g_yuk_n4_sd2_1_c = conjg (g_yuk_n4_sd2_1)
    g_yuk_n1_sl1_2 = (gcc * (mn_12 + ((sinthw * mn_11) / costhw)))
    g_yuk_n1_sl1_2_c = conjg (g_yuk_n1_sl1_2)
    g_yuk_n1_sl2_2 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_11)) / costhw)
    g_yuk_n1_sl2_2_c = conjg (g_yuk_n1_sl2_2)
    g_yuk_n1_su1_2 = (( - gcc) * (mn_12 + ((sinthw * mn_11) /  &
      (3.0_default * costhw))))
    g_yuk_n1_su1_2_c = conjg (g_yuk_n1_su1_2)
    g_yuk_n1_su2_2 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_11)) / costhw)
    g_yuk_n1_su2_2_c = conjg (g_yuk_n1_su2_2)
    g_yuk_n1_sd1_2 = (gcc * (mn_12 - ((sinthw * mn_11) /  &
      (costhw * 3.0_default))))
    g_yuk_n1_sd1_2_c = conjg (g_yuk_n1_sd1_2)
    g_yuk_n1_sd2_2 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_11)) / costhw)
    g_yuk_n1_sd2_2_c = conjg (g_yuk_n1_sd2_2)
    g_yuk_n2_sl1_2 = (gcc * (mn_22 + ((sinthw * mn_21) / costhw)))
    g_yuk_n2_sl1_2_c = conjg (g_yuk_n2_sl1_2)
    g_yuk_n2_sl2_2 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_21)) / costhw)
    g_yuk_n2_sl2_2_c = conjg (g_yuk_n2_sl2_2)
    g_yuk_n2_su1_2 = (( - gcc) * (mn_22 + ((sinthw * mn_21) /  &
      (3.0_default * costhw))))
    g_yuk_n2_su1_2_c = conjg (g_yuk_n2_su1_2)
    g_yuk_n2_su2_2 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_21)) / costhw)
    g_yuk_n2_su2_2_c = conjg (g_yuk_n2_su2_2)
    g_yuk_n2_sd1_2 = (gcc * (mn_22 - ((sinthw * mn_21) /  &
      (costhw * 3.0_default))))
    g_yuk_n2_sd1_2_c = conjg (g_yuk_n2_sd1_2)
    g_yuk_n2_sd2_2 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_21)) / costhw)
    g_yuk_n2_sd2_2_c = conjg (g_yuk_n2_sd2_2)
    g_yuk_n3_sl1_2 = (gcc * (mn_32 + ((sinthw * mn_31) / costhw)))
    g_yuk_n3_sl1_2_c = conjg (g_yuk_n3_sl1_2)
    g_yuk_n3_sl2_2 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_31)) / costhw)
    g_yuk_n3_sl2_2_c = conjg (g_yuk_n3_sl2_2)
    g_yuk_n3_su1_2 = (( - gcc) * (mn_32 + ((sinthw * mn_31) /  &
      (3.0_default * costhw))))
    g_yuk_n3_su1_2_c = conjg (g_yuk_n3_su1_2)
    g_yuk_n3_su2_2 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_31)) / costhw)
    g_yuk_n3_su2_2_c = conjg (g_yuk_n3_su2_2)
    g_yuk_n3_sd1_2 = (gcc * (mn_32 - ((sinthw * mn_31) /  &
      (costhw * 3.0_default))))
    g_yuk_n3_sd1_2_c = conjg (g_yuk_n3_sd1_2)
    g_yuk_n3_sd2_2 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_31)) / costhw)
    g_yuk_n3_sd2_2_c = conjg (g_yuk_n3_sd2_2)
    g_yuk_n4_sl1_2 = (gcc * (mn_42 + ((sinthw * mn_41) / costhw)))
    g_yuk_n4_sl1_2_c = conjg (g_yuk_n4_sl1_2)
    g_yuk_n4_sl2_2 = ((gcc * 2.0_default * q_lep * sinthw *  &
      conjg (mn_41)) / costhw)
    g_yuk_n4_sl2_2_c = conjg (g_yuk_n4_sl2_2)
    g_yuk_n4_su1_2 = (( - gcc) * (mn_42 + ((sinthw * mn_41) /  &
      (3.0_default * costhw))))
    g_yuk_n4_su1_2_c = conjg (g_yuk_n4_su1_2)
    g_yuk_n4_su2_2 = ((gcc * 2.0_default * q_up * sinthw *  &
      conjg (mn_41)) / costhw)
    g_yuk_n4_su2_2_c = conjg (g_yuk_n4_su2_2)
    g_yuk_n4_sd1_2 = (gcc * (mn_42 - ((sinthw * mn_41) /  &
      (costhw * 3.0_default))))
    g_yuk_n4_sd1_2_c = conjg (g_yuk_n4_sd1_2)
    g_yuk_n4_sd2_2 = ((gcc * 2.0_default * q_down * sinthw *  &
      conjg (mn_41)) / costhw)
    g_yuk_n4_sd2_2_c = conjg (g_yuk_n4_sd2_2)
end subroutine setup_parameters12
subroutine setup_parameters13 ()
    gncneu(1) = ((gz / 2.0_default) * ( &
      (2.0_default * 0.0_default * sin2thw) -  &
      (1.0_default / 2.0_default)))
    gncneu(2) = ((( - gz) / 2.0_default) *  &
      (1.0_default / 2.0_default))
    gnclep(1) = ((gz / 2.0_default) * ( &
      (2.0_default * (-1.0_default) * sin2thw) - ( -  &
      (1.0_default / 2.0_default))))
    gnclep(2) = ((( - gz) / 2.0_default) * ( -  &
      (1.0_default / 2.0_default)))
    gncup(1) = ((gz / 2.0_default) * ((2.0_default *  &
      (2.0_default / 3.0_default) * sin2thw) -  &
      (1.0_default / 2.0_default)))
    gncup(2) = ((( - gz) / 2.0_default) * (1.0_default / 2.0_default))
    gncdwn(1) = ((gz / 2.0_default) * ((2.0_default *  &
      ((-1.0_default) / 3.0_default) * sin2thw) - ( -  &
      (1.0_default / 2.0_default))))
    gncdwn(2) = ((( - gz) / 2.0_default) * ( -  &
      (1.0_default / 2.0_default)))
    g_yuk1_1_3(1) = ((gcc / mass(24)) * vckm_13 * (mass(2) / tanb))
    g_yuk1_1_3(2) = ((gcc / mass(24)) * vckm_13 * tanb * mass(5))
    g_yuk1_2_3(1) = ((gcc / mass(24)) * vckm_23 * (mass(4) / tanb))
    g_yuk1_2_3(2) = ((gcc / mass(24)) * vckm_23 * tanb * mass(5))
    g_yuk1_3_3(1) = ((gcc / mass(24)) * vckm_33 * (mass(6) / tanb))
    g_yuk1_3_3(2) = ((gcc / mass(24)) * vckm_33 * tanb * mass(5))
    g_yuk1_3_2(1) = ((gcc / mass(24)) * vckm_32 * (mass(6) / tanb))
    g_yuk1_3_2(2) = ((gcc / mass(24)) * vckm_32 * tanb * mass(3))
    g_yuk1_3_1(1) = ((gcc / mass(24)) * vckm_31 * (mass(6) / tanb))
    g_yuk1_3_1(2) = ((gcc / mass(24)) * vckm_31 * tanb * mass(1))
    g_yuk2_1_3(1) = conjg (g_yuk1_1_3(2))
    g_yuk2_1_3(2) = conjg (g_yuk1_1_3(1))
    g_yuk2_2_3(1) = conjg (g_yuk1_2_3(2))
    g_yuk2_2_3(2) = conjg (g_yuk1_2_3(1))
    g_yuk2_3_1(1) = conjg (g_yuk1_3_1(2))
    g_yuk2_3_1(2) = conjg (g_yuk1_3_1(1))
    g_yuk2_3_2(1) = conjg (g_yuk1_3_2(2))
    g_yuk2_3_2(2) = conjg (g_yuk1_3_2(1))
    g_yuk2_3_3(1) = conjg (g_yuk1_3_3(2))
    g_yuk2_3_3(2) = conjg (g_yuk1_3_3(1))
    gnzn_1_2(1) = (gz * vector0_12)
    gnzn_1_2(2) = (gz * axial0_12)
    gnzn_1_3(1) = (gz * vector0_13)
    gnzn_1_3(2) = (gz * axial0_13)
    gnzn_1_4(1) = (gz * vector0_14)
    gnzn_1_4(2) = (gz * axial0_14)
    gnzn_2_3(1) = (gz * vector0_23)
    gnzn_2_3(2) = (gz * axial0_23)
    gnzn_2_4(1) = (gz * vector0_24)
    gnzn_2_4(2) = (gz * axial0_24)
    gnzn_3_4(1) = (gz * vector0_34)
    gnzn_3_4(2) = (gz * axial0_34)
    gczc_1_1(1) = (gz * vp_11)
    gczc_1_1(2) = (gz * ap_11)
    gczc_1_2(1) = (gz * vp_12)
    gczc_1_2(2) = (gz * ap_12)
    gczc_2_1(1) = (gz * vp_21)
    gczc_2_1(2) = (gz * ap_21)
    gczc_2_2(1) = (gz * vp_22)
    gczc_2_2(2) = (gz * ap_22)
    gnwc_1_1(1) = (gcc * lnc_11)
    gnwc_1_1(2) = (gcc * rnc_11)
    g_nhc_1_1(1) = ((g / 2.0_default) * lnch_11)
    g_nhc_1_1(2) = ((g / 2.0_default) * rnch_11)
    gnwc_1_2(1) = (gcc * lnc_12)
    gnwc_1_2(2) = (gcc * rnc_12)
    g_nhc_1_2(1) = ((g / 2.0_default) * lnch_12)
    g_nhc_1_2(2) = ((g / 2.0_default) * rnch_12)
    gnwc_2_1(1) = (gcc * lnc_21)
    gnwc_2_1(2) = (gcc * rnc_21)
    g_nhc_2_1(1) = ((g / 2.0_default) * lnch_21)
    g_nhc_2_1(2) = ((g / 2.0_default) * rnch_21)
    gnwc_2_2(1) = (gcc * lnc_22)
    gnwc_2_2(2) = (gcc * rnc_22)
    g_nhc_2_2(1) = ((g / 2.0_default) * lnch_22)
    g_nhc_2_2(2) = ((g / 2.0_default) * rnch_22)
    gnwc_3_1(1) = (gcc * lnc_31)
    gnwc_3_1(2) = (gcc * rnc_31)
    g_nhc_3_1(1) = ((g / 2.0_default) * lnch_31)
    g_nhc_3_1(2) = ((g / 2.0_default) * rnch_31)
    gnwc_3_2(1) = (gcc * lnc_32)
    gnwc_3_2(2) = (gcc * rnc_32)
    g_nhc_3_2(1) = ((g / 2.0_default) * lnch_32)
    g_nhc_3_2(2) = ((g / 2.0_default) * rnch_32)
    gnwc_4_1(1) = (gcc * lnc_41)
    gnwc_4_1(2) = (gcc * rnc_41)
    g_nhc_4_1(1) = ((g / 2.0_default) * lnch_41)
    g_nhc_4_1(2) = ((g / 2.0_default) * rnch_41)
    gnwc_4_2(1) = (gcc * lnc_42)
    gnwc_4_2(2) = (gcc * rnc_42)
    g_nhc_4_2(1) = ((g / 2.0_default) * lnch_42)
    g_nhc_4_2(2) = ((g / 2.0_default) * rnch_42)
    gcwn_1_1(1) = (gcc * lcn_11)
    gcwn_1_1(2) = (gcc * rcn_11)
    g_chn_1_1(1) = ((g / 2.0_default) * conjg (rnch_11))
    g_chn_1_1(2) = ((g / 2.0_default) * conjg (lnch_11))
    gcwn_1_2(1) = (gcc * lcn_12)
    gcwn_1_2(2) = (gcc * rcn_12)
    g_chn_2_1(1) = ((g / 2.0_default) * conjg (rnch_21))
    g_chn_2_1(2) = ((g / 2.0_default) * conjg (lnch_21))
    gcwn_1_3(1) = (gcc * lcn_13)
    gcwn_1_3(2) = (gcc * rcn_13)
    g_chn_3_1(1) = ((g / 2.0_default) * conjg (rnch_31))
    g_chn_3_1(2) = ((g / 2.0_default) * conjg (lnch_31))
    gcwn_1_4(1) = (gcc * lcn_14)
    gcwn_1_4(2) = (gcc * rcn_14)
    g_chn_4_1(1) = ((g / 2.0_default) * conjg (rnch_41))
    g_chn_4_1(2) = ((g / 2.0_default) * conjg (lnch_41))
    gcwn_2_1(1) = (gcc * lcn_21)
    gcwn_2_1(2) = (gcc * rcn_21)
    g_chn_1_2(1) = ((g / 2.0_default) * conjg (rnch_12))
    g_chn_1_2(2) = ((g / 2.0_default) * conjg (lnch_12))
    gcwn_2_2(1) = (gcc * lcn_22)
    gcwn_2_2(2) = (gcc * rcn_22)
    g_chn_2_2(1) = ((g / 2.0_default) * conjg (rnch_22))
    g_chn_2_2(2) = ((g / 2.0_default) * conjg (lnch_22))
    gcwn_2_3(1) = (gcc * lcn_23)
    gcwn_2_3(2) = (gcc * rcn_23)
    g_chn_3_2(1) = ((g / 2.0_default) * conjg (rnch_32))
    g_chn_3_2(2) = ((g / 2.0_default) * conjg (lnch_32))
    gcwn_2_4(1) = (gcc * lcn_24)
    gcwn_2_4(2) = (gcc * rcn_24)
    g_chn_4_2(1) = ((g / 2.0_default) * conjg (rnch_42))
    g_chn_4_2(2) = ((g / 2.0_default) * conjg (lnch_42))
    gcicih1_1_1 = ((( - g ) / 2.0_default) * snnh1_11)
    gcicih2_1_1 = ((( - g ) / 2.0_default) * snnh2_11)
    gcicia_1_1 = ((( - g ) / 2.0_default) * pnna_11)
    gcicih1_1_2(1) = ((( - g ) / 2.0_default) * snnh1_12)
    gcicih1_1_2(2) = ((( - g ) / 2.0_default) * pnnh1_12)
    gcicih2_1_2(1) = ((( - g ) / 2.0_default) * snnh2_12)
    gcicih2_1_2(2) = ((( - g ) / 2.0_default) * pnnh2_12)
    gcicia_1_2(1) = ((( - g ) / 2.0_default) * snna_12)
    gcicia_1_2(2) = ((( - g ) / 2.0_default) * pnna_12)
    gcicih1_1_3(1) = ((( - g ) / 2.0_default) * snnh1_13)
    gcicih1_1_3(2) = ((( - g ) / 2.0_default) * pnnh1_13)
    gcicih2_1_3(1) = ((( - g ) / 2.0_default) * snnh2_13)
    gcicih2_1_3(2) = ((( - g ) / 2.0_default) * pnnh2_13)
    gcicia_1_3(1) = ((( - g ) / 2.0_default) * snna_13)
    gcicia_1_3(2) = ((( - g ) / 2.0_default) * pnna_13)
    gcicih1_1_4(1) = ((( - g ) / 2.0_default) * snnh1_14)
    gcicih1_1_4(2) = ((( - g ) / 2.0_default) * pnnh1_14)
    gcicih2_1_4(1) = ((( - g ) / 2.0_default) * snnh2_14)
    gcicih2_1_4(2) = ((( - g ) / 2.0_default) * pnnh2_14)
    gcicia_1_4(1) = ((( - g ) / 2.0_default) * snna_14)
    gcicia_1_4(2) = ((( - g ) / 2.0_default) * pnna_14)
    gcicih1_2_2 = ((( - g ) / 2.0_default) * snnh1_22)
    gcicih2_2_2 = ((( - g ) / 2.0_default) * snnh2_22)
    gcicia_2_2 = ((( - g ) / 2.0_default) * pnna_22)
    gcicih1_2_3(1) = ((( - g ) / 2.0_default) * snnh1_23)
    gcicih1_2_3(2) = ((( - g ) / 2.0_default) * pnnh1_23)
    gcicih2_2_3(1) = ((( - g ) / 2.0_default) * snnh2_23)
    gcicih2_2_3(2) = ((( - g ) / 2.0_default) * pnnh2_23)
 end subroutine setup_parameters13
 subroutine setup_parameters14 ()
!!! JR checked gch[x]h_[x]_[x]
    gcicia_2_3(1) = ((( - g ) / 2.0_default) * snna_23)
    gcicia_2_3(2) = ((( - g ) / 2.0_default) * pnna_23)
    gcicih1_2_4(1) = ((( - g ) / 2.0_default) * snnh1_24)
    gcicih1_2_4(2) = ((( - g ) / 2.0_default) * pnnh1_24)
    gcicih2_2_4(1) = ((( - g ) / 2.0_default) * snnh2_24)
    gcicih2_2_4(2) = ((( - g ) / 2.0_default) * pnnh2_24)
    gcicia_2_4(1) = ((( - g ) / 2.0_default) * snna_24)
    gcicia_2_4(2) = ((( - g ) / 2.0_default) * pnna_24)
    gcicih1_3_3 = ((( - g ) / 2.0_default) * snnh1_33)
    gcicih2_3_3 = ((( - g ) / 2.0_default) * snnh2_33)
    gcicia_3_3 = ((( - g ) / 2.0_default) * pnna_33)
    gcicih1_3_4(1) = ((( - g ) / 2.0_default) * snnh1_34)
    gcicih1_3_4(2) = ((( - g ) / 2.0_default) * pnnh1_34)
    gcicih2_3_4(1) = ((( - g ) / 2.0_default) * snnh2_34)
    gcicih2_3_4(2) = ((( - g ) / 2.0_default) * pnnh2_34)
    gcicia_3_4(1) = ((( - g ) / 2.0_default) * snna_34)
    gcicia_3_4(2) = ((( - g ) / 2.0_default) * pnna_34)
    gcicih1_4_4 = ((( - g ) / 2.0_default) * snnh1_44)
    gcicih2_4_4 = ((( - g ) / 2.0_default) * snnh2_44)
    gcicia_4_4 = ((( - g ) / 2.0_default) * pnna_44)
    gch1c_1_1 = (( - (g / sqrt (2.0_default))) * ((conjg (mu_11) *  &
      conjg (mv_12) * cosal) - (conjg (mu_12) * conjg (mv_11) * sinal)))
    gch2c_1_1 = (( - (g / sqrt (2.0_default))) * ((conjg (mu_12) *  &
      conjg (mv_11) * cosal) + (conjg (mu_11) * conjg (mv_12) * sinal)))
    gcac_1_1 = (imago * ( - (g /  &
      sqrt (2.0_default))) * ((mv_11 * mu_12 * sinbe) +  &
      (mv_12 * mu_11 * cosbe)))
    gch1c_1_2(1) = (( - gcc) * ((conjg (mu_11) *  &
      conjg (mv_22) * cosal) - (conjg (mu_12) * conjg (mv_21) * sinal)))
    gch1c_1_2(2) = (( - gcc) * ( &
      (mv_12 * mu_21 * cosal) - (mv_11 * mu_22 * sinal)))
    gch2c_1_2(1) = (( - gcc) * ((conjg (mu_12) *  &
      conjg (mv_21) * cosal) + (conjg (mu_11) * conjg (mv_22) * sinal)))
    gch2c_1_2(2) = (( - gcc) * ((mv_11 * mu_22 * cosal) &
      + (mv_12 * mu_21 * sinal)))
    gcac_1_2(1) = (imago * gcc * ((  &
      conjg (mu_12) * conjg (mv_21) * sinbe) + ( &
      conjg (mu_11) * conjg (mv_22) * cosbe)))
    gcac_1_2(2) = (( - imago) * gcc *  (( &
      mv_11 * mu_22 * sinbe) + (mv_12 * mu_21 * cosbe)))
    gch1c_2_1(1) = conjg (gch1c_1_2(2)) 
    gch1c_2_1(2) = conjg (gch1c_1_2(1))
    gch2c_2_1(1) = conjg (gch2c_1_2(2))
    gch2c_2_1(2) = conjg (gch2c_1_2(1))
    gcac_2_1(1) = conjg (gcac_1_2(2))
    gcac_2_1(2) = conjg (gcac_1_2(1))
    gch1c_2_2 = (( - (g / sqrt (2.0_default))) * ((conjg (mu_21) *  &
      conjg (mv_22) * cosal) - (conjg (mu_22) * conjg (mv_21) * sinal)))
    gch2c_2_2 = (( - (g / sqrt (2.0_default))) * ((conjg (mu_22) *  &
      conjg (mv_21) * cosal) + (conjg (mu_21) * conjg (mv_22) * sinal)))
    gcac_2_2 = (imago * ( - (g /  &
      sqrt (2.0_default))) * ((mv_21 * mu_22 * sinbe) +  &
      (mv_22 * mu_21 * cosbe)))
    g_yuk_ch1_sn1_3_c(1) = ((gcc * mass(15) * conjg (mu_12)) / (mass(24) &
      * cosbe))
    g_yuk_ch1_sn1_3_c(2) = ( - ((g * mv_11) / 2.0_default))
    g_yuk_ch1_sn1_3(1) = conjg (g_yuk_ch1_sn1_3_c(2))
    g_yuk_ch1_sn1_3(2) = conjg (g_yuk_ch1_sn1_3_c(1)) 
    g_yuk_ch2_sn1_3_c(1) = ((gcc * mass(15) * conjg (mu_22)) / (mass(24) &
      * cosbe))
    g_yuk_ch2_sn1_3_c(2) = ( - ((g * mv_21) / 2.0_default))
    g_yuk_ch2_sn1_3(1) = conjg (g_yuk_ch2_sn1_3_c(2))
    g_yuk_ch2_sn1_3(2) = conjg (g_yuk_ch2_sn1_3_c(1))
    g_yuk_ch1_sd1_1_3(1) = ((vckm_13 * gcc * mv_12 * mass(2) *  &
      conjg (mix_sd311)) / (mass(24) * sinbe))
    g_yuk_ch1_sd1_1_3(2) = (vckm_13 * gcc * (((conjg (mu_12) * mass(5) *  &
      conjg (mix_sd312)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd311))))
    g_yuk_ch1_sd1_1_3_c(1) = conjg (g_yuk_ch1_sd1_1_3(2))
    g_yuk_ch1_sd1_1_3_c(2) = conjg (g_yuk_ch1_sd1_1_3(1))
    g_yuk_ch1_su1_1_3(1) = (vckm_13 * gcc * (((conjg (mv_12) * mass(2) *  &
      conjg (mix_su112)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su111))))
    g_yuk_ch1_su1_1_3(2) = ((vckm_13 * gcc * mu_12 * mass(5) *  &
      conjg (mix_su111)) / (mass(24) * cosbe))
    g_yuk_ch1_su1_1_3_c(1) = conjg (g_yuk_ch1_su1_1_3(2))
    g_yuk_ch1_su1_1_3_c(2) = conjg (g_yuk_ch1_su1_1_3(1))
 end subroutine setup_parameters14
 subroutine setup_parameters15 ()
    g_yuk_ch1_sd1_2_3(1) = ((vckm_23 * gcc * mv_12 * mass(4) *  &
      conjg (mix_sd311)) / (mass(24) * sinbe))
    g_yuk_ch1_sd1_2_3(2) = (vckm_23 * gcc * (((conjg (mu_12) * mass(5) *  &
      conjg (mix_sd312)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd311))))
    g_yuk_ch1_sd1_2_3_c(1) = conjg (g_yuk_ch1_sd1_2_3(2))
    g_yuk_ch1_sd1_2_3_c(2) = conjg (g_yuk_ch1_sd1_2_3(1))
    g_yuk_ch1_su1_2_3(1) = (vckm_23 * gcc * (((conjg (mv_12) * mass(4) *  &
      conjg (mix_su212)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su211))))
    g_yuk_ch1_su1_2_3(2) = ((vckm_23 * gcc * mu_12 * mass(5) *  &
      conjg (mix_su211)) / (mass(24) * cosbe))
    g_yuk_ch1_su1_2_3_c(1) = conjg (g_yuk_ch1_su1_2_3(2))
    g_yuk_ch1_su1_2_3_c(2) = conjg (g_yuk_ch1_su1_2_3(1))
    g_yuk_ch1_sd1_3_3(1) = ((vckm_33 * gcc * mv_12 * mass(6) *  &
      conjg (mix_sd311)) / (mass(24) * sinbe))
    g_yuk_ch1_sd1_3_3(2) = (vckm_33 * gcc * (((conjg (mu_12) * mass(5) *  &
      conjg (mix_sd312)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd311))))
    g_yuk_ch1_sd1_3_3_c(1) = conjg (g_yuk_ch1_sd1_3_3(2))
    g_yuk_ch1_sd1_3_3_c(2) = conjg (g_yuk_ch1_sd1_3_3(1))
    g_yuk_ch1_su1_3_3(1) = (vckm_33 * gcc * (((conjg (mv_12) * mass(6) *  &
      conjg (mix_su312)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su311))))
    g_yuk_ch1_su1_3_3(2) = ((vckm_33 * gcc * mu_12 * mass(5) *  &
      conjg (mix_su311)) / (mass(24) * cosbe))
    g_yuk_ch1_su1_3_3_c(1) = conjg (g_yuk_ch1_su1_3_3(2))
    g_yuk_ch1_su1_3_3_c(2) = conjg (g_yuk_ch1_su1_3_3(1))
    g_yuk_ch1_sd1_3_2(1) = ((vckm_32 * gcc * mv_12 * mass(6) *  &
      conjg (mix_sd211)) / (mass(24) * sinbe))
    g_yuk_ch1_sd1_3_2(2) = (vckm_32 * gcc * (((conjg (mu_12) * mass(3) *  &
      conjg (mix_sd212)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd211))))
    g_yuk_ch1_sd1_3_2_c(1) = conjg (g_yuk_ch1_sd1_3_2(2))
    g_yuk_ch1_sd1_3_2_c(2) = conjg (g_yuk_ch1_sd1_3_2(1))
    g_yuk_ch1_su1_3_2(1) = (vckm_32 * gcc * (((conjg (mv_12) * mass(6) *  &
      conjg (mix_su312)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su311))))
    g_yuk_ch1_su1_3_2(2) = ((vckm_32 * gcc * mu_12 * mass(3) *  &
      conjg (mix_su311)) / (mass(24) * cosbe))
    g_yuk_ch1_su1_3_2_c(1) = conjg (g_yuk_ch1_su1_3_2(2))
    g_yuk_ch1_su1_3_2_c(2) = conjg (g_yuk_ch1_su1_3_2(1))
    g_yuk_ch1_sd1_3_1(1) = ((vckm_31 * gcc * mv_12 * mass(6) *  &
      conjg (mix_sd111)) / (mass(24) * sinbe))
    g_yuk_ch1_sd1_3_1(2) = (vckm_31 * gcc * (((conjg (mu_12) * mass(1) *  &
      conjg (mix_sd112)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd111))))
    g_yuk_ch1_sd1_3_1_c(1) = conjg (g_yuk_ch1_sd1_3_1(2))
    g_yuk_ch1_sd1_3_1_c(2) = conjg (g_yuk_ch1_sd1_3_1(1))
    g_yuk_ch1_su1_3_1(1) = (vckm_31 * gcc * (((conjg (mv_12) * mass(6) *  &
      conjg (mix_su312)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su311))))
    g_yuk_ch1_su1_3_1(2) = ((vckm_31 * gcc * mu_12 * mass(1) *  &
      conjg (mix_su311)) / (mass(24) * cosbe))
    g_yuk_ch1_su1_3_1_c(1) = conjg (g_yuk_ch1_su1_3_1(2))
    g_yuk_ch1_su1_3_1_c(2) = conjg (g_yuk_ch1_su1_3_1(1))
    g_yuk_ch1_sd2_1_3(1) = ((vckm_13 * gcc * mv_12 * mass(2) *  &
      conjg (mix_sd321)) / (mass(24) * sinbe))
    g_yuk_ch1_sd2_1_3(2) = (vckm_13 * gcc * (((conjg (mu_12) * mass(5) *  &
      conjg (mix_sd322)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd321))))
    g_yuk_ch1_sd2_1_3_c(1) = conjg (g_yuk_ch1_sd2_1_3(2))
    g_yuk_ch1_sd2_1_3_c(2) = conjg (g_yuk_ch1_sd2_1_3(1))
    g_yuk_ch1_su2_1_3(1) = (vckm_13 * gcc * (((conjg (mv_12) * mass(2) *  &
      conjg (mix_su122)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su121))))
    g_yuk_ch1_su2_1_3(2) = ((vckm_13 * gcc * mu_12 * mass(5) *  &
      conjg (mix_su121)) / (mass(24) * cosbe))
    g_yuk_ch1_su2_1_3_c(1) = conjg (g_yuk_ch1_su2_1_3(2))
    g_yuk_ch1_su2_1_3_c(2) = conjg (g_yuk_ch1_su2_1_3(1))
    g_yuk_ch1_sd2_2_3(1) = ((vckm_23 * gcc * mv_12 * mass(4) *  &
      conjg (mix_sd321)) / (mass(24) * sinbe))
    g_yuk_ch1_sd2_2_3(2) = (vckm_23 * gcc * (((conjg (mu_12) * mass(5) *  &
      conjg (mix_sd322)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd321))))
    g_yuk_ch1_sd2_2_3_c(1) = conjg (g_yuk_ch1_sd2_2_3(2))
    g_yuk_ch1_sd2_2_3_c(2) = conjg (g_yuk_ch1_sd2_2_3(1))
    g_yuk_ch1_su2_2_3(1) = (vckm_23 * gcc * (((conjg (mv_12) * mass(4) *  &
      conjg (mix_su222)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su221))))
    g_yuk_ch1_su2_2_3(2) = ((vckm_23 * gcc * mu_12 * mass(5) *  &
      conjg (mix_su221)) / (mass(24) * cosbe))
    g_yuk_ch1_su2_2_3_c(1) = conjg (g_yuk_ch1_su2_2_3(2))
    g_yuk_ch1_su2_2_3_c(2) = conjg (g_yuk_ch1_su2_2_3(1))
    g_yuk_ch1_sd2_3_3(1) = ((vckm_33 * gcc * mv_12 * mass(6) *  &
      conjg (mix_sd321)) / (mass(24) * sinbe))
    g_yuk_ch1_sd2_3_3(2) = (vckm_33 * gcc * (((conjg (mu_12) * mass(5) *  &
      conjg (mix_sd322)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd321))))
    g_yuk_ch1_sd2_3_3_c(1) = conjg (g_yuk_ch1_sd2_3_3(2))
    g_yuk_ch1_sd2_3_3_c(2) = conjg (g_yuk_ch1_sd2_3_3(1))
    g_yuk_ch1_su2_3_3(1) = (vckm_33 * gcc * (((conjg (mv_12) * mass(6) *  &
      conjg (mix_su322)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su321))))
    g_yuk_ch1_su2_3_3(2) = ((vckm_33 * gcc * mu_12 * mass(5) *  &
      conjg (mix_su321)) / (mass(24) * cosbe))
    g_yuk_ch1_su2_3_3_c(1) = conjg (g_yuk_ch1_su2_3_3(2))
    g_yuk_ch1_su2_3_3_c(2) = conjg (g_yuk_ch1_su2_3_3(1))
    g_yuk_ch1_sd2_3_2(1) = ((vckm_32 * gcc * mv_12 * mass(6) *  &
      conjg (mix_sd221)) / (mass(24) * sinbe))
    g_yuk_ch1_sd2_3_2(2) = (vckm_32 * gcc * (((conjg (mu_12) * mass(3) *  &
      conjg (mix_sd222)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd221))))
    g_yuk_ch1_sd2_3_2_c(1) = conjg (g_yuk_ch1_sd2_3_2(2))
    g_yuk_ch1_sd2_3_2_c(2) = conjg (g_yuk_ch1_sd2_3_2(1))
    g_yuk_ch1_su2_3_2(1) = (vckm_32 * gcc * (((conjg (mv_12) * mass(6) *  &
      conjg (mix_su322)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su321))))
    g_yuk_ch1_su2_3_2(2) = ((vckm_32 * gcc * mu_12 * mass(3) *  &
      conjg (mix_su321)) / (mass(24) * cosbe))
    g_yuk_ch1_su2_3_2_c(1) = conjg (g_yuk_ch1_su2_3_3(2))
    g_yuk_ch1_su2_3_2_c(2) = conjg (g_yuk_ch1_su2_3_3(1))
    g_yuk_ch1_sd2_3_1(1) = ((vckm_31 * gcc * mv_12 * mass(6) *  &
      conjg (mix_sd121)) / (mass(24) * sinbe))
    g_yuk_ch1_sd2_3_1(2) = (vckm_31 * gcc * (((conjg (mu_12) * mass(1) *  &
      conjg (mix_sd122)) / (mass(24) * cosbe)) - (conjg (mu_11) *  &
      sqrt (2.0_default) * conjg (mix_sd121))))
    g_yuk_ch1_sd2_3_1_c(1) = conjg (g_yuk_ch1_sd2_3_1(2))
    g_yuk_ch1_sd2_3_1_c(2) = conjg (g_yuk_ch1_sd2_3_1(1))
    g_yuk_ch1_su2_3_1(1) = (vckm_31 * gcc * (((conjg (mv_12) * mass(6) *  &
      conjg (mix_su322)) / (mass(24) * sinbe)) - (conjg (mv_11) *  &
      sqrt (2.0_default) * conjg (mix_su321))))
    g_yuk_ch1_su2_3_1(2) = ((vckm_31 * gcc * mu_12 * mass(1) *  &
      conjg (mix_su321)) / (mass(24) * cosbe))
    g_yuk_ch1_su2_3_1_c(1) = conjg (g_yuk_ch1_su2_3_1(2))
    g_yuk_ch1_su2_3_1_c(2) = conjg (g_yuk_ch1_su2_3_1(1))
    g_yuk_ch2_sd1_1_3(1) = ((vckm_13 * gcc * mv_22 * mass(2) *  &
      conjg (mix_sd311)) / (mass(24) * sinbe))
    g_yuk_ch2_sd1_1_3(2) = (vckm_13 * gcc * (((conjg (mu_22) * mass(5) *  &
      conjg (mix_sd312)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd311))))
    g_yuk_ch2_sd1_1_3_c(1) = conjg (g_yuk_ch2_sd1_1_3(2))
    g_yuk_ch2_sd1_1_3_c(2) = conjg (g_yuk_ch2_sd1_1_3(1))
    g_yuk_ch2_su1_1_3(1) = (vckm_13 * gcc * (((conjg (mv_22) * mass(2) *  &
      conjg (mix_su112)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su111))))
    g_yuk_ch2_su1_1_3(2) = ((vckm_13 * gcc * mu_22 * mass(5) *  &
      conjg (mix_su111)) / (mass(24) * cosbe))
    g_yuk_ch2_su1_1_3_c(1) = conjg (g_yuk_ch2_su1_1_3(2))
    g_yuk_ch2_su1_1_3_c(2) = conjg (g_yuk_ch2_su1_1_3(1))
    g_yuk_ch2_sd1_2_3(1) = ((vckm_23 * gcc * mv_22 * mass(4) *  &
      conjg (mix_sd311)) / (mass(24) * sinbe))
    g_yuk_ch2_sd1_2_3(2) = (vckm_23 * gcc * (((conjg (mu_22) * mass(5) *  &
      conjg (mix_sd312)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd311))))
    g_yuk_ch2_sd1_2_3_c(1) = conjg (g_yuk_ch2_sd1_2_3(2))
    g_yuk_ch2_sd1_2_3_c(2) = conjg (g_yuk_ch2_sd1_2_3(1))
    g_yuk_ch2_su1_2_3(1) = (vckm_23 * gcc * (((conjg (mv_22) * mass(4) *  &
      conjg (mix_su212)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su211))))
    g_yuk_ch2_su1_2_3(2) = ((vckm_23 * gcc * mu_22 * mass(5) *  &
      conjg (mix_su211)) / (mass(24) * cosbe))
    g_yuk_ch2_su1_2_3_c(1) = conjg (g_yuk_ch2_su1_2_3(2))
    g_yuk_ch2_su1_2_3_c(2) = conjg (g_yuk_ch2_su1_2_3(1))
    g_yuk_ch2_sd1_3_3(1) = ((vckm_33 * gcc * mv_22 * mass(6) *  &
      conjg (mix_sd311)) / (mass(24) * sinbe))
    g_yuk_ch2_sd1_3_3(2) = (vckm_33 * gcc * (((conjg (mu_22) * mass(5) *  &
      conjg (mix_sd312)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd311))))
    g_yuk_ch2_sd1_3_3_c(1) = conjg (g_yuk_ch2_sd1_3_3(2))
    g_yuk_ch2_sd1_3_3_c(2) = conjg (g_yuk_ch2_sd1_3_3(1))
    g_yuk_ch2_su1_3_3(1) = (vckm_33 * gcc * (((conjg (mv_22) * mass(6) *  &
      conjg (mix_su312)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su311))))
    g_yuk_ch2_su1_3_3(2) = ((vckm_33 * gcc * mu_22 * mass(5) *  &
      conjg (mix_su311)) / (mass(24) * cosbe))
    g_yuk_ch2_su1_3_3_c(1) = conjg (g_yuk_ch2_su1_3_3(2))
    g_yuk_ch2_su1_3_3_c(2) = conjg (g_yuk_ch2_su1_3_3(1))
    g_yuk_ch2_sd1_3_2(1) = ((vckm_32 * gcc * mv_22 * mass(6) *  &
      conjg (mix_sd211)) / (mass(24) * sinbe))
    g_yuk_ch2_sd1_3_2(2) = (vckm_32 * gcc * (((conjg (mu_22) * mass(3) *  &
      conjg (mix_sd212)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd211))))
    g_yuk_ch2_sd1_3_2_c(1) = conjg (g_yuk_ch2_sd1_3_2(2))
    g_yuk_ch2_sd1_3_2_c(2) = conjg (g_yuk_ch2_sd1_3_2(1))
    g_yuk_ch2_su1_3_2(1) = (vckm_32 * gcc * (((conjg (mv_22) * mass(6) *  &
      conjg (mix_su312)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su311))))
    g_yuk_ch2_su1_3_2(2) = ((vckm_32 * gcc * mu_22 * mass(3) *  &
      conjg (mix_su311)) / (mass(24) * cosbe))
    g_yuk_ch2_su1_3_2_c(1) = conjg (g_yuk_ch2_su1_3_2(2))
    g_yuk_ch2_su1_3_2_c(2) = conjg (g_yuk_ch2_su1_3_2(1))
    g_yuk_ch2_sd1_3_1(1) = ((vckm_31 * gcc * mv_22 * mass(6) *  &
      conjg (mix_sd111)) / (mass(24) * sinbe))
    g_yuk_ch2_sd1_3_1(2) = (vckm_31 * gcc * (((conjg (mu_22) * mass(1) *  &
      conjg (mix_sd112)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd111))))
    g_yuk_ch2_sd1_3_1_c(1) = conjg (g_yuk_ch2_sd1_3_1(2))
    g_yuk_ch2_sd1_3_1_c(2) = conjg (g_yuk_ch2_sd1_3_1(1))
    g_yuk_ch2_su1_3_1(1) = (vckm_31 * gcc * (((conjg (mv_22) * mass(6) *  &
      conjg (mix_su312)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su311))))
    g_yuk_ch2_su1_3_1(2) = ((vckm_31 * gcc * mu_22 * mass(1) *  &
      conjg (mix_su311)) / (mass(24) * cosbe))
    g_yuk_ch2_su1_3_1_c(1) = conjg (g_yuk_ch2_su1_3_1(2))
    g_yuk_ch2_su1_3_1_c(2) = conjg (g_yuk_ch2_su1_3_1(1))
    g_yuk_ch2_sd2_1_3(1) = ((vckm_13 * gcc * mv_22 * mass(2) *  &
      conjg (mix_sd321)) / (mass(24) * sinbe))
    g_yuk_ch2_sd2_1_3(2) = (vckm_13 * gcc * (((conjg (mu_22) * mass(5) *  &
      conjg (mix_sd322)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd321))))
    g_yuk_ch2_sd2_1_3_c(1) = conjg (g_yuk_ch2_sd2_1_3(2))
    g_yuk_ch2_sd2_1_3_c(2) = conjg (g_yuk_ch2_sd2_1_3(1))
    g_yuk_ch2_su2_1_3(1) = (vckm_13 * gcc * (((conjg (mv_22) * mass(2) *  &
      conjg (mix_su122)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su121))))
    g_yuk_ch2_su2_1_3(2) = ((vckm_13 * gcc * mu_22 * mass(5) *  &
      conjg (mix_su121)) / (mass(24) * cosbe))
    g_yuk_ch2_su2_1_3_c(1) = conjg (g_yuk_ch2_su2_1_3(2))
    g_yuk_ch2_su2_1_3_c(2) = conjg (g_yuk_ch2_su2_1_3(1))
    g_yuk_ch2_sd2_2_3(1) = ((vckm_23 * gcc * mv_22 * mass(4) *  &
      conjg (mix_sd321)) / (mass(24) * sinbe))
    g_yuk_ch2_sd2_2_3(2) = (vckm_23 * gcc * (((conjg (mu_22) * mass(5) *  &
      conjg (mix_sd322)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd321))))
    g_yuk_ch2_sd2_2_3_c(1) = conjg (g_yuk_ch2_sd2_2_3(2))
    g_yuk_ch2_sd2_2_3_c(2) = conjg (g_yuk_ch2_sd2_2_3(1))
    g_yuk_ch2_su2_2_3(1) = (vckm_23 * gcc * (((conjg (mv_22) * mass(4) *  &
      conjg (mix_su222)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su221))))
    g_yuk_ch2_su2_2_3(2) = ((vckm_23 * gcc * mu_22 * mass(5) *  &
      conjg (mix_su221)) / (mass(24) * cosbe))
    g_yuk_ch2_su2_2_3_c(1) = conjg (g_yuk_ch2_su2_2_3(2))
    g_yuk_ch2_su2_2_3_c(2) = conjg (g_yuk_ch2_su2_2_3(1))
    g_yuk_ch2_sd2_3_3(1) = ((vckm_33 * gcc * mv_22 * mass(6) *  &
      conjg (mix_sd321)) / (mass(24) * sinbe))
    g_yuk_ch2_sd2_3_3(2) = (vckm_33 * gcc * (((conjg (mu_22) * mass(5) *  &
      conjg (mix_sd322)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd321))))
    g_yuk_ch2_sd2_3_3_c(1) = conjg (g_yuk_ch2_sd2_3_3(2))
    g_yuk_ch2_sd2_3_3_c(2) = conjg (g_yuk_ch2_sd2_3_3(1))
    g_yuk_ch2_su2_3_3(1) = (vckm_33 * gcc * (((conjg (mv_22) * mass(6) *  &
      conjg (mix_su322)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su321))))
    g_yuk_ch2_su2_3_3(2) = ((vckm_33 * gcc * mu_22 * mass(5) *  &
      conjg (mix_su321)) / (mass(24) * cosbe))
    g_yuk_ch2_su2_3_3_c(1) = conjg (g_yuk_ch2_su2_3_3(2))
    g_yuk_ch2_su2_3_3_c(2) = conjg (g_yuk_ch2_su2_3_3(1))
    g_yuk_ch2_sd2_3_2(1) = ((vckm_32 * gcc * mv_22 * mass(6) *  &
      conjg (mix_sd221)) / (mass(24) * sinbe))
    g_yuk_ch2_sd2_3_2(2) = (vckm_32 * gcc * (((conjg (mu_22) * mass(3) *  &
      conjg (mix_sd222)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd221))))
    g_yuk_ch2_sd2_3_2_c(1) = conjg (g_yuk_ch2_sd2_3_2(2))
    g_yuk_ch2_sd2_3_2_c(2) = conjg (g_yuk_ch2_sd2_3_2(1))
    g_yuk_ch2_su2_3_2(1) = (vckm_32 * gcc * (((conjg (mv_22) * mass(6) *  &
      conjg (mix_su322)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su321))))
    g_yuk_ch2_su2_3_2(2) = ((vckm_32 * gcc * mu_22 * mass(3) *  &
      conjg (mix_su321)) / (mass(24) * cosbe))
    g_yuk_ch2_su2_3_2_c(1) = conjg (g_yuk_ch2_su2_3_2(2))
    g_yuk_ch2_su2_3_2_c(2) = conjg (g_yuk_ch2_su2_3_2(1))
    g_yuk_ch2_sd2_3_1(1) = ((vckm_31 * gcc * mv_22 * mass(6) *  &
      conjg (mix_sd121)) / (mass(24) * sinbe))
    g_yuk_ch2_sd2_3_1(2) = (vckm_31 * gcc * (((conjg (mu_22) * mass(1) *  &
      conjg (mix_sd122)) / (mass(24) * cosbe)) - (conjg (mu_21) *  &
      sqrt (2.0_default) * conjg (mix_sd121))))
    g_yuk_ch2_sd2_3_1_c(1) = conjg (g_yuk_ch2_sd2_3_1(2))
    g_yuk_ch2_sd2_3_1_c(2) = conjg (g_yuk_ch2_sd2_3_1(1))
    g_yuk_ch2_su2_3_1(1) = (vckm_31 * gcc * (((conjg (mv_22) * mass(6) *  &
      conjg (mix_su322)) / (mass(24) * sinbe)) - (conjg (mv_21) *  &
      sqrt (2.0_default) * conjg (mix_su321))))
    g_yuk_ch2_su2_3_1(2) = ((vckm_31 * gcc * mu_22 * mass(1) *  &
      conjg (mix_su321)) / (mass(24) * cosbe))
    g_yuk_ch2_su2_3_1_c(1) = conjg (g_yuk_ch2_su2_3_1(2))
    g_yuk_ch2_su2_3_1_c(2) = conjg (g_yuk_ch2_su2_3_1(1))
 end subroutine setup_parameters15
 subroutine setup_parameters16 ()
    g_yuk_n1_sl1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_11) * (sinthw / costhw) * mix_sl312) &
      + ((conjg (mn_13) * mass(15) * mix_sl311) / (mass(24) * cosbe)))))
    g_yuk_n1_sl1_3(2) = (gcc * ((1.0_default * (mn_12 + (1.0_default *  &
      (sinthw / costhw) * mn_11)) * mix_sl311) - ( &
      (mn_13 * mass(15) * mix_sl312) / (mass(24) * cosbe))))
    g_yuk_n1_sl1_3_c(1) = conjg (g_yuk_n1_sl1_3(2))
    g_yuk_n1_sl1_3_c(2) = conjg (g_yuk_n1_sl1_3(1))
    g_yuk_n1_su1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_11) *  &
      (sinthw / costhw) * mix_su312) + ((conjg (mn_14) * mass(6) * mix_su311) /  &
      (mass(24) * sinbe)))))
    g_yuk_n1_su1_3(2) = (gcc * (((-1.0_default) * (mn_12 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_11)) * mix_su311) - ( &
      (mn_14 * mass(6) * mix_su312) / (mass(24) * sinbe))))
    g_yuk_n1_su1_3_c(1) = conjg (g_yuk_n1_su1_3(2))
    g_yuk_n1_su1_3_c(2) = conjg (g_yuk_n1_su1_3(1))
    g_yuk_n1_sd1_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_11) *  &
      (sinthw / costhw) * mix_sd312) + ((conjg (mn_13) * mass(5) * mix_sd311) /  &
      (mass(24) * cosbe)))))
    g_yuk_n1_sd1_3(2) = (gcc * ((1.0_default * (mn_12 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_11)) * mix_sd311) - ( &
      (mn_13 * mass(5) * mix_sd312) / (mass(24) * cosbe))))
    g_yuk_n1_sd1_3_c(1) = conjg (g_yuk_n1_sd1_3(2))
    g_yuk_n1_sd1_3_c(2) = conjg (g_yuk_n1_sd1_3(1))
    g_yuk_n2_sl1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_21) * (sinthw / costhw) * mix_sl312) + ( &
      (conjg (mn_23) * mass(15) * mix_sl311) / (mass(24) * cosbe)))))
    g_yuk_n2_sl1_3(2) = (gcc * ((1.0_default * (mn_22 + (1.0_default *  &
      (sinthw / costhw) * mn_21)) * mix_sl311) - ( &
      (mn_23 * mass(15) * mix_sl312) / (mass(24) * cosbe))))
    g_yuk_n2_sl1_3_c(1) = conjg (g_yuk_n2_sl1_3(2))
    g_yuk_n2_sl1_3_c(2) = conjg (g_yuk_n2_sl1_3(1))
    g_yuk_n2_su1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_21) *  &
      (sinthw / costhw) * mix_su312) + ((conjg (mn_24) * mass(6) * mix_su311) /  &
      (mass(24) * sinbe)))))
    g_yuk_n2_su1_3(2) = (gcc * (((-1.0_default) * (mn_22 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_21)) * mix_su311) - ( &
      (mn_24 * mass(6) * mix_su312) / (mass(24) * sinbe))))
    g_yuk_n2_su1_3_c(1) = conjg (g_yuk_n2_su1_3(2))
    g_yuk_n2_su1_3_c(2) = conjg (g_yuk_n2_su1_3(1))
    g_yuk_n2_sd1_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_21) *  &
      (sinthw / costhw) * mix_sd312) + ((conjg (mn_23) * mass(5) * mix_sd311) /  &
      (mass(24) * cosbe)))))
    g_yuk_n2_sd1_3(2) = (gcc * ((1.0_default * (mn_22 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_21)) * mix_sd311) - ( &
      (mn_23 * mass(5) * mix_sd312) / (mass(24) * cosbe))))
    g_yuk_n2_sd1_3_c(1) = conjg (g_yuk_n2_sd1_3(2))
    g_yuk_n2_sd1_3_c(2) = conjg (g_yuk_n2_sd1_3(1))
    g_yuk_n3_sl1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_31) * (sinthw / costhw) * mix_sl312) + ( &
      (conjg (mn_33) * mass(15) * mix_sl311) / (mass(24) * cosbe)))))
    g_yuk_n3_sl1_3(2) = (gcc * ((1.0_default * (mn_32 + (1.0_default *  &
      (sinthw / costhw) * mn_31)) * mix_sl311) - ( &
      (mn_33 * mass(15) * mix_sl312) / (mass(24) * cosbe))))
    g_yuk_n3_sl1_3_c(1) = conjg (g_yuk_n3_sl1_3(2))
    g_yuk_n3_sl1_3_c(2) = conjg (g_yuk_n3_sl1_3(1))
    g_yuk_n3_su1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_31) *  &
      (sinthw / costhw) * mix_su312) + ((conjg (mn_34) * mass(6) * mix_su311) /  &
      (mass(24) * sinbe)))))
    g_yuk_n3_su1_3(2) = (gcc * (((-1.0_default) * (mn_32 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_31)) * mix_su311) - ( &
      (mn_34 * mass(6) * mix_su312) / (mass(24) * sinbe))))
    g_yuk_n3_su1_3_c(1) = conjg (g_yuk_n3_su1_3(2))
    g_yuk_n3_su1_3_c(2) = conjg (g_yuk_n3_su1_3(1))
    g_yuk_n3_sd1_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_31) *  &
      (sinthw / costhw) * mix_sd312) + ((conjg (mn_33) * mass(5) * mix_sd311) /  &
      (mass(24) * cosbe)))))
    g_yuk_n3_sd1_3(2) = (gcc * ((1.0_default * (mn_32 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_31)) * mix_sd311) - ( &
      (mn_33 * mass(5) * mix_sd312) / (mass(24) * cosbe))))
    g_yuk_n3_sd1_3_c(1) = conjg (g_yuk_n3_sd1_3(2))
    g_yuk_n3_sd1_3_c(2) = conjg (g_yuk_n3_sd1_3(1))
    g_yuk_n4_sl1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_41) * (sinthw / costhw) * mix_sl312) + ( &
      (conjg (mn_43) * mass(15) * mix_sl311) / (mass(24) * cosbe)))))
    g_yuk_n4_sl1_3(2) = (gcc * ((1.0_default * (mn_42 + (1.0_default *  &
      (sinthw / costhw) * mn_41)) * mix_sl311) - ( &
      (mn_43 * mass(15) * mix_sl312) / (mass(24) * cosbe))))
    g_yuk_n4_sl1_3_c(1) = conjg (g_yuk_n4_sl1_3(2))
    g_yuk_n4_sl1_3_c(2) = conjg (g_yuk_n4_sl1_3(1))
    g_yuk_n4_su1_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_41) *  &
      (sinthw / costhw) * mix_su312) + ((conjg (mn_44) * mass(6) * mix_su311) /  &
      (mass(24) * sinbe)))))
    g_yuk_n4_su1_3(2) = (gcc * (((-1.0_default) * (mn_42 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_41)) * mix_su311) - ( &
      (mn_44 * mass(6) * mix_su312) / (mass(24) * sinbe))))
    g_yuk_n4_su1_3_c(1) = conjg (g_yuk_n4_su1_3(2))
    g_yuk_n4_su1_3_c(2) = conjg (g_yuk_n4_su1_3(1))
    g_yuk_n4_sd1_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_41) *  &
      (sinthw / costhw) * mix_sd312) + ((conjg (mn_43) * mass(5) * mix_sd311) /  &
      (mass(24) * cosbe)))))
    g_yuk_n4_sd1_3(2) = (gcc * ((1.0_default * (mn_42 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_41)) * mix_sd311) - ( &
      (mn_43 * mass(5) * mix_sd312) / (mass(24) * cosbe))))
    g_yuk_n4_sd1_3_c(1) = conjg (g_yuk_n4_sd1_3(2))
    g_yuk_n4_sd1_3_c(2) = conjg (g_yuk_n4_sd1_3(1))
    g_yuk_n1_sl2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_11) * (sinthw / costhw) * mix_sl322) + ( &
      (conjg (mn_13) * mass(15) * mix_sl321) / (mass(24) * cosbe)))))
    g_yuk_n1_sl2_3(2) = (gcc * ((1.0_default * (mn_12 + (1.0_default *  &
      (sinthw / costhw) * mn_11)) * mix_sl321) - ( &
      (mn_13 * mass(15) * mix_sl322) / (mass(24) * cosbe))))
    g_yuk_n1_sl2_3_c(1) = conjg (g_yuk_n1_sl2_3(2))
    g_yuk_n1_sl2_3_c(2) = conjg (g_yuk_n1_sl2_3(1))
    g_yuk_n1_su2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_11) *  &
      (sinthw / costhw) * mix_su322) + ((conjg (mn_14) * mass(6) * mix_su321) /  &
      (mass(24) * sinbe)))))
    g_yuk_n1_su2_3(2) = (gcc * (((-1.0_default) * (mn_12 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_11)) * mix_su321) - ( &
      (mn_14 * mass(6) * mix_su322) / (mass(24) * sinbe))))
    g_yuk_n1_su2_3_c(1) = conjg (g_yuk_n1_su2_3(2))
    g_yuk_n1_su2_3_c(2) = conjg (g_yuk_n1_su2_3(1))
    g_yuk_n1_sd2_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_11) *  &
      (sinthw / costhw) * mix_sd322) + ((conjg (mn_13) * mass(5) * mix_sd321) /  &
      (mass(24) * cosbe)))))
    g_yuk_n1_sd2_3(2) = (gcc * ((1.0_default * (mn_12 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_11)) * mix_sd321) - ( &
      (mn_13 * mass(5) * mix_sd322) / (mass(24) * cosbe))))
    g_yuk_n1_sd2_3_c(1) = conjg (g_yuk_n1_sd2_3(2))
    g_yuk_n1_sd2_3_c(2) = conjg (g_yuk_n1_sd2_3(1))
    g_yuk_n2_sl2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_21) * (sinthw / costhw) * mix_sl322) + ( &
      (conjg (mn_23) * mass(15) * mix_sl321) / (mass(24) * cosbe)))))
    g_yuk_n2_sl2_3(2) = (gcc * ((1.0_default * (mn_22 + (1.0_default *  &
      (sinthw / costhw) * mn_21)) * mix_sl321) - ( &
      (mn_23 * mass(15) * mix_sl322) / (mass(24) * cosbe))))
    g_yuk_n2_sl2_3_c(1) = conjg (g_yuk_n2_sl2_3(2))
    g_yuk_n2_sl2_3_c(2) = conjg (g_yuk_n2_sl2_3(1))
    g_yuk_n2_su2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_21) *  &
      (sinthw / costhw) * mix_su322) + ((conjg (mn_24) * mass(6) * mix_su321) /  &
      (mass(24) * sinbe)))))
    g_yuk_n2_su2_3(2) = (gcc * (((-1.0_default) * (mn_22 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_21)) * mix_su321) - ( &
      (mn_24 * mass(6) * mix_su322) / (mass(24) * sinbe))))
    g_yuk_n2_su2_3_c(1) = conjg (g_yuk_n2_su2_3(2))
    g_yuk_n2_su2_3_c(2) = conjg (g_yuk_n2_su2_3(1))
    g_yuk_n2_sd2_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_21) *  &
      (sinthw / costhw) * mix_sd322) + ((conjg (mn_23) * mass(5) * mix_sd321) /  &
      (mass(24) * cosbe)))))
    g_yuk_n2_sd2_3(2) = (gcc * ((1.0_default * (mn_22 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_21)) * mix_sd321) - ( &
      (mn_23 * mass(5) * mix_sd322) / (mass(24) * cosbe))))
    g_yuk_n2_sd2_3_c(1) = conjg (g_yuk_n2_sd2_3(2))
    g_yuk_n2_sd2_3_c(2) = conjg (g_yuk_n2_sd2_3(1))
    g_yuk_n3_sl2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_31) * (sinthw / costhw) * mix_sl322) + ( &
      (conjg (mn_33) * mass(15) * mix_sl321) / (mass(24) * cosbe)))))
    g_yuk_n3_sl2_3(2) = (gcc * ((1.0_default * (mn_32 + (1.0_default *  &
      (sinthw / costhw) * mn_31)) * mix_sl321) - ( &
      (mn_33 * mass(15) * mix_sl322) / (mass(24) * cosbe))))
    g_yuk_n3_sl2_3_c(1) = conjg (g_yuk_n3_sl2_3(2))
    g_yuk_n3_sl2_3_c(2) = conjg (g_yuk_n3_sl2_3(1))
    g_yuk_n3_su2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_31) *  &
      (sinthw / costhw) * mix_su322) + ((conjg (mn_34) * mass(6) * mix_su321) /  &
      (mass(24) * sinbe)))))
    g_yuk_n3_su2_3(2) = (gcc * (((-1.0_default) * (mn_32 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_31)) * mix_su321) - ( &
      (mn_34 * mass(6) * mix_su322) / (mass(24) * sinbe))))
    g_yuk_n3_su2_3_c(1) = conjg (g_yuk_n3_su2_3(2))
    g_yuk_n3_su2_3_c(2) = conjg (g_yuk_n3_su2_3(1))
    g_yuk_n3_sd2_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_31) *  &
      (sinthw / costhw) * mix_sd322) + ((conjg (mn_33) * mass(5) * mix_sd321) /  &
      (mass(24) * cosbe)))))
    g_yuk_n3_sd2_3(2) = (gcc * ((1.0_default * (mn_32 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_31)) * mix_sd321) - ( &
      (mn_33 * mass(5) * mix_sd322) / (mass(24) * cosbe))))
    g_yuk_n3_sd2_3_c(1) = conjg (g_yuk_n3_sd2_3(2))
    g_yuk_n3_sd2_3_c(2) = conjg (g_yuk_n3_sd2_3(1))
    g_yuk_n4_sl2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_lep) * conjg (mn_41) * (sinthw / costhw) * mix_sl322) + ( &
      (conjg (mn_43) * mass(15) * mix_sl321) / (mass(24) * cosbe)))))
    g_yuk_n4_sl2_3(2) = (gcc * ((1.0_default * (mn_42 + (1.0_default *  &
      (sinthw / costhw) * mn_41)) * mix_sl321) - ( &
      (mn_43 * mass(15) * mix_sl322) / (mass(24) * cosbe))))
    g_yuk_n4_sl2_3_c(1) = conjg (g_yuk_n4_sl2_3(2))
    g_yuk_n4_sl2_3_c(2) = conjg (g_yuk_n4_sl2_3(1))
    g_yuk_n4_su2_3(1) = ( - (gcc * ((2.0_default * ( -  &
      q_up) * conjg (mn_41) *  &
      (sinthw / costhw) * mix_su322) + ((conjg (mn_44) * mass(6) * mix_su321) /  &
      (mass(24) * sinbe)))))
    g_yuk_n4_su2_3(2) = (gcc * (((-1.0_default) * (mn_42 + ( &
      (1.0_default / 3.0_default) *  &
      (sinthw / costhw) * mn_41)) * mix_su321) - ( &
      (mn_44 * mass(6) * mix_su322) / (mass(24) * sinbe))))
    g_yuk_n4_su2_3_c(1) = conjg (g_yuk_n4_su2_3(2))
    g_yuk_n4_su2_3_c(2) = conjg (g_yuk_n4_su2_3(1))
    g_yuk_n4_sd2_3(1) = ( - (gcc * ((2.0_default * ( - ( -  &
      (1.0_default / 3.0_default))) * conjg (mn_41) *  &
      (sinthw / costhw) * mix_sd322) + ((conjg (mn_43) * mass(5) * mix_sd321) /  &
      (mass(24) * cosbe)))))
    g_yuk_n4_sd2_3(2) = (gcc * ((1.0_default * (mn_42 + (( -  &
      (1.0_default / 3.0_default)) *  &
      (sinthw / costhw) * mn_41)) * mix_sd321) - ( &
      (mn_43 * mass(5) * mix_sd322) / (mass(24) * cosbe))))
    !!! For the adjoint color flow method these constants have to be
    !!! divided by a factor of sqrt(2). 
    g_yuk_n4_sd2_3_c(1) = conjg (g_yuk_n4_sd2_3(2))
    g_yuk_n4_sd2_3_c(2) = conjg (g_yuk_n4_sd2_3(1))
    !!! For the diagram-wise color calculation this has not to be
    !!! divided by an additional factor of sqrt(2) 
    g_yuk_gsu1_3(1) = ( - (mix_su312 * (gs / sqrt (2.0_default))))
    g_yuk_gsu1_3(2) = (mix_su311 * (gs / sqrt (2.0_default)))
    g_yuk_gsu1_3_c(1) = conjg (g_yuk_gsu1_3(2))
    g_yuk_gsu1_3_c(2) = conjg (g_yuk_gsu1_3(1))
    g_yuk_gsd1_3(1) = ( - (mix_sd312 * (gs / sqrt (2.0_default))))
    g_yuk_gsd1_3(2) = (mix_sd311 * (gs / sqrt (2.0_default)))
    g_yuk_gsd1_3_c(1) = conjg (g_yuk_gsd1_3(2))
    g_yuk_gsd1_3_c(2) = conjg (g_yuk_gsd1_3(1))
    g_yuk_gsu2_3(1) = ( - (mix_su322 * (gs / sqrt (2.0_default))))
    g_yuk_gsu2_3(2) = (mix_su321 * (gs / sqrt (2.0_default)))
    g_yuk_gsu2_3_c(1) = conjg (g_yuk_gsu2_3(2))
    g_yuk_gsu2_3_c(2) = conjg (g_yuk_gsu2_3(1))
    g_yuk_gsd2_3(1) = ( - (mix_sd322 * (gs / sqrt (2.0_default))))
    g_yuk_gsd2_3(2) = (mix_sd321 * (gs / sqrt (2.0_default)))
    g_yuk_gsd2_3_c(1) = conjg (g_yuk_gsd2_3(2))
    g_yuk_gsd2_3_c(2) = conjg (g_yuk_gsd2_3(1))
  end subroutine setup_parameters16
  subroutine setup_parameters17
   !!! This subroutine contains the gravitino parameters
   !!! ggrav = - m_planck / 4.0_default 
   ggrav = - imago / 4.0_default / m_planck
   ggrch1(1) = ggrav * conjg(mv_11)
   ggrch1(2) = ggrav * mu_11
   ggrch2(1) = ggrav * conjg(mv_21)
   ggrch2(2) = ggrav * mu_21
   ggrch1c = conjg(ggrch1)
   ggrch2c = conjg(ggrch2)
   ggrzneu1(1) = ggrav * (costhw * conjg(mn_12) - sinthw * conjg(mn_11)) 
   ggrzneu1(2) = conjg(ggrzneu1(1))
   ggrzneu2(1) = ggrav * (costhw * conjg(mn_22) - sinthw * conjg(mn_21)) 
   ggrzneu2(2) = conjg(ggrzneu2(1))
   ggrzneu3(1) = ggrav * (costhw * conjg(mn_32) - sinthw * conjg(mn_31)) 
   ggrzneu3(2) = conjg(ggrzneu3(1))
   ggrzneu4(1) = ggrav * (costhw * conjg(mn_42) - sinthw * conjg(mn_41)) 
   ggrzneu4(2) = conjg(ggrzneu4(1))
   ggraneu1(1) = ggrav * (sinthw * conjg(mn_12) + costhw * conjg(mn_11)) 
   ggraneu1(2) = conjg(ggraneu1(1))
   ggraneu2(1) = ggrav * (sinthw * conjg(mn_22) + costhw * conjg(mn_21)) 
   ggraneu2(2) = conjg(ggraneu2(1))
   ggraneu3(1) = ggrav * (sinthw * conjg(mn_32) + costhw * conjg(mn_31)) 
   ggraneu3(2) = conjg(ggraneu3(1))
   ggraneu4(1) = ggrav * (sinthw * conjg(mn_42) + costhw * conjg(mn_41)) 
   ggraneu4(2) = conjg(ggraneu4(1))
   ggr4neu1(1) = - ggrav * g * conjg(mn_12)
   ggr4neu1(1) = conjg(ggr4neu1(1))
   ggr4neu2(1) = - ggrav * g * conjg(mn_22)
   ggr4neu2(1) = conjg(ggr4neu2(1))
   ggr4neu3(1) = - ggrav * g * conjg(mn_32)
   ggr4neu3(1) = conjg(ggr4neu3(1))
   ggr4neu4(1) = - ggrav * g * conjg(mn_42)
   ggr4neu4(1) = conjg(ggr4neu4(1))
   ggr4ach1(1) = - ggrch1(1) * e
   ggr4ach1(2) = - ggrch1(2) * e
   ggr4ach2(1) = - ggrch2(1) * e
   ggr4ach2(2) = - ggrch2(2) * e
   ggr4ach1c = conjg(ggr4ach1)
   ggr4ach2c = conjg(ggr4ach2)
   ggr4zch1(1) = - ggrch1(1) * g * costhw
   ggr4zch1(2) = - ggrch1(2) * g * costhw
   ggr4zch2(1) = - ggrch2(1) * g * costhw
   ggr4zch2(2) = - ggrch2(2) * g * costhw
   ggr4zch1c = conjg(ggr4zch1)
   ggr4zch2c = conjg(ggr4zch2)
   ggravn = - imago / sqrt(2.0_default) / m_planck
   !!! ggravn = - 1.0_default / sqrt(2.0_default) / m_planck
   ggravu11(1) = ggravn * conjg(mix_su111)
   ggravu11(2) = ggravn * conjg(mix_su112)
   ggravu12(1) = ggravn * conjg(mix_su121)
   ggravu12(2) = ggravn * conjg(mix_su122)
   ggravu11c = conjg(ggravu11)
   ggravu12c = conjg(ggravu12)
   ggravu21(1) = ggravn * conjg(mix_su211)
   ggravu11(2) = ggravn * conjg(mix_su112)
   ggravu12(1) = ggravn * conjg(mix_su121)
   ggravu12(2) = ggravn * conjg(mix_su122)
   ggravu11c = conjg(ggravu11)
   ggravu12c = conjg(ggravu12)
   ggravu21(1) = ggravn * conjg(mix_su211)
   ggravu21(2) = ggravn * conjg(mix_su212)
   ggravu22(1) = ggravn * conjg(mix_su221)
   ggravu22(2) = ggravn * conjg(mix_su222)
   ggravu21c = conjg(ggravu21)
   ggravu22c = conjg(ggravu22)
   ggravu31(1) = ggravn * conjg(mix_su311)
   ggravu31(2) = ggravn * conjg(mix_su312)
   ggravu32(1) = ggravn * conjg(mix_su321)
   ggravu32(2) = ggravn * conjg(mix_su322)
   ggravu31c = conjg(ggravu31)
   ggravu32c = conjg(ggravu32)
   ggravd11(1) = ggravn * conjg(mix_sd111)
   ggravd11(2) = ggravn * conjg(mix_sd112)
   ggravd12(1) = ggravn * conjg(mix_sd121)
   ggravd12(2) = ggravn * conjg(mix_sd122)
   ggravd11c = conjg(ggravd11)
   ggravd12c = conjg(ggravd12)
   ggravd21(1) = ggravn * conjg(mix_sd211)
   ggravd21(2) = ggravn * conjg(mix_sd212)
   ggravd22(1) = ggravn * conjg(mix_sd221)
   ggravd22(2) = ggravn * conjg(mix_sd222)
   ggravd21c = conjg(ggravd21)
   ggravd22c = conjg(ggravd22)
   ggravd31(1) = ggravn * conjg(mix_sd311)
   ggravd31(2) = ggravn * conjg(mix_sd312)
   ggravd32(1) = ggravn * conjg(mix_sd321)
   ggravd32(2) = ggravn * conjg(mix_sd322)
   ggravd31c = conjg(ggravd31)
   ggravd32c = conjg(ggravd32)
   ggravl11(1) = ggravn * conjg(mix_sl111)
   ggravl11(2) = ggravn * conjg(mix_sl112)
   ggravl12(1) = ggravn * conjg(mix_sl121)
   ggravl12(2) = ggravn * conjg(mix_sl122)
   ggravl11c = conjg(ggravl11)
   ggravl12c = conjg(ggravl12)
   ggravl21(1) = ggravn * conjg(mix_sl211)
   ggravl21(2) = ggravn * conjg(mix_sl212)
   ggravl22(1) = ggravn * conjg(mix_sl221)
   ggravl22(2) = ggravn * conjg(mix_sl222)
   ggravl21c = conjg(ggravl21)
   ggravl22c = conjg(ggravl22)
   ggravl31(1) = ggravn * conjg(mix_sl311)
   ggravl31(2) = ggravn * conjg(mix_sl312)
   ggravl32(1) = ggravn * conjg(mix_sl321)
   ggravl32(2) = ggravn * conjg(mix_sl322)
   ggravl31c = conjg(ggravl31)
   ggravl32c = conjg(ggravl32)
   ggrhch1(1) = ggravn * (cosbe * conjg(mv_12) + sinbe * mu_12)  
   ggrhch1(2) = conjg(ggrhch1(1))
   ggrhch2(1) = ggravn * (cosbe * conjg(mv_22) + sinbe * mu_22)  
   ggrhch2(2) = conjg(ggrhch2(1))
   ggrhch1c = conjg(ggrhch1)
   ggrhch2c = conjg(ggrhch2)
   ggrh1neu1(1) = ggravn * (cosal * conjg(mn_14) - sinal * conjg(mn_13))/sqrt(2.0_default)
   ggrh1neu1(2) = conjg(ggrh1neu1(1))
   ggrh1neu2(1) = ggravn * (cosal * conjg(mn_24) - sinal * conjg(mn_23))/sqrt(2.0_default)
   ggrh1neu2(2) = conjg(ggrh1neu2(1))
   ggrh1neu3(1) = ggravn * (cosal * conjg(mn_34) - sinal * conjg(mn_33))/sqrt(2.0_default)
   ggrh1neu3(2) = conjg(ggrh1neu3(1))
   ggrh1neu4(1) = ggravn * (cosal * conjg(mn_44) - sinal * conjg(mn_43))/sqrt(2.0_default)
   ggrh1neu4(2) = conjg(ggrh1neu4(1))
   ggrh2neu1(1) = ggravn * (sinal * conjg(mn_14) + cosal * conjg(mn_13))/sqrt(2.0_default)
   ggrh2neu1(2) = conjg(ggrh2neu1(1))
   ggrh2neu2(1) = ggravn * (sinal * conjg(mn_24) + cosal * conjg(mn_23))/sqrt(2.0_default)
   ggrh2neu2(2) = conjg(ggrh2neu2(1))
   ggrh2neu3(1) = ggravn * (sinal * conjg(mn_34) + cosal * conjg(mn_33))/sqrt(2.0_default)
   ggrh2neu3(2) = conjg(ggrh2neu3(1))
   ggrh2neu4(1) = ggravn * (sinal * conjg(mn_44) + cosal * conjg(mn_43))/sqrt(2.0_default)
   ggrh2neu4(2) = conjg(ggrh2neu4(1))
   ggrh3neu1(1) = - imago * ggravn * (cosbe * conjg(mn_14) + sinbe * conjg(mn_13))/sqrt(2.0_default)
   ggrh3neu1(2) = conjg(ggrh3neu1(1))
   ggrh3neu2(1) = - imago * ggravn * (cosbe * conjg(mn_24) + sinbe * conjg(mn_23))/sqrt(2.0_default)
   ggrh3neu2(2) = conjg(ggrh3neu2(1))
   ggrh3neu3(1) = - imago * ggravn * (cosbe * conjg(mn_34) + sinbe * conjg(mn_33))/sqrt(2.0_default)
   ggrh3neu3(2) = conjg(ggrh3neu3(1))
   ggrh3neu4(1) = - imago * ggravn * (cosbe * conjg(mn_44) + sinbe * conjg(mn_43))/sqrt(2.0_default)
   ggrh3neu4(2) = conjg(ggrh3neu4(1))
   ggr4asl11(1) = - ggravn * qlep * conjg(mix_sl111)
   ggr4asl11(2) = ggravn * qlep * conjg(mix_sl112)
   ggr4asl12(1) = - ggravn * qlep * conjg(mix_sl121)
   ggr4asl12(2) = ggravn * qlep * conjg(mix_sl122)
   ggr4asl21(1) = - ggravn * qlep * conjg(mix_sl211)
   ggr4asl21(2) = ggravn * qlep * conjg(mix_sl212)
   ggr4asl22(1) = - ggravn * qlep * conjg(mix_sl221)
   ggr4asl22(2) = ggravn * qlep * conjg(mix_sl222)
   ggr4asl31(1) = - ggravn * qlep * conjg(mix_sl311)
   ggr4asl31(2) = ggravn * qlep * conjg(mix_sl312)
   ggr4asl32(1) = - ggravn * qlep * conjg(mix_sl321)
   ggr4asl32(2) = ggravn * qlep * conjg(mix_sl322)
   ggr4asl11c = - conjg(ggr4asl11) 
   ggr4asl12c = - conjg(ggr4asl12)
   ggr4asl21c = - conjg(ggr4asl21)
   ggr4asl22c = - conjg(ggr4asl22)
   ggr4asl31c = - conjg(ggr4asl31)
   ggr4asl32c = - conjg(ggr4asl32)
   ggr4asu11(1) = - ggravn * qup * conjg(mix_su111)
   ggr4asu11(2) = ggravn * qup * conjg(mix_su112)
   ggr4asu12(1) = - ggravn * qup * conjg(mix_su121)
   ggr4asu12(2) = ggravn * qup * conjg(mix_su122)
   ggr4asu21(1) = - ggravn * qup * conjg(mix_su211)
   ggr4asu21(2) = ggravn * qup * conjg(mix_su212)
   ggr4asu22(1) = - ggravn * qup * conjg(mix_su221)
   ggr4asu22(2) = ggravn * qup * conjg(mix_su222)
   ggr4asu31(1) = - ggravn * qup * conjg(mix_su311)
   ggr4asu31(2) = ggravn * qup * conjg(mix_su312)
   ggr4asu32(1) = - ggravn * qup * conjg(mix_su321)
   ggr4asu32(2) = ggravn * qup * conjg(mix_su322)
   ggr4asu11c = - conjg(ggr4asu11) 
   ggr4asu12c = - conjg(ggr4asu12)
   ggr4asu21c = - conjg(ggr4asu21)
   ggr4asu22c = - conjg(ggr4asu22)
   ggr4asu31c = - conjg(ggr4asu31)
   ggr4asu32c = - conjg(ggr4asu32)
   ggr4asd11(1) = - ggravn * qdwn * conjg(mix_sd111)
   ggr4asd11(2) = ggravn * qdwn * conjg(mix_sd112)
   ggr4asd12(1) = - ggravn * qdwn * conjg(mix_sd121)
   ggr4asd12(2) = ggravn * qdwn * conjg(mix_sd122)
   ggr4asd21(1) = - ggravn * qdwn * conjg(mix_sd211)
   ggr4asd21(2) = ggravn * qdwn * conjg(mix_sd212)
   ggr4asd22(1) = - ggravn * qdwn * conjg(mix_sd221)
   ggr4asd22(2) = ggravn * qdwn * conjg(mix_sd222)
   ggr4asd31(1) = - ggravn * qdwn * conjg(mix_sd311)
   ggr4asd31(2) = ggravn * qdwn * conjg(mix_sd312)
   ggr4asd32(1) = - ggravn * qdwn * conjg(mix_sd321)
   ggr4asd32(2) = ggravn * qdwn * conjg(mix_sd322)
   ggr4asd11c = - conjg(ggr4asd11) 
   ggr4asd12c = - conjg(ggr4asd12)
   ggr4asd21c = - conjg(ggr4asd21)
   ggr4asd22c = - conjg(ggr4asd22)
   ggr4asd31c = - conjg(ggr4asd31)
   ggr4asd32c = - conjg(ggr4asd32)
   gzll = gz * (- 1.0/2.0_default  - sin2thw * qelep)
   gzlr = - gz * sin2thw * qelep
   gzul = gz * (1.0/2.0_default - sin2thw * qeup)
   gzur = - gz * sin2thw * qeup
   gzdl = gz * (- 1.0/2.0_default  - sin2thw * qedwn)
   gzdr = - gz * sin2thw * qedwn   
   !!! Here the sign structure is different because of the sign
   !!! in qlep, qdwn, qup
   ggr4zsl11(1) = ggravn * gzll * conjg(mix_sl111) 
   ggr4zsl11(2) = - ggravn * gzlr * conjg(mix_sl112)
   ggr4zsl12(1) = ggravn * gzll * conjg(mix_sl121)
   ggr4zsl12(2) = - ggravn * gzlr * conjg(mix_sl122)
   ggr4zsl21(1) = ggravn * gzll * conjg(mix_sl211)
   ggr4zsl21(2) = - ggravn * gzlr * conjg(mix_sl212)
   ggr4zsl22(1) = ggravn * gzll * conjg(mix_sl221)
   ggr4zsl22(2) = - ggravn * gzlr * conjg(mix_sl222)
   ggr4zsl31(1) = ggravn * gzll * conjg(mix_sl311)
   ggr4zsl31(2) = - ggravn * gzlr * conjg(mix_sl312)
   ggr4zsl32(1) = ggravn * gzll * conjg(mix_sl321)
   ggr4zsl32(2) = - ggravn * gzlr * conjg(mix_sl322)
   ggr4zsl11c = - conjg(ggr4zsl11)  
   ggr4zsl12c = - conjg(ggr4zsl12)
   ggr4zsl21c = - conjg(ggr4zsl21)
   ggr4zsl22c = - conjg(ggr4zsl22)
   ggr4zsl31c = - conjg(ggr4zsl31)
   ggr4zsl32c = - conjg(ggr4zsl32)
   ggr4zsu11(1) = ggravn * gzul * conjg(mix_su111) 
   ggr4zsu11(2) = - ggravn * gzur * conjg(mix_su112)
   ggr4zsu12(1) = ggravn * gzul * conjg(mix_su121)
   ggr4zsu12(2) = - ggravn * gzur * conjg(mix_su122)
   ggr4zsu21(1) = ggravn * gzul * conjg(mix_su211)
   ggr4zsu21(2) = - ggravn * gzur * conjg(mix_su212)
   ggr4zsu22(1) = ggravn * gzul * conjg(mix_su221)
   ggr4zsu22(2) = - ggravn * gzur * conjg(mix_su222)
   ggr4zsu31(1) = ggravn * gzul * conjg(mix_su311)
   ggr4zsu31(2) = - ggravn * gzur * conjg(mix_su312)
   ggr4zsu32(1) = ggravn * gzul * conjg(mix_su321)
   ggr4zsu32(2) = - ggravn * gzur * conjg(mix_su322)
   ggr4zsu11c = - conjg(ggr4zsu11)  
   ggr4zsu12c = - conjg(ggr4zsu12)
   ggr4zsu21c = - conjg(ggr4zsu21)
   ggr4zsu22c = - conjg(ggr4zsu22)
   ggr4zsu31c = - conjg(ggr4zsu31)
   ggr4zsu32c = - conjg(ggr4zsu32)
   ggr4zsd11(1) = ggravn * gzdl * conjg(mix_sd111) 
   ggr4zsd11(2) = - ggravn * gzdr * conjg(mix_sd112)
   ggr4zsd12(1) = ggravn * gzdl * conjg(mix_sd121)
   ggr4zsd12(2) = - ggravn * gzdr * conjg(mix_sd122)
   ggr4zsd21(1) = ggravn * gzdl * conjg(mix_sd211)
   ggr4zsd21(2) = - ggravn * gzdr * conjg(mix_sd212)
   ggr4zsd22(1) = ggravn * gzdl * conjg(mix_sd221)
   ggr4zsd22(2) = - ggravn * gzdr * conjg(mix_sd222)
   ggr4zsd31(1) = ggravn * gzdl * conjg(mix_sd311)
   ggr4zsd31(2) = - ggravn * gzdr * conjg(mix_sd312)
   ggr4zsd32(1) = ggravn * gzdl * conjg(mix_sd321)
   ggr4zsd32(2) = - ggravn * gzdr * conjg(mix_sd322)
   ggr4zsd11c = - conjg(ggr4zsd11)  
   ggr4zsd12c = - conjg(ggr4zsd12)
   ggr4zsd21c = - conjg(ggr4zsd21)
   ggr4zsd22c = - conjg(ggr4zsd22)
   ggr4zsd31c = - conjg(ggr4zsd31)
   ggr4zsd32c = - conjg(ggr4zsd32)
   ggr4zsn  = ggravn * gz / 2.0_default
   ggr4zsnc = - conjg(ggr4zsn)
   ggr4wsn  = ggravn * g / sqrt(2.0_default)
   ggr4wsnc = - conjg(ggr4wsn)
   ggr4wsl11  = ggravn * g / sqrt(2.0_default) * conjg(mix_sl111) 
   ggr4wsl12  = ggravn * g / sqrt(2.0_default) * conjg(mix_sl121) 
   ggr4wsl21  = ggravn * g / sqrt(2.0_default) * conjg(mix_sl211) 
   ggr4wsl22  = ggravn * g / sqrt(2.0_default) * conjg(mix_sl221) 
   ggr4wsl31  = ggravn * g / sqrt(2.0_default) * conjg(mix_sl311) 
   ggr4wsl32  = ggravn * g / sqrt(2.0_default) * conjg(mix_sl321) 
   ggr4wsl11c = - conjg(ggr4wsl11)
   ggr4wsl12c = - conjg(ggr4wsl12)
   ggr4wsl21c = - conjg(ggr4wsl21)
   ggr4wsl22c = - conjg(ggr4wsl22)
   ggr4wsl31c = - conjg(ggr4wsl31)
   ggr4wsl32c = - conjg(ggr4wsl32)
   ggr4wsu11  = ggravn * g / sqrt(2.0_default) * conjg(mix_su111) 
   ggr4wsu12  = ggravn * g / sqrt(2.0_default) * conjg(mix_su121) 
   ggr4wsu21  = ggravn * g / sqrt(2.0_default) * conjg(mix_su211) 
   ggr4wsu22  = ggravn * g / sqrt(2.0_default) * conjg(mix_su221) 
   ggr4wsu31  = ggravn * g / sqrt(2.0_default) * conjg(mix_su311) 
   ggr4wsu32  = ggravn * g / sqrt(2.0_default) * conjg(mix_su321) 
   ggr4wsu11c = - conjg(ggr4wsu11)
   ggr4wsu12c = - conjg(ggr4wsu12)
   ggr4wsu21c = - conjg(ggr4wsu21)
   ggr4wsu22c = - conjg(ggr4wsu22)
   ggr4wsu31c = - conjg(ggr4wsu31)
   ggr4wsu32c = - conjg(ggr4wsu32)
   ggr4wsd11  = ggravn * g / sqrt(2.0_default) * conjg(mix_sd111) 
   ggr4wsd12  = ggravn * g / sqrt(2.0_default) * conjg(mix_sd121) 
   ggr4wsd21  = ggravn * g / sqrt(2.0_default) * conjg(mix_sd211) 
   ggr4wsd22  = ggravn * g / sqrt(2.0_default) * conjg(mix_sd221) 
   ggr4wsd31  = ggravn * g / sqrt(2.0_default) * conjg(mix_sd311) 
   ggr4wsd32  = ggravn * g / sqrt(2.0_default) * conjg(mix_sd321) 
   ggr4wsd11c = - conjg(ggr4wsd11)
   ggr4wsd12c = - conjg(ggr4wsd12)
   ggr4wsd21c = - conjg(ggr4wsd21)
   ggr4wsd22c = - conjg(ggr4wsd22)
   ggr4wsd31c = - conjg(ggr4wsd31)
   ggr4wsd32c = - conjg(ggr4wsd32)
   ggr4glsu11(1) = ggravn * gs * conjg(mix_su111)
   ggr4glsu11(2) = - ggravn * gs * conjg(mix_su112)
   ggr4glsu12(1) = ggravn * gs * conjg(mix_su121)
   ggr4glsu12(2) = - ggravn * gs * conjg(mix_su122)
   ggr4glsu21(1) = ggravn * gs * conjg(mix_su211)
   ggr4glsu21(2) = - ggravn * gs * conjg(mix_su212)
   ggr4glsu22(1) = ggravn * gs * conjg(mix_su221)
   ggr4glsu22(2) = - ggravn * gs * conjg(mix_su222)
   ggr4glsu31(1) = ggravn * gs * conjg(mix_su311)
   ggr4glsu31(2) = - ggravn * gs * conjg(mix_su312)
   ggr4glsu32(1) = ggravn * gs * conjg(mix_su321)
   ggr4glsu32(2) = - ggravn * gs * conjg(mix_su322)
   ggr4glsu11c = conjg(ggr4glsu11) 
   ggr4glsu12c = conjg(ggr4glsu12)
   ggr4glsu21c = conjg(ggr4glsu21)
   ggr4glsu22c = conjg(ggr4glsu22)
   ggr4glsu31c = conjg(ggr4glsu31)
   ggr4glsu32c = conjg(ggr4glsu32)
   ggr4glsd11(1) = ggravn * gs * conjg(mix_sd111)
   ggr4glsd11(2) = - ggravn * gs * conjg(mix_sd112)
   ggr4glsd12(1) = ggravn * gs * conjg(mix_sd121)
   ggr4glsd12(2) = - ggravn * gs * conjg(mix_sd122)
   ggr4glsd21(1) = ggravn * gs * conjg(mix_sd211)
   ggr4glsd21(2) = - ggravn * gs * conjg(mix_sd212)
   ggr4glsd22(1) = ggravn * gs * conjg(mix_sd221)
   ggr4glsd22(2) = - ggravn * gs * conjg(mix_sd222)
   ggr4glsd31(1) = ggravn * gs * conjg(mix_sd311)
   ggr4glsd31(2) = - ggravn * gs * conjg(mix_sd312)
   ggr4glsd32(1) = ggravn * gs * conjg(mix_sd321)
   ggr4glsd32(2) = - ggravn * gs * conjg(mix_sd322)
   ggr4glsd11c = conjg(ggr4glsd11) 
   ggr4glsd12c = conjg(ggr4glsd12)
   ggr4glsd21c = conjg(ggr4glsd21)
   ggr4glsd22c = conjg(ggr4glsd22)
   ggr4glsd31c = conjg(ggr4glsd31)
   ggr4glsd32c = conjg(ggr4glsd32)
   ggr4zh1_1   = gz * ggrh1neu1 / 2.0_default
   ggr4zh1_2   = gz * ggrh1neu2 / 2.0_default
   ggr4zh1_3   = gz * ggrh1neu3 / 2.0_default
   ggr4zh1_4   = gz * ggrh1neu4 / 2.0_default
   ggr4zh2_1   = gz * ggrh2neu1 / 2.0_default
   ggr4zh2_2   = gz * ggrh2neu2 / 2.0_default
   ggr4zh2_3   = gz * ggrh2neu3 / 2.0_default
   ggr4zh2_4   = gz * ggrh2neu4 / 2.0_default
   ggr4zh3_1   = gz * ggrh3neu1 / 2.0_default
   ggr4zh3_2   = gz * ggrh3neu2 / 2.0_default
   ggr4zh3_3   = gz * ggrh3neu3 / 2.0_default
   ggr4zh3_4   = gz * ggrh3neu4 / 2.0_default
   ggr4wh_1    = g * ggrhch1 / sqrt(2.0_default)
   ggr4wh_2    = g * ggrhch2 / sqrt(2.0_default)
   ggr4wh_1c   = g * ggrhch1c / sqrt(2.0_default)
   ggr4wh_2c   = g * ggrhch2c / sqrt(2.0_default)
   ggr4ha1(1)  = e * ggrhch1(1)
   ggr4ha1(2)  = - e * ggrhch1(2)
   ggr4ha2(1)  = e * ggrhch2(1)
   ggr4ha2(2)  = - e * ggrhch2(2)
   ggr4ha1c(1) = - e * ggrhch1c(1)
   ggr4ha1c(2) = e * ggrhch1c(2)
   ggr4ha2c(1) = - e * ggrhch2c(1)
   ggr4ha2c(2) = e * ggrhch2c(2)
   ggr4hz1     = gz * ggrhch1 / 2.0_default
   ggr4hz2     = gz * ggrhch2 / 2.0_default
   ggr4hz1c    = gz * ggrhch1c / 2.0_default
   ggr4hz2c    = gz * ggrhch2c / 2.0_default
  end subroutine setup_parameters17
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default * PI * alpha_s)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs
    gssq = (gs / sqrt (2.0_default))
    g_yuk_gsu1_3(1) = ( - (mix_su312 * (gs / sqrt (2.0_default))))
    g_yuk_gsu1_3(2) = (mix_su311 * (gs / sqrt (2.0_default)))
    g_yuk_gsu1_3_c(1) = conjg (g_yuk_gsu1_3(2))
    g_yuk_gsu1_3_c(2) = conjg (g_yuk_gsu1_3(1))
    g_yuk_gsd1_3(1) = ( - (mix_sd312 * (gs / sqrt (2.0_default))))
    g_yuk_gsd1_3(2) = (mix_sd311 * (gs / sqrt (2.0_default)))
    g_yuk_gsd1_3_c(1) = conjg (g_yuk_gsd1_3(2))
    g_yuk_gsd1_3_c(2) = conjg (g_yuk_gsd1_3(1))
    g_yuk_gsu2_3(1) = ( - (mix_su322 * (gs / sqrt (2.0_default))))
    g_yuk_gsu2_3(2) = (mix_su321 * (gs / sqrt (2.0_default)))
    g_yuk_gsu2_3_c(1) = conjg (g_yuk_gsu2_3(2))
    g_yuk_gsu2_3_c(2) = conjg (g_yuk_gsu2_3(1))
    g_yuk_gsd2_3(1) = ( - (mix_sd322 * (gs / sqrt (2.0_default))))
    g_yuk_gsd2_3(2) = (mix_sd321 * (gs / sqrt (2.0_default)))
    g_yuk_gsd2_3_c(1) = conjg (g_yuk_gsd2_3(2))
    g_yuk_gsd2_3_c(2) = conjg (g_yuk_gsd2_3(1))
    gglglsqsq = (gs**2)
    gglpsqsq = 2.0_default * e * gs / 3.0_default
    gglsu1su1_1 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su111 * conjg (mix_su111))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu2su2_1 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su121 * conjg (mix_su121))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu1su2_1 = (gz * gs * (1.0_default / 2.0_default) * mix_su111 *  &
      conjg (mix_su121))
    gglsu2su1_1 = (gz * gs * (1.0_default / 2.0_default) * mix_su121 *  &
      conjg (mix_su111))
    gglsd1sd1_1 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd111 * conjg (mix_sd111))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd2sd2_1 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd121 * conjg (mix_sd121))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd1sd2_1 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd111 * conjg (mix_sd121)))
    gglsd2sd1_1 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd121 * conjg (mix_sd111)))
    gglsu1su1_2 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su211 * conjg (mix_su211))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu2su2_2 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su221 * conjg (mix_su221))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu1su2_2 = (gz * gs * (1.0_default / 2.0_default) * mix_su211 *  &
      conjg (mix_su221))
    gglsu2su1_2 = (gz * gs * (1.0_default / 2.0_default) * mix_su221 *  &
      conjg (mix_su211))
    gglsd1sd1_2 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd211 * conjg (mix_sd211))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd2sd2_2 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd221 * conjg (mix_sd221))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd1sd2_2 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd211 * conjg (mix_sd221)))
    gglsd2sd1_2 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd221 * conjg (mix_sd211)))
    gglsu1su1_3 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su311 * conjg (mix_su311))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu2su2_3 = (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_su321 * conjg (mix_su321))) - (sin2thw *  &
      (2.0_default / 3.0_default))))
    gglsu1su2_3 = (gz * gs * (1.0_default / 2.0_default) * mix_su311 *  &
      conjg (mix_su321))
    gglsu2su1_3 = (gz * gs * (1.0_default / 2.0_default) * mix_su321 *  &
      conjg (mix_su311))
    gglsd1sd1_3 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd311 * conjg (mix_sd311))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd2sd2_3 = ( - (gz * gs * (((1.0_default / 2.0_default) *  &
      (mix_sd321 * conjg (mix_sd321))) - (sin2thw *  &
      (1.0_default / 3.0_default)))))
    gglsd1sd2_3 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd311 * conjg (mix_sd321)))
    gglsd2sd1_3 = ( - (gz * gs *  &
      (1.0_default / 2.0_default) * mix_sd321 * conjg (mix_sd311)))
    gglwsu1sd1_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su111) * mix_sd111)
    gglwsu2sd2_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su121) * mix_sd121)
    gglwsu1sd2_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su111) * mix_sd121)
    gglwsu2sd1_1_1 = (g * gs * sqrt (2.0_default) * vckm_11 *  &
      conjg (mix_su121) * mix_sd111)
    gglwsu1sd1_1_1_c = conjg (gglwsu1sd1_1_1) 
    gglwsu2sd2_1_1_c = conjg (gglwsu2sd2_1_1)
    gglwsu1sd2_1_1_c = conjg (gglwsu1sd2_1_1)
    gglwsu2sd1_1_1_c = conjg (gglwsu2sd1_1_1)
    gglwsu1sd1_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su111) * mix_sd211)
    gglwsu2sd2_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su121) * mix_sd221)
    gglwsu1sd2_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su111) * mix_sd221)
    gglwsu2sd1_1_2 = (g * gs * sqrt (2.0_default) * vckm_12 *  &
      conjg (mix_su121) * mix_sd211)
    gglwsu1sd1_1_2_c = conjg (gglwsu1sd1_1_2)
    gglwsu2sd2_1_2_c = conjg (gglwsu2sd2_1_2)
    gglwsu1sd2_1_2_c = conjg (gglwsu1sd2_1_2)
    gglwsu2sd1_1_2_c = conjg (gglwsu2sd1_1_2)
    gglwsu1sd1_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su111) * mix_sd311)
    gglwsu2sd2_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su121) * mix_sd321)
    gglwsu1sd2_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su111) * mix_sd321)
    gglwsu2sd1_1_3 = (g * gs * sqrt (2.0_default) * vckm_13 *  &
      conjg (mix_su121) * mix_sd311)
    gglwsu1sd1_1_3_c = conjg (gglwsu1sd1_1_3)  
    gglwsu2sd2_1_3_c = conjg (gglwsu2sd2_1_3)
    gglwsu1sd2_1_3_c = conjg (gglwsu1sd2_1_3)
    gglwsu2sd1_1_3_c = conjg (gglwsu2sd1_1_3)
    gglwsu1sd1_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su211) * mix_sd111)
    gglwsu2sd2_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su221) * mix_sd121)
    gglwsu1sd2_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su211) * mix_sd121)
    gglwsu2sd1_2_1 = (g * gs * sqrt (2.0_default) * vckm_21 *  &
      conjg (mix_su221) * mix_sd111)
    gglwsu1sd1_2_1_c = conjg (gglwsu1sd1_2_1) 
    gglwsu2sd2_2_1_c = conjg (gglwsu2sd2_2_1)
    gglwsu1sd2_2_1_c = conjg (gglwsu1sd2_2_1)
    gglwsu2sd1_2_1_c = conjg (gglwsu2sd1_2_1)
    gglwsu1sd1_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su211) * mix_sd211)
    gglwsu2sd2_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su221) * mix_sd221)
    gglwsu1sd2_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su211) * mix_sd221)
    gglwsu2sd1_2_2 = (g * gs * sqrt (2.0_default) * vckm_22 *  &
      conjg (mix_su221) * mix_sd211)
    gglwsu1sd1_2_2_c = conjg (gglwsu1sd1_2_2) 
    gglwsu2sd2_2_2_c = conjg (gglwsu2sd2_2_2)
    gglwsu1sd2_2_2_c = conjg (gglwsu1sd2_2_2)
    gglwsu2sd1_2_2_c = conjg (gglwsu2sd1_2_2)
    gglwsu1sd1_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su211) * mix_sd311)
    gglwsu2sd2_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su221) * mix_sd321)
    gglwsu1sd2_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su211) * mix_sd321)
    gglwsu2sd1_2_3 = (g * gs * sqrt (2.0_default) * vckm_23 *  &
      conjg (mix_su221) * mix_sd311)
    gglwsu1sd1_2_3_c = conjg (gglwsu1sd1_2_3)
    gglwsu2sd2_2_3_c = conjg (gglwsu2sd2_2_3)
    gglwsu1sd2_2_3_c = conjg (gglwsu1sd2_2_3)
    gglwsu2sd1_2_3_c = conjg (gglwsu2sd1_2_3)
    gglwsu1sd1_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su311) * mix_sd111)
    gglwsu2sd2_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su321) * mix_sd121)
    gglwsu1sd2_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su311) * mix_sd121)
    gglwsu2sd1_3_1 = (g * gs * sqrt (2.0_default) * vckm_31 *  &
      conjg (mix_su321) * mix_sd111)
    gglwsu1sd1_3_1_c = conjg (gglwsu1sd1_3_1)
    gglwsu2sd2_3_1_c = conjg (gglwsu2sd2_3_1)
    gglwsu1sd2_3_1_c = conjg (gglwsu1sd2_3_1)
    gglwsu2sd1_3_1_c = conjg (gglwsu2sd1_3_1)
    gglwsu1sd1_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su311) * mix_sd211)
    gglwsu2sd2_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su321) * mix_sd221)
    gglwsu1sd2_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su311) * mix_sd221)
    gglwsu2sd1_3_2 = (g * gs * sqrt (2.0_default) * vckm_32 *  &
      conjg (mix_su321) * mix_sd211)
    gglwsu1sd1_3_2_c = conjg (gglwsu1sd1_3_2)
    gglwsu2sd2_3_2_c = conjg (gglwsu2sd2_3_2)
    gglwsu1sd2_3_2_c = conjg (gglwsu1sd2_3_2)
    gglwsu2sd1_3_2_c = conjg (gglwsu2sd1_3_2)
    gglwsu1sd1_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su311) * mix_sd311)
    gglwsu2sd2_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su321) * mix_sd321)
    gglwsu1sd2_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su311) * mix_sd321)
    gglwsu2sd1_3_3 = (g * gs * sqrt (2.0_default) * vckm_33 *  &
      conjg (mix_su321) * mix_sd311)
    gglwsu1sd1_3_3_c = conjg (gglwsu1sd1_3_3) 
    gglwsu2sd2_3_3_c = conjg (gglwsu2sd2_3_3)
    gglwsu1sd2_3_3_c = conjg (gglwsu1sd2_3_3)
    gglwsu2sd1_3_3_c = conjg (gglwsu2sd1_3_3)
    ggr4glsu11(1) = ggravn * gs * conjg(mix_su111)
    ggr4glsu11(2) = - ggravn * gs * conjg(mix_su112)
    ggr4glsu12(1) = ggravn * gs * conjg(mix_su121)
    ggr4glsu12(2) = - ggravn * gs * conjg(mix_su122)
    ggr4glsu21(1) = ggravn * gs * conjg(mix_su211)
    ggr4glsu21(2) = - ggravn * gs * conjg(mix_su212)
    ggr4glsu22(1) = ggravn * gs * conjg(mix_su221)
    ggr4glsu22(2) = - ggravn * gs * conjg(mix_su222)
    ggr4glsu31(1) = ggravn * gs * conjg(mix_su311)
    ggr4glsu31(2) = - ggravn * gs * conjg(mix_su312)
    ggr4glsu32(1) = ggravn * gs * conjg(mix_su321)
    ggr4glsu32(2) = - ggravn * gs * conjg(mix_su322)
    ggr4glsu11c = conjg(ggr4glsu11) 
    ggr4glsu12c = conjg(ggr4glsu12)
    ggr4glsu21c = conjg(ggr4glsu21)
    ggr4glsu22c = conjg(ggr4glsu22)
    ggr4glsu31c = conjg(ggr4glsu31)
    ggr4glsu32c = conjg(ggr4glsu32)
    ggr4glsd11(1) = ggravn * gs * conjg(mix_sd111)
    ggr4glsd11(2) = - ggravn * gs * conjg(mix_sd112)
    ggr4glsd12(1) = ggravn * gs * conjg(mix_sd121)
    ggr4glsd12(2) = - ggravn * gs * conjg(mix_sd122)
    ggr4glsd21(1) = ggravn * gs * conjg(mix_sd211)
    ggr4glsd21(2) = - ggravn * gs * conjg(mix_sd212)
    ggr4glsd22(1) = ggravn * gs * conjg(mix_sd221)
    ggr4glsd22(2) = - ggravn * gs * conjg(mix_sd222)
    ggr4glsd31(1) = ggravn * gs * conjg(mix_sd311)
    ggr4glsd31(2) = - ggravn * gs * conjg(mix_sd312)
    ggr4glsd32(1) = ggravn * gs * conjg(mix_sd321)
    ggr4glsd32(2) = - ggravn * gs * conjg(mix_sd322)
    ggr4glsd11c = conjg(ggr4glsd11) 
    ggr4glsd12c = conjg(ggr4glsd12)
    ggr4glsd21c = conjg(ggr4glsd21)
    ggr4glsd22c = conjg(ggr4glsd22)
    ggr4glsd31c = conjg(ggr4glsd31)
    ggr4glsd32c = conjg(ggr4glsd32)
  end subroutine model_update_alpha_s
end module parameters_mssm_grav

