! parameters.PSSSM.f90
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
module parameters_psssm
  use kinds 
  use constants
  implicit none
  private
  public :: import_from_whizard, model_update_alpha_s
  real(kind=default), dimension(107), save, public :: mass = 0, width = 0
  real(kind=default), parameter, public :: GeV = 1.0_default
  real(kind=default), parameter, public :: MeV = GeV / 1000
  real(kind=default), parameter, public :: keV = MeV / 1000
  real(kind=default), parameter, public :: TeV = GeV * 1000
  real(kind=default), save, public :: &
       alpha = 1.0_default / 137.0359895_default, &
       sin2thw = 0.23124_default
  integer, save, public :: &
       sign1 = +1, sign2 = +1, sign3 = +1, sign4 = +1, sign5 = +1, &
       sign6 = +1, sign7 = +1, sign8 = +1, sign9 = +1, sign10 = +1, &
       sign11 = +1
  real(kind=default), save, public :: &
       sigch1 = +1, sigch2 = +1, sigch3 = +1, sigch4 = +1
  complex(kind=default), save, private :: vev
  real(kind=default), public, save :: sind = 0._default, & 
    cosd = 1._default, sinckm12 = 0.223_default, &
    sinckm13 = 0.004_default, sinckm23 = 0.04_default, &
    tana = 30._default, tanb = 30._default, as = 0._default
  real(kind=default), public, save :: & 
    sin2be, &
    cos2be, cosbe, sinbe, costhw, sinthw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: dimh0 = 9, dimA0 = 7, dimNeu = 11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!!!!!!!!!!!!
!!!Higgs Mixing
!!!!!!!!!!!!
  complex(kind=default), dimension(dimh0,dimh0), public, save :: mix_h0 = 0
  complex(kind=default), dimension(dimh0,dimh0), public, save :: mix_A0 = 0  
!!!!!!!!!!!!!!!!
!!!!!Neutralino-Mixing
!!!!!!!!!!!!!!
  complex(kind=default), dimension(dimNeu,dimNeu), public, save :: mix_neu = 0
  integer :: neu1, neu2
!!!!!!!!!!!!!
!!!Chargino Mixing
!!!!!!!!!!!!!
  complex(kind=default), dimension(2,2), public, save :: mix_charU = 0, mix_charV = 0
!!!!!!!!!!!!
!!!CKM-MAtrix
!!!!!!!!!!!!!!
  complex(kind=default), dimension(3,3), public, save :: vckm = 0
  real(kind=default), public, save :: q_lep, q_up, q_down 
  complex(kind=default), public, save :: gcc, qchar, qdwn, qup, qlep, & 
        gz, g, e, gs
  complex(kind=default), save, public :: xia = 1, xi0 = 1, xipm = 1
  complex(kind=default), dimension(2), public, save :: gncdwn
  complex(kind=default), dimension(2), public, save :: gncup
  complex(kind=default), dimension(2), public, save :: gnclep
  complex(kind=default), dimension(2), public, save :: gncneu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Couplings!!!
!! the Lorenz structure is encoded in the model3.ml file in term of
!! S,P,SP,SLR,VLR,.. (rather instructive!). If the type is one of 
!! SP,SLR,VLR,VA the targets.ml file creates an exrta 2-dim 
!! structure which appears as !first! rank of the correspnding 
!! array in this file. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!chargino fermion sfermion (slr,gen,char,sfm)
  complex(kind=default), dimension(2,3,2,1), public,  save :: &
    g_yuk_char_lsn, g_yuk_char_lsn_c 
  complex(kind=default), dimension(2,3,2,2), public,  save :: &
    g_yuk_char_nsl, g_yuk_char_nsl_c 
  complex(kind=default), dimension(2,3,3,2,2), public,  save :: &
    g_yuk_char_dsu, g_yuk_char_dsu_c, g_yuk_char_usd, g_yuk_char_usd_c  
!!!charged higgs to SLeptons (gauge basis)
 complex(kind=default), dimension(3),public, save :: & 
    g_hp_slLsnL, g_hp_slRsnL, g_hp_slLsnL_c, g_hp_slRsnL_c
!!!charged higgs to SLeptons (mass basis)
 complex(kind=default), dimension(3),public, save :: &
    g_hp_sl1sn1, g_hp_sl2sn1, g_hp_sl1sn1_c, g_hp_sl2sn1_c
!!!charged higgs to SQuarks  !!!(gauge basis)
 complex(kind=default), dimension(3,3),public, save :: & 
    g_hp_suLsdL, g_hp_suRsdL, g_hp_suLsdR, g_hp_suRsdR
!!!charged higgs to SQuarks (mass basis)
 complex(kind=default), dimension(3,3),public, save :: &
      g_hp_su1sd1, g_hp_su2sd1, g_hp_su1sd2, g_hp_su2sd2, &
      g_hp_su1sd1_c, g_hp_su2sd1_c, g_hp_su1sd2_c, g_hp_su2sd2_c
!!! Matrix structure coupling neutral Higgs to SQuarks dim(shiggs, gen) 
!!! final mass eigenstate coupling (rep gh1su1su2_3 ...)
  complex(kind=default), dimension(dimh0,3), public, save :: &
       g_h0_su1su1, g_h0_su2su2, g_h0_su1su2, g_h0_su2su1, &
       g_h0_sd1sd1, g_h0_sd2sd2, g_h0_sd1sd2, g_h0_sd2sd1, &
       g_h0_sl1sl1, g_h0_sl2sl2, g_h0_sl1sl2, g_h0_sl2sl1, & 
       g_h0_sn1sn1
  complex(kind=default), dimension(dimA0,3), public, save :: &
       g_A0_su1su1, g_A0_su2su2, g_A0_su1su2, g_A0_su2su1, &
       g_A0_sd1sd1, g_A0_sd2sd2, g_A0_sd1sd2, g_A0_sd2sd1, &
       g_A0_sl1sl1, g_A0_sl2sl2, g_A0_sl1sl2, g_A0_sl2sl1, &  
       g_A0_sn1sn1, g_A0_sn2sn2, g_A0_sn1sn2, g_A0_sn2sn1
!!! Matrix structure coupling neutral Higgs to SQuarks dim(shiggs, gen) 
!!! gauge eigenstate coupling (rep g_h1222susu ...) to be multiplied w/ mix_su
   complex(kind=default), dimension(dimh0,3), public, save :: &
        g_h0_suLsuL, g_h0_suRsuR, g_h0_suLsuR, g_h0_suRsuL, &
        g_h0_sdLsdL, g_h0_sdRsdR, g_h0_sdLsdR, g_h0_sdRsdL, &
        g_h0_slLslL, g_h0_slRslR, g_h0_slLslR, g_h0_slRslL, & 
        g_h0_snLsnL 
  complex(kind=default), dimension(dimA0,3), public, save :: &
       g_A0_suLsuR, g_A0_suRsuL, &
       g_A0_sdLsdR, g_A0_sdRsdL, &
       g_A0_slLslR, g_A0_slRslL 
  integer :: shiggs, phiggs
!!!Matrix structure coupling neutral Higgs to fermions dim (s/phiggs, gen)
  complex(kind=default), dimension(dimh0,3), public, save :: &
       g_yuk_h0_uu, g_yuk_h0_dd, g_yuk_h0_ll   
  complex(kind=default), dimension(dimA0,3), public, save :: &
       g_yuk_A0_uu, g_yuk_A0_dd, g_yuk_A0_ll   
  integer ::gen, gen1, gen2
!!!Charged Higgs to Quarks (SLR, gen1,gen2)
  complex(kind=default), dimension(2,3,3), public, save :: g_yuk_hp_ud, g_yuk_hm_du
!!!Charged Higgs to Leptons (gen)
  complex(kind=default), dimension(3), public, save :: g_yuk_hp_ln  
!!!3- VErtex Higgs Gauge
!!!Z/photon to charged Higgses
  complex(kind=default), public, save ::  g_Ahmhp, g_Zhmhp
!!!Z to (shiggs,phiggs)
  complex(kind=default), dimension(dimh0,dimA0), public, save :: g_Zh0A0 
!!!W to h+ , (shiggs)
  complex(kind=default), dimension(dimh0), public, save :: g_Whph0
!!!W to h+ , (phiggs)
  complex(kind=default), dimension(dimA0), public, save :: g_WhpA0
!!!ZZ to (shiggs)
  complex(kind=default), dimension(dimh0), public, save :: g_ZZh0
!!!WW to (shiggs)
  complex(kind=default), dimension(dimh0), public, save :: g_WWh0
!!!4- Vertex Higgs Gauge
!!!WW (shiggs1, shiggs2) , ZZ (shiggs1, shiggs2)  
  complex(kind=default), dimension(dimh0,dimh0), public, save :: g_WWh0h0 , g_ZZh0h0
!!!WZ hp (shiggs) , WA hp (shiggs)  
  complex(kind=default), dimension(dimh0), public, save :: g_ZWhph0 , g_AWhph0
!!!WW (phiggs1, phiggs2) , ZZ (phiggs1, phiggs2)  
  complex(kind=default), dimension(dimA0,dimA0), public, save :: g_WWA0A0 , g_ZZA0A0
!!!WZ hp (phiggs) , WA hp (phiggs)  
  complex(kind=default), dimension(dimA0), public, save :: g_ZWhpA0 , g_AWhpA0
  complex(kind=default), public, save :: g_ZZhphm, g_ZAhphm, & 
       g_AAhphm, g_WWhphm   
!!!Triple Higgs couplings
  complex(kind=default), dimension(dimh0,dimh0,dimh0) ,public, save :: g_h0h0h0
  complex(kind=default), dimension(dimh0,dimA0,dimA0) ,public, save :: g_h0A0A0
  complex(kind=default), dimension(dimh0) ,public, save :: g_h0hphm
  integer :: shiggs1, shiggs2, shiggs3, phiggs1, phiggs2
!!!Neutral Higgs to Neutralinos (SLR, neu1 , neu2, s/phiggs)
  complex(kind=default), dimension(2,dimNeu,dimNeu,dimh0) ,public, save :: g_neuneuh0  
  complex(kind=default), dimension(2,dimNeu,dimNeu,dimA0) ,public, save :: g_neuneuA0
!!!Neutral Higgs to Charginos (SLR, char1 , char2, s/phiggs)
  complex(kind=default), dimension(2,2,2,dimh0) ,public, save :: g_chchh0  
  complex(kind=default), dimension(2,2,2,dimA0) ,public, save :: g_chchA0
  integer :: ch1, ch2
!!!Chargino, charged Higgs, Neutralino (SLR, neu, char)
  complex(kind=default), dimension(2,dimNeu,2) ,public, save :: &
   g_neuhmchar, g_neuhmchar_c
!!! Neutralino - Fermion -Sfermion
!!!!!(SLR,gen,neu,sfm)
  complex(kind=default), dimension(2,3,dimNeu,2) ,public, save :: &
       g_yuk_neu_lsl, g_yuk_neu_usu, g_yuk_neu_dsd, &
       g_yuk_neu_lsl_c, g_yuk_neu_usu_c, g_yuk_neu_dsd_c, &
       g_yuk_neu_nsn, g_yuk_neu_nsn_c
!!! Gluino - Quark -SQuark
!!!!!(SLR,gen,sfm)
  complex(kind=default), dimension(2,3,2) ,public, save :: &
       g_yuk_gluino_usu, g_yuk_gluino_dsd, g_yuk_gluino_usu_c, g_yuk_gluino_dsd_c
!!! ZZ SFermions  (gen,sfm1,sfm2)
  complex(kind=default), dimension(3,2,2) ,public, save :: &
       g_zz_slsl, g_zz_snsn, g_zz_sdsd, g_zz_susu 
!!! GaGa SFermions  
  complex(kind=default), public, save :: &
       g_AA_slsl, g_AA_snsn, g_AA_sdsd, g_AA_susu 
!!! WW SFermions  (gen,sfm1,sfm2)
  complex(kind=default), dimension(3,2,2) ,public, save :: &
       g_ww_slsl, g_ww_snsn, g_ww_sdsd, g_ww_susu 
!!! ZGa SFermions  (gen,sfm1,sfm2)
  complex(kind=default), dimension(3,2,2) ,public, save :: &
       g_zA_slsl, g_zA_sdsd, g_zA_susu 
!!! W Ga SLeptons  (gen,sfm1)
  complex(kind=default), dimension(3,2) ,public, save :: &
       g_wA_slsn, g_wA_slsn_c 
!!! W Z SLeptons  (gen,sfm1)
  complex(kind=default), dimension(3,2) ,public, save :: &
       g_wz_slsn, g_wz_slsn_c 
!!! W Ga Squarks  (gen1,gen2,sfm1,sfm2)
  complex(kind=default), dimension(3,3,2,2) ,public, save :: &
       g_wA_susd, g_wA_susd_c 
!!! W Z Squarks  (gen1,gen2,sfm1,sfm2)
  complex(kind=default), dimension(3,3,2,2) ,public, save :: &
       g_wz_susd, g_wz_susd_c 
!!! Gluon W Squarks  (gen1,gen2,sfm1,sfm2)
  complex(kind=default), dimension(3,3,2,2) ,public, save :: &
       g_gw_susd, g_gw_susd_c 
!!! Gluon Z Squarks  (gen1,sfm1,sfm2)
  complex(kind=default), dimension(3,2,2) ,public, save :: &
       g_gz_susu, g_gz_sdsd
!!! Gluon A Squarks 
  complex(kind=default), public, save :: g_gA_sqsq
!!! Gluon Gluon Squarks 
  complex(kind=default), public, save :: g_gg_sqsq
!!! Gauge Leptoquarks
  complex(kind=default), public, save :: &
    g_zz_lqlq, g_zA_lqlq, g_AA_lqlq, g_zg_lqlq, g_Ag_lqlq, g_nlqc

  !!!SFermion mixing Matrix Style
  complex(kind=default), dimension(3,2,2), public, save :: &
    mix_sd, mix_su, mix_sl, mix_lq 
  
!!!!!W -slepton -sneutrino
  complex(kind=default), dimension(3,2), public, save :: g_wslsn, g_wslsn_c
!!!!!W -sup -sdown
  complex(kind=default), dimension(3,3,2,2), public, save :: g_wsusd, g_wsusd_c
!!!!!!!!!!!!
!!!Z to SFermions
!!!!!!!!!!!!
  complex(kind=default), dimension(3,1,1), public, save :: g_zsnsn
  complex(kind=default), dimension(3,2,2), public, save :: g_zslsl
  complex(kind=default), dimension(3,2,2), public, save :: g_zsusu
  complex(kind=default), dimension(3,2,2), public, save :: g_zsdsd
  complex(kind=default), dimension(3,2,2), public, save :: g_zlqlq
!!!!! E6 - Yukawa's & Stuff
  complex(kind=default), dimension(3,3,3) ,public, save :: &
  g_yuk_n_d = 0, g_yuk_n_h = 0, g_yuk_d_1 = 0, g_yuk_d_2 = 0, g_yuk_d_c = 0, g_yuk_e = 0, g_yuk_u = 0, g_yuk_nu = 0, g_yuk_d = 0
  complex(kind=default), dimension(2,3,9,3), public, save :: g_yuk_lq_s = 0, g_yuk_lq_p = 0
  complex(kind=default), dimension(2,2,3,3,dimNeu), public, save :: g_lq_neu = 0 
  complex(kind=default), dimension(2,2,3,3), public, save :: g_lq_gg = 0
  complex(kind=default), dimension(2,2,3,3,3), public, save :: g_lq_s = 0, g_lq_p = 0, g_lq_ssd = 0, g_lq_ssta = 0, &
   g_lq_ec_uc = 0,  g_lq_ec_uc_c = 0
  complex(kind=default), dimension(2,2,2,3,3,3), public, save :: g_lq_ssu = 0, g_lq_sst = 0
!!!2 glu 2 lq 
  complex(kind=default), public, save :: g_gg_lqlq = 0

!!! E6- vev-vector 
  complex(kind=default), dimension(9), public, save :: vevs = 0  
  real(kind=default), dimension(3,3), public, save :: delta3 = 0 
!!!NMSSM-parameters
  complex(kind=default), public, save ::  mu, lambda, &         
    A_lambda, k, A_k, r
  complex(kind=default), dimension(3), public, save :: al, au, ad
  complex(kind=default), public, save :: eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, etat, etae
  
  complex(kind=default), public, save :: eidelta, cosckm23, cosckm13, &
    cosckm12,gpzww, gppww, gzzww, gw4, igwww, igzww, iqw, igs, &
    gssq
!!!Charged Current to quarks
  complex(kind=default), dimension(3,3), public, save :: g_ccq, g_ccq_c 
!!!Neutralino, Chargino, W+/-
  complex(kind=default), dimension(2,2,dimNeu), public, save :: g_cwn, g_nwc
!!!Neutral current to Neutralino
  complex(kind=default), dimension(2,dimNeu,dimNeu), public, save :: g_zneuneu
!!!Neutral current to Chargino
  complex(kind=default), dimension(2,2,2), public, save :: g_zchch
   integer :: i,j,t,l,m,n,o,sfm1,sfm2
   complex(kind=default) :: sina, cosa

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(508), intent(in) :: par_array
    integer, intent(in) :: scheme
    type :: parameter_set
      real(default) :: gf
      real(default) :: mz
      real(default) :: mw
      real(default) :: wZ
      real(default) :: ww
      real(default) :: me
      real(default) :: mmu
      real(default) :: mtau
      real(default) :: mc
      real(default) :: ms
      real(default) :: mb
      real(default) :: mtop
      real(default) :: wtop
      real(default) :: alphas
      real(default) :: mtype
      real(default) :: m_zero
      real(default) :: m_half
      real(default) :: a0
      real(default) :: tanb
      real(default) :: sgn_mu
      real(default) :: lambda
      real(default) :: m_mes
      real(default) :: n_5
      real(default) :: c_grav
      real(default) :: m_grav
      real(default) :: yuk_d1
      real(default) :: yuk_d2
      real(default) :: yuk_dc
      real(default) :: yuk_nd
      real(default) :: ae_33
      real(default) :: au_33
      real(default) :: ad_33
      real(default) :: mh0_1
      real(default) :: wh0_1
      real(default) :: mh0_2
      real(default) :: wh0_2
      real(default) :: mh0_3
      real(default) :: wh0_3
      real(default) :: mh0_4
      real(default) :: wh0_4
      real(default) :: mh0_5
      real(default) :: wh0_5
      real(default) :: mh0_6
      real(default) :: wh0_6
      real(default) :: mh0_7
      real(default) :: wh0_7
      real(default) :: mh0_8
      real(default) :: wh0_8
      real(default) :: mh0_9
      real(default) :: wh0_9
      real(default) :: ma0_1
      real(default) :: wa0_1
      real(default) :: ma0_2
      real(default) :: wa0_2
      real(default) :: ma0_3
      real(default) :: wa0_3
      real(default) :: ma0_4
      real(default) :: wa0_4
      real(default) :: ma0_5
      real(default) :: wa0_5
      real(default) :: ma0_6
      real(default) :: wa0_6
      real(default) :: ma0_7
      real(default) :: wa0_7
      real(default) :: mhpm_1
      real(default) :: whpm_1
      real(default) :: mhpm_2
      real(default) :: whpm_2
      real(default) :: mhpm_3
      real(default) :: whpm_3
      real(default) :: mhpm_4
      real(default) :: whpm_4
      real(default) :: mhpm_5
      real(default) :: whpm_5
      real(default) :: tanb_h
      real(default) :: ls
      real(default) :: a_ls
      real(default) :: ks
      real(default) :: a_ks
      real(default) :: nmu
      real(default) :: al_h
      real(default) :: mixh0_11
      real(default) :: mixh0_12
      real(default) :: mixh0_13
      real(default) :: mixh0_14
      real(default) :: mixh0_15
      real(default) :: mixh0_16
      real(default) :: mixh0_17
      real(default) :: mixh0_18
      real(default) :: mixh0_19
      real(default) :: mixh0_21
      real(default) :: mixh0_22
      real(default) :: mixh0_23
      real(default) :: mixh0_24
      real(default) :: mixh0_25
      real(default) :: mixh0_26
      real(default) :: mixh0_27
      real(default) :: mixh0_28
      real(default) :: mixh0_29
      real(default) :: mixh0_31
      real(default) :: mixh0_32
      real(default) :: mixh0_33
      real(default) :: mixh0_34
      real(default) :: mixh0_35
      real(default) :: mixh0_36
      real(default) :: mixh0_37
      real(default) :: mixh0_38
      real(default) :: mixh0_39
      real(default) :: mixh0_41
      real(default) :: mixh0_42
      real(default) :: mixh0_43
      real(default) :: mixh0_44
      real(default) :: mixh0_45
      real(default) :: mixh0_46
      real(default) :: mixh0_47
      real(default) :: mixh0_48
      real(default) :: mixh0_49
      real(default) :: mixh0_51
      real(default) :: mixh0_52
      real(default) :: mixh0_53
      real(default) :: mixh0_54
      real(default) :: mixh0_55
      real(default) :: mixh0_56
      real(default) :: mixh0_57
      real(default) :: mixh0_58
      real(default) :: mixh0_59
      real(default) :: mixh0_61
      real(default) :: mixh0_62
      real(default) :: mixh0_63
      real(default) :: mixh0_64
      real(default) :: mixh0_65
      real(default) :: mixh0_66
      real(default) :: mixh0_67
      real(default) :: mixh0_68
      real(default) :: mixh0_69
      real(default) :: mixh0_71
      real(default) :: mixh0_72
      real(default) :: mixh0_73
      real(default) :: mixh0_74
      real(default) :: mixh0_75
      real(default) :: mixh0_76
      real(default) :: mixh0_77
      real(default) :: mixh0_78
      real(default) :: mixh0_79
      real(default) :: mixh0_81
      real(default) :: mixh0_82
      real(default) :: mixh0_83
      real(default) :: mixh0_84
      real(default) :: mixh0_85
      real(default) :: mixh0_86
      real(default) :: mixh0_87
      real(default) :: mixh0_88
      real(default) :: mixh0_89
      real(default) :: mixh0_91
      real(default) :: mixh0_92
      real(default) :: mixh0_93
      real(default) :: mixh0_94
      real(default) :: mixh0_95
      real(default) :: mixh0_96
      real(default) :: mixh0_97
      real(default) :: mixh0_98
      real(default) :: mixh0_99
      real(default) :: mixa0_11
      real(default) :: mixa0_12
      real(default) :: mixa0_13
      real(default) :: mixa0_14
      real(default) :: mixa0_15
      real(default) :: mixa0_16
      real(default) :: mixa0_17
      real(default) :: mixa0_18
      real(default) :: mixa0_19
      real(default) :: mixa0_21
      real(default) :: mixa0_22
      real(default) :: mixa0_23
      real(default) :: mixa0_24
      real(default) :: mixa0_25
      real(default) :: mixa0_26
      real(default) :: mixa0_27
      real(default) :: mixa0_28
      real(default) :: mixa0_29
      real(default) :: mixa0_31
      real(default) :: mixa0_32
      real(default) :: mixa0_33
      real(default) :: mixa0_34
      real(default) :: mixa0_35
      real(default) :: mixa0_36
      real(default) :: mixa0_37
      real(default) :: mixa0_38
      real(default) :: mixa0_39
      real(default) :: mixa0_41
      real(default) :: mixa0_42
      real(default) :: mixa0_43
      real(default) :: mixa0_44
      real(default) :: mixa0_45
      real(default) :: mixa0_46
      real(default) :: mixa0_47
      real(default) :: mixa0_48
      real(default) :: mixa0_49
      real(default) :: mixa0_51
      real(default) :: mixa0_52
      real(default) :: mixa0_53
      real(default) :: mixa0_54
      real(default) :: mixa0_55
      real(default) :: mixa0_56
      real(default) :: mixa0_57
      real(default) :: mixa0_58
      real(default) :: mixa0_59
      real(default) :: mixa0_61
      real(default) :: mixa0_62
      real(default) :: mixa0_63
      real(default) :: mixa0_64
      real(default) :: mixa0_65
      real(default) :: mixa0_66
      real(default) :: mixa0_67
      real(default) :: mixa0_68
      real(default) :: mixa0_69
      real(default) :: mixa0_71
      real(default) :: mixa0_72
      real(default) :: mixa0_73
      real(default) :: mixa0_74
      real(default) :: mixa0_75
      real(default) :: mixa0_76
      real(default) :: mixa0_77
      real(default) :: mixa0_78
      real(default) :: mixa0_79
      real(default) :: mixa0_81
      real(default) :: mixa0_82
      real(default) :: mixa0_83
      real(default) :: mixa0_84
      real(default) :: mixa0_85
      real(default) :: mixa0_86
      real(default) :: mixa0_87
      real(default) :: mixa0_88
      real(default) :: mixa0_89
      real(default) :: mixa0_91
      real(default) :: mixa0_92
      real(default) :: mixa0_93
      real(default) :: mixa0_94
      real(default) :: mixa0_95
      real(default) :: mixa0_96
      real(default) :: mixa0_97
      real(default) :: mixa0_98
      real(default) :: mixa0_99
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
      real(default) :: msnumu
      real(default) :: mstau1
      real(default) :: msntau
      real(default) :: mse2
      real(default) :: msmu2
      real(default) :: mstau2
      real(default) :: mlq_11
      real(default) :: mlq_12
      real(default) :: mlq_21
      real(default) :: mlq_22
      real(default) :: mlq_31
      real(default) :: mlq_32
      real(default) :: mlqino_1
      real(default) :: mlqino_2
      real(default) :: mlqino_3
      real(default) :: mgg
      real(default) :: mch1
      real(default) :: mch2
      real(default) :: mch3
      real(default) :: mch4
      real(default) :: mneu1
      real(default) :: mneu2
      real(default) :: mneu3
      real(default) :: mneu4
      real(default) :: mneu5
      real(default) :: mneu6
      real(default) :: mneu7
      real(default) :: mneu8
      real(default) :: mneu9
      real(default) :: mneu10
      real(default) :: mneu11
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
      real(default) :: wsnumu
      real(default) :: wstau1
      real(default) :: wsntau
      real(default) :: wse2
      real(default) :: wsmu2
      real(default) :: wstau2
      real(default) :: wgg
      real(default) :: wch1
      real(default) :: wch2
      real(default) :: wch3
      real(default) :: wch4
      real(default) :: wneu1
      real(default) :: wneu2
      real(default) :: wneu3
      real(default) :: wneu4
      real(default) :: wneu5
      real(default) :: wneu6
      real(default) :: wneu7
      real(default) :: wneu8
      real(default) :: wneu9
      real(default) :: wneu10
      real(default) :: wneu11
      real(default) :: wlq_11
      real(default) :: wlq_12
      real(default) :: wlq_21
      real(default) :: wlq_22
      real(default) :: wlq_31
      real(default) :: wlq_32
      real(default) :: wlqino_1
      real(default) :: wlqino_2
      real(default) :: wlqino_3
      real(default) :: milq1_11
      real(default) :: milq1_12
      real(default) :: milq1_21
      real(default) :: milq1_22
      real(default) :: milq2_11
      real(default) :: milq2_12
      real(default) :: milq2_21
      real(default) :: milq2_22
      real(default) :: milq3_11
      real(default) :: milq3_12
      real(default) :: milq3_21
      real(default) :: milq3_22
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
      real(default) :: mixn_01_01
      real(default) :: mixn_01_02
      real(default) :: mixn_01_03
      real(default) :: mixn_01_04
      real(default) :: mixn_01_05
      real(default) :: mixn_01_06
      real(default) :: mixn_01_07
      real(default) :: mixn_01_08
      real(default) :: mixn_01_09
      real(default) :: mixn_01_10
      real(default) :: mixn_01_11
      real(default) :: mixn_02_01
      real(default) :: mixn_02_02
      real(default) :: mixn_02_03
      real(default) :: mixn_02_04
      real(default) :: mixn_02_05
      real(default) :: mixn_02_06
      real(default) :: mixn_02_07
      real(default) :: mixn_02_08
      real(default) :: mixn_02_09
      real(default) :: mixn_02_10
      real(default) :: mixn_02_11
      real(default) :: mixn_03_01
      real(default) :: mixn_03_02
      real(default) :: mixn_03_03
      real(default) :: mixn_03_04
      real(default) :: mixn_03_05
      real(default) :: mixn_03_06
      real(default) :: mixn_03_07
      real(default) :: mixn_03_08
      real(default) :: mixn_03_09
      real(default) :: mixn_03_10
      real(default) :: mixn_03_11
      real(default) :: mixn_04_01
      real(default) :: mixn_04_02
      real(default) :: mixn_04_03
      real(default) :: mixn_04_04
      real(default) :: mixn_04_05
      real(default) :: mixn_04_06
      real(default) :: mixn_04_07
      real(default) :: mixn_04_08
      real(default) :: mixn_04_09
      real(default) :: mixn_04_10
      real(default) :: mixn_04_11
      real(default) :: mixn_05_01
      real(default) :: mixn_05_02
      real(default) :: mixn_05_03
      real(default) :: mixn_05_04
      real(default) :: mixn_05_05
      real(default) :: mixn_05_06
      real(default) :: mixn_05_07
      real(default) :: mixn_05_08
      real(default) :: mixn_05_09
      real(default) :: mixn_05_10
      real(default) :: mixn_05_11
      real(default) :: mixn_06_01
      real(default) :: mixn_06_02
      real(default) :: mixn_06_03
      real(default) :: mixn_06_04
      real(default) :: mixn_06_05
      real(default) :: mixn_06_06
      real(default) :: mixn_06_07
      real(default) :: mixn_06_08
      real(default) :: mixn_06_09
      real(default) :: mixn_06_10
      real(default) :: mixn_06_11
      real(default) :: mixn_07_01
      real(default) :: mixn_07_02
      real(default) :: mixn_07_03
      real(default) :: mixn_07_04
      real(default) :: mixn_07_05
      real(default) :: mixn_07_06
      real(default) :: mixn_07_07
      real(default) :: mixn_07_08
      real(default) :: mixn_07_09
      real(default) :: mixn_07_10
      real(default) :: mixn_07_11
      real(default) :: mixn_08_01
      real(default) :: mixn_08_02
      real(default) :: mixn_08_03
      real(default) :: mixn_08_04
      real(default) :: mixn_08_05
      real(default) :: mixn_08_06
      real(default) :: mixn_08_07
      real(default) :: mixn_08_08
      real(default) :: mixn_08_09
      real(default) :: mixn_08_10
      real(default) :: mixn_08_11
      real(default) :: mixn_09_01
      real(default) :: mixn_09_02
      real(default) :: mixn_09_03
      real(default) :: mixn_09_04
      real(default) :: mixn_09_05
      real(default) :: mixn_09_06
      real(default) :: mixn_09_07
      real(default) :: mixn_09_08
      real(default) :: mixn_09_09
      real(default) :: mixn_09_10
      real(default) :: mixn_09_11
      real(default) :: mixn_10_01
      real(default) :: mixn_10_02
      real(default) :: mixn_10_03
      real(default) :: mixn_10_04
      real(default) :: mixn_10_05
      real(default) :: mixn_10_06
      real(default) :: mixn_10_07
      real(default) :: mixn_10_08
      real(default) :: mixn_10_09
      real(default) :: mixn_10_10
      real(default) :: mixn_10_11
      real(default) :: mixn_11_01
      real(default) :: mixn_11_02
      real(default) :: mixn_11_03
      real(default) :: mixn_11_04
      real(default) :: mixn_11_05
      real(default) :: mixn_11_06
      real(default) :: mixn_11_07
      real(default) :: mixn_11_08
      real(default) :: mixn_11_09
      real(default) :: mixn_11_10
      real(default) :: mixn_11_11
      real(default) :: mu1_11
      real(default) :: mu1_12
      real(default) :: mu1_21
      real(default) :: mu1_22
      real(default) :: mv1_11
      real(default) :: mv1_12
      real(default) :: mv1_21
      real(default) :: mv1_22
      real(default) :: mu2_11
      real(default) :: mu2_12
      real(default) :: mu2_21
      real(default) :: mu2_22
      real(default) :: mv2_11
      real(default) :: mv2_12
      real(default) :: mv2_21
      real(default) :: mv2_22
      real(default) :: mu3_11
      real(default) :: mu3_12
      real(default) :: mu3_21
      real(default) :: mu3_22
      real(default) :: mv3_11
      real(default) :: mv3_12
      real(default) :: mv3_21
      real(default) :: mv3_22
      real(default) :: v
      real(default) :: cw
      real(default) :: sw
      real(default) :: ee
   end type parameter_set
   real(kind=default) :: sinthw, costhw, qelep, qeup, qedwn, v
   type(parameter_set) :: par    
    par%gf         = par_array(1)
    par%mz         = par_array(2)
    par%mw         = par_array(3)
    par%wZ         = par_array(4)
    par%ww         = par_array(5)
    par%me         = par_array(6)
    par%mmu        = par_array(7)
    par%mtau       = par_array(8)
    par%mc         = par_array(9)
    par%ms         = par_array(10)
    par%mb         = par_array(11)
    par%mtop       = par_array(12)
    par%wtop       = par_array(13)
    par%alphas     = par_array(14)
    par%mtype      = par_array(15)
    par%m_zero     = par_array(16)
    par%m_half     = par_array(17)
    par%a0         = par_array(18)
    par%tanb       = par_array(19)
    par%sgn_mu     = par_array(20)
    par%lambda     = par_array(21)
    par%m_mes      = par_array(22)
    par%n_5        = par_array(23)
    par%c_grav     = par_array(24)
    par%m_grav     = par_array(25)
    par%yuk_d1     = par_array(26)
    par%yuk_d2     = par_array(27)
    par%yuk_dc     = par_array(28)
    par%yuk_nd     = par_array(29)
    par%ae_33      = par_array(30)
    par%au_33      = par_array(31)
    par%ad_33      = par_array(32)
    par%mh0_1      = par_array(33)
    par%wh0_1      = par_array(34)
    par%mh0_2      = par_array(35)
    par%wh0_2      = par_array(36)
    par%mh0_3      = par_array(37)
    par%wh0_3      = par_array(38)
    par%mh0_4      = par_array(39)
    par%wh0_4      = par_array(40)
    par%mh0_5      = par_array(41)
    par%wh0_5      = par_array(42)
    par%mh0_6      = par_array(43)
    par%wh0_6      = par_array(44)
    par%mh0_7      = par_array(45)
    par%wh0_7      = par_array(46)
    par%mh0_8      = par_array(47)
    par%wh0_8      = par_array(48)
    par%mh0_9      = par_array(49)
    par%wh0_9      = par_array(50)
    par%ma0_1      = par_array(51)
    par%wa0_1      = par_array(52)
    par%ma0_2      = par_array(53)
    par%wa0_2      = par_array(54)
    par%ma0_3      = par_array(55)
    par%wa0_3      = par_array(56)
    par%ma0_4      = par_array(57)
    par%wa0_4      = par_array(58)
    par%ma0_5      = par_array(59)
    par%wa0_5      = par_array(60)
    par%ma0_6      = par_array(61)
    par%wa0_6      = par_array(62)
    par%ma0_7      = par_array(63)
    par%wa0_7      = par_array(64)
    par%mhpm_1     = par_array(65)
    par%whpm_1     = par_array(66)
    par%mhpm_2     = par_array(67)
    par%whpm_2     = par_array(68)
    par%mhpm_3     = par_array(69)
    par%whpm_3     = par_array(70)
    par%mhpm_4     = par_array(71)
    par%whpm_4     = par_array(72)
    par%mhpm_5     = par_array(73)
    par%whpm_5     = par_array(74)
    par%tanb_h     = par_array(75)
    par%ls         = par_array(76)
    par%a_ls       = par_array(77)
    par%ks         = par_array(78)
    par%a_ks       = par_array(79)
    par%nmu        = par_array(80)
    par%al_h       = par_array(81)
    par%mixh0_11   = par_array(82)
    par%mixh0_12   = par_array(83)
    par%mixh0_13   = par_array(84)
    par%mixh0_14   = par_array(85)
    par%mixh0_15   = par_array(86)
    par%mixh0_16   = par_array(87)
    par%mixh0_17   = par_array(88)
    par%mixh0_18   = par_array(89)
    par%mixh0_19   = par_array(90)
    par%mixh0_21   = par_array(91)
    par%mixh0_22   = par_array(92)
    par%mixh0_23   = par_array(93)
    par%mixh0_24   = par_array(94)
    par%mixh0_25   = par_array(95)
    par%mixh0_26   = par_array(96)
    par%mixh0_27   = par_array(97)
    par%mixh0_28   = par_array(98)
    par%mixh0_29   = par_array(99)
    par%mixh0_31   = par_array(100)
    par%mixh0_32   = par_array(101)
    par%mixh0_33   = par_array(102)
    par%mixh0_34   = par_array(103)
    par%mixh0_35   = par_array(104)
    par%mixh0_36   = par_array(105)
    par%mixh0_37   = par_array(106)
    par%mixh0_38   = par_array(107)
    par%mixh0_39   = par_array(108)
    par%mixh0_41   = par_array(109)
    par%mixh0_42   = par_array(110)
    par%mixh0_43   = par_array(111)
    par%mixh0_44   = par_array(112)
    par%mixh0_45   = par_array(113)
    par%mixh0_46   = par_array(114)
    par%mixh0_47   = par_array(115)
    par%mixh0_48   = par_array(116)
    par%mixh0_49   = par_array(117)
    par%mixh0_51   = par_array(118)
    par%mixh0_52   = par_array(119)
    par%mixh0_53   = par_array(120)
    par%mixh0_54   = par_array(121)
    par%mixh0_55   = par_array(122)
    par%mixh0_56   = par_array(123)
    par%mixh0_57   = par_array(124)
    par%mixh0_58   = par_array(125)
    par%mixh0_59   = par_array(126)
    par%mixh0_61   = par_array(127)
    par%mixh0_62   = par_array(128)
    par%mixh0_63   = par_array(129)
    par%mixh0_64   = par_array(130)
    par%mixh0_65   = par_array(131)
    par%mixh0_66   = par_array(132)
    par%mixh0_67   = par_array(133)
    par%mixh0_68   = par_array(134)
    par%mixh0_69   = par_array(135)
    par%mixh0_71   = par_array(136)
    par%mixh0_72   = par_array(137)
    par%mixh0_73   = par_array(138)
    par%mixh0_74   = par_array(139)
    par%mixh0_75   = par_array(140)
    par%mixh0_76   = par_array(141)
    par%mixh0_77   = par_array(142)
    par%mixh0_78   = par_array(143)
    par%mixh0_79   = par_array(144)
    par%mixh0_81   = par_array(145)
    par%mixh0_82   = par_array(146)
    par%mixh0_83   = par_array(147)
    par%mixh0_84   = par_array(148)
    par%mixh0_85   = par_array(149)
    par%mixh0_86   = par_array(150)
    par%mixh0_87   = par_array(151)
    par%mixh0_88   = par_array(152)
    par%mixh0_89   = par_array(153)
    par%mixh0_91   = par_array(154)
    par%mixh0_92   = par_array(155)
    par%mixh0_93   = par_array(156)
    par%mixh0_94   = par_array(157)
    par%mixh0_95   = par_array(158)
    par%mixh0_96   = par_array(159)
    par%mixh0_97   = par_array(160)
    par%mixh0_98   = par_array(161)
    par%mixh0_99   = par_array(162)
    par%mixa0_11   = par_array(163)
    par%mixa0_12   = par_array(164)
    par%mixa0_13   = par_array(165)
    par%mixa0_14   = par_array(166)
    par%mixa0_15   = par_array(167)
    par%mixa0_16   = par_array(168)
    par%mixa0_17   = par_array(169)
    par%mixa0_18   = par_array(170)
    par%mixa0_19   = par_array(171)
    par%mixa0_21   = par_array(172)
    par%mixa0_22   = par_array(173)
    par%mixa0_23   = par_array(174)
    par%mixa0_24   = par_array(175)
    par%mixa0_25   = par_array(176)
    par%mixa0_26   = par_array(177)
    par%mixa0_27   = par_array(178)
    par%mixa0_28   = par_array(179)
    par%mixa0_29   = par_array(180)
    par%mixa0_31   = par_array(181)
    par%mixa0_32   = par_array(182)
    par%mixa0_33   = par_array(183)
    par%mixa0_34   = par_array(184)
    par%mixa0_35   = par_array(185)
    par%mixa0_36   = par_array(186)
    par%mixa0_37   = par_array(187)
    par%mixa0_38   = par_array(188)
    par%mixa0_39   = par_array(189)
    par%mixa0_41   = par_array(190)
    par%mixa0_42   = par_array(191)
    par%mixa0_43   = par_array(192)
    par%mixa0_44   = par_array(193)
    par%mixa0_45   = par_array(194)
    par%mixa0_46   = par_array(195)
    par%mixa0_47   = par_array(196)
    par%mixa0_48   = par_array(197)
    par%mixa0_49   = par_array(198)
    par%mixa0_51   = par_array(199)
    par%mixa0_52   = par_array(200)
    par%mixa0_53   = par_array(201)
    par%mixa0_54   = par_array(202)
    par%mixa0_55   = par_array(203)
    par%mixa0_56   = par_array(204)
    par%mixa0_57   = par_array(205)
    par%mixa0_58   = par_array(206)
    par%mixa0_59   = par_array(207)
    par%mixa0_61   = par_array(208)
    par%mixa0_62   = par_array(209)
    par%mixa0_63   = par_array(210)
    par%mixa0_64   = par_array(211)
    par%mixa0_65   = par_array(212)
    par%mixa0_66   = par_array(213)
    par%mixa0_67   = par_array(214)
    par%mixa0_68   = par_array(215)
    par%mixa0_69   = par_array(216)
    par%mixa0_71   = par_array(217)
    par%mixa0_72   = par_array(218)
    par%mixa0_73   = par_array(219)
    par%mixa0_74   = par_array(220)
    par%mixa0_75   = par_array(221)
    par%mixa0_76   = par_array(222)
    par%mixa0_77   = par_array(223)
    par%mixa0_78   = par_array(224)
    par%mixa0_79   = par_array(225)
    par%mixa0_81   = par_array(226)
    par%mixa0_82   = par_array(227)
    par%mixa0_83   = par_array(228)
    par%mixa0_84   = par_array(229)
    par%mixa0_85   = par_array(230)
    par%mixa0_86   = par_array(231)
    par%mixa0_87   = par_array(232)
    par%mixa0_88   = par_array(233)
    par%mixa0_89   = par_array(234)
    par%mixa0_91   = par_array(235)
    par%mixa0_92   = par_array(236)
    par%mixa0_93   = par_array(237)
    par%mixa0_94   = par_array(238)
    par%mixa0_95   = par_array(239)
    par%mixa0_96   = par_array(240)
    par%mixa0_97   = par_array(241)
    par%mixa0_98   = par_array(242)
    par%mixa0_99   = par_array(243)
    par%msu1       = par_array(244)
    par%msd1       = par_array(245)
    par%msc1       = par_array(246)
    par%mss1       = par_array(247)
    par%mstop1     = par_array(248)
    par%msb1       = par_array(249)
    par%msu2       = par_array(250)
    par%msd2       = par_array(251)
    par%msc2       = par_array(252)
    par%mss2       = par_array(253)
    par%mstop2     = par_array(254)
    par%msb2       = par_array(255)
    par%mse1       = par_array(256)
    par%msne       = par_array(257)
    par%msmu1      = par_array(258)
    par%msnumu     = par_array(259)
    par%mstau1     = par_array(260)
    par%msntau     = par_array(261)
    par%mse2       = par_array(262)
    par%msmu2      = par_array(263)
    par%mstau2     = par_array(264)
    par%mlq_11     = par_array(265)
    par%mlq_12     = par_array(266)
    par%mlq_21     = par_array(267)
    par%mlq_22     = par_array(268)
    par%mlq_31     = par_array(269)
    par%mlq_32     = par_array(270)
    par%mlqino_1   = par_array(271)
    par%mlqino_2   = par_array(272)
    par%mlqino_3   = par_array(273)
    par%mgg        = par_array(274)
    par%mch1       = par_array(275)
    par%mch2       = par_array(276)
    par%mch3       = par_array(277)
    par%mch4       = par_array(278)
    par%mneu1      = par_array(279)
    par%mneu2      = par_array(280)
    par%mneu3      = par_array(281)
    par%mneu4      = par_array(282)
    par%mneu5      = par_array(283)
    par%mneu6      = par_array(284)
    par%mneu7      = par_array(285)
    par%mneu8      = par_array(286)
    par%mneu9      = par_array(287)
    par%mneu10     = par_array(288)
    par%mneu11     = par_array(289)
    par%wsu1       = par_array(290)
    par%wsd1       = par_array(291)
    par%wsc1       = par_array(292)
    par%wss1       = par_array(293)
    par%wstop1     = par_array(294)
    par%wsb1       = par_array(295)
    par%wsu2       = par_array(296)
    par%wsd2       = par_array(297)
    par%wsc2       = par_array(298)
    par%wss2       = par_array(299)
    par%wstop2     = par_array(300)
    par%wsb2       = par_array(301)
    par%wse1       = par_array(302)
    par%wsne       = par_array(303)
    par%wsmu1      = par_array(304)
    par%wsnumu     = par_array(305)
    par%wstau1     = par_array(306)
    par%wsntau     = par_array(307)
    par%wse2       = par_array(308)
    par%wsmu2      = par_array(309)
    par%wstau2     = par_array(310)
    par%wgg        = par_array(311)
    par%wch1       = par_array(312)
    par%wch2       = par_array(313)
    par%wch3       = par_array(314)
    par%wch4       = par_array(315)
    par%wneu1      = par_array(316)
    par%wneu2      = par_array(317)
    par%wneu3      = par_array(318)
    par%wneu4      = par_array(319)
    par%wneu5      = par_array(320)
    par%wneu6      = par_array(321)
    par%wneu7      = par_array(322)
    par%wneu8      = par_array(323)
    par%wneu9      = par_array(324)
    par%wneu10     = par_array(325)
    par%wneu11     = par_array(326)
    par%wlq_11     = par_array(327)
    par%wlq_12     = par_array(328)
    par%wlq_21     = par_array(329)
    par%wlq_22     = par_array(330)
    par%wlq_31     = par_array(331)
    par%wlq_32     = par_array(332)
    par%wlqino_1   = par_array(333)
    par%wlqino_2   = par_array(334)
    par%wlqino_3   = par_array(335)
    par%milq1_11   = par_array(336)
    par%milq1_12   = par_array(337)
    par%milq1_21   = par_array(338)
    par%milq1_22   = par_array(339)
    par%milq2_11   = par_array(340)
    par%milq2_12   = par_array(341)
    par%milq2_21   = par_array(342)
    par%milq2_22   = par_array(343)
    par%milq3_11   = par_array(344)
    par%milq3_12   = par_array(345)
    par%milq3_21   = par_array(346)
    par%milq3_22   = par_array(347)
    par%mt_11      = par_array(348)
    par%mt_12      = par_array(349)
    par%mt_21      = par_array(350)
    par%mt_22      = par_array(351)
    par%mb_11      = par_array(352)
    par%mb_12      = par_array(353)
    par%mb_21      = par_array(354)
    par%mb_22      = par_array(355)
    par%ml_11      = par_array(356)
    par%ml_12      = par_array(357)
    par%ml_21      = par_array(358)
    par%ml_22      = par_array(359)
    par%mixn_01_01 = par_array(360)
    par%mixn_01_02 = par_array(361)
    par%mixn_01_03 = par_array(362)
    par%mixn_01_04 = par_array(363)
    par%mixn_01_05 = par_array(364)
    par%mixn_01_06 = par_array(365)
    par%mixn_01_07 = par_array(366)
    par%mixn_01_08 = par_array(367)
    par%mixn_01_09 = par_array(368)
    par%mixn_01_10 = par_array(369)
    par%mixn_01_11 = par_array(370)
    par%mixn_02_01 = par_array(371)
    par%mixn_02_02 = par_array(372)
    par%mixn_02_03 = par_array(373)
    par%mixn_02_04 = par_array(374)
    par%mixn_02_05 = par_array(375)
    par%mixn_02_06 = par_array(376)
    par%mixn_02_07 = par_array(377)
    par%mixn_02_08 = par_array(378)
    par%mixn_02_09 = par_array(379)
    par%mixn_02_10 = par_array(380)
    par%mixn_02_11 = par_array(381)
    par%mixn_03_01 = par_array(382)
    par%mixn_03_02 = par_array(383)
    par%mixn_03_03 = par_array(384)
    par%mixn_03_04 = par_array(385)
    par%mixn_03_05 = par_array(386)
    par%mixn_03_06 = par_array(387)
    par%mixn_03_07 = par_array(388)
    par%mixn_03_08 = par_array(389)
    par%mixn_03_09 = par_array(390)
    par%mixn_03_10 = par_array(391)
    par%mixn_03_11 = par_array(392)
    par%mixn_04_01 = par_array(393)
    par%mixn_04_02 = par_array(394)
    par%mixn_04_03 = par_array(395)
    par%mixn_04_04 = par_array(396)
    par%mixn_04_05 = par_array(397)
    par%mixn_04_06 = par_array(398)
    par%mixn_04_07 = par_array(399)
    par%mixn_04_08 = par_array(400)
    par%mixn_04_09 = par_array(401)
    par%mixn_04_10 = par_array(402)
    par%mixn_04_11 = par_array(403)
    par%mixn_05_01 = par_array(404)
    par%mixn_05_02 = par_array(405)
    par%mixn_05_03 = par_array(406)
    par%mixn_05_04 = par_array(407)
    par%mixn_05_05 = par_array(408)
    par%mixn_05_06 = par_array(409)
    par%mixn_05_07 = par_array(410)
    par%mixn_05_08 = par_array(411)
    par%mixn_05_09 = par_array(412)
    par%mixn_05_10 = par_array(413)
    par%mixn_05_11 = par_array(414)
    par%mixn_06_01 = par_array(415)
    par%mixn_06_02 = par_array(416)
    par%mixn_06_03 = par_array(417)
    par%mixn_06_04 = par_array(418)
    par%mixn_06_05 = par_array(419)
    par%mixn_06_06 = par_array(420)
    par%mixn_06_07 = par_array(421)
    par%mixn_06_08 = par_array(422)
    par%mixn_06_09 = par_array(423)
    par%mixn_06_10 = par_array(424)
    par%mixn_06_11 = par_array(425)
    par%mixn_07_01 = par_array(426)
    par%mixn_07_02 = par_array(427)
    par%mixn_07_03 = par_array(428)
    par%mixn_07_04 = par_array(429)
    par%mixn_07_05 = par_array(430)
    par%mixn_07_06 = par_array(431)
    par%mixn_07_07 = par_array(432)
    par%mixn_07_08 = par_array(433)
    par%mixn_07_09 = par_array(434)
    par%mixn_07_10 = par_array(435)
    par%mixn_07_11 = par_array(436)
    par%mixn_08_01 = par_array(437)
    par%mixn_08_02 = par_array(438)
    par%mixn_08_03 = par_array(439)
    par%mixn_08_04 = par_array(440)
    par%mixn_08_05 = par_array(441)
    par%mixn_08_06 = par_array(442)
    par%mixn_08_07 = par_array(443)
    par%mixn_08_08 = par_array(444)
    par%mixn_08_09 = par_array(445)
    par%mixn_08_10 = par_array(446)
    par%mixn_08_11 = par_array(447)
    par%mixn_09_01 = par_array(448)
    par%mixn_09_02 = par_array(449)
    par%mixn_09_03 = par_array(450)
    par%mixn_09_04 = par_array(451)
    par%mixn_09_05 = par_array(452)
    par%mixn_09_06 = par_array(453)
    par%mixn_09_07 = par_array(454)
    par%mixn_09_08 = par_array(455)
    par%mixn_09_09 = par_array(456)
    par%mixn_09_10 = par_array(457)
    par%mixn_09_11 = par_array(458)
    par%mixn_10_01 = par_array(459)
    par%mixn_10_02 = par_array(460)
    par%mixn_10_03 = par_array(461)
    par%mixn_10_04 = par_array(462)
    par%mixn_10_05 = par_array(463)
    par%mixn_10_06 = par_array(464)
    par%mixn_10_07 = par_array(465)
    par%mixn_10_08 = par_array(466)
    par%mixn_10_09 = par_array(467)
    par%mixn_10_10 = par_array(468)
    par%mixn_10_11 = par_array(469)
    par%mixn_11_01 = par_array(470)
    par%mixn_11_02 = par_array(471)
    par%mixn_11_03 = par_array(472)
    par%mixn_11_04 = par_array(473)
    par%mixn_11_05 = par_array(474)
    par%mixn_11_06 = par_array(475)
    par%mixn_11_07 = par_array(476)
    par%mixn_11_08 = par_array(477)
    par%mixn_11_09 = par_array(478)
    par%mixn_11_10 = par_array(479)
    par%mixn_11_11 = par_array(480)
    par%mu1_11     = par_array(481)
    par%mu1_12     = par_array(482)
    par%mu1_21     = par_array(483)
    par%mu1_22     = par_array(484)
    par%mv1_11     = par_array(485)
    par%mv1_12     = par_array(486)
    par%mv1_21     = par_array(487)
    par%mv1_22     = par_array(488)
    par%mu2_11     = par_array(489)
    par%mu2_12     = par_array(490)
    par%mu2_21     = par_array(491)
    par%mu2_22     = par_array(492)
    par%mv2_11     = par_array(493)
    par%mv2_12     = par_array(494)
    par%mv2_21     = par_array(495)
    par%mv2_22     = par_array(496)
    par%mu3_11     = par_array(497)
    par%mu3_12     = par_array(498)
    par%mu3_21     = par_array(499)
    par%mu3_22     = par_array(500)
    par%mv3_11     = par_array(501)
    par%mv3_12     = par_array(502)
    par%mv3_21     = par_array(503)
    par%mv3_22     = par_array(504)
    par%v          = par_array(505)
    par%cw         = par_array(506)
    par%sw         = par_array(507)
    par%ee         = par_array(508)
    mass(1:99) = 0               
    width(1:99) = 0
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
    mass(25) = par%mh0_1
    width(25) = par%wh0_1
    mass(26) =  xi0 * mass(23)
    width(26) =  0
    mass(27) =  xipm * mass(24)
    width(27) =  0
    mass(35) = par%mh0_2
    width(35) = par%wh0_2
    mass(36) = par%mA0_1
    width(36) = par%wA0_1
    mass(37) = par%mhpm_1
    width(37) = par%whpm_1
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
    mass(56) = par%msnumu
    width(56) = par%wsnumu
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
    mass(69) = abs(par%mneu5)
    width(69) = par%wneu5
    mass(70) = abs(par%mch1)
    width(70) = par%wch1
    mass(71) = abs(par%mch2)
    width(71) = par%wch2
    mass(72) = par%mA0_2
    width(72) = par%wA0_2
    mass(73) = par%mh0_3
    width(73) = par%wh0_3
    mass(80) = par%mlq_11
    width(80) = par%wlq_11
    mass(81) = par%mlq_12
    width(81) = par%wlq_12
    mass(82) = par%mlq_21
    width(82) = par%wlq_21
    mass(83) = par%mlq_22
    width(83) = par%wlq_22
    mass(84) = par%mlq_31
    width(84) = par%wlq_31
    mass(85) = par%mlq_32
    width(85) = par%wlq_32
    mass(86) = par%mlqino_1
    width(86) = par%wlqino_1
    mass(87) = par%mlqino_2
    width(87) = par%wlqino_2
    mass(88) = par%mlqino_3
    width(88) = par%wlqino_3
    mass(89) = par%mA0_3
    width(89) = par%wA0_3
    mass(90) = par%mA0_4
    width(90) = par%wA0_4
    mass(91) = par%mA0_5
    width(91) = par%wA0_5
    mass(92) = par%mA0_6
    width(92) = par%wA0_6
    mass(93) = par%mA0_7
    width(93) = par%wA0_7
    mass(94) = par%mh0_4
    width(94) = par%wh0_4
    mass(95) = par%mh0_5
    width(95) = par%wh0_5
    mass(96) = par%mh0_6
    width(96) = par%wh0_6
    mass(97) = par%mh0_7
    width(97) = par%wh0_7
    mass(98) = par%mh0_8
    width(98) = par%wh0_8
    mass(99) = par%mh0_9
    width(99) = par%wh0_9
    mass(100) = abs(par%mneu6)
    width(100) = par%wneu6
    mass(101) = abs(par%mneu7)
    width(101) = par%wneu7
    mass(102) = abs(par%mneu8)
    width(102) = par%wneu8
    mass(103) = abs(par%mneu9)
    width(103) = par%wneu9
    mass(104) = abs(par%mneu10)
    width(104) = par%wneu10
    mass(105) = abs(par%mneu11)
    width(105) = par%wneu11
    mass(106) = abs(par%mch3)
    width(106) = par%wch3
    mass(107) = abs(par%mch4)
    width(107) = par%wch4
    sigch1    = sign (1._default, par%mch1) 
    sigch2    = sign (1._default, par%mch2)
    sigch3    = sign (1._default, par%mch3) 
    sigch4    = sign (1._default, par%mch4)
    sign1     = sign (1, int(par%mneu1))
    sign2     = sign (1, int(par%mneu2))
    sign3     = sign (1, int(par%mneu3))
    sign4     = sign (1, int(par%mneu4))
    sign5     = sign (1, int(par%mneu5))
    sign6     = sign (1, int(par%mneu6))
    sign7     = sign (1, int(par%mneu7))
    sign8     = sign (1, int(par%mneu8))
    sign9     = sign (1, int(par%mneu9))
    sign10    = sign (1, int(par%mneu10))
    sign11    = sign (1, int(par%mneu11))
    vckm(1,1)  = 1
    vckm(1,2)  = 0
    vckm(1,3)  = 0
    vckm(2,1)  = 0
    vckm(2,2)  = 1
    vckm(2,3)  = 0
    vckm(3,1)  = 0
    vckm(3,2)  = 0
    vckm(3,3)  = 1
    v = 2 * par%mW * par%sw / par%ee
    e = par%ee
    !!! This should not be color flow basis !!!
    as = par%alphas
    tanb = par%tanb_h
!!!Higgs Mixing
    sina = sin(par%al_h)
    cosa = cos(par%al_h)

    mix_h0(1,1) = par%mixh0_11
    mix_h0(1,2) = par%mixh0_12
    mix_h0(1,3) = par%mixh0_13
    mix_h0(1,4) = par%mixh0_14
    mix_h0(1,5) = par%mixh0_15
    mix_h0(1,6) = par%mixh0_16
    mix_h0(1,7) = par%mixh0_17
    mix_h0(1,8) = par%mixh0_18
    mix_h0(1,9) = par%mixh0_19
    mix_h0(2,1) = par%mixh0_21
    mix_h0(2,2) = par%mixh0_22
    mix_h0(2,3) = par%mixh0_23
    mix_h0(2,4) = par%mixh0_24
    mix_h0(2,5) = par%mixh0_25
    mix_h0(2,6) = par%mixh0_26
    mix_h0(2,7) = par%mixh0_27
    mix_h0(2,8) = par%mixh0_28
    mix_h0(2,9) = par%mixh0_29
    mix_h0(3,1) = par%mixh0_31
    mix_h0(3,2) = par%mixh0_32
    mix_h0(3,3) = par%mixh0_33
    mix_h0(3,4) = par%mixh0_34
    mix_h0(3,5) = par%mixh0_35
    mix_h0(3,6) = par%mixh0_36
    mix_h0(3,7) = par%mixh0_37
    mix_h0(3,8) = par%mixh0_38
    mix_h0(3,9) = par%mixh0_39
    mix_h0(4,1) = par%mixh0_41
    mix_h0(4,2) = par%mixh0_42
    mix_h0(4,3) = par%mixh0_43
    mix_h0(4,4) = par%mixh0_44
    mix_h0(4,5) = par%mixh0_45
    mix_h0(4,6) = par%mixh0_46
    mix_h0(4,7) = par%mixh0_47
    mix_h0(4,8) = par%mixh0_48
    mix_h0(4,9) = par%mixh0_49
    mix_h0(5,1) = par%mixh0_51
    mix_h0(5,2) = par%mixh0_52
    mix_h0(5,3) = par%mixh0_53
    mix_h0(5,4) = par%mixh0_54
    mix_h0(5,5) = par%mixh0_55
    mix_h0(5,6) = par%mixh0_56
    mix_h0(5,7) = par%mixh0_57
    mix_h0(5,8) = par%mixh0_58
    mix_h0(5,9) = par%mixh0_59
    mix_h0(6,1) = par%mixh0_61
    mix_h0(6,2) = par%mixh0_62
    mix_h0(6,3) = par%mixh0_63
    mix_h0(6,4) = par%mixh0_64
    mix_h0(6,5) = par%mixh0_65
    mix_h0(6,6) = par%mixh0_66
    mix_h0(6,7) = par%mixh0_67
    mix_h0(6,8) = par%mixh0_68
    mix_h0(6,9) = par%mixh0_69
    mix_h0(7,1) = par%mixh0_71
    mix_h0(7,2) = par%mixh0_72
    mix_h0(7,3) = par%mixh0_73
    mix_h0(7,4) = par%mixh0_74
    mix_h0(7,5) = par%mixh0_75
    mix_h0(7,6) = par%mixh0_76
    mix_h0(7,7) = par%mixh0_77
    mix_h0(7,8) = par%mixh0_78
    mix_h0(7,9) = par%mixh0_79
    mix_h0(8,1) = par%mixh0_81
    mix_h0(8,2) = par%mixh0_82
    mix_h0(8,3) = par%mixh0_83
    mix_h0(8,4) = par%mixh0_84
    mix_h0(8,5) = par%mixh0_85
    mix_h0(8,6) = par%mixh0_86
    mix_h0(8,7) = par%mixh0_87
    mix_h0(8,8) = par%mixh0_88
    mix_h0(8,9) = par%mixh0_89
    mix_h0(9,1) = par%mixh0_91
    mix_h0(9,2) = par%mixh0_92
    mix_h0(9,3) = par%mixh0_93
    mix_h0(9,4) = par%mixh0_94
    mix_h0(9,5) = par%mixh0_95
    mix_h0(9,6) = par%mixh0_96
    mix_h0(9,7) = par%mixh0_97
    mix_h0(9,8) = par%mixh0_98
    mix_h0(9,9) = par%mixh0_99


    mix_A0(1,1) = par%mixa0_11
    mix_A0(1,2) = par%mixa0_12
    mix_A0(1,3) = par%mixa0_13
    mix_A0(1,4) = par%mixa0_14
    mix_A0(1,5) = par%mixa0_15
    mix_A0(1,6) = par%mixa0_16
    mix_A0(1,7) = par%mixa0_17
    mix_A0(1,8) = par%mixa0_18
    mix_A0(1,9) = par%mixa0_19
    mix_A0(2,1) = par%mixa0_21
    mix_A0(2,2) = par%mixa0_22
    mix_A0(2,3) = par%mixa0_23
    mix_A0(2,4) = par%mixa0_24
    mix_A0(2,5) = par%mixa0_25
    mix_A0(2,6) = par%mixa0_26
    mix_A0(2,7) = par%mixa0_27
    mix_A0(2,8) = par%mixa0_28
    mix_A0(2,9) = par%mixa0_29
    mix_A0(3,1) = par%mixa0_31
    mix_A0(3,2) = par%mixa0_32
    mix_A0(3,3) = par%mixa0_33
    mix_A0(3,4) = par%mixa0_34
    mix_A0(3,5) = par%mixa0_35
    mix_A0(3,6) = par%mixa0_36
    mix_A0(3,7) = par%mixa0_37
    mix_A0(3,8) = par%mixa0_38
    mix_A0(3,9) = par%mixa0_39
    mix_A0(4,1) = par%mixa0_41
    mix_A0(4,2) = par%mixa0_42
    mix_A0(4,3) = par%mixa0_43
    mix_A0(4,4) = par%mixa0_44
    mix_A0(4,5) = par%mixa0_45
    mix_A0(4,6) = par%mixa0_46
    mix_A0(4,7) = par%mixa0_47
    mix_A0(4,8) = par%mixa0_48
    mix_A0(4,9) = par%mixa0_49
    mix_A0(5,1) = par%mixa0_51
    mix_A0(5,2) = par%mixa0_52
    mix_A0(5,3) = par%mixa0_53
    mix_A0(5,4) = par%mixa0_54
    mix_A0(5,5) = par%mixa0_55
    mix_A0(5,6) = par%mixa0_56
    mix_A0(5,7) = par%mixa0_57
    mix_A0(5,8) = par%mixa0_58
    mix_A0(5,9) = par%mixa0_59
    mix_A0(6,1) = par%mixa0_61
    mix_A0(6,2) = par%mixa0_62
    mix_A0(6,3) = par%mixa0_63
    mix_A0(6,4) = par%mixa0_64
    mix_A0(6,5) = par%mixa0_65
    mix_A0(6,6) = par%mixa0_66
    mix_A0(6,7) = par%mixa0_67
    mix_A0(6,8) = par%mixa0_68
    mix_A0(6,9) = par%mixa0_69
    mix_A0(7,1) = par%mixa0_71
    mix_A0(7,2) = par%mixa0_72
    mix_A0(7,3) = par%mixa0_73
    mix_A0(7,4) = par%mixa0_74
    mix_A0(7,5) = par%mixa0_75
    mix_A0(7,6) = par%mixa0_76
    mix_A0(7,7) = par%mixa0_77
    mix_A0(7,8) = par%mixa0_78
    mix_A0(7,9) = par%mixa0_79
    mix_A0(8,1) = par%mixa0_81
    mix_A0(8,2) = par%mixa0_82
    mix_A0(8,3) = par%mixa0_83
    mix_A0(8,4) = par%mixa0_84
    mix_A0(8,5) = par%mixa0_85
    mix_A0(8,6) = par%mixa0_86
    mix_A0(8,7) = par%mixa0_87
    mix_A0(8,8) = par%mixa0_88
    mix_A0(8,9) = par%mixa0_89
    mix_A0(9,1) = par%mixa0_91
    mix_A0(9,2) = par%mixa0_92
    mix_A0(9,3) = par%mixa0_93
    mix_A0(9,4) = par%mixa0_94
    mix_A0(9,5) = par%mixa0_95
    mix_A0(9,6) = par%mixa0_96
    mix_A0(9,7) = par%mixa0_97
    mix_A0(9,8) = par%mixa0_98
    mix_A0(9,9) = par%mixa0_99




 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
            
    select case (sign1)
       case (1)
          eta1 = (1.0_default,0.0_default)
       case (-1)
          eta1 = (0.0_default,1.0_default)
       case default 
          print *, 'sign1', sign1
          stop "parameters_PSSSM: No definite sign neutralino1"
    end select
    select case (sign2)
       case (1)
          eta2 = (1.0_default,0.0_default)
       case (-1)
          eta2 = (0.0_default,1.0_default)
       case default 
          print *, 'sign2', sign2
          stop "parameters_PSSSM: No definite sign neutralino2"
    end select
    select case (sign3)
       case (1)
          eta3 = (1.0_default,0.0_default)
       case (-1)
          eta3 = (0.0_default,1.0_default)
       case default 
          print *, 'sign3', sign3
          stop "parameters_PSSSM: No definite sign neutralino3"
    end select
    select case (sign4)
       case (1)
          eta4 = (1.0_default,0.0_default)
       case (-1)
          eta4 = (0.0_default,1.0_default)
       case default 
          print *, 'sign4', sign4
          stop "parameters_PSSSM: No definite sign neutralino4"
    end select
    select case (sign5)
       case (1)
          eta5 = (1.0_default,0.0_default)
       case (-1)
          eta5 = (0.0_default,1.0_default)
       case default 
          print *, 'sign5', sign5
          stop "parameters_PSSSM: No definite sign neutralino5"
    end select
    select case (sign6)
       case (1)
          eta6 = (1.0_default,0.0_default)
       case (-1)
          eta6 = (0.0_default,1.0_default)
       case default 
          print *, 'sign6', sign6
          stop "parameters_PSSSM: No definite sign neutralino6"
    end select
    select case (sign7)
       case (1)
          eta7 = (1.0_default,0.0_default)
       case (-1)
          eta7 = (0.0_default,1.0_default)
       case default 
          print *, 'sign7', sign7
          stop "parameters_PSSSM: No definite sign neutralino7"
    end select
    select case (sign8)
       case (1)
          eta8 = (1.0_default,0.0_default)
       case (-1)
          eta8 = (0.0_default,1.0_default)
       case default 
          print *, 'sign8', sign8
          stop "parameters_PSSSM: No definite sign neutralino8"
    end select
    select case (sign9)
       case (1)
          eta9 = (1.0_default,0.0_default)
       case (-1)
          eta9 = (0.0_default,1.0_default)
       case default 
          print *, 'sign9', sign9
          stop "parameters_PSSSM: No definite sign neutralino9"
    end select
    select case (sign10)
       case (1)
          etat = (1.0_default,0.0_default)
       case (-1)
          etat = (0.0_default,1.0_default)
       case default 
          print *, 'sign10', sign10
          stop "parameters_PSSSM: No definite sign neutralino10"
    end select
    select case (sign11)
       case (1)
          etae = (1.0_default,0.0_default)
       case (-1)
          etae = (0.0_default,1.0_default)
       case default 
          print *, 'sign11', sign11
          stop "parameters_PSSSM: No definite sign neutralino11"
    end select

    sinthw = par%sw
    sin2thw = sinthw**2
    costhw = par%cw
    qelep = - 1.0_default
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  MSSM Limit of the NMSSM, setting the superpotential terms
    !!!  (and corresponding A terms) equal to the mu term, and
    !!!  shrinking the Higgs mixing matrices to the MSSM case
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! lambda=0
    !!! k=0
    !!! A_k=0
    !!! A_lambda=0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! mix_A0(1,1) = sinbe
    !!! mix_A0(1,2) = cosbe
    !!! !!!
    !!! mix_h0(1,1) = -sina
    !!! mix_h0(1,2) = cosa
    !!! mix_h0(2,2) = sina
    !!! mix_h0(2,1) = cosa
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call setup_parameters1
    call setup_parameters2
    call setup_parameters3
    call setup_parameters4
    call setup_parameters5
    call setup_parameters6
    call setup_parameters7
    call setup_parameters8
    call setup_parameters9
contains
 subroutine setup_parameters1 ()
    ! e = sqrt ((4.0_default * PI * alpha)) *** This is predefined! ***
    g = e / sinthw
    gz = (g / costhw)
    !!! Color flow basis, divided by sqrt(2)
    gs = sqrt (2.0_default * PI * par%alphas)
    igs = (imago * gs)
    vev = ((2.0_default * mass(24)) / g)
    q_lep  = (- 1.0_default) 
    q_up   = 2.0_default / 3.0_default
    q_down = -( 1.0_default / 3.0_default)    
    qlep = - e * qelep   !!! This is the negative particle charge !!! 
    qup = - e * qeup     !!! This is the negative particle charge !!! 
    qdwn = - e * qedwn   !!! This is the negative particle charge !!! 
    qchar = ( - e)       !!! This is the negative particle charge !!! 
    gcc = g / (2.0_default * sqrt (2.0_default))
    gssq = (gs / sqrt (2.0_default))
    iqw = imago * e
    igzww = imago * g * costhw
    gw4 = (g**2)
    gzzww = ((g**2) * (costhw**2))
    gppww = (e**2)
    gpzww = (e * g * costhw) 
    sinbe = (tanb / sqrt ((1.0_default + (tanb**2))))
    cosbe = (1.0_default / sqrt ((1.0_default + (tanb**2))))
    eidelta = (cosd + (imago * sind))
    cos2be = ((cosbe**2) - (sinbe**2))
    sin2be = (2.0_default * cosbe * sinbe)

!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Neutralino Mixing
!!!!!!!!!!!!!!!!!!!!!
    mix_neu(1,1) = eta1 * par%mixn_01_01
    mix_neu(1,2) = eta1 * par%mixn_01_02
    mix_neu(1,3) = eta1 * par%mixn_01_03
    mix_neu(1,4) = eta1 * par%mixn_01_04
    mix_neu(1,5) = eta1 * par%mixn_01_05
    mix_neu(1,6) = eta1 * par%mixn_01_06
    mix_neu(1,7) = eta1 * par%mixn_01_07
    mix_neu(1,8) = eta1 * par%mixn_01_08
    mix_neu(1,9) = eta1 * par%mixn_01_09
    mix_neu(1,10) = eta1 * par%mixn_10_10
    mix_neu(1,11) = eta1 * par%mixn_10_11
    mix_neu(2,1) = eta2 * par%mixn_02_01
    mix_neu(2,2) = eta2 * par%mixn_02_02
    mix_neu(2,3) = eta2 * par%mixn_02_03
    mix_neu(2,4) = eta2 * par%mixn_02_04
    mix_neu(2,5) = eta2 * par%mixn_02_05
    mix_neu(2,6) = eta2 * par%mixn_02_06
    mix_neu(2,7) = eta2 * par%mixn_02_07
    mix_neu(2,8) = eta2 * par%mixn_02_08
    mix_neu(2,9) = eta2 * par%mixn_02_09
    mix_neu(2,10) = eta2 * par%mixn_02_10
    mix_neu(2,11) = eta2 * par%mixn_02_11
    mix_neu(3,1) = eta3 * par%mixn_03_01
    mix_neu(3,2) = eta3 * par%mixn_03_02
    mix_neu(3,3) = eta3 * par%mixn_03_03
    mix_neu(3,4) = eta3 * par%mixn_03_04
    mix_neu(3,5) = eta3 * par%mixn_03_05
    mix_neu(3,6) = eta3 * par%mixn_03_06
    mix_neu(3,7) = eta3 * par%mixn_03_07
    mix_neu(3,8) = eta3 * par%mixn_03_08
    mix_neu(3,9) = eta3 * par%mixn_03_09
    mix_neu(3,10) = eta3 * par%mixn_03_10
    mix_neu(3,11) = eta3 * par%mixn_03_11
    mix_neu(4,1) = eta4 * par%mixn_04_01
    mix_neu(4,2) = eta4 * par%mixn_04_02
    mix_neu(4,3) = eta4 * par%mixn_04_03
    mix_neu(4,4) = eta4 * par%mixn_04_04
    mix_neu(4,5) = eta4 * par%mixn_04_05
    mix_neu(4,6) = eta4 * par%mixn_04_06
    mix_neu(4,7) = eta4 * par%mixn_04_07
    mix_neu(4,8) = eta4 * par%mixn_04_08
    mix_neu(4,9) = eta4 * par%mixn_04_09
    mix_neu(4,10) = eta4 * par%mixn_04_10
    mix_neu(4,11) = eta4 * par%mixn_04_11
    mix_neu(5,1) = eta5 * par%mixn_05_01
    mix_neu(5,2) = eta5 * par%mixn_05_02
    mix_neu(5,3) = eta5 * par%mixn_05_03
    mix_neu(5,4) = eta5 * par%mixn_05_04
    mix_neu(5,5) = eta5 * par%mixn_05_05
    mix_neu(5,6) = eta5 * par%mixn_05_06
    mix_neu(5,7) = eta5 * par%mixn_05_07
    mix_neu(5,8) = eta5 * par%mixn_05_08
    mix_neu(5,9) = eta5 * par%mixn_05_09
    mix_neu(5,10) = eta5 * par%mixn_05_10
    mix_neu(5,11) = eta5 * par%mixn_05_11
    mix_neu(6,1) = eta6 * par%mixn_06_01
    mix_neu(6,2) = eta6 * par%mixn_06_02
    mix_neu(6,3) = eta6 * par%mixn_06_03
    mix_neu(6,4) = eta6 * par%mixn_06_04
    mix_neu(6,5) = eta6 * par%mixn_06_05
    mix_neu(6,6) = eta6 * par%mixn_06_06
    mix_neu(6,7) = eta6 * par%mixn_06_07
    mix_neu(6,8) = eta6 * par%mixn_06_08
    mix_neu(6,9) = eta6 * par%mixn_06_09
    mix_neu(6,10) = eta6 * par%mixn_06_10
    mix_neu(6,11) = eta6 * par%mixn_06_11
    mix_neu(7,1) = eta7 * par%mixn_07_01
    mix_neu(7,2) = eta7 * par%mixn_07_02
    mix_neu(7,3) = eta7 * par%mixn_07_03
    mix_neu(7,4) = eta7 * par%mixn_07_04
    mix_neu(7,5) = eta7 * par%mixn_07_05
    mix_neu(7,6) = eta7 * par%mixn_07_06
    mix_neu(7,7) = eta7 * par%mixn_07_07
    mix_neu(7,8) = eta7 * par%mixn_07_08
    mix_neu(7,9) = eta7 * par%mixn_07_09
    mix_neu(7,10) = eta7 * par%mixn_07_10
    mix_neu(7,11) = eta7 * par%mixn_07_11
    mix_neu(8,1) = eta8 * par%mixn_08_01
    mix_neu(8,2) = eta8 * par%mixn_08_02
    mix_neu(8,3) = eta8 * par%mixn_08_03
    mix_neu(8,4) = eta8 * par%mixn_08_04
    mix_neu(8,5) = eta8 * par%mixn_08_05
    mix_neu(8,6) = eta8 * par%mixn_08_06
    mix_neu(8,7) = eta8 * par%mixn_08_07
    mix_neu(8,8) = eta8 * par%mixn_08_08
    mix_neu(8,9) = eta8 * par%mixn_08_09
    mix_neu(8,10) = eta8 * par%mixn_08_10
    mix_neu(8,11) = eta8 * par%mixn_08_11
    mix_neu(9,1) = eta9 * par%mixn_09_01
    mix_neu(9,2) = eta9 * par%mixn_09_02
    mix_neu(9,3) = eta9 * par%mixn_09_03
    mix_neu(9,4) = eta9 * par%mixn_09_04
    mix_neu(9,5) = eta9 * par%mixn_09_05
    mix_neu(9,6) = eta9 * par%mixn_09_06
    mix_neu(9,7) = eta9 * par%mixn_09_07
    mix_neu(9,8) = eta9 * par%mixn_09_08
    mix_neu(9,9) = eta9 * par%mixn_09_09
    mix_neu(9,10) = eta9 * par%mixn_09_10
    mix_neu(9,11) = eta9 * par%mixn_09_11
    mix_neu(10,1) = etat * par%mixn_10_01
    mix_neu(10,2) = etat * par%mixn_10_02
    mix_neu(10,3) = etat * par%mixn_10_03
    mix_neu(10,4) = etat * par%mixn_10_04
    mix_neu(10,5) = etat * par%mixn_10_05
    mix_neu(10,6) = etat * par%mixn_10_06
    mix_neu(10,7) = etat * par%mixn_10_07
    mix_neu(10,8) = etat * par%mixn_10_08
    mix_neu(10,9) = etat * par%mixn_10_09
    mix_neu(10,10) = etat * par%mixn_10_10
    mix_neu(10,11) = etat * par%mixn_10_11
    mix_neu(11,1) = etae * par%mixn_11_01
    mix_neu(11,2) = etae * par%mixn_11_02
    mix_neu(11,3) = etae * par%mixn_11_03
    mix_neu(11,4) = etae * par%mixn_11_04
    mix_neu(11,5) = etae * par%mixn_11_05
    mix_neu(11,6) = etae * par%mixn_11_06
    mix_neu(11,7) = etae * par%mixn_11_07
    mix_neu(11,8) = etae * par%mixn_11_08
    mix_neu(11,9) = etae * par%mixn_11_09
    mix_neu(11,10) = etae * par%mixn_11_10
    mix_neu(11,11) = etae * par%mixn_11_11

    !!! Checked by JR !!! 
    mix_charU(1,1) = par%mu1_11                !!! Rotat. matrix containing phi_R
    mix_charU(1,2) = par%mu1_12                !!! Rotat. matrix containing phi_R
    mix_charU(2,1) = par%mu1_21                !!! Rotat. matrix containing phi_R
    mix_charU(2,2) = par%mu1_22                !!! Rotat. matrix containing phi_R
    mix_charV(1,1) = sigch1 * par%mv1_11       !!! Rotat. matrix containing phi_L
    mix_charV(1,2) = sigch1 * par%mv1_12       !!! Rotat. matrix containing phi_L
    mix_charV(2,1) = sigch2 * par%mv1_21       !!! Rotat. matrix containing phi_L
    mix_charV(2,2) = sigch2 * par%mv1_22       !!! Rotat. matrix containing phi_L
    al(1) = 0
    au(1) = 0
    ad(1) = 0
    al(2) = 0
    au(2) = 0
    ad(2) = 0
    !!! SLHA has a different sign for the trilinear scalar parameters
    al(3) = par%Ae_33
    au(3) = par%Au_33
    ad(3) = par%Ad_33
!!! NMSSM par's
    mu = par%nmu
    lambda = par%ls   
    A_lambda = par%a_ls  
    k = par%ks
    A_k = par%a_ks

    if (lambda == 0) then
    r = 0.0_default
    else
    r = mu/(vev * lambda)
    end if

!!!E6 - vev-vector

    do i=1,3
       delta3(i,i) = 1
       vevs(-2+3*i) = vev*cosbe
       vevs(-1+3*i) = vev*sinbe
       vevs(3*i) = vev*r
    end do

!!! E6 Yukawa 
!!! coupling Test !!!!!!!!!!!!!!!!!!!!!!!!!!!

    
    do i = 1,3
       do j = 1,3
          do t = 1,3      
!!! Symmetrisch !!!
!!!                g_yuk_n_d(i,i,i) = par%yuk_nd
!!!                g_yuk_n_h(i,i,i) = 1
!!!                g_yuk_d_1(i,i,i) = par%yuk_d1 !yuk lq ec uc
!!!                g_yuk_d_2(i,i,i) = par%yuk_d2
!!!                g_yuk_d_c(i,i,i) = par%yuk_dc !---"---"---
!!!                g_yuk_e(i,i,i) = 0.100914_default !mass(15)/(vev*cosbe)
!!!                g_yuk_u(i,i,i) = 0.889738_default  !mass(6)/(vev*sinbe)
!!!                g_yuk_d(i,i,i) = 0.139733_default !mass(5)/(vev*cosbe)
!!!                g_yuk_nu(i,i,i) = 0.01_default

!!! Antisymmetrisch !!!
       if( i /= j .and. i /= t .and. j /= t) then
          
          g_yuk_n_d(i,j,t) = par%yuk_nd
          g_yuk_n_h(i,j,t) = 1
          g_yuk_d_1(i,j,t) = par%yuk_d1 !yuk lq ec uc
          g_yuk_d_2(i,j,t) = par%yuk_d2
          g_yuk_d_c(i,j,t) = par%yuk_dc !---"---"---
        !!! careful: SM Yuk's should not be offdiag !!!
          g_yuk_e(i,j,t) = 0.100914_default !mass(15)/(vev*cosbe)
          g_yuk_u(i,j,t) = 0.889738_default  !mass(6)/(vev*sinbe)
          g_yuk_d(i,j,t) = 0.139733_default !mass(5)/(vev*cosbe)
          g_yuk_nu(i,j,t) = 0.01_default
 
       else
 
          g_yuk_n_d(i,j,t) = 0
          g_yuk_n_h(i,j,t) = 0
          g_yuk_d_1(i,j,t) = 0
          g_yuk_d_2(i,j,t) = 0
          g_yuk_d_c(i,j,t) = 0
          g_yuk_e(i,j,t) = 0
          g_yuk_u(i,j,t) = 0
          g_yuk_d(i,j,t) = 0
          g_yuk_nu(i,j,t) = 0
          
 
       end if
             



          end do
       end do
    end do

!!! SM Yukawas aus SPS1A   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!!! Block yu Q= 4.64241862e+02  
!!!    3  3     8.89738596e-01   # Yt(Q)MSSM DRbar
!!! Block yd Q= 4.64241862e+02  
!!!    3  3     1.39733096e-01   # Yb(Q)MSSM DRbar
!!! Block ye Q= 4.64241862e+02  
!!!    3  3     1.00914889e-01   # Ytau(Q)MSSM DRbar




!!! SFermion mixing matrices (FB)
    mix_sl = 0.0_default
    mix_su = 0.0_default
    mix_sd = 0.0_default
    
    do i = 1,2
       mix_sl(1:2,i,i) = 1.0_default
       mix_su(1:2,i,i) = 1.0_default
       mix_sd(1:2,i,i) = 1.0_default
    end do
        
    mix_sl(3,1,1) = par%ml_11
    mix_sl(3,1,2) = par%ml_12
    mix_sl(3,2,1) = par%ml_21
    mix_sl(3,2,2) = par%ml_22
    mix_su(3,1,1) = par%mt_11
    mix_su(3,1,2) = par%mt_12
    mix_su(3,2,1) = par%mt_21
    mix_su(3,2,2) = par%mt_22
    mix_sd(3,1,1) = par%mb_11
    mix_sd(3,1,2) = par%mb_12
    mix_sd(3,2,1) = par%mb_21
    mix_sd(3,2,2) = par%mb_22 

!!! Leptoquark mixing matrices 

   mix_lq = 0.0_default

   mix_lq(1,1,1) = par%milq1_11
   mix_lq(1,1,2) = par%milq1_12
   mix_lq(1,2,1) = par%milq1_21
   mix_lq(1,2,2) = par%milq1_22   
   mix_lq(2,1,1) = par%milq2_11
   mix_lq(2,1,2) = par%milq2_12
   mix_lq(2,2,1) = par%milq2_21
   mix_lq(2,2,2) = par%milq2_22
   mix_lq(3,1,1) = par%milq3_11
   mix_lq(3,1,2) = par%milq3_12
   mix_lq(3,2,1) = par%milq3_21
   mix_lq(3,2,2) = par%milq3_22    

 end subroutine setup_parameters1
 subroutine setup_parameters2 ()
!!!!!!!!!!!
!!!Neutral Current to Charged Higgs
!!!!!!!!!!!   
   g_Zhmhp = ((gz / 2.0_default) * (1.0_default -  &
        (2.0_default * sin2thw)))
   g_Ahmhp = e
   
!!!!!!!!!!!!!!!
!!!Charged Current to quarks
!!!!!!!!!!!!!!!!!
   do i = 1,3
      do j = 1,3
         g_ccq(i,j) = (gcc * vckm(i,j))
         g_ccq_c(i,j) = (gcc * conjg (vckm(i,j)))
      end do
   end do
!!!!!!!!!!!!!!!
!!!Charged Current to SQuarks
!!!!!!!!!!!!!!!!!
   do gen1 = 1,3
      do gen2 = 1,3
         do sfm1 = 1,2
            do sfm2 = 1,2
               g_wsusd(gen1,gen2,sfm1,sfm2) = &
                    ( - (gcc * 2.0_default * vckm(gen1,gen2) * &
                    conjg(mix_su(gen1,sfm1,1)) * mix_sd(gen2,sfm2,1)))
               
               g_wsusd_c(gen1,gen2,sfm1,sfm2) = conjg(g_wsusd(gen1,gen2,sfm1,sfm2))
            end do
         end do
      end do
   end do
 end subroutine setup_parameters2
 subroutine setup_parameters3 ()
!!!!!!!!!!!!!!!
!!!Neutral current to Sneutrinos
!!!!!!!!!!!!!!!!     
   g_zsnsn(1:3,1,1)  = (gz / 2.0_default)
!!!!!!!!!!!!!!!
!!!Neutral current to SLeptons
!!!!!!!!!!!!!!!!     
   
   do gen = 1,3
      do sfm1 = 1,2
         do sfm2 = 1,2
            g_zslsl(gen,sfm1,sfm2) = &
                 ((gz / 2.0_default) * ( mod(sfm1 + sfm2 + 1,2) *  &      
                 (2.0_default) * sin2thw  - (mix_sl(gen,sfm2,1) *  &
                 conjg (mix_sl(gen,sfm1,1)))))
            
            g_zsusu(gen,sfm1,sfm2) = &
                 ((gz / 2.0_default) * ((mix_su(gen,sfm1,1) * & 
                 conjg (mix_su(gen,sfm2,1))) - mod(sfm1 + sfm2 + 1,2) * &
                 ((4.0_default / 3.0_default) * sin2thw) ))
            
            g_zsdsd(gen,sfm1,sfm2) = & 
                 ((gz / 2.0_default) * (( mod(sfm1 + sfm2 + 1,2) * &
                 (2.0_default / 3.0_default) * sin2thw) - & 
                 (mix_sd(gen,sfm2,1) *  conjg (mix_sd(gen,sfm1,1)))))
          end do
       end do
    end do
!!!!!!!!!!!
!!!Neutral current to Leptoquarks
!!!!!!!!!!!

 g_zlqlq=0.0_default
 g_zlqlq(1:3,1,1) = (gz * sin2thw) / 3
 g_zlqlq(1:3,2,2) = (gz * sin2thw) / 3

!!! Neutral LQino Current
 g_nlqc = gz*(1.0_default/3.0_default)*sin2thw

!!!!!!!!!!!
!!!2 Gauge & GluonGauge to Leptoquarks
!!!!!!!!!!!

g_zz_lqlq = (gz**2) * (sin2thw**2) *  2.0 * ((q_down)**2)
g_zA_lqlq = e*gz * sin2thw * 2.0 * ((q_down)**2)
g_AA_lqlq = (e**2) * 2.0 * ((q_down)**2)
g_zg_lqlq = gz*gs * sin2thw * (q_down)
g_Ag_lqlq = e*gs * (q_down)


!!! 2 Glu 2 LQ

g_gg_lqlq = ((gs/sqrt (2.0_default))**2)

!!! LQ(ino)-Yuk's


   do i = 1,3
       do j = 1,3
          do t = 1,3      
             do l = 1,9
                g_lq_ec_uc(1,1,i,j,t) = g_yuk_d_1(i,t,j)
                g_lq_ec_uc(2,2,i,j,t) = -conjg(g_yuk_d_c(i,j,t))
                g_lq_ec_uc_c(1,2,i,j,t) = -g_yuk_d_c(i,j,t) 
                g_lq_ec_uc_c(2,1,i,j,t) = conjg(g_yuk_d_1(i,t,j)) 
 
                g_yuk_lq_s(1,j,l,t) = g_yuk_n_d(i,j,t)* &
                (1/sqrt(2.0_default))*mix_h0(3*i,l)
                g_yuk_lq_p(1,j,l,t) = g_yuk_n_d(i,j,t)* &
                (1/sqrt(2.0_default))*mix_h0(3*i,l)

             end do
          end do
       end do
   end do



!!! 3 scalar w/ LQ's: F-Terms w/ vev
 
    g_lq_ssu = 0.0_default
    g_lq_ssd = 0.0_default

    do i = 1,3
       do j = 1,3
          do t = 1,3
             do l = 1,3
                do m = 1,3              
                   
                   
                   g_lq_ssu(1,1,2,i,j,t) = g_lq_ssu(1,1,2,i,j,t) - vevs(-2+3*m)*g_yuk_d_1(i,t,l)*conjg (g_yuk_e(l,j,m))
                   g_lq_ssu(2,2,1,i,j,t) = g_lq_ssu(2,2,1,i,j,t) + vevs(-2+3*m)*g_yuk_d_c(i,l,t)*conjg (g_yuk_e(j,l,m))
                   g_lq_ssu(2,1,2,i,j,t) = g_lq_ssu(2,1,2,i,j,t) - vevs(-1+3*m)*g_yuk_d_c(i,j,l)*conjg (g_yuk_u(t,l,m))
                   g_lq_ssu(2,2,2,i,j,t) = g_lq_ssu(2,2,2,i,j,t) + vevs(3*l)*g_yuk_n_d(l,i,m)*conjg (g_yuk_d_1(m,t,j))
                   g_lq_ssu(1,1,1,i,j,t) = g_lq_ssu(1,1,1,i,j,t) + vevs(3*m)*g_yuk_n_d(m,l,i)*conjg (g_yuk_d_c(l,j,t))
                   g_lq_ssu(1,2,1,i,j,t) = g_lq_ssu(1,2,1,i,j,t) - vevs(-1+3*m)*g_yuk_d_1(i,j,l)*conjg (g_yuk_u(l,t,m))
                
                   g_lq_ssd(1,2,i,j,t) = g_lq_ssd(1,2,i,j,t) - vevs(-2+3*l)*g_yuk_d_2(i,j,m)*conjg (g_yuk_nu(m,l,t))
                   g_lq_ssd(2,2,i,j,t) = g_lq_ssd(2,2,i,j,t) - vevs(-2+3*m)*g_yuk_d_c(i,t,l)*conjg (g_yuk_d(j,l,m))
                   g_lq_ssd(1,1,i,j,t) = g_lq_ssd(1,1,i,j,t) + vevs(3*l)*g_yuk_n_d(l,m,i)*conjg (g_yuk_d_c(m,t,j))
                  
                end do
             end do
          end do
       end do
    end do
    
!!! g_lq_ssu/t (sfmLQ,sfmSel,sfmSUp,genLQ,genSel,genSUp)

    g_lq_sst = 0.0_default

    do i = 1,2
       do j = 1,2
          do m = 1,3              
             do n = 1,3
                do o = 1,3
                   
                   g_lq_sst(1,1,1,m,n,o) = g_lq_sst(1,1,1,m,n,o) + g_lq_ssu(1,i,j,m,n,o)*mix_su(o,j,1)*mix_sl(n,i,1)
                   g_lq_sst(1,2,1,m,n,o) = g_lq_sst(1,2,1,m,n,o) + g_lq_ssu(1,i,j,m,n,o)*mix_su(o,j,1)*mix_sl(n,i,2)
                   g_lq_sst(1,1,2,m,n,o) = g_lq_sst(1,1,2,m,n,o) + g_lq_ssu(1,i,j,m,n,o)*mix_su(o,j,2)*mix_sl(n,i,1)
                   g_lq_sst(1,2,2,m,n,o) = g_lq_sst(1,2,2,m,n,o) + g_lq_ssu(1,i,j,m,n,o)*mix_su(o,j,2)*mix_sl(n,i,2)
                   g_lq_sst(2,2,2,m,n,o) = g_lq_sst(2,2,2,m,n,o) + g_lq_ssu(2,i,j,m,n,o)*mix_su(o,j,2)*mix_sl(n,i,2)
                   g_lq_sst(2,2,1,m,n,o) = g_lq_sst(2,2,1,m,n,o) + g_lq_ssu(2,i,j,m,n,o)*mix_su(o,j,1)*mix_sl(n,i,2)
                   g_lq_sst(2,1,2,m,n,o) = g_lq_sst(2,1,2,m,n,o) + g_lq_ssu(2,i,j,m,n,o)*mix_su(o,j,2)*mix_sl(n,i,1)
                   g_lq_sst(2,1,1,m,n,o) = g_lq_sst(2,1,1,m,n,o) + g_lq_ssu(2,i,j,m,n,o)*mix_su(o,j,1)*mix_sl(n,i,1)

                end do
             end do
          end do
       end do
    end do


    g_lq_ssta = 0.0_default

    do i = 1,2
       do m = 1,3              
          do n = 1,3
             do o = 1,3
                             
                g_lq_ssta(1,1,m,n,o) = g_lq_ssta(1,1,m,n,o) + g_lq_ssd(1,i,m,n,o)*mix_sd(o,i,1)
                g_lq_ssta(1,2,m,n,o) = g_lq_ssta(1,2,m,n,o) + g_lq_ssd(1,i,m,n,o)*mix_sd(o,i,2)
                g_lq_ssta(2,1,m,n,o) = g_lq_ssta(2,1,m,n,o) + g_lq_ssd(2,i,m,n,o)*mix_sd(o,i,1)
                g_lq_ssta(2,2,m,n,o) = g_lq_ssta(2,2,m,n,o) + g_lq_ssd(2,i,m,n,o)*mix_sd(o,i,2)
                   
             end do
          end do
       end do
    end do





!!! LQ to Neutralino (slr,sfm,lq,lqino,neu)

    do i = 1,3
      do j = 1,3
        do t = 1,3
          do l = 1,11

                 g_lq_neu(1,2,j,t,l) = g_lq_neu(1,2,j,t,l) + g_yuk_n_d(i,j,t)*mix_neu((i+8),l)
                 g_lq_neu(2,2,j,t,l) = ((sqrt (2.0_default))/3)*(e/costhw)*delta3(j,t)*conjg (mix_neu(1,l))
                 g_lq_neu(1,1,j,t,l) = conjg (g_lq_neu(2,2,j,t,l))
                 g_lq_neu(2,1,j,t,l) = conjg (g_lq_neu(1,2,j,t,l))

          end do
        end do
      end do
    end do

!!! LQ(ino) to Gluino

    do i = 1,3
      do j = 1,3
                 g_lq_gg(1,:,i,j) = (sqrt(2.0_default)/3)*(e/costhw)
      end do
    end do



!!! LQ to SHiggs

    do i = 1,3
      do j = 1,3
        do t = 1,3
          do l = 1,3
            do m = 1,3              
              do n = 1,3
                      g_lq_s(1,1,l,m,n) = &
                      (1/sqrt(2.0_default)) * g_yuk_n_d(i,j,m) * &
                      conjg (g_yuk_n_d(t,j,l)) * ( vevs(3*t) * conjg(mix_h0((3*i),n)) + &
                      vevs(3*i) * mix_h0((3*t),n)) - delta3(l,m)*(((e / costhw)**2)/12) * &
                      ((mix_h0((-1+3*i),n) + conjg(mix_h0((-1+3*i),n))) * vevs(-1+3*i) - &
                       (mix_h0((-2+3*i),n) + conjg(mix_h0((-2+3*i),n))) * vevs(-2+3*i))

                      g_lq_s(2,2,m,j,n) = &
                      (1/sqrt(2.0_default)) * g_yuk_n_d(i,j,t) * &
                      conjg (g_yuk_n_d(l,m,t)) * ( vevs(3*l) * conjg(mix_h0((3*i),n)) + &
                      vevs(3*i) * mix_h0((3*l),n)) - delta3(j,m)*(((e / costhw)**2)/12) * &
                      ((mix_h0((-1+3*i),n) + conjg(mix_h0((-1+3*i),n))) * vevs(-1+3*i) - &
                       (mix_h0((-2+3*i),n) + conjg(mix_h0((-2+3*i),n))) * vevs(-2+3*i))

                      g_lq_s(2,1,j,t,n) = &
                      -(1/sqrt(2.0_default)) * conjg (g_yuk_n_d(i,j,t)) * & 
                      g_yuk_n_h(i,l,m) * ( vevs(-1+3*m) * mix_h0((-2+3*l),n) + &
                      vevs(-2+3*l) * mix_h0((-1+3*m),n) )
                      
                      g_lq_s(1,2,j,t,n) = g_lq_s(2,1,j,t,n)

              end do
            end do
          end do
        end do
      end do
    end do

!!! LQ to PHiggs

    do i = 1,3
      do j = 1,3
        do t = 1,3
          do l = 1,3
            do m = 1,3              
              do n = 1,3
                      g_lq_p(1,1,l,m,n) = &
                      (1/sqrt(2.0_default)) * g_yuk_n_d(i,j,m) * &
                      conjg (g_yuk_n_d(t,j,l)) * ( vevs(3*t) * conjg(mix_A0((3*i),n)) + &
                      vevs(3*i) * mix_A0((3*t),n)) - delta3(l,m)*(((e / costhw)**2)/12) * &
                      ((mix_A0((-1+3*i),n) + conjg(mix_A0((-1+3*i),n))) * vevs(-1+3*i) - &
                       (mix_A0((-2+3*i),n) + conjg(mix_A0((-2+3*i),n))) * vevs(-2+3*i))

                      g_lq_p(2,2,m,j,n) = &
                      (1/sqrt(2.0_default)) * g_yuk_n_d(i,j,t) * &
                      conjg (g_yuk_n_d(l,m,t)) * ( vevs(3*l) * conjg(mix_A0((3*i),n)) + &
                      vevs(3*i) * mix_A0((3*l),n)) - delta3(j,m)*(((e / costhw)**2)/12) * &
                      ((mix_A0((-1+3*i),n) + conjg(mix_A0((-1+3*i),n))) * vevs(-1+3*i) - &
                       (mix_A0((-2+3*i),n) + conjg(mix_A0((-2+3*i),n))) * vevs(-2+3*i))

                      g_lq_p(2,1,j,t,n) = &
                      -(1/sqrt(2.0_default)) * conjg (g_yuk_n_d(i,j,t)) * & 
                      g_yuk_n_h(i,l,m) * ( vevs(-1+3*m) * mix_A0((-2+3*l),n) + &
                      vevs(-2+3*l) * mix_A0((-1+3*m),n) )

                      g_lq_p(1,2,j,t,n) = g_lq_p(2,1,j,t,n) 


              end do
            end do
          end do
        end do
      end do
    end do


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!W to slepton sneutrino 
!!!!!!!!!!!!!!!!!!!!!!!!
    do gen = 1,3
       do i = 1,2
          g_wslsn(gen,i) = (gcc * 2.0_default * mix_sl(gen,i,1))
          g_wslsn_c(gen,i) = conjg(g_wslsn(gen,i))
       end do
    end do
    
!!!!!!!!!!!!!!!!!!!!!!
!!!Neutral current to Neutralinos
!!!!!!!!!!!!!!!!!!!!!!
    g_zneuneu=0
    do neu1 = 1,dimNeu
       g_zneuneu(2,neu1,neu1) =  &
            gz * real (((mix_neu(neu1,4) * conjg (mix_neu(neu1,4))) -  &
            (mix_neu(neu1,3) * conjg (mix_neu(neu1,3))))) / 2.0_default  
       do neu2 = neu1+1,dimNeu
          g_zneuneu(1,neu1,neu2) =  (&
               gz * cmplx (0.0_default, aimag ((mix_neu(neu1,4) * &
               conjg (mix_neu(neu2,4))) - (mix_neu(neu1,3) * &
               conjg (mix_neu(neu2,3)))), kind=default) / 2.0_default)
          g_zneuneu(2,neu1,neu2) =  (&
               gz * real (((mix_neu(neu1,4) * conjg (mix_neu(neu2,4))) -  &
               (mix_neu(neu1,3) * conjg (mix_neu(neu2,3))))) / 2.0_default)  
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!
!!!Neutral current to Charginos
!!!!!!!!!!!!!!!!!!!!!!
    do ch2 = 1,2
       do ch1 = 1,ch2
          g_zchch(1,ch1,ch2) =  &
               (gz*((((1.0_default - (2.0_default * sin2thw)) / 4.0_default) * &
               ((mix_charV(ch1,2) * conjg (mix_charV(ch2,2))) + (&
               conjg (mix_charU(ch1,2)) * mix_charU(ch2,2)))) + ((&
               (costhw**2) / 2.0_default) * ( (mix_charV(ch1,1) * &
               conjg (mix_charV(ch2,1))) + ( conjg (mix_charU(ch1,1)) * &
               mix_charU(ch2,1))))))
          
           g_zchch(2,ch1,ch2) =  &
                (gz*((((1.0_default - (2.0_default * sin2thw)) / 4.0_default) * &
                ((mix_charV(ch1,2) * conjg (mix_charV(ch2,2))) - (&
                conjg (mix_charU(ch1,2)) * mix_charU(ch2,2)))) + ((&
                (costhw**2) / 2.0_default) * ( (mix_charV(ch1,1) * &
                conjg (mix_charV(ch2,1))) - ( conjg (mix_charU(ch1,1)) * &
                mix_charU(ch2,1))))))
        end do
     end do
     g_zchch(1,2,1) = conjg(g_zchch(1,1,2))
     g_zchch(2,2,1) = conjg(g_zchch(2,1,2))
!!!!!!!!!!!!!!!!!!!!!!
!!! Z Z to SFermions 
!!!!!!!!!!!!!!!!!!!!!!
     do gen = 1,3
        do sfm1 = 1,2
           do sfm2 = 1,sfm1
              g_zz_slsl(gen,sfm1,sfm2) = &
                   (((gz**2) / 2.0_default) * ((1.0_default -  &
                   (4.0_default * sin2thw)) * mix_sl(gen,sfm1,1) * &
                   conjg (mix_sl(gen,sfm2,1)) + &
                   mod(sfm1+sfm2+1,2) * ((sin2thw**2) * 4.0_default)))
              g_zz_susu(gen,sfm1,sfm2) = &
                   (((gz**2) / 2.0_default) * ((1.0_default - (sin2thw *  &
                   (8.0_default / 3.0_default))) * mix_su(gen,sfm1,1) *  &
                   conjg (mix_su(gen,sfm2,1)) + &
                   mod(sfm1+sfm2+1,2) * ((sin2thw**2) * (2.0 * q_up)**2)))
              g_zz_sdsd(gen,sfm1,sfm2) = &
                   (((gz**2) / 2.0_default) * ((1.0_default - (sin2thw *  &
                   (4.0_default / 3.0_default))) * mix_sd(gen,sfm1,1) *  &
                   conjg (mix_sd(gen,sfm2,1)) + &
                   mod(sfm1+sfm2+1,2) * ((sin2thw**2) *  (2.0 * q_down)**2)))
           end do
        end do
     end do
     
     g_zz_slsl(:,1,2) = conjg(g_zz_slsl(:,2,1))
     g_zz_susu(:,1,2) = conjg(g_zz_susu(:,2,1))
     g_zz_sdsd(:,1,2) = conjg(g_zz_sdsd(:,2,1))     
     g_zz_snsn(1:3,1,1) = ((gz**2) / 2.0_default)
                
!!!!!!!!!!!!!!!
!!! Photon photon SFerm2
!!!!!!!!!!!!!!!
    g_AA_slsl = (2.0_default * (e**2))
    g_AA_susu = ((8.0_default / 9.0_default) * (e**2))
    g_AA_sdsd = ((2.0_default / 9.0_default) * (e**2))

!!!!!!!!!!!!!!!
!!!W W Sfermion
!!!!!!!!!!!!!!!
    do gen = 1,3
     do sfm1 = 1,2
      do sfm2 = 1, sfm1
       g_ww_slsl(gen,sfm1,sfm2) = &
          (((g**2) / 2.0_default) * mix_sl(gen,sfm1,1) * &
          conjg (mix_sl(gen,sfm2,1)))
       g_ww_sdsd(gen,sfm1,sfm2) = &
          (((g**2) / 2.0_default) * mix_sd(gen,sfm1,1) * &
          conjg (mix_sd(gen,sfm2,1)))
       g_ww_susu(gen,sfm1,sfm2) = &
          (((g**2) / 2.0_default) * mix_su(gen,sfm1,1) * &
          conjg (mix_su(gen,sfm2,1))) 
      end do
     end do
    end do
       g_ww_slsl(:,1,2) = conjg(g_ww_slsl(:,2,1))
       g_ww_susu(:,1,2) = conjg(g_ww_susu(:,2,1))
       g_ww_sdsd(:,1,2) = conjg(g_ww_sdsd(:,2,1))

       g_ww_snsn(1:3,1,1) = ((g**2) / 2.0_default)

!!!!!!!!!!!!!!!
!!!Z Photon Sfermions
!!!!!!!!!!!!!!!
    do gen = 1,3
     do sfm1 = 1,2
      do sfm2 = 1,2
       g_zA_slsl(gen,sfm1,sfm2) = &
         (e * gz  * &
         ((mix_sl(gen,sfm1,1) *  conjg (mix_sl(gen,sfm2,1))) - &          
         (mod(sfm1+sfm2+1,2) * sin2thw * 2.0_default )))       
       g_zA_sdsd(gen,sfm1,sfm2) = &
         (e * gz * (1.0_default / 3.0_default) * &
         ((mix_sd(gen,sfm1,1) *  conjg (mix_sd(gen,sfm2,1))) - &          
         (mod(sfm1+sfm2+1,2) * sin2thw * (2.0_default / 3.0_default))))
       g_zA_susu(gen,sfm1,sfm2) = &
         (e * gz * (2.0_default / 3.0_default) * &
         ((mix_su(gen,sfm1,1) * conjg (mix_su(gen,sfm2,1))) - &
         (mod(sfm1+sfm2+1,2) * sin2thw * (4.0_default / 3.0_default))))
      end do
     end do
    end do
 !      g_zA_slsl(:,1,2) = conjg(g_ww_slsl(:,2,1))
  !     g_zA_susu(:,1,2) = conjg(g_ww_susu(:,2,1))
   !    g_zA_sdsd(:,1,2) = conjg(g_ww_sdsd(:,2,1))

!!!!!!!!!!!!!!!
!!!W Photon SLeptons
!!!!!!!!!!!!!!!
    do gen = 1,3  
     do sfm1 = 1,2
        g_wA_slsn(gen,sfm1) = &
           ( - (e * 2.0_default * gcc * mix_sl(gen,sfm1,1)))          
     end do
    end do
        g_wA_slsn_c = conjg(g_wA_slsn)

!!!!!!!!!!!!!!!
!!!W Z SLeptons
!!!!!!!!!!!!!!!
    do gen = 1,3  
     do sfm1 = 1,2
      g_wz_slsn(gen,sfm1) = &
         gcc * gz * 2.0_default * sin2thw * mix_sl(gen,sfm1,1)          
     end do
    end do
        g_wz_slsn_c = conjg(g_wz_slsn)
!!!!!!!!!!!!!!!
!!!W Photon SQuarks
!!!!!!!!!!!!!!!
    do gen1 = 1,3
     do gen2 = 1,3
      do sfm1 = 1,2
       do sfm2 = 1,2
         g_wA_susd(gen1,gen2,sfm1,sfm2) = &
            ((gcc * e * (2.0_default / 3.0_default) * vckm(gen1,gen2) *  &
            conjg (mix_su(gen1,sfm1,1)) * mix_sd(gen2,sfm2,1)))
       end do
      end do
     end do
    end do
         g_wA_susd_c = conjg(g_wA_susd)

!!!!!!!!!!!!!!!
!!!W Z SQuarks
!!!!!!!!!!!!!!!
    do gen1 = 1,3
     do gen2 = 1,3
      do sfm1 = 1,2
       do sfm2 = 1,2
         g_wz_susd(gen1,gen2,sfm1,sfm2) = &
           ( - (gcc * gz * (2.0_default / 3.0_default) * &
           sin2thw * vckm(gen1,gen2) *  &
           conjg (mix_su(gen1,sfm1,1)) * mix_sd(gen2,sfm2,1)))
       end do
      end do
     end do
    end do
         g_wz_susd_c = conjg(g_wz_susd)

!!!!!!!!!!!!!!!
!!!Gluon W SQuarks
!!!!!!!!!!!!!!!
         do gen1 = 1,3
            do gen2 = 1,3
               do sfm1 = 1,2
                  do sfm2 = 1,2
                     g_gw_susd(gen1,gen2,sfm1,sfm2) = &
                          (g * gs * sqrt (2.0_default) * vckm(gen1,gen2) *  &
                          conjg (mix_su(gen1,sfm1,1)) * mix_sd(gen2,sfm2,1))
                  end do
               end do
            end do
         end do
         g_gw_susd_c = conjg(g_gw_susd)
!!!!!!!!!!!!!!!
!!!Gluon Z SQuarks
!!!!!!!!!!!!!!!
         do gen = 1,3
            do sfm1 = 1,2
               do sfm2 = 1,2
                  g_gz_susu(gen,sfm1,sfm2) = &
                       (gz * gs * (((1.0_default / 2.0_default) *  &
                       (mix_su(gen,sfm1,1) * conjg (mix_su(gen,sfm2,1)))) - &
                       mod(sfm1+sfm2+1,2)* (sin2thw *(q_up))))
                  g_gz_sdsd(gen,sfm1,sfm2) = &
                       (-(gz * gs * (((1.0_default / 2.0_default) *  &
                       (mix_sd(gen,sfm1,1) * conjg (mix_sd(gen,sfm2,1)))) + &
                       mod(sfm1+sfm2+1,2)* (sin2thw *(q_down)))))
               end do
            end do
         end do
!!!!!!!!!!!!!!!
!!!Glu Photon Squarks
!!!!!!!!!!!!!!!!
    g_gA_sqsq = 2.0_default * e * gs / 3.0_default
!!!!!!!!!!!!!!!
!!!Glu Glu Squarks
!!!!!!!!!!!!!!!!
    g_gg_sqsq = (gssq**2)    
!!!!!!!!!!!
!!!W to Chargino-Neuralino
!!!!!!!!!!!!!
    do ch1 = 1,2
       do neu1 = 1,dimNeu
          g_cwn(1,ch1,neu1) = gcc*((conjg (mix_neu(neu1,2)) * mix_charV(ch1,1) * sqrt (2.0)) - ( &
               conjg (mix_neu(neu1,4)) * mix_charV(ch1,2)))  
          g_cwn(2,ch1,neu1) = gcc*((mix_neu(neu1,2) * conjg (mix_charU(ch1,1)) * sqrt (2.0)) + (mix_neu(neu1,3) *  &
               conjg (mix_charU(ch1,2))))
          g_nwc(1,ch1,neu1) = conjg(g_cwn(1,ch1,neu1))
          g_nwc(2,ch1,neu1) = conjg(g_cwn(2,ch1,neu1))
       end do
    end do
 end subroutine setup_parameters3
subroutine setup_parameters4 ()
!!!!!!!!!!!!
!!!!Scalar Higgs coupling to chiral SQuarks a la Franke Fraas
!!!!!!!!!!!!
 do shiggs = 1,dimh0
   do gen = 1,3
   g_h0_suLsuL(shiggs, gen) = &
      (-g * (mass(2*gen) **2) / ( mass(24) * sinbe ) * mix_h0(shiggs,2) + &
      (g / 2.0_default ) * ( mass(23) / costhw) * (1.0_default - &
       2.0_default * q_up * (sinthw ** 2)) * (mix_h0(shiggs,2) * &
       sinbe - mix_h0(shiggs,1) * cosbe))
    g_h0_suRsuR(shiggs, gen) = &
       (-g * (mass(2*gen) **2) / ( mass(24) *sinbe ) * mix_h0(shiggs,2) + &
       g * mass(24) * q_up * ((sinthw / costhw)** 2) * (mix_h0(shiggs,2) * &
       sinbe -  mix_h0(shiggs,1) * cosbe))
    g_h0_suLsuR(shiggs, gen) = &
      (-g * mass(2*gen)  / (2.0_default * mass(24) * sinbe ) * &
      ( lambda *  vev * cosbe * mix_h0(shiggs,3)  - mu * &
      mix_h0(shiggs,1) + au(gen) * mix_h0(shiggs,2) ) )
    g_h0_suRsuL(shiggs, gen) = conjg (g_h0_suLsuR(shiggs, gen))

    g_h0_sdLsdL(shiggs, gen) = &
      (-g * (mass(2*gen-1) **2) / ( mass(24) * cosbe ) * mix_h0(shiggs,1) - &
      (g / 2.0_default ) * ( mass(23) / costhw) * (1.0_default + &
      2.0_default * q_down * (sinthw ** 2)) * (mix_h0(shiggs,2) * &
      sinbe - mix_h0(shiggs,1) * cosbe))
    g_h0_sdRsdR(shiggs, gen) = &
      (-g * (mass(2*gen-1) **2) / ( mass(24) * cosbe ) * mix_h0(shiggs,1) + &
      g * mass(24) * q_down * ((sinthw / costhw)** 2) * (mix_h0(shiggs,2) * &
      sinbe - mix_h0(shiggs,1) * cosbe))
    g_h0_sdLsdR(shiggs, gen) = &
      (-g * mass(2*gen-1) / (2.0_default * mass(24) * cosbe ) * & 
      (( lambda *  vev * sinbe * mix_h0(shiggs,3) - mu * &
      mix_h0(shiggs,2)) + ad(gen) * mix_h0(shiggs,1) ))
    g_h0_sdRsdL(shiggs, gen) = conjg (g_h0_sdLsdR(shiggs, gen))


    g_h0_slLslL(shiggs, gen) = &
      (-g * (mass(2*gen+9) **2) / ( mass(24) * cosbe ) * mix_h0(shiggs,1) - &
      (g / 2.0_default ) * ( mass(23) / costhw) * (1.0_default + &
      2.0_default * q_lep * (sinthw ** 2)) * (mix_h0(shiggs,2) * sinbe -&
       mix_h0(shiggs,1) * cosbe))
    g_h0_slRslR(shiggs, gen) = &
      (-g * (mass(2*gen+9) **2) / ( mass(24) * cosbe ) * mix_h0(shiggs,1) + &
      g * mass(24) * q_lep * ((sinthw / costhw)** 2) * (mix_h0(shiggs,2) * &
      sinbe - mix_h0(shiggs,1) * cosbe))
    g_h0_slLslR(shiggs, gen) = &
      (-g * mass(2*gen+9) / (2.0_default * mass(24) * cosbe ) * & 
      ( lambda *  vev * sinbe * mix_h0(shiggs,3) - mu * mix_h0(shiggs,2) + al(gen) * mix_h0(shiggs,1) ))
    g_h0_slRslL(shiggs, gen) = conjg (g_h0_slLslR(shiggs, gen))

   end do
 end do
   do shiggs = 1,3
    g_h0_snLsnL(shiggs,:) = &
      (gz * mass(23) * (1.0_default / 2.0_default) * &
      (mix_h0(shiggs,2) * sinbe - mix_h0(shiggs,1) * cosbe ))
   end do

!!!!!!!!!!!!
!!!!Axial Higgs coupling to chiral SQuarks a la Franke Fraas
!!!!!!!!!!!!
 do phiggs = 1,dimA0
   do gen = 1,3

    g_A0_suLsuR(phiggs, gen) = (- imago * g * mass(2*gen)  / (2.0_default * mass(24) * sinbe ) * & 
    ( lambda *  vev * cosbe * mix_A0(phiggs,3) - mu * mix_A0(phiggs,1)) - au(gen) * mix_A0(phiggs,2) )
    g_A0_suRsuL(phiggs, gen) = conjg (g_A0_suLsuR(phiggs, gen))
    
    g_A0_sdLsdR(phiggs, gen) = (- imago * g * mass(2*gen-1) / (2.0_default * mass(24) * cosbe ) * & 
    ( lambda *  vev * sinbe * mix_A0(phiggs,3) -mu * mix_A0(phiggs,2)) - ad(gen) * mix_A0(phiggs,1) )
    g_A0_sdRsdL(phiggs, gen) = conjg (g_A0_sdLsdR(phiggs, gen))

    g_A0_slLslR(phiggs, gen) = (- imago * g * mass(2*gen+9) / (2.0_default * mass(24) * cosbe ) * & 
    ( lambda *  vev * sinbe * mix_A0(phiggs,3) - mu * mix_A0(phiggs,2)) - al(gen) * mix_A0(phiggs,1) )
    g_A0_slRslL(phiggs, gen) = conjg (g_A0_slLslR(phiggs, gen))

   end do 
 end do
!!!!!!!!!!!!!!!!!!
!!!!Charged Higgs SLepton Sneutrino (L/R)
!!!!!!!!!!!!!!!!!!
 do gen = 1,3 
   g_hp_slLsnL(gen) = &
      ((g / (sqrt (2.0_default) * mass(24))) * (( &
      (mass(9 + 2*gen)**2) * tanb) - ((mass(24)**2) * sin2be)))
   g_hp_slRsnL(gen) = &
      (sqrt (2.0_default) * ((g * mass(9 + 2*gen) * ((conjg ( &
      al(gen)) * sinbe) + (mu * cosbe))) /  &
      (2.0_default * mass(24) * cosbe)))
 end do

 
end subroutine setup_parameters4
subroutine setup_parameters5 ()
!!!!!!!!!!!!!!!!!!!
!!!Charged Higgs to squarks (gauge) (gensu,gensd)
!!!!!!!!!!!!!!!!!!!    
  do gen1 = 1,3
     do gen2 = 1,3
        g_hp_suLsdL(gen1,gen2) = &
             ((g / (sqrt (2.0_default) * mass(24))) * vckm(gen1,gen2) * &
             ((-( (mass(24)**2) * sin2be)) + (((mass(2*gen2-1)**2) * tanb + & 
             ((mass(2*gen1)**2) / tanb)))))
        g_hp_suRsdR(gen1,gen2) = &
             ((sqrt (2.0_default) * g * mass(2*gen1) * mass(2*gen2-1) * &
             vckm(gen1,gen2)) / (mass(24) * sin2be))
        g_hp_suLsdR(gen1,gen2) = &
             (((g * mass(2*gen2-1)) / ( sqrt (2.0_default) * mass(24))) * &   
             vckm(gen1,gen2) * (mu + ( conjg(ad(gen2)) * tanb)))
        g_hp_suRsdL(gen1,gen2) = &
             (((g * mass(2*gen1)) / ( sqrt (2.0_default) * mass(24))) * &
             vckm(gen1,gen2) * (conjg(mu) + (au(gen1) / tanb)))
     end do
  end do
  
end subroutine setup_parameters5
!!!!!!!!!!!!
!!!gauge => masseigenstates
!!!!!!!!!!!!
 subroutine setup_parameters6 ()
!!!!!!!!!!!!
!!!Sfermions to Scalar Higgs 
!!!!!!!!!!!!
   do shiggs = 1,dimh0
      do gen = 1,3
!!!SUp
         g_h0_su1su1(shiggs,gen) = &
              (conjg (mix_su(gen,1,1)) * mix_su(gen,1,1) * g_h0_suLsuL(shiggs,gen) + &
              conjg (mix_su(gen,1,2)) * mix_su(gen,1,2) * g_h0_suRsuR(shiggs,gen) + &
              conjg (mix_su(gen,1,1)) * mix_su(gen,1,2) * g_h0_suLsuR(shiggs,gen) + &
              conjg (mix_su(gen,1,2)) * mix_su(gen,1,1) * g_h0_suRsuL(shiggs,gen))
         
         g_h0_su2su2(shiggs,gen) = &
              (conjg (mix_su(gen,2,1)) * mix_su(gen,2,1) * g_h0_suLsuL(shiggs,gen) + &
              conjg (mix_su(gen,2,2)) * mix_su(gen,2,2) * g_h0_suRsuR(shiggs,gen) + &
              conjg (mix_su(gen,2,1)) * mix_su(gen,2,2) * g_h0_suLsuR(shiggs,gen) + &
              conjg (mix_su(gen,2,2)) * mix_su(gen,2,1) * g_h0_suRsuL(shiggs,gen))
         
         g_h0_su1su2(shiggs,gen) = &
              (conjg (mix_su(gen,1,1)) * mix_su(gen,2,1) * g_h0_suLsuL(shiggs,gen) + &
              conjg (mix_su(gen,1,2)) * mix_su(gen,2,2) * g_h0_suRsuR(shiggs,gen) + &
              conjg (mix_su(gen,1,1)) * mix_su(gen,2,2) * g_h0_suLsuR(shiggs,gen) + &
              conjg (mix_su(gen,1,2)) * mix_su(gen,2,1) * g_h0_suRsuL(shiggs,gen))
         
         g_h0_su2su1(shiggs,gen) = conjg (g_h0_su1su2(shiggs,gen)) 
!!!SDown
         g_h0_sd1sd1(shiggs,gen) = &
              (conjg (mix_sd(gen,1,1)) * mix_sd(gen,1,1) * g_h0_sdLsdL(shiggs,gen) + &
              conjg (mix_sd(gen,1,2)) * mix_sd(gen,1,2) * g_h0_sdRsdR(shiggs,gen) + &
              conjg (mix_sd(gen,1,1)) * mix_sd(gen,1,2) * g_h0_sdLsdR(shiggs,gen) + &
              conjg (mix_sd(gen,1,2)) * mix_sd(gen,1,1) * g_h0_sdRsdL(shiggs,gen))
         
         g_h0_sd2sd2(shiggs,gen) = &
              (conjg (mix_sd(gen,2,1)) * mix_sd(gen,2,1) * g_h0_sdLsdL(shiggs,gen) + &
              conjg (mix_sd(gen,2,2)) * mix_sd(gen,2,2) * g_h0_sdRsdR(shiggs,gen) + &
              conjg (mix_sd(gen,2,1)) * mix_sd(gen,2,2) * g_h0_sdLsdR(shiggs,gen) + &
              conjg (mix_sd(gen,2,2)) * mix_sd(gen,2,1) * g_h0_sdRsdL(shiggs,gen))
         
         g_h0_sd1sd2(shiggs,gen) = &
              (conjg(mix_sd(gen,1,1)) * mix_sd(gen,2,1) * g_h0_sdLsdL(shiggs,gen) + &
              conjg(mix_sd(gen,1,2)) * mix_sd(gen,2,2) * g_h0_sdRsdR(shiggs,gen) + &
              conjg(mix_sd(gen,1,1)) * mix_sd(gen,2,2) * g_h0_sdLsdR(shiggs,gen) + &
              conjg(mix_sd(gen,1,2)) * mix_sd(gen,2,1) * g_h0_sdRsdL(shiggs,gen))
         
         g_h0_sd2sd1(shiggs,gen) = conjg(g_h0_sd1sd2(shiggs,gen)) 
!!!SLep
        g_h0_sl1sl1(shiggs,gen) = &
             (conjg (mix_sl(gen,1,1)) * mix_sl(gen,1,1) * g_h0_slLslL(shiggs,gen) + &
             conjg (mix_sl(gen,1,2)) * mix_sl(gen,1,2) * g_h0_slRslR(shiggs,gen) + &
             conjg (mix_sl(gen,1,1)) * mix_sl(gen,1,2) * g_h0_slLslR(shiggs,gen) + &
             conjg (mix_sl(gen,1,2)) * mix_sl(gen,1,1) * g_h0_slRslL(shiggs,gen))
        
        g_h0_sl2sl2(shiggs,gen) = &
             (conjg (mix_sl(gen,2,1)) * mix_sl(gen,2,1) * g_h0_slLslL(shiggs,gen) + &
             conjg (mix_sl(gen,2,2)) * mix_sl(gen,2,2) * g_h0_slRslR(shiggs,gen) + &
             conjg (mix_sl(gen,2,1)) * mix_sl(gen,2,2) * g_h0_slLslR(shiggs,gen) + &
             conjg (mix_sl(gen,2,2)) * mix_sl(gen,2,1) * g_h0_slRslL(shiggs,gen))
        
        g_h0_sl1sl2(shiggs,gen) = &
             (conjg (mix_sl(gen,1,1)) * mix_sl(gen,2,1) * g_h0_slLslL(shiggs,gen) + &
             conjg (mix_sl(gen,1,2)) * mix_sl(gen,2,2) * g_h0_slRslR(shiggs,gen) + &
             conjg (mix_sl(gen,1,1)) * mix_sl(gen,2,2) * g_h0_slLslR(shiggs,gen) + &
             conjg (mix_sl(gen,1,2)) * mix_sl(gen,2,1) * g_h0_slRslL(shiggs,gen))
        
        g_h0_sl2sl1(shiggs,gen) = conjg (g_h0_sl1sl2(shiggs,gen)) 
!!!SNeutrino
        g_h0_sn1sn1(shiggs,gen) =  g_h0_snLsnL(shiggs,gen)       
        
     end do
  end do
  
!!!!!!!!!!!!
!!!Sfermions to Axial Higgs 
!!!!!!!!!!!! 
  do phiggs = 1,dimA0
     do gen = 1,3  
!!!SUp
        g_A0_su1su1(phiggs, gen) = &
             (conjg (mix_su(gen,1,1)) * mix_su(gen,2,1) * g_A0_suLsuR(phiggs, gen) + &
             conjg (mix_su(gen,2,1)) * mix_su(gen,1,1) * g_A0_suRsuL(phiggs, gen))
        
        g_A0_su2su2(phiggs, gen) = &
             (conjg (mix_su(gen,1,2)) * mix_su(gen,2,2) * g_A0_suLsuR(phiggs, gen) + &
             conjg (mix_su(gen,2,2)) * mix_su(gen,1,2) * g_A0_suRsuL(phiggs, gen))
        
        g_A0_su1su2(phiggs, gen) = &
             (conjg (mix_su(gen,1,1)) * mix_su(gen,2,2) * g_A0_suLsuR(phiggs, gen) + &
             conjg (mix_su(gen,2,1)) * mix_su(gen,1,2) * g_A0_suRsuL(phiggs, gen))
        
        g_A0_su2su1(phiggs, gen) = conjg (g_A0_su1su2(phiggs, gen)) 
!!!SDown
        g_A0_sd1sd1(phiggs, gen) = &
             (conjg (mix_sd(gen,1,1)) * mix_sd(gen,2,1) * g_A0_sdLsdR(phiggs, gen) + &
             conjg (mix_sd(gen,2,1)) * mix_sd(gen,1,1) * g_A0_sdRsdL(phiggs, gen))
        
        g_A0_sd2sd2(phiggs, gen) = &
             (conjg (mix_sd(gen,1,2)) * mix_sd(gen,2,2) * g_A0_sdLsdR(phiggs, gen) + &
             conjg (mix_sd(gen,2,2)) * mix_sd(gen,1,2) * g_A0_sdRsdL(phiggs, gen))
        
        g_A0_sd1sd2(phiggs, gen) = &
             (conjg (mix_sd(gen,1,1)) * mix_sd(gen,2,2) * g_A0_sdLsdR(phiggs, gen) + &
             conjg (mix_sd(gen,2,1)) * mix_sd(gen,1,2) * g_A0_sdRsdL(phiggs, gen))
        
        g_A0_sd2sd1(phiggs, gen) = conjg (g_A0_sd1sd2(phiggs, gen)) 
!!!SLep
        g_A0_sl1sl1(phiggs, gen) = &
             (conjg (mix_sl(gen,1,1)) * mix_sl(gen,2,1) * g_A0_slLslR(phiggs, gen) + &
             conjg (mix_sl(gen,2,1)) * mix_sl(gen,1,1) * g_A0_slRslL(phiggs, gen))
        
        g_A0_sl2sl2(phiggs, gen) = &
             (conjg (mix_sl(gen,1,2)) * mix_sl(gen,2,2) * g_A0_slLslR(phiggs, gen) + &
             conjg (mix_sl(gen,2,2)) * mix_sl(gen,1,2) * g_A0_slRslL(phiggs, gen))

        g_A0_sl1sl2(phiggs, gen) = &
             (conjg (mix_sl(gen,1,1)) * mix_sl(gen,2,2) * g_A0_slLslR(phiggs, gen) + &
             conjg (mix_sl(gen,2,1)) * mix_sl(gen,1,2) * g_A0_slRslL(phiggs, gen))
        
        g_A0_sl2sl1(phiggs, gen) = conjg (g_A0_sl1sl2(phiggs, gen)) 
        
     end do
  end do
!!!!!!!!!!!!!!!!!
!!!Charged Higgs to SQuarks
!!!!!!!!!!!!!!!!!
  do gen1 = 1,3  
     do gen2 = 1,3  
        g_hp_su1sd1(gen1,gen2) = &
             ((conjg (mix_su(gen1,1,1)) * mix_sd(gen2,1,1) * &
             g_hp_suLsdL(gen1,gen2)) + ( &
             conjg (mix_su(gen1,1,2)) * mix_sd(gen2,1,2) * & 
             g_hp_suRsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,1,1)) * mix_sd(gen2,1,2) * &
             g_hp_suLsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,1,2)) * mix_sd(gen2,1,1) * &
             g_hp_suRsdL(gen1,gen2)))
        g_hp_su1sd2(gen1,gen2) = &
             ((conjg (mix_su(gen1,1,1)) * mix_sd(gen2,2,1) * &
             g_hp_suLsdL(gen1,gen2)) + ( &
             conjg (mix_su(gen1,1,2)) * mix_sd(gen2,2,2) * & 
             g_hp_suRsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,1,1)) * mix_sd(gen2,2,2) * &
             g_hp_suLsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,1,2)) * mix_sd(gen2,2,1) * &
             g_hp_suRsdL(gen1,gen2)))
        g_hp_su2sd1(gen1,gen2) = &
             ((conjg (mix_su(gen1,2,1)) * mix_sd(gen2,1,1) * &
             g_hp_suLsdL(gen1,gen2)) + ( &
             conjg (mix_su(gen1,2,2)) * mix_sd(gen2,1,2) * & 
             g_hp_suRsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,2,1)) * mix_sd(gen2,1,2) * &
             g_hp_suLsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,2,2)) * mix_sd(gen2,1,1) * &
             g_hp_suRsdL(gen1,gen2)))
        g_hp_su2sd2(gen1,gen2) = &
             ((conjg (mix_su(gen1,2,1)) * mix_sd(gen2,2,1) * &
             g_hp_suLsdL(gen1,gen2)) + ( &
             conjg (mix_su(gen1,2,2)) * mix_sd(gen2,2,2) * & 
             g_hp_suRsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,2,1)) * mix_sd(gen2,2,2) * &
             g_hp_suLsdR(gen1,gen2)) + ( &
             conjg (mix_su(gen1,2,2)) * mix_sd(gen2,2,1) * &
             g_hp_suRsdL(gen1,gen2)))
     end do
  end do
  g_hp_su1sd1_c = conjg(g_hp_su1sd1)
  g_hp_su1sd2_c = conjg(g_hp_su1sd2)
  g_hp_su2sd1_c = conjg(g_hp_su2sd1)
  g_hp_su2sd2_c = conjg(g_hp_su2sd2)
  
!!!!!!!!!!!!!!!!!!!!!!!
!!!Scalar Higgs to fermions
!!!!!!!!!!!!!!!!!!!!!!!
  
  do shiggs = 1,dimh0
     do gen = 1,3
        g_yuk_h0_uu(shiggs, gen) = (-g * mass(2 * gen) / &
             ( 2.0_default * mass(24) * sinbe ) * mix_h0(shiggs,2)) 
        
        g_yuk_h0_ll(shiggs, gen) = (-g * mass(9 + 2*gen) / &
             ( 2.0_default * mass(24) * cosbe ) * mix_h0(shiggs,1)) 
        
        g_yuk_h0_dd(shiggs, gen) = (-g * mass(2*gen -1) / &
             ( 2.0_default * mass(24) * cosbe ) * mix_h0(shiggs,1))    
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!
!!!Axial Higgs to fermions
!!!!!!!!!!!!!!!!!!!!!!!
  do phiggs = 1,dimA0
     do gen = 1,3
        g_yuk_A0_uu(phiggs, gen) = ( imago * g * mass(2 * gen) / &
             ( 2.0_default * mass(24) * sinbe ) * mix_A0(phiggs,2)) 
        
        g_yuk_A0_ll(phiggs, gen) = ( imago * g * mass(9 + 2*gen) / &
             ( 2.0_default * mass(24) * cosbe ) * mix_A0(phiggs,1) )
        
        g_yuk_A0_dd(phiggs, gen) = ( imago * g * mass(2*gen -1) / &
             ( 2.0_default * mass(24) * cosbe ) * mix_A0(phiggs,1) )  
     end do
  end do
  
!!!!!!!!!!!!!!!!!!!!!!!!
!!!!Charged Higgs Slepton Sneutrino (mass basis)
!!!!!!!!!!!!!!!!!!!!!!!!
  g_hp_sl1sn1(1:3) = ((conjg (mix_sl(3,1,1)) * g_hp_slLsnL(1:3)) + ( &
       conjg (mix_sl(3,1,2)) * g_hp_slRsnL(1:3)))
  g_hp_sl2sn1(1:3) = ((conjg (mix_sl(3,2,1)) * g_hp_slLsnL(1:3)) + ( &
       conjg (mix_sl(3,2,2)) * g_hp_slRsnL(1:3)))
  g_hp_sl1sn1_c = conjg(g_hp_sl1sn1)
  g_hp_sl2sn1_c = conjg(g_hp_sl2sn1)
  

!!!!!!!!!!!!!!!!!!!!!!!!
!!! Z-shiggs-phiggs
!!!!!!!!!!!!!!!!!!!!!!!!

  do shiggs = 1,dimh0
     do phiggs =1,dimA0
        g_Zh0A0(shiggs,phiggs) = ( - imago * gz / (2.0_default ) * &
             ( mix_h0(shiggs,1) * mix_A0(phiggs,1) - mix_h0(shiggs,2) * mix_A0(phiggs,2) )  ) 
     end do
  end do
  

!!!!!!!!!!!!!!!!!!!!!!!!
!!! W-h+-phiggs
!!!!!!!!!!!!!!!!!!!!!!!!
  do phiggs =1,dimA0   
     g_WhpA0(phiggs) = ( -imago * g / (2.0_default ) * &
          (sinbe * mix_A0(phiggs,1)  + cosbe * mix_A0(phiggs,2) ) )
  end do
  
!!!!!!!!!!!!!!!!!!!!!!!!
!!! Z-Z-shiggs, W-W-shiggs, W-h+-shiggs
!!!!!!!!!!!!!!!!!!!!!!!!
  do shiggs = 1,dimh0   
     g_ZZh0(shiggs) = (  g * mass(23) / costhw * &
        (cosbe * mix_h0(shiggs,1)  + sinbe * mix_h0(shiggs,2) ) )
     g_WWh0(shiggs) = (  g * mass(24)  * &
        (cosbe * mix_h0(shiggs,1)  + sinbe * mix_h0(shiggs,2) ) )
     g_Whph0(shiggs) = ( - g / (2.0_default ) * &
        (sinbe * mix_h0(shiggs,1)  - cosbe * mix_h0(shiggs,2) )) 
  end do
  
!!!!!!!!!!!!!!!!!!!
!!!WW-h0h0, ZZ-h0h0, ZW-hp-h0, AW-hp-h0
!!!!!!!!!!!!!!!!!!!!
  do shiggs1 = 1,dimh0
     do shiggs2 = 1,dimh0
        g_WWh0h0(shiggs1,shiggs2) = ( g**2 / (2.0_default) *(mix_h0(shiggs1,1) * & 
             mix_h0(shiggs2,1) + mix_h0(shiggs1,2) * mix_h0(shiggs2,2))) 
        g_ZZh0h0(shiggs1,shiggs2) = ( g**2 / (2.0_default * costhw**2) *(mix_h0(shiggs1,1) &
             * mix_h0(shiggs2,1) + mix_h0(shiggs1,2) * mix_h0(shiggs2,2))) 
     end do
     g_ZWhph0(shiggs1) = ( g**2 * sinthw**2 / (2.0_default * costhw) *(mix_h0(shiggs1,1) &
          * sinbe - mix_h0(shiggs1,2) * cosbe))
     g_AWhph0(shiggs1) = ( - e * g  / (2.0_default) *(mix_h0(shiggs1,1) * sinbe - & 
          mix_h0(shiggs1,2) * cosbe ))
  end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!WW-A0A0, ZZ-A0A0, ZW-hp-A0, AW-hp-A0
!!!!!!!!!!!!!!!!!!!!!!!!
  do phiggs1 = 1,dimA0
     do phiggs2 = 1,dimA0
        g_WWA0A0(phiggs1,phiggs2) = ( g**2 / (2.0_default) *(mix_A0(phiggs1,1) * &
             mix_A0(phiggs2,1) + mix_A0(phiggs1,2) * mix_A0(phiggs2,2))) 
        g_ZZA0A0(phiggs1,phiggs2) = ( g**2 / (2.0_default * costhw**2) *(mix_A0(phiggs1,1) &
             * mix_A0(phiggs2,1) + mix_A0(phiggs1,2) * mix_A0(phiggs2,2))) 
     end do
     g_ZWhpA0(phiggs1) = (imago * g**2 * sinthw**2 / (2.0_default * costhw) * & 
          (mix_A0(phiggs1,1) * sinbe + mix_A0(phiggs1,2) * cosbe))
     g_AWhpA0(phiggs1) = ( - imago * e * g  / (2.0_default) *(mix_A0(phiggs1,1) * &
          sinbe + mix_A0(phiggs1,2) * cosbe))
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!4-vertex: NEutral Gauge bosons to charged Higgses
!!!!!!!!!!!!!!!!!!!!!!!!!!  
    g_ZZhphm = (((gz**2) / 2.0_default) * &
            (((2.0_default * (costhw**2)) - 1.0_default)**2))
    g_AAhphm = (2.0_default * (e**2))
    g_ZAhphm = (e * gz * ((2.0_default * (costhw**2)) - &    
            1.0_default))
    g_WWhphm = ((g**2) / 2.0_default)
!!!!!!!!!!!!!!!!!!!!!!!!
!!!h0h0h0
!!!!!!!!!!!!!!!!!!!!!!!!
  do shiggs1 = 1,dimh0
   do shiggs2 = 1,shiggs1
    do shiggs3 = 1,shiggs2
     g_h0h0h0(shiggs1,shiggs2,shiggs3) = &
     -(3.0_default / 2.0_default) * (gz**2 / 2.0) *&
        ( cosbe * vev * mix_h0(shiggs1,1) * mix_h0(shiggs2,1) * &
        mix_h0(shiggs3,1) + sinbe * vev * mix_h0(shiggs1,2) * &
        mix_h0(shiggs2,2) * mix_h0(shiggs3,2))  + &
     (gz**2 / (2.0_default * 2.0)- sqrt(2.0) * lambda**2) *&
        cosbe * vev *( mix_h0(shiggs1,1) * mix_h0(shiggs2,2) * &
        mix_h0(shiggs3,2) + mix_h0(shiggs1,2) * mix_h0(shiggs2,1) *&
        mix_h0(shiggs3,2) +  mix_h0(shiggs1,2) * mix_h0(shiggs2,2) *&
        mix_h0(shiggs3,1)) + &
     (gz**2 / (2.0_default * 2.0)- sqrt(2.0) * lambda**2) *&
        sinbe * vev *( mix_h0(shiggs1,1) * mix_h0(shiggs2,1) * &
        mix_h0(shiggs3,2) + mix_h0(shiggs1,1) * mix_h0(shiggs2,2) * &
        mix_h0(shiggs3,1) +  mix_h0(shiggs1,2) * mix_h0(shiggs2,1) * &
        mix_h0(shiggs3,1)) + &
     sqrt(2.0) * ( k * lambda * sinbe - sqrt(2.0) * lambda**2 * cosbe) *&
        vev *( mix_h0(shiggs1,1) * mix_h0(shiggs2,3) * mix_h0(shiggs3,3) +&
        mix_h0(shiggs1,3) * mix_h0(shiggs2,1) * mix_h0(shiggs3,3) + &
        mix_h0(shiggs1,3) * mix_h0(shiggs2,3) * mix_h0(shiggs3,1)) + & 
     sqrt(2.0) * ( k * lambda * cosbe - sqrt(2.0) * lambda**2 * sinbe) *&
        vev *( mix_h0(shiggs1,2) * mix_h0(shiggs2,3) * mix_h0(shiggs3,3) +&
        mix_h0(shiggs1,3) * mix_h0(shiggs2,2) * mix_h0(shiggs3,3) +&
        mix_h0(shiggs1,3) * mix_h0(shiggs2,3) * mix_h0(shiggs3,2)) - &
     sqrt(2.0) * lambda**2 * r * vev *( mix_h0(shiggs1,1) *&
        mix_h0(shiggs2,1) * mix_h0(shiggs3,3) + mix_h0(shiggs1,1) * &
        mix_h0(shiggs2,3) * mix_h0(shiggs3,1) +  mix_h0(shiggs1,3) * &
        mix_h0(shiggs2,1) * mix_h0(shiggs3,1) +  mix_h0(shiggs1,1) * &
        mix_h0(shiggs2,1) * mix_h0(shiggs3,3) + mix_h0(shiggs1,1) * &
        mix_h0(shiggs2,3) * mix_h0(shiggs3,1) +  mix_h0(shiggs1,3) * &
        mix_h0(shiggs2,1) * mix_h0(shiggs3,1)) + &
     lambda* ((A_lambda / sqrt(2.0)) + sqrt(2.0) * k * r * vev ) * &
        ( mix_h0(shiggs1,1) * mix_h0(shiggs2,2) * mix_h0(shiggs3,3) + &
        mix_h0(shiggs1,1) * mix_h0(shiggs2,3) * mix_h0(shiggs3,2) + &
        mix_h0(shiggs1,2) * mix_h0(shiggs2,1) * mix_h0(shiggs3,3) + &
        mix_h0(shiggs1,2) * mix_h0(shiggs2,3) * mix_h0(shiggs3,1) + &
        mix_h0(shiggs1,3) * mix_h0(shiggs2,1) * mix_h0(shiggs3,2) + &
        mix_h0(shiggs1,3) * mix_h0(shiggs2,2) * mix_h0(shiggs3,1)) + &
     (sqrt(2.0) * k * A_k - 6.0_default * sqrt(2.0) * k**2 *r *vev) *&
        (mix_h0(shiggs1,3) * mix_h0(shiggs2,3) * mix_h0(shiggs3,3))
             
    end do
   end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!h0h0A0
!!!!!!!!!!!!!!!!!!!!!!!!
  do shiggs1 = 1,dimh0
   do phiggs1 = 1,dimA0
    do phiggs2 = 1,phiggs1
     g_h0A0A0(shiggs1,phiggs1,phiggs2) = &
       ( -(1.0_default / 2.0_default) * (gz**2 / 2.0) *( cosbe * vev * mix_h0(shiggs1,1) * &
         mix_A0(phiggs1,1) * mix_A0(phiggs2,1) + sinbe * vev * mix_h0(shiggs1,2) * & 
         mix_A0(phiggs1,2) * mix_A0(phiggs2,2)) + (gz**2 / (2.0_default * 2.0)- sqrt(2.0) * &
         lambda**2) * vev *( cosbe * mix_h0(shiggs1,1) * mix_A0(phiggs1,2) * mix_A0(phiggs2,2) + &
         sinbe * mix_h0(shiggs1,2) * mix_A0(phiggs1,1) * mix_A0(phiggs2,1)) - &
         sqrt(2.0) * (  k * lambda * cosbe + sqrt(2.0) * lambda**2 * sinbe) * &
         vev * mix_h0(shiggs1,2) * mix_A0(phiggs1,3) * mix_A0(phiggs2,3) - &
         sqrt(2.0) * ( k * lambda * sinbe + sqrt(2.0) * lambda**2 * cosbe) * &
         vev * mix_h0(shiggs1,1) * mix_A0(phiggs1,3) * mix_A0(phiggs2,3) - & 
         sqrt(2.0) * lambda**2 * r * vev * mix_h0(shiggs1,3) * ( mix_A0(phiggs1,1) * &
         mix_A0(phiggs2,1) + mix_A0(phiggs1,1) * mix_A0(phiggs2,1)) - &
         (sqrt(2.0) * k * A_k + 2.0_default * sqrt(2.0) * k**2 * r *vev ) * & 
         mix_h0(shiggs1,1) * mix_A0(phiggs1,3) * mix_A0(phiggs2,3) + &
         sqrt(2.0) * lambda * k * mix_h0(shiggs1,3) * vev *(cosbe * ( mix_A0(phiggs1,2) * &
         mix_A0(phiggs2,3) + mix_A0(phiggs1,3) * mix_A0(phiggs2,2) ) + sinbe * & 
         ( mix_A0(phiggs1,1) * mix_A0(phiggs2,3) + mix_A0(phiggs1,3) * mix_A0(phiggs2,1) ) ) + &
         lambda* ( sqrt(2.0) * k * r * vev - (A_lambda / sqrt(2.0)) ) * (mix_h0(shiggs1,1) * &
         (mix_A0(phiggs1,2) * mix_A0(phiggs2,3) + mix_A0(phiggs1,3) * mix_A0(phiggs2,2)) + & 
         mix_h0(shiggs1,2) * ( mix_A0(phiggs1,1) * mix_A0(phiggs2,3) + mix_A0(phiggs1,3) * & 
         mix_A0(phiggs2,1) ) ) - lambda* ( sqrt(2.0) * k * r * vev + (A_lambda / sqrt(2.0)) ) * &
         mix_h0(shiggs1,3) *  (mix_A0(phiggs1,1) * mix_A0(phiggs2,2) + mix_A0(phiggs1,2) * mix_A0(phiggs2,1)) )
    end do
   end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!h0h+h-
!!!!!!!!!!!!!!!!!!!!!!!!
  do shiggs1 = 1,dimh0
    g_h0hphm(shiggs1) = - g * mass(24) * (mix_h0(shiggs1,1) * cosbe + mix_h0(shiggs1,2) * sinbe )  - &
    (g * mass(23) / (2.0_default * costhw) ) * (mix_h0(shiggs1,2) * sinbe - mix_h0(shiggs1,1) * cosbe) * cos2be + &
    ( lambda**2 / sqrt(2.0) * vev ) * (mix_h0(shiggs1,2) * sinbe + mix_h0(shiggs1,1) * cosbe) * sin2be - &
    ( lambda / sqrt(2.0) ) * mix_h0(shiggs1,3) * ( (2 * k * r * vev + A_lambda ) * sin2be + 2.0 * lambda *r * vev  ) 
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Scalar Higgs to Neutralinos
!!!!!!!!!!!!!!!!!!!!!!!!!!
  do neu1 = 1,dimNeu
   do neu2 = 1,dimNeu
    do shiggs1 = 1,dimh0
       g_neuneuh0(1,neu1,neu2,shiggs1) = real(&
            mix_h0(shiggs1,2) * (   (g/2.0_default) * (( - sinthw / costhw ) * &
            (mix_neu(neu1,1) * (mix_neu(neu2,4)) + &
            mix_neu(neu2,1) * (mix_neu(neu1,4))) +&
            (mix_neu(neu1,2) * (mix_neu(neu2,4)) + &
            mix_neu(neu2,2) * (mix_neu(neu1,4)))) -&
            (lambda / sqrt(2.0) * & 
            (mix_neu(neu1,5) * (mix_neu(neu2,4)) + &
            mix_neu(neu2,5) * (mix_neu(neu1,4))))    ) - &
            mix_h0(shiggs1,1) * (   (g/2.0_default) * (( - sinthw / costhw ) * &
            (mix_neu(neu1,1) * (mix_neu(neu2,3)) + &
            mix_neu(neu2,1) * (mix_neu(neu1,3))) +&
            (mix_neu(neu1,2) * (mix_neu(neu2,3)) + &
            mix_neu(neu2,2) * (mix_neu(neu1,3)))) +&
            (lambda / sqrt(2.0) * & 
            (mix_neu(neu1,5) * (mix_neu(neu2,3)) + &
            mix_neu(neu2,5) * (mix_neu(neu1,3))))    ) + &
            mix_h0(shiggs1,3) * sqrt(2.0) * k * &
            (mix_neu(neu1,5) * (mix_neu(neu2,5)) + &
            mix_neu(neu2,5) * (mix_neu(neu1,5)))      )
       
       g_neuneuh0(2,neu1,neu2,shiggs1) =&
            cmplx(0.0, aimag(&
            mix_h0(shiggs1,2) * (   (g/2.0_default) * (( - sinthw / costhw ) * &
            (mix_neu(neu1,1) * (mix_neu(neu2,4)) + &
            mix_neu(neu2,1) * (mix_neu(neu1,4))) +&
            (mix_neu(neu1,2) * (mix_neu(neu2,4)) + &
            mix_neu(neu2,2) * (mix_neu(neu1,4)))) -&
            (lambda / sqrt(2.0) * & 
            (mix_neu(neu1,5) * (mix_neu(neu2,4)) + &
            mix_neu(neu2,5) * (mix_neu(neu1,4))))    ) - &
            mix_h0(shiggs1,1) * (   (g/2.0_default) * (( - sinthw / costhw ) * &
            (mix_neu(neu1,1) * (mix_neu(neu2,3)) + &
            mix_neu(neu2,1) * (mix_neu(neu1,3))) +&
            (mix_neu(neu1,2) * (mix_neu(neu2,3)) + &
            mix_neu(neu2,2) * (mix_neu(neu1,3)))) +&
            (lambda / sqrt(2.0) * & 
            (mix_neu(neu1,5) * (mix_neu(neu2,3)) + &
            mix_neu(neu2,5) * (mix_neu(neu1,3))))    ) + &
            mix_h0(shiggs1,3) * sqrt(2.0) * k * &
            (mix_neu(neu1,5) * (mix_neu(neu2,5)) + &
            mix_neu(neu2,5) * (mix_neu(neu1,5)))     ), kind=default )   
    end do
 end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Axial Higgs to Neutralinos
!!!!!!!!!!!!!!!!!!!!!!!!!!

do neu1 = 1,dimNeu
   do neu2 = 1,dimNeu
      do phiggs1 = 1,dimA0
         g_neuneuA0(1,neu1,neu2,phiggs1) = aimag( -(&
              mix_A0(phiggs1,2) * (   (g/2.0_default) * ((- sinthw / costhw ) * &
              (mix_neu(neu1,1) * (mix_neu(neu2,4)) + &
              mix_neu(neu2,1) * (mix_neu(neu1,4))) +&
              (mix_neu(neu1,2) * (mix_neu(neu2,4)) + &
              mix_neu(neu2,2) * (mix_neu(neu1,4)))) -&
              (lambda / sqrt(2.0) * & 
              (mix_neu(neu1,5) * (mix_neu(neu2,4)) + &
              mix_neu(neu2,5) * (mix_neu(neu1,4))))    ) - &
              mix_A0(phiggs1,1) * (   (g/2.0_default) * ((- sinthw / costhw ) * &
              (mix_neu(neu1,1) * (mix_neu(neu2,3)) + &
              mix_neu(neu2,1) * (mix_neu(neu1,3))) +&
              (mix_neu(neu1,2) * (mix_neu(neu2,3)) + &
              mix_neu(neu2,2) * (mix_neu(neu1,3)))) +&
              (lambda / sqrt(2.0) * & 
              (mix_neu(neu1,5) * (mix_neu(neu2,3)) + &
              mix_neu(neu2,5) * (mix_neu(neu1,3))))    ) + &
              mix_A0(phiggs1,3) * sqrt(2.0) * k * &
              (mix_neu(neu1,5) * (mix_neu(neu2,5)) + &
              mix_neu(neu2,5) * (mix_neu(neu1,5)))      ))
         
         g_neuneuA0(2,neu1,neu2,phiggs1) = &
              cmplx(0.0, real((&
              mix_A0(phiggs1,2) * (   (g/2.0_default) * ((- sinthw / costhw ) * &
              (mix_neu(neu1,1) * (mix_neu(neu2,4)) + &
              mix_neu(neu2,1) * (mix_neu(neu1,4))) +&
              (mix_neu(neu1,2) * (mix_neu(neu2,4)) + &
              mix_neu(neu2,2) * (mix_neu(neu1,4)))) -&
              (lambda / sqrt(2.0) * & 
              (mix_neu(neu1,5) * (mix_neu(neu2,4)) + &
              mix_neu(neu2,5) * (mix_neu(neu1,4))))    ) - &
              mix_A0(phiggs1,1) * (   (g/2.0_default) * ((- sinthw / costhw ) * &
              (mix_neu(neu1,1) * (mix_neu(neu2,3)) + &
              mix_neu(neu2,1) * (mix_neu(neu1,3))) +&
              (mix_neu(neu1,2) * (mix_neu(neu2,3)) + &
              mix_neu(neu2,2) * (mix_neu(neu1,3)))) +&
              (lambda / sqrt(2.0) * & 
              (mix_neu(neu1,5) * (mix_neu(neu2,3)) + &
              mix_neu(neu2,5) * (mix_neu(neu1,3))))    ) + &
              mix_A0(phiggs1,3) * sqrt(2.0) * k * &
              (mix_neu(neu1,5) * (mix_neu(neu2,5)) + &
              mix_neu(neu2,5) * (mix_neu(neu1,5)))      )),kind=default)   
      end do
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Scalar Higgs to Charginos
!!!!!!!!!!!!!!!!!!!!!!!!!!
do ch1 = 1,2
   do ch2 = 1,2
      do shiggs1 = 1,dimh0
         g_chchh0(1,ch1,ch2,shiggs1) = - conjg((gcc) * &
              (mix_h0(shiggs1,1) * mix_charU(ch1,2) * mix_charV(ch2,1) + mix_h0(shiggs1,2) * &
              mix_charU(ch1,1) * mix_charV(ch2,2)) - (lambda / sqrt(2.0)) * mix_h0(shiggs1,3) &
              * mix_charU(ch1,2) * mix_charV(ch2,2) )
      end do
   end do
end do
do ch1 = 1,2
   do ch2 = 1,2
      do shiggs1 = 1,dimh0         
         g_chchh0(2,ch1,ch2,shiggs1) = conjg(g_chchh0(1,ch2,ch1,shiggs1))
      end do
   end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Axial Higgs to Charginos
!!!!!!!!!!!!!!!!!!!!!!!!!!
do ch1 = 1,2
   do ch2 = 1,2
      do phiggs1 = 1,dimA0
         g_chchA0(1,ch1,ch2,phiggs1) = imago * conjg((gcc) * &
              (mix_A0(phiggs1,1) * mix_charU(ch1,2) * mix_charV(ch2,1) + mix_A0(phiggs1,2) * &
              mix_charU(ch1,1) * mix_charV(ch2,2)) + (lambda / sqrt(2.0)) * mix_A0(phiggs1,3) &
              * mix_charU(ch1,2) * mix_charV(ch2,2) )
      end do
   end do
end do
do ch1 = 1,2
   do ch2 = 1,2
      do phiggs1 = 1,dimA0
         g_chchA0(2,ch1,ch2,phiggs1) = conjg( g_chchA0(1,ch2,ch1,phiggs1))
      end do
   end do
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Neutralino, H+/- , Chargino
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do ch1 = 1,2
   do neu1 = 1,dimNeu
      
      g_neuhmchar(1,neu1,ch1) = (1.0_default / 2.0_default) * ( &
           g * cosbe * (( conjg(mix_neu(neu1,4)) ) * conjg(mix_charV(ch1,1)) + &
           ( 1.0_default / sqrt(2.0_default) ) * (( sinthw / costhw ) *&
           conjg( mix_neu(neu1,1)) + &
           conjg(mix_neu(neu1,2)) ) * conjg(mix_charV(ch1,2) )) - &
           lambda * sinbe * conjg(mix_neu(neu1,5)) * conjg(mix_charV(ch1,2))) 
      
      g_neuhmchar(2,neu1,ch1) = (1.0_default / 2.0_default) * ( &
           g * sinbe * (( mix_neu(neu1,3) ) * mix_charU(ch1,1) - &
           ( 1.0_default / sqrt(2.0) ) * (( sinthw / costhw ) * &
           mix_neu(neu1,1) + &
           mix_neu(neu1,2) ) * mix_charU(ch1,2) ) - &
           conjg(lambda) * cosbe * mix_neu(neu1,5) * mix_charU(ch1,2)) 
      
   end do
end do
g_neuhmchar_c(1,:,:) = conjg(g_neuhmchar(2,:,:))
g_neuhmchar_c(2,:,:) = conjg(g_neuhmchar(1,:,:))

end subroutine setup_parameters6
subroutine setup_parameters7 ()
!!!!!!!!!!!!!!!!
!!!Neutralino Neutrino Sneutrino
!!!!!!!!!!!!!!!!  
  do gen = 1,3
     do neu1 = 1,dimNeu
        g_yuk_neu_nsn(2,gen,neu1,1) = &
              (- g / (2.0_default * sqrt(2.0_default))) * &
              ((mix_neu(neu1,2) - &
              (sinthw / costhw) * mix_neu(neu1,1)) )   
     end do
  end do
  g_yuk_neu_nsn (2,:,:,2) = 0
  g_yuk_neu_nsn_c(1,:,:,:) = conjg(g_yuk_neu_nsn(2,:,:,:))
  
!!!!!!!!!!!!!!!!
!!!Neutralino Lepton SLepton
!!!!!!!!!!!!!!!!
   do gen = 1,3
      do neu1 = 1,dimNeu
         do sfm1 = 1,2
              g_yuk_neu_lsl(1,gen,neu1,sfm1) = &
                   ( - (gcc * ((2.0_default * ( - q_lep) * conjg (mix_neu(neu1,1) * &
                   (sinthw / costhw) * mix_sl(gen,sfm1,2)) + ((conjg (mix_neu(neu1,3)) * &
                   mass(9 + 2*gen) * mix_sl(gen,sfm1,1)) / (mass(24) * cosbe))))))
              g_yuk_neu_lsl(2,gen,neu1,sfm1) = &
                   (gcc * ((1.0_default * (mix_neu(neu1,2) + (1.0_default *  &
                   (sinthw / costhw) * mix_neu(neu1,1))) * mix_sl(gen,sfm1,1)) - ( &     
                   (mix_neu(neu1,3) * mass(9+2*gen) * mix_sl(gen,sfm1,2)) / (mass(24) * & 
                   cosbe))))
           end do
        end do
     end do
     g_yuk_neu_lsl_c(2,:,:,:) = conjg(g_yuk_neu_lsl(1,:,:,:))
     g_yuk_neu_lsl_c(1,:,:,:) = conjg(g_yuk_neu_lsl(2,:,:,:))     
!!!!!!!!!!!!!!!!
!!!Neutralino Up SUp
!!!!!!!!!!!!!!!!
     do gen = 1,3
        do neu1 = 1,dimNeu
           do sfm1 = 1,2
              g_yuk_neu_usu(1,gen,neu1,sfm1) = &
                   ( - (gcc * ((2.0_default * ( - q_up) * conjg (mix_neu(neu1,1) * &
                   (sinthw / costhw) * mix_su(gen,sfm1,2)) + ((conjg (mix_neu(neu1,4)) *&
                   mass(2*gen) * mix_su(gen,sfm1,1)) / (mass(24) * sinbe))))))
              
              g_yuk_neu_usu(2,gen,neu1,sfm1) = &
                   (gcc * ((-1.0_default * (mix_neu(neu1,2) + &          
                   ((1.0_default/3.0_default )*  &
                   (sinthw / costhw) * mix_neu(neu1,1))) * mix_su(gen,sfm1,1)) - ( &
                   (mix_neu(neu1,4) * mass(2*gen) * mix_su(gen,sfm1,2)) / (mass(24) * & 
                   sinbe))))
           end do
        end do
     end do
     g_yuk_neu_usu_c(2,:,:,:) = conjg(g_yuk_neu_usu(1,:,:,:))
     g_yuk_neu_usu_c(1,:,:,:) = conjg(g_yuk_neu_usu(2,:,:,:))
     
!!!!!!!!!!!!!!!!
!!!Neutralino Down SDown
!!!!!!!!!!!!!!!!


     do gen = 1,3
        do neu1 = 1,dimNeu
           do sfm1 = 1,2
              g_yuk_neu_dsd(1,gen,neu1,sfm1) = &
                   ( - (gcc * ((2.0_default * ( - q_down) * conjg (mix_neu(neu1,1) * &
                   (sinthw / costhw) * mix_sd(gen,sfm1,2)) + ((conjg (mix_neu(neu1,3)) * &
                   mass(2*gen-1) * mix_sd(gen,sfm1,1)) / (mass(24) * cosbe))))))
              g_yuk_neu_dsd(2,gen,neu1,sfm1) = &
                   (gcc * ((1.0_default * (mix_neu(neu1,2) + (- &
                   (1.0_default / 3.0_default ) *  &
                   (sinthw / costhw) * mix_neu(neu1,1))) * mix_sd(gen,sfm1,1)) - ( &
                   (mix_neu(neu1,3) * mass(2*gen -1) * mix_sd(gen,sfm1,2)) / (mass(24) * & 
                   cosbe))))
           end do
        end do
     end do
     g_yuk_neu_dsd_c(2,:,:,:) = conjg(g_yuk_neu_dsd(1,:,:,:))
     g_yuk_neu_dsd_c(1,:,:,:) = conjg(g_yuk_neu_dsd(2,:,:,:))
     
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
!!!!!!!!!!!!!!!!
!!!!!!Charged Higgs to Quarks
!!!!!!!!!!!!!!!!
     do gen1 = 1,3
        do gen2 = 1,3
           g_yuk_hp_ud(1,gen1,gen2) = ((gcc / mass(24)) * vckm(gen1,gen2) * (mass(2*gen1) / tanb))
           g_yuk_hp_ud(2,gen1,gen2) = ((gcc / mass(24)) * vckm(gen1,gen2) * tanb * mass(2*gen2-1))
   
           g_yuk_hm_du(1,gen1,gen2) = conjg(g_yuk_hp_ud(2,gen1,gen2))  
           g_yuk_hm_du(2,gen1,gen2) = conjg(g_yuk_hp_ud(1,gen1,gen2))
        end do
     end do
!!!!!!!!!!!!!!!!
!!!!!!Charged Higgs to Leptons
!!!!!!!!!!!!!!!!
     do gen = 1,3
        g_yuk_hp_ln(gen) = ((gcc / mass(24))  * (mass(9 + 2*(gen)) * tanb))
     end do
   end subroutine setup_parameters7
   subroutine setup_parameters8 ()
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Chargino Lepton Sneutrino
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do gen = 1,3
      do ch1 = 1,2
         g_yuk_char_lsn_c(2,gen,ch1,1) = &
              ( - ((g * mix_charV(ch1,1)) / 2.0_default))
         g_yuk_char_lsn_c(1,gen,ch1,1)=  &
              ((gcc * mass(9 + 2*gen) * conjg (mix_charU(ch1,2))) / &
              (mass(24)  * cosbe))
      end do
   end do
   g_yuk_char_lsn(2,:,:,1) = conjg(g_yuk_char_lsn_c(1,:,:,1))
   g_yuk_char_lsn(1,:,:,1) = conjg(g_yuk_char_lsn_c(2,:,:,1))
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Chargino Neutrino Slepton
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   g_yuk_char_nsl_c(1,:,:,:) = 0 
   
   do gen = 1,3
      do ch1 = 1,2
         do sfm1 = 1,2
            g_yuk_char_nsl_c(2,gen,ch1,sfm1) = &  
                 ((((( - g) / 2.0_default) * mix_charU(ch1,1)) * &
                 conjg (mix_sl(gen,sfm1,1))) +&
                 (((gcc * mass(9+2*gen) * mix_charU(ch1,2)) / &
                 (mass(24) * cosbe)) * conjg (mix_sl(gen,sfm1,2))))
         end do
      end do
   end do
   g_yuk_char_nsl(2,:,:,:) = conjg(g_yuk_char_nsl_c(1,:,:,:))
   g_yuk_char_nsl(1,:,:,:) = conjg(g_yuk_char_nsl_c(2,:,:,:))   
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Chargino Down SUp (slr,sgen,fgen,char,sfm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do gen1 = 1,3
      do gen2 = 1,3
         do ch1 = 1,2
            do sfm1 = 1,2
               g_yuk_char_dsu(1,gen1,gen2,ch1,sfm1) = &  
                    (vckm(gen1,gen2) * gcc * (((conjg (mix_charV(ch1,2)) * &
                    mass(2*gen2) * conjg (mix_su(gen1,sfm1,2))) / &
                    (mass(24) * sinbe)) - (conjg (mix_charV(ch1,1)) *  &
                    sqrt (2.0_default) * conjg (mix_su(gen1,sfm1,1)))))
               g_yuk_char_dsu(2,gen1,gen2,ch1,sfm1) = &  
                    ((vckm(gen1,gen2) * gcc * mix_charU(ch1,2) * mass(2*gen2-1) * &
                    conjg (mix_su(gen1,sfm1,1))) / (mass(24) * cosbe))
            end do
         end do
      end do
   end do
   g_yuk_char_dsu_c(2,:,:,:,:) = conjg(g_yuk_char_dsu(1,:,:,:,:))
   g_yuk_char_dsu_c(1,:,:,:,:) = conjg(g_yuk_char_dsu(2,:,:,:,:))
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!Chargino Up SDown (slr,sgen,fgen,char,sfm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do gen1 = 1,3
      do gen2 = 1,3
         do ch1 = 1,2
            do sfm1 = 1,2
               g_yuk_char_usd(1,gen1,gen2,ch1,sfm1) = &  
                    ((vckm(gen2,gen1) * gcc * mix_charV(ch1,2) * mass(2*gen2) *  &
                    conjg (mix_sd(gen1,sfm1,1))) / (mass(24) * sinbe))
               g_yuk_char_usd(2,gen1,gen2,ch1,sfm1) = &  
                    (vckm(gen2,gen1) * gcc * (((conjg (mix_charU(ch1,2)) * &
                    mass(2*gen1-1) * conjg (mix_sd(gen1,sfm1,2))) / &
                    (mass(24) * cosbe)) - (conjg (mix_charU(ch1,1)) *  &
                    sqrt (2.0_default) * conjg (mix_sd(gen1,sfm1,1)))))
            end do
         end do
      end do
   end do
   g_yuk_char_usd_c(2,:,:,:,:) = conjg(g_yuk_char_usd(1,:,:,:,:))
   g_yuk_char_usd_c(1,:,:,:,:) = conjg(g_yuk_char_usd(2,:,:,:,:))   
 end subroutine setup_parameters8
 subroutine setup_parameters9 ()
!!!!!!!!!!!!!!!
!!!!!!Gluino_Quark_SQuark
!!!!!!!!!!!!!!!
   do gen = 1,3
      do sfm1 = 1,2
         g_yuk_gluino_usu(1,gen,sfm1) = &
              ( - (mix_su(gen,sfm1,2) * (gs / sqrt (2.0_default))))
         g_yuk_gluino_usu(2,gen,sfm1) = & 
              (mix_su(gen,sfm1,1) * (gs / sqrt (2.0_default)))    
         
         g_yuk_gluino_dsd(1,gen,sfm1) = &
              ( - (mix_sd(gen,sfm1,2) * (gs / sqrt (2.0_default))))
         g_yuk_gluino_dsd(2,gen,sfm1) = & 
              (mix_sd(gen,sfm1,1) * (gs / sqrt (2.0_default)))    
      end do
   end do
   g_yuk_gluino_usu_c(1,:,:) =  g_yuk_gluino_usu(2,:,:) 
   g_yuk_gluino_usu_c(2,:,:) =  g_yuk_gluino_usu(1,:,:) 
   
   g_yuk_gluino_dsd_c(1,:,:) =  g_yuk_gluino_dsd(2,:,:) 
   g_yuk_gluino_dsd_c(2,:,:) =  g_yuk_gluino_dsd(1,:,:) 
 end subroutine setup_parameters9
end subroutine import_from_whizard

subroutine model_update_alpha_s (alpha_s)
  real(default), intent(in) :: alpha_s
  integer :: gen, gen1, gen2, sfm1, sfm2
  gs = sqrt (2.0_default * PI * alpha_s)
  igs = cmplx(0.0_default, 1.0_default, kind=default) * gs
  gssq = (gs / sqrt (2.0_default))
  !!!!!!!!!!!!!!!
  !!!Gluon W SQuarks
  !!!!!!!!!!!!!!!
           do gen1 = 1,3
              do gen2 = 1,3
                 do sfm1 = 1,2
                    do sfm2 = 1,2
                       g_gw_susd(gen1,gen2,sfm1,sfm2) = &
                            (g * gs * sqrt (2.0_default) * vckm(gen1,gen2) *  &
                            conjg (mix_su(gen1,sfm1,1)) * mix_sd(gen2,sfm2,1))
                    end do
                 end do
              end do
           end do
           g_gw_susd_c = conjg(g_gw_susd)
  !!!!!!!!!!!!!!!
  !!!Gluon Z SQuarks
  !!!!!!!!!!!!!!!
           do gen = 1,3
              do sfm1 = 1,2
                 do sfm2 = 1,2
                    g_gz_susu(gen,sfm1,sfm2) = &
                         (gz * gs * (((1.0_default / 2.0_default) *  &
                         (mix_su(gen,sfm1,1) * conjg (mix_su(gen,sfm2,1)))) - &
                         mod(sfm1+sfm2+1,2)* (sin2thw *(q_up))))
                    g_gz_sdsd(gen,sfm1,sfm2) = &
                         (-(gz * gs * (((1.0_default / 2.0_default) *  &
                         (mix_sd(gen,sfm1,1) * conjg (mix_sd(gen,sfm2,1)))) + &
                         mod(sfm1+sfm2+1,2)* (sin2thw *(q_down)))))
                 end do
              end do
           end do
  !!!!!!!!!!!!!!!
  !!!Glu Photon Squarks
  !!!!!!!!!!!!!!!!
      g_gA_sqsq = 2.0_default * e * gs / 3.0_default
  !!!!!!!!!!!!!!!
  !!!Glu Glu Squarks
  !!!!!!!!!!!!!!!!
      g_gg_sqsq = (gssq**2)    
  !!!!!!!!!!!!!!!
  !!!!!!Gluino_Quark_SQuark
  !!!!!!!!!!!!!!!
   do gen = 1,3
      do sfm1 = 1,2
         g_yuk_gluino_usu(1,gen,sfm1) = &
              ( - (mix_su(gen,sfm1,2) * (gs / sqrt (2.0_default))))
         g_yuk_gluino_usu(2,gen,sfm1) = & 
              (mix_su(gen,sfm1,1) * (gs / sqrt (2.0_default)))    
         
         g_yuk_gluino_dsd(1,gen,sfm1) = &
              ( - (mix_sd(gen,sfm1,2) * (gs / sqrt (2.0_default))))
         g_yuk_gluino_dsd(2,gen,sfm1) = & 
              (mix_sd(gen,sfm1,1) * (gs / sqrt (2.0_default)))    
      end do
   end do
   g_yuk_gluino_usu_c(1,:,:) =  g_yuk_gluino_usu(2,:,:) 
   g_yuk_gluino_usu_c(2,:,:) =  g_yuk_gluino_usu(1,:,:)    
   g_yuk_gluino_dsd_c(1,:,:) =  g_yuk_gluino_dsd(2,:,:) 
   g_yuk_gluino_dsd_c(2,:,:) =  g_yuk_gluino_dsd(1,:,:) 
end subroutine model_update_alpha_s
end module parameters_psssm

