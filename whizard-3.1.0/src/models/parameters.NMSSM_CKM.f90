! parameters.NMSSM_CKM.f90
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
module parameters_nmssm_ckm
  use kinds
  use constants
  implicit none
  private
  public :: import_from_whizard, model_update_alpha_s
  real(kind=default), dimension(73), save, public :: mass = 0, width = 0
  real(kind=default), parameter, public :: GeV = 1.0_default
  real(kind=default), parameter, public :: MeV = GeV / 1000
  real(kind=default), parameter, public :: keV = MeV / 1000
  real(kind=default), parameter, public :: TeV = GeV * 1000
  real(kind=default), save, public :: &
       alpha = 1.0_default / 137.0359895_default, &
       sin2thw = 0.23124_default
  integer, save, public :: &
       sign1 = +1, sign2 = +1, sign3 = +1, sign4 = +1, sign5 = +1
  real(kind=default), save, public :: &
       sigch1 = +1, sigch2 = +1 
  complex(kind=default), save, private :: vev
  real(kind=default), public, save :: sind = 0._default, & 
       cosd = 1._default, sinckm12 = 0.223_default, &
       sinckm13 = 0.004_default, sinckm23 = 0.04_default, &
       tana = 30._default, tanb = 30._default, as = 0._default
  real(kind=default), public, save :: e
  real(kind=default), public, save :: sin2be, &
       cos2be, cosbe, sinbe, costhw, sinthw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer, parameter :: dimh0 = 3, dimA0 = 2, dimNeu = 5 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!
!!!Higgs Mixing
!!!!!!!!!!!!
  real(kind=default), dimension(dimh0,dimh0), public, save :: mix_h0 = 0
  real(kind=default), dimension(2,3), public, save :: mix_A0 = 0  
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
       gz, g, gs
  complex(kind=default), save, public :: xia = 1, xi0 = 1, xipm = 1
  complex(kind=default), dimension(2), public, save :: gncdwn
  complex(kind=default), dimension(2), public, save :: gncup
  complex(kind=default), dimension(2), public, save :: gnclep
  complex(kind=default), dimension(2), public, save :: gncneu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Couplings!!!
!! the Lorenz structure is encoded in the Caml file in term of
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
!!!Matrix structure coupling neutral Higgs to SQuarks dim(shiggs, gen) final mass eigenstate coupling (rep gh1su1su2_3 ...)
  complex(kind=default), dimension(dimh0,3), public, save :: &
       g_h0_su1su1, g_h0_su2su2, g_h0_su1su2, g_h0_su2su1, &
       g_h0_sd1sd1, g_h0_sd2sd2, g_h0_sd1sd2, g_h0_sd2sd1, &
       g_h0_sl1sl1, g_h0_sl2sl2, g_h0_sl1sl2, g_h0_sl2sl1, & 
       g_h0_sn1sn1
  complex(kind=default), dimension(2,3), public, save :: &
       g_A0_su1su1, g_A0_su2su2, g_A0_su1su2, g_A0_su2su1, &
       g_A0_sd1sd1, g_A0_sd2sd2, g_A0_sd1sd2, g_A0_sd2sd1, &
       g_A0_sl1sl1, g_A0_sl2sl2, g_A0_sl1sl2, g_A0_sl2sl1, &  
       g_A0_sn1sn1, g_A0_sn2sn2, g_A0_sn1sn2, g_A0_sn2sn1
!!! Matrix structure coupling neutral Higgs to SQuarks 
!!! dim(shiggs, gen) gauge eigenstate coupling 
!!! (rep g_h1222susu ...) to be multiplied w/ mix_su
   complex(kind=default), dimension(dimh0,3), public, save :: &
        g_h0_suLsuL, g_h0_suRsuR, g_h0_suLsuR, g_h0_suRsuL, &
        g_h0_sdLsdL, g_h0_sdRsdR, g_h0_sdLsdR, g_h0_sdRsdL, &
        g_h0_slLslL, g_h0_slRslR, g_h0_slLslR, g_h0_slRslL, & 
        g_h0_snLsnL 
   complex(kind=default), dimension(2,3), public, save :: &
        g_A0_suLsuR, g_A0_suRsuL, &
        g_A0_sdLsdR, g_A0_sdRsdL, &
        g_A0_slLslR, g_A0_slRslL 
  integer :: shiggs, phiggs
!!!Matrix structure coupling neutral Higgs to fermions dim (s/phiggs, gen)
  complex(kind=default), dimension(dimh0,3), public, save :: &
       g_yuk_h0_uu, g_yuk_h0_dd, g_yuk_h0_ll   
  complex(kind=default), dimension(2,3), public, save :: &
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
  complex(kind=default), dimension(dimh0,2), public, save :: g_Zh0A0 
!!!W to h+ , (shiggs)
  complex(kind=default), dimension(dimh0), public, save :: g_Whph0
!!!W to h+ , (phiggs)
  complex(kind=default), dimension(2), public, save :: g_WhpA0
!!!ZZ to (shiggs)
  complex(kind=default), dimension(dimh0), public, save :: g_ZZh0
!!!WW to (shiggs)
  complex(kind=default), dimension(dimh0), public, save :: g_WWh0
!!!4- Vertex Higgs Gauge
!!!WW (shiggs1, shiggs2) , ZZ (shiggs1, shiggs2)  
  complex(kind=default), dimension(dimh0,dimh0), public, save :: &
       g_WWh0h0 , g_ZZh0h0
!!!WZ hp (shiggs) , WA hp (shiggs)  
  complex(kind=default), dimension(dimh0), public, save :: &
       g_ZWhph0 , g_AWhph0
!!!WW (phiggs1, phiggs2) , ZZ (phiggs1, phiggs2)  
  complex(kind=default), dimension(2,2), public, save :: &
       g_WWA0A0 , g_ZZA0A0
!!!WZ hp (phiggs) , WA hp (phiggs)  
  complex(kind=default), dimension(2), public, save :: g_ZWhpA0 , g_AWhpA0
  complex(kind=default), public, save :: g_ZZhphm, g_ZAhphm, & 
       g_AAhphm, g_WWhphm   
!!!Triple Higgs couplings
  complex(kind=default), dimension(dimh0,dimh0,dimh0) ,public, save :: g_h0h0h0
  complex(kind=default), dimension(dimh0,dimA0,dimA0) ,public, save :: g_h0A0A0
  complex(kind=default), dimension(dimh0) ,public, save :: g_h0hphm
  integer :: shiggs1, shiggs2, shiggs3, phiggs1, phiggs2
!!!Neutral Higgs to Neutralinos (SLR, neu1 , neu2, s/phiggs)
  complex(kind=default), dimension(2,dimNeu,dimNeu,dimh0) ,public, save :: &
       g_neuneuh0, g_neuneuh0_1  
  complex(kind=default), dimension(2,dimNeu,dimNeu,2) ,public, save :: &
       g_neuneuA0, g_neuneuA0_1
!!!Neutral Higgs to Charginos (SLR, char1 , char2, s/phiggs)
  complex(kind=default), dimension(2,2,2,dimh0) ,public, save :: g_chchh0  
  complex(kind=default), dimension(2,2,2,2) ,public, save :: g_chchA0
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

!!!SFermion mixing Matrix Style
  complex(kind=default), dimension(3,2,2), public, save :: &
    mix_sd, mix_su, mix_sl 
  
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
!!!NMSSM-parameters
  complex(kind=default), public, save ::  mu, lambda, &         
       A_lambda, k, A_k, r
  complex(kind=default), dimension(3), public, save :: al, au, ad
  complex(kind=default), public, save :: eta1, eta2, eta3, eta4, eta5
  
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
  integer :: i,j,sfm1,sfm2

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(177), intent(in) :: par_array
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
       real(default) :: mh0_1
       real(default) :: wh0_1
       real(default) :: mh0_2
       real(default) :: mh0_3
       real(default) :: mA0_1
       real(default) :: mA0_2
       real(default) :: mHpm 
       real(default) :: wh0_2
       real(default) :: wh0_3
       real(default) :: wHpm
       real(default) :: wA0_1
       real(default) :: wA0_2      
       real(default) :: tanb_h
       real(default) :: ls
       real(default) :: a_ls
       real(default) :: ks
       real(default) :: a_ks
       real(default) :: nmu
       real(default) :: mixh0_11
       real(default) :: mixh0_12
       real(default) :: mixh0_13
       real(default) :: mixh0_21
       real(default) :: mixh0_22
       real(default) :: mixh0_23
       real(default) :: mixh0_31
       real(default) :: mixh0_32
       real(default) :: mixh0_33
       real(default) :: mixa0_11
       real(default) :: mixa0_12
       real(default) :: mixa0_13
       real(default) :: mixa0_21
       real(default) :: mixa0_22
       real(default) :: mixa0_23
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
       real(default) :: mch1
       real(default) :: mch2
       real(default) :: mneu1
       real(default) :: mneu2
       real(default) :: mneu3
       real(default) :: mneu4
       real(default) :: mneu5
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
       real(default) :: wneu5
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
       real(default) :: mixn_11
       real(default) :: mixn_12
       real(default) :: mixn_13
       real(default) :: mixn_14
       real(default) :: mixn_15
       real(default) :: mixn_21
       real(default) :: mixn_22
       real(default) :: mixn_23
       real(default) :: mixn_24
       real(default) :: mixn_25
       real(default) :: mixn_31
       real(default) :: mixn_32
       real(default) :: mixn_33
       real(default) :: mixn_34
       real(default) :: mixn_35
       real(default) :: mixn_41
       real(default) :: mixn_42
       real(default) :: mixn_43
       real(default) :: mixn_44
       real(default) :: mixn_45
       real(default) :: mixn_51
       real(default) :: mixn_52
       real(default) :: mixn_53
       real(default) :: mixn_54
       real(default) :: mixn_55
       real(default) :: mu_11
       real(default) :: mu_12
       real(default) :: mu_21
       real(default) :: mu_22
       real(default) :: mv_11
       real(default) :: mv_12
       real(default) :: mv_21
       real(default) :: mv_22
       real(default) :: vckm11
       real(default) :: vckm12
       real(default) :: vckm13
       real(default) :: vckm21
       real(default) :: vckm22
       real(default) :: vckm23
       real(default) :: vckm31
       real(default) :: vckm32
       real(default) :: vckm33
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    real(kind=default) :: qelep, qeup, qedwn, v
    par%gf       = par_array(1)
    par%mZ       = par_array(2)
    par%wZ       = par_array(3)
    par%mW       = par_array(4)
    par%wW       = par_array(5)
    par%me       = par_array(6)
    par%mmu      = par_array(7)
    par%mtau     = par_array(8)
    par%ms       = par_array(9)
    par%mc       = par_array(10)
    par%mb       = par_array(11)
    par%mtop     = par_array(12)
    par%wtop     = par_array(13)
    par%alphas   = par_array(14)
    par%mtype    = par_array(15)
    par%m_zero   = par_array(16)
    par%m_half   = par_array(17)
    par%A0       = par_array(18)
    par%tanb     = par_array(19)
    par%sgn_mu   = par_array(20)
    par%lambda   = par_array(21)
    par%m_mes    = par_array(22)
    par%n5       = par_array(23)
    par%c_grav   = par_array(24)
    par%m_grav   = par_array(25)
    par%ae_33    = par_array(26)
    par%au_33    = par_array(27)
    par%ad_33    = par_array(28)
    par%mh0_1    = par_array(29)
    par%wh0_1    = par_array(30)
    par%mh0_2    = par_array(31)
    par%mh0_3    = par_array(32)
    par%mA0_1    = par_array(33)
    par%mA0_2    = par_array(34)
    par%mHpm     = par_array(35)
    par%wh0_2    = par_array(36)
    par%wh0_3    = par_array(37)
    par%wHpm     = par_array(38)
    par%wA0_1    = par_array(39)
    par%wA0_2    = par_array(40)
    par%tanb_h   = par_array(41)
    par%ls       = par_array(42)
    par%a_ls     = par_array(43)
    par%ks       = par_array(44)
    par%a_ks     = par_array(45)
    par%nmu      = par_array(46)
    par%mixh0_11 = par_array(47)
    par%mixh0_12 = par_array(48)
    par%mixh0_13 = par_array(49)
    par%mixh0_21 = par_array(50)
    par%mixh0_22 = par_array(51)
    par%mixh0_23 = par_array(52)
    par%mixh0_31 = par_array(53)
    par%mixh0_32 = par_array(54)
    par%mixh0_33 = par_array(55)
    par%mixa0_11 = par_array(56)
    par%mixa0_12 = par_array(57)
    par%mixa0_13 = par_array(58)
    par%mixa0_21 = par_array(59)
    par%mixa0_22 = par_array(60)
    par%mixa0_23 = par_array(61)
    par%msu1     = par_array(62)
    par%msd1     = par_array(63)
    par%msc1     = par_array(64)
    par%mss1     = par_array(65)
    par%mstop1   = par_array(66)
    par%msb1     = par_array(67)
    par%msu2     = par_array(68)
    par%msd2     = par_array(69)
    par%msc2     = par_array(70)
    par%mss2     = par_array(71)
    par%mstop2   = par_array(72)
    par%msb2     = par_array(73)
    par%mse1     = par_array(74)
    par%msne     = par_array(75)
    par%msmu1    = par_array(76)
    par%msnmu    = par_array(77)
    par%mstau1   = par_array(78)
    par%msntau   = par_array(79)
    par%mse2     = par_array(80)
    par%msmu2    = par_array(81)
    par%mstau2   = par_array(82)
    par%mgg      = par_array(83)
    par%mch1     = par_array(84)
    par%mch2     = par_array(85)
    par%mneu1    = par_array(86)
    par%mneu2    = par_array(87)
    par%mneu3    = par_array(88)
    par%mneu4    = par_array(89)
    par%mneu5    = par_array(90)
    par%wsu1     = par_array(91)
    par%wsd1     = par_array(92)
    par%wsc1     = par_array(93)
    par%wss1     = par_array(94)
    par%wstop1   = par_array(95)
    par%wsb1     = par_array(96)
    par%wsu2     = par_array(97)
    par%wsd2     = par_array(98)
    par%wsc2     = par_array(99)
    par%wss2     = par_array(100)
    par%wstop2   = par_array(101)
    par%wsb2     = par_array(102)
    par%wse1     = par_array(103)
    par%wsne     = par_array(104)
    par%wsmu1    = par_array(105)
    par%wsnmu    = par_array(106)
    par%wstau1   = par_array(107)
    par%wsntau   = par_array(108)
    par%wse2     = par_array(109)
    par%wsmu2    = par_array(110)
    par%wstau2   = par_array(111)
    par%wgg      = par_array(112)
    par%wch1     = par_array(113)
    par%wch2     = par_array(114)
    par%wneu1    = par_array(115)
    par%wneu2    = par_array(116)
    par%wneu3    = par_array(117)
    par%wneu4    = par_array(118)
    par%wneu5    = par_array(119)
    par%mt_11    = par_array(120)
    par%mt_12    = par_array(121)
    par%mt_21    = par_array(122)
    par%mt_22    = par_array(123)
    par%mb_11    = par_array(124)
    par%mb_12    = par_array(125)
    par%mb_21    = par_array(126)
    par%mb_22    = par_array(127)
    par%ml_11    = par_array(128)
    par%ml_12    = par_array(129)
    par%ml_21    = par_array(130)
    par%ml_22    = par_array(131)
    par%mixn_11  = par_array(132)
    par%mixn_12  = par_array(133)
    par%mixn_13  = par_array(134)
    par%mixn_14  = par_array(135)
    par%mixn_15  = par_array(136)
    par%mixn_21  = par_array(137)
    par%mixn_22  = par_array(138)
    par%mixn_23  = par_array(139)
    par%mixn_24  = par_array(140)
    par%mixn_25  = par_array(141)
    par%mixn_31  = par_array(142)
    par%mixn_32  = par_array(143)
    par%mixn_33  = par_array(144)
    par%mixn_34  = par_array(145)
    par%mixn_35  = par_array(146)
    par%mixn_41  = par_array(147)
    par%mixn_42  = par_array(148)
    par%mixn_43  = par_array(149)
    par%mixn_44  = par_array(150)
    par%mixn_45  = par_array(151)
    par%mixn_51  = par_array(152)
    par%mixn_52  = par_array(153)
    par%mixn_53  = par_array(154)
    par%mixn_54  = par_array(155)
    par%mixn_55  = par_array(156)
    par%mu_11    = par_array(157)
    par%mu_12    = par_array(158)
    par%mu_21    = par_array(159)
    par%mu_22    = par_array(160)
    par%mv_11    = par_array(161)
    par%mv_12    = par_array(162)
    par%mv_21    = par_array(163)
    par%mv_22    = par_array(164)
    par%vckm11   = par_array(165)
    par%vckm12   = par_array(166)
    par%vckm13   = par_array(167)
    par%vckm21   = par_array(168)
    par%vckm22   = par_array(169)
    par%vckm23   = par_array(170)
    par%vckm31   = par_array(171)
    par%vckm32   = par_array(172)
    par%vckm33   = par_array(173)
    par%v        = par_array(174)
    par%cw       = par_array(175)
    par%sw       = par_array(176)
    par%ee       = par_array(177)
    mass(1:73) = 0
    width(1:73) = 0
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
    mass(37) = par%mHpm
    width(37) = par%wHpm
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
    sigch1   = sign (1._default, par%mch1) 
    sigch2   = sign (1._default, par%mch2)
    sign1    = sign (1, int(par%mneu1))
    sign2    = sign (1, int(par%mneu2))
    sign3    = sign (1, int(par%mneu3))
    sign4    = sign (1, int(par%mneu4))
    sign5    = sign (1, int(par%mneu5))
    vckm(1,1)  = par%vckm11
    vckm(1,2)  = par%vckm12
    vckm(1,3)  = par%vckm13
    vckm(2,1)  = par%vckm21
    vckm(2,2)  = par%vckm22
    vckm(2,3)  = par%vckm23
    vckm(3,1)  = par%vckm31
    vckm(3,2)  = par%vckm32
    vckm(3,3)  = par%vckm33
    v = 2 * par%mW * par%sw / par%ee
    e = par%ee
    !!! This should not be the color flow basis
    as = par%alphas
    tanb = par%tanb_h

    !!! Higgs Mixing
    mix_h0(1,1) = par%mixh0_11
    mix_h0(1,2) = par%mixh0_12
    mix_h0(1,3) = par%mixh0_13
    mix_h0(2,1) = par%mixh0_21
    mix_h0(2,2) = par%mixh0_22
    mix_h0(2,3) = par%mixh0_23
    mix_h0(3,1) = par%mixh0_31
    mix_h0(3,2) = par%mixh0_32
    mix_h0(3,3) = par%mixh0_33

    mix_A0 = 0
    mix_A0(1,1) = par%mixa0_11
    mix_A0(1,2) = par%mixa0_12
    mix_A0(1,3) = par%mixa0_13
    mix_A0(2,1) = par%mixa0_21
    mix_A0(2,2) = par%mixa0_22
    mix_A0(2,3) = par%mixa0_23
    select case (sign1)
       case (1)
          eta1 = (1.0_default,0.0_default)
       case (-1)
          eta1 = (0.0_default,1.0_default)
       case default 
          print *, 'sign1', sign1
          stop "parameters_NMSSM: No definite sign neutralino1"
    end select
    select case (sign2)
       case (1)
          eta2 = (1.0_default,0.0_default)
       case (-1)
          eta2 = (0.0_default,1.0_default)
       case default 
          print *, 'sign2', sign2
          stop "parameters_NMSSM: No definite sign neutralino2"
    end select
    select case (sign3)
       case (1)
          eta3 = (1.0_default,0.0_default)
       case (-1)
          eta3 = (0.0_default,1.0_default)
       case default 
          print *, 'sign3', sign3
          stop "parameters_NMSSM: No definite sign neutralino3"
    end select
    select case (sign4)
       case (1)
          eta4 = (1.0_default,0.0_default)
       case (-1)
          eta4 = (0.0_default,1.0_default)
       case default 
          print *, 'sign4', sign4
          stop "parameters_NMSSM: No definite sign neutralino4"
    end select
    select case (sign5)
       case (1)
          eta5 = (1.0_default,0.0_default)
       case (-1)
          eta5 = (0.0_default,1.0_default)
       case default 
          print *, 'sign5', sign5
          stop "parameters_NMSSM: No definite sign neutralino5"
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
    !!!  Further rename the input Block in the SLHA-file from NMIX to NMNMIX 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! lambda=0
    !!! k=0
    !!! A_k=0
    !!! A_lambda=0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! mix_A0 = 0
    !!! mix_A0(1,1) = sinbe
    !!! mix_A0(1,2) = cosbe
    !!! !!!
    !!! mix_h0 = 0
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
    g = (e / sinthw)
    gz = (g / costhw)
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt (2.0_default * PI * par%alphas)
    igs = (imago * gs)
    vev = ((2.0_default * mass(24)) / g)
    q_lep  = (- 1.0_default) 
    q_up   = (2.0_default / 3.0_default)
    q_down = (- 1.0_default / 3.0_default)    
    qlep = - e * qelep   !!! This is the negative particle charge !!! 
    qup = - e * qeup     !!! This is the negative particle charge !!! 
    qdwn = - e * qedwn   !!! This is the negative particle charge !!! 
    qchar = ( - e)       !!! This is the negative particle charge !!! 
    gcc = (g / (2.0_default * sqrt (2.0_default)))
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
    mix_neu(1,1) = eta1 * par%mixn_11
    mix_neu(1,2) = eta1 * par%mixn_12
    mix_neu(1,3) = eta1 * par%mixn_13
    mix_neu(1,4) = eta1 * par%mixn_14
    mix_neu(1,5) = eta1 * par%mixn_15
    mix_neu(2,1) = eta2 * par%mixn_21
    mix_neu(2,2) = eta2 * par%mixn_22
    mix_neu(2,3) = eta2 * par%mixn_23
    mix_neu(2,4) = eta2 * par%mixn_24
    mix_neu(2,5) = eta2 * par%mixn_25
    mix_neu(3,1) = eta3 * par%mixn_31
    mix_neu(3,2) = eta3 * par%mixn_32
    mix_neu(3,3) = eta3 * par%mixn_33
    mix_neu(3,4) = eta3 * par%mixn_34
    mix_neu(3,5) = eta3 * par%mixn_35
    mix_neu(4,1) = eta4 * par%mixn_41
    mix_neu(4,2) = eta4 * par%mixn_42
    mix_neu(4,3) = eta4 * par%mixn_43
    mix_neu(4,4) = eta4 * par%mixn_44
    mix_neu(4,5) = eta4 * par%mixn_45
    mix_neu(5,1) = eta5 * par%mixn_51
    mix_neu(5,2) = eta5 * par%mixn_52
    mix_neu(5,3) = eta5 * par%mixn_53
    mix_neu(5,4) = eta5 * par%mixn_54
    mix_neu(5,5) = eta5 * par%mixn_55
    !!! Checked by JR !!! 
    mix_charU(1,1) = par%mu_11                !!! Rotat. matrix containing phi_R
    mix_charU(1,2) = par%mu_12                !!! Rotat. matrix containing phi_R
    mix_charU(2,1) = par%mu_21                !!! Rotat. matrix containing phi_R
    mix_charU(2,2) = par%mu_22                !!! Rotat. matrix containing phi_R
    mix_charV(1,1) = sigch1 * par%mv_11       !!! Rotat. matrix containing phi_L
    mix_charV(1,2) = sigch1 * par%mv_12       !!! Rotat. matrix containing phi_L
    mix_charV(2,1) = sigch2 * par%mv_21       !!! Rotat. matrix containing phi_L
    mix_charV(2,2) = sigch2 * par%mv_22       !!! Rotat. matrix containing phi_L
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
!!!NMSSM parameters
    mu = par%nmu
    lambda = par%ls   
    A_lambda = par%a_ls  
    k = par%ks
    A_k = par%a_ks
!!! ensuring that r==0 if lambda=0
    if (lambda == 0) then
    r = 0.0_default
    else
    r = mu/(vev * lambda)
    end if
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
       do neu2 = 1,dimNeu
          g_zneuneu(1,neu1,neu2) = &
               (gz/4.0_default) * ( mix_neu(neu1,4) * &
               conjg (mix_neu(neu2,4)) - mix_neu(neu1,3) * &
               conjg (mix_neu(neu2,3)))
          g_zneuneu(2,neu1,neu2) = -conjg(g_zneuneu(1,neu1,neu2))
       end do
    end do
!!!!!!!!!!!!!!!!!!!!!!
!!!Neutral current to Charginos
!!!!!!!!!!!!!!!!!!!!!!
 do ch2 = 1,2
    do ch1 = 1,ch2
       g_zchch(1,ch1,ch2) =  &
            (gz*((((1.0_default - (2.0_default * sin2thw)) / 4.0_default) * &
            ((mix_charV(ch1,2) * conjg (mix_charV(ch2,2))) + &
            (conjg (mix_charU(ch1,2)) * mix_charU(ch2,2)))) + &
            (((costhw**2) / 2.0_default) * ( (mix_charV(ch1,1) * &
            conjg (mix_charV(ch2,1))) + ( conjg (mix_charU(ch1,1)) * &
            mix_charU(ch2,1))))))
       
       g_zchch(2,ch1,ch2) =  &
            (gz*((((1.0_default - (2.0_default * sin2thw)) / 4.0_default) * &
            ((mix_charV(ch1,2) * conjg (mix_charV(ch2,2))) - &
            (conjg (mix_charU(ch1,2)) * mix_charU(ch2,2)))) + &
            (((costhw**2) / 2.0_default) * ( (mix_charV(ch1,1) * &
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
               (e * gz * ((mix_sl(gen,sfm1,1) *  conjg (mix_sl(gen,sfm2,1))) - &          
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
       g_cwn(1,ch1,neu1) = gcc*((conjg (mix_neu(neu1,2)) * &
            mix_charV(ch1,1) * sqrt (2.0)) - &
            (conjg (mix_neu(neu1,4)) * mix_charV(ch1,2)))  
       g_cwn(2,ch1,neu1) = gcc*((mix_neu(neu1,2) * &
            conjg (mix_charU(ch1,1)) * sqrt (2.0)) + &
            (mix_neu(neu1,3) * conjg (mix_charU(ch1,2))))
       g_nwc(1,ch1,neu1) = conjg(g_cwn(1,ch1,neu1))
       g_nwc(2,ch1,neu1) = conjg(g_cwn(2,ch1,neu1))
    end do
 end do
end subroutine setup_parameters3
subroutine setup_parameters4 ()
!!!!!!!!!!!!
!!!!Scalar Higgs coupling to chiral SQuarks,
!!!!SLeptons, and SNeutrinos a la Franke Fraas
!!!!!!!!!!!!
!!!!!!!second sfermion is cc in fr
!!!FB check 07/24/09
 do shiggs = 1,dimh0
    do gen = 1,3
       g_h0_suLsuL(shiggs, gen) = (-g * (mass(2*gen) **2) / ( mass(24) * &
            sinbe ) * mix_h0(shiggs,2) + (g / 2.0_default ) * ( &
            mass(23) / costhw) * (1.0_default - 2.0_default * &
            q_up * (sinthw ** 2)) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe))
       g_h0_suRsuR(shiggs, gen) = (-g * (mass(2*gen) **2) / ( mass(24) &
            *sinbe ) * mix_h0(shiggs,2) + g * mass(24) * q_up * &
            ((sinthw / costhw)** 2) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe))
       g_h0_suLsuR(shiggs, gen) = - (g * mass(2*gen)  / (2.0_default * &
            mass(24) * sinbe )) &
            * (- lambda * (vev / sqrt (2.0_default) )* cosbe * &
            mix_h0(shiggs,3)  - mu * mix_h0(shiggs,1) + au(gen) * &
            mix_h0(shiggs,2))
       g_h0_suRsuL(shiggs, gen) = conjg (g_h0_suLsuR(shiggs, gen))

       g_h0_suLsuL(shiggs, gen) = (-g * (mass(2*gen) **2) / ( mass(24) * &
            sinbe ) * mix_h0(shiggs,2) + (g / 2.0_default ) * ( &
            mass(23) / costhw) * (1.0_default - 2.0_default * &
            q_up * (sinthw ** 2)) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe))

       g_h0_sdLsdL(shiggs, gen) = (-g * (mass(2*gen-1) **2) / ( mass(24) &
            * cosbe ) * mix_h0(shiggs,1) - (g / 2.0_default ) * ( &
            mass(23) / costhw) * (1.0_default + 2.0_default * &
            q_down * (sinthw ** 2)) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe) )
       g_h0_sdRsdR(shiggs, gen) = (-g * (mass(2*gen-1) **2) / ( mass(24) &
            * cosbe ) * mix_h0(shiggs,1) + g * mass(24) * q_down * &
            ((sinthw / costhw)** 2) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe))

       g_h0_sdLsdR(shiggs, gen) =  - (g * mass(2*gen-1) / (2.0_default &
            * mass(24) * cosbe )) * &
            (- lambda * (vev/sqrt(2.0_default))*  sinbe * &
             mix_h0(shiggs,3) - mu * mix_h0(shiggs,2) + ad(gen) * &
             mix_h0(shiggs,1) )
       g_h0_sdRsdL(shiggs, gen) = conjg (g_h0_sdLsdR(shiggs, gen))

       g_h0_slLslL(shiggs, gen) = (-g * (mass(2*gen+9) **2) / ( mass(24) &
            * cosbe ) * mix_h0(shiggs,1) - (g / 2.0_default ) * ( &
            mass(23) / costhw) * (1.0_default + 2.0_default * &
            q_lep * (sinthw ** 2)) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe))
       g_h0_slRslR(shiggs, gen) = (-g * (mass(2*gen+9) **2) / ( mass(24) &
            * cosbe ) * mix_h0(shiggs,1) + g * mass(24) * q_lep * &
            ((sinthw / costhw)** 2) * (mix_h0(shiggs,2) * sinbe - &
            mix_h0(shiggs,1) * cosbe))
       g_h0_slLslR(shiggs, gen) = (-g * mass(2*gen+9) / (2.0_default &
            * mass(24) * cosbe ) * ( - lambda * ( vev/sqrt(2.0)) * sinbe * &
             mix_h0(shiggs,3) - mu * mix_h0(shiggs,2) + al(gen) * &
             mix_h0(shiggs,1) ))
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
!!!FB check 07/24/09
 do phiggs = 1,dimA0
    do gen = 1,3
       g_A0_suLsuR(phiggs, gen) = (- imago * g * mass(2*gen)  / &
            (2.0_default * mass(24) * sinbe )) * &
            ( lambda * ( vev / sqrt (2.0_default)) * &
            cosbe * mix_A0(phiggs,3) + mu * mix_A0(phiggs,1) + &
            au(gen) * mix_A0(phiggs,2) )
       g_A0_suRsuL(phiggs, gen) = conjg (g_A0_suLsuR(phiggs, gen))

       g_A0_sdLsdR(phiggs, gen) = (- imago * g * mass(2*gen-1)) / &
            (2.0_default * mass(24) * cosbe ) * ( lambda *  (vev/sqrt(2.0)) * &
            sinbe * mix_A0(phiggs,3) + mu * mix_A0(phiggs,2) + ad(gen) &
            * mix_A0(phiggs,1) )
       g_A0_sdRsdL(phiggs, gen) = conjg (g_A0_sdLsdR(phiggs, gen))

       g_A0_slLslR(phiggs, gen) = (- imago * g * mass(2*gen+9)) / &
            (2.0_default * mass(24) * cosbe ) * ( lambda * ( vev/sqrt(2.0)) * &
            sinbe * mix_A0(phiggs,3) + mu * mix_A0(phiggs,2) + &
            al(gen) * mix_A0(phiggs,1) )
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
      (sqrt (2.0_default) * ((g * mass(9 + 2*gen) * &
      ((conjg (al(gen)) * sinbe) + (mu * cosbe))) /  &
      (2.0_default * mass(24) * cosbe)))
 end do

 
end subroutine setup_parameters4
subroutine setup_parameters5 ()
!!!!!!!!!!!!!!!!!!!
!!!Charged Higgs to squarks (gauge) (gensu,gensd)
!!!!!!!!!!!!!!!!!!! 
  do gen1 = 1,3
     do gen2 = 1,3
        g_hp_suLsdL(gen1,gen2) = vckm(gen1,gen2) * &
            (g**2 * sinbe * cosbe * vev / sqrt(2.0_default)) &
           *( - 1.0_default  &
              + mass(2*gen2-1)**2 / (2.0_default*mass(24)**2 *cosbe**2) &
              + mass( 2*gen1 )**2 / (2.0_default*mass(24)**2 *sinbe**2) &
            )
        g_hp_suRsdR(gen1,gen2) = vckm(gen1,gen2) * &
            (g**2*vev*mass(2*gen2-1)*mass(2*gen1)/(mass(24)**2*sqrt(8.0))) &
           *(sinbe/cosbe + cosbe/sinbe)
        g_hp_suLsdR(gen1,gen2) = vckm(gen1,gen2) * &
            (g * mass(2*gen2-1) / (sqrt(2.0_default) * mass(24) )) &
           *(ad(gen2) * sinbe / cosbe + mu )
        g_hp_suRsdL(gen1,gen2) = vckm(gen1,gen2) * &
            (g * mass(2*gen1) / (sqrt(2.0_default) * mass(24) )) &
           *(au(gen1) * cosbe / sinbe + mu)
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
         g_h0_su1su1(shiggs,gen) = (conjg (mix_su(gen,1,1)) * &
              mix_su(gen,1,1) * g_h0_suLsuL(shiggs,gen) + &
              conjg(mix_su(gen,1,2)) * mix_su(gen,1,2) * &
              g_h0_suRsuR(shiggs,gen) + conjg (mix_su(gen,1,1)) * &
              mix_su(gen,1,2) * g_h0_suLsuR(shiggs,gen) + &
              conjg(mix_su(gen,1,2)) * mix_su(gen,1,1) * &
              g_h0_suRsuL(shiggs,gen))

         g_h0_su2su2(shiggs,gen) = (conjg (mix_su(gen,2,1)) * &
              mix_su(gen,2,1) * g_h0_suLsuL(shiggs,gen) + &
              conjg(mix_su(gen,2,2)) * mix_su(gen,2,2) * &
              g_h0_suRsuR(shiggs,gen) + conjg (mix_su(gen,2,1)) * &
              mix_su(gen,2,2) * g_h0_suLsuR(shiggs,gen) + &
              conjg(mix_su(gen,2,2)) * mix_su(gen,2,1) * &
              g_h0_suRsuL(shiggs,gen))

         g_h0_su1su2(shiggs,gen) = (conjg (mix_su(gen,1,1)) * &
              mix_su(gen,2,1) * g_h0_suLsuL(shiggs,gen) + &
              conjg(mix_su(gen,1,2)) * mix_su(gen,2,2) * &
              g_h0_suRsuR(shiggs,gen) + conjg (mix_su(gen,1,1)) * &
              mix_su(gen,2,2) * g_h0_suLsuR(shiggs,gen) + &
              conjg(mix_su(gen,1,2)) * mix_su(gen,2,1) * &
              g_h0_suRsuL(shiggs,gen))
         
         g_h0_su2su1(shiggs,gen) = conjg (g_h0_su1su2(shiggs,gen)) 
!!!SDown
         g_h0_sd1sd1(shiggs,gen) = (conjg(mix_sd(gen,1,1)) * &
              mix_sd(gen,1,1) * g_h0_sdLsdL(shiggs,gen) + &
              conjg(mix_sd(gen,1,2)) * mix_sd(gen,1,2) * &
              g_h0_sdRsdR(shiggs,gen) + conjg(mix_sd(gen,1,1)) * &
              mix_sd(gen,1,2) * g_h0_sdLsdR(shiggs,gen) + &
              conjg(mix_sd(gen,1,2)) * mix_sd(gen,1,1) * &
              g_h0_sdRsdL(shiggs,gen))

         g_h0_sd2sd2(shiggs,gen) = (conjg(mix_sd(gen,2,1)) * &
              mix_sd(gen,2,1) * g_h0_sdLsdL(shiggs,gen) + &
              conjg(mix_sd(gen,2,2)) * mix_sd(gen,2,2) * &
              g_h0_sdRsdR(shiggs,gen) + conjg (mix_sd(gen,2,1)) * &
              mix_sd(gen,2,2) * g_h0_sdLsdR(shiggs,gen) + &
              conjg(mix_sd(gen,2,2)) * mix_sd(gen,2,1) * &
              g_h0_sdRsdL(shiggs,gen))

         g_h0_sd1sd2(shiggs,gen) = (conjg(mix_sd(gen,1,1)) * &
              mix_sd(gen,2,1) * g_h0_sdLsdL(shiggs,gen) + &
              conjg(mix_sd(gen,1,2)) * mix_sd(gen,2,2) * &
              g_h0_sdRsdR(shiggs,gen) + conjg(mix_sd(gen,1,1)) *&
              mix_sd(gen,2,2) * g_h0_sdLsdR(shiggs,gen) + &
              conjg(mix_sd(gen,1,2)) * mix_sd(gen,2,1) * &
              g_h0_sdRsdL(shiggs,gen))
         
         g_h0_sd2sd1(shiggs,gen) = conjg(g_h0_sd1sd2(shiggs,gen)) 
!!!SLep                                                                                                                         
        g_h0_sl1sl1(shiggs,gen) = (conjg (mix_sl(gen,1,1)) * &
             mix_sl(gen,1,1) * g_h0_slLslL(shiggs,gen) + &
             conjg(mix_sl(gen,1,2)) * mix_sl(gen,1,2) * &
             g_h0_slRslR(shiggs,gen) + conjg (mix_sl(gen,1,1)) * &
             mix_sl(gen,1,2) * g_h0_slLslR(shiggs,gen) + &
             conjg(mix_sl(gen,1,2)) * mix_sl(gen,1,1) * &
             g_h0_slRslL(shiggs,gen))
	
        g_h0_sl2sl2(shiggs,gen) = (conjg (mix_sl(gen,2,1)) * &
             mix_sl(gen,2,1) * g_h0_slLslL(shiggs,gen) + &
             conjg(mix_sl(gen,2,2)) * mix_sl(gen,2,2) * &
             g_h0_slRslR(shiggs,gen) + conjg (mix_sl(gen,2,1)) * &
             mix_sl(gen,2,2) * g_h0_slLslR(shiggs,gen) + &
             conjg(mix_sl(gen,2,2)) * mix_sl(gen,2,1) * &
             g_h0_slRslL(shiggs,gen))
	
	g_h0_sl1sl2(shiggs,gen) = (conjg (mix_sl(gen,1,1)) * &
             mix_sl(gen,2,1) * g_h0_slLslL(shiggs,gen) + &
             conjg(mix_sl(gen,1,2)) * mix_sl(gen,2,2) * &
             g_h0_slRslR(shiggs,gen) + conjg (mix_sl(gen,1,1)) * &
             mix_sl(gen,2,2) * g_h0_slLslR(shiggs,gen) + &
             conjg(mix_sl(gen,1,2)) * mix_sl(gen,2,1) * &
             g_h0_slRslL(shiggs,gen))

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
             (conjg (mix_su(gen,1,2)) * mix_su(gen,1,1) * g_A0_suLsuR(phiggs, gen) + &
             conjg (mix_su(gen,1,1)) * mix_su(gen,1,2) * g_A0_suRsuL(phiggs, gen))
        
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
        g_hp_su1sd1(gen1,gen2) = ((conjg (mix_su(gen1,1,1)) * &
             mix_sd(gen2,1,1) * g_hp_suLsdL(gen1,gen2)) + &
             ( conjg(mix_su(gen1,1,2)) * mix_sd(gen2,1,2) * &
             g_hp_suRsdR(gen1,gen2)) + ( conjg (mix_su(gen1,1,1)) * &
             mix_sd(gen2,1,2) * g_hp_suLsdR(gen1,gen2)) + &
             ( conjg(mix_su(gen1,1,2)) * mix_sd(gen2,1,1) *&
             g_hp_suRsdL(gen1,gen2)))
        g_hp_su1sd2(gen1,gen2) = ((conjg (mix_su(gen1,1,1)) * &
             mix_sd(gen2,2,1) * g_hp_suLsdL(gen1,gen2)) + &
             ( conjg(mix_su(gen1,1,2)) * mix_sd(gen2,2,2) *&
             g_hp_suRsdR(gen1,gen2)) + ( conjg (mix_su(gen1,1,1)) * &
             mix_sd(gen2,2,2) * g_hp_suLsdR(gen1,gen2)) + &
             ( conjg(mix_su(gen1,1,2)) * mix_sd(gen2,2,1) * &
             g_hp_suRsdL(gen1,gen2)))
        g_hp_su2sd1(gen1,gen2) = ((conjg (mix_su(gen1,2,1)) * &
             mix_sd(gen2,1,1) * g_hp_suLsdL(gen1,gen2)) + &
             ( conjg(mix_su(gen1,2,2)) * mix_sd(gen2,1,2) * &
             g_hp_suRsdR(gen1,gen2)) + ( conjg (mix_su(gen1,2,1)) * &
             mix_sd(gen2,1,2) * g_hp_suLsdR(gen1,gen2)) + &
             ( conjg(mix_su(gen1,2,2)) * mix_sd(gen2,1,1) * &
             g_hp_suRsdL(gen1,gen2)))
        g_hp_su2sd2(gen1,gen2) = ((conjg (mix_su(gen1,2,1)) * &
             mix_sd(gen2,2,1) * g_hp_suLsdL(gen1,gen2)) + &
             ( conjg(mix_su(gen1,2,2)) * mix_sd(gen2,2,2) * &
             g_hp_suRsdR(gen1,gen2)) + ( conjg (mix_su(gen1,2,1)) * &
             mix_sd(gen2,2,2) * g_hp_suLsdR(gen1,gen2)) + &
             ( conjg(mix_su(gen1,2,2)) * mix_sd(gen2,2,1) * &
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
    do shiggs2 = 1,dimh0
       do shiggs3 = 1,dimh0
          !!! checked by fb with FeynRules 8/31/09
          g_h0h0h0(shiggs1,shiggs2,shiggs3) = (&
                -(3.0_default / 4.0_default) * (gz**2) *&
                ( cosbe * vev * mix_h0(shiggs1,1) * mix_h0(shiggs2,1) * &
                mix_h0(shiggs3,1) + sinbe * vev * mix_h0(shiggs1,2) * &
                mix_h0(shiggs2,2) * mix_h0(shiggs3,2))  + &
                 (gz**2 / (4.0_default ) - ( lambda**2)) *&
               (  cosbe * vev * &
               (  mix_h0(shiggs1,1) * mix_h0(shiggs2,2) * mix_h0(shiggs3,2)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,1) * mix_h0(shiggs3,2)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,2) * mix_h0(shiggs3,1) )&
                + sinbe * vev * &
               (  mix_h0(shiggs1,1) * mix_h0(shiggs2,1) * mix_h0(shiggs3,2)&
               +  mix_h0(shiggs1,1) * mix_h0(shiggs2,2) * mix_h0(shiggs3,1)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,1) * mix_h0(shiggs3,1) ) )&
                + vev * ( k * lambda * sinbe - lambda**2 * cosbe) *&
               (  mix_h0(shiggs1,1) * mix_h0(shiggs2,3) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,1) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,3) * mix_h0(shiggs3,1)  )&
                + vev * ( k * lambda * cosbe - lambda**2 * sinbe) *&
               (  mix_h0(shiggs1,2) * mix_h0(shiggs2,3) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,2) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,3) * mix_h0(shiggs3,2)  )&
                - (lambda * mu * sqrt(2.0_default) )* &
               (  mix_h0(shiggs1,1) * mix_h0(shiggs2,1) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,1) * mix_h0(shiggs2,3) * mix_h0(shiggs3,1)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,1) * mix_h0(shiggs3,1)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,2) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,3) * mix_h0(shiggs3,2)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,2) * mix_h0(shiggs3,2)  )&
               + ( lambda * A_lambda / sqrt (2.0_default) &
                + sqrt(2.0_default) * mu * k ) * &
               (  mix_h0(shiggs1,1) * mix_h0(shiggs2,2) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,1) * mix_h0(shiggs2,3) * mix_h0(shiggs3,2)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,1) * mix_h0(shiggs3,2)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,1) * mix_h0(shiggs3,3)&
               +  mix_h0(shiggs1,2) * mix_h0(shiggs2,3) * mix_h0(shiggs3,1)&
               +  mix_h0(shiggs1,3) * mix_h0(shiggs2,2) * mix_h0(shiggs3,1)  )&
                - ( k * A_k  * sqrt (2.0_default)  + &
                 6.0_default * sqrt (2.0_default) * k**2 * mu / lambda ) * &
                 mix_h0(shiggs1,3) * mix_h0(shiggs2,3) * mix_h0(shiggs3,3)  )
       end do
    end do
 end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!h0A0A0
!!!!!!!!!!!!!!!!!!!!!!!!
  do shiggs1 = 1,dimh0
     do phiggs1 = 1,dimA0
        do phiggs2 = 1,dimA0
           g_h0A0A0(shiggs1,phiggs1,phiggs2) = &
                (-1.0_default / 4.0_default) * gz**2 * vev  *&
                 ( cosbe * mix_h0(shiggs1,1) * mix_A0(phiggs1,1) * mix_A0(phiggs2,1)&
                 + sinbe * mix_h0(shiggs1,2) * mix_A0(phiggs1,2) * mix_A0(phiggs2,2)  )&
                 !!!
                 + (gz**2 / 4.0_default - lambda**2 ) * vev * &
                 ( cosbe * mix_h0(shiggs1,1) * mix_A0(phiggs1,2) * mix_A0(phiggs2,2)&
                 + sinbe * mix_h0(shiggs1,2) * mix_A0(phiggs1,1) * mix_A0(phiggs2,1)  )&
                !!!
                 - (  k * lambda * cosbe + lambda**2 * sinbe) * vev * &
                  mix_h0(shiggs1,2) * mix_A0(phiggs1,3) * mix_A0(phiggs2,3)&
                 - (  k * lambda * sinbe + lambda**2 * cosbe) * vev * &
                  mix_h0(shiggs1,1) * mix_A0(phiggs1,3) * mix_A0(phiggs2,3)&
                !!!
                 - lambda * mu * sqrt (2.0_default) *  mix_h0(shiggs1,3) * &
                          ( mix_A0(phiggs1,1) * mix_A0(phiggs2,1)  &
                                      + mix_A0(phiggs1,2) * mix_A0(phiggs2,2)  ) &
                !!!
                 - sqrt(2.0_default)* ( 2.0_default * k**2 * mu / lambda - k * A_k )&
                 * mix_h0(shiggs1,3) * mix_A0(phiggs1,3) * mix_A0(phiggs2,3)  &
                 !!!
                 +  lambda * k  * vev *  mix_h0(shiggs1,3) * &
                  ( cosbe * ( mix_A0(phiggs1,2) * mix_A0(phiggs2,3) &
                            + mix_A0(phiggs1,3) * mix_A0(phiggs2,2) ) &
                  + sinbe * ( mix_A0(phiggs1,1) * mix_A0(phiggs2,3) &
                            + mix_A0(phiggs1,3) * mix_A0(phiggs2,1) ) ) &
                !!!
                 + ( k * mu * sqrt(2.0_default)&
                  - lambda * A_lambda / sqrt(2.0_default)  )*&
                  (mix_h0(shiggs1,1) * (  mix_A0(phiggs1,2) * mix_A0(phiggs2,3) &
                                        + mix_A0(phiggs1,3) * mix_A0(phiggs2,2) ) &
                 + mix_h0(shiggs1,2) * (  mix_A0(phiggs1,1) * mix_A0(phiggs2,3) &
                                        + mix_A0(phiggs1,3) * mix_A0(phiggs2,1) ) ) &
                !!!
                 - (k * mu * sqrt(2.0_default) &
                 + lambda * A_lambda / sqrt(2.0_default)  )*&
                   mix_h0(shiggs1,3) * (  mix_A0(phiggs1,1) * mix_A0(phiggs2,2) &
                 + mix_A0(phiggs1,2) * mix_A0(phiggs2,1) )
        end do
     end do
  end do
!!!!!!!!!!!!!!!!!!!!!!!!
!!!h0h+h-
!!!!!!!!!!!!!!!!!!!!!!!!
  do shiggs1 = 1,dimh0
     g_h0hphm(shiggs1) = &
           - g * mass(24) * (mix_h0(shiggs1,1) * cosbe + mix_h0(shiggs1,2) * sinbe ) &
           -(g * mass(23) / (2.0_default * costhw) ) * (mix_h0(shiggs1,2) * sinbe &
           - mix_h0(shiggs1,1) * cosbe) * cos2be &
           + ( lambda**2 / 2.0_default * vev ) * &
             (mix_h0(shiggs1,2) * cosbe + mix_h0(shiggs1,1) * sinbe) * sin2be  &
           - ( lambda / sqrt(2.0_default) ) * mix_h0(shiggs1,3) * &
          ( (2.0_default * k * mu / lambda  + A_lambda ) * sin2be + 2.0 * mu )
  end do
!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Scalar Higgs to Neutralinos
!!!!!!!!!!!!!!!!!!!!!!!!!!
  do neu1 = 1,dimNeu
     do neu2 = 1,dimNeu
        do shiggs1 = 1,dimh0
	   !! corrected convention w.r.t. Franke/Fraas
           g_neuneuh0_1(2,neu1,neu2,shiggs1) = &
                 (1.0_default / sqrt (2.0_default))* ( &
                - g / sqrt (2.0_default) * mix_neu(neu1,2) * &
                 ( mix_neu(neu2,3) * mix_h0(shiggs1,1)  &
                 - mix_neu(neu2,4) * mix_h0(shiggs1,2)) &
                + g * ( sinthw/costhw ) / sqrt (2.0_default) * mix_neu(neu1,1) * &
                 ( mix_neu(neu2,3) * mix_h0(shiggs1,1)  &
                 - mix_neu(neu2,4) * mix_h0(shiggs1,2)) &
                + lambda * mix_h0(shiggs1,3) * mix_neu(neu1,3) * mix_neu(neu2,4) &
		+ lambda * mix_h0(shiggs1,1) * mix_neu(neu1,4) * mix_neu(neu2,5) &
                + lambda * mix_h0(shiggs1,2) * mix_neu(neu1,3) * mix_neu(neu2,5) &
                - k * mix_h0(shiggs1,3) * mix_neu(neu1,5) * mix_neu(neu2,5) &
		)
        end do
     end do
  end do

  do neu1 = 1,dimNeu
     do neu2 = 1,dimNeu
        do shiggs1 = 1,dimh0
           g_neuneuh0(2,neu1,neu2,shiggs1) = &
                0.5_default * ( g_neuneuh0_1(2,neu1,neu2,shiggs1) &
                                 + g_neuneuh0_1(2,neu2,neu1,shiggs1))
	end do
     end do
  end do

  g_neuneuh0(1,:,:,:) = conjg (g_neuneuh0(2,:,:,:) )

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Axial Higgs to Neutralinos
!!!!!!!!!!!!!!!!!!!!!!!!!!
  do neu1 = 1,dimNeu
     do neu2 = 1,dimNeu
        do phiggs1 = 1,dimA0
           !! corrected convention w.r.t. Franke/Fraas
           g_neuneuA0_1(1,neu1,neu2,phiggs1) = &
                 (imago / sqrt (2.0_default))* ( &
                - g / sqrt (2.0_default) * mix_neu(neu1,2) * &
                 ( mix_neu(neu2,3) * mix_A0(phiggs1,1)  &
                 - mix_neu(neu2,4) * mix_A0(phiggs1,2)) &
                + g * ( sinthw/costhw ) / sqrt (2.0_default) * mix_neu(neu1,1) * &
                 ( mix_neu(neu2,3) * mix_A0(phiggs1,1)  &
                 - mix_neu(neu2,4) * mix_A0(phiggs1,2)) &
                - lambda * mix_A0(phiggs1,3) * mix_neu(neu1,3) * mix_neu(neu2,4) &
                - lambda * mix_A0(phiggs1,1) * mix_neu(neu1,4) * mix_neu(neu2,5) &
                - lambda * mix_A0(phiggs1,2) * mix_neu(neu1,3) * mix_neu(neu2,5) &
                + k * mix_A0(phiggs1,3) * mix_neu(neu1,5) * mix_neu(neu2,5) &
                )
        end do
     end do
  end do

  do neu1 = 1,dimNeu
     do neu2 = 1,dimNeu
        do phiggs1 = 1,dimA0
           g_neuneuA0(1,neu1,neu2,phiggs1) = &
                0.5_default * conjg( g_neuneuA0_1(1,neu1,neu2,phiggs1) &
                                      + g_neuneuA0_1(1,neu2,neu1,phiggs1))
        end do
     end do
  end do

  g_neuneuA0(2,:,:,:) = conjg (g_neuneuA0(1,:,:,:) )

!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!Scalar Higgs to Charginos
!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ch1 = 1,2
     do ch2 = 1,2
        do shiggs1 = 1,dimh0
           g_chchh0(1,ch1,ch2,shiggs1) =  - gcc * &
                (  mix_h0(shiggs1,1) * mix_charU(ch1,2) * mix_charV(ch2,1) &
                + mix_h0(shiggs1,2) * mix_charU(ch1,1) * mix_charV(ch2,2) ) &
                - (lambda / ( sqrt (8.0_default ))) * mix_h0(shiggs1,3) * &
                mix_charU(ch1,2) * mix_charV(ch2,2)
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
           g_chchA0(1,ch1,ch2,phiggs1) = imago *( gcc * &
                (mix_A0(phiggs1,1) * mix_charU(ch1,2) * mix_charV(ch2,1) &
                + mix_A0(phiggs1,2) * mix_charU(ch1,1) * mix_charV(ch2,2))  &
                -(lambda / sqrt(8.0_default)) * &
                mix_A0(phiggs1,3) * mix_charU(ch1,2) * mix_charV(ch2,2) )
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
     conjg(mix_neu(neu1,2)) ) * conjg(mix_charV(ch1,2) )) + &
     lambda * sinbe * conjg(mix_neu(neu1,5)) * conjg(mix_charV(ch1,2))) 

   g_neuhmchar(2,neu1,ch1) = (1.0_default / 2.0_default) * ( &
     g * sinbe * (( mix_neu(neu1,3) ) * mix_charU(ch1,1) - &
     ( 1.0_default / sqrt(2.0_default) ) * (( sinthw / costhw ) * &
     mix_neu(neu1,1) + &
     mix_neu(neu1,2) ) * mix_charU(ch1,2) ) + &
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
  igs = cmplx(0.0_default,1.0_default,kind=default) * gs
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
end module parameters_nmssm_ckm

