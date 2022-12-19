!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

module parameters_thdm_ckm
  use kinds
  use constants
  implicit none
  private

  real(default), dimension(37), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs

  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, &
       ghbb, ghtt, ghcc, ghtautau, ghmm, &
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn, &
       gh1uu, gh2uu, gh3uu, gh1uc, gh2uc, gh3uc, gh1ut, gh2ut, gh3ut, gh1cu, &
       gh2cu, gh3cu, gh1cc, gh2cc, gh3cc, gh1ct, gh2ct, gh3ct, gh1tu, gh2tu, &
       gh3tu, gh1tc, gh2tc, gh3tc, gh1tt, gh2tt, gh3tt, &
       gh1dd, gh2dd, gh3dd, gh1ds, gh2ds, gh3ds, gh1db, gh2db, gh3db, gh1sd, &
       gh2sd, gh3sd, gh1ss, gh2ss, gh3ss, gh1sb, gh2sb, gh3sb, gh1bd, gh2bd, &
       gh3bd, gh1bs, gh2bs, gh3bs, gh1bb, gh2bb, gh3bb,&
       ghud, ghus, ghub, ghcd, ghcs, ghcb, ghtd, ghts, ghtb,&
       ghdu, ghdc, ghdt, ghsu, ghsc, ghst, ghbu, ghbc, ghbt, &
       gh1e1e1, gh2e1e1, gh3e1e1, gh1e1e2, gh2e1e2, gh3e1e2, &
       gh1e1e3, gh2e1e3, gh3e1e3, gh1e2e1, gh2e2e1, gh3e2e1, &
       gh1e2e2, gh2e2e2, gh3e2e2, gh1e2e3, gh2e2e3, gh3e2e3, &
       gh1e3e1, gh2e3e1, gh3e3e1, gh1e3e2, gh2e3e2, gh3e3e2, &
       gh1e3e3, gh2e3e3, gh3e3e3
  real(default), public :: cotanthw2, sinthw2, &
       R11, R12, R13, R21, R22, R23, R31, R32, R33, &
       lamb1, lamb2, lamb3, lamb4, lamb5R, lamb5I, lamb6R, &
       lamb6I, lamb7R, lamb7I
  complex(default), public :: gAHpHm, gZHpHm, gZh1h2, gZh1h3, gZh2h3, &
       gWpHmh1, gWpHmh2, gWpHmh3, gWmHph1, gWmHph2, gWmHph3, gh1ZZ, &
       gh2ZZ, gh3ZZ, gh1WpWm, gh2WpWm, gh3WpWm, ghhWpWm, ghhZZ, gHpHmAA, &
       gHpHmZZ, gHpHmAZ, gHpHmWpWm, gh1HpAWm, gh2HpAWm, gh3HpAWm, &
       gh1HpZWm, gh2HpZWm, gh3HpZWm, gh1HpAWmC, gh2HpAWmC, gh3HpAWmC, &
       gh1HpZWmC, gh2HpZWmC, gh3HpZWmC, gh1HpHm, gh2HpHm, gh3HpHm, &
       gh111, gh112, gh113, gh221, gh222, gh223, gh331, gh332, gh333, &
       gh123, gHpHmHpHm, gHpHm11, gHpHm12, gHpHm13, gHpHm22, gHpHm23, &
       gHpHm33, gh1111, gh1112, gh1113, gh1122, gh1123, gh1133, gh1222, &
       gh1223, gh1233, gh1333, gh2222, gh2223, gh2233, gh2333, gh3333, &
       lamb5,lamb6,lamb7, &
       Yu11, Yu12, Yu13, Yu21, Yu22, Yu23, Yu31, Yu32, Yu33, &
       Yd11, Yd12, Yd13, Yd21, Yd22, Yd23, Yd31, Yd32, Yd33, &
       Yl11, Yl12, Yl13, Yl21, Yl22, Yl23, Yl31, Yl32, Yl33, &
       V11, V12, V13, V21, V22, V23, V31, V32, V33, &
       ghe1n1, ghe1n2, ghe1n3, ghe2n1, ghe2n2, ghe2n3, ghe3n1, &
       ghe3n2, ghe3n3, ghn1e1, ghn1e2, ghn1e3, ghn2e1, ghn2e2, &
       ghn2e3, ghn3e1, ghn3e2, ghn3e3

  public :: import_from_whizard, model_update_alpha_s, init_parameters

contains

  subroutine init_parameters
    real(default), dimension(86) :: vars
    vars(1)  = 1.16639E-5   ! Fermi constant
    vars(2)  = 91.1882      ! Z-boson mass
    vars(3)  = 80.419       ! W-boson mass
    vars(4)  = 0.1178       ! Strong coupling constant (Z point)
    vars(5)  = 0.000510997  ! electron mass
    vars(6)  = 0.105658389  ! muon mass
    vars(7)  = 1.77705      ! tau-lepton mass
    vars(8)  = 0.095        ! s-quark mass
    vars(9)  = 1.2          ! c-quark mass
    vars(10) = 4.2          ! b-quark mass
    vars(11) = 173.1        ! t-quark mass
    vars(12) = 1.523        ! t-quark width
    vars(13) = 2.443        ! Z-boson width
    vars(14) = 2.049        ! W-boson width
    vars(15) = 0.000        ! anomaly Higgs couplings K factors
    vars(16) = 0.000        ! anomaly Higgs couplings K factors
    vars(17) = 0.000        ! anomaly Higgs couplings K factors
    vars(18) = 1.000        ! R_xi parameter for Z-boson
    vars(19) = 1.000        ! R_xi parameter for W-boson
    vars(20) = 125          ! Higgs h mass
    vars(21) = 0.004143     ! Higgs h width
    vars(22) = 500          ! Higgs H mass
    vars(23) = 1.000        ! Higgs H width
    vars(24) = 500          ! Higgs A mass
    vars(25) = 1.000        ! Higgs A width
    vars(26) = 500          ! Higgs charged mass
    vars(27) = 1.000        ! Higgs charged width
    vars(28) = 1.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(29) = 0.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(30) = 0.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(31) = 0.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(32) = 1.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(33) = 0.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(34) = 0.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(35) = 0.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(36) = 1.000        ! Mix. matrix: Si = Rji hj; hj= (h,H,A)
    vars(37) = 1.0          ! Higgs potential parameter
    vars(38) = 1.0          ! Higgs potential parameter
    vars(39) = 1.0          ! Higgs potential parameter
    vars(40) = 1.0          ! Higgs potential parameter
    vars(41) = 0.0          ! Higgs potential parameter
    vars(42) = 1.0          ! Higgs potential parameter
    vars(43) = 0.0          ! Higgs potential parameter
    vars(44) = 1.0          ! Higgs potential parameter
    vars(45) = 1.0          ! Higgs potential parameter
    vars(46) = 0.0          ! Higgs potential parameter
    vars(47) = 0.000   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(48) = 0.100   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(49) = 0.200   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(50) = 0.300   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(51) = 0.400   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(52) = 0.500   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(53) = 0.600   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(54) = 0.700   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(55) = 0.800   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(56) = 0.900   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(57) = 0.800   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(58) = 0.700   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(59) = 0.600   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(60) = 0.500   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(61) = 0.400   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(62) = 0.300   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(63) = 0.200   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(64) = 0.100   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(65) = 0.000   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(66) = 0.100   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(67) = 0.200   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(68) = 0.300   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(69) = 0.400   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(70) = 0.500   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(71) = 0.600   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(72) = 0.700   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(73) = 0.800   ! Yukawa matrix to the 0 vev Higgs doublet
    vars(74) = 1.000   ! CKM Matrix
    vars(75) = 0.000   ! CKM Matrix
    vars(76) = 0.000   ! CKM Matrix
    vars(77) = 0.000   ! CKM Matrix
    vars(78) = 1.000   ! CKM Matrix
    vars(79) = 0.000   ! CKM Matrix
    vars(80) = 0.000   ! CKM Matrix
    vars(81) = 0.000   ! CKM Matrix
    vars(82) = 1.000   ! CKM Matrix
    vars(83) = 1 / sqrt (sqrt (2.) * vars(1))      ! v (Higgs vev)
    vars(84) = vars(3) / vars(2)                   ! cos(theta-W)
    vars(85) = sqrt (1-vars(84)**2)                ! sin(theta-W)
    vars(86) = 2 * vars(85) * vars(3) / vars(83)   ! em-coupling (GF scheme)
    call import_from_whizard (vars)
  end subroutine init_parameters

  subroutine import_from_whizard (par_array)
    real(default), dimension(86), intent(in) :: par_array
    type :: parameter_set
       real(default) :: gf
       real(default) :: mZ
       real(default) :: mW
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
       real(default) :: khgaz
       real(default) :: khgaga
       real(default) :: khgg
       real(default) :: xi0
       real(default) :: xipm
       real(default) :: mh1
       real(default) :: wh1
       real(default) :: mh2
       real(default) :: wh2
       real(default) :: mh3
       real(default) :: wh3
       real(default) :: mHp
       real(default) :: wHp
       real(default) :: R_11
       real(default) :: R_12
       real(default) :: R_13
       real(default) :: R_21
       real(default) :: R_22
       real(default) :: R_23
       real(default) :: R_31
       real(default) :: R_32
       real(default) :: R_33
       real(default) :: lamb_1
       real(default) :: lamb_2
       real(default) :: lamb_3
       real(default) :: lamb_4
       real(default) :: lamb_5R
       real(default) :: lamb_5I
       real(default) :: lamb_6R
       real(default) :: lamb_6I
       real(default) :: lamb_7R
       real(default) :: lamb_7I
       real(default) :: Yu_11
       real(default) :: Yu_12
       real(default) :: Yu_13
       real(default) :: Yu_21
       real(default) :: Yu_22
       real(default) :: Yu_23
       real(default) :: Yu_31
       real(default) :: Yu_32
       real(default) :: Yu_33
       real(default) :: Yd_11
       real(default) :: Yd_12
       real(default) :: Yd_13
       real(default) :: Yd_21
       real(default) :: Yd_22
       real(default) :: Yd_23
       real(default) :: Yd_31
       real(default) :: Yd_32
       real(default) :: Yd_33
       real(default) :: Yl_11
       real(default) :: Yl_12
       real(default) :: Yl_13
       real(default) :: Yl_21
       real(default) :: Yl_22
       real(default) :: Yl_23
       real(default) :: Yl_31
       real(default) :: Yl_32
       real(default) :: Yl_33
       real(default) :: V_11
       real(default) :: V_12
       real(default) :: V_13
       real(default) :: V_21
       real(default) :: V_22
       real(default) :: V_23
       real(default) :: V_31
       real(default) :: V_32
       real(default) :: V_33
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
    par%gf      = par_array(1)
    par%mZ      = par_array(2)
    par%mW      = par_array(3)
    par%alphas  = par_array(4)
    par%me      = par_array(5)
    par%mmu     = par_array(6)
    par%mtau    = par_array(7)
    par%ms      = par_array(8)
    par%mc      = par_array(9)
    par%mb      = par_array(10)
    par%mtop    = par_array(11)
    par%wtop    = par_array(12)
    par%wZ      = par_array(13)
    par%wW      = par_array(14)
    par%khgaz   = par_array(15)
    par%khgaga  = par_array(16)
    par%khgg    = par_array(17)
    par%xi0     = par_array(18)
    par%xipm    = par_array(19)
    par%mh1     = par_array(20)
    par%wh1     = par_array(21)
    par%mh2     = par_array(22)
    par%wh2     = par_array(23)
    par%mh3     = par_array(24)
    par%wh3     = par_array(25)
    par%mHp     = par_array(26)
    par%wHp     = par_array(27)
    par%R_11    = par_array(28)
    par%R_12    = par_array(29)
    par%R_13    = par_array(30)
    par%R_21    = par_array(31)
    par%R_22    = par_array(32)
    par%R_23    = par_array(33)
    par%R_31    = par_array(34)
    par%R_32    = par_array(35)
    par%R_33    = par_array(36)
    par%lamb_1  = par_array(37)
    par%lamb_2  = par_array(38)
    par%lamb_3  = par_array(39)
    par%lamb_4  = par_array(40)
    par%lamb_5R = par_array(41)
    par%lamb_5I = par_array(42)
    par%lamb_6R = par_array(43)
    par%lamb_6I = par_array(44)
    par%lamb_7R = par_array(45)
    par%lamb_7I = par_array(46)
    par%Yu_11   = par_array(47)
    par%Yu_12   = par_array(48)
    par%Yu_13   = par_array(49)
    par%Yu_21   = par_array(50)
    par%Yu_22   = par_array(51)
    par%Yu_23   = par_array(52)
    par%Yu_31   = par_array(53)
    par%Yu_32   = par_array(54)
    par%Yu_33   = par_array(55)
    par%Yd_11   = par_array(56)
    par%Yd_12   = par_array(57)
    par%Yd_13   = par_array(58)
    par%Yd_21   = par_array(59)
    par%Yd_22   = par_array(60)
    par%Yd_23   = par_array(61)
    par%Yd_31   = par_array(62)
    par%Yd_32   = par_array(63)
    par%Yd_33   = par_array(64)
    par%Yl_11   = par_array(65)
    par%Yl_12   = par_array(66)
    par%Yl_13   = par_array(67)
    par%Yl_21   = par_array(68)
    par%Yl_22   = par_array(69)
    par%Yl_23   = par_array(70)
    par%Yl_31   = par_array(71)
    par%Yl_32   = par_array(72)
    par%Yl_33   = par_array(73)
    par%V_11    = par_array(74)
    par%V_12    = par_array(75)
    par%V_13    = par_array(76)
    par%V_21    = par_array(77)
    par%V_22    = par_array(78)
    par%V_23    = par_array(79)
    par%V_31    = par_array(80)
    par%V_32    = par_array(81)
    par%V_33    = par_array(82)
    par%v       = par_array(83)
    par%cw      = par_array(84)
    par%sw      = par_array(85)
    par%ee      = par_array(86)
    mass(1:37) = 0
    width(1:37) = 0
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
    mass(25) = par%mh1
    width(25) = par%wh1
    mass(35) = par%mh2
    width(35) = par%wh2
    mass(36) = par%mh3
    width(36) = par%wh3
    mass(27) = par%mHp
    width(27) = par%wHp
    mass(26) = par%xi0 * mass(23)
    width(26) = 0
    mass(37) = par%xipm * mass(24)
    width(37) = 0
    R11 = par%R_11
    R12 = par%R_12
    R13 = par%R_13
    R21 = par%R_21
    R22 = par%R_22
    R23 = par%R_23
    R31 = par%R_31
    R32 = par%R_32
    R33 = par%R_33
    Yu11 = par%Yu_11
    Yu12 = par%Yu_12
    Yu13 = par%Yu_13
    Yu21 = par%Yu_21
    Yu22 = par%Yu_22
    Yu23 = par%Yu_23
    Yu31 = par%Yu_31
    Yu32 = par%Yu_32
    Yu33 = par%Yu_33
    Yd11 = par%Yd_11
    Yd12 = par%Yd_12
    Yd13 = par%Yd_13
    Yd21 = par%Yd_21
    Yd22 = par%Yd_22
    Yd23 = par%Yd_23
    Yd31 = par%Yd_31
    Yd32 = par%Yd_32
    Yd33 = par%Yd_33
    Yl11 = par%Yl_11
    Yl12 = par%Yl_12
    Yl13 = par%Yl_13
    Yl21 = par%Yl_21
    Yl22 = par%Yl_22
    Yl23 = par%Yl_23
    Yl31 = par%Yl_31
    Yl32 = par%Yl_32
    Yl33 = par%Yl_33
    lamb1 = par%lamb_1
    lamb2 = par%lamb_2
    lamb3 = par%lamb_3
    lamb4 = par%lamb_4
    lamb5R = par%lamb_5R
    lamb5I = par%lamb_5I
    lamb6R = par%lamb_6R
    lamb6I = par%lamb_6I
    lamb7R = par%lamb_7R
    lamb7I = par%lamb_7I
    lamb5 = par%lamb_5R + imago*par%lamb_5I
    lamb6 = par%lamb_6R + imago*par%lamb_6I
    lamb7 = par%lamb_7R + imago*par%lamb_7I
    V11 = par%V_11
    V12 = par%V_12
    V13 = par%V_13
    V21 = par%V_21
    V22 = par%V_22
    V23 = par%V_23
    V31 = par%V_31
    V32 = par%V_32
    V33 = par%V_33
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
    iqw = (0,1)*qw
    gzww = g * costhw
    igzww = (0,1)*gzww
    gwww = g
    igwww = (0,1)*gwww
    gw4 = gwww**2
    gzzww = gzww**2
    gazww = gzww * qw
    gaaww = qw**2
    ghtt = - mass(6) / vev
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    ghmm = - mass(13) / vev
    !!!
    cotanthw2 = (1.0-tanthw**2)/(2.0*tanthw)
    sinthw2 = 2.0*sinthw*costhw
    !Gauge boson with 2 Higgs coupligns
    gAHpHm = e
    gZHpHm = e * cotanthw2
    gZh1h2 = - (0,1)*(e/sinthw2)*(R13*R22-R12*R23)/2.0
    gZh1h3 = -(0,1)*(e/sinthw2)*(R13*R32-R12*R33)/2.0
    gZh2h3 = -(0,1)*(e/sinthw2)*(R23*R32-R22*R33)/2.0
    gWpHmh1 = -(g/2.0)*(R12+(0,1)*R13)
    gWpHmh2 = -(g/2.0)*(R22+(0,1)*R23)
    gWpHmh3 = -(g/2.0)*(R32+(0,1)*R33)
    gWmHph1 = (g/2.0)*(R12-(0,1)*R13)
    gWmHph2 = (g/2.0)*(R22-(0,1)*R23)
    gWmHph3 = (g/2.0)*(R32-(0,1)*R33)
    !Higgs + Gauge boson + Gauge boson couplings
    gh1ZZ = (mass(23)**2) * R11 / vev
    gh2ZZ = (mass(23)**2) * R21 / vev
    gh3ZZ = (mass(23)**2) * R31 / vev
    gh1WpWm = 2.0*(mass(24)**2) * R11 / vev
    gh2WpWm = 2.0*(mass(24)**2) * R21 / vev
    gh3WpWm = 2.0*(mass(24)**2) * R31 / vev
    !Higgs+Higgs+gauge+gauge
    ghhWpWm = (mass(24)/vev)**2
    ghhZZ = ((mass(23)/vev)**2) / 2.0
    gHpHmAA = e**2
    gHpHmZZ = (e*cotanthw2)**2
    gHpHmAZ = 2.0*cotanthw2*e**2
    gHpHmWpWm = (g**2)/2.0
    gh1HpAWm = e*g*(R12-(0,1)*R13) / 2.0
    gh2HpAWm = e*g*(R22-(0,1)*R23) / 2.0
    gh3HpAWm = e*g*(R32-(0,1)*R33) / 2.0
    gh1HpZWm = - tanthw*gh1HpAWm
    gh2HpZWm = - tanthw*gh2HpAWm
    gh3HpZWm = - tanthw*gh3HpAWm
    gh1HpAWmC = conjg (gh1HpAWm)
    gh2HpAWmC = conjg (gh2HpAWm)
    gh3HpAWmC = conjg (gh3HpAWm)
    gh1HpZWmC = conjg (gh1HpZWm)
    gh2HpZWmC = conjg (gh2HpZWm)
    gh3HpZWmC = conjg (gh3HpZWm)
    !!! Cubic Higgs
    gh1HpHm = vev * (lamb3*R11 - R13*aimag(lamb7) &
         + R12*real(lamb7))
    gh2HpHm = vev * (lamb3*R21 - R23*aimag(lamb7) &
         + R22*real(lamb7))
    gh3HpHm = vev * (lamb3*R31 - R33*aimag(lamb7) + R32*real(lamb7))
    gh111 = (vev*(2*lamb1*R11**3 + R11*(-4*R12*R13*aimag(lamb5) &
         + R13**2*(lamb3 + lamb4 - 2*real(lamb5)) + &
         R12**2*(lamb3 + lamb4 + 2*real(lamb5))) + &
         R11**2*(-3*R13*aimag(lamb6) + 3*R12*real(lamb6)) + &
         (R12**2 + R13**2)*(-(R13*aimag(lamb7)) + R12*real(lamb7))))/2.0
    gh112 = (vev*(3*R11**2*(2*lamb1*R21 - R23*aimag(lamb6) &
         + R22*real(lamb6)) + 2*R11*(R13*(-2*R22*aimag(lamb5) &
         - 3*R21*aimag(lamb6) + R23*(lamb3 + lamb4 - &
         2*real(lamb5))) + R12*(-2*R23*aimag(lamb5) + R22*(lamb3 + lamb4 + &
         2*real(lamb5)) + 3*R21*real(lamb6))) + R13**2*(-3*R23*aimag(lamb7) + &
         R21*(lamb3 + lamb4 - 2*real(lamb5)) + R22*real(lamb7)) + &
         R12**2*(-(R23*aimag(lamb7)) + R21*(lamb3 + lamb4 + 2*real(lamb5)) + &
         3*R22*real(lamb7)) - 2*R12*R13*(2*R21*aimag(lamb5) + &
         R22*aimag(lamb7) - R23*real(lamb7))))/2.0
    gh113 = (vev*(3*R11**2*(2*lamb1*R31 - R33*aimag(lamb6) + &
         R32*real(lamb6)) + 2*R11*(R13*(-2*R32*aimag(lamb5) - &
         3*R31*aimag(lamb6) + R33*(lamb3 + lamb4 - &
         2*real(lamb5))) + R12*(-2*R33*aimag(lamb5) + R32*(lamb3 + lamb4 + &
         2*real(lamb5)) + 3*R31*real(lamb6))) + R13**2*(-3*R33*aimag(lamb7) + &
         R31*(lamb3 + lamb4 - 2*real(lamb5)) + R32*real(lamb7)) + &
         R12**2*(-(R33*aimag(lamb7)) + R31*(lamb3 + lamb4 + 2*real(lamb5)) + &
         3*R32*real(lamb7)) - 2*R12*R13*(2*R31*aimag(lamb5) + &
         R32*aimag(lamb7) - R33*real(lamb7))))/2.0
    gh221 = (vev*(R11*(6*lamb1*R21**2 - 4*R22*R23*aimag(lamb5) + &
         R23**2*(lamb3 + lamb4 - 2*real(lamb5)) + R22**2*(lamb3 + &
         lamb4 + 2*real(lamb5)) + R21*(-6*R23*aimag(lamb6) + &
         6*R22*real(lamb6))) - R13*(3*R21**2*aimag(lamb6) + &
         R22**2*aimag(lamb7) + 3*R23**2*aimag(lamb7) - &
         2*R21*(-2*R22*aimag(lamb5) + R23*(lamb3 + lamb4 - 2*real(lamb5))) &
         - 2*R22*R23*real(lamb7)) + R12*(-2*R22*R23*aimag(lamb7) + &
         2*R21*(-2*R23*aimag(lamb5) + R22*(lamb3 + &
         lamb4 + 2*real(lamb5))) + 3*R21**2*real(lamb6) + &
         3*R22**2*real(lamb7) + R23**2*real(lamb7))))/2.0
    gh222 = (vev*(2*lamb1*R21**3 + R21*(-4*R22*R23*aimag(lamb5) &
         + R23**2*(lamb3 + lamb4 - 2*real(lamb5)) + R22**2*(lamb3 + &
         lamb4 + 2*real(lamb5))) + R21**2*(-3*R23*aimag(lamb6) + &
         3*R22*real(lamb6)) + (R22**2 + R23**2)*(-(R23*aimag(lamb7)) &
         + R22*real(lamb7))))/2.0
    gh223 = (vev*(3*R21**2*(2*lamb1*R31 - R33*aimag(lamb6) + &
         R32*real(lamb6)) + 2*R21*(R23*(-2*R32*aimag(lamb5) - &
         3*R31*aimag(lamb6) + R33*(lamb3 + lamb4 - &
         2*real(lamb5))) + R22*(-2*R33*aimag(lamb5) + R32*(lamb3 + lamb4 + &
         2*real(lamb5)) + 3*R31*real(lamb6))) + R23**2*(-3*R33*aimag(lamb7) + &
         R31*(lamb3 + lamb4 - 2*real(lamb5)) + R32*real(lamb7)) + &
         R22**2*(-(R33*aimag(lamb7)) + R31*(lamb3 + lamb4 + 2*real(lamb5)) + &
         3*R32*real(lamb7)) - 2*R22*R23*(2*R31*aimag(lamb5) + &
         R32*aimag(lamb7) - R33*real(lamb7))))/2.0
    gh331 = (vev*(R11*(6*lamb1*R31**2 - 4*R32*R33*aimag(lamb5) + &
         R33**2*(lamb3 + lamb4 - 2*real(lamb5)) + R32**2*(lamb3 + lamb4 &
         + 2*real(lamb5)) + R31*(-6*R33*aimag(lamb6) + 6*R32*real(lamb6))) &
         - R13*(3*R31**2*aimag(lamb6) + R32**2*aimag(lamb7) + &
         3*R33**2*aimag(lamb7) - 2*R31*(-2*R32*aimag(lamb5) + &
         R33*(lamb3 + lamb4 - 2*real(lamb5))) - 2*R32*R33*real(lamb7)) + &
         R12*(-2*R32*R33*aimag(lamb7) + 2*R31*(-2*R33*aimag(lamb5) + &
         R32*(lamb3 + lamb4 + 2*real(lamb5))) + 3*R31**2*real(lamb6) + &
         3*R32**2*real(lamb7) + R33**2*real(lamb7))))/2.0
    gh332 = (vev*(R21*(6*lamb1*R31**2 - 4*R32*R33*aimag(lamb5) + &
         R33**2*(lamb3 + lamb4 - 2*real(lamb5)) + R32**2*(lamb3 + lamb4 &
         + 2*real(lamb5)) + R31*(-6*R33*aimag(lamb6) + 6*R32*real(lamb6))) &
         - R23*(3*R31**2*aimag(lamb6) + R32**2*aimag(lamb7) + &
         3*R33**2*aimag(lamb7) - 2*R31*(-2*R32*aimag(lamb5) + &
         R33*(lamb3 + lamb4 - 2*real(lamb5))) - 2*R32*R33*real(lamb7)) + &
         R22*(-2*R32*R33*aimag(lamb7) + 2*R31*(-2*R33*aimag(lamb5) + &
         R32*(lamb3 + lamb4 + 2*real(lamb5))) + 3*R31**2*real(lamb6) + &
         3*R32**2*real(lamb7) + R33**2*real(lamb7))))/2.0
    gh333 = (vev*(2*lamb1*R31**3 + R31*(-4*R32*R33*aimag(lamb5) + &
         R33**2*(lamb3 + lamb4 - 2*real(lamb5)) + R32**2*(lamb3 + &
         lamb4 + 2*real(lamb5))) + R31**2*(-3*R33*aimag(lamb6) + &
         3*R32*real(lamb6)) + (R32**2 + R33**2)*(-(R33*aimag(lamb7)) + &
         R32*real(lamb7))))/2.0
    gh123 = vev*(R11*(R23*(-2*R32*aimag(lamb5) - 3*R31*aimag(lamb6) + &
         R33*(lamb3 + lamb4 - 2*real(lamb5))) + R22*(-2*R33*aimag(lamb5) &
         + R32*(lamb3 + lamb4 + 2*real(lamb5)) + 3*R31*real(lamb6)) + &
         3*R21*(2*lamb1*R31 - R33*aimag(lamb6) + R32*real(lamb6))) + &
         R13*(R21*(-2*R32*aimag(lamb5) - 3*R31*aimag(lamb6) + R33*(lamb3 &
         + lamb4 - 2*real(lamb5))) + R23*(-3*R33*aimag(lamb7) + &
         R31*(lamb3 + lamb4 - 2*real(lamb5)) + &
         R32*real(lamb7)) + R22*(-2*R31*aimag(lamb5) - R32*aimag(lamb7) + &
         R33*real(lamb7))) + R12*(R21*(-2*R33*aimag(lamb5) + R32*(lamb3 + &
         lamb4 + 2*real(lamb5)) + 3*R31*real(lamb6)) + &
         R22*(-(R33*aimag(lamb7)) + R31*(lamb3 + lamb4 + 2*real(lamb5)) &
         + 3*R32*real(lamb7)) + R23*(-2*R31*aimag(lamb5) - R32*aimag(lamb7) &
         + R33*real(lamb7))))
    !!! Quartic charged Higgs
    gHpHmHpHm = lamb2
    gHpHm11 = (lamb3*R11**2)/2. + lamb7R*R11*R12 + lamb2*R12**2 - &
         lamb7I*R11*R13 + lamb2*R13**2
    gHpHm12 = lamb3*R11*R21 + lamb7R*R12*R21 - lamb7I*R13*R21 + &
         lamb7R*R11*R22 + 2*lamb2*R12*R22 - lamb7I*R11*R23 + 2*lamb2*R13*R23
    gHpHm13 = lamb3*R11*R31 + lamb7R*R12*R31 - lamb7I*R13*R31 + &
         lamb7R*R11*R32 + 2*lamb2*R12*R32 - lamb7I*R11*R33 + 2*lamb2*R13*R33
    gHpHm22 = (lamb3*R21**2)/2. + lamb7R*R21*R22 + lamb2*R22**2 - &
         lamb7I*R21*R23 + lamb2*R23**2
    gHpHm23 = lamb3*R21*R31 + lamb7R*R22*R31 - lamb7I*R23*R31 + &
         lamb7R*R21*R32 + 2*lamb2*R22*R32 - lamb7I*R21*R33 + 2*lamb2*R23*R33
    gHpHm33 = (lamb3*R31**2)/2. + lamb7R*R31*R32 + lamb2*R32**2 - &
         lamb7I*R31*R33 + lamb2*R33**2
    !!! Quartic neutral higgs
    gh1111 = (lamb1*R11**4 + R11**3*(2*lamb6R*R12 - 2*lamb6I*R13) + &
         2*R11*(lamb7R*R12 - lamb7I*R13)*(R12**2 + R13**2) + lamb2*(R12**2 + &
         R13**2)**2 + R11**2*((lamb3 + lamb4 + 2*lamb5R)*R12**2 - &
         4*lamb5I*R12*R13 + (lamb3 + lamb4 - 2*lamb5R)*R13**2))/4.0
    gh1112 = (R11**3*(2*lamb1*R21 + lamb6R*R22 - lamb6I*R23) + &
         (R12**2 + R13**2)*(lamb7R*R12*R21 - lamb7I*R13*R21 + &
         2*lamb2*R12*R22 + 2*lamb2*R13*R23) + R11**2*(R12*(3*lamb6R*R21 + &
         (lamb3 + lamb4 + 2*lamb5R)*R22 - 2*lamb5I*R23) + &
         R13*(-3*lamb6I*R21 - 2*lamb5I*R22 + (lamb3 + lamb4 - &
         2*lamb5R)*R23)) + R11*(R13**2*((lamb3 + lamb4 - &
         2*lamb5R)*R21 + lamb7R*R22 - 3*lamb7I*R23) + R12**2*((lamb3 + &
         lamb4 + 2*lamb5R)*R21 + 3*lamb7R*R22 - lamb7I*R23) - &
         2*R12*R13*(2*lamb5I*R21 + lamb7I*R22 - lamb7R*R23)))/2.0
    gh1113 = (R11**3*(2*lamb1*R31 + lamb6R*R32 - lamb6I*R33) + (R12**2 + &
         R13**2)*(lamb7R*R12*R31 - lamb7I*R13*R31 + 2*lamb2*R12*R32 + &
         2*lamb2*R13*R33) + R11**2*(R12*(3*lamb6R*R31 + (lamb3 + lamb4 + &
         2*lamb5R)*R32 - 2*lamb5I*R33) + R13*(-3*lamb6I*R31 - 2*lamb5I*R32 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33)) + R11*(R13**2*((lamb3 + lamb4 - &
         2*lamb5R)*R31 + lamb7R*R32 - 3*lamb7I*R33) + R12**2*((lamb3 + &
         lamb4 + 2*lamb5R)*R31 + 3*lamb7R*R32 - lamb7I*R33) - &
         2*R12*R13*(2*lamb5I*R31 + lamb7I*R32 - lamb7R*R33)))/2.0
    gh1122 = (4*R12*R13*(-(lamb5I*R21**2) - lamb7I*R21*R22 + lamb7R*R21*R23 + &
         2*lamb2*R22*R23) + R12**2*((lamb3 + lamb4 + 2*lamb5R)*R21**2 + &
         6*lamb7R*R21*R22 + 6*lamb2*R22**2 - 2*lamb7I*R21*R23 + &
         2*lamb2*R23**2) + R13**2*((lamb3 + lamb4 - 2*lamb5R)*R21**2 + &
         2*lamb7R*R21*R22 + 2*lamb2*R22**2 - 6*lamb7I*R21*R23 + &
         6*lamb2*R23**2) + R11**2*(6*lamb1*R21**2 + (lamb3 + lamb4 + &
         2*lamb5R)*R22**2 - 4*lamb5I*R22*R23 + (lamb3 + lamb4 - &
         2*lamb5R)*R23**2 + R21*(6*lamb6R*R22 - 6*lamb6I*R23)) + &
         2*R11*(R12*(3*lamb6R*R21**2 + 3*lamb7R*R22**2 - 2*lamb7I*R22*R23 &
         + lamb7R*R23**2 + 2*R21*((lamb3 + lamb4 + 2*lamb5R)*R22 &
         - 2*lamb5I*R23)) - R13*(3*lamb6I*R21**2 + &
         lamb7I*R22**2 - 2*lamb7R*R22*R23 + 3*lamb7I*R23**2 - &
         2*R21*(-2*lamb5I*R22 + (lamb3 + lamb4 - 2*lamb5R)*R23))))/4.0
    gh1123 = (2*R12*R13*(-2*lamb5I*R21*R31 - lamb7I*R22*R31 + &
         lamb7R*R23*R31 - lamb7I*R21*R32 + 2*lamb2*R23*R32 + lamb7R*R21*R33 &
         + 2*lamb2*R22*R33) + R11**2*(R22*(3*lamb6R*R31 + (lamb3 + lamb4 + &
         2*lamb5R)*R32 - 2*lamb5I*R33) + R23*(-3*lamb6I*R31 - &
         2*lamb5I*R32 + (lamb3 + lamb4 - 2*lamb5R)*R33) + &
         3*R21*(2*lamb1*R31 + lamb6R*R32 - lamb6I*R33)) + &
         R13**2*(R22*(lamb7R*R31 + 2*lamb2*R32) + R23*(-3*lamb7I*R31 &
         + 6*lamb2*R33) + R21*((lamb3 + lamb4 - 2*lamb5R)*R31 + &
         lamb7R*R32 - 3*lamb7I*R33)) + R12**2*(3*R22*(lamb7R*R31 + &
         2*lamb2*R32) + R23*(-(lamb7I*R31) + 2*lamb2*R33) + R21*((lamb3 &
         + lamb4 + 2*lamb5R)*R31 + 3*lamb7R*R32 - lamb7I*R33)) + &
         2*R11*(R13*(R21*(-3*lamb6I*R31 - 2*lamb5I*R32 + (lamb3 + lamb4 - &
         2*lamb5R)*R33) + R23*((lamb3 + lamb4 - 2*lamb5R)*R31 + lamb7R*R32 - &
         3*lamb7I*R33) + R22*(-2*lamb5I*R31 - lamb7I*R32 + lamb7R*R33)) + &
         R12*(R21*(3*lamb6R*R31 + (lamb3 + lamb4 + 2*lamb5R)*R32 - &
         2*lamb5I*R33) + R22*((lamb3 + lamb4 + 2*lamb5R)*R31 + &
         3*lamb7R*R32 - lamb7I*R33) + R23*(-2*lamb5I*R31 - &
         lamb7I*R32 + lamb7R*R33))))/2.0
    gh1133 = (4*R12*R13*(-(lamb5I*R31**2) - lamb7I*R31*R32 + lamb7R*R31*R33 + &
         2*lamb2*R32*R33) + R12**2*((lamb3 + lamb4 + 2*lamb5R)*R31**2 + &
         6*lamb7R*R31*R32 + 6*lamb2*R32**2 - 2*lamb7I*R31*R33 + &
         2*lamb2*R33**2) + R13**2*((lamb3 + lamb4 - 2*lamb5R)*R31**2 + &
         2*lamb7R*R31*R32 + 2*lamb2*R32**2 - 6*lamb7I*R31*R33 + &
         6*lamb2*R33**2) + R11**2*(6*lamb1*R31**2 + (lamb3 + lamb4 + &
         2*lamb5R)*R32**2 - 4*lamb5I*R32*R33 + (lamb3 + lamb4 - &
         2*lamb5R)*R33**2 + R31*(6*lamb6R*R32 - 6*lamb6I*R33)) + &
         2*R11*(R12*(3*lamb6R*R31**2 + 3*lamb7R*R32**2 - 2*lamb7I*R32*R33 &
         + lamb7R*R33**2 + 2*R31*((lamb3 + lamb4 + 2*lamb5R)*R32 - &
         2*lamb5I*R33)) - R13*(3*lamb6I*R31**2 + lamb7I*R32**2 - &
         2*lamb7R*R32*R33 + 3*lamb7I*R33**2 - 2*R31*(-2*lamb5I*R32 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33))))/4.0
    gh1222 = (R13*(-(lamb6I*R21**3) + (lamb3 + lamb4 - &
         2*lamb5R)*R21**2*R23 - 3*lamb7I*R21*R23**2 + 2*lamb2*R23**3 + &
         R22**2*(-(lamb7I*R21) + 2*lamb2*R23) + 2*R21*R22*(-(lamb5I*R21) &
         + lamb7R*R23)) + R11*(2*lamb1*R21**3 + R21**2*(3*lamb6R*R22 - &
         3*lamb6I*R23) + (lamb7R*R22 - lamb7I*R23)*(R22**2 + R23**2) + &
         R21*((lamb3 + lamb4 + 2*lamb5R)*R22**2 - 4*lamb5I*R22*R23 + &
         (lamb3 + lamb4 - 2*lamb5R)*R23**2)) + R12*(3*lamb7R*R21*R22**2 &
         + 2*lamb2*R22**3 + R22*((lamb3 + lamb4 + 2*lamb5R)*R21**2 &
         - 2*lamb7I*R21*R23 + 2*lamb2*R23**2) + R21*(lamb6R*R21**2 - &
         2*lamb5I*R21*R23 + lamb7R*R23**2)))/2.0
    gh1223 = (R13*(R22**2*(-(lamb7I*R31) + 2*lamb2*R33) + &
         R23**2*(-3*lamb7I*R31 + 6*lamb2*R33) + R21**2*(-3*lamb6I*R31 - &
         2*lamb5I*R32 + (lamb3 + lamb4 - 2*lamb5R)*R33) + 2*R21*R23*((lamb3 &
         + lamb4 - 2*lamb5R)*R31 + lamb7R*R32 - 3*lamb7I*R33) + &
         2*R22*(-2*lamb5I*R21*R31 + lamb7R*R23*R31 - lamb7I*R21*R32 + &
         2*lamb2*R23*R32 + lamb7R*R21*R33)) + R11*(3*R21**2*(2*lamb1*R31 &
         + lamb6R*R32 - lamb6I*R33) + R23**2*((lamb3 + &
         lamb4 - 2*lamb5R)*R31 + lamb7R*R32 - 3*lamb7I*R33) + &
         R22**2*((lamb3 + lamb4 + 2*lamb5R)*R31 + 3*lamb7R*R32 - &
         lamb7I*R33) - 2*R22*R23*(2*lamb5I*R31 + lamb7I*R32 - &
         lamb7R*R33) + 2*R21*(R22*(3*lamb6R*R31 + (lamb3 + lamb4 &
         + 2*lamb5R)*R32 - 2*lamb5I*R33) + R23*(-3*lamb6I*R31 - &
         2*lamb5I*R32 + (lamb3 + lamb4 - 2*lamb5R)*R33))) + &
         R12*(3*R22**2*(lamb7R*R31 + 2*lamb2*R32) + R23**2*(lamb7R*R31 &
         + 2*lamb2*R32) + R21**2*(3*lamb6R*R31 + (lamb3 + lamb4 + &
         2*lamb5R)*R32 - 2*lamb5I*R33) - 2*R21*R23*(2*lamb5I*R31 + &
         lamb7I*R32 - lamb7R*R33) + 2*R22*(R23*(-(lamb7I*R31) + &
         2*lamb2*R33) + R21*((lamb3 + lamb4 + 2*lamb5R)*R31 + &
         3*lamb7R*R32 - lamb7I*R33))))/2.0
    gh1233 = (R12*(2*R23*(-(lamb5I*R31**2) - lamb7I*R31*R32 + &
         lamb7R*R31*R33 + 2*lamb2*R32*R33) + R22*((lamb3 + lamb4 + &
         2*lamb5R)*R31**2 + 6*lamb7R*R31*R32 + 6*lamb2*R32**2 - &
         2*lamb7I*R31*R33 + 2*lamb2*R33**2) + R21*(3*lamb6R*R31**2 + &
         3*lamb7R*R32**2 - 2*lamb7I*R32*R33 + lamb7R*R33**2 + &
         2*R31*((lamb3 + lamb4 + 2*lamb5R)*R32 - 2*lamb5I*R33))) + &
         R13*(2*R22*(-(lamb5I*R31**2) - lamb7I*R31*R32 + lamb7R*R31*R33 + &
         2*lamb2*R32*R33) + R23*((lamb3 + lamb4 - 2*lamb5R)*R31**2 + &
         2*lamb7R*R31*R32 + 2*lamb2*R32**2 - 6*lamb7I*R31*R33 + &
         6*lamb2*R33**2) - R21*(3*lamb6I*R31**2 + lamb7I*R32**2 - &
         2*lamb7R*R32*R33 + 3*lamb7I*R33**2 - &
         2*R31*(-2*lamb5I*R32 + (lamb3 + lamb4 - 2*lamb5R)*R33))) + &
         R11*(R22*(3*lamb6R*R31**2 + 3*lamb7R*R32**2 - 2*lamb7I*R32*R33 + &
         lamb7R*R33**2 + 2*R31*((lamb3 + lamb4 + 2*lamb5R)*R32 - &
         2*lamb5I*R33)) - R23*(3*lamb6I*R31**2 + lamb7I*R32**2 - &
         2*lamb7R*R32*R33 + 3*lamb7I*R33**2 - 2*R31*(-2*lamb5I*R32 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33)) + R21*(6*lamb1*R31**2 + &
         (lamb3 + lamb4 + 2*lamb5R)*R32**2 - 4*lamb5I*R32*R33 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33**2 + R31*(6*lamb6R*R32 - &
         6*lamb6I*R33))))/2.
    gh1333 = R13*(-(lamb6I*R31**3) + (lamb3 + lamb4 - 2*lamb5R)*R31**2*R33 &
         - 3*lamb7I*R31*R33**2 + 2*lamb2*R33**3 + R32**2*(-(lamb7I*R31) &
         + 2*lamb2*R33) + 2*R31*R32*(-(lamb5I*R31) + lamb7R*R33)) + &
         R11*(2*lamb1*R31**3 + R31**2*(3*lamb6R*R32 - 3*lamb6I*R33) + &
         (lamb7R*R32 - lamb7I*R33)*(R32**2 + R33**2) + R31*((lamb3 + &
         lamb4 + 2*lamb5R)*R32**2 - 4*lamb5I*R32*R33 + (lamb3 + lamb4 &
         - 2*lamb5R)*R33**2)) + R12*(3*lamb7R*R31*R32**2 + &
         2*lamb2*R32**3 + R32*((lamb3 + lamb4 + 2*lamb5R)*R31**2 - &
         2*lamb7I*R31*R33 + 2*lamb2*R33**2) + R31*(lamb6R*R31**2 - &
         2*lamb5I*R31*R33 + lamb7R*R33**2))
    gh2222 = (lamb1*R21**4 + R21**3*(2*lamb6R*R22 - 2*lamb6I*R23) + &
         2*R21*(lamb7R*R22 - lamb7I*R23)*(R22**2 + R23**2) + lamb2*(R22**2 + &
         R23**2)**2 + R21**2*((lamb3 + lamb4 + 2*lamb5R)*R22**2 - &
         4*lamb5I*R22*R23 + (lamb3 + lamb4 - 2*lamb5R)*R23**2))/4.0
    gh2223 = (R21**3*(2*lamb1*R31 + lamb6R*R32 - lamb6I*R33) + &
         (R22**2 + R23**2)*(lamb7R*R22*R31 - lamb7I*R23*R31 + &
         2*lamb2*R22*R32 + 2*lamb2*R23*R33) + R21**2*(R22*(3*lamb6R*R31 &
         + (lamb3 + lamb4 + 2*lamb5R)*R32 - 2*lamb5I*R33) + &
         R23*(-3*lamb6I*R31 - 2*lamb5I*R32 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33)) + R21*(R23**2*((lamb3 + lamb4 - &
         2*lamb5R)*R31 + lamb7R*R32 - 3*lamb7I*R33) + R22**2*((lamb3 &
         + lamb4 + 2*lamb5R)*R31 + 3*lamb7R*R32 - lamb7I*R33) - &
         2*R22*R23*(2*lamb5I*R31 + lamb7I*R32 - lamb7R*R33)))/2.0
    gh2233 = (4*R22*R23*(-(lamb5I*R31**2) - lamb7I*R31*R32 + &
         lamb7R*R31*R33 + 2*lamb2*R32*R33) + R22**2*((lamb3 + &
         lamb4 + 2*lamb5R)*R31**2 + 6*lamb7R*R31*R32 + 6*lamb2*R32**2 &
         - 2*lamb7I*R31*R33 + 2*lamb2*R33**2) + R23**2*((lamb3 + lamb4 &
         - 2*lamb5R)*R31**2 + 2*lamb7R*R31*R32 + 2*lamb2*R32**2 - &
         6*lamb7I*R31*R33 + 6*lamb2*R33**2) + R21**2*(6*lamb1*R31**2 + &
         (lamb3 + lamb4 + 2*lamb5R)*R32**2 - 4*lamb5I*R32*R33 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33**2 + R31*(6*lamb6R*R32 - &
         6*lamb6I*R33)) + 2*R21*(R22*(3*lamb6R*R31**2 + 3*lamb7R*R32**2 - &
         2*lamb7I*R32*R33 + lamb7R*R33**2 + 2*R31*((lamb3 + lamb4 + &
         2*lamb5R)*R32 - 2*lamb5I*R33)) - R23*(3*lamb6I*R31**2 + &
         lamb7I*R32**2 - 2*lamb7R*R32*R33 + 3*lamb7I*R33**2 - &
         2*R31*(-2*lamb5I*R32 + (lamb3 + lamb4 - 2*lamb5R)*R33))))/4.0
    gh2333 = (R23*(-(lamb6I*R31**3) + (lamb3 + lamb4 - &
         2*lamb5R)*R31**2*R33 - 3*lamb7I*R31*R33**2 + 2*lamb2*R33**3 &
         + R32**2*(-(lamb7I*R31) + 2*lamb2*R33) + 2*R31*R32*(-(lamb5I*R31) &
         + lamb7R*R33)) + R21*(2*lamb1*R31**3 + R31**2*(3*lamb6R*R32 &
         - 3*lamb6I*R33) + (lamb7R*R32 - lamb7I*R33)*(R32**2 + R33**2) + &
         R31*((lamb3 + lamb4 + 2*lamb5R)*R32**2 - 4*lamb5I*R32*R33 + &
         (lamb3 + lamb4 - 2*lamb5R)*R33**2)) + R22*(3*lamb7R*R31*R32**2 &
         + 2*lamb2*R32**3 + R32*((lamb3 + lamb4 + 2*lamb5R)*R31**2 - &
         2*lamb7I*R31*R33 + 2*lamb2*R33**2) + R31*(lamb6R*R31**2 - &
         2*lamb5I*R31*R33 + lamb7R*R33**2)))/2.0
    gh3333 = (lamb1*R31**4 + R31**3*(2*lamb6R*R32 - 2*lamb6I*R33) + &
         2*R31*(lamb7R*R32 - lamb7I*R33)*(R32**2 + R33**2) + lamb2*(R32**2 + &
         R33**2)**2 + R31**2*((lamb3 + lamb4 + 2*lamb5R)*R32**2 - &
         4*lamb5I*R32*R33 + (lamb3 + lamb4 - 2*lamb5R)*R33**2))/4.
    !!! Yukawa couplings. Up-type quark, up-type quark
    gh1uu(1) = -((0,-1)*R13*(Yu11 - conjg(Yu11)) + &
         R12*(Yu11 + conjg(Yu11)) + 2*R11*mass(2))/(2.*vev)
    gh2uu(1) = -((0,-1)*R23*(Yu11 - conjg(Yu11)) + &
         R22*(Yu11 + conjg(Yu11)) + 2*R21*mass(2))/(2.*vev)
    gh3uu(1) = -((0,-1)*R33*(Yu11 - conjg(Yu11)) + &
         R32*(Yu11 + conjg(Yu11)) + 2*R31*mass(2))/(2.*vev)
    gh1uc(1) = (imago*R13*(Yu12 - conjg(Yu21)) - &
         R12*(Yu12 + conjg(Yu21)))/(2.*vev)
    gh2uc(1) = (imago*R23*(Yu12 - conjg(Yu21)) - &
         R22*(Yu12 + conjg(Yu21)))/(2.*vev)
    gh3uc(1) = (imago*R33*(Yu12 - conjg(Yu21)) - &
         R32*(Yu12 + conjg(Yu21)))/(2.*vev)
    gh1ut(1) = (imago*R13*(Yu13 - conjg(Yu31)) - &
         R12*(Yu13 + conjg(Yu31)))/(2.*vev)
    gh2ut(1) = (imago*R23*(Yu13 - conjg(Yu31)) - &
         R22*(Yu13 + conjg(Yu31)))/(2.*vev)
    gh3ut(1) = (imago*R33*(Yu13 - conjg(Yu31)) - &
         R32*(Yu13 + conjg(Yu31)))/(2.*vev)
    gh1cu(1) = -(imago*R13*(-Yu21 + conjg(Yu12)) + &
         R12*(Yu21 + conjg(Yu12)))/(2.*vev)
    gh2cu(1) = -(imago*R23*(-Yu21 + conjg(Yu12)) + &
         R22*(Yu21 + conjg(Yu12)))/(2.*vev)
    gh3cu(1) = -(imago*R33*(-Yu21 + conjg(Yu12)) + &
         R32*(Yu21 + conjg(Yu12)))/(2.*vev)
    gh1cc(1) = -((0,-1)*R13*(Yu22 - conjg(Yu22)) + &
         R12*(Yu22 + conjg(Yu22)) + 2*R11*mass(4))/(2.*vev)
    gh2cc(1) = -((0,-1)*R23*(Yu22 - conjg(Yu22)) + &
         R22*(Yu22 + conjg(Yu22)) + 2*R21*mass(4))/(2.*vev)
    gh3cc(1) = -((0,-1)*R33*(Yu22 - conjg(Yu22)) + &
         R32*(Yu22 + conjg(Yu22)) + 2*R31*mass(4))/(2.*vev)
    gh1ct(1) = (imago*R13*(Yu23 - conjg(Yu32)) - &
         R12*(Yu23 + conjg(Yu32)))/(2.*vev)
    gh2ct(1) = (imago*R23*(Yu23 - conjg(Yu32)) - &
         R22*(Yu23 + conjg(Yu32)))/(2.*vev)
    gh3ct(1) = (imago*R33*(Yu23 - conjg(Yu32)) - &
         R32*(Yu23 + conjg(Yu32)))/(2.*vev)
    gh1tu(1) = -(imago*R13*(-Yu31 + conjg(Yu13)) + &
         R12*(Yu31 + conjg(Yu13)))/(2.*vev)
    gh2tu(1) = -(imago*R23*(-Yu31 + conjg(Yu13)) + &
         R22*(Yu31 + conjg(Yu13)))/(2.*vev)
    gh3tu(1) = -(imago*R33*(-Yu31 + conjg(Yu13)) + &
         R32*(Yu31 + conjg(Yu13)))/(2.*vev)
    gh1tc(1) = -(imago*R13*(-Yu32 + conjg(Yu23)) + &
         R12*(Yu32 + conjg(Yu23)))/(2.*vev)
    gh2tc(1) = -(imago*R23*(-Yu32 + conjg(Yu23)) + &
         R22*(Yu32 + conjg(Yu23)))/(2.*vev)
    gh3tc(1) = -(imago*R33*(-Yu32 + conjg(Yu23)) + &
         R32*(Yu32 + conjg(Yu23)))/(2.*vev)
    gh1tt(1) = -((0,-1)*R13*(Yu33 - conjg(Yu33)) + &
         R12*(Yu33 + conjg(Yu33)) + 2*R11*mass(6))/(2.*vev)
    gh2tt(1) = -((0,-1)*R23*(Yu33 - conjg(Yu33)) + &
         R22*(Yu33 + conjg(Yu33)) + 2*R21*mass(6))/(2.*vev)
    gh3tt(1) = -((0,-1)*R33*(Yu33 - conjg(Yu33)) + &
         R32*(Yu33 + conjg(Yu33)) + 2*R31*mass(6))/(2.*vev)
    gh1uu(2) = (R12*(-Yu11 + conjg(Yu11)) + &
         imago*R13*(Yu11 + conjg(Yu11)))/(2.*vev)
    gh2uu(2) = (R22*(-Yu11 + conjg(Yu11)) + &
         imago*R23*(Yu11 + conjg(Yu11)))/(2.*vev)
    gh3uu(2) = (R32*(-Yu11 + conjg(Yu11)) + &
         imago*R33*(Yu11 + conjg(Yu11)))/(2.*vev)
    gh1uc(2) = (R12*(-Yu12 + conjg(Yu21)) + &
         imago*R13*(Yu12 + conjg(Yu21)))/(2.*vev)
    gh2uc(2) = (R22*(-Yu12 + conjg(Yu21)) + &
         imago*R23*(Yu12 + conjg(Yu21)))/(2.*vev)
    gh3uc(2) = (R32*(-Yu12 + conjg(Yu21)) + &
         imago*R33*(Yu12 + conjg(Yu21)))/(2.*vev)
    gh1ut(2) = (R12*(-Yu13 + conjg(Yu31)) + &
         imago*R13*(Yu13 + conjg(Yu31)))/(2.*vev)
    gh2ut(2) = (R22*(-Yu13 + conjg(Yu31)) + &
         imago*R23*(Yu13 + conjg(Yu31)))/(2.*vev)
    gh3ut(2) = (R32*(-Yu13 + conjg(Yu31)) + &
         imago*R33*(Yu13 + conjg(Yu31)))/(2.*vev)
    gh1cu(2) = (R12*(-Yu21 + conjg(Yu12)) + &
         imago*R13*(Yu21 + conjg(Yu12)))/(2.*vev)
    gh2cu(2) = (R22*(-Yu21 + conjg(Yu12)) + &
         imago*R23*(Yu21 + conjg(Yu12)))/(2.*vev)
    gh3cu(2) = (R32*(-Yu21 + conjg(Yu12)) + &
         imago*R33*(Yu21 + conjg(Yu12)))/(2.*vev)
    gh1cc(2) = (R12*(-Yu22 + conjg(Yu22)) + &
         imago*R13*(Yu22 + conjg(Yu22)))/(2.*vev)
    gh2cc(2) = (R22*(-Yu22 + conjg(Yu22)) + &
         imago*R23*(Yu22 + conjg(Yu22)))/(2.*vev)
    gh3cc(2) = (R32*(-Yu22 + conjg(Yu22)) + &
         imago*R33*(Yu22 + conjg(Yu22)))/(2.*vev)
    gh1ct(2) = (R12*(-Yu23 + conjg(Yu32)) + &
         imago*R13*(Yu23 + conjg(Yu32)))/(2.*vev)
    gh2ct(2) = (R22*(-Yu23 + conjg(Yu32)) + &
         imago*R23*(Yu23 + conjg(Yu32)))/(2.*vev)
    gh3ct(2) = (R32*(-Yu23 + conjg(Yu32)) + &
         imago*R33*(Yu23 + conjg(Yu32)))/(2.*vev)
    gh1tu(2) = (R12*(-Yu31 + conjg(Yu13)) + &
         imago*R13*(Yu31 + conjg(Yu13)))/(2.*vev)
    gh2tu(2) = (R22*(-Yu31 + conjg(Yu13)) + &
         imago*R23*(Yu31 + conjg(Yu13)))/(2.*vev)
    gh3tu(2) = (R32*(-Yu31 + conjg(Yu13)) + &
         imago*R33*(Yu31 + conjg(Yu13)))/(2.*vev)
    gh1tc(2) = (R12*(-Yu32 + conjg(Yu23)) + &
         imago*R13*(Yu32 + conjg(Yu23)))/(2.*vev)
    gh2tc(2) = (R22*(-Yu32 + conjg(Yu23)) + &
         imago*R23*(Yu32 + conjg(Yu23)))/(2.*vev)
    gh3tc(2) = (R32*(-Yu32 + conjg(Yu23)) + &
         imago*R33*(Yu32 + conjg(Yu23)))/(2.*vev)
    gh1tt(2) = (R12*(-Yu33 + conjg(Yu33)) + &
         imago*R13*(Yu33 + conjg(Yu33)))/(2.*vev)
    gh2tt(2) = (R22*(-Yu33 + conjg(Yu33)) + &
         imago*R23*(Yu33 + conjg(Yu33)))/(2.*vev)
    gh3tt(2) = (R32*(-Yu33 + conjg(Yu33)) + &
         imago*R33*(Yu33 + conjg(Yu33)))/(2.*vev)
    !!! Yukawa couplings. Down-type quark, down-type quark
    gh1dd(1) = -(((imago*R13*(Yd11 - conjg(Yd11)) + &
         R12*(Yd11 + conjg(Yd11)))/2. + R11*mass(1))/vev)
    gh2dd(1) = -(((imago*R23*(Yd11 - conjg(Yd11)) + &
         R22*(Yd11 + conjg(Yd11)))/2. + R21*mass(1))/vev)
    gh3dd(1) = -(((imago*R33*(Yd11 - conjg(Yd11)) + &
         R32*(Yd11 + conjg(Yd11)))/2. + R31*mass(1))/vev)
    gh1ds(1) = -(imago*R13*(Yd12 - conjg(Yd21)) + &
         R12*(Yd12 + conjg(Yd21)))/(2.*vev)
    gh2ds(1) = -(imago*R23*(Yd12 - conjg(Yd21)) + &
         R22*(Yd12 + conjg(Yd21)))/(2.*vev)
    gh3ds(1) = -(imago*R33*(Yd12 - conjg(Yd21)) + &
         R32*(Yd12 + conjg(Yd21)))/(2.*vev)
    gh1db(1) = -(imago*R13*(Yd13 - conjg(Yd31)) + &
         R12*(Yd13 + conjg(Yd31)))/(2.*vev)
    gh2db(1) = -(imago*R23*(Yd13 - conjg(Yd31)) + &
         R22*(Yd13 + conjg(Yd31)))/(2.*vev)
    gh3db(1) = -(imago*R33*(Yd13 - conjg(Yd31)) + &
         R32*(Yd13 + conjg(Yd31)))/(2.*vev)
    gh1sd(1) = (imago*R13*(-Yd21 + conjg(Yd12)) - &
         R12*(Yd21 + conjg(Yd12)))/(2.*vev)
    gh2sd(1) = (imago*R23*(-Yd21 + conjg(Yd12)) - &
         R22*(Yd21 + conjg(Yd12)))/(2.*vev)
    gh3sd(1) = (imago*R33*(-Yd21 + conjg(Yd12)) - &
         R32*(Yd21 + conjg(Yd12)))/(2.*vev)
    gh1ss(1) = -(((imago*R13*(Yd22 - conjg(Yd22)) + &
         R12*(Yd22 + conjg(Yd22)))/2. + R11*mass(3))/vev)
    gh2ss(1) = -(((imago*R23*(Yd22 - conjg(Yd22)) + &
         R22*(Yd22 + conjg(Yd22)))/2. + R21*mass(3))/vev)
    gh3ss(1) = -(((imago*R33*(Yd22 - conjg(Yd22)) + &
         R32*(Yd22 + conjg(Yd22)))/2. + R31*mass(3))/vev)
    gh1sb(1) = -(imago*R13*(Yd23 - conjg(Yd32)) + &
         R12*(Yd23 + conjg(Yd32)))/(2.*vev)
    gh2sb(1) = -(imago*R23*(Yd23 - conjg(Yd32)) + &
         R22*(Yd23 + conjg(Yd32)))/(2.*vev)
    gh3sb(1) = -(imago*R33*(Yd23 - conjg(Yd32)) + &
         R32*(Yd23 + conjg(Yd32)))/(2.*vev)
    gh1bd(1) = (imago*R13*(-Yd31 + conjg(Yd13)) - &
         R12*(Yd31 + conjg(Yd13)))/(2.*vev)
    gh2bd(1) = (imago*R23*(-Yd31 + conjg(Yd13)) - &
         R22*(Yd31 + conjg(Yd13)))/(2.*vev)
    gh3bd(1) = (imago*R33*(-Yd31 + conjg(Yd13)) - &
         R32*(Yd31 + conjg(Yd13)))/(2.*vev)
    gh1bs(1) = (imago*R13*(-Yd32 + conjg(Yd23)) - &
         R12*(Yd32 + conjg(Yd23)))/(2.*vev)
    gh2bs(1) = (imago*R23*(-Yd32 + conjg(Yd23)) - &
         R22*(Yd32 + conjg(Yd23)))/(2.*vev)
    gh3bs(1) = (imago*R33*(-Yd32 + conjg(Yd23)) - &
         R32*(Yd32 + conjg(Yd23)))/(2.*vev)
    gh1bb(1) = -(((imago*R13*(Yd33 - conjg(Yd33)) + &
         R12*(Yd33 + conjg(Yd33)))/2. + R11*mass(5))/vev)
    gh2bb(1) = -(((imago*R23*(Yd33 - conjg(Yd33)) + &
         R22*(Yd33 + conjg(Yd33)))/2. + R21*mass(5))/vev)
    gh3bb(1) = -(((imago*R33*(Yd33 - conjg(Yd33)) + &
         R32*(Yd33 + conjg(Yd33)))/2. + R31*mass(5))/vev)
    gh1dd(2) = (R12*(-Yd11 + conjg(Yd11)) - &
         imago*R13*(Yd11 + conjg(Yd11)))/(2.*vev)
    gh2dd(2) = (R22*(-Yd11 + conjg(Yd11)) - &
         imago*R23*(Yd11 + conjg(Yd11)))/(2.*vev)
    gh3dd(2) = (R32*(-Yd11 + conjg(Yd11)) - &
         imago*R33*(Yd11 + conjg(Yd11)))/(2.*vev)
    gh1ds(2) = (R12*(-Yd12 + conjg(Yd21)) - &
         imago*R13*(Yd12 + conjg(Yd21)))/(2.*vev)
    gh2ds(2) = (R22*(-Yd12 + conjg(Yd21)) - &
         imago*R23*(Yd12 + conjg(Yd21)))/(2.*vev)
    gh3ds(2) = (R32*(-Yd12 + conjg(Yd21)) - &
         imago*R33*(Yd12 + conjg(Yd21)))/(2.*vev)
    gh1db(2) = (R12*(-Yd13 + conjg(Yd31)) - &
         imago*R13*(Yd13 + conjg(Yd31)))/(2.*vev)
    gh2db(2) = (R22*(-Yd13 + conjg(Yd31)) - &
         imago*R23*(Yd13 + conjg(Yd31)))/(2.*vev)
    gh3db(2) = (R32*(-Yd13 + conjg(Yd31)) - &
         imago*R33*(Yd13 + conjg(Yd31)))/(2.*vev)
    gh1sd(2) = (R12*(-Yd21 + conjg(Yd12)) - &
         imago*R13*(Yd21 + conjg(Yd12)))/(2.*vev)
    gh2sd(2) = (R22*(-Yd21 + conjg(Yd12)) - &
         imago*R23*(Yd21 + conjg(Yd12)))/(2.*vev)
    gh3sd(2) = (R32*(-Yd21 + conjg(Yd12)) - &
         imago*R33*(Yd21 + conjg(Yd12)))/(2.*vev)
    gh1ss(2) = (R12*(-Yd22 + conjg(Yd22)) - &
         imago*R13*(Yd22 + conjg(Yd22)))/(2.*vev)
    gh2ss(2) = (R22*(-Yd22 + conjg(Yd22)) - &
         imago*R23*(Yd22 + conjg(Yd22)))/(2.*vev)
    gh3ss(2) = (R32*(-Yd22 + conjg(Yd22)) - &
         imago*R33*(Yd22 + conjg(Yd22)))/(2.*vev)
    gh1sb(2) = (R12*(-Yd23 + conjg(Yd32)) - &
         imago*R13*(Yd23 + conjg(Yd32)))/(2.*vev)
    gh2sb(2) = (R22*(-Yd23 + conjg(Yd32)) - &
         imago*R23*(Yd23 + conjg(Yd32)))/(2.*vev)
    gh3sb(2) = (R32*(-Yd23 + conjg(Yd32)) - &
         imago*R33*(Yd23 + conjg(Yd32)))/(2.*vev)
    gh1bd(2) = (R12*(-Yd31 + conjg(Yd13)) - &
         imago*R13*(Yd31 + conjg(Yd13)))/(2.*vev)
    gh2bd(2) = (R22*(-Yd31 + conjg(Yd13)) - &
         imago*R23*(Yd31 + conjg(Yd13)))/(2.*vev)
    gh3bd(2) = (R32*(-Yd31 + conjg(Yd13)) - &
         imago*R33*(Yd31 + conjg(Yd13)))/(2.*vev)
    gh1bs(2) = (R12*(-Yd32 + conjg(Yd23)) - &
         imago*R13*(Yd32 + conjg(Yd23)))/(2.*vev)
    gh2bs(2) = (R22*(-Yd32 + conjg(Yd23)) - &
         imago*R23*(Yd32 + conjg(Yd23)))/(2.*vev)
    gh3bs(2) = (R32*(-Yd32 + conjg(Yd23)) - &
         imago*R33*(Yd32 + conjg(Yd23)))/(2.*vev)
    gh1bb(2) = (R12*(-Yd33 + conjg(Yd33)) - &
         imago*R13*(Yd33 + conjg(Yd33)))/(2.*vev)
    gh2bb(2) = (R22*(-Yd33 + conjg(Yd33)) - &
         imago*R23*(Yd33 + conjg(Yd33)))/(2.*vev)
    gh3bb(2) = (R32*(-Yd33 + conjg(Yd33)) - &
         imago*R33*(Yd33 + conjg(Yd33)))/(2.*vev)
    !!! Yukawa couplings: up quark down quark Higgs Plus
    ghud(1) = -(sqrt2*(V12*Yd21 + V13*Yd31 + V11*(Yd11 - conjg(Yu11)) - &
         V21*conjg(Yu21) - V31*conjg(Yu31)))/(2.*vev)
    ghus(1) = -(sqrt2*(V11*Yd12 + V13*Yd32 + V12*(Yd22 - conjg(Yu11)) - &
         V22*conjg(Yu21) - V32*conjg(Yu31)))/(2.*vev)
    ghub(1) = -(sqrt2*(V11*Yd13 + V12*Yd23 + V13*Yd33 - V13*conjg(Yu11) - &
         V23*conjg(Yu21) - V33*conjg(Yu31)))/(2.*vev)
    ghcd(1) = -(sqrt2*(V22*Yd21 + V23*Yd31 - V11*conjg(Yu12) + V21*(Yd11 - &
         conjg(Yu22)) - V31*conjg(Yu32)))/(2.*vev)
    ghcs(1) = -(sqrt2*(V21*Yd12 + V23*Yd32 - V12*conjg(Yu12) + V22*(Yd22 - &
         conjg(Yu22)) - V32*conjg(Yu32)))/(2.*vev)
    ghcb(1) = -(sqrt2*(V21*Yd13 + V22*Yd23 + V23*Yd33 - V13*conjg(Yu12) - &
         V23*conjg(Yu22) - V33*conjg(Yu32)))/(2.*vev)
    ghtd(1) = -(sqrt2*(V32*Yd21 + V33*Yd31 - V11*conjg(Yu13) - &
         V21*conjg(Yu23) + V31*(Yd11 - conjg(Yu33))))/(2.*vev)
    ghts(1) = -(sqrt2*(V31*Yd12 + V33*Yd32 - V12*conjg(Yu13) - &
         V22*conjg(Yu23) + V32*(Yd22 - conjg(Yu33))))/(2.*vev)
    ghtb(1) = -(sqrt2*(V31*Yd13 + V32*Yd23 + V33*Yd33 - V13*conjg(Yu13) - &
         V23*conjg(Yu23) - V33*conjg(Yu33)))/(2.*vev)
    ghud(2) = -(sqrt2*(V12*Yd21 + V13*Yd31 + V11*(Yd11 + conjg(Yu11)) + &
         V21*conjg(Yu21) + V31*conjg(Yu31)))/(2.*vev)
    ghus(2) = -(sqrt2*(V11*Yd12 + V13*Yd32 + V12*(Yd22 + conjg(Yu11)) + &
         V22*conjg(Yu21) + V32*conjg(Yu31)))/(2.*vev)
    ghub(2) = -(sqrt2*(V11*Yd13 + V12*Yd23 + V13*Yd33 + V13*conjg(Yu11) + &
         V23*conjg(Yu21) + V33*conjg(Yu31)))/(2.*vev)
    ghcd(2) = -(sqrt2*(V22*Yd21 + V23*Yd31 + V11*conjg(Yu12) + V21*(Yd11 + &
         conjg(Yu22)) + V31*conjg(Yu32)))/(2.*vev)
    ghcs(2) = -(sqrt2*(V21*Yd12 + V23*Yd32 + V12*conjg(Yu12) + V22*(Yd22 + &
         conjg(Yu22)) + V32*conjg(Yu32)))/(2.*vev)
    ghcb(2) = -(sqrt2*(V21*Yd13 + V22*Yd23 + V23*Yd33 + V13*conjg(Yu12) + &
         V23*conjg(Yu22) + V33*conjg(Yu32)))/(2.*vev)
    ghtd(2) = -(sqrt2*(V32*Yd21 + V33*Yd31 + V11*conjg(Yu13) + &
         V21*conjg(Yu23) + V31*(Yd11 + conjg(Yu33))))/(2.*vev)
    ghts(2) = -(sqrt2*(V31*Yd12 + V33*Yd32 + V12*conjg(Yu13) + &
         V22*conjg(Yu23) + V32*(Yd22 + conjg(Yu33))))/(2.*vev)
    ghtb(2) = -(sqrt2*(V31*Yd13 + V32*Yd23 + V33*Yd33 + &
         V13*conjg(Yu13) + V23*conjg(Yu23) + V33*conjg(Yu33)))/(2.*vev)
    !!! Yukawa couplings: down quark up quark Higgs Minus
    ghdu(1) = (sqrt2*(Yu21*conjg(V12) + Yu31*conjg(V13) + &
         conjg(V11)*(Yu11 - conjg(Yd11)) - conjg(V21)*conjg(Yd21) - &
         conjg(V31)*conjg(Yd31)))/(2.*vev)
    ghdc(1) = (sqrt2*(Yu12*conjg(V11) + Yu32*conjg(V13) + conjg(V12) * &
         (Yu22 - conjg(Yd11)) - conjg(V22)*conjg(Yd21) - &
         conjg(V32)*conjg(Yd31)))/(2.*vev)
    ghdt(1) = (sqrt2*(Yu13*conjg(V11) + Yu23*conjg(V12) + Yu33*conjg(V13) - &
         conjg(V13)*conjg(Yd11) - conjg(V23)*conjg(Yd21) - conjg(V33) * &
         conjg(Yd31)))/(2.*vev)
    ghsu(1) = (sqrt2*(Yu21*conjg(V22) + Yu31*conjg(V23) - conjg(V11) * &
         conjg(Yd12) + conjg(V21)*(Yu11 - conjg(Yd22)) - conjg(V31) * &
         conjg(Yd32)))/(2.*vev)
    ghsc(1) = (sqrt2*(Yu12*conjg(V21) + Yu22*conjg(V22) + Yu32*conjg(V23) - &
         conjg(V12)*conjg(Yd12) - conjg(V22)*conjg(Yd22) - conjg(V32) * &
         conjg(Yd32)))/(2.*vev)
    ghst(1) = (sqrt2*(Yu13*conjg(V21) + Yu23*conjg(V22) + Yu33*conjg(V23) - &
         conjg(V13)*conjg(Yd12) - conjg(V23)*conjg(Yd22) - conjg(V33) * &
         conjg(Yd32)))/(2.*vev)
    ghbu(1) = (sqrt2*(Yu11*conjg(V31) + Yu21*conjg(V32) + Yu31*conjg(V33) - &
         conjg(V11)*conjg(Yd13) - conjg(V21)*conjg(Yd23) - conjg(V31) * &
         conjg(Yd33)))/(2.*vev)
    ghbc(1) = (sqrt2*(Yu12*conjg(V31) + Yu22*conjg(V32) + Yu32*conjg(V33) - &
         conjg(V12)*conjg(Yd13) - conjg(V22)*conjg(Yd23) - conjg(V32) * &
         conjg(Yd33)))/(2.*vev)
    ghbt(1) = (sqrt2*(Yu13*conjg(V31) + Yu23*conjg(V32) + Yu33*conjg(V33) - &
         conjg(V13)*conjg(Yd13) - conjg(V23)*conjg(Yd23) - conjg(V33) * &
         conjg(Yd33)))/(2.*vev)
    ghdu(2) = (sqrt2*(Yu21*conjg(V12) + Yu31*conjg(V13) + conjg(V11)*(Yu11 + &
         conjg(Yd11)) + conjg(V21)*conjg(Yd21) + conjg(V31) * &
         conjg(Yd31)))/(2.*vev)
    ghdc(2) = (sqrt2*(Yu12*conjg(V11) + Yu32*conjg(V13) + conjg(V12)*(Yu22 + &
         conjg(Yd11)) + conjg(V22)*conjg(Yd21) + conjg(V32) * &
         conjg(Yd31)))/(2.*vev)
    ghdt(2) = (sqrt2*(Yu13*conjg(V11) + Yu23*conjg(V12) + Yu33*conjg(V13) + &
         conjg(V13)*conjg(Yd11) + conjg(V23)*conjg(Yd21) + conjg(V33) * &
         conjg(Yd31)))/(2.*vev)
    ghsu(2) = (sqrt2*(Yu21*conjg(V22) + Yu31*conjg(V23) + conjg(V11) * &
         conjg(Yd12) + conjg(V21)*(Yu11 + conjg(Yd22)) + conjg(V31) * &
         conjg(Yd32)))/(2.*vev)
    ghsc(2) = (sqrt2*(Yu12*conjg(V21) + Yu22*conjg(V22) + Yu32*conjg(V23) + &
         conjg(V12)*conjg(Yd12) + conjg(V22)*conjg(Yd22) + conjg(V32) * &
         conjg(Yd32)))/(2.*vev)
    ghst(2) = (sqrt2*(Yu13*conjg(V21) + Yu23*conjg(V22) + Yu33*conjg(V23) + &
         conjg(V13)*conjg(Yd12) + conjg(V23)*conjg(Yd22) + conjg(V33) * &
         conjg(Yd32)))/(2.*vev)
    ghbu(2) = (sqrt2*(Yu11*conjg(V31) + Yu21*conjg(V32) + Yu31*conjg(V33) + &
         conjg(V11)*conjg(Yd13) + conjg(V21)*conjg(Yd23) + conjg(V31) * &
         conjg(Yd33)))/(2.*vev)
    ghbc(2) = (sqrt2*(Yu12*conjg(V31) + Yu22*conjg(V32) + Yu32*conjg(V33) + &
         conjg(V12)*conjg(Yd13) + conjg(V22)*conjg(Yd23) + conjg(V32) * &
         conjg(Yd33)))/(2.*vev)
    ghbt(2) = (sqrt2*(Yu13*conjg(V31) + Yu23*conjg(V32) + Yu33*conjg(V33) + &
         conjg(V13)*conjg(Yd13) + conjg(V23)*conjg(Yd23) + conjg(V33) * &
         conjg(Yd33)))/(2.*vev)
    !!! Yukawa couplings: Charged Higgs + Neutrino + Electron
    ghe1n1 = - (conjg(Yl11)/(sqrt2*vev))
    ghe1n2 = - (conjg(Yl21)/(sqrt2*vev))
    ghe1n3 = - (conjg(Yl31)/(sqrt2*vev))
    ghe2n1 = - (conjg(Yl12)/(sqrt2*vev))
    ghe2n2 = - (conjg(Yl22)/(sqrt2*vev))
    ghe2n3 = - (conjg(Yl32)/(sqrt2*vev))
    ghe3n1 = - (conjg(Yl13)/(sqrt2*vev))
    ghe3n2 = - (conjg(Yl23)/(sqrt2*vev))
    ghe3n3 = - (conjg(Yl33)/(sqrt2*vev))
    ghn1e1 = - (Yl11/(sqrt2*vev))
    ghn1e2 = - (Yl12/(sqrt2*vev))
    ghn1e3 = - (Yl13/(sqrt2*vev))
    ghn2e1 = - (Yl21/(sqrt2*vev))
    ghn2e2 = - (Yl22/(sqrt2*vev))
    ghn2e3 = - (Yl23/(sqrt2*vev))
    ghn3e1 = - (Yl31/(sqrt2*vev))
    ghn3e2 = - (Yl32/(sqrt2*vev))
    ghn3e3 = - (Yl33/(sqrt2*vev))
    gh1e1e1(1) = - (((imago*R13*(Yl11 - conjg(Yl11)) + &
         R12*(Yl11 + conjg(Yl11)))/2. + R11*mass(11))/vev)
    gh2e1e1(1) = -(((imago*R23*(Yl11 - conjg(Yl11)) + &
         R22*(Yl11 + conjg(Yl11)))/2. + R21*mass(11))/vev)
    gh3e1e1(1) = -(((imago*R33*(Yl11 - conjg(Yl11)) + &
         R32*(Yl11 + conjg(Yl11)))/2. + R31*mass(11))/vev)
    gh1e1e2(1) = -(imago*R13*(Yl12 - conjg(Yl21)) + &
         R12*(Yl12 + conjg(Yl21)))/(2.*vev)
    gh2e1e2(1) = -(imago*R23*(Yl12 - conjg(Yl21)) + &
         R22*(Yl12 + conjg(Yl21)))/(2.*vev)
    gh3e1e2(1) = -(imago*R33*(Yl12 - conjg(Yl21)) + &
         R32*(Yl12 + conjg(Yl21)))/(2.*vev)
    gh1e1e3(1) = -(imago*R13*(Yl13 - conjg(Yl31)) + &
         R12*(Yl13 + conjg(Yl31)))/(2.*vev)
    gh2e1e3(1) = -(imago*R23*(Yl13 - conjg(Yl31)) + &
         R22*(Yl13 + conjg(Yl31)))/(2.*vev)
    gh3e1e3(1) = -(imago*R33*(Yl13 - conjg(Yl31)) + &
         R32*(Yl13 + conjg(Yl31)))/(2.*vev)
    gh1e2e1(1) = (imago*R13*(-Yl21 + conjg(Yl12)) - &
         R12*(Yl21 + conjg(Yl12)))/(2.*vev)
    gh2e2e1(1) = (imago*R23*(-Yl21 + conjg(Yl12)) - &
         R22*(Yl21 + conjg(Yl12)))/(2.*vev)
    gh3e2e1(1) = (imago*R33*(-Yl21 + conjg(Yl12)) - &
         R32*(Yl21 + conjg(Yl12)))/(2.*vev)
    gh1e2e2(1) = -(((imago*R13*(Yl22 - conjg(Yl22)) + &
         R12*(Yl22 + conjg(Yl22)))/2. + R11*mass(13))/vev)
    gh2e2e2(1) = -(((imago*R23*(Yl22 - conjg(Yl22)) + &
         R22*(Yl22 + conjg(Yl22)))/2. + R21*mass(13))/vev)
    gh3e2e2(1) = -(((imago*R33*(Yl22 - conjg(Yl22)) + &
         R32*(Yl22 + conjg(Yl22)))/2. + R31*mass(13))/vev)
    gh1e2e3(1) = -(imago*R13*(Yl23 - conjg(Yl32)) + &
         R12*(Yl23 + conjg(Yl32)))/(2.*vev)
    gh2e2e3(1) = -(imago*R23*(Yl23 - conjg(Yl32)) + &
         R22*(Yl23 + conjg(Yl32)))/(2.*vev)
    gh3e2e3(1) = -(imago*R33*(Yl23 - conjg(Yl32)) + &
         R32*(Yl23 + conjg(Yl32)))/(2.*vev)
    gh1e3e1(1) = (imago*R13*(-Yl31 + conjg(Yl13)) - &
         R12*(Yl31 + conjg(Yl13)))/(2.*vev)
    gh2e3e1(1) = (imago*R23*(-Yl31 + conjg(Yl13)) - &
         R22*(Yl31 + conjg(Yl13)))/(2.*vev)
    gh3e3e1(1) = (imago*R33*(-Yl31 + conjg(Yl13)) - &
         R32*(Yl31 + conjg(Yl13)))/(2.*vev)
    gh1e3e2(1) = (imago*R13*(-Yl32 + conjg(Yl23)) - &
         R12*(Yl32 + conjg(Yl23)))/(2.*vev)
    gh2e3e2(1) = (imago*R23*(-Yl32 + conjg(Yl23)) - &
         R22*(Yl32 + conjg(Yl23)))/(2.*vev)
    gh3e3e2(1) = (imago*R33*(-Yl32 + conjg(Yl23)) - &
         R32*(Yl32 + conjg(Yl23)))/(2.*vev)
    gh1e3e3(1) = -(((imago*R13*(Yl33 - conjg(Yl33)) + &
         R12*(Yl33 + conjg(Yl33)))/2. + R11*mass(15))/vev)
    gh2e3e3(1) = -(((imago*R23*(Yl33 - conjg(Yl33)) + &
         R22*(Yl33 + conjg(Yl33)))/2. + R21*mass(15))/vev)
    gh3e3e3(1) = -(((imago*R33*(Yl33 - conjg(Yl33)) + &
         R32*(Yl33 + conjg(Yl33)))/2. + R31*mass(15))/vev)
    gh1e1e1(2) = (R12*(-Yl11 + conjg(Yl11)) - &
         imago*R13*(Yl11 + conjg(Yl11)))/(2.*vev)
    gh2e1e1(2) = (R22*(-Yl11 + conjg(Yl11)) - &
         imago*R23*(Yl11 + conjg(Yl11)))/(2.*vev)
    gh3e1e1(2) = (R32*(-Yl11 + conjg(Yl11)) - &
         imago*R33*(Yl11 + conjg(Yl11)))/(2.*vev)
    gh1e1e2(2) = (R12*(-Yl12 + conjg(Yl21)) - &
         imago*R13*(Yl12 + conjg(Yl21)))/(2.*vev)
    gh2e1e2(2) = (R22*(-Yl12 + conjg(Yl21)) - &
         imago*R23*(Yl12 + conjg(Yl21)))/(2.*vev)
    gh3e1e2(2) = (R32*(-Yl12 + conjg(Yl21)) - &
         imago*R33*(Yl12 + conjg(Yl21)))/(2.*vev)
    gh1e1e3(2) = (R12*(-Yl13 + conjg(Yl31)) - &
         imago*R13*(Yl13 + conjg(Yl31)))/(2.*vev)
    gh2e1e3(2) = (R22*(-Yl13 + conjg(Yl31)) - &
         imago*R23*(Yl13 + conjg(Yl31)))/(2.*vev)
    gh3e1e3(2) = (R32*(-Yl13 + conjg(Yl31)) - &
         imago*R33*(Yl13 + conjg(Yl31)))/(2.*vev)
    gh1e2e1(2) = (R12*(-Yl21 + conjg(Yl12)) - &
         imago*R13*(Yl21 + conjg(Yl12)))/(2.*vev)
    gh2e2e1(2) = (R22*(-Yl21 + conjg(Yl12)) - &
         imago*R23*(Yl21 + conjg(Yl12)))/(2.*vev)
    gh3e2e1(2) = (R32*(-Yl21 + conjg(Yl12)) - &
         imago*R33*(Yl21 + conjg(Yl12)))/(2.*vev)
    gh1e2e2(2) = (R12*(-Yl22 + conjg(Yl22)) - &
         imago*R13*(Yl22 + conjg(Yl22)))/(2.*vev)
    gh2e2e2(2) = (R22*(-Yl22 + conjg(Yl22)) - &
         imago*R23*(Yl22 + conjg(Yl22)))/(2.*vev)
    gh3e2e2(2) = (R32*(-Yl22 + conjg(Yl22)) - &
         imago*R33*(Yl22 + conjg(Yl22)))/(2.*vev)
    gh1e2e3(2) = (R12*(-Yl23 + conjg(Yl32)) - &
         imago*R13*(Yl23 + conjg(Yl32)))/(2.*vev)
    gh2e2e3(2) = (R22*(-Yl23 + conjg(Yl32)) - &
         imago*R23*(Yl23 + conjg(Yl32)))/(2.*vev)
    gh3e2e3(2) = (R32*(-Yl23 + conjg(Yl32)) - &
         imago*R33*(Yl23 + conjg(Yl32)))/(2.*vev)
    gh1e3e1(2) = (R12*(-Yl31 + conjg(Yl13)) - &
         imago*R13*(Yl31 + conjg(Yl13)))/(2.*vev)
    gh2e3e1(2) = (R22*(-Yl31 + conjg(Yl13)) - &
         imago*R23*(Yl31 + conjg(Yl13)))/(2.*vev)
    gh3e3e1(2) = (R32*(-Yl31 + conjg(Yl13)) - &
         imago*R33*(Yl31 + conjg(Yl13)))/(2.*vev)
    gh1e3e2(2) = (R12*(-Yl32 + conjg(Yl23)) - &
         imago*R13*(Yl32 + conjg(Yl23)))/(2.*vev)
    gh2e3e2(2) = (R22*(-Yl32 + conjg(Yl23)) - &
         imago*R23*(Yl32 + conjg(Yl23)))/(2.*vev)
    gh3e3e2(2) = (R32*(-Yl32 + conjg(Yl23)) - &
         imago*R33*(Yl32 + conjg(Yl23)))/(2.*vev)
    gh1e3e3(2) = (R12*(-Yl33 + conjg(Yl33)) - &
         imago*R13*(Yl33 + conjg(Yl33)))/(2.*vev)
    gh2e3e3(2) = (R22*(-Yl33 + conjg(Yl33)) - &
         imago*R23*(Yl33 + conjg(Yl33)))/(2.*vev)
    gh3e3e3(2) = (R32*(-Yl33 + conjg(Yl33)) - &
         imago*R33*(Yl33 + conjg(Yl33)))/(2.*vev)
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs
  end subroutine model_update_alpha_s
end module parameters_thdm_ckm
