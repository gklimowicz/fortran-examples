! parameters.Zprime.f90
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
module parameters_Zprime
  use kinds
  use constants
  implicit none
  private
  real(default), dimension(34), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, &
       ghgaga, ghgaz, ghgg, ghmm, & 		
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  complex(default), dimension(2), public :: gnchlep, gnchneu, gnchup, gnchdwn
  real(default) :: y_ql, y_ur, y_dr, y_ll, y_er
  real(default) :: bml_ql, bml_ur, bml_dr, bml_ll, bml_er
  real(default), public :: xi0 = 0, xipm = 0  
  real(default), public :: gyp, gbmlp
  real(default) :: glep_l, glep_r, gneu_l, gup_l, gup_r, gdwn_l, gdwn_r
  public :: import_from_whizard, model_update_alpha_s
contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(40), intent(in) :: par_array
    integer, intent(in) :: scheme
    type :: parameter_set
       real(default) :: GF
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
       real(default) :: mZH
       real(default) :: wZ
       real(default) :: wW
       real(default) :: wZH
       real(default) :: wH
       real(default) :: v_lep
       real(default) :: v_neu
       real(default) :: v_up
       real(default) :: v_dwn
       real(default) :: a_lep
       real(default) :: a_neu
       real(default) :: a_up
       real(default) :: a_dwn
       real(default) :: glepv
       real(default) :: gneuv
       real(default) :: gupv
       real(default) :: gdwnv
       real(default) :: glepa
       real(default) :: gneua
       real(default) :: gupa
       real(default) :: gdwna
       real(default) :: gyp
       real(default) :: gbmlp
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    real(default) :: e, g, qelep, qeup, qedwn, sin2thw, sinthw, costhw
    real(default) :: sinpsi, cospsi, atpsi
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
    par%mZH    = par_array(14)
    par%wZ     = par_array(15)
    par%wW     = par_array(16)
    par%wZH    = par_array(17)
    par%wH     = par_array(18)
    par%v_lep  = par_array(19)
    par%v_neu  = par_array(20)
    par%v_up   = par_array(21)
    par%v_dwn  = par_array(22)
    par%a_lep  = par_array(23)
    par%a_neu  = par_array(24)
    par%a_up   = par_array(25)
    par%a_dwn  = par_array(26)
    par%glepv  = par_array(27)
    par%gneuv  = par_array(28)
    par%gupv   = par_array(29)
    par%gdwnv  = par_array(30) 
    par%glepa  = par_array(31)
    par%gneua  = par_array(32)
    par%gupa   = par_array(33)
    par%gdwna  = par_array(34)
    par%gyp    = par_array(35)
    par%gbmlp  = par_array(36)
    par%v      = par_array(37)
    par%cw     = par_array(38)
    par%sw     = par_array(39)
    par%ee     = par_array(40)
    mass(1:34) = 0
    width(1:34) = 0
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
    mass(32) = par%mZH
    width(32) = par%wZH
    vev = 2 * par%mW * par%sw / par%ee 
    e = par%ee
    sinthw = par%sw
    sin2thw = sinthw**2
    costhw = par%cw
    qelep = - 1.0_default
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
    gyp = par%gyp
    gbmlp = par%gbmlp
    gcc = - g / 2 / sqrt (2.0_double)
    gncneu(1) = - g / 2 / costhw * ( + 0.5_double)
    gnclep(1) = - g / 2 / costhw * ( - 0.5_double - 2 * qelep * sin2thw)
    gncup(1)  = - g / 2 / costhw * ( + 0.5_double - 2 * qeup  * sin2thw)
    gncdwn(1) = - g / 2 / costhw * ( - 0.5_double - 2 * qedwn * sin2thw)
    gncneu(2) = - g / 2 / costhw * ( + 0.5_double)
    gnclep(2) = - g / 2 / costhw * ( - 0.5_double)
    gncup(2)  = - g / 2 / costhw * ( + 0.5_double)
    gncdwn(2) = - g / 2 / costhw * ( - 0.5_double)
    qlep = - e * qelep
    qup = - e * qeup
    qdwn = - e * qedwn
    qw = e
    iqw = (0,1)*qw
    gzww = g * costhw
    igzww = (0,1)*gzww
    gwww = g
    igwww = (0,1)*gwww
    ghww = mass(24) * g
    ghhww = g**2 / 2.0_default
    ghzz = mass(23) * g / costhw
    ghhzz = g**2 / 2.0_default / costhw**2
    ghtt = - mass(6) / vev
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    gh3 = - 3 * mass(25)**2 / vev
    gh4 = - 3 * mass(25)**2 / vev**2
    gw4 = gwww**2
    gzzww = gzww**2
    gazww = gzww * qw
    gaaww = qw**2
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs    
!!! Additional parameters to the Standard Model
!!! There are three independent possibilities to set the couplings of the heavy
!!! Z boson: multiples of the Z couplings and and explicit V and A coupling,
!!! also there is the option of a sequential anomaly-free Zprime.
    y_ql   =   1.0_default / 6.0_default
    y_ur   =   2.0_default / 3.0_default
    y_dr   = - 1.0_default / 3.0_default
    y_ll   = - 1.0_default / 2.0_default
    y_er   = - 1.0_default 
    bml_ql =   1.0_default / 3.0_default
    bml_ur =   1.0_default / 3.0_default
    bml_dr =   1.0_default / 3.0_default
    bml_ll = - 1.0_default
    bml_er = - 1.0_default     
    glep_l = gyp * y_ll + gbmlp * bml_ll 
    glep_r = gyp * y_er + gbmlp * bml_er 
    gneu_l = gyp * y_ll + gbmlp * bml_ll 
    gup_l  = gyp * y_ql + gbmlp * bml_ql 
    gup_r  = gyp * y_ur + gbmlp * bml_ur
    gdwn_l = gyp * y_ql + gbmlp * bml_ql
    gdwn_r = gyp * y_dr + gbmlp * bml_dr 
    gnchlep(1) = par%v_lep * gnclep(1) + par%glepv &
                 + (glep_r + glep_l) / 2.0_default
    gnchneu(1) = par%v_neu * gncneu(1) + par%gneuv &
                 + gneu_l / 2.0_default
    gnchup(1)  = par%v_up  * gncup(1)  + par%gupv &
                 + (gup_r + gup_l) / 2.0_default
    gnchdwn(1) = par%v_dwn * gncdwn(1) + par%gdwnv &
                 + (gdwn_r + gdwn_l) / 2.0_default
    gnchlep(2) = par%a_lep * gnclep(2) + par%glepa &
                 + (glep_r - glep_l) / 2.0_default                 
    gnchneu(2) = par%a_neu * gncneu(2) + par%gneua &
                 - gneu_l / 2.0_default                 
    gnchup(2)  = par%a_up  * gncup(2)  + par%gupa &
                 + (gup_r - gup_l) / 2.0_default                 
    gnchdwn(2) = par%a_dwn * gncdwn(2) + par%gdwna &
                 + (gdwn_r - gdwn_l) / 2.0_default                 
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs     
    !!! The Hgg coupling should not get a running alpha_s
  end subroutine model_update_alpha_s
end module parameters_Zprime
