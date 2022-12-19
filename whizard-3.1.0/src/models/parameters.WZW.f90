! parameters.WZW.f90
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
module parameters_wzw
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
  complex(default), public :: getagg, gpsigg, getaww, gpsiww, getazz, gpsizz
  complex(default), public :: getaaa, gpsiaa, getaaz, gpsiaz
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  complex(default), dimension(2), public :: gnchlep, gnchneu, gnchup, gnchdwn
  real(default) :: y_ql, y_ur, y_dr, y_ll, y_er
  real(default) :: bml_ql, bml_ur, bml_dr, bml_ll, bml_er
  real(default), public :: xi0 = 0, xipm = 0  
  real(default) :: glep_l, glep_r, gneu_l, gup_l, gup_r, gdwn_l, gdwn_r
  public :: import_from_whizard, model_update_alpha_s
contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(33), intent(in) :: par_array
    integer, intent(in) :: scheme
    type :: parameter_set
       real(default) :: GF
       real(default) :: mZ
       real(default) :: mW
       real(default) :: mH
       real(default) :: mPsi
       real(default) :: mEta
       real(default) :: alphas
       real(default) :: fpsi
       real(default) :: feta
       real(default) :: kpsig
       real(default) :: kpsiw
       real(default) :: kpsib
       real(default) :: ketag
       real(default) :: ketaw
       real(default) :: ketab
       real(default) :: lambda
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
       real(default) :: wPsi
       real(default) :: wEta
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
    par%mPsi   = par_array(5)
    par%mEta   = par_array(6)
    par%alphas = par_array(7)
    par%fPsi   = par_array(8)
    par%fEta   = par_array(9)
    par%kpsig  = par_array(10)
    par%kpsiw  = par_array(11)
    par%kpsib  = par_array(12)
    par%ketag  = par_array(13)
    par%ketaw  = par_array(14)
    par%ketab  = par_array(15)
    par%lambda = par_array(16)
    par%me     = par_array(17)
    par%mmu    = par_array(18)
    par%mtau   = par_array(19)    
    par%ms     = par_array(20)
    par%mc     = par_array(21)
    par%mb     = par_array(22)
    par%mtop   = par_array(23)
    par%wtop   = par_array(24)
    par%wZ     = par_array(25)
    par%wW     = par_array(26)
    par%wH     = par_array(27)
    par%wPsi   = par_array(28)
    par%wEta   = par_array(29)
    par%v      = par_array(30)
    par%cw     = par_array(31)
    par%sw     = par_array(32)
    par%ee     = par_array(33)
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
    mass(28) = par%mPsi
    width(28) = par%wPsi
    mass(29) = par%mEta
    width(29) = par%wEta    
    vev = 2 * par%mW * par%sw / par%ee 
    e = par%ee
    sinthw = par%sw
    sin2thw = sinthw**2
    costhw = par%cw
    qelep = - 1.0_default
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
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
    !!!
    getagg = par%ketag * par%alphas / PI / par%feta
    getaww = par%ketaw * g**2 / 4 / PI**2 / par%feta
    getazz = g**2 * (par%ketaw * costhw**2 + par%ketab * sinthw**4/costhw**2) / 4 / PI**2 / par%feta
    getaaa = e**2 * (par%ketaw + par%ketab) / 4 / PI**2 / par%feta
    getaaz = e * g * (par%ketaw * costhw - par%ketab * sinthw**2/costhw) / 4 / PI**2 / par%feta
    gpsigg = par%kpsig * par%alphas / PI / par%fpsi
    gpsiww = par%kpsiw * g**2 / 4 / PI**2 / par%fpsi
    gpsizz = g**2 * (par%kpsiw * costhw**2 + par%kpsib * sinthw**4/costhw**2) / 4 / PI**2 / par%fpsi
    gpsiaa = e**2 * (par%kpsiw + par%kpsib) / 4 / PI**2 / par%fpsi
    gpsiaz = e * g * (par%kpsiw * costhw - par%kpsib * sinthw**2/costhw) / 4 / PI**2 / par%fpsi    
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs     
    !!! The Hgg coupling should not get a running alpha_s
  end subroutine model_update_alpha_s
end module parameters_wzw
