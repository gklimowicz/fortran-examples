! parameters.SM_top.f90
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
module parameters_sm_top
  use kinds
  use constants
  use sm_physics !NODEP!
  implicit none
  private

  real(default), dimension(27), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs
  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn, qetop
  complex(default), public :: qlep, qup, qdwn, qtop, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, ghmm, &             
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gnclep, gncneu, gncup, gncdwn, gnctop

  public :: import_from_whizard, model_update_alpha_s

contains

  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(25), intent(in) :: par_array
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
    par%v      = par_array(22)
    par%cw     = par_array(23)
    par%sw     = par_array(24)
    par%ee     = par_array(25)
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
    vev = par%v
    e = par%ee
    sinthw = par%sw
    sin2thw = par%sw**2
    costhw = par%cw
    tanthw = sinthw/costhw
    qelep = - 1
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    qetop = - 4.0_default / 3.0_default
    g = e / sinthw
    gcc = - g / 2 / sqrt (2.0_double)
    gncneu(1) = - g / 2 / costhw * ( + 0.5_double)
    gnclep(1) = - g / 2 / costhw * ( - 0.5_double - 2 * qelep * sin2thw)
    gncup(1)  = - g / 2 / costhw * ( + 0.5_double - 2 * qeup  * sin2thw)
    gnctop(1) = - g / 2 / costhw * ( + 0.5_double - 2 * qetop  * sin2thw)
    gncdwn(1) = - g / 2 / costhw * ( - 0.5_double - 2 * qedwn * sin2thw)
    gncneu(2) = - g / 2 / costhw * ( + 0.5_double)
    gnclep(2) = - g / 2 / costhw * ( - 0.5_double)
    gncup(2)  = - g / 2 / costhw * ( + 0.5_double)
    gnctop(2) = - g / 2 / costhw * ( + 0.5_double)
    gncdwn(2) = - g / 2 / costhw * ( - 0.5_double)
    qlep = - e * qelep
    qup = - e * qeup
    qdwn = - e * qedwn
    qtop = - e * qetop
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
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs     
    !!! The Hgg coupling should not get a running alpha_s
  end subroutine model_update_alpha_s
end module parameters_sm_top
