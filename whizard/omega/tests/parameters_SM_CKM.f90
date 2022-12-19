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

module parameters_sm_ckm
  use kinds
  use constants
  implicit none
  private

  real(default), dimension(27), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs

  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, ghmm, &
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  complex(default), public :: &
       gccq11 = 0, gccq12 = 0, gccq13 = 0, gccq21 = 0, &
       gccq22 = 0, gccq23 = 0, gccq31 = 0, gccq32 = 0, gccq33 = 0
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn

  public :: import_from_whizard, model_update_alpha_s, init_parameters

contains

  subroutine init_parameters
    real(default), dimension(34) :: vars
    vars(1)  = 1.16639E-5   ! Fermi constant
    vars(2)  = 91.1882      ! Z-boson mass  
    vars(3)  = 80.419       ! W-boson mass  
    vars(4)  = 125          ! Higgs mass
    vars(5)  = 0.1178       ! Strong coupling constant (Z point)             
    vars(6)  = 0.000510997  ! electron mass                       
    vars(7)  = 0.105658389  ! muon mass                         
    vars(8)  = 1.77705      ! tau-lepton mass                               
    vars(9)  = 0.095        ! s-quark mass                         
    vars(10) = 1.2          ! c-quark mass                         
    vars(11) = 4.2          ! b-quark mass                         
    vars(12) = 171.9        ! t-quark mass                         
    vars(13) = 1.523        ! t-quark width
    vars(14) = 2.443        ! Z-boson width
    vars(15) = 2.049        ! W-boson width
    vars(16) = 0.004143     ! Higgs width  
    vars(17) = 0.97383      ! Vud
    vars(18) = 0.2272       ! Vus
    vars(19) = 0.00396      ! Vub
    vars(20) = -0.2271      ! Vcd
    vars(21) = 0.97296      ! Vcs
    vars(22) = 0.04221      ! Vcb
    vars(23) = 0.00814      ! Vtd
    vars(24) = -0.04161     ! Vts
    vars(25) = 0.99910      ! Vtb
    vars(26) = 0.100        ! anomaly Higgs couplings K factors
    vars(27) = 0.200        ! anomaly Higgs couplings K factors
    vars(28) = 0.300        ! anomaly Higgs couplings K factors
    vars(29) = 0.400        ! R_xi parameter for Z-boson
    vars(30) = 0.500        ! R_xi parameter for W-boson
    vars(31) = 1 / sqrt (sqrt (2.) * vars(1))    ! v (Higgs vev)
    vars(32) = vars(3) / vars(2)                 ! cos(theta-W)
    vars(33) = sqrt (1-vars(32)**2)              ! sin(theta-W)
    vars(34) = 2 * vars(33) * vars(3) / vars(31) ! em-coupling (GF scheme)
    call import_from_whizard (vars)
  end subroutine init_parameters

  subroutine import_from_whizard (par_array)
    real(default), dimension(34), intent(in) :: par_array
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
       real(default) :: vckm11
       real(default) :: vckm12
       real(default) :: vckm13
       real(default) :: vckm21
       real(default) :: vckm22
       real(default) :: vckm23
       real(default) :: vckm31
       real(default) :: vckm32
       real(default) :: vckm33
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
    par%vckm11 = par_array(17)
    par%vckm12 = par_array(18)
    par%vckm13 = par_array(19)
    par%vckm21 = par_array(20)
    par%vckm22 = par_array(21)
    par%vckm23 = par_array(22)
    par%vckm31 = par_array(23)
    par%vckm32 = par_array(24)
    par%vckm33 = par_array(25)
    par%khgaz  = par_array(26)
    par%khgaga = par_array(27)
    par%khgg   = par_array(28)
    par%xi0    = par_array(29)
    par%xipm   = par_array(30)
    par%v      = par_array(31)
    par%cw     = par_array(32)
    par%sw     = par_array(33)
    par%ee     = par_array(34)
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
    g = e / sinthw
    gcc = - g / 2 / sqrt (2.0_default)
    gccq11 = gcc * par%vckm11
    gccq12 = gcc * par%vckm12
    gccq13 = gcc * par%vckm13
    gccq21 = gcc * par%vckm21
    gccq22 = gcc * par%vckm22
    gccq23 = gcc * par%vckm23
    gccq31 = gcc * par%vckm31
    gccq32 = gcc * par%vckm32
    gccq33 = gcc * par%vckm33
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
end module parameters_sm_ckm
