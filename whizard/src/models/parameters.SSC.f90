! parameters.SSC.f90
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
module parameters_ssc
  use kinds
  use constants 
  implicit none
  private

  real(default), dimension(59), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs

  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, &
       ghgaga, ghgaz, ghgg, ghmm, & 		
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww, &
       gszz, gszzt, gsww, gswwt, gssww, gsszz, &
       gpnww, gpnzz, gpwz, gpww, &
       gpsnww, gpsnzz, gpsnhh, &
       gfww, gfzz, gfwwt, gfzzt, &
       gtnww, gtnzz, gtwz, gtww, &
       gtsnww, gtsnzz
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  real(default), public :: a4, a5, a6, a7, a10
  complex(default), public :: ig1a, ig1z, rg5a, rg5z, &
       ig1pkpg4a, ig1pkpg4z, ig1pkmg4a, ig1pkmg4z, &
       ig1mkpg4a, ig1mkpg4z, ig1mkmg4a, ig1mkmg4z, &
       ila, ilz, il5a, il5z, ik5a, ik5z, &
       ialww0, ialww2, ialzw0, ialzw1, ialzz, &
       alww0, alww2, alzw0, alzw1, alzz, &
       igdh4, gdh2w2, gdh2z2, gdhw2, gdhz2
  real(default), public :: lam_reg   
  real(default), public :: fudge_higgs, fudge_km, w_res, &
       amp00, amp02, amp11, amp20, amp22
  real(default), dimension(1:12), public :: gkm
  real(default), dimension(1:5), public :: mkm, wkm
  complex(default), public :: ghvva
  integer, public :: part_r
  logical, public :: unit_limit
  integer, private :: i

  public :: import_from_whizard, model_update_alpha_s

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(66), intent(in) :: par_array
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
       real(default) :: xi0
       real(default) :: xipm
       real(default) :: a4
       real(default) :: a5
       real(default) :: a6
       real(default) :: a7
       real(default) :: a10
       real(default) :: mkm_s
       real(default) :: mkm_p
       real(default) :: mkm_r
       real(default) :: mkm_f
       real(default) :: mkm_t
       real(default) :: gkm_s
       real(default) :: gkm_p
       real(default) :: gkm_r
       real(default) :: gkm_f
       real(default) :: gkm_t
       real(default) :: wkm_s
       real(default) :: wkm_p
       real(default) :: wkm_r
       real(default) :: wkm_f
       real(default) :: wkm_t
       real(default) :: g1a
       real(default) :: g1z
       real(default) :: g4a
       real(default) :: g4z
       real(default) :: g5a
       real(default) :: g5z
       real(default) :: ka
       real(default) :: kz
       real(default) :: la
       real(default) :: lz
       real(default) :: k5a
       real(default) :: k5z
       real(default) :: l5a
       real(default) :: l5z
       real(default) :: mreg
       real(default) :: fkm
       real(default) :: wres
       real(default) :: gkm_st
       real(default) :: gkm_pt
       real(default) :: gkm_rt
       real(default) :: gkm_ft
       real(default) :: gkm_tt
       real(default) :: fmixed
       real(default) :: fkappa
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
    par%xi0    = par_array(17)
    par%xipm   = par_array(18)
    par%a4     = par_array(19)         
    par%a5     = par_array(20)
    par%a6     = par_array(21)
    par%a7     = par_array(22)
    par%a10    = par_array(23)
    par%mkm_s  = par_array(24)
    par%mkm_p  = par_array(25)
    par%mkm_r  = par_array(26)
    par%mkm_f  = par_array(27)
    par%mkm_t  = par_array(28)
    par%gkm_s  = par_array(29)
    par%gkm_p  = par_array(30)
    par%gkm_r  = par_array(31)
    par%gkm_f  = par_array(32)
    par%gkm_t  = par_array(33)
    par%wkm_s  = par_array(34)
    par%wkm_p  = par_array(35)
    par%wkm_r  = par_array(36)
    par%wkm_f  = par_array(37)
    par%wkm_t  = par_array(38)
    par%g1a    = par_array(39)
    par%g1z    = par_array(40)
    par%g4a    = par_array(41)
    par%g4z    = par_array(42)
    par%g5a    = par_array(43)
    par%g5z    = par_array(44)
    par%ka     = par_array(45)
    par%kz     = par_array(46)
    par%la     = par_array(47)
    par%lz     = par_array(48)
    par%k5a    = par_array(49)
    par%k5z    = par_array(50)
    par%l5a    = par_array(51)
    par%l5z    = par_array(52)
    par%mreg   = par_array(53)
    par%fkm    = par_array(54)
    par%wres   = par_array(55)
    par%gkm_st = par_array(56)
    par%gkm_pt = par_array(57)
    par%gkm_rt = par_array(58)
    par%gkm_ft = par_array(59)
    par%gkm_tt = par_array(60)
    par%fmixed = par_array(61)
    par%fkappa = par_array(62)
    par%v      = par_array(63)
    par%cw     = par_array(64)
    par%sw     = par_array(65)
    par%ee     = par_array(66)
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
    mass(45) = par%mkm_s
    width(45) = par%wkm_s
    mass(46) = par%mkm_p
    width(46) = par%wkm_p
    mass(47) = par%mkm_p
    width(47) = par%wkm_p
    mass(48) = par%mkm_p
    width(48) = par%wkm_p
    mass(49) = par%mkm_p
    width(49) = par%wkm_p
    mass(52) = par%mkm_f
    width(52) = par%wkm_f
    mass(53) = par%mkm_t
    width(53) = par%wkm_t
    mass(54) = par%mkm_t
    width(54) = par%wkm_t
    mass(55) = par%mkm_t
    width(55) = par%wkm_t
    mass(59) = par%mkm_t
    width(59) = par%wkm_t
    mkm(1) = par%mkm_s
    mkm(2) = par%mkm_p
    mkm(3) = par%mkm_r
    mkm(4) = par%mkm_f
    mkm(5) = par%mkm_t
    gkm(1) = par%gkm_s
    gkm(2) = par%gkm_p
    gkm(3) = par%gkm_r
    gkm(4) = par%gkm_f
    gkm(5) = par%gkm_t
    gkm(6) = par%gkm_st
    gkm(7) = par%gkm_pt
    gkm(8) = par%gkm_rt
    gkm(9) = par%gkm_ft
    gkm(10) = par%gkm_tt
    gkm(11) = par%fmixed
    gkm(12) = par%fkappa
    wkm(1) = par%wkm_s
    wkm(2) = par%wkm_p
    wkm(3) = par%wkm_r
    wkm(4) = par%wkm_f
    wkm(5) = par%wkm_t
    e = par%ee
    sinthw = par%sw
    sin2thw = par%sw**2
    costhw = par%cw
    qelep = - 1
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
    w_res = par%wres
    vev = par%v
    do i=1,5
      if (w_res == 1 .and. wkm(i) == 0) then
        select case (i)
          case (1) !!! Scalar isosinglet
            wkm(1) = 3.*gkm(1)**2/32./Pi * mkm(1)**3/vev**2
            width(45) = wkm(1)
!!            write (*, "(1x,A,ES19.12)")  "Setting width: wkm_s =", wkm(1)
          case (2) !!! Scalar isoquintet
            wkm(2) = gkm(2)**2/64./Pi * mkm(2)**3/vev**2
            width(46) = wkm(2)
            width(47) = wkm(2)
            width(48) = wkm(2)
!!            write (*, "(1x,A,ES19.12)")  "Setting width: wkm_p =", wkm(2)
          case (3) !!! Vector isotriplet
            wkm(3) = gkm(3)**2/48./Pi * mkm(3)
!!            write (*, "(1x,A,ES19.12)")  "Setting width: wkm_r =", wkm(3)
          case (4) !!! Tensor isosinglet
            wkm(4) = gkm(4)**2/320./Pi * mkm(4)**3/vev**2
            width(52) = wkm(4)
 !!           write (*, "(1x,A,ES19.12)")  "Setting width: wkm_f =", wkm(4)
          case (5) !!! Tensor isoquintet
            wkm(5) = gkm(5)**2/1920./Pi * mkm(5)**3/vev**2
            width(53) = wkm(5)
            width(54) = wkm(5)
            width(55) = wkm(5)
!!            write (*, "(1x,A,ES19.12)")  "Setting width: wkm_t =", wkm(5)
        end select
      end if
    end do
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
    gsww = gkm(1) * mass(24) * g
    gszz = gkm(1) * mass(23) * g / costhw
    gswwt = gkm(6) * g**3 / mass(24) / (16.0 * PI) 
    gszzt = gkm(6) * g**3 / costhw**3 / mass(23) /(16.0 * PI)
    gpnww = - gkm(2) * mass(24) * g / 2 / sqrt(3.0_default)
    gpnzz = gkm(2) * mass(23) * g / costhw  / sqrt(3.0_default)
    gpsnww = 0
    gpsnzz = 0
    gpsnhh = 0
    gpwz = gkm(2) * mass(23) * g / 2
    gpww = gkm(2) * mass(24) * g / sqrt(2.0_default)
    gfww = gkm(4) * mass(24) * g / 2
    gfzz = gkm(4) * mass(23) * g / costhw / 2
    gfwwt = gkm(9) * g**3 / mass(24) / (32.0 * PI) 
    gfzzt = gkm(9) * g**3 / costhw**3 / mass(23) /(32.0 * PI)
    gtnww = - gkm(5) * mass(24) * g / 4 / sqrt(3.0_default)
    gtnzz = gkm(5) * mass(23) * g / costhw / 2 / sqrt(3.0_default)
    gtsnww = 0
    gtsnzz = 0
    gtwz = gkm(5) * mass(23) * g / 4
    gtww = gkm(5) * mass(24) * g / 2 / sqrt(2.0_default)
    gssww = 0
    gsszz = 0
    part_r = 1
    unit_limit = .false.
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs    
    a4 = par%a4
    a5 = par%a5
    a6 = par%a6
    a7 = par%a7
    a10 = par%a10
    lam_reg = par%mreg
    fudge_higgs = 0
    ghvva = 0
    fudge_km = par%fkm
    ig1a = iqw * par%g1a
    ig1z = igzww * par%g1z
    ig1pkpg4a = iqw   * (par%g1a + par%ka + par%g4a) / 2
    ig1pkpg4z = igzww * (par%g1z + par%kz + par%g4z) / 2
    ig1pkmg4a = iqw   * (par%g1a + par%ka - par%g4a) / 2
    ig1pkmg4z = igzww * (par%g1z + par%kz - par%g4z) / 2
    ig1mkpg4a = iqw   * (par%g1a - par%ka + par%g4a) / 2
    ig1mkpg4z = igzww * (par%g1z - par%kz + par%g4z) / 2
    ig1mkmg4a = iqw   * (par%g1a - par%ka - par%g4a) / 2
    ig1mkmg4z = igzww * (par%g1z - par%kz - par%g4z) / 2
    ila = iqw   * par%la / (mass(24)*mass(24))
    ilz = igzww * par%lz / (mass(24)*mass(24))
    rg5a = qw   * par%g5a
    rg5z = gzww * par%g5z
    ik5a = iqw   * par%k5a
    ik5z = igzww * par%k5z
    il5a = iqw   * par%l5a / (mass(24)*mass(24))
    il5z = igzww * par%l5z / (mass(24)*mass(24))
    alww0 = g**4 * (a4 + 2 * a5)
    alww2 = g**4 * 2 * a4
    alzw1 = g**4 / costhw**2 * (a4 + a6)
    alzw0 = g**4 / costhw**2 * 2 * (a5 + a7)
    alzz = g**4 / costhw**4 * 2 * (a4 + a5 + (a6+a7+a10)*2)
    ialww0 = g**2 * sqrt (-cmplx(a4 + 2 * a5, kind=default))
    ialww2 = g**2 * sqrt (-cmplx(2 * a4, kind=default))
    ialzw1 = g**2 / costhw * sqrt (-cmplx(a4 + a6, kind=default))
    ialzw0 = g**2 / costhw &
         & * sqrt (-cmplx(2 * (a5 + a7), kind=default))
    ialzz  = g**2 / (costhw*costhw) &
         & * sqrt (-cmplx(2 * (a4 + a5 + (a6+a7+a10)*2), &
         &                kind=default))
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs     
    !!! The Hgg coupling should not get a running alpha_s
  end subroutine model_update_alpha_s
end module parameters_ssc
