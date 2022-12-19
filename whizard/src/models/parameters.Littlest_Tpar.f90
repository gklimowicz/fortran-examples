! parameters.LittlestTpar.f90
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
module parameters_littlest_tpar
  use kinds
  use constants
  use sm_physics !NODEP!
  implicit none
  private
  complex(default), public :: &
       gh0ww, gh0zz, &
       gh0tt, gh0bb, gh0cc, gh0tautau, gh0mumu, &
       iga0tt, iga0bb, iga0cc, iga0tautau, iga0mumu, &
       gahh, gzhh, igzha, igzh0a 
  complex(default), dimension(2), public :: &
       ghptb, ghpcs, ghptaunu, ghpmunu 
  !!! Additional Littlest Higgs parameters 
  complex(default), public :: &
       ghwhw, ghwhwh, ghahah, ghzhz, ghzhah, &
       ghahz, ghaa, geaa, geaz, gegg 
  complex(default), public :: & 
       gpsiww, gpsiwhw, gpsizz, gpsizhzh, & 
       gpsizhz, gpsizah, gpsizhah, gpsiahah, & 
       gpsizw, gpsizwh, gpsiahw, gpsiahwh, &
       gpsizhw, gpsizhwh, gpsippww, gpsippwhw, &
       gpsippwhwh, gpsihw, gpsihwh, gpsi0w, & 
       gpsi0wh, gpsi1w, gpsi1wh, gpsippw, &
       gpsippwh 
  complex(default), public :: & 
       gpsihah, gpsi0ah, gahpsip, &
       gpsi1hz, gpsi1hzh, gpsi01z, gpsi01zh, &
       gzpsip, gzpsipp, gzhpsipp 
  complex(default), public :: &
       ghhaa, ghhwhw, ghhzhz, ghhahz, ghhzhah, &
       ghpsi0ww, ghpsi0whw, ghpsi0zz, ghpsi0zhzh, & 
       ghpsi0zhz, ghpsi0ahah, ghpsi0zah, ghpsi0zhah 
  complex(default), public :: &
      ghpsipwa, ghpsipwha, ghpsipwz, ghpsiwhz, &
      ghpsipwah, ghpsipwhah, ghpsipwzh, ghpsipwhzh, & 
      ghpsippww, ghpsippwhwh, ghpsippwhw, gpsi00zh, &
      gpsi00ah, gpsi00zhah, gpsi0pwa, gpsi0pwha, & 
      gpsi0pwz, gpsi0pwhz, gpsi0pwah, gpsi0pwhah, & 
      gpsi0pwzh, gpsi0pwhzh, gpsi0ppww, gpsi0ppwhwh, &
      gpsi0ppwhw, i_gpsi0pwa, i_gpsi0pwha, & 
      i_gpsi0pwz, i_gpsi0pwhz, i_gpsi0pwah, i_gpsi0pwhah, & 
      i_gpsi0pwzh, i_gpsi0pwhzh, i_gpsi0ppww, i_gpsi0ppwhwh, &
      i_gpsi0ppwhw 
  complex(default), public :: & 
      gpsippzz, gpsippzhzh, gpsippaz, gpsippaah, & 
      gpsippzah, gpsippwa, gpsippwha, gpsippwz, &
      gpsippwhz, gpsippwah, gpsippwhah, gpsippwzh, &
      gpsippwhzh, gpsicczz, gpsiccaz, gpsiccaah, &
      gpsicczzh, gpsiccazh, gpsicczah 
  complex(default), public :: &
      igahww, igzhww, igzwhw, igahwhwh, igzhwhwh, &
      igahwhw 
  complex(default), public :: &
      gwh4, gwhwhww, gwhwww, gwh3w, gwwaah, &
      gwwazh, gwwzzh, gwwzah, gwhwhaah, gwhwhazh, &
      gwhwhzzh, gwhwhzah, gwwzhah, gwhwhzhah, &
      gwhwzz, gwhwaz, gwhwaah, gwhwzah, gwhwzhzh, &
      gwhwzhah, gwhwazh, gwhwzzh 
  complex(default), public :: &
      qzup, gcch, gcctop, gccw, gccwh, &
      gnch, gztht, gzhtht, gah 
  complex(default), dimension(2), public :: &
      gnchup, gnchdwn, gnchneu, gnchlep, gahtt, &
      gahthth, ghtht, gpsipq2, gpsipq3, &
      ghhtht 
  complex(default), public :: &
      gahtht, ghthth, &
      gpsi0tt, gpsi0bb, gpsi0cc, gpsi0tautau, &
      gpsipl3, gpsi0tth, gpsi1tth, gpsipbth, &
      ghhtt, ghhthth 
  integer, parameter, public :: &
       n0 = 5, nloop = 2 
    real(default), dimension(38), public :: mass, width
    real(default) :: as
    complex(default), public :: gs, igs
    real(default), public :: e, g, gp, sinthw, costhw, sin2thw, tanthw, e_em
    real(default), public :: qelep, qeup, qedwn, vev, vevp
    real(default), public :: ttop, tbot, tch, ttau, tw
    real(default), public :: ltop, lbot, lc, ltau, lw
    real(default), public :: sint, sintp, sin2t, sin2tp, &
         cost, costp, cos2t, cos2tp
    real(default), public :: spsip, cpsip, spsi1, cpsi1, spsi0, cpsi0
    real(default), public :: t_fac, tp_fac, ttp_fac, c4s4
    real(default), public :: xzbp, xzwp, xh, xlam
    real(default), public :: f_vev, lam1, lam2
    real(default), public :: ye, yu
    complex(default), public :: qlep, qup, qdwn, gcc, qw, &
         gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
         ghbb, ghtt, ghcc, ghtautau, gh3, gh4, &
         ghgaga, ghgaz, ghgg, ghmm, & 		
         iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
    complex(default), dimension(2), public :: &
         gncneu, gnclep, gncup, gncdwn
    real(default), parameter :: xi0 = 0.0, xipm = 0.0 
    real(default), parameter :: &
         acc = 1.e-12_default
    real(default), parameter :: &
         asmz = .118_default

  public :: import_from_whizard, model_update_alpha_s

contains

  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(46), intent(in) :: par_array
    integer, intent(in) :: scheme
    !!! This corresponds to 1/alpha = 137.03598949333
    real(default), parameter :: &
         alpha = 1.0_default/137.03598949333_default
    real(default), parameter :: &
         asmz = 0.118_default
    real(default), parameter :: & 
         one = 1.0_default, two = 2.0_default, three = 3.0_default, &
         four = 4.0_default, five = 5.0_default
    type :: parameter_set
       real(default) :: gf
       real(default) :: mZ
       real(default) :: mW
       real(default) :: mH
       real(default) :: meta
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
       real(default) :: weta
       real(default) :: khgaz
       real(default) :: khgaga
       real(default) :: khgg
       real(default) :: mtoph
       real(default) :: mah
       real(default) :: mzh
       real(default) :: mwh
       real(default) :: mpsi0
       real(default) :: mpsi1
       real(default) :: mpsip
       real(default) :: mpsipp
       real(default) :: wtoph
       real(default) :: wah
       real(default) :: wzh
       real(default) :: wwh
       real(default) :: wpsi0
       real(default) :: wpsi1
       real(default) :: wpsip
       real(default) :: wpsipp
       real(default) :: st
       real(default) :: stp
       real(default) :: vp
       real(default) :: f
       real(default) :: lam1      
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
       par%gf     = par_array(1)
       par%mZ     = par_array(2)
       par%mW     = par_array(3)
       par%mH     = par_array(4)
       par%meta   = par_array(5)
       par%alphas = par_array(6)
       par%me     = par_array(7)
       par%mmu    = par_array(8)
       par%mtau   = par_array(9)
       par%ms     = par_array(10)
       par%mc     = par_array(11)
       par%mb     = par_array(12)
       par%mtop   = par_array(13)
       par%wtop   = par_array(14)
       par%wZ     = par_array(15)
       par%wW     = par_array(16)
       par%wH     = par_array(17)
       par%weta   = par_array(18)
       par%khgaz  = par_array(19)
       par%khgaga = par_array(20)
       par%khgg   = par_array(21)
       par%mtoph  = par_array(22)
       par%mah    = par_array(23)
       par%mzh    = par_array(24)
       par%mwh    = par_array(25)
       par%mpsi0  = par_array(26)
       par%mpsi1  = par_array(27)
       par%mpsip  = par_array(28)
       par%mpsipp = par_array(29)
       par%wtoph  = par_array(30)
       par%wah    = par_array(31)
       par%wzh    = par_array(32)
       par%wwh    = par_array(33)
       par%wpsi0  = par_array(34)
       par%wpsi1  = par_array(35)
       par%wpsip  = par_array(36)
       par%wpsipp = par_array(37)
       par%st     = par_array(38)
       par%stp    = par_array(39)
       par%vp     = par_array(40)
       par%f      = par_array(41)
       par%lam1   = par_array(42)
       par%v      = par_array(43)
       par%cw     = par_array(44)
       par%sw     = par_array(45) 
       par%ee     = par_array(46)
    e_em = sqrt(four * PI * alpha)
    mass(1:38) = 0
    width(1:38) = 0
    mass(3) = par%ms
    mass(4) = par%mc
    mass(5) = par%mb
    mass(6) = par%mtop
    width(6) = par%wtop
    mass(8) = par%mtoph
    width(8) = par%wtoph
    mass(11) = par%me
    mass(13) = par%mmu
    mass(15) = par%mtau
    mass(23) = par%mZ
    width(23) = par%wZ
    mass(24) = par%mW
    width(24) = par%wW
    mass(25) = par%mH
    width(25) = par%wH
    mass(26) =  xi0 * mass(23)
    width(26) =  0
    mass(27) =  xipm * mass(24)
    width(27) =  0
    mass(32) = par%mAH
    width(32) = par%wAH
    mass(33) = par%mZH
    width(33) = par%wZH
    mass(34) = par%mWH
    width(34) = par%wWH
    mass(35) = par%mpsi0
    width(35) = par%wpsi0
    mass(36) = par%mpsi1
    width(36) = par%wpsi1
    mass(37) = par%mpsip
    width(37) = par%wpsip
    mass(38) = par%mpsipp
    width(38) = par%wpsipp
    !!! This choice is responsible for anomaly cancellation:
    yu = - two / five
    ye = three / five
    ttop = four * mass(6)**2 / mass(25)**2
    tbot = four * mass(5)**2 / mass(25)**2
    tch  = four * mass(4)**2 / mass(25)**2
    ttau = four * mass(15)**2 / mass(25)**2
    tw   = four * mass(24)**2 / mass(25)**2  
    ltop = four * mass(6)**2 / mass(23)**2
    lbot = four * mass(5)**2 / mass(23)**2  
    lc   = four * mass(4)**2 / mass(23)**2
    ltau = four * mass(15)**2 / mass(23)**2
    lw   = four * mass(24)**2 / mass(23)**2
    vev = 2 * mass(24) * par%sw / par%ee 
    vevp = par%vp
    f_vev = par%f
    e = par%ee
    sinthw = par%sw
    sin2thw = sinthw**2
    costhw = par%cw
    tanthw = sinthw/costhw
    qelep = - one
    qeup = two / three
    qedwn = - one / three
    g = e / sinthw
    gp = e / costhw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Littlest Higgs special values
    sint = par%st
    sintp = par%stp
    sin2t = sint**2
    sin2tp = sintp**2
    cos2t = 1.0_default - sin2t
    cos2tp =  1.0_default - sin2tp
    cost = sqrt(cos2t)
    costp = sqrt(cos2tp)
    if (abs(sint) < 1.e-5 .or. abs(sint) > 1.-1.e-5 .or. &
         abs(sintp) < 1.e-5 .or. abs(sintp) > 1.-1.e-5) then
       write (6, "(A)") "This mixing is unphysical and leads to singularities"
       stop
    end if
    !!!
    lam1 = par%lam1
    lam2 = mass(8)/f_vev * (mass(6)/vev) / lam1
    xlam = lam1**2 / (lam1**2 + lam2**2)
    !!!
    xh = five*g*gp * sint*cost*sintp*costp * (cos2t*sin2tp + &
         sin2t*cos2tp) / (five*g**2*sin2tp*cos2tp - &
         gp**2*sin2t*cos2t) / two
    xzbp = - five/two*sintp*costp*(cos2tp-sin2tp)/sinthw
    xzwp = -sint*cost*(cos2t-sin2t)/two/costhw
    !!!
    spsi1 = sqrt(8.0_default)*vevp/sqrt(vev**2+8.0_default*vevp**2)
    cpsi1 = vev/sqrt(vev**2+8.0_default*vevp**2)
    spsip = two*vevp/sqrt(vev**2+four*vevp**2)
    cpsip = vev/sqrt(vev**2+four*vevp**2)
    spsi0 = sqrt(8.0_default)*vevp/vev
    !!!
    t_fac = (cos2t - sin2t) / two / sint / cost
    tp_fac = (cos2tp - sin2tp) / two / sintp / costp
    ttp_fac = (cos2t * sin2tp + sin2t * cos2tp) / sint / sintp &
         / cost / costp
    c4s4 = (cos2t**2 + sin2t**2) / two / cos2t / sin2t
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! In principle, we should allow here for deviations from 
    !!! the SM values 
    gcc = - g / 2 / sqrt (2.0_double)
    gncneu(1) = - g / two / costhw * ( + 0.5_double)
    gnclep(1) = - g / two / costhw * ( - 0.5_double - 2 * qelep * sin2thw)
    gncup(1)  = - g / two / costhw * ( + 0.5_double - 2 * qeup  * sin2thw)
    gncdwn(1) = - g / two / costhw * ( - 0.5_double - 2 * qedwn * sin2thw)
    gncneu(2) = - g / two / costhw * ( + 0.5_double)
    gnclep(2) = - g / two / costhw * ( - 0.5_double)
    gncup(2)  = - g / two / costhw * ( + 0.5_double)
    gncdwn(2) = - g / two / costhw * ( - 0.5_double)
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
    ghww = mass(24) * g * (one - vev**2/three/f_vev**2 + vev**2/f_vev**2 * &
         (cos2t - sin2t)**2/two - spsi0**2/two - two * sqrt(two) * &
         spsi0 * vevp/vev)
    !!! SM value !!! ghww = mass(24) * g    
    ghhww = g**2 / two
    ghzz = mass(23) * g / costhw * (one - vev**2/three/f_vev**2 - spsi0**2/two &
         + four*sqrt(two)*vevp/vev  - vev**2/f_vev**2/two * ((cos2t - sin2t)**2 &
         + five * (cos2tp - sin2tp)**2))
    !!! SM value !!! ghzz = mass(23) * g / costhw
    ghhzz = g**2 / two / costhw**2
    !!! ghtt = - mass(6) / vev !!! Will be corrected below
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    gh3 = - 3 * mass(25)**2 / vev
    gh4 = - 3 * mass(25)**2 / vev**2
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default * PI * par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Higgs anomaly couplings
    !!! SM LO loop factor (top,bottom,W)
    ghgaga = (-1._default) * alpha / vev / 2.0_default / PI * &
         (( 4.0_default * (fonehalf(ttop) + fonehalf(tch)) &
         + fonehalf(tbot)) / 3.0_default + fonehalf(ttau) + fone(tw)) &
         * sqrt(par%khgaga)
    !!! asymptotic limit:
    !!! ghgaga = (par%ee)**2 / vev / &
    !!!      9.0_default / pi**2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SM LO loop factor (only top and W)
    ghgaz = e * e_em / 8.0_default / PI**2 / vev *  &
          ( - 2.0_default + &
          16.0_default/3.0_default * sin2thw) * &
          (tri_i1(ttop,ltop) - tri_i2(ttop,ltop)) / costhw & 
          + ( - 1.0_default + &
          4.0_default/3.0_default * sin2thw) & 
          * (tri_i1(tbot,lbot) - tri_i2(tbot,lbot)) / costhw &
          + (-1.0_default + 4.0_default * sin2thw) &
          * (tri_i1(ttau,ltau) - tri_i2(ttau,ltau)) / costhw &
           - costhw * ( 4.0_default * (3.0_default - tanthw**2) * &
           tri_i2(tw,lw) + ((1 + 2.0_default/tw) * tanthw**2 - ( &
           5.0_default + 2.0_default/tw)) * tri_i1(tw,lw)) &
          /sinthw * sqrt(par%khgaz)
    !!! SM LO order loop factor with 
    !!! N(N)LO K factor = 2.1 (only top)
    !!! Limit of infinite top quark mass:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! We use par%gg because of sqrt(2) above
    ghgg = (-1._double) * par%alphas / vev / 4.0_default / PI * &
         (fonehalf(ttop) + fonehalf(tbot) + fonehalf(tch)) * &
         sqrt(par%khgg)
    !!! ghgg   = par%alphas / 3.0_default &
    !!!      / vev / pi * 2.1_default
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ghwhwh = - g**2 * vev / two
    ghahah = - gp**2 * vev / two
    ghwhw = t_fac * ghwhwh
    ghzhz = ghwhw / costhw
    ghzhah = - g * gp * vev * ttp_fac / four
    ghahz = - g * gp * tp_fac * vev / costhw / two
    !!!
    gpsiww = - g**2 * (spsi0*vev - sqrt(8.0_default)*vevp)/two
    gpsiwhw = - gpsiww * t_fac
    gpsizz = - g**2 * (spsi0*vev - four*sqrt(two)*vevp) / two / costhw**2
    gpsizhzh = g**2 * (spsi0*vev + sqrt(two)*vevp*t_fac) / two
    gpsizhz = - gpsizz * t_fac * costhw
    gpsizah = - gpsizz * tp_fac * costhw * gp / g
    gpsizhah = g*gp/four/sint/cost/sintp/costp * (vev*spsi0 * &
         (cos2t*sin2tp + sin2t*cos2tp) + sqrt(8.0_default) * vevp * &
         (cos2t - sin2t) * (cos2tp - sin2tp))  
    gpsiahah = gp**2 * (vev*spsi0 + tp_fac * sqrt(two) * vevp) / two
    !!!
    gpsizw = - g**2 / costhw * vevp 
    gpsizwh = - gpsizw * t_fac
    gpsiahw = - g*gp * tp_fac * (vev*spsip - four*vevp) / two
    gpsiahwh = - g*gp * vevp * (cos2t*cos2tp + sin2t*sin2tp) / two / &
         sint / cost / sintp / costp 
    gpsizhw = g**2 * vevp * t_fac
    gpsizhwh = - g**2 * vevp * c4s4
    !!!
    gpsippww = two * g**2 * vevp
    gpsippwhw = - gpsippww * t_fac
    gpsippwhwh = gpsippww * c4s4
    !!!
    gpsihw = - g/two * (sqrt(two) * spsi0 - spsip)
    gpsihwh = - gpsihw * t_fac
    gpsi0w = - g/sqrt(two)
    gpsi0wh = - gpsi0w * t_fac
    gpsi1w = imago * gpsi0w
    gpsi1wh = imago * gpsi0wh
    gpsippw = - g
    gpsippwh = g * t_fac
    !!!
    gpsihah = imago / two * gp * tp_fac * (spsi1 - two * spsi0)
    gpsi0ah = - imago * gp * tp_fac
    gahpsip = gp * tp_fac
    gpsi1hz = - imago * g / 2 / costhw * (spsi1 - two * spsi0)
    gpsi1hzh = imago * g / two * t_fac * (spsi1 - two * spsi0)
    gpsi01z = imago * g / costhw
    gpsi01zh = - imago * g * t_fac
    gzpsip = g / costhw * sin2thw
    gzpsipp = - g / costhw * (one - two * sin2thw)
    gzhpsipp = g * t_fac
    !!!
    ghhaa = - gp**2 / two
    ghhwhw = - g**2 / two * t_fac 
    ghhzhz = - g**2 / two / costhw * t_fac 
    ghhahz = - g*gp / 2 / costhw * tp_fac 
    ghhzhah = - g*gp / four * ttp_fac
    !!!
    ghpsi0ww = g**2 / two * spsi0
    ghpsi0whw =  - g**2 / two * spsi0 * t_fac
    ghpsi0zz = three * g**2 * spsi0 / two / costhw**2
    ghpsi0zhzh = g**2 / two * (one + t_fac**2) * spsi0 
    ghpsi0zhz =  - three * g**2 * t_fac / two / costhw * spsi0
    ghpsi0ahah = gp**2 / two * (one + tp_fac**2) * spsi0
    ghpsi0zah =  - three * g*gp * spsi0 * tp_fac / two / costhw
    ghpsi0zhah = g*gp * spsi0 / four * (ttp_fac + two**3 * t_fac * tp_fac)
    !!!
    ghpsipwa = - e*g * (spsip - sqrt(two) * spsi0) / two 
    ghpsipwha = - ghpsipwa * t_fac
    ghpsipwz = g**2 / costhw / two * (spsip * sin2thw - sqrt(two) * spsi0 * &
         (one + sin2thw))
    ghpsiwhz = - ghpsipwz * t_fac
    ghpsipwah = - g*gp * (spsip - two*sqrt(two)*spsi0) * tp_fac / two
    ghpsipwhah = - g*gp * (ttp_fac*spsip + & 
         four*sqrt(two)*t_fac*tp_fac*spsi0) / four
    ghpsipwzh = g**2 * t_fac * spsi0 / two 
    ghpsipwhzh = - g**2 * c4s4 / two * spsi0
    ghpsippww = sqrt(two) * g**2 * spsi0
    ghpsippwhwh = sqrt(two) * g**2 * c4s4 * spsi0
    ghpsippwhw = - sqrt(two) * g**2 * t_fac * spsi0
    gpsi00zh =  two * g**2 * t_fac**2
    gpsi00ah = two * gp**2 * tp_fac**2
    gpsi00zhah = two * g*gp * t_fac * tp_fac
    !!!
    gpsi0pwa = - e * g / sqrt(two)
    gpsi0pwha = - gpsi0pwa * t_fac 
    gpsi0pwz = - g**2 * (one + sin2thw) / costhw / sqrt(two)
    gpsi0pwhz = - gpsi0pwz * t_fac
    gpsi0pwah = sqrt(two) * g * gp * tp_fac 
    gpsi0pwhah = - gpsi0pwah * t_fac
    gpsi0pwzh = g**2 * t_fac / sqrt(two)
    gpsi0pwhzh = - g**2 * c4s4 / sqrt(two)
    gpsi0ppww = sqrt(two) * g**2
    gpsi0ppwhwh = gpsi0ppww * c4s4
    gpsi0ppwhw = - gpsi0ppww * t_fac
    i_gpsi0pwa = imago * gpsi0pwa
    i_gpsi0pwha = imago * gpsi0pwha
    i_gpsi0pwz = imago * gpsi0pwz
    i_gpsi0pwhz = imago * gpsi0pwhz
    i_gpsi0pwah = imago * gpsi0pwah
    i_gpsi0pwhah = imago * gpsi0pwhah
    i_gpsi0pwzh = imago * gpsi0pwzh
    i_gpsi0pwhzh = imago * gpsi0pwhzh
    i_gpsi0ppww = imago * gpsi0ppww
    i_gpsi0ppwhwh = imago * gpsi0ppwhwh
    i_gpsi0ppwhw = imago * gpsi0ppwhw
    !!!
    gpsippzz = two * g**2 / costhw**2 * sin2thw**2
    gpsippzhzh = - two * g**2 / four / sin2t / cos2t
    gpsippaz = - two * e * g / costhw * sin2thw
    gpsippaah = - two * e * gp * tp_fac
    gpsippzah = two * g * gp * tp_fac * sin2thw / costhw 
    !!!
    gpsippwa = three * e * g
    gpsippwha = - gpsippwa * t_fac
    gpsippwz = g**2 * (one - three * sin2thw) / costhw
    gpsippwhz = - gpsippwz * t_fac
    gpsippwah = two * g * gp * tp_fac
    gpsippwhah = - gpsippwah * t_fac
    gpsippwzh = - g**2 * t_fac / sqrt(two)
    gpsippwhzh = g**2 * c4s4 / sqrt(two)
    !!!
    gpsicczz = g**2 * (one - two * sin2thw)**2 / costhw**2 
    gpsiccaz = four * e * g * (one - two * sin2thw) / costhw
    gpsiccaah = - four * e * gp * tp_fac
    gpsicczzh = two * g**2 * t_fac * (one - two * sin2thw) / costhw
    gpsiccazh = four * e * g * t_fac
    gpsicczah = - two * g * gp * tp_fac * (one - two * sin2thw) / costhw
    !!! heavy triple gauge couplings
    igahww = - imago * g * costhw * xzbp * vev**2 / f_vev**2
    igzhww = - imago * g * (costhw*xzwp + sint*cost*(cos2t-sin2t)) * vev**2/f_vev**2
    igzwhw = imago * g * xzwp * vev**2 / f_vev**2
    igahwhwh = - imago * g * (two*t_fac*xh + costhw*xzbp) * vev**2 / f_vev**2
    igzhwhwh = - imago * g * two * t_fac
    igahwhw = imago * g * xh * vev**2 / f_vev**2
    !!!
    gwh4 = g**2 * (cos2t**3 + sin2t**3)/sin2t/cos2t 
    gwhwhww = g**2 / four
    gwhwww =  g**2 * sint * cost * (cos2t - sin2t) * vev**2 / four / f_vev**2
    gwh3w = - g**2 * t_fac
    !!!
    gwwaah = - g**2*sinthw*costhw*xzbp*vev**2/f_vev**2
    gwwazh = - g**2*sinthw*costhw*xzwp*vev**2/f_vev**2 + g**2*sinthw*sint*cost* &
         (cos2t-sin2t)/two * vev**2/f_vev**2
    gwwzzh = - g**2 * (costhw**2 - sin2thw) * xzwp * vev**2/f_vev**2
    gwwzah = - g**2 * costhw**2 * xzbp * vev**2/f_vev**2
    gwhwhaah = gwwaah  - g**2*sinthw*two*t_fac*xh*vev**2/f_vev**2 
    gwhwhazh = - g**2*sinthw*two*t_fac
    gwhwhzzh = - g**2*costhw*two*t_fac
    gwhwhzah = gwwzah - g**2*costhw*xh*two*t_fac*vev**2/f_vev**2
    gwwzhah = g**2 * xh * vev**2/f_vev**2
    gwhwhzhah = gwh4*xh*vev**2/f_vev**2 + g**2*costhw*xzbp**two*t_fac*vev**2/f_vev**2
    gwhwzz = g**2*two*costhw*xzwp*vev**2/f_vev**2
    gwhwaz = g**2*sinthw*xzwp*vev**2/two/f_vev**2 
    gwhwaah = g**2*sinthw*xh*vev**2/f_vev**2 
    gwhwzah = g**2*costhw*xh*vev**2/f_vev**2 
    gwhwzhzh = - g**2 * two * t_fac
    gwhwzhah = - g**2*vev**2/f_vev**2 * (xh*two*t_fac + costhw*xzbp)
    gwhwazh = g**2 * sinthw 
    gwhwzzh = g**2 * costhw 
    !!! 
    qzup = g / costhw * qeup * sin2thw
    gcch = - gcc * cost / sint
    gcctop = gcc * (one - (vev * xlam / f_vev)**2 / two)
    gccw = gcc * vev / f_vev * xlam
    gccwh = - gccw * cost / sint
    gnch = g * cost / four / sint
    gztht = - g * xlam * vev / four / costhw / f_vev
    gzhtht = g * xlam * (vev/f_vev) * cost / sint
    !!! Here we  give the formulae for the masses as a function of other
    !!! input parameter in the case that they do not want to be given by
    !!! the user
    ! mass(32) = sqrt((f_vev*gp/sintp/costp)**2/20_default - (gp*vev/two)**2 &
    !            + (g*vev/two/sint/cost)**2 * xh)
    ! mass(33) = sqrt((f_vev*g/two/sint/cost)**2 - (g*vev/two)**2 - &
    !            (gp*vev/two/sintp/costp)**2 * xh)
    ! mass(34) = sqrt((f_vev*g/two/sint/cost)**2 - (g*vev/two)**2)
    gah = gp / two / sintp / costp
    gnchup(1) = gah * (two*yu + 17.0_default/15. - five/6. * cos2tp)
    gnchup(2) = - gah * (one/five - cos2tp/two)
    gnchdwn(1) = gah * (two*yu + 11.0_default/15. + one/6. * cos2tp)
    gnchdwn(2) = gah * (one/five - cos2tp/two)    
    gnchneu(1) = gah * (ye - four/five + one/two * cos2tp)
    gnchneu(2) = - gah * (- ye + four/five - cos2tp/two)    
    gnchlep(1) = gah * (two*ye - 9.0_default/five + three/two * cos2tp)
    gnchlep(2) = gah * (one/five - cos2tp/two)    
    gahtht = gah / five * lam1 * lam2 / sqrt(lam1**2 + lam2**2)
    gahtt(1) = gnchup(1) - gah * xlam / five
    gahtt(2) = gnchup(2) + gah * xlam / five
    gahthth(1) = gah * (two * yu + 14.0_default/15. - four/three * &
         cos2tp + xlam / five)
    gahthth(2) = - gah * xlam / five
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Littlest Higgs Yukawa couplings
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ghtt = - mass(6) / vev * ( one + (vev/f_vev)**2 * xlam * (one + xlam))
    ghthth = - xlam * (one + xlam) * mass(8) * vev / f_vev**2
    ghtht(1) = - mass(6) / two / f_vev * (one + xlam)
    ghtht(2) = - mass(8) * xlam / two / f_vev
    gpsi0tt = - mass(6)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi0)
    gpsi0bb = - mass(5)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi0)
    gpsi0cc = - mass(4)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi0)
    gpsi0tautau = -mass(15)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi0)
    gpsi0tt = imago * mass(6)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi1)
    gpsi0bb = - imago * mass(5)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi1)
    gpsi0cc = imago * mass(4)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi1)
    gpsi0tautau = - imago * mass(15)/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi1)
    gpsipq2(1) = - mass(6)/sqrt(two)/vev * (vev/f_vev - two * spsip) 
    gpsipq2(2) = - mass(5)/sqrt(two)/vev * (vev/f_vev - two * spsip) 
    gpsipq3(1) = - mass(4)/sqrt(two)/vev * (vev/f_vev - two * spsip) 
    gpsipq3(2) = - mass(3)/sqrt(two)/vev * (vev/f_vev - two * spsip) 
    gpsipl3 = - mass(15)/two/sqrt(two)/vev * (vev/f_vev - two * spsip)
    gpsi0tth = - mass(6)/two/sqrt(two)/vev * (vev/f_vev - sqrt(two) * spsi0) &
         * lam1/lam2
    gpsi1tth = imago * gpsi0tth
    gpsipbth = - mass(6)/two/sqrt(two)/vev * (vev/f_vev - two * spsip) &
         * lam1/lam2
    !!!
    ghhtt = two*mass(6)/f_vev**2 * (one - two*f_vev*vevp/vev**2 - xlam/two)
    ghhthth = - lam1**2/mass(8)
    ghhtht(1) = ghhtt * lam1 / lam2
    ghhtht(2) = - mass(6) / vev / f_vev
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default * PI * alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs        
  end subroutine model_update_alpha_s
end module parameters_littlest_tpar
