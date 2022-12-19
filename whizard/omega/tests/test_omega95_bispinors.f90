!  omegalib.nw --
!
!  Copyright (C) 1999-2022 by
!      Wolfgang Kilian <kilian@physik.uni-siegen.de>
!      Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!      Juergen Reuter <juergen.reuter@desy.de>
!      with contributions from                                                                                                                                    
!      Fabian Bach <fabian.bach@t-online.de>                                                                                                                 
!      Bijan Chokoufe Nejad <bijan.chokoufe@desy.de>                                                                                                              
!      Christian Speckner <cnspeckn@googlemail.com>     
!
!  WHIZARD is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  WHIZARD is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_omega95_bispinors
  use kinds
  use omega95_bispinors
  use omega_vspinor_polarizations
  use omega_testtools
  implicit none
  integer :: i, j
  real(kind=default) :: m, pabs, qabs, tabs, zabs, w
  real(kind=default), dimension(4) :: r
  complex(kind=default) :: c_nil, c_one, c_two
  type(momentum) :: p, q, t, z, p_0
  type(vector) :: vp, vq, vt, vz
  type(vectorspinor) :: testv
  type(bispinor) :: vv
  logical :: passed
  call random_seed ()
  c_nil = 0.0_default
  c_one = 1.0_default
  c_two = 2.0_default
  w = 1.4142
  m = 13
  pabs = 42
  qabs = 137
  tabs = 84
  zabs = 3.1415
  p_0%t = m
  p_0%x = 0
  call random_momentum (p, pabs, m)
  call random_momentum (q, qabs, m)
  call random_momentum (t, tabs, m)
  call random_momentum (z, zabs, m)
  call random_number (r)
  do i = 1, 4
     testv%psi(1)%a(i) = (0.0_default, 0.0_default)
  end do
  do i = 2, 3
     do j = 1, 4
        testv%psi(i)%a(j) = cmplx (10.0_default * r(j), kind=default)
    end do
  end do
  testv%psi(4)%a(1) = (1.0_default, 0.0_default)
  testv%psi(4)%a(2) = (0.0_default, 2.0_default)
  testv%psi(4)%a(3) = (1.0_default, 0.0_default)
  testv%psi(4)%a(4) = (3.0_default, 0.0_default)
  vp = p
  vq = q
  vt = t
  vz = z
  passed = .true.
  vv%a(1) = (1.0_default, 0.0_default)
  vv%a(2) = (0.0_default, 2.0_default)
  vv%a(3) = (1.0_default, 0.0_default)
  vv%a(4) = (3.0_default, 0.0_default)
  vv = pr_psi(p, m, w, .false., vv)
  print *, "*** Checking the equations of motion ***:"
  call expect (abs(f_vf(c_one,vp,u(m,p,+1))-m*u(m,p,+1)), 0, "|[p-m]u(+)|=0", passed)
  call expect (abs(f_vf(c_one,vp,u(m,p,-1))-m*u(m,p,-1)), 0, "|[p-m]u(-)|=0", passed)
  call expect (abs(f_vf(c_one,vp,v(m,p,+1))+m*v(m,p,+1)), 0, "|[p+m]v(+)|=0", passed)
  call expect (abs(f_vf(c_one,vp,v(m,p,-1))+m*v(m,p,-1)), 0, "|[p+m]v(-)|=0", passed)
  print *, "*** Checking the equations of motion for negative masses***:"
  call expect (abs(f_vf(c_one,vp,u(-m,p,+1))+m*u(-m,p,+1)), 0, "|[p+m]u(+)|=0", passed)
  call expect (abs(f_vf(c_one,vp,u(-m,p,-1))+m*u(-m,p,-1)), 0, "|[p+m]u(-)|=0", passed)
  call expect (abs(f_vf(c_one,vp,v(-m,p,+1))-m*v(-m,p,+1)), 0, "|[p-m]v(+)|=0", passed)
  call expect (abs(f_vf(c_one,vp,v(-m,p,-1))-m*v(-m,p,-1)), 0, "|[p-m]v(-)|=0", passed)
  print *, "*** Checking the normalization ***:"
  call expect (s_ff(c_one,v(m,p,+1),u(m,p,+1)), +2*m, "ubar(+)*u(+)=+2m", passed)
  call expect (s_ff(c_one,v(m,p,-1),u(m,p,-1)), +2*m, "ubar(-)*u(-)=+2m", passed)
  call expect (s_ff(c_one,u(m,p,+1),v(m,p,+1)), -2*m, "vbar(+)*v(+)=-2m", passed)
  call expect (s_ff(c_one,u(m,p,-1),v(m,p,-1)), -2*m, "vbar(-)*v(-)=-2m", passed)
  call expect (s_ff(c_one,v(m,p,+1),v(m,p,+1)),    0, "ubar(+)*v(+)=0  ", passed)
  call expect (s_ff(c_one,v(m,p,-1),v(m,p,-1)),    0, "ubar(-)*v(-)=0  ", passed)
  call expect (s_ff(c_one,u(m,p,+1),u(m,p,+1)),    0, "vbar(+)*u(+)=0  ", passed)
  call expect (s_ff(c_one,u(m,p,-1),u(m,p,-1)),    0, "vbar(-)*u(-)=0  ", passed)
  print *, "*** Checking the normalization for negative masses***:"
  call expect (s_ff(c_one,v(-m,p,+1),u(-m,p,+1)), -2*m, "ubar(+)*u(+)=-2m", passed)
  call expect (s_ff(c_one,v(-m,p,-1),u(-m,p,-1)), -2*m, "ubar(-)*u(-)=-2m", passed)
  call expect (s_ff(c_one,u(-m,p,+1),v(-m,p,+1)), +2*m, "vbar(+)*v(+)=+2m", passed)
  call expect (s_ff(c_one,u(-m,p,-1),v(-m,p,-1)), +2*m, "vbar(-)*v(-)=+2m", passed)
  call expect (s_ff(c_one,v(-m,p,+1),v(-m,p,+1)),    0, "ubar(+)*v(+)=0  ", passed)
  call expect (s_ff(c_one,v(-m,p,-1),v(-m,p,-1)),    0, "ubar(-)*v(-)=0  ", passed)
  call expect (s_ff(c_one,u(-m,p,+1),u(-m,p,+1)),    0, "vbar(+)*u(+)=0  ", passed)
  call expect (s_ff(c_one,u(-m,p,-1),u(-m,p,-1)),    0, "vbar(-)*u(-)=0  ", passed)
  print *, "*** Checking the currents ***:"
  call expect (abs(v_ff(c_one,v(m,p,+1),u(m,p,+1))-2*vp), 0, "ubar(+).V.u(+)=2p", passed)
  call expect (abs(v_ff(c_one,v(m,p,-1),u(m,p,-1))-2*vp), 0, "ubar(-).V.u(-)=2p", passed)
  call expect (abs(v_ff(c_one,u(m,p,+1),v(m,p,+1))-2*vp), 0, "vbar(+).V.v(+)=2p", passed)
  call expect (abs(v_ff(c_one,u(m,p,-1),v(m,p,-1))-2*vp), 0, "vbar(-).V.v(-)=2p", passed)
  print *, "*** Checking the currents for negative masses***:"
  call expect (abs(v_ff(c_one,v(-m,p,+1),u(-m,p,+1))-2*vp), 0, "ubar(+).V.u(+)=2p", passed)
  call expect (abs(v_ff(c_one,v(-m,p,-1),u(-m,p,-1))-2*vp), 0, "ubar(-).V.u(-)=2p", passed)
  call expect (abs(v_ff(c_one,u(-m,p,+1),v(-m,p,+1))-2*vp), 0, "vbar(+).V.v(+)=2p", passed)
  call expect (abs(v_ff(c_one,u(-m,p,-1),v(-m,p,-1))-2*vp), 0, "vbar(-).V.v(-)=2p", passed)
  print *, "*** Checking current conservation ***:"
  call expect ((vp-vq)*v_ff(c_one,v(m,p,+1),u(m,q,+1)), 0, "d(ubar(+).V.u(+))=0", passed)
  call expect ((vp-vq)*v_ff(c_one,v(m,p,-1),u(m,q,-1)), 0, "d(ubar(-).V.u(-))=0", passed)
  call expect ((vp-vq)*v_ff(c_one,u(m,p,+1),v(m,q,+1)), 0, "d(vbar(+).V.v(+))=0", passed)
  call expect ((vp-vq)*v_ff(c_one,u(m,p,-1),v(m,q,-1)), 0, "d(vbar(-).V.v(-))=0", passed)
  print *, "*** Checking current conservation for negative masses***:"
  call expect ((vp-vq)*v_ff(c_one,v(-m,p,+1),u(-m,q,+1)), 0, "d(ubar(+).V.u(+))=0", passed)
  call expect ((vp-vq)*v_ff(c_one,v(-m,p,-1),u(-m,q,-1)), 0, "d(ubar(-).V.u(-))=0", passed)
  call expect ((vp-vq)*v_ff(c_one,u(-m,p,+1),v(-m,q,+1)), 0, "d(vbar(+).V.v(+))=0", passed)
  call expect ((vp-vq)*v_ff(c_one,u(-m,p,-1),v(-m,q,-1)), 0, "d(vbar(-).V.v(-))=0", passed)
  if (m == 0) then
     print *, "*** Checking axial current conservation ***:"
     call expect ((vp-vq)*a_ff(c_one,v(m,p,+1),u(m,q,+1)), 0, "d(ubar(+).A.u(+))=0", passed)
     call expect ((vp-vq)*a_ff(c_one,v(m,p,-1),u(m,q,-1)), 0, "d(ubar(-).A.u(-))=0", passed)
     call expect ((vp-vq)*a_ff(c_one,u(m,p,+1),v(m,q,+1)), 0, "d(vbar(+).A.v(+))=0", passed)
     call expect ((vp-vq)*a_ff(c_one,u(m,p,-1),v(m,q,-1)), 0, "d(vbar(-).A.v(-))=0", passed)
  end if
  print *, "*** Checking implementation of the sigma vertex funktions ***:"
  call expect ((vp*tvam_ff(c_one,c_nil,v(m,p,+1),u(m,q,+1),q) - (p*q-m**2)*(v(m,p,+1)*u(m,q,+1))), 0, &
               "p*[ubar(p,+).(Isigma*q).u(q,+)] - (p*q-m^2)*ubar(p,+).u(q,+) = 0", passed)
  call expect ((vp*tvam_ff(c_one,c_nil,v(m,p,-1),u(m,q,-1),q) - (p*q-m**2)*(v(m,p,-1)*u(m,q,-1))), 0, &
               "p*[ubar(p,-).(Isigma*q).u(q,-)] - (p*q-m^2)*ubar(p,-).u(q,-) = 0", passed)
  call expect ((vp*tvam_ff(c_one,c_nil,u(m,p,+1),v(m,q,+1),q) - (p*q-m**2)*(u(m,p,+1)*v(m,q,+1))), 0, &
               "p*[vbar(p,+).(Isigma*q).v(q,+)] - (p*q-m^2)*vbar(p,+).v(q,+) = 0", passed)
  call expect ((vp*tvam_ff(c_one,c_nil,u(m,p,-1),v(m,q,-1),q) - (p*q-m**2)*(u(m,p,-1)*v(m,q,-1))), 0, &
               "p*[vbar(p,-).(Isigma*q).v(q,-)] - (p*q-m^2)*vbar(p,-).v(q,-) = 0", passed)
  call expect ((v(m,p,+1)*f_tvamf(c_one,c_nil,vp,u(m,q,+1),q) - (p*q-m**2)*(v(m,p,+1)*u(m,q,+1))), 0, &
               "ubar(p,+).[p*(Isigma*q).u(q,+)] - (p*q-m^2)*ubar(p,+).u(q,+) = 0", passed)
  call expect ((v(m,p,-1)*f_tvamf(c_one,c_nil,vp,u(m,q,-1),q) - (p*q-m**2)*(v(m,p,-1)*u(m,q,-1))), 0, &
               "ubar(p,-).[p*(Isigma*q).u(q,-)] - (p*q-m^2)*ubar(p,-).u(q,-) = 0", passed)
  call expect ((u(m,p,+1)*f_tvamf(c_one,c_nil,vp,v(m,q,+1),q) - (p*q-m**2)*(u(m,p,+1)*v(m,q,+1))), 0, &
               "vbar(p,+).[p*(Isigma*q).v(q,+)] - (p*q-m^2)*vbar(p,+).v(q,+) = 0", passed)
  call expect ((u(m,p,-1)*f_tvamf(c_one,c_nil,vp,v(m,q,-1),q) - (p*q-m**2)*(u(m,p,-1)*v(m,q,-1))), 0, &
               "vbar(p,-).[p*(Isigma*q).v(q,-)] - (p*q-m^2)*vbar(p,-).v(q,-) = 0", passed)

  call expect ((vp*tvam_ff(c_nil,c_one,v(m,p,+1),u(m,q,+1),q) - (p*q+m**2)*p_ff(c_one,v(m,p,+1),u(m,q,+1))), 0, &
               "p*[ubar(p,+).(Isigma*q).g5.u(q,+)] - (p*q+m^2)*ubar(p,+).g5.u(q,+) = 0", passed)
  call expect ((vp*tvam_ff(c_nil,c_one,v(m,p,-1),u(m,q,-1),q) - (p*q+m**2)*p_ff(c_one,v(m,p,-1),u(m,q,-1))), 0, &
               "p*[ubar(p,-).(Isigma*q).g5.u(q,-)] - (p*q+m^2)*ubar(p,-).g5.u(q,-) = 0", passed)
  call expect ((vp*tvam_ff(c_nil,c_one,u(m,p,+1),v(m,q,+1),q) - (p*q+m**2)*p_ff(c_one,u(m,p,+1),v(m,q,+1))), 0, &
               "p*[vbar(p,+).(Isigma*q).g5.v(q,+)] - (p*q+m^2)*vbar(p,+).g5.v(q,+) = 0", passed)
  call expect ((vp*tvam_ff(c_nil,c_one,u(m,p,-1),v(m,q,-1),q) - (p*q+m**2)*p_ff(c_one,u(m,p,-1),v(m,q,-1))), 0, &
               "p*[vbar(p,-).(Isigma*q).g5.v(q,-)] - (p*q+m^2)*vbar(p,-).g5.v(q,-) = 0", passed)
  call expect ((v(m,p,+1)*f_tvamf(c_nil,c_one,vp,u(m,q,+1),q) - (p*q+m**2)*p_ff(c_one,v(m,p,+1),u(m,q,+1))), 0, &
               "p*[ubar(p,+).(Isigma*q).g5.u(q,+)] - (p*q+m^2)*ubar(p,+).g5.u(q,+) = 0", passed)
  call expect ((v(m,p,-1)*f_tvamf(c_nil,c_one,vp,u(m,q,-1),q) - (p*q+m**2)*p_ff(c_one,v(m,p,-1),u(m,q,-1))), 0, &
               "p*[ubar(p,-).(Isigma*q).g5.u(q,-)] - (p*q+m^2)*ubar(p,-).g5.u(q,-) = 0", passed)
  call expect ((u(m,p,+1)*f_tvamf(c_nil,c_one,vp,v(m,q,+1),q) - (p*q+m**2)*p_ff(c_one,u(m,p,+1),v(m,q,+1))), 0, &
               "p*[vbar(p,+).(Isigma*q).g5.v(q,+)] - (p*q+m^2)*vbar(p,+).g5.v(q,+) = 0", passed)
  call expect ((u(m,p,-1)*f_tvamf(c_nil,c_one,vp,v(m,q,-1),q) - (p*q+m**2)*p_ff(c_one,u(m,p,-1),v(m,q,-1))), 0, &
               "p*[vbar(p,-).(Isigma*q).g5.v(q,-)] - (p*q+m^2)*vbar(p,-).g5.v(q,-) = 0", passed)
  print *, "*** Checking polarization vectors: ***"
  call expect (conjg(eps(m,p, 1))*eps(m,p, 1), -1, "e( 1).e( 1)=-1", passed)
  call expect (conjg(eps(m,p, 1))*eps(m,p,-1),  0, "e( 1).e(-1)= 0", passed)
  call expect (conjg(eps(m,p,-1))*eps(m,p, 1),  0, "e(-1).e( 1)= 0", passed)
  call expect (conjg(eps(m,p,-1))*eps(m,p,-1), -1, "e(-1).e(-1)=-1", passed)
  call expect (                 p*eps(m,p, 1),  0, "    p.e( 1)= 0", passed)
  call expect (                 p*eps(m,p,-1),  0, "    p.e(-1)= 0", passed)
  if (m > 0) then
     call expect (conjg(eps(m,p, 1))*eps(m,p, 0),  0, "e( 1).e( 0)= 0", passed)
     call expect (conjg(eps(m,p, 0))*eps(m,p, 1),  0, "e( 0).e( 1)= 0", passed)
     call expect (conjg(eps(m,p, 0))*eps(m,p, 0), -1, "e( 0).e( 0)=-1", passed)
     call expect (conjg(eps(m,p, 0))*eps(m,p,-1),  0, "e( 0).e(-1)= 0", passed)
     call expect (conjg(eps(m,p,-1))*eps(m,p, 0),  0, "e(-1).e( 0)= 0", passed)
     call expect (                 p*eps(m,p, 0),  0, "    p.e( 0)= 0", passed)
  end if
  print *, "*** Checking polarization vectorspinors: ***"
  call expect (abs(p * ueps(m, p,  2)),  0, "p.ueps ( 2)= 0", passed)
  call expect (abs(p * ueps(m, p,  1)),  0, "p.ueps ( 1)= 0", passed)
  call expect (abs(p * ueps(m, p, -1)),  0, "p.ueps (-1)= 0", passed)
  call expect (abs(p * ueps(m, p, -2)),  0, "p.ueps (-2)= 0", passed)
  call expect (abs(p * veps(m, p,  2)),  0, "p.veps ( 2)= 0", passed)
  call expect (abs(p * veps(m, p,  1)),  0, "p.veps ( 1)= 0", passed)
  call expect (abs(p * veps(m, p, -1)),  0, "p.veps (-1)= 0", passed)
  call expect (abs(p * veps(m, p, -2)),  0, "p.veps (-2)= 0", passed)
  print *, "*** Checking polarization vectorspinors (neg. masses): ***"
  call expect (abs(p * ueps(-m, p,  2)),  0, "p.ueps ( 2)= 0", passed)
  call expect (abs(p * ueps(-m, p,  1)),  0, "p.ueps ( 1)= 0", passed)
  call expect (abs(p * ueps(-m, p, -1)),  0, "p.ueps (-1)= 0", passed)
  call expect (abs(p * ueps(-m, p, -2)),  0, "p.ueps (-2)= 0", passed)
  call expect (abs(p * veps(-m, p,  2)),  0, "p.veps ( 2)= 0", passed)
  call expect (abs(p * veps(-m, p,  1)),  0, "p.veps ( 1)= 0", passed)
  call expect (abs(p * veps(-m, p, -1)),  0, "p.veps (-1)= 0", passed)
  call expect (abs(p * veps(-m, p, -2)),  0, "p.veps (-2)= 0", passed)
  print *, "*** in the rest frame ***"
  call expect (abs(p_0 * ueps(m, p_0,  2)),  0, "p0.ueps ( 2)= 0", passed)
  call expect (abs(p_0 * ueps(m, p_0,  1)),  0, "p0.ueps ( 1)= 0", passed)
  call expect (abs(p_0 * ueps(m, p_0, -1)),  0, "p0.ueps (-1)= 0", passed)
  call expect (abs(p_0 * ueps(m, p_0, -2)),  0, "p0.ueps (-2)= 0", passed)
  call expect (abs(p_0 * veps(m, p_0,  2)),  0, "p0.veps ( 2)= 0", passed)
  call expect (abs(p_0 * veps(m, p_0,  1)),  0, "p0.veps ( 1)= 0", passed)
  call expect (abs(p_0 * veps(m, p_0, -1)),  0, "p0.veps (-1)= 0", passed)
  call expect (abs(p_0 * veps(m, p_0, -2)),  0, "p0.veps (-2)= 0", passed)
  print *, "*** in the rest frame (neg. masses) ***"
  call expect (abs(p_0 * ueps(-m, p_0,  2)),  0, "p0.ueps ( 2)= 0", passed)
  call expect (abs(p_0 * ueps(-m, p_0,  1)),  0, "p0.ueps ( 1)= 0", passed)
  call expect (abs(p_0 * ueps(-m, p_0, -1)),  0, "p0.ueps (-1)= 0", passed)
  call expect (abs(p_0 * ueps(-m, p_0, -2)),  0, "p0.ueps (-2)= 0", passed)
  call expect (abs(p_0 * veps(-m, p_0,  2)),  0, "p0.veps ( 2)= 0", passed)
  call expect (abs(p_0 * veps(-m, p_0,  1)),  0, "p0.veps ( 1)= 0", passed)
  call expect (abs(p_0 * veps(-m, p_0, -1)),  0, "p0.veps (-1)= 0", passed)
  call expect (abs(p_0 * veps(-m, p_0, -2)),  0, "p0.veps (-2)= 0", passed)
  print *, "*** Checking the irreducibility condition: ***"
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p,  2))),  0, "g.ueps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p,  1))),  0, "g.ueps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p, -1))),  0, "g.ueps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p, -2))),  0, "g.ueps (-2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p,  2))),  0, "g.veps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p,  1))),  0, "g.veps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p, -1))),  0, "g.veps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p, -2))),  0, "g.veps (-2)", passed)
  print *, "*** Checking the irreducibility condition (neg. masses): ***"
  call expect (abs(f_potgr (c_one, c_one, ueps(-m, p,  2))),  0, "g.ueps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(-m, p,  1))),  0, "g.ueps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(-m, p, -1))),  0, "g.ueps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(-m, p, -2))),  0, "g.ueps (-2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(-m, p,  2))),  0, "g.veps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(-m, p,  1))),  0, "g.veps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(-m, p, -1))),  0, "g.veps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(-m, p, -2))),  0, "g.veps (-2)", passed)
  print *, "*** in the rest frame ***"
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0,  2))),  0, "g.ueps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0,  1))),  0, "g.ueps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0, -1))),  0, "g.ueps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0, -2))),  0, "g.ueps (-2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0,  2))),  0, "g.veps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0,  1))),  0, "g.veps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0, -1))),  0, "g.veps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0, -2))),  0, "g.veps (-2)", passed)
  print *, "*** in the rest frame (neg. masses) ***"
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0,  2))),  0, "g.ueps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0,  1))),  0, "g.ueps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0, -1))),  0, "g.ueps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, ueps(m, p_0, -2))),  0, "g.ueps (-2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0,  2))),  0, "g.veps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0,  1))),  0, "g.veps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0, -1))),  0, "g.veps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, veps(m, p_0, -2))),  0, "g.veps (-2)", passed)
  print *, "*** Testing vectorspinor normalization ***"
  call expect (veps(m,p, 2)*ueps(m,p, 2), -2*m, "ueps( 2).ueps( 2)= -2m", passed)
  call expect (veps(m,p, 1)*ueps(m,p, 1), -2*m, "ueps( 1).ueps( 1)= -2m", passed)
  call expect (veps(m,p,-1)*ueps(m,p,-1), -2*m, "ueps(-1).ueps(-1)= -2m", passed)
  call expect (veps(m,p,-2)*ueps(m,p,-2), -2*m, "ueps(-2).ueps(-2)= -2m", passed)
  call expect (ueps(m,p, 2)*veps(m,p, 2),  2*m, "veps( 2).veps( 2)= +2m", passed)
  call expect (ueps(m,p, 1)*veps(m,p, 1),  2*m, "veps( 1).veps( 1)= +2m", passed)
  call expect (ueps(m,p,-1)*veps(m,p,-1),  2*m, "veps(-1).veps(-1)= +2m", passed)
  call expect (ueps(m,p,-2)*veps(m,p,-2),  2*m, "veps(-2).veps(-2)= +2m", passed)
  call expect (ueps(m,p, 2)*ueps(m,p, 2),    0, "ueps( 2).veps( 2)=   0", passed)
  call expect (ueps(m,p, 1)*ueps(m,p, 1),    0, "ueps( 1).veps( 1)=   0", passed)
  call expect (ueps(m,p,-1)*ueps(m,p,-1),    0, "ueps(-1).veps(-1)=   0", passed)
  call expect (ueps(m,p,-2)*ueps(m,p,-2),    0, "ueps(-2).veps(-2)=   0", passed)
  call expect (veps(m,p, 2)*veps(m,p, 2),    0, "veps( 2).ueps( 2)=   0", passed)
  call expect (veps(m,p, 1)*veps(m,p, 1),    0, "veps( 1).ueps( 1)=   0", passed)
  call expect (veps(m,p,-1)*veps(m,p,-1),    0, "veps(-1).ueps(-1)=   0", passed)
  call expect (veps(m,p,-2)*veps(m,p,-2),    0, "veps(-2).ueps(-2)=   0", passed)
  print *, "*** Testing vectorspinor normalization (neg. masses) ***"
  call expect (veps(-m,p, 2)*ueps(-m,p, 2), +2*m, "ueps( 2).ueps( 2)= +2m", passed)
  call expect (veps(-m,p, 1)*ueps(-m,p, 1), +2*m, "ueps( 1).ueps( 1)= +2m", passed)
  call expect (veps(-m,p,-1)*ueps(-m,p,-1), +2*m, "ueps(-1).ueps(-1)= +2m", passed)
  call expect (veps(-m,p,-2)*ueps(-m,p,-2), +2*m, "ueps(-2).ueps(-2)= +2m", passed)
  call expect (ueps(-m,p, 2)*veps(-m,p, 2), -2*m, "veps( 2).veps( 2)= -2m", passed)
  call expect (ueps(-m,p, 1)*veps(-m,p, 1), -2*m, "veps( 1).veps( 1)= -2m", passed)
  call expect (ueps(-m,p,-1)*veps(-m,p,-1), -2*m, "veps(-1).veps(-1)= -2m", passed)
  call expect (ueps(-m,p,-2)*veps(-m,p,-2), -2*m, "veps(-2).veps(-2)= -2m", passed)
  call expect (ueps(-m,p, 2)*ueps(-m,p, 2),    0, "ueps( 2).veps( 2)=   0", passed)
  call expect (ueps(-m,p, 1)*ueps(-m,p, 1),    0, "ueps( 1).veps( 1)=   0", passed)
  call expect (ueps(-m,p,-1)*ueps(-m,p,-1),    0, "ueps(-1).veps(-1)=   0", passed)
  call expect (ueps(-m,p,-2)*ueps(-m,p,-2),    0, "ueps(-2).veps(-2)=   0", passed)
  call expect (veps(-m,p, 2)*veps(-m,p, 2),    0, "veps( 2).ueps( 2)=   0", passed)
  call expect (veps(-m,p, 1)*veps(-m,p, 1),    0, "veps( 1).ueps( 1)=   0", passed)
  call expect (veps(-m,p,-1)*veps(-m,p,-1),    0, "veps(-1).ueps(-1)=   0", passed)
  call expect (veps(-m,p,-2)*veps(-m,p,-2),    0, "veps(-2).ueps(-2)=   0", passed)
  print *, "*** in the rest frame ***"
  call expect (veps(m,p_0, 2)*ueps(m,p_0, 2), -2*m, "ueps( 2).ueps( 2)= -2m", passed)
  call expect (veps(m,p_0, 1)*ueps(m,p_0, 1), -2*m, "ueps( 1).ueps( 1)= -2m", passed)
  call expect (veps(m,p_0,-1)*ueps(m,p_0,-1), -2*m, "ueps(-1).ueps(-1)= -2m", passed)
  call expect (veps(m,p_0,-2)*ueps(m,p_0,-2), -2*m, "ueps(-2).ueps(-2)= -2m", passed)
  call expect (ueps(m,p_0, 2)*veps(m,p_0, 2),  2*m, "veps( 2).veps( 2)= +2m", passed)
  call expect (ueps(m,p_0, 1)*veps(m,p_0, 1),  2*m, "veps( 1).veps( 1)= +2m", passed)
  call expect (ueps(m,p_0,-1)*veps(m,p_0,-1),  2*m, "veps(-1).veps(-1)= +2m", passed)
  call expect (ueps(m,p_0,-2)*veps(m,p_0,-2),  2*m, "veps(-2).veps(-2)= +2m", passed)
  call expect (ueps(m,p_0, 2)*ueps(m,p_0, 2),    0, "ueps( 2).veps( 2)=   0", passed)
  call expect (ueps(m,p_0, 1)*ueps(m,p_0, 1),    0, "ueps( 1).veps( 1)=   0", passed)
  call expect (ueps(m,p_0,-1)*ueps(m,p_0,-1),    0, "ueps(-1).veps(-1)=   0", passed)
  call expect (ueps(m,p_0,-2)*ueps(m,p_0,-2),    0, "ueps(-2).veps(-2)=   0", passed)
  call expect (veps(m,p_0, 2)*veps(m,p_0, 2),    0, "veps( 2).ueps( 2)=   0", passed)
  call expect (veps(m,p_0, 1)*veps(m,p_0, 1),    0, "veps( 1).ueps( 1)=   0", passed)
  call expect (veps(m,p_0,-1)*veps(m,p_0,-1),    0, "veps(-1).ueps(-1)=   0", passed)
  call expect (veps(m,p_0,-2)*veps(m,p_0,-2),    0, "veps(-2).ueps(-2)=   0", passed)
  print *, "*** in the rest frame (neg. masses) ***"
  call expect (veps(-m,p_0, 2)*ueps(-m,p_0, 2), +2*m, "ueps( 2).ueps( 2)= +2m", passed)
  call expect (veps(-m,p_0, 1)*ueps(-m,p_0, 1), +2*m, "ueps( 1).ueps( 1)= +2m", passed)
  call expect (veps(-m,p_0,-1)*ueps(-m,p_0,-1), +2*m, "ueps(-1).ueps(-1)= +2m", passed)
  call expect (veps(-m,p_0,-2)*ueps(-m,p_0,-2), +2*m, "ueps(-2).ueps(-2)= +2m", passed)
  call expect (ueps(-m,p_0, 2)*veps(-m,p_0, 2), -2*m, "veps( 2).veps( 2)= -2m", passed)
  call expect (ueps(-m,p_0, 1)*veps(-m,p_0, 1), -2*m, "veps( 1).veps( 1)= -2m", passed)
  call expect (ueps(-m,p_0,-1)*veps(-m,p_0,-1), -2*m, "veps(-1).veps(-1)= -2m", passed)
  call expect (ueps(-m,p_0,-2)*veps(-m,p_0,-2), -2*m, "veps(-2).veps(-2)= -2m", passed)
  call expect (ueps(-m,p_0, 2)*ueps(-m,p_0, 2),    0, "ueps( 2).veps( 2)=   0", passed)
  call expect (ueps(-m,p_0, 1)*ueps(-m,p_0, 1),    0, "ueps( 1).veps( 1)=   0", passed)
  call expect (ueps(-m,p_0,-1)*ueps(-m,p_0,-1),    0, "ueps(-1).veps(-1)=   0", passed)
  call expect (ueps(-m,p_0,-2)*ueps(-m,p_0,-2),    0, "ueps(-2).veps(-2)=   0", passed)
  call expect (veps(-m,p_0, 2)*veps(-m,p_0, 2),    0, "veps( 2).ueps( 2)=   0", passed)
  call expect (veps(-m,p_0, 1)*veps(-m,p_0, 1),    0, "veps( 1).ueps( 1)=   0", passed)
  call expect (veps(-m,p_0,-1)*veps(-m,p_0,-1),    0, "veps(-1).ueps(-1)=   0", passed)
  call expect (veps(-m,p_0,-2)*veps(-m,p_0,-2),    0, "veps(-2).ueps(-2)=   0", passed)
  print *, "*** Majorana properties of gravitino vertices: ***"
  call expect (abs(u (m,q,1) * f_sgr (c_one, c_one, ueps(m,p,2), t) + &
     ueps(m,p,2) * gr_sf(c_one,c_one,u(m,q,1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_sgr (c_one, c_one, ueps(m,p,2), t) + &
  !!!    ueps(m,p,2) * gr_sf(c_one,c_one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  !!! call expect (abs(u (m,q,1) * f_sgr (c_one, c_one, ueps(m,p,1), t) + &
  !!!    ueps(m,p,1) * gr_sf(c_one,c_one,u(m,q,1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_sgr (c_one, c_one, ueps(m,p,1), t) + &
  !!!    ueps(m,p,1) * gr_sf(c_one,c_one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  !!! call expect (abs(u (m,q,1) * f_sgr (c_one, c_one, ueps(m,p,-1), t) + &
  !!!    ueps(m,p,-1) * gr_sf(c_one,c_one,u(m,q,1),t)),  0, "f_sgr   + gr_sf       = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_sgr (c_one, c_one, ueps(m,p,-1), t) + &
  !!!    ueps(m,p,-1) * gr_sf(c_one,c_one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  !!! call expect (abs(u (m,q,1) * f_sgr (c_one, c_one, ueps(m,p,-2), t) + &
  !!!    ueps(m,p,-2) * gr_sf(c_one,c_one,u(m,q,1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_sgr (c_one, c_one, ueps(m,p,-2), t) + &
  !!!    ueps(m,p,-2) * gr_sf(c_one,c_one,u(m,q,-1),t)),  0, "f_sgr     + gr_sf     = 0", passed)
  call expect (abs(u (m,q,1) * f_slgr (c_one, c_one, ueps(m,p,2), t) + &
     ueps(m,p,2) * gr_slf(c_one,c_one,u(m,q,1),t)),  0, "f_slgr    + gr_slf    = 0", passed, threshold = 0.5_default)
  call expect (abs(u (m,q,1) * f_srgr (c_one, c_one, ueps(m,p,2), t) + &
     ueps(m,p,2) * gr_srf(c_one,c_one,u(m,q,1),t)),  0, "f_srgr    + gr_srf    = 0", passed, threshold = 0.5_default)
  call expect (abs(u (m,q,1) * f_slrgr (c_one, c_two, c_one, ueps(m,p,2), t) + &
     ueps(m,p,2) * gr_slrf(c_one,c_two,c_one,u(m,q,1),t)),  0, "f_slrgr   + gr_slrf   = 0", passed, threshold = 0.5_default)
  call expect (abs(u (m,q,1) * f_pgr (c_one, c_one, ueps(m,p,2), t) + &
     ueps(m,p,2) * gr_pf(c_one,c_one,u(m,q,1),t)),  0, "f_pgr     + gr_pf     = 0", passed, threshold = 0.5_default)
  call expect (abs(u (m,q,1) * f_vgr (c_one, vt, ueps(m,p,2), p+q) + &
     ueps(m,p,2) * gr_vf(c_one,vt,u(m,q,1),p+q)),  0, "f_vgr     + gr_vf = 0", passed, threshold = 0.5_default)
  call expect (abs(u (m,q,1) * f_vlrgr (c_one, c_two, vt, ueps(m,p,2), p+q) + &
     ueps(m,p,2) * gr_vlrf(c_one,c_two,vt,u(m,q,1),p+q)),  0, "f_vlrgr   + gr_vlrf   = 0", &
     passed, threshold = 0.5_default)
  !!! call expect (abs(u (m,q,-1) * f_vgr (c_one, vt, ueps(m,p,2), p+q) + &
  !!!    ueps(m,p,2) * gr_vf(c_one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  !!! call expect (abs(u (m,q,1) * f_vgr (c_one, vt, ueps(m,p,1), p+q) + &
  !!!    ueps(m,p,1) * gr_vf(c_one,vt,u(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_vgr (c_one, vt, ueps(m,p,1), p+q) + &
  !!!    ueps(m,p,1) * gr_vf(c_one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  !!! call expect (abs(u (m,q,1) * f_vgr (c_one, vt, ueps(m,p,-1), p+q) + &
  !!!    ueps(m,p,-1) * gr_vf(c_one,vt,u(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_vgr (c_one, vt, veps(m,p,-1), p+q) + &
  !!!    veps(m,p,-1) * gr_vf(c_one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  !!! call expect (abs(v (m,q,1) * f_vgr (c_one, vt, ueps(m,p,-2), p+q) + &
  !!!    ueps(m,p,-2) * gr_vf(c_one,vt,v(m,q,1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  !!! call expect (abs(u (m,q,-1) * f_vgr (c_one, vt, ueps(m,p,-2), p+q) + &
  !!!    ueps(m,p,-2) * gr_vf(c_one,vt,u(m,q,-1),p+q)),  0, "f_vgr     + gr_vf     = 0", passed)
  call expect (abs(s_grf (c_one, ueps(m,p,2), u(m,q,1),t) + &
     s_fgr(c_one,u(m,q,1),ueps(m,p,2),t)),  0, "s_grf     + s_fgr     = 0", passed)
  call expect (abs(sl_grf (c_one, ueps(m,p,2), u(m,q,1),t) + &
     sl_fgr(c_one,u(m,q,1),ueps(m,p,2),t)),  0, "sl_grf    + sl_fgr    = 0", passed)
  call expect (abs(sr_grf (c_one, ueps(m,p,2), u(m,q,1),t) + &
     sr_fgr(c_one,u(m,q,1),ueps(m,p,2),t)),  0, "sr_grf    + sr_fgr    = 0", passed)
  call expect (abs(slr_grf (c_one, c_two, ueps(m,p,2), u(m,q,1),t) + &
     slr_fgr(c_one,c_two,u(m,q,1),ueps(m,p,2),t)),  0, "slr_grf   + slr_fgr   = 0", passed)
  call expect (abs(p_grf (c_one, ueps(m,p,2), u(m,q,1),t) + &
     p_fgr(c_one,u(m,q,1),ueps(m,p,2),t)),  0, "p_grf     + p_fgr     = 0", passed)
  call expect (abs(v_grf (c_one, ueps(m,p,2), u(m,q,1),t) + &
     v_fgr(c_one,u(m,q,1),ueps(m,p,2),t)),  0, "v_grf     + v_fgr     = 0", passed)
  call expect (abs(vlr_grf (c_one, c_two, ueps(m,p,2), u(m,q,1),t) + &
     vlr_fgr(c_one,c_two,u(m,q,1),ueps(m,p,2),t)),  0, "vlr_grf   + vlr_fgr   = 0", passed)
  call expect (abs(u(m,p,1) * f_potgr (c_one,c_one,testv) - testv * gr_potf &
     (c_one,c_one,u (m,p,1))), 0, "f_potgr   - gr_potf   = 0", passed)
  call expect (abs (pot_fgr (c_one,u(m,p,1),testv) - pot_grf(c_one, &
     testv,u(m,p,1))), 0, "pot_fgr   - pot_grf   = 0", passed)
  call expect (abs(u(m,p,1) * f_s2gr (c_one,c_one,c_one,testv) - testv * gr_s2f &
     (c_one,c_one,c_one,u (m,p,1))), 0, "f_s2gr    - gr_s2f    = 0", passed)
  call expect (abs (s2_fgr (c_one,u(m,p,1),c_one,testv) - s2_grf(c_one, &
     testv,c_one,u(m,p,1))), 0, "s2_fgr    - s2_grf    = 0", passed)
  call expect (abs(u (m,q,1) * f_svgr (c_one, c_one, vt, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_svf(c_one,c_one,vt,u(m,q,1))),  0, "f_svgr    + gr_svf    = 0", passed)
  call expect (abs(u (m,q,1) * f_slvgr (c_one, c_one, vt, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_slvf(c_one,c_one,vt,u(m,q,1))),  0, "f_slvgr   + gr_slvf   = 0", passed)
  call expect (abs(u (m,q,1) * f_srvgr (c_one, c_one, vt, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_srvf(c_one,c_one,vt,u(m,q,1))),  0, "f_srvgr   + gr_srvf   = 0", passed)
  call expect (abs(u (m,q,1) * f_slrvgr (c_one, c_two, c_one, vt, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_slrvf(c_one,c_two,c_one,vt,u(m,q,1))),  0, "f_slrvgr  + gr_slrvf  = 0", passed)
  call expect (abs (sv1_fgr (c_one,u(m,p,1),vt,ueps(m,q,2)) + sv1_grf(c_one, &
     ueps(m,q,2),vt,u(m,p,1))), 0, "sv1_fgr   + sv1_grf   = 0", passed)
  call expect (abs (sv2_fgr (c_one,u(m,p,1),c_one,ueps(m,q,2)) + sv2_grf(c_one, &
     ueps(m,q,2),c_one,u(m,p,1))), 0, "sv2_fgr   + sv2_grf   = 0", passed)
  call expect (abs (slv1_fgr (c_one,u(m,p,1),vt,ueps(m,q,2)) + slv1_grf(c_one, &
     ueps(m,q,2),vt,u(m,p,1))), 0, "slv1_fgr  + slv1_grf  = 0", passed)
  call expect (abs (srv2_fgr (c_one,u(m,p,1),c_one,ueps(m,q,2)) + srv2_grf(c_one, &
     ueps(m,q,2),c_one,u(m,p,1))), 0, "srv2_fgr  + srv2_grf  = 0", passed)
  call expect (abs (slrv1_fgr (c_one,c_two,u(m,p,1),vt,ueps(m,q,2)) + slrv1_grf(c_one,c_two, &
     ueps(m,q,2),vt,u(m,p,1))), 0, "slrv1_fgr + slrv1_grf = 0", passed)
  call expect (abs (slrv2_fgr (c_one,c_two,u(m,p,1),c_one,ueps(m,q,2)) + slrv2_grf(c_one, &
     c_two,ueps(m,q,2),c_one,u(m,p,1))), 0, "slrv2_fgr + slrv2_grf = 0", passed)
  call expect (abs(u (m,q,1) * f_pvgr (c_one, c_one, vt, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_pvf(c_one,c_one,vt,u(m,q,1))),  0, "f_pvgr    + gr_pvf    = 0", passed)
  call expect (abs (pv1_fgr (c_one,u(m,p,1),vt,ueps(m,q,2)) + pv1_grf(c_one, &
     ueps(m,q,2),vt,u(m,p,1))), 0, "pv1_fgr   + pv1_grf   = 0", passed)
  call expect (abs (pv2_fgr (c_one,u(m,p,1),c_one,ueps(m,q,2)) + pv2_grf(c_one, &
     ueps(m,q,2),c_one,u(m,p,1))), 0, "pv2_fgr   + pv2_grf   = 0", passed)
  call expect (abs(u (m,q,1) * f_v2gr (c_one, vt, vz, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_v2f(c_one,vt,vz,u(m,q,1))),  0, "f_v2gr    + gr_v2f    = 0", passed)
  call expect (abs(u (m,q,1) * f_v2lrgr (c_one, c_two, vt, vz, ueps(m,p,2)) + &
     ueps(m,p,2) * gr_v2lrf(c_one,c_two,vt,vz,u(m,q,1))),  0, "f_v2lrgr  + gr_v2lrf  = 0", passed)
  call expect (abs (v2_fgr (c_one,u(m,p,1),vt,ueps(m,q,2)) + v2_grf(c_one, &
     ueps(m,q,2),vt,u(m,p,1))), 0, "v2_fgr    + v2_grf    = 0", passed)
  call expect (abs (v2lr_fgr (c_one,c_two,u(m,p,1),vt,ueps(m,q,2)) + v2lr_grf(c_one, c_two, &
     ueps(m,q,2),vt,u(m,p,1))), 0, "v2lr_fgr  + v2lr_grf  = 0", passed)
  print *, "*** Testing the gravitino propagator: ***"
  print *, "Transversality:"
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,testv))), 0, "p.pr.test", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,ueps(m,p,2)))),  0, "p.pr.ueps ( 2)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,ueps(m,p,1)))),  0, "p.pr.ueps ( 1)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,ueps(m,p,-1)))), 0, "p.pr.ueps (-1)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,ueps(m,p,-2)))), 0, "p.pr.ueps (-2)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,veps(m,p,2)))),  0, "p.pr.veps ( 2)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,veps(m,p,1)))),  0, "p.pr.veps ( 1)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,veps(m,p,-1)))), 0, "p.pr.veps (-1)", passed)
  call expect (abs(p * (cmplx (p*p - m**2, m*w, kind=default) * &
               pr_grav(p,m,w,veps(m,p,-2)))), 0, "p.pr.veps (-2)", passed)
  print *, "Irreducibility:"
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,testv)))), 0, "g.pr.test", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,ueps(m,p,2))))), 0, &
               "g.pr.ueps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,ueps(m,p,1))))), 0, &
               "g.pr.ueps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,ueps(m,p,-1))))), 0, &
               "g.pr.ueps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,ueps(m,p,-2))))), 0, &
               "g.pr.ueps (-2)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,veps(m,p,2))))), 0, &
               "g.pr.veps ( 2)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,veps(m,p,1))))), 0, &
               "g.pr.veps ( 1)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,veps(m,p,-1))))), 0, &
               "g.pr.veps (-1)", passed)
  call expect (abs(f_potgr (c_one, c_one, (cmplx (p*p - m**2, m*w, &
               kind=default) * pr_grav(p,m,w,veps(m,p,-2))))), 0, &
               "g.pr.veps (-2)", passed)
  if (.not. passed) then
    stop 1
  end if
end program test_omega95_bispinors
