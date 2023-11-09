! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
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
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module sf_aux_uti

  use kinds, only: default
  use numeric_utils, only: pacify
  use lorentz

  use sf_aux

  implicit none
  private

  public :: sf_aux_1
  public :: sf_aux_2
  public :: sf_aux_3
  public :: sf_aux_4

contains

  subroutine sf_aux_1 (u)
    integer, intent(in) :: u
    type(splitting_data_t) :: sd
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q, q0
    real(default) :: E, mk, mp, mq
    real(default) :: x, r1, r2, r1o, r2o
    real(default) :: k2, q0_2, q1_2, q2_2

    write (u, "(A)")  "* Test output: sf_aux_1"
    write (u, "(A)")  "*   Purpose: compute momentum splitting"
    write (u, "(A)")  "             (massless radiated particle)"
    write (u, "(A)")

    E = 1
    mk = 0.3_default
    mp = 0
    mq = mk

    k = vector4_moving (E, sqrt (E**2 - mk**2), 3)
    k2 = k ** 2;  call pacify (k2, 1e-10_default)

    x = 0.6_default
    r1 = 0.5_default
    r2 = 0.125_default

    write (u, "(A)")  "* (1) Non-collinear setup"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%sample_t (r1)
    call sd%sample_phi (r2)

    call sd%write (u)

    q = sd%split_momentum (k)
    q1_2 = q(1) ** 2;  call pacify (q1_2, 1e-10_default)
    q2_2 = q(2) ** 2;  call pacify (q2_2, 1e-10_default)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: s"
    write (u, "(2(1x,F11.8))")  sd%s, k2

    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  sd%t, q2_2

    write (u, "(A)")  "Compare: u"
    write (u, "(2(1x,F11.8))")  sd%u, q1_2

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  sd%x, energy (q(2)) / energy (k)

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  sd%xb, energy (q(1)) / energy (k)

    write (u, "(A)")
    write (u, "(A)")  "Extract: x, 1-x"
    write (u, "(2(1x,F11.8))")  sd%get_x (), sd%get_xb ()

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep energy)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_ENERGY)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q0_2 = q0(2) ** 2;  call pacify (q0_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q0_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%recover (k, q0, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t


    call sd%inverse_t (r1o)

    write (u, "(A)")  "Compare: r1"
    write (u, "(2(1x,F11.8))")  r1, r1o

    call sd%inverse_phi (r2o)

    write (u, "(A)")  "Compare: r2"
    write (u, "(2(1x,F11.8))")  r2, r2o

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep momentum)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_MOMENTUM)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q0_2 = q0(2) ** 2;  call pacify (q0_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q0_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%recover (k, q0, KEEP_MOMENTUM)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    call sd%inverse_t (r1o)

    write (u, "(A)")  "Compare: r1"
    write (u, "(2(1x,F11.8))")  r1, r1o

    call sd%inverse_phi (r2o)

    write (u, "(A)")  "Compare: r2"
    write (u, "(2(1x,F11.8))")  r2, r2o

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* (2) Collinear setup"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, 1 - x)

    call sd%write (u)

    q = sd%split_momentum (k)
    q1_2 = q(1) ** 2;  call pacify (q1_2, 1e-10_default)
    q2_2 = q(2) ** 2;  call pacify (q2_2, 1e-10_default)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: s"
    write (u, "(2(1x,F11.8))")  sd%s, k2

    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  sd%t, q2_2

    write (u, "(A)")  "Compare: u"
    write (u, "(2(1x,F11.8))")  sd%u, q1_2

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  sd%x, energy (q(2)) / energy (k)

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  sd%xb, energy (q(1)) / energy (k)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep energy)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_ENERGY)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q0_2 = q0(2) ** 2;  call pacify (q0_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q0_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%recover (k, q0, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep momentum)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_MOMENTUM)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q0_2 = q0(2) ** 2;  call pacify (q0_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q0_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%recover (k, q0, KEEP_MOMENTUM)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_aux_1"

  end subroutine sf_aux_1

  subroutine sf_aux_2 (u)
    integer, intent(in) :: u
    type(splitting_data_t) :: sd
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q, q0
    real(default) :: E, mk, mp, mq
    real(default) :: x, r1, r2, r1o, r2o
    real(default) :: k2, q02_2, q1_2, q2_2

    write (u, "(A)")  "* Test output: sf_aux_2"
    write (u, "(A)")  "*   Purpose: compute momentum splitting"
    write (u, "(A)")  "             (massless outgoing particle)"
    write (u, "(A)")

    E = 1
    mk = 0.3_default
    mp = mk
    mq = 0

    k = vector4_moving (E, sqrt (E**2 - mk**2), 3)
    k2 = k ** 2;  call pacify (k2, 1e-10_default)

    x = 0.6_default
    r1 = 0.5_default
    r2 = 0.125_default

    write (u, "(A)")  "* (1) Non-collinear setup"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%sample_t (r1)
    call sd%sample_phi (r2)

    call sd%write (u)

    q = sd%split_momentum (k)
    q1_2 = q(1) ** 2;  call pacify (q1_2, 1e-10_default)
    q2_2 = q(2) ** 2;  call pacify (q2_2, 1e-10_default)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: s"
    write (u, "(2(1x,F11.8))")  sd%s, k2

    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  sd%t, q2_2

    write (u, "(A)")  "Compare: u"
    write (u, "(2(1x,F11.8))")  sd%u, q1_2

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  sd%x, energy (q(2)) / energy (k)

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  sd%xb, energy (q(1)) / energy (k)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep energy)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_ENERGY)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t


    call sd%inverse_t (r1o)

    write (u, "(A)")  "Compare: r1"
    write (u, "(2(1x,F11.8))")  r1, r1o

    call sd%inverse_phi (r2o)

    write (u, "(A)")  "Compare: r2"
    write (u, "(2(1x,F11.8))")  r2, r2o

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep momentum)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_MOMENTUM)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_MOMENTUM)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    call sd%inverse_t (r1o)

    write (u, "(A)")  "Compare: r1"
    write (u, "(2(1x,F11.8))")  r1, r1o

    call sd%inverse_phi (r2o)

    write (u, "(A)")  "Compare: r2"
    write (u, "(2(1x,F11.8))")  r2, r2o

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* (2) Collinear setup"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, 1 - x)

    call sd%write (u)

    q = sd%split_momentum (k)
    q1_2 = q(1) ** 2;  call pacify (q1_2, 1e-10_default)
    q2_2 = q(2) ** 2;  call pacify (q2_2, 1e-10_default)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: s"
    write (u, "(2(1x,F11.8))")  sd%s, k2

    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  sd%t, q2_2

    write (u, "(A)")  "Compare: u"
    write (u, "(2(1x,F11.8))")  sd%u, q1_2

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  sd%x, energy (q(2)) / energy (k)

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  sd%xb, energy (q(1)) / energy (k)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep energy)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_ENERGY)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep momentum)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_MOMENTUM)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_MOMENTUM)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_aux_2"

  end subroutine sf_aux_2

  subroutine sf_aux_3 (u)
    integer, intent(in) :: u
    type(splitting_data_t) :: sd
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q, q0
    real(default) :: E, mk, mp, mq, qmin, qmax
    real(default) :: x, r1, r2, r1o, r2o
    real(default) :: k2, q02_2, q1_2, q2_2

    write (u, "(A)")  "* Test output: sf_aux_3"
    write (u, "(A)")  "*   Purpose: compute momentum splitting"
    write (u, "(A)")  "             (all massless, q cuts)"
    write (u, "(A)")

    E = 1
    mk = 0
    mp = 0
    mq = 0
    qmin = 1e-2_default
    qmax = 1e0_default

    k = vector4_moving (E, sqrt (E**2 - mk**2), 3)
    k2 = k ** 2;  call pacify (k2, 1e-10_default)

    x = 0.6_default
    r1 = 0.5_default
    r2 = 0.125_default

    write (u, "(A)")  "* (1) Non-collinear setup"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%sample_t (r1, t1 = - qmin ** 2, t0 = - qmax **2)
    call sd%sample_phi (r2)

    call sd%write (u)

    q = sd%split_momentum (k)
    q1_2 = q(1) ** 2;  call pacify (q1_2, 1e-10_default)
    q2_2 = q(2) ** 2;  call pacify (q2_2, 1e-10_default)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: s"
    write (u, "(2(1x,F11.8))")  sd%s, k2

    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  sd%t, q2_2

    write (u, "(A)")  "Compare: u"
    write (u, "(2(1x,F11.8))")  sd%u, q1_2

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  sd%x, energy (q(2)) / energy (k)

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  sd%xb, energy (q(1)) / energy (k)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep energy)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_ENERGY)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t


    call sd%inverse_t (r1o, t1 = - qmin ** 2, t0 = - qmax **2)

    write (u, "(A)")  "Compare: r1"
    write (u, "(2(1x,F11.8))")  r1, r1o

    call sd%inverse_phi (r2o)

    write (u, "(A)")  "Compare: r2"
    write (u, "(2(1x,F11.8))")  r2, r2o

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep momentum)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_MOMENTUM)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_MOMENTUM)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    call sd%inverse_t (r1o, t1 = - qmin ** 2, t0 = - qmax **2)

    write (u, "(A)")  "Compare: r1"
    write (u, "(2(1x,F11.8))")  r1, r1o

    call sd%inverse_phi (r2o)

    write (u, "(A)")  "Compare: r2"
    write (u, "(2(1x,F11.8))")  r2, r2o

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* (2) Collinear setup"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, 1 - x)

    call sd%write (u)

    q = sd%split_momentum (k)
    q1_2 = q(1) ** 2;  call pacify (q1_2, 1e-10_default)
    q2_2 = q(2) ** 2;  call pacify (q2_2, 1e-10_default)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: s"
    write (u, "(2(1x,F11.8))")  sd%s, k2

    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  sd%t, q2_2

    write (u, "(A)")  "Compare: u"
    write (u, "(2(1x,F11.8))")  sd%u, q1_2

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  sd%x, energy (q(2)) / energy (k)

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  sd%xb, energy (q(1)) / energy (k)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep energy)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_ENERGY)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Project on-shell (keep momentum)"

    q0 = q
    call on_shell (q0, [mp**2, mq**2], KEEP_MOMENTUM)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q0), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q0(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q0(2), u)
    write (u, "(A)")

    write (u, "(A)")  "Compare: mo^2"
    q02_2 = q0(2) ** 2;  call pacify (q02_2, 1e-10_default)
    write (u, "(2(1x,F11.8))")  sd%m2, q02_2
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momentum"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2)
    call sd%set_t_bounds (x, 1 - x)
    call sd%recover (k, q0, KEEP_MOMENTUM)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x
    write (u, "(A)")  "Compare: t"
    write (u, "(2(1x,F11.8))")  q2_2, sd%t

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_aux_3"

  end subroutine sf_aux_3

  subroutine sf_aux_4 (u)
    integer, intent(in) :: u
    type(splitting_data_t) :: sd
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E, mk, mp, mq, qmin, qmax
    real(default) :: x, xb

    write (u, "(A)")  "* Test output: sf_aux_4"
    write (u, "(A)")  "*   Purpose: compute massless collinear splitting near endpoint"

    E = 1
    mk = 0
    mp = 0
    mq = 0
    qmin = 1e-2_default
    qmax = 1e0_default

    k = vector4_moving (E, sqrt (E**2 - mk**2), 3)

    x = 0.1_default
    xb = 1 - x

    write (u, "(A)")
    write (u, "(A)")  "* (1) Collinear setup, moderate kinematics"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, xb)

    call sd%write (u)

    q = sd%split_momentum (k)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momenta"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, xb)
    call sd%recover (k, q, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  xb, sd%xb

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* (2) Close to x=0"
    write (u, "(A)")

    x = 1e-9_default
    xb = 1 - x

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, xb)

    call sd%write (u)

    q = sd%split_momentum (k)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momenta"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, xb)
    call sd%recover (k, q, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  xb, sd%xb

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* (3) Close to x=1"
    write (u, "(A)")

    xb = 1e-9_default
    x = 1 - xb

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, xb)

    call sd%write (u)

    q = sd%split_momentum (k)

    write (u, "(A)")
    write (u, "(A)")  "Incoming momentum k ="
    call vector4_write (k, u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum sum p + q ="
    call vector4_write (sum (q), u)
    write (u, "(A)")
    write (u, "(A)")  "Radiated momentum p ="
    call vector4_write (q(1), u)
    write (u, "(A)")
    write (u, "(A)")  "Outgoing momentum q ="
    call vector4_write (q(2), u)
    write (u, "(A)")

    write (u, "(A)")  "* Recover parameters from outgoing momenta"
    write (u, "(A)")

    call sd%init (k, mk**2, mp**2, mq**2, collinear = .true.)
    call sd%set_t_bounds (x, xb)
    call sd%recover (k, q, KEEP_ENERGY)

    write (u, "(A)")  "Compare: x"
    write (u, "(2(1x,F11.8))")  x, sd%x

    write (u, "(A)")  "Compare: 1-x"
    write (u, "(2(1x,F11.8))")  xb, sd%xb

    write (u, "(A)")
    call sd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_aux_4"

  end subroutine sf_aux_4


end module sf_aux_uti
