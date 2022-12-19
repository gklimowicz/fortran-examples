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

submodule (phs_fks) phs_fks_s

  use debug_master, only: debug_on
  use diagnostics
  use io_units, only: given_output_unit, free_unit
  use format_defs, only: FMT_17
  use format_utils, only: write_separator
  use physics_defs
  use flavors
  use pdg_arrays, only: is_colored
  use cascades
  use cascades2
  use ttv_formfactors, only: generate_on_shell_decay_threshold, m1s_to_mpole

  implicit none


  interface compute_beta
    module procedure compute_beta_massless
    module procedure compute_beta_massive
  end interface

  interface get_xi_max_fsr
    module procedure get_xi_max_fsr_massless
    module procedure get_xi_max_fsr_massive
  end interface


contains

  module subroutine isr_kinematics_write (isr, unit)
    class(isr_kinematics_t), intent(in) :: isr
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit); if (u < 0) return
    write (u,"(A)") "ISR kinematics: "
    write (u,"(A," // FMT_17 // ",1X)") "x(+):         ", isr%x(1)
    write (u,"(A," // FMT_17 // ",1X)") "x(-):         ", isr%x(2)
    write (u,"(A," // FMT_17 // ",1X)") "z(+):         ", isr%z(1)
    write (u,"(A," // FMT_17 // ",1X)") "z(-):         ", isr%z(2)
    write (u,"(A," // FMT_17 // ",1X)") "sqrts (Born): ", isr%sqrts_born
    if (allocated (isr%beam_energy)) then
       do i = 1, size (isr%beam_energy)
          write (u,"(A," // FMT_17 // ",1X)") "Beam energy:  ", &
               isr%beam_energy(i)
       end do
    end if
    write (u,"(A," // FMT_17 // ",1X)") "Fac. scale:   ", isr%fac_scale
    do i = 1, 2
       write (u,"(A," // FMT_17 // ",1X)") "Jacobian:     ", isr%jacobian(i)
    end do
    write (u,"(A,I0,1X)") "ISR mode:     ", isr%isr_mode
  end subroutine isr_kinematics_write

  module subroutine phs_point_set_init (phs_point_set, n_particles, n_phs)
    class(phs_point_set_t), intent(out) :: phs_point_set
    integer, intent(in) :: n_particles, n_phs
    integer :: i_phs
    allocate (phs_point_set%phs_point (n_phs))
    do i_phs = 1, n_phs
       phs_point_set%phs_point(i_phs) = n_particles
    end do
    phs_point_set%initialized = .true.
  end subroutine phs_point_set_init

  module subroutine phs_point_set_write (phs_point_set, i_phs, contributors, &
       unit, show_mass, testflag, check_conservation, ultra, n_in)
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in), optional :: i_phs
    integer, intent(in), dimension(:), optional :: contributors
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_mass
    logical, intent(in), optional :: testflag, ultra
    logical, intent(in), optional :: check_conservation
    integer, intent(in), optional :: n_in
    integer :: i, u
    type(vector4_t) :: p_sum
    u = given_output_unit (unit); if (u < 0) return
    if (present (i_phs)) then
       call phs_point_set%phs_point(i_phs)%write &
            (unit = u, show_mass = show_mass, testflag = testflag, &
            check_conservation = check_conservation, ultra = ultra, n_in = n_in)
    else
       do i = 1, size(phs_point_set%phs_point)
          call phs_point_set%phs_point(i)%write &
               (unit = u, show_mass = show_mass, testflag = testflag, &
               check_conservation = check_conservation, ultra = ultra,&
               n_in = n_in)
       end do
    end if
    if (present (contributors)) then
       if (debug_on)  call msg_debug &
            (D_SUBTRACTION, "Invariant masses for real emission: ")
       associate (pp => phs_point_set%phs_point(i_phs))
         p_sum = sum (pp, [contributors, size (pp)])
       end associate
       if (debug_active (D_SUBTRACTION)) &
            call vector4_write (p_sum, unit = unit, show_mass = show_mass, &
            testflag = testflag, ultra = ultra)
    end if
  end subroutine phs_point_set_write

  elemental module function phs_point_set_get_n_momenta &
       (phs_point_set, i_res) result (n)
    integer :: n
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_res
    n = size (phs_point_set%phs_point(i_res))
  end function phs_point_set_get_n_momenta

  pure module function phs_point_set_get_momenta &
       (phs_point_set, i_phs, n_in) result (p)
    type(vector4_t), dimension(:), allocatable :: p
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_phs
    integer, intent(in), optional :: n_in
    integer :: i
    if (present (n_in)) then
       p = phs_point_set%phs_point(i_phs)%select ([(i, i=1, n_in)])
    else
       p = phs_point_set%phs_point(i_phs)
    end if
  end function phs_point_set_get_momenta

  pure module function phs_point_set_get_momentum &
       (phs_point_set, i_phs, i_mom) result (p)
    type(vector4_t) :: p
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_phs, i_mom
    p = phs_point_set%phs_point(i_phs)%select (i_mom)
  end function phs_point_set_get_momentum

  pure module function phs_point_set_get_energy &
       (phs_point_set, i_phs, i_mom) result (E)
    real(default) :: E
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_phs, i_mom
    E = energy (phs_point_set%phs_point(i_phs)%select (i_mom))
  end function phs_point_set_get_energy

  module function phs_point_set_get_sqrts &
       (phs_point_set, i_phs) result (sqrts)
    real(default) :: sqrts
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_phs
    sqrts = sqrt (phs_point_set%phs_point(i_phs)%get_msq ([1,2]))
  end function phs_point_set_get_sqrts

  module subroutine phs_point_set_set_momenta_p (phs_point_set, i_phs, p)
    class(phs_point_set_t), intent(inout) :: phs_point_set
    integer, intent(in) :: i_phs
    type(vector4_t), intent(in), dimension(:) :: p
    phs_point_set%phs_point(i_phs) = p
  end subroutine phs_point_set_set_momenta_p

  module subroutine phs_point_set_set_momenta_phs_point &
       (phs_point_set, i_phs, p)
    class(phs_point_set_t), intent(inout) :: phs_point_set
    integer, intent(in) :: i_phs
    type(phs_point_t), intent(in) :: p
    phs_point_set%phs_point(i_phs) = p
  end subroutine phs_point_set_set_momenta_phs_point

  module function phs_point_set_get_n_particles &
       (phs_point_set, i) result (n_particles)
    integer :: n_particles
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in), optional :: i
    integer :: j
    j = 1; if (present (i)) j = i
    n_particles = size (phs_point_set%phs_point(j))
  end function phs_point_set_get_n_particles

  module function phs_point_set_get_n_phs (phs_point_set) result (n_phs)
    integer :: n_phs
    class(phs_point_set_t), intent(in) :: phs_point_set
    n_phs = size (phs_point_set%phs_point)
  end function phs_point_set_get_n_phs

  module function phs_point_set_get_invariant_mass &
       (phs_point_set, i_phs, i_part) result (m2)
    real(default) :: m2
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_phs
    integer, intent(in), dimension(:) :: i_part
    m2 = phs_point_set%phs_point(i_phs)%get_msq (i_part)
  end function phs_point_set_get_invariant_mass

  module subroutine phs_point_set_write_phs_point (phs_point_set, i_phs, &
       unit, show_mass, testflag, check_conservation, ultra, n_in)
    class(phs_point_set_t), intent(in) :: phs_point_set
    integer, intent(in) :: i_phs
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: show_mass
    logical, intent(in), optional :: testflag, ultra
    logical, intent(in), optional :: check_conservation
    integer, intent(in), optional :: n_in
    call phs_point_set%phs_point(i_phs)%write (unit, show_mass, testflag, &
       check_conservation, ultra, n_in)
  end subroutine phs_point_set_write_phs_point

  module subroutine phs_point_set_final (phs_point_set)
    class(phs_point_set_t), intent(inout) :: phs_point_set
    integer :: i
    deallocate (phs_point_set%phs_point)
    phs_point_set%initialized = .false.
  end subroutine phs_point_set_final

  module subroutine real_kinematics_init (r, n_tot, n_phs, n_alr, n_contr)
    class(real_kinematics_t), intent(inout) :: r
    integer, intent(in) :: n_tot, n_phs, n_alr, n_contr
    allocate (r%xi_max (n_phs))
    allocate (r%y (n_phs))
    allocate (r%y_soft (n_phs))
    call r%p_born_cms%init (n_tot - 1, 1)
    call r%p_born_lab%init (n_tot - 1, 1)
    call r%p_real_cms%init (n_tot, n_phs)
    call r%p_real_lab%init (n_tot, n_phs)
    allocate (r%jac (n_phs), r%jac_rand (n_phs))
    allocate (r%alr_to_i_phs (n_alr))
    allocate (r%xi_ref_momenta (n_contr))
    r%alr_to_i_phs = 0
    r%xi_tilde = zero; r%xi_mismatch = zero
    r%xi_max = zero
    r%y = zero; r%y_mismatch = zero
    r%y_soft = zero
    r%phi = zero
    r%cms_energy2 = zero
    r%xi_ref_momenta = vector4_null
    r%jac_mismatch = one
    r%jac_rand = one
  end subroutine real_kinematics_init

  module subroutine real_kinematics_init_onshell (r, n_tot, n_phs)
    class(real_kinematics_t), intent(inout) :: r
    integer, intent(in) :: n_tot, n_phs
    call r%p_born_onshell%init (n_tot - 1, 1)
    call r%p_real_onshell(1)%init (n_tot, n_phs)
    call r%p_real_onshell(2)%init (n_tot, n_phs)
  end subroutine real_kinematics_init_onshell

  module subroutine real_kinematics_write (r, unit)
    class(real_kinematics_t), intent(in) :: r
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit); if (u < 0) return
    write (u,"(A)") "Real kinematics: "
    write (u,"(A," // FMT_17 // ",1X)") "xi_tilde: ", r%xi_tilde
    write (u,"(A," // FMT_17 // ",1X)") "phi: ", r%phi
    do i = 1, size (r%xi_max)
       write (u,"(A,I1,1X)") "i_phs: ", i
       write (u,"(A," // FMT_17 // ",1X)") "xi_max: ", r%xi_max(i)
       write (u,"(A," // FMT_17 // ",1X)") "y: ", r%y(i)
       write (u,"(A," // FMT_17 // ",1X)") "jac_rand: ", r%jac_rand(i)
       write (u,"(A," // FMT_17 // ",1X)") "y_soft: ", r%y_soft(i)
    end do
    write (u, "(A)") "Born Momenta: "
    write (u, "(A)") "CMS: "
    call r%p_born_cms%write (unit = u)
    write (u, "(A)") "Lab: "
    call r%p_born_lab%write (unit = u)
    write (u, "(A)") "Real Momenta: "
    write (u, "(A)") "CMS: "
    call r%p_real_cms%write (unit = u)
    write (u, "(A)") "Lab: "
    call r%p_real_lab%write (unit = u)
  end subroutine real_kinematics_write

  module function get_boost_for_threshold_projection &
       (p, sqrts, mtop) result (L)
    type(lorentz_transformation_t) :: L
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(in) :: sqrts, mtop
    type(vector4_t) :: p_tmp
    type(vector3_t) :: dir
    real(default) :: scale_factor, arg
    p_tmp = p(THR_POS_WP) + p(THR_POS_B)
    arg = sqrts**2 - four * mtop**2
    if (arg > zero) then
       scale_factor = sqrt (arg) / two
    else
       scale_factor = tiny_07*1000
    end if
    dir = scale_factor * create_unit_vector (p_tmp)
    p_tmp = [sqrts / two, dir%p]
    L = boost (p_tmp, mtop)
  end function get_boost_for_threshold_projection

  function get_generation_phi (p_born, p_real, emitter, i_gluon) result (phi)
    real(default) :: phi
    type(vector4_t), intent(in), dimension(:) :: p_born, p_real
    integer, intent(in) :: emitter, i_gluon
    type(vector4_t) :: p1, p2, pp
    type(lorentz_transformation_t) :: rot_to_gluon, rot_to_z
    type(vector3_t) :: dir, z
    real(default) :: cpsi
    pp = p_real(emitter) + p_real(i_gluon)
    cpsi = (space_part_norm (pp)**2 - space_part_norm (p_real(emitter))**2 &
           + space_part_norm (p_real(i_gluon))**2) / &
           (two * space_part_norm (pp) * space_part_norm (p_real(i_gluon)))
    dir = create_orthogonal (space_part (p_born(emitter)))
    rot_to_gluon = rotation (cpsi, sqrt (one - cpsi**2), dir)
    pp = rot_to_gluon * p_born(emitter)
    z%p = [0._default, 0._default, 1._default]
    rot_to_z = rotation_to_2nd &
         (space_part (p_born(emitter)) / space_part_norm (p_born(emitter)), z)
    p1 = rot_to_z * pp / space_part_norm (pp)
    p2 = rot_to_z * p_real(i_gluon)
    phi = azimuthal_distance (p1, p2)
    if (phi < zero) phi = twopi - abs(phi)
  end function get_generation_phi

  module subroutine real_kinematics_apply_threshold_projection_real &
       (r, i_phs, mtop, L_to_cms, invert)
    class(real_kinematics_t), intent(inout) :: r
    integer, intent(in) :: i_phs
    real(default), intent(in) :: mtop
    type(lorentz_transformation_t), intent(in), dimension(:) :: L_to_cms
    logical, intent(in) :: invert
    integer :: leg, other_leg
    type(vector4_t), dimension(:), allocatable :: p_real
    type(vector4_t), dimension(:), allocatable :: p_real_onshell
    type(vector4_t), dimension(4) :: k_tmp
    type(vector4_t), dimension(4) :: k_decay_onshell_real
    type(vector4_t), dimension(3) :: k_decay_onshell_born
    do leg = 1, 2
       other_leg = 3 - leg
       p_real = r%p_real_cms%phs_point(i_phs)
       allocate (p_real_onshell (size (p_real)))
       p_real_onshell(1:2) = p_real(1:2)
       k_tmp(1) = p_real(7)
       k_tmp(2) = p_real(ass_quark(leg))
       k_tmp(3) = p_real(ass_boson(leg))
       k_tmp(4) = [mtop, zero, zero, zero]
       call generate_on_shell_decay_threshold (k_tmp(1:3), &
            k_tmp(4), k_decay_onshell_real (2:4))
       k_decay_onshell_real (1) = k_tmp(4)
       k_tmp(1) = p_real(ass_quark(other_leg))
       k_tmp(2) = p_real(ass_boson(other_leg))
       k_decay_onshell_born = create_two_particle_decay (mtop**2, k_tmp(1), k_tmp(2))
       p_real_onshell(THR_POS_GLUON) = L_to_cms(leg) * k_decay_onshell_real (2)
       p_real_onshell(ass_quark(leg)) = L_to_cms(leg) * k_decay_onshell_real(3)
       p_real_onshell(ass_boson(leg)) = L_to_cms(leg) * k_decay_onshell_real(4)
       p_real_onshell(ass_quark(other_leg)) = L_to_cms(leg) * k_decay_onshell_born (2)
       p_real_onshell(ass_boson(other_leg)) = L_to_cms(leg) * k_decay_onshell_born (3)
       if (invert) then
          call vector4_invert_direction (p_real_onshell (ass_quark(other_leg)))
          call vector4_invert_direction (p_real_onshell (ass_boson(other_leg)))
       end if
       r%p_real_onshell(leg)%phs_point(i_phs) = p_real_onshell
       deallocate (p_real_onshell)
    end do
  end subroutine real_kinematics_apply_threshold_projection_real

  module subroutine threshold_projection_born &
       (mtop, L_to_cms, p_in, p_onshell)
    real(default), intent(in) :: mtop
    type(lorentz_transformation_t), intent(in) :: L_to_cms
    type(vector4_t), intent(in), dimension(:) :: p_in
    type(vector4_t), intent(out), dimension(:) :: p_onshell
    type(vector4_t), dimension(3) :: k_decay_onshell
    type(vector4_t) :: p_tmp_1, p_tmp_2
    type(lorentz_transformation_t) :: L_to_cms_inv
    p_onshell(1:2) = p_in(1:2)
    L_to_cms_inv = inverse (L_to_cms)
    p_tmp_1 = L_to_cms_inv * p_in(THR_POS_B)
    p_tmp_2 = L_to_cms_inv * p_in(THR_POS_WP)
    k_decay_onshell = create_two_particle_decay (mtop**2, &
         p_tmp_1, p_tmp_2)
    p_onshell([THR_POS_B, THR_POS_WP]) = k_decay_onshell([2, 3])
    p_tmp_1 = L_to_cms * p_in(THR_POS_BBAR)
    p_tmp_2 = L_to_cms * p_in(THR_POS_WM)
    k_decay_onshell = create_two_particle_decay (mtop**2, &
         p_tmp_1, p_tmp_2)
    p_onshell([THR_POS_BBAR, THR_POS_WM]) = k_decay_onshell([2, 3])
    p_onshell([THR_POS_WP, THR_POS_B]) = L_to_cms * p_onshell([THR_POS_WP, THR_POS_B])
    p_onshell([THR_POS_WM, THR_POS_BBAR]) = L_to_cms_inv * p_onshell([THR_POS_WM, THR_POS_BBAR])
  end subroutine threshold_projection_born

  pure module subroutine compute_dalitz_bounds &
       (q0, m2, mrec2, z1, z2, k0_rec_max)
    real(default), intent(in) :: q0, m2, mrec2
    real(default), intent(out) :: z1, z2, k0_rec_max
    k0_rec_max = (q0**2 - m2 + mrec2) / (two * q0)
    z1 = (k0_rec_max + sqrt(k0_rec_max**2 - mrec2)) / q0
    z2 = (k0_rec_max - sqrt(k0_rec_max**2 - mrec2)) / q0
  end subroutine compute_dalitz_bounds

  module function real_kinematics_kt2 &
       (real_kinematics, i_phs, emitter, kt2_type, xi, y) result (kt2)
    real(default) :: kt2
    class(real_kinematics_t), intent(in) :: real_kinematics
    integer, intent(in) :: emitter, i_phs, kt2_type
    real(default), intent(in), optional :: xi, y
    real(default) :: xii, yy
    real(default) :: q, E_em, z, z1, z2, m2, mrec2, k0_rec_max
    type(vector4_t) :: p_emitter
    if (present (y)) then
       yy = y
    else
       yy = real_kinematics%y (i_phs)
    end if
    if (present (xi)) then
       xii = xi
    else
       xii = real_kinematics%xi_tilde * real_kinematics%xi_max (i_phs)
    end if
    select case (kt2_type)
    case (UBF_FSR_SIMPLE)
       kt2 = real_kinematics%cms_energy2 / two * xii**2 * (1 - yy)
    case (UBF_FSR_MASSIVE)
       q = sqrt (real_kinematics%cms_energy2)
       p_emitter = real_kinematics%p_born_cms%phs_point(1)%select (emitter)
       mrec2 = (q - p_emitter%p(0))**2 - sum (p_emitter%p(1:3)**2)
       m2 = p_emitter**2
       E_em = energy (p_emitter)
       call compute_dalitz_bounds (q, m2, mrec2, z1, z2, k0_rec_max)
       z = z2 - (z2 - z1) * (one + yy) / two
       kt2 = xii**2 * q**3 * (one - z) / &
          (two * E_em - z * xii * q)
    case (UBF_FSR_MASSLESS_RECOIL)
       kt2 = real_kinematics%cms_energy2 / two * xii**2 * (1 - yy**2) / two
    case (UBF_ISR)
       kt2 = real_kinematics%cms_energy2 / four * xii**2 * (1 - yy**2) / (1 - xii)
    case default
       kt2 = zero
       call msg_bug ("kt2_type must be set to a known value")
    end select
  end function real_kinematics_kt2

  module subroutine real_kinematics_final (real_kin)
    class(real_kinematics_t), intent(inout) :: real_kin
    if (allocated (real_kin%xi_max)) deallocate (real_kin%xi_max)
    if (allocated (real_kin%y)) deallocate (real_kin%y)
    if (allocated (real_kin%alr_to_i_phs)) deallocate (real_kin%alr_to_i_phs)
    if (allocated (real_kin%jac_rand)) deallocate (real_kin%jac_rand)
    if (allocated (real_kin%y_soft)) deallocate (real_kin%y_soft)
    if (allocated (real_kin%xi_ref_momenta))  &
         deallocate (real_kin%xi_ref_momenta)
    call real_kin%p_born_cms%final (); call real_kin%p_born_lab%final ()
    call real_kin%p_real_cms%final (); call real_kin%p_real_lab%final ()
  end subroutine real_kinematics_final

  module subroutine fks_config_clear_phase_space (phs_config)
    class(phs_fks_config_t), intent(inout) :: phs_config
  end subroutine fks_config_clear_phase_space

  module subroutine phs_fks_config_write (object, unit, include_id)
    class(phs_fks_config_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: include_id
    integer :: u
    u = given_output_unit (unit)
    call object%phs_wood_config_t%write (u)
    write (u, "(3x,A,I0)")  "NLO mode  = ", object%mode
    write (u, "(3x,A,L1)")  "2->1 proc = ", object%born_2_to_1
    write (u, "(3x,A,L1)")  "Dalitz    = ", object%make_dalitz_plot
    write (u, "(A,A)") "Extra Born md5sum: ", object%md5sum_born_config
  end subroutine phs_fks_config_write

  module subroutine phs_fks_config_set_mode (phs_config, mode)
    class(phs_fks_config_t), intent(inout) :: phs_config
    integer, intent(in) :: mode
    select case (mode)
    case (NLO_REAL, NLO_MISMATCH)
       phs_config%mode = PHS_MODE_ADDITIONAL_PARTICLE
    case (NLO_DGLAP)
       phs_config%mode = PHS_MODE_COLLINEAR_REMNANT
    end select
  end subroutine phs_fks_config_set_mode

  module subroutine phs_fks_config_configure (phs_config, sqrts, &
       sqrts_fixed, lab_is_cm, azimuthal_dependence, rebuild, &
       ignore_mismatch, nlo_type, subdir)
    class(phs_fks_config_t), intent(inout) :: phs_config
    real(default), intent(in) :: sqrts
    logical, intent(in), optional :: sqrts_fixed
    logical, intent(in), optional :: lab_is_cm
    logical, intent(in), optional :: azimuthal_dependence
    logical, intent(in), optional :: rebuild
    logical, intent(in), optional :: ignore_mismatch
    integer, intent(in), optional :: nlo_type
    type(string_t), intent(in), optional :: subdir
    if (present (nlo_type))  phs_config%nlo_type = nlo_type
    if (.not. phs_config%is_combined_integration) then
       select case (phs_config%mode)
       case (PHS_MODE_ADDITIONAL_PARTICLE)
          phs_config%n_par = phs_config%n_par + 3
          if (phs_config%nlo_type == NLO_REAL .and. phs_config%n_out == 2) then
             phs_config%born_2_to_1 = .true.
          end if
       case (PHS_MODE_COLLINEAR_REMNANT)
          phs_config%n_par = phs_config%n_par + 1
       end select
    end if
    call phs_config%compute_md5sum ()
  end subroutine phs_fks_config_configure

  module subroutine phs_fks_config_startup_message (phs_config, unit)
    class(phs_fks_config_t), intent(in) :: phs_config
    integer, intent(in), optional :: unit
    call phs_config%phs_wood_config_t%startup_message (unit)
  end subroutine phs_fks_config_startup_message

  module subroutine phs_fks_config_generate_phase_space_extra (phs_config)
    class(phs_fks_config_t), intent(inout) :: phs_config
    integer :: off_shell, extra_off_shell
    type(flavor_t), dimension(:,:), allocatable :: flv_born
    integer :: i, j
    integer :: n_state, n_flv_born
    integer :: unit_fds
    logical :: valid
    type(string_t) :: file_name
    logical :: file_exists
    if (phs_config%use_cascades2) then
       allocate (phs_config%feyngraph_set)
    else
       allocate (phs_config%cascade_set)
    end if
    n_flv_born = size (phs_config%flv, 1) - 1
    n_state = size (phs_config%flv, 2)
    allocate (flv_born (n_flv_born, n_state))
    do i = 1, n_flv_born
       do j = 1, n_state
          flv_born(i, j) = phs_config%flv(i, j)
       end do
    end do
    if (phs_config%use_cascades2) then
       file_name = char (phs_config%id) // ".fds"
       inquire (file=char (file_name), exist=file_exists)
       if (.not. file_exists) call msg_fatal &
            ("The O'Mega input file " // char (file_name) // &
            " does not exist. " // "Please make sure that the " // &
            "variable ?omega_write_phs_output has been set correctly.")
       unit_fds = free_unit ()
       open (unit=unit_fds, file=char(file_name), status='old', action='read')
    end if
    off_shell = phs_config%par%off_shell
    do extra_off_shell = 0, max (n_flv_born - 2, 0)
       phs_config%par%off_shell = off_shell + extra_off_shell
       if (phs_config%use_cascades2) then
          call feyngraph_set_generate (phs_config%feyngraph_set, &
               phs_config%model, phs_config%n_in, phs_config%n_out - 1, &
               flv_born, phs_config%par, phs_config%fatal_beam_decay, unit_fds, &
               phs_config%vis_channels)
          if (feyngraph_set_is_valid (phs_config%feyngraph_set)) exit
       else
          call cascade_set_generate (phs_config%cascade_set, &
               phs_config%model, phs_config%n_in, phs_config%n_out - 1, &
               flv_born, phs_config%par, phs_config%fatal_beam_decay)
          if (cascade_set_is_valid (phs_config%cascade_set)) exit
       end if
    end do
    if (phs_config%use_cascades2) then
       close (unit_fds)
       valid = feyngraph_set_is_valid (phs_config%feyngraph_set)
    else
       valid = cascade_set_is_valid (phs_config%cascade_set)
    end if
    if (.not. valid) &
       call msg_fatal ("Resonance extraction: Phase space generation failed")
  end subroutine phs_fks_config_generate_phase_space_extra

  module subroutine phs_fks_config_set_born_config (phs_config, phs_cfg_born)
    class(phs_fks_config_t), intent(inout) :: phs_config
    type(phs_wood_config_t), intent(in), target :: phs_cfg_born
    if (debug_on) &
         call msg_debug (D_PHASESPACE, "phs_fks_config_set_born_config")
    phs_config%forest = phs_cfg_born%forest
    phs_config%n_channel = phs_cfg_born%n_channel
    allocate (phs_config%channel (phs_config%n_channel))
    phs_config%channel = phs_cfg_born%channel
    phs_config%n_par = phs_cfg_born%n_par
    phs_config%n_state = phs_cfg_born%n_state
    phs_config%sqrts = phs_cfg_born%sqrts
    phs_config%par = phs_cfg_born%par
    phs_config%sqrts_fixed = phs_cfg_born%sqrts_fixed
    phs_config%azimuthal_dependence = phs_cfg_born%azimuthal_dependence
    phs_config%provides_chains = phs_cfg_born%provides_chains
    phs_config%lab_is_cm = phs_cfg_born%lab_is_cm
    phs_config%vis_channels = phs_cfg_born%vis_channels
    phs_config%provides_equivalences = phs_cfg_born%provides_equivalences
    allocate (phs_config%chain (size (phs_cfg_born%chain)))
    phs_config%chain = phs_cfg_born%chain
    phs_config%model => phs_cfg_born%model
    phs_config%use_cascades2 = phs_cfg_born%use_cascades2
    if (allocated (phs_cfg_born%cascade_set)) then
       allocate (phs_config%cascade_set)
       phs_config%cascade_set = phs_cfg_born%cascade_set
    end if
    if (allocated (phs_cfg_born%feyngraph_set)) then
       allocate (phs_config%feyngraph_set)
       phs_config%feyngraph_set = phs_cfg_born%feyngraph_set
    end if
    phs_config%md5sum_born_config = phs_cfg_born%md5sum_phs_config
  end subroutine phs_fks_config_set_born_config

  module function phs_fks_config_get_resonance_histories &
       (phs_config) result (resonance_histories)
    type(resonance_history_t), dimension(:), allocatable :: resonance_histories
    class(phs_fks_config_t), intent(inout) :: phs_config
    if (allocated (phs_config%cascade_set)) then
       call cascade_set_get_resonance_histories (phs_config%cascade_set, &
            n_filter = 2, res_hists = resonance_histories)
    else if (allocated (phs_config%feyngraph_set)) then
       call feyngraph_set_get_resonance_histories (phs_config%feyngraph_set, &
            n_filter = 2, res_hists = resonance_histories)
    else
       if (debug_on) call msg_debug (D_PHASESPACE, "Have to rebuild phase space for resonance histories")
       call phs_config%generate_phase_space_extra ()
       if (phs_config%use_cascades2) then
          call feyngraph_set_get_resonance_histories &
               (phs_config%feyngraph_set, n_filter = 2, res_hists = resonance_histories)
       else
          call cascade_set_get_resonance_histories &
               (phs_config%cascade_set, n_filter = 2, res_hists = resonance_histories)
       end if
    end if
  end function phs_fks_config_get_resonance_histories

  module subroutine dalitz_plot_init (plot, unit, filename, inverse)
    class(dalitz_plot_t), intent(inout) :: plot
    integer, intent(in) :: unit
    type(string_t), intent(in) :: filename
    logical, intent(in) :: inverse
    plot%active = .true.
    plot%unit = unit
    plot%inverse = inverse
    open (plot%unit, file = char (filename), action = "write")
  end subroutine dalitz_plot_init

  module subroutine dalitz_plot_write_header (plot)
    class(dalitz_plot_t), intent(in) :: plot
    write (plot%unit, "(A36)") "### Dalitz plot generated by WHIZARD"
    if (plot%inverse) then
       write (plot%unit, "(A10,1x,A4)") "### k0_n+1", "k0_n"
    else
       write (plot%unit, "(A8,1x,A6)") "### k0_n", "k0_n+1"
    end if
  end subroutine dalitz_plot_write_header

  module subroutine dalitz_plot_register (plot, k0_n, k0_np1)
    class(dalitz_plot_t), intent(in) :: plot
    real(default), intent(in) :: k0_n, k0_np1
    if (plot%inverse) then
       write (plot%unit, "(F8.4,1X,F8.4)") k0_np1, k0_n
    else
       write (plot%unit, "(F8.4,1X,F8.4)") k0_np1, k0_n
    end if
  end subroutine dalitz_plot_register

  module subroutine dalitz_plot_final (plot)
    class(dalitz_plot_t), intent(inout) :: plot
    logical :: opened
    plot%active = .false.
    plot%inverse = .false.
    if (plot%unit >= 0) then
       inquire (unit = plot%unit, opened = opened)
       if (opened) close (plot%unit)
    end if
    plot%filename = var_str ('')
    plot%unit = -1
  end subroutine dalitz_plot_final

  module function check_scalar_products (p) result (valid)
    logical :: valid
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), parameter :: tolerance = 1E-7_default
    !!! for testing EW singularities: tolerance = 5E-9_default
    integer :: i, j
    valid = .true.
    do i = 1, size (p)
       do j = i, size (p)
          if (i /= j) then
             if (abs(p(i) * p(j)) < tolerance) then
                valid = .false.
                exit
             end if
          end if
       end do
    end do
  end function check_scalar_products

  module subroutine phs_fks_generator_connect_kinematics &
       (generator, isr_kinematics, real_kinematics, massive_phsp)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(isr_kinematics_t), intent(in), pointer :: isr_kinematics
    type(real_kinematics_t), intent(in), pointer :: real_kinematics
    logical, intent(in) :: massive_phsp
    generator%real_kinematics => real_kinematics
    generator%isr_kinematics => isr_kinematics
    generator%massive_phsp = massive_phsp
  end subroutine phs_fks_generator_connect_kinematics

  module subroutine phs_fks_generator_compute_isr_kinematics &
       (generator, r, p_in)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: r
    type(vector4_t), dimension(2), intent(in), optional :: p_in
    integer :: em
    type(vector4_t), dimension(2) :: p

    if (present (p_in)) then
       p = p_in
    else
       p = generator%real_kinematics%p_born_lab%phs_point(1)%select ([1,2])
    end if

    associate (isr_kinematics => generator%isr_kinematics)
       do em = 1, 2
          isr_kinematics%x(em) = p(em)%p(0) / isr_kinematics%beam_energy(em)
          isr_kinematics%z(em) = one - (one - isr_kinematics%x(em)) * r
          isr_kinematics%jacobian(em) =  one - isr_kinematics%x(em)
       end do
       isr_kinematics%sqrts_born = (p(1) + p(2))**1
    end associate
  end subroutine phs_fks_generator_compute_isr_kinematics

  module subroutine phs_fks_generator_final (generator)
    class(phs_fks_generator_t), intent(inout) :: generator
    if (allocated (generator%emitters)) deallocate (generator%emitters)
    if (associated (generator%real_kinematics)) &
         nullify (generator%real_kinematics)
    if (associated (generator%isr_kinematics)) &
         nullify (generator%isr_kinematics)
    if (allocated (generator%m2)) deallocate (generator%m2)
    generator%massive_phsp = .false.
    if (allocated (generator%is_massive)) deallocate (generator%is_massive)
    generator%singular_jacobian = .false.
    generator%i_fsr_first = -1
    if (allocated (generator%resonance_contributors)) &
           deallocate (generator%resonance_contributors)
    generator%mode = GEN_REAL_PHASE_SPACE
  end subroutine phs_fks_generator_final

  module subroutine phs_identifier_init_from_emitter (phs_id, emitter)
    class(phs_identifier_t), intent(out) :: phs_id
    integer, intent(in) :: emitter
    phs_id%emitter = emitter
  end subroutine phs_identifier_init_from_emitter

  module subroutine phs_identifier_init_from_emitter_and_contributors &
     (phs_id, emitter, contributors)
     class(phs_identifier_t), intent(out) :: phs_id
     integer, intent(in) :: emitter
     integer, intent(in), dimension(:) :: contributors
     allocate (phs_id%contributors (size (contributors)))
     phs_id%contributors = contributors
     phs_id%emitter = emitter
  end subroutine phs_identifier_init_from_emitter_and_contributors

  module function phs_identifier_check &
       (phs_id, emitter, contributors) result (check)
    logical :: check
    class(phs_identifier_t), intent(in) :: phs_id
    integer, intent(in) :: emitter
    integer, intent(in), dimension(:), optional :: contributors
    check = phs_id%emitter == emitter
    if (present (contributors)) then
       if (.not. allocated (phs_id%contributors)) &
          call msg_fatal ("Phs identifier: contributors not allocated!")
       check = check .and. all (phs_id%contributors == contributors)
    end if
  end function phs_identifier_check

  module subroutine phs_identifier_write (phs_id, unit)
    class(phs_identifier_t), intent(in) :: phs_id
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit); if (u < 0) return
    write (u, '(A)') 'phs_identifier: '
    write (u, '(A,1X,I1)') 'Emitter: ', phs_id%emitter
    if (allocated (phs_id%contributors)) then
       write (u, '(A)', advance = 'no') 'Resonance contributors: '
       do i = 1, size (phs_id%contributors)
          write (u, '(I1,1X)', advance = 'no') phs_id%contributors(i)
       end do
    else
       write (u, '(A)') 'No Contributors allocated'
    end if
  end subroutine phs_identifier_write

  module subroutine check_for_phs_identifier &
       (phs_id, n_in, emitter, contributors, phs_exist, i_phs)
    type(phs_identifier_t), intent(in), dimension(:) :: phs_id
    integer, intent(in) :: n_in, emitter
    integer, intent(in), dimension(:), optional :: contributors
    logical, intent(out) :: phs_exist
    integer, intent(out) :: i_phs
    integer :: i
    phs_exist = .false.
    i_phs = -1
    do i = 1, size (phs_id)
       if (phs_id(i)%emitter < 0) then
          i_phs = i
          exit
       end if
       phs_exist = phs_id(i)%emitter == emitter
       if (present (contributors)) &
            phs_exist = phs_exist .and. &
            all (phs_id(i)%contributors == contributors)
       if (phs_exist) then
          i_phs = i
          exit
       end if
    end do
  end subroutine check_for_phs_identifier

  module subroutine phs_fks_write (object, unit, verbose)
    class(phs_fks_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u, i, n_id
    u = given_output_unit (unit)
    call object%base_write ()
    n_id = size (object%phs_identifiers)
    if (n_id == 0) then
       write (u, "(A)") "No phs identifiers allocated! "
    else
       do i = 1, n_id
          call object%phs_identifiers(i)%write (u)
       end do
    end if
  end subroutine phs_fks_write

  module subroutine phs_fks_init (phs, phs_config)
    class(phs_fks_t), intent(out) :: phs
    class(phs_config_t), intent(in), target :: phs_config

    call phs%base_init (phs_config)
    select type (phs_config)
    type is (phs_fks_config_t)
       phs%config => phs_config
       phs%forest = phs_config%forest
    end select

    select type (phs)
    type is (phs_fks_t)
      select type (phs_config)
      type is (phs_fks_config_t)
         phs%mode = phs_config%mode
      end select

      select case (phs%mode)
      case (PHS_MODE_ADDITIONAL_PARTICLE)
         phs%n_r_born = phs%config%n_par - 3
      case (PHS_MODE_COLLINEAR_REMNANT)
         phs%n_r_born = phs%config%n_par - 1
      end select
    end select
  end subroutine phs_fks_init

  module subroutine phs_fks_compute_flux (phs)
    class(phs_fks_t), intent(inout) :: phs
    call phs%compute_base_flux ()
    select type (config => phs%config)
    type is (phs_fks_config_t)
       if (config%born_2_to_1) then
          phs%flux = conv * twopi &
               / (2 * config%sqrts ** 2 * phs%m_out(1) ** 2)
       end if
    end select
  end subroutine phs_fks_compute_flux

  module subroutine phs_fks_allocate_momenta (phs, phs_config, data_is_born)
    class(phs_fks_t), intent(inout) :: phs
    class(phs_config_t), intent(in) :: phs_config
    logical, intent(in) :: data_is_born
    integer :: n_out_born
    allocate (phs%p_born (phs_config%n_in))
    allocate (phs%p_real (phs_config%n_in))
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       if (data_is_born) then
          n_out_born = phs_config%n_out
       else
          n_out_born = phs_config%n_out - 1
       end if
       allocate (phs%q_born (n_out_born))
       allocate (phs%q_real (n_out_born + 1))
       allocate (phs%p_born_tot (phs_config%n_in + n_out_born))
    end select
  end subroutine phs_fks_allocate_momenta

  module subroutine phs_fks_evaluate_selected_channel (phs, c_in, r_in)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    real(default), intent(in), dimension(:) :: r_in
    integer :: n_in

    call phs%phs_wood_t%evaluate_selected_channel (c_in, r_in)
    phs%r(:,c_in) = r_in

    phs%q_defined = phs%phs_wood_t%q_defined
    if (.not. phs%q_defined) return

    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       n_in = phs%config%n_in
       phs%p_born = phs%phs_wood_t%p
       phs%q_born = phs%phs_wood_t%q
       phs%p_born_tot (1: n_in) = phs%p_born
       phs%p_born_tot (n_in + 1 :) = phs%q_born
       call phs%set_reference_frames (.true.)
       call phs%set_isr_kinematics (.true.)
    case (PHS_MODE_COLLINEAR_REMNANT)
       call phs%compute_isr_kinematics (r_in(phs%n_r_born + 1))
       phs%r_isr = r_in(phs%n_r_born + 1)
    end select
  end subroutine phs_fks_evaluate_selected_channel

  module subroutine phs_fks_evaluate_other_channels (phs, c_in)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: c_in
    call phs%phs_wood_t%evaluate_other_channels (c_in)
    phs%r_defined = .true.
  end subroutine phs_fks_evaluate_other_channels

  module subroutine phs_fks_get_mcpar (phs, c, r)
    class(phs_fks_t), intent(in) :: phs
    integer, intent(in) :: c
    real(default), dimension(:), intent(out) :: r
    r(1 : phs%n_r_born) = phs%r(1 : phs%n_r_born,c)
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       r(phs%n_r_born + 1 :) = phs%r_real
    case (PHS_MODE_COLLINEAR_REMNANT)
       r(phs%n_r_born + 1 :) = phs%r_isr
    end select
  end subroutine phs_fks_get_mcpar

  module subroutine phs_fks_set_beam_energy (phs)
    class(phs_fks_t), intent(inout) :: phs
    call phs%generator%set_sqrts_hat (phs%config%sqrts)
  end subroutine phs_fks_set_beam_energy

  module subroutine phs_fks_set_emitters (phs, emitters)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in), dimension(:), allocatable :: emitters
    call phs%generator%set_emitters (emitters)
  end subroutine phs_fks_set_emitters

  module subroutine phs_fks_set_momenta (phs, p)
    class(phs_fks_t), intent(inout) :: phs
    type(vector4_t), intent(in), dimension(:) :: p
    integer :: n_in, n_tot_born
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       n_in = phs%config%n_in; n_tot_born = phs%config%n_tot - 1
       phs%p_born = p(1 : n_in)
       phs%q_born = p(n_in + 1 : n_tot_born)
       phs%p_born_tot = p
    end select
  end subroutine phs_fks_set_momenta

  module subroutine phs_fks_setup_masses (phs, n_tot)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: n_tot
    call phs%generator%setup_masses (n_tot)
  end subroutine phs_fks_setup_masses

  module subroutine phs_fks_get_born_momenta (phs, p)
    class(phs_fks_t), intent(inout) :: phs
    type(vector4_t), intent(out), dimension(:) :: p
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       p(1 : phs%config%n_in) = phs%p_born
       p(phs%config%n_in + 1 :) = phs%q_born
    case (PHS_MODE_COLLINEAR_REMNANT)
       p(1:phs%config%n_in) = phs%phs_wood_t%p
       p(phs%config%n_in + 1 : ) = phs%phs_wood_t%q
    end select
    if (.not. phs%config%lab_is_cm)  p = phs%lt_cm_to_lab * p
  end subroutine phs_fks_get_born_momenta

  module subroutine phs_fks_get_outgoing_momenta (phs, q)
    class(phs_fks_t), intent(in) :: phs
    type(vector4_t), intent(out), dimension(:) :: q
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       q = phs%q_real
    case (PHS_MODE_COLLINEAR_REMNANT)
       q = phs%phs_wood_t%q
    end select
  end subroutine phs_fks_get_outgoing_momenta

  module subroutine phs_fks_get_incoming_momenta (phs, p)
    class(phs_fks_t), intent(in) :: phs
    type(vector4_t), intent(inout), dimension(:), allocatable :: p
    p = phs%p_real
  end subroutine phs_fks_get_incoming_momenta

  module subroutine phs_fks_set_isr_kinematics (phs, requires_boost)
    class(phs_fks_t), intent(inout) :: phs
    logical, intent(in) :: requires_boost
    type(vector4_t), dimension(2) :: p
    if (phs%generator%isr_kinematics%isr_mode == SQRTS_VAR) then
       if (requires_boost) then
          p = phs%lt_cm_to_lab &
               * phs%generator%real_kinematics%p_born_cms%phs_point(1)%select ([1,2])
       else
          p = phs%generator%real_kinematics%p_born_lab%phs_point(1)%select ([1,2])
       end if
       call phs%generator%set_isr_kinematics (p)
    end if
  end subroutine phs_fks_set_isr_kinematics

  module subroutine phs_fks_generate_radiation_variables &
       (phs, r_in, threshold)
    class(phs_fks_t), intent(inout) :: phs
    real(default), intent(in), dimension(:) :: r_in
    logical, intent(in) :: threshold
    type(vector4_t), dimension(:), allocatable :: p_born
    if (size (r_in) /= 3) call msg_fatal &
         ("Real kinematics need to be generated using three random numbers!")
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       allocate (p_born (size (phs%p_born_tot)))
       if (threshold) then
          p_born = phs%get_onshell_projected_momenta ()
       else
          p_born = phs%p_born_tot
       if (.not. phs%lab_is_cm ()) &
            p_born = inverse (phs%lt_cm_to_lab) * p_born
       end if
       call phs%generator%generate_radiation_variables &
            (r_in, p_born, phs%phs_identifiers, threshold)
       phs%r_real = r_in
    end select
  end subroutine phs_fks_generate_radiation_variables

  module subroutine phs_fks_compute_xi_ref_momenta (phs, p_in, contributors)
    class(phs_fks_t), intent(inout) :: phs
    type(vector4_t), intent(in), dimension(:), optional :: p_in
    type(resonance_contributors_t), intent(in), dimension(:), optional :: &
         contributors
    if (phs%mode == PHS_MODE_ADDITIONAL_PARTICLE) then
       if (present (p_in)) then
          call phs%generator%compute_xi_ref_momenta (p_in, contributors)
       else
          call phs%generator%compute_xi_ref_momenta &
               (phs%p_born_tot, contributors)
       end if
    end if
  end subroutine phs_fks_compute_xi_ref_momenta

  module subroutine phs_fks_compute_xi_ref_momenta_threshold (phs)
    class(phs_fks_t), intent(inout) :: phs
    select case (phs%mode)
    case (PHS_MODE_ADDITIONAL_PARTICLE)
       call phs%generator%compute_xi_ref_momenta_threshold &
            (phs%get_onshell_projected_momenta ())
    end select
  end subroutine phs_fks_compute_xi_ref_momenta_threshold

  module subroutine phs_fks_compute_cms_energy (phs)
    class(phs_fks_t), intent(inout) :: phs
    if (phs%mode == PHS_MODE_ADDITIONAL_PARTICLE) &
         call phs%generator%compute_cms_energy (phs%p_born_tot)
  end subroutine phs_fks_compute_cms_energy

  module subroutine phs_fks_set_reference_frames (phs, is_cms)
    class(phs_fks_t), intent(inout) :: phs
    logical, intent(in) :: is_cms
    associate (real_kinematics => phs%generator%real_kinematics)
       if (phs%config%lab_is_cm) then
          real_kinematics%p_born_cms%phs_point(1) = phs%p_born_tot
          real_kinematics%p_born_lab%phs_point(1) = phs%p_born_tot
       else
          if (is_cms) then
             real_kinematics%p_born_cms%phs_point(1) &
                  = phs%p_born_tot
             real_kinematics%p_born_lab%phs_point(1) &
                  = phs%lt_cm_to_lab * phs%p_born_tot
          else
             real_kinematics%p_born_cms%phs_point(1) &
                  = inverse (phs%lt_cm_to_lab) * phs%p_born_tot
             real_kinematics%p_born_lab%phs_point(1) &
                  = phs%p_born_tot
          end if
       end if
    end associate
  end subroutine phs_fks_set_reference_frames

  module function phs_fks_i_phs_is_isr (phs, i_phs) result (is_isr)
    logical :: is_isr
    class(phs_fks_t), intent(in) :: phs
    integer, intent(in) :: i_phs
    is_isr = phs%phs_identifiers(i_phs)%emitter <= phs%generator%n_in
  end function phs_fks_i_phs_is_isr

  module subroutine phs_fks_generator_generate_fsr_default &
       (generator, emitter, i_phs, &
       p_born, p_real, xi_y_phi, no_jacobians)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: emitter, i_phs
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    real(default), intent(in), dimension(3), optional :: xi_y_phi
    logical, intent(in), optional :: no_jacobians
    real(default) :: q0

    call generator%generate_fsr_in (p_born, p_real)
    q0 = sum (p_born(1:generator%n_in))**1
    generator%i_fsr_first = generator%n_in + 1
    call generator%generate_fsr_out (emitter, i_phs, p_born, p_real, q0, &
         xi_y_phi = xi_y_phi, no_jacobians = no_jacobians)
    if (debug_active (D_PHASESPACE)) then
       call vector4_check_momentum_conservation (p_real, generator%n_in, &
           rel_smallness = 1000 * tiny_07, abs_smallness = tiny_07)
    end if
  end subroutine phs_fks_generator_generate_fsr_default

  module subroutine phs_fks_generator_generate_fsr_resonances (generator, &
       emitter, i_phs, i_con, p_born, p_real, xi_y_phi, no_jacobians)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: emitter, i_phs
    integer, intent(in) :: i_con
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(inout), dimension(:) :: p_real
    real(default), intent(in), dimension(3), optional :: xi_y_phi
    logical, intent(in), optional :: no_jacobians
    integer, dimension(:), allocatable :: resonance_list
    integer, dimension(size(p_born)) :: inv_resonance_list
    type(vector4_t), dimension(:), allocatable :: p_tmp_born
    type(vector4_t), dimension(:), allocatable :: p_tmp_real
    type(vector4_t) :: p_resonance
    real(default) :: q0
    integer :: i, j, nlegborn, nlegreal
    integer :: i_emitter
    type(lorentz_transformation_t) :: boost_to_resonance
    integer :: n_resonant_particles
    if (debug_on) call msg_debug2 (D_PHASESPACE, "phs_fks_generator_generate_fsr_resonances")
    nlegborn = size (p_born); nlegreal = nlegborn + 1
    allocate (resonance_list (size (generator%resonance_contributors(i_con)%c)))
    resonance_list = generator%resonance_contributors(i_con)%c
    n_resonant_particles = size (resonance_list)

    if (.not. any (resonance_list == emitter)) then
       call msg_fatal ("Emitter must be included in the resonance list!")
    else
       do i = 1, n_resonant_particles
          if (resonance_list (i) == emitter) i_emitter = i
       end do
    end if

    inv_resonance_list = &
       create_inverse_resonance_list (nlegborn, resonance_list)

    allocate (p_tmp_born (n_resonant_particles))
    allocate (p_tmp_real (n_resonant_particles + 1))
    p_tmp_born = vector4_null
    p_tmp_real = vector4_null
    j = 1
    do i = 1, n_resonant_particles
       p_tmp_born(j) = p_born(resonance_list(i))
       j = j + 1
    end do

    call generator%generate_fsr_in (p_born, p_real)

    p_resonance = generator%real_kinematics%xi_ref_momenta(i_con)
    q0 = p_resonance**1

    boost_to_resonance = inverse (boost (p_resonance, q0))
    p_tmp_born = boost_to_resonance * p_tmp_born

    generator%i_fsr_first = 1
    call generator%generate_fsr_out (emitter, i_phs, p_tmp_born, p_tmp_real, &
         q0, i_emitter, xi_y_phi)
    p_tmp_real = inverse (boost_to_resonance) * p_tmp_real

    do i = generator%n_in + 1, nlegborn
       if (any (resonance_list == i)) then
          p_real(i) = p_tmp_real(inv_resonance_list (i))
       else
          p_real(i) = p_born (i)
       end if
    end do
    p_real(nlegreal) = p_tmp_real (n_resonant_particles + 1)

    if (debug_active (D_PHASESPACE)) then
       call vector4_check_momentum_conservation (p_real, generator%n_in, &
            rel_smallness = 1000 * tiny_07, abs_smallness = tiny_07)
    end if

  contains

    function create_inverse_resonance_list (nlegborn, resonance_list) &
       result (inv_resonance_list)
       integer, intent(in) :: nlegborn
       integer, intent(in), dimension(:) :: resonance_list
       integer, dimension(nlegborn) :: inv_resonance_list
       integer :: i, j
       inv_resonance_list = 0
       j = 1
       do i = 1, nlegborn
          if (any (i == resonance_list)) then
             inv_resonance_list (i) = j
             j = j + 1
          end if
       end do
    end function create_inverse_resonance_list

    function boosted_energy () result (E)
      real(default) :: E
      type(vector4_t) :: p_boost
      p_boost = boost_to_resonance * p_resonance
      E = p_boost%p(0)
    end function boosted_energy
  end subroutine phs_fks_generator_generate_fsr_resonances

  module subroutine phs_fks_generator_generate_fsr_threshold (generator, &
       emitter, i_phs, p_born, p_real, xi_y_phi)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: emitter, i_phs
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(inout), dimension(:) :: p_real
    real(default), intent(in), dimension(3), optional :: xi_y_phi
    type(vector4_t), dimension(2) :: p_tmp_born
    type(vector4_t), dimension(3) :: p_tmp_real
    integer :: nlegborn, nlegreal
    type(vector4_t) :: p_top
    real(default) :: q0
    type(lorentz_transformation_t) :: boost_to_top
    integer :: leg, other_leg
    real(default) :: sqrts, mtop
    if (debug_on)  call msg_debug2 &
         (D_PHASESPACE, "phs_fks_generator_generate_fsr_resonances")
    nlegborn = size (p_born); nlegreal = nlegborn + 1

    leg = thr_leg(emitter); other_leg = 3 - leg

    p_tmp_born(1) = p_born (ass_boson(leg))
    p_tmp_born(2) = p_born (ass_quark(leg))

    call generator%generate_fsr_in (p_born, p_real)

    p_top = generator%real_kinematics%xi_ref_momenta(leg)

    q0 = p_top**1
    sqrts = two * p_born(1)%p(0)
    mtop = m1s_to_mpole (sqrts)
    if (sqrts**2 - four * mtop**2 > zero) then
       boost_to_top = inverse (boost (p_top, q0))
    else
       boost_to_top = identity
    end if
    p_tmp_born = boost_to_top * p_tmp_born

    generator%i_fsr_first = 1
    call generator%generate_fsr_out (emitter, i_phs, p_tmp_born, &
         p_tmp_real, q0, 2, xi_y_phi)
    p_tmp_real = inverse (boost_to_top) * p_tmp_real

    p_real(ass_boson(leg)) = p_tmp_real(1)
    p_real(ass_quark(leg)) = p_tmp_real(2)
    p_real(ass_boson(other_leg)) = p_born(ass_boson(other_leg))
    p_real(ass_quark(other_leg)) = p_born(ass_quark(other_leg))
    p_real(THR_POS_GLUON) = p_tmp_real(3)

  end subroutine phs_fks_generator_generate_fsr_threshold

  module subroutine phs_fks_generator_generate_fsr_in &
       (generator, p_born, p_real)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    integer :: i
    do i = 1, generator%n_in
       p_real(i) = p_born(i)
    end do
  end subroutine phs_fks_generator_generate_fsr_in

  module subroutine phs_fks_generator_generate_fsr_out (generator, &
      emitter, i_phs, p_born, p_real, q0, p_emitter_index, &
      xi_y_phi, no_jacobians)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: emitter, i_phs
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(inout), dimension(:) :: p_real
    real(default), intent(in) :: q0
    integer, intent(in), optional :: p_emitter_index
    real(default), intent(in), dimension(3), optional :: xi_y_phi
    logical, intent(in), optional :: no_jacobians
    real(default) :: xi, y, phi
    integer :: nlegborn, nlegreal
    real(default) :: uk_np1, uk_n
    real(default) :: uk_rec, k_rec0
    type(vector3_t) :: k_n_born, k
    real(default) :: uk_n_born, uk, k2, k0_n
    real(default) :: cpsi, beta
    type(vector3_t) :: vec, vec_orth
    type(lorentz_transformation_t) :: rot
    integer :: i, p_em
    logical :: compute_jac
    p_em = emitter; if (present (p_emitter_index)) p_em = p_emitter_index
    compute_jac = .true.
    if (present (no_jacobians)) compute_jac = .not. no_jacobians
    if (generator%i_fsr_first < 0) &
       call msg_fatal ("FSR generator is called for outgoing particles but "&
            &"i_fsr_first is not set!")

    if (present (xi_y_phi)) then
       xi = xi_y_phi(I_XI)
       y = xi_y_phi(I_Y)
       phi = xi_y_phi(I_PHI)
    else
       associate (rad_var => generator%real_kinematics)
          xi = rad_var%xi_tilde
          if (rad_var%supply_xi_max) xi = xi * rad_var%xi_max(i_phs)
          y = rad_var%y(i_phs)
          phi = rad_var%phi
       end associate
    end if

    nlegborn = size (p_born)
    nlegreal = nlegborn + 1
    generator%E_gluon = q0 * xi / two
    uk_np1 = generator%E_gluon
    k_n_born = p_born(p_em)%p(1:3)
    uk_n_born = k_n_born**1

    generator%mrec2 = (q0 - p_born(p_em)%p(0))**2 &
         - space_part_norm(p_born(p_em))**2
    if (generator%is_massive(emitter)) then
       call generator%compute_emitter_kinematics (y, emitter, &
            i_phs, q0, k0_n, uk_n, uk, compute_jac)
    else
       call generator%compute_emitter_kinematics (y, q0, uk_n, uk)
       generator%real_kinematics%y_soft(i_phs) = y
       k0_n = uk_n
    end if

    if (debug_on)  call msg_debug2 &
         (D_PHASESPACE, "phs_fks_generator_generate_fsr_out")
    call debug_input_values ()

    vec = uk_n / uk_n_born * k_n_born
    vec_orth = create_orthogonal (vec)
    p_real(p_em)%p(0) = k0_n
    p_real(p_em)%p(1:3) = vec%p(1:3)
    cpsi = (uk_n**2 + uk**2 - uk_np1**2) / (two * uk_n * uk)
    !!! This is to catch the case where cpsi = 1, but numerically
    !!! turns out to be slightly larger than 1.
    call check_cpsi_bound (cpsi)
    rot = rotation (cpsi, - sqrt (one - cpsi**2), vec_orth)
    p_real(p_em) = rot * p_real(p_em)
    vec = uk_np1 / uk_n_born * k_n_born
    vec_orth = create_orthogonal (vec)
    p_real(nlegreal)%p(0) = uk_np1
    p_real(nlegreal)%p(1:3) = vec%p(1:3)
    cpsi = (uk_np1**2 + uk**2 - uk_n**2) / (two * uk_np1 * uk)
    call check_cpsi_bound (cpsi)
    rot = rotation (cpsi, sqrt (one - cpsi**2), vec_orth)
    p_real(nlegreal) = rot * p_real(nlegreal)
    call construct_recoiling_momenta ()
    if (compute_jac) call compute_jacobians ()

  contains

  subroutine debug_input_values ()
    if (debug2_active (D_PHASESPACE)) then
       call generator%write ()
       print *, 'emitter =    ', emitter
       print *, 'p_born:'
       call vector4_write_set (p_born)
       print *, 'p_real:'
       call vector4_write_set (p_real)
       print *, 'q0 =    ', q0
       if (present(p_emitter_index)) then
          print *, 'p_emitter_index =    ', p_emitter_index
       else
          print *, 'p_emitter_index not given'
       end if
    end if
  end subroutine debug_input_values

  subroutine check_cpsi_bound (cpsi)
    real(default), intent(inout) :: cpsi
    if (cpsi > one) then
       cpsi = one
    else if (cpsi < -one) then
       cpsi = - one
    end if
  end subroutine check_cpsi_bound

  subroutine construct_recoiling_momenta ()
    type(lorentz_transformation_t) :: lambda
    k_rec0 = q0 - p_real(p_em)%p(0) - p_real(nlegreal)%p(0)
    if (k_rec0**2 > generator%mrec2) then
       uk_rec = sqrt (k_rec0**2 - generator%mrec2)
    else
       uk_rec = 0
    end if
    if (generator%is_massive(emitter)) then
       beta = compute_beta (q0**2, k_rec0, uk_rec, &
              p_born(p_em)%p(0), uk_n_born)
    else
       beta = compute_beta (q0**2, k_rec0, uk_rec)
    end if
    k = p_real(p_em)%p(1:3) + p_real(nlegreal)%p(1:3)
    vec%p(1:3) = one / uk * k%p(1:3)
    lambda = boost (beta / sqrt(one - beta**2), vec)
    do i = generator%i_fsr_first, nlegborn
      if (i /= p_em) then
         p_real(i) = lambda * p_born(i)
      end if
    end do
    vec%p(1:3) = p_born(p_em)%p(1:3) / uk_n_born
    rot = rotation (cos(phi), sin(phi), vec)
    p_real(nlegreal) = rot * p_real(nlegreal)
    p_real(p_em) = rot * p_real(p_em)
  end subroutine construct_recoiling_momenta

  subroutine compute_jacobians ()
    associate (jac => generator%real_kinematics%jac(i_phs))
       if (generator%is_massive(emitter)) then
          jac%jac(1) = jac%jac(1) * four / q0 / uk_n_born / xi
       else
          k2 = two * uk_n * uk_np1* (one - y)
          jac%jac(1) = uk_n**2 / uk_n_born / (uk_n - k2 / (two * q0))
       end if
       jac%jac(2) = one
       jac%jac(3) = one - xi / two * q0 / uk_n_born
    end associate
  end subroutine compute_jacobians


  end subroutine phs_fks_generator_generate_fsr_out

  module subroutine phs_fks_generate_fsr_in (phs)
    class(phs_fks_t), intent(inout) :: phs
    type(vector4_t), dimension(:), allocatable :: p
    p = phs%generator%real_kinematics%p_born_lab%get_momenta &
         (1, phs%generator%n_in)
  end subroutine phs_fks_generate_fsr_in

  module subroutine phs_fks_generate_fsr (phs, emitter, i_phs, p_real, &
       i_con, xi_y_phi, no_jacobians)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: emitter, i_phs
    type(vector4_t), intent(out), dimension(:) :: p_real
    integer, intent(in), optional :: i_con
    real(default), intent(in), dimension(3), optional :: xi_y_phi
    logical, intent(in), optional :: no_jacobians
    type(vector4_t), dimension(:), allocatable :: p
    associate (generator => phs%generator)
      p = generator%real_kinematics%p_born_cms%phs_point(1)
      generator%real_kinematics%supply_xi_max = .true.
      if (present (i_con)) then
         call generator%generate_fsr (emitter, i_phs, i_con, p, p_real, &
              xi_y_phi, no_jacobians)
      else
         call generator%generate_fsr (emitter, i_phs, p, p_real, &
              xi_y_phi, no_jacobians)
      end if
      generator%real_kinematics%p_real_cms%phs_point(i_phs) = p_real
      if (.not. phs%config%lab_is_cm)  p_real = phs%lt_cm_to_lab * p_real
      generator%real_kinematics%p_real_lab%phs_point(i_phs) = p_real
    end associate
  end subroutine phs_fks_generate_fsr

  pure module function phs_fks_get_onshell_projected_momenta (phs) result (p)
    type(vector4_t), dimension(:), allocatable :: p
    class(phs_fks_t), intent(in) :: phs
    p = phs%generator%real_kinematics%p_born_onshell%phs_point(1)
  end function phs_fks_get_onshell_projected_momenta

  module subroutine phs_fks_generate_fsr_threshold &
       (phs, emitter, i_phs, p_real)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: emitter, i_phs
    type(vector4_t), intent(inout), dimension(:), optional :: p_real
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: pp
    integer :: leg
    associate (generator => phs%generator)
       generator%real_kinematics%supply_xi_max = .true.
       allocate (p_born (1 : generator%real_kinematics%p_born_cms%get_n_particles()))
       p_born = generator%real_kinematics%p_born_onshell%get_momenta (1)
       allocate (pp (size (p_born) + 1))
       call generator%generate_fsr_threshold (emitter, i_phs, p_born, pp)
       leg = thr_leg (emitter)
       call generator%real_kinematics%p_real_onshell(leg)%set_momenta (i_phs, pp)
       if (present (p_real))  p_real = pp
    end associate
  end subroutine phs_fks_generate_fsr_threshold

  module subroutine phs_fks_compute_xi_max_internal (phs, p, threshold)
    class(phs_fks_t), intent(inout) :: phs
    type(vector4_t), intent(in), dimension(:) :: p
    logical, intent(in) :: threshold
    integer :: i_phs, i_con, emitter
    do i_phs = 1, size (phs%phs_identifiers)
       associate (phs_id => phs%phs_identifiers(i_phs), generator => phs%generator)
          emitter = phs_id%emitter
          if (threshold) then
             call generator%compute_xi_max (emitter, i_phs, p, &
                  generator%real_kinematics%xi_max(i_phs), i_con = thr_leg(emitter))
          else if (allocated (phs_id%contributors)) then
             do i_con = 1, size (phs_id%contributors)
                call generator%compute_xi_max (emitter, i_phs, p, &
                     generator%real_kinematics%xi_max(i_phs), i_con = 1)
             end do
          else
             call generator%compute_xi_max (emitter, i_phs, p, &
                  generator%real_kinematics%xi_max(i_phs))
          end if
       end associate
    end do
  end subroutine phs_fks_compute_xi_max_internal

  module subroutine phs_fks_compute_xi_max_with_output &
       (phs, emitter, i_phs, y, p, xi_max)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: i_phs, emitter
    real(default), intent(in) :: y
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(out) :: xi_max
    call phs%generator%compute_xi_max (emitter, i_phs, p, xi_max, y_in = y)
  end subroutine phs_fks_compute_xi_max_with_output

  module subroutine phs_fks_generator_compute_emitter_kinematics_massless &
       (generator, y, q0, uk_em, uk)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: y, q0
    real(default), intent(out) :: uk_em, uk
    real(default) :: k0_np1, q2

    k0_np1 = generator%E_gluon
    q2 = q0**2

    uk_em = (q2 - generator%mrec2 - two * q0 * k0_np1) / &
         (two * (q0 - k0_np1 * (one - y)))
    uk = sqrt (uk_em**2 + k0_np1**2 + two * uk_em * k0_np1 * y)
  end subroutine phs_fks_generator_compute_emitter_kinematics_massless

  module subroutine phs_fks_generator_compute_emitter_kinematics_massive &
       (generator, y, em, i_phs, q0, k0_em, uk_em, uk, compute_jac)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: y
    integer, intent(in) :: em, i_phs
    real(default), intent(in) :: q0
    real(default), intent(inout) :: k0_em, uk_em, uk
    logical, intent(in) :: compute_jac
    real(default) :: k0_np1, q2, mrec2, m2
    real(default) :: k0_rec_max, k0_em_max, k0_rec, uk_rec
    real(default) :: z, z1, z2

    k0_np1 = generator%E_gluon
    q2 = q0**2
    mrec2 = generator%mrec2
    m2 = generator%m2(em)

    k0_rec_max = (q2 - m2 + mrec2) / (two * q0)
    k0_em_max = (q2 + m2 - mrec2)  /(two * q0)
    z1 = (k0_rec_max + sqrt (k0_rec_max**2 - mrec2)) / q0
    z2 = (k0_rec_max - sqrt (k0_rec_max**2 - mrec2)) / q0
    z = z2 - (z2 - z1) * (one + y) / two
    k0_em = k0_em_max - k0_np1 * z
    k0_rec = q0 - k0_np1 - k0_em
    uk_em = sqrt(k0_em**2 - m2)
    uk_rec = sqrt(k0_rec**2 - mrec2)
    uk = uk_rec
    if (compute_jac) &
         generator%real_kinematics%jac(i_phs)%jac = q0 * (z1 - z2) / four * k0_np1
    generator%real_kinematics%y_soft(i_phs) = &
       (two * q2 * z - q2 - mrec2 + m2) / (sqrt(k0_em_max**2 - m2) * q0) / two
  end subroutine phs_fks_generator_compute_emitter_kinematics_massive

  function recompute_xi_max (q0, mrec2, m2, y) result (xi_max)
    real(default) :: xi_max
    real(default), intent(in) :: q0, mrec2, m2, y
    real(default) :: q2, k0_np1_max, k0_rec_max
    real(default) :: z1, z2, z
    q2 = q0**2
    k0_rec_max = (q2 - m2 + mrec2) / (two * q0)
    z1 = (k0_rec_max + sqrt (k0_rec_max**2 - mrec2)) / q0
    z2 = (k0_rec_max - sqrt (k0_rec_max**2 - mrec2)) / q0
    z = z2 - (z2 - z1) * (one + y) / 2
    k0_np1_max = - (q2 * z**2 - two * q0 * k0_rec_max * z + mrec2) / (two * q0 * z * (one - z))
    xi_max = two * k0_np1_max / q0
  end function recompute_xi_max

  function compute_beta_massless (q2, k0_rec, uk_rec) result (beta)
    real(default), intent(in) :: q2, k0_rec, uk_rec
    real(default) :: beta
    beta = (q2 - (k0_rec + uk_rec)**2) / (q2 + (k0_rec + uk_rec)**2)
  end function compute_beta_massless

  function compute_beta_massive (q2, k0_rec, uk_rec, &
     k0_em_born, uk_em_born) result (beta)
    real(default), intent(in) :: q2, k0_rec, uk_rec
    real(default), intent(in) :: k0_em_born, uk_em_born
    real(default) :: beta
    real(default) :: k0_rec_born, uk_rec_born, alpha
    k0_rec_born = sqrt(q2) - k0_em_born
    uk_rec_born = uk_em_born
    alpha = (k0_rec + uk_rec) / (k0_rec_born + uk_rec_born)
    beta = (one - alpha**2) / (one + alpha**2)
  end function compute_beta_massive

  pure function get_xi_max_fsr_massless (p_born, q0, emitter) result (xi_max)
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: q0
    integer, intent(in) :: emitter
    real(default) :: xi_max
    real(default) :: uk_n_born
    uk_n_born = space_part_norm (p_born(emitter))
    xi_max = two * uk_n_born / q0
  end function get_xi_max_fsr_massless

  pure function get_xi_max_fsr_massive (p_born, q0, emitter, m2, y) result (xi_max)
    real(default) :: xi_max
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: q0
    integer, intent(in) :: emitter
    real(default), intent(in) :: m2, y
    real(default) :: mrec2
    real(default) :: k0_rec_max
    real(default) :: z, z1, z2
    real(default) :: k0_np1_max
    associate (p => p_born(emitter)%p)
       mrec2 = (q0 - p(0))**2 - p(1)**2 - p(2)**2 - p(3)**2
    end associate
    call compute_dalitz_bounds (q0, m2, mrec2, z1, z2, k0_rec_max)
    z = z2 - (z2 - z1) * (one + y) / two
    k0_np1_max = - (q0**2 * z**2 - two * q0 * k0_rec_max * z + mrec2) &
       / (two * q0 * z * (one - z))
    xi_max = two * k0_np1_max / q0
  end function get_xi_max_fsr_massive

  module function get_xi_max_isr (xb, y) result (xi_max)
    real(default) :: xi_max
    real(default), dimension(2), intent(in) :: xb
    real(default), intent(in) :: y
    xi_max = one - max (xi_max_isr_plus (xb(I_PLUS), y), xi_max_isr_minus (xb(I_MINUS), y))
  contains
    function xi_max_isr_plus (x, y)
      real(default) :: xi_max_isr_plus
      real(default), intent(in) :: x, y
      real(default) :: deno
      deno = sqrt ((one + x**2)**2 * (one - y)**2 + 16 * y * x**2) + (one - y) * (1 - x**2)
      xi_max_isr_plus = two * (one + y) * x**2 / deno
    end function xi_max_isr_plus

    function xi_max_isr_minus (x, y)
      real(default) :: xi_max_isr_minus
      real(default), intent(in) :: x, y
      real(default) :: deno
      deno = sqrt ((one + x**2)**2 * (one + y)**2 - 16 * y * x**2) + (one + y) * (1 - x**2)
      xi_max_isr_minus = two * (one - y) * x**2 / deno
    end function xi_max_isr_minus
  end function get_xi_max_isr

  recursive function get_xi_max_isr_decay (p) result (xi_max)
     real(default) :: xi_max
     type(vector4_t), dimension(:), intent(in) :: p
     integer :: n_tot
     type(vector4_t), dimension(:), allocatable :: p_dec_new
     n_tot = size (p)
     if (n_tot == 3) then
        xi_max = xi_max_one_to_two (p(1), p(2), p(3))
     else
        allocate (p_dec_new (n_tot - 1))
        p_dec_new(1) = sum (p (3 : ))
        p_dec_new(2 : n_tot - 1) = p (3 : n_tot)
        xi_max = min (xi_max_one_to_two (p(1), p(2), sum(p(3 : ))), &
             get_xi_max_isr_decay (p_dec_new))
     end if
  contains
    function xi_max_one_to_two (p_in, p_out1, p_out2) result (xi_max)
      real(default) :: xi_max
      type(vector4_t), intent(in) :: p_in, p_out1, p_out2
      real(default) :: m_in, m_out1, m_out2
      m_in = p_in**1
      m_out1 = p_out1**1; m_out2 = p_out2**1
      xi_max = one - (m_out1 + m_out2)**2 / m_in**2
    end function xi_max_one_to_two
  end function get_xi_max_isr_decay

  module subroutine phs_fks_generate_isr (phs, i_phs, p_real)
    class(phs_fks_t), intent(inout) :: phs
    integer, intent(in) :: i_phs
    type(vector4_t), intent(out), dimension(:) :: p_real
    type(vector4_t) :: p0, p1
    type(lorentz_transformation_t) :: lt
    real(default) :: sqrts_hat
    type(vector4_t), dimension(:), allocatable :: p_work

    associate (generator => phs%generator)
       select case (generator%n_in)
       case (1)
          p_work = generator%real_kinematics%p_born_cms%phs_point(1)
          call generator%generate_isr_fixed_beam_energy (i_phs, p_work, p_real)
          phs%config%lab_is_cm = .true.
       case (2)
          select case (generator%isr_kinematics%isr_mode)
          case (SQRTS_FIXED)
             p_work = generator%real_kinematics%p_born_cms%phs_point(1)
             call generator%generate_isr_fixed_beam_energy (i_phs, p_work, p_real)
          case (SQRTS_VAR)
             p_work = generator%real_kinematics%p_born_lab%phs_point(1)
             call generator%generate_isr (i_phs, p_work, p_real)
          end select
       end select
       generator%real_kinematics%p_real_lab%phs_point(i_phs) = p_real
       if (.not. phs%config%lab_is_cm) then
          sqrts_hat = (p_real(1) + p_real(2))**1
          p0 = p_real(1) + p_real(2)
          lt = boost (p0, sqrts_hat)
          p1 = inverse(lt) * p_real(1)
          lt = lt * rotation_to_2nd (3, space_part (p1))
          phs%generator%real_kinematics%p_real_cms%phs_point(i_phs) &
               = inverse (lt) * p_real
       else
          phs%generator%real_kinematics%p_real_cms%phs_point(i_phs) &
               = p_real
       end if
     end associate
  end subroutine phs_fks_generate_isr

  module subroutine phs_fks_generator_generate_isr_fixed_beam_energy &
       (generator, i_phs, p_born, p_real)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: i_phs
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    real(default) :: xi_max, xi, y, phi
    integer :: nlegborn, nlegreal, i
    real(default) :: k0_np1
    real(default) :: msq_in
    type(vector4_t) :: p_virt
    real(default) :: jac_real

    associate (rad_var => generator%real_kinematics)
      xi_max = rad_var%xi_max(i_phs)
      xi = rad_var%xi_tilde * xi_max
      y = rad_var%y(i_phs)
      phi = rad_var%phi
      rad_var%y_soft(i_phs) = y
    end associate

    nlegborn = size (p_born)
    nlegreal = nlegborn + 1

    msq_in = sum (p_born(1:generator%n_in))**2
    generator%real_kinematics%jac(i_phs)%jac = one

    p_real(1) = p_born(1)
    if (generator%n_in > 1) p_real(2) = p_born(2)
    k0_np1 = zero
    do i = 1, generator%n_in
       k0_np1 = k0_np1 + p_real(i)%p(0) * xi / two
    end do
    p_real(nlegreal)%p(0) = k0_np1
    p_real(nlegreal)%p(1) = k0_np1 * sqrt(one - y**2) * sin(phi)
    p_real(nlegreal)%p(2) = k0_np1 * sqrt(one - y**2) * cos(phi)
    p_real(nlegreal)%p(3) = k0_np1 * y

    p_virt = sum (p_real(1:generator%n_in)) - p_real(nlegreal)

    jac_real = one
    call generate_on_shell_decay (p_virt, &
         p_born(generator%n_in + 1 : nlegborn), &
         p_real(generator%n_in + 1 : nlegreal - 1), 1, msq_in, jac_real)

    associate (jac => generator%real_kinematics%jac(i_phs))
      jac%jac(1) = jac_real
      jac%jac(2) = one
    end associate

  end subroutine phs_fks_generator_generate_isr_fixed_beam_energy

  module subroutine phs_fks_generator_generate_isr_factorized &
       (generator, i_phs, emitter, p_born, p_real)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: i_phs, emitter
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    type(vector4_t), dimension(3) :: p_tmp_born
    type(vector4_t), dimension(4) :: p_tmp_real
    type(vector4_t) :: p_top
    type(lorentz_transformation_t) :: boost_to_rest_frame
    integer, parameter :: nlegreal = 7
    !!! Factorized phase space so far only required for ee -> bwbw

    p_tmp_born = vector4_null; p_tmp_real = vector4_null
    p_real(1:2) = p_born(1:2)
    if (emitter == THR_POS_B) then
       p_top = p_born (THR_POS_WP) + p_born (THR_POS_B)
       p_tmp_born(2) = p_born (THR_POS_WP)
       p_tmp_born(3) = p_born (THR_POS_B)
    else if (emitter == THR_POS_BBAR) then
       p_top = p_born (THR_POS_WM) + p_born (THR_POS_BBAR)
       p_tmp_born(2) = p_born (THR_POS_WM)
       p_tmp_born(3) = p_born (THR_POS_BBAR)
    else
       call msg_fatal ("Threshold computation requires emitters to be at position 5 and 6 " // &
          "Please check if your process specification fulfills this requirement.")
    end if
    p_tmp_born (1) = p_top
    boost_to_rest_frame = inverse (boost (p_top, p_top**1))
    p_tmp_born = boost_to_rest_frame * p_tmp_born
    call generator%compute_xi_max_isr_factorized (i_phs, p_tmp_born)
    call generator%generate_isr_fixed_beam_energy &
         (i_phs, p_tmp_born, p_tmp_real)
    p_tmp_real = inverse (boost_to_rest_frame) * p_tmp_real
    if (emitter == THR_POS_B) then
       p_real(THR_POS_WP) = p_tmp_real(2)
       p_real(THR_POS_B) = p_tmp_real(3)
       p_real(THR_POS_WM) = p_born(THR_POS_WM)
       p_real(THR_POS_BBAR) = p_born(THR_POS_BBAR)
    !!! Exception has been handled above
    else
       p_real(THR_POS_WM) = p_tmp_real(2)
       p_real(THR_POS_BBAR) = p_tmp_real(3)
       p_real(THR_POS_WP) = p_born(THR_POS_WP)
       p_real(THR_POS_B) = p_born(THR_POS_B)
    end if
    p_real(nlegreal) = p_tmp_real(4)
  end subroutine phs_fks_generator_generate_isr_factorized

  module subroutine phs_fks_generator_generate_isr &
       (generator, i_phs, p_born, p_real)
    !!! Important: Import Born momenta in the lab frame
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: i_phs
    type(vector4_t), intent(in) , dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    real(default) :: xi_max, xi_tilde, xi, y, phi
    integer :: nlegborn, nlegreal
    real(default) :: sqrts_real
    real(default) :: k0_np1
    type(lorentz_transformation_t) :: &
         lambda_transv, lambda_longit, lambda_longit_inv
    real(default) :: x_plus, x_minus, xb_plus, xb_minus
    real(default) :: onemy, onepy
    integer :: i
    real(default) :: xi_plus, xi_minus
    real(default) :: beta_gamma
    type(vector3_t) :: beta_vec

    associate (rad_var => generator%real_kinematics)
       xi_max = rad_var%xi_max(i_phs)
       xi_tilde = rad_var%xi_tilde
       xi = xi_tilde * xi_max
       y = rad_var%y(i_phs)
       onemy = one - y; onepy = one + y
       phi = rad_var%phi
       rad_var%y_soft(i_phs) = y
    end associate

    nlegborn = size (p_born)
    nlegreal = nlegborn + 1
    generator%isr_kinematics%sqrts_born = (p_born(1) + p_born(2))**1

    !!! Initial state real momenta
    xb_plus = generator%isr_kinematics%x(I_PLUS)
    xb_minus = generator%isr_kinematics%x(I_MINUS)
    x_plus = xb_plus / sqrt(one - xi) * sqrt ((two - xi * onemy) / &
         (two - xi * onepy))
    x_minus = xb_minus / sqrt(one - xi) * sqrt ((two - xi * onepy) / &
         (two - xi * onemy))
    xi_plus = xi_tilde * (one - xb_plus)
    xi_minus = xi_tilde * (one - xb_minus)
    p_real(I_PLUS) = x_plus / xb_plus * p_born(I_PLUS)
    p_real(I_MINUS) = x_minus / xb_minus * p_born(I_MINUS)

    !!! Fraction of momentum fractions in a collinear splitting
    generator%isr_kinematics%z(I_PLUS) = (one - xi_plus)
    generator%isr_kinematics%z(I_MINUS) = (one - xi_minus)

    !!! Create radiation momentum in the real CMS
    sqrts_real = generator%isr_kinematics%sqrts_born / sqrt (one - xi)
    k0_np1 = sqrts_real * xi / two
    p_real(nlegreal)%p(0) = k0_np1
    p_real(nlegreal)%p(1) = k0_np1 * sqrt (one - y**2) * sin(phi)
    p_real(nlegreal)%p(2) = k0_np1 * sqrt (one - y**2) * cos(phi)
    p_real(nlegreal)%p(3) = k0_np1 * y

    !!! Boosts the radiation from real CMS to the real LAB frame
    call get_boost_parameters (p_real, beta_gamma, beta_vec)
    lambda_longit = create_longitudinal_boost &
         (beta_gamma, beta_vec, inverse = .true.)
    p_real(nlegreal) = lambda_longit * p_real(nlegreal)

    call get_boost_parameters (p_born, beta_gamma, beta_vec)
    lambda_longit = create_longitudinal_boost &
         (beta_gamma, beta_vec, inverse = .false.)
    forall (i = 3 : nlegborn)  p_real(i) = lambda_longit * p_born(i)

    lambda_transv = create_transversal_boost (p_real(nlegreal), xi, sqrts_real)
    forall (i = 3 : nlegborn)  p_real(i) = lambda_transv * p_real(i)

    lambda_longit_inv = create_longitudinal_boost &
         (beta_gamma, beta_vec, inverse = .true.)
    forall (i = 3 : nlegborn)  p_real(i) = lambda_longit_inv * p_real(i)

    !!! Compute Jacobians
    associate (jac => generator%real_kinematics%jac(i_phs))
       !!! Additional 1 / (1 - xi) factor because in the real jacobian,
       !!! there is s_real in the numerator
       !!! We also have to adapt the flux factor, which is 1/(2s_real) for
       !!! the real component
       !!! The reweighting factor is s_born / s_real, cancelling the
       !!! (1-xi) factor from above
       jac%jac(1) = one / (one - xi)
       jac%jac(2) = one
       jac%jac(3) = one / (one - xi_plus)**2
       jac%jac(4) = one / (one - xi_minus)**2
    end associate
  contains
    subroutine get_boost_parameters (p, beta_gamma, beta_vec)
       type(vector4_t), intent(in), dimension(:) :: p
       real(default), intent(out) :: beta_gamma
       type(vector3_t), intent(out) :: beta_vec
       beta_vec = (p(1)%p(1:3) + p(2)%p(1:3)) / (p(1)%p(0) + p(2)%p(0))
       beta_gamma = beta_vec**1 / sqrt (one - beta_vec**2)
       beta_vec = beta_vec / beta_vec**1
    end subroutine get_boost_parameters

    function create_longitudinal_boost &
         (beta_gamma, beta_vec, inverse) result (lambda)
       real(default), intent(in) :: beta_gamma
       type(vector3_t), intent(in) :: beta_vec
       logical, intent(in) :: inverse
       type(lorentz_transformation_t) :: lambda
       if (inverse) then
          lambda = boost (beta_gamma, beta_vec)
       else
          lambda = boost (-beta_gamma, beta_vec)
       end if
    end function create_longitudinal_boost

    function create_transversal_boost (p_rad, xi, sqrts_real) result (lambda)
       type(vector4_t), intent(in) :: p_rad
       real(default), intent(in) :: xi, sqrts_real
       type(lorentz_transformation_t) :: lambda
       type(vector3_t) :: vec_transverse
       real(default) :: pt2, beta, beta_gamma
       pt2 = transverse_part (p_rad)**2
       beta = one / sqrt (one + sqrts_real**2 * (one - xi) / pt2)
       beta_gamma = beta / sqrt (one - beta**2)
       vec_transverse%p(1:2) = p_rad%p(1:2)
       vec_transverse%p(3) = zero
       vec_transverse = normalize (vec_transverse)
       lambda = boost (-beta_gamma, vec_transverse)
    end function create_transversal_boost
  end subroutine phs_fks_generator_generate_isr

  module subroutine phs_fks_generator_set_sqrts_hat (generator, sqrts)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: sqrts
    generator%sqrts = sqrts
  end subroutine phs_fks_generator_set_sqrts_hat

  module subroutine phs_fks_generator_set_emitters (generator, emitters)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in), dimension(:), allocatable ::  emitters
    allocate (generator%emitters (size (emitters)))
    generator%emitters = emitters
  end subroutine phs_fks_generator_set_emitters

  module subroutine phs_fks_generator_setup_masses (generator, n_tot)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: n_tot
    if (.not. allocated (generator%m2)) then
       allocate (generator%is_massive (n_tot))
       allocate (generator%m2 (n_tot))
       generator%is_massive = .false.
       generator%m2 = zero
    end if
  end subroutine phs_fks_generator_setup_masses

  module subroutine phs_fks_generator_set_xi_and_y_bounds &
       (generator, fks_xi_min, fks_y_max)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in), optional :: fks_xi_min, fks_y_max
    real(default) :: xi_min, y_max
    xi_min = zero; y_max = one
    if (present (fks_xi_min)) xi_min = fks_xi_min
    if (present (fks_y_max)) y_max = fks_y_max
    generator%xi_min = min (one, max (xi_min, tiny_07))
    generator%y_max = min (abs (y_max), one)
  end subroutine phs_fks_generator_set_xi_and_y_bounds

  module subroutine phs_fks_generator_set_isr_kinematics (generator, p)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(vector4_t), dimension(2), intent(in) :: p
    if (allocated (generator%isr_kinematics%beam_energy)) then
       select case (size (generator%isr_kinematics%beam_energy))
       case (1)
          generator%isr_kinematics%x = p%p(0) / &
               generator%isr_kinematics%beam_energy(1)
       case (2)
          generator%isr_kinematics%x = p%p(0) / &
               generator%isr_kinematics%beam_energy
       end select
    else
       generator%isr_kinematics%x = 0
    end if
  end subroutine phs_fks_generator_set_isr_kinematics

  module subroutine phs_fks_generator_generate_radiation_variables &
       (generator, r_in, p_born, phs_identifiers, threshold)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in), dimension(:) :: r_in
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
    logical, intent(in), optional :: threshold

    associate (rad_var => generator%real_kinematics)
       rad_var%phi = r_in (I_PHI) * twopi
       select case (generator%mode)
       case (GEN_REAL_PHASE_SPACE)
          rad_var%jac_rand = twopi
          call generator%compute_y_real_phs (r_in(I_Y), p_born, phs_identifiers, &
               rad_var%jac_rand, rad_var%y, threshold)
       case (GEN_SOFT_MISMATCH)
          rad_var%jac_mismatch = twopi
          call generator%compute_y_mismatch (r_in(I_Y), rad_var%jac_mismatch, &
               rad_var%y_mismatch, rad_var%y_soft)
       case default
          call generator%compute_y_test (rad_var%y)
       end select
       call generator%compute_xi_tilde (r_in(I_XI))
       call generator%set_masses (p_born, phs_identifiers)
    end associate
  end subroutine phs_fks_generator_generate_radiation_variables

  module subroutine phs_fks_generator_compute_xi_ref_momenta &
       (generator, p_born, resonance_contributors)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(resonance_contributors_t), intent(in), dimension(:), optional &
       :: resonance_contributors
    integer :: i_con, n_contributors
    if (present (resonance_contributors)) then
       n_contributors = size (resonance_contributors)
       if (.not. allocated (generator%resonance_contributors)) &
          allocate (generator%resonance_contributors (n_contributors))
       do i_con = 1, n_contributors
          generator%real_kinematics%xi_ref_momenta(i_con) = &
               get_resonance_momentum (p_born, resonance_contributors(i_con)%c)
          generator%resonance_contributors(i_con) = &
               resonance_contributors(i_con)
       end do
    else
       generator%real_kinematics%xi_ref_momenta(1) = &
            sum (p_born(1:generator%n_in))
    end if
  end subroutine phs_fks_generator_compute_xi_ref_momenta

  module subroutine phs_fks_generator_compute_xi_ref_momenta_threshold &
       (generator, p_born)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(vector4_t), intent(in), dimension(:) :: p_born
    generator%real_kinematics%xi_ref_momenta(1) = &
         p_born(THR_POS_WP) + p_born(THR_POS_B)
    generator%real_kinematics%xi_ref_momenta(2) = &
         p_born(THR_POS_WM) + p_born(THR_POS_BBAR)
  end subroutine phs_fks_generator_compute_xi_ref_momenta_threshold

  module subroutine phs_fks_generator_compute_cms_energy (generator, p_born)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t) :: p_sum
    p_sum = sum (p_born (1 : generator%n_in))
    generator%real_kinematics%cms_energy2 = p_sum**2
  end subroutine phs_fks_generator_compute_cms_energy

  module subroutine phs_fks_generator_compute_xi_max (generator, emitter, &
       i_phs, p, xi_max, i_con, y_in)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: i_phs, emitter
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(out) :: xi_max
    integer, intent(in), optional :: i_con
    real(default), intent(in), optional :: y_in
    real(default) :: q0
    type(vector4_t), dimension(:), allocatable :: pp, pp_decay
    type(vector4_t) :: p_res
    type(lorentz_transformation_t) :: L_to_resonance
    real(default) :: y
    if (.not. any (generator%emitters == emitter)) return
    allocate (pp (size (p)))
    associate (rad_var => generator%real_kinematics)
       if (present (i_con)) then
          q0 = rad_var%xi_ref_momenta(i_con)**1
       else
          q0 = energy (sum (p(1:generator%n_in)))
       end if
       if (present (y_in)) then
          y = y_in
       else
          y = rad_var%y(i_phs)
       end if
       if (present (i_con)) then
          p_res = rad_var%xi_ref_momenta(i_con)
          L_to_resonance = inverse (boost (p_res, q0))
          pp = L_to_resonance * p
       else
          pp = p
       end if
       if (emitter <= generator%n_in) then
          select case (generator%isr_kinematics%isr_mode)
          case (SQRTS_FIXED)
             if (generator%n_in > 1) then
                allocate (pp_decay (size (pp) - 1))
             else
                allocate (pp_decay (size (pp)))
             end if
             pp_decay (1) = sum (pp(1:generator%n_in))
             pp_decay (2 : ) = pp (generator%n_in + 1 : )
             xi_max = get_xi_max_isr_decay (pp_decay)
             deallocate (pp_decay)
          case (SQRTS_VAR)
             xi_max = get_xi_max_isr (generator%isr_kinematics%x, y)
          end select
       else
          if (generator%is_massive(emitter)) then
             xi_max = get_xi_max_fsr (pp, q0, emitter, generator%m2(emitter), y)
          else
             xi_max = get_xi_max_fsr (pp, q0, emitter)
          end if
       end if
       deallocate (pp)
    end associate
  end subroutine phs_fks_generator_compute_xi_max

  module subroutine phs_fks_generator_compute_xi_max_isr_factorized &
     (generator, i_phs, p)
    class(phs_fks_generator_t), intent(inout) :: generator
    integer, intent(in) :: i_phs
    type(vector4_t), intent(in), dimension(:) :: p
    generator%real_kinematics%xi_max(i_phs) = get_xi_max_isr_decay (p)
  end subroutine phs_fks_generator_compute_xi_max_isr_factorized

  module subroutine phs_fks_generator_set_masses &
       (generator, p, phs_identifiers)
    class(phs_fks_generator_t), intent(inout) :: generator
    type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
    type(vector4_t), intent(in), dimension(:) :: p
    integer :: emitter, i_phs
    do i_phs = 1, size (phs_identifiers)
       emitter = phs_identifiers(i_phs)%emitter
       if (any (generator%emitters == emitter) .and. emitter > 0) then
          if (generator%is_massive (emitter) .and. emitter > generator%n_in) &
             generator%m2(emitter) = p(emitter)**2
       end if
    end do
  end subroutine phs_fks_generator_set_masses

  module subroutine compute_y_from_emitter (r_y, p, n_in, emitter, &
       massive, y_max, jac_rand, y, contributors, threshold)
    real(default), intent(in) :: r_y
    type(vector4_t), intent(in), dimension(:) :: p
    integer, intent(in) :: n_in
    integer, intent(in) :: emitter
    logical, intent(in) :: massive
    real(default), intent(in) :: y_max
    real(default), intent(inout) :: jac_rand
    real(default), intent(out) :: y
    integer, intent(in), dimension(:), allocatable, optional :: contributors
    logical, intent(in), optional :: threshold
    logical :: thr, resonance
    type(vector4_t) :: p_res, p_em
    real(default) :: q0
    type(lorentz_transformation_t) :: boost_to_resonance
    integer :: i
    real(default) :: beta, one_m_beta, one_p_beta
    thr = .false.; if (present (threshold)) thr = threshold
    p_res = vector4_null
    if (present (contributors)) then
       resonance = allocated (contributors)
    else
       resonance = .false.
    end if
    if (massive) then
       if (resonance) then
          do i = 1, size (contributors)
             p_res = p_res + p(contributors(i))
          end do
       else if (thr) then
          p_res = p(ass_boson(thr_leg(emitter))) + p(ass_quark(thr_leg(emitter)))
       else
          p_res = sum (p(1:n_in))
       end if
       q0 = p_res**1
       boost_to_resonance = inverse (boost (p_res, q0))
       p_em = boost_to_resonance * p(emitter)
       beta = beta_emitter (q0, p_em)
       one_m_beta = one - beta
       one_p_beta = one + beta
       y = one / beta * (one - one_p_beta * &
              exp ( - r_y * log(one_p_beta / one_m_beta)))
       jac_rand = jac_rand * &
              (one - beta * y) * log(one_p_beta / one_m_beta) / beta
    else
       y = (one - two * r_y) * y_max
       jac_rand = jac_rand * 3 * (one - y**2) * y_max
       y = 1.5_default * (y - y**3 / 3)
    end if
  end subroutine compute_y_from_emitter

  module subroutine phs_fks_generator_compute_y_real_phs &
       (generator, r_y, p, phs_identifiers, &
       jac_rand, y, threshold)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: r_y
    type(vector4_t), intent(in), dimension(:) :: p
    type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
    real(default), intent(inout), dimension(:) :: jac_rand
    real(default), intent(out), dimension(:) :: y
    logical, intent(in), optional :: threshold
    real(default) :: beta, one_p_beta, one_m_beta
    type(lorentz_transformation_t) :: boost_to_resonance
    real(default) :: q0
    type(vector4_t) :: p_res, p_em
    integer :: i, i_phs, emitter
    logical :: thr
    logical :: construct_massive_fsr
    construct_massive_fsr = .false.
    thr = .false.; if (present (threshold)) thr = threshold
    do i_phs = 1, size (phs_identifiers)
       emitter = phs_identifiers(i_phs)%emitter
       !!! We need this additional check because of decay phase spaces
       !!! t -> bW has a massive emitter at position 1, which should
       !!! not be treated here.
       construct_massive_fsr = emitter > generator%n_in
       if (construct_massive_fsr) construct_massive_fsr = &
            construct_massive_fsr .and. generator%is_massive (emitter)
       call compute_y_from_emitter (r_y, p, generator%n_in, &
            emitter, construct_massive_fsr, &
            generator%y_max, jac_rand(i_phs), y(i_phs), &
            phs_identifiers(i_phs)%contributors, threshold)
    end do
  end subroutine phs_fks_generator_compute_y_real_phs

  module subroutine phs_fks_generator_compute_y_mismatch &
       (generator, r_y, jac_rand, y, y_soft)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: r_y
    real(default), intent(inout) :: jac_rand
    real(default), intent(out) :: y
    real(default), intent(out), dimension(:) :: y_soft
    y = (one - two * r_y) * generator%y_max
    jac_rand = jac_rand * 3 * (one - y**2) * generator%y_max
    y = 1.5_default * (y - y**3 / 3)
    y_soft = y
  end subroutine phs_fks_generator_compute_y_mismatch

  module subroutine phs_fks_generator_compute_y_test (generator, y)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(out), dimension(:):: y
    select case (generator%mode)
    case (GEN_SOFT_LIMIT_TEST)
       y = y_test_soft
    case (GEN_COLL_LIMIT_TEST)
       y = y_test_coll
    case (GEN_ANTI_COLL_LIMIT_TEST)
       y = - y_test_coll
    case (GEN_SOFT_COLL_LIMIT_TEST)
       y = y_test_coll
    case (GEN_SOFT_ANTI_COLL_LIMIT_TEST)
       y = - y_test_coll
    end select
  end subroutine phs_fks_generator_compute_y_test

  pure module function beta_emitter (q0, p) result (beta)
    real(default), intent(in) :: q0
    type(vector4_t), intent(in) :: p
    real(default) :: beta
    real(default) :: m2, mrec2, k0_max
    m2 = p**2
    mrec2 = (q0 - p%p(0))**2 - p%p(1)**2 - p%p(2)**2 - p%p(3)**2
    k0_max = (q0**2 - mrec2 + m2) / (two * q0)
    beta = sqrt(one - m2 / k0_max**2)
  end function beta_emitter

  pure module subroutine phs_fks_generator_compute_xi_tilde (generator, r)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: r
    real(default) :: deno
    associate (rad_var => generator%real_kinematics)
       select case (generator%mode)
       case (GEN_REAL_PHASE_SPACE)
          if (generator%singular_jacobian) then
             rad_var%xi_tilde = (one - generator%xi_min) - (one - r)**2 * &
                (one - two * generator%xi_min)
             rad_var%jac_rand = rad_var%jac_rand * two * (one - r) * &
                (one - two * generator%xi_min)
          else
             rad_var%xi_tilde = generator%xi_min + r * (one - generator%xi_min)
             rad_var%jac_rand = rad_var%jac_rand * (one - generator%xi_min)
          end if
       case (GEN_SOFT_MISMATCH)
          deno = one - r
          if (deno < tiny_13) deno = tiny_13
          rad_var%xi_mismatch = generator%xi_min + r / deno
          rad_var%jac_mismatch = rad_var%jac_mismatch / deno**2
       case (GEN_SOFT_LIMIT_TEST)
          rad_var%xi_tilde = r * two * xi_tilde_test_soft
          rad_var%jac_rand = two * xi_tilde_test_soft
       case (GEN_COLL_LIMIT_TEST)
          rad_var%xi_tilde = xi_tilde_test_coll
          rad_var%jac_rand = xi_tilde_test_coll
       case (GEN_ANTI_COLL_LIMIT_TEST)
          rad_var%xi_tilde = xi_tilde_test_coll
          rad_var%jac_rand = xi_tilde_test_coll
       case (GEN_SOFT_COLL_LIMIT_TEST)
          rad_var%xi_tilde = r * two * xi_tilde_test_soft
          rad_var%jac_rand = two * xi_tilde_test_soft
       case (GEN_SOFT_ANTI_COLL_LIMIT_TEST)
          rad_var%xi_tilde = r * two * xi_tilde_test_soft
          rad_var%jac_rand = two * xi_tilde_test_soft
       end select
    end associate
  end subroutine phs_fks_generator_compute_xi_tilde

  module subroutine phs_fks_generator_prepare_generation (generator, &
       r_in, i_phs, emitter, p_born, phs_identifiers, contributors, i_con)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), dimension(3), intent(in) :: r_in
    integer, intent(in) :: i_phs, emitter
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(phs_identifier_t), intent(in), dimension(:) :: phs_identifiers
    type(resonance_contributors_t), intent(in), dimension(:), optional :: &
         contributors
    integer, intent(in), optional :: i_con
    call generator%generate_radiation_variables (r_in, p_born, phs_identifiers)
    call generator%compute_xi_ref_momenta &
         (generator%real_kinematics%p_born_lab%phs_point(1)%get (), &
         contributors)
    call generator%compute_xi_max (emitter, i_phs, p_born, &
         generator%real_kinematics%xi_max(i_phs), i_con = i_con)
  end subroutine phs_fks_generator_prepare_generation

  module subroutine phs_fks_generator_generate_fsr_from_xi_and_y &
       (generator, xi, y, &
     phi, emitter, i_phs, p_born, p_real)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: xi, y, phi
    integer, intent(in) :: emitter, i_phs
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    associate (rad_var => generator%real_kinematics)
       rad_var%supply_xi_max = .false.
       rad_var%xi_tilde = xi
       rad_var%y(i_phs) = y
       rad_var%phi = phi
    end associate
    call generator%set_sqrts_hat (p_born(1)%p(0) + p_born(2)%p(0))
    call generator%generate_fsr (emitter, i_phs, p_born, p_real)
  end subroutine phs_fks_generator_generate_fsr_from_xi_and_y

  module subroutine phs_fks_generator_generate_isr_from_xi_and_y &
       (generator, xi, xi_max, y, phi, i_phs, p_born, p_real)
    class(phs_fks_generator_t), intent(inout) :: generator
    real(default), intent(in) :: xi, xi_max, y, phi
    integer, intent(in) :: i_phs
    type(vector4_t), intent(in), dimension(:) :: p_born
    type(vector4_t), intent(out), dimension(:) :: p_real
    associate (rad_var => generator%real_kinematics)
       rad_var%xi_max(i_phs) = xi_max
       rad_var%xi_tilde = xi / xi_max
       rad_var%y(i_phs) = y
       rad_var%phi = phi
    end associate
    call generator%generate_isr (i_phs, p_born, p_real)
  end subroutine phs_fks_generator_generate_isr_from_xi_and_y

  pure module subroutine phs_fks_generator_get_radiation_variables &
       (generator, i_phs, xi, y, phi)
    class(phs_fks_generator_t), intent(in) :: generator
    integer, intent(in) :: i_phs
    real(default), intent(out) :: xi, y
    real(default), intent(out), optional :: phi
    associate (rad_var => generator%real_kinematics)
       xi = rad_var%xi_max(i_phs) * rad_var%xi_tilde
       y = rad_var%y(i_phs)
       if (present (phi)) phi = rad_var%phi
    end associate
  end subroutine phs_fks_generator_get_radiation_variables

  module subroutine phs_fks_generator_write (generator, unit)
    class(phs_fks_generator_t), intent(in) :: generator
    integer, intent(in), optional :: unit
    integer :: u
    type(string_t) :: massive_phsp
    u = given_output_unit (unit); if (u < 0) return
    if (generator%massive_phsp) then
       massive_phsp = " massive "
    else
       massive_phsp = " massless "
    end if
    write (u, "(A)") char ("This is a generator for a" &
         // massive_phsp // "phase space")
    if (associated (generator%real_kinematics)) then
       call generator%real_kinematics%write ()
    else
       write (u, "(A)") "Warning: There are no real " // &
            "kinematics associated with this generator"
    end if
    call write_separator (u)
    write (u, "(A," // FMT_17 // ",1X)") "sqrts     : ", generator%sqrts
    write (u, "(A," // FMT_17 // ",1X)") "E_gluon   : ", generator%E_gluon
    write (u, "(A," // FMT_17 // ",1X)") "mrec2     : ", generator%mrec2
  end subroutine phs_fks_generator_write

  module subroutine phs_fks_compute_isr_kinematics (phs, r)
    class(phs_fks_t), intent(inout) :: phs
    real(default), intent(in) :: r
    if (.not. phs%config%lab_is_cm) then
       call phs%generator%compute_isr_kinematics &
            (r, phs%lt_cm_to_lab * phs%phs_wood_t%p)
    else
       call phs%generator%compute_isr_kinematics (r, phs%phs_wood_t%p)
    end if
  end subroutine phs_fks_compute_isr_kinematics

  module subroutine phs_fks_final (object)
    class(phs_fks_t), intent(inout) :: object
    call object%forest%final ()
    call object%generator%final ()
  end subroutine phs_fks_final

  subroutine filter_particles_from_resonances &
       (res_hist, exclusion_list, &
    model, res_hist_filtered)
    type(resonance_history_t), intent(in), dimension(:) :: res_hist
    type(string_t), intent(in), dimension(:) :: exclusion_list
    type(model_t), intent(in) :: model
    type(resonance_history_t), intent(out), dimension(:), allocatable :: &
         res_hist_filtered
    integer :: i_hist, i_flv, i_new, n_orig
    logical, dimension(size (res_hist)) :: to_filter
    type(flavor_t) :: flv
    to_filter = .false.
    n_orig = size (res_hist)
    do i_flv = 1, size (exclusion_list)
       call flv%init (exclusion_list (i_flv), model)
       do i_hist = 1, size (res_hist)
          if (res_hist(i_hist)%has_flavor (flv)) to_filter (i_hist) = .true.
       end do
    end do
    allocate (res_hist_filtered (n_orig - count (to_filter)))
    i_new = 1
    do i_hist = 1, size (res_hist)
       if (.not. to_filter (i_hist)) then
          res_hist_filtered (i_new) = res_hist (i_hist)
          i_new = i_new + 1
       end if
    end do
  end subroutine filter_particles_from_resonances

  subroutine clean_resonance_histories &
       (res_hist, n_in, flv, res_hist_clean, success)
    type(resonance_history_t), intent(in), dimension(:) :: res_hist
    integer, intent(in) :: n_in
    integer, intent(in), dimension(:) :: flv
    type(resonance_history_t), intent(out), dimension(:), allocatable :: &
         res_hist_clean
    logical, intent(out) :: success
    integer :: i_hist
    type(resonance_history_t), dimension(:), allocatable :: &
         res_hist_colored, res_hist_contracted

    if (debug_on) call msg_debug (D_SUBTRACTION, "resonance_mapping_init")
    if (debug_active (D_SUBTRACTION)) then
       call msg_debug (D_SUBTRACTION, "Original resonances:")
       do i_hist = 1, size(res_hist)
          call res_hist(i_hist)%write ()
       end do
    end if

    call remove_uncolored_resonances ()
    call contract_resonances (res_hist_colored, res_hist_contracted)
    call remove_subresonances (res_hist_contracted, res_hist_clean)
    !!! Here, we are still not sure whether we actually would rather use
    !!! call remove_multiple_resonances (res_hist_contracted, res_hist_clean)
    if (debug_active (D_SUBTRACTION)) then
       call msg_debug (D_SUBTRACTION, "Resonances after removing uncolored and duplicates: ")
       do i_hist = 1, size (res_hist_clean)
          call res_hist_clean(i_hist)%write ()
       end do
    end if
    if (size (res_hist_clean) == 0) then
       call msg_warning ("No resonances found. Proceed in usual FKS mode.")
       success = .false.
    else
       success = .true.
    end if

  contains
    subroutine remove_uncolored_resonances ()
      type(resonance_history_t), dimension(:), allocatable :: res_hist_tmp
      integer :: n_hist, nleg_out, n_removed
      integer :: i_res, i_hist
      n_hist = size (res_hist)
      nleg_out = size (flv) - n_in
      allocate (res_hist_tmp (n_hist))
      allocate (res_hist_colored (n_hist))
      do i_hist = 1, n_hist
         res_hist_tmp(i_hist) = res_hist(i_hist)
         call res_hist_tmp(i_hist)%add_offset (n_in)
         n_removed = 0
         do i_res = 1, res_hist_tmp(i_hist)%n_resonances
            associate (resonance => res_hist_tmp(i_hist)%resonances(i_res - n_removed))
               if (.not. any (is_colored (flv (resonance%contributors%c))) &
                  .or. size (resonance%contributors%c) == nleg_out) then
                     call res_hist_tmp(i_hist)%remove_resonance (i_res - n_removed)
                     n_removed = n_removed + 1
               end if
            end associate
         end do
         if (allocated (res_hist_tmp(i_hist)%resonances)) then
            if (any (res_hist_colored == res_hist_tmp(i_hist))) then
               cycle
            else
               do i_res = 1, res_hist_tmp(i_hist)%n_resonances
                  associate (resonance => res_hist_tmp(i_hist)%resonances(i_res))
                     call res_hist_colored(i_hist)%add_resonance (resonance)
                  end associate
               end do
            end if
         end if
      end do
    end subroutine remove_uncolored_resonances

    subroutine contract_resonances (res_history_in, res_history_out)
      type(resonance_history_t), intent(in), dimension(:) :: res_history_in
      type(resonance_history_t), intent(out), dimension(:), allocatable :: res_history_out
      logical, dimension(:), allocatable :: i_non_zero
      integer :: n_hist_non_zero, n_hist
      integer :: i_hist_new
      n_hist = size (res_history_in); n_hist_non_zero = 0
      allocate (i_non_zero (n_hist))
      i_non_zero = .false.
      do i_hist = 1, n_hist
         if (res_history_in(i_hist)%n_resonances /= 0) then
            n_hist_non_zero = n_hist_non_zero + 1
            i_non_zero(i_hist) = .true.
         end if
      end do
      allocate (res_history_out (n_hist_non_zero))
      i_hist_new = 1
      do i_hist = 1, n_hist
         if (i_non_zero (i_hist)) then
            res_history_out (i_hist_new) = res_history_in (i_hist)
            i_hist_new = i_hist_new + 1
         end if
      end do
    end subroutine contract_resonances

    subroutine remove_subresonances (res_history_in, res_history_out)
      type(resonance_history_t), intent(in), dimension(:) :: res_history_in
      type(resonance_history_t), intent(out), dimension(:), allocatable :: res_history_out
      logical, dimension(:), allocatable :: i_non_sub_res
      integer :: n_hist, n_hist_non_sub_res
      integer :: i_hist1, i_hist2
      logical :: is_not_subres
      n_hist = size (res_history_in); n_hist_non_sub_res = 0
      allocate (i_non_sub_res (n_hist)); i_non_sub_res = .false.
      do i_hist1 = 1, n_hist
         is_not_subres = .true.
         do i_hist2 = 1, n_hist
            if (i_hist1 == i_hist2) cycle
            is_not_subres = is_not_subres .and. &
               .not.(res_history_in(i_hist2) .contains. res_history_in(i_hist1))
         end do
         if (is_not_subres) then
            n_hist_non_sub_res = n_hist_non_sub_res + 1
            i_non_sub_res (i_hist1) = .true.
         end if
      end do

      allocate (res_history_out (n_hist_non_sub_res))
      i_hist2 = 1
      do i_hist1 = 1, n_hist
         if (i_non_sub_res (i_hist1)) then
            res_history_out (i_hist2) = res_history_in (i_hist1)
            i_hist2 = i_hist2 + 1
         end if
      end do
    end subroutine remove_subresonances

    subroutine remove_multiple_resonances (res_history_in, res_history_out)
      type(resonance_history_t), intent(in), dimension(:) :: res_history_in
      type(resonance_history_t), intent(out), dimension(:), allocatable :: res_history_out
      integer :: n_hist, n_hist_single
      logical, dimension(:), allocatable :: i_hist_single
      integer :: i_hist, j
      n_hist = size (res_history_in)
      n_hist_single = 0
      allocate (i_hist_single (n_hist)); i_hist_single = .false.
      do i_hist = 1, n_hist
         if (res_history_in(i_hist)%n_resonances == 1) then
            n_hist_single = n_hist_single + 1
            i_hist_single(i_hist) = .true.
         end if
      end do

      allocate (res_history_out (n_hist_single))
      j = 1
      do i_hist = 1, n_hist
         if (i_hist_single(i_hist)) then
            res_history_out(j) = res_history_in(i_hist)
            j = j + 1
         end if
      end do
    end subroutine remove_multiple_resonances
  end subroutine clean_resonance_histories

  module subroutine get_filtered_resonance_histories &
       (phs_config, n_in, flv_state, model, excluded_resonances, &
        resonance_histories_filtered, success)
    type(phs_fks_config_t), intent(inout) :: phs_config
    integer, intent(in) :: n_in
    integer, intent(in), dimension(:,:), allocatable :: flv_state
    type(model_t), intent(in) :: model
    type(string_t), intent(in), dimension(:), allocatable :: &
         excluded_resonances
    type(resonance_history_t), intent(out), dimension(:), &
       allocatable :: resonance_histories_filtered
    logical, intent(out) :: success
    type(resonance_history_t), dimension(:), allocatable :: resonance_histories
    type(resonance_history_t), dimension(:), allocatable :: &
       resonance_histories_clean!, resonance_histories_filtered
    allocate (resonance_histories (size (phs_config%get_resonance_histories ())))
    resonance_histories = phs_config%get_resonance_histories ()
    call clean_resonance_histories (resonance_histories, &
         n_in, flv_state (:,1), resonance_histories_clean, success)
    if (success .and. allocated (excluded_resonances)) then
       call filter_particles_from_resonances (resonance_histories_clean, &
            excluded_resonances, model, resonance_histories_filtered)
    else
       allocate (resonance_histories_filtered (size (resonance_histories_clean)))
       resonance_histories_filtered = resonance_histories_clean
    end if
  end subroutine get_filtered_resonance_histories


end submodule phs_fks_s

