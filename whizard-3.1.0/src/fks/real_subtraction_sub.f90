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

submodule (real_subtraction) real_subtraction_s

  use debug_master, only: debug_on
  use io_units
  use format_defs, only: FMT_15
  use string_utils
  use constants
  use numeric_utils
  use diagnostics
  use pdg_arrays
  use sm_physics
  use models
  use ttv_formfactors, only: m1s_to_mpole

  implicit none

contains

  function nlo_purpose (purpose) result (of_purpose)
    type(string_t) :: of_purpose
    integer, intent(in) :: purpose
    select case (purpose)
    case (INTEGRATION)
       of_purpose = var_str ("Integration")
    case (FIXED_ORDER_EVENTS)
       of_purpose = var_str ("Fixed order NLO events")
    case default
       of_purpose = var_str ("Undefined!")
    end select
  end function nlo_purpose

  module subroutine soft_subtraction_init (sub_soft, reg_data)
    class(soft_subtraction_t), intent(inout) :: sub_soft
    type(region_data_t), intent(in), target :: reg_data
    sub_soft%reg_data => reg_data
    allocate (sub_soft%momentum_matrix (reg_data%n_legs_born, &
         reg_data%n_legs_born))
  end subroutine soft_subtraction_init

  module function soft_subtraction_requires_boost &
       (sub_soft, sqrts) result (requires_boost)
    logical :: requires_boost
    class(soft_subtraction_t), intent(in) :: sub_soft
    real(default), intent(in) :: sqrts
    real(default) :: mtop
    logical :: above_threshold
    if (sub_soft%factorization_mode == FACTORIZATION_THRESHOLD) then
       mtop = m1s_to_mpole (sqrts)
       above_threshold = sqrts**2 - four * mtop**2 > zero
    else
       above_threshold = .false.
    end if
    requires_boost = sub_soft%use_resonance_mappings .or. above_threshold
  end function soft_subtraction_requires_boost

  module subroutine soft_subtraction_create_softvec_fsr &
       (sub_soft, p_born, y, phi, emitter, xi_ref_momentum)
    class(soft_subtraction_t), intent(inout) :: sub_soft
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: y, phi
    integer, intent(in) :: emitter
    type(vector4_t), intent(in) :: xi_ref_momentum
    type(vector3_t) :: dir
    type(vector4_t) :: p_em
    type(lorentz_transformation_t) :: rot
    type(lorentz_transformation_t) :: boost_to_rest_frame
    logical :: requires_boost
    associate (p_soft => sub_soft%p_soft)
       p_soft%p(0) = one
       requires_boost = sub_soft%requires_boost (two * p_born(1)%p(0))
       if (requires_boost) then
          boost_to_rest_frame = inverse (boost (xi_ref_momentum, xi_ref_momentum**1))
          p_em = boost_to_rest_frame * p_born(emitter)
       else
          p_em = p_born(emitter)
       end if
       p_soft%p(1:3) = p_em%p(1:3) / space_part_norm (p_em)
       dir = create_orthogonal (space_part (p_em))
       rot = rotation (y, sqrt(one - y**2), dir)
       p_soft = rot * p_soft
       if (.not. vanishes (phi)) then
         dir = space_part (p_em) / space_part_norm (p_em)
         rot = rotation (cos(phi), sin(phi), dir)
         p_soft = rot * p_soft
       end if
       if (requires_boost) p_soft = inverse (boost_to_rest_frame) * p_soft
    end associate
  end subroutine soft_subtraction_create_softvec_fsr

  module subroutine soft_subtraction_create_softvec_isr (sub_soft, y, phi)
    class(soft_subtraction_t), intent(inout) :: sub_soft
    real(default), intent(in) :: y, phi
    real(default) :: sin_theta
    sin_theta = sqrt(one - y**2)
    associate (p => sub_soft%p_soft%p)
       p(0) = one
       p(1) = sin_theta * sin(phi)
       p(2) = sin_theta * cos(phi)
       p(3) = y
    end associate
  end subroutine soft_subtraction_create_softvec_isr

  module subroutine soft_subtraction_create_softvec_mismatch &
       (sub_soft, E, y, phi, p_em)
    class(soft_subtraction_t), intent(inout) :: sub_soft
    real(default), intent(in) :: E, phi, y
    type(vector4_t), intent(in) :: p_em
    real(default) :: sin_theta
    type(lorentz_transformation_t) :: rot_em_off_3_axis
    sin_theta = sqrt (one - y**2)
    associate (p => sub_soft%p_soft%p)
       p(0) = E
       p(1) = E * sin_theta * sin(phi)
       p(2) = E * sin_theta * cos(phi)
       p(3) = E * y
    end associate
    rot_em_off_3_axis = rotation_to_2nd (3, space_part (p_em))
    sub_soft%p_soft = rot_em_off_3_axis * sub_soft%p_soft
  end subroutine soft_subtraction_create_softvec_mismatch

  module function soft_subtraction_compute (sub_soft, p_born, &
       born_ij, y, q2, alpha_coupling, alr, emitter, i_res) result (sqme)
    real(default) :: sqme
    class(soft_subtraction_t), intent(inout) :: sub_soft
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in), dimension(:,:) :: born_ij
    real(default), intent(in) :: y
    real(default), intent(in) :: q2, alpha_coupling
    integer, intent(in) :: alr, emitter, i_res
    real(default) :: s_alpha_soft
    real(default) :: kb
    real(default) :: xi2_factor

    if (.not. vector_set_is_cms (p_born, sub_soft%reg_data%n_in)) then
       call vector4_write_set (p_born, show_mass = .true., &
              check_conservation = .true.)
       call msg_fatal ("Soft subtraction: phase space point must be in CMS")
    end if
    if (debug2_active (D_SUBTRACTION)) then
       select case (char (sub_soft%reg_data%regions(alr)%nlo_correction_type))
       case ("QCD")
          print *, 'Compute soft subtraction using alpha_s = ', alpha_coupling
       case ("EW")
          print *, 'Compute soft subtraction using alpha_qed = ', alpha_coupling
       end select
    end if

    s_alpha_soft = sub_soft%reg_data%get_svalue_soft (p_born, &
         sub_soft%p_soft, alr, emitter, i_res)
    if (s_alpha_soft > one + tiny_07) call msg_fatal ("s_alpha_soft > 1!")
    if (debug2_active (D_SUBTRACTION)) &
         call msg_print_color ('s_alpha_soft', s_alpha_soft, COL_YELLOW)
    select case (sub_soft%factorization_mode)
    case (NO_FACTORIZATION)
       kb = sub_soft%evaluate_factorization_default (p_born, born_ij)
    case (FACTORIZATION_THRESHOLD)
       kb = sub_soft%evaluate_factorization_threshold (thr_leg(emitter), p_born, born_ij)
    end select
    if (debug_on) call msg_debug2 (D_SUBTRACTION, 'KB', kb)
    sqme = four * pi * alpha_coupling * s_alpha_soft * kb
    if (sub_soft%xi2_expanded) then
       xi2_factor = four / q2
    else
       xi2_factor = one
    end if
    if (emitter <= sub_soft%reg_data%n_in) then
       sqme = xi2_factor * (one - y**2) * sqme
    else
       sqme = xi2_factor * (one - y) * sqme
    end if
    if (sub_soft%reg_data%regions(alr)%double_fsr) sqme = sqme * two
  end function soft_subtraction_compute

  module function soft_subtraction_evaluate_factorization_default &
       (sub_soft, p, born_ij) result (kb)
    real(default) :: kb
    class(soft_subtraction_t), intent(inout) :: sub_soft
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(in), dimension(:,:) :: born_ij
    integer :: i, j
    kb = zero
    call sub_soft%compute_momentum_matrix (p)
    do i = 1, size (p)
       do j = 1, size (p)
          kb = kb + sub_soft%momentum_matrix (i, j) * born_ij (i, j)
       end do
    end do
  end function soft_subtraction_evaluate_factorization_default

  module subroutine soft_subtraction_compute_momentum_matrix &
       (sub_soft, p_born)
    class(soft_subtraction_t), intent(inout) :: sub_soft
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default) :: num, deno1, deno2
    integer :: i, j
    do i = 1, sub_soft%reg_data%n_legs_born
       do j = 1, sub_soft%reg_data%n_legs_born
          if (i <= j) then
             num = p_born(i) * p_born(j)
             deno1 = p_born(i) * sub_soft%p_soft
             deno2 = p_born(j) * sub_soft%p_soft
             sub_soft%momentum_matrix(i, j) = num / (deno1 * deno2)
          else
             !!! momentum matrix is symmetric.
            sub_soft%momentum_matrix(i, j) = sub_soft%momentum_matrix(j, i)
          end if
       end do
    end do
  end subroutine soft_subtraction_compute_momentum_matrix

  module function soft_subtraction_evaluate_factorization_threshold &
       (sub_soft, leg, p_born, born_ij) result (kb)
    real(default) :: kb
    class(soft_subtraction_t), intent(inout) :: sub_soft
    integer, intent(in) :: leg
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in), dimension(:,:) :: born_ij
    type(vector4_t), dimension(4) :: p
    p = get_threshold_momenta (p_born)
    kb = evaluate_leg_pair (ASSOCIATED_LEG_PAIR (leg))
    if (debug2_active (D_SUBTRACTION))  call show_debug ()

  contains

    function evaluate_leg_pair (i_start) result (kbb)
      real(default) :: kbb
      integer, intent(in) :: i_start
      integer :: i1, i2
      real(default) :: numerator, deno1, deno2
      kbb = zero
      do i1 = i_start, i_start + 1
         do i2 = i_start, i_start + 1
            numerator = p(i1) * p(i2)
            deno1 = p(i1) * sub_soft%p_soft
            deno2 = p(i2) * sub_soft%p_soft
            kbb = kbb +  numerator * born_ij (i1, i2) / deno1 / deno2
         end do
      end do
      if (debug2_active (D_SUBTRACTION)) then
         do i1 = i_start, i_start + 1
            do i2 = i_start, i_start + 1
               call msg_print_color('i1', i1, COL_PEACH)
               call msg_print_color('i2', i2, COL_PEACH)
               call msg_print_color('born_ij (i1,i2)', born_ij (i1,i2), &
                    COL_PINK)
               print *, 'Top momentum: ', p(1)%p
            end do
         end do
      end if
    end function evaluate_leg_pair

    subroutine show_debug ()
      integer :: i
      call msg_print_color &
           ('soft_subtraction_evaluate_factorization_threshold', COL_GREEN)
      do i = 1, 4
         print *, 'sqrt(p(i)**2) =    ', sqrt(p(i)**2)
      end do
    end subroutine show_debug

  end function soft_subtraction_evaluate_factorization_threshold

  module function soft_subtraction_i_xi_ref &
       (sub_soft, alr, i_phs) result (i_xi_ref)
    integer :: i_xi_ref
    class(soft_subtraction_t), intent(in) :: sub_soft
    integer, intent(in) :: alr, i_phs
    if (sub_soft%use_resonance_mappings) then
       i_xi_ref = sub_soft%reg_data%alr_to_i_contributor (alr)
    else if (sub_soft%factorization_mode == FACTORIZATION_THRESHOLD) then
       i_xi_ref = i_phs
    else
       i_xi_ref = 1
    end if
  end function soft_subtraction_i_xi_ref

  module subroutine soft_subtraction_final (sub_soft)
    class(soft_subtraction_t), intent(inout) :: sub_soft
    if (associated (sub_soft%reg_data)) nullify (sub_soft%reg_data)
    if (allocated (sub_soft%momentum_matrix)) &
         deallocate (sub_soft%momentum_matrix)
  end subroutine soft_subtraction_final

  module subroutine soft_mismatch_init (soft_mismatch, reg_data, &
       real_kinematics, factorization_mode)
    class(soft_mismatch_t), intent(inout) :: soft_mismatch
    type(region_data_t), intent(in), target :: reg_data
    type(real_kinematics_t), intent(in), target :: real_kinematics
    integer, intent(in) :: factorization_mode
    soft_mismatch%reg_data => reg_data
    allocate (soft_mismatch%sqme_born (reg_data%n_flv_born))
    allocate (soft_mismatch%sqme_born_color_c (reg_data%n_legs_born, &
         reg_data%n_legs_born, reg_data%n_flv_born))
    allocate (soft_mismatch%sqme_born_charge_c (reg_data%n_legs_born, &
         reg_data%n_legs_born, reg_data%n_flv_born))
    call soft_mismatch%sub_soft%init (reg_data)
    soft_mismatch%sub_soft%xi2_expanded = .false.
    soft_mismatch%real_kinematics => real_kinematics
    soft_mismatch%sub_soft%factorization_mode = factorization_mode
  end subroutine soft_mismatch_init

  module function soft_mismatch_evaluate &
       (soft_mismatch, alpha_s) result (sqme_mismatch)
    real(default) :: sqme_mismatch
    class(soft_mismatch_t), intent(inout) :: soft_mismatch
    real(default), intent(in) :: alpha_s
    integer :: alr, i_born, emitter, i_res, i_phs, i_con
    real(default) :: xi, y, q2, s
    real(default) :: E_gluon
    type(vector4_t) :: p_em
    real(default) :: sqme_alr, sqme_soft
    type(vector4_t), dimension(:), allocatable :: p_born
    sqme_mismatch = zero
    associate (real_kinematics => soft_mismatch%real_kinematics)
       xi = real_kinematics%xi_mismatch
       y = real_kinematics%y_mismatch
       s = real_kinematics%cms_energy2
       E_gluon = sqrt (s) * xi / two

       if (debug_active (D_MISMATCH)) then
          print *, 'Evaluating soft mismatch: '
          print *, 'Phase space: '
          call vector4_write_set (real_kinematics%p_born_cms%get_momenta(1), &
               show_mass = .true.)
          print *, 'xi: ', xi, 'y: ', y, 's: ', s, 'E_gluon: ', E_gluon
       end if

       allocate (p_born (soft_mismatch%reg_data%n_legs_born))

       do alr = 1, soft_mismatch%reg_data%n_regions

           i_phs = real_kinematics%alr_to_i_phs (alr)
           if (soft_mismatch%reg_data%has_pseudo_isr ()) then
              i_con = 1
              p_born = &
                   soft_mismatch%real_kinematics%p_born_onshell%get_momenta(1)
           else
              i_con = soft_mismatch%reg_data%alr_to_i_contributor (alr)
              p_born = soft_mismatch%real_kinematics%p_born_cms%get_momenta(1)
           end if
           q2 = real_kinematics%xi_ref_momenta(i_con)**2
           emitter = soft_mismatch%reg_data%regions(alr)%emitter
           p_em = p_born (emitter)
           i_res = soft_mismatch%reg_data%regions(alr)%i_res
           i_born = soft_mismatch%reg_data%regions(alr)%uborn_index

           call print_debug_alr ()

           call soft_mismatch%sub_soft%create_softvec_mismatch &
                (E_gluon, y, real_kinematics%phi, p_em)
           if (debug_active (D_MISMATCH)) &
                print *, 'Created soft vector: ', &
                soft_mismatch%sub_soft%p_soft%p

           select type (fks_mapping => soft_mismatch%reg_data%fks_mapping)
           type is (fks_mapping_resonances_t)
              call fks_mapping%set_resonance_momentum &
                   (real_kinematics%xi_ref_momenta(i_con))
           end select

           sqme_soft = soft_mismatch%sub_soft%compute &
                (p_born, soft_mismatch%sqme_born_color_c(:,:,i_born), y, &
                q2, alpha_s, alr, emitter, i_res)

           sqme_alr = soft_mismatch%compute (alr, xi, y, p_em, &
                real_kinematics%xi_ref_momenta(i_con), &
                soft_mismatch%sub_soft%p_soft, &
                soft_mismatch%sqme_born(i_born), sqme_soft, &
                alpha_s, s)

           if (debug_on) call msg_debug (D_MISMATCH, 'sqme_alr: ', sqme_alr)
           sqme_mismatch = sqme_mismatch + sqme_alr

       end do
    end associate
  contains
    subroutine print_debug_alr ()
      if (debug_active (D_MISMATCH)) then
         print *, 'alr: ', alr
         print *, 'i_phs: ', i_phs, 'i_con: ', i_con, 'i_res: ', i_res
         print *, 'emitter: ', emitter, 'i_born: ', i_born
         print *, 'emitter momentum: ', p_em%p
         print *, 'resonance momentum: ', &
            soft_mismatch%real_kinematics%xi_ref_momenta(i_con)%p
         print *, 'q2: ', q2
      end if
    end subroutine print_debug_alr
  end function soft_mismatch_evaluate

  module function soft_mismatch_compute &
       (soft_mismatch, alr, xi, y, p_em, p_res, p_soft, &
     sqme_born, sqme_soft, alpha_s, s) result (sqme_mismatch)
    real(default) :: sqme_mismatch
    class(soft_mismatch_t), intent(in) :: soft_mismatch
    integer, intent(in) :: alr
    real(default), intent(in) :: xi, y
    type(vector4_t), intent(in) :: p_em, p_res, p_soft
    real(default), intent(in) :: sqme_born, sqme_soft
    real(default), intent(in) :: alpha_s, s
    real(default) :: q2, expo, sm1, sm2, jacobian

    q2 = p_res**2
    expo = - two * p_soft * p_res / q2
    !!! Divide by 1 - y to factor out the corresponding
    !!! factor in the soft matrix element
    sm1 = sqme_soft / (one - y) * ( exp(expo) - exp(- xi) )
    if (debug_on)  call msg_debug2 &
         (D_MISMATCH, 'sqme_soft in mismatch ', sqme_soft)

    sm2 = zero
    if (soft_mismatch%reg_data%regions(alr)%has_collinear_divergence ()) then
       expo = - two * p_em * p_res / q2 * &
          p_soft%p(0) / p_em%p(0)
       sm2 = 32 * pi * alpha_s * cf / (s * xi**2) * sqme_born * &
          ( exp(expo) - exp(- xi) ) / (one - y)
    end if

    jacobian = soft_mismatch%real_kinematics%jac_mismatch * &
         s * xi / (8 * twopi3)
    sqme_mismatch = (sm1 - sm2) * jacobian

  end function soft_mismatch_compute

  module subroutine soft_mismatch_final (soft_mismatch)
    class(soft_mismatch_t), intent(inout) :: soft_mismatch
    call soft_mismatch%sub_soft%final ()
    if (associated (soft_mismatch%reg_data)) nullify (soft_mismatch%reg_data)
    if (allocated (soft_mismatch%sqme_born)) &
         deallocate (soft_mismatch%sqme_born)
    if (allocated (soft_mismatch%sqme_born_color_c)) &
         deallocate (soft_mismatch%sqme_born_color_c)
    if (allocated (soft_mismatch%sqme_born_charge_c)) &
         deallocate (soft_mismatch%sqme_born_charge_c)
    if (associated (soft_mismatch%real_kinematics)) &
         nullify (soft_mismatch%real_kinematics)
  end subroutine soft_mismatch_final

  module subroutine coll_subtraction_init (coll_sub, n_alr, n_in)
    class(coll_subtraction_t), intent(inout) :: coll_sub
    integer, intent(in) :: n_alr, n_in
    coll_sub%n_in = n_in
    coll_sub%n_alr = n_alr
  end subroutine coll_subtraction_init

  module subroutine coll_subtraction_set_parameters (coll_sub, CA, CF, TR)
    class(coll_subtraction_t), intent(inout) :: coll_sub
    real(default), intent(in) :: CA, CF, TR
    coll_sub%CA = CA
    coll_sub%CF = CF
    coll_sub%TR = TR
  end subroutine coll_subtraction_set_parameters

  module function coll_subtraction_compute_fsr (coll_sub, emitter, &
       flst, p_res, p_born, sqme_born, mom_times_sqme_spin_c, &
       xi, alpha_coupling, double_fsr) result (sqme)
    real(default) :: sqme
    class(coll_subtraction_t), intent(in) :: coll_sub
    integer, intent(in) :: emitter
    integer, dimension(:), intent(in) :: flst
    type(vector4_t), intent(in) :: p_res
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: sqme_born, mom_times_sqme_spin_c
    real(default), intent(in) :: xi, alpha_coupling
    logical, intent(in) :: double_fsr
    real(default) :: q0, z, p0, z_o_xi, onemz
    integer :: nlegs, flv_em, flv_rad
    nlegs = size (flst)
    flv_rad = flst(nlegs); flv_em = flst(emitter)
    q0 = p_res**1
    p0 = p_res * p_born(emitter) / q0
    !!! Here, z corresponds to 1-z in the formulas of arXiv:1002.2581;
    !!! the integrand is symmetric under this variable change
    z_o_xi = q0 / (two * p0)
    z = xi * z_o_xi; onemz = one - z
    if (is_gluon (flv_em) .and. is_gluon (flv_rad)) then
       sqme = coll_sub%CA * ( two * ( z / onemz * xi + onemz / z_o_xi ) * sqme_born &
            + four * xi * z * onemz * mom_times_sqme_spin_c )
    else if (is_fermion (flv_em) .and. is_fermion (flv_rad)) then
       sqme = coll_sub%TR * xi * (sqme_born - four * z * onemz * mom_times_sqme_spin_c)
    else if (is_fermion (flv_em) .and. is_massless_vector (flv_rad)) then
       sqme = sqme_born * coll_sub%CF * (one + onemz**2) / z_o_xi
    else
       sqme = zero
    end if
    sqme = sqme / (p0**2 * onemz * z_o_xi)
    sqme = sqme * four * pi * alpha_coupling
    if (double_fsr) sqme = sqme * onemz * two
  end function coll_subtraction_compute_fsr

  module function coll_subtraction_compute_isr &
     (coll_sub, emitter, flst, p_born, sqme_born, mom_times_sqme_spin_c, &
     xi, alpha_coupling, isr_mode) result (sqme)
    real(default) :: sqme
    class(coll_subtraction_t), intent(in) :: coll_sub
    integer, intent(in) :: emitter
    integer, dimension(:), intent(in) :: flst
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: sqme_born
    real(default), intent(in) :: mom_times_sqme_spin_c
    real(default), intent(in) :: xi, alpha_coupling
    integer, intent(in) :: isr_mode
    real(default) :: z, onemz, p02
    integer :: nlegs, flv_em, flv_rad
    !!! p_born must be in lab frame.
    nlegs = size (flst)
    flv_rad = flst(nlegs); flv_em = flst(emitter)
    !!! No need to pay attention to n_in = 1, because this case always has a
    !!! massive initial-state particle and thus no collinear divergence.
    p02 = p_born(1)%p(0) * p_born(2)%p(0) / two
    z = one - xi; onemz = xi
    if (is_massless_vector (flv_em) .and. is_massless_vector (flv_rad)) then
       sqme = coll_sub%CA * (two * (z + z * onemz**2) * sqme_born + four * onemz**2 &
          / z * mom_times_sqme_spin_c)
    else if (is_fermion (flv_em) .and. is_massless_vector (flv_rad)) then
       sqme = coll_sub%CF * (one + z**2) * sqme_born
    else if (is_fermion (flv_em) .and. is_fermion (flv_rad)) then
       sqme = coll_sub%CF * (z * onemz * sqme_born + four * onemz**2 / z * mom_times_sqme_spin_c)
    else if (is_massless_vector (flv_em) .and. is_fermion (flv_rad)) then
       sqme = coll_sub%TR * (z**2 + onemz**2) * onemz * sqme_born
    else
       sqme = zero
    end if
    if (isr_mode == SQRTS_VAR) then
       sqme = sqme / p02 * z
    else
       !!! We have no idea why this seems to work as there should be no factor
       !!! of z for the fixed-beam settings. This should definitely be understood in the
       !!! future!
       sqme = sqme / p02 / z
    end if
    sqme = sqme * four * pi * alpha_coupling
  end function coll_subtraction_compute_isr

  module subroutine coll_subtraction_final (sub_coll)
    class(coll_subtraction_t), intent(inout) :: sub_coll
    sub_coll%use_resonance_mappings = .false.
  end subroutine coll_subtraction_final

  module subroutine real_subtraction_init (rsub, reg_data, settings)
    class(real_subtraction_t), intent(inout), target :: rsub
    type(region_data_t), intent(in), target :: reg_data
    type(nlo_settings_t), intent(in), target :: settings
    integer :: alr
    if (debug_on) call msg_debug (D_SUBTRACTION, "real_subtraction_init")
    if (debug_on) call msg_debug (D_SUBTRACTION, "n_in", reg_data%n_in)
    if (debug_on) call msg_debug &
         (D_SUBTRACTION, "nlegs_born", reg_data%n_legs_born)
    if (debug_on) call msg_debug &
         (D_SUBTRACTION, "nlegs_real", reg_data%n_legs_real)
    if (debug_on) call msg_debug &
         (D_SUBTRACTION, "reg_data%n_regions", reg_data%n_regions)
    if (debug2_active (D_SUBTRACTION))  call reg_data%write ()
    rsub%reg_data => reg_data
    allocate (rsub%sqme_born (reg_data%n_flv_born))
    rsub%sqme_born = zero
    allocate (rsub%sf_factors (reg_data%n_regions, 0:reg_data%n_in))
    rsub%sf_factors = zero
    allocate (rsub%sqme_real_arr (reg_data%n_regions))
    rsub%sqme_real_arr = zero
    allocate (rsub%sqme_born_color_c &
         (reg_data%n_legs_born, reg_data%n_legs_born, reg_data%n_flv_born))
    rsub%sqme_born_color_c = zero
    allocate (rsub%sqme_born_charge_c &
         (reg_data%n_legs_born, reg_data%n_legs_born, reg_data%n_flv_born))
    rsub%sqme_born_charge_c = zero
    allocate (rsub%sqme_real_non_sub (reg_data%n_flv_real, reg_data%n_phs))
    rsub%sqme_real_non_sub = zero
    allocate (rsub%sc_required (reg_data%n_regions))
    do alr = 1, reg_data%n_regions
       rsub%sc_required(alr) = reg_data%regions(alr)%sc_required
    end do
    if (rsub%requires_spin_correlations ()) then
       allocate (rsub%sqme_born_spin_c &
            (1:3, 1:3, reg_data%n_legs_born, reg_data%n_flv_born))
       rsub%sqme_born_spin_c = zero
    end if
    call rsub%sub_soft%init (reg_data)
    call rsub%sub_coll%init (reg_data%n_regions, reg_data%n_in)
    rsub%settings => settings
    rsub%sub_soft%use_resonance_mappings = settings%use_resonance_mappings
    rsub%sub_coll%use_resonance_mappings = settings%use_resonance_mappings
    rsub%sub_soft%factorization_mode = settings%factorization_mode
  end subroutine real_subtraction_init

  module subroutine real_subtraction_set_real_kinematics (rsub, real_kinematics)
    class(real_subtraction_t), intent(inout) :: rsub
    type(real_kinematics_t), intent(in), target :: real_kinematics
    rsub%real_kinematics => real_kinematics
  end subroutine real_subtraction_set_real_kinematics

  module subroutine real_subtraction_set_isr_kinematics (rsub, fractions)
    class(real_subtraction_t), intent(inout) :: rsub
    type(isr_kinematics_t), intent(in), target :: fractions
    rsub%isr_kinematics => fractions
  end subroutine real_subtraction_set_isr_kinematics

  module function real_subtraction_get_i_res (rsub, alr) result (i_res)
    integer :: i_res
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr
    select type (fks_mapping => rsub%reg_data%fks_mapping)
    type is (fks_mapping_resonances_t)
       i_res = fks_mapping%res_map%alr_to_i_res (alr)
    class default
       i_res = 0
    end select
  end function real_subtraction_get_i_res

  module subroutine real_subtraction_compute (rsub, emitter, &
       i_phs, alpha_s, alpha_qed, separate_alrs, sqme)
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: emitter, i_phs
    logical, intent(in) :: separate_alrs
    real(default), intent(inout), dimension(:) :: sqme
    real(default), dimension(:), allocatable :: sqme_alr_arr
    real(default), intent(in) :: alpha_s, alpha_qed
    real(default) :: sqme_alr, alpha_coupling
    integer :: alr, i_con, i_res, this_emitter
    logical :: same_emitter
    logical, dimension(:), allocatable :: alr_evaluated
    allocate (alr_evaluated(rsub%reg_data%n_regions))
    allocate (sqme_alr_arr(rsub%reg_data%n_regions))
    sqme_alr_arr = zero
    alr_evaluated = .false.
    do alr = 1, rsub%reg_data%n_regions
       if (.not. alr_evaluated(rsub%reg_data%regions(alr)%eqv_index)) then
          if (allocated (rsub%selected_alr)) then
             if (.not. any (rsub%selected_alr == alr)) cycle
          end if
          sqme_alr = zero
          if (emitter > rsub%isr_kinematics%n_in) then
             same_emitter = emitter == rsub%reg_data%regions(alr)%emitter
          else
             same_emitter = rsub%reg_data%regions(alr)%emitter <= rsub%isr_kinematics%n_in
          end if
          select case (char(rsub%reg_data%regions(alr)%nlo_correction_type))
          case ("QCD")
             alpha_coupling = alpha_s
          case ("EW")
             alpha_coupling = alpha_qed
          end select
          if (same_emitter .and. i_phs == rsub%real_kinematics%alr_to_i_phs (alr)) then
             i_res = rsub%get_i_res (alr)
             this_emitter = rsub%reg_data%regions(alr)%emitter
             sqme_alr = rsub%evaluate_emitter_region (alr, this_emitter, i_phs, &
                  i_res, alpha_coupling)
             i_con = rsub%get_i_contributor (alr)
             sqme_alr = sqme_alr * rsub%get_phs_factor (i_con)
          end if
          sqme_alr_arr(alr) = sqme_alr_arr(alr) + sqme_alr
          if (.not. (debug_active (D_SUBTRACTION) .or. debug2_active (D_SUBTRACTION))) then
             if (.not. allocated (rsub%selected_alr)) &
                  alr_evaluated(rsub%reg_data%regions(alr)%eqv_index) = .true.
          end if
       else
          sqme_alr_arr(alr) = sqme_alr_arr(rsub%reg_data%regions(alr)%eqv_index)
       end if
       if (rsub%radiation_event .and. sqme_alr_arr(alr) /= zero) then
          rsub%sqme_real_arr(alr) = sqme_alr_arr(alr)
       end if
       if (separate_alrs) then
          sqme(alr) = sqme(alr) + sqme_alr_arr(alr)
       else
          sqme(1) = sqme(1) + sqme_alr_arr(alr)
       end if
    end do
    if (debug_on) then
       if (debug2_active (D_SUBTRACTION)) call check_s_alpha_consistency ()
    end if
  contains
    subroutine check_s_alpha_consistency ()
      real(default) :: sum_s_alpha, sum_s_alpha_soft
      integer :: i_ftuple
      if (debug_on) call msg_debug2 (D_SUBTRACTION, "Check consistency of s_alpha: ")
      do alr = 1, rsub%reg_data%n_regions
         sum_s_alpha = rsub%sum_up_s_alpha(alr, i_phs)
         call msg_debug2 (D_SUBTRACTION, 'sum_s_alpha', sum_s_alpha)
         if (.not. nearly_equal(sum_s_alpha, one)) then
            call msg_bug ("The sum of all S functions should be equal to one!")
         end if
         sum_s_alpha_soft = rsub%sum_up_s_alpha_soft(alr, i_phs)
         call msg_debug2 (D_SUBTRACTION, 'sum_s_alpha_soft', sum_s_alpha_soft)
         if (.not. nearly_equal(sum_s_alpha_soft, one)) then
            call msg_bug ("The sum of all soft S functions should be equal to one!")
         end if
      end do
    end subroutine check_s_alpha_consistency
  end subroutine real_subtraction_compute

  module function real_subtraction_evaluate_emitter_region (rsub, alr, &
       emitter, i_phs, i_res, alpha_coupling) result (sqme)
    real(default) :: sqme
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, emitter, i_phs, i_res
    real(default), intent(in) :: alpha_coupling
    if (emitter <= rsub%isr_kinematics%n_in) then
       sqme = rsub%evaluate_region_isr &
            (alr, emitter, i_phs, i_res, alpha_coupling)
    else
       select type (fks_mapping => rsub%reg_data%fks_mapping)
       type is (fks_mapping_resonances_t)
          call fks_mapping%set_resonance_momenta &
               (rsub%real_kinematics%xi_ref_momenta)
       end select
       sqme = rsub%evaluate_region_fsr (alr, emitter, i_phs, i_res, alpha_coupling)
    end if
  end function real_subtraction_evaluate_emitter_region

  module function real_subtraction_sum_up_s_alpha &
       (rsub, alr, i_phs) result (sum_s_alpha)
    real(default) :: sum_s_alpha
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, i_phs
    type(vector4_t), dimension(:), allocatable :: p_real
    integer :: i_res, i_ftuple, i1, i2
    allocate (p_real (rsub%reg_data%n_legs_real))
    if (rsub%reg_data%has_pseudo_isr ()) then
       p_real = rsub%real_kinematics%p_real_onshell(i_phs)%get_momenta (i_phs)
    else
       p_real = rsub%real_kinematics%p_real_cms%get_momenta (i_phs)
    end if
    i_res = rsub%get_i_res (alr)
    sum_s_alpha = zero
    do i_ftuple = 1, rsub%reg_data%regions(alr)%nregions
       call rsub%reg_data%regions(alr)%ftuples(i_ftuple)%get (i1, i2)
       sum_s_alpha = sum_s_alpha + rsub%reg_data%get_svalue (p_real, alr, i1, i2, i_res)
    end do
  end function real_subtraction_sum_up_s_alpha

  module function real_subtraction_sum_up_s_alpha_soft &
       (rsub, alr, i_phs) result (sum_s_alpha_soft)
    real(default) :: sum_s_alpha_soft
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, i_phs
    type(vector4_t), dimension(:), allocatable :: p_born
    integer :: i_res, i_ftuple, i1, i2, emitter, nlegs
    allocate (p_born (rsub%reg_data%n_legs_born))
    if (rsub%reg_data%has_pseudo_isr ()) then
       p_born = rsub%real_kinematics%p_born_onshell%get_momenta (1)
    else
       p_born = rsub%real_kinematics%p_born_cms%get_momenta (1)
    end if
    i_res = rsub%get_i_res (alr)
    emitter = rsub%reg_data%regions(alr)%emitter
    associate (r => rsub%real_kinematics)
       if (emitter > rsub%sub_soft%reg_data%n_in) then
          call rsub%sub_soft%create_softvec_fsr (p_born, r%y_soft(i_phs), r%phi, &
               emitter, r%xi_ref_momenta(rsub%sub_soft%i_xi_ref (alr, i_phs)))
       else
          call rsub%sub_soft%create_softvec_isr (r%y_soft(i_phs), r%phi)
       end if
    end associate
    nlegs = rsub%reg_data%n_legs_real
    sum_s_alpha_soft = zero
    do i_ftuple = 1, rsub%reg_data%regions(alr)%nregions
       call rsub%reg_data%regions(alr)%ftuples(i_ftuple)%get (i1, i2)
       if (i2 == nlegs) then
          sum_s_alpha_soft = sum_s_alpha_soft + rsub%reg_data%get_svalue_soft &
            (p_born, rsub%sub_soft%p_soft, alr, i1, i_res)
       end if
    end do
  end function real_subtraction_sum_up_s_alpha_soft

  module function real_subtraction_evaluate_region_fsr (rsub, alr, &
       emitter, i_phs, i_res, alpha_coupling) result (sqme_tot)
    real(default) :: sqme_tot
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, emitter, i_phs, i_res
    real(default), intent(in) :: alpha_coupling
    real(default) :: sqme_rad, sqme_soft, sqme_coll, sqme_cs, sqme_remn
    sqme_rad = zero; sqme_soft = zero; sqme_coll = zero
    sqme_cs = zero; sqme_remn = zero
    associate (region => rsub%reg_data%regions(alr), &
         template => rsub%settings%fks_template)
      if (rsub%radiation_event) then
         sqme_rad = rsub%sqme_real_non_sub &
                     (rsub%reg_data%get_matrix_element_index (alr), i_phs)
         call evaluate_fks_factors (sqme_rad, rsub%reg_data, &
              rsub%real_kinematics, alr, i_phs, emitter, i_res)
         call apply_kinematic_factors_radiation (sqme_rad, &
              rsub%real_kinematics, i_phs)
      end if
      if (rsub%subtraction_event .and. .not. (rsub%subtraction_deactivated &
            .or. region%nlo_correction_type == "none")) then
         if (debug2_active (D_SUBTRACTION)) then
            print *, "[real_subtraction_evaluate_region_fsr]"
            print *, "xi: ", rsub%real_kinematics%xi_max(i_phs) * &
                 rsub%real_kinematics%xi_tilde
            print *, "y: ", rsub%real_kinematics%y(i_phs)
         end if
         call rsub%evaluate_subtraction_terms_fsr (alr, emitter, i_phs, &
              i_res, alpha_coupling, sqme_soft, sqme_coll, sqme_cs)
         call apply_kinematic_factors_subtraction_fsr (sqme_soft, &
              sqme_coll, sqme_cs, rsub%real_kinematics, i_phs)
         associate (symm_factor_fs => &
              rsub%reg_data%born_to_real_symm_factor_fs (alr))
            sqme_soft = sqme_soft * symm_factor_fs
            sqme_coll = sqme_coll * symm_factor_fs
            sqme_cs = sqme_cs * symm_factor_fs
         end associate
         sqme_remn = compute_sqme_remnant_fsr (sqme_soft, sqme_cs, &
              rsub%real_kinematics%xi_max(i_phs), template%xi_cut, &
              rsub%real_kinematics%xi_tilde)
         sqme_tot = - sqme_soft - sqme_coll + sqme_cs + sqme_remn
      else
         sqme_tot = sqme_rad
      end if
      sqme_tot = sqme_tot * rsub%real_kinematics%jac_rand(i_phs)
      sqme_tot = sqme_tot * rsub%reg_data%regions(alr)%mult
    end associate
    if (debug_active (D_SUBTRACTION) .and. .not. &
         debug2_active (D_SUBTRACTION)) then
       call real_subtraction_register_debug_sqme (rsub, alr, emitter, &
            i_phs, sqme_rad, sqme_soft, sqme_coll=sqme_coll, sqme_cs=sqme_cs)
    else if (debug2_active (D_SUBTRACTION)) then
       call write_computation_status_fsr ()
    end if
  contains
    function compute_sqme_remnant_fsr (sqme_soft, sqme_cs, xi_max, xi_cut, xi_tilde) result (sqme_remn)
      real(default) :: sqme_remn
      real(default), intent(in) :: sqme_soft, sqme_cs, xi_max, xi_cut, xi_tilde
      if (debug_on) call msg_debug (D_SUBTRACTION, "compute_sqme_remnant_fsr")
      sqme_remn = (sqme_soft - sqme_cs) * log (xi_max) * xi_tilde / xi_cut
    end function compute_sqme_remnant_fsr

    subroutine write_computation_status_fsr (passed, total, region_type, full)
      integer, intent(in), optional :: passed, total
      character(*), intent(in), optional :: region_type
      integer :: i_born
      integer :: u
      real(default) :: xi
      logical :: yorn
      logical, intent(in), optional :: full
      yorn = .true.
      if (present (full)) yorn = full
      if (debug_on)  call msg_debug &
           (D_SUBTRACTION, "real_subtraction_evaluate_region_fsr")
      u = given_output_unit (); if (u < 0) return
      i_born = rsub%reg_data%regions(alr)%uborn_index
      xi = rsub%real_kinematics%xi_max (i_phs) * rsub%real_kinematics%xi_tilde
      write (u,'(A,I2)') 'rsub%purpose: ', rsub%purpose
      write (u,'(A,I4)') 'alr: ', alr
      write (u,'(A,I3)') 'emitter: ', emitter
      write (u,'(A,I3)') 'i_phs: ', i_phs
      write (u,'(A,F6.4)') 'xi_max: ', rsub%real_kinematics%xi_max (i_phs)
      write (u,'(A,F6.4)') 'xi_cut: ', rsub%real_kinematics%xi_max(i_phs) * &
           rsub%settings%fks_template%xi_cut
      write (u,'(A,F6.4,2X,A,F6.4)') 'xi: ', xi, 'y: ', &
           rsub%real_kinematics%y (i_phs)
      if (yorn) then
         write (u,'(A,ES16.9)')  'sqme_born: ', rsub%sqme_born(i_born)
         write (u,'(A,ES16.9)')  'sqme_real: ', sqme_rad
         write (u,'(A,ES16.9)')  'sqme_soft: ', sqme_soft
         write (u,'(A,ES16.9)')  'sqme_coll: ', sqme_coll
         write (u,'(A,ES16.9)')  'sqme_coll-soft: ', sqme_cs
         write (u,'(A,ES16.9)')  'sqme_remn: ', sqme_remn
         write (u,'(A,ES16.9)')  'sqme_tot: ', sqme_tot
         if (present (passed) .and. present (total) .and. &
              present (region_type)) &
              write (u,'(A)') char (str (passed) // " of " // str (total) // &
              " " // region_type // " points passed in total")
      end if
      write (u,'(A,ES16.9)')  'jacobian - real: ', &
           rsub%real_kinematics%jac(i_phs)%jac(1)
      write (u,'(A,ES16.9)')  'jacobian - soft: ', &
           rsub%real_kinematics%jac(i_phs)%jac(2)
      write (u,'(A,ES16.9)')  'jacobian - coll: ', &
           rsub%real_kinematics%jac(i_phs)%jac(3)
    end subroutine write_computation_status_fsr
  end function real_subtraction_evaluate_region_fsr

  subroutine real_subtraction_register_debug_sqme (rsub, alr, emitter, i_phs, sqme_rad, sqme_soft,&
      sqme_coll, sqme_cs, sqme_coll_plus, sqme_coll_minus, sqme_cs_plus, sqme_cs_minus)
    class(real_subtraction_t), intent(in) :: rsub
    integer, intent(in) :: alr, emitter, i_phs
    real(default), intent(in) :: sqme_rad, sqme_soft
    real(default), intent(in), optional :: sqme_coll, sqme_cs, sqme_coll_plus, sqme_coll_minus, sqme_cs_plus, sqme_cs_minus
    real(default), dimension(:), allocatable, save :: sqme_rad_store
    logical :: is_soft, is_collinear_plus, is_collinear_minus, is_fsr
    real(default), parameter :: soft_threshold = 0.001_default
    real(default), parameter :: coll_threshold = 0.99_default
    real(default), parameter :: rel_smallness = 0.01_default
    real(default) :: sqme_dummy, this_sqme_rad, y, xi_tilde
    logical, dimension(:), allocatable, save :: count_alr

    if (.not. allocated (sqme_rad_store)) then
       allocate (sqme_rad_store (rsub%reg_data%n_regions))
       sqme_rad_store = zero
    end if
    if (rsub%radiation_event) then
       sqme_rad_store(alr) = sqme_rad
    else
       if (.not. allocated (count_alr)) then
          allocate (count_alr (rsub%reg_data%n_regions))
          count_alr = .false.
       end if

       if (is_massless_vector (rsub%reg_data%regions(alr)%flst_real%flst(rsub%reg_data%n_legs_real))) then
          xi_tilde = rsub%real_kinematics%xi_tilde
          is_soft = xi_tilde < soft_threshold
       else
          is_soft = .false.
       end if
       y = rsub%real_kinematics%y(i_phs)
       is_collinear_plus = y > coll_threshold .and. &
            rsub%reg_data%regions(alr)%has_collinear_divergence()
       is_collinear_minus = -y > coll_threshold .and. &
            rsub%reg_data%regions(alr)%has_collinear_divergence()

       is_fsr = emitter > rsub%isr_kinematics%n_in

       if (is_fsr) then
          if (.not. present(sqme_coll) .or. .not. present(sqme_cs)) &
             call msg_error ("real_subtraction_register_debug_sqme: Wrong arguments for FSR")
       else
          if (.not. present(sqme_coll_plus) .or. .not. present(sqme_coll_minus) &
               .or. .not. present(sqme_cs_plus) .or. .not. present(sqme_cs_minus)) &
             call msg_error ("real_subtraction_register_debug_sqme: Wrong arguments for ISR")
       end if

       this_sqme_rad = sqme_rad_store(alr)
       if (is_soft .and. .not. is_collinear_plus .and. .not. is_collinear_minus) then
          if ( .not. nearly_equal (this_sqme_rad, sqme_soft, &
               abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
             call msg_print_color (char ("Soft MEs do not match in region " // str (alr)), COL_RED)
          else
             call msg_print_color (char ("sqme_soft OK in region " // str (alr)), COL_GREEN)
          end if
          print *, 'this_sqme_rad, sqme_soft =    ', this_sqme_rad, sqme_soft
       end if

       if (is_collinear_plus .and. .not. is_soft) then
          if (is_fsr) then
             if ( .not. nearly_equal (this_sqme_rad, sqme_coll, &
                 abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
               call msg_print_color (char ("Collinear MEs do not match in region " // str (alr)), COL_RED)
             else
               call msg_print_color (char ("sqme_coll OK in region " // str (alr)), COL_GREEN)
             end if
             print *, 'this_sqme_rad, sqme_coll =    ', this_sqme_rad, sqme_coll
          else
             if ( .not. nearly_equal (this_sqme_rad, sqme_coll_plus, &
                  abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
             call msg_print_color (char ("Collinear MEs do not match in region " // str (alr)), COL_RED)
             else
                call msg_print_color (char ("sqme_coll_plus OK in region " // str (alr)), COL_GREEN)
             end if
             print *, 'this_sqme_rad, sqme_coll_plus =    ', this_sqme_rad, sqme_coll_plus
          end if
       end if

       if (is_collinear_minus .and. .not. is_soft) then
          if (.not. is_fsr) then
             if ( .not. nearly_equal (this_sqme_rad, sqme_coll_minus, &
                  abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
                call msg_print_color (char ("Collinear MEs do not match in region " // str (alr)), COL_RED)
             else
                call msg_print_color (char ("sqme_coll_minus OK in region " // str (alr)), COL_GREEN)
             end if
             print *, 'this_sqme_rad, sqme_coll_minus =    ', this_sqme_rad, sqme_coll_minus
          end if
       end if

       if (is_soft .and. is_collinear_plus) then
          if (is_fsr) then
             if ( .not. nearly_equal (this_sqme_rad, sqme_cs, &
                  abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
                call msg_print_color (char ("Soft-collinear MEs do not match in region " // str (alr)), COL_RED)
             else
                call msg_print_color (char ("sqme_cs OK in region " // str (alr)), COL_GREEN)
             end if
             print *, 'this_sqme_rad, sqme_cs =    ', this_sqme_rad, sqme_cs
          else
             if ( .not. nearly_equal (this_sqme_rad, sqme_cs_plus, &
                  abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
                call msg_print_color (char ("Soft-collinear MEs do not match in region " // str (alr)), COL_RED)
             else
                call msg_print_color (char ("sqme_cs_plus OK in region " // str (alr)), COL_GREEN)
             end if
             print *, 'this_sqme_rad, sqme_cs_plus =    ', this_sqme_rad, sqme_cs_plus
          end if
       end if

       if (is_soft .and. is_collinear_minus) then
          if (.not. is_fsr) then
             if ( .not. nearly_equal (this_sqme_rad, sqme_cs_minus, &
                  abs_smallness=tiny(1._default), rel_smallness=rel_smallness)) then
                call msg_print_color (char ("Soft-collinear MEs do not match in region " // str (alr)), COL_RED)
             else
                call msg_print_color (char ("sqme_cs_minus OK in region " // str (alr)), COL_GREEN)
             end if
             print *, 'this_sqme_rad, sqme_cs_minus =    ', this_sqme_rad, sqme_cs_minus
          end if
       end if

       count_alr (alr) = .true.
       if (all (count_alr)) then
          deallocate (count_alr)
          deallocate (sqme_rad_store)
       end if
    end if
  end subroutine real_subtraction_register_debug_sqme

  module function real_subtraction_evaluate_region_isr (rsub, alr, &
       emitter, i_phs, i_res, alpha_coupling) result (sqme_tot)
    real(default) :: sqme_tot
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, emitter, i_phs, i_res
    real(default), intent(in) :: alpha_coupling
    real(default) :: sqme_rad, sqme_soft, sqme_coll_plus, sqme_coll_minus
    real(default) :: sqme_cs_plus, sqme_cs_minus
    real(default) :: sqme_remn
    sqme_rad = zero; sqme_soft = zero;
    sqme_coll_plus = zero; sqme_coll_minus = zero
    sqme_cs_plus = zero; sqme_cs_minus = zero
    sqme_remn = zero
    associate (region => rsub%reg_data%regions(alr), template => rsub%settings%fks_template)
      if (rsub%radiation_event) then
         sqme_rad = rsub%sqme_real_non_sub (rsub%reg_data%get_matrix_element_index (alr), i_phs)
         call evaluate_fks_factors (sqme_rad, rsub%reg_data, rsub%real_kinematics, &
              alr, i_phs, emitter, i_res)
         call apply_kinematic_factors_radiation (sqme_rad, rsub%real_kinematics, i_phs)
      end if
      if (rsub%subtraction_event .and. .not. (rsub%subtraction_deactivated &
            .or. region%nlo_correction_type == "none")) then
         call rsub%evaluate_subtraction_terms_isr (alr, emitter, i_phs, i_res, alpha_coupling, &
              sqme_soft, sqme_coll_plus, sqme_coll_minus, sqme_cs_plus, sqme_cs_minus)
         call apply_kinematic_factors_subtraction_isr (sqme_soft, sqme_coll_plus, &
              sqme_coll_minus, sqme_cs_plus, sqme_cs_minus, rsub%real_kinematics, i_phs)
         associate (symm_factor_fs => rsub%reg_data%born_to_real_symm_factor_fs (alr))
            sqme_soft = sqme_soft * symm_factor_fs
            sqme_coll_plus = sqme_coll_plus * symm_factor_fs
            sqme_coll_minus = sqme_coll_minus * symm_factor_fs
            sqme_cs_plus = sqme_cs_plus * symm_factor_fs
            sqme_cs_minus = sqme_cs_minus * symm_factor_fs
         end associate
         sqme_remn = compute_sqme_remnant_isr (rsub%isr_kinematics%isr_mode, &
              sqme_soft, sqme_cs_plus, sqme_cs_minus, &
              rsub%isr_kinematics, rsub%real_kinematics, i_phs, template%xi_cut)
         sqme_tot = - sqme_soft - sqme_coll_plus - sqme_coll_minus &
              + sqme_cs_plus + sqme_cs_minus + sqme_remn
      else
         sqme_tot = sqme_rad
      end if
    end associate
    sqme_tot = sqme_tot * rsub%real_kinematics%jac_rand (i_phs)
    sqme_tot = sqme_tot * rsub%reg_data%regions(alr)%mult
    if (debug_active (D_SUBTRACTION) .and. .not. debug2_active (D_SUBTRACTION)) then
       call real_subtraction_register_debug_sqme (rsub, alr, emitter, i_phs, sqme_rad,&
            sqme_soft, sqme_coll_plus=sqme_coll_plus, sqme_coll_minus=sqme_coll_minus,&
            sqme_cs_plus=sqme_cs_plus, sqme_cs_minus=sqme_cs_minus)
    else if (debug2_active (D_SUBTRACTION)) then
       call write_computation_status_isr ()
    end if
  contains
    function compute_sqme_remnant_isr (isr_mode, sqme_soft, sqme_cs_plus, sqme_cs_minus, &
       isr_kinematics, real_kinematics, i_phs, xi_cut) result (sqme_remn)
      real(default) :: sqme_remn
      integer, intent(in) :: isr_mode
      real(default), intent(in) :: sqme_soft, sqme_cs_plus, sqme_cs_minus
      type(isr_kinematics_t), intent(in) :: isr_kinematics
      type(real_kinematics_t), intent(in) :: real_kinematics
      integer, intent(in) :: i_phs
      real(default), intent(in) :: xi_cut
      real(default) :: xi_tilde, xi_max, xi_max_plus, xi_max_minus, xb_plus, xb_minus
      xi_max = real_kinematics%xi_max (i_phs)
      xi_tilde = real_kinematics%xi_tilde
      select case (isr_mode)
      case (SQRTS_VAR)
         xb_plus = isr_kinematics%x(I_PLUS)
         xb_minus = isr_kinematics%x(I_MINUS)
         xi_max_plus = one - xb_plus
         xi_max_minus = one - xb_minus
      case (SQRTS_FIXED)
         xi_max_plus = real_kinematics%xi_max (i_phs)
         xi_max_minus = real_kinematics%xi_max (i_phs)
      end select
      sqme_remn = log (xi_max) * xi_tilde * sqme_soft &
           - log (xi_max_plus) * xi_tilde * sqme_cs_plus &
           - log (xi_max_minus) * xi_tilde * sqme_cs_minus
      sqme_remn = sqme_remn / xi_cut
    end function compute_sqme_remnant_isr

    subroutine write_computation_status_isr (unit)
       integer, intent(in), optional :: unit
       integer :: i_born
       integer :: u
       real(default) :: xi
       u = given_output_unit (unit); if (u < 0) return
       i_born = rsub%reg_data%regions(alr)%uborn_index
       xi = rsub%real_kinematics%xi_max (i_phs) * rsub%real_kinematics%xi_tilde
       write (u,'(A,I4)') 'alr: ', alr
       write (u,'(A,I2)') 'emitter: ', emitter
       write (u,'(A,F4.2)') 'xi_max: ', rsub%real_kinematics%xi_max (i_phs)
       print *, 'xi: ', xi, 'y: ', rsub%real_kinematics%y (i_phs)
       print *, 'xb1: ', rsub%isr_kinematics%x(1), 'xb2: ', rsub%isr_kinematics%x(2)
       print *, 'random jacobian: ', rsub%real_kinematics%jac_rand (i_phs)
       write (u,'(A,ES16.9)')  'sqme_born: ', rsub%sqme_born(i_born)
       write (u,'(A,ES16.9)')  'sqme_real: ', sqme_rad
       write (u,'(A,ES16.9)')  'sqme_soft: ', sqme_soft
       write (u,'(A,ES16.9)')  'sqme_coll_plus: ', sqme_coll_plus
       write (u,'(A,ES16.9)')  'sqme_coll_minus: ', sqme_coll_minus
       write (u,'(A,ES16.9)')  'sqme_cs_plus: ', sqme_cs_plus
       write (u,'(A,ES16.9)')  'sqme_cs_minus: ', sqme_cs_minus
       write (u,'(A,ES16.9)')  'sqme_remn: ', sqme_remn
       write (u,'(A,ES16.9)')  'sqme_tot: ', sqme_tot
       write (u,'(A,ES16.9)')  'jacobian - real: ', rsub%real_kinematics%jac(i_phs)%jac(1)
       write (u,'(A,ES16.9)')  'jacobian - soft: ', rsub%real_kinematics%jac(i_phs)%jac(2)
       write (u,'(A,ES16.9)')  'jacobian - collplus: ', rsub%real_kinematics%jac(i_phs)%jac(3)
       write (u,'(A,ES16.9)')  'jacobian - collminus: ', rsub%real_kinematics%jac(i_phs)%jac(4)
    end subroutine write_computation_status_isr
  end function real_subtraction_evaluate_region_isr

  module subroutine real_subtraction_evaluate_subtraction_terms_fsr &
       (rsub, alr, emitter, i_phs, i_res, alpha_coupling, sqme_soft, &
        sqme_coll, sqme_cs)
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, emitter, i_phs, i_res
    real(default), intent(in) :: alpha_coupling
    real(default), intent(out) :: sqme_soft, sqme_coll, sqme_cs
    if (debug_on)  call msg_debug &
         (D_SUBTRACTION, "real_subtraction_evaluate_subtraction_terms_fsr")
    sqme_soft = zero; sqme_coll = zero; sqme_cs = zero
    associate (xi_tilde => rsub%real_kinematics%xi_tilde, &
         y => rsub%real_kinematics%y(i_phs), &
         template => rsub%settings%fks_template)
      if (template%xi_cut > xi_tilde) &
           sqme_soft = rsub%compute_sub_soft &
                         (alr, emitter, i_phs, i_res, alpha_coupling)
      if (y - 1 + template%delta_o > 0) &
           sqme_coll = rsub%compute_sub_coll &
                         (alr, emitter, i_phs, alpha_coupling)
      if (template%xi_cut > xi_tilde .and. y - 1 + template%delta_o > 0) &
           sqme_cs = rsub%compute_sub_coll_soft &
                         (alr, emitter, i_phs, alpha_coupling)
      if (debug2_active (D_SUBTRACTION)) then
         print *, "FSR Cutoff:"
         print *, "sub_soft: ", &
              template%xi_cut > xi_tilde, "(ME: ", sqme_soft, ")"
         print *, "sub_coll: ", &
              (y - 1 + template%delta_o) > 0, "(ME: ", sqme_coll, ")"
         print *, "sub_coll_soft: ", &
              template%xi_cut > xi_tilde .and. (y - 1 + template%delta_o) > 0, &
              "(ME: ", sqme_cs, ")"
      end if
    end associate
  end subroutine real_subtraction_evaluate_subtraction_terms_fsr

  subroutine evaluate_fks_factors (sqme, reg_data, real_kinematics, &
      alr, i_phs, emitter, i_res)
    real(default), intent(inout) :: sqme
    type(region_data_t), intent(inout) :: reg_data
    type(real_kinematics_t), intent(in), target :: real_kinematics
    integer, intent(in) :: alr, i_phs, emitter, i_res
    real(default) :: s_alpha
    type(phs_point_set_t), pointer :: p_real => null ()
    if (reg_data%has_pseudo_isr ()) then
       p_real => real_kinematics%p_real_onshell (i_phs)
    else
       p_real => real_kinematics%p_real_cms
    end if
    s_alpha = reg_data%get_svalue (p_real%get_momenta(i_phs), alr, emitter, i_res)
    if (debug2_active (D_SUBTRACTION)) call msg_print_color('s_alpha', s_alpha, COL_YELLOW)
    if (s_alpha > one + tiny_07) call msg_fatal ("s_alpha > 1!")
    sqme = sqme * s_alpha
    associate (region => reg_data%regions(alr))
       if (emitter > reg_data%n_in) then
          if (debug2_active (D_SUBTRACTION)) &
               print *, 'Double FSR: ', region%double_fsr_factor (p_real%get_momenta(i_phs))
          sqme = sqme * region%double_fsr_factor (p_real%get_momenta(i_phs))
       end if
    end associate
  end subroutine evaluate_fks_factors

  subroutine apply_kinematic_factors_radiation (sqme, &
       real_kinematics, i_phs)
    real(default), intent(inout) :: sqme
    type(real_kinematics_t), intent(in) :: real_kinematics
    integer, intent(in) :: i_phs
    real(default) :: xi, xi_tilde, s_b
    xi_tilde = real_kinematics%xi_tilde
    xi = xi_tilde * real_kinematics%xi_max (i_phs)
    sqme = sqme * xi**2 / xi_tilde * real_kinematics%jac(i_phs)%jac(1)
  end subroutine apply_kinematic_factors_radiation

  subroutine apply_kinematic_factors_subtraction_fsr &
     (sqme_soft, sqme_coll, sqme_cs, real_kinematics, i_phs)
    real(default), intent(inout) :: sqme_soft, sqme_coll, sqme_cs
    type(real_kinematics_t), intent(in) :: real_kinematics
    integer, intent(in) :: i_phs
    real(default) :: xi_tilde, onemy
    xi_tilde = real_kinematics%xi_tilde
    onemy = one - real_kinematics%y(i_phs)
    sqme_soft = sqme_soft / onemy / xi_tilde
    sqme_coll = sqme_coll / onemy / xi_tilde
    sqme_cs = sqme_cs / onemy / xi_tilde
    associate (jac => real_kinematics%jac(i_phs)%jac)
       sqme_soft = sqme_soft * jac(2)
       sqme_coll = sqme_coll * jac(3)
       sqme_cs = sqme_cs * jac(2)
    end associate
  end subroutine apply_kinematic_factors_subtraction_fsr

  subroutine apply_kinematic_factors_subtraction_isr &
     (sqme_soft, sqme_coll_plus, sqme_coll_minus, sqme_cs_plus, &
      sqme_cs_minus, real_kinematics, i_phs)
    real(default), intent(inout) :: sqme_soft, sqme_coll_plus, sqme_coll_minus
    real(default), intent(inout) :: sqme_cs_plus, sqme_cs_minus
    type(real_kinematics_t), intent(in) :: real_kinematics
    integer, intent(in) :: i_phs
    real(default) :: xi_tilde, y, onemy, onepy
    xi_tilde = real_kinematics%xi_tilde
    y = real_kinematics%y (i_phs)
    onemy = one - y; onepy = one + y
    associate (jac => real_kinematics%jac(i_phs)%jac)
       sqme_soft = sqme_soft / (one - y**2) / xi_tilde * jac(2)
       sqme_coll_plus = sqme_coll_plus / onemy / xi_tilde / two * jac(3)
       sqme_coll_minus = sqme_coll_minus / onepy / xi_tilde / two * jac(4)
       sqme_cs_plus = sqme_cs_plus / onemy / xi_tilde / two * jac(2)
       sqme_cs_minus = sqme_cs_minus / onepy / xi_tilde / two * jac(2)
    end associate
  end subroutine apply_kinematic_factors_subtraction_isr

  module subroutine real_subtraction_evaluate_subtraction_terms_isr (rsub, &
      alr, emitter, i_phs, i_res, alpha_coupling, sqme_soft, sqme_coll_plus, &
      sqme_coll_minus, sqme_cs_plus, sqme_cs_minus)
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, emitter, i_phs, i_res
    real(default), intent(in) :: alpha_coupling
    real(default), intent(out) :: sqme_soft
    real(default), intent(out) :: sqme_coll_plus, sqme_coll_minus
    real(default), intent(out) :: sqme_cs_plus, sqme_cs_minus
    sqme_coll_plus = zero; sqme_cs_plus = zero
    sqme_coll_minus = zero; sqme_cs_minus = zero
    associate (xi_tilde => rsub%real_kinematics%xi_tilde, &
         y => rsub%real_kinematics%y(i_phs), &
         template => rsub%settings%fks_template)
      if (template%xi_cut > xi_tilde) &
           sqme_soft = rsub%compute_sub_soft &
           (alr, emitter, i_phs, i_res, alpha_coupling)
      if (emitter /= 2) then
         if (y - 1 + template%delta_i > 0) then
            sqme_coll_plus = &
                 rsub%compute_sub_coll (alr, 1, i_phs, alpha_coupling)
            if (template%xi_cut > xi_tilde) then
               sqme_cs_plus = &
                    rsub%compute_sub_coll_soft (alr, 1, i_phs, alpha_coupling)
            end if
         end if
      end if
      if (emitter /= 1) then
         if (-y - 1 + template%delta_i > 0) then
            sqme_coll_minus = &
                 rsub%compute_sub_coll (alr, 2, i_phs, alpha_coupling)
            if (template%xi_cut > xi_tilde) then
               sqme_cs_minus = &
                    rsub%compute_sub_coll_soft (alr, 2, i_phs, alpha_coupling)
            end if
         end if
      end if
      if (debug2_active (D_SUBTRACTION)) then
         print *, "ISR Cutoff:"
         print *, "y: ", y
         print *, "delta_i: ", template%delta_i
         print *, "emitter: ", emitter
         print *, "sub_soft: ", &
              template%xi_cut > xi_tilde, "(ME: ", sqme_soft, ")"
         print *, "sub_coll_plus: ", &
              (y - 1 + template%delta_i) > 0, "(ME: ", sqme_coll_plus, ")"
         print *, "sub_coll_minus: ", &
              (-y - 1 + template%delta_i) > 0, "(ME: ", sqme_coll_minus, ")"
         print *, "sub_coll_soft_plus: ", template%xi_cut > xi_tilde .and. &
              (y - 1 + template%delta_i) > 0, "(ME: ", sqme_cs_plus, ")"
         print *, "sub_coll_soft_minus: ", template%xi_cut > xi_tilde .and. &
              (-y - 1 + template%delta_i) > 0, "(ME: ", sqme_cs_minus, ")"
      end if
    end associate
  end subroutine real_subtraction_evaluate_subtraction_terms_isr

  module function real_subtraction_get_phs_factor (rsub, i_con) result (factor)
    real(default) :: factor
    class(real_subtraction_t), intent(in) :: rsub
    integer, intent(in) :: i_con
    real(default) :: s
    s = rsub%real_kinematics%xi_ref_momenta (i_con)**2
    factor = s / (8 * twopi3)
  end function real_subtraction_get_phs_factor

  module function real_subtraction_get_i_contributor (rsub, alr) result (i_con)
    integer :: i_con
    class(real_subtraction_t), intent(in) :: rsub
    integer, intent(in) :: alr
    if (allocated (rsub%reg_data%alr_to_i_contributor)) then
       i_con = rsub%reg_data%alr_to_i_contributor (alr)
    else
       i_con = 1
    end if
  end function real_subtraction_get_i_contributor

  module function real_subtraction_compute_sub_soft (rsub, alr, emitter, &
       i_phs, i_res, alpha_coupling) result (sqme_soft)
    real(default) :: sqme_soft
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, emitter, i_phs, i_res
    real(default), intent(in) :: alpha_coupling
    integer :: i_xi_ref, i_born
    real(default) :: q2, sf_factor
    type(vector4_t), dimension(:), allocatable :: p_born
    associate (real_kinematics => rsub%real_kinematics, &
            nlo_corr_type => rsub%reg_data%regions(alr)%nlo_correction_type, &
            sregion => rsub%reg_data%regions(alr))
      sqme_soft = zero
      if (sregion%has_soft_divergence ()) then
         i_xi_ref = rsub%sub_soft%i_xi_ref (alr, i_phs)
         q2 = real_kinematics%xi_ref_momenta (i_xi_ref)**2
         allocate (p_born (rsub%reg_data%n_legs_born))
         if (rsub%reg_data%has_pseudo_isr ()) then
            p_born = real_kinematics%p_born_onshell%get_momenta (1)
         else
            p_born = real_kinematics%p_born_cms%get_momenta (1)
         end if
         if (emitter > rsub%sub_soft%reg_data%n_in) then
            call rsub%sub_soft%create_softvec_fsr &
                 (p_born, real_kinematics%y_soft(i_phs), &
                 real_kinematics%phi, emitter, &
                 real_kinematics%xi_ref_momenta(i_xi_ref))
            sf_factor = one
         else
            call rsub%sub_soft%create_softvec_isr &
                 (real_kinematics%y_soft(i_phs), real_kinematics%phi)
            sf_factor = rsub%sf_factors(alr, 0)
         end if
         i_born = sregion%uborn_index
         select case (char (nlo_corr_type))
         case ("QCD")
            sqme_soft = rsub%sub_soft%compute &
            (p_born, rsub%sqme_born_color_c(:,:,i_born) * &
            sf_factor, real_kinematics%y(i_phs), &
            q2, alpha_coupling, alr, emitter, i_res)
         case ("EW")
            sqme_soft = rsub%sub_soft%compute &
            (p_born, rsub%sqme_born_charge_c(:,:,i_born) * &
            sf_factor, real_kinematics%y(i_phs), &
            q2, alpha_coupling, alr, emitter, i_res)
         end select
      end if
    end associate
    if (debug2_active (D_SUBTRACTION)) call check_soft_vector ()
  contains
    subroutine check_soft_vector ()
      !!! p_soft = p_gluon / E_gluon only in the soft limit
      !!! This check only has to be passed for ISR or for FSR if ?test_soft_limit = true is set.
      type(vector4_t) :: p_gluon
      if (debug_on) call msg_debug2 (D_SUBTRACTION, "Compare soft vector: ")
      print *, 'p_soft: ', rsub%sub_soft%p_soft%p
      print *, 'Normalized gluon momentum: '
      if (rsub%reg_data%has_pseudo_isr ()) then
         p_gluon = rsub%real_kinematics%p_real_onshell(thr_leg(emitter))%get_momentum &
              (i_phs, rsub%reg_data%n_legs_real)
      else
         p_gluon = rsub%real_kinematics%p_real_cms%get_momentum &
              (i_phs, rsub%reg_data%n_legs_real)
      end if
      call vector4_write (p_gluon / p_gluon%p(0), show_mass = .true.)
    end subroutine check_soft_vector
  end function real_subtraction_compute_sub_soft

  module function real_subtraction_get_spin_correlation_term &
       (rsub, alr, i_born, emitter) result (mom_times_sqme)
    real(default) :: mom_times_sqme
    class(real_subtraction_t), intent(in) :: rsub
    integer, intent(in) :: alr, i_born, emitter
    real(default), dimension(0:3) :: k_perp
    integer :: mu, nu
    if (rsub%sc_required(alr)) then
       if (debug2_active(D_SUBTRACTION)) call check_me_consistency ()
       associate (real_kin => rsub%real_kinematics)
         if (emitter > rsub%reg_data%n_in) then
            k_perp = real_subtraction_compute_k_perp_fsr ( &
                 real_kin%p_born_lab%get_momentum(1, emitter), &
                 rsub%real_kinematics%phi)
         else
            k_perp = real_subtraction_compute_k_perp_isr ( &
                 real_kin%p_born_lab%get_momentum(1, emitter), &
                 rsub%real_kinematics%phi)
         end if
       end associate
       mom_times_sqme = zero
       do mu = 1, 3
          do nu = 1, 3
             mom_times_sqme = mom_times_sqme + &
                  k_perp(mu) * k_perp(nu) * rsub%sqme_born_spin_c (mu, nu, emitter, i_born)
          end do
       end do
    else
       mom_times_sqme = zero
    end if
  contains
    subroutine check_me_consistency ()
      real(default) ::  sqme_sum
      if (debug_on) call msg_debug2 (D_SUBTRACTION, "Spin-correlation: Consistency check")
      sqme_sum = &
               - rsub%sqme_born_spin_c(1,1,emitter,i_born) &
               - rsub%sqme_born_spin_c(2,2,emitter,i_born) &
               - rsub%sqme_born_spin_c(3,3,emitter,i_born)
      if (.not. nearly_equal (sqme_sum, -rsub%sqme_born(i_born))) then
         print *, 'Spin-correlated matrix elements are not consistent: '
         print *, 'emitter: ', emitter
         print *, 'g^{mu,nu} B_{mu,nu}: ', -sqme_sum
         print *, 'all Born matrix elements: ', rsub%sqme_born
         call msg_fatal ("FAIL")
      else
         call msg_print_color ("Success", COL_GREEN)
      end if
    end subroutine check_me_consistency
  end function real_subtraction_get_spin_correlation_term

  module function real_subtraction_compute_k_perp_fsr &
       (p, phi) result (k_perp_fsr)
    real(default), dimension(0:3) :: k_perp_fsr
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: phi
    type(vector4_t) :: k
    type(vector3_t) :: vec
    type(lorentz_transformation_t) :: rot
    vec = p%p(1:3) / p%p(0)
    k%p(0) = zero
    k%p(1) = p%p(1); k%p(2) = p%p(2)
    k%p(3) = - (p%p(1)**2 + p%p(2)**2) / p%p(3)
    rot = rotation (cos(phi), sin(phi), vec)
    k = rot * k
    k%p(1:3) = k%p(1:3) / space_part_norm (k)
    k_perp_fsr = k%p
  end function real_subtraction_compute_k_perp_fsr

  module function real_subtraction_compute_k_perp_isr &
       (p, phi) result (k_perp_isr)
    real(default), dimension(0:3) :: k_perp_isr
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: phi
    k_perp_isr(0) = zero
    k_perp_isr(1) = sin(phi)
    k_perp_isr(2) = cos(phi)
    k_perp_isr(3) = zero
  end function real_subtraction_compute_k_perp_isr

  module function real_subtraction_compute_sub_coll &
       (rsub, alr, em, i_phs, alpha_coupling) result (sqme_coll)
    real(default) :: sqme_coll
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, em, i_phs
    real(default), intent(in) :: alpha_coupling
    real(default) :: xi, xi_max
    real(default) :: mom_times_sqme_spin_c, sqme_born_coll
    real(default) :: N_col
    integer :: i_con, i
    real(default) :: pfr
    associate (sregion => rsub%reg_data%regions(alr))
      sqme_coll = zero
      sqme_born_coll = zero
      N_col = 1
      if (sregion%has_collinear_divergence ()) then
         xi = rsub%real_kinematics%xi_tilde * rsub%real_kinematics%xi_max(i_phs)
         if (rsub%sub_coll%use_resonance_mappings) then
            i_con = rsub%reg_data%alr_to_i_contributor (alr)
         else
            i_con = 1
         end if
         mom_times_sqme_spin_c = &
              rsub%get_spin_correlation_term (alr, sregion%uborn_index, em)
         if (rsub%reg_data%nlo_correction_type == "EW" .and. &
             sregion%nlo_correction_type == "QCD" .and. &
             qcd_ew_interferences (sregion%flst_uborn%flst)) then
            do i = 1, size (sregion%flst_uborn%flst)
               if (is_quark (sregion%flst_uborn%flst (i))) then
                  sqme_born_coll = &
                       -rsub%sqme_born_color_c (i, i, sregion%uborn_index)/CF
                  exit
               end if
            end do
         else
            sqme_born_coll = rsub%sqme_born(sregion%uborn_index)
         end if
         if (em <= rsub%sub_coll%n_in) then
            select case (rsub%isr_kinematics%isr_mode)
            case (SQRTS_FIXED)
               xi_max = rsub%real_kinematics%xi_max(i_phs)
            case (SQRTS_VAR)
               xi_max = one - rsub%isr_kinematics%x(em)
            end select
            xi = rsub%real_kinematics%xi_tilde * xi_max
            if (sregion%nlo_correction_type == "QCD") then
               call rsub%sub_coll%set_parameters (CA = CA, CF = CF, TR = TR)
            else if (sregion%nlo_correction_type == "EW") then
               if (is_quark (sregion%flst_real%flst(size(sregion%flst_real%flst)))) N_col = 3
               call rsub%sub_coll%set_parameters (CA = zero, &
                    CF = sregion%flst_real%charge(em)**2, &
                    TR = N_col*sregion%flst_real%charge(size(sregion%flst_real%flst))**2)
            end if
            sqme_coll = rsub%sub_coll%compute_isr (em, sregion%flst_real%flst, &
                 rsub%real_kinematics%p_born_lab%phs_point(1)%get (), &
                 sqme_born_coll * rsub%sf_factors(alr, em), &
                 mom_times_sqme_spin_c * rsub%sf_factors(alr, em), &
                 xi, alpha_coupling, rsub%isr_kinematics%isr_mode)
         else
            if (sregion%nlo_correction_type == "QCD") then
               call rsub%sub_coll%set_parameters (CA = CA, CF = CF, TR = TR)
            else if (sregion%nlo_correction_type == "EW") then
               if (is_quark (sregion%flst_real%flst(sregion%emitter))) N_col = 3
               call rsub%sub_coll%set_parameters (CA = zero, &
                    CF = sregion%flst_real%charge(sregion%emitter)**2, &
                    TR = N_col*sregion%flst_real%charge(sregion%emitter)**2)
            end if
            sqme_coll = rsub%sub_coll%compute_fsr (sregion%emitter, &
                 sregion%flst_real%flst, &
                 rsub%real_kinematics%xi_ref_momenta (i_con), &
                 rsub%real_kinematics%p_born_lab%get_momenta(1), &
                 sqme_born_coll, &
                 mom_times_sqme_spin_c, &
                 xi, alpha_coupling, sregion%double_fsr)
            if (rsub%sub_coll%use_resonance_mappings) then
               select type (fks_mapping => rsub%reg_data%fks_mapping)
               type is (fks_mapping_resonances_t)
                  pfr = fks_mapping%get_resonance_weight (alr, &
                       rsub%real_kinematics%p_born_cms%get_momenta(1))
               end select
               sqme_coll = sqme_coll * pfr
            end if
         end if
      end if
    end associate
  end function real_subtraction_compute_sub_coll

  module function real_subtraction_compute_sub_coll_soft &
       (rsub, alr, em, i_phs, alpha_coupling) result (sqme_cs)
    real(default) :: sqme_cs
    class(real_subtraction_t), intent(inout) :: rsub
    integer, intent(in) :: alr, em, i_phs
    real(default), intent(in) :: alpha_coupling
    real(default) :: mom_times_sqme_spin_c, sqme_born_coll
    real(default) :: N_col
    integer :: i_con, em_pdf, i
    associate (sregion => rsub%reg_data%regions(alr))
      sqme_cs = zero
      sqme_born_coll = zero
      N_col = 1
      if (sregion%has_collinear_divergence ()) then
         if (rsub%sub_coll%use_resonance_mappings) then
            i_con = rsub%reg_data%alr_to_i_contributor (alr)
         else
            i_con = 1
         end if
         mom_times_sqme_spin_c = rsub%get_spin_correlation_term (alr, sregion%uborn_index, em)
         if (rsub%reg_data%nlo_correction_type == "EW" .and. &
             sregion%nlo_correction_type == "QCD" .and. &
             qcd_ew_interferences (sregion%flst_uborn%flst)) then
            do i = 1, size (sregion%flst_uborn%flst)
               if (is_quark (sregion%flst_uborn%flst (i))) then
                  sqme_born_coll = -rsub%sqme_born_color_c (i, i, sregion%uborn_index)/CF
                  exit
               end if
            end do
         else
            sqme_born_coll = rsub%sqme_born(sregion%uborn_index)
         end if
         if (em <= rsub%sub_coll%n_in) then
            em_pdf = 0
            if (sregion%nlo_correction_type == "QCD") then
               call rsub%sub_coll%set_parameters (CA = CA, CF = CF, TR = TR)
            else if (sregion%nlo_correction_type == "EW") then
               if (is_quark (sregion%flst_real%flst(size(sregion%flst_real%flst)))) N_col = 3
               call rsub%sub_coll%set_parameters (CA = zero, &
                    CF = sregion%flst_real%charge(em)**2, &
                    TR = N_col*sregion%flst_real%charge(size(sregion%flst_real%flst))**2)
            end if
            sqme_cs = rsub%sub_coll%compute_isr (em, sregion%flst_real%flst, &
                 rsub%real_kinematics%p_born_lab%phs_point(1)%get (), &
                 sqme_born_coll * rsub%sf_factors(alr, em_pdf), &
                 mom_times_sqme_spin_c * rsub%sf_factors(alr, em_pdf), &
                 zero, alpha_coupling, rsub%isr_kinematics%isr_mode)
         else
            if (sregion%nlo_correction_type == "QCD") then
               call rsub%sub_coll%set_parameters (CA = CA, CF = CF, TR = TR)
            else if (sregion%nlo_correction_type == "EW") then
               if (is_quark (sregion%flst_real%flst(sregion%emitter))) N_col = 3
               call rsub%sub_coll%set_parameters (CA = zero, &
                    CF = sregion%flst_real%charge(sregion%emitter)**2, &
                    TR = N_col*sregion%flst_real%charge(sregion%emitter)**2)
            end if
            sqme_cs = rsub%sub_coll%compute_fsr (sregion%emitter, sregion%flst_real%flst, &
                 rsub%real_kinematics%xi_ref_momenta(i_con), &
                 rsub%real_kinematics%p_born_lab%phs_point(1)%get (), &
                 sqme_born_coll, &
                 mom_times_sqme_spin_c, &
                 zero, alpha_coupling, sregion%double_fsr)
         end if
      end if
    end associate
  end function real_subtraction_compute_sub_coll_soft

  module function real_subtraction_requires_spin_correlations &
       (rsub) result (val)
    logical :: val
    class(real_subtraction_t), intent(in) :: rsub
    val = any (rsub%sc_required)
  end function real_subtraction_requires_spin_correlations

  module subroutine real_subtraction_final (rsub)
    class(real_subtraction_t), intent(inout) :: rsub
    call rsub%sub_soft%final ()
    call rsub%sub_coll%final ()
    !!! Finalization of region data is done in pcm_nlo_final
    if (associated (rsub%reg_data)) nullify (rsub%reg_data)
    !!! Finalization of real kinematics is done in pcm_instance_nlo_final
    if (associated (rsub%real_kinematics)) nullify (rsub%real_kinematics)
    if (associated (rsub%isr_kinematics)) nullify (rsub%isr_kinematics)
    if (allocated (rsub%sqme_real_non_sub)) deallocate (rsub%sqme_real_non_sub)
    if (allocated (rsub%sqme_born)) deallocate (rsub%sqme_born)
    if (allocated (rsub%sf_factors)) deallocate (rsub%sf_factors)
    if (allocated (rsub%sqme_real_arr)) deallocate (rsub%sqme_real_arr)
    if (allocated (rsub%sqme_born_color_c)) deallocate (rsub%sqme_born_color_c)
    if (allocated (rsub%sqme_born_charge_c)) deallocate (rsub%sqme_born_charge_c)
    if (allocated (rsub%sc_required)) deallocate (rsub%sc_required)
    if (allocated (rsub%selected_alr)) deallocate (rsub%selected_alr)
  end subroutine real_subtraction_final

  module function powheg_damping_simple_get_f (partition, p) result (f)
    real(default) :: f
    class(powheg_damping_simple_t), intent(in) :: partition
    type(vector4_t), intent(in), dimension(:) :: p
    !!! real(default) :: pt2
    f = 1
    call msg_bug ("Simple damping currently not available")
    !!! TODO (cw-2017-03-01) Compute pt2 from emitter)
    !!! f = partition%h2 / (pt2 + partition%h2)
  end function powheg_damping_simple_get_f

  module subroutine powheg_damping_simple_init (partition, scale, reg_data)
    class(powheg_damping_simple_t), intent(out) :: partition
    real(default), intent(in) :: scale
    type(region_data_t), intent(in) :: reg_data
    partition%h2 = scale**2
  end subroutine powheg_damping_simple_init

  module subroutine powheg_damping_simple_write (partition, unit)
    class(powheg_damping_simple_t), intent(in) :: partition
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit); if (u < 0) return
    write (u, "(1x,A)") "Powheg damping simple: "
    write (u, "(1x,A, "// FMT_15 // ")") "scale h2: ", partition%h2
  end subroutine powheg_damping_simple_write

  module subroutine real_partition_fixed_order_init (partition, scale, reg_data)
    class(real_partition_fixed_order_t), intent(out) :: partition
    real(default), intent(in) :: scale
    type(region_data_t), intent(in) :: reg_data
  end subroutine real_partition_fixed_order_init

  module subroutine real_partition_fixed_order_write (partition, unit)
    class(real_partition_fixed_order_t), intent(in) :: partition
    integer, intent(in), optional :: unit
  end subroutine real_partition_fixed_order_write

  module function real_partition_fixed_order_get_f (partition, p) result (f)
    real(default) :: f
    class(real_partition_fixed_order_t), intent(in) :: partition
    type(vector4_t), intent(in), dimension(:) :: p
    integer :: i, em
    f = zero
    PAIRS: do i = 1, size (partition%fks_pairs)
       associate (ii => partition%fks_pairs(i)%ireg)
          if (ii(1) == 0) then
            IS: do em = 1, 2
               if ((p(em) + p(ii(2)))**1 < p(em)**1 + p(ii(2))**1 + &
                    partition%scale) then
                  f = one
                  exit PAIRS
               end if
            end do IS
          else
             if ((p(ii(1)) + p(ii(2)))**1 < p(ii(1))**1 + p(ii(2))**1 + &
                  partition%scale) then
                f = one
                exit PAIRS
             end if
          end if
       end associate
    end do PAIRS
  end function real_partition_fixed_order_get_f


end submodule real_subtraction_s

