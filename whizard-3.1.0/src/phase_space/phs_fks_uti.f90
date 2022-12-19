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

module phs_fks_uti

  use kinds, only: default
  use format_utils, only: write_separator, pac_fmt
  use format_defs, only: FMT_15, FMT_19
  use numeric_utils, only: nearly_equal
  use constants, only: tiny_07, zero, one, two
  use lorentz
  use phs_points, only: assignment(=)

  use physics_defs, only: THR_POS_B, THR_POS_BBAR, THR_POS_WP, THR_POS_WM, THR_POS_GLUON
  use physics_defs, only: thr_leg

  use resonances, only: resonance_contributors_t
  use phs_fks

  implicit none
  private

  public :: phs_fks_generator_1
  public :: phs_fks_generator_2
  public :: phs_fks_generator_3
  public :: phs_fks_generator_4
  public :: phs_fks_generator_5
  public :: phs_fks_generator_6
  public :: phs_fks_generator_7

contains

  subroutine phs_fks_generator_1 (u)
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: p_real
    integer :: emitter, i_phs
    real(default) :: x1, x2, x3
    real(default), parameter :: sqrts = 250.0_default
    type(phs_identifier_t), dimension(2) :: phs_identifiers
    write (u, "(A)") "* Test output: phs_fks_generator_1"
    write (u, "(A)") "* Purpose: Create massless fsr phase space"
    write (u, "(A)")

    allocate (p_born (4))
    p_born(1)%p(0) = 125.0_default
    p_born(1)%p(1:2) = 0.0_default
    p_born(1)%p(3) = 125.0_default
    p_born(2)%p(0) = 125.0_default
    p_born(2)%p(1:2) = 0.0_default
    p_born(2)%p(3) = -125.0_default
    p_born(3)%p(0) = 125.0_default
    p_born(3)%p(1) = -39.5618_default
    p_born(3)%p(2) = -20.0791_default
    p_born(3)%p(3) = -114.6957_default
    p_born(4)%p(0) = 125.0_default
    p_born(4)%p(1:3) = -p_born(3)%p(1:3)

    allocate (generator%isr_kinematics)
    generator%n_in = 2
    generator%isr_kinematics%isr_mode = SQRTS_FIXED
    call generator%set_xi_and_y_bounds ()

    call generator%set_sqrts_hat (sqrts)

    write (u, "(A)") "* Use four-particle phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    write (u, "(A)") "***********************"
    write (u, "(A)")

    x1 = 0.5_default; x2 = 0.25_default; x3 = 0.75_default
    write (u, "(A)" ) "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,1X,A,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (4, 2, 2, 1)

    allocate (generator%emitters (2))
    generator%emitters(1) = 3; generator%emitters(2) = 4
    allocate (generator%m2 (4))
    generator%m2 = zero
    allocate (generator%is_massive (4))
    generator%is_massive(1:2) = .false.
    generator%is_massive(3:4) = .true.
    phs_identifiers(1)%emitter = 3
    phs_identifiers(2)%emitter = 4
    call generator%compute_xi_ref_momenta (p_born)
    call generator%generate_radiation_variables ([x1,x2,x3], p_born, phs_identifiers)
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       call generator%compute_xi_max (emitter, i_phs, p_born, &
            generator%real_kinematics%xi_max(i_phs))
    end do
    write (u, "(A)")  &
         "* With these, the following radiation variables have been produced:"
    associate (rad_var => generator%real_kinematics)
      write (u, "(A,F3.2)") "xi_tilde: ", rad_var%xi_tilde
      write (u, "(A,F3.2)") "y: " , rad_var%y(1)
      write (u, "(A,F3.2)") "phi: ", rad_var%phi
    end associate
    call write_separator (u)
    write (u, "(A)") "Produce real momenta: "
    i_phs = 1; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter

    allocate (p_real (5))
    call generator%generate_fsr (emitter, i_phs, p_born, p_real)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    call write_separator (u)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_1"

  end subroutine phs_fks_generator_1

  subroutine phs_fks_generator_2 (u)
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: p_real
    integer :: emitter, i_phs
    real(default) :: x1, x2, x3
    real(default), parameter :: sqrts_hadronic = 250.0_default
    type(phs_identifier_t), dimension(2) :: phs_identifiers
    write (u, "(A)") "* Test output: phs_fks_generator_2"
    write (u, "(A)") "* Purpose: Create massless ISR phase space"
    write (u, "(A)")


    allocate (p_born (4))
    p_born(1)%p(0) = 114.661_default
    p_born(1)%p(1:2) = 0.0_default
    p_born(1)%p(3) = 114.661_default
    p_born(2)%p(0) = 121.784_default
    p_born(2)%p(1:2) = 0.0_default
    p_born(2)%p(3) = -121.784_default
    p_born(3)%p(0) = 115.148_default
    p_born(3)%p(1) = -46.250_default
    p_born(3)%p(2) = -37.711_default
    p_born(3)%p(3) = 98.478_default
    p_born(4)%p(0) = 121.296_default
    p_born(4)%p(1:2) = -p_born(3)%p(1:2)
    p_born(4)%p(3) = -105.601_default

    phs_identifiers(1)%emitter = 1
    phs_identifiers(2)%emitter = 2

    allocate (generator%emitters (2))
    allocate (generator%isr_kinematics)
    generator%emitters(1) = 1; generator%emitters(2) = 2
    generator%sqrts = sqrts_hadronic
    allocate (generator%isr_kinematics%beam_energy(2))
    generator%isr_kinematics%beam_energy = sqrts_hadronic / two
    call generator%set_sqrts_hat (sqrts_hadronic)
    call generator%set_isr_kinematics (p_born)
    generator%n_in = 2
    generator%isr_kinematics%isr_mode = SQRTS_VAR
    call generator%set_xi_and_y_bounds ()
    write (u, "(A)") "* Use four-particle phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    write (u, "(A)") "***********************"
    write (u, "(A)")

    x1=0.5_default; x2=0.25_default; x3=0.65_default
    write (u, "(A)" ) "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,1X,A,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (4, 2, 2, 1)
    call generator%real_kinematics%p_born_lab%set_momenta (1, p_born)

    allocate (generator%m2 (2))
    generator%m2(1) = 0._default; generator%m2(2) = 0._default
    allocate (generator%is_massive (4))
    generator%is_massive = .false.
    call generator%generate_radiation_variables ([x1,x2,x3], p_born, phs_identifiers)
    call generator%compute_xi_ref_momenta (p_born)
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       call generator%compute_xi_max (emitter, i_phs, p_born, &
            generator%real_kinematics%xi_max(i_phs))
    end do
    write (u, "(A)")  &
         "* With these, the following radiation variables have been produced:"
    associate (rad_var => generator%real_kinematics)
      write (u, "(A,F3.2)") "xi_tilde: ", rad_var%xi_tilde
      write (u, "(A,F3.2)") "y: " , rad_var%y(1)
      write (u, "(A,F3.2)") "phi: ", rad_var%phi
    end associate
    write (u, "(A)") "Initial-state momentum fractions: "
    associate (xb => generator%isr_kinematics%x)
       write (u, "(A,F3.2)") "x_born_plus: ", xb(1)
       write (u, "(A,F3.2)") "x_born_minus: ", xb(2)
    end associate
    call write_separator (u)
    write (u, "(A)") "Produce real momenta: "
    i_phs = 1; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    allocate (p_real(5))
    call generator%generate_isr (i_phs, p_born, p_real)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    call write_separator (u)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_2"

  end subroutine phs_fks_generator_2

  subroutine phs_fks_generator_3 (u)
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: p_real
    real(default) :: x1, x2, x3
    real(default) :: mB, mW, mT
    integer :: i, emitter, i_phs
    type(phs_identifier_t), dimension(2) :: phs_identifiers

    write (u, "(A)") "* Test output: phs_fks_generator_3"
    write (u, "(A)") "* Puropse: Create real phase space for particle decays"
    write (u, "(A)")

    allocate (p_born(3))
    p_born(1)%p(0) = 172._default
    p_born(1)%p(1) = 0._default
    p_born(1)%p(2) = 0._default
    p_born(1)%p(3) = 0._default
    p_born(2)%p(0) = 104.72866679_default
    p_born(2)%p(1) = 45.028053213_default
    p_born(2)%p(2) = 29.450337581_default
    p_born(2)%p(3) = -5.910229156_default
    p_born(3)%p(0) = 67.271333209_default
    p_born(3)%p(1:3) = -p_born(2)%p(1:3)

    generator%n_in = 1
    allocate (generator%isr_kinematics)
    generator%isr_kinematics%isr_mode = SQRTS_FIXED
    call generator%set_xi_and_y_bounds ()

    mB = 4.2_default
    mW = 80.376_default
    mT = 172._default

    generator%sqrts = mT

    write (u, "(A)") "* Use three-particle phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    write (u, "(A)") "**********************"
    write (u, "(A)")

    x1 = 0.5_default; x2 = 0.25_default; x3 = 0.6_default
    write (u, "(A)") "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,A,1X,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (3, 2, 2, 1)
    call generator%real_kinematics%p_born_lab%set_momenta (1, p_born)

    allocate (generator%emitters(2))
    generator%emitters(1) = 1
    generator%emitters(2) = 3
    allocate (generator%m2 (3), generator%is_massive(3))
    generator%m2(1) = mT**2
    generator%m2(2) = mW**2
    generator%m2(3) = mB**2
    generator%is_massive = .true.
    phs_identifiers(1)%emitter = 1
    phs_identifiers(2)%emitter = 3

    call generator%generate_radiation_variables ([x1,x2,x3], p_born, phs_identifiers)
    call generator%compute_xi_ref_momenta (p_born)
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       call generator%compute_xi_max (emitter, i_phs, p_born, &
            generator%real_kinematics%xi_max(i_phs))
    end do

    write (u, "(A)") &
       "* With these, the following radiation variables have been produced: "
    associate (rad_var => generator%real_kinematics)
      write (u, "(A,F4.2)") "xi_tilde: ", rad_var%xi_tilde
      do i = 1, 2
         write (u, "(A,I1,A,F5.2)") "i: ", i, "y: " , rad_var%y(i)
      end do
      write (u, "(A,F4.2)") "phi: ", rad_var%phi
    end associate

    call write_separator (u)
    write (u, "(A)") "Produce real momenta via initial-state emission: "
    i_phs = 1; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    allocate (p_real (4))
    call generator%generate_isr_fixed_beam_energy (i_phs, p_born, p_real)
    call pacify (p_real, 1E-6_default)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    call write_separator(u)
    write (u, "(A)") "Produce real momenta via final-state emisson: "
    i_phs = 2; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    call generator%generate_fsr (emitter, i_phs, p_born, p_real)
    call pacify (p_real, 1E-6_default)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_3"

  end subroutine phs_fks_generator_3

  subroutine phs_fks_generator_4 (u)
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: p_real
    integer, dimension(:), allocatable :: emitters
    integer, dimension(:,:), allocatable :: resonance_lists
    type(resonance_contributors_t), dimension(2) :: alr_contributors
    real(default) :: x1, x2, x3
    real(default), parameter :: sqrts = 250.0_default
    integer, parameter :: nlegborn = 6
    integer :: i_phs, i_con, emitter
    real(default) :: m_inv_born, m_inv_real
    character(len=7) :: fmt
    type(phs_identifier_t), dimension(2) :: phs_identifiers

    call pac_fmt (fmt, FMT_19, FMT_15, .true.)

    write (u, "(A)") "* Test output: phs_fks_generator_4"
    write (u, "(A)") "* Purpose: Create FSR phase space with fixed resonances"
    write (u, "(A)")

    allocate (p_born (nlegborn))
    p_born(1)%p(0) = 250._default
    p_born(1)%p(1) = 0._default
    p_born(1)%p(2) = 0._default
    p_born(1)%p(3) = 250._default
    p_born(2)%p(0) = 250._default
    p_born(2)%p(1) = 0._default
    p_born(2)%p(2) = 0._default
    p_born(2)%p(3) = -250._default
    p_born(3)%p(0) = 145.91184486_default
    p_born(3)%p(1) = 50.39727589_default
    p_born(3)%p(2) = 86.74156041_default
    p_born(3)%p(3) = -69.03608748_default
    p_born(4)%p(0) = 208.1064784_default
    p_born(4)%p(1) = -44.07610020_default
    p_born(4)%p(2) = -186.34264578_default
    p_born(4)%p(3) = 13.48038407_default
    p_born(5)%p(0) = 26.25614471_default
    p_born(5)%p(1) = -25.12258068_default
    p_born(5)%p(2) = -1.09540228_default
    p_born(5)%p(3) = -6.27703505_default
    p_born(6)%p(0) = 119.72553196_default
    p_born(6)%p(1) = 18.80140499_default
    p_born(6)%p(2) = 100.69648766_default
    p_born(6)%p(3) = 61.83273846_default

    allocate (generator%isr_kinematics)
    generator%n_in = 2
    generator%isr_kinematics%isr_mode = SQRTS_FIXED
    call generator%set_xi_and_y_bounds ()

    call generator%set_sqrts_hat (sqrts)

    write (u, "(A)") "* Test process: e+ e- -> W+ W- b b~"
    write (u, "(A)") "* Resonance pairs: (3,5) and (4,6)"
    write (u, "(A)") "* Use four-particle phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    write (u, "(A)") "******************************"
    write (u, "(A)")

    x1 = 0.5_default; x2 = 0.25_default; x3 = 0.75_default
    write (u, "(A)") "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,1X,A,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (nlegborn, 2, 2, 2)

    allocate (generator%emitters (2))
    generator%emitters(1) = 5; generator%emitters(2) = 6
    allocate (generator%m2 (nlegborn))
    generator%m2 = p_born**2
    allocate (generator%is_massive (nlegborn))
    generator%is_massive (1:2) = .false.
    generator%is_massive (3:6) = .true.

    phs_identifiers(1)%emitter = 5
    phs_identifiers(2)%emitter = 6
    do i_phs = 1, 2
       allocate (phs_identifiers(i_phs)%contributors (2))
    end do
    allocate (resonance_lists (2, 2))
    resonance_lists (1,:) = [3,5]
    resonance_lists (2,:) = [4,6]
    !!! Here is obviously some redundance. Surely we can improve on this.
    do i_phs = 1, 2
       phs_identifiers(i_phs)%contributors = resonance_lists(i_phs,:)
    end do
    do i_con = 1, 2
       allocate (alr_contributors(i_con)%c (size (resonance_lists(i_con,:))))
       alr_contributors(i_con)%c = resonance_lists(i_con,:)
    end do
    call generator%generate_radiation_variables &
       ([x1, x2, x3], p_born, phs_identifiers)

    allocate (p_real(nlegborn + 1))
    call generator%compute_xi_ref_momenta (p_born, alr_contributors)
    !!! Keep the distinction between i_phs and i_con because in general,
    !!! they are not the same.
    do i_phs = 1, 2
       i_con = i_phs
       emitter = phs_identifiers(i_phs)%emitter
       write (u, "(A,I1,1X,A,I1,A,I1,A)") &
            "* Generate FSR phase space for emitter ", emitter, &
            "and resonance pair (",  resonance_lists (i_con, 1), ",", &
            resonance_lists (i_con, 2), ")"
       call generator%compute_xi_max (emitter, i_phs, p_born, &
            generator%real_kinematics%xi_max(i_phs), i_con = i_con)
       call generator%generate_fsr (emitter, i_phs, i_con, p_born, p_real)
       call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
       call write_separator(u)
       write (u, "(A)") "* Check if resonance masses are conserved: "
       m_inv_born = compute_resonance_mass (p_born, resonance_lists (i_con,:))
       m_inv_real = compute_resonance_mass (p_real, resonance_lists (i_con,:), 7)
       write (u, "(A,1X, " // fmt // ")") "m_inv_born = ", m_inv_born
       write (u, "(A,1X, " // fmt // ")") "m_inv_real = ", m_inv_real
       if (abs (m_inv_born - m_inv_real) < tiny_07) then
          write (u, "(A)") " Success! "
       else
          write (u, "(A)") " Failure! "
       end if
       call write_separator(u)
       call write_separator(u)
    end do
    deallocate (p_real)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_4"
  end subroutine phs_fks_generator_4

  subroutine phs_fks_generator_5 (u)
    use ttv_formfactors, only: init_parameters
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born, pb1
    type(vector4_t), dimension(:), allocatable :: p_born_onshell, pb1_os
    type(vector4_t), dimension(:), allocatable :: p_real
    real(default) :: x1, x2, x3
    real(default) :: mB, mW, mtop, mcheck
    integer :: i, emitter, i_phs
    type(phs_identifier_t), dimension(2) :: phs_identifiers
    type(lorentz_transformation_t) :: L_to_cms
    real(default), parameter :: sqrts = 360._default
    real(default), parameter :: momentum_tolerance = 1E-10_default
    real(default) :: mpole, gam_out

    write (u, "(A)") "* Test output: phs_fks_generator_5"
    write (u, "(A)") "* Puropse: Perform threshold on-shell projection of "
    write (u, "(A)") "*          Born momenta and create a real phase-space "
    write (u, "(A)") "*          point from those. "
    write (u, "(A)")

    allocate (p_born(6), p_born_onshell(6))
    p_born(1)%p(0) = sqrts / two
    p_born(1)%p(1:2) = zero
    p_born(1)%p(3) = sqrts / two
    p_born(2)%p(0) = sqrts / two
    p_born(2)%p(1:2) = zero
    p_born(2)%p(3) = -sqrts / two
    p_born(3)%p(0) = 117.1179139230_default
    p_born(3)%p(1) = 56.91215483880_default
    p_born(3)%p(2) = -40.02386013017_default
    p_born(3)%p(3) = -49.07634310496_default
    p_born(4)%p(0) = 98.91904548743_default
    p_born(4)%p(1) = 56.02241403836_default
    p_born(4)%p(2) = -8.302977504723_default
    p_born(4)%p(3) = -10.50293716131_default
    p_born(5)%p(0) = 62.25884689208_default
    p_born(5)%p(1) = -60.00786540278_default
    p_born(5)%p(2) = 4.753602375910_default
    p_born(5)%p(3) = 15.32916731546_default
    p_born(6)%p(0) = 81.70419369751_default
    p_born(6)%p(1) = -52.92670347439_default
    p_born(6)%p(2) = 43.57323525898_default
    p_born(6)%p(3) = 44.25011295081_default

    generator%n_in = 2
    allocate (generator%isr_kinematics)
    generator%isr_kinematics%isr_mode = SQRTS_FIXED
    call generator%set_xi_and_y_bounds ()

    mB = 4.2_default
    mW = 80.376_default
    mtop = 172._default

    generator%sqrts = sqrts

    !!! Dummy-initialization of the threshold model because generate_fsr_threshold
    !!! uses m1s_to_mpole to determine if it is above or below threshold.
    call init_parameters (mpole, gam_out, mtop, one, one / 1.5_default, 125._default, &
         0.47_default, 0.118_default, 91._default, 80._default, 4.2_default, &
         one, one, one, one, zero, zero, zero, zero, zero, zero, .false., zero)

    write (u, "(A)") "* Use four-particle phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    call vector4_check_momentum_conservation &
         (p_born, 2, unit = u, abs_smallness = momentum_tolerance, verbose = .true.)
    write (u, "(A)") "**********************"
    write (u, "(A)")

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (7, 2, 2, 2)
    call generator%real_kinematics%init_onshell (7, 2)
    generator%real_kinematics%p_born_cms%phs_point(1) = p_born

    write (u, "(A)") "Get boost projection system -> CMS: "
    L_to_cms = get_boost_for_threshold_projection (p_born, sqrts, mtop)
    call L_to_cms%write (u, testflag = .true., ultra = .true.)
    write (u, "(A)") "**********************"
    write (u, "(A)")

    write (u, "(A)") "* Perform onshell-projection:"
    pb1 = generator%real_kinematics%p_born_cms%phs_point(1)
    call threshold_projection_born (mtop, L_to_cms, pb1, p_born_onshell)
    generator%real_kinematics%p_born_onshell%phs_point(1) = p_born_onshell

    call generator%real_kinematics%p_born_onshell%write &
        (1, unit = u, testflag = .true., ultra = .true.)

    pb1_os = generator%real_kinematics%p_born_onshell%phs_point(1)
    call check_phsp (pb1_os, 0)

    allocate (generator%emitters (2))
    generator%emitters(1) = THR_POS_B; generator%emitters(2) = THR_POS_BBAR

    allocate (generator%m2 (6), generator%is_massive(6))
    generator%m2 = p_born**2
    generator%is_massive (1:2) = .false.
    generator%is_massive (3:6) = .true.

    phs_identifiers(1)%emitter = THR_POS_B
    phs_identifiers(2)%emitter = THR_POS_BBAR

    x1 = 0.5_default; x2 = 0.25_default; x3 = 0.6_default
    write (u, "(A)") "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,A,1X,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3


    call generator%generate_radiation_variables ([x1,x2,x3], p_born_onshell, phs_identifiers)
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       call generator%compute_xi_ref_momenta_threshold (p_born_onshell)
       call generator%compute_xi_max (emitter, i_phs, p_born_onshell, &
            generator%real_kinematics%xi_max(i_phs), i_con = thr_leg(emitter))
    end do
    write (u, "(A)") &
       "* With these, the following radiation variables have been produced: "
    associate (rad_var => generator%real_kinematics)
      write (u, "(A,F4.2)") "xi_tilde: ", rad_var%xi_tilde
      write (u, "(A)") "xi_max: "
      write (u, "(2F5.2)") rad_var%xi_max(1), rad_var%xi_max(2)
      write (u, "(A)") "y: "
      write (u, "(2F5.2)") rad_var%y(1), rad_var%y(2)
      write (u, "(A,F4.2)") "phi: ", rad_var%phi
    end associate

    call write_separator (u)
    write (u, "(A)") "* Produce real momenta from on-shell phase space: "
    allocate (p_real(7))
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       write (u, "(A,I1)") "emitter: ", emitter
       call generator%generate_fsr_threshold (emitter, i_phs, p_born_onshell, p_real)
       call check_phsp (p_real, emitter)
    end do

    call write_separator(u)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_5"

  contains
    subroutine check_phsp (p, emitter)
       type(vector4_t), intent(inout), dimension(:) :: p
       integer, intent(in) :: emitter
       type(vector4_t) :: pp
       real(default) :: E_tot
       logical :: check
       write (u, "(A)") "* Check momentum conservation: "
       call vector4_check_momentum_conservation &
            (p, 2, unit = u, abs_smallness = momentum_tolerance, verbose = .true.)
       write (u, "(A)") "* Check invariant masses: "
       write (u, "(A)", advance = "no") "inv(W+, b, gl): "
       pp = p(THR_POS_WP) + p(THR_POS_B)
       if (emitter == THR_POS_B) pp = pp + p(THR_POS_GLUON)
       if (nearly_equal (pp**1, mtop)) then
          write (u, "(A)") "CHECK"
       else
          write (u, "(A,F7.3)") "FAIL: ", pp**1
       end if
       write (u, "(A)", advance = "no") "inv(W-, bbar): "
       pp = p(THR_POS_WM) + p(THR_POS_BBAR)
       if (emitter == THR_POS_BBAR) pp = pp + p(THR_POS_GLUON)
       if (nearly_equal (pp**1, mtop)) then
          write (u, "(A)") "CHECK"
       else
          write (u, "(A,F7.3)") "FAIL: ", pp**1
       end if
       write (u, "(A)") "* Sum of energies equal to sqrts?"
       E_tot = sum(p(1:2)%p(0)); check = nearly_equal (E_tot, sqrts)
       write (u, "(A,L1)") "Initial state: ", check
       if (.not. check) write (u, "(A,F7.3)") "E_tot: ", E_tot
       if (emitter > 0) then
          E_tot = sum(p(3:7)%p(0))
       else
          E_tot = sum(p(3:6)%p(0))
       end if
       check = nearly_equal (E_tot, sqrts)
       write (u, "(A,L1)") "Final state  : ", check
       if (.not. check) write (u, "(A,F7.3)") "E_tot: ", E_tot
       call pacify (p, 1E-6_default)
       call vector4_write_set (p, u, testflag = .true., ultra = .true.)

    end subroutine check_phsp
  end subroutine phs_fks_generator_5

  subroutine phs_fks_generator_6 (u)
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: p_real
    real(default) :: x1, x2, x3
    real(default) :: mB, mW, mT
    integer :: i, emitter, i_phs
    type(phs_identifier_t), dimension(2) :: phs_identifiers

    write (u, "(A)") "* Test output: phs_fks_generator_6"
    write (u, "(A)") "* Puropse: Create real phase space for particle decays"
    write (u, "(A)")

    allocate (p_born(4))
    p_born(1)%p(0) = 173.1_default
    p_born(1)%p(1) = zero
    p_born(1)%p(2) = zero
    p_born(1)%p(3) = zero
    p_born(2)%p(0) = 68.17074462929_default
    p_born(2)%p(1) = -37.32578717617_default
    p_born(2)%p(2) = 30.99675959336_default
    p_born(2)%p(3) = -47.70321718398_default
    p_born(3)%p(0) = 65.26639312326_default
    p_born(3)%p(1) = -1.362927648502_default
    p_born(3)%p(2) = -33.25327150840_default
    p_born(3)%p(3) = 56.14324922494_default
    p_born(4)%p(0) = 39.66286224745_default
    p_born(4)%p(1) = 38.68871482467_default
    p_born(4)%p(2) = 2.256511915049_default
    p_born(4)%p(3) = -8.440032040958_default

    generator%n_in = 1
    allocate (generator%isr_kinematics)
    generator%isr_kinematics%isr_mode = SQRTS_FIXED
    call generator%set_xi_and_y_bounds ()

    mB = 4.2_default
    mW = 80.376_default
    mT = 173.1_default

    generator%sqrts = mT

    write (u, "(A)") "* Use four-particle phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    write (u, "(A)") "**********************"
    write (u, "(A)")

    x1=0.5_default; x2=0.25_default; x3=0.6_default
    write (u, "(A)") "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,A,1X,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (3, 2, 2, 1)
    call generator%real_kinematics%p_born_lab%set_momenta (1, p_born)

    allocate (generator%emitters(2))
    generator%emitters(1) = 1
    generator%emitters(2) = 2
    allocate (generator%m2 (4), generator%is_massive(4))
    generator%m2(1) = mT**2
    generator%m2(2) = mB**2
    generator%m2(3) = zero
    generator%m2(4) = zero
    generator%is_massive(1:2) = .true.
    generator%is_massive(3:4) = .false.
    phs_identifiers(1)%emitter = 1
    phs_identifiers(2)%emitter = 2

    call generator%generate_radiation_variables ([x1,x2,x3], p_born, phs_identifiers)
    call generator%compute_xi_ref_momenta (p_born)
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       call generator%compute_xi_max (emitter, i_phs, p_born, &
            generator%real_kinematics%xi_max(i_phs))
    end do

    write (u, "(A)") &
       "* With these, the following radiation variables have been produced: "
    associate (rad_var => generator%real_kinematics)
      write (u, "(A,F4.2)") "xi_tilde: ", rad_var%xi_tilde
      do i = 1, 2
         write (u, "(A,I1,A,F5.2)") "i: ", i, "y: " , rad_var%y(i)
      end do
      write (u, "(A,F4.2)") "phi: ", rad_var%phi
    end associate

    call write_separator (u)
    write (u, "(A)") "Produce real momenta via initial-state emission: "
    i_phs = 1; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    allocate (p_real(5))
    call generator%generate_isr_fixed_beam_energy (i_phs, p_born, p_real)
    call pacify (p_real, 1E-6_default)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    call write_separator(u)
    write (u, "(A)") "Produce real momenta via final-state emisson: "
    i_phs = 2; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    call generator%generate_fsr (emitter, i_phs, p_born, p_real)
    call pacify (p_real, 1E-6_default)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_6"

  end subroutine phs_fks_generator_6

  subroutine phs_fks_generator_7 (u)
    integer, intent(in) :: u
    type(phs_fks_generator_t) :: generator
    type(vector4_t), dimension(:), allocatable :: p_born
    type(vector4_t), dimension(:), allocatable :: p_real
    real(default) :: x1, x2, x3
    integer :: i, emitter, i_phs
    type(phs_identifier_t), dimension(2) :: phs_identifiers
    real(default), parameter :: sqrts = 1000.0_default

    write (u, "(A)") "* Test output: phs_fks_generator_7"
    write (u, "(A)") "* Puropse: Create real phase space for scattering ISR"
    write (u, "(A)") "*          keeping the beam energy fixed."
    write (u, "(A)")

    allocate (p_born(4))
    p_born(1)%p(0) = 500._default
    p_born(1)%p(1) = 0._default
    p_born(1)%p(2) = 0._default
    p_born(1)%p(3) = 500._default
    p_born(2)%p(0) = 500._default
    p_born(2)%p(1) = 0._default
    p_born(2)%p(2) = 0._default
    p_born(2)%p(3) = -500._default
    p_born(3)%p(0) = 500._default
    p_born(3)%p(1) = 11.275563070_default
    p_born(3)%p(2) = -13.588797663_default
    p_born(3)%p(3) = 486.93070588_default
    p_born(4)%p(0) = 500._default
    p_born(4)%p(1:3) = -p_born(3)%p(1:3)

    phs_identifiers(1)%emitter = 1
    phs_identifiers(2)%emitter = 2

    allocate (generator%emitters(2))
    generator%n_in = 2
    allocate (generator%isr_kinematics)
    generator%isr_kinematics%isr_mode = SQRTS_FIXED
    call generator%set_xi_and_y_bounds ()
    generator%emitters(1) = 1; generator%emitters(2) = 2
    generator%sqrts = sqrts

    write (u, "(A)") "* Use 2 -> 2 phase space containing: "
    call vector4_write_set (p_born, u, testflag = .true., ultra = .true.)
    write (u, "(A)") "**********************"
    write (u, "(A)")

    x1 = 0.5_default; x2 = 0.25_default; x3 = 0.6_default
    write (u, "(A)") "* Use random numbers: "
    write (u, "(A,F3.2,1X,A,F3.2,A,1X,F3.2)") &
       "x1: ", x1, "x2: ", x2, "x3: ", x3

    allocate (generator%real_kinematics)
    call generator%real_kinematics%init (4, 2, 2, 1)
    call generator%real_kinematics%p_born_lab%set_momenta (1, p_born)

    allocate (generator%m2 (4))
    generator%m2 = 0._default
    allocate (generator%is_massive(4))
    generator%is_massive = .false.
    call generator%generate_radiation_variables ([x1,x2,x3], p_born, phs_identifiers)
    call generator%compute_xi_ref_momenta (p_born)
    do i_phs = 1, 2
       emitter = phs_identifiers(i_phs)%emitter
       call generator%compute_xi_max (emitter, i_phs, p_born, &
            generator%real_kinematics%xi_max(i_phs))
    end do

    write (u, "(A)") &
       "* With these, the following radiation variables have been produced: "
    associate (rad_var => generator%real_kinematics)
       write (u, "(A,F4.2)") "xi_tilde: ", rad_var%xi_tilde
       do i = 1, 2
          write (u, "(A,I1,A,F5.2)") "i: ", i, "y: " , rad_var%y(i)
       end do
       write (u, "(A,F4.2)") "phi: ", rad_var%phi
    end associate

    call write_separator (u)
    write (u, "(A)") "Produce real momenta via initial-state emission: "
    i_phs = 1; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    allocate (p_real(5))
    call generator%generate_isr_fixed_beam_energy (i_phs, p_born, p_real)
    call pacify (p_real, 1E-6_default)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    call write_separator(u)
    i_phs = 2; emitter = phs_identifiers(i_phs)%emitter
    write (u, "(A,I1)") "emitter: ", emitter
    call generator%generate_isr_fixed_beam_energy (i_phs, p_born, p_real)
    call pacify (p_real, 1E-6_default)
    call vector4_write_set (p_real, u, testflag = .true., ultra = .true.)
    write (u, "(A)")
    write (u, "(A)") "* Test output end: phs_fks_generator_7"

  end subroutine phs_fks_generator_7


end module phs_fks_uti
