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

submodule (ckkw_matching) ckkw_matching_s

  use debug_master, only: debug_on
  use io_units
  use format_utils, only: write_separator
  use diagnostics
  use physics_defs
  use shower_core

  implicit none

contains

  module subroutine ckkw_matching_settings_init (settings, var_list)
    class(ckkw_matching_settings_t), intent(out) :: settings
    type(var_list_t), intent(in) :: var_list
    settings%alphaS = 1.0_default
    settings%Qmin = 1.0_default
    settings%n_max_jets = 3
  end subroutine ckkw_matching_settings_init

  module subroutine ckkw_matching_settings_write (settings, unit)
    class(ckkw_matching_settings_t), intent(in) :: settings
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)")  "CKKW matching settings:"
    call write_separator (u)
    write (u, "(3x,A,1x,ES19.12)") &
         "alphaS       = ", settings%alphaS
    write (u, "(3x,A,1x,ES19.12)") &
         "Qmin         = ", settings%Qmin
    write (u, "(3x,A,1x,I0)") &
         "n_max_jets   = ", settings%n_max_jets
  end subroutine ckkw_matching_settings_write

  module subroutine ckkw_pseudo_shower_weights_init (weights)
    class(ckkw_pseudo_shower_weights_t), intent(out) :: weights
    weights%alphaS = zero
  end subroutine ckkw_pseudo_shower_weights_init

  module subroutine ckkw_pseudo_shower_weights_write (weights, unit)
    class(ckkw_pseudo_shower_weights_t), intent(in) :: weights
    integer, intent(in), optional :: unit
    integer :: s, i, u
    u = given_output_unit (unit); if (u < 0) return
    s = size (weights%weights)
    write (u, "(1x,A)")  "CKKW (pseudo) shower weights: "
    do i = 1, s
       write (u, "(3x,I0,2(ES19.12))")  i, weights%weights(i), &
            weights%weights_by_type(i,:)
    end do
    write (u, "(3x,A,1x,I0)")  "alphaS =", weights%alphaS
  end subroutine ckkw_pseudo_shower_weights_write

  pure module subroutine ckkw_pseudo_shower_weights_fake (weights, particle_set)
    class(ckkw_pseudo_shower_weights_t), intent(inout) :: weights
    type(particle_set_t), intent(in) :: particle_set
    integer :: i, j, n
    type(vector4_t) :: momentum
    n = 2**particle_set%n_tot
    if (allocated (weights%weights)) then
       deallocate (weights%weights)
    end if
    allocate (weights%weights (1:n))
    do i = 1, n
       momentum = vector4_null
       do j = 1, particle_set%n_tot
          if (btest (i,j-1)) then
             momentum = momentum + particle_set%prt(j)%p
          end if
       end do
       if (momentum**1 > 0.0) then
          weights%weights(i) = 1.0 / (momentum**2)
       end if
    end do
    ! equally distribute the weights by type
    if (allocated (weights%weights_by_type)) then
       deallocate (weights%weights_by_type)
    end if
    allocate (weights%weights_by_type (1:n, 0:4))
    do i = 1, n
       do j = 0, 4
          weights%weights_by_type(i,j) = 0.2 * weights%weights(i)
       end do
    end do
  end subroutine ckkw_pseudo_shower_weights_fake

  module subroutine ckkw_matching_init (matching, var_list, process_name)
    class(ckkw_matching_t), intent(out) :: matching
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_name
    if (debug_on) call msg_debug (D_MATCHING, "matching_init")
    call matching%settings%init (var_list)
    matching%process_name = process_name
  end subroutine ckkw_matching_init

  module subroutine ckkw_matching_write (matching, unit)
    class(ckkw_matching_t), intent(in) :: matching
    integer, intent(in), optional :: unit
    call matching%settings%write (unit)
    call matching%weights%write (unit)
  end subroutine ckkw_matching_write

  module function ckkw_matching_get_method (matching) result (method)
    type(string_t) :: method
    class(ckkw_matching_t), intent(in) :: matching
    method = matching_method (MATCH_CKKW)
  end function ckkw_matching_get_method

  module subroutine ckkw_matching_before_shower &
       (matching, particle_set, vetoed)
    class(ckkw_matching_t), intent(inout) :: matching
    type(particle_set_t), intent(inout) :: particle_set
    logical, intent(out) :: vetoed
    call matching%weights%init ()
    call matching%weights%fake (particle_set)
    select type (shower => matching%shower)
    type is (shower_t)
       call ckkw_matching_apply (shower%partons, &
            matching%settings, &
            matching%weights, matching%rng, vetoed)
    class default
       call msg_bug ("CKKW matching only works with WHIZARD shower.")
    end select
  end subroutine ckkw_matching_before_shower

  module subroutine ckkw_matching_apply &
       (partons, settings, weights, rng, vetoed)
    type(parton_pointer_t), dimension(:), intent(inout), allocatable :: &
         partons
    type(ckkw_matching_settings_t), intent(in) :: settings
    type(ckkw_pseudo_shower_weights_t), intent(in) :: weights
    class(rng_t), intent(inout), allocatable :: rng
    logical, intent(out) :: vetoed

    real(default), dimension(:), allocatable :: scales
    real(double) :: weight, sf
    real(default) :: rand
    integer :: i, n_partons

    if (signal_is_pending ()) return
    weight = one

    n_partons = size (partons)

    do i = 1, n_partons
       call partons(i)%p%write ()
    end do

    !!! the pseudo parton shower is already simulated by shower_add_interaction
    !!! get the respective clustering scales
    allocate (scales (1:n_partons))
    do i = 1, n_partons
       if (.not. associated (partons(i)%p)) cycle
       if (partons(i)%p%type == INTERNAL) then
          scales(i) = two * min (partons(i)%p%child1%momentum%p(0),  &
                                 partons(i)%p%child2%momentum%p(0))**2 * &
               (1.0 - (space_part (partons(i)%p%child1%momentum) * &
                space_part (partons(i)%p%child2%momentum)) / &
               (space_part (partons(i)%p%child1%momentum)**1 * &
                space_part (partons(i)%p%child2%momentum)**1))
          scales(i) = sqrt (scales(i))
          partons(i)%p%ckkwscale = scales(i)
          print *, scales(i)
       end if
    end do

    print *, " scales finished"
    !!! if (highest multiplicity) -> reweight with PDF(mu_F) / PDF(mu_cut)
    do i = 1, n_partons
       call partons(i)%p%write ()
    end do

    !!! Reweight and possibly veto the whole event

    !!! calculate the relative alpha_S weight

    !! calculate the Sudakov weights for internal lines
    !! calculate the Sudakov weights for external lines
    do i = 1, n_partons
       if (signal_is_pending ()) return
       if (.not. associated (partons(i)%p)) cycle
       if (partons(i)%p%type == INTERNAL) then
          !!! get type
          !!! check that all particles involved are colored
          if ((partons(i)%p%is_colored () .or. &
               partons(i)%p%ckkwtype > 0) .and. &
               (partons(i)%p%child1%is_colored () .or. &
               partons(i)%p%child1%ckkwtype > 0) .and. &
               (partons(i)%p%child1%is_colored () .or. &
               partons(i)%p%child1%ckkwtype > 0)) then
             print *, "reweight with alphaS(" , partons(i)%p%ckkwscale, &
                  ") for particle ", partons(i)%p%nr
             if (partons(i)%p%belongstoFSR) then
                print *, "FSR"
                weight = weight * D_alpha_s_fsr (partons(i)%p%ckkwscale**2, &
                     partons(i)%p%settings) / settings%alphas
             else
                print *, "ISR"
                weight = weight * &
                     D_alpha_s_isr (partons(i)%p%ckkwscale**2, &
                     partons(i)%p%settings) / settings%alphas
             end if
          else
             print *, "no reweight with alphaS for ", partons(i)%p%nr
          end if
          if (partons(i)%p%child1%type == INTERNAL) then
             print *, "internal line from ", &
                  partons(i)%p%child1%ckkwscale, &
                  " to ", partons(i)%p%ckkwscale, &
                  " for type ", partons(i)%p%child1%ckkwtype
             if (partons(i)%p%child1%ckkwtype == 0) then
                sf = 1.0
             else if (partons(i)%p%child1%ckkwtype == 1) then
                sf = SudakovQ (partons(i)%p%child1%ckkwscale, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFQ = ", sf
             else if (partons(i)%p%child1%ckkwtype == 2) then
                sf = SudakovG (partons(i)%p%child1%ckkwscale, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFG = ", sf
             else
                print *, "SUSY not yet implemented"
             end if
             weight = weight * min (one, sf)
          else
             print *, "external line from ", settings%Qmin, &
                  partons(i)%p%ckkwscale
             if (partons(i)%p%child1%is_quark ()) then
                sf = SudakovQ (settings%Qmin, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFQ = ", sf
             else if (partons(i)%p%child1%is_gluon ()) then
                sf = SudakovG (settings%Qmin, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFG = ", sf
             else
                print *, "not yet implemented (", &
                     partons(i)%p%child2%type, ")"
                sf = one
             end if
             weight = weight * min (one, sf)
          end if
          if (partons(i)%p%child2%type == INTERNAL) then
             print *, "internal line from ", partons(i)%p%child2%ckkwscale, &
                  " to ", partons(i)%p%ckkwscale, &
                  " for type ", partons(i)%p%child2%ckkwtype
             if (partons(i)%p%child2%ckkwtype == 0) then
                sf = 1.0
             else if (partons(i)%p%child2%ckkwtype == 1) then
                sf = SudakovQ (partons(i)%p%child2%ckkwscale, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFQ = ", sf
             else if (partons(i)%p%child2%ckkwtype == 2) then
                sf = SudakovG (partons(i)%p%child2%ckkwscale, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFG = ", sf
             else
                print *, "SUSY not yet implemented"
             end if
             weight = weight * min (one, sf)
          else
             print *, "external line from ", settings%Qmin, &
                  partons(i)%p%ckkwscale
             if (partons(i)%p%child2%is_quark ()) then
                sf = SudakovQ (settings%Qmin, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFQ = ", sf
             else if (partons(i)%p%child2%is_gluon ()) then
                sf = SudakovG (settings%Qmin, &
                     partons(i)%p%ckkwscale, &
                     partons(i)%p%settings, .true., rng)
                print *, "SFG = ", sf
             else
                print *, "not yet implemented (", &
                     partons(i)%p%child2%type, ")"
                sf = one
             end if
             weight = weight * min (one, sf)
          end if
       end if
    end do

    call rng%generate (rand)

    print *, "final weight: ", weight

    !!!!!!! WRONG
    vetoed = .false.
    ! vetoed = (rand > weight)
    if (vetoed) then
       return
    end if

    !!! finally perform the parton shower
    !!! veto emissions that are too hard

    deallocate (scales)
  end subroutine ckkw_matching_apply

  module subroutine ckkw_matching_after_shower (matching, particle_set, vetoed)
    class(ckkw_matching_t), intent(inout) :: matching
    type(particle_set_t), intent(inout) :: particle_set
    logical, intent(out) :: vetoed
    vetoed = .false.
  end subroutine ckkw_matching_after_shower

  function GammaQ (smallq, largeq, settings, fsr) result (gamma)
    real(default), intent(in) :: smallq, largeq
    type(shower_settings_t), intent(in) :: settings
    logical, intent(in) :: fsr
    real(default) :: gamma
    gamma = (8._default / three) / (pi * smallq)
    gamma = gamma * (log(largeq / smallq) - 0.75)
    if (fsr) then
       gamma = gamma * D_alpha_s_fsr (smallq**2, settings)
    else
       gamma = gamma * D_alpha_s_isr (smallq**2, settings)
    end if
  end function GammaQ

  function GammaG (smallq, largeq, settings, fsr) result (gamma)
    real(default), intent(in) :: smallq, largeq
    type(shower_settings_t), intent(in) :: settings
    logical, intent(in) :: fsr
    real(default) :: gamma
    gamma = 6._default / (pi * smallq)
    gamma = gamma *( log(largeq / smallq) - 11.0 / 12.0)
    if (fsr) then
       gamma = gamma * D_alpha_s_fsr (smallq**2, settings)
    else
       gamma = gamma * D_alpha_s_isr (smallq**2, settings)
    end if
  end function GammaG

  function GammaF (smallq, settings, fsr) result (gamma)
    real(default), intent(in) :: smallq
    type(shower_settings_t), intent(in) :: settings
    logical, intent(in) :: fsr
    real(default) :: gamma
    gamma = number_of_flavors (smallq, settings%max_n_flavors, &
         settings%min_virtuality) / (three * pi * smallq)
    if (fsr) then
       gamma = gamma * D_alpha_s_fsr (smallq**2, settings)
    else
       gamma = gamma * D_alpha_s_isr (smallq**2, settings)
    end if
  end function GammaF

  function SudakovQ (Q1, Q, settings, fsr, rng) result (sf)
    real(default), intent(in) :: Q1, Q
    type(shower_settings_t), intent(in) :: settings
    class(rng_t), intent(inout), allocatable :: rng
    logical, intent(in) :: fsr
    real(default) :: sf
    real(default) :: integral
    integer, parameter :: NTRIES = 100
    integer :: i
    real(default) :: rand
    integral = zero
    do i = 1, NTRIES
       call rng%generate (rand)
       integral = integral + GammaQ (Q1 + rand * (Q - Q1), Q, settings, fsr)
    end do
    integral = integral / NTRIES
    sf = exp (-integral)
  end function SudakovQ

  function SudakovG (Q1, Q, settings, fsr, rng) result (sf)
    real(default), intent(in) :: Q1, Q
    type(shower_settings_t), intent(in) :: settings
    logical, intent(in) :: fsr
    real(default) :: sf
    real(default) :: integral
    class(rng_t), intent(inout), allocatable :: rng
    integer, parameter :: NTRIES = 100
    integer :: i
    real(default) :: rand
    integral = zero
    do i = 1, NTRIES
       call rng%generate (rand)
       integral = integral + &
            GammaG (Q1 + rand * (Q - Q1), Q, settings, fsr) + &
            GammaF (Q1 + rand * (Q - Q1), settings, fsr)
    end do
    integral = integral / NTRIES
    sf = exp (-integral)
  end function SudakovG


end submodule ckkw_matching_s

