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

submodule (mci_base) mci_base_s

  use io_units
  use format_utils, only: pac_fmt
  use format_defs, only: FMT_14, FMT_17
  use diagnostics

  implicit none

contains

  module subroutine mci_final (object)
    class(mci_t), intent(inout) :: object
    if (allocated (object%rng))  call object%rng%final ()
  end subroutine mci_final

  module subroutine mci_write (object, unit, pacify, md5sum_version)
    class(mci_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: pacify
    logical, intent(in), optional :: md5sum_version
    logical :: md5sum_ver
    integer :: u, i, j
    character(len=7) :: fmt
    call pac_fmt (fmt, FMT_17, FMT_14, pacify)
    u = given_output_unit (unit)
    md5sum_ver = .false.
    if (present (md5sum_version))  md5sum_ver = md5sum_version
    if (object%use_timer .and. .not. md5sum_ver) then
       write (u, "(2x)", advance="no")
       call object%timer%write (u)
    end if
    if (object%integral_known) then
       write (u, "(3x,A," // fmt // ")") &
            "Integral             = ", object%integral
    end if
    if (object%error_known) then
       write (u, "(3x,A," // fmt // ")") &
            "Error                = ", object%error
    end if
    if (object%efficiency_known) then
       write (u, "(3x,A," // fmt // ")")  &
            "Efficiency           = ", object%efficiency
    end if
    write (u, "(3x,A,I0)")  "Number of channels   = ", object%n_channel
    write (u, "(3x,A,I0)")  "Number of dimensions = ", object%n_dim
    if (object%n_chain > 0) then
       write (u, "(3x,A,I0)")  "Number of chains     = ", object%n_chain
       write (u, "(3x,A)")  "Chains:"
       do i = 1, object%n_chain
          write (u, "(5x,I0,':')", advance = "no")  i
          do j = 1, object%n_channel
             if (object%chain(j) == i) &
                  write (u, "(1x,I0)", advance = "no")  j
          end do
          write (u, "(A)")
       end do
    end if
  end subroutine mci_write

  module subroutine mci_startup_message (mci, unit, n_calls)
    class(mci_t), intent(in) :: mci
    integer, intent(in), optional :: unit, n_calls
    if (mci%n_chain > 0) then
       write (msg_buffer, "(A,3(1x,I0,1x,A))") &
            "Integrator:", mci%n_chain, "chains,", &
            mci%n_channel, "channels,", &
            mci%n_dim, "dimensions"
    else
       write (msg_buffer, "(A,3(1x,I0,1x,A))") &
            "Integrator:", &
            mci%n_channel, "channels,", &
            mci%n_dim, "dimensions"
    end if
    call msg_message (unit = unit)
  end subroutine mci_startup_message

  module subroutine mci_record_index (mci, i_mci)
    class(mci_t), intent(inout) :: mci
    integer, intent(in) :: i_mci
  end subroutine mci_record_index

  module subroutine mci_set_dimensions (mci, n_dim, n_channel)
    class(mci_t), intent(inout) :: mci
    integer, intent(in) :: n_dim
    integer, intent(in) :: n_channel
    mci%n_dim = n_dim
    mci%n_channel = n_channel
  end subroutine mci_set_dimensions

  module subroutine mci_declare_chains (mci, chain)
    class(mci_t), intent(inout) :: mci
    integer, dimension(:), intent(in) :: chain
    allocate (mci%chain (size (chain)))
    mci%n_chain = maxval (chain)
    allocate (mci%chain_weights (mci%n_chain), source = 0._default)
    mci%chain = chain
  end subroutine mci_declare_chains

  module subroutine mci_collect_chain_weights (mci, weight)
    class(mci_t), intent(inout) :: mci
    real(default), dimension(:), intent(in) :: weight
    integer :: i, c
    if (allocated (mci%chain)) then
       mci%chain_weights = 0
       do i = 1, size (mci%chain)
          c = mci%chain(i)
          mci%chain_weights(c) = mci%chain_weights(c) + weight(i)
       end do
    end if
  end subroutine mci_collect_chain_weights

  module function mci_has_chains (mci) result (flag)
    class(mci_t), intent(in) :: mci
    logical :: flag
    flag = allocated (mci%chain)
  end function mci_has_chains

  module subroutine mci_write_chain_weights (mci, unit)
    class(mci_t), intent(in) :: mci
    integer, intent(in), optional :: unit
    integer :: u, i, n, n_digits
    character(4) :: ifmt
    u = given_output_unit (unit)
    if (allocated (mci%chain_weights)) then
       write (u, "(1x,A)")  "Weights of channel chains (groves):"
       n_digits = 0
       n = size (mci%chain_weights)
       do while (n > 0)
          n = n / 10
          n_digits = n_digits + 1
       end do
       write (ifmt, "(A1,I1)") "I", n_digits
       do i = 1, size (mci%chain_weights)
          write (u, "(3x," // ifmt // ",F13.10)")  i, mci%chain_weights(i)
       end do
    end if
  end subroutine mci_write_chain_weights

  module subroutine mci_set_md5sum (mci, md5sum)
    class(mci_t), intent(inout) :: mci
    character(32), intent(in) :: md5sum
    mci%md5sum = md5sum
  end subroutine mci_set_md5sum

  module subroutine mci_add_pass (mci, adapt_grids, adapt_weights, final_pass)
    class(mci_t), intent(inout) :: mci
    logical, intent(in), optional :: adapt_grids
    logical, intent(in), optional :: adapt_weights
    logical, intent(in), optional :: final_pass
  end subroutine mci_add_pass

  module subroutine mci_import_rng (mci, rng)
    class(mci_t), intent(inout) :: mci
    class(rng_t), intent(inout), allocatable :: rng
    call move_alloc (rng, mci%rng)
  end subroutine mci_import_rng

  module subroutine mci_set_timer (mci, active)
    class(mci_t), intent(inout) :: mci
    logical, intent(in) :: active
    mci%use_timer = active
  end subroutine mci_set_timer

  module subroutine mci_start_timer (mci)
    class(mci_t), intent(inout) :: mci
    if (mci%use_timer)  call mci%timer%start ()
  end subroutine mci_start_timer

  module subroutine mci_stop_timer (mci)
    class(mci_t), intent(inout) :: mci
    if (mci%use_timer)  call mci%timer%stop ()
  end subroutine mci_stop_timer

  module subroutine mci_sampler_test (mci, sampler, n_calls)
    class(mci_t), intent(inout) :: mci
    class(mci_sampler_t), intent(inout), target :: sampler
    integer, intent(in) :: n_calls
    real(default), dimension(:), allocatable :: x_in, f
    real(default), dimension(:,:), allocatable :: x_out
    real(default) :: val
    integer :: i, c
    allocate (x_in (mci%n_dim))
    allocate (f (mci%n_channel))
    allocate (x_out (mci%n_dim, mci%n_channel))
    do i = 1, n_calls
       c = mod (i, mci%n_channel) + 1
       call mci%rng%generate_array (x_in)
       call sampler%evaluate (c, x_in, val, x_out, f)
    end do
  end subroutine mci_sampler_test

  module subroutine mci_pacify (object, efficiency_reset, error_reset)
    class(mci_t), intent(inout) :: object
    logical, intent(in), optional :: efficiency_reset, error_reset
  end subroutine mci_pacify

  module function mci_get_integral (mci) result (integral)
    class(mci_t), intent(in) :: mci
    real(default) :: integral
    if (mci%integral_known) then
       integral = mci%integral
    else
       call msg_bug ("The integral is unknown. This is presumably a" // &
            "WHIZARD bug.")
    end if
  end function mci_get_integral

  module function mci_get_error (mci) result (error)
    class(mci_t), intent(in) :: mci
    real(default) :: error
    if (mci%error_known) then
       error = mci%error
    else
       error = 0
    end if
  end function mci_get_error

  module function mci_get_efficiency (mci) result (efficiency)
    class(mci_t), intent(in) :: mci
    real(default) :: efficiency
    if (mci%efficiency_known) then
       efficiency = mci%efficiency
    else
       efficiency = 0
    end if
  end function mci_get_efficiency

  module function mci_get_time (mci) result (time)
    class(mci_t), intent(in) :: mci
    real(default) :: time
    if (mci%use_timer) then
       time = mci%timer
    else
       time = 0
    end if
  end function mci_get_time

  pure module function mci_get_md5sum (mci) result (md5sum)
    class(mci_t), intent(in) :: mci
    character(32) :: md5sum
    md5sum = mci%md5sum
  end function mci_get_md5sum

  module subroutine mci_instance_base_init (mci_instance, mci)
    class(mci_instance_t), intent(out) :: mci_instance
    class(mci_t), intent(in), target :: mci
    allocate (mci_instance%w (mci%n_channel))
    allocate (mci_instance%f (mci%n_channel))
    allocate (mci_instance%x (mci%n_dim, mci%n_channel))
    if (mci%n_channel > 0) then
       call mci_instance%set_channel_weights &
            (spread (1._default, dim=1, ncopies=mci%n_channel))
    end if
    mci_instance%f = 0
    mci_instance%x = 0
  end subroutine mci_instance_base_init

  module subroutine mci_instance_set_channel_weights &
       (mci_instance, weights, sum_non_zero)
    class(mci_instance_t), intent(inout) :: mci_instance
    real(default), dimension(:), intent(in) :: weights
    logical, intent(out), optional :: sum_non_zero
    real(default) :: wsum
    wsum = sum (weights)
    if (wsum /= 0) then
       mci_instance%w = weights / wsum
       if (present (sum_non_zero)) sum_non_zero = .true.
    else
       if (present (sum_non_zero)) sum_non_zero = .false.
       call msg_warning ("MC sampler initialization:&
            & sum of channel weights is zero")
    end if
  end subroutine mci_instance_set_channel_weights

  module subroutine mci_instance_evaluate (mci, sampler, c, x)
    class(mci_instance_t), intent(inout) :: mci
    class(mci_sampler_t), intent(inout) :: sampler
    integer, intent(in) :: c
    real(default), dimension(:), intent(in) :: x
    real(default) :: val
    call sampler%evaluate (c, x, val, mci%x, mci%f)
    mci%valid = sampler%is_valid ()
    if (mci%valid) then
       call mci%compute_weight (c)
       call mci%record_integrand (val)
    end if
  end subroutine mci_instance_evaluate

  module subroutine mci_instance_fetch (mci, sampler, c)
    class(mci_instance_t), intent(inout) :: mci
    class(mci_sampler_t), intent(in) :: sampler
    integer, intent(in) :: c
    real(default) :: val
    mci%valid = sampler%is_valid ()
    if (mci%valid) then
       call sampler%fetch (val, mci%x, mci%f)
       call mci%compute_weight (c)
       call mci%record_integrand (val)
    end if
  end subroutine mci_instance_fetch

  module function mci_instance_get_value (mci) result (value)
    class(mci_instance_t), intent(in) :: mci
    real(default) :: value
    if (mci%valid) then
       value = mci%integrand * mci%mci_weight
    else
       value = 0
    end if
  end function mci_instance_get_value

  module function mci_instance_get_n_event_dropped (mci) result (n_dropped)
    class(mci_instance_t), intent(in) :: mci
    integer :: n_dropped
    n_dropped = mci%n_dropped
  end function mci_instance_get_n_event_dropped

  module subroutine mci_instance_reset_n_event_dropped (mci)
    class(mci_instance_t), intent(inout) :: mci
    mci%n_dropped = 0
  end subroutine mci_instance_reset_n_event_dropped

  module subroutine mci_instance_record_event_dropped (mci)
    class(mci_instance_t), intent(inout) :: mci
    mci%n_dropped = mci%n_dropped + 1
  end subroutine mci_instance_record_event_dropped

  module subroutine mci_state_write (object, unit)
    class(mci_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "MCI state:"
    write (u, "(3x,A,I0)")  "Channel   = ", object%selected_channel
    write (u, "(3x,A,999(1x,F12.10))")  "x (in)    =", object%x_in
    write (u, "(3x,A,ES19.12)")  "Integrand = ", object%val
  end subroutine mci_state_write

  module subroutine mci_instance_store (mci, state)
    class(mci_instance_t), intent(in) :: mci
    class(mci_state_t), intent(out) :: state
    state%selected_channel = mci%selected_channel
    allocate (state%x_in (size (mci%x, 1)))
    state%x_in = mci%x(:,mci%selected_channel)
    state%val = mci%integrand
  end subroutine mci_instance_store

  module subroutine mci_instance_recall (mci, sampler, state)
    class(mci_instance_t), intent(inout) :: mci
    class(mci_sampler_t), intent(inout) :: sampler
    class(mci_state_t), intent(in) :: state
    if (size (state%x_in) == size (mci%x, 1) &
         .and. state%selected_channel <= size (mci%x, 2)) then
       call sampler%rebuild (state%selected_channel, &
            state%x_in, state%val, mci%x, mci%f)
       call mci%compute_weight (state%selected_channel)
       call mci%record_integrand (state%val)
    else
       call msg_fatal ("Recalling event: mismatch in channel or dimension")
    end if
  end subroutine mci_instance_recall


end submodule mci_base_s

