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

module muli
  use, intrinsic :: iso_fortran_env
  use kinds, only: default

  implicit none
  private



  public :: muli_t

  type :: muli_t     
     real(default) :: GeV2_scale_cutoff
     logical :: initialized = .false.
     integer, dimension(4) :: flow
   contains
     generic :: initialize => muli_initialize
     procedure :: muli_initialize
     procedure :: apply_initial_interaction => muli_apply_initial_interaction
     procedure :: restart => muli_restart
     procedure :: is_initialized => muli_is_initialized
     procedure :: generate_gev2_pt2 => muli_generate_gev2_pt2
     procedure :: generate_partons => muli_generate_partons
     procedure :: generate_color_flows => muli_generate_color_flows
     procedure :: get_color_flow => muli_get_color_flow
     procedure :: get_color_correlations => muli_get_color_correlations
     procedure :: replace_parton => muli_replace_parton
     procedure :: generate_next_scale => muli_generate_next_scale
  end type muli_t


contains

  subroutine muli_initialize (this, GeV2_scale_cutoff, gev2_s, &
       muli_dir, random_seed)
    class(muli_t), intent(out) :: this
    real(kind=default), intent(in) :: gev2_s, GeV2_scale_cutoff
    character(*), intent(in) :: muli_dir
    integer, intent(in), optional :: random_seed
    this%initialized = .true.
  end subroutine muli_initialize

  subroutine muli_apply_initial_interaction (this, GeV2_s, &
       x1, x2, pdg_f1, pdg_f2, n1, n2)
    class(muli_t), intent(inout) :: this
    real(default), intent(in) :: Gev2_s, x1, x2
    integer, intent(in):: pdg_f1, pdg_f2, n1, n2
  end subroutine muli_apply_initial_interaction

  subroutine muli_restart (this)
    class(muli_t), intent(inout) :: this
  end subroutine muli_restart

  elemental function muli_is_initialized (this) result (res)
    logical :: res
    class(muli_t), intent(in) :: this
    res = this%initialized
  end function muli_is_initialized

  subroutine muli_generate_gev2_pt2 (this, gev2_start_scale, gev2_new_scale)
    class(muli_t), intent(inout) :: this
    real(kind=default), intent(in) :: gev2_start_scale
    real(kind=default), intent(out) :: gev2_new_scale
    call this%generate_next_scale ()
    gev2_new_scale = 1
  end subroutine muli_generate_gev2_pt2

  subroutine muli_generate_partons (this, n1, n2, x_proton_1, x_proton_2, &
       pdg_f1, pdg_f2, pdg_f3, pdg_f4)
    class(muli_t), intent(inout) :: this
    integer, intent(in) :: n1, n2
    real(kind=default), intent(out) :: x_proton_1, x_proton_2
    integer, intent(out) :: pdg_f1, pdg_f2, pdg_f3, pdg_f4
    integer, dimension(4) :: pdg_f
    call this%generate_color_flows ()
    pdg_f = 0
    pdg_f1 = 0
    pdg_f2 = 0
    pdg_f3 = 0
    pdg_f4 = 0
  end subroutine muli_generate_partons

  subroutine muli_generate_color_flows (this)
    class(muli_t), intent(inout) :: this
    integer, dimension(4) :: flow
  end subroutine muli_generate_color_flows

  pure function muli_get_color_flow (this) result (flow)
    class(muli_t), intent(in) :: this
    integer, dimension(4) :: flow
    flow = this%flow
  end function muli_get_color_flow

  subroutine muli_get_color_correlations &
       (this, start_index, final_index, flow)
    class(muli_t), intent(in) :: this
    integer, intent(in) :: start_index
    integer, intent(out) :: final_index
    integer, dimension(2,4), intent(out) :: flow
    integer :: pos, f_end, f_beginning
    final_index = start_index
    flow = reshape([0,0,0,0,0,0,0,0],[2,4])
  end subroutine muli_get_color_correlations

  subroutine muli_replace_parton &
       (this, proton_id, old_id, new_id, pdg_f, x_proton, gev_scale)
    class(muli_t), intent(inout) :: this
    integer, intent(in) :: proton_id, old_id, new_id, pdg_f
    real(kind=default), intent(in) :: x_proton, gev_scale
  end subroutine muli_replace_parton

  subroutine muli_generate_next_scale (this, integrand_kind)
    class(muli_t), intent(inout) :: this
    integer, intent(in), optional :: integrand_kind
  end subroutine muli_generate_next_scale


end module muli
