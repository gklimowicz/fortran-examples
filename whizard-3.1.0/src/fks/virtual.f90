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

module virtual

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use pdg_arrays
  use models
  use model_data, only: model_data_t, field_data_t
  use lorentz
  use nlo_data, only: nlo_settings_t
  use fks_regions

  implicit none
  private

  public :: virtual_t

  type :: virtual_t
     type(nlo_settings_t) :: settings
     real(default), dimension(:,:,:), allocatable :: gamma_0, gamma_p, c_flv
     real(default) :: fac_scale
     real(default), allocatable :: ren_scale
     real(default), allocatable :: es_scale2
     integer, dimension(:), allocatable :: n_is_neutrinos
     integer :: n_in, n_legs, n_flv
     logical :: bad_point = .false.
     type(string_t) :: selection
     real(default), dimension(:), allocatable :: sqme_born
     real(default), dimension(:), allocatable :: sqme_virt_fin
     real(default), dimension(:,:,:), allocatable :: sqme_color_c
     real(default), dimension(:,:,:), allocatable :: sqme_charge_c
     logical :: has_pdfs = .false.
   contains
     procedure :: init => virtual_init
     procedure :: init_constants => virtual_init_constants
     procedure :: set_ren_scale => virtual_set_ren_scale
     procedure :: set_fac_scale => virtual_set_fac_scale
     procedure :: set_ellis_sexton_scale => virtual_set_ellis_sexton_scale
     procedure :: evaluate => virtual_evaluate
     procedure :: compute_eikonals => virtual_compute_eikonals
     procedure :: compute_eikonals_threshold => virtual_compute_eikonals_threshold
     procedure :: set_bad_point => virtual_set_bad_point
     procedure :: evaluate_initial_state => virtual_evaluate_initial_state
     procedure :: compute_collinear_contribution &
        => virtual_compute_collinear_contribution
     procedure :: compute_massive_self_eikonals => &
          virtual_compute_massive_self_eikonals
     procedure :: final => virtual_final
  end type virtual_t


  interface
    module subroutine virtual_init &
         (virt, flv_born, n_in, settings, model, has_pdfs)
      class(virtual_t), intent(inout) :: virt
      integer, intent(in), dimension(:,:) :: flv_born
      integer, intent(in) :: n_in
      type(nlo_settings_t), intent(in) :: settings
      class(model_data_t), intent(in) :: model
      logical, intent(in) :: has_pdfs
    end subroutine virtual_init
    module subroutine virtual_init_constants (virt, flv_born, nf_input, model)
      class(virtual_t), intent(inout) :: virt
      integer, intent(in), dimension(:,:) :: flv_born
      integer, intent(in) :: nf_input
      type(string_t), dimension(2) :: corr_type
      class(model_data_t), intent(in) :: model
    end subroutine virtual_init_constants
    module subroutine virtual_set_ren_scale (virt, ren_scale)
      class(virtual_t), intent(inout) :: virt
      real(default), allocatable, intent(in) :: ren_scale
    end subroutine virtual_set_ren_scale
    module subroutine virtual_set_fac_scale (virt, p, fac_scale)
      class(virtual_t), intent(inout) :: virt
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in), optional :: fac_scale
    end subroutine virtual_set_fac_scale
    module subroutine virtual_set_ellis_sexton_scale (virt, Q)
      class(virtual_t), intent(inout) :: virt
      real(default), allocatable, intent(in) :: Q
    end subroutine virtual_set_ellis_sexton_scale
    module subroutine virtual_evaluate (virt, reg_data, alpha_coupling, &
         p_born, separate_uborns, sqme_virt)
      class(virtual_t), intent(inout) :: virt
      type(region_data_t), intent(in) :: reg_data
      real(default), dimension(2), intent(in) :: alpha_coupling
      type(vector4_t), intent(in), dimension(:)  :: p_born
      logical, intent(in) :: separate_uborns
      real(default), dimension(:), intent(inout) :: sqme_virt
    end subroutine virtual_evaluate
    module subroutine virtual_compute_eikonals (virtual, i_flv, i_corr, &
         p_born, s_o_Q2, reg_data, BI)
      class(virtual_t), intent(inout) :: virtual
      integer, intent(in) :: i_flv, i_corr
      type(vector4_t), intent(in), dimension(:)  :: p_born
      real(default), intent(in) :: s_o_Q2
      type(region_data_t), intent(in) :: reg_data
      real(default), intent(inout), dimension(:) :: BI
    end subroutine virtual_compute_eikonals
    module subroutine virtual_compute_eikonals_threshold (virtual, i_flv, &
           p_born, s_o_Q2, QB, BI)
      class(virtual_t), intent(in) :: virtual
      integer, intent(in) :: i_flv
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in) :: s_o_Q2
      real(default), intent(inout), dimension(:) :: QB
      real(default), intent(inout), dimension(:) :: BI
    end subroutine virtual_compute_eikonals_threshold
    module subroutine virtual_set_bad_point (virt, value)
      class(virtual_t), intent(inout) :: virt
      logical, intent(in) :: value
    end subroutine virtual_set_bad_point
    module subroutine virtual_evaluate_initial_state &
         (virt, i_flv, i_corr, reg_data, QB)
      class(virtual_t), intent(inout) :: virt
      type(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: i_flv, i_corr
      real(default), intent(inout), dimension(:) :: QB
    end subroutine virtual_evaluate_initial_state
    module subroutine virtual_compute_collinear_contribution &
         (virt, i_flv, i_corr, p_born, sqrts, reg_data, QB)
      class(virtual_t), intent(inout) :: virt
      integer, intent(in) :: i_flv, i_corr
      type(vector4_t), dimension(:), intent(in) :: p_born
      real(default), intent(in) :: sqrts
      type(region_data_t), intent(in) :: reg_data
      real(default), intent(inout), dimension(:) :: QB
    end subroutine virtual_compute_collinear_contribution
    module subroutine virtual_compute_massive_self_eikonals (virt, &
         i_flv, i_corr, p_born, s_over_Q2, reg_data,  QB)
      class(virtual_t), intent(inout) :: virt
      integer, intent(in) :: i_flv, i_corr
      type(vector4_t), intent(in), dimension(:) :: p_born
      real(default), intent(in) :: s_over_Q2
      type(region_data_t), intent(in) :: reg_data
      real(default), intent(inout), dimension(:) :: QB
    end subroutine virtual_compute_massive_self_eikonals
    module subroutine virtual_final (virtual)
      class(virtual_t), intent(inout) :: virtual
    end subroutine virtual_final
  end interface

end module virtual
