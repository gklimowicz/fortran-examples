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

module kinematics

  use kinds, only: default
  use lorentz
  use physics_defs
  use sf_base
  use phs_base
  use fks_regions
  use mci_base
  use process_config
  use process_mci
  use pcm_base, only: pcm_t, pcm_workspace_t
  use pcm, only: pcm_nlo_t, pcm_nlo_workspace_t

  implicit none
  private

  public :: kinematics_t

  type :: kinematics_t
     integer :: n_in = 0
     integer :: n_channel = 0
     integer :: selected_channel = 0
     type(sf_chain_instance_t), pointer :: sf_chain => null ()
     class(phs_t), pointer :: phs => null ()
     real(default), dimension(:), pointer :: f => null ()
     real(default) :: phs_factor
     logical :: sf_chain_allocated = .false.
     logical :: phs_allocated = .false.
     logical :: f_allocated = .false.
     integer :: emitter = -1
     integer :: i_phs = 0
     integer :: i_con = 0
     logical :: only_cm_frame = .false.
     logical :: new_seed = .true.
     logical :: threshold = .false.
   contains
     procedure :: write => kinematics_write
     procedure :: final => kinematics_final
     procedure :: configure => kinematics_configure
     procedure :: set_nlo_info => kinematics_set_nlo_info
     procedure :: set_threshold => kinematics_set_threshold
     procedure :: init_sf_chain => kinematics_init_sf_chain
     procedure :: init_phs => kinematics_init_phs
     procedure :: evaluate_radiation_kinematics => &
          kinematics_evaluate_radiation_kinematics
     procedure :: generate_fsr_in => kinematics_generate_fsr_in
     procedure :: compute_xi_ref_momenta => kinematics_compute_xi_ref_momenta
     procedure :: compute_selected_channel => kinematics_compute_selected_channel
     procedure :: redo_sf_chain => kinematics_redo_sf_chain
     procedure :: compute_other_channels => kinematics_compute_other_channels
     procedure :: get_incoming_momenta => kinematics_get_incoming_momenta
     procedure :: get_boost_to_lab => kinematics_get_boost_to_lab
     procedure :: get_boost_to_cms => kinematics_get_boost_to_cms
     procedure :: recover_mcpar => kinematics_recover_mcpar
     procedure :: recover_sfchain => kinematics_recover_sfchain
     procedure :: get_mcpar => kinematics_get_mcpar
     procedure :: evaluate_sf_chain => kinematics_evaluate_sf_chain
     procedure :: return_beam_momenta => kinematics_return_beam_momenta
     procedure :: lab_is_cm => kinematics_lab_is_cm
     procedure :: modify_momenta_for_subtraction => &
          kinematics_modify_momenta_for_subtraction
     procedure :: threshold_projection => kinematics_threshold_projection
     procedure :: evaluate_radiation => kinematics_evaluate_radiation
  end type kinematics_t


  interface
    module subroutine kinematics_write (object, unit)
      class(kinematics_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine kinematics_write
    module subroutine kinematics_final (object)
      class(kinematics_t), intent(inout) :: object
    end subroutine kinematics_final
    module subroutine kinematics_configure (kin, pcm, pcm_work, &
         sf_chain, beam_config, phs_config, nlo_type, is_i_sub)
      class(kinematics_t), intent(out) :: kin
      class(pcm_t), intent(inout) :: pcm
      class(pcm_workspace_t), intent(in) :: pcm_work
      type(sf_chain_t), intent(in), target :: sf_chain
      type(process_beam_config_t), intent(in), target :: beam_config
      class(phs_config_t), intent(in), target :: phs_config
      integer, intent(in) :: nlo_type
      logical, intent(in) :: is_i_sub
    end subroutine kinematics_configure
    module subroutine kinematics_set_nlo_info (k, nlo_type)
      class(kinematics_t), intent(inout) :: k
      integer, intent(in) :: nlo_type
    end subroutine kinematics_set_nlo_info
    module subroutine kinematics_set_threshold (kin, factorization_mode)
      class(kinematics_t), intent(inout) :: kin
      integer, intent(in) :: factorization_mode
    end subroutine kinematics_set_threshold
    module subroutine kinematics_init_sf_chain &
         (k, sf_chain, config, extended_sf)
      class(kinematics_t), intent(inout) :: k
      type(sf_chain_t), intent(in), target :: sf_chain
      type(process_beam_config_t), intent(in) :: config
      logical, intent(in), optional :: extended_sf
    end subroutine kinematics_init_sf_chain
    module subroutine kinematics_init_phs (k, config)
      class(kinematics_t), intent(inout) :: k
      class(phs_config_t), intent(in), target :: config
    end subroutine kinematics_init_phs
    module subroutine kinematics_evaluate_radiation_kinematics (k, r_in)
      class(kinematics_t), intent(inout) :: k
      real(default), intent(in), dimension(:) :: r_in
    end subroutine kinematics_evaluate_radiation_kinematics
    module subroutine kinematics_generate_fsr_in (kin)
      class(kinematics_t), intent(inout) :: kin
    end subroutine kinematics_generate_fsr_in
    module subroutine kinematics_compute_xi_ref_momenta (k, reg_data, nlo_type)
      class(kinematics_t), intent(inout) :: k
      type(region_data_t), intent(in) :: reg_data
      integer, intent(in) :: nlo_type
    end subroutine kinematics_compute_xi_ref_momenta
    module subroutine kinematics_compute_selected_channel &
         (k, mci_work, phs_channel, p, success)
      class(kinematics_t), intent(inout) :: k
      type(mci_work_t), intent(in) :: mci_work
      integer, intent(in) :: phs_channel
      type(vector4_t), dimension(:), intent(out) :: p
      logical, intent(out) :: success
    end subroutine kinematics_compute_selected_channel
    module subroutine kinematics_redo_sf_chain (kin, mci_work, phs_channel)
      class(kinematics_t), intent(inout) :: kin
      type(mci_work_t), intent(in) :: mci_work
      integer, intent(in) :: phs_channel
    end subroutine kinematics_redo_sf_chain
    module subroutine kinematics_compute_other_channels &
         (k, mci_work, phs_channel)
      class(kinematics_t), intent(inout) :: k
      type(mci_work_t), intent(in) :: mci_work
      integer, intent(in) :: phs_channel
    end subroutine kinematics_compute_other_channels
    module subroutine kinematics_get_incoming_momenta (k, p)
      class(kinematics_t), intent(in) :: k
      type(vector4_t), dimension(:), intent(out) :: p
    end subroutine kinematics_get_incoming_momenta
    module function kinematics_get_boost_to_lab (kin) result (lt)
      type(lorentz_transformation_t) :: lt
      class(kinematics_t), intent(in) :: kin
    end function kinematics_get_boost_to_lab
    module function kinematics_get_boost_to_cms (kin) result (lt)
      type(lorentz_transformation_t) :: lt
      class(kinematics_t), intent(in) :: kin
    end function kinematics_get_boost_to_cms
    module subroutine kinematics_recover_mcpar (k, mci_work, phs_channel, p)
      class(kinematics_t), intent(inout) :: k
      type(mci_work_t), intent(inout) :: mci_work
      integer, intent(in) :: phs_channel
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine kinematics_recover_mcpar
    module subroutine kinematics_recover_sfchain (k, channel, p)
      class(kinematics_t), intent(inout) :: k
      integer, intent(in) :: channel
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine kinematics_recover_sfchain
    module subroutine kinematics_get_mcpar (k, phs_channel, r)
      class(kinematics_t), intent(in) :: k
      integer, intent(in) :: phs_channel
      real(default), dimension(:), intent(out) :: r
    end subroutine kinematics_get_mcpar
    module subroutine kinematics_evaluate_sf_chain &
         (k, fac_scale, negative_sf, sf_rescale)
      class(kinematics_t), intent(inout) :: k
      real(default), intent(in) :: fac_scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(inout), optional :: sf_rescale
    end subroutine kinematics_evaluate_sf_chain
    module subroutine kinematics_return_beam_momenta (k)
      class(kinematics_t), intent(in) :: k
    end subroutine kinematics_return_beam_momenta
    module function kinematics_lab_is_cm (k) result (lab_is_cm)
       logical :: lab_is_cm
       class(kinematics_t), intent(in) :: k
    end function kinematics_lab_is_cm
    module subroutine kinematics_modify_momenta_for_subtraction (k, p_in, p_out)
      class(kinematics_t), intent(inout) :: k
      type(vector4_t), intent(in), dimension(:) :: p_in
      type(vector4_t), intent(out), dimension(:), allocatable :: p_out
    end subroutine kinematics_modify_momenta_for_subtraction
    module subroutine kinematics_threshold_projection (k, pcm_work, nlo_type)
      class(kinematics_t), intent(inout) :: k
      type(pcm_nlo_workspace_t), intent(inout) :: pcm_work
      integer, intent(in) :: nlo_type
    end subroutine kinematics_threshold_projection
    module subroutine kinematics_evaluate_radiation (k, p_in, p_out, success)
      class(kinematics_t), intent(inout) :: k
      type(vector4_t), intent(in), dimension(:) :: p_in
      type(vector4_t), intent(out), dimension(:), allocatable :: p_out
      logical, intent(out) :: success
    end subroutine kinematics_evaluate_radiation
  end interface

end module kinematics
