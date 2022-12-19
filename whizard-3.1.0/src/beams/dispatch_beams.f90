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

module dispatch_beams

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use constants, only: PI, one
  use numeric_utils, only: vanishes
  use diagnostics
  use os_interface, only: os_data_t
  use variables, only: var_list_t
  use rng_base, only: rng_factory_t
  use pdg_arrays
  use model_data, only: model_data_t
  use flavors, only: flavor_t
  use physics_defs, only: PHOTON
  use physics_defs, only: MZ_REF, ME_REF, ALPHA_QCD_MZ_REF, ALPHA_QED_ME_REF
  use sm_qcd, only: qcd_t, alpha_qcd_fixed_t, alpha_qcd_from_scale_t
  use sm_qcd, only: alpha_qcd_from_lambda_t
  use sm_qed, only: qed_t, alpha_qed_fixed_t, alpha_qed_from_scale_t
  use beam_structures
  use dispatch_rng, only: dispatch_rng_factory
  use dispatch_rng, only: update_rng_seed_in_var_list
  use sf_base
  use sf_mappings
  use sf_isr
  use sf_epa
  use sf_ewa
  use sf_escan
  use sf_gaussian
  use sf_beam_events
  use sf_circe1
  use sf_circe2
  use sf_pdf_builtin
  use sf_lhapdf

  implicit none
  private

  public :: sf_prop_t
  public :: dispatch_sf_data
  public :: dispatch_sf_data_extra
  public :: strfun_mode
  public :: dispatch_sf_config
  public :: dispatch_qcd
  public :: dispatch_qed

  type :: sf_prop_t
     real(default), dimension(2) :: isr_eps = 1
  end type sf_prop_t


  procedure (dispatch_sf_data), pointer :: &
       dispatch_sf_data_extra => null ()

  interface
    module function strfun_mode (name) result (n)
      type(string_t), intent(in) :: name
      integer :: n
    end function strfun_mode
    module subroutine dispatch_sf_config (sf_config, sf_prop, beam_structure, &
           var_list, var_list_global, model, os_data, sqrts, pdg_prc)
      type(sf_config_t), dimension(:), allocatable, intent(out) :: sf_config
      type(sf_prop_t), intent(out) :: sf_prop
      type(beam_structure_t), intent(inout) :: beam_structure
      type(var_list_t), intent(in) :: var_list
      type(var_list_t), intent(inout) :: var_list_global
      class(model_data_t), target, intent(in) :: model
      type(os_data_t), intent(in) :: os_data
      real(default), intent(in) :: sqrts
      class(sf_data_t), allocatable :: sf_data
      type(beam_structure_t) :: beam_structure_tmp
      type(pdg_array_t), dimension(:,:), intent(in) :: pdg_prc
      type(string_t), dimension(:), allocatable :: prt_in
      type(pdg_array_t), dimension(:), allocatable :: pdg_in
    end subroutine dispatch_sf_config
  end interface

contains

  subroutine dispatch_sf_data (data, sf_method, i_beam, sf_prop, &
       var_list, var_list_global, model, &
       os_data, sqrts, pdg_in, pdg_prc, polarized)
    class(sf_data_t), allocatable, intent(inout) :: data
    type(string_t), intent(in) :: sf_method
    integer, dimension(:), intent(in) :: i_beam
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_in
    type(pdg_array_t), dimension(:,:), intent(in) :: pdg_prc
    type(sf_prop_t), intent(inout) :: sf_prop
    type(var_list_t), intent(in) :: var_list
    type(var_list_t), intent(inout) :: var_list_global
    integer :: next_rng_seed
    class(model_data_t), target, intent(in) :: model
    type(os_data_t), intent(in) :: os_data
    real(default), intent(in) :: sqrts
    logical, intent(in) :: polarized
    type(pdg_array_t), dimension(:), allocatable :: pdg_out
    real(default) :: isr_alpha, isr_q_max, isr_mass
    integer :: isr_order
    logical :: isr_recoil, isr_keep_energy
    real(default) :: epa_alpha, epa_x_min, epa_q_min, epa_q_max, epa_mass
    logical :: epa_recoil, epa_keep_energy
    integer :: epa_int_mode
    type(string_t) :: epa_mode
    real(default) :: ewa_x_min, ewa_pt_max, ewa_mass
    logical :: ewa_recoil, ewa_keep_energy
    type(pdg_array_t), dimension(:), allocatable :: pdg_prc1
    integer :: ewa_id
    type(string_t) :: pdf_name
    type(string_t) :: lhapdf_dir, lhapdf_file
    type(string_t), dimension(13) :: lhapdf_photon_sets
    integer :: lhapdf_member, lhapdf_photon_scheme
    logical :: hoppet_b_matching
    class(rng_factory_t), allocatable :: rng_factory
    logical :: circe1_photon1, circe1_photon2, circe1_generate, &
         circe1_with_radiation
    real(default) :: circe1_sqrts, circe1_eps
    integer :: circe1_version, circe1_chattiness, &
         circe1_revision
    character(6) :: circe1_accelerator
    logical :: circe2_polarized
    type(string_t) :: circe2_design, circe2_file
    real(default), dimension(2) :: gaussian_spread
    logical :: beam_events_warn_eof
    type(string_t) :: beam_events_dir, beam_events_file
    logical :: escan_normalize
    integer :: i
    lhapdf_photon_sets = [var_str ("DOG0.LHgrid"), var_str ("DOG1.LHgrid"), &
         var_str ("DGG.LHgrid"), var_str ("LACG.LHgrid"), &
         var_str ("GSG0.LHgrid"), var_str ("GSG1.LHgrid"), &
         var_str ("GSG960.LHgrid"), var_str ("GSG961.LHgrid"), &
         var_str ("GRVG0.LHgrid"), var_str ("GRVG1.LHgrid"), &
         var_str ("ACFGPG.LHgrid"), var_str ("WHITG.LHgrid"), &
         var_str ("SASG.LHgrid")]
    select case (char (sf_method))
    case ("pdf_builtin")
       allocate (pdf_builtin_data_t :: data)
       select type (data)
       type is (pdf_builtin_data_t)
          pdf_name = &
               var_list%get_sval (var_str ("$pdf_builtin_set"))
          hoppet_b_matching = &
               var_list%get_lval (var_str ("?hoppet_b_matching"))
          call data%init ( &
               model, pdg_in(i_beam(1)), &
               name = pdf_name, &
               path = os_data%pdf_builtin_datapath, &
               hoppet_b_matching = hoppet_b_matching)
       end select
    case ("pdf_builtin_photon")
       call msg_fatal ("Currently, there are no photon PDFs built into WHIZARD,", &
            [var_str ("for the photon content inside a proton or neutron use"), &
             var_str ("the 'lhapdf_photon' structure function.")])
    case ("lhapdf")
       allocate (lhapdf_data_t :: data)
       if (pdg_in(i_beam(1))%get (1) == PHOTON) then
          call msg_fatal ("The 'lhapdf' structure is intended only for protons and", &
               [var_str ("pions, please use 'lhapdf_photon' for photon beams.")])
       end if
       lhapdf_dir = &
            var_list%get_sval (var_str ("$lhapdf_dir"))
       lhapdf_file = &
            var_list%get_sval (var_str ("$lhapdf_file"))
       lhapdf_member = &
            var_list%get_ival (var_str ("lhapdf_member"))
       lhapdf_photon_scheme = &
            var_list%get_ival (var_str ("lhapdf_photon_scheme"))
       hoppet_b_matching = &
            var_list%get_lval (var_str ("?hoppet_b_matching"))
       select type (data)
       type is (lhapdf_data_t)
          call data%init &
               (model, pdg_in(i_beam(1)), &
                lhapdf_dir, lhapdf_file, lhapdf_member, &
                lhapdf_photon_scheme, hoppet_b_matching)
       end select
    case ("lhapdf_photon")
       allocate (lhapdf_data_t :: data)
       if (pdg_in(i_beam(1))%get_length () /= 1 .or. &
            pdg_in(i_beam(1))%get (1) /= PHOTON) then
          call msg_fatal ("The 'lhapdf_photon' structure function is exclusively for", &
               [var_str ("photon PDFs, i.e. for photons as beam particles")])
       end if
       lhapdf_dir = &
            var_list%get_sval (var_str ("$lhapdf_dir"))
       lhapdf_file = &
            var_list%get_sval (var_str ("$lhapdf_photon_file"))
       lhapdf_member = &
            var_list%get_ival (var_str ("lhapdf_member"))
       lhapdf_photon_scheme = &
            var_list%get_ival (var_str ("lhapdf_photon_scheme"))
       if (.not. any (lhapdf_photon_sets == lhapdf_file)) then
          call msg_fatal ("This PDF set is not supported or not " // &
               "intended for photon beams.")
       end if
       select type (data)
       type is (lhapdf_data_t)
          call data%init &
               (model, pdg_in(i_beam(1)), &
                lhapdf_dir, lhapdf_file, lhapdf_member, &
                lhapdf_photon_scheme)
       end select
    case ("isr")
       allocate (isr_data_t :: data)
       isr_alpha = &
            var_list%get_rval (var_str ("isr_alpha"))
       if (vanishes (isr_alpha)) then
          isr_alpha = (var_list%get_rval (var_str ("ee"))) &
               ** 2 / (4 * PI)
       end if
       isr_q_max = &
            var_list%get_rval (var_str ("isr_q_max"))
       if (vanishes (isr_q_max)) then
          isr_q_max = sqrts
       end if
       isr_mass   = var_list%get_rval (var_str ("isr_mass"))
       isr_order  = var_list%get_ival (var_str ("isr_order"))
       isr_recoil = var_list%get_lval (var_str ("?isr_recoil"))
       isr_keep_energy = var_list%get_lval (var_str ("?isr_keep_energy"))
       select type (data)
       type is (isr_data_t)
          call data%init &
               (model, pdg_in (i_beam(1)), isr_alpha, isr_q_max, &
               isr_mass, isr_order, recoil = isr_recoil, keep_energy = &
               isr_keep_energy)
          call data%check ()
          sf_prop%isr_eps(i_beam(1)) = data%get_eps ()
       end select
    case ("epa")
       allocate (epa_data_t :: data)
       epa_mode = var_list%get_sval (var_str ("$epa_mode"))
       epa_int_mode = 0
       epa_alpha = var_list%get_rval (var_str ("epa_alpha"))
       if (vanishes (epa_alpha)) then
          epa_alpha = (var_list%get_rval (var_str ("ee"))) &
               ** 2 / (4 * PI)
       end if
       epa_x_min = var_list%get_rval (var_str ("epa_x_min"))
       epa_q_min = var_list%get_rval (var_str ("epa_q_min"))
       epa_q_max = var_list%get_rval (var_str ("epa_q_max"))
       if (vanishes (epa_q_max)) then
          epa_q_max = sqrts
       end if
       select case (char (epa_mode))
       case ("default", "Budnev_617")
          epa_int_mode = 0
       case ("Budnev_616e")
          epa_int_mode = 1
       case ("log_power")
          epa_int_mode = 2
          epa_q_max = sqrts
       case ("log_simple")
          epa_int_mode = 3
          epa_q_max = sqrts
       case ("log")
          epa_int_mode = 4
          epa_q_max = sqrts
       case default
          call msg_fatal ("EPA: unsupported EPA mode; please choose " // &
               "'default', 'Budnev_616', 'Budnev_616e', 'log_power', " // &
               "'log_simple', or 'log'")
       end select
       epa_mass   = var_list%get_rval (var_str ("epa_mass"))
       epa_recoil = var_list%get_lval (var_str ("?epa_recoil"))
       epa_keep_energy = var_list%get_lval (var_str ("?epa_keep_energy"))
       select type (data)
       type is (epa_data_t)
          call data%init &
               (model, epa_int_mode, pdg_in (i_beam(1)), epa_alpha, &
               epa_x_min, epa_q_min, epa_q_max, epa_mass, &
               recoil = epa_recoil, keep_energy = epa_keep_energy)
          call data%check ()
       end select
    case ("ewa")
       allocate (ewa_data_t :: data)
       allocate (pdg_prc1 (size (pdg_prc, 2)))
       pdg_prc1 = pdg_prc(i_beam(1),:)
       if (any (pdg_prc1%get_length () /= 1) &
            .or. any (pdg_prc1 /= pdg_prc1(1))) then
          call msg_fatal &
               ("EWA: process incoming particle (W/Z) must be unique")
       end if
       ewa_id = abs (pdg_prc1(1)%get (1))
       ewa_x_min = var_list%get_rval (var_str ("ewa_x_min"))
       ewa_pt_max = var_list%get_rval (var_str ("ewa_pt_max"))
       if (vanishes (ewa_pt_max)) then
          ewa_pt_max = sqrts
       end if
       ewa_mass = var_list%get_rval (var_str ("ewa_mass"))
       ewa_recoil = var_list%get_lval (&
            var_str ("?ewa_recoil"))
       ewa_keep_energy = var_list%get_lval (&
            var_str ("?ewa_keep_energy"))
       select type (data)
       type is (ewa_data_t)
          call data%init &
               (model, pdg_in (i_beam(1)), ewa_x_min, &
               ewa_pt_max, sqrts, ewa_recoil, &
               ewa_keep_energy, ewa_mass)
          call data%set_id (ewa_id)
          call data%check ()
       end select
    case ("circe1")
       allocate (circe1_data_t :: data)
       select type (data)
       type is (circe1_data_t)
          circe1_photon1 = &
               var_list%get_lval (var_str ("?circe1_photon1"))
          circe1_photon2 = &
               var_list%get_lval (var_str ("?circe1_photon2"))
          circe1_sqrts = &
               var_list%get_rval (var_str ("circe1_sqrts"))
          circe1_eps = &
               var_list%get_rval (var_str ("circe1_eps"))
          if (circe1_sqrts <= 0)  circe1_sqrts = sqrts
          circe1_generate = &
               var_list%get_lval (var_str ("?circe1_generate"))
          circe1_version = &
               var_list%get_ival (var_str ("circe1_ver"))
          circe1_revision = &
               var_list%get_ival (var_str ("circe1_rev"))
          circe1_accelerator = &
               char (var_list%get_sval (var_str ("$circe1_acc")))
          circe1_chattiness = &
               var_list%get_ival (var_str ("circe1_chat"))
          circe1_with_radiation = &
               var_list%get_lval (var_str ("?circe1_with_radiation"))
          call data%init (model, pdg_in, circe1_sqrts, circe1_eps, &
               [circe1_photon1, circe1_photon2], &
               circe1_version, circe1_revision, circe1_accelerator, &
               circe1_chattiness, circe1_with_radiation)
          if (circe1_generate) then
             call msg_message ("CIRCE1: activating generator mode")
             call dispatch_rng_factory &
                  (rng_factory, var_list_global, next_rng_seed)
             call update_rng_seed_in_var_list (var_list_global, next_rng_seed)
             call data%set_generator_mode (rng_factory)
          end if
       end select
    case ("circe2")
       allocate (circe2_data_t :: data)
       select type (data)
       type is (circe2_data_t)
          circe2_polarized = &
               var_list%get_lval (var_str ("?circe2_polarized"))
          circe2_file = &
               var_list%get_sval (var_str ("$circe2_file"))
          circe2_design = &
               var_list%get_sval (var_str ("$circe2_design"))
          call data%init (os_data, model, pdg_in, sqrts, &
               circe2_polarized, polarized, circe2_file, circe2_design)
          call msg_message ("CIRCE2: activating generator mode")
          call dispatch_rng_factory &
               (rng_factory, var_list_global, next_rng_seed)
          call update_rng_seed_in_var_list (var_list_global, next_rng_seed)
          call data%set_generator_mode (rng_factory)
       end select
    case ("gaussian")
       allocate (gaussian_data_t :: data)
       select type (data)
       type is (gaussian_data_t)
          gaussian_spread = &
               [var_list%get_rval (var_str ("gaussian_spread1")), &
               var_list%get_rval (var_str ("gaussian_spread2"))]
          call dispatch_rng_factory &
               (rng_factory, var_list_global, next_rng_seed)
          call update_rng_seed_in_var_list (var_list_global, next_rng_seed)
          call data%init (model, pdg_in, gaussian_spread, rng_factory)
       end select
    case ("beam_events")
       allocate (beam_events_data_t :: data)
       select type (data)
       type is (beam_events_data_t)
          beam_events_dir = os_data%whizard_beamsimpath
          beam_events_file = var_list%get_sval (&
               var_str ("$beam_events_file"))
          beam_events_warn_eof = var_list%get_lval (&
               var_str ("?beam_events_warn_eof"))
          call data%init (model, pdg_in, &
                  beam_events_dir, beam_events_file, beam_events_warn_eof)
       end select
    case ("energy_scan")
       escan_normalize = &
            var_list%get_lval (var_str ("?energy_scan_normalize"))
       allocate (escan_data_t :: data)
       select type (data)
       type is (escan_data_t)
          if (escan_normalize) then
             call data%init (model, pdg_in)
          else
             call data%init (model, pdg_in, sqrts)
          end if
       end select
    case default
       if (associated (dispatch_sf_data_extra)) then
          call dispatch_sf_data_extra (data, sf_method, i_beam, &
               sf_prop, var_list, var_list_global, model, os_data, sqrts, pdg_in, &
               pdg_prc, polarized)
       end if
       if (.not. allocated (data)) then
          call msg_fatal ("Structure function '" &
               // char (sf_method) // "' not implemented")
       end if
    end select
    if (allocated (data)) then
       allocate (pdg_out (size (pdg_prc, 1)))
       call data%get_pdg_out (pdg_out)
       do i = 1, size (i_beam)
          pdg_in(i_beam(i)) = pdg_out(i)
       end do
    end if
  end subroutine dispatch_sf_data

  subroutine dispatch_qcd (qcd, var_list, os_data)
    type(qcd_t), intent(inout) :: qcd
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    logical :: fixed, from_mz, from_pdf_builtin, from_lhapdf, from_lambda_qcd
    real(default) :: mz, alpha_val, lambda
    integer :: nf, order, lhapdf_member
    type(string_t) :: pdfset, lhapdf_dir, lhapdf_file
    call unpack_variables ()
    if (allocated (qcd%alpha))  deallocate (qcd%alpha)
    if (from_lhapdf .and. from_pdf_builtin) then
        call msg_fatal (" Mixing alphas evolution",  &
             [var_str (" from LHAPDF and builtin PDF is not permitted")])
    end if
    select case (count ([from_mz, from_pdf_builtin, from_lhapdf, from_lambda_qcd]))
    case (0)
       if (fixed) then
          allocate (alpha_qcd_fixed_t :: qcd%alpha)
       else
          call msg_fatal ("QCD alpha: no calculation mode set")
       end if
    case (2:)
       call msg_fatal ("QCD alpha: calculation mode is ambiguous")
    case (1)
       if (fixed) then
          call msg_fatal ("QCD alpha: use '?alphas_is_fixed = false' for " // &
               "running alphas")
       else if (from_mz) then
          allocate (alpha_qcd_from_scale_t :: qcd%alpha)
       else if (from_pdf_builtin) then
          allocate (alpha_qcd_pdf_builtin_t :: qcd%alpha)
       else if (from_lhapdf) then
          allocate (alpha_qcd_lhapdf_t :: qcd%alpha)
       else if (from_lambda_qcd) then
          allocate (alpha_qcd_from_lambda_t :: qcd%alpha)
       end if
       call msg_message ("QCD alpha: using a running strong coupling")
    end select
    call init_alpha ()
    qcd%n_f = var_list%get_ival (var_str ("alphas_nf"))
  contains
    subroutine unpack_variables ()
      fixed = var_list%get_lval (var_str ("?alphas_is_fixed"))
      from_mz = var_list%get_lval (var_str ("?alphas_from_mz"))
      from_pdf_builtin = &
           var_list%get_lval (var_str ("?alphas_from_pdf_builtin"))
      from_lhapdf = &
           var_list%get_lval (var_str ("?alphas_from_lhapdf"))
      from_lambda_qcd = &
           var_list%get_lval (var_str ("?alphas_from_lambda_qcd"))
      pdfset = var_list%get_sval (var_str ("$pdf_builtin_set"))
      lambda = var_list%get_rval (var_str ("lambda_qcd"))
      nf = var_list%get_ival (var_str ("alphas_nf"))
      order = var_list%get_ival (var_str ("alphas_order"))
      lhapdf_dir = var_list%get_sval (var_str ("$lhapdf_dir"))
      lhapdf_file = var_list%get_sval (var_str ("$lhapdf_file"))
      lhapdf_member = var_list%get_ival (var_str ("lhapdf_member"))
      if (var_list%contains (var_str ("mZ"))) then
         mz = var_list%get_rval (var_str ("mZ"))
      else
         mz = MZ_REF
      end if
      if (var_list%contains (var_str ("alphas"))) then
         alpha_val = var_list%get_rval (var_str ("alphas"))
      else
         alpha_val = ALPHA_QCD_MZ_REF
      end if
    end subroutine unpack_variables

    subroutine init_alpha ()
      select type (alpha => qcd%alpha)
      type is (alpha_qcd_fixed_t)
         alpha%val = alpha_val
      type is (alpha_qcd_from_scale_t)
         alpha%mu_ref = mz
         alpha%ref = alpha_val
         alpha%order = order
         alpha%nf = nf
      type is (alpha_qcd_from_lambda_t)
         alpha%lambda = lambda
         alpha%order = order
         alpha%nf = nf
      type is (alpha_qcd_pdf_builtin_t)
         call alpha%init (pdfset, &
              os_data%pdf_builtin_datapath)
      type is (alpha_qcd_lhapdf_t)
         call alpha%init (lhapdf_file, lhapdf_member, lhapdf_dir)
      end select
    end subroutine init_alpha

  end subroutine dispatch_qcd

  subroutine dispatch_qed (qed, var_list)
    type(qed_t), intent(inout) :: qed
    type(var_list_t), intent(in) :: var_list
    logical :: fixed, from_me, analytic
    real(default) :: me, alpha_val
    integer :: nf, nlep, order
    call unpack_variables ()
    if (allocated (qed%alpha))  deallocate (qed%alpha)
    select case (count ([from_me]))
    case (0)
       if (fixed) then
          allocate (alpha_qed_fixed_t :: qed%alpha)
       else
          call msg_fatal ("QED alpha: no calculation mode set")
       end if
    case (2:)
       call msg_fatal ("QED alpha: calculation mode is ambiguous")
    case (1)
       if (fixed) then
          call msg_fatal ("QED alpha: use '?alphas_is_fixed = false' for " // &
               "running alpha")
       else if (from_me) then
          allocate (alpha_qed_from_scale_t :: qed%alpha)
       end if
       call msg_message ("QED alpha: using a running electromagnetic coupling")
    end select
    call init_alpha ()
    if (var_list%get_ival (var_str ("alpha_nf")) == -1) then
       qed%n_f = var_list%get_ival (var_str ("alphas_nf"))
    else
       qed%n_f = var_list%get_ival (var_str ("alpha_nf"))
    end if
    qed%n_lep = var_list%get_ival (var_str ("alpha_nlep"))
  contains
    subroutine unpack_variables ()
      fixed = var_list%get_lval (var_str ("?alpha_is_fixed"))
      from_me = var_list%get_lval (var_str ("?alpha_from_me"))
      if (var_list%get_ival (var_str ("alpha_nf")) == -1) then
         nf = var_list%get_ival (var_str ("alphas_nf"))
      else
         nf = var_list%get_ival (var_str ("alpha_nf"))
      end if
      analytic = var_list%get_lval (var_str ("?alpha_evolve_analytic"))
      nlep = var_list%get_ival (var_str ("alpha_nlep"))
      order = var_list%get_ival (var_str ("alpha_order"))
      if (var_list%contains (var_str ("me"))) then
         me = var_list%get_rval (var_str ("me"))
      else
         me = ME_REF
      end if
      if (var_list%contains (var_str ("alpha_em_i"))) then
         alpha_val = one / var_list%get_rval (var_str ("alpha_em_i"))
      else
         alpha_val = ALPHA_QED_ME_REF
      end if
    end subroutine unpack_variables

    subroutine init_alpha ()
      select type (alpha => qed%alpha)
      type is (alpha_qed_fixed_t)
         alpha%val = alpha_val
      type is (alpha_qed_from_scale_t)
         alpha%mu_ref = me
         alpha%ref = alpha_val
         alpha%order = order
         alpha%nf = nf
         alpha%nlep = nlep
         alpha%analytic = analytic
      end select
    end subroutine init_alpha

  end subroutine dispatch_qed


end module dispatch_beams
