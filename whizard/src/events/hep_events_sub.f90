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

submodule (hep_events) hep_events_s

  use system_dependencies, only: HEPMC2_AVAILABLE
  use system_dependencies, only: HEPMC3_AVAILABLE
  use diagnostics
  use numeric_utils
  use flavors
  use colors
  use helicities
  use subevents, only: PRT_BEAM, PRT_INCOMING, PRT_OUTGOING
  use subevents, only: PRT_UNDEFINED
  use subevents, only: PRT_VIRTUAL, PRT_RESONANT, PRT_BEAM_REMNANT
  use hep_common

  implicit none

contains

  module subroutine hepeup_from_event &
       (event, keep_beams, keep_remnants, process_index)
    class(generic_event_t), intent(in), target :: event
    logical, intent(in), optional :: keep_beams
    logical, intent(in), optional :: keep_remnants
    integer, intent(in), optional :: process_index
    type(particle_set_t), pointer :: particle_set
    real(default) :: scale, alpha_qcd
    if (event%has_valid_particle_set ()) then
       particle_set => event%get_particle_set_ptr ()
       call hepeup_from_particle_set (particle_set, keep_beams, keep_remnants)
       if (present (process_index)) then
          call hepeup_set_event_parameters (proc_id = process_index)
       end if
       scale = event%get_fac_scale ()
       if (.not. vanishes (scale)) then
          call hepeup_set_event_parameters (scale = scale)
       end if
       alpha_qcd = event%get_alpha_s ()
       if (.not. vanishes (alpha_qcd)) then
          call hepeup_set_event_parameters (alpha_qcd = alpha_qcd)
       end if
       if (event%weight_prc_is_known ()) then
          call hepeup_set_event_parameters (weight = event%get_weight_prc ())
       end if
    else
       call msg_bug ("HEPEUP: event incomplete")
    end if
  end subroutine hepeup_from_event

  module subroutine hepeup_to_event &
       (event, fallback_model, process_index, recover_beams, &
       use_alpha_s, use_scale)
    class(generic_event_t), intent(inout), target :: event
    class(model_data_t), intent(in), target :: fallback_model
    integer, intent(out), optional :: process_index
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: use_alpha_s
    logical, intent(in), optional :: use_scale
    class(model_data_t), pointer :: model
    real(default) :: weight, scale, alpha_qcd
    type(particle_set_t) :: particle_set
    model => event%get_model_ptr ()
    call hepeup_to_particle_set &
         (particle_set, recover_beams, model, fallback_model)
    call event%set_hard_particle_set (particle_set)
    call particle_set%final ()
    if (present (process_index)) then
       call hepeup_get_event_parameters (proc_id = process_index)
    end if
    call hepeup_get_event_parameters (weight = weight, &
         scale = scale, alpha_qcd = alpha_qcd)
    call event%set_weight_ref (weight)
    if (present (use_alpha_s)) then
       if (use_alpha_s .and. alpha_qcd > 0) &
            call event%set_alpha_qcd_forced (alpha_qcd)
    end if
    if (present (use_scale)) then
       if (use_scale .and. scale > 0) &
            call event%set_scale_forced (scale)
    end if
  end subroutine hepeup_to_event

  module subroutine hepevt_from_event  &
         (event, process_index, i_evt, keep_beams, keep_remnants, &
          ensure_order, fill_hepev4)
    class(generic_event_t), intent(in), target :: event
    integer, intent(in), optional :: i_evt, process_index
    logical, intent(in), optional :: keep_beams
    logical, intent(in), optional :: keep_remnants
    logical, intent(in), optional :: ensure_order
    logical, intent(in), optional :: fill_hepev4
    type(particle_set_t), pointer :: particle_set
    real(default) :: alpha_qcd, scale
    if (event%has_valid_particle_set ()) then
       particle_set => event%get_particle_set_ptr ()
       call hepevt_from_particle_set (particle_set, keep_beams, &
            keep_remnants, ensure_order, fill_hepev4)
       if (present (process_index)) then
          call hepevt_set_event_parameters (proc_id = process_index)
       end if
       if (event%weight_prc_is_known ()) then
          call hepevt_set_event_parameters (weight = event%get_weight_prc ())
       end if
       if (event%sqme_prc_is_known ()) then
          call hepevt_set_event_parameters &
               (function_value = event%get_sqme_prc ())
       end if
       scale = event%get_fac_scale ()
       if (.not. vanishes (scale)) then
          call hepevt_set_event_parameters (scale = scale)
       end if
       alpha_qcd = event%get_alpha_s ()
       if (.not. vanishes (alpha_qcd)) then
          call hepevt_set_event_parameters (alpha_qcd = alpha_qcd)
       end if
       if (present (i_evt)) then
          call hepevt_set_event_parameters (i_evt = i_evt)
       else if (event%has_index ()) then
          call hepevt_set_event_parameters (i_evt = event%get_index ())
       else
          call hepevt_set_event_parameters (i_evt = 0)
       end if
    else
       call msg_bug ("HEPEVT: event incomplete")
    end if
  end subroutine hepevt_from_event

  subroutine particle_to_hepmc (prt, hprt)
    type(particle_t), intent(in) :: prt
    type(hepmc_particle_t), intent(out) :: hprt
    integer :: hepmc_status
    select case (prt%get_status ())
    case (PRT_UNDEFINED)
       hepmc_status = 0
    case (PRT_OUTGOING)
       hepmc_status = 1
    case (PRT_BEAM)
       hepmc_status = 4
    case (PRT_RESONANT)
       hepmc_status = 2
    case (PRT_BEAM_REMNANT)
       if (prt%get_n_children () == 0) then
          hepmc_status = 1
       else
          hepmc_status = 3
       end if
    case default
       hepmc_status = 3
    end select
    call hepmc_particle_init (hprt, &
         prt%get_momentum (), prt%get_pdg (), &
         hepmc_status)
    if (HEPMC2_AVAILABLE) then
       call hepmc_particle_set_color (hprt, prt%get_color ())
       select case (prt%get_polarization_status ())
       case (PRT_DEFINITE_HELICITY)
          call hepmc_particle_set_polarization (hprt, &
               prt%get_helicity ())
       case (PRT_GENERIC_POLARIZATION)
          call hepmc_particle_set_polarization (hprt, &
               prt%get_polarization ())
       end select
    end if
  end subroutine particle_to_hepmc

  module subroutine hepmc_event_from_particle_set &
         (evt, particle_set, cross_section, error, color)
    type(hepmc_event_t), intent(inout) :: evt
    type(particle_set_t), intent(in) :: particle_set
    real(default), intent(in), optional :: cross_section, error
    logical, intent(in), optional :: color
    type(hepmc_vertex_t), dimension(:), allocatable :: v
    type(hepmc_particle_t), dimension(:), allocatable :: hprt
    type(hepmc_particle_t), dimension(2) :: hbeam
    type(vector4_t), dimension(:), allocatable :: vtx
    logical, dimension(:), allocatable :: is_beam
    integer, dimension(:), allocatable :: v_from, v_to
    integer :: n_vertices, n_tot, i
    logical :: write_color
    write_color = .false.
    if (present (color))  write_color = color
    n_tot = particle_set%get_n_tot ()
    allocate (v_from (n_tot), v_to (n_tot))
    call particle_set%assign_vertices (v_from, v_to, n_vertices)
    allocate (hprt (n_tot))
    allocate (vtx (n_vertices))
    vtx = vector4_null
    do i = 1, n_tot
       if (v_to(i) /= 0 .or. v_from(i) /= 0) then
          call particle_to_hepmc (particle_set%prt(i), hprt(i))
          if (v_from(i) /= 0) then
             vtx(v_from(i)) = particle_set%prt(i)%get_vertex ()
          end if
       end if
    end do
    if (present (cross_section) .and. present(error)) &
       call hepmc_event_set_cross_section (evt, cross_section, error)
    allocate (v (n_vertices))
    do i = 1, n_vertices
       call hepmc_vertex_init (v(i), vtx(i))
       call hepmc_event_add_vertex (evt, v(i))
    end do
    allocate (is_beam (n_tot))
    is_beam = particle_set%prt(1:n_tot)%get_status () == PRT_BEAM
    if (.not. any (is_beam)) then
       is_beam = particle_set%prt(1:n_tot)%get_status () == PRT_INCOMING
    end if
    if (count (is_beam) == 2) then
       hbeam = pack (hprt, is_beam)
       call hepmc_event_set_beam_particles (evt, hbeam(1), hbeam(2))
    end if
    do i = 1, n_tot
       if (v_to(i) /= 0) then
          call hepmc_vertex_add_particle_in (v(v_to(i)), hprt(i))
       end if
    end do
    do i = 1, n_tot
       if (v_from(i) /= 0) then
          call hepmc_vertex_add_particle_out (v(v_from(i)), hprt(i))
       end if
    end do
    FIND_SIGNAL_PROCESS: do i = 1, n_tot
       if (particle_set%prt(i)%get_status () == PRT_INCOMING) then
          call hepmc_event_set_signal_process_vertex (evt, v(v_to(i)))
          exit FIND_SIGNAL_PROCESS
       end if
    end do FIND_SIGNAL_PROCESS
    if (HEPMC3_AVAILABLE) then
       do i = 1, n_tot
          if (write_color) then
             call hepmc_particle_set_color (hprt(i), &
                  particle_set%prt(i)%get_color ())
          end if
          select case (particle_set%prt(i)%get_polarization_status ())
          case (PRT_DEFINITE_HELICITY)
             call hepmc_particle_set_polarization (hprt(i), &
                  particle_set%prt(i)%get_helicity ())
          case (PRT_GENERIC_POLARIZATION)
             call hepmc_particle_set_polarization (hprt(i), &
                  particle_set%prt(i)%get_polarization ())
          end select
       end do
    end if
  end subroutine hepmc_event_from_particle_set

  subroutine particle_from_hepmc_particle &
       (prt, hprt, model, fallback_model, polarization, barcode)
    type(particle_t), intent(out) :: prt
    type(hepmc_particle_t), intent(in) :: hprt
    type(model_data_t), intent(in), target :: model
    type(model_data_t), intent(in), target :: fallback_model
    type(hepmc_vertex_t) :: vtx
    integer, intent(in) :: polarization
    integer, dimension(:), intent(in) :: barcode
    type(hepmc_polarization_t) :: hpol
    type(flavor_t) :: flv
    type(color_t) :: col
    type(helicity_t) :: hel
    type(polarization_t) :: pol
    type(vector4_t) :: vertex
    integer :: n_parents, n_children
    integer, dimension(:), allocatable :: &
         parent_barcode, child_barcode, parent, child
    integer :: i
    select case (hepmc_particle_get_status (hprt))
    case (1);  call prt%set_status (PRT_OUTGOING)
    case (2);  call prt%set_status (PRT_RESONANT)
    case (3);  call prt%set_status (PRT_VIRTUAL)
    end select
    if (hepmc_particle_is_beam (hprt)) call prt%set_status (PRT_BEAM)
    call flv%init (hepmc_particle_get_pdg (hprt), model, fallback_model)
    call col%init (hepmc_particle_get_color (hprt))
    call prt%set_flavor (flv)
    call prt%set_color (col)
    call prt%set_polarization (polarization)
    select case (polarization)
    case (PRT_DEFINITE_HELICITY)
       hpol = hepmc_particle_get_polarization (hprt)
       call hepmc_polarization_to_hel (hpol, prt%get_flv (), hel)
       call prt%set_helicity (hel)
       call hepmc_polarization_final (hpol)
    case (PRT_GENERIC_POLARIZATION)
       hpol = hepmc_particle_get_polarization (hprt)
       call hepmc_polarization_to_pol (hpol, prt%get_flv (), pol)
       call prt%set_pol (pol)
       call hepmc_polarization_final (hpol)
    end select
    call prt%set_momentum (hepmc_particle_get_momentum (hprt), &
         hepmc_particle_get_mass_squared (hprt))
    n_parents  = hepmc_particle_get_n_parents  (hprt)
    n_children = hepmc_particle_get_n_children (hprt)
    if (HEPMC2_AVAILABLE) then
       allocate (parent_barcode (n_parents),  parent (n_parents))
       allocate (child_barcode  (n_children), child  (n_children))
       parent_barcode = hepmc_particle_get_parent_barcodes (hprt)
       child_barcode  = hepmc_particle_get_child_barcodes  (hprt)
       do i = 1, size (barcode)
          where (parent_barcode == barcode(i))  parent = i
          where (child_barcode  == barcode(i))  child  = i
       end do
       call prt%set_parents (parent)
       call prt%set_children (child)
    else if (HEPMC3_AVAILABLE) then
       allocate (parent_barcode (n_parents),  parent (n_parents))
       allocate (child_barcode  (n_children), child  (n_children))
       parent_barcode = hepmc_particle_get_parent_barcodes (hprt)
       child_barcode  = hepmc_particle_get_child_barcodes  (hprt)
       do i = 1, size (barcode)
          where (parent_barcode == barcode(i))  parent = i
          where (child_barcode  == barcode(i))  child  = i
       end do
       call prt%set_parents (parent)
       call prt%set_children (child)
    end if
    if (prt%get_status () == PRT_VIRTUAL .and. n_parents == 0) &
         call prt%set_status (PRT_INCOMING)
    if (HEPMC2_AVAILABLE) then
       vtx = hepmc_particle_get_decay_vertex (hprt)
       if (hepmc_vertex_is_valid (vtx)) then
          vertex = hepmc_vertex_to_vertex (vtx)
          if (vertex /= vector4_null)  call prt%set_vertex (vertex)
       end if
    end if
  end subroutine particle_from_hepmc_particle

  module subroutine hepmc_event_to_particle_set &
       (particle_set, evt, model, fallback_model, polarization)
    type(particle_set_t), intent(inout), target :: particle_set
    type(hepmc_event_t), intent(in) :: evt
    class(model_data_t), intent(in), target :: model, fallback_model
    integer, intent(in) :: polarization
    type(hepmc_event_particle_iterator_t) :: it
    type(hepmc_vertex_t) :: v
    type(hepmc_vertex_particle_in_iterator_t) :: v_it
    type(hepmc_particle_t) :: prt
    integer, dimension(:), allocatable :: barcode, n_parents
    integer :: n_tot, n_beam, i, bc
    n_tot = hepmc_event_get_n_particles(evt)
    allocate (barcode (n_tot))
    if (HEPMC2_AVAILABLE) then
       call hepmc_event_particle_iterator_init (it, evt)
       do i = 1, n_tot
          barcode(i) = hepmc_particle_get_barcode &
               (hepmc_event_particle_iterator_get (it))
          call hepmc_event_particle_iterator_advance (it)
       end do
       allocate (particle_set%prt (n_tot))
       call hepmc_event_particle_iterator_reset (it)
       do i = 1, n_tot
          prt = hepmc_event_particle_iterator_get (it)
          call particle_from_hepmc_particle (particle_set%prt(i), &
               prt, model, fallback_model, polarization, barcode)
          call hepmc_event_particle_iterator_advance (it)
       end do
       call hepmc_event_particle_iterator_final (it)
       v = hepmc_event_get_signal_process_vertex (evt)
       if (hepmc_vertex_is_valid (v)) then
          call hepmc_vertex_particle_in_iterator_init (v_it, v)
          do while (hepmc_vertex_particle_in_iterator_is_valid (v_it))
             prt = hepmc_vertex_particle_in_iterator_get (v_it)
             bc = hepmc_particle_get_barcode &
                  (hepmc_vertex_particle_in_iterator_get (v_it))
             do i = 1, size(barcode)
                if (bc == barcode(i))  &
                     call particle_set%prt(i)%set_status (PRT_INCOMING)
             end do
             call hepmc_vertex_particle_in_iterator_advance (v_it)
          end do
          call hepmc_vertex_particle_in_iterator_final (v_it)
       end if
    else if (HEPMC3_AVAILABLE) then
       allocate (particle_set%prt (n_tot))
       do i = 1, n_tot
          barcode(i) = hepmc_particle_get_barcode &
               (hepmc_event_get_nth_particle (evt, i))
       end do
       do i = 1, n_tot
          prt = hepmc_event_get_nth_particle (evt, i)
          call particle_from_hepmc_particle (particle_set%prt(i), &
               prt, model, fallback_model, polarization, barcode)
       end do
    end if
    do i = 1, n_tot
       if (particle_set%prt(i)%get_status () == PRT_VIRTUAL &
            .and. particle_set%prt(i)%get_n_children () == 0) &
            call particle_set%prt(i)%set_status (PRT_OUTGOING)
    end do
    if (HEPMC3_AVAILABLE) then
       n_beam = hepmc_event_get_n_beams (evt)
       do i = 1, n_beam
          bc = hepmc_event_get_nth_beam (evt, i)
          if (.not. particle_set%prt(bc)%get_status () == PRT_INCOMING) &
               call particle_set%prt(bc)%set_status (PRT_BEAM)
       end do
       do i = 1, n_tot
          if (particle_set%prt(i)%get_status () == PRT_VIRTUAL) then
             n_parents = particle_set%prt(i)%get_parents ()
             if (all &
                  (particle_set%prt(n_parents)%get_status () == PRT_BEAM)) then
                call particle_set%prt(i)%set_status (PRT_INCOMING)
             end if
          end if
       end do
    end if
    particle_set%n_tot = n_tot
    particle_set%n_beam  = &
         count (particle_set%prt%get_status () == PRT_BEAM)
    particle_set%n_in  = &
         count (particle_set%prt%get_status () == PRT_INCOMING)
    particle_set%n_out = &
         count (particle_set%prt%get_status () == PRT_OUTGOING)
    particle_set%n_vir = &
         particle_set%n_tot - particle_set%n_in - particle_set%n_out
  end subroutine hepmc_event_to_particle_set

  module subroutine hepmc_to_event &
       (event, hepmc_event, fallback_model, process_index, &
       recover_beams, use_alpha_s, use_scale)
    class(generic_event_t), intent(inout), target :: event
    type(hepmc_event_t), intent(inout) :: hepmc_event
    class(model_data_t), intent(in), target :: fallback_model
    integer, intent(out), optional :: process_index
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: use_alpha_s
    logical, intent(in), optional :: use_scale
    class(model_data_t), pointer :: model
    real(default) :: scale, alpha_qcd
    type(particle_set_t) :: particle_set
    model => event%get_model_ptr ()
    call event%set_index (hepmc_event_get_event_index (hepmc_event))
    call hepmc_event_to_particle_set (particle_set, &
         hepmc_event, model, fallback_model, PRT_DEFINITE_HELICITY)
    call event%set_hard_particle_set (particle_set)
    call particle_set%final ()
    call event%set_weight_ref (1._default)
    alpha_qcd = hepmc_event_get_alpha_qcd (hepmc_event)
    scale = hepmc_event_get_scale (hepmc_event)
    if (present (use_alpha_s)) then
       if (use_alpha_s .and. alpha_qcd > 0) &
            call event%set_alpha_qcd_forced (alpha_qcd)
    end if
    if (present (use_scale)) then
       if (use_scale .and. scale > 0) &
            call event%set_scale_forced (scale)
    end if
  end subroutine hepmc_to_event

  module subroutine particle_to_lcio (prt, lprt)
    type(particle_t), intent(in) :: prt
    type(lcio_particle_t), intent(out) :: lprt
    integer :: lcio_status
    type(vector4_t) :: vtx
    select case (prt%get_status ())
    case (PRT_UNDEFINED)
       lcio_status = 0
    case (PRT_OUTGOING)
       lcio_status = 1
    case (PRT_BEAM_REMNANT)
       if (prt%get_n_children () == 0) then
          lcio_status = 1
       else
          lcio_status = 3
       end if
    case (PRT_BEAM)
       lcio_status = 4
    case (PRT_RESONANT)
       lcio_status = 2
    case default
       lcio_status = 3
    end select
    call lcio_particle_init (lprt, &
         prt%get_momentum (), &
         prt%get_pdg (), &
         prt%flv%get_charge (), &
         lcio_status)
    call lcio_particle_set_color (lprt, prt%get_color ())
    vtx = prt%get_vertex ()
    call lcio_particle_set_vtx (lprt, space_part (vtx))
    call lcio_particle_set_t (lprt, vtx%p(0))
    select case (prt%get_polarization_status ())
    case (PRT_DEFINITE_HELICITY)
       call lcio_polarization_init (lprt, prt%get_helicity ())
    case (PRT_GENERIC_POLARIZATION)
       call lcio_polarization_init (lprt, prt%get_polarization ())
    end select
  end subroutine particle_to_lcio

  module subroutine particle_from_lcio_particle &
     (prt, lprt, model, fallback_model, daughters, parents, polarization)
    type(particle_t), intent(out) :: prt
    type(lcio_particle_t), intent(in) :: lprt
    type(model_data_t), intent(in), target :: model
    type(model_data_t), intent(in), target :: fallback_model
    integer, dimension(:), intent(in) :: daughters, parents
    integer, intent(in) :: polarization
    type(vector4_t) :: vtx4
    type(flavor_t) :: flv
    type(color_t) :: col
    type(helicity_t) :: hel
    type(polarization_t) :: pol
    select case (lcio_particle_get_status (lprt))
    case (1);  call prt%set_status (PRT_OUTGOING)
    case (2);  call prt%set_status (PRT_RESONANT)
    case (3)
       select case (size (parents))
       case (0)
          call prt%set_status (PRT_INCOMING)
       case default
          call prt%set_status (PRT_VIRTUAL)
       end select
    case (4);  call prt%set_status (PRT_BEAM)
    end select
    call flv%init (lcio_particle_get_pdg (lprt), model, fallback_model)
    call col%init (lcio_particle_get_flow (lprt))
    if (flv%is_beam_remnant ())  call prt%set_status (PRT_BEAM_REMNANT)
    call prt%set_flavor (flv)
    call prt%set_color (col)
    call prt%set_polarization (polarization)
    select case (polarization)
    case (PRT_DEFINITE_HELICITY)
       call lcio_particle_to_hel (lprt, prt%get_flv (), hel)
       call prt%set_helicity (hel)
    case (PRT_GENERIC_POLARIZATION)
       call lcio_particle_to_pol (lprt, prt%get_flv (), pol)
       call prt%set_pol (pol)
    end select
    call prt%set_momentum (lcio_particle_get_momentum (lprt), &
         lcio_particle_get_mass_squared (lprt))
    call prt%set_parents (parents)
    call prt%set_children (daughters)
    vtx4 = vector4_moving (lcio_particle_get_time (lprt), &
         lcio_particle_get_vertex (lprt))
    if (vtx4 /= vector4_null)  call prt%set_vertex (vtx4)
  end subroutine particle_from_lcio_particle

  module subroutine lcio_event_from_particle_set (evt, particle_set)
    type(lcio_event_t), intent(inout) :: evt
    type(particle_set_t), intent(in) :: particle_set
    type(lcio_particle_t), dimension(:), allocatable :: lprt
    type(particle_set_t), target :: pset_filtered
    integer, dimension(:), allocatable :: parent
    integer :: n_tot, i, j, n_beam, n_parents, type, beam_count

    call particle_set%filter_particles ( pset_filtered, &
         real_parents = .true. , keep_beams = .true. , keep_virtuals = .false.)
    n_tot = pset_filtered%n_tot
    n_beam = count (pset_filtered%prt%get_status () == PRT_BEAM)
    if (n_beam == 0) then
       type = PRT_INCOMING
    else
       type = PRT_BEAM
    end if
    beam_count = 0
    allocate (lprt (n_tot))
    do i = 1, n_tot
       call particle_to_lcio (pset_filtered%prt(i), lprt(i))
       n_parents = pset_filtered%prt(i)%get_n_parents ()
       if (n_parents /= 0) then
          allocate (parent (n_parents))
          parent = pset_filtered%prt(i)%get_parents ()
          do j = 1, n_parents
             call lcio_particle_set_parent  (lprt(i), lprt(parent(j)))
          end do
          deallocate (parent)
       end if
       if (pset_filtered%prt(i)%get_status () == type) then
          beam_count = beam_count + 1
          call lcio_event_set_beam &
               (evt, pset_filtered%prt(i)%get_pdg (), beam_count)
       end if
       call lcio_particle_add_to_evt_coll (lprt(i), evt)
    end do
    call lcio_event_add_coll (evt)
  end subroutine lcio_event_from_particle_set

  module subroutine lcio_event_to_particle_set &
       (particle_set, evt, model, fallback_model, polarization)
    type(particle_set_t), intent(inout), target :: particle_set
    type(lcio_event_t), intent(in) :: evt
    class(model_data_t), intent(in), target :: model, fallback_model
    integer, intent(in) :: polarization
    type(lcio_particle_t) :: prt
    integer, dimension(:), allocatable :: parents, daughters
    integer :: n_tot, i, j, n_parents, n_children
    n_tot = lcio_event_get_n_tot (evt)
    allocate (particle_set%prt (n_tot))
    do i = 1, n_tot
       prt = lcio_event_get_particle (evt, i-1)
       n_parents = lcio_particle_get_n_parents (prt)
       n_children = lcio_particle_get_n_children (prt)
       allocate (daughters (n_children))
       allocate (parents (n_parents))
       if (n_children > 0) then
          do j = 1, n_children
             daughters(j) = lcio_get_n_children (evt,i,j)
          end do
       end if
       if (n_parents > 0) then
          do j = 1, n_parents
             parents(j) = lcio_get_n_parents (evt,i,j)
          end do
       end if
       call particle_from_lcio_particle (particle_set%prt(i), &
            prt, model, fallback_model, &
            daughters, parents, polarization)
       deallocate (daughters, parents)
    end do
    do i = 1, n_tot
       if (particle_set%prt(i)%get_status () == PRT_VIRTUAL) then
          CHECK_BEAM: do j = 1, particle_set%prt(i)%get_n_parents ()
             if (particle_set%prt(j)%get_status () == PRT_BEAM) &
                  call particle_set%prt(i)%set_status (PRT_INCOMING)
             exit CHECK_BEAM
          end do CHECK_BEAM
       end if
    end do
    particle_set%n_tot = n_tot
    particle_set%n_beam  = &
         count (particle_set%prt%get_status () == PRT_BEAM)
    particle_set%n_in  = &
         count (particle_set%prt%get_status () == PRT_INCOMING)
    particle_set%n_out = &
         count (particle_set%prt%get_status () == PRT_OUTGOING)
    particle_set%n_vir = &
         particle_set%n_tot - particle_set%n_in - particle_set%n_out
  end subroutine lcio_event_to_particle_set

  module subroutine lcio_to_event &
       (event, lcio_event, fallback_model, process_index, recover_beams, &
       use_alpha_s, use_scale)
    class(generic_event_t), intent(inout), target :: event
    type(lcio_event_t), intent(inout) :: lcio_event
    class(model_data_t), intent(in), target :: fallback_model
    integer, intent(out), optional :: process_index
    logical, intent(in), optional :: recover_beams
    logical, intent(in), optional :: use_alpha_s
    logical, intent(in), optional :: use_scale
    class(model_data_t), pointer :: model
    real(default) :: scale, alpha_qcd
    type(particle_set_t) :: particle_set
    model => event%get_model_ptr ()
    call lcio_event_to_particle_set (particle_set, &
         lcio_event, model, fallback_model, PRT_DEFINITE_HELICITY)
    call event%set_hard_particle_set (particle_set)
    call particle_set%final ()
    call event%set_weight_ref (1._default)
    alpha_qcd = lcio_event_get_alphas (lcio_event)
    scale = lcio_event_get_scaleval (lcio_event)
    if (present (use_alpha_s)) then
       if (use_alpha_s .and. alpha_qcd > 0) &
            call event%set_alpha_qcd_forced (alpha_qcd)
    end if
    if (present (use_scale)) then
       if (use_scale .and. scale > 0) &
            call event%set_scale_forced (scale)
    end if
  end subroutine lcio_to_event


end submodule hep_events_s

