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

submodule (shower_core) shower_core_s

  use io_units
  use constants
  use format_utils, only: write_separator
  use numeric_utils
  use diagnostics
  use physics_defs
  use sm_physics
  use hep_common
  use flavors
  use colors
  use subevents

  implicit none

  real(default), save :: alphasxpdfmax = 12._default


contains

  module subroutine shower_init &
       (shower, settings, taudec_settings, pdf_data, os_data)
    class(shower_t), intent(out) :: shower
    type(shower_settings_t), intent(in) :: settings
    type(taudec_settings_t), intent(in) :: taudec_settings
    type(pdf_data_t), intent(in) :: pdf_data
    type(os_data_t), intent(in) :: os_data
    if (debug_on) call msg_debug (D_SHOWER, "shower_init")
    shower%settings = settings
    shower%taudec_settings = taudec_settings
    shower%os_data = os_data
    call shower%pdf_data%init (pdf_data)
    shower%name = "WHIZARD internal"
    call shower%write_msg ()
  end subroutine shower_init

  module subroutine shower_prepare_new_event (shower, fac_scale, alpha_s)
    class(shower_t), intent(inout) :: shower
    real(default), intent(in) :: fac_scale, alpha_s
    call shower%cleanup ()
    shower%next_free_nr = 1
    shower%next_color_nr = 1
    if (debug_active (D_SHOWER)) then
       if (allocated (shower%interactions)) then
          call msg_bug ("Shower: creating new shower while old one " // &
               "is still associated (interactions)")
       end if
       if (allocated (shower%partons)) then
          call msg_bug ("Shower: creating new shower while old one " // &
               "is still associated (partons)")
       end if
    end if
    treat_light_quarks_massless = .true.
    treat_duscb_quarks_massless = .false.
    shower%valid = .true.
  end subroutine shower_prepare_new_event

  module subroutine shower_activate_multiple_interactions (shower)
    class(shower_t), intent(inout) :: shower
    if (shower%mi%is_initialized ()) then
       call shower%mi%restart ()
    else
       call shower%mi%initialize (&
            GeV2_scale_cutoff=shower%settings%min_virtuality, &
            GeV2_s=shower_interaction_get_s &
            (shower%interactions(1)%i), &
            muli_dir=char(shower%os_data%whizard_mulipath))
    end if
    call shower%mi%apply_initial_interaction ( &
         GeV2_s=shower_interaction_get_s(shower%interactions(1)%i), &
         x1=shower%interactions(1)%i%partons(1)%p%parent%x, &
         x2=shower%interactions(1)%i%partons(2)%p%parent%x, &
         pdg_f1=shower%interactions(1)%i%partons(1)%p%parent%type, &
         pdg_f2=shower%interactions(1)%i%partons(2)%p%parent%type, &
         n1=shower%interactions(1)%i%partons(1)%p%parent%nr, &
         n2=shower%interactions(1)%i%partons(2)%p%parent%nr)
  end subroutine shower_activate_multiple_interactions

  module subroutine shower_import_particle_set (shower, particle_set)
    class(shower_t), target, intent(inout) :: shower
    type(particle_set_t), intent(in) :: particle_set
    !integer, dimension(:), allocatable :: connections
    type(parton_t), dimension(:), allocatable, target, save :: partons, hadrons
    type(parton_pointer_t), dimension(:), allocatable :: &
         parton_pointers
    integer :: n_beam, n_in, n_out, n_tot
    integer :: i, j, nr, max_color_nr
    if (debug_on) call msg_debug (D_SHOWER, 'shower_import_particle_set')
    call count_and_allocate ()
    call setup_hadrons_from_particle_set ()
    call setup_partons_from_particle_set ()
    call shower%update_max_color_nr (1 + max_color_nr)
    call shower%add_interaction_2ton (parton_pointers)
    if (shower%settings%muli_active) then
       call shower%activate_multiple_interactions ()
    end if
    if (debug_on) call msg_debug2 (D_SHOWER, 'shower%write() after shower_import_particle_set')
    if (debug2_active (D_SHOWER)) then
       call shower%write ()
    end if
  contains
    subroutine count_and_allocate ()
      max_color_nr = 0
      n_beam = particle_set%get_n_beam ()
      n_in = particle_set%get_n_in ()
      n_out = particle_set%get_n_out ()
      n_tot = particle_set%get_n_tot ()
      if (allocated (partons))  deallocate (partons)
      allocate (partons (n_in + n_out))
      allocate (parton_pointers (n_in+n_out))
    end subroutine count_and_allocate

    subroutine setup_hadrons_from_particle_set ()
      j = 0
      !!! !!! !!! Workaround for Portland 16.1 compiler bug
      !!! if (n_beam > 0 .and. all (particle_set%prt(1:2)%flv%get_pdg_abs () > TAU)) then
      if (n_beam > 0 .and. particle_set%prt(1)%flv%get_pdg_abs () > TAU .and. &
           particle_set%prt(2)%flv%get_pdg_abs () > TAU) then
         if (debug_on) call msg_debug (D_SHOWER, 'Copy hadrons from particle_set to hadrons')
         if (.not. allocated (hadrons))  allocate (hadrons (1:2))
         do i = 1, n_tot
            if (particle_set%prt(i)%status == PRT_BEAM) then
               j = j + 1
               nr = shower%get_next_free_nr ()
               hadrons(j) = parton_of_particle (particle_set%prt(i), nr)
               hadrons(j)%settings => shower%settings
               max_color_nr = max (max_color_nr, abs(hadrons(j)%c1), &
                                                 abs(hadrons(j)%c2))
         end if
         end do
      end if
    end subroutine setup_hadrons_from_particle_set

    subroutine setup_partons_from_particle_set ()
      integer, dimension(1) :: parent
      j = 0
      if (debug_on) call msg_debug (D_SHOWER, "Copy partons from particle_set to partons")
      do i = 1, n_tot
         if (particle_set%prt(i)%get_status () == PRT_INCOMING .or. &
             particle_set%prt(i)%get_status () == PRT_OUTGOING) then
            j = j + 1
            nr = shower%get_next_free_nr ()
            partons(j) = parton_of_particle (particle_set%prt(i), nr)
            partons(j)%settings => shower%settings
            parton_pointers(j)%p => partons(j)
            max_color_nr = max (max_color_nr, abs (partons(j)%c1), &
                                              abs (partons(j)%c2))
            if (particle_set%prt(i)%get_status () == PRT_INCOMING .and. &
                particle_set%prt(i)%get_n_parents () == 1 .and. &
                allocated (hadrons)) then
               parent = particle_set%prt(i)%get_parents ()
               partons(j)%initial => hadrons (parent(1))
               partons(j)%x = space_part_norm (partons(j)%momentum) / &
                              space_part_norm (partons(j)%initial%momentum)
                  end if
         end if
      end do
    end subroutine setup_partons_from_particle_set

  end subroutine shower_import_particle_set

  module subroutine shower_generate_emissions &
       (shower, valid, number_of_emissions)
    class(shower_t), intent(inout), target :: shower
    logical, intent(out) :: valid
    integer, optional, intent(in) :: number_of_emissions

    type(parton_t), dimension(:), allocatable, target :: partons
    type(parton_pointer_t), dimension(:), allocatable :: &
         parton_pointers
    real(default) :: mi_scale, ps_scale, shat, phi
    type(parton_pointer_t) :: temppp
    integer :: i, j, k
    integer :: n_int, max_color_nr
    integer, dimension(2,4) :: color_corr
    if (debug_on) call msg_debug (D_SHOWER, "shower_generate_emissions")
    if (shower%settings%isr_active) then
       if (debug_on) call msg_debug (D_SHOWER, "Generate ISR with FSR")
       i = 0
       BRANCHINGS: do
          i = i + 1
          if (signal_is_pending ()) return
          if (shower%settings%muli_active) then
          call shower%mi%generate_gev2_pt2 &
                  (shower%get_ISR_scale (), mi_scale)
          else
             mi_scale = 0.0
          end if

          !!! Shower: debugging
          !!! shower%generate_next_isr_branching returns a pointer to
          !!! the parton with the next ISR-branching, this parton's
          !!! scale is the scale of the next branching
          ! temppp=shower%generate_next_isr_branching_veto ()
          temppp = shower%generate_next_isr_branching ()

          if (.not. associated (temppp%p) .and. &
            mi_scale < shower%settings%min_virtuality) then
             exit BRANCHINGS
          end if
          !!! check if branching or interaction occurs next
          if (associated (temppp%p)) then
             ps_scale = abs(temppp%p%t)
          else
             ps_scale = 0._default
          end if
          if (mi_scale > ps_scale) then
             !!! discard branching evolution lower than mi_scale
             call shower%set_max_ISR_scale (mi_scale)
             if (associated (temppp%p)) &
                  call temppp%p%set_simulated (.false.)

             !!! execute new interaction
             deallocate (partons)
             deallocate (parton_pointers)
             allocate (partons(1:4))
             allocate (parton_pointers(1:4))
             do j = 1, 4
                partons(j)%nr = shower%get_next_free_nr ()
                partons(j)%belongstointeraction = .true.
                parton_pointers(j)%p => partons(j)
             end do
          call shower%mi%generate_partons (partons(1)%nr, partons(2)%nr, &
                  partons(1)%x, partons(2)%x, &
                  partons(1)%type, partons(2)%type, &
                  partons(3)%type, partons(4)%type)
             !!! calculate momenta
             shat = partons(1)%x *partons(2)%x * &
                  shower_interaction_get_s (shower%interactions(1)%i)
             partons(1)%momentum = [0.5_default * sqrt(shat), &
                  zero, zero, 0.5_default*sqrt(shat)]
             partons(2)%momentum = [0.5_default * sqrt(shat), &
                  zero, zero, -0.5_default*sqrt(shat)]
             call parton_set_initial (partons(1), &
                  shower%interactions(1)%i%partons(1)%p%initial)
             call parton_set_initial (partons(2), &
                  shower%interactions(1)%i%partons(2)%p%initial)
             partons(1)%belongstoFSR = .false.
             partons(2)%belongstoFSR = .false.
             !!! calculate color connection
          call shower%mi%get_color_correlations &
                  (shower%get_next_color_nr (), &
                  max_color_nr,color_corr)
          call shower%update_max_color_nr (max_color_nr)

             partons(1)%c1 = color_corr(1,1)
             partons(1)%c2 = color_corr(2,1)
             partons(2)%c1 = color_corr(1,2)
             partons(2)%c2 = color_corr(2,2)
             partons(3)%c1 = color_corr(1,3)
             partons(3)%c2 = color_corr(2,3)
             partons(4)%c1 = color_corr(1,4)
             partons(4)%c2 = color_corr(2,4)

             call shower%rng%generate (phi)
             phi = 2 * pi * phi
             partons(3)%momentum = [0.5_default*sqrt(shat), &
                  sqrt(mi_scale)*cos(phi), &
                  sqrt(mi_scale)*sin(phi), &
                  sqrt(0.25_default*shat - mi_scale)]
             partons(4)%momentum = [ 0.5_default*sqrt(shat), &
                  -sqrt(mi_scale)*cos(phi), &
                  -sqrt(mi_scale)*sin(phi), &
                  -sqrt(0.25_default*shat - mi_scale)]
             partons(3)%belongstoFSR = .true.
             partons(4)%belongstoFSR = .true.

             call shower%add_interaction_2ton (parton_pointers)
             n_int = size (shower%interactions)
             do k = 1, 2
             call shower%mi%replace_parton &
                  (shower%interactions(n_int)%i%partons(k)%p%initial%nr, &
                   shower%interactions(n_int)%i%partons(k)%p%nr, &
                   shower%interactions(n_int)%i%partons(k)%p%parent%nr, &
                   shower%interactions(n_int)%i%partons(k)%p%type, &
                   shower%interactions(n_int)%i%partons(k)%p%x, &
                   mi_scale)
             end do
          else
          !!! execute the next branching 'found' in the previous step
             call shower%execute_next_isr_branching (temppp)
             if (shower%settings%muli_active) then
             call shower%mi%replace_parton (temppp%p%initial%nr, &
                     temppp%p%child1%nr, temppp%p%nr, &
                     temppp%p%type, temppp%p%x, ps_scale)
             end if

          end if
       end do BRANCHINGS

       call shower%generate_fsr_for_isr_partons ()
    else
       if (signal_is_pending ()) return
       if (debug_on) call msg_debug (D_SHOWER, "Generate FSR without ISR")
       call shower%simulate_no_isr_shower ()
    end if

    !!! some bookkeeping, needed after the shower is done
    call shower%boost_to_labframe ()
    call shower%generate_primordial_kt ()
    call shower%update_beamremnants ()

    if (shower%settings%fsr_active) then
       do i = 1, size (shower%interactions)
          if (signal_is_pending ()) return
          call shower%interaction_generate_fsr_2ton &
               (shower%interactions(i)%i)
       end do
    else
       call shower%simulate_no_fsr_shower ()
    end if
    if (debug_on) call msg_debug (D_SHOWER, "Shower finished:")
    if (debug_active (D_SHOWER))  call shower%write ()

    valid = shower%valid
    !!! clean-up muli: we should finalize the muli pdf sets
    !!!      when _all_ runs are done. Not after every event if possible
    ! call shower%mi%finalize()
  end subroutine shower_generate_emissions

  module subroutine shower_make_particle_set &
       (shower, particle_set, model, model_hadrons)
    class(shower_t), intent(in) :: shower
    type(particle_set_t), intent(inout) :: particle_set
    class(model_data_t), intent(in), target :: model
    class(model_data_t), intent(in), target :: model_hadrons
    call shower%combine_with_particle_set (particle_set, model, &
         model_hadrons)
    if (shower%settings%hadronization_active) then
       call shower%converttopythia ()
    end if
  end subroutine shower_make_particle_set

  module subroutine shower_add_interaction_2ton (shower, partons)
    class(shower_t), intent(inout) :: shower
    type(parton_pointer_t), intent(in), dimension(:), allocatable :: partons
    !type(ckkw_pseudo_shower_weights_t), intent(in) :: ckkw_pseudo_weights

    integer :: n_partons, n_out
    integer :: i, j, imin, jmin
    real(default) :: y, ymin
    !real(default) :: w, wmax
    !real(default) :: random, sum
    type(parton_pointer_t), dimension(:), allocatable :: new_partons
    type(parton_t), pointer :: prt
    integer :: n_int
    type(shower_interaction_pointer_t), dimension(:), allocatable :: temp
    type(vector4_t) :: prtmomentum, childmomentum
    logical :: isr_is_possible_and_allowed
    type(lorentz_transformation_t) :: L

    if (signal_is_pending ()) return
    if (debug_on) call msg_debug (D_SHOWER, "Add interaction2toN")
    n_partons = size (partons)
    n_out = n_partons - 2
    if (n_out < 2) then
       call msg_bug &
            ("Shower core: trying to add a 2-> (something<2) interaction")
    end if

    isr_is_possible_and_allowed = (associated (partons(1)%p%initial) &
         .and. associated (partons(2)%p%initial)) .and. &
         shower%settings%isr_active
    if (debug_on) call msg_debug (D_SHOWER, "isr_is_possible_and_allowed", &
         isr_is_possible_and_allowed)

    if (associated (partons(1)%p%initial) .and. &
         partons(1)%p%is_quark ()) then
       if (partons(1)%p%momentum%p(0) < &
            two * partons(1)%p%mass()) then
          if (abs(partons(1)%p%type) < 2) then
             treat_light_quarks_massless = .true.
          else
             treat_duscb_quarks_massless = .true.
          end if
       end if
    end if
    if (associated (partons(2)%p%initial) .and. &
         partons(2)%p%is_quark ()) then
       if (partons(2)%p%momentum%p(0) < &
            two * partons(2)%p%mass()) then
          if (abs(partons(2)%p%type) < 2) then
             treat_light_quarks_massless = .true.
          else
             treat_duscb_quarks_massless = .true.
          end if
       end if
    end if

    if (allocated (shower%interactions)) then
       n_int = size (shower%interactions) + 1
    else
       n_int = 1
    end if
    allocate (temp (1:n_int))
    do i = 1, n_int - 1
       allocate (temp(i)%i)
       temp(i)%i = shower%interactions(i)%i
    end do
    allocate (temp(n_int)%i)
    allocate (temp(n_int)%i%partons(1:n_partons))
    do i = 1, n_partons
       allocate (temp(n_int)%i%partons(i)%p)
       call parton_copy (partons(i)%p, temp(n_int)%i%partons(i)%p)
    end do
    if (allocated (shower%interactions)) deallocate(shower%interactions)
    allocate (shower%interactions(1:n_int))
    do i = 1, n_int
       shower%interactions(i)%i => temp(i)%i
    end do
    deallocate (temp)

    if (associated (shower%interactions(n_int)%i%partons(1)%p%initial))  &
         call shower%interactions(n_int)%i%partons(1)%p%initial%set_simulated ()
    if (associated (shower%interactions(n_int)%i%partons(2)%p%initial))  &
         call shower%interactions(n_int)%i%partons(2)%p%initial%set_simulated ()
    if (isr_is_possible_and_allowed) then
       !!! boost to the CMFrame of the incoming partons
       L = boost (-(shower%interactions(n_int)%i%partons(1)%p%momentum + &
            shower%interactions(n_int)%i%partons(2)%p%momentum), &
            (shower%interactions(n_int)%i%partons(1)%p%momentum + &
            shower%interactions(n_int)%i%partons(2)%p%momentum)**1 )
       do i = 1, n_partons
          call parton_apply_lorentztrafo &
               (shower%interactions(n_int)%i%partons(i)%p, L)
       end do
    end if
    do i = 1, size (partons)
       if (signal_is_pending ()) return
       !!! partons are marked as belonging to the hard interaction
       shower%interactions(n_int)%i%partons(i)%p%belongstointeraction &
            = .true.
       shower%interactions(n_int)%i%partons(i)%p%belongstoFSR = i > 2
       shower%interactions(n_int)%i%partons(i)%p%interactionnr = n_int
       !!! include a 2^(i - 1) number as a label for the ckkw clustering
       shower%interactions(n_int)%i%partons(i)%p%ckkwlabel = 2**(i - 1)
    end do

    if (allocated (shower%partons)) then
       allocate (new_partons(1:size(shower%partons) + &
            size(shower%interactions(n_int)%i%partons)))
       do i = 1, size (shower%partons)
          new_partons(i)%p => shower%partons(i)%p
       end do
       do i = 1, size (shower%interactions(n_int)%i%partons)
          new_partons(size(shower%partons) + i)%p => &
               shower%interactions(n_int)%i%partons(i)%p
       end do
       deallocate (shower%partons)
    else
       allocate (new_partons(1:size(shower%interactions(n_int)%i%partons)))
       do i = 1, size (partons)
          new_partons(i)%p => shower%interactions(n_int)%i%partons(i)%p
       end do
    end if
    allocate (shower%partons(1:size (new_partons)))
    do i = 1, size (new_partons)
       shower%partons(i)%p => new_partons(i)%p
    end do
    deallocate (new_partons)

    if (isr_is_possible_and_allowed) then
       if (shower%settings%isr_pt_ordered) then
          call shower_prepare_for_simulate_isr_pt &
               (shower, shower%interactions(size (shower%interactions))%i)
       else
          call shower_prepare_for_simulate_isr_ana_test &
               (shower, shower%interactions(n_int)%i%partons(1)%p, &
               shower%interactions(n_int)%i%partons(2)%p)
       end if
    end if

    !!! generate pseudo parton shower history and add all partons to
    !!!          shower%partons-array
    !!! TODO initial -> initial + final branchings ??
    allocate (new_partons(1:(n_partons - 2)))
    do i = 1, size (new_partons)
       nullify (new_partons(i)%p)
    end do
    do i = 1, size (new_partons)
       new_partons(i)%p => shower%interactions(n_int)%i%partons(i + 2)%p
    end do
    imin = 0
    jmin = 0

    ! TODO: (bcn 2015-04-24) make this a clustering step of the matching
!    if (allocated (ckkw_pseudo_weights%weights)) then
!       !<Perform clustering using the CKKW weights>>
!    else
       CLUSTERING: do
          !!! search for the partons to be clustered together
          ymin = zero
          OUTER: do i = 1, size (new_partons)
             INNER: do j = i + 1, size (new_partons)
                !!! calculate the jet measure
                if (.not.associated (new_partons(i)%p)) cycle INNER
                if (.not.associated (new_partons(j)%p)) cycle INNER
                !if (.not. shower_clustering_allowed &
                     !(shower, new_partons, i,j)) &
                     !cycle inner
                !!! Durham jet-measure ! don't care about constants
                y = min (new_partons(i)%p%momentum%p(0), &
                     new_partons(j)%p%momentum%p(0)) * &
                     (one - enclosed_angle_ct &
                     (new_partons(i)%p%momentum, &
                     new_partons(j)%p%momentum))
                if (y < ymin .or. vanishes(ymin)) then
                   ymin = y
                   imin = i
                   jmin = j
                end if
             end do INNER
          end do OUTER
          if (ymin > zero) then
             call shower%add_parent (new_partons(imin)%p)
             call parton_set_child &
                  (new_partons(imin)%p%parent, new_partons(jmin)%p, 2)
             call parton_set_parent &
                  (new_partons(jmin)%p, new_partons(imin)%p%parent)
             prt => new_partons(imin)%p%parent
             prt%nr = shower_get_next_free_nr (shower)
             prt%type = INTERNAL

             prt%momentum = new_partons(imin)%p%momentum + &
                  new_partons(jmin)%p%momentum
             prt%t = prt%momentum**2
             !!! TODO -> calculate costheta and store it for
             !!!       later use in generate_ps

             if (space_part_norm(prt%momentum) > tiny_10) then
                prtmomentum = prt%momentum
                childmomentum = prt%child1%momentum
                prtmomentum = boost (- prt%get_beta () / sqrt(one - &
                     (prt%get_beta ())**2), space_part(prt%momentum) / &
                     space_part_norm(prt%momentum)) * prtmomentum
                childmomentum = boost (- prt%get_beta() / &
                     sqrt(one - &
                     (prt%get_beta ())**2), space_part(prt%momentum) / &
                     space_part_norm(prt%momentum)) * childmomentum
                prt%costheta = enclosed_angle_ct (prtmomentum, childmomentum)
             else
                prt%costheta = - one
             end if

             prt%belongstointeraction = .true.
             nullify (new_partons(imin)%p)
             nullify (new_partons(jmin)%p)
             new_partons(imin)%p => prt
          else
             exit CLUSTERING
          end if
       end do CLUSTERING
!    end if

    !!! set the FSR starting scale for all partons
    do i = 1, size (new_partons)
       !!! the imaginary mother is the only parton remaining in new_partons
       if (.not. associated (new_partons(i)%p)) cycle
       call set_starting_scale (new_partons(i)%p, &
            get_starting_scale (new_partons(i)%p))
       exit
    end do

  contains

    recursive subroutine transfer_pointers (destiny, start, prt)
      type(parton_pointer_t), dimension(:), allocatable :: destiny
      integer, intent(inout) :: start
      type(parton_t), pointer :: prt
      destiny(start)%p => prt
      start = start + 1
      if (associated (prt%child1)) then
         call transfer_pointers (destiny, start, prt%child1)
      end if
      if (associated (prt%child2)) then
         call transfer_pointers (destiny, start, prt%child2)
      end if
    end subroutine transfer_pointers

    recursive function get_starting_scale (prt) result (scale)
      type(parton_t), pointer :: prt
      real(default) :: scale
      scale = huge (scale)
      if (associated (prt%child1) .and. associated (prt%child2)) then
         scale = min(scale, prt%t)
      end if
      if (associated (prt%child1)) then
         scale = min (scale, get_starting_scale (prt%child1))
      end if
      if (associated (prt%child2)) then
         scale = min (scale, get_starting_scale (prt%child2))
      end if
    end function get_starting_scale

    recursive subroutine set_starting_scale (prt, scale)
      type(parton_t), pointer :: prt
      real(default) :: scale
      if (prt%type /= INTERNAL) then
         if (scale > prt%settings%min_virtuality + prt%mass_squared ()) then
            prt%t = scale
         else
            prt%t = prt%mass_squared ()
            call prt%set_simulated ()
         end if
      end if
      if (associated (prt%child1)) then
         call set_starting_scale (prt%child1, scale)
      end if
      if (associated (prt%child2)) then
         call set_starting_scale (prt%child2, scale)
      end if
    end subroutine set_starting_scale

  end subroutine shower_add_interaction_2ton

  module subroutine shower_simulate_no_isr_shower (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i, j
    type(parton_t), pointer :: prt
    if (debug_on) call msg_debug (D_SHOWER, "shower_simulate_no_isr_shower")
    do i = 1, size (shower%interactions)
       do j = 1, 2
          prt => shower%interactions(i)%i%partons(j)%p
          if (associated (prt%initial)) then
             !!! for virtuality ordered: remove unneeded partons
             if (associated (prt%parent)) then
                if (.not. prt%parent%is_proton ()) then
                   if (associated (prt%parent%parent)) then
                      if (.not. prt%parent%is_proton ()) then
                         call shower_remove_parton_from_partons &
                              (shower, prt%parent%parent)
                      end if
                   end if
                   call shower_remove_parton_from_partons &
                        (shower, prt%parent)
                end if
             end if
             call parton_set_parent (prt, prt%initial)
             call parton_set_child (prt%initial, prt, 1)
             if (associated (prt%initial%child2)) then
                call shower_remove_parton_from_partons &
                     (shower,prt%initial%child2)
                deallocate (prt%initial%child2)
             end if
             call shower%add_child (prt%initial, 2)
          end if
       end do
    end do
  end subroutine shower_simulate_no_isr_shower

  module subroutine shower_simulate_no_fsr_shower (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i, j
    type(parton_t), pointer :: prt
    do i = 1, size (shower%interactions)
       do j = 3, size (shower%interactions(i)%i%partons)
          prt => shower%interactions(i)%i%partons(j)%p
          call prt%set_simulated ()
          prt%scale = zero
          prt%t = prt%mass_squared ()
       end do
    end do
  end subroutine shower_simulate_no_fsr_shower

  subroutine swap_pointers (prtp1, prtp2)
    type(parton_pointer_t), intent(inout) :: prtp1, prtp2
    type(parton_pointer_t) :: prtptemp
    prtptemp%p => prtp1%p
    prtp1%p => prtp2%p
    prtp2%p => prtptemp%p
  end subroutine swap_pointers

  recursive subroutine shower_remove_parton_from_partons (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), pointer :: prt
    integer :: i
    if (.not. prt%belongstoFSR .and. associated (prt%child2)) then
       call shower_remove_parton_from_partons_recursive (shower, prt%child2)
    end if
    do i = 1, size (shower%partons)
       if (associated (shower%partons(i)%p, prt)) then
          shower%partons(i)%p => null()
          ! TODO: (bcn 2015-05-05) memory leak here? no deallocation?
          exit
       end if
       if (debug_active (D_SHOWER)) then
          if (i == size (shower%partons)) then
             call msg_bug ("shower_remove_parton_from_partons: parton&
                  &to be removed not found")
          end if
       end if
    end do
  end subroutine shower_remove_parton_from_partons

  recursive subroutine shower_remove_parton_from_partons_recursive (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), pointer :: prt
    if (associated (prt%child1)) then
       call shower_remove_parton_from_partons_recursive (shower, prt%child1)
       deallocate (prt%child1)
    end if
    if (associated (prt%child2)) then
       call shower_remove_parton_from_partons_recursive (shower, prt%child2)
       deallocate (prt%child2)
    end if
    call shower_remove_parton_from_partons (shower, prt)
  end subroutine shower_remove_parton_from_partons_recursive

  module subroutine shower_sort_partons (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i, j, maxsort, size_partons
    logical :: changed
    if (debug_on) call msg_debug2 (D_SHOWER, "shower_sort_partons")
    if (.not. allocated (shower%partons)) return
    size_partons = size (shower%partons)
    maxsort = 0
    do i = 1, size_partons
       if (associated (shower%partons(i)%p))  maxsort = i
    end do
    if (signal_is_pending ()) return
    size_partons = size (shower%partons)
    if (size_partons <= 1) return
    do i = 1, maxsort
       if (.not. associated (shower%partons(i)%p)) cycle
       if (.not. shower%settings%isr_pt_ordered) then
          !!! set unsimulated ISR partons to be "typeless" to prevent
          !!!     influences from "wrong" masses
          if (.not. shower%partons(i)%p%belongstoFSR .and. &
              .not. shower%partons(i)%p%simulated .and. &
              .not. shower%partons(i)%p%belongstointeraction) then
             shower%partons(i)%p%type = 0
          end if
       end if
    end do
    if (signal_is_pending ()) return
    !!! Just a Bubblesort
    !!! Different algorithms needed for t-ordered and pt^2-ordered shower
    !!! Pt-ordered:
    if (shower%settings%isr_pt_ordered) then
       OUTERDO_PT: do i = 1, maxsort - 1
          changed = .false.
          INNERDO_PT: do j = 1, maxsort - i
             if (.not. associated (shower%partons(j + 1)%p)) cycle
             if (.not. associated (shower%partons(j)%p)) then
                !!! change if j + 1 ist assoaciated and j is not
                call swap_pointers (shower%partons(j), shower%partons(j + 1))
                changed = .true.
             else if (shower%partons(j)%p%scale < &
                  shower%partons(j + 1)%p%scale) then
                call swap_pointers (shower%partons(j), shower%partons(j + 1))
                changed = .true.
             else if (nearly_equal(shower%partons(j)%p%scale, &
                  shower%partons(j + 1)%p%scale)) then
                if (shower%partons(j)%p%nr > shower%partons(j + 1)%p%nr) then
                   call swap_pointers (shower%partons(j), shower%partons(j + 1))
                   changed = .true.
                end if
             end if
          end do INNERDO_PT
          if (.not. changed) exit OUTERDO_PT
       end do outerdo_pt
    !!! |t|-ordered
    else
       OUTERDO_T: do i = 1, maxsort - 1
          changed = .false.
          INNERDO_T: do j = 1, maxsort - i
             if (.not. associated (shower%partons(j + 1)%p)) cycle
             if (.not. associated (shower%partons(j)%p)) then
                !!! change if j+1 is associated and j isn't
                call swap_pointers (shower%partons(j), shower%partons(j + 1))
                changed = .true.
             else if (.not. shower%partons(j)%p%belongstointeraction .and. &
                  shower%partons(j + 1)%p%belongstointeraction) then
                !!! move partons belonging to the interaction to the front
                call swap_pointers (shower%partons(j), shower%partons(j + 1))
                changed = .true.
             else if (.not. shower%partons(j)%p%belongstointeraction .and. &
                  .not. shower%partons(j + 1)%p%belongstointeraction ) then
                if (abs (shower%partons(j)%p%t) - &
                     shower%partons(j)%p%mass_squared () < &
                     abs(shower%partons(j + 1)%p%t) - &
                     shower%partons(j + 1)%p%mass_squared ()) then
                   call swap_pointers (shower%partons(j), shower%partons(j + 1))
                   changed = .true.
                else
                   if (nearly_equal(abs (shower%partons(j)%p%t) - &
                        shower%partons(j)%p%mass_squared (), &
                        abs(shower%partons(j + 1)%p%t) - &
                        shower%partons(j + 1)%p%mass_squared ())) then
                      if (shower%partons(j)%p%nr > &
                           shower%partons(j + 1)%p%nr) then
                         call swap_pointers (shower%partons(j), &
                              shower%partons(j + 1))
                         changed = .true.
                      end if
                   end if
                end if
             end if
          end do INNERDO_T
          if (.not. changed) exit OUTERDO_T
       end do OUTERDO_T
    end if
  end subroutine shower_sort_partons

  module subroutine shower_cleanup (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i
    if (allocated (shower%interactions)) then
    do i = 1, size (shower%interactions)
       if (allocated (shower%interactions(i)%i%partons)) &
            deallocate (shower%interactions(i)%i%partons)
       deallocate (shower%interactions(i)%i)
    end do
    deallocate (shower%interactions)
    end if
    if (allocated (shower%partons))  deallocate (shower%partons)
  end subroutine shower_cleanup

  module function shower_get_next_free_nr (shower) result (next_number)
    class(shower_t), intent(inout) :: shower
    integer :: next_number
    next_number = shower%next_free_nr
    shower%next_free_nr = shower%next_free_nr + 1
  end function shower_get_next_free_nr

  pure module subroutine shower_update_max_color_nr (shower, index)
    class(shower_t), intent(inout) :: shower
    integer, intent(in) :: index
    if (index > shower%next_color_nr) then
       shower%next_color_nr = index
    end if
  end subroutine shower_update_max_color_nr

  module function shower_get_next_color_nr (shower) result (next_color)
    class(shower_t), intent(inout) :: shower
    integer :: next_color
    next_color = shower%next_color_nr
    shower%next_color_nr = shower%next_color_nr + 1
  end function shower_get_next_color_nr

  subroutine shower_enlarge_partons_array (shower, custom_length)
    type(shower_t), intent(inout) :: shower
    integer, intent(in), optional :: custom_length
    integer :: i, length, oldlength
    type(parton_pointer_t), dimension(:), allocatable :: tmp_partons
    if (debug_on) call msg_debug (D_SHOWER, "shower_enlarge_partons_array")
    if (present(custom_length)) then
       length = custom_length
    else
       length = 10
    end if
    if (debug_active (D_SHOWER)) then
       if (length < 1) then
          call msg_bug ("Shower: no parton_pointers added in shower%partons")
       end if
    end if
    if (allocated (shower%partons)) then
       oldlength = size (shower%partons)
       allocate (tmp_partons(1:oldlength))
       do i = 1, oldlength
          tmp_partons(i)%p => shower%partons(i)%p
       end do
       deallocate (shower%partons)
    else
       oldlength = 0
    end if
    allocate (shower%partons(1:oldlength + length))
    do i = 1, oldlength
       shower%partons(i)%p => tmp_partons(i)%p
    end do
    do i = oldlength + 1, oldlength + length
       shower%partons(i)%p => null()
    end do
  end subroutine shower_enlarge_partons_array

  module subroutine shower_add_child (shower, prt, child)
    class(shower_t), intent(inout) :: shower
    type(parton_t), pointer :: prt
    integer, intent(in) :: child
    integer :: i, lastfree
    type(parton_pointer_t) :: newprt
    if (child /= 1 .and. child /= 2) then
       call msg_bug ("Shower: Adding child in nonexisting place")
    end if
    allocate (newprt%p)
    newprt%p%nr = shower%get_next_free_nr ()
    !!! add new parton as child
    if (child == 1) then
       prt%child1 => newprt%p
    else
       prt%child2 => newprt%p
    end if
    newprt%p%parent => prt
    if (associated (prt%settings)) then
       newprt%p%settings => prt%settings
    end if
    newprt%p%interactionnr = prt%interactionnr
    !!! add new parton to shower%partons list
    if (associated (shower%partons (size(shower%partons))%p)) then
       call shower_enlarge_partons_array (shower)
    end if
    !!! find last free pointer and let it point to the new parton
    lastfree = 0
    do i = size (shower%partons), 1, -1
       if (.not. associated (shower%partons(i)%p)) then
          lastfree = i
       end if
    end do
    if (lastfree == 0) then
       call msg_bug ("Shower: no free pointers found")
    end if
    shower%partons(lastfree)%p => newprt%p
  end subroutine shower_add_child

  module subroutine shower_add_parent (shower, prt)
    class(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout), target :: prt
    integer :: i, lastfree
    type(parton_pointer_t) :: newprt
    if (debug_on)  call msg_debug2 &
         (D_SHOWER, "shower_add_parent: for parton nr", prt%nr)
    allocate (newprt%p)
    newprt%p%nr = shower%get_next_free_nr ()
    !!! add new parton as parent
    newprt%p%child1 => prt
    prt%parent => newprt%p
    if (associated (prt%settings)) then
       newprt%p%settings => prt%settings
    end if
    newprt%p%interactionnr = prt%interactionnr
    !!! add new parton to shower%partons list
    if (.not. allocated (shower%partons) .or. &
         associated (shower%partons(size(shower%partons))%p)) then
       call shower_enlarge_partons_array (shower)
    end if
    !!! find last free pointer and let it point to the new parton
    lastfree = 0
    do i = size(shower%partons), 1, -1
       if (.not. associated (shower%partons(i)%p)) then
          lastfree = i
       end if
    end do
    if (debug_active (D_SHOWER)) then
       if (lastfree == 0) then
          call msg_bug ("Shower: no free pointers found")
       end if
    end if
    shower%partons(lastfree)%p => newprt%p
  end subroutine shower_add_parent

  pure function shower_get_total_momentum (shower) result (mom)
    type(shower_t), intent(in) :: shower
    type(vector4_t) :: mom
    integer :: i
    if (.not. allocated (shower%partons)) return
    mom = vector4_null
    do i = 1, size (shower%partons)
       if (.not. associated (shower%partons(i)%p)) cycle
       if (shower%partons(i)%p%is_final ()) then
          mom = mom + shower%partons(i)%p%momentum
       end if
    end do
  end function shower_get_total_momentum

  module function shower_get_nr_of_partons (shower, mine, &
       include_remnants, no_hard_prts, only_colored) result (nr)
    class(shower_t), intent(in) :: shower
    real(default), intent(in), optional :: mine
    logical, intent(in), optional :: include_remnants, no_hard_prts, &
         only_colored
    logical :: no_hard, only_col, include_rem
    integer :: nr, i
    nr = 0
    no_hard = .false.; if (present (no_hard_prts)) &
         no_hard = no_hard_prts
    only_col = .false.; if (present (only_colored)) &
         only_col = only_colored
    include_rem = .true.; if (present (include_remnants)) &
         include_rem = include_remnants
    do i = 1, size (shower%partons)
       if (.not. associated (shower%partons(i)%p)) cycle
       associate (prt => shower%partons(i)%p)
          if (.not. prt%is_final ()) cycle
          if (present (only_colored)) then
             if (only_col) then
                if (.not. prt%is_colored ()) cycle
             else
                if (prt%is_colored ()) cycle
             end if
          end if
          if (no_hard) then
             if (shower%partons(i)%p%belongstointeraction) cycle
          end if
          if (.not. include_rem) then
             if (prt%type == BEAM_REMNANT) cycle
          end if
          if (present(mine)) then
             if (prt%momentum%p(0) < mine) cycle
          end if
          nr = nr + 1
       end associate
    end do
  end function shower_get_nr_of_partons

  function shower_get_nr_of_final_colored_ME_partons (shower) result (nr)
    type(shower_t), intent(in) :: shower
    integer :: nr
    integer :: i, j
    type(parton_t), pointer :: prt
    nr = 0
    do i = 1, size (shower%interactions)
       do j = 1, size (shower%interactions(i)%i%partons)
          prt => shower%interactions(i)%i%partons(j)%p
          if (.not. associated (prt)) cycle
          if (.not. prt%is_colored ()) cycle
          if (prt%belongstointeraction .and. prt%belongstoFSR .and. &
               (prt%type /= INTERNAL)) then
             nr = nr +1
          end if
       end do
    end do
  end function shower_get_nr_of_final_colored_ME_partons

  module subroutine shower_get_final_colored_ME_momenta (shower, momenta)
    class(shower_t), intent(in) :: shower
    type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
    type(parton_pointer_t), dimension(:), allocatable :: partons
    integer :: i, j, index, s
    type(parton_t), pointer :: prt
    s = shower_get_nr_of_final_colored_ME_partons (shower)
    if (s == 0) return
    allocate (partons(1:s))
    allocate (momenta(1:s))
    index = 0
    do i = 1, size (shower%interactions)
       do j = 1, size (shower%interactions(i)%i%partons)
          prt => shower%interactions(i)%i%partons(j)%p
          if (.not. associated (prt)) cycle
          if (.not. prt%is_colored ()) cycle
          if (prt%belongstointeraction .and. prt%belongstoFSR .and. &
               (prt%type /= INTERNAL)) then
             index = index + 1
             partons(index)%p => prt
          end if
       end do
    end do
    do i = 1, s   ! pointers forbid array notation
       momenta(i) = partons(i)%p%momentum
    end do
  end subroutine shower_get_final_colored_ME_momenta

  recursive function interaction_fsr_is_finished_for_parton &
       (prt) result (finished)
    type(parton_t), intent(in) :: prt
    logical :: finished
    if (prt%belongstoFSR) then
       !!! FSR partons
       if (associated (prt%child1)) then
          finished = interaction_fsr_is_finished_for_parton (prt%child1) &
               .and. interaction_fsr_is_finished_for_parton (prt%child2)
       else
          finished = prt%t <= prt%mass_squared ()
       end if
    else
       !!! search for emitted timelike partons in ISR shower
       if (.not. associated (prt%initial)) then
          !!! no inital -> no ISR
          finished = .true.
       else if (.not. associated (prt%parent)) then
          finished = .false.
       else
          if (.not. prt%parent%is_proton ()) then
             if (associated (prt%child2)) then
                finished = interaction_fsr_is_finished_for_parton (prt%parent) .and. &
                     interaction_fsr_is_finished_for_parton (prt%child2)
             else
                finished = interaction_fsr_is_finished_for_parton (prt%parent)
             end if
          else
             if (associated (prt%child2)) then
                finished = interaction_fsr_is_finished_for_parton (prt%child2)
             else
                !!! only second partons can come here -> if that happens FSR
                !!!     evolution is not existing
                finished = .true.
             end if
          end if
       end if
    end if
  end function interaction_fsr_is_finished_for_parton

  function interaction_fsr_is_finished (interaction) result (finished)
    type(shower_interaction_t), intent(in) :: interaction
    logical :: finished
    integer :: i
    finished = .true.
    if (.not. allocated (interaction%partons)) return
    do i = 1, size (interaction%partons)
       if (.not. interaction_fsr_is_finished_for_parton &
            (interaction%partons(i)%p)) then
          finished = .false.
          exit
       end if
    end do
  end function interaction_fsr_is_finished

  module function shower_interaction_get_s (interaction) result (s)
    type(shower_interaction_t), intent(in) :: interaction
    real(default) :: s
    s = (interaction%partons(1)%p%initial%momentum + &
         interaction%partons(2)%p%initial%momentum)**2
  end function shower_interaction_get_s

  function shower_fsr_is_finished (shower) result (finished)
    type(shower_t), intent(in) :: shower
    logical :: finished
    integer :: i
    finished = .true.
    if (.not. allocated (shower%interactions)) return
    do i = 1, size(shower%interactions)
       if (.not. interaction_fsr_is_finished (shower%interactions(i)%i)) then
          finished = .false.
          exit
       end if
    end do
  end function shower_fsr_is_finished

  function shower_isr_is_finished (shower) result (finished)
    type(shower_t), intent(in) :: shower
    logical :: finished
    integer :: i
    type(parton_t), pointer :: prt
    finished = .true.
    if (.not.allocated (shower%partons)) return
    do i = 1, size (shower%partons)
       if (.not. associated (shower%partons(i)%p)) cycle
       prt => shower%partons(i)%p
       if (shower%settings%isr_pt_ordered) then
          if (.not. prt%belongstoFSR .and. .not. prt%simulated &
              .and. prt%scale > zero) then
             finished = .false.
             exit
          end if
       else
          if (.not. prt%belongstoFSR .and. .not. prt%simulated &
              .and. prt%t < zero) then
             finished = .false.
             exit
          end if
       end if
    end do
  end function shower_isr_is_finished

  subroutine interaction_find_partons_nearest_to_hadron &
       (interaction, prt1, prt2, isr_pt_ordered)
    type(shower_interaction_t), intent(in) :: interaction
    type(parton_t), pointer :: prt1, prt2
    logical, intent(in) :: isr_pt_ordered
    prt1 => null ()
    prt2 => null ()
    prt1 => interaction%partons(1)%p
    do
       if (associated (prt1%parent)) then
          if (prt1%parent%is_proton ()) then
             exit
          else if ((.not. isr_pt_ordered .and. .not. prt1%parent%simulated) &
               .or. (isr_pt_ordered .and. .not. prt1%simulated)) then
             exit
          else
             prt1 => prt1%parent
          end if
       else
          exit
       end if
    end do
    prt2 => interaction%partons(2)%p
    do
       if (associated (prt2%parent)) then
          if (prt2%parent%is_proton ()) then
             exit
          else if ((.not. isr_pt_ordered .and. .not. prt2%parent%simulated) &
               .or. (isr_pt_ordered .and. .not. prt2%simulated)) then
             exit
          else
             prt2 => prt2%parent
          end if
       else
          exit
       end if
    end do
  end subroutine interaction_find_partons_nearest_to_hadron

  module subroutine shower_update_beamremnants (shower)
    class(shower_t), intent(inout) :: shower
    type(parton_t), pointer :: hadron, remnant
    integer :: i
    real(default) :: random
    !!! only proton in first interaction !!?
    !!! currently only first beam-remnant will be updated
    do i = 1,2
       if (associated (shower%interactions(1)%i%partons(i)%p%initial)) then
          hadron => shower%interactions(1)%i%partons(i)%p%initial
       else
          cycle
       end if
       remnant => hadron%child2
       if (associated (remnant)) then
          remnant%momentum = hadron%momentum - hadron%child1%momentum
       end if
       !!! generate flavor of the beam-remnant if beam was proton
       if (abs (hadron%type) == PROTON .and. associated (hadron%child1)) then
          if (hadron%child1%is_quark ()) then
             !!! decide if valence (u,d) or sea quark (s,c,b)
             if ((abs (hadron%child1%type) <= 2) .and. &
                  (hadron%type * hadron%child1%type > zero)) then
                !!! valence quark
                if (abs (hadron%child1%type) == 1) then
                   !!! if d then remaining diquark is uu_1
                   remnant%type = sign (UU1, hadron%type)
                else
                   call shower%rng%generate (random)
                   !!! if u then remaining diquark is ud_0 or ud_1
                   if (random < 0.75_default) then
                      remnant%type = sign (UD0, hadron%type)
                   else
                      remnant%type = sign (UD1, hadron%type)
                   end if
                end if
                remnant%c1 = hadron%child1%c2
                remnant%c2 = hadron%child1%c1
             else if ((hadron%type * hadron%child1%type) < zero) then
                !!! antiquark
                if (.not. associated (remnant%child1)) then
                   call shower%add_child (remnant, 1)
                end if
                if (.not. associated (remnant%child2)) then
                   call shower%add_child (remnant, 2)
                end if
                call shower%rng%generate (random)
                if (random < 0.6666_default) then
                   !!! 2/3 into udq + u
                   if (abs (hadron%child1%type) == 1) then
                      remnant%child1%type = sign (NEUTRON, hadron%type)
                   else if (abs (hadron%child1%type) == 2) then
                      remnant%child1%type = sign (PROTON, hadron%type)
                   else if (abs (hadron%child1%type) == 3) then
                      remnant%child1%type = sign (SIGMA0, hadron%type)
                   else if (abs (hadron%child1%type) == 4) then
                      remnant%child1%type = sign (SIGMACPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 5) then
                      remnant%child1%type = sign (SIGMAB0, hadron%type)
                   end if
                   remnant%child2%type = sign (2, hadron%type)
                else
                   !!! 1/3 into uuq + d
                   if (abs (hadron%child1%type) == 1) then
                      remnant%child1%type = sign (PROTON, hadron%type)
                   else if (abs (hadron%child1%type) == 2) then
                      remnant%child1%type = sign (DELTAPLUSPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 3) then
                      remnant%child1%type = sign (SIGMAPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 4) then
                      remnant%child1%type = sign (SIGMACPLUSPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 5) then
                      remnant%child1%type = sign (SIGMABPLUS, hadron%type)
                   end if
                   remnant%child2%type = sign (1, hadron%type)
                end if
                remnant%c1 = hadron%child1%c2
                remnant%c2 = hadron%child1%c1
                remnant%child1%c1 = 0
                remnant%child1%c2 = 0
                remnant%child2%c1 = remnant%c1
                remnant%child2%c2 = remnant%c2
             else
                !!! sea quark
                if (.not. associated (remnant%child1)) then
                   call shower%add_child (remnant, 1)
                end if
                if (.not. associated (remnant%child2)) then
                   call shower%add_child (remnant, 2)
                end if
                call shower%rng%generate (random)
                if (random < 0.5_default) then
                   !!! 1/2 into usbar + ud_0
                   if (abs (hadron%child1%type) == 3) then
                      remnant%child1%type = sign (KPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 4) then
                      remnant%child1%type = sign (D0, hadron%type)
                   else if (abs (hadron%child1%type) == 5) then
                      remnant%child1%type = sign (BPLUS, hadron%type)
                   end if
                   remnant%child2%type = sign (UD0, hadron%type)
                else if (random < 0.6666_default) then
                   !!! 1/6 into usbar + ud_1
                   if (abs (hadron%child1%type) == 3) then
                      remnant%child1%type = sign (KPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 4) then
                      remnant%child1%type = sign (D0, hadron%type)
                   else if (abs (hadron%child1%type) == 5) then
                      remnant%child1%type = sign (BPLUS, hadron%type)
                   end if
                   remnant%child2%type = sign (UD1, hadron%type)
                else
                   !!! 1/3 into dsbar + uu_1
                   if (abs (hadron%child1%type) == 3) then
                      remnant%child1%type = sign (K0, hadron%type)
                   else if (abs (hadron%child1%type) == 4) then
                      remnant%child1%type = sign (DPLUS, hadron%type)
                   else if (abs (hadron%child1%type) == 5) then
                      remnant%child1%type = sign (B0, hadron%type)
                   end if
                   remnant%child2%type = sign (UU1, hadron%type)
                end if
                remnant%c1 = hadron%child1%c2
                remnant%c2 = hadron%child1%c1
                remnant%child1%c1 = 0
                remnant%child1%c2 = 0
                remnant%child2%c1 = remnant%c1
                remnant%child2%c2 = remnant%c2
             end if
          else if (hadron%child1%is_gluon ()) then
             if (.not.associated (remnant%child1)) then
                call shower%add_child (remnant, 1)
             end if
             if (.not.associated (remnant%child2)) then
                call shower%add_child (remnant, 2)
             end if
             call shower%rng%generate (random)
             if (random < 0.5_default) then
                !!! 1/2 into u + ud_0
                remnant%child1%type = sign (2, hadron%type)
                remnant%child2%type = sign (UD0, hadron%type)
             else if (random < 0.6666_default) then
                !!! 1/6 into u + ud_1
                remnant%child1%type = sign (2, hadron%type)
                remnant%child2%type = sign (UD1, hadron%type)
             else
                !!! 1/3 into d + uu_1
                remnant%child1%type = sign (1, hadron%type)
                remnant%child2%type = sign (UU1, hadron%type)
             end if
             remnant%c1 = hadron%child1%c2
             remnant%c2 = hadron%child1%c1
             if (hadron%type > 0) then
                remnant%child1%c1 = remnant%c1
                remnant%child2%c2 = remnant%c2
             else
                remnant%child1%c2 = remnant%c2
                remnant%child2%c1 = remnant%c1
             end if
          end if
          remnant%initial => hadron
          if (associated (remnant%child1)) then
             remnant%child1%initial => hadron
             remnant%child2%initial => hadron
             !!! don't care about on-shellness for now
             remnant%child1%momentum = 0.5_default * remnant%momentum
             remnant%child2%momentum = 0.5_default * remnant%momentum
             !!! but care about on-shellness for baryons
             if (mod (remnant%child1%type, 100) >= 10) then
                !!! check if the third quark is set -> meson or baryon
                remnant%child1%t = remnant%child1%mass_squared ()
                remnant%child1%momentum = [remnant%child1%momentum%p(0), &
                     (remnant%child1%momentum%p(1:3) / &
                      remnant%child1%momentum%p(1:3)**1) * &
                     sqrt (remnant%child1%momentum%p(0)**2 - remnant%child1%t)]
                remnant%child2%momentum = remnant%momentum &
                     - remnant%child1%momentum
             end if
          end if
       end if
    end do
  end subroutine shower_update_beamremnants

  subroutine interaction_apply_lorentztrafo (interaction, L)
    type(shower_interaction_t), intent(inout) :: interaction
    type(lorentz_transformation_t), intent(in) :: L
    type(parton_t), pointer :: prt
    integer :: i
    !!! ISR part
    do i = 1,2
       prt => interaction%partons(i)%p
       !!! loop over ancestors
       MOTHERS: do
          !!! boost parton
          call parton_apply_lorentztrafo (prt, L)
          if (associated (prt%child2)) then
             !!! boost emitted timelike parton (and daughters)
             call parton_apply_lorentztrafo_recursive (prt%child2, L)
          end if
          if (associated (prt%parent)) then
             if (.not. prt%parent%is_proton ()) then
                prt => prt%parent
             else
                exit
             end if
          else
             exit
          end if
       end do MOTHERS
    end do
    !!! FSR part
    if (associated (interaction%partons(3)%p%parent)) then
       !!! pseudo Parton-Shower histora has been generated -> find
       !!!     mother and go on from there recursively
       prt => interaction%partons(3)%p
       do while (associated (prt%parent))
          prt => prt%parent
       end do
       call parton_apply_lorentztrafo_recursive (prt, L)
    else
       do i = 3, size (interaction%partons)
          call parton_apply_lorentztrafo (interaction%partons(i)%p, L)
       end do
    end if
  end subroutine interaction_apply_lorentztrafo

  subroutine shower_apply_lorentztrafo (shower, L)
    type(shower_t), intent(inout) :: shower
    type(lorentz_transformation_t), intent(in) :: L
    integer :: i
    do i = 1, size (shower%interactions)
       call interaction_apply_lorentztrafo (shower%interactions(i)%i, L)
    end do
  end subroutine shower_apply_lorentztrafo

  subroutine interaction_boost_to_CMframe (interaction, isr_pt_ordered)
    type(shower_interaction_t), intent(inout) :: interaction
    logical, intent(in) :: isr_pt_ordered
    type(vector4_t) :: beta
    type(parton_t), pointer :: prt1, prt2
    call interaction_find_partons_nearest_to_hadron &
         (interaction, prt1, prt2, isr_pt_ordered)
    beta = prt1%momentum + prt2%momentum
    beta = beta / beta%p(0)
    if (debug_active (D_SHOWER)) then
       if (beta**2 > one) then
          call msg_error ("Shower: boost to CM frame: beta > 1")
          return
       end if
    end if
    if (space_part(beta)**2 > tiny_13) then
       call interaction_apply_lorentztrafo (interaction, &
            boost(space_part(beta)**1 / &
              sqrt (one - space_part(beta)**2), -direction(beta)))
    end if
  end subroutine interaction_boost_to_CMframe

  module subroutine shower_boost_to_CMframe (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i
    do i = 1, size (shower%interactions)
       call interaction_boost_to_CMframe &
            (shower%interactions(i)%i, shower%settings%isr_pt_ordered)
    end do
    ! TODO: (bcn 2015-03-23) this shouldnt be here !
    call shower%update_beamremnants ()
  end subroutine shower_boost_to_CMframe

  module subroutine shower_boost_to_labframe (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i
    do i = 1, size (shower%interactions)
       call interaction_boost_to_labframe &
            (shower%interactions(i)%i, shower%settings%isr_pt_ordered)
    end do
  end subroutine shower_boost_to_labframe

  subroutine interaction_boost_to_labframe (interaction, isr_pt_ordered)
    type(shower_interaction_t), intent(inout) :: interaction
    logical, intent(in) :: isr_pt_ordered
    type(parton_t), pointer :: prt1, prt2
    type(vector3_t) :: beta
    call interaction_find_partons_nearest_to_hadron &
         (interaction, prt1, prt2, isr_pt_ordered)
    if (.not. associated (prt1%initial) .or. .not. &
         associated (prt2%initial)) then
       return
    end if
    !!! transform partons to overall labframe.
    beta = vector3_canonical(3) * &
         ((prt1%x * prt2%momentum%p(0) - &
           prt2%x * prt1%momentum%p(0)) / &
          (prt1%x * prt2%momentum%p(3) - &
           prt2%x * prt1%momentum%p(3)))
    if (beta**1 > tiny_10) &
         call interaction_apply_lorentztrafo (interaction, &
             boost (beta**1 / sqrt(one - beta**2), -direction(beta)))
  end subroutine interaction_boost_to_labframe

  subroutine interaction_rotate_to_z (interaction, isr_pt_ordered)
    type(shower_interaction_t), intent(inout) :: interaction
    logical, intent(in) :: isr_pt_ordered
    type(parton_t), pointer :: prt1, prt2
    call interaction_find_partons_nearest_to_hadron &
         (interaction, prt1, prt2, isr_pt_ordered)
    if (associated (prt1%initial)) then
       call interaction_apply_lorentztrafo (interaction, &
            rotation_to_2nd (space_part (prt1%momentum), &
            vector3_canonical(3) * sign (one, &
            prt1%initial%momentum%p(3))))
    end if
  end subroutine interaction_rotate_to_z

  module subroutine shower_rotate_to_z (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i
    do i = 1, size (shower%interactions)
       call interaction_rotate_to_z &
            (shower%interactions(i)%i, shower%settings%isr_pt_ordered)
    end do
    ! TODO: (bcn 2015-03-23) this shouldnt be here !
    call shower%update_beamremnants ()
  end subroutine shower_rotate_to_z

  subroutine interaction_generate_primordial_kt &
       (interaction, primordial_kt_width, primordial_kt_cutoff, rng)
    type(shower_interaction_t), intent(inout) :: interaction
    real(default), intent(in) :: primordial_kt_width, primordial_kt_cutoff
    class(rng_t), intent(inout), allocatable :: rng
    type(parton_t), pointer :: had1, had2
    type(vector4_t) :: momenta(2)
    type(vector3_t) :: beta
    real(default) :: pt (2), phi(2)
    real(default) :: shat
    real(default) :: btheta, bphi
    integer :: i
    if (vanishes (primordial_kt_width))  return
    if (.not. associated (interaction%partons(1)%p%initial) .or. &
        .not. associated (interaction%partons(2)%p%initial)) then
       return
    end if
    had1 => interaction%partons(1)%p%initial
    had2 => interaction%partons(2)%p%initial
    !!! copy momenta and energy
    momenta(1) = had1%child1%momentum
    momenta(2) = had2%child1%momentum
    GENERATE_PT_PHI: do i = 1, 2
       !!! generate transverse momentum and phi
       GENERATE_PT: do
          call rng%generate (pt (i))
          pt(i) = primordial_kt_width * sqrt(-log(pt(i)))
          if (pt(i) < primordial_kt_cutoff) exit
       end do GENERATE_PT
       call rng%generate (phi (i))
       phi(i) = twopi * phi(i)
    end do GENERATE_PT_PHI
    !!! adjust momenta
    shat = (momenta(1) + momenta(2))**2
    momenta(1) = [momenta(1)%p(0), &
         pt(1) * cos(phi(1)), &
         pt(1) * sin(phi(1)), &
         momenta(1)%p(3)]
    momenta(2) = [momenta(2)%p(0), &
         pt(2) * cos(phi(2)), &
         pt(2) * sin(phi(2)), &
         momenta(2)%p(3)]
    beta = [momenta(1)%p(1) + momenta(2)%p(1), &
         momenta(1)%p(2) + momenta(2)%p(2), zero] / sqrt(shat)
    momenta(1) = boost (beta**1 / sqrt(one - beta**2), -direction(beta)) &
         * momenta(1)
    bphi = azimuthal_angle (momenta(1))
    btheta = polar_angle (momenta(1))
    call interaction_apply_lorentztrafo (interaction, &
         rotation (cos(bphi), sin(bphi), 3) * rotation(cos(btheta), &
         sin(btheta), 2) * rotation(cos(-bphi), sin(-bphi), 3))
    call interaction_apply_lorentztrafo (interaction, &
         boost (beta**1 / sqrt(one - beta**2), -direction(beta)))
  end subroutine interaction_generate_primordial_kt

  module subroutine shower_generate_primordial_kt (shower)
    class(shower_t), intent(inout) :: shower
    integer :: i
    do i = 1, size (shower%interactions)
       call interaction_generate_primordial_kt (shower%interactions(i)%i, &
            shower%settings%isr_primordial_kt_width, &
            shower%settings%isr_primordial_kt_cutoff, shower%rng)
    end do
    ! TODO: (bcn 2015-03-23) this shouldnt be here !
    call shower%update_beamremnants ()
  end subroutine shower_generate_primordial_kt

  subroutine interaction_write (interaction, unit)
    type(shower_interaction_t), intent(in) :: interaction
    integer, intent(in), optional :: unit
    integer :: i, u
    u = given_output_unit (unit); if (u < 0) return
    if (associated (interaction%partons(1)%p)) then
       if (associated (interaction%partons(1)%p%initial)) &
            call interaction%partons(1)%p%initial%write (u)
    end if
    if (associated (interaction%partons(2)%p)) then
       if (associated (interaction%partons(2)%p%initial)) &
            call interaction%partons(2)%p%initial%write (u)
    end if
    if (allocated (interaction%partons)) then
       do i = 1, size (interaction%partons)
          call interaction%partons(i)%p%write (u)
       end do
    end if
    write (u, "(A)")
  end subroutine interaction_write

  module subroutine shower_write (shower, unit)
    class(shower_t), intent(in) :: shower
    integer, intent(in), optional :: unit
    integer :: i, u
    u = given_output_unit (unit); if (u < 0) return
    write (u, "(1x,A)") "------------------------------"
    write (u, "(1x,A)") "WHIZARD internal parton shower"
    write (u, "(1x,A)") "------------------------------"
    call shower%pdf_data%write (u)
    if (size (shower%interactions) > 0) then
       write (u, "(3x,A)") "Interactions: "
       do i = 1, size (shower%interactions)
          write (u, "(4x,A,I0)") "Interaction number ", i
          if (.not. associated (shower%interactions(i)%i)) then
             call msg_fatal ("Shower: missing interaction in shower")
          end if
          call interaction_write (shower%interactions(i)%i, u)
       end do
    else
       write (u, "(3x,A)") "[no interactions in shower]"
    end if
    write (u, "(A)")
    if (allocated (shower%partons)) then
       write (u, "(5x,A)")  "Partons:"
       do i = 1, size (shower%partons)
          if (associated (shower%partons(i)%p)) then
             call shower%partons(i)%p%write (u)
             if (i < size (shower%partons)) then
                if (associated (shower%partons(i + 1)%p)) then
                   if (shower%partons(i)%p%belongstointeraction .and. &
                        .not. shower%partons(i + 1)%p%belongstointeraction) then
                      call write_separator (u)
                   end if
                end if
             end if
          end if
       end do
    else
       write (u, "(5x,A)")  "[no partons in shower]"
    end if
    write (u, "(4x,A)")  "Total Momentum: "
    call vector4_write (shower_get_total_momentum (shower))
    write (u, "(1x,A,L1)") "ISR finished: ", shower_isr_is_finished (shower)
    write (u, "(1x,A,L1)") "FSR finished: ", shower_fsr_is_finished (shower)
  end subroutine shower_write

  module subroutine shower_combine_with_particle_set (shower, particle_set, &
       model_in, model_hadrons)
    class(shower_t), intent(in) :: shower
    type(particle_set_t), intent(inout) :: particle_set
    class(model_data_t), intent(in), target :: model_in
    class(model_data_t), intent(in), target :: model_hadrons
    type(particle_t), dimension(:), allocatable :: particles
    integer, dimension(:), allocatable :: hard_colored_ids, &
         shower_partons_ids, incoming_ids, outgoing_ids
    class(model_data_t), pointer :: model
    logical, dimension(:), allocatable :: hard_colored_mask
    integer :: n_shower_partons, n_remnants, i, j
    integer :: n_in, n_out, n_beam, n_tot_old
    if (signal_is_pending ()) return
    if (debug_on) call msg_debug (D_SHOWER, "shower_combine_with_particle_set")
    if (debug_on) call msg_debug (D_SHOWER, "Particle set before replacing")
    if (debug_active (D_SHOWER))  &
         call particle_set%write (summary=.true., compressed=.true.)

    n_shower_partons = shower%get_nr_of_partons (only_colored = &
                                     .true., no_hard_prts = .true.)
    n_remnants = shower%get_nr_of_partons (only_colored = .false., &
                                           no_hard_prts = .true.)
    if (n_shower_partons + n_remnants > 0) then
       call particle_set%without_hadronic_remnants &
            (particles, n_tot_old, n_shower_partons + n_remnants)
       call count_and_allocate ()
       call replace_outgoings ()
       call set_hard_colored_as_resonant_parents_for_shower ()
       call add_to_pset (n_tot_old, .true.)
       call add_to_pset (n_tot_old + n_remnants, .false.)
       call particle_set%replace (particles)
    end if

    if (debug_on) call msg_debug (D_SHOWER, 'Particle set after replacing')
    if (debug_active (D_SHOWER))  &
         call particle_set%write (summary=.true., compressed=.true.)

  contains

    subroutine count_and_allocate ()
      n_beam = particle_set%get_n_beam ()
      n_in = particle_set%get_n_in ()
      n_out = particle_set%get_n_out ()
      allocate (hard_colored_mask (size (particles)))
      hard_colored_mask = (particles%get_status () == PRT_INCOMING .or. &
           particles%get_status () == PRT_OUTGOING) .and. &
           particles%is_colored ()
      allocate (hard_colored_ids (count (hard_colored_mask)))
      hard_colored_ids = pack ([(i, i=1, size (particles))], hard_colored_mask)
      allocate (shower_partons_ids (n_shower_partons))
      shower_partons_ids = [(n_tot_old + n_remnants + i, i=1, n_shower_partons)]
      allocate (incoming_ids(n_in))
      incoming_ids = [(n_beam +  i, i=1, n_in)]
      allocate (outgoing_ids(n_out))
      outgoing_ids = [(n_tot_old - n_out  + i, i=1, n_out )]
      if (debug_active (D_SHOWER)) then
          print *, 'n_remnants =    ', n_remnants
          print *, 'n_shower_partons =    ', n_shower_partons
          print *, 'n_tot_old =    ', n_tot_old
          print *, 'n_beam =    ', n_beam
          print *, 'n_in, n_out =    ', n_in, n_out
      end if
    end subroutine count_and_allocate

    subroutine replace_outgoings ()
      do i = 1, size (shower%interactions)
         if (i > 1) then
            call msg_bug ('shower_combine_with_particle_set assumes 1 interaction')
         end if
         associate (interaction => shower%interactions(i)%i)
            do j = 3, size (interaction%partons)
               if (associated (interaction%partons(j)%p)) then
                  call replace_parton_in_particles (j, interaction%partons(j)%p)
               end if
            end do
         end associate
      end do
    end subroutine replace_outgoings

    subroutine replace_parton_in_particles (j, prt)
      integer, intent(in) :: j
      type(parton_t), intent(in) :: prt
      integer :: idx
      if (j <= 2) then
         idx = n_beam + j
      else
         idx = n_tot_old - n_out - n_in + j
      end if
      call particles(idx)%set_momentum (prt%momentum)
    end subroutine replace_parton_in_particles

    subroutine set_hard_colored_as_resonant_parents_for_shower ()
      do i = 1, n_tot_old
         if (hard_colored_mask (i)) then
            if (has_splitted (i)) then
              call particles(i)%add_children (shower_partons_ids)
              if (particles(i)%get_status () == PRT_OUTGOING) then
                 call particles(i)%set_status (PRT_RESONANT)
              end if
            end if
         end if
      end do
    end subroutine set_hard_colored_as_resonant_parents_for_shower

    function has_splitted (i) result (splitted)
      logical :: splitted
      integer, intent(in) :: i
      splitted = .false.
      do j = 1, size (shower%partons)
         if (.not. associated (shower%partons(j)%p)) cycle
         if (particles(i)%flv%get_pdg () == shower%partons(j)%p%type) then
            if (all (nearly_equal (particles(i)%p%p, &
                                   shower%partons(j)%p%momentum%p))) then
               splitted = shower%partons(j)%p%is_branched ()
            end if
         end if
      end do
    end function has_splitted
    subroutine add_to_pset (offset, remnants)
      integer, intent(in) :: offset
      logical, intent(in) :: remnants
      integer :: i, j
      j = offset
      do i = 1, size (shower%partons)
         if (.not. associated (shower%partons(i)%p)) cycle
         associate (prt => shower%partons(i)%p)
            if (.not. prt%is_final () .or. &
                prt%belongstointeraction) cycle
            if (remnants) then
               if (prt%is_colored ()) cycle
            else
               if (.not. (prt%is_colored ())) cycle
            end if
            j = j + 1
            call find_model (model, prt%type, model_in, model_hadrons)
            particles (j) = prt%to_particle (model)
            if (remnants) then
               call particles(j)%set_parents ([prt%initial%nr])
               call particles(prt%initial%nr)%add_child (j)
            else
               call particles(j)%set_parents (hard_colored_ids)
            end if
         end associate
      end do
    end subroutine add_to_pset


  end subroutine shower_combine_with_particle_set

  module subroutine shower_write_lhef (shower, unit)
    class(shower_t), intent(in) :: shower
    integer, intent(in), optional :: unit
    integer :: u
    integer :: i
    integer :: c1, c2
    u = given_output_unit (unit);  if (u < 0)  return
    write(u,'(A)') '<LesHouchesEvents version="1.0">'
    write(u,'(A)') '<-- not a complete lhe file - just one event -->'
    write(u,'(A)') '<event>'
    write(u, *) 2 + shower%get_nr_of_partons (), 1, 1.0, 1.0, 1.0, 1.0
    !!! write incoming partons
    do i = 1, 2
       if (abs (shower%partons(i)%p%type) < 1000) then
          c1 = 0
          c2 = 0
          if (shower%partons(i)%p%is_colored ()) then
             if (shower%partons(i)%p%c1 /= 0) c1 = 500 + shower%partons(i)%p%c1
             if (shower%partons(i)%p%c2 /= 0) c2 = 500 + shower%partons(i)%p%c2
          end if
          write (u,*) shower%partons(i)%p%type, -1, 0, 0, c1, c2, &
               shower%partons(i)%p%momentum%p(1), &
               shower%partons(i)%p%momentum%p(2), &
               shower%partons(i)%p%momentum%p(3), &
               shower%partons(i)%p%momentum%p(0), &
               shower%partons(i)%p%momentum**2, zero, 9.0
       else
          write (u,*) shower%partons(i)%p%type, -9, 0, 0, 0, 0, &
               shower%partons(i)%p%momentum%p(1), &
               shower%partons(i)%p%momentum%p(2), &
               shower%partons(i)%p%momentum%p(3), &
               shower%partons(i)%p%momentum%p(0), &
               shower%partons(i)%p%momentum**2, zero, 9.0
       end if
    end do
    !!! write outgoing partons
    do i = 3, size (shower%partons)
       if (.not. associated (shower%partons(i)%p)) cycle
       if (.not. shower%partons(i)%p%is_final ()) cycle
       c1 = 0
       c2 = 0
       if (shower%partons(i)%p%is_colored ()) then
          if (shower%partons(i)%p%c1 /= 0) c1 = 500 + shower%partons(i)%p%c1
          if (shower%partons(i)%p%c2 /= 0) c2 = 500 + shower%partons(i)%p%c2
       end if
       write (u,*) shower%partons(i)%p%type, 1, 1, 2, c1, c2, &
            shower%partons(i)%p%momentum%p(1), &
            shower%partons(i)%p%momentum%p(2), &
            shower%partons(i)%p%momentum%p(3), &
            shower%partons(i)%p%momentum%p(0), &
            shower%partons(i)%p%momentum**2, zero, 9.0
    end do
    write(u,'(A)') '</event>'
    write(u,'(A)') '</LesHouchesEvents>'
  end subroutine shower_write_lhef

  subroutine shower_replace_parent_by_hadron (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout), target :: prt
    type(parton_t), pointer :: remnant => null()
    if (associated (prt%parent)) then
       call shower_remove_parton_from_partons (shower, prt%parent)
       deallocate (prt%parent)
    end if
    if (.not. associated (prt%initial%child2)) then
       call shower%add_child (prt%initial, 2)
    end if
    prt%parent => prt%initial
    prt%parent%child1 => prt
    ! make other child to be a beam-remnant
    remnant => prt%initial%child2
    remnant%type = BEAM_REMNANT
    remnant%momentum = prt%parent%momentum - prt%momentum
    remnant%x = one - prt%x
    remnant%parent => prt%initial
    remnant%t = zero
  end subroutine shower_replace_parent_by_hadron

  subroutine shower_get_first_ISR_scale_for_parton (shower, prt, tmax)
    type(shower_t), intent(inout), target :: shower
    type(parton_t), intent(inout), target :: prt
    real(default), intent(in), optional :: tmax
    real(default) :: t, tstep, random, integral, temp1
    real(default) :: temprand
    if (present(tmax)) then
       t = max (max (-shower%settings%isr_tscalefactor * prt%momentum%p(0)**2, &
            -abs(tmax)), prt%t)
    else
       t = max (-shower%settings%isr_tscalefactor * prt%momentum%p(0)**2, prt%t)
    end if
    call shower%rng%generate (random)
    random = -twopi * log(random)
    !!! compare Integral and log(random) instead of random and exp(-Integral)
    integral = zero
    call prt%set_simulated (.false.)
    do
       call shower%rng%generate (temprand)
       tstep = max (abs (0.01_default * t) * temprand, 0.1_default * &
            shower%settings%min_virtuality)
       if (t + 0.5_default * tstep > - shower%settings%min_virtuality) then
          prt%t = prt%mass_squared ()
          call prt%set_simulated ()
          exit
       end if
       prt%t = t + 0.5_default * tstep
       temp1 = integral_over_z_simple (prt, (random - integral) / tstep)
       integral = integral + tstep * temp1
       if (integral > random) then
          prt%t = t + 0.5_default * tstep
          exit
       end if
       t = t + tstep
    end do
    if (prt%t > - shower%settings%min_virtuality) then
       call shower_replace_parent_by_hadron (shower, prt)
    end if

  contains

    function integral_over_z_simple (prt, final) result (integral)
      type(parton_t), intent(inout) :: prt
      real(default), intent(in) :: final
      real(default), volatile :: integral

      real(default), parameter :: zstepfactor = one
      real(default), parameter :: zstepmin = 0.0001_default
      real(default) :: z, zstep, minz, maxz
      real(default) :: pdfsum
      integer :: quark, d_nf

      integral = zero
      if (debug2_active (D_SHOWER)) then
         print *, "D: integral_over_z_simple: t = ", prt%t
      end if
      minz = prt%x
      ! maxz = maxzz(shat, s, shower%settings%isr_z_cutoff, shower%settings%isr_minenergy)
      maxz = shower%settings%isr_z_cutoff
      z = minz
      d_nf = shower%settings%max_n_flavors
      !!! TODO -> Adapt zstep to structure of divergencies
      if (prt%child1%is_gluon ()) then
         !!! gluon coming from g->gg
         do
            call shower%rng%generate (temprand)
            zstep = max(zstepmin, temprand * zstepfactor * z * (one - z))
            zstep = min(zstep, maxz - z)
            integral = integral + zstep * (D_alpha_s_isr ((one - &
                 (z + 0.5_default * zstep)) * abs(prt%t), &
                 shower%settings) / (abs(prt%t))) * &
                 P_ggg (z + 0.5_default * zstep) * &
                 shower%get_pdf (prt%initial%type, &
                 prt%x / (z + 0.5_default * zstep), abs(prt%t), GLUON)
            if (integral > final) then
               exit
            end if
            z = z + zstep
            if (z >= maxz) then
               exit
            end if
         end do
         !!! gluon coming from q->qg  ! correctly implemented yet?
         if (integral < final) then
            z = minz
            do
               call shower%rng%generate (temprand)
               zstep = max(zstepmin, temprand * zstepfactor * z * (one - z))
               zstep = min(zstep, maxz - z)
               pdfsum = zero
               do quark = -d_nf, d_nf
                  if (quark == 0) cycle
                  pdfsum = pdfsum + shower%get_pdf (prt%initial%type, &
                       prt%x / (z + 0.5_default * zstep), abs(prt%t), quark)
               end do
               integral = integral + zstep * (D_alpha_s_isr &
                    ((z + 0.5_default * zstep) * abs(prt%t), &
                    shower%settings) / (abs(prt%t))) * &
                    P_qqg (one - (z + 0.5_default * zstep)) * pdfsum
               if (integral > final) then
                  exit
               end if
               z = z + zstep
               if (z >= maxz) then
                  exit
               end if
            end do
         end if
      else if (prt%child1%is_quark ()) then
         !!! quark coming from q->qg
         do
            call shower%rng%generate(temprand)
            zstep = max(zstepmin, temprand * zstepfactor * z * (one - z))
            zstep = min(zstep, maxz - z)
            integral = integral + zstep * (D_alpha_s_isr ((one - &
                 (z + 0.5_default * zstep)) * abs(prt%t), &
                 shower%settings) / (abs(prt%t))) * &
                 P_qqg (z + 0.5_default * zstep) * &
                 shower%get_pdf (prt%initial%type, &
                 prt%x / (z + 0.5_default * zstep), abs(prt%t), prt%type)
            if (integral > final) then
               exit
            end if
            z = z + zstep
            if (z >= maxz) then
               exit
            end if
         end do
         !!! quark coming from g->qqbar
         if (integral < final) then
            z = minz
            do
               call shower%rng%generate (temprand)
               zstep = max(zstepmin, temprand * zstepfactor * z*(one - z))
               zstep = min(zstep, maxz - z)
               integral = integral + zstep * (D_alpha_s_isr &
                    ((one - (z + 0.5_default * zstep)) * abs(prt%t), &
                    shower%settings) / (abs(prt%t))) * &
                    P_gqq (z + 0.5_default * zstep) * &
                    shower%get_pdf (prt%initial%type, &
                    prt%x / (z + 0.5_default * zstep), abs(prt%t), GLUON)
               if (integral > final) then
                  exit
               end if
               z = z + zstep
               if (z >= maxz) then
                  exit
               end if
            end do
         end if

      end if
      integral = integral / shower%get_pdf (prt%initial%type, prt%x, &
           abs(prt%t), prt%type)
    end function integral_over_z_simple

  end subroutine shower_get_first_ISR_scale_for_parton

  subroutine shower_prepare_for_simulate_isr_pt (shower, interaction)
    type(shower_t), intent(inout) :: shower
    type(shower_interaction_t), intent(inout) :: interaction
    real(default) :: s
    s = (interaction%partons(1)%p%momentum + &
         interaction%partons(2)%p%momentum)**2
    interaction%partons(1)%p%scale = shower%settings%isr_tscalefactor * 0.25_default * s
    interaction%partons(2)%p%scale = shower%settings%isr_tscalefactor * 0.25_default * s
  end subroutine shower_prepare_for_simulate_isr_pt

  subroutine shower_prepare_for_simulate_isr_ana_test (shower, prt1, prt2)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout), target :: prt1, prt2
    type(parton_t), pointer :: prt, prta, prtb
    real(default) ::  scale, factor, E
    integer :: i
    if (.not. associated (prt1%initial) .or. .not. associated (prt2%initial)) then
       return
    end if
    scale = - (prt1%momentum + prt2%momentum) ** 2
    call prt1%set_simulated ()
    call prt2%set_simulated ()
    call shower%add_parent (prt1)
    call shower%add_parent (prt2)
    factor = sqrt (energy (prt1%momentum)**2 - scale) / &
         space_part_norm(prt1%momentum)
    prt1%parent%type = prt1%type
    prt1%parent%z = one
    prt1%parent%momentum = prt1%momentum
    prt1%parent%t = scale
    prt1%parent%x = prt1%x
    prt1%parent%initial => prt1%initial
    prt1%parent%belongstoFSR = .false.
    prt1%parent%c1 = prt1%c1
    prt1%parent%c2 = prt1%c2

    prt2%parent%type= prt2%type
    prt2%parent%z = one
    prt2%parent%momentum = prt2%momentum
    prt2%parent%t = scale
    prt2%parent%x = prt2%x
    prt2%parent%initial => prt2%initial
    prt2%parent%belongstoFSR = .false.
    prt2%parent%c1 = prt2%c1
    prt2%parent%c2 = prt2%c2

    call shower_get_first_ISR_scale_for_parton (shower, prt1%parent)
    call shower_get_first_ISR_scale_for_parton (shower, prt2%parent)

    !!! redistribute energy among first partons
    prta => prt1%parent
    prtb => prt2%parent

    E = energy (prt1%momentum + prt2%momentum)
    prta%momentum%p(0) = (E**2 - prtb%t + prta%t) / (two * E)
    prtb%momentum%p(0) = E - prta%momentum%p(0)

    call prt1%parent%set_simulated ()
    call prt2%parent%set_simulated ()
    !!! rescale momenta
    do i = 1, 2
       if (i == 1) then
          prt => prt1%parent
       else
          prt => prt2%parent
       end if
       factor = sqrt (energy (prt%momentum)**2 - prt%t) &
            / space_part_norm (prt%momentum)
       prt%momentum = vector4_moving (energy (prt%momentum), &
            factor * space_part (prt%momentum))
    end do

    if (prt1%parent%t < zero) then
       call shower%add_parent (prt1%parent)
       prt1%parent%parent%momentum = prt1%parent%momentum
       prt1%parent%parent%t = prt1%parent%t
       prt1%parent%parent%x = prt1%parent%x
       prt1%parent%parent%initial => prt1%parent%initial
       prt1%parent%parent%belongstoFSR = .false.
       call shower%add_child (prt1%parent%parent, 2)
    end if

    if (prt2%parent%t < zero) then
       call shower%add_parent (prt2%parent)
       prt2%parent%parent%momentum = prt2%parent%momentum
       prt2%parent%parent%t = prt2%parent%t
       prt2%parent%parent%x = prt2%parent%x
       prt2%parent%parent%initial => prt2%parent%initial
       prt2%parent%parent%belongstoFSR = .false.
       call shower%add_child (prt2%parent%parent, 2)
    end if

  end subroutine shower_prepare_for_simulate_isr_ana_test

  subroutine shower_add_children_of_emitted_timelike_parton (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), pointer :: prt

    if (prt%t > prt%mass_squared () + shower%settings%min_virtuality) then
       if (prt%is_quark ()) then
          !!! q -> qg
          call shower%add_child (prt, 1)
          prt%child1%type = prt%type
          prt%child1%momentum%p(0) = prt%z * prt%momentum%p(0)
          prt%child1%t = prt%t
          call shower%add_child (prt, 2)
          prt%child2%type = GLUON
          prt%child2%momentum%p(0) = (one - prt%z) * prt%momentum%p(0)
          prt%child2%t = prt%t
       else
          if (int (prt%x) > 0) then
             call shower%add_child (prt, 1)
             prt%child1%type = int (prt%x)
             prt%child1%momentum%p(0) = prt%z * prt%momentum%p(0)
             prt%child1%t = prt%t
             call shower%add_child (prt, 2)
             prt%child2%type = -int (prt%x)
             prt%child2%momentum%p(0) = (one - prt%z) * prt%momentum%p(0)
             prt%child2%t= prt%t
          else
             call shower%add_child (prt, 1)
             prt%child1%type = GLUON
             prt%child1%momentum%p(0) = prt%z * prt%momentum%p(0)
             prt%child1%t = prt%t
             call shower%add_child (prt, 2)
             prt%child2%type = GLUON
             prt%child2%momentum%p(0) = (one - prt%z) * prt%momentum%p(0)
             prt%child2%t = prt%t
          end if
       end if
    end if
  end subroutine shower_add_children_of_emitted_timelike_parton

  subroutine shower_simulate_children_ana (shower,prt)
    type(shower_t), intent(inout), target :: shower
    type(parton_t), intent(inout) :: prt
    real(default), dimension(1:2) :: random, integral
    integer, dimension(1:2) :: gtoqq
    integer :: daughter
    type(parton_t), pointer :: daughterprt
    integer :: n_loop

    if (signal_is_pending ()) return
    if (debug2_active (D_SHOWER))  &
         print *, "D: shower_simulate_children_ana: for parton " , prt%nr
    gtoqq = 0
    if (.not. associated (prt%child1) .or. .not. associated (prt%child2)) then
       call msg_error ("Shower: error in simulate_children_ana: no children.")
       return
    end if

    if (HADRON_REMNANT <= abs (prt%type) .and. abs (prt%type) <= HADRON_REMNANT_OCTET) then
       !!! prt is beam-remnant
       call prt%set_simulated ()
       return
    end if

    !!! check if partons are "internal" -> fixed scale
    if (prt%child1%type == INTERNAL) then
       call prt%child1%set_simulated ()
    end if
    if (prt%child2%type == INTERNAL) then
       call prt%child2%set_simulated ()
    end if

    integral = zero

    !!! impose constraints by angular ordering -> cf. (26) of Gaining analytic control
    !!! check if no branchings are possible
    if (.not. prt%child1%simulated) then
       prt%child1%t = min (prt%child1%t, &
            0.5_default * prt%child1%momentum%p(0)**2 * (one - &
            prt%get_costheta ()))
       if (.not. associated (prt%child1%settings)) &
            prt%child1%settings => shower%settings
       if (min (prt%child1%t, prt%child1%momentum%p(0)**2) < &
            prt%child1%mass_squared () + &
            prt%child1%settings%min_virtuality) then
          prt%child1%t = prt%child1%mass_squared ()
          call prt%child1%set_simulated ()
       end if
    end if
    if (.not. prt%child2%simulated) then
       prt%child2%t = min (prt%child2%t, &
            0.5_default * prt%child2%momentum%p(0)**2 * (one - &
            prt%get_costheta ()))
       if (.not. associated (prt%child2%settings)) &
            prt%child2%settings => shower%settings
       if (min (prt%child2%t, prt%child2%momentum%p(0)**2) < &
            prt%child2%mass_squared () + &
            prt%child2%settings%min_virtuality) then
          prt%child2%t = prt%child2%mass_squared ()
          call prt%child2%set_simulated ()
       end if
    end if

    call shower%rng%generate (random)

    n_loop = 0
    do
       if (signal_is_pending ()) return
       n_loop = n_loop + 1
       if (n_loop > 900) then
          !!! try with massless quarks
          treat_duscb_quarks_massless = .true.
       end if
       if (n_loop > 1000) then
          call msg_message ("simulate_children_ana failed for parton ", prt%nr)
          call msg_warning ("too many loops in simulate_children_ana")
          call shower%write ()
          shower%valid = .false.
          return
       end if

       !!! check if a branching in the range t(i) to t(i) - tstep(i) occurs
       if (.not. prt%child1%simulated) then
          call parton_simulate_stept &
               (prt%child1, shower%rng, integral(1), random(1), gtoqq(1))
       end if
       if (.not. prt%child2%simulated) then
          call parton_simulate_stept &
               (prt%child2, shower%rng, integral(2), random(2), gtoqq(2))
       end if

       if (prt%child1%simulated .and. prt%child2%simulated) then
          if (sqrt (prt%t) <= sqrt (prt%child1%t) + sqrt (prt%child2%t)) then
            !!! virtuality : t - m**2 (assuming it's not fixed)
            if (prt%child1%type == INTERNAL .and. prt%child2%type == INTERNAL) then
               call msg_fatal &
                    ("Shower: both partons fixed, but momentum not conserved")
            else if (prt%child1%type == INTERNAL) then
               !!! reset child2
               call prt%child2%set_simulated (.false.)
               prt%child2%t = min (prt%child1%t, (sqrt (prt%t) - &
                    sqrt (prt%child1%t))**2)
               integral(2) = zero
               call shower%rng%generate (random(2))
            else if (prt%child2%type == INTERNAL) then
               ! reset child1
               call prt%child1%set_simulated (.false.)
               prt%child1%t = min (prt%child2%t, (sqrt (prt%t) - &
                    sqrt (prt%child2%t))**2)
               integral(1) = zero
               call shower%rng%generate (random(1))
            else if (prt%child1%t - prt%child1%mass_squared () > &
                 prt%child2%t - prt%child2%mass_squared ()) then
               !!! reset child2
               call prt%child2%set_simulated (.false.)
               prt%child2%t = min (prt%child1%t, (sqrt (prt%t) - &
                    sqrt (prt%child1%t))**2)
               integral(2) = zero
               call shower%rng%generate (random(2))
            else
               !!! reset child1 ! TODO choose child according to their t
               call prt%child1%set_simulated (.false.)
               prt%child1%t = min (prt%child2%t, (sqrt (prt%t) - &
                    sqrt (prt%child2%t))**2)
               integral(1) = zero
               call shower%rng%generate (random(1))
            end if
          else
             exit
          end if
       end if
    end do

    call parton_apply_costheta (prt, shower%rng)

    if (.not. associated (prt%child1%settings)) &
       prt%child1%settings => shower%settings
    if (.not. associated (prt%child2%settings)) &
       prt%child2%settings => shower%settings
    do daughter = 1, 2
       if (signal_is_pending ()) return
       if (daughter == 1) then
          daughterprt => prt%child1
       else
          daughterprt => prt%child2
       end if
       if (daughterprt%t < daughterprt%mass_squared () + &
            daughterprt%settings%min_virtuality) then
          cycle
       end if
       if (.not. (daughterprt%is_quark () .or. daughterprt%is_gluon ())) then
          cycle
       end if
       if (daughterprt%is_quark ()) then
          !!! q -> qg
          call shower%add_child (daughterprt, 1)
          daughterprt%child1%type = daughterprt%type
          daughterprt%child1%momentum%p(0) = daughterprt%z * &
               daughterprt%momentum%p(0)
          daughterprt%child1%t = daughterprt%t
          call shower%add_child (daughterprt, 2)
          daughterprt%child2%type = GLUON
          daughterprt%child2%momentum%p(0) = (one - daughterprt%z) * &
               daughterprt%momentum%p(0)
          daughterprt%child2%t = daughterprt%t
       else if (daughterprt%is_gluon ()) then
          if (gtoqq(daughter) > 0) then
             call shower%add_child (daughterprt, 1)
             daughterprt%child1%type = gtoqq (daughter)
             daughterprt%child1%momentum%p(0) = &
                  daughterprt%z * daughterprt%momentum%p(0)
             daughterprt%child1%t = daughterprt%t
             call shower%add_child (daughterprt, 2)
             daughterprt%child2%type = - gtoqq (daughter)
             daughterprt%child2%momentum%p(0) = (one - &
                  daughterprt%z) * daughterprt%momentum%p(0)
             daughterprt%child2%t = daughterprt%t
          else
             call shower%add_child (daughterprt, 1)
             daughterprt%child1%type = GLUON
             daughterprt%child1%momentum%p(0) = &
                  daughterprt%z * daughterprt%momentum%p(0)
             daughterprt%child1%t = daughterprt%t
             call shower%add_child (daughterprt, 2)
             daughterprt%child2%type = GLUON
             daughterprt%child2%momentum%p(0) = (one - &
                  daughterprt%z) * daughterprt%momentum%p(0)
             daughterprt%child2%t = daughterprt%t
          end if
       end if
    end do
    call shower_parton_update_color_connections (shower, prt)
  end subroutine shower_simulate_children_ana

  subroutine shower_isr_step_pt (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), target, intent(inout) :: prt
    type(parton_t), pointer :: otherprt

    real(default) :: scale, scalestep
    real(default), volatile :: integral
    real(default) :: random, factor
    real(default) :: temprand1, temprand2

    otherprt => shower%find_recoiler (prt)

    scale = prt%scale
    call shower%rng%generate (temprand1)
    call shower%rng%generate (temprand2)
    scalestep = max (abs (scalefactor1 * scale) * temprand1, &
         scalefactor2 * temprand2 * D_Min_scale)
    call shower%rng%generate (random)
    random = - twopi * log(random)
    integral = zero

    if (scale - 0.5_default * scalestep < D_Min_scale) then
       !!! close enough to cut-off scale -> ignore
       prt%scale = zero
       prt%t = prt%mass_squared ()
       call prt%set_simulated ()
    else
       prt%scale = scale - 0.5_default * scalestep
       factor = scalestep * (D_alpha_s_isr (prt%scale, &
            shower%settings) / (prt%scale * &
            shower%get_pdf (prt%initial%type, prt%x, prt%scale, prt%type)))
       integral = integral + factor * integral_over_z_isr_pt &
            (prt, otherprt, (random - integral) / factor)
       if (integral > random) then
          !!! prt%scale set above and prt%z set in integral_over_z_isr_pt
          call prt%set_simulated ()
          prt%t = - prt%scale / (one - prt%z)
       else
          prt%scale = scale - scalestep
       end if
    end if

  contains

    function integral_over_z_isr_pt (prt, otherprt, final) &
         result (integral)
      type(parton_t), intent(inout) :: prt, otherprt
      real(default), intent(in) :: final
      real(default), volatile :: integral
      real(default) :: mbr, r
      real(default) :: zmin, zmax, z, zstep
      integer :: n_bin
      integer, parameter :: n_total_bins = 100
      real(default) :: quarkpdfsum
      real(default) :: temprand
      integer :: quark, d_nf

      quarkpdfsum = zero
      d_nf = shower%settings%max_n_flavors
      if (debug2_active (D_SHOWER)) then
         print *, "D: integral_over_z_isr_pt: for scale = ", prt%scale
      end if

      integral = zero
      mbr = (prt%momentum + otherprt%momentum)**1
      zmin = prt%x
      zmax = min (one - (sqrt (prt%scale) / mbr) * &
           (sqrt(one + 0.25_default * prt%scale / mbr**2) - &
           0.25_default * sqrt(prt%scale) / mbr), shower%settings%isr_z_cutoff)
      zstep = (zmax - zmin) / n_total_bins

      if (debug_active (D_SHOWER)) then
         if (zmin > zmax) then
            call msg_bug(" error in integral_over_z_isr_pt: zmin > zmax ")
            integral = zero
         end if
      end if

      !!! divide the range [zmin:zmax] in n_total_bins
      BINS: do n_bin = 1, n_total_bins
         z = zmin + zstep * (n_bin - 0.5_default)
         !!! z-value in the middle of the bin

         if (prt%is_gluon ()) then
            QUARKS: do quark = -d_nf, d_nf
               if (quark == 0) cycle quarks
               quarkpdfsum = quarkpdfsum + shower%get_pdf &
                    (prt%initial%type, prt%x / z, prt%scale, quark)
            end do QUARKS
            !!! g -> gg or q -> gq
            integral = integral + (zstep / z) * ((P_ggg (z) + &
                 P_ggg (one - z)) * shower%get_pdf (prt%initial%type, &
                 prt%x / z, prt%scale, GLUON) + P_qqg (one - z) * quarkpdfsum)
         else if (prt%is_quark ()) then
            !!! q -> qg or g -> qq
            integral = integral + (zstep / z) * ( P_qqg (z) * &
                 shower%get_pdf (prt%initial%type, prt%x / z, prt%scale, &
                 prt%type) + &
                 P_gqq(z) * shower%get_pdf (prt%initial%type, prt%x / z, &
                 prt%scale, GLUON))
         else
            ! call msg_fatal ("Bug neither quark nor gluon in" &
            !           // " integral_over_z_isr_pt")
         end if
         if (integral > final) then
            prt%z = z
            call shower%rng%generate (temprand)
            !!! decide type of father partons
              if (prt%is_gluon ()) then
                 if (temprand > (P_qqg (one - z) * quarkpdfsum) / &
                      ((P_ggg (z) + P_ggg (one - z)) * shower%get_pdf &
                      (prt%initial%type, prt%x / z, prt%scale, GLUON) &
                      + P_qqg (one - z) * quarkpdfsum)) then
                    !!! gluon => gluon + gluon
                    prt%aux_pt = GLUON
                 else
                    !!! quark => quark + gluon
                    !!! decide which quark flavor the parent is
                    r = temprand * quarkpdfsum
                    WHICH_QUARK: do quark = -d_nf, d_nf
                       if (quark == 0) cycle WHICH_QUARK
                       if (r > quarkpdfsum - shower%get_pdf (prt%initial%type, &
                            prt%x / z, prt%scale, quark)) then
                          prt%aux_pt = quark
                          exit WHICH_QUARK
                       else
                          quarkpdfsum = quarkpdfsum - shower%get_pdf &
                               (prt%initial%type, prt%x / z, prt%scale, quark)
                       end if
                    end do WHICH_QUARK
                 end if

              else if (prt%is_quark ()) then
                 if (temprand > (P_qqg (z) * shower%get_pdf (prt%initial%type, &
                      prt%x / z, prt%scale, prt%type)) / &
                      (P_qqg (z) * shower%get_pdf (prt%initial%type, prt%x / z, &
                      prt%scale, prt%type) + &
                      P_gqq (z) * shower%get_pdf (prt%initial%type, prt%x / z, &
                      prt%scale, GLUON))) then
                    !!! gluon => quark + antiquark
                    prt%aux_pt = GLUON
                 else
                    !!! quark => quark + gluon
                    prt%aux_pt = prt%type
                 end if
              end if
              exit BINS
           end if
        end do BINS
      end function integral_over_z_isr_pt
    end subroutine shower_isr_step_pt

  module function shower_generate_next_isr_branching_veto &
       (shower) result (next_brancher)
    class(shower_t), intent(inout) :: shower
    type(parton_pointer_t) :: next_brancher
    integer :: i
    type(parton_t),  pointer :: prt
    real(default) :: random
    !!! pointers to branchable partons
    type(parton_pointer_t), dimension(:), allocatable :: partons
    integer :: n_partons
    real(default) :: weight
    real(default) :: temp1, temp2, temp3, E3

    if (signal_is_pending ()) return

    if (shower%settings%isr_pt_ordered) then
       next_brancher = shower%generate_next_isr_branching ()
       return
    end if
    next_brancher%p => null()
    !!! branchable partons
    n_partons = 0
    do i = 1,size (shower%partons)
       prt => shower%partons(i)%p
       if (.not. associated (prt)) cycle
       if (prt%belongstoFSR) cycle
       if (prt%is_final ()) cycle
       if (.not. prt%belongstoFSR .and. prt%simulated) cycle
       n_partons = n_partons + 1
    end do
    if (n_partons == 0) then
       return
    end if
    allocate (partons(1:n_partons))
    n_partons = 1
    do i = 1, size (shower%partons)
       prt => shower%partons(i)%p
       if (.not. associated (prt)) cycle
       if (prt%belongstoFSR) cycle
       if (prt%is_final ()) cycle
       if (.not. prt%belongstoFSR .and. prt%simulated) cycle
       partons(n_partons)%p => shower%partons(i)%p
       n_partons = n_partons + 1
    end do
    !!! generate initial trial scales
    do i = 1, size (partons)
       if (signal_is_pending ()) return
       call generate_next_trial_scale (partons(i)%p)
    end do

    do
       !!! search for parton with the highest trial scale
       prt => partons(1)%p
       do i = 1, size (partons)
          if (prt%t >= zero) cycle
          if (abs (partons(i)%p%t) > abs (prt%t)) then
             prt => partons(i)%p
          end if
       end do

       if (prt%t >= zero) then
          next_brancher%p => null()
          exit
       end if
       !!! generate trial z and type of mother prt
       call generate_trial_z_and_typ (prt)

       !!! weight with pdf and alpha_s
       temp1 = (D_alpha_s_isr ((one - prt%z) * abs(prt%t), &
            shower%settings) / sqrt (alphasxpdfmax))
       temp2 = shower%get_xpdf (prt%initial%type, prt%x, prt%t, &
            prt%type) / sqrt (alphasxpdfmax)
       temp3 = shower%get_xpdf (prt%initial%type, prt%child1%x, prt%child1%t, &
            prt%child1%type) / &
            shower%get_xpdf (prt%initial%type, prt%child1%x, prt%t, &
            prt%child1%type)
       ! TODO: (bcn 2015-02-19) ???
       if (temp1 * temp2 * temp3 > one) then
          print *, "weights:", temp1, temp2, temp3
       end if
       weight = (D_alpha_s_isr ((one - prt%z) * abs(prt%t), &
            shower%settings)) * &
            shower%get_xpdf (prt%initial%type, prt%x, prt%t, prt%type) * &
            shower%get_xpdf (prt%initial%type, prt%child1%x, prt%child1%t, &
            prt%child1%type) / &
            shower%get_xpdf &
            (prt%initial%type, prt%child1%x, prt%t, prt%child1%type)
       if (weight > alphasxpdfmax) then
          print *, "Setting alphasxpdfmax from ", alphasxpdfmax, " to ", weight
          alphasxpdfmax = weight
       end if
       weight = weight / alphasxpdfmax
       call shower%rng%generate (random)
       if (weight < random) then
          !!! discard branching
          call generate_next_trial_scale (prt)
          cycle
       end if
       !!! branching accepted so far
       !!! generate emitted parton
       prt%child2%t = abs(prt%t)
       prt%child2%momentum%p(0) = sqrt (abs(prt%t))
       if (shower%settings%isr_only_onshell_emitted_partons) then
          prt%child2%t = prt%child2%mass_squared ()
       else
          call prt%child2%next_t_ana (shower%rng)
       end if

       if (thetabar (prt, shower%find_recoiler (prt), &
            shower%settings%isr_angular_ordered, E3)) then
          prt%momentum%p(0) = E3
          prt%child2%momentum%p(0) = E3 - prt%child1%momentum%p(0)

          !!! found branching
          call prt%generate_ps_ini (shower%rng)
          next_brancher%p => prt
          call prt%set_simulated ()
          exit
       else
          call generate_next_trial_scale (prt)
          cycle
       end if
    end do
    if (.not. associated (next_brancher%p)) then
       !!! no further branching found -> all partons emitted by hadron
       print *, "--all partons emitted by hadrons---"
       do i = 1, size(partons)
          call shower_replace_parent_by_hadron (shower, partons(i)%p%child1)
       end do
    end if
    !!! some bookkeeping
    call shower%sort_partons ()
    ! call shower%boost_to_CMframe ()        ! really necessary?
    ! call shower%rotate_to_z ()             ! really necessary?
  contains

    subroutine generate_next_trial_scale (prt)
      type(parton_t), pointer, intent(inout) :: prt
      real(default) :: random, F
      real(default) :: zmax = 0.99_default !! ??
      call shower%rng%generate (random)
      F = one   !!! TODO
      F = alphasxpdfmax / (two * pi)
      if (prt%child1%is_quark ()) then
         F = F * (integral_over_P_gqq (prt%child1%x, zmax) + &
                  integral_over_P_qqg (prt%child1%x, zmax))
      else if (prt%child1%is_gluon ()) then
         F = F * (integral_over_P_ggg (prt%child1%x, zmax) + &
              two * shower%settings%max_n_flavors * &
              integral_over_P_qqg (one - zmax, one - prt%child1%x))
      else
         call msg_bug("neither quark nor gluon in generate_next_trial_scale")
      end if
      F = F / shower%get_xpdf (prt%child1%initial%type, prt%child1%x, &
           prt%child1%t, prt%child1%type)
      prt%t = prt%t * random**(one / F)
      if (abs (prt%t) - prt%mass_squared () < &
           prt%settings%min_virtuality) then
         prt%t = prt%mass_squared ()
      end if
    end subroutine generate_next_trial_scale

    subroutine generate_trial_z_and_typ (prt)
      type(parton_t), pointer, intent(inout) :: prt
      real(default) :: random
      real(default) :: z, zstep, zmin, integral
      real(default) :: zmax = 0.99_default !! ??
      if (debug_on) call msg_debug (D_SHOWER, "generate_trial_z_and_typ")
      call shower%rng%generate (random)
      integral = zero
      !!! decide which branching a->bc occurs
      if (prt%child1%is_quark ()) then
         if (random < integral_over_P_qqg (prt%child1%x, zmax) / &
              (integral_over_P_qqg (prt%child1%x, zmax) + &
               integral_over_P_gqq (prt%child1%x, zmax))) then
            prt%type = prt%child1%type
            prt%child2%type = GLUON
            integral = integral_over_P_qqg (prt%child1%x, zmax)
         else
            prt%type = GLUON
            prt%child2%type = - prt%child1%type
            integral = integral_over_P_gqq (prt%child1%x, zmax)
         end if
      else if (prt%child1%is_gluon ()) then
         if (random < integral_over_P_ggg (prt%child1%x, zmax) / &
              (integral_over_P_ggg (prt%child1%x, zmax) + two * &
              shower%settings%max_n_flavors * &
              integral_over_P_qqg (one - zmax, &
              one - prt%child1%x))) then
            prt%type = GLUON
            prt%child2%type = GLUON
            integral = integral_over_P_ggg (prt%child1%x, zmax)
         else
            call shower%rng%generate (random)
            prt%type = 1 + floor(random * shower%settings%max_n_flavors)
            call shower%rng%generate (random)
            if (random > 0.5_default) prt%type = - prt%type
            prt%child2%type = prt%type
            integral = integral_over_P_qqg (one - zmax, &
                 one - prt%child1%x)
         end if
      else
         call msg_bug("neither quark nor gluon in generate_next_trial_scale")
      end if
      !!! generate the z-value
      !!! z between prt%child1%x and zmax
      ! prt%z = one - random * (one - prt%child1%x)      ! TODO

      call shower%rng%generate (random)
      zmin = prt%child1%x
      zstep = max(0.1_default, 0.5_default * (zmax - zmin))
      z = zmin
      if (zmin > zmax) then
         print *, " zmin = ", zmin, " zmax = ", zmax
         call msg_fatal ("Shower: zmin greater than zmax")
      end if
      !!! procedure pointers would be helpful here
      if (prt%is_quark () .and. prt%child1%is_quark ()) then
         do
            zstep = min(zstep, 0.5_default * (zmax - z))
            if (abs(zstep) < 0.00001) exit
            if (integral_over_P_qqg (zmin, z) < random * integral) then
               if (integral_over_P_qqg (zmin, min(z + zstep, zmax)) &
                    < random * integral) then
                  z = min (z + zstep, zmax)
                  cycle
               else
                  zstep = zstep * 0.5_default
                  cycle
               end if
            end if
         end do
      else if (prt%is_quark () .and. prt%child1%is_gluon ()) then
         do
            zstep = min(zstep, 0.5_default * (zmax - z))
            if (abs(zstep) < 0.00001) exit
            if (integral_over_P_qqg (zmin, z) < random * integral) then
               if (integral_over_P_qqg (zmin, min(z + zstep, zmax)) &
                    < random * integral) then
                  z = min(z + zstep, zmax)
                  cycle
               else
                  zstep = zstep * 0.5_default
                  cycle
               end if
            end if
         end do
      else if (prt%is_gluon () .and. prt%child1%is_quark ()) then
         do
            zstep = min(zstep, 0.5_default * (zmax - z))
            if (abs (zstep) < 0.00001) exit
            if (integral_over_P_gqq (zmin, z) < random * integral) then
               if (integral_over_P_gqq (zmin, min(z + zstep, zmax)) &
                    < random * integral) then
                  z = min (z + zstep, zmax)
                  cycle
               else
                  zstep = zstep * 0.5_default
                  cycle
               end if
            end if
         end do
      else if (prt%is_gluon () .and. prt%child1%is_gluon ()) then
         do
            zstep = min(zstep, 0.5_default * (zmax - z))
            if (abs (zstep) < 0.00001) exit
            if (integral_over_P_ggg (zmin, z) < random * integral) then
               if (integral_over_P_ggg (zmin, min(z + zstep, zmax)) &
                    < random * integral) then
                  z = min(z + zstep, zmax)
                  cycle
               else
                  zstep = zstep * 0.5_default
                  cycle
               end if
            end if
         end do
      else
      end if
      prt%z = z
      prt%x = prt%child1%x / prt%z
    end subroutine generate_trial_z_and_typ
  end function shower_generate_next_isr_branching_veto

  module function shower_find_recoiler (shower, prt) result(recoiler)
    class(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout), target :: prt
    type(parton_t), pointer :: recoiler
    type(parton_t), pointer :: otherprt1, otherprt2
    integer :: n_int
    otherprt1 => null()
    otherprt2 => null()
    DO_INTERACTIONS: do n_int = 1, size(shower%interactions)
       otherprt1 => shower%interactions(n_int)%i%partons(1)%p
       otherprt2 => shower%interactions(n_int)%i%partons(2)%p
       PARTON1: do
          if (associated (otherprt1%parent)) then
             if (.not. otherprt1%parent%is_proton () .and. &
                  otherprt1%parent%simulated) then
                otherprt1 => otherprt1%parent
                if (associated (otherprt1, prt)) then
                   exit PARTON1
                end if
             else
                exit PARTON1
             end if
          else
             exit PARTON1
          end if
       end do PARTON1
       PARTON2: do
          if (associated (otherprt2%parent)) then
             if (.not. otherprt2%parent%is_proton () .and. &
                  otherprt2%parent%simulated) then
                otherprt2 => otherprt2%parent
                if (associated (otherprt2, prt)) then
                   exit PARTON2
                end if
             else
                exit PARTON2
             end if
          else
             exit PARTON2
          end if
       end do PARTON2

       if (associated (otherprt1, prt) .or. associated (otherprt2, prt)) then
          exit DO_INTERACTIONS
       end if
       if (associated (otherprt1%parent, prt) .or. &
            associated (otherprt2%parent, prt)) then
          exit DO_INTERACTIONS
       end if
    end do DO_INTERACTIONS

    recoiler => null()
    if (associated (otherprt1%parent, prt)) then
       recoiler => otherprt2
    else if (associated (otherprt2%parent, prt)) then
       recoiler => otherprt1
    else if (associated (otherprt1, prt)) then
       recoiler => otherprt2
    else if (associated (otherprt2, prt)) then
       recoiler => otherprt1
    else
       call shower%write ()
       call prt%write ()
       call msg_error ("shower_find_recoiler: no otherparton found")
    end if
  end function shower_find_recoiler

  subroutine shower_isr_step (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), target, intent(inout) :: prt
    type(parton_t), pointer :: otherprt => null()
    real(default) :: t, tstep
    real(default), volatile :: integral
    real(default) :: random
    real(default) :: temprand1, temprand2
    otherprt => shower%find_recoiler (prt)
    ! if (.not. otherprt%child1%belongstointeraction) then
    !    otherprt => otherprt%child1
    ! end if

    if (signal_is_pending ()) return
    t = max(prt%t, prt%child1%t)
    call shower%rng%generate (random)
    ! compare Integral and log(random) instead of random and exp(-Integral)
    random = - twopi * log(random)
    integral = zero
    call shower%rng%generate (temprand1)
    call shower%rng%generate (temprand2)
    tstep = max (abs (0.02_default * t) * temprand1, &
         0.02_default * temprand2 * shower%settings%min_virtuality)
    if (t + 0.5_default * tstep > - shower%settings%min_virtuality) then
       prt%t = prt%mass_squared ()
       call prt%set_simulated ()
    else
       prt%t = t + 0.5_default * tstep
       integral = integral + tstep * &
            integral_over_z_isr (shower, prt, otherprt,(random - integral) / tstep)
       if (integral > random) then
          prt%t = t + 0.5_default * tstep
          prt%x = prt%child1%x / prt%z
          call prt%set_simulated ()
       else
          prt%t = t + tstep
       end if
    end if
  end subroutine shower_isr_step

  function integral_over_z_isr (shower, prt, otherprt, final) result (integral)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout) :: prt, otherprt
    real(default), intent(in) :: final
    !!! !!! !!! volatile argument: gfortran 7 aggressive optimization (#809)
    real(default), volatile :: integral
    real(default) :: minz, maxz, shat,s
    integer :: quark

    !!! calculate shat -> s of parton-parton system
    shat = (otherprt%momentum + prt%child1%momentum)**2
    !!! calculate s -> s of hadron-hadron system
    s = (otherprt%initial%momentum + prt%initial%momentum)**2
    integral = zero
    minz = prt%child1%x
    maxz = maxzz (shat, s, shower%settings%isr_z_cutoff, &
         shower%settings%isr_minenergy)

    !!! for gluon
    if (prt%child1%is_gluon ()) then
       !!! 1: g->gg
       prt%type = GLUON
       prt%child2%type = GLUON
         prt%child2%t = abs(prt%t)
         call integral_over_z_part_isr &
              (shower, prt, otherprt, shat, minz, maxz, integral, final)
         if (integral > final) then
            return
         else
            !!! 2: q->gq
            do quark = - shower%settings%max_n_flavors, &
                    shower%settings%max_n_flavors
               if (quark == 0) cycle
               prt%type = quark
               prt%child2%type = quark
               prt%child2%t = abs(prt%t)
               call integral_over_z_part_isr &
                    (shower, prt, otherprt, shat, minz, maxz, integral, final)
               if (integral > final) then
                  return
               end if
            end do
         end if
      else if (prt%child1%is_quark ()) then
         !!! 1: q->qg
         prt%type = prt%child1%type
         prt%child2%type = GLUON
         prt%child2%t = abs(prt%t)
         call integral_over_z_part_isr &
              (shower, prt,otherprt, shat, minz, maxz, integral, final)
         if (integral > final) then
            return
         else
            !!! 2: g->qqbar
            prt%type = GLUON
            prt%child2%type = -prt%child1%type
            prt%child2%t = abs(prt%t)
            call integral_over_z_part_isr &
                 (shower, prt,otherprt, shat, minz, maxz, integral, final)
         end if
      end if
    end function integral_over_z_isr

  subroutine integral_over_z_part_isr &
       (shower, prt, otherprt, shat ,minz, maxz, retvalue, final)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout) :: prt, otherprt
    real(default), intent(in) :: shat, minz, maxz, final
    real(default), intent(inout) :: retvalue
    real(default) :: z, zstep
    real(default) :: r1,r3,s1,s3
    real(default) :: pdf_divisor
    real(default) :: temprand
    real(default), parameter :: zstepfactor = 0.1_default
    real(default), parameter :: zstepmin = 0.0001_default
    if (debug_on) call msg_debug2 (D_SHOWER, "integral_over_z_part_isr")
    if (signal_is_pending ()) return
    pdf_divisor = shower%get_pdf &
         (prt%initial%type, prt%child1%x, prt%t, prt%child1%type)
    z = minz
    s1 = shat + abs(otherprt%t) + abs(prt%child1%t)
    r1 = sqrt (s1**2 - four * abs(otherprt%t * prt%child1%t))
    ZLOOP: do
       if (signal_is_pending ()) return
       if (z >= maxz) then
          exit
       end if
       call shower%rng%generate (temprand)
       if (prt%child1%is_gluon ()) then
          if (prt%is_gluon ()) then
             !!! g-> gg -> divergencies at z->0 and z->1
             zstep = max(zstepmin, temprand * zstepfactor * z * (one - z))
          else
             !!! q-> gq -> divergencies at z->0
             zstep = max(zstepmin, temprand * zstepfactor * (one - z))
          end if
       else
          if (prt%is_gluon ()) then
             !!! g-> qqbar -> no divergencies
             zstep = max(zstepmin, temprand * zstepfactor)
          else
             !!! q-> qg -> divergencies at z->1
             zstep = max(zstepmin, temprand * zstepfactor * (one - z))
          end if
       end if
       zstep = min(zstep, maxz - z)
       prt%z = z + 0.5_default * zstep
       s3 = shat / prt%z + abs(otherprt%t) + abs(prt%t)
       r3 = sqrt (s3**2 - four * abs(otherprt%t * prt%t))
       !!! TODO: WHY is this if needed?
       if (abs(otherprt%t) > eps0) then
          prt%child2%t = min ((s1 * s3 - r1 * r3) / &
               (two * abs(otherprt%t)) - abs(prt%child1%t) - &
               abs(prt%t), abs(prt%child1%t))
       else
          prt%child2%t = abs(prt%child1%t)
       end if
       do
          prt%child2%momentum%p(0) = sqrt (abs(prt%child2%t))
          if (shower%settings%isr_only_onshell_emitted_partons) then
             prt%child2%t = prt%child2%mass_squared ()
          else
             call prt%child2%next_t_ana (shower%rng)
          end if
          !!! take limits by recoiler into account
          prt%momentum%p(0) = (shat / prt%z + &
               abs(otherprt%t) - abs(prt%child1%t) - &
               prt%child2%t) / (two * sqrt(shat))
          prt%child2%momentum%p(0) = &
               prt%momentum%p(0) - prt%child1%momentum%p(0)
          !!! check if E and t of prt%child2 are consistent
          if (prt%child2%momentum%p(0)**2 < prt%child2%t &
               .and. prt%child2%t > prt%child2%mass_squared ()) then
             !!! E is too small to have p_T^2 = E^2 - t > 0
             !!!      -> cycle to find another solution
             cycle
          else
             !!! E is big enough -> exit
             exit
          end if
       end do
       if (thetabar (prt, otherprt, shower%settings%isr_angular_ordered) &
            .and. pdf_divisor > zero &
            .and. prt%child2%momentum%p(0) > zero) then
          retvalue = retvalue + (zstep / prt%z) * &
               (D_alpha_s_isr ((one - prt%z) * prt%t, &
               shower%settings) * &
               P_prt_to_child1 (prt) * &
               shower%get_pdf (prt%initial%type, prt%child1%x / prt%z, &
               prt%t, prt%type)) / (abs(prt%t) * pdf_divisor)
       end if
       if (retvalue > final) then
          exit
       else
          z = z + zstep
       end if
    end do ZLOOP
  end subroutine integral_over_z_part_isr

  module function shower_generate_next_isr_branching &
       (shower) result (next_brancher)
    class(shower_t), intent(inout) :: shower
    type(parton_pointer_t) :: next_brancher
    integer i, index
    type(parton_t),  pointer :: prt
    next_brancher%p => null()
    do
       if (signal_is_pending ()) return
       if (shower_isr_is_finished (shower)) exit
       !!! find mother with highest |t| or pt to be simulated
       index = 0
       call shower%sort_partons ()
       do i = 1,size (shower%partons)
          prt => shower%partons(i)%p
          if (.not. associated (prt)) cycle
          if (.not. shower%settings%isr_pt_ordered) then
             if (prt%belongstointeraction) cycle
          end if
          if (prt%belongstoFSR) cycle
          if (prt%is_final ()) cycle
          if (.not. prt%belongstoFSR .and. prt%simulated) cycle
          index = i
          exit
       end do
       if (debug_active (D_SHOWER)) then
          if (index == 0) then
             call msg_fatal(" no branchable partons found")
          end if
       end if

       prt => shower%partons(index)%p

       !!! ISR simulation
       if (shower%settings%isr_pt_ordered) then
          call shower_isr_step_pt (shower, prt)
       else
          call shower_isr_step (shower, prt)
       end if
       if (prt%simulated) then
          if (prt%t < zero) then
             next_brancher%p => prt
             if (.not. shower%settings%isr_pt_ordered) &
                  call prt%generate_ps_ini (shower%rng)
             exit
          else
             if (.not. shower%settings%isr_pt_ordered) then
                call shower_replace_parent_by_hadron (shower, prt%child1)
             else
                call shower_replace_parent_by_hadron (shower, prt)
             end if
          end if
       end if
    end do

    !!! some bookkeeping
    call shower%sort_partons ()
    call shower%boost_to_CMframe ()        !!! really necessary?
    call shower%rotate_to_z ()             !!! really necessary?
  end function shower_generate_next_isr_branching

  module subroutine shower_generate_fsr_for_partons_emitted_in_ISR (shower)
    class(shower_t), intent(inout) :: shower
    integer :: n_int, i
    type(parton_t), pointer :: prt
    if (shower%settings%isr_only_onshell_emitted_partons) return
    if (debug_on)  call msg_debug &
         (D_SHOWER, "shower_generate_fsr_for_partons_emitted_in_ISR")
    INTERACTIONS_LOOP: do n_int = 1, size (shower%interactions)
       INCOMING_PARTONS_LOOP: do i = 1, 2
          if (signal_is_pending ()) return
          prt => shower%interactions(n_int)%i%partons(i)%p
          PARENT_PARTONS_LOOP: do
             if (associated (prt%parent)) then
                if (.not. prt%parent%is_proton ()) then
                   prt => prt%parent
                else
                   exit
                end if
             else
                exit
             end if
             if (associated (prt%child2)) then
                if (prt%child2%is_branched ()) then
                   call shower_parton_generate_fsr (shower, prt%child2)
                end if
             else
                ! call msg_fatal ("Shower: no child2 associated?")
             end if
          end do PARENT_PARTONS_LOOP
       end do INCOMING_PARTONS_LOOP
    end do INTERACTIONS_LOOP
  end subroutine shower_generate_fsr_for_partons_emitted_in_ISR

  module subroutine shower_execute_next_isr_branching (shower, prtp)
    class(shower_t), intent(inout) :: shower
    type(parton_pointer_t), intent(inout) :: prtp
    type(parton_t), pointer :: prt, otherprt
    type(parton_t), pointer :: prta, prtb, prtc, prtr
    real(default) :: mbr
    real(default) :: phirand
    if (debug_on) call msg_debug (D_SHOWER, "shower_execute_next_isr_branching")
    if (.not. associated (prtp%p)) then
       call msg_fatal ("Shower: prtp not associated")
    end if

    prt => prtp%p

    if ((.not. shower%settings%isr_pt_ordered .and. &
         prt%t > - shower%settings%min_virtuality) .or. &
         (shower%settings%isr_pt_ordered .and. prt%scale < D_Min_scale)) then
       call msg_error ("Shower: no branching to be executed.")
    end if

    otherprt => shower%find_recoiler (prt)
    if (shower%settings%isr_pt_ordered) then
       !!! get the recoiler
       otherprt => shower%find_recoiler (prt)
       if (associated (otherprt%parent)) then
          !!! Why only for pt ordered
          if (.not. otherprt%parent%is_proton () .and. &
               shower%settings%isr_pt_ordered) otherprt => otherprt%parent
       end if
       if (.not. associated (prt%parent)) then
          call shower%add_parent (prt)
       end if
       prt%parent%belongstoFSR = .false.
       if (.not. associated (prt%parent%child2)) then
          call shower%add_child (prt%parent, 2)
       end if

       prta => prt%parent          !!! new parton a with branching a->bc
       prtb => prt                 !!! former parton
       prtc => prt%parent%child2   !!! emitted parton
       prtr => otherprt            !!! recoiler

       mbr = (prtb%momentum + prtr%momentum)**1

       !!! 1. assume you are in the restframe
       !!! 2. rotate by random phi
       call shower%rng%generate (phirand)
       phirand = twopi * phirand
       call shower_apply_lorentztrafo (shower, &
            rotation(cos(phirand), sin(phirand),vector3_canonical(3)))
       !!! 3. Put the b off-shell
       !!! and
       !!! 4. construct the massless a
       !!! and the parton (eventually emitted by a)

       !!! generate the flavor of the parent (prta)
       if (prtb%aux_pt /= 0) prta%type = prtb%aux_pt
       if (prtb%is_quark ()) then
          if (prta%type == prtb%type) then
             !!! (anti)-quark -> (anti-)quark + gluon
             prta%type = prtb%type   ! quarks have same flavor
             prtc%type = GLUON        ! emitted gluon
          else
             !!! gluon -> quark + antiquark
             prta%type = GLUON
             prtc%type = - prtb%type
          end if
       else if (prtb%is_gluon ()) then
          prta%type = GLUON
          prtc%type = GLUON
       else
          ! STOP "Bug in shower_execute_next_branching: neither quark nor gluon"
       end if

       prta%initial => prtb%initial
       prta%belongstoFSR = .false.
       prta%scale = prtb%scale
       prta%x = prtb%x / prtb%z

       prtb%momentum = vector4_moving ((mbr**2 + prtb%t) / (two * mbr), &
            vector3_canonical(3) * &
            sign ((mbr**2 - prtb%t) / (two * mbr), &
            prtb%momentum%p(3)))
       prtr%momentum = vector4_moving ((mbr**2 - prtb%t) / (two * mbr), &
            vector3_canonical(3) * &
            sign( (mbr**2 - prtb%t) / (two * mbr), &
            prtr%momentum%p(3)))

       prta%momentum = vector4_moving ((0.5_default / mbr) * &
            ((mbr**2 / prtb%z) + prtb%t - prtc%mass_squared ()), &
            vector3_null)
       prta%momentum = vector4_moving (prta%momentum%p(0), &
            vector3_canonical(3) * &
            (0.5_default / prtb%momentum%p(3)) * &
            ((mbr**2 / prtb%z) - two &
            * prtr%momentum%p(0) * prta%momentum%p(0) ) )
       if (prta%momentum%p(0)**2 - prta%momentum%p(3)**2 - &
            prtc%mass_squared () > zero) then
          !!! This SHOULD be always fulfilled???
          prta%momentum = vector4_moving (prta%momentum%p(0), &
               vector3_moving([sqrt (prta%momentum%p(0)**2 - &
               prta%momentum%p(3)**2 - &
               prtc%mass_squared ()), zero, &
               prta%momentum%p(3)]))
       end if
       prtc%momentum = prta%momentum - prtb%momentum

       !!! 5. rotate to have a along z-axis
       call shower%boost_to_CMframe ()
       call shower%rotate_to_z ()
       !!! 6. rotate back in phi
       call shower_apply_lorentztrafo (shower, rotation &
            (cos(-phirand), sin(-phirand), vector3_canonical(3)))
    else
       if (prt%child2%t > prt%child2%mass_squared ()) then
          call shower_add_children_of_emitted_timelike_parton &
               (shower, prt%child2)
          call prt%child2%set_simulated ()
       end if

       call shower%add_parent (prt)
       call shower%add_child (prt%parent, 2)

       prt%parent%momentum = prt%momentum
       prt%parent%t = prt%t
       prt%parent%x = prt%x
       prt%parent%initial => prt%initial
       prt%parent%belongstoFSR = .false.

       prta => prt
       prtb => prt%child1
       prtc => prt%child2
    end if
    if (signal_is_pending ()) return
    if (shower%settings%isr_pt_ordered) then
       call prt%parent%generate_ps_ini (shower%rng)
    else
       call prt%generate_ps_ini (shower%rng)
    end if

    !!! add color connections
    if (prtb%is_quark ()) then

       if (prta%type == prtb%type) then
          if (prtb%type > 0) then
             !!! quark -> quark + gluon
             prtc%c2 = prtb%c1
             prtc%c1 = shower%get_next_color_nr ()
             prta%c1 = prtc%c1
          else
             !!! antiquark -> antiquark + gluon
             prtc%c1 = prtb%c2
             prtc%c2 = shower%get_next_color_nr ()
             prta%c2 = prtc%c2
          end if
       else
          !!! gluon -> quark + antiquark
          if (prtb%type > 0) then
             !!! gluon -> quark + antiquark
             prta%c1 = prtb%c1
             prtc%c1 = 0
             prtc%c2 = shower%get_next_color_nr ()
             prta%c2 = prtc%c2
          else
             !!! gluon -> antiquark + quark
             prta%c2 = prtb%c2
             prtc%c1 = shower%get_next_color_nr ()
             prtc%c2 = 0
             prta%c1 = prtc%c1
            end if
         end if
      else if (prtb%is_gluon ()) then
         if (prta%is_gluon ()) then
            !!! g -> gg
            prtc%c2 = prtb%c1
            prtc%c1 = shower%get_next_color_nr ()
            prta%c1 = prtc%c1
            prta%c2 = prtb%c2
         else if (prta%is_quark ()) then
            if (prta%type > 0) then
               prta%c1 = prtb%c1
               prta%c2 = 0
               prtc%c1 = prtb%c2
               prtc%c2 = 0
            else
               prta%c1 = 0
               prta%c2 = prtb%c2
               prtc%c1 = 0
               prtc%c2 = prtb%c1
            end if
         end if
      end if

      call shower%sort_partons ()
      call shower%boost_to_CMframe ()
      call shower%rotate_to_z ()

    end subroutine shower_execute_next_isr_branching

  subroutine shower_remove_parents_and_stuff (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout), target :: prt
    type(parton_t), pointer :: actprt, nextprt
    nextprt => prt%parent
    actprt => null()
    !!! remove children of emitted timelike parton
    if (associated (prt%child2)) then
       if (associated (prt%child2%child1)) then
          call shower_remove_parton_from_partons_recursive &
               (shower, prt%child2%child1)
       end if
       prt%child2%child1 => null()
       if (associated (prt%child2%child2)) then
          call shower_remove_parton_from_partons_recursive &
               (shower, prt%child2%child2)
       end if
       prt%child2%child2 => null()
    end if
    do
       actprt => nextprt
       if (.not. associated (actprt)) then
          exit
       else if (actprt%is_proton ()) then
          !!! remove beam-remnant
          call shower_remove_parton_from_partons (shower, actprt%child2)
          exit
       end if
       if (associated (actprt%parent)) then
          nextprt => actprt%parent
       else
          nextprt => null()
       end if
       call shower_remove_parton_from_partons_recursive &
            (shower, actprt%child2)
       call shower_remove_parton_from_partons (shower, actprt)

    end do
    prt%parent=>null()

  end subroutine shower_remove_parents_and_stuff

  module function shower_get_ISR_scale (shower) result (scale)
    class(shower_t), intent(in) :: shower
    real(default) :: scale
    type(parton_t), pointer :: prt1, prt2
    integer :: i
    scale = zero
    do i = 1, size (shower%interactions)
       call interaction_find_partons_nearest_to_hadron &
            (shower%interactions(i)%i, prt1, prt2, &
            shower%settings%isr_pt_ordered)
       if (.not. prt1%simulated .and. abs(prt1%scale) > scale) &
            scale = abs(prt1%scale)
       if (.not. prt1%simulated .and. abs(prt2%scale) > scale) &
            scale = abs(prt2%scale)
    end do
  end function shower_get_ISR_scale

  module subroutine shower_set_max_isr_scale (shower, newscale)
    class(shower_t), intent(inout) :: shower
    real(default), intent(in) :: newscale
    real(default) :: scale
    type(parton_t), pointer :: prt
    integer :: i,j
    if (debug_on) call msg_debug (D_SHOWER, "shower_set_max_isr_scale: newscale", &
         newscale)
    if (shower%settings%isr_pt_ordered) then
       scale = newscale
    else
       scale = - abs (newscale)
    end if

    INTERACTIONS: do i = 1, size (shower%interactions)
       PARTONS: do j = 1, 2
          prt => shower%interactions(i)%i%partons(j)%p
          do
             if (.not. shower%settings%isr_pt_ordered) then
                if (prt%belongstointeraction) prt => prt%parent
             end if
               if (prt%t < scale) then
                  if (associated (prt%parent)) then
                     prt => prt%parent
                  else
                     exit   !!! unresolved prt found
                  end if
               else
                  exit   !!! prt with scale above newscale found
               end if
            end do
            if (.not. shower%settings%isr_pt_ordered) then
               if (prt%child1%belongstointeraction .or. &
                    prt%is_proton ()) then
                  !!! don't reset scales of "first" spacelike partons
                  !!!    in virtuality ordered shower or hadrons
                  cycle
               end if
            else
               if (prt%is_proton ()) then
                  !!! don't reset scales of hadrons
                  cycle
               end if
            end if
            if (shower%settings%isr_pt_ordered) then
               prt%scale = scale
            else
               prt%t = scale
            end if
            call prt%set_simulated (.false.)
            call shower_remove_parents_and_stuff (shower, prt)
         end do PARTONS
      end do INTERACTIONS
    end subroutine shower_set_max_isr_scale

  module subroutine shower_interaction_generate_fsr_2ton (shower, interaction)
    class(shower_t), intent(inout) :: shower
    type(shower_interaction_t), intent(inout) :: interaction
    type(parton_t), pointer :: prt
    prt => interaction%partons(3)%p
    do
       if (.not. associated (prt%parent)) exit
       prt => prt%parent
    end do
    call shower_parton_generate_fsr (shower, prt)
    call shower_parton_update_color_connections (shower, prt)
  end subroutine shower_interaction_generate_fsr_2ton

  subroutine shower_parton_generate_fsr (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout), target :: prt
    type(parton_pointer_t), dimension(:), allocatable :: partons
    logical :: single_emission = .false.
    if (debug_on) call msg_debug (D_SHOWER, "shower_parton_generate_fsr")
    if (signal_is_pending ()) return
    if (debug_active (D_SHOWER)) then
       if (.not. prt%is_branched ()) then
          call msg_error ("shower_parton_generate_fsr: parton not branched")
          return
       end if
       if (prt%child1%simulated .or. &
           prt%child2%simulated) then
          print *, "children already simulated for parton ", prt%nr
          return
       end if
    end if
    allocate (partons(1))
    partons(1)%p => prt
    if (single_emission) then
       call shower%parton_pointer_array_generate_fsr (partons, partons)
    else
       call shower%parton_pointer_array_generate_fsr_recursive (partons)
    end if
  end subroutine shower_parton_generate_fsr

  recursive module subroutine shower_parton_pointer_array_generate_fsr_recursive &
       (shower, partons)
    class(shower_t), intent(inout) :: shower
    type(parton_pointer_t), dimension(:), allocatable, intent(inout) :: &
         partons
    type(parton_pointer_t), dimension(:), allocatable :: partons_new
    if (debug_on) call msg_debug (D_SHOWER, "shower_parton_pointer_array_generate_fsr_recursive")
    if (signal_is_pending ()) return
    if (size (partons) == 0) return
    call shower%parton_pointer_array_generate_fsr (partons, partons_new)
    call shower%parton_pointer_array_generate_fsr_recursive (partons_new)
  end subroutine shower_parton_pointer_array_generate_fsr_recursive

  module subroutine shower_parton_pointer_array_generate_fsr &
       (shower, partons, partons_new)
    class(shower_t), intent(inout) :: shower
    type(parton_pointer_t), dimension(:), allocatable, intent(inout) :: &
         partons
    type(parton_pointer_t), dimension(:), allocatable, intent(out) :: &
         partons_new
    integer :: i, size_partons, size_partons_new
    if (debug_on) call msg_debug (D_SHOWER, "shower_parton_pointer_array_generate_fsr")
    !!! Simulate highest/first parton
    call shower_simulate_children_ana (shower, partons(1)%p)
    !!! check for new daughters to be included in new_partons
    size_partons = size (partons)
    size_partons_new = size_partons - 1   !!! partons(1) not needed anymore
    if (partons(1)%p%child1%is_branched ()) &
         size_partons_new = size_partons_new + 1
    if (partons(1)%p%child2%is_branched ()) &
         size_partons_new = size_partons_new + 1

    allocate (partons_new (1:size_partons_new))

    if (size_partons > 1) then
       do i = 2, size_partons
          partons_new (i - 1)%p => partons(i)%p
       end do
    end if
    if (partons(1)%p%child1%is_branched ()) &
         partons_new (size_partons)%p => partons(1)%p%child1
    if (partons(1)%p%child2%is_branched ()) then
    !!! check if child1 is already included
       if (size_partons_new == size_partons) then
          partons_new (size_partons)%p => partons(1)%p%child2
       else if (size_partons_new == size_partons + 1) then
          partons_new (size_partons + 1)%p => partons(1)%p%child2
       else
          call msg_fatal ("Shower: wrong sizes in" &
               // "shower_parton_pointer_array_generate_fsr")
       end if
    end if
    deallocate (partons)

  end subroutine shower_parton_pointer_array_generate_fsr

  recursive subroutine shower_parton_update_color_connections &
       (shower, prt)
    type(shower_t), intent(inout) :: shower
    type(parton_t), intent(inout) :: prt
    real(default) :: temprand
    if (.not. associated (prt%child1) .or. &
        .not. associated (prt%child2)) return

    if (signal_is_pending ()) return
    if (prt%is_gluon ()) then
       if (prt%child1%is_quark ()) then
          !!! give the quark the colorpartner and the antiquark
          !!!     the anticolorpartner
          if (prt%child1%type > 0) then
             !!! child1 is quark, child2 is antiquark
             prt%child1%c1 = prt%c1
             prt%child2%c2 = prt%c2
          else
             !!! child1 is antiquark, child2 is quark
             prt%child1%c2 = prt%c2
             prt%child2%c1 = prt%c1
          end if
       else
          !!! g -> gg splitting -> random choosing of partners
          call shower%rng%generate (temprand)
          if (temprand > 0.5_default) then
             prt%child1%c1 = prt%c1
             prt%child1%c2 = shower%get_next_color_nr ()
             prt%child2%c1 = prt%child1%c2
             prt%child2%c2 = prt%c2
          else
             prt%child1%c2 = prt%c2
             prt%child2%c1 = prt%c1
             prt%child2%c2 = shower%get_next_color_nr ()
             prt%child1%c1 = prt%child2%c2
          end if
       end if
    else if (prt%is_quark ()) then
       if (prt%child1%is_quark ()) then
          if (prt%child1%type > 0) then
             !!! q -> q + g
             prt%child2%c1 = prt%c1
             prt%child2%c2 = shower%get_next_color_nr ()
             prt%child1%c1 = prt%child2%c2
          else
             !!! qbar -> qbar + g
             prt%child2%c2 = prt%c2
             prt%child2%c1 = shower%get_next_color_nr ()
             prt%child1%c2 = prt%child2%c1
          end if
       else
          if (prt%child2%type > 0) then
             !!! q -> g + q
             prt%child1%c1 = prt%c1
             prt%child1%c2 = shower%get_next_color_nr ()
             prt%child2%c1 = prt%child1%c2
          else
             !!! qbar -> g + qbar
             prt%child1%c2 = prt%c2
             prt%child1%c1 = shower%get_next_color_nr ()
             prt%child2%c2 = prt%child1%c1
          end if
       end if
    end if

    call shower_parton_update_color_connections (shower, prt%child1)
    call shower_parton_update_color_connections (shower, prt%child2)
  end subroutine shower_parton_update_color_connections

  module function shower_get_pdf (shower, mother, x, Q2, daughter) result (pdf)
    class(shower_t), intent(inout), target :: shower
    integer, intent(in) :: mother, daughter
    real(default), intent(in) :: x, Q2
    real(default) :: pdf
    real(double), save :: f(-6:6) = 0._double
    real(double), save :: lastx, lastQ2 = 0._double
    pdf = zero
    if (debug_active (D_SHOWER)) then
       if (abs (mother) /= PROTON) then
          if (debug_on) call msg_debug (D_SHOWER, "mother", mother)
          call msg_fatal ("Shower: pdf only implemented for (anti-)proton")
       end if
       if (.not. (abs (daughter) >= 1 .and. abs (daughter) <= 6 .or. &
                  daughter == GLUON)) then
          if (debug_on) call msg_debug (D_SHOWER, "daughter", daughter)
          call msg_fatal ("Shower: error in pdf, unknown daughter")
       end if
    end if
    if (x > zero .and. x < one) then
       if ((dble(Q2) - lastQ2) > eps0 .or. (dble(x) - lastx) > eps0) then
          call shower%pdf_data%evolve &
               (dble(x), sqrt (abs (dble(Q2))), f)
       end if
       if (abs (daughter) >= 1 .and. abs (daughter) <= 6) then
          pdf = max (f(daughter * sign (1,mother)), tiny_10)
       else
          pdf = max (f(0), tiny_10)
       end if
    end if
    lastQ2 = dble(Q2)
    lastx  = dble(x)
    if (x > eps0) then
       pdf = pdf / x
    end if
  end function shower_get_pdf

  module function shower_get_xpdf (shower, mother, x, Q2, daughter) result (pdf)
    class(shower_t), intent(inout), target :: shower
    integer, intent(in) :: mother, daughter
    real(default), intent(in) :: x, Q2
    real(default) :: pdf
    real(double), save :: f(-6:6) = 0._double
    real(double), save :: lastx, lastQ2 = 0._double
    pdf = zero
    if (debug_active (D_SHOWER)) then
       if (abs (mother) /= PROTON) then
          if (debug_on) call msg_debug (D_SHOWER, "mother", mother)
          call msg_fatal ("Shower: pdf only implemented for (anti-)proton")
       end if
       if (.not. (abs (daughter) >= 1 .and. abs (daughter) <= 6 .or. &
                  daughter == GLUON)) then
          if (debug_on) call msg_debug (D_SHOWER, "daughter", daughter)
          call msg_fatal ("Shower: error in pdf, unknown daughter")
       end if
    end if
    if (x > zero .and. x < one) then
       if ((dble(Q2) - lastQ2) > eps0 .or. (dble(x) - lastx) > eps0) then
          call shower%pdf_data%evolve &
               (dble(x), sqrt (abs (dble(Q2))), f)
       end if
       if (abs (daughter) >= 1 .and. abs (daughter) <= 6) then
          pdf = max (f(daughter * sign (1,mother)), tiny_10)
       else
          pdf = max (f(0), tiny_10)
       end if
    end if
    lastQ2 = dble(Q2)
    lastx  = dble(x)
  end function shower_get_xpdf

  module subroutine shower_converttopythia (shower)
    class(shower_t), intent(in) :: shower
    integer :: N, NPAD, K
    real(double) :: P, V
    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    SAVE /PYJETS/
    type(parton_t), pointer :: pp, ppparent
    integer :: i
    K = 0
    do i = 1, 2
       !!! get history of the event
       pp => shower%interactions(1)%i%partons(i)%p
       !!! add these partons to the event record
       if (associated (pp%initial)) then
          !!! add hadrons
          K(i,1) = 21
          K(i,2) = pp%initial%type
          K(i,3) = 0
          P(i,1:5) = pp%initial%momentum_to_pythia6 ()
          !!! add partons emitted by the hadron
          ppparent => pp
          do while (associated (ppparent%parent))
             if (ppparent%parent%is_proton ()) then
                exit
             else
                ppparent => ppparent%parent
             end if
          end do
          K(i+2,1) = 21
          K(i+2,2) = ppparent%type
          K(i+2,3) = i
          P(i+2,1:5) = ppparent%momentum_to_pythia6 ()
          !!! add partons in the initial state of the ME
          K(i+4,1) = 21
          K(i+4,2) = pp%type
          K(i+4,3) = i
          P(i+4,1:5) = pp%momentum_to_pythia6 ()
       else
          !!! for e+e- without ISR all entries are the same
          K(i,1) = 21
          K(i,2) = pp%type
          K(i,3) = 0
          P(i,1:5) = pp%momentum_to_pythia6 ()
          P(i+2,:) = P(1,:)
          K(i+2,:) = K(1,:)
          K(i+2,3) = i
          P(i+4,:) = P(1,:)
          K(i+4,:) = K(1,:)
          K(i+4,3) = i
          P(i+4,5) = 0.
       end if
    end do
    N = 6
    !!! create intermediate (fake) Z-Boson
    !K(7,1) = 21
    !K(7,2) = 23
    !K(7,3) = 0
    !P(7,1:4) = P(5,1:4) + P(6,1:4)
    !P(7,5) = P(7,4)**2 - P(7,3)**2 - P(7,2)**2 - P(7,1)**2
    !N = 7
    !!! include partons in the final state of the hard matrix element
    do i = 1, size (shower%interactions(1)%i%partons) - 2
       !!! get partons that are in the final state of the hard matrix element
       pp => shower%interactions(1)%i%partons(2+i)%p
       !!! add these partons to the event record
       K(7+I,1) = 21
       K(7+I,2) = pp%type
       K(7+I,3) = 7
       P(7+I,1:5) = pp%momentum_to_pythia6 ()
       !N = 7 + I
       N = 6 + I
    end do
    !!! include "Z" (again)
    !N = N + 1
    !K(N,1) = 11
    !K(N,2) = 23
    !K(N,3) = 7
    !P(N,1:5) = P(7,1:5)
    !nz = N
    !!! include partons from the final state of the parton shower
    call shower_transfer_final_partons_to_pythia (shower, 8)
    !!! set "children" of "Z"
    !K(nz,4) = 11
    !K(nz,5) = N

    !!! be sure to remove the next partons (=first obsolete partons)
    !!! otherwise they might be interpreted as thrust information
    K(N+1:N+3,1:3) = 0
  end subroutine shower_converttopythia

  subroutine shower_transfer_final_partons_to_pythia (shower, first)
    integer :: N, NPAD, K
    real(double) :: P, V
    COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
    SAVE /PYJETS/
    type(shower_t), intent(in) :: shower
    integer, intent(in) :: first
    type(parton_t), pointer :: prt
    integer :: i, j, n_finals
    type(parton_t), dimension(:), allocatable :: final_partons
    type(parton_t) :: temp_parton
    integer :: minindex, maxindex

    prt => null()

    !!! get total number of final partons
    n_finals = 0
    do i = 1, size (shower%partons)
       if (.not. associated (shower%partons(i)%p)) cycle
       prt => shower%partons(i)%p
       if (.not. prt%belongstoFSR) cycle
       if (associated (prt%child1)) cycle
       n_finals = n_finals + 1
    end do

    allocate (final_partons(1:n_finals))
    j = 1
    do i = 1, size (shower%partons)
       if (.not. associated (shower%partons(i)%p)) cycle
       prt => shower%partons(i)%p
       if (.not. prt%belongstoFSR) cycle
       if (associated (prt%child1)) cycle
       final_partons(j) = shower%partons(i)%p
       j = j + 1
    end do

    !!! move quark to front as beginning of color string
    minindex = 1
    maxindex = size (final_partons)
    FIND_Q: do i = minindex, maxindex
       if (final_partons(i)%type >= 1 .and. final_partons(i)%type <= 6) then
          temp_parton = final_partons(minindex)
          final_partons(minindex) = final_partons(i)
          final_partons(i) = temp_parton
          exit FIND_Q
       end if
    end do FIND_Q

    !!! sort so that connected partons are next to each other, don't care about zeros
    do i = 1, size (final_partons)
       !!! ensure that final_partnons begins with a color (not an anticolor)
       if (final_partons(i)%c1 > 0 .and. final_partons(i)%c2 == 0) then
          if (i == 1) then
             exit
          else
             temp_parton = final_partons(1)
             final_partons(1) = final_partons(i)
             final_partons(i) = temp_parton
             exit
          end if
       end if
    end do

    do i = 1, size (final_partons) - 1
       !!! search for color partner and move it to i + 1
       PARTNERS: do j = i + 1, size (final_partons)
          if (final_partons(j)%c2 == final_partons(i)%c1) exit PARTNERS
       end do PARTNERS
       if (j > size (final_partons)) then
          print *, "no color connected parton found" !WRONG???
          print *, "particle: ", final_partons(i)%nr, " index: ", &
               final_partons(i)%c1
          exit
       end if
       temp_parton = final_partons(i + 1)
       final_partons(i + 1) = final_partons(j)
       final_partons(j) = temp_parton
    end do

    !!! transfering partons
    do i = 1, size (final_partons)
       prt = final_partons(i)
       N = N + 1
       K(N,1) = 2
       if (prt%c1 == 0) K(N,1) = 1       !!! end of color string
       K(N,2) = prt%type
       !K(N,3) = first
       K(N,3) = 0
       K(N,4) = 0
       K(N,5) = 0
       P(N,1:5) = prt%momentum_to_pythia6()
    end do
    deallocate (final_partons)
  end subroutine shower_transfer_final_partons_to_pythia


end submodule shower_core_s

