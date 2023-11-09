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

submodule (hep_common) hep_common_s

  use io_units
  use diagnostics
  use numeric_utils
  use format_utils, only: refmt_tiny
  use physics_defs, only: pb_per_fb
  use physics_defs, only: HADRON_REMNANT
  use physics_defs, only: HADRON_REMNANT_SINGLET
  use physics_defs, only: HADRON_REMNANT_TRIPLET
  use physics_defs, only: HADRON_REMNANT_OCTET
  use xml
  use flavors
  use colors
  use subevents, only: PRT_BEAM, PRT_INCOMING, PRT_OUTGOING
  use subevents, only: PRT_UNDEFINED
  use subevents, only: PRT_VIRTUAL, PRT_RESONANT, PRT_BEAM_REMNANT

  implicit none

contains

  module subroutine heprup_init &
       (beam_pdg, beam_energy, n_processes, unweighted, negative_weights)
    integer, dimension(2), intent(in) :: beam_pdg
    real(default), dimension(2), intent(in) :: beam_energy
    integer, intent(in) :: n_processes
    logical, intent(in) :: unweighted
    logical, intent(in) :: negative_weights
    IDBMUP = beam_pdg
    EBMUP = beam_energy
    PDFGUP = -1
    PDFSUP = -1
    if (unweighted) then
       IDWTUP = 3
    else
       IDWTUP = 4
    end if
    if (negative_weights)  IDWTUP = - IDWTUP
    NPRUP = n_processes
  end subroutine heprup_init

  module subroutine assure_heprup (pset)
    type(particle_set_t), intent(in) :: pset
    integer :: i, num_id
    integer, parameter :: min_processes = 10
    num_id = 1
    if (LPRUP (num_id) /= 0)  return
    call heprup_init ( &
         [pset%prt(1)%get_pdg (), pset%prt(2)%get_pdg ()] , &
         [pset%prt(1)%p%p(0), pset%prt(2)%p%p(0)], &
         num_id, .false., .false.)
    do i = 1, (num_id / min_processes + 1) * min_processes
       call heprup_set_process_parameters (i = i, process_id = &
            i, cross_section = 1._default, error = 1._default)
    end do
  end subroutine assure_heprup

  module subroutine combine_lhef_with_particle_set &
       (particle_set, u, model_in, model_hadrons)
    type(particle_set_t), intent(inout) :: particle_set
    integer, intent(in) :: u
    class(model_data_t), intent(in), target :: model_in
    class(model_data_t), intent(in), target :: model_hadrons
    type(flavor_t) :: flv
    type(color_t) :: col
    class(model_data_t), pointer :: model
    type(particle_t), dimension(:), allocatable :: prt_tmp, prt
    integer :: i, j
    type(vector4_t) :: mom, d_mom
    integer, PARAMETER :: MAXLEN=200
    character(len=maxlen) :: string
    integer :: ibeg, n_tot, n_entries
    integer, dimension(:), allocatable :: relations, mothers, tbd
    INTEGER :: NUP,IDPRUP,IDUP,ISTUP
    real(kind=double) :: XWGTUP,SCALUP,AQEDUP,AQCDUP,VTIMUP,SPINUP
    integer :: MOTHUP(1:2), ICOLUP(1:2)
    real(kind=double) :: PUP(1:5)
    real(kind=default) :: pup_dum(1:5)
    character(len=5) :: buffer
    character(len=6) :: strfmt
    logical :: not_found
    logical :: debug_lhef = .false.
    STRFMT='(A000)'
    WRITE (STRFMT(3:5),'(I3)') MAXLEN

    if (debug_lhef)  call particle_set%write ()

    rewind (u)

    do
       read (u,*, END=501, ERR=502) STRING
       IBEG = 0
       do
          if (signal_is_pending ()) return
          IBEG = IBEG + 1
          ! Allow indentation.
          IF (STRING (IBEG:IBEG) .EQ. ' ' .and. IBEG < MAXLEN-6) cycle
          exit
       end do
       IF (string(IBEG:IBEG+6) /= '<event>' .and. &
            string(IBEG:IBEG+6) /= '<event ') cycle
       exit
    end do
    !!! Read first line of event info -> number of entries
    read (u, *, END=503, ERR=504) NUP, IDPRUP, XWGTUP, SCALUP, AQEDUP, AQCDUP
    n_tot = particle_set%get_n_tot ()
    allocate (prt_tmp (1:n_tot+NUP))
    allocate (relations (1:NUP), mothers (1:NUP), tbd(1:NUP))
    do i = 1, n_tot
       if (signal_is_pending ()) return
       prt_tmp (i) = particle_set%get_particle (i)
    end do
    !!! transfer particles from lhef to particle_set
    !!!...Read NUP subsequent lines with information on each particle.
    n_entries = 1
    mothers = 0
    relations = 0
    PARTICLE_LOOP: do I = 1, NUP
       read (u,*, END=200, ERR=505) IDUP, ISTUP, MOTHUP(1), MOTHUP(2), &
            ICOLUP(1), ICOLUP(2), (PUP (J),J=1,5), VTIMUP, SPINUP
       if (model_in%test_field (IDUP)) then
          model => model_in
       else if (model_hadrons%test_field (IDUP)) then
          model => model_hadrons
       else
          write (buffer, "(I5)") IDUP
          call msg_error ("Parton " // buffer // &
               " found neither in given model file nor in SM_hadrons")
          return
       end if
       if (debug_lhef) then
          print *, "IDUP, ISTUP, MOTHUP, PUP = ", IDUP, ISTUP, MOTHUP(1), &
             MOTHUP(2), PUP
       end if
       call flv%init (IDUP, model)
       if (IABS(IDUP) == 2212 .or. IABS(IDUP) == 2112) then
          ! PYTHIA sometimes sets color indices for protons and neutrons (?)
          ICOLUP (1) = 0
          ICOLUP (2) = 0
       end if
       call col%init_col_acl (ICOLUP (1), ICOLUP (2))
       !!! Settings for unpolarized particles
       ! particle_set%prt (oldsize+i)%hel = ??
       ! particle_set%prt (oldsize+i)%pol = ??
       if (MOTHUP(1) /= 0) then
          mothers(i) = MOTHUP(1)
       end if
       pup_dum = PUP
       if (pup_dum(4) < 1E-10_default)  cycle
       mom = vector4_moving (pup_dum (4), &
            vector3_moving ([pup_dum (1), pup_dum (2), pup_dum (3)]))
       not_found = .true.
       SCAN_PARTICLES: do j = 1, n_tot
          d_mom = prt_tmp(j)%get_momentum ()
          if (all (nearly_equal &
               (mom%p, d_mom%p, abs_smallness = 1.E-4_default)) .and. &
                (prt_tmp(j)%get_pdg () == IDUP)) then
             if (.not. prt_tmp(j)%get_status () == PRT_BEAM .or. &
                  .not. prt_tmp(j)%get_status () == PRT_BEAM_REMNANT) &
                  relations(i) = j
                  not_found = .false.
          end if
       end do SCAN_PARTICLES
       if (not_found) then
          if (debug_lhef) &
             print *, "Not found: adding particle"
          call prt_tmp(n_tot+n_entries)%set_flavor (flv)
          call prt_tmp(n_tot+n_entries)%set_color (col)
          call prt_tmp(n_tot+n_entries)%set_momentum (mom)
          if (MOTHUP(1) /= 0) then
             if (relations(MOTHUP(1)) /= 0) then
                call prt_tmp(n_tot+n_entries)%set_parents &
                     ([relations(MOTHUP(1))])
                call prt_tmp(relations(MOTHUP(1)))%add_child (n_tot+n_entries)
                if (prt_tmp(relations(MOTHUP(1)))%get_status () &
                     == PRT_OUTGOING) &
                     call prt_tmp(relations(MOTHUP(1)))%reset_status &
                     (PRT_VIRTUAL)
             end if
          end if
          call prt_tmp(n_tot+n_entries)%set_status (PRT_OUTGOING)
          if (debug_lhef) call prt_tmp(n_tot+n_entries)%write ()
          n_entries = n_entries + 1
       end if
    end do PARTICLE_LOOP
    do i = 1, n_tot
       if (prt_tmp(i)%get_status () == PRT_OUTGOING .and. &
           prt_tmp(i)%get_n_children () /= 0) then
                call prt_tmp(i)%reset_status (PRT_VIRTUAL)
       end if
    end do

    allocate (prt (1:n_tot+n_entries-1))
    prt = prt_tmp (1:n_tot+n_entries-1)
    ! transfer to particle_set
    call particle_set%replace (prt)
    deallocate (prt, prt_tmp)

    if (debug_lhef) then
       call particle_set%write ()
       print *, "combine_lhef_with_particle_set"
       ! stop
    end if

200 continue
    return

501 write(*,*) "READING LHEF failed 501"
    return
502 write(*,*) "READING LHEF failed 502"
    return
503 write(*,*) "READING LHEF failed 503"
    return
504 write(*,*) "READING LHEF failed 504"
    return
505 write(*,*) "READING LHEF failed 505"
    return
  end subroutine combine_lhef_with_particle_set

  module subroutine w2p_write_lhef_event (unit)
    integer, intent(in) :: unit
    type(xml_tag_t), allocatable :: tag_lhef, tag_head, tag_init, &
         tag_event, tag_gen_n, tag_gen_v
    if (debug_on) call msg_debug (D_EVENTS, "w2p_write_lhef_event")
    allocate (tag_lhef, tag_head, tag_init, tag_event, &
         tag_gen_n, tag_gen_v)
    call tag_lhef%init (var_str ("LesHouchesEvents"), &
         [xml_attribute (var_str ("version"), var_str ("1.0"))], .true.)
    call tag_head%init (var_str ("header"), .true.)
    call tag_init%init (var_str ("init"), .true.)
    call tag_event%init (var_str ("event"), .true.)
    call tag_gen_n%init (var_str ("generator_name"), .true.)
    call tag_gen_v%init (var_str ("generator_version"), .true.)
    call tag_lhef%write (unit); write (unit, *)
    call tag_head%write (unit); write (unit, *)
    write (unit, "(2x)", advance = "no")
    call tag_gen_n%write (var_str ("WHIZARD"), unit)
    write (unit, *)
    write (unit, "(2x)", advance = "no")
    call tag_gen_v%write (var_str ("3.1.0"), unit)
    write (unit, *)
    call tag_head%close (unit); write (unit, *)
    call tag_init%write (unit); write (unit, *)
    call heprup_write_lhef (unit)
    call tag_init%close (unit); write (unit, *)
    call tag_event%write (unit); write (unit, *)
    call hepeup_write_lhef (unit)
    call tag_event%close (unit); write (unit, *)
    call tag_lhef%close (unit); write (unit, *)
    deallocate (tag_lhef, tag_head, tag_init, tag_event, &
         tag_gen_n, tag_gen_v)
  end subroutine w2p_write_lhef_event

  module subroutine heprup_get_run_parameters &
       (beam_pdg, beam_energy, n_processes, unweighted, negative_weights)
    integer, dimension(2), intent(out), optional :: beam_pdg
    real(default), dimension(2), intent(out), optional :: beam_energy
    integer, intent(out), optional :: n_processes
    logical, intent(out), optional :: unweighted
    logical, intent(out), optional :: negative_weights
    if (present (beam_pdg))  beam_pdg = IDBMUP
    if (present (beam_energy))  beam_energy = EBMUP
    if (present (n_processes))  n_processes = NPRUP
    if (present (unweighted)) then
       select case (abs (IDWTUP))
       case (3)
          unweighted = .true.
       case (4)
          unweighted = .false.
       case (1,2)  !!! not supported by WHIZARD
          unweighted = .false.
       case default
          call msg_fatal ("HEPRUP: unsupported IDWTUP value")
       end select
    end if
    if (present (negative_weights)) then
       negative_weights = IDWTUP < 0
    end if
  end subroutine heprup_get_run_parameters

  module subroutine heprup_set_lhapdf_id (i_beam, pdf_id)
    integer, intent(in) :: i_beam, pdf_id
    PDFGUP(i_beam) = 0
    PDFSUP(i_beam) = pdf_id
  end subroutine heprup_set_lhapdf_id

  module subroutine heprup_set_process_parameters &
       (i, process_id, cross_section, error, max_weight, is_width)
    integer, intent(in) :: i, process_id
    real(default), intent(in), optional :: cross_section, error, max_weight
    logical, intent(in), optional :: is_width
    logical :: is_w
    is_w = .false.
    if (present (is_width))  is_w = is_width
    LPRUP(i) = process_id
    if (present (cross_section)) then
       if (is_w) then
          XSECUP(i) = cross_section
       else
          XSECUP(i) = cross_section * pb_per_fb
       end if
    else
       XSECUP(i) = 0
    end if
    if (present (error)) then
       if (is_w) then
          XERRUP(i) = error
       else
          XERRUP(i) = error * pb_per_fb
       end if
    else
       XERRUP(i) = 0
    end if
    select case (IDWTUP)
    case (3);  XMAXUP(i) = 1
    case (4)
       if (present (max_weight)) then
          if (is_w) then
             XMAXUP(i) = max_weight
          else
             XMAXUP(i) = max_weight * pb_per_fb
          end if
       else
          XMAXUP(i) = 0
       end if
    end select
  end subroutine heprup_set_process_parameters

  module subroutine heprup_get_process_parameters  &
       (i, process_id, cross_section, error, max_weight, is_width)
    integer, intent(in) :: i
    integer, intent(out), optional :: process_id
    real(default), intent(out), optional :: cross_section, error, max_weight
    logical, intent(in), optional :: is_width
    logical :: is_w
    is_w = .false.
    if (present (is_width))  is_w = is_width
    if (present (process_id))  process_id = LPRUP(i)
    if (present (cross_section)) then
       if (is_w) then
          cross_section = XSECUP(i)
       else
          cross_section = XSECUP(i) / pb_per_fb
       end if
    end if
    if (present (error)) then
       if (is_w) then
          error = XERRUP(i)
       else
          error = XERRUP(i) / pb_per_fb
       end if
    end if
    if (present (max_weight)) then
       select case (IDWTUP)
       case (3)
          max_weight = 1
       case (4)
          max_weight = XMAXUP(i) / pb_per_fb
       case (1,2)   !!! not supported by WHIZARD
          max_weight = 0
       case default
          call msg_fatal ("HEPRUP: unsupported IDWTUP value")
       end select
    end if
  end subroutine heprup_get_process_parameters

  module subroutine heprup_write_verbose (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)")  "HEPRUP Common Block"
    write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "IDBMUP", IDBMUP, &
         "PDG code of beams"
    write (u, "(3x,A6,' = ',G12.5,1x,G12.5,8x,A)")  "EBMUP ", EBMUP, &
         "Energy of beams in GeV"
    write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "PDFGUP", PDFGUP, &
         "PDF author group [-1 = undefined]"
    write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "PDFSUP", PDFSUP, &
         "PDF set ID       [-1 = undefined]"
    write (u, "(3x,A6,' = ',I9,3x,1x,9x,3x,8x,A)")  "IDWTUP", IDWTUP, &
         "LHA code for event weight mode"
    write (u, "(3x,A6,' = ',I9,3x,1x,9x,3x,8x,A)")  "NPRUP ", NPRUP, &
         "Number of user subprocesses"
    do i = 1, NPRUP
       write (u, "(1x,A,I0)")  "Subprocess #", i
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "XSECUP", XSECUP(i), &
            "Cross section in pb"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "XERRUP", XERRUP(i), &
            "Cross section error in pb"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "XMAXUP", XMAXUP(i), &
            "Maximum event weight (cf. IDWTUP)"
       write (u, "(3x,A6,' = ',I9,3x,1x,12x,8x,A)")  "LPRUP ", LPRUP(i), &
            "Subprocess ID"
    end do
  end subroutine heprup_write_verbose

  module subroutine heprup_write_lhef (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(2(1x,I0),2(1x,ES17.10),6(1x,I0))") &
         IDBMUP, EBMUP, PDFGUP, PDFSUP, IDWTUP, NPRUP
    do i = 1, NPRUP
       write (u, "(3(1x,ES17.10),1x,I0)") &
            XSECUP(i), XERRUP(i), XMAXUP(i), LPRUP(i)
    end do
  end subroutine heprup_write_lhef

  module subroutine heprup_write_ascii (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(2(1x,I0),2(1x,ES17.10),6(1x,I0))") &
         IDBMUP, EBMUP, PDFGUP, PDFSUP, IDWTUP, NPRUP
    do i = 1, NPRUP
       write (u, "(3(1x,ES17.10),1x,I0)") &
            XSECUP(i), XERRUP(i), XMAXUP(i), LPRUP(i)
    end do
  end subroutine heprup_write_ascii

  module subroutine heprup_read_lhef (u)
    integer, intent(in) :: u
    integer :: i
    read (u, *) &
         IDBMUP, EBMUP, PDFGUP, PDFSUP, IDWTUP, NPRUP
    do i = 1, NPRUP
       read (u, *) &
            XSECUP(i), XERRUP(i), XMAXUP(i), LPRUP(i)
    end do
  end subroutine heprup_read_lhef

  module subroutine hepeup_init (n_tot)
    integer, intent(in) :: n_tot
    NUP = n_tot
    IDPRUP = 0
    XWGTUP = 1
    SCALUP = -1
    AQEDUP = -1
    AQCDUP = -1
  end subroutine hepeup_init

  module subroutine hepeup_set_event_parameters &
       (proc_id, weight, scale, alpha_qed, alpha_qcd)
    integer, intent(in), optional :: proc_id
    real(default), intent(in), optional :: &
         weight, scale, alpha_qed, alpha_qcd
    if (present (proc_id))   IDPRUP = proc_id
    if (present (weight))    XWGTUP = weight
    if (present (scale))     SCALUP = scale
    if (present (alpha_qed)) AQEDUP = alpha_qed
    if (present (alpha_qcd)) AQCDUP = alpha_qcd
  end subroutine hepeup_set_event_parameters

  module subroutine hepeup_get_event_parameters &
       (proc_id, weight, scale, alpha_qed, alpha_qcd)
    integer, intent(out), optional :: proc_id
    real(default), intent(out), optional :: &
         weight, scale, alpha_qed, alpha_qcd
    if (present (proc_id))   proc_id   = IDPRUP
    if (present (weight))    weight    = XWGTUP
    if (present (scale))     scale     = SCALUP
    if (present (alpha_qed)) alpha_qed = AQEDUP
    if (present (alpha_qcd)) alpha_qcd = AQCDUP
  end subroutine hepeup_get_event_parameters

  module subroutine hepeup_set_particle (i, pdg, status, parent, col, p, m2)
    integer, intent(in) :: i
    integer, intent(in) :: pdg, status
    integer, dimension(:), intent(in) :: parent
    type(vector4_t), intent(in) :: p
    integer, dimension(2), intent(in) :: col
    real(default), intent(in) :: m2
    if (i > MAXNUP) then
       call msg_error (arr=[ &
            var_str ("Too many particles in HEPEUP common block. " // &
                            "If this happened "), &
            var_str ("during event output, your events will be " // &
                            "invalid; please consider "), &
            var_str ("switching to a modern event format like HEPMC. " // &
                            "If you are not "), &
            var_str ("using an old, HEPEUP based format and " // &
                            "nevertheless get this error,"), &
            var_str ("please notify the WHIZARD developers,") ])
       return
    end if
    IDUP(i) = pdg
    select case (status)
    case (PRT_BEAM);         ISTUP(i) = -9
    case (PRT_INCOMING);     ISTUP(i) = -1
    case (PRT_BEAM_REMNANT); ISTUP(i) =  3
    case (PRT_OUTGOING);     ISTUP(i) =  1
    case (PRT_RESONANT);     ISTUP(i) =  2
    case (PRT_VIRTUAL);      ISTUP(i) =  3
    case default;            ISTUP(i) =  0
    end select
    select case (size (parent))
    case (0);      MOTHUP(:,i) = 0
    case (1);      MOTHUP(1,i) = parent(1); MOTHUP(2,i) = 0
    case default;  MOTHUP(:,i) = [ parent(1), parent(size (parent)) ]
    end select
    if (col(1) > 0) then
       ICOLUP(1,i) = 500 + col(1)
    else
       ICOLUP(1,i) = 0
    end if
    if (col(2) > 0) then
       ICOLUP(2,i) = 500 + col(2)
    else
       ICOLUP(2,i) = 0
    end if
    PUP(1:3,i) = refmt_tiny (vector3_get_components (space_part (p)))
    PUP(4,i) = refmt_tiny (energy (p))
    PUP(5,i) = refmt_tiny (sign (sqrt (abs (m2)), m2))
    VTIMUP(i) = 0
    SPINUP(i) = 9
  end subroutine hepeup_set_particle

  module subroutine hepeup_set_particle_lifetime (i, lifetime)
    integer, intent(in) :: i
    real(default), intent(in) :: lifetime
    VTIMUP(i) = lifetime
  end subroutine hepeup_set_particle_lifetime

  module subroutine hepeup_set_particle_spin_pol (i, p, pol, p_mother)
    integer, intent(in) :: i
    type(vector4_t), intent(in) :: p
    type(polarization_t), intent(in) :: pol
    type(vector4_t), intent(in) :: p_mother
    type(vector3_t) :: s3, p3
    type(vector4_t) :: s4
    s3 = vector3_moving (pol%get_axis ())
    p3 = space_part (p)
    s4 = rotation_to_2nd (3, p3) * vector4_moving (0._default, s3)
    SPINUP(i) = enclosed_angle_ct (s4, p_mother)
  end subroutine hepeup_set_particle_spin_pol

  module subroutine hepeup_get_particle (i, pdg, status, parent, col, p, m2)
    integer, intent(in) :: i
    integer, intent(out), optional :: pdg, status
    integer, dimension(:), intent(out), optional :: parent
    type(vector4_t), intent(out), optional :: p
    integer, dimension(2), intent(out), optional :: col
    real(default), dimension(5,MAXNUP) :: pup_def
    real(default), intent(out), optional :: m2
    if (present (pdg))  pdg = IDUP(i)
    if (present (status)) then
       select case (ISTUP(i))
       case (-9);  status = PRT_BEAM
       case (-1);  status = PRT_INCOMING
       case (1);   status = PRT_OUTGOING
       case (2);   status = PRT_RESONANT
       case (3);
          select case (abs (IDUP(i)))
          case (HADRON_REMNANT, HADRON_REMNANT_SINGLET, &
               HADRON_REMNANT_TRIPLET, HADRON_REMNANT_OCTET)
             status = PRT_BEAM_REMNANT
          case default
             status = PRT_VIRTUAL
          end select
       case default
          status = PRT_UNDEFINED
       end select
    end if
    if (present (parent)) then
       select case (size (parent))
       case (0)
       case (1);    parent(1) = MOTHUP(1,i)
       case (2);    parent = MOTHUP(:,i)
       end select
    end if
    if (present (col)) then
       col = ICOLUP(:,i)
    end if
    if (present (p)) then
       pup_def = PUP
       p = vector4_moving (pup_def(4,i), vector3_moving (pup_def(1:3,i)))
    end if
    if (present (m2)) then
       m2 = sign (PUP(5,i) ** 2, PUP(5,i))
    end if
  end subroutine hepeup_get_particle

  module subroutine hepevt_init (n_tot, n_out)
    integer, intent(in) :: n_tot, n_out
    NHEP              = n_tot
    NEVHEP            = 0
    idruplh           = 0
    hepevt_n_out      = n_out
    hepevt_n_remnants = 0
    hepevt_weight     = 1
    eventweightlh     = 1
    hepevt_function_value = 0
    hepevt_function_ratio = 1
    alphaqcdlh        = -1
    alphaqedlh        = -1
    scalelh           = -1
  end subroutine hepevt_init

  module subroutine hepevt_set_event_parameters &
       (proc_id, weight, function_value, function_ratio, &
       alpha_qcd, alpha_qed, scale, i_evt)
    integer, intent(in), optional :: proc_id
    integer, intent(in), optional :: i_evt
    real(default), intent(in), optional :: weight, function_value, &
       function_ratio, alpha_qcd, alpha_qed, scale
    if (present (proc_id))  idruplh = proc_id
    if (present (i_evt))  NEVHEP = i_evt
    if (present (weight)) then
       hepevt_weight = weight
       eventweightlh = weight
    end if
    if (present (function_value)) hepevt_function_value = &
         function_value
    if (present (function_ratio)) hepevt_function_ratio = &
         function_ratio
    if (present (alpha_qcd))  alphaqcdlh = alpha_qcd
    if (present (alpha_qed))  alphaqedlh = alpha_qed
    if (present (scale))  scalelh(1) = scale
    if (present (i_evt))  NEVHEP = i_evt
  end subroutine hepevt_set_event_parameters

  module subroutine hepevt_set_particle &
       (i, pdg, status, parent, child, p, m2, hel, vtx, &
       col, pol_status, pol, fill_hepev4)
    integer, intent(in) :: i
    integer, intent(in) :: pdg, status
    integer, dimension(:), intent(in) :: parent
    integer, dimension(:), intent(in) :: child
    logical, intent(in), optional :: fill_hepev4
    type(vector4_t), intent(in) :: p
    real(default), intent(in) :: m2
    integer, dimension(2), intent(in) :: col
    integer, intent(in) :: pol_status
    integer, intent(in) :: hel
    type(polarization_t), intent(in), optional :: pol
    type(vector4_t), intent(in) :: vtx
    logical :: hepev4
    hepev4 = .false.; if (present (fill_hepev4))  hepev4 = fill_hepev4
    IDHEP(i) = pdg
    select case (status)
      case (PRT_BEAM);      ISTHEP(i) = 2
      case (PRT_INCOMING);  ISTHEP(i) = 2
      case (PRT_OUTGOING);  ISTHEP(i) = 1
      case (PRT_VIRTUAL);   ISTHEP(i) = 2
      case (PRT_RESONANT);  ISTHEP(i) = 2
      case default;         ISTHEP(i) = 0
    end select
    select case (size (parent))
    case (0);      JMOHEP(:,i) = 0
    case (1);      JMOHEP(1,i) = parent(1); JMOHEP(2,i) = 0
    case default;  JMOHEP(:,i) = [ parent(1), parent(size (parent)) ]
    end select
    select case (size (child))
    case (0);      JDAHEP(:,i) = 0
    case (1);      JDAHEP(:,i) = child(1)
    case default;  JDAHEP(:,i) = [ child(1), child(size (child)) ]
    end select
    PHEP(1:3,i) = refmt_tiny (vector3_get_components (space_part (p)))
    PHEP(4,i) = refmt_tiny (energy (p))
    PHEP(5,i) = refmt_tiny (sign (sqrt (abs (m2)), m2))
    VHEP(1:3,i) = vtx%p(1:3)
    VHEP(4,i) = vtx%p(0)
    hepevt_pol(i) = hel
    if (hepev4) then
       if (col(1) > 0) then
          icolorflowlh(1,i) = 500 + col(1)
       else
          icolorflowlh(1,i) = 0
       end if
       if (col(2) > 0) then
          icolorflowlh(2,i) = 500 + col(2)
       else
          icolorflowlh(2,i) = 0
       end if
       if (present (pol) .and. &
            pol_status == PRT_GENERIC_POLARIZATION) then
          if (pol%is_polarized ()) &
             spinlh(:,i) = pol%get_axis ()
       else
          spinlh(:,i) = zero
          spinlh(3,i) = hel
       end if
    end if
  end subroutine hepevt_set_particle

  module subroutine hepevt_write_verbose (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)")  "HEPEVT Common Block"
    write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)")  "NEVHEP", NEVHEP, &
         "Event number"
    write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)")  "NHEP  ", NHEP, &
         "Number of particles in event"
    do i = 1, NHEP
       write (u, "(1x,A,I0)")  "Particle #", i
       write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)", advance="no") &
            "ISTHEP", ISTHEP(i), "Status code: "
       select case (ISTHEP(i))
       case ( 0);  write (u, "(A)")  "null entry"
       case ( 1);  write (u, "(A)")  "outgoing"
       case ( 2);  write (u, "(A)")  "decayed"
       case ( 3);  write (u, "(A)")  "documentation"
       case (4:10);  write (u, "(A)")  "[unspecified]"
       case (11:200);  write (u, "(A)")  "[model-specific]"
       case (201:);  write (u, "(A)")  "[user-defined]"
       case default;  write (u, "(A)")  "[undefined]"
       end select
       write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)")  "IDHEP ", IDHEP(i), &
            "PDG code of particle"
       write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "JMOHEP", JMOHEP(:,i), &
            "Index of first/second mother"
       write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "JDAHEP", JDAHEP(:,i), &
            "Index of first/last daughter"
       write (u, "(3x,A6,' = ',ES12.5,1x,ES12.5,8x,A)")  "PHEP12", &
            PHEP(1:2,i), "Transversal momentum (x/y) in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "PHEP3 ", PHEP(3,i), &
            "Longitudinal momentum (z) in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "PHEP4 ", PHEP(4,i), &
            "Energy in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "PHEP5 ", PHEP(5,i), &
            "Invariant mass in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,ES12.5,8x,A)")  "VHEP12", VHEP(1:2,i), &
            "Transversal displacement (xy) in mm"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "VHEP3 ", VHEP(3,i), &
            "Longitudinal displacement (z) in mm"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "VHEP4 ", VHEP(4,i), &
            "Production time in mm"
    end do
  end subroutine hepevt_write_verbose

  module subroutine hepeup_write_verbose (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(A)")  "HEPEUP Common Block"
    write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)")  "NUP   ", NUP, &
         "Number of particles in event"
    write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)")  "IDPRUP", IDPRUP, &
         "Subprocess ID"
    write (u, "(3x,A6,' = ',ES12.5,1x,20x,A)")  "XWGTUP", XWGTUP, &
         "Event weight"
    write (u, "(3x,A6,' = ',ES12.5,1x,20x,A)")  "SCALUP", SCALUP, &
         "Event energy scale in GeV"
    write (u, "(3x,A6,' = ',ES12.5,1x,20x,A)")  "AQEDUP", AQEDUP, &
         "QED coupling [-1 = undefined]"
    write (u, "(3x,A6,' = ',ES12.5,1x,20x,A)")  "AQCDUP", AQCDUP, &
         "QCD coupling [-1 = undefined]"
    do i = 1, NUP
       write (u, "(1x,A,I0)")  "Particle #", i
       write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)")  "IDUP  ", IDUP(i), &
            "PDG code of particle"
       write (u, "(3x,A6,' = ',I9,3x,1x,20x,A)", advance="no") &
            "ISTUP ", ISTUP(i), "Status code: "
       select case (ISTUP(i))
       case (-1);  write (u, "(A)")  "incoming"
       case ( 1);  write (u, "(A)")  "outgoing"
       case (-2);  write (u, "(A)")  "spacelike"
       case ( 2);  write (u, "(A)")  "resonance"
       case ( 3);  write (u, "(A)")  "resonance (doc)"
       case (-9);  write (u, "(A)")  "beam"
       case default;  write (u, "(A)")  "[undefined]"
       end select
       write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "MOTHUP", MOTHUP(:,i), &
            "Index of first/last mother"
       write (u, "(3x,A6,' = ',I9,3x,1x,I9,3x,8x,A)")  "ICOLUP", ICOLUP(:,i), &
            "Color/anticolor flow index"
       write (u, "(3x,A6,' = ',ES12.5,1x,ES12.5,8x,A)")  "PUP1/2", PUP(1:2,i), &
            "Transversal momentum (x/y) in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "PUP3  ", PUP(3,i), &
            "Longitudinal momentum (z) in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "PUP4  ", PUP(4,i), &
            "Energy in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "PUP5  ", PUP(5,i), &
            "Invariant mass in GeV"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "VTIMUP", VTIMUP(i), &
            "Invariant lifetime in mm"
       write (u, "(3x,A6,' = ',ES12.5,1x,12x,8x,A)")  "SPINUP", SPINUP(i), &
            "cos(spin angle) [9 = undefined]"
    end do
  end subroutine hepeup_write_verbose

  module subroutine hepeup_write_lhef (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    if (debug_on) call msg_debug (D_EVENTS, "hepeup_write_lhef")
    if (debug_on) call msg_debug2 (D_EVENTS, "ID IST MOTH ICOL P VTIM SPIN")
    write (u, "(2(1x,I0),4(1x,ES17.10))") &
         NUP, IDPRUP, XWGTUP, SCALUP, AQEDUP, AQCDUP
    do i = 1, NUP
       write (u, "(6(1x,I0),7(1x,ES17.10))") &
            IDUP(i), ISTUP(i), MOTHUP(:,i), ICOLUP(:,i), &
            PUP(:,i), VTIMUP(i), SPINUP(i)
       if (debug2_active (D_EVENTS)) then
          write (msg_buffer, "(6(1x,I0),7(1x,ES17.10))") &
               IDUP(i), ISTUP(i), MOTHUP(:,i), ICOLUP(:,i), &
               PUP(:,i), VTIMUP(i), SPINUP(i)
          call msg_message ()
       end if
    end do
  end subroutine hepeup_write_lhef

  module subroutine hepeup_write_lha (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    integer, dimension(MAXNUP) :: spin_up
    spin_up = int(SPINUP)
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(2(1x,I5),1x,ES17.10,3(1x,ES13.6))") &
         NUP, IDPRUP, XWGTUP, SCALUP, AQEDUP, AQCDUP
    write (u, "(500(1x,I5))") IDUP(:NUP)
    write (u, "(500(1x,I5))") MOTHUP(1,:NUP)
    write (u, "(500(1x,I5))") MOTHUP(2,:NUP)
    write (u, "(500(1x,I5))") ICOLUP(1,:NUP)
    write (u, "(500(1x,I5))") ICOLUP(2,:NUP)
    write (u, "(500(1x,I5))") ISTUP(:NUP)
    write (u, "(500(1x,I5))") spin_up(:NUP)
    do i = 1, NUP
            write (u, "(1x,I5,4(1x,ES17.10))") i, PUP([ 4,1,2,3 ], i)
    end do

  end subroutine hepeup_write_lha

  module subroutine hepevt_write_hepevt (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3(1x,I0),(1x,ES17.10))") &
         NHEP, hepevt_n_out, hepevt_n_remnants, hepevt_weight
    do i = 1, NHEP
       write (u, "(7(1x,I0))") &
            ISTHEP(i), IDHEP(i), JMOHEP(:,i), JDAHEP(:,i), hepevt_pol(i)
       write (u, "(5(1x,ES17.10))") PHEP(:,i)
       write (u, "(5(1x,ES17.10))") VHEP(:,i), 0.d0
    end do
  end subroutine hepevt_write_hepevt

  module subroutine hepevt_write_ascii (unit, long)
    integer, intent(in), optional :: unit
    logical, intent(in) :: long
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3(1x,I0),(1x,ES17.10))") &
         NHEP, hepevt_n_out, hepevt_n_remnants, hepevt_weight
    do i = 1, NHEP
       if (ISTHEP(i) /= 1)  cycle
       write (u, "(2(1x,I0))") IDHEP(i), hepevt_pol(i)
       write (u, "(5(1x,ES17.10))") PHEP(:,i)
    end do
    if (long) then
       write (u, "(2(1x,ES17.10))") &
            hepevt_function_value, hepevt_function_ratio
    end if
  end subroutine hepevt_write_ascii

  module subroutine hepevt_write_athena (unit)
    integer, intent(in), optional :: unit
    integer :: u, i, num_event
    num_event = 0
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(2(1x,I0))") NEVHEP, NHEP
    do i = 1, NHEP
       write (u, "(7(1x,I0))") &
            i, ISTHEP(i), IDHEP(i), JMOHEP(:,i), JDAHEP(:,i)
       write (u, "(5(1x,ES17.10))") PHEP(:,i)
       write (u, "(5(1x,ES17.10))") VHEP(1:4,i)
    end do
  end subroutine hepevt_write_athena

  module subroutine hepevt_write_mokka (unit)
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3(1x,I0),(1x,ES17.10))") &
         NHEP, hepevt_n_out, hepevt_n_remnants, hepevt_weight
    do i = 1, NHEP
       write (u, "(4(1x,I0),4(1x,ES17.10))") &
            ISTHEP(i), IDHEP(i), JDAHEP(1,i), JDAHEP(2,i), &
            PHEP(1:3,i), PHEP(5,i)
    end do
  end subroutine hepevt_write_mokka

  module subroutine hepeup_read_lhef (u)
    integer, intent(in) :: u
    integer :: i
    read (u, *) &
         NUP, IDPRUP, XWGTUP, SCALUP, AQEDUP, AQCDUP
    do i = 1, NUP
       read (u, *) &
            IDUP(i), ISTUP(i), MOTHUP(:,i), ICOLUP(:,i), &
            PUP(:,i), VTIMUP(i), SPINUP(i)
    end do
  end subroutine hepeup_read_lhef

  module subroutine hepeup_from_particle_set (pset_in, &
     keep_beams, keep_remnants, tauola_convention)
    type(particle_set_t), intent(in) :: pset_in
    type(particle_set_t), target :: pset
    logical, intent(in), optional :: keep_beams
    logical, intent(in), optional :: keep_remnants
    logical, intent(in), optional :: tauola_convention
    integer :: i, n_parents, status, n_tot
    integer, dimension(1) :: i_mother
    logical :: kr, tc
    kr = .true.;  if (present (keep_remnants))  kr = keep_remnants
    tc = .false.;  if (present (tauola_convention))  tc = tauola_convention
    call pset_in%filter_particles (pset, real_parents = .true. , &
       keep_beams = keep_beams, keep_virtuals = .false.)
    n_tot = pset%get_n_tot ()
    call hepeup_init (n_tot)
    do i = 1, n_tot
       associate (prt => pset%prt(i))
         status = prt%get_status ()
         if (kr .and. status == PRT_BEAM_REMNANT &
                .and. prt%get_n_children () == 0) &
              status = PRT_OUTGOING
         call hepeup_set_particle (i, &
              prt%get_pdg (), &
              status, &
              prt%get_parents (), &
              prt%get_color (), &
              prt%get_momentum (), &
              prt%get_p2 ())
         n_parents = prt%get_n_parents ()
         call hepeup_set_particle_lifetime (i, &
              prt%get_lifetime ())
         if (.not. tc) then
            if (n_parents == 1) then
               i_mother = prt%get_parents ()
               select case (prt%get_polarization_status ())
               case (PRT_GENERIC_POLARIZATION)
                  call hepeup_set_particle_spin (i, &
                       prt%get_momentum (), &
                       prt%get_polarization (), &
                       pset%prt(i_mother(1))%get_momentum ())
               end select
            end if
         else
            select case (prt%get_polarization_status ())
            case (PRT_DEFINITE_HELICITY)
              SPINUP(i) = prt%get_helicity()
            end select
         end if
       end associate
    end do
  end subroutine hepeup_from_particle_set

  module subroutine hepeup_to_particle_set &
       (particle_set, recover_beams, model, alt_model)
    type(particle_set_t), intent(inout), target :: particle_set
    logical, intent(in), optional :: recover_beams
    class(model_data_t), intent(in), target :: model, alt_model
    type(particle_t), dimension(:), allocatable :: prt
    integer, dimension(2) :: parent
    integer, dimension(:), allocatable :: child
    integer :: i, j, k, pdg, status
    type(flavor_t) :: flv
    type(color_t) :: col
    integer, dimension(2) :: c
    type(vector4_t) :: p
    real(default) :: p2
    logical :: reconstruct
    integer :: off
    if (present (recover_beams)) then
       reconstruct = recover_beams .and. .not. all (ISTUP(1:2) == PRT_BEAM)
    else
       reconstruct = .false.
    end if
    if (reconstruct) then
       off = 4
    else
       off = 0
    end if
    allocate (prt (NUP + off), child (NUP + off))
    do i = 1, NUP
       k = i + off
       call hepeup_get_particle (i, pdg, status, col = c, p = p, m2 = p2)
       call flv%init (pdg, model, alt_model)
       call prt(k)%set_flavor (flv)
       call prt(k)%reset_status (status)
       call col%init (c)
       call prt(k)%set_color (col)
       call prt(k)%set_momentum (p, p2)
       where (MOTHUP(:,i) /= 0)
          parent = MOTHUP(:,i) + off
       elsewhere
          parent = 0
       end where
       call prt(k)%set_parents (parent)
       child = [(j, j = 1 + off, NUP + off)]
       where (MOTHUP(1,:NUP) /= i .and. MOTHUP(2,:NUP) /= i)  child = 0
       call prt(k)%set_children (child)
    end do
    if (reconstruct) then
       do k = 1, 2
          call prt(k)%reset_status (PRT_BEAM)
          call prt(k)%set_children ([k+2,k+4])
       end do
       do k = 3, 4
          call prt(k)%reset_status (PRT_BEAM_REMNANT)
          call prt(k)%set_parents ([k-2])
       end do
       do k = 5, 6
          call prt(k)%set_parents ([k-4])
       end do
    end if
    call particle_set%replace (prt)
  end subroutine hepeup_to_particle_set

  module subroutine hepevt_from_particle_set &
       (particle_set, keep_beams, keep_remnants, ensure_order, fill_hepev4)
    type(particle_set_t), intent(in) :: particle_set
    type(particle_set_t), target :: pset_hepevt, pset_tmp
    logical, intent(in), optional :: keep_beams
    logical, intent(in), optional :: keep_remnants
    logical, intent(in), optional :: ensure_order
    logical, intent(in), optional :: fill_hepev4
    integer :: i, status, n_tot
    logical :: activate_remnants, ensure
    activate_remnants = .true.
    if (present (keep_remnants))  activate_remnants = keep_remnants
    ensure = .false.
    if (present (ensure_order))  ensure = ensure_order
    call particle_set%filter_particles (pset_tmp, real_parents = .true., &
       keep_virtuals = .false., keep_beams = keep_beams)
    if (ensure) then
       call pset_tmp%to_hepevt_form (pset_hepevt)
    else
       pset_hepevt = pset_tmp
    end if
    n_tot = pset_hepevt%get_n_tot ()
    call hepevt_init (n_tot, pset_hepevt%get_n_out ())
    do i = 1, n_tot
       associate (prt => pset_hepevt%prt(i))
         status = prt%get_status ()
         if (activate_remnants &
              .and. status == PRT_BEAM_REMNANT &
              .and. prt%get_n_children () == 0) &
              status = PRT_OUTGOING
         select case (prt%get_polarization_status ())
         case (PRT_GENERIC_POLARIZATION)
            call hepevt_set_particle (i, &
                 prt%get_pdg (), status, &
                 prt%get_parents (), &
                 prt%get_children (), &
                 prt%get_momentum (), &
                 prt%get_p2 (), &
                 prt%get_helicity (), &
                 prt%get_vertex (), &
                 prt%get_color (), &
                 prt%get_polarization_status (), &
                 pol = prt%get_polarization (), &
                 fill_hepev4 = fill_hepev4)
         case default
            call hepevt_set_particle (i, &
                 prt%get_pdg (), status, &
                 prt%get_parents (), &
                 prt%get_children (), &
                 prt%get_momentum (), &
                 prt%get_p2 (), &
                 prt%get_helicity (), &
                 prt%get_vertex (), &
                 prt%get_color (), &
                 prt%get_polarization_status (), &
                 fill_hepev4 = fill_hepev4)
         end select
       end associate
    end do
    call pset_hepevt%final ()
  end subroutine hepevt_from_particle_set


end submodule hep_common_s

