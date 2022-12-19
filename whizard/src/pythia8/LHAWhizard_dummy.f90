function new_whizard_lha () bind(C) result (whizard_lha)
  use iso_c_binding !NODEP!
  type(c_ptr) :: whizard_lha
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function new_whizard_lha

subroutine lhaup_whizard_delete (cptr) bind(C)
  use iso_c_binding !NODEP!
  ! Attribute value cannot have intent(inout).
  type(c_ptr), value :: cptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine lhaup_whizard_delete


function lhaup_whizard_set_init (cptr, beam_pdg, beam_energy, n_processes, unweighted, negative_weights) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int), dimension(2), intent(in) :: beam_pdg
  real(c_double), dimension(2), intent(in) :: beam_energy
  integer(c_int), intent(in), value :: n_processes
  logical(c_bool), intent(in), value :: unweighted, negative_weights
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function lhaup_whizard_set_init

function lhaup_whizard_set_process_parameters (cptr, process_id, cross_section, error, max_weight) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int), intent(in), value :: process_id
  real(c_double), intent(in), value :: cross_section, error, max_weight
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function lhaup_whizard_set_process_parameters

subroutine lhaup_whizard_list_init (cptr) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine lhaup_whizard_list_init

subroutine lhaup_whizard_list_event (cptr) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine lhaup_whizard_list_event

subroutine lhaup_whizard_set_event_process &
     (cptr, process_id, scale, alpha_qcd, alpha_qed, weight) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int), intent(in), value :: process_id
  real(c_double), intent(in), value :: scale, alpha_qcd, alpha_qed, weight
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine lhaup_whizard_set_event_process


function lhaup_whizard_set_event (cptr, process_id, n_particles, particle_set) bind(C) result (flag)
  use iso_c_binding !NODEP!
  use whizard_lha, only: lha_particle_t
  type(c_ptr), value :: cptr
  integer(c_int), intent(in), value :: process_id
  integer(c_int), intent(in), value :: n_particles
  ! IMPORTANT NOTE: Assumed-size array has to be defined by *.
  type(lha_particle_t), dimension(*), intent(in) :: particle_set
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** LHAupWhizard: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function lhaup_whizard_set_event
