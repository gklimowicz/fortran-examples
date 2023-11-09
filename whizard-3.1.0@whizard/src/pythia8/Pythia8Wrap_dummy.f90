! Tell the caller that this is NOT the true Pythia8 library
logical(c_bool) function pythia8_available () bind(C)
  use iso_c_binding !NODEP!
  pythia8_available = .false.
end function pythia8_available

function new_pythia8 (filename, subrun) bind(C) result (pythia)
  use iso_c_binding !NODEP!
  character(len=1, kind=c_char), intent(in) :: filename
  integer(c_int), intent(in) :: subrun
  type(c_ptr) :: pythia
  pythia = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function new_pythia8

subroutine pythia8_delete (pythia) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: pythia
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pythia8_delete

function pythia8_set_lhaup_ptr (cptr, whizard_lha) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  type(c_ptr), intent(in), value :: whizard_lha
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_set_lhaup_ptr

function pythia8_set_rndm_engine_ptr (pythia, whizard_rndm) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: pythia
  type(c_ptr), intent(in), value :: whizard_rndm
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_set_rndm_engine_ptr

function pythia8_read_file (cptr, filename, subrun) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  character(kind=c_char), dimension(*), intent(in) :: filename
  integer(c_int), intent(in), value :: subrun
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_read_file

function pythia8_read_string (cptr, str) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  character(kind=c_char), dimension(*), intent(in) :: str
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_read_string

function pythia8_init (cptr) bind(C) result (flag)
  use iso_c_binding !NODEP!
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_init

function pythia8_next (cptr) bind(C) result (flag)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  logical(c_bool) :: flag
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_next

subroutine pythia8_list_lha_event (cptr) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pythia8_list_lha_event

subroutine pythia8_list_event (cptr) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pythia8_list_event

function pythia8_get_event_size (cptr) bind(C) result(n)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int) :: n
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_event_size

function pythia8_get_single_event (cptr, index) bind(C) result (particle)
  use iso_c_binding !NODEP!
  use whizard_lha, only: lha_particle_t
  type(c_ptr), value :: cptr
  integer(c_int), value :: index
  type(lha_particle_t) :: particle
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_single_event


function pythia8_get_particle_status (cptr, index) bind(C) result (status)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  integer(c_int), intent(in), value :: index
  integer(c_int) :: status
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_particle_status

function pythia8_get_particle_id (cptr, index) bind(C) result (status)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  integer(c_int), intent(in), value :: index
  integer(c_int) :: status
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_particle_id

subroutine pythia8_get_particle_momentum (cptr, index, momentum) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  integer(c_int), intent(in), value :: index
  real(c_double), dimension(*), intent(out) :: momentum
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pythia8_get_particle_momentum

function pythia8_get_n_mothers (cptr) bind(C) result (n_mothers)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int) :: n_mothers
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_n_mothers

function pythia8_get_mother_array (cptr, i_prt, index) bind(C) result (mother)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int), value, intent(in) :: i_prt
  integer(c_int), value, intent(in) :: index
  integer(c_int) :: mother
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_mother_array

function pythia8_get_n_daughters (cptr) bind(C) result (n_daughters)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int) :: n_daughters
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop

end function pythia8_get_n_daughters

function pythia8_get_daughter_array (cptr, i_prt, index) bind(C) result (daughter)
  use iso_c_binding !NODEP!
  type(c_ptr), value :: cptr
  integer(c_int), value, intent(in) :: i_prt
  integer(c_int), value, intent(in) :: index
  integer(c_int) :: daughter
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_daughter_array

function pythia8_get_status_hepmc (cptr, i_prt) bind(C) result (status)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  integer(c_int), value, intent(in) :: i_prt
  integer(c_int) :: status
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pythia8_get_status_hepmc

subroutine pythia8_get_decay_vertex (cptr, i_prt, time, vertex) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  integer(c_int), value, intent(in) :: i_prt
  real(c_double), intent(out) :: time
  real(c_double), dimension(3), intent(out) :: vertex
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pythia8_get_decay_vertex

subroutine pythia8_get_production_vertex (cptr, i_prt, time, vertex) bind(C)
  use iso_c_binding !NODEP!
  type(c_ptr), intent(in), value :: cptr
  integer(c_int), value, intent(in) :: i_prt
  real(c_double), intent(out) :: time
  real(c_double), dimension(3), intent(out) :: vertex
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Pythia8: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pythia8_get_production_vertex
