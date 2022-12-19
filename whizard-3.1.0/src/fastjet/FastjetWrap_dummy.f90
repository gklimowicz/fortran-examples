! Tell the caller that this is NOT the true Fastjet library
logical(c_bool) function fastjet_available () bind(C)
  use iso_c_binding
  fastjet_available = .false.
end function fastjet_available

! Dummy implementations for FastJet C wrapper functions
subroutine fastjet_print_banner () bind (C)
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine fastjet_print_banner

function new_pseudojet (px, py, pz, e) bind (C) result (j)
  use iso_c_binding
  real(c_double), intent(in), value :: px, py, pz, e
  type(c_ptr) :: j
  j = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function new_pseudojet

subroutine pseudojet_delete (j) bind (C)
  use iso_c_binding
  type(c_ptr), value :: j
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pseudojet_delete

function pseudojet_get_e (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_e

function pseudojet_get_px (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_px

function pseudojet_get_py (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_py

function pseudojet_get_pz (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_pz

function pseudojet_get_perp (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_perp

function pseudojet_get_rap (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_rap

function pseudojet_get_phi (j) bind (C) result (p)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  real(c_double) :: p
  p = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_phi

function pseudojet_get_constituents (j) bind (C) result (cv)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j
  type(c_ptr) :: cv
  cv = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_get_constituents

function pseudojet_contains (j, p) bind (C) result (flag)
  use iso_c_binding
  type(c_ptr), intent(in), value :: j, p
  logical(c_bool) :: flag
  flag = .false.
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_contains

function new_pseudojet_vector (j, n) bind (C) result (jv)
  use iso_c_binding
  type(c_ptr), dimension(*), intent(in) :: j
  integer(c_int), intent(in), value :: n
  type(c_ptr) :: jv
  jv = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function new_pseudojet_vector

subroutine pseudojet_vector_delete (jv) bind (C)
  use iso_c_binding
  type(c_ptr), value :: jv
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine pseudojet_vector_delete

function pseudojet_vector_get_size (jv) bind (C) result (n)
  use iso_c_binding
  type(c_ptr), intent(in), value :: jv
  integer(c_int) :: n
  n = 0
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_vector_get_size

function pseudojet_vector_get_jet (jv, i) bind (C) result (cptr)
  use iso_c_binding
  type(c_ptr), intent(in), value :: jv
  integer(c_int), intent(in), value :: i
  type(c_ptr) :: cptr
  cptr = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_vector_get_jet

function pseudojet_vector_sorted_by_pt (jets) bind (C) result (sorted_jets)
  use iso_c_binding
  type(c_ptr), intent(in), value :: jets
  type(c_ptr) :: sorted_jets
  sorted_jets = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function pseudojet_vector_sorted_by_pt

function new_jet_definition (jet_alg, r, jet_ycut) bind (C) result (jet_def)
  use iso_c_binding
  integer(c_int), intent(in), value :: jet_alg
  real(c_double), intent(in), value :: jet_ycut
  real(c_double), intent(in), value :: r
  type(c_ptr) :: jet_def
  jet_def = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function new_jet_definition

subroutine jet_definition_delete (jet_def) bind (C)
  use iso_c_binding
  type(c_ptr), value :: jet_def
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine jet_definition_delete

function jet_definition_get_description (jet_def) bind (C) result (str)
  use iso_c_binding
  type(c_ptr), intent(in), value :: jet_def
  type(c_ptr) :: str
  str = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function jet_definition_get_description

function new_cluster_sequence (jv, jet_def) bind (C) result (cs)
  use iso_c_binding
  type(c_ptr), intent(in), value :: jv
  type(c_ptr), intent(in), value :: jet_def
  type(c_ptr) :: cs
  cs = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function new_cluster_sequence

subroutine cluster_sequence_delete (cs) bind (C)
  use iso_c_binding
  type(c_ptr), value :: cs
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine cluster_sequence_delete

function cluster_sequence_get_inclusive_jets (cs) bind (C) result (jets)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cs
  type(c_ptr) :: jets
  jets = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function cluster_sequence_get_inclusive_jets

function cluster_sequence_get_exclusive_jets (cs, dcut) bind (C) result (jets)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cs
  real(c_double), intent(in) :: dcut
  type(c_ptr) :: jets
  jets = c_null_ptr
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function cluster_sequence_get_exclusive_jets

function cluster_sequence_get_jet_indices (cs, jv) bind (C) result (idx)
  use iso_c_binding
  type(c_ptr), intent(in), value :: cs, jv
  type(c_ptr) :: idx
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function cluster_sequence_get_jet_indices

function int_vector_get (iv, i) bind (C) result (j)
  use iso_c_binding
  type(c_ptr), intent(in), value :: iv
  integer(c_int), intent(in), value :: i
  integer(c_int) :: j
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function int_vector_get

subroutine int_vector_delete (iv) bind (C)
  use iso_c_binding
  type(c_ptr), value :: iv
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** FastJet: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine int_vector_delete
