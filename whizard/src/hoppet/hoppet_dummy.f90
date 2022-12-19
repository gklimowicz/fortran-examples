! Dummy replacement routines

subroutine InitForWhizard (pdf_builtin, pdf, pdf_id)
  use lhapdf
  logical, intent(in) :: pdf_builtin
  type(lhapdf_pdf_t), intent(inout), optional :: pdf
  integer, intent(in), optional :: pdf_id
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** HOPPET: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine InitForWhizard

subroutine EvalForWhizard (x, q, f)
  double precision, intent(in)  :: x, q
  double precision, intent(out) :: f(-6:6)
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** HOPPET: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine EvalForWhizard
