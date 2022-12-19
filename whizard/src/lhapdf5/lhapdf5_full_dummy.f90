! Dummy replacement routines for the case that LHAPDF is fully absent.

subroutine InitPDFsetM (set, file)
  integer, intent(in) :: set
  character(*), intent(in) :: file
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine InitPDFsetM

subroutine InitPDFM (set, mem)
  integer, intent(in) :: set, mem
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine InitPDFM

subroutine numberPDFM (set, n_members)
  integer, intent(in) :: set
  integer, intent(out) :: n_members
  n_members = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine numberPDFM

subroutine evolvePDF (x, q, ff)
  double precision, intent(in) :: x, q
  double precision, dimension(-6:6), intent(out) :: ff
  ff = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine evolvePDF

subroutine evolvePDFM (set, x, q, ff)
  integer, intent(in) :: set
  double precision, intent(in) :: x, q
  double precision, dimension(-6:6), intent(out) :: ff
  ff = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine evolvePDFM

subroutine evolvePDFpM (set, x, q, s, scheme, ff)
  integer, intent(in) :: set
  double precision, intent(in) :: x, q, s
  integer, intent(in) :: scheme
  double precision, dimension(-6:6), intent(out) :: ff
  ff = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine evolvePDFpM

subroutine GetXminM (set, mem, xmin)
  integer, intent(in) :: set, mem
  double precision, intent(out) :: xmin
  xmin = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetXminM

subroutine GetXmin (mem, xmin)
  integer, intent(in) :: mem
  double precision, intent(out) :: xmin
  xmin = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetXmin

subroutine GetXmaxM (set, mem, xmax)
  integer, intent(in) :: set, mem
  double precision, intent(out) :: xmax
  xmax = 1
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetXmaxM

subroutine GetXmax (mem, xmax)
  integer, intent(in) :: mem
  double precision, intent(out) :: xmax
  xmax = 1
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetXmax

subroutine GetQ2minM (set, mem, q2min)
  integer, intent(in) :: set, mem
  double precision, intent(out) :: q2min
  q2min = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetQ2minM

subroutine GetQ2min (mem, q2min)
  integer, intent(in) :: mem
  double precision, intent(out) :: q2min
  q2min = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetQ2min

subroutine GetQ2maxM (set, mem, q2max)
  integer, intent(in) :: set, mem
  double precision, intent(out) :: q2max
  q2max = huge (1.d0)
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetQ2maxM

subroutine GetQ2max (mem, q2max)
  integer, intent(in) :: mem
  double precision, intent(out) :: q2max
  q2max = huge (1.d0)
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetQ2max

double precision function alphasPDF (Q)
  double precision, intent(in) :: Q
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end function alphasPDF

subroutine evolvePDFphoton (x, q, ff, fphot)
  double precision, intent(in) :: x, q
  double precision, dimension(-6:6), intent(out) :: ff
  double precision, intent(out) :: fphot
  ff = 0
  fphot = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine evolvePDFphoton

subroutine evolvePDFphotonM (set, x, q, ff, fphot)
  integer, intent(in) :: set
  double precision, intent(in) :: x, q
  double precision, dimension(-6:6), intent(out) :: ff
  double precision, intent(out) :: fphot
  ff = 0
  fphot = 0
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine evolvePDFphotonM

subroutine GetQmass (nf, mass)
  integer, intent(in) :: nf
  double precision, intent(out) :: mass
  mass = 0 
  write (0, "(A)")  "*************************************************************"
  write (0, "(A)")  "*** LHAPDF: Error: library not linked, WHIZARD terminates ***"
  write (0, "(A)")  "*************************************************************"
  stop
end subroutine GetQmass
