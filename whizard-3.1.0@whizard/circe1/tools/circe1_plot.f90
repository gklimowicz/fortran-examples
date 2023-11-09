program circe1_plot
  use kinds
  use circe1 

  implicit none

  real(kind=double) :: xmin, xmax, y, roots
  integer :: xory, nstep, p1, p2, acc, ver, rev, i
  real(kind=double) :: x, logx, d
  read *, xory, xmin, xmax, nstep, y, p1, p2, roots, acc, ver, rev
  call circes (0d0, 0d0, roots, acc, ver, rev, 0)
  do i = 0, nstep
     logx = log (xmin) + i * log (xmax/xmin) / nstep
     x = exp (logx)
     d = 0d0
     if (xory .eq. 1) then
        if (p1 .eq. C1_PHOTON) then
           d = circe (x, y, p1, p2)
        else
           d = circe (1d0 - x, y, p1, p2)
        end if
     else if (xory .eq. 2) then
        if (p1 .eq. C1_PHOTON) then
           d = circe (y, x, p1, p2)
        else
           d = circe (y, 1d0 - x, p1, p2)
        end if
     end if
     if (d .gt. 1d-4) print *, x, d
  end do
end program circe1_plot
