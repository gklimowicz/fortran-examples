subroutine dekay (kto, hx)
  integer, intent(in) :: kto
  double precision, dimension(4), intent(in) :: hx
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine dekay

subroutine dexay (kto, pol)
  integer, intent(in) :: kto
  double precision, dimension(4), intent(in) :: pol
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine dexay

subroutine initdk (mode, keypol)
  integer, intent(in) :: mode, keypol
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine initdk

subroutine inimas (mode, keypol)
  integer, intent(in) :: mode, keypol
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine inimas

subroutine iniphx (xk00)
  double precision, intent(in) :: xk00
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine iniphx

subroutine inietc (jakk1, jakk2, itd, ifpho)
  integer, intent(in) :: jakk1, jakk2, itd, ifpho
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine inietc

subroutine phoini ()
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine phoini

function wthiggs (ifpseudo, hh1, hh2)
  double precision, dimension(4), intent(in) :: hh1, hh2
  logical, intent(in) :: ifpseudo
  double precision :: wthiggs
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end function wthiggs
  
subroutine taupi0 (mode, jak, ion)
  integer, intent(in) :: mode, jak
  integer, dimension(3), intent(in) :: ion
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine taupi0

subroutine photos (id)
  integer, intent(in) :: id
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine photos

subroutine ranmar (rvec, lenv)
  double precision, dimension(lenv) :: rvec
  integer, intent(in) :: lenv
  write (0, "(A)")  "**************************************************************"
  write (0, "(A)")  "*** Error: TAUOLA has not been enabled, WHIZARD terminates ***"
  write (0, "(A)")  "**************************************************************"
  stop
end subroutine ranmar
