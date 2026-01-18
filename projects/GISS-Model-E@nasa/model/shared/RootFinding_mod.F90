!------------------------------------------------------------------------------
Module RootFinding_mod
!------------------------------------------------------------------------------
!@sum Module that contains methods to find roots of a function
!@auth SSSO ASTG
  implicit none
  private
  public NewtonMethod

  integer, parameter :: DP = selected_real_kind(14)
  INTEGER, PARAMETER :: maxNumIterations=100

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  function NewtonMethod(func, a, b, accuracy, ierr) result(root)
!@sum Use Newton-Raphson + bisection to solve Func(x)=0 bracketed
!@+ between a and b
    REAL(kind=DP), INTENT(IN) :: a, b, accuracy
    integer, intent(out) :: ierr
    real(kind=DP) :: root 
    interface
      SUBROUTINE func(x, fval, fderiv)
        REAL(8), INTENT(IN) :: x
        REAL(8), INTENT(OUT) :: fval, fderiv
      end SUBROUTINE func
    end interface
!
    integer :: j
    REAL(kind=DP) :: df, dx, dxold, f
    REAL(kind=DP) :: fb, fa, temp, xb, xa

    ierr = 0
    call func(a, fa, df)
    call func(b, fb, df)
    if (fa*fb > 0.0d0) then
      ierr = 1
      print*, 'Error: root must be bracketed in Newton-Raphson'
    end if
    if (fa == 0.0d0) then
      root = a
      return
    else if (fb == 0.0d0) then
      root = b
      return
    else if (fa < 0.0d0) then
      xa = a
      xb = b
    else
      xb = a
      xa = b
    end if
    root = .5*(a+b)   ! initial guess
    dxold = abs(b-a)
    dx = dxold
    call func(root, f, df)
    do j = 1, maxNumIterations
      ! use bisection 
      if ( ((root-xb)*df-f)*((root-xa)*df-f) > 0.0 .or. &
        abs(2.0*f) > abs(dxold*df) ) then
        dxold = dx
        dx = 0.5*(xb-xa)
        root = xa + dx
        if (xa == root) return
      else  ! take a Newton step
        dxold = dx
        dx = f/df
        temp = root
        root = root - dx
        if (temp == root) return
      end if
      if (abs(dx) < accuracy) return
      call func(root, f, df)
      if (f < 0.0d0) then
        xa = root
      else
        xb = root
      end if
    end do
    ierr = 1
    print*, 'Error: Newton-Raphson exceeding maximum iterations'

  END function NewtonMethod

end Module RootFinding_mod
