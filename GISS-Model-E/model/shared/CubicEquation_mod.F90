!@sum this file contains a solver for cubic equation

module CubicEquation_mod
  implicit none

contains

  subroutine cubicroot(a,b,c,d,x,n) 
!@sum solves cubic equation: a x^3 + b x^2 + c x + d = 0 
!@+   using Cardano formula (see Korn, Korn, Mathematical Handbook)
!@auth Igor Aleinov from solution by Cardano in
    implicit none
    real*8,intent(in) :: a,b,c,d  ! coefficients of cubic
    real*8, intent(out) :: x(:)   ! results ( 0-3 roots )
    integer, intent(out) :: n     ! number of roots
    real*8 :: x0,x1,x2
    real*8 :: a0,a1,a2,Q1,R1,D1
    real*8, parameter :: EPS0 = 1.d-8 ! 1.d-15
    real*8, parameter :: one3rd = 1.d0/3.d0
    real*8 :: arg, S, T
    complex*16 :: ST

    !print *,"cubicroot:",a,b,c,d

    if (abs(a) < (abs(b)+abs(c)+abs(d))*EPS0 ) then
      if (abs(b) < (abs(c)+abs(d))*EPS0) then
        if (abs(c) < abs(d)*EPS0) then
          write(*,*) "Internal Error in Cardano: no solution."
          stop
        endif
        x0 = -d/c
        x(1) = x0
        !write(*,*) "Cardano: returning",x0
        n = 1
      else
        !write(*,*) "What's this?"
        D1 = c*c - 4.d0*b*d

        if (D1 > 0.d0) then
          Q1 = sqrt(D1)
          x0 = (-c + Q1) / (2.d0 * b)
          x1 = (-c - Q1) / (2.d0 * b)
          !return
          n = 2
        else if (D1.eq.0.) then
          x0 = -c / (2.d0 * b)
          x1 = x0
          n = 1
        else 
          x0 = -c /(2.d0 *b)
          x1 = sqrt(-D1) / (2.d0* b)
          n = 0
        end if
      end if
      !print *,"CX1",x0,x1
      !x = max(x0,x1)
      x(1) = x0
      x(2) = x1
    else
      a2 = b/a
      a1 = c/a
      a0 = d/a
      Q1 = (3.d0 * a1 - a2*a2 ) / 9.d0
      R1 = (9.d0 * a2 * a1 - 27.d0 * a0 - 2.d0 * a2*a2*a2) /54.d0
      D1 = Q1*Q1*Q1 + R1*R1
      !write(677,*) "Q", D1
      !write(*,*) "abcda2a1a0Q1R1D1",a,b,c,d,a2,a1,a0,Q1,R1,D1
      if (D1 > 0.d0) then       !* only one real root *!
        !write(*,*) "One real root."
        arg = R1 + sqrt(D1)
        S = sign(1.d0, arg) * (abs(arg)**one3rd)
        arg = R1 - sqrt(D1)
        T = sign(1.d0, arg) * (abs(arg)**one3rd)
        x0 = -a2/3.d0 + S + T
        x1 = -a2/3.d0 - (S+T)*0.5d0
        x2 = sqrt(3.d0) * (S-T)*0.5d0
        !print *,"CX2",x0,x1,x2
        n = 1
      else if (D1.eq.0.) then !* two roots coincide * *!
        !write(*,*) "Two roots coincide."
        S = sign(1.d0, R1) * (abs(R1)**one3rd)
        x0 = -a2/3.d0 + 2.d0*S
        x1 = -a2/3.d0 - S
        x2 = x1
        !print *,"CX3",x0,x1,x2
        n =2
      else                    !* three different real roots *!
        !call CRtCube( R1, sqrt(-D1), S, T)
        !write(*,*) "Three different real roots. a2R1D1ST",a2,R1,D1,S,T
        ST = ( cmplx(R1, sqrt(-D1),kind(1.d0)) )**one3rd
        S = real (ST)
        T = aimag(ST)
        x0 = -a2/3.d0 + 2.d0*S
        x1 = -a2/3.d0 - S + sqrt(3.d0)*T
        x2 = -a2/3.d0 - S - sqrt(3.d0)*T
        !print *,"CX4",x0,x1,x2
        n = 3
      end if
      !x = max(x0,x1,x2)
      !x = x2
      x(1) = x0
      x(2) = x1
      x(3) = x2
    end if
  end subroutine cubicroot

end module CubicEquation_mod

