module MathematicalConstants_mod
   implicit none
   public

  real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi
  real*8,parameter :: TWOPI = 2d0*PI           !@param twopi 2*pi
  real*8,parameter :: RADIAN = PI/180d0        !@param radian pi/180
!@param zero,one 0 and 1 for occasional use as arguments
  real*8,parameter :: zero = 0d0, one=1d0
!@param rt2,byrt2   sqrt(2), 1/sqrt(2)
  real*8,parameter :: rt2 = 1.4142135623730950d0
  real*8,parameter :: byrt2 = 1./rt2
!@param rt3,byrt3   sqrt(3), 1/sqrt(3)
  real*8,parameter :: rt3 = 1.7320508075688772d0
  real*8,parameter :: byrt3 = 1./rt3
!@param rt12,byrt12   sqrt(12), 1/sqrt(12)
  real*8,parameter :: rt12 = 3.4641016151377546d0
  real*8,parameter :: byrt12 = 1./rt12
  real*8,parameter :: by3 =1./3d0  !@param by3  1/3
  real*8,parameter :: by6 =1./6d0  !@param by6  1/6
  real*8,parameter :: by9 =1./9d0  !@param by9  1/9
  real*8,parameter :: by12=1./12d0 !@param by12 1/12

end module MATHEMATICALCONSTANTS_MOD
