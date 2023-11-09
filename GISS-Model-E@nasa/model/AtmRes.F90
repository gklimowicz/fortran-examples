!@sum Resolution file wrapper
!@auth SSSO Development Team
module Resolution
  use HorizontalRes
  use VerticalRes
  implicit none
  public      
end module Resolution

module threeD_mass_unfinished
! legacy parameters that still appear in a couple of places in the model
  use verticalres, only : plbot,ls1_nominal
  real*8, parameter :: psf=plbot(1),ptop=plbot(ls1_nominal),psfmpt=psf-ptop
end module threeD_mass_unfinished
