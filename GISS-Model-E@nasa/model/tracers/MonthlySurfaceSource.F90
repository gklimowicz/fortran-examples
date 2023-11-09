module MonthlySurfaceSource_mod
!!$  implicit none
!!$  private
!!$
!!$  public :: MonthlySurfaceSource
!!$
!!$  type MonthlySurfaceSource
!!$    character(len=?) :: fileName
!!$    real*8, allocatable :: bufferA(:,:)
!!$    real*8, allocatable :: bufferB(:,:)
!!$  contains
!!$    procedure :: get
!!$    procedure :: open
!!$    procedure :: close
!!$    procedure :: read
!!$  end type MonthlySurfaceSource
!!$
!!$contains
!!$
!!$  
!!$
end module MonthlySurfaceSource_mod
