!---------------------------------------------------------------------
!
! AnniversaryDate is a simple data structure that contains integer
! month and date to represent recurring events such as solstices and
! equinoctes.  It is an extension of the (empty) AbstractTimeStamp
! class, which provides print() methods.
!
! The structure is merely provided as a convenience for calendar
! diagnostics.
!
!---------------------------------------------------------------------

module AnniversaryDate_mod
  use AbstractTimeStamp_mod
  implicit none
  private

  public :: AnniversaryDate

  type, extends(AbstractTimeStamp) :: AnniversaryDate
    integer :: month
    integer :: date
  contains
    procedure :: toString
  end type AnniversaryDate

contains


  ! Convert data to human-readable string
  subroutine toString(this, string)
    class (AnniversaryDate), intent(in) :: this
    character(len=*), intent(out) :: string

    write(string,'(i2.2,"-",i2.2)') this%month, this%date

  end subroutine toString


end module AnniversaryDate_mod
