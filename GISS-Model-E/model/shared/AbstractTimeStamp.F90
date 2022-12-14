module AbstractTimeStamp_mod
  use StringUtilities_mod, only: toLowerCase
  implicit none
  private

  public :: AbstractTimeStamp

  type, abstract :: AbstractTimeStamp
  contains
    procedure(toString), deferred :: toString
    procedure :: print_unit
    procedure :: print_stdout
    generic :: print => print_unit, print_stdout
  end type AbstractTimeStamp

  abstract interface

    subroutine toString(this, string)
      import AbstractTimeStamp
      class (AbstractTimeStamp), intent(in) :: this
      character(len=*), intent(out) :: string
    end subroutine toString

  end interface

contains


  subroutine print_stdout(this) 
    use iso_fortran_env, only: OUTPUT_UNIT
    class (AbstractTimeStamp), intent(in) :: this
    call this%print(OUTPUT_UNIT)
  end subroutine print_stdout


  subroutine print_unit(this, unit) 
    class (AbstractTimeStamp), intent(in) :: this
    integer, intent(in) :: unit
    
    character(len=80) :: string

    call this%toString(string)
    write(unit,'(a)') trim(string)

  end subroutine print_unit


end module AbstractTimeStamp_mod


#define TYPE AbstractTimeStamp
#define USE_MODULE AbstractTimeStamp_mod
#define TYPE_NAME AbstractTimeStamp
#define HAS_PRINT

#include "AssociativeArrayTemplate.h"

#define VALUE_TYPE AbstractTimeStamp
#define ASSOCIATIVE_ARRAY_TYPE AbstractTimeStampAssociativeArray
#undef ITERATOR_TYPE
#define HASH_TYPE AbstractTimeStampHashMap

#include "HashMapTemplate.h"
