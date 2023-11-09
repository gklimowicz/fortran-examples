!@sum this file will contain some system calls which are trivial to
!@+   implement in C, but are missing or are too awkward in Fortran

module SystemTools
  !use iso_c_binding, only : c_ptr, c_int

  implicit none
  private

  public stFileList, stLinkStatus

  interface
     function c_field_list(dir,list,nl,ls) &
          bind(C,name='c_field_list')
       use iso_c_binding, only : c_ptr, c_int
       integer (c_int) :: c_field_list
       type (c_ptr), value :: dir
       type (c_ptr), value :: list
       integer (c_int), value :: nl
       integer (c_int), value :: ls
     end function c_field_list
  end interface

  interface
     function c_link_status(path) &
          bind(C,name='c_link_status')
       use iso_c_binding, only : c_ptr, c_int
       integer (c_int) :: c_link_status
       type (c_ptr), value :: path
     end function c_link_status
  end interface

contains

  subroutine stFileList(dir,list,num)
    use iso_c_binding, only : c_ptr, c_int, C_NULL_CHAR, c_loc
    character(len=*), intent(in) :: dir
    character(len=*), dimension(:), intent(out), target :: list
    integer, intent(out) :: num
    !---
    integer, parameter :: MAX_LEN=80
    character(len=len(dir)+1), target :: dir0
    !character(len=MAX_LEN), dimension(size(list)), target :: list_loc
    !character(len=:), dimension(:), allocatable, target :: list_loc
    !character(len=1), dimension(:,:), allocatable, target :: list_loc
    character(len=len(list(1))+1), dimension(size(list)), target :: list_loc
!    character(len=1), dimension(size(list)), target :: list_loc
    type(c_ptr), dimension(size(list)), target :: c_list_ptrs
    integer (c_int) :: nl
    integer (c_int) :: ls
    !integer (c_int) :: retcode
    integer :: n
    !character(len=len(list(1))+1) :: tmp_str

    nl = size(list)
    ls = len(list(1))

    !if ( nl > MAX_LEN ) call stop_model("stFileList: increase MAX_LEN",255)
    !allocate( list_loc(ls+1,nl) )

    dir0 = trim(dir)//C_NULL_CHAR

    do n = 1, nl
      c_list_ptrs(n) = c_loc(list_loc(n)(1:1))
      !c_list_ptrs(n) = c_loc(list_loc(1,n))
    end do

    num=c_field_list(c_loc(dir0(1:1)), c_loc(c_list_ptrs), nl, ls)

    !write(0,*) "stFileList: num=" , num

    select case(num)
    case(-1)
      call stop_model("stFileList: could not open directory",255)
    case(-2)
      call stop_model("stFileList: file names too long",255)
    case(-3)
      call stop_model("stFileList: list too long",255)
    case(0:)
      do n = 1, num
        list(n) = list_loc(n)(:scan(list_loc(n),C_NULL_CHAR)-1)
        !tmp_str = transfer( list_loc(:,n), tmp_str)
        !list(n) = tmp_str(:scan(tmp_str,C_NULL_CHAR)-1)
      end do
    case default
      call stop_model("stFileList: unknown error",255)
    end select

  end subroutine stFileList


  subroutine stLinkStatus(path, status)
!@sum returns the status of the file specified by "path" (if "path" is a
!@+ symbolic link it is followed to the actual file):
!@+ 1 - regular file
!@+ 2 - directory
!@+ 0 - file exists, but is nether regular file or directory
!@+ -1 - file doesn't exist or reading error
    use iso_c_binding, only : C_NULL_CHAR, c_loc
    character*(*) :: path
    integer status
    !---
    character(len=len(path)+1), target :: path0

    path0 = trim(path)//C_NULL_CHAR
    status = c_link_status( c_loc(path0(1:1)) )

  end subroutine stLinkStatus

end module SystemTools
