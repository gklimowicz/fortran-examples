! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module xml_uti

  use iso_varying_string, string_t => varying_string
  use io_units

  use xml

  implicit none
  private

  public :: xml_1
  public :: xml_2
  public :: xml_3
  public :: xml_4

contains

  subroutine xml_1 (u)
    integer, intent(in) :: u
    type(xml_tag_t), allocatable :: tag
    integer :: u_tmp
    type(cstream_t) :: cstream
    logical :: success

    write (u, "(A)")  "* Test output: xml_1"
    write (u, "(A)")  "*   Purpose: write and read tag"
    write (u, "(A)")

    write (u, "(A)")  "* Empty tag"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, *)
    write (u, "(A)")  "* Tag with preceding blank lines"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"))
    write (u_tmp, *)
    write (u_tmp, "(A)")  "    "
    write (u_tmp, "(2x)", advance = "no")
    call tag%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, *)
    write (u, "(A)")  "* Tag with preceding comments"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"))
    write (u_tmp, "(A)")  "<!-- comment -->"
    write (u_tmp, *)
    write (u_tmp, "(A)")  "<!-- multiline"
    write (u_tmp, "(A)")  "     comment -->"
    call tag%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    close (u_tmp)
    deallocate (tag)

    call cstream%final ()

    write (u, *)
    write (u, "(A)")  "* Tag with name mismatch"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("wrongname"))
    call tag%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: xml_1"

  end subroutine xml_1

  subroutine xml_2 (u)
    integer, intent(in) :: u
    type(xml_tag_t), allocatable :: tag1, tag2
    integer :: u_tmp
    type(cstream_t) :: cstream
    logical :: success

    write (u, "(A)")  "* Test output: xml_2"
    write (u, "(A)")  "*   Purpose: handle optional tag"
    write (u, "(A)")

    write (u, "(A)")  "* Optional tag present"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag1)
    call tag1%init (var_str ("option"))
    call tag1%write (u_tmp)
    write (u_tmp, *)
    allocate (tag2)
    call tag2%init (var_str ("tagname"))
    call tag2%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag1, tag2)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag1)
    call tag1%init (var_str ("option"))
    call tag1%read (cstream, success)
    call tag1%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    write (u, *)
    allocate (tag2)
    call tag2%init (var_str ("tagname"))
    call tag2%read (cstream, success)
    call tag2%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    deallocate (tag1, tag2)
    close (u_tmp)
    call cstream%final ()

    write (u, *)
    write (u, "(A)")  "* Optional tag absent"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag2)
    call tag2%init (var_str ("tagname"))
    call tag2%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag2)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag1)
    call tag1%init (var_str ("option"))
    call tag1%read (cstream, success)
    call tag1%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    write (u, *)
    allocate (tag2)
    call tag2%init (var_str ("tagname"))
    call tag2%read (cstream, success)
    call tag2%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    deallocate (tag1, tag2)
    close (u_tmp)
    call cstream%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: xml_2"

  end subroutine xml_2

  subroutine xml_3 (u)
    integer, intent(in) :: u
    type(xml_tag_t), allocatable :: tag
    integer :: u_tmp
    type(cstream_t) :: cstream
    logical :: success, closing
    type(string_t) :: content

    write (u, "(A)")  "* Test output: xml_3"
    write (u, "(A)")  "*   Purpose: handle tag with content"
    write (u, "(A)")

    write (u, "(A)")  "* Tag without content"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%write (u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    write (u, "(A,L1)")  "content = ", tag%has_content
    write (u, *)
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, "(A)")  "* Tag with content"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"), has_content = .true.)
    call tag%write (var_str ("Content text"), u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%read_content (cstream, content, closing)
    call tag%write (u)
    write (u, "(A)", advance = "no")  char (content)
    call tag%close (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    write (u, "(A,L1)")  "content = ", tag%has_content
    write (u, "(A,L1)")  "closing = ", closing
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, *)
    write (u, "(A)")  "* Tag with multiline content"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"), has_content = .true.)
    call tag%write (u_tmp)
    write (u_tmp, *)
    write (u_tmp, "(A)")  "Line 1"
    write (u_tmp, "(A)")  "Line 2"
    call tag%close (u_tmp)
    write (u_tmp, *)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"))
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    do
       call tag%read_content (cstream, content, closing)
       if (closing)  exit
       write (u, "(A)")  char (content)
    end do
    call tag%close (u)
    write (u, *)
    write (u, "(A,L1)")  "success = ", success
    write (u, "(A,L1)")  "content = ", tag%has_content
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: xml_3"

  end subroutine xml_3

  subroutine xml_4 (u)
    integer, intent(in) :: u
    type(xml_tag_t), allocatable :: tag
    integer :: u_tmp
    type(cstream_t) :: cstream
    logical :: success

    write (u, "(A)")  "* Test output: xml_4"
    write (u, "(A)")  "*   Purpose: handle tag with attributes"
    write (u, "(A)")

    write (u, "(A)")  "* Tag with one mandatory and one optional attribute,"
    write (u, "(A)")  "* unknown attribute ignored"
    write (u, *)

    u_tmp = free_unit ()
    open (u_tmp, status = "scratch", action = "readwrite")

    allocate (tag)
    call tag%init (var_str ("tagname"), &
         [xml_attribute (var_str ("a1"), var_str ("foo")), &
          xml_attribute (var_str ("a3"), var_str ("gee"))])
    call tag%write (u_tmp)
    deallocate (tag)

    call show (u_tmp, u)

    write (u, *)
    write (u, "(A)")  "Result from read:"
    call cstream%init (u_tmp)
    allocate (tag)
    call tag%init (var_str ("tagname"), &
         [xml_attribute (var_str ("a1")), &
          xml_attribute (var_str ("a2"), var_str ("bar"))])
    call tag%read (cstream, success)
    call tag%write (u)
    write (u, *)
    deallocate (tag)
    close (u_tmp)
    call cstream%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: xml_4"

  end subroutine xml_4


  subroutine show (u_tmp, u)
    integer, intent(in) :: u_tmp, u
    character (80) :: buffer
    integer :: iostat
    write (u, "(A)")  "File content:"
    rewind (u_tmp)
    do
       read (u_tmp, "(A)", iostat = iostat)  buffer
       if (iostat /= 0)  exit
       write (u, "(A)")  trim (buffer)
    end do
    rewind (u_tmp)
  end subroutine show


end module xml_uti
