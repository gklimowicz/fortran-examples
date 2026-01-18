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

submodule (xml) xml_s

  use io_units
  use system_defs, only: BLANK, TAB
  use diagnostics
  use ifiles

  implicit none

contains

  module subroutine cstream_final (stream)
    class(cstream_t), intent(inout) :: stream
    stream%cache_is_empty = .true.
    call stream%stream_t%final ()
  end subroutine cstream_final

  module subroutine cstream_get_record (cstream, string, iostat)
    class(cstream_t), intent(inout) :: cstream
    type(string_t), intent(out) :: string
    integer, intent(out) :: iostat
    if (cstream%cache_is_empty) then
       call stream_get_record (cstream%stream_t, string, iostat)
    else
       string = cstream%cache
       cstream%cache_is_empty = .true.
       iostat = 0
    end if
  end subroutine cstream_get_record

  module subroutine cstream_revert_record (cstream, string)
    class(cstream_t), intent(inout) :: cstream
    type(string_t), intent(in) :: string
    if (cstream%cache_is_empty) then
       cstream%cache = string
       cstream%cache_is_empty = .false.
    else
       call msg_bug ("CStream: attempt to revert twice")
    end if
  end subroutine cstream_revert_record

  module subroutine attribute_write (object, unit)
    class(attribute_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A,'=')", advance = "no")  char (object%name)
    if (object%known) then
       write (u, "(A,A,A)", advance = "no")  '"', char (object%value), '"'
    else
       write (u, "('?')", advance = "no")
    end if
  end subroutine attribute_write

  module function xml_attribute (name, value) result (attribute)
    type(string_t), intent(in) :: name
    type(string_t), intent(in), optional :: value
    type(attribute_t) :: attribute
    attribute%name = name
    if (present (value)) then
       attribute%value = value
       attribute%known = .true.
    else
       attribute%known = .false.
    end if
  end function xml_attribute

  module subroutine attribute_set_value (attribute, value)
    class(attribute_t), intent(inout) :: attribute
    type(string_t), intent(in) :: value
    attribute%value = value
    attribute%known = .true.
  end subroutine attribute_set_value

  module function attribute_get_value (attribute) result (value)
    class(attribute_t), intent(in) :: attribute
    type(string_t) :: value
    if (attribute%known) then
       value = attribute%value
    else
       value = "?"
    end if
  end function attribute_get_value

  module subroutine xml_tag_init_no_attributes (tag, name, has_content)
    class(xml_tag_t), intent(out) :: tag
    type(string_t), intent(in) :: name
    logical, intent(in), optional :: has_content
    tag%name = name
    allocate (tag%attribute (0))
    if (present (has_content))  tag%has_content = has_content
  end subroutine xml_tag_init_no_attributes

  module subroutine xml_tag_init_with_attributes (tag, name, attribute, has_content)
    class(xml_tag_t), intent(out) :: tag
    type(string_t), intent(in) :: name
    type(attribute_t), dimension(:), intent(in) :: attribute
    logical, intent(in), optional :: has_content
    tag%name = name
    allocate (tag%attribute (size (attribute)))
    tag%attribute = attribute
    if (present (has_content))  tag%has_content = has_content
  end subroutine xml_tag_init_with_attributes

  module subroutine xml_tag_set_attribute (tag, i, value)
    class(xml_tag_t), intent(inout) :: tag
    integer, intent(in) :: i
    type(string_t), intent(in) :: value
    call tag%attribute(i)%set_value (value)
  end subroutine xml_tag_set_attribute

  module function xml_tag_get_attribute (tag, i) result (value)
    class(xml_tag_t), intent(in) :: tag
    integer, intent(in) :: i
    type(string_t) :: value
    value = tag%attribute(i)%get_value ()
  end function xml_tag_get_attribute

  module subroutine xml_tag_write (tag, unit)
    class(xml_tag_t), intent(in) :: tag
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "('<',A)", advance = "no")  char (tag%name)
    do i = 1, size (tag%attribute)
       write (u, "(1x)", advance = "no")
       call tag%attribute(i)%write (u)
    end do
    if (tag%has_content) then
       write (u, "('>')", advance = "no")
    else
       write (u, "(' />')", advance = "no")
    end if
  end subroutine xml_tag_write

  module subroutine xml_tag_close (tag, unit)
    class(xml_tag_t), intent(in) :: tag
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "('</',A,'>')", advance = "no")  char (tag%name)
  end subroutine xml_tag_close

  module subroutine xml_tag_write_with_content (tag, content, unit)
    class(xml_tag_t), intent(in) :: tag
    type(string_t), intent(in) :: content
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    call tag%write (u)
    write (u, "(A)", advance = "no")  char (content)
    call tag%close (u)
  end subroutine xml_tag_write_with_content

  module subroutine xml_tag_read (tag, cstream, success)
    class(xml_tag_t), intent(inout) :: tag
    type(cstream_t), intent(inout) :: cstream
    logical, intent(out) :: success
    type(string_t) :: string
    integer :: iostat, p1, p2
    character(2), parameter :: WS = BLANK // TAB
    logical :: done

    ! Skip comments and blank lines
    FIND_NON_COMMENT: do
       FIND_NONEMPTY_RECORD: do
          call cstream%get_record (string, iostat)
          if (iostat /= 0)  call err_io ()
          p1 = verify (string, WS)
          if (p1 > 0)  exit FIND_NONEMPTY_RECORD
       end do FIND_NONEMPTY_RECORD

       ! Look for comment beginning
       p2 = p1 + 3
       if (extract (string, p1, p2) /= "<!--")  exit FIND_NON_COMMENT

       ! Look for comment end, then restart
       string = extract (string, p2 + 1)
       FIND_COMMENT_END: do
          do p1 = 1, len (string) - 2
             p2 = p1 + 2
             if (extract (string, p1, p2) == "-->") then

                ! Return trailing text to the stream
                string = extract (string, p2 + 1)
                if (string /= "")  call cstream%revert_record (string)
                exit FIND_COMMENT_END

             end if
          end do
          call cstream%get_record (string, iostat)
          if (iostat /= 0)  call err_io ()
       end do FIND_COMMENT_END
    end do FIND_NON_COMMENT

    ! Look for opening <
    p2 = p1
    if (extract (string, p1, p2) /= "<") then
       call cstream%revert_record (string)
       success = .false.;  return
    else

       ! Look for tag name
       string = extract (string, p2 + 1)
       p1 = verify (string, WS);  if (p1 == 0)  call err_incomplete ()
       p2 = p1 + len (tag%name) - 1
       if (extract (string, p1, p2) /= tag%name) then
          call cstream%revert_record ("<" // string)
          success = .false.;  return
       else

          ! Look for attributes
          string = extract (string, p2 + 1)
          READ_ATTRIBUTES: do
             call tag%read_attribute (string, done)
             if (done)  exit READ_ATTRIBUTES
          end do READ_ATTRIBUTES

          ! Look for closing >
          p1 = verify (string, WS);  if (p1 == 0)  call err_incomplete ()
          p2 = p1
          if (extract (string, p1, p1) == ">") then
             tag%has_content = .true.
          else

             ! Look for closing />
             p2 = p1 + 1
             if (extract (string, p1, p2) /= "/>")  call err_incomplete ()
          end if

          ! Return trailing text to the stream
          string = extract (string, p2 + 1)
          if (string /= "")  call cstream%revert_record (string)
          success = .true.

       end if
    end if

  contains

    subroutine err_io ()
      select case (iostat)
      case (:-1)
         call msg_fatal ("XML: Error reading tag '" // char (tag%name) &
              // "': end of file")
      case (1:)
         call msg_fatal ("XML: Error reading tag '" // char (tag%name) &
              // "': I/O error")
      end select
      success = .false.
    end subroutine err_io

    subroutine err_incomplete ()
      call msg_fatal ("XML: Error reading tag '" // char (tag%name) &
           // "': tag incomplete")
      success = .false.
    end subroutine err_incomplete

  end subroutine xml_tag_read

  module subroutine xml_tag_read_attribute (tag, string, done)
    class(xml_tag_t), intent(inout) :: tag
    type(string_t), intent(inout) :: string
    logical, intent(out) :: done
    character(2), parameter :: WS = BLANK // TAB
    type(string_t) :: name, value
    integer :: p1, p2, i

    p1 = verify (string, WS);  if (p1 == 0)  call err ()
    p2 = p1

    ! Look for first terminating '>' or '/>'
    if (extract (string, p1, p2) == ">") then
       done = .true.
    else
       p2 = p1 + 1
       if (extract (string, p1, p2) == "/>") then
          done = .true.
       else

          ! Look for '='
          p2 = scan (string, '=')
          if (p2 == 0)  call err ()
          name = trim (extract (string, p1, p2 - 1))

          ! Look for '"'
          string = extract (string, p2 + 1)
          p1 = verify (string, WS);  if (p1 == 0)  call err ()
          p2 = p1
          if (extract (string, p1, p2) /= '"')  call err ()

          ! Look for matching '"' and get value
          string = extract (string, p2 + 1)
          p1 = 1
          p2 = scan (string, '"')
          if (p2 == 0)  call err ()
          value = extract (string, p1, p2 - 1)

          SCAN_KNOWN_ATTRIBUTES: do i = 1, size (tag%attribute)
             if (name == tag%attribute(i)%name) then
                call tag%attribute(i)%set_value (value)
                exit SCAN_KNOWN_ATTRIBUTES
             end if
          end do SCAN_KNOWN_ATTRIBUTES

          string = extract (string, p2 + 1)
          done = .false.
       end if
    end if

  contains

    subroutine err ()
      call msg_fatal ("XML: Error reading attributes of '" // char (tag%name) &
           // "': syntax error")
    end subroutine err

  end subroutine xml_tag_read_attribute

  module subroutine xml_tag_read_content (tag, cstream, content, closing)
    class(xml_tag_t), intent(in) :: tag
    type(cstream_t), intent(inout) :: cstream
    type(string_t), intent(out) :: content
    type(string_t) :: string
    logical, intent(out) :: closing
    integer :: iostat
    integer :: p0, p1, p2
    character(2), parameter :: WS = BLANK // TAB
    call cstream%get_record (content, iostat)
    if (iostat /= 0)  call err_io ()
    closing = .false.
    FIND_CLOSING: do p0 = 1, len (content) - 1

       ! Look for terminating </
       p1 = p0
       p2 = p1 + 1
       if (extract (content, p1, p2) == "</") then

          ! Look for closing tag name
          string = extract (content, p2 + 1)
          p1 = verify (string, WS);  if (p1 == 0)  call err_incomplete ()
          p2 = p1 + len (tag%name) - 1
          if (extract (string, p1, p2) == tag%name) then

             ! Tag name matches: look for final >
             string = extract (string, p2 + 1)
             p1 = verify (string, WS);  if (p1 == 0)  call err_incomplete ()
             p2 = p1
             if (extract (string, p1, p2) /= ">")  call err_incomplete ()

             ! Return trailing text to the stream
             string = extract (string, p2 + 1)
             if (string /= "")  call cstream%revert_record (string)
             content = extract (content, 1, p0 -1)
             closing = .true.
             exit FIND_CLOSING

          end if
       end if
    end do FIND_CLOSING

  contains

    subroutine err_io ()
      select case (iostat)
      case (:-1)
         call msg_fatal ("XML: Error reading content of '" // char (tag%name) &
              // "': end of file")
      case (1:)
         call msg_fatal ("XML: Error reading content of '" // char (tag%name) &
              // "': I/O error")
      end select
      closing = .false.
    end subroutine err_io

    subroutine err_incomplete ()
      call msg_fatal ("XML: Error reading content '" // char (tag%name) &
           // "': closing tag incomplete")
      closing = .false.
    end subroutine err_incomplete

  end subroutine xml_tag_read_content


end submodule xml_s

