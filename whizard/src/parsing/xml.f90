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

module xml

  use iso_varying_string, string_t => varying_string
  use lexers

  implicit none
  private

  public :: cstream_t
  public :: xml_attribute
  public :: xml_tag_t

  type, extends (stream_t) :: cstream_t
     logical :: cache_is_empty = .true.
     type(string_t) :: cache
   contains
     procedure :: final => cstream_final
     procedure :: get_record => cstream_get_record
     procedure :: revert_record => cstream_revert_record
  end type cstream_t

  type :: attribute_t
     type(string_t) :: name
     type(string_t) :: value
     logical :: known = .false.
   contains
     procedure :: write => attribute_write
     procedure :: set_value => attribute_set_value
     procedure :: get_value => attribute_get_value
  end type attribute_t

  type :: xml_tag_t
     type(string_t) :: name
     type(attribute_t), dimension(:), allocatable :: attribute
     logical :: has_content = .false.
   contains
     generic :: init => init_no_attributes
     procedure :: init_no_attributes => xml_tag_init_no_attributes
     generic :: init => init_with_attributes
     procedure :: init_with_attributes => xml_tag_init_with_attributes
     procedure :: set_attribute => xml_tag_set_attribute
     procedure :: get_attribute => xml_tag_get_attribute
     generic :: write => write_without_content
     procedure :: write_without_content => xml_tag_write
     procedure :: close => xml_tag_close
     generic :: write => write_with_content
     procedure :: write_with_content => xml_tag_write_with_content
     procedure :: read => xml_tag_read
     procedure :: read_attribute => xml_tag_read_attribute
     procedure :: read_content => xml_tag_read_content
  end type xml_tag_t


  interface
    module subroutine cstream_final (stream)
      class(cstream_t), intent(inout) :: stream
    end subroutine cstream_final
    module subroutine cstream_get_record (cstream, string, iostat)
      class(cstream_t), intent(inout) :: cstream
      type(string_t), intent(out) :: string
      integer, intent(out) :: iostat
    end subroutine cstream_get_record
    module subroutine cstream_revert_record (cstream, string)
      class(cstream_t), intent(inout) :: cstream
      type(string_t), intent(in) :: string
    end subroutine cstream_revert_record
    module subroutine attribute_write (object, unit)
      class(attribute_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine attribute_write
    module function xml_attribute (name, value) result (attribute)
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: value
      type(attribute_t) :: attribute
    end function xml_attribute
    module subroutine attribute_set_value (attribute, value)
      class(attribute_t), intent(inout) :: attribute
      type(string_t), intent(in) :: value
    end subroutine attribute_set_value
    module function attribute_get_value (attribute) result (value)
      class(attribute_t), intent(in) :: attribute
      type(string_t) :: value
    end function attribute_get_value
    module subroutine xml_tag_init_no_attributes (tag, name, has_content)
      class(xml_tag_t), intent(out) :: tag
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: has_content
    end subroutine xml_tag_init_no_attributes
    module subroutine xml_tag_init_with_attributes (tag, name, attribute, has_content)
      class(xml_tag_t), intent(out) :: tag
      type(string_t), intent(in) :: name
      type(attribute_t), dimension(:), intent(in) :: attribute
      logical, intent(in), optional :: has_content
    end subroutine xml_tag_init_with_attributes
    module subroutine xml_tag_set_attribute (tag, i, value)
      class(xml_tag_t), intent(inout) :: tag
      integer, intent(in) :: i
      type(string_t), intent(in) :: value
    end subroutine xml_tag_set_attribute
    module function xml_tag_get_attribute (tag, i) result (value)
      class(xml_tag_t), intent(in) :: tag
      integer, intent(in) :: i
      type(string_t) :: value
    end function xml_tag_get_attribute
    module subroutine xml_tag_write (tag, unit)
      class(xml_tag_t), intent(in) :: tag
      integer, intent(in), optional :: unit
    end subroutine xml_tag_write
    module subroutine xml_tag_close (tag, unit)
      class(xml_tag_t), intent(in) :: tag
      integer, intent(in), optional :: unit
    end subroutine xml_tag_close
    module subroutine xml_tag_write_with_content (tag, content, unit)
      class(xml_tag_t), intent(in) :: tag
      type(string_t), intent(in) :: content
      integer, intent(in), optional :: unit
    end subroutine xml_tag_write_with_content
    module subroutine xml_tag_read (tag, cstream, success)
      class(xml_tag_t), intent(inout) :: tag
      type(cstream_t), intent(inout) :: cstream
      logical, intent(out) :: success
    end subroutine xml_tag_read
    module subroutine xml_tag_read_attribute (tag, string, done)
      class(xml_tag_t), intent(inout) :: tag
      type(string_t), intent(inout) :: string
      logical, intent(out) :: done
    end subroutine xml_tag_read_attribute
    module subroutine xml_tag_read_content (tag, cstream, content, closing)
      class(xml_tag_t), intent(in) :: tag
      type(cstream_t), intent(inout) :: cstream
      type(string_t), intent(out) :: content
      logical, intent(out) :: closing
    end subroutine xml_tag_read_content
  end interface

end module xml
