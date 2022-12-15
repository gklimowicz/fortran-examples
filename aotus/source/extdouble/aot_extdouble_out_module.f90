! Copyright (c) 2012-2013, 2016-2017 Harald Klimach <harald@klimachs.de>
!
! Parts of this file were written by Harald Klimach for
! German Research School of Simulation Sciences and University of Siegen.
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
! OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.
! **************************************************************************** !

!> A module to produce Lua scripts with nested tables.
!!
!! This module eases the output of readable Lua scripts.
!! It takes care of indentation with nested tables, and provides a concise
!! interface to output Fortran data into Lua tables.
!! Therefore this module is somehow the counter-part to the reading functions,
!! however, it is almost completely independent and relies purely on Fortran
!! output methods. Thus this module could stand alone, along with the
!! flu_kinds_module without the Lua library.
module aot_extdouble_out_module
  use aot_out_general_module, only: aot_out_type, aot_out_open, aot_out_close, &
    &                               aot_out_open_table, aot_out_close_table, &
    &                               aot_out_breakline
  use aot_extdouble_top_module, only: xdble_k

  implicit none

  public :: aot_out_val

  !> Put Fortran intrinsic types into the script.
  !!
  !! Scalar values and one-dimensional arrays are supported.
  !! Here we add support for extdouble precision.
  !! NOTE however, that the used format will only be in double precision, as
  !! Lua does not provide higher accuracy right now anyway.
  interface aot_out_val
    ! scalars
    module procedure aot_out_val_extdouble

    ! arrays
    module procedure aot_out_val_arr_extdouble
  end interface

  private

contains

  !>  Put extdouble variables into the Lua script.
  !!
  !! The value is passed in with val, optionally you can assign a name to it
  !! with the vname argument. If it should be put on the same line as the
  !! previous entry, you have to set advance_previous=.false.
  subroutine aot_out_val_extdouble(put_conf, val, vname, advance_previous)
    !------------------------------------------------------------------------
    type(aot_out_type), intent(inout)  :: put_conf
    character(len=*), optional, intent(in) :: vname
    logical, optional, intent(in) :: advance_previous
    real(kind=xdble_k), intent(in) :: val
    !------------------------------------------------------------------------
    character(len=3) :: adv_string
    !------------------------------------------------------------------------

    if (put_conf%level > 0) then
      ! Leave the advancing to the next entry in the table.
      adv_string = 'no'
    else
      ! Not within a table, finalize the global definition with a newline.
      adv_string = 'yes'
    end if

    call aot_out_breakline(put_conf, advance_previous)

    if (present(vname)) then
      write(put_conf%outunit, fmt="(a,EN24.15)", advance=adv_string) &
        & trim(vname)//" = ", val
    else
      write(put_conf%outunit, fmt="(EN24.15)", advance=adv_string) val
    end if

  end subroutine aot_out_val_extdouble
! *****************************************************************************!


! *****************************************************************************!
  !> This is a vectorized version of the value output.
  !!
  !! It takes a one-dimensional array and puts it into a table. The parameters
  !! have the usual meanings, as in the scalar routines, however and additional
  !! argument (max_per_line) allows the specification of the number of elements
  !! that might be put onto a single line.
  !! The first entry will be placed into the same line as the opening brace, and
  !! the closing brace will be put on the same line, as the last entry.
  subroutine aot_out_val_arr_extdouble(put_conf, val, vname, advance_previous, &
    &                                   max_per_line)
    !------------------------------------------------------------------------
    !> Lua script to write the array into.
    type(aot_out_type), intent(inout)  :: put_conf

    !> Name for this array
    character(len=*), optional, intent(in) :: vname

    !> Actual data to write into the script
    real(kind=xdble_k), intent(in) :: val(:)

    !> Flag if this array should be put on the same line as the last entry of
    !! the parent table.
    logical, optional, intent(in) :: advance_previous

    !> Maximal number of entries to put into a single line.
    !! Defaults to 3.
    integer, optional, intent(in) :: max_per_line
    !------------------------------------------------------------------------
    integer :: i
    integer :: nVals
    integer :: mpl
    logical :: bline
    !------------------------------------------------------------------------

    if (present(max_per_line)) then
      mpl = max_per_line
    else
      mpl = 3
    end if

    ! Opening the table(subtable for array actually)
    call aot_out_open_table(put_conf, vname, &
      &                     advance_previous = advance_previous)

    nVals = size(val)
    if (nVals > 0) then
      ! Always put the first entry on the same line as the opening brace.
      call aot_out_val(put_conf, val(1), advance_previous = .false.)

      do i=2,nVals
        ! Output each entry and break the line after mpl entries on a line.
        bline = (mod(i-1, mpl) == 0)
        call aot_out_val(put_conf, val(i), advance_previous = bline)
      end do
    end if

    ! Always put the closing brace on the same line as the last entry.
    call aot_out_close_table(put_conf, advance_previous = .false.)

  end subroutine aot_out_val_arr_extdouble
! *****************************************************************************!

end module aot_extdouble_out_module
