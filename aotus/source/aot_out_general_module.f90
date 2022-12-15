! Copyright (c) 2012-2016, 2019 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2016 Kannan Masilamani <kannan.masilamani@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach for
! German Research School of Simulation Sciences and University of Siegen.
! Parts of this file were written by Kannan Masilamani for University of Siegen.
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

!> Collection of general operations required for the output of Lua scripts.
module aot_out_general_module

  implicit none

  private

  public :: aot_out_type
  public :: aot_out_open
  public :: aot_out_close
  public :: aot_out_open_table
  public :: aot_out_close_table
  public :: aot_out_breakline
  public :: aot_out_toChunk

  !> This type provides the internal representation of the opened Lua script.
  !!
  !! It is used to keep track of the state in the script internally.
  type aot_out_type
    integer :: outunit !! Unit to write to
    integer :: indent !! Indentation level (number of spaces)
    integer :: stack(100) !! Number of entries on each level
    integer :: level !! Current nesting level in tables
    logical :: externalOpen !! Flag if file opened outside the aot_out scope
    integer :: in_step !! Number of spaces for each indentation level
  end type


contains


! **************************************************************************** !
  !> Open the file to write to and return a handle (put_conf) to it.
  !!
  !! This will overwrite the given file, if it already exists.
  !! Either filename of outUnit has to be specified, use outUnit to write to a
  !! pre-connected file.
  !! If both are given, the file will be opened and connected to a new unit,
  !! outUnit is ignored in this case.
  subroutine aot_out_open(put_conf, filename, outUnit, indentation, outstat)
    !------------------------------------------------------------------------
    type(aot_out_type), intent(out) :: put_conf !! Handle for the file
    character(len=*), optional, intent(in) :: filename !! File to open
    integer, optional, intent(in) :: outUnit !! Pre-connected unit to write to
    integer, optional, intent(in) :: indentation !! Spacer per indentation level

    !> IO status of the open operation for the given filename or an indication
    !! whether the given outUnit is actually connected to an open file.
    !!
    !! This returns 0 if the the returned unit has properly been properly
    !! connected to the file.
    integer, optional, intent(out) :: outstat
    !------------------------------------------------------------------------
    integer :: iError
    logical :: isOpen
    !------------------------------------------------------------------------

    if (present(indentation)) then
      put_conf%in_step = indentation
    else
      put_conf%in_step = 4
    end if

    if (present(filename)) then
      put_conf%outunit = newunit()
      open(unit = put_conf%outunit, file = trim(filename), action = 'write', &
        &  status='replace', recl = 360, iostat=iError)
      put_conf%externalOpen = .false.
    else if (present(outUnit)) then
      inquire(unit=outUnit, opened=isOpen)
      if (isOpen) then
        iError = 0
      else
        iError = -10
      end if
      put_conf%externalOpen = .true.
      put_conf%outunit = outUnit
    end if

    if (present(outstat)) outstat = iError

    put_conf%indent = 0
    put_conf%stack(:) = 0
    put_conf%level = 0

  end subroutine aot_out_open
! **************************************************************************** !


! **************************************************************************** !
  !>  Close the opened script again.
  !!
  !! This will close the file, if the data was not written to a pre-connected
  !! unit (that is the file for the script was opened in the aot_out_open).
  subroutine aot_out_close(put_conf)
    !------------------------------------------------------------------------
    type(aot_out_type), intent(inout)  :: put_conf
    !------------------------------------------------------------------------
    if ( .not. put_conf%externalOpen ) close( put_conf%outunit )
  end subroutine aot_out_close
! **************************************************************************** !


! **************************************************************************** !
  !> Start a new table to write to.
  !!
  !! You can give the table a name with the tname argument.
  !! If the table definition should NOT start on a new line, you have to pass
  !! in an advance_previous = .false.
  subroutine aot_out_open_table(put_conf, tname, advance_previous)
    !------------------------------------------------------------------------
    type(aot_out_type), intent(inout)  :: put_conf
    character(len=*), optional, intent(in) :: tname
    logical, optional, intent(in) :: advance_previous
    !------------------------------------------------------------------------

    call aot_out_breakline(put_conf, advance_previous)

    if (present(tname)) then
      write(put_conf%outunit, fmt='(a)', advance='no') trim(tname)//' = {'
    else
      write(put_conf%outunit, fmt='(a)', advance='no') '{'
    end if

    put_conf%level = put_conf%level + 1
    put_conf%indent = put_conf%indent + put_conf%in_step

  end subroutine aot_out_open_table
! **************************************************************************** !


! **************************************************************************** !
  !>  Close the current table.
  !!
  !! The table on the current table is closed with a curly bracket.
  !! If this bracket should be put to the same line as the last entry of the
  !! table, you have to set advance_previous = .false.
  subroutine aot_out_close_table(put_conf, advance_previous)
    !------------------------------------------------------------------------
    type(aot_out_type), intent(inout)  :: put_conf
    logical, optional, intent(in) :: advance_previous
    !------------------------------------------------------------------------
    logical :: loc_adv_prev
    character(len=max(put_conf%indent-put_conf%in_step,0)) :: indent
    character(len=3) :: adv_string
    !------------------------------------------------------------------------

    indent = ''
    adv_string = 'yes'

    if (present(advance_previous)) then
      loc_adv_prev = advance_previous
    else
      loc_adv_prev = .true.
    end if

    put_conf%indent = max(put_conf%indent - put_conf%in_step, 0)
    put_conf%stack(put_conf%level) = 0
    put_conf%level = max(put_conf%level - 1, 0)

    if (put_conf%level > 0) then
      ! Do not advance, to let the next entry append the separator to the line.
      adv_string = 'no'
    end if

    ! Close last entry without separator.
    if (loc_adv_prev) then
      ! Closing brace should be on new line.
      write(put_conf%outunit,*) ''
      write(put_conf%outunit, fmt="(a)", advance=adv_string) indent//'}'
    else
      ! Closing brace on same line as last entry.
      write(put_conf%outunit, fmt="(a)", advance=adv_string) ' }'
    end if

  end subroutine aot_out_close_table
! **************************************************************************** !




! **************************************************************************** !
  !> This subroutine takes care of the proper linebreaking in Lua-Tables.
  !!
  !! It takes care of a proper line-continuation, depending on the optional
  !! advance_previous flag and increases the count of elements in the current
  !! table.
  !! The default is to put each entry on a new line, if it should be on the
  !! same line advance_previous = .false. has to be set.
  subroutine aot_out_breakline(put_conf, advance_previous)
    type(aot_out_type), intent(inout)  :: put_conf
    logical, optional, intent(in) :: advance_previous

    character(len=put_conf%indent) :: indent
    character :: sep
    logical :: loc_adv_prev

    indent = ''
    if (present(advance_previous)) then
      loc_adv_prev = advance_previous
    else
      loc_adv_prev = .true.
    end if

    lev_if: if (put_conf%level > 0) then

      if (put_conf%stack(put_conf%level) > 0) then
        ! Use the separator to close the previous entry.
        sep = ','
      else
        ! First entry, nothing to separate yet.
        sep = ''
      end if

      if (loc_adv_prev) then
        write(put_conf%outunit, fmt='(a)') trim(sep)
        write(put_conf%outunit, fmt='(a)', advance='no') indent
      else
        write(put_conf%outunit, fmt='(a)', advance='no') trim(sep)//" "
      end if

      put_conf%stack(put_conf%level) = put_conf%stack(put_conf%level) + 1

    else if (put_conf%level .eq. 0)then
      write(put_conf%outunit, fmt='(a)', advance='no') " "
    end if lev_if

  end subroutine aot_out_breakline
! **************************************************************************** !


! **************************************************************************** !
  !> This subroutine converts information written in outunit to string
  subroutine aot_out_toChunk(out_conf, chunk, ErrCode, ErrString)
    type(aot_out_type), intent(in)  :: out_conf

    !> String with Lua code to load.
    character(len=*), intent(out) :: chunk

    !> Error code returned by Lua during loading or executing the file.
    !!
    !! This optional parameter might be used to react on errors in the calling
    !! side. If neither ErrCode nor ErrString are given, this subroutine will
    !! stop the program execution and print the error message
    integer, intent(out), optional :: ErrCode

    !> Error description
    !!
    !! This optional argument holds the error message in case something
    !! went wrong. It can be used to provide some feedback to the user in the
    !! calling routine. If neither ErrCode nor ErrString are provided,
    !! this subroutine will print the error message and stop program execution.
    character(len=*), intent(out), optional :: ErrString

    logical :: stop_on_error
    integer :: error
    integer :: chunk_len, chunk_left, read_len
    character(len=320) :: err_string
    logical :: unitOpened
    integer :: read_stat
    character(len=320) :: chunk_line

    stop_on_error = .not.(present(ErrString) .or. present(ErrCode))
    error = 0
    err_string = ''

    ! length of chunk
    chunk_len = len(chunk)

    inquire(unit=out_conf%outunit, opened=unitOpened)
    if (unitOpened) then
      chunk = ''
      chunk_left = chunk_len
      rewind(out_conf%outunit)
      do
        read(out_conf%outunit,'(a)', iostat=read_stat) chunk_line
        read_len = len(trim(chunk_line))
        if (read_stat /= 0) then
          if (read_stat > 0) then
            error = read_stat
            err_string = 'Error reading out conf unit'
          end if
          exit ! exit reading
        end if
        if (chunk_left >= read_len) then
          chunk_left = chunk_left - len(trim(chunk))
          chunk =  trim(chunk)//new_line('x')//trim(chunk_line)
        else
          error = 2
          err_string = 'Reached limit of output string length'
          exit
        end if
      end do
    else
      error = 1
      err_string = 'Output conf unit is not opened'
    end if

    if (present(ErrCode)) then
      ErrCode = error
    end if

    if (present(ErrString)) then
      ErrString = err_string
    end if

    if (error /= 0) then
      if (stop_on_error) then
        write(*,*) 'From aot_out_toChunk: '//trim(err_string)
        STOP
      end if
    end if
  end subroutine aot_out_toChunk
! **************************************************************************** !


! **************************************************************************** !
  !> Helper function to provide new unit, as long as F2008 newunit argument
  !! in open statement is not commonly available.
  !!
  !! To be used right in front of the open statement like this:
  !!  myUnit = newunit()
  !!  open(myUnit, ...)
  function newunit() result(nu)
    integer :: nu

    logical :: connected

    nu = 21
    inquire(unit=nu, opened=connected)
    do while(connected)
      nu = nu + 1
      inquire(unit=nu, opened=connected)
    end do
  end function newunit
! **************************************************************************** !


end module aot_out_general_module
