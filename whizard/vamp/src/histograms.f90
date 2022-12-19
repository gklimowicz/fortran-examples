! histograms.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module histograms
  use kinds
  use utils, only: find_free_unit
  implicit none
  private
  public :: create_histogram
  public :: fill_histogram
  public :: delete_histogram
  public :: write_histogram
  private :: create_histogram1, create_histogram2
  private :: fill_histogram1, fill_histogram2s, fill_histogram2v
  private :: delete_histogram1, delete_histogram2
  private :: write_histogram1, write_histogram2
  !WK! public :: write_histogram1_unit
  private :: midpoint
  private :: midpoint1, midpoint2
  interface create_histogram
     module procedure create_histogram1, create_histogram2
  end interface
  interface fill_histogram
     module procedure fill_histogram1, fill_histogram2s, fill_histogram2v
  end interface
  interface delete_histogram
     module procedure delete_histogram1, delete_histogram2
  end interface
  interface write_histogram
     module procedure write_histogram1, write_histogram2
     module procedure write_histogram1_unit
  end interface
  interface midpoint
     module procedure midpoint1, midpoint2
  end interface
  integer, parameter, private :: N_BINS_DEFAULT = 10
  type, public :: histogram
     private
     integer :: n_bins
     real(kind=default) :: x_min, x_max
     real(kind=default), dimension(:), pointer :: bins => null ()
     real(kind=default), dimension(:), pointer :: bins2 => null ()
     real(kind=default), dimension(:), pointer :: bins3 => null ()
  end type histogram
  type, public :: histogram2
     private
     integer, dimension(2) :: n_bins
     real(kind=default), dimension(2) :: x_min, x_max
     real(kind=default), dimension(:,:), pointer :: bins => null ()
     real(kind=default), dimension(:,:), pointer :: bins2 => null ()
  end type histogram2
contains
  elemental subroutine create_histogram1 (h, x_min, x_max, nb)
    type(histogram), intent(out) :: h
    real(kind=default), intent(in) :: x_min, x_max
    integer, intent(in), optional :: nb
    if (present (nb)) then
       h%n_bins = nb
    else
       h%n_bins = N_BINS_DEFAULT
    end if
    h%x_min = x_min
    h%x_max = x_max
    allocate (h%bins(0:h%n_bins+1), h%bins2(0:h%n_bins+1))
    h%bins = 0
    h%bins2 = 0
    allocate (h%bins3(0:h%n_bins+1))
    h%bins3 = 0
  end subroutine create_histogram1
  pure subroutine create_histogram2 (h, x_min, x_max, nb)
    type(histogram2), intent(out) :: h
    real(kind=default), dimension(:), intent(in) :: x_min, x_max
    integer, intent(in), dimension(:), optional :: nb
    if (present (nb)) then
       h%n_bins = nb
    else
       h%n_bins = N_BINS_DEFAULT
    end if
    h%x_min = x_min
    h%x_max = x_max
    allocate (h%bins(0:h%n_bins(1)+1,0:h%n_bins(1)+1), &
              h%bins2(0:h%n_bins(2)+1,0:h%n_bins(2)+1))
    h%bins = 0
    h%bins2 = 0
  end subroutine create_histogram2
  elemental subroutine fill_histogram1 (h, x, weight, excess)
    type(histogram), intent(inout) :: h
    real(kind=default), intent(in) :: x
    real(kind=default), intent(in), optional :: weight
    real(kind=default), intent(in), optional :: excess
    integer :: i
    if (x < h%x_min) then
       i = 0
    else if (x > h%x_max) then
       i = h%n_bins + 1
    else
       i = 1 + h%n_bins * (x - h%x_min) / (h%x_max - h%x_min)
  !WK! i = min (max (i, 0), h%n_bins + 1)
    end if
    if (present (weight)) then
       h%bins(i) = h%bins(i) + weight
       h%bins2(i) = h%bins2(i) + weight*weight
    else
       h%bins(i) = h%bins(i) + 1
       h%bins2(i) = h%bins2(i) + 1
    end if
    if (present (excess)) h%bins3(i) = h%bins3(i) + excess
  end subroutine fill_histogram1
  elemental subroutine fill_histogram2s (h, x1, x2, weight)
    type(histogram2), intent(inout) :: h
    real(kind=default), intent(in) :: x1, x2
    real(kind=default), intent(in), optional :: weight
    call fill_histogram2v (h, (/ x1, x2 /), weight)
  end subroutine fill_histogram2s
  pure subroutine fill_histogram2v (h, x, weight)
    type(histogram2), intent(inout) :: h
    real(kind=default), dimension(:), intent(in) :: x
    real(kind=default), intent(in), optional :: weight
    integer, dimension(2) :: i
    i = 1 + h%n_bins * (x - h%x_min) / (h%x_max - h%x_min)
    i = min (max (i, 0), h%n_bins + 1)
    if (present (weight)) then
       h%bins(i(1),i(2)) = h%bins(i(1),i(2)) + weight
       h%bins2(i(1),i(2)) = h%bins2(i(1),i(2)) + weight*weight
    else
       h%bins(i(1),i(2)) = h%bins(i(1),i(2)) + 1
       h%bins2(i(1),i(2)) = h%bins2(i(1),i(2)) + 1
    end if
  end subroutine fill_histogram2v
  elemental subroutine delete_histogram1 (h)
    type(histogram), intent(inout) :: h
    deallocate (h%bins, h%bins2)
    deallocate (h%bins3)
  end subroutine delete_histogram1
  elemental subroutine delete_histogram2 (h)
    type(histogram2), intent(inout) :: h
    deallocate (h%bins, h%bins2)
  end subroutine delete_histogram2
  subroutine write_histogram1 (h, name, over)
    type(histogram), intent(in) :: h
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: over
    integer :: i, iounit
    if (present (name)) then
       call find_free_unit (iounit)
       if (iounit > 0) then
          open (unit = iounit, action = "write", status = "replace", &
                file = name)
          if (present (over)) then
             if (over) then
                write (unit = iounit, fmt = *) &
                  "underflow", h%bins(0), sqrt (h%bins2(0))
             end if
          end if
          do i = 1, h%n_bins
             write (unit = iounit, fmt = *) &
                  midpoint (h, i), h%bins(i), sqrt (h%bins2(i))
          end do
          if (present (over)) then
             if (over) then
                write (unit = iounit, fmt = *) &
                  "overflow", h%bins(h%n_bins+1), &
                  sqrt (h%bins2(h%n_bins+1))
             end if
          end if
          close (unit = iounit)
       else
          print *, "write_histogram: Can't find a free unit!"
       end if
    else
       if (present (over)) then
          if (over) then
             print *, "underflow", h%bins(0), sqrt (h%bins2(0))
          end if
       end if
       do i = 1, h%n_bins
          print *, midpoint (h, i), h%bins(i), sqrt (h%bins2(i))
       end do
       if (present (over)) then
          if (over) then
             print *, "overflow", h%bins(h%n_bins+1), &
                      sqrt (h%bins2(h%n_bins+1))
          end if
       end if
    end if
  end subroutine write_histogram1
  subroutine write_histogram1_unit (h, iounit, over, show_excess)
    type(histogram), intent(in) :: h
    integer, intent(in) :: iounit
    logical, intent(in), optional :: over, show_excess
    integer :: i
    logical :: show_exc
    show_exc = .false.; if (present(show_excess)) show_exc = show_excess
    if (present (over)) then
       if (over) then
          if (show_exc) then
             write (unit = iounit, fmt = 1) &
                  "underflow", h%bins(0), sqrt (h%bins2(0)), h%bins3(0)
          else
             write (unit = iounit, fmt = 1) &
                  "underflow", h%bins(0), sqrt (h%bins2(0))
          end if
       end if
    end if
    do i = 1, h%n_bins
       if (show_exc) then
          write (unit = iounit, fmt = 1) &
               midpoint (h, i), h%bins(i), sqrt (h%bins2(i)), h%bins3(i)
       else
          write (unit = iounit, fmt = 1) &
               midpoint (h, i), h%bins(i), sqrt (h%bins2(i))
       end if
    end do
    if (present (over)) then
       if (over) then
          if (show_exc) then
             write (unit = iounit, fmt = 1) &
                  "overflow", h%bins(h%n_bins+1), &
                  sqrt (h%bins2(h%n_bins+1)), &
                  h%bins3(h%n_bins+1)
          else
             write (unit = iounit, fmt = 1) &
                  "overflow", h%bins(h%n_bins+1), &
                  sqrt (h%bins2(h%n_bins+1))
          end if
       end if
    end if
  1 format (1x,4(G16.9,2x))
  end subroutine write_histogram1_unit
  elemental function midpoint1 (h, bin) result (x)
    type(histogram), intent(in) :: h
    integer, intent(in) :: bin
    real(kind=default) :: x
    x = h%x_min + (h%x_max - h%x_min) * (bin - 0.5) / h%n_bins
  end function midpoint1
  elemental function midpoint2 (h, bin, d) result (x)
    type(histogram2), intent(in) :: h
    integer, intent(in) :: bin, d
    real(kind=default) :: x
    x = h%x_min(d) + (h%x_max(d) - h%x_min(d)) * (bin - 0.5) / h%n_bins(d)
  end function midpoint2
  subroutine write_histogram2 (h, name, over)
    type(histogram2), intent(in) :: h
    character(len=*), intent(in), optional :: name
    logical, intent(in), optional :: over
    integer :: i1, i2, iounit
    if (present (name)) then
       call find_free_unit (iounit)
       if (iounit > 0) then
          open (unit = iounit, action = "write", status = "replace", &
                file = name)
          if (present (over)) then
             if (over) then
                write (unit = iounit, fmt = *) &
                     "double underflow", h%bins(0,0), sqrt (h%bins2(0,0))
                do i2 = 1, h%n_bins(2)
                   write (unit = iounit, fmt = *) &
                        "x1 underflow", midpoint (h, i2, 2), &
                        h%bins(0,i2), sqrt (h%bins2(0,i2))
                end do
                do i1 = 1, h%n_bins(1)
                   write (unit = iounit, fmt = *) &
                        "x2 underflow", midpoint (h, i1, 1), &
                        h%bins(i1,0), sqrt (h%bins2(i1,0))
                end do
             end if
          end if
          do i1 = 1, h%n_bins(1)
             do i2 = 1, h%n_bins(2)
                write (unit = iounit, fmt = *) &
                     midpoint (h, i1, 1), midpoint (h, i2, 2), &
                     h%bins(i1,i2), sqrt (h%bins2(i1,i2))
             end do
          end do
          if (present (over)) then
             if (over) then
                do i2 = 1, h%n_bins(2)
                   write (unit = iounit, fmt = *) &
                        "x1 overflow", midpoint (h, i2, 2), &
                        h%bins(h%n_bins(1)+1,i2), &
                        sqrt (h%bins2(h%n_bins(1)+1,i2))
                end do
                do i1 = 1, h%n_bins(1)
                   write (unit = iounit, fmt = *) &
                        "x2 overflow", midpoint (h, i1, 1), &
                        h%bins(i1,h%n_bins(2)+1), &
                        sqrt (h%bins2(i1,h%n_bins(2)+1))
                end do
                write (unit = iounit, fmt = *) "double overflow", &
                     h%bins(h%n_bins(1)+1,h%n_bins(2)+1), &
                     sqrt (h%bins2(h%n_bins(1)+1,h%n_bins(2)+1))
             end if
          end if
          close (unit = iounit)
       else
          print *, "write_histogram: Can't find a free unit!"
       end if
    else
       if (present (over)) then
          if (over) then
             print *, "double underflow", h%bins(0,0), sqrt (h%bins2(0,0))
             do i2 = 1, h%n_bins(2)
                print *, "x1 underflow", midpoint (h, i2, 2), &
                     h%bins(0,i2), sqrt (h%bins2(0,i2))
             end do
             do i1 = 1, h%n_bins(1)
                print *, "x2 underflow", midpoint (h, i1, 1), &
                     h%bins(i1,0), sqrt (h%bins2(i1,0))
             end do
          end if
       end if
       do i1 = 1, h%n_bins(1)
          do i2 = 1, h%n_bins(2)
             print *, midpoint (h, i1, 1), midpoint (h, i2, 2), &
                  h%bins(i1,i2), sqrt (h%bins2(i1,i2))
          end do
       end do
       if (present (over)) then
          if (over) then
             do i2 = 1, h%n_bins(2)
                print *, "x1 overflow", midpoint (h, i2, 2), &
                     h%bins(h%n_bins(1)+1,i2), &
                     sqrt (h%bins2(h%n_bins(1)+1,i2))
             end do
             do i1 = 1, h%n_bins(1)
                print *, "x2 overflow", midpoint (h, i1, 1), &
                     h%bins(i1,h%n_bins(2)+1), &
                     sqrt (h%bins2(i1,h%n_bins(2)+1))
             end do
             print *, "double overflow", &
                  h%bins(h%n_bins(1)+1,h%n_bins(2)+1), &
                  sqrt (h%bins2(h%n_bins(1)+1,h%n_bins(2)+1))
          end if
       end if
    end if
  end subroutine write_histogram2
end module histograms
