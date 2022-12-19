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

module sf_mappings_uti

  use kinds, only: default
  use format_defs, only: FMT_11, FMT_12, FMT_13, FMT_14, FMT_15, FMT_16

  use sf_mappings

  implicit none
  private

  public :: sf_mappings_1
  public :: sf_mappings_2
  public :: sf_mappings_3
  public :: sf_mappings_4
  public :: sf_mappings_5
  public :: sf_mappings_6
  public :: sf_mappings_7
  public :: sf_mappings_8
  public :: sf_mappings_9
  public :: sf_mappings_10
  public :: sf_mappings_11
  public :: sf_mappings_12
  public :: sf_mappings_13
  public :: sf_mappings_14
  public :: sf_mappings_15
  public :: sf_mappings_16

contains

  subroutine sf_mappings_1 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p

    write (u, "(A)")  "* Test output: sf_mappings_1"
    write (u, "(A)")  "*   Purpose: probe standard mapping"
    write (u, "(A)")

    allocate (sf_s_mapping_t :: mapping)
    select type (mapping)
    type is (sf_s_mapping_t)
       call mapping%init ()
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.1):"
    p = [0.1_default, 0.1_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)
    allocate (sf_s_mapping_t :: mapping)
    select type (mapping)
    type is (sf_s_mapping_t)
       call mapping%init (power=2._default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    write (u, *)
    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.1):"
    p = [0.1_default, 0.1_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)")  "I =", mapping%integral (100000)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_1"

  end subroutine sf_mappings_1

  subroutine sf_mappings_2 (u)
    integer, intent(in) :: u
    type(sf_channel_t), dimension(:), allocatable :: channel
    integer :: c

    write (u, "(A)")  "* Test output: sf_mappings_2"
    write (u, "(A)")  "*   Purpose: construct and display &
         &mapping-channel objects"
    write (u, "(A)")

    call allocate_sf_channels (channel, n_channel = 8, n_strfun = 2)
    call channel(2)%activate_mapping ([1])
    call channel(3)%set_s_mapping ([1,2])
    call channel(4)%set_s_mapping ([1,2], power=2._default)
    call channel(5)%set_res_mapping ([1,2], m = 0.5_default, w = 0.1_default, single = .false.)
    call channel(6)%set_os_mapping ([1,2], m = 0.5_default, single = .false.)
    call channel(7)%set_res_mapping ([1], m = 0.5_default, w = 0.1_default, single = .true.)
    call channel(8)%set_os_mapping ([1], m = 0.5_default, single = .true.)

    call channel(3)%set_par_index (1, 1)
    call channel(3)%set_par_index (2, 4)

    call channel(4)%set_par_index (1, 1)
    call channel(4)%set_par_index (2, 4)

    call channel(5)%set_par_index (1, 1)
    call channel(5)%set_par_index (2, 3)

    call channel(6)%set_par_index (1, 1)
    call channel(6)%set_par_index (2, 2)

    call channel(7)%set_par_index (1, 1)

    call channel(8)%set_par_index (1, 1)

    do c = 1, size (channel)
       write (u, "(I0,':')", advance="no")  c
       call channel(c)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_2"

  end subroutine sf_mappings_2

  subroutine sf_mappings_3 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p

    write (u, "(A)")  "* Test output: sf_mappings_3"
    write (u, "(A)")  "*   Purpose: probe resonance pair mapping"
    write (u, "(A)")

    allocate (sf_res_mapping_t :: mapping)
    select type (mapping)
    type is (sf_res_mapping_t)
       call mapping%init (0.5_default, 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.1):"
    p = [0.1_default, 0.1_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_3"

  end subroutine sf_mappings_3

  subroutine sf_mappings_4 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p

    write (u, "(A)")  "* Test output: sf_mappings_4"
    write (u, "(A)")  "*   Purpose: probe on-shell pair mapping"
    write (u, "(A)")

    allocate (sf_os_mapping_t :: mapping)
    select type (mapping)
    type is (sf_os_mapping_t)
       call mapping%init (0.5_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.1):"
    p = [0._default, 0.1_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0,1.0):"
    p = [0._default, 1.0_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_4"

  end subroutine sf_mappings_4

  subroutine sf_mappings_5 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p

    write (u, "(A)")  "* Test output: sf_mappings_5"
    write (u, "(A)")  "*   Purpose: probe endpoint pair mapping"
    write (u, "(A)")

    allocate (sf_ep_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ep_mapping_t)
       call mapping%init ()
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_5"

  end subroutine sf_mappings_5

  subroutine sf_mappings_6 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p

    write (u, "(A)")  "* Test output: sf_mappings_6"
    write (u, "(A)")  "*   Purpose: probe endpoint resonant mapping"
    write (u, "(A)")

    allocate (sf_epr_mapping_t :: mapping)
    select type (mapping)
    type is (sf_epr_mapping_t)
       call mapping%init (a = 1._default, m = 0.5_default, w = 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Same mapping without resonance:"
    write (u, "(A)")

    allocate (sf_epr_mapping_t :: mapping)
    select type (mapping)
    type is (sf_epr_mapping_t)
       call mapping%init (a = 1._default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_6"

  end subroutine sf_mappings_6

  subroutine sf_mappings_7 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p

    write (u, "(A)")  "* Test output: sf_mappings_7"
    write (u, "(A)")  "*   Purpose: probe endpoint on-shell mapping"
    write (u, "(A)")

    allocate (sf_epo_mapping_t :: mapping)
    select type (mapping)
    type is (sf_epo_mapping_t)
       call mapping%init (a = 1._default, m = 0.5_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0):"
    p = [0._default, 0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1,0.5):"
    p = [0.1_default, 0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_7"

  end subroutine sf_mappings_7

  subroutine sf_mappings_8 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p, pb

    write (u, "(A)")  "* Test output: sf_mappings_8"
    write (u, "(A)")  "*   Purpose: probe power pair mapping"
    write (u, "(A)")

    allocate (sf_ip_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ip_mapping_t)
       call mapping%init (eps = 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.5):"
    p = [0._default, 0.5_default]
    pb= [1._default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    pb= [0.5_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9,0.5):"
    p = [0.9_default, 0.5_default]
    pb= [0.1_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    pb= [0.3_default, 0.8_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.8):"
    p = [0.7_default, 0.8_default]
    pb= [0.3_default, 0.2_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.99,0.02):"
    p = [0.99_default, 0.02_default]
    pb= [0.01_default, 0.98_default]
    call mapping%check (u, p, pb, FMT_14, FMT_12)

    write (u, *)
    write (u, "(A)")  "Probe at (0.99,0.98):"
    p = [0.99_default, 0.98_default]
    pb= [0.01_default, 0.02_default]
    call mapping%check (u, p, pb, FMT_14, FMT_12)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_8"

  end subroutine sf_mappings_8

  subroutine sf_mappings_9 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p, pb

    write (u, "(A)")  "* Test output: sf_mappings_9"
    write (u, "(A)")  "*   Purpose: probe power resonant pair mapping"
    write (u, "(A)")

    allocate (sf_ipr_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ipr_mapping_t)
       call mapping%init (eps = 0.1_default, m = 0.5_default, w = 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.5):"
    p = [0._default, 0.5_default]
    pb= [1._default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    pb= [0.5_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9,0.5):"
    p = [0.9_default, 0.5_default]
    pb= [0.1_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    pb= [0.3_default, 0.8_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.8):"
    p = [0.7_default, 0.8_default]
    pb= [0.3_default, 0.2_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9999,0.02):"
    p = [0.9999_default, 0.02_default]
    pb= [0.0001_default, 0.98_default]
    call mapping%check (u, p, pb, FMT_11, FMT_12)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9999,0.98):"
    p = [0.9999_default, 0.98_default]
    pb= [0.0001_default, 0.02_default]
    call mapping%check (u, p, pb, FMT_11, FMT_12)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Same mapping without resonance:"
    write (u, "(A)")

    allocate (sf_ipr_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ipr_mapping_t)
       call mapping%init (eps = 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.5):"
    p = [0._default, 0.5_default]
    pb= [1._default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.5,0.5):"
    p = [0.5_default, 0.5_default]
    pb= [0.5_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9,0.5):"
    p = [0.9_default, 0.5_default]
    pb= [0.1_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.2):"
    p = [0.7_default, 0.2_default]
    pb= [0.3_default, 0.8_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7,0.8):"
    p = [0.7_default, 0.8_default]
    pb= [0.3_default, 0.2_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_9"

  end subroutine sf_mappings_9

  subroutine sf_mappings_10 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(2) :: p, pb

    write (u, "(A)")  "* Test output: sf_mappings_10"
    write (u, "(A)")  "*   Purpose: probe power on-shell mapping"
    write (u, "(A)")

    allocate (sf_ipo_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ipo_mapping_t)
       call mapping%init (eps = 0.1_default, m = 0.5_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.5):"
    p = [0._default, 0.5_default]
    pb= [1._default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.02):"
    p = [0._default, 0.02_default]
    pb= [1._default, 0.98_default]
    call mapping%check (u, p, pb, FMT_15, FMT_12)

    write (u, *)
    write (u, "(A)")  "Probe at (0,0.98):"
    p = [0._default, 0.98_default]
    pb= [1._default, 0.02_default]
    call mapping%check (u, p, pb, FMT_15, FMT_12)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_10"

  end subroutine sf_mappings_10

  subroutine sf_mappings_11 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(4) :: p, pb

    write (u, "(A)")  "* Test output: sf_mappings_11"
    write (u, "(A)")  "*   Purpose: probe power pair mapping"
    write (u, "(A)")

    allocate (sf_ei_mapping_t :: mapping)
    select type (mapping)
    type is (sf_ei_mapping_t)
       call mapping%init (eps = 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
       call mapping%set_index (3, 3)
       call mapping%set_index (4, 4)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0.5, 0.5, 0.5, 0.5):"
    p = [0.5_default, 0.5_default, 0.5_default, 0.5_default]
    pb= [0.5_default, 0.5_default, 0.5_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7, 0.2, 0.4, 0.8):"
    p = [0.7_default, 0.2_default, 0.4_default, 0.8_default]
    pb= [0.3_default, 0.8_default, 0.6_default, 0.2_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9, 0.06, 0.95, 0.1):"
    p = [0.9_default, 0.06_default, 0.95_default, 0.1_default]
    pb= [0.1_default, 0.94_default, 0.05_default, 0.9_default]
    call mapping%check (u, p, pb, FMT_13, FMT_12)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_11"

  end subroutine sf_mappings_11

  subroutine sf_mappings_12 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(4) :: p, pb

    write (u, "(A)")  "* Test output: sf_mappings_12"
    write (u, "(A)")  "*   Purpose: probe resonant combined mapping"
    write (u, "(A)")

    allocate (sf_eir_mapping_t :: mapping)
    select type (mapping)
    type is (sf_eir_mapping_t)
       call mapping%init (a = 1._default, &
            eps = 0.1_default, m = 0.5_default, w = 0.1_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
       call mapping%set_index (3, 3)
       call mapping%set_index (4, 4)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0.5, 0.5, 0.5, 0.5):"
    p = [0.5_default, 0.5_default, 0.5_default, 0.5_default]
    pb= [0.5_default, 0.5_default, 0.5_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7, 0.2, 0.4, 0.8):"
    p = [0.7_default, 0.2_default, 0.4_default, 0.8_default]
    pb= [0.3_default, 0.8_default, 0.6_default, 0.2_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9, 0.06, 0.95, 0.1):"
    p = [0.9_default, 0.06_default, 0.95_default, 0.1_default]
    pb= [0.1_default, 0.94_default, 0.05_default, 0.9_default]
    call mapping%check (u, p, pb, FMT_15, FMT_12)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_12"

  end subroutine sf_mappings_12

  subroutine sf_mappings_13 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(4) :: p, pb

    write (u, "(A)")  "* Test output: sf_mappings_13"
    write (u, "(A)")  "*   Purpose: probe on-shell combined mapping"
    write (u, "(A)")

    allocate (sf_eio_mapping_t :: mapping)
    select type (mapping)
    type is (sf_eio_mapping_t)
       call mapping%init (a = 1._default, eps = 0.1_default, m = 0.5_default)
       call mapping%set_index (1, 1)
       call mapping%set_index (2, 2)
       call mapping%set_index (3, 3)
       call mapping%set_index (4, 4)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0.5, 0.5, 0.5, 0.5):"
    p = [0.5_default, 0.5_default, 0.5_default, 0.5_default]
    pb= [0.5_default, 0.5_default, 0.5_default, 0.5_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.7, 0.2, 0.4, 0.8):"
    p = [0.7_default, 0.2_default, 0.4_default, 0.8_default]
    pb= [0.3_default, 0.8_default, 0.6_default, 0.2_default]
    call mapping%check (u, p, pb, FMT_16)

    write (u, *)
    write (u, "(A)")  "Probe at (0.9, 0.06, 0.95, 0.1):"
    p = [0.9_default, 0.06_default, 0.95_default, 0.1_default]
    pb= [0.1_default, 0.94_default, 0.05_default, 0.9_default]
    call mapping%check (u, p, pb, FMT_14, FMT_12)

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_13"

  end subroutine sf_mappings_13

  subroutine sf_mappings_14 (u)
    integer, intent(in) :: u
    real(default), dimension(2) :: p2, r2
    real(default), dimension(1) :: p1, r1
    real(default) :: f, x_free, m2

    write (u, "(A)")  "* Test output: sf_mappings_14"
    write (u, "(A)")  "*   Purpose: probe rescaling in os mapping"
    write (u, "(A)")

    x_free = 0.9_default
    m2 = 0.5_default

    write (u, "(A)")  "* Two parameters"
    write (u, "(A)")

    p2 = [0.1_default, 0.2_default]

    call map_on_shell (r2, f, p2, -log (m2), x_free)

    write (u, "(A,9(1x," // FMT_14 // "))")  "p =", p2
    write (u, "(A,9(1x," // FMT_14 // "))")  "r =", r2
    write (u, "(A,9(1x," // FMT_14 // "))")  "f =", f
    write (u, "(A,9(1x," // FMT_14 // "))")  "*r=", x_free * product (r2)

    write (u, *)

    call map_on_shell_inverse (r2, f, p2, -log (m2), x_free)

    write (u, "(A,9(1x," // FMT_14 // "))")  "p =", p2
    write (u, "(A,9(1x," // FMT_14 // "))")  "r =", r2
    write (u, "(A,9(1x," // FMT_14 // "))")  "f =", f
    write (u, "(A,9(1x," // FMT_14 // "))")  "*r=", x_free * product (r2)

    write (u, "(A)")
    write (u, "(A)")  "* One parameter"
    write (u, "(A)")

    p1 = [0.1_default]

    call map_on_shell_single (r1, f, p1, -log (m2), x_free)

    write (u, "(A,9(1x," // FMT_14 // "))")  "p =", p1
    write (u, "(A,9(1x," // FMT_14 // "))")  "r =", r1
    write (u, "(A,9(1x," // FMT_14 // "))")  "f =", f
    write (u, "(A,9(1x," // FMT_14 // "))")  "*r=", x_free * product (r1)

    write (u, *)

    call map_on_shell_single_inverse (r1, f, p1, -log (m2), x_free)

    write (u, "(A,9(1x," // FMT_14 // "))")  "p =", p1
    write (u, "(A,9(1x," // FMT_14 // "))")  "r =", r1
    write (u, "(A,9(1x," // FMT_14 // "))")  "f =", f
    write (u, "(A,9(1x," // FMT_14 // "))")  "*r=", x_free * product (r1)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_14"

  end subroutine sf_mappings_14

  subroutine sf_mappings_15 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(1) :: p

    write (u, "(A)")  "* Test output: sf_mappings_15"
    write (u, "(A)")  "*   Purpose: probe resonance single mapping"
    write (u, "(A)")

    allocate (sf_res_mapping_single_t :: mapping)
    select type (mapping)
    type is (sf_res_mapping_single_t)
       call mapping%init (0.5_default, 0.1_default)
       call mapping%set_index (1, 1)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0):"
    p = [0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5):"
    p = [0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.1):"
    p = [0.1_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_15"

  end subroutine sf_mappings_15

  subroutine sf_mappings_16 (u)
    integer, intent(in) :: u
    class(sf_mapping_t), allocatable :: mapping
    real(default), dimension(1) :: p

    write (u, "(A)")  "* Test output: sf_mappings_16"
    write (u, "(A)")  "*   Purpose: probe on-shell single mapping"
    write (u, "(A)")

    allocate (sf_os_mapping_single_t :: mapping)
    select type (mapping)
    type is (sf_os_mapping_single_t)
       call mapping%init (0.5_default)
       call mapping%set_index (1, 1)
    end select

    call mapping%write (u)

    write (u, *)
    write (u, "(A)")  "Probe at (0):"
    p = [0._default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Probe at (0.5):"
    p = [0.5_default]
    call mapping%check (u, p, 1-p, "F7.5")

    write (u, *)
    write (u, "(A)")  "Compute integral:"
    write (u, "(3x,A,1x,F7.5)") "I =", mapping%integral (100000)

    deallocate (mapping)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_mappings_16"

  end subroutine sf_mappings_16


end module sf_mappings_uti
