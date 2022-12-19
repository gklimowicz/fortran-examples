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

module phs_forests_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use io_units
  use format_defs, only: FMT_12
  use lorentz
  use flavors
  use interactions
  use model_data
  use mappings
  use phs_base
  use resonances, only: resonance_history_set_t

  use phs_forests

  implicit none
  private

  public :: phs_forest_1
  public :: phs_forest_2

contains

  subroutine phs_forest_1 (u)
    use os_interface
    integer, intent(in) :: u
    type(phs_forest_t) :: forest
    type(phs_channel_t), dimension(:), allocatable :: channel
    type(model_data_t), target :: model
    type(string_t) :: process_id
    type(flavor_t), dimension(5) :: flv
    type(string_t) :: filename
    type(interaction_t) :: int
    integer :: unit_fix
    type(mapping_defaults_t) :: mapping_defaults
    logical :: found_process, ok
    integer :: n_channel, ch, i
    logical, dimension(4) :: active = .true.
    real(default) :: sqrts = 1000
    real(default), dimension(5,4) :: x
    real(default), dimension(4) :: factor
    real(default) :: volume

    write (u, "(A)")  "* Test output: PHS forest"
    write (u, "(A)")  "*   Purpose: test PHS forest routines"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Create phase-space file 'phs_forest_test.phs'"
    write (u, "(A)")

    call flv%init ([11, -11, 11, -11, 22], model)
    unit_fix = free_unit ()
    open (file="phs_forest_test.phs", unit=unit_fix, action="write")
    write (unit_fix, *) "process foo"
    write (unit_fix, *) 'md5sum_process    = "6ABA33BC2927925D0F073B1C1170780A"'
    write (unit_fix, *) 'md5sum_model_par  = "1A0B151EE6E2DEB92D880320355A3EAB"'
    write (unit_fix, *) 'md5sum_phs_config = "B6A8877058809A8BDD54753CDAB83ACE"'
    write (unit_fix, *) "sqrts         =    100.00000000000000"
    write (unit_fix, *) "m_threshold_s =    50.000000000000000"
    write (unit_fix, *) "m_threshold_t =    100.00000000000000"
    write (unit_fix, *) "off_shell =            2"
    write (unit_fix, *) "t_channel =            6"
    write (unit_fix, *) "keep_nonresonant =  F"
    write (unit_fix, *) ""
    write (unit_fix, *) "  grove"
    write (unit_fix, *) "    tree 3 7"
    write (unit_fix, *) "      map 3 s_channel 23"
    write (unit_fix, *) "    tree 5 7"
    write (unit_fix, *) "    tree 6 7"
    write (unit_fix, *) "  grove"
    write (unit_fix, *) "    tree 9 11"
    write (unit_fix, *) "      map 9 t_channel 22"
    close (unit_fix)

    write (u, "(A)")
    write (u, "(A)")  "* Read phase-space file 'phs_forest_test.phs'"

    call syntax_phs_forest_init ()
    process_id = "foo"
    filename = "phs_forest_test.phs"
    call forest%read (filename, process_id, 2, 3, model, found_process)

    write (u, "(A)")
    write (u, "(A)")  "* Set parameters, flavors, equiv, momenta"
    write (u, "(A)")

    call forest%set_flavors (flv)
    call forest%set_parameters (mapping_defaults, .false.)
    call forest%setup_prt_combinations ()
    call forest%set_equivalences ()
    call int%basic_init (2, 0, 3)
    call int%set_momentum &
         (vector4_moving (500._default, 500._default, 3), 1)
    call int%set_momentum &
         (vector4_moving (500._default,-500._default, 3), 2)
    call forest%set_prt_in (int)
    n_channel = 2
    x = 0
    x(:,n_channel) = [0.3, 0.4, 0.1, 0.9, 0.6]
    write (u, "(A)")  "   Input values:"
    write (u, "(3x,5(1x," // FMT_12 // "))")  x(:,n_channel)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluating phase space"

    call forest%evaluate_selected_channel (n_channel, active, sqrts, &
         x, factor, volume, ok)
    call forest%evaluate_other_channels (n_channel, active, sqrts, &
         x, factor, combine=.true.)
    call forest%get_prt_out (int)
    write (u, "(A)")  "   Output values:"
    do ch = 1, 4
       write (u, "(3x,5(1x," // FMT_12 // "))")  x(:,ch)
    end do
    call int%basic_write (u)
    write (u, "(A)")  "   Factors:"
    write (u, "(3x,5(1x," // FMT_12 // "))")  factor
    write (u, "(A)")  "   Volume:"
    write (u, "(3x,5(1x," // FMT_12 // "))")  volume
    call forest%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Compute equivalences"

    n_channel = 4
    allocate (channel (n_channel))
    call forest%get_equivalences (channel, .true.)
    do i = 1, n_channel
       write (u, "(1x,I0,':')", advance = "no")  ch
       call channel(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()
    call forest%final ()
    call syntax_phs_forest_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_forest_1"

  end subroutine phs_forest_1

  subroutine phs_forest_2 (u)
    use os_interface
    integer, intent(in) :: u
    integer :: unit_fix
    type(phs_forest_t) :: forest
    type(model_data_t), target :: model
    type(string_t) :: process_id
    type(string_t) :: filename
    logical :: found_process
    type(resonance_history_set_t) :: res_set
    integer :: i

    write (u, "(A)")  "* Test output: phs_forest_2"
    write (u, "(A)")  "*   Purpose: test PHS forest routines"
    write (u, "(A)")

    write (u, "(A)")  "* Reading model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Create phase-space file 'phs_forest_2.phs'"
    write (u, "(A)")

    unit_fix = free_unit ()
    open (file="phs_forest_2.phs", unit=unit_fix, action="write")
    write (unit_fix, *) "process foo"
    write (unit_fix, *) 'md5sum_process    = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"'
    write (unit_fix, *) 'md5sum_model_par  = "1A0B151EE6E2DEB92D880320355A3EAB"'
    write (unit_fix, *) 'md5sum_phs_config = "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"'
    write (unit_fix, *) "sqrts         =    100.00000000000000"
    write (unit_fix, *) "m_threshold_s =    50.000000000000000"
    write (unit_fix, *) "m_threshold_t =    100.00000000000000"
    write (unit_fix, *) "off_shell =            2"
    write (unit_fix, *) "t_channel =            6"
    write (unit_fix, *) "keep_nonresonant =  F"
    write (unit_fix, *) ""
    write (unit_fix, *) "  grove"
    write (unit_fix, *) "    tree 3 7"
    write (unit_fix, *) "    tree 3 7"
    write (unit_fix, *) "      map 3 s_channel -24"
    write (unit_fix, *) "    tree 5 7"
    write (unit_fix, *) "    tree 3 7"
    write (unit_fix, *) "      map 3 s_channel -24"
    write (unit_fix, *) "      map 7 s_channel 23"
    write (unit_fix, *) "    tree 5 7"
    write (unit_fix, *) "      map 7 s_channel 25"
    write (unit_fix, *) "    tree 3 11"
    write (unit_fix, *) "      map 3 s_channel -24"
    close (unit_fix)

    write (u, "(A)")  "* Read phase-space file 'phs_forest_2.phs'"

    call syntax_phs_forest_init ()
    process_id = "foo"
    filename = "phs_forest_2.phs"
    call forest%read (filename, process_id, 2, 3, model, found_process)

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonance history set"
    write (u, "(A)")

    call forest%extract_resonance_history_set (res_set)
    call res_set%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()
    call forest%final ()
    call syntax_phs_forest_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_forest_2"

  end subroutine phs_forest_2


end module phs_forests_uti
