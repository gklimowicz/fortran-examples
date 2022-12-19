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

submodule (prc_template_me) prc_template_me_s

  use io_units
  use system_defs, only: TAB
  use diagnostics
  use flavors

  implicit none

contains

  module function template_me_def_type_string () result (string)
    type(string_t) :: string
    string = "template"
  end function template_me_def_type_string

  module subroutine template_me_def_write (object, unit)
    class(template_me_def_t), intent(in) :: object
    integer, intent(in) :: unit
    select type (writer => object%writer)
    type is (template_me_writer_t)
       call writer%write (unit)
    end select
  end subroutine template_me_def_write

  module subroutine template_me_def_read (object, unit)
    class(template_me_def_t), intent(out) :: object
    integer, intent(in) :: unit
    call msg_bug &
         ("WHIZARD template process definition: input not supported (yet)")
  end subroutine template_me_def_read

  module function template_me_def_needs_code () result (flag)
    logical :: flag
    flag = .true.
  end function template_me_def_needs_code

  module subroutine template_me_def_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (5))
    features = [ &
         var_str ("init"), &
         var_str ("update_alpha_s"), &
         var_str ("is_allowed"), &
         var_str ("new_event"), &
         var_str ("get_amplitude")]
  end subroutine template_me_def_get_features

  module subroutine template_me_def_connect (def, lib_driver, i, proc_driver)
    class(template_me_def_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
    integer(c_int) :: pid, fid
    type(c_funptr) :: fptr
    select type (proc_driver)
    type is  (template_me_driver_t)
       pid = i
       fid = 1
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%init)
       fid = 2
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%update_alpha_s)
       fid = 3
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%is_allowed)
       fid = 4
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%new_event)
       fid = 5
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%get_amplitude)
    end select
  end subroutine template_me_def_connect

  module function template_me_writer_type_name () result (string)
    type(string_t) :: string
    string = "template"
  end function template_me_writer_type_name

  module function template_me_writer_get_module_name (id) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: id
    name = "tpr_" // id
  end function template_me_writer_get_module_name

  module subroutine template_me_writer_write (object, unit)
    class(template_me_writer_t), intent(in) :: object
    integer, intent(in) :: unit
    integer :: i, j
    write (unit, "(5x,A,I0)") "# incoming part. = ", object%n_in
    write (unit, "(7x,A)", advance="no") &
                              "   Initial state: "
    do i = 1, object%n_in - 1
       write (unit, "(1x,A)", advance="no") char (object%prt_in(i))
    end do
    write (unit, "(1x,A)") char (object%prt_in(object%n_in))
    write (unit, "(5x,A,I0)") "# outgoing part. = ", object%n_out
    write (unit, "(7x,A)", advance="no") &
                              "   Final state:   "
    do j = 1, object%n_out - 1
       write (unit, "(1x,A)", advance="no") char (object%prt_out(j))
    end do
    write (unit, "(1x,A)") char (object%prt_out(object%n_out))
    write (unit, "(5x,A,I0)") "# part. (total) = ", object%n_tot
  end subroutine template_me_writer_write

  module subroutine template_me_writer_init (writer, model, &
       prt_in, prt_out, unity)
    class(template_me_writer_t), intent(out) :: writer
    class(model_data_t), intent(in), target :: model
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    logical, intent(in) :: unity
    writer%model => model
    writer%model_name = model%get_name ()
    writer%n_in = size (prt_in)
    writer%n_out = size (prt_out)
    writer%n_tot = size (prt_in) + size (prt_out)
    allocate (writer%prt_in (size (prt_in)), source = prt_in)
    allocate (writer%prt_out (size (prt_out)), source = prt_out)
    writer%unity = unity
  end subroutine template_me_writer_init

  module subroutine template_me_write_makefile_code &
       (writer, unit, id, os_data, verbose, testflag)
    class(template_me_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    write (unit, "(5A)")  "SOURCES += ", char (id), ".f90"
    write (unit, "(5A)")  "OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  "clean-", char (id), ":"
    if (verbose) then
       write (unit, "(5A)")  TAB, "rm -f tpr_", char (id), ".mod"
       write (unit, "(5A)")  TAB, "rm -f ", char (id), ".lo"
    else
       write (unit, "(5A)")  TAB // '@echo  "  RM        ', &
            trim (char (id)), '.mod"'
       write (unit, "(5A)")  TAB, "@rm -f tpr_", char (id), ".mod"
       write (unit, "(5A)")  TAB // '@echo  "  RM        ', &
            trim (char (id)), '.lo"'
       write (unit, "(5A)")  TAB, "@rm -f ", char (id), ".lo"
    end if
    write (unit, "(5A)")  "CLEAN_SOURCES += ", char (id), ".f90"
    write (unit, "(5A)")  "CLEAN_OBJECTS += tpr_", char (id), ".mod"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  char (id), ".lo: ", char (id), ".f90"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"
  end subroutine template_me_write_makefile_code

  module subroutine template_me_write_source_code (writer, id)
    class(template_me_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    integer, dimension(writer%n_in) :: prt_in, mult_in, col_in
    type(flavor_t), dimension(1:writer%n_in) :: flv_in
    integer, dimension(writer%n_out) :: prt_out, mult_out
    integer, dimension(writer%n_tot) :: prt, mult
    integer, dimension(:,:), allocatable :: sxxx
    integer :: dummy, status
    type(flavor_t), dimension(1:writer%n_out) :: flv_out
    type(string_t) :: proc_str, comment_str, col_str
    integer :: u, i, j
    integer :: hel, hel_in, hel_out, fac, factor, col_fac
    type(string_t) :: filename
    comment_str = ""
    do i = 1, writer%n_in
       comment_str = comment_str // writer%prt_in(i) // " "
    end do
    do j = 1, writer%n_out
       comment_str = comment_str // writer%prt_out(j) // " "
    end do
    do i = 1, writer%n_in
       prt_in(i) = writer%model%get_pdg (writer%prt_in(i))
       call flv_in(i)%init (prt_in(i), writer%model)
       mult_in(i) = flv_in(i)%get_multiplicity ()
       col_in(i) = abs (flv_in(i)%get_color_type ())
       mult(i) = mult_in(i)
       end do
    do j = 1, writer%n_out
       prt_out(j) = writer%model%get_pdg (writer%prt_out(j))
       call flv_out(j)%init (prt_out(j), writer%model)
       mult_out(j) = flv_out(j)%get_multiplicity ()
       mult(writer%n_in + j) = mult_out(j)
       end do
    prt(1:writer%n_in) = prt_in(1:writer%n_in)
    prt(writer%n_in+1:writer%n_tot) = prt_out(1:writer%n_out)
    proc_str = converter (prt)
    hel_in = product (mult_in)
    hel_out = product (mult_out)
    col_fac = product (col_in)
    hel = hel_in * hel_out
    fac = hel
    dummy = 1
    factor = 1
    if (writer%n_out >= 3) then
       do i = 3, writer%n_out
          factor = factor * (i - 2) * (i - 1)
       end do
    end if
    factor = factor * col_fac
    allocate (sxxx(1:hel,1:writer%n_tot))
    call create_spin_table (dummy,hel,fac,mult,sxxx)
    call msg_message ("Writing test matrix element for process '" &
         // char (id) // "'")
    filename = id // ".f90"
    u = free_unit ()
    open (unit=u, file=char(filename), action="write")
    write (u, "(A)") "! File generated automatically by WHIZARD"
    write (u, "(A)") "!                                        "
    write (u, "(A)") "! Note that irresp. of what you demanded WHIZARD"
    write (u, "(A)") "! treats this as colorless process       "
    write (u, "(A)") "!                                        "
    write (u, "(A)") "module tpr_" // char(id)
    write (u, "(A)") "                                         "
    write (u, "(A)") "  use kinds"
    write (u, "(A)") "  use omega_color, OCF => omega_color_factor"
    write (u, "(A)") "                                         "
    write (u, "(A)") "  implicit none"
    write (u, "(A)") "  private"
    write (u, "(A)") "                                         "
    write (u, "(A)") "  public :: md5sum"
    write (u, "(A)") "  public :: number_particles_in, number_particles_out"
    write (u, "(A)") "  public :: number_spin_states, spin_states"
    write (u, "(A)") "  public :: number_flavor_states, flavor_states"
    write (u, "(A)") "  public :: number_color_flows, color_flows"
    write (u, "(A)") "  public :: number_color_indices, number_color_factors, &"
    write (u, "(A)") "     color_factors, color_sum, openmp_supported"
    write (u, "(A)") "  public :: init, final, update_alpha_s"
    write (u, "(A)") "                                         "
    write (u, "(A)") "  public :: new_event, is_allowed, get_amplitude"
    write (u, "(A)") "       "
    write (u, "(A)") "  real(default), parameter :: &"
    write (u, "(A)") "       & conv = 0.38937966e12_default"
    write (u, "(A)") "       "
    write (u, "(A)") "  real(default), parameter :: &"
    write (u, "(A)") "       & pi = 3.1415926535897932384626433832795028841972_default"
    write (u, "(A)") "       "
    write (u, "(A)") "  real(default), parameter :: &"
    if (writer%unity) then
       write (u, "(A)") "                   & const = 1"
    else
       write (u, "(A,1x,I0,A)") "       & const = (16 * pi / conv) * " &
          // "(16 * pi**2)**(", writer%n_out, "-2) "
    end if
    write (u, "(A)") "       "
    write (u, "(A,1x,I0)") "  integer, parameter, private :: n_prt =  ", &
       writer%n_tot
    write (u, "(A,1x,I0)") "  integer, parameter, private :: n_in = ", &
       writer%n_in
    write (u, "(A,1x,I0)") "  integer, parameter, private :: n_out = ", &
       writer%n_out
    write (u, "(A)") "  integer, parameter, private :: n_cflow = 1"
    write (u, "(A)") "  integer, parameter, private :: n_cindex = 2"
    write (u, "(A)") "  !!! We ignore tensor products and take only one flavor state."
    write (u, "(A)") "  integer, parameter, private :: n_flv = 1"
    write (u, "(A,1x,I0)") "  integer, parameter, private :: n_hel = ", hel
    write (u, "(A)") "                                           "
    write (u, "(A)") "  logical, parameter, private :: T = .true."
    write (u, "(A)") "  logical, parameter, private :: F = .false."
    write (u, "(A)") "                                           "
    do i = 1, hel
       write (u, "(A)") "  integer, dimension(n_prt), parameter, private :: &"
       write (u, "(A)") "    " // s_conv(i) // " = [ " // &
            char(converter(sxxx(i,1:writer%n_tot))) // " ]"
    end do
    write (u, "(A)") "  integer, dimension(n_prt,n_hel), parameter, private :: table_spin_states = &"
    write (u, "(A)") "    reshape ( [ & "
    do i = 1, hel-1
       write (u, "(A)") "                 " // s_conv(i) // ", & "
    end do
    write (u, "(A)") "                 " // s_conv(hel) // " & "
    write (u, "(A)") "              ], [ n_prt, n_hel ] )"
    write (u, "(A)") "                                                 "
    write (u, "(A)") "  integer, dimension(n_prt), parameter, private :: &"
    write (u, "(A)") "    f0001 = [ " // char(proc_str) // " ]   !  " // char(comment_str)
    write (u, "(A)") "  integer, dimension(n_prt,n_flv), parameter, private :: table_flavor_states = &"
    write (u, "(A)") "    reshape ( [ f0001 ], [ n_prt, n_flv ] )"
    write (u, "(A)") "                                                 "
    write (u, "(A)") "  integer, dimension(n_cindex, n_prt), parameter, private :: &"
    !!! This produces non-matching color flows, better keep it completely colorless
    ! write (u, "(A)") "    c0001 = reshape ( [ " // char (dummy_colorizer (flv_in)) // &
    !                          " " // &
    select case (writer%n_in)
    case (1)
       col_str = "0,0,"
    case (2)
       col_str = "0,0,0,0,"
    end select
    write (u, "(A)") "    c0001 = reshape ( [" // char (col_str) // &
      (repeat ("0,0, ", writer%n_out-1)) // "0,0 ], " // " [ n_cindex, n_prt ] )"
    write (u, "(A)") "  integer, dimension(n_cindex, n_prt, n_cflow), parameter, private :: &"
    write (u, "(A)") "  table_color_flows = reshape ( [ c0001 ], [ n_cindex, n_prt, n_cflow ] )"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  logical, dimension(n_prt), parameter, private :: & "
    write (u, "(A)") "    g0001 = [ "  // (repeat ("F, ", writer%n_tot-1)) // "F ] "
    write (u, "(A)") "  logical, dimension(n_prt, n_cflow), parameter, private " &
         // ":: table_ghost_flags = &"
    write (u, "(A)") "    reshape ( [ g0001 ], [ n_prt, n_cflow ] )"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  integer, parameter, private :: n_cfactors = 1"
    write (u, "(A)") "  type(OCF), dimension(n_cfactors), parameter, private :: &"
    write (u, "(A)") "    table_color_factors = [  OCF(1,1,+1._default) ]"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  logical, dimension(n_flv), parameter, private :: a0001 = [ T ]"
    write (u, "(A)") "  logical, dimension(n_flv, n_cflow), parameter, private :: &"
    write (u, "(A)") "    flv_col_is_allowed = reshape ( [ a0001 ], [ n_flv, n_cflow ] )"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  complex(default), dimension (n_flv, n_hel, n_cflow), private, save :: amp"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  logical, dimension(n_hel), private, save :: hel_is_allowed = T"
    write (u, "(A)") "                                           "
    write (u, "(A)") "contains"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function md5sum ()"
    write (u, "(A)") "    character(len=32) :: md5sum"
    write (u, "(A)") "    ! DON'T EVEN THINK of modifying the following line!"
    write (u, "(A)") "    md5sum = """ // writer%md5sum // """"
    write (u, "(A)") "  end function md5sum"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  subroutine init (par, scheme)"
    write (u, "(A)") "    real(default), dimension(*), intent(in) :: par"
    write (u, "(A)") "    integer, intent(in) :: scheme"
    write (u, "(A)") "  end subroutine init"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  subroutine final ()"
    write (u, "(A)") "  end subroutine final"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  subroutine update_alpha_s (alpha_s)"
    write (u, "(A)") "    real(default), intent(in) :: alpha_s"
    write (u, "(A)") "  end subroutine update_alpha_s"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_particles_in () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = n_in"
    write (u, "(A)") "  end function number_particles_in"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_particles_out () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = n_out"
    write (u, "(A)") "  end function number_particles_out"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_spin_states () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = size (table_spin_states, dim=2)"
    write (u, "(A)") "  end function number_spin_states"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure subroutine spin_states (a)"
    write (u, "(A)") "    integer, dimension(:,:), intent(out) :: a"
    write (u, "(A)") "    a = table_spin_states"
    write (u, "(A)") "  end subroutine spin_states"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_flavor_states () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = 1"
    write (u, "(A)") "  end function number_flavor_states"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure subroutine flavor_states (a)"
    write (u, "(A)") "    integer, dimension(:,:), intent(out) :: a"
    write (u, "(A)") "    a = table_flavor_states"
    write (u, "(A)") "  end subroutine flavor_states"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_color_indices () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = size(table_color_flows, dim=1)"
    write (u, "(A)") "  end function number_color_indices"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure subroutine color_factors (cf)"
    write (u, "(A)") "    type(OCF), dimension(:), intent(out) :: cf"
    write (u, "(A)") "    cf = table_color_factors"
    write (u, "(A)") "  end subroutine color_factors"
    write (u, "(A)") "                                           "
    !pure unless OpenMP
    !write (u, "(A)") "  pure function color_sum (flv, hel) result (amp2)"
    write (u, "(A)") "  function color_sum (flv, hel) result (amp2)"
    write (u, "(A)") "    integer, intent(in) :: flv, hel"
    write (u, "(A)") "    real(kind=default) :: amp2"
    write (u, "(A)") "    amp2 = real (omega_color_sum (flv, hel, amp, table_color_factors))"
    write (u, "(A)") "  end function color_sum"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_color_flows () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = size (table_color_flows, dim=3)"
    write (u, "(A)") "  end function number_color_flows"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure subroutine color_flows (a, g)"
    write (u, "(A)") "    integer, dimension(:,:,:), intent(out) :: a"
    write (u, "(A)") "    logical, dimension(:,:), intent(out) :: g"
    write (u, "(A)") "    a = table_color_flows"
    write (u, "(A)") "    g = table_ghost_flags"
    write (u, "(A)") "  end subroutine color_flows"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function number_color_factors () result (n)"
    write (u, "(A)") "    integer :: n"
    write (u, "(A)") "    n = size (table_color_factors)"
    write (u, "(A)") "  end function number_color_factors"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function openmp_supported () result (status)"
    write (u, "(A)") "    logical :: status"
    write (u, "(A)") "    status = .false."
    write (u, "(A)") "  end function openmp_supported"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  subroutine new_event (p)"
    write (u, "(A)") "    real(default), dimension(0:3,*), intent(in) :: p"
    write (u, "(A)") "    call calculate_amplitudes (amp, p)"
    write (u, "(A)") "  end subroutine new_event"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function is_allowed (flv, hel, col) result (yorn)"
    write (u, "(A)") "    logical :: yorn"
    write (u, "(A)") "    integer, intent(in) :: flv, hel, col"
    write (u, "(A)") "    yorn = hel_is_allowed(hel) .and. flv_col_is_allowed(flv,col)"
    write (u, "(A)") "  end function is_allowed"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure function get_amplitude (flv, hel, col) result (amp_result)"
    write (u, "(A)") "    complex(default) :: amp_result"
    write (u, "(A)") "    integer, intent(in) :: flv, hel, col"
    write (u, "(A)") "    amp_result = amp (flv, hel, col)"
    write (u, "(A)") "  end function get_amplitude"
    write (u, "(A)") "                                           "
    write (u, "(A)") "  pure subroutine calculate_amplitudes (amp, k)"
    write (u, "(A)") "    complex(default), dimension(:,:,:), intent(out) :: amp"
    write (u, "(A)") "    real(default), dimension(0:3,*), intent(in) :: k"
    write (u, "(A)") "    real(default) :: fac"
    write (u, "(A)") "    integer :: i"
    write (u, "(A)") "    ! We give all helicities the same weight!"
    if (writer%unity) then
       write (u, "(A,1x,I0,1x,A)") "    fac = ", col_fac
       write (u, "(A)") "    amp = const * sqrt(fac)"
    else
       write (u, "(A,1x,I0,1x,A)") "    fac = ", factor
       write (u, "(A)") "    amp = sqrt((2 * (k(0,1)*k(0,2) &"
       write (u, "(A,1x,I0,A)") "         - dot_product (k(1:,1), k(1:,2)))) ** (3-", &
                                  writer%n_out, ")) * sqrt(const * fac)"
    end if
    write (u, "(A,1x,I0,A)") "    amp = amp / sqrt(", hel_out, "._default)"
    write (u, "(A)") "  end subroutine calculate_amplitudes"
    write (u, "(A)") "                                           "
    write (u, "(A)") "end module tpr_" // char(id)
    close (u, iostat=status)
    deallocate (sxxx)
  contains
    function s_conv (num) result (chrt)
      integer, intent(in) :: num
      character(len=10) :: chrt
      write (chrt, "(I10)") num
      chrt = trim(adjustl(chrt))
      if (num < 10) then
         chrt = "s000" // chrt
      else if (num < 100) then
         chrt = "s00" // chrt
      else if (num < 1000) then
         chrt = "s0" // chrt
      else
         chrt = "s" // chrt
      end if
    end function s_conv
    function converter (flv) result (str)
      integer, dimension(:), intent(in) :: flv
      type(string_t) :: str
      character(len=150), dimension(size(flv)) :: chrt
      integer :: i
      str = ""
      do i = 1, size(flv) - 1
         write (chrt(i), "(I10)") flv(i)
         str = str // var_str(trim(adjustl(chrt(i)))) // ", "
      end do
      write (chrt(size(flv)), "(I10)") flv(size(flv))
      str = str // trim(adjustl(chrt(size(flv))))
    end function converter
    integer function sj (j,m)
      integer, intent(in) :: j, m
      if (((j == 1) .and. (m == 1)) .or. &
          ((j == 2) .and. (m == 2)) .or. &
          ((j == 3) .and. (m == 3)) .or. &
          ((j == 4) .and. (m == 3)) .or. &
          ((j == 5) .and. (m == 4))) then
         sj = 1
      else if (((j == 2) .and. (m == 1)) .or. &
          ((j == 3) .and. (m == 1)) .or. &
          ((j == 4) .and. (m == 2)) .or. &
          ((j == 5) .and. (m == 2))) then
         sj = -1
      else if (((j == 3) .and. (m == 2)) .or. &
          ((j == 5) .and. (m == 3))) then
         sj = 0
      else if (((j == 4) .and. (m == 1)) .or. &
          ((j == 5) .and. (m == 1))) then
         sj = -2
      else if (((j == 4) .and. (m == 4)) .or. &
          ((j == 5) .and. (m == 5))) then
         sj = 2
      else
         call msg_fatal ("template_me_write_source_code: Wrong spin type")
      end if
    end function sj
    recursive subroutine create_spin_table (index, nhel, fac, mult, inta)
      integer, intent(inout) :: index, fac
      integer, intent(in) :: nhel
      integer, dimension(:), intent(in) :: mult
      integer, dimension(nhel,size(mult)), intent(out) :: inta
      integer :: j
      if (index > size(mult)) return
      fac = fac / mult(index)
      do j = 1, nhel
         inta(j,index) = sj (mult(index),mod(((j-1)/fac),mult(index))+1)
      end do
      index = index + 1
      call create_spin_table (index, nhel, fac, mult, inta)
    end subroutine create_spin_table
    function dummy_colorizer (flv) result (str)
      type(flavor_t), dimension(:), intent(in) :: flv
      type(string_t) :: str
      integer :: i, k
      str = ""
      k = 0
      do i = 1, size(flv)
         k = k + 1
         select case (flv(i)%get_color_type ())
         case (1,-1)
            str = str // "0,0, "
         case (3)
            str = str // int2string(k) // ",0, "
         case (-3)
            str = str // "0," // int2string(-k) // ", "
         case (8)
            str = str // int2string(k) // "," // int2string(-k-1) // ", "
            k = k + 1
         case default
            call msg_error ("Color type not supported.")
         end select
      end do
      str = adjustl(trim(str))
    end function dummy_colorizer
  end subroutine template_me_write_source_code

  module subroutine template_me_before_compile (writer, id)
    class(template_me_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine template_me_before_compile

  module subroutine template_me_after_compile (writer, id)
    class(template_me_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine template_me_after_compile

  module function template_me_writer_get_procname (feature) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: feature
    select case (char (feature))
    case ("n_in");   name = "number_particles_in"
    case ("n_out");  name = "number_particles_out"
    case ("n_flv");  name = "number_flavor_states"
    case ("n_hel");  name = "number_spin_states"
    case ("n_col");  name = "number_color_flows"
    case ("n_cin");  name = "number_color_indices"
    case ("n_cf");   name = "number_color_factors"
    case ("flv_state");  name = "flavor_states"
    case ("hel_state");  name = "spin_states"
    case ("col_state");  name = "color_flows"
    case default
       name = feature
    end select
  end function template_me_writer_get_procname

  module subroutine template_me_write_interface (writer, unit, id, feature)
    class(template_me_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(string_t), intent(in) :: feature
    type(string_t) :: name
    name = writer%get_c_procname (id, feature)
    write (unit, "(2x,9A)")  "interface"
    select case (char (feature))
    case ("init")
       write (unit, "(5x,9A)")  "subroutine ", char (name), &
            " (par, scheme) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), dimension(*), &
            &intent(in) :: par"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: scheme"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("update_alpha_s")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " (alpha_s) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), intent(in) :: alpha_s"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("is_allowed")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " &
            &(flv, hel, col, flag) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(7x,9A)")  "logical(c_bool), intent(out) :: flag"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("new_event")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " (p) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), dimension(0:3,*), &
            &intent(in) :: p"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("get_amplitude")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " &
            &(flv, hel, col, amp) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(7x,9A)")  "complex(c_default_complex), intent(out) &
            &:: amp"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    end select
    write (unit, "(2x,9A)")  "end interface"
  end subroutine template_me_write_interface

  module subroutine template_me_write_wrapper (writer, unit, id, feature)
    class(template_me_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    type(string_t) :: name
    name = writer%get_c_procname (id, feature)
    write (unit, *)
    select case (char (feature))
    case ("init")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (par, scheme) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use tpr_", char (id)
       write (unit, "(2x,9A)")  "real(c_default_float), dimension(*), &
            &intent(in) :: par"
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: scheme"
       if (c_default_float == default .and. c_int == kind(1)) then
          write (unit, "(2x,9A)")  "call ", char (feature), " (par, scheme)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("update_alpha_s")
       write (unit, "(9A)")  "subroutine ", char (name), " (alpha_s) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use tpr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), intent(in) &
               &:: alpha_s"
          write (unit, "(2x,9A)")  "call ", char (feature), " (alpha_s)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("is_allowed")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (flv, hel, col, flag) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use tpr_", char (id)
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(2x,9A)")  "logical(c_bool), intent(out) :: flag"
       write (unit, "(2x,9A)")  "flag = ", char (feature), &
            " (int (flv), int (hel), int (col))"
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("new_event")
       write (unit, "(9A)")  "subroutine ", char (name), " (p) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use tpr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), dimension(0:3,*), &
               &intent(in) :: p"
          write (unit, "(2x,9A)")  "call ", char (feature), " (p)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("get_amplitude")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (flv, hel, col, amp) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use tpr_", char (id)
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(2x,9A)")  "complex(c_default_complex), intent(out) &
            &:: amp"
       write (unit, "(2x,9A)")  "amp = ", char (feature), &
            " (int (flv), int (hel), int (col))"
       write (unit, "(9A)")  "end subroutine ", char (name)
    end select
  end subroutine template_me_write_wrapper

  module function template_me_driver_type_name () result (string)
    type(string_t) :: string
    string = "template"
  end function template_me_driver_type_name

  module subroutine template_me_state_write (object, unit)
    class(template_me_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,L1)")  "Template ME state: new kinematics = ", &
         object%new_kinematics
  end subroutine template_me_state_write

  module subroutine template_me_state_reset_new_kinematics (object)
    class(template_me_state_t), intent(inout) :: object
  end subroutine template_me_state_reset_new_kinematics

  module subroutine prc_template_me_write (object, unit)
    class(prc_template_me_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(3x,A)", advance="no")  "Template process core:"
    if (object%data_known) then
       write (u, "(1x,A)")  char (object%data%id)
    else
       write (u, "(1x,A)")  "[undefined]"
    end if
    if (allocated (object%par)) then
       write (u, "(3x,A)")  "Parameter array:"
       do i = 1, size (object%par)
          write (u, "(5x,I0,1x,ES17.10)")  i, object%par(i)
       end do
    end if
  end subroutine prc_template_me_write

  module subroutine prc_template_me_write_name (object, unit)
    class(prc_template_me_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: template"
  end subroutine prc_template_me_write_name

  module subroutine prc_template_me_set_parameters (prc_template_me, model)
    class(prc_template_me_t), intent(inout) :: prc_template_me
    class(model_data_t), intent(in), target, optional :: model
    if (present (model)) then
       if (.not. allocated (prc_template_me%par)) &
            allocate (prc_template_me%par (model%get_n_real ()))
       call model%real_parameters_to_array (prc_template_me%par)
       prc_template_me%scheme = model%get_scheme_num ()
    end if
  end subroutine prc_template_me_set_parameters

  module subroutine prc_template_me_init (object, def, lib, id, i_component)
    class(prc_template_me_t), intent(inout) :: object
    class(prc_core_def_t), intent(in), target :: def
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    call object%base_init (def, lib, id, i_component)
    call object%activate_parameters ()
  end subroutine prc_template_me_init

  module subroutine prc_template_me_activate_parameters (object)
    class (prc_template_me_t), intent(inout) :: object
    if (allocated (object%driver)) then
       if (allocated (object%par)) then
          select type (driver => object%driver)
          type is (template_me_driver_t)
             if (associated (driver%init)) then
                call driver%init (object%par, object%scheme)
             end if
          end select
       else
          call msg_bug ("prc_template_me_activate: parameter set is not allocated")
       end if
    else
       call msg_bug ("prc_template_me_activate: driver is not allocated")
    end if
  end subroutine prc_template_me_activate_parameters

  module function prc_template_me_is_allowed &
      (object, i_term, f, h, c) result (flag)
    class(prc_template_me_t), intent(in) :: object
    integer, intent(in) :: i_term, f, h, c
    logical :: flag
    logical(c_bool) :: cflag
    select type (driver => object%driver)
    type is (template_me_driver_t)
       call driver%is_allowed (f, h, c, cflag)
       flag = cflag
    end select
  end function prc_template_me_is_allowed

  module subroutine prc_template_me_compute_hard_kinematics &
       (object, p_seed, i_term, int_hard, core_state)
    class(prc_template_me_t), intent(in) :: object
    type(vector4_t), dimension(:), intent(in) :: p_seed
    integer, intent(in) :: i_term
    type(interaction_t), intent(inout) :: int_hard
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    call int_hard%set_momenta (p_seed)
  end subroutine prc_template_me_compute_hard_kinematics

  module subroutine prc_template_me_compute_eff_kinematics &
       (object, i_term, int_hard, int_eff, core_state)
    class(prc_template_me_t), intent(in) :: object
    integer, intent(in) :: i_term
    type(interaction_t), intent(in) :: int_hard
    type(interaction_t), intent(inout) :: int_eff
    class(prc_core_state_t), intent(inout), allocatable :: core_state
  end subroutine prc_template_me_compute_eff_kinematics

  module function prc_template_me_compute_amplitude &
       (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
       core_state)  result (amp)
    class(prc_template_me_t), intent(in) :: object
    integer, intent(in) :: j
    type(vector4_t), dimension(:), intent(in) :: p
    integer, intent(in) :: f, h, c
    real(default), intent(in) :: fac_scale, ren_scale
    real(default), intent(in), allocatable :: alpha_qcd_forced
    class(prc_core_state_t), intent(inout), allocatable, optional :: core_state
    complex(default) :: amp
    integer :: n_tot, i
    real(c_default_float), dimension(:,:), allocatable :: parray
    complex(c_default_complex) :: camp
    logical :: new_event
    select type (driver => object%driver)
    type is (template_me_driver_t)
       new_event = .true.
       if (present (core_state)) then
          if (allocated (core_state)) then
             select type (core_state)
             type is (template_me_state_t)
                new_event = core_state%new_kinematics
                core_state%new_kinematics = .false.
             end select
          end if
       end if
       if (new_event) then
          n_tot = object%data%n_in + object%data%n_out
          allocate (parray (0:3, n_tot))
          forall (i = 1:n_tot)  parray(:,i) = vector4_get_components (p(i))
          call driver%new_event (parray)
       end if
       if (object%is_allowed (1, f, h, c)) then
          call driver%get_amplitude &
               (int (f, c_int), int (h, c_int), int (c, c_int), camp)
          amp = camp
       else
          amp = 0
       end if
    end select
  end function prc_template_me_compute_amplitude


end submodule prc_template_me_s

