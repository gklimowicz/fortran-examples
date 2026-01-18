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

submodule (analysis) analysis_s

  use io_units
  use format_utils, only: quote_underscore, tex_format
  use system_defs, only: TAB
  use diagnostics
  use ifiles

  implicit none

contains

  module subroutine graph_options_init (graph_options)
    class(graph_options_t), intent(out) :: graph_options
    graph_options%id = ""
    graph_options%title = ""
    graph_options%description = ""
    graph_options%x_label = ""
    graph_options%y_label = ""
    graph_options%gmlcode_bg = ""
    graph_options%gmlcode_fg = ""
  end subroutine graph_options_init

  module subroutine graph_options_set (graph_options, id, &
       title, description, x_label, y_label, width_mm, height_mm, &
       x_log, y_log, x_min, x_max, y_min, y_max, &
       gmlcode_bg, gmlcode_fg)
    class(graph_options_t), intent(inout) :: graph_options
    type(string_t), intent(in), optional :: id
    type(string_t), intent(in), optional :: title
    type(string_t), intent(in), optional :: description
    type(string_t), intent(in), optional :: x_label, y_label
    integer, intent(in), optional :: width_mm, height_mm
    logical, intent(in), optional :: x_log, y_log
    real(default), intent(in), optional :: x_min, x_max, y_min, y_max
    type(string_t), intent(in), optional :: gmlcode_bg, gmlcode_fg
    if (present (id))  graph_options%id = id
    if (present (title))  graph_options%title = title
    if (present (description))  graph_options%description = description
    if (present (x_label))  graph_options%x_label = x_label
    if (present (y_label))  graph_options%y_label = y_label
    if (present (width_mm))   graph_options%width_mm  = width_mm
    if (present (height_mm))  graph_options%height_mm = height_mm
    if (present (x_log))  graph_options%x_log = x_log
    if (present (y_log))  graph_options%y_log = y_log
    if (present (x_min))  graph_options%x_min = x_min
    if (present (x_max))  graph_options%x_max = x_max
    if (present (y_min))  graph_options%y_min = y_min
    if (present (y_max))  graph_options%y_max = y_max
    if (present (x_min))  graph_options%x_min_set = .true.
    if (present (x_max))  graph_options%x_max_set = .true.
    if (present (y_min))  graph_options%y_min_set = .true.
    if (present (y_max))  graph_options%y_max_set = .true.
    if (present (gmlcode_bg))  graph_options%gmlcode_bg = gmlcode_bg
    if (present (gmlcode_fg))  graph_options%gmlcode_fg = gmlcode_fg
  end subroutine graph_options_set

  module subroutine graph_options_write (gro, unit)
    class(graph_options_t), intent(in) :: gro
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
1   format (A,1x,'"',A,'"')
2   format (A,1x,L1)
3   format (A,1x,ES19.12)
4   format (A,1x,I0)
5   format (A,1x,'[undefined]')
    write (u, 1)  "title       =", char (gro%title)
    write (u, 1)  "description =", char (gro%description)
    write (u, 1)  "x_label     =", char (gro%x_label)
    write (u, 1)  "y_label     =", char (gro%y_label)
    write (u, 2)  "x_log       =", gro%x_log
    write (u, 2)  "y_log       =", gro%y_log
    if (gro%x_min_set) then
       write (u, 3)  "x_min       =", gro%x_min
    else
       write (u, 5)  "x_min       ="
    end if
    if (gro%x_max_set) then
       write (u, 3)  "x_max       =", gro%x_max
    else
       write (u, 5)  "x_max       ="
    end if
    if (gro%y_min_set) then
       write (u, 3)  "y_min       =", gro%y_min
    else
       write (u, 5)  "y_min       ="
    end if
    if (gro%y_max_set) then
       write (u, 3)  "y_max       =", gro%y_max
    else
       write (u, 5)  "y_max       ="
    end if
    write (u, 4)  "width_mm    =", gro%width_mm
    write (u, 4)  "height_mm   =", gro%height_mm
    write (u, 1)  "gmlcode_bg  =", char (gro%gmlcode_bg)
    write (u, 1)  "gmlcode_fg  =", char (gro%gmlcode_fg)
  end subroutine graph_options_write

  subroutine graph_options_write_tex_header (gro, unit)
    type(graph_options_t), intent(in) :: gro
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    if (gro%title /= "") then
       write (u, "(A)")
       write (u, "(A)")  "\section{" // char (gro%title) // "}"
    else
       write (u, "(A)")  "\section{" // char (quote_underscore (gro%id)) // "}"
    end if
    if (gro%description /= "") then
       write (u, "(A)")  char (gro%description)
       write (u, *)
       write (u, "(A)")  "\vspace*{\baselineskip}"
    end if
    write (u, "(A)")  "\vspace*{\baselineskip}"
    write (u, "(A)")  "\unitlength 1mm"
    write (u, "(A,I0,',',I0,A)")  &
         "\begin{gmlgraph*}(", &
         gro%width_mm, gro%height_mm, &
         ")[dat]"
  end subroutine graph_options_write_tex_header

  subroutine graph_options_write_tex_footer (gro, unit)
    type(graph_options_t), intent(in) :: gro
    integer, intent(in), optional :: unit
    integer :: u, width, height
    width = gro%width_mm - 10
    height = gro%height_mm - 10
    u = given_output_unit (unit)
    write (u, "(A)")  "  begingmleps ""Whizard-Logo.eps"";"
    write (u, "(A,I0,A,I0,A)")  &
         "    base := (", width, "*unitlength,", height, "*unitlength);"
    write (u, "(A)")  "    height := 9.6*unitlength;"
    write (u, "(A)")  "    width := 11.2*unitlength;"
    write (u, "(A)")  "  endgmleps;"
    write (u, "(A)")  "\end{gmlgraph*}"
  end subroutine graph_options_write_tex_footer

  function graph_options_get_id (gro) result (id)
    type(string_t) :: id
    type(graph_options_t), intent(in) :: gro
    id = gro%id
  end function graph_options_get_id

  function graph_options_get_gml_setup (gro) result (cmd)
    type(string_t) :: cmd
    type(graph_options_t), intent(in) :: gro
    type(string_t) :: x_str, y_str
    if (gro%x_log) then
       x_str = "log"
    else
       x_str = "linear"
    end if
    if (gro%y_log) then
       y_str = "log"
    else
       y_str = "linear"
    end if
    cmd = "setup (" // x_str // ", " // y_str // ");"
  end function graph_options_get_gml_setup

  function graph_options_get_gml_x_label (gro) result (cmd)
    type(string_t) :: cmd
    type(graph_options_t), intent(in) :: gro
    cmd = 'label.bot (<' // '<' // gro%x_label // '>' // '>, out);'
  end function graph_options_get_gml_x_label

  function graph_options_get_gml_y_label (gro) result (cmd)
    type(string_t) :: cmd
    type(graph_options_t), intent(in) :: gro
    cmd = 'label.ulft (<' // '<' // gro%y_label // '>' // '>, out);'
  end function graph_options_get_gml_y_label

  function graph_options_get_gml_graphrange &
       (gro, x_min, x_max, y_min, y_max) result (cmd)
    type(string_t) :: cmd
    type(graph_options_t), intent(in) :: gro
    real(default), intent(in), optional :: x_min, x_max, y_min, y_max
    type(string_t) :: x_min_str, x_max_str, y_min_str, y_max_str
    character(*), parameter :: fmt = "(ES15.8)"
    if (gro%x_min_set) then
       x_min_str = "#" // trim (adjustl (real2string (gro%x_min, fmt)))
    else if (present (x_min)) then
       x_min_str = "#" // trim (adjustl (real2string (x_min, fmt)))
    else
       x_min_str = "??"
    end if
    if (gro%x_max_set) then
       x_max_str = "#" // trim (adjustl (real2string (gro%x_max, fmt)))
    else if (present (x_max)) then
       x_max_str = "#" // trim (adjustl (real2string (x_max, fmt)))
    else
       x_max_str = "??"
    end if
    if (gro%y_min_set) then
       y_min_str = "#" // trim (adjustl (real2string (gro%y_min, fmt)))
    else if (present (y_min)) then
       y_min_str = "#" // trim (adjustl (real2string (y_min, fmt)))
    else
       y_min_str = "??"
    end if
    if (gro%y_max_set) then
       y_max_str = "#" // trim (adjustl (real2string (gro%y_max, fmt)))
    else if (present (y_max)) then
       y_max_str = "#" // trim (adjustl (real2string (y_max, fmt)))
    else
       y_max_str = "??"
    end if
    cmd = "graphrange (" // x_min_str // ", " // y_min_str // "), " &
         // "(" // x_max_str // ", " // y_max_str // ");"
  end function graph_options_get_gml_graphrange

  function graph_options_get_gml_bg_command (gro) result (cmd)
    type(string_t) :: cmd
    type(graph_options_t), intent(in) :: gro
    cmd = gro%gmlcode_bg
  end function graph_options_get_gml_bg_command

  function graph_options_get_gml_fg_command (gro) result (cmd)
    type(string_t) :: cmd
    type(graph_options_t), intent(in) :: gro
    cmd = gro%gmlcode_fg
  end function graph_options_get_gml_fg_command

  subroutine graph_options_get_header (pl, header, comment)
    type(graph_options_t), intent(in) :: pl
    type(ifile_t), intent(inout) :: header
    type(string_t), intent(in), optional :: comment
    type(string_t) :: c
    if (present (comment)) then
       c = comment
    else
       c = ""
    end if
    call ifile_append (header, &
         c // "ID: " // pl%id)
    call ifile_append (header, &
         c // "title: " // pl%title)
    call ifile_append (header, &
         c // "description: " // pl%description)
    call ifile_append (header, &
         c // "x axis label: " // pl%x_label)
    call ifile_append (header, &
         c // "y axis label: " // pl%y_label)
  end subroutine graph_options_get_header

  module subroutine drawing_options_write (dro, unit)
    class(drawing_options_t), intent(in) :: dro
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
1   format (A,1x,'"',A,'"')
2   format (A,1x,L1)
    write (u, 2)  "with_hbars  =", dro%with_hbars
    write (u, 2)  "with_base   =", dro%with_base
    write (u, 2)  "piecewise   =", dro%piecewise
    write (u, 2)  "fill        =", dro%fill
    write (u, 2)  "draw        =", dro%draw
    write (u, 2)  "err         =", dro%err
    write (u, 2)  "symbols     =", dro%symbols
    write (u, 1)  "fill_options=", char (dro%fill_options)
    write (u, 1)  "draw_options=", char (dro%draw_options)
    write (u, 1)  "err_options =", char (dro%err_options)
    write (u, 1)  "symbol      =", char (dro%symbol)
    write (u, 1)  "gmlcode_bg  =", char (dro%gmlcode_bg)
    write (u, 1)  "gmlcode_fg  =", char (dro%gmlcode_fg)
  end subroutine drawing_options_write

  module subroutine drawing_options_init_histogram (dro)
    class(drawing_options_t), intent(out) :: dro
    dro%dataset = "dat"
    dro%with_hbars = .true.
    dro%with_base = .true.
    dro%piecewise = .true.
    dro%fill = .true.
    dro%draw = .true.
    dro%fill_options = "withcolor col.default"
    dro%draw_options = ""
    dro%err_options = ""
    dro%symbol = "fshape(circle scaled 1mm)()"
    dro%gmlcode_bg = ""
    dro%gmlcode_fg = ""
  end subroutine drawing_options_init_histogram

  module subroutine drawing_options_init_plot (dro)
    class(drawing_options_t), intent(out) :: dro
    dro%dataset = "dat"
    dro%draw = .true.
    dro%fill_options = "withcolor col.default"
    dro%draw_options = ""
    dro%err_options = ""
    dro%symbol = "fshape(circle scaled 1mm)()"
    dro%gmlcode_bg = ""
    dro%gmlcode_fg = ""
  end subroutine drawing_options_init_plot

  module subroutine drawing_options_set (dro, dataset, &
       with_hbars, with_base, piecewise, fill, draw, err, symbols, &
       fill_options, draw_options, err_options, symbol, &
       gmlcode_bg, gmlcode_fg)
    class(drawing_options_t), intent(inout) :: dro
    type(string_t), intent(in), optional :: dataset
    logical, intent(in), optional :: with_hbars, with_base, piecewise
    logical, intent(in), optional :: fill, draw, err, symbols
    type(string_t), intent(in), optional :: fill_options, draw_options
    type(string_t), intent(in), optional :: err_options, symbol
    type(string_t), intent(in), optional :: gmlcode_bg, gmlcode_fg
    if (present (dataset))  dro%dataset = dataset
    if (present (with_hbars))  dro%with_hbars = with_hbars
    if (present (with_base))  dro%with_base = with_base
    if (present (piecewise))  dro%piecewise = piecewise
    if (present (fill))  dro%fill = fill
    if (present (draw))  dro%draw = draw
    if (present (err))  dro%err = err
    if (present (symbols))  dro%symbols = symbols
    if (present (fill_options))  dro%fill_options = fill_options
    if (present (draw_options))  dro%draw_options = draw_options
    if (present (err_options))  dro%err_options = err_options
    if (present (symbol))  dro%symbol = symbol
    if (present (gmlcode_bg))  dro%gmlcode_bg = gmlcode_bg
    if (present (gmlcode_fg))  dro%gmlcode_fg = gmlcode_fg
  end subroutine drawing_options_set

  function drawing_options_get_calc_command (dro) result (cmd)
    type(string_t) :: cmd
    type(drawing_options_t), intent(in) :: dro
    if (dro%with_base) then
       cmd = "calculate " // dro%dataset // ".base (" // dro%dataset // ") " &
            // "(x, #0);"
    else
       cmd = ""
    end if
  end function drawing_options_get_calc_command

  function drawing_options_get_draw_command (dro) result (cmd)
    type(string_t) :: cmd
    type(drawing_options_t), intent(in) :: dro
    if (dro%fill) then
       cmd = "fill"
    else if (dro%draw) then
       cmd = "draw"
    else
       cmd = ""
    end if
    if (dro%fill .or. dro%draw) then
       if (dro%piecewise)  cmd = cmd // " piecewise"
       if (dro%draw .and. dro%with_base)  cmd = cmd // " cyclic"
       cmd = cmd // " from (" // dro%dataset
       if (dro%with_base) then
          if (dro%piecewise) then
             cmd = cmd // ", " // dro%dataset // ".base/\"  ! "
          else
             cmd = cmd // " ~ " // dro%dataset // ".base\"  ! "
          end if
       end if
       cmd = cmd // ")"
       if (dro%fill) then
          cmd = cmd // " " // dro%fill_options
          if (dro%draw)  cmd = cmd // " outlined"
       end if
       if (dro%draw)  cmd = cmd // " " // dro%draw_options
       cmd = cmd // ";"
    end if
  end function drawing_options_get_draw_command

  function drawing_options_get_err_command (dro) result (cmd)
    type(string_t) :: cmd
    type(drawing_options_t), intent(in) :: dro
    if (dro%err) then
       cmd = "draw piecewise " &
            // "from (" // dro%dataset // ".err)" &
            // " " // dro%err_options // ";"
    else
       cmd = ""
    end if
  end function drawing_options_get_err_command

  function drawing_options_get_symb_command (dro) result (cmd)
    type(string_t) :: cmd
    type(drawing_options_t), intent(in) :: dro
    if (dro%symbols) then
       cmd = "phantom" &
            // " from (" // dro%dataset // ")" &
            // " withsymbol (" // dro%symbol // ");"
    else
       cmd = ""
    end if
  end function drawing_options_get_symb_command

  function drawing_options_get_gml_bg_command (dro) result (cmd)
    type(string_t) :: cmd
    type(drawing_options_t), intent(in) :: dro
    cmd = dro%gmlcode_bg
  end function drawing_options_get_gml_bg_command

  function drawing_options_get_gml_fg_command (dro) result (cmd)
    type(string_t) :: cmd
    type(drawing_options_t), intent(in) :: dro
    cmd = dro%gmlcode_fg
  end function drawing_options_get_gml_fg_command

  subroutine observable_init (obs, obs_label, obs_unit, graph_options)
    type(observable_t), intent(out) :: obs
    type(string_t), intent(in), optional :: obs_label, obs_unit
    type(graph_options_t), intent(in), optional :: graph_options
    if (present (obs_label)) then
       obs%obs_label = obs_label
    else
       obs%obs_label = ""
    end if
    if (present (obs_unit)) then
       obs%obs_unit = obs_unit
    else
       obs%obs_unit = ""
    end if
    if (present (graph_options)) then
       obs%graph_options = graph_options
    else
       call obs%graph_options%init ()
    end if
  end subroutine observable_init

  subroutine observable_clear (obs)
    type(observable_t), intent(inout) :: obs
    obs%sum_values = 0
    obs%sum_squared_values = 0
    obs%sum_weights = 0
    obs%sum_squared_weights = 0
    obs%count = 0
  end subroutine observable_clear

  module subroutine observable_record_value_unweighted (obs, value, success)
    type(observable_t), intent(inout) :: obs
    real(default), intent(in) :: value
    logical, intent(out), optional :: success
    obs%sum_values = obs%sum_values + value
    obs%sum_squared_values = obs%sum_squared_values + value**2
    obs%sum_weights = obs%sum_weights + 1
    obs%sum_squared_weights = obs%sum_squared_weights + 1
    obs%count = obs%count + 1
    if (present (success))  success = .true.
  end subroutine observable_record_value_unweighted

  module subroutine observable_record_value_weighted (obs, value, weight, success)
    type(observable_t), intent(inout) :: obs
    real(default), intent(in) :: value, weight
    logical, intent(out), optional :: success
    obs%sum_values = obs%sum_values + value * weight
    obs%sum_squared_values = obs%sum_squared_values + value**2 * weight
    obs%sum_weights = obs%sum_weights + weight
    obs%sum_squared_weights = obs%sum_squared_weights + weight**2
    obs%count = obs%count + 1
    if (present (success))  success = .true.
  end subroutine observable_record_value_weighted

  function observable_get_n_entries (obs) result (n)
    integer :: n
    type(observable_t), intent(in) :: obs
    n = obs%count
  end function observable_get_n_entries

  function observable_get_average (obs) result (avg)
    real(default) :: avg
    type(observable_t), intent(in) :: obs
    if (obs%sum_weights /= 0) then
       avg = obs%sum_values / obs%sum_weights
    else
       avg = 0
    end if
  end function observable_get_average

  function observable_get_error (obs) result (err)
    real(default) :: err
    type(observable_t), intent(in) :: obs
    real(default) :: var, n
    if (obs%sum_weights /= 0) then
       select case (obs%count)
       case (0:1)
          err = 0
       case default
          n = obs%count
          var = obs%sum_squared_values / obs%sum_weights &
                - (obs%sum_values / obs%sum_weights) ** 2
          err = sqrt (max (var, 0._default) / (n - 1))
       end select
    else
       err = 0
    end if
  end function observable_get_error

  function observable_get_label (obs, wl, wu) result (string)
    type(string_t) :: string
    type(observable_t), intent(in) :: obs
    logical, intent(in) :: wl, wu
    type(string_t) :: obs_label, obs_unit
    if (wl) then
       if (obs%obs_label /= "") then
          obs_label = obs%obs_label
       else
          obs_label = "\textrm{Observable}"
       end if
    else
       obs_label = ""
    end if
    if (wu) then
       if (obs%obs_unit /= "") then
          if (wl) then
             obs_unit = "\;[" // obs%obs_unit // "]"
          else
             obs_unit = obs%obs_unit
          end if
       else
          obs_unit = ""
       end if
    else
       obs_unit = ""
    end if
    string = obs_label // obs_unit
  end function observable_get_label

  subroutine observable_write (obs, unit)
    type(observable_t), intent(in) :: obs
    integer, intent(in), optional :: unit
    real(default) :: avg, err, relerr
    integer :: n
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    avg = observable_get_average (obs)
    err = observable_get_error (obs)
    if (avg /= 0) then
       relerr = err / abs (avg)
    else
       relerr = 0
    end if
    n = observable_get_n_entries (obs)
    if (obs%graph_options%title /= "") then
       write (u, "(A,1x,3A)") &
            "title       =", '"', char (obs%graph_options%title), '"'
    end if
    if (obs%graph_options%title /= "") then
       write (u, "(A,1x,3A)") &
            "description =", '"', char (obs%graph_options%description), '"'
    end if
    write (u, "(A,1x," // HISTOGRAM_DATA_FORMAT // ")", advance = "no") &
         "average     =", avg
    call write_unit ()
    write (u, "(A,1x," // HISTOGRAM_DATA_FORMAT // ")", advance = "no") &
         "error[abs]  =", err
    call write_unit ()
    write (u, "(A,1x," // HISTOGRAM_DATA_FORMAT // ")") &
         "error[rel]  =", relerr
    write (u, "(A,1x,I0)") &
         "n_entries   =", n
  contains
    subroutine write_unit ()
      if (obs%obs_unit /= "") then
         write (u, "(1x,A)")  char (obs%obs_unit)
      else
         write (u, *)
      end if
    end subroutine write_unit
  end subroutine observable_write

  subroutine observable_write_driver (obs, unit, write_heading)
    type(observable_t), intent(in) :: obs
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: write_heading
    real(default) :: avg, err
    integer :: n_digits
    logical :: heading
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    heading = .true.;  if (present (write_heading))  heading = write_heading
    avg = observable_get_average (obs)
    err = observable_get_error (obs)
    if (avg /= 0 .and. err /= 0) then
       n_digits = max (2, 2 - int (log10 (abs (err / real (avg, default)))))
    else if (avg /= 0) then
       n_digits = 100
    else
       n_digits = 1
    end if
    if (heading) then
       write (u, "(A)")
       if (obs%graph_options%title /= "") then
          write (u, "(A)")  "\section{" // char (obs%graph_options%title) &
               // "}"
       else
          write (u, "(A)")  "\section{Observable}"
       end if
       if (obs%graph_options%description /= "") then
          write (u, "(A)")  char (obs%graph_options%description)
          write (u, *)
       end if
       write (u, "(A)")  "\begin{flushleft}"
    end if
    write (u, "(A)", advance="no")  "  $\langle{" ! $ sign
    write (u, "(A)", advance="no")  char (observable_get_label (obs, wl=.true., wu=.false.))
    write (u, "(A)", advance="no")  "}\rangle = "
    write (u, "(A)", advance="no")  char (tex_format (avg, n_digits))
    write (u, "(A)", advance="no")  "\pm"
    write (u, "(A)", advance="no")  char (tex_format (err, 2))
    write (u, "(A)", advance="no")  "\;{"
    write (u, "(A)", advance="no")  char (observable_get_label (obs, wl=.false., wu=.true.))
    write (u, "(A)")  "}"
    write (u, "(A)", advance="no")  "     \quad[n_{\text{entries}} = "
    write (u, "(I0)",advance="no")  observable_get_n_entries (obs)
    write (u, "(A)")  "]$"          ! $ fool Emacs' noweb mode
    if (heading) then
       write (u, "(A)") "\end{flushleft}"
    end if
  end subroutine observable_write_driver

  subroutine bin_init (bin, midpoint, width)
    type(bin_t), intent(out) :: bin
    real(default), intent(in) :: midpoint, width
    bin%midpoint = midpoint
    bin%width = width
  end subroutine bin_init

  elemental subroutine bin_clear (bin)
    type(bin_t), intent(inout) :: bin
    bin%sum_weights = 0
    bin%sum_squared_weights = 0
    bin%sum_excess_weights = 0
    bin%count = 0
  end subroutine bin_clear

  subroutine bin_record_value (bin, normalize, weight, excess)
    type(bin_t), intent(inout) :: bin
    logical, intent(in) :: normalize
    real(default), intent(in) :: weight
    real(default), intent(in), optional :: excess
    real(default) :: w, e
    if (normalize) then
       if (bin%width /= 0) then
          w = weight / bin%width
          if (present (excess))  e = excess / bin%width
       else
          w = 0
          if (present (excess))  e = 0
       end if
    else
       w = weight
       if (present (excess))  e = excess
    end if
    bin%sum_weights = bin%sum_weights + w
    bin%sum_squared_weights = bin%sum_squared_weights + w ** 2
    if (present (excess)) &
         bin%sum_excess_weights = bin%sum_excess_weights + abs (e)
    bin%count = bin%count + 1
  end subroutine bin_record_value

  function bin_get_midpoint (bin) result (x)
    real(default) :: x
    type(bin_t), intent(in) :: bin
    x = bin%midpoint
  end function bin_get_midpoint

  function bin_get_width (bin) result (w)
    real(default) :: w
    type(bin_t), intent(in) :: bin
    w = bin%width
  end function bin_get_width

  function bin_get_n_entries (bin) result (n)
    integer :: n
    type(bin_t), intent(in) :: bin
    n = bin%count
  end function bin_get_n_entries

  function bin_get_sum (bin) result (s)
    real(default) :: s
    type(bin_t), intent(in) :: bin
    s = bin%sum_weights
  end function bin_get_sum

  function bin_get_error (bin) result (err)
    real(default) :: err
    type(bin_t), intent(in) :: bin
    err = sqrt (bin%sum_squared_weights)
  end function bin_get_error

  function bin_get_excess (bin) result (excess)
    real(default) :: excess
    type(bin_t), intent(in) :: bin
    excess = bin%sum_excess_weights
  end function bin_get_excess

  subroutine bin_write_header (unit)
    integer, intent(in), optional :: unit
    character(120) :: buffer
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (buffer, "(A,4(1x," //HISTOGRAM_HEAD_FORMAT // "),2x,A)") &
         "#", "bin midpoint", "value    ", "error    ", &
         "excess     ", "n"
    write (u, "(A)")  trim (buffer)
  end subroutine bin_write_header

  subroutine bin_write (bin, unit)
    type(bin_t), intent(in) :: bin
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,4(1x," // HISTOGRAM_DATA_FORMAT // "),2x,I0)") &
         bin_get_midpoint (bin), &
         bin_get_sum (bin), &
         bin_get_error (bin), &
         bin_get_excess (bin), &
         bin_get_n_entries (bin)
  end subroutine bin_write

  module subroutine histogram_init_n_bins (h, id, &
       lower_bound, upper_bound, n_bins, normalize_bins, &
       obs_label, obs_unit, graph_options, drawing_options)
    type(histogram_t), intent(out) :: h
    type(string_t), intent(in) :: id
    real(default), intent(in) :: lower_bound, upper_bound
    integer, intent(in) :: n_bins
    logical, intent(in) :: normalize_bins
    type(string_t), intent(in), optional :: obs_label, obs_unit
    type(graph_options_t), intent(in), optional :: graph_options
    type(drawing_options_t), intent(in), optional :: drawing_options
    real(default) :: bin_width
    integer :: i
    call observable_init (h%obs_within_bounds, obs_label, obs_unit)
    call observable_init (h%obs, obs_label, obs_unit)
    h%lower_bound = lower_bound
    h%upper_bound = upper_bound
    h%n_bins = max (n_bins, 1)
    h%width = h%upper_bound - h%lower_bound
    h%normalize_bins = normalize_bins
    bin_width = h%width / h%n_bins
    allocate (h%bin (h%n_bins))
    call bin_init (h%underflow, h%lower_bound, 0._default)
    do i = 1, h%n_bins
       call bin_init (h%bin(i), &
            h%lower_bound - bin_width/2 + i * bin_width, bin_width)
    end do
    call bin_init (h%overflow, h%upper_bound, 0._default)
    if (present (graph_options)) then
       h%graph_options = graph_options
    else
       call h%graph_options%init ()
    end if
    call graph_options_set (h%graph_options, id = id)
    if (present (drawing_options)) then
       h%drawing_options = drawing_options
    else
       call h%drawing_options%init_histogram ()
    end if
  end subroutine histogram_init_n_bins

  module subroutine histogram_init_bin_width (h, id, &
       lower_bound, upper_bound, bin_width, normalize_bins, &
       obs_label, obs_unit, graph_options, drawing_options)
    type(histogram_t), intent(out) :: h
    type(string_t), intent(in) :: id
    real(default), intent(in) :: lower_bound, upper_bound, bin_width
    logical, intent(in) :: normalize_bins
    type(string_t), intent(in), optional :: obs_label, obs_unit
    type(graph_options_t), intent(in), optional :: graph_options
    type(drawing_options_t), intent(in), optional :: drawing_options
    integer :: n_bins
    if (bin_width /= 0) then
       n_bins = nint ((upper_bound - lower_bound) / bin_width)
    else
       n_bins = 1
    end if
    call histogram_init_n_bins (h, id, &
         lower_bound, upper_bound, n_bins, normalize_bins, &
         obs_label, obs_unit, graph_options, drawing_options)
  end subroutine histogram_init_bin_width

  subroutine histogram_init_histogram (h, h_in, drawing_options)
    type(histogram_t), intent(out) :: h
    type(histogram_t), intent(in) :: h_in
    type(drawing_options_t), intent(in), optional :: drawing_options
    h = h_in
    if (present (drawing_options)) then
       h%drawing_options = drawing_options
    end if
  end subroutine histogram_init_histogram

  subroutine histogram_clear (h)
    type(histogram_t), intent(inout) :: h
    call observable_clear (h%obs)
    call observable_clear (h%obs_within_bounds)
    call bin_clear (h%underflow)
    if (allocated (h%bin))  call bin_clear (h%bin)
    call bin_clear (h%overflow)
  end subroutine histogram_clear

  subroutine histogram_record_value_unweighted (h, value, excess, success)
    type(histogram_t), intent(inout) :: h
    real(default), intent(in) :: value
    real(default), intent(in), optional :: excess
    logical, intent(out), optional :: success
    integer :: i_bin
    call observable_record_value (h%obs, value)
    if (h%width /= 0) then
       i_bin = floor (((value - h%lower_bound) / h%width) * h%n_bins) + 1
    else
       i_bin = 0
    end if
    if (i_bin <= 0) then
       call bin_record_value (h%underflow, .false., 1._default, excess)
       if (present (success))  success = .false.
    else if (i_bin <= h%n_bins) then
       call observable_record_value (h%obs_within_bounds, value)
       call bin_record_value &
            (h%bin(i_bin), h%normalize_bins, 1._default, excess)
       if (present (success))  success = .true.
    else
       call bin_record_value (h%overflow, .false., 1._default, excess)
       if (present (success))  success = .false.
    end if
  end subroutine histogram_record_value_unweighted

  subroutine histogram_record_value_weighted (h, value, weight, success)
    type(histogram_t), intent(inout) :: h
    real(default), intent(in) :: value, weight
    logical, intent(out), optional :: success
    integer :: i_bin
    call observable_record_value (h%obs, value, weight)
    if (h%width /= 0) then
       i_bin = floor (((value - h%lower_bound) / h%width) * h%n_bins) + 1
    else
       i_bin = 0
    end if
    if (i_bin <= 0) then
       call bin_record_value (h%underflow, .false., weight)
       if (present (success))  success = .false.
    else if (i_bin <= h%n_bins) then
       call observable_record_value (h%obs_within_bounds, value, weight)
       call bin_record_value (h%bin(i_bin), h%normalize_bins, weight)
       if (present (success))  success = .true.
    else
       call bin_record_value (h%overflow, .false., weight)
       if (present (success))  success = .false.
    end if
  end subroutine histogram_record_value_weighted

  function histogram_get_n_entries (h) result (n)
    integer :: n
    type(histogram_t), intent(in) :: h
    n = observable_get_n_entries (h%obs)
  end function histogram_get_n_entries

  function histogram_get_average (h) result (avg)
    real(default) :: avg
    type(histogram_t), intent(in) :: h
    avg = observable_get_average (h%obs)
  end function histogram_get_average

  function histogram_get_error (h) result (err)
    real(default) :: err
    type(histogram_t), intent(in) :: h
    err = observable_get_error (h%obs)
  end function histogram_get_error

  function histogram_get_n_entries_within_bounds (h) result (n)
    integer :: n
    type(histogram_t), intent(in) :: h
    n = observable_get_n_entries (h%obs_within_bounds)
  end function histogram_get_n_entries_within_bounds

  function histogram_get_average_within_bounds (h) result (avg)
    real(default) :: avg
    type(histogram_t), intent(in) :: h
    avg = observable_get_average (h%obs_within_bounds)
  end function histogram_get_average_within_bounds

  function histogram_get_error_within_bounds (h) result (err)
    real(default) :: err
    type(histogram_t), intent(in) :: h
    err = observable_get_error (h%obs_within_bounds)
  end function histogram_get_error_within_bounds

  function histogram_get_n_bins (h) result (n)
    type(histogram_t), intent(in) :: h
    integer :: n
    n = h%n_bins
  end function histogram_get_n_bins

  function histogram_get_n_entries_for_bin (h, i) result (n)
    integer :: n
    type(histogram_t), intent(in) :: h
    integer, intent(in) :: i
    if (i <= 0) then
       n = bin_get_n_entries (h%underflow)
    else if (i <= h%n_bins) then
       n = bin_get_n_entries (h%bin(i))
    else
       n = bin_get_n_entries (h%overflow)
    end if
  end function histogram_get_n_entries_for_bin

  function histogram_get_sum_for_bin (h, i) result (avg)
    real(default) :: avg
    type(histogram_t), intent(in) :: h
    integer, intent(in) :: i
    if (i <= 0) then
       avg = bin_get_sum (h%underflow)
    else if (i <= h%n_bins) then
       avg = bin_get_sum (h%bin(i))
    else
       avg = bin_get_sum (h%overflow)
    end if
  end function histogram_get_sum_for_bin

  function histogram_get_error_for_bin (h, i) result (err)
    real(default) :: err
    type(histogram_t), intent(in) :: h
    integer, intent(in) :: i
    if (i <= 0) then
       err = bin_get_error (h%underflow)
    else if (i <= h%n_bins) then
       err = bin_get_error (h%bin(i))
    else
       err = bin_get_error (h%overflow)
    end if
  end function histogram_get_error_for_bin

  function histogram_get_excess_for_bin (h, i) result (err)
    real(default) :: err
    type(histogram_t), intent(in) :: h
    integer, intent(in) :: i
    if (i <= 0) then
       err = bin_get_excess (h%underflow)
    else if (i <= h%n_bins) then
       err = bin_get_excess (h%bin(i))
    else
       err = bin_get_excess (h%overflow)
    end if
  end function histogram_get_excess_for_bin

  function histogram_get_graph_options_ptr (h) result (ptr)
    type(graph_options_t), pointer :: ptr
    type(histogram_t), intent(in), target :: h
    ptr => h%graph_options
  end function histogram_get_graph_options_ptr

  function histogram_get_drawing_options_ptr (h) result (ptr)
    type(drawing_options_t), pointer :: ptr
    type(histogram_t), intent(in), target :: h
    ptr => h%drawing_options
  end function histogram_get_drawing_options_ptr

  subroutine histogram_write (h, unit)
    type(histogram_t), intent(in) :: h
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    call bin_write_header (u)
    if (allocated (h%bin)) then
       do i = 1, h%n_bins
          call bin_write (h%bin(i), u)
       end do
    end if
    write (u, "(A)")
    write (u, "(A,1x,A)")  "#", "Underflow:"
    call bin_write (h%underflow, u)
    write (u, "(A)")
    write (u, "(A,1x,A)")  "#", "Overflow:"
    call bin_write (h%overflow, u)
    write (u, "(A)")
    write (u, "(A,1x,A)")  "#", "Summary: data within bounds"
    call observable_write (h%obs_within_bounds, u)
    write (u, "(A)")
    write (u, "(A,1x,A)")  "#", "Summary: all data"
    call observable_write (h%obs, u)
    write (u, "(A)")
  end subroutine histogram_write

  subroutine histogram_write_gml_reader (h, filename, unit)
    type(histogram_t), intent(in) :: h
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: unit
    character(*), parameter :: fmt = "(ES15.8)"
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(2x,A)")  'fromfile "' // char (filename) // '":'
    write (u, "(4x,A)")  'key "# Histogram:";'
    write (u, "(4x,A)")  'dx := #' &
         // real2char (h%width / h%n_bins / 2, fmt) // ';'
    write (u, "(4x,A)")  'for i withinblock:'
    write (u, "(6x,A)")  'get x, y, y.d, y.n, y.e;'
    if (h%drawing_options%with_hbars) then
       write (u, "(6x,A)")  'plot (' // char (h%drawing_options%dataset) &
            // ') (x,y) hbar dx;'
    else
       write (u, "(6x,A)")  'plot (' // char (h%drawing_options%dataset) &
            // ') (x,y);'
    end if
    if (h%drawing_options%err) then
       write (u, "(6x,A)")  'plot (' // char (h%drawing_options%dataset) &
         // '.err) ' &
            // '(x,y) vbar y.d;'
    end if
    !!! Future excess options for plots
    ! write (u, "(6x,A)")  'if show_excess: ' // &
    !            & 'plot(dat.e)(x, y plus y.e) hbar dx; fi'
    write (u, "(4x,A)")  'endfor'
    write (u, "(2x,A)")  'endfrom'
  end subroutine histogram_write_gml_reader

  subroutine histogram_write_gml_driver (h, filename, unit)
    type(histogram_t), intent(in) :: h
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: unit
    type(string_t) :: calc_cmd, bg_cmd, draw_cmd, err_cmd, symb_cmd, fg_cmd
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    call graph_options_write_tex_header (h%graph_options, unit)
    write (u, "(2x,A)")  char (graph_options_get_gml_setup (h%graph_options))
    write (u, "(2x,A)")  char (graph_options_get_gml_graphrange &
         (h%graph_options, x_min=h%lower_bound, x_max=h%upper_bound))
    call histogram_write_gml_reader (h, filename, unit)
    calc_cmd = drawing_options_get_calc_command (h%drawing_options)
    if (calc_cmd /= "")  write (u, "(2x,A)")  char (calc_cmd)
    bg_cmd = drawing_options_get_gml_bg_command (h%drawing_options)
    if (bg_cmd /= "")  write (u, "(2x,A)")  char (bg_cmd)
    draw_cmd = drawing_options_get_draw_command (h%drawing_options)
    if (draw_cmd /= "")  write (u, "(2x,A)")  char (draw_cmd)
    err_cmd = drawing_options_get_err_command (h%drawing_options)
    if (err_cmd /= "")  write (u, "(2x,A)")  char (err_cmd)
    symb_cmd = drawing_options_get_symb_command (h%drawing_options)
    if (symb_cmd /= "")  write (u, "(2x,A)")  char (symb_cmd)
    fg_cmd = drawing_options_get_gml_fg_command (h%drawing_options)
    if (fg_cmd /= "")  write (u, "(2x,A)")  char (fg_cmd)
    write (u, "(2x,A)")  char (graph_options_get_gml_x_label (h%graph_options))
    write (u, "(2x,A)")  char (graph_options_get_gml_y_label (h%graph_options))
    call graph_options_write_tex_footer (h%graph_options, unit)
    write (u, "(A)") "\vspace*{2\baselineskip}"
    write (u, "(A)") "\begin{flushleft}"
    write (u, "(A)") "\textbf{Data within bounds:} \\"
    call observable_write_driver (h%obs_within_bounds, unit, &
                                  write_heading=.false.)
    write (u, "(A)") "\\[0.5\baselineskip]"
    write (u, "(A)") "\textbf{All data:} \\"
    call observable_write_driver (h%obs, unit, write_heading=.false.)
    write (u, "(A)") "\end{flushleft}"
  end subroutine histogram_write_gml_driver

  subroutine histogram_get_header (h, header, comment)
    type(histogram_t), intent(in) :: h
    type(ifile_t), intent(inout) :: header
    type(string_t), intent(in), optional :: comment
    type(string_t) :: c
    if (present (comment)) then
       c = comment
    else
       c = ""
    end if
    call ifile_append (header, c // "WHIZARD histogram data")
    call graph_options_get_header (h%graph_options, header, comment)
    call ifile_append (header, &
         c // "range: " // real2string (h%lower_bound) &
         // " - " // real2string (h%upper_bound))
    call ifile_append (header, &
         c // "counts total: " &
         // int2char (histogram_get_n_entries_within_bounds (h)))
    call ifile_append (header, &
         c // "total average: " &
         // real2string (histogram_get_average_within_bounds (h)) // " +- " &
         // real2string (histogram_get_error_within_bounds (h)))
  end subroutine histogram_get_header

  module subroutine point_init_contents (point, x, y, yerr, xerr)
    type(point_t), intent(out) :: point
    real(default), intent(in) :: x, y
    real(default), intent(in), optional :: yerr, xerr
    point%x = x
    point%y = y
    if (present (yerr))  point%yerr = yerr
    if (present (xerr))  point%xerr = xerr
  end subroutine point_init_contents

  module subroutine point_init_point (point, point_in)
    type(point_t), intent(out) :: point
    type(point_t), intent(in) :: point_in
    point%x = point_in%x
    point%y = point_in%y
    point%yerr = point_in%yerr
    point%xerr = point_in%xerr
  end subroutine point_init_point

  function point_get_x (point) result (x)
    real(default) :: x
    type(point_t), intent(in) :: point
    x = point%x
  end function point_get_x

  function point_get_y (point) result (y)
    real(default) :: y
    type(point_t), intent(in) :: point
    y = point%y
  end function point_get_y

  function point_get_xerr (point) result (xerr)
    real(default) :: xerr
    type(point_t), intent(in) :: point
    xerr = point%xerr
  end function point_get_xerr

  function point_get_yerr (point) result (yerr)
    real(default) :: yerr
    type(point_t), intent(in) :: point
    yerr = point%yerr
  end function point_get_yerr

  subroutine point_write_header (unit)
    integer, intent(in) :: unit
    character(120) :: buffer
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (buffer, "(A,4(1x," // HISTOGRAM_HEAD_FORMAT // "))") &
         "#", "x       ", "y       ", "yerr     ", "xerr     "
    write (u, "(A)")  trim (buffer)
  end subroutine point_write_header

  subroutine point_write (point, unit)
    type(point_t), intent(in) :: point
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,4(1x," // HISTOGRAM_DATA_FORMAT // "))") &
         point_get_x (point), &
         point_get_y (point), &
         point_get_yerr (point), &
         point_get_xerr (point)
  end subroutine point_write

  module subroutine plot_init_empty (p, id, graph_options, drawing_options)
    type(plot_t), intent(out) :: p
    type(string_t), intent(in) :: id
    type(graph_options_t), intent(in), optional :: graph_options
    type(drawing_options_t), intent(in), optional :: drawing_options
    if (present (graph_options)) then
       p%graph_options = graph_options
    else
       call p%graph_options%init ()
    end if
    call p%graph_options%set (id = id)
    if (present (drawing_options)) then
       p%drawing_options = drawing_options
    else
       call p%drawing_options%init_plot ()
    end if
  end subroutine plot_init_empty

  module subroutine plot_init_plot (p, p_in, drawing_options)
    type(plot_t), intent(out) :: p
    type(plot_t), intent(in) :: p_in
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(point_t), pointer :: current, new
    current => p_in%first
    do while (associated (current))
       allocate (new)
       call point_init (new, current)
       if (associated (p%last)) then
          p%last%next => new
       else
          p%first => new
       end if
       p%last => new
       current => current%next
    end do
    p%count = p_in%count
    p%graph_options = p_in%graph_options
    if (present (drawing_options)) then
       p%drawing_options = drawing_options
    else
       p%drawing_options = p_in%drawing_options
    end if
  end subroutine plot_init_plot

  subroutine plot_final (plot)
    type(plot_t), intent(inout) :: plot
    type(point_t), pointer :: current
    do while (associated (plot%first))
       current => plot%first
       plot%first => current%next
       deallocate (current)
    end do
    plot%last => null ()
  end subroutine plot_final

  subroutine plot_clear (plot)
    type(plot_t), intent(inout) :: plot
    plot%count = 0
    call plot_final (plot)
  end subroutine plot_clear

  subroutine plot_record_value (plot, x, y, yerr, xerr, success)
    type(plot_t), intent(inout) :: plot
    real(default), intent(in) :: x, y
    real(default), intent(in), optional :: yerr, xerr
    logical, intent(out), optional :: success
    type(point_t), pointer :: point
    plot%count = plot%count + 1
    allocate (point)
    call point_init (point, x, y, yerr, xerr)
    if (associated (plot%first)) then
       plot%last%next => point
    else
       plot%first => point
    end if
    plot%last => point
    if (present (success))  success = .true.
  end subroutine plot_record_value

  function plot_get_n_entries (plot) result (n)
    integer :: n
    type(plot_t), intent(in) :: plot
    n = plot%count
  end function plot_get_n_entries

  function plot_get_graph_options_ptr (p) result (ptr)
    type(graph_options_t), pointer :: ptr
    type(plot_t), intent(in), target :: p
    ptr => p%graph_options
  end function plot_get_graph_options_ptr

  function plot_get_drawing_options_ptr (p) result (ptr)
    type(drawing_options_t), pointer :: ptr
    type(plot_t), intent(in), target :: p
    ptr => p%drawing_options
  end function plot_get_drawing_options_ptr

  subroutine plot_write (plot, unit)
    type(plot_t), intent(in) :: plot
    integer, intent(in), optional :: unit
    type(point_t), pointer :: point
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    call point_write_header (u)
    point => plot%first
    do while (associated (point))
       call point_write (point, unit)
       point => point%next
    end do
    write (u, *)
    write (u, "(A,1x,A)")  "#", "Summary:"
    write (u, "(A,1x,I0)") &
         "n_entries =", plot_get_n_entries (plot)
    write (u, *)
  end subroutine plot_write

  subroutine plot_write_gml_reader (p, filename, unit)
    type(plot_t), intent(in) :: p
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(2x,A)")  'fromfile "' // char (filename) // '":'
    write (u, "(4x,A)")  'key "# Plot:";'
    write (u, "(4x,A)")  'for i withinblock:'
    write (u, "(6x,A)")  'get x, y, y.err, x.err;'
    write (u, "(6x,A)")  'plot (' // char (p%drawing_options%dataset) &
         // ') (x,y);'
    if (p%drawing_options%err) then
       write (u, "(6x,A)")  'plot (' // char (p%drawing_options%dataset) &
            // '.err) (x,y) vbar y.err hbar x.err;'
    end if
    write (u, "(4x,A)")  'endfor'
    write (u, "(2x,A)")  'endfrom'
  end subroutine plot_write_gml_reader

  subroutine plot_write_gml_driver (p, filename, unit)
    type(plot_t), intent(in) :: p
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: unit
    type(string_t) :: calc_cmd, bg_cmd, draw_cmd, err_cmd, symb_cmd, fg_cmd
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    call graph_options_write_tex_header (p%graph_options, unit)
    write (u, "(2x,A)") &
         char (graph_options_get_gml_setup (p%graph_options))
    write (u, "(2x,A)") &
         char (graph_options_get_gml_graphrange (p%graph_options))
    call plot_write_gml_reader (p, filename, unit)
    calc_cmd = drawing_options_get_calc_command (p%drawing_options)
    if (calc_cmd /= "")  write (u, "(2x,A)")  char (calc_cmd)
    bg_cmd = drawing_options_get_gml_bg_command (p%drawing_options)
    if (bg_cmd /= "")  write (u, "(2x,A)")  char (bg_cmd)
    draw_cmd = drawing_options_get_draw_command (p%drawing_options)
    if (draw_cmd /= "")  write (u, "(2x,A)")  char (draw_cmd)
    err_cmd = drawing_options_get_err_command (p%drawing_options)
    if (err_cmd /= "")  write (u, "(2x,A)")  char (err_cmd)
    symb_cmd = drawing_options_get_symb_command (p%drawing_options)
    if (symb_cmd /= "")  write (u, "(2x,A)")  char (symb_cmd)
    fg_cmd = drawing_options_get_gml_fg_command (p%drawing_options)
    if (fg_cmd /= "")  write (u, "(2x,A)")  char (fg_cmd)
    write (u, "(2x,A)")  char (graph_options_get_gml_x_label (p%graph_options))
    write (u, "(2x,A)")  char (graph_options_get_gml_y_label (p%graph_options))
    call graph_options_write_tex_footer (p%graph_options, unit)
  end subroutine plot_write_gml_driver

  subroutine plot_get_header (plot, header, comment)
    type(plot_t), intent(in) :: plot
    type(ifile_t), intent(inout) :: header
    type(string_t), intent(in), optional :: comment
    type(string_t) :: c
    if (present (comment)) then
       c = comment
    else
       c = ""
    end if
    call ifile_append (header, c // "WHIZARD plot data")
    call graph_options_get_header (plot%graph_options, header, comment)
    call ifile_append (header, &
         c // "number of points: " &
         // int2char (plot_get_n_entries (plot)))
  end subroutine plot_get_header

  subroutine graph_element_final (el)
    type(graph_element_t), intent(inout) :: el
    select case (el%type)
    case (AN_HISTOGRAM)
       deallocate (el%h)
    case (AN_PLOT)
       call plot_final (el%p)
       deallocate (el%p)
    end select
    el%type = AN_UNDEFINED
  end subroutine graph_element_final

  function graph_element_get_n_entries (el) result (n)
    integer :: n
    type(graph_element_t), intent(in) :: el
    select case (el%type)
    case (AN_HISTOGRAM);  n = histogram_get_n_entries (el%h)
    case (AN_PLOT);       n = plot_get_n_entries (el%p)
    case default;         n = 0
    end select
  end function graph_element_get_n_entries

  function graph_element_get_graph_options_ptr (el) result (ptr)
    type(graph_options_t), pointer :: ptr
    type(graph_element_t), intent(in) :: el
    select case (el%type)
    case (AN_HISTOGRAM);  ptr => histogram_get_graph_options_ptr (el%h)
    case (AN_PLOT);       ptr => plot_get_graph_options_ptr (el%p)
    case default;         ptr => null ()
    end select
  end function graph_element_get_graph_options_ptr

  function graph_element_get_drawing_options_ptr (el) result (ptr)
    type(drawing_options_t), pointer :: ptr
    type(graph_element_t), intent(in) :: el
    select case (el%type)
    case (AN_HISTOGRAM);  ptr => histogram_get_drawing_options_ptr (el%h)
    case (AN_PLOT);       ptr => plot_get_drawing_options_ptr (el%p)
    case default;         ptr => null ()
    end select
  end function graph_element_get_drawing_options_ptr

  subroutine graph_element_write (el, unit)
    type(graph_element_t), intent(in) :: el
    integer, intent(in), optional :: unit
    type(graph_options_t), pointer :: gro
    type(string_t) :: id
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    gro => graph_element_get_graph_options_ptr (el)
    id = graph_options_get_id (gro)
    write (u, "(A,A)")  '#', repeat ("-", 78)
    select case (el%type)
    case (AN_HISTOGRAM)
       write (u, "(A)", advance="no")  "# Histogram: "
       write (u, "(1x,A)")  char (id)
       call histogram_write (el%h, unit)
    case (AN_PLOT)
       write (u, "(A)", advance="no")  "# Plot: "
       write (u, "(1x,A)")  char (id)
       call plot_write (el%p, unit)
    end select
  end subroutine graph_element_write

  subroutine graph_element_write_gml_reader (el, filename, unit)
    type(graph_element_t), intent(in) :: el
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: unit
    select case (el%type)
    case (AN_HISTOGRAM);  call histogram_write_gml_reader (el%h, filename, unit)
    case (AN_PLOT);       call plot_write_gml_reader (el%p, filename, unit)
    end select
  end subroutine graph_element_write_gml_reader

  subroutine graph_init (g, id, n_elements, graph_options)
    type(graph_t), intent(out) :: g
    type(string_t), intent(in) :: id
    integer, intent(in) :: n_elements
    type(graph_options_t), intent(in), optional :: graph_options
    allocate (g%el (n_elements))
    if (present (graph_options)) then
       g%graph_options = graph_options
    else
       call g%graph_options%init ()
    end if
    call g%graph_options%set (id = id)
  end subroutine graph_init

  subroutine graph_insert_histogram (g, i, h, drawing_options)
    type(graph_t), intent(inout), target :: g
    integer, intent(in) :: i
    type(histogram_t), intent(in) :: h
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(graph_options_t), pointer :: gro
    type(drawing_options_t), pointer :: dro
    type(string_t) :: id
    g%el(i)%type = AN_HISTOGRAM
    allocate (g%el(i)%h)
    call histogram_init_histogram (g%el(i)%h, h, drawing_options)
    gro => histogram_get_graph_options_ptr (g%el(i)%h)
    dro => histogram_get_drawing_options_ptr (g%el(i)%h)
    id = graph_options_get_id (gro)
    call dro%set (dataset = "dat." // id)
  end subroutine graph_insert_histogram

  subroutine graph_insert_plot (g, i, p, drawing_options)
    type(graph_t), intent(inout) :: g
    integer, intent(in) :: i
    type(plot_t), intent(in) :: p
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(graph_options_t), pointer :: gro
    type(drawing_options_t), pointer :: dro
    type(string_t) :: id
    g%el(i)%type = AN_PLOT
    allocate (g%el(i)%p)
    call plot_init_plot (g%el(i)%p, p, drawing_options)
    gro => plot_get_graph_options_ptr (g%el(i)%p)
    dro => plot_get_drawing_options_ptr (g%el(i)%p)
    id = graph_options_get_id (gro)
    call dro%set (dataset = "dat." // id)
  end subroutine graph_insert_plot

  subroutine graph_final (g)
    type(graph_t), intent(inout) :: g
    integer :: i
    do i = 1, size (g%el)
       call graph_element_final (g%el(i))
    end do
    deallocate (g%el)
  end subroutine graph_final

  function graph_get_n_elements (graph) result (n)
    integer :: n
    type(graph_t), intent(in) :: graph
    n = size (graph%el)
  end function graph_get_n_elements

  function graph_get_drawing_options_ptr (g, i) result (ptr)
    type(drawing_options_t), pointer :: ptr
    type(graph_t), intent(in), target :: g
    integer, intent(in) :: i
    ptr => graph_element_get_drawing_options_ptr (g%el(i))
  end function graph_get_drawing_options_ptr

  subroutine graph_write (graph, unit)
    type(graph_t), intent(in) :: graph
    integer, intent(in), optional :: unit
    integer :: i
    do i = 1, size (graph%el)
       call graph_element_write (graph%el(i), unit)
    end do
  end subroutine graph_write

  subroutine graph_write_gml_driver (g, filename, unit)
    type(graph_t), intent(in) :: g
    type(string_t), intent(in) :: filename
    type(string_t) :: calc_cmd, bg_cmd, draw_cmd, err_cmd, symb_cmd, fg_cmd
    integer, intent(in), optional :: unit
    type(drawing_options_t), pointer :: dro
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    call graph_options_write_tex_header (g%graph_options, unit)
    write (u, "(2x,A)") &
         char (graph_options_get_gml_setup (g%graph_options))
    write (u, "(2x,A)") &
         char (graph_options_get_gml_graphrange (g%graph_options))
    do i = 1, size (g%el)
       call graph_element_write_gml_reader (g%el(i), filename, unit)
       calc_cmd = drawing_options_get_calc_command &
                       (graph_element_get_drawing_options_ptr (g%el(i)))
       if (calc_cmd /= "")  write (u, "(2x,A)")  char (calc_cmd)
    end do
    bg_cmd = graph_options_get_gml_bg_command (g%graph_options)
    if (bg_cmd /= "")  write (u, "(2x,A)")  char (bg_cmd)
    do i = 1, size (g%el)
       dro => graph_element_get_drawing_options_ptr (g%el(i))
       bg_cmd = drawing_options_get_gml_bg_command (dro)
       if (bg_cmd /= "")  write (u, "(2x,A)")  char (bg_cmd)
       draw_cmd = drawing_options_get_draw_command (dro)
       if (draw_cmd /= "")  write (u, "(2x,A)")  char (draw_cmd)
       err_cmd = drawing_options_get_err_command (dro)
       if (err_cmd /= "")  write (u, "(2x,A)")  char (err_cmd)
       symb_cmd = drawing_options_get_symb_command (dro)
       if (symb_cmd /= "")  write (u, "(2x,A)")  char (symb_cmd)
       fg_cmd = drawing_options_get_gml_fg_command (dro)
       if (fg_cmd /= "")  write (u, "(2x,A)")  char (fg_cmd)
    end do
    fg_cmd = graph_options_get_gml_fg_command (g%graph_options)
    if (fg_cmd /= "")  write (u, "(2x,A)")  char (fg_cmd)
    write (u, "(2x,A)")  char (graph_options_get_gml_x_label (g%graph_options))
    write (u, "(2x,A)")  char (graph_options_get_gml_y_label (g%graph_options))
    call graph_options_write_tex_footer (g%graph_options, unit)
  end subroutine graph_write_gml_driver

  subroutine graph_get_header (graph, header, comment)
    type(graph_t), intent(in) :: graph
    type(ifile_t), intent(inout) :: header
    type(string_t), intent(in), optional :: comment
    type(string_t) :: c
    if (present (comment)) then
       c = comment
    else
       c = ""
    end if
    call ifile_append (header, c // "WHIZARD graph data")
    call graph_options_get_header (graph%graph_options, header, comment)
    call ifile_append (header, &
         c // "number of graph elements: " &
         // int2char (graph_get_n_elements (graph)))
  end subroutine graph_get_header

  subroutine analysis_object_init (obj, id, type)
    type(analysis_object_t), intent(out) :: obj
    type(string_t), intent(in) :: id
    integer, intent(in) :: type
    obj%id = id
    obj%type = type
    select case (obj%type)
    case (AN_OBSERVABLE);  allocate (obj%obs)
    case (AN_HISTOGRAM);   allocate (obj%h)
    case (AN_PLOT);        allocate (obj%p)
    case (AN_GRAPH);       allocate (obj%g)
    end select
  end subroutine analysis_object_init

  subroutine analysis_object_final (obj)
    type(analysis_object_t), intent(inout) :: obj
    select case (obj%type)
    case (AN_OBSERVABLE)
       deallocate (obj%obs)
    case (AN_HISTOGRAM)
       deallocate (obj%h)
    case (AN_PLOT)
       call plot_final (obj%p)
       deallocate (obj%p)
    case (AN_GRAPH)
       call graph_final (obj%g)
       deallocate (obj%g)
    end select
    obj%type = AN_UNDEFINED
  end subroutine analysis_object_final

  subroutine analysis_object_clear (obj)
    type(analysis_object_t), intent(inout) :: obj
    select case (obj%type)
    case (AN_OBSERVABLE)
       call observable_clear (obj%obs)
    case (AN_HISTOGRAM)
       call histogram_clear (obj%h)
    case (AN_PLOT)
       call plot_clear (obj%p)
    end select
  end subroutine analysis_object_clear

  subroutine analysis_object_record_data (obj, &
       x, y, yerr, xerr, weight, excess, success)
    type(analysis_object_t), intent(inout) :: obj
    real(default), intent(in) :: x
    real(default), intent(in), optional :: y, yerr, xerr, weight, excess
    logical, intent(out), optional :: success
    select case (obj%type)
    case (AN_OBSERVABLE)
       if (present (weight)) then
          call observable_record_value_weighted (obj%obs, x, weight, success)
       else
          call observable_record_value_unweighted (obj%obs, x, success)
       end if
    case (AN_HISTOGRAM)
       if (present (weight)) then
          call histogram_record_value_weighted (obj%h, x, weight, success)
       else
          call histogram_record_value_unweighted (obj%h, x, excess, success)
       end if
    case (AN_PLOT)
       if (present (y)) then
          call plot_record_value (obj%p, x, y, yerr, xerr, success)
       else
          if (present (success))  success = .false.
       end if
    case default
       if (present (success))  success = .false.
    end select
  end subroutine analysis_object_record_data

  subroutine analysis_object_set_next_ptr (obj, next)
    type(analysis_object_t), intent(inout) :: obj
    type(analysis_object_t), pointer :: next
    obj%next => next
  end subroutine analysis_object_set_next_ptr

  function analysis_object_get_next_ptr (obj) result (next)
    type(analysis_object_t), pointer :: next
    type(analysis_object_t), intent(in) :: obj
    next => obj%next
  end function analysis_object_get_next_ptr

  function analysis_object_get_n_elements (obj) result (n)
    integer :: n
    type(analysis_object_t), intent(in) :: obj
    select case (obj%type)
    case (AN_HISTOGRAM)
       n = 1
    case (AN_PLOT)
       n = 1
    case (AN_GRAPH)
       n = graph_get_n_elements (obj%g)
    case default
       n = 0
    end select
  end function analysis_object_get_n_elements

  function analysis_object_get_n_entries (obj, within_bounds) result (n)
    integer :: n
    type(analysis_object_t), intent(in) :: obj
    logical, intent(in), optional :: within_bounds
    logical :: wb
    select case (obj%type)
    case (AN_OBSERVABLE)
       n = observable_get_n_entries (obj%obs)
    case (AN_HISTOGRAM)
       wb = .false.;  if (present (within_bounds)) wb = within_bounds
       if (wb) then
          n = histogram_get_n_entries_within_bounds (obj%h)
       else
          n = histogram_get_n_entries (obj%h)
       end if
    case (AN_PLOT)
       n = plot_get_n_entries (obj%p)
    case default
       n = 0
    end select
  end function analysis_object_get_n_entries

  function analysis_object_get_average (obj, within_bounds) result (avg)
    real(default) :: avg
    type(analysis_object_t), intent(in) :: obj
    logical, intent(in), optional :: within_bounds
    logical :: wb
    select case (obj%type)
    case (AN_OBSERVABLE)
       avg = observable_get_average (obj%obs)
    case (AN_HISTOGRAM)
       wb = .false.;  if (present (within_bounds)) wb = within_bounds
       if (wb) then
          avg = histogram_get_average_within_bounds (obj%h)
       else
          avg = histogram_get_average (obj%h)
       end if
    case default
       avg = 0
    end select
  end function analysis_object_get_average

  function analysis_object_get_error (obj, within_bounds) result (err)
    real(default) :: err
    type(analysis_object_t), intent(in) :: obj
    logical, intent(in), optional :: within_bounds
    logical :: wb
    select case (obj%type)
    case (AN_OBSERVABLE)
       err = observable_get_error (obj%obs)
    case (AN_HISTOGRAM)
       wb = .false.;  if (present (within_bounds)) wb = within_bounds
       if (wb) then
          err = histogram_get_error_within_bounds (obj%h)
       else
          err = histogram_get_error (obj%h)
       end if
    case default
       err = 0
    end select
  end function analysis_object_get_error

  function analysis_object_get_observable_ptr (obj) result (obs)
    type(observable_t), pointer :: obs
    type(analysis_object_t), intent(in) :: obj
    select case (obj%type)
    case (AN_OBSERVABLE);  obs => obj%obs
    case default;          obs => null ()
    end select
  end function analysis_object_get_observable_ptr

  function analysis_object_get_histogram_ptr (obj) result (h)
    type(histogram_t), pointer :: h
    type(analysis_object_t), intent(in) :: obj
    select case (obj%type)
    case (AN_HISTOGRAM);  h => obj%h
    case default;         h => null ()
    end select
  end function analysis_object_get_histogram_ptr

  function analysis_object_get_plot_ptr (obj) result (plot)
    type(plot_t), pointer :: plot
    type(analysis_object_t), intent(in) :: obj
    select case (obj%type)
    case (AN_PLOT);  plot => obj%p
    case default;    plot => null ()
    end select
  end function analysis_object_get_plot_ptr

  function analysis_object_get_graph_ptr (obj) result (g)
    type(graph_t), pointer :: g
    type(analysis_object_t), intent(in) :: obj
    select case (obj%type)
    case (AN_GRAPH);  g => obj%g
    case default;     g => null ()
    end select
  end function analysis_object_get_graph_ptr

  function analysis_object_has_plot (obj) result (flag)
    logical :: flag
    type(analysis_object_t), intent(in) :: obj
    select case (obj%type)
    case (AN_HISTOGRAM);  flag = .true.
    case (AN_PLOT);       flag = .true.
    case (AN_GRAPH);      flag = .true.
    case default;         flag = .false.
    end select
  end function analysis_object_has_plot

  subroutine analysis_object_write (obj, unit, verbose)
    type(analysis_object_t), intent(in) :: obj
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    logical :: verb
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    verb = .false.;  if (present (verbose))  verb = verbose
    write (u, "(A)")  repeat ("#", 79)
    select case (obj%type)
    case (AN_OBSERVABLE)
       write (u, "(A)", advance="no")  "# Observable:"
    case (AN_HISTOGRAM)
       write (u, "(A)", advance="no")  "# Histogram: "
    case (AN_PLOT)
       write (u, "(A)", advance="no")  "# Plot: "
    case (AN_GRAPH)
       write (u, "(A)", advance="no")  "# Graph: "
    case default
       write (u, "(A)") "# [undefined analysis object]"
       return
    end select
    write (u, "(1x,A)")  char (obj%id)
    select case (obj%type)
    case (AN_OBSERVABLE)
       call observable_write (obj%obs, unit)
    case (AN_HISTOGRAM)
       if (verb) then
          call obj%h%graph_options%write (unit)
          write (u, *)
          call obj%h%drawing_options%write (unit)
          write (u, *)
       end if
       call histogram_write (obj%h, unit)
    case (AN_PLOT)
       if (verb) then
          call obj%p%graph_options%write (unit)
          write (u, *)
          call obj%p%drawing_options%write (unit)
          write (u, *)
       end if
       call plot_write (obj%p, unit)
    case (AN_GRAPH)
       call graph_write (obj%g, unit)
    end select
  end subroutine analysis_object_write

  subroutine analysis_object_write_driver (obj, filename, unit)
    type(analysis_object_t), intent(in) :: obj
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: unit
    select case (obj%type)
    case (AN_OBSERVABLE)
       call observable_write_driver (obj%obs, unit)
    case (AN_HISTOGRAM)
       call histogram_write_gml_driver (obj%h, filename, unit)
    case (AN_PLOT)
       call plot_write_gml_driver (obj%p, filename, unit)
    case (AN_GRAPH)
       call graph_write_gml_driver (obj%g, filename, unit)
    end select
  end subroutine analysis_object_write_driver

  subroutine analysis_object_get_header (obj, header, comment)
    type(analysis_object_t), intent(in) :: obj
    type(ifile_t), intent(inout) :: header
    type(string_t), intent(in), optional :: comment
    select case (obj%type)
    case (AN_HISTOGRAM)
       call histogram_get_header (obj%h, header, comment)
    case (AN_PLOT)
       call plot_get_header (obj%p, header, comment)
    end select
  end subroutine analysis_object_get_header

  subroutine analysis_iterator_init (iterator, object)
    type(analysis_iterator_t), intent(out) :: iterator
    type(analysis_object_t), intent(in), target :: object
    iterator%object => object
    if (associated (iterator%object)) then
       iterator%type = iterator%object%type
       select case (iterator%type)
       case (AN_PLOT)
          iterator%point => iterator%object%p%first
       end select
    end if
  end subroutine analysis_iterator_init

  function analysis_iterator_is_valid (iterator) result (valid)
    logical :: valid
    type(analysis_iterator_t), intent(in) :: iterator
    if (associated (iterator%object)) then
       select case (iterator%type)
       case (AN_HISTOGRAM)
          valid = iterator%index <= histogram_get_n_bins (iterator%object%h)
       case (AN_PLOT)
          valid = associated (iterator%point)
       case default
          valid = .false.
       end select
    else
       valid = .false.
    end if
  end function analysis_iterator_is_valid

  subroutine analysis_iterator_advance (iterator)
    type(analysis_iterator_t), intent(inout) :: iterator
    if (associated (iterator%object)) then
       select case (iterator%type)
       case (AN_PLOT)
          iterator%point => iterator%point%next
       end select
       iterator%index = iterator%index + 1
    end if
  end subroutine analysis_iterator_advance

  function analysis_iterator_get_type (iterator) result (type)
    integer :: type
    type(analysis_iterator_t), intent(in) :: iterator
    type = iterator%type
  end function analysis_iterator_get_type

  subroutine analysis_iterator_get_data (iterator, &
       x, y, yerr, xerr, width, excess, index, n_total)
    type(analysis_iterator_t), intent(in) :: iterator
    real(default), intent(out), optional :: x, y, yerr, xerr, width, excess
    integer, intent(out), optional :: index, n_total
    select case (iterator%type)
    case (AN_HISTOGRAM)
       if (present (x)) &
            x = bin_get_midpoint (iterator%object%h%bin(iterator%index))
       if (present (y)) &
            y = bin_get_sum (iterator%object%h%bin(iterator%index))
       if (present (yerr)) &
            yerr = bin_get_error (iterator%object%h%bin(iterator%index))
       if (present (xerr)) &
            call invalid ("histogram", "xerr")
       if (present (width)) &
            width = bin_get_width (iterator%object%h%bin(iterator%index))
       if (present (excess)) &
            excess = bin_get_excess (iterator%object%h%bin(iterator%index))
       if (present (index)) &
            index = iterator%index
       if (present (n_total)) &
            n_total = histogram_get_n_bins (iterator%object%h)
    case (AN_PLOT)
       if (present (x)) &
            x = point_get_x (iterator%point)
       if (present (y)) &
            y = point_get_y (iterator%point)
       if (present (yerr)) &
            yerr = point_get_yerr (iterator%point)
       if (present (xerr)) &
            xerr = point_get_xerr (iterator%point)
       if (present (width)) &
            call invalid ("plot", "width")
       if (present (excess)) &
            call invalid ("plot", "excess")
       if (present (index)) &
            index = iterator%index
       if (present (n_total)) &
            n_total = plot_get_n_entries (iterator%object%p)
    case default
       call msg_bug ("analysis_iterator_get_data: called " &
            // "for unsupported analysis object type")
    end select
  contains
    subroutine invalid (typestr, objstr)
      character(*), intent(in) :: typestr, objstr
      call msg_bug ("analysis_iterator_get_data: attempt to get '" &
           // objstr // "' for type '" // typestr // "'")
    end subroutine invalid
  end subroutine analysis_iterator_get_data

  module subroutine analysis_final ()
    type(analysis_object_t), pointer :: current
    do while (associated (analysis_store%first))
       current => analysis_store%first
       analysis_store%first => current%next
       call analysis_object_final (current)
    end do
    analysis_store%last => null ()
  end subroutine analysis_final

  subroutine analysis_store_append_object (id, type)
    type(string_t), intent(in) :: id
    integer, intent(in) :: type
    type(analysis_object_t), pointer :: obj
    allocate (obj)
    call analysis_object_init (obj, id, type)
    if (associated (analysis_store%last)) then
       analysis_store%last%next => obj
    else
       analysis_store%first => obj
    end if
    analysis_store%last => obj
  end subroutine analysis_store_append_object

  function analysis_store_get_object_ptr (id) result (obj)
    type(string_t), intent(in) :: id
    type(analysis_object_t), pointer :: obj
    obj => analysis_store%first
    do while (associated (obj))
       if (obj%id == id)  return
       obj => obj%next
    end do
  end function analysis_store_get_object_ptr

  subroutine analysis_store_init_object (id, type, obj)
    type(string_t), intent(in) :: id
    integer, intent(in) :: type
    type(analysis_object_t), pointer :: obj, next
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       next => analysis_object_get_next_ptr (obj)
       call analysis_object_final (obj)
       call analysis_object_init (obj, id, type)
       call analysis_object_set_next_ptr (obj, next)
    else
       call analysis_store_append_object (id, type)
       obj => analysis_store%last
    end if
  end subroutine analysis_store_init_object

  module function analysis_store_get_object_type (id) result (type)
    type(string_t), intent(in) :: id
    integer :: type
    type(analysis_object_t), pointer :: object
    object => analysis_store_get_object_ptr (id)
    if (associated (object)) then
       type = object%type
    else
       type = AN_UNDEFINED
    end if
  end function analysis_store_get_object_type

  function analysis_store_get_n_objects () result (n)
    integer :: n
    type(analysis_object_t), pointer :: current
    n = 0
    current => analysis_store%first
    do while (associated (current))
       n = n + 1
       current => current%next
    end do
  end function analysis_store_get_n_objects

  module subroutine analysis_store_get_ids (id)
    type(string_t), dimension(:), allocatable, intent(out) :: id
    type(analysis_object_t), pointer :: current
    integer :: i
    allocate (id (analysis_store_get_n_objects()))
    i = 0
    current => analysis_store%first
    do while (associated (current))
       i = i + 1
       id(i) = current%id
       current => current%next
    end do
  end subroutine analysis_store_get_ids

  subroutine analysis_store_write_driver_all (filename_data, unit)
    type(string_t), intent(in) :: filename_data
    integer, intent(in), optional :: unit
    type(analysis_object_t), pointer :: obj
    call analysis_store_write_driver_header (unit)
    obj => analysis_store%first
    do while (associated (obj))
       call analysis_object_write_driver (obj, filename_data, unit)
       obj => obj%next
    end do
    call analysis_store_write_driver_footer (unit)
  end subroutine analysis_store_write_driver_all

  subroutine analysis_store_write_driver_obj (filename_data, id, unit)
    type(string_t), intent(in) :: filename_data
    type(string_t), dimension(:), intent(in) :: id
    integer, intent(in), optional :: unit
    type(analysis_object_t), pointer :: obj
    integer :: i
    call analysis_store_write_driver_header (unit)
    do i = 1, size (id)
       obj => analysis_store_get_object_ptr (id(i))
       if (associated (obj))  &
            call analysis_object_write_driver (obj, filename_data, unit)
    end do
    call analysis_store_write_driver_footer (unit)
  end subroutine analysis_store_write_driver_obj

  subroutine analysis_store_write_driver_header (unit)
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, '(A)') "\documentclass[12pt]{article}"
    write (u, *)
    write (u, '(A)') "\usepackage{gamelan}"
    write (u, '(A)') "\usepackage{amsmath}"
    write (u, '(A)') "\usepackage{ifpdf}"
    write (u, '(A)') "\ifpdf"
    write (u, '(A)') "   \DeclareGraphicsRule{*}{mps}{*}{}"
    write (u, '(A)') "\else"
    write (u, '(A)') "   \DeclareGraphicsRule{*}{eps}{*}{}"
    write (u, '(A)') "\fi"
    write (u, *)
    write (u, '(A)') "\begin{document}"
    write (u, '(A)') "\begin{gmlfile}"
    write (u, *)
    write (u, '(A)') "\begin{gmlcode}"
    write (u, '(A)') "  color col.default, col.excess;"
    write (u, '(A)') "  col.default = 0.9white;"
    write (u, '(A)') "  col.excess  = red;"
    write (u, '(A)') "  boolean show_excess;"
    !!! Future excess options for plots
    ! if (mcs(1)%plot_excess .and. mcs(1)%unweighted) then
    !    write (u, '(A)') "  show_excess = true;"
    ! else
    write (u, '(A)') "  show_excess = false;"
    ! end if
    write (u, '(A)') "\end{gmlcode}"
    write (u, *)
  end subroutine analysis_store_write_driver_header

  subroutine analysis_store_write_driver_footer (unit)
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write(u, *)
    write(u, '(A)') "\end{gmlfile}"
    write(u, '(A)') "\end{document}"
  end subroutine analysis_store_write_driver_footer

  module subroutine analysis_init_observable (id, obs_label, obs_unit, graph_options)
    type(string_t), intent(in) :: id
    type(string_t), intent(in), optional :: obs_label, obs_unit
    type(graph_options_t), intent(in), optional :: graph_options
    type(analysis_object_t), pointer :: obj
    type(observable_t), pointer :: obs
    call analysis_store_init_object (id, AN_OBSERVABLE, obj)
    obs => analysis_object_get_observable_ptr (obj)
    call observable_init (obs, obs_label, obs_unit, graph_options)
  end subroutine analysis_init_observable

  module subroutine analysis_init_histogram_n_bins &
       (id, lower_bound, upper_bound, n_bins, normalize_bins, &
        obs_label, obs_unit, graph_options, drawing_options)
    type(string_t), intent(in) :: id
    real(default), intent(in) :: lower_bound, upper_bound
    integer, intent(in) :: n_bins
    logical, intent(in) :: normalize_bins
    type(string_t), intent(in), optional :: obs_label, obs_unit
    type(graph_options_t), intent(in), optional :: graph_options
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(analysis_object_t), pointer :: obj
    type(histogram_t), pointer :: h
    call analysis_store_init_object (id, AN_HISTOGRAM, obj)
    h => analysis_object_get_histogram_ptr (obj)
    call histogram_init (h, id, &
         lower_bound, upper_bound, n_bins, normalize_bins, &
         obs_label, obs_unit, graph_options, drawing_options)
  end subroutine analysis_init_histogram_n_bins

  module subroutine analysis_init_histogram_bin_width &
       (id, lower_bound, upper_bound, bin_width, normalize_bins, &
        obs_label, obs_unit, graph_options, drawing_options)
    type(string_t), intent(in) :: id
    real(default), intent(in) :: lower_bound, upper_bound, bin_width
    logical, intent(in) :: normalize_bins
    type(string_t), intent(in), optional :: obs_label, obs_unit
    type(graph_options_t), intent(in), optional :: graph_options
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(analysis_object_t), pointer :: obj
    type(histogram_t), pointer :: h
    call analysis_store_init_object (id, AN_HISTOGRAM, obj)
    h => analysis_object_get_histogram_ptr (obj)
    call histogram_init (h, id, &
         lower_bound, upper_bound, bin_width, normalize_bins, &
         obs_label, obs_unit, graph_options, drawing_options)
  end subroutine analysis_init_histogram_bin_width

  module subroutine analysis_init_plot (id, graph_options, drawing_options)
    type(string_t), intent(in) :: id
    type(graph_options_t), intent(in), optional :: graph_options
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(analysis_object_t), pointer :: obj
    type(plot_t), pointer :: plot
    call analysis_store_init_object (id, AN_PLOT, obj)
    plot => analysis_object_get_plot_ptr (obj)
    call plot_init (plot, id, graph_options, drawing_options)
  end subroutine analysis_init_plot

  module subroutine analysis_init_graph (id, n_elements, graph_options)
    type(string_t), intent(in) :: id
    integer, intent(in) :: n_elements
    type(graph_options_t), intent(in), optional :: graph_options
    type(analysis_object_t), pointer :: obj
    type(graph_t), pointer :: graph
    call analysis_store_init_object (id, AN_GRAPH, obj)
    graph => analysis_object_get_graph_ptr (obj)
    call graph_init (graph, id, n_elements, graph_options)
  end subroutine analysis_init_graph

  module subroutine analysis_store_clear_obj (id)
    type(string_t), intent(in) :: id
    type(analysis_object_t), pointer :: obj
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       call analysis_object_clear (obj)
    end if
  end subroutine analysis_store_clear_obj

  module subroutine analysis_store_clear_all ()
    type(analysis_object_t), pointer :: obj
    obj => analysis_store%first
    do while (associated (obj))
       call analysis_object_clear (obj)
       obj => obj%next
    end do
  end subroutine analysis_store_clear_all

  module subroutine analysis_record_data (id, x, y, yerr, xerr, &
       weight, excess, success, exist)
    type(string_t), intent(in) :: id
    real(default), intent(in) :: x
    real(default), intent(in), optional :: y, yerr, xerr, weight, excess
    logical, intent(out), optional :: success, exist
    type(analysis_object_t), pointer :: obj
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       call analysis_object_record_data (obj, x, y, yerr, xerr, &
            weight, excess, success)
       if (present (exist))  exist = .true.
    else
       if (present (success))  success = .false.
       if (present (exist))  exist = .false.
    end if
  end subroutine analysis_record_data

  module subroutine analysis_fill_graph (id, i, id_in, drawing_options)
    type(string_t), intent(in) :: id
    integer, intent(in) :: i
    type(string_t), intent(in) :: id_in
    type(drawing_options_t), intent(in), optional :: drawing_options
    type(analysis_object_t), pointer :: obj
    type(graph_t), pointer :: g
    type(histogram_t), pointer :: h
    type(plot_t), pointer :: p
    obj => analysis_store_get_object_ptr (id)
    g => analysis_object_get_graph_ptr (obj)
    obj => analysis_store_get_object_ptr (id_in)
    if (associated (obj)) then
       select case (obj%type)
       case (AN_HISTOGRAM)
          h => analysis_object_get_histogram_ptr (obj)
          call graph_insert_histogram (g, i, h, drawing_options)
       case (AN_PLOT)
          p => analysis_object_get_plot_ptr (obj)
          call graph_insert_plot (g, i, p, drawing_options)
       case default
          call msg_error ("Graph '" // char (id) // "': Element '" &
               // char (id_in) // "' is neither histogram nor plot.")
       end select
    else
       call msg_error ("Graph '" // char (id) // "': Element '" &
               // char (id_in) // "' is undefined.")
    end if
  end subroutine analysis_fill_graph

  module function analysis_exists (id) result (flag)
    type(string_t), intent(in) :: id
    logical :: flag
    type(analysis_object_t), pointer :: obj
    flag = .true.
    obj => analysis_store%first
    do while (associated (obj))
       if (obj%id == id)  return
       obj => obj%next
    end do
    flag = .false.
  end function analysis_exists

  module function analysis_get_n_elements (id) result (n)
    integer :: n
    type(string_t), intent(in) :: id
    type(analysis_object_t), pointer :: obj
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       n = analysis_object_get_n_elements (obj)
    else
       n = 0
    end if
  end function analysis_get_n_elements

  module function analysis_get_n_entries (id, within_bounds) result (n)
    integer :: n
    type(string_t), intent(in) :: id
    logical, intent(in), optional :: within_bounds
    type(analysis_object_t), pointer :: obj
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       n = analysis_object_get_n_entries (obj, within_bounds)
    else
       n = 0
    end if
  end function analysis_get_n_entries

  module function analysis_get_average (id, within_bounds) result (avg)
    real(default) :: avg
    type(string_t), intent(in) :: id
    type(analysis_object_t), pointer :: obj
    logical, intent(in), optional :: within_bounds
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       avg = analysis_object_get_average (obj, within_bounds)
    else
       avg = 0
    end if
  end function analysis_get_average

  module function analysis_get_error (id, within_bounds) result (err)
    real(default) :: err
    type(string_t), intent(in) :: id
    type(analysis_object_t), pointer :: obj
    logical, intent(in), optional :: within_bounds
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       err = analysis_object_get_error (obj, within_bounds)
    else
       err = 0
    end if
  end function analysis_get_error

  module function analysis_has_plots_any () result (flag)
    logical :: flag
    type(analysis_object_t), pointer :: obj
    flag = .false.
    obj => analysis_store%first
    do while (associated (obj))
       flag = analysis_object_has_plot (obj)
       if (flag)  return
    end do
  end function analysis_has_plots_any

  module function analysis_has_plots_obj (id) result (flag)
    logical :: flag
    type(string_t), dimension(:), intent(in) :: id
    type(analysis_object_t), pointer :: obj
    integer :: i
    flag = .false.
    do i = 1, size (id)
       obj => analysis_store_get_object_ptr (id(i))
       if (associated (obj)) then
          flag = analysis_object_has_plot (obj)
          if (flag)  return
       end if
    end do
  end function analysis_has_plots_obj

  subroutine analysis_init_iterator (id, iterator)
    type(string_t), intent(in) :: id
    type(analysis_iterator_t), intent(out) :: iterator
    type(analysis_object_t), pointer :: obj
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj))  call analysis_iterator_init (iterator, obj)
  end subroutine analysis_init_iterator

  module subroutine analysis_write_object (id, unit, verbose)
    type(string_t), intent(in) :: id
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    type(analysis_object_t), pointer :: obj
    obj => analysis_store_get_object_ptr (id)
    if (associated (obj)) then
       call analysis_object_write (obj, unit, verbose)
    else
       call msg_error ("Analysis object '" // char (id) // "' not found")
    end if
  end subroutine analysis_write_object

  module subroutine analysis_write_all (unit, verbose)
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    type(analysis_object_t), pointer :: obj
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    obj => analysis_store%first
    do while (associated (obj))
       call analysis_object_write (obj, unit, verbose)
       obj => obj%next
    end do
  end subroutine analysis_write_all

  module subroutine analysis_write_driver (filename_data, id, unit)
    type(string_t), intent(in) :: filename_data
    type(string_t), dimension(:), intent(in), optional :: id
    integer, intent(in), optional :: unit
    if (present (id)) then
       call analysis_store_write_driver_obj (filename_data, id, unit)
    else
       call analysis_store_write_driver_all (filename_data, unit)
    end if
  end subroutine analysis_write_driver

  module subroutine analysis_compile_tex (file, has_gmlcode, os_data)
    type(string_t), intent(in) :: file
    logical, intent(in) :: has_gmlcode
    type(os_data_t), intent(in) :: os_data
    integer :: status
    if (os_data%event_analysis_ps) then
       call os_system_call ("make compile " // os_data%makeflags // " -f " // &
            char (file) // "_ana.makefile", status)
       if (status /= 0) then
          call msg_error ("Unable to compile analysis output file")
       end if
    else
       call msg_warning ("Skipping results display because " &
            // "latex/mpost/dvips is not available")
    end if
  end subroutine analysis_compile_tex

  subroutine analysis_get_header (id, header, comment)
    type(string_t), intent(in) :: id
    type(ifile_t), intent(inout) :: header
    type(string_t), intent(in), optional :: comment
    type(analysis_object_t), pointer :: object
    object => analysis_store_get_object_ptr (id)
    if (associated (object)) then
       call analysis_object_get_header (object, header, comment)
    end if
  end subroutine analysis_get_header

  module subroutine analysis_write_makefile (filename, unit, has_gmlcode, os_data)
    type(string_t), intent(in) :: filename
    integer, intent(in) :: unit
    logical, intent(in) :: has_gmlcode
    type(os_data_t), intent(in) :: os_data
    write (unit, "(3A)")  "# WHIZARD: Makefile for analysis '", &
         char (filename), "'"
    write (unit, "(A)")  "# Automatically generated file, do not edit"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "# LaTeX setup"
    write (unit, "(A)")  "LATEX = " // char (os_data%latex)
    write (unit, "(A)")  "MPOST = " // char (os_data%mpost)
    write (unit, "(A)")  "GML = " // char (os_data%gml)
    write (unit, "(A)")  "DVIPS = " // char (os_data%dvips)
    write (unit, "(A)")  "PS2PDF = " // char (os_data%ps2pdf)
    write (unit, "(A)")  'TEX_FLAGS = "$$TEXINPUTS:' // &
         char(os_data%whizard_texpath) // '"'
    write (unit, "(A)")  'MP_FLAGS  = "$$MPINPUTS:' // &
         char(os_data%whizard_texpath) // '"'
    write (unit, "(A)")  ""
    write (unit, "(5A)")  "TEX_SOURCES = ", char (filename), ".tex"
    if (os_data%event_analysis_pdf) then
       write (unit, "(5A)")  "TEX_OBJECTS = ", char (filename), ".pdf"
    else
       write (unit, "(5A)")  "TEX_OBJECTS = ", char (filename), ".ps"
    end if
    if (os_data%event_analysis_ps) then
       if (os_data%event_analysis_pdf) then
          write (unit, "(5A)")  char (filename), ".pdf: ", &
               char (filename), ".tex"
       else
          write (unit, "(5A)")  char (filename), ".ps: ", &
               char (filename), ".tex"
       end if
       write (unit, "(5A)")  TAB, "-TEXINPUTS=$(TEX_FLAGS) $(LATEX) " // &
            char (filename) // ".tex"
       if (has_gmlcode) then
          write (unit, "(5A)")  TAB, "$(GML) " // char (filename)
          write (unit, "(5A)")  TAB, "TEXINPUTS=$(TEX_FLAGS) $(LATEX) " // &
            char (filename) // ".tex"
       end if
       write (unit, "(5A)")  TAB, "$(DVIPS) -o " // char (filename) // ".ps " // &
            char (filename) // ".dvi"
       if (os_data%event_analysis_pdf) then
          write (unit, "(5A)")  TAB, "$(PS2PDF) " // char (filename) // ".ps"
       end if
    end if
    write (unit, "(A)")
    write (unit, "(A)")  "compile: $(TEX_OBJECTS)"
    write (unit, "(A)")  ".PHONY: compile"
    write (unit, "(A)")
    write (unit, "(5A)")  "CLEAN_OBJECTS = ",  char (filename), ".aux"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".log"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".dvi"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".out"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".[1-9]"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".[1-9][0-9]"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".[1-9][0-9][0-9]"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".t[1-9]"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".t[1-9][0-9]"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".t[1-9][0-9][0-9]"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".ltp"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".mp"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".mpx"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".dvi"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".ps"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (filename), ".pdf"
    write (unit, "(A)")
    write (unit, "(A)")  "# Generic cleanup targets"
    write (unit, "(A)")  "clean-objects:"
    write (unit, "(A)")  TAB // "rm -f $(CLEAN_OBJECTS)"
    write (unit, "(A)")  ""
    write (unit, "(A)")  "clean: clean-objects"
    write (unit, "(A)")  ".PHONY: clean"
  end subroutine analysis_write_makefile


end submodule analysis_s

