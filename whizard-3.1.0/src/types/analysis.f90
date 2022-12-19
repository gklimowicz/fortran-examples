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

module analysis

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface

  implicit none
  private

  public :: graph_options_t
  public :: drawing_options_t
  public :: AN_UNDEFINED, AN_HISTOGRAM, AN_OBSERVABLE, AN_PLOT, AN_GRAPH

  public :: analysis_final
  public :: analysis_store_get_object_type
  public :: analysis_store_get_ids
  public :: analysis_init_observable
  public :: analysis_init_histogram
  public :: analysis_init_plot
  public :: analysis_init_graph
  public :: analysis_clear
  public :: analysis_record_data
  public :: analysis_fill_graph
  public :: analysis_exists
  public :: analysis_get_n_elements
  public :: analysis_get_n_entries
  public :: analysis_get_average
  public :: analysis_get_error

  public :: analysis_has_plots
  public :: analysis_write
  public :: analysis_write_driver
  public :: analysis_compile_tex
  public :: analysis_write_makefile

  character(*), parameter, public :: HISTOGRAM_HEAD_FORMAT = "1x,A15,3x"
  character(*), parameter, public :: HISTOGRAM_INTG_FORMAT = "3x,I9,3x"
  character(*), parameter, public :: HISTOGRAM_DATA_FORMAT = "ES19.12"

  integer, parameter :: AN_UNDEFINED = 0
  integer, parameter :: AN_OBSERVABLE = 1
  integer, parameter :: AN_HISTOGRAM = 2
  integer, parameter :: AN_PLOT = 3
  integer, parameter :: AN_GRAPH = 4

  type :: graph_options_t
     private
     type(string_t) :: id
     type(string_t) :: title
     type(string_t) :: description
     type(string_t) :: x_label
     type(string_t) :: y_label
     integer :: width_mm = 130
     integer :: height_mm = 90
     logical :: x_log = .false.
     logical :: y_log = .false.
     real(default) :: x_min = 0
     real(default) :: x_max = 1
     real(default) :: y_min = 0
     real(default) :: y_max = 1
     logical :: x_min_set = .false.
     logical :: x_max_set = .false.
     logical :: y_min_set = .false.
     logical :: y_max_set = .false.
     type(string_t) :: gmlcode_bg
     type(string_t) :: gmlcode_fg
   contains
     procedure :: init => graph_options_init
     procedure :: set => graph_options_set
     procedure :: write => graph_options_write
  end type graph_options_t

  type :: drawing_options_t
     type(string_t) :: dataset
     logical :: with_hbars = .false.
     logical :: with_base = .false.
     logical :: piecewise = .false.
     logical :: fill = .false.
     logical :: draw = .false.
     logical :: err = .false.
     logical :: symbols = .false.
     type(string_t) :: fill_options
     type(string_t) :: draw_options
     type(string_t) :: err_options
     type(string_t) :: symbol
     type(string_t) :: gmlcode_bg
     type(string_t) :: gmlcode_fg
   contains
     procedure :: write => drawing_options_write
     procedure :: init_histogram => drawing_options_init_histogram
     procedure :: init_plot => drawing_options_init_plot
     procedure :: set => drawing_options_set
  end type drawing_options_t

  type :: observable_t
     private
     real(default) :: sum_values = 0
     real(default) :: sum_squared_values = 0
     real(default) :: sum_weights = 0
     real(default) :: sum_squared_weights = 0
     integer :: count = 0
     type(string_t) :: obs_label
     type(string_t) :: obs_unit
     type(graph_options_t) :: graph_options
  end type observable_t

  type :: bin_t
     private
     real(default) :: midpoint = 0
     real(default) :: width = 0
     real(default) :: sum_weights = 0
     real(default) :: sum_squared_weights = 0
     real(default) :: sum_excess_weights = 0
     integer :: count = 0
  end type bin_t

  type :: histogram_t
     private
     real(default) :: lower_bound = 0
     real(default) :: upper_bound = 0
     real(default) :: width = 0
     integer :: n_bins = 0
     logical :: normalize_bins = .false.
     type(observable_t) :: obs
     type(observable_t) :: obs_within_bounds
     type(bin_t) :: underflow
     type(bin_t), dimension(:), allocatable :: bin
     type(bin_t) :: overflow
     type(graph_options_t) :: graph_options
     type(drawing_options_t) :: drawing_options
  end type histogram_t

  type :: point_t
     private
     real(default) :: x = 0
     real(default) :: y = 0
     real(default) :: yerr = 0
     real(default) :: xerr = 0
     type(point_t), pointer :: next => null ()
  end type point_t

  type :: plot_t
     private
     type(point_t), pointer :: first => null ()
     type(point_t), pointer :: last => null ()
     integer :: count = 0
     type(graph_options_t) :: graph_options
     type(drawing_options_t) :: drawing_options
  end type plot_t

  type :: graph_element_t
     private
     integer :: type = AN_UNDEFINED
     type(histogram_t), pointer :: h => null ()
     type(plot_t), pointer :: p => null ()
  end type graph_element_t

  type :: graph_t
     private
     type(graph_element_t), dimension(:), allocatable :: el
     type(graph_options_t) :: graph_options
  end type graph_t

  type :: analysis_object_t
     private
     type(string_t) :: id
     integer :: type = AN_UNDEFINED
     type(observable_t), pointer :: obs => null ()
     type(histogram_t), pointer :: h => null ()
     type(plot_t), pointer :: p => null ()
     type(graph_t), pointer :: g => null ()
     type(analysis_object_t), pointer :: next => null ()
  end type analysis_object_t

  type :: analysis_iterator_t
    private
    integer :: type = AN_UNDEFINED
    type(analysis_object_t), pointer :: object => null ()
    integer :: index = 1
    type(point_t), pointer :: point => null ()
  end type

  type :: analysis_store_t
     private
     type(analysis_object_t), pointer :: first => null ()
     type(analysis_object_t), pointer :: last => null ()
  end type analysis_store_t


  interface observable_record_value
     module procedure observable_record_value_unweighted
     module procedure observable_record_value_weighted
  end interface

  interface histogram_init
     module procedure histogram_init_n_bins
     module procedure histogram_init_bin_width
  end interface

  interface point_init
    module procedure point_init_contents
    module procedure point_init_point
  end interface
  interface plot_init
     module procedure plot_init_empty
     module procedure plot_init_plot
  end interface
  interface analysis_init_histogram
     module procedure analysis_init_histogram_n_bins
     module procedure analysis_init_histogram_bin_width
  end interface

  interface analysis_clear
     module procedure analysis_store_clear_obj
     module procedure analysis_store_clear_all
  end interface

  interface analysis_has_plots
     module procedure analysis_has_plots_any
     module procedure analysis_has_plots_obj
  end interface

  interface analysis_write
     module procedure analysis_write_object
     module procedure analysis_write_all
  end interface


  type(analysis_store_t), save :: analysis_store


  interface
    module subroutine graph_options_init (graph_options)
      class(graph_options_t), intent(out) :: graph_options
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
    end subroutine graph_options_set
    module subroutine graph_options_write (gro, unit)
      class(graph_options_t), intent(in) :: gro
      integer, intent(in), optional :: unit
    end subroutine graph_options_write
    module subroutine drawing_options_write (dro, unit)
      class(drawing_options_t), intent(in) :: dro
      integer, intent(in), optional :: unit
    end subroutine drawing_options_write
    module subroutine drawing_options_init_histogram (dro)
      class(drawing_options_t), intent(out) :: dro
    end subroutine drawing_options_init_histogram
    module subroutine drawing_options_init_plot (dro)
      class(drawing_options_t), intent(out) :: dro
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
    end subroutine drawing_options_set
    module subroutine observable_record_value_unweighted (obs, value, success)
      type(observable_t), intent(inout) :: obs
      real(default), intent(in) :: value
      logical, intent(out), optional :: success
    end subroutine observable_record_value_unweighted
    module subroutine observable_record_value_weighted (obs, value, weight, success)
      type(observable_t), intent(inout) :: obs
      real(default), intent(in) :: value, weight
      logical, intent(out), optional :: success
    end subroutine observable_record_value_weighted
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
    end subroutine histogram_init_bin_width
    module subroutine point_init_contents (point, x, y, yerr, xerr)
      type(point_t), intent(out) :: point
      real(default), intent(in) :: x, y
      real(default), intent(in), optional :: yerr, xerr
    end subroutine point_init_contents
    module subroutine point_init_point (point, point_in)
      type(point_t), intent(out) :: point
      type(point_t), intent(in) :: point_in
    end subroutine point_init_point
    module subroutine plot_init_empty (p, id, graph_options, drawing_options)
      type(plot_t), intent(out) :: p
      type(string_t), intent(in) :: id
      type(graph_options_t), intent(in), optional :: graph_options
      type(drawing_options_t), intent(in), optional :: drawing_options
    end subroutine plot_init_empty
    module subroutine plot_init_plot (p, p_in, drawing_options)
      type(plot_t), intent(out) :: p
      type(plot_t), intent(in) :: p_in
      type(drawing_options_t), intent(in), optional :: drawing_options
    end subroutine plot_init_plot
    module subroutine analysis_final ()
    end subroutine analysis_final
    module function analysis_store_get_object_type (id) result (type)
      type(string_t), intent(in) :: id
      integer :: type
    end function analysis_store_get_object_type
    module subroutine analysis_store_get_ids (id)
      type(string_t), dimension(:), allocatable, intent(out) :: id
    end subroutine analysis_store_get_ids
    module subroutine analysis_init_observable (id, obs_label, obs_unit, graph_options)
      type(string_t), intent(in) :: id
      type(string_t), intent(in), optional :: obs_label, obs_unit
      type(graph_options_t), intent(in), optional :: graph_options
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
    end subroutine analysis_init_histogram_bin_width
    module subroutine analysis_init_plot (id, graph_options, drawing_options)
      type(string_t), intent(in) :: id
      type(graph_options_t), intent(in), optional :: graph_options
      type(drawing_options_t), intent(in), optional :: drawing_options
    end subroutine analysis_init_plot
    module subroutine analysis_init_graph (id, n_elements, graph_options)
      type(string_t), intent(in) :: id
      integer, intent(in) :: n_elements
      type(graph_options_t), intent(in), optional :: graph_options
    end subroutine analysis_init_graph
    module subroutine analysis_store_clear_obj (id)
      type(string_t), intent(in) :: id
    end subroutine analysis_store_clear_obj
    module subroutine analysis_store_clear_all ()
    end subroutine analysis_store_clear_all
    module subroutine analysis_record_data (id, x, y, yerr, xerr, &
         weight, excess, success, exist)
      type(string_t), intent(in) :: id
      real(default), intent(in) :: x
      real(default), intent(in), optional :: y, yerr, xerr, weight, excess
      logical, intent(out), optional :: success, exist
    end subroutine analysis_record_data
    module subroutine analysis_fill_graph (id, i, id_in, drawing_options)
      type(string_t), intent(in) :: id
      integer, intent(in) :: i
      type(string_t), intent(in) :: id_in
      type(drawing_options_t), intent(in), optional :: drawing_options
    end subroutine analysis_fill_graph
    module function analysis_exists (id) result (flag)
      type(string_t), intent(in) :: id
      logical :: flag
    end function analysis_exists
    module function analysis_get_n_elements (id) result (n)
      integer :: n
      type(string_t), intent(in) :: id
    end function analysis_get_n_elements
    module function analysis_get_n_entries (id, within_bounds) result (n)
      integer :: n
      type(string_t), intent(in) :: id
      logical, intent(in), optional :: within_bounds
    end function analysis_get_n_entries
    module function analysis_get_average (id, within_bounds) result (avg)
      real(default) :: avg
      type(string_t), intent(in) :: id
      logical, intent(in), optional :: within_bounds
    end function analysis_get_average
    module function analysis_get_error (id, within_bounds) result (err)
      real(default) :: err
      type(string_t), intent(in) :: id
      logical, intent(in), optional :: within_bounds
    end function analysis_get_error
    module function analysis_has_plots_any () result (flag)
      logical :: flag
    end function analysis_has_plots_any
    module function analysis_has_plots_obj (id) result (flag)
      logical :: flag
      type(string_t), dimension(:), intent(in) :: id
    end function analysis_has_plots_obj
    module subroutine analysis_write_object (id, unit, verbose)
      type(string_t), intent(in) :: id
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine analysis_write_object
    module subroutine analysis_write_all (unit, verbose)
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine analysis_write_all
    module subroutine analysis_write_driver (filename_data, id, unit)
      type(string_t), intent(in) :: filename_data
      type(string_t), dimension(:), intent(in), optional :: id
      integer, intent(in), optional :: unit
    end subroutine analysis_write_driver
    module subroutine analysis_compile_tex (file, has_gmlcode, os_data)
      type(string_t), intent(in) :: file
      logical, intent(in) :: has_gmlcode
      type(os_data_t), intent(in) :: os_data
    end subroutine analysis_compile_tex
    module subroutine analysis_write_makefile (filename, unit, has_gmlcode, os_data)
      type(string_t), intent(in) :: filename
      integer, intent(in) :: unit
      logical, intent(in) :: has_gmlcode
      type(os_data_t), intent(in) :: os_data
    end subroutine analysis_write_makefile
  end interface

end module analysis
