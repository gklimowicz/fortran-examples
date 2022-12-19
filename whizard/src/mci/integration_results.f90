module integration_results

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use mci_base

  implicit none
  private

  public :: integration_entry_t
  public :: integration_results_t
  public :: integration_results_write_driver
  public :: integration_results_compile_driver

  integer, parameter, public :: PRC_UNKNOWN = 0
  integer, parameter, public :: PRC_DECAY = 1
  integer, parameter, public :: PRC_SCATTERING = 2

  integer, parameter :: RESULTS_CHUNK_SIZE = 10

  real(default), parameter, public :: INTEGRATION_ERROR_TOLERANCE = 1e-10
  real, parameter, public :: GML_MIN_RANGE_RATIO = 0.02

  type :: integration_entry_t
     private
     integer :: process_type = PRC_UNKNOWN
     integer :: pass = 0
     integer :: it = 0
     integer :: n_it = 0
     integer :: n_calls = 0
     integer :: n_calls_valid = 0
     logical :: improved = .false.
     real(default) :: integral = 0
     real(default) :: error = 0
     real(default) :: efficiency = 0
     real(default) :: efficiency_pos = 0
     real(default) :: efficiency_neg = 0
     real(default) :: chi2 = 0
     real(default), dimension(:), allocatable :: chain_weights
   contains
     procedure :: get_pass => integration_entry_get_pass
     procedure :: get_n_calls => integration_entry_get_n_calls
     procedure :: get_n_calls_valid => integration_entry_get_n_calls_valid
     procedure :: get_integral => integration_entry_get_integral
     procedure :: get_error => integration_entry_get_error
     procedure :: get_rel_error => integration_entry_get_relative_error
     procedure :: get_accuracy => integration_entry_get_accuracy
     procedure :: get_efficiency => integration_entry_get_efficiency
     procedure :: get_efficiency_pos => integration_entry_get_efficiency_pos
     procedure :: get_efficiency_neg => integration_entry_get_efficiency_neg
     procedure :: get_chi2 => integration_entry_get_chi2
     procedure :: has_improved => integration_entry_has_improved
     procedure :: get_n_groves => integration_entry_get_n_groves
     procedure :: write => integration_entry_write
     procedure :: write_verbose => integration_entry_write_verbose
     procedure :: read => integration_entry_read
     procedure :: write_chain_weights => integration_entry_write_chain_weights
  end type integration_entry_t

  type, extends (mci_results_t) :: integration_results_t
     private
     integer :: process_type = PRC_UNKNOWN
     integer :: current_pass = 0
     integer :: n_pass = 0
     integer :: n_it = 0
     logical :: screen = .false.
     integer :: unit = 0
     integer :: verbosity = 0
     real(default) :: error_threshold = 0
     type(integration_entry_t), dimension(:), allocatable :: entry
     type(integration_entry_t), dimension(:), allocatable :: average
   contains
     procedure :: init => integration_results_init
     procedure :: set_verbosity => integration_results_set_verbosity
     procedure :: set_error_threshold => integration_results_set_error_threshold
     procedure :: write => integration_results_write
     procedure :: write_verbose => integration_results_write_verbose
     procedure :: write_chain_weights => &
          integration_results_write_chain_weights
     procedure :: read => integration_results_read
     procedure, private :: write_header
     procedure, private :: write_hline
     procedure, private :: write_dline
     procedure :: display_init => integration_results_display_init
     procedure :: display_current => integration_results_display_current
     procedure :: display_pass => integration_results_display_pass
     procedure :: display_final => integration_results_display_final
     procedure :: expand => integration_results_expand
     procedure :: new_pass => integration_results_new_pass
     procedure :: append => integration_results_append
     procedure :: record_simple => integration_results_record_simple
     procedure :: record_extended => integration_results_record_extended
     procedure :: exist => integration_results_exist
     procedure :: get_entry => results_get_entry
     procedure :: get_n_calls => integration_results_get_n_calls
     procedure :: get_integral => integration_results_get_integral
     procedure :: get_error => integration_results_get_error
     procedure :: get_accuracy => integration_results_get_accuracy
     procedure :: get_chi2 => integration_results_get_chi2
     procedure :: get_efficiency => integration_results_get_efficiency
     procedure :: pacify => integration_results_pacify
     procedure :: record_correction => integration_results_record_correction
  end type integration_results_t


  interface integration_entry_t
     module procedure integration_entry_init
  end interface integration_entry_t


  interface
    module function integration_entry_init (process_type, pass,&
         & it, n_it, n_calls, n_calls_valid, improved, integral, error,&
         & efficiency, efficiency_pos, efficiency_neg, chi2, chain_weights)&
         & result (entry)
      type(integration_entry_t) :: entry
      integer, intent(in) :: process_type, pass, it, n_it, &
           n_calls, n_calls_valid
      logical, intent(in) :: improved
      real(default), intent(in) :: integral, error, efficiency, &
           efficiency_pos, efficiency_neg
      real(default), intent(in), optional :: chi2
      real(default), dimension(:), intent(in), optional :: chain_weights
    end function integration_entry_init
    elemental module function integration_entry_get_pass (entry) result (n)
      integer :: n
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_pass
    elemental module function integration_entry_get_n_calls (entry) result (n)
      integer :: n
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_n_calls
    elemental module function integration_entry_get_n_calls_valid &
         (entry) result (n)
      integer :: n
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_n_calls_valid
    elemental module function integration_entry_get_integral (entry) result (int)
      real(default) :: int
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_integral
    elemental module function integration_entry_get_error (entry) result (err)
      real(default) :: err
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_error
    elemental module function integration_entry_get_relative_error &
         (entry) result (err)
      real(default) :: err
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_relative_error
    elemental module function integration_entry_get_accuracy &
         (entry) result (acc)
      real(default) :: acc
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_accuracy
    elemental module function accuracy (integral, error, n_calls) result (acc)
      real(default) :: acc
      real(default), intent(in) :: integral, error
      integer, intent(in) :: n_calls
    end function accuracy
    elemental module function integration_entry_get_efficiency &
         (entry) result (eff)
      real(default) :: eff
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_efficiency
    elemental module function integration_entry_get_efficiency_pos &
         (entry) result (eff)
      real(default) :: eff
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_efficiency_pos
    elemental module function integration_entry_get_efficiency_neg &
         (entry) result (eff)
      real(default) :: eff
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_efficiency_neg
    elemental module function integration_entry_get_chi2 (entry) result (chi2)
      real(default) :: chi2
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_chi2
    elemental module function integration_entry_has_improved &
         (entry) result (flag)
      logical :: flag
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_has_improved
    elemental module function integration_entry_get_n_groves &
         (entry) result (n_groves)
      integer :: n_groves
      class(integration_entry_t), intent(in) :: entry
    end function integration_entry_get_n_groves
    module subroutine integration_entry_write (entry, unit, verbosity, suppress)
      class(integration_entry_t), intent(in) :: entry
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: verbosity
      logical, intent(in), optional :: suppress
    end subroutine integration_entry_write
    module subroutine integration_entry_write_verbose (entry, unit)
      class(integration_entry_t), intent(in) :: entry
      integer, intent(in) :: unit
    end subroutine integration_entry_write_verbose
    module subroutine integration_entry_read (entry, unit)
      class(integration_entry_t), intent(out) :: entry
      integer, intent(in) :: unit
    end subroutine integration_entry_read
    module subroutine integration_entry_write_chain_weights (entry, unit)
      class(integration_entry_t), intent(in) :: entry
      integer, intent(in), optional :: unit
    end subroutine integration_entry_write_chain_weights
    module subroutine integration_results_init (results, process_type)
      class(integration_results_t), intent(out) :: results
      integer, intent(in) :: process_type
    end subroutine integration_results_init
    module subroutine integration_results_set_verbosity (results, verbosity)
      class(integration_results_t), intent(inout) :: results
      integer, intent(in) :: verbosity
    end subroutine integration_results_set_verbosity
    module subroutine integration_results_set_error_threshold &
         (results, error_threshold)
      class(integration_results_t), intent(inout) :: results
      real(default), intent(in) :: error_threshold
    end subroutine integration_results_set_error_threshold
    module subroutine integration_results_write (object, unit, suppress)
      class(integration_results_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: suppress
    end subroutine integration_results_write
    module subroutine integration_results_write_verbose (object, unit)
      class(integration_results_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine integration_results_write_verbose
    module subroutine integration_results_write_chain_weights (results, unit)
      class(integration_results_t), intent(in) :: results
      integer, intent(in), optional :: unit
    end subroutine integration_results_write_chain_weights
    module subroutine integration_results_read (results, unit)
      class(integration_results_t), intent(out) :: results
      integer, intent(in) :: unit
    end subroutine integration_results_read
    module subroutine write_header (results, unit, logfile)
      class(integration_results_t), intent(in) :: results
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: logfile
    end subroutine write_header
    module subroutine write_hline (results, unit)
      class(integration_results_t), intent(in) :: results
      integer, intent(in), optional :: unit
    end subroutine write_hline
    module subroutine write_dline (results, unit)
      class(integration_results_t), intent(in) :: results
      integer, intent(in), optional :: unit
    end subroutine write_dline
    module subroutine integration_results_display_init &
         (results, screen, unit)
      class(integration_results_t), intent(inout) :: results
      logical, intent(in) :: screen
      integer, intent(in), optional :: unit
    end subroutine integration_results_display_init
    module subroutine integration_results_display_current (results, pacify)
      class(integration_results_t), intent(in) :: results
      logical, intent(in), optional :: pacify
    end subroutine integration_results_display_current
    module subroutine integration_results_display_pass (results, pacify)
      class(integration_results_t), intent(in) :: results
      logical, intent(in), optional :: pacify
    end subroutine integration_results_display_pass
    module subroutine integration_results_display_final (results)
      class(integration_results_t), intent(inout) :: results
    end subroutine integration_results_display_final
    module subroutine integration_results_expand (results)
      class(integration_results_t), intent(inout) :: results
    end subroutine integration_results_expand
    module subroutine integration_results_new_pass (results)
      class(integration_results_t), intent(inout) :: results
    end subroutine integration_results_new_pass
    module subroutine integration_results_append (results, &
         n_it, n_calls, n_calls_valid, &
         integral, error, efficiency, efficiency_pos, efficiency_neg, &
         chain_weights)
      class(integration_results_t), intent(inout) :: results
      integer, intent(in) :: n_it, n_calls, n_calls_valid
      real(default), intent(in) :: integral, error, efficiency, &
           efficiency_pos, efficiency_neg
      real(default), dimension(:), intent(in), optional :: chain_weights
    end subroutine integration_results_append
    module subroutine integration_results_record_simple &
         (object, n_it, n_calls, integral, error, efficiency, &
          chain_weights, suppress)
      class(integration_results_t), intent(inout) :: object
      integer, intent(in) :: n_it, n_calls
      real(default), intent(in) :: integral, error, efficiency
      real(default), dimension(:), intent(in), optional :: chain_weights
      logical, intent(in), optional :: suppress
    end subroutine integration_results_record_simple
    module subroutine integration_results_record_extended (object, n_it, &
         n_calls, n_calls_valid, integral, error, efficiency, efficiency_pos, &
         efficiency_neg, chain_weights, suppress)
      class(integration_results_t), intent(inout) :: object
      integer, intent(in) :: n_it, n_calls, n_calls_valid
      real(default), intent(in) :: integral, error, efficiency, &
           efficiency_pos, efficiency_neg
      real(default), dimension(:), intent(in), optional :: chain_weights
      logical, intent(in), optional :: suppress
    end subroutine integration_results_record_extended
    module function integration_results_exist (results) result (flag)
      logical :: flag
      class(integration_results_t), intent(in) :: results
    end function integration_results_exist
    module function results_get_entry (results, last, it, pass) result (entry)
      class(integration_results_t), intent(in) :: results
      type(integration_entry_t) :: entry
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function results_get_entry
    module function integration_results_get_n_calls (results, last, it, pass) &
         result (n_calls)
      class(integration_results_t), intent(in), target :: results
      integer :: n_calls
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function integration_results_get_n_calls
    module function integration_results_get_integral (results, last, it, pass) &
         result (integral)
      class(integration_results_t), intent(in), target :: results
      real(default) :: integral
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function integration_results_get_integral
    module function integration_results_get_error (results, last, it, pass) &
         result (error)
      class(integration_results_t), intent(in), target :: results
      real(default) :: error
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function integration_results_get_error
    module function integration_results_get_accuracy (results, last, it, pass) &
         result (accuracy)
      class(integration_results_t), intent(in), target :: results
      real(default) :: accuracy
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function integration_results_get_accuracy
    module function integration_results_get_chi2 (results, last, it, pass) &
         result (chi2)
      class(integration_results_t), intent(in), target :: results
      real(default) :: chi2
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function integration_results_get_chi2
    module function integration_results_get_efficiency &
         (results, last, it, pass) result (efficiency)
      class(integration_results_t), intent(in), target :: results
      real(default) :: efficiency
      logical, intent(in), optional :: last
      integer, intent(in), optional :: it, pass
    end function integration_results_get_efficiency
    module subroutine integration_results_pacify (results, efficiency_reset)
      class(integration_results_t), intent(inout) :: results
      logical, intent(in), optional :: efficiency_reset
    end subroutine integration_results_pacify
    module subroutine integration_results_record_correction (object, corr, err)
      class(integration_results_t), intent(inout) :: object
      real(default), intent(in) :: corr, err
    end  subroutine integration_results_record_correction
    module subroutine integration_results_write_driver &
         (results, filename, eff_reset)
      type(integration_results_t), intent(inout) :: results
      type(string_t), intent(in) :: filename
      logical, intent(in), optional :: eff_reset
    end subroutine integration_results_write_driver
    module subroutine integration_results_compile_driver &
         (results, filename, os_data)
      type(integration_results_t), intent(in) :: results
      type(string_t), intent(in) :: filename
      type(os_data_t), intent(in) :: os_data
    end subroutine integration_results_compile_driver
  end interface

end module integration_results
