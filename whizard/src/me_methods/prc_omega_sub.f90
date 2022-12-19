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

submodule (prc_omega) prc_omega_s

  use io_units
  use system_defs, only: TAB

  implicit none

contains

  module function omega_def_type_string () result (string)
    type(string_t) :: string
    string = "omega"
  end function omega_def_type_string

  module subroutine omega_def_write (object, unit)
    class(omega_def_t), intent(in) :: object
    integer, intent(in) :: unit
    select type (writer => object%writer)
    class is (omega_writer_t)
       call writer%write (unit)
    end select
  end subroutine omega_def_write

  module subroutine omega_def_read (object, unit)
    class(omega_def_t), intent(out) :: object
    integer, intent(in) :: unit
    call msg_bug ("O'Mega process definition: input not supported yet")
  end subroutine omega_def_read

  module function omega_def_needs_code () result (flag)
    logical :: flag
    flag = .true.
  end function omega_def_needs_code

  module subroutine omega_def_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (6))
    features = [ &
         var_str ("init"), &
         var_str ("update_alpha_s"), &
         var_str ("reset_helicity_selection"), &
         var_str ("is_allowed"), &
         var_str ("new_event"), &
         var_str ("get_amplitude")]
  end subroutine omega_def_get_features

  module subroutine omega_def_connect (def, lib_driver, i, proc_driver)
    class(omega_def_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
    integer(c_int) :: pid, fid
    type(c_funptr) :: fptr
    select type (proc_driver)
    type is  (omega_driver_t)
       pid = i
       fid = 1
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%init)
       fid = 2
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%update_alpha_s)
       fid = 3
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%reset_helicity_selection)
       fid = 4
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%is_allowed)
       fid = 5
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%new_event)
       fid = 6
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%get_amplitude)
    end select
  end subroutine omega_def_connect

  module function omega_omega_writer_type_name () result (string)
    type(string_t) :: string
    string = "omega"
  end function omega_omega_writer_type_name

  module function omega_ufo_writer_type_name () result (string)
    type(string_t) :: string
    string = "omega/UFO"
  end function omega_ufo_writer_type_name

  module function omega_ovm_writer_type_name () result (string)
    type(string_t) :: string
    string = "ovm"
  end function omega_ovm_writer_type_name

  module function omega_writer_get_module_name (id) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: id
    name = "opr_" // id
  end function omega_writer_get_module_name

  module subroutine omega_writer_write (object, unit)
    class(omega_writer_t), intent(in) :: object
    integer, intent(in) :: unit
    write (unit, "(5x,A,A)")  "Model name        = ", &
         '"' // char (object%model_name) // '"'
    write (unit, "(5x,A,A)")  "Mode string       = ", &
         '"' // char (object%process_mode) // '"'
    write (unit, "(5x,A,A)")  "Process string    = ", &
         '"' // char (object%process_string) // '"'
    write (unit, "(5x,A,A)")  "Restrictions      = ", &
         '"' // char (object%restrictions) // '"'
    write (unit, "(5x,A,L1)")  "OpenMP support    = ", object%openmp_support
    write (unit, "(5x,A,L1)")  "Report progress   = ", object%report_progress
    ! write (unit, "(5x,A,L1)")  "Write phs output  = ", object%write_phs_output
    write (unit, "(5x,A,A)")  "Extra options     = ", &
         '"' // char (object%extra_options) // '"'
    write (unit, "(5x,A,L1)")  "Write diagrams    = ", object%diags
    write (unit, "(5x,A,L1)")  "Write color diag. = ", object%diags_color
    write (unit, "(5x,A,L1)")  "Complex Mass S.   = ", &
         object%complex_mass_scheme
  end subroutine omega_writer_write

  module subroutine omega_writer_init (writer, model_name, prt_in, prt_out, &
       ufo_path, restrictions, cms_scheme, openmp_support, &
       report_progress, write_phs_output, extra_options, diags, diags_color)
    class(omega_writer_t), intent(out) :: writer
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in
    type(string_t), dimension(:), intent(in) :: prt_out
    type(string_t), intent(in), optional :: ufo_path
    type(string_t), intent(in), optional :: restrictions
    logical, intent(in), optional :: cms_scheme
    logical, intent(in), optional :: openmp_support
    logical, intent(in), optional :: report_progress
    logical, intent(in), optional :: write_phs_output
    type(string_t), intent(in), optional :: extra_options
    logical, intent(in), optional :: diags, diags_color
    integer :: i
    writer%model_name = model_name
    select type (writer)
    type is (omega_ufo_writer_t)
       if (present (ufo_path)) then
          writer%ufo_path = ufo_path
       else
          call msg_fatal ("O'Mega: UFO model option is selected, but UFO model path is unset")
       end if
    end select
    if (present (restrictions)) then
       writer%restrictions = restrictions
    else
       writer%restrictions = ""
    end if
    if (present (cms_scheme))  writer%complex_mass_scheme = cms_scheme
    if (present (openmp_support))  writer%openmp_support = openmp_support
    if (present (report_progress))  writer%report_progress = report_progress
    if (present (write_phs_output)) writer%write_phs_output = write_phs_output
    if (present (extra_options)) then
       writer%extra_options = " " // extra_options
    else
       writer%extra_options = ""
    end if
    if (present (diags))  writer%diags = diags
    if (present (diags_color))  writer%diags_color = diags_color
    select case (size (prt_in))
    case (1);  writer%process_mode = " -decay"
    case (2);  writer%process_mode = " -scatter"
    end select
    associate (s => writer%process_string)
      s = " '"
      do i = 1, size (prt_in)
         if (i > 1)  s = s // " "
         s = s // prt_in(i)
      end do
      s = s // " ->"
      do i = 1, size (prt_out)
         s = s // " " // prt_out(i)
      end do
      s = s // "'"
    end associate
  end subroutine omega_writer_init

  module subroutine omega_write_makefile_code &
       (writer, unit, id, os_data, verbose, testflag)
    class(omega_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    type(string_t) :: omega_binary, omega_path
    type(string_t) :: restrictions_string
    type(string_t) :: openmp_string
    type(string_t) :: kmatrix_string
    type(string_t) :: progress_string
    type(string_t) :: diagrams_string
    type(string_t) :: cms_string
    type(string_t) :: write_phs_output_string
    type(string_t) :: parameter_module
    logical :: escape_hyperref
    escape_hyperref = .false.
    if (present (testflag))  escape_hyperref = testflag
    select type (writer)
    type is (omega_omega_writer_t)
       omega_binary = "omega_" // writer%model_name // ".opt"
    type is (omega_ufo_writer_t)
       omega_binary = "omega_UFO.opt"
    type is (omega_ovm_writer_t)
       select case (char (writer%model_name))
       case ("SM", "SM_CKM", "SM_Higgs", "THDM", "THDM_CKM", &
             "HSExt", "QED", "QCD", "Zprime")
       case default
          call msg_fatal ("The model " // char (writer%model_name) &
               // " is not available for the O'Mega VM.")
       end select
       omega_binary = "omega_" // writer%model_name // "_VM.opt"
    end select
    omega_path = os_data%whizard_omega_binpath // "/" // omega_binary
    if (.not. verbose)  omega_path = "@" // omega_path
    if (writer%restrictions /= "") then
       restrictions_string = " -cascade '" // writer%restrictions // "'"
    else
       restrictions_string = ""
    end if
    if (writer%openmp_support) then
       openmp_string = " -target:openmp"
    else
       openmp_string = ""
    end if
    if (writer%report_progress) then
       progress_string = " -fusion:progress"
    else
       progress_string = ""
    end if
    if (writer%diags) then
       if (writer%diags_color) then
          diagrams_string = " -diagrams:C " // char(id) // &
               "_diags -diagrams_LaTeX"
       else
          diagrams_string = " -diagrams " // char(id) // &
               "_diags -diagrams_LaTeX"
       end if
    else
       if (writer%diags_color) then
          diagrams_string = " -diagrams:c " // char(id) // &
               "_diags -diagrams_LaTeX"
       else
          diagrams_string = ""
       end if
    end if
    if (writer%complex_mass_scheme) then
       cms_string = " -model:cms_width"
    else
       cms_string = ""
    end if
    if (writer%write_phs_output) then
       write_phs_output_string = " -phase_space " // char (id) // ".fds"
    else
       write_phs_output_string = ""
    endif
    select case (char (writer%model_name))
    case ("SM_rx", "SSC", "NoH_rx", "AltH")
       kmatrix_string = " -target:kmatrix_2_write"
    case ("SSC_2", "SSC_AltT", "SM_ul")
       kmatrix_string = " -target:kmatrix_write"
    case default
       kmatrix_string = ""
    end select
    write (unit, "(5A)")  "SOURCES += ", char (id), ".f90"
    select type (writer)
    type is (omega_ovm_writer_t)
       write (unit, "(5A)")  "SOURCES += ", char (id), ".hbc"
    end select
    if (writer%diags .or. writer%diags_color) then
       write (unit, "(5A)")  "TEX_SOURCES += ", char (id), "_diags.tex"
       if (os_data%event_analysis_pdf) then
          write (unit, "(5A)")  "TEX_OBJECTS += ", char (id), "_diags.pdf"
       else
          write (unit, "(5A)")  "TEX_OBJECTS += ", char (id), "_diags.ps"
       end if
    end if
    write (unit, "(5A)")  "OBJECTS += ", char (id), ".lo"
    select type (writer)
    type is (omega_omega_writer_t)
       write (unit, "(5A)")  char (id), ".f90:"
       if (.not. verbose) then
          write (unit, "(5A)")  TAB // '@echo  "  OMEGA     ', trim (char (id)), '.f90"'
       end if
       write (unit, "(99A)")  TAB, char (omega_path), &
            " -o ", char (id), ".f90", &
            " -target:whizard", &
            " -target:parameter_module parameters_", char (writer%model_name), &
            " -target:module opr_", char (id), &
            " -target:md5sum '", writer%md5sum, "'", &
            char (cms_string), &
            char (openmp_string), &
            char (progress_string), &
            char (kmatrix_string), &
            char (writer%process_mode), char (writer%process_string), &
            char (restrictions_string), char (diagrams_string), &
            char (writer%extra_options), char (write_phs_output_string)
    type is (omega_ufo_writer_t)
       parameter_module = char (id) // "_par_" // replace (char (writer%model_name), &
            "-", "_", every=.true.)
       write (unit, "(5A)")  char (id), ".f90: ", char (parameter_module), ".lo"
       if (.not. verbose) then
          write (unit, "(5A)")  TAB // '@echo  "  OMEGA[UFO]', trim (char (id)), '.f90"'
       end if
       write (unit, "(99A)")  TAB, char (omega_path), &
            " -o ", char (id), ".f90", &
            " -model:UFO_dir ", &
            char (writer%ufo_path), "/", char (writer%model_name), &
            " -model:exec", &
            " -target:whizard", &
            " -target:parameter_module ", char (parameter_module), &
            " -target:module opr_", char (id), &
            " -target:md5sum '", writer%md5sum, "'", &
            char (cms_string), &
            char (openmp_string), &
            char (progress_string), &
            char (kmatrix_string), &
            char (writer%process_mode), char (writer%process_string), &
            char (restrictions_string), char (diagrams_string), &
            char (writer%extra_options), char (write_phs_output_string)
       write (unit, "(5A)") "SOURCES += ", char (parameter_module), ".f90"
       write (unit, "(5A)") "OBJECTS += ", char (parameter_module), ".lo"
       write (unit, "(5A)")  char (parameter_module), ".f90:"
       write (unit, "(99A)")  TAB, char (omega_path), &
            " -model:UFO_dir ", &
            char (writer%ufo_path), "/", char (writer%model_name), &
            " -model:exec", &
            " -target:parameter_module ", char (parameter_module), &
            " -params", &
            " -o $@"
       write (unit, "(5A)")  char (parameter_module), ".lo: ", char (parameter_module), ".f90"
       if (.not. verbose) then
          write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
       end if
       write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"
    type is (omega_ovm_writer_t)
       write (unit, "(5A)")  char (id), ".hbc:"
       write (unit, "(99A)")  TAB, char (omega_path), &
            " -o ", char (id), ".hbc", &
            char (progress_string), &
            char (cms_string), &
            char (writer%process_mode), char (writer%process_string), &
            char (restrictions_string), char (diagrams_string), &
            char (writer%extra_options), char (write_phs_output_string)
       write (unit, "(5A)")  char (id), ".f90:"
       if (.not. verbose) then
          write (unit, "(5A)")  TAB // '@echo  "  OVM       ', trim (char (id)), '.f90"'
       end if
       write (unit, "(99A)")  TAB, char (omega_path), &
            " -o ", char (id), ".f90 -params", &
            " -target:whizard ", &
            " -target:bytecode_file ", char (id), ".hbc", &
            " -target:wrapper_module opr_", char (id), &
            " -target:parameter_module_external parameters_", &
            char (writer%model_name), &
            " -target:md5sum '", writer%md5sum, "'", &
            char (openmp_string)
    end select
    if (writer%diags .or. writer%diags_color) &
       write (unit, "(5A)")  char (id), "_diags.tex: ", char (id), ".f90"
    write (unit, "(5A)")  "clean-", char (id), ":"
    if (verbose) then
       write (unit, "(5A)")  TAB, "rm -f ", char (id), ".f90"
       write (unit, "(5A)")  TAB, "rm -f opr_", char (id), ".mod"
       write (unit, "(5A)")  TAB, "rm -f ", char (id), ".lo"
    else
       write (unit, "(5A)")  TAB // '@echo  "  RM        ', &
            trim (char (id)), '.f90,.mod,.lo"'
       write (unit, "(5A)")  TAB, "@rm -f ", char (id), ".f90"
       write (unit, "(5A)")  TAB, "@rm -f opr_", char (id), ".mod"
       write (unit, "(5A)")  TAB, "@rm -f ", char (id), ".lo"
    end if
    write (unit, "(5A)")  "CLEAN_SOURCES += ", char (id), ".f90"
    select type (writer)
    type is (omega_ufo_writer_t)
       write (unit, "(5A)")  "CLEAN_SOURCES += ", char (writer%model_name), ".mdl"
       write (unit, "(5A)")  "CLEAN_SOURCES += ", char (parameter_module), ".f90"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (parameter_module), ".mod"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (parameter_module), ".lo"
    type is (omega_ovm_writer_t)
       write (unit, "(5A)")  "CLEAN_SOURCES += ", char (id), ".hbc"
    end select
    if (writer%diags .or. writer%diags_color) then
       write (unit, "(5A)")  "CLEAN_SOURCES += ", char (id), "_diags.tex"
    end if
    write (unit, "(5A)")  "CLEAN_OBJECTS += opr_", char (id), ".mod"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), ".lo"
    if (writer%diags .or. writer%diags_color) then
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.aux"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.log"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.dvi"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.toc"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.out"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.[1-9]"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.[1-9][0-9]"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.[1-9][0-9][0-9]"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.t[1-9]"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.t[1-9][0-9]"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.t[1-9][0-9][0-9]"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.mp"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags-fmf.log"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.dvi"
       write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.ps"
       if (os_data%event_analysis_pdf) &
            write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_diags.pdf"
    end if
    write (unit, "(5A)")  char (id), ".lo: ", char (id), ".f90"
    write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"
    if (writer%diags .or. writer%diags_color) then
       if (os_data%event_analysis_ps) then
          if (os_data%event_analysis_pdf) then
             write (unit, "(5A)")  char (id), "_diags.pdf: ", char (id), "_diags.tex"
          else
             write (unit, "(5A)")  char (id), "_diags.ps: ", char (id), "_diags.tex"
          end if
          if (escape_hyperref) then
             if (verbose) then
                write (unit, "(5A)")  TAB, "-cat ", char (id), "_diags.tex | \"
             else
                write (unit, "(5A)")  TAB // '@echo  "  HYPERREF  ', &
                     trim (char (id)) // '_diags.tex"'
                write (unit, "(5A)")  TAB, "@cat ", char (id), "_diags.tex | \"
             end if
             write (unit, "(5A)")  TAB, "   sed -e" // &
                  "'s/\\usepackage\[colorlinks\]{hyperref}.*/%\\usepackage" // &
                  "\[colorlinks\]{hyperref}/' > \"
             write (unit, "(5A)")  TAB, "   ", char (id), "_diags.tex.tmp"
             if (verbose) then
                write (unit, "(5A)")  TAB, "mv -f ", char (id), "_diags.tex.tmp \"
             else
                write (unit, "(5A)")  TAB, "@mv -f ", char (id), "_diags.tex.tmp \"
             end if
             write (unit, "(5A)")  TAB, "   ", char (id), "_diags.tex"
          end if
          if (verbose) then
             write (unit, "(5A)")  TAB, "-TEXINPUTS=$(TEX_FLAGS) $(LATEX) " // &
                  char (id) // "_diags.tex"
             write (unit, "(5A)")  TAB, "MPINPUTS=$(MP_FLAGS) $(MPOST) " // &
                  char (id) // "_diags-fmf.mp"
             write (unit, "(5A)")  TAB, "TEXINPUTS=$(TEX_FLAGS) $(LATEX) " // &
                  char (id) // "_diags.tex"
             write (unit, "(5A)")  TAB, "$(DVIPS) -o " // char (id) // "_diags.ps " // &
                  char (id) // "_diags.dvi"
          else
             write (unit, "(5A)")  TAB // '@echo  "  LATEX     ', &
                  trim (char (id)) // '_diags.tex"'
             write (unit, "(5A)")  TAB, "@TEXINPUTS=$(TEX_FLAGS) $(LATEX) " // &
                  char (id) // "_diags.tex > /dev/null"
             write (unit, "(5A)")  TAB // '@echo  "  METAPOST  ', &
                  trim (char (id)) // '_diags-fmf.mp"'
             write (unit, "(5A)")  TAB, "@MPINPUTS=$(MP_FLAGS) $(MPOST) " // &
                  char (id) // "_diags-fmf.mp  > /dev/null"
             write (unit, "(5A)")  TAB // '@echo  "  LATEX     ', &
                  trim (char (id)) // '_diags.tex"'
             write (unit, "(5A)")  TAB, "@TEXINPUTS=$(TEX_FLAGS) $(LATEX) " // &
                  char (id) // "_diags.tex > /dev/null"
             write (unit, "(5A)")  TAB // '@echo  "  DVIPS     ', &
                  trim (char (id)) // '_diags.dvi"'
             write (unit, "(5A)")  TAB, "@$(DVIPS) -q -o " // char (id) &
                  // "_diags.ps " // char (id) // "_diags.dvi"
          end if
          if (os_data%event_analysis_pdf) then
             if (verbose) then
                write (unit, "(5A)")  TAB, "$(PS2PDF) " // char (id) // "_diags.ps"
             else
                write (unit, "(5A)")  TAB // '@echo  "  PS2PDF    ', &
                     trim (char (id)) // '_diags.ps"'
                write (unit, "(5A)")  TAB, "@$(PS2PDF) " // char (id) // "_diags.ps"
             end if
          end if
       end if
    end if
  end subroutine omega_write_makefile_code

  module subroutine omega_write_source_code (writer, id)
    class(omega_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine omega_write_source_code

  module subroutine omega_before_compile (writer, id)
    class(omega_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine omega_before_compile

  module subroutine omega_after_compile (writer, id)
    class(omega_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
  end subroutine omega_after_compile

  module function omega_writer_get_procname (feature) result (name)
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
  end function omega_writer_get_procname

  module subroutine omega_write_interface (writer, unit, id, feature)
    class(omega_writer_t), intent(in) :: writer
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
    case ("reset_helicity_selection")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " &
            &(threshold, cutoff) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), intent(in) :: threshold"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: cutoff"
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
  end subroutine omega_write_interface

  module subroutine omega_write_wrapper (writer, unit, id, feature)
    class(omega_writer_t), intent(in) :: writer
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
       write (unit, "(2x,9A)")  "use opr_", char (id)
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
       write (unit, "(2x,9A)")  "use opr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), intent(in) &
               &:: alpha_s"
          write (unit, "(2x,9A)")  "call ", char (feature), " (alpha_s)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("reset_helicity_selection")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (threshold, cutoff) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), intent(in) &
               &:: threshold"
          write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: cutoff"
          write (unit, "(2x,9A)")  "call ", char (feature), &
               " (threshold, int (cutoff))"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("is_allowed")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (flv, hel, col, flag) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(2x,9A)")  "logical(c_bool), intent(out) :: flag"
       write (unit, "(2x,9A)")  "flag = ", char (feature), &
            " (int (flv), int (hel), int (col))"
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("new_event")
       write (unit, "(9A)")  "subroutine ", char (name), " (p) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
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
       write (unit, "(2x,9A)")  "use opr_", char (id)
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(2x,9A)")  "complex(c_default_complex), intent(out) &
            &:: amp"
       write (unit, "(2x,9A)")  "amp = ", char (feature), &
            " (int (flv), int (hel), int (col))"
       write (unit, "(9A)")  "end subroutine ", char (name)
    end select
  end subroutine omega_write_wrapper

  module function omega_driver_type_name () result (string)
    type(string_t) :: string
    string = "omega"
  end function omega_driver_type_name

  module subroutine omega_state_write (object, unit)
    class(omega_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,L1)")  "O'Mega state: new kinematics = ", &
         object%new_kinematics
  end subroutine omega_state_write

  module subroutine omega_state_reset_new_kinematics (object)
    class(omega_state_t), intent(inout) :: object
    object%new_kinematics = .true.
  end subroutine omega_state_reset_new_kinematics

  module subroutine prc_omega_write (object, unit)
    class(prc_omega_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(3x,A)", advance="no")  "O'Mega process core:"
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
    call object%helicity_selection%write (u)
    call object%qcd%write (u)
    call object%qed%write (u)
  end subroutine prc_omega_write

  module subroutine prc_omega_write_name (object, unit)
    class(prc_omega_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: O'Mega"
  end subroutine prc_omega_write_name

  module subroutine prc_omega_init (object, def, lib, id, i_component)
    class(prc_omega_t), intent(inout) :: object
    class(prc_core_def_t), intent(in), target :: def
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    call object%base_init (def, lib, id, i_component)
    call object%activate_parameters ()
  end subroutine prc_omega_init

  module subroutine prc_omega_activate_parameters (object)
    class (prc_omega_t), intent(inout) :: object
    if (allocated (object%driver)) then
       if (allocated (object%par)) then
          select type (driver => object%driver)
          type is (omega_driver_t)
             if (associated (driver%init)) then
                call driver%init (object%par, object%scheme)
             end if
          end select
       else
          call msg_bug ("prc_omega_activate: parameter set is not allocated")
       end if
       call object%reset_helicity_selection ()
    else
       call msg_bug ("prc_omega_activate: driver is not allocated")
    end if
  end subroutine prc_omega_activate_parameters

  module function prc_omega_is_allowed (object, i_term, f, h, c) result (flag)
    class(prc_omega_t), intent(in) :: object
    integer, intent(in) :: i_term, f, h, c
    logical :: flag
    logical(c_bool) :: cflag
    select type (driver => object%driver)
    type is (omega_driver_t)
       call driver%is_allowed (f, h, c, cflag)
       flag = cflag
    end select
  end function prc_omega_is_allowed

  module subroutine prc_omega_compute_hard_kinematics &
       (object, p_seed, i_term, int_hard, core_state)
    class(prc_omega_t), intent(in) :: object
    type(vector4_t), dimension(:), intent(in) :: p_seed
    integer, intent(in) :: i_term
    type(interaction_t), intent(inout) :: int_hard
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    call int_hard%set_momenta (p_seed)
  end subroutine prc_omega_compute_hard_kinematics

  module subroutine prc_omega_compute_eff_kinematics &
       (object, i_term, int_hard, int_eff, core_state)
    class(prc_omega_t), intent(in) :: object
    integer, intent(in) :: i_term
    type(interaction_t), intent(in) :: int_hard
    type(interaction_t), intent(inout) :: int_eff
    class(prc_core_state_t), intent(inout), allocatable :: core_state
  end subroutine prc_omega_compute_eff_kinematics

  module subroutine prc_omega_reset_helicity_selection (object)
    class(prc_omega_t), intent(inout) :: object
    select type (driver => object%driver)
    type is (omega_driver_t)
       if (associated (driver%reset_helicity_selection)) then
          if (object%helicity_selection%active) then
             call driver%reset_helicity_selection &
                  (real (object%helicity_selection%threshold, &
                  c_default_float), &
                  int (object%helicity_selection%cutoff, c_int))
          else
             call driver%reset_helicity_selection &
                  (0._c_default_float, 0_c_int)
          end if
       end if
    end select
  end subroutine prc_omega_reset_helicity_selection

  module function prc_omega_compute_amplitude &
       (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
       core_state)  result (amp)
    class(prc_omega_t), intent(in) :: object
    integer, intent(in) :: j
    type(vector4_t), dimension(:), intent(in) :: p
    integer, intent(in) :: f, h, c
    real(default), intent(in) :: fac_scale, ren_scale
    real(default), intent(in), allocatable :: alpha_qcd_forced
    class(prc_core_state_t), intent(inout), allocatable, optional :: core_state
    real(default) :: alpha_qcd
    complex(default) :: amp
    integer :: n_tot, i
    real(c_default_float), dimension(:,:), allocatable :: parray
    complex(c_default_complex) :: camp
    logical :: new_event
    select type (driver => object%driver)
    type is (omega_driver_t)
       new_event = .true.
       if (present (core_state)) then
          if (allocated (core_state)) then
             select type (core_state)
             type is (omega_state_t)
                new_event = core_state%new_kinematics
                core_state%new_kinematics = .false.
             end select
          end if
       end if
       if (new_event) then
          if (allocated (object%qcd%alpha)) then
             if (allocated (alpha_qcd_forced)) then
                alpha_qcd = alpha_qcd_forced
             else
                alpha_qcd = object%qcd%alpha%get (ren_scale)
             end if
             call driver%update_alpha_s (alpha_qcd)
             if (present (core_state)) then
                if (allocated (core_state)) then
                   select type (core_state)
                   type is (omega_state_t)
                      core_state%alpha_qcd = alpha_qcd
                   end select
                end if
             end if
          end if
          n_tot = object%data%get_n_tot ()
          allocate (parray (0:3, n_tot))
          do i = 1, n_tot
             parray(:,i) = vector4_get_components (p(i))
          end do
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
  end function prc_omega_compute_amplitude

  module function prc_omega_get_alpha_s &
       (object, core_state) result (alpha_qcd)
    class(prc_omega_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(default) :: alpha_qcd
    alpha_qcd = -1
    if (allocated (object%qcd%alpha) .and. allocated (core_state)) then
       select type (core_state)
       type is (omega_state_t)
          alpha_qcd = core_state%alpha_qcd
       end select
    end if
  end function prc_omega_get_alpha_s

  module function prc_omega_get_alpha_qed &
       (object, core_state) result (alpha_qed)
    class(prc_omega_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(default) :: alpha_qed
    alpha_qed = -1
    if (allocated (object%qed%alpha) .and. allocated (core_state)) then
       select type (core_state)
       type is (omega_state_t)
          alpha_qed = core_state%alpha_qed
       end select
    end if
  end function prc_omega_get_alpha_qed


end submodule prc_omega_s

