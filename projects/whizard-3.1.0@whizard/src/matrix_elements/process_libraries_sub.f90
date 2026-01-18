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

submodule (process_libraries) process_libraries_s

  use io_units
  use diagnostics
  use md5

  implicit none

contains

  module subroutine strip_equation_lhs (buffer)
    character(*), intent(inout) :: buffer
    type(string_t) :: string, prefix
    string = buffer
    call split (string, prefix, "=")
    buffer = string
  end subroutine strip_equation_lhs

  module subroutine process_component_def_write (object, unit)
    class(process_component_def_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(3x,A,A)")  "Component ID        = ", char (object%basename)
    write (u, "(3x,A,L1)") "Initial component   = ", object%initial
    write (u, "(3x,A,I0,1x,I0,1x,I0)") "N (in, out, tot)    = ", &
         object%n_in, object%n_out, object%n_tot
    write (u, "(3x,A)", advance="no") "Particle content    = "
    if (allocated (object%prt_in)) then
       call prt_spec_write (object%prt_in, u, advance="no")
    else
       write (u, "(A)", advance="no")  "[undefined]"
    end if
    write (u, "(A)", advance="no") " => "
    if (allocated (object%prt_out)) then
       call prt_spec_write (object%prt_out, u, advance="no")
    else
       write (u, "(A)", advance="no")  "[undefined]"
    end if
    write (u, "(A)")
    if (object%method /= "") then
       write (u, "(3x,A,A)")  "Method              = ", &
            char (object%method)
    else
       write (u, "(3x,A)")  "Method              = [undefined]"
    end if
    if (allocated (object%core_def)) then
       write (u, "(3x,A,A)")  "Process variant     = ", &
            char (object%core_def%type_string ())
       call object%core_def%write (u)
    else
       write (u, "(3x,A)")  "Process variant     = [undefined]"
    end if
    write (u, "(3x,A,A,A)") "MD5 sum (def)       = '", object%md5sum, "'"
  end subroutine process_component_def_write

  module subroutine process_component_def_read (component, unit, core_def_templates)
    class(process_component_def_t), intent(out) :: component
    integer, intent(in) :: unit
    type(prc_template_t), dimension(:), intent(in) :: core_def_templates
    character(80) :: buffer
    type(string_t) :: var_buffer, prefix, in_state, out_state
    type(string_t) :: variant_type

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    component%basename = trim (adjustl (buffer))

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    read (buffer, *)  component%initial

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    read (buffer, *)  component%n_in, component%n_out, component%n_tot

    call get (unit, var_buffer)
    call split (var_buffer, prefix, "=")   ! keeps 'in => out'
    call split (var_buffer, prefix, "=")   ! actually: separator is '=>'

    in_state = prefix
    if (component%n_in > 0) then
       call prt_spec_read (component%prt_in, in_state)
    end if

    out_state = extract (var_buffer, 2)
    if (component%n_out > 0) then
       call prt_spec_read (component%prt_out, out_state)
    end if

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    component%method = trim (adjustl (buffer))
    if (component%method == "[undefined]") &
         component%method = ""

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    variant_type = trim (adjustl (buffer))
    call allocate_core_def &
         (core_def_templates, variant_type, component%core_def)
    if (allocated (component%core_def)) then
       call component%core_def%read (unit)
    end if

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    read (buffer(3:34), "(A32)")  component%md5sum

  end subroutine process_component_def_read

  module subroutine process_component_def_show (object, unit)
    class(process_component_def_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(6x,A)", advance="no")  char (object%basename)
    if (.not. object%initial) &
         write (u, "('*')", advance="no")
    write (u, "(':',1x)", advance="no")
    if (allocated (object%prt_in)) then
       call prt_spec_write (object%prt_in, u, advance="no")
    else
       write (u, "(A)", advance="no")  "[undefined]"
    end if
    write (u, "(A)", advance="no") " => "
    if (allocated (object%prt_out)) then
       call prt_spec_write (object%prt_out, u, advance="no")
    else
       write (u, "(A)", advance="no")  "[undefined]"
    end if
    if (object%method /= "") then
       write (u, "(2x,'[',A,']')")  char (object%method)
    else
       write (u, *)
    end if
  end subroutine process_component_def_show

  module subroutine process_component_def_compute_md5sum (component, model)
    class(process_component_def_t), intent(inout) :: component
    class(model_data_t), intent(in), optional, target :: model
    integer :: u
    component%md5sum = ""
    u = free_unit ()
    open (u, status = "scratch", action = "readwrite")
    if (present (model))  write (u, "(A32)")  model%get_md5sum ()
    call component%write (u)
    rewind (u)
    component%md5sum = md5sum (u)
    close (u)
    if (allocated (component%core_def)) then
       call component%core_def%set_md5sum (component%md5sum)
    end if
  end subroutine process_component_def_compute_md5sum

  module function process_component_def_get_def_type_string (component) result (type_string)
    type(string_t) :: type_string
    class(process_component_def_t), intent(in) :: component
    type_string = component%core_def%type_string ()
  end function process_component_def_get_def_type_string

  module subroutine process_component_def_allocate_driver (component, driver)
    class(process_component_def_t), intent(in) :: component
    class(prc_core_driver_t), intent(out), allocatable :: driver
    if (allocated (component%core_def)) then
       call component%core_def%allocate_driver (driver, component%basename)
    end if
  end subroutine process_component_def_allocate_driver

  module function process_component_def_needs_code (component) result (flag)
    class(process_component_def_t), intent(in) :: component
    logical :: flag
    flag = component%core_def%needs_code ()
  end function process_component_def_needs_code

  module function process_component_def_get_writer_ptr (component) result (writer)
    class(process_component_def_t), intent(in), target :: component
    class(prc_writer_t), pointer :: writer
    writer => component%core_def%writer
  end function process_component_def_get_writer_ptr

  module function process_component_def_get_features (component) result (features)
    class(process_component_def_t), intent(in) :: component
    type(string_t), dimension(:), allocatable :: features
    call component%core_def%get_features (features)
  end function process_component_def_get_features

  module subroutine process_component_def_connect &
       (component, lib_driver, i, proc_driver)
    class(process_component_def_t), intent(in) :: component
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
    select type (proc_driver)
    class is (process_driver_internal_t)
       !!! Nothing to do
    class default
       call component%core_def%connect (lib_driver, i, proc_driver)
    end select
  end subroutine process_component_def_connect

  module function process_component_get_core_def_ptr (component) result (ptr)
    class(process_component_def_t), intent(in), target :: component
    class(prc_core_def_t), pointer :: ptr
    ptr => component%core_def
  end function process_component_get_core_def_ptr

  module function process_component_def_get_n_in (component) result (n_in)
    class(process_component_def_t), intent(in) :: component
    integer :: n_in
    n_in = component%n_in
  end function process_component_def_get_n_in

  module function process_component_def_get_n_out (component) result (n_out)
    class(process_component_def_t), intent(in) :: component
    integer :: n_out
    n_out = component%n_out
  end function process_component_def_get_n_out

  module function process_component_def_get_n_tot (component) result (n_tot)
    class(process_component_def_t), intent(in) :: component
    integer :: n_tot
    n_tot = component%n_tot
  end function process_component_def_get_n_tot

  module subroutine process_component_def_get_prt_in (component, prt)
    class(process_component_def_t), intent(in) :: component
    type(string_t), dimension(:), intent(out), allocatable :: prt
    integer :: i
    allocate (prt (component%n_in))
    do i = 1, component%n_in
       prt(i) = component%prt_in(i)%to_string ()
    end do
  end subroutine process_component_def_get_prt_in

  module subroutine process_component_def_get_prt_out (component, prt)
    class(process_component_def_t), intent(in) :: component
    type(string_t), dimension(:), intent(out), allocatable :: prt
    integer :: i
    allocate (prt (component%n_out))
    do i = 1, component%n_out
       prt(i) = component%prt_out(i)%to_string ()
    end do
  end subroutine process_component_def_get_prt_out

  module function process_component_def_get_prt_spec_in (component) result (prt)
    class(process_component_def_t), intent(in) :: component
    type(prt_spec_t), dimension(:), allocatable :: prt
    allocate (prt (component%n_in))
    prt(:) = component%prt_in(:)
  end function process_component_def_get_prt_spec_in

  module function process_component_def_get_prt_spec_out (component) result (prt)
    class(process_component_def_t), intent(in) :: component
    type(prt_spec_t), dimension(:), allocatable :: prt
    allocate (prt (component%n_out))
    prt(:) = component%prt_out(:)
  end function process_component_def_get_prt_spec_out

  module subroutine process_component_def_get_pdg_in (component, model, pdg)
    class(process_component_def_t), intent(in) :: component
    class(model_data_t), intent(in), target :: model
    integer, intent(out), dimension(:) :: pdg
    integer :: i
    do i = 1, size (pdg)
       pdg(i) = model%get_pdg (component%prt_in(i)%to_string ())
    end do
  end subroutine process_component_def_get_pdg_in

  pure module function process_component_def_get_md5sum (component) result (md5sum)
    class(process_component_def_t), intent(in) :: component
    character(32) :: md5sum
    md5sum = component%md5sum
  end function process_component_def_get_md5sum

  elemental module function process_component_def_get_nlo_type &
       (component) result (nlo_type)
    integer :: nlo_type
    class(process_component_def_t), intent(in) :: component
    nlo_type = component%nlo_type
  end function process_component_def_get_nlo_type

  elemental module function process_component_def_get_associated_born &
       (component) result (i_born)
    integer :: i_born
    class(process_component_def_t), intent(in) :: component
    i_born = component%associated_components(ASSOCIATED_BORN)
  end function process_component_def_get_associated_born

  elemental module function process_component_def_get_associated_real_fin &
       (component) result (i_rfin)
    integer :: i_rfin
    class(process_component_def_t), intent(in) :: component
    i_rfin = component%associated_components(ASSOCIATED_REAL_FIN)
  end function process_component_def_get_associated_real_fin

  elemental module function process_component_def_get_associated_real_sing &
       (component) result (i_rsing)
    integer :: i_rsing
    class(process_component_def_t), intent(in) :: component
    i_rsing = component%associated_components(ASSOCIATED_REAL_SING)
  end function process_component_def_get_associated_real_sing

  elemental module function process_component_def_get_associated_subtraction &
       (component) result (i_sub)
    integer :: i_sub
    class(process_component_def_t), intent(in) :: component
    i_sub = component%associated_components(ASSOCIATED_SUB)
  end function process_component_def_get_associated_subtraction

  elemental module function process_component_def_can_be_integrated &
       (component) result (active)
    logical :: active
    class(process_component_def_t), intent(in) :: component
    active = component%active
  end function process_component_def_can_be_integrated

  module function process_component_def_get_association_list &
       (component, i_skip_in) result (list)
    integer, dimension(:), allocatable :: list
    class(process_component_def_t), intent(in) :: component
    integer, intent(in), optional :: i_skip_in
    integer :: i, j, n, i_skip
    logical :: valid
    i_skip = 0; if (present (i_skip_in)) i_skip = i_skip_in
    n = count (component%associated_components /= 0) - 1
    if (i_skip > 0) then
       if (component%associated_components(i_skip) > 0) n = n - 1
    end if
    allocate (list (n))
    j = 1
    do i = 1, size(component%associated_components)
       valid = component%associated_components(i) /= 0 &
               .and. i /= ASSOCIATED_SUB .and. i /= i_skip
       if (valid) then
          list(j) = component%associated_components(i)
          j = j + 1
       end if
    end do
  end function process_component_def_get_association_list

  module function process_component_def_get_associated_real &
       (component) result (i_real)
    integer :: i_real
    class(process_component_def_t), intent(in) :: component
    i_real = component%associated_components(ASSOCIATED_REAL)
  end function process_component_def_get_associated_real

  elemental module function process_component_def_get_me_method (component) result (method)
    type(string_t) :: method
    class(process_component_def_t), intent(in) :: component
    method = component%method
  end function process_component_def_get_me_method

  module function process_component_def_get_fixed_emitter (component) result (emitter)
    integer :: emitter
    class(process_component_def_t), intent(in) :: component
    emitter = component%fixed_emitter
  end function process_component_def_get_fixed_emitter

  pure module subroutine process_component_def_get_coupling_powers &
       (component, alpha_power, alphas_power)
    class(process_component_def_t), intent(in) :: component
    integer, intent(out) :: alpha_power, alphas_power
    alpha_power = component%alpha_power
    alphas_power = component%alphas_power
  end subroutine process_component_def_get_coupling_powers

  module subroutine process_def_write (object, unit)
    class(process_def_t), intent(in) :: object
    integer, intent(in) :: unit
    integer :: i
    write (unit, "(1x,A,A,A)") "ID = '", char (object%id), "'"
    if (object%num_id /= 0) &
         write (unit, "(1x,A,I0)")  "ID(num) = ", object%num_id
    select case (object%n_in)
    case (1);  write (unit, "(1x,A)")  "Decay"
    case (2);  write (unit, "(1x,A)")  "Scattering"
    case default
       write (unit, "(1x,A)")  "[Undefined process]"
       return
    end select
    if (object%model_name /= "") then
       write (unit, "(1x,A,A)")  "Model = ", char (object%model_name)
    else
       write (unit, "(1x,A)")  "Model = [undefined]"
    end if
    write (unit, "(1x,A,I0)")  "Initially defined component(s) = ", &
         object%n_initial
    write (unit, "(1x,A,I0)")  "Extra generated component(s)   = ", &
         object%n_extra
    if (object%requires_resonances) then
       ! This line has to matched with the reader below!
       write (unit, "(1x,A,I0)")  "Resonant subprocesses required"
    end if
    write (unit, "(1x,A,A,A)") "MD5 sum   = '", object%md5sum, "'"
    if (allocated (object%initial)) then
       do i = 1, size (object%initial)
          write (unit, "(1x,A,I0)")  "Component #", i
          call object%initial(i)%write (unit)
       end do
    end if
    if (allocated (object%extra)) then
       do i = 1, size (object%extra)
          write (unit, "(1x,A,I0)")  "Component #", object%n_initial + i
          call object%extra(i)%write (unit)
       end do
    end if
  end subroutine process_def_write

  module subroutine process_def_read (object, unit, core_def_templates)
    class(process_def_t), intent(out) :: object
    integer, intent(in) :: unit
    type(prc_template_t), dimension(:), intent(in) :: core_def_templates
    integer :: i, i1, i2
    character(80) :: buffer, ref
    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    i1 = scan (buffer, "'")
    i2 = scan (buffer, "'", back=.true.)
    if (i2 > i1) then
       object%id = buffer(i1+1:i2-1)
    else
       object%id = ""
    end if

    read (unit, "(A)")  buffer
    select case (buffer(2:11))
    case ("Decay     "); object%n_in = 1
    case ("Scattering"); object%n_in = 2
    case default
       return
    end select

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    object%model_name = trim (adjustl (buffer))
    if (object%model_name == "[undefined]")  object%model_name = ""

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    read (buffer, *)  object%n_initial

    read (unit, "(A)")  buffer
    call strip_equation_lhs (buffer)
    read (buffer, *)  object%n_extra

    read (unit, "(A)")  buffer
    if (buffer(1:9) == " Resonant") then
       object%requires_resonances = .true.
       read (unit, "(A)")  buffer
    else
       object%requires_resonances = .false.
    end if

    call strip_equation_lhs (buffer)
    read (buffer(3:34), "(A32)")  object%md5sum

    if (object%n_initial > 0) then
       allocate (object%initial (object%n_initial))
       do i = 1, object%n_initial
          read (unit, "(A)")  buffer
          write (ref, "(1x,A,I0)")  "Component #", i
          if (buffer /= ref)  return                ! Wrong component header
          call object%initial(i)%read (unit, core_def_templates)
       end do
    end if

  end subroutine process_def_read

  module subroutine process_def_show (object, unit)
    class(process_def_t), intent(in) :: object
    integer, intent(in) :: unit
    integer :: i
    write (unit, "(4x,A)", advance="no") char (object%id)
    if (object%num_id /= 0) &
         write (unit, "(1x,'(',I0,')')", advance="no")  object%num_id
    if (object%model_name /= "") &
         write (unit, "(1x,'[',A,']')", advance="no")  char (object%model_name)
    if (object%requires_resonances) then
       write (unit, "(1x,A)", advance="no")  "[+ resonant subprocesses]"
    end if
    write (unit, *)
    if (allocated (object%initial)) then
       do i = 1, size (object%initial)
          call object%initial(i)%show (unit)
       end do
    end if
    if (allocated (object%extra)) then
       do i = 1, size (object%extra)
          call object%extra(i)%show (unit)
       end do
    end if
  end subroutine process_def_show

  module subroutine process_def_init (def, id, &
       model, model_name, n_in, n_components, num_id, &
       nlo_process, negative_sf, requires_resonances)
    class(process_def_t), intent(out) :: def
    type(string_t), intent(in), optional :: id
    class(model_data_t), intent(in), optional, target :: model
    type(string_t), intent(in), optional :: model_name
    integer, intent(in), optional :: n_in
    integer, intent(in), optional :: n_components
    integer, intent(in), optional :: num_id
    logical, intent(in), optional :: nlo_process
    logical, intent(in), optional :: negative_sf
    logical, intent(in), optional :: requires_resonances
    character(16) :: suffix
    integer :: i
    if (present (id)) then
       def%id = id
    else
       def%id = ""
    end if
    if (present (num_id)) then
       def%num_id = num_id
    end if
    if (present (model)) then
       def%model => model
       def%model_name = model%get_name ()
    else
       def%model => null ()
       if (present (model_name)) then
          def%model_name = model_name
       else
          def%model_name = ""
       end if
    end if
    if (present (n_in))  def%n_in = n_in
    if (present (n_components)) then
       def%n_initial = n_components
       allocate (def%initial (n_components))
    end if
    if (present (nlo_process)) then
       def%nlo_process = nlo_process
    end if
    if (present (negative_sf)) then
       def%negative_sf = negative_sf
    end if
    if (present (requires_resonances)) then
       def%requires_resonances = requires_resonances
    end if
    def%initial%initial = .true.
    def%initial%method     = ""
    do i = 1, def%n_initial
       write (suffix, "(A,I0)")  "_i", i
       def%initial(i)%basename = def%id // trim (suffix)
    end do
    def%initial%description = ""
  end subroutine process_def_init

  module subroutine process_def_set_model_name (def, model_name)
    class(process_def_t), intent(inout) :: def
    type(string_t), intent(in) :: model_name
    def%model_name = model_name
  end subroutine process_def_set_model_name

  module subroutine process_def_import_component (def, &
       i, n_out, prt_in, prt_out, method, variant, &
       nlo_type, can_be_integrated)
    class(process_def_t), intent(inout) :: def
    integer, intent(in) :: i
    integer, intent(in), optional :: n_out
    type(prt_spec_t), dimension(:), intent(in), optional :: prt_in
    type(prt_spec_t), dimension(:), intent(in), optional :: prt_out
    type(string_t), intent(in), optional :: method
    integer, intent(in), optional :: nlo_type
    logical, intent(in), optional :: can_be_integrated
    type(string_t) :: nlo_type_string
    class(prc_core_def_t), &
         intent(inout), allocatable, optional :: variant
    integer :: p
    associate (comp => def%initial(i))
      if (present (n_out)) then
         comp%n_in  = def%n_in
         comp%n_out = n_out
         comp%n_tot = def%n_in + n_out
      end if
      if (present (prt_in)) then
         allocate (comp%prt_in (size (prt_in)))
         comp%prt_in = prt_in
      end if
      if (present (prt_out)) then
         allocate (comp%prt_out (size (prt_out)))
         comp%prt_out = prt_out
      end if
      if (present (method))  comp%method = method
      if (present (variant)) then
         call move_alloc (variant, comp%core_def)
      end if
      if (present (nlo_type)) then
        comp%nlo_type = nlo_type
      end if
      if (present (can_be_integrated)) then
         comp%active = can_be_integrated
      else
         comp%active = .true.
      end if
      if (allocated (comp%prt_in) .and. allocated (comp%prt_out)) then
         associate (d => comp%description)
           d = ""
           do p = 1, size (prt_in)
              if (p > 1)  d = d // ", "
              d = d // comp%prt_in(p)%to_string ()
           end do
           d = d // " => "
           do p = 1, size (prt_out)
              if (p > 1)  d = d // ", "
              d = d // comp%prt_out(p)%to_string ()
           end do
           if (comp%method /= "") then
              if ((def%nlo_process .and. .not. comp%active) .or. &
                   comp%nlo_type == NLO_SUBTRACTION) then
                 d = d // " [inactive]"
              else
                 d = d // " [" // comp%method // "]"
              end if
           end if
           nlo_type_string = component_status (comp%nlo_type)
           if (nlo_type_string /= "born") then
             d = d // ", [" // nlo_type_string // "]"
           end if
         end associate
      end if
    end associate
  end subroutine process_def_import_component

  module function process_def_get_n_components (def) result (n)
    class(process_def_t), intent(in) :: def
    integer :: n
    n = size (def%initial)
  end function process_def_get_n_components

  module subroutine process_def_set_fixed_emitter (def, i, emitter)
    class(process_def_t), intent(inout) :: def
    integer, intent(in) :: i, emitter
    def%initial(i)%fixed_emitter = emitter
  end subroutine process_def_set_fixed_emitter

  module subroutine process_def_set_coupling_powers (def, alpha_power, alphas_power)
    class(process_def_t), intent(inout) :: def
    integer, intent(in) :: alpha_power, alphas_power
    def%initial(1)%alpha_power = alpha_power
    def%initial(1)%alphas_power = alphas_power
  end subroutine process_def_set_coupling_powers

  module subroutine process_def_set_associated_components (def, i, &
       i_list, remnant, real_finite, mismatch)
    class(process_def_t), intent(inout) :: def
    logical, intent(in) :: remnant, real_finite, mismatch
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: i_list
    integer :: add_index
    add_index = 0
    associate (comp => def%initial(i)%associated_components)
       comp(ASSOCIATED_BORN) = i_list(1)
       comp(ASSOCIATED_REAL) = i_list(2)
       comp(ASSOCIATED_VIRT) = i_list(3)
       comp(ASSOCIATED_SUB) = i_list(4)
       if (remnant) then
          comp(ASSOCIATED_PDF) = i_list(5)
          add_index = add_index + 1
       end if
       if (real_finite) then
          comp(ASSOCIATED_REAL_FIN) = i_list(5+add_index)
          add_index = add_index + 1
       end if
       if (mismatch) then
          !!! incomplete
       end if
    end associate
  end subroutine process_def_set_associated_components

  module subroutine process_def_compute_md5sum (def, model)
    class(process_def_t), intent(inout) :: def
    class(model_data_t), intent(in), optional, target :: model
    integer :: i
    type(string_t) :: buffer
    buffer = def%model_name
    do i = 1, def%n_initial
       call def%initial(i)%compute_md5sum (model)
       buffer = buffer // def%initial(i)%md5sum
    end do
    do i = 1, def%n_extra
       call def%extra(i)%compute_md5sum (model)
       buffer = buffer // def%initial(i)%md5sum
    end do
    def%md5sum = md5sum (char (buffer))
  end subroutine process_def_compute_md5sum

  module function process_def_get_md5sum (def, i_component) result (md5sum)
    class(process_def_t), intent(in) :: def
    integer, intent(in), optional :: i_component
    character(32) :: md5sum
    if (present (i_component)) then
       md5sum = def%initial(i_component)%md5sum
    else
       md5sum = def%md5sum
    end if
  end function process_def_get_md5sum

  module function process_def_get_core_def_ptr (def, i_component) result (ptr)
    class(process_def_t), intent(in), target :: def
    integer, intent(in) :: i_component
    class(prc_core_def_t), pointer :: ptr
    ptr => def%initial(i_component)%get_core_def_ptr ()
  end function process_def_get_core_def_ptr

  module function process_def_needs_code (def, i_component) result (flag)
    class(process_def_t), intent(in) :: def
    integer, intent(in) :: i_component
    logical :: flag
    flag = def%initial(i_component)%needs_code ()
  end function process_def_needs_code

  module subroutine process_def_get_pdg_in_1 (def, pdg)
    class(process_def_t), intent(in), target :: def
    integer, dimension(:), intent(out) :: pdg
    call def%initial(1)%get_pdg_in (def%model, pdg)
  end subroutine process_def_get_pdg_in_1

  elemental module function process_def_is_nlo (def) result (flag)
    logical :: flag
    class(process_def_t), intent(in) :: def
    flag = def%nlo_process
  end function process_def_is_nlo

  elemental module function process_def_get_nlo_type (def, i_component) result (nlo_type)
    integer :: nlo_type
    class(process_def_t), intent(in) :: def
    integer, intent(in) :: i_component
    nlo_type = def%initial(i_component)%nlo_type
  end function process_def_get_nlo_type

  elemental module function process_def_get_negative_sf (def) result (neg_sf)
    logical :: neg_sf
    class(process_def_t), intent(in) :: def
    neg_sf = def%negative_sf
  end function process_def_get_negative_sf

  module function process_def_get_n_in (def) result (n_in)
    class(process_def_t), intent(in) :: def
    integer :: n_in
    n_in = def%n_in
  end function process_def_get_n_in

  module function process_def_get_component_def_ptr (def, i) result (component)
    type(process_component_def_t), pointer :: component
    class(process_def_t), intent(in), target :: def
    integer, intent(in) :: i
    if (i <= def%n_initial) then
       component => def%initial(i)
    else
       component => null ()
    end if
  end function process_def_get_component_def_ptr

  module subroutine process_def_list_final (list)
    class(process_def_list_t), intent(inout) :: list
    type(process_def_entry_t), pointer :: current
    nullify (list%last)
    do while (associated (list%first))
       current => list%first
       list%first => current%next
       deallocate (current)
    end do
  end subroutine process_def_list_final

  module subroutine process_def_list_write (object, unit, libpath)
    class(process_def_list_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: libpath
    type(process_def_entry_t), pointer :: entry
    integer :: i, u
    u = given_output_unit (unit)
    if (associated (object%first)) then
       i = 1
       entry => object%first
       do while (associated (entry))
          write (u, "(1x,A,I0,A)")  "Process #", i, ":"
          call entry%write (u)
          i = i + 1
          entry => entry%next
          if (associated (entry))  write (u, *)
       end do
    else
       write (u, "(1x,A)")  "Process definition list: [empty]"
    end if
  end subroutine process_def_list_write

  module subroutine process_def_list_show (object, unit)
    class(process_def_list_t), intent(in) :: object
    integer, intent(in), optional :: unit
    type(process_def_entry_t), pointer :: entry
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%first)) then
       write (u, "(2x,A)")  "Processes:"
       entry => object%first
       do while (associated (entry))
          call entry%show (u)
          entry => entry%next
       end do
    else
       write (u, "(2x,A)")  "Processes: [empty]"
    end if
  end subroutine process_def_list_show

  module subroutine process_def_list_read (object, unit, core_def_templates)
    class(process_def_list_t), intent(out) :: object
    integer, intent(in) :: unit
    type(prc_template_t), dimension(:), intent(in) :: core_def_templates
    type(process_def_entry_t), pointer :: entry
    character(80) :: buffer, ref
    integer :: i
    read (unit, "(A)")  buffer
    write (ref, "(1x,A)")  "Process definition list: [empty]"
    if (buffer == ref)  return         ! OK: empty library
    backspace (unit)
    READ_ENTRIES: do i = 1, huge(0)-1
       if (i > 1) read (unit, *, end=1)
       read (unit, "(A)")  buffer

       write (ref, "(1x,A,I0,A)")  "Process #", i, ":"
       if (buffer /= ref)  return      ! Wrong process header: done.
       allocate (entry)
       call entry%read (unit, core_def_templates)
       call object%append (entry)
    end do READ_ENTRIES
1   continue                           ! EOF: done
  end subroutine process_def_list_read

  module subroutine process_def_list_append (list, entry)
    class(process_def_list_t), intent(inout) :: list
    type(process_def_entry_t), intent(inout), pointer :: entry
    if (list%contains (entry%id)) then
       call msg_fatal ("Recording process: '" // char (entry%id) &
            // "' has already been defined")
    end if
    if (associated (list%first)) then
       list%last%next => entry
    else
       list%first => entry
    end if
    list%last => entry
    entry => null ()
  end subroutine process_def_list_append

  module function process_def_list_get_n_processes (list) result (n)
    integer :: n
    class(process_def_list_t), intent(in) :: list
    type(process_def_entry_t), pointer :: current
    n = 0
    current => list%first
    do while (associated (current))
       n = n + 1
       current => current%next
    end do
  end function process_def_list_get_n_processes

  module subroutine process_def_list_get_process_id_list (list, id)
    class(process_def_list_t), intent(in) :: list
    type(string_t), dimension(:), allocatable, intent(out) :: id
    type(process_def_entry_t), pointer :: current
    integer :: i
    allocate (id (list%get_n_processes ()))
    i = 0
    current => list%first
    do while (associated (current))
       i = i + 1
       id(i) = current%id
       current => current%next
    end do
  end subroutine process_def_list_get_process_id_list

  module subroutine process_def_list_get_process_id_req_resonant (list, id)
    class(process_def_list_t), intent(in) :: list
    type(string_t), dimension(:), allocatable, intent(out) :: id
    type(process_def_entry_t), pointer :: current
    integer :: i
    allocate (id (list%get_n_processes ()))
    i = 0
    current => list%first
    do while (associated (current))
       if (current%requires_resonances) then
          i = i + 1
          id(i) = current%id
       end if
       current => current%next
    end do
    id = id(1:i)
  end subroutine process_def_list_get_process_id_req_resonant

  module function process_def_list_get_process_def_ptr (list, id) result (entry)
    type(process_def_entry_t), pointer :: entry
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(process_def_entry_t), pointer :: current
    current => list%first
    do while (associated (current))
       if (id == current%id)  exit
       current => current%next
    end do
    entry => current
  end function process_def_list_get_process_def_ptr

  module function process_def_list_contains (list, id) result (flag)
    logical :: flag
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(process_def_entry_t), pointer :: current
    current => list%get_process_def_ptr (id)
    flag = associated (current)
  end function process_def_list_contains

  module function process_def_list_get_entry_index (list, id) result (n)
    integer :: n
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(process_def_entry_t), pointer :: current
    n = 0
    current => list%first
    do while (associated (current))
       n = n + 1
       if (id == current%id) then
          return
       end if
       current => current%next
    end do
    n = 0
  end function process_def_list_get_entry_index

  module function process_def_list_get_num_id (list, id) result (num_id)
    integer :: num_id
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(process_def_entry_t), pointer :: current
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       num_id = current%num_id
    else
       num_id = 0
    end if
  end function process_def_list_get_num_id

  module function process_def_list_get_model_name (list, id) result (model_name)
    type(string_t) :: model_name
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(process_def_entry_t), pointer :: current
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       model_name = current%model_name
    else
       model_name = ""
    end if
  end function process_def_list_get_model_name

  module function process_def_list_get_n_in (list, id) result (n)
    integer :: n
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(process_def_entry_t), pointer :: current
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       n = current%n_in
    else
       n = 0
    end if
  end function process_def_list_get_n_in

  module subroutine process_def_list_get_pdg_in_1 (list, id, pdg)
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    integer, dimension(:), intent(out) :: pdg
    type(process_def_entry_t), pointer :: current
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       call current%get_pdg_in_1 (pdg)
    else
       pdg = 0
    end if
  end subroutine process_def_list_get_pdg_in_1

  module subroutine process_def_list_get_component_list (list, id, cid)
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(string_t), dimension(:), allocatable, intent(out) :: cid
    type(process_def_entry_t), pointer :: current
    integer :: i, n
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       allocate (cid (current%n_initial + current%n_extra))
       do i = 1, current%n_initial
          cid(i) = current%initial(i)%basename
       end do
       n = current%n_initial
       do i = 1, current%n_extra
          cid(n + i) = current%extra(i)%basename
       end do
    end if
  end subroutine process_def_list_get_component_list

  module subroutine process_def_list_get_component_description_list &
       (list, id, description)
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    type(string_t), dimension(:), allocatable, intent(out) :: description
    type(process_def_entry_t), pointer :: current
    integer :: i, n
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       allocate (description (current%n_initial + current%n_extra))
       do i = 1, current%n_initial
          description(i) = current%initial(i)%description
       end do
       n = current%n_initial
       do i = 1, current%n_extra
          description(n + i) = current%extra(i)%description
       end do
    end if
  end subroutine process_def_list_get_component_description_list

  module function process_def_list_req_resonant (list, id) result (flag)
    class(process_def_list_t), intent(in) :: list
    type(string_t), intent(in) :: id
    logical :: flag
    type(process_def_entry_t), pointer :: current
    current => list%get_process_def_ptr (id)
    if (associated (current)) then
       flag = current%requires_resonances
    else
       flag = .false.
    end if
  end function process_def_list_req_resonant

  module function process_library_entry_to_string (object) result (string)
    type(string_t) :: string
    class(process_library_entry_t), intent(in) :: object
    character(32) :: buffer
    string = "[" // STATUS_LETTER(object%status) // "]"
    select case (object%status)
    case (STAT_UNKNOWN)
    case default
       if (associated (object%def)) then
          write (buffer, "(I0)")  object%i_component
          string = string // " " // object%def%id // "." // trim (buffer)
       end if
       if (object%i_external /= 0) then
          write (buffer, "(I0)")  object%i_external
          string = string // " = ext:" // trim (buffer)
       else
          string = string // " = int"
       end if
       if (allocated (object%driver)) then
          string = string // " (" // object%driver%type_name () // ")"
       end if
    end select
  end function process_library_entry_to_string

  module subroutine process_library_entry_init (object, &
       status, def, i_component, i_external, driver_template)
    class(process_library_entry_t), intent(out) :: object
    integer, intent(in) :: status
    type(process_def_t), target, intent(in) :: def
    integer, intent(in) :: i_component
    integer, intent(in) :: i_external
    class(prc_core_driver_t), intent(inout), allocatable, optional &
         :: driver_template
    object%status = status
    object%def => def
    object%i_component = i_component
    object%i_external = i_external
    if (present (driver_template)) then
       call move_alloc (driver_template, object%driver)
    end if
  end subroutine process_library_entry_init

  module subroutine process_library_entry_connect (entry, lib_driver, i)
    class(process_library_entry_t), intent(inout) :: entry
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    call entry%def%initial(entry%i_component)%connect &
         (lib_driver, i, entry%driver)
  end subroutine process_library_entry_connect

  module subroutine process_library_write (object, unit, libpath)
    class(process_library_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: libpath
    integer :: i, u
    u = given_output_unit (unit)
    write (u, "(1x,A,A)")  "Process library: ", char (object%basename)
    write (u, "(3x,A,L1)")   "external        = ", object%external
    write (u, "(3x,A,L1)")   "makefile exists = ", object%makefile_exists
    write (u, "(3x,A,L1)")   "driver exists   = ", object%driver_exists
    write (u, "(3x,A,A1)")   "code status     = ", &
         STATUS_LETTER (object%status)
    write (u, *)
    if (allocated (object%entry)) then
       write (u, "(1x,A)", advance="no")  "Process library entries:"
       write (u, "(1x,I0)")  object%n_entries
       do i = 1, size (object%entry)
          write (u, "(1x,A,I0,A,A)")  "Entry #", i, ": ", &
               char (object%entry(i)%to_string ())
       end do
       write (u, *)
    end if
    if (object%external) then
       call object%driver%write (u, libpath)
       write (u, *)
    end if
    call object%process_def_list_t%write (u)
  end subroutine process_library_write

  module subroutine process_library_show (object, unit)
    class(process_library_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(A,A)")  "Process library: ", char (object%basename)
    write (u, "(2x,A,L1)")   "external        = ", object%external
    if (object%static) then
       write (u, "(2x,A,L1)")   "static          = ", .true.
    else
       write (u, "(2x,A,L1)")   "makefile exists = ", object%makefile_exists
       write (u, "(2x,A,L1)")   "driver exists   = ", object%driver_exists
    end if
    write (u, "(2x,A,A1)", advance="no")   "code status     = "
    select case (object%status)
    case (STAT_UNKNOWN);    write (u, "(A)")  "[unknown]"
    case (STAT_OPEN);       write (u, "(A)")  "open"
    case (STAT_CONFIGURED); write (u, "(A)")  "configured"
    case (STAT_SOURCE);     write (u, "(A)")  "source code exists"
    case (STAT_COMPILED);   write (u, "(A)")  "compiled"
    case (STAT_LINKED);     write (u, "(A)")  "linked"
    case (STAT_ACTIVE);     write (u, "(A)")  "active"
    end select
    call object%process_def_list_t%show (u)
  end subroutine process_library_show

  module subroutine process_library_init (lib, basename)
    class(process_library_t), intent(out) :: lib
    type(string_t), intent(in) :: basename
    lib%basename = basename
    lib%status = STAT_OPEN
    call msg_message ("Process library '" // char (basename) &
         // "': initialized")
  end subroutine process_library_init

  module subroutine process_library_init_static (lib, basename)
    class(process_library_t), intent(out) :: lib
    type(string_t), intent(in) :: basename
    lib%basename = basename
    lib%status = STAT_OPEN
    lib%static = .true.
    call msg_message ("Static process library '" // char (basename) &
         // "': initialized")
  end subroutine process_library_init_static

  module subroutine process_library_configure (lib, os_data)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    type(process_def_entry_t), pointer :: def_entry
    integer :: n_entries, n_external, i_entry, i_external
    type(string_t) :: model_name
    integer :: i_component

    n_entries = 0
    n_external = 0
    if (allocated (lib%entry))  deallocate (lib%entry)

    def_entry => lib%first
    do while (associated (def_entry))
       do i_component = 1, def_entry%n_initial
          n_entries = n_entries + 1
          if (def_entry%initial(i_component)%needs_code ()) then
             n_external = n_external + 1
             lib%external = .true.
          end if
       end do
       def_entry => def_entry%next
    end do

    call lib%allocate_entries (n_entries)

    i_entry = 0
    i_external = 0
    def_entry => lib%first
    do while (associated (def_entry))
       do i_component = 1, def_entry%n_initial
          i_entry = i_entry + 1
          associate (lib_entry => lib%entry(i_entry))
            lib_entry%status = STAT_CONFIGURED
            lib_entry%def => def_entry%process_def_t
            lib_entry%i_component = i_component
            if (def_entry%initial(i_component)%needs_code ()) then
               i_external = i_external + 1
               lib_entry%i_external = i_external
            end if
            call def_entry%initial(i_component)%allocate_driver &
                 (lib_entry%driver)
          end associate
       end do
       def_entry => def_entry%next
    end do

    call dispatch_prclib_driver (lib%driver, &
         lib%basename, lib%get_modellibs_ldflags (os_data))
    call lib%driver%init (n_external)
    do i_entry = 1, n_entries
       associate (lib_entry => lib%entry(i_entry))
         i_component = lib_entry%i_component
         model_name = lib_entry%def%model_name
         associate (def => lib_entry%def%initial(i_component))
           if (def%needs_code ()) then
              call lib%driver%set_record (lib_entry%i_external, &
                   def%basename, &
                   model_name, &
                   def%get_features (), def%get_writer_ptr ())
           end if
         end associate
       end associate
    end do

    if (lib%static) then
       if (lib%n_entries /= 0)  lib%entry%status = STAT_LINKED
       lib%status = STAT_LINKED
    else if (lib%external) then
       where (lib%entry%i_external == 0)  lib%entry%status = STAT_LINKED
       lib%status = STAT_CONFIGURED
       lib%makefile_exists = .false.
       lib%driver_exists = .false.
    else
       if (lib%n_entries /= 0)  lib%entry%status = STAT_LINKED
       lib%status = STAT_LINKED
    end if
  end subroutine process_library_configure

  module subroutine process_library_allocate_entries (lib, n_entries)
    class(process_library_t), intent(inout) :: lib
    integer, intent(in) :: n_entries
    lib%n_entries = n_entries
    allocate (lib%entry (n_entries))
  end subroutine process_library_allocate_entries

  module subroutine process_library_init_entry (lib, i, &
       status, def, i_component, i_external, driver_template)
    class(process_library_t), intent(inout) :: lib
    integer, intent(in) :: i
    integer, intent(in) :: status
    type(process_def_t), target, intent(in) :: def
    integer, intent(in) :: i_component
    integer, intent(in) :: i_external
    class(prc_core_driver_t), intent(inout), allocatable, optional &
         :: driver_template
    call lib%entry(i)%init (status, def, i_component, i_external, &
         driver_template)
  end subroutine process_library_init_entry

  module subroutine process_library_compute_md5sum (lib, model)
    class(process_library_t), intent(inout) :: lib
    class(model_data_t), intent(in), optional, target :: model
    type(process_def_entry_t), pointer :: def_entry
    type(string_t) :: buffer
    buffer = lib%basename
    def_entry => lib%first
    do while (associated (def_entry))
       call def_entry%compute_md5sum (model)
       buffer = buffer // def_entry%md5sum
       def_entry => def_entry%next
    end do
    lib%md5sum = md5sum (char (buffer))
    call lib%driver%set_md5sum (lib%md5sum)
  end subroutine process_library_compute_md5sum

  module subroutine process_library_write_makefile &
       (lib, os_data, force, verbose, testflag, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: force, verbose
    logical, intent(in), optional :: testflag
    type(string_t), intent(in), optional :: workspace
    character(32) :: md5sum_file
    logical :: generate
    integer :: unit
    if (lib%external .and. .not. lib%static) then
       generate = .true.
       if (.not. force) then
          md5sum_file = lib%driver%get_md5sum_makefile (workspace)
          if (lib%md5sum == md5sum_file) then
             call msg_message ("Process library '" // char (lib%basename) &
                  // "': keeping makefile")
             generate = .false.
          end if
       end if
       if (generate) then
          call msg_message ("Process library '" // char (lib%basename) &
               // "': writing makefile")
          unit = free_unit ()
          open (unit, &
               file = char (workspace_prefix (workspace) &
               &            // lib%driver%basename // ".makefile"), &
               status="replace", action="write")
          call lib%driver%generate_makefile (unit, os_data, verbose, testflag)
          close (unit)
       end if
       lib%makefile_exists = .true.
    end if
  end subroutine process_library_write_makefile

  module subroutine process_library_write_driver (lib, force, workspace)
    class(process_library_t), intent(inout) :: lib
    logical, intent(in) :: force
    type(string_t), intent(in), optional :: workspace
    character(32) :: md5sum_file
    logical :: generate
    integer :: unit
    if (lib%external .and. .not. lib%static) then
       generate = .true.
       if (.not. force) then
          md5sum_file = lib%driver%get_md5sum_driver (workspace)
          if (lib%md5sum == md5sum_file) then
             call msg_message ("Process library '" // char (lib%basename) &
                  // "': keeping driver")
             generate = .false.
          end if
       end if
       if (generate) then
          call msg_message ("Process library '" // char (lib%basename) &
               // "': writing driver")
          unit = free_unit ()
          open (unit, &
               file = char (workspace_prefix (workspace) &
               &            // lib%driver%basename // ".f90"), &
               status="replace", action="write")
          call lib%driver%generate_driver_code (unit)
          close (unit)
       end if
       lib%driver_exists = .true.
    end if
  end subroutine process_library_write_driver

  module subroutine process_library_update_status (lib, os_data, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in), optional :: workspace
    character(32) :: md5sum_file
    integer :: i, i_external, i_component
    if (lib%external) then
       select case (lib%status)
       case (STAT_CONFIGURED:STAT_LINKED)
          call lib%driver%load (os_data, noerror=.true., workspace=workspace)
       end select
       if (lib%driver%loaded) then
          md5sum_file = lib%driver%get_md5sum (0)
          if (lib%md5sum == md5sum_file) then
             call lib%load_entries ()
             lib%entry%status = STAT_ACTIVE
             lib%status = STAT_ACTIVE
             call msg_message ("Process library '" // char (lib%basename) &
                  // "': active")
          else
             do i = 1, lib%n_entries
                associate (entry => lib%entry(i))
                  i_external = entry%i_external
                  i_component = entry%i_component
                  if (i_external /= 0) then
                     md5sum_file = lib%driver%get_md5sum (i_external)
                     if (entry%def%get_md5sum (i_component) == md5sum_file) then
                        entry%status = STAT_COMPILED
                     else
                        entry%status = STAT_CONFIGURED
                     end if
                  end if
                end associate
             end do
             call lib%driver%unload ()
             lib%status = STAT_CONFIGURED
          end if
       end if
       select case (lib%status)
       case (STAT_CONFIGURED)
          do i = 1, lib%n_entries
             associate (entry => lib%entry(i))
               i_external = entry%i_external
               i_component = entry%i_component
               if (i_external /= 0) then
                  select case (entry%status)
                  case (STAT_CONFIGURED)
                     md5sum_file = lib%driver%get_md5sum_source &
                          (i_external, workspace)
                     if (entry%def%get_md5sum (i_component) == md5sum_file) then
                        entry%status = STAT_SOURCE
                     end if
                  end select
               end if
             end associate
          end do
          if (all (lib%entry%status >= STAT_SOURCE)) then
             md5sum_file = lib%driver%get_md5sum_driver (workspace)
             if (lib%md5sum == md5sum_file) then
                lib%status = STAT_SOURCE
             end if
          end if
       end select
    end if
  end subroutine process_library_update_status

  module subroutine process_library_make_source &
       (lib, os_data, keep_old_source, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    logical, intent(in), optional :: keep_old_source
    type(string_t), intent(in), optional :: workspace
    logical :: keep_old
    integer :: i, i_external
    keep_old = .false.
    if (present (keep_old_source))  keep_old = keep_old_source
    if (lib%external .and. .not. lib%static) then
       select case (lib%status)
       case (STAT_CONFIGURED)
          if (keep_old) then
             call msg_message ("Process library '" // char (lib%basename) &
                  // "': keeping source code")
          else
             call msg_message ("Process library '" // char (lib%basename) &
                  // "': creating source code")
             do i = 1, size (lib%entry)
                associate (entry => lib%entry(i))
                  i_external = entry%i_external
                  if (i_external /= 0 &
                       .and. lib%entry(i)%status == STAT_CONFIGURED) then
                     call lib%driver%clean_proc &
                          (i_external, os_data, workspace)
                  end if
                end associate
                if (signal_is_pending ())  return
             end do
             call lib%driver%make_source (os_data, workspace)
          end if
          lib%status = STAT_SOURCE
          where (lib%entry%i_external /= 0 &
               .and. lib%entry%status == STAT_CONFIGURED)
             lib%entry%status = STAT_SOURCE
          end where
          lib%status = STAT_SOURCE
       end select
    end if
  end subroutine process_library_make_source

  module subroutine process_library_make_compile &
       (lib, os_data, keep_old_source, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    logical, intent(in), optional :: keep_old_source
    type(string_t), intent(in), optional :: workspace
    if (lib%external .and. .not. lib%static) then
       select case (lib%status)
       case (STAT_CONFIGURED)
          call lib%make_source (os_data, keep_old_source, workspace)
       end select
       if (signal_is_pending ())  return
       select case (lib%status)
       case (STAT_SOURCE)
          call msg_message ("Process library '" // char (lib%basename) &
               // "': compiling sources")
          call lib%driver%make_compile (os_data, workspace)
          where (lib%entry%i_external /= 0 &
               .and. lib%entry%status == STAT_SOURCE)
             lib%entry%status = STAT_COMPILED
          end where
          lib%status = STAT_COMPILED
       end select
    end if
  end subroutine process_library_make_compile

  module subroutine process_library_make_link &
       (lib, os_data, keep_old_source, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    logical, intent(in), optional :: keep_old_source
    type(string_t), intent(in), optional :: workspace
    if (lib%external .and. .not. lib%static) then
       select case (lib%status)
       case (STAT_CONFIGURED:STAT_SOURCE)
          call lib%make_compile (os_data, keep_old_source, workspace)
       end select
       if (signal_is_pending ())  return
       select case (lib%status)
       case (STAT_COMPILED)
          call msg_message ("Process library '" // char (lib%basename) &
               // "': linking")
          call lib%driver%make_link (os_data, workspace)
          lib%entry%status = STAT_LINKED
          lib%status = STAT_LINKED
       end select
    end if
  end subroutine process_library_make_link

  module subroutine process_library_load (lib, os_data, keep_old_source, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    logical, intent(in), optional :: keep_old_source
    type(string_t), intent(in), optional :: workspace
    select case (lib%status)
    case (STAT_CONFIGURED:STAT_COMPILED)
       call lib%make_link (os_data, keep_old_source, workspace)
    end select
    if (signal_is_pending ())  return
    select case (lib%status)
    case (STAT_LINKED)
       if (lib%external) then
          call msg_message ("Process library '" // char (lib%basename) &
               // "': loading")
          call lib%driver%load (os_data, workspace=workspace)
          call lib%load_entries ()
       end if
       lib%entry%status = STAT_ACTIVE
       lib%status = STAT_ACTIVE
    end select
  end subroutine process_library_load

  module subroutine process_library_load_entries (lib)
    class(process_library_t), intent(inout) :: lib
    integer :: i
    do i = 1, size (lib%entry)
       associate (entry => lib%entry(i))
         if (entry%i_external /= 0) then
            call entry%connect (lib%driver, entry%i_external)
         end if
       end associate
    end do
  end subroutine process_library_load_entries

  module subroutine process_library_unload (lib)
    class(process_library_t), intent(inout) :: lib
    select case (lib%status)
    case (STAT_ACTIVE)
       if (lib%external) then
          call msg_message ("Process library '" // char (lib%basename) &
               // "': unloading")
          call lib%driver%unload ()
       end if
       lib%entry%status = STAT_LINKED
       lib%status = STAT_LINKED
    end select
  end subroutine process_library_unload

  module subroutine process_library_clean (lib, os_data, distclean, workspace)
    class(process_library_t), intent(inout) :: lib
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: distclean
    type(string_t), intent(in), optional :: workspace
    call lib%unload ()
    if (lib%external .and. .not. lib%static) then
       call msg_message ("Process library '" // char (lib%basename) &
            // "': removing old files")
       if (distclean) then
          call lib%driver%distclean (os_data, workspace)
       else
          call lib%driver%clean (os_data, workspace)
       end if
    end if
    where (lib%entry%i_external /= 0)
       lib%entry%status = STAT_CONFIGURED
    elsewhere
       lib%entry%status = STAT_LINKED
    end where
    if (lib%external) then
       lib%status = STAT_CONFIGURED
    else
       lib%status = STAT_LINKED
    end if
  end subroutine process_library_clean

  module subroutine process_library_open (lib)
    class(process_library_t), intent(inout) :: lib
    select case (lib%status)
    case (STAT_OPEN)
    case default
       call lib%unload ()
       if (.not. lib%static) then
          lib%entry%status = STAT_OPEN
          lib%status = STAT_OPEN
          if (lib%external)  lib%update_counter = lib%update_counter + 1
          call msg_message ("Process library '" // char (lib%basename) &
               // "': open")
       else
          call msg_error ("Static process library '" // char (lib%basename) &
               // "': processes can't be appended")
       end if
    end select
  end subroutine process_library_open

  module function process_library_get_name (lib) result (name)
    class(process_library_t), intent(in) :: lib
    type(string_t) :: name
    name = lib%basename
  end function process_library_get_name

  module function process_library_is_active (lib) result (flag)
    logical :: flag
    class(process_library_t), intent(in) :: lib
    flag = lib%status == STAT_ACTIVE
  end function process_library_is_active

  module function process_library_get_status (lib, i) result (status)
    class(process_library_t), intent(in) :: lib
    integer, intent(in), optional :: i
    integer :: status
    if (present (i)) then
       status = lib%entry(i)%status
    else
       status = lib%status
    end if
  end function process_library_get_status

  module function process_library_get_update_counter (lib) result (counter)
    class(process_library_t), intent(in) :: lib
    integer :: counter
    counter = lib%update_counter
  end function process_library_get_update_counter

  module subroutine process_library_set_status (lib, status, entries)
    class(process_library_t), intent(inout) :: lib
    integer, intent(in) :: status
    logical, intent(in), optional :: entries
    lib%status = status
    if (present (entries)) then
       if (entries)  lib%entry%status = status
    end if
  end subroutine process_library_set_status

  module function process_library_is_loaded (lib) result (flag)
    class(process_library_t), intent(in) :: lib
    logical :: flag
    flag = lib%driver%loaded
  end function process_library_is_loaded

  module subroutine process_library_entry_fill_constants (entry, driver, data)
    class(process_library_entry_t), intent(in) :: entry
    class(prclib_driver_t), intent(in) :: driver
    type(process_constants_t), intent(out) :: data
    integer :: i
    if (entry%i_external /= 0) then
       i = entry%i_external
       data%id         = driver%get_process_id (i)
       data%model_name = driver%get_model_name (i)
       data%md5sum     = driver%get_md5sum (i)
       data%openmp_supported = driver%get_openmp_status (i)
       data%n_in  = driver%get_n_in  (i)
       data%n_out = driver%get_n_out (i)
       data%n_flv = driver%get_n_flv (i)
       data%n_hel = driver%get_n_hel (i)
       data%n_col = driver%get_n_col (i)
       data%n_cin = driver%get_n_cin (i)
       data%n_cf  = driver%get_n_cf  (i)
       call driver%set_flv_state (i, data%flv_state)
       call driver%set_hel_state (i, data%hel_state)
       call driver%set_col_state (i, data%col_state, data%ghost_flag)
       call driver%set_color_factors (i, data%color_factors, data%cf_index)
    else
       select type (proc_driver => entry%driver)
       class is (process_driver_internal_t)
          call proc_driver%fill_constants (data)
       end select
    end if
  end subroutine process_library_entry_fill_constants

  module subroutine process_library_fill_constants (lib, id, i_component, data)
    class(process_library_t), intent(in) :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    type(process_constants_t), intent(out) :: data
    integer :: i
    do i = 1, size (lib%entry)
       associate (entry => lib%entry(i))
          if (entry%def%id == id .and. entry%i_component == i_component) then
             call entry%fill_constants (lib%driver, data)
             return
          end if
       end associate
    end do
  end subroutine process_library_fill_constants

  module subroutine process_library_connect_process &
       (lib, id, i_component, data, proc_driver)
    class(process_library_t), intent(in) :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    type(process_constants_t), intent(out) :: data
    class(prc_core_driver_t), allocatable, intent(out) :: proc_driver
    integer :: i
    do i = 1, size (lib%entry)
       associate (entry => lib%entry(i))
         if (entry%def%id == id .and. entry%i_component == i_component) then
            call entry%fill_constants (lib%driver, data)
            allocate (proc_driver, source = entry%driver)
            return
         end if
       end associate
    end do
    call msg_fatal ("Process library '" // char (lib%basename) &
               // "': process '" // char (id) // "' not found")
  end subroutine process_library_connect_process

  module subroutine process_library_test_transfer_md5sum (lib, r, e, c)
    class(process_library_t), intent(inout) :: lib
    integer, intent(in) :: r, e, c
    associate (writer => lib%driver%record(r)%writer)
       writer%md5sum = lib%entry(e)%def%get_md5sum (c)
    end associate
  end subroutine process_library_test_transfer_md5sum

  module function process_library_get_nlo_type (lib, id, i_component) result (nlo_type)
    integer :: nlo_type
    class(process_library_t), intent(in) :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    integer :: i
    do i = 1, size (lib%entry)
       if (lib%entry(i)%def%id == id .and. lib%entry(i)%i_component == i_component) then
          nlo_type = lib%entry(i)%def%get_nlo_type (i_component)
          exit
       end if
    end do
  end function process_library_get_nlo_type

  module function process_library_get_modellibs_ldflags (prc_lib, os_data) result (flags)
    class(process_library_t), intent(in) :: prc_lib
    type(os_data_t), intent(in) :: os_data
    type(string_t) :: flags
    type(string_t), dimension(:), allocatable :: models
    type(string_t) :: modelname, modellib, modellib_full
    logical :: exist
    integer :: i, j, mi
    flags = " -lomega"
    if ((.not. os_data%use_testfiles) .and. &
               os_dir_exist (os_data%whizard_models_libpath_local)) &
                 flags = flags // " -L" // os_data%whizard_models_libpath_local
    flags = flags // " -L" // os_data%whizard_models_libpath
    allocate (models(prc_lib%n_entries + 1))
    models = ""
    mi = 1
    if (allocated (prc_lib%entry)) then
       SCAN: do i = 1, prc_lib%n_entries
          if (associated (prc_lib%entry(i)%def)) then
             if (prc_lib%entry(i)%def%model_name /= "") then
                modelname = prc_lib%entry(i)%def%model_name
             else
                cycle SCAN
             end if
          else
             cycle SCAN
          end if
          do j = 1, mi
             if (models(mi) == modelname) cycle SCAN
          end do
          models(mi) = modelname
          mi = mi + 1
          if (os_data%use_libtool) then
             modellib = "libparameters_" // modelname // ".la"
          else
             modellib = "libparameters_" // modelname // ".a"
          end if
          exist = .false.
          if (.not. os_data%use_testfiles) then
             modellib_full = os_data%whizard_models_libpath_local &
                  // "/" // modellib
             inquire (file=char (modellib_full), exist=exist)
          end if
          if (.not. exist) then
             modellib_full = os_data%whizard_models_libpath &
                  // "/" // modellib
             inquire (file=char (modellib_full), exist=exist)
          end if
          if (exist) flags = flags // " -lparameters_" // modelname
       end do SCAN
    end if
    deallocate (models)
    flags = flags // " -lwhizard"
  end function process_library_get_modellibs_ldflags

  module function process_library_get_static_modelname (prc_lib, os_data) result (name)
    class(process_library_t), intent(in) :: prc_lib
    type(os_data_t), intent(in) :: os_data
    type(string_t) :: name
    type(string_t), dimension(:), allocatable :: models
    type(string_t) :: modelname, modellib, modellib_full
    logical :: exist
    integer :: i, j, mi
    name = ""
    allocate (models(prc_lib%n_entries + 1))
    models = ""
    mi = 1
    if (allocated (prc_lib%entry)) then
       SCAN: do i = 1, prc_lib%n_entries
          if (associated (prc_lib%entry(i)%def)) then
             if (prc_lib%entry(i)%def%model_name /= "") then
                modelname = prc_lib%entry(i)%def%model_name
             else
                cycle SCAN
             end if
          else
             cycle SCAN
          end if
          do j = 1, mi
             if (models(mi) == modelname) cycle SCAN
          end do
          models(mi) = modelname
          mi = mi + 1
          modellib = "libparameters_" // modelname // ".a"
          exist = .false.
          if (.not. os_data%use_testfiles) then
             modellib_full = os_data%whizard_models_libpath_local &
                  // "/" // modellib
             inquire (file=char (modellib_full), exist=exist)
          end if
          if (.not. exist) then
             modellib_full = os_data%whizard_models_libpath &
                  // "/" // modellib
             inquire (file=char (modellib_full), exist=exist)
          end if
          if (exist) name = name // " " // modellib_full
       end do SCAN
    end if
    deallocate (models)
  end function process_library_get_static_modelname


end submodule process_libraries_s

