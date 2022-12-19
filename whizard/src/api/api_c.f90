subroutine whizard_create (whizard_handle) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_loc  !NODEP!
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(inout) :: whizard_handle
  type(whizard_api_t), pointer :: whizard

  allocate (whizard)
  whizard_handle = c_loc (whizard)

end subroutine whizard_create

subroutine whizard_option (whizard_handle, key, value) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(inout) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: key
  character(c_char), dimension(*), intent(in) :: value

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%option (string_c2f (key), string_c2f (value))

end subroutine whizard_option

subroutine whizard_init (whizard_handle) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%init ()

end subroutine whizard_init

subroutine whizard_final (whizard_handle) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%final ()
  deallocate (whizard)

end subroutine whizard_final


subroutine whizard_set_double (whizard_handle, var, value) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  real(c_double), intent(in), value :: value

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%set_var (string_c2f (var), real (value, real64))

end subroutine whizard_set_double

subroutine whizard_set_int (whizard_handle, var, value) bind (C)
  use iso_fortran_env, only: int32  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  integer(c_int), intent(in), value :: value

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%set_var (string_c2f (var), int (value, int32))

end subroutine whizard_set_int

subroutine whizard_set_bool (whizard_handle, var, value) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  integer(c_int), intent(in), value :: value

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%set_var (string_c2f (var), value /= 0)

end subroutine whizard_set_bool

subroutine whizard_set_char (whizard_handle, var, value) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  character(c_char), dimension(*), intent(in) :: value

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%set_var (string_c2f (var), string_c2f (value))

end subroutine whizard_set_char

function whizard_get_double (whizard_handle, var, value) result (stat) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  real(c_double), intent(inout) :: value
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  logical :: known
  real(real64) :: v

  call c_f_pointer (whizard_handle, whizard)
  call whizard%get_var (string_c2f (var), v, known)
  if (known) then
     value = v
     stat = 0
  else
     value = 0
     stat = 1
  end if

end function whizard_get_double

function whizard_get_int (whizard_handle, var, value) result (stat) bind (C)
  use iso_fortran_env, only: int32  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  integer(c_int), intent(inout) :: value
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  logical :: known
  integer(int32) :: v

  call c_f_pointer (whizard_handle, whizard)
  call whizard%get_var (string_c2f (var), v, known)
  if (known) then
     value = v
     stat = 0
  else
     value = 0
     stat = 1
  end if

end function whizard_get_int

function whizard_get_bool (whizard_handle, var, value) result (stat) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  integer(c_int), intent(inout) :: value
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  logical :: known
  logical :: v

  call c_f_pointer (whizard_handle, whizard)
  call whizard%get_var (string_c2f (var), v, known)
  if (known) then
     if (v) then
        value = 1
     else
        value = 0
     end if
     stat = 0
  else
     value = 0
     stat = 1
  end if

end function whizard_get_bool

function whizard_get_char (whizard_handle, var, value, strlen) result (stat) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use string_utils, only: strcpy_f2c
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  character(c_char), dimension(*), intent(inout) :: value
  integer(c_int), value :: strlen
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  logical :: known
  character(:), allocatable :: v

  call c_f_pointer (whizard_handle, whizard)
  call whizard%get_var_character (string_c2f (var), v, known, strlen-1)
  if (known) then
     call strcpy_f2c (v, value)
     stat = 0
  else
     call strcpy_f2c ("", value)
     stat = 1
  end if

end function whizard_get_char

function whizard_get_char_len (whizard_handle, var) result (strlen) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: var
  integer(c_int) :: strlen

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  strlen = whizard%get_var_character_length (string_c2f (var)) + 1

end function whizard_get_char_len

function whizard_flv_string (whizard_handle, pdg, fstr, strlen) result (stat) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: strcpy_f2c
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  integer(c_int), value :: pdg
  character(c_char), dimension(*), intent(inout) :: fstr
  integer(c_int), value :: strlen
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  character(:), allocatable :: v

  call c_f_pointer (whizard_handle, whizard)
  v = whizard%flv_string (pdg)
  if (len (v) < strlen) then
     call strcpy_f2c (v, fstr)
     stat = 0
  else
     call strcpy_f2c ("", fstr)
     stat = 1
  end if

end function whizard_flv_string

function whizard_flv_string_len (whizard_handle, pdg) result (strlen) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  integer(c_int), value :: pdg
  integer(c_int) :: strlen

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  strlen = len (whizard%flv_string (pdg)) + 1

end function whizard_flv_string_len

function whizard_flv_array_string (whizard_handle, pdg, nf, fstr, strlen) result (stat) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: strcpy_f2c
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  integer(c_int), dimension(*), intent(in) :: pdg
  integer(c_int), value :: nf
  character(c_char), dimension(*), intent(inout) :: fstr
  integer(c_int), value :: strlen
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  character(:), allocatable :: v

  call c_f_pointer (whizard_handle, whizard)
  v = whizard%flv_string (pdg(1:nf))
  if (len (v) < strlen) then
     call strcpy_f2c (v, fstr)
     stat = 0
  else
     call strcpy_f2c ("", fstr)
     stat = 1
  end if

end function whizard_flv_array_string

function whizard_flv_array_string_len (whizard_handle, pdg, nf) result (strlen) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  integer(c_int), dimension(*), intent(in) :: pdg
  integer(c_int), value :: nf
  integer(c_int) :: strlen

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  strlen = len (whizard%flv_string (pdg(1:nf))) + 1

end function whizard_flv_array_string_len

subroutine whizard_command (whizard_handle, cmd) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: cmd

  type(whizard_api_t), pointer :: whizard

  call c_f_pointer (whizard_handle, whizard)
  call whizard%command (string_c2f (cmd))

end subroutine whizard_command

function whizard_get_integration_result (whizard_handle, proc_id, integral, error) result (stat) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: proc_id
  real(c_double), intent(inout) :: integral
  real(c_double), intent(inout) :: error
  integer(c_int) :: stat

  type(whizard_api_t), pointer :: whizard
  logical :: known
  real(real64) :: int, err

  call c_f_pointer (whizard_handle, whizard)
  call whizard%get_integration_result (string_c2f (proc_id), int, err, known)
  if (known) then
     integral = int
     error = err
     stat = 0
  else
     integral = 0
     error = 0
     stat = 1
  end if

end function whizard_get_integration_result

subroutine whizard_new_sample (whizard_handle, name, sample_handle) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use iso_c_binding, only: c_loc  !NODEP!
  use string_utils, only: string_c2f
  use api, only: whizard_api_t
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(in) :: whizard_handle
  character(c_char), dimension(*), intent(in) :: name
  type(c_ptr), intent(inout) :: sample_handle

  type(whizard_api_t), pointer :: whizard
  type(simulation_api_t), pointer :: sample

  call c_f_pointer (whizard_handle, whizard)
  allocate (sample)
  sample_handle = c_loc (sample)
  call whizard%new_sample (string_c2f (name), sample)

end subroutine whizard_new_sample

subroutine whizard_sample_open (sample_handle, it_begin, it_end) bind (C)
  use iso_fortran_env, only: int32  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  type(c_ptr), intent(inout) :: sample_handle
  integer(c_int), intent(inout) :: it_begin
  integer(c_int), intent(inout) :: it_end

  type(simulation_api_t), pointer :: sample
  integer(int32) :: it_b, it_e

  call c_f_pointer (sample_handle, sample)
  call sample%open (it_b, it_e)
  it_begin = it_b
  it_end = it_e

end subroutine whizard_sample_open

subroutine whizard_sample_next_event (sample_handle) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  type(c_ptr), intent(inout) :: sample_handle

  type(simulation_api_t), pointer :: sample

  call c_f_pointer (sample_handle, sample)
  call sample%next_event ()

end subroutine whizard_sample_next_event

subroutine whizard_sample_close (sample_handle) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  type(c_ptr), intent(inout) :: sample_handle

  type(simulation_api_t), pointer :: sample

  call c_f_pointer (sample_handle, sample)
  call sample%close ()
  deallocate (sample)

end subroutine whizard_sample_close

function whizard_sample_next_event_hepmc (sample_handle) result (evt) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use hepmc_interface, only: hepmc_event_t
  use hepmc_interface, only: hepmc_event_get_c_ptr
  use api, only: simulation_api_t

  type(c_ptr), intent(inout) :: sample_handle
  type(c_ptr) :: evt

  type(simulation_api_t), pointer :: sample
  type(hepmc_event_t) :: hepmc_event

  call c_f_pointer (sample_handle, sample)
  call sample%next_event (hepmc_event)
  evt = hepmc_event_get_c_ptr (hepmc_event)

end function whizard_sample_next_event_hepmc

function whizard_sample_next_event_lcio (sample_handle) result (evt) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use lcio_interface, only: lcio_event_t
  use lcio_interface, only: lcio_event_get_c_ptr
  use api, only: simulation_api_t

  type(c_ptr), intent(inout) :: sample_handle
  type(c_ptr) :: evt

  type(simulation_api_t), pointer :: sample
  type(lcio_event_t) :: lcio_event

  call c_f_pointer (sample_handle, sample)
  call sample%next_event (lcio_event)
  evt = lcio_event_get_c_ptr (lcio_event)

end function whizard_sample_next_event_lcio

subroutine whizard_sample_get_event_index (sample_handle, idx) bind (C)
  use iso_fortran_env, only: int32  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  integer(c_int), intent(inout) :: idx

  type(simulation_api_t), pointer :: sample
  integer(int32) :: i

  call c_f_pointer (sample_handle, sample)
  call sample%get_event_index (i)
  idx = i

end subroutine whizard_sample_get_event_index

subroutine whizard_sample_get_process_index (sample_handle, i_proc) bind (C)
  use iso_fortran_env, only: int32  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  integer(c_int), intent(inout) :: i_proc

  type(simulation_api_t), pointer :: sample
  integer(int32) :: i

  call c_f_pointer (sample_handle, sample)
  call sample%get_process_index (i)
  i_proc = i

end subroutine whizard_sample_get_process_index

subroutine whizard_sample_get_process_id (sample_handle, proc_id, strlen) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: strcpy_f2c
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  character(c_char), dimension(*), intent(inout) :: proc_id
  integer(c_int), value :: strlen

  type(simulation_api_t), pointer :: sample
  character(:), allocatable :: p_id

  call c_f_pointer (sample_handle, sample)
  call sample%get_process_id (p_id)
  if (len (p_id) < strlen) then
     call strcpy_f2c (p_id, proc_id)
  else
     call strcpy_f2c (p_id(1:strlen-1), proc_id)
  end if

end subroutine whizard_sample_get_process_id

function whizard_sample_get_process_id_len (sample_handle) result (strlen) bind (C)
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  integer(c_int) :: strlen

  type(simulation_api_t), pointer :: sample
  character(:), allocatable :: p_id

  call c_f_pointer (sample_handle, sample)
  call sample%get_process_id (p_id)
  strlen = len (p_id)

end function whizard_sample_get_process_id_len

subroutine whizard_sample_get_fac_scale (sample_handle, f_scale) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  real(c_double), intent(inout) :: f_scale

  type(simulation_api_t), pointer :: sample
  real(real64) :: s

  call c_f_pointer (sample_handle, sample)
  call sample%get_fac_scale (s)
  f_scale = s

end subroutine whizard_sample_get_fac_scale

subroutine whizard_sample_get_alpha_s (sample_handle, alpha_s) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  real(c_double), intent(inout) :: alpha_s

  type(simulation_api_t), pointer :: sample
  real(real64) :: a

  call c_f_pointer (sample_handle, sample)
  call sample%get_alpha_s (a)
  alpha_s = a

end subroutine whizard_sample_get_alpha_s

subroutine whizard_sample_get_weight (sample_handle, weight) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  real(c_double), intent(inout) :: weight

  type(simulation_api_t), pointer :: sample
  real(real64) :: w

  call c_f_pointer (sample_handle, sample)
  call sample%get_weight (w)
  weight = w

end subroutine whizard_sample_get_weight

subroutine whizard_sample_get_sqme (sample_handle, sqme) bind (C)
  use iso_fortran_env, only: real64  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_double  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use api, only: simulation_api_t

  implicit none

  type(c_ptr), intent(inout) :: sample_handle
  real(c_double), intent(inout) :: sqme

  type(simulation_api_t), pointer :: sample
  real(real64) :: s

  call c_f_pointer (sample_handle, sample)
  call sample%get_sqme (s)
  sqme = s

end subroutine whizard_sample_get_sqme


!!! ************************************************************
!!! Below: procedures used for intrinsic tests only
subroutine whizard_ut_setup (ut_name, u_log, results_handle) bind (C)
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_loc  !NODEP!
  use string_utils, only: string_c2f
  use io_units, only: free_unit
  use unit_tests, only: test_results_t

  implicit none

  character(c_char), dimension(*), intent(in) :: ut_name
  integer(c_int), intent(out) :: u_log
  type(c_ptr), intent(inout) :: results_handle

  type(test_results_t), pointer :: results

  allocate (results)
  results_handle = c_loc (results)

  u_log = free_unit ()
  open (unit = u_log, &
       file = "whizard_check." // string_c2f (ut_name) // ".log", &
       action = "write", status = "replace")

end subroutine whizard_ut_setup

subroutine whizard_ut_wrapup (u_log, results_handle) bind (C)
  use iso_fortran_env, only: output_unit  !NODEP!
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use unit_tests, only: test_results_t

  implicit none

  integer(c_int), value :: u_log
  type(c_ptr), intent(inout) :: results_handle

  type(test_results_t), pointer :: results

  call c_f_pointer (results_handle, results)
  call results%wrapup (output_unit)

  close (u_log)

end subroutine whizard_ut_wrapup

function whizard_ut_get_n_pass (results_handle) result (n_pass) bind (C)
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use unit_tests, only: test_results_t

  implicit none

  type(c_ptr), intent(inout) :: results_handle
  integer(c_int) :: n_pass

  type(test_results_t), pointer :: results

  call c_f_pointer (results_handle, results)
  n_pass = results%get_n_pass ()

end function whizard_ut_get_n_pass

function whizard_ut_get_n_fail (results_handle) result (n_fail) bind (C)
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use unit_tests, only: test_results_t

  implicit none

  type(c_ptr), intent(inout) :: results_handle
  integer(c_int) :: n_fail

  type(test_results_t), pointer :: results

  call c_f_pointer (results_handle, results)
  n_fail = results%get_n_fail ()

end function whizard_ut_get_n_fail

function whizard_ut_get_n_total (results_handle) result (n_total) bind (C)
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use unit_tests, only: test_results_t

  implicit none

  type(c_ptr), intent(inout) :: results_handle
  integer(c_int) :: n_total

  type(test_results_t), pointer :: results

  call c_f_pointer (results_handle, results)
  n_total = results%get_n_total ()

end function whizard_ut_get_n_total

subroutine whizard_ut_start (u_log, name) bind (C)
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use string_utils, only: string_c2f
  use unit_tests, only: start_test

  implicit none

  integer(c_int), value :: u_log
  character(c_char), dimension(*), intent(in) :: name

  call start_test (int (u_log), string_c2f (name))

end subroutine whizard_ut_start

subroutine whizard_ut_end (u_log, name, description, results_handle) bind (C)
  use iso_c_binding, only: c_int  !NODEP!
  use iso_c_binding, only: c_ptr  !NODEP!
  use iso_c_binding, only: c_char  !NODEP!
  use iso_c_binding, only: c_f_pointer  !NODEP!
  use string_utils, only: string_c2f
  use unit_tests, only: compare_test_results
  use unit_tests, only: test_results_t

  implicit none

  integer(c_int), value :: u_log
  character(c_char), dimension(*), intent(in) :: name
  character(c_char), dimension(*), intent(in) :: description
  type(c_ptr), intent(inout) :: results_handle

  type(test_results_t), pointer :: results

  integer :: u_test
  logical :: success

  open (newunit = u_test, file = string_c2f (name) // ".out", &
       action = "read", status = "old")
  call compare_test_results (u_test, int (u_log), string_c2f (name), success)
  call c_f_pointer (results_handle, results)
  call results%add (string_c2f (name), string_c2f (description), success)

end subroutine whizard_ut_end

