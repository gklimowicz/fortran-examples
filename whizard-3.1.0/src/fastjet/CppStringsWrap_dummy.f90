! Dummy implementations for C++ string library wrapper functions
subroutine cpp_str_delete (str) bind (C)
  use iso_c_binding
  type(c_ptr), value :: str
end subroutine cpp_str_delete

function cpp_str_length (str) bind (C) result (length)
  use iso_c_binding
  type(c_ptr), intent(in), value :: str
  integer(c_int) :: length
  length = 0
end function cpp_str_length

function cpp_str_get (str, i) bind (C) result (c)
  use iso_c_binding
  type(c_ptr), intent(in), value :: str
  integer(c_int), intent(in), value :: i
  character(c_char) :: c
  c = achar (0)
end function cpp_str_get
