! tao52_random_numbers.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module tao52_random_numbers
  use kinds
  implicit none
  integer, parameter, private :: tao_i32 = selected_int_kind (9)
  integer, parameter, private :: tao_r64 = selected_real_kind (15)
  private :: generate
  private :: seed_static, seed_state, seed_raw_state
  private :: seed_stateless
  private :: create_state_from_seed, create_raw_state_from_seed, &
             create_state_from_state, create_raw_state_from_state, &
             create_state_from_raw_state, create_raw_state_from_raw_st
  private :: destroy_state, destroy_raw_state
  public :: assignment(=)
  private :: copy_state, copy_raw_state, &
             copy_raw_state_to_state, copy_state_to_raw_state
  private :: write_state_unit, write_state_name
  private :: write_raw_state_unit, write_raw_state_name
  private :: read_state_unit, read_state_name
  private :: read_raw_state_unit, read_raw_state_name
  private :: find_free_unit
  public :: tao_random_marshal
  private :: marshal_state, marshal_raw_state
  public :: tao_random_marshal_size
  private :: marshal_state_size, marshal_raw_state_size
  public :: tao_random_unmarshal
  private :: unmarshal_state, unmarshal_raw_state
  public :: tao_random_number
  public :: tao_random_seed
  public :: tao_random_create
  public :: tao_random_destroy
  public :: tao_random_copy
  public :: tao_random_read
  public :: tao_random_write
  public :: tao_random_flush
  ! public :: tao_random_luxury
  public :: tao_random_test
  private :: write_state_array
  private :: read_state_array
  private :: real_stateless, real_array_stateless
  private :: real_static, real_state, real_array_static, real_array_state
  interface tao_random_seed
     module procedure seed_static, seed_state, seed_raw_state
  end interface
  interface tao_random_create
     module procedure create_state_from_seed, create_raw_state_from_seed, &
                      create_state_from_state, create_raw_state_from_state, &
                      create_state_from_raw_state, create_raw_state_from_raw_st
  end interface
  interface tao_random_destroy
     module procedure destroy_state, destroy_raw_state
  end interface
  interface tao_random_copy
     module procedure copy_state, copy_raw_state, &
                      copy_raw_state_to_state, copy_state_to_raw_state
  end interface
  interface assignment(=)
     module procedure copy_state, copy_raw_state, &
                      copy_raw_state_to_state, copy_state_to_raw_state
  end interface
  interface tao_random_write
     module procedure &
          write_state_unit, write_state_name, &
          write_raw_state_unit, write_raw_state_name
  end interface
  interface tao_random_read
     module procedure &
          read_state_unit, read_state_name, &
          read_raw_state_unit, read_raw_state_name
  end interface
  interface tao_random_marshal_size
     module procedure marshal_state_size, marshal_raw_state_size
  end interface
  interface tao_random_marshal
     module procedure marshal_state, marshal_raw_state
  end interface
  interface tao_random_unmarshal
     module procedure unmarshal_state, unmarshal_raw_state
  end interface
  interface tao_random_number
     module procedure real_static, real_state, real_array_static, real_array_state
  end interface
  integer, parameter, private :: K = 100, L = 37
  integer, parameter, private :: DEFAULT_BUFFER_SIZE = 1009
  integer, parameter, private :: MIN_UNIT = 11, MAX_UNIT = 99
  real(kind=tao_r64), parameter, private :: M = 1.0_tao_r64
  real(kind=tao_r64), dimension(K), save, private :: s_state
  logical, save, private :: s_virginal = .true.
  real(kind=tao_r64), dimension(DEFAULT_BUFFER_SIZE), save, private :: s_buffer
  integer, save, private :: s_buffer_end = size (s_buffer)
  integer, save, private :: s_last = size (s_buffer)
  type, public :: tao_random_raw_state
     private
     real(kind=tao_r64), dimension(K) :: x
  end type tao_random_raw_state
  type, public :: tao_random_state
     private
     type(tao_random_raw_state) :: state
     real(kind=tao_r64), dimension(:), pointer :: buffer => null ()
     integer :: buffer_end, last
  end type tao_random_state
contains
  subroutine seed_static (seed)
    integer, optional, intent(in) :: seed
    call seed_stateless (s_state, seed)
    s_virginal = .false.
    s_last = size (s_buffer)
  end subroutine seed_static
  elemental subroutine seed_raw_state (s, seed)
    type(tao_random_raw_state), intent(inout) :: s
    integer, optional, intent(in) :: seed
    call seed_stateless (s%x, seed)
  end subroutine seed_raw_state
  elemental subroutine seed_state (s, seed)
    type(tao_random_state), intent(inout) :: s
    integer, optional, intent(in) :: seed
    call seed_raw_state (s%state, seed)
    s%last = size (s%buffer)
  end subroutine seed_state
  elemental subroutine create_state_from_seed (s, seed, buffer_size)
    type(tao_random_state), intent(out) :: s
    integer, intent(in) :: seed
    integer, intent(in), optional :: buffer_size
    call create_raw_state_from_seed (s%state, seed)
    if (present (buffer_size)) then
       s%buffer_end = max (buffer_size, K)
    else
       s%buffer_end = DEFAULT_BUFFER_SIZE
    end if
    allocate (s%buffer(s%buffer_end))
    call tao_random_flush (s)
  end subroutine create_state_from_seed
  elemental subroutine create_state_from_state (s, state)
    type(tao_random_state), intent(out) :: s
    type(tao_random_state), intent(in) :: state
    call create_raw_state_from_raw_st (s%state, state%state)
    allocate (s%buffer(size(state%buffer)))
    call tao_random_copy (s, state)
  end subroutine create_state_from_state
  elemental subroutine create_state_from_raw_state &
       (s, raw_state, buffer_size)
    type(tao_random_state), intent(out) :: s
    type(tao_random_raw_state), intent(in) :: raw_state
    integer, intent(in), optional :: buffer_size
    call create_raw_state_from_raw_st (s%state, raw_state)
    if (present (buffer_size)) then
       s%buffer_end = max (buffer_size, K)
    else
       s%buffer_end = DEFAULT_BUFFER_SIZE
    end if
    allocate (s%buffer(s%buffer_end))
    call tao_random_flush (s)
  end subroutine create_state_from_raw_state
  elemental subroutine create_raw_state_from_seed (s, seed)
    type(tao_random_raw_state), intent(out) :: s
    integer, intent(in) :: seed
    call seed_raw_state (s, seed)
  end subroutine create_raw_state_from_seed
  elemental subroutine create_raw_state_from_state (s, state)
    type(tao_random_raw_state), intent(out) :: s
    type(tao_random_state), intent(in) :: state
    call copy_state_to_raw_state (s, state)
  end subroutine create_raw_state_from_state
  elemental subroutine create_raw_state_from_raw_st (s, raw_state)
    type(tao_random_raw_state), intent(out) :: s
    type(tao_random_raw_state), intent(in) :: raw_state
    call copy_raw_state (s, raw_state)
  end subroutine create_raw_state_from_raw_st
  elemental subroutine destroy_state (s)
    type(tao_random_state), intent(inout) :: s
    deallocate (s%buffer)
  end subroutine destroy_state
  elemental subroutine destroy_raw_state (s)
    type(tao_random_raw_state), intent(inout) :: s
  end subroutine destroy_raw_state
  elemental subroutine copy_state (lhs, rhs)
    type(tao_random_state), intent(inout) :: lhs
    type(tao_random_state), intent(in) :: rhs
    call copy_raw_state (lhs%state, rhs%state)
    if (size (lhs%buffer) /= size (rhs%buffer)) then
       deallocate (lhs%buffer)
       allocate (lhs%buffer(size(rhs%buffer)))
    end if
    lhs%buffer = rhs%buffer
    lhs%buffer_end = rhs%buffer_end
    lhs%last = rhs%last
  end subroutine copy_state
  elemental subroutine copy_raw_state (lhs, rhs)
    type(tao_random_raw_state), intent(out) :: lhs
    type(tao_random_raw_state), intent(in) :: rhs
    lhs%x = rhs%x
  end subroutine copy_raw_state
  elemental subroutine copy_raw_state_to_state (lhs, rhs)
    type(tao_random_state), intent(inout) :: lhs
    type(tao_random_raw_state), intent(in) :: rhs
    call copy_raw_state (lhs%state, rhs)
    call tao_random_flush (lhs)
  end subroutine copy_raw_state_to_state
  elemental subroutine copy_state_to_raw_state (lhs, rhs)
    type(tao_random_raw_state), intent(out) :: lhs
    type(tao_random_state), intent(in) :: rhs
    call copy_raw_state (lhs, rhs%state)
  end subroutine copy_state_to_raw_state
  elemental subroutine tao_random_flush (s)
    type(tao_random_state), intent(inout) :: s
    s%last = size (s%buffer)
  end subroutine tao_random_flush
  subroutine write_state_unit (s, unit)
    type(tao_random_state), intent(in) :: s
    integer, intent(in) :: unit
    write (unit = unit, fmt = *) "BEGIN TAO_RANDOM_STATE"
    call write_raw_state_unit (s%state, unit)
    write (unit = unit, fmt = "(2(1x,a16,1x,i10/),1x,a16,1x,i10)") &
         "BUFFER_SIZE", size (s%buffer), &
         "BUFFER_END", s%buffer_end, &
         "LAST", s%last
    write (unit = unit, fmt = *) "BEGIN BUFFER"
    call write_state_array (s%buffer, unit)
    write (unit = unit, fmt = *) "END BUFFER"
    write (unit = unit, fmt = *) "END TAO_RANDOM_STATE"
  end subroutine write_state_unit
  subroutine read_state_unit (s, unit)
    type(tao_random_state), intent(inout) :: s
    integer, intent(in) :: unit
    integer :: buffer_size
    read (unit = unit, fmt = *)
    call read_raw_state_unit (s%state, unit)
    read (unit = unit, fmt = "(2(1x,16x,1x,i10/),1x,16x,1x,i10)") &
         buffer_size, s%buffer_end, s%last
    read (unit = unit, fmt = *)
    if (buffer_size /= size (s%buffer)) then
       deallocate (s%buffer)
       allocate (s%buffer(buffer_size))
    end if
    call read_state_array (s%buffer, unit)
    read (unit = unit, fmt = *)
    read (unit = unit, fmt = *)
  end subroutine read_state_unit
  subroutine write_raw_state_unit (s, unit)
    type(tao_random_raw_state), intent(in) :: s
    integer, intent(in) :: unit
    write (unit = unit, fmt = *) "BEGIN TAO_RANDOM_RAW_STATE"
    call write_state_array (s%x, unit)
    write (unit = unit, fmt = *) "END TAO_RANDOM_RAW_STATE"
  end subroutine write_raw_state_unit
  subroutine read_raw_state_unit (s, unit)
    type(tao_random_raw_state), intent(inout) :: s
    integer, intent(in) :: unit
    read (unit = unit, fmt = *)
    call read_state_array (s%x, unit)
    read (unit = unit, fmt = *)
  end subroutine read_raw_state_unit
  subroutine find_free_unit (u, iostat)
    integer, intent(out) :: u
    integer, intent(out), optional :: iostat
    logical :: exists, is_open
    integer :: i, status
    do i = MIN_UNIT, MAX_UNIT
       inquire (unit = i, exist = exists, opened = is_open, &
                iostat = status)
       if (status == 0) then
          if (exists .and. .not. is_open) then
             u = i
             if (present (iostat)) then
                iostat = 0
             end if
             return
          end if
       end if
    end do
    if (present (iostat)) then
       iostat = -1
    end if
    u = -1
  end subroutine find_free_unit
  subroutine write_state_name (s, name)
    type(tao_random_state), intent(in) :: s
    character(len=*), intent(in) :: name
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", file = name)
    call write_state_unit (s, unit)
    close (unit = unit)
  end subroutine write_state_name
  subroutine write_raw_state_name (s, name)
    type(tao_random_raw_state), intent(in) :: s
    character(len=*), intent(in) :: name
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", file = name)
    call write_raw_state_unit (s, unit)
    close (unit = unit)
  end subroutine write_raw_state_name
  subroutine read_state_name (s, name)
    type(tao_random_state), intent(inout) :: s
    character(len=*), intent(in) :: name
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", file = name)
    call read_state_unit (s, unit)
    close (unit = unit)
  end subroutine read_state_name
  subroutine read_raw_state_name (s, name)
    type(tao_random_raw_state), intent(inout) :: s
    character(len=*), intent(in) :: name
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", file = name)
    call read_raw_state_unit (s, unit)
    close (unit = unit)
  end subroutine read_raw_state_name
  pure subroutine luxury_stateless &
       (buffer_size, buffer_end, last, consumption)
    integer, intent(in) :: buffer_size
    integer, intent(inout) :: buffer_end
    integer, intent(inout) :: last
    integer, intent(in) :: consumption
    if (consumption >= 1 .and. consumption <= buffer_size) then
       buffer_end = consumption
       last = min (last, buffer_end)
    else
       !!! print *, "tao_random_luxury: ", "invalid consumption ", &
       !!!      consumption, ", not in [ 1,", buffer_size, "]."
       buffer_end = buffer_size
    end if
  end subroutine luxury_stateless
  elemental subroutine luxury_state (s)
    type(tao_random_state), intent(inout) :: s
    call luxury_state_integer (s, size (s%buffer))
  end subroutine luxury_state
  elemental subroutine luxury_state_integer (s, consumption)
    type(tao_random_state), intent(inout) :: s
    integer, intent(in) :: consumption
    call luxury_stateless (size (s%buffer), s%buffer_end, s%last, consumption)
  end subroutine luxury_state_integer
  elemental subroutine luxury_state_real (s, consumption)
    type(tao_random_state), intent(inout) :: s
    real(kind=default), intent(in) :: consumption
    call luxury_state_integer (s, int (consumption * size (s%buffer)))
  end subroutine luxury_state_real
  subroutine luxury_static ()
    if (s_virginal) then
       call tao_random_seed ()
    end if
    call luxury_static_integer (size (s_buffer))
  end subroutine luxury_static
  subroutine luxury_static_integer (consumption)
    integer, intent(in) :: consumption
    if (s_virginal) then
       call tao_random_seed ()
    end if
    call luxury_stateless (size (s_buffer), s_buffer_end, s_last, consumption)
  end subroutine luxury_static_integer
  subroutine luxury_static_real (consumption)
    real(kind=default), intent(in) :: consumption
    if (s_virginal) then
       call tao_random_seed ()
    end if
    call luxury_static_integer (int (consumption * size (s_buffer)))
  end subroutine luxury_static_real
  pure subroutine generate (a, state)
    real(kind=tao_r64), dimension(:), intent(inout) :: a
    real(kind=tao_r64), dimension(:), intent(inout) :: state
    integer :: j, n
    n = size (a)
    a(1:K) = state(1:K)
    do j = K+1, n
       a(j) = modulo (a(j-K) + a(j-L), M)
    end do
    state(1:L) = modulo (a(n+1-K:n+L-K) + a(n+1-L:n), M)
    do j = L+1, K
       state(j) = modulo (a(n+j-K) + state(j-L), M)
    end do
  end subroutine generate
  pure subroutine seed_stateless (state, seed)
    real(kind=tao_r64), dimension(:), intent(out) :: state
    integer, optional, intent(in) :: seed
    integer, parameter :: DEFAULT_SEED = 0
    integer, parameter :: MAX_SEED = 2**30 - 3
    integer, parameter :: TT = 70
    real(kind=tao_r64), parameter :: ULP = 2.0_tao_r64**(-52)
    real(kind=tao_r64), dimension(2*K-1) :: x
    real(kind=tao_r64) :: ss
    integer :: seed_value, t, s, j
    if (present (seed)) then
       seed_value = modulo (seed, MAX_SEED + 1)
    else
       seed_value = DEFAULT_SEED
    end if
    ss = 2*ULP * (seed_value + 2)
    do j = 1, K
       x(j) = ss
       ss = 2*ss
       if (ss >= 1) then
          ss = ss - 1 + 2*ULP
       end if
    end do
    x(K+1:2*K-1) = 0.0
    x(2) = x(2) + ULP
    s = seed_value
    t = TT - 1
    do
       x(3:2*K-1:2) = x(2:K)
       x(2:2*K-2:2) = 0
       do j = 2*K-1, K+1, -1
          x(j-(K-L)) = modulo (x(j-(K-L)) + x(j), M)
          x(j-K) = modulo (x(j-K) + x(j), M)
       end do
       if (modulo (s, 2) == 1) THEN
          x(2:K+1) = x(1:K)
          x(1) = x(K+1)
          x(L+1) = modulo (x(L+1) + x(K+1), M)
       end if
       if (s /= 0) then
          s = s / 2
       else
          t = t - 1
       end if
       if (t <= 0) then
          exit
       end if
    end do
    state(1:K-L) = x(L+1:K)
    state(K-L+1:K) = x(1:L)
    do j = 1, 10
       call generate (x, state)
    end do
  end subroutine seed_stateless
  subroutine write_state_array (a, unit)
    real(kind=tao_r64), dimension(:), intent(in) :: a
    integer, intent(in) :: unit
    integer :: i
    do i = 1, size (a)
       write (unit = unit, fmt = "(1x,i10,1x,f30.0)") i, 2.0_tao_r64**52 * a(i)
    end do
  end subroutine write_state_array
  subroutine read_state_array (a, unit)
    real(kind=tao_r64), dimension(:), intent(inout) :: a
    integer, intent(in) :: unit
    real(kind=tao_r64) :: x
    integer :: i, idum
    do i = 1, size (a)
       read (unit = unit, fmt = *) idum, x
       a(i) = 2.0_tao_r64**(-52) * x
    end do
  end subroutine read_state_array
  pure subroutine marshal_state (s, ibuf, dbuf)
    type(tao_random_state), intent(in) :: s
    integer, dimension(:), intent(inout) :: ibuf
    real(kind=tao_r64), dimension(:), intent(inout) :: dbuf
    integer :: buf_size
    buf_size = size (s%buffer)
    ibuf(1) = s%buffer_end
    ibuf(2) = s%last
    ibuf(3) = buf_size
    dbuf(1:buf_size) = s%buffer
    call marshal_raw_state (s%state, ibuf(4:), dbuf(buf_size+1:))
  end subroutine marshal_state
  pure subroutine marshal_state_size (s, iwords, dwords)
    type(tao_random_state), intent(in) :: s
    integer, intent(out) :: iwords, dwords
    call marshal_raw_state_size (s%state, iwords, dwords)
    iwords = iwords + 3
    dwords = dwords + size(s%buffer)
  end subroutine marshal_state_size
  pure subroutine unmarshal_state (s, ibuf, dbuf)
    type(tao_random_state), intent(inout) :: s
    integer, dimension(:), intent(in) :: ibuf
    real(kind=tao_r64), dimension(:), intent(in) :: dbuf
    integer :: buf_size
    s%buffer_end = ibuf(1)
    s%last = ibuf(2)
    buf_size = ibuf(3)
    s%buffer = dbuf(1:buf_size)
    call unmarshal_raw_state (s%state, ibuf(4:), dbuf(buf_size+1:))
  end subroutine unmarshal_state
  pure subroutine marshal_raw_state (s, ibuf, dbuf)
    type(tao_random_raw_state), intent(in) :: s
    integer, dimension(:), intent(inout) :: ibuf
    real(kind=tao_r64), dimension(:), intent(inout) :: dbuf
    ibuf(1) = size (s%x)
    dbuf(1:size(s%x)) = s%x
  end subroutine marshal_raw_state
  pure subroutine marshal_raw_state_size (s, iwords, dwords)
    type(tao_random_raw_state), intent(in) :: s
    integer, intent(out) :: iwords, dwords
    iwords = 1
    dwords = size (s%x)
  end subroutine marshal_raw_state_size
  pure subroutine unmarshal_raw_state (s, ibuf, dbuf)
    type(tao_random_raw_state), intent(inout) :: s
    integer, dimension(:), intent(in) :: ibuf
    real(kind=tao_r64), dimension(:), intent(in) :: dbuf
    integer :: buf_size
    buf_size = ibuf(1)
    s%x = dbuf(1:buf_size)
  end subroutine unmarshal_raw_state
  pure subroutine real_stateless (state, buffer, buffer_end, last, r)
    real(kind=tao_r64), dimension(:), intent(inout) :: state, buffer
    integer, intent(in) :: buffer_end
    integer, intent(inout) :: last
    real(kind=default), intent(out) :: r
    integer, parameter :: NORM = 1
    last = last + 1
    if (last > buffer_end) then
       call generate (buffer, state)
       last = 1
    end if
    r = NORM * buffer(last) 
  end subroutine real_stateless
  pure subroutine real_array_stateless &
       (state, buffer, buffer_end, last, v, num)
    real(kind=tao_r64), dimension(:), intent(inout) :: state, buffer
    integer, intent(in) :: buffer_end
    integer, intent(inout) :: last
    real(kind=default), dimension(:), intent(out) :: v
    integer, optional, intent(in) :: num
    integer, parameter :: NORM = 1
    integer :: nu, done, todo, chunk
    if (present (num)) then
       nu = num
    else
       nu = size (v)
    end if
    if (last >= buffer_end) then
       call generate (buffer, state)
       last = 0
    end if
    done = 0
    todo = nu
    chunk = min (todo, buffer_end - last)
    v(1:chunk) = NORM * buffer(last+1:last+chunk)
    do
       last = last + chunk
       done = done + chunk
       todo = todo - chunk
       chunk = min (todo, buffer_end)
       if (chunk <= 0) then
          exit
       end if
       call generate (buffer, state)
       last = 0
       v(done+1:done+chunk) = NORM * buffer(1:chunk)
    end do
  end subroutine real_array_stateless
  elemental subroutine real_state (s, r)
    type(tao_random_state), intent(inout) :: s
    real(kind=default), intent(out) :: r
    call real_stateless (s%state%x, s%buffer, s%buffer_end, s%last, r)
  end subroutine real_state
  pure subroutine real_array_state (s, v, num)
    type(tao_random_state), intent(inout) :: s
    real(kind=default), dimension(:), intent(out) :: v
    integer, optional, intent(in) :: num
    call real_array_stateless &
         (s%state%x, s%buffer, s%buffer_end, s%last, v, num)
  end subroutine real_array_state
  subroutine real_static (r)
    real(kind=default), intent(out) :: r
    if (s_virginal) then
       call tao_random_seed ()
    end if
    call real_stateless (s_state, s_buffer, s_buffer_end, s_last, r)
  end subroutine real_static
  subroutine real_array_static (v, num)
    real(kind=default), dimension(:), intent(out) :: v
    integer, optional, intent(in) :: num
    if (s_virginal) then
       call tao_random_seed ()
    end if
    call real_array_stateless &
         (s_state, s_buffer, s_buffer_end, s_last, v, num)
  end subroutine real_array_static
  subroutine tao_random_test (name)
    character(len=*), optional, intent(in) :: name
    character(len=*), parameter :: &
         OK = "(1x,f22.20,' is ok.')", &
         NOT_OK = "(1x,f22.20,' is not ok, (A_2027082 ',f22.20,')!')"
    integer, parameter :: &
         SEED = 310952, &
         N = 2009, M = 1009, &
         N_SHORT = 1984
    real(kind=default), parameter :: &
         A_2027082 = 0.36410514377569680455_tao_r64
    real(kind=default), dimension(N) :: a
    type(tao_random_state) :: s, t
    integer, dimension(:), allocatable :: ibuf
    real(kind=tao_r64), dimension(:), allocatable :: dbuf
    integer :: i, ibuf_size, dbuf_size
    print *, "testing the 52-bit tao_random_numbers ..."
    ! call tao_random_luxury ()
    call tao_random_seed (SEED)
    do i = 1, N+1
       call tao_random_number (a, M)
    end do
    if (a(1) == A_2027082) then
       print OK, a(1)
    else
       print NOT_OK, a(1), A_2027082
       stop 1
    end if
    call tao_random_seed (SEED)
    do i = 1, M+1
       call tao_random_number (a)
    end do
    if (a(1) == A_2027082) then
       print OK, a(1)
    else
       print NOT_OK, a(1), A_2027082
       stop 1
    end if
    print *, "testing the stateless stuff ..."
    call tao_random_create (s, SEED)
    do i = 1, N_SHORT
       call tao_random_number (s, a, M)
    end do
    call tao_random_create (t, s)
    do i = 1, N+1 - N_SHORT
       call tao_random_number (s, a, M)
    end do
    if (a(1) == A_2027082) then
       print OK, a(1)
    else
       print NOT_OK, a(1), A_2027082
       stop 1
    end if
    do i = 1, N+1 - N_SHORT
       call tao_random_number (t, a, M)
    end do
    if (a(1) == A_2027082) then
       print OK, a(1)
    else
       print NOT_OK, a(1), A_2027082
       stop 1
    end if
    if (present (name)) then
       print *, "testing I/O ..."
       call tao_random_seed (s, SEED)
       do i = 1, N_SHORT
          call tao_random_number (s, a, M)
       end do
       call tao_random_write (s, name)
       do i = 1, N+1 - N_SHORT
          call tao_random_number (s, a, M)
       end do
       if (a(1) == A_2027082) then
          print OK, a(1)
       else
          print NOT_OK, a(1), A_2027082
          stop 1
       end if
       call tao_random_read (s, name)
       do i = 1, N+1 - N_SHORT
          call tao_random_number (s, a, M)
       end do
       if (a(1) == A_2027082) then
          print OK, a(1)
       else
          print NOT_OK, a(1), A_2027082
          stop 1
       end if
    end if
    print *, "testing marshaling/unmarshaling ..."
    call tao_random_seed (s, SEED)
    do i = 1, N_SHORT
       call tao_random_number (s, a, M)
    end do
    call tao_random_marshal_size (s, ibuf_size, dbuf_size)
    allocate (ibuf(ibuf_size), dbuf(dbuf_size))
    call tao_random_marshal (s, ibuf, dbuf)
    do i = 1, N+1 - N_SHORT
       call tao_random_number (s, a, M)
    end do
    if (a(1) == A_2027082) then
       print OK, a(1)
    else
       print NOT_OK, a(1), A_2027082
       stop 1
    end if
    call tao_random_unmarshal (s, ibuf, dbuf)
    do i = 1, N+1 - N_SHORT
       call tao_random_number (s, a, M)
    end do
    if (a(1) == A_2027082) then
       print OK, a(1)
    else
       print NOT_OK, a(1), A_2027082
       stop 1
    end if
  end subroutine tao_random_test
end module tao52_random_numbers
