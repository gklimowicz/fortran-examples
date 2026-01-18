program circe2_generate_program
  use kinds
  use circe2
  use tao_random_objects 
  implicit none
  type(circe2_state) :: c2s
  type(rng_tao), save :: rng
  character(len=1024) :: filename, design, buffer
  integer :: status, nevents, seed
  real(kind=default) :: roots
  real(kind=default), dimension(2) :: x
  integer :: i, ierror
  call get_command_argument (1, value = filename, status = status)
  if (status /= 0) filename = ""
  call get_command_argument (2, value = design, status = status)
  if (status /= 0) design = ""
  if (filename == "" .or. design == "") then
     print *, "usage: circe2_generate filename design [roots] [#events] [seed]"
     stop
  end if
  call get_command_argument (3, value = buffer, status = status)
  if (status == 0) then
     read (buffer, *, iostat = status) roots
     if (status /= 0) roots = 500
  else
     roots = 500
  end if
  call get_command_argument (4, value = buffer, status = status)
  if (status == 0) then
     read (buffer, *, iostat = status) nevents
     if (status /= 0) nevents = 1000
  else
     nevents = 1000
  end if
  call get_command_argument (5, value = buffer, status = status)
  if (status == 0) then
     read (buffer, *, iostat = status) seed
     if (status == 0) then
        call random2_seed (rng, seed)
     else
        call random2_seed (rng)
     end if
  else
     call random2_seed (rng)
  end if
  call circe2_load (c2s, trim(filename), trim(design), roots, ierror)
  if (ierror /= 0) then
     print *, "circe2_generate: failed to load design ", trim(design), &
          " for ", real (roots, kind=single), &
          " GeV from  ", trim(filename)
     stop
  end if
  do i = 1, nevents
     call circe2_generate (c2s, rng, x, [11, -11], [0, 0])
     write (*, '(F12.10,1X,F12.10)') x
  end do
  contains
     subroutine random2_seed (rng, seed)
       class(rng_tao), intent(inout) :: rng
       integer, intent(in), optional:: seed
       integer, dimension(8) :: date_time
       integer :: seed_value
       if (present (seed)) then
          seed_value = seed
       else
          call date_and_time (values = date_time)
          seed_value = product (date_time)
       endif
       call rng%init (seed_value)
     end subroutine random2_seed
end program circe2_generate_program
