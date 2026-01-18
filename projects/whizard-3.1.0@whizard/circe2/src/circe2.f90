! circe2.f90 -- correlated beam spectra for linear colliders
! Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!
! Circe2 is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! Circe2 is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!-----------------------------------------------------------------------
module circe2
  use kinds
  implicit none
  private
  integer, parameter :: POLAVG = 1, POLHEL = 2, POLGEN = 3
  type circe2_division
    real(kind=default), dimension(:), allocatable :: x
    integer, dimension(:), allocatable :: map
    real(kind=default), dimension(:), allocatable :: y
    real(kind=default), dimension(:), allocatable :: alpha, xi, eta, a, b
  end type circe2_division
  type circe2_channel
    real(kind=default), dimension(:), allocatable :: wgt
    type(circe2_division), dimension(2) :: d
    real(kind=default), dimension(:,:), allocatable :: val
    logical :: triang
    integer, dimension(2) :: pid, pol
    real(kind=default) :: lumi
  end type circe2_channel
  type circe2_state
    type(circe2_channel), dimension(:), allocatable :: ch
    real(kind=default), dimension(:), allocatable :: cwgt
    integer :: polspt
  end type circe2_state
  public :: circe2_state
  public :: rng_type
  type, abstract :: rng_type
    contains
      procedure(rng_generate), deferred :: generate
  end type rng_type
  abstract interface
     subroutine rng_generate (rng_obj, u)
       import :: rng_type, default
       class(rng_type), intent(inout) :: rng_obj
       real(kind=default), intent(out) :: u
     end subroutine rng_generate
  end interface
  public :: circe2_generate
  interface circe2_generate
     module procedure circe2_generate_ph
  end interface circe2_generate
  interface circe2_generate
     module procedure circe2_generate_channel
  end interface circe2_generate
  public :: circe2_choose_channel
  interface circe2_choose_channel
     module procedure circe2_choose_channel
  end interface circe2_choose_channel
  public :: circe2_generate_polarized
  interface circe2_generate_polarized
     module procedure circe2_generate_polarized
  end interface circe2_generate_polarized
  public :: circe2_luminosity
  public :: circe2_distribution
  interface circe2_distribution
     module procedure circe2_distribution_ph
  end interface circe2_distribution
  interface circe2_distribution
     module procedure circe2_distribution_channel
  end interface circe2_distribution
  public :: circe2_density_matrix
  public :: circe2_load
  integer, parameter, public :: &
       EOK = 0, EFILE = -1, EMATCH = -2, EFORMT = -3, ESIZE = -4
contains
  subroutine circe2_generate_ph (c2s, rng, y, p, h)
    type(circe2_state), intent(in) :: c2s
    class(rng_type), intent(inout) :: rng
    real(kind=default), dimension(:), intent(out) :: y
    integer, dimension(:), intent(in) :: p
    integer, dimension(:), intent(in) :: h
    integer :: i, ic
    ic = 0
    if ((c2s%polspt == POLAVG .or. c2s%polspt == POLGEN) .and. any (h /= 0)) then
       write (*, '(2A)') 'circe2: current beam description ', &
            'supports only polarization averages'
    else if (c2s%polspt == POLHEL .and. any (h == 0)) then
       write (*, '(2A)') 'circe2: polarization averages ', &
            'not supported by current beam description'
    else
       do i = 1, size (c2s%ch)
          if (all (p == c2s%ch(i)%pid .and. h == c2s%ch(i)%pol)) then
             ic = i
          end if
       end do
    end if
    if (ic <= 0) then
       write (*, '(A,2I4,A,2I3)') &
            'circe2: no channel for particles', p, &
            ' and polarizations', h
       y = - huge (y)
       return
    end if
    call circe2_generate_channel (c2s%ch(ic), rng, y)
  end subroutine circe2_generate_ph
  !-----------------------------------------------------------------------
  subroutine circe2_generate_channel (ch, rng, y)
    type(circe2_channel), intent(in) :: ch
    class(rng_type), intent(inout) :: rng
    real(kind=default), dimension(:), intent(out) :: y
    integer :: i, d, ibot, itop
    integer, dimension(2) :: ib
    real(kind=default), dimension(2) :: x, v
    real(kind=default) :: u, tmp
    call rng%generate (u)
    ibot = 0
    itop = ubound (ch%wgt, dim=1)
    do
       if (itop <= ibot + 1) then
          i = ibot + 1
          exit
       else
          i = (ibot + itop) / 2
          if (u < ch%wgt(i)) then
             itop = i
          else
             ibot = i
          end if
       end if
    end do
    ib(2) = 1 + (i - 1) / ubound (ch%d(1)%x, dim=1)
    ib(1) = i - (ib(2) - 1) * ubound (ch%d(1)%x, dim=1)
    call rng%generate (v(1))
    call rng%generate (v(2))
    do d = 1, 2
      x(d) = ch%d(d)%x(ib(d))*v(d) + ch%d(d)%x(ib(d)-1)*(1-v(d))
    end do
    y = circe2_map (ch%d, x, ib)
    if (ch%triang) then
       y(2) = y(1) * y(2)
       call rng%generate (u)
       if (2*u >= 1) then
          tmp = y(1)
          y(1) = y(2)
          y(2) = tmp
       end if
    end if
  end subroutine circe2_generate_channel
  !-----------------------------------------------------------------------
  elemental function circe2_map (d, x, b) result (y)
     type(circe2_division), intent(in) :: d
     real(kind=default), intent(in) :: x
     integer, intent(in) :: b
     real(kind=default) :: y
     real(kind=default) :: z
     select case (d%map(b))
     case (0)
        y = x
     case (1)
        z = d%a(b) * (x - d%xi(b))
        if (abs (z) <= tiny (z)) then
           z = abs (z)
        end if
        y = z**d%alpha(b) / d%b(b) + d%eta(b)
     case (2)
        y = d%a(b) * tan (d%a(b)*(x-d%xi(b)) / d%b(b)**2) + d%eta(b)
     case default
        y = - huge (y)
     end select
  end function circe2_map
  elemental function circe2_jacobian (d, y, b) result (j)
     type(circe2_division), intent(in) :: d
     real(kind=default), intent(in) :: y
     integer, intent(in) :: b
     real(kind=default) :: j
     select case (d%map(b))
     case (0)
        j = 1
     case (1)
        j = d%b(b) / (d%a(b)*d%alpha(b)) &
          * (d%b(b)*(y-d%eta(b)))**(1/d%alpha(b)-1)
     case (2)
        j = d%b(b)**2 / ((y-d%eta(b))**2 + d%a(b)**2)
     case default
        j = - huge (j)
     end select
  end function circe2_jacobian
  subroutine circe2_choose_channel (c2s, rng, p, h)
    type(circe2_state), intent(in) :: c2s
    class(rng_type), intent(inout) :: rng
    integer, dimension(:), intent(out) :: p, h
    integer :: ic, ibot, itop
    real(kind=default) :: u
    call rng%generate (u)
    ibot = 0
    itop = size (c2s%ch)
    do
       if (itop <= ibot + 1) then
          ic = ibot + 1
          p = c2s%ch(ic)%pid
          h = c2s%ch(ic)%pol
          return
       else
          ic = (ibot + itop) / 2
          if (u < c2s%cwgt(ic)) then
             itop = ic
          else
             ibot = ic
          end if
       end if
    end do
    write (*, '(A)') 'circe2: internal error'
    stop
  end subroutine circe2_choose_channel
  subroutine circe2_generate_polarized (c2s, rng, p, pol, x)
    type(circe2_state), intent(in) :: c2s
    class(rng_type), intent(inout) :: rng
    integer, dimension(:), intent(out) :: p
    real(kind=default), intent(out) :: pol(0:3,0:3)
    real(kind=default), dimension(:), intent(out) :: x
    integer, dimension(2) :: h
    integer :: i1, i2
    real(kind=default) :: pol00
    call circe2_choose_channel (c2s, rng, p, h)
    call circe2_generate (c2s, rng, x, p, h)
    call circe2_density_matrix (c2s, pol, p, x)
    pol00 = pol(0,0)
    do i1 = 0, 3
       do i2 = 0, 3
          pol(i1,i2) = pol(i1,i2) / pol00
       end do
    end do
  end subroutine circe2_generate_polarized
  function circe2_luminosity (c2s, p, h)
    type(circe2_state), intent(in) :: c2s
    integer, dimension(:), intent(in) :: p
    integer, dimension(:), intent(in) :: h
    real(kind=default) :: circe2_luminosity
    integer :: ic
    circe2_luminosity = 0
    do ic = 1, size (c2s%ch)
       if (       all (p == c2s%ch(ic)%pid .or. p == 0) &
            .and. all (h == c2s%ch(ic)%pol .or. h == 0)) then
          circe2_luminosity = circe2_luminosity + c2s%ch(ic)%lumi
       end if
    end do
  end function circe2_luminosity
  !-----------------------------------------------------------------------
  function circe2_distribution_ph (c2s, p, h, yy)
    type(circe2_state), intent(in) :: c2s
    integer, dimension(:), intent(in) :: p
    real(kind=default), dimension(:), intent(in)  :: yy
    integer, dimension(:), intent(in) :: h
    real(kind=default) :: circe2_distribution_ph
    integer :: i, ic
    ic = 0
    if ((c2s%polspt == POLAVG .or. c2s%polspt == POLGEN) .and. any (h /= 0)) then
       write (*, '(2A)') 'circe2: current beam description ', &
            'supports only polarization averages'
    else if (c2s%polspt == POLHEL .and. any (h == 0)) then
       write (*, '(2A)') 'circe2: polarization averages ', &
            'not supported by current beam description'
    else
       do i = 1, size (c2s%ch)
          if (all (p == c2s%ch(i)%pid .and. h == c2s%ch(i)%pol)) then
             ic = i
          end if
       end do
    end if
    if (ic <= 0) then
       circe2_distribution_ph = 0
    else
       circe2_distribution_ph = circe2_distribution_channel (c2s%ch(ic), yy)
    end if
  end function circe2_distribution_ph
  !-----------------------------------------------------------------------
  function circe2_distribution_channel (ch, yy)
    type(circe2_channel), intent(in) :: ch
    real(kind=default), dimension(:), intent(in)  :: yy
    real(kind=default) :: circe2_distribution_channel
    real(kind=default), dimension(2) :: y
    integer :: d, ibot, itop
    integer, dimension(2) :: ib
    if (ch%triang) then
       y(1) = maxval (yy)
       y(2) = minval (yy) / y(1)
    else
       y = yy
    end if
    if (      y(1) < ch%d(1)%y(0) &
         .or. y(1) > ch%d(1)%y(ubound (ch%d(1)%y, dim=1)) &
         .or. y(2) < ch%d(2)%y(0) &
         .or. y(2) > ch%d(2)%y(ubound (ch%d(2)%y, dim=1))) then
       circe2_distribution_channel = 0
       return
    end if
    do d = 1, 2
       ibot = 0
       itop = ubound (ch%d(d)%x, dim=1)
       search: do
          if (itop <= ibot + 1) then
             ib(d) = ibot + 1
             exit search
          else
             ib(d) = (ibot + itop) / 2
             if (y(d) < ch%d(d)%y(ib(d))) then
                itop = ib(d)
             else
                ibot = ib(d)
             end if
          end if
       end do search
    end do
    circe2_distribution_channel = &
        ch%val(ib(1),ib(2)) * product (circe2_jacobian (ch%d, y, ib))
    if (ch%triang) then
       circe2_distribution_channel = circe2_distribution_channel / y(1)
    end if
  end function circe2_distribution_channel
  !-----------------------------------------------------------------------
  subroutine circe2_density_matrix (c2s, pol, p, x)
    type(circe2_state), intent(in) :: c2s
    real(kind=default), dimension(0:,0:), intent(out) :: pol
    integer, dimension(:), intent(in) :: p
    real(kind=default), dimension(:), intent(in) :: x
    if (c2s%polspt /= POLGEN) then
       write (*, '(2A)') 'circe2: current beam ', &
            'description supports no density matrices'
       return
    end if
    print *, 'circe2: circe2_density_matrix not implemented yet!'
    if (p(1) < p(2) .and. x(1) < x(2)) then
       ! nonsense test to suppress warning
    end if
    pol = 0
  end subroutine circe2_density_matrix
  !-----------------------------------------------------------------------
  subroutine circe2_load (c2s, file, design, roots, ierror)
    type(circe2_state), intent(out) :: c2s
    character(len=*), intent(in) :: file, design
    real(kind=default), intent(in) :: roots
    integer, intent(out) :: ierror
    character(len=72) :: buffer, fdesgn, fpolsp
    real(kind=default) :: froots
    integer :: lun, loaded, prefix
    logical match
    integer :: ic, nc
    integer :: status
    logical exists, isopen
    scan: do lun = 10, 99
       inquire (unit = lun, exist = exists, opened = isopen, iostat = status)
       if (status == 0 .and. exists .and. .not.isopen) exit scan
    end do scan
    if (lun > 99) lun = -1
    if (lun < 0) then
       write (*, '(A)') 'circe2_load: no free unit'
       ierror = ESIZE
       return
    end if
    loaded = 0
    open (unit = lun, file = file, status = 'old', iostat = status)
    if (status /= 0) then
       write (*, '(2A)') 'circe2_load: can''t open ', file
       ierror = EFILE
       return
    end if
    if (ierror .gt. 0) then
       write (*, '(2A)') 'circe2_load: ', 'Version 3.1.0'                         
    end if
    prefix = index (design, '*') - 1
    do
       find_circe2: do
          skip_comments: do
             read (lun, '(A)', iostat = status) buffer
             if (status /= 0) then
                close (unit = lun)
                if (loaded > 0) then           
                   ierror = EOK
                else
                   ierror = EMATCH
                end if
                return
             else
                if (buffer(1:6) == 'CIRCE2') then
                   exit find_circe2
                else if (buffer(1:1) == '!') then
                   if (ierror > 0) then
                      write (*, '(A)') buffer
                   end if
                else
                   exit skip_comments
                end if
              end if
           end do skip_comments
           write (*, '(A)') 'circe2_load: invalid file'
           ierror = EFORMT
           return
        end do find_circe2
       if (buffer(8:15) == 'FORMAT#1') then
          read (lun, *)
          read (lun, *) fdesgn, froots
          match = .false.
          if (fdesgn == design) then
             match = .true.
          else if (prefix == 0) then
             match = .true.
          else if (prefix .gt. 0) then
             if (fdesgn(1:min(prefix,len(fdesgn))) &
                  == design(1:min(prefix,len(design)))) then
                match = .true.
             end if
          end if
          if (match .and. abs (froots - roots) <= 1) then
             read (lun, *) 
             read (lun, *) nc, fpolsp
             allocate (c2s%ch(nc), c2s%cwgt(0:nc))
             if (fpolsp(1:1)=='a' .or. fpolsp(1:1)=='A') then
                c2s%polspt = POLAVG
             else if (fpolsp(1:1)=='h' .or. fpolsp(1:1)=='H') then
                c2s%polspt = POLHEL
             else if (fpolsp(1:1)=='d' .or. fpolsp(1:1)=='D') then
                c2s%polspt = POLGEN
             else
                write (*, '(A,I5)') 'circe2_load: invalid polarization support: ', fpolsp
                ierror = EFORMT
                return
             end if
             c2s%cwgt(0) = 0
             do ic = 1, nc
                call circe2_load_channel (c2s%ch(ic), c2s%polspt, lun, ierror)
                c2s%cwgt(ic) = c2s%cwgt(ic-1) + c2s%ch(ic)%lumi
             end do
             c2s%cwgt = c2s%cwgt / c2s%cwgt(nc)
             loaded = loaded + 1
          else
             skip_data: do
                read (lun, *) buffer
                if (buffer(1:6) == 'ECRIC2') then
                   exit skip_data
                end if
             end do skip_data
             cycle
          end if
       else
          write (*, '(2A)') 'circe2_load: invalid format: ', buffer(8:72)
          ierror = EFORMT
          return
       end if
       read (lun, '(A)') buffer
       if (buffer(1:6) /= 'ECRIC2') then
          write (*, '(A)') 'circe2_load: invalid file'
          ierror = EFORMT
          return
       end if
    end do
  end subroutine circe2_load
  !-----------------------------------------------------------------------
  subroutine circe2_load_channel (ch, polspt, lun, ierror)
    type(circe2_channel), intent(out) :: ch
    integer, intent(in) :: polspt, lun
    integer, intent(out) :: ierror
    integer :: d, i, ib
    integer :: i1, i2
    integer, dimension(2) :: nb
    real(kind=default) :: w
    read (lun, *)
    read (lun, *) ch%pid(1), ch%pol(1), ch%pid(2), ch%pol(2), ch%lumi
    if (polspt == POLAVG .and. any (ch%pol /= 0)) then
       write (*, '(A)') 'circe2_load: expecting averaged polarization'
       ierror = EFORMT
       return
    else if (polspt == POLHEL .and. any (ch%pol == 0)) then
       write (*, '(A)') 'circe2_load: expecting helicities'
       ierror = EFORMT
       return
    else if (polspt == POLGEN) then
       write (*, '(A)') 'circe2_load: general polarizations not supported yet'
       ierror = EFORMT
       return
    else if (polspt == POLGEN .and. any (ch%pol /= 0)) then
       write (*, '(A)') 'circe2_load: expecting pol = 0'
       ierror = EFORMT
       return
    end if
    read (lun, *)
    read (lun, *) nb, ch%triang
    do d = 1, 2
       read (lun, *)
       allocate (ch%d(d)%x(0:nb(d)), ch%d(d)%y(0:nb(d)))
       allocate (ch%d(d)%map(nb(d)), ch%d(d)%alpha(nb(d)))
       allocate (ch%d(d)%xi(nb(d)), ch%d(d)%eta(nb(d)))
       allocate (ch%d(d)%a(nb(d)), ch%d(d)%b(nb(d)))
       read (lun, *) ch%d(d)%x(0)
       do ib = 1, nb(d)
          read (lun, *) ch%d(d)%x(ib), ch%d(d)%map(ib), &
               ch%d(d)%alpha(ib), ch%d(d)%xi(ib), ch%d(d)%eta(ib), &
               ch%d(d)%a(ib), ch%d(d)%b(ib)
          if (ch%d(d)%map(ib) < 0 .or. ch%d(d)%map(ib) > 2) then
             write (*, '(A,I3)') 'circe2_load: invalid map: ', ch%d(d)%map(ib)
             ierror = EFORMT
             return
          end if
       end do
    end do
    do d = 1, 2
       do i = 0, ubound (ch%d(d)%x, dim=1)
          ch%d(d)%y(i) = circe2_map (ch%d(d), ch%d(d)%x(i), max (i, 1))
       end do
    end do
    read (lun, *)
    allocate (ch%wgt(0:product(nb)), ch%val(nb(1),nb(2)))
    ch%wgt(0) = 0
    do i = 1, ubound (ch%wgt, dim=1)
       read (lun, *) w
       ch%wgt(i) = ch%wgt(i-1) + w
       i2 = 1 + (i - 1) / ubound (ch%d(1)%x, dim=1)
       i1 = i - (i2 - 1) * ubound (ch%d(1)%x, dim=1)
       ch%val(i1,i2) = w &
            / (  (ch%d(1)%x(i1) - ch%d(1)%x(i1-1)) &
               * (ch%d(2)%x(i2) - ch%d(2)%x(i2-1)))
    end do
    ch%wgt(ubound (ch%wgt, dim=1)) = 1
  end subroutine circe2_load_channel
end module circe2
