  module int_routines
    use kinds

    implicit none

  public :: gethst
  public :: puthst
  public :: putcol   

contains

  subroutine gethst (tag, ndata, x, f, x1, x2, pwr, scale)
    character(len=2) :: tag
    integer :: ndata
    real(kind=double), dimension(2,NDATA,NDATA) :: x
    real(kind=double), dimension(0:NDATA+1,0:NDATA+1) :: f
    real(kind=double), dimension(NDATA) :: x1, x2
    real(kind=double) :: pwr, scale, s, t, tmp
    real(kind=double), dimension(3) :: dum
    integer :: i, j
        open (10, file = 'lumidiff-'//tag//'.dat')
        read (10, *) pwr
        s = 0d0
        do i = 1, ndata
           do j = 1, ndata
              read (10, *) x(1,i,j), x(2,i,j), f(i,j)
              s = s + f(i,j)
           end do   
        end do   
        t = s
        do i = 1, NDATA
           read (10, *) dum(1), f(i,0), dum(2), f(i,NDATA+1), dum(3)
           t = t + f(i,0) + f(i,NDATA+1)
        end do
        do i = 1, NDATA
           read (10, *) dum(1), f(0,i), dum(2), f(NDATA+1,i), dum(3)
           t = t + f(0,i) + f(NDATA+1,i)
        end do
        read (10, *) f(0,0), dum(1), f(0,NDATA+1), dum(2)
        t = t + f(0,0) + f(0,NDATA+1)
        read (10, *) f(NDATA+1,0), dum(1), f(NDATA+1,NDATA+1), dum(2)
        t = t + f(NDATA+1,0) + f(NDATA+1,NDATA+1)
        close (10)
        do i = 1, ndata
           x1(i) = x(1,i,1)
           x2(i) = x(2,1,i)
           do j = 2, ndata
              if (x(1,i,j) .ne. x1(i)) then
                 print *, 'error: grid not rectangular: ', &
                         x(1,i,j), ' != ', x1(i)
                 stop
              endif
              if (x(2,j,i) .ne. x2(i)) then
                 print *, 'error: grid not rectangular: ', &
                         x(2,j,i), ' != ', x2(i)
                 stop
              endif
           end do   
        end do   
          if (scale .eq. 0d0) scale = t
        do i = 1, NDATA
           do j = 1, NDATA
              f(i,j) = dble(NDATA*NDATA) * f(i,j) / scale
           end do   
           f(0,i) = dble(NDATA) * f(0,i) / scale
           f(NDATA+1,i) = dble(NDATA) * f(NDATA+1,i) / scale
           f(i,0) = dble(NDATA) * f(i,0) / scale
           f(i,NDATA+1) = dble(NDATA) * f(i,NDATA+1) / scale
        end do   
        f(0,0) = f(0,0) / scale
        f(NDATA+1,0) = f(NDATA+1,0) / scale
        f(0,NDATA+1) = f(0,NDATA+1) / scale
        f(NDATA+1,NDATA+1) = f(NDATA+1,NDATA+1) / scale
        if (tag(1:1) .eq. 'e') then
           do i = 1, NDATA
              x1(i) = max (0d0, min (1d0, 1d0-x1(i)))
           end do   
           do i = 0, NDATA+1
              tmp = f(NDATA+1,i) 
              f(NDATA+1,i) = f(0,i)
              f(0,i) = tmp
           end do   
        end if
        if (tag(2:2) .eq. 'e') then
           do i = 1, NDATA
              x2(i) = max (0d0, min (1d0, 1d0-x2(i)))
           end do   
           do i = 0, NDATA+1
              tmp = f(i,NDATA+1)
              f(i,NDATA+1) = f(i,0)
              f(i,0) = tmp
           end do   
        end if      
        do i = 1, NDATA
           x1(i) = min (1d0, (x1(i)**(1d0/pwr) + 0.5d0/dble(NDATA))**pwr)
           x2(i) = min (1d0, (x2(i)**(1d0/pwr) + 0.5d0/dble(NDATA))**pwr)
        end do   
  end subroutine gethst

  subroutine puthst (file, tag, NDATA, f, x1, x2)
    integer, intent(in) :: file, NDATA
    character(len=2), intent(in) :: tag
    real(kind=double), dimension(NDATA) :: x1, x2
    real(kind=double), dimension(0:NDATA+1,0:NDATA+1) :: f
    integer :: i
    do i = 0, NDATA+1
       write (file, 1000) 'xa2h'//tag, i, NDATA+1
       call putcol (file, NDATA+2, f(0,i))
    end do
    write (file, 1001) 'xa21'//tag, NDATA
    call putcol (file, NDATA, x1)
    write (file, 1001) 'xa22'//tag, NDATA
    call putcol (file, NDATA, x2)
1000 format (6X, 'data (', A6, '(i,', I3, ',@ENERGY@,@ACC@,1),i=0,', I3, ') /')
1001 format (6X, 'data (', A6, '(i', ',@ENERGY@,@ACC@,1),i=1,', I3, ') /')
  end subroutine puthst

  subroutine putcol (file, NDATA, x)
    integer, intent(in) :: file, NDATA
    character(len=80) :: fmt
    real(kind=double), dimension(NDATA) :: x
    integer :: i, j
    do i = 1, NDATA - 5, 5
       write (file, 1001) (x(j),j=i,i+4)
1001    format (5X, '$', 5(E12.6,','))
    end do    
      if (i .eq. NDATA) then
         write (file, 1002) x(i)
 1002    format (5X, '$', E12.6, '/')
      else
         write (fmt, 1010) NDATA - i
 1010    format ('(5X,''$'',', I2, '(E12.6,'',''),E12.6,''/'')')
         write (file, fmt) (x(j),j=i,NDATA)
      end if
  end subroutine putcol
   
       
  end module int_routines

  program interp
    use kinds
    use int_routines
    
    integer, parameter :: NDATA = 20
    double precision x1(NDATA), x2(NDATA), f(0:NDATA+1,0:NDATA+1)
    double precision x(2,NDATA,NDATA), pwr, scale
    character(len=2), dimension(4) :: tags = &
       (/ 'ee', 'eg', 'ge', 'gg' /)
    integer :: i
    open (11, file = 'Parameters')
    scale = 0d0
    do i = 1, 4
       call gethst (tags(i), NDATA, x, f, x1, x2, pwr, scale)
       call puthst (11, tags(i), NDATA, f, x1, x2)
    end do
    write (11, 1000) pwr
1000 format (6X, 'data xa2pwr(@ENERGY@,@ACC@,1) / ', E12.6, ' /')
    close (11)
  end program interp
