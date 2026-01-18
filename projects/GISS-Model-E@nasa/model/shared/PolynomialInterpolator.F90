!------------------------------------------------------------------------------
module PolynomialInterpolator_mod
!------------------------------------------------------------------------------
!@sum Module that contains methods to perform polynomial interpolation
!@auth SSSO ASTG
  implicit none
  private

  public interpolator2D
  public interpolator3D
  
  type interpolator
    private
    real*8  :: errorEstimate
    contains
      procedure :: getErrorEstimate
      procedure :: setErrorEstimate
  end type interpolator

  type, extends(interpolator) :: interpolator2D
    private
    real*8, allocatable, dimension(:) :: xCoordinates
    real*8, allocatable, dimension(:) :: yCoordinates
    real*8, allocatable, dimension(:,:) :: tabulatedValues     
    real*8  :: x, y                        
    contains
      procedure :: interpolate2dlin
      procedure :: interpolate2dcub
  end type interpolator2D

  type, extends(interpolator) :: interpolator3D
    private
    real*8, allocatable, dimension(:) :: xCoordinates
    real*8, allocatable, dimension(:) :: yCoordinates
    real*8, allocatable, dimension(:) :: zCoordinates
    real*8, allocatable, dimension(:,:,:) :: tabulatedValues     
    real*8  :: x, y, z                     
    contains
      procedure :: interpolate3dlin
      procedure :: interpolate3dcub
  end type interpolator3D

  interface interpolator2D
    module procedure newInterpolator2D
  end interface

  interface interpolator3D
    module procedure newInterpolator3D
  end interface

  integer, parameter :: NLIN = 2 ! need two points to make a line
  integer, parameter :: NCUB = 4 ! need four points for a cubic polynomial
  integer, parameter :: UNDEFINED = -1000. 

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

  function newInterpolator2D(xCoordinates, yCoordinates, tabulatedValues) &
    result(i2D)
!@sum Construct 2D interpolation object
    type (interpolator2D) :: i2D
    real*8, dimension(:)   :: xCoordinates, yCoordinates
    real*8, dimension(:,:) :: tabulatedValues   

    allocate(i2D%xCoordinates(size(xCoordinates)))
    allocate(i2D%yCoordinates(size(yCoordinates)))
    allocate(i2D%tabulatedValues(size(tabulatedValues,1), &
                                 size(tabulatedValues,2)))
    i2D%xCoordinates = xCoordinates
    i2D%yCoordinates = yCoordinates
    i2D%tabulatedValues = tabulatedValues

  end function newInterpolator2D

  function newInterpolator3D(xCoordinates, yCoordinates, zCoordinates, &
    tabulatedValues) result(i3D)
!@sum Construct 3D interpolation object
    type (interpolator3D)    :: i3D 
    real*8, dimension(:)     :: xCoordinates, yCoordinates, zCoordinates
    real*8, dimension(:,:,:) :: tabulatedValues

    allocate(i3D%xCoordinates(size(xCoordinates)))
    allocate(i3D%yCoordinates(size(yCoordinates)))
    allocate(i3D%zCoordinates(size(zCoordinates)))
    allocate(i3D%tabulatedValues(size(tabulatedValues,1), & 
                                 size(tabulatedValues,2), &
                                 size(tabulatedValues,3)))
    i3D%xCoordinates = xCoordinates
    i3D%yCoordinates = yCoordinates
    i3D%zCoordinates = zCoordinates
    i3D%tabulatedValues = tabulatedValues

  end function newInterpolator3D

  function interpolate2dlin(this, x) result(y)
!@sum 2D linear interpolation
    class (interpolator2D), intent(in) :: this
    real*8 :: y             ! interpolant
    real*8 :: x(2)
!
    integer i, j
    integer iTemp, jTemp, iLoc, jLoc
    real*8 ym(NLIN), yn(NLIN)
    real*8 xTemp(NLIN), yTemp(NLIN)
    real*8 dy

    iLoc = bracketAt(this%xCoordinates, x(1))
    jLoc = bracketAt(this%yCoordinates, x(2))

    do i = 1, NLIN
      iTemp = i
      do j = 1, NLIN
        jTemp = j
        yn(j) = this%tabulatedValues(iLoc+jTemp-1, jLoc+iTemp-1)
        xTemp(j) = this%xCoordinates(iLoc+jTemp-1)
      enddo
      if (yn(1) == UNDEFINED) then
        ym(i) = UNDEFINED
        yTemp(i) = this%yCoordinates(jLoc+iTemp-1)
      else
        call interpolate(xTemp, yn, x(1), ym(i), dy)
        yTemp(i) = this%yCoordinates(jLoc+iTemp-1)
      endif
    enddo
    if (ym(1) == UNDEFINED) then
      y = UNDEFINED
    else
      call interpolate(yTemp, ym, x(2), y, dy)
    endif
    call this%setErrorEstimate(dy)

  end function interpolate2dlin

  function interpolate2dcub(this, x) result(y)
!@sum 2D cubic interpolation
    class (interpolator2D), intent(in) :: this
    real*8 :: y             ! interpolant
    real*8 :: x(2)
!
    integer i, j
    integer iTemp, jTemp, iLoc, jLoc
    real*8 ym(NCUB), yn(NCUB)
    real*8 xTemp(NCUB), yTemp(NCUB)
    real*8 dy

    iLoc = bracketAt(this%xCoordinates, x(1))
    jLoc = bracketAt(this%yCoordinates, x(2))

    do i = 1, NCUB
      iTemp = i-1
      do j =1, NCUB
        jTemp = j-1
        yn(j) = this%tabulatedValues(iLoc+jTemp-1, jLoc+iTemp-1)
        xTemp(j) = this%xCoordinates(iLoc+jTemp-1)
      enddo
      if ((yn(1) == UNDEFINED) .or. (yn(2) == UNDEFINED)) then
        ym(i) = UNDEFINED
        yTemp(i) = this%yCoordinates(jLoc+iTemp-1)
      else
        call interpolate(xTemp, yn, x(1), ym(i), dy)
        yTemp(i) = this%yCoordinates(jLoc+iTemp-1)
      endif
    enddo
    if ((ym(1) == UNDEFINED) .or. (ym(2) == UNDEFINED)) then
      y = UNDEFINED
    else
      call interpolate(yTemp, ym, x(2), y, dy)
    endif
    call this%setErrorEstimate(dy)

  end function interpolate2dcub

  function interpolate3dlin(this, x) result(y) 
!@sum 3D linear interpolation
    class (interpolator3D), intent(in) :: this
    real*8 :: y             ! interpolant
    real*8 :: x(3)
!
    integer i, j, k
    integer iTemp, jTemp, kTemp, iLoc, jLoc, kLoc
    real*8 ym(NLIN), yn(NLIN), yo(NLIN) 
    real*8 xTemp(NLIN), yTemp(NLIN), zTemp(NLIN) 
    real*8 dy

    iLoc = bracketAt(this%xCoordinates, x(1))
    jLoc = bracketAt(this%yCoordinates, x(2))
    kLoc = bracketAt(this%zCoordinates, x(3))

    do i = 1, NLIN
      kTemp = i
      zTemp(i) = this%zCoordinates(kLoc+kTemp-1)
      do k = 1, NLIN
        iTemp = k
        yTemp(k) = this%yCoordinates(jLoc+iTemp-1)  
        do j = 1, NLIN
          jTemp = j
          xTemp(j) = this%xCoordinates(iLoc+jTemp-1)
          yn(j) = this%tabulatedValues(iLoc+jTemp-1, jLoc+iTemp-1, kLoc+kTemp-1)
        enddo
        if (yn(1) == UNDEFINED) then
          ym(k) = UNDEFINED
        else
          call interpolate(xTemp, yn, x(1), ym(k), dy)
        endif
      enddo
      if (ym(1) == UNDEFINED) then
        yo(i) = UNDEFINED
      else
        call interpolate(yTemp, ym, x(2), yo(i), dy)
      endif
    enddo
    if (yo(2) == UNDEFINED)  then
      y = UNDEFINED
    else
      call interpolate(zTemp, yo, x(3), y, dy)
    endif
    call this%setErrorEstimate(dy)

  end function interpolate3dlin

  function interpolate3dcub(this, x) result(y)
!@sum 3D cubic interpolation
    class (interpolator3D), intent(in) :: this
    real*8 :: y             ! interpolant
    real*8 :: x(3)
!
    integer i, j, k
    integer iTemp, jTemp, kTemp, iLoc, jLoc, kLoc
    real*8 ym(NCUB), yn(NCUB), yo(NCUB) 
    real*8 xTemp(NCUB), yTemp(NCUB), zTemp(NCUB) 
    real*8 dy
 
    iLoc = bracketAt(this%xCoordinates, x(1))
    jLoc = bracketAt(this%yCoordinates, x(2))
    kLoc = bracketAt(this%zCoordinates, x(3))

    do i = 1, NCUB
      kTemp = i-1
      zTemp(i) = this%zCoordinates(kLoc+kTemp-1)
      do k = 1, NCUB
        iTemp = k-1
        yTemp(k) = this%yCoordinates(jLoc+iTemp-1)
        do j = 1,NCUB
          jTemp = j-1
          xTemp(j) = this%xCoordinates(iLoc+jTemp-1)
          yn(j) = this%tabulatedValues(iLoc+jTemp-1, jLoc+iTemp-1, kLoc+kTemp-1)
        enddo
        if ((yn(1) == UNDEFINED) .or. (yn(2) == UNDEFINED)) then
          ym(k) = UNDEFINED
        else
          call interpolate(xTemp, yn, x(1), ym(k), dy)
        endif
      enddo
      if ((ym(1) == UNDEFINED) .or. (ym(2) == UNDEFINED)) then
        yo(i) = UNDEFINED
      else
        call interpolate(yTemp, ym, x(2), yo(i), dy)
      endif
    enddo
    if ((yo(3) == UNDEFINED) .or. (yo(4) == UNDEFINED))  then
      y = UNDEFINED
    else
      call interpolate(zTemp, yo, x(3), y, dy)
    endif
    call this%setErrorEstimate(dy)

  end function interpolate3dcub

  subroutine setErrorEstimate(this, error)
    class (interpolator) :: this
    real*8,intent(in) :: error

    select type(this)
    type is (interpolator)
    class is (interpolator2d)
    class is (interpolator3d)
      this%errorEstimate = error
    class default
    end select

  end subroutine setErrorEstimate

  function getErrorEstimate(this) result (error)
    class (interpolator) :: this
    real*8 :: error
    error = this%errorEstimate
  end function getErrorEstimate

  subroutine interpolate(xValues, yValues, x, y, dy)
!@sum This is an implementation of Neville's algorithm used for polynomial
!@+ interpolation (e.g. http://en.wikipedia.org/wiki/Neville's_algorithm) 
    real*8 :: xValues(:)
    real*8 :: yValues(:)
    real*8, intent(in)  :: x      ! interpolated point
    real*8, intent(out) :: y      ! interpolant
    real*8, intent(out) :: dy     ! error
!
    integer i, j, n           ! indices
    integer ix                ! ix is closest index in data to x   
    real*8 xerr1, xerr2, xerr ! errors in x
    real*8 abDiff, dif, difx  ! differences
    real*8, allocatable, dimension(:) :: a, b       ! corrections

    n = size(xValues)
    allocate(a(n), b(n))

    ix = 1
    dif = abs(x-xValues(1))

    ! first iteration: initialize corrections and ix 
    do i = 1, n
      difx = abs(x-xValues(i))
      if (difx < dif) then
        ix = i 
        dif = difx
      end if
      a(i) = yValues(i)
      b(i) = yValues(i)
    end do

    y = yValues(ix)       ! first guess
    ix = ix - 1

    ! second iteration: recur to interpolant
    do j = 1, n-1
      do i = 1, n-j
        xerr1 = xValues(i) - x
        xerr2 = xValues(i+j) - x
        xerr  = xerr1 - xerr2
        abDiff = a(i+1) - b(i)
        ! exit if input x's are equal (within round off)
        if(xerr == 0.d0) call stop_model('failure in interpolate',255)
        xerr = abDiff/xerr
        a(i) = xerr1*xerr
        b(i) = xerr2*xerr
      end do
      ! update y using corrections
      if (2*ix < n-j) then
        dy = a(ix+1)
      else
        dy = b(ix)
        ix = ix - 1
      end if
      y = y + dy 
    end do

  end subroutine interpolate

  function bracketAt(values, x) result(i)
!@sum locates parameters of integration in lookup table
    real*8 :: values(:) ! data
    real*8, intent(in)  :: x     ! point to be bracketed
    integer :: i                 ! bracketed index
!
    integer il, im, iu ! lower, middle and upper limits
    integer  :: n     ! data size

    n = size(values)
    il = 0
    iu = n + 1
    do while (iu - il > 1)
      im = (iu + il) / 2
      if( (values(n) > values(1)) .eqv. (x > values(im)) )then
        il = im
      else
        iu = im
      end if
    end do
    i = il
    
  end function bracketAt

end module PolynomialInterpolator_mod
