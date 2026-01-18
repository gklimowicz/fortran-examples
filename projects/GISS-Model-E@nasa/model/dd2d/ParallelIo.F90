module ParallelIo_mod
!@sum Convenience wrapper around pario_nc.f
!@+
!@+ Allows the same function call to be used to define, write or read
!@+ a NetCDF variable.  This saves programming effort and eliminates
!@+ problems of consistency between subroutines that define, write or
!@+ read a set of arrays.

   use dd2d_utils, only : dist_grid
   use pario
   implicit none
   private

   public :: ParallelIo ! type and constructor
   public :: doVar      ! procedure and method

   type ParallelIo
      private
      type (dist_grid), pointer :: grid
      integer :: fileId
      character(len=:), allocatable :: compareString
   contains
      generic :: doVar => doVar_2d_real64
      generic :: doVar => doVar_3d_real64
      generic :: doVar => doVar_4d_real64
      generic :: doVar => doVar_5d_real64
      generic :: doVar => doVar_2d_integer
      generic :: doVar => doVar_3d_integer
      generic :: doVar => doVar_4d_integer
      generic :: doVar => doVar_2d_logical
      procedure :: doVar_2d_real64
      procedure :: doVar_3d_real64
      procedure :: doVar_4d_real64
      procedure :: doVar_5d_real64
      procedure :: doVar_2d_integer
      procedure :: doVar_3d_integer
      procedure :: doVar_4d_integer
      procedure :: doVar_2d_logical
   end type ParallelIo

   ! This type is used to force keyword association
   ! for optional arguments in argument lists.
   ! it is intentionally private
   type UnusableArgument
   end type UnusableArgument


   interface doVar


   !@sum subroutine doVar(handle, action, arr, varinfo, ...
   !@+     r4_on_disk, jdim, no_xdim, record, record1)
   !@+
   !@+ All-in-one subroutine that will defin, read or write an array in
   !@+ NetCDF, depending on the contents of action.  Equivalent to
   !@+ defvar/read_dist_/write_dist in pario_nc.f
   !@+
   !@var class (ParallelIo), intent(in) :: handle
   !@+     Instance of ParallelIo constructed (see Parallelio() below)
   !@+
   !@var action :: character(len=*)
   !@+    'define': Define a variable
   !@+    'read_dist': Read a distributed array
   !@+    'write_dist': Write a distributed array
   !@+    NOTE: There is no equivalent to read_data() or write_data() in pario_nc.f
   !@+
   !@var <type> :: arr(:,:,...)
   !@+     Array or scalar to write.
   !@+     <type> may be real*8, integer or logical
   !@+
   !@+     NOTE: If the desired type/dimension implementation of this
   !@+           interface does not yet exist, it should be added.
   !@+
   !@var character(*) :: varinfo
   !@+     String defining name of variable and its dimensions to
   !@+     define in NetCDF.
   !@+     Example: 't(im,jm,lm)'
   !@+
   !@+ ***** Keyword args: the remaining must be specfied by keyword, eg jdim=17
   !@+
   !@var logical, intent(in), optional :: r4_on_disk
   !@+     Used only when action=='define'
   !@+     Indicates the defined real variable should be a 4-byte float
   !@+     even though the passed array is 8-byte (which is what happens
   !@+     writing out diagnostic acc files).
   !@+
   !@var integer, intent(in), optional :: jdim = 2
   !@+     Used only when action=='write_dist'
   !@+     Specifies the index (starting from 1) of the LAST horizontal
   !@+     dimension.  If not specified, jdim=2; correct for model arrays
   !@+     like T(i,j,l).  To write an array dimensioned T(l,i,j) set jdim=3.
   !@+
   !@var logical, intent(in), optional :: no_xdim = .false.
   !@+     Used only when action=='read_dist'
   !@+     (WARNING: Negative logic; let has_xdim = .not. no_xdim)
   !@+     ?????
   !@+
   !@var integer, intent(in), optional :: record,record1
   !@+     Used only when action=='read_dist'
   !@+     ????????
      module procedure doVar_2d_real64
      module procedure doVar_3d_real64
      module procedure doVar_4d_real64
      module procedure doVar_5d_real64
      module procedure doVar_2d_integer
      module procedure doVar_3d_integer
      module procedure doVar_4d_integer
      module procedure doVar_2d_logical
   end interface doVar

   interface ParallelIo
   !@sum type(ParallelIo) ParallelIo(grid, fid)
   !@+ --------------------------------------------------
   !@+
   !@+ Constructor of ParallelIO object.
   !@+
   !@+ EXAMPLE:
   !@+     type (ParallelIo) :: handle
   !@+     handle = ParallelIo(grid, fid)
   !@+
   !@var type(dist_grid), intent(in) :: grid
   !@+     The grid on which the array exists.
   !@+
   !@var integer :: fid
   !@+     Open file handle to write to (obtained via par_open())
   !@+     
      module procedure newParallelIo
   end interface ParallelIo


contains


   function newParallelIo(grid, fid, compareString) result(handle)
      type (ParallelIo) :: handle
      type (dist_grid), target :: grid
      integer, intent(in) :: fid
      character(len=*), optional, intent(in) :: compareString
      
      handle%grid => grid
      handle%fileId = fid
      if (present(compareString)) then
         handle%compareString = compareString
      end if

   end function newParallelIo

#define _DECLARE_(arr) real(kind=real64), dimension(:,:) :: arr
#define _PROC_NAME_ doVar_2d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) real(kind=real64), dimension(:,:,:) :: arr
#define _PROC_NAME_ doVar_3d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) real(kind=real64), dimension(:,:,:,:) :: arr
#define _PROC_NAME_ doVar_4d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) real(kind=real64), dimension(:,:,:,:,:) :: arr
#define _PROC_NAME_ doVar_5d_real64
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) integer, dimension(:,:) :: arr
#define _PROC_NAME_ doVar_2d_integer
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) integer, dimension(:,:,:) :: arr
#define _PROC_NAME_ doVar_3d_integer
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) integer, dimension(:,:,:,:) :: arr
#define _PROC_NAME_ doVar_4d_integer
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_

#define _DECLARE_(arr) logical, dimension(:,:) :: arr
#define _PROC_NAME_ doVar_2d_logical
#include "do_generic.inc"
#undef _PROC_NAME_
#undef _DECLARE_


   ! Extract variable name from varInfo string
   ! Such strings are all of the form    <varName>[(dims)]
   function getVariableNameFrom(varInfo) result(variableName)
      character(len=:), allocatable :: variableName
      character(len=*), intent(in) :: varInfo

      integer :: idx
      
      idx = index(varInfo, '(')
      if (idx == 0) idx = len_trim(varInfo)
      variableName = trim(adjustl(varInfo(1:idx-1)))
      
   end function getVariableNameFrom


end module ParallelIo_mod

