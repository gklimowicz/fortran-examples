dnl 
dnl This file provides m4 macros for automatic generation
dnl of f90 interfaces
dnl
dnl Author: I. Aleinov
dnl Email: ialeinov@giss.nasa.gov
dnl
dnl                Description
dnl
dnl The fortran file which wants to utilize these macros should 
dnl start with including this file:
dnl
dnl include(f90_interface.m4)
dnl
dnl It should then define fortran types and array ranks it wants
dnl to loop over. For example, to loop over integer and real*8 types
dnl and over array ranks 1,2,3,4 one defines
dnl
dnl define(`m4_type_list',`integer,real*8')
dnl define(`m4_rank_list',`1,2,3,4')
dnl 
dnl The main function to use to loop over specified types and ranks
dnl is m4_loop(). One typically should pu it around the code
dnl taht needs to be duplicated. For example,
dnl 
dnl m4_loop(`
dnl      some code
dnl ')
dnl 
dnl will duplicate "some code" for each type and rank that we specified
dnl (2*4 = 8 times in our example).
dnl 
dnl The following macro definitions are set inside the m4_loop:
dnl 
dnl Loop variables:
dnl m4_type - current type ( abbreviated to
dnl                          i - integer
dnl                          r4 - real*4
dnl                          r8 - real*8
dnl                          c - character
dnl                          l - logical
dnl 
dnl m4_rank - current runk (integer number)
dnl 
dnl Macros for array descriptions:
dnl 
dnl m4_typedecl - fortran declaration for current type, i.e.
dnl               integer,real*4,real*8,character or logical
dnl m4_ext - "extent" (expression :,:,:,... corresponding to current rank)
dnl m4_dims(dims) - "dimensions" (expression of a kind:
dnl              dims(1,1):dims(2,1), dims(1,2):dims(2,2), ...
dnl              where dims(2,:) is a pre-defined array with lower
dnl              and upper bounds
dnl m4_shape(a) - same as above but will use lower and upper bounds
dnl              of an existion array "a" (i.e. will create array of the
dnl              same "shape")
dnl 
dnl "Suffix":
dnl m4_suffix - expression of a kind "m4_type()_m4_rank()" which can
dnl              be appended to identifiers to create unique names
dnl              for particular type and kind (for example, when one
dnl              defines module procedures for an interface block).
dnl 
dnl Explicit DO loops:
dnl
dnl m4_do(dims) - start DO loop over all indices, lover and upper bounds
dnl               are described by the array dims (defined as above)
dnl m4_enddo    - "END DO" for m4_do(dims)
dnl m4_ind      - generalized DO index, ie. the expression
dnl               i1,i2,i3,... corresponding to current rank
dnl 
dnl The following macros are "hacks" to deal with rank=0, i.e. when
dnl one wants to create an interface which handles both arrays and
dnl scalars:
dnl 
dnl m4_brkts(code) - will expand to "(code)" if "code" is not empty.
dnl                Otherwise it will expand to an empty string.
dnl                For rank=0 m4_ext, m4_dims, m4_shape, m4_ind
dnl                will become an empty string. So one has to use
dnl                them with m4_brkts() to omit brackets for scalars.
dnl m4_arr_only(m4_rank,code) - will expand to "code" if m4_rank
dnl                is non-zero. Otherwise it will expand to an empty 
dnl                line. Useful when some code has to be omitted
dnl                for scalars. For example, one doesn't need DO 
dnl                indices for scalars:
dnl                m4_arr_only(m4_rank,integer :: m4_ind)
dnl 
dnl Other macros (may become useful in some special cases):
dnl 
dnl m4_loop_types(`code') - same as m4_loop(`code') but loop only
dnl                over types
dnl m4_loop_ranks(`code') - same as m4_loop(`code') but loop only
dnl                over ranks
dnl 
dnl m4_ext_r, m4_dims_r, m4_shape_r, m4_do_r, m4_enddo_r, m4_ind_r
dnl                - same as above but take an extra argument (the first)
dnl                which specifies the rank. These may be useful if one
dnl                needs to deal with arrays of a rank which differs
dnl                from m4_rank (for example rank=m4_rank-1)
dnl 
dnl m4_foreach(`index', (list),`code') - generic loop. Will loop over
dnl                all variables in the "list" (coma-separated) repeating
dnl                "code" with macro "index" being expanded to the current
dnl                variable from the "list".
dnl 
dnl 
dnl               Examples:
dnl 
dnl Here is a most trivial example of something useful. It creates
dnl an interface for integers and real*8 for ranks 1,2,6 . The function
dnl takes "a" as input, allocates temporary array "c" (of the same
dnl shape), performs some trivial computations and returns result in "b";
dnl 
dnl include(f90_interface.m4)
dnl define(`m4_rank_list',`1,2,6')
dnl define(`m4_type_list',`integer,real*8')
dnl 
dnl       module aaa
dnl 
dnl       interface mycopy
dnl         m4_loop(`
dnl         module procedure mycopy_`'m4_suffix
dnl         ')
dnl       end interface mycopy
dnl 
dnl       contains
dnl 
dnl       m4_loop(`
dnl       subroutine mycopy_`'m4_suffix`'(a,b)
dnl       m4_typedecl, intent(in) :: a(m4_ext)
dnl       m4_typedecl, intent(out) :: b(m4_ext)
dnl       m4_typedecl :: c(m4_shape(a))
dnl       c = a
dnl       c = c + a
dnl       b = c
dnl       end subroutine mycopy_`'m4_suffix
dnl       ')
dnl 
dnl       end module aaa
dnl 
divert(-1)

dnl define default lists of types and ranks
dnl (can be redefined in user's file)

define(`m4_type_list',`integer,real*8')
define(`m4_rank_list',`1,3')

dnl define some useful functions

define(`_for',`ifelse($#,0,``$0'',`ifelse(eval($2<=$3),1,
    `pushdef(`$1',$2)$4`'popdef(`$1')$0(`$1',incr($2),$3,`$4')')')')

define(`m4_foreach', `pushdef(`$1')__foreach($@)popdef(`$1')')
dnl define(`_foreach', `pushdef(`$1')__foreach($@)popdef(`$1')')
define(`_arg1', `$1')
define(`__foreach', `ifelse(`$2', `()', `',
  	`define(`$1', _arg1$2)$3`'$0(`$1', (shift$2), `$3')')')

dnl define lookup tables and typical patterns

define(`stype[integer]',`i')
define(`stype[real*4]',`r4')
define(`stype[real*8]',`r8')
define(`stype[character]',`c')
define(`stype[logical]',`l')

define(`decl[i]',`integer')
define(`decl[r4]',`real*4')
define(`decl[r8]',`real*8')
define(`decl[c]',`character')
define(`decl[l]',`logical')

define(`m4_ext_r',`_for(`x',1,eval(`$1'),`:ifelse(x,eval(`$1'),,`,')')')

define(`m4_dims_r',`_for(`x',1,eval(`$1'),`
     &  $2(1,x):$2(2,x)ifelse(x,eval(`$1'),,`,')')')

dnl define(`m4_shape_r',`_for(`x',1,$1,`
dnl      &  size($2,x)ifelse(x,$1,,`,')')')
define(`m4_shape_r',`_for(`x',1,eval(`$1'),`
     &  lbound($2,x):ubound($2,x)ifelse(x,eval(`$1'),,`,')')')

define(`m4_do_r',`_for(`x',1,eval(`$1'),`
      do i`'x=$2(1,x),$2(2,x)')')

define(`m4_enddo_r',`_for(`x',1,eval(`$1'),`
      enddo')')

define(`m4_ind_r',`_for(`x',1,eval(`$1'),`i`'x`'ifelse(x,eval(`$1'),,`,')')')

define(`m4_brkts',`ifelse(`$1',`',`',`($@)')')

define(`m4_arr_only',`ifelse(eval(`$1'),`0',`',`shift($@)')')

dnl redefine type list to short type list
define(`m4_stype_list',
 `shift(m4_foreach(`_t', (m4_type_list), `,defn(stype[_t])'))')
dnl  define(`m4_stype_list',
dnl  `m4_foreach(`_t', (m4_type_list), `,YYYdefn(stype[_t])')')

dnl define main loops over types and ranks

define(`m4_loop_types',`m4_foreach(`m4_type', (m4_stype_list),`$1')')

define(`m4_loop_ranks',`m4_foreach(`m4_rank', (m4_rank_list),`$1')')

define(`m4_loop',`m4_loop_types(`m4_loop_ranks(`$1')')')

dnl define shorter definitions for default rank
define(`m4_ext',`m4_ext_r(m4_rank)')dnl
define(`m4_dims',`m4_dims_r(m4_rank,$1)')dnl
define(`m4_shape',`m4_shape_r(m4_rank,$1)')dnl
define(`m4_do',`m4_do_r(m4_rank,$1)')dnl
define(`m4_enddo',`m4_enddo_r(m4_rank)')dnl
define(`m4_ind',`m4_ind_r(m4_rank)')dnl

define(`m4_typedecl',`defn(decl[m4_type])')dnl
define(`m4_suffix',`m4_type()_`'m4_rank()')dnl

divert(0)
dnl
dnl The following lines will be included at the top of automatically
dnl generated file to prevent users from doing something silly
dnl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!                   DO NOT EDIT THIS FILE !!!                  !
!                                                              !
! This file is generated automatically by m4 preprocessor.     !
! If you need to modyfy it plese edit corresponding *.m4f file.!
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

