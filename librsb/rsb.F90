!> @file.
!! @brief Header file automatically generated from <rsb.h>, offering ISO-C-BINDING interfaces to <rsb.h>'s functions.
!! Defines \c MODULE \c rsb.
!! For examples of usage, see Fortran examples in \ref rsb_doc_examples.
!! The official documentation is that of <rsb.h>.
!! Make sure you are using a modern Fortran compiler.
      
!DEC$IF .NOT. DEFINED (RSB_FORTRAN_HEADER)
!DEC$DEFINE RSB_FORTRAN_HEADER
      
      MODULE rsb
         USE ISO_C_BINDING, ONLY: C_INT,C_INT64_T,C_PTR,C_NULL_PTR,C_SIGNED_CHAR
      
#ifdef RSB_WANT_LONG_IDX_TYPE
      INTEGER,PARAMETER :: RSB_IDX_KIND=8
#define C_RSB_INT_KND_ C_INT64_T
#else
      INTEGER,PARAMETER :: RSB_IDX_KIND=4
#define C_RSB_INT_KND_ C_INT
#endif
      
!> ISO C BINDING interface to ::rsb_strerror_r.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_strerror_r&
        &(errval,buf,buflen)&
        &BIND(c,NAME = 'rsb_strerror_r')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: errval
       CHARACTER(C_CHAR), DIMENSION(*) :: buf ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       INTEGER(C_SIZE_T), VALUE  :: buflen
       END FUNCTION rsb_strerror_r
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_perror.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_perror&
        &(stream,errval)&
        &BIND(c,NAME = 'rsb_perror')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: stream ! A numerical type
       INTEGER(C_INT), VALUE  :: errval
       END FUNCTION rsb_perror
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_lib_init.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_lib_init&
        &(iop)&
        &BIND(c,NAME = 'rsb_lib_init')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: iop ! C_NULL_PTR is a safe value. Please consult the rsb.h documentation for other options.
       END FUNCTION rsb_lib_init
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_lib_reinit.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_lib_reinit&
        &(iop)&
        &BIND(c,NAME = 'rsb_lib_reinit')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: iop ! C_NULL_PTR is a safe value. Please consult the rsb.h documentation for other options.
       END FUNCTION rsb_lib_reinit
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_lib_set_opt_str.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_lib_set_opt_str&
        &(opnp,opvp)&
        &BIND(c,NAME = 'rsb_lib_set_opt_str')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: opnp ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       CHARACTER(C_CHAR), DIMENSION(*) :: opvp ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       END FUNCTION rsb_lib_set_opt_str
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_lib_set_opt.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_lib_set_opt&
        &(iof,iop)&
        &BIND(c,NAME = 'rsb_lib_set_opt')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: iof
       TYPE(C_PTR),VALUE :: iop ! C_NULL_PTR is a safe value. Please consult the rsb.h documentation for other options.
       END FUNCTION rsb_lib_set_opt
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_lib_get_opt.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_lib_get_opt&
        &(iof,iop)&
        &BIND(c,NAME = 'rsb_lib_get_opt')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: iof
       TYPE(C_PTR),VALUE :: iop ! C_NULL_PTR is a safe value. Please consult the rsb.h documentation for other options.
       END FUNCTION rsb_lib_get_opt
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_lib_exit.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_lib_exit&
        &(iop)&
        &BIND(c,NAME = 'rsb_lib_exit')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: iop ! C_NULL_PTR is a safe value. Please consult the rsb.h documentation for other options.
       END FUNCTION rsb_lib_exit
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_coo_begin.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_alloc_from_coo_begin&
        &(nnzA,typecode,nrA,ncA,flagsA,errvalp)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_coo_begin')
       USE ISO_C_BINDING
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_mtx_alloc_from_coo_begin
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_coo_end.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_alloc_from_coo_end&
        &(mtxApp)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_coo_end')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxApp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       END FUNCTION rsb_mtx_alloc_from_coo_end
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_csr_const.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_alloc_from_csr_const&
        &(VA,RP,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_csr_const')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: RP ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_INT), VALUE  :: brA
       INTEGER(C_INT), VALUE  :: bcA
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_mtx_alloc_from_csr_const
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_csc_const.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_alloc_from_csc_const&
        &(VA,IA,CP,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_csc_const')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR), VALUE  :: CP ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_INT), VALUE  :: brA
       INTEGER(C_INT), VALUE  :: bcA
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_mtx_alloc_from_csc_const
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_csr_inplace.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_alloc_from_csr_inplace&
        &(VA,RP,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_csr_inplace')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: RP ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_INT), VALUE  :: brA
       INTEGER(C_INT), VALUE  :: bcA
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_mtx_alloc_from_csr_inplace
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_coo_const.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_alloc_from_coo_const&
        &(VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_coo_const')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_INT), VALUE  :: brA
       INTEGER(C_INT), VALUE  :: bcA
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_mtx_alloc_from_coo_const
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_alloc_from_coo_inplace.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_alloc_from_coo_inplace&
        &(VA,IA,JA,nnzA,typecode,nrA,ncA,brA,bcA,flagsA,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_mtx_alloc_from_coo_inplace')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_INT), VALUE  :: brA
       INTEGER(C_INT), VALUE  :: bcA
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_mtx_alloc_from_coo_inplace
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_clone.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_clone&
        &(mtxBpp,typecode,transA,alphap,mtxAp,flags)&
        &BIND(c,NAME = 'rsb_mtx_clone')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxBpp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_clone
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_free.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_mtx_free&
        &(mtxAp)&
        &BIND(c,NAME = 'rsb_mtx_free')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       END FUNCTION rsb_mtx_free
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_nrm.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_nrm&
        &(mtxAp,Np,flags)&
        &BIND(c,NAME = 'rsb_mtx_get_nrm')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: Np ! A single variable of the same numerical type of the matrix.
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_extff_t
       END FUNCTION rsb_mtx_get_nrm
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_vec.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_vec&
        &(mtxAp,Dp,flags)&
        &BIND(c,NAME = 'rsb_mtx_get_vec')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: Dp ! A single variable of the same numerical type of the matrix.
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_extff_t
       END FUNCTION rsb_mtx_get_vec
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_rndr.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_rndr&
        &(filename,mtxAp,pmWidth,pmHeight,rflags)&
        &BIND(c,NAME = 'rsb_mtx_rndr')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: pmWidth
       INTEGER(C_RSB_INT_KND_), VALUE  :: pmHeight
       INTEGER(C_INT), VALUE  :: rflags !> ISO C BINDING interface to ::rsb_marf_t
       END FUNCTION rsb_mtx_rndr
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_file_mtx_rndr.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_file_mtx_rndr&
        &(pmp,filename,pmlWidth,pmWidth,pmHeight,rflags)&
        &BIND(c,NAME = 'rsb_file_mtx_rndr')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: pmp ! A numerical type
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       INTEGER(C_RSB_INT_KND_), VALUE  :: pmlWidth
       INTEGER(C_RSB_INT_KND_), VALUE  :: pmWidth
       INTEGER(C_RSB_INT_KND_), VALUE  :: pmHeight
       INTEGER(C_INT), VALUE  :: rflags !> ISO C BINDING interface to ::rsb_marf_t
       END FUNCTION rsb_file_mtx_rndr
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_spmv.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_spmv&
        &(transA,alphap,mtxAp,Xp,incX,betap,Yp,incY)&
        &BIND(c,NAME = 'rsb_spmv')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: Xp ! A single variable of the same numerical type of the matrix.
       INTEGER(C_RSB_INT_KND_), VALUE  :: incX
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: Yp ! A single variable of the same numerical type of the matrix.
       INTEGER(C_RSB_INT_KND_), VALUE  :: incY
       END FUNCTION rsb_spmv
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_spmm.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_spmm&
        &(transA,alphap,mtxAp,nrhs,order,Bp,ldB,betap,Cp,ldC)&
        &BIND(c,NAME = 'rsb_spmm')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrhs
       INTEGER(C_INT), VALUE  :: order !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR), VALUE  :: Bp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldB
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: Cp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldC
       END FUNCTION rsb_spmm
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_spsv.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_spsv&
        &(transT,alphap,mtxTp,Xp,incX,Yp,incY)&
        &BIND(c,NAME = 'rsb_spsv')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: transT
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxTp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: Xp ! A single variable of the same numerical type of the matrix.
       INTEGER(C_RSB_INT_KND_), VALUE  :: incX
       TYPE(C_PTR),VALUE :: Yp ! A single variable of the same numerical type of the matrix.
       INTEGER(C_RSB_INT_KND_), VALUE  :: incY
       END FUNCTION rsb_spsv
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_spsm.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_spsm&
        &(transT,alphap,mtxTp,nrhs,order,betap,Bp,ldB,Cp,ldC)&
        &BIND(c,NAME = 'rsb_spsm')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: transT
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxTp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrhs
       INTEGER(C_INT), VALUE  :: order !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: Bp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldB
       TYPE(C_PTR), VALUE  :: Cp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldC
       END FUNCTION rsb_spsm
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_add_to_dense.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_add_to_dense&
        &(alphap,mtxAp,ldB,nrB,ncB,rowmajorB,Bp)&
        &BIND(c,NAME = 'rsb_mtx_add_to_dense')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldB
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrB
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncB
       INTEGER(C_INT), VALUE  :: rowmajorB !> ISO C BINDING interface to ::rsb_bool_t
       TYPE(C_PTR), VALUE  :: Bp ! A numerical type
       END FUNCTION rsb_mtx_add_to_dense
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_sppsp.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_sppsp&
        &(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_sppsp')
       USE ISO_C_BINDING
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: transB
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxBp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_sppsp
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_spmsp.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_spmsp&
        &(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,errvalp&
      &)&
        &BIND(c,NAME = 'rsb_spmsp')
       USE ISO_C_BINDING
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: transB
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxBp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_spmsp
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_spmsp_to_dense.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_spmsp_to_dense&
        &(typecode,transA,alphap,mtxAp,transB,betap,mtxBp,ldC&
      &,nrC,ncC,rowmajorC,Cp)&
        &BIND(c,NAME = 'rsb_spmsp_to_dense')
       USE ISO_C_BINDING
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: transB
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxBp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldC
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrC
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncC
       INTEGER(C_INT), VALUE  :: rowmajorC !> ISO C BINDING interface to ::rsb_bool_t
       TYPE(C_PTR), VALUE  :: Cp ! A numerical type
       END FUNCTION rsb_spmsp_to_dense
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_switch_to_coo.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_switch_to_coo&
        &(mtxAp,VAp,IAp,JAp,flags)&
        &BIND(c,NAME = 'rsb_mtx_switch_to_coo')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VAp ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IAp ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JAp ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_switch_to_coo
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_switch_to_csr.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_switch_to_csr&
        &(mtxAp,VAp,IAp,JAp,flags)&
        &BIND(c,NAME = 'rsb_mtx_switch_to_csr')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VAp ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IAp ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JAp ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_switch_to_csr
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_coo.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_coo&
        &(mtxAp,VA,IA,JA,flags)&
        &BIND(c,NAME = 'rsb_mtx_get_coo')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_get_coo
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_csr.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_csr&
        &(typecode,mtxAp,VA,RP,JA,flags)&
        &BIND(c,NAME = 'rsb_mtx_get_csr')
       USE ISO_C_BINDING
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: RP ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_get_csr
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_rows_sparse.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_rows_sparse&
        &(transA,alphap,mtxAp,VA,IA,JA,frA,lrA,rnzp,flags)&
        &BIND(c,NAME = 'rsb_mtx_get_rows_sparse')
       USE ISO_C_BINDING
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: frA
       INTEGER(C_RSB_INT_KND_), VALUE  :: lrA
       TYPE(C_PTR), VALUE  :: rnzp ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_get_rows_sparse
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_coo_block.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_coo_block&
        &(mtxAp,VA,IA,JA,frA,lrA,fcA,lcA,IREN,JREN,rnzp,flags&
      &)&
        &BIND(c,NAME = 'rsb_mtx_get_coo_block')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: frA
       INTEGER(C_RSB_INT_KND_), VALUE  :: lrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: fcA
       INTEGER(C_RSB_INT_KND_), VALUE  :: lcA
       TYPE(C_PTR), VALUE  :: IREN ! INTEGER(C_INT)
       TYPE(C_PTR), VALUE  :: JREN ! INTEGER(C_INT)
       TYPE(C_PTR), VALUE  :: rnzp ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_get_coo_block
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_info.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_info&
        &(mtxAp,miflags,minfop)&
        &BIND(c,NAME = 'rsb_mtx_get_info')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: miflags
       TYPE(C_PTR), VALUE  :: minfop ! A numerical type
       END FUNCTION rsb_mtx_get_info
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_info_str.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_info_str&
        &(mtxAp,mis,minfop,buflen)&
        &BIND(c,NAME = 'rsb_mtx_get_info_str')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       CHARACTER(C_CHAR), DIMENSION(*) :: mis ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       TYPE(C_PTR), VALUE  :: minfop ! A numerical type
       INTEGER(C_SIZE_T), VALUE  :: buflen
       END FUNCTION rsb_mtx_get_info_str
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_upd_vals.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_upd_vals&
        &(mtxAp,elop_flags,omegap)&
        &BIND(c,NAME = 'rsb_mtx_upd_vals')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: elop_flags
       TYPE(C_PTR), VALUE  :: omegap ! A numerical type
       END FUNCTION rsb_mtx_upd_vals
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_prec.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_prec&
        &(opdp,mtxAp,prec_flags,ipdp)&
        &BIND(c,NAME = 'rsb_mtx_get_prec')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: opdp ! A numerical type
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_INT), VALUE  :: prec_flags
       TYPE(C_PTR), VALUE  :: ipdp ! A numerical type
       END FUNCTION rsb_mtx_get_prec
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_set_vals.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_set_vals&
        &(mtxAp,VA,IA,JA,nnz,flags)&
        &BIND(c,NAME = 'rsb_mtx_set_vals')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnz
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_set_vals
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_mtx_get_vals.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_mtx_get_vals&
        &(mtxAp,VA,IA,JA,nnz,flags)&
        &BIND(c,NAME = 'rsb_mtx_get_vals')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnz
       INTEGER(C_INT), VALUE  :: flags !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_mtx_get_vals
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_tune_spmm.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_tune_spmm&
        &(mtxOpp,sfp,tnp,maxr,maxt,transA,alphap,mtxAp,nrhs,order&
      &,Bp,ldB,betap,Cp,ldC)&
        &BIND(c,NAME = 'rsb_tune_spmm')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxOpp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR), VALUE  :: sfp ! REAL*8
       TYPE(C_PTR), VALUE  :: tnp ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: maxr
       REAL(C_DOUBLE), VALUE  :: maxt
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrhs
       INTEGER(C_INT), VALUE  :: order !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR), VALUE  :: Bp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldB
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: Cp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldC
       END FUNCTION rsb_tune_spmm
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_tune_spsm.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_tune_spsm&
        &(mtxOpp,sfp,tnp,maxr,maxt,transA,alphap,mtxAp,nrhs,order&
      &,Bp,ldB,betap,Cp,ldC)&
        &BIND(c,NAME = 'rsb_tune_spsm')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxOpp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       TYPE(C_PTR), VALUE  :: sfp ! REAL*8
       TYPE(C_PTR), VALUE  :: tnp ! INTEGER(C_INT)
       INTEGER(C_INT), VALUE  :: maxr
       REAL(C_DOUBLE), VALUE  :: maxt
       INTEGER(C_INT), VALUE  :: transA
       TYPE(C_PTR),VALUE :: alphap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrhs
       INTEGER(C_INT), VALUE  :: order !> ISO C BINDING interface to ::rsb_flags_t
       TYPE(C_PTR), VALUE  :: Bp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldB
       TYPE(C_PTR),VALUE :: betap ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: Cp ! A numerical type
       INTEGER(C_RSB_INT_KND_), VALUE  :: ldC
       END FUNCTION rsb_tune_spsm
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_psblas_trans_to_rsb_trans.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_psblas_trans_to_rsb_trans&
        &(psbtrans)&
        &BIND(c,NAME = 'rsb_psblas_trans_to_rsb_trans')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), VALUE  :: psbtrans
       END FUNCTION rsb_psblas_trans_to_rsb_trans
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_file_mtx_save.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_file_mtx_save&
        &(mtxAp,filename)&
        &BIND(c,NAME = 'rsb_file_mtx_save')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: mtxAp ! A matrix pointer variable: (TYPE(C_PTR),TARGET)
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       END FUNCTION rsb_file_mtx_save
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_file_mtx_load.
      INTERFACE
       TYPE(C_PTR) FUNCTION &
        &rsb_file_mtx_load&
        &(filename,flagsA,typecode,errvalp)&
        &BIND(c,NAME = 'rsb_file_mtx_load')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       TYPE(C_PTR),VALUE :: errvalp ! INTEGER(C_INT)
       END FUNCTION rsb_file_mtx_load
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_file_vec_load.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_file_vec_load&
        &(filename,typecode,Yp,yvlp)&
        &BIND(c,NAME = 'rsb_file_vec_load')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       TYPE(C_PTR),VALUE :: Yp ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR), VALUE  :: yvlp ! INTEGER(C_INT)
       END FUNCTION rsb_file_vec_load
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_file_vec_save.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_file_vec_save&
        &(filename,typecode,Yp,yvl)&
        &BIND(c,NAME = 'rsb_file_vec_save')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       TYPE(C_PTR),VALUE :: Yp ! A single variable of the same numerical type of the matrix.
       INTEGER(C_RSB_INT_KND_), VALUE  :: yvl
       END FUNCTION rsb_file_vec_save
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_file_mtx_get_dims.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_file_mtx_get_dims&
        &(filename,nrp,ncp,nzp,flagsp)&
        &BIND(c,NAME = 'rsb_file_mtx_get_dims')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       TYPE(C_PTR), VALUE  :: nrp ! INTEGER(C_INT)
       TYPE(C_PTR), VALUE  :: ncp ! INTEGER(C_INT)
       TYPE(C_PTR), VALUE  :: nzp ! INTEGER(C_INT)
       TYPE(C_PTR), VALUE  :: flagsp ! INTEGER(C_INT)
       END FUNCTION rsb_file_mtx_get_dims
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_coo_sort.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_coo_sort&
        &(VA,IA,JA,nnzA,nrA,ncA,typecode,flagsA)&
        &BIND(c,NAME = 'rsb_coo_sort')
       USE ISO_C_BINDING
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_coo_sort
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_coo_cleanup.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_coo_cleanup&
        &(nnzp,VA,IA,JA,nnzA,nrA,ncA,typecode,flagsA)&
        &BIND(c,NAME = 'rsb_coo_cleanup')
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE  :: nnzp ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: VA ! A single variable of the same numerical type of the matrix.
       TYPE(C_PTR),VALUE :: IA ! INTEGER(C_INT)
       TYPE(C_PTR),VALUE :: JA ! INTEGER(C_INT)
       INTEGER(C_RSB_INT_KND_), VALUE  :: nnzA
       INTEGER(C_RSB_INT_KND_), VALUE  :: nrA
       INTEGER(C_RSB_INT_KND_), VALUE  :: ncA
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       INTEGER(C_INT), VALUE  :: flagsA !> ISO C BINDING interface to ::rsb_flags_t
       END FUNCTION rsb_coo_cleanup
      END INTERFACE
      
!> ISO C BINDING interface to ::rsb_time.
      INTERFACE
       REAL(C_DOUBLE) FUNCTION &
        &rsb_time&
        &()&
        &BIND(c,NAME = 'rsb_time')
       USE ISO_C_BINDING
       END FUNCTION rsb_time
      END INTERFACE
      
!DEC$ENDIF
      
      
!> ISO C BINDING interface to ::rsb_blas_file_mtx_load.
      INTERFACE
       INTEGER(C_INT) FUNCTION &
        &rsb_blas_file_mtx_load&
        &(filename,typecode)&
        &BIND(c,NAME = 'rsb_blas_file_mtx_load')
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*) :: filename ! C interoperable string, e.g. 'filename'//C_NULL_CHAR or C_NULL_CHAR .
       INTEGER(C_SIGNED_CHAR), VALUE  :: typecode
       END FUNCTION rsb_blas_file_mtx_load
      END INTERFACE
! Error values 
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_NO_ERROR&
            & = -INT(Z"0000",C_INT) !< See #RSB_ERR_NO_ERROR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_GENERIC_ERROR&
            & = -INT(Z"0001",C_INT) !< See #RSB_ERR_GENERIC_ERROR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_UNSUPPORTED_OPERATION&
            & = -INT(Z"0002",C_INT) !< See #RSB_ERR_UNSUPPORTED_OPERATION.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_UNSUPPORTED_TYPE&
            & = -INT(Z"0004",C_INT) !< See #RSB_ERR_UNSUPPORTED_TYPE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_UNSUPPORTED_FORMAT&
            & = -INT(Z"0008",C_INT) !< See #RSB_ERR_UNSUPPORTED_FORMAT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_INTERNAL_ERROR&
            & = -INT(Z"0010",C_INT) !< See #RSB_ERR_INTERNAL_ERROR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_BADARGS&
            & = -INT(Z"0020",C_INT) !< See #RSB_ERR_BADARGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_ENOMEM&
            & = -INT(Z"0040",C_INT) !< See #RSB_ERR_ENOMEM.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_UNIMPLEMENTED_YET&
            & = -INT(Z"0100",C_INT) !< See #RSB_ERR_UNIMPLEMENTED_YET.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_LIMITS&
            & = -INT(Z"0200",C_INT) !< See #RSB_ERR_LIMITS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_UNSUPPORTED_FEATURE&
            & = -INT(Z"0400",C_INT) !< See #RSB_ERR_UNSUPPORTED_FEATURE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_NO_USER_CONFIGURATION&
            & = -INT(Z"0800",C_INT) !< See #RSB_ERR_NO_USER_CONFIGURATION.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_CORRUPT_INPUT_DATA&
            & = -INT(Z"01000",C_INT) !< See #RSB_ERR_CORRUPT_INPUT_DATA.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_FAILED_MEMHIER_DETECTION&
            & = -INT(Z"02000",C_INT) !< See #RSB_ERR_FAILED_MEMHIER_DETECTION.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS&
            & = -INT(Z"04000",C_INT) !< See #RSB_ERR_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT&
            & = -INT(Z"08000",C_INT) !< See #RSB_ERR_NO_STREAM_OUTPUT_CONFIGURED_OUT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_INVALID_NUMERICAL_DATA&
            & = -INT(Z"010000",C_INT) !< See #RSB_ERR_INVALID_NUMERICAL_DATA.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_MEMORY_LEAK&
            & = -INT(Z"020000",C_INT) !< See #RSB_ERR_MEMORY_LEAK.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ERR_ELEMENT_NOT_FOUND&
            & = -INT(Z"040000000",C_INT) !< See #RSB_ERR_ELEMENT_NOT_FOUND.
! Matrix flags values 
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_NOFLAGS&
            & = INT(Z"0000000",C_INT) !< See #RSB_FLAG_NOFLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_FORTRAN_INDICES_INTERFACE&
            & = INT(Z"0000001",C_INT) !< See #RSB_FLAG_FORTRAN_INDICES_INTERFACE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_C_INDICES_INTERFACE&
            & = INT(Z"0000000",C_INT) !< See #RSB_FLAG_C_INDICES_INTERFACE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_USE_HALFWORD_INDICES&
            & = INT(Z"0000002",C_INT) !< See #RSB_FLAG_USE_HALFWORD_INDICES.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_WANT_ROW_MAJOR_ORDER&
            & = INT(Z"0000000",C_INT) !< See #RSB_FLAG_WANT_ROW_MAJOR_ORDER.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_WANT_COLUMN_MAJOR_ORDER&
            & = INT(Z"04000000",C_INT) !< See #RSB_FLAG_WANT_COLUMN_MAJOR_ORDER.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_SORTED_INPUT&
            & = INT(Z"0000004",C_INT) !< See #RSB_FLAG_SORTED_INPUT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_TRIANGULAR&
            & = INT(Z"0000008",C_INT) !< See #RSB_FLAG_TRIANGULAR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_LOWER&
            & = INT(Z"0000010",C_INT) !< See #RSB_FLAG_LOWER.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_UPPER&
            & = INT(Z"0000020",C_INT) !< See #RSB_FLAG_UPPER.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_UNIT_DIAG_IMPLICIT&
            & = INT(Z"0000040",C_INT) !< See #RSB_FLAG_UNIT_DIAG_IMPLICIT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_WANT_COO_STORAGE&
            & = INT(Z"0000100",C_INT) !< See #RSB_FLAG_WANT_COO_STORAGE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DUPLICATES_KEEP_LAST&
            & = INT(Z"0000000",C_INT) !< See #RSB_FLAG_DUPLICATES_KEEP_LAST.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DUPLICATES_DEFAULT_HANDLE&
            & = INT(Z"0000000",C_INT) !< See #RSB_FLAG_DUPLICATES_DEFAULT_HANDLE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DUPLICATES_SUM&
            & = INT(Z"0000200",C_INT) !< See #RSB_FLAG_DUPLICATES_SUM.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DISCARD_ZEROS&
            & = INT(Z"0000400",C_INT) !< See #RSB_FLAG_DISCARD_ZEROS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_QUAD_PARTITIONING&
            & = INT(Z"0002000",C_INT) !< See #RSB_FLAG_QUAD_PARTITIONING.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_WANT_BCSS_STORAGE&
            & = INT(Z"0004000",C_INT) !< See #RSB_FLAG_WANT_BCSS_STORAGE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS&
            & = INT(Z"0040000",C_INT) !< See #RSB_FLAG_ASSEMBLED_IN_COO_ARRAYS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT&
            & = INT(Z"0080000",C_INT) !< See #RSB_FLAG_EXPERIMENTAL_IN_PLACE_PERMUTATION_SORT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_SYMMETRIC&
            & = INT(Z"0400000",C_INT) !< See #RSB_FLAG_SYMMETRIC.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_HERMITIAN&
            & = INT(Z"0800000",C_INT) !< See #RSB_FLAG_HERMITIAN.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS&
            & = INT(Z"01000000",C_INT) !< See #RSB_FLAG_RECURSIVE_MORE_LEAVES_THAN_THREADS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG&
            & = INT(Z"08000000",C_INT) !< See #RSB_FLAG_RECURSIVE_SUBDIVIDE_MORE_ON_DIAG.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS&
            & = INT(Z"040000000",C_INT) !< See #RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_USE_CSR_RESERVED&
            & = INT(Z"0200000",C_INT) !< See #RSB_FLAG_USE_CSR_RESERVED.
! Composite flags 
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DEFAULT_STORAGE_FLAGS&
            & = (RSB_FLAG_WANT_BCSS_STORAGE+&
            &RSB_FLAG_WANT_COO_STORAGE) !< See #RSB_FLAG_DEFAULT_STORAGE_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS&
            & = RSB_FLAG_WANT_COO_STORAGE  !< See #RSB_FLAG_DEFAULT_COO_MATRIX_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS&
            & = RSB_FLAG_WANT_BCSS_STORAGE !< See #RSB_FLAG_DEFAULT_CSR_MATRIX_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS&
            & = (RSB_FLAG_QUAD_PARTITIONING+&
            &RSB_FLAG_USE_HALFWORD_INDICES+&
            &RSB_FLAG_WANT_COO_STORAGE+&
            &RSB_FLAG_WANT_BCSS_STORAGE) !< See #RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DEFAULT_MATRIX_FLAGS&
            & = RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS !< See #RSB_FLAG_DEFAULT_MATRIX_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_IDENTICAL_FLAGS&
            & = RSB_FLAG_NOFLAGS !< See #RSB_FLAG_IDENTICAL_FLAGS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_LOWER_HERMITIAN&
            & = (RSB_FLAG_HERMITIAN + RSB_FLAG_LOWER) !< See #RSB_FLAG_LOWER_HERMITIAN.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_UPPER_HERMITIAN&
            & = (RSB_FLAG_HERMITIAN + RSB_FLAG_UPPER) !< See #RSB_FLAG_UPPER_HERMITIAN.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_LOWER_TRIANGULAR&
            & = (RSB_FLAG_TRIANGULAR + RSB_FLAG_LOWER) !< See #RSB_FLAG_LOWER_TRIANGULAR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_UPPER_TRIANGULAR&
            & = (RSB_FLAG_TRIANGULAR + RSB_FLAG_UPPER) !< See #RSB_FLAG_UPPER_TRIANGULAR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_LOWER_SYMMETRIC&
            & = (RSB_FLAG_SYMMETRIC + RSB_FLAG_LOWER) !< See #RSB_FLAG_LOWER_SYMMETRIC.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_DIAGONAL&
            & = (RSB_FLAG_UPPER_TRIANGULAR + RSB_FLAG_LOWER_TRIANGULAR) !< See #RSB_FLAG_DIAGONAL.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_UPPER_SYMMETRIC&
            & = (RSB_FLAG_SYMMETRIC + RSB_FLAG_UPPER) !< See #RSB_FLAG_UPPER_SYMMETRIC.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_USE_HALFWORD_INDICES_CSR&
            & = (RSB_FLAG_USE_HALFWORD_INDICES+&
            &RSB_FLAG_USE_CSR_RESERVED) !< See #RSB_FLAG_USE_HALFWORD_INDICES_CSR.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_USE_HALFWORD_INDICES_COO&
            & = (RSB_FLAG_USE_HALFWORD_INDICES+&
            &RSB_FLAG_WANT_COO_STORAGE) !< See #RSB_FLAG_USE_HALFWORD_INDICES_COO.
      INTEGER(C_INT),PARAMETER&
            &::RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES&
            & = (RSB_FLAG_USE_HALFWORD_INDICES_COO+&
            &RSB_FLAG_USE_HALFWORD_INDICES_CSR) !< See #RSB_FLAG_MUTUALLY_EXCLUSIVE_SWITCHES.
! Transposition constants 
      INTEGER(C_INT),PARAMETER::RSB_TRANSPOSITION_N=INT(Z"04E",C_INT) 
      INTEGER(C_INT),PARAMETER::RSB_TRANSPOSITION_T=INT(Z"054",C_INT) 
      INTEGER(C_INT),PARAMETER::RSB_TRANSPOSITION_C=INT(Z"043",C_INT) 
! Numerical types constants 
      INTEGER(C_SIGNED_CHAR),PARAMETER&
            &::RSB_NUMERICAL_TYPE_SAME_TYPE=1 
      INTEGER(C_SIGNED_CHAR),PARAMETER&
            &::RSB_NUMERICAL_TYPE_INT=ICHAR('I') 
      INTEGER(C_SIGNED_CHAR),PARAMETER&
            &::RSB_NUMERICAL_TYPE_DOUBLE=ICHAR('D') 
      INTEGER(C_SIGNED_CHAR),PARAMETER&
            &::RSB_NUMERICAL_TYPE_FLOAT=ICHAR('S') 
      INTEGER(C_SIGNED_CHAR),PARAMETER&
            &::RSB_NUMERICAL_TYPE_FLOAT_COMPLEX=ICHAR('C') 
      INTEGER(C_SIGNED_CHAR),PARAMETER&
            &::RSB_NUMERICAL_TYPE_DOUBLE_COMPLEX=ICHAR('Z') 
! Other enumerations constants 
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_VERBOSE_INIT&
            &=INT(Z"0000001",C_INT)  !< See #RSB_IO_WANT_VERBOSE_INIT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_VERBOSE_EXIT&
            &=INT(Z"0000002",C_INT)  !< See #RSB_IO_WANT_VERBOSE_EXIT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_OUTPUT_STREAM&
            &=INT(Z"0000003",C_INT)  !< See #RSB_IO_WANT_OUTPUT_STREAM.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_SORT_METHOD&
            &=INT(Z"0000004",C_INT)  !< See #RSB_IO_WANT_SORT_METHOD.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_CACHE_BLOCKING_METHOD&
            &=INT(Z"0000005",C_INT)  !< See #RSB_IO_WANT_CACHE_BLOCKING_METHOD.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_SUBDIVISION_MULTIPLIER&
            &=INT(Z"0000006",C_INT)  !< See #RSB_IO_WANT_SUBDIVISION_MULTIPLIER.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_VERBOSE_ERRORS&
            &=INT(Z"0000007",C_INT)  !< See #RSB_IO_WANT_VERBOSE_ERRORS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_BOUNDED_BOX_COMPUTATION&
            &=INT(Z"0000008",C_INT)  !< See #RSB_IO_WANT_BOUNDED_BOX_COMPUTATION.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_EXECUTING_THREADS&
            &=INT(Z"0000009",C_INT)  !< See #RSB_IO_WANT_EXECUTING_THREADS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE&
            &=INT(Z"0000010",C_INT)  !< See #RSB_IO_WANT_EXTRA_VERBOSE_INTERFACE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING&
            &=INT(Z"0000011",C_INT)  !< See #RSB_IO_WANT_MEMORY_HIERARCHY_INFO_STRING.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_IS_INITIALIZED_MARKER&
            &=INT(Z"0000012",C_INT)  !< See #RSB_IO_WANT_IS_INITIALIZED_MARKER.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_MEM_ALLOC_CNT&
            &=INT(Z"0000013",C_INT)  !< See #RSB_IO_WANT_MEM_ALLOC_CNT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_MEM_ALLOC_TOT&
            &=INT(Z"0000014",C_INT)  !< See #RSB_IO_WANT_MEM_ALLOC_TOT.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_LEAF_LEVEL_MULTIVEC&
            &=INT(Z"0000015",C_INT)  !< See #RSB_IO_WANT_LEAF_LEVEL_MULTIVEC.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_MAX_MEMORY_ALLOCATIONS&
            &=INT(Z"0000016",C_INT)  !< See #RSB_IO_WANT_MAX_MEMORY_ALLOCATIONS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_MAX_MEMORY_ALLOCATED&
            &=INT(Z"0000017",C_INT)  !< See #RSB_IO_WANT_MAX_MEMORY_ALLOCATED.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_LIBRSB_ETIME&
            &=INT(Z"0000018",C_INT)  !< See #RSB_IO_WANT_LIBRSB_ETIME.
      INTEGER(C_INT),PARAMETER&
            &::RSB_IO_WANT_VERBOSE_TUNING&
            &=INT(Z"0000019",C_INT)  !< See #RSB_IO_WANT_VERBOSE_TUNING.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_NORM_ONE&
            &=INT(Z"000001001",C_INT)  !< See #RSB_EXTF_NORM_ONE.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_NORM_TWO&
            &=INT(Z"000001002",C_INT)  !< See #RSB_EXTF_NORM_TWO.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_NORM_INF&
            &=INT(Z"000001003",C_INT)  !< See #RSB_EXTF_NORM_INF.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_SUMS_ROW&
            &=INT(Z"000001004",C_INT)  !< See #RSB_EXTF_SUMS_ROW.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_SUMS_COL&
            &=INT(Z"000001005",C_INT)  !< See #RSB_EXTF_SUMS_COL.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_ASUMS_ROW&
            &=INT(Z"000001006",C_INT)  !< See #RSB_EXTF_ASUMS_ROW.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_ASUMS_COL&
            &=INT(Z"000001007",C_INT)  !< See #RSB_EXTF_ASUMS_COL.
      INTEGER(C_INT),PARAMETER&
            &::RSB_EXTF_DIAG&
            &=INT(Z"000000004",C_INT)  !< See #RSB_EXTF_DIAG.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MARF_RGB&
            &=INT(Z"000000001",C_INT)  !< See #RSB_MARF_RGB.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MARF_EPS_S&
            &=INT(Z"000000010",C_INT)  !< See #RSB_MARF_EPS_S.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MARF_EPS_B&
            &=INT(Z"000000020",C_INT)  !< See #RSB_MARF_EPS_B.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MARF_EPS&
            &=INT(Z"000000030",C_INT)  !< See #RSB_MARF_EPS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MARF_EPS_L&
            &=INT(Z"000000070",C_INT)  !< See #RSB_MARF_EPS_L.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T&
            &=INT(Z"000000001",C_INT)  !< See #RSB_MIF_INDEX_STORAGE_IN_BYTES__TO__SIZE_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T&
            &=INT(Z"000000002",C_INT)  !< See #RSB_MIF_INDEX_STORAGE_IN_BYTES_PER_NNZ__TO__RSB_REAL_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T&
            &=INT(Z"000000004",C_INT)  !< See #RSB_MIF_MATRIX_ROWS__TO__RSB_COO_INDEX_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T&
            &=INT(Z"000000008",C_INT)  !< See #RSB_MIF_MATRIX_COLS__TO__RSB_COO_INDEX_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T&
            &=INT(Z"000000010",C_INT)  !< See #RSB_MIF_MATRIX_NNZ__TO__RSB_NNZ_INDEX_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_TOTAL_SIZE__TO__SIZE_T&
            &=INT(Z"000000020",C_INT)  !< See #RSB_MIF_TOTAL_SIZE__TO__SIZE_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T&
            &=INT(Z"000000040",C_INT)  !< See #RSB_MIF_MATRIX_FLAGS__TO__RSB_FLAGS_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T&
            &=INT(Z"000000080",C_INT)  !< See #RSB_MIF_MATRIX_TYPECODE__TO__RSB_TYPE_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_MATRIX_INFO__TO__CHAR_P&
            &=INT(Z"000000100",C_INT)  !< See #RSB_MIF_MATRIX_INFO__TO__CHAR_P.
      INTEGER(C_INT),PARAMETER&
            &::RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T&
            &=INT(Z"000000200",C_INT)  !< See #RSB_MIF_LEAVES_COUNT__TO__RSB_BLK_INDEX_T.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_MUL&
            &=INT(Z"000000001",C_INT)  !< See #RSB_ELOPF_MUL.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_DIV&
            &=INT(Z"000000002",C_INT)  !< See #RSB_ELOPF_DIV.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_POW&
            &=INT(Z"000000004",C_INT)  !< See #RSB_ELOPF_POW.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_NEG&
            &=INT(Z"000000008",C_INT)  !< See #RSB_ELOPF_NEG.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_SCALE_ROWS&
            &=INT(Z"000000010",C_INT)  !< See #RSB_ELOPF_SCALE_ROWS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_SCALE_COLS&
            &=INT(Z"000000020",C_INT)  !< See #RSB_ELOPF_SCALE_COLS.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_SCALE_ROWS_REAL&
            &=INT(Z"000000040",C_INT)  !< See #RSB_ELOPF_SCALE_ROWS_REAL.
      INTEGER(C_INT),PARAMETER&
            &::RSB_ELOPF_SCALE_COLS_REAL&
            &=INT(Z"000000080",C_INT)  !< See #RSB_ELOPF_SCALE_COLS_REAL.
      INTEGER(C_INT),PARAMETER&
            &::RSB_PRECF_ILU&
            &0=INT(Z"000000001",C_INT)  !< See #RSB_PRECF_ILU0.
      TYPE(C_PTR),PARAMETER&
            &::RSB_NULL_INIT_OPTIONS&
            &=C_NULL_PTR  !< See #RSB_NULL_INIT_OPTIONS.
      TYPE(C_PTR),PARAMETER&
            &::RSB_NULL_EXIT_OPTIONS&
            &=C_NULL_PTR  !< See #RSB_NULL_EXIT_OPTIONS.
      END MODULE rsb
