/*
Copyright (C) 2008-2021 Michele Martone

This file is part of librsb.

librsb is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

librsb is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details.

You should have received a copy of the GNU Lesser General Public
License along with librsb; see the file COPYING.
If not, see <http://www.gnu.org/licenses/>.

*/
/* @cond INNERDOC */
/*
 * @author Michele Martone
 */
#ifndef RSB_ERR_H_INCLUDED
#define RSB_ERR_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include "rsb_common.h"

#define RSB_ERRM_E_MTXAP	"Supplied NULL  matrix structure pointer !"
#define RSB_ERRM_E_MTXAPP	"Supplied NULL  matrix structure pointer pointer !"
#define RSB_ERRM_E_VIJ	"Supplied NULL  VA, IA, JA arrays !"
#define RSB_ERRM_E_VIJP	"Supplied NULL  VAP, IAP, JAP arrays !"
#define RSB_ERRM_CNHEAF	"Cannot reuse arrays while passing them as const."
#define RSB_ERRM_EM	"!\n"
#define RSB_ERRM_NULL_VA "Supplied NULL VA array pointer!\n"
#define RSB_ERRM_BCE	"Internal bounds computing error!\n"
#define RSB_ERRM_ZSM	"WARNING : zero sized malloc !\n"
#define RSB_ERRM_BFEANS	"Bad flags: RSB_FLAG_EXTERNALLY_ALLOCATED_ARRAYS are not supported here!\n"
#define RSB_ERRM_BFSAH	"Bad flags: both RSB_FLAG_SYMMETRIC and RSB_FLAG_HERMITIAN?!\n"
#define RSB_ERRM_SEOWS	"some error occurred while shuffling\n"
#define RSB_ERRM_ALSMINT	"a leaf submatrix is non terminal!\n"
#define RSB_ERRM_ANLSMIT	"a non leaf submatrix is terminal!\n"
#define RSB_ERRM_NT	"non-terminal matrix as input!\n"
#define RSB_ERRM_FAOTAFS	"failed allocation of temporary array for swap"
#define RSB_ERRM_ES	""
#define RSB_ERRM_NON_SQUARE	"matrix non-square\n"
#define RSB_ERRM_NPSFF	"NULL pointer supplied for filename."
#define RSB_ERRM_CBAEM	"Cannot build and empty matrix.\n"
#define RSB_ERRM_CMOINIY	"Column Major order is not implemented yet.\n"
#define RSB_ERRM_MDNFARTS	"matrix does not fit as RCSR (too sparse)!\n"
#define RSB_ERRM_FYRYNS	"fatal : your rsb_coo_idx_t type is not supported.."
#define RSB_ERRM_WTC		"error : specified type code is not valid\n"
#define RSB_ERRM_COVMUINS	"control of virtual memory usage is not supported\n"
#define RSB_ERRM_FCOVMU		"failed control of virtual memory usage\n"
#define RSB_ERRM_PFTM		"Problems finalizing the matrix!\n"
#define RSB_ERRM_WUF		"WARNING: unfinished code!\n"
#define RSB_ERRM_PAL		"probable allocation problem\n"
#define RSB_ERRM_NAOL		"no array of leaves ?\n"
#define RSB_ERRM_NNTC		"no nonzeros to clone ?\n"
#define RSB_ERRM_NDIANN		"no diagonal implicit and no nonzeros ?\n"
#define RSB_ERRM_INZCL		"invalid nonzeroes count left! did you pass only duplicates and zeros ?\n"
#define RSB_ERRM_DICBU		"diagonal-implicit matrix can't be updated on the diagonal!\n"
#define RSB_ERRM_CP		"cleanup problem\n"
#define RSB_ERRM_BNCS		"blocking not correctly specified!\n"
#define RSB_ERRM_CIIAUF		"Clique insertion is an unfinished functionality!\n"
#define RSB_ERRM_FMMTDT		"failed matrix multiplication to dense test\n"
#define RSB_ERRM_FMATD		"failed matrix add to dense\n"
#define RSB_ERRM_FMM		"failed matrix multiplication\n"
#define RSB_ERRM_FMC		"failed matrix cloning\n"
#define RSB_ERRM_CMINBC		"cloned matrix is not built correctly\n"
#define RSB_ERRM_FCMS		"Failed computing matrix sum.\n"
#define RSB_ERRM_FMATDBC	"failed matrix add to dense basic checksum\n"
#define RSB_ERRM_FYCITINS	"fatal : your rsb_coo_idx_t type is not supported.."
#define RSB_ERRM_WOPSTASA	"WARNING : overflow possible. Switching to another sort algorithm.\n"
#define RSB_ERRM_SLIINS		"error : seems like input is not sorted\n"
#define RSB_ERRM_EWEFCFD	"error while estimating fillin (corrupt fillin data?)\n"
#define RSB_ERRM_EWEFCTD	"error while estimating fillin (corrupt timing data?)\n"
#define RSB_ERRM_ERROR		"error\n"
#define RSB_ERRM_ESIIB		"error : sorting input is bad\n"
#define RSB_ERRM_SLSIB		"error : seems like sorting is bugged\n"
#define RSB_ERRM_UL		"error : seems like you did not initialize librsb!\n"
#define RSB_ERRM_SIL		"error initializing library!\n"
#define RSB_ERRM_SILTC          "Error initializing the library! Nevertheless continuing just to print configuration info (something might be wrong though, and this can help diagnose the problem). Program will return an error code anyway.\n"
#define RSB_ERRM_EDNC		"error during nonzeros compacting\n"
#define RSB_ERRM_ZCFRWAEC	"zero compacting function returned with an error code\n"
#define RSB_ERRM_EQRPF		"error reading performance file.\n"
#define RSB_ERRM_ELMPF		"error loading memory performance file.\n"
#define RSB_ERRM_AE		"allocation error\n"
#define RSB_ERRM_IE		"internal error ?\n"
#define RSB_ERRM_NL		"\n"
#define RSB_ERRM_MBE		"matrix build error!\n"
#define RSB_ERRM_BM		"bad mtxAp:"
#define RSB_ERRM_BS		"bad submatrix: "
#define RSB_ERRM_TS		"a problem occurred in triangular solve!\n"
#define RSB_ERRM_BOH_TRI	"triangular solve requires either a lower, or an upper triangle!\n"
#define RSB_ERRM_MV		"a problem occurred in sparse matrix-vector product!\n"
#define RSB_ERRM_NOLP		"timer-based profiling has not been enabled at configure time!\n"
#define RSB_ERRM_IMNIP		"input matrix is not an in place one!\n"
#define RSB_ERRM_DNSAMIWAFCB		"Did not supply a matrix initiated with rsb_mtx_alloc_from_coo_begin!\n"
#define RSB_ERRM_IPEWIEM	"Internal problem encounted when initiating an empty matrix\n"
#define RSB_ERRM_NO_XDR	"No XDR configured in: binary matrix I/O disabled.\n"
#define RSB_ERRM_MASM	"problems assembling / converting matrix.\n"
#define RSB_ERRM_BADDIM	"bad matrix dimensions\n"

/* RSB_ERR_* messages */
#define RSB_ERRM_GENERIC_ERROR "An unspecified error occurred."
#define RSB_ERRM_UNSUPPORTED_OPERATION "The user requested an operation which is not supported (e.g.: was opted out at build time)."
#define RSB_ERRM_UNSUPPORTED_TYPE "The user requested to use a type which is not supported (e.g.: was opted out at build time)."
#define RSB_ERRM_UNSUPPORTED_FORMAT "The user requested to use a matrix storage format which is not supported (e.g.: was opted out at build time)."
#define RSB_ERRM_INTERNAL_ERROR "An error occurred which is not apparently caused by a user's fault (internal error)."
#define RSB_ERRM_BADARGS "The user supplied corrupt or inconsistent data as argument."
#define RSB_ERRM_ENOMEM	"There is not enough dynamical memory to perform the requested operation."
#define RSB_ERRM_UNIMPLEMENTED_YET "The requested operation was not implemented yet in this code revision."
#define RSB_ERRM_LIMITS "The requested operation could not be executed, or index overflow will happen."
#define RSB_ERRM_NO_USER_CONFIGURATION "A file containing user set configuration was not present."
#define RSB_ERRM_CORRUPT_INPUT_DATA "User-supplied data (e.g.: from file) was corrupt."
#define RSB_ERRM_FAILED_MEMHIER_DETECTION "Memory hierarchy info failed to be detected. You can bypass this by setting a meaningful RSB_USER_SET_MEM_HIERARCHY_INFO environment variable."
#define RSB_ERRM_COULD_NOT_HONOUR_EXTERNALLY_ALLOCATION_FLAGS "User gave flags for an inplace constructor in a copy-based routine."
#define RSB_ERRM_UNSUPPORTED_FEATURE "The requested feature is not available because it was opted out or not configured at build time."
#define RSB_ERRM_NO_STREAM_OUTPUT_CONFIGURED_OUT "Output to stream feature has been disabled at configure time."
#define RSB_ERRM_INVALID_NUMERICAL_DATA "User gave some input with invalid numerical data."
#define RSB_ERRM_MEMORY_LEAK "Probable memory leak (user did not deallocate librsb structures before calling rsb_lib_exit())."
#define RSB_ERRM_ELEMENT_NOT_FOUND "Element not found by rsb_mtx_get_vals() or rsb_mtx_set_vals()."

const char *rsb__get_errstr_ptr(rsb_err_t errval);
rsb_err_t rsb__do_perror(FILE *stream, rsb_err_t errval);
/*const rsb_char_t * rsb__do_strerror_r(rsb_err_t errval, rsb_char_t * buf, size_t buflen);*/
rsb_err_t rsb__do_strerror_r(rsb_err_t errval, rsb_char_t * buf, size_t buflen);

#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif /* RSB_ERR_H_INCLUDED */
/* @endcond */
