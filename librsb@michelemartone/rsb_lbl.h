/*                                                                                                                            

Copyright (C) 2008-2015 Michele Martone

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
/** @file
 *  @brief	Macros for linked block formats (and more)
 *  @author Michele Martone
 * */

#ifndef RSB_LBL_H_INCLUDED
#define RSB_LBL_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef RSB_FLAG_WANT_LINKED_STORAGE

/**
 * \brief A local matrix block info structure.
 * */
struct rsb_block_tail_t
{
	/* TODO : USE SHORT INT'S */
/*	void *next_block_on_row;
	void *next_block_on_col;*/
	/* seems like using integers is faster ! */
#if RSB_BLOCK_SMALL_INDICES
	/* 
	 * DANGER  :
	 * WARNING : MISSING CHECKS IN THE CODE FOR BIGGER BLOCKS !
	 * DANGER  :
	 */
	#define intrablock_int unsigned char
	#define interblock_int unsigned short int
	#define index_int      unsigned short int
#else
/*	#define intrablock_int short unsigned int*/
	#define intrablock_int unsigned int	/* FIXME : temporary (for debugging!) */
	#define interblock_int unsigned int
	#define index_int      unsigned int
#endif /* RSB_BLOCK_SMALL_INDICES */
	intrablock_int block_rows;
	intrablock_int block_columns;
	interblock_int block_row;
	interblock_int block_column;
	index_int base_column;
	index_int base_row;
#if 0
	void(*block_spmv_f_p[])(const void*,const void*,void*,int,int) /* block multiplication function pointer */ 
	rsb_flags_t flags;
#endif
	/* only for test purposes right now  */
	/* FIXME : and data alignment ? we should not break it ! */
#ifdef RSB_WANT_BLOCK_TRAILING_STRUCT_QUICK
	void(*spmv_fp)(const double*,const double*,double*,int,int); /* sample spmv protorype */
#if RSB_BLOCK_SMALL_INDICES
	short int foo;	/* dummy 16 bits bits for 64 bits alignment purposes : FIXME : temporary */
#else
	int foo;	/* dummy 32 bits bits for 64 bits alignment purposes : FIXME : temporary */
#endif
#endif /* RSB_WANT_BLOCK_TRAILING_STRUCT_QUICK */
};


#define RSB_BLOCK_EXTRA_BYTES (sizeof(struct rsb_block_tail_t))

#define RSB_BLOCK_TRAILING_STRUCT_TRANSPOSE(bt) \
	{\
	intrablock_int ib_tmp; \
	index_int ii_tmp; \
	ib_tmp=bt->block_columns;bt->block_columns=bt->block_rows;bt->block_rows=ib_tmp; \
	ib_tmp=bt->block_column ;bt->block_column =bt->block_row ;bt->block_row =ib_tmp; \
	ii_tmp=bt->base_column;  bt->base_column  =bt->base_row  ;bt->base_row  =ii_tmp; \
	/*printf("br:%d\n",(bt)->block_row);*/ \
	/*printf("bc:%d\n",(bt)->block_column);*/ \
	/*printf("%d %d\n",block_rows_,block_columns_);*/ \
	}

/* trailing struct, full */
#define RSB_BLOCK_TRAILING_STRUCT_SET_(bt,block_row_,block_column_,block_rows_,block_columns_,base_row_,base_column_) \
	{\
	{(bt)->block_row =(block_row_ );} \
	{(bt)->block_rows=(block_rows_);} \
	{(bt)->block_column=(block_column_);} \
	{(bt)->block_columns=(block_columns_);} \
	{(bt)->base_row =(base_row_ );} \
	{(bt)->base_column=(base_column_);} \
	/*printf("br:%d\n",(bt)->block_row);*/ \
	/*printf("bc:%d\n",(bt)->block_column);*/ \
	/*printf("%d %d\n",block_rows_,block_columns_);*/ \
	}

/*	{(bt)->flags=0x0;}*/

/* still unsupported */
#ifdef RSB_WANT_BLOCK_TRAILING_STRUCT_QUICK
#define RSB_BLOCK_TRAILING_STRUCT_SET(bt,block_row,block_column,block_rows,block_columns,base_row,base_column) \
	{RSB_BLOCK_TRAILING_STRUCT_SET_(((struct rsb_block_tail_t *)(bt)),(block_row),(block_column),(block_rows),(block_columns),(base_row),(base_column)) \
	{((struct rsb_block_tail_t *)(bt))->spmv_fp=(RSB_double_spmv((block_rows),(block_columns)));} }
#else
#define RSB_BLOCK_TRAILING_STRUCT_SET(bt,block_row,block_column,block_rows,block_columns,base_row,base_column) \
	{RSB_BLOCK_TRAILING_STRUCT_SET_(((struct rsb_block_tail_t *)(bt)),(block_row),(block_column),(block_rows),(block_columns),(base_row),(base_column))}
#endif /* RSB_WANT_BLOCK_TRAILING_STRUCT_QUICK */

#define RSB_BLOCK_TRAILING_STRUCT_NEXT(m,bt) \
	(struct rsb_block_tail_t*)(((char*)(bt))+((m)->el_size)*((rsb_coo_idx_t)(bt)->block_rows)*((rsb_coo_idx_t)(bt)->block_columns) + RSB_BLOCK_EXTRA_BYTES)




/* * This macro returns the block offset (in term of bytes) of the k^th block of matrix m */
#define RSB_BLOCK_OFFSET(m,k) (((m)->indptr[(k)]*(m)->el_size) + (((m)->flags&RSB_FLAG_WANT_LINKED_STORAGE)?(k+1)*(RSB_BLOCK_EXTRA_BYTES):0 ))

/* * This macro returns the count of needed bytes for the blocks array. */
#define RSB_TOTAL_BLOCK_BYTES(matrix,options) \
	((RSB_TOTAL_BLOCK_ELEMENT_COUNT(matrix))*(matrix)->el_size+((matrix)->block_count+1)*(((matrix)->flags&RSB_FLAG_WANT_LINKED_STORAGE)?sizeof(struct rsb_block_tail_t):0))




/* the extra structure is placed at the block first byte */
#define RSB_BLOCK_TRAILING_STRUCT_GET(m,k)  \
	((struct rsb_block_tail_t*)((char*)(RSB_BLOCK_ADDRESS(m,k))-sizeof(struct rsb_block_tail_t)))

/* new */
/*#define _VBR_BLOCK_TRAILING_STRUCT_GET(m,k)  \
	(((m)->flags & RSB_FLAG_WANT_LINKED_STORAGE)?(RSB_BLOCK_TRAILING_STRUCT_GET(m,k)):-1) 
*/






/*	*	*	*	*	*	*	*	*	*	*	*/

/*
 * The following three macros should be used for scanning a whole matrix.
 * They are intended to be used in a situation when performance is important, when
 * specialized kernels are called on each matrix sub block.
 * */
#if 0
#ifdef RSB_WANT_BLOCK_TRAILING_STRUCT
/* Implementation for the trailing structures enhanced variant. */
/* WARNING : UNTESTED */
#define	RSB_GET_NEXT_BLOCK_POINTER(BP,M,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	/*										\
	 * *input*									\
	 * M		should be a valid rsb_mtx_t structure pointer		\
	 * *output*									\
	 * ROWVAR	will be set to the base row    of this block			\
	 * COLVAR	will be set to the base column of this block			\
	 * BLOCKROWSVAR	will be set to the rows   count of this block			\
	 * BLOCKCOLSVAR	will be set to the column count of this block			\
	 * BP		 will be set to the current block pointer			\
	 * */										\
	_bt = RSB_BLOCK_TRAILING_STRUCT_NEXT((M),_bt);					\
	_lastk=_k;									\
	(BLOCKROWVAR)=_bt->block_row;							\
	(BLOCKCOLUMNVAR)=_bt->block_column;						\
	(BLOCKROWSVAR)=_bt->block_rows;							\
	(BLOCKCOLSVAR)=_bt->block_columns;						\
	(ROWVAR)=_bt->base_row;								\
	(COLVAR)=_bt->base_column;							\
	(BP)=((char*)_bt)+sizeof(struct rsb_block_tail_t);						\
	_k++;			/* for the future macro calls */			\
	;
/* WARNING : UNTESTED */
#define RSB_GET_FIRST_BLOCK_POINTER(BP,M,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)	\
	int /*_i=0,_j=0,*/_k=0,_lastk=0;						\
	const struct rsb_block_tail_t *_bt;							\
	_bt = (struct rsb_block_tail_t*)RSB_BLOCK_TRAILING_STRUCT_GET((M),0) ;					\
	(BLOCKROWVAR)=_bt->block_row;							\
	(BLOCKCOLUMNVAR)=_bt->block_column;						\
	(BLOCKROWSVAR)=_bt->block_rows;							\
	(BLOCKCOLSVAR)=_bt->block_columns;						\
	(ROWVAR)=_bt->base_row;								\
	(COLVAR)=_bt->base_column;							\
	(BP)=((char*)_bt)+sizeof(struct rsb_block_tail_t);						\
	/*int _lasti=0;*/											\
	/*int _lastj=0;*/											\
	/*RSB_GET_NEXT_BLOCK_POINTER(BP,M,ROWVAR,COLVAR,BLOCKROWSVAR,BLOCKCOLSVAR,BLOCKROWVAR,BLOCKCOLUMNVAR)*/

/* WARNING : UNTESTED */
#define RSB_GOT_LAST_BLOCK_POINTER(M)	( _k >= (M)->block_count )
#endif	/* ifdef RSB_WANT_BLOCK_TRAILING_STRUCT */
#endif	/* 0 */


#else	/* RSB_FLAG_WANT_LINKED_STORAGE */

/* * This macro returns the block offset (in term of bytes) of the k^th block of matrix m */
#define RSB_BLOCK_OFFSET(m,k) (((m)->indptr[(k)]*(m)->el_size))

/* * This macro returns the count of needed bytes for the blocks array. */
#define RSB_TOTAL_BLOCK_BYTES(matrix,options) \
	((RSB_TOTAL_BLOCK_ELEMENT_COUNT(matrix))*(matrix)->el_size)

#define RSB_TOTAL_BLOCK_ELEMENT_COUNT(matrix) ((matrix)->element_count)

/*	((RSB_TOTAL_BLOCK_ELEMENT_COUNT(matrix))*(matrix)->el_size+((matrix)->block_count)*(RSB_BLOCK_EXTRA_BYTES))*/


#endif	/* RSB_FLAG_WANT_LINKED_STORAGE */

/* * This macro returns the block address of the k^th block of matrix m */
#define RSB_BLOCK_ADDRESS(m,k) (((char*)((m)->VA))+RSB_BLOCK_OFFSET((m),(k)))



#ifdef __cplusplus
}
#endif  /* __cplusplus */

#endif	/* RSB_LBL_H_INCLUDED */
/* @endcond */
