/*                                                                                                                            

Copyright (C) 2008-2019 Michele Martone

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
/* @cond INNERDOC  */
/**
 * @file
 * @author Michele Martone
 * @brief
 * Auxiliary functionalities file.
 */

#include <stdlib.h>
#include "rsb.h"
#include "rsb_internals.h"

/* void qsort(void *base, size_t nmemb, size_t size,
	int(*compar)(const void *, const void *)); */

RSB_INTERNALS_COMMON_HEAD_DECLS
#if 0
/* 20121001 old code, declaring static vars */
/**
 * \ingroup gr_internals
 * \brief An internal use only structure.
 * */
static struct nzinfo
{
	int i;	/* index */
	int n;	/* value */
	rsb_flags_t flags;
};
#define BLOCK_START 1

const struct nzinfo *g_sort_array;

static int rsb_compar_nzinfo(const void * ap, const void * bp)
{
	int a=*(int*)ap;
	int b=*(int*)bp;
	return ( g_sort_array[a].n - g_sort_array[b].n );
}

static void rsb_qsort_nzinfo(const struct nzinfo*nnzpr, const struct nzinfo*nnzpc, rsb_coo_idx_t *nnzpr_s, rsb_coo_idx_t *nnzpc_s, size_t rows, size_t columns)
{
	/* sort by smallest */
	g_sort_array=nnzpr;
	qsort( nnzpr_s , (size_t)rows   , sizeof(rsb_coo_idx_t), &rsb_compar_nzinfo );
	g_sort_array=nnzpc;
	qsort( nnzpc_s , (size_t)columns, sizeof(rsb_coo_idx_t), &rsb_compar_nzinfo );
	g_sort_array=NULL;
}

#define MEDBLOCKSIZE	8		/* FIXME */
#define MAXBLOCKSIZE	16

static int rsb_do_partition(struct nzinfo * nnzpx,rsb_blk_idx_t * p_x,rsb_blk_idx_t *X_b, rsb_blk_idx_t maxk, int blocksize)
{
	int k;
	rsb_coo_idx_t last_i=0;

	
	*X_b=0;
	if( ! ( nnzpx[0].flags & BLOCK_START ) || nnzpx[0].i!=0 )
	{
		p_x[*X_b]=0;/* if not here, this will be performed in the next loop */
		(*X_b)++;
	}

	for(k=0;k<maxk;++k)
	{
		if(
			( nnzpx[k].flags & BLOCK_START) || 
			( nnzpx[k].i-last_i >= blocksize) ||
			( k>0 && k+1<maxk && nnzpx[k-1].n>nnzpx[k].n && 	/* local minumum */
				nnzpx[k+1].n>nnzpx[k].n && nnzpx[k].i-last_i /* >= MEDBLOCKSIZE */ ) 
		)
		{
			last_i=nnzpx[k].i;
			p_x[*X_b]=last_i;
			(*X_b)++;
			RSB_DO_FLAG_ADD(nnzpx[k].flags,BLOCK_START);	/* we mark a new block */
		}
	}
//	p_x[*X_b]=0;
//	for(k=0;k<*X_b;++k) p_x[k+1]=p_x[k+1]+p_x[k];
	p_x[*X_b]=maxk;
	return 0;
}

/* 
 * This implementation partitions a a matrix specified through its nonzero coordinates into a vbr partitioning
 * */
int rsb__util_nnz2vbr(const rsb_coo_idx_t * IA, const rsb_coo_idx_t * JA, const rsb_nnz_idx_t nnz, const rsb_blk_idx_t rows, const rsb_blk_idx_t columns, rsb_blk_idx_t **p_rp, rsb_blk_idx_t **p_cp, rsb_blk_idx_t *M_b, rsb_blk_idx_t *K_b, int blockrowsize, int blockcolumnsize)
{
/**
	FIXME : UNTESTED

*/

/* Example driver program:
#include <stdlib.h>
#include <rsb_partition.h>


int main()
{
	const int rows=10,columns=10;
	const int IA[]={0,1,1,2,3,4,5,6,7,8,9};
	const int JA[]={0,2,1,2,3,4,5,6,7,8,9};
	const int NNZ=11;
	int *p_rp;
	int *p_cp;
	return nnz2vbr(IA,JA,NNZ,rows,columns,&p_r,&p_c);
}
*/

	struct nzinfo *nnzpc=NULL  ,*nnzpr=NULL  ;
	rsb_coo_idx_t *nnzpc_s,*nnzpr_s;
	int k;
	rsb_blk_idx_t *p_r,*p_c;
	if(!M_b||!K_b||!IA||!JA||nnz<1||rows<1||columns<1||!p_rp||!p_cp)return -1;

	nnzpr=calloc(rows   ,sizeof(struct nzinfo));
	nnzpc=calloc(columns,sizeof(struct nzinfo));
	nnzpr_s=calloc(rows   ,sizeof(rsb_coo_idx_t));
	nnzpc_s=calloc(columns,sizeof(rsb_coo_idx_t));
	p_r=calloc(1+rows,sizeof(rsb_coo_idx_t));
	p_c=calloc(1+columns,sizeof(rsb_coo_idx_t));

	/* FIXME */
	if(!nnzpr_s)goto err;
	if(!nnzpc_s)goto err;
	if(!nnzpr)goto err;
	if(!nnzpc)goto err;
//	if(!p_r)goto err;
//	if(!p_c)goto err;

	/* check */
	for(k=0;k<nnz;++k) if(IA[k]>=rows   )goto err;
	for(k=0;k<nnz;++k) if(JA[k]>=columns)goto err;
	for(k=0;k<nnz;++k) if(IA[k]< 0   )goto err;
	for(k=0;k<nnz;++k) if(JA[k]< 0   )goto err;

	/* we count the nonzeros on each column and on each row */
	for(k=0;k<nnz;++k) nnzpr[IA[k]].i=IA[k], nnzpc[JA[k]].i=JA[k];
	for(k=0;k<nnz;++k) nnzpr[IA[k]].n++, nnzpc[JA[k]].n++;
	for(k=0;k<rows   ;++k) nnzpr_s[k]=k;
	for(k=0;k<columns;++k) nnzpc_s[k]=k;

	/* we sort the index arrays on the basis of the count value */
	rsb_qsort_nzinfo(nnzpr, nnzpc, nnzpr_s, nnzpc_s, rows, columns);
//	for(k=0;k<columns;++k) RSB_STDERR("(%d -> %d)\n", nnzpc[nnzpc_s[k]].i, nnzpc[nnzpc_s[k]].n) ; RSB_STDERR("---\n");
//	for(k=0;k<rows   ;++k) RSB_STDERR("(%d -> %d)\n", nnzpr[nnzpr_s[k]].i, nnzpr[nnzpr_s[k]].n) ; RSB_STDERR("---\n");
//	for(k=0;k<columns;++k) RSB_STDERR("%d ", nnzpc[nnzpc_s[k]].n) ; RSB_STDERR("\n");
//	for(k=0;k<rows   ;++k) RSB_STDERR("%d ", nnzpr[nnzpr_s[k]].n) ; RSB_STDERR("\n");

	
//	RSB_STDERR("c=[0"); for(k=0;k<columns;++k) RSB_STDERR(",%d", nnzpc[k].n) ; RSB_STDERR("];\n");
//	RSB_STDERR("r=[0"); for(k=0;k<rows   ;++k) RSB_STDERR(",%d", nnzpr[k].n) ; RSB_STDERR("];\n");

	/* now we determine blocks dimensions, using least populated columns/rows as split points */

	/* first heuristic : partition tightly, then expand */
	for(k=0;k<rows   ;++k) p_r[k]=k;
	for(k=0;k<columns;++k) p_c[k]=k;
	
	/* we sorted the arrays with per-column and per-row nnz counts */

	
	/* we start with a minimal partitioning (bogus) */
	*M_b= (rows    / ( MAXBLOCKSIZE/4 ));
	*K_b= (columns / ( MAXBLOCKSIZE/4 ));
	/* we use the minimal column and row indices as first block starting indices */
	for(k=0;k<*M_b;++k) { RSB_DO_FLAG_ADD(nnzpr[nnzpr_s[k]].flags,BLOCK_START); }
	for(k=0;k<*K_b;++k) { RSB_DO_FLAG_ADD(nnzpc[nnzpc_s[k]].flags,BLOCK_START); }
//	for(k=0;k<*M_b;++k) RSB_STDERR("break at %d\n",nnzpr[nnzpr_s[k]].i);

	/* 
	 * we create (right now, with no furhter strategy, but we would like to use local minimality criteria)
	 * intermediate blocks of average size
	 * */

	if(blockrowsize   <1)blockrowsize   =MEDBLOCKSIZE;/* FIXME */
	if(blockcolumnsize<1)blockcolumnsize=MEDBLOCKSIZE;/* FIXME */

	*M_b=0;
	*K_b=0;
	rsb_do_partition(nnzpr,p_r,M_b,rows   ,blockrowsize);
	rsb_do_partition(nnzpc,p_c,K_b,columns,blockcolumnsize);
//	RSB_STDERR("M_b:%d\n",*M_b);
//	RSB_STDERR("K_b:%d\n",*K_b);


//	for(k=0;k<*K_b;++k) p_c[k]=p_c[k+1]-p_c[k];
//	for(k=0;k<*M_b;++k) p_r[k]=p_r[k+1]-p_r[k];

//	for(k=0;k<M_b   ;++k) RSB_STDERR(" %d ", p_r[k]);RSB_STDERR("\n");
//	for(k=0;k<K_b   ;++k) RSB_STDERR(" %d ", p_c[k]);RSB_STDERR("\n");

	if(nnzpr_s) free(nnzpr_s);
	if(nnzpc_s) free(nnzpc_s);
	if(nnzpr) free(nnzpr);
	if(nnzpc) free(nnzpc);
//	if(p_r) free(p_r);
//	if(p_c) free(p_c);

	/* FIXME : MISSING REALLOC FOR P_R, P_C ! */
	*p_rp=p_r,
	*p_cp=p_c;

	/* FIXME : FREE MEMORY ! */
	return 0;
err:
	if(nnzpr_s) free(nnzpr_s);
	if(nnzpc_s) free(nnzpc_s);
	if(nnzpr) free(nnzpr);
	if(nnzpc) free(nnzpc);
//	if(p_r) free(p_r);
//	if(p_c) free(p_c);
	/* FIXME : FREE MEMORY ! */
	return -1;
}
#endif

rsb_bool_t rsb__should_rejoin_small_leaf(rsb_nnz_idx_t nnz, rsb_nnz_idx_t mk, rsb_nnz_idx_t uk, rsb_nnz_idx_t lk, rsb_type_t typecode)
{
	/**
	 *
	 */
	return
		(
		(( nnz-lk>0 && nnz-lk < ((rsb_global_session_handle.min_leaf_matrix_bytes)/RSB_SIZEOF(typecode))) ||
		( lk-mk>0 && lk-mk < ((rsb_global_session_handle.min_leaf_matrix_bytes)/RSB_SIZEOF(typecode))) ||
		( mk-uk>0 && mk-uk < ((rsb_global_session_handle.min_leaf_matrix_bytes)/RSB_SIZEOF(typecode))) ||
		( uk>0 && uk < ((rsb_global_session_handle.min_leaf_matrix_bytes)/RSB_SIZEOF(typecode))))
		&&
		(( nnz-lk>0 && nnz-lk > ((rsb_global_session_handle.avg_leaf_matrix_bytes)/RSB_SIZEOF(typecode))) ||
		( lk-mk>0 && lk-mk > ((rsb_global_session_handle.avg_leaf_matrix_bytes)/RSB_SIZEOF(typecode))) ||
		( mk-uk>0 && mk-uk > ((rsb_global_session_handle.avg_leaf_matrix_bytes)/RSB_SIZEOF(typecode))) ||
		( uk>0 && uk > ((rsb_global_session_handle.avg_leaf_matrix_bytes)/RSB_SIZEOF(typecode))))
		);
}

/* @endcond */
