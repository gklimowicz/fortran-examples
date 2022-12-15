dnl
dnl
dnl	RSB_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_NAME(TYPE,BLOCKORIENTED)
dnl	-------------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function dispatcher function name.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_NAME',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
rsb__do_mergesort`_'blockoriented`'dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_ARGS(TYPE,BLOCKORIENTED)
dnl	-------------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function dispatcher function arguments.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_ARGS',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
dnl
(
	rsb_coo_idx_t *iarray,
	rsb_coo_idx_t *jarray,
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *biarray,
',`')dnl
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *bjarray,'
,`')dnl
	void *array,
	rsb_nnz_idx_t length, 
ifelse(blockoriented,`BCSR',`	rsb_coo_idx_t mb,
',`')dnl
ifelse(blockoriented,`BCSR',`	rsb_coo_idx_t kb,
',`')dnl
	rsb_coo_idx_t *iresult,
	rsb_coo_idx_t *jresult,
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *biresult,
',`')dnl
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *bjresult,
',`')dnl
	void *result,
	rsb_type_t type)dnl
dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_PROTOTYPE(TYPE,BLOCKORIENTED)
dnl	------------------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function dispatcher function prototype.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_PROTOTYPE',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
rsb_err_t RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_NAME(mtype,blockoriented)dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_ARGS(mtype,blockoriented)dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER(TYPE,BLOCKORIENTED)
dnl	--------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function dispatcher.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER',`dnl
dnl
pushdef(`types',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_DISPATCHER_PROTOTYPE(,blockoriented)
dnl
{
	/*!
	 * \ingroup gr_util
	 *	This function will sort the non zero elements of a sparse ifelse(blockoriented,`VBR',`blocked ')`'matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 * 	\param length  the input  arrays length
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  mtype array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
ifelse(blockoriented,`CSR',`dnl
	 *	FIXME : UNDOCUMENTED
')dnl
ifelse(blockoriented,`VBR',`dnl
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
')dnl
ifelse(blockoriented,`BCSR',`
	 *	FIXME : UNDOCUMENTED
',`')dnl
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */


dnl
pushdef(`blockoriented',$2)dnl
pushdef(`args',`iarray, jarray,
ifelse(blockoriented,`VBR',`	biarray,bjarray,',`')dnl
ifelse(blockoriented,`BCSR',`	mb,kb,',`')dnl
array, length,
iresult, jresult,
ifelse(blockoriented,`VBR',`	biresult,bjresult,',`')dnl
result')dnl
dnl
foreach(`mtype',RSB_M4_TYPES,`dnl
dnl
	if(`type' == RSB_M4_NUMERICAL_TYPE_PREPROCESSOR_SYMBOL(mtype))
	return RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_NAME(mtype,blockoriented)(args);dnl

	else
dnl
')dnl
dnl
popdef(`blockoriented')dnl
popdef(`args')dnl
dnl
	return RSB_ERR_UNSUPPORTED_TYPE;
}
dnl
popdef(`types')dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_NAME(TYPE,BLOCKORIENTED)
dnl	--------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates merge function name.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_NAME',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
rsb_do_merge_`'RSB_M4_CHOPSPACES(mtype)`_'blockoriented`'dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_ARGS(TYPE,BLOCKORIENTED)
dnl	--------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates merge function arguments.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_ARGS',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
(
		const rsb_coo_idx_t* RSB_M4_RESTRICT ileft, const rsb_coo_idx_t* RSB_M4_RESTRICT iright,  rsb_coo_idx_t*RSB_M4_RESTRICT iresult,
		const rsb_coo_idx_t* RSB_M4_RESTRICT jleft, const rsb_coo_idx_t* RSB_M4_RESTRICT jright,  rsb_coo_idx_t*RSB_M4_RESTRICT jresult,
ifelse(blockoriented,`VBR',`const rsb_coo_idx_t * RSB_M4_RESTRICT bileft, const rsb_coo_idx_t * RSB_M4_RESTRICT biright, rsb_coo_idx_t * RSB_M4_RESTRICT biresult,',`')dnl
ifelse(blockoriented,`VBR',`const rsb_coo_idx_t * RSB_M4_RESTRICT bjleft, const rsb_coo_idx_t * RSB_M4_RESTRICT bjright, rsb_coo_idx_t * RSB_M4_RESTRICT bjresult,',`')dnl
ifelse(blockoriented,`BCSR',`const rsb_coo_idx_t mb, const rsb_coo_idx_t kb,',`')dnl
		const mtype* left, const mtype* RSB_M4_RESTRICT right,  mtype* RSB_M4_RESTRICT result,
dnl		rsb_coo_idx_t left_index, rsb_coo_idx_t right_index, rsb_coo_idx_t result_index, 
dnl		rsb_coo_idx_t left_mod, rsb_coo_idx_t right_mod, rsb_coo_idx_t result_mod )dnl
		rsb_nnz_idx_t left_length,
		rsb_nnz_idx_t right_length )dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_PROTOTYPE(TYPE,BLOCKORIENTED)
dnl	-------------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates merge function prototype.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_PROTOTYPE',`
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
void RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_NAME(mtype,blockoriented)dnl
RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_ARGS(mtype,blockoriented)
dnl
popdef(`blockoriented')dnl
popdef(`mtype')
dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION(TYPE,BLOCKORIENTED)
dnl	---------------------------------------------------------------
dnl	Expands to the mergesort on coordinates merge function function.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION',`
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_PROTOTYPE(mtype,blockoriented)
{
	/*!
	 * \ingroup gr_util
	 * The merge function for our blockoriented matrix coefficients sorting.
	 *
	 * NOTE :This function is the mergesort bottleneck.
	 */
dnl	rsb_nnz_idx_t left_length=left_mod;
dnl	rsb_nnz_idx_t right_length = right_mod;
	register int left_index=0, right_index=0, result_index=0;
	
	/*
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
	          |<-  length ----->|
dnl	<--------------- mod ------------------------->
	          ^- index
	 */

dnl #define LEFT_ADVANCE		left_index =( left_index+1)% left_mod; left_length-- ;
dnl #define RIGHT_ADVANCE		right_index=(right_index+1)%right_mod; right_length--;
dnl #define RESULT_ADVANCE		result_index =( result_index+1)% result_mod;

#define LEFT_ADVANCE		left_index =( left_index+1); left_length-- ;
#define RIGHT_ADVANCE		right_index=(right_index+1); right_length--;
#define RESULT_ADVANCE		result_index =( result_index+1);

ifelse(blockoriented,`VBR',`dnl
#define RESULT_APPEND(IEL,JEL,BIEL,BJEL,EL)	\
',`dnl
#define RESULT_APPEND(IEL,JEL,EL)	\
')dnl
	iresult[result_index]=(IEL);  \
	jresult[result_index]=(JEL);  \
	result[result_index]=( EL);  \
ifelse(blockoriented,`VBR',`dnl
	biresult[result_index]=(BIEL);  \
',`')dnl
ifelse(blockoriented,`VBR',`dnl
	bjresult[result_index]=(BJEL);  \
',`')dnl
	RESULT_ADVANCE;

#define LRESULT_APPEND	\
	iresult[result_index]=ileft[left_index];\
	jresult[result_index]=jleft[left_index];\
ifelse(blockoriented,`VBR',`dnl
	biresult[result_index]=bileft[left_index];  \
',`')dnl
ifelse(blockoriented,`VBR',`dnl
	bjresult[result_index]=bjleft[left_index];  \
',`')dnl
	result[ result_index]= left[left_index];\
	RESULT_ADVANCE; \
	LEFT_ADVANCE;

#define RRESULT_APPEND	\
	iresult[result_index]=iright[right_index];\
	jresult[result_index]=jright[right_index];\
ifelse(blockoriented,`VBR',`dnl
	biresult[result_index]=biright[right_index];  \
',`')dnl
ifelse(blockoriented,`VBR',`dnl
	bjresult[result_index]=bjright[right_index];  \
',`')dnl
	 result[result_index]= right[right_index];\
	RESULT_ADVANCE; \
	RIGHT_ADVANCE; 

	while( left_length > 0 && right_length > 0)
	if(
	ifelse(blockoriented,`VBR',`
		bileft[left_index] < biright[right_index] ||
		(	bileft[left_index] == biright[right_index]	&&
			bjleft[left_index] <= bjright[right_index]	)
		)
	')dnl
	ifelse(blockoriented,`CSR',`
		ileft[left_index] < iright[right_index] ||
		(	ileft[left_index] == iright[right_index]	&&
			jleft[left_index] <= jright[right_index]	)
		)
	')dnl
	ifelse(blockoriented,`BCSR',`
		ileft[left_index]/mb < iright[right_index]/mb ||
		(	ileft[left_index]/mb == iright[right_index]/mb	&&
			jleft[left_index]/kb <= jright[right_index]/kb	)
		)
	')dnl
	ifelse(blockoriented,`PACK',`
		ileft[left_index].spmv_uauz < iright[right_index].v )
	')dnl
	{
		LRESULT_APPEND
	}
	else
	{
		RRESULT_APPEND
	}

	while( left_length  > 0 )
	{
		LRESULT_APPEND
	}
	while( right_length  > 0 )
	{
		RRESULT_APPEND
	}
#undef LEFT_ADVANCE
#undef RIGHT_ADVANCE
#undef RESULT_ADVANCE
#undef RESULT_APPEND
#undef LRESULT_APPEND
#undef RRESULT_APPEND

}
dnl
popdef(`blockoriented')dnl
popdef(`type')dnl
dnl
')
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_NAME(TYPE,BLOCKORIENTED)
dnl	--------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function name.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_NAME',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
dnl rsb_do_mergesort_`'mtype`'dnl
rsb_do_mergesort_`'RSB_M4_CHOPSPACES(mtype)`_'blockoriented`'dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_ARGS(TYPE,BLOCKORIENTED)
dnl	--------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function arguments.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_ARGS',`dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
(
	rsb_coo_idx_t *RSB_M4_RESTRICT iarray,
	rsb_coo_idx_t *RSB_M4_RESTRICT jarray,
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *RSB_M4_RESTRICT biarray,
',`')dnl
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *RSB_M4_RESTRICT bjarray,
',`')dnl
ifelse(blockoriented,`BCSR',`	rsb_coo_idx_t mb, rsb_coo_idx_t kb,
',`')dnl
	mtype *array,
	rsb_nnz_idx_t length, 
	rsb_coo_idx_t *RSB_M4_RESTRICT iresult,
	rsb_coo_idx_t *RSB_M4_RESTRICT jresult,
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *RSB_M4_RESTRICT biresult,
',`')dnl
ifelse(blockoriented,`VBR',`	rsb_coo_idx_t *RSB_M4_RESTRICT bjresult,
',`')dnl
	mtype *RSB_M4_RESTRICT result)
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_PROTOTYPE(TYPE,BLOCKORIENTED)
dnl	-------------------------------------------------------------------
dnl	Expands to the mergesort on coordinates function prototype.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_PROTOTYPE',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
rsb_err_t RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_NAME(`mtype',`blockoriented')dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_ARGS(`mtype',`blockoriented')dnl
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
')dnl
dnl
dnl
dnl
dnl	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION(TYPE,BLOCKORIENTED)
dnl	---------------------------------------------------------
dnl	Expands to the block oriented mergesort on coordinates function
dnl	for VBR partitioning.
dnl	The blocking is specified as a function argument.
dnl
define(`RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION',`dnl
dnl
pushdef(`mtype',$1)dnl
pushdef(`blockoriented',$2)dnl
dnl
RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_PROTOTYPE(mtype,blockoriented)
dnl
{
	/*!
	 *	\ingroup gr_util
	 *	This function will sort the non zero elements of a sparse blocked
	 *      mtype matrix.
	 *      It will read row and column indices arrays, the values array,
	 *	and will sort them in separate output arrays.
	 *
	 *	NOTE : This function could be optimized.
	 *
	 * 	\param iarray  the input  row    indices array
	 * 	\param jarray  the input  column indices array
	 * 	\param array   the input  mtype array
	 * 	\param iresult the output row    indices array
	 * 	\param jresult the output column indices array
	 * 	\param result  the output values array
ifelse(blockoriented,`VBR',`dnl
	 * 	\param biarray  the input  block row    indices array
	 * 	\param bjarray  the input  block column indices array
	 * 	\param biresult the output block row    indices array
	 * 	\param bjresult the output block column indices array
')dnl
	 *	Will sort  thethree arrays (iarray, jarray, array) following the 
	 *	criteria :
	 *
	 * 	(ia1,ja1)<=(ia2,ja2) iff (ia1<ia2) or ( (ia1==ia2) and (ja1<ja2) )
	 * 	i.e.: C (row major) ordering
	 */

	rsb_nnz_idx_t middle;
	rsb_coo_idx_t so=sizeof(rsb_coo_idx_t);

	rsb_coo_idx_t * ileft  ;
	rsb_coo_idx_t * iright ;
ifelse(blockoriented,`VBR',`
	rsb_coo_idx_t * bileft  ;
	rsb_coo_idx_t * biright ;
',`')dnl
	rsb_coo_idx_t * jleft  ;
	rsb_coo_idx_t * jright ;
ifelse(blockoriented,`VBR',`
	rsb_coo_idx_t * bjleft  ;
	rsb_coo_idx_t * bjright ;
',`')dnl
ifelse(RSB_M4_WANT_OMP,`1',`dnl
	size_t tn=0;
dnl	size_t nt;
')dnl
	mtype * left  ;
	mtype * right ;
	
#define LIMIT 1
	if(length==LIMIT)
	{
		*iresult = *iarray;
		*jresult = *jarray;
ifelse(blockoriented,`VBR',`
		*biresult = *biarray;
		*bjresult = *bjarray;
',`')dnl
		*(mtype*)result = *(mtype*)array;
	}
	if(length<=LIMIT) return RSB_ERR_NO_ERROR;
#undef LIMIT
	middle = length/2;

ifelse(blockoriented,`VBR',`
	bileft  = biarray;
	bjleft  = bjarray;
	biright = biarray+middle;
	bjright = bjarray+middle;
',`')dnl
	left  = array;
	right  = array+middle;
	ileft  = iarray;
	jleft  = jarray;
	iright = iarray+middle;
	jright = jarray+middle;

ifelse(`0',`1',`dnl 20121016 
ifelse(RSB_M4_WANT_OMP,`1',`dnl
`#'dnl
dnl       pragma omp parallel num_threads(rsb_global_session_handle.rsb_g_threads)
       pragma omp parallel
	/*	FIXME : warning : experimental */
	{
	tn = omp_get_thread_num();
	nt = omp_get_num_threads();
	if(tn==0)
')dnl
',`dnl
/* 20121016 commented out omp usage because broke serial compilation  */
	{
')dnl
	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_NAME(mtype,blockoriented)
	( ileft, jleft,
ifelse(blockoriented,`VBR',`dnl
		bileft, bjleft,
',`')dnl
ifelse(blockoriented,`BCSR',` mb, kb, ',`')dnl
		left,   middle,
	        iresult  ,       jresult,
ifelse(blockoriented,`VBR',`dnl
		biresult, bjresult,
',`')dnl
		result         );

ifelse(RSB_M4_WANT_OMP,`1',`dnl
	if(tn==1)
')dnl
	RSB_M4_MERGESORT_ON_COORDINATES_FUNCTION_NAME(mtype,blockoriented)
	(iright, jright,
ifelse(blockoriented,`VBR',`dnl
		biright, bjright,
',`')dnl
ifelse(blockoriented,`BCSR',` mb, kb, ',`')dnl
		right, length-middle,  iresult+middle  ,jresult+middle,
ifelse(blockoriented,`VBR',`dnl
		biresult+middle, bjresult+middle,
',`')dnl
	((mtype*)result)+middle  );
ifelse(RSB_M4_WANT_OMP,`1',`dnl
	}
',`dnl
	}
')dnl

	RSB_MEMCPY(ileft ,iresult       ,so*middle);
	RSB_MEMCPY(jleft ,jresult       ,so*middle);
ifelse(blockoriented,`VBR',`
	RSB_MEMCPY(bileft ,biresult       ,so*middle);
	RSB_MEMCPY(bjleft ,bjresult       ,so*middle);
',`')dnl
	RSB_MEMCPY(  left, result       ,sizeof(mtype)*middle);
	RSB_MEMCPY(iright,iresult+middle,so*(length-middle));
	RSB_MEMCPY(jright,jresult+middle,so*(length-middle));
ifelse(blockoriented,`VBR',`
	RSB_MEMCPY(biright ,biresult+middle       ,so*(length-middle));
	RSB_MEMCPY(bjright ,bjresult+middle       ,so*(length-middle));
',`')dnl
	RSB_MEMCPY( right, ((mtype*)result)+middle ,sizeof(mtype)*(length-middle));

	RSB_M4_MERGESORT_ON_COORDINATES_MERGE_FUNCTION_NAME(mtype,blockoriented)dnl
		(
			ileft,iright,iresult,
			jleft,jright,jresult,
ifelse(blockoriented,`BCSR',`	mb,kb,')
ifelse(blockoriented,`VBR',`dnl
			bileft, biright, biresult,
',`')dnl
ifelse(blockoriented,`VBR',`dnl
			bjleft, bjright, bjresult,
',`')dnl
			left, right, result,
dnl			0,0,0,
dnl			middle,length-middle,length
			middle,length-middle
			);
	return RSB_ERR_NO_ERROR;	/* ! */
}
dnl
popdef(`blockoriented')dnl
popdef(`mtype')dnl
dnl
')
dnl
dnl
dnl
