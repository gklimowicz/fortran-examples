dnl
dnl
dnl	@author: Michele Martone
dnl	@brief
dnl	A cache estimator code  (EXPERIMENTAL, FIXME)
dnl
dnl	FIXME: INCOMPLETE
dnl
include(`rsb_misc.m4')dnl
include(`do_unroll.m4')dnl
dnl
dnl
dnl
/* @cond INNERDOC */
/*!
 * @file
 * @author Michele Martone
 * @brief L1 cache probing code (OBSOLETE)
 */
#include "rsb_sys.h"	/* rsb__aligned_malloc */

#include <stdio.h>
#include <stdlib.h>

typedef size_t rsb_word_t ;/* FIXME */

define(`RSB_M4_CACHE_SCAN',`dnl
pushdef(`FORCE_EVICTION',ifelse($1,1,1,0))dnl
pushdef(`args',`$1')dnl
pushdef(`want_what',$2)dnl
dnl
ifelse(want_what,`function_identifier',`k_wordops_'FORCE_EVICTION)`'dnl
dnl
ifelse(want_what,`function_declaration',`dnl
int $0(args,`function_identifier')dnl
$0(args,`function_args');dnl
')dnl
dnl
ifelse(want_what,`function_definition',`dnl
static inline size_t $0(args,`function_identifier')dnl
$0(args,`function_args')dnl
$0(args,`function_body')dnl
')dnl
dnl
ifelse(want_what,`function_args',`dnl
(rsb_word_t *p)dnl
')dnl
dnl
ifelse(want_what,`function_body',`dnl
dnl
{
dnl pushdef(`WORDS',1024)dnl
pushdef(`CACHE_SIZE',eval(8192))dnl
pushdef(`WORD_SIZE',8)dnl
pushdef(`WORDS',eval(CACHE_SIZE/WORD_SIZE))dnl
pushdef(`CACHE_LINE_SIZE',eval(64))dnl
pushdef(`CACHE_WAYS',4)dnl
pushdef(`CACHE_SETS',eval(CACHE_SIZE/(CACHE_LINE_SIZE*CACHE_WAYS)))dnl
pushdef(`WORDS_PER_CACHE_LINE',eval(CACHE_LINE_SIZE/WORD_SIZE))dnl	words per cache line
pushdef(`CACHE_SET_OFFSET_WORDS',eval(CACHE_WAYS*WORDS_PER_CACHE_LINE))dnl
pushdef(`CACHE_SET_OFFSET',eval(CACHE_SIZE/CACHE_SETS))dnl
	/*
	 * Code for a cache size of CACHE_SIZE bytes, WORD_SIZE bytes sized words,
	 * CACHE_WAYS-way associativity, with CACHE_LINE_SIZE bytes sized cache lines,
	 * each line thus fitting WORDS_PER_CACHE_LINE words,
	 * for a total of CACHE_SETS cache sets, distanced CACHE_SET_OFFSET_WORDS words (CACHE_SET_OFFSET bytes) each.
	 *
	 */
dnl	/* CSI th cache set */
dnl	/* CWI th cache line touched */
	/* FIXME : NESTED LOOPS IS BUGGY */
	RSB_M4_SIMPLE_UNROLL(`CWI',`0',`eval(CACHE_WAYS)',`dnl
	RSB_M4_SIMPLE_UNROLL(`CSI',`0',`eval(CACHE_SETS)',`dnl
	RSB_M4_SIMPLE_UNROLL(`CLI',`0',`eval(WORDS_PER_CACHE_LINE)',`
	p[CLI + CWI*CACHE_SET_OFFSET_WORDS + CSI*WORDS_PER_CACHE_LINE]*=-1;dnl
	/* CLI^th word of CSI^th cache set of CWI^th associativity way */
	')
	')
	/* each cache set has been loaded with a minimum of cache misses */
	')
	return eval(CACHE_SETS*WORDS_PER_CACHE_LINE*CACHE_WAYS);
popdef(`CACHE_SET_OFFSET')dnl
popdef(`CACHE_SET_OFFSET_WORDS')dnl
popdef(`WORDS_PER_CACHE_LINE')dnl	words per cache line
popdef(`CACHE_SETS')dnl
popdef(`CACHE_WAYS')dnl
popdef(`CACHE_LINE_SIZE')dnl
popdef(`WORDS')dnl
popdef(`WORD_SIZE')dnl
popdef(`CACHE_SIZE')dnl
	/* FIXME : should follow nothing : CSI CLI CWI */
}
')dnl
popdef(`args')dnl
popdef(`want_what')dnl
popdef(`FORCE_EVICTION')dnl
dnl
dnl
')dnl
dnl

RSB_M4_CACHE_SCAN(1,`function_definition')
RSB_M4_CACHE_SCAN(0,`function_definition')

int main()
{
	size_t i,j=0,it,times=100000;
	rsb_word_t * p = NULL;
	size_t N,KW,K=1024,W;
	double t,bt;
	size_t ops;

	for(i=1;i<7;++i)
	{
		N=(1<<i)*K;			/* bytes */
		KW=N/(K*sizeof(rsb_word_t));	/* kilowords */
		W=N/(   sizeof(rsb_word_t));	/* words */

		//p = rsb__aligned_malloc(N,N);	/* we want aligned bytes*/
		p = rsb__aligned_malloc(((1<<8) * K),N);	/* we want aligned bytes*/
		if(!p)goto err;

		ops=0;				/* op count reset */
		t = - rsb_time();			/* clock reset */
		for(it=0;it<times;++it)	
/*			for(j=0;j<KW;++j)	*//* we process one kiloword at a time */
				ops += RSB_M4_CACHE_SCAN(0,`function_identifier')( p+j*K );
		t += rsb_time();

		RSB_STDOUT("%10d times, %10d Kwords == %10d bytes : %10lg secs : %10lg ops per sec\n",times,KW,N,t,((double)(ops))/t);

		ops=0;				/* op count reset */
		bt = - rsb_time();			/* clock reset */
		for(it=0;it<times;++it)	
/*			for(j=0;j<KW;++j)	*//* we process one kiloword at a time */
				ops += RSB_M4_CACHE_SCAN(1,`function_identifier')( p+j*K );
		bt += rsb_time();

		RSB_STDOUT("%10d times, %10d Kwords == %10d bytes : %10lg secs : %10lg ops per sec\n",times,KW,N,bt,((double)(ops))/bt);
		RSB_STDOUT("ratio = %lg\n",bt/t);
		RSB_STDOUT("\n");
		if(p){rsb__free(p);p=NULL;}
	}
	
	if(p)free(p);
	return 0;
err:
	return -1;
}

/* @endcond */
