# Copyright (C) 2008-2021 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

#
## octave-based tester : generates a program with hard coded matrices, right hand sides, and result vectors.
# Assumes that octave/matlab computations are correct.

# Requires octave-4.4 or newer.
#
# NOTE : we have scarce control on the nonzero density!
# NOTE : this tester is largely superseded by ./rsbench -Q
# NOTE : apparent problems could arise in case of integer overflows due to extra-big matrices and ill matrix building parameters (not so unlikely).
# TODO :
# * a failure/success counter, and more verbose errors
# * someday should split in pieces (separate files) the octave based tester and lots of code.
# * s/printf/RSB_ERROR/g
# * should run some test on non square matrices

source("./sbtg.m")

function dump_transposition(a,br,bc)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	ts="int";
	printf("/* type is %s */\n",ts);
	printf("\n");
	s=sparse(a);
	nz=nnz(a);
	m=size(a,1);
	k=size(a,2);
	I=zeros(nz,1);
	JI=zeros(nz,1);
	SVA=zeros(nz,1);
	l=1;
	for i=1:m;
	for j=1:k;
		if(s(i,j)!=0)
			I(l)=i-1;
			JI(l)=j-1;
			SVA(l)=s(i,j);
			++l;
		endif
	endfor
	endfor
	its="rsb_coo_idx_t";
	printf("\n%s NIA[]=",its);
	printf("%s",dump_c_vec(I));
	printf("\n%s NJA[]=",its);
	printf("%s",dump_c_vec(JI));
	printf("\n%s NSVA[]=",ts);
	printf("%s",dump_c_vec(SVA*0));
	printf("\n%s SVA[]=",ts);
	printf("%s",dump_c_vec(SVA));
	printf("\n");
end 

function dump_scale(a,br,bc,transi)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	SV=zeros(m+extra,1);
	SV(1:m)=linspace(1,m,m);
	ts="int";
	printf("/* type is %s */\n",ts);
	printf("\n");
	s=sparse(a);
	nz=nnz(a);
	m=size(a,1);
	k=size(a,2);
	I=zeros(nz,1);
	JI=zeros(nz,1);
	SVA=zeros(nz,1);
	l=1;
	for i=1:m;
	for j=1:k;
		if(s(i,j)!=0)
			I(l)=i-1;
			JI(l)=j-1;
			if transi>1
				SVA(l)=s(i,j)*SV(j);
			else
				SVA(l)=s(i,j)*SV(i);
			endif
			++l;
		endif
	endfor
	endfor
	its="rsb_coo_idx_t";
	printf("\n%s NIA[]=",its);
	printf("%s",dump_c_vec(I));
	printf("\n%s NJA[]=",its);
	printf("%s",dump_c_vec(JI));
	printf("\n%s SV[]=",ts);
	printf("%s",dump_c_vec(SV));
# the following will not work with integers :)
#	printf("\n%s ISV[]=",ts);
#	dump_c_vec(1.0/SV)
	printf("\n%s NSVA[]=",ts);
	printf("%s",dump_c_vec(SVA*0));
	printf("\n%s SVA[]=",ts);
	printf("%s",dump_c_vec(SVA));
	printf("\n");
end 

function dump_negation(a,br,bc)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	dx=zeros(k+extra,1);
	d =zeros(k+extra,1);
	ts="int";
	printf("/* type is %s */\n",ts);
	printf("\n");
	s=sparse(a);
	nz=nnz(a);
	m=size(a,1);
	k=size(a,2);
	I=zeros(nz,1);
	JI=zeros(nz,1);
	V=zeros(nz,1);
	l=1;
	for i=1:m;
	for j=1:k;
		if(s(i,j)!=0)
			I(l)=i-1;
			JI(l)=j-1;
			V(l)=s(i,j);
			++l;
		endif
	endfor
	endfor
	V=-V;
	its="rsb_coo_idx_t";
	printf("\n%s NIA[]=",its);
	printf("%s",dump_c_vec(I));
	printf("\n%s NJA[]=",its);
	printf("%s",dump_c_vec(JI));
	printf("\n/*const*/ %s NVA[]=",ts);
	printf("%s",dump_c_vec(V));
	printf("\n%s NNVA[]=",ts);
	printf("%s",dump_c_vec(V*0));
	printf("\n");
end 

function dump_infty_norm(a,br,bc,transi)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	dx=zeros(k+extra,1);
	d =zeros(k+extra,1);
	ts="int";
	printf("/* type is %s */\n",ts);
	printf("\n");
	if transi==2
		printf("const %s in= %d;\n",ts,norm(a',Inf));
	else
		printf("const %s in= %d;\n",ts,norm(a,Inf));
	endif
	printf("\n");
	printf("%s inx=0;",ts);
end 

function dump_getdiag(a,br,bc)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	dx=zeros(k+extra,1);
	d =zeros(k+extra,1);
	ts="int";
	printf("/* type is %s */\n",ts);
	printf("\n");
	printf("/* the result diagonal vector */\n",ts);
	printf("const %s d[]=",ts);
	printf("%s",dump_c_vec(diag(a)))
	printf("\n");
	printf("/* the vector which will store the result */\n",ts);
	printf("%s dx[]=",ts);
	printf("%s",dump_c_vec(dx))
end 

function dump_getrow(a,br,bc)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	r1x=zeros(k+extra,1);
	r1 =zeros(k+extra,1);
	ts="int";
	printf("/* type is %s */\n",ts);
	printf("\n");
	printf("/* the result vector */\n",ts);
	for i=1:m
	printf("const %s r%d[]=",ts,i);
		r1 (1:k)=a(i,1:k);
	printf("%s",dump_c_vec(r1))
	printf("\n");
	endfor
	printf("/* the vector which will store the result */\n",ts);
	printf("%s r1x[]=",ts);
	printf("%s",dump_c_vec(r1x))
end 

function dump_spmv(alpha,a,br,bc,transi)
	dump_spmm(alpha,a,br,bc,1,transi);
end 

function dump_spmm(alpha,a,br,bc,nrhs,transi)
	extra=br+bc; # padding far more than necessary: will do no harm
	m=size(a,1);
	k=size(a,2);
	b=zeros(k+extra,nrhs);
	b(1:k,:)=ones(k,nrhs);
	for nrhsi=1:nrhs
		b(1:k,nrhsi)*=nrhsi;
	endfor
	x=zeros(m+extra,nrhs);
	r=zeros(m+extra,nrhs);
	if transi>1
		x(1:m,:)=alpha*a'*b(1:k,:);
	else
		x(1:m,:)=alpha*a *b(1:k,:);
	endif
	ts="int";
	if m<20 
		printf("/* \n");
		a
		printf("*/ \n");
	endif
	printf("/* type is %s */\n",ts);
	printf("/* matrix in coo */\n",ts);
	printf("\n");
	printf("/* the vector which will store the result */\n",ts);
	printf("%s X[]=",ts);
	printf("%s",dump_c_vec(r));
	printf("\n");
	printf("/* the result vector */\n",ts);
	printf("const %s R[]=",ts);
	printf("%s",dump_c_vec(x));
	printf("\n");
	printf("/* the right hand side vector */\n",ts);
	printf("const %s B[]=",ts);
	printf("%s",dump_c_vec(b));
end 

global aops="";
global oops="";
global ops="";
global op="";
global main=0;

if nargin == 0
	# default unrolls
	rua=cua=linspace(1,4,4);
elseif nargin == 1
	# same unrolls
	rua=eval(["[",cell2mat(argv()(1)),"]"]);
	cua=eval(["[",cell2mat(argv()(1)),"]"]);
endif

if nargin >=2
	# different unrolls
	rua=eval(["[",cell2mat(argv()(1)),"]"]);
	cua=eval(["[",cell2mat(argv()(2)),"]"]);
endif


if nargin >= 5
	op=argv()(5){1};
	if strcmp(op,"main")
		main=1;
		op="";
	endif
endif

if nargin >= 4
#	ops=eval(["[",cell2mat(argv()(3)),"]"]);
	ops=argv()(4){1};
	oops=ops;
	#if want_op("other")
#	ops=strcat(ops,",");
	ops=char(strsplit(ops,","));
endif

if nargin >= 3
	aops=argv()(3){1};
	aops=char(strsplit(aops,","));
endif
global extra_ops=("transposition,getdiag,getrow,");


function a=want_op(o)
	global oops;
	global ops;
	global op;
	global main;
	#a=strfind(ops,o)
	#a=0;
	#a=(strmatch(op,o))
	a=0;
	a=strfind(strcat(oops,","),o);
	if a
		a=a(1);
	else
		a=0;
	endif
#	oops,op,o,a
#	for i=1:size(ops,1)
#		if strfind(ops(i,:),op)
#			a+=1;
#		endif
#	endfor
end


printf("/**\n@file\n@brief integer type coverage testing program\n\n*/\n");
printf("#include \"rsb.h\"\n");
printf("#include \"rsb_internals.h\"\n");
printf("#include <string.h>\n");
printf("/*\n");
ops
rua
cua
op
main
printf("*/\n");

if main
	printf("%s",rsb_octave_doc_c_header);
	printf("%s\n","RSB_INTERNALS_COMMON_HEAD_DECLS");
	printf("int main%s()\n{\n",op);
	printf("rsb_err_t errval=RSB_ERR_NO_ERROR;int fi,Ri,Ti,Ci;int octave_failed_tests=0;\n");
	for oi=1:size(ops,1)
		printf("int main_%s(void);\n",ops(oi,:)); # avoid implicit declaration
		printf("if((octave_failed_tests+=main_%s()))\n\tgoto err;\n",ops(oi,:));
	endfor
	if main
		oops=strcat("extra_ops",",",oops);
	endif
#	printf("int main()\n{return 0;}\n");
#	printf("return 0;\n");
#	printf("err: return -1;\n");
#	printf(";}\n");
#	quit;
else
	printf("int main_%s()\n{\n",op);
	printf("rsb_err_t errval=RSB_ERR_NO_ERROR;int fi,Ri,Ti,Ci;int octave_failed_tests=0;\n");
	#if !strfind(ops,op)
#	if !strcmp(ops,strcat(op," "))
#		printf("return 0;}\n");
#	endif
#	ops, op
	if want_op(op)==0
		printf("return 0;\n");
#		printf("}\n"); quit; #
	endif
endif


# We should split the test program in pieces: compilation of this sequence 
# of subprograms takes more than the sum of individual pieces.
	if want_here_op("spsv")
		# FIXME : should build arrays with all supported format flags, for each operation.
printf("rsb_flags_t Cflagsa[]={  0 }; \n");
	else
printf("rsb_flags_t Cflagsa[]={ /*RSB_FLAG_WANT_COLUMN_MAJOR_ORDER ,*/ 0 }; \n");
	endif
printf("rsb_flags_t Rflagsa[]={ RSB_FLAG_QUAD_PARTITIONING, 0 }; \n");
#printf("rsb_flags_t Tflagsa[]=RSB_ROWS_TRANSPOSITIONS_ARRAY; \n");
printf("rsb_flags_t flagsa[]={ \n");

#printf("#ifdef RSB_FLAG_DEFAULT_STORAGE_FLAGS\n RSB_FLAG_DEFAULT_STORAGE_FLAGS, \n#endif /* RSB_FLAG_DEFAULT_STORAGE_FLAGS */\n");
#printf("RSB_FLAG_DEFAULT_STORAGE_FLAGS\n");
printf("RSB_FLAG_DEFAULT_RSB_MATRIX_FLAGS \n");
printf("#ifdef RSB_MATRIX_STORAGE_BCSC \n RSB_FLAG_WANT_BCSS_STORAGE, \n#endif /*RSB_MATRIX_STORAGE_BCSC*/\n");
printf("#ifdef RSB_MATRIX_STORAGE_LR   \n RSB_FLAG_WANT_LINKED_STORAGE, \n#endif /* RSB_MATRIX_STORAGE_LR */   \n");
printf("#ifdef RSB_MATRIX_STORAGE_VBR  \n  0, \n#endif /* RSB_MATRIX_STORAGE_VBR */   \n");
printf("};\n");
flagsa=[
	# FIXME: and what about halfword switches ?
	"RSB_FLAG_EXPERIMENTAL_SWITCH_TO_COO";
	"RSB_FLAG_WANT_BCSS_STORAGE";
	"RSB_FLAG_WANT_LINKED_STORAGE ";
	"0"; ];
printf("if(rsb_lib_init(RSB_NULL_INIT_OPTIONS))\n\tgoto err;\n");
#flagsa=[ "RSB_FLAG_WANT_BCSS_STORAGE"; "RSB_FLAG_WANT_BCSS_STORAGE|RSB_FLAG_WANT_COLUMN_MAJOR_ORDER "; ];

###############################################################################
printf("{\n");
#for fi=1:size(flagsa,1)
#printf("for(Ti=0;Ti<sizeof(Tflagsa)/sizeof(rsb_trans_t);++Ti)\n");
printf("for(Ri=0;Ri<sizeof(Rflagsa)/sizeof(rsb_flags_t);++Ri)\n");
printf("for(Ci=0;Ci<sizeof(Cflagsa)/sizeof(rsb_flags_t);++Ci)\n");
printf("for(fi=0;fi<sizeof(flagsa)/sizeof(rsb_flags_t);++fi)\n");
#printf("if(Tflagsa[Ti]!=RSB_INVALID_TRANS)\n");
	if want_here_op("spsv")
#printf("if(Tflagsa[Ti]!=RSB_INVALID_TRANS)\n");
	endif
printf("{\n");
#flags=flagsa(fi,:);
#printf("rsb_flags_t flags=%s; /*:) */\n",flags);
printf("rsb_flags_t flags=flagsa[fi] | Rflagsa[Ri] | Cflagsa[Ci];\n");
#printf("rsb_flags_t trans=Tflagsa[Ti];\n");
printf("rsb_flags_t trans=RSB_NUMERICAL_TYPE_INVALID_TYPE ;\n");
#printf("flags|=RSB_FLAG_SHOULD_DEBUG;\n"); # for heavy debugging only
#printf("printf(\"%%d\\n\",flags);\n"); # print flags


l=1;u=20;
l=1;u=8;
l=1;u=10;
#l=1;u=4;
#l=5;u=8;
#l=3;u=3;
#l=1;u=1;
#for transi=1:3 
	if want_here_op("spsv")
		printf("flags |= RSB_FLAG_LOWER_TRIANGULAR;\n");
	endif
for transi=1:2 # this is integer testing, dude ..
	if want_here_op("spsv") && transi==2 ;  continue ; endif
	if transi==1
		printf("trans=RSB_TRANSPOSITION_N;\n");trans="RSB_TRANSPOSITION_N";
	elseif transi==2
		printf("trans=RSB_TRANSPOSITION_T;\n");trans="RSB_TRANSPOSITION_T";
	elseif transi==3
		printf("trans=RSB_TRANSPOSITION_C;\n");
	endif
for n=l:u
# TODO : should adapt to available blockings 
for bri=1:length(rua)
for bci=1:length(cua)
	br=rua(bri);
	bc=cua(bci);
	incx=1;incy=1;
	bis="";
	if br * bc != 1 ;  bis=sprintf(" blocked %d x %d",br,bc) ; endif # blocking info string
	fprintf(stderr,"creating test matrix %d/%d%s\n",n,u,bis);
	#fprintf(stderr,"creating test matrix %d/%d blocked %d x %d\n",n,u,br,bc);
	printf("{\n");

	its="rsb_coo_idx_t";
	ts="int";
	mti=1;
	mdi=1;

	printf("struct rsb_mtx_t *mtxAp=NULL;\n");
	printf("rsb_err_t errval=RSB_ERR_NO_ERROR;\n");
	printf("rsb_flags_t typecode=RSB_NUMERICAL_TYPE_INT;\n");
	printf("int Mb=%d,Kb=%d; /* plain unblocked */\n",br,bc);
	printf("char buf[RSB_CONST_MATRIX_IMPLEMENTATION_CODE_STRING_MAX_LENGTH];/* Flawfinder: ignore */\n");
	
	a=gen_test_matrix(op,n);
#	if transi>=2 ; a=a'; endif
#	if transi==3 ; a=conj(a); endif
	ts="int";
	printf("%s\n",dump_c_coo(a,ts,"c"))
	printf("\nif((mtxAp=rsb_mtx_alloc_from_coo_const(VA,IA,JA,nnz,typecode,nr,nc,Mb,Kb,flags,&errval))==NULL)\ngoto err;");
	printf("\nprintf(\"testing matrix of type %%s\\t\",rsb__sprint_matrix_implementation_code2(mtxAp,buf,flags));\n");
	printf("\n%s","rsb__fprint_matrix_implementation_code(mtxAp,buf,flags,stdout);\n");

	if want_here_op("v")
	printf("/* begin spmv_uaua test */\n");
	dump_spmv(1,a,br,bc,transi)
#	printf("\nif((mtxAp=rsb_allocate_bcsr_sparse_matrix(VA,IA,JA,nnz,typecode,nr,nc,Mb,Kb))==NULL)\ngoto err;");
	printf("\n\tif((errval=rsb_spmv(mtxAp,B,X))                    !=RSB_ERR_NO_ERROR)\n\t\tgoto err;\n");
	printf("\n#error internal error: wrong code generator parameters!");
	printf("\nif(memcmp(X,R,sizeof(%s)*nr))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(X,R,nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"spmv_uaua test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"spmv_uaua test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("/* end spmv_uaua test */\n");
	endif

#
	# FIXME : this test misses the += part of this operation!
	extra=br+bc; # padding far more than necessary: will do no harm
	if want_here_op("vadd")
	printf("/* begin spmvadd test */\n");
	dump_spmv(1,a,br,bc,transi) # yes
#	printf("\nif((mtxAp=rsb_allocate_bcsr_sparse_matrix(VA,IA,JA,nnz,typecode,nr,nc,Mb,Kb))==NULL)\ngoto err;");
#	printf("\nif((errval=rsb_fill_with_zeros(X,mtxAp->typecode,nr+%d,1))                    !=RSB_ERR_NO_ERROR)\ngoto err;",extra);
	printf("\nif((errval=rsb_spmv_add(mtxAp,B,X))                    !=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(memcmp(X,R,sizeof(%s)*nr))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(X,R,nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"spmvadd test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"spmvadd test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("/* end spmvadd test */\n");
	endif
#
	# FIXME : this test misses the -= part of this operation!
	if want_here_op("vsub")
	printf("/* begin spmvsub est */\n");
	dump_spmv(-1,a,br,bc,transi) # yes
	printf("\nif((errval=rsb_spmv_sub(mtxAp,B,X))                    !=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(memcmp(X,R,sizeof(%s)*nr))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(X,R,nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"spmvsub test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"spmvsub test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("/* end spmvsub test */\n");
	endif
#
	if want_here_op("m")
	# not implemented yet
	#printf("return 0;\n");

	printf("/* begin spmm test */\n");
	dump_spmm(1,a,br,bc,2,transi)
	printf("\nif((errval=rsb_m(mtxAp,B,X,k+%d+%d,nr+%d+%d,2))                    !=RSB_ERR_NO_ERROR)\ngoto err;",br,bc,br,bc);
	printf("\nif(memcmp(X,R,sizeof(%s)*nr*2))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(X,R,2*nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"spmm test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"spmm test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("/* end spmm test */\n");

	endif
#
	if want_here_op("getrow")
	printf("/* FIXME: skipping getrow test because of obsolete/superseded rsb__get_row_dense() ! */\n");
	printf("if(0)\n");
	printf("{\n");
	printf("/* begin getrow test */\n");
	dump_getrow(a,br,bc)
	for i=1:n
		printf("\nif((errval=rsb__get_row_dense(mtxAp,r1x,%d-1))                    !=RSB_ERR_NO_ERROR)\ngoto err;",i);
		printf("\nif(memcmp(r1x,r%d,sizeof(%s)*nc))",i,ts);
		printf("\n{\n");
		printf("\nif((errval=rsb__debug_print_vectors_diff(r1x,r%d,nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;",i);
		printf("\nRSB_OCTAVE_ERROR(\"getrow %d test matrix %d/%d blocked %d x %d is not ok\\n\");\n",i,n,u,br,bc);
		printf("\n}else printf(\"getrow %d test matrix %d/%d blocked %d x %d is ok\\n\");\n",i,n,u,br,bc);
	endfor
	printf("/* end getrow test */\n");
	printf("}\n");
	endif
#
	if want_here_op("getdiag")
	printf("/* begin getdiag test */\n");
	printf("{\n");
	dump_getdiag(a,br,bc)
	printf("\nif((errval = rsb__dodo_getdiag /*rsb_do_getdiag <- crashes */(mtxAp,dx))                    !=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(memcmp(dx,d,sizeof(%s)*nr))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(dx,d,nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"diag test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"diag test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("/* end getdiag test */\n");
	printf("}\n");
	endif
#
	mvec=0;
	lang="c";
	if want_here_op("spsv")
	printf("/* begin spsv test */\n");
	printf("{\n");
#	v=[-4,-1,0,1,4]; # THIS IS THE FULL, RIGHTEOUS TEST
#	v=[-4,0,4];	# YES, BUT THIS IS JUST FASTER :)
	v=[-4,4];
#	v=[1];
#	v=[-1]; # only 1 is possible in this (integer) case
	av=v;
	bv=v;
	for ai=1:length(av)
	printf("{\n");
		alpha=av(ai);
#		beta=bv(bi);
		beta=0;
		printf("%s",dump_spsv(a,mdi,br,bc,alpha,mvec,transi,ts,incx,incy,lang));
#		printf("\nint alpha=%d,beta=%d;",alpha,beta);
		printf("\nint alpha=%d;",alpha);
		printf("\n\tif((errval=rsb_spsv(trans,&alpha,mtxAp,y,1,y,1))  !=RSB_ERR_NO_ERROR)\n\t\tgoto err;\n");
		printf(check_spsv(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy));
#		printf("\nif(memcmp(y,cy,sizeof(%s)*nr))",ts);
#		printf("\n{\n");
#		printf("\nif((errval=rsb__debug_print_vectors_diff(y,cy,nr,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
#		printf("\nRSB_OCTAVE_ERROR(\"spsv test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
#		printf("\n}else printf(\"spsv test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("} \n");
	endfor
	printf("/* end spsv test */\n");
	printf("} \n");
	endif
#
	if want_here_op("csmm") ||  want_here_op("spmv")
	printf("/* begin spmv test */\n");
	printf("{\n");
#	v=[-4,-1,0,1,4]; # THIS IS THE FULL, RIGHTEOUS TEST
	v=[-4,0,4];	# YES, BUT THIS IS JUST FASTER :)
#	v=[-4,4];
	av=v;
	bv=v;
	for ai=1:length(av)
	for bi=1:length(bv)
	printf("{\n");
		alpha=av(ai);
		beta =bv(bi);
		printf(dump_csmm(a,mti,mdi,br,bc,alpha,mvec,beta,transi,ts,incx,incy,lang));
		printf("\nint alpha=%d,beta=%d;",alpha,beta);
#		printf("\nif((errval=rsb_csmm(mtxAp,x,y,&alpha,&beta))  !=RSB_ERR_NO_ERROR\ngoto err;");
		printf("\n\tif((errval=rsb_spmv(trans,&alpha,mtxAp,x,1,&beta,y,1))  !=RSB_ERR_NO_ERROR)\n\t\tgoto err;\n");
		printf(check_csmm(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy));
	printf("} \n");
	endfor
	endfor
	printf("/* end spmv test */\n");
	printf("} \n");
	endif
#
	if want_here_op("infty_norm")
	printf("/* begin infty_norm test */\n");
	printf("/*");
	trans
	a
	printf("*/");
	dump_infty_norm(a,br,bc,transi)
	printf("\nif((errval = rsb__do_infinity_norm(mtxAp,&inx,RSB_BOOL_FALSE,trans))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(in!=inx)");
	printf("\n{\n");
	printf("\nprintf(\"infty_norm : should be %%d, not %%d !\\n\",in,inx);\n");
	printf("\nRSB_OCTAVE_ERROR(\"infty_norm test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"infty_norm test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("/* end infty_norm test */\n");
	endif
#
	if want_here_op("negation")
#
	printf("/* begin negation test */\n");
	printf("{");
	printf("/*");
	-a
	printf("*/");
	dump_negation(a,br,bc)
	printf("\nif((errval=rsb_coo_sort(   NVA , NIA, NJA, nnz, nr, nc, typecode, mtxAp->flags ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval=rsb_negation(mtxAp))                    !=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval = rsb_mtx_get_coo(mtxAp, NNVA , NIA, NJA, RSB_FLAG_C_INDICES_INTERFACE ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval=rsb_coo_sort(  NNVA , NIA, NJA, nnz, nr, nc, typecode, mtxAp->flags ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(memcmp(NNVA,NVA,sizeof(%s)*nnz))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(NVA,NNVA,nnz,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"negation test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"negation test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("}");
	printf("\nif((errval=rsb_negation(mtxAp))                    !=RSB_ERR_NO_ERROR)\ngoto err;/*negating back :) */");
	printf("/* end negation test */\n");
#
	endif
#
	if false
	if want_here_op("transposition")
	printf("if(0)/* FIXME: rsb_sym_transpose is still not mature enough */");
	printf("{");
	printf("/*");
	a
	printf("*/");
	printf("/* begin transposition test */\n");
	dump_transposition(a,br,bc)
	printf("\nif((errval=rsb_coo_sort(   SVA , NIA, NJA, nnz, nr, nc, typecode, mtxAp->flags ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\n errval=rsb_sym_transpose(mtxAp);");
	printf("\n if(errval!=RSB_ERR_UNIMPLEMENTED_YET ){");
	printf("\n if(errval!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval = rsb_mtx_get_coo(mtxAp, NSVA , NJA, NIA, RSB_FLAG_C_INDICES_INTERFACE ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval=rsb_coo_sort(   NSVA , NJA, NIA, nnz, nr, nc, typecode, mtxAp->flags ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(memcmp(NSVA,SVA,sizeof(%s)*nnz))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(SVA,NSVA,nnz,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"transpose test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"transpose test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
	printf("\n errval=rsb_sym_transpose(mtxAp);/* transposing back */");
	printf("\n if(errval!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\n }");
	printf("\n else");
	printf("\n { printf(\"skipping transposition test: unsupported\\n\");}");
	printf("/* end transpose test */\n");
	printf("}");
	endif
	endif
#
	if want_here_op("scale")
	printf("{");
	printf("/*");
	a
	printf("*/");
	printf("/* begin scale test */\n");
	printf("/* (since scaling does modify the matrix, it should be the last op here.. ) */\n");
	dump_scale(a,br,bc,transi)
	printf("\nif((errval=rsb_coo_sort(   SVA , NIA, NJA, nnz, nr, nc, typecode, mtxAp->flags ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval=rsb__do_scal(mtxAp,SV,trans))                 !=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval = rsb_mtx_get_coo(mtxAp, NSVA , NIA, NJA, RSB_FLAG_C_INDICES_INTERFACE ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif((errval=rsb_coo_sort(   NSVA , NIA, NJA, nnz, nr, nc, typecode, mtxAp->flags ))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nif(memcmp(NSVA,SVA,sizeof(%s)*nnz))",ts);
	printf("\n{\n");
	printf("\nif((errval=rsb__debug_print_vectors_diff(SVA,NSVA,nnz,typecode,1,1,0))!=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("\nRSB_OCTAVE_ERROR(\"scale test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc);
	printf("\n}else printf(\"scale test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc);
#	printf("\nif((errval=rsb_scal(mtxAp,ISV))                 !=RSB_ERR_NO_ERROR)\ngoto err;");
	printf("/* end scale test */\n");
	printf("}");
	endif
	printf("\n rsb_mtx_free(mtxAp);");
#	printf("goto end;\n");
#
	printf("\n}\n");
	# done
end
end
end
end
printf("\n}/*fi*/\n");
printf("}\n");
#end	# fi
#printf("\n}/*Ri*/\n");

printf("goto end; /* support compiler happyness worldwide */\n");
printf("end:;\n");
printf("if( errval!=RSB_ERR_NO_ERROR)\ngoto err;\n");
printf("if((errval=rsb_lib_exit(RSB_NULL_EXIT_OPTIONS))!=RSB_ERR_NO_ERROR)\ngoto err;\n");
if main
	printf("\nif(octave_failed_tests)RSB_INFO(\"ERROR: %%d failed tests (please file a bug report with this output).\\n\",octave_failed_tests);else RSB_INFO(\"all tests passed.\\n\");\n");
	printf("return octave_failed_tests?-1:0;\nerr:rsb_perror(NULL,errval);return -1;\n}\n");
else
	printf("return octave_failed_tests;\n  err:rsb_perror(NULL,errval);return -1;\n \n");
	printf("                              goto ferr/*only here as a token use of ferr*/;ferr:rsb_perror(NULL,errval);return -1;\n}\n");
endif

