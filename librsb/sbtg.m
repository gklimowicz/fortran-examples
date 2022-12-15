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
# (Octave based) Sparse BLAS Tester Generating code
#
global rsb_octave_doc_f_header="!> @cond INNERDOC\n!> @file\n!> @author Michele Martone \n!> @brief This file is part of the Octave based test suite for librsb\n";
global rsb_octave_doc_f_footer="!> @endcond\n";
global rsb_octave_doc_c_header="/*!\n * @file\n * @author Michele Martone \n * @brief This file is part of the Octave based test suite for librsb\n */\n";
global matrix_diagonal=['e','i']; # explicit or implicit, FIXME
#global matrix_diagonal=['e']; # explicit or implicit, FIXME
#
#global matrix_types_array=['g'];
#global matrix_types_array=['u'];
#global matrix_types_array=['l'];
#global matrix_types_array=['g','s'];
#global matrix_types_array=['g','l'];
#global matrix_types_array=['u','l'];
#global matrix_types_array=['u','g'];
#global matrix_types_array=['g','u','l','s','h'];
global matrix_types_array=['g','u','l','s','h'];
#
global blas_ctor_funcs=['r']; # 'r'=random ;'c','b','v'
#
global blas_trans_codes_array=['n','t','c'];
#global blas_trans_codes_array=['n'];
#global blas_trans_codes_array=['c'];
#
global findent="      ";
global ftab="  ";
#
#global alpha_array=[1];
#global alpha_array=[-1];
#global alpha_array=[1,-1];
global alpha_array=[3,1,-1,-3];
#global alpha_array=[4,1,-1,-4];
#
#global beta_array=[1,-1]; # FIXME: beta != 1 is unsupported by the sparse blas
global beta_array=[1,0]; # FIXME: beta != 1 is unsupported by the sparse blas
#
global incx_array=[1,2];
#global incx_array=[1];
#global incx_array=[2];
#
global incy_array=[1,2];
#global incy_array=[1];
#
#global beta_array=[4,1,-1,-4];
#
#global blas_type_codes_array=['s','d','c','z'];
#global blas_type_codes_array=['s'];
sbtgprepend = getenv("abs_builddir");
if length (sbtgprepend ) > 0
	sbtgprepend = [sbtgprepend,"/"];
end
source([sbtgprepend,"sbtg-types.m"]);
#
#global max_random_matrix_size_a=[1,2,3,4];
#global max_random_matrix_size_a=[1,2,3,4,5];
#
global max_random_matrix_size_a=[2];
#
global blas_op_codes_num=3;
#global blas_op_codes_num=2;
#global blas_op_codes_num=1;
#
#rand("seed",0);
#
global sbtg_random_fixed_seq=0; # if 1, the generated file will be the same on multiple runs
global sbtg_random_init_seed=[0];
global sbtg_want_only_full_matrices=0; # if 1, the generated matrices will be full

function res=iscomplextype(tc)
	res=(tc=='c' || tc=='z');
end

function res=sbtg_rand(n)
	global sbtg_random_fixed_seq;
	global sbtg_random_init_seed;

	if sbtg_random_fixed_seq
		rand("state",sbtg_random_init_seed);
		sbtg_random_fixed_seq=0;
	endif
	res=rand(n);
end 

function res=dump_c_vec(v,lang)
	if nargin < 2
		lang="c";
	endif
	if lang == "c"
		res="{ ";
		vv=vec(v);
		if(size(vv,1)>0)
		  if(iscomplex(v))
			for i=1:size(vv,1)-1; res=sprintf("%s%d+%d*I, ",res,real(vv(i)),imag(vv(i))); endfor
			i=size(vv,1);res=sprintf("%s%d+%d*I ",res,real(vv(i)),imag(vv(i)));
		  else
			for i=1:size(vv,1)-1; res=sprintf("%s%d, ",res,vv(i)); endfor
			i=size(vv,1);res=sprintf("%s%d ",res,vv(i)); 
		  endif
		else
			res=sprintf("%s-1 /* a dummy value */",res); 
		endif
		res=sprintf("%s};",res);
	else
		res="(/";
		vv=vec(v);
		vlen=size(vv,1)-1;
		if(iscomplex(v))
			for i=1:vlen; res=sprintf("%s(%d.e0,%d.e0), ",res,real(vv(i)),imag(vv(i))); endfor
		else
			for i=1:vlen; res=sprintf("%s%d, ",res,vv(i)); endfor
		endif
		if size(vv,1)>0
			if(iscomplex(v))
				res=sprintf("%s(%d,%d)",res,real(vv(size(vv,1))),imag(vv(size(vv,1))));
			else
				res=sprintf("%s%d",res,vv(size(vv,1)));
			endif
		endif
		res=sprintf("%s/)",res);
	endif
end

function res=dump_vec(v,type,id,comment,lang,mvec)
	if lang == "c"
		res=sprintf("%s %s[]=%s",type,id,dump_c_vec(v,lang));
	else
		if(nargin<6)
			mvec=0;
		endif
		if size(v,1)==0
			vv=[-999999];
		else
			vv=v;
		endif
		if size(v,1)>2
			bs="&\n          &";
		else
			bs="";
		endif
		if(mvec > 0)
			sc=sprintf("%d",size(vv,1));
			c_vec=["(reshape(",dump_c_vec(vv,lang),",(/",sc,",1/)))"];
			sd=",1";
		else
			c_vec=dump_c_vec(vv,lang);
			sd="";
		endif
		res=sprintf("%s :: %s(%d%s)=%s%s",type,id,size(vv,1),sd,bs,c_vec);
		if size(v,1)==0
			res=[res," ! fortran does not support empty arrays"];
		endif
	endif
	if(nargin>=4)
		res=sprintf("%s%s",res,comment);
	endif
end

function res=decl_var(lang,type,id,val)
	if lang == "c"
		res=sprintf("%s %s=%d;\n",type,id,val);
	else
		res=sprintf("%s :: %s=%d\n",type,id,val);
	endif
end

function res=dump_c_coo(a,ts,lang)
	#global ts;
	global findent;
	#
	s=sparse(a);
	nz=nnz(a);
	nr=size(a,1);
	nc=size(a,2);
	I=zeros(nz,1);
	JI=zeros(nz,1);
	V=zeros(nz,1);
	l=1;
	for i=1:nr;
	for j=1:nc;
		if(s(i,j)!=0)
			I(l)=i-1;
			JI(l)=j-1;
			V(l)=s(i,j);
			++l;
		endif
	endfor
	endfor
	res="";
	if lang == "c"
#		its="rsb_coo_idx_t";
		its="int";
		cts="rsb_coo_idx_t";
		nts="rsb_nnz_idx_t";
#		res=sprintf("%s	const rsb_nnz_idx_t nnz=%d;\n",res,nz);
		res=sprintf("%s	const int nnz=%d;\n",res,nz);
		res=sprintf("%s	const %s nr=%d;\n",res,cts,nr);
		res=sprintf("%s	const %s nc=%d;\n",res,cts,nc);
		res=sprintf("%s\tconst %s\n\tconst %s\n\tconst %s",res,dump_vec(I,cts,"IA","",lang),dump_vec(JI,cts,"JA","",lang),dump_vec(V,ts,"VA","",lang));
	else
	# FIXME: WRITE ME
		its="INTEGER(KIND=RSB_IDX_KIND)";
		res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","nnz",nz)];
		res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","nr",nr)];
		res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","nc",nc)];
		res=[sprintf("%s%s%s\n%s%s\n%s%s",res,findent,dump_vec(I+1,its,"IA","",lang),findent,dump_vec(JI+1,its,"JA","",lang),findent,dump_vec(V,ts,"VA","",lang))];
	endif
end

function a=want_here_op(o)
	global op;
	global main;
	global oops;
	global extra_ops;
	if main
#		a=want_op(o)
		a=0;
		a=strfind(strcat(extra_ops,","),o);
		if a
			a=a(1);
		else
			a=0;
		endif
	else
		a=(strcmp(op,o));
	endif
end

function d=mydiagf(m)
	d=m-(tril(m,-1)+triu(m,1));
end

function s=gen_test_matrix(op,n,tc,mti,mdi)
	global matrix_types_array;
	global matrix_diagonal;
	global sbtg_want_only_full_matrices;
	# generate a sparse random matrix sized n
	m=factorial(4);
#	s=mod(ceil(sbtg_rand(n)*m),4);
	s=zeros(n);
	s+=1*ceil(sbtg_rand(n)>.7);
	s+=2*ceil(sbtg_rand(n)>.8);
	s+=3*ceil(sbtg_rand(n)>.9);
	if sbtg_want_only_full_matrices ; s+=1; endif
	s(1,1)=1; # fix to avoid empty matrices
	if (nargin>=3 && iscomplextype(tc)) # complex type: adding a random imaginary part
		si=gen_test_matrix(op,n);
	#	si=si*i;
		si=(si+transpose(si))*i; # more chances to get nonzeros
		s=si+s;
	endif
#	if want_here_op("spsv") ||  want_here_op("ussv") 
	if nargin>=4 # matrix type index
		mt=matrix_types_array(mti);
	else
		mt='g'; # normal
	endif
	if (length(op)>=4 && (op(1:4)=="spsv" ||  op(1:4)=="ussv")) || (mt != 'g')
		if (nargin == 3 && iscomplextype(tc)) # complex type: real diagonal
			for k=1:n
				s(k,k)=real(s(k,k));
			endfor
		endif
	   if mt=='u'
		   s=triu(s);
	   else
		   s=tril(s);
           endif
	   s=(s-mydiagf(s))+speye(n);
	endif
	if nargin>=5
		md=matrix_diagonal(mdi);
	else
		md='e';
	endif
	if md=='i' # diagonal implicit (1)
		for k=1:n
			s(k,k)=0;
		endfor
	endif
end

function res=check_message(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy,tc,op,lang)
	global blas_trans_codes_array;
	global matrix_diagonal;
	global blas_type_codes_array;
	global matrix_types_array;
	mt=matrix_types_array(mti);
	md=matrix_diagonal(mdi);
	trans=blas_trans_codes_array(transi);
	bs="";
	if nargin >= 13
	if lang == "f"
		bs="&\n          &";
	endif
	endif
	res=sprintf("type=%s dims=%dx%d sym=%s diag=%s blocks=%s%dx%d %s alpha=%2d beta=%2d incx=%d incy=%d trans=%c",
		tc,size(a,1),size(a,2),mt,md,bs,br,bc,op,alpha,beta,incx,incy,trans);
end

function res=check_csmm(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy,tc,op,lang,mvec)
	global findent;
	global ftab;
	global matrix_diagonal;
	if nargin<14
		mvec = 0;
	end
	#
	md=matrix_diagonal(mdi);
	ts="int";
	res="";
	n=size(a,1);
	if nargin <= 6+2+2
		typecode="typecode";
	else
		#typecode="RSB_BLAS_C_HANDLER_TO_C_TYPECODE(A)";
		typecode=["'",toupper(tc),"'"];
	endif
	if nargin < 7+2+2
		tc='d';
	endif
	if nargin < 8+2+2
		op="spmv";
	endif
	if nargin < 11+2
		lang="c";
	endif
	if(mvec > 0)
		sd="i,1";
	else
		sd="i";
	end
	if lang == "c"
	cm=[check_message(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy,tc,op,lang)];
	nok=[cm," is not ok"];
	ok=[cm," is ok"];
	res=[res,sprintf("	if( (errval = rsb__do_are_same(y,cy,nr,%s,%d,%d)) != RSB_ERR_NO_ERROR )",typecode,incy,incy),"{ "];
	res=[res," rsb__debug_print_vectors_diff(y,cy,nr,",typecode,sprintf(",%d,%d,RSB_VECTORS_DIFF_DISPLAY_N",incy,incy),"); "];
#	res=[res,"		RSB_ERROR(\"",nok,"\\n\");\n"];
#	res=[res,"		goto err;\n	}\n		else printf(\"",ok,"\\n\");\n"];
	res=[res," goto ferr; }\n		else printf(\"",ok,"\\n\");\n"];
#	res=[res,sprintf("\nif(memcmp(y,cy,sizeof(%s)*nr))",ts),"\n{\n"];
#	res=[res,"\t	if(( errval = rsb__debug_print_vectors_diff(y,cy,nr,",typecode,")) != RSB_ERR_NO_ERROR)\n		goto err;"];
#	res=[res,sprintf("\nRSB_OCTAVE_ERROR(\"spmv test matrix %d/%d blocked %d x %d is not ok\\n\");\n",n,u,br,bc)];
#	res=[res,printf("\n}else printf(\"spmv test matrix %d/%d blocked %d x %d is ok\\n\");\n",n,u,br,bc)];
#	res=[res,sprintf("\nRSB_OCTAVE_ERROR(\"spmv test matrix %d blocked %d x %d is not ok\\n\");\n",n,br,bc)];
#	res=[res,sprintf("		RSB_ERROR(\"%s test matrix %d blocked %d x %d is not ok\\n\");\n",op,n,br,bc)];
#	res=[res,sprintf("		goto err;\n	}\n	else\n		printf(\"%s test matrix %d blocked %d x %d is ok\\n\");\n",op,n,br,bc) ];
#	res=[res,sprintf("\nprintf(\"spmv test\\n\");\n")];
	elseif lang == "f"
		nok=["LCS,\" is not ok\""];
		ok=["LCS,\" is ok\""];
		res=[res,findent,sprintf("DO i=1,%d\n",size(a,1))]; # FIXME: with transA and non square matrices, this breaks
		res=[res,findent,ftab,"IF(y(",sd,").NE.cy(",sd,"))PRINT*,",nok,"\n"];
		res=[res,findent,ftab,"IF(y(",sd,").NE.cy(",sd,"))errval=-1\n"];
		res=[res,findent,ftab,"IF(y(",sd,").NE.cy(",sd,"))GOTO 9997\n"];
#		res=[res,findent,findent,"IF(y(i",sd,").NE.cy(i))THEN\n"];
#		res=[res,findent,findent,findent,"errval =-1\n"];
#		res=[res,findent,findent,findent,"GOTO 9997\n"];
#		res=[res,findent,findent,"ENDIF\n"];
		res=[res,findent,"ENDDO\n"];
		res=[res,findent,"PRINT*,",ok,"\n"];
#		res=[res,"rsb__debug_print_vectors_diff(y,cy,nr,",typecode,");\n"];
#		res=[res,"RSB_ERROR(\"",check_message(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy,tc,op)," is not ok\\n\");\n"];
#		res=[res,sprintf("		goto err;\n	}\n		else printf(\"%s is ok\\n\");\n",check_message(a,mti,mdi,br,bc,alpha,beta,transi,incx,incy,tc,op) ) ];
	elseif lang == "p"
		res=[res,findent,sprintf("DO i=1,%d\n",size(a,1))]; # FIXME: with transA and non square matrices, this breaks
#		res=[res,findent,findent,"IF(y(",sd,").NE.cy(",sd,"))PRINT*,\"",nok,"\"\n"];
		res=[res,findent,findent,"IF(y(",sd,").NE.cy(",sd,"))PRINT*,\"results mismatch:\",y,\"instead of\",cy","\n"];
		res=[res,findent,findent,"IF(y(",sd,").NE.cy(",sd,"))info=-1\n"];
		res=[res,findent,findent,"IF(y(",sd,").NE.cy(",sd,"))GOTO 9996\n"];
		res=[res,findent,"ENDDO\n"];
#		res=[res,findent,"PRINT*,\"",ok,"\"\n"];
	endif
end 

function res=check_spsv(a,mti,mdi,br,bc,alpha,beta_,transi,incx,incy,tc,op,lang,mvec)
	if nargin==6+2+2
		res=check_csmm(a,mti,mdi,br,bc,alpha,beta_,transi,incx,incy);
	elseif nargin==7+2+2
		res=check_csmm(a,mti,mdi,br,bc,alpha,beta_,transi,incx,incy,tc);
	elseif nargin==8+2+2
		res=check_csmm(a,mti,mdi,br,bc,alpha,beta_,transi,incx,incy,tc,op);
	elseif nargin==8+2+2+1
		res=check_csmm(a,mti,mdi,br,bc,alpha,beta_,transi,incx,incy,tc,op,lang);
	elseif nargin==8+2+2+1+1
		res=check_csmm(a,mti,mdi,br,bc,alpha,beta_,transi,incx,incy,tc,op,lang,mvec);
	endif
end 

function sx=stride_apply(x,incx)
	# we use zeros as sentinel values, but we could use some big values or nan's as well (e.g.: 99999)
	sx=zeros(size(x,1)*incx,1);
	for i=1:size(x,1)
		sx((i-1)*incx+1)=x(i);
	end
end

function res=dump_csmm(a,mti,mdi,br,bc,alpha,mvec,beta,transi,ts,incx,incy,lang)
	global findent;
	global matrix_types_array;
	global matrix_diagonal;
	#
	extra=br+bc; # padding far more than necessary: will do no harm
	extra=0; # FIXME: for now it's ok
	#
	mt=matrix_types_array(mti);
	md=matrix_diagonal(mdi);
	#
	#incx=1;incy=1;
	nr=size(a,1);
	nc=size(a,2);
	x=zeros(nc+extra,1);
	x(1:nc)=ones(nc,1);
	y=zeros(nr+extra,1);
	y(1:nr)=ones(nc,1)*3;
	cy=zeros(nr+extra,1);
	# we obtain a matrix from its packed version
	if mt == 's'
		aa=a+transpose(a);
		for k=1:nr;aa(k,k)-=a(k,k); end
	elseif mt == 'h'
		aa=a+a';
		for k=1:nr;aa(k,k)-=a(k,k); end
	else
		aa=a;
	endif
	if md=='i'
		aa=aa+speye(nr);
	endif
	# we compute the vectors
	if transi==3
		cy(1:nr)=beta*y(1:nr)+alpha*aa'*x(1:nc);
		transc='H';
	elseif transi==2
		cy(1:nr)=beta*y(1:nr)+alpha*transpose(aa)*x(1:nc);
		transc='T';
	elseif transi==1
		cy(1:nr)=beta*y(1:nr)+alpha*aa *x(1:nc);
		transc='1';
	endif
	res="";
	scy=stride_apply(cy,incy);
	sy=stride_apply(y,incy);
	sx=stride_apply(x,incx);
	if lang == "c"
		res=[res,sprintf("		/* x: %d */\n",length(sx))];
		res=[res,"\tconst ",dump_vec(sx,ts,"x","/* reference x */\n",lang)];
		res=[res,"\tconst ",dump_vec(scy,ts,"cy","/* reference cy after */\n",lang)];
		sbcy=sy;
#		res=[res,"\tconst ",dump_vec(sbcy,ts,"bcy","/* reference bcy before */\n",lang)];
#		res=[res,"\t",dump_vec(sy,ts,"y","/* y */\n",lang),"\n"];
#		res=[res,"\t",sprintf("rsb__memcpy(y,bcy,(%d*%d+%d)*sizeof(%s)); /* because y will get overwritten otherwise */\n",incy,nr,extra,ts)];
		res=[res,"\t",dump_vec(sbcy,ts,"y","/* y */\n",lang),"\n"];
		res=sprintf("%s\t\n	const char*lsc=\"System and hardcoded solution: y' <- y + %d A^%c * x \\n\"%s",res,alpha,transc,print_matrix(aa,"k"));
		res=sprintf("%s\t%s",res,print_matrix(cy,"k","y'"));
		res=sprintf("%s\t%s",res,print_matrix(y,"k","y"));
		res=sprintf("%s\t%s;",res,print_matrix(x,"k","x"));
	else
		res=[res,findent,dump_vec(sx,ts,"x","! reference x \n",lang,mvec)];
		res=[res,findent,dump_vec(scy,ts,"cy","! reference cy after \n",lang,mvec)];
		sbcy=sy;
#		res=[res,findent,dump_vec(sbcy,ts,"bcy","! reference bcy before \n",lang,mvec)];
#		res=[res,findent,dump_vec(sy,ts,"y","! y will be overwritten\n",lang,mvec),"\n"];
#		res=[res,findent,"y=bcy\n"];
		res=[res,findent,dump_vec(sbcy,ts,"y","! y will be overwritten\n",lang,mvec),"\n"];
	endif
end 

function res=dump_spsv(a,mdi,br,bc,alpha,mvec,transi,ts,incx,incy,lang)
	global findent;
	global matrix_diagonal;
	#
	md=matrix_diagonal(mdi);
	#
	# FIXME: problems with incy != 1
	beta=0;
	extra=br+bc; # padding far more than necessary: will do no harm
	extra=0; # FIXME
	nr=size(a,1);
	nc=size(a,2);
	x=zeros(nc+extra,1);
	x(1:nc)=ones(nc,1);
	y=zeros(nr+extra,1);
	y(1:nr)=ones(nc,1)*0;
	cy=zeros(nr+extra,1);
#	cy(1:nr)=beta*y(1:nr)+alpha*a*x(1:nc);
	if md=='i'
		aa=a+speye(nr);
	else
		aa=a;
	endif
	if(transi==1)
		cy(1:nr)=beta*y(1:nr)+aa *x(1:nc);
		transc='1';
	elseif(transi==2)
		cy(1:nr)=beta*y(1:nr)+transpose(aa)*x(1:nc);
		transc='T';
	elseif(transi==3)
		cy(1:nr)=beta*y(1:nr)+aa'*x(1:nc);
		transc='H';
	endif
	tmp=cy;cy=x*alpha*alpha;x=tmp*alpha;
#	if nr<20 
#		printf("/* \n");
#		a
#		printf("*/ \n");
#	endif
	scy=stride_apply(cy,incy);
	sy=stride_apply(y,incy);
	sx=stride_apply(x,incx);
	res="";
	if lang == "c"
		res=[res,sprintf("/* type is %s */\n",ts),"\n"];
	#	res=[res,"\tconst ",dump_vec(sx,ts,"x","/* reference x */\n",lang)]; # for ot-spsv.c
		res=[res,"\tconst ",dump_vec(scy,ts,"cy","/* reference cy after */\n",lang)];
	#	sbcy=sy;
	#	res=[res,"\tconst ",dump_vec(sbcy,ts,"bcy","/* reference y before */\n",lang)];
	#	res=[res,"\t",dump_vec(sy,ts,"y","/* y */\n",lang),"\n"];
		res=[res,"\t",dump_vec(sx,ts,"y","/* y */\n",lang),"\n"];
		res=sprintf("%s\t\n	const char*lsc=\"System and hardcoded solution: y' <- %d A^-%c * y \\n\"%s",res,alpha,transc,print_matrix(aa,"k"));
		res=sprintf("%s\t%s",res,print_matrix(cy,"k","y"));
		res=sprintf("%s\t%s;\n",res,print_matrix(x,"k","y'"));
	#	res=[res,"\t",sprintf("rsb__memcpy(y,bcy,%d*sizeof(%s)); /* because y will get overwritten otherwise */\n",nr+extra,ts)];
	#	res=[res,"\t",sprintf("rsb__memcpy(y,x,%d*sizeof(%s)); /* because y will get overwritten otherwise */\n",nr+extra,ts)];
	#	res=[res,"\t",sprintf("rsb__memcpy(y,x,(%d*%d+%d)*sizeof(%s)); /* because y will get overwritten otherwise */\n",incy,nr,extra,ts)];
	else
		res=[res,findent,dump_vec(sx,ts,"x","! reference x \n",lang,mvec)];
		res=[res,findent,dump_vec(scy,ts,"cy","! reference cy after \n",lang,mvec)];
		sbcy=sy;
	#	res=[res,findent,dump_vec(sbcy,ts,"bcy","! reference bcy before \n",lang,mvec)];
		res=[res,findent,dump_vec(sy,ts,"y","! y \n",lang,mvec),"\n"];
		res=[res,findent,"y=x\n"];
	end 
end 

function res=print_matrix(a,lang,lab)
	ms="";
	if nargin<3
		lab="A";
	endif
#	if nargin<2
#		lang="c";
#	endif
	if lang == "c"
		hc="";
		co=["/*\n ",lab," = \n"];
		cc="*/\n";
		el="\n";
	elseif lang == "k"
		hc="";
		co=["\" ",lab," = \\n"];
		cc="\"";
		el="\\n";
	else # f
		hc="!";
		co=["! ",lab," =\n"];
		cc="";
		el="\n";
	endif
	res=co;
	for i=1:size(a,1)
	ms=[ms,hc];
	for j=1:size(a,2)
		if(iscomplex(a))
			ms=sprintf("%s %d+%di",ms,real(a(i,j)),imag(a(i,j)));
		else
			ms=sprintf("%s %d",ms,a(i,j));
		endif
	end
		ms=sprintf("%s%s",ms,el);
	end
	res=sprintf("%s%s%s",res,ms,cc);
end

function trc=blas_trans_char_array(tri)
	if tri==1
	trc="'n'";
	elseif tri==2
	trc="'t'";
	elseif tri==3
	trc="'c'";
	else
	trrc="'error :-)'";
	end
end

function tr=blas_trans_array(tri)
	if tri==1
	tr="blas_no_trans";
	elseif tri==2
	tr="blas_trans";
	elseif tri==3
	tr="blas_conj_trans ";
	end
end

function op=blas_types_array(tc,lang)
	if nargin == 1
		lang="c";
	endif
	if(lang=="c")
		if tc=='s'
		op="float";
		elseif tc=='d'
		op="double";
		elseif tc=='c'
		op="float complex";
		elseif tc=='z'
		op="double complex";
		end
	else
		if tc=='s'
		op="REAL*4";
		elseif tc=='d'
		op="REAL*8";
		elseif tc=='c'
		op="COMPLEX*8";
		elseif tc=='z'
		op="COMPLEX*16";
		end
	end
end

function s=pointer_symbol_if_type(tc)
	if tc=='s'
	s="";
	elseif tc=='d'
	s="";
	elseif tc=='c'
	s="&";
	elseif tc=='z'
	s="&";
	end
end

function op=blas_op_codes_array(oi)
	if oi==1
	op="usmv";
	elseif oi==2
	op="ussv";
	elseif oi==3
	op="usmm";
	else
	# error
	op="????"
	end
end

function res=blas_tester_function(what,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct)
	global blas_trans_codes_array;
	global matrix_types_array;
	global blas_type_codes_array;
	global findent;
	global matrix_diagonal;
	#
	n=rms;
	op=blas_op_codes_array(oi);
	bc=tc;#blas_type_codes_array(ti);
	tr=blas_trans_array(tri);
	#tc=blas_type_codes_array(ti);
	tn=blas_types_array(tc,lang);
	btc=blas_trans_codes_array(tri);
	trc=blas_trans_char_array(tri);
	mt=matrix_types_array(mti);
	md=matrix_diagonal(mdi);
	#
	res="";
	id="";
	#
	#
	if what=="id  "
		if alpha<1
			alphas=sprintf("nr%d",-alpha);
		else
			alphas=sprintf("p%d", alpha);
		endif
		if beta<1
			betas=sprintf("nr%d",-beta);
		else
			betas=sprintf("p%d", beta);
		endif
		res=sprintf("t%s_s%s_d%s_%s_%d_%c_a%s_b%s_ix%d_iy%d",tc,mt,md,op,rms,btc,alphas,betas,incx,incy);
	else
		id=blas_tester_function("id  ",lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy);
	endif
	#
	if what=="comm"
		res=sprintf("/* op:%s; type:%s; trans:%s kind:%s; diag:%s */", blas_op_codes_array(oi),bc,btc,mt,md);
	endif
	#
	if lang=="c"
	#
	if what=="CALL"
		res=sprintf("%s()",blas_tester_function("id  ",lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy));
	endif
	#
	if what=="decl"
		#bih="blas_invalid_handle";
		bih="-1";
		res=sprintf("static rsb_err_t %s(void)\n{\n\t%s\n",
		blas_tester_function("id  ",lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy),
		blas_tester_function("comm",lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy)
		);
		a=gen_test_matrix(op,n,tc,mti,mdi);
#		res=[res,"	rsb_err_t errval = RSB_ERR_NO_ERROR;\n"];
		res=[res,"	rsb_err_t errval = RSB_BLAS_ERROR;\n"];
		res=[res,"\tblas_sparse_matrix A = ",bih,";\n"];
		res=[res,sprintf("\tconst enum blas_trans_type transT=%s;\n",tr)];
		res=[res,sprintf("\tconst int incy=%d;\n",incy)];
		if op=="usmv"
			res=[res,sprintf("\tconst int incx=%d;\n",incx)];
		elseif op=="ussv"
			res=[res,sprintf("\tconst int incx=%d;\n",incx)];
		elseif op=="usmm"
			res=[res,sprintf("\tconst int ldX=%d;\n",n)]; # note: this only OK for square
			res=[res,sprintf("\tconst int ldY=%d;\n",n)]; # note: this only OK for square
		endif
		res=[res,sprintf("\tconst %s alpha=%d;\n",tn,alpha)];
		res=sprintf("%s\t%s",res,print_matrix(a,"c"));
		res=sprintf("%s\t/* declaration of VA,IA,JA */\n %s",res,dump_c_coo(a,tn,"c"));
		if ct=='v'
			res=[res,"	",dump_vec(ones(n,1),"rsb_coo_idx_t","K","// K\n",lang)];
			res=[res,"	",dump_vec(ones(n,1),"rsb_coo_idx_t","L","// L\n",lang)];
		endif
		mvec=1;
		if op=="usmv"
			res=[res,dump_csmm(a,mti,mdi,1,1,alpha,mvec,beta,tri,tn,incx,incy,lang)];
		elseif op=="usmm"
			res=[res,dump_csmm(a,mti,mdi,1,1,alpha,mvec,beta,tri,tn,incx,incy,lang)];
		else
			res=[res,dump_spsv(a,mdi,1,1,alpha,mvec,tri,tn,incx,incy,lang)];
		endif

		res=sprintf("%s\t\n\t%s%s%s%s%s\n",res,"if(!RSB_BLAS_SUPPORTED_TYPE('",bc,"')){printf(\"type=",bc," unsupported: skipping test.\\n\");errval=RSB_ERR_UNSUPPORTED_TYPE;goto err;}");
		res=sprintf("%s\t%s%s%s\n",res,"if((nnz == 0 ) && !RSB_BLAS_SUPPORT_EMPTY ){ printf(\"empty matrices are unsupported: skipping test.\\n\");errval=RSB_ERR_UNSUPPORTED_TYPE;goto err;}\n");
		if ct=='c'
			res=sprintf("%s\t%s%s%s\n",res,"A = BLAS_",bc,"uscr_begin(nr,nc);");
		endif
		if ct=='b'
			res=sprintf("%s\t%s%s%s\n",res,"A = BLAS_",bc,"uscr_block_begin(1,1,nr,nc);");
		endif
		if ct=='v'
			res=sprintf("%s\t%s%s%s\n",res,"A = BLAS_",bc,"uscr_variable_block_begin(nr,nc,K,L);");
		endif
#		gotoferrlabel="goto ferr;";
		gotoferrlabel="{RSB_ERROR(\"!\\n\");goto ferr;}";
		gotoferrlabel_pah="{RSB_ERROR(\"uscr_begin() gave %d!\\n\",A);goto ferr;}";
		gotoferrlabel_ussp="{RSB_ERROR(\"ussp() gave %d!\\n\",A);goto ferr;}";
		gotoferrlabel_inse="{RSB_ERROR(\"uscr_insert_entries() gave %d!\\n\",A);goto ferr;}";
		gotoferrlabel_end="{RSB_ERROR(\"uscr_end() gave %d!\\n\",A);goto ferr;}";
#		res=sprintf("%s\t%s%s\n",res,"if( A == blas_invalid_handle )",gotoferrlabel_pah);
		res=sprintf("%s\t%s%s\n",res,"if( A == -1 )",gotoferrlabel_pah);
		if op=="ussv"
			if mt=='u'
				res=[res,"	if( BLAS_ussp(A,blas_upper_triangular) != RSB_BLAS_NO_ERROR ) "];
			elseif mt=='l'
				res=[res,"	if( BLAS_ussp(A,blas_lower_triangular) != RSB_BLAS_NO_ERROR ) "];
			else
				res=[res,"	if( BLAS_ussp(A,blas_lower_triangular) != RSB_BLAS_NO_ERROR ) "];
			endif
			res=[res," ",gotoferrlabel_ussp,"\n"];
			#res=[res,"\tBLAS_ussp(A,blas_upper_triangular);\n"]; # NEW
		endif
		if md=='i'
			res=[res,"	if( BLAS_ussp(A,blas_unit_diag) != RSB_BLAS_NO_ERROR ) "];
			res=[res," ",gotoferrlabel_ussp,"\n"];
		endif
		if mt=='s'
			res=[res,"	if( BLAS_ussp(A,blas_lower_symmetric) != RSB_BLAS_NO_ERROR ) "];
			res=[res," ",gotoferrlabel_ussp,"\n"];
		endif
		if mt=='h'
			res=[res,"	if( BLAS_ussp(A,blas_lower_hermitian) != RSB_BLAS_NO_ERROR ) "];
			res=[res," ",gotoferrlabel_ussp,"\n"];
		endif
		res=sprintf("%s\t%s%s%s%s\n",res,"if( BLAS_",bc,"uscr_insert_entries(A,nnz,VA,IA,JA) != RSB_BLAS_NO_ERROR)",gotoferrlabel_inse);
		res=[res,sprintf("\t%s%s%s%s\n","if( BLAS_",bc,"uscr_end(A) != RSB_BLAS_NO_ERROR )",gotoferrlabel_end)];
	
		if op=="usmv"
			res=[res,sprintf("\tif( BLAS_%s%s(transT,%salpha,A,x,incx,y,incy) != RSB_BLAS_NO_ERROR )%s\n",bc,op,pointer_symbol_if_type(tc),gotoferrlabel)];
			res=[res,check_csmm(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang)];
		elseif op=="usmm"
			res=[res,sprintf("\tif( BLAS_%s%s(blas_colmajor,transT,1,%salpha,A,x,ldX,y,ldY) != RSB_BLAS_NO_ERROR )%s\n",bc,op,pointer_symbol_if_type(tc),gotoferrlabel)];
			res=[res,check_csmm(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang)];
		elseif op=="ussv"
			res=[res,sprintf("\tif( BLAS_%s%s(transT,%salpha,A,y,incx) != RSB_BLAS_NO_ERROR )%s\n",bc,op,pointer_symbol_if_type(tc),gotoferrlabel)];
			res=[res,check_spsv(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang)];
		endif
		res=[res,sprintf("%s%s\n","\n\tif( BLAS_usds(A) != RSB_BLAS_NO_ERROR )",gotoferrlabel)];
		cm=[check_message(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op)];
		nok=[cm," is not ok"];
#		res=[res,"	goto ok;\n"];
		res=[res,"	errval = RSB_ERR_NO_ERROR;\n"];
		res=[res,"ferr:\t"];
		res=[res,"if(errval != RSB_ERR_NO_ERROR){ "];
		res=[res," RSB_ERROR(\"",nok,"\\n\"); "];
		res=[res," RSB_ERRLN(lsc); "];
		res=[res," RSB_ERROR(\"Computed solution: y'=\\n\"); "];
		typecode=["'",toupper(tc),"'"];
		res=[res,"	rsb_sbtc_print_vec(y,nr*incy,",typecode,"); }\n"];
		res=[res,"err:\t"];
#		res=[res,sprintf("	%s","return RSB_ERR_NO_ERROR;\n\n	BLAS_usds(A);\n")];
		res=[res,sprintf("return errval;\n")];
#		res=[res,sprintf("ok:	return RSB_ERR_NO_ERROR;\n")];
		res=[res,sprintf("}\n")];
	endif
	#
	elseif lang=="f"
	#
		if(rand()>0.5)
			bpfx=sprintf("BLAS_%s",bc);
			bspfx=sprintf("BLAS_%s",bc);
			bipfx=bc;
		else
			bpfx="";
			bspfx=bc;
			bipfx="";
		endif
	#
		dprv=rand()*3;
		bdpfx="";
		if(dprv>1)
			bdpfx=sprintf("%s",bc);
		endif
		if(dprv>2)
			bdpfx=sprintf("BLAS_%s",bc);
		endif
	#
		if what=="CALL"
			res=[res,findent,"",id,""];
		endif
		if what=="decl"
			res=[res,findent,"SUBROUTINE ",id,"(errval)\n"];
			res=[res,findent,"USE blas_sparse\n"];
			res=[res,findent,"USE rsb ! for RSB_IDX_KIND\n"];
			res=[res,findent,"IMPLICIT NONE\n"];
			res=[res,findent,"INTEGER::errval,istat=0,i\n"];
			res=[res,findent,"INTEGER(KIND=RSB_IDX_KIND)::mone=-1,pone=+1\n"];
			res=[res,findent,"INTEGER::A\n"];
			res=[res,findent,"INTEGER::transT=",tr,"\n"];
			res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","incx",incx)];
			if op=="usmv"
			res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","incy",incy)];
			elseif op=="usmm"
			res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","incy",incy)];
			res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","ldX",n)];
			res=[res,findent,decl_var(lang,"INTEGER(KIND=RSB_IDX_KIND)","ldY",n)];
			endif
			res=[res,findent,decl_var(lang,tn,"alpha",alpha)];
			a=gen_test_matrix(op,n,tc,mti,mdi);
			cm=[check_message(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang)];
      			res=[res,findent,"CHARACTER(LEN=*),PARAMETER :: LCS = \"&\n",findent," &",cm,"\"","\n"];
			res=sprintf("%s%s\n",res,print_matrix(a,lang));
			res=sprintf("%s%s! declaration of VA,IA,JA \n%s\n",res,findent,dump_c_coo(a,tn,lang));
			if ct=='v'
				res=[res,findent,dump_vec(ones(n,1),"INTEGER(KIND=RSB_IDX_KIND)","K","! K\n",lang)];
				res=[res,findent,dump_vec(ones(n,1),"INTEGER(KIND=RSB_IDX_KIND)","L","! L\n",lang)];
			endif
			if( ( rand()*2>1 && incx==1 && incy==1 ) || op == "usmm" )
				mvec=1;
			else
				mvec=0;
			endif
			if op=="usmv"
				res=[res,dump_csmm(a,mti,mdi,1,1,alpha,mvec,beta,tri,tn,incx,incy,lang)];
			elseif op=="usmm"
				res=[res,dump_csmm(a,mti,mdi,1,1,alpha,mvec,beta,tri,tn,incx,incy,lang)];
			else
				res=[res,dump_spsv(a,mdi,1,1,alpha,mvec,tri,tn,incx,incy,lang)];
			endif
			res=[res,findent,"errval=0\n"];
			if ct=='c'
				res=[res,findent,sprintf("CALL %suscr_begin(nr,nc,A,errval)\n",bspfx)];
			endif
			if ct=='b'
				res=[res,findent,sprintf("CALL %suscr_block_begin(pone,pone,nr,nc,A,errval)\n",bspfx)];
			endif
			if ct=='v'
				res=[res,findent,sprintf("CALL %suscr_variable_block_begin(nr,nc,K,L,A,errval)\n",bspfx)];
			endif
			res=[res,findent,"IF(errval.NE.0)GOTO 9999\n"];
###############################################################################
			if op=="ussv"
				if mt=='u'
					res=[res,findent,"CALL ussp(A,blas_upper_triangular,istat)\n"]; # NEW
				elseif mt=='l'
					res=[res,findent,"CALL ussp(A,blas_lower_triangular,istat)\n"]; # NEW
				else
					res=[res,findent,"CALL ussp(A,blas_lower_triangular,istat)\n"]; # NEW
				endif
			endif
			if md=='i'
				res=[res,findent,"CALL ussp(A,blas_unit_diag,istat)\n"];
			endif
			if mt=='s'
				res=[res,findent,"CALL ussp(A,blas_lower_symmetric,istat)\n"];
			endif
			if mt=='h'
				res=[res,findent,"CALL ussp(A,blas_lower_hermitian,istat)\n"];
			endif
			res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
			if(rand()*3>1)
				res=[res,findent,sprintf("CALL %suscr_insert_row(A,mone,nnz,VA,JA,istat)\n",bpfx)];
				res=[res,findent,"IF(istat.EQ.0)THEN\n",findent,"  istat=-1;GOTO 9997;\n",findent,"ELSE\n",findent,"  istat=0;\n",findent,"ENDIF\n"];
			elseif(rand()*2>1)
				res=[res,findent,sprintf("CALL %suscr_insert_col(A,mone,nnz,VA,IA,istat)\n",bpfx)];
				res=[res,findent,"IF(istat.EQ.0)THEN\n",findent,"  istat=-1;GOTO 9997;\n",findent,"ELSE\n",findent,"  istat=0;\n",findent,"ENDIF\n"];
			endif
			if(rand()*2>1)
				res=[res,findent,sprintf("CALL %suscr_insert_entries(A,nnz,VA,IA,JA,istat)\n",bpfx)];
				res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
			else
				for i=1:nnz(a)
					if ct=='v'
						res=[res,findent,sprintf("CALL %suscr_insert_block(A,VA(%d:),pone,pone,IA(%d),JA(%d),istat)\n",bipfx,i,i,i)];
					else
					if(rand()*2>1)
						res=[res,findent,sprintf("CALL %suscr_insert_entry(A,VA(%d),IA(%d),JA(%d),istat)\n",bpfx,i,i,i)];
					else
						res=[res,findent,sprintf("CALL %suscr_insert_clique(A,pone,pone,VA(%d:),pone,pone,IA(%d:),JA(%d:),istat)\n",bspfx,i,i,i)];
					endif
					endif
					res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
				endfor
			endif
###############################################################################
			res=[res,findent,sprintf("CALL %suscr_end(A,istat)\n",bdpfx)];
			res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
			if op=="usmv"
				if(mvec>0)
					res=[res,findent,sprintf("CALL %susmm(blas_colmajor,transT,pone,alpha,A,x,nr,y,nc,istat)\n",bspfx),""];
				else
					res=[res,findent,sprintf("CALL %s%s(transT,alpha,A,x,incx,y,incy,istat)\n",bpfx,op),""];
				endif
				res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
				res=[res,check_csmm(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang,mvec)];
			elseif op=="usmm"
				res=[res,findent,sprintf("CALL %s%s(blas_colmajor,transT,1,alpha,A,x,ldX,y,ldY,istat)\n",bpfx,op),""];
				res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
				res=[res,check_csmm(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang,mvec)];
			elseif op=="ussv"
				if(mvec>0)
					res=[res,findent,sprintf("CALL %sussm(blas_colmajor,transT,pone,alpha,A,y,nr,istat)\n",bspfx),""];
				else
					res=[res,findent,sprintf("CALL %s%s(transT,alpha,A,y,incx,istat)\n",bpfx,op)];
				endif
				res=[res,findent,"IF(istat.NE.0)GOTO 9997\n"];
				res=[res,check_spsv(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang,mvec)];
			endif
#			res=[res,findent,"GOTO 9998\n"];
#			res=[res,"9997",findent,"errval=-1\n"];
			res=[res,"9997",findent,"IF(istat.NE.0)errval=-1\n"];
#			res=[res,"9998",findent,"CONTINUE\n"];
			res=[res,findent,"CALL usds(A,istat)\n"];
			res=[res,findent,"IF(istat.NE.0)errval=-1\n"];
			res=[res,"9999",findent,"CONTINUE\n"];
			res=[res,findent,"end SUBROUTINE ",id," \n"];
		endif
	#
	elseif lang=="p"
	#
		
		if what=="CALL"
			res=[res,findent,"",id,""];
		endif
		if what=="decl"
			a=gen_test_matrix(op,n,tc,mti,mdi);
			cm=[check_message(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op)];
			#pro=[cm," has problems"];
			nok=[cm," is not ok"];
			pro=nok;
			ok=[cm," is ok"];
			#
			res=[res,findent,"SUBROUTINE ",id,"(errval,afmt,ictxt)\n"];
			res=[res,findent,"USE psb_base_mod\n"];
			res=[res,findent,"IMPLICIT NONE\n"];
			res=[res,findent,"CHARACTER(LEN=*) :: afmt\n"];
#			res=[res,findent,"CHARACTER(LEN=psb_fidasize_) :: afmt\n"];
			res=[res,findent,"TYPE(psb_",bc,"spmat_type) :: a\n"];
			res=[res,findent,"TYPE(psb_desc_type)   :: desc_a\n"];
			res=[res,findent,"INTEGER            :: ictxt, iam=-1, np=-1\n"];
			res=[res,findent,"INTEGER            :: info=-1\n"];
#			res=[res,findent,"INTEGER   :: idim\n"];
			res=[res,findent,"\n"];
			res=[res,findent,"INTEGER::errval,istat=0,i\n"];
			res=[res,findent,"CHARACTER::transA=",trc,"\n"];
			res=[res,findent,decl_var(lang,"INTEGER","incx",incx)];
			if op=="usmv"
			res=[res,findent,decl_var(lang,"INTEGER","incy",incy)];
			elseif op=="usmm"
			res=[res,findent,decl_var(lang,"INTEGER","incy",incy)];
			endif
			res=[res,findent,decl_var(lang,tn,"alpha",alpha)];
			res=[res,findent,decl_var(lang,tn,"beta",beta)];
			res=sprintf("%s%s\n",res,print_matrix(a,lang));
			res=sprintf("%s%s! declaration of VA,IA,JA \n%s\n",res,findent,dump_c_coo(a,tn,lang));
			mvec=0;
			if op=="usmv"
				res=[res,dump_csmm(a,mti,mdi,1,1,alpha,mvec,beta,tri,tn,incx,incy,lang)];
			elseif op=="usmm"
				res=[res,dump_csmm(a,mti,mdi,1,1,alpha,mvec,beta,tri,tn,incx,incy,lang)];
			else
				res=[res,dump_spsv(a,mdi,1,1,alpha,mvec,tri,tn,incx,incy,lang)];
			endif
			res=[res,findent,"errval=0\n"];
#			res=[res,findent,"afmt='CSR'\n"];
#			res=[res,findent,"CALL psb_init(ictxt)\n"];
			res=[res,findent,"CALL psb_info(ictxt,iam,np)\n"];
			res=[res,findent,"IF(iam<0)THEN\n"];
			res=[res,findent,findent,"info=-1\n"];
			res=[res,findent,findent,"GOTO 9999\n"];
			res=[res,findent,"ENDIF\n"];
			res=[res,findent,"CALL psb_barrier(ictxt)\n"];
			res=[res,findent,"CALL psb_cdall(ictxt,desc_a,info,nl=nr)\n"];
			res=[res,findent,"IF (info .NE. 0)GOTO 9996\n"];
			res=[res,findent,"CALL psb_spall(a,desc_a,info,nnz=nnz)\n"];
			res=[res,findent,"IF (info .NE. 0)GOTO 9996\n"];
			if op=="ussv"
			res=[res,findent,"a%descra='TLN'\n"];
			endif
			res=[res,findent,"CALL psb_barrier(ictxt)\n"];
			res=[res,findent,"CALL psb_spins(nnz,IA,JA,VA,a,desc_a,info)\n"];
			res=[res,findent,"IF (info .NE. 0)GOTO 9996\n"];
			res=[res,findent,"CALL psb_cdasb(desc_a,info)\n"];
			res=[res,findent,"IF (info .NE. 0)GOTO 9996\n"];
			res=[res,findent,"CALL psb_spasb(a,desc_a,info,dupl=psb_dupl_err_,afmt=afmt)\n"];
			res=[res,findent,"IF(info.NE.0)PRINT *,\"matrix assembly failed\"\n"];
			res=[res,findent,"IF(info.NE.0)GOTO 9996\n"];
			res=[res,findent,"\n"];
#			res=[res,findent,sprintf("CALL %suscr_begin(nr,nc,A,errval)\n",bc)];
#			res=[res,findent,"IF(errval.NE.0)GOTO 9999\n"];
#			res=[res,findent,"CALL uscr_insert_entries(A,nnz,VA,IA,JA,istat)\n"];
#			res=[res,findent,"IF(istat.NE.0)GOTO 9996\n"];
#			res=[res,findent,"CALL uscr_end(A,istat)\n"];
#			res=[res,findent,"IF(istat.NE.0)GOTO 9996\n"];
			if op=="usmv"
				res=[res,findent,sprintf("CALL psb_spmm(alpha,A,x,beta,y,desc_a,info,transA)\n",op),""];
				res=[res,findent,"IF(info.NE.0)PRINT *,\"psb_spmm failed\"\n"];
				res=[res,findent,"IF(info.NE.0)GOTO 9996\n"];
				res=[res,check_csmm(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang)];
			elseif op=="usmm"
				;
			elseif op=="ussv"
#				res=[res,findent,"x(:)=y(:)\n"];
#				res=[res,findent,"beta=0\n"];
#				res=[res,findent,sprintf("CALL psb_spsm(alpha,A,x,beta,y,desc_a,info,transA,diag=x)\n",op)];
				res=[res,findent,sprintf("CALL psb_spsm(alpha,A,x,beta,y,desc_a,info,transA)\n",op)];
				res=[res,findent,"IF(info.NE.0)PRINT *,\"psb_spsm failed\"\n"];
				res=[res,findent,"IF(info.NE.0)GOTO 9996\n"];
				res=[res,check_spsv(a,mti,mdi,1,1,alpha,beta,tri,incx,incy,tc,op,lang)];
			endif
#			res=[res,findent,"GOTO 9998\n"];
#			res=[res,"9996",findent,"errval=-1\n"];
			res=[res,"9996",findent,"CONTINUE\n"];
			res=[res,"",findent,"IF(info .NE. 0)errval=errval+1\n"];
			res=[res,"",findent,"CALL psb_spfree(a,desc_a,info)\n"];
			res=[res,findent,"IF (info .NE. 0)GOTO 9997\n"];
			res=[res,"9997",findent,"CONTINUE\n"];
			res=[res,"",findent,"IF(info .NE. 0)errval=errval+1\n"];
			res=[res,"",findent,"CALL psb_cdfree(desc_a,info)\n"];
			res=[res,findent,"IF (info .NE. 0)GOTO 9998\n"];
			res=[res,"9998",findent,"CONTINUE\n"];
			res=[res,"",findent,"IF(info .NE. 0)errval=errval+1\n"];
#			res=[res,"9998",findent,"CONTINUE\n"];
#			res=[res,findent,"CALL usds(A,istat)\n"];
#			res=[res,findent,"IF(istat.NE.0)errval=-1\n"];
#			res=[res,"9999",findent,"CONTINUE\n"];
#			res=sprintf("%s%s",res,findent,"CALL psb_exit(ictxt)\n");
			res=[res,"9999",findent,"CONTINUE\n"];
			res=[res,"",findent,"IF(info .NE. 0)errval=errval+1\n"];
			res=[res,findent,findent,"IF(errval.NE.0)PRINT*,\"",pro,"\"\n"];
			res=[res,findent,findent,"IF(errval.EQ.0)PRINT*,\"",ok,"\"\n"];
			res=[res,findent,"END SUBROUTINE ",id," \n"];
		endif
	end
end

function all_test(lang,what)
global blas_op_codes_num;
global max_random_matrix_size_a;
global blas_type_codes_array;
global blas_trans_codes_array;
global alpha_array;
global beta_array;
global incx_array;
global incy_array;
global findent;
global matrix_types_array;
global matrix_diagonal;
global blas_ctor_funcs;
#
if what=="decl"
	what_here="decl";
else
	what_here="CALL";
endif
#
if what_here=="CALL"
if lang=="c"
	res="";
	res=[res,sprintf("	printf(\"Test print BEGIN\\n\");\n")];
	for ti=1:length(blas_type_codes_array)
		tc=blas_type_codes_array(ti);
		tn=blas_types_array(tc,lang);
		res=sprintf("%s\t%s"  ,res,"{\n");
		res=sprintf("%s\t const %s v [] = {1,2,3};\n%s",res,tn);
		res=sprintf("%s\t rsb_sbtc_print_vec(v,sizeof(v)/sizeof(v[0]),toupper('%s'));\n",res,tc);
		res=sprintf("%s\t%s"  ,res,"}\n");
	endfor
	res=[res,sprintf("	printf(\"Test print END\\n\");\n")];
	printf("%s",res);
endif
endif
#
for cfi=1:length(blas_ctor_funcs)
#for mdi=2:length(matrix_diagonal)
for mdi=1:length(matrix_diagonal)
for mti=1:length(matrix_types_array)
for ti=1:length(blas_type_codes_array)
tc=blas_type_codes_array(ti);
for oi=1:blas_op_codes_num
#for oi=2:2
#for oi=1:1
for rmsi=1:length(max_random_matrix_size_a)
rms=max_random_matrix_size_a(rmsi);
for alphai=1:length(alpha_array)
for betai=1:length(beta_array)
for incxi=1:length(incx_array)
for incyi=1:length(incy_array)
incy=incy_array(incyi);
incx=incx_array(incxi);
for tri=1:length(blas_trans_codes_array)
res="";
beta=beta_array(betai);
alpha=alpha_array(alphai);
op=blas_op_codes_array(oi);
if ( op == "usmm" && ( beta != 1 || incx != 1 || incy != 1  ) ) continue; end # only incx==1 incy==1 make sense for usmm
mt=matrix_types_array(mti);
md=matrix_diagonal(mdi);
ct=blas_ctor_funcs(cfi);
ctrv=rand()*3;
if(ctrv>1)
	ct='c';
else
	ct='b';
endif
if(ctrv>2)
	ct='v';
endif
#
#op
#mt
if !xor( (mt=='l' || mt=='u'),( strcmp(op,"usmv")==0) ) # FIXME: for some reason, op=="usmv" is not as op=="ussv"
##if true
#if (((mt=='l') || (mt=='u')) && ( op == "ussv")) || (((mt!='l') && (mt!='u')) && ( op != "ussv")) 
#	res=[res,sprintf("/* %s %s */\n",op,mt)];
#if (mt=='l' || mt=='u' ) && op=="ussv"
#if 1
#if (mt=='l' || mt=='u') && (op != "ussv") ; continue ; endif
#
if !(lang=="p" && (incx!=1 || incy!=1 || (beta!=0 && op=="ussv") || mt!='g' || md!='e'))
#
if (!(op=="ussv" && lang!="p"  && (beta!=1 || incy!=incx))) && (!(op=="usmv" && lang!="p"  && (beta==0 ))) && (!(op=="ussv" && (mt=='s' || mt=='h' || mt=='g')))
#if !(op=="ussv" && lang!="p"  && (beta!=1 || tri>1 || incy!=incx)) && !(op=="usmv" && lang!="p"  && (beta==0 ))
#
if lang=="c"
#
if what_here=="CALL"
	fid=blas_tester_function(what_here,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct);
	res=sprintf("%s\t errval = %s;\n",res,fid);
#	res=sprintf("%s\tif( errval != RSB_ERR_NO_ERROR )++failed;else++passed;\n",res);
res=sprintf("%s\tif( errval== RSB_ERR_NO_ERROR )++passed;else{if(errval==RSB_ERR_UNSUPPORTED_TYPE)++skipped,errval=RSB_ERR_NO_ERROR ;else++failed;}\n",res);
	res=sprintf("%s\tif( errval != RSB_ERR_NO_ERROR )RSB_ERROR(\"%s failed!\\n\");\n",res,fid);
else
	res=sprintf("%s",sprintf("%s\t%s",res,blas_tester_function(what_here,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct)));
	res=sprintf("%s\n",res);
endif
#
#
elseif lang=="f"
#
if what_here=="CALL"
	res=[res,findent,"CALL ",blas_tester_function(what_here,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct),"(errval)\n"];
	res=[res,findent,"IF(errval.LT.0)failed=failed+1\n"];
	res=[res,findent,"IF(errval.EQ.0)passed=passed+1\n"];
	res=[res,findent,"\n"];
else
#	res=[res,"! declaration ... \n"];
	res=[res,"! \n"];
	res=[res,"",blas_tester_function(what_here,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct),""];
#	res=[res,findent,"! TODO: still unimplemented \n"];
endif
#
elseif lang=="p"
#
if what_here=="CALL"
	res=[res,findent,"CALL ",blas_tester_function(what_here,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct),"(errval,afmt,ictxt)\n"];
	res=[res,findent,"IF(errval.NE.0)failed=failed+1\n"];
	res=[res,findent,"IF(errval.EQ.0)passed=passed+1\n"];
	res=[res,findent,"errval=0\n"];
	res=[res,findent,"\n"];
else
#	res=[res,"! declaration ... \n"];
	res=[res,"! \n"];
	res=[res,"",blas_tester_function(what_here,lang,mti,mdi,tc,oi,rms,tri,alpha,beta,incx,incy,ct),""];
#	res=[res,findent,"! TODO: still unimplemented \n"];
endif
#
endif
	printf("%s",res);
#
end
end
end
end
end
end
end
end
end
end
end
end
end
end
#
#
end # end all_test function

function lt=rsb_octave_license(lang);
#
pre="";
if lang == "f" 
	pre="! ";
end 
lt="";
fd=fopen("rsb_license_header.inc","r");
while (ll=fgets(fd,1024)) != -1 ;
lt=sprintf("%s%s%s",lt,pre,ll);
endwhile;
fclose(fd);
#
end 


