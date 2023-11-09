/*	Solution of a system of linear algebraic equations using,
	Gaussian elimination, Crout's triangular decomposition,
	Crout's decomposition with iterative refinement,
	   This will work only if long double has higher precision than double
	Cholesky's triangular decomposition for a real symmetric matrix,
	Singular value decomposition */

/*	CHOLSK will work only for symmetric matrices, while other functions
	should work for a general matrix */

#include <stdio.h>
#include <math.h>

int cholsk(int n, int num, double * a, double *x, double *det, int nd,
	int *iflg);
int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg);
int crouth(int n, int num, double *a, double *b, double *x, double *det,
	int *idet, int inc[], int lj, double reps, int *iflg, double wk[]);
int gauelm(int n, int num, double *a, double *x, double * det,
       int inc[], int lj, int * iflg);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);

main()
{
	int i,i1,it,j,n, num, iflg, idet, lj, inc[20];
	double r1, det, reps, x[2][20], wk[60], am[20][20], bm[20][20], sigma[20];

/*	Example 3.5 : Hilbert matrix */

	num=1;
	lj=20;
	for(i1=0; i1<99; ++i1) {
		printf("type n = no. of equations         (quits when n<=0)\n");
		scanf("%d", &n);
		if(n<=0) return 0;
		printf("type it=1/2/3/4/5  for gauelm/crout/crouth/cholsk/svd\n");
		scanf("%d", &it);

/*	Generating the n x n Hilbert matrix and the right hand side */
		for(i=0; i<n; ++i) {
			r1=0.0;
			for(j=0; j<n; ++j) {
				am[i][j]=1.0/(i+j+1.0);
				r1=r1+am[i][j];
			}
			x[0][i]=r1;
		}

		iflg=0;
		switch(it) {
			case 1: {
				i=gauelm(n,num,&am[0][0],&x[0][0],&det,inc,lj,&iflg);
				printf("gauelm:     %d    n = %d    det = %e \n",i,n,det);
				break;
			}
			case 2: {
				i=crout(n,num,&am[0][0],&x[0][0],&det,&idet,inc,lj,&iflg);
				printf("crout:     %d    n = %d    det = %e  %d \n",i,n,det,idet);
				break;
			}
			case 3: {
				reps=1.e-14;
				i=crouth(n,num,&am[0][0],&bm[0][0],&x[0][0],&det,&idet,inc,lj,reps,&iflg,wk);
				printf("crouth:     %d    n = %d    det = %e  %d \n",i,n,det,idet);
				printf("error estimate = %e \n",wk[0]);
				break;
			}
			case 4: {
				i=cholsk(n,num,&am[0][0],&x[0][0],&det,lj,&iflg);
				printf("cholsk :     %d    n = %d    det = %e \n",i,n,det);
				printf("error estimate = %e \n",wk[0]);
				break;
			}
			default : {
				i=svd(n,n,&am[0][0],&bm[0][0],sigma,lj,lj);
				printf("svd:     %d    n = %d    singular values =  \n",i,n);
				for(j=0; j<n; ++j) printf(" %e ",sigma[j]);
				printf("\n");
				reps=1.e-14;
				i=svdevl(n,n,&am[0][0],&bm[0][0],sigma,lj,lj,&x[0][0],reps);
				break;
			}
		}
		printf(" solution :");
		for(j=0; j<n; j++) printf(" %f ",x[0][j]);
		printf(" \n");
	}
        return;
}




/*     Solution of a system of linear equations with real symmetric
            positive definite matrix using Cholesky decomposition

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
        	equations to be solved
	A : (input/output) The matrix of coefficient of size ND*N
        	A[i][j] is the coefficient of x_j in ith equation
		at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size ND*NUM)
        	X[j][i] is the ith element of jth right hand side
		at output it will contain the solutions
	DET : (output) The determinant of the matrix
	ND : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter which determines the type of computation
		required.
        	If IFLG<=0, both elimination and solution are calculated
	    		and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
			decomposition should have been calculated earlier
		
	Error status is returned by the value of the function CHOLSK.
		0 value implies successful execution
		103 implies (N<=0 or N>ND)
		123 implies some pivot turned out to be zero

	Required functions : None
*/
 
#include <math.h>

int cholsk(int n, int num, double * a, double *x, double *det, int nd,
	int *iflg)

{
	int i,j,k,l;
	double sum;

	if(n<=0 || n>nd) return 103;

	if(*iflg < 2) {
/*     Perform elimination  */
		*det = 1.0;
		for(k=0; k<n; ++k) {
			for(i=0; i<=k-1; ++i) {
				if( a[i*nd+i] == 0.0) return 123;
				sum= a[k*nd+i];
				for(j=0; j<=i-1; ++j) sum=sum-a[i*nd+j]*a[k*nd+j];
				a[k*nd+i]=sum/a[i*nd+i];
			}

			sum = a[k*nd+k];
			for(j=0; j<=k-1; ++j) sum=sum-a[k*nd+j]*a[k*nd+j];
			if(sum<= 0.0) return 123;
			a[k*nd+k]=sqrt(sum);
			*det = (*det)*sum;
		}

		if(*iflg == 1) {*iflg = 2; return 0;}
		*iflg = 2;
	}
/*     Solution for the NUM different right-hand sides */
	for(j=0; j<num; ++j) {
/*     Forward substitution  */
		x[j*nd]= x[j*nd]/a[0];
		for(k=1; k<n; ++k) {
			for(l=0; l<=k-1; ++l) x[j*nd+k]=x[j*nd+k]-a[k*nd+l]*x[j*nd+l];
			x[j*nd+k]=x[j*nd+k]/a[k*nd+k];
		}
/*     back-substitution */
		x[j*nd+n-1]=x[j*nd+n-1]/a[(n-1)*nd+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*nd+k]=x[j*nd+k]-x[j*nd+l]*a[l*nd+k];
			x[j*nd+k] = x[j*nd+k]/a[k*nd+k];
		}
	}
	return 0;
}



 
/*	Solution of a system of linear equations using Crout's algorithm 
	with partial pivoting

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input/output) The matrix of coefficient of size LJ*N
	        A[i][j] is the coefficient of x_j in ith equation
			at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
			at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
			interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter which determines the type of computation
		required.
		If IFLG<=0, both elimination and solution are calculated
			and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function CROUT
		0 value implies successful execution
		102 implies (N<0 or N>LJ) 
		122 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg)

{
	int i,j,k,km,l;
	double r1,t1,r2,d1;
	double *wk;

	if(n<=0 || n>lj) return 102;

	if((*iflg)<2)
/*	Perform LU decomposition  */
	{wk=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i){
		r1=0.0;
		for(j=0; j<n; ++j){
			if(fabs(a[i*lj+j]) > r1) r1=fabs(a[i*lj+j]);}
		wk[i]=r1;
		if(r1==0.0) {free(wk); return 122;}
	}

	*det=1.0; *idet=0;
	for(k=0; k<n; ++k){
		r1=0.0; km=k;
/*	Generate the Kth column of L  */
		for (l=k; l<n; ++l){
			d1=a[l*lj+k];
			for(i=0; i<=k-1; ++i) d1=d1-a[l*lj+i]*a[i*lj+k];
			a[l*lj+k]=d1;
/*	Finding the pivot  */
			r2=fabs(d1/wk[l]);
			if(r2>r1) {r1=r2; km=l;}
		}

		inc[k]=km;
/*	Interchange the rows if needed  */
		if(km != k) {
			*det=-(*det);
			for(l=0; l<n; ++l) {
				t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;
			}
			t1=wk[k]; wk[k]=wk[km]; wk[km]=t1;
		}

		*det=*det*a[k*lj+k];
		if(a[k*lj+k]==0.0) {free(wk); return 122;}

		if(*det != 0.0) {
/*	Scale the value of the determinant DET */
			while(fabs(*det)>32.0) {
				*det=(*det)*0.03125; *idet=*idet+5;}

			while(fabs(*det)<0.03125) {
				*det=(*det)*32.0; *idet=*idet-5;}
		}

/*	Generate the Kth row of U */
		for(l=k+1; l<n; ++l) {
			d1=a[k*lj+l];
			for(i=0; i<=k-1; ++i) d1=d1-a[k*lj+i]*a[i*lj+l];
			a[k*lj+l]=d1/a[k*lj+k];
		}
	}
	free(wk);

	if(*iflg==1) {*iflg=2; return 0;}
	*iflg=2;
	}

/*	Solution for NUM different right-hand sides  */
	for(j=0; j<num; ++j){
/*	Forward substitution  */
		for(k=0; k<n; ++k) {
			if(k != inc[k]) {
				t1=x[j*lj+k];
				x[j*lj+k]=x[j*lj+inc[k]];
				x[j*lj+inc[k]]=t1;
			}

			d1=x[j*lj+k];
			for(l=0; l<=k-1; ++l) d1=d1-a[k*lj+l]*x[j*lj+l];
			x[j*lj+k]=d1/a[k*lj+k];
		}

/*	Back-substitution  */
		for(k=n-2; k>=0; --k){
			d1=x[j*lj+k];
			for(l=n-1; l>=k+1; --l) d1=d1-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=d1;
		}
	}
	return 0;
}



/*	Solution of a system of linear equations using Crout's algorithm 
            with iterative refinement

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input) The matrix of coefficient of size LJ*N
	         A[i][j] is the coefficient of x_j in ith equation
	B : (output) Array of size LJ*N containing triangular decomposition of matrix A
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
			at output it will contain the solutions
	DET, IDET : (output) The determinant of the matrix = DET*2**IDET
	INC : (output) Integer array of length N containing information about
			interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	REPS : (input) Required relative precision in solution
	IFLG : (input) Integer parameter which determines the type of computation
		required.
		If IFLG<=0, both elimination and solution are calculated
			and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
	WK : (output) Array of length NUM containing the estimated errors.
		WK[I] will contain estimated error in the Ith solution.
		
	Error status is returned by the value of the function CROUTH.
		0 value implies successful execution
		11 implies that iterative refinement did not converge
		102 implies (N<=0 or N>LJ) 
		122 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : CROUT
*/

#include <math.h>
#include <stdlib.h>

int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg);

int crouth(int n, int num, double *a, double *b, double *x, double *det,
	int *idet, int inc[], int lj, double reps, int *iflg, double wk[])

{
	int i,j,k,it,iflg1,ier,num1, nitr=10;
	double r1,r2,rp1,rp;
/* If long double is supported, use that for residuals */
	long double d1,d2;
	double *wk1, *wk2;

	if(n<=0 || n>lj) return 102;

	if(*iflg<2) {
/*	Preserving the matrix for calculating the residuals */
		for(i=0; i<n; ++i) {
			for(j=0; j<n; ++j) b[j*lj+i]= a[j*lj+i];
		}
		iflg1=1;
/*	Perform LU decomposition using crout */
		ier=crout(n,num,b,x,det,idet,inc,lj,&iflg1);
		if(ier>100) return ier;
	}

	if(*iflg==1) {*iflg=2; return 0;}
	*iflg=2;
	num1=1;
	wk1=(double *) calloc((size_t) n, sizeof(double));
	wk2=(double *) calloc((size_t) n, sizeof(double));

/*	Solving the systems with NUM different right-hand sides */
	for(j=0; j<num; ++j) {
		for(i=0; i<n; ++i) {
			wk1[i]= x[j*lj+i];
/*	Preserving the RHS for calculating residuals */
			wk2[i]= x[j*lj+i];
			x[j*lj+i]=0.0;
		}

		rp1=0.0;
/*	The iterative refinement */
		for(it=1; it<= nitr ; ++it) {
			ier=crout(n,num1,b,wk1,det,idet,inc,lj,&iflg1);
			if(ier>100) {free(wk2); free(wk1); return ier;}
			r1=0.0; r2=0.0;
			for(i=0; i<n; ++i) {
				if(fabs(wk1[i])>r1) r1=fabs(wk1[i]);
				x[j*lj+i]=x[j*lj+i]+wk1[i];
				if(fabs(x[j*lj+i])>r2) r2=fabs(x[j*lj+i]);
			}

/*	The error estimate  */
			if(it==2) wk[j]=r1;
			if(r2==0.0) goto nextit;
			rp=r1/r2;
			if(rp<reps) goto nextit;
			if(rp>rp1 && it>1) {ier=11; goto nextit;}
			rp1=rp;

/*	Calculating the residue  */
			for(i=0; i<n; ++i) {
				d1= wk2[i];
				for(k=0; k<n; ++k) {
					d2=a[i*lj+k];
					d1=d1-d2*x[j*lj+k];
				}
				wk1[i]=d1;
			}
		}
		ier=11;
nextit: ;
	}
	free(wk2); free(wk1);
	return ier;
}





/*	Solution of a system of linear equations using Gaussian elimination
	with partial pivoting

	N : (input) Number of equations to be solved
	NUM : (input) Number of different sets (each with N equations) of
	         equations to be solved
	A : (input/output) The matrix of coefficient of size LJ*N
	        A[i][j] is the coefficient of x_j in ith equation
	     	at output it will contain the triangular decomposition
	X : (input/output) The matrix containing right hand sides (size LJ*NUM)
	        X[j][i] is the ith element of jth right hand side
	     	at output it will contain the solutions
	DET : (output) The determinant of the matrix
	INC : (output) Integer array of length N containing information about
		interchanges performed during elimination
	LJ : (input) Second dimension of arrays A and X in calling function
	IFLG : (input) Integer parameter to specify the type of computation required
		If IFLG<=0, both elimination and solution are
			done and IFLG is set to 2
		If IFLG=1, only elimination is done and IFLG is set to 2
		If IFLG>=2 only solution is calculated, the triangular
		    decomposition should have been calculated earlier
		
	Error status is returned by the value of the function GAUELM.
		0 value implies successful execution
		101 implies (N<=0 or N>LJ) 
		121 implies some pivot turned out to be zero and hence
			matrix must be nearly singular

	Required functions : None
*/

#include <math.h>

int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg)

{
int i,j,k,km,l;
double r1,t1;

	if(n<=0 || n>lj) return 101;
 
	if((*iflg)<2) {
/*	Perform elimination  */
		*det=1.0;
		for(k=0; k<n-1; ++k) {
/*	Find the maximum element in the Kth column  */
			r1=0.0; km=k;
			for(l=k; l<n; ++l)
				if(fabs(a[l*lj+k])>r1) {r1=fabs(a[l*lj+k]); km=l;}

			inc[k]=km;
			if(km != k) {
/*	Interchange the rows if needed  */
				for(l=k; l<n; ++l) 
				{t1=a[k*lj+l]; a[k*lj+l]=a[km*lj+l]; a[km*lj+l]=t1;}
				*det=-(*det);
			}

			*det=(*det)*a[k*lj+k];
			if(a[k*lj+k]==0) return 121;

			for(l=k+1; l<n; ++l) {
				a[l*lj+k]=a[l*lj+k]/a[k*lj+k];
				for(i=k+1; i<n; ++i) a[l*lj+i]=a[l*lj+i]-a[l*lj+k]*a[k*lj+i];
			}
		}
		*det=(*det)*a[(n-1)*lj+n-1];
		inc[n-1]=n-1;
		if(a[(n-1)*lj+n-1]==0) return 121;

		if((*iflg)==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solution for the num different right-hand sides  */
	for(j=0; j<num; ++j) {
/*	forward-substitution  */
		for(k=0; k<n-1; ++k) {
			if(k != inc[k])
			{t1=x[j*lj+k]; x[j*lj+k]=x[j*lj+inc[k]]; x[j*lj+inc[k]]=t1;}
			for(l=k+1; l<n; ++l) x[j*lj+l]=x[j*lj+l]-a[l*lj+k]*x[j*lj+k];
		}

/*	back-substitution  */

		x[j*lj+n-1]=x[j*lj+n-1]/a[(n-1)*lj+n-1];
		for(k=n-2; k>=0; --k) {
			for(l=n-1; l>=k+1; --l) x[j*lj+k]=x[j*lj+k]-x[j*lj+l]*a[k*lj+l];
			x[j*lj+k]=x[j*lj+k]/a[k*lj+k];
		}
	}
	return 0;
}



/*	To calculate the Singular Value Decomposition of a matrix A=U D Vtranspose

	N : (input) Number of variables
	M : (input) Number of equations
	A : (input/output) Matrix of coefficients of size LA*M
		After execution it will contain the matrix U
	V : (output) The matrix V of size LV*N
	SIGMA : (output) Array of length N, containing the singular values
	LA : (input) Actual value of second dimension of A in the calling function
	LV : (input) Actual value of second dimension of V in the calling function
		
	Error status is returned by the value of the function SVD.
		0 value implies successful execution
		12 QR iteration failed to converge to required accuracy
		105 implies N<=0, N>LV, M<=0, N>LA, N>M

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv)

{
	int i,j,k,l,itr,ier, itmax=30;
	double f, g, h, rmax, s, r1, r2, c, x, y, z, aeps, eps=1.e-16;
	double *e;

	if(n>m || n<=0 || m<=0 || n>la || n>lv) return 105;
	ier=0;

/*	Reduction to Bidiagonal form using Householder transformations */
	g=0.0; rmax=0.0;
	e=(double *) calloc((size_t) n, sizeof(double));

	for(i=0; i<n; ++i) {
/*	Off-diagonal elements of bidiagonal form  */
		e[i]=g;
		s=0.0;
		for(j=i; j<m; ++j) s=s+a[i+j*la]*a[i+j*la];
		if(s <= 0.0) {
/*	transformation not required */
			g=0.0;
		}
		else {
			f= a[i+i*la];
			g=sqrt(s);
			if(f>=0.0) g=-g;
			h=f*g-s;
			a[i+i*la] = f-g;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]= a[j+k*la]+f*a[i+k*la];
			}
		}

/*	Diagonal elements of bidiagonal form  */
		sigma[i]=g;
		s=0.0;
		for(j=i+1; j<n; ++j) s=s+a[j+i*la]*a[j+i*la];

		if(s<= 0.0) g=0.0;
		else {
			f= a[i*la+(i+1)];
			g=sqrt(s);
			if(f>= 0.0) g=-g;
			h=f*g-s;
			a[i*la+(i+1)]=f-g;
			for(j=i+1; j<n; ++j) e[j]=a[j+i*la]/h;

			for(j=i+1; j<m; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+j*la]*a[k+i*la];
				for(k=i+1; k<n; ++k) a[k+j*la] = a[k+j*la]+s*e[k];
			}
		}
		r1=fabs(sigma[i])+fabs(e[i]);
		if(r1 > rmax) rmax=r1;
	}

/*	Accumulation of right hand transformation in array V */
	for(i=n-1; i>=0; --i) {
		if(g != 0.0) {
			h=g*a[i*la+(i+1)];
			for(j=i+1; j<n; ++j) v[i+j*lv]=a[j+i*la]/h;

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<n; ++k) s=s+a[k+i*la]*v[j+k*lv];
				for(k=i+1; k<n; ++k) v[j+k*lv]=v[j+k*lv]+s*v[i+k*lv];
			}
		}

		for(j=i+1; j<n; ++j) {
			v[j+i*lv]=0.0; v[i+j*lv]=0.0;
		}
		v[i+i*lv]=1;
		g= e[i];
	}

/*	Accumulation of left hand transformation overwritten on matrix A */
	for(i=n-1; i>=0; --i) {
		g=sigma[i];
		for(j=i+1; j<n; ++j) a[j+i*la]=0.0;
		if(g != 0.0) {
			h=g*a[i+i*la];

			for(j=i+1; j<n; ++j) {
				s=0.0;
				for(k=i+1; k<m; ++k) s=s+a[i+k*la]*a[j+k*la];
				f=s/h;
				for(k=i; k<m; ++k) a[j+k*la]=a[j+k*la]+f*a[i+k*la];
			}

			for(j=i; j<m; ++j) a[i+j*la]=a[i+j*la]/g;
		}
		else {
			for(j=i; j<m; ++j) a[i+j*la]=0.0;
		}
		a[i+i*la] = a[i+i*la]+1;
	}

/*	Diagonalisation of the bidiagonal form */
	aeps=eps*rmax;
/*	Loop over the singular values */
	for(k=n-1; k>=0; --k) {
/*	The QR transformation */
		for(itr=1; itr<=itmax; ++itr) {

/*	Test for splitting */
			for(l=k; l>=0; --l) {
				if(fabs(e[l]) < aeps) goto split;
				if(fabs(sigma[l-1]) < aeps) break;
			}

/*	cancellation of E[L] if L>1  */
			c=0.0; s=1.0;
			for(i=l; i<=k; ++i) {
				f=s*e[i];
				e[i] = c*e[i];
				if(fabs(f) < aeps) goto split;
				g=sigma[i];
				sigma[i]=sqrt(f*f+g*g);
				c=g/sigma[i];
				s=-f/sigma[i];

				for(j=0; j<m; ++j) {
					r1= a[j*la+(l-1)];
					r2= a[i+j*la];
					a[j*la+(l-1)]=r1*c+r2*s;
					a[i+j*la]=c*r2-s*r1;
				}
			}

split:			z=sigma[k];
			if(l == k) {
/*	QR iteration has converged */
				if(z < 0.0) {
					sigma[k] = -z;
					for(j=0; j<n; ++j) v[k+j*lv]=-v[k+j*lv];
				}
				break;
			}

			if(itr==itmax) {ier=12; break;}
					
/*	calculating shift from bottom 2x2 minor */
			x=sigma[l];
			y=sigma[k-1];
			g=e[k-1];
			h=e[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.*h*y);
			g=sqrt(1.+f*f);
			if(f < 0.0) g=-g;
			f=((x-z)*(x+z)+h*(y/(f+g)-h))/x;

/*	next QR transformation */
			c=1.0; s=1.0;
/*	Given's rotation  */
			for(i=l+1; i<=k; ++i) {
				g=e[i];
				y=sigma[i];
				h=s*g;
				g=c*g;
				e[i-1]=sqrt(f*f+h*h);
				c=f/e[i-1];
				s=h/e[i-1];
				f=c*x+s*g;
				g=c*g-s*x;
				h=s*y;
				y=c*y;

				for(j=0; j<n; ++j) {
					x=v[j*lv+(i-1)];
					z=v[i+j*lv];
					v[j*lv+(i-1)]=c*x+s*z;
					v[i+j*lv]=c*z-s*x;
				}

				sigma[i-1]=sqrt(f*f+h*h);
				if(sigma[i-1] != 0.0) {
					c=f/sigma[i-1];
					s=h/sigma[i-1];
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for(j=0; j<m; ++j) {
					y= a[j*la+(i-1)];
					z= a[i+j*la];
					a[j*la+(i-1)] = c*y+s*z;
					a[i+j*la] = c*z-s*y;
				}
			}

			e[l]=0.0;
			e[k]=f;
			sigma[k]=x;
		}
	}
	free(e);
	return ier;
}



/*	To evaluate the solution of a system of linear equations using SVD

	N : (input) Number of variables
	M : (input) Number of equations
	U : (input) Array of size LU*M containing the left-hand transformation
	V : (input) Array of size LV*N containing the right-hand transformation
	SIGMA : (input) Array of size N containing the singular values
	LU : (input) Second dimension of array U in the calling function
	LV : (input) Second dimension of array V in the calling function
	B : (input/output) Array of length M containing the RHS
		after execution it will contain the solution
	REPS : (input) Relative accuracy. All singular values < REPS*(Max of singular values)
		will be reduced to zero

	The returned value is always zero

	Required functions : None
*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps)

{
	int i,j;
	double smax, aeps, s;
	double *wk;

/*	Finding the largest singular value */
	smax=0.0;
	for(i=0; i<n; ++i)
		if(sigma[i] > smax) smax=sigma[i];

	aeps=smax*reps;
	wk=(double *)calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i) {
		s=0.0;
/*	Only SIGMA[I] > AEPS contribute to the solution */
		if(sigma[i] > aeps) {
			for(j=0; j<m; ++j) s=s+b[j]*u[i+j*lu];
			s=s/sigma[i];
		}
		wk[i]=s;
	}

	for(i=0; i<n; ++i) {
		s=0.0;
		for(j=0; j<n; ++j) s=s+v[j+i*lv]*wk[j];
		b[i]=s;
	}
	free(wk);
	return 0;
}
