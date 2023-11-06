/*	To solve a system of linear equations arising from finite difference
	approximation of ordinary differential equations

	N : (input) Number of mesh points used.
	M : (input) Number of first order differential equations in the system
	ML : (input) Number of boundary conditions at the first boundary t=T[0]
	A : (input/output) Array of length (M+ML)*2M*N containing
		the matrix of equations. After execution it will contain
		the triangular decomposition
	IFLG : (input/output) Integer variable used as a flag to decide
		the nature of computation.
		If IFLG=0 the triangular decomposition, determinant and
			solution of equations is computed. IFLG is set to 2
		If IFLG=1 only the triangular decomposition and determinant
			are calculated and IFLG is set to 2
		If IFLG=2 then it is assumed that triangular decomposition
			is already done and available in A and only the solution
			of system of equations is solved
	DET : (output) Scaled value of determinant of the matrix
	IDET : (output) exponent of determinant, the value of determinant
		is DET*2**IDET
	INC : (input/output) Integer array of length M*N containing the
		information about interchanges used by Gaussian elimination.
		It is calculated if IFLG=0 or 1, while for IFLG=2 it must
		be supplied from previous calculation.
	X : (input/output) Array of length M*N containing the right hand
		side of equation at input. After execution it will be overwritten
		by the solution if IFLG=0 or 2. For IFLG=1 it is not used.
		
	Error status is returned by the value of the function GAUBLK.
		0 value implies successful execution
		735 implies that the matrix is singular

	Required functions : None
*/

#include <math.h>

int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
	int *idet, int *inc, double *x)

{
	int i,j,k,mr,mc,kmax,ki,kj,kk,m1,m2;
	double rmax,r1,at,xt,d1;

	m1=m+ml;
	m2=m1*2*m;

	if(*iflg<=1) {
		*idet=0; *det=1.0;
		mr=m1;		/*	The number of rows in each block of the matrix */
		mc=2*m;		/*	The number of columns in each block of the matrix*/

		for(j=0; j<n; ++j) {
			if(j==n-1) {mr=m; mc=m;}

			kk=mr-1; if(m<kk) kk=m;
			for(k=0; k<kk; ++k) {
				rmax=fabs(a[k+k*m1+j*m2]);
				kmax=k;
/*		Find the pivot */
				for(ki=k+1; ki<mr; ++ki) {
					r1=fabs(a[ki+k*m1+j*m2]);
					if(r1>rmax) {rmax=r1; kmax=ki;}
				}
				inc[k+j*m]=kmax;

				if(kmax != k) {
/*	exchange rows K and KMAX */
					*det= -(*det);
					for(ki=k; ki<mc; ++ki) {
						at=a[k+ki*m1+j*m2];
						a[k+ki*m1+j*m2]=a[kmax+ki*m1+j*m2];
						a[kmax+ki*m1+j*m2]=at;
					}
				}

				*det=(*det)*a[k+k*m1+j*m2];
/*		If the pivot is zero, then quit */
				if(a[k+k*m1+j*m2]==0.0) return 735;

/*	Gaussian elimination */
				for(ki=k+1; ki<mr; ++ki) {
					a[ki+k*m1+j*m2]=a[ki+k*m1+j*m2]/a[k+k*m1+j*m2];
					for(kj=k+1; kj<mc; ++kj) a[ki+kj*m1+j*m2] =
						a[ki+kj*m1+j*m2]-a[ki+k*m1+j*m2]*a[k+kj*m1+j*m2];
				}
			}

			if(*det != 0.0) {
/*	Scale the determinant if necessary */
				while(fabs(*det)>32.0) {
					*det=(*det)*0.03125;
					*idet=(*idet)+5;
				}

				while(fabs(*det)<0.03125) {
					*det=(*det)*32.0;
					*idet=(*idet)-5;
				}
			}

/*	Copy the overlapping elements into the next block */
			if(j<n-1) {
				for(k=0; k<ml; ++k) {
					for(ki=0; ki<m; ++ki) a[k+ki*m1+(j+1)*m2]=a[k+m+(ki+m)*m1+j*m2];
				}
			}
		}

		inc[n*m-1]=m-1;
		*det=(*det)*a[m-1+(m-1)*m1+(n-1)*m2];
		if(a[m-1+(m-1)*m1+(n-1)*m2]==0.0) return 735;

		if(*iflg==1) {*iflg=2; return 0;}
		*iflg=2;
	}

/*	Solve the system of linear equations */
	mr=m1;
	for(j=0; j<n; ++j) {
		if(j==n-1) mr=m;
		kj=mr-1; if(m<kj) kj=m;
		for(k=0; k<kj; ++k) {
			kk=inc[k+j*m];
			if(k != kk) {
/*	exchange the corresponding elements of RHS */
				xt=x[k+j*m];
				x[k+j*m]=x[kk+j*m];
				x[kk+j*m]=xt;
			}
			
/*	Gaussian elimination */
			for(ki=k+1; ki<mr; ++ki) x[ki+j*m]=x[ki+j*m]-a[ki+k*m1+j*m2]*x[k+j*m];
		}
	}

/*	back-substitution */
	mc=m;
	for(j=n-1; j>=0; --j) {
		for(k=m-1; k>=0; --k) {
			d1=x[k+j*m];
			for(ki=mc-1; ki>=k+1; --ki) d1=d1-x[ki+j*m]*a[k+ki*m1+j*m2];
			x[k+j*m]=d1/a[k+k*m1+j*m2];
		}
		mc=2*m;
	}
	return 0;
}
