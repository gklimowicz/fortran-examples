/*	Linear least square fit using SVD in 1 dimension */

#include <stdio.h>
#include <math.h>

double fun(double x);
void phi(int n, double x[], double f[]);
double rangau(double *seed);
int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);

double seed;

main()
{
	int i,i1,j,nx,nt,m, id, iflg, ier,np,nd;
	double hh, x[100], f[200], ax[90000],ay[10000],reps, chi, a[20][20],y[200],fy[200];
	double ef[200],c[200],cov[100000];

/*	Example 10.3 */

	seed=11; np=100; id=1; reps=1.e-7;


	printf("type nx= no. of pts,  m = No. of basis functions\n");
	scanf(" %d %d", &nx, &m);

	hh=1.0/(nx-1.0);
	nd=1;
	for(i=0; i<nx; ++i) {
		x[i]=i*hh;
	}
	for(i=0; i<nx; ++i) {
		f[i]=fun(x[i]);
		ef[i]=1.e-3;
	}

	j=llsq(nx,m,nd,x,id,f,ef,c,ax,ay,np,np,y,fy,phi,reps,&chi,cov);
	printf(" \n ier = %d   no. of pts = %d    no. of basis functions = %d  , chisq = %e \n", j,nx,m,chi);
	printf("   Coefficients : ");
	for(i=0; i<m; ++i) printf(" %e ",c[i]);
	printf(" \n");

	return;
}

double fun(double x)

{
	double fx;
	fx=((231*x*x-315)*x*x+105)*x*x-5+rangau(&seed)*1.e-3;
	return fx;

}



/*  To calculate the set of basis functions in 1 dimension, these are x**i */

void phi(int n, double *x, double f[])

{
	int i,j;
	double fx,f1;

	for(i=0; i<n; ++i) {
		f1=1.0;
		if(i>0) f1=pow(*x,(double ) i);
		f[i]=f1;
	}
	return;
}



/*	Linear least squares fit in K dimensions

	N : (input) Number of data points to be fitted
	M : (input) Number of basis functions to be used
	K : (input) Number of dimensions
	X : (input) Array of length IX*N containing the coordinates
		of points at which function values is available
		X[j][i] is the ith coordinate of jth point
		The points may have arbitrary distribution in K-space
	IX : (input) Second dimension of X in the calling function, IX >= K 
	F : (input) Array of length N containing the function values
		F[I] should be function value at X[I][0],...,X[I][K-1]
	EF : (input) Array of length N containing the estimated error in F[I]. 
	A : (output) Array of length N containing the fitted coefficients
		Note that although the number of coefficients is M, the
		rest of array is used as scratch space
	U : (output) Array of length IU*N containing the matrix U of SVD
		of the design matrix
	V : (output) Array of length IV*M containing the matrix V of SVD
		of the design matrix
	IU : (input) Second dimension of U in the calling function (IU>=M)
	IV : (input) Second dimension of V in the calling function (IV>=M)
	SIGMA : (output) Array of length M containing the singular values
		of the design matrix
	Y : (output) Array of length N containing the values of fitted
		function at each of the tabular points
	PHI : (input) Name of function to calculate the basis functions
		at any given point
	REPS : (input) Required accuracy for solution of equations using SVD
		singular values less than REPS times maximum will be set to zero
	CHISQ : (output) The value of Chi square at minimum
        COV : (output) Array of length IV*M containing the covariance
                matrix of the fitted parameters. COV(I,I) will be the
                variance in A(I).
		
	Error status is returned by the value of the function LLSQ.
		0 value implies successful execution
		606 implies that M>N, M<=0, N<=0, or K>IX
		607 implies that EF[I]<=0 for some I
		No calculations are done in both these cases
		Other values may be set by SVD

	Function PHI(M,X,Y) must be supplied by the user to calculate
	the required basis functions. M is the number of basis functions,
	X is an array of length K containing the coordinates of point
	at which the basis function needs to be calculated. Y is an
	array of length M containing the computed basis functions at X

	Required functions : SVD, SVDEVL, PHI

	THE ARGUMENTS OF THIS FUNCTION HAVE CHANGED FROM THE EARLIER VERSION.
	NOW THERE IS AN ADDITIONAL ARGUMENT COV TO CALCULATE THE COVARIANCE
	MATRIX.

*/

#include <math.h>
#include <stdlib.h>

int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
	int lv, double b[], double reps);
int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);

int llsq(int n, int m, int k, double *x, int ix, double f[], double ef[],
	double a[], double *u, double *v, int iu, int iv, double sigma[],
	double y[], void (* phi) (int , double * , double * ), double reps,
	double *chisq, double *cov)

{
	int i,j,ik,ier;
	double s1,sigmax;
	double *wk;

	if(m>n || m<=0 || n<=0 || k>ix) return 606;
	
	wk=(double *) calloc((size_t) m, sizeof(double));
/*	Setting up the design matrix and the RHS */
	for(i=0; i<n; ++i) {
		if(ef[i]<=0.0) {free(wk); return 607;}
		a[i]=f[i]/ef[i];
		phi(m,&x[i*ix],wk);
		for(j=0; j<m; ++j) u[j+i*iu]=wk[j]/ef[i];
	}

	ier=svd(m,n,u,v,sigma,iu,iv);
	if(ier>100) {free(wk); return ier;}

/*	Calculate the least squares solution */
	ier=svdevl(m,n,u,v,sigma,iu,iv,a,reps);
 
/*	Computing the \chi^2 from fitted coefficients */
	*chisq=0.0;
	for(i=0; i<n; ++i) {
		phi(m,&x[i*ix],wk);
		s1=0.0;
		for(j=0; j<m; ++j) s1=s1+a[j]*wk[j];
		y[i]=s1;
		s1=(f[i]-y[i])/ef[i];
		*chisq=(*chisq)+s1*s1;
	}
 
/*	Computing the covariance matrix for fitted coefficients */
	sigmax=0.0;
	for(i=0; i<m; ++i)
		if(sigma[i]>sigmax) sigmax=sigma[i];
	for(i=0; i<m; ++i) {
		for(j=0; j<=i; ++j) {
			cov[j+i*iv]=0.0;
			for(ik=0; ik<m; ++ik) if(sigma[ik]>reps*sigmax)
			cov[j+i*iv]=cov[j+i*iv]+v[ik+j*iv]*v[ik+i*iv]/(sigma[ik]*sigma[ik]);
			cov[i+j*iv]=cov[j+i*iv];
		}
	}


	free(wk);
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
	wk=(double *) calloc((size_t) n, sizeof(double));
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



/*	To generate random numbers with Gaussian probability distribution
	It generates random numbers with zero mean and variance of 1.
	
	SEED : (input/output) real seed, it should be positive and
		less than AN. It is updated by the function and should
		not be modified between two calls, unless a fresh
		sequence is required

	Required functions : None
	
	THE ARGUMENT OF THIS FUNCTION HAS CHANGED AS COMPARED TO EARLIER
	VERSION AS THE SEED IS NOW DOUBLE INSTEAD OF INT.

*/

#include <math.h>


double rangau(double *seed)

{
	int n2;
	double am=2147483648.0, a=45875.0, ac=453816693.0, an=2147483647.0, r1, rn;

	rn=a*(*seed)+ac; n2=rn/am; r1=rn-am*n2;
	if(*seed==0.0) *seed=0.1;
	rn=sqrt(2.0*log(an/(*seed)))*cos(2.0*M_PI*r1/an);
	*seed=r1;
	return rn;
}

