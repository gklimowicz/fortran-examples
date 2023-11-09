/*	to calculate weights and abscissas of Gauss-Laguerre quadrature formulas */

#include <stdio.h>
#include <math.h>

int tql2(double *z, int n, int iz, double d[], double e[], double reps);
int gausrc(int n, double w[], double ab[], double *cof, double ri0);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double xl, xu,  ri0, a[50], w[50], cof[50][3];

/*	For Laguerre polynomial */
	ri0=1.0;

	for(i1=0; i1<99; ++i1) {
		printf("type  nuse = No. of points in quadrature formula \n");
		printf("                        (quits when n<=0)\n");
		scanf(" %d", &nuse);
		if(nuse<=0) return 0;

		for(i=1; i<=nuse; ++i) {
/*	For Laguerre polynomial */
			cof[i-1][0]=1.0;
			cof[i-1][1]=-(2*i-1);
			cof[i-1][2]=(i-1.0)*(i-1.0);
		}
		i=gausrc(nuse,w,a,&cof[0][0],ri0);
		printf(" ier = %d  n = %d  \n", i,nuse);
		printf("  Abscissas             Weights  \n");
		for(i=0; i<nuse; ++i) printf(" %20.13e %20.13e \n",a[i],w[i]);
	}
	return;
}


/*	To calculate weights and abscissas of a quadrature formula with
	specified weight function when recurrence relation for orthogonal
	polynomial is known.

	N : (input) Number of points in the required quadrature formula
	W : (output) Array of length N, which will contain the weights
	AB : (output) Array of length N containing the abscissas
	COF : (input) Array of length 3*N containing the coefficients of
		the recurrence relation for orthogonal polynomials
		P_i(x)=(COF[3*i]*x+COF[3*i+1])P_{i-1}(x) - COF[3*i+2]*P_{i-2}(x)
	RI0 : (input) The integral of weight function over the required interval.
		
	Error status is returned by the value of the function GAUSRC.
		0 value implies successful execution
		302 implies N<=0
		321 implies that some coefficient becomes imaginary
			during calculations.
		In both these cases calculations are abandoned.
		Other values may be set by TQL2

     Required functions : TQL2
*/
 
#include <math.h>
#include <stdlib.h>

int tql2(double *z, int n, int iz, double d[], double e[], double reps);

int gausrc(int n, double w[], double ab[], double *cof, double ri0)

{
	int i,j,ier;
	double r1,reps=1.0e-15;
	double *wk,*d,*e;

	if(n<=0) return 302;
	wk=(double *) calloc((size_t) (n*n), sizeof(double));
	d=(double *) calloc((size_t) n, sizeof(double));
	e=(double *) calloc((size_t) n, sizeof(double));

/*     Calculate the coefficients of symmetric tridiagonal matrix */
	for(i=0; i<n; ++i) {
		d[i]=-cof[3*i+1]/cof[3*i];
		if(i<n-1) {
			r1=cof[3*i+5]/(cof[3*i]*cof[3*i+3]);
			if(r1>=0.0) e[i+1]=sqrt(r1);
			else {free(e); free(d); free(wk); return 321;}
		}
		for(j=0; j<n; ++j) wk[j+i*n]=0.0;
		wk[i+i*n]=1.0;
	}

/*     Find eigenvalues and eigenvectors of the tridiagonal matrix */
	ier=tql2(wk,n,n,d,e,reps);
	if(ier>0) {free(e); free(d); free(wk); return ier;}

/*     Calculate the abscissas and weights */
	for(i=0; i<n; ++i) {
		ab[i]=d[i];
		w[i]=wk[i]*wk[i]*ri0;
	}
	free(e); free(d); free(wk);
	return 0;
}


/*	To find eigenvalues and eigenvectors of Z T Z^T using QL algorithm
	where T is a symmetric tridiagonal matrix and Z is an orthogonal matrix.
	If Z is the transformation matrix to reduce original matrix to
	tridiagonal matrix, it will calculate the eigenvectors of original matrix

	Z : (input/output) Array of length IZ*N which should contain
		the transformation matrix required to reduce original real
		symmetric matrix to tridiagonal form. To find eigenvectors
		of a symmetric tridiagonal matrix, set Z to unit matrix.
		After execution Z will contain the eigenvector of the original
		matrix Z T Z^T. Z[i][j] should contain the ith component of
		jth eigenvector
	N : (input) Order of matrix
	IZ : (input) The second dimension of array Z as declared in the
		calling function. (IZ>=N)
	D : (input/output) Array of length N, containing the diagonal
		elements of the tridiagonal matrix, D[i]=T[i][i]
		After execution it will contain the eigenvalues of the matrix
	E : (input/output) Array of length N containing the off-diagonal
		elements of the tridiagonal matrix, E[i]=T[i][i+1]=T[i+1][i]
		It is used as scratch space and its contents will be destroyed
		during execution.
	REPS : (input) Required tolerance, it should be of order of machine
		accuracy
		
	Error status is returned by the value of the function TQL2.
		0 value implies successful execution
		108 implies that N<=1 or N>IZ, in which case no
			calculations are performed
		143 implies that the QL algorithm failed to converge
			for some eigenvalue, the calculations are abandoned

	Required functions : None
*/

#include <math.h>

int tql2(double *z, int n, int iz, double d[], double e[], double reps)

{
	int i,j,k,l,m,it, nit=30;
	double b,f,g,h,p,r,c,s;

	if(n<=0 || n>iz) return 108;

	for(i=1; i<n; ++i) e[i-1]=e[i];
	e[n-1]=0.0;
	b=0.0; f=0.0;

	for(l=0; l<n; ++l) {
		h=reps*(fabs(d[l])+fabs(e[l]));
		if(b<h) b=h;
/*	Look for small off-diagonal elements */
		for(i=l; i<n; ++i) {
			m=i;
			if(fabs(e[m])<=b) break;
		}
		if(m==l) {
/*	one eigenvalue is isolated */
			d[l]=d[l]+f;
			goto nextit;
		}

/*	Loop for QL transformation  */
		for(it=1; it<=nit; ++it) {
/*	Find shift */
			g=d[l];
			p=(d[l+1]-g)/(2.0*e[l]);
			r=sqrt(p*p+1.0);
			if(p<0.0) r=-r;
			d[l]=e[l]/(p+r);
			h=g-d[l];
			for(i=l+1; i<n; ++i) d[i]=d[i]-h;
			f=f+h;

/*	The QL transformation */
			p=d[m];
			c=1.0; s=0.0;
/*	Given's rotations */
			for(i=m-1; i>=l; --i) {
				g=c*e[i];
				h=c*p;
				if(fabs(p)>=fabs(e[i])) {
					c=e[i]/p;
					r=sqrt(c*c+1.0);
					e[i+1]=s*p*r;
					s=c/r;
					c=1.0/r;
				}
				else {
					c=p/e[i];
					r=sqrt(c*c+1.0);
					e[i+1]=s*e[i]*r;
					s=1.0/r;
					c=c/r;
				}
				p=c*d[i]-s*g;
				d[i+1]=h+s*(c*g+s*d[i]);

/*	Transforming the eigenvectors */
				for(k=0; k<n; ++k) {
					h=z[i+1+k*iz];
					z[i+1+k*iz]=s*z[i+k*iz]+c*h;
					z[i+k*iz]=c*z[i+k*iz]-s*h;
				}
			}
			e[l]=s*p;
			d[l]=c*p;
			if(fabs(e[l])<=b) {
/*	One eigenvalue is isolated */
				d[l]=d[l]+f;
				goto nextit;
			}
		}
/*	QL iteration fails to converge */
		return 143;

nextit: b=b;
	}

/*	Sort eigenvalues in ascending order by straight selection */
	for(i=0; i<n-1; ++i) {
		k=i;
		p=d[i];
		for(j=i+1; j<n; ++j) {
			if(d[j]<p) {k=j; p=d[j];}
		}
		if(k != i) {
/*	exchange the eigenvalues and eigenvectors */
			d[k]=d[i];
			d[i]=p;
			for(j=0; j<n; ++j) {
				p=z[i+j*iz]; z[i+j*iz]=z[k+j*iz]; z[k+j*iz]=p;
			}
		}
	}
	return 0;
}
