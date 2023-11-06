/*	Eigenvalue problem for a complex Hermitian matrix */

#include <stdio.h>
#include <math.h>

int tred2(double *a, int n, int ia, double d[], double e[]);
int tql2(double *z, int n, int iz, double d[], double e[], double reps);
int herevp(double *ar, double *ai, int n, int ia, double ei[], double *vr,
	double *vi, int iv, double reps);

main()
{
	int i,i1,j,n,m,m1,m2, id, iflg, ier,np,nmax;
	double hh, d[100],e[100], vr[4][4],vi[4][4],p,ei[20],eps,reps,p0;
/*	The real part of matrix */
	double ar[4][4]={7,3,1,-1, 3,7,1,-1, 1,1,7,-3, -1,-1,-3,7};
/*	The imaginary part of matrix */
	double ai[4][4]={0,0,2,2, 0,0,-2,-2, -2,2,0,0, -2,2,0,0};

	id=4; n=4; np=0; reps=1.e-14;
	printf(" The matrix is :\n");
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) printf(" %e %e   ",ar[i][j],ai[i][j]);
		printf("\n");
	}

		i=herevp(&ar[0][0],&ai[0][0],n,id,ei,&vr[0][0],&vi[0][0],id,reps);
		printf(" i= %d  n =  %d \n  eigenvalue      eigenvector\n",i,n);
		for(i=0; i<n; ++i) {
			printf(" %e    ",ei[i]);
			for(j=0; j<n; ++j) printf(" %e %e   ",vr[j][i],vi[j][i]);
			printf(" \n");
		}

	return;
}




/*	To find eigenvalues and eigenvectors of a complex Hermitian matrix

	AR : (input) Array of length IA*N containing the real part of the matrix
	AI : (input) Array of length IA*N containing the imaginary part of the matrix
	N : (input) Order of the matrix
	IA : (input) Second dimension of arrays AR, AI as declared in the calling function
	EI : (output) Array of length N containing the eigenvalues which
		should be real.
	VR : (output) Array of length IV*N containing the real part of eigenvectors
		VR[i][j] is the ith component of jth eigenvector
	VI : (output) Array of length IV*N containing the imaginary part of
		eigenvectors VI[i][j] is the ith component of jth eigenvector
	IV : (input) Second dimension of array VR, VI as declared in the calling function
	REPS : (input) Required tolerance, should be of the order of machine
		accuracy
		
	Error status is returned by the value of the function HEREVP.
		0 value implies successful execution
		111 implies that N<=1 or N>IA or N>IZ
			in which case no calculations are done
		Other values may be set by TRED2 or TQL2

	Required functions : TRED2, TQL2
*/

#include <math.h>
#include <stdlib.h>

int tred2(double *a, int n, int ia, double d[], double e[]);
int tql2(double *z, int n, int iz, double d[], double e[], double reps);

int herevp(double *ar, double *ai, int n, int ia, double ei[], double *vr,
	double *vi, int iv, double reps)

{
	int i,j,i1,n2,ier;
	double *wk,*d,*e;

	if(n<=1 || n>ia || n>iv) return 111;

	wk=(double *) calloc((size_t) (4*n*n), sizeof(double));
	d=(double *) calloc((size_t) (2*n), sizeof(double));
	e=(double *) calloc((size_t) (2*n), sizeof(double));
	n2=2*n;

/*	Setup the 2N*2N real matrix */
	for(i=0; i<n; ++i) {
		for(j=0; j<n; ++j) {
			wk[j+i*n2]=ar[j+i*ia];
			wk[j+n+(i+n)*n2]=ar[j+i*ia];
			wk[j+(i+n)*n2]=ai[j+i*ia];
			wk[j+n+i*n2]=-ai[j+i*ia];
		}
	}

/*	To reduce the 2N*2N matrix to tridiagonal form */
	ier=tred2(wk,n2,n2,d,e);
	if(ier>100) {free(e); free(d); free(wk); return ier;}

/*	Find eigenvalues and eigenvectors of tridiagonal matrix */
	ier=tql2(wk,n2,n2,d,e,reps);
	if(ier>100) {free(e); free(d); free(wk); return ier;}
 
/*	Since all eigenvalues are repeated and sorted in ascending order
	pick alternate eigenvalues and eigenvectors */
	for(i=0; i<n; ++i) {
		i1=2*i;
		ei[i]=d[i1];
		for(j=0; j<n; ++j) {
			vr[i+j*iv]=wk[i1+j*n2];
			vi[i+j*iv]=wk[i1+(j+n)*n2];
		}
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



/*	Reduction of real symmetric matrix to tridiagonal form using
	Householder's method

	A : (input/output) Array of length IA*N containing the matrix
		elements. After execution it will be overwritten by the
		transformation matrix
	N : (input) Order of matrix
	IA : (input) The second dimension of A as specified in the calling function
	D : (output) Diagonal elements of the transformed tridiagonal matrix
		D[I] would contain A[I][I].
	E : (output) Off-diagonal elements of the transformed tridiagonal matrix
		E[I+1] would contain A[I][I+1]=A[I+1][I]
		
	Error status is returned by the value of the function TRED2.
		0 value implies successful execution
		107 implies that N<=1 or N>IA, in which case no calculations are done

	Required functions : None
*/

#include <math.h>

int tred2(double *a, int n, int ia, double d[], double e[])

{
	int i,j,k;
	double f,g,h,hh, reps=1.e-300;

	if(n<=1 || n>ia) return 107;

	for(i=n-1; i>=1; --i) {
		f=a[i-1+i*ia];
		g=0.0;
		for(k=0; k<=i-2; ++k) g=g+a[k+i*ia]*a[k+i*ia];
		h=g+f*f;
		if(g<= reps) {
/*	Skip the transformation */
			e[i]=f;
			h=0.0;
		}

		else {
			g=sqrt(h);
			if(f>0.0) g=-g;
			e[i]=g;
			h=h-f*g;
			a[i-1+i*ia]=f-g;
			f=0.0;

			for(j=0; j<=i-1; ++j) {
/*	Elements of u_i/H_i */
				a[i+j*ia]=a[j+i*ia]/h;
				g=0.0;
/*	Form elements of A_iu_i */
				for(k=0; k<=j; ++k) g=g+a[k+j*ia]*a[k+i*ia];
				for(k=j+1; k<=i-1; ++k) g=g+a[j+k*ia]*a[k+i*ia];
/*	Components of p_i */
				e[j]=g/h;
				f=f+g*a[i+j*ia];
			}

/*	calculate u_i^Tp_i/2H_i */
			hh=0.5*f/h;
			for(j=0; j<=i-1; ++j) {
				f=a[j+i*ia];
				g=e[j]-hh*f;
/*	Elements of q_i */
				e[j]=g;
				for(k=0; k<=j; ++k) a[k+j*ia]=a[k+j*ia]-f*e[k]-g*a[k+i*ia];
			}
		}
		d[i]=h;
	}

	d[0]=0.0; e[0]=0.0;
/*	accumulation of transformation matrix Q */
	for(i=0; i<n; ++i) {
		if(d[i] != 0.0) {
			for(j=0; j<=i-1; ++j) {
				g=0.0;
				for(k=0; k<=i-1; ++k) g=g+a[k+i*ia]*a[j+k*ia];
				for(k=0; k<=i-1; ++k) a[j+k*ia]=a[j+k*ia]-g*a[i+k*ia];
			}
		}
		d[i]=a[i+i*ia];
		a[i+i*ia]=1.0;
		for(j=0; j<=i-1; ++j) {a[j+i*ia]=0.0; a[i+j*ia]=0.0;}
	}
	return 0;
}
