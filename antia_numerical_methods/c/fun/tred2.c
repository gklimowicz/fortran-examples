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
