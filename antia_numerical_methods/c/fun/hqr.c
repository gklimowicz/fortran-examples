/*	To find all eigenvalues of an upper-Hessenberg matrix using QR algorithm

	H : (input/output) Array of length IH*NN containing the matrix
		elements. The contents of array are destroyed during execution
	NN : (input) Order of matrix
	IH : (input) Second dimension of array H as declared in the calling
		function
	ER : (output) Array of length NN containing the real part of
		the eigenvalues.
	EI : (output) Array of length NN containing the imaginary part of
		the eigenvalues.
	REPS : (input) Required relative accuracy in eigenvalues. It should
		be of the order of machine accuracy
		
	Error status is returned by the value of the function HQR.
		0 value implies successful execution
		114 implies that N<=1 or N>IH, in which case no
			calculations are done
		145 implies that QR iteration failed to converge at
			some stage and calculations are abandoned.

	Required functions : None
*/

#include <math.h>

int hqr(double *h, int nn, int ih, double er[], double ei[], double reps)

{
	int i,j,k,l,m,n,it, nit=30;
	double t,x,w,y,z,p,d,r,s,u;

	t=0.0; n=nn-1;
	if(nn<=1 || nn>ih) return 114;

/*	If all eigenvalues are found, then quit */
nexti: if(n<0) return 0;

/*	Loop for double QR iteration */
	for(it=1; it<=nit; ++it) {
/*	Search for a small subdiagonal element */
		for(l=n; l>=1; --l) {
			if(fabs(h[l-1+l*ih]) <= reps*(fabs(h[l-1+(l-1)*ih])+fabs(h[l+l*ih]))) break;
		}
		x=h[n+n*ih];

		if(l==n) {
/*	One eigenvalue is isolated */
			er[n]=x+t; ei[n]=0.0;
			n=n-1;
			goto nexti;
		}

		y=h[n-1+(n-1)*ih];
		w=h[n-1+n*ih]*h[n+(n-1)*ih];
		if(l == n-1) {
/*	A pair of eigenvalues are isolated */
			p=(y-x)/2.0;
			d=p*p+w;
			y=sqrt(fabs(d));
			x=x+t;
			if(d>0.0) {
/*	Pair of real eigenvalues */
				if(p<0.0) y=-y;
				y=p+y;
				er[n-1]=x+y; er[n]=x-w/y;
				ei[n-1]=0.0; ei[n]=0.0;
			}
			else {
/*	Pair of complex conjugate eigenvalues */
				er[n]=x+p; er[n-1]=er[n];
				ei[n-1]=y; ei[n]=-y;
			}
			n=n-2;
			goto nexti;
		}

		if(it==10 || it==20) {
/*	Apply special shifts */
			t=t+x;
			for(i=0; i<=n; ++i) h[i+i*ih]=h[i+i*ih]-x;
			s=fabs(h[n-1+n*ih])+fabs(h[n-2+(n-1)*ih]);
			x=0.75*s;
			y=x;
			w=-0.4375*s*s;
		}

/*	Search for two consecutive small subdiagonal elements */
		for(m=n-2; m>=l; --m) {
			z=h[m+m*ih];
			r=x-z; s=y-z;
			p=(r*s-w)/h[m+(m+1)*ih]+h[m+1+m*ih];
			u=h[m+1+(m+1)*ih]-z-r-s;
			r=h[m+1+(m+2)*ih];
			s=fabs(p)+fabs(u)+fabs(r);
			p=p/s; u=u/s; r=r/s;
			if(m==l) break;
			if(fabs(h[m-1+m*ih])*(fabs(u)+fabs(r))<= reps*fabs(p)*
				(fabs(h[m-1+(m-1)*ih])+fabs(z)+fabs(h[m+1+(m+1)*ih]))) break;
		}
		for(i=m+3; i<=n; ++i) {h[i-3+i*ih]=0.0; h[i-2+i*ih]=0.0;}
		h[m+(m+2)*ih]=0.0;

/*	Double QR transform for rows L to N and columns M to N */
		for(k=m; k<=n-1; ++k) {
			if(k != m) {
				p=h[k-1+k*ih];
				u=h[k-1+(k+1)*ih];
				r=0.0;
				if(k != n-1) r=h[k-1+(k+2)*ih];
				x=fabs(p)+fabs(u)+fabs(r);
				if(x != 0.0) {p=p/x; u=u/x; r=r/x;}
			}

			if(x != 0.0) {
				s=sqrt(p*p+u*u+r*r);
				if(p<0.0) s=-s;
				if(k != m) h[k-1+k*ih]=-s*x;
				else if(l != m) h[k-1+k*ih]=-h[k-1+k*ih];
				p=p+s; x=p/s; y=u/s; z=r/s; u=u/p; r=r/p;

/*	Row modification */
				for(j=k; j<=n; ++j) {
					p=h[j+k*ih]+u*h[j+(k+1)*ih];
					if(k != n-1) {
						p=p+r*h[j+(k+2)*ih];
						h[j+(k+2)*ih]=h[j+(k+2)*ih]-p*z;
					}
					h[j+(k+1)*ih]=h[j+(k+1)*ih]-p*y;
					h[j+k*ih]=h[j+k*ih]-p*x;
				}

/*	Column modification */
				j=k+3; if(n<j) j=n;
				for(i=l; i<=j; ++i) {
					p=x*h[k+i*ih]+y*h[k+1+i*ih];
					if(k != n-1) {
						p=p+z*h[k+2+i*ih];
						h[k+2+i*ih]=h[k+2+i*ih]-p*r;
					}
					h[k+1+i*ih]=h[k+1+i*ih]-p*u;
					h[k+i*ih]=h[k+i*ih]-p;
				}
			}
		}
	}

	return 145;
}
