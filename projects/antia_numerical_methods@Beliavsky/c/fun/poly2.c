/*	To calculate polynomial interpolation in 2 dimensions on a
	rectangular mesh of points

	(XB1,XB2) : (input) is the point at which interpolation is required
	X1 : (input) Array of length N1 containing the abscissas
	X2 : (input) Array of length N2 containing the abscissas
	F : (input) Array of length NDIM*N2 containing the function values
		F[J][I]=f(X1(I),X2(J))
	NDIM : (input) Second dimension of array F as specified in calling function
	N1 : (input) Length of array X1, i.e Number of points along first dimension
	N2 : (input) Length of array X2, i.e Number of points along second dimension
	NP1 : (input) Number of points to be used for interpolation along X1
	NP2 : (input) Number of points to be used for interpolation along X2
	FB : (output) Interpolated value of function at (XB1,XB2)
		
	Error status is returned by the value of the function POLY2.
		0 value implies successful execution
		206 implies N1 > NDIM, in which case no calculations
			are done

	Required functions : DIVDIF0, NEARST
*/

#include <math.h>

int nearst(double xb, double x[], int ntab);
int divdif0(double xb, double x[], double f[], int *nuse, int ntab,
		double fb[], double aeps, int iflg, int *if1);

int poly2(double xb1, double xb2, double x1[], double x2[], double *f, int ndim,
	int n1, int n2, int np1, int np2, double *fb)

{
	int i,j,k,next,in,ip,nuse,iflg,if1,nit,ier, nmax=10;
	double err,px,reps, xn[11], xd[11], fb1[11];
	
	if(n1>ndim) return 206;
	ier=0;
	reps=0.0;
/*	Find the point nearest to XB1 in X1 */
	next=nearst(xb2,x2,n2);

	nuse=np1; iflg=0;
	i=divdif0(xb1,x1,(f+next*ndim),&nuse,n1,fb1,reps,iflg,&if1);

/*	Set IFLG=1 so that next time DIVDIF0 does not try to locate
	the point again in the table. */
	iflg=1;
	*fb=fb1[nuse];
	xd[1]=*fb;
	xn[1]=x2[next];
	px=1.0;

	ip=next; in=next;

/*	The number of points to be used along X1 */
	nit=np2; if(nmax<nit) nit=nmax; if(n2<nit) nit=n2;
	if(np2<1) {
		nit=4; if(nmax<nit) nit=nmax; if(n2<nit) nit=n2;
	}

/*	Calculate the successive interpolation polynomials */
	for(j=2; j<=nit; ++j) {
	
/*	Find the next nearest point in X1 */
		if(in<=0 ) {
			ip=ip+1; next=ip;
		}
		else if(ip >= n2-1) {
			in=in-1; next=in;
		}
		else if(fabs(xb2-x2[ip+1]) < fabs(xb2-x2[in-1]) ) {
			ip=ip+1; next=ip;
		}
		else {
			in=in-1; next=in;
		}

/*	Interpolate along X2 to calculate function value at (x1[next], xb2) */
		nuse=np1;
		i=divdif0(xb1,x1,(f+next*ndim),&nuse,n1,fb1,reps,iflg,&if1);
		xd[j]=fb1[nuse];
		xn[j]=x2[next];

/*	Calculate the divided differences for interpolation in X1 */
		for(k=j-1; k>=1; --k) xd[k]=(xd[k+1]-xd[k])/(xn[j]-xn[k]);

		px=px*(xb2-xn[j-1]);
		err=xd[1]*px;
		*fb=*fb+err;

	}
	return 0;
}
