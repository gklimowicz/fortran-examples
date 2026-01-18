/*	Interpolation using Newton's divided difference formula

	XB : (input) value of x at which interpolation is required
	X : (input) Array of length NTAB containing x values
	F : (input) Array of length NTAB containing function values
		F[I] is the tabulated function value at X[I].
	NUSE : (input/output) Number of points to be used for interpolation
		After execution it will contain the number actually used
	NTAB : (input) Number of points in the table
	FB : (output) Array containing interpolated values
	       FB[I] should contain interpolation using I points
	       FB[NUSE] should be the final value
	AEPS : (input) Required accuracy
	DFB : (output) First derivative of interpolating function at XB
	DDFB : (output) Second derivative of interpolating function at XB
		
	Error status is returned by the value of the function DIVDIF.
		0 value implies successful execution
		21 implies NUSE<1, in which case it is set to MIN(6,NTAB)
		22 implies NUSE>NTAB or NMAX, in which case it is reduced
		23 implies interpolation has not converged to specified accuracy

	Required functions : NEARST
*/

#include <math.h>

int nearst(double xb, double x[], int ntab);

int divdif(double xb, double x[], double f[], int *nuse, int ntab,
		double fb[], double aeps, double *dfb, double *ddfb)

{
	int i,j,k,next,in,ip,nit,ier, nmax=10;
	double err,px,dpx,ddpx,xn[11],xd[11];

/*	Find the nearest point */

	next=nearst(xb,x,ntab);
	fb[1]=f[next];
	xd[1]=f[next];
	xn[1]=x[next];
	ier=0;
	px=1.0;

/*	Initialisation for the derivatives */

	*dfb=0.0; *ddfb=0.0;
	dpx=0.0; ddpx=0.0;

/*	Points between IN and IP are used for interpolation */

	ip=next; in=next;

/*	Maximum number of points to be used for interpolation */
	nit=*nuse; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
	if(*nuse>nmax || *nuse>ntab) ier=22;
	if(*nuse<1) {
		ier=21;
		nit=6; if(nmax<nit) nit=nmax; if(ntab<nit) nit=ntab;
	}
	*nuse=1;
		
/*	Calculate successive interpolation polynomial */
	for(j=2; j<=nit; ++j) {

/*	Choose the next nearest point to XB */
		if(in<=0 ) {
			ip=ip+1; next=ip;
		}
		else if(ip >= ntab-1) {
			in=in-1; next=in;
		}
		else if(fabs(xb-x[ip+1]) < fabs(xb-x[in-1]) ) {
			ip=ip+1; next=ip;
		}
		else {
			in=in-1; next=in;
		}

/*	Calculating the divided differences */
		xd[j]=f[next];
		xn[j]=x[next];
		for(k=j-1; k>=1; --k) xd[k]=(xd[k+1]-xd[k])/(xn[j]-xn[k]);

/*	Calculating the derivatives */
		ddpx=ddpx*(xb-xn[j-1])+2.*dpx;
		dpx=dpx*(xb-xn[j-1])+px;
		*dfb = *dfb+dpx*xd[1];
		*ddfb = *ddfb+ddpx*xd[1];

		px=px*(xb-xn[j-1]);
		err=xd[1]*px;
		fb[j]=fb[j-1]+err;
		*nuse=j;

		if(fabs(err) < aeps) return ier;
	}
	return 23;
}
