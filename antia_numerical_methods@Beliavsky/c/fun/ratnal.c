/*	To calculate rational function interpolation

	XB : (input) x value at which interpolation is required
	X : (input) Array of length NTAB containing the abscissas
	F : (input) Array of length NTAB containing the function values at X[I]
	NUSE : (input/output) Number of points to be used for interpolation
		after execution NUSE will contain the number of points actually used
	NTAB : (input) Number of points in the table
	FB : (output) The interpolated value at x=XB
	AEPS : (input) The required accuracy in interpolation
		
	Error status is returned by the value of the function RATNAL.
		0 value implies successful execution
		21 implies NUSE <2 in which case it is increased
		22 implies NUSE>NTAB or NMAX, in which case it is reduced
		23 implies interpolation did not converge to specified accuracy
		205 implies denominator turns out to be zero and
			interpolation cannot proceed

	Required functions : NEARST
*/

#include <math.h>

int nearst(double xb, double x[], int ntab);

int ratnal(double xb, double x[], double f[], int *nuse, int ntab,
	double *fb, double aeps)

{
	int i,j,ii,next,in,ip,nit,ier, nmax=10;
	double w,rx,fac,xn[10],d[10],c[10];

/*	To find the entry nearest to xb */
	next=nearst(xb,x,ntab);
	*fb=f[next];
	d[0]=f[next];
	c[0]=d[0];
	xn[0]=x[next];
	ier=0;

	ip=next; in=next;

/*	Maximum number of points to be used for interpolation */
	nit=*nuse; if(ntab<nit) nit=ntab; if(nmax<nit) nit=nmax;
	if(*nuse>nmax || *nuse>ntab) ier=22;
	if(*nuse<1) {
		ier=21;
		nit=6; if(ntab<nit) nit=ntab; if(nmax<nit) nit=nmax;
	}
	*nuse=1;

/*	If XB coincides with a tabular point, then return */
	if(xb == xn[0]) return 0;

/*	Calculate the successive rational function interpolation */
	for(j=1; j<nit; ++j) {
/*	Choosing the next nearest point */
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

		xn[j]=x[next];
		c[j]=f[next];
		d[j]=f[next];

/*	Using the recurrences to calculate the differences C[I] and D[I] */
		for(ii=j-1; ii>=0; --ii) {
			w=c[ii+1]-d[ii];
			rx=(xb-xn[ii])/(xb-xn[j]);
			fac=rx*d[ii]-c[ii+1];

			if(fac == 0) return 205;

			fac=w/fac;
			c[ii]=rx*d[ii]*fac;
			d[ii]=c[ii+1]*fac;
		}

		*fb = *fb + c[0];
		*nuse=j+1;
		if(fabs(c[0])<aeps) return ier;
	}
	return 23;
}
