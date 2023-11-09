/*	To compute spherical harmonic Y_lm(COS(THETA),PHI)
	Version of YLM with argument as COS(THETA) instead of THETA

	L : (input) Degree of spherical harmonic
	M : (input) Azimuthal order of spherical harmonic
	X, PHI : (input) Real variables specifying the values of cos(theta)
		and phi at which the spherical harmonic needs to be evaluated
	Y : (output) Array of length 2 containing the complex value of
		spherical harmonic.

	Required functions : PLM
*/

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979324

void plm(int l, int m, double x, double p[]);

void ylm_x(int l, int m, double x, double phi, double y[])

{
	int i,mm;
	double clm;
	double *p;

	y[0]=0.0; y[1]=0.0;
	if(l<0 || abs(m)>l) return;

	p=(double *) calloc((size_t) (l+2), sizeof(double));
	mm=abs(m);
 
/*	To use THETA instead of THETA in argument uncomment  this line 
	x=cos(theta); */
	plm(l,mm,x,p);

	clm=(2*l+1.0)/(4.*PI);
	for(i=l-mm+1; i<=l+mm; ++i) clm=clm/i;
	clm=sqrt(clm);

	if(mm - 2*(mm/2) == 1 && m>=0) clm=-clm;
	y[0]=clm*p[l]*cos(m*phi);
	y[1]=clm*p[l]*sin(m*phi);
	free(p);
	return;
}
