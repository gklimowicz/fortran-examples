/*	To calculate spherical harmonic */

#include <stdio.h>
#include <math.h>

void plm(int l, int m, double x, double p[]);
void ylm(int l, int m, double theta, double phi, double y[]);
void ylm_x(int l, int m, double x, double phi, double y[]);

main()
{
	int i,i1,l,m,n;
	double x,theta,phi,f1,f2,y1[2],y[2];

	for(i1=0; i1<99; ++i1) {
		printf("type l, m, theta, phi       (quits when theta=-9999)\n");
		scanf(" %d %d %le %le",&l,&m,&theta,&phi);
		if(theta==-9999) return 0;

		ylm(l,m,theta,phi,y1);
		x=cos(theta);
		ylm_x(l,m,x,phi,y);
		printf(" l = %d    m = %d    theta = %e   phi = %e \n",l,m,theta,phi);
		printf(" x = %e   Ylm(x) = %e %e    %e %e \n",x,y[0],y[1],y1[0],y1[1]);
	}
	return;
}




/*	To calculate the associated Legendre functions P_lm(X)

	L,M : (input) Order of function (L>=0), abs(M)<=L
	X : (input) Argument at which the value of polynomial is required
	P : (output) Array of length L+1, which will contain the
		calculated values of polynomials. P[j] will contain
		P_jM(X) for j>=M

	Required functions : None
*/

#include <math.h>

void plm(int l, int m, double x, double p[])

{
	int i,n,mm;
	double pm,rm;

	if(l<0 || abs(m)>l ) return;
 
/*	First compute P_MM(x) */
	mm=abs(m);
	pm=1.0;
	for(i=1; i<2*mm; i += 2) pm=pm*i;

	if(m<0) {
/*	Modify the normalisation factor */
		for(i=1; i<=2*mm; ++i) pm=pm/i;
		if(mm - 2*(mm/2) == 1) pm=-pm;
	}

	rm=mm/2.0;
	if(mm == 0) p[mm]=1.0;
	else if(fabs(x)<1.0) p[mm]=pm*pow(1.0-x*x, rm);
	else p[mm]=0.0;

/*	Use the recurrence relation to compute P_nM(x) */
	p[mm+1]=x*p[mm]*(2*mm+1)/(mm+1-m);
	for(n=mm+2; n<=l; ++n) p[n]=((2*n-1)*x*p[n-1]-(n-1+m)*p[n-2])/(n-m);

	return;
}




/*	To compute spherical harmonic Y_lm(THETA,PHI)

	L : (input) Degree of spherical harmonic
	M : (input) Azimuthal order of spherical harmonic
	THETA, PHI : (input) Real variables specifying the angular coordinates 
		at which the spherical harmonic needs to be evaluated
	Y : (output) Array of length 2 containing the complex value of
		spherical harmonic.

	Required functions : PLM
*/

#include <math.h>
#include <stdlib.h>

#define PI 3.14159265358979324

void plm(int l, int m, double x, double p[]);

void ylm(int l, int m, double theta, double phi, double y[])

{
	int i,mm;
	double x,clm;
	double *p;

	y[0]=0.0; y[1]=0.0;
	if(l<0 || abs(m)>l) return;

	p=(double *) calloc((size_t) (l+2), sizeof(double));
	mm=abs(m);
 
/*	To use X instead of THETA in argument comment out this line */
	x=cos(theta);
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
