/*	To search for complex zeros in a rectangular region in complex plane

	RX1 : (input) Lower limit of real part for the region to be searched
	RX2 : (input) Upper limit of real part for the region to be searched
	RY1 : (input) Lower limit of imaginary part for the region to be searched
	RY2 : (input) Upper limit of imaginary part for the region to be searched
		The region with real part between RX1 & RX2 and 
		imaginary part between RY1 & RY2 will be searched
	NX : (input/output) number of points to be used along real axis
	NY : (input/output) number of points to be used along imaginary axis
	CFUN : (input) Name of the function routine to calculate the complex
		function. Function CFUN(CZ,CF) must be supplied by the user.
		CZ and CF are both complex variables (array of length 2 containing
		the real and imaginary parts). CF should be the value of complex
		function at CZ.

	The zero will be located near the point where all 4 quadrants
	in function values meet in the figure that is produced.

	The returned value is always zero.

	Required functions : CFUN
*/

#include <stdio.h>
#include <math.h>

int search(double rx1, double rx2, double ry1, double ry2, int *nx,
	int *ny, void cfun(double * , double * ))

{
	int i,j, iq[41], imax=41;
	double hx,hy,y,rx,ry, x[41],cz[2],cf[2];

	if(rx1==rx2 || ry1==ry2) return 0;
	if(*nx > imax) *nx=imax;
	if(*nx<=1) *nx=21;
	if(*ny<=1) *ny=21;

	hx=(rx2-rx1)/(*nx-1);
	hy=(ry2-ry1)/(*ny-1);
	for(i=0; i<(*nx); ++i) x[i]=rx1+hx*i;

	for(i=0; i<(*ny); ++i) {
		y=ry2-hy*i;
		for(j=0; j<(*nx); ++j) {
			cz[0]=x[j]; cz[1]=y; cfun(cz,cf);

/*	Determine the quadrant in which the function value lies */
			rx=cf[0]; ry=cf[1];
			if(rx>=0.0 && ry>= 0.0) iq[j]=1;
			else if(rx<0.0 && ry>= 0.0) iq[j]=2;
			else if(rx<0.0 && ry< 0.0) iq[j]=3;
			else iq[j]=4;
		}
		printf("%12.3e  ",y);
		for(j=0; j<(*nx); ++j) printf("%3d",iq[j]);
		printf("\n");
	}

	printf("\n   ");
	for(i=0; i<(*nx); i += 5) printf("%15.3e",x[i]);
	printf("\n");
	return 0;
}
