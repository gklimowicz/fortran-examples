/*	To calculate the discrete Fourier Transform
	Normal evaluation of the sum, valid for any number of points
	Should be used when the number of points is small

	N : (input) Number of points
	CG : (input) Array of length 2*N containing the data points
	CF : (output) Array of length 2*N which will contain the
		Fourier transform of CG
		Both CG and CF are complex numbers and hence the arrays should
		contain the real and imaginary parts.
	IFLG : (input) Flag to decide whether to calculate forward or inverse
		transform. If IFLG>=0 then Fourier transform is calculated
		IF IFLG<0 then inverse Fourier transform is calculated
		
	Error status is returned by the value of the function DFT.
		0 value implies successful execution
		611 implies that N<2 and no calculations are done

	Required functions : None
*/

#include <math.h>

#define PI 3.14159265358979324

int dft(int n, double cg[][2], double cf[][2], int iflg)

{
	int i,iw,j;
	double cs[2],cw;

	if(n<2) return 611;

	if(iflg>=0) iw=1;
	else iw=-1;

/*	Loop for Fourier transform */
	for(i=0; i<n; ++i) {
		cw=2.0*iw*i*PI/n;
		cs[0]=0.0; cs[1]=0.0;
		for(j=0; j<n; ++j) {
			cs[0]=cs[0]+cg[j][0]*cos(j*cw)-cg[j][1]*sin(j*cw);
			cs[1]=cs[1]+cg[j][1]*cos(j*cw)+cg[j][0]*sin(j*cw);
		}
		cf[i][0]=cs[0]; cf[i][1]=cs[1];
	}
	return 0;
}
