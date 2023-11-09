/*	To calculate the discrete Fourier Transform of real data using FFT algorithm

	N : (input) Number of points, which must be a power of 2
	CG : (input/output) Array of length 2*(N/2) containing the data points
		After execution it will contain the Fourier transform of CG
		In the calling function this array may be treated as a
		real array of length N, though the Fourier transform
		will be complex.
	IFLG : (input) Flag to decide whether to calculate forward or inverse
		transform. If IFLG>=0 then Fourier transform is calculated
		IF IFLG<0 then inverse Fourier transform is calculated
		
	Error status is returned by the value of the function FFTR.
		0 value implies successful execution
		611 implies that N<4, no calculations are done
		631 implies that N is not a power of 2, in this case
			contents of CG will be destroyed but will not
			contain the Fourier transform.

	Required functions : FFT, CDIV
*/

#include <math.h>

#define PI 3.14159265358979324

void cdiv(double c1[2], double c2[2], double cr[2]);

int fftr(int n, double cg[][2], int iflg)

{
	int i,i1,nn,ier;
	double r,gr,gi,th,cw[2],cf[2],c1[2],c2[2],cwf[2],ct[2];

	nn=n/2; th=PI/nn;
	if(iflg>=0) {
		cw[0]=cos(th); cw[1]=sin(th);
		cf[0]=0.0; cf[1]=-0.5;
/*	Calculate the DFT of complex array of length N/2 */
		ier=fft(nn,cg,iflg);
		if(ier>0) return ier;
	}
	else {
		cf[0]=0.0; cf[1]=0.5;
		cw[0]=cos(th); cw[1]=-sin(th);
	}

/*	Rearranging the DFT */
	cwf[0]=cw[0]; cwf[1]=cw[1];
	for(i=1; i<=nn/2; ++i) {
		i1=nn-i;
		ct[0]=cf[0]*cwf[0]-cf[1]*cwf[1]; ct[1]=cf[0]*cwf[1]+cf[1]*cwf[0];
		c1[0]=0.5*(cg[i][0]+cg[i1][0])+ct[0]*(cg[i][0]-cg[i1][0])-
				ct[1]*(cg[i][1]+cg[i1][1]);
		c1[1]=0.5*(cg[i][1]-cg[i1][1])+ct[0]*(cg[i][1]+cg[i1][1])+
				ct[1]*(cg[i][0]-cg[i1][0]);
		cdiv(cf,cwf,c2);
		r=0.5*(cg[i1][0]+cg[i][0])-c2[0]*(cg[i1][0]-cg[i][0])+
				c2[1]*(cg[i1][1]+cg[i][1]);
		cg[i1][1]=0.5*(cg[i1][1]-cg[i][1])-c2[0]*(cg[i1][1]+cg[i][1])-
				c2[1]*(cg[i1][0]-cg[i][0]);
		cg[i1][0]=r;
		cg[i][0]=c1[0]; cg[i][1]=c1[1];
		r=cw[0]*cwf[0]-cw[1]*cwf[1]; cwf[1]=cw[0]*cwf[1]+cw[1]*cwf[0];
		cwf[0]=r;
	}

/*	The end points */
	gr=cg[0][0]; gi=cg[0][1];
	if(iflg>=0) {cg[0][0]=gr+gi; cg[0][1]=gr-gi;}
	else {
		cg[0][0]=0.5*(gr+gi); cg[0][1]=0.5*(gr-gi);
/*	Calculate the inverse DFT */
		ier=fft(nn,cg,iflg);
	}
	return ier;
}
