/*	To calculate the discrete Fourier Transform using FFT algorithm

	N : (input) Number of points, which must be a power of 2
	CG : (input/output) Array of length 2*N containing the data points
		After execution it will contain the Fourier transform of CG
		CG is complex and hence the array should
		contain the real and imaginary parts.
	IFLG : (input) Flag to decide whether to calculate forward or inverse
		transform. If IFLG>=0 then Fourier transform is calculated
		IF IFLG<0 then inverse Fourier transform is calculated
		
	Error status is returned by the value of the function FFT.
		0 value implies successful execution
		611 implies that N<2, no calculations are done
		631 implies that N is not a power of 2, in this case
			contents of CG will be destroyed but will not
			contain the Fourier transform.

	Required functions : None
*/

#include <math.h>

#define PI 3.14159265358979324

int fft(int n, double cg[][2], int iflg)

{
	int i,j,m,j0,k0,iw,i1,jr;
	double r,ct[2],th,cwj[2],cwf[2];

	if(n<2) return 611;

/*	Bit reversal */
	j=0;
	for(i=0; i<n; ++i) {
		if(j>i) {
/*	exchange CG[I] with CG[J] */
			r=cg[i][0]; cg[i][0]=cg[j][0]; cg[j][0]=r;
			r=cg[i][1]; cg[i][1]=cg[j][1]; cg[j][1]=r;
		}
		m=n/2;
		while(m>=1 && j>=m) {j=j-m; m=m/2;}
/*	J-1 is the bit reverse of I */
		j=j+m;
	}

	j0=1; k0=n/2;
	th=PI/k0;
	if(iflg>=0) iw=1;
	else iw=-1;
	cwf[0]=-1; cwf[1]=0.0;

/*	Main loop for FFT executed Log_2(N) times */
	do {
		cwj[0]=1.0; cwj[1]=0.0;
/*	Inner loop over all elements */
		for(jr=0; jr<j0; ++jr) {
			for(i=jr; i<n; i += 2*j0) {
				i1=i+j0;
				ct[0]=cg[i1][0]*cwj[0]-cg[i1][1]*cwj[1];
				ct[1]=cg[i1][0]*cwj[1]+cg[i1][1]*cwj[0];
				cg[i1][0]=cg[i][0]-ct[0]; cg[i1][1]=cg[i][1]-ct[1];
				cg[i][0]=cg[i][0]+ct[0]; cg[i][1]=cg[i][1]+ct[1];
			}
			r=cwj[0]*cwf[0]-cwj[1]*cwf[1];
			cwj[1]=cwj[0]*cwf[1]+cwj[1]*cwf[0];
			cwj[0]=r;
		}

		j0=2*j0; k0=k0/2;
		cwf[0]=cos(iw*k0*th); cwf[1]=sin(iw*k0*th);
	} while(j0<n);
	if(j0==n) return 0;
	else return 631;
}
