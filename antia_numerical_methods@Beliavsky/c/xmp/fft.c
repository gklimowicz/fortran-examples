/*	Discrete Fourier transform of complex data */

#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979324

int fft(int n, double cg[][2], int iflg);
int dft(int n, double cg[][2], double cf[][2], int iflg);

main()
{
	int i,i1,j,n, id, iflg, ier,np;
	double hh,x,a,d1,d2,cg[2000][2],f[2000],cf1[2000][2],cf[2000][2];

/*	Example 10.5 */

	np=100; iflg=1;

/*	a is the coefficient in function used to generate data points */
	for(i1=0; i1<99; ++i1) {
		printf("type n = No. of pts,   a         (quits when n<=0)\n");
		scanf(" %d %le", &n, &a);
		if(n<=0) return 0;

		hh=1.0/n; iflg=1;
		for(i=0; i<n; ++i) {
			x=i*hh;
			cg[i][0]=sin(2.*PI*a*x);
			cg[i][1]=0.0;
			f[i]=cg[i][0];
		}

		i=dft(n,cg,cf,iflg);
		i=fft(n,cg,iflg); 
		printf(" ier = %d  n = %d    a = %e \n",i,n,a);
		printf("    Fourier transform (using fft)   Fourier transform (using dft)\n");
		for(i=0; i<n; ++i) printf(" %d   %e %e    %e %e \n",i,cg[i][0],cg[i][1],cf[i][0],cf[i][1]);

/*	Take inverse transform and compare with input data */
		iflg=-1;
		i=fft(n,cg,iflg); 
		i=dft(n,cf,cf1,iflg); 

		d1=0.0; d2=0.0;
		for(i=0; i<n; ++i) {
/*	The inverse transform should be divided by n before comparing */
			if(fabs(cg[i][0]/n-f[i]) > d1) d1=fabs(cg[i][0]/n-f[i]);
			if(fabs(cf1[i][0]/n-f[i]) > d2) d2=fabs(cf1[i][0]/n-f[i]);
		}
		printf(" maximum difference :  with fft = %e    with dft = %e\n",d1,d2);
	}
	return;
}



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
