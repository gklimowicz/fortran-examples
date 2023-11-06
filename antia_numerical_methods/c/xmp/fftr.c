/*	Fast Fourier transform of real data */

#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979324

void cdiv(double c1[2], double c2[2], double cr[2]);

int fftr(int n, double cg[][2], int iflg);

main()
{
	int i,i1,j,n, id, iflg, ier,np;
	double hh,x,a,d1,d2,cg[2000],f[2000];

/*	Example 10.5 */

	iflg=1;

/*	a is the coefficient in function used to generate data points */
	for(i1=0; i1<99; ++i1) {
		printf("type n = No. of pts,   a         (quits when n<=0)\n");
		scanf(" %d %le", &n, &a);
		if(n<=0) return 0;

		hh=1.0/n; iflg=1;
		for(i=0; i<n; ++i) {
			x=i*hh;
			cg[i]=sin(2.*PI*a*x);
			f[i]=cg[i];
		}

		i=fftr(n,&cg[0],iflg); 
		printf(" ier = %d  n = %d    a = %e \n",i,n,a);
		printf("    Fourier transform \n");
		for(i=0; i<n; i=i+2) printf(" %e %e \n",cg[i],cg[i+1]);

/*	Take inverse transform and compare with input data */
		iflg=-1;
		i=fftr(n,&cg[0],iflg); 

		d1=0.0;
		for(i=0; i<n; ++i) {
/*	The inverse transform should be divided by n/2 before comparing */
			if(fabs(cg[i]*2/n-f[i]) > d1) d1=fabs(cg[i]*2/n-f[i]);
		}
		printf(" maximum difference  = %e \n",d1);
	}
	return;
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


/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	Complex division

	C1 : (input) Array of length 2 containing the numerator 
	C2 : (input) Array of length 2 containing the denominator
	CR : (output) Array of length 2 containing the result (C1/C2)
		CR can be same as  C1 or C2

	Required functions : None

*/

#include <math.h>

void cdiv(double c1[2], double c2[2], double cr[2])

{
	double r,r1,den;

	if(fabs(c2[0])>fabs(c2[1])) {
		r=c2[1]/c2[0];
		den=c2[0]+c2[1]*r;
/*	To avoid overwriting on c1 if c1 and cr are same */
		r1=(c1[0]+c1[1]*r)/den;
		cr[1]=(c1[1]-c1[0]*r)/den;
		cr[0]=r1;
	}
	else {
		r=c2[0]/c2[1];
		den=c2[0]*r+c2[1];
		r1=(c1[0]*r+c1[1])/den;
		cr[1]=(c1[1]*r-c1[0])/den;
		cr[0]=r1;
	}
	return;
}
