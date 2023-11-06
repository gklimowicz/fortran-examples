/*	Fast Fourier transform of real data in n dimensions */

#include <stdio.h>
#include <math.h>

#define PI 3.14159265358979324

int fftn(int nd, int nn[], double cg[][2], int iflg);
double cabs(double cx[2]);

main()
{
	int i,i1,j,k,n, id, iflg, ier,nn[9];
	double hh,x,y,z,a,d1,cg[34000][2],f[34000];

	id=3; iflg=1;

	for(i1=0; i1<99; ++i1) {
		printf("type nn[0],...,nn[id-1] = No. of points along each axis\n");
		printf("                        (quits when nn[0]<=0)\n");
		for(i=0; i<id; ++i) scanf(" %d", &nn[i]);
		if(nn[0]<=0) return 0;

		printf("type A = coef in sin and cos terms \n");
		scanf(" %le", &a);
		hh=1.0/nn[0];
		for(i=0; i<nn[0]; ++i) {
			x=i*hh;
			for(j=0; j<nn[1]; ++j) {
				y=j*hh;
				for(k=0; k<nn[2]; ++k) {
					z=k*hh;
					n=i+j*nn[0]+k*nn[0]*nn[1];
					cg[n][0]=sin(2.*PI*a*y)*cos(2.*PI*a*x)*cos(2.*PI*a*z);
					cg[n][1]=0.0;
					f[n]=cg[n][0];
				}
			}
		}

		iflg=1;
		i=fftn(id,nn,cg,iflg);
		printf(" ier = %d    no. of points = %d %d %d    a = %e\n",i,nn[0],nn[1],nn[2],a);
		n=nn[0]*nn[1]*nn[2];
		for(i=0; i<n; ++i) if(cabs(&cg[i][0]) > 1.e-4) printf(" %d %e %e\n",i,cg[i][0],cg[i][1]);

/*	Take the inverse transform and compare with input data */
		iflg=-1;
		i=fftn(id,nn,cg,iflg);
		printf(" ier = %d  \n",i);
		d1=0.0;
		for(i=0; i<n; ++i) {
			if(fabs(cg[i][0]/n-f[i])+fabs(cg[i][1]) > d1) d1=fabs(cg[i][0]/n-f[i])+fabs(cg[i][1]);
		}
	printf(" maximum difference = %e \n",d1);
	}
	return;
}



/*	To calculate the discrete Fourier Transform using FFT algorithm in n dimensions

	ND : (input) Number of dimensions
	NN : (input) Integer array of length ND containing the number of points
		along each dimension. NN[I] is the number of points along
		the Ith dimension, which should be a power of 2
	CG : (input/output) Array of length 2*NN[0]*NN[1]*...*NN[ND-1]
		containing the data points. The dimensions of CG in the calling
		function must exactly match the number of points, e.g.
		CG[NN[ND-1]]...[NN[1]][NN[0]][2]
		The real and imaginary part of data must be stored in consecutive
		locations.
		After execution it will contain the Fourier transform of CG
	IFLG : (input) Flag to decide whether to calculate forward or inverse
		transform. If IFLG>=0 then Fourier transform is calculated
		IF IFLG<0 then inverse Fourier transform is calculated
		
	Error status is returned by the value of the function FFTN.
		0 value implies successful execution
		631 implies that at-least one of NN[I] is not a power of 2,
			in this case contents of CG will be destroyed but will not
			contain the Fourier transform.

	Required functions : None
*/

#include <math.h>

#define PI 3.14159265358979324

int fftn(int nd, int nn[], double cg[][2], int iflg)

{
	int i,j,i1,i2,id,j0,j2,k0,iw,m,n,npr,npr1,ntot,ir,jr,jr0;
	double ct[2],r,th,cwf[2],cwj[2];

	ntot=1;
	for(i=0; i<nd; ++i) ntot=ntot*nn[i];
	
	npr1=1;
	if(iflg>=0) iw=1;
	else iw=-1;

/*	Loop over each dimension */
	for(id=0; id<nd; ++id) {
		n=nn[id];
		npr=npr1;
		npr1=npr*n;

/*	Loop for bit reversal */
		j=0;
		for(i=0; i<npr1; i += npr) {
			if(j>i) {
				for(i1=i; i1<ntot; i1 += npr1) {
					for(i2=i1; i2<=i1+npr-1; ++i2) {
						j2=i2+j-i;
						r=cg[i2][0]; cg[i2][0]=cg[j2][0]; cg[j2][0]=r;
						r=cg[i2][1]; cg[i2][1]=cg[j2][1]; cg[j2][1]=r;
					}
				}
			}
			m=npr1/2;
			while(m>=npr && j>=m) {j=j-m; m=m/2;}
			j=j+m;
		}

		j0=1; k0=n/2;
		th=PI/k0;
		cwf[0]=-1.0; cwf[1]=0.0;

/*	Loop for FFT calculation */
		do {
			cwj[0]=1.0; cwj[1]=0.0;
			for(jr=0; jr<j0; ++jr) {
				jr0=jr*npr;
				for(ir=jr0; ir<ntot; ir += 2*j0*npr) {
					for(i=ir; i<=ir+npr-1; ++i) {
						i1=i+j0*npr;
						ct[0]=cg[i1][0]*cwj[0]-cg[i1][1]*cwj[1];
						ct[1]=cg[i1][0]*cwj[1]+cg[i1][1]*cwj[0];
						cg[i1][0]=cg[i][0]-ct[0]; cg[i1][1]=cg[i][1]-ct[1];
						cg[i][0]=cg[i][0]+ct[0]; cg[i][1]=cg[i][1]+ct[1];
					}
				}
				r=cwj[0]*cwf[0]-cwj[1]*cwf[1];
				cwj[1]=cwj[0]*cwf[1]+cwj[1]*cwf[0];
				cwj[0]=r;
			}

			j0=2*j0; k0=k0/2;
			if(j0==n) break;
			if(j0>n || k0==0) return 631;

			cwf[0]=cos(iw*k0*th); cwf[1]=sin(iw*k0*th);
		} while(j0<n);
	}
	return 0;
}



/* Utility functions for complex arithmetic, complex numbers are
	represented as double array of length 2, with first element as
	real part and second element as imaginary part.

*/

/*	For finding absolute value of a complex number

	CX : (input) Array of length 2 containing the complex number
	CABS will be the absolute value of CX

	Required functions : None

*/

#include <math.h>

double cabs(double cx[2])

{
	double r;

	if(fabs(cx[0])>fabs(cx[1])) {
		r=cx[1]/cx[0];
		return fabs(cx[0])*sqrt(1+r*r);
	}
	else if(cx[1] != 0.0) {
		r=cx[0]/cx[1];
		return fabs(cx[1])*sqrt(1+r*r);
	}
	else return 0.0;

}
