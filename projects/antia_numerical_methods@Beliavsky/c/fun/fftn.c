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
