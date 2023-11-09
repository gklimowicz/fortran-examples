/*	Multiple integration over a hyper-rectangle in n-dimensions
	using compound monomial rules with given number of points

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions, N>0 and N<=NMAX (=50)
	M : (input) Integer specifying the formula to be used
		M can be 1,3 or 5, otherwise no calculations are done
		M=1 selects 1-point formula of degree 1
		M=3 selects 2N-point formula of degree 3 due to Stroud
		M=5 selects (2N*N+1)-point formula of degree 5
	IND : (input) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise no calculations are done
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RI : (output) The calculated value of the integral
	NUM : (output) Number of function evaluations used by STROUD
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1000000)
		
	Error status is returned by the value of the function STROUD.
		0 value implies successful execution
		308 implies that number of points exceeded MAXPT and
			no calculations are done
		309 implies N<1 or N>=NMAX, in which case no calculations are done
		310 implies M is not 1,3 or 5 or IND[J]<1 for some J
			in which case no calculations are done

	Function F(N,X) must be provided by the user
	
	Required functions :  F
*/

#include <math.h>

#define PI 3.14159265358979324e0

int stroud(double a[], double b[], int n, int m, int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt)

{
	int i,j,k,mpt,in,ier,ip[50], nmax=50, maxpts=1100000;
	double r,s1,s2,a0,a1,a2,an,xi,xi2,xi3;
	double x[50],x3[50][50],h[50],xa[50],wt[50];

	*ri=0.0; *num=0;
	if(n>nmax || n<1) return 309;

/*	Calculate the number of function evaluations required */
	if(m == 1) mpt=1;
	else if(m == 3) mpt=2*n;
	else if(m == 5) mpt=2*n*n+1;
	else return 310;

	*num=mpt;
	for(i=0; i<n; ++i) {
		if(ind[i]<1) return 310;
		*num=(*num)*ind[i];
	}
	if(*maxpt <1) *maxpt=maxpts;
	if(*num > (*maxpt)) return 308;

/*	Constants for the (2N*N+1)-point formula of degree 5 */
	xi=sqrt(0.6e0);
	a0=(25*n*n-115*n+162.0)/162.0;
	a1=(70.0-25*n)/162.0;
	a2=25.0/324.0;

/*	Abscissas for the 2N-point formula of degree 3 */
	if(m == 3) {
		xi3=sqrt(2.0/3.0);
		xi2=1.0/sqrt(3.0);
		for(i=0; i<n; ++i) {
			for(j=0; j<n-1; j=j+2) {
				an=(j+1)*(i+1)*PI/n;
				x3[j][i]=xi3*cos(an);
				x3[j+1][i]=xi3*sin(an);
			}
/*	When N is odd */
			if(n != 2*(n/2)) {
				if(i == 2*(i/2)) x3[n-1][i]=-xi2;
				else x3[n-1][i]=xi2;
			}
		}
	}

	for(i=0; i<n; ++i) {
		ip[i]=1;
		h[i]=(b[i]-a[i])/(2*ind[i]);
/*	For abscissas of (2N*N+1)-point formula of degree 5 */
		wt[i]=h[i]*xi;
	}

/*	loop for the sum over all subintervals */
	k=n-1;
	do {
		for(in=k; in>=0; --in) xa[in]=a[in]+(2*ip[in]-1)*h[in];

/*	Generalised midpoint rule */
		if(m == 1) r=f(n,xa);

		else if(m==3) {
/*	Stroud's 2N-point rule of degree 3 */
			r=0.0;
			for(i=0; i<n; ++i) {
				for(j=0; j<n; ++j) x[j]=xa[j]+x3[j][i]*h[j];
				r=r+f(n,x);
				for(j=0; j<n; ++j) x[j]=xa[j]-x3[j][i]*h[j];
				r=r+f(n,x);
			}
			r=r/(2*n);
		}

		else if(m == 5) {
/*	(2N*N+1)-point rule of degree 5 */
			r=a0*f(n,xa);
			s1=0.0; s2=0.0;
			for(i=0; i<n; ++i) x[i]=xa[i];
			for(i=0; i<n; ++i) {
				x[i]=x[i]+wt[i];
				s1=s1+f(n,x);
				x[i]=xa[i]-wt[i];
				s1=s1+f(n,x);
				x[i]=xa[i];
				for(j=i+1; j<n; ++j) {
					x[i]=xa[i]+wt[i];
					x[j]=xa[j]+wt[j];
					s2=s2+f(n,x);
					x[j]=xa[j]-wt[j];
					s2=s2+f(n,x);
					x[i]=xa[i]-wt[i];
					s2=s2+f(n,x);
					x[j]=xa[j]+wt[j];
					s2=s2+f(n,x);
					x[j]=xa[j];
					x[i]=xa[i];
				}
			}
			r=r+a1*s1+a2*s2;
		}

		*ri=(*ri)+r;
		k=0;
		while(k <= n-1) {
			if(ip[k]>= ind[k]) {
				ip[k]=1;
				k=k+1;
			}
			else {
				ip[k]=ip[k]+1;
				break;
			}
		}
	} while (k<n);

/*	If all directions are exhausted, compute the value of integral */
	for(i=0; i<n; ++i) *ri=2.*(*ri)*h[i];
	return 0;
}
