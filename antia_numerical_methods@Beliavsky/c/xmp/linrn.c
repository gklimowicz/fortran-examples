/*	Linear interpolation in n-dimensions */

#include <stdio.h>
#include <math.h>
int linrn(int n, double xb[], double x[], double f[], int np[], double *fb,
	int ndim[], int nxd);
int locate(double xb, double x[], int np);
double fun(double x1, double x2, double x3, double x4, double x5);

main()
{
	int i,i1,i2,i3,i4,j,n,nxd,n1[10],ndim[10];
	double fb, h, f0, xb[10], x[10][20], f[9][9][9][9][9];

	ndim[0]=ndim[1]=ndim[2]=ndim[3]=ndim[4]=9;
	n1[0]=n1[1]=n1[2]=n1[3]=n1[4]=5;
	h=0.1;
	n=5;
	nxd=20;

/*	generating the table of values in four dimensions with uniform
	spacing of h=0.1 and containing n=5 points along each axes
	spanning the interval [0,0.4]    */

	for(i=0; i<n; ++i) {
		for(j=0; j<n1[0]; ++j) x[i][j]=j*h;
	}
	for(i=0; i<n1[0]; ++i) {
		for(i1=0; i1<n1[1]; ++i1) {
			for(i2=0; i2<n1[2]; ++i2) {
				for(i3=0; i3<n1[3]; ++i3) {
					for(i4=0; i4<n1[4]; ++i4) {
						f[i4][i3][i2][i1][i]=fun(x[0][i], x[1][i1], x[2][i2], x[3][i3], x[4][i4]);}}}}}

	for(i1=0; i1<99; ++i1) {
		printf("type xb(i), i=0,4  the coordinates of the required point\n");
		printf("             (quits when xb[0]<-100)\n");
		scanf("%le %le %le %le %le", &xb[0], &xb[1], &xb[2], &xb[3], &xb[4]);
		if(xb[0]<-100) return 0;
		i=linrn(n,xb, &x[0][0], &f[0][0][0][0][0], n1, &fb, ndim, nxd);
		printf(" ier = %d  x = %e %e %e %e %e \n", i,xb[0],xb[1],xb[2],xb[3],xb[4]);
		f0=fun(xb[0],xb[1],xb[2],xb[3],xb[4]);
		printf("  interpolated value = %e    exact value = %e \n", fb,f0);

	}
	return;
}

/*	sample function to generate the table for interpolation 
    since this is linear function interpolation should give exact
	value apart from roundoff  */

double fun(double x1, double x2, double x3, double x4, double x5)

{
	double f;

	f=(1.+x1)*(2.+x2)*(3.+x3)*(x4+x5);
	return f;
}
	


/*	Linear interpolation in N dimensions

	N : (input) Number of variables
	XB : (input) Array of length N containing the coordinates
		of the point at which interpolation is required
	X : (input) Array of length NXD*N containing the abscissas
	F : (input) Array of dimension exactly F[NDIM[N-1]]...[NDIM[0]]
		containing the function values
		F[IN_1]...[I1][I0]=f(X[0][I0],X[1][I1],...,X[N-1][IN_1])
	NP : (input) Integer array of length N containing the number of tabular points
		NP[I] is the number of tabular point along Ith dimension
	FB : (output) The interpolated value
	NDIM : (input) Integer array of length N, containing the
		dimension of F as declared in calling function
	NXD : (input) Second dimension of array X, NXD >= max(NP[i])
		
	Error status is returned by the value of the function LINRN.
		0 value implies successful execution
		207 implies NP[I]>NDIM[I] or NP[I]<2 for some I

	Required functions : LOCATE
*/

#include <math.h>
#include <stdlib.h>

int locate(double xb, double x[], int np);

int linrn(int n, double xb[], double x[], double f[], int np[], double *fb,
	int ndim[], int nxd)

{
	int i,j,k, nd, low, ndp, index;
	double term;
	int *jn;
	double *hn;

	nd=0;
	jn=(int *) calloc((size_t) (2*n), sizeof(int));
	hn=(double *) calloc((size_t) n, sizeof(double));

/*	Locate the point within a hypercube in the table */

	for(i=0; i<n; ++i) {
		if(np[i]>ndim[i] || np[i]>nxd || np[i]<2) {free(hn); free(jn); return 207;}

		low=locate(xb[i],&x[nd],np[i]);
		jn[i+n]=low;
		hn[i]=(xb[i]-x[low+nd])/(x[low+nd+1]-x[low+nd]);
		nd=nd+nxd;
		jn[i]=0;
	}

	*fb=0.0;
	do {
		index=jn[0]+jn[n];
		ndp=ndim[0];
		for(i=1; i<n; ++i) {
			index=index+(jn[i]+jn[i+n])*ndp;
			ndp=ndp*ndim[i];
		}
	
/*	F(IN[0]+IN[N], IN[1]+IN[1+N], ... ,IN[N-1]+IN[2N-1]) */
		term=f[index];

		for(i=0; i<n; ++i) {
			if(jn[i] == 0) term=term*(1-hn[i]);
			else term=term*hn[i];
		}

		*fb = *fb + term;

/*	select the next point */

		k=0;
		while(k <= n-1) {
			if(jn[k] >= 1) {
				jn[k]=0;
				k=k+1;
			}
			else {
				jn[k]=jn[k]+1;
				break;
			}
		}
	} while(k<n);

	free(hn); free(jn);
	return 0;
}


/*	To locate a given point between two points of an ordered table

	XB : (input) The point which needs to be located
	X : (input) Array of length NP containing the ordered table
	NP : (input) Length of table
	LOCATE should give the value such that XB is between X[LOCATE] and 
		X[LOCATE+1]

	Required functions : None
*/

#include <math.h>

int locate(double xb, double x[], int np)

{
	int low, igh, mid; 

	low=0;
	igh=np-1;
	if( (xb<x[low]) == (xb<x[igh]) ) {
		if(fabs(xb-x[low]) > fabs(xb-x[igh])) low=igh-1;
	}
	else {

/*	If the point is within the range of table locate it by bisection */
		while(igh-low > 1) {
			mid=(low+igh)/2;
			if( (xb<x[mid]) == (xb<x[low])) low=mid;
			else igh=mid;
		}
	}
	return low;
}
