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
