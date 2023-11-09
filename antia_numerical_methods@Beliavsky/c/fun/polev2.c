/*	Evaluating the fitted polynomial and its derivatives at any value
	of x using known coefficients of orthogonal polynomials in 2 dimensions
	Should be used to evaluate the polynomial using coefficients calculated
	by POLFIT2.

	NX : (input) Degree of polynomial in X
	NY : (input) Degree of polynomial in Y
	AX : (input) Array of length LA*2 containing the coefficients
		alpha and beta for orthogonal polynomials in X
		AX[0][I] contains alpha and AX[1][I] contains beta
	AY : (input) Array of length LA*2 containing the coefficients
		alpha and beta for orthogonal polynomials in Y
		AY[0][I] contains alpha and AY[1][I] contains beta
		The arrays AX and AY can be calculated using POLFIT2
	LA : (input) Second dimension of arrays AX, AY and WT in the calling
		function. LA > MAX(NX,NY)
	WT : (input) Array of length LA*(MY+1) containing the coefficients
		of the fit. WT[J][I] is the coefficient of
		PHI_I(X)PSI_J(Y), where PHI_I and PSI_J are orthogonal
		polynomials in X and Y
	X0,Y0 : (input) Coordinates of the point at which polynomial
		needs to be evaluated
	F : (output) Calculated value of the fitted polynomial at (X0,Y0)
	DFX : (output) First derivative  dF/dX at X0,Y0
	DFY : (output) First derivative dF/dY at X0,Y0
	DFXX : (output) Second derivative d^2F/dXdX at X0,Y0
	DFXY : (output) Second derivative d^2F/dXdY at X0,Y0
	DFYY : (output) Second derivative d^2F/dYdY at X0,Y0
		
	Error status is returned by the value of the function POLEV2.
		0 value implies successful execution
		604 implies that LA<=MAX(NX,NY), in which case no
			calculations are done.
	
	Required functions : POLORT
*/

#include <math.h>
#include <stdlib.h>

int polort(int m, double alp[], double beta[], double x, double f[],
	double df[], double ddf[]);

int polev2(int nx, int ny, double *ax, double *ay, int la, double *wt,
	double x0, double y0, double *f, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy)

{
	int i,j,nk;
	double *wkx,*wky;

	nk=nx+1; if(ny>nx) nk=ny+1;
	if(nk-1>la) return 604;

	wkx=(double *) calloc((size_t) (3*nk),sizeof(double));
	wky=(double *) calloc((size_t) (3*nk),sizeof(double));

/*	Calculate the orthogonal polynomials along each dimension */
	i=polort(nx,ax,&ax[la],x0,wkx,&wkx[nk],&wkx[2*nk]);
	i=polort(ny,ay,&ay[la],y0,wky,&wky[nk],&wky[2*nk]);
 
/*	Calculate the fitted polynomial and its derivatives */
	*f=0.0; *dfx=0.0; *dfy=0.0;
	*dfxx=0.0; *dfxy=0.0; *dfyy=0.0;
	for(i=0; i<=nx; ++i) {
		for(j=0; j<=ny; ++j) {
			*f=(*f)+wt[i+j*la]*wkx[i]*wky[j];
			*dfx=(*dfx)+wt[i+j*la]*wkx[nk+i]*wky[j];
			*dfy=(*dfy)+wt[i+j*la]*wkx[i]*wky[nk+j];
			*dfxx=(*dfxx)+wt[i+j*la]*wkx[2*nk+i]*wky[j];
			*dfyy=(*dfyy)+wt[i+j*la]*wkx[i]*wky[2*nk+j];
			*dfxy=(*dfxy)+wt[i+j*la]*wkx[nk+i]*wky[nk+j];
		}
	}
	free(wky); free(wkx);
	return 0;
}
