/*	Least squares polynomial fit using orthogonal polynomials in 2 dimensions
	Weights are assumed to be equal for all points

	NX : (input) Number of data points along X direction
	NY : (input) Number of data points along Y direction
	X,Y : (input) Arrays of length NX,NY containing the coordinates
		of tabular points
	F : (input) Array of length LA*NY containing the function values
		F[j][i] is the value at X[I],Y[J]
	AX : (output) Array of length IC*3 containing information about
		fit along X direction. AX[0][I], AX[1][I], AX[2][I] will
		respectively contain the coefficients, alpha, beta, gamma
	AY : (output) Array of length IC*3 containing information about
		fit along Y direction. AY[0][I], AY[1][I], AY[2][I] will
		respectively contain the coefficients, alpha, beta, gamma
	LA : (input) Second dimension of arrays  F, FY as declared
		in the calling function. LA >= MAX(NX,NY)
	C : (output) Array of length IC*(MY+1) containing the fitted
		coefficients of product of orthogonal polynomials in X & Y
	IC : (input) Second dimension of arrays C, AX, AY as declared in the calling
		function. IC > MAX(MX,MY)+1
	MX : (input) Required degree of polynomial in X
	MY : (input) Required degree of polynomial in Y
	FY : (output) Array of length LA*NY containing the values of
		fitted function at each tabular point	
	CHISQ : (output) the Chi square value for the fit
		
	Error status is returned by the value of the function POLFIT2.
		0 value implies successful execution
		602 implies that IC<MX+1 or IC<MY+1 or LA<NX or LA<NY
		603 implies that NX<MX+1 or MX<0 or NY < MY+1 or MY<0
		In both these cases no calculations are done
		Other values may be set by POLFIT1
	
	The fitted polynomial can be calculated at any value of x using POLEV2

	Required functions : POLFIT1, POLEV2, POLORT
*/

#include <math.h>
#include <stdlib.h>

int polfit1(int n, int m, int num, double x[], double *f, double w[],
	double *a, double alp[], double beta[], double gam[]);
int polev2(int nx, int ny, double *ax, double *ay, int la, double *wt,
	double x0, double y0, double *f, double *dfx, double *dfy,
	double *dfxx, double *dfxy, double *dfyy);

int polfit2(int nx, int ny, double x[], double y[], double *f, double *ax,
	double *ay, int la, double *c, int ic, int mx, int my, double *fy,
	double *chisq)

{
	int i,j,lj,ln1,m,num,ier;
	double dfx,dfy,dfxx,dfxy,dfyy,r1;
	double *wk,*wk1, *aw;

	if(nx>la || ny>la || ic<mx+1 || ic<my+1) return 602;

	if(mx<0 || my<0 || nx<=mx || ny<=my) return 603;

/*     Set the weights to 1 */
	lj=nx; if(lj<ny) lj=ny;
	aw=(double *) calloc((size_t) lj, sizeof(double));
	for(i=0; i<lj; ++i) aw[i]=1.0;

/*     Set up the RHS for calculating the fits along y-axis */
	wk=(double *) calloc((size_t) (nx*ny), sizeof(double));
	wk1=(double *) calloc((size_t) (nx*(my+1)), sizeof(double));
	lj=ny;
	for(i=0; i<nx; ++i) {
		for(j=0; j<ny; ++j) wk[j+i*lj]=f[i+j*la];
	}

	m=my+1;
	num=nx;
	ier=polfit1(ny,my,num,y,wk,aw,wk1,ay,&ay[ic],&ay[2*ic]);
	if(ier>100) { free(wk1); free(wk); free(aw); return ier;}

/*     Set up the RHS for calculating the fits along x-axis */
	for(j=0; j<m; ++j) {
		for(i=0; i<nx; ++i) wk[i+j*nx]=wk1[j+i*m];
	}
	num=m;
	ier=polfit1(nx,mx,num,x,wk,aw,wk1,ax,&ax[ic],&ax[2*ic]);
	if(ier>100) {free(wk1); free(wk); free(aw); return ier;}

/*	Store the calculated coefficients in array C */
	m=mx+1;
	for(i=0; i<=my; ++i) {
		for(j=0; j<m; ++j) c[j+ic*i]=wk1[j+i*m];
	}
 
/*     Calculate the CHI square */
	*chisq=0.0;
	for(i=0; i<ny; ++i) {
		for(j=0; j<nx; ++j) {
			ier=polev2(mx,my,ax,ay,ic,c,x[j],y[i],&fy[j+i*la],
				&dfx,&dfy,&dfxx,&dfxy,&dfyy);
			r1=f[j+i*la]-fy[j+i*la];
			*chisq=(*chisq)+r1*r1;
		}
	}
	free(wk1); free(wk); free(aw);
	return 0;
}
