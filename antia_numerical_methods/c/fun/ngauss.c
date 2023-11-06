/*	Multiple integration over a hyper-rectangle in n-dimensions
	using product Gauss-Legendre formulas with given number of points

	A : (input) Array of length N containing the lower limit
		along each dimension
	B : (input) Array of length N containing the upper limit
		along each dimension
	N : (input) The number of dimensions, N>0
	M : (input) Integer array of length N specifying the formula
		to be used along each dimension. M[J]-point formula will
		be used along Jth dimension, M[J] should be 2,4,8,16 or 32
		otherwise no calculations are done
	IND : (input) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise no calculations are done
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RI : (output) The calculated value of the integral
	NUM : (output) Number of function evaluations used
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
		
	Error status is returned by the value of the function NGAUSS.
		0 value implies successful execution
		305 implies N<1, in which case no calculations are done
		306 implies M[J] is not 2,4,8,16 or 32 or IND[J]<1 for some J
			in which case no calculations are done
		307 implies that number of points exceeded MAXPT and
			no calculations are done
	
	Required functions : F
*/

#include <math.h>
#include <stdlib.h>

#define MAXPTS 1100000

int ngauss(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt)

{
	int i,j,no,j1,j2,m1,m2,k,ier, maxpts=1100000;
	double x1,h1;
	int *ip, *npt;
	double *h, *xa, *wt;

/*	Weights and abscissas for Gauss-Legendre quadrature.
	For N-point formula W[K]=W[N-K-1] and X[K]=-X[N-K-1]
		For K=0,1,...,N/2-1. Hence only half points are tabulated.
	For 2-point W[0]; 4-point W[1], W[2]; 8-point W[3],...,W[6];
	16-point W[7],...,W[14]; 32-point W[15],...,W[30] are the
	weights corresponding to abscissas X[I].
*/

	double w[31] = {1.0e0,
      0.34785484513745385737e0, 0.65214515486254614263e0,
      0.10122853629037625915e0, 0.22238103445337447054e0,
      0.31370664587788728734e0, 0.36268378337836198297e0,
      0.02715245941175409485e0, 0.06225352393864789286e0,
      0.09515851168249278481e0, 0.12462897125553387205e0,
      0.14959598881657673208e0, 0.16915651939500253819e0,
      0.18260341504492358887e0, 0.18945061045506849629e0,
      0.00701861000947009660e0, 0.01627439473090567061e0,
      0.02539206530926205945e0, 0.03427386291302143310e0,
      0.04283589802222668066e0, 0.05099805926237617620e0,
      0.05868409347853554714e0, 0.06582222277636184684e0,
      0.07234579410884850623e0, 0.07819389578707030647e0,
      0.08331192422694675522e0, 0.08765209300440381114e0,
      0.09117387869576388471e0, 0.09384439908080456564e0,
      0.09563872007927485942e0, 0.09654008851472780057e0};

	double x[31] = {0.57735026918962576451e0,
      0.86113631159405257522e0, 0.33998104358485626480e0,
      0.96028985649753623168e0, 0.79666647741362673959e0,
      0.52553240991632898582e0, 0.18343464249564980494e0,
      0.98940093499164993260e0, 0.94457502307323257608e0,
      0.86563120238783174388e0, 0.75540440835500303390e0,
      0.61787624440264374845e0, 0.45801677765722738634e0,
      0.28160355077925891323e0, 0.09501250983763744019e0,
      0.99726386184948156354e0, 0.98561151154526833540e0,
      0.96476225558750643077e0, 0.93490607593773968917e0,
      0.89632115576605212397e0, 0.84936761373256997013e0,
      0.79448379596794240696e0, 0.73218211874028968039e0,
      0.66304426693021520098e0, 0.58771575724076232904e0,
      0.50689990893222939002e0, 0.42135127613063534536e0,
      0.33186860228212764978e0, 0.23928736225213707454e0,
      0.14447196158279649349e0, 0.04830766568773831623e0};

	*ri=0.0; *num=0;
	if(n<1) return 305;
	if(*maxpt < 1) *maxpt=maxpts;

/*	calculate the number of function evaluations required */
	*num=m[0]*ind[0];
	for(i=1; i<n; ++i) *num=(*num)*m[i]*ind[i];
	if(*num > (*maxpt)) return 307;

/*	Initialisation */
	ip=(int *) calloc((size_t) (n+2), sizeof(int));
	npt=(int *) calloc((size_t) (n+2), sizeof(int));
	h=(double *) calloc((size_t) (n+2), sizeof(double));
	xa=(double *) calloc((size_t) (n+2), sizeof(double));
	wt=(double *) calloc((size_t) (n+2), sizeof(double));
	ier=0;
	for(i=0; i<n; ++i) {
		ip[i]=0;
		npt[i]=m[i]*ind[i]-1;
		if((m[i]!=2) && (m[i]!=4) && (m[i]!=8) && (m[i]!=16) && (m[i]!=32)) ier=306;
		if(ind[i] < 1) ier=306;
		h[i]=(b[i]-a[i])/(2*ind[i]);
	}
	if(ier!=0) {free(wt); free(xa); free(h); free(npt); free(ip); return ier;}
	for(i=n; i<n+1; ++i) {h[i]=1.0; wt[i]=1.0;}

/*	Loop for sum over N dimensions */
	k=n-1;
	do {
		for(i=k; i>=0; --i) {
			m2=m[i]/2;
/*	The abscissas are X(NO),...,X(NO+M2-1) */
			no=m2-1;
			h1=h[i];
			j1=ip[i]/m[i];
			j2=ip[i]-j1*m[i];
/*	Use the (J2)th point in (J1)th subinterval */
			x1=a[i]+(2*j1+1)*h1;
			if(j2<m2) {
/*	For the first M2 abscissas */
				xa[i]=x1+h1*x[no+j2];
				wt[i]=w[no+j2]*wt[i+1];
			}
			else if(j2-m2 < m2) {
/*	For the next M2 abscissas */
				xa[i]=x1-h1*x[no+j2-m2];
				wt[i]=w[no+j2-m2]*wt[i+1];
			}
			else {
/*	For Gaussian formula with odd number of points use the abscissa at x=0 */
				xa[i]=x1;
				wt[i]=w[no+m2]*wt[i+1];
			}
		}

/*	Add the new point to running sum */
		*ri=(*ri)+wt[0]*f(n,xa);
		k=0;
		while(k <= n-1) {
			if(ip[k] >= npt[k]) {
				ip[k]=0;
				k=k+1;
			}
			else {
				ip[k]=ip[k]+1;
				break;
			}
		}
	} while(k<n);

/*	If all points are exhausted compute the value of integral */
	for(i=0; i<n; ++i) *ri=(*ri)*h[i];
	free(wt); free(xa); free(h); free(npt); free(ip);
	return 0;
}
