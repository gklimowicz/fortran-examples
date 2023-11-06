/*	To integrate a function over finite interval using Gauss-Legendre formulas

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	B : (input) The upper limit
	NP : (input/output) The Gauss-Legendre formula to be used, (NP=2,4,8,16,32)
		For other values of NP it will be set to 8.
		It will use composite NP-point formula
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	NPT : (output) Number of function evaluations used by the function
	FUN : (input) Name of the function to calculate the integrand
		
	Error status is returned by the value of the function GAUSS.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy
		36 implies that NP was not 2,4,8,16 or 32. In which case
			it is set to 8.

		FUNCTION FUN(X) must be supplied by the user.

	Required functions : FUN
*/

#include <math.h>

int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ))

{
	int i,j,k,n,no,np2,ier;
	double r1,dx,a1,s1,at,bt, nmax=9;

/*	Weights and abscissas for Gauss-Legendre quadrature.
	For N-point formula W[K]=W[N-K-1] and X[K]=-X[N-K-1]
		For K=0,1,...,N/2-1. Hence only half points are tabulated.
	For 2-point W[0]; 4-point W[1], W[2]; 8-point W[3],...,W[6];
	16-point W[7],...,W[14]; 32-point W[15],...,W[30] are the
	weights corresponding to abscissas X[I].
*/

	double w[31]={1.0e0,
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

	n=1;
	dx=b-a;
	ier=0;
	*rint=0.0;
	*npt=0;

	if(*np==2) no=0;
	else if(*np==4) no=1;
	else if(*np==8) no=3;
	else if(*np==16) no=7;
	else if(*np==32) no=15;
	else {
/*	If NP-point formula is not available use NP=8 */
		*np=8; no=3; ier=36;
	}
/*	X[NO],...,X[NP2] are the abscissas for the formula */
	np2=no+(*np)/2-1;

/*	Subdivide the interval until convergence */
	for(i=0; i<nmax; ++i) {
		r1=0.0;
		for(j=0; j<n; ++j) {
			a1=a+j*dx;
			at=dx/2.;
			bt=a1+dx/2.;

/*	To reduce roundoff errors sum over each subinterval is evaluated separately */
			s1=0.0;
			for(k=no;  k<=np2; ++k) s1=s1+w[k]*(fun(at*x[k]+bt)+fun(bt-at*x[k]));
			r1=r1+s1;
		}
		r1=r1*dx/2.;

/*	convergence check */
		*dif=r1-(*rint);
		*rint=r1;
		*npt=(*npt)+n*(*np);
		s1=reps*fabs(*rint); if(aeps>s1) s1=aeps;
		if(i>0 && fabs(*dif) < s1) return ier;
		dx=dx/2.;
		n=n*2;
	}
	return 30;
}
