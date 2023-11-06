/*	To integrate    f(x)*log(x)  over  (0,a] */

#include <stdio.h>
#include <math.h>

double fun(double x);
double f2(double x);
int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));
int gaulog(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt);
int gaulg2(double *rint, double a, double *a1, double reps, double aeps,
double *dif, double (*f) (double ), double (*f2) (double ), int *np);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double a1, xu, rint, reps, aeps, df;

/*	Exercise 6.2 (I5) */

	aeps=1.e-19; reps=1.e-14;
	for(i1=0; i1<99; ++i1) {
		printf("type  xu = upper limit,   a1 = point at which integral is to be split \n");
		printf("                      (quits when a1> 2*a)\n");
		scanf(" %le %le", &xu, &a1);
		if(a1> 2*xu) return 0;

		i=gaulg2(&rint,xu,&a1,reps,aeps,&df,fun,f2,&nuse);
		printf(" ier = %d  no. of function evaluations = %d  upper limit = %e   a1 = %e \n", i,nuse,xu,a1);
		printf(" integral = %e   estimated error = %e  \n", rint, df);
	}
	return;
}

/* the integrand */

double fun(double x)

{
	return -log(x)*sin(x);
}

/* the integrand/log(1/x)  for gaulog */

double f2(double x)

{
	return sin(x);
}

 

/*     To integrate a function with logarithmic singularity over (0,A]
     using a combination of Gaussian formulas

	RINT : (output) Calculated value of the integral
	A : (input) The upper limit
	A1 : (input/output) The point at which integral has to be broken
		A1 will be adjusted by the function.
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	F2 : (input) Name of the function to calculate F(X)/LOG(X)
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function GAULG2.
		0 value implies successful execution
		31 implies specified accuracy was not achieved by GAUSS over [A1,A]
		32 implies specified accuracy was not achieved by GAULOG
		34 implies specified accuracy was not achieved by GAUSS over (0,A1]
		In case of multiple failures second digit will be sum
			of these values.
		In all cases DIF will contain the estimated accuracy

     	FUNCTION F(X) and F2(X) must be supplied by the user.

	Required functions : GAUSS, GAULOG, F, F2
*/

#include <math.h>

int gaulog(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt);
int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));


int gaulg2(double *rint, double a, double *a1, double reps, double aeps,
	double *dif, double (*f) (double ), double (*f2) (double ), int *np)

{
	int npt,npt1,npt2,ier,ier1,ier2;
	double r1,r2,r3,a2,dif1,dif2,t1, amn=0.01;

	r1=0.0; r2=0.0; r3=0.0; a2=0.0; *dif=0.0; dif2=0.0;
	*np=0; npt=0; npt2=0; ier1=0; ier2=0;
	if( *a1 > a) *a1=a;
	if( *a1<= 0.0) *a1=a;
 
/*     Evaluate the integral over (0,A1] */
	do {
		ier=gaulog(&r1,*a1,aeps,reps,&dif1,f2,&npt1);
		*np=(*np)+npt1;
		if(ier == 0) break;

/*     If GAULOG fails decrease A1 */
		t1= (*a1);
		*a1=(*a1)/2.;
	} while (*a1 > amn);
 
	if(ier>0) {ier1=2; *a1=t1;}

/*     Evaluate the integral over [A1,A] */
	npt1=16;
	ier=0;
	if(a-(*a1) > aeps) ier=gauss(&r2,*a1,a,&npt1,reps,aeps,dif,&npt,f);
	if(ier>0) ier=1;

/*     Evaluate the regular part over [0,A1] */
	if((*a1) != 1) ier2=gauss(&r3,a2,*a1,&npt1,reps,aeps,&dif2,&npt2,f2);
	if(ier2>0) ier2=4;
	*rint=r1+r2-r3*log(*a1);
	*dif=fabs(*dif)+fabs(dif1)+fabs(dif2);
	*np = (*np)+npt+npt2;
	ier=ier+ier1+ier2;
	if(ier>0) ier=ier+30;
	return ier;
}



/*     To integrate a function with logarithmic singularity Gaussian formulas

	RINT : (output) Calculated value of the integral
	A : (input) The upper limit
	AEPS : (input) The required absolute accuracy
	REPS : (input) The required relative accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the
		integrand (divided by LOG(A/X))
	NPT : (output) Number of function evaluations used
		
	Error status is returned by the value of the function GAULOG.
		0 value implies successful execution
     		30 implies specified accuracy was not achieved
     		DIF will contain the estimated accuracy

	Function F(X) must be supplied by the user
	Note that function calculates integral of F(X)*LOG(A/X)

	Required functions : F
*/

#include <math.h>

int gaulog(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt)

{
	int i,j,n;
	double r1,s1,r2;

/*     Weights and abscissas for Gaussian formula for logarithmic singularity
     W[N-2],...,W[2N-3], are the weights for N-point rule and
     X[N-2],...,X[2N-3], the corresponding abscissas
     Weights and abscissas are available for N=2,4,8,16
*/
 
	double x[30] = {1.120088061669761830e-1, 6.022769081187381028e-1,
			4.144848019938322080e-2, 2.452749143206022519e-1,
                        5.561654535602758372e-1, 8.489823945329851746e-1,
                        1.332024416089246501e-2, 7.975042901389493841e-2,
                        1.978710293261880538e-1, 3.541539943519094197e-1,
                        5.294585752349172777e-1, 7.018145299390999638e-1,
                        8.493793204411066760e-1, 9.533264500563597888e-1,
                        3.897834487115909095e-3, 2.302894561687320045e-2,
                        5.828039830624031972e-2, 1.086783650910538817e-1,
                        1.726094549098437244e-1, 2.479370544705782363e-1,
                        3.320945491299168705e-1, 4.221839105819483085e-1,
                        5.150824733814623250e-1, 6.075561204477284747e-1,
                        6.963756532282138523e-1, 7.784325658732652431e-1,
                        8.508502697153909688e-1, 9.110868572222718348e-1,
                        9.570255717035421226e-1, 9.870478002479844660e-1};

	double w[30] = {7.185393190303844407e-1,2.814606809696155593e-1,
                        3.834640681451351249e-1, 3.868753177747626273e-1,
                        1.904351269501424154e-1, 3.922548712995983245e-2,
                        1.644166047280028868e-1, 2.375256100233060205e-1,
                        2.268419844319191264e-1, 1.757540790060702450e-1,
                        1.129240302467590519e-1, 5.787221071778207240e-2,
                        2.097907374213297804e-2, 3.686407104027619013e-3,
                        6.079171004359114509e-2, 1.029156775175820228e-1,
                        1.223556620460090919e-1, 1.275692469370159323e-1,
                        1.230135746000709083e-1, 1.118472448554857552e-1,
                        9.659638515212439849e-2, 7.935666435147320573e-2,
                        6.185049458196527197e-2, 4.543524650772672381e-2,
                        3.109897475158184829e-2, 1.945976592736087029e-2,
                        1.077625496320554213e-2, 4.972542890087649610e-3,
                        1.678201110051197249e-3, 2.823537646684367889e-4};
 
/*     The 2-point formula */
	r1=(f(a*x[0])*w[0]+f(a*x[1])*w[1])*a;
	*npt=2;
	n=2;
 
/*     Use higher order formula until convergence */
	for(j=2; j<=4; ++j) {
		n=n*2;
		r2=0.0;
		for(i=n-2; i<2*n-2; ++i) r2=r2+f(x[i]*a)*w[i];
		r2=r2*a;
 
		*npt=(*npt)+n;
		*dif=r2-r1;
		*rint=r2;
		s1=reps*fabs(*rint); if(aeps>s1) s1=aeps;
		if(fabs(*dif) < s1) return 0;
		r1=r2;
	}
 
	return 30;
}



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
