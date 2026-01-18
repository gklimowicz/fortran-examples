/*	To integrate  f(x)*exp(-x) over semi infinite interval */

#include <stdio.h>
#include <math.h>

double fun(double x);
double f2(double x);
int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));
int gaulag(double *rint, double a, double *a1, double reps, double aeps,
	double *dif, double (*f) (double ), double (*f2) (double ), int *np);
int lagure(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt);

main()
{
	int i,i1,j,nuse, id, iflg, ier,np;
	double xl, a1, rint, reps, aeps, df, fb[20];

/*	Example 6.10 */

	aeps=1.e-13; reps=1.e-12;
	for(i1=0; i1<99; ++i1) {
		printf("type  xl = lower limit,   a1 = point at which integral is to be split \n");
		printf("                           (quits when a1<xl/2)\n");
		scanf(" %le %le", &xl, &a1);
		if(a1<xl/2) return 0;

		i=gaulag(&rint,xl,&a1,reps,aeps,&df,fun,f2,&nuse);
		printf(" ier = %d  no. of function evaluations = %d  xl = %e   a1 = %e \n", i,nuse,xl,a1);
		printf(" integral = %e   estimated error = %e  \n", rint, df);
	}
	return;
}

/*  The integrand  */

double fun(double x)

{
	return x/(1.+exp(x));
}

/*  The integrand*exp(x)  for lagure  */

double f2(double x)

{
/*	return 1./sqrt((x+1)*(b+x)); */
	return x/(1+exp(-x));
}

 

/*	To integrate a function over semi-infinite interval using a
	combination of Gauss-Legendre and Gauss-Laguerre formulas

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	A1 : (input/output) The point at which integral has to be broken
		A1 will be adjusted by the function.
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the integrand
	F2 : (input) Name of the function to calculate F(X)*EXP(X)
	NP : (output) Number of function evaluations used
		
	Error status is returned by the value of the function GAULAG.
		0 value implies successful execution
		30 implies specified accuracy was not achieved by GAUSS
		37 implies specified accuracy was not achieved by LAGURE
		38 implies specified accuracy was not achieved by both
			GAUSS and LAGURE
		In all cases DIF will contain the estimated accuracy

	Functions F(X) and F2(X) must be supplied by the user.

	Required functions : GAUSS, LAGURE, F, F2
*/

#include <math.h>

int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
	double *dif, int *npt, double (*fun) (double ));
int lagure(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt);

int gaulag(double *rint, double a, double *a1, double reps, double aeps,
	double *dif, double (*f) (double ), double (*f2) (double ), int *np)

{
	int npt,npt1,ier,ier1;
	double r1,r2,t1,dif1,a2, amax=50.0;

	r1=0.0; r2=0.0; *dif=0.0;
	*np=0; npt=0; ier1=0;
	if( (*a1)< a) *a1=a;

/*	To calculate integral over [A1,Infinity) */
	do {
		ier=lagure(&r1,*a1,aeps,reps,&dif1,f2,&npt1);
		*np=(*np)+npt1;
		if(ier == 0) break;

/*	If LAGURE fails then increase A1 */
		t1= (*a1);
		a2=(*a1)*2.; if(*a1+2. > a2) a2=(*a1)+2.;
		*a1=a2;
	} while (*a1 < amax);

	if(ier >0) {ier1=37; *a1=t1;}

/*	To calculate integral over [A,A1] */
	npt1=16;
	ier=0;
	if((*a1)-a > aeps) ier=gauss(&r2,a,*a1,&npt1,reps,aeps,dif,&npt,f);
	*rint=r1+r2;
	*dif=fabs(*dif)+fabs(dif1);
	*np=(*np)+npt;
	ier=ier+ier1;
	if(ier>ier1 && ier1>0) ier=38;
	return ier;
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



/*	To integrate a function over semi-infinite interval using Gauss-Laguerre formulas

	RINT : (output) Calculated value of the integral
	A : (input) The lower limit
	AEPS : (input) The required absolute accuracy
	REPS : (input) The required relative accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	F : (input) Name of the function to calculate the
		integrand (multiplied by EXP(X))
	NPT : (output) Number of function evaluations used
		
	Error status is returned by the value of the function LAGURE.
		0 value implies successful execution
		30 implies specified accuracy was not achieved
			DIF will contain the estimated accuracy

	Function F(X) must be supplied by the user

	Required functions : F
*/

#include <math.h>

int lagure(double *rint, double a, double aeps, double reps, double *dif,
	double (*f) (double ), int *npt)

{
	int i,j,n;
	double r1,r2,s1,expa;

/*	Weights and abscissas for Gauss-Laguerre quadrature
	W[N-2],...,W[2N-3], are the weights for N-point rule and
	X[N-2],...,X[2N-3], the corresponding abscissas
	Weights and abscissas are available for N=2,4,8,16,32
*/

	double x[62]={0.585786437627e0,     3.414213562373e0,
     .322547689619392312e0,1.74576110115834658e0,4.53662029692112798e0,
     9.39507091230113313e0,.170279632305101000e0,.903701776799379912e0,
     2.25108662986613069e0,4.26670017028765879e0,7.04590540239346570e0,
     10.7585160101809952e0,15.7406786412780046e0,22.8631317368892641e0,
     .087649410478927840e0,.462696328915080832e0,1.14105777483122686e0,
     2.12928364509838062e0,3.43708663389320665e0,5.07801861454976791e0,
     7.07033853504823413e0,9.43831433639193878e0,12.2142233688661587e0,
     15.4415273687816171e0,19.1801568567531349e0,23.5159056939919085e0,
     28.5787297428821404e0,34.5833987022866258e0,41.9404526476883326e0,
     51.7011603395433184e0,0.0444893658332670184e0,.234526109519618537e0,
     .576884629301886426e0,1.07244875381781763e0,1.72240877644464544e0,
     2.52833670642579488e0,3.49221327302199449e0,4.61645676974976739e0,
     5.90395850417424395e0,7.35812673318624111e0,8.98294092421259610e0,
     10.7830186325399721e0,12.7636979867427251e0,14.9311397555225573e0,
     17.2924543367153148e0,19.8558609403360547e0,22.6308890131967745e0,
     25.6286360224592478e0,28.8621018163234747e0,32.3466291539647370e0,
     36.1004948057519738e0,40.1457197715394415e0,44.5092079957549380e0,
     49.2243949873086392e0,54.3337213333969073e0,59.8925091621340182e0,
     65.9753772879350528e0,72.6876280906627086e0,80.1874469779135231e0,
     88.7353404178923987e0,98.8295428682839726e0,111.751398097937695e0};

	double w[62] = { 0.853553390593e0,     0.146446609407e0,
     .603154104341633602e0,.357418692437799687e0,.38887908515005384e-1,
     .53929470556132745e-3,.369188589341637530e0,.418786780814342956e0,
     .175794986637171806e0,.33343492261215651e-1,.27945362352256725e-2,
     .90765087733582131e-4,.84857467162725315e-6,.10480011748715104e-8,
     .206151714957800994e0,.331057854950884166e0,.265795777644214153e0,
     .136296934296377540e0,.47328928694125219e-1,.11299900080339453e-1,
     .18490709435263109e-2,.20427191530827846e-3,.14844586873981299e-4,
     .68283193308711996e-6,.18810248410796732e-7,.28623502429738816e-9,
     .2127079033224103e-11,.6297967002517868e-14,.5050473700035513e-17,
     .4161462370372855e-21,0.109218341952384971e0,0.210443107938813234e0,
     .235213229669848005e0,.195903335972881043e0,.129983786286071761e0,
     .70578623865717442e-1,.31760912509175070e-1,.11918214834838557e-1,
     .37388162946115248e-2,.98080330661495513e-3,.21486491880136419e-3,
     .39203419679879472e-4,.59345416128686329e-5,.74164045786675522e-6,
     .76045678791207815e-7,.63506022266258067e-8,.42813829710409289e-9,
     .2305899491891336e-10,.9799379288727094e-12,.3237801657729266e-13,
     .8171823443420719e-15,.1542133833393823e-16,.2119792290163619e-18,
     .2054429673788045e-20,.1346982586637395e-22,.5661294130397359e-25,
     .1418560545463037e-27,.1913375494454224e-30,.1192248760098222e-33,
     .2671511219240137e-37,.1338616942106256e-41,.4510536193898974e-47};


	expa=exp(-a);
/*	The 2-point formula */
	r1=(f(x[0]+a)*w[0]+f(x[1]+a)*w[1])*expa;
	*npt=2;
	n=2;

/*	Use higher order formula until convergence */
	for(j=2; j<=5; ++j) {
		n=n*2;
		r2=0.0;
		for(i=n-2; i<2*n-2; ++i) r2=r2+f(x[i]+a)*w[i];
		r2=r2*expa;
		*npt=(*npt)+n;
		*dif=r2-r1;
		*rint=r2;
		s1=reps*fabs(*rint); if(aeps>s1) s1=aeps;
		if(fabs(*dif) < s1) return 0;
		r1=r2;
	}

	return 30;
}
