/*	Multiple integration over a hyper-rectangle in n-dimensions
	using product Gauss-Legendre formulas

	A : (input) Array of length N containing the lower limit along each dimension
	B : (input) Array of length N containing the upper limit along each dimension
	N : (input) The number of dimensions
	M : (input/output) Integer array of length N specifying the formula
		to be used along each dimension. M[J]-point formula will
		be used along Jth dimension, M[J] should be 2,4,8,16 or 32
		otherwise it will be set to a default value of 2
	IND : (input/output) Integer array of length N specifying the number
		of subintervals to be used along each dimension. IND[J]>0
		otherwise it will be set to a default value of 1
	F : (input) Name of the function to calculate the integrand
		Function(N,X) should calculate the integrand, where N is the
		number of dimensions and X is an array of length N containing
		the coordinates of the point where integrand is to be calculated
	RINT : (output) The calculated value of the integral
	REPS : (input) The required relative accuracy
	AEPS : (input) The required absolute accuracy
		The estimated error should be less than max(AEPS,REPS*fabs(RINT))
	DIF : (output) estimated (absolute) error achieved by the function
	NUM : (output) Number of function evaluations used
	MAXPT : (input/output) Maximum number of function evaluations permitted
		If MAXPT <1 it is set to a default value of MAXPTS (=1100000)
		
	Error status is returned by the value of the function MULINT.
		0 value implies successful execution
		39 implies specified accuracy was not achieved in
			which case DIF will contain the estimated accuracy
		305 implies N<1 and no calculations are done
		307 implies that number of points exceeded MAXPT in
			first attempt and no approximation of RINT is calculated

	Function F(N,X) must be supplied by the user
	
	Required functions : NGAUSS, F
*/

#include <math.h>

int ngauss(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *ri, int *num, int *maxpt);

int mulint(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
	double *rint, double reps, double aeps, double *dif, int *num, int *maxpt)

{
	int i,j,i1,m1,no,qc,ier, maxpts=1100000;
	double rint1,r1;

/*	Set M[I] and IND[I] to default values if they are unacceptable */
	for(i=0; i<n; ++i) {
		if((m[i]!=2) && (m[i]!=4) && (m[i]!=8) && (m[i]!=16) && (m[i]!=32)) m[i]=2;
		if(ind[i] <1) ind[i]=1;
	}
	if(*maxpt <1) *maxpt=maxpts;
	*num=0;

/*	Evaluate the integral */
	ier=ngauss(a,b,n,m,ind,f,rint,&no,maxpt);
	*num=(*num)+no;
	if(ier > 100) return ier;

/*	Iteration to check and improve the accuracy of integral */
	for(i=0; i<10; ++i) {
		qc=1;
		*dif=0.0;

/*	Check for convergence along each dimension by doubling the
	number of points at which function is evaluated */
		for(j=0; j<n; ++j) {
			if(*num > *maxpt) return 39;
			m1=m[j]; i1=ind[j];
/*	If M[J]<32 then double the order of formula */
			if(m[j]<32) m[j]=2*m[j];

/*	otherwise double the number of subintervals */
			else ind[j]=2*ind[j];

			ier=ngauss(a,b,n,m,ind,f,&rint1,&no,maxpt);
			*num=(*num)+no;
			*dif=(*dif)+fabs(*rint-rint1);
			if(ier > 100) return 39;

			r1=reps*fabs(rint1); if(aeps>r1) r1=aeps;
			if(fabs(rint1-*rint) < r1) {m[j]=m1; ind[j]=i1;}
			else {*rint=rint1; qc=0;}
		}

/*	If satisfactory accuracy is achieved for all dimensions then return */
		if(qc==1) return 0;

/*	If the number of function evaluations exceeds maxpt then quit */
		if(*num> (*maxpt)) return 39;
	}
	return 39;
}
