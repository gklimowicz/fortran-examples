/*	To evaluate a polynomial and its derivatives at any point

	N : (input) Degree of the polynomial
	A : (input) Array of length N+1 containing the coefficients
		of the polynomial. A[0] is the constant term and A[N] is
		the coefficient of X**N
	X : (input) The value of x at which the polynomial is to be evaluated
	ND : (input) Number of derivatives required. The function always
		calculates the first derivative irrespective of value of ND
	PD : (output) Array of length PD containing the calculated
		values of the derivatives. PD[I-1] will contain the Ith derivative
	The polynomial value is returned through POLD

	Required functions : None
*/

double pold(int n, double a[], double x, int nd, double pd[])

{
	int i,j;
	double pf;

	pf=a[n];
	for(i=0; i<nd; ++i) pd[i]=0.0;

	for(j=n-1; j>=0; --j) {

		for(i=nd-1; i>=1; --i) pd[i]=pd[i]*x+(i+1)*pd[i-1];
		pd[0]=pd[0]*x+pf;
		pf=pf*x+a[j];
	}
	return pf;
}
