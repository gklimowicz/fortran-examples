/*	To evaluate a rational function at any point

	M : (input) Degree of numerator
	K : (input) Degree of denominator
	A : (input) Array of length M+1 containing the coefficients
		of the polynomial in numerator. A[0] is the constant term
		while A[M] is the coefficient of X**M.
	B : (input) Array of length K+1 containing the coefficients
		of the polynomial in denominator. B[0] is the constant term
		while B[K] is the coefficient of X**K.
	X : (input) The value of X at which the rational function needs
		to be evaluated
	The value of rational function will be returned as RMK

	Required functions : None
*/

double rmk(int m, int k, double a[], double b[], double x)

{
	int i;
	double den, rnu;

	rnu=a[m];
	for(i=m-1; i>=0; --i) rnu=rnu*x+a[i];

	den=b[k];
	for(i=k-1; i>=0; --i) den=den*x+b[i];

	return rnu/den;
}
