/*	To calculate Legendre polynomial P_l(X)

	L : (input) Order of polynomial (L>=0)
	X : (input) Argument at which the value of polynomial is required
	P : (output) Array of length L+1, which will contain the
		calculated values of polynomials. P[j] will contain P_j(X)

	Required functions : None
*/

void pleg(int l, double x, double p[])

{
	int n;
	
	if(l<0) return;
	p[0]=1.0;
	p[1]=x;
	for(n=2; n<=l; ++n) p[n]=((2*n-1)*x*p[n-1]-(n-1)*p[n-2])/n;
	return;
}
