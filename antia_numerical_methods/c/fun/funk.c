/*	Function to calculate the integrand for calculating
	PSI(I,X) as required by function FREDCO.

	Function FKER(X,T) is the kernel K(x,t) and PHI(I,T) is the Ith
	basis function, phi_i(t). The argument X and I are passed through
	global variables XFRED and IFRED respectively.

	Functions FKER(X,T) and PHI(I,T) must be supplied by the user

	Required functions : FKER, PHI
*/

/*	To pass parameters from function FREDCO put these statements before fredco */
double XFRED;
int IFRED;

double funk(double t)

{
	return fker(XFRED,t)*phi(IFRED,t);
}
