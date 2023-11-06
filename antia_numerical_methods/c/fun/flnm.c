/*	Function to calculate the function value and derivative
	as required for line search

	FCN : (input) Name of function to calculate the required function
	X : (input) Parameter along the line to specify the point where
		function evaluation is required
	DF : (output) First derivative of function along the line at X
	V : (input/output) Array of length 3N, first N elements specify the
		direction of line search. Next N elements will contain the
		coordinates of the point at which function is evaluated,
		while the last N elements contain the gradient vector at the point
	X0 : (input) Array of length N, containing the coordinates
		of starting point for line search
	N : (input) Number of variables in the function to be minimised
	NUM : (input/output) Integer variable to keep count of function evaluations

	Function FCN(N,X,FX,G) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X and G is the gradient vector. X and G
		are arrays of length N.

	Required functions : FCN
*/

double flnm(void (*fcn) (int , double * , double * , double * ), double x,
	double *df, double v[], double x0[], int n, int *num)

{
	int i,n2;
	double fl;

	*num=(*num)+1;
	n2=2*n;

/*	The coordinates of the required point */
	for(i=0; i<n; ++i) v[n+i]=x0[i]+v[i]*x;

	fcn(n,&v[n],&fl,&v[n2]);

/*	The first derivative along the search direction */
	*df=0.0;
	for(i=0; i<n; ++i) *df=(*df)+v[i]*v[n2+i];
	return fl;
}
