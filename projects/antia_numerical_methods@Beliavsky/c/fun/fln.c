/*	Function to calculate the function value as required for
	line search without derivatives

	FCN : (input) Name of function to calculate the required function
	X : (input) Parameter along the line to specifying the point where
		function evaluation is required
	V : (input/output) Array of length 2N, first N elements specify the
		direction of line search. After execution next N elements will contain
		the coordinates of the point at which function is evaluated.
	X0 : (input) Array of length N, containing the coordinates
		of the starting point for line search
	N : (input) Number of variables in the function to be minimised
	NUM : (input/output) Integer variable to keep count of function evaluations

	Function FCN(N,X,FX) to calculate the required function, must be supplied
		by the user. Here N is the number of variables, FX is the
		function value at X. X is an array of length N.

	Required functions : FCN
*/

double fln(void fcn(int , double * , double * ), double x, double v[],
	double x0[], int n, int *num)

{
	int i;
	double fl;

	*num=(*num)+1;
/*	coordinates of the required points */
	for(i=0; i<n; ++i) x0[n+i]=x0[i]+v[i]*x;
	fcn(n,&x0[n],&fl);
	return fl;
}
