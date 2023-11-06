/*	To calculate coefficients of cubic spline interpolation with
		not-a-knot boundary conditions

	X : (input) Array of length N containing x values
	F : (input) Array of length N containing values of function at X[I]
		F[I] is the tabulated function value at X[I].
	N : (input) Length of table X, F
	C : (output) Array of length 3*N containing the spline coefficients
		
	Error status is returned by the value of the function SPLINE.
		0 value implies successful execution
		201 implies that N<2

	Required functions : None
*/

#include <math.h>

double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);

int spline(double x[], double f[], int n, double c[][3])

{
	int i,j;
	double g, c1, cn, div12, div01;

	if(n<2) return 201;
	else if(n == 2) {
/*	Use linear interpolation */
		c[0][0]=(f[1]-f[0])/(x[1]-x[0]);
		c[0][1]=0.0;
		c[0][2]=0.0;
		return 0;
	}
	else if(n == 3) {
/*	Use quadratic interpolation */
		div01=(f[1]-f[0])/(x[1]-x[0]);
		div12=(f[2]-f[1])/(x[2]-x[1]);
		c[0][2]=0.0;
		c[1][2]=0.0;
		c[0][1]=(div12-div01)/(x[2]-x[0]);
		c[1][1]=c[0][1];
		c[0][0]=div01+c[0][1]*(x[0]-x[1]);
		c[1][0]=div12+c[0][1]*(x[1]-x[2]);
	        return 0;
	}
	else {
/*	Use cubic splines 

	Setting up the coefficients of tridiagonal matrix */
		c[n-1][2]=(f[n-1]-f[n-2])/(x[n-1]-x[n-2]);
		for(i=n-2; i>=1; --i) {
			c[i][2]=(f[i]-f[i-1])/(x[i]-x[i-1]);
			c[i][1]=2.*(x[i+1]-x[i-1]);
/*	The right hand sides */
			c[i][0]=3.*(c[i][2]*(x[i+1]-x[i])+c[i+1][2]*(x[i]-x[i-1]));
		}

/*	The not-a-knot boundary conditions */
		c1=x[2]-x[0];
		c[0][1]=x[2]-x[1];
		c[0][0]=c[1][2]*c[0][1]*(2.*c1+x[1]-x[0])+c[2][2]*(x[1]-x[0])*(x[1]-x[0]);
		c[0][0]=c[0][0]/c1;
		cn=x[n-1]-x[n-3];
		c[n-1][1]=x[n-2]-x[n-3];
		c[n-1][0]=c[n-1][2]*c[n-1][1]*(2.*cn+x[n-1]-x[n-2]);
		c[n-1][0]=(c[n-1][0]+c[n-2][2]*(x[n-1]-x[n-2])*(x[n-1]-x[n-2]))/cn;
/*	Solving the equation by Gaussian elimination */
		g=(x[2]-x[1])/c[0][1];
		c[1][1]=c[1][1]-g*c1;
		c[1][0]=c[1][0]-g*c[0][0];
		for(j=1; j<n-2; ++j) {
			g=(x[j+2]-x[j+1])/c[j][1];
			c[j+1][1]=c[j+1][1]-g*(x[j]-x[j-1]);
			c[j+1][0]=c[j+1][0]-g*c[j][0];
		}
		g=cn/c[n-2][1];
		c[n-1][1]=c[n-1][1]-g*(x[n-2]-x[n-3]);
		c[n-1][0]=c[n-1][0]-g*c[n-2][0];


/*	The back-substitution */
		c[n-1][0]=c[n-1][0]/c[n-1][1];
		for(i=n-2; i>=1; --i) c[i][0]=(c[i][0]-c[i+1][0]*(x[i]-x[i-1]))/c[i][1];
		c[0][0]=(c[0][0]-c[1][0]*c1)/c[0][1];

/*	Calculating the coefficients of cubic spline */
		for(i=0; i<n-1; ++i) {
			c[i][1]=(3.*c[i+1][2]-2.*c[i][0]-c[i+1][0])/(x[i+1]-x[i]);
			c[i][2]=(c[i][0]+c[i+1][0]-2.*c[i+1][2])/((x[i+1]-x[i])*(x[i+1]-x[i]));
		}
/*	Set the coefficients for interval beyond X(N) using continuity
	of second derivative, although they may not be used. */
		c[n-1][1]=c[n-1][1]+3*(x[n-1]-x[n-2])*c[n-2][2];
		c[n-1][2]=0.0;
		return 0;
	}
}
