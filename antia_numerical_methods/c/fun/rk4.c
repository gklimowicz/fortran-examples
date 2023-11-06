/*	To perform one step of integration of ordinary differential equations
	using a fourth-order Runge Kutta method 

	N : (input) Number of first order differential equations to be solved
	T : (input) Initial value of independent variable t, at
		which initial values are specified. This value is not updated.
	Y0 : (input) Array of length N containing the initial
		values of variables.
	DY0 : (input) Array of length N containing the derivatives
		of Y at the initial point Y0
	H : (input) The step size to be used for integration
	Y1 : (output) Array of length N containing the solution at t=T+H
	DIF : (input) Name of function to calculate the right hand side
		of differential equation y'=f(t,y)

	Function DIF(T,N,Y,DY) must be supplied by the user to specify
	the differential equation. T is the value of independent variable,
	N is the number of variables, Y is an array of length N containing
	the values of variables. DY is an array of length N which should
	contain the calculated values of derivatives at (T,Y).

	Required functions : DIF
*/	

#include <math.h>
#include <stdlib.h>

void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
	void dif(double , int , double * , double * ))

{
	int i;
	double h2,t1;
	double *wk1, *wk2;

	h2=h/2.;
	wk1=(double *) calloc((size_t) n, sizeof(double));
	wk2=(double *) calloc((size_t) n, sizeof(double));
	for(i=0; i<n; ++i) {
		y1[i]=h2*dy0[i];
/*	y_0+0.5k_1 */
		wk1[i]=y0[i]+y1[i];
	}

	t1=t+h2;
	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i) {
		y1[i]=y1[i]+h*wk2[i];
/*	y_0+0.5k_2 */
		wk1[i]=y0[i]+h2*wk2[i];
	}

	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i) {
		y1[i]=y1[i]+h*wk2[i];
/*	y_0+k_3 */
		wk1[i]=y0[i]+h*wk2[i];
	}

	t1=t+h;
	dif(t1,n,wk1,wk2);
	for(i=0; i<n; ++i)  y1[i]=y0[i]+(y1[i]+h2*wk2[i])/3.0;
	free(wk2); free(wk1);
	return;
}
	
