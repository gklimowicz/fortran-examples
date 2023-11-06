/*	Function to transform from hyper-spherical coordinates to
	Cartesian coordinates.
	It can be used with any of the functions for multiple integration
	for integration over hyper-spherical shells

	N : (input) Number of dimensions
	X : (input) Array of length N, giving hyper-spherical coordinates 
		First coordinate is radial distance r
	Y : (output) Array of length N giving the transformed Cartesian
		coordinates, which is passed on to the Function FUNSPH(N,Y)
		for calculating the required function.

	The Function FUNSPH(N,Y) must be supplied by the user. The name of
	this function is not passed as argument and hence it has to have
	the same name as occurring here.
		
	Required functions : FUNSPH
*/

#include <math.h>
#include <stdlib.h>

double funsph(int n, double y[]);

double sphnd(int n, double x[])

{
	int i;
	double t1,p1,s1;
	double *y;

/*	The Cartesian coordinates are given by
	y(1)=x(1)*COS(x(2))
	y(2)=x(1)*SIN(x(2))*COS(x(3))
	y(3)=x(1)*SIN(x(2))*SIN(x(3))*COS(x(4))
	.........
	y(n-1)=x(1)*SIN(x(2))*SIN(x(3))*...*SIN(x(n-1))*COS(x(n))
	y(n)=x(1)*SIN(x(2))*SIN(x(3))*...*SIN(x(n))

	P1 is the volume element in hyper-spherical coordinates
*/
	y=(double *) calloc((size_t) n, sizeof(double));
	t1=x[0];
	p1=1.0;
	for(i=0; i<n-1; ++i) {
		y[i]=t1*cos(x[i+1]);
		p1=p1*t1;
		t1=t1*sin(x[i+1]);
	}
	y[n-1]=t1;

/*	FUNSPH should calculate the required integrand using Cartesian coordinates */
	s1=p1*funsph(n,y);
	free(y);
	return s1;
}
