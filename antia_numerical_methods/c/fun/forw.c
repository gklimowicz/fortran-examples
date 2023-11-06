/*	To solve the forward problem corresponding to a linear inverse
	problem. This function may be used to generate artificial data
	for testing inversion techniques.

	NP : (input) Number of points used in defining the kernels.
	NM : (input) Number of data points in the inverse problem
		which should be same as the number of kernels that are
		supplied in array RKER.
	R : (input) Array of length NP containing the coordinates
		of points at which kernels are available.
	RKER : (input) Array of length IK*NP containing the kernels
		for the inverse problem. RKER[j][i] should contain the
		value at R[j] for the ith kernel.
	IK : (input) The second dimension of RKER as declared in the calling
		function
	DI : (output) Array of length NM containing the calculated
		data points using the kernel.
	F : (input/output) Array of length NP containing the function
		value at points in R. F[I] should contain the function
		value at R[I]. If IFLG=0, the function values are
		calculated using user supplied function FUN,
		otherwise, these values must be supplied while calling
		the function.
	FUN : (input) Name of function routine to calculate the given
		function. This is used only if IFLG=0, otherwise the
		function values are to be supplied in array F.
	IFLG : (input/output) Integer parameter used as a flag to decide
		the type of computation required.
		If IFLG=0, then the function values are calculated using
			a user supplied function FUN. These values are stored
			in array F and IFLG is set to 1 so that next time
			the values need not be calculated.
		For other values of IFLG the function values must be
			supplied in array F.
		
	Error status is returned by the value of the function FORW.
		0 value implies successful execution
		711 implies that IK<NM and no calculations are done

	Function FUN(X) must be supplied by the user

	Required functions : FUN
*/

#include <math.h>

int forw(int np, int nm, double r[], double *rker, int ik, double di[],
	double f[], double fun(double ), int *iflg)

{
	int i,ir;
	double h,s1;

	if(ik<nm) return 711;

	if(*iflg==0) {
/*     Calculate the function value using supplied routine */
		for(i=0; i<np; ++i) f[i]=fun(r[i]);
		*iflg=1;
	}
 
/*     Calculate the integrals */
	for(i=0; i<nm; ++i) {
		s1=0.0;
		h=(r[1]-r[0])/2.0;
		for(ir=0; ir<np; ++ir) {
			s1=s1+h*f[ir]*rker[i+ir*ik];
			h=(r[ir+2]-r[ir])/2.0;
			if(ir==np-2) h=(r[ir+1]-r[ir])/2.0;
		}
		di[i]=s1;
	}
	return 0;
}
