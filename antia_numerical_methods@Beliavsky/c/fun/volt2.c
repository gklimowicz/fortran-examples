/*	To solve nonlinear Volterra equation of the second kind using
	quadrature method with Simpson's rule

	Integral[K(x,t,f(x)) dt] over [0,x] = f(x)+g(x)

	N : (input) Number of points at which the solution is required
	A : (input) Lower limit of the integral
	H : (input) Uniform spacing to be used between abscissas
	F : (output) Array of length N containing the calculated solution
		F[I] is the value of solution at X[I].
	X : (output) Array of length N containing the abscissas used
		in the trapezoidal rule. The spacing is assumed to be uniform
	FG : (input) Name of the function used to calculate the right
		hand side g(X). 
	FKER : (input) Name of the function used to calculate the
		kernel K(x,t,f), note that in this case the "kernel" includes
		the unknown solution also, because the equation is nonlinear.
	REPS : (input) The required accuracy to which nonlinear equation is
		solved using fixed point iteration. It will not control the
		truncation error which is determined by the mesh spacing.
		
	Error status is returned by the value of the function VOLT2.
		0 value implies successful execution
		712 implies N<3, in which case no calculations are done.
		752 implies that fixed point iteration failed at some stage.
			No further calculations are done

	Functions FG(X) and FKER(X,T,F) must be supplied by the user.
		FG is the right hand side function g(x) and FKER(X,T,F) is 
		the kernel (the integrand in integral equation). 

	Required functions : FG, FKER
*/

#include <math.h>

int volt2(int n, double a, double h, double f[], double x[],
	double fg(double ), double fker(double , double , double), double reps)

{
	int i,j, nit=30;
	double g2,fk1,f21,f22,f23,f24,fi,f0,f1,r1,r2, eta=1.0e-30;

	if(n<3) return 712;

/*	Generating the starting values, F[0] and F[1] */
	x[0]=a;
	f[0]=-fg(a);
	g2=fg(a+h);
	fk1=fker(a+h,a,f[0]);
	f21=-g2+h*fk1;
	f22=-g2+0.5*h*(fk1+fker(a+h,a+h,f21));
	f23=0.5*(f[0]+f22);
	f24=-fg(a+0.5*h)+0.25*h*(fker(a+0.5*h,a,f[0])+fker(a+0.5*h,a+0.5*h,f23));
	x[1]=a+h;
	f[1]=-g2+h*(fk1+4.*fker(a+h,a+0.5*h,f24)+fker(a+h,a+h,f22))/6.0;

/*	Continuing the solution */
	for(i=2; i<n; i +=2) {
		x[i]=a+i*h;
/*	Use Simpson's 1/3 rule for odd I */
		fi=fker(x[i],a,f[0]);
		for(j=1; j<=i-2; j +=2) {
			fi=fi+4.*fker(x[i],x[j],f[j]);
			fi=fi+2.*fker(x[i],x[j+1],f[j+1]);
		}
		fi=(fi+4.*fker(x[i],x[i-1],f[i-1]))*h/3.0-fg(x[i]);

/*	The predictor */
		f0=2.0*f[i-1]-f[i-2];
/*	Iteration on the corrector */
		for(j=1; j<=nit; ++j) {
			f1=fi+h*fker(x[i],x[i],f0)/3.;
			r2=fabs(f1-f0)/(fabs(f1)+fabs(f1-f[i-1])+eta);
			if(r2<reps) break;
			f0=f1;
		}
		if(j>nit) return 752;

		f[i]=f1;
		if(i>=n-1) return 0;

/*	For even I+1 use Simpson's 1/3 and 3/8 rule */
		x[i+1]=a+(i+1)*h;
		fi=fker(x[i+1],a,f[0]);
/*	Use 1/3 rule on the first I-3 points */
		for(j=1; j<=i-4; j +=2) {
			fi=fi+4.*fker(x[i+1],x[j],f[j]);
			fi=fi+2.*fker(x[i+1],x[j+1],f[j+1]);
		}
		if(i>2) {
			fi=(fi+4.0*fker(x[i+1],x[i-3],f[i-3]))/3.0;
			fi=fi+17.0*fker(x[i+1],x[i-2],f[i-2])/24.0;
		}
		else fi=3.0*fi/8.0;		 /*	If I+1=4 use only 3/8 rule */
		fi=h*(fi+9.*(fker(x[i+1],x[i-1],f[i-1])+fker(x[i+1],x[i],f[i]))/8.)-fg(x[i+1]);

/*	The predictor */
		f0=2.0*f[i]-f[i-1];
/*	Iteration on the corrector */
		for(j=1; j<=nit; ++j) {
			f1=fi+3.*h*fker(x[i+1],x[i+1],f0)/8.0;
			r2=fabs(f1-f0)/(fabs(f1)+fabs(f1-f[i])+eta);
			if(r2<reps) break;
			f0=f1;
		}
		if(j>nit) return 752;

		f[i+1]=f1;
	}
	return 0;
}
