/*	To solve linear Volterra equation using Trapezoidal rule */

#include <stdio.h>
#include <math.h>

int volt(int n, double a, double h, double f[], double x[],
	double fg(double ), double fker(double , double ), int it);

double fg(double x);
double fker(double x, double t);

double al;

main()
{
	int i,i1,j,n,m, id, iflg, ier,np,nmax;
	double hh, x[500], f[500],tn,erc,reps,t0,fi;

/*	Example 13.7 */

	t0=0.0;

/*	it=-1 for smoothed solution, for it=1 no smoothing is applied */
	for(i1=0; i1<99; ++i1) {
		printf("type n=no. of steps,  it=+1/-1,   hh=step size \n");
		printf("                       (quits when n<=0)\n");
		scanf(" %d %d %le",&n,&id,&hh);
		if(n<=0) return 0;

		i=volt(n,t0,hh,f,x,fg,fker,id);
		printf(" ier = %d    no. of steps =  %d    it = %d    step size = %e\n",i,n,id,hh);
		printf("   x       solution\n");
		for(i=0; i<n; i=i+5) printf(" %e %e \n",x[i],f[i]);

	}
	return;
}


/*	The right hand side function */

double fg(double x)

{
	return sinh(x);
}


/*	The kernel */

double fker(double x, double t)

{
		return exp(x-t);
}



/*	To solve linear Volterra equations using quadrature method with the
	trapezoidal rule

	N : (input) Number of points at which the solution is required
	A : (input) Lower limit of the integral
	H : (input) Uniform spacing to be used between abscissas
	F : (output) Array of length N containing the calculated solution
		F[I] is the value of solution at X[I].
	X : (output) Array of length N containing the abscissas used
		in the trapezoidal rule. The spacing is assumed to be uniform
	FG : (input) Name of the function used to calculate the right
		hand side g(X). 
	FKER : (input) Name of the function used to calculate the kernel K(x,t)
	IT : (input) Integer variable to specify the type of integral equation
		If IT=1 Volterra equation of the first kind is solved
		If IT=-1 Volterra equation of the first kind is solved and
			computed values are smoothed as explained in Section 12.8
		If IT=2 Volterra equation of the second kind is solved
		
	Error status is returned by the value of the function VOLT.
		0 value implies successful execution
		712 implies N<3, in which case no calculations are done.
		751 implies that denominator is zero at some stage
			No further calculations are done

	Functions FG(X) and FKER(X,T) must be supplied by the user.
		FG is the right hand side function g(x) and FKER(X,T) is 
		the kernel. 

	Required functions : FG, FKER
*/

#include <math.h>


int volt(int n, double a, double h, double f[], double x[],
	double fg(double ), double fker(double , double ), int it)

{
	int i,j;
	double di,fi,fi1;

	if(n<3) return 712;

	x[0]=a;
	if(it==2) f[0]=-fg(a);	 /*	Starting value for equations of second kind */
	else {
		di=h*fker(a,a);
		if(di==0.0) return 751;
/*	Starting value for equations of first kind */
		f[0]=(fg(a+h)-fg(a-h))/(2.*di);
	}

/*	Continue the solution using the trapezoidal rule */
	for(i=1; i<n; ++i) {
		x[i]=a+i*h;
		fi=0.5*fker(x[i],a)*f[0];
		for(j=1; j<=i-1; ++j) fi=fi+fker(x[i],x[j])*f[j];
		fi=fi*h-fg(x[i]);

		di=-h*fker(x[i],x[i])/2.0;
		if(it==2) di=di+1.0;
		if(di==0.0) return 751;
		f[i]=fi/di;
	}

	if(it==-1) {
/*	Apply smoothing for F(2),...,F(N-1) */
		fi=0.25*(f[0]+2.*f[1]+f[2]);
		for(i=1; i<n-2; ++i) {
			fi1=0.25*(f[i]+2.*f[i+1]+f[i+2]);
			f[i]=fi;
			fi=fi1;
		}
		f[n-2]=fi;
	}
	return 0;
}
