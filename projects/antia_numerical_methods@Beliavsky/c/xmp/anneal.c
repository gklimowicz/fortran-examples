/*	Minimisation using simulated annealing */

#include <stdio.h>
#include <math.h>

FILE *FP;

void change(double temp, int *idum, int ip[]);
double fcn(int ip[]);
double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier);

int spline(double x[], double f[], int n, double c[][3]);
double ranf(int *iseed);

double ERR[600],X0[600],F0[600],X[50];

main()
{
	int i,j,id,idum,nch,ich,it,ic,n, ip[20],ip1[20],ip0[20];
	double res,resmin,res0,temp,temp0,expi;

/*	Example 8.7 */

/*     The number of knots to be fixed */
	n=12;
	resmin=1.0;

/*	The output file: The first 50 lines are the input data
	The next 13 lines contain the selected knots
	The last 501 lines contain the interpolated values using selected knots */

	FP=fopen("titan.out","w");

/*   Try 2 different sequence of random numbers, -ID is the seed
     Only the last solution is written out in the file */

      for(id=1; id<3; ++id) {

/*     Initial distribution of knots */
		  ip[0]=1;
		  ip[n-1]=48;
		  for(i=1; i<n-1; ++i) {ip[i]=4*i; ip1[i]=ip[i];}
		  res0=fcn(ip);
		  resmin=res0;
		  temp=res0; if(temp<1.0) temp=1.0;
		  temp0=temp;
		  nch=500;
		  idum=-id;
		
/*       The loop for annealing  */
		  for(it=0; it<400; ++it) {
			  ich=0;

/*	      Perform 10000 trials at each temperature */
			  for(ic=0; ic<10000; ++ic) {
				  for(i=0; i<n; ++i) ip0[i]=ip[i];
				  change(temp,&idum,ip);
				  res=fcn(ip);
				  if(res<resmin) {
/*			Store the minimum value */
					  for(i=0; i<n; ++i) ip1[i]=ip[i];
					  resmin=res;
				  }
				  expi=2.0;
				  if(fabs((res0-res)/temp)<500.0) expi=exp((res0-res)/temp);
				  else if(res>res0) expi=0.0;

				  if(res<res0 || ranf(&idum)<expi) {
/*			Accept the new point */
					  res0=res;
					  ich=ich+1;
				  }
				  else {
					  for(i=0; i<n; ++i) ip[i]=ip0[i];
				  }
				  if(ich>=nch) break;
			  }
/*         Reduce the temperature slowly */
			  temp=temp*0.98;

			  if(res0-resmin > 5*temp || ich==0) {
				  for(i=0; i<n; ++i) ip[i]=ip1[i];
				  res0=resmin;
			  }
 
/*         The convergence criterion */
			  if(ich==0 && temp<0.04*res0) break;
			  if(temp<1.e-5*res0) break;
		  }

/*       Write out the solution */
		  res=fcn(ip1);
		  printf("random no. seed = %d    function value at minimum = %e\n",-id,res);
		  printf(" serial no. of knots for minimum error :\n");
		  for(i=1; i<n-1; ++i) printf(" %d ",ip1[i]);
		  printf("\n");
	  }

/*     write the solution in file titan.out */
	  fprintf(FP,"#   the selected knots :\n");
	  for(i=0; i<n; ++i) fprintf(FP," %f \n",X[ip1[i]]);
	  fprintf(FP,"#    X     interpolated value   error \n");
	  for(i=0; i<501; ++i) fprintf(FP," %e %e %e \n",X0[i],F0[i],ERR[i]);
	  return;
}
 
 

 
/*	To evaluate the cubic spline interpolant at a specified point

	XB : (input) point at which interpolation is required
	N : (input) Number of points in the table
	X : (input) Array of length N, containing the abscissas
	F : (input) Array of length N, containing the function values at X[I]
	C : (input) Array of length 3*N containing the spline coefficients
		which should have been calculated using SPLINE
	DFB : (output) First derivative of spline at x=XB
	DDFB : (output) Second derivative of spline at x=XB
	IER : (output) error parameter, IER=0 if execution is successful
		IER=24 implies XB is outside the range of table on higher side
		IER=25 implies XB is outside the range of table on lower side
		IER=201 implies N<2
	SPLEVL will be the interpolated value at x=XB

	Required functions : None
*/

#include <math.h>

double splevl(double xb, int n, double x[], double f[], double c[][3],
		double *dfb, double *ddfb, int *ier)

{
	int i,j,igh,nigh, mid;
	double dx, r1;
	static int low=-1;
	
	if(n<2) {*ier=201; return 0.0;}

	*ier=0;
/*	If the previous value of LOW is inadmissible, set the range to (0,N-1) */
	if(low<0 || low>=n-1) {low=0; igh=n-1;}
	else igh=low+1;

	while((xb<x[low] && xb<x[igh]) || (xb>x[low] && xb>x[igh])) {
/*	Extend the range */
		if((xb>x[low]) == (x[n-1]>x[0])) {
/*	Extend the range on higher side */
			if(igh >= n-1) {*ier=24; low=n-2; break;}
			else {
				nigh=igh+2*(igh-low); if(n-1 < nigh) nigh=n-1;
				low=igh; igh=nigh;
			}
		}

		else {
/*	Extend the range on lower side */
			if(low <= 0) {*ier=25; igh=low+1; break;}
			else {
				nigh=low;
				low=low-2*(igh-low); if(low<0) low=0;
				igh=nigh;
			}
		}
	}


/*	Once the point is bracketed between two tabular points locate it by bisection */
	while((igh-low > 1) && (xb != x[low])) {
		mid=(low+igh)/2;
		if((xb<= x[mid]) == (xb<= x[low])) low=mid;
		else igh=mid;
	}

	dx=xb-x[low];
	r1=((c[low][2]*dx+c[low][1])*dx+c[low][0])*dx+f[low];
	*dfb=(3.0*c[low][2]*dx+2.*c[low][1])*dx+c[low][0];
	*ddfb=6.*c[low][2]*dx+2.*c[low][1];
	return r1;
}


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



/*	To generate uniformly distributed random numbers in interval (0,1)

	ISEED : (input/output) is an integer value used as the seed
		It should be initialised to negative value before first call
		and should not be modified between successive calls.

	Required functions : None
*/

#include <math.h>


double ranf(int *iseed)

{
	int m1=714025, ia1=1366, ic1=150889;
	int m2=214326, ia2=3613, ic2=45289;
	int m3=139968, ia3=3877, ic3=29573;
	int i,j, ish=43;
	double r1;

	static double rm1,rm2,ran[43];
	static int iflg=0, is1,is2,is3;
	
/*	Initialise on first call or when ISEED<0 */
	if(*iseed < 0 || iflg == 0) {
		iflg=1;
		rm1=1.0/m1;
		rm2=1.0/m2;

/*	Seeds for the three random number generators */
		is1=-(*iseed); is1=is1-m1*(is1/m1);
		is2=ia1*is1+ic1; is2=is2-m1*(is2/m1);
		is3=ia2*is2+ic2; is3=is3-m2*(is3/m2);
		*iseed=1;

/*	Store ISH random numbers in the array RAN */
		for(j=0; j<ish; ++j) {
			is1=ia1*is1+ic1; is1=is1-m1*(is1/m1);
			is2=ia2*is2+ic2; is2=is2-m2*(is2/m2);
			ran[j]=(is1+is2*rm2)*rm1;
		}
	}

	is1=ia1*is1+ic1; is1=is1-m1*(is1/m1);
	is2=ia2*is2+ic2; is2=is2-m2*(is2/m2);
	is3=ia3*is3+ic3; is3=is3-m3*(is3/m3);

/*	Select a random entry from RAN and store a new number in its place */
	i=(ish*is3)/m3;
	r1=ran[i];
	ran[i]=(is1+is2*rm2)*rm1;
	return r1;
}

 

double fcn(int ip[])

{
	static int i,j,np,ipas=-999,n0,n,ier;
	static double h,xb,dfb,ddfb,errmax,fb,df, c[50][3],xp[20],fp[20];

/*     To calculate the function to be minimised.
     The function value is the maximum error in spline approximation
     using the selected knots.
     IP is an array specifying the serial no. of knots to be used */

/*	The input data for the Titanium problem */
	double f[49] = {0.644,0.622,0.638,0.649,0.652,0.639,0.646,
            0.657,0.652,0.655,0.644,0.663,0.663,0.668,
            0.676,0.676,0.686,0.679,0.678,0.683,0.694,
            0.699,0.710,0.730,0.763,0.812,0.907,1.044,
            1.336,1.881,2.169,2.075,1.598,1.211,0.916,
            0.746,0.672,0.627,0.615,0.607,0.606,0.609,
            0.603,0.601,0.603,0.601,0.611,0.601,0.608};
 
 
	if(ipas==-999) {
/*	During the first call write out the input data in file titan.out */
		np=49;
		fprintf(FP,"#     X        f(X)\n");
		for(i=0; i<np; ++i) {
			X[i]=585+10*i;
			fprintf(FP," %e %e\n",X[i],f[i]);
		}

/*	Calculate the cubic spline interpolation at 501 points using full table */
		i=spline(X,f,np,c);
		n0=501;
		h=(X[np-1]-X[0])/(n0-1);
		for(i=0; i<n0; ++i) {
			xb=X[0]+h*i;
			X0[i]=xb;
			F0[i]=splevl(xb,np,X,f,c,&dfb,&ddfb,&ier);
		}
		ipas=1;
	}

	n=12;
/*	Setup the knots using array IP */
	for(i=0; i<n; ++i) {
		xp[i]=X[ip[i]];
		fp[i]=f[ip[i]];
	}
 
/*	Calculate the cubic spline interpolation using only N chosen knots
	and compare with full result to calculate the maximum error */
	i=spline(xp,fp,n,c);
	errmax=0.0;
	for(i=0; i<n0; ++i) {
		xb=X0[i];
		fb=splevl(xb,n,xp,fp,c,&dfb,&ddfb,&ier);
		df=fabs(fb-F0[i]);
		ERR[i]=fb-F0[i];
		if(df>errmax) errmax=df;
	}
	return errmax;
}



/*	To change the set of knots selected for interpolation */

void change(double temp, int *idum, int ip[])

{
	int i,j,n,ip1;

	n=12;
	if(temp>0.05) {
		for(i=0; i<n; ++i) {
 
/*	Shift the index of all knots by at most 1 randomly and
	maintain the ordering of knots */
			ip1=ip[i]+3*ranf(idum)-1;
			if(i>0) {
				if(ip1<= ip[i-1]+1) ip1=ip[i-1]+2;
			}
			else {
				if(ip1<0) ip1=0;
			}
			if(ip1>50+2*i-2*n) ip1=50+2*i-2*n;
			ip[i]=ip1;
		}
	} else {

/*	If the temperature is less than 0.05 then shift only one knot at a time */
		i=n*ranf(idum);
		if(ranf(idum)<0.5) ip1=ip[i]-1;
		else ip1=ip[i]+1;

		if(i>0) {
			if(ip1<=ip[i-1]) ip1=ip[i-1]+1;
		} else {
			if(ip1<0) ip1=0;
		}

		if(i<n-1) {
			if(ip1>= ip[i+1]) ip1=ip[i+1]-1;
		} else {
			if(ip1>48) ip1=48;
		}
		ip[i]=ip1;
	}
	return;
}
