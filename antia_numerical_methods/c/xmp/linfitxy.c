/*	Straight line fit to data when there are errors in both x and y */

#include <stdio.h>
#include <math.h>

double ran1(double *seed);
int linfitxy(int n, double *x, double *y, double sigx, double sigy,
	double rho, double *xi, double *yi, double *a, double *b,
	double *chi);

main()
{
	int i,n,ier;
	double seed,r,sx,sy,rxy,c,rm,chi;
	double x[100],y[100],xi[100],yi[100];

	seed=9;
	n=100;
/*	Generate the data */
	for(i=0; i<100; ++i) {
		r=ran1(&seed);
		x[i]=i*0.01+0.001*(r-0.5);
/*	There is no correlation as independent random numbers are used */
		r=ran1(&seed);
		y[i]=i*0.005+0.01*(r-0.5);
	}
/*	Set the error and correlation	*/
	sx=3.e-4; sy=3.e-3; rxy=0;
	ier=linfitxy(n,x,y,sx,sy,rxy,xi,yi,&c,&rm,&chi);
	printf(" ier = %d  slope = %f  intercept = %f \n ",ier,rm,c);
	printf(" \n      x    fitted-x       y    fitted-y \n");
	for(i=0; i<100; ++i) 
		printf(" %f  %f  %f  %f \n",x[i],xi[i],y[i],yi[i]);
}

/*	Least squares straight line fit when there is error in both x and y
        To fit equation of form y=a+b*x

        N : (input) Number of data points
        X : (input) Real array of length N containing the x values
        Y : (input) Real array of length N containing the y values
        SIGX : (input) estimated error in x values, assumed to be the
                same for all points
        SIGY : (input) estimated error in y values, assumed to be the
                same for all points
        RHO : (input) estimated correlation between errors in x and y
                assumed to be the same for all points
        XI : (output) Reall array of length N which will give the fitted x values
        YI : (output) Reall array of length N which will give the fitted y values
        A : (output) fitted value of the intercept
        B : (output) fitted value of the slope
        CHI : (output) value of chi^2 at the minimum

	Error status is returned by the value of the function LINFITXY.
        	0 value implies successful execution
                617 implies that discriminant of quadratic equation
                       is negative and calculations are aborted

	Required functions : none
*/

#include <math.h>

int linfitxy(int n, double *x, double *y, double sigx, double sigy,
	double rho, double *xi, double *yi, double *a, double *b,
	double *chi)
{
	int i;
	double r,s,sx,sy,sxy,sxx,syy,c0,c1,c2,del,b1,a1,t,t1,chi1;
	double x1,x2,y1,y2;

	r=rho*sigy/sigx; s=sigy*sigy/(sigx*sigx);
	sx=0.0; sy=0.0; sxy=0.0; sxx=0.0; syy=0.0;
	for(i=0; i<n; ++i){
		sx=sx+x[i]; sy=sy+y[i];
		sxy=sxy+x[i]*y[i]; sxx=sxx+x[i]*x[i]; syy=syy+y[i]*y[i];
	}
	sx=sx/n; sy=sy/n; sxy=sxy/n; sxx=sxx/n; syy=syy/n;

/*	Find the quadratic in the slope and solve it */
	c2=sigx*sigx*(sxy-sx*sy)+rho*sigx*sigy*(sx*sx-sxx);
	c1=sigx*sigx*(sy*sy-syy)+sigy*sigy*(sxx-sx*sx);
	c0=-rho*sigx*sigy*(sy*sy-syy)-sigy*sigy*(sxy-sx*sy);
	del=c1*c1-4*c0*c2;
	if(del<0.0) return 617;
	del=sqrt(del);
	if(c1>0.0) {b1=(-c1-del)/(2*c2);}
	else b1=(-c1+del)/(2*c2);
	*b=c0/(c2*b1);
	a1=sy-b1*sx;
	*a=sy-(*b)*sx;
	t1=(b1*r-s)/(b1-r);
	t=((*b)*r-s)/((*b)-r);

/*	Choose the solution with smaller chi^2 */
	chi1=0.0; *chi=0.0;
	for(i=0; i<n; ++i){
		x1=(t1*x[i]-y[i]+a1)/(t1-b1);
		x2=(t*x[i]-y[i]+(*a))/(t-(*b));
		y1=a1+b1*x1; y2=(*a)+(*b)*x2;
		chi1=chi1+(x[i]-x1)*(x[i]-x1)/(sigx*sigx)-2*rho*(x[i]-x1)*(y[i]-y1)
			/(sigx*sigy)+(y[i]-y1)*(y[i]-y1)/(sigy*sigy);
		*chi=*chi+(x[i]-x2)*(x[i]-x2)/(sigx*sigx)-2*rho*(x[i]-x2)*(y[i]-y2)
			/(sigx*sigy)+(y[i]-y2)*(y[i]-y2)/(sigy*sigy);
	}

	if(chi1< (*chi)){
		*chi=chi1; *a=a1; *b=b1; t=t1;
	}

	for(i=0; i<n; ++i){
		xi[i]=(t*x[i]-y[i]+(*a))/(t-(*b));
		yi[i]=(*a)+(*b)*xi[i];
	}
	return 0;
}


/*	To generate uniformly distributed random numbers in interval (0,1)

	SEED : (input/output) is a real value used as the seed
		It should be positive during initial call and
		should not be modified between different calls

	Required functions : None
*/

#include <math.h>

double ran1(double *seed)

{
	double am=2147483648e0, a=45875e0, ac=453816693e0, an=2147483647e0;

	*seed=fmod((*seed)*a+ac,am);
	return (*seed)/an;
}
	


