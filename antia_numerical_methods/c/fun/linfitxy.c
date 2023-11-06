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


