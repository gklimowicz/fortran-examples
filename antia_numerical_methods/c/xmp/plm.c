/*	To calculate associated Legendre function */

#include <stdio.h>
#include <math.h>

void pleg(int l, double x, double p[]);
void plm(int l, int m, double x, double p[]);

main()
{
	int i,i1,l,m,n;
	double x,theta,phi,f1,f2,p[5000],p1[5000];

	for(i1=0; i1<99; ++i1) {
		printf("type l, m, x,        (quits when x=-9999)\n");
		scanf(" %d %d %le",&l,&m,&x);
		if(x==-9999) return 0;

		pleg(l,x,p);
		plm(l,m,x,p1);
		printf(" l = %d    m = %d   x = %e\n  Pl(x) = %e    Plm(x) = %e\n",l,m,x,p[l],p1[l]);
	}
	return;
}




/*	To calculate Legendre polynomial P_l(X)

	L : (input) Order of polynomial (L>=0)
	X : (input) Argument at which the value of polynomial is required
	P : (output) Array of length L+1, which will contain the
		calculated values of polynomials. P[j] will contain P_j(X)

	Required functions : None
*/

void pleg(int l, double x, double p[])

{
	int n;
	
	if(l<0) return;
	p[0]=1.0;
	p[1]=x;
	for(n=2; n<=l; ++n) p[n]=((2*n-1)*x*p[n-1]-(n-1)*p[n-2])/n;
	return;
}



/*	To calculate the associated Legendre functions P_lm(X)

	L,M : (input) Order of function (L>=0), abs(M)<=L
	X : (input) Argument at which the value of polynomial is required
	P : (output) Array of length L+1, which will contain the
		calculated values of polynomials. P[j] will contain
		P_jM(X) for j>=M

	Required functions : None
*/

#include <math.h>

void plm(int l, int m, double x, double p[])

{
	int i,n,mm;
	double pm,rm;

	if(l<0 || abs(m)>l ) return;
 
/*	First compute P_MM(x) */
	mm=abs(m);
	pm=1.0;
	for(i=1; i<2*mm; i += 2) pm=pm*i;

	if(m<0) {
/*	Modify the normalisation factor */
		for(i=1; i<=2*mm; ++i) pm=pm/i;
		if(mm - 2*(mm/2) == 1) pm=-pm;
	}

	rm=mm/2.0;
	if(mm == 0) p[mm]=1.0;
	else if(fabs(x)<1.0) p[mm]=pm*pow(1.0-x*x, rm);
	else p[mm]=0.0;

/*	Use the recurrence relation to compute P_nM(x) */
	p[mm+1]=x*p[mm]*(2*mm+1)/(mm+1-m);
	for(n=mm+2; n<=l; ++n) p[n]=((2*n-1)*x*p[n-1]-(n-1+m)*p[n-2])/(n-m);

	return;
}
