/*	To calculate Bessel functions of various types of integral order */

#include <stdio.h>
#include <math.h>

double bj0(double xg);
double bj1(double xg);
void bjn(int n, double xb, double bj[]);
double by0(double xg);
double by1(double xg);
void byn(int n, double xb, double bj[]);
void bjy0(double xb, double *bj0, double *by0);
void bjy1(double xb, double *bj1, double *by1);
void sphbjn(int n, double xb, double bj[]);

double bi0(double xg);
double bi1(double xg);
void bin(int n, double xb, double bj[]);
double bk0(double xg);
double bk1(double xg);
void bkn(int n, double xb, double bj[]);

main()
{
	int i,i1,j,k,n;
	double x,f1,f2,g1,g2,bj[2000];

	for(i1=0; i1<99; ++i1) {
		printf("type x,n        (quits when x=-9999)\n");
		scanf(" %le %d", &x, &n);
		if(x==-9999) return 0;

/*	Bessel function of first kind */
		f1=bj0(x);
		f2=bj1(x);
		bjn(n,x,bj);
		k=abs(n);
		printf(" x = %e    n = %d   J0 = %e    J1 = %e \n",x,n,f1,f2);
		printf("bjn:   J0 = %e     J1 = %e    Jn = %e \n",bj[0],bj[1],bj[k]);

/*	Bessel function of first and second kind */
		bjy0(x, &f1, &g1);
		bjy1(x, &f2, &g2);
		printf(" J0= %e   Y0= %e   J1 = %e   Y1 = %e\n",f1,g1,f2,g2);

/*	Bessel function of second kind */
		f1=by0(x);
		f2=by1(x);
		byn(n,x,bj);
		printf(" Y0 = %e    Y1 = %e \n",f1,f2);
		printf("byn:   Y0 = %e    Y1 = %e    Yn = %e\n",bj[0],bj[1],bj[k]);

/*	Spherical Bessel function */
		sphbjn(n,x,bj);
		printf(" j0= %e    j1 = %e     jn = %e\n",bj[0],bj[1],bj[k]);

/*	Modified Bessel function of first kind */
		f1=bi0(x);
		f2=bi1(x);
		bin(n,x,bj);
		printf("  I0= %e     I1 = %e\n",f1,f2);
		printf("bin:    I0= %e     I1 = %e     In = %e\n",bj[0],bj[1],bj[k]);

/*	Modified Bessel function of second kind */
		f1=bk0(x);
		f2=bk1(x);
		bkn(n,x,bj);
		printf("  K0= %e     K1 =%e\n",f1,f2);
		printf("bkn:    K0= %e     K1 =%e     Kn = %e\n",bj[0],bj[1],bj[k]);

	}
	return;
}



/*	To calculate the Bessel function of order zero for real argument 

	Required functions : None
*/
 
#include <math.h>

#define PI 3.14159265358979324

double bj0(double x)

{
	double y,fn,fd,f,p0,q0;
 
/*	Coefficients of rational function approximations */
	double a[8] = {1.293686560051304152e-02,  8.573459862295151747e-05,
                       3.854769244046149702e-07,  1.308534328117880493e-09,
                       3.512360907188715842e-12,  7.512575042421009221e-15,
                       1.229302278444845702e-17,  1.311883486088925264e-20};
	double b[8] = {9.999999999999999878e-01, -2.370631343994868513e-01,
                       1.247651819849453565e-02, -2.529374255010058573e-04,
                       2.411267406461247155e-06, -1.159484705672466498e-08,
                       2.730546745501229851e-11, -2.517936655103065990e-14};

	double a1[5] = { 8.911849018950665793e+01,
                         2.078818787053760203e+03,  1.366258799766718466e+04,
                         1.800383785973922830e+04,  4.923440494847201509e+02};
	double b1[5] = { 9.999999999999999908e-01,
                         8.904817768950681616e+01,  2.072664795311476688e+03,
                         1.352584337655551999e+04,  1.723138433448795889e+04};

	double a2[5] = { 1.046366195300779895e+02,  2.980727727381642723e+03,
                         2.570418404044668245e+04,  5.291014161889741749e+04,
                         9.497228391199055149e+03};
	double b2[5] = {-1.249999999999999534e-01, -1.300633525376058087e+01,
                        -3.651542590150084093e+02, -3.016744074771875522e+03,
                        -5.251679479249748063e+03};
 
	if(fabs(x)<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=(((((((a[7]*y+a[6])*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=fn/fd;
	}

	else {
/*	Use rational function approximations for P_0 and Q_0 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p0=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q0=fn/(fd*fabs(x));
 
		y=fabs(x);
		f=sqrt(2.0/(PI*y))*(p0*cos(y-PI/4.0)-q0*sin(y-PI/4.0));
	}

	return f;
}


/*	To calculate the Bessel function of order one for real argument 

	Required functions : None
*/
 
#include <math.h>

#define PI 3.14159265358979324

double bj1(double x)

{
	double y,fn,fd,f,p1,q1,x2;
 
/*	Coefficients of rational function approximations */
	double a[7] = {1.156878510067059849e-02,  6.749406787503859073e-05,
                       2.614560837317535451e-07,  7.408815126464007290e-10,
                       1.577637796406197189e-12,  2.432945305413635549e-15,
                       2.257446839754248784e-18};
	double b[8] = {5.000000000000000074e-01, -5.671560744966475746e-02,
                       1.914864631812345532e-03, -2.821407888958585592e-05,
                       2.103168789640803591e-07, -8.322474383730280556e-10,
                       1.678871778708754849e-12, -1.372424374400306547e-15};
	double a1[5] = { 8.659888261699365129e+01,
                         1.932665751369749084e+03,  1.172714583536277145e+04,
                         1.256737699073784218e+04, -6.147124347503755010e+02};
	double b1[5] = { 1.000000000000000011e+00,
                         8.671607011699346720e+01,  1.942669862370300601e+03,
                         1.194181952104744095e+04,  1.371467864585746530e+04};
	double a2[5] = { 1.021472573795463627e+02,  2.807865400111916226e+03,
                         2.280402060738415865e+04,  4.121116954504273053e+04,
                         3.501974669280301705e+03};
	double b2[5] = { 3.749999999999999461e-01,  3.820268245483084309e+01,
                         1.042753017477090289e+03,  8.289951986135169400e+03,
                         1.371889615877945967e+04};
 
	if(fabs(x)<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=x*fn/fd;
	}

	else {
/*	Use rational function approximations for P_1 and Q_1 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p1=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q1=fn/(fd*fabs(x));
 
		y=fabs(x);
		x2=y-0.75*PI;
		f=sqrt(2.0/(PI*y))*(p1*cos(x2)-q1*sin(x2));
		if(x<0.0) f=-f;
	}

	return f;
}


/*	To calculate the Bessel function of integral order for real argument

	N : (input) Order of Bessel function required, N may be negative
		or positive. For N=0,1 use BJ0 and BJ1 respectively.
	XB : (input) Argument at which the value is required
	BJ : (output) Array of length at least 
		abs(N)+16+MAX(25,5*sqrt(N))
		which will contain the value of Bessel function of order
		0,1,...,abs(N). BJ[i] will contain Bessel function of order i
		or -i (if N<0)
		Remaining elements of array are used as scratch space

	Required functions : BJ0, BJ1
*/

#include <math.h>

#include <stdio.h>

double bj0(double x);
double bj1(double x);

void bjn(int n, double xb, double bj[])

{
	int i,j,n1,na;
/*	REPS should be less than the machine accuracy */
	double x,s,xa,t0,t, reps=1.e-17;

	x=fabs(xb);
	na=abs(n);
	if(xb==0.0) {
		bj[0]=1.0;
		for(i=1; i<=na; ++i) bj[i]=0.0;
		return;
	}

	else if(na<x || na<2) {
/*	Use the recurrence relation in the forward direction  */
		bj[0]=bj0(x);
		bj[1]=bj1(x);
		for(i=2; i<=na; ++i) bj[i]=2*(i-1)*bj[i-1]/x-bj[i-2];
	}

	else if(x<=4.0) {
/*	Use series expansion to calculate  bj[na], bj[na-1] */
		xa=x*x/4.0;
		t0=x/2.0;
		for(i=2; i<=na; ++i) {
			t0=t0*x/(2.0*i);
			if(i>=na-1) {
				t=t0; s=t0;
				for(j=1; j<51; ++j) {
					t=-t*xa/(j*(j+i));
					s=s+t;
					if(fabs(t)<fabs(s*reps)) break;
				}
				bj[i]=s;
			}
		}
		for(i=na-1; i>=1; --i) bj[i-1]=2.*i*bj[i]/x-bj[i+1];
	}

	else {
/*	Use the recurrence relation in the backward direction  */
		n1=5.*sqrt((double) na); if(n1<25) n1=25;
		if(x<na/6.0) n1=n1/log(na*0.5/x);
		n1=n1+na;
		if(na<x+15.0) n1=n1+15;
		if(n1-na<15) n1=na+15;
		bj[n1]=0.0; bj[n1-1]=1.0;
		for(i=n1-1; i>=1; --i) bj[i-1]=2.*i*bj[i]/x-bj[i+1];

		s=bj[0];
		for(i=2; i<n1; i += 2) s=s+2.*bj[i];

		for(i=0; i<=na; ++i) bj[i]=bj[i]/s;
/*	If fabs(bj(na+1))<1/reps, then the required accuracy may not be achieved
	hence printout an error message */
		if(fabs(bj[na+1]*reps)< 1.0) printf(" bjn failed at n = %d  x = %e s = %e \n", n,x,bj[na+1]);
	}

	if((n<0) != (xb<0.0)) {
		for(i=1; i<=na; i += 2) bj[i]= -bj[i];
	}
	return;
}


/*	To calculate the Bessel function of first and second kind of order zero
	for real positive argument
	For XB<=0 the function of second kind is not defined and function
	will return zero value  without any error message or flag

	XB : (input) Argument at which the functions need to be evaluated
	BJ0 : (output) Calculated value of Bessel function of first kind
	BY0 : (output) Calculated value of Bessel function of second kind

	Required functions : None
*/

#include <math.h>

#define PI 3.14159265358979324
#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

void bjy0(double xb, double *bj0, double *by0)

{
	double x,y,fn,fd,f,p0,q0;
 
/*	Coefficients of rational function approximations */
	double a[8] = {1.293686560051304152e-02,  8.573459862295151747e-05,
                       3.854769244046149702e-07,  1.308534328117880493e-09,
                       3.512360907188715842e-12,  7.512575042421009221e-15,
                       1.229302278444845702e-17,  1.311883486088925264e-20};
	double b[8] = {9.999999999999999878e-01, -2.370631343994868513e-01,
                       1.247651819849453565e-02, -2.529374255010058573e-04,
                       2.411267406461247155e-06, -1.159484705672466498e-08,
                       2.730546745501229851e-11, -2.517936655103065990e-14};

	double a0[7] = {1.089079731266387424e-02,  5.954632605213292419e-05,
                        2.150109922530480401e-07,  5.641082188778387960e-10,
                        1.102341761343675716e-12,  1.539990321465010920e-15,
                        1.263081729204829936e-18};
	double b0[8] = {2.500000000000000006e-01, -2.071480067183403591e-02,
                        5.553511120900719150e-04, -6.804373640943337406e-06,
                        4.346149688717712144e-08, -1.505635199331021665e-10,
                        2.703193275976574669e-13, -1.993047807317608951e-16};

	double a1[5] = { 8.911849018950665793e+01,
                         2.078818787053760203e+03,  1.366258799766718466e+04,
                         1.800383785973922830e+04,  4.923440494847201509e+02};
	double b1[5] = { 9.999999999999999908e-01,
                         8.904817768950681616e+01,  2.072664795311476688e+03,
                         1.352584337655551999e+04,  1.723138433448795889e+04};

	double a2[5] = { 1.046366195300779895e+02,  2.980727727381642723e+03,
                         2.570418404044668245e+04,  5.291014161889741749e+04,
                         9.497228391199055149e+03};
	double b2[5] = {-1.249999999999999534e-01, -1.300633525376058087e+01,
                        -3.651542590150084093e+02, -3.016744074771875522e+03,
                        -5.251679479249748063e+03};
 
	x=fabs(xb);

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=(((((((a[7]*y+a[6])*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		*bj0=fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		if(x>0.0) *by0=2.*(*bj0*(log(x/2)+EUGAM)+y*fn/fd)/PI;
	}

	else {
/*	Use rational function approximations for P_0 and Q_0 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p0=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q0=fn/(fd*fabs(x));
 
		y=fabs(x);
		*by0=sqrt(2.0/(PI*y))*(p0*sin(y-PI/4.0)+q0*cos(y-PI/4.0));
		*bj0=sqrt(2.0/(PI*y))*(p0*cos(y-PI/4.0)-q0*sin(y-PI/4.0)); 
	}

	if(xb<=0.0) *by0=0.0;
	return;
}



/*	To calculate the Bessel function of first and second kind of order one
	for real positive argument
	For XB<=0 the function of second kind is not defined and function
	will return zero value  without any error message or flag.

	XB : (input) Argument at which the functions need to be evaluated
	BJ1 : (output) Calculated value of Bessel function of first kind
	BY1 : (output) Calculated value of Bessel function of second kind

	Required functions : None
*/
 
#include <math.h>

#define PI 3.14159265358979324
#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

void bjy1(double xb, double *bj1, double *by1)

{
	double x,y,fn,fd,f,p1,q1,x2;
 
/*	Coefficients of rational function approximations */
	double a[7] = {1.156878510067059849e-02,  6.749406787503859073e-05,
                       2.614560837317535451e-07,  7.408815126464007290e-10,
                       1.577637796406197189e-12,  2.432945305413635549e-15,
                       2.257446839754248784e-18};
	double b[8] = {5.000000000000000074e-01, -5.671560744966475746e-02,
                       1.914864631812345532e-03, -2.821407888958585592e-05,
                       2.103168789640803591e-07, -8.322474383730280556e-10,
                       1.678871778708754849e-12, -1.372424374400306547e-15};

	double a0[7] = {1.186963690270342970e-02,  7.123839029323002736e-05,
                        2.848196759538669928e-07,  8.365303089083305885e-10,
                        1.857096246589275584e-12,  3.012506935004947491e-15,
                        2.996951174746838817e-18};
	double b0[8] = {2.500000000000000299e-01, -7.515759077432437273e-02,
                        3.430771992327672576e-03, -6.022315614557372919e-05,
                        5.067136874996839630e-07, -2.197514674456554803e-09,
                        4.768619679411702014e-12, -4.139491442515065355e-15};
    
	double a1[5] = { 8.659888261699365129e+01,
                         1.932665751369749084e+03,  1.172714583536277145e+04,
                         1.256737699073784218e+04, -6.147124347503755010e+02};
	double b1[5] = { 1.000000000000000011e+00,
                         8.671607011699346720e+01,  1.942669862370300601e+03,
                         1.194181952104744095e+04,  1.371467864585746530e+04};

	double a2[5] = { 1.021472573795463627e+02,  2.807865400111916226e+03,
                         2.280402060738415865e+04,  4.121116954504273053e+04,
                         3.501974669280301705e+03};
	double b2[5] = { 3.749999999999999461e-01,  3.820268245483084309e+01,
                         1.042753017477090289e+03,  8.289951986135169400e+03,
                         1.371889615877945967e+04};
 
	x=fabs(xb);

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		*bj1=x*fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		if(x>0.0) *by1=2.*(*bj1*(log(x/2)+EUGAM)-1.0/x-x*fn/fd)/PI;
	}

	else {
/*	Use rational function approximations for P_1 and Q_1 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p1=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q1=fn/(fd*fabs(x));
 
		y=fabs(x);
		x2=y-0.75*PI;
		*by1=sqrt(2.0/(PI*y))*(p1*sin(x2)+q1*cos(x2));
		*bj1=sqrt(2.0/(PI*y))*(p1*cos(x2)-q1*sin(x2));
	}

	if(xb<=0.0) *by1=0.0;
	return;
}


/*	To calculate the modified Bessel function of first kind of order
	zero for real argument

	Required functions : None
*/

#include <math.h>

double bi0(double x)

{
	double xa,y,fn,fd,f;

/*	Coefficients of rational function approximations */
	double a[8] = {-1.212778758454003780e-02,  7.395384394332586467e-05,
                       -2.981475020389602859e-07,  8.749000589372218583e-10,
                       -1.925066925538111917e-12,  3.116611043626875576e-15,
                       -3.403826199284281585e-18,  1.919794284150553073e-21};
	double b[8] = { 9.999999999999999887e-01,  2.378722124154600447e-01,
                        1.266700694780821586e-02,  2.627214102531807755e-04,
                        2.599754169250266946e-06,  1.322628469861000483e-08,
                        3.380151377715495065e-11,  3.502398414704805956e-14};
 
	double a1[9] = {-4.368454479954000936e+01,  8.662328585154270157e+02,
                        -9.865641284565305256e+03,  6.940667740497879456e+04,
                        -3.022161789164256804e+05,  7.741268251525742446e+05,
                        -1.016413121409283393e+06,  5.111938862294702143e+05,
                        -5.159799972015923803e+04};
	double b1[10] = { 3.989422804014326451e-01, -1.737774413557676433e+01,
                          3.434265110737608756e+02, -3.893820477243345847e+03,
                          2.722034716144037839e+04, -1.173589448296665233e+05,
                          2.954547271385017668e+05, -3.737076509558829360e+05,
                          1.685151343593986274e+05, -9.954986502601715062e+03};
 
	xa=fabs(x);
	if(xa<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=(((((((a[7]*y+a[6])*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=fn/fd;
	}

	else {
/*	Use rational function approximation to the asymptotic form */
		y=1./xa;
        fn=((((((((b1[9]*y+b1[8])*y+b1[7])*y+b1[6])*y+b1[5])*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
        fd=((((((((a1[8]*y+a1[7])*y+a1[6])*y+a1[5])*y+a1[4])*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=exp(xa)*fn/(fd*sqrt(xa));
	}

	return f;
}


/*	To calculate the modified Bessel function of first kind of order
	one for real argument

	Required functions : None
*/
 
#include <math.h>

double bi1(double x)

{
	double xa,y,fn,fd,f;

/*	Coefficients of rational function approximations */
	double a[7] = {-1.070229423222827071e-02,  5.648033546364466297e-05,
                       -1.918252346392549379e-07,  4.558004979539079070e-10,
                       -7.621816709677626962e-13,  8.342527256774311367e-16,
                       -4.621507477175875649e-19};
	double b[8] = { 4.999999999999999914e-01,  5.714885288388592026e-02,
                        1.963513444884162262e-03,  2.981702267497806226e-05,
                        2.328548099504502597e-07,  9.862614413072000150e-10,
                        2.192079167003372805e-12,  2.056894175269540721e-15};
 
	double a1[9] = {-4.354982033508663805e+01,  8.607755018195304592e+02,
                        -9.764617274576599078e+03,  6.833415275089609019e+04,
                        -2.952159059017153952e+05,  7.462070230806506945e+05,
                        -9.526843724015979184e+05,  4.423280833889137276e+05,
                        -2.993689102971858935e+04};
	double b1[10] = { 3.989422804014327116e-01, -1.752346799070278231e+01,
                          3.498681897994668814e+02, -4.022298493285640672e+03,
                          2.868368405169782754e+04, -1.275734365963396269e+05,
                          3.390186096946305841e+05, -4.802175709866444585e+05,
                          2.931298028513980846e+05, -5.264175913253416995e+04};
 
	xa=fabs(x);
	if(xa<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=x*fn/fd;
	}

	else {
/*	Use rational function approximation to the asymptotic form */
		y=1./xa;
        fn=((((((((b1[9]*y+b1[8])*y+b1[7])*y+b1[6])*y+b1[5])*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
        fd=((((((((a1[8]*y+a1[7])*y+a1[6])*y+a1[5])*y+a1[4])*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=exp(xa)*fn/(fd*sqrt(xa));
		if(x<0.0) f= -f;
	}

	return f;
}

 
/*	To calculate the modified Bessel function of first kind of
		integral order for real argument

	N : (input) Order of Bessel function required, N must be positive
	XB : (input) Argument at which the value is required
	BI : (output) Array of length at least N+16+MAX(25, 5*sqrt(N))
		which will contain the value of Bessel function of order
		0,1,...,n. BI[i] will contain Bessel function of order i
		Remaining elements of array are used as scratch space

	Required functions : BI0, BI1
*/
#include <stdio.h>
#include <math.h>

double bi0(double x);
double bi1(double x);

void bin(int n, double xb, double bi[])

{
	int i,j,n1,na;
/*	REPS should be less than the machine accuracy */
	double x,s,t,t0,xa, reps=1.e-17;

	x=fabs(xb);
	na=abs(n);
	if(xb==0.0) {
		bi[0]=1.0;
		for(i=1; i<=na; ++i) bi[i]=0.0;
		return;
	}

	else if(na<x-10 || na<2) {
/*	Use the recurrence relation in the forward direction  */
		bi[0]=bi0(x);
		bi[1]=bi1(x);
		for(i=2; i<=na; ++i) bi[i]=-2*(i-1)*bi[i-1]/x+bi[i-2];
	}

	else if(x<=4.0) {
/*	Use series expansion to calculate  bi[na], bi[na-1] */
		xa=x*x/4.0;
		t0=x/2.0;
		for(i=2; i<=na; ++i) {
			t0=t0*x/(2.0*i);
			if(i>=na-1) {
				t=t0; s=t0;
				for(j=1; j<51; ++j) {
					t=t*xa/(j*(j+i));
					s=s+t;
					if(fabs(t)<fabs(s*reps)) break;
				}
				bi[i]=s;
			}
		}
		for(i=na-1; i>=1; --i) bi[i-1]=2.*i*bi[i]/x+bi[i+1];
	}


	else {
/*	Use the recurrence relation in the backward direction  */
		n1=5.*sqrt((double) na); if(n1<25) n1=25;
		if(x<na/6.0) n1=n1/log(na*0.5/x);
		n1=n1+na;
		if(na<x+15.0) n1=n1+15;
		if(n1-na<15) n1=na+15;
		bi[n1]=0.0; bi[n1-1]=1.0;
		for(i=n1-1; i>=1; --i) bi[i-1]=2.*i*bi[i]/x+bi[i+1];

		s=bi[0]/bi0(x);
		for(i=0; i<=na; ++i) bi[i]=bi[i]/s;
/*	If fabs(bi[na+1])<1/reps, then the required accuracy may not be achieved
	hence printout an error message */
		if(fabs(bi[na+1]*reps)< 1.0) printf(" bin failed at n = %d  x = %e  s = %e \n", n,x,bi[na+1]);
	}

	if(xb<0.0) {
		for(i=1; i<=na; i += 2) bi[i]= -bi[i];
	}
	return;
}

 


/*	To calculate the modified Bessel function of second kind of order
	zero for real argument

	Required functions : None
*/
 
#include <math.h>

#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

double bk0(double x)

{
	double y,fn,fd,f,bi0;
 
/*	Coefficients of rational function approximations */
	double a[7] = {-1.011123277211023982e-02,  5.027323349741487160e-05,
                       -1.604132766768703653e-07,  3.571438627971404230e-10,
                       -5.582724661115911120e-13,  5.702487045740630357e-16,
                       -2.945917638250090849e-19};
	double b[8] = { 2.499999999999999993e-01,  2.090969180697244468e-02,
                        5.713038828706075545e-04,  7.220998182565402322e-06,
                        4.832471102493292948e-08,  1.789925692526897035e-10,
                        3.530871144986696274e-13,  2.972558105712627660e-16};
 
	double a0[8] = {-1.212778758454003780e-02,  7.395384394332586467e-05,
                        -2.981475020389602859e-07,  8.749000589372218583e-10,
                        -1.925066925538111917e-12,  3.116611043626875576e-15,
                        -3.403826199284281585e-18,  1.919794284150553073e-21};
	double b0[8] = { 9.999999999999999887e-01,  2.378722124154600447e-01,
                         1.266700694780821586e-02,  2.627214102531807755e-04,
                         2.599754169250266946e-06,  1.322628469861000483e-08,
                         3.380151377715495065e-11,  3.502398414704805956e-14};
 
	double a1[5] = { 1.134095488162070337e+01,  4.133208417436337182e+01,
                         5.704583485964969346e+01,  2.690522251864472800e+01,
                         2.767771031868579433e+00};
	double b1[6] = { 1.253314137315500251e+00,  1.405711481662802554e+01,
                         5.011348668524983318e+01,  6.592870149979143648e+01,
                         2.752549950796738039e+01,  1.796256223905248475e+00};
 
	if(x <= 0.0) return 0.0;

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=(((((((a0[7]*y+a0[6])*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		bi0=fn/fd;

		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=y*fn/fd-bi0*(EUGAM+log(x/2));
	}

	else {
/*	Use rational function approximations for P_0 and Q_0 */
		y=1./x;
		fn=((((b1[5]*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=exp(-x)*fn/(fd*sqrt(x));
	}

	return f;
}


/*	To calculate the modified Bessel function of second kind of order
	one for real argument

	Required functions : None
*/
 
#include <math.h>

#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

double bk1(double x)

{
	double y,fn,fd,f,bi1;
 
/*	Coefficients of rational function approximations */
	double a[7] = {-1.070229423222827071e-02,  5.648033546364466297e-05,
                       -1.918252346392549379e-07,  4.558004979539079070e-10,
                       -7.621816709677626962e-13,  8.342527256774311367e-16,
                       -4.621507477175875649e-19};
	double b[8] = { 4.999999999999999914e-01,  5.714885288388592026e-02,
                        1.963513444884162262e-03,  2.981702267497806226e-05,
                        2.328548099504502597e-07,  9.862614413072000150e-10,
                        2.192079167003372805e-12,  2.056894175269540721e-15};
 
	double a0[7] = { -1.097271232519001047e-02,  5.948570919145243152e-05,
                         -2.079707688602524007e-07,  5.097946369393536825e-10,
                         -8.813373771643053620e-13,  9.993295069392185250e-16,
                         -5.743938665570093767e-19};
	double b0[8] = { 2.499999999999999663e-01,  7.538182191870271423e-02,
                         3.497906054669938529e-03,  6.321709197334740028e-05,
                         5.569209160486120478e-07,  2.585145364373725340e-09,
                         6.185681870407510042e-12,  6.175435537988524458e-15};
 
	double a1[5] = { 9.546087640477785968e+00,  2.842458429353285707e+01,
                         3.039341689594721892e+01,  9.857685145186666128e+00,
                         4.376813258161116587e-01};
	double b1[6] = { 1.253314137315500259e+00,  1.243423939735687111e+01,
                         3.996465306538935273e+01,  5.017830258820480688e+01,
                         2.351074680068346432e+01,  2.993271886379231231e+00};
 
	if(x <= 0.0) return 0.0;

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		bi1=x*fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		f= -x*fn/fd+bi1*(EUGAM+log(x/2))+1./x;
	}

	else {
/*	Use rational function approximations for P_0 and Q_0 */
		y=1./x;
		fn=((((b1[5]*y+b1[4])*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		f=exp(-x)*fn/(fd*sqrt(x));
	}

	return f;
}
 


/*	To calculate the modified Bessel function of second kind of
		integral order for real argument
		The function value is not defined for x<=0 and no calculations
		are done in that case but no warning is issued

	N : (input) Order of Bessel function required, N must be positive
	X : (input) Argument at which the value is required
	BK : (output) Array of length at least abs(N)+1
		which will contain the value of Bessel function of order
		0,1,...,N. BK[i] will contain Bessel function of order i

	Required functions : BK0, BK1
*/

#include <math.h>

double bk0(double x);
double bk1(double x);

void bkn(int n, double x, double bk[])

{
	int i;
 
	bk[0]=bk0(x);
	bk[1]=bk1(x);
	if(x<= 0.0) return;
/*	Use the recurrence relation in Forward direction */
	for(i=2; i<=abs(n); ++i) bk[i]=2*(i-1)*bk[i-1]/x+bk[i-2];
	return;
}
 


/*	To calculate the Bessel function of second kind of order zero
	for real positive argument
	For X<=0 the function is not defined and function will return
	a zero value without any error message or flag.

	Required functions : None
*/
 
#include <math.h>

#define PI 3.14159265358979324
#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

double by0(double x)

{
	double y,fn,fd,f,p0,q0,bj0;
 
/*	Coefficients of rational function approximations */
	double a[8] = {1.293686560051304152e-02,  8.573459862295151747e-05,
                       3.854769244046149702e-07,  1.308534328117880493e-09,
                       3.512360907188715842e-12,  7.512575042421009221e-15,
                       1.229302278444845702e-17,  1.311883486088925264e-20};
	double b[8] = {9.999999999999999878e-01, -2.370631343994868513e-01,
                       1.247651819849453565e-02, -2.529374255010058573e-04,
                       2.411267406461247155e-06, -1.159484705672466498e-08,
                       2.730546745501229851e-11, -2.517936655103065990e-14};

	double a0[7] = {1.089079731266387424e-02,  5.954632605213292419e-05,
                        2.150109922530480401e-07,  5.641082188778387960e-10,
                        1.102341761343675716e-12,  1.539990321465010920e-15,
                        1.263081729204829936e-18};
	double b0[8] = {2.500000000000000006e-01, -2.071480067183403591e-02,
                        5.553511120900719150e-04, -6.804373640943337406e-06,
                        4.346149688717712144e-08, -1.505635199331021665e-10,
                        2.703193275976574669e-13, -1.993047807317608951e-16};

	double a1[5] = { 8.911849018950665793e+01,
                         2.078818787053760203e+03,  1.366258799766718466e+04,
                         1.800383785973922830e+04,  4.923440494847201509e+02};
	double b1[5] = { 9.999999999999999908e-01,
                         8.904817768950681616e+01,  2.072664795311476688e+03,
                         1.352584337655551999e+04,  1.723138433448795889e+04};

	double a2[5] = { 1.046366195300779895e+02,  2.980727727381642723e+03,
                         2.570418404044668245e+04,  5.291014161889741749e+04,
                         9.497228391199055149e+03};
	double b2[5] = {-1.249999999999999534e-01, -1.300633525376058087e+01,
                        -3.651542590150084093e+02, -3.016744074771875522e+03,
                        -5.251679479249748063e+03};
 
	if(x <= 0.0) return 0.0;

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=(((((((a[7]*y+a[6])*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		bj0=fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		f=2.*(bj0*(log(x/2)+EUGAM)+y*fn/fd)/PI;
	}

	else {
/*	Use rational function approximations for P_0 and Q_0 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p0=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q0=fn/(fd*fabs(x));
 
		y=fabs(x);
		f=sqrt(2.0/(PI*y))*(p0*sin(y-PI/4.0)+q0*cos(y-PI/4.0));
/*		bj0=sqrt(2.0/(PI*y))*(p0*cos(y-PI/4.0)-q0*sin(y-PI/4.0)); */
	}

	return f;
}
 


/*	To calculate the Bessel function of second kind of order one
	for real positive argument
	For X<=0 the function is not defined and function will return
	a zero value without any error message or flag.

	Required functions : None
*/
 
#include <math.h>

#define PI 3.14159265358979324
#define EUGAM 0.5772156649015328606    /*	the Euler's constant */

double by1(double x)

{
	double y,fn,fd,f,p1,q1,x2,bj1;
 
/*	Coefficients of rational function approximations */
	double a[7] = {1.156878510067059849e-02,  6.749406787503859073e-05,
                      2.614560837317535451e-07,  7.408815126464007290e-10,
                      1.577637796406197189e-12,  2.432945305413635549e-15,
                      2.257446839754248784e-18};
	double b[8] = {5.000000000000000074e-01, -5.671560744966475746e-02,
                      1.914864631812345532e-03, -2.821407888958585592e-05,
                      2.103168789640803591e-07, -8.322474383730280556e-10,
                      1.678871778708754849e-12, -1.372424374400306547e-15};

	double a0[7] = {1.186963690270342970e-02,  7.123839029323002736e-05,
                        2.848196759538669928e-07,  8.365303089083305885e-10,
                        1.857096246589275584e-12,  3.012506935004947491e-15,
                        2.996951174746838817e-18};
	double b0[8] = {2.500000000000000299e-01, -7.515759077432437273e-02,
                        3.430771992327672576e-03, -6.022315614557372919e-05,
                        5.067136874996839630e-07, -2.197514674456554803e-09,
                        4.768619679411702014e-12, -4.139491442515065355e-15};

	double a1[5] = { 8.659888261699365129e+01,
                         1.932665751369749084e+03,  1.172714583536277145e+04,
                         1.256737699073784218e+04, -6.147124347503755010e+02};
	double b1[5] = { 1.000000000000000011e+00,
                         8.671607011699346720e+01,  1.942669862370300601e+03,
                         1.194181952104744095e+04,  1.371467864585746530e+04};

	double a2[5] = { 1.021472573795463627e+02,  2.807865400111916226e+03,
                         2.280402060738415865e+04,  4.121116954504273053e+04,
                         3.501974669280301705e+03};
	double b2[5] = { 3.749999999999999461e-01,  3.820268245483084309e+01,
                         1.042753017477090289e+03,  8.289951986135169400e+03,
                         1.371889615877945967e+04};
 
	if(x<=0.0) return 0.0;

	if(x<8.0) {
/*	Use rational function approximation */
		y=x*x;
		fn=((((((b[7]*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((a[6]*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		bj1=x*fn/fd;

		fn=((((((b0[7]*y+b0[6])*y+b0[5])*y+b0[4])*y+b0[3])*y+b0[2])*y+b0[1])*y+b0[0];
		fd=((((((a0[6]*y+a0[5])*y+a0[4])*y+a0[3])*y+a0[2])*y+a0[1])*y+a0[0])*y+1;
		f=2.*(bj1*(log(x/2)+EUGAM)-1.0/x-x*fn/fd)/PI;
	}

	else {
/*	Use rational function approximations for P_1 and Q_1 */
		y=1./(x*x);
		fn=(((b1[4]*y+b1[3])*y+b1[2])*y+b1[1])*y+b1[0];
		fd=((((a1[4]*y+a1[3])*y+a1[2])*y+a1[1])*y+a1[0])*y+1;
		p1=fn/fd;

		fn=(((b2[4]*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((a2[4]*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		q1=fn/(fd*fabs(x));
 
		y=fabs(x);
		x2=y-0.75*PI;
		f=sqrt(2.0/(PI*y))*(p1*sin(x2)+q1*cos(x2));
/*		bj1=sqrt(2.0/(PI*y))*(p1*cos(x2)-q1*sin(x2)); */
	}

	return f;
}



/*	To calculate the Bessel function of second kind of integral order
		for real argument

	N : (input) Order of Bessel function required, N may be negative
		or positive. For N=0,1 use BY0 and BY1 respectively.
	X : (input) Argument at which the value is required
	BY : (output) Array of length at least abs(N)+1
		which will contain the value of Bessel function of order
		0,1,...,ABS(N). BY[I] will contain Bessel function of order I
		or -I (if N<0).

	Required functions : BY0, BY1
*/

#include <math.h>

double by0(double x);
double by1(double x);

void byn(int n, double x, double by[])

{
	int i;
 
	by[0]=by0(x);
	by[1]=by1(x);
	if(x<= 0.0) return;
/*	Use the recurrence relation in Forward direction */
	for(i=2; i<=abs(n); ++i) by[i]=2*(i-1)*by[i-1]/x-by[i-2];

	if(n<0) {
		for(i=1; i<=abs(n); i +=2) by[i]= -by[i];
	}
	return;
}
 


/*	To calculate the spherical Bessel function of integral order
	(j_n(x)=Sqrt(PI/(2x))*J_{n+1/2}(x)) for a real argument

	N : (input) Order of Bessel function required, N may be negative
		or positive. 
	XB : (input) Argument at which the value is required
	BJ : (output) Array of length at least 
		ABS(N)+16+MAX(25,5*SQRT(N))
		which will contain the value of Bessel function of order
		0,1,...,ABS(N). BJ[I] will contain Bessel function of order I
		or -I (if N<0)
		Remaining elements of array are used as scratch space

	Required functions : None
*/

#include <math.h>
#include <stdio.h>

void sphbjn(int n, double xb, double bj[])

{
	int i,j,n1;
/*	REPS should be less than the machine accuracy */
	double s,x,xa,t,t0, reps=1.e-17; 

	x=fabs(xb);
	if(xb==0.0) {
		bj[0]=1.0;
		for(i=1; i<=abs(n); ++i) bj[i]=0.0;
		return;
	}
	if(n>0) {
		if(n<x || n<2) {
/*	Use the recurrence relation in the forward direction  */
 
			bj[0]=sin(x)/x;
			bj[1]=sin(x)/(x*x)-cos(x)/x;
			for(i=2; i<=n; ++i) bj[i]=(2*i-1)*bj[i-1]/x-bj[i-2];
		}

		else if(x<=4.0) {
/*	Use series expansion to calculate  bj[n], bj[n-1] */
			xa=x*x/4.0;
			t0=1.0;
			for(i=1; i<=n; ++i) {
				t0=t0*x/(2.0*i+1);
				if(i>=n-1) {
					t=t0; s=t0;
					for(j=1; j<51; ++j) {
						t=-t*xa/(j*(j+i+0.5));
						s=s+t;
						if(fabs(t)<fabs(s*reps)) break;
					}
					bj[i]=s;
				}
			}
			for(i=n-1; i>=1; --i) bj[i-1]=(2.*i+1)*bj[i]/x-bj[i+1];
		}

		else {
/*	Use the recurrence relation in the backward direction */
			n1=5.*sqrt((double) n); if(n1<25) n1=25;
			if(x<n/6.0) n1=n1/log(n*0.5/x);
			n1=n1+n;
			if(n<x+15.0) n1=n1+15;
			if(n1-n<15) n1=n+15;
			bj[n1]=0.0;
			bj[n1-1]=1.0;
			for(i=n1-1; i>=1; --i) bj[i-1]=(2.0*i+1)*bj[i]/x-bj[i+1];

			s=bj[0]*x/sin(x);
			for(i=0; i<=n; ++i) bj[i]=bj[i]/s;

/*	If ABS(bj[n+2])<1./reps, then the required accuracy may not be achieved */
			if(fabs(bj[n+2]*reps)<1.0) printf(" sphbjn failed at n = %d  x = %e s = %e \n",n,x,bj[n+2]);
		}
	}

	else {
/*	For negative N use the recurrence relation in the forward direction  */
		bj[0]=sin(x)/x;
		bj[1]=cos(x)/x;
		for(i=2; i<=abs(n); ++i) bj[i]=(3-2*i)*bj[i-1]/x-bj[i-2];
	}

	if(xb<0.0) {
		for(i=1; i<=abs(n); i +=2) bj[i]=-bj[i];
	}
	return;
}
