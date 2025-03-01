/*	To calculate the Fermi integral for k=5/2 for real argument

	Required functions : None
*/

#include <math.h>

double ferm25(double x)

{
	double y,f,fn,fd;

/*	Coefficients of rational function approximations */
	double a[9] = { 2.571894150057574794e+00,  2.583310822362159194e+00,
                        1.297464965411635611e+00,  3.478951634872978662e-01,
                        4.953536822731509544e-02,  3.534319266081264567e-03,
                        1.114660054008188913e-04,  1.196397340523770481e-06,
                        2.108653816687462637e-09};
	double b[9] = { 3.323350970447843045e+00,  8.253561418549546218e+00,
                        7.900830547729824064e+00,  3.709901839066293724e+00,
                        9.037486960062833711e-01,  1.118605468265125269e-01,
                        6.453871643985184629e-03,  1.439223289367074978e-04,
                        7.851354785439239006e-07};
 
	double a1[9] = { 3.143350423552492924e-01,  1.185925514982136793e-01,
                         1.507118078811659113e-02,  1.570482413916096774e-03,
                        -5.595169056590362619e-05, -3.021530267629797655e-06,
                         2.935920822791014998e-08, -3.139429619713453590e-10,
                         1.862374323772158779e-12};
	double b1[10] = { 3.082574389659515354e+00,  3.851032874650932823e+00,
                          2.542778697192629395e+00,  1.123138726193956412e+00,
                          3.568215902085525392e-01,  8.321216695003629528e-02,
                          1.368008186484687443e-02,  1.461373661943904411e-03,
                          4.840746803030768445e-05, -6.684259889599044773e-06};
 
	double a2[7] = {-1.010224365188074094e+02,  6.300193167257369564e+04,
                        -2.611825140904720814e+06,  6.096057209422390797e+08,
                        -1.021314648072490935e+10,  7.043216447315414049e+11,
                        -6.058065587728390704e+11};
	double b2[8] = { 2.857142857142856463e-01, -2.475121812396447356e+01,
                         1.758866516117776321e+04, -4.875066917911561303e+05,
                         1.636558505721404325e+08, -4.202416018195836072e+08,
                         1.613941002540396296e+11,  2.688699848072963231e+12};
 
	if(x<2.0) {
		y=exp(x);
		fn=(((((((b[8]*y+b[7])*y+b[6])*y+b[5])*y+b[4])*y+b[3])*y+b[2])*y+b[1])*y+b[0];
		fd=((((((((a[8]*y+a[7])*y+a[6])*y+a[5])*y+a[4])*y+a[3])*y+a[2])*y+a[1])*y+a[0])*y+1;
		f=y*fn/fd;
	}

	else if(x<10.0) {
		fn=((((((((b1[9]*x+b1[8])*x+b1[7])*x+b1[6])*x+b1[5])*x+b1[4])*x+b1[3])*x+b1[2])*x+b1[1])*x+b1[0];
		fd=((((((((a1[8]*x+a1[7])*x+a1[6])*x+a1[5])*x+a1[4])*x+a1[3])*x+a1[2])*x+a1[1])*x+a1[0])*x+1;
		f=fn/fd;
	}
	else {
		y=1.0/(x*x);
		fn=((((((b2[7]*y+b2[6])*y+b2[5])*y+b2[4])*y+b2[3])*y+b2[2])*y+b2[1])*y+b2[0];
		fd=((((((a2[6]*y+a2[5])*y+a2[4])*y+a2[3])*y+a2[2])*y+a2[1])*y+a2[0])*y+1;
		f=pow(x,3.5)*fn/fd;
	}

	return f;
}
