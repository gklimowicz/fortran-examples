/*	header file for all functions included in Appendix C of
	Numerical Methods for Scientists and Engineers, Third Edition,
	by H. M. Antia
*/


int adams(int n, double *y, double *dy, void dif(double , int , double * , double * ),
        double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
        int ijm3, int ijm4, double *wk);


int adi(double x0, double xn, double y0, double yn, int *kn, int nx, int ny,
        double x[], double y[], double *u, int iu,
        void cof(double , double , double * , double * , double * , double * ,
                double * , double * ),
        void bc(int , double , double , double * , double * , double * ),
        double *el, double *eu, double aeps, int *nit);


int adm(double *t, double dt, double x0, double xn, double y0, double yn,
        int nt, int nx, int ny, double x[], double y[], double *u, int iu,
        void cof(double , double , double , double * , double * , double * ,
                double * , double * , double * ),
        double  bc(int , double , double , double ),
        double fic(double , double , double ), int *iflg);


int adpint(double *rint, double xl, double xu, double reps, double aeps,
        double *dif, double (*f) (double ), int *npt, int *nmax);


int balanc(double *a, int n, int ia, double b, int *low, int *igh, double d[]);


int balbak(int n, int low, int igh, double *cz, int m, int iz, double d[]);


int balbak_l(int n, int low, int igh, double *cz, int m, int iz, double d[]);


double betai(double a, double b, double x);

double betap(double a, double b, double x);

double betcon(double a, double b, double x);

double betcon1(double a, double b, double x);

double betcon(double a, double b, double x);

int bfgs(int n, double x[], double *f, double g[], double *h, int *num,
        double reps, double aeps, void fcn(int , double * , double * , double * ));

double bi0(double x);


double bi1(double x);


void bin(int n, double xb, double bi[]);


int bisect(double *x, double *xl, double *xu, int nb, double (*f) (double ));


double bj0(double x);


double bj1(double x);


void bjn(int n, double xb, double bj[]);


void bjy0(double xb, double *bj0, double *by0);


void bjy1(double xb, double *bj1, double *by1);


double bk0(double x);


double bk1(double x);


void bkn(int n, double x, double bk[]);


int brackm(double *a, double *b, double *x, double (*f) (double ));


int brent(double *a, double *b, double *x, double reps, double aeps,
        double (*f) (double ));


int brentm(double *a, double *b, double *x, double *fx, double reps,
        double aeps, double (*f) (double ));


int broydn(void fcn(int , double * , double * ), int np, double x[],
        double f[], double reps, double aeps);


double bspev2(int nx, int ny, double *x, double *y, int k, int nderiv,
        double *wt, int iw, double x0, double y0, double *dfx, double *dfy,
        double *dfxx, double *dfxy, double *dfyy, int *ier);


double bspevl(int n, double *x, int k, int nderiv, double wt[], double x0,
                double *df, double *ddf, int *ier);


double bspevn(int n, int nk[], double *x, int nxd, int k, double *wt,
        double x0[], int *ier);


double bspevn1(int n, int nk[], double *x, int nxd, int k, double *wt,
        double x0[], double df[], int *ier);


double bspevn2(int n, int nk[], double *x, int nxd, int k, double *wt,
        double x0[], double df[], double *ddf, int *ier);


int bspfit(int n, double x[], double f[], double ef[], int k, double *a,
	int la, double *v, int iv, double sigma[], double c[], double xf[],
	int no, double y[], int *iflg, double reps, double rlm, int ide,
	double *chisq, double *cov);


int bspfit2(int nx, int ny, double x[], double y[], double *f, int k,
        double *ax, double *ay, int la, double *vx, double *vy, int iv,
        double sigmax[], double sigmay[], double *c, double xf[], double yf[],
        int mx, int my, double *fy, double reps, double rlm, int ide, double *chisq);


int bspfitn(int n, int nk[], double *x, double *f, int k, double *a, int la,
        double *v, int iv, double *sigma, double *c, double *xf, int mk[],
        double *fy, double reps, double rlm, double ide, double *chisq);


int bspfitw2(int nx, int ny, double x[], double y[], double *f, double *ef,
        int k, double *a, int la, double *v, int iv, double sigma[], double *c,
        int ic, double xf[], double yf[], int mx, int my, double *fy,
        double reps, double rlm, int ide, double *chisq);


int bspfitwn(int n, int nk[], double *x, double *f, double *ef, int k,
        double *a, int la, double *v, int iv, double sigma[], double *c,
        double *xf, int mk[], double *fy, double reps, double rlm, int ide, double *chisq);


int bspint(int n, double x[], double f[], int k, double *a, int la, double c[],
        double *xf, int *no, int *iflg, int inc[]);


int bspint2(int nx, int ny, double x[], double y[], double *f, int k, double *ax,
        double *ay, int la, double *c, double *xf, double *yf, int *mx, int *my,
        int *iflg, int intx[], int inty[]);


int bspintn(int n, int nk[], double *x, int nxd, double *f, int k, double *ax,
        double *c, double *xf, int mk[], int *intx);


int bsplin(double *x, int nx, int k, double xb, int nderiv, double b[],
                double db[], double ddb[], int *left);


int bspode(int nk, int k, int m, int ml, double par[], double *x,
        double *a, double t[], int n, double tx[],
        void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
        void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
        int iflag, double reps);


double bspqd(int n, double *x, int k, double wt[], double xl, double xu, int *ier);


double bspqd2(int nx, int ny, double *x, double *y, int k, double *wt,
        int iw, double xl, double xu, double yl, double yu, int *ier);


double bspqdn(int n, int nk[], double *x, int nxd, int k, double *wt,
        double xl[], double xu[], int *ier);


double by0(double x);


double by1(double x);


void byn(int n, double x, double by[]);


double cabs(double cx[2]);


double cassum(double (*term)(int ), int n);


double cassum_a(double term[ ], int n);


int cauchy(double *ri, double a, double b, double c, double reps, double aeps,
        double *dif, double (*f) (double ), double (*funp) (double ), int *npt);


void cdiv(double c1[2], double c2[2], double cr[2]);


int chebap(int m, int k, double a[], double c[]);


int chebcf(int n, double c[], double p[], int iflg);


int chebex(int n, double c[], double fun(double ));


int cholsk(int n, int num, double * a, double *x, double *det, int nd,
        int *iflg);


int crank(double *t, double dt, double x0, double xn, int nt, int nx,
        double x[], double u[],
        void cof(double , double , double * , double * , double * , double * ),
        void bc(double , double , double , double * , double * , double * ,
                double * , double * , double * ),
        double fic(double , double ), int *iflg);


int crout(int n, int num, double *a, double *x, double *det, int *idet,
	int inc[], int lj, int *iflg);


int crouth(int n, int num, double *a, double *b, double *x, double *det,
        int *idet, int inc[], int lj, double reps, int *iflg, double wk[]);


void csqrt(double cx[2], double cr[2]);


int davidm(double *x1, double *x2, double *f2, double *d2f, double reps,
        double aeps, double (*f) (double , double *));


int davidn(void (*fcn) (int ,double * , double * ,double *), int np, double x[],
        double f[], double reps, double aeps);


int davidn_b(void (*fcn) (int ,double * , double * ), int np, double x[],
        double f[], double reps, double aeps);


double dawson(double x);


int dft(int n, double cg[][2], double cf[][2], int iflg);


int divdif(double xb, double x[], double f[], int *nuse, int ntab,
                double fb[], double aeps, double *dfb, double *ddfb);


int divdif0(double xb, double x[], double f[], int *nuse, int ntab,
                double fb[], double aeps, int iflg, int *if1);


double drvt(double a, int id, double hh0, double reps, double aeps,
        double (*f)(double ), int *ier);


int elmhes(double *a, int n, int ia, int low, int igh, int inc[]);


int epsiln(double *ri, double a, double b, double reps, double aeps,
        double *dif, int *n, double (*fun)(double ));


int equids(double a[], double b[], int n, int npt, double (*f) (int , double * ),
        double *s1, double *s2, double reps, double aeps, double *dif, int *np);


double erf(double x);


double erfc(double x);


int euler(int n, int *m1, int *m2, double a0, double reps, double aeps,
        double *dif, int *n1, int *n2, double *sum, double (*term) (int i));


int extp(int n, double y[], double dy[], void dif(double , int , double * , double * ),
        double *h, double *t0, double tn, double reps, int *nstep, int *nmax,
        int iflg);

double fbeta(double x);

int fdm(int n, int m, int ml, double par[], double *x, double *xc, double t[],
        void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
        void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
        int iflag, double reps);


double ferm05(double x);


double ferm15(double x);


double ferm25(double x);


double fermm05(double x);


int fft(int n, double cg[][2], int iflg);


int fftn(int nd, int nn[], double cg[][2], int iflg);


int fftr(int n, double cg[][2], int iflg);


int filon(double *ri, double xl, double xu, double rk, int qsin, double reps,
        double aeps, double *dif, int *n, double (*fun) (double ));


double fln(void fcn(int , double * , double * ), double x, double v[],
        double x0[], int n, int *num);


double flnm(void (*fcn) (int , double * , double * , double * ), double x,
        double *df, double v[], double x0[], int n, int *num);


double fm(double x);


int forw(int np, int nm, double r[], double *rker, int ik, double di[],
        double f[], double fun(double ), int *iflg);


int fred(int m, double a, double b, double wt[], double x[], double f[],
        double fc[], double fg(double ), double fker(double , double ),
        double *ei, int iq, int it, double reps);


int fredco(int n, double a, double b, double f[], double x[], double reps,
        double aeps, int iq, int it);


double funk(double t);



double gamma(double xg);

double gammaln(double xg);


double gammap(double a, double x);


int gaublk(int n, int m, int ml, double *a, int *iflg, double *det,
        int *idet, int *inc, double *x);


int gaubnd(int n, int kb, int num, double *a, double *x, double *det,
                int *idet, int inc[], int lj, int *iflg);


int gaucb1(double *rint, double a, double b, double reps, double aeps,
        double *dif, int *npt, double (*fun) (double));


int gaucb2(double *rint, double a, double b, double reps, double aeps,
        double *dif, int *npt, double (*fun) (double));


int gaucby(double *rint, double a, double b, double reps, double aeps,
        double *dif, int *npt, double (*fun) (double));


int gauelm(int n, int num, double * a, double * x, double * det,
       int inc[], int lj, int * iflg);


int gauher(int n, double w[], double a[]);


int gaujac(int n, double alp, double beta, double w[], double a[]);


int gaulag(double *rint, double a, double *a1, double reps, double aeps,
        double *dif, double (*f) (double ), double (*f2) (double ), int *np);


int gauleg(int n, double w[], double a[]);


int gaulg2(double *rint, double a, double *a1, double reps, double aeps,
        double *dif, double (*f) (double ), double (*f2) (double ), int *np);


int gaulog(double *rint, double a, double aeps, double reps, double *dif,
        double (*f) (double ), int *npt);


int gaus16(double *ri, double a, double b, double *dif, int *n,
        double (*f)(double ));


int gausq(double *rint, double a, double aeps, double reps, double *dif,
        int *npt, double (*fun) (double ) );


int gausq2(double *rint, double a, double *a1, double reps, double aeps,
        double *dif, double (*f) (double ), double (*f2) (double ), int *np);


int gausrc(int n, double w[], double ab[], double *cof, double ri0);


int gauss(double *rint, double a, double b, int *np, double reps, double aeps,
        double *dif, int *npt, double (*fun) (double ));


int gauswt(int n, double w[], double ab[], double (*fmom) (int ), int qgaus);


int gear(int n, double *y, double *dy, void dif(double , int , double * , double * ),
        double h, double t, double reps, int *nstp, int ij, int ijm1, int ijm2,
        int ijm3, int ijm4, int *iflag, double *wk1, double *wk);


int gevp(int n, int m, int ml, double par[], double *x, double *xc,
        double t[], double *e0,
        void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
        void bcs(int , int , double * , double * , double * , double , double , double * , double * ),
        void eqnd(int , int , int , double * , double * , double * , double ),
        void bcsd(int , int , double * , double * ,  double , double ),
        int iflag, double reps, double el, double eu);


int golden(double *a, double *b, double *x, double *fx, double reps,
        double aeps, double (*f) (double ));


int herevp(double *ar, double *ai, int n, int ia, double ei[], double *vr,
        double *vi, int iv, double reps);


int hermit(double *rint, double aeps, double reps, double *dif,
        double (*f) (double ), int *npt);


int hqr(double *h, int nn, int ih, double er[], double ei[], double reps);


int invit(double *a, int m, int ia, double *p, double u[], int iflg,
        double *ei, double *erc, double reps, int *nit);


int invit_l(double *a, int m, int ia, double *p, double u[], int iflg,
        double *ei, double *erc, double reps, int *nit);


int iranbin(double *seed, int n, double p, double c[]);


int iranpoi(double *seed, double rmu, double p[]);


int kronrd(double *ri, double a, double b, double *dif, int *n,
        double (*f)(double ));


int lagitr(int n, double a[], double cxi[2]);


int lagure(double *rint, double a, double aeps, double reps, double *dif,
        double (*f) (double ), int *npt);


int lagurw(int n, double alp, double w[], double a[]);


int lapinv(int n, double t[], double f[], void cfs(double * , double * ),
        double alpha, double *reps);


int lax(int n, double *t, double dt, double x0, double xn, int nt, int nx,
        double x[], double *u, int iu,
        void flux(int , double , double , double * , double * ),
        void bc(int , int , double , double , int *, double * ),
        void fic(int , double , double , double * ), int *iflg);


void lines(double t, int n, double u[], double du[]);


int linfitxy(int n, double *x, double *y, double sigx, double sigy,
	double rho, double *xi, double *yi, double *a, double *b,
	double *chi);


int linl1(int m, int n, double a[], double f[], double *g, int ig,
        double eps, double *esum);


int linmin(double *x1, double *x2, double *f1, double *df1, double reps,
        double aeps, void (*f) (int , double * , double * , double * ),
        double v[], double xi[], int n, int *num);


int linmnf(double *x0, double *x1, double *f0, double reps, double aeps,
        void (*f) (int , double * , double * ), double v[], double xi[],
        int n, int *num);


int linrn(int n, double xb[], double x[], double f[], int np[], double *fb,
        int ndim[], int nxd);


int llsq(int n, int m, int k, double *x, int ix, double f[], double ef[],
	double a[], double *u, double *v, int iu, int iv, double sigma[],
	double y[], void (* phi) (int , double * , double * ), double reps,
	double *chisq, double *cov);


int locate(double xb, double x[], int np);


int matinv(int n, int ia, double *a, double *ai);


int mcarlo(double a[], double b[], int n, int npt, double (*f) (int , double * ),
        double *ri, double reps, double aeps, double *err, int *np);


int minmax(int m, int k, int n, double a[], double x[], double f[],
        double eps, double *emax);


int mstep(int n, double *y, double *dy, void dif(double , int , double * , double * ),
        double *h, double *t0, double tn, double yf[], double reps, int *nstp,
        int *nmax, int *iflg, int ist, double *wk);


void muler2(double cx1[], double cx2[], double cx3[], double reps, double aeps,
        int *ier, double cf[], double cx[], int ix, int nz, double czero[][2], double rmax);


int mulint(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
        double *rint, double reps, double aeps, double *dif, int *num, int *maxpt);


int muller(double cx1[], double cx2[], double cx3[], double reps, double aeps,
        void cf(double * , double * ), int nz, double czero[][2], double rmax);


int nearst(double xb, double x[], int ntab);


int newrap(double x0, double xl, double xu, double *x, double reps,
        double aeps, double (*fun) (double , double * ));


int newton(void fcn(int , double * ,double *, double * ), int np,
        double x[], double f[], double reps, double aeps);


int ngauss(double a[], double b[], int n, int m[], int ind[], double (*f) (int , double * ),
        double *ri, int *num, int *maxpt);


void nllsq(int n, double a[], double *f, double g[]);


void nllsq_f(int n, double a[], double *f);


int nminf(int n, double x[], double *f, int *num, double reps, double aeps,
        void (*fcn) (int , double * , double * ));


int pade(int m, int k, double a[], double c[]);


double pcor(int n, double xx);


void pleg(int l, double x, double p[]);


void plm(int l, int m, double x, double p[]);


double pold(int n, double a[], double x, int nd, double pd[]);


int polev2(int nx, int ny, double *ax, double *ay, int la, double *wt,
        double x0, double y0, double *f, double *dfx, double *dfy,
        double *dfxx, double *dfxy, double *dfyy);


int polevl(int m, double a[], double alp[], double beta[], double x,
        double *f, double *df, double *ddf);


int polevn(int n, int nk[], double *ax, int la, double *wt, double x0[],
        double *f);


int polevn1(int n, int nk[], double *ax, int la, double *wt, double x0[],
        double *f, double df[]);


int polevn2(int n, int nk[], double *ax, int la, double *wt, double x0[],
        double *f, double df[], double *ddf);


int polfit(int n, int m, double x[], double f[], double sig[], double a[],
	double alp[], double beta[], double y[], double h[], double gam[]);


int polfit1(int n, int m, int num, double x[], double *f, double w[],
        double *a, double alp[], double beta[], double gam[]);


int polfit2(int nx, int ny, double x[], double y[], double *f, double *ax,
        double *ay, int la, double *c, int ic, int mx, int my, double *fy,
        double *chisq);


int polfitn(int n, int nk[], double *x, double *f, double *ax, int la,
        double *c, int mk[], double *fy, double *chisq);


int polort(int m, double alp[], double beta[], double x, double f[],
        double df[], double ddf[]);


int poly2(double xb1, double xb2, double x1[], double x2[], double *f, int ndim,
        int n1, int n2, int np1, int np2, double *fb);


int polyl1(int m, int n, double a[], double x[], double f[], double eps,
        double *esum);


int polyr(int n, double a[], double rr[], double ri[], int qrefin);


double ran1(double *seed);


double ranf(int *iseed);


double rangau(double *seed);


int ratnal(double xb, double x[], double f[], int *nuse, int ntab,
        double *fb, double aeps);


int remes(int m, int k, int *n, double xl, double xu, double a[], double x[],
        double f[], double ex[], int *ie, double *emax, double eps, double epsm,
        int iflg);


void rk2(int n, double t, double y0[], double dy0[], double h, double y1[],
        void dif(double , int , double * , double * ));


void rk4(int n, double t, double y0[], double dy0[], double h, double y1[],
        void dif(double , int , double * , double * ));


int rkm(int n, double y[], double dy[], void dif(double , int , double * , double *),
        double *h, double *t0, double tn, double reps, int *nstep, int *nmax);


int rkm_2(int n, double y[], double dy[], void dif(double , int , double * , double *),
        double *h, double *t0, double tn, double reps, int *nstep, int *nmax);


int rls(int nk, double xo[], int k, int nr, double r[], double *rker,
        int ik, double *ac, int nm, int ns, double alp, int ide, double di[],
        double de[], double df[], double f[], double b[], int *iflg,
        double reps, double *chisq, double *sumd, double *a, double *av,
        int iv, double sigma[], int nsim, double fe[]);


double rmk(int m, int k, double a[], double b[], double x);


double rmk1(int m, int k, double a[], double b[], double x);


double rmkd(int m, int k, double a[], double b[], double x, double *df);


double rmkd1(int m, int k, double a[], double b[], double x, double *df);


int rombrg(double *ri, double a, double b, double gi[], double reps,
        double aeps, double *dif, int *n, double (*fun)(double ));


double round(double x, int n, double b);


int search(double rx1, double rx2, double ry1, double ry2, int *nx,
        int *ny, void cfun(double * , double * ));


int secan_2(double x0, double xl, double xu, double *x, double reps,
        double aeps, double (*fun) (double , int * ));


void secani(double x0, double xl, double xu, double *x, double f, int jf,
        double reps,double aeps, int *ier);


int secant(double x0, double xl, double xu, double *x, double reps,
        double aeps, double (*fun) (double ));


int setmat(int n, int m, int ml, double *a, double *bc, double *x,
        double *xc, double t[], double par[],
        void eqn(int , int , int , double * , double * , double * , double * , double * , double ),
        void bcs(int , int , double * , double * , double * , double , double , double * , double * ));

void shsort(double a[], int n);

int simpl1(double *a, int ia, int n, int m, int id[], int iv[], double aeps);


int simplx(double *a, int ia, int n, int m1, int m2, int m3, double x[],
        double *f, double aeps);


int simpx(double *a, int ia, int n, int m, int nv, int qf, int id[],
        int iv[], double aeps);


int simson(double *ri, double xl, double xu, double reps, double aeps,
        double *dif, int *n, double (*fun)(double ));


int smooth(int ntab, double x[], double f[], double c[][3], int np,
        double xp[], double fp[]);


int sor(double x0, double xn, double y0, double yn, int nx, int ny,
        double x[], double y[], double *u, int iu,
        void cof(double , double , double * , double * , double * ,
                double * , double * , double * , double *),
        double bc(int , double , double ), double *omega, double aeps, int *nit);


void sphbjn(int n, double xb, double bj[]);


double sphnd(int n, double x[]);


double splevl(double xb, int n, double x[], double f[], double c[][3],
                double *dfb, double *ddfb, int *ier);


int spline(double x[], double f[], int n, double c[][3]);


int splint(double a, double b, double *sint, double *tint, int n, double x[],
        double f[], double c[][3]);


int strint(double a[], double b[], int n, int *m, int ind[], double (*f) (int , double * ),
        double *rint, double reps, double aeps, double *dif, int *num, int *maxpt);


int stroud(double a[], double b[], int n, int m, int ind[], double (*f) (int , double * ),
        double *ri, int *num, int *maxpt);


int strt4(int n, double *y, double *dy, void dif(double , int , double * , double * ),
        double *h, double *t, double reps, int iflg, double tstep, int *nstp,
        double *wk);


int sturm(double e[], double d[], int n, int m1, int m2, double el[],
        double eu[], int *num, double reps);


int svd(int n, int m, double *a, double *v, double sigma[], int la, int lv);


int svdevl(int n, int m, double *u, double *v, double sigma[], int lu,
        int lv, double b[], double reps);


int tinvit(double e[], double d[], int n, double el, double eu, double *ei,
        double ev[], double reps, int iflg, int *num);


int tql2(double *z, int n, int iz, double d[], double e[], double reps);


int trbak(double *a, int ia, int n, double *z, int iz, int nz);


int tred2(double *a, int n, int ia, double d[], double e[]);


int tridia(double e[], double d[], int n, int m1, int m2, double ei[],
        double eps1, double reps, double *ev, int iv);


int volt(int n, double a, double h, double f[], double x[],
        double fg(double ), double fker(double , double ), int it);


int volt2(int n, double a, double h, double f[], double x[],
        double fg(double ), double fker(double , double , double), double reps);


void ylm(int l, int m, double theta, double phi, double y[]);


void ylm_x(int l, int m, double x, double phi, double y[]);


int zroot(int n, double cx[][2], double czero[][2], int *nz, double reps,
        double aeps, double rmax, void cf(double * , double * ));


int zroot2(int n, double cx[][2], double czero[][2], int *nz, double reps,
        double aeps, double rmax, void cf(double * , double * , int * ));




