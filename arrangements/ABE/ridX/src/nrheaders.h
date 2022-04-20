#define NRANSI
#define FREE_ARG char*
#define NR_END 1

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))
static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

int *ivector_ridX(long nl, long nh);
double *dvector_ridX(long nl, long nh);
double **dmatrix_ridX(long nrl, long nrh, long ncl, long nch);
unsigned long *lvector_ridX(long nl, long nh);
void free_dvector_ridX(double *v, long nl, long nh);
void free_dmatrix_ridX(double **m, long nrl, long nrh, long ncl, long nch);
void free_lvector_ridX(unsigned long *v, long nl, long nh);
void nrerror_ridX(char error_text[]);
void banbks_ridX(double **a, unsigned long n, int m1, int m2, double **al,
        unsigned long indx[], double b[]);
void bandec_ridX(double **a, unsigned long n, int m1, int m2, double **al,
        unsigned long indx[], double *d);
void banmul_ridX(double **a, unsigned long n, int m1, int m2, double x[], double b[]);
void lubksb_ridX(double **a, int n, int *indx, double b[]);
void ludcmp_ridX(double **a, int n, int *indx, double *d);

void dsvdcmp_ridX(double **a, int m, int n, double w[], double **v);
void dsvbksb_ridX(double **u, double w[], double **v, int m, int n, double b[], double x[]);
double dpythag_ridX(double a, double b);
