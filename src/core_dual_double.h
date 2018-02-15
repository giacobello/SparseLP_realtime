# 1 "core_dual.h"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "core_dual.h"



# 1 "ftype.h" 1
# 178 "ftype.h"
extern "C" void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy);
extern "C" void dcopy_(const int *n, const double *x, const int *incx, double *y, const int *incy);
extern "C" double ddot_(const int *n, const double *x, const int *incx, const double *y, const int *incy);
extern "C" void dscal_(const int *n, const double *a, double *x, const int *incx);
extern "C" double dnrm2_(const int *n, const double *x, const int *incx);
extern "C" double dasum_(const int *n, const double *x, const int *incx);
extern "C" int idamax_(const int *n, const double *x, const int *incx);

extern "C" void dgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const double *alpha, const double *a, const int *lda, const double *b, const int *ldb,
           const double *beta, double *c, const int *ldc);

extern "C" void dgemv_(const char *trans, const int *m, const int *n, const double *alpha,
           const double *a, const int *lda, const double *x, const int *incx,
           const double *beta, double *y, const int *incy);

extern "C" void dgesdd_( const char* jobz, const int* m, const int* n, double* a,
             const int* lda, double* s, double* u, const int* ldu,
             double* vt, const int* ldvt, double* work,
             const int* lwork, int* iwork, int* info );

extern "C" void dpotrf_( const char* uplo, const int* n, double* a, const int* lda,
             int* info );

extern "C" void dpotrs_( const char* uplo, const int* n, const int* nrhs,
             const double* a, const int* lda, double* b,
             const int* ldb, int* info );
# 5 "core_dual.h" 2

typedef struct{
 double *y;
 double *Y;
 int m;
 int n;
 double gamma;
} isignal_double;

typedef struct{

 double epsilon;
 int verbose;
} settings_double;

typedef struct{
 double *a;
 int *kp;
 int *status;
} solution_double;

typedef struct{
 double *c;
 double *xk;
 double *sk;
 double *lk;
 double *b;
 double *rb;
 double *rc;
 double *tn;
 double *skm1;
 double *xkskm1;
 double *d1;
 double *d2;
 double *T;
 double *Kp;
 double *Kpp;
 double *tN;
 double *tN2;
 double *ra;
 double *tm;
 double *dxa;
 double *dsa;
 double *rs;
 double *dss;
 double *dxs;

} variables_double;


int read_problem(char *s, isignal_double *sig);

void free_problem( isignal_double sig);

void cpdip_slp_core( isignal_double sig, settings_double set, variables_double var, solution_double sol);

void multiply_A( isignal_double sig, double *x, double *y, double *tn);

void multiply_At( isignal_double sig, double *y, double *x);

void init_variables( isignal_double sig, settings_double set, variables_double *var, solution_double *sol);

void reset_problem( isignal_double sig, variables_double *var);

void reset_variables( isignal_double sig, variables_double *var);

void free_variables( variables_double var, solution_double sol);

void init_problem( int m, int n, isignal_double *sig);
