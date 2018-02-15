# 1 "tools.h"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "tools.h"



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
# 5 "tools.h" 2

void print_vector( int size, double *a, char *s);

void print_matrix( int size_m, int size_n, double *X, char *s);

double tools_stepsize( int size, double *dx, double *x);

double tools_stepsize_scale( int size, double scale, double *dx, double *x);

void tools_axpy( int size, double alpha, double *x, double *y);

void tools_copy( int size, double *x, double *y);

void tools_scal( int size, double alpha, double *x);

double tools_asum( int size, double *x);

double tools_sum( int size, double *x);

void tools_add( int size, double alpha, double *x);

void tools_set( int size, double alpha, double *x);

double tools_dot( int size, double *x, double *y);

double tools_nrm2( int size, double *x);

void tools_m1( int size, double *x, double *y);

void tools_ea( int size, double *x, double *y);

double tools_inf( int size, double *x);

void tools_hp( int size, double *x, double *y, double *z);

void tools_md( int size_m, int size_n, double *X, double *d, double *Z);

void tools_mdt( int size_m, int size_n, double *X, double *d, double *Z);

void tools_ad( int size, double *d, double *X);

void tools_ads( int size, double alpha, double *X);

void tools_gemvp( char NT, int m, int n, double *A, double *x, double *b);

void tools_gemv (const char transa, const int m, const int n,
         const double alpha, const double *a, const int lda,
         const double *x, const int incx, const double beta,
         double *y, const int incy);

void tools_gemm(const char transa, const char transb, const int m,
        const int n, const int k, const double alpha, const double *a,
        const int lda, const double *b, const int ldb,
        const double beta, double *c, const int ldc);

void tools_gesdd( const char jobz, const int m, const int n, double* a,
             const int lda, double* s, double* u, const int ldu,
             double* vt, const int ldvt, double* work,
         const int lwork, int* iwork, int *info);

void tools_potrf( const char uplo, const int n, double* a, const int lda, int* info);

void tools_potrs( const char uplo, const int n, const int nrhs,
             const double* a, const int lda, double* b,
         const int ldb, int* info);
