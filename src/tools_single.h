# 1 "tools.h"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "tools.h"



# 1 "ftype.h" 1
# 84 "ftype.h"
extern "C" void saxpy_(const int *n, const float *alpha, const float *x, const int *incx, float *y, const int *incy);
extern "C" void scopy_(const int *n, const float *x, const int *incx, float *y, const int *incy);
extern "C" float sdot_(const int *n, const float *x, const int *incx, const float *y, const int *incy);
extern "C" void sscal_(const int *n, const float *a, float *x, const int *incx);
extern "C" float snrm2_(const int *n, const float *x, const int *incx);
extern "C" float sasum_(const int *n, const float *x, const int *incx);
extern "C" int isamax_(const int *n, const float *x, const int *incx);

extern "C" void sgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k,
           const float *alpha, const float *a, const int *lda, const float *b, const int *ldb,
           const float *beta, float *c, const int *ldc);

extern "C" void sgemv_(const char *trans, const int *m, const int *n, const float *alpha,
           const float *a, const int *lda, const float *x, const int *incx,
           const float *beta, float *y, const int *incy);

extern "C" void sgesdd_( const char* jobz, const int* m, const int* n, float* a,
             const int* lda, float* s, float* u, const int* ldu,
             float* vt, const int* ldvt, float* work,
             const int* lwork, int* iwork, int* info );

extern "C" void spotrf_( const char* uplo, const int* n, float* a, const int* lda,
             int* info );

extern "C" void spotrs_( const char* uplo, const int* n, const int* nrhs,
             const float* a, const int* lda, float* b,
             const int* ldb, int* info );
# 5 "tools.h" 2

void print_vector( int size, float *a, char *s);

void print_matrix( int size_m, int size_n, float *X, char *s);

float tools_stepsize( int size, float *dx, float *x);

float tools_stepsize_scale( int size, float scale, float *dx, float *x);

void tools_axpy( int size, float alpha, float *x, float *y);

void tools_copy( int size, float *x, float *y);

void tools_scal( int size, float alpha, float *x);

float tools_asum( int size, float *x);

float tools_sum( int size, float *x);

void tools_add( int size, float alpha, float *x);

void tools_set( int size, float alpha, float *x);

float tools_dot( int size, float *x, float *y);

float tools_nrm2( int size, float *x);

void tools_m1( int size, float *x, float *y);

void tools_ea( int size, float *x, float *y);

float tools_inf( int size, float *x);

void tools_hp( int size, float *x, float *y, float *z);

void tools_md( int size_m, int size_n, float *X, float *d, float *Z);

void tools_mdt( int size_m, int size_n, float *X, float *d, float *Z);

void tools_ad( int size, float *d, float *X);

void tools_ads( int size, float alpha, float *X);

void tools_gemvp( char NT, int m, int n, float *A, float *x, float *b);

void tools_gemv (const char transa, const int m, const int n,
         const float alpha, const float *a, const int lda,
         const float *x, const int incx, const float beta,
         float *y, const int incy);

void tools_gemm(const char transa, const char transb, const int m,
        const int n, const int k, const float alpha, const float *a,
        const int lda, const float *b, const int ldb,
        const float beta, float *c, const int ldc);

void tools_gesdd( const char jobz, const int m, const int n, float* a,
             const int lda, float* s, float* u, const int ldu,
             float* vt, const int ldvt, float* work,
         const int lwork, int* iwork, int *info);

void tools_potrf( const char uplo, const int n, float* a, const int lda, int* info);

void tools_potrs( const char uplo, const int n, const int nrhs,
             const float* a, const int lda, float* b,
         const int ldb, int* info);
