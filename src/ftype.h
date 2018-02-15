#ifndef __FTYPE_H__
#define __FTYPE_H__

#ifdef DOUBLE
#define FTYPE double
#define FNAME double
#endif

#ifdef SINGLE
#define FTYPE float
#define FNAME single
#endif

#define PASTER(x,y) x ## _ ## y
#define EVALUATOR(x,y)  PASTER(x,y)
#define FFTYPE(fun) EVALUATOR(fun, FNAME)

#ifdef MATLAB
#include <mex.h>
/*
#include <blas.h>
#include <lapack.h>

#define AXPY daxpy
#define COPY dcopy
#define DOT ddot 
#define SCAL dscal
#define NRM2 dnrm2
#define ASUM dasum
#define AMAX idamax
#define GEMM dgemm
#define GEMV dgemv
#define GESDD dgesdd
#define POTRF dpotrf
#define POTRS dpotrs
*/

#define INT mwSignedIndex
#define INTA mwSignedIndex

#define CONST 

#else
#define INT int
#define INTA int
#define CONST const 
#endif

#ifdef __cplusplus
#define EXTERN extern "C"
#else
#define EXTERN extern
#endif


#ifdef SINGLE

#if ( defined(_WIN32) || defined(__WIN32__) || defined(PC) || defined(__WINDOWS__) || defined(MKL) )
#define AXPY saxpy
#define COPY scopy
#define DOT sdot 
#define SCAL sscal
#define NRM2 snrm2
#define ASUM sasum
#define AMAX isamax
#define GEMM sgemm
#define GEMV sgemv
#define GESDD sgesdd
#define POTRF spotrf
#define POTRS spotrs

EXTERN void saxpy(CONST INT *n, CONST float *alpha, CONST float *x, CONST INT *incx, float *y, CONST INT *incy);
EXTERN void scopy(CONST INT *n, CONST float *x, CONST INT *incx, float *y, CONST INT *incy);
EXTERN float sdot(CONST  INT *n, CONST float *x, CONST INT *incx, CONST float *y, CONST INT *incy);
EXTERN void sscal(CONST INT *n, CONST float *a, float *x, CONST INT *incx);
EXTERN float snrm2(CONST INT *n, CONST float *x, CONST INT *incx);
EXTERN float sasum(CONST INT *n, CONST float *x, CONST INT *incx);
EXTERN INT isamax(CONST INT *n, CONST float *x, CONST INT *incx);

EXTERN void sgemm(CONST char *transa, CONST char *transb, CONST INT *m, CONST INT *n, CONST INT *k,
           CONST float *alpha, CONST float *a, CONST INT *lda, CONST float *b, CONST INT *ldb,
           CONST float *beta, float *c, CONST INT *ldc);

EXTERN void sgemv(CONST char *trans, CONST INT *m, CONST INT *n, CONST float *alpha,
           CONST float *a, CONST INT *lda, CONST float *x, CONST INT *incx,
           CONST float *beta, float *y, CONST INT *incy);

EXTERN void sgesdd( CONST char* jobz, CONST INT* m, CONST INT* n, float* a, 
             CONST INT* lda, float* s, float* u, CONST INT* ldu,
             float* vt, CONST INT* ldvt, float* work,
             CONST INT* lwork, INTA* iwork, INTA* info );

EXTERN void spotrf( CONST char* uplo, CONST INT* n, float* a, CONST INT* lda, 
             INTA* info );

EXTERN void spotrs( CONST char* uplo, CONST INT* n, CONST INT* nrhs, 
             CONST float* a, CONST INT* lda, float* b,
             CONST INT* ldb, INTA* info );

#else 

#define AXPY saxpy_
#define COPY scopy_
#define DOT sdot_ 
#define SCAL sscal_
#define NRM2 snrm2_
#define ASUM sasum_
#define AMAX isamax_
#define GEMM sgemm_
#define GEMV sgemv_
#define GESDD sgesdd_
#define POTRF spotrf_
#define POTRS spotrs_

EXTERN void saxpy_(CONST INT *n, CONST float *alpha, CONST float *x, CONST INT *incx, float *y, CONST INT *incy);
EXTERN void scopy_(CONST INT *n, CONST float *x, CONST INT *incx, float *y, CONST INT *incy);
EXTERN float sdot_(CONST  INT *n, CONST float *x, CONST INT *incx, CONST float *y, CONST INT *incy);
EXTERN void sscal_(CONST INT *n, CONST float *a, float *x, CONST INT *incx);
EXTERN float snrm2_(CONST INT *n, CONST float *x, CONST INT *incx);
EXTERN float sasum_(CONST INT *n, CONST float *x, CONST INT *incx);
EXTERN INT isamax_(CONST INT *n, CONST float *x, CONST INT *incx);
 
EXTERN void sgemm_(CONST char *transa, CONST char *transb, CONST INT *m, CONST INT *n, CONST INT *k,
           CONST float *alpha, CONST float *a, CONST INT *lda, CONST float *b, CONST INT *ldb,
           CONST float *beta, float *c, CONST INT *ldc);

EXTERN void sgemv_(CONST char *trans, CONST INT *m, CONST INT *n, CONST float *alpha,
           CONST float *a, CONST INT *lda, CONST float *x, CONST INT *incx,
           CONST float *beta, float *y, CONST INT *incy);

EXTERN void sgesdd_( CONST char* jobz, CONST INT* m, CONST INT* n, float* a, 
             CONST INT* lda, float* s, float* u, CONST INT* ldu,
             float* vt, CONST INT* ldvt, float* work,
             CONST INT* lwork, INTA* iwork, INTA* info );

EXTERN void spotrf_( CONST char* uplo, CONST INT* n, float* a, CONST INT* lda, 
             INTA* info );

EXTERN void spotrs_( CONST char* uplo, CONST INT* n, CONST INT* nrhs, 
             CONST float* a, CONST INT* lda, float* b,
             CONST INT* ldb, INTA* info );

#endif

#endif //SINGLE

#ifdef DOUBLE

#if (defined(_WIN32) || defined(__WIN32__) || defined(PC) || defined(__WINDOWS__) || defined(MKL) )

#define AXPY daxpy
#define COPY dcopy
#define DOT ddot 
#define SCAL dscal
#define NRM2 dnrm2
#define ASUM dasum
#define AMAX idamax
#define GEMM dgemm
#define GEMV dgemv
#define GESDD dgesdd
#define POTRF dpotrf
#define POTRS dpotrs


EXTERN void daxpy(CONST INT *n, CONST double *alpha, CONST double *x, CONST INT *incx, double *y, CONST INT *incy);
EXTERN void dcopy(CONST INT *n, CONST double *x, CONST INT *incx, double *y, CONST INT *incy);
EXTERN double ddot(CONST  INT *n, CONST double *x, CONST INT *incx, CONST double *y, CONST INT *incy);
EXTERN void dscal(CONST INT *n, CONST double *a, double *x, CONST INT *incx);
EXTERN double dnrm2(CONST INT *n, CONST double *x, CONST INT *incx);
EXTERN double dasum(CONST INT *n, CONST double *x, CONST INT *incx);
EXTERN INT idamax(CONST INT *n, CONST double *x, CONST INT *incx);

EXTERN void dgemm(CONST char *transa, CONST char *transb, CONST INT *m, CONST INT *n, CONST INT *k,
           CONST double *alpha, CONST double *a, CONST INT *lda, CONST double *b, CONST INT *ldb,
           CONST double *beta, double *c, CONST INT *ldc);

EXTERN void dgemv(CONST char *trans, CONST INT *m, CONST INT *n, CONST double *alpha,
           CONST double *a, CONST INT *lda, CONST double *x, CONST INT *incx,
           CONST double *beta, double *y, CONST INT *incy);

EXTERN void dgesdd( CONST char* jobz, CONST INT* m, CONST INT* n, double* a, 
             CONST INT* lda, double* s, double* u, CONST INT* ldu,
             double* vt, CONST INT* ldvt, double* work,
             CONST INT* lwork, INTA* iwork, INTA* info );

EXTERN void dpotrf( CONST char* uplo, CONST INT* n, double* a, CONST INT* lda, 
             INTA* info );

EXTERN void dpotrs( CONST char* uplo, CONST INT* n, CONST INT* nrhs, 
             CONST double* a, CONST INT* lda, double* b,
             CONST INT* ldb, INTA* info );



#else 
#define AXPY daxpy_
#define COPY dcopy_
#define DOT ddot_ 
#define SCAL dscal_
#define NRM2 dnrm2_
#define ASUM dasum_
#define AMAX idamax_
#define GEMM dgemm_
#define GEMV dgemv_
#define GESDD dgesdd_
#define POTRF dpotrf_
#define POTRS dpotrs_

EXTERN void daxpy_(CONST INT *n, CONST double *alpha, CONST double *x, CONST INT *incx, double *y, CONST INT *incy);
EXTERN void dcopy_(CONST INT *n, CONST double *x, CONST INT *incx, double *y, CONST INT *incy);
EXTERN double ddot_(CONST  INT *n, CONST double *x, CONST INT *incx, CONST double *y, CONST INT *incy);
EXTERN void dscal_(CONST INT *n, CONST double *a, double *x, CONST INT *incx);
EXTERN double dnrm2_(CONST INT *n, CONST double *x, CONST INT *incx);
EXTERN double dasum_(CONST INT *n, CONST double *x, CONST INT *incx);
EXTERN INT idamax_(CONST INT *n, CONST double *x, CONST INT *incx);

EXTERN void dgemm_(CONST char *transa, CONST char *transb, CONST INT *m, CONST INT *n, CONST INT *k,
           CONST double *alpha, CONST double *a, CONST INT *lda, CONST double *b, CONST INT *ldb,
           CONST double *beta, double *c, CONST INT *ldc);

EXTERN void dgemv_(CONST char *trans, CONST INT *m, CONST INT *n, CONST double *alpha,
           CONST double *a, CONST INT *lda, CONST double *x, CONST INT *incx,
           CONST double *beta, double *y, CONST INT *incy);

EXTERN void dgesdd_( CONST char* jobz, CONST INT* m, CONST INT* n, double* a, 
             CONST INT* lda, double* s, double* u, CONST INT* ldu,
             double* vt, CONST INT* ldvt, double* work,
             CONST INT* lwork, INTA* iwork, INTA* info );

EXTERN void dpotrf_( CONST char* uplo, CONST INT* n, double* a, CONST INT* lda, 
             INTA* info );

EXTERN void dpotrs_( CONST char* uplo, CONST INT* n, CONST INT* nrhs, 
             CONST double* a, CONST INT* lda, double* b,
             CONST INT* ldb, INTA* info );

#endif

#endif //DOUBLE

#endif
