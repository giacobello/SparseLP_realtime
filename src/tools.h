#ifndef __TOOLS_H__
#define __TOOLS_H__

#include "ftype.h"

void print_vector( int size, FTYPE *a, char *s);

void print_matrix( int size_m, int size_n, FTYPE *X, char *s);

FTYPE tools_stepsize( INT size, FTYPE *dx, FTYPE *x);

FTYPE tools_stepsize_scale( INT size, FTYPE scale, FTYPE *dx, FTYPE *x);

void tools_axpy( INT size, FTYPE alpha, FTYPE *x, FTYPE *y);

void tools_copy( INT size, FTYPE *x, FTYPE *y);

void tools_scal( INT size, FTYPE alpha, FTYPE *x);

FTYPE tools_asum( INT size, FTYPE *x);

FTYPE tools_sum( INT size, FTYPE *x);

void tools_add( INT size, FTYPE alpha, FTYPE *x);

void tools_set( INT size, FTYPE alpha, FTYPE *x);

FTYPE tools_dot( INT size, FTYPE *x, FTYPE *y);

FTYPE tools_nrm2( INT size, FTYPE *x);

void tools_m1( INT size, FTYPE *x, FTYPE *y);

void tools_ea( INT size, FTYPE *x, FTYPE *y);

FTYPE tools_inf( INT size, FTYPE *x);

void tools_hp( INT size, FTYPE *x, FTYPE *y, FTYPE *z);

void tools_md( INT size_m, INT size_n, FTYPE *X, FTYPE *d, FTYPE *Z);

void tools_mdt( INT size_m, INT size_n, FTYPE *X, FTYPE *d, FTYPE *Z);

void tools_ad( INT size, FTYPE *d, FTYPE *X);

void tools_ads( INT size, FTYPE alpha, FTYPE *X);

void tools_gemvp( char NT, INT m, INT n, FTYPE *A, FTYPE *x, FTYPE *b);

void tools_gemv (CONST char transa, CONST INT m, CONST INT n, 
								 CONST FTYPE alpha, CONST FTYPE *a, CONST INT lda, 
								 CONST FTYPE *x, CONST INT incx, CONST FTYPE beta, 
								 FTYPE *y, CONST INT incy);

void tools_gemm(CONST char transa, CONST char transb, CONST INT m, 
								CONST INT n, CONST INT k,	CONST FTYPE alpha, CONST FTYPE *a, 
								CONST INT lda, CONST FTYPE *b, CONST INT ldb,
								CONST FTYPE beta, FTYPE *c, CONST INT ldc);

void tools_gesdd( CONST char jobz, CONST INT m, CONST INT n, FTYPE* a, 
             CONST INT lda, FTYPE* s, FTYPE* u, CONST INT ldu,
             FTYPE* vt, CONST INT ldvt, FTYPE* work,
									CONST INT lwork, INTA* iwork, INTA *info);

void tools_potrf( CONST char uplo, CONST INT n, FTYPE* a, CONST INT lda, INTA* info);

void tools_potrs( CONST char uplo, CONST INT n, CONST INT nrhs, 
             CONST FTYPE* a, CONST INT lda, FTYPE* b,
									CONST INT ldb, INTA* info);


#endif
