# 1 "core_primal.h"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "core_primal.h"



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
# 5 "core_primal.h" 2


typedef struct{
 float *y;
 float *Y;
 int m;
 int n;
 float gamma;
} isignal_single;

typedef struct{

 float epsilon;
 int verbose;
} settings_single;

typedef struct{
 float *a;
 int *kp;
 int *status;
} solution_single;

typedef struct{
 float *uk;
 float *vk;
 float *g;
 float *gs;
 float *b;
 float *sk1;
 float *sk2;
 float *zk1;
 float *zk2;
 float *rh1;
 float *rh2;
 float *p1;
 float *p2;
 float *d1;
 float *d2;
 float *sk1m1;
 float *sk2m1;
 float *rs1;
 float *rs2;
 float *r1;
 float *r2;
 float *rc1;
 float *rc2;
 float *rhs;
 float *tM;
 float *tM2;
 float *d;
 float *d1pd2m1;
 float *zk1m1;
 float *zk2m1;
 float *dza1;
 float *dza2;
 float *dsa1;
 float *dsa2;
 float *dzs1;
 float *dzs2;
 float *dss1;
 float *dss2;
 float *dv;
 float *Kp;
 float *Kpp;
 float *T;

 int L;
 float *work;
 int *iwork;
 float *s;
 float *U;
 float *Vt;
 float *t1;
 float *t2;
 float *t3;

} variables_single;

int read_problem(char *s, isignal_single *sig);

void free_problem( isignal_single sig);

void cpdip_slp_core( isignal_single sig, settings_single set, variables_single var, solution_single sol);

void multiply_A( isignal_single sig, float *x, float *y);

void multiply_At( isignal_single sig, float *y, float *x);

void init_variables( isignal_single sig, settings_single set, variables_single *var, solution_single *sol);

void reset_problem( isignal_single sig, variables_single *var);

void reset_variables( isignal_single sig, variables_single *var);

void free_variables( variables_single var, solution_single sol);

void init_problem( int m, int n, isignal_single *sig);
