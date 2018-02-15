#ifndef __CORE_DUAL_H__
#define __CORE_DUAL_H__

#include "ftype.h"

typedef struct{
	FTYPE *y;
	FTYPE *Y;
	int m;
	int n;
	FTYPE gamma;
} FFTYPE(isignal);

typedef struct{

	FTYPE epsilon;
	int verbose;
} FFTYPE(settings);

typedef struct{
	FTYPE *a;
	int *kp;
	int *status;
} FFTYPE(solution);

typedef struct{
	FTYPE *c;
	FTYPE *xk;
	FTYPE *sk;
	FTYPE *lk;
	FTYPE *b;
	FTYPE *rb;
	FTYPE *rc;
	FTYPE *tn;
	FTYPE *skm1;
	FTYPE *xkskm1;
	FTYPE *d1;
	FTYPE *d2;
	FTYPE *T;
	FTYPE *Kp;
	FTYPE *Kpp;
	FTYPE *tN;
	FTYPE *tN2;
	FTYPE *ra;
	FTYPE *tm;
	FTYPE *dxa;
	FTYPE *dsa;
	FTYPE *rs;
	FTYPE *dss;
	FTYPE *dxs;

} FFTYPE(variables);


int read_problem(char *s, FFTYPE(isignal) *sig);

void free_problem( FFTYPE(isignal) sig);

void cpdip_slp_core( FFTYPE(isignal) sig, FFTYPE(settings) set,	FFTYPE(variables) var, FFTYPE(solution) sol);

void multiply_A( FFTYPE(isignal) sig, FTYPE *x, FTYPE *y, FTYPE *tn);

void multiply_At( FFTYPE(isignal) sig, FTYPE *y, FTYPE *x);

void init_variables( FFTYPE(isignal) sig, FFTYPE(settings) set, FFTYPE(variables) *var, FFTYPE(solution) *sol);

void reset_problem( FFTYPE(isignal) sig, FFTYPE(variables) *var);

void reset_variables( FFTYPE(isignal) sig, FFTYPE(variables) *var);

void free_variables( FFTYPE(variables) var, FFTYPE(solution) sol);

void init_problem( int m, int n, FFTYPE(isignal) *sig);

#endif
