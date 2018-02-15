#ifndef __CORE_PRIMAL_H__
#define __CORE_PRIMAL_H__

#include "ftype.h"


typedef struct{
	FTYPE *y;
	FTYPE *Y;
	INT m;
	INT n;
	FTYPE gamma;
} FFTYPE(isignal);

typedef struct{

	FTYPE epsilon;
	INT verbose;
} FFTYPE(settings);

typedef struct{
	FTYPE *a;
	INT *kp;
	INT *status;
} FFTYPE(solution);

typedef struct{
	FTYPE *uk;
	FTYPE *vk;
	FTYPE *g;
	FTYPE *gs;
	FTYPE *b;
	FTYPE *sk1;
	FTYPE *sk2;
	FTYPE *zk1;
	FTYPE *zk2;
	FTYPE *rh1;
	FTYPE *rh2;
	FTYPE *p1;
	FTYPE *p2;
	FTYPE *d1;
	FTYPE *d2;
	FTYPE *sk1m1;
	FTYPE *sk2m1;
	FTYPE *rs1;
	FTYPE *rs2;
	FTYPE *r1;
	FTYPE *r2;
	FTYPE *rc1;
	FTYPE *rc2;
	FTYPE *rhs;
	FTYPE *tM;
	FTYPE *tM2;
	FTYPE *d;
	FTYPE *d1pd2m1;
	FTYPE *zk1m1;
	FTYPE *zk2m1;
	FTYPE *dza1;
	FTYPE *dza2;
	FTYPE *dsa1;
	FTYPE *dsa2;
	FTYPE *dzs1;
	FTYPE *dzs2;
	FTYPE *dss1;
	FTYPE *dss2;
	FTYPE *dv;
	FTYPE *Kp;
	FTYPE *Kpp;
	FTYPE *T;

	INT L;
	FTYPE *work;
	INTA *iwork;
	FTYPE *s;
	FTYPE *U;
	FTYPE *Vt;
	FTYPE *t1;
	FTYPE *t2;
	FTYPE *t3;

} FFTYPE(variables);

int read_problem(char *s, FFTYPE(isignal) *sig);

void free_problem( FFTYPE(isignal) sig);

void cpdip_slp_core( FFTYPE(isignal) sig, FFTYPE(settings) set,	FFTYPE(variables) var, FFTYPE(solution) sol);

void multiply_A( FFTYPE(isignal) sig, FTYPE *x, FTYPE *y);

void multiply_At( FFTYPE(isignal) sig, FTYPE *y, FTYPE *x);

void init_variables( FFTYPE(isignal) sig, FFTYPE(settings) set, FFTYPE(variables) *var, FFTYPE(solution) *sol);

void reset_problem( FFTYPE(isignal) sig, FFTYPE(variables) *var);

void reset_variables( FFTYPE(isignal) sig, FFTYPE(variables) *var);

void free_variables( FFTYPE(variables) var, FFTYPE(solution) sol);

void init_problem( INT m, INT n, FFTYPE(isignal) *sig);

#endif
