#ifndef __FCORE_H__
#define __FCORE_H__

#if defined(SD)
#ifdef PRIMAL

#include "core_primal_single.h"
#include "tools_single.h"

#include "core_primal_double.h"
#include "tools_double.h"

/* copy variables from single to double precision to start from the new position*/
void copy_variables( isignal_double *sig_double, variables_single *var_single, variables_double *var_double){
	int i, M;

	M = sig_double->m + sig_double->n;

	for( i = 0 ; i < M ; ++i ){
		/* copy sk */
		var_double->sk1[i] = (double) var_single->sk1[i];
		var_double->sk2[i] = (double) var_single->sk2[i];

		/* copy zk */
		var_double->zk1[i] = (double) var_single->zk1[i];
		var_double->zk2[i] = (double) var_single->zk2[i];

		/* copy vk */
		var_double->vk[i] = (double) var_single->vk[i];
	}

	/* copy uk */
	for( i = 0 ; i < sig_double->n ; ++i )
		var_double->uk[i] = (double) var_single->uk[i];

	/* set b vector */
	for( i = 0 ; i < sig_double->m ; ++i)
		var_double->b[i] = -sig_double->y[i];
}

#endif /* end primal */
#ifdef DUAL 

#include "core_dual_single.h"
#include "tools_single.h"

#include "core_dual_double.h"
#include "tools_double.h"

/* copy variables from single to double precision to start from the new position*/
void copy_variables( isignal_double *sig_double, variables_single *var_single, variables_double *var_double){
	int i, M;

	M = sig_double->m + sig_double->n;

	/* copy xk */
	for( i = 0 ; i < 2*M ; ++i )
		var_double->xk[i] = (double) var_single->xk[i];

	/* copy sk */
	for( i = 0 ; i < 2*M ; ++i )
		var_double->sk[i] = (double) var_single->sk[i];

	/* copy lk */
	for( i = 0 ; i < sig_double->m ; ++i)
		var_double->lk[i] = (double) var_single->lk[i];

}

#endif 

/* copy the problem from single to double precision*/
void copy_problem( isignal_single *sig_single, isignal_double *sig_double){
	int i;

	sig_double->m = sig_single->m;
	sig_double->n = sig_single->n;
	sig_double->gamma = (double) sig_single->gamma;

	for( i = 0 ; i < sig_single->m ; ++i)
		sig_double->y[i] = (double) sig_single->y[i];

	for( i = 0 ; i < (sig_single->m)*(sig_single->n) ; ++i)
		sig_double->Y[i] = (double) sig_single->Y[i];
}

/* copy the problem from double to single precision*/
void copy_problem( isignal_double *sig_double, isignal_single *sig_single){
	int i;

	sig_single->m = sig_double->m;
	sig_single->n = sig_double->n;
	sig_single->gamma = (float) sig_double->gamma;

	for( i = 0 ; i < sig_single->m ; ++i)
		sig_single->y[i] = (float) sig_double->y[i];

	for( i = 0 ; i < (sig_single->m)*(sig_single->n) ; ++i)
		sig_single->Y[i] = (float) sig_double->Y[i];
}

#else /* else SINGLE && DOUBLE */

#ifdef PRIMAL
#include "core_primal.h"
#include "tools.h"
#endif

#ifdef DUAL
#include "core_dual.h"
#include "tools.h"
#endif

#endif /* end SINGLE && DOUBLE */

#endif /* __FCORE_H__ */
