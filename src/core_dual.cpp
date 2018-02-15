/* COPYRIGHT:
   2012 Tobias L. Jensen

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE?2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
   implied. See the License for the specific language governing
   permissions and limitations under the License.
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "ftype.h"

#if defined(SD) && defined(SINGLE)
#include "core_dual_single.h"
#include "tools_single.h"
#elif defined(SD )&& defined(DOUBLE)
#include "core_dual_double.h"
#include "tools_double.h"
#else
#include "tools.h"
#include "core_dual.h"
#endif

#define SOLVED 1
#define FAILED 0

#define MAXIT 30

#define COUNT_MAX 5

/* Solves the sparse linear prediction problem using the dual algorithm */
void cpdip_slp_core( FFTYPE(isignal) sig, FFTYPE(settings) set,	FFTYPE(variables) var, FFTYPE(solution) sol)
{
	int k,M,N;
	
	INTA info;
	int count;

	FTYPE mu, pinfeas, dinfeas, relgap, cxk;
	FTYPE axa, asa, mu_aff, sigma;
	FTYPE etak, ap, ad;

	M = sig.m;
	N = 2*sig.m+2*sig.n;

	*sol.status = FAILED; /* Set initial */

	if( set.verbose ){
		printf("-- Sparse Linear Prediction (dual) ---- \n");
		printf("Ite.     pinfeas     dinfeas     rel_gap\n");
	}

	/* main iteration */
	for( k = 1 ; k <= MAXIT ; ++k ){
		/* ************************************************** */
		/* 1. compute residual and evaluate stopping criteria */

		/*
		print_vector( N, var.xk, " xk ");
		print_vector( M, var.lk, " lk ");
		print_vector( N, var.sk, " sk ");
		*/

		/* 	rb = A*xk - b; */
		multiply_A( sig, var.xk, var.rb, var.tn); 
		tools_axpy( M, -1.0, var.b, var.rb);

		/* rc = At*lk + sk - c; */
		multiply_At( sig, var.lk, var.rc);
		tools_axpy( N, 1.0, var.sk, var.rc);
		tools_axpy( N, -1.0, var.c, var.rc);

		mu = tools_dot( N, var.xk,var.sk)/N;

		pinfeas = tools_nrm2( M, var.rb)/(1.0+tools_nrm2( M, var.rb));
		dinfeas = tools_nrm2( N, var.rc)/(1.0+tools_nrm2( N, var.c));
		cxk = tools_dot( N, var.c, var.xk);
		relgap = (cxk-tools_dot( M, var.b, var.lk))/
			(1.0+fabs(cxk));

		if( set.verbose ){
			printf("%2d      %.2e     %.2e    %.2e\n",k,pinfeas,dinfeas,relgap);
		}

		/* Stop if sufficiently accurate */
		if( relgap < set.epsilon && pinfeas < set.epsilon && dinfeas < set.epsilon ){
			*sol.status = SOLVED;
			break;
		}

		/* ************************************************** */
		/* 2. calculate the affine scaling                    */

		/* precalculate skm1 = 1./xk*/
		tools_m1( N, var.sk, var.skm1);
		
		/* xkskm1 = xk.*skm1;*/
		tools_hp( N, var.xk, var.skm1, var.xkskm1);

		/* d1 = xkskm1(1:m)+xkskm1(m+1:2*m); */
		tools_copy( M, var.xkskm1, var.d1);
		tools_axpy( M, 1.0, &(var.xkskm1[sig.m]), var.d1);

		/* d2 = xkskm1(2*m+1:2*m+n)+xkskm1(2*m+n+1:end); */
		tools_copy( sig.n, &(var.xkskm1[2*sig.m]), var.d2);
		tools_axpy( sig.n, 1.0, &(var.xkskm1[2*sig.m+sig.n]),var.d2);
		
		/* Kp = Y*diag(d2)*Y'+diag(d1);*/
		tools_md( sig.m, sig.n, sig.Y, var.d2, var.T);

		tools_gemm( 'N', 'T', sig.m, sig.m, sig.n, 
								 1.0, var.T, sig.m, sig.Y, sig.m, 0.0, var.Kp, sig.m);

		tools_ad( M, var.d1, var.Kp);

		/* Cholesky factorization, Kp lower triangular part is the C, upper is trash */
    #if defined(SD) && defined(SINGLE)
		tools_potrf( 'L', sig.m, var.Kp, sig.m, &info);
		if( info != 0){
			if( set.verbose )
				printf("info = %d\n", info);

			*sol.status = FAILED;
			*sol.kp = k;
			break;
		}
		#else
		tools_copy( sig.m*sig.m, var.Kp, var.Kpp);

		tools_potrf( 'L', sig.m, var.Kp, sig.m, &info);

		if( info != 0 && set.verbose )
				printf("info = %d\n", info);

		count = 0;
		while( info != 0 && count < COUNT_MAX){
			
			tools_ads( sig.m, set.epsilon, var.Kpp);
			tools_copy( sig.m*sig.m, var.Kpp, var.Kp);
			tools_potrf( 'L', sig.m, var.Kp, sig.m, &info);
			if( set.verbose )
				printf("info = %d\n", info);

			count++;
		}
		
		if( info != 0 && count == COUNT_MAX ){
			 *sol.status = FAILED;
			 *sol.kp = k;
			 break;
		}
		#endif

		/* Compute right hand side 	
			 ra = (-rb-Af((xkskm1).*rc-xk,Y,m,n)); */
		tools_hp( N, var.xkskm1, var.rc, var.tN);
		tools_axpy( N, -1.0, var.xk, var.tN); 
		multiply_A( sig, var.tN, var.ra, var.tn); 
		tools_axpy( sig.m, 1.0, var.rb, var.ra);
		tools_scal( sig.m, -1.0, var.ra);

		/* Computes solution*/
		tools_potrs( 'L', sig.m, 1, var.Kp, sig.m, var.ra, sig.m, &info);

		/* dsa = -rc-Atf(dla,Y); */
		multiply_At( sig, var.ra, var.dsa);
		tools_scal( N, -1.0, var.dsa );
		tools_axpy( N, -1.0, var.rc, var.dsa);

		/* dxa = -xk-(xkskm1).*dsa;*/
		tools_hp( N, var.xkskm1, var.dsa, var.dxa);
		tools_scal( N, -1.0, var.dxa );
		tools_axpy( N, -1.0, var.xk, var.dxa);

		
		/* Calculate the primal-dual stepsize */
		axa = tools_stepsize( N, var.dxa, var.xk);
		asa = tools_stepsize( N, var.dsa, var.sk);


		/* and obtain the affine scalingx */
		/* mu_aff = (xk+axa*dxa)'*(sk+asa*dsa)/N; */
		tools_copy( N, var.xk, var.tN);
		tools_axpy( N, axa, var.dxa, var.tN );

		tools_copy( N, var.sk, var.tN2);
		tools_axpy( N, asa, var.dsa, var.tN2 );
		
		mu_aff = tools_dot( N, var.tN, var.tN2)/N;

		/* 	sigma = (mu_aff/mu)^3; */
		sigma = mu_aff/mu; sigma = sigma*sigma*sigma;


		/* **************************************** */
		/* 3. Calculate the search direction        */

		/* rs = (-rb-Af((xkskm1).*rc-xk-skm1.*(dxa.*dsa - sigma*mu),Y,m,n));*/

		tools_hp( N, var.dxa, var.dsa, var.tN2);
		tools_add( N, -sigma*mu, var.tN2);
		tools_hp( N, var.skm1, var.tN2, var.tN);
		tools_axpy( N, 1.0, var.xk, var.tN);

		tools_hp( N, var.xkskm1, var.rc, var.tN2);
		tools_axpy( N, -1.0, var.tN2, var.tN);
		multiply_A( sig, var.tN, var.rs, var.tN);
		tools_axpy( M, -1.0, var.rb, var.rs);

		/* Computes solution*/
		tools_potrs( 'L', sig.m, 1, var.Kp, sig.m, var.rs, sig.m, &info);

		/* dss = -rc-A'*dls; */
		multiply_At( sig, var.rs, var.dss);
		tools_scal( N, -1.0, var.dss );
		tools_axpy( N, -1.0, var.rc, var.dss);

		/* dxs = -xk-skm1.*(dxa.*dsa - sigma*mu)-(xkskm1).*dss; */
		tools_hp( N, var.dxa, var.dsa, var.tN2);
		tools_add( N, -sigma*mu, var.tN2);
		tools_hp( N, var.skm1, var.tN2, var.tN);
		tools_axpy( N, 1.0, var.xk, var.tN);
		tools_hp( N, var.xkskm1, var.dss, var.dxs);
		tools_axpy( N, 1.0, var.tN, var.dxs);
		tools_scal( N, -1.0, var.dxs);

		/* 	Calculate the primal-dual stepsize and take the step */
		etak = 1.0-0.1/((FTYPE)k*k*k);

		ap = tools_stepsize_scale( N, etak, var.dxs, var.xk);

		ad = tools_stepsize_scale( N, etak, var.dss, var.sk);

		/* 	xk = xk+ap*dxs;
				lk = lk+ad*dls;
				sk = sk+ad*dss; */
		tools_axpy( N, ap, var.dxs, var.xk);
		tools_axpy( M, ad, var.rs, var.lk); //drs = dls
		tools_axpy( N, ad, var.dss, var.sk); //drs = dls

	}

	/* Obtain the original variables */
	/* 	a = (-xk(2*m+1:2*m+n)+xk(2*m+n+1:2*m+2*n))*0.5; */
	tools_copy( sig.n, &(var.xk[2*sig.m]), sol.a);
	tools_scal( sig.n, -1.0, sol.a);
	tools_axpy( sig.n, 1.0, &(var.xk[2*sig.m+sig.n]), sol.a);
	tools_scal( sig.n, 0.5, sol.a);

	*sol.kp = k;
}

/* Multiplication with standard form A, y = A*x*/
void multiply_A( FFTYPE(isignal) sig, FTYPE *x, FTYPE *y, FTYPE *tn){
	int one = 1;

	/* A = [-I I -X X] */
	tools_copy( sig.m, &(x[sig.m]), y);
	tools_axpy( sig.m, -1.0, x, y);

	tools_copy( sig.n, &(x[2*sig.m+sig.n]), tn);
	tools_axpy( sig.n, -1.0, &(x[2*sig.m]), tn);

	/*          transform  lenY lenX alpha  a  lda  X  incX  beta  Y, incY */
	tools_gemv( 'N', sig.m, sig.n, one, sig.Y, sig.m, tn, one, 1.0, y, one);

}

/* Multiplication with standard form AT, x = A'*y */
void multiply_At( FFTYPE(isignal) sig, FTYPE *y, FTYPE *x){
	int one = 1;

	/* AT = [-I
		       I 
           -X'
           X'] */

	tools_copy( sig.m, y, &x[0]);
	tools_scal( sig.m, -1.0, &x[0]);
	tools_copy( sig.m, y, &x[sig.m]);
	
	/*         transform  lenY lenX alpha  a  lda  X  incX  beta  Y, incY */
	tools_gemv( 'T', sig.m, sig.n, one, sig.Y, sig.m, y, one, 0.0, &x[2*sig.m+sig.n], one);
	tools_copy( sig.n, &x[2*sig.m+sig.n], &x[2*sig.m]);
	tools_scal( sig.n, -1.0, &x[2*sig.m]);
}


/* Free allocated variables */
void free_variables( FFTYPE(variables) var, FFTYPE(solution) sol){
	free(sol.a);
	free(sol.kp);
	free(sol.status);

	free(var.c);
	free(var.xk);
	free(var.sk);
	free(var.lk);
	free(var.b);
	free(var.rb);
	free(var.rc);
	free(var.tn);
	free(var.skm1);
	free(var.xkskm1);
	free(var.d1);
	free(var.d2);
	free(var.T);
	free(var.Kp);
	free(var.Kpp);
	free(var.ra);
	free(var.tN);
	free(var.tN2);
	free(var.tm);
	free(var.dxa);
	free(var.dsa);
	free(var.rs);
	free(var.dss);
	free(var.dxs);

}

/* Allocates and initializes the variables */
void init_variables( FFTYPE(isignal) sig, FFTYPE(settings) set, FFTYPE(variables) *var, FFTYPE(solution) *sol){
	
	/* allocate memory for solution */
	sol->a = (FTYPE*) malloc(sig.n*sizeof(FTYPE));
	sol->kp = (int*) malloc(sizeof(FTYPE));
	sol->status = (int*) malloc(sizeof(FTYPE));

	/* allocate memory for variables */
	var->c = (FTYPE*) malloc( (2*sig.m+2*sig.n)*sizeof(FTYPE));

	var->b = (FTYPE*) malloc(sig.m*sizeof(FTYPE));

	var->xk = (FTYPE*) malloc( (2*sig.m+2*sig.n)*sizeof(FTYPE));
	
	var->sk = (FTYPE*) malloc( (2*sig.m+2*sig.n)*sizeof(FTYPE));

	var->lk = (FTYPE*) malloc( sig.m*sizeof(FTYPE));

	var->rb = (FTYPE*) malloc( sig.m*sizeof(FTYPE));

	var->rc = (FTYPE*) malloc( (2*sig.m+2*sig.n)*sizeof(FTYPE));

	var->tn = (FTYPE*) malloc( sig.n*sizeof(FTYPE));

	var->skm1 = (FTYPE*) malloc( (2*sig.m+2*sig.n)*sizeof(FTYPE));

	var->xkskm1 = (FTYPE*) malloc( (2*sig.m+2*sig.n)*sizeof(FTYPE));

	var->d1 = (FTYPE*) malloc( sig.m*sizeof(FTYPE));

	var->d2 = (FTYPE*) malloc( sig.n*sizeof(FTYPE));

	var->T = (FTYPE*) malloc( sig.n*sig.m*sizeof(FTYPE));

	var->Kp = (FTYPE*) malloc( sig.m*sig.m*sizeof(FTYPE));

	var->Kpp = (FTYPE*) malloc( sig.m*sig.m*sizeof(FTYPE));

	var->tN = (FTYPE*) malloc( (2*sig.m + 2*sig.n)*sizeof(FTYPE));

	var->tN2 = (FTYPE*) malloc( (2*sig.m + 2*sig.n)*sizeof(FTYPE));

	var->tm = (FTYPE*) malloc( sig.m*sizeof(FTYPE));

	var->ra = (FTYPE*) malloc( sig.m*sizeof(FTYPE));

	var->dsa = (FTYPE*) malloc( (2*sig.m + 2*sig.n)*sizeof(FTYPE));

	var->dxa = (FTYPE*) malloc( (2*sig.m + 2*sig.n)*sizeof(FTYPE));

	var->rs = (FTYPE*) malloc( (sig.m)*sizeof(FTYPE));

	var->dss = (FTYPE*) malloc( (2*sig.m + 2*sig.n)*sizeof(FTYPE));

	var->dxs = (FTYPE*) malloc( (2*sig.m + 2*sig.n)*sizeof(FTYPE));

}

void reset_problem( FFTYPE(isignal) sig, FFTYPE(variables) *var){
	int i;
	
	for( i = 0 ; i < 2*sig.m ; ++i)
		var->c[i] = 1.0;

	for( i = 2*sig.m ; i < 2*sig.m+2*sig.n ; ++i)
		var->c[i] = sig.gamma;

	for( i = 0 ; i < sig.m ; ++i)
		var->b[i] = 2*sig.y[i];

}

void reset_variables( FFTYPE(isignal) sig, FFTYPE(variables) *var){
	int i;
	
	for( i = 0 ; i < 2*sig.m+2*sig.n ; ++i )
		var->xk[i] = 1.0;

	for( i = 0 ; i < 2*sig.m+2*sig.n ; ++i )
		var->sk[i] = 1.0;

	for( i = 0 ; i < sig.m ; ++i)
		var->lk[i] = var->b[i];

}
