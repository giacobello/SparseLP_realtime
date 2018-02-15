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
#include "core_primal_single.h"
#include "tools_single.h"
#elif defined(SD) && defined(DOUBLE)
#include "core_primal_double.h"
#include "tools_double.h"
#else
#include "tools.h"
#include "core_primal.h"
#endif

#define SOLVED 1
#define FAILED 0

#define MAXIT 30

#define COUNT_MAX 5

#define MIN(A,B) (A < B ? A : B)
#define MAX(A,B) (A < B ? B : A)


/* Solves the sparse linear prediction problem using the dual algorithm */
void cpdip_slp_core( FFTYPE(isignal) sig, FFTYPE(settings) set,	FFTYPE(variables) var, FFTYPE(solution) sol)
{
	INT k,M,N;
	
	INTA info;
	INT count;

	FTYPE mu, pinfeas, dinfeas, relgap, cxk;
	FTYPE asa, aza, mu_aff, sigma;
	FTYPE etak, ap, ad;

	M = sig.n+sig.m;
	N = 2*M;

	*sol.status = FAILED; /* Set initial */

	if( set.verbose ){
		printf("- Sparse Linear Prediction (primal) -\n");
		printf("Ite.     pinfeas   dinfeas    rel_gap\n");
	}

	/* main iteration */
	for( k = 1 ; k <= MAXIT ; ++k ){
		/* ************************************************** */
		/* 1. compute residual and evaluate stopping criteria */

		/* 		g = A*uk;
					rh1 = sk1 + g - vk - b;
					rh2 = sk2 - g - vk + b;
		*/

		if( k == 1)
			multiply_A( sig, var.uk, var.g);
		else
			tools_axpy( M, ap, var.gs, var.g);

		tools_copy( M, var.sk1, var.rh1);
		tools_axpy( M, 1.0, var.g, var.rh1);
		tools_axpy( M, -1.0, var.vk, var.rh1);
		tools_axpy( M, -1.0, var.b, var.rh1);

		tools_copy( M, var.sk2, var.rh2);
		tools_axpy( M, -1.0, var.g, var.rh2);
		tools_axpy( M, -1.0, var.vk, var.rh2);
		tools_axpy( M, 1.0, var.b, var.rh2);


		/* mu = (sk1'*zk1+sk2'*zk2)/N;*/
		mu = (tools_dot( M, var.sk1,var.zk1)+tools_dot( M, var.sk2,var.zk2))/N;
		
		/* pinfeas = norm([rh1;rh2])/(one+sqrt(2)*norm(b)); */
	
		pinfeas = sqrt( (FTYPE)(pow(tools_nrm2( M, var.rh1),2) + pow(tools_nrm2( M, var.rh2),2)) )
			/ (1.0+1.4142*tools_nrm2( M, var.b));

		/* relgap = (sum(vk)+b'*(zk1-zk2))/(one+abs(sum(vk))); */

		cxk = tools_sum( M, var.vk);
		tools_copy( M, var.zk1, var.tM);
		tools_axpy( M, -1.0, var.zk2, var.tM);
		relgap = (cxk + tools_dot( M, var.b, var.tM))/
			(1.0+fabs(cxk));

		/* rc1 = -Atf(zk1-zk2,m,Y,gamma); */
		/* rc2 = zk1 + zk2 - 1; */
		multiply_At(sig, var.tM, var.rc1);
		tools_scal( sig.n, -1.0, var.rc1);

		tools_copy( M, var.zk1,var.rc2);
		tools_axpy( M, 1.0, var.zk2, var.rc2);
		tools_add( M, -1.0, var.rc2 );

		/* dinfeas = norm([rc1;rc2])/(one+sqrt(M)); */

		dinfeas = sqrt( (FTYPE) ( pow(tools_nrm2( sig.n, var.rc1),2) + pow(tools_nrm2( M, var.rc2),2) ) )
			/ (1.0+sqrt((FTYPE)M));

		if( set.verbose ){
			printf("%2d     %.2e    %.2e    %.2e\n", k, pinfeas, dinfeas, relgap);
		}

		/* Stop if sufficiently accurate */
		if( relgap < set.epsilon && pinfeas < set.epsilon){
			*sol.status = SOLVED;
			break;
		}

		/* ************************************************** */
		/* 2. calculate the affine scaling                    */

		/* precalculate skm1 = 1./sk*/
		tools_m1( M, var.sk1, var.sk1m1);
		tools_m1( M, var.sk2, var.sk2m1);

		tools_m1( M, var.zk1, var.zk1m1);
		tools_m1( M, var.zk2, var.zk2m1);

		/* 	rs1 = sk1.*zk1;
				rs2 = sk2.*zk2; */

		tools_hp( M, var.sk1, var.zk1, var.rs1);
		tools_hp( M, var.sk2, var.zk2, var.rs2);

		/*d1 = zk1./sk1;
			d2 = zk2./sk2;*/

		tools_hp( M, var.zk1, var.sk1m1, var.d1);
		tools_hp( M, var.zk2, var.sk2m1, var.d2);

		/*p1 = (zk1.*rh1-rs1)./sk1;
			p2 = (zk2.*rh2-rs2)./sk2; */

		tools_hp( M, var.zk1, var.rh1, var.tM);
		tools_axpy( M, -1.0, var.rs1, var.tM);
		tools_hp( M, var.tM, var.sk1m1, var.p1);

		tools_hp( M, var.zk2, var.rh2, var.tM);
		tools_axpy( M, -1.0, var.rs2, var.tM);
		tools_hp( M, var.tM, var.sk2m1, var.p2);

		/* r1 = - Atf(p1-p2,Y,gamma,m);
			 r2 = rc2+p1+p2;*/
		/*tools_copy( M, var.p2, var.tM);
		tools_axpy( M, -1.0, var.p1, var.tM);
		multiply_At( sig, var.tM, var.r1);
		tools_axpy( sig.n, 1.0, var.rc1, var.r1);*/

		tools_copy( M, var.p1, var.r2);
		tools_axpy( M, 1.0, var.p2, var.r2);
		tools_axpy( M, 1.0, var.rc2, var.r2);

		/* 	d = 4*(d1.*d2./(d1+d2));*/
		tools_copy( M, var.d1, var.tM);
		tools_axpy( M, 1.0, var.d2, var.tM);
		tools_m1( M, var.tM, var.d1pd2m1);
		tools_hp( M, var.d1pd2m1, var.d1, var.tM);
		tools_hp( M, var.tM, var.d2, var.d);
		tools_scal( M, 4.0, var.d);

		/* Kp = Y*diag(dt)*Y'+gamma^2diag(dp);*/
		tools_mdt( sig.m, sig.n, sig.Y, var.d, var.T);

		tools_gemm( 'T', 'N', sig.n, sig.n, sig.m, 
								 1.0, sig.Y, sig.m, var.T, sig.m, 0.0, var.Kp, sig.n);

		tools_scal( sig.n, pow(sig.gamma,2), &(var.d[sig.m]) );
		tools_ad( sig.n, &(var.d[sig.m]), var.Kp);

		/* Cholesky factorization, Kp lower triangular part is the C, upper is trash */
    #if defined(SD) && defined(SINGLE)
		tools_potrf( 'L', sig.n, var.Kp, sig.n, &info);
		if( info != 0){
			if( set.verbose )
				printf("potrf: info = %d\n", info);

			*sol.status = FAILED;
			*sol.kp = k;
			break;
		}
		#else
		/* Copy matrix if refactorization is needed */
		tools_copy( sig.n*sig.n, var.Kp, var.Kpp);
		
		tools_potrf( 'L', sig.n, var.Kp, sig.n, &info);



		if( info != 0 && set.verbose )
				printf("potrf: info = %d\n", info);

		count = 0;
		while( info != 0 && count < COUNT_MAX){
			if( set.verbose )
				printf("potrf: info = %d\n", info);
			
			tools_ads( sig.n, 10.0, var.Kpp);
			tools_copy( sig.n*sig.n, var.Kpp, var.Kp);
			tools_potrf( 'L', sig.n, var.Kp, sig.n, &info);
			count++;
		}
		
		if( info != 0 && count == COUNT_MAX ){
			 *sol.status = FAILED;
			 *sol.kp = k;
			 break;
		}
		#endif

		
		/* Compute right hand side 	
			 rhs = r1-Atf( ((d2-d1)./(d1+d2)).*r2,Y,gamma,m);
		*/

		tools_copy( M, var.d2, var.tM);
		tools_axpy( M, -1.0, var.d1, var.tM);
		tools_hp( M, var.tM, var.d1pd2m1, var.tM2);
		tools_hp( M, var.tM2, var.r2, var.tM);
		tools_axpy( M, -1.0, var.p2, var.tM);
		tools_axpy( M, 1.0, var.p1, var.tM);
		multiply_At( sig, var.tM, var.rhs);

		tools_scal( sig.n, -1.0, var.rhs);
		tools_axpy( sig.n, 1.0, var.rc1, var.rhs);


		/* Computes solution*/
		tools_potrs( 'L', sig.n, 1, var.Kp, sig.n, var.rhs, sig.n, &info);
		


		/* 	dv = -(d2-d1).*gs./(d1+d2)+r2./(d1+d2); */
		multiply_A( sig, var.rhs, var.gs);
		tools_copy( M, var.d1, var.tM);
		tools_axpy( M, -1.0, var.d2, var.tM);
		tools_hp( M, var.tM, var.gs, var.tM2);
		tools_axpy( M, 1.0, var.r2, var.tM2);
		tools_hp( M, var.tM2, var.d1pd2m1, var.dv);

		/* dza1 = -(zk1./sk1).*(-rh1+rs1./zk1-g+dv);
			 dza2 = -(zk2./sk2).*(-rh2+rs2./zk2+g+dv);*/
		tools_hp( M, var.rs1, var.zk1m1, var.tM);
		tools_axpy( M, -1.0, var.gs, var.tM); 
		tools_axpy( M, 1.0, var.dv, var.tM);
		tools_axpy( M, -1.0, var.rh1, var.tM);
		tools_hp( M, var.d1, var.tM, var.dza1);
		tools_scal( M, -1.0, var.dza1);
		
		tools_hp( M, var.rs2, var.zk2m1, var.tM);
		tools_axpy( M, 1.0, var.gs, var.tM); 
		tools_axpy( M, 1.0, var.dv, var.tM);
		tools_axpy( M, -1.0, var.rh2, var.tM);
		tools_hp( M, var.d2, var.tM, var.dza2);
		tools_scal( M, -1.0, var.dza2);

		/* dsa1 = (-rs1-sk1.*dza1)./zk1;
			 dsa2 = (-rs2-sk2.*dza2)./zk2; */

		tools_hp( M, var.sk1, var.dza1, var.tM);
		tools_axpy( M, 1.0, var.rs1, var.tM );
		tools_hp( M, var.tM, var.zk1m1, var.dsa1 );
		tools_scal( M, -1.0, var.dsa1);

		tools_hp( M, var.sk2, var.dza2, var.tM);
		tools_axpy( M, 1.0, var.rs2, var.tM );
		tools_hp( M, var.tM, var.zk2m1, var.dsa2 );
		tools_scal( M, -1.0, var.dsa2);

		/* Calculate the primal-dual stepsize */
		asa = tools_stepsize( M, var.dsa1, var.sk1);
		asa = MIN(tools_stepsize( M, var.dsa2, var.sk2),asa);

		aza = tools_stepsize( M, var.dza1, var.zk1);
		aza = MIN(tools_stepsize( M, var.dza2, var.zk2),asa);

		/* and obtain the affine scalingx */
		/* 	mu_aff = ((sk1+asa*dsa1)'*(zk1+aza*dza1)+(sk2+asa*dsa2)'*(zk2+aza*dza2))/N;*/
		tools_copy( M, var.sk1, var.tM);
		tools_axpy( M, asa, var.dsa1, var.tM );
		
		tools_copy( M, var.zk1, var.tM2);
		tools_axpy( M, aza, var.dza1, var.tM2 );
		
		mu_aff = tools_dot( M, var.tM, var.tM2);

		tools_copy( M, var.sk2, var.tM);
		tools_axpy( M, asa, var.dsa2, var.tM );
		
		tools_copy( M, var.zk2, var.tM2);
		tools_axpy( M, aza, var.dza2, var.tM2 );
		
		mu_aff = (mu_aff+tools_dot( M, var.tM, var.tM2))/N;

		/* 	sigma = (mu_aff/mu)^3; */
		sigma = pow(mu_aff/mu,3);

		/* **************************************** */
		/* 3. Calculate the search direction        */

		/* 	rs1 = sk1.*zk1+dsa1.*dza1-sigma*mu;
				rs2 = sk2.*zk2+dsa2.*dza2-sigma*mu;
		*/
		tools_hp( M, var.sk1, var.zk1, var.rs1);
		tools_hp( M, var.dsa1, var.dza1, var.tM);
		tools_axpy( M, 1.0, var.tM, var.rs1);
		tools_add( M, -sigma*mu, var.rs1);

		tools_hp( M, var.sk2, var.zk2, var.rs2);
		tools_hp( M, var.dsa2, var.dza2, var.tM);
		tools_axpy( M, 1.0, var.tM, var.rs2);
		tools_add( M, -sigma*mu, var.rs2);
		
		/*p1 = (zk1.*rh1-rs1)./sk1;
			p2 = (zk2.*rh2-rs2)./sk2; */

		tools_hp( M, var.zk1, var.rh1, var.tM);
		tools_axpy( M, -1.0, var.rs1, var.tM);
		tools_hp( M, var.tM, var.sk1m1, var.p1);

		tools_hp( M, var.zk2, var.rh2, var.tM);
		tools_axpy( M, -1.0, var.rs2, var.tM);
		tools_hp( M, var.tM, var.sk2m1, var.p2);

		/* r1 = - Atf(p1-p2,Y,gamma,m);
			 r2 = p1+p2;*/
		/*tools_copy( M, var.p2, var.tM);
		tools_axpy( M, -1.0, var.p1, var.tM);
		multiply_At( sig, var.tM, var.r1);
		tools_axpy( sig.n, 1.0, var.rc1, var.r1);*/

		tools_copy( M, var.p1, var.r2);
		tools_axpy( M, 1.0, var.p2, var.r2);
		tools_axpy( M, 1.0, var.rc2, var.r2);
		/* 
			 rhs = r1-Atf( ((d2-d1)./(d1+d2)).*r2,Y,gamma,m);
		*/
		tools_copy( M, var.d2, var.tM);
		tools_axpy( M, -1.0, var.d1, var.tM);
		tools_hp( M, var.tM, var.d1pd2m1, var.tM2);
		tools_hp( M, var.tM2, var.r2, var.tM);
		tools_axpy( M, -1.0, var.p2, var.tM);
		tools_axpy( M, 1.0, var.p1, var.tM);
		multiply_At( sig, var.tM, var.rhs);
		tools_scal( sig.n, -1.0, var.rhs);
		tools_axpy( sig.n, 1.0, var.rc1, var.rhs);


		/* Computes solution*/
		tools_potrs( 'L', sig.n, 1, var.Kp, sig.n, var.rhs, sig.n, &info);

		/* 	dvs = -(d2-d1).*gs./(d1+d2)+r2./(d1+d2); */
		multiply_A( sig, var.rhs, var.gs);
		tools_copy( M, var.d1, var.tM);
		tools_axpy( M, -1.0, var.d2, var.tM);
		tools_hp( M, var.tM, var.gs, var.tM2);
		tools_axpy( M, 1.0, var.r2, var.tM2);
		tools_hp( M, var.tM2, var.d1pd2m1, var.dv);

		/* dzs1 = -(zk1./sk1).*(-rh1+rs1./zk1-g+dvs);
			 dzs2 = -(zk2./sk2).*(-rh2+rs2./zk2+g+dvs);*/
		tools_hp( M, var.rs1, var.zk1m1, var.tM);
		tools_axpy( M, -1.0, var.gs, var.tM); 
		tools_axpy( M, 1.0, var.dv, var.tM);
		tools_axpy( M, -1.0, var.rh1, var.tM);
		tools_hp( M, var.d1, var.tM, var.dzs1);
		tools_scal( M, -1.0, var.dzs1);

		tools_hp( M, var.rs2, var.zk2m1, var.tM);
		tools_axpy( M, 1.0, var.gs, var.tM); 
		tools_axpy( M, 1.0, var.dv, var.tM);
		tools_axpy( M, -1.0, var.rh2, var.tM);
		tools_hp( M, var.d2, var.tM, var.dzs2);
		tools_scal( M, -1.0, var.dzs2);

		/* 	dss1 = (-rs1-sk1.*dzs1)./zk1;
				dss2 = (-rs2-sk2.*dzs2)./zk2; */

		tools_hp( M, var.sk1, var.dzs1, var.tM);
		tools_axpy( M, 1.0, var.rs1, var.tM );
		tools_hp( M, var.tM, var.zk1m1, var.dss1 );
		tools_scal( M, -1.0, var.dss1);

		tools_hp( M, var.sk2, var.dzs2, var.tM);
		tools_axpy( M, 1.0, var.rs2, var.tM );
		tools_hp( M, var.tM, var.zk2m1, var.dss2 );
		tools_scal( M, -1.0, var.dss2);


		/* 	Calculate the primal-dual stepsize and take the step */
		etak = 1.0-0.1/((FTYPE)k*k*k);

		ap = tools_stepsize_scale( M, etak, var.dss1, var.sk1);
		ap = MIN( ap, tools_stepsize_scale( M, etak, var.dss2, var.sk2));

		ad = tools_stepsize_scale( M, etak, var.dzs1, var.zk1);
		ad = MIN( ad, tools_stepsize_scale( M, etak, var.dzs2, var.zk2));

		/* 		xk = xk+ap*dxs;
					sk = sk+ap*dss;
					zk = zk+ad*dzs;
		*/
		tools_axpy( sig.n, ap, var.rhs, var.uk);
		tools_axpy( M, ap, var.dv, var.vk);
		
		tools_axpy( M, ap, var.dss1, var.sk1);
		tools_axpy( M, ap, var.dss2, var.sk2);
		
		tools_axpy( M, ad, var.dzs1, var.zk1);
		tools_axpy( M, ad, var.dzs2, var.zk2);

	}

	/* Obtain the original variables */
	tools_copy( sig.n, var.uk, sol.a);
	*sol.kp = k;
}

/* Multiplication with standard form A, y = A*x*/
void multiply_A( FFTYPE(isignal) sig, FTYPE *x, FTYPE *y){
	INT one = 1;

	/* A = [-X; 
		 gamma I] */

	/*      transform lenY   lenX alpha  A     lda    X  incX beta Y, incY */
	tools_gemv( 'N', sig.m, sig.n, 1.0, sig.Y, sig.m, x, one, 0.0, y, one);

	/* tools_gemv( 'N', sig.m, sig.n, sig.Y, x, y); */
	tools_scal( sig.m, -1.0, y);
	tools_copy( sig.n, x, &y[sig.m]);
	tools_scal( sig.n, sig.gamma, &y[sig.m]);

}

/* Multiplication with standard form AT, x = A'*y */
void multiply_At( FFTYPE(isignal) sig, FTYPE *y, FTYPE *x){
	INT one = 1;

	/* AT = [-X' gamma*I] */

	tools_copy( sig.n, &y[sig.m], x);

	/*      transform lenY   lenX alpha  A     lda    X  incX   beta      Y, incY */
	tools_gemv( 'T', sig.m, sig.n, 1.0, sig.Y, sig.m, y, one, -sig.gamma, x, one);

	tools_scal( sig.n, -1.0, x);

	/* tools_gemv( 'T', sig.m, sig.n, sig.Y, y, x); */
	/* tools_scal( sig.n, -1.0, x); */
	/* tools_axpy( sig.n, sig.gamma, &y[sig.m], x); */
}



/* Free allocated variables */
void free_variables( FFTYPE(variables) var, FFTYPE(solution) sol){
	free(sol.a);
	free(sol.kp);
	free(sol.status);

	free(var.uk);
	free(var.vk);
	free(var.dv);
	free(var.g);
	free(var.gs);
	free(var.b);
	free(var.sk1);
	free(var.sk2);
	free(var.zk1);
	free(var.zk2);
	free(var.rh1);
	free(var.rh2);
	free(var.p1);
	free(var.p2);
	free(var.d1);
	free(var.d2);
	free(var.sk1m1);
	free(var.sk2m1);
	free(var.rs1);
	free(var.rs2);
	free(var.r1);
	free(var.r2);
	free(var.rc1);
	free(var.rc2);
	free(var.rhs);
	free(var.tM);
	free(var.tM2);
	free(var.d);
	free(var.d1pd2m1);
	free(var.zk1m1);
	free(var.zk2m1);
	free(var.dza1);
	free(var.dza2);
	free(var.dsa1);
	free(var.dsa2);
	free(var.dzs1);
	free(var.dzs2);
	free(var.dss1);
	free(var.dss2);	
	free(var.Kp);
	free(var.Kpp);
	free(var.T);

	/* free variables used for initialization and SVD */
	free(var.work);
	free(var.iwork);
	free(var.s);
	free(var.U);
	free(var.Vt);
	free(var.t1);
	free(var.t2);
	free(var.t3);
}

/* Allocates and initializes the variables */
void init_variables( FFTYPE(isignal) sig, FFTYPE(settings) set, FFTYPE(variables) *var, FFTYPE(solution) *sol){
	
	INT M;

	M = sig.n+sig.m;

	/* allocate memory for solution */
	sol->a = (FTYPE*) malloc(sig.n*sizeof(FTYPE));
	sol->kp = (INT*) malloc(sizeof(FTYPE));
	sol->status = (INT*) malloc(sizeof(FTYPE));

	/* allocate memory for variables */
	var->b = (FTYPE*) malloc(M*sizeof(FTYPE));

	var->uk = (FTYPE*) malloc( sig.n*sizeof(FTYPE));
	var->vk = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->g = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->gs = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->sk1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->sk2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->zk1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->zk2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->rh1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->rh2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->p1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->p2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->d1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->d2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->sk1m1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->sk2m1 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->p1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->p2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->rs1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->rs2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->rhs = (FTYPE*) malloc( sig.n*sizeof(FTYPE));
	var->r1 = (FTYPE*) malloc( sig.n*sizeof(FTYPE));
	var->r2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->rc1 = (FTYPE*) malloc( sig.n*sizeof(FTYPE));
	var->rc2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->tM = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->tM2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->d = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->d1pd2m1 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->dv = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->zk1m1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->zk2m1 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->dza1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->dza2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->dsa1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->dsa2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->dzs1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->dzs2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->dss1 = (FTYPE*) malloc( M*sizeof(FTYPE));
	var->dss2 = (FTYPE*) malloc( M*sizeof(FTYPE));

	var->Kp = (FTYPE*) malloc( sig.n*sig.n*sizeof(FTYPE));

	var->Kpp = (FTYPE*) malloc( sig.n*sig.n*sizeof(FTYPE));

	var->T = (FTYPE*) malloc( sig.n*sig.m*sizeof(FTYPE));

	/* For initialization */
	var->L = 3*MIN(sig.m,sig.n) + MAX( MAX(sig.m,sig.n), 4*MIN(sig.m,sig.n)*MIN(sig.m,sig.n)+3*MIN(sig.m,sig.n)+MAX(sig.m,sig.n) );

	var->work = (FTYPE*) malloc( var->L*sizeof(FTYPE));
  var->iwork = (INTA*) malloc( 8*MIN(sig.m,sig.n)*sizeof(INT));

	var->s = (FTYPE*) malloc( sig.n*sizeof(FTYPE));

	var->U = (FTYPE*) malloc( sig.m*sig.n*sizeof(FTYPE));
	var->Vt = (FTYPE*) malloc( sig.n*sig.n*sizeof(FTYPE));

	var->t1 = (FTYPE*) malloc( sig.n*sizeof(FTYPE));
	var->t2 = (FTYPE*) malloc( sig.n*sizeof(FTYPE));
	var->t3 = (FTYPE*) malloc( sig.n*sizeof(FTYPE));

}

/* (re)set problem parameters to the initial values. */
void reset_problem( FFTYPE(isignal) sig, FFTYPE(variables) *var){
	INT i;
	INT M;

	M = sig.n + sig.m;

	/* set b vector */
	for( i = 0 ; i < sig.m ; ++i)
		var->b[i] = -sig.y[i];
	for( i = sig.m ; i < M ; ++i)
		var->b[i] = 0.0;
}

/* (re)set variables to the initial values. */
void reset_variables( FFTYPE(isignal) sig, FFTYPE(variables) *var){
	INTA info;
	INT M, i, one = 1;;

	M = sig.n + sig.m;

	/*
	tools_gemm ( 'T', 'N', sig.n, sig.n, sig.m, 
							1.0, sig.Y, sig.m, sig.Y, sig.m, 0.0, var->Kp, sig.n);
	
	tools_ads( sig.n, sig.gamma*sig.gamma, var->Kp);

	tools_potrf( 'L', sig.n, var->Kp, sig.n, &info);
		
	prINTf("info potrf %d\n",info);
	
	tools_gemv ( 'T', sig.m, sig.n, one, sig.Y, sig.m, sig.y, one, 0, var->uk, one);
	
	tools_potrs ( 'L', sig.n, 1, var->Kp, sig.n, var->uk, sig.n, &info);
	printf("info potrs %d\n",info);
	*/

	/* calculate Tikhonov regularisation via the SVD */
	/* u0 = A\b; */		
	/* Copy Y to a T matrix to preserve Y*/
	
	for( i = 0 ; i < (sig.m)*(sig.n) ; ++i) var->T[i] = sig.Y[i];
			
	tools_gesdd( 'S', sig.m, sig.n, var->T, sig.m, var->s, var->U, sig.m, 
							 var->Vt, sig.n , var->work, var->L, var->iwork, &info);
		
	if( info != 0 ){
		printf(" gesdd failed in reset_variables with code %d\n",info);
		return;
	}
	
	tools_gemv( 'T', sig.m, sig.n, 1.0, var->U, sig.m, sig.y, one, 0.0, var->t1, one);
	
	tools_hp( sig.n, var->s, var->s, var->t2);
	tools_add( sig.n, sig.gamma*sig.gamma, var->t2);
	tools_m1( sig.n, var->t2, var->t3);
	tools_hp( sig.n, var->t3, var->s, var->t2);
	tools_hp( sig.n, var->t1, var->t2, var->t3);
	
	tools_gemv( 'T', sig.n, sig.n, 1.0, var->Vt, sig.n, var->t3, one, 0.0, var->uk, one);

		
	/* rls = Af(x0,Y,gamma) - b; */
	multiply_A( sig, var->uk, var->tM);
	tools_axpy( M, -1.0, var->b, var->tM);

	/* v0 = abs(rls)+1e-1; */
	tools_ea( M, var->tM, var->vk);
	tools_add( M, 0.1, var->vk);

	/* z = (1-1e-1)*rls/norm(rls,'inf');*/
	tools_scal( M, 0.9/tools_inf(M, var->tM),var->tM);

	/* 
		 zk1 = (1+z)/2;
		 zk2 = (1-z)/2;
		 sk1 = ones(length(zk1),1)*1e-2;
		 sk2 = ones(length(zk2),1)*1e-2;
	*/

	tools_copy( M, var->tM, var->zk1);
	tools_add( M, 1.0 , var->zk1);
	tools_scal( M, 0.5 , var->zk1);

	tools_copy( M, var->tM, var->zk2);
	tools_add( M, -1.0 , var->zk2);
	tools_scal( M, -0.5 , var->zk2);

	tools_set( M, 0.01, var->sk1);
	tools_set( M, 0.01, var->sk2);

}
