/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2012 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file examples/exampleLP.cpp
 *	\author Hans Joachim Ferreau
 *	\version 3.0beta
 *	\date 2008-2009
 *
 *	Very simple example for solving a LP sequence using qpOASES.
 */

/* Modified to be used for benchmarking */
/* for the sparse linear prediction */
/* Tobias Lindstrøm Jensen, Aalborg University, 2012, tlj@es.aau.dk */

#include <qpOASES.hpp>
#include <sys/time.h>
#include <stdlib.h>

USING_NAMESPACE_QPOASES

int load_problem(char *s, real_t *A, real_t *g, real_t *lb, real_t *ub, real_t *lbA, real_t *ubA);
int write_solution(char *fn, int n, real_t *alpha, int num_it, double time);

void print_matrix( int size_m, int size_n, real_t *X, char *s){
	int i,j;

	printf("%s\n",s);
	
	for( i = 0 ; i < size_m ; ++i){
		for( j = 0 ; j < size_n ; ++j)
			printf("[%d,%d] %.4f  ",i,j,X[i+j*size_m]);
		printf("\n");
	}
}

void print_vector( int size, real_t *a, char *s){
	int i;

	printf("%s\n",s);
	
	for( i = 0 ; i < size ; ++i)
		printf("%d %.4f\n",i,a[i]);
}

/** Example for qpOASES main function solving LPs. */
int main(int argc, char **argv)
{

	if( argc != 6 ){
		printf("Number of input should be 5: problem_name solution_name no_of_frames m n\n");
		return 0;
	}

	char s [100];
	sprintf(s,"../data/%s",argv[1]);
	char v [100];
	sprintf(v,"../data/%s",argv[2]);
	
	
	char f [100];
	int no_of_frames;
	sscanf(argv[3],"%d",&no_of_frames);

	int m,n;
	sscanf(argv[4],"%d",&m);
	sscanf(argv[5],"%d",&n);

	int rpt = 100;

	int i=1,j;
	int r;
	int num_it;

	int nWSR = 100000;

	double l1, l2;
	struct timeval tim;

	/* Setting up QProblem object with zero Hessian matrix. */
	#ifdef SR
	QProblem slp( n+m, m+m, HST_ZERO);
#else
	QProblem slp( n+m+n, m+m+n+n, HST_ZERO);
#endif

	Options options;
	// 	options.setToMPC();
	slp.setOptions( options );
	slp.setPrintLevel( PL_NONE );

	#ifdef SR
	real_t* alpha = (real_t*) malloc(sizeof(real_t)*(n+m+n)); 
	real_t* A = (real_t*) malloc(sizeof(real_t)*(n+m+n)*(m+m+n+n)); 
	real_t* g = (real_t*) malloc(sizeof(real_t)*(n+m+n));
	real_t* lb = (real_t*) malloc(sizeof(real_t)*(n+m+n));
	real_t* ub = (real_t*) malloc(sizeof(real_t)*(n+m+n));
	real_t* lbA = (real_t*) malloc(sizeof(real_t)*(m+m+n+n));
	real_t* ubA = (real_t*) malloc(sizeof(real_t)*(m+m+n+n));
	#else
	real_t* alpha = (real_t*) malloc(sizeof(real_t)*(n+m)); 
	real_t* A = (real_t*) malloc(sizeof(real_t)*(n+m)*(m+m)); 
	real_t* g = (real_t*) malloc(sizeof(real_t)*(n+m));
	real_t* lb = (real_t*) malloc(sizeof(real_t)*(n+m));
	real_t* ub = (real_t*) malloc(sizeof(real_t)*(n+m));
	real_t* lbA = (real_t*) malloc(sizeof(real_t)*(m+m));
	real_t* ubA = (real_t*) malloc(sizeof(real_t)*(m+m));
	#endif

	/* run over all problems in the signal */
	for( i = 1 ; i <= no_of_frames ; i++ ){

		sprintf(f,"%s_%d",s,i);
		//sprintf(f,"%s",s,i);
		if( load_problem(f, A, g, lb, ub, lbA, ubA) ){
		
			gettimeofday(&tim, NULL);
			l1 = tim.tv_sec + (tim.tv_usec/1000000.0);			
			for( r = 0 ; r < rpt; ++r){
				slp.init( 0, g, A, lb, ub, lbA, ubA, nWSR);
				slp.getPrimalSolution( alpha );
			}

			gettimeofday(&tim, NULL);
			l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
			
			sprintf(f,"%s_%d",v,i);
			write_solution( f, n, alpha, 0, (l2-l1)/rpt);
		}
		else 
			return 0;
	}
		return 1;
}

int load_problem(char *s, real_t *A, real_t *g, real_t *lb, real_t *ub, real_t *lbA, real_t *ubA) {
	FILE *f;
	int i, j, m, n;
	double gamma;
	double r;
	real_t l = 1e6;

	// open the file
	f = fopen(s,"r");
	if( f != NULL ){
		// first read m, n and gamma
		if(fscanf(f, "%d", &(m)) != 1)
			return 0;
		if(fscanf(f, "%d", &(n)) != 1)
			return 0;
		if(fscanf(f, "%lf", &gamma) != 1)
			 return 0;

		double X[m*n];
		double x[m];

		// read in X and x
		for( i = 0 ; i < m*n ; ++i ){
			if(fscanf(f,"%lf", &r) != 1) return 0; else X[i]=(double)r;
			 }
		for( i = 0 ; i < m ; ++i ){
			if(fscanf(f,"%lf", &r) != 1) return 0; else x[i]=(double)r;
		}
		// close the file
		fclose(f);

		#ifdef SR
		// write data to the form qpOASES requires
		//real_t g[4+6] = { 0,0,0,0, 1,1,1,1,1,1};
		for( i = 0 ; i < n ; ++i )
			g[i] = 0.0;
		for( i = n ; i < n+m ; ++i )
			g[i] = 1.0;

		//real_t lb[4+6] = {-l, -l, -l, -l, -l, -l, -l, -l, -l, -l };
		for( i = 0 ; i < n+m ; ++i )
			lb[i] = -l;
		//real_t ub[4+6] = { l,  l,  l,  l,  l,  l,  l,  l,  l,  l };
		for( i = 0 ; i < n+m+n ; ++i )
			ub[i] = l;
		//real_t lbA[6+6] = {-l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l };
		for( i = 0 ; i < m+m ; ++i )
			lbA[i] = -l;
		///*real_t ubA[6+6] = { 0.0050, 0.0136, 0.0297, 0.0533, 0.0800, 0.1005, -0.0050, -0.0136, -0.0297, -0.0533, -0.0800, -0.1005};*/
		for( i = 0 ; i < m ; ++i )
			ubA[i] = -x[i];
		for( i = m ; i < m+m ; ++i )
			ubA[i] = x[i-m];
		
		//	A = 
		//    -X -I 
		//		X -I
		//	in row major order

		for( j = 0 ; j < m ; ++j )
			for( i = 0 ; i < n ; ++i )
				A[j*(n+m)+i] = -X[j+i*m];

		for( j = 0 ; j < m ; ++j )
			for( i = 0 ; i < n ; ++i )
				A[(j+m)*(n+m)+i] = X[j+i*m];

		for( j = 0 ; j < m ; ++j )
			A[j*(n+m)+(j+n)] = -1.0;

		for( j = 0 ; j < m ; ++j )
			A[(j+m)*(n+m)+(j+n)] = -1.0;

		#else
		// write data to the form qpOASES requires
		//real_t g[4+6+4] = { 0,0,0,0, 1,1,1,1,1,1, gamma,gamma,gamma,gamma};
		for( i = 0 ; i < n ; ++i )
			g[i] = 0.0;
		for( i = n ; i < n+m ; ++i )
			g[i] = 1.0;
		for( i = n+m ; i < n+m+n ; ++i )
			g[i] = gamma;

		//real_t lb[4+6+4] = {-l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l};
		for( i = 0 ; i < n+m+n ; ++i )
			lb[i] = -l;
		//real_t ub[4+6+4] = { l,  l,  l,  l,  l,  l,  l,  l,  l,  l,  l,  l,  l,  l};
		for( i = 0 ; i < n+m+n ; ++i )
			ub[i] = l;
		//real_t lbA[6+6+4+4] = {-l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l, -l};
		for( i = 0 ; i < m+m+n+n ; ++i )
			lbA[i] = -l;
		///*real_t ubA[6+6+4+4] = { 0.0050, 0.0136, 0.0297, 0.0533, 0.0800, 0.1005, -0.0050, -0.0136, -0.0297, -0.0533, -0.0800, -0.1005, 0, 0, 0, 0, 0, 0, 0, 0};*/
		for( i = 0 ; i < m ; ++i )
			ubA[i] = -x[i];
		for( i = m ; i < m+m ; ++i )
			ubA[i] = x[i-m];
		for( i = m+m ; i < m+m+n+n ; ++i ) 
			ubA[i] = 0.0;

		
		//	A = 
		//    -X -I 0 
		//		X -I 0
		//		I  0 -I 
		//		-I 0 -I
		//	in row major order

		for( j = 0 ; j < m ; ++j )
			for( i = 0 ; i < n ; ++i )
				A[j*(n+m+n)+i] = -X[j+i*m];

		for( j = 0 ; j < m ; ++j )
			for( i = 0 ; i < n ; ++i )
				A[(j+m)*(n+m+n)+i] = X[j+i*m];

		for( j = 0 ; j < n ; ++j )
			A[(j+m+m)*(n+m+n)+j] = 1.0;

		for( j = 0 ; j < n ; ++j )
			A[(j+m+m+n)*(n+m+n)+j] = -1.0;

		for( j = 0 ; j < m ; ++j )
			A[j*(n+m+n)+(j+n)] = -1.0;

		for( j = 0 ; j < m ; ++j )
			A[(j+m)*(n+m+n)+(j+n)] = -1.0;

		for( j = 0 ; j < n ; ++j )
			A[(j+m+m)*(n+m+n)+(j+n+m)] = -1.0;

		for( j = 0 ; j < n ; ++j )
			A[(j+m+m+n)*(n+m+n)+(j+n+m)] = -1.0;

		#endif
		return 1;
	}
	else 
		printf("Could not load file %s\n",s);
		return 0;
}


/*writes the solution and timing to the file f*/
int write_solution(char *fn, int n, real_t *alpha, int num_it, double time){
	FILE *f;
	int i;

	f = fopen(fn,"w");
	if( f != NULL){
		fprintf( f, "%d\n", num_it);
		fprintf( f, "%d\n", 1);
		fprintf( f, "%lf\n", (double) time);
		for( i = 0 ; i < n ; i++)
			fprintf( f, "%lf\n", (double) alpha[i]);
		
		fclose(f);
		return 1;
	}
	else{
		printf("Could not open file %s in write_solution\n",fn);
		return 0;
	}
}

/*
 *	end of file
 */
