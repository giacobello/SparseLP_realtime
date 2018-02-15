/* Produced by CVXGEN, 2012-11-08 05:03:14 -0800.  */
/* CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com. */
/* The code in this file is Copyright (C) 2006-2011 Jacob Mattingley. */
/* CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial */
/* applications without prior written permission from Jacob Mattingley. */

/* Modification of the testsolver.c to be used for benchmarking */
/* for the sparse linear prediction */
/* Tobias Lindstrøm Jensen, Aalborg University, 2012, tlj@es.aau.dk */

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>


#include "solver.h"

Vars vars;
Params params;
Workspace work;
Settings settings;

int main(int argc, char **argv) {

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

	int i=1;
	int r;
	int num_it;

	double l1, l2;
	struct timeval tim;

  set_defaults();
  setup_indexing();

  settings.verbose = 0;
	
	/* run over all problems in the signal */
	for( i = 1 ; i <= no_of_frames ; i++ ){

		//printf("Run problem instance %d\n",i);

		/* get signal */
		sprintf(f,"%s_%d",s,i);
		if( load_problem(f) ){

			gettimeofday(&tim, NULL);
			l1 = tim.tv_sec + (tim.tv_usec/1000000.0);			
			for( r = 0 ; r < rpt; ++r){
				num_it = solve();
			}
			gettimeofday(&tim, NULL);
			l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
				
			sprintf(f,"%s_%d",v,i);
			write_solution( f, n, vars.alpha, num_it, (l2-l1)/rpt);
		}
		else
			return 0;
	}
	return 0;
}

int load_problem(char *s) {
	FILE *f;
	int i,m,n;
	double r;

	// open the file
	f = fopen(s,"r");
	if( f != NULL ){
		// first read m, n and gamma
		if(fscanf(f, "%d", &(m)) != 1)
			return 0;
		if(fscanf(f, "%d", &(n)) != 1)
			return 0;
		if(fscanf(f, "%lf", &r) != 1)
			 return 0;

		params.gamma[0] = (double)r;

	
		// read in X and x
		for( i = 0 ; i < m*n ; ++i ){
			if(fscanf(f,"%lf", &r) != 1) return 0; else params.X[i]=(double)r;
			 }
		for( i = 0 ; i < m ; ++i ){
			if(fscanf(f,"%lf", &r) != 1) return 0; else params.x[i]=(double)r;

		}
		// close the file
		fclose(f);
		
		return 1;
	}
	else 
		printf("Could not load file %s\n",s);
		return 0;
}

/*writes the solution and timing to the file f*/
int write_solution(char *fn, int n, double *alpha, int num_it, double time){
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
