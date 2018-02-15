#include <stdlib.h>
#include <stdio.h>

#include "ftype.h"

#include "fcore.h"

#include <sys/time.h>
int main(int argc, char **argv)
{ 
	FFTYPE(isignal) sig;
	FFTYPE(settings) set;
	FFTYPE(solution) sol;
	FFTYPE(variables) var;

	char *s = "../data/default_problem";
	int i;
	double l1, l2;
	double t1, t2, t=0.0;
	struct timeval tim;
	int rep = 2000;

	/* get signal */
	if( read_problem(s,&sig) ){
	
		/* get settings */	
		#ifdef SINGLE
		set.epsilon = 1e-3;
		#endif
		#ifdef DOUBLE
		set.epsilon = 1e-6;
		#endif
		set.verbose = 0;
		
		init_variables( sig, set, &var, &sol);

		/* warm up */
		for( i = 0 ; i < 100; ++i){
			reset_problem( sig, &var);
			reset_variables( sig, &var);
			cpdip_slp_core( sig, set, var, sol);
			if( *sol.kp >= 29)
				printf("Warning. Did not converge. Iteration %d",*sol.kp);
		}

		/* run timing */
		gettimeofday(&tim, NULL);
		l1 = tim.tv_sec + (tim.tv_usec/1000000.0);
		for( i = 0 ; i < rep; ++i){
			gettimeofday(&tim, NULL);
			t1 = tim.tv_sec + (tim.tv_usec/1000000.0);
			reset_problem( sig, &var);
			reset_variables( sig, &var);
			cpdip_slp_core( sig, set, var, sol);
			gettimeofday(&tim, NULL);
			t2 = tim.tv_sec + (tim.tv_usec/1000000.0);
			t += t2-t1;
			if( *sol.kp >= 29)
				printf("Warning. Did not converge. Iteration %d",*sol.kp);
		}
		gettimeofday(&tim, NULL);
		l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
		printf("Average time : %2.5f [s]\n", (l2-l1)/rep);
		printf("Average time via individual times: %2.5f [s]\n", t/rep);

		free_variables( var, sol);
		free_problem( sig);

		return 0;
	}
	else
		return 1;
}
