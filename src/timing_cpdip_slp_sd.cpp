#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "fcore.h"

int main(int argc, char **argv)
{
 	isignal_single sig_single;
	settings_single set_single;
	solution_single sol_single;
	variables_single var_single;

	isignal_double sig_double;
	settings_double set_double;
	solution_double sol_double;
	variables_double var_double;

	char *s = "../data/default_problem";
	int i;
	double l1, l2;
	double t1, t2, t=0.0;
	struct timeval tim;
	int rep = 2000;

	/* get signal */
	if( read_problem(s,&sig_single) ){
	
		/* get settings */	
		set_single.epsilon = 1e-3;
		set_single.verbose = 0;

		set_double.epsilon = 1e-6;
		set_double.verbose = 0;

		init_problem( sig_single.m, sig_single.n, &sig_double);
		copy_problem( &sig_single, &sig_double);
		
		init_variables( sig_single, set_single, &var_single, &sol_single);
		init_variables( sig_double, set_double, &var_double, &sol_double);		

		/* warm up */
		for( i = 0 ; i < 100; ++i){

			/*reset and run single */
			reset_problem( sig_single, &var_single);
			reset_variables( sig_single, &var_single);
			cpdip_slp_core( sig_single, set_single, var_single, sol_single);

			/*reset and run double */
			reset_problem( sig_double, &var_double);
			copy_variables( &sig_double, &var_single, &var_double);
			cpdip_slp_core( sig_double, set_double, var_double, sol_double);

			if( *sol_double.kp >= 29)
				printf("Warning. Did not converge. Iteration %d\n",*sol_double.kp);
		}

		/* run timing */
		gettimeofday(&tim, NULL);
		l1 = tim.tv_sec + (tim.tv_usec/1000000.0);

		for( i = 0 ; i < rep; ++i){
			gettimeofday(&tim, NULL);
			t1 = tim.tv_sec + (tim.tv_usec/1000000.0);

			/*reset and run single */
			reset_variables( sig_single, &var_single);
			cpdip_slp_core( sig_single, set_single, var_single, sol_single);

			/*reset and run double */
			copy_variables( &sig_double, &var_single, &var_double);
			//reset_variables( sig_double, &var_double);
			cpdip_slp_core( sig_double, set_double, var_double, sol_double);

			gettimeofday(&tim, NULL);
			t2 = tim.tv_sec + (tim.tv_usec/1000000.0);
			t += t2-t1;

			if( *sol_double.kp >= 29)
				printf("Warning. Did not converge. Iteration %d\n",*sol_double.kp);
		}


		gettimeofday(&tim, NULL);
		l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
		printf("Average time : %2.4f [s]\n", (l2-l1)/rep);
		printf("Average time via individual times: %2.4f [s]\n", t/rep);

		/* free variables */
		free_variables( var_single, sol_single);
		free_variables( var_double, sol_double);

		free_problem( sig_single);
		free_problem( sig_double);

		return 0;
	}
	else
		return 1;
}
