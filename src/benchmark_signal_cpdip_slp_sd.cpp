#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

#include "fcore.h"

int write_solution(char *fn, int n, solution_double sol, double time);

int main(int argc, char **argv){
 	isignal_single sig_single;
	settings_single set_single;
	solution_single sol_single;
	variables_single var_single;

	isignal_double sig_double;
	settings_double set_double;
	solution_double sol_double;
	variables_double var_double;

	/* get settings */	
	set_single.epsilon = 1e-3;
	set_single.verbose = 0;
	
	set_double.epsilon = 1e-6;
	set_double.verbose = 0;

	if( argc != 4 ){
		printf("Number of input should be 3: problem_name solution_name no_of_frames\n");
		return 0;
	}

	char s [100];
	sprintf(s,"../data/%s",argv[1]);
	char v [100];
	sprintf(v,"../data/%s",argv[2]);
	
	
	char f [100];
	int no_of_frames;
	sscanf(argv[3],"%d",&no_of_frames);

	int rpt = 100;

	int i=1;
	int r;
	double l1, l2;
	struct timeval tim;


	/* run over all problems in the signal */
	for( i = 1 ; i <= no_of_frames ; i++ ){

		//printf("Run problem instance %d\n",i);
		
		/* get signal */
		sprintf(f,"%s_%d",s,i);
		if( read_problem(f,&sig_single) ){
			
			init_problem( sig_single.m, sig_single.n, &sig_double);
			copy_problem( &sig_single, &sig_double);
			
			init_variables( sig_single, set_single, &var_single, &sol_single);
			init_variables( sig_double, set_double, &var_double, &sol_double);		
			
			gettimeofday(&tim, NULL);
			l1 = tim.tv_sec + (tim.tv_usec/1000000.0);			
			for( r = 0 ; r < rpt; ++r){
				
				/*reset and run single */
				reset_problem( sig_single, &var_single);
				reset_variables( sig_single, &var_single);
				cpdip_slp_core( sig_single, set_single, var_single, sol_single);
				
				/*copy and run double */
				reset_problem( sig_double, &var_double);
				copy_variables( &sig_double, &var_single, &var_double);
				cpdip_slp_core( sig_double, set_double, var_double, sol_double);

				if( *sol_double.kp >= 29)
					printf("Warning. Did not converge. Iteration %d, problem %d\n",*sol_double.kp,i);
			}
			gettimeofday(&tim, NULL);
			l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
			
			sprintf(f,"%s_%d",v,i);
			write_solution( f, sig_double.n, sol_double, (l2-l1)/rpt);
			//printf("Average time : %2.5f [s]\n", (l2-l1)/rpt);
			
			/* free variables */
			free_variables( var_single, sol_single);
			free_variables( var_double, sol_double);
			
			free_problem( sig_single);
			free_problem( sig_double);
			
		}//end read problem
		else
			return 1;
		
	}//end loop problem
	return 0;
}
