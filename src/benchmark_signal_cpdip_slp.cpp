#include <stdlib.h>
#include <stdio.h>

#include "ftype.h"

#include "fcore.h"

#include <sys/time.h>

int write_solution(char *fn, int n, FFTYPE(solution) sol, FTYPE time);

int main(int argc, char **argv)
{ 
	FFTYPE(isignal) sig;
	FFTYPE(settings) set;
	FFTYPE(solution) sol;
	FFTYPE(variables) var;

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

	/* get settings */	
	set.epsilon = 1e-6;
	set.verbose = 0;

	/* init variables based on the first problem */
	sprintf(f,"%s_%d",s,i);
	if( read_problem(f,&sig) ) {
		init_variables( sig, set, &var, &sol);

		/* run over all problems in the signal */
		for( i = 1 ; i <= no_of_frames ; i++ ){

			//printf("Run problem instance %d\n",i);

			/* get signal */
			sprintf(f,"%s_%d",s,i);
			if( read_problem(f,&sig) ){

				gettimeofday(&tim, NULL);
				l1 = tim.tv_sec + (tim.tv_usec/1000000.0);			
				for( r = 0 ; r < rpt; ++r){

					reset_problem( sig, &var);
					reset_variables( sig, &var);
					cpdip_slp_core( sig, set, var, sol);

					if( *sol.kp >= 29)
						printf("Warning. Did not converge. Iteration %d",*sol.kp);
				}
				gettimeofday(&tim, NULL);
				l2 = tim.tv_sec + (tim.tv_usec/1000000.0);
				
				sprintf(f,"%s_%d",v,i);
				write_solution( f, sig.n, sol, (l2-l1)/rpt);
				//printf("Average time : %2.5f [s]\n", (l2-l1)/rpt);
				free_problem( sig);
			}
			else
				return 0;
		}
		free_variables( var, sol);
	}
	else 
		return 0;
}
