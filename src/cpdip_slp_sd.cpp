#include <cstdlib>
#include <cstdio>

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

	/* get signal */
	if( read_problem(s,&sig_single) ){
	
		/* get settings */	
		set_single.epsilon = 1e-3;
		set_single.verbose = 1;

		set_double.epsilon = 1e-6;
		set_double.verbose = 1;

		init_problem( sig_single.m, sig_single.n, &sig_double);
		copy_problem( &sig_single, &sig_double);
		
		init_variables( sig_single, set_single, &var_single, &sol_single);
		init_variables( sig_double, set_double, &var_double, &sol_double);		

		reset_problem( sig_single, &var_single);
		reset_problem( sig_double, &var_double);

		/*reset and run single */
		reset_variables( sig_single, &var_single);
		cpdip_slp_core( sig_single, set_single, var_single, sol_single);

		/*reset and run double */
		copy_variables( &sig_double, &var_single, &var_double);
		//reset_variables( sig_double, &var_double);
		cpdip_slp_core( sig_double, set_double, var_double, sol_double);
		
		/* print solution */
		//print_vector( sig_double.n, sol_double.a, "sol_double" );

		free_variables( var_single, sol_single);
		free_variables( var_double, sol_double);

		free_problem( sig_single);
		free_problem( sig_double);

		return 1;
	}
	else
		return 0;
}
