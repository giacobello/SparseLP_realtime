#include <stdlib.h>

#include "ftype.h"
#include "fcore.h"

int main(int argc, const char **argv)
{ 
	FFTYPE(isignal) sig;
	FFTYPE(settings) set;
	FFTYPE(solution) sol;
	char* s = "../data/default_problem";
	FFTYPE(variables) var;

	/* get signal */
	if( read_problem(s,&sig) ){
	
		/* get settings */	
		#ifdef SINGLE
		set.epsilon = 1e-3;
		#endif
		#ifdef DOUBLE
		set.epsilon = 1e-6;
		#endif

		set.verbose = 1;
		
		init_variables( sig, set, &var, &sol);

		reset_problem( sig, &var);

		reset_variables( sig, &var);

		/* run the algorithm */
		cpdip_slp_core( sig, set, var, sol);

		free_variables( var, sol);

		free_problem( sig);

		return 0;
	}
	else
		return 1;
}
