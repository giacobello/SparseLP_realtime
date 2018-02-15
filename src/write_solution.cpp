#include <stdio.h>
#include <stdlib.h>

#include "ftype.h"

#ifdef PRIMAL

#if defined(SD) && defined(SINGLE)
#include "core_primal_single.h"
#elif defined(SD) && defined(DOUBLE)
#include "core_primal_double.h"
#else
#include "core_primal.h"
#endif

#else

#if defined(SD) && defined(SINGLE)
#include "core_dual_single.h"
#elif  defined(SD) && defined(DOUBLE)
#include "core_dual_double.h"
#else
#include "core_dual.h"
#endif

#endif

/*writes the solution and timing to the file f*/
int write_solution(char *fn, int n, FFTYPE(solution) sol, FTYPE time){
	FILE *f;
	int i;

	f = fopen(fn,"w");
	if( f != NULL){
		fprintf( f, "%d\n", *sol.kp);
		fprintf( f, "%d\n", *sol.status);
		fprintf( f, "%lf\n", (double) time);
		for( i = 0 ; i < n ; i++)
			fprintf( f, "%lf\n", (double) sol.a[i]);
		
		fclose(f);
		return 1;
	}
	else{
		printf("Could not open file %s in write_solution\n",fn);
		return 0;
	}
}
