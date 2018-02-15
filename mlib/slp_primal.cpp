#include <mex.h>

#include "ftype.h"
#include "fcore.h"

void copy_solution( int n, solution_double sol_double, mxArray *plhs[] );

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

	mxArray *Xm;
	int prec_id, verbose;

	isignal_double sig_double;
	settings_double set_double;
	solution_double sol_double;
	variables_double var_double;

	// Verify the number of inputs
	if(nrhs != 7)
		printf("Should contain 7 input parameters but has %i\n",nrhs);

	// Pass the data structures
	Xm = (mxArray*) prhs[0]; 
	sig_double.m = mxGetM(Xm); sig_double.n = mxGetN(Xm);

	sig_double.Y = mxGetPr(Xm);
	sig_double.y = mxGetPr( (mxArray*) prhs[1] ); 

	sig_double.gamma = (double) mxGetScalar( (mxArray*) prhs[2] );
	prec_id = (int) mxGetScalar( (mxArray*) prhs[3] );
	set_double.epsilon = (double) mxGetScalar( (mxArray*) prhs[4] );

	//set_single.epsilon = (double) mxGetScalar( (mxArray*) prhs[5] );
	verbose = (int) mxGetScalar( (mxArray*) prhs[6] );
	set_double.verbose = verbose; //set_single.verbose = verbose; 

	// Allocate
	init_variables( sig_double, set_double, &var_double, &sol_double);

	// Reset to default and calculate initial point 
	reset_problem( sig_double, &var_double);
	reset_variables( sig_double, &var_double);

	//// Run the main algorithm
	cpdip_slp_core( sig_double, set_double, var_double, sol_double);

	//// Allocate memory and assign output pointer
	plhs[0] = mxCreateDoubleMatrix(sig_double.n, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);

	////Copy solution to output
	copy_solution( sig_double.n, sol_double, plhs);

	//// Free variables
	free_variables( var_double, sol_double);

}// exit


void copy_solution( int n, solution_double sol_double, mxArray *plhs[] ){
	int i;
	double *a = mxGetPr( plhs[0] );
	double *s = mxGetPr( plhs[1] );
	double *kp = mxGetPr( plhs[2] );

	for( i = 0 ; i < n ; i++ )
		a[i] = sol_double.a[i];

	s[0] = (double) *sol_double.status;
	kp[0] = (double) *sol_double.kp;
}
