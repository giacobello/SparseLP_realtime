CC = g++
CFLAGS = -Wall -march=native -Os -g

#MKLROOT = /opt/intel/composer_xe_2011_sp1.10.319
MKLROOT = /opt/intel/composer_xe_2011_sp1.11.339

#Dynamic linking
#MKL = -I$(MKLROOT)/mkl/include/ -lm -L$(MKLROOT)/mkl/lib/intel64 -L$(MKLROOT)/compiler/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread

#Static linking
MKL =  -Wl,--start-group  $(MKLROOT)/mkl/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/mkl/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -L$(MKLROOT)/compiler/lib/intel64 -liomp5 -lpthread -lm

#LINK = -lblas -llapack
LINK = -DMKL $(MKL)

CVXGEN_DIR = /home/tlj/Dropbox/SparSig_local/SparseLP/cvxgen5010/
CVXGEN_DIR_BLAS = /home/tlj/Dropbox/SparSig_local/SparseLP/cvxgen5010_blas/

primal_single: tools.cpp tools.h ftype.h core_primal.h core_primal.cpp
	$(CC) -DSD -DSINGLE -E tools.h > tools_single.h
	$(CC) -DSD -DSINGLE -E tools.cpp > tools_single.cpp
	$(CC) -DSD -DSINGLE -E core_primal.h > core_primal_single.h
	$(CC) -DSD -DSINGLE -E core_primal.cpp > core_primal_single.cpp

primal_double: tools.cpp tools.h ftype.h core_primal.h core_primal.cpp
	$(CC) -DSD -DDOUBLE -E tools.h > tools_double.h
	$(CC) -DSD -DDOUBLE -E tools.cpp > tools_double.cpp
	$(CC) -DSD -DDOUBLE -E core_primal.h > core_primal_double.h
	$(CC) -DSD -DDOUBLE -E core_primal.cpp > core_primal_double.cpp

dual_single: tools.cpp tools.h ftype.h core_dual.h core_dual.cpp
	$(CC) -DSD -DSINGLE -E tools.h > tools_single.h
	$(CC) -DSD -DSINGLE -E tools.cpp > tools_single.cpp
	$(CC) -DSD -DSINGLE -E core_dual.h > core_dual_single.h
	$(CC) -DSD -DSINGLE -E core_dual.cpp > core_dual_single.cpp

dual_double: tools.cpp tools.h ftype.h core_dual.h core_dual.cpp
	$(CC) -DSD -DDOUBLE -E tools.h > tools_double.h
	$(CC) -DSD -DDOUBLE -E tools.cpp > tools_double.cpp
	$(CC) -DSD -DDOUBLE -E core_dual.h > core_dual_double.h
	$(CC) -DSD -DDOUBLE -E core_dual.cpp > core_dual_double.cpp

read_problem: read_problem.cpp write_solution.cpp
	$(CC) -DSD -DDOUBLE -E read_problem.cpp > read_problem_double.cpp
	$(CC) -DSD -DSINGLE -E read_problem.cpp > read_problem_single.cpp
	$(CC) -DSD -DDOUBLE -E write_solution.cpp > write_solution_double.cpp
	$(CC) -DSD -DSINGLE -E write_solution.cpp > write_solution_single.cpp

slp_primal_sd: Makefile primal_single primal_double read_problem cpdip_slp_sd.cpp fcore.h
	$(CC) $(CFLAGS) -DSD -DPRIMAL -o slp_primal_sd cpdip_slp_sd.cpp read_problem_single.cpp read_problem_double.cpp tools_single.cpp tools_double.cpp core_primal_single.cpp core_primal_double.cpp  $(LINK) 

slp_dual_sd: Makefile dual_single dual_double read_problem cpdip_slp_sd.cpp fcore.h
	$(CC) $(CFLAGS)  -DSD -DDUAL -o slp_dual_sd cpdip_slp_sd.cpp read_problem_single.cpp read_problem_double.cpp tools_single.cpp tools_double.cpp core_dual_single.cpp core_dual_double.cpp  $(LINK) 

slp_primal_single: Makefile read_problem.cpp cpdip_slp.cpp core_primal.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DPRIMAL -DSINGLE -o slp_primal_single cpdip_slp.cpp read_problem.cpp tools.cpp core_primal.cpp $(LINK)

slp_dual_single: Makefile read_problem.cpp cpdip_slp.cpp core_dual.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DDUAL -DSINGLE -o slp_dual_single cpdip_slp.cpp read_problem.cpp tools.cpp core_dual.cpp $(LINK)

slp_primal_double: Makefile read_problem.cpp cpdip_slp.cpp core_primal.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DPRIMAL -DDOUBLE -o slp_primal_double cpdip_slp.cpp read_problem.cpp tools.cpp core_primal.cpp $(LINK)

slp_dual_double: Makefile read_problem.cpp cpdip_slp.cpp core_dual.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DDUAL -DDOUBLE -o slp_dual_double cpdip_slp.cpp read_problem.cpp tools.cpp core_dual.cpp $(LINK)

timing_primal_single: Makefile read_problem.cpp timing_cpdip_slp.cpp core_primal.cpp fcore.h ftype.h
	$(CC) $(CFLAGS) -DPRIMAL -DSINGLE -o timing_primal_single timing_cpdip_slp.cpp read_problem.cpp tools.cpp core_primal.cpp $(LINK)

timing_dual_single: Makefile read_problem.cpp timing_cpdip_slp.cpp core_dual.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DDUAL -DSINGLE -o timing_dual_single timing_cpdip_slp.cpp read_problem.cpp tools.cpp core_dual.cpp $(LINK)

timing_primal_double: Makefile read_problem.cpp timing_cpdip_slp.cpp fcore.h ftype.h
	$(CC) $(CFLAGS) -DPRIMAL -DDOUBLE -o timing_primal_double timing_cpdip_slp.cpp read_problem.cpp tools.cpp core_primal.cpp $(LINK)

timing_dual_double: Makefile read_problem.cpp timing_cpdip_slp.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DDUAL -DDOUBLE -o timing_dual_double timing_cpdip_slp.cpp read_problem.cpp tools.cpp core_dual.cpp $(LINK)

timing_primal_sd: Makefile primal_single primal_double read_problem timing_cpdip_slp_sd.cpp fcore.h
	$(CC) $(CFLAGS) -g  -DPRIMAL -DSD -o timing_primal_sd timing_cpdip_slp_sd.cpp read_problem_single.cpp read_problem_double.cpp tools_single.cpp tools_double.cpp core_primal_single.cpp core_primal_double.cpp  $(LINK) 

timing_dual_sd: Makefile dual_single dual_double read_problem timing_cpdip_slp_sd.cpp fcore.h
	$(CC) $(CFLAGS) -g  -DDUAL -DSD -o timing_dual_sd timing_cpdip_slp_sd.cpp read_problem_single.cpp read_problem_double.cpp tools_single.cpp tools_double.cpp core_dual_single.cpp core_dual_double.cpp  $(LINK) 

benchmark_signal_primal_double: Makefile read_problem.cpp benchmark_signal_cpdip_slp.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DPRIMAL -DDOUBLE -o benchmark_signal_primal_double benchmark_signal_cpdip_slp.cpp read_problem.cpp write_solution.cpp tools.cpp core_primal.cpp $(LINK)

benchmark_signal_dual_double: Makefile read_problem.cpp benchmark_signal_cpdip_slp.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DDUAL -DDOUBLE -o benchmark_signal_dual_double benchmark_signal_cpdip_slp.cpp read_problem.cpp write_solution.cpp tools.cpp core_dual.cpp $(LINK)

benchmark_signal_primal_sd: Makefile primal_single primal_double read_problem benchmark_signal_cpdip_slp_sd.cpp fcore.h ftype.h
	$(CC) $(CFLAGS) -DPRIMAL -DSD -o benchmark_signal_primal_sd benchmark_signal_cpdip_slp_sd.cpp read_problem_double.cpp read_problem_single.cpp write_solution_double.cpp tools_single.cpp tools_double.cpp core_primal_single.cpp core_primal_double.cpp $(LINK)

benchmark_signal_dual_sd: Makefile dual_single dual_double read_problem benchmark_signal_cpdip_slp_sd.cpp fcore.h ftype.h
	$(CC) $(CFLAGS)  -DDUAL -DSD -o benchmark_signal_dual_sd benchmark_signal_cpdip_slp_sd.cpp read_problem_double.cpp read_problem_single.cpp write_solution_double.cpp tools_single.cpp tools_double.cpp core_dual_single.cpp core_dual_double.cpp $(LINK)

benchmark_signal_cvxgen: Makefile benchmark_signal_cvxgen.c
	gcc $(CFLAGS) -I$(CVXGEN_DIR) -c benchmark_signal_cvxgen.c -o benchmark_signal_cvxgen.o

	gcc $(CFLAGS) -I$(CVXGEN_DIR) $(CVXGEN_DIR)solver.h $(CVXGEN_DIR)solver.o $(CVXGEN_DIR)matrix_support.o $(CVXGEN_DIR)ldl.o $(CVXGEN_DIR)util.o benchmark_signal_cvxgen.o -lm -o benchmark_signal_cvxgen

benchmark_signal_cvxgen_blas: Makefile benchmark_signal_cvxgen_blas.cpp
	$(CC) $(CFLAGS) -I$(CVXGEN_DIR_BLAS) -c benchmark_signal_cvxgen_blas.cpp -o benchmark_signal_cvxgen_blas.o

	$(CC) $(CFLAGS) -I$(CVXGEN_DIR_BLAS) $(CVXGEN_DIR_BLAS)solver.h $(CVXGEN_DIR_BLAS)solver.o $(CVXGEN_DIR_BLAS)matrix_support.o $(CVXGEN_DIR_BLAS)ldl.o $(CVXGEN_DIR_BLAS)util.o benchmark_signal_cvxgen_blas.o -lm -o benchmark_signal_cvxgen_blas $(LINK)

benchmark: benchmark_signal_primal_double benchmark_signal_dual_double benchmark_signal_primal_sd benchmark_signal_dual_sd benchmark_signal_cvxgen

clean:
	rm -f *.o
	find . -perm +100 -type f -delete
	rm -f *.*~