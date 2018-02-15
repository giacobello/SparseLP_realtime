# Call like
# make benchmark_small P=problem5010 NF=487 M=50 N=10

OVERRIDING=default

STOP = dropbox stop
START = dropbox start
NICE = nice --20 #standard user can set to -20

MATLAB = /home/tlj/Programs/R2011b/bin/matlab

benchmark_small: benchmark_primal benchmark_dual benchmark_primal_sd benchmark_dual_sd benchmark_cvxgen benchmark_Matlab 

benchmark_large: benchmark_primal benchmark_dual benchmark_primal_sd benchmark_dual_sd benchmark_Matlab

benchmark_Matlab:
	$(STOP)
	@echo "benchmark_signal('$(P)',$(NF))" > matlab_command.m
	$(NICE) $(MATLAB) -nodesktop < matlab_command.m
	$(START)

benchmark_primal:
	$(STOP)
	$(NICE) ../src/benchmark_signal_primal_double $(P) sol_$(P)_primal $(NF)
	$(START)

benchmark_dual:
	$(STOP)
	$(NICE) ../src/benchmark_signal_dual_double $(P) sol_$(P)_dual $(NF)
	$(START)	

benchmark_primal_sd:
	$(STOP)
	$(NICE) ../src/benchmark_signal_primal_sd $(P) sol_$(P)_primal_sd $(NF)
	$(START)

benchmark_dual_sd:
	$(STOP)
	$(NICE) ../src/benchmark_signal_dual_sd $(P) sol_$(P)_dual_sd $(NF)
	$(START)

benchmark_cvxgen:
	$(STOP)
	$(NICE) ../src/benchmark_signal_cvxgen $(P) sol_$(P)_cvxgen $(NF) $(M) $(N)
	$(START)

benchmark_cvxgen_blas:
	$(STOP)
	$(NICE) ../src/benchmark_signal_cvxgen_blas $(P) sol_$(P)_cvxgen_blas $(NF) $(M) $(N)
	$(START)