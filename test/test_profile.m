% Solves the Sparse Linear Prediction 
% problem and observe the profiling

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, June 6, 2012

n = 160;
k = 100;
N = 10;
gamma = 10.0;
rpt = 500;

% warm up

for i = 1:rpt
	[B,A] = butter(N,0.2);
	xs = filter(B,A,randn(n+2*k,1));
	X = toeplitz(xs(k:n+k-1),xs(k:-1:1));
	x = xs(k+1:n+k);
	
	a = pdip_slp_primal(x,X,gamma,1e-6,0);
end

% run
profile on
for i = 1:rpt
	[B,A] = butter(N,0.2);
	xs = filter(B,A,randn(n+2*k,1));
	X = toeplitz(xs(k:n+k-1),xs(k:-1:1));
	x = xs(k+1:n+k);
	
	a = pdip_slp_primal(x,X,gamma,1e-6,0);
end
profile off


profile report
