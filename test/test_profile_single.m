% Solves the Sparse Linear Prediction 
% problem and observe the profiling

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, June 6, 2012

n = 160;
k = 100;
N = 10;
gamma = 4.5;
rpt = 100;

% warm up
for i = 1:rpt
	[B,A] = butter(N,0.2);
	xs = filter(B,A,randn(n+2*k,1));
	X = toeplitz(xs(k:n+k-1),xs(k:-1:1));
	x = xs(k+1:n+k);
	
	x_single = single(x);
	X_single = single(X);
	gamma_single = single(gamma);

	a = pdip_slp(x_single,X_single,gamma_single,1e-2,0);
end

% run
profile on
for i = 1:rpt
	[B,A] = butter(N,0.2);
	xs = filter(B,A,randn(n+2*k,1));
	X = toeplitz(xs(k:n+k-1),xs(k:-1:1));
	x = xs(k+1:n+k);
	
	x_single = single(x);
	X_single = single(X);
	gamma_single = single(gamma);

	a = pdip_slp(x_single,X_single,gamma_single,1e-2,0);

end
profile off

profile report
