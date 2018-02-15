% Solves the Sparse Linear Prediction 
% using two different methods
% and compare objective

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, June 6, 2012

n = 6;
k = 4;
N = 10;
gamma = 0.1;

rpt = 1000;

RandStream.setDefaultStream(RandStream('mt19937ar','seed',1004397));
[B,A] = butter(N,0.2);
xs = filter(B,A,randn(n+2*k,1));
X = toeplitz(xs(k:n+k-1),xs(k:-1:1));
x = xs(k+1:n+k);
xp = xs(k:n+k-1);
f = @(a) norm(x-X*a,1)+gamma*norm(a,1);

keyboard
% --- CVX ---
for ii=1:rpt/10
cvx_begin quiet
variable a_cvx(k)
minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
cvx_end
end

tic 
for ii=1:rpt/10
cvx_begin quiet
variable a_cvx(k)
minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
cvx_end
end
t_cvx = toc;

% --- Double precision ---
for ii=1:rpt
	a_cip_double = pdip_slp(x,X,gamma,1e-7,0);
end

tic
for ii=1:rpt
a_cip_double = pdip_slp(x,X,gamma,1e-7,0);
end
t_cip_double = toc;

% --- Single precision ---
x_single = single(x);
X_single = single(X);
gamma_single = single(gamma);

for ii=1:rpt
	a_cip_single = pdip_slp_dual(x_single,X_single,gamma_single,1e-3,0);
end

tic
for ii=1:rpt
a_cip_single = pdip_slp_dual(x_single,X_single,gamma_single,1e-3,0);
end

t_cip_single = toc;

fprintf('\nComparison:\n')
fprintf('t_cvx = %.4f   t_cip_double = %.4f t_cip_single = %.4f \n',t_cvx/rpt, ...
		t_cip_double/rpt, t_cip_single/rpt);

fprintf('f(a_cvx) = %.5f  f(a_cip_double) = %.5f   f(a_cip_single) = %.5f\n',f(a_cvx),f(a_cip_double),f(a_cip_single));
