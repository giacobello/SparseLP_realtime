function becnhmark_signal( problem, no_of_frames)
% Benchmark several methods for computing
% sparse linear prediction using Matlab
% based methods.
%
% problem = string, name of problem
% no_of_frames = integer, number of frames to process
%
% Methods:
% Mosek (Mosek "overwrites" linprog in Matlab)
% CVX+SeDuMi
% pdip_slp_primal
% pdip_slp_dual
%
% with warm-up. Saves the result in a file. 
% Use plot_benchmark_signal.m to plot the results
%
% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, November 7, 2012

addpath ../mlib
addpath ../data_process

problem_name = sprintf('../data/%s', problem);
f = @(alpha,x,X,gamma) norm(x-X*alpha,1)+gamma*norm(alpha,1);

[X x gamma] = readProblem(sprintf('%s_%d',problem_name,1));
nt = size(X,2);

% Benchmark settings
rpt = 100;
rpt_CVX = 100;

epsilon = 1e-6;

% For results
timing_Mosek = zeros(no_of_frames,1);
timing_CVX_SeDuMi = zeros(no_of_frames,1);
timing_pdip_slp_primal = zeros(no_of_frames,1);
timing_pdip_slp_dual = zeros(no_of_frames,1);

result_Mosek = zeros(nt,no_of_frames);
result_CVX_SeDuMi = zeros(nt,no_of_frames);
result_pdip_slp_primal = zeros(nt,no_of_frames);
result_pdip_slp_dual = zeros(nt,no_of_frames);

obj_Mosek = zeros(no_of_frames,1);
obj_CVX_SeDuMi = zeros(no_of_frames,1);
obj_pdip_slp_primal = zeros(no_of_frames,1);
obj_pdip_slp_dual = zeros(no_of_frames,1);

% ------------ pdip_slp_primal ------------
fprintf('Runs for Mosek\n');
% Process ones for warm-up
for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));
	[m,n] = size(X);

	A = [-X;gamma*eye(n)]; 
	b = [-x;zeros(n,1)];
	
	%G = sparse([A, -eye(m+n); -A, -eye(m+n)]);
	G = [A, -eye(m+n); -A, -eye(m+n)];
	h = [b;-b];
	c = [zeros(n,1);ones(m+n,1)];
	
	xt = linprog(c,G,h);
	alpha = xt(1:n);
end

for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));
	[m,n] = size(X);
	
	t = tic;
	for r = 1:rpt
			A = [-X;gamma*eye(n)]; 
			b = [-x;zeros(n,1)];
				
				G = [A, -eye(m+n); -A, -eye(m+n)];
				h = [b;-b];
				c = [zeros(n,1);ones(m+n,1)];
				
				xt = linprog(c,G,h);
				alpha = xt(1:n);
	end
	time = toc(t);
	
	timing_Mosek(ii) = time/rpt;
	result_Mosek(:,ii) = alpha;
	obj_Mosek(ii) = f(alpha,x,X,gamma);
	fprintf('Run procent done %2.2f\r',100*ii/no_of_frames);
end

% ------------ pdip_slp_primal ------------
fprintf('Runs for pdip_slp_primal\n');
% Process ones for warm-up
for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));
	[m,n] = size(X);
		
	alpha = pdip_slp_primal(x,X,gamma,epsilon,0);
end

for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));
	[m,n] = size(X);
		
	t = tic;
	for r = 1:rpt
		alpha = pdip_slp_primal(x,X,gamma,epsilon,0);
	end
	time = toc(t);
	
	timing_pdip_slp_primal(ii) = time/rpt;
	result_pdip_slp_primal(:,ii) = alpha;
	obj_pdip_slp_primal(ii) = f(alpha,x,X,gamma);
	fprintf('Run procent done %2.2f\r',100*ii/no_of_frames);
end

% ------------ pdip_slp_dual ------------
fprintf('Runs for pdip_slp_dual\n');
% Process ones for warm-up
for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));	
	[m,n] = size(X);
	
	alpha = pdip_slp_dual(x,X,gamma,epsilon,0);
end

for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));

	t = tic;
	for r = 1:rpt
		alpha = pdip_slp_dual(x,X,gamma,epsilon,0);
	end
	time = toc(t);
	
	timing_pdip_slp_dual(ii) = time/rpt;
	result_pdip_slp_dual(:,ii) = alpha;
	obj_pdip_slp_dual(ii) = f(alpha,x,X,gamma);
	fprintf('Run procent done %2.2f\r',100*ii/no_of_frames);
end

% ------------ CVX+SeDuMi ------------
fprintf('Runs for CVX+SeDuMi\n');
% Process ones for warm-up
for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));
	[m,n] = size(X);
	
	cvx_begin
	cvx_quiet(true)
	variable a_cvx(n)
	minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
	cvx_end
end

for ii = 1:no_of_frames
	[X x gamma] = readProblem(sprintf('%s_%d',problem_name,ii));
	[m,n] = size(X);
	
	t = tic;
	for r = 1:rpt_CVX
		cvx_begin
		cvx_quiet(true)
		variable a_cvx(n)
		minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
		cvx_end
	end
	time = toc(t);
	
	timing_CVX_SeDuMi(ii) = time/rpt_CVX;
	result_CVX_SeDuMi(:,ii) = a_cvx;
	obj_CVX_SeDuMi(ii) = f(a_cvx,x,X,gamma);
	fprintf('Run procent done %2.2f\r',100*ii/no_of_frames);
end


save(sprintf('../data/sol_%s_Matlab',problem),'n','problem','problem_name','no_of_frames','rpt','epsilon','timing_Mosek','timing_CVX_SeDuMi','timing_pdip_slp_primal','timing_pdip_slp_dual','result_Mosek','result_CVX_SeDuMi','result_pdip_slp_primal','result_pdip_slp_dual','obj_Mosek','obj_CVX_SeDuMi','obj_pdip_slp_primal','obj_pdip_slp_dual');
