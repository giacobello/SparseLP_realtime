% Solves the Sparse Linear Prediction 
% using two different methods
% and measures the timing
%
% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, Jan 15, 2013


addpath ../mlib
addpath ../data_process

problem = '../data/problem260100_44';
[X x gamma] = readProblem(problem);
[m,n] = size(X);

rep = 2000;

% Warm up
for r=1:rep
	[a_matlab info_matlab]= pdip_slp_primal(x,X,gamma,1e-6,0);
end

%Measure
tic;
for r=1:rep
	[a_matlab info_matlab]= pdip_slp_primal(x,X,gamma,1e-6,0);
end
t_matlab=toc;

% Warm up
for r=1:rep
	[a_mex info_mex] = slp(x,X,gamma);
end

%Measure
tic;
for r=1:rep
	[a_mex info_mex] = slp(x,X,gamma);
end
t_mex=toc;

fprintf('Timings for problem %s\n',problem)
fprintf('Matlab %f [s]\n',t_matlab/rep);
fprintf('Mex %f [s]\n',t_mex/rep);
