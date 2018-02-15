% Solves a small instances of the 
% sparse linear prediction problem
% 
% Mainly used for debugging/verfication 
% of the matlab / c++ -implementation

% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, June 6, 2012


addpath ../mlib
addpath ../data_process

[X x gamma] = readProblem('../data/default_problem');
[m,n] = size(X);

fprintf('-------------- Matlab ----------------\n')
[a_matlab opt2]= pdip_slp_primal(x,X,gamma,1e-6,1);
fprintf('\n');

fprintf('-------------- MEX C++ ---------------\n')
set.verbose = 1;set.type = 'primal';
[a_mex opt]= slp(x,X,gamma,set);
fprintf('\n');

%fprintf('-------------- CVX -------------------\n')
%cvx_begin
%variables a_cvx(n)
%minimize norm(x-X*a_cvx,1)+gamma*norm(a_cvx,1)
%cvx_end
%fprintf('\n');
%f = @(a,X,x,gamma) norm(x-X*a,1)+gamma*norm(a,1);

fprintf('--------- Solutions ------------------\n')
fprintf('    Matlab      Mex  ')
[a_matlab a_mex]