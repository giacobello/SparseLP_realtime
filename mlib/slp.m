function [alpha, info] = slp( x, X, gamma, varargin)
% slp - Sparse Linear Prediction
%
% [alpha, info] = slp( x, X, gamma, varargin)
%
% This function solves the problem
% 
% minimize ||x - X alpha||_1 + gamma ||alpha||_1
%
% for alpha. 
%
% To use this function you will need to install Math Kernel 
% Library (MKL) and compile the algorithms using the script 
% make_linux.m and make_windows.m depending on your system. 
%
% INPUTS:
%   X     : A real matrix size m x n
%   x     : A real vector size m
%   gamma : A positive (or zero) scalar
%
%settings : (Optional) A struct with some of the following fields
%
%      type: String 'primal' (default) or 'dual'. The algorithm
%           will use a primal or dual based method of solving the
%           the problem. In general the smaller n compared to m,
%           the more favorable is the primal method, and vice versa.
%           
%      prec: String 'double' (default) or 'single-double'. Double
%           will run the entire algorithm in double precision
%           where-as single-double will start the algorithm in single 
%           precision and shift to double-precision when an epsilon
%           accuracy of ep_single=1e-3 (default) is reached
%     
%       eps: Accuracy of the solution. Default 1e-6.
%
%eps_single: Accuracy where we shift from single to double precision
%           if prec == 'single-double'. Note that if eps_single <= eps
%           then the algorithm will run only in single precision.
%
% verbose  : Writes the iteration on screen if true. Default false.
%
% OUTPUTS:
%   alpha : Epsilon accuracy solution of optimization problem. 
%  
%   info  : A struct with the following fiels
%   
%    solved: True if solved to the requested accuracy, false otherwise
%    k     : Number of iterations
%
% AUTHOR:
%   Post-doc Tobias L. Jensen, Aalborg University, Denmark.
%     E-mail: tlj@es.aau.dk
%
% VERSION HISTORY:
%   0.1 [20-DEC-2012] Written the help part.
%
% COPYRIGHT:
%   2012 Tobias L. Jensen
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE?2.0
%
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
%   implied. See the License for the specific language governing
%   permissions and limitations under the License.


if length(varargin) == 1
	settings = varargin{1};
elseif length(varargin) > 1
	fprintf('slp: incorrect number of inputs. See slp.m\n');
	return
end

settings.def = 0;
% Set default values
if( ~isfield(settings,'type') )
	settings.type = 'primal';
end

if( ~isfield(settings,'prec') )
	settings.prec = 'double';
end

if( ~isfield(settings,'eps') )
	settings.eps = 1e-6;
end

if( ~isfield(settings,'eps_single') )
	settings.eps_single = 1e-3;
end

if( ~isfield(settings,'verbose') )
	settings.verbose = false;
end

% Validate input
if strcmp(settings.prec,'double')
	prec_id = 1;
elseif strcmp(settings.prec,'single-double')
	prec_id = 2;
else
	fprintf('No such prec %s. See help slp.\n',settings.prec);
	return
end

if size(X,1) ~= length(x)
	fprintf('Size of X and x does not match. See help slp.\n');
	return
end

% Fork of for primal and dual based algorithm
if strcmp(settings.type,'primal')
	% Primal algorithm
	
	[alpha, info.solved, info.k] = slp_primal( X, x, gamma, prec_id, ...
			settings.eps, settings.eps_single, settings.verbose);
	
elseif strcmp(settings.type,'dual')
	%Dual algorithm
	
	[alpha, info.solved, info.k] = slp_dual( X, x, gamma, prec_id, ...
			settings.eps, settings.eps_single, settings.verbose);
	
else
	fprintf('No such type %s. See help slp.\n',settings.type);
	return
end

