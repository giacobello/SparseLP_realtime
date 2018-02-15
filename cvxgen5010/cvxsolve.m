% Produced by CVXGEN, 2012-11-12 06:39:09 -0800.
% CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: cvxsolve.m.
% Description: Solution file, via cvx, for use with sample.m.

function [vars, status] = cvxsolve(params, settings)

X = params.X;
gamma = params.gamma;
x = params.x;

cvx_begin
  % Caution: automatically generated by cvxgen. May be incorrect.
  variable t(50, 1);
  variable s(10, 1);
  variable alpha(10, 1);

  minimize(sum(t) + gamma*sum(s));
  subject to
    -s <= alpha;
    alpha <= s;
    -t <= x - X*alpha;
    x - X*alpha <= t;
cvx_end

vars.alpha = alpha;
vars.s = s;
vars.t = t;
status.cvx_status = cvx_status;

% Provide a drop-in replacement for csolve.
status.optval = cvx_optval;
status.converged = strcmp(cvx_status, 'Solved');
