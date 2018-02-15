% csolve  Solves a custom quadratic program very rapidly.
%
% [vars, status] = csolve(params, settings)
%
% solves the convex optimization problem
%
%   minimize(sum(t) + gamma*sum(s))
%   subject to
%     -s <= alpha
%     alpha <= s
%     -t <= x - X*alpha
%     x - X*alpha <= t
%
% with variables
%    alpha  10 x 1
%        s  10 x 1
%        t  50 x 1
%
% and parameters
%        X  50 x 10
%    gamma   1 x 1    positive
%        x  50 x 1
%
% Note:
%   - Check status.converged, which will be 1 if optimization succeeded.
%   - You don't have to specify settings if you don't want to.
%   - To hide output, use settings.verbose = 0.
%   - To change iterations, use settings.max_iters = 20.
%   - You may wish to compare with cvxsolve to check the solver is correct.
%
% Specify params.X, ..., params.x, then run
%   [vars, status] = csolve(params, settings)


% Produced by CVXGEN, 2012-11-12 06:39:10 -0800.
% CVXGEN is Copyright (C) 2006-2011 Jacob Mattingley, jem@cvxgen.com.
% The code in this file is Copyright (C) 2006-2011 Jacob Mattingley.
% CVXGEN, or solvers produced by CVXGEN, cannot be used for commercial
% applications without prior written permission from Jacob Mattingley.

% Filename: csolve.m.
% Description: Help file for the Matlab solver interface.
