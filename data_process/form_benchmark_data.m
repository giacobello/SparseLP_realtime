function form_benchmark_data( N, K, varargin)
% This function will form some test data, and insert it into the
% the data/ folder with the same data structure as used for 
% loading single problems in the C/C++ implementation
%
% Input:
%   N: frame size
%   K: predictor order
% Optional:
%   gamma: regularization parameter
%
% Tobias Lindstrøm Jensen, 2012, Aalborg University
% tlj@es.aau.dk

% Order and size
%N = 160;
%K = 100;

% Toeplitz (problem) sizes
m = K + N;
n = K;

%M = 256, D = 96;
%W = sin(pi*([0:D-1]'+0.5)/2/D);
%W = [W;ones(M-2*D,1);flipud(W);];

% Offset
offset = 1000;

% Gamma (regularization)
if length(varargin) == 1
	gamma = varargin{1};
else
	gamma = 0.1;
end

% Filename
filename = '../data/fa.wav'

[xs fs nbits] = wavread(filename);

nFrames = floor((length(xs)-offset)/N);

for k = 1:nFrames
	
	s = [0 xs((1+(k-1)*N):k*N)' zeros(1,K-1)];
	x = [s(2:end) 0]';
	X = toeplitz(s,zeros(K,1)');
	
	sn = sprintf('../data/problem%d%d_%d',m,n,k);
	writeProblem(sn,X,x,gamma);
	fprintf('%1.2f procent done. Written to %s\r', 100*k/nFrames,sn);
end

fprintf('Last index is k=%d\n',nFrames);
