function [X x gamma] = readProblem(filename)

data = dlmread(filename,'\t');
m = data(1);
n = data(2);
gamma = data(3);
X = data(4:4+m*n-1); X = reshape(X,m,n);
x = data(4+m*n:end);