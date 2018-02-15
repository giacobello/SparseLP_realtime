function [k status time alpha] = readSolution(filename)

data = dlmread(filename,'\t');
k = data(1);
status = data(2);
time = data(3);
alpha = data(4:end);