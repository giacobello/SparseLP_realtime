function writeProblem(filename,X,x,gamma)

[m,n] = size(X);

dlmwrite(filename,m,'precision','%d','delimiter','\t');
dlmwrite(filename,n,'-append','precision','%d','delimiter','\t');
dlmwrite(filename,gamma,'-append','precision','%.8f','delimiter','\t');
dlmwrite(filename,X(:),'-append','precision','%.8f','delimiter','\t');
dlmwrite(filename,x(:),'-append','precision','%.8f','delimiter','\t');