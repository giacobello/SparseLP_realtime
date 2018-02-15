function [a opt] = pdip_slp_dual(y,Y,gamma,epsilon,verbose,varargin)
% Solves
% minimize   ||y-Ya||_1 + gamma ||a||_1
%
% to epsilon accuracy
%
% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, June 6, 2012

solver = 1;
COUNT_MAX = 5;

[m,n]=size(Y);

%Form new matrices

cls = superiorfloat(y);

one = ones(1,1,cls);

c = [ones(2*m,1);gamma*ones(2*n,1)];


%Am = [-eye(m) eye(m) -Y Y];
%A = @(x) -x(1:m)+x(m+1:2*m) + Y*(-x(2*m+1:2*m+n)+x(2*m+n+1:end));
%At = @(z) Atf(z,Y);

b = 2 * y;
% such that we have the primal and dual problems
% minimize c'*x
% subject to A*x=b, x>=0
%
% maximize b'*lambda
% subject to A'*lambda + s = c, s>=0
%

N = one*(2*m+2*n);
M = one*(m);

if solver == 1
	opts1.LT = true; 
	opts2.LT = true; opts2.TRANSA = true;
elseif solver == 2
	opts1.UT = true; opts1.TRANSA = true;
	opts2.UT = true; 

else
	fprintf('No such solver\n');return
end

if verbose
	fprintf('Ite.    pinfeas.      dinfeas       rel_gap\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Obtain initial settings
% As suggested by Nocedal & Wright
%xt = Am'*((Am*Am')\b);
%lt = (Am*Am')\(Am*c);
%st = c-Am'*lt;
%
%dx = max(-(3/2)*min(xt),0);
%ds = max(-(3/2)*min(st),0);
%
%xh = xt + dx*ones(size(xt));
%sh = st + ds*ones(size(st));
%
%dxh = 0.5* (xh'*sh)/(sum(sh));
%dsh = 0.5* (xh'*sh)/(sum(xh));
%
%xk = xh + dxh*ones(size(xh));
%sk = sh + dsh*ones(size(sh));
%lk = lt;

% Simple initial setting
xk = ones(N,1);
sk = ones(N,1);
lk = b;

for k=1:100
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%1. compute residual and evaluate stopping criteria
	%rb = A*xk - b;
	%rc = At*lk + sk -c;

	rb = Af(xk,Y,m,n) - b;
	
	
	if k ==  1
		Atlk = Atf(lk,Y); %First iteration calc via a function handle
	else
		Atlk =  Atlk+ad*Atdls; %Compute via precomputed values
	end
	
	rc = Atlk + sk - c;	
	
	mu = xk'*sk/N;

	pinfeas = norm(rb)/(one+norm(b));
	dinfeas = norm(rc)/(one+norm(c));
	relgap = (c'*xk-b'*lk)/(one+abs(c'*xk));

	if verbose
		fprintf('%2d       %.2e      %.2e      %.2e\n',k,pinfeas,dinfeas,relgap);
	end
	
	if max([relgap,pinfeas,dinfeas]) <= epsilon
		break;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%2. Calculate the affine scaling
	
	% KKT system
	%K = [zeros(N,N) A' eye(N);
	%	   A  zeros(M,M) zeros(M,N);
	%	   diag(sk) zeros(N,M) diag(xk)];
	%ra = [-rc;-rb;-xk.*sk];
	%da = K\ra;
	%dxa = da(1:N);
	%dsa = da(N+M+1:end);
	%dla = da(N+1:N+M);
	
	% Solved by
	%Ds = diag(xk./sk);
	%Kp = (A*Ds*A');

	skm1 = one./sk; % precalculate 
	
	%D1 = diag(xk(1:m).*skm1(1:m))+diag(xk(m+1:2*m).*skm1(m+1:2*m));
	%D2 = diag(xk(2*m+1:2*m+n).*skm1(2*m+1:2*m+n))+diag(xk(2*m+n+1:end).*skm1(2*m+n+1:end));
	%Kp = D1 + Y*D2*Y';
	xkskm1 = xk.*skm1;
	d1 = xkskm1(1:m)+xkskm1(m+1:2*m); %xk(1:m).*skm1(1:m)+xk(m+1:2*m).*skm1(m+1:2*m);
	d2 = xkskm1(2*m+1:2*m+n)+xkskm1(2*m+n+1:end);

	%dla = Kp\(-rb-A*diag(xk./sk)*rc+A*diag(1./sk)*(xk.*sk));

	%ra = (-rb-A*((xk.*skm1).*rc-xk));
	ra = (-rb-Af((xkskm1).*rc-xk,Y,m,n));


	if solver == 1
		Kpp = bsxfun(@times,d2,Y');
		Kp = Y*Kpp;

		Kp(1:M+1:end) = Kp(1:M+1:end)' + d1; %adds to the diagonal

		[C,p] = chol(Kp,'lower');
	
		count = 0;
		while p ~= 0
			if verbose
				fprintf('Cholesky failed in iteration %d. Add regularization\n',k);
			end
			
			Kp(1:M+1:end) = Kp(1:M+1:end)'+epsilon; %adds regularization to the diagonal
			[C,p] = chol(Kp,'lower');
			
			count = count + 1;
		
			if count > COUNT_MAX
				a = (-xk(2*m+1:2*m+n)+xk(2*m+n+1:2*m+2*n))*0.5;
				opt.k = k;
				opt.status = 'failed';
				return
			end
		
		end


		%v = C'\ra;
		%dla = C\v;

	elseif solver == 2
		Kpp = bsxfun(@times,sqrt(d2),Y');	
		Kp = [diag(sqrt(d1)); Kpp];
	
		Cp = qr(Kp,0);
		C = Cp(1:m,1:m);
	end

	v = linsolve(C,ra,opts1);
	dla = linsolve(C,v,opts2);

	%dsa = -rc-A'*dla;
	dsa = -rc-Atf(dla,Y);
	dxa = -xk-(xkskm1).*dsa;


	% Calculate the primal-dual stepsize
	I = find(dxa<0);
	axa = min([1.0;-xk(I)./dxa(I)]);
    
	I = find(dsa<0);
	asa = min([1.0;-sk(I)./dsa(I)]);
	
	mu_aff = (xk+axa*dxa)'*(sk+asa*dsa)/N;
	
	% and obtain the affine scalingx
	sigma = (mu_aff/mu)^3;


	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%3. Calculate the search direction
	% KKT system 
	%rs = [-rc;-rb;-xk.*sk-dxa.*dsa + sigma*mu];
	%ds = K\rs;
	%dxs = ds(1:N);
	%dss = ds(N+M+1:end);
	%dls = ds(N+1:N+M);
	
	%Solved by
	%dls = Kp\(-rb-A*diag(xk./sk)*rc+A*diag(1./sk)*(xk.*sk+dxa.*dsa - sigma*mu));
	rs = (-rb-Af((xkskm1).*rc-xk-skm1.*(dxa.*dsa - sigma*mu),Y,m,n));
	%v = C'\rs;
	%dls = C\v;

	v = linsolve(C,rs,opts1);
	dls = linsolve(C,v,opts2);

	%dss = -rc-A'*dls;
	Atdls = Atf(dls,Y);
	dss = -rc-Atdls;
	dxs = -xk-skm1.*(dxa.*dsa - sigma*mu)-(xkskm1).*dss;
	
	etak = one-0.1/k.^3; %0.99
	
	% Calculate the primal-dual stepsize and take the step
	I = find(dxs<0);
	ap = min([1.0;-etak*xk(I)./dxs(I)]);
    
	I = find(dss<0);
	ad = min([1.0;-etak*sk(I)./dss(I)]);

	xk = xk+ap*dxs;
	lk = lk+ad*dls;
	sk = sk+ad*dss;

end

% Obtain the original variables
a = (-xk(2*m+1:2*m+n)+xk(2*m+n+1:2*m+2*n))*0.5;
opt.k = k;
opt.status='solved';


function s = Atf(z,Y)
p = Y'*z;
s = [-z;z;-p;p];

function y = Af(x,Y,m,n)
 y =-x(1:m)+x(m+1:2*m) + Y*(-x(2*m+1:2*m+n)+x(2*m+n+1:end));
