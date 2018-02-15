function [a opt] = pdip_slp_primal(y,Y,gamma,epsilon,verbose,varargin)
% Solves
%
% minimize   ||y-Ya||_1 + gamma ||a||_1
%
% to epsilon accuracy
%
% Tobias Lindstroem Jensen
% tlj@es.aau.dk, Aalborg University, June 6, 2012

COUNT_MAX = 5;

[m,n]=size(Y);

%Form new matrices
%A = [-Y; gamma*eye(n)];
b = [-y;zeros(n,1)];

N = 2*n+2*m;

%G = [A -eye(M);-A -eye(M)];
%h = [b;-b];
%c = [zeros(N,1);ones(M,1)];

opts1.LT = true; 
opts2.LT = true; opts2.TRANSA = true;

cls = superiorfloat(y);
one = ones(1,1,cls);

if verbose
	fprintf('Ite.      pinfeas.     dinfeas.      rel_gap\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Obtain initial settings

% Simple initial setting
%x0 = A\b;

[U S V] = svd(Y);
x0 = V*inv(S'*S+diag(ones(size(S,2),1)*gamma^2))*S'*U'*y;

rls = Af(x0,Y,gamma) - b;

y0 = abs(rls)+1e-1;

%xk = [x0;y0];
z = (1-1e-1)*rls/norm(rls,'inf');

%zk = [(1+z)/2;(1-z)/2];
%sk = ones(length(h),1)*1e-2;

zk1 = (1+z)/2;
zk2 = (1-z)/2;
sk1 = ones(length(zk1),1)*1e-2;
sk2 = ones(length(zk2),1)*1e-2;

uk = x0;
vk = y0;

for k=1:100
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%1. compute residual and evaluate stopping criteria
	%uk = xk(1:n);
	%vk = xk(n+1:end);
	%sk1 = sk(1:M);
	%sk2 = sk(M+1:end);
	%zk1 = zk(1:M);
	%zk2 = zk(M+1:end);
	
	%rc = -G'*zk - c;
	%rh = sk + G*xk - h;
	

	rc1 = -Atf(zk1-zk2,Y,gamma,m);
	rc2 = zk1 + zk2 - 1;
	
	if k ==  1
		g = Af(uk,Y,gamma); %First iteration calc via a function handle
	else
		g = g+ap*gs; %Compute via precomputed values
	end

	rh1 = sk1 + g - vk - b;
	rh2 = sk2 - g - vk + b;

	%mu = sk'* zk/N;
	mu = (sk1'*zk1+sk2'*zk2)/N;
	
	%pinfeas = norm(rh)/(one+norm(h))
	%dinfeas = norm(rc)/(one+norm(c))
	%relgap = (c'*xk+h'*zk)/(one+abs(c'*xk))
	pinfeas = norm([rh1;rh2])/(one+sqrt(2)*norm(b));
	dinfeas = norm([rc1;rc2])/(one+sqrt(2*(m+n)));
	relgap = (sum(vk)+b'*(zk1-zk2))/(one+abs(sum(vk)));

	if verbose
		%fprintf('%2d       %.2e      %.2e ...
				%.2e\n',k,pinfeas,dinfeas,relgap);
		fprintf('%2d       %.2e      %.2e      %.2e\n',k,pinfeas,dinfeas,relgap);
	end


	if max([relgap,pinfeas]) <= epsilon
		break;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%2. Calculate the affine scaling
	%rs = sk.*zk;
	
	rs1 = sk1.*zk1;
	rs2 = sk2.*zk2;
	
	d1 = zk1./sk1;
	d2 = zk2./sk2;
	
	% KKT system
	%K = G'*diag(zk./sk)*G;
	%rhs = (rc-G'*diag(zk./sk)*rh+G'*(rs./sk));
	%dxa = K\rhs;
	%dxa = [du;dv];
	
	%r1 = rc1 - Atf((zk1.*rh1-rs1)./sk1 - (zk2.*rh2-rs2)./sk2,m,Y,gamma);
	%r2 = rc2 + ((zk1.*rh1-rs1)./sk1 + (zk2.*rh2-rs2)./sk2);

	p1 = (zk1.*rh1-rs1)./sk1;
	p2 = (zk2.*rh2-rs2)./sk2;
	%r1 = - Atf(p1-p2,Y,gamma,m);
	r2 = rc2+p1+p2;
	
	%D = diag(4*(d1.*d2./(d1+d2)));
	%K = (A'*D*A)
	d = 4*(d1.*d2./(d1+d2));
	
	%Kp = bsxfun(@times,d,A);
	%K = A'*Kp;

	Kp = bsxfun(@times,d(1:m),Y);
	
	K = Y'*Kp; % + diag(gamma^2*d(m+1:end));

	K(1:n+1:end) = K(1:n+1:end)'+gamma^2*d(m+1:end);

	[C,p] = chol(K,'lower');

	count = 0;
	while p ~= 0
		if verbose
			fprintf('Cholesky failed in iteration %d. Add regularization\n',k);
		end

		K(1:n+1:end) = K(1:n+1:end)'+epsilon; %adds regularization to the diagonal
		[C,p] = chol(K,'lower');
			
		count = count + 1;
		
		if count > COUNT_MAX
			a = uk;
			opt.k = k;
			opt.status = 'failed';
			return
		end
		
	end

	%du = K\(r1-A'*(((d2-d1)./(d1+d2)).*r2));
	rhs = rc1-Atf((p1-p2)+((d2-d1)./(d1+d2)).*r2,Y,gamma,m);
	
	v = linsolve(C,rhs,opts1);
	du = linsolve(C,v,opts2);

	gs = Af(du,Y,gamma);
	dv = -(d2-d1).*gs./(d1+d2)+r2./(d1+d2);

	%dza = -diag(zk./sk)*(-rh+rs./zk-G*dxa);
	dza1 = -(zk1./sk1).*(-rh1+rs1./zk1-gs+dv);
	dza2 = -(zk2./sk2).*(-rh2+rs2./zk2+gs+dv);
	
	%dsa = (-rs-sk.*dza)./zk;
	dsa1 = (-rs1-sk1.*dza1)./zk1;
	dsa2 = (-rs2-sk2.*dza2)./zk2;
	
	% Calculate the primal-dual stepsize
	%I = find(dsa<0);
	%asa = min([1.0;-sk(I)./dsa(I)]);
  
	%I = find(dza<0);
	%aza = min([1.0;-zk(I)./dza(I)]);


	I1 = find(dsa1<0);
	I2 = find(dsa2<0);
	asa = min([1.0;-sk1(I1)./dsa1(I1);-sk2(I2)./dsa2(I2)]);
  
	I1 = find(dza1<0);
	I2 = find(dza2<0);
	aza = min([1.0;-zk1(I1)./dza1(I1);-zk2(I2)./dza2(I2)]);

	%mu_aff = (sk+asa*dsa)'*(zk+aza*dza)/N;
	mu_aff = ((sk1+asa*dsa1)'*(zk1+aza*dza1)+(sk2+asa*dsa2)'*(zk2+aza*dza2))/N;

	% and obtain the affine scalingx
	sigma = (mu_aff/mu)^3;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%3. Calculate the search direction
	%rs = sk.*zk+dsa.*dza-sigma*mu;

	rs1 = sk1.*zk1+dsa1.*dza1-sigma*mu;
	rs2 = sk2.*zk2+dsa2.*dza2-sigma*mu;
	
	%r1 = rc1 - Atf( (zk1.*rh1-rs1)./sk1 - (zk2.*rh2-rs2)./sk2,m,Y,gamma);
	%r2 = rc2 + ((zk1.*rh1-rs1)./sk1 + (zk2.*rh2-rs2)./sk2);

	%r1 = - Atf( (zk1.*rh1-rs1)./sk1 - (zk2.*rh2-rs2)./sk2,Y,gamma,m);
	r2 = rc2 + ((zk1.*rh1-rs1)./sk1 + (zk2.*rh2-rs2)./sk2);

	%dxs = K\(rc-G'*diag(zk./sk)*rh+G'*(rs./sk));

	%dus = K\(r1-A'*(((d2-d1)./(d1+d2)).*r2));
	rhs = rc1 - Atf(((zk1.*rh1-rs1)./sk1 - (zk2.*rh2-rs2)./sk2)+ ((d2-d1)./(d1+d2)).*r2,Y,gamma,m);
	v = linsolve(C,rhs ,opts1);
	dus = linsolve(C,v,opts2);

	gs = Af(dus,Y,gamma);
	dvs = -(d2-d1).*gs./(d1+d2)+r2./(d1+d2);

	%dzs = -diag(zk./sk)*(-rh+rs./zk-G*dxs);
	dzs1 = -zk1./sk1.*(-rh1+rs1./zk1-gs+dvs);
	dzs2 = -zk2./sk2.*(-rh2+rs2./zk2+gs+dvs);
	
	%dss = (-rs-sk.*dzs)./zk;
	dss1 = (-rs1-sk1.*dzs1)./zk1;
	dss2 = (-rs2-sk2.*dzs2)./zk2;

	etak = one-0.1/k.^3; %0.99
	
	% Calculate the primal-dual stepsize and take the step
	%I = find(dss<0);
	%ap = min([1.0;-etak*sk(I)./dss(I)]); 
	%I = find(dzs<0);
	%ad = min([1.0;-etak*zk(I)./dzs(I)]);

	I1 = find(dss1<0);
	I2 = find(dss2<0);
	ap = min([1.0;-etak*sk1(I1)./dss1(I1);-etak*sk2(I2)./dss2(I2)]);
    
	I1 = find(dzs1<0);
	I2 = find(dzs2<0);
	ad = min([1.0;-etak*zk1(I1)./dzs1(I1);-etak*zk2(I2)./dzs2(I2)]);

	%xk = xk+ap*dxs;
	uk = uk+ap*dus;
	vk = vk+ap*dvs;
	
	%sk = sk+ap*dss;
	sk1 = sk1+ap*dss1;
	sk2 = sk2+ap*dss2;
	
	%zk = zk+ad*dzs;
	zk1 = zk1+ad*dzs1;
	zk2 = zk2+ad*dzs2;

end

% Obtain the original variables
a = uk;
opt.k = k;
opt.status = 'solved';

function s = Atf(z,Y,gamma,m)
s = -Y'*z(1:m)+gamma*z(m+1:end);


function y = Af(z,Y,gamma)
y = [-Y*z;gamma*z];

