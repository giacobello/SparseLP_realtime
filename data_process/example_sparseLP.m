clear; close all; clc;
% this script is used to compare the various method for linear prediction:
% for each method (AC2,CV2,AC1,CV1,R11(gamma=1),R101(gamma=0.1),R21(gamma=1),R201(gamma=0.1),B1,B2,HU,HUR(gamma=1)) the output is:
% 1. coefficient vectors
% 2. stability check
% 3. residual
% results are saved in: 
% mixO40F320.mat
% faO40F320.mat

%% file reading

[y fs1 bq]=wavread('../data/fa.wav');
fs=8000; %new sampling frequency
if fs1~=fs
    xo=resample(y,fs,fs1)'; %resampling from fs1 to 8000Hz 
    else xo=y;
end
frame=320; %%frame length (multiple of 4)
nframes=floor(length(xo)/frame);
xo=xo(1:nframes*frame); %resizing input speech

%% removal of low frequencies components

fc=50; %cut-off frequency (Hz)
ord=2; %order
rp=20; %ripple amplitude
[b,a] = cheby2(ord,rp,fc/fs,'high');
x=filter(b,a,xo);

%% LPC analysis (autocorrelation method)
order=40; %order of the LPC filter
tetaAC2=zeros(order+1,nframes);
stab_checkAC2=zeros(1,nframes);
resAC2=zeros(size(x));
for k=1:nframes
    [r,lags]=xcorr(x((1+(k-1)*frame):k*frame)); %autocorrelation (biased)
    acm=zeros(order+1); %creation of autoccrealtion matrix (can be replaced by toeplitz(r(frame:frame+order)) )
    for i=0:order
        acm=acm+diag(ones(1,order+1-i)*r(frame+i),i);
    end
    for i=0:order-1
        acm=acm+diag(ones(1,order-i)*r(frame+i+1),-i-1);
    end
    rc=r(frame+1:frame+order);
    acmc=acm(1:order,1:order);
    tetaAC2(:,k)=[1; -inv(acmc)*rc]; %calculation of A(z) parameters
    stab_checkAC2(k)=sum(abs(roots(tetaAC2(:,k)'))'>1); %stability check --> poles inside unit circle
    resAC2((1+(k-1)*frame):k*frame)=filter(tetaAC2(:,k)',1,x((1+(k-1)*frame):k*frame));
end

%% 1-norm minimization regularized (N1=1 and N2=N+K), gamma=1

%order=40; %order of the LPC filter
tetaR11=zeros(order+1,nframes);
stab_checkR11=zeros(1,nframes);
resR11=zeros(size(x));
for k=1:nframes    
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    gamma=1;
		keyboard
		cvx_begin
        variable a1(order)
        minimize norm(Y*a1+y,1)+gamma*norm(a1,1)
    cvx_end
    tetaR11(:,k)=[1;a1];
    stab_checkR11(k)=sum(abs(roots(tetaR11(:,k)'))'>1);    
    resR11((1+(k-1)*frame):k*frame)=filter(tetaR11(:,k)',1,x((1+(k-1)*frame):k*frame));
end


%% LPC analysis (covariance method)
%order=40; %order of the LPC filter
tetaCV2=zeros(order+1,nframes);
stab_checkCV2=zeros(1,nframes);
resCV2=zeros(size(x));
for k=1:nframes
    s=x((1+(k-1)*frame):k*frame);
    y=s(order+1:frame);
    Y=toeplitz(s(order:frame-1),s(order:-1:1)');
    tetaCV2(:,k)=[1; -(Y'*Y)\(Y'*y)];
    stab_checkCV2(k)=sum(abs(roots(tetaCV2(:,k)'))'>1);
    resCV2((1+(k-1)*frame):k*frame)=filter(tetaCV2(:,k)',1,x((1+(k-1)*frame):k*frame));
end


%% 1-norm minimization with N1=1 and N2=N+K (same as autocorrelation)

%order=40; %order of the LPC filter
tetaAC1=zeros(order+1,nframes);
stab_checkAC1=zeros(1,nframes);
resAC1=zeros(size(x));
for k=1:nframes    
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    cvx_begin
        variable a1(order)
        minimize norm(Y*a1+y,1)
    cvx_end
    tetaAC1(:,k)=[1;a1];
    stab_checkAC1(k)=sum(abs(roots(tetaAC1(:,k)'))'>1);    
    resAC1((1+(k-1)*frame):k*frame)=filter(tetaAC1(:,k)',1,x((1+(k-1)*frame):k*frame));
end


%% 1-norm minimization with N1=K+1 and N2=N (same as covariance)

%order=40; %order of the LPC filter
tetaCV1=zeros(order+1,nframes);
stab_checkCV1=zeros(1,nframes);
resCV1=zeros(size(x));
for k=1:nframes
    s=x((1+(k-1)*frame):k*frame);
    y=s(order+1:frame);
    Y=toeplitz(s(order:frame-1),s(order:-1:1)');
    cvx_begin
        variable a1(order)
        minimize norm(Y*a1+y,1)
    cvx_end
    tetaCV1(:,k)=[1;a1]';
    stab_checkCV1(k)=sum(abs(roots(tetaCV1(:,k)'))'>1);    
    resCV1((1+(k-1)*frame):k*frame)=filter(tetaCV1(:,k)',1,x((1+(k-1)*frame):k*frame));
end


%% 1-norm minimization regularized (N1=1 and N2=N+K), gamma=0.1

%order=40; %order of the LPC filter
tetaR101=zeros(order+1,nframes);
stab_checkR101=zeros(1,nframes);
resR101=zeros(size(x));
for k=1:nframes    
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    gamma=0.1;
    cvx_begin
        variable a1(order)
        minimize norm(Y*a1+y,1)+gamma*norm(a1,1)
    cvx_end
    tetaR101(:,k)=[1;a1];
    stab_checkR101(k)=sum(abs(roots(tetaR101(:,k)'))'>1);    
    resR101((1+(k-1)*frame):k*frame)=filter(tetaR101(:,k)',1,x((1+(k-1)*frame):k*frame));
end

%% 2-norm minimization regularized (N1=1 and N2=N+K), gamma=1

%order=40; %order of the LPC filter
tetaR21=zeros(order+1,nframes);
stab_checkR21=zeros(1,nframes);
resR21=zeros(size(x));
for k=1:nframes    
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    gamma=1;
    cvx_begin
        variable a1(order)
        minimize norm(Y*a1+y,2)+gamma*norm(a1,1)
    cvx_end
    tetaR21(:,k)=[1;a1];
    stab_checkR21(k)=sum(abs(roots(tetaR21(:,k)'))'>1);    
    resR21((1+(k-1)*frame):k*frame)=filter(tetaR21(:,k)',1,x((1+(k-1)*frame):k*frame));
end

%% 2-norm minimization regularized (N1=1 and N2=N+K), gamma=0.1

%order=40; %order of the LPC filter
tetaR201=zeros(order+1,nframes);
stab_checkR201=zeros(1,nframes);
resR201=zeros(size(x));
for k=1:nframes    
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    gamma=0.1;
    cvx_begin
        variable a1(order)
        minimize norm(Y*a1+y,2)+gamma*norm(a1,1)
    cvx_end
    tetaR201(:,k)=[1;a1];
    stab_checkR201(k)=sum(abs(roots(tetaR201(:,k)'))'>1);    
    resR201((1+(k-1)*frame):k*frame)=filter(tetaR201(:,k)',1,x((1+(k-1)*frame):k*frame));
end


%% burg method based on the 1-norm ans 2-norm
%order=40;
tetaB1=zeros(order+1,nframes);
stab_checkB1=zeros(1,nframes);
resB1=zeros(size(x));
tetaB2=zeros(order+1,nframes);
stab_checkB2=zeros(1,nframes);
resB2=zeros(size(x));
for k=1:nframes
    tetaB1(:,k)=arburgABS(x((1+(k-1)*frame):k*frame),order)';
    stab_checkB1(k)=sum(abs(roots(tetaB1(:,k)'))'>1);
    resB1((1+(k-1)*frame):k*frame)=filter(tetaB1(:,k)',1,x((1+(k-1)*frame):k*frame));
    tetaB2(:,k)=arburg(x((1+(k-1)*frame):k*frame),order)';
    stab_checkB2(k)=sum(abs(roots(tetaB2(:,k)'))'>1);
    resB2((1+(k-1)*frame):k*frame)=filter(tetaB2(:,k)',1,x((1+(k-1)*frame):k*frame));
end

%% minimization with HUBER cost function
%order=40; %order of the LPC filter
tetaHU=zeros(order+1,nframes);
stab_checkHU=zeros(1,nframes);
resHU=zeros(size(x));
for k=1:nframes
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    cvx_begin
        variable a1(order)
        minimize sum(huber(Y*a1+y))
    cvx_end
    tetaHU(:,k)=[1;a1]';
    stab_checkHU(k)=sum(abs(roots(tetaHU(:,k)'))'>1);
    resHU((1+(k-1)*frame):k*frame)=filter(tetaHU(:,k)',1,x((1+(k-1)*frame):k*frame));
end

%% minimization with HUBER cost function regularized
%order=40; %order of the LPC filter
tetaHUR=zeros(order+1,nframes);
stab_checkHUR=zeros(1,nframes);
resHUR=zeros(size(x));
for k=1:nframes
    s=[0 x((1+(k-1)*frame):k*frame)' zeros(1,order-1)];
    y=[s(2:end) 0]';
    Y=toeplitz(s,s(order:-1:1)');
    Y(1:order,1:order)=tril(Y(1:order,1:order));
    gamma=1;
    cvx_begin
        variable a1(order)
        minimize sum(huber(Y*a1+y))+gamma*norm(a1,1)
    cvx_end
    tetaHUR(:,k)=[1;a1]';
    stab_checkHUR(k)=sum(abs(roots(tetaHUR(:,k)'))'>1);
    resHUR((1+(k-1)*frame):k*frame)=filter(tetaHUR(:,k)',1,x((1+(k-1)*frame):k*frame));
end

%% figures

figure; 
plot(resAC2); hold on; plot(resCV2,'g'); hold on; plot(resAC1,'r'); hold on; plot(resCV1,'m'); hold on; plot(resB1,'k'); hold on; plot(resB2,'c'); hold on; plot(resHU,'y');
legend('resAC2','resCV2','resAC1','resCV1','resB1','resB2','resHU');
grid on;

figure; 
plot(resR11); hold on; plot(resR101,'g'); hold on; plot(resR21,'m'); hold on; plot(resR201,'k'); hold on; plot(resHUR,'r')
legend('resR11','resR101','resR21','resR201','resHUR');
grid on;









