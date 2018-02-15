clear; clc; close all;

% reading file

[SinX fs nbits]=wavread('../data/fa.wav');

x=SinX(:,1);
 
% parameters

L=80; D=88; M=256;
S=M/2+1;

N=floor(length(x)/(L+D))-5;

quadwin=0;

% sin-cos window
wu=sin(2*pi*((0:D-1)+0.5)./(D*4));
wm=ones(1,L);
wd=cos(2*pi*((0:D-1)+0.5)./(D*4));
w=[wu wm wd]';
if quadwin==1
    w=ones(L+D+D,1);
end
% matrices to store data

x_mat=zeros(M,N);
X_mat=zeros(M/2+1,N);
X_mod=zeros(M/2+1,N);
x_mod=zeros(M,N);


%% segmentation & windowing

for n=1:N
      
        t1=1+(L+D)*(n-1);
        t2=M+(L+D)*(n-1);
        
        x_mat(:,n)=x(t1:t2).*w;
        
               
    % FFT
    
        X=fft(x_mat(:,n),M);
      
        X_mat(:,n)=X(1:M/2+1);       
    
        
        
end


%% modification

X_mod(:,1)=X_mat(:,1);
X_mem=zeros(size(X_mat(:,1)));

for n=1:N

    
        X_mod(:,n)=X_mat(:,n);
                                
        % synthesis    
        XX=[X_mod(:,n); flipud(conj(X_mod(2:M/2,n)))];
        
       
        x_mod(:,n)=real(ifft(XX));

end

        


%% reconstruction
        
xrec=zeros(size(x));
xmufrec=zeros(size(x));

if quadwin==1
    wu=sin(2*pi*((0:D-1)+0.5)./(D*4));
    wm=ones(1,L);
    wd=cos(2*pi*((0:D-1)+0.5)./(D*4));
    w=[wu wm wd]'.^2;
end

for n=1:N
    

    
        t1=1+(L+D)*(n-1);
        t2=M+(L+D)*(n-1);
    
        xrec(t1:t2)=xrec(t1:t2)+x_mod(:,n).*w;
				%xmufrec(t1:t2)=xmufrec(t1:t2)+x_muf(:,n).*w;
     
    
end


figure; subplot(211); plot(x); subplot(212); plot(xrec,'r'); 

