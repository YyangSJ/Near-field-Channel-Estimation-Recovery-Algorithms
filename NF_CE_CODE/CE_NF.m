clear all;
close all;
warning off;
rng(2024) 

f=10e9; % system frequency
lambda=3e8/f; % wavelength
d=lambda/2; % element spacing
N=100; % number of antennas
Q=100; % number of measurements
P=8; % number of channel paths
r=3:1:30;  % distance samples in [3,30]
theta=-pi/2:0.1:pi/2; %angle samples in [-pi/2,pi/2]
realization=200;
SNR=[0,10,20];
for reali=1:realization
    alpha=1/sqrt(2)*(normrnd(0,1,P,1)+1i*normrnd(0,1,P,1));
    angle_rnd=unifrnd(-pi/2,pi/2,P,1); 
    r_rnd=unifrnd(3,30,P,1); 
    % angle_rnd=theta(randperm(length(theta),P)); % on-grid parameter
    % r_rnd=r(randperm(length(r),P)); % on-grid parameter
    G_size=length(theta)*length(r);
    G=zeros(N,G_size);
    fla=0;
    for ang=1:length(theta)
        for rr=1:length(r)
            fla=fla+1;
            G(:,fla)=AM(N,lambda,d,r(rr),theta(ang));  % dictionary G generation
        end
    end
    h=zeros(N,1);
    for p=1:P
        h=h+ alpha(p)*AM(N,lambda,d,r_rnd(p),angle_rnd(p));  % channel generation
    end
    % W=normrnd(0,1,N,Q)+1i*normrnd(0,1,N,Q);
    % W=W./abs(W);
    % W = randi([0, 1], N, Q);
    % W(W == 0) = -1;
    W=eye(N,Q);  % measurement matrix, a simple case is considered.
    n=1/sqrt(2)*(normrnd(0,1,N,1)+1i*normrnd(0,1,N,1)); % AWGN
    for snr=1:length(SNR)
        pow=10^(SNR(snr)/10);
        n=sqrt(1/pow)*n;
        y=W'*h+W'*n;
        
        %% sparse recovery algorithms, convex/Bayesian methods are setting to maximum 200 iterations
        Iter_max=200;
        [~,x_omp]=cs_somp(y,W'*G,G_size,P); % OMP
        [supp_ols,x_ols]=OLS_MMV(y,W'*G,P); % OLS
        x_FISTA=FISTA(W'*G,y,zeros(G_size,1),1.05,1.01,0.15,Iter_max,1, 0); % FISTA
        x_ADMM=ADMM(W'*G,y,0.7,0.15,Iter_max); % ADMM
        x_mfvsbl=MFV_SBL(y,W'*G,(W'*G)'*W'*G,Iter_max); % SAVE-SBL
        x_fmfsbl=FMFSBL(W'*G,y,Iter_max); %FMFSBL
        x_vsbl=V_SBL(y,W'*G,Iter_max); % VSBL
        x_sbl=SBL(y,W'*G,Iter_max);  % SBL
        x_laomp=LAOMP_f(y,W'*G,P,P); % LAOMP
        
        
        nmse_omp(reali,snr)=norm(G*x_omp-h)^2/norm(h)^2
        nmse_ols(reali,snr)=norm(G(:,supp_ols)*x_ols-h)^2/norm(h)^2
        nmse_FISTA(reali,snr)=norm(G*x_FISTA(:,Iter_max)-h)^2/norm(h)^2
        nmse_ADMM(reali,snr)=norm(G*x_ADMM(:,Iter_max)-h)^2/norm(h)^2
        nmse_mfvsbl(reali,snr)=norm(G*x_mfvsbl-h)^2/norm(h)^2
        nmse_fmfsbl(reali,snr)=norm(G*x_fmfsbl-h)^2/norm(h)^2
        nmse_vsbl(reali,snr)=norm(G*x_vsbl-h)^2/norm(h)^2
        nmse_sbl(reali,snr)=norm(G*x_sbl-h)^2/norm(h)^2
        nmse_laomp(reali,snr)=norm(G*x_laomp-h)^2/norm(h)^2
    end
    
    mean_nmse_omp=mean(nmse_omp,1);
    mean_nmse_ols=mean(nmse_ols,1);
    mean_nmse_FISTA=mean(nmse_FISTA,1);
    mean_nmse_ADMM=mean(nmse_ADMM,1);
    mean_nmse_mfvsbl=mean(nmse_mfvsbl,1);
    mean_nmse_fmfsbl=mean(nmse_fmfsbl,1);
    mean_nmse_vsbl=mean(nmse_vsbl,1);
    mean_nmse_sbl=mean(nmse_sbl,1);
    mean_nmse_laomp=mean(nmse_laomp,1);
    
    
end

%% spherical-wave array manifold
function g= AM(N,lambda,d,r,theta)
g=zeros(N,1);
for n=1:N
    r_n=sqrt(r^2-2*r*(n-1)*d*sin(theta)+(n-1)^2*d^2);
    g(n)=r/r_n*exp(-1i*2*pi/lambda*(r_n-r));
end
end


