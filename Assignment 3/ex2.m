clc;
clear;
close all;

% Repeat Project 3 for 50 realizations of input and output and work with
% the mean values

%data length
N = 2048;
%variables for cumulants
K = 32;
M = 64;
L3 = 20;
%order of impulse response
q = 5;
qSub = q-2;
qSup = q+3;
%confidence interval
vcrit = 1.96;
lconfbound = -vcrit/sqrt(N);
upconfbound = vcrit/sqrt(N);

snr = (30:-5:-5);
n = 50;

hEst = zeros(n,q+1);
hSub = zeros(n,qSub+1);
hSup = zeros(n,qSup+1);
nrmseEst = zeros(1,n);
nrmseEstSub = zeros(1,n);
nrmseEstSup = zeros(1,n);
nrmseEstY = zeros(n,length(snr));
for j=1:n
    % construct 50 inputs kai outputs
    [X,V] = Signals_xv();
    
    % skewness of V
    SK = skewness(V');
    
    % cumulants
    p = reshape(X,M,K);
    [~,~,cumX,~] = bisp3cum(p,M,L3,'n','u'); 
    close
    
    hEst(j,:) = GiannnakisFormula(q,cumX);
    hSub(j,:) = GiannnakisFormula(qSub,cumX);
    hSup(j,:) = GiannnakisFormula(qSup,cumX);
    
    [nrmseEst(j),xEst] = myFun(hEst(j,:),v,N,X);
    [nrmseEstSub(j),xEstSub] = myFun(hSub(j,:),v,N,X);
    [nrmseEstSup(j),xEstSup] = myFun(hSup(j,:),v,N,X);
  
    for i=1:8
        y = awgn(X,snr(i),'measured');
        r = reshape(y,M,K);
        [~,~,cumY,~] = bisp3cum(r,M,L3,'n','u');
        close
        hEstY = GiannnakisFormula(q,cumY);
        [nrmseEstY(j,i),yEst] = myFun(hEstY,v,N,y);
        
    end
end

skMean = median(SK);

hEstMean = median(hEst);
hSubMean = median(hSub);
HSupMean = median(hSup);

nrmseXmean = median(nrmseEst);
nrmseXsubMean = median(nrmseEstSub);
nrmseXsupMean = median(nrmseEstSup);
nrmseX = [nrmseXmean; nrmseXsubMean; nrmseXsupMean]';

figure;
boxplot(nrmseX,{'q','qsub','qsup'})
title('Compare NRSME for different order q')

% nrsme of y versus snr
nrmseYmean = median(nrmseEstY);
nrmseYstd = std(nrmseEstY);
% 95% confidence interval
S = size(nrmseEstY,1);
nrsmeYs = nrmseYstd/sqrt(S);
ConfInt_up95 = bsxfun(@times,nrsmeYs,upconfbound);
ConfInt_low95 = bsxfun(@times,nrsmeYs,lconfbound);

%plot nrmse+-std vs snr
figure
plot(snr,nrmseYmean);
hold on
plot(snr,nrmseYmean+nrmseYstd,'r--');
hold on
plot(snr,nrmseYmean-nrmseYstd,'k--');
title('Mean NRMSE of Y +-std versus SNR range')
xlabel('SNR')
ylabel('NRMSE')
hold off

%plot nrmse+confidence interval vs snr
figure
plot(snr,nrmseYmean);
hold on
plot(snr,nrmseYmean+ConfInt_up95,'r--');
hold on
plot(snr,nrmseYmean+ConfInt_low95,'k--');
title('Mean NRMSE of Y +confidence interval versus SNR range')
xlabel('SNR')
ylabel('NRMSE')

function h = GiannnakisFormula(q,c3)
h = NaN(1,length(q)+1);
for k=0:q
    h(k+1) = c3(k+21,q+21)/c3(21,q+21);
end
end

function [nrmse,xEst] = myFun(h,v,N,x)
xEst = conv(h,v);
xEst = xEst(1:N);
dif = 0;
for k=1:N
    dif = dif + (xEst(k)-x(k))^2;
end
rmse = sqrt(dif/N);
nrmse = rmse/(max(x)-min(x));

end



