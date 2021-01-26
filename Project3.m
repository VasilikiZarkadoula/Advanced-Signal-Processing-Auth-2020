% Giannakis Formula
clc;
clear;
close all;

% load real discrete signal x[k]
load('SignalX.mat');

% skewness of white non-Gaussian noise
sk = skewness(v);

% 3rd order cumulants of x[k] using the indirect method 
K = 32;
M = 64;
L3 = 20;
p = reshape(x,M,K);

[C3,~,c3,~] = bisp3cum(p,M,L3,'n','u');

% Graphic display
maxlag = L3;
samprate = M;
lagindex=-maxlag:maxlag;
lag=lagindex/samprate;
figure;
imagesc(lag,lag,c3);
axis xy
grid
%3d
figure;
axX=-20:20;
axY=-20:20;
surf(axX,axY,c3);

% estimate the impulse response of the MA system using Giannakis' formula
q = 5;
hEst = GiannnakisFormula(q,c3);

% estimate the impulse response of the MA system using Giannakis' formula
% a) consider sub-estimation of the order q
qSub = q-2;
hSub  =GiannnakisFormula(qSub,c3);
% b) consider sup-estimation of the order q
qSup = q+3;
hSup = GiannnakisFormula(qSup,c3);

% estimate the MA-q system output and compute the nrmse
[nrmse,xEst] = myFun(hEst,v,N,x);

% estimate the MA-q system output and compute the nrmse for the cases of
% hSub and hSup
[nrmseSub,xSubEst] = myFun(hSub,v,N,x);
[nrmseSup,xsup_est]=myFun(hSup,v,N,x);

NRMSE = [nrmse nrmseSub nrmseSup];

% repeat the above but instead of x[k] use the noise contaminated output
% yi[k]
snr = (30:-5:-5);

y = zeros(length(snr),N);
hEstY = zeros(length(snr),q+1);
nrmseY = zeros(1,length(snr));
for i = 1:length(snr)
   y(i,:) = awgn(x,snr(i),'measured');
   figure;
   p = reshape(y(i,:),M,K);
   [~,~,c3i,~] = bisp3cum(p,M,L3,'n','u');
   hEstY(i,:) = GiannnakisFormula(q,c3i);
   set(0,'DefaultFigureVisible','off');
   [nrmseY(i),yEsti] = myFun(hEstY(i,:),v,N,y(i,:));
end

set(0,'DefaultFigureVisible','on')
figure;
plot(snr,nrmseY);
title('NRMSE of y versus SNR range')
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

figure;
plot(x,'blue');
hold on
plot(xEst,'red');
legend('Origin','Estimated')
end




