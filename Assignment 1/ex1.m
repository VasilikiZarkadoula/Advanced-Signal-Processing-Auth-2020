% Advanced Signal Processing
% Vasiliki Zarkadoula
clc
clear
close all

% data
N = 8192;   % data length
lamda = [0.12 0.3 0.42 0.19 0.17 0.36];
omega = 2*pi*lamda;
rng shuffle
a = 0;
b = 2*pi;
phi(1) = (b-a).*rand+a;
phi(2) = (b-a).*rand+a;
phi(3) = phi(1)+phi(2);
phi(4) = (b-a).*rand+a;
phi(5) = (b-a).*rand+a;
phi(6) = phi(4)+phi(5);

%--------------------------------------------------------------------------

% 1.

% Construct real discrete process X

X=zeros(N,1);
for k=0:N-1
    X(k+1)=0;
    for j=1:6
        X(k+1)=X(k+1)+cos(omega(j)*k+phi(j));
    end
   
end
plot(X); 
title('Real discrete process X')

%--------------------------------------------------------------------------

% 2.

% Estimate the power spectrum, use 128 shiftings for autocorrelation

m1 = mean(X);       
m2 = ACF(X,128);    
c2 = m2-m1^2;       %covariance
C2 = fft(c2);       %power spectrum
C2 = fftshift(C2);
fs = 1;
n = length(C2);
x = (-(n-1)/2:(n-1)/2)*(fs/n);
y = abs(C2).^2/n;
figure;
plot(x, y);         
xlabel('Frequency')
ylabel('Power')
title('Power Spectrum')

%--------------------------------------------------------------------------
       
% 3.

% Estimate the bispectrum using indirect and direct method

M = 256;
K = 32;
L = 64;

% reshape X[k] from(N*1)vector to (M*K) array
Y = reshape(X,M,K);

% a)indirect method
figure;
C3a1 = bisp3cum(Y,M,L,'n','u');   %rectangular window

figure;
C3a2 = bisp3cum(Y,M,L,'pa','u');  %parzen window

% b)direct method
figure;
C3b = bispecd (X,M,1,M,0);
maxlag = L;
samprate = M;
lagindex = -maxlag:maxlag;
freq = lagindex/maxlag/2*samprate;
imagesc(freq,freq,abs(C3b));
axis xy
grid
colormap gray
