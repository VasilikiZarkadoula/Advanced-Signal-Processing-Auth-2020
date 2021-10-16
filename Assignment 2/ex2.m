% Advanced Signal Processing
% Vasiliki Zarkadoula
clc
clear
close all

% Repeat exercise 1 but now take into account 50 realizations of X and 
% compare the mean values of the estimated power spectrums and bispectrums

% data
N = 8192;   % data length
lamda = [0.12 0.3 0.42 0.19 0.17 0.36];
omega = 2*pi*lamda;
M = 512;
L2 = 127;
G = 50;

% construct 50 realizations of X
m1 = zeros(1,G);
m2 = zeros(L2+1,G);
powerSpectrum = zeros(L2+1,G);
bispectrum = zeros(M,M,G);
X = zeros(N,1);
A = [];
for i=1:G
    rng shuffle
    a = 0;
    b = 2*pi;
    phi(1) = (b-a).*rand+a;
    phi(2) = (b-a).*rand+a;
    phi(3) = phi(1)+phi(2);
    phi(4) = (b-a).*rand+a;
    phi(5) = (b-a).*rand+a;
    phi(6) = phi(4)+phi(5);
    for k=1:N
        X(k)=0;
        for j=1:6
            X(k)=X(k)+cos(omega(j)*k+phi(j));
        end
    end
    m1(1,i) = mean(X);
    m2(:,i)=ACF(X,L2+1);
    powerSpectrum(:,i) = fft(m2(:,i)-m1(1,i)^2);
    bispectrum(:,:,i) = bispecd (X,M,1,M,0);
    A = [A X];
    
end

% mean +- std of Power Spectrum 
powerSpectrum = fftshift(powerSpectrum);
powerSpectrumMean = mean(abs(powerSpectrum),2);
poweSpectrumStd = std(powerSpectrum,[],2);

% mean + std of Bispectrum
bispectrumMean = mean(abs(bispectrum),3);
bispectrumStd = std(bispectrum,[],3);

% graphic display 
figure;
plot(powerSpectrumMean,'red')
hold on
plot((powerSpectrumMean+poweSpectrumStd),'blue')
hold on
plot((powerSpectrumMean-poweSpectrumStd),'green')
title('Power Spectrum Mean Values+-std')
xlabel('Samples')

fshift = (-L2/2:L2/2)/L2;
figure;
plot(fshift,powerSpectrumMean,'red')
hold on
plot(fshift,(powerSpectrumMean+poweSpectrumStd),'blue')
hold on
plot(fshift,(powerSpectrumMean-poweSpectrumStd),'green')
title('Power Spectrum Mean Values+-std')
xlabel('f[HZ]')

f1 = [-(M/2-1):M/2]/(M+1);
f2 = f1;
figure;
contour(f1,f2,(bispectrumMean),4); grid on
hold on;
contour(f1,f2,bispectrumMean+bispectrumStd,4); grid on
hold on;
contour(f1,f2,bispectrumMean-bispectrumStd,4);
title('Bispectrum via Direct Method  -  Mean Values +- Std')
xlabel('f1[HZ]'), ylabel('f2[HZ]')
