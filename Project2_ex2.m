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


%BISPECD function - HOSA toolbox
function [Bspec,waxis] = bispecd (y,  nfft, wind, nsamp, overlap)
[ly, nrecs] = size(y);
    if (ly == 1) y = y(:);  ly = nrecs; nrecs = 1; end

    if (exist('nfft') ~= 1)            nfft = 128; end
    if (exist('overlap') ~= 1)      overlap = 50;  end
    overlap = min(99,max(overlap,0));
    if (nrecs > 1)                  overlap =  0;  end
    if (exist('nsamp') ~= 1)          nsamp = 0;   end
    if (nrecs > 1)                    nsamp = ly;  end

    if (nrecs == 1 & nsamp <= 0)
       nsamp = fix(ly/ (8 - 7 * overlap/100));
    end
    if (nfft  < nsamp)   nfft = 2^nextpow2(nsamp); end

    overlap  = fix(nsamp * overlap / 100);             % added 2/14
    nadvance = nsamp - overlap;
    nrecs    = fix ( (ly*nrecs - overlap) / nadvance);


% ------------------- create the 2-D window -------------------------
  if (exist('wind') ~= 1) wind = 5; end
  [m,n] = size(wind);
  window = wind;
  if (max(m,n) == 1)     % scalar: wind is size of Rao-Gabr window
     winsize = wind;
     if (winsize < 0) winsize = 5; end        % the window length L
     winsize = winsize - rem(winsize,2) + 1;  % make it odd
     if (winsize > 1)
        mwind   = fix (nfft/winsize);            % the scale parameter M
        lby2    = (winsize - 1)/2;

        theta  = -lby2:lby2;
        opwind = ones(winsize,1) * (theta .^2);       % w(m,n)=m^2
        opwind = opwind + opwind' + theta' * theta;   % m^2 + n^2 + mn
        opwind = 1 - (2*mwind/nfft)^2 * opwind;       %
        hex    = ones(winsize,1) * theta;             % m
        hex    = abs(hex) + abs(hex') + abs(hex+hex');
        hex    = (hex < winsize);
        opwind = opwind .* hex;
        opwind = opwind * (4 * mwind^2) / (7 * pi^2) ;
     else
        opwind = 1;
     end

  elseif (min(m,n) == 1)  % 1-D window passed: convert to 2-D
     window = window(:);
     if (any(imag(window) ~= 0))
        disp(['1-D window has imaginary components: window ignored'])
        window = 1;
     end
     if (any(window < 0))
        disp(['1-D window has negative components: window ignored'])
        window = 1;
     end
     lwind  = length(window);
     windf  = [window(lwind:-1:2); window];    % the full symmetric 1-D
     window = [window; zeros(lwind-1,1)];
     opwind = (windf * windf')      ...
              .* hankel(flipud(window), window); % w(m)w(n)w(m+n)
     winsize = length(window);

  else                    % 2-D window passed: use directly
    winsize = m;
    if (m ~= n)
       disp('2-D window is not square: window ignored')
       window = 1;
       winsize = m;
    end
    if (rem(m,2) == 0)
       disp('2-D window does not have odd length: window ignored')
       window = 1;
       winsize = m;
    end
    opwind  = window;
  end

% ---------------- accumulate triple products ----------------------

    Bspec    = zeros(nfft,nfft);

    mask = hankel([1:nfft],[nfft,1:nfft-1] );   % the hankel mask (faster)
    locseg = [1:nsamp]';
    for krec = 1:nrecs
        xseg   = y(locseg);
        Xf     = fft(xseg-mean(xseg), nfft)/nsamp;
        CXf    = conj(Xf);
        Bspec  = Bspec + (Xf * Xf.') .* ...
	         reshape(CXf(mask), nfft, nfft);
        locseg = locseg + nadvance;
    end

    Bspec = fftshift(Bspec)/(nrecs);



% ----------------- frequency-domain smoothing ------------------------

  if (winsize > 1)
      lby2 = (winsize-1)/2;
      Bspec = conv2(Bspec,opwind);
      Bspec = Bspec(lby2+1:lby2+nfft,lby2+1:lby2+nfft);
  end
% ------------ contour plot of magnitude bispectum --------------------

   if (rem(nfft,2) == 0)
       waxis = [-nfft/2:(nfft/2-1)]'/nfft;
   else
       waxis = [-(nfft-1)/2:(nfft-1)/2]'/nfft;
   end

   hold off, clf
%  contour(abs(Bspec),4,waxis,waxis),grid
   contour(waxis,waxis,abs(Bspec),4),grid on 
   title('Bispectrum estimated via the direct (FFT) method')
   xlabel('f1'), ylabel('f2')
   set(gcf,'Name','Hosa BISPECD')
   return
end

    
    

     


