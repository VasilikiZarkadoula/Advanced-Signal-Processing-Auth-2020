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

%--------------------------------------------------------------------------

%BISPECD function from HOSA toolbox
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



