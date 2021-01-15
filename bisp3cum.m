function [bisp,freq,cum,lag]=bisp3cum(signal,samprate,maxlag,window,scale)

%	BISP3CUM Auto bispectrum/3rd order cumulant
%
%	[bisp,freq,cum,lag]=bisp3cum(signal,samprate,maxlag,window,scale)
%
%	The maxlag*2+1 x maxlag*2+1 element auto bispectrum and 3rd order cumulant matrices
%	and maxlag*2+1 element frequency and lag vectors are computed from the signal
%	matrix containing samples in rows and records in columns, signal sample rate and
%	maximum lag scalars, and lag window function and scale strings.
%
%	If unspecified, the signal matrix is entered after the prompt from the keyboard,
%	and the default assignments samprate=1 and maxlag=0 are used.  The window and scale
%	strings specify lag window and scale matrix computation, according to:
%
%	window = 'none', 'n', or unspecified does not compute a window
%	       = 'uniform' or 'u' computes the uniform hexagonal window
%	       = 'sasaki' or 's' computes the sasaki window
%	       = 'priestley' or 'p' computes the priestley window
%	       = 'parzen' or 'pa' computes the parzen window
%	       = 'hamming' or 'h' computes the hamming window
%	       = 'gaussian' or 'g' computes the gaussian distribution window
%	       = 'daniell' or 'd' computes the daniell window
%
%	scale  = 'biased', 'b', or unspecified computes the biased estimate
%	       = 'unbiased' or 'u' computes the unbiased estimate

%	Implemented using MATLAB 5.3.1 and additional functions:
%
%	mat=toep(column,row)
%	wind=lagwind(lag,window)
%
%	Implementation:
%
%	cum(k,l) = sum_{n=0}^{N-1} conj(signal(n))*signal(n+k)*signal(n+l)/N
%
%	k,l = {-maxlag,...,-1,0,1,...,maxlag}, n = {0,1,...,N-1}
%
%	bisp=fftshift(fft2(ifftshift(cum.*wind)))
%
%	Example:
%
%	» [b,f,c,l]=bisp3cum([1-i -1+i],1,1)
%
%	b =
%
%	  -5.1962 - 5.1962i        0            -0.0000 - 0.0000i
%	        0                  0                  0          
%	  -0.0000 - 0.0000i   0.0000 + 0.0000i   5.1962 + 5.1962i
%
%	f =
%
%	   -0.5000         0    0.5000
%
%	c =
%
%	  -1.0000 + 1.0000i   1.0000 - 1.0000i        0          
%	   1.0000 - 1.0000i        0            -1.0000 + 1.0000i
%	        0            -1.0000 + 1.0000i   1.0000 - 1.0000i
%
%	l =
%
%	    -1     0     1
%
%	References:
%
%	C. L. Nikias, A. P. Petropulu, Higher-Order Spectra Analysis:  A Nonlinear Signal
%	Processing Framework, PTR Prentice Hall, Englewood Cliffs, NJ, 1993.
%
%	T. S. Rao, M. M. Gabr, An Introduction to Bispectral Analysis and Bilinear Time
%	Series Models, Lecture Notes in Statistics, Volume 24, D. Brillinger, S. Fienberg,
%	J. Gani, J. Hartigan, K. Krickeberg, Editors, Springer-Verlag, New York, NY, 1984.
%
%	Copyright (c) 2000
%	Tom McMurray
%	mcmurray@teamcmi.com

%	assign default input parameters

if ~nargin
   signal=input('enter signal matrix or return for 0 outputs\n');
   if isempty(signal)
      bisp=0;
      freq=0;
      cum=0;
      lag=0;
      return
   end
end
if nargin<2
   samprate=1;
end
if nargin<3
   maxlag=0;
end
if nargin<4
   window='n';
end
if nargin<5
   scale='b';
end

%	while signal is unsupported, enter supported signal or return for 0 outputs

while isempty(signal)|~isnumeric(signal)|~all(all(isfinite(signal)))|ndims(signal)>2
   signal=input(['signal is empty, nonnumeric, nonfinite, or > 2 '...
         'dimensional:\nenter finite matrix or return for 0 outputs\n']);
   if isempty(signal)
      bisp=0;
      freq=0;
      cum=0;
      lag=0;
      return
   end
end
[sample,record]=size(signal);

%	if signal is a row vector, modify to a column vector

if sample==1
   sample=record;
   record=1;
   signal=signal(:);
end

%	while samprate is unsupported, enter supported samprate or return for samprate=1

while isempty(samprate)
   samprate=input(['signal sample rate is empty:\nenter finite positive scalar or '...
         'return for signal sample rate = 1\n']);
   if isempty(samprate)
      samprate=1;
   end
end
while samprate<=0|~isnumeric(samprate)|~isfinite(samprate)|length(samprate)~=1
   samprate=input(['signal sample rate = ' num2str(samprate(:).') ' <= 0, '...
         'nonnumeric, nonfinite, or nonscalar:\nenter finite positive scalar or '...
         'return for signal sample rate = 1\n']);
   if isempty(samprate)
      samprate=1;
   end
end
sample1=sample-1;
strsample1=num2str(sample1);

%	while maxlag is unsupported, enter supported maxlag or return for maxlag=0

while isempty(maxlag)
   maxlag=input(['maximum lag is empty:\nenter integer scalar >= 0, <= signal '...
         'sample length - 1 = ' strsample1 ', or return for maximum lag = 0\n']);
   if isempty(maxlag)
      maxlag=0;
   end
end
while maxlag<0|maxlag>sample1|rem(maxlag,1)|~isnumeric(maxlag)|~isfinite(maxlag)...
      |length(maxlag)~=1
   maxlag=input(['maximum lag = ' num2str(maxlag(:).') ' < 0, > signal sample '...
         'length - 1 = ' strsample1 ', noninteger, nonnumeric, nonfinite, or '...
         'nonscalar:\nenter integer scalar >= 0, <= ' strsample1 ', or return for '...
         'maximum lag = 0\n']);
   if isempty(maxlag)
      maxlag=0;
   end
end

%	compute lag vector

lagindex=-maxlag:maxlag;
lag=lagindex/samprate;

%	if maxlag, compute freq vector

if maxlag
   freq=lagindex/maxlag/2*samprate;
   
%	else, freq=0 and specify no window/biased estimate computation
   
else
   freq=0;
   window='n';
   scale='b';
end
clear lagindex

%	resolve window

window=lower(num2str(window));
windowarr={'none' 'uniform' 'sasaki' 'priestley' 'parzen' 'hamming' 'gaussian'...
      'daniell'};
windowind=strmatch(window,windowarr);

%	while window is unsupported, enter supported window or return for window='n'

while isempty(windowind)|isempty(window)
   window=lower(input(['window = ' window(:).' ' unresolved:\nenter none, '...
         'uniform, sasaki, priestley, parzen, hamming, gaussian, daniell, or '...
         'return for window = none\n'],'s'));
   if isempty(window)
      window='n';
   end
   windowind=strmatch(window,windowarr);
end

%	if window is unique, assign supported window

if length(windowind)==1
   window=windowarr{windowind};
   
%	else window=='p', window='priestley'
   
else
   window='priestley';
end
clear windowarr

%	resolve scale

scale=lower(num2str(scale));
scalearr={'biased' 'unbiased'};
scaleind=strmatch(scale,scalearr);

%	while scale is unsupported, enter supported scale or return for scale='b'

while isempty(scaleind)|isempty(scale)
   scale=lower(input(['scale = ' scale(:).' ' unresolved:\nenter biased, '...
         'unbiased, or return for scale = biased\n'],'s'));
   if isempty(scale)
      scale='b';
   end
   scaleind=strmatch(scale,scalearr);
end
scale=scalearr{scaleind};
scale=scale(1);
clear scalearr

%	generate constants

maxlag1=maxlag+1;
maxlag2=maxlag*2;
maxlag21=maxlag2+1;
samp1ind=sample:-1:1;
samlsamind=sample-maxlag:sample;
ml1samind=maxlag1:sample;
ml211ind=maxlag21:-1:1;
zeros1maxlag=zeros(1,maxlag);
zerosmaxlag1=zeros1maxlag(:);
onesmaxlag211=ones(maxlag21,1);
strmaxlag21=num2str(maxlag21);

%	subtract mean from signal

meansig=mean(signal);
signal=signal-meansig(ones(sample,1),:);

%	initialize cumulant matrix

cum=zeros(maxlag21);

%	signal record cumulant computation loop

tic
for k=1:record
   time=cputime;
   sig=signal(:,k);
   trflsig=sig(samp1ind)';
   toepsig=toep([sig(samlsamind);zerosmaxlag1],[conj(trflsig(ml1samind)) ...
         zeros1maxlag]);
   
%	compute cumulant
   
   cum=cum+toepsig.*trflsig(onesmaxlag211,:)*toepsig.';
   disp(['record ' num2str(k) ':  time = ' num2str(cputime-time) ' seconds'])
end
cum=cum/record;
clear samp1ind samlsamind ml1samind zerosmaxlag1 sig trflsig toepsig

%	if scale=='b', compute biased cumulant

if scale=='b'
   cum=cum/sample;
   
%	else, compute unbiased cumulant

else
   scalmat=zeros(maxlag1);
   for k=1:maxlag1
      maxlag1k=maxlag1-k;
      scalmat(k,k:maxlag1)=repmat(sample-maxlag1k,1,maxlag1k+1);
   end
   scalmat=scalmat+triu(scalmat,1).';
   samplemaxlag1=sample-maxlag1;
   maxlag1ind=maxlag:-1:1;
   scalmat=[scalmat [toep((samplemaxlag1:sample-2).',...
            samplemaxlag1:-1:sample-maxlag2);scalmat(maxlag1,maxlag1ind)]];
   scalmat=[scalmat;scalmat(maxlag1ind,ml211ind)];
   scalmat(scalmat<1)=1;
   cum=cum./scalmat;
   clear scalmat maxlag1ind
end
time=num2str(toc);
disp(' ')
disp([strmaxlag21 ' x ' strmaxlag21 ' element cumulant computed in ' time ' seconds'])

%	generate lag window function

if window(1)=='n'
   wind=1;
else
   wind=lagwind(maxlag1,window);
   
%	generate 2d even window function
   
   windeven=[wind(maxlag1:-1:2) wind];
   windeven=windeven(onesmaxlag211,:);
   
%	0 pad window function
   
   wind=[wind zeros1maxlag];

%	generate 2d window function
   
   wind=toep(wind(:),[wind(1) zeros(1,maxlag2)]);
   wind=wind+tril(wind,-1).';
   wind=wind(ml211ind,:).*windeven.*windeven.';
   clear windeven
end
clear ml211ind zeros1maxlag onesmaxlag211

%	compute bispectrum

tic
bisp=fftshift(fft2(ifftshift(cum.*wind)));
time=num2str(toc);
disp(' ')
disp([strmaxlag21 ' x ' strmaxlag21 ' element bispectrum computed in ' time...
      ' seconds'])
disp(' ')
clear wind

%	plot mean signal

meansig=mean(signal,2);
realsig=isreal(meansig);
subplot 221
if realsig
   plot((0:sample1)/samprate,meansig)
else
   plot((0:sample1)/samprate,abs(meansig))
end
fontsize='\fontsize{8}';
seconds='(\its \rm)';
title([fontsize 'Averaged Signal'])
xlabel([fontsize 'Time ' seconds])
ylabel([fontsize 'Signal (\itV \rm)'])
grid
clear meansig

%	plot 3rd order cumulant

subplot 222
if realsig
   imagesc(lag,lag,cum)
else
   imagesc(lag,lag,abs(cum))
end
title([fontsize '3^{rd} Order Cumulant (\itV \rm^{3} )'])
xlabel([fontsize 'Lag \tau_{0} ' seconds])
ylabel([fontsize 'Lag \tau_{1} ' seconds])
axis xy
grid
colorbar

%	plot bispectrum

subplot 223
imagesc(freq,freq,abs(bisp))
imagesc(freq,freq,abs(bisp))
title([fontsize 'Bispectrum Magnitude (\itV \rm^{3}/\itHz \rm^{2} )'])
hertz='(\itHz \rm)';
fxlabel=[fontsize 'Frequency f_{0} ' hertz];
xlabel(fxlabel)
fylabel=[fontsize 'Frequency f_{1} ' hertz];
ylabel(fylabel)
axis xy
grid
colorbar
subplot 224
imagesc(freq,freq,angle(bisp)*180/pi)
title([fontsize 'Bispectrum Phase (\itdeg \rm)'])
xlabel(fxlabel)
ylabel(fylabel)
axis xy
grid
colorbar
colormap gray