function wind=lagwind(lag,window)

%	LAGWIND Lag window function
%
%	wind=lagwind(lag,window)
%
%	The lag element window vector is computed from the lag scalar and window function
%	string.  If unspecified, the lag scalar is entered after the prompt from the
%	keyboard.  The window string specifies lag window vector computation, according to:
%
%	window = 'uniform', 'u', or unspecified computes the uniform window
%	       = 'sasaki' or 's' computes the sasaki window
%	       = 'priestley' or 'p' computes the priestley window
%	       = 'parzen' or 'pa' computes the parzen window
%	       = 'hamming' or 'h' computes the hamming window
%	       = 'gaussian' or 'g' computes the gaussian distribution window
%	       = 'daniell' or 'd' computes the daniell window

%	Implemented using MATLAB 5.3.1
%
%	Example:
%
%	» w=lagwind(4,'parzen')
%
%	w =
%
%	    1.0000    0.5556    0.0741         0
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
   lag=input('enter lag integer scalar > 1 or return for wind = 1\n');
   if isempty(lag)|lag==1
      wind=1;
      return
   end
end
if nargin<2
   window='uniform';
end

%	while lag is unsupported, enter supported lag or return for wind=1

while isempty(lag)
   lag=input(['lag is empty:\nenter integer scalar > 1 or return for wind = 1\n']);
   if isempty(lag)|lag==1
      wind=1;
      return
   end
end
while lag<1|rem(lag,1)|~isnumeric(lag)|~isfinite(lag)|length(lag)~=1
   lag=input(['lag = ' num2str(lag(:).') ' < 1, noninteger, nonnumeric, '...
         'nonfinite, or nonscalar:\nenter integer scalar > 1 or return for wind = '...
         '1\n']);
   if isempty(lag)|lag==1
      wind=1;
      return
   end
end

%	if lag==1, wind=1

if lag==1
   wind=1;
   return
end

%	resolve window

window=lower(num2str(window));
windowarr={'uniform' 'sasaki' 'priestley' 'parzen' 'hamming' 'gaussian' 'daniell'};
windowind=strmatch(window,windowarr);

%	while window is unsupported, enter supported window or return for window='uniform'

while isempty(windowind)|isempty(window)
   window=lower(input(['window = ' window(:).' ' unresolved:\nenter uniform, '...
         'sasaki, priestley, parzen, hamming, gaussian, daniell, or return for '...
         'window = uniform\n'],'s'));
   if isempty(window)
      window='uniform';
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
lag1=lag-1;

%	compute lag window vector

switch window
   
%	uniform window

case 'uniform'
   wind=ones(1,lag);
   
%	sasaki window
   
case 'sasaki'
   windlag=(0:lag1)/lag1;
   wind=sin(pi*windlag)/pi+cos(pi*windlag).*(1-windlag);
   
%	priestley window
   
case 'priestley'
   windlag=(1:lag1)/lag1;
   wind=[1 (sin(pi*windlag)/pi./windlag-cos(pi*windlag))*3/pi/pi./windlag./windlag];
   
%	parzen window
   
case 'parzen'
   fixlag12=fix(lag1/2);
   fixlag121=fixlag12+1;
   windlag0=(0:fixlag12)/lag1;
   windlag1=1-(fixlag121:lag1)/lag1;
   wind(1:fixlag121)=1-(1-windlag0).*windlag0.*windlag0*6;
   wind(fixlag121+1:lag)=windlag1.*windlag1.*windlag1*2;
   
%	hamming window
   
case 'hamming'
   wind=0.54+0.46*cos(pi*(0:lag1)/lag1);
   
%	0.5 mean, 0.125 standard deviation, truncated gaussian distribution window
   
case 'gaussian'
   wind=[1 erfc(((1:lag1-1)/lag1-0.5)*8/sqrt(2))/2 0];
   
%	daniell window
   
case 'daniell'
   windlag=(1:lag1)/lag1;
   wind=[1 sin(pi*windlag)/pi./windlag];
   
%	end switch window
   
end