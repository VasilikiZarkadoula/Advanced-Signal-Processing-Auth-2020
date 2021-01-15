% ACF - Estimate AutoCorrelation Function
% >> r = acf(y,p)
%
% Inputs:
% y - series to compute acf for, Nx1 vector
% p - total number of lags, real positive scalar integer < N
%
% Output:
% r: A px1 vector containing autocorrelations
%    (First lag computed is lag 1. Lag 0 not computed)
%
%
% Examples:
% >> acf(randn(100,1), 10)
%
% >> plot(acf(filter(ones(10,1),1,randn(1e6,1)),100));
%
% ===============================
% Barry Cardiff
% University College Dublin (UCD)
% email: barry.cardiff@ucd.ie
% ===============================
%
% Algorithm:
% The def of the ACF is:
%               r(m) = E[conj(y(n)-E[y])*(y(n+m)-E[y])]/Var[y]
% as the mean does not effect the variance, the mean can be removed from
% the signal simplying the calculation:
% step 1)   y <--- y - mean(y)
%           then we have:
%           r(m) = E[conj(y(n))*(y(n+m))]/Var[y]
%
% step 2)   The Expectation is replaced by a flat averag, and the variance
%           is estimated using matlab var() function
% 
%
% Version History:
% v1: 26th June 2017
% Based on ACF downloaded from matworks' file exchange by Calvin Price,
% but I've made it work for complex data, made it moe efficient (removed a
% loop), removed the ploting stuff, and improved / fixed argument checking.
%
function r = acf(y,p)
    if nargin > 2
        error('acf usage error, see help.\n');
    end

    if ~isvector(y)
        error('Input series y must be a vector.\n');
    end
    N = length(y);

    if ~(isscalar(p) && (p<N) && isreal(p) && (p>0) && (floor(p)==p) )
        error('Input number of lags p must be a real positive integer < N.\n');
    end

    % remove mean as it does not impact the acf, and make colunm vector
    y = y(:) - mean(y);
    v = var(y); % estimate the variance
    
    % Loop over all lags
    r = nan(p,1); % preload for speed
    for m = 1:p
       r(m) = (1/(N-m))*((y(1:(end-m)))'*y((1+m):end))/v;
    end

return
