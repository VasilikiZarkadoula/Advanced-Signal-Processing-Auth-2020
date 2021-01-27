clc;
clear;

% Construct real discrete signal x[k]

N = 2048;

h = [1; 0.93; 0.85; 0.72; 0.59; -0.1];      % impulse response 
v = exprnd(1,[1,N]);     % white non-Gaussian noise


% convolution
x=zeros(1,N);
for k=1:N
   for i=0:length(h)-1
       if k>i
          x(k)= x(k) + h(i+1)*v(k-i);
       end
   end
end

filename = 'SignalX.mat';
save ('SignalX.mat');