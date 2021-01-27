function [x,v,h,N] = Signals_xv()

N=2048;
h=[1; 0.93; 0.85; 0.72; 0.59; -0.1];
rng shuffle
v = exprnd(1,[1,N]);
x=zeros(1,N);
for k=1:N
        for i=0:5
            if k>i
                x(k)= x(k)+h(i+1)*v(k-i);
            end
        end
end

end