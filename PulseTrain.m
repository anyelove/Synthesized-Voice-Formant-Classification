function x=PulseTrain(N,P,jitter,shimmer)
% Pulse Train Generator
% N : Size Vector in samples
% P : Period of Pulse in samples
% jitter : in percentage of period
% shimmer : in percentage of amplitude
% script: Victor Espinoza

x=[1 zeros(1,2*N)];
m=0:floor(N/P)-1;
for ii=m
    Prand=round(P-P*jitter*randn(1,1)/100);
    x(Prand+P*ii+1)=1*(1-randn(1,1)*shimmer/100);
end
x=x(1:N);
