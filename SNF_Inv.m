function [x, b, a]=SNF_Inv(x,fc,fs,BW,g)

zp = vtractz(fc,BW,fs);                 % Get z-domain poles from fp & b
[b,a] = zp2tf([],zp,1);
if g==1
  b=sum(a);
else
  b=g*sum(a);
end
x=filter(b,a,x);