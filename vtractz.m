function zp = vtractz(p,bw,fs)

% This MATLAB function produces a vector of z-domain poles corresponding
% to the input pole frequencies and bandwidths, assuming an all-pole
% model of the vocal tract.  It takes vector input of pole frequencies
% and bandwidths (in Hz), and a scalar sampling frequency fs, and
% outputs a vector containing the z-domain poles for the vocal tract.
%
% USAGE:   zp = vtractz(p,bw,fs)
%
% Input:   p = column vector of the system pole frequencies (Hz)
%          bw = column vector of the system pole bandwidths (Hz)
%          fs = sampling frequency (Hz)
%
% Output:  zp = vocal tract model z-domain poles
%
% Harold Cheyne
% 7 December 2000
% revised 14 December 2000
% 
% eg.
% zp = vtractz(1000,500,8000)
% 
% zp =
% 
%    0.5814 + 0.5814i
%    0.5814 - 0.5814i

% The z-domain poles are of the form a+jb, where a and b are between 0 and 1.
% And each pole needs to appear in a pair, a+jb and a-jb.
n = length(p);										% get the number of poles
zp = zeros(2.*n,1);								% initialize the output vector
B = pi.*bw./fs;									% Express the bandwidths in radians
P = 2.*pi.*p./fs;									% Express the pole frequencies in radians
r = -cos(B)+2-sqrt((cos(B)-1).*(cos(B)-3));	% Calculate the z-plane pole magnitudes
for m = 1:n										% For each input pole frequency, calculate
   if p(m) < fs/4									% the real part, which is positive for poles
      a = r(m)./sqrt(1+tan(P(m)).^2);		% in the 1st quadrant and negative for poles
   else a = -r(m)./sqrt(1+tan(P(m)).^2);	% in the 2nd quadrant.
   end
   b = a.*tan(P(m));								% Calculate the imaginary part.
   zp(2*m-1) = a + 1i.*b;
   zp(2*m) = a - 1i.*b;
end
