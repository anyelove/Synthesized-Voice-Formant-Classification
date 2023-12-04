function g=GlottalPulse(N1,N2,N3)
%------------------------------------------------------------
% Glottal Pulse 
%------------------------------------------------------------
% script: Victor Espinoza, 2013 . Angelo Morales, 2022.
%------------------------------------------------------------
% Ref.:
% Rosenberg 1971; Effect of Glottal Pulse Shape on the Quality of Natural Vowels, JASA. 
% Bell Labs
%------------------------------------------------------------
% N1 : opening slope in samples
% N2 : closing slope in samples
% N3 : pulse shape index in Rosenberg Paper
%------------------------------------------------------------

% N1=25; % TP in Rosenberg paper
% N2=10; % TN-TP in Rosenberg paper
% N3= 4; % Pulse shape index in Rosenberg Paper a=1,b=2,c=3,d=4,e=5,f=6

a = 1; % amplitude of function
g = zeros(1 , N1 + N2 + 1); %length of pulse

for t = 1:1:N1 + 1
    if N3 == 1
        g(1,t) = a*((t-1)/N1);
    elseif N3 == 2
        g(1,t) = a*(3*((t-1)/N1)^2 - 2*((t-1)/N1)^3);
    elseif N3 == 3
        g(1,t) = (a/2)*(1-cos(pi*(t-1)/N1));
    elseif N3 == 4
        g(1,t) = (a/2)*(1-cos(pi*(t-1)/N1));
    elseif N3 == 5
        g(1,t) = a*(sin((pi/2)*(t-1)/N1));
    else
        g(1,t) = (3*a/2)*(t/N1);
        if g(1,t) > 1
            g(1,t) = 1;
        else
            g(1,t) = g(1,t);
        end
    end
end

for t = N1 + 2:1:N1 + N2 + 1
    if N3 == 1
        g(1,t) = a*(1-(((t-1)-N1)/N2));
    elseif N3 == 2
        g(1,t) = a*(1-(((t-1)-N1)/N2)^2);
    elseif N3 == 3
        g(1,t) = a*(cos((pi/2)*((t-1)-N1)/N2));
    elseif N3 == 4
        g(1,t) = (a/2)*(1+cos(pi*((t-1)-N1)/N2));
    elseif N3 == 5
        g(1,t) = a*(cos((pi/2)*((t-1)-N1)/N2));
    else
        g(1,t) = (3*a/2)*(1-(((t-1)-N1)/N2));
        if g(1,t) > 1
            g(1,t) = 1;
        else
            g(1,t) = g(1,t);
        end
    end
end
