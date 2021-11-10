function seis_tap = Haning_kai(seis)
% Function for applying a cosine taper on spesified percentage of each end 
% of the input signal.
%
% Input:
%       seis = the signal to be tapered

%
% Output:
%       seis_tap = the tapered signal
%
% Written by Kai 
% 

[m ,n]= size(seis); % length of the seismogram

% The first part of the costap(t) function:
theta=pi*(0:n-1)/(n-1);
costap = 0.5*(cos(theta)+1);
seis_tap = seis.*costap;   


