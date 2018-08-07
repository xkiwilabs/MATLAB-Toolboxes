function [meanRP, sdRP, rvRP, radians] = hilbertphase(x1, x2)
%************************************************************************************
%   HILBERTPHASE calculates instantaneous relative phase angle between
%   two time-series. Returns the instantaneous relative phase time-series in 
%   radians (using Hilbert transform), the circular mean and SD of the
%   relative phase in radians, and the mean resultant vector (rvRP = 0 to 1)
%   
%
%   Syntax:
%   [meanRP, sdRP, rvRP, radians] = hilbertphase(x1, x2)                                                                      
%
%   BY: Michael J. Richardson (michael.richardson@uc.edu), 2004     
%       Updated to include circular stats in 2009.
%
%-------------------------------------------------------------------------------------
 
%% Center Data
x1=x1-mean(x1);
x2=x2-mean(x2);
 
%% Using Hilbert Transform to Calculate Relative Phase
h1 = hilbert(x1);
h2 = hilbert(x2);
h1 = imag(h1);
h2 = imag(h2);
num = h2.*x1 - x2.*h1;
denom = x2.*x1 + h2.*h1;
radians = atan2(num,denom);
 
%% Get Circular Stats
wgh = ones(size(radians));
wr = wgh'*exp(1i*radians);
meanRP = angle(wr);
rvRP = abs(wr)/sum(wgh);
sdRP = sqrt(2*(1-rvRP));
 
%% End of function
return
%---------------------------------------------------------------------------------------
