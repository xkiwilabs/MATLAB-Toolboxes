function [meanRP, sdRP, rvRP, radians] = discretephase(x1, x2, samplerate, DistFLT, AmpFLT)
%************************************************************************************
%   DISCRETEPHASE calculates discrete instantaneous relative phase angle 
%   between two time-series (using Hilbert transform), with time-series 1 as 
%   the referent. Returns the discrete relative phase time-series in radians 
%   as well as the circular mean and SD of the discrete relative phase in 
%   radians and the mean resultant vector (rvRP = 0 to 1)
%
%   Syntax:
%   [meanRP, sdRP, rvRP, radians] = discretephase(x1, x2)                                                                      
%
%   BY: Michael J. Richardson (michael.richardson@uc.edu), 2007     
%       Updated to include circular stats in 2011
%
%-------------------------------------------------------------------------------------
 
%% Center Data
x1=x1-mean(x1);
x2=x2-mean(x2);
 
%% Get Peak locations for first time-Series (i.e., x1);
[~, locs] = findpeaks(x1,'MinPeakDistance',floor(samplerate*DistFLT),'MinPeakHeight',(max(x1))*AmpFLT);
 
%% Using Hilbert Transform to Calculate Phase Angles for Second time-series
hbt = hilbert(x2);
ihbt = imag(hbt);
radians = atan2(ihbt,x2);
 
%% Calculate Phase Angle of second time-series at peak locations of first time-series.
radians = radians(locs);
 
%% Get Circular Stats
wgh = ones(size(radians));
wr = wgh'*exp(1i*radians);
meanRP = angle(wr);
rvRP = abs(wr)/sum(wgh);
sdRP = sqrt(2*(1-rvRP));
 
%% End of function
return
%---------------------------------------------------------------------------------------


