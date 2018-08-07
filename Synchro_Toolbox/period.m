function [meanPeriod, sdPeriod, peaks, pkLocs] = period(data, samplerate, DistFLT, AmpFLT)
%************************************************************************************
%   PERIOD calculates the periods for rhythmic time-series data.
%
%   Returns the peaks and peaks locations for rhythmic time-series data,
%   as well as the mean period and SD period in seconds.
%
%   Syntax:
%   [meanPeriod, sdPeriod, peaks, pkLocs] = period(data, 100, .5, .3)                                                                    
%
%   BY: Michael J. Richardson (michael.richardson@uc.edu), 2004, updated 2009 
%
%-------------------------------------------------------------------------------------
 
%% Center Data
x=data(:,1)-mean(data(:,1));
 
%% Determine Peaks and Peak Locations
[peaks, pkLocs] = findpeaks(x,'MinPeakDistance',floor(samplerate*DistFLT),'MinPeakHeight',(max(x))*AmpFLT);
 
%% Get Period Stats
meanPeriod = mean(diff(pkLocs))/samplerate;
sdPeriod = std(diff(pkLocs))/samplerate;
 
%% End of function
return
%---------------------------------------------------------------------------------------