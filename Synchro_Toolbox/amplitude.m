function [meanAmp, sdAmp, peaks, pkLocs, valleys, vLocs] = amplitude(data, samplerate, DistFLT, AmpFLT)
%************************************************************************************
%   PERIOD    Calculates the amplitude for rhythmic time-series data.
%
%   Returns the peaks and valleys and peaks and valley locations for a 
%   rhythmic positional time-series, as well as the mean and SD amplitude.
%
%   Syntax:
%   [meanAmp, sdAmp, peaks, pkLocs, valleys, vLocs] = amplitude(data, 100, .5, .3)                                                                    
%
%   BY: Michael J. Richardson (michael.richardson@uc.edu), 2004, updated 2009     
%
%-------------------------------------------------------------------------------------
 
%% Center Data
x=data(:,1)-mean(data(:,1));
 
%% Determine Peaks and Peak Locations
[peaks, pkLocs] = findpeaks(x,'MinPeakDistance',floor(samplerate*DistFLT),'MinPeakHeight',(max(x))*AmpFLT);
 
%% Determine Valleys and Valley Locations
[valleys, vLocs] = findpeaks(x.*-1,'MinPeakDistance',floor(samplerate*DistFLT),'MinPeakHeight',(max(x))*AmpFLT);
valleys = valleys.*-1;
 
%% Determine appropriate data length for Amp calculations
pvLength = length(peaks);
if pvLength > length(valleys) 
    pvLength = length(valleys);
end
 
%% Calculate mean and SD
meanAmp = mean(peaks(1:pvLength)-valleys(1:pvLength));
sdAmp = std(peaks(1:pvLength)-valleys(1:pvLength));
 
%% End of function
return
%---------------------------------------------------------------------------------------