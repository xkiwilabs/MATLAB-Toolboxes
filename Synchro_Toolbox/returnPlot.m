function returnPlot (filename, samplerate)
%**************************************************************************
%   RETURNPLOT calculates the discrete relative phase between two times
%   series and plots the n:m relative phase on return plot. 
%   Great for determining if coordination is polyrhythmic (i.e., 1:2, 2:3, etc).
%   
%   User needs to specify:
%       filename        : data file to open; should be 2-column txt or csv file
%       samplerate      : sample rate of the time series
%
%   Syntax:
%   returnPlot(filename, samplerate)
%   
%   Examples:
%       >> returnPlot('ExData_TwoToOne.csv', 100);
%       >> returnPlot('ExData_ThreeToTwo.csv', 100);
%   
%   Michael J. Richardson (2005)
%   Last Updated 2009, 2013.
%
%**************************************************************************
 
%% Define Fixed Parameters
linearDetrend = 1;         % 0=no; 1= perform linear detrend (good idea if drift in data)
peakDistance = 0.4;        % 0.5 second minimum period
peakAmp = .3;              % 20% of max amplitude
filterCutoff = 20;         % cutoff frequency for filter (Hz)
rad2deg = 360/(2*pi);      % for converting radians to degrees
 
 
%% Load Data from file
x_data = load(filename);  % should be a 2-column txt or csv file
x1 = x_data(:,1);       % 1 indicates the first columns of data
x2 = x_data(:,2);       % 2 indicates the second column of data
    
 
%% Filter Data using 2nd Order Low-Pass Butterworth Filter
[weight_b,weight_a] = butter(2,filterCutoff/(samplerate/2));
x1 = filtfilt(weight_b,weight_a,x1);
x2 = filtfilt(weight_b,weight_a,x2);
 
 
%% Linear detrend data
if linearDetrend == 1
    x1 = detrend(x1);
    x2 = detrend(x2);
end
 
%% Normalize Data
x1 = x1-mean(x1);
x2 = x2-mean(x2);  
 
 
%% Get Peaks
[~, ~, ~, pLocs1] = period(x1, samplerate, peakDistance, peakAmp);
[~, ~, ~, pLocs2] = period(x2, samplerate, peakDistance, peakAmp);
 
 
%% Get Relative Phase Time Series x1:x2 and x2:x1
[~, ~, ~, radians] = discretephase(x1, x2, samplerate, peakDistance, peakAmp);         
angles1 = radians*rad2deg;
[~, ~, ~, radians] = discretephase(x2, x1, samplerate, peakDistance, peakAmp);         
angles2 = radians*rad2deg;
 
 
%% Define time vector for plotting
delta_t = 1/samplerate; % for time array for plotting
data_len = length(x1);
t = 1:data_len;
t = t*delta_t; 
 
 
%% Plot Data
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/2 scrsz(4)/3]);
h = axes('Position', [0 0 1 1], 'Visible', 'off');
 
% Plot time-series with peaks
axes('Position',[.07 .575 .4 .35]);
hold on;
plot(t, x1,'r-');
plot(pLocs1*delta_t,x1(pLocs1),'r+', 'MarkerSize', 6);
plot(t, x2,'b-');
plot(pLocs2*delta_t,x2(pLocs2),'b+', 'MarkerSize', 6);
xlabel('time (s)');
ylabel('Data');
hold off;
 
% Plot discrete relative phase angles
axes('Position',[.07 .1 .4 .35]);
hold on;
plot(pLocs1*delta_t, angles1,'or');
plot(pLocs2*delta_t, angles2,'ob');
ylim([-200 200]);
set(gca,'YTick',[-180 0 180])
xlabel('time (s)');
ylabel('DRP');
hold off;
 
%Plot relative phase return plot
axes('Position',[.55 .1 .4 .825]);
hold on;
plot(angles1(1:end-1), angles1(2:end), 'or', 'MarkerSize', 4);
plot(angles2(1:end-1), angles2(2:end), 'ob', 'MarkerSize', 4);
xlim([-200 200]);
set(gca,'XTick',[-180 -90 0 90 180])
ylim([-200 200]);
set(gca,'YTick',[-180 -90 0 90 180])
xlabel('RPt');
ylabel('RPt+1');
hold off;
 
%% end of function
return
 


