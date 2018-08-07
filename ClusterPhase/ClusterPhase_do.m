function [GRPrhoM INDrhoM INDrpM TSrhoGRP TSrpIND] = ClusterPhase_do(TSfilename, TSnumber, TSfsamp, TSlsamp, TSsamplerate, plotflag)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   ClusterPhase_do.m
%
%   [GRPrhoM INDrhoM INDrpM TSrhoGRP TSrpIND] = ClusterPhase_do(TSfilename, TSnumber, TSfsamp, TSlsamp, TSsamplerate, plotflag)
%
%   Input:
%       TSfilename    : data file
%       TSnumber      : number of time series
%       TSsamplerate  : sample rate of the time series
%       TSfsamp       : first data point in time series used
%       TSlsamp       : last data point in time series used
%       plotflag      : do plots (0=no, 1=yes)
%
%   Output:
%       GRPrhoM        : mean group rho (0 to 1; 1 = perfect sync)
%       INDrhoM        : mean rho for each TS to group (0 to 1; 1 = perfect sync)
%       INDrpM         : mean Relative Phase for each TS to group cluster phase 
%       TSrhoGRP        : group rho time-series
%       TSrpIND         : relative phase time-series for each individual TS to cluster phase
%
%   Example:
%       [GRPrhoM INDrhoM INDrpM TSrhoGRP TSrpIND] = ClusterPhase_do('G201EO1.txt', 6, 1, 7200, 120, 1);
%
%   BY (2008):
%   Michael J Richardson (Univeristy of Cincinnati) & Till D. Frank (UCONN) 
%   
%   UPDATED (2011):
%   Michael J Richardson (Univeristy of Cincinnati)
%
%   References:
%   [1]  Frank, T. D., & Richardson, M. J. (2010). On a test statistic for 
%        the Kuramoto order parameter of synchronization: with an illustration 
%        for group synchronization during rocking chairs.
%
%   [2]  Richardson,M.J., Garcia, R., Frank, T. D., Gregor, M., & 
%        Marsh,K. L. (2010). Measuring Group Synchrony: A Cluster-Phase Method 
%        for Analyzing Multivariate Movement Time-Series 
%
%   Code Contact & References:
%        michael.richardson@uc.edu
%        http://homepages.uc.edu/~richamo/
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%% Set fixed parameters
filterfreq = 10;
close all;


%% load time-series (TS)
%**************************************************************************
data = load(TSfilename);

for nts=1:TSnumber
    ts_data(:,nts) = data(TSfsamp:TSlsamp,nts);
end

TSlength = length(ts_data(:,1));
delta_t = 1/TSsamplerate;
t = (1:TSlength)*delta_t;


%% Downsample, Filter and normalize data
%**************************************************************************
%linear detrend data to remove drift (chiar moving slightly during trial)
for nts=1:TSnumber
    ts_data(:,nts) = detrend(ts_data(:,nts));
end

%normlaize data
for nts=1:TSnumber
    ts_data(:,nts) = zscore(ts_data(:,nts));
end

%filter data
for nts=1:TSnumber
    [weight_b,weight_a] = butter(2,filterfreq/(TSsamplerate/2));
    ts_data(:,nts) = filtfilt(weight_b,weight_a,ts_data(:,nts));
end


%% Compute phase for each TS using Hilbert transform
%**************************************************************************
TSphase = zeros(TSlength-1,TSnumber);
for k=1:TSnumber
    hrp = hilbert(ts_data(:,k));
    for n=1:TSlength-1
        TSphase(n,k)=atan2(real(hrp(n)),imag(hrp(n)));
    end
    TSphase(:,k)=unwrap(TSphase(:,k));
end


%% Compute mean running (Cluster) phase
%**************************************************************************
clusterphase = zeros(1,TSlength-1);
for n=1:TSlength-1
    ztot=complex(0,0);
    for k=1:TSnumber
        z=exp(1i*TSphase(n,k));
        ztot=ztot+z;
    end
    ztot=ztot/TSnumber;
    clusterphase(n)=angle(ztot);
end
clusterphase = unwrap(clusterphase);


%% Compute relative phases between phase of TS and cluster phase
%**************************************************************************
TSrpIND=zeros(TSlength-1,TSnumber);
INDrpM = zeros(TSnumber,1);
INDrhoM = zeros(TSnumber,1);
for k=1:TSnumber
    ztot=complex(0,0);
    for n=1:TSlength-1
        z=exp(1i*(TSphase(n,k)-clusterphase(n)));
        TSrpIND(n,k) = z;
        ztot=ztot+z;
    end
    TSrpIND(:,k) = angle(TSrpIND(:,k))*360/(2*pi); % convert radian to degrees
    ztot=ztot/(TSlength-1);
    INDrpM(k) = angle(ztot);
    INDrhoM(k) = abs(ztot);
end
TSRPM = INDrpM;
INDrpM = (INDrpM(:,1)./(2*pi)*360); % convert radian to degrees
disp(' ');
disp('Mean relative phases of individuals to cluster phase')
disp(INDrpM(:,1).');
disp('Averaged degree of synchronization of individuals (Rho = 1-circular variance)')
disp(INDrhoM(:,1).');


%% Compute cluster amplitude rhotot in rotation frame
%**************************************************************************
TSrhoGRP=zeros(TSlength-1,1);
for n=1:TSlength-1
    ztot=complex(0,0);
    for k=1:TSnumber
        z=exp(1i*(TSphase(n,k)-clusterphase(n)-TSRPM(k)));
        ztot=ztot+z;
    end
    ztot=ztot/TSnumber;
    TSrhoGRP(n)=abs(ztot);
end

GRPrhoM = mean(TSrhoGRP);
disp('Averaged degree of synchronization of the group')
disp(GRPrhoM);


%% Do Plot
%**************************************************************************
% plot data for time-series (separeted on graph for display purposes)
scrsz = get(0,'ScreenSize');
h = figure('Position',[scrsz(3)/3 scrsz(4)/3 scrsz(3)/2 scrsz(4)/2]);

if plotflag == 1
    subplot(3,1,1);
    hold on;
    tmpdata = zeros(TSlength,TSnumber);
    for nts=1:TSnumber
        tmpdata(:,nts) = (ts_data(:,nts)+(nts*4));
    end
    plot(t(1:TSlength), tmpdata);
    xlabel('Time','fontsize',10)
    ylabel('RAW Data','fontsize',10)
    hold off;
end

% plot individ-cluster relative phase
if plotflag == 1
    subplot(3,1,2);
    set(gca,'fontsize',10)
    plot(t(1:length(t)-1), TSrpIND)
    xlabel('Time','fontsize',10)
    ylabel('IND-Clust Relative Phase','fontsize',10)
    xlim([0 max(t)]);
    ylim([-185 185]);
end
get(0, 'CurrentFigure');
str(1) = {['Mean RPs: ', sprintf('%.2f  ', INDrpM(:,1))]};
text(0, 220, 0, str, 'FontSize', 10, 'Color', 'k');


% plot group-cluster amplitiude (rho) timeseries
if plotflag == 1
    subplot(3,1,3);
    set(gca,'fontsize',10)
    plot(t(1:length(t)-1), TSrhoGRP)
    xlabel('Time','fontsize',10)
    ylabel('GRP-Clust Amplitude','fontsize',10)
    ylim([0 1]);
    xlim([0 max(t)]);
end
get(0, 'CurrentFigure');
str(1) = {['Mean GRP Rho: ', sprintf('%.3f  ', GRPrhoM), '  Mean IND Rhos: ', sprintf('%.3f  ', INDrhoM(:,1))]};
text(0, -.4, 0, str, 'FontSize', 10, 'Color', 'k');	


%%
return;
%**************************************************************************
%**************************************************************************
%**************************************************************************
