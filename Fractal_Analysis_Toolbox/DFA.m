function DFA(filename, output)
%
%   DFA - Detrended Fluctuation Analysis
%   
%   NOTES:  DFA can be sued for both fGn and fBm time-seres data.
%           
%           DFA works very well for fGn-like time series data even when the 
%           data are relatively short (256 points).
%           
%           DFA may be problematic for classifying low values of H when
%           the process is fBm (under-diffusive fBm with H from 0 to about 0.2) because
%           they produce negative H estimates! One soluation is to
%           differentiate your fBm data prior to analysis (i.e., convert it to fGn).
%
%   Examples:
%       DFA('fGn1.txt', 1);
%       DFA('fBm1.txt', 1);
%
%   Adpated and updated code from an uncountable number of authors over the years:
%   Jay Holden, Michael J. Richardson, R. C. Schmidt, Charles Coey, Nikita Kuznetsov
%   Jing Hu, Jianbo Gao
%--------------------------------------------------------------------------

%% Load Data
dataRaw = load(filename);
dataRaw = dataRaw(:,1)';

%% Setup Analysis Parameters and Data
dataLength = length(dataRaw);
minscale = 2;
maxscale = floor(log2(dataLength));
step = 1;

%% Center and Integrate Data
dataINT = dataRaw - mean(dataRaw);
dataINT = cumsum(dataINT-mean(dataINT));


%% Divide the integrated time series into boxes of equal length for analysis
%  This is where the DFA action is.
k = 1;
for i = minscale : step : maxscale
    
    l = floor(2^i);
    
    b_block = floor(dataLength/l);
    clear z;
    
    for j = 1 : b_block
        zx = 1 : l;
        zy = dataINT((j - 1)*l + 1 : j*l);
        sx = sum(zx);
        sy = sum(zy);
        sxy = sum(zx.*zy);
        sx2 = sum(zx.^2);
        x_bar = sx/l;
        y_bar = sy/l;
        slope = (l*sxy-sx*sy)/(l*sx2-sx*sx);
        zz = slope*(zx - x_bar) + y_bar;
        z((j - 1)*l + 1 : j*l) = (dataINT((j - 1)*l + 1 : j*l) - zz).^2;
    end
    
    dfa_data(k) = (sum(z)/length(z))^(1/2);
    k = k + 1;
end

%% Organize DFA results
DFAresult(1, :) = minscale : step : maxscale;
DFAresult(2, :) = log2(dfa_data);

%% Plot Results
figure;
subplot(1,3,1:2); % Plot data
plot(dataRaw, '-b');
xlim([1 length(dataRaw)]);

subplot(1,3,3);
plot(DFAresult(1,:),DFAresult(2,:),'.-'); hold on;
ylabel('Log_2(F(n))')
xlabel('Log_2(Window size, n)')
p = polyfit(DFAresult(1,:),DFAresult(2,:),1);
plot(DFAresult(1,:),polyval(p,DFAresult(1,:)),'r','linewidth',2); hold off;

%% Get Scaling
alpha_r = round(p(1)*100)/100;
if alpha_r > 1;
    hurst = alpha_r-1;
    fprintf('%s   Hurst (fBm): %.2f\n\n', filename, alpha_r-1);
else
    hurst = alpha_r;
    fprintf('%s   Hurst (fGn): %.2f\n\n', filename, alpha_r); 
end

%% Output to File
if output == 1;
    fid=fopen('DFA_Stats.csv','a');
    fprintf(fid,'%s,%.4f,%.4f\n',filename, alpha_r, hurst);
    fclose(fid);
end
    
    
    
    
    
   