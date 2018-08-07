function AutoCorrel(filename, lag)
%
%   AutoCorrel - Calculates Auto Correlation Function
%   
%   Examples:
%       AutoCorrel('fGn1.txt', 128);
%       AutoCorrel('fGn2.txt', 128);
%       AutoCorrel('fGn3.txt', 128);
%--------------------------------------------------------------------------

%% Load Data
dataRaw = load(filename);
data = dataRaw(:,1);
    
%% Z-score transform
data_z = zscore(data);

%% ACF calculation
c = zeros(lag,1);
c(1) = sum((data_z-mean(data_z)).*(data_z-mean(data_z)))./length(data_z);
for k = 1:lag
    c(k+1) = sum((data_z(1:end-k)-mean(data_z)).*(data_z(1+k:end)-mean(data_z)))./length(data_z);
end
r=c./repmat(c(1),length(c),1);

%% Plot the result
figure;
subplot(1,2,1);
plot(data,'k.-','linewidth',1);
xlabel('Time');
ylabel('Amplitude, A(t)');
xlim([0 length(data)]);

subplot(1,2,2);
stem(0:lag,r,'k','linewidth',1);
xlabel('Lag');
ylabel('Autocorrelation');
xlim([0 lag]);

    
    
    
    
    
   