function SDAstats = SDA(filename, diffTS, nfreq, percent_spec, output)
%
%   SDA - Standardized Dispersion Analysis
%   
%   Examples:
%       SDA('fGn1.txt', 0, 128, 50, 1);
%       SDA('fGn2.txt', 0, 128, 50, 1);
%       SDA('fGn3.txt', 0, 128, 50, 1);
%
%   Adpated and updated code from an uncountable number of authors over the years:
%   Jay Holden, Michael J. Richardson, R. C. Schmidt, Charles Coey, Nikita Kuznetsov
%--------------------------------------------------------------------------

%% Load and Check Data
x_data = load(filename);
x = x_data(:,1);
len_p2 = log2(length(x)); 

if length(x) < nfreq*4
   disp(' ERROR: data to short for nFreq...try again'); 
   return;
end

%% Diff TS
if diffTS == 1
   x = diff(x); 
end


%% Z-Score Normalization
if abs(mean(x))>=1e-6 || std(x,1)~=1;
    mu=mean(x); sigma=std(x,1);
    x = (x - mu)./ sigma;
end

    
%% Begin Standardized Dispersion
dispersion=zeros(len_p2,1); % This initalizes the dispersion vector

for j=1:(len_p2-1)
    % computes bin size for each scale
    bin_size=length(x)/pow2(j);
    k = fix(length(x)/bin_size);
    start=1;
    segment_means=zeros(k,1);
    
    for i=1:k  %Compute mean for each bin
        segment_means(i,1)=mean(x(start:start+bin_size-1,1));
        start=start+bin_size;
    end
    
    %Compute std (population formula) across the bins
    dispersion(j,1)=std(segment_means,1);
end

%Get the standard deviation of all points, (bins of 1 data
%point). Should equal 1, so round to emiminate e-15 output
dispersion(len_p2,1)=round(std(x,1)); 
           
           
%% Finish up
%create the bin size vector 
bin_sz=(0:len_p2-1)';

% flip vector since it was done from large to small samples and
% is ploted from small to large samples.
disper=flipud(log2(dispersion));
sdaoutput=[bin_sz disper];
srlength = round(length(sdaoutput)*(percent_spec/100)); 


%% Plot Results
figure;
subplot(1,2,1); % Plot data
plot(x, '-b');
xlim([1 length(x)]);

subplot(1,2,2);
plot(sdaoutput(:,1),sdaoutput(:,2),'bo-','linewidth',2)
hold on;
xlabel('Log_2(Window Size)')
ylabel('Log_2(SD Dispersion)')
p = polyfit(sdaoutput(1:srlength,1),sdaoutput(1:srlength,2),1);
plot(sdaoutput(1:srlength,1),polyval(p,sdaoutput(1:srlength,1),1),'r','linewidth',2);
hold off;

%% Get final SDA statistics
pfit=round(p(1)*100)/100;
SDAstats = [round(pfit*100)/100, round((pfit+1)*100)/100, round((1-pfit)*100)/100];
fprintf('%s      Slope: %.3f    Hurst: %.3f     FD: %.3f\n\n', filename, SDAstats(1), SDAstats(2), SDAstats(3)) 


%% Output to File
if output == 1;
    fid=fopen('SDA_Stats.csv','a');
    fprintf(fid,'%s,',filename);
    fprintf(fid,'%.4f,%.4f,%.4f\r\n',SDAstats(1), SDAstats(2), SDAstats(3));
    fclose(fid);
end
    
    
    
    
    
   