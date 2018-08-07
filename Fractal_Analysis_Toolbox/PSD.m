function PSD(filename, diffTS, nfreq, percent_spec, output)
%
%   PSD - Power Spectral Density Analysis
%	Compute power spectrum using Welch's method from a single time-series
%	of a length of a power of two.
%   
%   Examples:
%       PSD('fGn1.txt', 0, 128, 75, 0);
%       PSD('fGn2.txt', 0, 128, 75, 0);
%       PSD('fGn3.txt', 0, 128, 75, 0);
%       PSD('fBm1.txt', 0, 128, 75, 1);
%       PSD('fBm2.txt', 0, 128, 75, 1);
%
%   Adpated and updated code from an uncountable number of authors over the years:
%   Jay Holden, Michael J. Richardson, R. C. Schmidt, Charles Coey, Nikita Kuznetsov
%--------------------------------------------------------------------------

%% Load and Check Data
x_data = load(filename);
x = x_data(:,1);

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

    
%% Power Spectrum Density
[P,F]=pwelch(x,triang(2*nfreq),nfreq,2*nfreq,1,'onesided');
F(1)=[];
F(nfreq)=[];
P(1,:)=[];
P(nfreq,:)=[];
Plog=log10(P);
Flog=log10(F);

spec_len = (length(P)+1);
quart_spec = spec_len/4;

switch (percent_spec)
    case 75
        pSCL = spec_len - quart_spec;
        powers = Plog(1:pSCL);
        freqs = Flog(1:pSCL);
    case 50
        pSCL = spec_len - (2*quart_spec);
        powers = Plog(1:pSCL);
        freqs = Flog(1:pSCL);
    case 25
        pSCL = spec_len - (3*quart_spec);
        powers = Plog(1:pSCL);
        freqs = Flog(1:pSCL);
    otherwise
        powers = Plog;
        freqs = Flog;
end

p = polyfit(freqs,powers,1);
slope = p(:,1);
alpha = -1*(slope);

fprintf('%s      Slope: %.2f    Alpha: %.2f\n\n', filename, slope, alpha)  

%% Plot Results
figure;
subplot(2,2,1:2); % Plot data
plot(x, '-b');
xlim([1 length(x)]);

subplot(2,2,3);
plot(F,P,'.-');
ylabel('Power');
xlabel('Frequency');

subplot(2,2,4);
plot(Flog,Plog,'.-'); hold on;
ylabel('Log_1_0(Power)');
xlabel('Log_1_0(Frequency)');
plot(freqs,polyval(p,freqs),'r','linewidth',2); hold off;

%% Output to File
if output == 1;
    fid=fopen('PSD_Stats.csv','a');
    fprintf(fid,'%s,',filename);
    fprintf(fid,'%d,%d,%1.4f\r\n',nfreq,percent_spec,alpha);
    fclose(fid);
end
    
    
    
    
    
   