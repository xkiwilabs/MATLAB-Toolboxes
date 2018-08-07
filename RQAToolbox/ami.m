function amiF = ami(filename,min_lag,max_lag)

%% *************************************************************************************%	
%	ami_do:
%
%	Syntax:
%		amiF = ami(file_name,min_lag,max_lag)
% 
%	Example:
%		amiF = ami('datafile.txt',1,100);
%
%	BY: C. A. Coey and Michael Richardson (03/11/2013)
%
%   Based on example code  by Alexandros Leontitsis 
%***************************************************************************************%

%% Load Data
x_data = load(filename);
x = x_data(:,1);
len = length(x);


%% Create Lag Vector
if max_lag <= (len/2 - 1)
    if min_lag < max_lag
        i = min_lag;
        j = 1;
        while i <= max_lag
            lag(j) = i;
            j = j+1;
            i = i+1;
        end;
    else
        fprintf('error - maximum lag not greater than minimum lag\n');
        fprintf('default lag vector used (0 - 50)\n');
        lag = (0:50);
    end;

else
    fprintf('error - maximum lag exceeds recommendation\n');
    fprintf('maximum lag set to n/2-1\n');
    lag = (0:(len/2 - 1));
end;

%% Compute Average Mutual Information
x=x-min(x);
x=x/max(x);
h = waitbar(0, 'Processing AMI...');
for i=1:length(lag)
   % Define the number of bins
   k=floor(1+log2(len-lag(i))+0.5);
   
   % If the time series has no variance then the MAI is 0
   if var(x,1)==0
      v(i)=0;
   else
      v(i)=0;
      for k1=1:k
         for k2=1:k
             
            ppp=find((k1-1)/k<x(1:len-lag(i)) & x(1:len-lag(i))<=k1/k ...
               & (k2-1)/k<x(1+lag(i):len) & x(1+lag(i):len)<=k2/k);
            ppp=length(ppp);
            px1=find((k1-1)/k<x(1:len-lag(i)) & x(1:len-lag(i))<=k1/k);
            px2=find((k2-1)/k<x(1+lag(i):len) & x(1+lag(i):len)<=k2/k);
            if ppp>0
               ppp=ppp/(len-lag(i));
               px1=length(px1)/(len-lag(i));
               px2=length(px2)/(len-lag(i));
               v(i)=v(i)+ppp*log2(ppp/px1/px2);
            end
         end
      end
   end
   waitbar(i/length(lag), h);
end
delete(h);
amiF = [lag' v'];

%% Plot AMI Function
figure;
plot(amiF(:,1), amiF(:,2));
xlabel('TLag');
ylabel('AMI');

return;
