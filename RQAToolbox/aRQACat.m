function aRQACat(file_name, doPlots, doStatsFile)
% **********************************************************************************************
%   Performs Categorical Auto RQA on data in specified file
%
%   Inputs:
%           filename   : data file; must be txt or csv file with no headers
%           doPlots    : 1=yes, 0=no
%           doStatsFile: 1=yes, 0=no
%
%   Example:
%           aRQACat('Elvis.txt', 1, 1)
%
%   BY: Michael J. Richardson 2009, 2015, UC
%       
%   NOTE: Using some sub-routines developed by Bruce Kay and Michael Richardson
%   at UCONN from 2003 to 2005.
%
%   Contact: michael.richardson@uc.edu
%--------------------------------------------------------------------------------------------

%% Fixed Parameters: Change if needed.
tmin = 1; % do NOT remove line of identity for auto recurrence
minl = 2; % minimum number of points for a line; effects %DET and MeanLine
radius = 0.001; %approx zero.

%% Load data form file; should be two column file
Data_raw = load(file_name);
DataX = Data_raw(:,1);
   
%% Perform RQA stuff
ds  = xRQA_dist(DataX,DataX,1,1);
[td, rs, ~, err_code] = xRQACat_Stats(ds.d,radius,tmin,minl);

if err_code == 0        
    fprintf('%s      REC: %2.3f      DET: %3.3f       MxL: %5.2f\n', file_name, rs.perc_recur, rs.perc_determ, rs.maxl_found);    
else
    fprintf('%s      REC: 0.000        DET: 0.000       MxL: 0.000\n', file_name);  
end;
    
%% Do Plots
if doPlots == 1
    close all;
    scrsz = get(0,'ScreenSize');
    figure('Position',[scrsz(3)/4 scrsz(4)/4 scrsz(3)/3 scrsz(4)/2]);
    h = axes('Position', [0 0 1 1], 'Visible', 'off');

    axes('Position',[.375 .35 .58 .6], 'FontSize', 8)
    image(rot90(td*25)); % invert color
    xlabel('X(i)','Interpreter','none', 'FontSize', 10);
    ylabel('X(j)','Interpreter','none', 'FontSize', 10);
    set(gca,'XTick',[ ]);
    set(gca,'YTick',[ ]);

    axes('Position',[.375 .1 .58 .15], 'FontSize', 8)
    plot(1:length(DataX), DataX, 'bo');
    xlim([1 length(DataX)]);

    axes('Position',[.09 .35 .15 .6], 'FontSize', 8)
    plot(DataX, 1:length(DataX), 'bo');
    ylim([1 length(DataX)]);

    set(gcf, 'CurrentAxes', h);
    str(1) = {['%REC = ', sprintf('%.2f',rs.perc_recur)]};
    text(.05, 0.24, str, 'FontSize', 10, 'Color', 'b');
    str(1) = {['%DET = ', sprintf('%.2f',rs.perc_determ)]};
    text(.05, .20, str, 'FontSize', 10, 'Color', 'b');
    str(1) = {['MAXLINE = ', sprintf('%.0f',rs.maxl_found)]};
    text(.05, .16, str, 'FontSize', 10, 'Color', 'b');
    str(1) = {['MEANLINE = ', sprintf('%.0f',rs.llmnsd(1))]};
    text(.05, .12, str, 'FontSize', 10, 'Color', 'b');
    str(1) = {['ENTROPY = ', sprintf('%.2f',rs.entropy(1))]};
    text(.05, .08, str, 'FontSize', 10, 'Color', 'b');

end;

%% Print Output Stats
if doStatsFile > 0
    fid = fopen('xRQACat_Stats.csv','a');
    fprintf(fid,'%s,',file_name);
    if err_code == 0 
        fprintf(fid,'%f,%f,%f,%f,%f\n',rs.perc_recur, rs.perc_determ, rs.maxl_found, rs.llmnsd(1), rs.entropy(1));
    else
        fprintf(fid,'0.0,0.0,0.0,0.0,0.0\n');
    end
    fclose(fid);
end

return; % end of function



function [ds,err_code] = xRQA_dist(a,b,dim,lag)
%*****************************************************************************************
% FUNCTION [ds,err_code] = xRQA_dist(a,b,dim,lag);
%
% Compute distances between all points of two vectors,
% which can be embedded using time-lags. (The Webber method of PSR for xRQA).
% For recurrence plots.
%
% Input:
%   a, b   -- data vector (1XN or NX1)
%   dim -- embedding dimension (scalar)
%   lag -- time lag in samples (scalar)
%
% Output:
%   ds  -- structure containing the input parameters dim & lag
%          and the distance matrix (ds.d). (Type the structure name
%          to see the list.)
%       Send ds.d to xrqa_do to do the thresholding, etc.
%
% Note: any errors in the input arguments causes the function
% to print an error message to the screen and return dist = 0 (scalar).
%-----------------------------------------------------------------------------------

%% Check input arguments
% Make sure the data vectors are column and of the same length
if size(a,1) < size(a,2)
    a = a';
end

if size(1) < size(2)
    b = b';
end

if size(a,2) ~= 1
    fprintf('Input a data must be a vector, not a matrix\n');
    err_code = 1;
    ds = 0;
    return;
end

if size(b,2) ~= 1
    fprintf('Input b data must be a vector, not a matrix\n');
    err_code = 2;
    ds = 0;
    return;
end

if size(a,1) ~= size(b,1)
    fprintf('Vectors are of different length\n');
    err_code = 3;
    ds = 0;
    return
end

n = size(a,1);

if ~all(size(dim) == [1 1])
    fprintf('Please use a scalar embedding dimension\n');
    err_code = 4;
    ds = 0;
    return;
end

if mod(dim,1) ~= 0
    fprintf('Please use integer embedding dimension\n');
    err_code = 5;
    ds = 0;
    return;
end

if dim <= 0
    fprintf('Please use an embedding dimension >= 1\n');
    err_code = 6;
    ds = 0;
    return;
end

if ~all(size(lag) == [1 1])
    fprintf('Please use a scalar timelag\n');
    err_code = 7;
    ds = 0;
    return;
end

if mod(lag,1) ~= 0
    fprintf('Please use integer timelag\n');
    err_code = 8;
    ds = 0;
    return;
end

if lag <= 0
    fprintf('Please use a time-lag >= 1\n');
    err_code = 9;
    ds = 0;
    return
end

n2 = n - lag*(dim-1);
if n2 <= 0
    fprintf('Not enough data for these embedding parameters\n');
    err_code = 10;
    ds = 0;
    return
end

err_code = 0;

%% Initialize needed matrices
dist = zeros(n2,n2,'single');
v = zeros(n2,dim,'single');


%% Calculate Distance Matrix
if dim > 1
    % Embed if requested
    emb_a=zeros(n2,dim,'single');
    emb_b=zeros(n2,dim,'single');
    for k=1:dim
        emb_a(:,k) = a((1:n2)+lag*(k-1));
        emb_b(:,k) = b((1:n2)+lag*(k-1));
    end
    
    % Compute embedded case
    for i = 1:n2
        for k = 1:dim
            v(:,k) = emb_a(i,k) - emb_b(:,k);
        end
        v = v.^2;
        dist(i,:) = sqrt(sum(v(:,1:dim)'));
    end
else
    
    % Compute non-embedded case
    for i = 1:n2
        dist(i,:) = (abs(a(i) - b))';
    end
    
end

ds = struct('dim',dim','lag',lag,'d',dist);

return;

function [td, rs, mats, err_code] = xRQACat_Stats(d,rad,tmin,minl)
%% Threshold the distance matrix =========================================
td = xRQACat_radius(d,rad,tmin); %Threshold the distance matrix
if td == -1
    fprintf('Error in thresholding ');
    err_code = 1;
    rs = 0;
    mats = 0;
    return
end

%% Find all lines & compute TRENDs ============================================================
[ll,maxl_poss,npts,trend1,trend2] = xRQA_line(td,tmin);
if ll == -1
    fprintf('\nError in line counting\n');
    err_code = 2;
    rs = 0;
    mats = 0;
    return
end

%% Do statistics on the line length distribution ==============================
%fprintf('Histogramming the line lengths...\n');
[lh,llmnsd] = xRQA_histlines(ll,minl); %Compute histogram of line lengths
if lh == -1
    fprintf('Error in creating histograms\n');
    err_code = 3;
    rs = 0;
    mats = 0;
    return
end
%lh is histogram of lines >= minl only
%llmnsd = mean, SD, and N of lines >= minl only

%% ENTROPY of LL distribution ==============
if lh ~= 0
    entropy = xRQA_entropy(double(lh(:,2)),maxl_poss-minl+1); %Compute entropy of line length histogram
    if entropy == -1
        fprintf('Error in computing entropy\n');
        err_code = 4;
        rs = 0;
        mats = 0;
        return
    end
else
    entropy = [0 0];
end

%% Percent recurrence (%RECUR) =============================== 
recur = sum(ll);
perc_rec = 100*recur/npts;

%% DETER, MAXLINE, algorithmic complexity of LL distribution ===============
if lh ~= 0
    perc_determ = 100*sum(double(lh(:,1)).*double(lh(:,2)))/recur; % %DETER; double-precision arithmetic required
    maxl_found = max(lh(:,1)); % Maximum diagonal line found (MAXLINE)
else
    perc_determ = 0;
    maxl_found = 0;
end

%% Output to structures ===============================================================================
rs = struct('tmin',tmin, 'minl',minl, ...
    'perc_recur',perc_rec, 'perc_determ',perc_determ, 'npts',npts, ...
    'entropy',entropy, 'maxl_poss',maxl_poss, 'maxl_found',maxl_found, ...
    'trend1',trend1, 'trend2',trend2, 'llmnsd',llmnsd);

mats = struct('tmin',tmin, 'minl',minl, ...
    'td',td, 'll',ll, 'lh',lh );

err_code = 0;

return


function thrd = xRQACat_radius(dist,rad,tmin)
%*****************************************************************************************
% FUNCTION thrd = xRQACat_radius(dist,rad,tmin);
%% Check input arguments
s=size(dist);
if s(1) ~= s(2)
    fprintf('Distance matrix must be square\n');
    thrd = -1;
    return;
end
if s(1) == 1
    fprintf('Distance matrix has only one element!\n');
    thrd = -1;
    return;
end

s=size(rad);
if ~all(s == [1 1])
    fprintf('Please use a scalar threshold\n');
    thrd = -1;
    return;
end
if rad <= 0
    fprintf('Please use a threshold > 0\n');
    thrd = -1;
    return;
end

s=size(tmin);
if ~all(s == [1 1])
    fprintf('Please use a scalar tmin parameter\n');
    thrd = -1;
    return;
end
if mod(tmin,1) ~= 0
    fprintf('Please use integer tmin\n');
    thrd = -1;
    return;
end
if tmin < 0
    fprintf('Please use a tmin >= 0\n');
    thrd = -1;
    return;
end

%% Nonvectorized algorithm: actually faster almost all of the time for this single-precision algorithm
ldist = length(dist);
thrd = zeros(ldist,ldist,'int8');
for i = 1:ldist
    for j = 1:ldist
        if dist(i,j) <= rad
            thrd(i,j) = 1;
        end
    end
end

%% Set distances within tmin as non-recurrent
if tmin ~= 0
    lt = length(thrd);
    for i = 0:tmin-1
        for j = 1:lt-i
            thrd(j,j+i) = 0;
            thrd(j+i,j) = 0;
        end
    end
end
return


function [ll,maxl_poss,npts,trend1,trend2] = xRQA_line(thrd,tmin)
%*****************************************************************************************
% FUNCTION [ll,maxl_poss,npts,trend1,trend2] = xRQA_line(thrd,tmin);
%
% Single-precision version of xrqa_line1.m; works only with Matlab R14+
%
% Find all diagonal lines (parallel to the main diagonal) and their lengths in a thresholded distance matrix.
% Do so for ENTIRE recurrence plots (compare xrqa_line2_single). Also computes TRENDs.
%
% Input:
%   thrd -- square thresholded distance matrix having 1's and 0's at least on the upper diagonal
%   tmin -- minimum recurrence time (retain main diagonal if == 0); sets the Theiler window
%   plotopt -- plotting option (TREND: recurrence as function of distance from main diagonal) (1=yes, other = no)
%
% Output:
%   ll -- vector containing the lengths of all the lines found, even those >= 1 and < minl
%   max_poss -- maximum possible line length
%   npts -- total number of possible points in the recurrence plot, excluding the Theiler window
%   trend1 -- Webber's Trend measure for one triangle of the threshold matrix
%   trend2 -- " for the other triangle of the threshold matrix
%
% Note: any errors in the input argument causes the function
% to print an error message and return all outputs = 0.
%-----------------------------------------------------------------------------------

%% Check input matrix
s=size(thrd);
if s(1) ~= s(2)
    fprintf('Thresholded distance matrix must be square\n');
    ll = 0;
    trend1 = 0;
    trend2 = 0;
    maxl_poss = 0;
    npts = 0;
    return;
end

if s(1) == 1
    fprintf('Thresholded distance matrix has only one element!\n');
    ll = 0;
    trend1 = 0;
    trend2 = 0;
    maxl_poss = 0;
    npts = 0;
    return;
end

mx = max(max(thrd));
if mx > 1
    fprintf('Please use a thresholded distance matrix\n');
    ll = 0;
    trend1 = 0;
    trend2 = 0;
    maxl_poss = 0;
    npts = 0;
    return;
end
n = s(1);

s=size(tmin);
if ~all(s == [1 1])
    fprintf('Please use a scalar tmin parameter\n');
    ll = 0;
    trend1 = 0;
    trend2 = 0;
    maxl_poss = 0;
    npts = 0;
    return;
end

if mod(tmin,1) ~= 0
    fprintf('Please use integer tmin\n');
    ll = 0;
    trend1 = 0;
    trend2 = 0;
    maxl_poss = 0;
    npts = 0;
    return;
end

if tmin < 0
    fprintf('Please use a tmin >= 0\n');
    ll = 0;
    trend1 = 0;
    trend2 = 0;
    maxl_poss = 0;
    npts = 0;
    return;
end

% Compute; very slow if lots of long lines, but at least it works correctly
% The following "brute force" method seems inefficient, but it's the fastest thing I've tried.
% See count_lines.m & xrqa_lines1_single_b.m for another algorithm.

possnumll = ceil((n^2)/2);
recur = zeros(2*n-1,2,'single');
%recur = zeros(2*n-1,2);

ll = zeros(possnumll,1,'int16'); %line lengths found; this seems to be a reasonable guess as to the max possible # of found lines,
%but we may need to change this.
nlines = 0;
for i = 1:2*n-1 %Does all diagonals, for cross-recurrence case
    d = diag(thrd,i-n);
    ld = length(d);
    recur(i,1) = ld;
    recur(i,2) = 0; %Total recurrence along ith diagonal

    j = 1;
    while j <= ld

        if d(j) == 1
            nlines = nlines + 1;
            ll(nlines) = 1;
            recur(i,2) = recur(i,2)+1;

            k=j+1;
            while k <= ld

                if d(k) == 1
                    ll(nlines) = ll(nlines) + 1;
                    recur(i,2) = recur(i,2)+1;
                    k = k + 1;
                else
                    break
                end

            end
            j=k+1;

        else
            j=j+1;
        end
    end
end

% Shrink ll to actual # of lines found
if nlines > possnumll
    fprintf('***Please increase possnumll in xrqa_line1.m!***\n');
    fprintf('Initialized line length vector size exceeded, performance greatly reduced!\n');
end

ll = ll(1:nlines);

% Compute trend
% Acc. to Weber, = slope of line-of-best fit through percentrecurrence
% as a function of displacement from main diagonal, excluding the last
% ten percent of the range, in units of percent recurrence (locally
% measured) per 1000 points ("points"=distance from main diagonal).
% Do it for the upper & lower diagonal matrices and average

mid = find(recur(:,1) == max(recur(:,1)));

recur(:,2) = recur(:,2)./recur(:,1);

% lower:
first = tmin;
last = n-1;
indx = (first:last)';

p = polyfit(indx,100*recur(mid-tmin:-1:1,2),1);
trend1 = 1000*p(1);
% upper:
first = mid+tmin;
last = 2*n-1;

p = polyfit(indx,100*recur(first:last,2),1);
trend2 = 1000*p(1);

lmd = length(thrd);
maxl_poss = lmd - tmin;
if tmin == 0
    npts = lmd^2;
else %Total number of points in the recurrence plot, excluding the diagonals taken up by tmin
    npts = lmd^2 - lmd - 2*lmd*(tmin-1) + tmin*(tmin-1);
end

return

function [linehist,linestats] = xRQA_histlines(llengths, minl)
%*****************************************************************************************
% FUNCTION [linehist,linestats] = xRQA_histlines(llengths, minl);
%
% Inputs:
%   llengths -- vector of line lengths
%   minl -- minimum line length to consider
%	 plotopt  -- 1 = do plots of line length histogram
%			subplot 1: linear scaling on both axes
%			subplot 2: log scaling on both axes
%
% Output:
%   linehist -- NX2 matrix for doing histogram of line lengths a la Webber
%               (N=# of distinct line lengths found)
%               column 1 = line length
%               column 2 = frequency
%               Contains only the line lengths >= minl
%   linestats -- mean, SD, and N of line lengths >= minl
%
% Note: any errors in the input arguments causes the function
% to print an error and return linehist = -1 (scalar) and linestats = [-1 -1 -1].
%-----------------------------------------------------------------------------------

%% Check input arguments ===================================================
s = size(llengths);
if s(1) < s(2)
    llengths=llengths';
end
s=size(llengths);

if s(2) ~= 1
    fprintf('Input data must be a vector, not a matrix\n');
    linehist = single(-1);
    linestats = [-1 -1 -1];
    return
end

s=size(minl);
if ~all(s == [1 1])
    fprintf('Please use a scalar for min line length\n');
    linehist = single(-1);
    linestats = [-1 -1 -1];
    return;
end

if mod(minl,1) ~= 0
    fprintf('Please use integer min line length\n');
    linehist = single(-1);
    linestats = [-1 -1 -1];
    return;
end

if minl <= 0
    minl = 1;
end

%% Begin computations ===========================================================

ll=length(llengths);
if min(llengths) == max(llengths)
    %  fprintf('Only one line length found\n');
    if llengths(1) >= minl
        linehist = single([llengths(1) ll]);
        linestats = [double(llengths(1)) 0 1];
    else
        linehist = single(0);
        linestats = [0 0 0];
    end
    return;
end
llengths = sort(llengths);

%% Find first observed line length >= minimum line length
x2 = [];
for i = 1:ll
    if llengths(i) >= minl
        x2 = llengths(i);
        i_first = i;
        break;
    end
end

%% Perform Calculations
if ~isempty(x2)
    %Do stats on the line lengths that are greater than minl
    ll_mean = mean(double(llengths(i_first:end)));
    ll_std = std(double(llengths(i_first:end)));
    ll_n = length(llengths(i_first:end));
    linestats = [ll_mean ll_std ll_n];

    %Do the histogram
    max_num_bins = llengths(end)-llengths(1)+1;
    x3 = zeros(max_num_bins,1,'single');
    lh2 = zeros(max_num_bins,1,'single');
    x3(1) = x2;
    nbins = 1;
    lh2(1) = 1;

    % Proceed from there
    for j = i+1:ll

        if llengths(j) == x3(nbins)
            lh2(nbins) = lh2(nbins)+1;
        else
            nbins = nbins + 1;
            x3(nbins)=llengths(j);
            lh2(nbins)=1;
        end

    end

    x3 = x3(1:nbins);
    lh2 = lh2(1:nbins);
    linehist =  [x3 lh2];

else
    linehist = single(0);
    linestats = [0 0 0];
end

return

function entropy = xRQA_entropy(distr,nstates)
%*****************************************************************************************
% FUNCTION entropy = xrqa_entropy(distr,nstates);
%
% Compute Shannon entropy of a distribution
%
% Input
%   distr -- distr vector, either column or row.  All elements
%            are considered to be part of the distribution.
%   nstates -- number of possible states
%
% Output
%   entropy -- 1X2 vector containing shannon entropy &
%              max entropy possible given nstates - shannon entropy =
%              information according to Layzer 1988
%
% Example:
% my_entropy = xrqa_entropy(my_distr, 234);
%
% where my_distr is your distribution vector. 
%
% Note: any errors in the input arguments causes the function
% to print an error and return entropy = -1 (scalar).
%-----------------------------------------------------------------------------------

%% Check input argument
s=size(distr);
if s(1) < s(2)
    distr=distr';
end

s=size(distr);
if s(2) ~= 1
    fprintf('Input data must be a vector, not a matrix\n');
    entropy = -1;
    return
end

s=size(nstates);
if ~all(s == [1 1])
    fprintf('Please use a scalar for number of states\n');
    entropy = -1;
    return;
end

if mod(nstates,1) ~= 0
    fprintf('Please enter integer number of states\n');
    entropy = -1;
    return;
end

if nstates <= 0
    fprintf('Number of states must be > 0\n');
    entropy = -1;
    return;
end


%% Compute
distr2 = distr/sum(distr);
ld = length(distr2);
sumdist = 0;
for i = 1:ld
    if distr2(i) ~= 0
        sumdist = sumdist - distr2(i)*log2(distr2(i));
    end
end
maxen = log2(nstates);
entropy = [sumdist maxen-sumdist];

return

