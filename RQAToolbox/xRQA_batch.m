function xRQA_batch()
% **********************************************************************************************
%   function xRQA_batch()
%
%   Performs a batch Cross Recurrence Quantification Analysis
%
%   M.J. Richardson 7/2004, updated 2009, 2015
%
%-------------------------------------------------------------------------------------
fprintf('Processes xRqa on txt files\n');

% Get User Input for Parameters and File Settings
dirtxt_prim = input('Enter primary .txt or .csv file pattern (wildcards, but NO EXTENSION): ','s');
ext1 = input('Enter file extension (e.g. .txt or .csv): ','s');
norm = input('Normalize data (0=none, 1=unit interval, 2=zscore, 3=center): ');
Edim = input('Embedding dimen: ');
Tlag = input('Time lag: ');
rescale = input('Rescale option (1=mean 2=max):   ');
radius = input('Radius: ');
doStatsFile = input('Output stats to file (1=yes, 0=no):' );

if isempty(strfind(dirtxt_prim,'.'))
    dirtxt_prim = strcat(dirtxt_prim,ext1);
end
dir_in_prim=dir(dirtxt_prim);
nfiles_prim = length(dir_in_prim);

if nfiles_prim==0
    fprintf('No %s files found\n',ext1);
    return   
end

fns_prim = sortfiles(dir_in_prim);

fprintf('\n\n%d file(s) to process... ... ...\n',nfiles_prim);

for i = 1:nfiles_prim
    file_name = fns_prim(i,:);
    file_name = deblank(file_name);
    
    xRQA(file_name, norm, Edim, Tlag, rescale, radius, 0, doStatsFile);
    
end

fprintf('\nDone!\n');
return;


%**************************************************************************
function sorted_names = sortfiles(direc);
%function sorted_names = sortfiles(direc);
%Sort a list of file names
%M.J. Richardson 2004
ld = length(direc);
[fn{1:ld,1}] = deal(direc.name);
sorted_names = char(fn);
sorted_names = sortrows(sorted_names);
return