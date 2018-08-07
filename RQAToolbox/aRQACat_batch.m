function aRQACat_batch()
% **********************************************************************************************
%   function aRQACat_batch()
%
%   Performs a batch Categorical Auto Recurrence Quantification Analysis
%
%   M.J. Richardson 2009, updated 2015
%
%-------------------------------------------------------------------------------------
fprintf('Processes xRqaCat on txt or csv files\n');

dirtxt_prim = input('Enter primary .txt file pattern (wildcards, but NO EXTENSION): ','s');
ext1 = input('Enter file extension (e.g. .txt or .csv): ','s');
doStatsFile = input('Output stats to file (1=yes, 0=no):' );

if isempty(findstr(dirtxt_prim,'.'))
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
    aRQACat(file_name, 0, doStatsFile);
end

fprintf('\nDone\n');
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