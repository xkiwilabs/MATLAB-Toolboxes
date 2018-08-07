function synchroBatch()
%**************************************************************************
%   SYNCHROBATCH calculates and outputs all of the above IRP or DRP measures 
%   on a batch of selected continuous two column .txt or .csv data files in
%   the current directory.
%
%   Syntax:
%   	synchroBatch();
%
%   Michael J. Richardson 2009
%   Last Updated 2013.
%--------------------------------------------------------------------------

%% Get User Input for Parameter settings
dirtxt = input('Enter file pattern (wildcards (i.e., myFiles*), but NO EXTENSION): ','s');
ext = input('Enter file extension (e.g. .txt or .csv): ','s');
samplerate = input('Sample rate (Hz): ');
phaseType = input('relative phase type (1 = IRP, 2 = DRP): ');
phaseMode = input('relative phase mode (1 = inphase, 2 = antiphase): ');
plotFlag = input('display plots (0 = no,  1 = yes): ');
outFlag = input('output stats file (0 = none,  1 = stats + 0-180 dist, 2 = stats + 0-360 dist): ');

%% Find and sort files defined by dirtxt.ext in current directory
if isempty(findstr(dirtxt,'.'))
    dirtxt = strcat(dirtxt,ext);
end
dir_in=dir(dirtxt);
nfiles = length(dir_in);
if nfiles==0
    fprintf('No %s files found\n',ext);
    return   
end;
ld = length(dir_in);
[fn{1:ld,1}] = deal(dir_in.name);
fns = char(fn);
fns = sortrows(fns);

fprintf('\n\n%d file(s) to process... ... ...\n',nfiles);

%% Process loop for all files defined by  dirtxt.ext in current directory 
for i = 1:nfiles 
    filename = fns(i,:);
    filename = deblank(filename);
    fprintf('processing...  %s\n',filename); 
    if phaseType == 2
        synchroDRP(filename, samplerate, phaseMode, plotFlag, outFlag);
    else
        synchroIRP(filename, samplerate, phaseMode, plotFlag, outFlag);
    end
    fprintf('\n');
end

%% End of Function
return;