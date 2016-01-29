function NS_PreprocessStimData_streaming(rawDataDir,WritePathBase,fileNos,varargin)
% NS_PreprocessStimData_streaming() organizes raw electrical stimulation 
% data according to stimulation patterns. A folder for each pattern is 
% created, and data surrounding current pulses for all repetitions are 
% saved in a new file for each movie chunk / amplitude
%       INPUTS: 
%           rawDataDir: points to the raw data directory
%           WritePathBase: Points to the directory of the analysis output
%           fileNos: vector containing the dataruns to preproecss.  Appends these numbers to 'data0--'
%       optional:
%           system: 'stim512' or 'stim64'. Default is 'stim512'
%           traceLength: length of trace to save after each pulse (in samples)
% At 20 samples/millisecond, 100 samples = 5 milliseconds. Default is 100
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Set up defaults for optional parameters
system = 'stim512';  %stim512 or stim64
traceLength = 100; % length of trace to save after each pulse (in samples)
% At 20 samples/millisecond, 100 samples = 5 milliseconds

nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn','Unexpected number of arguments');
    throw(err);
end

% Read the optional input arguments
for j=1:(nbin/2)
    if ~ischar(varargin{j*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{j*2-1})
        case 'system'
            system = varargin{j*2};
        case 'tracelength'
            system = varargin{j*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

if ~strcmp(rawDataDir(end),filesep)
    rawDataDir = [rawDataDir filesep];
end
disp(rawDataDir); 

% Points to the directory of the output.
if ~strcmp(WritePathBase(end),filesep)
    WritePathBase = [WritePathBase filesep];
end

% Choose system
switch system
    case 'stim512'
        numElectrodes = 512; 
    case 'stim64'
        numElectrodes = 64; 
end

% Sets some parameter values depending on what system is in use.
NS_GlobalConstants = NS_GenerateGlobalConstants(numElectrodes);

for i = fileNos
    % Set the number of zeros in the file name.
    if i < 10
        FileName = ['00' num2str(i)];
    elseif i < 100
        FileName = ['0' num2str(i)];
    else
        FileName = num2str(i);
    end
    
    rawDataPath = [rawDataDir 'data' FileName];
    WritePath = [WritePathBase 'data' FileName];
    disp(['analyzing ' rawDataPath]);
    % Determine if WritePath exists, if not then make it.
    if ~exist(WritePath, 'file')
        mkdir(WritePath)
        disp(['Created new directory: ' WritePath])
    end

    if ~exist([WritePath filesep 'pattern_files'], 'file')
        mkdir([WritePath filesep 'pattern_files'])
        disp(['Created new directory: ' WritePath filesep 'pattern_files'])
    end

    if ~exist([WritePath filesep 'status_files'], 'file')
        mkdir([WritePath filesep 'status_files'])
        disp(['Created new directory: ' WritePath filesep 'status_files'])
    end
    
    %%%%%% standard %%%%%%
    %     All repetitions of a particular pattern, at a particular amplitude, inside a particular movie chunk.
    %     For example: Not appropriate for moving bar.
    %     Generates the "p m" files (the pattern-movie files), the
    %     stimulus_files, and the status_files
    
    NS_SaveRespForMovieAllPatternsAllChannels_stream(rawDataPath,WritePath,...
        NS_GlobalConstants,traceLength);
    
end
