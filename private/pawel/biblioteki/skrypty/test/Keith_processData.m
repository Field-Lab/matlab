%{
File: 101014 - Data processing V3\processData.m
Version: V4 - 16/11/10 15:13

This script removes the artifacts from the raw data and creates a new bin 
file with the data in the original format, so that the data with artifacts 
substracted can be processed in vision.
The first 10 pulses of data are skipped because the artifact takes a few 
pulses to converge to a stable shape and adding these pulses of data could
produce false positives.

Useful variables:
    - dataFolder: where the original data file can be found

%}

function processData(dataFolder,varargin)
% This function removes the artifacts from raw MEA recordings. 
% It processes the data stored in Vision-compressed .bin files that can be 
% found in dataFolder, creates a new subfolder 'Processed' and writes the
% new .bin file adequately compressed there.
% In order to remove the artifacts from the data file, it needs to have
% access to the artifact shapes, so it first tries to load them if they
% have been computed before, if they haven't they are computed and stored
% in a 'avTrace[nnn].mat' file, where [nnn] is the reference number for the
% data being processed.
%
% Parameters:
%   - dataFolder: directory of the raw data
%
% dataFolder can be followed by parameter/value pairs to specify additional
% properties. Acceptable parameters are:
%
%   - newDataFolder: if specified, the script will create a new this data
%   folder if it does not exist and will write the new bin file there
%
%   - pulseDuration: if not specified, it tries to read it from the raw
%   data file header comment. If not specified and it is not in the
%   comment, an error will occur.
%
%   - useTrigger: boolean indicating if the trigger recording should be
%   used to detect when pulses occured. If 0, then it is assumed there is
%   exactly one pulse per half second, if 1 then in every half second of
%   data an artifact is removed only if a trigger pulse is found. By
%   default, useTrigger is set as false.
%
%   - useSmoothing: boolean indicating if you wish to smooth the beginning
%   and ending of the pulse or not to remove potential leftover artifacts
%
%   - visionPath: path to the .jar vision archive. Set acceptable default
%   for your settings in "Reading the input Arguments" cell.
%
%   - upSampling: by default 1, it should be increased for finer
%   positionning of the artifact. A value of n leads to positionning of the
%   artifact to the 1/n sample
%
%   - nTrials: number of trials you wish to use for computation of the
%   artifact, by default 200.
%
%   - Fexperiment: frequency of the light pulse, by default set to 2Hz
%
%   - nSkipAverage: number of pulses skipped at the beginning of the data for
%   avTrace computation
%
%   - nSkip: number of pulses skipped at the beginning of the artifact
%   removal
%
% Example:
%   processData('C:\Data\data001'); 
%   processData('C:\Data\data001','pulseDuration',4); 
%   processData('C:\Data\data001','useTrigger',false); 
%   processData('C:\Data\data001','pulseDuration',4,...
%                                 'useTrigger',false); 
%
% Returns: []

%% Reading the input arguments 

% Making sure dataFolder ends by '/'
if dataFolder(end:end)~='/'
    dataFolder = [dataFolder '/'];
end

% Getting the data reference
strPos = strfind(dataFolder,'data');
strPos = strPos(end);
dataReference = dataFolder((strPos+4):(strPos+6));

% Setting default value for the optional input arguments
pulseDuration = -1;
useTrigger = 0;
useSmoothing = 0;
visionPath = '/Users/keithmathieson/analysis/vision8/Vision.jar';
upSampling = 1;
nTrials = 200;
Fexperiment = 2;
nSkipAverage = 10;
newDataFolder = [dataFolder 'Processed/'];
nSkip = 20;

% Checking the optional parameters
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'newdatafolder'
            newDataFolder = varargin{kk*2};
            if newDataFolder(end:end)~='/'
                 newDataFolder = [newDataFolder '/'];
            end
        case 'pulseduration'
            pulseDuration = varargin{kk*2};
        case 'usetrigger'
            useTrigger = varargin{kk*2};
        case 'usesmoothing'
            useSmoothing = varargin{kk*2};
        case 'visionpath'
            visionPath = varargin{kk*2};
        case 'upsampling'
            upSampling = varargin{kk*2};
        case 'ntrials'
            nTrials = varargin{kk*2};
        case 'fexperiment'
            Fexperiment = varargin{kk*2};
        case 'nskipaverage'
            nSkipAverage = varargin{kk*2};
        case 'nSkip'
            nSkip = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

clear nbin kk

%% Linking to the bin files

% Linking to the jar archive
javaaddpath(visionPath);

% Creating the new directory for the processed file if it doesn't exist
% already
mkdir(newDataFolder);

% Path to the new file
full_path_new=[newDataFolder dataReference '.bin'];

% Linking to the files 
rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(dataFolder);
header = rawFile.getHeader();
newFile = edu.ucsc.neurobiology.vision.io.ModifyRawDataFile(full_path_new,header);

%% Initialization of the artifact
ah = artifact(dataFolder,...
                'pulseDuration',pulseDuration,...
                'useTrigger',useTrigger,...
                'useSmoothing',useSmoothing,...
                'visionPath',visionPath,...
                'upSampling',upSampling,...
                'nTrials',nTrials,...
                'Fexperiment',Fexperiment,...
                'nSkip',nSkipAverage);

clear pulseDuration useTrigger useSmoothing visionPath upSampling nTrials Fexperiment nSkipAverage

%% Computing the artifact

try
    ah.load()
catch loadError
    ah.compute()
end

%% Removing the artifact

% Waitbar and counter
startSample = nSkip*ah.nSamplesPulse; 

nn = nSkip; 
nMax = (header.getNumberOfSamples()/ah.nSamplesPulse);
wh = waitbar(0,'Processing and writing the data');

% Processing the raw data
while (true)
    try
        % Getting all the data including trigger
        rawData = double(rawFile.getData(startSample,ah.nSamplesPulse));        
    catch javaEndOfFileError
        break;
    end

    % Removing the artifacts
    procData = ah.processData(rawData);

    % Writing the data
    newFile.writeRawData(startSample,procData);

    % Incrementing the sample counter
    startSample = startSample + ah.nSamplesPulse;

    % Waitbar update
    nn = nn+1;
    waitbar(nn/nMax,wh);
end

close(wh);


