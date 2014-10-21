%{
File: analyzeSpikes.m
Version: v5.01 - 04/22/2013

This function computes the sliding PSTH for different neurons. 
The spikes are detected and sorted by the Vision software.

Useful parameters: 
    - dataFolder
    - imageFormat: by default, eps, could be any format supported by Matlab
    - visionPath: change the default value if you don't want to specify it 
    every time
    - rawDataFolder: directory where the original data can be found

Make sure there is only one neuron file in the dataFolder folder, otherwise 
weird things could happen.

%}

function analyzeSpikesGratings(dataFolder,logFilePath,varargin)
% This function computes the PTSH of all the neurons found in a MEA
% experiment recording. The function takes as input a folder containing the
% Vision neuron file and creates a new subfolder in which it saves the
% histograms in the format specified.
% The function also saves the information relative to each neuron in a
% Matlab structure named neuronData[nnn].mat, where [nnn] is the dta
% reference. 
% This information is saved in dataFolder.
% 
% Parameters:
%   - dataFolder: the directory containing the neuron file to be analyzed
%   - logFilePath: the path to the experiment log file. If specified, instead
%   of plotting all the experiments in separate figures, the program will
%   try to combine them in relevant ways (same power varying pulse width,
%   same pulse width varying power).
%
% dataFolder can be followed by parameter/value pairs to specify additional
% properties. Acceptable parameters are:
%   - visionPath: path to the .jar vision archive. Set acceptable default
%   for your settings in "Reading the input Arguments" cell.
%   - rawDataFolder: if specified it should point towards the unprocessed
%   data. 
%   - binsize: if specified the script does not choose the bin size anymore
%   but uses what the user specified. This bin size should be specified in
%   number of samples, not milliseconds or seconds. Good bin size: 5ms =
%   100 samples.
%
% Returns: []

%% Reading the input arguments

% Making sure dataFolder ends by '\' or '/', whichever is right
if dataFolder(end:end)~=filesep
    dataFolder = [dataFolder filesep];
end

% Setting default values for the optional input arguments
if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    visionPath = 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end
rawDataFolder = '';
overrideBinSize = true;
binSize = 100;
outputPath = [dataFolder 'statistics'];
neuronDataFolder = dataFolder;
ttlsPerPulse = 1;

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
        case 'visionpath'
            visionPath = varargin{kk*2};
        case 'rawdatafolder'
            rawDataFolder = varargin{kk*2};
        case 'neurondatafolder'
            neuronDataFolder = varargin{kk*2};
        case 'binsize'
            overrideBinSize = true;
            binSize = varargin{kk*2};
        case 'ttlsperpulse'
            ttlsPerPulse = varargin{kk*2};
        case 'outputpath'
            outputPath = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

if neuronDataFolder(end:end)~=filesep
    neuronDataFolder = [neuronDataFolder filesep];
end

clear nbin

%% Linking to the bin and neuron file

% Finding the path to the neuron file in the data folder
contentsNeuronDataFolder = dir(neuronDataFolder);

for kk=1:length(contentsNeuronDataFolder)
    isNeuronFile = strfind(contentsNeuronDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsNeuronDataFolder(kk).name,'.neurons-raw');
    if isNeuronFile
        if isNeuronRawFile
        else
            neuron_path =  [neuronDataFolder contentsNeuronDataFolder(kk).name];
            break;
        end
    end
end

if ~exist('edu/ucsc/neurobiology/vision/io/NeuronFile','class')
    javaaddpath(visionPath);
end

neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);
binFile = edu.ucsc.neurobiology.vision.io.RawDataFile(dataFolder);       
if ~isempty(rawDataFolder)
    rawBinFile = edu.ucsc.neurobiology.vision.io.RawDataFile(dataFolder);
end
header = binFile.getHeader();

clear isNeuronFile

%% Reading the log file and loading the pulse times, if necessary

M = readLogFileGratings(logFilePath);
FExperiments = M.frequency;
nExperiments = length(FExperiments);
stimType = M.stimulus_type;

Fs = header.getSamplingFrequency();     % Sampling frequency

allPulseTimes = loadAllPulsesTimeFromNeuronFile(neuronFile, nExperiments, Fs);

%% Basic file information

% Creating the folder where the neuron data will be saved
if ~exist(outputPath,'dir')
    mkdir(outputPath);
end

% Closing the bin files
binFile.close;      
if ~isempty(rawDataFolder)
    rawBinFile.close;
end

%% Processing all the neurons

neuronList = neuronFile.getIDList();

for kk=1:length(neuronList)
    % Creating the neuron information instance
    neuronData = neuronInformation('neuronID',neuronList(kk),...
                                   'rawSpikeTimes',neuronFile.getSpikeTimes(neuronList(kk)),...
                                   'pulseTimeInformation',allPulseTimes,...
                                   'Fexperiment',FExperiments,...
                                   'stimulationType',stimType,...
                                   'Fs',Fs,...
                                   'ttlsPerPulse',ttlsPerPulse,...
                                   'electrode',neuronFile.getElectrode(neuronList(kk)));
    
    % Processing the data
    if overrideBinSize
        neuronData.computeFiringInformation(binSize);
    else
        neuronData.computeFiringInformation();
    end
    
    % Saving the neuron information
    neuronData.save(outputPath);
end

end % analyzeSpikes

function allPTStruct = loadAllPulsesTimeFromNeuronFile(neuronFile, nExperiments, Fs)
% This function loads all the TTL pulses information into a matlab
% structure.

% More than two seconds between TTL pulses: assume next experiment.
thresholdTTLDiff = 2*Fs; 

rawPT = double(neuronFile.getTTLTimes());
timeDiffTTL = diff(rawPT);
posStartExp = [0; find(timeDiffTTL>thresholdTTLDiff)] + 1;
% Checking we're seeing the right number of experiments
if length(posStartExp)~=nExperiments
    err = MException('MATLAB:TTLPulses:InconsistentTTLs',...
                sprintf('%s: %d stimuli expected, %d seen',...
                'TTL pulses information does not appear to be consistent with the number of stimuli',...
                nExperiments,length(posStartExp)));
            throw(err);
end

allPTStruct = repmat(struct('experimentID',[],'startSample',[],...
    'TTLTimes',[],'nSamplesExperiment',[]),1,nExperiments);
for kk=1:nExperiments
    allPTStruct(kk).experimentID = kk;
    allPTStruct(kk).startSample = rawPT(posStartExp(kk)) - 100;
    if kk<nExperiments
        allPTStruct(kk).TTLTimes = rawPT(posStartExp(kk):(posStartExp(kk+1)-1)) ...
            - allPTStruct(kk).startSample;
    else
        allPTStruct(kk).TTLTimes = rawPT(posStartExp(kk):end) ...
            - allPTStruct(kk).startSample;
    end
    % Here: assume there is at least 1 second between experiments.
    allPTStruct(kk).nSamplesExperiment = allPTStruct(kk).TTLTimes(end) + Fs;
end

end