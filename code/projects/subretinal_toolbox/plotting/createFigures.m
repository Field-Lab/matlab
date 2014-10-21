%{
File: createFigures.m
Version: v4.05 - 03/19/2013

This function plots the PSTH for all specified neurons. If a neuron list is
not specified it plots the PSTH for all neurons.

Useful parameters: 
    - dataFolder
    - imageFormat: by default, eps, could be any format supported by Matlab
    - logFilePath: path to the experiment log file 
%}

function createFigures(statsFolder,figuresFolder,varargin)
% This function plots the sliding PTSH of all the neurons found in a MEA
% experiment recording.
% 
% Parameters:
%   - statsFolder: the directory containing the neuron file to be analyzed
%
% statsFolder can be followed by parameter/value pairs to specify additional
% properties. Acceptable parameters are:
%   - imageFormat: a string specifying which image format should be used to
%   save the images. By default images are saved as eps file, other formats
%   could be 'jpg', 'bmp', 'png',...
%   - logFile: the path to the experiment log file. If specified, instead
%   of plotting all the experiments in separate figures, the program will
%   try to combine them in relevant ways (same power varying pulse width,
%   same pulse width varying power).
%   - neuronList: if specified, the figures are plotted only for the
%   neurons in this list. This list should take the form of a vector of
%   neuron IDs.
%   - combineStim: by default true, if true the script will try to match
%   parameters together when creating the plots, otherwise it will create
%   one plot per stimulus
%
% Returns: []

%% Reading the input arguments

% Making sure dataFolder ends by '\' or '/', whichever is right
if statsFolder(end:end)~=filesep
    statsFolder = [statsFolder filesep];
end

% Setting default values for the optional input arguments
imageFormat = 'jpg';
logFilePath = '';
useNeuronList = false;
powerToIrradianceFactor = 1;
combineStim = true;

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
        case 'imageformat'
            imageFormat = varargin{kk*2};
        case 'logfile'
            logFilePath = varargin{kk*2};
        case 'neuronlist'
            useNeuronList = true;
            neuronList = varargin{kk*2};
        case 'powertoirradiancefactor'
            powerToIrradianceFactor = varargin{kk*2};
        case 'combinestim'
            combineStim = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

clear nbin

%% Creating the figures folder if it doesn't already exist
if ~exist(figuresFolder,'dir')
    mkdir(figuresFolder);
end

%% Reading the log file, if it was specified

if ~isempty(logFilePath)
    [M,paramNames] = readLogFile(logFilePath);
    
    % Cleaning up M: removing parameters we don't care about, such as
    % dataset and stimulus number, or number of pulses
    irrelevantParams = {'Experiment duration','Stimulus','Start Time','Number of pulses','Image projected'};
    for kk=1:length(irrelevantParams)
        pos = find(ismember(paramNames, irrelevantParams{kk})==1);
        M(:,pos) = [];
        paramNames(pos) = [];
    end
    
    M(:,ismember(paramNames,'Power')) = M(:,ismember...
        (paramNames,'Power'))*powerToIrradianceFactor;
end

clear irrelevantParams kk pos

%% Getting the list of neurons to read

contentsStatsFolder = dir(statsFolder);
neuronNames = struct('name','');
nNeurons = 0;

if useNeuronList
    for kk=1:length(neuronList)
        nNeurons = nNeurons + 1;
        neuronNames(nNeurons).name = ['neuron' num2str(neuronList(kk)) '.mat'];
    end
else
    for kk=1:length(contentsStatsFolder)
        if strfind(contentsStatsFolder(kk).name,'.mat')
            nNeurons = nNeurons + 1;
            neuronNames(nNeurons).name = contentsStatsFolder(kk).name;
        end
    end
end

%% Processing all the neurons

for kk=1:nNeurons
    % Loading
    load([statsFolder neuronNames(kk).name]);

    if isempty(logFilePath)||~combineStim
        % Plotting all the figures separately
        obj.saveFigures(figuresFolder,'imageFormat',imageFormat);
    else
        % Plotting combined figures

        % Selecting the columns we care about: power, pulse width, pattern
        relevantParams = {'Pulse Duration','Power','Frequency'};
        unitParams = {'ms','mW','Hz'};

        % For each parameter...
        for ll=1:length(relevantParams)
%             combineExperiments.varyingParam = strrep(relevantParams{ll},' ','');
            combineExperiments.varyingParam = relevantParams{ll};
            combineExperiments.unitParam = unitParams{ll};

            % Finding the relevant combinations for each of these parameters
            pos = find(ismember(paramNames, relevantParams{ll})==1);
            [combinations, varyingParamValues] = findCombinations(M,pos);

            for mm=1:length(combinations)
                % For each of these combinations, wrapping it in a structure
                % adapted to the saveFigures method of the neuronInformation class
                currentCombination = combinations{mm};
                
                % If there is only one element in the combination, it's not
                % really a combination and we don't plot it
                if length(currentCombination)>1
                    combineExperiments.expID = currentCombination;
                    combineExperiments.valuesParam = varyingParamValues{mm};

                    % Then: creating the figure and saving it with a meaningful tag
                    tag = '';
                    for nn=1:length(paramNames)
                        if ~strcmp(paramNames{nn},relevantParams{ll})
                            tag = strcat(tag, paramNames{nn}, '_', ...
                                num2str(M(currentCombination(1),nn)), '_');
                        end
                    end
                    tag(end:end) = '';
                    tag = strrep(tag,'.',',');


                    obj.saveFigures(figuresFolder,'imageFormat',imageFormat,...
                        'combineExperiments',combineExperiments,...
                        'nameTag',tag,'displayimage','off');
                end
            end
        end
    end
end

end % createFigures
