function computeActivationData(statsFolder,logfilePath,resultFolder,...
    powerScalingFactor,powerToIrradianceFactor,ratType,varargin)
% computeActivationData(statsFolder,logfilePath,resultFolder,
%    powerScalingFactor,powerToIrradianceFactor,ratType,varargin)
%
% This function computes the activation data from the neuron information
% data. 
%
% Parameters:
%   - statsFolder: path to the statistics folder created by analyzeSpikes
%   - logfilePath: path to the logfile corresponding to the experiment
%   - resultFolder: output path
%   - powerScalingFactor: used to normalize power in logfiles
%   - powerToIrradianceFactor: max irradiance
%   - RatType: 'WT' or 'RCS', if results not satisfying start and end of
%   integration should be changed in the function. 
%
% Optional parameters, specified by key/value pairs:
%   - 'createTextFiles': boolean, false by default. When set to true the
%   activation data will be saved as a .mat and .txt file instead of .mat
%   only.
%
% Returns: []
%
% Todo: add start and end of integration as optional parameters
%   
% Version: v4.04 - 05/08/2011
%

%%

% Setting default parameters
createTextFile = false;              % save data as mat and txt or mat only

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
        case 'createtextfiles'
            createTextFile = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Making sure the folders ends by '\' or '/', whichever is right
if statsFolder(end:end)~=filesep
    statsFolder = [statsFolder filesep];
end
if resultFolder(end:end)~=filesep
    resultFolder = [resultFolder filesep];
end

%% Finding the neuron names and the number of experiments, power levels and PW

contentsStatsFolder = dir(statsFolder);
neuronNames = struct('name','');
nNeurons = 0;
for kk=1:length(contentsStatsFolder)
    if strfind(contentsStatsFolder(kk).name,'.mat')
        nNeurons = nNeurons + 1;
        neuronNames(nNeurons).name = contentsStatsFolder(kk).name;
    end
end

if ~isempty(logfilePath)
    [M,paramNames] = readLogFile(logfilePath);
end
powerLevels = M(:,ismember(paramNames,'Power'));
pulseWidths = M(:,ismember(paramNames,'Pulse Duration'));
fStim = M(:,ismember(paramNames,'Frequency'));
imProjected = M(:,ismember(paramNames,'Image projected'));
nExperiments = length(powerLevels);

if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

%% Processing each neuron data
labels = {'Power@Logfile-mW','Power@Array-mW','Irradiance-mW/mm^2',...
    'Duration','Frequency','Image projected','nSpikes','variance'};
for kk=1:nNeurons
    % Loading
    load([statsFolder neuronNames(kk).name]);
    
    % Creating the data structure
    activationData.labels = labels;
    activationData.data = zeros(nExperiments,length(labels));
    
    % Computing the number of spikes elicited for each experiment
    for ll=1:nExperiments
        switch ratType
            case 'WT'
                % Determining where we're looking for spikes and baseline
                % rate
                startTimeStim = 100;
                endTimeStim = min(220,floor(1000/fStim(ll)/5)*5);            % Original value: 40
                startTimeBaseline = 220;                                     % Original value: 125
                endTimeBaseline = min(340,floor(1000/fStim(ll)/5)*5);       % Original value: 155
                
                % Getting the firing information
                [nSpikesStim, var1] = getNumberOfSpikes(startTimeStim,...
                    endTimeStim,obj.PSTHs(ll));
                [nSpikesBaseline, var2] = getNumberOfSpikes(startTimeBaseline,...
                    endTimeBaseline,obj.PSTHs(ll));
                
                % If both intervals are of different lengths, we rescale
                % the baseline
                if (endTimeBaseline-startTimeBaseline)~=(endTimeStim-startTimeStim)
                    if (endTimeBaseline-startTimeBaseline)>0
                        nSpikesBaseline = nSpikesBaseline*(endTimeStim-...
                            startTimeStim)/(endTimeBaseline-startTimeBaseline);
                    end
                end
                
            case 'RCS'
                % Determining where we're looking for spikes and baseline
                % rate
                startTimeStim = 5;
                endTimeStim = min(80,floor(1000/fStim(ll)/5)*5);            % Original value: 80
                startTimeBaseline = 90;                                    % Original value: 125
                endTimeBaseline = min(165,floor(1000/fStim(ll)/5)*5);       % Original value: 245
                
                % Getting the firing information
                [nSpikesStim, var1] = getNumberOfSpikes(startTimeStim,...
                    endTimeStim,obj.PSTHs(ll));
                [nSpikesBaseline, var2] = getNumberOfSpikes(startTimeBaseline,...
                    endTimeBaseline,obj.PSTHs(ll));
                
                % If both intervals are of different lengths, we rescale
                % the baseline
                if (endTimeBaseline-startTimeBaseline)~=(endTimeStim-startTimeStim)
                    if (endTimeBaseline-startTimeBaseline)>0
                        nSpikesBaseline = nSpikesBaseline*(endTimeStim-...
                            startTimeStim)/(endTimeBaseline-startTimeBaseline);
                    end
                end
                
        end
        
        activationData.data(ll,1) = powerLevels(ll);
        activationData.data(ll,2) = powerLevels(ll)*powerScalingFactor;
        activationData.data(ll,3) = powerLevels(ll)*powerScalingFactor*...
            powerToIrradianceFactor;
        activationData.data(ll,4) = pulseWidths(ll);
        activationData.data(ll,5) = fStim(ll);
        activationData.data(ll,6) = imProjected(ll);
        activationData.data(ll,7) = nSpikesStim - nSpikesBaseline;
        activationData.data(ll,8) = (var1 + var2);
        activationData.data(ll,9) = abs((nSpikesStim - nSpikesBaseline)/sqrt(var2)); % Significance of stimulation
        if isnan(activationData.data(ll,9))
            activationData.data(ll,9) = Inf;
        end
    end
    
    % Saving the results
    save([resultFolder neuronNames(kk).name(1:(end-4)) '_act.mat'],...
        'activationData');
    
    % Also writing a text file if required
    if createTextFile
        writeToTextFile(activationData.data,neuronNames(kk).name(1:(end-4)),labels,...
            [resultFolder neuronNames(kk).name(1:(end-4)) '_act.txt']);
    end
end

end %powerStr

function [nSpikes, var] = getNumberOfSpikes(startTime,endTime,PSTH)
% startTime, endTime specified in ms. 
% Assuming bin size of 5 ms in the PSTH, Fs = 20,000 Hz

x = (floor(startTime/5)+1):(floor(endTime/5));
nSpikes = sum(PSTH.data(x));
var = sum(PSTH.variances(x))/PSTH.numberOfPulses; % Square of SEM

end % getNumberOfSpikes

function writeToTextFile(activationData, neuronName, labels, fileName)
% Writes the stimulation data to a text file, sorting by pulseWidth

nParams = length(labels);
activationData = sortrows(activationData,find(ismember(labels,'Duration')));

fid = fopen(fileName,'w+');

fprintf(fid,[neuronName '\r\n']);
for kk=1:nParams
    fprintf(fid, [labels{kk} '\t']);
end
fprintf(fid, '\r\n');
fprintf(fid,[repmat('%4.4f\t',1,nParams) '\r\n'],activationData.');

fclose(fid);

end % writeToTextFile