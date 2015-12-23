function computeActivationData(statsFolder,logfilePath,resultFolder,...
    powerScalingFactor,powerToIrradianceFactor,varargin)
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
%
% Optional parameters, specified by key/value pairs:
%   - 'createTextFiles': boolean, false by default. When set to true the
%   activation data will be saved as a .mat and .txt file instead of .mat
%   only.
%   - 'pValue': significance of stimulated bins. By default 0.01.
%   - 'createClassification': boolean, by default true. A classification
%   file is created with neurons above threshols in the 'Stimulated'
%   category, neurons without stimulation in 'No stimulation' and neurons
%   stimulated but under threshold in 'Weakly stimulated'.
%   - positiveOnly: boolean, default is false. If set to false will count
%   suppression as stimulation, otherwise will only count actual spikes
%   elicited above baseline.
%
% Returns: []
%
% Todo: add start and end of integration as optional parameters
%   
% Version: v5.00 - 03/21/2013
%

%%

% Setting default parameters
createTextFile = false;              % save data as mat and txt or mat only
pThreshold = 0.01;
createClassification = true;
positiveOnly = false;

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
        case 'pvalue'
            pThreshold = varargin{kk*2};
        case 'positiveonly'
            positiveOnly = varargin{kk*2};
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

%% If write to text file and old file, remove the old one

if createClassification && exist(fullfile(resultFolder, 'allNeurons.classification'),'file')
    delete(fullfile(resultFolder, 'allNeurons.classification'));
%     delete(fullfile(resultFolder, 'allNeurons_stimVsSupp.classification'));
end

%% Processing each neuron data
labels = {'Power@Logfile-mW','Power@Array-mW','Irradiance-mW/mm^2',...
    'Duration','Frequency','Image projected','nSpikes', 'SEM','nSpikes+','nSpikes-'};
for kk=1:nNeurons
    % Loading
    load([statsFolder neuronNames(kk).name]);
    
    % Creating the data structure
    activationData.labels = labels;
    activationData.data = zeros(nExperiments,length(labels));
    
    % Computing the number of spikes elicited for each experiment
    for ll=1:nExperiments
        % Getting the baseline statistics. For now done from last 40%,
        % but it would be better to have spontaneous data every time. 
        % startTimeBaseline and endTimeBaseline have to be rounded to the
        % nearest multiple of 5.
        startTimeBaseline = round(0.6/fStim(ll)*1000/5)*5;
        endTimeBaseline = floor(1/fStim(ll)*1000/5)*5;
        [nSpikesBackground, semBackground, nTrials, nBins] = getNumberOfSpikes(...
            startTimeBaseline, endTimeBaseline, obj.PSTHs(ll));
        nSpikesExpected = nSpikesBackground/nBins;
        % Estimate from SEM is biased, hence the correcting factor.
        % Other assumption made is that what happens in a bin is
        % independent from what happens in others, which is an
        % approximation but should be close to true.
        varSpikesExpected = nTrials/(nTrials-1)*semBackground*nTrials/nBins;
        
        % Going over all the bins and finding the statistically significant
        % ones, adding their contribution to the number of spikes elicited.
        nSpikesElicited = 0;
        nPosSpikes = 0;
        nNegSpikes = 0;
        varSpikesElicited = 0;
        for mm=3:length(obj.PSTHs(ll).data)-1 % First bins are usually blanked or artifacts
            nSpikesBin = obj.PSTHs(ll).data(mm);
            varBin = obj.PSTHs(ll).variances(mm);
            
            tValue = abs((nSpikesBin-nSpikesExpected)/...
                sqrt(varSpikesExpected/nTrials + varBin/nTrials));
            pValue = 1-tcdf(tValue, nTrials-1);
            
            if pValue < pThreshold
                if nSpikesBin>nSpikesExpected % Stimulation
                    nSpikesElicited = nSpikesElicited + (nSpikesBin - nSpikesExpected);
                    nPosSpikes = nPosSpikes + (nSpikesBin - nSpikesExpected);
                else % Suppression
                    if ~positiveOnly
                        nSpikesElicited = nSpikesElicited - (nSpikesBin - nSpikesExpected);
                        nNegSpikes = nNegSpikes + (nSpikesBin - nSpikesExpected);
                    end
                end
                varSpikesElicited = varSpikesElicited + varBin;
            end
        end
        
        % Updating the activationData structure used to return the results
        activationData.data(ll,1) = powerLevels(ll);
        activationData.data(ll,2) = powerLevels(ll)*powerScalingFactor;
        activationData.data(ll,3) = powerLevels(ll)*powerScalingFactor*...
            powerToIrradianceFactor;
        activationData.data(ll,4) = pulseWidths(ll);
        activationData.data(ll,5) = fStim(ll);
        activationData.data(ll,6) = imProjected(ll);
        activationData.data(ll,7) = nSpikesElicited;
        activationData.data(ll,8) = varSpikesElicited/nTrials; % Note: this is the SEM
        activationData.data(ll,9) = nPosSpikes;
        activationData.data(ll,10) = nNegSpikes;
    end
    
    % Saving the results
    save([resultFolder neuronNames(kk).name(1:(end-4)) '_act.mat'],...
        'activationData');
    
    % Also writing a text file if required
    if createTextFile
        writeToTextFile(activationData.data,neuronNames(kk).name(1:(end-4)),labels,...
            [resultFolder neuronNames(kk).name(1:(end-4)) '_act.txt']);
    end
    
    % Add the neuron to the classification file if required
    if createClassification
        writeToStimulationClassificationFile(activationData.data(:,7),neuronNames(kk).name(7:(end-4)),...
            fullfile(resultFolder, 'allNeurons.classification'));
%         writeToStimVsSuppFile(activationData.data(:,9),activationData.data(:,10),...
%             neuronNames(kk).name(7:(end-4)),fullfile(resultFolder, 'allNeurons_stimVsSupp.classification'));
    end
end

end %powerStr

function [nSpikes, sem, nTrials, nBins] = getNumberOfSpikes(startTime,endTime,PSTH)
% startTime, endTime specified in ms. 
% Assuming bin size of 5 ms in the PSTH, Fs = 20,000 Hz
% Returns the standard error of the mean.

x = (floor(startTime/5)+1):(floor(endTime/5));
nSpikes = sum(PSTH.data(x));
sem = sum(PSTH.variances(x))/PSTH.numberOfPulses; % Square of SEM. Addition because Poisson distr. assumption
nTrials = PSTH.numberOfPulses;
nBins = length(x);

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

function writeToStimulationClassificationFile(nSpikesElicited, neuronName, fileName)
% Writes the neuron classification to the specified file. Compares the
% number of elicited spikes to a fixed threshold.

thresholdStim = 0.5;

fid = fopen(fileName, 'a');
maxNumSpikes = max(nSpikesElicited);
if maxNumSpikes==0
    fprintf(fid, '%s  All/No stimulation\n', neuronName);
else
    if maxNumSpikes<thresholdStim
        fprintf(fid, '%s  All/Weakly stimulated\n', neuronName);
    else
        fprintf(fid, '%s  All/Stimulated\n', neuronName);
    end
end

fclose(fid);

        
end % writeToStimulationClassificationFile

% function writeToStimVsSuppFile(nPosSpikes, nNegSpikes, neuronName, fileName)
% % Writes the neuron classification to the specified file
% 
% fid = fopen(fileName, 'a');
% 
% [~,ixx] = max(nPosSpikes-nNegSpikes);
% 
% fprintf(fid, '%s  %2.4f  %2.4f\n', neuronName, nPosSpikes(ixx), nNegSpikes(ixx));
% fclose(fid);
%         
% end % writeToStimVsSuppFile