function [goodIDs, duplicateIDs, duplicatePairs, duplicateGroups]...
    = removeDuplicates(dataFolder,varargin)
% [GOODIDS, DUPLICATEIDS, PAIRS, GROUPS] = removeDuplicates(DATAFOLDER,...)
%
% This function identifies the duplicates in a neuron file and returns a
% list of unique neurons, duplicates and duplicate groups. 
% Duplicate identification is done by comparing EIs.  The criterion used
% can be either a threshold on the l2 norm of the difference between neuron
% EIs, or a chi squared test on that same difference. 
%
% Parameters:
%   - dataFolder:  root folder in which a neuron and an EI file can be
%   found.
%
% Optional parameters should be specified in pairs (tag/value):
%   - 'useChiSquared': boolean, false by default.  If set to false l2 norm
%   difference between EIs is used as similarity criterion, otherwise it is
%   chi squared hypothesis test.
%   - 'threshold': double, similarity thresold for EIs comparison.  For the
%   l2 norm difference it should be on the order of a couple of hundreds.
%   For the chi squared test, the default of 0.05 is a conservative value. 
%   For the l2 norm, the default is 5*10^3, conservative value.
%   - 'detectThreshold': boolean, set to true by default.  When true, the
%   algorithm tries to detect an optimal threshold for duplicates
%   identification (the optional 'threshold' parameter above).  If
%   'threshold' is specified detectThreshold should either be left
%   unspecified or set to false.  If 'threshold' is specified and
%   'detectThreshold' is set to true, the last parameter specified will
%   win.
%   - 'useRawNeurons': default is false. When true the neurons are
%   identified in the neurons-raw file, otherwise they're identified in the
%   neurons file.
%   - 'neuronList': when specified, the method will identify duplicates
%   amongst the neurons in this list.
%
% Returns:
%   - goodIDs: a sorted list of unique neurons
%   - duplicateIDs: a sorted list of neuron duplicates
%   - duplicatePairs: an undirected graph of neuron duplicates in array
%   list.  The third column is the "similarity" value between both neurons
%   (l2 norm or reduced chi square)
%   - duplicateGroups: a structure grouping the neurons by duplicate
%   groups. The structure has 3 fields: groupNumber, the number of the
%   duplicate group; bestNeuron:  the ID of the neuron with highest spike
%   count in the group;  allDuplicates: the ID of all the neurons that were
%   identified as duplicates.
%
% Version: v4.04 - 25/07/2011

%% Parameters

% Making sure dataFolder ends by '\' or '/', whichever is right
if dataFolder(end:end)~=filesep
    dataFolder = [dataFolder filesep];
end

% Setting default values
if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    visionPath = 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end

useNeuronRawFile = false;
neuronIDs = [];

useChiSquared = false;
autoThreshold = true;
l2_threshold = 5*10^3;
chi2_threshold = 0.05;

% Checking the optional parameters
nbin = length(varargin);
if mod(nbin,2)==1
    err = MException('MATLAB:InvArgIn', ...
        'Unexpected number of arguments');
    throw(err);
end

% Reading the optional input arguments
tempThresh = 0;
for kk=1:(nbin/2)
    if ~ischar(varargin{kk*2-1})
        err = MException('MATLAB:InvArgIn',...
            'Unexpected additional property');
        throw(err);
    end
    
    switch lower(varargin{kk*2-1})
        case 'usechisquared'
            useChiSquared = varargin{kk*2};
        case 'detectthreshold'
            autoThreshold = varargin{kk*2};
        case 'visionpath'
            visionPath = varargin{kk*2};
        case 'threshold'
            tempThresh = varargin{kk*2};
            autoThreshold = false;
        case 'userawneurons'
            useNeuronRawFile = varargin{kk*2};
        case 'neuronlist'
            neuronIDs = varargin{kk*2};
        otherwise
            err = MException('MATLAB:InvArgIn',...
                'Unknown parameter specified');
            throw(err);
    end
end

% Way of setting the threshold is kind of clumsy... 
if tempThresh>0
    if useChiSquared
        chi2_threshold = tempThresh;
    else
        l2_threshold = tempThresh;
    end
end

%% Linking to the files, setting threshold

% Finding the path to the neuron file in the data folder
contentsDataFolder = dir(dataFolder);

for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isEIFile = strfind(contentsDataFolder(kk).name,'.ei');
%     isParamsFile = strfind(contentsDataFolder(kk).name,'.params');
    if isNeuronFile
        if isNeuronRawFile
            if useNeuronRawFile
                neuron_path =  [dataFolder contentsDataFolder(kk).name];
            end
        else
            if (~useNeuronRawFile)&&...
               (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuron_path =  [dataFolder contentsDataFolder(kk).name];
            end
        end
    end
    if isEIFile
        ei_path = [dataFolder contentsDataFolder(kk).name];
    end
%     if isParamsFile
%         params_path = [dataFolder contentsDataFolder(kk).name];
%     end
end

if ~exist('edu/ucsc/neurobiology/vision/io/RawDataFile','class')
    javaaddpath(visionPath);
end
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(ei_path);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuron_path);
% paramFile = edu.ucsc.neurobiology.vision.io.ParametersFile(param_path);

if useChiSquared
    threshold = chi2_threshold;
else
    threshold = l2_threshold;
end

%% Reading the EI and the neuron file

if isempty(neuronIDs)
    neuronIDs = neuronFile.getIDList();
end
nNeurons = length(neuronIDs);

currentEI = eiFile.getImage(neuronIDs(1));
% Last frame of EI is completely blank: ignored
currentEI = currentEI(:,:,1:end);     
allEIs = zeros(nNeurons,size(currentEI,2));
allVariances = zeros(size(allEIs));

% Loading all the EIs in memory
for kk=1:nNeurons
    currentEI = eiFile.getImage(neuronIDs(kk));
    currentEI_volt = squeeze(currentEI(1,:,:)).';
    currentEI_var = squeeze(currentEI(2,:,:)).';
    
    [EI2D,pos] = max(currentEI_volt);
    var2D = currentEI_var((0:size(currentEI_volt,2)-1)*...
        size(currentEI_volt,1)+pos);
    
    allEIs(kk,:) = EI2D.';
    allVariances(kk,:) = var2D.';
end

%% Computing the EI distance matrix and the chi-square matrix

l2dist = zeros(nNeurons);
chi2dist = zeros(size(l2dist));

for kk=1:nNeurons
    currentEI_kk = allEIs(kk,:);
    currentEI_kk_var = allVariances(kk,:);
    nPoints = numel(currentEI_kk);
    
    for ll=(kk+1):nNeurons
        currentEI_ll = allEIs(ll,:);
        currentEI_ll_var = allVariances(ll,:);
        
        EIDiff = currentEI_kk-currentEI_ll;
        
        % Computing the l2 distance between the EIs of the two neurons
        l2dist(kk,ll) = sum(sum(EIDiff.^2));
        l2dist(ll,kk) = l2dist(kk,ll);
        
        % Regularized chi2 test to take into account the variances as well
        % Note: there should always be 0 varianve in the first column, it
        % should be taken out.
        % This is why divide by nPoints - 2
        aux_var = (EIDiff).^2./(currentEI_kk_var+currentEI_ll_var);
        aux_var = aux_var(2:end);
        
        chi2dist(kk,ll) = sum(aux_var)/(nPoints-2);
        chi2dist(ll,kk) = chi2dist(kk,ll);
    end
end
% If there was a NaN it's due to a zero variance, very unlikely for an
% offdiagonal term
chi2dist(isnan(chi2dist)) = 0;


%% Building the output matrices - either with chi squared or l2 distance

% Building the adjacency matrix for the neuron duplicates: two neurons
% connected in this graph are duplicates
if useChiSquared
    metric = chi2dist;
else
    metric = l2dist;
end

% If we're automatically detecting the threshold, we do it here
if autoThreshold
    allDist = reshape(metric,numel(metric),1);
    allDist = sort(allDist);
    nIdentical = sum(allDist==0);
    [dummyVar,thresholdPos] = max(diff(allDist((nIdentical+1):round(numel(allDist)/2)))); 
    thresholdPos = thresholdPos + nIdentical;
    threshold = allDist(thresholdPos);
end

duplicatePairsMatrix = (metric<=threshold);

% Formatting this undirected graph in an array list form
duplicatePairs = [];
for kk=1:(nNeurons-1)
    for ll=(kk+1):nNeurons
        if duplicatePairsMatrix(kk,ll)
            duplicatePairs = [duplicatePairs; ...
                [double(neuronIDs(kk)) double(neuronIDs(ll)) metric(kk,ll)]];
        end
    end
end

% Building the duplicate groups
duplicateGroups = struct('groupNumber',{},'bestNeuron',{}','allDuplicates',{});
nDuplicateGroups = 0;
visitedNeurons = [-1,-1];
for kk=1:nNeurons
    % For the group we're currently looking at, checking if we've already
    % seen one of the neurons and associating a duplicate group number
    % accordingly
    currentDuplicateGroup = double(neuronIDs(duplicatePairsMatrix(:,kk)));
    commonNeurons = intersect(visitedNeurons(:,1),currentDuplicateGroup);
    if ~isempty(commonNeurons)
        currentDuplicateGroupNumber = visitedNeurons(visitedNeurons(:,1)...
            ==commonNeurons(1),2);
        currentDuplicateGroupNumber = currentDuplicateGroupNumber(1);
    else
        nDuplicateGroups = nDuplicateGroups + 1;
        currentDuplicateGroupNumber = nDuplicateGroups;
    end
    % Affecting the group of neurons we're currently looking at
    duplicateGroups(currentDuplicateGroupNumber).groupNumber = ...
        currentDuplicateGroupNumber;
    duplicateGroups(currentDuplicateGroupNumber).allDuplicates = ...
        [duplicateGroups(currentDuplicateGroupNumber).allDuplicates;...
        currentDuplicateGroup];
    % Adding the neurons in this group to the visited neurons list
    visitedNeurons = [visitedNeurons; [currentDuplicateGroup,...
        currentDuplicateGroupNumber*ones(length(currentDuplicateGroup),1)]];
end

% We might have put in some neurons several times in the duplicate group:
% cleaning this now and selecting the best neuron in the group
goodIDs = zeros(nDuplicateGroups,1);
duplicateIDs = [];
for kk=1:nDuplicateGroups
    duplicateGroups(kk).allDuplicates = unique(...
        duplicateGroups(kk).allDuplicates);
    
    % Finding the best neuron in the group
    maxSpikeCount = 0;
    for ll=1:length(duplicateGroups(kk).allDuplicates)
        thisSpikeCount = neuronFile.getSpikeCount(...
            duplicateGroups(kk).allDuplicates(ll));
        if thisSpikeCount>maxSpikeCount
            maxSpikeCount = thisSpikeCount;
            bestNeuron = duplicateGroups(kk).allDuplicates(ll);
        end
    end
    duplicateGroups(kk).bestNeuron = bestNeuron;
    
    % Building the duplicate matrices we will return
    goodIDs(kk) = duplicateGroups(kk).bestNeuron;
    currentDuplicates = duplicateGroups(kk).allDuplicates;
    currentDuplicates(currentDuplicates==duplicateGroups(kk).bestNeuron) = [];
    if ~isempty(currentDuplicates)
        duplicateIDs = [duplicateIDs; currentDuplicates];
    end
end
goodIDs = sort(goodIDs);
duplicateIDs = sort(unique(duplicateIDs));

end % removeDuplicates
