function neuronPairsRefVsNew = crossIdentifyNeuronIDs(dataFolderRef,dataFolderNew,...
                                refNeuronIDs, newNeuronIDs)
% crossIdentifyNeuronIDs(dataFolderRef,dataFolderNew,...
%                        newNeuronIDs,refNeuronIDs)
%
% This function compares the neurons in a data set (the one in
% dataFolderNew) with the neurons in a reference data set contained in
% dataFolderRef. 
% It returns a nNeurons*2 matrix, where nNeurons is the number of neurons
% in the reference dataset.  The first column in this matrix corresponds to
% the neuronID of a neuron in the new dataset, and the second column
% corresponds to the neuron in the reference dataset that looks the most
% like it.  
% Similarities between neurons are measured using the l2 norm of the
% difference between EIs, the smaller this quantity is the better. 
%
% Parameters:
%    - dataFolderRef: reference dataset, we want to find out which neurons
%    in the new dataset correspond to which neurons in this dataset (these
%    neurons will be the first column of the output).
%    - dataFolderNew: new dataset, we want to find its neuron in the
%    reference dataset (first column of the output)
%
% Returns:
%    - neuronPairsRefVsNew: a nNeurons*2 matrix representing the most
%    likely pairing of neurons from the ref dataset into the new dataset.
%    The likelyhood of pairing is estimated from the EIs.
%
% Version: 27/02/2012

%% Parameters

nSim = 3;

if dataFolderRef(end:end)~=filesep
    dataFolderRef = [dataFolderRef filesep];
end
if dataFolderNew(end:end)~=filesep
    dataFolderNew = [dataFolderNew filesep];
end

if isunix
    visionPath = '/home/ggoetz/Research/Vision/Vision815/Vision.jar';
else
    visionPath = 'C:\Users\ggoetz\Research\Vision\Vision815\Vision.jar';
end

%% Linking to the files

if ~exist('edu/ucsc/neurobiology/vision/io/NeuronFile','class')
    javaaddpath(visionPath);
end

% Finding the reference neuron and ei files
contentsDataFolder = dir(dataFolderRef);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isEIFile = strfind(contentsDataFolder(kk).name,'.ei');
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuronRef_path =  [dataFolderRef contentsDataFolder(kk).name];
                break;
            end
        end
    end
    if isEIFile
        eiRef_path = [dataFolderRef contentsDataFolder(kk).name];
    end
end
neuronRefFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronRef_path);
eiRefFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiRef_path);

% Finding the new neuron and ei files
contentsDataFolder = dir(dataFolderNew);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isEIFile = strfind(contentsDataFolder(kk).name,'.ei');
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuronNew_path =  [dataFolderNew contentsDataFolder(kk).name];
                break;
            end
        end
    end
    if isEIFile
        eiNew_path = [dataFolderNew contentsDataFolder(kk).name];
    end
end
neuronNewFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronNew_path);
eiNewFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiNew_path);

%% Reading the EIs 

if ~exist('newNeuronIDs','var')||isempty(newNeuronIDs)
    newNeuronIDs = neuronNewFile.getIDList();
end
if ~exist('refNeuronIDs','var')||isempty(refNeuronIDs)
    refNeuronIDs = neuronRefFile.getIDList();
end

% Making sure the neuron ID lists are column vectors
if size(newNeuronIDs,1)<size(newNeuronIDs,2)
    newNeuronIDs = newNeuronIDs.';
end
if size(refNeuronIDs,1)<size(refNeuronIDs,2)
    refNeuronIDs = refNeuronIDs.';
end

nNeuronsRef = length(refNeuronIDs);
nNeuronsNew = length(newNeuronIDs);

% Initialization of the EI storage matrices
currentEI = eiRefFile.getImage(refNeuronIDs(1));
% Last frame of EI is completely blank: ignored
currentEI = currentEI(:,:,1:end);     
allEIsRef = zeros(nNeuronsRef,size(currentEI,2));
allVariancesRef = zeros(size(allEIsRef));
allEIsNew = zeros(nNeuronsNew,size(currentEI,2));
allVariancesNew = zeros(size(allEIsNew));

% Loading all the EIs in memory, first the new dataset then the ref one
for kk=1:nNeuronsNew
    currentEI = eiNewFile.getImage(newNeuronIDs(kk));
    currentEI_volt = squeeze(currentEI(1,:,:)).';
    currentEI_var = squeeze(currentEI(2,:,:)).';
    
    [EI2D,pos] = max(currentEI_volt);
    var2D = currentEI_var((0:size(currentEI_volt,2)-1)*...
        size(currentEI_volt,1)+pos);
    
    allEIsNew(kk,:) = EI2D.';
    allVariancesNew(kk,:) = var2D.';
end
for kk=1:nNeuronsRef
    currentEI = eiRefFile.getImage(refNeuronIDs(kk));
    currentEI_volt = squeeze(currentEI(1,:,:)).';
    currentEI_var = squeeze(currentEI(2,:,:)).';
    
    [EI2D,pos] = max(currentEI_volt);
    var2D = currentEI_var((0:size(currentEI_volt,2)-1)*...
        size(currentEI_volt,1)+pos);
    
    allEIsRef(kk,:) = EI2D.';
    allVariancesRef(kk,:) = var2D.';
end


%% Comparing the EIs and building the neuronPairs matrix

neuronPairsRefVsNew = zeros(length(refNeuronIDs),nSim+1);

for kk=1:nNeuronsRef
    % Computing the l2 distance between EIs
    l2dist = zeros(1,nNeuronsNew);
    currentEI_kk = allEIsRef(kk,:);
    for ll = 1:nNeuronsNew
        currentEI_ll = allEIsNew(ll,:);
        EIDiff = currentEI_kk - currentEI_ll;
        l2dist(ll) = sum(sum(EIDiff.^2));
    end
    
    % Findind the best matches
    [~,pos] = sort(l2dist);
    
    % Storing the nearest neighbor into the neuronPairs matrix
    neuronPairsRefVsNew(kk,:) = [refNeuronIDs(kk), newNeuronIDs(pos(1:nSim)).'];
end

end % crossCheckDatasets