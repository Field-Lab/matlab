function neuronPairsRefVsNew = crossIdentifyNeuronIDs_NB(dataFolderRef, dataFolderNew,...
                                refNeuronIDs, newNeuronIDs, autoClassify)
% crossIdentifyNeuronIDs(dataFolderRef,dataFolderNew,...
%                        refNeuronIDs,newNeuronIDs)
%
% This function compares the neurons in a data set (the one in
% dataFolderNew) with the neurons in a reference data set contained in
% dataFolderRef. 
% It returns a nNeurons*2 matrix, where nNeurons is the number of neurons
% in the reference dataset.  The first column in this matrix corresponds to
% the neuronIDs of neurons in the reference dataset, and the second column
% corresponds to the neurons in the new dataset that looks the most
% like the reference neurons. 
% Those matches are the argmin of the l2 norm of the difference between 
% the reference neuron EI and EIs of neurons in the new dataset.
%
% Parameters:
%    - dataFolderRef: reference dataset, we want to find out which neurons
%    in the new dataset correspond to which neurons in this dataset (these
%    reference neurons will be in the first column of the output).
%    dataFolderRef should be a path to a folder containing both an EI and
%    a neurons file (neurons-raw files will be ignored).
%    - dataFolderNew: new dataset, we want to match its neurons with those 
%    of the reference dataset. dataFolderNew should be a path to a folder
%    containing both an EI and a neurons file (neurons-raw files will be 
%    ignored).
%
% Optional:
%    - refNeuronsIDs, newNeuronIDs: if those two lists are specified, only
%    the neurons in this list will be considered for pairing.
%
% Returns:
%    - neuronPairsRefVsNew: a nNeurons*2 matrix representing the most
%    likely pairing of neurons from the ref dataset into the new dataset.
%    The likelyhood of pairing is estimated from the EIs.
%
% Version: 27/02/2012

%% Parameters

% If you want to output more than one match candidate for each neuron,
% you can increase nSim
nSim = 3; 
threshold = 3*15000;

if dataFolderRef(end:end)~=filesep
    dataFolderRef = [dataFolderRef filesep];
end
if dataFolderNew(end:end)~=filesep
    dataFolderNew = [dataFolderNew filesep];
end

if isunix
    visionPath = '/Lab/Development/vision7/Vision.jar';
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
    if autoClassify
            isParamFile = strfind(contentsDataFolder(kk).name,'.params');
    end
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuronRef_path =  [dataFolderRef contentsDataFolder(kk).name];                
            end
        end
    end
    if isEIFile
        eiRef_path = [dataFolderRef contentsDataFolder(kk).name];
    end
    if autoClassify       
        if isParamFile
            paramRef_path = [dataFolderRef contentsDataFolder(kk).name];
        end
        
    end
end
neuronRefFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronRef_path);
eiRefFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiRef_path);
if autoClassify
    paramRefFile = edu.ucsc.neurobiology.vision.io.ParametersFile(paramRef_path);
end
% Finding the new neuron and ei files
contentsDataFolder = dir(dataFolderNew);
for kk=1:length(contentsDataFolder)
    isNeuronFile = strfind(contentsDataFolder(kk).name,'.neurons');
    isNeuronRawFile = strfind(contentsDataFolder(kk).name,'.neurons-raw');
    isEIFile = strfind(contentsDataFolder(kk).name,'.ei');
    if autoClassify
        isParamFile = strfind(contentsDataFolder(kk).name,'.params');
    end
    if isNeuronFile
        if isNeuronRawFile
        else
            if (isNeuronFile+length('.neurons')-1)==length(contentsDataFolder(kk).name)
                neuronNew_path =  [dataFolderNew contentsDataFolder(kk).name];
            end
        end
    end
    if isEIFile
        eiNew_path = [dataFolderNew contentsDataFolder(kk).name];
    end
    if autoClassify
        if isParamFile
            paramNew_path = [dataFolderNew contentsDataFolder(kk).name];
        end
    end
end
neuronNewFile = edu.ucsc.neurobiology.vision.io.NeuronFile(neuronNew_path);
eiNewFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiNew_path);
if autoClassify
    paramNewFile = edu.ucsc.neurobiology.vision.io.ParametersFile(paramNew_path);

end

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
counter = 0;
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
    if l2dist(pos(1)) > threshold
        neuronPairsRefVsNew(kk,2) = nan;
        counter = counter +1;
    end
end

if autoClassify
    for i = 1:length(neuronPairsRefVsNew)
        if ~isnan(neuronPairsRefVsNew(i,2))
            class = paramRefFile.getCell(neuronPairsRefVsNew(i,1), 'classID');
            string = char(class.toString);
            begin = strfind(string, 'All');
            paramNewFile.setCell(neuronPairsRefVsNew(i,2), 'classID', string(begin:end));
        end
    end
end

match_rate = (nNeuronsRef - counter)/nNeuronsRef
mean(l2dist)
if autoClassify
    paramNewFile.close(1);
    paramRefFile.close(1);
end

end % crossCheckDatasets