% How many cells need to be analyzed? 
% For an example dataset (2015-08-17-4) there are 357 neurons

% For all electrodes, how many cells pass by the stimulating electrode. 
% Load ei file .. 
% Testing with an example dataset
eiFilePath = '/Volumes/Analysis/2015-08-17-4/data001/data001.ei';
pathToAnalysisData = '/Volumes/Analysis/2015-08-17-4/data004/'; %elec stim data

% [eiFile, eiPath] = uigetfile('/Volumes/Analysis/*.ei');
% eiFilePath = fullfile(eiPath,eiFile);

[eiM,nIdList] = convertEiFileToMatrix(eiFilePath);
[xC,yC] = getElectrodeCoords512();

threshold = 7; 
cellsToAnalyze = cell(512,1); 
for e = 1:512
    electrodeNo = e;
    cellIDs = [];
    for n = 1:1:size(nIdList,1)
        tempID  = nIdList(n);          % gets the id number of the ith neuron
        ei = squeeze(eiM(:,:,n))';
        eiAmps = max(ei)-min(ei);
        if eiAmps(electrodeNo) > threshold
            cellIDs = cat(1,cellIDs,tempID);
        end
    end
    cellsToAnalyze{e} = cellIDs; 
end

%% How many patterns are there to analyze?
totalNo = 0; 
for e = 1:512
    totalNo = length(cellsToAnalyze{e}) + totalNo; 
end
fprintf('Total number of pattern-cell pairs is %0.0f\n',totalNo); 
%%
currentPair = 0; 
for e = 1:232
    currentPair = length(cellsToAnalyze{e}) + currentPair; 
end
%%
for e = 233:512
    cellIDs = cellsToAnalyze{e} ; 
    patternNo = e; 
    disp( '**************************************************'); 
    disp(['****** analyzing pattern no. ' num2str(e) ' *****************']); 
    disp( '**************************************************'); 
    
    for n = 1:length(cellIDs)
        currentPair = 1+currentPair; 
        nId = cellIDs(n);
        elecStimSorting(pathToAnalysisData,eiFilePath,patternNo,nId);
        disp(['finished analyzing neuron ' num2str(nId) ' for pattern '...
            num2str(patternNo) ' : ' num2str(currentPair) ' of ' num2str(totalNo)]); 
%         disp(['finished analyzing ' num2str(currentPair) ' of ' num2str(totalNo)]); 
    end
end
%% Plotting ..
electrodeNo = 50; 
figure; set(gcf,'Color',[1 1 1]); hold on;
% title(sprintf('Stimulating Electrode %0.0f',electrodeNo)); 
cellIds = cellsToAnalyze{electrodeNo}; 
colors = jet(length(cellIds));
legendEntry = zeros(1,length(cellIds)); 
for nn = 1:length(cellIds); 
    ei = squeeze(eiM(:,:,find(nIdList == cellIds(nn))))';
    eiAmps = max(ei)-min(ei);
    threshVals = find(eiAmps>threshold); 
    legendEntry(nn) = scatter(xC(threshVals),yC(threshVals),eiAmps(threshVals),colors(nn,:),'filled');
end
scatter(xC(electrodeNo),yC(electrodeNo),500,[0 0 0],'filled');
legend(legendEntry,num2str(cellIds)); 
axis image; axis off;
title(sprintf('stimulating electrode %0.0f, threshold %0.0f',...
    electrodeNo,threshold));
