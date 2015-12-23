function altClusterPlots(main)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



filePath = get(main.filePath, 'String');
residLim = [str2double(get(main.residualLimitLow, 'String')) str2double(get(main.residualLimitHigh, 'String'))];
currentMovie = str2double(get(main.movieNo, 'String'));






%load data
if get(main.analyzeModeButton, 'Value') %analysis mode
    fileName = get(main.fileName, 'String');
    temp = load([filePath fileName]);
    elecResp = temp.elecResp;
    
    movieInd = find(elecResp.stimInfo.movieNos == currentMovie);
    
    if exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
        dataTraces=NS_ReadPreprocessedData([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo)],...
            '', 0, elecResp.stimInfo.patternNo, currentMovie, 99999);
    elseif exist([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie)], 'file')
        dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
            currentMovie, 99999);
    else
        warnH = warndlg([elecResp.names.data_path filesep 'p' num2str(elecResp.stimInfo.patternNo) '_m' num2str(currentMovie) ' couldn''t be found']);
        uiwait(warnH)
        return
    end
    centerChannel = elecResp.cells.recElec;
    
    eiNorm = norm(elecResp.cells.mainEI(centerChannel,:)); %used to draw a threshold line in linkage plot
    clusterElecs = getCluster(centerChannel);
    eiNormSurrElecs = norm(reshape(elecResp.cells.mainEI(clusterElecs,:)', 1, []));
    
    
    if ~isempty(elecResp.analysis.type{movieInd})
        succBin = elecResp.analysis.latencies{movieInd}~=0;
    else
        succBin = [];
    end
    
    residLim = elecResp.analysis.details.residCalcWindow{movieInd};
    
    if isempty(residLim)
        residLim = [str2double(get(main.residualLimitLow, 'String')) str2double(get(main.residualLimitHigh, 'String'))];
    end
    
    
else %explore mode
    patternNoStr = get(main.patternNo, 'String');
    patternNo = str2double(patternNoStr);
    centerChannel = str2double(get(main.centerElec, 'String'));
    residLim = [str2double(get(main.residualLimitLow, 'String')) str2double(get(main.residualLimitHigh, 'String'))];
    
    %determine whether data has been clustered temporarily
    try
        if isnan(patternNo);
            succBin = main.tempAnalysisNoSave;
        else
            succBin = main.tempAnalysis{patternNo, currentMovie};
        end
    catch
        succBin = [];
    end
    succBin = logical(succBin);
    
    if exist([filePath filesep 'p' patternNoStr filesep 'p' patternNoStr '_m' num2str(currentMovie)], 'file')
        dataTraces=NS_ReadPreprocessedData([filePath filesep 'p' patternNoStr], '', 0, patternNo, currentMovie, 99999);
    elseif exist([filePath filesep 'p' patternNoStr '_m' num2str(currentMovie)], 'file')
        dataTraces=NS_ReadPreprocessedData(filePath, '', 0, patternNo, currentMovie, 99999);
    end
end

% determine which channels to use
clusterElecs = getCluster(centerChannel);
nChannels = length(clusterElecs);

%if in explore mode OR analyzed but without all successes, determine which
%ei channels are checked and only use those channels
if ~get(main.analyzeModeButton, 'Value') || (~all(succBin) && ~isempty(elecResp.analysis.type))
    for ii = nChannels:-1:1
        if ~get(main.aCheckBox{1,ii}, 'Value') %if it's unchecked
            clusterElecs(ii) = [];
        end
    end
end

validData1Elec = squeeze(dataTraces(:, centerChannel, residLim(1):residLim(2)));

validDataSurrElecs = dataTraces(:, clusterElecs, residLim(1):residLim(2));
validDataSurrElecs = permute(validDataSurrElecs,[1 3 2]);
validDataSurrElecs = reshape(validDataSurrElecs, size(validDataSurrElecs,1), []);

%linkage analysis
dist = pdist(validData1Elec);
distSurrElecs = pdist(validDataSurrElecs);
linkSing = linkage(dist, 'average');
linkSingSurrElecs = linkage(distSurrElecs, 'average');

%PCA
[~, PCs] = princomp(validData1Elec);
[~, PCsSurrElecs] = princomp(validDataSurrElecs);

%% plotting
figure('position', [100 200 800 800])

if ~isempty(succBin)
    %dendrogram single electrode data
    axes('position', [0.1 0.55 0.35 0.4]);
    [hBranches, ~, perm] = dendrogram(linkSing, 0, 'labels', char(zeros(size(validData1Elec,1),1)));
        
    %determine which branches contain only successful leaves
    succLeaves = find(succBin);
    succBranches = [];
    for ii = 1:size(linkSing,1)
        if (any(succLeaves == linkSing(ii,1)) || any(succBranches+size(validData1Elec,1) == linkSing(ii,1))) && ...
                (any(succLeaves == linkSing(ii,2)) || any(succBranches+size(validData1Elec,1) == linkSing(ii,2))) %%branches are number starting from nLeaves+1
            succBranches = [succBranches ii]; %numbered as they are in hBranches (without the +nLeaves)
        end
    end
        
    hold on
    for ii = 1:length(succBin)
        if succBin(ii)
            [x,~] = ind2sub(size(linkSing(:,1:2)), find(linkSing(:,1:2)==ii));
            plot([find(perm == ii) find(perm == ii)], [0 linkSing(x,3)], 'r-')
        end
    end
    set(hBranches(succBranches), 'color', 'r')
    if get(main.analyzeModeButton, 'Value') %analysis mode
        plot([1 size(validData1Elec,1)], [eiNorm eiNorm], 'k')
    end
    hold off
    title('center electrode only')
    ylabel('mean Euclidean distance')

    
    %dendrogram with surrounding electrode data
    axes('position', [0.1 0.05 0.35 0.4]);
    [hBranches, ~, perm] = dendrogram(linkSingSurrElecs, 0, 'labels', char(zeros(size(validData1Elec,1),1)));

    %determine which branches contain only successful leaves
    succBranches = [];
    for ii = 1:size(linkSingSurrElecs,1)
        if (any(succLeaves == linkSingSurrElecs(ii,1)) || any(succBranches+size(validData1Elec,1) == linkSingSurrElecs(ii,1))) && ...
                (any(succLeaves == linkSingSurrElecs(ii,2)) || any(succBranches+size(validData1Elec,1) == linkSingSurrElecs(ii,2))) %%branches are number starting from nLeaves+1
            succBranches = [succBranches ii]; %numbered as they are in hBranches (without the +nLeaves)
        end
    end
    
    hold on
    for ii = 1:length(succBin)
        if succBin(ii)
            [x,~] = ind2sub(size(linkSingSurrElecs(:,1:2)), find(linkSingSurrElecs(:,1:2)==ii));
            plot([find(perm == ii) find(perm == ii)], [0 linkSingSurrElecs(x,3)], 'r-')
        end
    end
    set(hBranches(succBranches), 'color', 'r')

    if get(main.analyzeModeButton, 'Value') %analysis mode
        plot([1 size(validData1Elec,1)], [eiNormSurrElecs eiNormSurrElecs], 'k')
    end
    hold off
    title(['electrodes ' num2str(clusterElecs)])
    ylabel('mean Euclidean distance')

else
    axes('position', [0.1 0.55 0.35 0.4]);
    dendrogram(linkSing, 0, 'labels', char(zeros(size(validData1Elec,1),1))); %0 signifies to display all branches
    if get(main.analyzeModeButton, 'Value') %analysis mode
        hold on
        plot([1 size(validData1Elec,1)], [eiNorm eiNorm], 'k')
        hold off
    end
    title('center electrode only')
    ylabel('mean Euclidean distance')
    
    
    axes('position', [0.1 0.05 0.35 0.4]);
    dendrogram(linkSingSurrElecs, 0, 'labels', char(zeros(size(validData1Elec,1),1)));
    if get(main.analyzeModeButton, 'Value') %analysis mode
        hold on
        plot([1 size(validData1Elec,1)], [eiNormSurrElecs eiNormSurrElecs], 'k')
        hold off
    end
    title(['electrodes ' num2str(clusterElecs)])
    ylabel('mean Euclidean distance')
end

if ~isempty(succBin)
    axes('position', [0.55 0.55 0.4 0.4]);
    hold on
    plot(PCs(succBin, 1), PCs(succBin, 2), 'r.')
    plot(PCs(~succBin, 1), PCs(~succBin, 2), 'k.')
    hold off
    title('center electrode only')
    xlabel('PC1'); ylabel('PC2')
    
    axes('position', [0.55 0.05 0.4 0.4]);
    hold on
    plot(PCsSurrElecs(succBin, 1), PCsSurrElecs(succBin, 2), 'r.')
    plot(PCsSurrElecs(~succBin, 1), PCsSurrElecs(~succBin, 2), 'k.')
    hold off
    title(['electrodes ' num2str(clusterElecs)])    
    xlabel('PC1'); ylabel('PC2')
else
    axes('position', [0.55 0.55 0.4 0.4]);
    plot(PCs(:,1), PCs(:,2), 'k.')
    title('center electrode only')
    xlabel('PC1'); ylabel('PC2')
    
    axes('position', [0.55 0.05 0.4 0.4]);
    plot(PCsSurrElecs(:, 1), PCsSurrElecs(:, 2), 'k.')
    title(['electrodes ' num2str(clusterElecs)])
    xlabel('PC1'); ylabel('PC2')
end


