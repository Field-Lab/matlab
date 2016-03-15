% Input list of cells to analyze: 
% Get a list of ON parasol cells from 2015-11-09-3. 
% fPath = '/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/';
% visionDataPath = '/Volumes/Analysis/2015-11-09-3/data000/data000'; 

fPath = '/Volumes/Analysis/2015-10-06-3/data001-data002-autosort-template000/';
visionDataPath = '/Volumes/Analysis/2015-10-06-3/data000/data000'; 

eipath = [visionDataPath '.ei']; 
datarun = load_data(visionDataPath);  
datarun = load_params(datarun); 
datarun = load_sta(datarun); 
datarun = load_ei(datarun,'all'); 
fillcolors = [0.275 0.094, 0.471;
    0.471 0.094 0.29;
    0.29 0.471 0.094;
    0.094, 0.471, 0.275];

celltype = datarun.cell_types{2}.name; 
target_cell_group = datarun.cell_types{1}.cell_ids; % List of cell IDs of the group to target. 
group_to_analyze_against = [datarun.cell_types{1}.cell_ids datarun.cell_types{2}.cell_ids]; % target_cell_group; %
% [allVals, actProb, actAmp] = genActThreshSpatialMaps_auto(fPath,eipath,cell_ids); 
%%
disp('Find pattern-cell pairs to analyze (patternNo = stim electrode)..'); 
[eiM,nIdList] = convertEiFileToMatrix([visionDataPath '.ei']);
threshold =10; 
sortingThreshold = 30; 
cellsToAnalyze = cell(512,1); 
for e = 1:512
    electrodeNo = e;
    cellIDs = [];
    for n = 1:1:size(group_to_analyze_against,2)
        tempID  = group_to_analyze_against(n);          % gets the id number of the ith neuron
        
        ei = squeeze(eiM(:,:,find(tempID == nIdList)))';
        eiMin = abs(min(ei));   
        if eiMin(electrodeNo) > threshold && (max(eiMin(:)) > sortingThreshold)
            
            minIdx = find(ei(:,electrodeNo) == min(ei(:,electrodeNo)));
            max1 = max(ei(1:minIdx,electrodeNo));
            max2 = max(ei(minIdx:end,electrodeNo));
            if max1>2
%                 plot(ei(:,electrodeNo),'r'); title(['AXON max first phase is ' num2str(max1)])
            else % Only includes somatic activation
%                 plot(ei(:,electrodeNo),'b'); title(['SOMA max first phase is ' num2str(max1)])
                cellIDs = cat(1,cellIDs,tempID);
            end
            
        end
    end
    cellsToAnalyze{e} = cellIDs; 
    fprintf('*');
    if mod(e,80) == 0
        fprintf('electrode %0.0f \n',e); 
    end
end
%% Gets the probabilities of activation of all the cells in the group to analyze against.
for p = 1:512
    cellIDs = cellsToAnalyze{p};
    patternNo = p; 
    [actCells,probs] = plotActivatedCells(fPath,eipath,...
        patternNo,'cellIDs',cellIDs,'probability',0.5,'plotsoff',1);
    allCells{p} = actCells;
    allprobs{p} = probs;
    disp(sprintf('finished pattern %0.0f\n',p));
end
%% For each cell, find the probability of activation
maxRatio = zeros(size(target_cell_group,2),1); 
for n= 1:1:size(target_cell_group,2)
    probability_by_electrode = zeros(512,1); 
    probabilityRatio_by_electrode = zeros(512,1); 
    for e=1:512
        if ismember(target_cell_group(n),allCells{e})
            probability = allprobs{e} ; 
            idx = find(target_cell_group(n)==allCells{e}); 
            probability_by_electrode(e) = probability(idx);
            try
                otherProb = max(probability(setdiff(1:length(probability),idx))); 
                if otherProb < 0.01; otherProb=0.01; end
                probabilityRatio_by_electrode(e) = probability_by_electrode(e)-otherProb; 
            catch
                probabilityRatio_by_electrode(e) = 5;
            end  
        end
    end
    maxRatio(n) = max(probabilityRatio_by_electrode); 
end

colors = jet(100);
colorScale = [-0.5 0.5];
figure; hold on;
for n = 1:1:size(target_cell_group,2)
    if maxRatio(n)>colorScale(2)
        color_idx = size(colors,1);
    else
        color_idx = ceil(size(colors,1)*(maxRatio(n)-colorScale(1))/diff(colorScale)); 
    end
    if color_idx == 0; color_idx = 1; end
    plot_rf_summaries(datarun, target_cell_group(n), 'clear', false, 'label', false,...
        'plot_fits', true, 'fit_color', 'k','fill_color',colors(color_idx,:));
end
axis image; axis off; 
colorbar; colormap jet; 
caxis(colorScale)
title('ON parasol activation probability - other cell activation probability'); 

%% Now for each cell, only look at the electrode with the maximum recorded EI signal

maxRatioEI = zeros(size(target_cell_group,2),1); 
for n= 1:1:size(target_cell_group,2)
    probability_by_electrode = zeros(512,1); 
    probabilityRatio_by_electrode = zeros(512,1); 
    tempID  = target_cell_group(n);          % gets the id number of the ith neuron     
    ei = squeeze(eiM(:,:,find(tempID == nIdList)))';
        [eiMin,minIdx] = min(min(ei));
    for e=minIdx
        if ismember(target_cell_group(n),allCells{e})
            probability = allprobs{e} ; 
            idx = find(target_cell_group(n)==allCells{e}); 
            probability_by_electrode(e) = probability(idx);
            try
                otherProb = max(probability(setdiff(1:length(probability),idx))); 
                if otherProb < 0.1; otherProb=0.1; end
                probabilityRatio_by_electrode(e) = probability_by_electrode(e)-otherProb; 
            catch
                probabilityRatio_by_electrode(e) = 5;
            end  
        end
    end
    maxRatioEI(n) = max(probabilityRatio_by_electrode); 
end

colors = jet(100);
colorScale = [-0.5 0.5];
figure; hold on;
for n = 1:1:size(target_cell_group,2)
    if maxRatioEI(n)>colorScale(2)
        color_idx = size(colors,1);
    else
        color_idx = ceil(size(colors,1)*(maxRatioEI(n)-colorScale(1))/diff(colorScale)); 
    end
    if color_idx == 0; color_idx = 1; end
    plot_rf_summaries(datarun, target_cell_group(n), 'clear', false, 'label', false,...
        'plot_fits', true, 'fit_color', 'k','fill_color',colors(color_idx,:));
end
axis image; axis off; 
colorbar; colormap jet; 
caxis(colorScale)
title('ON parasol activation probability - other cell activation probability EIs only'); 