function analysis_eiAmps_vs_actThresh()
% The purpose of this analysis is to understand whether the EI can be used
% to predict the best electrodes for electrical stimulation.
%%
% Load the plot showing spatial sensitivity of seven cells in a single prep
[xc,yc] = getElectrodeCoords512(); 
codebase_path = matlab_code_path; 
cmap = load([codebase_path 'code/projects/electrical_stim/'...
    'resources/redtealcmap.mat']); 

% cells ids from /Volumes/Analysis/2012-09-24-3/data000/data000.ei
dataPath = '/Volumes/Analysis/2012-09-24-3/data000/data000'; 
cellIds = [3457 4058 4656 2796 5123 5748 6439];

% Analyze all electrodes that recorded voltages over some threshold. 
% Load EI. 
datarun  = load_data(dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, cellIds,'keep_java_ei','false');
ei_thresh = 10; 
alleis = zeros(512,70); 
figure(100); 
allthresholds = zeros(1,512); 
shape={'o','+','*','x','s','d','h'};
colors = lines(length(cellIds)); 
%%
for n = 1:length(cellIds)
    thresh_quad1 = zeros(1,512);
    thresh_quad2 = zeros(1,512);
    thresh_quad3 = zeros(1,512);
    thresh_quad4 = zeros(1,512);
    neuronId = cellIds(n);
    ei = datarun.ei.eis{get_cell_indices(datarun,neuronId)};
    alleis = alleis + ei(:,1:70); 
    eiAmps = max(ei,[],2) - min(ei,[],2);
    if n == 1
        plotCoords = true;
    else
        plotCoords = false;
    end; 
    eiContour_wLinFit(eiAmps,'linFitThresh',5,'figureNum',100,...
        'plotCoords',plotCoords);
    
    % Find ei amps greater than a particular threshold
    elecs = find(eiAmps > ei_thresh); 
   
    % Determine the quadrant for each electrode in a given cell .  
    if ~isempty(find(xc(elecs) < min(xc/2), 1))
        fpath = '/Volumes/Analysis/2012-09-24-3/data003/';
        [~, thresh_quad1] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1,'suppressPlots',1);
    end

    if ~isempty(find(xc(elecs) > min(xc/2) & xc(elecs) < 0, 1));
        fpath = '/Volumes/Analysis/2012-09-24-3/data004/';
        [~, thresh_quad2] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1,'suppressPlots',1);
    end
    
    if ~isempty(find(xc(elecs) < max(xc/2) & xc(elecs) > 0, 1))
        fpath = '/Volumes/Analysis/2012-09-24-3/data005/';
        [~, thresh_quad3] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1, 'suppressPlots',1);

    end
    if ~isempty(find(xc(elecs) > max(xc/2), 1))
        fpath = '/Volumes/Analysis/2012-09-24-3/data006/';
        [~, thresh_quad4] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1, 'suppressPlots',1);
    end
    
    thresholds = thresh_quad1 + thresh_quad2 + thresh_quad3 + thresh_quad4;
    thresholds(thresholds>6) = 0;
    thresholds(thresholds<0) = 0;
    idx = find(thresholds);
    % Plot the EI amplitude vs the activation threshold
    figure(600);
    hold on; scatter(eiAmps(idx),thresholds(idx),50,'marker',shape{n});
    [sortedAmps, sortingOrder] = sort(eiAmps(idx));
    thresh_subset = thresholds(idx);
    thresh_subset = thresh_subset(sortingOrder);
    stopIdx = find(sortedAmps<26,1,'last');
    % Calculate regression.
    X = [ones(length(sortedAmps(1:stopIdx)),1) sortedAmps(1:stopIdx)];
    b = X\thresh_subset(1:stopIdx)';
    yCalc1 = X*b;
    figure(500); 
    hold on;
    plot(sortedAmps(1:stopIdx),yCalc1,'Color',colors(n,:));
    plot(sortedAmps,thresh_subset,[shape{n}],'Color',colors(n,:),...
        'MarkerSize',16);
    xlabel('EI amplitude (mV)'); ylabel('activation threshold (uA)');
    figure(100);
    hold on; scatter(xc(idx),yc(idx),150, thresholds(idx),'filled');
    allthresholds = allthresholds + thresholds;
end
figure(100); 
colormap(cmap.m); caxis([0.5 4.5]);
colorbar;
title('50% activation thresholds, human sorting'); 



%% Plot all the thresholds grouped as somatic or axonal EI spike. 
figure(10); 
for ii = 1:512
    if allthresholds(ii)
        figure(10); cla; 
        plot(alleis(ii,:)); 
        figure(20);  hold on; 
        h=scatter(max(alleis(ii,:))-min(alleis(ii,:)),allthresholds(ii),100,'filled');
        ButtonName = questdlg('What part of the cell?', ...
            'Spatial location', ...
            'axon', 'soma', 'axon');
        switch ButtonName,
            case 'axon',
                h.CData = [0 1 1];
                h.MarkerEdgeColor = [0 0 0]; 
                h.Marker = '+';
            case 'soma',
                h.CData = [0 0 1]; 
        end % switch
    end
end
xlabel('EI amplitude (uV)')
ylabel('actvation thresholds (uA)')
title(sprintf(['ei amplitude vs activation thresholds grouped by \n'...
    'axonal stim (black +) and somatic stim (blue o)'])); 

