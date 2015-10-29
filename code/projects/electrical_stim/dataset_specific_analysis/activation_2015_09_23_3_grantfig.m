% Plot a group of cells from 2015-09-23-3

dataPath =  '/Volumes/Analysis/2015-09-23-3/data005/data005';
datarun = load_data(dataPath);
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');
positions = datarun.ei.position;

ONPcellIDs = [2662,2808,2178]; 
OFFPcellIDs = [2417,2332,2791]; 
ONMcellIDs = [2776,2402,2446,2536];
OFFMcellIDs = [2672,2538];
largecellID = 2765; 
figure; 
h1= plot_rf_summaries(datarun, ONPcellIDs, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b');
h2= plot_rf_summaries(datarun, OFFPcellIDs, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k');
h3= plot_rf_summaries(datarun, ONMcellIDs, 'clear', false, 'label', true, 'label_color', 'm', 'plot_fits', true, 'fit_color', 'm');
h4= plot_rf_summaries(datarun, OFFMcellIDs, 'clear', false, 'label', true, 'label_color', 'g', 'plot_fits', true, 'fit_color', 'g');
h5= plot_rf_summaries(datarun, largecellID, 'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y');
axis image; axis off;
    % legend([h1 h2 h3 h4 h5],'ON parasol','OFF parasol','ON midget','OFF midget','large')
title(sprintf('ON parasol (blue) OFF parasol (black)\n ON midget (pink) OFF midget (green)\n stimulating electrode %0.0f',electrodeNo));
 
elecStimPath = '/Volumes/Analysis/2015-04-09-2/data003/';


%% Visualize EIs
[xc,yc] = getElectrodeCoords512(); 
codebase_path = matlab_code_path; 
cmap = load([codebase_path 'code/projects/electrical_stim/'...
    'resources/redtealcmap.mat']); 

cellIds = [ONPcellIDs OFFPcellIDs ONMcellIDs OFFMcellIDs largecellID]; 

ei_thresh = 10; 
%%
pats = [185   186   187   177   178   179   180   169   170   171   161   162   163   164   153   154   155   156];
figure(400); scatter(xc(pats),yc(pats),10,'k','filled'); axis image; 
figure(100); 
figure(200); scatter(xc,yc,3,'k','filled');    axis off; axis image;

for n = 1:length(cellIds)
    thresh_quad1 = zeros(1,512);
    thresh_quad2 = zeros(1,512);
    thresh_quad3 = zeros(1,512);
    thresh_quad4 = zeros(1,512);
    neuronId = cellIds(n);
    ei = datarun.ei.eis{get_cell_indices(datarun,neuronId)};
    eiAmps = max(ei,[],2) - min(ei,[],2);
    if n == 1
        plotCoords = true;
    else
        plotCoords = false;
    end; 
     
    eiContour_wLinFit(eiAmps,'linFitThresh',5,'figureNum',0,'plotCoords',plotCoords,'numContourLevels',2);
    hold on; scatter(xc,yc,eiAmps,'filled'); title(num2str(neuronId)); 
    eiContour_wLinFit(eiAmps,'linFitThresh',5,'figureNum',300,'plotCoords',plotCoords,'numContourLevels',2);
    % Find ei amps greater than a particular threshold
    elecs = find(eiAmps > ei_thresh); 
    figure(200); 
    hold on; scatter(xc(elecs),yc(elecs),eiAmps(elecs),'r','filled'); 
 
    % Determine the quadrant for each electrode in a given cell .  
    
    if ~isempty(find(xc(elecs) < min(xc/2), 1))
        fpath = '/Volumes/Analysis/2015-09-23-3/data001/';
        [~, thresh_quad1] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1,'suppressPlots',1);
    end

    if ~isempty(find(xc(elecs) > min(xc/2) & xc(elecs) < 0, 1));
        fpath = '/Volumes/Analysis/2015-09-23-3/data002/';
        [~, thresh_quad2] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1,'suppressPlots',1);
    end
    
    if ~isempty(find(xc(elecs) < max(xc/2) & xc(elecs) > 0, 1))
        fpath = '/Volumes/Analysis/2015-09-23-3/data003/';
        [~, thresh_quad3] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1, 'suppressPlots',1);

    end
    if ~isempty(find(xc(elecs) > max(xc/2), 1))
        fpath = '/Volumes/Analysis/2015-09-23-3/data004/';
        [~, thresh_quad4] = genActThreshSpatialMaps(fpath, neuronId, ...
            'dispCheckedElecs',1, 'suppressPlots',1);

    end
   
    thresholds = thresh_quad1 + thresh_quad2 + thresh_quad3 + thresh_quad4;
    
    thresholds(thresholds>6) = 0; 
    thresholds(thresholds<0) = 0; 
    idx = find(thresholds);
    figure(100); 
    hold on; scatter(xc(idx),yc(idx),150, thresholds(idx),'filled');
  
    
end
colormap(cmap.m); caxis([0.5 4.5]);
colorbar;
title('50% activation thresholds, human sorting'); 
%%
% 1. Check to make sure that the elecResp files for all the relevant
% electrodes are made & checked.
% 2. Make sure that the automated spike sorting has been done for each of
% the cells & patterns
% 3. Generate comparison plots. 

%%
if ~exist('spatialSensitivity','var')
    spatialSensAlg = load(['/Volumes/Lab/Projects/electrical_stim/GM-sorting-'...
        'validation/2015-06-algorithm-results/spatialSensitivity']);
    spatialSensitivity = spatialSensAlg.Output; clear spatialSensAlg;
end
%%

patternNos = zeros(size(spatialSensitivity)); 
thresh_hum = zeros(size(spatialSensitivity)); 
thresh_alg = zeros(size(spatialSensitivity)); 
algThreshElectrodes = zeros(size(spatialSensitivity)); 
for p=1:size(spatialSensitivity,2)
    
    alg_output = spatialSensitivity(p);
    [thresholdHum, thresholdAlg] = fitToErfOutputAndHuman(alg_output);
    patternNos(p) = alg_output.stimInfo.patternNo;
    algThreshElectrodes(p)= alg_output.stimInfo.patternNo;
    if (thresholdHum > 4.5) || (thresholdHum < 0)
         thresh_hum(p) = 0; 
    else
        thresh_hum(p) = thresholdHum;
    end
    if (thresholdAlg > 4.5) || (thresholdAlg < 0)
        thresh_alg(p) = 0;
    else
        thresh_alg(p) = thresholdAlg;
    end
    if p<100 && p>50
        % Display results if the patterns are not well matched.
        if abs(thresh_hum(p) - thresh_alg(p)) > 0.1
            displayResults_compareOutput(alg_output,1,0)
            suptitle([alg_output.path ' n:' num2str(alg_output.neuronInfo.neuronIds) ...
                ' pattern: ' num2str(alg_output.stimInfo.patternNo)]);
        end
    end
end


figure; %subplot(1,3,1); 
line(0:4,0:4,'Color','k'); hold on; 
scatter(thresh_hum, thresh_alg, 50, [1  .271  0],'filled');
xlabel('thresholds (\muA) using human analysis'); 
ylabel('thresholds (\muA) using GM algorithm'); 
title(sprintf('50%% activation thresholds, %0.0f stim patterns',length(thresh_hum))); 

%% Display threshold contour plot based on 
figure(300); hold on; 
idx = find(thresh_alg); 
scatter(xc(algThreshElectrodes(idx)),yc(algThreshElectrodes(idx)),150,...
    thresh_alg(idx),'filled');
colormap(cmap.m); caxis([0.5 4.5]);
colorbar;
title('50% activation thresholds, algorithm sorting'); 
