% spatial sensitivity plots comparing algorithm and human

% % cells ids from /Volumes/Analysis/2015-04-14-0/data000/data000.ei
% cellIds = [1 392 633 917 1308 1699 2345 3049 3466 3947 4655 6873];
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
%%
figure(100); 
figure(200); scatter(xc,yc,3,'k','filled');    axis off; axis image;
line([0 0],[min(yc) max(yc)],'Color','k');
line([min(xc)/2 min(xc)/2],[min(yc) max(yc)],'Color','k');
line([max(xc)/2 max(xc)/2],[min(yc) max(yc)],'Color','k');
allthresholds = zeros(1,512); 
shape={'o','+','*','x','s','d','h'};
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
    eiContour_wPolyFit(eiAmps,'ei_thresh',ei_thresh,'figureNum',100,'plotCoords',plotCoords);
    eiContour_wPolyFit(eiAmps,'ei_thresh',ei_thresh,'figureNum',300,'plotCoords',plotCoords);
    % Find ei amps greater than a particular threshold
    elecs = find(eiAmps > ei_thresh); 
    figure(200); 
    hold on; scatter(xc(elecs),yc(elecs),eiAmps(elecs),'r','filled'); 
 
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
    figure(100); 
    hold on; scatter(xc(idx),yc(idx),150, thresholds(idx),'filled');
    
   allthresholds = allthresholds + thresholds; 
end
figure(100); 
colormap(cmap.m); caxis([0.5 4.5]);
colorbar;
title('50% activation thresholds, human sorting'); 
% 1. Check to make sure that the elecResp files for all the relevant
% electrodes are made & checked.
% 2. Make sure that the automated spike sorting has been done for each of
% the cells & patterns
% 3. Generate comparison plots. 
%% Show same plot where the bundle threshold is lower
load(['/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/'...
    'hand_sorted_bundle_thresholds/axonBundleThresholds_byPattern_2012_09_24_3_data001.mat']); 
idx = find(allthresholds);
idx_good = find(allthresholds < bundle_thresholds_2012_09_24');
idx_good = intersect(idx_good,idx);
figure(300); 
scatter(xc(idx_good),yc(idx_good),150, allthresholds(idx_good),'filled');
colormap(cmap.m); caxis([0.5 4.5]);
colorbar;
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
%% Unchecked code.
% subplot(1,3,2); line(0:4,0:4,'Color','k'); hold on; 
% scatter(thresh_hum, thresh_alg, 50,[0  .231  1], 'filled');
% xlabel('thresholds (\muA) using human analysis'); 
% ylabel('thresholds (\muA) using GM algorithm'); 
% title(sprintf('50%% activation thresholds, %0.0f stim patterns\n recording electrode = stimulating electrode',length(thresh_hum))); 
% 
% 
% % Calculate percent accuracy
% idx_fnD = find(thresh_humD - thresh_algD > 0.1);
% idx_fpD = find(thresh_algD - thresh_humD > 0.1);
% false_pos = length(find(thresh_humD(idx_fpD)==0));
% false_neg = length(find(thresh_algD(idx_fnD)==0));
% total_trials = length(thresh_humD);
% correct_trials = length(find(abs(thresh_humD - thresh_algD)<0.1));
% false_negD = false_neg/total_trials; 
% false_posD = false_pos/total_trials; 
% agreementD = correct_trials/total_trials; 
% 
% idx_fnS = find(thresh_hum - thresh_alg > 0.1);
% idx_fpS = find(thresh_alg - thresh_hum > 0.1);
% false_pos = length(find(thresh_hum(idx_fpS)==0));
% false_neg = length(find(thresh_alg(idx_fnS)==0));
% total_trials = length(thresh_hum);
% correct_trials = length(find(abs(thresh_hum - thresh_alg)<0.1));
% false_negS = false_neg/total_trials; 
% false_posS = false_pos/total_trials; 
% agreementS = correct_trials/total_trials; 
% 
% all_valsD = [agreementD; false_negD; false_posD ];
% all_valsS = [agreementS; false_negS; false_posS ];
% errorsD = [1.96*sqrt(1/total_trials*agreementD*(1-agreementD));...
%     1.96*sqrt(1/total_trials*false_negD*(1-false_negD));...
%     1.96*sqrt(1/total_trials*false_negD*(1-false_negD))];
% errorsS = [1.96*sqrt(1/total_trials*agreementS*(1-agreementS));...
%     1.96*sqrt(1/total_trials*false_negS*(1-false_negS));...
%     1.96*sqrt(1/total_trials*false_negS*(1-false_negS))];
% 
% subplot(1,3,3); 
% groupnames = {'agreement';'false negative';'false positive'};
% bw_colormap = [1  .271  0;
%     0  .231  1];
% gridstatus = 'y';
% bw_legend = {'rec elec != stim elec','rec elec = stim elec'}; 
% barweb([all_valsD all_valsS], [errorsD errorsS], [], ...
%     groupnames, 'Threshold analysis of algorithm v. human', ...
%     'error type', 'proportional agreement', bw_colormap, gridstatus, bw_legend);