%% Script with code snippets to generate the figures for 2016 NIH electrical stimulation grant.  

%% Figure showing bidirectional propagation + target eis. 

% see eifigs_2015_10_06_6.m

[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',197,'plotElecWaveforms',[238 230 231 232 259 264 276 105 93 86 83 79 76 72 64]);
[rawdata, amps0] = generateEiFromStimPattern('/Volumes/Analysis/2015-10-06-6/data001/',...
    90,'movieNo',221,'plotElecWaveforms',[238 230 231 232 259 264 276 105 93 86 83 79 76 72 64]);

% Contour plot of n1338
[eiM,neuronIdList] = convertEiFileToMatrix('/Volumes/Analysis/2015-10-06-6/data002/data002.ei');
ei = eiM(:,:,find(neuronIdList==1338)); 
eiamps = max(ei,[],2) - min(ei,[],2); 

eiContour_wLinFit(eiamps)
[xc,yc]=getElectrodeCoords512(); 
elec_waveforms = [230 231 232 259 264 93 83 79 76 72]; 
for ii = 1:length(elec_waveforms)
    text(xc(elec_waveforms(ii)), yc(elec_waveforms(ii)),...
        num2str(elec_waveforms(ii)),'HorizontalAlignment','center');
end
figure; 
for ii = 1:length(elec_waveforms); 
    subplot(2,5,ii)
    plot((0:size(ei,2)-1)/20,ei(elec_waveforms(ii),:)); 
    title([num2str(ii) ' : ' num2str(elec_waveforms(ii))]);
    xlim([0 1.5]); ylim([-100 100]);
end

% Contour plot of n1338 + other eis
[eiM,neuronIdList] = convertEiFileToMatrix('/Volumes/Analysis/2015-10-06-6/data002/data002.ei');
nlist = [1338 2988 4279]; 
plotCoords= true; 
for n = 1:length(nlist)
    ei = eiM(:,:,find(neuronIdList==nlist(n)));
    eiamps = max(ei,[],2) - min(ei,[],2);
    if n>1
        plotCoords = false; 
    end
    eiContour_wLinFit(eiamps,'figureNum',100,'showSoma',false,'plotCoords',plotCoords); 
end
 

%% Show the bidirectional voltage amps.

load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/Safe1193.mat');
idx = find(Safe1193); 
figure; 
scatter(xc(idx),yc(idx),500,Safe1193(idx),'filled'); 
axis image; axis off; colormap jet; colorbar; 
title('Safe zone 2015-11-09-3'); 

[xc,yc]=getElectrodeCoords512();
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_3.mat')
axonBundleThresholds_2015_10_6_3(isnan(axonBundleThresholds_2015_10_6_3)) = 4; 
zoomElecs = intersect(find(xc>-590),find(xc<200));
figure; scatter(xc(zoomElecs),yc(zoomElecs),100,axonBundleThresholds_2015_10_6_3(zoomElecs),'filled'); colormap(flipud(jet)); colorbar; axis image; axis off; 
figure; scatter(xc(zoomElecs),yc(zoomElecs),150./axonBundleThresholds_2015_10_6_3(zoomElecs),axonBundleThresholds_2015_10_6_3(zoomElecs),'filled'); colormap(flipud(jet)); colorbar; axis image; axis off; 

figure; scatter(xc,yc,10,'k'); hold on; scatter(xc(zoomElecs),yc(zoomElecs),20,'filled');
text(xc(135),yc(135),num2str(135),'HorizontalAlignment','center'); title('Electrodes of interest'); 
pathToAnalysisData = '/Volumes/Analysis/2015-10-06-3/data001-data002/'; 
%%
figure; hold on; 
axes_handle = gca; 
figure; axes_handle2 = gca; 
scatter(xc,yc,10,'k'); 
for p = 1:length(zoomElecs);
    patternNo = zoomElecs(p);
    movieNos = findMovieNos(pathToAnalysisData,patternNo); 
    [rawData, amplitudes] = generateEiFromStimPattern(pathToAnalysisData, patternNo,'movieNo',movieNos(end),'suppressPlots',true);
    %  polynomial fit
    axes(axes_handle2); cla(axes_handle2); scatter(xc,yc,amplitudes); axis image; axis off; 
    [curr_axon_x, curr_axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(amplitudes,'axonBundleFit',true);   
    axes(axes_handle2); hold on; plot(curr_axon_x,curr_axon_y,'k','LineWidth',2); 
    toplotornot = questdlg('plot?','plot','yes','no','yes'); 
    switch toplotornot
        case 'yes'
            axes(axes_handle); plot(curr_axon_x,curr_axon_y,'k','LineWidth',2); 
    end
    disp([num2str(p) ':' num2str(patternNo)]); 
end

axis image; axis off; 
text(xc(135),yc(135),num2str(135),'Color','b','HorizontalAlignment','center');


for ii = 103
    figure;
    patternNo = zoomElecs(ii);
    movieNos = findMovieNos(pathToAnalysisData,patternNo);
    [rawData, amplitudes] = generateEiFromStimPattern(pathToAnalysisData, patternNo,'movieNo',movieNos(end),'suppressPlots',true);
    [curr_axon_x, curr_axon_y, ~, soma_x, soma_y] = weighted_axon_poly_reg(amplitudes,'axonBundleFit',true);
    scatter(xc,yc,amplitudes); axis image; axis off;
    hold on; plot(curr_axon_x,curr_axon_y,'k','LineWidth',2);
end

 tubulin = imread('/Volumes/Data/2015-10-06-3/Imaging/confocal/MIP_tiled_hyperstack-tubulin.tif');
 smaller = tubulin(540:6310,2930:7288); clear tubulin; 
 
hsize = 100; % 241=60 um separation of electrodes.
sigma = 0.46;
h = fspecial('log',hsize,sigma);
filtered = imfilter(smaller(:,:,1),h);
figure; imshow(filtered); title(sprintf('laplacian of a gaussian hsize %0.0f sigma %0.2f',hsize,sigma)) ;

%% Figure FFF (A) RGC and bundle threshold comparison histogram for two 
% preps. (B) scatter plot showing lowest RGC somatic threshold vs. 
% bundle threshold across all electrodes, for two preps
% load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_3.mat')
% Preparation 1: 2015-11-09-3 (high density recording)
% load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/data_files/Safe1193.mat');
% fPath = '/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/';
% visionDataPath = '/Volumes/Analysis/2015-11-09-3/data000/data000'; 

% Preparation 2: 2015-10-06-3 (lower density recording)
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/hand_sorted_bundle_thresholds/axonBundleThresholds_2015_10_06_3.mat')
fPath = '/Volumes/Analysis/2015-10-06-3/data001-data002-autosort-template000/';
visionDataPath = '/Volumes/Analysis/2015-10-06-3/data000/data000'; 
safe = axonBundleThresholds_2015_10_6_3; 
 % Get a list of somatic activation for every electrode. 
disp('Find pattern-cell pairs to analyze (patternNo = stim electrode)..'); 
[eiM,nIdList] = convertEiFileToMatrix([visionDataPath '.ei']);
cellsToAnalyze = cell(512,1); 
sortingThreshold = 45; 
threshold = 10; 
for e = 1:512
    electrodeNo = e;
    cellIDs = [];
    for n = 1:1:size(nIdList,1)
        tempID  = nIdList(n);          % gets the id number of the ith neuron
        ei = squeeze(eiM(:,:,n))';
        eiMin = abs(min(ei));   
        if eiMin(electrodeNo) > threshold && (max(eiMin(:)) > sortingThreshold) 
            minIdx = find(ei(:,electrodeNo) == min(ei(:,electrodeNo)));
            max1 = max(ei(1:minIdx,electrodeNo));
            max2 = max(ei(minIdx:end,electrodeNo));
            if max1>2
%                 plot(ei(:,electrodeNo),'r'); title(['AXON max first phase is ' num2str(max1)])
            else
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
%% Save all the lowest thresholds by electode. 
ii = 1; 
for p = [1:512];
    patternNo = p;
    cellIDs = cellsToAnalyze{p};
    bundleThresh = safe(patternNo); 
    for n = 1:length(cellIDs)
        cellID = cellIDs(n); 
        fname = sprintf('elecRespAuto_n%0.0f_p%0.0f.mat',cellID,patternNo);
        try load(fullfile(fPath,fname))
            [thresh,~, ~]  = getThresholdFromAuto(elecRespAuto);
            allthresholds(ii) = thresh; 
            allbundlethresh(ii) = bundleThresh; 
            ii = ii+1; 
        catch
            disp([fullfile(fPath,fname) ' does not exist']);
        end
    end
    disp(sprintf('finished pattern %0.0f\n',p));
end
idx = find(allthresholds>0.1 & allthresholds<5); 
figure; scatter(allthresholds(idx),allbundlethresh(idx)/1.1,'k'); 
xlabel('somatic thresholds \muA'); ylabel('bundle safe zone \uA'); 
hold on; line(0:5,0:5,'Color','k','LineStyle','--');
xlim([0 4])
%%
% unpackcell ids
allcellIDs = [];
for p=1:512
    allcellIDs = [allcellIDs cellsToAnalyze{p}'];
end
[C,IA,IC] = unique(allcellIDs);
bestCellElecBundleThresh = zeros(size(C));
bestCellThresh=zeros(size(C));
for ii = 1:size(C,2)
    these_thresh = allthresholds(find(IC == ii));
    these_bundle_thresh = allbundlethresh(find(IC == ii));
    these_thresh(these_thresh<0) = 5; 
    these_thresh(these_thresh>5) = 5; 
    these_bundle_thresh(these_bundle_thresh==0) = 5; 
    these_bundle_thresh(isnan(these_bundle_thresh)) = 5; 
    thresh_diffs = these_bundle_thresh-these_thresh; 
    idx = find(max(thresh_diffs) == thresh_diffs,1)
    bestCellElecBundleThresh(ii) = these_bundle_thresh(idx); 
    bestCellThresh(ii) = these_thresh(idx); 
end
figure; scatter(bestCellThresh,bestCellElecBundleThresh/1.1,'k'); 
xlabel('somatic thresholds \muA'); ylabel('bundle thresholds \muA'); 
hold on; line(0:5,0:5,'Color','k','LineStyle','--'); 
ylim([0 4]); 
bestCellThresh(bestCellThresh==5) = []; 
figure; 
h1 = histfit(bestCellThresh,13,'kernel');
h1(1).FaceAlpha = 0.5;
h1(1).FaceColor = [1 0.2 0]; 
h1(2).Color = h1(1).FaceColor;
bestCellElecBundleThresh(bestCellElecBundleThresh==5)=[];
hold on;
h2 = histfit(bestCellElecBundleThresh,8,'kernel');
% sz = safe;
% sz(sz==0)=[]; 
%  h2 = histfit(sz,7,'kernel')
h2(1).FaceAlpha = 0.4;
h2(1).FaceColor = [0 0 0]; 
h2(2).Color = h2(1).FaceColor;
xlabel('current amplitude \muA'); 
ylabel('number of cells measured')
legend([h1(2) h2(2)],'soma threshold','safe zone'); 
title('2015-11-09-3 or 2015-10-06-3'); 

%% Figure HHH
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/AlgorithmOutputs_2015_10_6_6.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/HandAnalysis_2015_10_6_6.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/Electrodes_2015_10_6_6.mat')
[xc,yc] = getElectrodeCoords512(); 
figure; 
subplot(2,1,1); 
scatter(xc(Electrodes_2015_10_6_6),yc(Electrodes_2015_10_6_6),200./AlgorithmOutputs_2015_10_6_6,AlgorithmOutputs_2015_10_6_6,'filled'); 
colormap(flipud(jet)); colorbar; title('Algorithm Output 2015-10-06-6'); 
caxis([0.5 2.5]); axis image; axis off; 
subplot(2,1,2); 
scatter(xc(Electrodes_2015_10_6_6),yc(Electrodes_2015_10_6_6),200./HandAnalysis_2015_10_6_6,HandAnalysis_2015_10_6_6,'filled'); 
colormap(flipud(jet)); colorbar; title('Hand Analysis 2015-10-6-6 '); 
caxis([0.5 2.5]); axis image; axis off; 

load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/AlgorithmOutputs_2015_10_6_3.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/HandAnalysis_2015_10_6_3.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/Electrodes_2015_10_6_3.mat')
[xc,yc] = getElectrodeCoords512(); 
figure; 
subplot(2,1,1); 
scatter(xc(Electrodes_2015_10_6_3),yc(Electrodes_2015_10_6_3),200./AlgorithmOutputs_2015_10_6_3,AlgorithmOutputs_2015_10_6_3,'filled'); 
colormap(flipud(jet)); colorbar; title('Algorithm Output 2015-10-06-3'); 
caxis([0.5 2.5]); axis image; axis off; 
subplot(2,1,2); 
scatter(xc(Electrodes_2015_10_6_3),yc(Electrodes_2015_10_6_3),200./HandAnalysis_2015_10_6_3,HandAnalysis_2015_10_6_3,'filled'); 
colormap(flipud(jet)); colorbar; title('Hand Analysis 2015-10-6-3 '); 
caxis([0.5 2.5]); axis image; axis off; 

load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/AlgorithmOutputs_2015_09_23_2.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/HandAnalysis_2015_09_23_2.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/Electrodes_2015_09_23_2.mat')
[xc,yc] = getElectrodeCoords512(); 
figure; 
subplot(2,1,1); 
scatter(xc(Electrodes_2015_09_23_2),yc(Electrodes_2015_09_23_2),200./...
    AlgorithmOutputs_2015_09_23_2,AlgorithmOutputs_2015_09_23_2,'filled'); 
colormap(flipud(jet)); colorbar; title('Algorithm Output 2015-09-23-2'); 
caxis([0.2 2]); axis image; axis off; 
subplot(2,1,2); 
scatter(xc(Electrodes_2015_09_23_2),yc(Electrodes_2015_09_23_2),200./...
    HandAnalysis_2015_09_23_2,HandAnalysis_2015_09_23_2,'filled'); 
colormap(flipud(jet)); colorbar; title('Hand Analysis 2015-09-23-2 '); 
caxis([0.2 2]); axis image; axis off; 

load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/AlgorithmOutputs_2012_09_24_3.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/HandAnalysis_2012_09_24_3.mat')
load('/Volumes/Lab/Projects/electrical_stim/bundle-safe-zones/AlgorithmOutputs/Electrodes_2012_09_24_3.mat')
[xc,yc] = getElectrodeCoords512(); 
figure; 
subplot(2,1,1); 
scatter(xc(Electrodes_2012_09_24_3),yc(Electrodes_2012_09_24_3),200./...
    AlgorithmOutputs_2012_09_24_3,AlgorithmOutputs_2012_09_24_3,'filled'); 
colormap(flipud(jet)); colorbar; title('Algorithm Output 2012-09-24-3'); 
caxis([0.8 2.5]); axis image; axis off; 
subplot(2,1,2); 
scatter(xc(Electrodes_2012_09_24_3),yc(Electrodes_2012_09_24_3),200./...
    HandAnalysis_2012_09_24_3,HandAnalysis_2012_09_24_3,'filled'); 
colormap(flipud(jet)); colorbar; title('Hand Analysis 2012-09-24-3 '); 
caxis([0.8 2.5]); axis image; axis off; 

%% Figure LLL response curves for two cells, manual vs automated. 
load('/Volumes/Lab/Users/grosberg/Gonzalo_sorting_test_data/OutputAll.mat')
kk = 551; 
alg_output = OutputAll(kk);
gm_path = alg_output.path;
ii = find(gm_path==filesep,3,'last');
my_path = fullfile('/Volumes/Analysis/',gm_path(ii(1):ii(3)));
alg_output.path = my_path;
compareAlgHum_dataplots(alg_output, 1);

% For complete scatter plots, see spikesorting_alg_vs_human_710ampseries.m
%% Figure TTT Raphe
raphe_plots(); 
 %% Figure MMM Correlation between EI amplitude and sensitivity
analysis_eiAmps_vs_actThresh()

%% Dancers figure. 
dataPath = '/Volumes/Analysis/2012-09-24-3/data000/data000'; 
cellIds = [3457 4058 4656 2796 5123 5748 6439];
 % Analyze all electrodes that recorded voltages over some threshold. 
% Load EI. 
datarun  = load_data(dataPath);
datarun  = load_neurons(datarun);
datarun  = load_sta(datarun, 'load_sta', 'all');
datarun  = load_params(datarun);
datarun  = load_ei(datarun, cellIds,'keep_java_ei','false');
figure; 
for n = 1:length(cellIds)
    neuronId = cellIds(n);
    ei = datarun.ei.eis{get_cell_indices(datarun,neuronId)};
    eiAmps = max(ei,[],2) - min(ei,[],2);
    [curve_x, curve_y, ~, soma_x, soma_y, ~, ~] = weighted_axon_poly_reg(eiAmps);
    hold on; plot(curve_x,curve_y,'k'); 
end

%% Figure 9 using EIs to predict the best stimulating electrode. 

% Get a list of ON parasol cells from 2015-11-09-3. 
fPath = '/Volumes/Analysis/2015-11-09-3/data001-data002-autosort/';
visionDataPath = '/Volumes/Analysis/2015-11-09-3/data000/data000'; 
eipath = [visionDataPath '.ei']; 
% fPath = '/Volumes/Analysis/2015-10-06-3/data001-data002-autosort-template000/';
% visionDataPath = '/Volumes/Analysis/2015-10-06-3/data000/data000'; 

datarun = load_data(visionDataPath);  
datarun = load_params(datarun); 
datarun = load_sta(datarun); 
datarun = load_ei(datarun,'all'); 
sortingThreshold = 40; 
fillcolors = [0.275 0.094, 0.471;
    0.471 0.094 0.29;
    0.29 0.471 0.094;
    0.094, 0.471, 0.275];

celltype = datarun.cell_types{1}.name; 
cell_ids = datarun.cell_types{1}.cell_ids; 
 [allVals, actProb, actAmp] = genActThreshSpatialMaps_auto(fPath,eipath,cell_ids); 