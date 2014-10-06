
clear all


%same for 2011-06-24-5, 2011-07-14-0 and 2013-05-28-3
startPadding = 0.15;  endPadding = 0.05; 
barColor = 'white';

% (subset of) options related to algorithm
% run_opt.trial_num = 80; % > 0 %example trial for plotting
% run_opt.vel = 0:1:500; % velocity values to calculate motion signal for when generating motion signal curves
%
% run_opt.cell_type = 'ON parasol'; % on/off parasol, on/off midget

% general use parameters
%bar_speed = 240; %in pixels/second



%%%% dataset info for 2011-06-24-5
if 1
    %analysisPathBase = '/Analysis/2011-06-24-5/';
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/';
    
    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_mb = [2 33 183 332 214 783]; %from data001-from-data000.neurons
    cell_ids_er = [1 33 183 332 214 784]; %from data006.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat
    cell_ids_rf = [1 32 229 396 531 737]; %from rf_dataset
    
    % electrical stimulation parameters.
    movieChunk = 30;
    EstimData = 'data015';

    visChosenBoth = [80 77]; % which of the visual response trials to replicate
    visChosen = visChosenBoth(1);
    
    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data001-from-data000/data001-from-data000.neurons'];
end

%%%% dataset info for 2011-07-14-0
if 0
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/';
    
    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_mb = [213 200 348 467 543 663 783 901]; %from mbVisNeuronsFile
    cell_ids_er = [48  187 304 407 498 618 783 918]; %from data005.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat
    cell_ids_rf = [47  187 303 473 498 713 784 903]; %from rf_dataset
    
    % electrical stimulation parameters.
    movieChunk = 33; %30/33 or 37/40
    EstimData = 'data013';
    
    visChosenBoth = [1 48]; %1 (m30, m37) or 48 (m33, m40); which of the visual response trials to replicate
    visChosen = visChosenBoth(2);
    
    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data001/data001.neurons'];
end

%%%% dataset info for 2013-05-28-3/data012 (gray background)
if 0
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/';
    
    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_mb = [76  138 377 437 783 917]; %from mbVisNeuronsFile
    cell_ids_er = [166 138 349 437 784 916]; %from data010.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat
    cell_ids_rf = [167 137 350 439 722 932]; %from rf_dataset (closest in time to MB vis stim)
    
    % electrical stimulation parameters.
    movieChunk = 41;
    EstimData = 'data012';
    
    visChosenBoth = [46 53];% which of the visual response trials to replicate
    visChosen = visChosenBoth(1);
    
    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data003-gw/data003-gw.neurons'];
end

%%%% dataset info for 2013-05-28-3/data013 (no light)
if 0
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/';
    
    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_mb = [76  138 377 437 783 917]; %from mbVisNeuronsFile
    cell_ids_er = [166 138 349 437 784 916]; %from data010.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat
    cell_ids_rf = [167 137 350 439 722 932]; %from rf_dataset (closest in time to MB vis stim)
    
    % electrical stimulation parameters.
    movieChunk = 41;
    EstimData = 'data013';
    
    visChosenBoth = [46 53];% which of the visual response trials to replicate
    visChosen = visChosenBoth(1);
    
    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data003-gw/data003-gw.neurons'];
end


nCells = length(cell_ids_mb);

% for labeling figures
date_piece = analysisPathBase(strfind(analysisPathBase,'/Analysis/')+10:end-1);



%% load analysis results from file

output_save_path_estim = [analysisPathBase EstimData filesep 'm' num2str(movieChunk) '_mb_decode_out.mat'];
output_save_path_vis = [analysisPathBase EstimData filesep 'vis_mb_decode_out.mat'];

try %load from server
    load(output_save_path_estim)
catch
    load(output_save_path_estim(strfind(output_save_path_estim,'/Analysis'):end))
end
mb_decode_estim = mb_decode_out; clear mb_decode_out


try %load from server
    load(output_save_path_vis)
catch
    load(output_save_path_vis(strfind(output_save_path_vis,'/Analysis'):end))
end
mb_decode_vis = mb_decode_out; clear mb_decode_out


%% load spike times (in response to electrical and visual stimulation)

bar_speed = mb_decode_estim.bar_speed;
spikeTimesAll_estim = mb_decode_estim.spikeTimesAll;
spikeTimesAll_vis =   mb_decode_vis.spikeTimesAll;

nTrials_estim = length(spikeTimesAll_estim);
nTrials_vis = length(spikeTimesAll_vis);

cell_x_pos = mb_decode_estim.cell_x_pos; %should be same for estim and vis
[~, cellOrder] = sort(cell_x_pos);


run_opt = mb_decode_estim.run_opt; %should be the same for vis and estim

%% load mbSpikes directly to get electrically-elicited flags

try
    load([analysisPathBase EstimData filesep 'm' num2str(movieChunk) '_mbSpikes.mat']);
catch
    load([analysisPathBase(strfind(analysisPathBase,'/Analysis'):end) EstimData filesep 'm' num2str(movieChunk) '_mbSpikes.mat'])
end

elicitedBin = cell(nTrials_estim,1);
%sort moving bar spikes according to cell_ids_er (so that cell_x_pos indeces match up);
for ii = 1:nCells
    cellInd = find([mbSpikes.cellIDVis] == cell_ids_er(ii));
    mbSpikesAll(ii) = mbSpikes(cellInd); %#ok<SAGROW>
    for jj = 1:nTrials_estim
        elicitedBin{jj}{ii} = mbSpikesAll(ii).elicited{jj};
    end
end

clear mbSpikesAll mbSpikes

%% load pulse times

try %load from server
    pattern_application_times = extract_stim_occurrence_timing(EstimData(end-2:end), analysisPathBase, nCells, movieChunk);
catch
    pattern_application_times = extract_stim_occurrence_timing(EstimData(end-2:end),...
        analysisPathBase(strfind(analysisPathBase,'/Analysis'):end), nCells, movieChunk);
end

%% get length of moving bar -- in future this should be extracted from MovingBarVisSpikes structure

try %load from server
    [~, intervals] = loadMovingBarResponses(mbVisNeuronsFile, cell_ids_mb, barColor, 'trialPadding', [startPadding endPadding]);
catch
    [~, intervals] = loadMovingBarResponses(mbVisNeuronsFile(strfind(mbVisNeuronsFile,'/Analysis'):end), cell_ids_mb, barColor, 'trialPadding', [startPadding endPadding]);
end


%% verify that visChosenBoth is correct

trial_length = mean(cellfun(@diff, intervals)); %in seconds

% check for concistency with visChosen
[~, order] = spike_train_median_calculator(spikeTimesAll_vis, 2, 1, 10, 'both', trial_length, 'cellIDs', cell_ids_mb);

if ~isequal(order(1:2), visChosenBoth)
    keyboard
else
    disp('all good')
end


%% figure 1 panel a: visual stim with mosaic

% bar vs. rf positions verified by comparison with Vision display

datarun = mb_decode_vis.datarun;

datarun.stimulus.x_start = 200;
datarun.stimulus.x_end = 440;
datarun.stimulus.y_start = 120;
datarun.stimulus.y_end = 360;

figure; hold on
plot_rf_summaries(datarun, cell_ids_rf, 'coordinates', 'monitor');

fill([168 168 200 200], [80 400 400 80], 'k'); %initial bar position

title(date_piece)
%plot([200 440 440 200 200], [120 120 360 360 120]) %outline of WN region, for checking RF positions against Vision display


%% figure 1 panel b: target visually-elicited spike sequence and single
% (median-distance) electrical stim response

%determine which estim response is median-distance from target visual
%response

costParam = 10; %in ms
costParam = costParam/1000; %in seconds; the distance that a spike moves that is equivalent (in cost) to deleting and inserting spike
costParam = 2/costParam; %cost per second to move spike

targetTrain = spikeTimesAll_vis{visChosen};

costs = zeros(1,nTrials_estim);
for ii = 1:nTrials_estim %trial being compared
    pairCost = 0;
    for kk = 1:nCells %calculate pairwise cost for each cell and sum
        d = spkd(spikeTimesAll_estim{ii}{kk}, targetTrain{kk}, costParam);
        
        % normType = 'both'
        thisCost = d/(length(spikeTimesAll_estim{ii}{kk}) + length(targetTrain{kk})); %normalize by total number of spikes in both spike trains
        
        pairCost = pairCost + thisCost; %sum over cells for ith vs. jth trains
    end
    costs(ii) = pairCost;
end

costSorted = sort(costs);
medianCost = costSorted(find(costSorted >= median(costs),1));
medianTrial = find(costs == medianCost);

clear costSorted medianCost


%make figure

spotSize = 3;

figure
for ii = 1:nCells
    axes('position', [0.1 (nCells - ii + 1)/(nCells+2) 0.8 1/(nCells+3)])
    hold on
    cellInd = cellOrder(ii);
    
    % target visual spike train
    plot([0 trial_length], [2 2], '-', 'color', [0.8 0.8 0.8])
    for jj = 1:length(targetTrain{cellInd})
        plot(targetTrain{cellInd}(jj), 2, 'o', 'MarkerSize', spotSize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0 0])
    end
    
%     % pulse application times
%     for jj = 1:length(pattern_application_times{cellInd,2})
%         plot(pattern_application_times{cellInd,2}(jj)/20000, 1.5, 'o', 'MarkerSize', spotSize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0 0 1])
%     end
%     
    % median distance electrically elicited spike train
    for jj = 1:length(spikeTimesAll_estim{medianTrial}{cellInd})
        if elicitedBin{medianTrial}{cellInd}(jj) %electrically elicited spike
            plot(spikeTimesAll_estim{medianTrial}{cellInd}(jj), 1, 'o', 'MarkerSize', spotSize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [1 0 0])
        else %spontaneous spike or spike from cross-talk
            plot(spikeTimesAll_estim{medianTrial}{cellInd}(jj), 1, 'o', 'MarkerSize', spotSize, 'MarkerEdgeColor', 'none', 'MarkerFaceColor', [0.5 0.5 0.5])
        end
    end
    
    set(gca, 'xlim', [0 trial_length], 'ylim', [0 3], 'ytick', [])
    if ii == 1
        title([date_piece filesep EstimData  ', vis trial ' num2str(visChosen) ', estim trial ' num2str(medianTrial) ', m' num2str(movieChunk)])
    end
    if ii == nCells
        xlabel('time (s)')
    else
        set(gca, 'xtick', [])
    end
    ylabel(['cell' num2str(cell_ids_er(cellInd))])
    hold off
end


%% plot rasters of visually and electrically elicited spikes

figure

%axes for title
axes('position', [0.1 (nCells+1)/(nCells+2) 0.8 1/(nCells+20)])
title([date_piece filesep EstimData  ', vis trial ' num2str(visChosen) ', m' num2str(movieChunk)])
axis off

for i = 1:nCells    
    %visual stimulation
    axes('position', [0.1 (nCells - i + 1)/(nCells+2) 0.35 1/(nCells+3)])
    hold on
    cellInd = cellOrder(i);
    for k = 1:nTrials_vis
        for j = 1:length(spikeTimesAll_vis{k}{cellInd})
            if k == visChosen; %highlight trial being recreated
                plot(spikeTimesAll_vis{k}{cellInd}(j), k-0.5, '.', 'markerFaceColor', 'none', 'markerEdgeColor', [0 0 0])
            else
                plot(spikeTimesAll_vis{k}{cellInd}(j)*[1 1], [k-1 k], '-', 'LineWidth', 1, 'color', [0 0 0])
            end
        end
    end
    hold off
    
    set(gca, 'yLim', [0 nTrials_vis], 'xlim', [0 trial_length])
    if i == nCells
        xlabel('time (s)')
    else
        set(gca, 'xtick', [])
    end
    ylabel(['cell' num2str(cell_ids_mb(cellInd))])
    
    % electrical stimulation
    axes('position', [0.55 (nCells - i + 1)/(nCells+2) 0.35 1/(nCells+3)])
    hold on
    cellInd = cellOrder(i);
    for k = 1:nTrials_estim
        for j = 1:length(spikeTimesAll_estim{k}{cellInd})
            plot(spikeTimesAll_estim{k}{cellInd}(j)*[1 1], [k-1 k], '-', 'LineWidth', 1, 'color', [1 0 0])
        end
    end
    hold off
    
    set(gca, 'yLim', [0 nTrials_estim], 'xlim', [0 trial_length])
    if i == nCells
        xlabel('time (s)')
    else
        set(gca, 'xtick', [])
    end
    ylabel(['cell' num2str(cell_ids_mb(cellInd))])    
end


%% plot motion signal curve for a single trial
% 
% sig_strength_curve_1_trial = sig_strength_curve{run_opt.trial_num};
% vel_est = (run_opt.vel(sig_strength_curve_1_trial==max(sig_strength_curve_1_trial)));
% 
% figure; hold on
% plot(bar_speed*[1 1], [0 max(sig_strength_curve_1_trial)*1.2], 'k--')
% plot(run_opt.vel, sig_strength_curve_1_trial)
% xlabel('speed tuning')
% ylabel('signal strength')
% 
% clear sig_strength_curve
% 

%% motion signal curves and estimated velocities


vel_ests_brute_vis   = mb_decode_vis.vel_ests_brute;
vel_ests_brute_estim = mb_decode_estim.vel_ests_brute;

vel_ests_eff_vis     = mb_decode_vis.vel_ests_eff;
vel_ests_eff_estim   = mb_decode_estim.vel_ests_eff;

sig_strength_curve_vis   = mb_decode_vis.sig_strength_curve;
sig_strength_curve_estim = mb_decode_estim.sig_strength_curve;


%% plot motion signal curves on top of each other

figure; hold on
plot(bar_speed*[1 1], [0 max(sig_strength_curve_vis{jj})*1.5], 'k--')
for jj = 1:nTrials_vis
    plot(run_opt.vel, sig_strength_curve_vis{jj}, 'color', [0 0 0])
    plot(vel_ests_brute_vis(jj), max(sig_strength_curve_vis{jj}), '*', 'markerEdgeColor', [0 0 0]);
    if visChosen == jj
        plot(run_opt.vel, sig_strength_curve_vis{jj}, 'color', [0 0 1], 'linewidth', 2)
    end
    %plot(vel_ests_brute_vis(jj)*[1 1], [0 max(sig_strength_curve_vis{jj})], 'k-')
end

for jj = 1:nTrials_estim
    plot(run_opt.vel, sig_strength_curve_estim{jj}, 'color', [1 0 0])
    plot(vel_ests_brute_estim(jj), max(sig_strength_curve_estim{jj}), '*', 'markerEdgeColor', [1 0 0]);
    %plot(vel_ests_brute_estim(jj)*[1 1], [0 max(sig_strength_curve_estim{jj})], 'r-')
end

xlabel('speed (pixels/second)')
ylabel('net motion signal')
title(['speed tuning curves' 10 date_piece filesep EstimData  ', vis trial ' num2str(visChosen) ', m' num2str(movieChunk)])

%% plot histograms of velocity estimates


histCenters = 150:300;
histBinWidth = diff(histCenters(1:2)); %assumes uniform bin widths
histCounts_vis   = hist(vel_ests_brute_vis, histCenters);
histCounts_estim = hist(vel_ests_brute_estim, histCenters);

hist_vis   = 0;
hist_estim = 0;
hist_x = histCenters(1)-0.5*histBinWidth;
for ii = 1:length(histCenters)
    hist_vis =   [hist_vis   histCounts_vis(ii)*[1 1]];
    hist_estim = [hist_estim histCounts_estim(ii)*[1 1]];
    hist_x =     [hist_x histCenters(ii) + histBinWidth*[-0.5 0.5]];
end
hist_vis   = [hist_vis 0];
hist_estim = [hist_estim 0];
hist_x     = [hist_x histCenters(end) + histBinWidth*0.5];

visChosenEst = vel_ests_brute_vis(visChosen);


figure; hold on
fill(hist_x, hist_vis, [0 0 0])
fill(hist_x, -hist_estim, [1 0 0])
plot(bar_speed*[1 1], [-10 10], 'k--')

% mark velocity estimate for replicated visual trial
plot(visChosenEst, histCounts_vis(histCenters==visChosenEst)+1, 'kv')

set(gca, 'xlim', [histCenters(1)-1 histCenters(end)+1])

title(['velocity estimates (brute force)' 10 date_piece filesep EstimData  ', vis trial ' num2str(visChosen) ', m' num2str(movieChunk)])


%% more efficient version (search-based) to estimate velocity for each trial
% 
% %%
% % plot histogram of velocity estimates
% figure
% histCenters = 200:300;
% hist(mb_decode_out.vel_ests_eff, histCenters);
% title(['speed estimates (search-based)'])



