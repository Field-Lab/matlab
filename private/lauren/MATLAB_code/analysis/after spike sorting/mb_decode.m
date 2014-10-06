% Script to calculate motion signals across different cell populations
% This script calculates the motion signals for a data run displays several
% nice guis for visualization.  Can set lots of options for different
% behavior. 
%
% adapted from: Marvin Thielk 2013
% mthielk@salk.edu
% adapted for electrical stimulation analysis by Lauren Jepson 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vary the options set in run_opt to turn different behaviors on and off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load -> reloads datavel_ests_brute
% remote -> if true uses data on the network, otherwise, local data
% data_run -> which data run to load and analyze
% raster -> creates a plot where you can vary the cell number and see a
%           raster for all the trials to demonstrate the reliability of
%           signals
% trial_raster -> creates a plot where you can vary the trial number and
%                 see all the spiking cells arranged by x position of their
%                 receptive field
% trial_raster_shift -> creates a plot where you can vary the trial number
%                       and manually shift the cells spike timing to see
%                       them line up.
% manual_speed_tuning -> creates a plot that demonstrates the algorithm and
%                        allows you see the gaussians line up as you vary
%                        the velocity.  This is for pairwise motion signal.
% velocity_lim -> max velocity that auto_speed_tuning automatically
%                 calculates
% auto_speed_tuning -> for a pair of neurons calculates the motion signal
%                      at different velocities and plots it.  Useful to
%                      show the pairwise shape of the speed tuning cuve
% tau -> tuning parameter that dictates the width of the gaussians used
% pop_speed_tuning -> plots the tuning curve summed over all possible pairs
%                     of neurons. Parallelized.
% tol -> determines the various tolerances for the integration and
%        optimizations used.
% savefig -> if true the script tries to save the figures that take longer
%            to produce in a figs/data_set folder.  The folder must already
%            exist or it will throw an error.
% trial_num -> Which trial to use for pop_speed_tuning.  
%                   0 < trial_num <= ~50
% trial_estimate -> creates a histogram of the fminunc determined velocity
%                   estimate of all the trials.  Parallelized.
% trial_estimate_start -> the velocity to use as the initial guess in
%                         fminunc to estimate the velocity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% parameters that shouldn't need to change

clear all
addpath(genpath('/snle/home/snl-e/matlab-standard/private/mthielk'))

% options related to algorithm
run_opt.tau = .01; % tuning parameter (in seconds)
run_opt.tol = 1e-3; % used in pop_motion_signal (and fminunc)
run_opt.trial_estimate_start = 200;
run_opt.trial_num = 20; % > 0
run_opt.vel = 0:1:500; % velocity values to calculate motion signal for when generating motion signal curves

run_opt.cell_type = 'ON parasol'; % on/off parasol, on/off midget

mb_decode_out.run_opt = run_opt;

%% general use parameters

bar_speed = 240; %in pixels/second
stim_type = 'electrical'; % visual or electrical



%%%% dataset info for 2011-06-24-5
if 0
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2011-06-24-5/';
    rf_dataset = 'data000'; %used only to determine x-position of RF centers

    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_rf = [1 32 229 396 531 737]; %from rf_dataset
    cell_ids_mb = [2 33 183 332 214 783]; %from data001-from-data000.neurons
    cell_ids_er = [1 33 183 332 214 784]; %from data006.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat

    % electrical stimulation parameters.
    movieChunk = 30; %30 or 37
    EstimData = 'data015';

    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data001-from-data000/data001-from-data000.neurons'];
    startPadding = 0.15;
    endPadding = 0.05;
    barColor = 'white';
end

%%%% dataset info for 2011-07-14-0
if 0
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2011-07-14-0/';
    rf_dataset = 'data000'; %used only to determine x-position of RF centers

    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_rf = [47  187 303 473 498 713 784 903]; %from rf_dataset
    cell_ids_mb = [213 200 348 467 543 663 783 901]; %from mbVisNeuronsFile
    cell_ids_er = [48  187 304 407 498 618 783 918]; %from data005.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat

    % electrical stimulation parameters.
    movieChunk = 40; %30/33 or 37/40
    EstimData = 'data013';

    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data001/data001.neurons'];
    startPadding = 0.15; %check
    endPadding = 0.05; %check
    barColor = 'white';
end

%%%% dataset info for 2013-05-28-3/data012 (gray background)
if 1
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/';
    rf_dataset = 'data004'; %used only to determine x-position of RF centers (closest in time to MB vis stim)

    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_rf = [167 137 350 439 722 932]; %from rf_dataset (closest in time to MB vis stim)
    cell_ids_mb = [76  138 377 437 783 917]; %from mbVisNeuronsFile
    cell_ids_er = [166 138 349 437 784 916]; %from data010.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat

    % electrical stimulation parameters.
    movieChunk = 41;
    EstimData = 'data012';

    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data003-gw/data003-gw.neurons'];
    startPadding = 0.15; %check
    endPadding = 0.05; %check
    barColor = 'white';
end

%%%% dataset info for 2013-05-28-3/data013 (no light)
if 0
    analysisPathBase = '/snle/lab/Experiments/Array/Analysis/2013-05-28-3/';
    rf_dataset = 'data004'; %used only to determine x-position of RF centers (closest in time to MB vis stim)

    % cell ids, ordered according to electrical stimulation pattern number
    cell_ids_rf = [167 137 350 439 722 932]; %from rf_dataset (closest in time to MB vis stim)
    cell_ids_mb = [76  138 377 437 783 917]; %from mbVisNeuronsFile
    cell_ids_er = [166 138 349 437 784 916]; %from data010.neurons (used in elecResp files) -- should match values in *MovingBarSpikes.mat

    % electrical stimulation parameters.
    movieChunk = 41;
    EstimData = 'data013';

    % visual stimulation parameters
    mbVisNeuronsFile = [analysisPathBase 'data003-gw/data003-gw.neurons'];
    startPadding = 0.15; %check
    endPadding = 0.05; %check
    barColor = 'white';
end


nCells = length(cell_ids_rf);

if strcmpi(stim_type, 'electrical')
    output_save_path = [analysisPathBase EstimData filesep 'm' num2str(movieChunk) '_mb_decode_out.mat'];
else
    output_save_path = [analysisPathBase EstimData filesep 'vis_mb_decode_out.mat'];
end

%% get length of moving bar -- in future this should be extracted from MovingBarVisSpikes structure

[~, intervals] = loadMovingBarResponses(mbVisNeuronsFile, cell_ids_mb, barColor, 'trialPadding', [startPadding endPadding]);

trial_length = mean(cellfun(@diff, intervals)); %in seconds

%% load spike times (in response to electrical OR visual stimulation) and
% format for pop_motion_signal

if strcmpi(stim_type, 'electrical')
    % load structure of electrical stimulation response times from disk
    % load([analysisPathBase EstimData filesep 'm' num2str(movieChunk) '_mbSpikes.mat']);
    load([analysisPathBase EstimData filesep 'm' num2str(movieChunk) '_mbSpikes.mat']);

    
    %sort moving bar spikes according to cell_ids_er (so that cell_x_pos indeces match up);
    for ii = 1:length(mbSpikes)
        cellInd = find([mbSpikes.cellIDVis] == cell_ids_er(ii));
        mbSpikesAll(ii) = mbSpikes(cellInd); %#ok<SAGROW>
    end
    
    clear mbSpikes cellInd
    
    %reorganize data structure
    spikeTimesAll = cell(length(mbSpikesAll(1).spikeTimes),1);
    for jj = 1:length(mbSpikesAll(1).spikeTimes)
        spikeTimesAll{jj} = cell(nCells,1);
        for ii = 1:nCells
            spikeTimesAll{jj}{ii} = reshape(mbSpikesAll(ii).spikeTimes{jj}, [], 1)/20000;
        end
    end
    clear mbSpikesAll
elseif strcmpi(stim_type, 'visual')
    spikeTimesAll = loadMovingBarResponses(mbVisNeuronsFile, cell_ids_mb, barColor, 'trialPadding', [startPadding endPadding]);
else
    error('unrecognized stimulus type (must be ''electrical'' or ''visual''')
end

% spikes = cell(nCells,1);
% for ii = 1:nCells
%     spikes{ii} = reshape(mbSpikesAll(ii).spikeTimes{run_opt.trial_num}, [], 1)/20000; %in seconds
% end

nTrials = length(spikeTimesAll);

mb_decode_out.bar_speed = bar_speed;
mb_decode_out.stim_type = stim_type;
mb_decode_out.spikeTimesAll = spikeTimesAll;

%% load RF center positions from vision files
% note:
%
% datarun.stas.fits is in stixel coordinates wrt (1,1) at CENTER of TOP left
% stixel, and is positive in downward and rightward direction (as shown in
% Vision interactive)
%
% datarun.vision.sta_fits is in stixel coordinates wrt (0,0) at LOWER left CORNER of
% bottom left stixel, and is positive in upward and rightward directions
% (as shown in Vision interactive)

datarun.names.rrs_params_path =  [analysisPathBase rf_dataset filesep rf_dataset '.params'];
datarun.names.rrs_neurons_path = [analysisPathBase rf_dataset filesep rf_dataset '.neurons'];
datarun.names.rrs_sta_path =     [analysisPathBase rf_dataset filesep rf_dataset '.sta'];
datarun.names.rrs_globals_path = [analysisPathBase rf_dataset filesep rf_dataset '.globals'];

opt=struct('verbose', true,'load_params', true,'load_neurons', true, 'load_sta', true);
datarun=load_data(datarun,opt);
%datarun = load_index(datarun);

% gets x positions of RF centers in stixels from left EDGE of vis stim
cellInds =  get_cell_indices(datarun,cell_ids_rf);
cell_x_pos = cellfun(@(x) x.mean(1), datarun.vision.sta_fits(cellInds));
clear cellInds

%convert from stixels to pixels
if datarun.stimulus.stixel_width ~= 8 % check because sometimes it is guessed
    warn(['stixel width appears to be ' num2str(datarun.stimulus.stixel_width) '--make sure this is correct'])
end
cell_x_pos = cell_x_pos*datarun.stimulus.stixel_width;

mb_decode_out.datarun = datarun;
mb_decode_out.cell_x_pos = cell_x_pos;

%% plot rasters to make sure what you loaded looks right

figure
[~, cellOrder] = sort(cell_x_pos);
for i = 1:nCells
    axes('position', [0.1 (nCells - i + 1)/(nCells+1) 0.8 1/(nCells+2)])
    hold on
    cellInd = cellOrder(i);
    for k = 1:nTrials
        for j = 1:length(spikeTimesAll{k}{cellInd})
            if strcmpi(stim_type, 'visual') && k == run_opt.trial_num;
                plot(spikeTimesAll{k}{cellInd}(j), k-0.5, '.', 'markerFaceColor', 'none', 'markerEdgeColor', [0.8 0 0])
            else
                plot([spikeTimesAll{k}{cellInd}(j) spikeTimesAll{k}{cellInd}(j)], [k-1 k], '-', 'LineWidth', 1, 'color', [0.8 0 0])

            end
        end
    end
    hold off
    
    set(gca, 'yLim', [0 nTrials], 'xlim', [0 trial_length])
    if i == nCells
        xlabel('time (s)')
    else
        set(gca, 'xtick', [])
    end
    ylabel(['cell' num2str(cell_ids_mb(cellInd))])
end

%% run Marvin's algorithm for a single trial at a single velocity

spikes = spikeTimesAll{run_opt.trial_num};

velocity = 200; %in stixels/s?

indices_both = 1:nCells; %because spikes and cell_x_pos already only contain cells of interest in specific order
trigger = 0; %spike times are already in reference to trial start time

tic
if strcmp('7.10.0.499 (R2010a)', version) %old version -- doesn't have integral.m built in, so use alternative version of martin's code
    sig_strength = pop_motion_signal_no_integral(velocity, spikes, indices_both, indices_both, cell_x_pos, trigger, trial_length, run_opt.tau) %#ok<NOPTS>
else
    sig_strength = pop_motion_signal(velocity, spikes, indices_both, indices_both, cell_x_pos, trigger, trial_length, run_opt.tau) %#ok<NOPTS>
end
toc

%% generate motion signal curve for a single trial

sig_strength_curve = zeros(size(run_opt.vel));

tic
if strcmp('7.10.0.499 (R2010a)', version) %old version -- doesn't have integral.m built in, so use alternative version of martin's code
    for ii = 1:length(run_opt.vel)
        sig_strength_curve(ii) = pop_motion_signal_no_integral(run_opt.vel(ii), spikes, indices_both, indices_both, cell_x_pos, trigger, trial_length, run_opt.tau);
    end
else
    for ii = 1:length(run_opt.vel)
        sig_strength_curve(ii) = pop_motion_signal(run_opt.vel(ii), spikes, indices_both, indices_both, cell_x_pos, trigger, trial_length, run_opt.tau);
    end
end
toc

vel_est = (run_opt.vel(sig_strength_curve==max(sig_strength_curve)));

figure; hold on
plot(bar_speed*[1 1], [0 max(sig_strength_curve)*1.2], 'k--')
plot(run_opt.vel, sig_strength_curve)
xlabel('speed tuning')
ylabel('signal strength')

clear sig_strength_curve


%% estimate velocity for single trial

options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 30, 'LargeScale', 'off');

if strcmp('7.10.0.499 (R2010a)', version) %old version -- doesn't have integral.m built in, so use alternative version of martin's code
    v_estimate = fminunc(@(v) -pop_motion_signal_no_integral(v, spikes, indices_both,...
        indices_both, cell_x_pos, trigger, trial_length, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
else
    v_estimate = fminunc(@(v) -pop_motion_signal(v, spikes, indices_both,...
        indices_both, cell_x_pos, trigger, trial_length, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
end
fprintf('for trial %d, the estimated speed was %d', run_opt.trial_num, v_estimate)

%% run brute-force algorithm (calculates motion signal for each velocity and takes max) to estimate velocity for each trial

vel_ests_brute = zeros(nTrials,1);
sig_strength_curve = cell(nTrials,1);

for jj = 1:nTrials
    jj
    sig_strength_curve{jj} = zeros(size(run_opt.vel));
    
    theseSpikes = spikeTimesAll{jj};
    
    tic
    if strcmp('7.10.0.499 (R2010a)', version) %old version -- doesn't have integral.m built in, so use alternative version of martin's code
        for ii = 1:length(run_opt.vel)
            sig_strength_curve{jj}(ii) = pop_motion_signal_no_integral(run_opt.vel(ii), theseSpikes, indices_both, indices_both, cell_x_pos, trigger, trial_length, run_opt.tau);
        end
    else
        for ii = 1:length(run_opt.vel)
            sig_strength_curve{jj}(ii) = pop_motion_signal(run_opt.vel(ii), theseSpikes, indices_both, indices_both, cell_x_pos, trigger, trial_length, run_opt.tau);
        end
    end
    toc
    
    vel_ests_brute(jj) = (run_opt.vel(sig_strength_curve{jj}==max(sig_strength_curve{jj})));
end

mb_decode_out.vel_ests_brute = vel_ests_brute;
mb_decode_out.sig_strength_curve = sig_strength_curve;


%% plot motion signal curves on top of each other

lineColors = lines(nTrials);
figure; hold on
plot(bar_speed*[1 1], [0 max(sig_strength_curve{jj})*1.5], 'k--')
for jj = 1:nTrials
    plot(run_opt.vel, sig_strength_curve{jj}, 'color', lineColors(jj,:))
    plot(vel_ests_brute(jj), max(sig_strength_curve{jj}), '*', 'markerEdgeColor', lineColors(jj,:));
end
xlabel('speed (pixels/second)')
ylabel('net motion signal')
title(['speed tuning curve of each trial (' stim_type ' stimulation)'])

% plot histogram
figure
histCenters = 200:300;
hist(vel_ests_brute, histCenters);
title(['speed estimates (brute force, ' stim_type ' stimulation)'])


%% run more efficient version (search-based) to estimate velocity for each trial

vel_ests_eff = zeros(nTrials,1);

options = optimset('Display', 'iter', 'TolFun', run_opt.tol , 'MaxFunEvals', 30, 'LargeScale', 'off');
for jj = 1:nTrials
    theseSpikes = spikeTimesAll{jj};
    
    tic
    if strcmp('7.10.0.499 (R2010a)', version) %old version -- doesn't have integral.m built in, so use alternative version of martin's code
        vel_ests_eff(jj) = fminunc(@(v) -pop_motion_signal_no_integral(v, theseSpikes, indices_both,...
            indices_both, cell_x_pos, trigger, trial_length, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
    else
        vel_ests_eff(jj) = fminunc(@(v) -pop_motion_signal(v, theseSpikes, indices_both,...
            indices_both, cell_x_pos, trigger, trial_length, run_opt.tau, run_opt.tol*.1), run_opt.trial_estimate_start, options);
    end
    fprintf('for trial %d, the estimated speed was %d', run_opt.trial_num, vel_ests_eff(jj))
    toc
end

mb_decode_out.vel_ests_eff = vel_ests_eff;

%%
% plot histogram
figure
histCenters = 200:300;
hist(vel_ests_eff, histCenters);
title(['speed estimates (search-based, ' stim_type ' stimulation)'])

%% save results to file

% get rid of java objects (can't be saved in this way)
mb_decode_out.datarun = rmfield(mb_decode_out.datarun, 'globals');
mb_decode_out.datarun.stas = rmfield(mb_decode_out.datarun.stas, 'java_sta');

save(output_save_path, 'mb_decode_out')



