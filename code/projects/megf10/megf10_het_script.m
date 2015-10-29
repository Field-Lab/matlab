%% load data
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/Classification/
addpath /Volumes/lab/Experiments/Array/Shared/sravi/Tchaikowsky/sravi/matlab/DS' cell analysis'/
cd ~/Desktop/DS_code/
opt = struct('load_params', 1,'load_neurons', 1, 'load_ei', 1);

% mgf10 HET Animal
datamb{1} = load_data('/Volumes/lab/Experiments/Array/Analysis/2015-06-09-0/data003/data003', opt);
datamb{1}.names.stimulus_path = '/Volumes/lab/Experiments/Array/Analysis/2015-06-09-0/stimuli/s03.mat';
datamb{1} = load_stim_matlab(datamb{1});


%% classify DS vs non-DS cells (using moving bars)
  
i = 1; % which datarun to use for classification
[NumSpikesCell, MaxRate, StimComb] = get_spikescellstim_mb_gdf(datamb{i},datamb{i}.cell_ids, datamb{1}.triggers(end)+10, 1);
ds_struct = mbcellanalysis(NumSpikesCell, StimComb);

a = 1; b = 2; % which parameters to use for classification
Taylor1 = ds_struct.MAG{a,1}./((sum(ds_struct.RHO{a,1},2))');
Taylor2 = ds_struct.MAG{b,1}./((sum(ds_struct.RHO{b,1},2))');

figure(1); clf;
plot(log(Taylor1), log(Taylor2), 'o')
title('ndf0')
xlabel('TP 1')
ylabel('TP 2')
hold on

[x, y] = ginput;
plot(x, y);
IN = inpolygon(log(Taylor1), log(Taylor2), x, y);
[~, I] = find(IN == 1);
id_init = datamb{i}.cell_ids(I);

[C ia ib] = intersect(id_init, datamb{i}.cell_ids);
idx = ones(length(datamb{i}.cell_ids),1);
idx(ib) = 2; %initializing ds cells to cluster 2, everything else cluster 1

% Runs Gaussian mixture model to determine clusters based on chosen cells.
% close all;
% X = [];
% N = [];
% p = [];
% X(:,1) = log(Taylor1)';
% X(:,2) = log(Taylor2)';
% [idx N p] = clustering_analysis_plots(X, 0,1, 2, 0, 1, 0, 0, 0,0, idx);
% 
ds_id = [];
ds_id = datamb{i}.cell_ids(idx==2);
nonds_id = datamb{i}.cell_ids(idx==1);

%%
% get number of directions
num_conditions = length(datamb{1}.stimulus.trial_list);
% extract rates for low speed trials
low_speed = 2;
slow_indices = find(StimComb(:,2) == low_speed);

[direction_list, sorted_dir_indices] = sort(StimComb(slow_indices,3), 'ascend');
% sort indices into a ascending angles 
num_reps = datamb{1}.stimulus.repetitions;
stim_duration = 10.8;

tuning_struct = [];
analyzed_cell_num = 1;
for rgc = 1:length(ds_id);

    temp_index = get_cell_indices(datamb{1}, ds_id(rgc));
    
    figure(1); clf;

    for dr = 1:length(direction_list)
        temp_triggers = [];
        for rep = 1:num_reps
            temp_trigger = slow_indices(sorted_dir_indices(dr)) + rep + ((rep-1) * num_conditions);
            temp_triggers = [temp_triggers; temp_trigger];
        end
        temp_epochs = get_raster(datamb{1}.spikes{temp_index}, datamb{1}.triggers(temp_triggers), 'stop', stim_duration,'plot', false);   
        mb_rasters(dr).rasters = temp_epochs;
        mb_rasters(dr).direction = direction_list(dr);

        subplot(3,4,dr)
        plot_raster(temp_epochs, 0, stim_duration)
    end

    raster_cut_times = [];
    [raster_cut_times,y] = ginput;
    
    if isempty(raster_cut_times)
        continue
    end

    for dr = 1:length(direction_list)
        temp_raster = mb_rasters(dr).rasters;
        spike_counter_one = 0;
        spike_counter_two = 0;
        for trial = 1:size(temp_raster,1)
            temp_count = length(find(temp_raster{trial} < raster_cut_times(dr)));
            spike_counter_one = temp_count + spike_counter_one;
            temp_count = length(find(temp_raster{trial} > raster_cut_times(dr)));
            spike_counter_two = temp_count + spike_counter_two;       
        end
        first_half_counter(dr) = spike_counter_one;
        second_half_counter(dr) = spike_counter_two;
    end


    ds_cell(rgc).first_half_counter = first_half_counter;
    ds_cell(rgc).second_half_counter = second_half_counter;


    % plot tuning curves for each phase of response
    figure(1); clf;
    subplot(2,1,1)
    polar([direction_list' 360] * (pi/180), [second_half_counter second_half_counter(1)], 'r')
    hold on
    polar([direction_list' 360] * (pi/180), [first_half_counter first_half_counter(1)], 'k')
    hold off

   % plot tuning curve
    subplot(2,1,2)
    plot(direction_list, second_half_counter ./ max(second_half_counter), 'ro')
    hold on
    plot(direction_list, first_half_counter ./ max(first_half_counter), 'ko')
    hold off
    drawnow
    pause

    temp_width_one = circ_std(direction_list * (pi/180), first_half_counter', 30*pi/180) *180/pi;
    temp_width_two = circ_std(direction_list * (pi/180), second_half_counter', 30*pi/180) *180/pi;
    
    first_phase_width(analyzed_cell_num) = temp_width_one;
    second_phase_width(analyzed_cell_num) = temp_width_two;
    
    temp_dir = circ_mean(direction_list * (pi/180), first_half_counter') *180/pi;
    if temp_dir < 0
        temp_dir = 360 - temp_dir;
    end
    first_phase_dir(analyzed_cell_num) = temp_dir;
    
    temp_dir = circ_mean(direction_list * (pi/180), second_half_counter') *180/pi;
    if temp_dir < 0
        temp_dir = 360 - temp_dir;
    end
    second_phase_dir(analyzed_cell_num) = temp_dir;
    
    analyzed_cell_num = analyzed_cell_num + 1;
end


figure(1); clf;
subplot(1,2,1)
plot(first_phase_width, second_phase_width, 'ko')
hold on
plot([0:100],[0:100], 'k')
axis([0 100 0 100])
axis square
xlabel('first phase tuning width')
ylabel('second phase tuning width')
hold off

subplot(1,2,2)
plot(first_phase_dir, second_phase_dir, 'ko')
hold on
plot([0:360],[0:360], 'k')
axis([0 360 0 360])
axis square
xlabel('first phase direction')
ylabel('second phase direction')
hold off
