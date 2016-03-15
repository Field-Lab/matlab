% add parasol spikes
% Simulation to determine how much contamination of parasol spikes is
% needed for a blue/green cell to be affected

tic
% Information about run to be analysized
datarun.names.rrs_neurons_path='/Volumes/Analysis/2007-09-18-3/data000-mg/data000-cr2/data000-cr2.neurons';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-09-18-3/data000-mg/data000-cr2/data000-cr2.sta';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';

% First cell is parasol, second cell is blue/green
cell_specification = [1651,1652];
cell_type = 'SBC';

% Parce the name of the datarun for saving purposes
slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataset(to_replace) = '-';

% Load the data
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true,'load_all',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

%Load the movie
[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);
[mvi] = load_movie(mdf_file, triggers);


% check that stas have been loaded into datarun
if ~isfield(datarun.stas, 'stas')
    error('STAs are not contained in datarun, see LOAD_STA.M')
end

% initialize_output
if ~isfield(datarun, 'matlab')
    datarun = setfield(datarun, 'matlab', []);
elseif ~isfield(datarun.matlab, 'sta_fits')
    temp_cell = cell(length(datarun.cell_ids), 1);
    datarun.matlab.sta_fits = temp_cell;
end

% initialization for saving information
cell_indices = get_cell_indices(datarun, cell_specification);
num_rgcs = length(cell_indices);
parameters= zeros(num_rgcs,21);
variables = {'cell_specification', 'cell_indices', 'center_point_x', 'center_point_y', 'center_sd_x', 'center_sd_y', 'center_rotation_angle', 'color_weight_a', 'color_weight_b', 'color_weight_c', 'x_dim', 'y_dim', 'surround_sd_scale', 'surround_amp_scale', 'scale_one', 'scale_two', 'tau_one', 'tau_two','n_filters', 'frame_number', 'rmse'};
information{1} = dataset;
information{2} = variables;
sta = cell(2,1);
sig_stixels_sum = cell(2,1);
% Set this parameter based on timecourses
num_frames = 20; % both have to be run with the name number of frames
for i = 1:num_rgcs
    figure
    hb= tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]);
    if i == 1
        label = [cell_type, ' Spikes'];
    else
        label = 'Blue/Green Spikes';
    end
    
    
    spikes=datarun.spikes{cell_indices(i)};
    % Cut down parasol spikes because simulation takes too long
%     if length(spikes) > 20000
%         spikes = spikes(1:20000);
%     end
    % spikes are in s. convert to ms
    spikes=round(spikes*1000);
    % compute the STA in matlab (not using java information)
    [sta{i}, timecourse, sig_stixels] = compute_only_sta(datarun, mdf_file, num_frames, spikes, 0);
    sig_stixels_sum{i} = sig_stixels;
    for a = 1:3
        axes(hb(a)); 
        frame = 11;
        imagesc(squeeze(sta{i}(:,:,:, frame+a)));
        title({['STA Frame: ' num2str(frame+a)] ; label})
        axis image
        axis off
    end
  
    
        
    % Plot the time course from the significant stixels
    h = plot_time_course_(timecourse, 'colors', ['rgb']', 'foa', 0)
    title({'Time Course'; label});
end
% Find the peak frame of the STA
[junk,start_index] = max(sum(reshape(sta{2}.^2,[],size(sta{2},4)),1));

% Percent contamination
perc = [0,.05,.10,.20,.30,.40,.50,.75];
sta_perc = cell(size(perc,2),1);

% Get the spikes from both the b/g cell and the parasol
spikes1 = datarun.spikes{cell_indices(1)};
spikes2 = datarun.spikes{cell_indices(2)};
total_bg = length(spikes2);

 
    
for percentage = 1:length(perc)
%     percentage = 8;
    total_bg = 4250;
     figure
    ha = tight_subplot(2,3,[.01 .03],[.1 .01],[.01 .01]);
%     total_bg = 1500;
    % compute how many parasol spikes to add
    total_spikes = round(total_bg/(1-perc(percentage)));
    
    parasol_to_include = total_spikes - total_bg;
    % randomly select parasol spikes to add avoiding the first 100
    % spikes because the piece could be unstable
    p = randperm(size(spikes1,1)-100, parasol_to_include)+100; % +100 to account for spikes at the beginning of the train
    spikes=datarun.spikes{cell_indices(1)}(p);
    
    %Combine b/g and comtamination spikes into one variables
    spikes = [spikes; spikes2(100:4349)];
    
    spikes=round(spikes*1000);
    % Compute new sta with the b/g spikes and some parasol spikes
    [combined_sta, timecourse, sig_stixels] = compute_only_sta(datarun, mdf_file, num_frames, spikes, 0); % 0 indicates no plotting
    % Comput significant stixels into order to determine the timecourse
    %[sig_stixels] = significant_stixels(combined_sta);
    %[timecourse, params] = time_course_from_sta(combined_sta, sig_stixels);
    % Plot the important frames of the STA
  for a = 1:6
    axes(ha(a)); 
    frame = 11;
    if a < 4
        imagesc(squeeze(sta{2}(:,:,:, frame+a)));
        title({['STA Frame: ' num2str(frame+a)] ; 'Blue/Green Spikes'})
        axis image
        axis off
    else
        frame = 8;
        imagesc(squeeze(combined_sta(:,:,:, frame+a)));
        title({['STA Frame: ' num2str(frame+a)]; sprintf(['%d Percent ' cell_type ' Spikes'], perc(percentage)*100)});
           axis image
        axis off

    end
  end
  
    % Plot the time course from the significant stixels
    h = plot_time_course_(timecourse, 'colors', ['rgb']', 'foa', 0)
    title({'Time Course'; sprintf(['%d Percent ' cell_type ' Spikes'], perc(percentage)*100)});
   
    % Save resulting STA for quantification
    sta_perc{percentage} = combined_sta;
end
toc


% Per sig_stixel differences
% non blue/green cell type
diff_sig = zeros(sum(sum(sig_stixels_sum{1})), 3, num_frames, length(perc));
diff_bg = zeros(sum(sum(sig_stixels_sum{2})), 3, num_frames, length(perc));

for b = 1:length(perc)
    full_indexing = repmat(full(sig_stixels_sum{1}), [1,1,3,num_frames]);
    sig_para =  sta{1}(full_indexing);
        sig_para = reshape(sig_para, [sum(sum(sig_stixels_sum{1})), 3, num_frames]);

    sig_contam = sta_perc{b}(full_indexing);
    sig_contam = reshape(sig_contam, [sum(sum(sig_stixels_sum{1})), 3, num_frames]);
    diff_sig(:,:,:,b) = sig_para - sig_contam;
end
mean_para_change = squeeze(mean(diff_sig,1));
% select_frame = 12;
% select_color = 2;
% if select_color == 2
%     color_str = 'Green';
% elseif select_color == 3
%     color_str = 'Blue';
% else
%     color_str = 'Red';
% end
% h = figure; 
% plot(perc*100, squeeze(mean_para_change(select_color, select_frame,:)))
% set(h, 'Position', [680 460 650 650])
% %  set(h, {'color'}, num2cell(jet(sum(sum(sig_stixels_sum{2}))), 2));
% xlabel('Percent of contaminated pixels')
% ylabel(['Mean(',  cell_type , ' sign.stixels - contaminated sign. stixels)'])
% title({['Contamination with ', cell_type, ' Cells']; dataset; [num2str(cell_specification(1)), ' and ' num2str(cell_specification(2))]; ['Frame ',num2str(select_frame), ' ', color_str]})


figure; plot(perc*100, squeeze(diff_sig(:,2,12,:)))
xlabel('Percent of contaminated pixels')
ylabel('Contaminated sign. stixels - Parasol sign.stixels')
% legend(legendCell)

% Blue/green
for c = 1:length(perc)
    full_indexing = repmat(full(sig_stixels_sum{2}), [1,1,3,num_frames]);
    sig_bg =  sta{2}(full_indexing);
        sig_bg = reshape(sig_bg, [sum(sum(sig_stixels_sum{2})), 3, num_frames]);

    sig_contam = sta_perc{c}(full_indexing);
    sig_contam = reshape(sig_contam, [sum(sum(sig_stixels_sum{2})), 3, num_frames]);
    diff_bg(:,:,:,c) = sig_contam - sig_bg;
    
end
mean_bg_change = squeeze(mean(diff_bg,1));


select_frame = 12;
select_color = 1;
if select_color == 2
    color_str = 'Green';
elseif select_color == 3
    color_str = 'Blue';
else
    color_str = 'Red';
end

g = figure; 
plot(perc*100, squeeze(mean_bg_change(1, select_frame,:)), 'r')
hold on
plot(perc*100, squeeze(mean_bg_change(1, select_frame+1,:)), 'r:')
plot(perc*100, squeeze(mean_bg_change(1, select_frame+2,:)), 'r--')

plot(perc*100, squeeze(mean_bg_change(2, select_frame,:)), 'g')
plot(perc*100, squeeze(mean_bg_change(3, select_frame,:)), 'b')


plot(perc*100, squeeze(mean_bg_change(2, select_frame+1,:)), 'g:')
plot(perc*100, squeeze(mean_bg_change(3, select_frame+1,:)), 'b:')


plot(perc*100, squeeze(mean_bg_change(2, select_frame+2,:)), 'g--')
plot(perc*100, squeeze(mean_bg_change(3, select_frame+2,:)), 'b--')
set(g, 'Position', [680 460 650 650])
xlabel('Percent of contaminated pixels')
ylabel('Mean(Contaminated sign. stixels - Blue/Green sign.stixels)')
title({['Contamination with ', cell_type, ' Cells']; dataset; [num2str(cell_specification(1)), ' and ' num2str(cell_specification(2))]})
legend(['Frame ' num2str(select_frame)], ['Frame ' num2str(select_frame+1)], ['Frame ' num2str(select_frame+2)])