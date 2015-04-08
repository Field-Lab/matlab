%% --------------- Produce STAs from focal white noise -------------------
%% Function: Generate STAs for cells that were stimulated with the Voronoi stimulus

%% How to use: 
%     1) List data file name where neurons and params are.
%     2) Enter the Voronoi stimulus stixel size
%     3) Enter the cells (from Vision) in the same order as they were
%     entered to generate the mask in Voronoi_stimulus
%     4) Enter movie file

%% Potential problems:
%     1) Could have a problem parsing movie file correctly for multiple
%     cells (limited testing
%     2) Binning problems with matlab STA code
    

%% Inputs
% file_name : start with piece number and continue until the .neurons etc files (eg '2006-06-06-2/data003/data003')
% stixel_size : Voronoi stixel size
% cell_specification : cell IDs from Vision
% mdf_file : white noise xml 



%% Results
% STA : Computed in matlab not java, binning problems
% timecourse: From STA sig stixels
% significant stixels : standard parameters

%% Author 
% Colleen Rhoades (rhoades@stanford.edu)
% April 7, 2015
clear

%% -----------------INPUTS---------------------------------
file_name = '2006-06-06-2/data003/data003';
stixel_size = 16;
cell_specification = [1487,5462]; % put cells in same order as Voronoi_stimulus where the mask was generated
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-16-4-0.48-33333.xml';

%% --------------------------- Load Data -----------------------------

datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];

slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataset(to_replace) = '-';
num_frames = 30; % both have to be run with the name number of frames

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits', 0, 'load_sta', 0, 'load_sta_params', 0, 'load_all',true);
opt.load_sta_params.save_rf = 0;
opt.load_sta_params.frames = 1:30; % have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);
%% --------------------------- Get Spikes and Movie -----------------------------
% Find out indices for desired cell type
cellID = get_cell_indices(datarun, cell_specification);
for cell = 1:length(cellID)
    spikes=datarun.spikes{cellID(cell)};
    spikes = spikes*1000;
    triggers=datarun.triggers; %onsets of the stimulus presentation
    
    [mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
        triggers, 1,2);
    
    [mvi] = load_movie(mdf_file, triggers);
    

    
%% --------------------------- Compute STA of one cell -----------------------------

    [sta, timecourse, sig_stixels] = compute_only_sta(datarun, mdf_file, num_frames, spikes, 1, cell, length(cellID));
end
%% --------------------------- Fitting code below, STOP HERE -----------------------------

return
%%%%%%% ----------STOP HERE FOR JUST STA WITHOUT FIT--------------------
cell_indices = get_cell_indices(datarun, cell_specification);
num_rgcs = length(cell_indices);
parameters= zeros(num_rgcs,21);
variables = {'cell_specification', 'cell_indices', 'center_point_x', 'center_point_y', 'center_sd_x', 'center_sd_y', 'center_rotation_angle', 'color_weight_a', 'color_weight_b', 'color_weight_c', 'x_dim', 'y_dim', 'surround_sd_scale', 'surround_amp_scale', 'scale_one', 'scale_two', 'tau_one', 'tau_two','n_filters', 'frame_number', 'rmse'};
information{1} = dataset;
information{2} = variables;
for rgc = 1:num_rgcs
    
    fprintf('fitting the STA for cell %d... \n', datarun.cell_ids(cell_indices(rgc)))
    
    temp_sta = sta;
    
    % fit_surround_sd_scale is necessary for any fitting to occur
    %     temp_fit_params = fit_sta(temp_sta, 'fit_n_filters', true, 'initial_n_filters', 10, 'initial_scale_one',0.15,'initial_scale_two',-0.2,'initial_tau_one',5,'initial_tau_two',5, 'fit_surround_sd_scale', true, 'fit_surround', true);
    
    
    
    figure
    [temp_fit_params, sta, sta_temp, sig_stixels] = fit_sta(temp_sta, 'fit_n_filters', true, 'fit_surround_sd_scale', false, 'fit_surround', false, 'initial_n_filters', 8, 'interpolate', false, 'frame_number', num_frames, 'num_colors',1);
    % Plot result
    %  plot_sta_fit(sta_temp, temp_fit_params.fit_params, temp_fit_params.fixed_params, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, sig_stixels, 'off');
    
    
    % Add raw data to plot
    % hold on
    %         real_stix = significant_stixels(temp_sta);
    %       biggestBlob = ExtractNLargestBlobs(full(real_stix), 1);
    %     real_stix = biggestBlob;
    %     tc = time_course_from_sta(temp_sta, real_stix);
    %     norm_factor = max(abs(reshape(tc, 1, [])));
    %     tc = tc ./ norm_factor;
    %     hold on
    %     subplot(2,1,2)
    %         plot(tc(:,1), 'r', 'linewidth', 2)
    %         plot(tc(:,2), 'g', 'linewidth', 2)
    %         plot(tc(:,3), 'b', 'linewidth', 2)
    
    suptitle({dataset; sprintf('Fit for cell %d', datarun.cell_ids(cell_indices(rgc)))})
    
    
    if isempty(temp_fit_params)
        temp_id = datarun.cell_ids(cell_indices(rgc));
        warn_message = ['cell ',num2str(temp_id), ' has no sig stixels and no fit'];
        warning(warn_message)
    end
    
    datarun.matlab.sta_fits{cell_indices(rgc)} = temp_fit_params;
    
    
    
    
    parameters(rgc,1) = cell_specification(rgc);
    parameters(rgc,2) = cell_indices(rgc);
    parameters(rgc,3) = datarun.matlab.sta_fits{cell_indices(rgc)}.center_point_x;
    parameters(rgc,4) = datarun.matlab.sta_fits{cell_indices(rgc)}.center_point_y;
    parameters(rgc,5) = datarun.matlab.sta_fits{cell_indices(rgc)}.center_sd_x;
    parameters(rgc,6) = datarun.matlab.sta_fits{cell_indices(rgc)}.center_sd_y;
    parameters(rgc,7) = datarun.matlab.sta_fits{cell_indices(rgc)}.center_rotation_angle;
    parameters(rgc,8) = datarun.matlab.sta_fits{cell_indices(rgc)}.color_weight_a;
    parameters(rgc,9) = datarun.matlab.sta_fits{cell_indices(rgc)}.color_weight_b;
    parameters(rgc,10) = datarun.matlab.sta_fits{cell_indices(rgc)}.color_weight_c;
    parameters(rgc,11) = datarun.matlab.sta_fits{cell_indices(rgc)}.x_dim;
    parameters(rgc,12) = datarun.matlab.sta_fits{cell_indices(rgc)}.y_dim;
    parameters(rgc,13) = datarun.matlab.sta_fits{cell_indices(rgc)}.surround_sd_scale;
    parameters(rgc,14) = datarun.matlab.sta_fits{cell_indices(rgc)}.surround_amp_scale;
    parameters(rgc,15) = datarun.matlab.sta_fits{cell_indices(rgc)}.scale_one;
    parameters(rgc,16) = datarun.matlab.sta_fits{cell_indices(rgc)}.scale_two;
    parameters(rgc,17) = datarun.matlab.sta_fits{cell_indices(rgc)}.tau_one;
    parameters(rgc,18) = datarun.matlab.sta_fits{cell_indices(rgc)}.tau_two;
    parameters(rgc,19) = datarun.matlab.sta_fits{cell_indices(rgc)}.n_filters;
    parameters(rgc,20) = datarun.matlab.sta_fits{cell_indices(rgc)}.frame_number;
    parameters(rgc,21) = datarun.matlab.sta_fits{cell_indices(rgc)}.rmse;
    
    
    
    
    information{3} = parameters;
    information{4} = datarun.matlab.sta_fits;
    % save([dataset,'-', cell_type{1}], 'information')
    
end

figure
for q = 1:num_rgcs
    temp_sta = datarun.stas.stas{cell_indices(q)};
    temp_fit_params = datarun.matlab.sta_fits{cell_indices(q)};
    fit_tc{q} =  plot_fit_timecourses(temp_sta, temp_fit_params.fit_indices, temp_fit_params.fixed_indices, temp_fit_params.fit_params, temp_fit_params.fixed_params, sig_stixels, 1); % plot_raw = 1
    hold on
end
title({['Fits of ' num2str(length(cell_specification)), ' ', cell_type{1}, ' Cells']; dataset})
xlabel('Time')
ylabel('Amplitude')



%% Look at the results
figure
suptitle({information{1}; [num2str(size(information{3},1)), cell_type{1}, ' cells']})

%scale 1
scale_one = information{3}(:,15);
subplot(3,2,1)
hist(scale_one);
title('Scale One');

scale_two = information{3}(:,16);
subplot(3,2,2) ; hist(scale_two);
title('Scale Two');

tau_one = information{3}(:,17);
subplot(3,2,3) ; hist(tau_one);
title('Tau One');

tau_two = information{3}(:,18);
subplot(3,2,4) ; hist(tau_two);
title('Tau Two');

n = information{3}(:,19);
subplot(3,2,5) ; hist(n);
title('N Filters');

area = information{3}(:,5).*information{3}(:,6)*pi;
subplot(3,2,6) ; hist(area);
title('RF Size');



%% Plot mosaic
figure
for q = 1:num_rgcs
    temp_sta = datarun.stas.stas{cell_indices(q)};
    temp_fit_params = datarun.matlab.sta_fits{cell_indices(q)};
    fit_indices = temp_fit_params.fit_indices;
    fixed_indices = temp_fit_params.fixed_indices;
    fit_params = temp_fit_params.fit_params;
    fixed_params = temp_fit_params.fixed_params;
    
    all_params(fit_indices) = fit_params;
    all_params(fixed_indices) = fixed_params;
    
    % get sta fit
    sta_fit = sta_fit_function(all_params);
    % spatial fit
    
    
    hold on
    %     temp_rf = rf_from_sta(sta);
    %     imagesc(norm_image(temp_rf))
    %     hold on
    h(q) = plot_spatial_sd(all_params);
    set(h(q), 'DisplayName', 'Fitting Code');
    drawnow
end
axis equal

set(gca, 'xlim', [0, datarun.stimulus.field_width]);
set(gca, 'ylim', [0, datarun.stimulus.field_height]);
set(gca,'YDir','reverse');
%% compare mosaic to that from vision
% datarun=load_data('2010-09-24-0/data001-nwpca-cr');
% datarun=load_sta(datarun);
% datarun=load_params(datarun);
% datarun=load_neurons(datarun);
% datarun=set_polarities(datarun);
hold on
h(rgc+1) = plot_rf_summaries(datarun,{cell_type_index}, 'plot_fits',1,'label',0, 'foa', -1, 'clear', 0); %looking at the data type 8 is on large
set(h(rgc+1), 'DisplayName', 'Vision')
legend(h(end-1:end));
title({['STA Fitting Code:  ', cell_type{1}] ; 'Clean Spike Sorting'; dataset})



area  = parameters(:,5) .* parameters(:,6) .* pi;
diameters = 2*[parameters(:,5) , parameters(:,6)];
column1_larger = find(diameters(:,1)>diameters(:,2));
column2_larger = find(diameters(:,1)<=diameters(:,2));
max_diameters = [diameters(column1_larger,1); diameters(column2_larger,2)];
disp('Mean RF diameter')
mean(max_diameters)
%% Calculate zero crossing
zero_crossing = zeros(size(fit_tc,2),1);
for i = 1:size(fit_tc,2)
    cell_ = fit_tc{i}(:,2)*16.65; % units of ms
    
    [x0,y0] = intersections(1:num_frames,cell_,1:num_frames, zeros(num_frames,1)); % find zero crossing
    zero_crossing(i) = x0(1);
end
information{5} = zero_crossing;
save([dataset,'-', cell_type{1}], 'information')
