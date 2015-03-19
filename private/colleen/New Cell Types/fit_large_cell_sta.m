clear
file_name = '2008-08-27-6/data009-mg/data009';
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];

datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-10-1-0.48-11111.xml';
interpolate = false;
% cell_specification = [502,860,1024,1130,2076,2361,2618,2705,3022,3172,3213,3559,4022,4071,4238,4774,4852,5496,6518,6533,6860,7279,7671];
cell_type = {'on parasol'};
slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
to_replace = strfind(dataset, '/');
dataset(to_replace) = '-';
num_frames = 20; % both have to be run with the name number of frames

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames = 11:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

cell_type_index= zeros(1,size(cell_type,2));
for num_cell_types = 1:size(cell_type,2)
    for i = 1:size(datarun.cell_types,2)
        right_cell_type = strcmpi(datarun.cell_types{i}.name, cell_type{num_cell_types}); % case insensitive
        if right_cell_type == 1;
            cell_type_index(num_cell_types) = i;
            break
        end
        cell_type_index(num_cell_types) = 0;% couldn't find the right cell type
    end
    
end

cell_specification = datarun.cell_types{cell_type_index}.cell_ids;
% cell_specification = [6241];

% Find out indices for desired cell type


triggers=datarun.triggers; %onsets of the stimulus presentation

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


cell_indices = get_cell_indices(datarun, cell_specification);
num_rgcs = length(cell_indices);
parameters= zeros(num_rgcs,21);
variables = {'cell_specification', 'cell_indices', 'center_point_x', 'center_point_y', 'center_sd_x', 'center_sd_y', 'center_rotation_angle', 'color_weight_a', 'color_weight_b', 'color_weight_c', 'x_dim', 'y_dim', 'surround_sd_scale', 'surround_amp_scale', 'scale_one', 'scale_two', 'tau_one', 'tau_two','n_filters', 'frame_number', 'rmse'};
information{1} = dataset;
information{2} = variables;
for rgc = 1:num_rgcs
    
    fprintf('fitting the STA for cell %d... \n', datarun.cell_ids(cell_indices(rgc)))
    
    temp_sta = datarun.stas.stas{cell_indices(rgc)};

    % fit_surround_sd_scale is necessary for any fitting to occur
%     temp_fit_params = fit_sta(temp_sta, 'fit_n_filters', true, 'initial_n_filters', 10, 'initial_scale_one',0.15,'initial_scale_two',-0.2,'initial_tau_one',5,'initial_tau_two',5, 'fit_surround_sd_scale', true, 'fit_surround', true);



                    
[temp_fit_params, sta, sta_temp, sig_stixels] = fit_sta(temp_sta, 'fit_n_filters', false, 'fit_surround_sd_scale', false, 'fit_surround', false, 'initial_n_filters', 8, 'interpolate', false, 'frame_number', num_frames);
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
