datarun.names.rrs_neurons_path='/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.neurons';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.sta';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-80x60.xml';

slashes = strfind(datarun.names.rrs_neurons_path, '/');
dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
% num_frames = 30; % both have to be run with the name number of frames

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true,'load_all',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

cell_specification = [502,860,1024,1130,2076,2361,2618,2705,3172,3213,3559,4022,4071,4238,4774,4852,5496,6518,6533,6860,7279,7671];

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
    %temp_fit_params = fit_sta(temp_sta, 'fit_n_filters', true, 'initial_n_filters', 6, 'initial_scale_one',0.25,'initial_scale_two',-0.1,'initial_tau_one',3.15,'initial_tau_two',6, 'fit_surround_sd_scale', true, 'fit_surround', true);
    temp_fit_params = fit_sta(temp_sta, 'fit_n_filters', true, 'fit_surround_sd_scale', true, 'fit_surround', true);
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
save([dataset,'/large-cells'], 'information')
end


%% Look at the results
figure
suptitle({information{1}; [num2str(size(information{3},1)) ' large cells']})

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

