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

% loop over cells and fit
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
    
end



%datarun = compute_sta_fits(datarun, cell_specification, 'verbose', true);
cell_indices = get_cell_indices(datarun, cell_specification);
%output_matrix = make_Gaussian_two_d('center_point_x', datarun.matlab.sta_fits{cell_indices}.center_point_x, 'center_point_y', datarun.matlab.sta_fits{cell_indices}.center_point_y, 'rotation_angle', datarun.matlab.sta_fits{cell_indices}.center_rotation_angle, 'amp_scale', datarun.matlab.sta_fits{cell_indices}.surround_amp_scale, 'sd_x', datarun.matlab.sta_fits{cell_indices}.center_sd_x, 'sd_y',datarun.matlab.sta_fits{cell_indices}.center_sd_y, 'x_dim', datarun.matlab.sta_fits{cell_indices}.x_dim, 'y_dim', datarun.matlab.sta_fits{cell_indices}.y_dim);

