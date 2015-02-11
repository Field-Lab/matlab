datarun.names.rrs_neurons_path='/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.neurons';
datarun.names.rrs_sta_path = '/Volumes/Analysis/2010-09-24-0/data001-nwpca/data001-nwpca.sta';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-80x60.xml';
% num_frames = 30; % both have to be run with the name number of frames

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true,'load_all',true);
datarun=load_data(datarun,opt);

triggers=datarun.triggers; %onsets of the stimulus presentation

[mov,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

[mvi] = load_movie(mdf_file, triggers);

cell_specification = [139];
datarun = compute_sta_fits(datarun, cell_specification, 'verbose', true);
cell_indices = get_cell_indices(datarun, cell_specification);
output_matrix = make_Gaussian_two_d('center_point_x', datarun.matlab.sta_fits{cell_indices}.center_point_x, 'center_point_y', datarun.matlab.sta_fits{cell_indices}.center_point_y, 'rotation_angle', datarun.matlab.sta_fits{cell_indices}.center_rotation_angle, 'amp_scale', datarun.matlab.sta_fits{cell_indices}.surround_amp_scale, 'sd_x', datarun.matlab.sta_fits{cell_indices}.center_sd_x, 'sd_y',datarun.matlab.sta_fits{cell_indices}.center_sd_y, 'x_dim', datarun.matlab.sta_fits{cell_indices}.x_dim, 'y_dim', datarun.matlab.sta_fits{cell_indices}.y_dim)

