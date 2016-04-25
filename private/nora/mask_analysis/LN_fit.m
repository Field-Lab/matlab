clear 
params_201602176
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

n_subgroups = length(cell_idx);
GS = 1;
sta_scale = 4;
colors = colormap(winter(6));


%% load and organize repeats 

% load reg repeats
for i_reg = 1:n_reg
    datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
    datarun{i_reg} = load_neurons(datarun{i_reg});
    reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 1);
    clear datarun
end

% load masking repeats
datarun = load_data([ Analysis_Path masking '/' masking]);
datarun = load_neurons(datarun);
mask_data = interleaved_data_prep(datarun, 1100, n_masks*2*30, 'cell_spec', cells_masking, 'visual_check', 0);
clear datarun
% separate each condition
idx = 1:30;
for count = 1:(n_masks*2)
    condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
end

%%
% load classification run 
datarun_class = load_data(class);
datarun_class = load_params(datarun_class);
datarun_class = load_sta(datarun_class);


%% load movie
i_chunk = 1;
load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
NSmovie = movie_chunk;
i_chunk = 2;

% mask movie
while size(NSmovie,3) < 1200
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NSmovie = cat(3,NSmovie, movie_chunk);
    if size(NSmovie,3) > 1200
        NSmovie = NSmovie(:,:,1:1200);
    end
    i_chunk = i_chunk + 1;
end
NSmovie = NSmovie - 64;
NSmovie = imresize(NSmovie, 1/sta_scale, 'box');
NSmovie = permute(NSmovie, [2 1 3]);


%%
% plot the PSTH for each condition for each cell

for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        disp(cells_orig(i_cell))
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        [reg, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0, 'smoothing', 10);
        
        % get the spatial filted from the STA fit
        [params,~,~] = fit_sta(datarun_class.stas.stas{master_idx}, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true); 
        ideal_center = make_Gaussian_two_d('center_point_x', params.center_point_x, 'center_point_y', params.center_point_y, 'sd_x', params.center_sd_x, 'sd_y', params.center_sd_y, 'rotation_angle', params.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
        ideal_surround = make_Gaussian_two_d('center_point_x', params.center_point_x, 'center_point_y', params.center_point_y, 'sd_x', params.surround_sd_scale*params.center_sd_x, 'sd_y', params.surround_sd_scale*params.center_sd_y, 'rotation_angle', params.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
        ideal_center = flip(ideal_center,1); 
        ideal_center = flip(ideal_center,2); 
        ideal_surround = flip(ideal_surround,1); 
        ideal_surround = flip(ideal_surround,2);


        
        for i = 1:length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            if s==2
                load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_allcells_sigma' num2str(s) '.mat'])
            else
                load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat']);
            end
           
            % getting the right mask movie
            comp_mask = mod(mask+1,2);  
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
            comp_mask = repmat(comp_mask, 1, 1, 1200);
            mask_movie = NSmovie .* mask;
            
            % filtering
            mask_gen_signal_center = squeeze(convn(mask_movie, ideal_center, 'valid')); 
            mask_gen_signal_surround = squeeze(convn(mask_movie, ideal_surround, 'valid')); 
            reg_gen_signal_center = squeeze(convn(NSmovie, ideal_center, 'valid')); 
            reg_gen_signal_surround = squeeze(convn(NSmovie, ideal_surround, 'valid')); 
            comp_gen_signal_center = squeeze(convn(NSmovie.*comp_mask, ideal_center, 'valid')); 
            comp_gen_signal_surround = squeeze(convn(NSmovie.*comp_mask, ideal_surround, 'valid')); 
            
            % responses
            mask_PSTH = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            comp = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
                        
            
            responses = reg(30:end);
            b_center_init = [params.scale_one params.tau_one params.n_one_filters params.scale_two params.tau_two params.n_two_filters];
            b_surround_init = [params.scale_one params.tau_one params.n_one_filters params.scale_two params.tau_two params.n_two_filters];
            b_NL_init = [1 1];
            [fit_params, fval] = fminsearch(@(b) model_func(b, reg_gen_signal_center, reg_gen_signal_surround, responses), [b_center_init, b_surround_init, b_NL_init]);


        end
    end
end 
close all

%%

for i_cell = 1:6
    disp(cells_orig(i_cell))
    master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
    [params{i_cell},filter{i_cell},~] = fit_sta(datarun_class.stas.stas{master_idx}, 'fit_surround', true, 'fit_surround_sd_scale', true, 'fit_surround_amp_scale', true);
end