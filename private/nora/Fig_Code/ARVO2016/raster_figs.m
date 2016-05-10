clear
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 6 4])
GS = 1;
sta_scale = 4;
params_list

%% load movie
if GS
    i_chunk = 1;
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NSmovie = movie_chunk;
    interval(i_chunk) = size(movie_chunk,3);
    i_chunk = 2;
    
    % mask movie
    while size(NSmovie,3) < 1200
        load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
        interval(i_chunk) = size(movie_chunk,3);
        NSmovie = cat(3,NSmovie, movie_chunk);
        if size(NSmovie,3) > 1200
            NSmovie = NSmovie(:,:,1:1200);
        end
        i_chunk = i_chunk + 1;
    end
    NSmovie = (NSmovie - 64)/255;
    NSmovie = imresize(NSmovie, 1/sta_scale, 'box');
    NSmovie = permute(NSmovie, [2 1 3]);
    
end

interval_frame = cumsum(interval);


%%
e = 5; c = 5; SG = 2; 
for exp = e
    
    eval(exps{exp})
    n_reg = length(reg);
    mkdir(fig_save);
    n_subgroups = length(cell_idx);

    % load reg repeats
    for i_reg = 1:n_reg
        datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
        datarun{i_reg} = load_neurons(datarun{i_reg});
        reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 1);
        clear datarun
    end
    reg_data{i_reg}.testspikes = reg_data{i_reg}.testspikes([1:6 8:end], :);
    
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
    
    % load classification run
    datarun_class = load_data(class);
    datarun_class = load_params(datarun_class);

    glm_files = dir(['/Volumes/Lab/Users/Nora/GLMFits_masking/' piece '/NSEM_full_opt/*NSEM_optmodel.mat']);
    if strcmp(piece, '2016-04-21-1')
        glm_fit_idx = [6 3 5 1 2 4];
    end
    surround_offset = 0;
    %%
    % 1 = center, 2 = surround
    % 1 = reg, 2 = spot, 3 = gap
    for subgroup = SG
        for i_cell = c
            disp(cells_orig(i_cell))
            master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
            [PSTH{1}, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0, 'smoothing', 10);
                if strcmp(piece, '2016-04-21-1')
                    load(['/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1/NSEM_full_opt/' glm_files(glm_fit_idx(i_cell)).name]);
                else
                    load(['/Volumes/Lab/Users/Nora/GLMFits_masking/' piece '/NSEM_full_opt/' num2str(cells_orig(i_cell)) 'NSEM_optmodel.mat']);
                end
                
                sta_fit = opt_model.sta_fit.params;
                
                spatial{1} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.center_sd_x, 'sd_y',  sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
                spatial{2} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.surround_sd_scale*sta_fit.center_sd_x, 'sd_y',sta_fit.surround_sd_scale*sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
                
                for i = 1:2
                    spatial{i} = flip(spatial{i},1);
                    spatial{i} = flip(spatial{i},2);
                end
                if strcmp(cell_type{subgroup}, 'on')
                    time_course{1}= -opt_model.new_time;
                else
                    time_course{1}= opt_model.new_time;
                end

                time_course{2} = [time_course{1}(1:(end-surround_offset)); zeros(surround_offset,1)];
                
                for i = 1
                    s = sigmas(mask_conditions{subgroup}(i));
                    if strcmp(piece, '2016-04-21-1')
                        if s<8
                            load(['/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_allcells_sigma' num2str(s) '.mat'])
                        else
                            load(['/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat'])
                        end
                    else
                        if s==2
                            load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_allcells_sigma' num2str(s) '.mat'])
                        else
                            load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat']);
                        end
                    end
                    
                    % organize the masks
                    comp_mask = mod(mask+1,2);
                    mask = imresize(mask, 1/sta_scale, 'box');
                    mask = repmat(mask, 1, 1, 1200);
                    comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
                    comp_mask = repmat(comp_mask, 1, 1, 1200);
                    movie{1} = NSmovie;
                    movie{2} = NSmovie .* mask;
                    movie{3} = NSmovie.*comp_mask;
                    
                    % find the generator signals
                    for i_sta_part = 1:2
                        for i_movie_type = 1:3
                            spatial_GS = squeeze(convn(movie{i_movie_type}, spatial{i_sta_part}, 'valid'));
                            gen_signal{i_movie_type}{i_sta_part} = conv(spatial_GS, time_course{i_sta_part}, 'valid');
                        end
                    end
                    
                    PSTH{2} = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
                    PSTH{3} = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);

                    model_arch = 'y~b1/(b2+exp(b3*(x1+b4*x2)))';
                    init = [opt_model.model.Coefficients.Estimate; -sta_fit.surround_amp_scale];
                    
                    close all
                    
                    for i_movie_type = 1:3
                        % get the correlation and the generator signal for
                        % later
                        data_matrix = [gen_signal{i_movie_type}{1}(1:1166) gen_signal{i_movie_type}{2}(1:1166)];
                        responses = PSTH{i_movie_type}(30:end);
                        model = fitnlm(data_matrix, responses, model_arch, init);

                        %raster plotting
                        firing_rate = predict(model, data_matrix);
                        switch i_movie_type
                            case 1
                                IDP_plot_raster(reg_data{1},i_cell, firing_rate);
                            case 2
                                IDP_plot_raster(condition{mask_conditions{subgroup}(i)},i_cell, firing_rate);
                        end
                        
                        axis off
                    end
                    
                end
        end
        
    end

    
end

%% 
figure(1); ylim([1 58]); exportfig(gcf, ['/Users/Nora/Desktop/cell' num2str(i_cell) '_raster_full.eps'], 'Bounds', 'loose', 'Color', 'rgb')
figure(2); ylim([2 59]);  exportfig(gcf, ['/Users/Nora/Desktop/cell' num2str(i_cell) 'sigma' num2str(i) '_raster_mask.eps'], 'Bounds', 'loose', 'Color', 'rgb')

%% also plog rasters for background
% NB 2016-04-25

%%{
% choose cell and experiment 
exp_nm = '2013-08-19-6';
cell_savename = 'OFFPar_1101';
cell_no = 23;

% basic info
GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';%
GLMType.map_type = 'mapPRJ';
GLMType.fit_type = 'NSEM';
tstim = .00832750;
bins_per_frame = 10;

% load Alex's spikes
alex_fit_dir = '/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/rk1_MU_PS_noCP_timekernelCONEMODEL_stimnonlin_log_powerraise/standardparams/PS_netinhibitory_domainconstrain_COB/NSEM_mapPRJ/';
load([alex_fit_dir exp_nm '/' cell_savename '.mat'])
figure; plotrasters(fittedGLM.xvalperformance, fittedGLM, 'raster_length', 6)
exportfig(gcf, ['/Users/Nora/Desktop/AKH_NSEM_cell_raster_full.eps'], 'Bounds', 'loose', 'Color', 'rgb')
alex_fit_dir = '/Volumes/Lab/Users/akheitman/NSEM_Home/GLM_Output_Analysis/rk1_MU_PS_noCP_timekernelCONEMODEL_stimnonlin_log_powerraise/standardparams/PS_netinhibitory_domainconstrain_COB/WN_mapPRJ/';
load([alex_fit_dir exp_nm '/' cell_savename '.mat'])
figure; plotrasters(fittedGLM.xvalperformance, fittedGLM, 'raster_length', 6)     
exportfig(gcf, ['/Users/Nora/Desktop/SKH_WN_cell_raster_full.eps'], 'Bounds', 'loose', 'Color', 'rgb')
%}