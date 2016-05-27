clear
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 12 4])
GS = 1;
sta_scale = 4;
colors = colormap(winter(6));
for i = 1:5
    Corr_ON{i} = [];
    Corr_OFF{i} = [];
end

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
params_list

for i_exp = [2 3 5]
    eval(exps{i_exp}) % 2 3 5
    n_reg = length(reg);
    mkdir(fig_save);
    n_subgroups = length(cell_idx);
    
    %% load and organize repeats
    
    % load reg repeats
    for i_reg = 1:n_reg
        datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
        datarun{i_reg} = load_neurons(datarun{i_reg});
        reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 0);
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
    % avg RF
    %datarun_class = load_sta(datarun_class);
     
    %%
    % plot the PSTH for each condition for each cell
    
    glm_files = dir(['/Volumes/Lab/Users/Nora/GLMFits/GLMFits_masking/' piece '/NSEM_full_opt/*NSEM_optmodel.mat']);
    if strcmp(piece, '2016-04-21-1')
        glm_fit_idx = [6 3 5 1 2 4];
    end
    surround_offset = 0;
    %%
    % 1 = center, 2 = surround
    % 1 = reg, 2 = spot, 3 = gap
    for subgroup = 1:n_subgroups
        for i_cell = cell_idx{subgroup}
            disp(cells_orig(i_cell))
            master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));

            if strcmp(piece, '2016-04-21-1')
                load(['/Volumes/Lab/Users/Nora/GLMFits/GLMFits_masking/2016-04-21-1/NSEM_full_opt/' glm_files(glm_fit_idx(i_cell)).name]);
            else
                load(['/Volumes/Lab/Users/Nora/GLMFits/GLMFits_masking/' piece '/NSEM_full_opt/' num2str(cells_orig(i_cell)) 'NSEM_optmodel.mat']);
            end
            
            sta_fit = opt_model.sta_fit.params;
            
            spatial{1} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.center_sd_x, 'sd_y',  sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
            spatial{2} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.surround_sd_scale*sta_fit.center_sd_x, 'sd_y',sta_fit.surround_sd_scale*sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
            for i = 1:2
                spatial{i} = flip(spatial{i},1);
                spatial{i} = flip(spatial{i},2);
            end
            spatial_final = spatial{1}-sta_fit.surround_amp_scale*spatial{2};
            clear spatial
            
            if strcmp(cell_type{subgroup}, 'on')
                time_course= -opt_model.new_time;
            else
                time_course= opt_model.new_time;
            end
            
            i = 4; s = sigmas(mask_conditions{subgroup}(i));
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
            
            
            PSTH = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 1);
            
            % organize the masks
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            movie = NSmovie .* mask;
            
            % find the generator signals
            spatial_GS = squeeze(convn(movie, spatial_final, 'valid'));
            gen_signal = conv(spatial_GS, time_course, 'valid');
            
            model_arch = 'y~b1/(b2+exp(b3*(x1)))';
            init = [opt_model.model.Coefficients.Estimate];
            
            data_matrix = gen_signal(1:1166);
            responses = PSTH(30:end);
            model = fitnlm(data_matrix, responses, model_arch, init);
            %model_PSTH = predict(model, data_matrix);
            %temp = corrcoef(responses, model_PSTH);
            %temp = temp(2,1);
            temp = model.Rsquared.Ordinary;
            thresh = 0;
            if strcmp(cell_type{subgroup}, 'on') && temp>thresh
                Corr_ON{i_exp} = [Corr_ON{i_exp} temp];
            elseif temp>thresh
                Corr_OFF{i_exp} = [Corr_OFF{i_exp} temp];
            end
            clear temp
  
        end
    end
    
end
%%
%{
default_colors = get(gca,'ColorOrder');
exp_means = [mean(Corr_ON) mean(Corr_OFF)];
exp_std = [std(Corr_ON)/sqrt(length(Corr_ON)) std(Corr_OFF)/sqrt(length(Corr_OFF))];
bar(1:2, exp_means, 'FaceColor', default_colors(1,:))
hold on; 
errorbar(exp_means, exp_std, '.k', 'LineWidth', 5)
set(gca, 'XTickLabels', {'ON', 'OFF'});
ylabel('Model Performance (R^2)')
%}
%%
figure; hold on
exp = [2 3 5];
for i=1:3
   plot((i+0.15)*ones(size(Corr_ON{exp(i)})),Corr_ON{exp(i)},'.', 'Color', default_colors(1,:), 'MarkerSize', 30)
   plot((i-0.15)*ones(size(Corr_OFF{exp(i)})),Corr_OFF{exp(i)},'.', 'Color', default_colors(2,:), 'MarkerSize', 30)
end
ylim([0.3 0.9])
set(gca,'XTick',[0.85 1.15 1.85 2.15 2.85 3.15], 'XTickLabel', {'ON', 'OFF'});
%legend({'On', 'Off'})
ylabel('Model Performance (R^2')