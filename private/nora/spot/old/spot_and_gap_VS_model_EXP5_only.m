clear
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 12 4])
GS = 1;
sta_scale = 4;
colors = colormap(winter(6));
mean_cell_corr = cell(10,1);
mean_cell_surr = cell(10,1);
model_mean_cell_corr = cell(10,1);
model_mean_cell_surr = cell(10,1);

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

i_exp = [5];
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
    reg_data_first.testspikes = reg_data{i_reg}.testspikes(1:15,:);
    reg_data_second.testspikes = reg_data{i_reg}.testspikes(16:30,:);
    clear datarun
end

% load masking repeats
datarun = load_data([ Analysis_Path masking '/' masking]);
datarun = load_neurons(datarun);
mask_data = interleaved_data_prep(datarun, 1100, n_masks*2*30, 'cell_spec', cells_masking, 'visual_check', 0);
monitor_refresh = 100/median(diff(datarun.triggers));
clear datarun
% separate each condition
idx = 1:2:29;
for count = 1:(n_masks*2)
    condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
end

%% Separate each condition into halves to calculate the noise floors 
idx = 2:2:30;
for count = 1:(n_masks*2)
    second_half_condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
end

%% load classification run
datarun_class = load_data(class);
datarun_class = load_params(datarun_class);
% avg RF
%datarun_class = load_sta(datarun_class);

%% Load up the cells in question
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
        disp(i_cell)
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        
        % Use if comparing to full field rather than to spot size 10 
        %PSTH{1} = IDP_plot_PSTH(reg_data_first,i_cell, 'color', 0, 'smoothing', 1);
        %PSTH_compare = IDP_plot_PSTH(reg_data_second,i_cell, 'color', 0, 'smoothing', 1);
        
        if strcmp(piece, '2016-04-21-1')
            load(['/Volumes/Lab/Users/Nora/GLMFits/GLMFits_masking/2016-04-21-1/NSEM_full_opt/' glm_files(glm_fit_idx(i_cell)).name]);
        else
            load(['/Volumes/Lab/Users/Nora/GLMFits/GLMFits_masking/' piece '/NSEM_full_opt/' num2str(cells_orig(i_cell)) 'NSEM_optmodel.mat']);
        end
        
        % Create the center and surround filters from the STA Fit
        sta_fit = opt_model.sta_fit.params;
        spatial{1} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.center_sd_x, 'sd_y',  sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
        spatial{2} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.surround_sd_scale*sta_fit.center_sd_x, 'sd_y',sta_fit.surround_sd_scale*sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
        
        % flip the spatial part in preparation for the convolution
        for i = 1:2
            spatial{i} = flip(spatial{i},1);
            spatial{i} = flip(spatial{i},2);
        end
        
        % time course from GLM / NL-time iterative fit
        % Fix the timecourse for polarity of cell
        if strcmp(cell_type{subgroup}, 'on')
            time_course{1}= -opt_model.new_time;
        else
            time_course{1}= opt_model.new_time;
        end
        
        % time course from vision
        %time_course{1}=flip(mean([datarun_class.vision.timecourses(master_idx).r ...
        %    datarun_class.vision.timecourses(master_idx).g ...
        %    datarun_class.vision.timecourses(master_idx).b],2));
        
        % Allow for delay in surround 
        % time_course{2} = [time_course{1}(1:(end-surround_offset)); zeros(surround_offset,1)];
        
        for i = 1:length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            
            % load up the masks
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
            
            % save the PSTH for each condition to compare to the last condition later 
            PSTH{2}{i} = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 1);
            if s==10
                second_half_PSTH = IDP_plot_PSTH(second_half_condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 1);
            end
            
            % organize the masks and movies
            comp_mask = mod(mask+1,2);
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
            comp_mask = repmat(comp_mask, 1, 1, 1200);
            movie{1} = NSmovie;
            movie{2} = NSmovie .* mask;
            movie{3} = NSmovie.*comp_mask;
            
            % find the generator signals
            
            % Fixed CS weights
            %%{
            for i_movie_type = 2
                spatial_GS = squeeze(convn(movie{i_movie_type}, spatial{1}-sta_fit.surround_amp_scale*spatial{2}, 'valid'));
                gen_signal{i_movie_type} = conv(spatial_GS, time_course{1}, 'valid');
            end
            model_arch = 'y~b1/(b2+exp(b3*x1))';

            % Refit center-surround weights
            %{
                for i_sta_part = 1:2
                    for i_movie_type = 1:3
                        spatial_GS = squeeze(convn(movie{i_movie_type}, spatial{i_sta_part}, 'valid'));
                        gen_signal{i_movie_type}{i_sta_part} = conv(spatial_GS, time_course{i_sta_part}, 'valid');
                    end
                end
                model_arch = 'y~b1/(b2+exp(b3*(x1+b4*x2)))';
                init = [opt_model.model.Coefficients.Estimate; -sta_fit.surround_amp_scale];
            %}
            
            for i_movie_type = 2
                % "Fixed" model
                %{
                    x1 = gen_signal{i_movie_type}{1};
                    x2 = gen_signal{i_movie_type}{2};
                    x = x1+init(4)*x2;
                    model_PSTH{i_movie_type} = predict(opt_model.model, x);
                %}
                
                % Refit center-surround weights
                %data_matrix = [gen_signal{i_movie_type}{1}(1:1166) gen_signal{i_movie_type}{2}(1:1166)];
                
                % Fixed CS weights
                data_matrix = gen_signal{i_movie_type}(1:1166);
                
                % Fit the model
                responses = PSTH{i_movie_type}{i}(30:end);
                model = fitnlm(data_matrix, responses, model_arch, init);
                
                % use the firing rate 
                model_PSTH{i_movie_type}{i}= predict(model, data_matrix);
                % use poisson spiking
                %{
                firing_rate = predict(model, data_matrix);
                logical_sim = Poisson_spiking(firing_rate, 15, 10, monitor_refresh);
                model_PSTH{i_movie_type}{i} = imresize(sum(logical_sim), [1 1166], 'box')*monitor_refresh/1.5;
                %}
                
                % Model performance
                Corr{i_exp}(i_cell,i,i_movie_type) = model.Rsquared.Ordinary;
                
                % Calculate the Noise Floor
                if s == 10
                    responses = PSTH{2}{i}(30:end);
                    model = fitnlm(data_matrix, responses, model_arch, init);
                    % use the firing rate 
                    model_second_half_PSTH= predict(model, data_matrix);
                    % use poisson spiking
                    %{
                    firing_rate= predict(model, data_matrix);
                    logical_sim = Poisson_spiking(firing_rate, 15, 10, monitor_refresh);
                    model_second_half_PSTH = imresize(sum(logical_sim), [1 1166], 'box')*monitor_refresh/1.5;
                    %}
                end
            end
            
            %model_surround_struct(i_cell,i) = var(model_PSTH{3});
        end
        for i = 1:(length(mask_conditions{subgroup}))
            NSEM_Corr(i_cell,i) = err(second_half_PSTH, PSTH{2}{i});
            model_NSEM_Corr(i_cell,i) = err(model_second_half_PSTH, model_PSTH{2}{i});
        end
    end
end

%% Plotting

close all
default_colors = get(gca,'ColorOrder');
all_sigmas = [2 4 6 8 10];
mean_corr = mean(NSEM_Corr([1:3 5:6],:));
noise_floor = mean_corr(end);
mean_corr = mean_corr - noise_floor;
std_corr =  std(NSEM_Corr([1:3 5:6],:))/sqrt(5);
figure(1);
error_bars_y = [std_corr+mean_corr flip(mean_corr-std_corr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*[1 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

%%{
mean_corr = mean(model_NSEM_Corr([1:3 5:6],:));
noise_floor = mean_corr(end);
mean_corr = mean_corr - noise_floor;
std_corr =  std(model_NSEM_Corr([1:3 5:6],:))/sqrt(5);
figure(2); hold on;
error_bars_y = [std_corr+mean_corr flip(mean_corr-std_corr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*default_colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])
%}

%% Plotting ON V OFF

close all
default_colors = get(gca,'ColorOrder');
all_sigmas = [2 4 6 8 10];
mean_corr_ON = mean(NSEM_Corr(1:3,:));
noise_floor = mean_corr(end);
mean_corr = mean_corr - noise_floor;
std_corr =  std(NSEM_Corr(1:3,:))/sqrt(5);
figure(1);
error_bars_y = [std_corr+mean_corr flip(mean_corr-std_corr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*[1 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

%%{
mean_corr = mean(model_NSEM_Corr([1:3 5:6],:));
noise_floor = mean_corr(end);
mean_corr = mean_corr - noise_floor;
std_corr =  std(model_NSEM_Corr([1:3 5:6],:))/sqrt(5);
figure(2); hold on;
error_bars_y = [std_corr+mean_corr flip(mean_corr-std_corr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*default_colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])
%}