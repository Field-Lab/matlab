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

for i_exp = [5]
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
            i_cell
            disp(cells_orig(i_cell))
            master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
            [PSTH{1}, time] = IDP_plot_PSTH(reg_data{n_reg},i_cell, 'color', 0, 'smoothing', 1);
            %%{
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
            if strcmp(cell_type{subgroup}, 'on')
                time_course{1}= -opt_model.new_time;
            else
                time_course{1}= opt_model.new_time;
            end
            
            %time_course{1}=flip(mean([datarun_class.vision.timecourses(master_idx).r ...
            %    datarun_class.vision.timecourses(master_idx).g ...
            %    datarun_class.vision.timecourses(master_idx).b],2));
            time_course{2} = [time_course{1}(1:(end-surround_offset)); zeros(surround_offset,1)];
            %}
            for i = 1:length(mask_conditions{subgroup})
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
                
                % spot size analysis code
                PSTH{2} = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 1);
                PSTH{3} = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 1);
                NSEM_Corr(i_cell,i) = err(PSTH{1}, PSTH{2});
                surround_struct(i_cell,i) = var(PSTH{3});
                
                %%{
                
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

                model_arch = 'y~b1/(b2+exp(b3*(x1+b4*x2)))';
                init = [opt_model.model.Coefficients.Estimate; -sta_fit.surround_amp_scale];
                for i_movie_type = 1:3
                    % get the correlation and the generator signal for
                    % later
                    x1 = gen_signal{i_movie_type}{1};
                    x2 = gen_signal{i_movie_type}{2};
                    x = x1+init(4)*x2;
                    model_PSTH{i_movie_type} = predict(opt_model.model, x);
                    
                    %data_matrix = [gen_signal{i_movie_type}{1}(1:1166) gen_signal{i_movie_type}{2}(1:1166)];
                    %responses = PSTH{i_movie_type}(30:end);
                    %model = fitnlm(data_matrix, responses, model_arch, init);
                    %model_PSTH{i_movie_type} = predict(model, data_matrix);
                    %Corr{i_exp}(i_cell, i,i_movie_type) = model.Rsquared.Ordinary;
                    
                    
                    %total_gen_signal{i_movie_type} = model.Coefficients.Estimate(3)*(gen_signal{i_movie_type}{1}(1:1166)+model.Coefficients.Estimate(4)*gen_signal{i_movie_type}{2}(1:1166));
                    %Corr{exp}(i_cell, i,i_movie_type) = model.Rsquared.Ordinary;
                    %{
                    total_gen_signal{i_movie_type} = -model.Coefficients.Estimate(3)*(gen_signal{i_movie_type}{1}(1:1166)+model.Coefficients.Estimate(4)*gen_signal{i_movie_type}{2}(1:1166));
                    PSTH{2} = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
                    PSTH{3} = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
                    model_arch = 'y~b1/(b2+exp(b3*(x1+b4*x2)))';
                    init = [opt_model.model.Coefficients.Estimate; -sta_fit.surround_amp_scale];
                    firing_rate = predict(model, data_matrix);
                    switch i_movie_type
                        case 1
                            IDP_plot_raster(reg_data{1},i_cell, firing_rate);
                        case 2
                            IDP_plot_raster(condition{mask_conditions{subgroup}(i)},i_cell, firing_rate);
                        case 3
                            IDP_plot_raster(condition{comp_conditions{subgroup}(i)},i_cell, firing_rate);
                    end
                   
                    axis off
                    % compare model performance
                    %sum((responses-predict(model, data_matrix)).^2);
                    %}
                    
                    %exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_raster_' num2str(i_movie_type)], 'Bounds', 'loose', 'Color', 'rgb')
                    %close all
                end
                model_NSEM_Corr(i_cell,i) = err(model_PSTH{1}, model_PSTH{2});
                model_surround_struct(i_cell,i) = var(model_PSTH{3});
                % compare model performance
                %sum((responses-predict(model, data_matrix)).^2);
                %exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_raster_' num2str(i_movie_type)], 'Bounds', 'loose', 'Color', 'rgb')
                %close all
                %}
            end
            
            
            %{
                % Generator signals correlated vs anticorrelated
            start = 1;
            for i_interval = 1:24
               comp_firing(i_interval) = sum(PSTH{3}(start:interval_frame(i_interval)-1));
               mask_firing(i_interval) = sum(PSTH{2}(start:interval_frame(i_interval)-1));
               try
                   mask_gen(i_interval) = sum(total_gen_signal{2}((start:interval_frame(i_interval)-1)-30));
                   comp_gen(i_interval) = sum(total_gen_signal{3}((start:interval_frame(i_interval)-1)-30));
               catch
                   mask_gen(i_interval) = 0;
                   comp_gen(i_interval) = 0;
               end
               start = interval_frame(i_interval);
            end
            TCF = [TCF; comp_firing];
            TMF = [TMF; mask_firing];
            TMG = [TMG; mask_gen];
            TCG = [TCG; comp_gen];
            %}
            %{
                default_colors = get(gca,'ColorOrder');
                figure;
            plot(time(30:end),PSTH{3}(30:end))
            hold on;
            plot(time(30:end), 10*total_gen_signal{2}(1:1166))
            xlim([1 3])
            axis off
            plot(time(30:end),PSTH{1}(30:end), 'Color', default_colors(4,:))
            set(gcf,'Position', [2 2 4 4])
            exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_compVSmaskGS'], 'Bounds', 'loose', 'Color', 'rgb')
            %close all
            
            figure;
            hold on;
            plot(time(30:end),PSTH{2}(30:end))
            plot(time(30:end), 10*total_gen_signal{3}(1:1166))
            plot(time(30:end),PSTH{1}(30:end), 'Color', default_colors(4,:))
            xlim([1 3])
            axis off
            set(gcf,'Position', [2 2 4 4])
            exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_maskVScompGS'], 'Bounds', 'loose', 'Color', 'rgb')
            pause()
            close all
            %}
            
        end
    end
        sigma = sigmas(mask_conditions{subgroup});
    
    for i = 1:length(cells_masking)
        MSE =  NSEM_Corr(i,:);
        surr = surround_struct(i,:);
        if ~any(diff(MSE)>20) && ~any(diff(surr)>100)

            figure(1); hold on; plot(sigma, MSE, 'Color', [1 1 1]*0.75); hold on
            figure(2); hold on; plot(sigma, surr,  'Color', [1 1 1]*0.75); hold on
            if strcmp(piece, '2016-04-21-1') && i == 1
                example_Cell = {sigma, MSE, surr};
            end
            for s = 1:length(sigma)
                mean_cell_corr{sigma(s)} = [mean_cell_corr{sigma(s)} MSE(s)];
                mean_cell_surr{sigma(s)} = [mean_cell_surr{sigma(s)} surr(s)];
            end

        end
    end
    for i = 1:length(cells_masking)
        MSE =  model_NSEM_Corr(i,:);
        surr = model_surround_struct(i,:);
        if ~any(diff(MSE)>20) && ~any(diff(surr)>100)

            figure(3); hold on; plot(sigma, MSE, 'Color', [1 1 1]*0.75); hold on
            figure(4); hold on; plot(sigma, surr,  'Color', [1 1 1]*0.75); hold on
            if strcmp(piece, '2016-04-21-1') && i == 1
                example_Cell = {sigma, MSE, surr};
            end
            for s = 1:length(sigma)
                model_mean_cell_corr{sigma(s)} = [model_mean_cell_corr{sigma(s)} MSE(s)];
                model_mean_cell_surr{sigma(s)} = [model_mean_cell_surr{sigma(s)} surr(s)];
            end

        end
    end
end


%%
close all
default_colors = get(gca,'ColorOrder');
all_sigmas = [2 4 5 6 8 10];
for s = 1:6
   mean_corr(s) = mean(mean_cell_corr{all_sigmas(s)}); 
   mean_surr(s) = mean(mean_cell_surr{all_sigmas(s)}); 
   std_corr(s) = std(mean_cell_corr{all_sigmas(s)})/sqrt(length(mean_cell_corr{all_sigmas(s)})); 
   std_surr(s) = std(mean_cell_surr{all_sigmas(s)})/sqrt(length(mean_cell_surr{all_sigmas(s)})); 
end
figure(1); 
error_bars_y = [std_corr+mean_corr flip(mean_corr-std_corr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*[1 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

figure(2); 
error_bars_y = [std_surr+mean_surr flip(mean_surr-std_surr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*[1 1 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_surr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

for s = 1:6
   mean_corr(s) = mean(model_mean_cell_corr{all_sigmas(s)}); 
   mean_surr(s) = mean(model_mean_cell_surr{all_sigmas(s)}); 
   std_corr(s) = std(model_mean_cell_corr{all_sigmas(s)})/sqrt(length(model_mean_cell_corr{all_sigmas(s)})); 
   std_surr(s) = std(model_mean_cell_surr{all_sigmas(s)})/sqrt(length(model_mean_cell_surr{all_sigmas(s)})); 
end
figure(1); hold on;
error_bars_y = [std_corr+mean_corr flip(mean_corr-std_corr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*default_colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

figure(2); hold on;
error_bars_y = [mean_surr+std_surr flip(mean_surr-std_surr)];
error_bars_x = [all_sigmas flip(all_sigmas)];
fill(error_bars_x, error_bars_y, 0.75*default_colors(1,:), 'FaceAlpha', 0.25, 'EdgeColor', 'none');
hold on; plot(all_sigmas, mean_surr, 'k', 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

%% Plot Correlation across spot size
%{
figure; hold on;
plot([2 4 5 6 12], [squeeze(Corr{2}(:,:,2)) squeeze(Corr{2}(:,1,1))]','-', 'Color', [1 1 1]*0.75);
plot([2 4 5 6 12], [squeeze(Corr{3}(:,:,2)) squeeze(Corr{3}(:,1,1))]','-', 'Color', [1 1 1]*0.75);
plot([2 4 6 8 10 12], [squeeze(Corr{5}([1:3 5:6],1:5,2)) squeeze(Corr{5}([1:3 5:6],1,1))]','-', 'Color', [1 1 1]*0.75);
%i_cell = 5; plot([2 4 5 6 12], [squeeze(Corr{2}(i_cell,:,2)) squeeze(Corr{2}(i_cell,1,1))]','-', 'Color', default_colors(1,:), 'LineWidth', 2);

clear mean_corr
for sigma = [2 4 6]
    mean_corr(1:3) = mean(cat(1, Corr{2}(:,[1 2 4],2),Corr{3}(:,[1 2 4],2),Corr{5}([1:3 5:6],[1 2 3],2)));
    %mean_corr(4:5) = mean(Corr{5}([1:3 5:6], [4:5], 2));
    mean_corr(4) = mean(cat(1, Corr{2}(:,1,1),Corr{3}(:,1,1),Corr{5}([1:3 5:6],1,1)));
end

plot([2 4 6 12], mean_corr, 'k.-', 'LineWidth', 2, 'MarkerSize', 10)
set(gcf, 'Position', [2 2 5 6])
set(gca, 'XTick', [2 4 6 8 10 12])
set(gca, 'XTickLabels', {'2', '4', '6', '8', '10', 'Full'})
set(gca, 'YTick', 0.4:0.2:0.8)
%exportfig(gcf, [fig_save '/Corr_across_spot_size.eps'], 'Bounds', 'loose', 'Color', 'rgb')
%}