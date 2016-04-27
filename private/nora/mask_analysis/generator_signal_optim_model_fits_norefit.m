clear
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 12 4])
GS = 1;
sta_scale = 4;
colors = colormap(winter(6));

TCF = [];
TMG = [];
TMF = [];
TCG = [];

params_list
eval(exps{end})
n_reg = length(reg);
mkdir(fig_save);
n_subgroups = length(cell_idx);
close all

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
% avg RF
%datarun_class = load_sta(datarun_class);


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
    NSmovie = NSmovie - 64;
    NSmovie = imresize(NSmovie, 1/sta_scale, 'box');
    NSmovie = permute(NSmovie, [2 1 3]);
    
end

interval_frame = cumsum(interval);

%%
% plot the PSTH for each condition for each cell

glm_files = dir('/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1/NSEM_full_opt/*NSEM_optmodel.mat');
glm_fit_idx = [6 3 5 1 2 4];
%%
% 1 = center, 2 = surround
% 1 = reg, 2 = spot, 3 = gap
for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        disp(cells_orig(i_cell))
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        [PSTH{1}, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0, 'smoothing', 10);
        
        load(['/Volumes/Lab/Users/Nora/GLMFits_masking/2016-04-21-1/NSEM_full_opt/' glm_files(glm_fit_idx(i_cell)).name]);
        spatial = opt_model.orig_glm.linearfilters.Stimulus.space_rk1;
        ROIcoord.xvals = opt_model.orig_glm.linearfilters.Stimulus.x_coord;
        ROIcoord.yvals = opt_model.orig_glm.linearfilters.Stimulus.y_coord;
        time_course= opt_model.new_time;
        %spatial = flip(spatial,1);
        %spatial = flip(spatial,2);
        time_course = flip(time_course);
        
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
            
            % organize the mask movies
            comp_mask = mod(mask+1,2);
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
            comp_mask = repmat(comp_mask, 1, 1, 1200);
            movie{1} = NSmovie;
            movie{2} = NSmovie .* mask;
            movie{3} = NSmovie.*comp_mask;
            
            for i = 1:3
                movie{i} = movie{i}(ROIcoord.xvals ,ROIcoord.yvals ,:);
            end
                
                % find the generator signals
            for i_movie_type = 1:3
                gen_signal{i_movie_type} = conv(squeeze(convn(movie{i_movie_type}, spatial, 'valid')), time_course, 'full');
            end
            
            PSTH{2} = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            PSTH{3} = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            
            for i_movie_type = 1:3
                model = fitnlm(gen_signal{i_movie_type}(1:1166), PSTH{i_movie_type}(30:end), 'y~b1./(b2+exp(b3*x1))', opt_model.model.Coefficients.Estimate)
                figure; plot(predict(model, gen_signal{i_movie_type}(1:1166))); hold on; plot(PSTH{i_movie_type}(30:end)); pause()
                Corr(i_cell, i,i_movie_type, :) = model.Rsquared.Ordinary;
            end

        end
    end
end

            
            
            
            %{
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
            plot(time(30:end),comp(30:end))
            hold on;
            plot(time(30:end), mask_gen_signal(1:1166))
            xlim([1 9])
            axis off
            exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_compVSmaskGS'], 'Bounds', 'loose', 'Color', 'rgb')
            close all
            
            plot(time(30:end),reg(30:end))
            hold on;
            plot(time(30:end),mask_PSTH(30:end))
            plot(time(30:end), comp_gen_signal(1:1166))
            xlim([1 9])
            axis off
            exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_maskVScompGS'], 'Bounds', 'loose', 'Color', 'rgb')
            %}
    
    %{
    plot([2 4 6 8 10], Corr([1:3 5:6],:,2)', 'Color', [1 1 1]*0.75)
    hold on;
    plot([2 4 6 8 10], mean(Corr([1:3 5:6],:,2)), 'Color', [1 1 1]*0, 'LineWidth', 2)
    xlabel('Spot Size (sigmas)')
    ylabel('Model Performance (R^2)')
    %}
 

    %{ 
for i_movie_type = 1:3
                
                data_matrix = [gen_signal{i_movie_type}{1}(1:1166) gen_signal{i_movie_type}{2}(1:1166)];
                responses = PSTH{i_movie_type}(30:end);
                model = fitnlm(data_matrix, responses, model_arch, init);
                %total_gen_signal{i_movie_type} = model.Coefficients.Estimate(3)*(gen_signal{i_movie_type}{1}(1:1166)+model.Coefficients.Estimate(3)*gen_signal{i_movie_type}{2}(1:1166));
                
                %{
                    %raster plotting
                    plot(predict(model, data_matrix)); hold on; plot(responses);
                    %pause()
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
                %}
                
                % compare model performance
                Corr(i_cell, i,i_movie_type, :) = model.Rsquared.Ordinary;
                
                
                %exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) 'sigma' num2str(i) '_raster_' num2str(i_movie_type)], 'Bounds', 'loose', 'Color', 'rgb')
                close all
                %}
