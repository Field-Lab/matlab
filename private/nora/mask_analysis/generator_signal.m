clear 
params_201602178
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

n_subgroups = length(cell_idx);
GS = 1;
sta_scale = 4;
default_colors = get(gca,'ColorOrder');


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

%% load movie
if GS
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
    
end


%%
% plot the PSTH for each condition for each cell

for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        [reg, time] = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0);
        %sta = squeeze(sum(datarun_class.stas.stas{master_idx},3));
        %%{
        sta_fit = datarun_class.vision.sta_fits{master_idx}; 
        ideal_center = make_Gaussian_two_d('center_point_x', sta_fit.mean(1), 'center_point_y', 40 -sta_fit.mean(2), 'sd_x', sta_fit.sd(1), 'sd_y', sta_fit.sd(2), 'rotation_angle', sta_fit.angle, 'x_dim', 80, 'y_dim', 40);
        ideal_surround = make_Gaussian_two_d('center_point_x', sta_fit.mean(1), 'center_point_y', 40-sta_fit.mean(2), 'sd_x', 2*sta_fit.sd(1), 'sd_y',2*sta_fit.sd(2), 'rotation_angle', sta_fit.angle, 'x_dim', 80, 'y_dim', 40);
        time_course=mean([datarun_class.vision.timecourses(master_idx).r ...
            datarun_class.vision.timecourses(master_idx).g ...
            datarun_class.vision.timecourses(master_idx).b],2);
        sta_center = reshape((ideal_center(:))*time_course', [40 80 30]);
        sta_surround = reshape((ideal_surround(:))*time_course', [40 80 30]);
        %}
        sta_center = flip(sta_center,1);  
        sta_center = flip(sta_center,2); 
        sta_center = flip(sta_center,3); 
        sta_surround = flip(sta_surround,1);  
        sta_surround = flip(sta_surround,2); 
        sta_surround = flip(sta_surround,3); 
        for i = 1:length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            if s==2
                load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_allcells_sigma' num2str(s) '.mat'])
            else
                load(['/Volumes/Data/2016-02-17-1/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma' num2str(s) '.mat']);
            end
            comp_mask = mod(mask+1,2);  
            mask = imresize(mask, 1/sta_scale, 'box');
            mask = repmat(mask, 1, 1, 1200);
            comp_mask = imresize(comp_mask, 1/sta_scale, 'box');
            comp_mask = repmat(comp_mask, 1, 1, 1200);
            mask_movie = NSmovie .* mask;
            mask_gen_signal_center = squeeze(convn(mask_movie, sta_center, 'valid')); 
            mask_gen_signal_surround = squeeze(convn(mask_movie, sta_surround, 'valid')); 
            comp_gen_signal_center = squeeze(convn(NSmovie.*comp_mask, sta_center, 'valid')); 
            comp_gen_signal_surround = squeeze(convn(NSmovie.*comp_mask, sta_surround, 'valid')); 
            mask = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0);
            comp = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0);
            close all
            figure(1); hold on
            plot(time(1:1166), mask(30:end), 'Color', default_colors(1,:));
            plot(time(1:1171), mask_gen_signal_center-0.25*mask_gen_signal_surround, 'Color', default_colors(2,:)); 
            xlim([0 8])
            %exportfig(gcf, [fig_save '/mask_gen_' 'cell' num2str(i_cell) '_sigma_' num2str(sigmas(mask_conditions{subgroup}(i)))], 'Bounds', 'loose', 'Color', 'rgb')

            figure(2); hold on
            plot(time(1:1166), comp(30:end), 'Color', default_colors(1,:));
            plot(time(1:1171), comp_gen_signal_center-0.25*comp_gen_signal_surround, 'Color', default_colors(2,:));
            xlim([0 8])
            %exportfig(gcf, [fig_save '/comp_gen_' 'cell' num2str(i_cell) '_sigma_' num2str(sigmas(mask_conditions{subgroup}(i)))], 'Bounds', 'loose', 'Color', 'rgb')

            figure(3); hold on
            plot(time(1:1166), reg(30:end), 'k'); 
            plot(time(1:1166), comp(30:end), 'Color', default_colors(3,:));
            plot(time(1:1171), mask_gen_signal_center-0.25*mask_gen_signal_surround, 'Color', default_colors(1,:)); 
            xlim([0 8])
            %exportfig(gcf, [fig_save '/mask_suppress_' 'cell' num2str(i_cell) '_sigma_' num2str(sigmas(mask_conditions{subgroup}(i)))], 'Bounds', 'loose', 'Color', 'rgb')
             pause()
            close all
        end
    end
end 
close all