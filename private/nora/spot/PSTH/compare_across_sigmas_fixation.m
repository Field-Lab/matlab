clear
total_mask = [];
total_reg = [];
celltype_idx = [];

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 8 4])

params_list
total_firing = [];
total_reg_firing = [];

%%
for exp = exps
    disp(exp{1})
    eval(exp{1})
    n_reg = length(reg);
    mkdir(fig_save);
    
    
    n_subgroups = length(cell_idx);

%%
%%{
% load reg repeats
for i_reg = 1:n_reg
    datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
    datarun{i_reg} = load_neurons(datarun{i_reg});
    reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 1);
    clear datarun
end
%}

%{
    % plot reg repeats for stability check
    for i_cell = 1:size(reg_data{1}.testspikes,2)
        figure; hold on
        for i_reg = 1:n_reg
            IDP_plot_PSTH(reg_data{i_reg}, i_cell, 'color', i_reg);
        end
        title(['cell ' num2str(i_cell)]);
    end
%}
close all

%%
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
movie = [];
for i = 1:24
    load(['/Volumes/Lab/Users/Nora/new_stim_nora/NSEM intervals/matfiles/movie_chunk_' num2str(i) '.mat'])
    interval(i) = size(movie_chunk, 3);
    movie = cat(3, movie, movie_chunk);
end
default_colors = get(gca,'ColorOrder');
interval_frame = cumsum(interval);

%%
% plot the PSTH for each condition for each cell
start = 1;
for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        reg_psth = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0, 'smoothing', 10);
        for i = 1:length(mask_conditions{subgroup})
            mask = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0, 'smoothing', 10);
            for fixation = 1:23
                Firing(i_cell, fixation, i) = max(mask(start:interval_frame(fixation)));
                reg_firing(i_cell, fixation) =  max(reg_psth(start:interval_frame(fixation)));
                start = interval_frame(fixation);
            end
            start = 1;
        end
    end
end

total_firing = [total_firing; Firing(:,:,1)];
total_reg_firing = [total_reg_firing; reg_firing];


end

plot(total_firing, total_reg_firing, 'k.')
hold on; plot([0 300], [0 300])

%{
%%
idx = -10:10;
outer_idx = -100:100;
for i_cell = cell_idx{1}
    datarun_class = load_data([ '/Volumes/Analysis/' piece '/streamed/' class '/' class]);
    datarun_class = load_params(datarun_class);
end
for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        index = get_cell_indices(datarun_class, cells_orig(i_cell));
        center = round(datarun_class.vision.sta_fits{index}.mean);
        x_ind = center(1)*4+idx;
        x_ind = x_ind(x_ind>0 & x_ind<320);
        y_ind = (40 - center(2))*4+idx;
        y_ind = y_ind(y_ind>0 & y_ind<160);
        outer_x_ind = center(1)*4+outer_idx;
        outer_x_ind = outer_x_ind(outer_x_ind>0 & outer_x_ind<=320);
        outer_y_ind = (40 - center(2))*4+outer_idx;
        outer_y_ind = outer_y_ind(outer_y_ind>0 & outer_y_ind<=160);
        start = 1;
        for fixation = 1:23
            image = movie(:,:,start+1); 
            image(x_ind, y_ind) = 0;
            image = image(outer_x_ind, outer_y_ind);
            fixation_contrast(fixation) = std(image(:));
            start = interval_frame(fixation);
        end
        plot(diff_by_fix(i_cell, :), fixation_contrast, '.');
        pause()
    end
end
%}
