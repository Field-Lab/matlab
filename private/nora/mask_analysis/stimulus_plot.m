clear 
params_201602171
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 8 4])

n_subgroups = length(cell_idx);
GS = 0;
sta_scale = 4;

%%
colors = [0    0    1;...
    0    0    0.3;...
          0.6    0    0.6;...
          1    0    0];
      colors = flip(colors, 1);
%figure; hold on; for i=1:4; plot([0 1], [0 i], 'Color', colors(i,:)); end
%%

the_chosen_one =2;
sub_TCO = 1;

%%
% load classification run 
datarun_class = load_data(class);
datarun_class = load_params(datarun_class);
datarun_class = load_sta(datarun_class, 'load_sta', cells_orig);

%% load movie
if 1
    for i_chunk = 7
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NSmovie = movie_chunk;
    clear movie_chunk
    NSmovie = NSmovie - 64;
    NSmovie = permute(NSmovie, [2 1 3]);
    end
end

%%

for subgroup = sub_TCO%1:n_subgroups
    for i_cell = the_chosen_one%cell_idx{subgroup}
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        sta = squeeze(sum(datarun_class.stas.stas{master_idx},3));
        [U,~,~] = svd(reshape(sta, [40*80, 30]));
        sta = reshape(U(:,2), 40, 80);
        imagesc(sta)
        for i = [1 2 4]%length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            the_fit = datarun_class.stas.fits{master_idx};
            ctr = the_fit.mean;
            rad = s*the_fit.sd;
            hold on; [X,Y] = drawEllipse_upsampled([ctr rad the_fit.angle]);
            hold on; plot(X,Y, 'Color', colors(i,:), 'LineWidth', 2)
        end
        colormap gray
        axis image
        axis off
        xlim([31 51]);ylim([2 22])
       exportfig(gcf, [fig_save '/sta_' 'cell' num2str(i_cell)], 'Bounds', 'loose', 'Color', 'rgb')
        %close all
        
    end
    figure;
    load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma4.mat']);
    %load('/Volumes/Data/2016-04-21-1/Visual/2016-04-21-1_NJB_Masks/Maskin_allcells_sigma4.mat');  

    mask_frame = NSmovie(:,:,2).*mask;
    imagesc(mask_frame)
    axis image
    colormap gray
    axis off
    caxis([0 255]-64)
    xlim([31 51]*4);ylim([2 22]*4)
    exportfig(gcf, [fig_save '/spot_plot_sigma4'], 'Bounds', 'loose', 'Color', 'rgb')
    
    figure;
    comp_mask = mod(mask+1,2);
    mask_frame = NSmovie(:,:,2).*comp_mask;
    imagesc(mask_frame)
    axis image
    colormap gray
    axis off
    caxis([0 255]-64)
    xlim([31 51]*4);ylim([2 22]*4)
   exportfig(gcf, [fig_save '/gap_plot_sigma4'], 'Bounds', 'loose', 'Color', 'rgb')
end

pause();
close all

%%
for i_reg = 1:n_reg
    datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
    datarun{i_reg} = load_neurons(datarun{i_reg});
    reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 0);
    clear datarun
end
reg_data{1}.testspikes = reg_data{1}.testspikes([1:6 8:end],:);
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
%close all

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
% plot the PSTH for each condition for each cell
set(0, 'defaultFigurePosition', [2 2 10 3])
for subgroup = sub_TCO%:2
    for i_cell = the_chosen_one%cell_idx{subgroup} 
        for i = [1 2 4]%length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i)); 
            figure(1); hold on
            IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', colors(i,:), 'smoothing', 20);
            figure(2); hold on
            IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', colors(i,:), 'smoothing', 20);
        end
        exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) '_mask'], 'Bounds', 'loose', 'Color', 'rgb')

        figure(1); hold on; IDP_plot_PSTH(reg_data{1}, i_cell, 'color', colors(3,:), 'plot_offset', 0, 'smoothing', 20);
        axis([-0.1 8 -250 250])
        axis off
        hold on; plot([-0.01 -0.01], [0 100], 'k', 'LineWidth', 3)
        exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) '_mask'], 'Bounds', 'loose', 'Color', 'rgb')
        figure(2); hold on; IDP_plot_PSTH(reg_data{1}, i_cell, 'color', colors(3,:), 'plot_offset', 0, 'smoothing', 20, 'negative', 1);
        axis([-0.1 8 -250 250])
        hold on; plot([-0.01 -0.01], [0 100], 'k', 'LineWidth', 3)
        axis off
        exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) '_comp'], 'Bounds', 'loose', 'Color', 'rgb')
        
        disp(i_cell);%pause();
        figure(1); clf; figure(2); clf;
    end
end
%close all