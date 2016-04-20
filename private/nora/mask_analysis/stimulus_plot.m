clear 
params_201602178
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

n_subgroups = length(cell_idx);
GS = 0;
sta_scale = 4;
default_colors = get(gca,'ColorOrder');


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

for subgroup = 1:n_subgroups
    for i_cell = cell_idx{subgroup}
        master_idx = get_cell_indices(datarun_class, cells_orig(i_cell));
        sta = squeeze(sum(datarun_class.stas.stas{master_idx},3));
        [U,~,~] = svd(reshape(sta, [40*80, 30]));
        sta = reshape(U(:,2), 40, 80);
        imagesc(sta)
        for i = 1:length(mask_conditions{subgroup})
            s = sigmas(mask_conditions{subgroup}(i));
            the_fit = datarun_class.stas.fits{master_idx};
            ctr = the_fit.mean;
            rad = s*the_fit.sd;
            hold on; [X,Y] = drawEllipse_upsampled([ctr rad the_fit.angle]);
            hold on; plot(X,Y, 'k')
        end
        colormap gray
        axis image
        axis off
        exportfig(gcf, [fig_save '/mask_plot_' 'cell' num2str(i_cell)], 'Bounds', 'loose', 'Color', 'rgb')
        close all
        
    end
    load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells' num2str(subgroup) '_sigma4.mat']);
    mask_frame = NSmovie(:,:,2).*mask;
    imagesc(mask_frame)
    axis image
    colormap gray
    axis off
    caxis([0 255]-64)
    exportfig(gcf, [fig_save '/spot_plot_sigma4'], 'Bounds', 'loose', 'Color', 'rgb')
    comp_mask = mod(mask+1,2);
    mask_frame = NSmovie(:,:,2).*comp_mask;
    imagesc(mask_frame)
    axis image
    colormap gray
    axis off
    caxis([0 255]-64)
    exportfig(gcf, [fig_save '/gap_plot_sigma4'], 'Bounds', 'loose', 'Color', 'rgb')
end
close all