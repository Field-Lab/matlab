clear 
params_201602176_A
n_reg = length(reg);
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 4 4])

n_subgroups = length(cell_idx);
GS = 0;
sta_scale = 4;

%%
datarun_class = load_data(class);
datarun_class = load_params(datarun_class);
datarun_class = load_neurons(datarun_class);

   %%
%{
movie = load_java_movie(datarun_class, '/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-119.5.xml', datarun_class.triggers);

for i = 1:3
frame = movie.stimulus.java_movie.getFrame(i).getBuffer;
whole_Frame = reshape(frame, 40, 80, 3);
image(whole_Frame(1:10, 1:10, :));
axis square
axis off
%exportfig(gcf, ['/Volumes/Lab/Users/Nora/white_noise_stim_frame' num2str(i) '.eps'], 'Bounds', 'tight', 'Color', 'rgb', 'Renderer', 'zbuffer')
end
%}
%%
cell= 2; 
datarun_class = load_sta(datarun_class, 'load_sta', cells_orig(cell), 'save_rf', true);
plot_rf(datarun_class, cells_orig(cell), 'fit', false);

axis([28 53 15 40])
axis off
%exportfig(gcf, ['ex_sta.eps'], 'Bounds', 'tight', 'Color', 'rgb', 'Renderer', 'zbuffer')

load(['/Volumes/Lab/Users/Nora/GLMFits_masking/' piece '/NSEM_full_opt/' num2str(cells_orig(cell)) 'NSEM_optmodel.mat']);
sta_fit = opt_model.sta_fit.params;


spatial{1} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', 1.4403, 'sd_y',  1.1435, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);
spatial{2} = make_Gaussian_two_d('center_point_x', sta_fit.center_point_x, 'center_point_y', sta_fit.center_point_y, 'sd_x', sta_fit.surround_sd_scale*sta_fit.center_sd_x, 'sd_y',sta_fit.surround_sd_scale*sta_fit.center_sd_y, 'rotation_angle', sta_fit.center_rotation_angle, 'x_dim', 80, 'y_dim', 40);

spatial_total = spatial{1} - sta_fit.surround_amp_scale*spatial{2};
spatial_total = cat(3, spatial_total*sta_fit.color_weight_a, spatial_total*sta_fit.color_weight_b, spatial_total*sta_fit.color_weight_c);
%%

figure;
image(1.3*(norm_image(spatial_total)-0.1))
axis([28 53 15 40]) 
axis off

%figure; surf(spatial{1} - sta_fit.surround_amp_scale*spatial{2})

%
 for s = [2 4 6 ]
            the_fit = datarun_class.stas.fits{datarun_class.cell_ids == cells_orig(cell)};
            ctr = the_fit.mean;
            rad = s*the_fit.sd;
            hold on; [X,Y] = drawEllipse_upsampled([ctr rad the_fit.angle]);
            hold on; plot(X,Y, 'k', 'LineWidth', 2)
 end
        
 exportfig(gcf, ['ex_stafit_with_ellipse.eps'], 'Bounds', 'loose', 'Color', 'rgb', 'Renderer', 'zbuffer')
 
 %%
 figure;
 load(['/Volumes/Data/' piece '/Visual/masks/Maskin/Maskin_cells1_sigma4.mat']);
 mask(1:70,:) = 0;
 mask_frame = NSmovie(:,:,2).*mask;
 imagesc(mask_frame)
 axis image
 colormap gray
 axis off
 caxis([0 255]-64)
axis([28 53 15 40]*4)
%exportfig(gcf, ['ex_sta_spot_plot.eps'], 'Bounds', 'loose', 'Color', 'rgb')
 
 mask = mod(2, mask);
 mask_frame = NSmovie(:,:,2).*mask;
 imagesc(mask_frame)
 axis image
 colormap gray
 axis off
 caxis([0 255]-64)
axis([28 53 15 40]*4)
%exportfig(gcf, ['ex_sta_gap_plot.eps'], 'Bounds', 'loose', 'Color', 'rgb')

 mask_frame = NSmovie(:,:,2);
 imagesc(mask_frame)
 axis image
 colormap gray
 axis off
 caxis([0 255]-64)
axis([28 53 15 40]*4)
%exportfig(gcf, ['ex_sta_FF.eps'], 'Bounds', 'loose', 'Color', 'rgb')


     