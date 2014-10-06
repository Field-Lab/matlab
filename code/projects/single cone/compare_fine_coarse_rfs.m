% Compare RFs in coarse and fine stixel runs


% DATASET SELECTION

switch 1
    case 1 % plantain
        % data003  RGB 1-48-0.48-11111  (320x320) NDF 0.0 7200sec
        % data008 RGB bin 5-16-0.48-11111-64x64 ndf 0.0 xxsec
        fine_data_path = [single_cone_path 'saved/plantain.mat'];
        coarse_data_spec = {'2008-08-27-5/data008/data008/data008'};
        stixel_ratio = 5;

    case 2 % apricot (8x8)
        % data000 Binary RGB 8-16-0.48-11111 ndf 0.0 1800 s
        % data005 Binary RGB 1-32-0.48-11111 ndf 0.0 7200 s
        fine_data_path = [single_cone_path 'saved/apricot.mat'];  %   168 ON m    349 OFF m
        coarse_data_spec = {'2009-04-13-5/data000/data000'};      %   286 ON m     64 OFF m
        stixel_ratio = 8;
        
    case 3 % apple (4x4)
        % data010 Binary RGB 4-16-0.48-11111 ndf 0.0 rpe ~2400 s
        %data013 Binary RGB 1-64-0.48-33333 320x320 ndf 0.0 5082 s
        fine_data_path = [single_cone_path 'saved/apple-13.mat'];
        coarse_data_spec = {'2010-03-05-2/data010/data010'};
        stixel_ratio = 4;
        
end


% PARAMETERS

cell_spec_fine = {3}; % which cells to compute average RF
cell_spec_coarse = cell_spec_fine;
cell_spec_recon = cell_spec_fine;
max_amp = .9; % saturate to show surround
scale_up = 5; % for display only, so pdf looks good
plot_polarity = 1;






% load coarse RFs
if ~exist('datarun_coarse','var')
    % plantain
    datarun_coarse = load_data(coarse_data_spec{:});
    datarun_coarse = load_params(datarun_coarse);
    datarun_coarse = load_sta(datarun_coarse,'load_sta',[]);
    datarun_coarse = get_sta_summaries(datarun_coarse,'all');
end



% load up fine RFs, and bin them
if ~exist('datarun_fine','var')
    tic

    load(fine_data_path)
    datarun_fine = datarun;
    clear datarun;

    cell_ids = get_cell_ids(datarun_fine,'all');

    for cc=1:length(cell_ids)
        cell_id = cell_ids(cc);
        cell_index = get_cell_indices(datarun_fine,cell_id);

        % get center
        center = datarun_fine.stas.rf_coms{cell_index};

        % change to new coordinates
        new_center = 1 + (center - (stixel_ratio+1)/2) / stixel_ratio;
        datarun_fine.stas.rf_coms{cell_index} = new_center;
        
        
        % bin the RF
        rf = matrix_rebinned(get_rf(datarun_fine,cell_id),stixel_ratio);

        % plot binned RF with cetner
        if 0
            figure(3);clf;imagesc(norm_image(rf));axis image;
            hold on; plot(new_center(1),new_center(2),'or')
            title(sprintf('cell id %d, index %d',cell_id,cell_index));pause
        end
        
        % put it in datarun
        datarun_fine.stas.rfs{cell_index} = rf;
    end

    toc

end


% GENERATE STAS BY RECONSTRUCTING CONE WEIGHTS
if ~exist('datarun_recon','var')
    datarun_recon = datarun_fine;
    
    tic 
    
    % load up single cone data
    datarun_recon = import_single_cone_data(datarun_recon);
    load([single_cone_path datarun_recon.names.nickname '/Wc.mat'])
    datarun_recon.cones.mosaic = make_mosaic_struct(datarun_recon.cones.centers);
    
    cell_ids = get_cell_ids(datarun_recon,'all');

    for cc=1:length(cell_ids)
        cell_id = cell_ids(cc);
        cell_index = get_cell_indices(datarun_recon,cell_id);
        
        % bin the RF
        rf_recon = Wc*datarun_recon.cones.weights(:,cell_index);
        rf_recon = reshape(rf_recon,datarun_recon.stimulus.field_height,datarun_recon.stimulus.field_width,3);
        rf = matrix_rebinned(rf_recon,stixel_ratio);
        
        % plot binned RF with cetner
        if 0
            figure(3);clf;imagesc(norm_image(rf));axis image;
            hold on; plot(new_center(1),new_center(2),'or')
            title(sprintf('cell id %d, index %d',cell_id,cell_index));pause
        end
        
        % put it in datarun
        datarun_recon.stas.rfs{cell_index} = rf;
    end

    toc
    
    
    
end


% PLOT COMPARISON

figure(1);clf;
T = @(x)x;


% fine image 
[average_rf_fine,xdata,ydata] = compute_average_rf(datarun_fine,cell_spec_fine);
average_rf_fine = clip_rgb_image(average_rf_fine, max_amp);
subplot(221);imagesc(matrix_scaled_up(norm_image(plot_polarity * average_rf_fine),scale_up),'xdata',T(xdata),'ydata',T(ydata));axis image;
title(sprintf('fine (%d RGCs)',length(get_cell_indices(datarun_fine,cell_spec_fine))))


% coarse image
[average_rf_coarse,xdata,ydata] = compute_average_rf(datarun_coarse,cell_spec_coarse);
average_rf_coarse = clip_rgb_image(average_rf_coarse, max_amp);
subplot(222);imagesc(matrix_scaled_up(norm_image(plot_polarity * average_rf_coarse),scale_up),'xdata',T(xdata),'ydata',T(ydata));axis image;
title(sprintf('coarse (%d RGCs)',length(get_cell_indices(datarun_coarse,cell_spec_coarse))))


% reconstructed image 
[average_rf_recon,xdata,ydata] = compute_average_rf(datarun_recon,cell_spec_recon);
average_rf_recon = clip_rgb_image(average_rf_recon, max_amp);
subplot(223);imagesc(matrix_scaled_up(norm_image(plot_polarity * average_rf_recon),scale_up),'xdata',T(xdata),'ydata',T(ydata));axis image;
title(sprintf('reconstructed (%d RGCs)',length(get_cell_indices(datarun_recon,cell_spec_recon))))



% compute fine profile
[x, y] = get_rf_profiles(datarun_fine,cell_spec_fine,'radius',10);
[average_x_fine, average_y_fine] = curve_from_binning(x,sum(y,2), 'bin_edges',0:.1:10);
average_y_fine = plot_polarity * average_y_fine;

% compute coarse profile
[x, y] = get_rf_profiles(datarun_coarse,cell_spec_coarse,'radius',10);
[average_x_coarse, average_y_coarse] = curve_from_binning(x,sum(y,2), 'bin_edges',0:.1:10);
average_y_coarse = plot_polarity * average_y_coarse;

% compute reconstructed profile
[x, y] = get_rf_profiles(datarun_recon,cell_spec_recon,'radius',10);
[average_x_recon, average_y_recon] = curve_from_binning(x,sum(y,2), 'bin_edges',0:.1:10);
average_y_recon = plot_polarity * average_y_recon;

% compute cone weight profile

% center
[center_x, center_y] = get_rf_cone_profiles(datarun_recon,cell_spec_recon,'radius',stixel_ratio*10,...
    'center_type', @(d,cc)d.cones.rf_fits{cc}.center,...
    'selection',struct('thresh', 0.10,'radius', [0 inf], 'polarity', 1,'contiguity', true,'scale', 3.0));
% surround
[surround_x, surround_y] = get_rf_cone_profiles(datarun_recon,cell_spec_recon,'radius',stixel_ratio*10,...
    'center_type', @(d,cc)d.cones.rf_fits{cc}.center,...
    ...'selection',struct('thresh', 0.05,'radius', [0 8], 'polarity', -1,'contiguity', false,'scale', 3.0));
    'selection',struct('thresh', 0.00,'radius', [0 8], 'polarity', 0,'contiguity', false,'scale', 3.0));
% bin
center_x = center_x/stixel_ratio;
surround_x = surround_x/stixel_ratio;
[average_x_center, average_y_center] = curve_from_binning(center_x,center_y, 'bin_edges',0:.1:8);
[average_x_surround, average_y_surround] = curve_from_binning(surround_x,surround_y, 'bin_edges',0:.1:8);
% polarity
average_y_center = plot_polarity * average_y_center;
average_y_surround = plot_polarity * average_y_surround;


% normalize
average_y_coarse = average_y_coarse / std(average_y_coarse);
average_y_fine = average_y_fine / std(average_y_fine);
average_y_recon = average_y_recon / std(average_y_recon);


% plot overlay
subplot(224);
switch 2
    case 1 % everything
        plot([flipud(-average_x_fine); average_x_fine], [flipud(average_y_fine); average_y_fine],'g',...
            [flipud(-average_x_coarse); average_x_coarse], [flipud(average_y_coarse); average_y_coarse],'r',...
            [flipud(-average_x_recon); average_x_recon], [flipud(average_y_recon); average_y_recon],'k',...
            [flipud(-average_x_center); average_x_center], [flipud(average_y_center); average_y_center],'c',...
            [flipud(-average_x_surround); average_x_surround], [flipud(average_y_surround); average_y_surround],'b')
        title('fine = green, coarse = red, reconstructed = black, center = cyan, surround = blue')
    case 2 % coarse & reconstructed
        plot([flipud(-average_x_coarse); average_x_coarse], [flipud(average_y_coarse); average_y_coarse],'r',...
            [flipud(-average_x_recon); average_x_recon], [flipud(average_y_recon); average_y_recon],'b')
        title('coarse = red, reconstructed = blue')
end

