

% load images
if 1
    % array map
    array_map.im = imread('/snle/home/gauthier2/Desktop/imaging/methods/alignment graphics/512.gif','gif');

    % array montage (cones AND electrodes)
    array_montage.im = imread('/snle/lab/Experiments/Array/Analysis/2006-06-12-0/images/alive-montage-10x-complete-flat.tiff','tif');

    % alignment image, electrodes
    align_32p_6x_elec.im = imread('/snle/data/2006-06-12-0/slidebook/Image 1 000.tif','tif');

    % alignment image, cones
    align_32p_6x_cones.im = imread('/snle/data/2006-06-12-0/slidebook/Image 1 001.tif','tif');

    % rf alignment matrices
    align_32p.im = imread('/snle/home/gauthier/Desktop/imaging/methods/alignment graphics/photographic mapping 32.tif','tif');
    align_10p.im = imread('/snle/home/gauthier/Desktop/imaging/methods/alignment graphics/photographic mapping 10.tif','tif');
    align_2p.im = imread('/snle/home/gauthier/Desktop/imaging/methods/alignment graphics/photographic mapping 2.tif','tif');

    % 10x array image
    %array_10x.im = imread('/snle/data/2006-06-12-0/slidebook/Image 12 000.tif','tif');

    % 10x cone image
    %cone_10x.im = imread('/snle/data/2006-06-12-0/slidebook/Image 12 001.tif','tif');

end


% make tforms
if 1
    % make tform to array map from array montage
    input_points = ...
        [1.723306480648065e+03     2.989104410441041e+02
        3.707124212421241e+02     2.299311431143105e+02
        3.227268226822681e+02     3.004098559855986e+03
        1.675320882088209e+03     3.076076957695770e+03];
    base_points = ...
        [1.501614317019722e+02     2.093323593864135e+02
        1.374083272461651e+02     5.679196493791089e+02
        8.830898466033600e+02     5.679196493791088e+02
        8.958429510591671e+02     2.085821767713661e+02];
    tf_array_map_from_array_montage = cp2tform(input_points, base_points, 'projective');


    % make tform to array_montage from align_6x_elec
    base_points =...
        [1.699517921146953e+03     1.371665770609319e+03
        1.251668458781362e+03     1.410479390681003e+03
        8.978673835125450e+02     1.134305555555555e+03
        7.067849462365592e+02     1.758309139784946e+03
        6.067652329749104e+02     2.339020609318996e+03
        1.411401433691756e+03     2.398733870967742e+03];
    input_points =...
        [7.780000000000001e+02     3.679999999999999e+02
        5.130000000000000e+02     3.900000000000000e+02
        3.020000000000000e+02     2.260000000000000e+02
        1.910000000000000e+02     5.960000000000000e+02
        1.300000000000000e+02     9.390000000000000e+02
        6.080000000000000e+02     9.740000000000000e+02];
    tf_array_montage_from_align_32p_6x_elec = cp2tform(input_points, base_points, 'projective');
    
    
    % make tform to align_32p_6x_cones from align_32p
    input_points =...
        [3.522500000000000e+02     2.567500000000000e+02
        3.202500000000000e+02     2.242500000000000e+02
        2.562500000000000e+02     2.242500000000000e+02
        2.247500000000000e+02     2.567500000000000e+02
        2.567500000000000e+02     2.887500000000000e+02
        2.567500000000000e+02     3.205000000000000e+02
        3.527500000000000e+02     3.205000000000000e+02
        3.842500000000000e+02     2.882500000000000e+02
        3.842500000000000e+02     1.927500000000000e+02
        3.842500000000000e+02     1.282500000000000e+02
        3.522500000000000e+02     1.602500000000000e+02
        2.887500000000000e+02     9.675000000000000e+01
        2.242500000000000e+02     1.607500000000000e+02];
    base_points =...
        [8.460000000000000e+02     3.170000000000000e+02
        6.820000000000000e+02     4.889999999999999e+02
        6.790000000000000e+02     8.320000000000000e+02
        8.410000000000000e+02     1.005000000000000e+03
        1.006000000000000e+03     8.320000000000000e+02
        1.169000000000000e+03     8.330000000000000e+02
        1.173000000000000e+03     3.150000000000000e+02
        1.010000000000000e+03     1.430000000000000e+02
        5.240000000000000e+02     1.420000000000000e+02
        2.109999999999999e+02     1.380000000000000e+02
        3.639999999999999e+02     3.140000000000000e+02
        4.799999999999989e+01     6.550000000000000e+02
        3.629999999999999e+02     1.004000000000000e+03];
    tf_align_32p_6x_cones_from_align_32p = cp2tform(input_points, base_points, 'projective');


    % make tform to array_montage from align_32p
    input_points =...
        [3.522500000000000e+02     2.567500000000000e+02
        3.202500000000000e+02     2.242500000000000e+02
        2.562500000000000e+02     2.242500000000000e+02
        2.247500000000000e+02     2.567500000000000e+02
        2.567500000000000e+02     2.887500000000000e+02
        2.567500000000000e+02     3.205000000000000e+02
        3.527500000000000e+02     3.205000000000000e+02
        3.842500000000000e+02     2.882500000000000e+02
        3.842500000000000e+02     1.927500000000000e+02
        3.842500000000000e+02     1.282500000000000e+02
        3.522500000000000e+02     1.602500000000000e+02
        2.887500000000000e+02     9.675000000000000e+01
        2.242500000000000e+02     1.607500000000000e+02];
    base_points = tformfwd(tf_align_32p_6x_cones_from_align_32p,input_points);
    base_points = tformfwd(tf_array_montage_from_align_32p_6x_elec,base_points);
    tf_array_montage_from_align_32p = cp2tform(input_points, base_points, 'projective');

end


% compare
if 1
    %[spatial_sensitivity,all_sig_stixels,spatial_cell_ids] =...
    %    compute_spatial_sensitivity(datarun, [31 33 34 197 256 273 274 468 602 634], spat_sens_params);
    
    % plot single RF
    if 1
        [rf_tf,xdata,ydata] = imtransform(rf, tf_array_montage_from_align_32p);
        vdata = [1 size(rf,1)];
        udata = [1 size(rf,2)];
        [temp,xdata,ydata] = imtransform(matrix_scaled_up(rf,4), tf_array_montage_from_align_32p,...
            'udata',udata,'vdata',vdata,'XYScale',2);
    end
    
    figure(4);clf;colormap gray
    subplot(121); imagesc(xdata,ydata,norm_image(temp)); axis image
    subplot(122); imagesc(array_montage.im(:,:,1)); axis image
    linkaxes([subplot(121) subplot(122)])
    
    % overlay cones
    if 1
        plot_centers = tformfwd(tf_array_montage_from_align_32p,datarun.cones.centers);
        subplot(122); hold on
        plot(plot_centers(:,1),plot_centers(:,2),'.r','MarkerSize',10)
    end
end

