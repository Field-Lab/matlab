% the 1-sigma fits over the top of this
cell_type = 4;
cell_indices = get_cell_indices(datarun, {cell_type});
num_rgcs = length(cell_indices);
%datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);
 
% get connectivity and purities
% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   
     
center_connectivity = mosaic_weights .* selection;                                       
% get elliptical fits from cone weights
 

[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.0,...
                                            'radius', [0 4], 'polarity', 0,...
                                            'contiguity', false,'scale', 3.0,...
                                            'remove_cones', 'U');   
full_connectivity = mosaic_weights .* selection;                                       
% get elliptical fits from cone weights
 
fit_params.fit_center = true;
fit_params.fit_center_sd = true;
fit_params.fit_center_scale = true;
fit_params.fit_surround_scale = false;
fit_params.fit_surround_sd_scale = true;
 
which_to_fit = [...
    fit_params.fit_center...
    fit_params.fit_center_sd...
    fit_params.fit_center_scale...
    fit_params.fit_surround_scale...
    fit_params.fit_surround_sd_scale...
    ];
 
fit_indices = find(which_to_fit);
fix_indices = find(~which_to_fit);
 
 
for rgc = 1:num_rgcs
    center_cone_indices = find(center_connectivity(:,rgc));
    center_cone_locations = datarun.cones.centers(center_cone_indices,:);
    % get distances between cones
    cone_dists = ipdm(center_cone_locations);
    max_dist = max(max(cone_dists));
    dist_index = find(cone_dists == max_dist,1);
    [cone_A, cone_B] = ind2sub(size(cone_dists), dist_index);
    cone_A_location = center_cone_locations(cone_A, :);
    cone_B_location = center_cone_locations(cone_B, :);
    axis_vector = cone_A_location - cone_B_location;
    [init_theta, junk_rho] = cart2pol(axis_vector(1), axis_vector(2));
    
    cone_indices = find(full_connectivity(:,rgc));
    cone_locations = datarun.cones.centers(cone_indices,:);
    cone_weights = full_connectivity(cone_indices,rgc);
    cone_weights = cone_weights ./ max(cone_weights);

    
    circular_fit = datarun.cones.rf_fits{cell_indices(rgc)};
    
    clear x;
    % set starting conditions of fit
    %x(1) = circular_fit.center(1);
    %x(2) = circular_fit.center(2);
    x(1) = 1;
    x(2) = circular_fit.center_radius;
    x(3) = circular_fit.center_radius;
    x(4) = init_theta;
    surround_sd_scale = 2;
    surround_scale = 0;
    center_point = circular_fit.center;
    
    % fminsearch
%     fit_fun = @(x)spatial_diff_of_Gaussians_err(x,surround_sd_scale, surround_scale, cone_locations,cone_weights);
%     fit_options = optimset('MaxIter', 10000, 'MaxFunEval', 4000);
%     fit_params = fminsearch(fit_fun, x,fit_options)
%     
 
    % fmincon
    fit_fun = @(x)spatial_diff_of_Gaussians_err(x,center_point,surround_sd_scale, surround_scale, cone_locations,cone_weights);
    fit_options = optimset('MaxIter', 10000, 'MaxFunEval', 4000);
    lb = [0 5 5 -pi/2];
    ub = [inf 1000 1000 pi/2];
    fit_options = optimset('Algorithm', 'interior-point', 'MaxIter', 10000, 'MaxFunEval', 4000);
    [fit_params, fval] = fmincon(fit_fun, x,[],[],[],[],lb,ub,[],fit_options)
 
    
    % plot the fit and the weights to check
    figure(1); clf;
    sig_cones = Wc(:,cone_indices);
    sig_cones_image = sig_cones * cone_weights;
    sig_cones_image = reshape(sig_cones_image, [320, 640, 3]);
    image(norm_image(sig_cones_image))
    axis image; hold on
    drawEllipse(center_point(1),center_point(2),sqrt(fit_params(2)),sqrt(fit_params(3)),fit_params(4));
    window_size = 30;
    begin_x = circular_fit.center(1) - window_size;
    end_x = circular_fit.center(1) + window_size;
    begin_y = circular_fit.center(2) - window_size;
    end_y = circular_fit.center(2) + window_size;
    axis([begin_x end_x begin_y end_y])
    plot([cone_A_location(1), cone_B_location(1)], [cone_A_location(2), cone_B_location(2)], 'k')
    pause

    
    
end
    
    
    % lsqcurvefit
%     fit_fun = @(x, surround_sd_scale, surround_scale, cone_locations)spatial_diff_of_Gaussians(x,surround_sd_scale, surround_scale, cone_locations);
%     lb = [0 0 0 1 1 0];
%     ub = [inf inf inf 1000 1000 inf];
%     fit_params = lsqcurvefit(fit_fun, x, cone_locations, cone_weights, lb, ub);
%     
 
    
    
                                        
                                        
