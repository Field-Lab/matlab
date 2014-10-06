[LMS_paths, LMS_names] = get_LMS_paths('high');

num_datasets = length(LMS_paths);
verbose = 0;

mosaic_counter = 0;

for dataset = 1:num_datasets
    clear datarun new_datarun sim_datarun
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = import_single_cone_data(datarun, LMS_names{dataset});
    
    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);

        if num_on_midgets > 20
        on_midget_flag = true;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > 20
        off_midget_flag = true;
    else
        off_midget_flag = false;
    end
    
    if on_midget_flag && off_midget_flag
        cell_types = [3,4];
    elseif on_midget_flag && ~off_midget_flag
        cell_types = 3;
    elseif ~on_midget_flag && off_midget_flag
        cell_type = 4;
    else
        cell_types = [];
    end
    
    for tp = 1:length(cell_types)
        mosaic_counter = mosaic_counter + 1;
        cell_type = cell_types(tp);


        temp_cell_indices = get_cell_indices(datarun, {cell_type});

        ss_params.thresh = 3.0;
        distance_threshold = 20;
        temp_cell_indices = get_cell_indices(datarun, {cell_type});
        for RGC = 1:length(temp_cell_indices);
            RGC
            temp_sta = get_sta(datarun, datarun.cell_ids(temp_cell_indices(RGC)));
            new_datarun.stas.stas{RGC} = temp_sta;
            datarun.stas.rf_coms{temp_cell_indices(RGC)} = rf_com(temp_sta,'ss_params',ss_params, 'distance_thresh', distance_threshold, 'positive', true);
            datarun.stas;
        end




       % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.05, 'radius', [0 inf], 'polarity', 1,...
                                                    'contiguity', true, 'scale', 3.0, 'remove_cones', 'SU');   
        connectivity = mosaic_weights .* selection;
        new_datarun = extras.new_datarun;
        [num_cones, num_RGCs] = size(connectivity);   
        
        for RGC = 1:num_RGCs
            temp_cone_indices = find(connectivity(:,RGC) > 0);
            temp_cone_weights = connectivity(temp_cone_indices, RGC);
            temp_cone_locations = new_datarun.cones.centers(temp_cone_indices, :);

            % get peak cone location
            [max_weight, max_index] = max(temp_cone_weights);
            max_index = temp_cone_indices(max_index);
            temp_fit = datarun.cones.rf_fits{temp_cell_indices(RGC)};

            
            % set center
            center = datarun.stas.rf_coms{temp_cell_indices(RGC)};
            %center = temp_fit.center;
            %center = new_datarun.cones.centers(max_index,:);

            params(1) = temp_fit.center_radius; %x_sd
            params(2) = temp_fit.center_radius; %y_sd
            if params(1) < 2  
               params(1) = 5;
               params(2) = 5;
            end

            params(3) = temp_fit.center_scale; % amplitude scale
            params(4) = 0; % angle

            [fit_params, fval, exitflag, options] = fminsearch(@(params) gauss_2d(params, center, temp_cone_weights, temp_cone_locations), params, struct('MaxFunEvals', 2000));
%             fval
%             options.funcCount
            fit_params

        % plot reconstructed RF
            % get reconstruction
            cone_rf = cone_rf_reconstructed([datarun.stimulus.field_height datarun.stimulus.field_width],...
                temp_cone_weights,temp_cone_locations);
            % plot it
            figure(1);clf;
            image(norm_image(cone_rf)); axis image; drawnow
            % add center point
            hold on; plot(center(1),center(2),'r.')


       % plot fit
            % generate points to plot a circle at one center sigma

            %sx = fit_params(1);
            %sy = fit_params(2);
            %w = fit_params(4);
            
            mx = center(1);
            my = center(2);
            sx = temp_fit.center_radius;
            sy = temp_fit.center_radius;
            w = 0;


            % circle samples
            circle_samples = 0:0.05:2*pi;
            x_circle = cos(circle_samples);
            y_circle = sin(circle_samples);

            % rotate by angle and stretch
            R = rmatrix2(-w / (2*pi) * 360);
            L = [sx 0; 0 sy];
            X = R * L * [x_circle; y_circle];
            X(:,end+1) = X(:,1);

            X(1,:)=(X(1,:)+mx);
            X(2,:)=(abs((X(2,:)+my)));

            plot(X(1,:), X(2,:), 'r');
            window = 30;
            axis([center(1)-window center(1)+window center(2)-window center(2)+window])
            pause
        end
   
 
    end
end



