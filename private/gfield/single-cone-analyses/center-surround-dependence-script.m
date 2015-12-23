
% script for comparing the surround strength with the center strength and correlating this with purity

% get data names and paths
[LMS_paths, LMS_names, ] = get_LMS_paths('high');

num_datasets = length(LMS_paths);
verbose = 0;
mosaic_counter = 0;

for dataset = 1:num_datasets
    clear datarun 
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

%         ss_params.thresh = 2.5;
%         distance_threshold = 40;
%         temp_cell_indices = get_cell_indices(datarun, {cell_type});
%         for RGC = 1:length(temp_cell_indices);
%             RGC
%             temp_sta = get_sta(datarun, datarun.cell_ids(temp_cell_indices(RGC)));
%             new_datarun.stas.stas{RGC} = temp_sta;
%             datarun.stas.rf_coms{temp_cell_indices(RGC)} = rf_com(temp_sta,'ss_params',ss_params, 'distance_thresh', distance_threshold, 'positive', true);
%             datarun.stas;
%         end
%         % new fits
%         datarun = fit_cone_rfs(datarun, {cell_type},'verbose', true, 'show_fit_params', false, 'foa_profile', 1, 'foa_2d', 2);
%         temp_filename = [LMS_names{dataset},'-',num2str(tp)];
%         save(temp_filename, 'datarun')

          temp_filename = [LMS_names{dataset},'-',num2str(tp)];
          load(temp_filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract center connectivity
        [center_weights, center_selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 2.5,...
                                                    'radius', [0 2.5], 'polarity', 1,...
                                                    'contiguity', false, 'scale', 3.0, 'remove_cones', 'SU');   

        center_connectivity = center_weights .* center_selection;
        new_datarun = extras.new_datarun;
        [num_cones, num_RGCs] = size(center_connectivity);                                     

        % compute purity of cone mosaic
        purity_indices = compute_opponency_index(center_connectivity, new_datarun.cones.types);
        purity_sd = std(purity_indices);
        mosaic_purities(mosaic_counter) = purity_sd;

        total_center_vols = sum(center_connectivity, 1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract surround connectivity
        [surround_weights, surround_selection, surround_extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 2.5,...
                                                    'radius', [0 4], 'polarity', -1,...
                                                    'contiguity', false, 'scale', 5.0, 'remove_cones', 'SU');   

        surround_connectivity = surround_weights .* surround_selection;
        [num_cones, num_RGCs] = size(surround_weights);            
                         
        % compute purity of cone mosaic
        surround_purity_indices = compute_opponency_index(surround_connectivity, new_datarun.cones.types);
        
        % remove NaNs
        nan_indices = isnan(surround_purity_indices);
        non_nans = find(nan_indices == 0);
        surround_purity_indices = surround_purity_indices(non_nans);
        surround_purity_sd = std(surround_purity_indices);
        surround_mosaic_purities(mosaic_counter) = surround_purity_sd;

        temp_nans = length(find(nan_indices == 1));
        num_nans(mosaic_counter) = temp_nans

        % surround convergences
        temp_converge = zeros(1, num_RGCs);
        for rgc = 1:num_RGCs
            temp_converg(rgc) = length(find(abs(surround_connectivity(:,rgc)) > 0));
        end
        surround_convergs(mosaic_counter) = mean(temp_converg);
        %hist(temp_converg, [0:5:150]) 


        % remove same cell from center pool
        purity_indices = purity_indices(non_nans);
        purity_sd = std(purity_indices);
        mosaic_purities(mosaic_counter) = purity_sd;       

        total_surround_vols = sum(surround_connectivity, 1);
        total_surround_vols = total_surround_vols(non_nans);
        total_center_vols = total_center_vols(non_nans);

%         figure(1)
%         plot(purity_indices, surround_purity_indices, 'ko')
%         axis([-1.5 1.5 -1.5 1.5])
%         axis square
%         xlabel('center purity')
%         ylabel('surround_purity')
         rel_surround = -1.*total_surround_vols./total_center_vols;
         relative_surround_strength(mosaic_counter) = mean(rel_surround)
% 
%         figure(2)
% %        hist(rel_surround, [0:0.05:1.5])
%         plot(rel_surround, abs(purity_indices), 'ko')
%         xlabel('surround_strength')
%         ylabel('center_purity')


%         % compute surround strength based of profile fit
%         temp_indices = get_cell_indices(datarun, {cell_type});
%         center_area = zeros(1, length(temp_indices));
%         surround_area = zeros(1, length(temp_indices));
%         for rgc = 1:length(temp_indices)
%             temp_fit = datarun.cones.rf_fits{temp_indices(rgc)};
%             center_area(rgc) = temp_fit.center_radius^2 .* temp_fit.center_scale * pi;
%             surround_area(rgc) = temp_fit.surround_radius^2 .* temp_fit.surround_scale * pi;
%         end
%         surround_strength = surround_area ./ center_area;
% 
%         fit_surround_strength(mosaic_counter) = surround_strength;
%         figure(3)
%         plot(surround_strength, abs(purity_indices), 'ko')

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % permute cone labels and get distribution of PIs
        num_iters = 100;
        for iter = 1:num_iters
            rand_cone_indices = randperm(length(new_datarun.cones.types));
            permuted_cone_PIs = compute_opponency_index(center_connectivity, new_datarun.cones.types(rand_cone_indices));
            permuted_cone_PIs = permuted_cone_PIs(non_nans);
            permuted_cone_PI_SD(iter) = std(permuted_cone_PIs);
        end
        permuted_cone_purity(mosaic_counter) = mean(permuted_cone_PI_SD);
        permuted_cone_purity_error(mosaic_counter) = std(permuted_cone_PI_SD);


    end

end

figure(5)
relative_purity = mosaic_purities./permuted_cone_purity;
plot(relative_surround_strength, relative_purity, 'k.')


datanames = cell(1, length(mosaic_purities));
for ii = 1:length(mosaic_purities)
    if mod(ii,2)
        suffix = '1';
    else
        suffix = '2';
    end
    datanames{ii} = [LMS_names{round(ii/2)},'-',suffix]
end


figure(2)
clf
hold on
for ii = 1:length(datanames)
    text(relative_surround_strength(ii), relative_purity(ii), datanames{ii})
end
axis([0 0.4 0.5 1.7])








        