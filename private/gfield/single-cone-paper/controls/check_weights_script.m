%%%% check for a bias in weights to see if higher L or M cone weights can explain the purity effect
%%%% The Lennie hypothesis


[LMS_paths, LMS_names, cone_paths] = get_LMS_paths('high');

num_datasets = length(LMS_paths);
verbose = 0;

mosaic_counter = 0;
cell_cutoff = 5;
num_cone_mosaics = 64;
use_cone_files = false;
clumped_flag = false;

for dataset = 1:num_datasets
    clear datarun new_datarun sim_datarun cell_types
    datarun = load_data(LMS_paths{dataset});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    %datarun = import_single_cone_data(datarun, LMS_names{dataset});
    datarun = import_single_cone_data(datarun, cone_paths{dataset});    

    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

    num_on_midgets = length(datarun.cell_types{3}.cell_ids);
    num_off_midgets = length(datarun.cell_types{4}.cell_ids);

    if num_on_midgets > cell_cutoff
        on_midget_flag = true;
    else
        on_midget_flag = false;
    end
    if num_off_midgets > cell_cutoff
        off_midget_flag = true;
    else
        off_midget_flag = false;
    end
    
    if on_midget_flag && off_midget_flag
        cell_types = [3,4];
    elseif on_midget_flag && ~off_midget_flag
        cell_types = 3;
    elseif ~on_midget_flag && off_midget_flag
        cell_types = 4;
    else
        cell_types = [];
    end
    
    
    for tp = 1:length(cell_types)
        mosaic_counter = mosaic_counter + 1;
        cell_type = cell_types(tp);
        % extract connectivity
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type}, 'thresh', 0.1, 'radius', [0 inf], 'polarity', 1,...
                                                    'contiguity', true, 'scale', 3.0, 'remove_cones', 'SU');   
                                     


        connectivity = mosaic_weights .* selection;
        new_datarun = extras.new_datarun;
        [num_cones, num_RGCs] = size(connectivity);

        L_indices = find(new_datarun.cones.types == 'L');
        M_indices = find(new_datarun.cones.types == 'M');
 

        L_weights = [];
        M_weights = [];
        for RGC = 1:num_RGCs
            % normalize cone weight
            connectivity(:,RGC) = connectivity(:,RGC) ./ sum(connectivity(:,RGC));
            temp_cone_indices = find(connectivity(:,RGC) > 0);
            
            temp_L_indices = intersect(L_indices, temp_cone_indices);
            temp_M_indices = intersect(M_indices, temp_cone_indices);
                
            L_weights = [L_weights; connectivity(temp_L_indices, RGC)];
            M_weights = [M_weights; connectivity(temp_M_indices, RGC)];

        end

        L_mean(mosaic_counter) = mean(L_weights);
        M_mean(mosaic_counter) = mean(M_weights);

        % interate cone ids and recomp the  means
        num_iter = 100;
        rand('twister', 11111);
        temp_L_means = zeros(1,num_iter);
        temp_M_means = temp_L_means;
        for iter = 1:num_iter

            rand_cone_indices = randperm(num_cones);
            new_datarun.cones.types = new_datarun.cones.types(rand_cone_indices);

            L_indices = find(new_datarun.cones.types == 'L');
            M_indices = find(new_datarun.cones.types == 'M');

            L_weights = [];
            M_weights = [];
            for RGC = 1:num_RGCs
                temp_cone_indices = find(connectivity(:,RGC) > 0);

                temp_L_indices = intersect(L_indices, temp_cone_indices);
                temp_M_indices = intersect(M_indices, temp_cone_indices);

                L_weights = [L_weights; connectivity(temp_L_indices, RGC)];
                M_weights = [M_weights; connectivity(temp_M_indices, RGC)];
            end

            temp_L_means(iter) = mean(L_weights);
            temp_M_means(iter) = mean(M_weights);
    
         end
   
        perm_L_mean(mosaic_counter) = mean(temp_L_means);
        perm_M_mean(mosaic_counter) = mean(temp_M_means);
        perm_L_sd(mosaic_counter) = std(perm_L_mean);
        perm_M_sd(mosaic_counter) = std(perm_M_mean);
        

%         figure(1)
%         bins = 0:0.01:0.4;
%         [L_hist, hist_bins] = hist(L_weights, bins);
%         [M_hist, hist_bins] = hist(M_weights, bins);
%         subplot(2,1,1)
%         xlabel('L')
%         bar(hist_bins, L_hist)
%         subplot(2,1,2)
%         bar(hist_bins, M_hist)
%         prep = [LMS_names{dataset},'-',num2str(cell_types(tp))];
%         title(prep)
%         xlabel('M')
        %pause
        %temp_path = ['~/Desktop/bias/',prep];
        %print(1, temp_path, '-dpdf')

    end
end


figure(1)
clf; hold on
errorbar(L_mean, perm_L_mean, perm_L_sd, 'ko')
plot([0 0.1], [0 0.1], 'k')
xlabel('data')
ylabel('permuted cones')
title('l cones')

figure(2)
clf; hold on
errorbar(M_mean, perm_M_mean, perm_M_sd, 'ko')
plot([0 0.1], [0 0.1], 'k')
xlabel('data')
ylabel('permuted cones')
title('m cones')




