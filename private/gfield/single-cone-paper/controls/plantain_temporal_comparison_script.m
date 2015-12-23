datarunA = load_data('2008-08-27-5', 'rf-3-plantain');
datarunA = load_index(datarunA);
datarunA = load_sta(datarunA, 'load_sta', []);
datarunA = load_params(datarunA);
datarunA = import_single_cone_data(datarunA,'2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard');
datarunA.cones.mosaic = make_mosaic_struct(datarunA.cones.centers);


datarunB = load_data('2008-08-27-5', 'rf-1-gf');
datarunB = load_index(datarunB);
datarunB = load_sta(datarunB, 'load_sta', []);
datarunB = load_params(datarunB);
datarunB = import_single_cone_data(datarunB, '2008-08-27-5_data001_data001-bayes-msf_70.00--standard');
datarunB.cones.mosaic = make_mosaic_struct(datarunB.cones.centers);

a_cone_colors = [1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 0 0 0];

plot_cone_mosaic(datarunA, 'fig_or_axes', 1, 'cone_size', 8)

plot_cone_mosaic(datarunA, 'fig_or_axes', 1, 'cone_size', 8, 'type_colors',a_cone_colors)
print(1, '~/Desktop/cone-mosaic-a.pdf','-dpdf')

plot_cone_mosaic(datarunB, 'fig_or_axes', 2, 'cone_size', 4)
print(2, '~/Desktop/cone-mosaic-b.pdf','-dpdf')



%%
% analyze the purity for the interval 1 and interval 4 1x1 stimulation in plantain

datarun = datarunB;

% set analysis parameters
binarize_weights_flag = false; % if true, weights are binarized
use_cone_files = false; % if true, cone files are loaded from server
clumped_flag = false; % if true, clumped cone mosaics are loaded, 'use_cone_files', must also be true
rf_center_threshold = 0.1; % defines threshold for what is RF center
num_permuted_cone_mosaics = 100;  % num
confidence_interval = 0.8413; % 0.8413 is 1-sigma
min_convergence = 4; % each rgc must contact at least this number of cones
verbose = 0;
cell_types = [3,4]; % analyze both the on and off midget cells

mosaic_counter = 0;

% analyze mosaic
for tp = 1:length(cell_types)
    mosaic_counter = mosaic_counter +1;
    cell_type = cell_types(tp);

    % extract connectivity
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                                'thresh', rf_center_threshold,...
                                                'radius', [0 inf], 'polarity', 1,...
                                                'contiguity', true,'scale', 3.0,...
                                                'remove_cones', 'S');   

    % determine wither weights should be binarized
    if binarize_weights_flag
        connectivity = double(selection); % binarize weights
    else
        connectivity = mosaic_weights .* selection; % keep weights continuous valued
    end

    % the new datarun has S cone excluded from all fields of datarun.cone 
    new_datarun = extras.new_datarun;

    [num_cones, num_RGCs] = size(connectivity);


    % get pointers to sampled cones
    % purpose: this ensures that on cones sampled by the mosaic are permuted
    % this ensures that L:M cone ratio of the sampled cones is the same per and post permutation
    summed_selection = sum(selection, 2);
    cone_indices = find(summed_selection > 0);
    clear cone_pointers
    cone_pointers(1:num_cones) = false;
    cone_pointers(cone_indices) = true;

    % normalize cone weights
    for RGC = 1:num_RGCs
        connectivity(:,RGC) = connectivity(:,RGC) ./ sum(connectivity(:,RGC));
    end

    % get list of pointers to rgcs with sufficiently high convergence
    keeper_list = false(num_RGCs,1);
    for cc = 1:num_RGCs
        if length(find(selection(:,cc))) >= min_convergence;
            keeper_list(cc) = true;
        end
    end

    % compute purity of cone mosaic
    purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
    purity_sd = std(purity_indices(keeper_list));
    mosaic_purities(mosaic_counter) = purity_sd;

    % determines wheter pre-calculated permutations are loaded from the server,
    % or whether permutations are calculated online.
    if use_cone_files  
        % load cone files
        for cone_mosaic = 1:num_permuted_cone_mosaics
            clear cleaned_cone_types temp_types cone_types
            suffix_string = num2str(11110+cone_mosaic);
            if clumped_flag
                str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations-no-off/clumped/',LMS_names{dataset},'-',suffix_string,'.txt'];
            else
                str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations-no-off/permuted/',LMS_names{dataset},'-',suffix_string,'.txt'];  
            end
            cones_file = dlmread([str2load],'\t');
            cone_types = cones_file(:,4);
            S_indices = find(cone_types == 3); % s cones
            U_indices = find(cone_types == 0); % untyped cones
            keep_cone_indices = setdiff([1:length(cone_types)], [S_indices; U_indices]);    

            temp_types = cone_types(keep_cone_indices);
            cleaned_cone_types(temp_types == 1,1) = 'L';
            cleaned_cone_types(temp_types == 2,1) = 'M';

            permuted_cone_PIs = compute_opponency_index(connectivity, cleaned_cone_types);
            %finite_indices = isfinite(permuted_cone_PIs);
            permuted_cone_PI_SD(cone_mosaic) = std(permuted_cone_PIs(keeper_list));
        end
    else  
        % don't load cone files
        % permute cone labels within code -- clumping not supported
        for iter = 1:num_permuted_cone_mosaics
            new_cones = new_datarun.cones.types(cone_pointers);
            new_cone_num = length(new_cones);
            rand_cone_indices = randperm(new_cone_num);
            permuted_cone_PIs = compute_opponency_index(connectivity(cone_pointers,:), new_cones(rand_cone_indices));

%            finite_indices = isfinite(permuted_cone_PIs);
            permuted_cone_PI_SD(iter) = std(permuted_cone_PIs(keeper_list));
        end
    end
    permuted_cone_purity(mosaic_counter) = mean(permuted_cone_PI_SD);
    permuted_cone_purity_error(mosaic_counter) = norminv(confidence_interval,0,1) * std(permuted_cone_PI_SD);

end

mosaic_purities
permuted_cone_purity
permuted_cone_purity_error

figure
clf
hold on
errorbar(mosaic_purities, permuted_cone_purity, permuted_cone_purity_error,'ko')
hold on
plot([0 0.5], [0 0.5], 'k')
axis([0 0.5 0 0.5])
axis square
xlabel('data')
ylabel('permute')
title('purity')
hold off
%print(1, '/snle/home/gfield/Desktop/permutation-4.pdf','-dpdf')







