%% peach

%% BW data
BW_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data003/data003';

BW_path = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data012/data012';


% load BW data
BW_datarun = load_data(BW_path);
BW_datarun = load_params(BW_datarun,struct('verbose',1));  
BW_datarun = load_sta(BW_datarun,'load_sta',[]);
BW_datarun = set_polarities(BW_datarun);


%% RGB Data
% peach
RGB_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
RGB_cone_path = 'peach';

RGB_path = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
RGB_cone_path = 'apple';



% load RGB data
RGB_datarun = load_data(RGB_path);
RGB_datarun = load_params(RGB_datarun,struct('verbose',1));  
RGB_datarun = load_sta(RGB_datarun,'load_sta',[]);
RGB_datarun = set_polarities(RGB_datarun);
RGB_datarun = import_single_cone_data(RGB_datarun, RGB_cone_path);    
RGB_datarun.cones.mosaic = make_mosaic_struct(RGB_datarun.cones.centers);

% load cone RFs
cd /snle/lab/Experiments/Array/Shared/one/apple
load Wc

% Convert Wc into "colorless" cone kernels 
BW_cone_rfs = zeros(size(Wc,1)./3,size(Wc,2));
for cn = 1:size(Wc,2)
    tmp_cn = Wc(:,cn);
    tmp_cn = reshape(full(tmp_cn),[320,320,3]);
    summed_cn_colors = squeeze(sum(tmp_cn,3));
    reshaped_cn = reshape(summed_cn_colors,[],1);
    reshaped_cn = reshaped_cn ./ norm(reshaped_cn);
    BW_cone_rfs(:,cn) = reshaped_cn;
end
    
BW_cone_rfs = sparse(BW_cone_rfs);

%% Fit the BW STAs with the cones identified from RGB Data

% GET CONE WEIGHTS IN EACH RGC

num_cones = size(Wc,2);

% get list of cell indices
cell_indices = get_cell_indices(BW_datarun,'all');

% note when it started
fprintf('\nComputing cone weights in %d RFs',length(cell_indices));
start_time_regress = clock; 

% initialize
cone_weights = zeros(num_cones,length(cell_indices));

% go through list of cells
for cc = 1:length(cell_indices)
    
    fprintf('.')

    % get summary frame
    rf = get_rf(BW_datarun,BW_datarun.cell_ids(cell_indices(cc)));

    if isempty(rf)
        exit_flag = 1
        continue
    end

    % reshape for the regression
    rf = reshape(rf,[],1);

    % put in units of SNR
    rf = rf / robust_std(rf);


    % regress to get the cone weights
    cone_weights(:,cc) = BW_cone_rfs\rf;

end

% display how long it took
fprintf('\n   done (%0.1f seconds)\n',etime(clock,start_time_regress));


%% use RGB cones to fit BW STAs

BW_datarun.cones = RGB_datarun.cones;
BW_datarun.cones.weights = cone_weights;
BW_datarun.cones = rmfield(BW_datarun.cones, 'rf_fits');


cone_info.cone_weights = cone_weights;
cone_info.cone_centers = RGB_datarun.cones.centers;
cone_info.cell_ids = BW_datarun.cell_ids;

BW_datarun = get_sta_summaries(BW_datarun, {3,4}, 'keep_stas', false)

BW_datarun = fit_cone_rfs(BW_datarun,{3,4},'foa_profile',[],'fit_radius',150,'verbose',1,'foa_2d',4);


%% use BW cones to fit RGB STAs

% RGB_datarun.cones.center = BW_datarun.cones;
% BW_datarun.cones.weights = cone_weights;
% BW_datarun.cones = rmfield(BW_datarun.cones, 'rf_fits');
% 
% 
% cone_info.cone_weights = cone_weights;
% cone_info.cone_centers = RGB_datarun.cones.centers;
% cone_info.cell_ids = BW_datarun.cell_ids;
% 
% BW_datarun = get_sta_summaries(BW_datarun, {3,4}, 'keep_stas', false)
% 
% BW_datarun = fit_cone_rfs(BW_datarun,{3,4},'foa_profile',[],'fit_radius',150,'verbose',1,'foa_2d',4);
% 

%%

use_cone_files = true;
clumped_flag = true;


%%
temp_cell_type = 3;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(BW_datarun, {temp_cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'S');   

connectivity = mosaic_weights .* selection; % keep weights continuous valued


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

% compute purity of cone mosaic
purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
finited_indices = isfinite(purity_indices);
purity_sd = std(purity_indices(finited_indices));

% determines wheter pre-calculated permutations are loaded from the server,
% or whether permutations are calculated online.
if use_cone_files  
    % load cone files
    for cone_mosaic = 1:100
        clear cleaned_cone_types temp_types cone_types
        suffix_string = num2str(11110+cone_mosaic);
        if clumped_flag
            str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/clumped/',RGB_cone_path,'-',suffix_string,'.txt'];
        else
            str2load = ['/snle/lab/Experiments/Array/Shared/one/simulations/permuted/',RGB_cone_path,'-',suffix_string,'.txt'];  
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
        finite_indices = isfinite(permuted_cone_PIs);
        permuted_cone_PI_SD(cone_mosaic) = std(permuted_cone_PIs(finite_indices));
    end
else  
    % don't load cone files
    % permute cone labels within code -- clumping not supported
    for iter = 1:num_permuted_cone_mosaics
        new_cones = new_datarun.cones.types(cone_pointers);
        new_cone_num = length(new_cones);
        rand_cone_indices = randperm(new_cone_num);
        permuted_cone_PIs = compute_opponency_index(connectivity(cone_pointers,:), new_cones(rand_cone_indices));

        finite_indices = isfinite(permuted_cone_PIs);
        permuted_cone_PI_SD(iter) = std(permuted_cone_PIs(finite_indices));
    end
end

purity_sd
mean(permuted_cone_PI_SD)
std(permuted_cone_PI_SD)
hold on
errorbar(purity_sd, mean(permuted_cone_PI_SD), std(permuted_cone_PI_SD), 'go')




