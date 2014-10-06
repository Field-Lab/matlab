% s-cone efficiency script
% computes fraction of total input to parasol and midget cells is provided by S cones.
%% new cone finding

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';

% blueberry
path_and_name{2,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{2,2} = 'blueberry';

% apricot
path_and_name{3,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{3,2} = 'apricot';

% peach
path_and_name{4,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{4,2} = 'peach';


%kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi';

% grapes
path_and_name{1,1} = '/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = 'grapes';

% apple
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
path_and_name{1,2} = 'apple';

%% old cone finding

% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain-old';

% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = 'blueberry-old';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{1,2} = 'apricot-old';

%kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi-old';

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = 'peach-old';

% grapes
path_and_name{1,1} = '/jacob/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = 'grapes-old';

% apple
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data013/data013';
path_and_name{1,2} = 'apple-old';

%%

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


% get S cone centers
[s_weights, s_selection, s_extras] = select_cone_weights(datarun, {5,6}, 'thresh', 3,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false, 'scale', 20.0, 'remove_cones', 'LMU');   
s_cone_centers = s_extras.new_datarun.cones.centers;

% or just grab all S cones
s_cone_indices = find(datarun.cones.types == 'S');
% s_cone_centers = datarun.cones.centers(s_cone_indices,:);

cell_types = [1,2,3,4];
mosaic_s_cone_encounters = zeros(length(cell_types), 1);
mosaic_s_cone_hits = zeros(length(cell_types), 1);
mosaic_s_cone_efficiency = zeros(length(cell_types), 1);

thresholds = 3;
radius = 3;

s_cone_fraction = zeros(4,1);
s_cone_efficiency = zeros(4,1);
plot_counter = 0;
figure
for thresh = 1:length(thresholds)

    weight_fraction = thresholds(thresh);

    for ctype = 1:length(cell_types)
        
        cell_indices = get_cell_indices(datarun, {ctype});
        num_rgcs = length(cell_indices);
        %cell_radii = zeros(num_rgcs, 1);

        [all_weights, all_selection, extras] = select_cone_weights(datarun, {ctype}, 'thresh', 0,...
                                                'radius', [0 radius], 'polarity', 0,...
                                                'contiguity', false, 'remove_cones', 'U');     
   
        [temp_s_weights, temp_s_selection, extras] = select_cone_weights(datarun, {ctype}, 'thresh', 0,...
                                                'radius', [0 radius], 'polarity', 0,...
                                                'contiguity', false, 'remove_cones', 'LMU');    

        all_connectivity = all_weights .* all_selection;
        temp_s_connectivity = temp_s_weights .* temp_s_selection;

        % get_number of cones for each cell of given type
        all_cone_counter = sum(all_selection,1); % total number of cones in each cell
        s_cone_counter = sum(temp_s_selection,1); % total s cones in each cell

        s_cone_present_fraction = s_cone_counter ./ all_cone_counter;


        plot_counter = plot_counter +1;
        subplot(4,2,plot_counter)
        [weight_vals, hist_bins] = hist(reshape(all_weights,[],1),50);
        plot(hist_bins,weight_vals)

        plot_counter = plot_counter +1;
        subplot(4,2,plot_counter)
        [s_weight_vals, hist_bins] = hist(reshape(temp_s_weights,[],1),hist_bins);
        plot(hist_bins, s_weight_vals)

        length(find(temp_s_connectivity))

        % NOW GET CONES WITH WEIGHTS
        [all_weights, all_selection, extras] = select_cone_weights(datarun, {ctype}, 'thresh', thresholds,...
                                            'radius', [0 radius], 'polarity', 0,...
                                            'contiguity', false, 'remove_cones', 'U');     
   
        [temp_s_weights, temp_s_selection, extras] = select_cone_weights(datarun, {ctype}, 'thresh', thresholds,...
                                                'radius', [0 radius], 'polarity', 0,...
                                                'contiguity', false, 'remove_cones', 'LMU');    
       

        lm_cone_hits = sum(all_selection);
        s_cones_hit = sum(temp_s_selection);
        
        efficiency = s_cones_hit ./ lm_cone_hits;


        % remove NaNs
        s_cone_present_fraction(isnan(s_cone_present_fraction)) = 0;
        efficiency(isnan(efficiency)) = 0;

        mean_fraction= mean(s_cone_present_fraction);
        mean_efficiency = mean(efficiency);

        s_cone_fraction(ctype) = mean_fraction;
        s_cone_efficiency(ctype) = mean_efficiency;

    end

end

s_cone_fraction
s_cone_efficiency

%% compute total fraction of S cone input to each cell type

num_dset = size(path_and_name,1);
s_input = cell(4,1);

for dset = 1:num_dset

    % load data
    datarun = load_data(path_and_name{dset,1});
    datarun = load_params(datarun,struct('verbose',1));  
    datarun = load_sta(datarun,'load_sta',[]);
    datarun = set_polarities(datarun);
    datarun = import_single_cone_data(datarun, path_and_name{dset,2});    
    datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

    cell_types = [1,2,3,4];
    mosaic_s_cone_encounters = zeros(length(cell_types), 1);
    mosaic_s_cone_hits = zeros(length(cell_types), 1);
    mosaic_s_cone_efficiency = zeros(length(cell_types), 1);

    %thresholds = 3;
    radius = 1.5;
    s_density_threshold = 0.06;

    fraction_s_input = zeros(4,1);
    figure
    for ctype = 1:length(cell_types)

        cell_indices = get_cell_indices(datarun, {ctype});
        num_rgcs = length(cell_indices);
        %cell_radii = zeros(num_rgcs, 1);

        [all_weights, all_selection, extras] = select_cone_weights(datarun, {ctype}, 'thresh', 0,...
                                                'radius', [0 radius], 'polarity', 1,...
                                                'contiguity', false);     

        [temp_s_weights, temp_s_selection, extras] = select_cone_weights(datarun, {ctype}, 'thresh', 0,...
                                                'radius', [0 radius], 'polarity', 0,...
                                                'contiguity', false, 'remove_cones', 'LMU');    

        % only keep cellls that exceed a threshold S cone density in their region
        all_cone_counter = sum(all_selection,1); % total number of cones in each cell
        s_cone_counter = sum(temp_s_selection,1); % total s cones in each cell

        s_cone_present_fraction = s_cone_counter ./ all_cone_counter;

        keeper_cells = find(s_cone_present_fraction >= s_density_threshold);


        all_connectivity = all_weights .* all_selection;
        s_connectivity = temp_s_weights .* temp_s_selection;

        total_cone_input = sum(all_connectivity,1);
        total_s_input = sum(s_connectivity,1);

        temp_fraction_s_input = total_s_input(keeper_cells) ./ total_cone_input(keeper_cells);

        % remove NaNs
        temp_fraction_s_input(isnan(temp_fraction_s_input)) = 0;

        subplot(2,2,ctype)
        hist_bins = [-0.02:0.001:0.2];
        [temp_vals, hist_bins] = hist(temp_fraction_s_input,hist_bins);
        bar(hist_bins, temp_vals)
        title([datarun.cell_types{ctype}.name, ' (',num2str(length(keeper_cells)),')'])
        axis([-0.05 0.2 0 10])

        fraction_s_input(ctype) = mean(temp_fraction_s_input);

        temp_input = s_input{ctype};
        s_input{ctype} = [temp_input, temp_fraction_s_input];

    end
fraction_s_input
end

%%
figure
hist_bins = [-0.02:0.001:0.2];
mean_s_input = zeros(4,1);
std_s_input = zeros(4,1);
skew_s_input = zeros(4,1);
for ct = 1:4
    subplot(2,2,ct)
    [temp_val, hist_bins] = hist(s_input{ct},hist_bins);
    bar(hist_bins, temp_val)
    title([datarun.cell_types{ctype}.name, ' (',num2str(length(s_input{ct})),')'])
    axis([-0.05 0.2 0 20])
    mean_s_input(ct) = mean(s_input{ct});
    std_s_input(ct) = std(s_input{ct});
    skew_s_input(ct) = skewness(s_input{ct});
end

mean_s_input
std_s_input
skew_s_input








