function [RGC_rf_map, RGC_rf_weights, extras] = cone_to_RGC_connectivity(cone_locations, RGC_locations, varargin)

p = inputParser;

p.addRequired('cone_locations', @isnumeric)
p.addRequired('RGC_locations', @isnumeric)
p.addParamValue('center_amplitude', 10, @isnumeric)
p.addParamValue('center_space_constant', 2, @isnumeric)
p.addParamValue('center_weight_sd_factor', 0.25, @isnumeric)
p.addParamValue('surround', false, @islogical)
p.addParamValue('surround_weight_sd_factor', 0.25, @isnumeric)
p.addParamValue('surround_amplitude', 3, @isnumeric)
p.addParamValue('surround_space_constant', 0.5, @isnumeric)
p.addParamValue('sharing_probability', 0.1, @isnumeric)
p.addParamValue('seed', [], @isnumeric)

p.parse(cone_locations, RGC_locations, varargin{:});

center_amplitude = p.Results.center_amplitude;
center_space_constant = p.Results.center_space_constant;
center_weight_sd_factor = p.Results.center_weight_sd_factor;
surround = p.Results.surround;
surround_weight_sd_factor = p.Results.surround_weight_sd_factor;
surround_amplitude = p.Results.surround_amplitude;
surround_space_constant = p.Results.surround_space_constant;
seed = p.Results.seed;
sharing_probability = p.Results.sharing_probability;


% seed random number generator for reproducibility
if ~isempty(seed)
    rand('twister', 3333);
end

num_RGCs = length(RGC_locations(:,1));
num_cones = length(cone_locations(:,1));
RGC_rf_map = cell(num_RGCs, 1); % initialize cell array of connections
RGC_rf_weights = cell(num_RGCs, 1); % initialize cell array of weights
unsampled_cones = 0; % initialize counter

for cone = 1:num_cones
    clear weights
    temp_location = cone_locations(cone,:);
    temp_location = repmat(temp_location, num_RGCs, 1);
    temp_distances = sqrt((temp_location(:,1) - RGC_locations(:,1)).^2 + (temp_location(:,2) - RGC_locations(:,2)).^2);
    connection_potential = center_amplitude * exp(-(center_space_constant.^-1 * temp_distances).^2);    
    [weights, RGC_indices] = sort(connection_potential, 'descend');
    signif_indices = find(weights > 1);
    if isempty(signif_indices)
        unsampled_cones = unsampled_cones + 1;
    end
    if length(signif_indices) > 1
        rand_sample = rand(1);
        if rand_sample > sharing_probability
            signif_indices = signif_indices(1);
            weights = abs(weights(1) + normrnd(0,weights(1)*center_weight_sd_factor));  % add random strength to sampling weights
        else
            signif_indices = signif_indices(1:2);
            weights = abs(weights(1:2) + normrnd(0,weights(1)*center_weight_sd_factor, 2,1));
        end
    end
    
    for cll = 1:length(signif_indices);
        temp_rf = RGC_rf_map{RGC_indices(signif_indices(cll))};
        temp_rf = [temp_rf, cone];
        RGC_rf_map{RGC_indices(signif_indices(cll))} = temp_rf;
        temp_weights = RGC_rf_weights{RGC_indices(signif_indices(cll))};
        temp_weights = [temp_weights, weights(cll)];
        RGC_rf_weights{RGC_indices(signif_indices(cll))} = temp_weights;
    end

    
    if surround
        surround_connection_potential = surround_amplitude * exp(-(surround_space_constant * temp_distances).^2);
        [surround_weights, RGC_indices] = sort(surround_connection_potential, 'descend');
        surround_signif_indices = find(surround_weights > 0.1);
        excluded_center_indices = find(surround_signif_indices ~= signif_indices);
        surround_signif_indices = surround_signif_indices(excluded_center_indices);
        RGC_indices = RGC_indices(surround_signif_indices);
        surround_weights = surround_weights(surround_signif_indices);
        surround_weights = surround_weights + normnd(zeros(length(surround_weights)), surround_weights.*surround_weight_sd_factor);

        for cll = 1:length(surround_signif_indices);
            temp_rf = RGC_rf_map{RGC_indices(surround_signif_indices(cll))};
            temp_rf = [temp_rf, cone];
            RGC_rf_map{RGC_indices(surround_signif_indices(cll))} = temp_rf;
            temp_weights = RGC_rf_weights{RGC_indicies(surround_signif_indices(cll))};
            temp_weights = [temp_weights, -1*abs(surround_weights(cll))];
            RGC_rf_weights{RGC_indices(surround_signif_indices(cll))} = temp_weights;
        end
    end
end
  
extras.unsampled_cones = unsampled_cones;

