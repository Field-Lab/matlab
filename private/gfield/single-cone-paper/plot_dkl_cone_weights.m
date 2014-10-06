function OIs = plot_dkl_cone_weights(datarun, cell_spec, varargin)

% plot_dkl_cone_weights     This function plots the cone weights to a cell in DKL weight space
%
% usage:  plot_dkl_cone_weights(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct with field specifying X
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field X
%
%
% optional parameters, their default values, and what they specify:
%
%
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot.  if -1, plot in current.
% method            time_course             information to base cone weights on
%                                           'time_course' - time course of response -- for LMS stimuli only!
%                                           'cones' - use cone weight information
%                                           'area_around_peak', a window around the peak
%
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not passed.
%
%
%       passed to 'select_cone_weights'
%
%           thresh          2.5     	threshold to include cone weights
%           radius          [0 6]     	radius for select cone weights
%           polarity         0          polarity of cones to use
%           contiguity      false       contiguity condition
%           scale           3.0         
%           remove_cones      'U'       remove "undefined" cones
%
%
%
%
%
% date author
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa', []);
p.addParamValue('method','time_course', @ischar);
p.addParamValue('permute_cones', false, @islogical); 

% parameters to be passed on
%    select_cone_weights function
p.addParamValue('thresh', 2.5, @isnumeric);
p.addParamValue('radius', [0 6], @isnumeric);
p.addParamValue('polarity', 0, @isnumeric);
p.addParamValue('contiguity', false, @islogical');
p.addParamValue('scale', 3.0);
p.addParamValue('remove_cones', 'U', @ischar);
p.addParamValue('window_size', 2, @isnumeric);
p.addParamValue('S_contrast_factor', 1, @isnumeric);
p.addParamValue('sigma_cutoff', 2, @isnumeric);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
selection_params.thresh = params.thresh;
selection_params.radius = params.radius;
selection_params.polarity = params.polarity;
selection_params.contiguity = params.contiguity;
selection_params.scale = params.scale;
selection_params.remove_cones = params.remove_cones;

% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa);
axes(plot_axes)


% BODY OF FUNCTION

% get cell indices and number of RGCs
cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);

% extract polarities of cells
polarities = datarun.stas.polarities;

% initalize matrix that will be filled in during loop
weighted_triplets = zeros(num_rgcs, 3);

if strcmp(params.method, 'cones')

    % extract connectivity
    selection_params
    [weights, selection, extras] = select_cone_weights(datarun, cell_spec, selection_params); 
    connectivity = weights .* selection;
    
    % extract cone types
    cone_types = extras.new_datarun.cones.types;

    % randomly permute cone labels
    if params.permute_cones
        num_cones = length(weights(:,1));
        shuffled_indices = randperm(num_cones);
        cone_types = cone_types(shuffled_indices);
    end

    % get indicies to L, M, and S cones.
    L_cone_indices = find(cone_types == 'L');
    M_cone_indices = find(cone_types == 'M');
    S_cone_indices = find(cone_types == 'S');

    % initialize variables that will be filled in during loop
    %temp_opponent_indices = zeros(num_rgcs, 1);

    for rgc = 1:num_rgcs
        temp_cone_indices = find(connectivity(:,rgc));

        temp_L_indices = intersect(temp_cone_indices, L_cone_indices);
        temp_M_indices = intersect(temp_cone_indices, M_cone_indices);
        temp_S_indices = intersect(temp_cone_indices, S_cone_indices);

        total_L = sum(connectivity(temp_L_indices, rgc) * polarities{cell_indices(rgc)});
        total_M = sum(connectivity(temp_M_indices, rgc) * polarities{cell_indices(rgc)});
        total_S = sum(connectivity(temp_S_indices, rgc) * polarities{cell_indices(rgc)});

        total_input = abs(total_L) + abs(total_M) + abs(total_S);

        norm_L = total_L ./ total_input;
        norm_M = total_M ./ total_input;
        norm_S = total_S ./ total_input;

        weighted_triplets(rgc,:) = [norm_L, norm_M, norm_S];

    %     if (norm_M * norm_L) < 0 && (norm_M * polarities{cell_indices(rgc)}) > 0
    %         temp_opponent_indices(rgc) = datarun.cell_ids(cell_indices(rgc));
    %     end

    end

    %tmp_indices = find(temp_opponent_indices);
    %opponent_ids = temp_opponent_indices(tmp_indices)
end


if strcmp(params.method, 'time_course')
    for rgc = 1:num_rgcs
        temp_tc = datarun.stas.time_courses{cell_indices(rgc)};


        if isempty(temp_tc)
            error(['Cell ', num2str(datarun.cell_ids(cell_indices(rgc))),' has no time course'])
        end
        
        L_tc = temp_tc(:,1) - temp_tc(1,1);
        M_tc = temp_tc(:,2) - temp_tc(1,2);
        S_tc = temp_tc(:,3) - temp_tc(1,3);
    
        % make TC template
        L_max = max(L_tc);
        L_min = min(L_tc);
        if L_max < abs(L_min)
            L_tc_temp = L_tc .*-1;
        else
            L_tc_temp = L_tc;
        end
        M_max = max(M_tc);
        M_min = min(M_tc);
        if M_max < abs(M_min)
            M_tc_temp = M_tc .* -1;
        else
            M_tc_temp = M_tc;
        end
        S_max = max(S_tc);
        S_min = max(S_tc);
        if S_max < abs(S_min)
            S_tc_temp = S_tc *-1;
        else
            S_tc_temp = S_tc;
        end

        template = (L_tc_temp + M_tc_temp + S_tc_temp) ./ 3;
        
        L_proj = dot(L_tc, template);
        M_proj = dot(M_tc, template);
        S_proj = dot(S_tc, template);

        total_proj = abs(L_proj) + abs(M_proj) + abs(S_proj);
        
        norm_L = L_proj ./ total_proj;
        norm_M = M_proj ./ total_proj;
        norm_S = S_proj ./ total_proj;

        weighted_triplets(rgc,:) = [norm_L, norm_M, norm_S];
    end
end

if strcmp(params.method, 'area_around_peak')
    for rgc = 1:num_rgcs
        temp_rf = datarun.stas.rfs{cell_indices(rgc)};

%         imagesc(norm_image(temp_rf))
 %       rgc
% pause
        window_size = params.window_size;
        polarity = datarun.stas.polarities{cell_indices(rgc)};

        % get peak stixel in RF
        [red_peak, red_peak_index] = max(reshape(abs(squeeze(temp_rf(:,:,1))), 1, []));
        [green_peak, green_peak_index] = max(reshape(abs(squeeze(temp_rf(:,:,2))), 1, []));
        
        % check that green and red peaks are in same location
        % otherwise choose point with greater peak
        if green_peak > red_peak
            peak_index = green_peak_index;
        else
            peak_index = red_peak_index;
        end

        [peak_rw peak_cl] = ind2sub(size(squeeze(temp_rf(:,:,1))), peak_index);

        rw = peak_rw-window_size:1:peak_rw+window_size;
        cl = peak_cl-window_size:1:peak_cl+window_size;

        % ensure rows and column indices are within the bounds of the stimulus
        rw(rw<1) = 1;
        cl(cl<1) = 1;
        rw(rw>datarun.stimulus.field_height) = datarun.stimulus.field_height;
        cl(cl>datarun.stimulus.field_width) = datarun.stimulus.field_width;
        rw = unique(rw);
        cl = unique(cl);


        windowed_rf = temp_rf(rw, cl, :);
%imagesc(norm_image(windowed_rf))
%pause

        % scale the blue channel by the S cone contrast factor (something a different S cone contrast was used)
        windowed_rf = windowed_rf(:,:,3) ./ params.S_contrast_factor;
        glbl_extreme = max(max(max(abs(windowed_rf))));

        % get a robust sigma for each channel
        red_std = robust_std(reshape(squeeze(temp_rf(:,:,1))./glbl_extreme, 1, [])) .* params.sigma_cutoff;
        green_std = robust_std(reshape(squeeze(temp_rf(:,:,2))./glbl_extreme, 1, [])) .* params.sigma_cutoff;
        blue_std = robust_std(reshape(squeeze(temp_rf(:,:,3))./glbl_extreme./params.S_contrast_factor, 1, [])) .* params.sigma_cutoff;


        % normalize the r,g,and b values by the global extreme
        red_rf = squeeze(temp_rf(rw, cl, 1)) ./ glbl_extreme;
        green_rf = squeeze(temp_rf(rw, cl, 2)) ./ glbl_extreme;
        blue_rf = squeeze(temp_rf(rw, cl, 3)) ./ glbl_extreme ./ params.S_contrast_factor;

        % get indices to values that excede sigma cutoff
        red_signal_indices = find(abs(red_rf) > red_std);
        green_signal_indices = find(abs(green_rf) > green_std);
        blue_signal_indices = find(abs(blue_rf) > blue_std);

        red_sum = sum(sum(red_rf(red_signal_indices)));
        green_sum = sum(sum(green_rf(green_signal_indices)));
        blue_sum = sum(sum(blue_rf(blue_signal_indices)));


        red = sum(sum(abs(red_rf(red_signal_indices))));
        green = sum(sum(abs(green_rf(green_signal_indices))));


        total_input = abs(red_sum) + abs(green_sum) + abs(blue_sum);
        norm_L = red_sum ./ total_input .* polarity;
        norm_M = green_sum ./ total_input .* polarity;
        norm_S = blue_sum ./ total_input .* polarity;

        weighted_triplets(rgc,:) = [norm_L, norm_M, norm_S];

        
        OIs(rgc) = (red_sum - green_sum) ./ (red + green);
    end
end

% figure(11)        
% hist(weighted_triplets(:, 1))
% figure(12)
% hist(weighted_triplets(:,2))
% figure(13)
% hist(weighted_triplets(:,3))

hold on
plot(weighted_triplets(:,2), weighted_triplets(:,1), 'ko')
plot([0 0], [-1 1], 'k',[-1 1], [0 0], 'k') % axes
plot([0 -1], [1, 0], 'k', [0 1], [-1 0], 'k') % opponent diagonals
plot([1 0], [0, 1], 'k', [-1 0], [0, -1], 'k') %non-opponent diagonals

% s cone axis
plot([0, 1/sqrt(4)], [0, 1/sqrt(4)], 'k')
plot([-0.05, 0.05] + [1/4, 1/4], [0.05, -0.05] + [1/4, 1/4], 'k')
xlabel('M cone')
ylabel('L cone')
axis square
hold off









