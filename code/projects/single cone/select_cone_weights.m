function [weights, selection, extras] = select_cone_weights(datarun,cell_spec, varargin)
% select_cone_weights     choose cone weights based on various criterion
%
% usage:  [weights,selection,extras] = select_cone_weights(datarun,cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     weights - NxL double matrix of cone weights
%            selection - NxL binary matrix of which weights were selected
%               extras - extra information.  
%                           outputs a new datarun.cones structure with cones
%                           removed as specified by remove_cones
%
% where N = # cones, L = # cells
%
% optional params, their default values, and what they specify:
%
% thresh      	0.1           threshold based on weight amplitude
%                           if length==1, params.thresh is a lower bound, and
%                           there is no upper bound
%                           if any value of params.thresh is strictly greater than 1,
%                               the threshold is in units of SNR range, where
%                               noise amplitude is computed as robust_std(all weights)
%                           otherwise
%                               the threshold is a fraction of the largest weight
%                           e.g. thresh==5 means all weights > 5*noise sigma
%                                thresh==0.2 means all weights > 0.2*largest weight
%                                thresh==[1 5] means all weights between 1 and 5 noise sigmas
%
% radius        [0 Inf]     range in units of center sigmas
%
% polarity     	1           which polarity to return
%                           	1 - positive only
%                              -1 - negative only
%                               0 - either
%
%
% 
% contiguity    true        enforces selected cones to be contiguous by
%                           clustering cones around the strongest cone, cutoff=scale*(median nnd)  
%
% scale         3           maximum separation between cones, expressed in number of median nnd, 
%                           to be considered contiguous
%
%
% remove_cones  []          a character string that defines cones to be
%                           excluded from the output weights and selection matrix
%                           i.e. 'SU' will remove 'S' and 'U' cones
%
% 2009-02 gauthier
% edits: 2009-07 gdf
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% forked parameters

% get criterion type
p.addParamValue('thresh', 0.1,@(x)length(x)<=2);
p.addParamValue('radius', [0 Inf],@(x)length(x)==2);
p.addParamValue('polarity', 0, @(x)numel(x)==1 & any(x==[-1 0 1]));
p.addParamValue('contiguity', true, @islogical);
p.addParamValue('scale', 3, @isnumeric);
p.addParamValue('remove_cones', [], @ischar)

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% BODY OF THE FUNCTION


% parse threshold
if any(params.thresh > 1)
    thresh_type = 'snr';
else
    thresh_type = 'peak';
end
% establish upper bound, if nonexistant
if length(params.thresh) == 1, params.thresh = [params.thresh Inf]; end


% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);


% remove specified cones
if ~isempty(params.remove_cones)
    for tp = 1:length(params.remove_cones)
        temp_type(tp) = params.remove_cones(tp);
        temp_exclude_cone_indices = find(datarun.cones.types == temp_type(tp));
        keep_cone_indices = setdiff([1:length(datarun.cones.types)], temp_exclude_cone_indices);
        
        datarun.cones.centers = datarun.cones.centers(keep_cone_indices,:);
        datarun.cones.types = datarun.cones.types(keep_cone_indices);
        datarun.cones.rgb = datarun.cones.rgb(keep_cone_indices,:);
        datarun.cones.likelihoods = datarun.cones.likelihoods(keep_cone_indices);
        datarun.cones.types_em = datarun.cones.types_em(keep_cone_indices);
        datarun.cones.types_kmeans = datarun.cones.types_kmeans(keep_cone_indices);
        datarun.cones.roi = datarun.cones.roi(keep_cone_indices);
        datarun.cones.weights = datarun.cones.weights(keep_cone_indices,:);
    end
end


% initialize output variables
all_weights = datarun.cones.weights(:,cell_indices);
all_selections = false(size(datarun.cones.weights,1),length(cell_indices));


% loop through cells
for cc = 1:length(cell_indices)

    % get cell index
    cell_index = cell_indices(cc);
    

    % get weights
    weights = datarun.cones.weights(:,cell_index);

    % intially, set all weights to be selected, and subsequently pare them down based on the criteria
    selection = true(size(weights));


    % threshold
    switch thresh_type
        case 'snr'
            % threshold based on SNR
            % identify noise
            noise_sigma = robust_std(weights);
            % keep only weights inside the range
            selection = selection & (abs(weights)>=noise_sigma*params.thresh(1)) & ...
                (abs(weights)<=noise_sigma*params.thresh(2));
        case 'peak'
            % identify maximum weight
            max_weight = max(weights);
            % threshold based on largest weight
            selection = selection & (abs(weights)>=params.thresh(1)*max_weight) & (abs(weights)<=params.thresh(2)*max_weight);
            
    end

    % only check if range is not [0 Inf]
    if ~all(params.radius == [0 Inf])
        % identify rf center and center sigma, using the single cone RF fit
        center = datarun.cones.rf_fits{cell_index}.center;
        radius = datarun.cones.rf_fits{cell_index}.center_radius;
        % get distances from RF center
        distances = sqrt(sum(    ( datarun.cones.centers - repmat(center,size(datarun.cones.centers,1),1) ).^2   ,2));
        % keep only weights inside the range
        selection = selection & (distances>=radius*params.radius(1)) & (distances<=radius*params.radius(2));
    end

    % polarity
    if params.polarity ~= 0  % if 0, keep all weights
        % otherwise, remove weights of the undesired sign
        selection = selection & (weights*params.polarity > 0);
    end

    % enter saved information into output structures
    all_selections(:,cc) = selection;
    
end


% if contiguity
if params.contiguity
    all_selections = enforce_cone_contiguity(datarun, all_weights, all_selections, params.scale);
end


% store the fits only for the specified RGCs
if isfield(datarun.cones, 'rf_fits')
    for RGC = 1:length(cell_indices);
        temp_rf_fits{RGC} = datarun.cones.rf_fits{cell_indices(RGC)};
    end
    datarun.cones.rf_fits = temp_rf_fits;
end


weights = all_weights;
selection = all_selections;
extras.new_datarun.cones = datarun.cones;

