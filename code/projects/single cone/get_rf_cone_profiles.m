function [all_x, all_y] = get_rf_cone_profiles(datarun,cell_spec, varargin)
% get_rf_cone_profiles     get cone weights as a function of distance from center, aggregated across many RFs
%
% usage:  [all_x, all_y] = get_rf_cone_profiles(datarun,cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells
%            varargin - struct or list of optional parameters (see below)
%
% outputs:      all_x - Nx1 matrix, distance of each cone from the RF center
%               all_y - Nx1 matrix, weight of each cone
%
% NOTE: If by_cone_type is true, all_x and all_y are cell arrays of length 4,
%       where all_x{ii}, all_y{ii} are the distances and weights for a single cone type.
%       The order is L,M,S,U.
%
%
% optional params, their default values, and what they specify:
%
% radius            Inf        	region in which to get average profile
% center_type       'com'     	type of center point to use
% by_cone_type      false   	group results by cone type?
% normalize         'fit'       how to normalize the weights before returning.  the weights are divided by:
%                                   'fit' - (center_scale - surround_scale) for the fit stored in datarun.cones.rf_fits
%                                   'max' - the largest cone weight
%                                   'new fit' - center_scale of a single gaussian fit to the cell
%                                   'none' - don't normalize
% selection         struct      struct of parameters to pass to select_cone_weights
%                                   if an empty struct, use all cones
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius', Inf);
p.addParamValue('center_type','com');
p.addParamValue('by_cone_type',false);
p.addParamValue('normalize', 'fit', @(x)any(strcmpi(x,{'none','max','fit','new fit'})));
p.addParamValue('selection',struct);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% initialize variables to store profiles
if ~params.by_cone_type
    all_x = [];
    all_y = [];
else
    all_x = cell(4,1);
    all_y = cell(4,1);
end

% go through each cell
for cc =1:length(cell_indices)
    
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % get center point
    rf_ctr = rf_center(datarun,cell_id,params.center_type);
    
    % skip if there's no center point
    if isempty(rf_ctr)
        continue
    end
    
    
    % get profile for this cell
    
    % if no selection criteria, use all cones 
    if isempty(params.selection)
        [x,y,cone_indices] = rf_cone_profile(datarun.cones.weights(:,cell_index),...
            datarun.cones.centers,rf_ctr,'radius',params.radius);
    else
        % otherwise use only specified cones

        % select cones
        [weights, selection] = select_cone_weights(datarun,cell_id,params.selection);

        [x,y,cone_indices] = rf_cone_profile(weights(selection),...
            datarun.cones.centers(selection,:),rf_ctr,'radius',params.radius);
    end
    
    % normalize by y amplitude
    switch params.normalize
        case 'max'
            % divide by max
            y = y/max(y);
            
        case 'fit'
            % divide by max of current fit
            
            if isempty(datarun.cones.rf_fits{cell_index})
                continue
            end
            
            % get center scale - surround scale of current fit
            current_fit_scale = datarun.cones.rf_fits{cell_index}.center_scale - ...
                datarun.cones.rf_fits{cell_index}.surround_scale;
            
            % divide by this
            y = y/current_fit_scale;
            
        case 'new fit'
            % fit a single gaussian
            fit_params = fit_profile(x,y,'center_radius',20);
            
            % divide by the amplitude
            y = y/(fit_params.center_scale);
            
        case 'none'
        otherwise
            error('Normalization type ''%s'' not recognized.',params.normalize)
    end
    
    % add to the growing accumulation of profiles

    % if no color specificity
    if ~params.by_cone_type
        % include all cones indiscriminately
        all_x = [all_x; x];
        all_y = [all_y; y];
    else
        % otherwise, group by cone color
        
        % get list of which cones are which color
        L_cone_indices = datarun.cones.types(cone_indices)=='L';
        M_cone_indices = datarun.cones.types(cone_indices)=='M';
        S_cone_indices = datarun.cones.types(cone_indices)=='S';
        U_cone_indices = datarun.cones.types(cone_indices)=='U';
        
        % add to the accumulating list of numbers
        all_x{1} = [all_x{1}; x(L_cone_indices)];
        all_x{2} = [all_x{2}; x(M_cone_indices)];
        all_x{3} = [all_x{3}; x(S_cone_indices)];
        all_x{4} = [all_x{4}; x(U_cone_indices)];
        
        all_y{1} = [all_y{1}; y(L_cone_indices)];
        all_y{2} = [all_y{2}; y(M_cone_indices)];
        all_y{3} = [all_y{3}; y(S_cone_indices)];
        all_y{4} = [all_y{4}; y(U_cone_indices)];
        
    end

end


