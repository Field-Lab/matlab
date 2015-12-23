function [all_x, all_y] = get_rf_profiles(datarun,cell_spec, varargin)
% get_rf_profiles     get cone weights as a function of distance from center, aggregated across many RFs
%
% usage:  [all_x, all_y] = get_rf_profiles(datarun,cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells
%            varargin - struct or list of optional parameters (see below)
%
% outputs:      all_x - Nx1 matrix, distance of each cone from the RF center
%               all_y - NxC matrix, profile strength.  C = size(rf,3).
%
%
%
% optional params, their default values, and what they specify:
%
% radius            Inf        	region in which to get average profile
% center_type       'com'     	type of center point to use
% normalize         'std'       how to normalize the stixels before returning.  the stixels are divided by:
%                                   'none' - don't normalize
%                                   'max' - the largest stixel
%                                   'std' - standard deviation of the stixel
%
%
% gauthier 2008-12
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius', Inf);
p.addParamValue('center_type','com');
p.addParamValue('normalize', 'std', @(x)any(strcmpi(x,{'none','max','std'})));

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% initialize variables to store profiles
all_x = [];
all_y = [];

% go through each cell
for cc =1:length(cell_indices)
    
    cell_index = cell_indices(cc);
    
    % get center point
    rf_ctr = rf_center(datarun,datarun.cell_ids(cell_index),params.center_type);
    
    % get rf
    rf = get_rf(datarun,datarun.cell_ids(cell_index));
    
    % skip if there's no center point, or no rf
    if isempty(rf_ctr)
        continue
    end
    
    % get profile
    [x,y] = rf_profile(rf,rf_ctr,'radius',params.radius);
    
    % normalize by y amplitude
    switch params.normalize
        case 'max'
            % divide by max
            y = y/max(max(y));
            
        case 'std'
            % divide by std
            y = y/std(reshape(y,[],1));
            
        case 'none'
        otherwise
            error('Normalization type ''%s'' not recognized.',params.normalize)
    end
    
    % add to the growing accumulation of profiles
    all_x = [all_x; x];
    all_y = [all_y; y];

end


