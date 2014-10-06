function net_inputs = get_net_cone_input(datarun,cell_spec, varargin)
% get_net_cone_input     get sum of all cone weights, sorted by type, within a specified radius
%
% usage:  net_inputs = get_net_cone_input(datarun,cell_spec, varargin)
%
% arguments:    datarun - datarun struct
%             cell_spec - which cells to use
%              varargin - struct or list of optional parameters (see below)
%
% outputs:     net_inputs - Nx3 matrix of the sum of cone weights for each cell (rows) and cone type (columns: L,M,S)
%
%
% optional params, their default values, and what they specify:
%
% radius_type       'pixels'            the units of the radius
%                                           'pixels'
%                                           'center'
%                                           'surround'
% radius            Inf               	radius in which to gather cone weights
%
%
%
% parameters passed on to get_rf_cone_profiles, their name when passed, and what they do.
%       normalize       normalize       how to normalize the amplitude of cone weights
%       center_type     center_type     type of center point around which to get cones
%
%
% gauthier 2008-10
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('radius', Inf);
p.addParamValue('radius_type', 'pixels',@(x)any(strcmpi(x,{'pixels','center','surround'})));

% parameters to be passed on to get_rf_cone_profiles
p.addParamValue('normalize', 'default value');
p.addParamValue('center_type', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% make parameters to pass
profile_params = make_struct_to_pass(p.Results,{'normalize','normalize','center_type','center_type'});



% BODY OF THE FUNCTION


% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% initialize output variable
net_inputs = zeros(length(cell_indices),3);

for cc = 1:length(cell_indices)

    cell_index = cell_indices(cc);
    
    % don't get input if no fit exists
    if isempty(datarun.cones.rf_fits{cell_index})
        continue
    end

    % get radius in pixels
    switch params.radius_type

        case 'pixels'
            radius = params.radius;
        
        case {'center','surround'}
            % don't get input if no fit exists
            if isempty(datarun.cones.rf_fits{cell_index})
                continue
            end
            
            switch params.radius_type
                case 'center'
                    radius = params.radius * abs(datarun.cones.rf_fits{cell_index}.center_radius);

                case 'surround'
                    radius = params.radius * datarun.cones.rf_fits{cell_index}.surround_radius;
            end

        otherwise
            error('radius type ''%s'' not recognized',params.radius_type)
    end

    % get cone weights within the specified radius
    [junk, weight_array] = get_rf_cone_profiles(datarun,datarun.cell_ids(cell_index),...
        'radius',radius,'by_cone_type',true,profile_params);

    % load weights of each cone type
    for tt = 1:3
        net_inputs(cc,tt) = sum(weight_array{tt});
    end
end




