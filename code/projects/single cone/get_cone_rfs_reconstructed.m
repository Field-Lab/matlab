function rfs = get_cone_rfs_reconstructed(datarun,cell_spec, varargin)
% get_cone_rfs_reconstructed     return cell array of reconstructed cone RFs
%
% usage:  rfs = get_cone_rfs_reconstructed(datarun,cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     rfs - cell array of matrices
%                       each matrix is YxXxC, where C = 1 or 3 (see below)
%
%
% optional params, their default values, and what they specify:
%
%  by_type      true        color cones according to their type
%                               if true, output matrix has 3 colors
%                               if false, output matrix has 1 color
%  Wc           []          cone weight matrix
%                               if empty, cones are approximated as a single pixel
%
%
% 2008-11 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('by_type', true);
p.addParamValue('Wc', []);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

rfs = cell(length(cell_indices),1);

% get reconstructed RF from each
for cc = 1:length(cell_indices)
   
    % get cone weights
    the_weights = datarun.cones.weights(:,cell_indices(cc));
    
    if isempty(params.Wc)
        % call cone_rf_reconstructed to approximate each cone as a single pixel
        
        % get parameter to pass to cone_rf_reconstructed
        if params.by_type
            cone_types = datarun.cones.types;
        else
            cone_types = [];
        end

        % get single cone RF
        rfs{cc} = cone_rf_reconstructed([datarun.stimulus.field_height datarun.stimulus.field_width],...
            the_weights,datarun.cones.centers,'cone_types',cone_types);
    
    else
        % use Wc to reconsruct
        rfs{cc} = reshape(params.Wc*the_weights,datarun.stimulus.field_height,datarun.stimulus.field_width,[]);
        
    end
end

