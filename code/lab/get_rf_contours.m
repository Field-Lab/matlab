function datarun = get_rf_contours(datarun, cell_specification, contour_levels, varargin)
% get_rf_contours     Compute STA contours for several cells and save to dataset
%
% usage:  datarun = get_rf_contours(datarun, cell_specification, contour_levels, params)
%
% arguments:          datarun - datarun struct with field dataset.stas.rfs
%          cell_specification - see get_cell_indices for options
%              contour_levels - levels at which to compute contour
%                      params - struct of optional parameters (see below)
%
% outputs:            datarun - datarun struct with fields:
%
%   datarun.contours{}{}    or another name
%
%
% optional fields in params, their default values, and what they specify:
%
% verbose           false               display info
% norm_params       []                  type of normalization to apply, see 'rf_normalized' for options
% save_cont_params  []                  save parameters to datarun.stas.(save_cont_params)
%                                           if empty, params are not saved.
%                                           this struct contains normalization parameters,
%                                           contour levels, and name of summary frame
% save_name         'rf_contours'          name of field in which to save contours
% rfs               'summaries_filt'    which summaries to use for computing contours
%              
%
%   gauthier 2008-03
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('norm_params', struct);
p.addParamValue('save_cont_params', []);
p.addParamValue('save_name', 'rf_contours');
p.addParamValue('rfs', 'summaries_filt');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;





% get list of cells
[cell_nums, cell_type_name, cell_type_num] = get_cell_indices(datarun, cell_specification);


% display how many cells will be loaded
if params.verbose;
    fprintf('\nCompute contours %s for %d %s cells ', num2str(contour_levels), length(cell_nums), cell_type_name);
    start_loading = clock; % note when it started
end


% go through through list of cells
for cell_num = cell_nums
    
    % show a tick of progress
    if params.verbose;fprintf('.');end
    
    %[filt_summary, filt_params] = rf_filtered(datarun.stas.rfs{cell_num}, struct('filt_type',params.filt_type));
    
    % get summary
    summary = datarun.stas.(params.rfs){cell_num};
    
    % normalize summary frame
    [norm_summary,extras,norm_params] = rf_normalized(summary, params.norm_params);
    
    % compute contours of normalized summary
    contour_polygons = rf_contours(norm_summary, contour_levels);
    
    % save to datarun in the specified field
    datarun.stas.(params.save_name){cell_num} = contour_polygons;
    
    % save parameters which defined the contours, if desired
    if ~isempty(params.save_cont_params)
        
        % note the normalization parameters
        cont_params = norm_params;
        
        % note which contour levels
        cont_params.levels = contour_levels;
        
        % note which summary frame was used
        cont_params.rfs = params.rfs;
        
        % save in dataset
        datarun.stas.(params.save_cont_params){cell_num} = cont_params;
    end
end

% if not all cells have contours, ensure datarun.stas.(params.save_name) is the right length
if length(datarun.stas.(params.save_name)) < length(datarun.cell_ids)
    datarun.stas.(params.save_name){length(datarun.cell_ids)} = [];
end

% display how long it took
if params.verbose;
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end
