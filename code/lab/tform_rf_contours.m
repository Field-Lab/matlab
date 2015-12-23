function datarun = tform_rf_contours(datarun, cell_spec, tform, varargin)
% TFORM_RF_CONTOURS
% usage : datarun = tform_rf_contours(datarun, cell_spec, tform, [opts])
%
% opts:     contours_field
%           contour_indices
%
% phli 2010-04
%

opts = inputParser;
opts.addParamValue('contours_field', 'rf_contours');
opts.addParamValue('output_field', []);
opts.addParamValue('output_suffix', []);
opts.addParamValue('contour_indices', []);
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.output_suffix)
    opts.output_suffix = 'tform';
end

if isempty(opts.output_field)
    opts.output_field = [opts.contours_field '_' opts.output_suffix];
end



[cell_nums, cell_type_name, cell_type_num] = get_cell_indices(datarun, cell_spec);

for cell_num = cell_nums
    % Build up tform if needed
    if ~isstruct(tform)
	    rf_tform = build_rf_tform(datarun, datarun.cell_ids(cell_num), tform);
	else
	    rf_tform = tform;
    end

    if isempty(opts.contour_indices)
        contour_indices = 1:length(datarun.stas.(opts.contours_field){cell_num});
    else
        contour_indices = opts.contour_indices;
    end

    
	% Apply tform to contours
    for i = contour_indices
        contour = datarun.stas.(opts.contours_field){cell_num}{i};
        contour = tform_polygon_struct(contour, rf_tform);
        datarun.stas.(opts.output_field){cell_num}{i} = contour;
    end
end