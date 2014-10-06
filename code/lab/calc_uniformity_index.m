function ui = calc_uniformity_index(datarun, cell_spec, varargin)
% CALC_UNIFORMITY_INDEX   Calculate uniformity index for given mosaic
% (see Gauthier et al. 2009, DOI:10.1371/journal.pbio.1000063)
%
% usage: ui = calc_uniformity_index(datarun, cell_spec, opts)
%
% opts: thresh              []
%       min_thresh          0.05
%       verbose             false
%       plot                false
%       tform               []
%       simplify_contours   false
%
% 2010-02 phli
%

opts = inputParser;
opts.addParamValue('contours_field', 'rf_contours');
opts.addParamValue('thresh', []);
opts.addParamValue('min_thresh', 0.05);
opts.addParamValue('verbose', false);
opts.addParamValue('plot', false);
opts.addParamValue('tform', []);
opts.addParamValue('simplify_contours', false);
opts.parse(varargin{:});
opts = opts.Results;

if opts.verbose
    disp('Calculating uniformity index...')
end

if opts.thresh
    if opts.thresh < opts.min_thresh
        opts.thresh = opts.min_thresh;
    end
    
    if opts.verbose
        disp(['  Thresh: ' num2str(opts.thresh)]);
    end
    
    % Recalculate contours at new threshold level
    datarun = get_rf_contours(datarun, cell_spec, opts.thresh, 'norm_params', struct('method', 'peak'));
end

if ~isempty(opts.tform)
    % Should be fine to overwrite contours_field with transformed version
    % since this version of datarun is not returned
    datarun = tform_rf_contours(datarun, cell_spec, opts.tform, 'contours_field', opts.contours_field, 'output_field', opts.contours_field);
end

if opts.simplify_contours
    datarun = simplify_rf_contours(datarun, cell_spec);
end



[cell_nums, cell_type_name, cell_type_num] = get_cell_indices(datarun, cell_spec); %#ok<ASGLU> % cell_type_name is not used

contours = datarun.stas.(opts.contours_field);
num_contours = length(contours{cell_nums(1)}); % Bit of a hack...

for i = 1:num_contours
    % Initialize ui (needed in case contours turn out to be empty and
    % verbose output is on)
    ui = 0;
    
    % Initialize union polygon with first contour
    union_polygon = contours{cell_nums(1)}{i};

    % Polygon for storing overlaps starts empty
    overlap_polygon = [];

    % Settings for PolygonClip
    diff_mode = 0;
    intersect_mode = 1;
    union_mode = 3;

    for j = 2:numel(cell_nums)
        rf = contours{cell_nums(j)}{i};
        if isempty_polygon_struct(rf)
            continue;
        end

        % First find overlap of current rf with all those that came before
        curr_overlap_polygon = PolygonClip(union_polygon, rf, intersect_mode);

        % If there was any overlap, union it to the overlap total
        if ~isempty(curr_overlap_polygon)
            if ~isempty(overlap_polygon)
                overlap_polygon = PolygonClip(overlap_polygon, curr_overlap_polygon, union_mode);
            else
                overlap_polygon = curr_overlap_polygon;
            end
        end

        % Now add current rf to the union total and repeat
        union_polygon = PolygonClip(union_polygon, rf, union_mode);
    end

    % Subtract the overlap total from the union total to get area covered
    % by only one cell
    if ~isempty_polygon_struct(overlap_polygon)
        ui_polygon = PolygonClip(union_polygon, overlap_polygon, diff_mode);
    end
    
    % Intersect with ROI
    if ~isempty_polygon_struct(union_polygon)
        roi = datarun.stas.delaunays{cell_type_num}.roi;
        ui_polygon = PolygonClip(ui_polygon, roi, intersect_mode);
        
        ui(i) = calc_polygon_struct_area(ui_polygon) / calc_polygon_struct_area(roi);
        
        if opts.plot
            figure;
            plot_polygon_struct(union_polygon, 'b');
            hold on;
            plot_polygon_struct(overlap_polygon, 'r');
        end     
    end
end
if opts.verbose
    disp(['  UI: ' num2str(ui)]);
    disp(' ');
end
