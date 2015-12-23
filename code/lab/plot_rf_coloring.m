function outmosaic = plot_rf_coloring(datarun, cell_spec, varargin)
% PLOT_RF_COLORING    Do a topocoloring of mosaic
%
% usage: mosaic = plot_rf_coloring(datarun, cell_spec, opts)
%
% opts: skip            []                      Cell ids within cell_spec to skip
%       rfs             'rfs'                   Name of field in datarun.stas holding the RFs
%       norm_params     {'method', 'peak'}      Params to pass to RF_NORMALIZED
%       clip            0.25                    Threshold level below which to clip normalized RFs to zero
%
%       tform           []                      Affine transformation matrix to pass to BUILD_RF_TFORM.  A 
%                                               scalar will be interpreted as degrees clockwise rotation.  
%                                               Can also use keywords 'xreflect' and 'yreflect'
%
%       ifactor         2                       Number of times to consecutively 2D interpolate
%       imethod         'cubic'                 Interpolation method for interp2
%
% 2010-02 phli
%

% Sequence of colors to use for topocoloring; in practice only 5-6 are
% generally needed; can add to this if we ever need more.
colors = [ 
    1 0 0;       % Red
    0 1 0;       % Green
    0 0 1;       % Blue
    0 1 1;       % Cyan
    1 0 1;       % Purple
    1 1 0;       % Yellow
    1 1 1;       % White/Grey
];


opts = inputParser;
opts.addParamValue('skip', []);
opts.addParamValue('rfs', 'rfs');
opts.addParamValue('norm_params', {'method', 'peak'});
opts.addParamValue('clip', 0.25);       % Threshold level below which to clip normalized rfs to zero
opts.addParamValue('tform', []);        % Affine transformation matrix to pass to BUILD_RF_TFORM
opts.addParamValue('ifactor', 2);       % Number of times to consecutively 2D interpolate the results
opts.addParamValue('imethod', 'cubic'); % Interpolation methods for interp2
opts.addParamValue('label', false);     % Show cell_id labels on each RF
opts.addParamValue('label_size', 10);
opts.parse(varargin{:});
opts = opts.Results;


[cell_type_num, cell_type_name] = get_cell_type_nums(datarun, cell_spec);
cell_ids = datarun.cell_types{cell_type_num}.cell_ids;

if ~isfield(datarun.stas, 'delaunays') || cell_type_num > length(datarun.stas.delaunays) || isempty(datarun.stas.delaunays{cell_type_num})
    datarun = do_delaunay_tri(datarun, cell_spec);
end
dt = datarun.stas.delaunays{cell_type_num};

if ~isfield(dt, 'coloring')
    % Get a topocoloring of the mosaic through greedy color picking on Delaunay Triangulation
    dt.coloring = delaunay_color(dt);
end


% Look for skip data in datarun if not given in args
if isempty(opts.skip) && isfield(datarun.stas, 'skip') && (numel(datarun.stas.skip) >= cell_type_num)
    opts.skip = datarun.stas.skip{cell_type_num};
    disp(['Skipping: ' num2str(datarun.stas.skip{cell_type_num})]);
end


% Build up the mosaic image by overlaying each topocolored rf
rf_size = size(get_rf(datarun, cell_ids(1)));
mosaic = zeros(rf_size(1), rf_size(2), 3);
for i = 1:numel(cell_ids)
    if any(opts.skip == cell_ids(i))
        continue;
    end
    
    if dt.coloring(i) == 0
        continue;
    end
    
    rf = get_rf(datarun, cell_ids(i), 'where', opts.rfs);
    rf = sum(rf, 3); % Sum rf over RGB dimension, if any
    %rf = sqrt(sum(rf.^2, 3)); % Get vector norm of each pixel's RGB values
    
    % Normalize rf
    rf = rf_normalized(rf, opts.norm_params{:});
    
    % Clip rf at threshold level
    rf(rf < opts.clip) = 0;
    
    % Renormalize
    rf(rf > 0) = rf(rf > 0) -  min(rf(rf > 0));
    rf = rf ./ max(rf(:));

    % Apply image transformation, if any
    if ~isempty(opts.tform)
        rf = imtransform(rf, build_rf_tform(datarun, cell_ids(i), opts.tform), 'XData', [1 size(rf, 2)], 'YData', [1 size(rf, 1)]);
    end
    
    % Set rf color and add to mosaic image
    color = colors(dt.coloring(i), :);
    for j = 1:3
        mosaic(:,:,j) = mosaic(:,:,j) + (rf .* color(j));
    end
end

if opts.ifactor > 0
    % Interpolate mosaic image
    mosaici        = interp2(mosaic(:,:,1), opts.ifactor, opts.imethod);
    mosaici(:,:,2) = interp2(mosaic(:,:,2), opts.ifactor, opts.imethod);
    mosaici(:,:,3) = interp2(mosaic(:,:,3), opts.ifactor, opts.imethod);
    mosaic = mosaici;
end

% Normalize mosaic
mosaic = mosaic -  min(mosaic(:));
mosaic = mosaic ./ max(mosaic(:));

image(mosaic);
axis equal;
title([strrep(datarun.names.short_name, '_', ' ') ': ' cell_type_name ' ' num2str(opts.tform)]);




% Show cell_id labels?
if opts.label
    hold on;
    
    % If we interpolated the images, we have to hack the coordinates to
    % scale them up the right number of factors
    input_scale = 1 / (2.^opts.ifactor);
    plot_rf_labels(datarun, cell_spec, 'input', 'sta scaled', 'scale', input_scale, 'label_color', 'white', 'label_size', opts.label_size);
end



% To facilitate saving figures automatically
savename = ['rfcolor_' datarun.names.short_name '_ct' num2str(cell_type_num)];
if ~isempty(opts.tform)
    savename = [savename '_' num2str(opts.tform)];
end
set_savename(savename);



if nargout > 0
    outmosaic = mosaic;
end
