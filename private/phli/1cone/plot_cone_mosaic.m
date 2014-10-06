function plot_axes = plot_cone_mosaic(datarun, varargin)
% plot_cone_mosaic     Plot the mosaic of cones
%
% usage:  plot_cone_mosaic(datarun, varargin)
%
% arguments:  datarun - datarun struct with fields stimulus.stixel_height, width
%              varargin - struct of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% fig_or_axes       0               which figure or axes to plot in.  if 0, make new.
% cone_size         15              size of the cones to plot
% bg_color          [.5 .5 .5]      RGB triplet of the background color
%                                       if empty, don't plot a background image
% cone_colors       []              Nx3 matrix, N = number of cones, specifying color triplet for each cone
%                                       if empty, use cone_roi.  if cone_roi is empty, use type_colors
% cone_roi          []              boolean vector or list numbers specifying cones in ROI to use / highlight
% roi_highlight     true            if true, plot roi cones blue, and others black
%                                   if false, only plot roi cones and plot as specified by cone colors
%
% type_colors       4x3 matrix of awesome cone colors
%                                   4x3 matrix of colors, each row specifying the color for L, M, S, or U cones.
%                                       only used if cone_colors is empty.
% clear             true            clear axes before plotting
% drawnow           true            draw after completing?
% cone_circles      false           uses drawCircle to draw cones as circles of radius specified by cone_radius
% cone_size         []              Either marker size or radius if circles
% label             false           Print cone ids
% label_color       'black'         
% label_size        10
%
% 2008-10 gauthier
% 2012-07 phli, cleanup, can treat individual cones more flexibly if desired
%


p = inputParser;

p.addParamValue('fig_or_axes', 0);
p.addParamValue('clear', true);
p.addParamValue('drawnow', false, @islogical);
p.addParamValue('scale', [1 1]);
p.addParamValue('ticks', false);

p.addParamValue('bounds', []);
p.addParamValue('cone_roi', []);
p.addParamValue('roi_highlight', true, @islogical);

p.addParamValue('bg_color',[1 1 1]*0.5);
p.addParamValue('cone_size', []);
p.addParamValue('cone_colors', []);
p.addParamValue('type_colors', [1 0 0; 0 1 0; 0 0 1; 0 0 0]);
p.addParamValue('cone_circles', false, @islogical);
p.addParamValue('label', false);
p.addParamValue('label_color', 'black');
p.addParamValue('label_size' , 10);

p.parse(varargin{:});
params = p.Results;


ncones = length(datarun.cones.types);

% Convert list of numbers to boolean vector if needed
if ~islogical(params.cone_roi)
    boolvect = false(ncones,1);
    boolvect(params.cone_roi) = true;
    params.cone_roi = boolvect;
end

if isempty(params.cone_size)
    if params.cone_circles, params.cone_size = 1; 
    else                    params.cone_size = 15; end
end

% Functional approach
if params.cone_circles, cone_func = @circle_cones;
else                    cone_func = @marker_cones; end


% Set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes,params.clear);
if ~isempty(params.bg_color), set(plot_axes, 'Color', params.bg_color); end
if ~params.ticks, set(plot_axes, 'box', 'on', 'xtick', [], 'ytick', []); end
hold on


% Highlights; this is a bit futzy for the sake of backward compatibility
if ~isempty(params.cone_roi) && params.roi_highlight  && ncones>length(params.cone_roi)
    % If highlight is true, then we actually want to include all cones in
    % roi but just set colors differently.
    params.cone_colors = zeros(ncones, 3);
    params.cone_colors(params.cone_roi,3) = 1;
    params.cone_roi = true(ncones,1);
elseif isempty(params.cone_roi)    
    % Although highlight defaults to true, if the cone_roi is left empty
    % then assume that highlight is not actually wanted.  This is for
    % backward compatibility.
    params.cone_roi = true(ncones,1);
else
    % plot all cones with specified color - ath
    params.cone_colors = repmat(params.cone_colors,ncones,1);
    params.cone_roi = true(ncones,1);
end


% Bounds provided?
if ~isempty(params.bounds)
    inbounds = params.bounds(1) <= datarun.cones.centers(:,1) & datarun.cones.centers(:,1) <= params.bounds(2) & ...
               params.bounds(3) <= datarun.cones.centers(:,2) & datarun.cones.centers(:,2) <= params.bounds(4);
    params.cone_roi = params.cone_roi & inbounds;
end


% Put variables into shorter names
cone_centers = datarun.cones.centers(params.cone_roi,:);
cone_centers(:,1) = cone_centers(:,1) .* params.scale(1);
cone_centers(:,2) = cone_centers(:,2) .* params.scale(2);
cone_types = datarun.cones.types(params.cone_roi);



% Plot cones in efficient blocks, or plot each cone individually for flexibility?
if isempty(params.cone_colors) && length(params.cone_size) == 1
    % Plot in blocks
    cone_func(cone_centers(cone_types == 'L',1), cone_centers(cone_types == 'L',2), params.cone_size, 'Color', params.type_colors(1,:));
    cone_func(cone_centers(cone_types == 'M',1), cone_centers(cone_types == 'M',2), params.cone_size, 'Color', params.type_colors(2,:));
    cone_func(cone_centers(cone_types == 'S',1), cone_centers(cone_types == 'S',2), params.cone_size, 'Color', params.type_colors(3,:));
    cone_func(cone_centers(cone_types == 'U',1), cone_centers(cone_types == 'U',2), params.cone_size, 'Color', params.type_colors(4,:));
else
    for cn = 1:length(cone_centers(:,1))
        
        cone_func(cone_centers(cn,1), cone_centers(cn,2), params.cone_size, 'Color', params.cone_colors(mod(cn-1,size(params.cone_colors,1))+1,:));
    end
end


% Label cones?
if params.label
    for c = find(params.cone_roi(:)')
        text('String', num2str(c), 'Position', datarun.cones.centers(c,:), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', params.label_color, ...
            'FontSize', params.label_size);
    end
end


% Fix axes
axis equal ij;
width  = datarun.stimulus.field_width  * params.scale(1);
height = datarun.stimulus.field_height * params.scale(2);
axis([0, width, 0, height]);

if params.drawnow, drawnow; end
if nargout < 1, clear plot_axes; end



function out = marker_cones(x, y, marker_size, varargin)
out = plot(x, y, '.', 'MarkerSize', marker_size, varargin{:});

function out = circle_cones(x, y, r, varargin)
x = x(:);
y = y(:);
r = repmat(r, size(x));
out = drawCircle(x, y, r, varargin{:});