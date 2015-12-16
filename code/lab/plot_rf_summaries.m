function plot_axes = plot_rf_summaries(datarun, cell_specification, varargin)
% plot_rf_summaries     Plot collection of RFs in a variety of ways
%
%
%  NOTE: THIS FUNCTION IS STILL UNDER CONSTRUCTION -- NOT ALL FEATURES HAVE BEEN IMPLEMENTED
%
%
%
% usage:  plot_axes = plot_rf_summaries(datarun,cell_specification, varargin)
%
% arguments:             datarun - datarun struct
%             cell_specification - see get_cell_indices for options
%                       varargin - struct or list of optional parameters (see below)
%
% outputs:     plot_axes - handle to the figure or axes of the plot
%
%
% optional fields in params, their default values, and what they specify:
%
%   general parameters
%
% foa               -1          figure or axes to plot in, see set_up_fig_or_axes for options
% bg_color          []         title(['trial order, based on distance metric (bottom to top), normalization type = ' normType]) matlab color spec, color of background.
%                                   if empty, no bg is plotted
% clear             true        clear previous contents of axes
% skip              []          list of cell_ids to skip plotting; useful
%                                   in conjunction with a {1} style cell_specification
% array             false       plot array outline.  this requires datarun.piece.corners to be non-empty.
% array_label       false       plot electrode IDs of corner electrodes
% array_color       [1 0 0]     specifies color of array outline
%
% ***** STA contour parameters *****
%
% plot_contours         false           plot contours
% contours_field       	'rf_contours' 	name of field in dataset.stas in which contours are stored
% contour_colors       	'kr'            colors in which to plot the contours
%                                           for numbers, use this form:  [1 0 0; 0 1 0; 0 0 1];
% contour_fill          'r'             color to fill contours
% contour_alpha         0.5             transparency for contour fills
% contour_widths        10              width to plot contours
% rf_tform              []              transformation to apply to each rf; this is passed to
%                                           tform_rf_contours; refer there for options
%
%
% ***** center point parameters *****
%
% plot_centers      false       plot center points
% center_type       []          type of center point to plot
%                                   see rf_center for options
% label             false       label each center point with the cell ID
% label_color       'k'         color to plot the IDs.  can be standard matlab color spec,e.g. [.5 .5 .5]
% label_size        10          point size to plot the cell IDs
%
%
% ***** gaussian STA fit parameters *****
%
% plot_fits         false       plot gaussian STA fits
% fit_color         'k'         color to plot the fits.  can be standard matlab color spec,e.g. [.5 .5 .5]
% fit_width         1           thickness of plotted fits
% fill_color        'none'      color to fill fits with.  can be standard matlab color spec
%
%
% ***** coordinate parameters *****
%
% coordinates       'sta'       what coordinate space to plot in
%                                   'sta' - coordinates of the STA
%                                   'monitor' - stimulus monitor coordinates
%                                   'rect' - coordinates of provided rectangle
% rect              []          rectangle in which to plot
%                                   should be a cell array, {coord, rect_spec}
%                                   coord can be 'sta' or 'monitor', telling which coordinates the rect_spec are in
%                                   rect_spec should be a vector: [x_center y_center width height angle]
%
%
%
%
%   NOTE: if plot_centers, plot_fits, and plot_contours are all false,
%       the highest level description will be plotted.
%
%
%
% This function plots RFs in all possible ways, from simple center points
% to multiple levels of contours.  The user can specify which type of plot
% to make, along with various parameters (e.g. color) for each plot.  Each
% option can be specified independently, allowing for multiple types of
% information to be shown at once (e.g. contours and gaussian fit outlines).
% If no specification is provided, the function defaults to plotting the most
% descriptive information available.
%
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa', -1);
p.addParamValue('bg_color', []);
p.addParamValue('clear', true);
p.addParamValue('skip', []);
p.addParamValue('scale', 1, @isnumeric);

p.addParamValue('plot_contours', false);
p.addParamValue('contours_field', 'rf_contours');
p.addParamValue('contour_colors', []);
p.addParamValue('contour_widths', 1);
p.addParamValue('contour_fill', 'r');
p.addParamValue('contour_alpha', 0.5);

p.addParamValue('plot_centers', false);
p.addParamValue('center_type', []);
p.addParamValue('label', false);
p.addParamValue('label_color', 'k');
p.addParamValue('label_size', 10);

p.addParamValue('plot_fits', false);
p.addParamValue('fit_color', 'k');
p.addParamValue('fit_width', 1);
p.addParamValue('fill_color', 'none');

p.addParamValue('coordinates', 'sta');
p.addParamValue('rect', []);

p.addParamValue('mirror', false); % Backwards compatibility...
p.addParamValue('rf_tform', []);

p.addParamValue('array', false);
p.addParamValue('array_label', false);
p.addParamValue('array_color', [1 0 0]);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% Backwards compatibility...
if params.mirror && isempty(params.rf_tform)
    params.rf_tform = 'xreflect';
end


% if all plots are false, search for the most sophisticated kind of information available
if ~params.plot_centers && ~params.plot_fits && ~params.plot_contours

    % plot fits if they exist...
    if isfield(datarun,'stas') && isfield(datarun.stas,'fits')
        params.plot_fits = true;
    else
        % otherwise, just center points
        params.plot_centers = true;
        params.label = true;
    end
end



% SET UP AXES FOR PLOTTING

% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa,params.clear);
axes(plot_axes)


% ensure hold is on
hold on

coord_tform = coordinate_transform(datarun,params.coordinates);

% set output coordinates
switch params.coordinates
    case 'sta'
        bg_size = [datarun.stimulus.field_width datarun.stimulus.field_height] * params.scale;
    case 'monitor'
        bg_size = [datarun.stimulus.monitor_x datarun.stimulus.monitor_y] * params.scale;
    case 'rect'
        error('not yet implemented')
end 

% compute axes size
%pixel_height = datarun.stimulus.field_height * datarun.stimulus.stixel_height;
%pixel_width = datarun.stimulus.field_width * datarun.stimulus.stixel_width;

% set scale factor
% switch 2
%     case 1 % plot in pixel space
%         scale_x = datarun.stimulus.stixel_width;
%         scale_y = datarun.stimulus.stixel_height;
%     case 2 % plot in stixel space
%         scale_x = 1;
%         scale_y = 1;
% end
        

% plot color in background, if desired
if ~isempty(params.bg_color)
    %image_data = repmat(reshape(params.bg_color,1,1,3),datarun.stimulus.field_height,datarun.stimulus.field_width);
    image_data = repmat(reshape(params.bg_color,1,1,3),...  color
        bg_size(2),...    y size
        bg_size(1)); %    x size
    image('Parent',plot_axes,'CData',image_data)
end

% set axes size
%set(plot_axes,'XLim',[1 datarun.stimulus.field_width],'YLim',[1 datarun.stimulus.field_height],'DataAspectRatio',[1 1 1])
set(plot_axes,'XLim',[1 bg_size(1)],'YLim',[1 bg_size(2)],'DataAspectRatio',[1 1 1])
axis ij



% get list of cells to plot
[cell_indices, cell_type_name, cell_type_num] = get_cell_indices(datarun, cell_specification);



if params.plot_contours
    % Apply transforms; OVERWRITES DATARUN!!!  This must be changed if
    % datarun is to be returned!!
    if params.rf_tform
        datarun = tform_rf_contours(datarun, cell_specification, params.rf_tform, 'contours_field', params.contours_field, 'output_field', params.contours_field);
    end
    datarun = tform_rf_contours(datarun, cell_specification, coord_tform, 'contours_field', params.contours_field, 'output_field', params.contours_field);

    
    % loop through cells
    for cell_index = cell_indices(~ismember(datarun.cell_ids(cell_indices), params.skip))
        
        % get contours
        cnt_poly = datarun.stas.(params.contours_field){cell_index};
        num_contours = length(cnt_poly);
                
        % if no contour colors were specified, generate a gradient of colors
        if isempty(params.contour_colors)
            grad_min = [1 1 1]*0.4;
            grad_max = [1 0 0];
            params.contour_colors = dotted_line([grad_min;grad_max], num_contours);
        end

        % if a list of colors was specified, ensure they are a column vector
        if ischar(params.contour_colors) && size(params.contour_colors,1) == 1
            params.contour_colors = params.contour_colors';
        end

        % set up list of colors as a cell array
        for cc = 1:size(params.contour_colors, 1)
            contour_colors{cc} = params.contour_colors(cc,:);
        end

        % go through each threshold
        for tt = 1:length(cnt_poly)
            % get color
            col = contour_colors{mod(tt-1,cc) + 1};
            
            % get width
            wid = params.contour_widths(mod(tt-1,length(params.contour_widths)) + 1);

            % Plot it
            plot_polygon_struct(cnt_poly{tt}, 'facecolor', params.contour_fill, 'bgcolor', 'w', 'alpha', params.contour_alpha, 'linecolor', col, 'LineWidth', wid);
        end

    end
    
end



% PLOT FITS

if params.plot_fits
    
    % 2D Gaussian fits
    
    % go through each cell
    for cell_index = cell_indices(~ismember(datarun.cell_ids(cell_indices), params.skip))

        % get the fit
        the_fit = datarun.stas.fits{cell_index};

        % skip if doesn't exist
        if isempty(the_fit);continue;end

        % get center
        ctr = the_fit.mean * params.scale;

        % get center radius
        rad = the_fit.sd * params.scale * 1.5;
 
        % get points of an ellipse with these parameters
        [X,Y] = drawEllipse([ctr rad the_fit.angle]);
        
        % if any are NaN, skip
        if any(isnan([X Y]));continue;end
        
        % transform to desired coordinates
        [X, Y] = tformfwd(coord_tform, X, Y);

        % plot the points and fill
        
        if ~strcmpi(params.fill_color, 'none')
            fill(X,Y,params.fill_color)
        end
        plot(X,Y,'Color',params.fit_color, 'LineWidth', params.fit_width)
    end
    
    % cone fits
    if 0
        % plot gaussian fit to each cell
        for cell_index = cell_indices(~ismember(datarun.cell_ids(cell_indices), params.skip))
            % get the fit
            the_fit = datarun.cones.rf_fits{cell_index};

            % skip if doesn't exist
            if isempty(the_fit);continue;end

            % get center
            ctr = (the_fit.center) .* [scale_x scale_y];

            % get center radius
            rad = the_fit.center_radius * scale_x;

            % get points of an ellipse with these parameters
            [X,Y] = drawEllipse([ctr rad rad 0]);

            % transform to desired coordinates
            [X, Y] = tformfwd(coord_tform, X, Y);

            % plot the points
            plot(X,Y,'b')
        end
    end
    
end



if params.plot_centers    
end


if params.label
    plot_rf_labels(datarun, cell_specification, 'tform', coord_tform,...
                'center_type', params.center_type, 'skip', params.skip,...
                'label_color', params.label_color, 'label_size',...
                 params.label_size, 'axes', plot_axes, 'scale', params.scale);
end


% PLOT ARRAY
if params.array
    switch params.coordinates
        case 'monitor'
            plot_array(datarun,params.array_color,[],params.array_label)
        case 'sta'
            T = coordinate_transform(datarun,'monitor');
            plot_array(datarun,params.array_color,fliptform(T),params.array_label);
        otherwise
            error('coordinate system ''%s'' not recognized',params.coordinates)
    end
        
end



% For use with SAVE_FIGURES
if params.plot_contours
    prefix = 'rfcont_';
else
    prefix = 'rfsummaries_';
end
savename = [prefix datarun.names.short_name '_ct' num2str(cell_type_num)];
set_savename(savename);


% don't return handle if not requested
if nargout == 0
    clear plot_axes
end
