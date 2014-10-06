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
% cone_roi          []              binary vector of cones in ROI
%                                       if empty, plot as normal
%                                       if non-empty, plot ROI cones in blue, other in black
% roi_highlight     true            if true, plot roi cones blue, and others black
%                                   if false, only plot roi cones and plot as specified by cone colors
%                
% type_colors       4x3 matrix of awesome cone colors
%                                   4x3 matrix of colors, each row specifying the color for L, M, S, or U cones.
%                                       only used if cone_colors is empty.
% clear             true            clear axes before plotting
% drawnow           true            draw after completing?
% cone_circles      false           uses drawCircle to draw cones as circles of radius specified by cone_radius
% cone_radius       1               radius used to draw cones as circles
%
% 2008-10 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('fig_or_axes', 0);
p.addParamValue('cone_size', 15);
p.addParamValue('bg_color',[1 1 1]*0.5);
p.addParamValue('type_colors',...
    [1 0 0; 0 1 0; 0 0 1; 0 0 0]);
    %[225/255 75/255 83/255;
    %116/255 181/255 73/255;
    %100/255 134/255 228/255;
    %182/255 158/255 61/255]);
    %[219/255 57/255 65/255;
    %98/255 167/255 55/255;
    %81/255 116/255 222/255;
    %150/255 150/255 150/255]);
p.addParamValue('cone_roi', []);
p.addParamValue('cone_colors', []);
p.addParamValue('clear', true);
p.addParamValue('drawnow', false, @islogical);
p.addParamValue('roi_highlight', true, @islogical);
p.addParamValue('scale', [1 1]);
p.addParamValue('cone_circles', false, @islogical);
p.addParamValue('cone_radius', 1, @isnumeric);
p.addParamValue('label', false);
p.addParamValue('label_color', 'black');
p.addParamValue('label_size' , 10);
p.addParamValue('cones', []);
p.addParamValue('bounds', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% BODY OF THE FUNCTION

% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes,params.clear);

% Cones provided by user, or bounds provided, or use all?
if ~isempty(params.bounds)
    inbounds = params.bounds(1) <= datarun.cones.centers(:,1) & datarun.cones.centers(:,1) <= params.bounds(2) & ...
               params.bounds(3) <= datarun.cones.centers(:,2) & datarun.cones.centers(:,2) <= params.bounds(4);
    params.cones = find(inbounds);
end
if isempty(params.cones)
    params.cones = 1:size(datarun.cones.centers,1);
end

% put variables into shorter names
cone_centers = datarun.cones.centers(params.cones,:);
cone_centers(:,1) = cone_centers(:,1) .* params.scale(1);
cone_centers(:,2) = cone_centers(:,2) .* params.scale(2);
cone_types = datarun.cones.types(params.cones);


% Set background color
if ~isempty(params.bg_color)
    set(plot_axes, 'Color', params.bg_color);
end

% Set axes, and turn hold on
axes(plot_axes);
width  = datarun.stimulus.field_width  * params.scale(1);
height = datarun.stimulus.field_height * params.scale(2);
axis([0, width, 0, height]);
axis equal;
axis ij;
set(plot_axes, 'box', 'on', 'xtick', [], 'ytick', []);
hold on

% if cone colors are not provided, then plot based on the presence of a ROI
if isempty(params.cone_colors)
    if isempty(params.cone_roi)
        
        % plot each cone
        if params.cone_circles % plot cones as circles
            for cn = 1:length(cone_centers(:,1))
                drawCircle(cone_centers(cone_types == 'L',1), cone_centers(cone_types == 'L',2), params.cone_radius, 'Color', params.type_colors(1,:))
                drawCircle(cone_centers(cone_types == 'M',1), cone_centers(cone_types == 'M',2), params.cone_radius, 'Color', params.type_colors(2,:))
                drawCircle(cone_centers(cone_types == 'S',1), cone_centers(cone_types == 'S',2), params.cone_radius, 'Color', params.type_colors(3,:))
                drawCircle(cone_centers(cone_types == 'U',1), cone_centers(cone_types == 'U',2), params.cone_radius, 'Color', params.type_colors(4,:))
            end
            
        else % plot cones as points
            plot(cone_centers(cone_types=='L',1),cone_centers(cone_types=='L',2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(1,:))
            plot(cone_centers(cone_types=='M',1),cone_centers(cone_types=='M',2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(2,:))
            plot(cone_centers(cone_types=='S',1),cone_centers(cone_types=='S',2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(3,:))
            plot(cone_centers(cone_types=='U',1),cone_centers(cone_types=='U',2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(4,:))
        end
        
        
    else
        if params.roi_highlight
            % plot cone ROI
            check = 1;
            plot(cone_centers(params.cone_roi,1),cone_centers(params.cone_roi,2),'.','MarkerSize',params.cone_size,'Color',[0 1 1])
            plot(cone_centers(~params.cone_roi,1),cone_centers(~params.cone_roi,2),'.','MarkerSize',params.cone_size,'Color',[0 0 0])
        else
            roi_indices = find(params.cone_roi);
            L_indices = intersect(roi_indices, find(cone_types == 'L'));
            M_indices = intersect(roi_indices, find(cone_types == 'M'));
            S_indices = intersect(roi_indices, find(cone_types == 'S'));
            U_indices = intersect(roi_indices, find(cone_types == 'U'));

            % plot each cone
            if params.cone_circles % plot cones as circles
                for cn = 1:length(cone_centers(:,1))
                    drawCircle(cone_centers(cone_types == 'L',1), cone_centers(cone_types == 'L',2), params.cone_radius, 'Color', params.type_colors(1,:))
                    drawCircle(cone_centers(cone_types == 'M',1), cone_centers(cone_types == 'M',2), params.cone_radius, 'Color', params.type_colors(2,:))
                    drawCircle(cone_centers(cone_types == 'S',1), cone_centers(cone_types == 'S',2), params.cone_radius, 'Color', params.type_colors(3,:))
                    drawCircle(cone_centers(cone_types == 'U',1), cone_centers(cone_types == 'U',2), params.cone_radius, 'Color', params.type_colors(4,:))
                end
            else
                plot(cone_centers(L_indices,1),cone_centers(L_indices,2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(1,:))           
                plot(cone_centers(M_indices,1),cone_centers(M_indices,2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(2,:))           
                plot(cone_centers(S_indices,1),cone_centers(S_indices,2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(3,:))           
                plot(cone_centers(U_indices,1),cone_centers(U_indices,2),'.','MarkerSize',params.cone_size,'Color',params.type_colors(4,:))           
            end
      
        end
    end
else
    % plot each cone
    if params.cone_circles % plot cones as circles
        for cn = 1:length(cone_centers(:,1))
            drawCircle(cone_centers(cn,1), cone_centers(cn,2), params.cone_radius, 'Color', params.cone_colors(mod(cn-1,size(params.cone_colors,1))+1,:))
        end
    end
   % if params.cone_colors is non-empty, plot each cone a different color 
   for cc = 1:size(cone_centers,1)
        plot(cone_centers(cc,1),cone_centers(cc,2),'.','MarkerSize',params.cone_size,...
            'Color',params.cone_colors(mod(cc-1,size(params.cone_colors,1))+1,:))
   end
end


% draw, if desired
if params.drawnow
    drawnow
end

if params.label
    for c = params.cones
        text('String', num2str(c), 'Position', datarun.cones.centers(c,:), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', params.label_color, ...
            'FontSize', params.label_size);
    end
    
    % draw, if desired
    if params.drawnow
        drawnow
    end
end

% don't spew the plot axes if they're not requested
if nargout < 1
    clear plot_axes
end
