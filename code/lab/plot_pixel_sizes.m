function plot_pixel_sizes(datarun, varargin)
% plot_pixel_sizes     plot pixel size at various points on the monitor
%
% usage:  result = my_function(arg1, <params>)
%
% arguments:     datarun - datarun struct
%                <params> - struct or list of optional parameters (see below)
%
%
%
% optional params, their default values, and what they specify:
%
% figure        0        	figure to plot in. if 0, make new figure.
%                            	if -1, plot in current.
%
%
% 2010-06  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure', 0);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up figure
set_up_fig_or_axes(params.figure);
plot_fig = gcf;
figure(plot_fig);clf;
xlim([1 datarun.stimulus.monitor_x])
ylim([1 datarun.stimulus.monitor_y])
axis image
set(gca,'ydir','reverse')
hold on
title('pixel height and width at various points on the monitor')

% plot the array
plot_array(datarun,[1 1 1]*0.7);

% choose points to plot
comparison_points = [datarun.piece.corners(2:end,:); mean(datarun.piece.corners(2:end,:),1)];

% for each point...
for cc = 1:size(comparison_points,1)
    
    % get pixel size
    [junk,pixel_height,pixel_width] = compute_pixel_size(datarun,'monitor_location',comparison_points(cc,:));
    
    % plot a point
    plot(comparison_points(cc,1),comparison_points(cc,2),'.r')
    
    % show the height and width
    text(comparison_points(cc,1),comparison_points(cc,2),...
        sprintf('h: %0.2f\nw: %0.2f',pixel_height,pixel_width),...
        'Color',[0 0 0],'FontSize',12,'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    
end
