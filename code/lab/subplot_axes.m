function plot_axes = subplot_axes(fig,rect,pad_x,pad_y,x,y)
% generates coordinates for axes which are located in the specified
% rectangle of the specified figure.  can return either the axis boundaries
% or actually make the axes and return the handles to them.
%
% 
% usage:  axes = subplot_axes(fig,rect,pad_x,pad_y,x,y)
%
% arguments:      fig - which figure to use (can be 0, see below)
%                rect - subset of the figure in which to make a grid of subplots
%                           standard matlab format for specifying a rectangle:
%                           [y_start x_start height width]
%        pad_x, pad_y - how much space to put between subplots
%                           units are portion of the subplot axis size
%                           e.g. pax_x = 0.5 means the margins are half as wide
%                                   as the plots themselves
%                 x,y - number of subplots across and down
%
%
% if fig is 0, it returns a cell array of axis positions
% otherwise, it returns a cell array of axis handles in the specified
% figure


% DOCUMENTATION SHOULD BE IMPROVED


% compute size of each axis
% a_width = (rect(3) - (x-1)*pad_x ) / x;
% a_height = (rect(4) - (y-1)*pad_y ) / y;
a_width = rect(3) / ( x + (x)*pad_x );
a_height = rect(4) / ( y + (y)*pad_y );

% if 0 or NaN input, return empty
if x ==0 || y == 0 || isnan(x) || isnan(y)
    plot_axes=cell(1);
    return
end

% initialize plot_axes
plot_axes = cell(x*y,1);


aa = 1;
for yy = y:-1:1
    for xx = 1:x
        begin_x = rect(1) + pad_x*a_width/2 + (xx-1) * (a_width*(1+pad_x));
        begin_y = rect(2) + pad_y*a_height/2 + (yy-1) * (a_height*(1+pad_y));
        plot_axes{aa} = [begin_x begin_y a_width a_height];
        aa = aa+1;
    end
end

if 1
    figure(fig);
    for aa =1:length(plot_axes)
        plot_axes{aa} = subplot('Position',plot_axes{aa});
    end
end

