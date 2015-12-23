function [plot_axes,figure_num] = set_up_fig_or_axes(fig_or_axes,clear_it)
% set_up_fig_or_axes     Clear figure or axes and return handle to the axes.
%
% usage:  [plot_axes,figure_num] = set_up_fig_or_axes(fig_or_axes,clear_it)
%
% arguments:   fig_or_axes - axis or figure handle
%                               If 0, make new figure.
%                               If -1, return current axes.
%                               If empty, does nothing and returns empty.  This is
%                                    usually the signal not to plot at all.
%                 clear_it - boolean, clear axes?
%                               does not apply if fig_or_axes = -1
%
% outputs:     plot_axes - 
%             figure_num - 
%
%  if only one argument is provided, clear_it is set to true
%
% 2008-01  gauthier
% 2009-04  gauthier, added -1 option
%


 
% if empty, return empty
if isempty(fig_or_axes)
    plot_axes = [];
    return
end

% if only one argument is provided, clear_it is set to true
if ~exist('clear_it','var')
    clear_it = true;
end

% if -1, return current axes
if fig_or_axes == -1
    plot_axes = gca;
    
    if clear_it
        cla
    end
    
    return
end



% if an integer, assume it's a figure number
if mod(fig_or_axes,1) == 0

    % if 0, make new figure
    if fig_or_axes == 0
        fig_or_axes = figure;
    end
    
    % clear the figure
    figure(fig_or_axes);
    
    % return axes
    plot_axes = gca;
    
else
    % otherwise, assume that axes were given
    
    % clear them
    if clear_it
        cla(fig_or_axes)
    end
        
    % return them
    plot_axes = fig_or_axes;
end

if nargout>1
    figure_num = get(plot_axes,'Parent');
end
