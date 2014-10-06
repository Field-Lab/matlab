function plot_axes = plot_rf_cone_profile(datarun,cell_id, varargin)
% plot_rf_cone_profile     Plot cone weights as a function of distance from the center of the R
%
% usage:  plot_axes = plot_rf_cone_profile(datarun,cell_id, varargin)
%
% arguments:     datarun - datarun struct
%                cell_id - which cell to plot
%               varargin - struct or list of optional parameters (see below)
%
% outputs:     plot_axes - where it was plotted
%
%
% optional params, their default values, and what they specify:
%
% by_cone_type      true          	distinguish L, M, S, and U contributions
% fig_or_axes       []            	figure or axes to plot in. if 0, make new figure. if empty, don't plot
%
%
% parameters passed on to other functions. if not specified by the user, these parameters are not passed.
%
%       passed to 'get_rf_cone_profiles'
%
%         	center_type       	center point of the RF
%           radius              radius in which to plot points
%
%
%       passed to 'curve_from_binning'
%
%           bin_edges           how to bin cone weights
%           average_y           how to combine cone weights in a bin
%
%
%
% gauthier 2008-10
%
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('fig_or_axes', 0);
p.addParamValue('by_cone_type', true);

% parameters to pass
%       to get_rf_cone_profiles
p.addParamValue('center_type', 'default value');
p.addParamValue('bin_edges', 'default value');
%       to curve_from_binning
p.addParamValue('average_y', 'default value');
p.addParamValue('radius', 'default value');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;

% generate structs to pass on
profiles_params = make_struct_to_pass(p.Results,{'center_type','center_type','radius','radius'});
curve_params = make_struct_to_pass(p.Results,{'bin_edges','bin_edges','average_y','average_y'});



% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% get cone weights as a function of distance
[all_x, all_y] = get_rf_cone_profiles(datarun,cell_id,'by_cone_type',params.by_cone_type,profiles_params);

% interpolate a curve
[average_x, average_y] = curve_from_binning(all_x,all_y,'bin_edges',0:2:100,curve_params);


% plot the curve(s)

% add a zero line
%plot([min(params.bin_edges) max(params.bin_edges)],[0 0],'-','Color',[.7 .7 .7])

% if no cone colors were given
if ~params.by_cone_type
    % plot everything in black
    plot(average_x,average_y,'k')
else
    % otherwise plot by color
    hold on
    colors = ['rgbk']';
    for col = 1:4
        plot(average_x{col},average_y{col},['-' colors(col)])
    end
end

cell_type_name = datarun.cell_types{find_cell_types(datarun,cell_id)}.name;
title(sprintf('cell id %d  (%s)',cell_id,cell_type_name))


% don't return if not requested
if nargout<1
    clear plot_axes
end


drawnow


