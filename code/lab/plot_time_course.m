function plot_time_course(datarun,cell_id, varargin)
% plot_time_course     plot RF, marks, and time course 
%
% usage:  plot_time_course(datarun,cell_id, varargin)
%
% arguments:     datarun - datarun struct
%                cell_id - which cell
%               varargin - struct or list of optional parameters (see below)
%
% optional params, their default values, and what they specify:
%
% figure            0              	figure to plot in. if 0, make new figure. if -1, use current figure
% clear             true            clear figure before plotting
%
%
% 2009-06 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure', 0);
p.addParamValue('clear', true, @islogical);
p.addParamValue('normalize', false, @islogical);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;


% BODY OF THE FUNCTION

% set up plot axes
set_up_fig_or_axes(params.figure);

% clear if desired
if params.clear; clf; end

% get cell index
cell_index = get_cell_indices(datarun,cell_id);

% get cell type name
cell_type_number = find_cell_types(datarun, cell_id);
if cell_type_number > 0
    cell_type = datarun.cell_types{cell_type_number}.name;
else
    cell_type = 'unknown type';
end

% make title
title(sprintf('cell id %d (%s)',cell_id,cell_type))

% plot time_course, if it exists
time_course = datarun.stas.time_courses{cell_index};
if ~isempty(time_course)
    if params.normalize
        time_course = time_course ./ norm(time_course);
    end
    plot_time_course_(time_course)
    xlabel('frame number')
    ylabel('contrast')
end

