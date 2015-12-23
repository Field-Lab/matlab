function plot_rfs(datarun, cell_specification, varargin)
% PLOT_SUMMARY     Plot the summary frame of a list of cells
%
% usage:  datarun = my_function(datarun, arg1, params)
%
% arguments:             datarun - datarun struct
%             cell_specification - see get_cell_indices for options
%                         params - struct of optional parameters (see below)
%
% outputs:     result - result of computation
%
%
% optional fields in params, their default values, and what they specify:
%
% figure        0         	figure to plot in.  if empty, creates new figure.
% overlay      	true        overlay RF summaries
% 
%
% 2009-04  gauthier
% 2010-02  phli, changed callback; no longer blocks execution to plot
%

p = inputParser;
p.addParamValue('figure', 0);
p.addParamValue('overlay', true);
p.KeepUnmatched = true;
p.parse(varargin{:});
params = p.Results;

if params.figure == 0
    params.figure = figure();
end
figure(params.figure);
params.ax = gca();

% get list of cell numbers
cell_indices = get_cell_indices(datarun, cell_specification);

% parameters of the loop
start_index=1; index_min=1; index_max=length(cell_indices);

% set up figure with scroll bar
slider = make_loop_slider_list(start_index, index_min, index_max, {@slider_plot, datarun, cell_specification, cell_indices, params, p.Unmatched});

% Initial plot call
slider_plot(slider, [], datarun, cell_specification, cell_indices, params, p.Unmatched);



function slider_plot(handle, event, datarun, cell_specification, cell_indices, params, plotargs) %#ok<INUSL> % EVENT is not used

% identify which panel this is
k = round(get(handle, 'Value'));

% get cell index and cell id
cell_index = cell_indices(k);
cell_id = datarun.cell_ids(cell_index);

plot_rf(datarun, cell_id, plotargs, 'foa', params.ax);

% show other cell IDs
if params.overlay
    plot_rf_summaries(datarun, cell_specification, struct('clear', false, 'foa', params.ax, 'label', true))
end