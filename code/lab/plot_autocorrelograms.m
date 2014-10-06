function plot_autocorrelograms(datarun, cell_spec, varargin)
% plot_autocorrelograms     plot the autocorrelogram for cell(s)
%
% usage: plot_autocorrelograms(datarun, cell_spec, <params>)
%
% arguments:  datarun - datarun struct 
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    plots a figure
%
%
% optional parameters, their default values, and what they specify:
%
%
% foa               []               figure or axes to plot in. if 0, make new figure. if empty, don't plot.  if -1, plot in current.
% clear_fig         'true'           clear the figure
% line_color         'k'             line color of the autocorrelation
% normalize          false           normalize the so the peak amplitude is
%                                    one
% xlabel            []               label x axis
% ylabel            []               label y axis
% title             []               title of figure
%
% 2011-02 gdf
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('foa', []);
p.addParamValue('clear_fig', true, @islogical);
p.addParamValue('line_color', 'k');
p.addParamValue('normalize', false, @islogical)
p.addParamValue('title', []);
p.addParamValue('xlabel', []);
p.addParamValue('ylabel', []);

p.parse(varargin{:});


% BODY OF FUNCTION

% check that autocorrelation is a field in datarun
if ~isfield(datarun, 'autocorrelation')
    error('Autocorrelation is not a field, see get_autocorrelation')
end

% get cell_indices and num cells
cell_indices = get_cell_indices(datarun, cell_spec);
num_cells = length(cell_indices);

% initialize figure
plot_axes = set_up_fig_or_axes(p.Results.foa,p.Results.clear_fig);
axes(plot_axes)
hold on

for cc = 1:num_cells
    probabilities = datarun.autocorrelation{cell_indices(cc)}.probabilities;
    bins = datarun.autocorrelation{cell_indices(cc)}.bins;
    if p.Results.normalize
        plot(bins, probabilities ./ max(probabilities), 'Color', p.Results.line_color);
    else
        plot(bins, probabilities, 'Color', p.Results.line_color);
    end
end

title(p.Results.title)
xlabel(p.Results.xlabel)
ylabel(p.Results.ylabel)



    