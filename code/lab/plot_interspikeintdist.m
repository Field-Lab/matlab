function plot_interspikeintdist(datarun, cell_spec, varargin)
% plot_interspikeintdist     plot the interspikeinterval probability distribution for cell(s)
%
% usage: plot_interspikeintdist(datarun, cell_spec, <params>)
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
% line_color         'k'             line color of the isi distribution
% normalize          false           normalize the so the peak amplitude is
%                                    one
% xlabel            []               label x axis
% ylabel            []               label y axis
% title             []               title of figure
%
% 2011-02 gdf
% 2014-01 sravi modified initial autocorrelation function to isi


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

% check that interspikeinterval is a field in datarun
if ~isfield(datarun, 'interspikeinterval')
    error('Interspikeinterval is not a field, see get_interspikeinterval')
end

% get cell_indices and num cells
cell_indices = get_cell_indices(datarun, cell_spec);
num_cells = length(cell_indices);

% initialize figure
plot_axes = set_up_fig_or_axes(p.Results.foa,p.Results.clear_fig);
axes(plot_axes)
hold on

for cc = 1:num_cells
    probabilities = datarun.interspikeinterval{cell_indices(cc)}.probabilities;
    bins = datarun.interspikeinterval{cell_indices(cc)}.bins;
    if p.Results.normalize
        plot(bins, probabilities ./ max(probabilities), 'Color', p.Results.line_color);
    else
        plot(bins, probabilities, 'Color', p.Results.line_color);
    end
end

title(p.Results.title)
xlabel(p.Results.xlabel)
ylabel(p.Results.ylabel)