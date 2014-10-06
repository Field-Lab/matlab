function datarun = get_cone_rf_contours(datarun, cell_spec, varargin)
% get_cone_rf_contours     Compute contour levels for single cone rfs
%
% usage:  datarun = get_cone_rf_contours(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in
%
%   datarun.stas.contours
%
%   by default.  see options below
%
%
% optional parameters, their default values, and what they specify:
%
%
% verbose           false               show output
% fig_or_axes       []                  figure or axes to plot in. if 0, make new figure. if empty, don't plot
% foo               'bar'               how to activate the rings
%                                           'bar' - activate on site
%                                           'bore' - activate remotely
%
%
%
% 2008-11 gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('name', 'contours');
p.addParamValue('verbose', false);
p.addParamValue('fig_or_axes', []);
p.addParamValue('foo','bar', @(x)any(strcmpi(x,{'bar','bore'})));


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);

% get cell indices
cell_indices = get_cell_indices(datarun,cell_spec);

% show output
if params.verbose
    fprintf('\nComputing contours for %d cells...',length(cell_indices));
    start_time = clock; % note when it started
end


plot_cone_mosaic(datarun,'fig_or_axes',20); hold on

% loop through cells
for cc = 1:length(cell_indices)
    
    % get cell index and id
    cell_index = cell_indices(cc);
    cell_id = datarun.cell_ids(cell_index);
    
    % get cone weights
    the_weights = datarun.cones.weights(:,cell_indices(cc));
    
    % find strong weihghts
    strong_wts = the_weights> 5*robust_std(the_weights);
    
    if sum(strong_wts)<=2;continue;end
    
    % get center points
    x = datarun.cones.centers(strong_wts,1);
    y = datarun.cones.centers(strong_wts,2);
    
    % compute convex hull
    K = convhull(x,y);
    
    hull = [x(K) y(K)];
    
    plot(gca,hull(:,1),hull(:,2),'k')
end


% display how long it took
if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_time));
end





