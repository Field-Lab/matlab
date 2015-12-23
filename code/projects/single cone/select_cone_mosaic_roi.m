function datarun = select_cone_mosaic_roi(datarun, varargin)
% select_cone_mosaic_roi     select a region of interest in the cones
%
% usage:  datarun = select_cone_mosaic_roi(datarun, <params>)
%
% arguments:  datarun - datarun struct
%            varargin - struct or list of optional parameters (see below)
%
% outputs:     datarun - datarun struct, with datarun.cones.roi set to the selected ROI
%
%
% optional params, their default values, and what they specify:
%
% fig_or_axes       0        	figure or axes to plot in. if 0, make new figure.
% bg                []          image to plot in the background
%                                   if empty, a plot of the cones is used
%
%
% 2009-01  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('fig_or_axes', 0, @(x)~isempty(x));
p.addParamValue('bg', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.fig_or_axes);
axes(plot_axes)

% plot background image
if ~isempty(params.bg)
    imagesc(params.bg);axis image
else
    plot_cone_mosaic(datarun,'fig_or_axes',plot_axes)
end

% prompt user to select region
[junk,xi,yi] = roipoly;

% identify which points are inside the polygon
cone_roi = inpolygon(datarun.cones.centers(:,1),datarun.cones.centers(:,2),xi,yi);

% save to datarun
datarun.cones.roi = cone_roi;

