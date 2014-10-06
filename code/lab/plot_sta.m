function plot_sta(datarun,cell_id,varargin)
% plot_sta     plot an STA in an interactive slider window
%
% usage:  plot_sta(datarun, cell_id, varargin)
%
% arguments:  datarun - datarun struct with field datarun.stas.stas
%             cell_id - which cell id to plot
%            varargin - struct or list of optional parameters (see below)
%
%
% optional params, their default values, and what they specify:
%
% figure          	0            	figure or axes where to plot it. 0 for new figure, -1 in current axes
% array             false           plot array outline?  taken from datarun.piece.corners
%
%
% 2009-06  gauthier
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('figure', 0);%, @(x)~isempty(x)&&round(x)==x);
p.addParamValue('array', false);


% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% get array corners, if desired
if params.array
    if isfield(datarun,'piece') && isfield(datarun.piece,'corners')
        T = coordinate_transform(datarun,'monitor');
        corners = tforminv(T,datarun.piece.corners);
    else
        corners = [];
        warning('could not plot array outline because datarun.piece.corners does not exist.  see compute_monitor_to_array_transformation.')
    end
else
    corners = [];
end

% set title

% cell ID
title_text = sprintf('cell id %d',cell_id);

% cell type
ct = find_cell_types(datarun,datarun.cell_ids);
if ct; title_text = [title_text ' (' datarun.cell_types{ct(1)}.name ')']; end

% spike count
if isfield(datarun,'spikes'); title_text = sprintf('%s, %d spikes',title_text,length(datarun.spikes{get_cell_indices(datarun,cell_id)}));end

% plot sta with title
plot_sta_(get_sta(datarun,cell_id),'figure',params.figure,'overlay',corners,...
    'prefix',[title_text ', '])


