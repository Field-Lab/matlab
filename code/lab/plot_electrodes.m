function h = plot_electrodes(positions, varargin)
% PLOT_ELECTRODES    Plot electrode positions, especially color groups
% usage: plot_electrodes(positions, [groups, plot_ei_options])
%
% Typical usage would be: plot_electrodes(positions, {[1:100] [101:200]}),
% which would plot all electrodes as filled circles, with the first 100
% electrodes one color, the second 100 another color, and the rest grey.
%
% Additional options are passed on to plot_ei_.  Useful example: 
%   'label', true, 'scale', 0.1
%
% 2010-06 phli
%

numelec = size(positions, 1);
fake_ei = ones(numelec, 1);
frame = 1;

if ~isempty(varargin) && iscell(varargin{1})
    % Assume electrode groupings were given
    groups = varargin{1};
    varargin = varargin(2:end);
    
    % Get RGB triples for electrode colors according to groupings, for
    % colors not included in a group, default to light grey.
    colors = grouped_colors(groups, numelec, [0.9 0.9 0.9]);

    h = plot_ei_(fake_ei, positions, frame, 'elec_colors', colors, 'max_scale', 0.8, varargin{:});
else
    % Just pass everything to plot_ei_
    h = plot_ei_(fake_ei, positions, frame, varargin{:});
end