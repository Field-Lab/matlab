function plot_snl_(gen_signal,spikes,varargin)
% plot_snl_     plot SNL, optionally with fit
%
%    Three things are plotted vs. the generator signal
%           1) spike counts (plus small amount of noise for visualization)
%           2) mean spike count
%           3) fit to mean spike count
%
%
% usage:  plot_snl_(gen_signal,spikes, <params>)
%
% arguments:     gen_signal - N-length vector of generator signal values
%                    spikes - N-length vector of spike times
%                  <params> - struct or list of optional parameters (see below)
%
% outputs:     none
%
%
% optional params, their default values, and what they specify:
%
% verbose           false          	show output
% foa               []            	figure or axes to plot in. if 0, make new figure.
%                                    	if empty, don't plot.  if -1, plot in current.
% fit               []              struct specifying SNL fit
%                                       if empty, don't plot fit
% fig_title         'unknown cell'
%                                   title of figure
%
%
% 2009-09  gauthier
%
%


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('foa', []);
p.addParamValue('fit', []);
p.addParamValue('fig_title', 'unknown cell');

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



% BODY OF THE FUNCTION


% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa);

if isempty(plot_axes)
    error('plot axes not specified')
end

% turn hold on
hold(plot_axes,'on')

% plot spikes
plot(gen_signal,spikes+0.05*randn(length(spikes),1),'.','Color',[1 1 1]*.8,'Parent',plot_axes);

% plot mean
[X,Y] = curve_from_binning(gen_signal,spikes,'average_y','mean','average_x','mean');
plot(X,Y,'k-+','Parent',plot_axes)

% plot fit
if ~isempty(params.fit)
    % generate fit
    switch params.fit.type
        case 'exp'
            fit_fcn = @(x)exp(params.fit.b + params.fit.a*x);
            fit_name = sprintf('fit = exp( %0.2f + %0.2f * g )',params.fit.b,params.fit.a);
        otherwise
            error('fit type ''%s'' not recognized',params.fit.type)
    end
        
    % plot it
    fitx = linspace(min(gen_signal), max(gen_signal), 100);
    plot(fitx,fit_fcn(fitx),'r-','Parent',plot_axes)
    % append to title
    params.fig_title = [params.fig_title ', ' fit_name];
end

% add title
title(plot_axes,params.fig_title)

