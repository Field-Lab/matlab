function [best_thresh, maximized_ui] = calc_best_ui_thresh(datarun, cell_spec, init_thresh, min_thresh, thresh_tolerance, varargin)
% CALC_BEST_UI_THRESH    Iteratively search for the contour threshold level that gives the highest uniformity index
% (see Gauthier et al. 2009, DOI:10.1371/journal.pbio.1000063)
%
% usage: [best_thresh, maximized_ui] = calc_best_ui_thresh(datarun, cell_spec, init_thresh, min_thresh, thresh_tolerance, [opts])
%
% opts: tform               []
%       simplify_contours   false
%       contours_field      'rf_contours'
%       verbose             true
%       plot                false
%
% 2010-02 phli
%

opts = inputParser;
opts.addParamValue('tform', []);
opts.addParamValue('simplify_contours', false);
opts.addParamValue('contours_field', 'rf_contours');
opts.addParamValue('verbose', true);
opts.addParamValue('plot', false);
opts.parse(varargin{:});
opts = opts.Results;

% Make function to fminsearch ui
ui_func = @(thresh) -calc_uniformity_index(datarun, cell_spec, 'thresh', thresh, 'min_thresh', min_thresh, 'contours_field', opts.contours_field, 'tform', opts.tform, 'verbose', opts.verbose, 'plot', opts.plot);

% Search for the minimum
[best_thresh, inverse_maximized_ui] = fminsearch(ui_func, init_thresh, struct('TolX', thresh_tolerance));

maximized_ui = -inverse_maximized_ui;