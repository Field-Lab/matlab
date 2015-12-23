function datarun = maximize_uniformity_index(datarun, cell_spec, init_thresh, min_thresh, thresh_tolerance, varargin)
% MAXIMIZE_UNIFORMITY_INDEX    Iteratively calculate contour threshold that maximizes UI
% usage: datarun = maximize_uniformity_index(datarun, cell_spec, init_thresh, min_thresh, thresh_tolerance, [opts])
%
% OPTS are passed to CALC_BEST_UI_THRESH; see help there for syntax.
%
% See CALC_BEST_UI_THRESH for most argument info.  Puts outputs on DATARUN
% under DATARUN.STAS.MOSAICS{cell_type_num}
% 
% TODO: Merge with the mosaics struct?
%
% phli 2010-03
%

cell_type_num = get_cell_type_nums(datarun, cell_spec);

[best_thresh, maximized_ui] = calc_best_ui_thresh(datarun, cell_spec, init_thresh, min_thresh, thresh_tolerance, varargin{:});

datarun.stas.mosaics{cell_type_num}.best_thresh = best_thresh;
datarun.stas.mosaics{cell_type_num}.max_ui = maximized_ui;