function plot_rf_labels(datarun, cell_spec, varargin)
% PLOT_RF_LABELS    Plot cell id labels at rf center points
% usage: plot_rf_labels(datarun, cell_spec, [opts])
%
% Abstracted from plot_rf_summaries
%
% ToDo: * Auto adjust axes after plotting so that labels actually show up!
%       * Allow arbitrary labels
%
% phli 2010-04
%


opts = inputParser;

opts.addParamValue('tform' , []);
opts.addParamValue('input' , 'sta');
opts.addParamValue('output', 'sta');
opts.addParamValue('scale' , 1);
opts.addParamValue('skip', []);

opts.addParamValue('center_type', 'com');
opts.addParamValue('label_color', 'black');
opts.addParamValue('label_size' , 10);
opts.addParamValue('axes'       , gca);

opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.tform)
    opts.tform = coordinate_transform(datarun, opts.output, 'input', opts.input, 'scale', opts.scale);
end

% Save previous hold setting
old_hold = get(opts.axes, 'NextPlot');



cell_indices = get_cell_indices(datarun, cell_spec);
for cell_num = cell_indices(~ismember(datarun.cell_ids(cell_indices), opts.skip))
    
    center_point = rf_center(datarun, datarun.cell_ids(cell_num), opts.center_type);

    % Skip if doesn't exist
    if isempty(center_point); continue; end

    % Transform to desired coordinates
    [X, Y] = tformfwd(opts.tform, center_point(1), center_point(2));
    
    % plot cell ID there
    text('Parent', opts.axes, 'String', num2str(datarun.cell_ids(cell_num)), 'Position', [X Y] * opts.scale,  ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'color', opts.label_color, ...
        'FontSize', opts.label_size)
    hold(opts.axes, 'on');
end



% Restore previous hold setting
set(opts.axes, 'NextPlot', old_hold);
