function [templatex,templatey] = build_contrast_response_template(dataruns, stimstructs, cellsets, varargin)
% BUILD_CONTRAST_RESPONSE_TEMPLATE
% 
% DEPRECATED; should move towards using more general BUILD_RESPONSE_TEMPLATE
%
% usage: [templatex,templatey] = build_contrast_response_template(dataruns, stimstructs, cellsets, varargin)
%
% 2011 phli
%

opts = inputParser;
opts.addRequired('stop');
opts.addRequired('bin_width');
opts.addParamValue('intensities', []);
opts.addParamValue('start', 0);
opts.KeepUnmatched = true;
opts.parse(varargin{:});
unmatched = opts.Unmatched;
opts = opts.Results;

if ~isempty(opts.intensities)
    intensities = opts.intensities;
else
    giopts = keepfields(unmatched, 'intensity_predicate', 'with_blanks', '-force');
    intensities = get_intensities(stimstructs, giopts);
end


for i = 1:length(dataruns)
    if iscell(dataruns), datarun = dataruns{i};
    else datarun = dataruns(i); end
    
    if iscell(stimstructs), stimstruct = stimstructs{i};
    else stimstruct = stimstructs(i); end
    
    % Collect triggers for all intensities
    clear triggers; % Get rid of version from previous loop
    for j = 1:length(intensities)
        triggers{j,1} = datarun.triggers(stimstruct.rgbs(1:length(datarun.triggers),1) == intensities(j));
    end
    triggers = cell2mat(triggers);
    
    if iscell(cellsets)
        cells = cellsets{i};
        if iscell(cells)
            cells = cell2mat(select(cells, @(e)(~isempty(e))));
        end
    elseif length(dataruns) == 1, 
        cells = cellsets;
    else
        error('Not enough cellsets provided'); 
    end
    
    indices = get_cell_indices(datarun, cells);
    spikes = cell2mat(datarun.spikes(indices));
    
    % Get spike times in window relative to triggers
    for j = 1:length(triggers);
        time_diffs = spikes - triggers(j);
        in_window = time_diffs >= opts.start & time_diffs <= opts.stop;
        raster_times{j,1} = time_diffs(in_window);
    end
    d_raster_times{i,1} = cell2mat(raster_times);
end
d_raster_times = cell2mat(d_raster_times);


templatex = opts.start:opts.bin_width:opts.stop;
templatey = histc(d_raster_times, templatex);

% Subtract baseline if possible
if opts.start < 0
    baseline = mean(templatey(templatex <= 0));
    templatey = templatey - baseline;
end

% Normalize template
templatey = templatey ./ sqrt(templatey' * templatey);


if nargout == 0
    plot(templatex, templatey);
    clear templatex templatey;
end