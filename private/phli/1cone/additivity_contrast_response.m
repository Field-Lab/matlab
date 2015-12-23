function [responses,templatex,templateyon,templateyoff] = additivity_contrast_response(dataruns, stimstructs, cellsets, varargin)
% ADDITIVITY_CONTRAST_RESPONSE
%
% Deprecated; should move towards more general methods
% BUILD_RESPONSE_TEMPLATE and ...
%

opts = inputParser;
opts.addRequired('stop');
opts.addRequired('bin_width');
opts.addParamValue('start', 0);
opts.addParamValue('template_dataruns', 1:length(dataruns));
opts.addParamValue('cr_dataruns', 1);
opts.parse(varargin{:});
opts = opts.Results;


% Build positive and negative response templates from selected dataruns
t_dataruns    = dataruns(opts.template_dataruns);
t_stimstructs = stimstructs(opts.template_dataruns);
t_cellsets    = cellsets(opts.template_dataruns);
[templatex,templateyon]  = build_contrast_response_template(t_dataruns, t_stimstructs, t_cellsets, opts.stop, opts.bin_width, 'start', opts.start, 'intensity_predicate', @(i) (i > 0));
[templatex,templateyoff] = build_contrast_response_template(t_dataruns, t_stimstructs, t_cellsets, opts.stop, opts.bin_width, 'start', opts.start, 'intensity_predicate', @(i) (i < 0));


% Calculate matched filter response values for selected dataruns
dataruns    = dataruns(opts.cr_dataruns);
stimstructs = stimstructs(opts.cr_dataruns);
cellsets    = cellsets(opts.cr_dataruns);
intensities = get_intensities(stimstructs);
for i = 1:length(dataruns)
    if iscell(dataruns) datarun = dataruns{i};
    else datarun = dataruns(i); end

    cells = cellsets{i};
    
    % New format is to leave cells as cellstruct that may have empties
    empties = cell2mat(collect(cells, @isempty));
    indices = zeros(size(empties));
    indices(~empties) = get_cell_indices(datarun, cell2mat(cells(~empties)));
    for ci = 1:length(indices)            
        cellid = indices(ci);
        if ~cellid, continue; end
        
        spikes = cell2mat(datarun.spikes(cellid));
        if iscell(stimstructs) stimstruct = stimstructs{i};
        else stimstruct = stimstructs(i); end
        intensities = get_intensities(stimstruct);

        for j = 1:length(intensities)
            intensity = intensities(j);
            matching_rgbs = stimstruct.rgbs(1:length(datarun.triggers),1) == intensity;
            if intensity >= 0, templatey = templateyon;
            else templatey = templateyoff; end
            
            for k = 1:length(stimstruct.mapindices)
                map = stimstruct.mapindices(k);
                matching_maps = stimstruct.maps(1:length(datarun.triggers)) == map;
                
                % Just copied logic for pulling PSTH spikes from raster_phli.m
                triggers = datarun.triggers(matching_rgbs & matching_maps');
                for t = 1:length(triggers);
                    time_diffs = spikes - triggers(t);
                    in_window = time_diffs >= opts.start & time_diffs <= opts.stop;
                    raster_times{t,1} = time_diffs(in_window);
                end
                d_raster_times = cell2mat(raster_times);
                
                response = histc(d_raster_times, templatex);
                
                % Subtract baseline if possible
                if opts.start < 0
                    baseline = mean(response(templatex <= 0));
                    response = response - baseline;
                end
                
                % Inner product / matched filter
                responses{i}(j,k,ci) = response(:)' * templatey;
            end
        end
    end

    if nargout == 0
        % Plot!
    end
end


if length(responses) == 1, responses = responses{1}; end
if nargout == 0, clear responses; end