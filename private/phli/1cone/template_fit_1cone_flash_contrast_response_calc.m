function datarun = template_fit_1cone_flash_contrast_response_calc(datarun, rgcfieldname, templaterun, templateindex, varargin)

% Clipping half wave rectification
clipfitfunc = @(s, timecourse_template) (max([zeros(1,length(timecourse_template)); timecourse_template.*s]));


opts = inputParser();
opts.addParamValue('hist_start', -0.25);
opts.addParamValue('hist_end', 0.75);
opts.addParamValue('fitfunc', clipfitfunc);
opts.parse(varargin{:});
opts = opts.Results;


stimstruct = datarun.stimulus;
rgcids = datarun.(rgcfieldname);

% Futz to get hist_bins centered around zero
interval = templaterun.stimulus.interval;
monitor_refresh = templaterun.stimulus.refresh_period / interval / 1000;
hist_bins = [fliplr(0:-monitor_refresh:opts.hist_start) monitor_refresh:monitor_refresh:opts.hist_end];

% Build response rasters and PSTHs
triggers = datarun.triggers(1:2:end);
for cid = rgcids;
    if isempty(cid{1}), continue; end
    cnum = get_cell_indices(datarun, cid{1});
    
    for i = 1:length(stimstruct.uintensities)
        intensities = stimstruct.uintensities{i};
        for j = 1:length(intensities)
            intensity = intensities(j);
            if intensity == 0
                trials = stimstruct.urgb.blanks;
            else
                trials = stimstruct.urgb.singles & stimstruct.urgb.intensities(i,:) == intensity;
            end
            trials = cell2mat(stimstruct.urgbi(trials)');
            
            % If run was interrupted, we won't have all the trials
            trials = trials(trials <= length(triggers));
            
            [res, ~, histx, hist] = rasterphli(datarun, cid{1}, triggers(trials), 'hist', true, 'hist_bins', hist_bins);
            datarun.irasters{cnum}{i}(j).res   = res;
            datarun.irasters{cnum}{i}(j).histx = histx;
            datarun.irasters{cnum}{i}(j).hist  = hist;
        end
    end
end


% Build timecourse response templates from earlier white noise run
% x scaled to fundamental monitor rate
template = fliplr(templaterun.vision.timecourse_templates{templateindex}');
[timecourse_template, timecourse_templatex] = flash_timecourse_template(template, templaterun, stimstruct);

% Calculate matched filter response (with half wave rectification)
for cid = rgcids
    if isempty(cid{1}), continue; end
    cnum = get_cell_indices(datarun, cid{1});
    
    irasters = datarun.irasters{cnum};
    for cone = 1:length(irasters)
        crasters = irasters{cone};
        intensities = stimstruct.uintensities{cone};
        
        for j = 1:length(intensities)
            intensity = intensities(j);
            
            rasterstruct = crasters(j);
            psth = rasterstruct.hist;
            psthx = rasterstruct.histx;
            
            baseline = mean(psth(psthx < 0));
            
            % Line up x axes (should already be same scale just shifted)
            psth = interp1(psthx, psth, timecourse_templatex, 'linear', NaN);
            temp = timecourse_template(~isnan(psth));
            tempx = timecourse_templatex(~isnan(psth));
            psth = psth(~isnan(psth));
            
            % Fit
            polarity = intensity / abs(intensity);
            if intensity == 0, polarity = -1; end
            [clipfit, resnorm] = lsqcurvefit(clipfitfunc, 1, temp.*polarity + baseline, psth);
            
            % Plot
            figure(cid{1}+10000);
            ax(j) = sanesubplot(length(irasters), length(intensities).*2, [cone j]);
            plot(tempx, [psth; temp.*polarity.*clipfit]);
            title(['b' num2str(baseline)]);
            
            datarun.irasters{cnum}{cone}(j).clipfit = clipfit .* polarity .* -1;
            datarun.irasters{cnum}{cone}(j).fitresnorm = resnorm;
        end
        
        setaxesy(ax);
        clear ax
    end
end

% Calc and plot xscale cone weight fits
crsx = cell2mat(stimstruct.uintensities');
for cid = rgcids
    if isempty(cid{1}), continue; end
    cnum = get_cell_indices(datarun, cid{1});
    
    irasters = datarun.irasters{cnum};
    crs = cell2mat(cellfun(@(s) ([s.clipfit]), irasters, 'UniformOutput', false)');
    
    figure(cid{1}+10000);
    sanesubplot(1, 2, [1 2]);
    p = normcdfxscale(crs, crsx, 'plot', true, 'rectfit', false);
end