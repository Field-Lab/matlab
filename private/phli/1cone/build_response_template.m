function [templatex,templatey] = build_response_template(dataruns, triggersets, rgcsets, varargin)
% BUILD_RESPONSE_TEMPLATE
% usage: [templatex,templatey] = build_response_template(dataruns, stimstructs, opts)
%
% OPTS are mostly options for RASTERPHLI, but a few are used here.
%
% Triggers can be selected for the desired trials before being passed here,
% or else they can be selected using OPTS.STIMSTRUCTS and OPTS.PREDICATES.
% If these are to be used then the triggers passed in should be equal in
% number to the number of trials.  OPTS.STIMSTRUCTS is a struct or cell 
% array with whatever information about each trial is needed.
% OPTS.PREDICATES is a cell array of function handles, one for each
% datarun.  The PREDICATE is run on each element of STIMSTRUCTS and should
% return true or false to indicate which trials are selected.  PREDICATES
% can also be a single function handle to be used for every datarun.
%
% 2012-02 phli, derived from BUILD_CONTRAST_RESPONSE_TEMPLATE but designed to be more general
%

opts = inputParser;
opts.addParamValue('stop', []);
opts.addParamValue('hist_bin', 0.01);
opts.addParamValue('start', -0.1);
opts.addParamValue('stimstructs', []);
opts.addParamValue('predicates', {});
opts.KeepUnmatched = true;
opts.parse(varargin{:});
rasteropts = opts.Unmatched;
opts = opts.Results;


rasteropts.start    = opts.start;
rasteropts.stop     = opts.stop;
tasteropts.hist_bin = opts.hist_bin;


d_raster_times = cell(length(dataruns),1);
for i = 1:length(dataruns)
    if iscell(dataruns), datarun = dataruns{i};
    else datarun = dataruns(i); end    
    triggers = triggersets{i};
    rgcs = rgcsets{i};
    
    % If RGCS are still in cell array form, convert them to num array; may need to have blanks removed.
    if iscell(rgcs), rgcs = cell2mat(rgcs); end

    % Use stimstructs and predicates to select from among the triggers given (or else just use all triggers given)
    if ~isempty(opts.stimstructs)
        if iscell(opts.stimstructs), stimstruct = opts.stimstructs{i};
        else stimstruct = opts.stimstructs(i); end
        
        if iscell(opts.predicates), predicate = opts.predicates{i};
        else predicate = opts.predicates; end
        
        selected_trials = cell2mat(collect(stimstruct, predicate));
        selected_trials = selected_trials(1:length(triggers)); % Drop selected trials that were not actually run (determined from having more selected trials than triggers)
        triggers = triggers(selected_trials);
    end
    
    res = rasterphli(datarun, rgcs, triggers, rasteropts, 'plot', false);
    if ~isempty(res)
        d_raster_times{i,1} = res(:,1);
    end
end
d_raster_times = cell2mat(d_raster_times);


templatex = opts.start:opts.hist_bin:opts.stop;
templatey = hist(d_raster_times, templatex);

% FIXME: Ugly; can handle above instead?
% End bins are not full width
diffx = diff(templatex);
half_bin_width = diffx(1) / 2;
first_bin_width = half_bin_width + templatex(1)   - opts.start;
last_bin_width  = half_bin_width - templatex(end) + opts.stop;
templatey(1)   = templatey(1)   * 2 * half_bin_width / first_bin_width;
templatey(end) = templatey(end) * 2 * half_bin_width / last_bin_width;



% Subtract baseline
baseline = mean(templatey(templatex <= 0));
templatey = templatey - baseline;

% Normalize template
templatey = templatey ./ sqrt(templatey * templatey');


if nargout == 0
    plot(templatex, templatey);
    clear templatex templatey;
end