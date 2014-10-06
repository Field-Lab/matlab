function hists = collect_spikehists(rasterrun, rgc, stimparams, varargin)
% Flexibly iterate through various stim parameters and collect hists

opts = inputParser();
opts.addParamValue('triggers', []);
opts.addParamValue('histx', []);
opts.addParamValue('start', 0);
opts.addParamValue('stop', []);
opts.addParamValue('binwidth', 0.05);
opts.addParamValue('selectedvals', []);
opts.parse(varargin{:});
opts = opts.Results;

% Load stimulus if necessary
if ~isfield(rasterrun,          'stimulus'), rasterrun = read_stim_lisp_output(rasterrun);             end
if ~isfield(rasterrun.stimulus, 'urgb'),     rasterrun.stimulus = parse_stim_rgbs(rasterrun.stimulus); end

% Fill in default opts
if isempty(opts.triggers), opts.triggers = rasterrun.triggers; end
if isempty(opts.stop)
    if ~isempty(opts.histx), opts.stop = opts.histx(end); 
    else                     opts.stop = min(diff(opts.triggers)); end
end
if isempty(opts.histx), opts.histx = opts.start:opts.binwidth:opts.stop; end

% Initialize data for recursion
trigselect = true(1,length(opts.triggers));

% Recursively build up nested cellarrays of results
hists = collect_spikehists_(rasterrun, rgc, stimparams, trigselect, opts);

% Unpack nested recursion results
hists = cellnest2matrix(hists);


function hists = collect_spikehists_(rasterrun, rgc, stimparams, trigselect, opts)
% Not a real recursion endpoint; just in case the starting point is useless
if isempty(stimparams)
    hists = {};
    return
end

% Get the stimparam and range of values to use
stimparam = stimparams{1};
if ~isempty(opts.selectedvals) && isfield(opts.selectedvals, stimparam)
    stimvals = opts.selectedvals.(stimparam);
else
    % urgbs treated specially
    if strcmp(stimparam, 'urgbs')
        stimvals = 1:length(rasterrun.stimulus.urgbs);
    else
        stimvals = unique(rasterrun.stimulus.(stimparam)(trigselect));
    end
end

% Go through the stimvals, filter selected triggers on each val, then
% either get the hist if this is the last stimparam, or else recurse deeper
nvals = length(stimvals);
for i = 1:nvals
    stimval = stimvals(i);

    % Update the trigger selection that will be passed down
    newtrigselect = trigselect;
    % urbs treated specially
    if strcmp(stimparam, 'urgbs')
        urgbselect = false(1,length(trigselect));
        urgbselect(rasterrun.stimulus.urgbi{stimval}) = true;
        newtrigselect = trigselect & urgbselect;
    else
        newtrigselect = trigselect & rasterrun.stimulus.(stimparam) == stimval;
    end

    % If last stimparam, end recursion
    if length(stimparams) == 1
        % Last level so get the data!
        res = rasterphli(rasterrun, rgc, opts.triggers(newtrigselect), 'start', opts.start, 'stop', opts.stop, 'hist', false);
        hists{i} = histc(res(:,1), opts.histx);
    else
        % Recurse
        hists{i} = collect_spikehists_(rasterrun, rgc, stimparams(2:end), newtrigselect, opts);
    end
end