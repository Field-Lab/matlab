function hists = collect_grating_spikehists(rasterrun, rgcs, triggers, varargin)

opts = inputParser();
opts.addParamValue('spatialperiods', []);
opts.addParamValue('phaserelperiod', []); % Define spatial phases to use as fractions of spatial period?
opts.addParamValue('urgbs', []);
opts.addParamValue('temporalperiods', []);
opts.addParamValue('histx', []);
opts.addParamValue('start', 0);
opts.addParamValue('stop', []);
opts.addParamValue('binwidth', 0.05);
opts.parse(varargin{:});
opts = opts.Results;

% Load stimulus if necessary
if ~isfield(rasterrun,          'stimulus'), rasterrun = read_stim_lisp_output(rasterrun);             end
if ~isfield(rasterrun.stimulus, 'urgb'),     rasterrun.stimulus = parse_stim_rgbs(rasterrun.stimulus); end

hists = cell(length(rgcs),1);
for i = 1:length(rgcs)
    rgc = rgcs(i);
    if iscell(rgc), rgc = rgc{1}; end
    if isempty(rgc), continue; end
    
    hists{i} = gratings_analysis_rgc(rasterrun, rgc, triggers, opts);
end



function hists = gratings_analysis_rgc(rasterrun, rgc, triggers, opts)
if isempty(opts.stop), opts.stop = min(diff(triggers)); end
histx = opts.histx;
if isempty(histx), histx = opts.start:opts.binwidth:opts.stop; end

stimstruct = rasterrun.stimulus;

spatialperiods = opts.spatialperiods;
if isempty(spatialperiods), spatialperiods = unique(stimstruct.spatialperiods); end
if isempty(opts.urgbs), opts.urgbs = 1:length(stimstruct.urgbs); end
if isempty(opts.temporalperiods), opts.temporalperiods = unique(stimstruct.temporalperiods); end

for i = 1:length(spatialperiods)
    spatialperiod = spatialperiods(i);
    periodtrials = stimstruct.spatialperiods == spatialperiod;
    
    if ~isempty(opts.phaserelperiod)
        phases = spatialperiod * opts.phaserelperiod;
    else
        phases = unique(stimstruct.spatialphases(periodtrials));
    end
    for j = 1:length(phases)
        phase = phases(j);
        periodphasetrials = periodtrials & stimstruct.spatialphases == phase;
        
        for k = opts.urgbs
            
            for l = 1:length(opts.temporalperiods)
                temporalperiod = opts.temporalperiods(l);
                currtrials = periodphasetrials & stimstruct.temporalperiods == temporalperiod;
                currtrials = intersect(find(currtrials), stimstruct.urgbi{k});
                
                res = rasterphli(rasterrun, rgc, triggers(currtrials), 'start', opts.start, 'stop', opts.stop, 'hist', false);
                hists{j,i,k,l} = histc(res(:,1), histx);
            end
        end
    end
end