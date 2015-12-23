function save_raster_data(conerun, rasterrun, conergcs, rasterrgcs, savedir, varargin)
opts = inputParser();
opts.addParamValue('triggers', rasterrun.triggers(1:2:end));
opts.addParamValue('marksthresh', 3.5);
opts.addParamValue('stablerun', []);
opts.addParamValue('stablergcs', []);
opts.parse(varargin{:});
opts = opts.Results;

if ~exist(savedir, 'dir'), mkdir(savedir); end

for i = 1:length(rasterrgcs)
    rgc = rasterrgcs(i);
    if iscell(rgc), rgc = rgc{1}; end
    if isempty(rgc), continue; end
    
    conergc = conergcs(i);
    filename = ['ud_' build_1cone_filename(conerun, rasterrun, conergc)];
    savepath = fullfile(savedir, filename);
    fprintf('Processing %s\n', savepath);
    
    conergcind = get_cell_indices(conerun, conergc);
    polarity = conerun.stas.polarities{conergcind};
    [rax histx histy res] = ...
        ww_raster_plot(rasterrun, rasterrun.stimulus, rgc, polarity, 'triggers', opts.triggers, 'plot', false);
    
    simpconerun = setsimplemarks(conerun, conergc, opts.marksthresh);
    rfstimweight.hairyrf  = calc_stim_rf_weights(conerun, conergc, rasterrun.stimulus.mapims{1});
    rfstimweight.simplerf = calc_stim_rf_weights(simpconerun, conergc, rasterrun.stimulus.mapims{1});
    
    if ~isempty(opts.stablerun) && ~isempty(opts.stablergcs{i})
        stablerun = opts.stablerun;
        stablergc = opts.stablergcs{i};
        simpstablerun = setsimplemarks(stablerun, stablergc, opts.marksthresh);
        rfstimweight.stabilitycheck.hairyrf = calc_stim_rf_weights(stablerun, stablergc, rasterrun.stimulus.mapims{1});
        rfstimweight.stabilitycheck.simplerf = calc_stim_rf_weights(simpstablerun, stablergc, rasterrun.stimulus.mapims{1});
    end
    
    save(savepath, 'histx', 'histy', 'res', 'polarity', 'rfstimweight');
end