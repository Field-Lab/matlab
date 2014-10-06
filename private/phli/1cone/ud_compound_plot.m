function outstructs = ud_compound_plot(conerun, rasterrun, conergc, varargin)

opts = inputParser();
opts.addParamValue('stabilityrun', []);
opts.addParamValue('colors', []);
opts.addParamValue('cmf', @jet);
opts.addParamValue('triggers', rasterrun.triggers);
opts.addParamValue('rasteropts', {});
opts.addParamValue('rfopts', {});
opts.parse(varargin{:});
opts = opts.Results;

conergcnum = conerun.cell_nums(conergc);

if ~isfield(rasterrun, 'rgc')
    if ~isfield(rasterrun, 'map')
        rasterrun.map = map_ei(conerun, rasterrun, 'master_cell_type', conergc);
    end
    rasterrun.rgc = rasterrun.map{conergcnum};
end

[ud_urgbs ud_colors] = ud_rasterurgbs(rasterrun.stimulus, 'colors', opts.colors, 'cmf', opts.cmf);

% Check polarity, transpose rasters if OFF
if offcell(conerun, conergcnum), ud_urgbs  = ud_urgbs'; end


% Pad to leave room for RFs
rastersize = size(ud_urgbs,1);
rfsize = floor(rastersize/2);
ud_urgbs  = [ud_urgbs,  cell(rastersize,rfsize)];
ud_colors = [ud_colors, cell(rastersize,rfsize)];

% Plot rasters
opts.rasteropts = {'hist_color' ud_colors opts.rasteropts{:}};
outstructs = urgb_raster_subplot(rasterrun, rasterrun.rgc, opts.triggers, ud_urgbs, opts.rasteropts{:});


% Plot RF
opts.rfopts = {'fit' false 'colors' opts.colors 'cmf' opts.cmf opts.rfopts{:}};
width = rastersize+rfsize;
sanesubplot(rastersize, width, {1:rfsize rastersize+(1:rfsize)});
plot_rf_stimmap(conerun, conergc, rasterrun.stimulus.mapnycpoly{1}, opts.rfopts{:});

% Plot stability?
if ~isempty(opts.stabilityrun)
    if ~isfield(opts.stabilityrun, 'rgc')
        if ~isfield(opts.stabilityrun, 'map')
            opts.stabilityrun.map = map_ei(conerun, opts.stabilityrun, 'master_cell_type', conergc);
        end
        opts.stabilityrun.rgc = opts.stabilityrun.map{conergcnum};
    end
    
    sanesubplot(rastersize, width, {rastersize-rfsize+(1:rfsize) rastersize+(1:rfsize)});
    plot_rf_stimmap(opts.stabilityrun, opts.stabilityrun.rgc, rasterrun.stimulus.mapnycpoly{1}, opts.rfopts{:});
end

function bool = offcell(datarun, cellnum)
bool = false;
if isfield(datarun.stas, 'polarities') && ~isempty(datarun.stas.polarities{cellnum})
    bool = datarun.stas.polarities{cellnum} == -1;
end