function save_udraster_data(udstruct, savedir, varargin)
conerun = udstruct.conerun;
rasterrun = udstruct.wwrun;
conergcs = conerun.wwrgcs;
rasterrgcs = rasterrun.wwrgcs;

if isfield(udstruct, 'stablerun')
    stablerun = udstruct.stablerun;
    stablergcs = stablerun.wwrgcs;
    save_raster_data(conerun, rasterrun, conergcs, rasterrgcs, savedir, 'stablerun', stablerun, 'stablergcs', stablergcs, varargin{:});
else
    save_raster_data(conerun, rasterrun, conergcs, rasterrgcs, savedir, varargin{:});
end