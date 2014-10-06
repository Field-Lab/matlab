function save_crraster_data(conerun, rasterrun, savedir, varargin)
rasterrun = read_stim_lisp_output(rasterrun);
rasterrun.stimulus = parse_stim_rgbs(rasterrun.stimulus);
save_raster_data(conerun, rasterrun, conerun.rgcs, rasterrun.rgcs, savedir, varargin{:});