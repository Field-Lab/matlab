datarun = load_data('/Volumes/Analysis/2015-08-17-5/d01-29-norefit/data014/data014');
datarun = load_params(datarun);

matfiles = dir('/Volumes/Lab/Users/Nora/amacrine_rasters/C*.mat');

for file = 1% :length(matfiles)
    load(['/Volumes/Lab/Users/Nora/amacrine_rasters/' matfiles(file).name])
    [xx, yy] = rasterplot(nmrasters, 40, 30000);
end

res_spikes_plot(testmovie, firing_rate)
