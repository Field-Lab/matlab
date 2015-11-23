clear
run_number = {'19', '29', '30', '31', '32', '35', '37', '38', '20', '23', '25', '26', '28', '39'};
neurons = [530,1103,2378,2644,4043,5299,6788,918,3931,4132,5075,5506,6830,7203,532,1006,1010,1597,1880,2421,3398,3559,4658,5342,5493,7087,7145,81,1175,3394,4657,4850,4864,5090,5361,6323,7131];
spikes = zeros(size(run_number, 2), length(neurons));
for r = 1:size(run_number, 2)
    clearvars datarun
    datarun.names.rrs_neurons_path = ['/Volumes/Analysis/2015-09-23-7/d19-39/data0', run_number{r}, '-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039/data0', run_number{r},'-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039.neurons'];
    datarun.names.rrs_params_path = ['/Volumes/Analysis/2015-09-23-7/d19-39/data0', run_number{r},'-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039/data0', run_number{r},'-from-data019_data020_data021_data022_data023_data024_data025_data026_data027_data028_data029_data030_data031_data032_data033_data034_data035_data036_data037_data038_data039.params'];
    % datarun.names.rrs_sta_path = '/Volumes/Analysis/2007-03-01-3/d00_03_06_07_08_09-norefit/data000-from-data000_data003_data006_data007_data008_data009/data000-from-data000_data003_data006_data007_data008_data009.sta';
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_sta', 0, 'load_all',0);

    datarun=load_data(datarun, opt);
    indicies = get_cell_indices(datarun, neurons);
    for i = 1:length(indicies)
        spikes(r,i) = size(datarun.spikes{indicies(i)},1);
    end
end

% sta = datarun.stas.stas{indicies};
 