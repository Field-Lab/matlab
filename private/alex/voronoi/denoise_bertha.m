local_path = '/Volumes/Analysis/';
datarun = load_data([local_path, '2016-02-17-4/d00-05-norefit/data001/data001']);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

[inputs, refresh, duration] = get_wn_movie_ath(datarun, 'BW-2-6-0.48-11111-400x300-60.35.xml');
p = mean(inputs,3);
q = mean(inputs(:,:,1:18000),3);

save([local_path, '/2016-02-17-4/noise_BW_2_6_11111.mat'], 'p', 'q')
