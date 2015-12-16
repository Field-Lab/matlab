
% load data
datarun = load_data('/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data005/data005');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun);
% get raw inputs from white noie movie
[inputs, refresh, ~] = get_wn_movie_ath(datarun, 'BW-20-10-0.48-11111-16x16.xml');
inputs = logical((inputs+0.48)/0.96);

% prepare sta of a cell (choose frame with max STA - it is only to select significant stixels. Can use any alternative code for this.)
sta = datarun.stas.stas{5};
sta = squeeze(sta);
tmp = sta(:,:,26);
a = robust_std(tmp(:)); % threshold in units of std
tmp_stixels = find(abs(tmp)>a*3); % so that we don't care about polarity. Might include a couple of surround stixels though.

% select other stixels for test if nonlinearity will rise at the end
% tmp_stixels = 1:15;

% prepare spikes and inputs
spikes = datarun.spikes{5};
train_spikes = ceil(spikes*1000/refresh);
train_inputs = inputs(tmp_stixels,:)*0.96-0.48;

% set sta params
sta_params.length = 15;
sta_params.offset = 0;
fraction = 0.9;

[mean_unbiased_sta, gensig_bins, mean_nonlinearity, stix_sign]=unbiased_STA_coarse(train_inputs, train_spikes, fraction, sta_params);

figure
plot(mean_unbiased_sta')

figure
plot(mean_nonlinearity)
