
datarun = load_data('/Volumes/Analysis/2012-09-24-5/d03-06-07-norefit/data003/data003');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun,'load_sta',[],'keep_java_sta',true);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = load_cones_ath(datarun,'d03');

datarun = conepreprocess_wnm(datarun, 'cone_data_ind', 'bayes');
conepreprocess_save(datarun, 'cone_data_ind', 'data', 'date', '2012-09-24-5');

