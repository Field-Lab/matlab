
% load data from java platform into matlab
datarun = load_data('/Volumes/salk-transfer/Array/Analysis/2008-08-27-5/data008/data008/data008');
datarun = load_neurons(datarun);
datarun = load_sta(datarun, 'load_sta', 'all');
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all');



datarun = get_sta_summaries(datarun, 'all');





