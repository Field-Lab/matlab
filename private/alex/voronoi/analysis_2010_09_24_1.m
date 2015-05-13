%% load stuff
datarun = load_data('/Volumes/Analysis/2010-09-24-1/data006/data006');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2010-09-24-1/data021/data021');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);


datarun2 = load_data('/Volumes/Analysis/2010-09-24-1/data035/data035');
datarun2 = load_params(datarun2,'verbose',1);
datarun2 = load_sta(datarun2);
datarun2 = set_polarities(datarun2);
datarun2 = load_neurons(datarun2);

datarun3 = load_data('/Volumes/Analysis/2010-09-24-1/data002/data002');
datarun3 = load_params(datarun3,'verbose',1);
datarun3 = load_sta(datarun3);
datarun3 = set_polarities(datarun3);
datarun3 = load_neurons(datarun3);

datarun4 = load_data('/Volumes/Analysis/2010-09-24-1/data005/data005');
datarun4 = load_params(datarun4,'verbose',1);
datarun4 = load_sta(datarun4);
datarun4 = set_polarities(datarun4);
datarun4 = load_neurons(datarun4);

%% find stable cells

figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun1, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun2, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'g')
plot_rf_summaries(datarun3, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'k')

figure
plot_rf_summaries(datarun3, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun4, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun4, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun2, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


figure
plot_rf_summaries(datarun1, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun2, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')
