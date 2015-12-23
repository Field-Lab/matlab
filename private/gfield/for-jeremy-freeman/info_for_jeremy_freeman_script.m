%% define various quantitites
cone_weight_params.thresh = 0.03;
cone_weight_params.radius = [0 4];
cone_weight_params.polarity = 0;
cone_weight_params.contiguity = false;
cone_weight_params.scale = 3.0;

make_figure_flag = true;
matrix_scale_factor = 10;

dataname = 'plantain';

%% load data
[datarun, cone_info] = load_data_and_cones(dataname, 'sta_summaries', false);

key_info.cone_inputs = cone_info.cone_inputs; %generate signal for each cone (1870 cones in all)
key_info.spike_rate = cone_info.spike_rate; %the number of spikes in each time bin (a bin = 1 monitor frame)
key_info.rgc_ids = datarun.cell_ids; %numbers identify each RGC
key_info.cone_weights = datarun.cones.weights; %the weight from each cone to each RGC
key_info.cone_types = datarun.cones.types; %identifies the cones as L,M, or S.
key_info.cone_locations = datarun.cones.centers; % x-y coordinates for each cone.
key_info.rgc_locations = datarun.stas.rf_coms; % x-y coordinates for the center of mass of each RGC RF
key_info.cell_types = datarun.cell_types; % identifies the type of each RGC (i.e. ON-midget cell)


%%


% set path to java code
visionPath = '/snle/lab/Applications/java-volatile/vision/vision.app/Contents/Resources/Java/Vision.jar';
javaaddpath(visionPath)


% load data from java files into matlab
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2011-10-25-9/data000/data000');
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_sta(datarun,'load_sta', 'all','save_sta', false);

% compute the RFs and time courses for parasols, midgets, and sbc
datarun = get_sta_summaries(datarun, {1,2,3,4,5});





