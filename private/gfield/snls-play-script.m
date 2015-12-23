% Script for computing Static Nonlinearities for 1x1 data

if ~exist('datarun','var')

    switch 2
        case 1;  datarun = load_data('2008-08-26-2','rf-1-blueberry');
        case 2;  datarun = load_data('2008-08-27-0','rf-1-peach');
        case 3;  datarun = load_data('2008-05-13-3','rf-6-kiwi');
    end

    % load movie_xml_path & other info
    datarun = load_index(datarun);
    
    % load spikes times and trigger times
    datarun = load_neurons(datarun);
    
    % load java object of movie
    datarun = load_java_movie(datarun); 
    
    % load cone weights matrix
    load([single_cone_path datarun.names.nickname '/Wc.mat'])
end



% note start and end times, and set time offset
switch datarun.names.nickname
    case 'blueberry'
        start_time = 6400;
        end_time = 9550;
        time_offset = -1;
    case 'peach'
        start_time = 0;
        end_time = 2399;
        time_offset = -1;
    case 'kiwi'
        start_time = 0;
        end_time = 7440;
        time_offset = -1;
    otherwise
        error('start time and end time not set')
end




datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
%datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
%datarun = import_single_cone_data(datarun, '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_10.00');

cell_type = 2;
datarun = get_sta_summaries(datarun, {cell_type}, 'keep_stas', false);

plot_rf_portraits(datarun, {cell_type}, 'plot_radius', 30, 'figure', 10)


cell_id = 3334;
cell_index = get_cell_indices(datarun, cell_id);

% define range over which to compute SNLs
start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;


% set up cell array of STA
strfs = datarun.stas.stas(cell_index);

% set up spikes
% compute spike rate at all times
spike_rate = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);
% store spikes in the relevant region
spikes = spike_rate;



[fit_params,gen_signals,params] = compute_snls(strfs,spikes,datarun.stimulus.java_movie,...
                                    'start_stim', start_stim, 'end_stim', end_stim,...
                                    'num_frames', size(strfs{1}, 4));



datarun = get_snls(datarun, cell_id, 'start_stim', start_stim, 'end_stim', end_stim./10);


snl_struct = datarun.stas.snls{cell_index};

plot_snl_(snl_struct.gen_signal, snl_struct.spikes, 'fit', snl_struct.fit_params, 'foa', 1)












