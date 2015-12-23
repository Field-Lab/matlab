% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{1,2} = '2009-04-13-5_data005-s3600-s7200_data005_data005-bayes-msf_15.00';

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_10.00';

% kiwi
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-05-13-3/data006/data006';
path_and_name{1,2} = 'kiwi';

% blueberry
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-26-2/data001-s6369-s9551/data001-s6369-s9551';
path_and_name{1,2} = 'blueberry';

% grapes
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2007-03-27-2/data014/data014/data014';
path_and_name{1,2} = '2007-03-27-2_data014_data014_data014-bayes-msf_25.00';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);




cell_counter = 18;
cell_id = datarun.cell_types{1}.cell_ids(cell_counter);
cell_index = get_cell_indices(datarun, cell_id);
datarun = get_sta_summaries(datarun, cell_id, 'keep_stas', false);

%%%%%%%%%%%%%
wdw = 50;
threshold = 0.15;
fig_num = 1;

plot_rf(datarun, cell_id, 'scale', 1, 'polarity', true, 'foa', fig_num)
hold on

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', threshold, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0, 'remove_cones', 'U');   
centers = extras.new_datarun.cones.centers(selection,:);
%plot(centers(:,1), centers(:,2), 'k.', 'MarkerSize', 8)
plot_cone_mosaic(datarun, 'fig_or_axes', fig_num, 'bg_color', [], 'clear', false, 'cone_size', 6)
temp_com = datarun.stas.rf_coms{cell_index};
axis([temp_com(1) - wdw, temp_com(1) + wdw, temp_com(2) - wdw, temp_com(2) + wdw])
hold off
%%%%%%%%%%%%%

%%%%%%%%%%%%%
threshold = 0.1;
fig_num = 2;

plot_rf(datarun, cell_id, 'scale', 1, 'polarity', true, 'foa', fig_num)
hold on

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', threshold, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0, 'remove_cones', 'U');   
centers = extras.new_datarun.cones.centers(selection,:);
plot(centers(:,1), centers(:,2), 'k.', 'MarkerSize', 8)
axis([temp_com(1) - wdw, temp_com(1) + wdw, temp_com(2) - wdw, temp_com(2) + wdw])
hold off
%%%%%%%%%%%%%

%%%%%%%%%%%%%
threshold = 0.05;
fig_num = 3;

plot_rf(datarun, cell_id, 'scale', 1, 'polarity', true, 'foa', fig_num)
hold on

[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', threshold, 'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true, 'scale', 3.0, 'remove_cones', 'U');   
centers = extras.new_datarun.cones.centers(selection,:);
plot(centers(:,1), centers(:,2), 'k.', 'MarkerSize', 8)
axis([temp_com(1) - wdw, temp_com(1) + wdw, temp_com(2) - wdw, temp_com(2) + wdw])
hold off
%%%%%%%%%%%%%




