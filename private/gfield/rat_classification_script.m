path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2010-03-31-0/data001/data001';

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun = get_sta_fits_from_vision(datarun);
datarun = get_autocorrelations(datarun, 'all','bin_size', 0.002);

stixel_params.thresh = 4;
datarun = get_sta_summaries(datarun, 'all', 'marks_params', stixel_params);

num_cell_types = length(datarun.cell_types);
for ct = 1:num_cell_types
    ct
    datarun.cell_types{ct}.name
end

ct_index_list = 8:1:30;

for ct = 1:length(ct_index_list)
   

    % plot the RF outlines
    plot_rf_summaries(datarun, {ct_index_list(ct)}, 'plot_fits', true, 'coordinates', 'sta', 'foa', 1)

    % get info for time course
    temp_indices = get_cell_indices(datarun, {ct_index_list(ct)});
    num_cells = length(temp_indices);
    tc_length = length(datarun.stas.time_courses{temp_indices(1)});
    tc_matrix = zeros(num_cells, tc_length);
    t_step = (1000./ (datarun.stimulus.monitor_refresh./datarun.stimulus.interval));
    t_steps = -(tc_length-1)*t_step:t_step:0;

    % plot the time course
    figure(2); clf;
    hold on
    for rgc = 1:num_cells
        plot(t_steps,datarun.stas.time_courses{temp_indices(rgc)} ./ std(datarun.stas.time_courses{temp_indices(rgc)}), 'k')
        title(datarun.cell_types{ct_index_list(ct)}.name)
    end

    plot_autocorrelograms(datarun, {ct_index_list(ct)}, 'foa', 3, 'clear_fig', true, 'normalize', true)
    
    % print figures
    ct_name = datarun.cell_types{ct_index_list(ct)}.name;
    print(1, ['~/Desktop/rat-classification/',ct_name,'-fig1.pdf',], '-dpdf')
    print(2, ['~/Desktop/rat-classification/',ct_name,'-fig2.pdf',], '-dpdf')
    print(3, ['~/Desktop/rat-classification/',ct_name,'-fig3.pdf',], '-dpdf')

end

%%
