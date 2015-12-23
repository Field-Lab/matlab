%% DEFINE DATA RUN TO ANALYZE

clear all

data_name = 'apricot';

plot_matrix_scale_factor = 10;
cone_weight_params.thresh = 0.1;
cone_weight_params.radius = [0 inf];
cone_weight_params.polarity = 1;
cone_weight_params.contiguity = true;
cone_weight_params.scale = 3.0;

save_flag = true;

set(0,'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Helvetica')


%% LOAD DATA
switch data_name

    case 'apple';  datarun = load_data('2010-03-05-2', 'rf-13-apple-gf');
                fruit_name = 'apple/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/apple/'];               
    case 'peach';  datarun = load_data('2008-08-27-0','rf-1-peach');
                fruit_name = 'peach/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/peach/'];                
    case 'plantain';  datarun = load_data('2008-08-27-5', 'rf-3-plantain');
                fruit_name = 'plantain/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/plantain/'];    
    case 'blueberry';  datarun = load_data('2008-08-26-2','rf-1-blueberry');
                fruit_name = 'blueberry/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/blueberry/'];   
    case 'kiwi';  datarun = load_data('2008-05-13-3','rf-6-kiwi');
                fruit_name = 'kiwi/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/kiwi/'];  
    case 'apricot';  datarun = load_data('2009-04-13-5', 'rf-5-apricot');
                fruit_name = 'apricot/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/apricot/'];     
    case 'grapes';  datarun = load_data('2007-03-27-2', 'rf-14-grapes');
                fruit_name = 'grapes/';
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 25;
                parasol_window_size = 50;
                save_path = ['~/Desktop/grapes/'];  

    otherwise
        disp('unknow data_name')
end


% load movie_xml_path & other info
datarun = load_index(datarun);

% load spikes times, trigger times, params, stas and cone info
datarun = load_neurons(datarun);
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, datarun.names.nickname);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
datarun = get_sta_summaries(datarun, {1,2,3,4}, 'keep_stas', false);



%% DETERMINE CELLS TO ANALYZE, PLOT AND PRINT RESULTS

cell_types = {1,2,3,4};
cell_indices = get_cell_indices(datarun,cell_types);


for cc = 1:length(cell_indices)

    temp_cell_id = datarun.cell_ids(cell_indices(cc));
    % -- Determine window size for plotting --
    if ismember(temp_cell_id, datarun.cell_types{3}.cell_ids) || ismember(temp_cell_id, datarun.cell_types{4}.cell_ids)
        window_size = midget_window_size;
    else
        window_size = parasol_window_size;
    end
  
    % --- Prepare figure ---
    figure(1); clf;

     
    % get the portrait for the RF
    temp_rf = get_rf(datarun, temp_cell_id, 'polarity', true);
    image(matrix_scaled_up(norm_image(temp_rf),plot_matrix_scale_factor))
    axis image; axis off
    temp_com = datarun.stas.rf_coms{cell_indices(cc),:};
    bg_window_x = (temp_com(1) - window_size) * plot_matrix_scale_factor;
    ed_window_x = (temp_com(1) + window_size) * plot_matrix_scale_factor;
    bg_window_y = (temp_com(2) - window_size) * plot_matrix_scale_factor;
    ed_window_y = (temp_com(2) + window_size) * plot_matrix_scale_factor;
    axis([bg_window_x ed_window_x bg_window_y ed_window_y])
    title_text = ['Cell ',num2str(temp_cell_id),' RF', temp_cell_id];
    title(title_text)

    

    % save figures to disk
    if save_flag
       disp('writing cell data to disk') 
       print(1, [save_path, num2str(temp_cell_id), '.pdf'], '-dpdf')
    end
        
        
end

