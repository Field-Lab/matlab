
% load data from plantain
[datarun, cone_info] = load_data_and_cones('plantain', 'sta_summaries', true);


% choose RGCs
datarun = get_sta_fits_from_vision(datarun);
plot_rf_summaries(datarun, {4}, 'foa', 1, 'plot_fits', true, 'array', false, 'label', true)


off_midget_list = [3317 3197 2628 2656 3167];

plot_rf_summaries(datarun, off_midget_list, 'foa', 1, 'clear', false, 'plot_fits', true,...
                    'array', false, 'label', true, 'fit_color', 'r')


%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINE DATA RUN TO ANALYZE

clear all

data_name = 'plantain';

plot_matrix_scale_factor = 1;
cone_weight_params.thresh = 0.1;
cone_weight_params.radius = [0 inf];
cone_weight_params.polarity = 1;
cone_weight_params.contiguity = true;
cone_weight_params.scale = 3.0;
significance_iters = 50;

save_flag = true;

set(0,'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Helvetica')


%% LOAD DATA
switch data_name

    case 'apple';  datarun = load_data('2010-03-05-2', 'rf-13-apple-gf');
                fruit_name = 'apple/';
                num_frames = 76284;
                num_cones = 2287;
                num_rgcs = 295;
                tc_length = 6;
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/apple/'];
                
    case 'peach';  datarun = load_data('2008-08-27-0','rf-1-peach');
                fruit_name = 'peach/';
                num_frames = 48015;
                num_cones = 2107;
                num_rgcs = 328;
                tc_length = 7;
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/peach/'];                

    case 'plantain';  datarun = load_data('2008-08-27-5', 'rf-3-plantain');
                fruit_name = 'plantain/';
                num_frames = 144103;
                num_cones = 1870;
                num_rgcs = 376;
                tc_length = 7;
                data_path = [single_cone_path,fruit_name];
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/plantain/'];     
    case 'blueberry';  datarun = load_data('2008-08-26-2','rf-1-blueberry');
                fruit_name = 'blueberry/';
                data_path = [single_cone_path,fruit_name];
                num_frames = 63046;
                num_cones = 1382;
                num_rgcs = 226;
                tc_length = 6;
                midget_window_size = 22;
                parasol_window_size = 50;
                save_path = ['~/Desktop/blueberry/'];
    case 'apricot';  datarun = load_data('2009-04-13-5', 'rf-5-apricot');
                fruit_name = 'apricot/';
                data_path = [single_cone_path,fruit_name];
                num_frames = 108078;
                num_cones = 6197;
                num_rgcs = 837;
                tc_length = 7;
                midget_window_size = 15;
                parasol_window_size = 40;
                save_path = ['~/Desktop/apricot/'];
    case 'kiwi';  datarun = load_data('2008-05-13-3','rf-6-kiwi');
                fruit_name = 'kiwi/';
                data_path = [single_cone_path,fruit_name];
                num_frames = 111679;
                num_cones = 1941;
                num_rgcs = 216;
                tc_length = 6;
                midget_window_size = 20;
                parasol_window_size = 45;
                save_path = ['~/Desktop/kiwi/'];                    
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

% load cone weights matrix
load([single_cone_path datarun.names.nickname '/Wc.mat'])
  
% --- read cone generator signals from disk ---
disp('Reading file cone_input.bin...')
file_name = [data_path, 'cone_input.bin'];
fid = fopen(file_name, 'r');
cone_inputs = fread(fid, [num_frames, num_cones], 'float32', 0,'b');
fclose(fid);

% --- read spike_rate file from disk ---
disp('Reading file spike_count.bin...')
file_name = [data_path, 'spike_counts.bin'];
fid = fopen(file_name, 'r');
spike_rate = fread(fid, [num_frames, num_rgcs], 'float32', 0,'b');
spike_rate = spike_rate';
fclose(fid);

% --- read the time courses for the RGCs ---
disp('Reading time course file...')
file_name = [data_path, 'rgc_tcs.bin'];
time_course_fid = fopen(file_name, 'r');
rgc_tcs = fread(time_course_fid, [tc_length, num_rgcs],'float32', 0, 'b');
rgc_tcs = rgc_tcs';
fclose(time_course_fid);

% -- load RGC cell IDs ---
disp('reading cell IDs...')
file_name = [data_path,'rgc_IDs.txt'];
rgc_ids = dlmread(file_name, '\t');

% --- load cone IDs ---
disp('reading cone IDs...')
file_name = [data_path, 'cone_IDs.txt'];
cone_ids = dlmread(file_name, '\t');

%% ORGANIZE SPIKE TIMES

% go through each cell and make spike_time cell from spike_rate matrix
for cc = 1:length(rgc_ids)

    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate(cc,:))
        spike_times_ = [spike_times_, find( spike_rate(cc,:) > (nn-1) )];
    end
    
    % put into storage variables
    spike_times{cc} = sort(spike_times_);
    
end

%% DETERMINE CELLS TO ANALYZE, PLOT AND PRINT RESULTS

cell_types = off_midget_list;
cell_indices = get_cell_indices(datarun,cell_types);
width = datarun.stimulus.field_width;
height = datarun.stimulus.field_height;


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
    px = ceil(sqrt(9));
    py = ceil(9)/px;
    plot_axes = subplot_axes(1,[0 0 1 .95],0.05,0.15,px,py);

     
    % get the portrait for the RF
    temp_rf = get_rf(datarun, temp_cell_id, 'polarity', true);
    temp_com = round(datarun.stas.rf_coms{cell_indices(cc),:});
    bg_window_x = (temp_com(1) - window_size); %* plot_matrix_scale_factor;
    ed_window_x = (temp_com(1) + window_size); %* plot_matrix_scale_factor;
    bg_window_y = (temp_com(2) - window_size); %* plot_matrix_scale_factor;
    ed_window_y = (temp_com(2) + window_size); %* plot_matrix_scale_factor;
    % ensure that sta image can fit in window, otherwise skip cell.
    if ed_window_y > height || ed_window_x > width || bg_window_x < 0 || bg_window_y < 0
        continue
    end
    temp_rf = temp_rf(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    axes(plot_axes{1})
    image(matrix_scaled_up(norm_image(temp_rf),plot_matrix_scale_factor))
    axis image; axis off
    title_text = ['Cell ',num2str(temp_cell_id),' STA', temp_cell_id];
    title(title_text)


    
    
    % --- GET STA AND STC IN CONE SPACE ---
    
    
    
    
    % determine significant cones for the given rgc
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, temp_cell_id,cone_weight_params);
    temp_cone_indices = find(selection);
    
    if length(temp_cone_indices) < 7;
        continue
    end

    % get time course from cell of interest, zero the mean and make unit
    % lenth
    temp_tc = rgc_tcs(cell_indices(cc),:);
    temp_tc = temp_tc - mean(temp_tc);
    temp_tc = temp_tc ./ norm(temp_tc);

    % flip the impulse response to estimate the temporal filter
    reverse_indices = tc_length:-1:1;
    impulse_filter = temp_tc(reverse_indices);

    % filter the cone projected stimulus by the RGC time course
    filtered_cone_inputs = filter(impulse_filter, 1, cone_inputs(:,temp_cone_indices));

   

    % --- COMSTRUCT THE STE AND COMPUTE STA AND STC ---

    
    % get the covariance matrix of the stimulus (will be subtracted from the STC)
    stimulus_cov = cov(filtered_cone_inputs);

    % CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
    STE = filtered_cone_inputs(spike_times{cell_indices(cc)},:);


    % --- COMPUTE STA ---
    temp_polarity = datarun.stas.polarities{cell_indices(cc)};
    STA = mean(STE);
    STA = STA ./ max(STA) ./2;
    STA = STA * temp_polarity;


    % --- PLOT STA ---
    axes(plot_axes{2})
    cone_STA = Wc(:,temp_cone_indices) * STA';
    cone_rf = reshape(cone_STA, [height, width, 3]);
    cone_rf = cone_rf(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
    imagesc(norm_image(matrix_scaled_up(cone_rf, plot_matrix_scale_factor)))
    axis off; axis image
    title('cone space from STA')

    
    % --- COMPUTE STC ---
    STA_ = mean(STE)';
    z_STE = STE';
    z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
    z_filtered_cone_inputs = filtered_cone_inputs';
    z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean

    STC = cov(z_STE') - cov(z_filtered_cone_inputs');



    cones_stas{cc} = STA;    
    cones_stcs{cc} = cov(STE);
    stim_cov{cc} = cov(filtered_cone_inputs);
    stim_mean{cc} = mean(filtered_cone_inputs);
    selected_cones{cc} = selection;
    
    %     % --- FACTORIZE THE STC ---
%     [PCs, eig_vals_matrix] = eig(STC);
%     eig_vals = diag(eig_vals_matrix);
%     [eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
%     PCs = PCs(:, sorted_eig_indices);
% 
%     % get indices to meaningful dimensions
%     keep_indices = find(abs(eig_vals) > 1e-10); % note this number might need to be changed depending on data
%     % store eigvals;
%     eig_vals = eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
%     PCs = PCs(:,keep_indices);
% 
%     % --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
%     axes(plot_axes{3})
%     plot(eig_vals, 'ko')
%     
%     % --- PLOT EIGENVECTORS OF THE STC with increased variance---
%     num_dims = 6;
%     for dm = 1:num_dims
%         image_PC = Wc(:,temp_cone_indices) * PCs(:,dm);
%         image_PC = reshape(image_PC, [height, width, 3]);
%         image_PC = image_PC(bg_window_y:ed_window_y,bg_window_x:ed_window_x,:);
%         axes(plot_axes{3+dm})
%         image(norm_image(matrix_scaled_up(image_PC, plot_matrix_scale_factor)))
%         axis off; axis square
%         title(['D',num2str(dm)])
%         drawnow
%     end

    
end

%%
temp_cell_indices = get_cell_indices(datarun, off_midget_list);

for rgc = 1:length(temp_cell_indices)
    spike_list{rgc} = datarun.spikes{temp_cell_indices(rgc)};
%    sta_list{rgc} = datarun.stas.rfs{temp_cell_indices(rgc)};
    rf_coms{rgc} = datarun.stas.rf_coms{temp_cell_indices(rgc)};
end
    

datarun_for_liam.duration = datarun.duration;
datarun_for_liam.cones = datarun.cones;
datarun_for_liam.selected_cones = selected_cones;
datarun_for_liam.cell_ids = off_midget_list;
datarun_for_liam.rf_coms = rf_coms;
datarun_for_liam.spikes = spike_list;
datarun_for_liam.cone_stas = cones_stas;
datarun_for_liam.cone_stcs = cones_stcs;
datarun_for_liam.stim_cov = stim_cov;
datarun_for_liam.stim_mean = stim_mean;











                
                