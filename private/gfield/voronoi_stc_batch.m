%% Analysis Parameters

cell_types = {1,2,3,4};
template_length = 8;
midget_window_size = 20;
parasol_window_size = 40;

save_flag = true;
save_path = ['~/Desktop/apple/'];

set(0,'DefaultAxesFontSize', 10, 'DefaultAxesFontName', 'Helvetica')


%%
% load dataset
datarun = load_data('2010-03-05-2','rf-14-apple-gf');

% load movie_xml_path & other info
datarun = load_index(datarun);

% load spikes times and trigger times
datarun = load_neurons(datarun);

% load java object of movie
datarun = load_java_movie(datarun); 


cd /jacob/snle/data/2010-03-05-2/
cone_map = load('cone_map-data-014.txt');


%-----------------------------------------
% GET INFORMATION ABOUT STIMULUS


start_time = 0;
end_time = 6879;
time_offset = 0;

% refresh time
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;

% compute start and stop times in terms of stimulus number
start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
stims = start_stim:end_stim;

% note stimulus size
field_width = datarun.stimulus.java_movie.getWidth;
field_height = datarun.stimulus.java_movie.getHeight;


% CONSTRUCT THE FULL STIMULUS
movie_path = '/snle/acquisition/movie-xml2/BW-1-8-0.48-11111-1517x1.xml';
[mov,height,width,duration,refresh] = get_movie(movie_path, datarun.triggers, length(stims));
mov = squeeze(mov(:,:,1,:));

mov = mov - 0.5;


%-----------------------------------------
% BIN SPIKES ACCORDING TO STIMULUS FRAMES


% get list of cells
cell_indices = get_cell_indices(datarun,'all');

% initialize storage variables
spike_rate = zeros(length(cell_indices),length(stims));
spike_times = cell(length(cell_indices),1);


% go through each cell
for cc = 1:length(cell_indices)

    cell_index = cell_indices(cc);

    % bin up spikes for entire duration
    spike_rate_ = histc(datarun.spikes{cell_index},datarun.triggers(1):refresh_time:datarun.stimulus.java_movie.size*refresh_time);

    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    %spike_rate_ = circshift(spike_rate_(start_stim:end_stim),time_offset);
    spike_rate_ = spike_rate_(start_stim:end_stim);


    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate_)
        spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
    end
    
    % put into storage variables
    spike_rate(cc,:) = spike_rate_;
    spike_times{cc} = sort(spike_times_);
    
end

%-----------------------------------------
% GET TIME COURSES AND SIGNIFICANT PIXELS (MARKS) FOR CELLS


% LOAD INFORMATION
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta','all');
datarun = set_polarities(datarun);

% DEFINE SET OF CELLS TO PROCESS -- USUALLY OFF PARASOLS {2}
cell_indices = get_cell_indices(datarun, cell_types);
num_cells = length(cell_indices);

% GET SIGINFICANT PIXELS AND TIME COURSES FOR CELLS
for cll = 1:num_cells
    tmp_sta = datarun.stas.stas{cell_indices(cll)};
    tmp_sig_stixels = significant_stixels(tmp_sta, 'thresh', 5);
    if isempty(find(tmp_sig_stixels,1))
        tmp_sig_stixels = significant_stixels(tmp_sta, 'thresh', 3.5);
        warning('lowering significance threshold of cel %d', cll)
    end
    if isempty(find(tmp_sig_stixels,1))
        warning('no significant pixels for cell %d', cll)
    end
    tmp_tc = time_course_from_sta(tmp_sta, tmp_sig_stixels);
    datarun.stas.time_courses{cell_indices(cll)} = tmp_tc;
    datarun.stas.marks{cell_indices(cll)} = tmp_sig_stixels;
end   

%-----------------------------------------
%% define cell types of interest to print figures

cell_indices = get_cell_indices(datarun,cell_types);
num_rgcs = length(cell_indices);


for cc = 1:num_rgcs

    temp_index = cell_indices(cc);

    % PARAMETERS
    % get relevant parameters of STC
    num_patches = numel(datarun.stas.stas{1}) ./ size(datarun.stas.stas{1}, 4);

    % --- SPIKE TIMES FOR CELL OF INTEREST---
    % get spike times for cell of interest
    cell_spike_times = spike_times{temp_index};
    % find all frames that are atleast "template_length" after the beginning of the stimulus.
    cleared_frames = find(cell_spike_times >= template_length);
    cell_spike_times = cell_spike_times(cleared_frames);
    num_spikes = length(cell_spike_times);


    % --- MAKE TEMPORAL KERNEL (IMPULSE RESPONSE) FROM TIME COURSE FOR CELL OF INTEREST ---
    % get TC for cell of interest
    temp_tc = datarun.stas.time_courses{temp_index};
    num_frames = length(temp_tc);
    % cut our region of interest from TC
    temp_tc = temp_tc(num_frames - template_length+1:num_frames);
    %impulse_response = temp_tc;
    impulse_response = temp_tc - mean(temp_tc); % ZERO MEAN
    impulse_response = impulse_response ./ norm(impulse_response); % UNIT LENGTH

    % --- SIGNIFICANT PIXELS FOR CELL OF INTEREST ---
    marks = datarun.stas.marks{temp_index};
    mark_indices = find(full(marks));
    num_marks = length(mark_indices);

    if num_marks < 6
        continue
    end
        

    %-----------------------------------------
    % COMSTRUCT THE STE AND COMPUTE STA AND STC

    % filter by the time course of the cell
    reverse_indices = template_length:-1:1;
    reverse_impulse = impulse_response(reverse_indices);
    filtered_patches = filter(reverse_impulse, 1, mov(mark_indices,:)');

    % get the covariance matrix of the stimulus (will be subtracted from the STC)
    cov_patches = cov(filtered_patches);

    % CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
    STE = filtered_patches(cell_spike_times,:);

    % --- COMPUTE STA ---
    temp_polarity = datarun.stas.polarities{temp_index};
    STA = mean(STE);
    STA = STA ./ max(STA) ./2;
    STA = STA * temp_polarity;
    % MAKE MAP FROM STIMULUS PIXELS TO VORONOI PATCHES

    % --- COMPUTE STC ---
    STA_ = mean(STE)';
    % project the STA out of the STE
    z_STE = STE';
    z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
    % project the STA out of the stimulus
    z_filtered_patches = filtered_patches';
    z_filtered_patches = z_filtered_patches - STA_*(STA_'*z_filtered_patches)./ norm(STA_)^2;  % project out the mean
 
    STC = cov(z_STE') - cov(z_filtered_patches');

    % --- FACTORIZE THE STC ---
    [PCs, eig_vals_matrix] = eig(STC);
    eig_vals = diag(eig_vals_matrix);
    [eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
    PCs = PCs(:, sorted_eig_indices);

    %% PLOTTING
    
    temp_cell_id = datarun.cell_ids(temp_index);
    
    if ismember(temp_cell_id, datarun.cell_types{3}.cell_ids) || ismember(temp_cell_id, datarun.cell_types{4}.cell_ids)
        window_size = midget_window_size;
    else
        window_size = parasol_window_size;
    end
        
     
    % --- PLOT STA ---
    figure(1); clf;
    px = ceil(sqrt(9));
    py = ceil(9)/px;
    plot_axes = subplot_axes(1,[0 0 1 .95],0.05,0.15,px,py);

    
    axes(plot_axes{2})
    temp_cone_map = zeros(size(cone_map))+0.5;
    for mark = 1:num_marks
        map_indices = find(cone_map == mark_indices(mark));
        temp_cone_map(map_indices) = STA(mark) + 0.5;
    end
    image(matrix_scaled_up(repmat(temp_cone_map, [1 1 3]), 5))
    axis off
    axis image
    title('STA from marks')

    % find center of STA
    temp_sta_indices = find(temp_cone_map ~= 0.5);
    [i,j] = ind2sub(size(temp_cone_map), temp_sta_indices);
    mark_logic = false(size(temp_cone_map));
    mark_logic(i,j) = true;
    %temp_com = getfield(regionprops(mark_logic, 'WeightedCentroid'),'WeightedCentroid');
    temp_com = getfield(regionprops(mark_logic, 'Centroid'),'Centroid');

    begin_window_x = (temp_com(1) - window_size) *5;
    end_window_x = (temp_com(1) + window_size) * 5;
    begin_window_y = (temp_com(2) - window_size) * 5;
    end_window_y = (temp_com(2) + window_size) * 5;
    axis([begin_window_x end_window_x begin_window_y end_window_y])

    
    % --- PLOT FULL STA ---
    axes(plot_axes{1})
    temp_cone_map = zeros(size(cone_map))+0.5;
    all_filtered_patches = filter(reverse_impulse, 1, mov');
    full_STE = all_filtered_patches(cell_spike_times,:);
    full_STA = mean(full_STE);
    full_STA = full_STA ./ max(full_STA) ./2;
    full_STA = full_STA * temp_polarity;
    for cn = 1:length(full_STA)
        map_indices = find(cone_map == cn);
        temp_cone_map(map_indices) = full_STA(cn) + 0.5;
    end
    image(matrix_scaled_up(repmat(temp_cone_map, [1 1 3]), 5))
    axis off
    axis image
    title(['cell ', num2str(temp_cell_id), ' : STA'])
    axis([begin_window_x end_window_x begin_window_y end_window_y])

        
    % --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
    axes(plot_axes{3})
    plot(eig_vals, 'ko')
    title('PV spectrum')
    
    % --- PLOT EIGENVECTORS OF THE STC ---
    num_dims = 3;
    for dm = 1:num_dims
        temp_map = zeros(size(cone_map)) + 0.5;
        temp_pc = PCs(:,dm) ./ (2 * abs(ext(PCs(:,dm)))) + 0.5; 
        for mark = 1:num_marks
            map_indices = find(cone_map == mark_indices(mark));
            temp_map(map_indices) = temp_pc(mark);
        end
        axes(plot_axes{dm+3})
        image(matrix_scaled_up(repmat(temp_map, [1,1,3]), 5))
        axis off; axis image; hold on
        axis([begin_window_x end_window_x begin_window_y end_window_y])
        title_text = ['D',num2str((dm))];
        title(title_text)
    end
    
    num_dims = 3;
    dm_counter = 0;
    for dm = num_marks-num_dims+1:num_marks
        temp_map = zeros(size(cone_map)) + 0.5;
        temp_pc = PCs(:,dm) ./ (2 * abs(ext(PCs(:,dm)))) + 0.5; 
        for mark = 1:num_marks
            map_indices = find(cone_map == mark_indices(mark));
            temp_map(map_indices) = temp_pc(mark);
        end
        dm_counter = dm_counter+1;
        axes(plot_axes{dm_counter+6})
        image(matrix_scaled_up(repmat(temp_map, [1,1,3]),5))
        axis off; axis image; hold on
        axis([begin_window_x end_window_x begin_window_y end_window_y])
        title_text = ['D',num2str((num_marks-num_dims+dm_counter))];
        title(title_text)
    end
    
    

    % save figures to disk
    if save_flag
       disp('writing cell data to disk') 
       print(1, [save_path, num2str(temp_cell_id), '.pdf'], '-dpdf')
    end
        
    

end
