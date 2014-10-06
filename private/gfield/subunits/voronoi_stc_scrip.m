% load dataset
datarun = load_data('2010-03-05-2','rf-14-apple-gf');
datarun = load_data('2010-09-24-1', 'rf-34-orange-gf');

% load movie_xml_path & other info
datarun = load_index(datarun);

% load spikes times and trigger times
datarun = load_neurons(datarun);

% load java object of movie
datarun = load_java_movie(datarun); 

% specify location of cone map
cd /braid/snle/data/2010-03-05-2/ % apple
cd /braid/snle/data/2010-03-05-2/cone' maps data006'/
cone_map = load('full_cone_map.txt');


%-----------------------------------------
% GET INFORMATION ABOUT STIMULUS


start_time = 0;
end_time = 3599;
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
movie_path = '/snle/acquisition/movie-xml2/RGB-gaussian-1-6-0.48-11111-2023x1.xml';
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
cell_type = 'all';
cell_indices = get_cell_indices(datarun, cell_type);
num_cells = length(cell_indices);

% GET SIGINFICANT PIXELS AND TIME COURSES FOR CELLS
for cll = 1:num_cells
    tmp_sta = datarun.stas.stas{cell_indices(cll)};
    tmp_sig_stixels = significant_stixels(tmp_sta, 'thresh', 5);
    tmp_tc = time_course_from_sta(tmp_sta, tmp_sig_stixels);
    datarun.stas.time_courses{cell_indices(cll)} = tmp_tc;
    datarun.stas.marks{cell_indices(cll)} = tmp_sig_stixels;
end   

%-----------------------------------------
% DEFINE A CELL OF INTEREST AMONG THE CELLS ABOVE, I.E. A PARTICULAR OFF PARASOL
%cell_of_interest = 3;
%apple
cell_id = 826; % off parasol orange
cell_id = 1366;  % off parasol
cell_id = 977;  % off midget apple
cell_id = 1156; % off midget 


cell_list = get_cell_indices(datarun, cell_id);
temp_index = cell_list;
%temp_index = cell_list(cell_of_interest);
% SET THE LENGTH OF THE TEMPORAL KERNAL -- CAN BE SHORTER THAT THE FULL TIME COURSE

% PARAMETERS
% get relevant parameters of STC
num_patches = numel(datarun.stas.stas{1}) ./ size(datarun.stas.stas{1}, 4);
template_length = 8;

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
temp_cone_map = zeros(size(cone_map))+0.5;

% --- PLOT STA ---
figure(101)
for mark = 1:num_marks
    map_indices = find(cone_map == mark_indices(mark));
    temp_cone_map(map_indices) = STA(mark) + 0.5;
end
image(repmat(temp_cone_map, [1 1 3]))


% --- COMPUTE STC ---
STA_ = mean(STE)';
z_STE = STE';
z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
z_filtered_patches = filtered_patches';
z_filtered_patches = z_filtered_patches - STA_*(STA_'*z_filtered_patches)./ norm(STA_)^2;  % project out the mean

% % just subtract the STA from the STE and the filtered_pathes
% z_STE = STE';
% z_STE = z_STE - repmat(STA', 1, size(z_STE,2));
% z_filtered_patches = filtered_patches';
% z_filtered_patches = z_filtered_patches - repmat(STA', 1, size(z_filtered_patches,2));


STC = cov(z_STE') - cov(z_filtered_patches');

% --- FACTORIZE THE STC ---
[PCs, eig_vals_matrix] = eig(STC);
eig_vals = diag(eig_vals_matrix);
[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
PCs = PCs(:, sorted_eig_indices);

% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
figure(8); clf;
plot(eig_vals, 'ko')


% --- PLOT EIGENVECTORS OF THE STC ---
num_dims = 6;
for dm = 1:num_dims
    temp_map = zeros(size(cone_map)) + 0.5;
    temp_pc = PCs(:,dm) ./ (2 * abs(ext(PCs(:,dm)))) + 0.5; 
    for mark = 1:num_marks
        map_indices = find(cone_map == mark_indices(mark));
        temp_map(map_indices) = temp_pc(mark);
    end
    figure(dm)
    image(repmat(temp_map, [1,1,3]))
    axis off
end



num_iters = 50;
time_offsets = ceil(size(filtered_patches,1)*rand(num_iters,1));
shifted_eig_vals = zeros(num_iters, num_marks);
for iter = 1:num_iters
    
    s_spike_times = 1+mod(cell_spike_times+time_offsets(iter), size(filtered_patches,1));
    s_STE = filtered_patches(s_spike_times,:);   
    s_STA_ = mean(s_STE);
    s_STE = s_STE - (s_STE*s_STA_') * s_STA_ ./ norm(s_STA_)^2;  % project out the mean
    s_STC = cov(s_STE) - cov_patches;

    % --- FACTORIZE ---
    [s_PCs, s_eig_val_matrix] = eig(s_STC);
    s_eig_vals = diag(s_eig_val_matrix);
    [s_eig_vals, s_sorted_eig_indices] = sort(s_eig_vals, 'descend');
    s_PCs = s_PCs(:, s_sorted_eig_indices);  

    % store eigvals;
    shifted_eig_vals(iter,:) = s_eig_vals ./ abs(sum(s_eig_vals));
end

num_sds = 3;
mean_shifted_egvals = mean(shifted_eig_vals);
sd_shifted_egvals = std(shifted_eig_vals);
high_shifted_egvals = mean_shifted_egvals + num_sds*sd_shifted_egvals;
low_shifted_egvals = mean_shifted_egvals - num_sds*sd_shifted_egvals;

figure(102); clf;
plot(high_shifted_egvals, 'r-')
hold on
plot(low_shifted_egvals, 'r-')
plot(eig_vals./abs(sum(eig_vals)), 'ko')
axis([0 40 -0.1 0.1])
hold off

% --- PLOT SHIFTED EIGENVECTORS OF THE STC ---
num_dims = 6;
for dm = 1:num_dims
    temp_map = zeros(size(cone_map)) + 0.5;
    temp_pc = s_PCs(:,dm) ./ (2 * abs(ext(s_PCs(:,dm)))) + 0.5; 
    for mark = 1:num_marks
        map_indices = find(cone_map == mark_indices(mark));
        temp_map(map_indices) = temp_pc(mark);
    end
    figure(dm)
    image(repmat(temp_map, [1,1,3]))
end





% -------------------------------------------
% --- CLUSTERING ---

% PROJECT DATA THROUGH PCS;
num_dims = 32;
weights = STE * PCs(:,1:num_dims); 

[Xid, VQs] = kmeans(weights, num_dims, 'MaxIter', 200); 
kernels = PCs(:,1:num_dims) * VQs';

% --- PLOT FACTORS OF THE STC ---
for dm = 1:num_dims
    temp_map = zeros(size(cone_map)) + 0.5;
    temp_pc = kernels(:,dm) ./ (2 * abs(ext(kernels(:,dm)))) + 0.5; 
    for mark = 1:num_marks
        map_indices = find(cone_map == mark_indices(mark));
        temp_map(map_indices) = temp_pc(mark);
    end
    figure(dm)
    imagesc(repmat(temp_map, [1,1,3]))
end    



% -------------------------------------------
% --- NNMF ---

[W, H] = nnmf(STE'+1.45,10);
PCs = W;

num_dims = 10
for dm = 1:num_dims
        temp_pc = PCs(:,dm);
    for mark = 1:num_marks
        map_indices = find(cone_map == mark_indices(mark));
        temp_map(map_indices) = temp_pc(mark);
    end
    figure(dm)
    imagesc(temp_map)
    colormap gray
end




