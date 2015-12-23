% For a datarun with single cone RFs, compute the input to each cone.  
%
% set time_offset to change the time offset between spikes and stimuli.  Usually -1 is best.
%
% S = # stimuli
% C = # cones
% N = # cell ids
%
% Results are stored in three variables:
%
% cone_inputs - S x C matrix, each entry a cone activation value
% spike_times - N x S matrix, each entry number of spikes during that stimulus
% spike_rate  - length N cell, each entry a list of spike times (in units of stimulus number)
%
% cone_inputs are saved to a text file located at "save_path"
%

% test_blur = true
write_flag = false;

% LOAD DATARUN

%%
if ~exist('datarun','var')

    switch 3
        case 1;  datarun = load_data('2008-08-26-2','rf-1-blueberry');
        case 2;  datarun = load_data('2008-08-27-0','rf-1-peach');
        case 3;  datarun = load_data('2008-05-13-3','rf-6-kiwi');
        case 4;  datarun = load_data('2008-08-27-5', 'rf-3-plantain');
        case 5;  datarun = load_data('2009-04-13-5', 'rf-5-apricot');
        case 6;  datarun = load_data('2007-03-27-2', 'rf-14-grapes');
        case 7;  datarun = load_data('2010-03-05-2', 'rf-13-apple-gf');
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

% Set path to write saved data
save_path = [single_cone_path,datarun.names.nickname,'/cone_stimulus.txt'];


%%

% SET PARAMETERS

% note start and end times, and set time offset
switch datarun.names.nickname
    case 'blueberry'
        start_time = 6400;
        end_time = 9550;
    case 'peach'
        start_time = 0;
        end_time = 2399;
    case 'kiwi'
        start_time = 0;
        end_time = 7440;
    case 'plantain'
        start_time = 0;
        end_time = 7200;
    case 'apricot'
        start_time = 3600;
        end_time = 7200;
    case 'apple-13'
        start_time = 0;
        end_time = 5082;
    otherwise
        error('start time and end time not set')
end


% refresh time
refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;

% compute start and stop times in terms of stimulus number
start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
%end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
end_stim = datarun.stimulus.java_movie.size;
stims = start_stim:end_stim;

% note stimulus size
field_width = datarun.stimulus.java_movie.getWidth;
field_height = datarun.stimulus.java_movie.getHeight;


datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = import_single_cone_data(datarun, datarun.names.nickname);
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
datarun = get_sta_summaries(datarun, 'all', 'keep_stas', false);


%%
% IDENTIFY SPIKES IN STIMULUS BINS


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

    % take spikes just from the relevant subset of stimulus
%    spike_rate_ = spike_rate_(start_stim:end_stim);

    % translate to spike times (with duplicates for multiple spikes per time bin)
    spike_times_ = [];
    for nn = 1:max(spike_rate_)
        spike_times_ = [spike_times_; find( spike_rate_ > (nn-1) )];
    end

    
    % put into storage variables
    spike_rate(cc,:) = spike_rate_';
    spike_times{cc} = sort(spike_times_);
    
end

%%
fprintf('Computing cone input in frames %d to %d...',start_stim,end_stim)
T=text_waitbar;
start_time_ = clock; % note when it started

% initialize storage variables
cone_inputs = zeros(length(stims),size(Wc,2));


% cycle through each stimulus
for ss = 1:length(stims)

    T=text_waitbar(T,ss/length(stims) - 0.01);

    % note which stim
    this_stim = stims(ss);

    % get new frame
    STAFrame = datarun.stimulus.java_movie.getFrame(this_stim-1);
    new_frame = permute(reshape(STAFrame.getBuffer,3,field_width,field_height),[3 2 1]) - .5;

    %new_frame = new_frame - 0.5;
    new_frame = reshape(new_frame,[],1);

    % convert to cone space
    cone_inputs(ss,:) = full(Wc'*double(new_frame));

end

if 0
    % choose directory to write files
    single_cone_folder = '2011-06-30-0_rf-3-bayes-msf_200.00-BW-2-4';
    cd([single_cone_path,single_cone_folder])

    % write cone generator signals to file
    disp('writing cone generator signals...')
    write_cone_inputs = cone_inputs; 
    size(write_cone_inputs)
    cone_fid = fopen('cone_input.bin','w+', 'b');
    fwrite(cone_fid, write_cone_inputs, 'float32');
    fclose(cone_fid);
    %clear write_cone_inputs

    %fid = fopen('cone_input.bin', 'r');
    %test_out = fread(fid, size(write_cone_inputs), 'float32', 0,'b');
    %test_out = fread(fid, 4, 'float32', 0,'b');
    %fclose(fid);
    
    % write spike times to file
    disp('writing spike counts in each frame...')
    write_spike_rate = spike_rate';
    size(write_spike_rate)
    spike_fid = fopen('spike_counts.bin', 'w+', 'b');
    fwrite(spike_fid, write_spike_rate, 'float32');
    fclose(spike_fid);
    %clear write_spike_rate
    
    %fid = fopen('spike_counts.bin', 'r');
    %spike_out = fread(fid, size(spike_rate), 'float32', 0,'b');
    %test_out = fread(fid, 4, 'float32', 0,'b');
    %fclose(fid);
    
    % write RGC time courses
    num_frames = length(datarun.stas.time_courses{1});
    template_length = num_frames;
    disp('writing RGC time courses...')
    rgc_tcs = zeros(length(cell_indices), length(datarun.stas.time_courses{cell_indices(1)}));
    for cc = 1:length(cell_indices)
        temp_tc = datarun.stas.time_courses{cell_indices(cc)};
        if size(temp_tc,2) == 3 % if rgb stimulus
            temp_tc = sum(temp_tc, 2);
        elseif size(temp_tc,2) ==1 % if bw stimulus
            % cut out region of interest from TC
            temp_tc = temp_tc(num_frames - template_length+1:num_frames);
        elseif isempty(temp_tc)
            warning_message = ['cell ',cc,' has an empty TC, so significant marks'];
            warning(warning_message)
            temp_tc = zeros(1,size(rgc_tcs,2));
        else
            error('number of color channels does not equal either 1 or 3')
        end
        rgc_tcs(cc,:) = temp_tc';
    end
    write_rgc_tcs = rgc_tcs';
    size(write_rgc_tcs)
    time_course_fid = fopen('rgc_tcs.bin', 'w+', 'b');
    fwrite(time_course_fid, write_rgc_tcs, 'float32');
    fclose(time_course_fid);
    
    %time_course_fid = fopen('rgc_tcs.bin');
    %check_tcs = fread(time_course_fid, size(rgc_tcs),'float32', 0, 'b');
    %check_tcs = fread(time_course_fid, 4,'float32', 0, 'b');
    %fclose(time_course_fid);

    % write out cell_ids
    disp('writing cell IDs...')
    write_rgc_ids = datarun.cell_ids';
    tmp_save_path = [single_cone_path,single_cone_folder,'/rgc_IDs','.txt']
    dlmwrite(tmp_save_path, write_rgc_ids, '\t')
    %cell_ids_fid = fopen('rgc_IDs.bin', 'w+', 'b');
    %fwrite(cell_ids_fid, write_cell_ids, 'float32');
    %fclose(cell_ids_fid);

    % write out cone_ids
    disp('writing cone IDs...')
    write_cone_ids = 1:length(datarun.cones.types);
    write_cone_ids = write_cone_ids';
    tmp_save_path = [single_cone_path,single_cone_folder,'/cone_IDs','.txt']
    dlmwrite(tmp_save_path, write_cone_ids, '\t')
    %cone_ids_fid = fopen('cone_IDs.bin', 'w+', 'b');
    %fwrite(cone_ids_fid, write_cone_ids, 'float32');
    %fclose(cone_ids_fid);
    
end

%%
template_length = 7;

% apple
cell_id = 826; % off parasol 648
cell_id = 947;  % on parasol
cell_id = 677; % on midget 841
cell_id = 1156; % off midget 8111
% peach
cell_id = 3334; % a nice off parasol from peach
%plantain
cell_id = 1742; % off midget
cell_id = 887; % off parasol

temp_cell_index = get_cell_indices(datarun, cell_id);

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);
                                        
temp_cone_indices = find(selection);

% get time course from cell of interest
temp_tc = datarun.stas.time_courses{temp_cell_index};
num_frames = size(temp_tc,1);

if size(temp_tc,2) == 3
    temp_tc = temp_tc(num_frames - template_length+1:num_frames,:);
    temp_tc = sum(temp_tc, 2);
elseif size(temp_tc,2) ==1
    % cut out region of interest from TC
    temp_tc = temp_tc(num_frames - template_length+1:num_frames);
else
    error('number of color channels does not equal either 1 or 3')
end

temp_tc = temp_tc - mean(temp_tc);
temp_tc = temp_tc ./ norm(temp_tc);

% flip the impulse response to estimate the temporal filter
reverse_indices = template_length:-1:1;
impulse_filter = temp_tc(reverse_indices);
%impulse_filter = temp_tc;

filtered_cone_inputs = filter(impulse_filter, 1, cone_inputs(:,temp_cone_indices));

   

%-----------------------------------------
% COMSTRUCT THE STE AND COMPUTE STA AND STC

% get the covariance matrix of the stimulus (will be subtracted from the STC)
stimulus_cov = cov(filtered_cone_inputs);

% CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
STE = filtered_cone_inputs(spike_times{temp_cell_index},:);

% get image frame info
cell_com = datarun.cones.rf_fits{temp_cell_index}.center;
window = 15;
frame_begin_x = cell_com(2) - window;
frame_end_x = cell_com(2) + window;
frame_begin_y = cell_com(1) - window;
frame_end_y = cell_com(1) + window;

% --- COMPUTE STA ---
temp_polarity = datarun.stas.polarities{temp_cell_index};
STA = mean(STE);
STA = STA ./ max(STA) ./2;
STA = STA * temp_polarity;

matrix_scale_factor = 10;

% --- PLOT STA ---
figure(101)
image_STA = Wc(:,temp_cone_indices) * STA';
imagesc(norm_image(matrix_scaled_up(reshape(image_STA, [320 320 3]),matrix_scale_factor)))
axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
axis off
axis square
print(101, '~/Desktop/sta.pdf', '-dpdf')

% --- COMPUTE STC ---
STA_ = mean(STE)';
z_STE = STE';
z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
%STE = STE - (STE*STA_') * STA_ ./ norm(STA_)^2;  % project out the mean
z_filtered_cone_inputs = filtered_cone_inputs';
z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean

STC = cov(z_STE') - cov(z_filtered_cone_inputs');

% --- FACTORIZE THE STC ---
[PCs, eig_vals_matrix] = eig(STC);
eig_vals = diag(eig_vals_matrix);
[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
PCs = PCs(:, sorted_eig_indices);

% get indices to meaningful dimensions
keep_indices = find(abs(eig_vals) > 1e-10); % note this number might need to be changed depending on data
% store eigvals;
eig_vals = eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
PCs = PCs(:,keep_indices);

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = 4;
for dm = 1:num_dims
    image_PC = Wc(:,temp_cone_indices) * PCs(:,dm);
    figure(dm)
    image(norm_image(matrix_scaled_up(reshape(image_PC, [320,320,3]),10)))
    axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off
    axis square
    save_path = ['~/Desktop/PC',num2str(dm),'.pdf'];
    print(dm, save_path, '-dpdf')

end
% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
figure(num_dims + 1); clf;
plot(eig_vals, 'ko')



%% significance of STCs by iterative method (Schwartz et al. 2006)


num_sig_dims = 0;

num_iters = 1;
time_offsets = ceil(10000*rand(num_iters,1));
shifted_eig_vals = zeros(num_iters, length(temp_cone_indices) - num_sig_dims - 1);
for iter = 1:num_iters
    
    % bin up spikes for entire duration
    s_spike_times = 1+mod(spike_times{temp_cell_index}+time_offsets(iter), size(filtered_cone_inputs,1));
    s_STE = filtered_cone_inputs(s_spike_times,:); 
    
    s_STA_ = mean(s_STE)';
    z_STE = s_STE';
    z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
    z_filtered_cone_inputs = filtered_cone_inputs';
    z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean
    for dm = 1:num_sig_dims
        z_STE = z_STE - PCs(:,dm)*(PCs(:,dm)'*z_STE)./ norm(PCs(:,dm))^2;
        z_filtered_cone_inputs = z_filtered_cone_inputs - PCs(:,dm)*(PCs(:,dm)'*z_filtered_cone_inputs)./ norm(PCs(:,dm))^2;
    end
    s_STC = cov(z_STE') - cov(z_filtered_cone_inputs');
    
    % --- FACTORIZE ---
    [s_PCs, s_eig_val_matrix] = eig(s_STC);
    s_eig_vals = diag(s_eig_val_matrix);
    [s_eig_vals, s_sorted_eig_indices] = sort(s_eig_vals, 'descend');
    s_PCs = s_PCs(:, s_sorted_eig_indices);  

    % get indices to meaningful dimensions
    keep_indices = find(abs(s_eig_vals) > 1e-10); % note this number might need to be changed depending on data
    % store eigvals;
    shifted_eig_vals(iter,:) = s_eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
end

num_sds = 3;
mean_shifted_egvals = mean(shifted_eig_vals);
sd_shifted_egvals = std(shifted_eig_vals);
high_shifted_egvals = mean_shifted_egvals + num_sds*sd_shifted_egvals;
low_shifted_egvals = mean_shifted_egvals - num_sds*sd_shifted_egvals;

figure(103); clf;
max_eig_val = max(eig_vals);
plot(eig_vals./max_eig_val, 'k.', 'MarkerSize', 20)
hold on
temp_counter = (num_sig_dims+1):1:length(eig_vals);
plot(temp_counter,high_shifted_egvals./max_eig_val, 'r-', 'LineWidth', 2)
plot(temp_counter,low_shifted_egvals./max_eig_val, 'r-', 'LineWidth', 2)
axis([0, size(PCs,1)+1, -0.5 1.1])
title(['PC spectrum ', num2str(num_sig_dims)])
hold off
print(103, '~/Desktop/spectrum.pdf', '-dpdf')

% --- PLOT SHIFTED EIGENVECTORS OF THE STC ---
for dm = 1:4
    image_shifted_PCs = Wc(:,temp_cone_indices) * s_PCs(:,dm);
    figure(dm+50)
    image(norm_image(matrix_scaled_up(reshape(image_shifted_PCs, [320,320,3]),matrix_scale_factor)))
    axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off
    axis square
    save_path = ['~/Desktop/PC_s',num2str(dm),'.pdf'];
    print(dm+50, save_path, '-dpdf')

end


% -------------------------------------------
% --- CLUSTERING ---


% PROJECT DATA THROUGH PCS;
num_dims = 15;
sub_space = [STA' PCs(:,1:num_dims)];
weights = filtered_cone_inputs * sub_space; 

figure
plot(weights(1:10000,1), weights(1:10000,2), 'k.')


[Xid, VQs] = kmeans(weights, num_dims+1, 'MaxIter', 1000); 
kernels = sub_space * VQs';

% --- PLOT FACTORS OF THE STC ---
for dm = 1:4
    image_shifted_PCs = Wc(:,temp_cone_indices) * kernels(:,dm);
    figure(dm+50)
    image(norm_image(reshape(image_shifted_PCs, [320,320,3])))
    axis([frame_begin_y frame_end_y frame_begin_x frame_end_x])
end    


% -------------------------------------------
% --- NNMF ---

[W, H] = nnmf(weights'+0.7,9);
kernels = sub_space * W;

% --- PLOT FACTORS OF THE STC ---
for dm = 1:num_dims+1
    image_shifted_PCs = Wc(:,temp_cone_indices) * kernels(:,dm);
    figure(dm+50)
    image(norm_image(reshape(image_shifted_PCs, [320,320,3])))
    axis([frame_begin_y frame_end_y frame_begin_x frame_end_x])

end   

