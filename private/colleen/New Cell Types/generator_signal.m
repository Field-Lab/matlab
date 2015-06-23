%generator signal
clear

date='2012-08-09-3';
concatname='data002';
stim_time = 1800; % seconds
% Wrong Movie Information
file_name = [date, '/', concatname, '/', concatname];
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-1-0.48-11111.xml';

cell_type = {'OFF parasol'};
num_frames = 30; % both have to be run with the name number of frames

cell_specification = 586; %OFF parasol

%% END OF INPUT
folder = cell_type{1};
% file path to save pictures
filepath=['/Users/colleen/Desktop/GeneratorSignal/',date,'/',concatname,'/'];
if ~exist([filepath,folder],'dir')
    mkdir([filepath,folder]);
end

% Wrong Movie
datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];


%% Load Data1
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_sta_params', 1, 'load_all',1);
opt.load_sta_params.save_rf = 1;
opt.load_sta_params.frames =1:30% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

[cell_numbers] = get_cell_indices( datarun, cell_specification );
sta = datarun.stas.stas{cell_numbers};
[sig_stixels] = significant_stixels(sta);
[timecourse, params] = time_course_from_sta(sta, sig_stixels);
bounds = autozoom_to_fit(datarun, cell_numbers)
bounds = round(bounds);

% sta = sta(bounds(3):bounds(4), bounds(1):bounds(2), :,:);
% sta = squeeze(sta(6,5,2,:));
triggers = datarun.triggers;
[mov,height,width,duration,refresh] = get_movie_ath(mdf_file, triggers, 1,2);

[movie] = load_movie(mdf_file, triggers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note stimulus size
field_width = movie.getWidth;
field_height = movie.getHeight;

% identify numbers of things
start_stim = 0;
end_stim = 216000;
num_frames = 30;
num_stims = end_stim - start_stim;
num_cells = 1;


sta = squeeze(sta(bounds(3):bounds(4),bounds(1):bounds(2),:,:));
% stavect = sta;
% stavect_shaped = reshape(stavect, height*width*1, []);
% ind = 1:size(stavect_shaped,1);
% ind = reshape(ind, height, width, []);
% ind = permute(ind, [3 2 1]);
% ind = reshape(ind, [], 1);
% buf = stavect_shaped(ind,:);
% strfs{1} = buf;


% Initialize
diag_gen_signals = cell(num_cells,1);
for cc = 1:num_cells
    diag_gen_signals{cc} = zeros(num_stims, num_frames);
end

% % Step through movie frame by frame
% % The gen_signal for each frame of the STRF is calculated separately for
% % each movie frame here (diag_gen_signals).  The diag_gen_signals matrix is
% % then summed along its diagonals below to calculate the total gen_signal
% % for the whole STRF for each movie frame.
% %for ss = start_stim:end_stim
% for ss = 1:num_stims
%     
%     stimbuf = movie.getFrame(ss-1+start_stim).getBuffer();
%     stimbuf = stimbuf(2:3:end);
% %     if ~opts.isRGB, stimbuf = stimbuf(1:3:end); end
% 
%     % Get generator signal values for each cell
%     for cc = 1:num_cells
%         diag_gen_signals{cc}(ss,:) = double(stimbuf')*strfs{cc};
%     end
% end

stim = zeros(field_height, field_width, 3);
for ss = 1:num_stims
    
    stimbuf = movie.getFrame(ss-1+start_stim).getBuffer();
     stim(:,:,1) = reshape(stimbuf(1:3:end),width,height)'; % store the three color channels
      stim(:,:,2) = reshape(stimbuf(2:3:end),width,height)';
   stim(:,:,3) = reshape(stimbuf(3:3:end),width,height)';
         stim_red = stim(bounds(3):bounds(4), bounds(1):bounds(2),:);
%     stimbuf = stimbuf(2:3:end);
%     if ~opts.isRGB, stimbuf = stimbuf(1:3:end); end

    % Get generator signal values for each cell
    for cc = 1:num_cells
        stim_mult_frames = repmat(stim_red,[1,1,1,30]);
        multiplication = double(stim_mult_frames).*double(sta);
        temp = sum(multiplication,3);
        temp = sum(temp,2);
        temp = squeeze(sum(temp,1));
        diag_gen_signals{cc}(ss,:) = temp;

%         diag_gen_signals{cc}(ss,:) = double(stimbuf')*strfs{cc};
    end
    if mod(ss,10000) == 0
        disp(ss)
    end
end


% Sum across strf frames 
% The result of the above is for each cell a series of vectors giving the
% generator signal across the whole movie for each separate frame of the
% STRF.  We want to sum across STRF frames for the total generator signal,
% which means summing along diagonals.  conv2 with an identity matrix is an
% efficient way to get Matlab to sum down all the diagonals.
gen_signals = zeros(num_stims - num_frames + 1, num_cells);
diag_sum_kernel = eye(num_frames);
size(diag_gen_signals{cc})
size(diag_sum_kernel)
for cc = 1:num_cells
    gen_signals(:,cc) = conv2(diag_gen_signals{cc}, diag_sum_kernel, 'valid');
end


% Unbias gen_signals
% The frames from Java movie are biased by +0.5.  This could be corrected
% in the main loop, but more efficient to compensate here.
for cc = 1:num_cells
    gen_signals(:,cc) = gen_signals(:,cc) - 0.5*sum(sta(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movie = zeros(height, width, 3, 120*stim_time);
% for i =1:120*stim_time
%     F= round(mvi.getFrame(i).getBuffer);
%     movie(:,:,1, i) = round(reshape(F(1:3:end),width,height)'-0.5); % store the three color channels
%     movie(:,:,2, i) = round(reshape(F(2:3:end),width,height)'-0.5);
%     movie(:,:,3, i) = round(reshape(F(3:3:end),width,height)'-0.5);
%     div= mod(i,10000);
%     if div== 0 
%         disp(i)
%     end
%     
% end
% movie = movie(bounds(3):bounds(4), bounds(1):bounds(2), :,:);
% counter = 1;
% generator = zeros(1, 120*stim_time/30);
% % sta = norm(sta);
% for j = 1:30:120*stim_time
%     
% %     stimulus_chunk = movie(:,:,:,j:j+29);
%     stimulus_chunk = movie(j:j+29);
% 
%     multiplication = stimulus_chunk(:).*sta(:);
%     generator(counter) = sum(multiplication(:));
%     counter = counter +1;
% end
% 

spikes = datarun.spikes{cell_numbers};
time = linspace(0, 1800, 216000);

[N, edges] = histcounts(spikes,linspace(time(30), time(1)+time(end), 216000-28));

data = [gen_signals, N'];

data_sorted=  sortrows(data);
range = data_sorted(end,1) - data_sorted(1,1);

num_bins = 40;
%% bins even in number of spikes
spikes_per_bin = length(spikes)/num_bins;

spike_count = 0;
bin_counter = 1;
start_iter = 1;
while bin_counter < num_bins
    for i = start_iter:size(gen_signals,1)
        spike_count = spike_count + data_sorted(i,2);
        if spike_count >= spikes_per_bin

            bin_edges(bin_counter) = data_sorted(i,1);
            total_spike(bin_counter) = spike_count;
            bin_counter = bin_counter + 1;
            spike_count= 0;
            break
        end
    end
    start_iter = i+1
end

bin_centers = zeros(1,39);
bin_centers(1) = (data_sorted(1,1)+bin_edges(1))/2;
bin_centers(2:39) = (bin_edges(2:end) + bin_edges(1:end-1))/2;
bin_edges = [data_sorted(1,1), bin_edges];
width_bins = bin_edges(2:end) - bin_edges(1:end-1);
FR = total_spike./width_bins(1:39);

figure;
plot(bin_centers,FR)
%% bins even in time
% bin_width = range/num_bins;
% 
% bins = data_sorted(1,1):bin_width:data_sorted(end,1);
% [P] = histcounts(data_sorted(:,1), bins)
% 
% count = 0;
% FR = zeros(1, num_bins);
% for i = 1:length(P)-1
%     start = P(i);
%     FR(i) = mean(data_sorted((count+1):(count+start), 2));
%     count = count + start;
% 
% end
% 
% center = (bins(2:end) - bins(1:end-1))./2 + bins(1:end-1);