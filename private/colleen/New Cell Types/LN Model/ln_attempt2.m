close all
clear all


%% Fit the static nonlinearity


date='2016-02-17-1/data022-data028';
concatname='data026';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-8-1-0.24-11111-119.5.xml';
bw = 1;
cell_specification = [4];
frames = 1:30;
frame_width = 80;
frame_height = 40;
interval  =1;
num_bins = 10;
num_repeats = 30;
time = 10;

file_name = [date, '/', concatname,'/', concatname];

datarun.names.rrs_neurons_path=['/Volumes/Analysis/', file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', file_name, '.sta'];


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_all',true);
datarun=load_data(datarun,opt);

datarun = get_sta_summaries(datarun, cell_specification, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
    'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
    'thresh',4,'robust_std_method',1));

[inputs, ~, ~] = get_wn_movie_ath(datarun, mdf_file, bw);
% inputs = inputs(2:3:end,:); % just the green part of the movie
inputs = double(inputs);
inputs(inputs ==-1) = 0.02;
inputs(inputs ==1) = 0.98;
ind = get_cell_indices(datarun, cell_specification);
spikes_fitting=  datarun.spikes{ind};
triggers = datarun.triggers;

% upsampled_num_frames = length(triggers)*100/interval;
% upsample_factor = interval; % should be interval
% frames_needed = kron(1:upsampled_num_frames, ones(1,upsample_factor));
% frame_spacing = zeros(1, size(frames_needed,2));
% 
% 
% for i= 1:length(triggers)-1
%     spacing = linspace(triggers(i), triggers(i+1),101);
%     frame_spacing(1, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1); %% assume triggers every 100 frames
% end


sta = double(datarun.stas.stas{ind});
% sta = sta(:,:,2,:);
sig_stixels = datarun.stas.marks{ind};
sig_stixels = sig_stixels(:); % Vectorize

% identify STA frames to use
frames = 1:30;
nframes = length(frames);
num_stims = size(inputs,2);
% Expand across color channels
sig_stixels = repmat(sig_stixels, size(sta,3), 1);

% Get desired STA frames, vectorize them
staframes = reshape(sta(:,:,:,frames), [], nframes);

% Sparsify according to sigstix
numnonzero = sum(sig_stixels);
strfs = spalloc(size(staframes,1), nframes, numnonzero*nframes);
for j = 1:length(frames)
    strfs(sig_stixels,j) = staframes(sig_stixels,j);
end

strfs = stavect2visbuf(strfs, frame_width, frame_height, 0);
diag_gen_signals = zeros(num_stims, nframes);
for n = 1:num_stims
    diag_gen_signals(n,:)  = inputs(:, n)'*strfs;
end
diag_sum_kernel = eye(nframes);
gen_signals = conv2(diag_gen_signals, diag_sum_kernel, 'valid');

gen_signals =gen_signals -0.5* sum(strfs(:));
refresh_time =  sum(diff(triggers))/((length(triggers)-1)*100)*interval;% assumes 20 samples/sec and 100 frames/TTL

spike_rate = histc(spikes_fitting,triggers(1):refresh_time:size(inputs,2)*refresh_time+triggers(1));
spikes_fitting = spike_rate(1:end-1);

start_stim  = 1;
end_stim = size(inputs,2);
first_stim = start_stim + nframes - 1; % first_stim is the earliest stimulus frame for which a generator signal value can be computed
stim_frames = first_stim:end_stim;
spikes_fitting = spikes_fitting(stim_frames,:);

data = [gen_signals, spikes_fitting];
% Sort the data from lowest generator signal to highest. Column 2 is the
% corresponding spike count for each generator signal
data_sorted = sortrows(data);

num_of_gen_signals = size(data_sorted,1);
gen_signals_per_bin = num_of_gen_signals/num_bins;


bin_edges = [1:ceil(gen_signals_per_bin):size(data_sorted,1), size(data_sorted,1)];
binned_gen_signals=nan(num_bins,1);
mean_spikes=nan(num_bins,1);

for b = 1:length(bin_edges)-1
%     binned_gen_signals(b,1) =data_sorted(bin_edges(b)); % mean(data_sorted(bin_edges(i:i+1),1));
    
    prob_spikes(b,1) = sum(data_sorted(bin_edges(b):bin_edges(b+1),2))/sum(data_sorted(:,2));
    bin_center(b,1) = mean(data_sorted(bin_edges(b:b+1),1));
end
% 
% ft = fittype( 'a/(1+exp(-b*x))', 'independent', 'x', 'dependent', 'y' );
% options = fitoptions( 'Method', 'NonlinearLeastSquares' );
% options.Display = 'Off';
% options.StartPoint = [10 10];
% 
% f = fit(bin_center, prob_spikes, ft, options); % f(x) = a*exp(b*x)

f = fit(bin_center(1:end-1), prob_spikes(1:end-1), 'exp1'); % f(x) = a*exp(b*x)
figure; plot(f,bin_center, prob_spikes);
figure; plot(bin_center, prob_spikes);

coeffval = coeffvalues(f);



% coeffval = load('coeffval.mat');


%% Plot WN repeats



date='2016-02-17-1/data022-data028';
concatname='data027';
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-8-1-0.24-22222-119.5.xml';
bw = 1;
cell_specification = [4];
frames = 1:30;
frame_width = 80;
frame_height = 40;
interval  =1;
num_bins = 10;



datarun_repeats = load_data(fullfile(server_path(),date,concatname, concatname));
datarun_repeats = load_params(datarun_repeats,'verbose',1);
datarun_repeats = load_neurons(datarun_repeats);
% datarun_repeats = load_sta(datarun_repeats);

% datarun = get_sta_summaries(datarun, cell_specification, ...
%     'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
%     'marks_params',struct( ...
%     'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
%     'thresh',4,'robust_std_method',1));


triggers = datarun_repeats.triggers;
ind = get_cell_indices(datarun_repeats,cell_specification);

spikes = datarun_repeats.spikes{ind};
% trigger_last = 1;
% nm_trigs = 1;
% while triggers(trigger_last) < time*(num_repeats-1)
%     a = find(triggers> triggers(trigger_last) + time);
%     if isempty(a)
%         break;
%     end
%     
%     nm_trigs = [nm_trigs; a(1)];
%     trigger_last = a(1);
% end

nm_trigs(:,1)=[1; find(diff(triggers)>0.84)+1];

rasters = [];
cnt = 0;
beg_points = nm_trigs(:,1);
for j=1:num_repeats
    if j == num_repeats
        tmp=spikes(spikes>triggers(beg_points(j)))...
            - triggers(beg_points(j));
        rasters=[rasters tmp' + time*cnt];
        cnt = cnt+1;
    else
        
        tmp=spikes(spikes>triggers(beg_points(j)) & spikes<triggers(beg_points(j+1)))...
            - triggers(beg_points(j));
        rasters=[rasters tmp' + time*cnt];
        cnt = cnt+1;
    end
end



fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
set(fig,'color','white','position',[82 242 1785 856]);

% plot NSEM
h=subplot('position',[0.05 0.65, 0.9,0.3]);
rasterplot(rasters,num_repeats ,time,h) % convert from s to ms





%% Predict spikes in the first 30 seconds of WN run
nframes = length(frames);

[inputs, ~, ~] = get_wn_movie_ath(datarun, mdf_file, bw);

% Grab the right frames
start_stim = floor(triggers(1)*120)
end_stim = start_stim+ time*120-1;

inputs = inputs(:, start_stim:end_stim);
num_stims = size(inputs,2);

inputs = double(inputs);
% put it in the right range
inputs(inputs ==-1) = 0.02;
inputs(inputs ==1) = 0.98;

% sta = double(datarun.stas.stas{ind});
% sta = sta(:,:,2,:); % only use green channel
% sig_stixels = datarun.stas.marks{ind};
% sig_stixels = sig_stixels(:); % Vectorize
% sig_stixels = repmat(sig_stixels, size(sta,3), 1);
% staframes = reshape(sta(:,:,:,frames), [], nframes);
% 
% numnonzero = sum(sig_stixels);
% strfs = spalloc(size(staframes,1), nframes, numnonzero*nframes);
% for j = 1:nframes
%     strfs(sig_stixels,j) = staframes(sig_stixels,j);
% end
% strfs = stavect2visbuf(strfs, frame_width, frame_height, ~bw);

diag_gen_signals = zeros(num_stims, nframes);
for n = 1:num_stims
    diag_gen_signals(n,:)  = inputs(:, n)'*strfs;
end
diag_sum_kernel = eye(nframes);
gen_signals = conv2(diag_gen_signals, diag_sum_kernel, 'valid');
gen_signals =gen_signals -0.5* sum(strfs(:));

gen_signals = [zeros(29,1);gen_signals];
gen_signals = gen_signals;

%  gen_signals_finer=repmat(gen_signals,1,4)';
%   gen_signals_finer= gen_signals_finer(:);


% gen_signals_finer = resizem(gen_signals, 16*length(gen_signals), 'nearest');
% coeffval = coeffval.coeffval;
% predict_prob = coeffval(1)/(1+exp(-coeffval(2)*gen_signals)); % compare to binned_spikes  % not a FR, is a probability of a spike

predict_prob = coeffval(1)*exp(gen_signals*coeffval(2)); % compare to binned_spikes  % not a FR, is a probability of a spike

num_trials = 19;
for i = 1:num_trials
    gen_spikes(i,:) = rand(length(predict_prob),1) < predict_prob;
end

figure; imagesc(gen_spikes)
new_spikes = [];
spike_bins = 0:1/120:time-1/120; % fix later
for i =1 :length(spike_bins)
    new_spikes = [new_spikes, (spike_bins(i)+time*find(gen_spikes(:,i) == 1) - time)'];
end


new_spikes = sort(new_spikes, 'ascend');
figure(fig);
h=subplot('position',[0.05 0.35, 0.9,0.3]);
figure;
rasterplot(new_spikes, 19 ,time,h) % convert from ms to seconds
