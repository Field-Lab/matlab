clear
close all
dbstop if error
dataparam.date='2016-02-17-1/';
dataparam.concatname='data000';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-119.5.xml';
num_bins = 10;

summary = figure;
% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
    dataparam.cell_specification = [1456];
end


slashes = strfind(dataparam.mdf_file, '/');
dashes = strfind(dataparam.mdf_file ,'-');
dashes_end = dashes(dashes >slashes(end));

if strcmp(dataparam.mdf_file(slashes(end)+1:slashes(end)+3), 'RGB')
    num_colors =3;
else
    num_colors = 1;
end
dataparam.stixel_size = str2num(dataparam.mdf_file(dashes_end(1)+1:dashes_end(2)-1));
dataparam.interval = str2num(dataparam.mdf_file(dashes_end(2)+1:dashes_end(3)-1));
if length(dashes_end)>=5
    dataparam.seed = str2num(dataparam.mdf_file(dashes_end(4)+1:dashes_end(5)-1));
else
    dataparam.seed = dataparam.mdf_file(dashes_end(4)+1:end);
    dataparam.seed = str2num(dataparam.seed(1:strfind(dataparam.seed, '.xml')-1))
    
end

if isempty(strfind(dataparam.mdf_file, '119.5')) && isempty(strfind(dataparam.mdf_file, '60.35'))
    dataparam.refresh_rate = 120;
else
    dataparam.refresh_rate = dataparam.mdf_file(dashes_end(end)+1:end);
    dataparam.refresh_rate = str2num(dataparam.refresh_rate(1:strfind(dataparam.refresh_rate, '.xml')-1))
end

if length(dashes_end) > 5
    
    dataparam.x_dim = 800;
    dataparam.y_dim = 600;
else
    dataparam.x_dim = 640;
    dataparam.y_dim = 320;
end

dataparam.num_of_interval = 50; % number of * in progress bar
frame_width = dataparam.x_dim/dataparam.stixel_size;
frame_height = dataparam.y_dim/dataparam.stixel_size;
stixels_per_frame = frame_width*frame_height;





dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
dataparam.cell_type = {'all'};

%% END OF INPUT
% load movie

if num_colors == 3
    bw = 0;
else
    bw= 1;
end


datarun.names.rrs_neurons_path=['/Volumes/Analysis/', dataparam.file_name, '.neurons'];
datarun.names.rrs_params_path=['/Volumes/Analysis/', dataparam.file_name, '.params'];
datarun.names.rrs_sta_path = ['/Volumes/Analysis/', dataparam.file_name, '.sta'];


opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
datarun=load_data(datarun,opt);

datarun = get_sta_summaries(datarun, dataparam.cell_specification, ...
    'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
    'marks_params',struct( ...
    'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
    'thresh',4,'robust_std_method',1));


[inputs, ~, ~] = get_wn_movie_ath(datarun, dataparam.mdf_file, bw);

inputs = double(inputs);
inputs(inputs ==-1) = 0.02;
inputs(inputs ==1) = 0.98;

cell_ids = get_cell_indices(datarun, dataparam.cell_specification);

for j = 1:length(cell_ids)
    
    
    spikes{j}=  datarun.spikes{cell_ids(j)};
    
end


%% account for dropped frames
triggers = datarun.triggers;
difference = diff(triggers);
mean_val = mean(difference);
high_point_trigger = find(difference>(mean_val*1.2)==1);
high_point_trigger = sort(high_point_trigger,'descend');
for i = 1:length(high_point_trigger)
    triggers = [datarun.triggers(1:high_point_trigger(i)); datarun.triggers((high_point_trigger(i)+1):end) - mean(diff(datarun.triggers(1:high_point_trigger(i))))];
    %1496*100/refresh (refresh is 2)+1 : 392*100/2 (refresh is 2)+50
    
    % inputs = [inputs(:,1:149601), inputs(:,(149601 + (100/dataparam.interval)-1):end)];
    for j = 1:size(spikes,2)
        ind = find(spikes{j} > datarun.triggers(high_point_trigger(i)) & spikes{j} <= (datarun.triggers(high_point_trigger(i)+1)+num_frames*mean(diff(datarun.triggers(1:high_point_trigger(i))))/(100/dataparam.interval))); % rate = 120 if refresh =1 or 60 if refresh = 2
        if ~isempty(ind)
            spikes{j}(ind(end)+1:end) = spikes{j}(ind(end)+1:end) - mean(diff(datarun.triggers(1:high_point_trigger(i))));
            spikes{j} = [spikes{j}(1:ind(1)-1); spikes{j}(ind(end)+1:end)];
        end
    end
    
end




upsampled_num_frames = length(triggers)*100/dataparam.interval;
upsample_factor = dataparam.interval; % should be interval
frames_needed = kron(1:upsampled_num_frames, ones(1,upsample_factor));
frame_spacing = zeros(1, size(frames_needed,2));


for i= 1:length(triggers)-1
    spacing = linspace(triggers(i), triggers(i+1),101);
    frame_spacing(1, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1); %% assume triggers every 100 frames
end

% binned_spikes = zeros(size(spikes,2), size(frames_needed,2)-1);
% for j = 1:size(spikes,2)
%     for i = 1:size(frames_needed,2)-1
%         binned_spikes(j,i) = sum(spikes{j} >= frame_spacing(1,i) & spikes{j} < frame_spacing(1,i+1));
%     end
% end


for i =1:length(cell_ids)
    sta = double(datarun.stas.stas{cell_ids(i)});
    sig_stixels = datarun.stas.marks{cell_ids(i)};
    sig_stixels = sig_stixels(:); % Vectorize
    
    % identify STA frames to use
    frames = 1:30;
    nframes = length(frames);
    num_frames = 30;
    num_stims = size(inputs,2);
    % Expand across color channels
    sig_stixels = repmat(sig_stixels, size(sta,3), 1);
    
    % Get desired STA frames, vectorize them
    staframes = reshape(sta(:,:,:,frames), [], nframes);
    
    % Sparsify according to sigstix
    numnonzero = sum(sig_stixels);
    strfs{i} = spalloc(size(staframes,1), nframes, numnonzero*nframes);
    for j = 1:length(frames)
        strfs{i}(sig_stixels,j) = staframes(sig_stixels,j);
    end
    
    strfs{i} = stavect2visbuf(strfs{i}, frame_width, frame_height, num_colors ==3);
    diag_gen_signals = zeros(num_stims, num_frames);
    for n = 1:num_stims
        diag_gen_signals(n,:)  = double(inputs(:, n))'*strfs{i};
    end
    
    % gen_signals = zeros(num_stims-num_frames + 1, 1);
    diag_sum_kernel = eye(num_frames);
    gen_signals{i} = conv2(diag_gen_signals, diag_sum_kernel, 'valid');
    
    gen_signals{i} =gen_signals{i} - 0.5* sum(strfs{i}(:));
    refresh_time =  sum(diff(triggers))/((length(triggers)-1)*100)*dataparam.interval;% assumes 20 samples/sec and 100 frames/TTL
   
   spike_rate = histc(datarun.spikes{cell_ids(i)},triggers(1):refresh_time:size(inputs,2)*refresh_time+triggers(1));
    one_cell_spikes= spike_rate(1:end-1);
    start_stim  = 1;
    end_stim = size(inputs,2);
    first_stim = start_stim + nframes - 1; % first_stim is the earliest stimulus frame for which a generator signal value can be computed
    stim_frames = first_stim:end_stim;
    one_cell_spikes = one_cell_spikes(stim_frames,:);
    
    data = [gen_signals{i}, one_cell_spikes];
% Sort the data from lowest generator signal to highest. Column 2 is the
% corresponding spike count for each generator signal
data_sorted = sortrows(data);

num_of_gen_signals = size(data_sorted,1);
gen_signals_per_bin = num_of_gen_signals/num_bins;


bin_edges = [1:ceil(gen_signals_per_bin):size(data_sorted,1), size(data_sorted,1)];
binned_gen_signals=nan(num_bins,1);
mean_spikes=nan(num_bins,1);

for b = 1:length(bin_edges)-1
    binned_gen_signals(b,1) =data_sorted(bin_edges(b)); % mean(data_sorted(bin_edges(i:i+1),1));
   
    %     binned_spikes(i) = sum(data_sorted(bin_edges(i):bin_edges(i+1),2))
    prob_spikes(b,1) = sum(data_sorted(bin_edges(b):bin_edges(b+1),2))/sum(data_sorted(:,2));
    bin_center(b,1) = mean(data_sorted(bin_edges(b):bin_edges(b+1),1));
end
% figure;
% plot(binned_gen_signals, mean_spikes, 'o-')
% xlabel('Generator Signal')
% ylabel('Spike Rate (spikes/bin)')
% title({[date, ' ', concatname]; cell_type{1}; ['Cell ' num2str(cell_specification)]})

%     
% %     binned_spikes = binned_spikes(end-length(gen_signals{i})+1:end);
%     [Y, edges] = quantileranks(gen_signals{i},num_bins);
% %     
%     for k =1:num_bins
%         mean_GS(k) = mean(gen_signals{i}(Y == k));
%         spikes_in_range = one_cell_spikes(Y == k);
%         prob(k) = sum(spikes_in_range)./sum(one_cell_spikes);
%     end
% %     
%     
    f = fit(bin_center, prob_spikes, 'exp1'); % f(x) = a*exp(b*x)
    figure; plot(f,bin_center, prob_spikes)
    title({[dataparam.date, ' ', dataparam.concatname]; ['Cell ', num2str(dataparam.cell_specification(i))]})
    xlabel('mean generator signal in bin')
    ylabel('Probability')
    coeffval = coeffvalues(f);
    gen_spread = linspace(-0.6, 0.6, 0.01);
    
    exp_summary = coeffval(1).*exp(gen_spread*coeffval(2)); % compare to binned_spikes  % not a FR, is a probability of a spike
    figure(summary)
    hold on
    plot(gen_spread, exp_summary)
end

GS_predict = GS_fit;
% should be GS predict in a real situation

num_trials = 19;
for i = 1:num_trials
    gen_spikes(:,i) = rand(length(FR_predict),1) < FR_predict;
end

binned_spikes = binned_spikes(1:1000);
x = 1:length(binned_spikes);

gen_spikes = gen_spikes(1:1000,:);
figure; plot([x(binned_spikes >=1);x(binned_spikes >=1)], [1.21*ones(size(x(binned_spikes >=1)));1.215*ones(size(x(binned_spikes >=1)))], 'k')
hold on
for i = 1:num_trials
    plot([x(gen_spikes(:,i) >=1);x(gen_spikes(:,i) >=1)], [(1+0.01*i)*ones(size(x(gen_spikes(:,i) == 1)));(1+0.01*i+0.008)*ones(size(x(gen_spikes(:,i) == 1)))], 'r')
end
%
%     figure; plot(binned_spikes)
%     hold on
%     plot(FR, 'r')
