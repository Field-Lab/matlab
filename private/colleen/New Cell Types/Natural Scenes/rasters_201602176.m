close all
clear all
%% ------- INPUTS -------
% % predict_on_nsem(date, concatname, cell_specification, num_repeats, movie_length,rate)
% date='2016-02-17-6';
% concatname='data010';
% 
% cell_specification = [34];
% num_repeats = 30;
% movie_length = 10000; % in ms
% rate = 119.5;
% 
% 
% %% ------ END OF INPUTS --------
% % file path to save pictures
% filepath=['/Users/colleen/Desktop/NSEM_Rasters/',date,'/',concatname,'/'];
% if ~exist(filepath,'dir')
%     mkdir(filepath);
% end
% 
% % load NSEM runs
% 
% datarun = load_data(fullfile(server_path(),date,concatname, concatname));
% datarun = load_params(datarun,'verbose',1);
% datarun = load_neurons(datarun);
% nmrun = datarun;
% 
% 
% % triggers NSEM
% nm_trigs(:,1)=[1; find(diff(nmrun.triggers)>0.9)+1];
% 
% 
% for i=1:length(cell_specification)
%     
%     % Vision ID of the cell
%     visionID = cell_specification(i);
%     ind = get_cell_indices(nmrun,visionID);
%     
%     
%     % NSEM processing
%     nmrasters = [];
%     cnt = 0; % 1st trial number for each NDF
%     
%     spikes = nmrun.spikes{ind};
%     trigs = nmrun.triggers;
%     beg_points = nm_trigs(:,1);
%     for j=1:num_repeats
%         if j == num_repeats
%             tmp=spikes(spikes>trigs(beg_points(j)))...
%                 - trigs(beg_points(j));
%             nmrasters=[nmrasters tmp'*1000 + movie_length*cnt];
%             cnt = cnt+1;
%         else
%             
%             tmp=spikes(spikes>trigs(beg_points(j)) & spikes<trigs(beg_points(j+1)))...
%                 - trigs(beg_points(j));
%             nmrasters=[nmrasters tmp'*1000 + movie_length*cnt];
%             cnt = cnt+1;
%         end
%     end
% end
% 
% 
% % plot stuff
% fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
% set(fig,'color','white','position',[82 242 1785 856]);
% 
% % plot NSEM
% h=subplot('position',[0.05 0.65, 0.9,0.3]);
% rasterplot(nmrasters/1000,(num_repeats) ,movie_length/1000,h) % convert from ms to seconds
% 
% total_computed = 0;
% n_frames = rate*10;
% 
% i_chunk = 1;
% NS_movie = zeros(320,160,n_frames);
% while total_computed < n_frames
%     load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
%     NS_movie(:,:,(total_computed+1):(total_computed+size(movie_chunk,3))) = movie_chunk;
%     interval(i_chunk) = size(movie_chunk,3);
%     total_computed = size(movie_chunk,3)+total_computed;
%     i_chunk = i_chunk + 1;
% end
% 
% interval_frame = cumsum(interval)/rate;
% 
% % interval_frame = interval_frame(interval_frame < movie_length/1000);
% % interval_frame = [0, interval_frame];
% % for k=1:length(interval_frame)
% %     line([0,0]+interval_frame(k),[0,(num_repeats)*2],'color','b','linewidth',0.3)
% % end
% 
% axis([0 movie_length/1000 0 (num_repeats)*2])
% ylabel('Spikes for each repeat')
% xlabel('time (sec)')
% title([date ', cell ',int2str(visionID)])
% 
% % [psth, bins] = get_psth(spike_times, trial_begin_times, varargin)
% 
% 
% %
% %     % save figure
% % %     print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));
% % export_fig(sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]), '-pdf')
% %
% %         close(fig)
% 
% 
% 
% %
% %make block_frames 1200-100
% % prepped_data = interleaved_data_prep(datarun, 1100, num_repeats, 'cell_spec', cell_specification, 'visual_check', 0);
% % [spike_rate] = IDP_plot_PSTH(prepped_data, 1, [1 1 0], 0);
% 
% 
% 
% % plot stuff
% % figure(fig)
% % h=subplot('position',[0.05 0.35, 0.9,0.3]);
% % plot(spike_rate)


% predict_on_nsem(date, concatname, cell_specification, num_repeats, movie_length,rate)
date='2016-02-17-6';
concatname='data000';

cell_specification = [36];
num_repeats = 30;
movie_length = 10000; % in ms
rate = 119.5;


%% ------ END OF INPUTS --------
% file path to save pictures
filepath=['/Users/colleen/Desktop/NSEM_Rasters/',date,'/',concatname,'/'];
if ~exist(filepath,'dir')
    mkdir(filepath);
end

% load NSEM runs

datarun = load_data(fullfile(server_path(),date,concatname, concatname));
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);
nmrun = datarun;


% triggers NSEM
nm_trigs(:,1)=[1; find(diff(nmrun.triggers)>0.9)+1];
num_repeats = 1;
movie_length = 30;
for i=1:length(cell_specification)
    
    % Vision ID of the cell
    visionID = cell_specification(i);
    ind = get_cell_indices(nmrun,visionID);
    
    
    % NSEM processing
    nmrasters = [];
    cnt = 0; % 1st trial number for each NDF
    
    spikes = nmrun.spikes{ind};
    nmrasters=[nmrasters spikes(spikes>datarun.triggers(1) && spikes<30+datarun.triggers(1))'*1000];
 
    
end


% plot stuff
fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
set(fig,'color','white','position',[82 242 1785 856]);

% plot NSEM
h=subplot('position',[0.05 0.65, 0.9,0.3]);
rasterplot(nmrasters/1000,(num_repeats) ,movie_length/1000,h) % convert from ms to seconds

total_computed = 0;
n_frames = rate*10;

i_chunk = 1;
NS_movie = zeros(320,160,n_frames);
while total_computed < n_frames
    load(['/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSinterval/matfiles/movie_chunk_' num2str(i_chunk) '.mat']);
    NS_movie(:,:,(total_computed+1):(total_computed+size(movie_chunk,3))) = movie_chunk;
    interval(i_chunk) = size(movie_chunk,3);
    total_computed = size(movie_chunk,3)+total_computed;
    i_chunk = i_chunk + 1;
end

interval_frame = cumsum(interval)/rate;

% interval_frame = interval_frame(interval_frame < movie_length/1000);
% interval_frame = [0, interval_frame];
% for k=1:length(interval_frame)
%     line([0,0]+interval_frame(k),[0,(num_repeats)*2],'color','b','linewidth',0.3)
% end

axis([0 movie_length/1000 0 (num_repeats)*2])
ylabel('Spikes for each repeat')
xlabel('time (sec)')
title([date ', cell ',int2str(visionID)])

% [psth, bins] = get_psth(spike_times, trial_begin_times, varargin)


%
%     % save figure
% %     print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]));
% export_fig(sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(visionID)]), '-pdf')
%
%         close(fig)



%
%make block_frames 1200-100
% prepped_data = interleaved_data_prep(datarun, 1100, num_repeats, 'cell_spec', cell_specification, 'visual_check', 0);
% [spike_rate] = IDP_plot_PSTH(prepped_data, 1, [1 1 0], 0);



% plot stuff
% figure(fig)
% h=subplot('position',[0.05 0.35, 0.9,0.3]);
% plot(spike_rate)



%% 
%%%%%%%%%%%%%%%%%%%%%%%% FIT LN MODEL %%%%%%%%%%%%%%%%%%%%%%%% 

clear datarun

dataparam.date='2016-02-17-6/';
dataparam.concatname{1}='data000';
dataparam.mdf_file{1}='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-119.5.xml';
% dataparam.concatname{2}='data028';
% dataparam.mdf_file{2}='/Volumes/Analysis/stimuli/white-noise-xml/RGB-8-2-0.48-11111-119.5.xml';
select_cells = 1;
if select_cells == 1
    dataparam.cell_specification{1} = [36];
end
num_bins = 10;

summary = figure;

for num_datarun = 1:length(dataparam.concatname)
    clear datarun
    color = [0 0 1; 0 1 0];
    % list specific cell (1), or run for a whole cell type (0)
    
    
    
    slashes = strfind(dataparam.mdf_file{num_datarun}, '/');
    dashes = strfind(dataparam.mdf_file{num_datarun} ,'-');
    dashes_end = dashes(dashes >slashes(end));
    
    if strcmp(dataparam.mdf_file{num_datarun}(slashes(end)+1:slashes(end)+3), 'RGB')
        num_colors =3;
    else
        num_colors = 1;
    end
    dataparam.stixel_size = str2num(dataparam.mdf_file{num_datarun}(dashes_end(1)+1:dashes_end(2)-1));
    dataparam.interval = str2num(dataparam.mdf_file{num_datarun}(dashes_end(2)+1:dashes_end(3)-1));
    if length(dashes_end)>=5
        dataparam.seed = str2num(dataparam.mdf_file{num_datarun}(dashes_end(4)+1:dashes_end(5)-1));
    else
        dataparam.seed = dataparam.mdf_file{num_datarun}(dashes_end(4)+1:end);
        dataparam.seed = str2num(dataparam.seed(1:strfind(dataparam.seed, '.xml')-1))
        
    end
    
    if isempty(strfind(dataparam.mdf_file{num_datarun}, '119.5')) && isempty(strfind(dataparam.mdf_file{num_datarun}, '60.35'))
        dataparam.refresh_rate = 120;
    else
        dataparam.refresh_rate = dataparam.mdf_file{num_datarun}(dashes_end(end)+1:end);
        dataparam.refresh_rate = str2num(dataparam.refresh_rate(1:strfind(dataparam.refresh_rate, '.xml')-1))
    end
    
    if length(dashes_end) > 5
        
        dataparam.x_dim = 800;
        dataparam.y_dim = 600;
    else
        dataparam.x_dim = 640;
        dataparam.y_dim = 320;
    end
    
    frame_width = dataparam.x_dim/dataparam.stixel_size;
    frame_height = dataparam.y_dim/dataparam.stixel_size;
    stixels_per_frame = frame_width*frame_height;
    
    
    
    
    
    dataparam.file_name = [dataparam.date, '/', dataparam.concatname{num_datarun},'/', dataparam.concatname{num_datarun}];
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
    
    datarun = get_sta_summaries(datarun, dataparam.cell_specification{num_datarun}, ...
        'verbose',0,'keep_stas',0,'keep_rfs',1,'fig_or_axes',[],...
        'marks_params',struct( ...
        'strength','vector length', 'filter', fspecial('gauss',15,0.7), ...
        'thresh',4,'robust_std_method',1));
    
%     loaded = load('/Volumes/Lab/Users/bhaishahster/pc2016_02_17_6_analysis_fits/OFF_large_smooth_trial_resp_data.mat');
%     inputs = loaded.condMov{1};
    [inputs, ~, ~] = get_wn_movie_ath(datarun, dataparam.mdf_file{num_datarun}, bw);
%     inputs = inputs(2:3:end, :);
    inputs = double(inputs);
    inputs(inputs ==-1) = 0.02;
    inputs(inputs ==1) = 0.98;
    
    cell_ids = get_cell_indices(datarun, dataparam.cell_specification{num_datarun});
    
    for j = 1:length(cell_ids)
        
        
        spikes_new{j}=  datarun.spikes{cell_ids(j)};
        
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
        for j = 1:size(spikes_new,2)
            ind = find(spikes_new{j} > datarun.triggers(high_point_trigger(i)) & spikes_new{j} <= (datarun.triggers(high_point_trigger(i)+1)+num_frames*mean(diff(datarun.triggers(1:high_point_trigger(i))))/(100/dataparam.interval))); % rate = 120 if refresh =1 or 60 if refresh = 2
            if ~isempty(ind)
                spikes_new{j}(ind(end)+1:end) = spikes_new{j}(ind(end)+1:end) - mean(diff(datarun.triggers(1:high_point_trigger(i))));
                spikes_new{j} = [spikes_new{j}(1:ind(1)-1); spikes_new{j}(ind(end)+1:end)];
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
%         sta = sta(:,:,2,:);
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
        
        strfs{i} = stavect2visbuf(strfs{i}, frame_width, frame_height, 1);
        diag_gen_signals = zeros(num_stims, num_frames);
        for n = 1:num_stims
            diag_gen_signals(n,:)  = double(inputs(:, n))'*strfs{i};
        end
        
        % gen_signals = zeros(num_stims-num_frames + 1, 1);
        diag_sum_kernel = eye(num_frames);
        gen_signals{i} = conv2(diag_gen_signals, diag_sum_kernel, 'valid');
        
                gen_signals{i} =gen_signals{i} -0.5* sum(strfs{i}(:));
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
        f = fit(bin_center, prob_spikes, 'exp1'); % f(x) = exp(b*x+a)
        figure; plot(f,bin_center, prob_spikes)
        title({[dataparam.date, ' ', dataparam.concatname{num_datarun}]; ['Cell ', num2str(dataparam.cell_specification{num_datarun}(i))]})
        xlabel('mean generator signal in bin')
        ylabel('Probability')
        coeffval = coeffvalues(f);
        gen_spread = linspace(-0.6, 0.6, 1000);
        
        exp_summary = coeffval(1)*exp(gen_spread*coeffval(2)); % compare to binned_spikes  % not a FR, is a probability of a spike
        figure(summary)
        hold on
        plot(gen_spread, exp_summary, 'Color', color(num_datarun,:))
    end
end
figure(summary)
xlabel('generator signal')
ylabel('Probability of spike')
% title([dataparam.date dataparam.concatname{1}, ' ', dataparam.concatname{2}])


sta = double(datarun.stas.stas{cell_ids});
%         sta = sta(:,:,2,:);
%         for t = 1:30
%         upsamp_image(:,:,t) = imresize(sta(:,:,2,t), 4, 'nearest');
%         end

sig_stixels = datarun.stas.marks{cell_ids};


%         upsamp_sig_stixels = imresize(temp, 4, 'nearest');

%        upsamp_sig_stixels = upsamp_sig_stixels(:); % Vectorize

% identify STA frames to use
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


% Expand across color channels
%         upsamp_sig_stixels = repmat(upsamp_sig_stixels, size(upsamp_image,3), 1);



strfs{i} = stavect2visbuf(strfs{i}, frame_width, frame_height, 1);




% upsamp_image = imresize(image(:,:,2), 2, 'nearest')
% [inputs, refresh, duration] = get_wn_movie_ath(datarun, dataparam.mdf_file{1}, 0);
inputs = inputs(:,1:30*120);
%     loaded = load('/Volumes/Lab/Users/bhaishahster/pc2016_02_17_6_analysis_fits/OFF_large_smooth_trial_resp_data.mat');
%     inputs = loaded.condMov{1};

% inputs = double(inputs);
% inputs(inputs == -0.48) =0.02;
% inputs(inputs == 0.48) = 0.98;
%     inputs(inputs ==-0.2412) = 0.01;
%     inputs(inputs ==0.2412) = 0.4924;

% inputs = NS_movie;



%         inputs = permute(inputs,[2 1 3]);

% inputs_new = reshape(inputs, size(inputs,1)*size(inputs,2),size(inputs,3));
for n = 1:size(inputs,2)
    diag_gen_signals_NS(n,:)  = double(inputs(:, n))'*strfs{1};
end

% gen_signals = zeros(num_stims-num_frames + 1, 1);
diag_sum_kernel = eye(num_frames);
gen_signals{i} = conv2(diag_gen_signals_NS, diag_sum_kernel, 'valid');

        gen_signals{i} =gen_signals{i} - 0.5* sum(strfs{i}(:));

        
predict_prob = coeffval(1)*exp(gen_signals{i}*coeffval(2)); % compare to binned_spikes  % not a FR, is a probability of a spike




%
% GS_predict = GS_fit;
% should be GS predict in a real situation

num_trials = 19;
for i = 1:num_trials
    gen_spikes(i,:) = rand(length(predict_prob),1) < predict_prob;
end

new_spikes = [];
spike_bins = 1/120*30:1/120:10; % fix later
for i =1 :length(spike_bins)
    new_spikes = [new_spikes, (spike_bins(i)+10*find(gen_spikes(:,i) == 1) - 10)'];
end

new_spikes = sort(new_spikes, 'ascend');
h=subplot('position',[0.05 0.35, 0.9,0.3]);

rasterplot(new_spikes/1000,19 ,10/1000,h) % convert from ms to seconds


%
% binned_spikes = binned_spikes(1:1000);
% x = 1:length(binned_spikes);
%
% gen_spikes = gen_spikes(1:1000,:);
% figure; plot([x(binned_spikes >=1);x(binned_spikes >=1)], [1.21*ones(size(x(binned_spikes >=1)));1.215*ones(size(x(binned_spikes >=1)))], 'k')
% hold on
% for i = 1:num_trials
%     plot([x(gen_spikes(:,i) >=1);x(gen_spikes(:,i) >=1)], [(1+0.01*i)*ones(size(x(gen_spikes(:,i) == 1)));(1+0.01*i+0.008)*ones(size(x(gen_spikes(:,i) == 1)))], 'r')
% end
%
%     figure; plot(binned_spikes)
%     hold on
%     plot(FR, 'r')
