clear
close all
dbstop if error
dataparam.date='2016-02-17-1/';
dataparam.concatname='data004';
dataparam.mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/RGB-16-1-0.48-11111-119.5.xml';
dataparam.stixel_size = 16;
dataparam.seed = 11111;
dataparam.interval = 1;
dataparam.refresh_rate = 119.5;
fitparam.num_frames = 30;
dataparam.x_dim = 640;
dataparam.y_dim = 320;
frame_width = 640/dataparam.stixel_size;
frame_height = 320/dataparam.stixel_size;
stixels_per_frame = frame_width*frame_height;
num_colors =3;
num_bins = 100;


dataparam.file_name = [dataparam.date, '/', dataparam.concatname,'/', dataparam.concatname];
% dataparam.save_path = ['/Volumes/Lab/Users/crhoades/JitterMovie/', dataparam.date, '/', dataparam.concatname, '/'];

% list specific cell (1), or run for a whole cell type (0)
select_cells = 1;
if select_cells == 1
 dataparam.cell_specification = [183];   
end
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


% slashes = strfind(datarun.names.rrs_neurons_path, '/');
% dataset = datarun.names.rrs_neurons_path(slashes(3)+1:slashes(5)-1);
% to_replace = strfind(dataset, '/');
% dataparam.dataset(to_replace) = '-';

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',true, 'load_sta', 1, 'load_all',true);
opt.load_sta_params.save_rf = 1;
% opt.load_sta_params.frames = 1:fitparam.num_frames;% have to input as a vector list of frames, not the number of frames total, counting backwards
datarun=load_data(datarun,opt);

[inputs, ~, ~] = get_wn_movie_ath(datarun, dataparam.mdf_file, bw);



cell_ids = get_cell_indices(datarun, dataparam.cell_specification);
spikes = datarun.spikes{cell_ids};
%% account for dropped frames

%% 2016-02-17-1 data004
high_point_trigger = 1496;
triggers = [datarun.triggers(1:high_point_trigger); datarun.triggers((high_point_trigger+1):end) - mean(diff(datarun.triggers(1:high_point_trigger)))];
 %1496*100/refresh (refresh is 2)+1 : 392*100/2 (refresh is 2)+50

% inputs = [inputs(:,1:149601), inputs(:,(149601 + (100/dataparam.interval)-1):end)];
    ind = find(spikes > datarun.triggers(high_point_trigger) & spikes <= (datarun.triggers(high_point_trigger+1)+fitparam.num_frames*mean(diff(datarun.triggers(1:high_point_trigger)))/(100/dataparam.interval))); % rate = 120 if refresh =1 or 60 if refresh = 2
    if ~isempty(ind)
      spikes(ind(end)+1:end) = spikes(ind(end)+1:end) - mean(diff(datarun.triggers(1:high_point_trigger)));
      spikes = [spikes(1:ind(1)-1); spikes(ind(end)+1:end)];
    end
%%

% 
% 
% image_width = dataparam.x_dim/dataparam.stixel_size;
% image_height = dataparam.y_dim/dataparam.stixel_size;
% 
% real_frame = zeros(image_width, image_height, num_colors, size(inputs,2));
% if num_colors == 3
%     real_frame(:,:,1,1) = reshape(inputs(1:3:image_width*image_height*3)',image_width, image_height);
%     real_frame(:,:,2,1) = reshape(inputs(2:3:image_width*image_height*3)',image_width, image_height);
%     real_frame(:,:,3,1) = reshape(inputs(3:3:image_width*image_height*3)',image_width, image_height);
% else
%     real_frame(:,:,1,1) = reshape(inputs(1:image_width*image_height)',image_width, image_height);
% 
% end
% 
% 
% pointer = image_width*image_height*num_colors+1;
% %     pointer = 2+25+2;
% i =2;
% 
% %%%% REMOVE THE +6 WHEN NOT 2016-02-17-6 JITTER
% while pointer+image_height*image_width*num_colors-1<size(inputs,2)*size(inputs,1)
%     temp = inputs(pointer:pointer+image_height*image_width*num_colors-1);
%     if num_colors == 3
%         real_frame(:,:,1,i) = reshape(temp(1:3:end), image_width, image_height);
%         real_frame(:,:,2,i) = reshape(temp(2:3:end), image_width, image_height);
%         real_frame(:,:,3,i) = reshape(temp(3:3:end), image_width, image_height);
%     else
%         real_frame(:,:,1,i) = reshape(temp, image_width, image_height);
%     end
% 
%     pointer = pointer+image_height*image_width*num_colors;
%     i = i+1;
% end
% 


upsampled_num_frames = length(triggers)*100/dataparam.interval; 
upsample_factor = dataparam.interval; % should be interval
frames_needed = kron(1:upsampled_num_frames, ones(1,upsample_factor));
frame_spacing = zeros(1, size(frames_needed,2));


for i= 1:length(triggers)-1
    spacing = linspace(triggers(i), triggers(i+1),101);
    frame_spacing(1, (i-1)*100+1:(i-1)*100+100)= spacing(1:end-1); %% assume triggers every 100 frames
end

binned_spikes = zeros(size(spikes,2), size(frames_needed,2)-1);
for j = 1:size(spikes,2)
    for i = 1:size(frames_needed,2)-1
        binned_spikes(j,i) = sum(spikes >= frame_spacing(1,i) & spikes < frame_spacing(1,i+1));
    end
end



sta = double(datarun.stas.stas{cell_ids});
[sig_stixels] = significant_stixels(sta);
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
    strfs = spalloc(size(staframes,1), nframes, numnonzero*nframes);
    for i = 1:length(frames)
        strfs(sig_stixels,i) = staframes(sig_stixels,i);
    end
    
strf = stavect2visbuf(strfs, 40, 20, 1);
diag_gen_signals = zeros(num_stims, num_frames);
for i = 1:num_stims
    diag_gen_signals(i,:)  = double(inputs(:, i))'*strfs;
end

% gen_signals = zeros(num_stims-num_frames + 1, 1);
diag_sum_kernel = eye(num_frames);
gen_signals = conv2(diag_gen_siganls, diag_sum_kernel, 'valid');
% 
% 
% [~, m] = max(sta(:));
% [x,y,z,w] = ind2sub(size(sta), m);
% 
% sta_red = sta(x-2:x+2, y-2:y+2, :,:);
% real_frame_red = real_frame(x-2:x+2, y-2:y+2, :,:);
% 
% [sig_stixels] = significant_stixels(sta);
% 
%  sig_stixels = sig_stixels(:); % Vectorize
% 
%     % expand across color channels
%  sig_stixels = repmat(sig_stixels, size(sta,3), 1);
%     
%     % Get desired STA frames, vectorize them
%     staframes = reshape(sta(:,:,:,1:30), [], 30);
%         numnonzero = sum(sig_stixels);
% frames = 1:30;
% strfs = staframes(sig_stixels,:);
% %  strfs = spalloc(size(staframes,1), 30, numnonzero*30);
%     for i = 1:length(frames)
%         strfs(sig_stixels,i) = staframes(sig_stixels,i);
%     end   
%     
%     movie = inputs(sig_stixels,:);
%     
    
% sta_red = permute(sta_red, [2 1 3 4]);
% 
% 
% sta_red = flip(sta_red, 4);
% sta_red = flip(sta_red, 1);
% sta_red = flip(sta_red, 2);
% sta_red = flip(sta_red, 3);
% 
%     GS_fit = squeeze(convn(real_frame_red, sta_red,'valid'));
% 
% 
%     GS_fit = squeeze(convn(movie, strfs,'valid'));
    binned_spikes = binned_spikes(end-length(GS_fit)+1:end);
    [Y, edges] = quantileranks(GS_fit,num_bins);
    
    for i =1:num_bins
        mean_GS_fit(i) = mean(GS_fit(Y == i));
        spikes_in_range = binned_spikes(Y == i);
        prob(i) = sum(spikes_in_range)./sum(binned_spikes);
    end

    
   f = fit(mean_GS_fit', prob', 'exp1'); % f(x) = a*exp(b*x)
    figure; plot(f,mean_GS_fit, prob)
        xlabel('mean generator signal in bin')
    ylabel('Probability')
    
    coeffval = coeffvalues(f);
    GS_predict = GS_fit;
    % should be GS predict in a real situation
    FR_predict = (coeffval(1).*exp(GS_predict*coeffval(2)))/max((coeffval(1).*exp(GS_predict*coeffval(2)))); % compare to binned_spikes  % not a FR, is a probability of a spike
    
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
    