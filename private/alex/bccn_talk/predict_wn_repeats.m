%% %%%%%%%%%%%%%%%%%%%%%%%%%% NDF2 and NDF1 %%%%%%%%%%%%%%%%%%%%%%%%%%

% NDF2
% data015 BW-16-6-0.48-11111 900s
% data016 NSEM repeats
% data017 BW-8-6-0.48-11111 repeats
% data018 BW-8-6-0.48-11111 1200s

% NDF1
% data019 RGB-16-6-0.48-11111 900s
% data020 NSEM repeats
% data021 BW-8-6-0.48-11111 repeats
% data022 BW-8-6-0.48-11111 1200s

path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';

% ndf2 and 1
for_use = ['data017'; 'data021'];

datarun = cell(2,1);
for i=1:2
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_neurons(datarun{i});
end


for_sta = ['data018'; 'data022'];
starun = cell(2,1);
for i=1:2
    starun{i} = load_data(fullfile(path2data, for_sta(i,:), for_sta(i,:)));
    starun{i} = load_params(starun{i},'verbose',1);
    starun{i} = load_neurons(starun{i});
    starun{i} = load_sta(starun{i}, 'load_sta', 'all');
    starun{i} = set_polarities(starun{i});
end


[inputs{1}, refresh(1), duration(1)] = get_wn_movie_ath(starun{1}, 'BW-8-6-0.48-11111-40x40.xml');
[inputs{2}, refresh(2), duration(2)] = get_wn_movie_ath(starun{2}, 'BW-8-6-0.48-11111-40x40.xml');

ndfs = ['2', '1'];
pos = [960 741 522 303];

[inputs_rep{1}, refresh_rep(1), duration_rep(1)] = get_wn_movie_ath(datarun{1}, 'BW-8-6-0.48-11111-40x40.xml');
[inputs_rep{2}, refresh_rep(2), duration_rep(2)] = get_wn_movie_ath(datarun{2}, 'BW-8-6-0.48-11111-40x40.xml');


sta_params.length = 20;
sta_params.offset = 0;
fraction = 0.9;


%% get parameters of LN model
cell_type = 1;
cellIDs = get_cell_indices(starun{1}, {cell_type});

cols = 'rb';


cellID = 151;
close all

for run_number = 1:2
    
    spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/(refresh(run_number))); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spikes (spikes<602) = []; % to omit first 30s - not to overfit!
    spikes = spikes-601;
    pol = starun{run_number}.stas.polarities{cellID};
    
    sta = reshape(squeeze(starun{run_number}.stas.stas{cellID}),40*40,30);
    sta = pol*sta(:,11:end);
%     figure
%     plot(sta')
    [~,max_loc] = find(sta==max(sta(:)));
    a = robust_std(sta(:,max_loc))*3.5;
    my_stixels = find(sta(:,max_loc)>a);
    
    % [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs, spikes, fraction, sta_params);
    [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(inputs{run_number}(my_stixels,602:end), spikes, fraction, sta_params);
    
    for k=1:length(nonlinearity)
        if isnan(nonlinearity(k))
            nonlinearity(k) = nonlinearity(k-1);
        end
    end
    %
%     figure
%     plot(unbiased_sta')
%     
%     figure
%     plot(gensig_bins(1:end-1)+diff(gensig_bins)/2,nonlinearity, '*-')
%     
    %% get mean response for WN repeats
    
    %%% EVERYTHING IN FRAMES
    % beginnings of repeats
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    
    % spike rate continuous across all trials
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    
    % mean spike rate across trials
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_sr=mean_sr/length(trial_begins);
    
    % figure
    % plot(mean_sr(sta_params.length:end))
    
    %% predict response with LN model
    
    % get raw inputs
    
    one_repeat_end = find(diff(datarun{run_number}.triggers)>0.84,1);
    real_dur = datarun{run_number}.triggers(one_repeat_end)*1000/refresh_rep(run_number); % duration of 1 repeat in frames
    inputs_rep_tmp = inputs_rep{run_number}(my_stixels,1:floor(real_dur));
    
    
    % filter raw inputs with STA - get generator signal
    my_filtered_inputs=zeros(size(inputs_rep_tmp,1),size(inputs_rep_tmp,2)-size(unbiased_sta,2)+1);
    for i=1:size(inputs_rep_tmp,1)
        my_filtered_inputs(i,:)=conv(inputs_rep_tmp(i,:),unbiased_sta(i,:),'valid');
    end
    gen_sig=sum(my_filtered_inputs,1);
    gen_sig=gen_sig/max(gen_sig);
    
    % nonlinearity(end) = nonlinearity(end-1);
    % pass generator signal through nonlinearity
    predicted_rate = zeros(size(gen_sig));
    for j=1:size(gensig_bins,2)-1
        a = find(gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1));
        predicted_rate(a)= nonlinearity(j);
    end
    %
    % figure
    % plot(predicted_rate)
    
    
    tmp = mean_sr(sta_params.length:end-18)';
    tmp = tmp/mean(tmp);
    tmp1= predicted_rate/mean(predicted_rate);
    sse = sum((tmp1-tmp).^2);
    sst = sum((tmp-mean(tmp)).^2);
    r2 = 1 - sse/sst;
    
    if run_number ==1
        tmp_pos = [1,2];
        tt = 739;
    else
        tmp_pos = [4,3];
        tt = 293;
    end
    
    figure
    set(gcf, 'Position', [-1916 tt 1912 366])
    plot(tmp, cols(run_number))
    hold on
    plot(tmp1, 'k')
    title(['r2=', num2str(r2),',  NDF ', ndfs(run_number), ', cell ', int2str(cellID)])
    
    
    figure(3)
    set(gcf, 'Position', [-1917  1 1912 243]);
    hold on
    plot(tmp, cols(run_number))
end


%% predict ndf2
cell_type = 4;
cellIDs = get_cell_indices(starun{1}, {cell_type});

run_number = 1;

clear r2_4by4model 
cell_for_use = [];
cnt = 1;
for cellID = cellIDs
    
    spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/(refresh(run_number))); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spikes (spikes<602) = []; % to omit first 30s - not to overfit!
    spikes = spikes-601;
    pol = starun{run_number}.stas.polarities{cellID};
    
    sta = reshape(squeeze(starun{run_number}.stas.stas{cellID}),40*40,30);
    sta = pol*sta(:,11:end);
    [~,max_loc] = find(sta==max(sta(:)));
    a = robust_std(sta(:,max_loc(1)))*3.5;
    my_stixels = find(sta(:,max_loc(1))>a);
    
    if ~isempty(my_stixels)
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(inputs{run_number}(my_stixels,602:end), spikes, fraction, sta_params);
        
        for k=1:length(nonlinearity)
            if isnan(nonlinearity(k))
                nonlinearity(k) = nonlinearity(k-1);
            end
        end
        
        %%% EVERYTHING IN FRAMES
        trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
        spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        spike_rate=zeros(max(diff(trial_begins))*20,1);
        spikes(spikes>max(diff(trial_begins))*20)=[];
        while ~isempty(spikes)
            [c, ia, ic] = unique(spikes);
            spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
            spikes(ia)=[];
        end
        mean_sr=0;
        for i=1:length(trial_begins)
            mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
        end
        mean_ndf4=mean_sr/length(trial_begins);
        mean_ndf4_norm = mean_ndf4/mean(mean_ndf4);
        
        % get raw inputs
        
        one_repeat_end = find(diff(datarun{run_number}.triggers)>0.84,1);
        real_dur = datarun{run_number}.triggers(one_repeat_end)*1000/refresh_rep(run_number); % duration of 1 repeat in frames
        inputs_rep_tmp = inputs_rep{run_number}(my_stixels,1:floor(real_dur));
        
        
        % filter raw inputs with STA - get generator signal
        my_filtered_inputs=zeros(size(inputs_rep_tmp,1),size(inputs_rep_tmp,2)-size(unbiased_sta,2)+1);
        for i=1:size(inputs_rep_tmp,1)
            my_filtered_inputs(i,:)=conv(inputs_rep_tmp(i,:),unbiased_sta(i,:),'valid');
        end
        gen_sig=sum(my_filtered_inputs,1);
        gen_sig=gen_sig/max(gen_sig);
        
        % pass generator signal through nonlinearity
        predicted_rate = zeros(size(gen_sig));
        for j=1:size(gensig_bins,2)-1
            a = find(gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1));
            predicted_rate(a)= nonlinearity(j);
        end
        
        
        mean_ndf4_norm_cut = mean_ndf4_norm(sta_params.length:end-18)';
        
        predicted_norm= predicted_rate/mean(predicted_rate);
        sse = sum((predicted_norm-mean_ndf4_norm_cut).^2);
        sst = sum((mean_ndf4_norm_cut-mean(mean_ndf4_norm_cut)).^2);
        r2_4by4model(cnt) = 1 - sse/sst;
        cell_for_use = [cell_for_use cellID];
        cnt=cnt+1;
    end
    
end

clear r2_4by3 r2_3by4
cnt = 1;
for cellID = cell_for_use
    
    run_number = 1;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf4=mean_sr/length(trial_begins);
    mean_ndf4_norm = mean_ndf4/mean(mean_ndf4);
    
    
    run_number = 2;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf3=mean_sr/length(trial_begins);
    mean_ndf3_norm = mean_ndf3/mean(mean_ndf3);
    
%     figure
%     plot(mean_ndf3_norm)
%     hold on
%     plot(mean_ndf4_norm)
        
    % how much of ndf4 variance is explained by ndf3
    sse = sum((mean_ndf4_norm-mean_ndf3_norm).^2);
    sst = sum((mean_ndf3_norm-mean(mean_ndf3_norm)).^2);
    r2_3by4(cnt) = 1 - sse/sst;
    
    % how much of ndf3 variance is explained by ndf4
    sse = sum((mean_ndf3_norm-mean_ndf4_norm).^2);
    sst = sum((mean_ndf4_norm-mean(mean_ndf4_norm)).^2);
    r2_4by3(cnt) = 1 - sse/sst;
    
    cnt = cnt+1;
end


% 
% figure
% % plot(r2_4by3)
% hold on
% plot(r2_4by4model)

figure
hold on
% ON parasols
plot(r2_4by3, r2_4by4model,'b*')

% on midgets
plot(r2_4by3, r2_4by4model,'+r')

% off parasols
plot(r2_4by3, r2_4by4model,'vk')

% off midgets
plot(r2_4by3, r2_4by4model,'^g')

axis([-1 1 -1 1])
line([-1,1], [0,0])
line([0,0], [-1,1])
line([-1,1], [-1,1])
xlabel('wrong light level')
ylabel('model')


%% predict ndf1
cell_type = 4;
cellIDs = get_cell_indices(starun{1}, {cell_type});

run_number = 2;

clear r2_1by1model 
cell_for_use = [];
cnt = 1;
for cellID = cellIDs
    
    spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/(refresh(run_number))); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spikes (spikes<602) = []; % to omit first 30s - not to overfit!
    spikes = spikes-601;
    pol = starun{run_number}.stas.polarities{cellID};
    
    sta = reshape(squeeze(starun{run_number}.stas.stas{cellID}),40*40,30);
    sta = pol*sta(:,11:end);
    [~,max_loc] = find(sta==max(sta(:)));
    a = robust_std(sta(:,max_loc(1)))*3.5;
    my_stixels = find(sta(:,max_loc(1))>a);
    
    if ~isempty(my_stixels)
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(inputs{run_number}(my_stixels,602:end), spikes, fraction, sta_params);
        
        for k=1:length(nonlinearity)
            if isnan(nonlinearity(k))
                nonlinearity(k) = nonlinearity(k-1);
            end
        end
        
        %%% EVERYTHING IN FRAMES
        trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
        spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        spike_rate=zeros(max(diff(trial_begins))*20,1);
        spikes(spikes>max(diff(trial_begins))*20)=[];
        while ~isempty(spikes)
            [c, ia, ic] = unique(spikes);
            spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
            spikes(ia)=[];
        end
        mean_sr=0;
        for i=1:length(trial_begins)
            mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
        end
        mean_ndf1=mean_sr/length(trial_begins);
        mean_ndf1_norm = mean_ndf1/mean(mean_ndf1);
        
        % get raw inputs
        
        one_repeat_end = find(diff(datarun{run_number}.triggers)>0.84,1);
        real_dur = datarun{run_number}.triggers(one_repeat_end)*1000/refresh_rep(run_number); % duration of 1 repeat in frames
        inputs_rep_tmp = inputs_rep{run_number}(my_stixels,1:floor(real_dur));
        
        
        % filter raw inputs with STA - get generator signal
        my_filtered_inputs=zeros(size(inputs_rep_tmp,1),size(inputs_rep_tmp,2)-size(unbiased_sta,2)+1);
        for i=1:size(inputs_rep_tmp,1)
            my_filtered_inputs(i,:)=conv(inputs_rep_tmp(i,:),unbiased_sta(i,:),'valid');
        end
        gen_sig=sum(my_filtered_inputs,1);
        gen_sig=gen_sig/max(gen_sig);
        
        % pass generator signal through nonlinearity
        predicted_rate = zeros(size(gen_sig));
        for j=1:size(gensig_bins,2)-1
            a = find(gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1));
            predicted_rate(a)= nonlinearity(j);
        end
        
        
        mean_ndf1_norm_cut = mean_ndf1_norm(sta_params.length:end-18)';
        
        predicted_norm= predicted_rate/mean(predicted_rate);
        sse = sum((predicted_norm-mean_ndf1_norm_cut).^2);
        sst = sum((mean_ndf1_norm_cut-mean(mean_ndf1_norm_cut)).^2);
        r2_1by1model(cnt) = 1 - sse/sst;
        cell_for_use = [cell_for_use cellID];
        cnt=cnt+1;
    end
    
end

clear r2_2by1 r2_1by2
cnt = 1;
for cellID = cell_for_use
    
    run_number = 1;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf2=mean_sr/length(trial_begins);
    mean_ndf2_norm = mean_ndf2/mean(mean_ndf2);
    
    
    run_number = 2;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf1=mean_sr/length(trial_begins);
    mean_ndf1_norm = mean_ndf1/mean(mean_ndf1);
    
%     figure
%     plot(mean_ndf3_norm)
%     hold on
%     plot(mean_ndf4_norm)
        
    % how much of ndf1 variance is explained by ndf2
    sse = sum((mean_ndf2_norm-mean_ndf1_norm).^2);
    sst = sum((mean_ndf1_norm-mean(mean_ndf1_norm)).^2);
    r2_1by2(cnt) = 1 - sse/sst;
    
    % how much of ndf2 variance is explained by ndf1
    sse = sum((mean_ndf1_norm-mean_ndf2_norm).^2);
    sst = sum((mean_ndf2_norm-mean(mean_ndf2_norm)).^2);
    r2_2by1(cnt) = 1 - sse/sst;
    
    cnt = cnt+1;
end


% 
% figure
% % plot(r2_4by3)
% hold on
% plot(r2_4by4model)

figure
hold on
% ON parasols
plot(r2_1by2, r2_1by1model,'b*')

% on midgets
plot(r2_1by2, r2_1by1model,'+r')

% off parasols
plot(r2_1by2, r2_1by1model,'vk')

% off midgets
plot(r2_1by2, r2_1by1model,'^g')

axis([-3 1 -3 1])
line([-3,3], [0,0])
line([0,0], [-3,3])
line([-3,3], [-3,3])
xlabel('wrong light level')
ylabel('model')



%% %%%%%%%%%%%%%%%%%%%%%%%%%% NDF4 and NDF3 %%%%%%%%%%%%%%%%%%%%%%%%%%

% NDF4
% data002 BW-20-8-0.48-11111 1400s
% data003 BW-20-8-0.48-11111 repeats
% data004 BW-16-8-0.48-11111 repeats
% data005 NSEM repeats
% data006 BW-16-8-0.48-11111 1200s

% NDF3
% data007 BW-16-8-0.48-11111 1800s
% data008 BW-16-8-0.48-11111 repeats
% data009 BW-10-6-0.48-11111 repeats
% data010 NSEM repeats

path2data = '/Volumes/Analysis/2015-08-17-1/d01-29-norefit/';

% ndf4 and 3
for_use = ['data004'; 'data008'];
ndfs = ['4', '3'];
pos = [960 741 522 303];


datarun = cell(2,1);
for i=1:2
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_neurons(datarun{i});
end


for_sta = ['data006'; 'data007'];
starun = cell(2,1);
for i=1:2
    starun{i} = load_data(fullfile(path2data, for_sta(i,:), for_sta(i,:)));
    starun{i} = load_params(starun{i},'verbose',1);
    starun{i} = load_neurons(starun{i});
    starun{i} = load_sta(starun{i}, 'load_sta', 'all');
    starun{i} = set_polarities(starun{i});
end



[inputs{1}, refresh(1), duration(1)] = get_wn_movie_ath(starun{1}, 'BW-16-8-0.48-11111-20x20.xml');
[inputs{2}, refresh(2), duration(2)] = get_wn_movie_ath(starun{2}, 'BW-16-8-0.48-11111-20x20.xml');

[inputs_rep{1}, refresh_rep(1), duration_rep(1)] = get_wn_movie_ath(datarun{1}, 'BW-16-8-0.48-11111-20x20.xml');

[inputs_rep{2}, refresh_rep(2), duration_rep(2)] = get_wn_movie_ath(datarun{2}, 'BW-16-8-0.48-11111-20x20.xml');


%% get parameters of LN model
cell_type = 2;
cellIDs = get_cell_indices(starun{1}, {cell_type});
 

sta_params.length = 20;
sta_params.offset = 0;
fraction = 0.9;
cols = 'rb';


cellID = 145;
close all

for run_number = 1:2
    
    spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/(refresh(run_number))); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spikes (spikes<602) = []; % to omit first 30s - not to overfit!
    spikes = spikes-601;
    pol = starun{run_number}.stas.polarities{cellID};
    
    sta = reshape(squeeze(starun{run_number}.stas.stas{cellID}),20*20,30);
    sta = pol*sta(:,11:end);
%     figure
%     plot(sta')
    a = robust_std(sta(:,16))*3.5;
    my_stixels = find(sta(:,16)>a);
    
    % [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs, spikes, fraction, sta_params);
    [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(inputs{run_number}(my_stixels,602:end), spikes, fraction, sta_params);
    
    for k=1:length(nonlinearity)
        if isnan(nonlinearity(k))
            nonlinearity(k) = nonlinearity(k-1);
        end
    end
    %
    % figure
    % plot(unbiased_sta')
    %
    % figure
    % plot(gensig_bins(1:end-1)+diff(gensig_bins)/2,nonlinearity, '*-')
    
    %% get mean response for WN repeats
    
    %%% EVERYTHING IN FRAMES
    % beginnings of repeats
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    
    % spike rate continuous across all trials
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    
    % mean spike rate across trials
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_sr=mean_sr/length(trial_begins);
    
    % figure
    % plot(mean_sr(sta_params.length:end))
    
    %% predict response with LN model
    
    % get raw inputs
    
    one_repeat_end = find(diff(datarun{run_number}.triggers)>0.84,1);
    
    real_dur = datarun{run_number}.triggers(one_repeat_end)*1000/refresh_rep(run_number); % duration of 1 repeat in frames
    
    inputs_rep_tmp = inputs_rep{run_number}(my_stixels,1:floor(real_dur));
    
    
    % filter raw inputs with STA - get generator signal
    my_filtered_inputs=zeros(size(inputs_rep_tmp,1),size(inputs_rep_tmp,2)-size(unbiased_sta,2)+1);
    for i=1:size(inputs_rep_tmp,1)
        my_filtered_inputs(i,:)=conv(inputs_rep_tmp(i,:),unbiased_sta(i,:),'valid');
    end
    gen_sig=sum(my_filtered_inputs,1);
    gen_sig=gen_sig/max(gen_sig);
    
    % nonlinearity(end) = nonlinearity(end-1);
    % pass generator signal through nonlinearity
    predicted_rate = zeros(size(gen_sig));
    for j=1:size(gensig_bins,2)-1
        a = find(gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1));
        predicted_rate(a)= nonlinearity(j);
    end
    %
    % figure
    % plot(predicted_rate)
    
    
    tmp = mean_sr(sta_params.length:end-14)';
    tmp = tmp/mean(tmp);
    tmp1= predicted_rate/mean(predicted_rate);
    sse = sum((tmp1-tmp).^2);
    sst = sum((tmp-mean(tmp)).^2);
    r2 = 1 - sse/sst;
    
    if run_number ==1
        tmp_pos = [1,2];
        tt = 739;
    else
        tmp_pos = [4,3];
        tt = 293;
    end
    
    figure
    set(gcf, 'Position', [-1916 tt 1912 366])
    plot(tmp, cols(run_number))
    hold on
    plot(tmp1, 'k')
    title(['r2=', num2str(r2),',  NDF ', ndfs(run_number), ', cell ', int2str(cellID)])
    
    
    figure(3)
    set(gcf, 'Position', [-1917  1 1912 243]);
    hold on
    plot(tmp, cols(run_number))
end



%
%
% fr = predicted_rate; % spikes per frames
% dt = 0.83/2; % frame
% nBins = length(fr); % 10 ms spike train
% nTrials = 20; % number of simulations
% spikeMat = rand(nTrials, nBins) < repmat(fr*dt,20,1);
%
%
% figure
% set(gcf, 'Position', [-1915 pos(tmp_pos(1)) 1878 145])
% hold off
% for j=1:14
%     tmp = find(spikeMat(j,:));
%     plot(tmp,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
%     hold on
% end
% plot(predicted_rate/max(predicted_rate))
% title(['predicted, NDF ', ndfs(run_number), ', cell ', int2str(cellID)])
%
%
% % spike rate continuous across all trials
% spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh);
% spikes(spikes<sta_params.length-sta_params.offset)=[];
%
% splitSpikes=cell(14,1);
% for j=6:19
%     tmp=spikes(spikes>=trial_begins(j)+1 & spikes<trial_begins(j+1))...
%         - trial_begins(j);
%     tmp(tmp<20)=[];
%     splitSpikes{j-5}=tmp-20;
% end
% figure
% set(gcf, 'Position', [-1915 pos(tmp_pos(2)) 1878 145])
% hold off
% for j=1:14
%     plot(splitSpikes{j},zeros(length(splitSpikes{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
%     hold on
% end
% tmp = mean_sr(sta_params.length:end-14)';
% plot(tmp/max(tmp))
% title(['real, NDF ', ndfs(run_number), ', cell ', int2str(cellID)])
%
%
%
%


%% predict ndf4
cell_type = 4;
cellIDs = get_cell_indices(starun{1}, {cell_type});

run_number = 1;

clear r2_4by4model 
cell_for_use = [];
cnt = 1;
for cellID = cellIDs
    
    spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/(refresh(run_number))); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spikes (spikes<602) = []; % to omit first 30s - not to overfit!
    spikes = spikes-601;
    pol = starun{run_number}.stas.polarities{cellID};
    
    sta = reshape(squeeze(starun{run_number}.stas.stas{cellID}),20*20,30);
    sta = pol*sta(:,11:end);
    [~,max_loc] = find(sta==max(sta(:)));
    a = robust_std(sta(:,max_loc(1)))*3.5;
    my_stixels = find(sta(:,max_loc(1))>a);
    
    if ~isempty(my_stixels)
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(inputs{run_number}(my_stixels,602:end), spikes, fraction, sta_params);
        
        for k=1:length(nonlinearity)
            if isnan(nonlinearity(k))
                nonlinearity(k) = nonlinearity(k-1);
            end
        end
        
        %%% EVERYTHING IN FRAMES
        trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
        spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        spike_rate=zeros(max(diff(trial_begins))*20,1);
        spikes(spikes>max(diff(trial_begins))*20)=[];
        while ~isempty(spikes)
            [c, ia, ic] = unique(spikes);
            spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
            spikes(ia)=[];
        end
        mean_sr=0;
        for i=1:length(trial_begins)
            mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
        end
        mean_ndf4=mean_sr/length(trial_begins);
        mean_ndf4_norm = mean_ndf4/mean(mean_ndf4);
        
        % get raw inputs
        
        one_repeat_end = find(diff(datarun{run_number}.triggers)>0.84,1);
        real_dur = datarun{run_number}.triggers(one_repeat_end)*1000/refresh_rep(run_number); % duration of 1 repeat in frames
        inputs_rep_tmp = inputs_rep{run_number}(my_stixels,1:floor(real_dur));
        
        
        % filter raw inputs with STA - get generator signal
        my_filtered_inputs=zeros(size(inputs_rep_tmp,1),size(inputs_rep_tmp,2)-size(unbiased_sta,2)+1);
        for i=1:size(inputs_rep_tmp,1)
            my_filtered_inputs(i,:)=conv(inputs_rep_tmp(i,:),unbiased_sta(i,:),'valid');
        end
        gen_sig=sum(my_filtered_inputs,1);
        gen_sig=gen_sig/max(gen_sig);
        
        % pass generator signal through nonlinearity
        predicted_rate = zeros(size(gen_sig));
        for j=1:size(gensig_bins,2)-1
            a = find(gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1));
            predicted_rate(a)= nonlinearity(j);
        end
        
        
        mean_ndf4_norm_cut = mean_ndf4_norm(sta_params.length:end-14)';
        
        predicted_norm= predicted_rate/mean(predicted_rate);
        sse = sum((predicted_norm-mean_ndf4_norm_cut).^2);
        sst = sum((mean_ndf4_norm_cut-mean(mean_ndf4_norm_cut)).^2);
        r2_4by4model(cnt) = 1 - sse/sst;
        cell_for_use = [cell_for_use cellID];
        cnt=cnt+1;
    end
    
end

clear r2_4by3 r2_3by4
cnt = 1;
for cellID = cell_for_use
    
    run_number = 1;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf4=mean_sr/length(trial_begins);
    mean_ndf4_norm = mean_ndf4/mean(mean_ndf4);
    
    
    run_number = 2;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf3=mean_sr/length(trial_begins);
    mean_ndf3_norm = mean_ndf3/mean(mean_ndf3);
    
%     figure
%     plot(mean_ndf3_norm)
%     hold on
%     plot(mean_ndf4_norm)
        
    % how much of ndf4 variance is explained by ndf3
    sse = sum((mean_ndf4_norm-mean_ndf3_norm).^2);
    sst = sum((mean_ndf3_norm-mean(mean_ndf3_norm)).^2);
    r2_3by4(cnt) = 1 - sse/sst;
    
    % how much of ndf3 variance is explained by ndf4
    sse = sum((mean_ndf3_norm-mean_ndf4_norm).^2);
    sst = sum((mean_ndf4_norm-mean(mean_ndf4_norm)).^2);
    r2_4by3(cnt) = 1 - sse/sst;
    
    cnt = cnt+1;
end

figure
hold on
% ON parasols
plot(r2_4by3, r2_4by4model,'b*')

% on midgets
plot(r2_4by3, r2_4by4model,'+r')

% off parasols
plot(r2_4by3, r2_4by4model,'vk')

% off midgets
plot(r2_4by3, r2_4by4model,'^g')

axis([-1 1 -1 1])
line([-1,1], [0,0])
line([0,0], [-1,1])
line([-1,1], [-1,1])
xlabel('wrong light level')
ylabel('model')



%% predict ndf3
cell_type = 4;
cellIDs = get_cell_indices(starun{1}, {cell_type});

run_number = 2;

clear r2_3by3model 
cell_for_use = [];
cnt = 1;
for cellID = cellIDs
    
    spikes=ceil((starun{run_number}.spikes{cellID}-starun{run_number}.triggers(1))*1000/(refresh(run_number))); % spikes in frames
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spikes (spikes<602) = []; % to omit first 30s - not to overfit!
    spikes = spikes-601;
    pol = starun{run_number}.stas.polarities{cellID};
    
    sta = reshape(squeeze(starun{run_number}.stas.stas{cellID}),20*20,30);
    sta = pol*sta(:,11:end);
    [~,max_loc] = find(sta==max(sta(:)));
    a = robust_std(sta(:,max_loc(1)))*3.5;
    my_stixels = find(sta(:,max_loc(1))>a);
    
    if ~isempty(my_stixels)
        
        [unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA_coarse(inputs{run_number}(my_stixels,602:end), spikes, fraction, sta_params);
        
        for k=1:length(nonlinearity)
            if isnan(nonlinearity(k))
                nonlinearity(k) = nonlinearity(k-1);
            end
        end
        
        %%% EVERYTHING IN FRAMES
        trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
        spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
        spikes(spikes<sta_params.length-sta_params.offset)=[];
        spike_rate=zeros(max(diff(trial_begins))*20,1);
        spikes(spikes>max(diff(trial_begins))*20)=[];
        while ~isempty(spikes)
            [c, ia, ic] = unique(spikes);
            spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
            spikes(ia)=[];
        end
        mean_sr=0;
        for i=1:length(trial_begins)
            mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
        end
        mean_ndf3=mean_sr/length(trial_begins);
        mean_ndf3_norm = mean_ndf3/mean(mean_ndf3);
        
        % get raw inputs
        
        one_repeat_end = find(diff(datarun{run_number}.triggers)>0.84,1);
        real_dur = datarun{run_number}.triggers(one_repeat_end)*1000/refresh_rep(run_number); % duration of 1 repeat in frames
        inputs_rep_tmp = inputs_rep{run_number}(my_stixels,1:floor(real_dur));
        
        
        % filter raw inputs with STA - get generator signal
        my_filtered_inputs=zeros(size(inputs_rep_tmp,1),size(inputs_rep_tmp,2)-size(unbiased_sta,2)+1);
        for i=1:size(inputs_rep_tmp,1)
            my_filtered_inputs(i,:)=conv(inputs_rep_tmp(i,:),unbiased_sta(i,:),'valid');
        end
        gen_sig=sum(my_filtered_inputs,1);
        gen_sig=gen_sig/max(gen_sig);
        
        % pass generator signal through nonlinearity
        predicted_rate = zeros(size(gen_sig));
        for j=1:size(gensig_bins,2)-1
            a = find(gen_sig >= gensig_bins(j) & gen_sig < gensig_bins(j+1));
            predicted_rate(a)= nonlinearity(j);
        end
        
        
        mean_ndf3_norm_cut = mean_ndf3_norm(sta_params.length:end-14)';
        
        predicted_norm= predicted_rate/mean(predicted_rate);
        sse = sum((predicted_norm-mean_ndf3_norm_cut).^2);
        sst = sum((mean_ndf3_norm_cut-mean(mean_ndf3_norm_cut)).^2);
        r2_3by3model(cnt) = 1 - sse/sst;
        cell_for_use = [cell_for_use cellID];
        cnt=cnt+1;
    end
    
end

clear r2_4by3 r2_3by4
cnt = 1;
for cellID = cell_for_use
    
    run_number = 1;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf4=mean_sr/length(trial_begins);
    mean_ndf4_norm = mean_ndf4/mean(mean_ndf4);
    
    
    run_number = 2;
    
    trial_begins=round([0; datarun{run_number}.triggers(find(diff(datarun{run_number}.triggers)>0.84)+1)]*1000/refresh(run_number));
    spikes=ceil((datarun{run_number}.spikes{cellID}-datarun{run_number}.triggers(1))*1000/refresh(run_number));
    spikes(spikes<sta_params.length-sta_params.offset)=[];
    spike_rate=zeros(max(diff(trial_begins))*20,1);
    spikes(spikes>max(diff(trial_begins))*20)=[];
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    mean_sr=0;
    for i=1:length(trial_begins)
        mean_sr = mean_sr + spike_rate(trial_begins(i)+1:trial_begins(i)+max(diff(trial_begins)));
    end
    mean_ndf3=mean_sr/length(trial_begins);
    mean_ndf3_norm = mean_ndf3/mean(mean_ndf3);
    
%     figure
%     plot(mean_ndf3_norm)
%     hold on
%     plot(mean_ndf4_norm)
        
    % how much of ndf4 variance is explained by ndf3
    sse = sum((mean_ndf4_norm-mean_ndf3_norm).^2);
    sst = sum((mean_ndf3_norm-mean(mean_ndf3_norm)).^2);
    r2_3by4(cnt) = 1 - sse/sst;
    
    % how much of ndf3 variance is explained by ndf4
    sse = sum((mean_ndf3_norm-mean_ndf4_norm).^2);
    sst = sum((mean_ndf4_norm-mean(mean_ndf4_norm)).^2);
    r2_4by3(cnt) = 1 - sse/sst;
    
    cnt = cnt+1;
end

figure
hold on
% ON parasols
plot(r2_3by4, r2_3by3model,'b*')

% on midgets
plot(r2_3by4, r2_3by3model,'+r')

% off parasols
plot(r2_3by4, r2_3by3model,'vk')

% off midgets
plot(r2_3by4, r2_3by3model,'^g')

axis([-1 1 -1 1])
line([-1,1], [0,0])
line([0,0], [-1,1])
line([-1,1], [-1,1])
xlabel('wrong light level')
ylabel('model')

