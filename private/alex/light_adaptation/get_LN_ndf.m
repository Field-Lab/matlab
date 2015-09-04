function [predicted_spike_rate_acc, spike_rate_acc, gen_sig_acc, r2, cntstx_acc, cell_ind] = get_LN_ndf(data, scale, mean_rgb)


downsize=1;
offset = 0;
int_method='spline';
sta_length=15;
sigstix_threshold = 4;
movie_duration = 31000; % ms,single repetition
nm_repeats = 20;


mvpath='/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/';
my_movie=zeros(320,320,3600);
cnt=1;
for i=1:30
    load([mvpath, 'movie_chunk_', int2str(i), '.mat']);
    movie=movie(81:240,:,:);
    movie=imresize(movie,2,'nearest');
    my_movie(:,:,cnt:cnt+119)=movie;
    cnt=cnt+120;
end
my_movie=my_movie/255-0.5;
clear mvpath movie

for i=1
    switch data
        case '004'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            frames=8;
            ndf=4;
            nmdata='009';
        case '008'
            wn_movie_name = 'BW-10-8-0.48-11111-32x32.xml';
            frames=8;
            ndf=4;
            nmdata='009';
        case '011'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            frames=8;
            ndf=3;
            nmdata='012';
        case '014'
            wn_movie_name = 'BW-8-8-0.48-11111-40x40.xml';
            frames=8;
            ndf=3;
            nmdata='012';
        case '015'
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            frames=6;
            ndf=2;
            nmdata='016';
        case '018'
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            frames=6;
            ndf=2;
            nmdata='016';
        case '019'
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            frames=6;
            ndf=1;
            nmdata='020';
        case '022'
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            frames=6;
            ndf=1;
            nmdata='020';
        case '023'
            wn_movie_name = 'BW-16-4-0.48-11111-20x20.xml';
            frames=4;
            ndf=0;
            nmdata='024';
        case '027'
            wn_movie_name = 'BW-4-2-0.48-11111-80x80.xml';
            frames=2;
            ndf=0;
            nmdata='024';
    end
    starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'/data',data]);
    starun = load_params(starun,'verbose',1);
    starun = set_polarities(starun);
    starun = load_neurons(starun);
    starun = load_sta(starun,'load_sta','all','keep_java_sta',true);
    
    nmrun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',nmdata,'/data',nmdata]);
    nmrun = load_params(nmrun,'verbose',1);
    nmrun = load_neurons(nmrun);
end
clear nmdata

[inputs, refresh, duration] = get_wn_movie_ath(starun, wn_movie_name);
clear wn_movie_name

scaled_movie=zeros(size(inputs,1), 3600);
stix_size=320/sqrt(size(inputs,1));
cnt=1;
for i=1:stix_size:320
    for j=1:stix_size:320
        tmp=reshape(my_movie(i:i+stix_size-1,j:j+stix_size-1, :), stix_size*stix_size, 3600);
        scaled_movie(cnt,:)=mean(tmp);
        cnt=cnt+1;
    end
end
clear stix_size

% figure
tc=zeros(length(starun.cell_types),length(starun.vision.timecourses(1).g));
for j=1:size(tc,1)
    cnt=1;
   % subplot(3,4,j)
    for i=starun.cell_types{j}.cell_ids
        tmp=find(starun.cell_ids==i);
        if ~isempty(starun.vision.timecourses(tmp).g)
            if mean_rgb % for rod-dominated RGB runs, take average of green and blue channels
                tc(j,:)=tc(j,:)+mean([starun.vision.timecourses(tmp).b'; starun.vision.timecourses(tmp).g']);
            elseif  max(starun.vision.timecourses(tmp).g)<max(starun.vision.timecourses(tmp).b) % for blue cells, RGB run                
                tc(j,:)=tc(j,:)+starun.vision.timecourses(tmp).b';
            else % BW run, all cells, or RGB run, green cells
                %if sum(starun.vision.timecourses(tmp).g == starun.vision.timecourses(tmp).b)
                tc(j,:)=tc(j,:)+starun.vision.timecourses(tmp).g';
            end
%              plot(starun.vision.timecourses(tmp).g','g')
%              hold on
%              plot(starun.vision.timecourses(tmp).b','b')
            cnt=cnt+1;
        end
    end
    tc(j,:)=tc(j,:)'/cnt;
  
end

nm_triggers=nmrun.triggers;
movie_refresh=refresh/frames;%1/120*1000; %ms

predicted_spike_rate_acc = [];
spike_rate_acc = [];
gen_sig_acc = [];
r2=zeros(length(starun.cell_ids),2);
cntstx_acc = r2;
cell_ind = zeros(length(starun.cell_ids),1);

cell_cnt=1;
for cellID = starun.cell_ids
    datarunID=find(starun.cell_ids==cellID);
    
    [~, my_type] = find_cell_type(starun, cellID);
    
    %****** PART 1 calculate linear filter from STA ******
    % sta of a chosen cell
    sta = double(squeeze(starun.stas.stas{datarunID}));
    
    cell_ind(cell_cnt) = datarunID;
    
    if ndims(sta)==4
        if mean_rgb
            sta=squeeze(mean(sta,3));
        elseif max(starun.vision.timecourses(datarunID).g)<max(starun.vision.timecourses(datarunID).b)
            sta=squeeze(sta(:,:,3,:)); % take blue channel
        else
            sta=squeeze(sta(:,:,2,:)); % take green channel
        end
    end
    
    % get weights by convolving sta of a chosen cell with mean time course
    temp_sta = reshape(sta, size(sta,1)*size(sta,2),size(sta,3));
    weights=sum(repmat(tc(my_type,:),size(sta,1)*size(sta,2),1) .* temp_sta, 2);
    weights = weights/max(weights);
    
    % threshold weights
    weight_threshold=robust_std(weights)*sigstix_threshold;
    inds=find(weights>weight_threshold);
    cntstx = 0;
    while length(inds)<2
        cntstx=cntstx+0.5;
        weight_threshold=robust_std(weights)*(sigstix_threshold-cntstx);
        inds=find(weights>weight_threshold);
    end
    
    cntstx_acc(cell_cnt,1) = cntstx;
    cntstx_acc(cell_cnt,2) = length(inds);
    
    my_inputs = inputs(inds, :); % to save memory and time
    inp= zeros(length(inds),size(inputs,2)*downsize);
    % down to a frame: instead of refresh of 66.6192, use refresh/8 = 8.3274
    for i=1:downsize
        inp(:,i:downsize:end) = my_inputs;
    end
    my_inputs=inp;
        
    spikes=ceil((starun.spikes{datarunID}-starun.triggers(1))*1000/(refresh/downsize)); % spikes in ms
    spikes(spikes<sta_length-offset)=[];
    
    spike_rate=zeros(duration*downsize,1);
    my_sta=zeros(length(inds),sta_length);
    
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        for j=1:sta_length
            my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(my_inputs(:,spikes(ia)-sta_length+j+offset),2);
        end
        spikes(ia)=[];
    end
    my_sta=my_sta/sum(spike_rate);    

    
    % convolve raw input with linear filter
    my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(my_sta,2)+1);
    for i=1:length(inds)
        my_filtered_inputs(i,:)=conv(my_inputs(i,:),my_sta(i,:),'valid');
    end
    
    gen_sig=sum(my_filtered_inputs);
    gen_sig=gen_sig/max(gen_sig);
    
    spike_rate = spike_rate(sta_length:end);
    
    bins = (max(gen_sig)-min(gen_sig))/100;
    my_bins = min(gen_sig):bins:max(gen_sig)-bins;
    
    % calculate nonlinearity
    my_nl=zeros(size(my_bins));
    cnt=1;
    for i=my_bins
        tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
        tmp(tmp<1)=[];
        tmp(tmp>length(spike_rate))=[];
        my_nl(cnt)=mean(spike_rate(tmp)); % find mean firing rate after this certain value
        cnt=cnt+1;
    end
    
    figure
    plot(my_nl/2)
    
    % calculate predicted FR
    cnt=1;
    predicted_rate=zeros(size(spike_rate));
    for i=min(gen_sig):bins:max(gen_sig)-bins
        tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
        tmp(tmp<1)=[];
        tmp(tmp>length(spike_rate))=[];
        predicted_rate(tmp)=my_nl(cnt);
        cnt=cnt+1;
    end
    sse = sum((predicted_rate-spike_rate).^2);
    sst = sum((spike_rate-mean(spike_rate)).^2);
    r2(cell_cnt,1) = 1 - sse/sst;
    
    figure
    plot(predicted_rate)
    hold on
    plot(spike_rate)
    
    %****** PART 2 calculate gen signal from NSEM using STA ******
    
    spikes=ceil((nmrun.spikes{datarunID}-nmrun.triggers(1))*1000/movie_refresh); % spikes in ms
    spikes(spikes<sta_length-offset)=[];
    
    spike_rate=zeros(movie_duration*nm_repeats,1);
    
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end
    
    
    my_begs=round([0; nm_triggers(find(diff(nm_triggers)>0.9)+1)]*1000/movie_refresh);
    acc_spike_rate=0;
    
    for i=1:length(my_begs)-1
        acc_spike_rate = acc_spike_rate + spike_rate(my_begs(i)+1:my_begs(i)+3600);
    end
    acc_spike_rate=acc_spike_rate/(length(my_begs)-1);
    
    % interpolate linear filter to match movie frames
    if downsize==1
        tt=[];
        for i=1:size(my_sta,1)
            tt=[tt; interp1(1:8:8*size(my_sta,2), my_sta(i,:),1:8*size(my_sta,2),int_method)]; % spline, linear
        end
        my_sta=tt;
    end
    

    
    % convolve raw input with linear filter
    my_inputs=scaled_movie(inds, :);
    my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(my_sta,2)+1);
    for i=1:length(inds)
        my_filtered_inputs(i,:)=conv(my_inputs(i,:),my_sta(i,:),'valid');
    end
    
    gen_sig=sum(my_filtered_inputs);
    gen_sig=gen_sig/max(gen_sig);

    
    gen_sig_acc = [gen_sig_acc; gen_sig];
    
    acc_spike_rate = acc_spike_rate(size(my_sta,2):end);
    
    bins = (max(gen_sig)-min(gen_sig))/100;
    my_bins = min(gen_sig):bins:max(gen_sig)-bins;
    
    % calculate nonlinearity
    my_nl=zeros(size(my_bins));
    cnt=1;
    for i=my_bins
        tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
        tmp(tmp<1)=[];
        tmp(tmp>length(acc_spike_rate))=[];
        my_nl(cnt)=mean(acc_spike_rate(tmp)); % find mean firing rate after this certain value
        cnt=cnt+1;
    end
    
    % calculate predicted FR
    cnt=1;
    predicted_rate=zeros(size(acc_spike_rate));
    for i=min(gen_sig):bins:max(gen_sig)-bins
        tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
        tmp(tmp<1)=[];
        tmp(tmp>length(acc_spike_rate))=[];
        predicted_rate(tmp)=my_nl(cnt);
        cnt=cnt+1;
    end
    sse = sum((predicted_rate-acc_spike_rate).^2);
    sst = sum((acc_spike_rate-mean(acc_spike_rate)).^2);
    r2(cell_cnt,2) = 1 - sse/sst;
    
    predicted_spike_rate_acc = [predicted_spike_rate_acc predicted_rate];
    spike_rate_acc = [spike_rate_acc acc_spike_rate];
    
    cell_cnt=cell_cnt+1;
end
